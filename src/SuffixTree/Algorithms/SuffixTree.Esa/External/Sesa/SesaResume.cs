using System.Diagnostics;
using System.Runtime.InteropServices;
using Microsoft.Win32.SafeHandles;

namespace SuffixTree.Esa.External;

/// <summary>
/// The recording and the verification of the stages of one build with <see cref="SesaBuildOptions.Resume"/> (work package 1b.9e, ADR-004 section 4.3 "Crash safety and atomicity"). <see cref="Load"/> reads
/// <c>state.json</c> of an earlier run, checks it against this run's input and layout and then, stage by stage, the outputs that each completed stage vouches for (the ranges of the container and the scratch files: length
/// and CRC-32C); the stages up to the first one that does not verify are <em>kept</em> and are not run again. <see cref="Prepare"/> removes what the stages that are not kept left behind. Every stage that runs then
/// records itself with <c>Commit</c>: its outputs are forced to storage (when <see cref="SesaBuildOptions.Fsync"/>), their CRC-32C is computed, and only then the new state is written (a temporary file, forced to storage,
/// renamed over the old one), so that a state never names an output that is not on the disk. Nothing here is on the path of a build without Resume.
/// </summary>
internal sealed class SesaResumeRun
{
    private const string RanksFile = "ranks.bin";
    private const string BwtFile = "ph6/bwt.bin";

    private readonly SesaBuildOptions _options;
    private readonly string _tmp;
    private readonly string _scratch;
    private readonly PackedTextInfo _info;
    private readonly SesaPlan _plan;
    private readonly bool _psi;
    private readonly ChildCodec _codec;
    private readonly IoCounters _counters;
    private readonly DifferenceCover? _cover;
    private readonly List<SesaCheckpointRow> _rows = [];
    private SesaResumeState _state;
    private List<SesaStageRecord> _skipped = [];
    private int _kept;
    private int _onDiskStages;
    private bool _onDisk;
    private int _writes;
    private long _bytes;
    private double _checkpointSeconds;
    private ClearResult _cleared;
    private double _clearSeconds;

    private SesaResumeRun(SesaBuildOptions options, string tmp, string scratch, PackedTextInfo info, SesaPlan plan, ChildCodec codec, IoCounters counters, DifferenceCover? cover)
    {
        _options = options;
        _tmp = tmp;
        _scratch = scratch;
        _info = info;
        _plan = plan;
        _psi = options.IncludePsi;
        _codec = codec;
        _counters = counters;
        _cover = cover;
        _state = new SesaResumeState
        {
            N = info.N,
            Bits = info.Bits,
            Sigma = info.Sigma,
            Sha256 = info.Sha256,
            Letters = [.. info.Letters.Select(l => (int)l)],
            Counts = [.. info.Counts],
            Psi = _psi,
            ChildCodec = codec == ChildCodec.Tagged ? "tagged" : "packed",
            ReservedLength = plan.ReservedLength,
        };
    }

    /// <summary>The difference cover of this run (null when it cannot be made: the build reports that when it tries).</summary>
    internal DifferenceCover? Cover => _cover;

    /// <summary>A state file was found (whether or not it was usable).</summary>
    internal bool StateFound { get; private set; }

    /// <summary>Stages of an earlier run were verified and are not run again.</summary>
    internal bool Accepted => _skipped.Count > 0;

    /// <summary>Why the state was not used, or why it was cut short; empty when everything it held was used.</summary>
    internal string Reason { get; private set; } = string.Empty;

    /// <summary>Seconds that <see cref="Load"/> spent reading and verifying.</summary>
    internal double VerifySeconds { get; private set; }

    /// <summary>The names of the stages that are kept, in order.</summary>
    internal IReadOnlyList<string> Skipped => [.. _skipped.Select(s => s.Name)];

    /// <summary>Bytes that <see cref="Prepare"/> frees (partial files of the stages that were not completed, the part of the container behind its reserved length): they count towards the free space of the pre-flight.</summary>
    internal long ReclaimableBytes { get; private set; }

    /// <summary>Bytes of the scratch files that the kept stages vouch for (they exist already and count against the scratch that the rest of the build needs).</summary>
    internal long ExistingScratchBytes { get; private set; }

    /// <summary>A state file that this run trusts or wrote is on the disk, so a failure of the build leaves its temporary file and scratch for the next call.</summary>
    internal bool Keeps => _onDisk;

    /// <summary>Seconds spent recording stages so far (fsync, CRC, state file).</summary>
    internal double CheckpointSeconds => _checkpointSeconds;

    /// <summary>
    /// The ranges of the container that the kept stages vouch for, each with a CRC-32C that <see cref="Load"/> has just verified (work package 1b.9e.1). Every other byte of the temporary container is, for this run, nothing that an earlier
    /// run or anyone else can be believed about: <see cref="MappedOutputFile.OpenExisting"/> makes it zero, as it is in a new file, before any stage runs.
    /// </summary>
    internal IReadOnlyList<(long Offset, long Length)> KeptContainerRanges() => [.. _skipped.SelectMany(s => s.Checks).Where(c => c.InContainer).Select(c => (c.Offset, c.Length))];

    /// <summary>Notes what the normalisation of the reopened container did and how long it took (for the report).</summary>
    internal void NoteCleared(ClearResult cleared, double seconds)
    {
        _cleared = cleared;
        _clearSeconds = seconds;
    }

    /// <summary>The stage was completed by an earlier run and its outputs verified: this run does not run it.</summary>
    internal bool Has(string stage) => _skipped.Any(s => s.Name == stage);

    internal SesaStageRecord Record(string stage) => _skipped.First(s => s.Name == stage);

    /// <summary>The scratch file that a kept stage vouches for (false while it is not recorded or not kept): a failure must not delete it.</summary>
    internal bool Vouches(string relativePath) => _state.Stages.Any(s => s.Checks.Any(c => !c.InContainer && c.Role == relativePath));

    /// <summary>The report for the result.</summary>
    internal SesaResumeReport Report() => new()
    {
        StateFound = StateFound,
        Accepted = Accepted,
        Reason = Reason,
        Skipped = Skipped,
        VerifySeconds = VerifySeconds,
        Checkpoints = [.. _rows],
        StateWrites = _writes,
        StateBytes = _bytes,
        ClearScannedBytes = _cleared.ScannedBytes,
        ClearedBytes = _cleared.ClearedBytes,
        ClearSeconds = _clearSeconds,
    };

    private string StatePath => Path.Combine(_scratch, SesaResumeFormat.FileName);

    // ------------------------------------------------------------------ reading and verifying an earlier run

    /// <summary>Reads and verifies the state of an earlier run of this build; never throws because of what it finds (any damage is a reason to start clean or to redo a stage).</summary>
    internal static SesaResumeRun Load(SesaBuildOptions options, string tmp, string scratch, PackedTextInfo info, SesaPlan plan, ChildCodec codec, IoCounters counters)
    {
        DifferenceCover? cover = null;
        try
        {
            cover = SaBuilder.CoverFor(options.Period, options.ExactPeriod);
        }
        catch (Exception e) when (e is ArgumentException or InvalidOperationException)
        {
            // the build reports it when it makes the cover; a state of the ranks is then simply not reusable
        }
        var run = new SesaResumeRun(options, tmp, scratch, info, plan, codec, counters, cover);
        long t0 = Stopwatch.GetTimestamp();
        run.Verify();
        run.VerifySeconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        return run;
    }

    private void Verify()
    {
        StateFound = File.Exists(StatePath) || Directory.Exists(StatePath);
        if (!StateFound)
        {
            Reason = "no state file";
            return;
        }
        byte[] bytes;
        try
        {
            bytes = File.ReadAllBytes(StatePath);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            Reason = "the state file cannot be read";
            return;
        }
        if (!SesaResumeFormat.TryParse(bytes, out SesaResumeState? found) || found is null)
        {
            Reason = "the state file is damaged, cut, or of another format";
            return;
        }
        string? mismatch = HeaderMismatch(found);
        if (mismatch is not null)
        {
            Reason = mismatch;
            return;
        }
        if (found.Stages.Count == 0)
        {
            Reason = "the state records no completed stage";
            return;
        }
        long length;
        try
        {
            var fi = new FileInfo(_tmp);
            length = fi.Exists ? fi.Length : -1;
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            length = -1;
        }
        if (length < _plan.ReservedLength)
        {
            Reason = length < 0 ? "the temporary file of the earlier run is gone" : $"the temporary file is shorter ({length} bytes) than the container was created ({_plan.ReservedLength} bytes)";
            return;
        }
        _state = found;
        _onDisk = true;
        _onDiskStages = found.Stages.Count;
        int kept = 0;
        foreach (SesaStageRecord st in found.Stages)
        {
            string? why = VerifyStage(st, length);
            if (why is not null)
            {
                Reason = $"stage {st.Name}: {why}; it and the stages after it are done again";
                break;
            }
            kept++;
        }
        _kept = kept;
        _skipped = [.. found.Stages.Take(kept)];
        if (kept == 0)
        {
            _state = _state with { Stages = [] };
            _onDisk = false;
            _onDiskStages = 0;
            return;
        }
        ReclaimableBytes = Has(SesaStages.Lcp) ? 0 : Math.Max(0, length - _plan.ReservedLength);
        foreach (SesaFileCheck f in _skipped.SelectMany(s => s.Checks).Where(c => !c.InContainer)) ExistingScratchBytes += f.Length;
        ReclaimableBytes += PartialScratchBytes();
    }

    private string? HeaderMismatch(SesaResumeState s)
    {
        if (s.Writer != SesaResumeFormat.WriterId) return "the state was written by another writer";
        if (s.N != _info.N || s.Bits != _info.Bits || s.Sigma != _info.Sigma || s.Sha256 != _info.Sha256) return "the state is of another input (length or SHA-256 differ)";
        if (!s.Letters.SequenceEqual(_info.Letters.Select(l => (int)l)) || !s.Counts.SequenceEqual(_info.Counts)) return "the state is of another input (alphabet or symbol counts differ)";
        if (s.Psi != _psi) return "the state is of a build with another PSI option";
        if (s.ChildCodec != (_codec == ChildCodec.Tagged ? "tagged" : "packed")) return "the state is of a build with another CHILD codec";
        if (s.ReservedLength != _plan.ReservedLength) return "the state is of another layout of the container";
        return null;
    }

    private List<SesaFileCheck> Expected(SesaStageRecord st, out string? why)
    {
        why = null;
        var list = new List<SesaFileCheck>();
        switch (st.Name)
        {
            case SesaStages.Packed:
                list.Add(new SesaFileCheck("text", true, _plan.TextOffset, _plan.TextLength, 0));
                list.Add(new SesaFileCheck("alphabet", true, _plan.AlphabetOffset, _plan.AlphabetLength, 0));
                list.Add(new SesaFileCheck("counts", true, _plan.CountsOffset, _plan.CountsLength, 0));
                break;
            case SesaStages.Ranks:
                if (_cover is null) { why = "the difference cover of this run cannot be made"; break; }
                if (!st.Values.TryGetValue("period", out long period) || period != _cover.Period || st.Value("kind") != (long)_cover.Kind || st.Value("residues") != _cover.Residues.Length || st.Value("digest") != CoverDigest(_cover))
                {
                    why = "the rank file was built for another difference cover (another period)";
                    break;
                }
                list.Add(new SesaFileCheck(RanksFile, false, 0, 4 * new SampleGeometry(_info.N, _cover).Slots, 0));
                break;
            case SesaStages.Sa:
                list.Add(new SesaFileCheck("sa", true, _plan.SaOffset, _plan.SaLength, 0));
                break;
            case SesaStages.Lcp:
                {
                    long kind = st.Value("kind"), overflow = st.Value("overflow");
                    if (kind is < 0 or > (long)LcpKind.U32 || overflow < 0 || (overflow > 0 && kind == (long)LcpKind.U32)) { why = "the width or the escapes of the LCP are impossible"; break; }
                    SesaTail tail = SesaLayout.PlaceTail(_plan, (LcpKind)kind, overflow);
                    long nodes = st.Value("node_count"), internalNodes = st.Value("internal_count");
                    if (st.Value("final_length") != tail.FinalLength) { why = "the recorded final length of the container does not follow from the width and the escapes"; break; }
                    long expectedNodes = _info.N == 0 ? 1 : internalNodes + _info.N + 2; // the scan of the empty text counts the root only
                    if (nodes != expectedNodes || internalNodes < 0) { why = "NODE_COUNT is not 1 + I + n + 1"; break; }
                    if (st.Value("lrs_depth") < 0 || st.Value("lrs_depth") > _info.N || st.Value("max_lcp") < 0 || st.Value("max_lcp") > _info.N || st.Value("lrs_rank") < -1 || st.Value("lrs_rank") > _info.N) { why = "the statistics of the scan are out of range"; break; }
                    list.Add(new SesaFileCheck("child", true, _plan.ChildOffset, _plan.ChildLength, 0));
                    list.Add(new SesaFileCheck("lcp", true, _plan.LcpOffset, tail.LcpLength, 0));
                    if (overflow > 0) list.Add(new SesaFileCheck("overflow", true, tail.OverflowOffset, tail.OverflowLength, 0));
                    break;
                }
            case SesaStages.Bwt:
                if (!PrimaryRowPossible(st)) { why = "the primary row is impossible"; break; }
                list.Add(new SesaFileCheck(BwtFile, false, 0, _info.N + 1, 0));
                break;
            case SesaStages.Psi:
                if (!PrimaryRowPossible(st)) { why = "the primary row is impossible"; break; }
                list.Add(new SesaFileCheck("psi", true, _plan.PsiOffset, _plan.PsiLength, 0));
                break;
            default:
                why = "an unknown stage";
                break;
        }
        return list;
    }

    private bool PrimaryRowPossible(SesaStageRecord st) => st.Values.TryGetValue("primary_row", out long row) && row >= 0 && row <= _info.N;

    /// <summary>Null when the stage's outputs are what its record says; else the reason in words.</summary>
    private string? VerifyStage(SesaStageRecord st, long containerLength)
    {
        try
        {
            List<SesaFileCheck> expected = Expected(st, out string? why);
            if (why is not null) return why;
            if (st.Checks.Count != expected.Count) return "the record lists other outputs than this layout has";
            for (int i = 0; i < expected.Count; i++)
            {
                SesaFileCheck a = st.Checks[i], b = expected[i];
                if (a.Role != b.Role || a.InContainer != b.InContainer || a.Offset != b.Offset || a.Length != b.Length) return $"the record of {a.Role} is not where or as long as this layout has it";
            }
            if (st.Name == SesaStages.Lcp && containerLength != st.Value("final_length")) return $"the container has {containerLength} bytes, the stage grew it to {st.Value("final_length")}";
            var ranges = new List<(int Kind, long Offset, long Length)>();
            foreach (SesaFileCheck c in st.Checks.Where(c => c.InContainer)) ranges.Add((0, c.Offset, c.Length));
            if (ranges.Count > 0)
            {
                (uint Crc, SesaCrcRow Row)[] sums = SesaChecksums.Compute(_tmp, ranges);
                int at = 0;
                foreach (SesaFileCheck c in st.Checks.Where(c => c.InContainer))
                {
                    if (sums[at++].Crc != c.Crc) return $"the {c.Role} section has another CRC-32C than when the stage finished";
                }
            }
            foreach (SesaFileCheck c in st.Checks.Where(c => !c.InContainer))
            {
                string path = Path.Combine(_scratch, c.Role);
                var fi = new FileInfo(path);
                if (!fi.Exists) return $"the scratch file {c.Role} is gone";
                if (fi.Length != c.Length) return $"the scratch file {c.Role} has {fi.Length} bytes, the stage wrote {c.Length}";
                if (SesaChecksums.Compute(path, [(0, 0, c.Length)])[0].Crc != c.Crc) return $"the scratch file {c.Role} has another CRC-32C than when the stage finished";
            }
            return null;
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException or InvalidDataException or KeyNotFoundException or ArgumentException)
        {
            return "its outputs cannot be read: " + e.Message;
        }
    }

    private long PartialScratchBytes()
    {
        long total = 0;
        try
        {
            foreach (string file in Directory.EnumerateFiles(_scratch, "*", SearchOption.AllDirectories))
            {
                if (IsKept(file)) continue;
                total += new FileInfo(file).Length;
            }
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // an estimate that only makes the pre-flight stricter
        }
        return total;
    }

    private bool IsKept(string fullPath)
    {
        string relative = Path.GetRelativePath(_scratch, fullPath).Replace('\\', '/');
        if (relative == SesaResumeFormat.FileName) return true;
        return _skipped.Any(s => s.Checks.Any(c => !c.InContainer && c.Role == relative));
    }

    /// <summary>Removes what the stages that are not kept left in the scratch (partial files, directories of PH4, a half-written state) and rewrites the state when it names stages that are not kept.</summary>
    internal void Prepare()
    {
        if (!Accepted) return;
        foreach (string dir in Directory.EnumerateDirectories(_scratch).ToArray())
        {
            string rel = Path.GetRelativePath(_scratch, dir).Replace('\\', '/');
            bool holdsKept = _skipped.Any(s => s.Checks.Any(c => !c.InContainer && c.Role.StartsWith(rel + "/", StringComparison.Ordinal)));
            if (!holdsKept) Directory.Delete(dir, true);
            else
            {
                foreach (string file in Directory.EnumerateFiles(dir, "*", SearchOption.AllDirectories).ToArray())
                {
                    if (!IsKept(file)) File.Delete(file);
                }
            }
        }
        foreach (string file in Directory.EnumerateFiles(_scratch).ToArray())
        {
            if (!IsKept(file)) File.Delete(file);
        }
        if (_kept < _onDiskStages)
        {
            _state = _state.Truncate(_kept);
            WriteState(_state);
            _onDiskStages = _kept;
        }
    }

    // ------------------------------------------------------------------ recording a stage

    /// <summary>The PH0 stage: the container is complete up to TEXT, ALPHABET and COUNTS.</summary>
    internal void CommitPacked(MappedOutputFile container) => Commit(SesaStages.Packed, container,
        [("text", true, _plan.TextOffset, _plan.TextLength), ("alphabet", true, _plan.AlphabetOffset, _plan.AlphabetLength), ("counts", true, _plan.CountsOffset, _plan.CountsLength)], new Dictionary<string, long>());

    /// <summary>The PH1 stage: the rank file of the difference-cover sample.</summary>
    internal void CommitRanks(DifferenceCover cover)
    {
        long bytes = new FileInfo(Path.Combine(_scratch, RanksFile)).Length;
        Commit(SesaStages.Ranks, null, [(RanksFile, false, 0, bytes)], new Dictionary<string, long>
        {
            ["period"] = cover.Period,
            ["kind"] = (long)cover.Kind,
            ["residues"] = cover.Residues.Length,
            ["digest"] = CoverDigest(cover),
        });
    }

    /// <summary>The PH2 / PH3 stage: the SA section.</summary>
    internal void CommitSa(MappedOutputFile container) => Commit(SesaStages.Sa, container, [("sa", true, _plan.SaOffset, _plan.SaLength)], new Dictionary<string, long>());

    /// <summary>The PH4 / PH5 stage: the container has grown, the LCP, the overflow list and CHILD are written.</summary>
    internal void CommitLcp(MappedOutputFile container, SesaTail tail, ChildScanResult scan)
    {
        var outputs = new List<(string, bool, long, long)>
        {
            ("child", true, _plan.ChildOffset, _plan.ChildLength),
            ("lcp", true, _plan.LcpOffset, tail.LcpLength),
        };
        if (tail.OverflowCount > 0) outputs.Add(("overflow", true, tail.OverflowOffset, tail.OverflowLength));
        Commit(SesaStages.Lcp, container, outputs, new Dictionary<string, long>
        {
            ["kind"] = (long)tail.Kind,
            ["overflow"] = tail.OverflowCount,
            ["final_length"] = tail.FinalLength,
            ["node_count"] = scan.NodeCount,
            ["internal_count"] = scan.InternalCount,
            ["lrs_depth"] = scan.LrsDepth,
            ["lrs_rank"] = scan.LrsDepth == 0 ? -1 : scan.LrsRank,
            ["max_lcp"] = scan.MaxLcp,
        });
    }

    /// <summary>The PH6 stage a: the BWT file (scratch).</summary>
    internal void CommitBwt(long primaryRow) => Commit(SesaStages.Bwt, null, [(BwtFile, false, 0, _info.N + 1)], new Dictionary<string, long> { ["primary_row"] = primaryRow });

    /// <summary>The PH6 stage b: the PSI section (the BWT file is not needed any more).</summary>
    internal void CommitPsi(MappedOutputFile container, long primaryRow) => Commit(SesaStages.Psi, container, [("psi", true, _plan.PsiOffset, _plan.PsiLength)], new Dictionary<string, long> { ["primary_row"] = primaryRow });

    /// <summary>The primary row that a kept BWT or PSI stage recorded.</summary>
    internal long PrimaryRow() => (Has(SesaStages.Psi) ? Record(SesaStages.Psi) : Record(SesaStages.Bwt)).Value("primary_row");

    private void Commit(string stage, MappedOutputFile? container, IReadOnlyList<(string Role, bool InContainer, long Offset, long Length)> outputs, Dictionary<string, long> values)
    {
        long t0 = Stopwatch.GetTimestamp();
        if (container is not null) container.Complete(_options.Fsync, _counters);
        foreach ((string role, bool inContainer, _, _) in outputs)
        {
            if (inContainer || !_options.Fsync) continue;
            using SafeFileHandle handle = File.OpenHandle(Path.Combine(_scratch, role), FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite);
            RandomAccess.FlushToDisk(handle);
            _counters.AddFlush();
        }
        long t1 = Stopwatch.GetTimestamp();

        var checks = new List<SesaFileCheck>();
        long crcBytes = 0;
        var inside = outputs.Where(o => o.InContainer).ToList();
        if (inside.Count > 0)
        {
            (uint Crc, SesaCrcRow Row)[] sums = SesaChecksums.Compute(_tmp, [.. inside.Select(o => (0, o.Offset, o.Length))]);
            for (int i = 0; i < inside.Count; i++) checks.Add(new SesaFileCheck(inside[i].Role, true, inside[i].Offset, inside[i].Length, sums[i].Crc));
            crcBytes += inside.Sum(o => o.Length);
        }
        foreach ((string role, _, _, long length) in outputs.Where(o => !o.InContainer))
        {
            checks.Add(new SesaFileCheck(role, false, 0, length, SesaChecksums.Compute(Path.Combine(_scratch, role), [(0, 0, length)])[0].Crc));
            crcBytes += length;
        }
        long t2 = Stopwatch.GetTimestamp();

        var record = new SesaStageRecord(stage, checks, values);
        SesaResumeState next = _state.Truncate(_kept).With(record);
        WriteState(next);
        _state = next;
        _kept = next.Stages.Count;
        _onDiskStages = _kept;
        long t3 = Stopwatch.GetTimestamp();
        _rows.Add(new SesaCheckpointRow(stage, Seconds(t0, t1), Seconds(t1, t2), Seconds(t2, t3), crcBytes));
        _checkpointSeconds += Seconds(t0, t3);
    }

    private void WriteState(SesaResumeState state)
    {
        byte[] bytes = SesaResumeFormat.Serialize(state);
        string temp = Path.Combine(_scratch, SesaResumeFormat.TempName);
        using (var fs = new FileStream(temp, FileMode.Create, FileAccess.Write, FileShare.None, 1))
        {
            fs.Write(bytes);
            fs.Flush(_options.Fsync);
        }
        _counters.AddWrite(bytes.Length);
        if (_options.Fsync) _counters.AddFlush();
        File.Move(temp, StatePath, overwrite: true);
        _writes++;
        _bytes += bytes.Length;
        _onDisk = true;
    }

    private static double Seconds(long from, long to) => (to - from) / (double)Stopwatch.Frequency;

    /// <summary>A number that identifies the cover (the period and the residues): the CRC-32C of them as little-endian 32-bit integers.</summary>
    internal static long CoverDigest(DifferenceCover cover)
    {
        int[] words = [cover.Period, .. cover.Residues];
        return Crc32C.Compute(MemoryMarshal.AsBytes(words.AsSpan()));
    }
}
