using System.Diagnostics;
using System.Globalization;
using System.Runtime.InteropServices;
using System.Text;

namespace SuffixTree.Esa.External;

/// <summary>The names of the points at which <see cref="SesaBuildOptions.Hook"/> is called, in the order of the build (tests make the build fail, or wait, at each).</summary>
internal static class SesaHooks
{
    internal const string Start = "start";
    internal const string Scanned = "scanned";
    internal const string Created = "created";
    internal const string Packed = "packed";
    internal const string SaBuilt = "sa-built";
    internal const string LcpPlaced = "lcp-placed";
    internal const string LcpChildBuilt = "lcp-child-built";
    internal const string BwtBuilt = "bwt-built";
    internal const string PsiBuilt = "psi-built";
    internal const string InfoWritten = "info-written";
    internal const string Checksummed = "checksummed";
    internal const string HeaderWritten = "header-written";
    internal const string HeaderSynced = "header-synced";
    internal const string MagicWritten = "magic-written";

    /// <summary>Every hook point of a build with PSI, in order (a build without PSI skips <see cref="BwtBuilt"/> and <see cref="PsiBuilt"/>).</summary>
    internal static readonly string[] All =
    [
        Start, Scanned, Created, Packed, SaBuilt, LcpPlaced, LcpChildBuilt, BwtBuilt, PsiBuilt, InfoWritten, Checksummed, HeaderWritten, HeaderSynced, MagicWritten,
    ];
}

/// <summary>Explicit values for the sizes of PH6 (null = the default).</summary>
internal sealed record SesaPh6Tuning
{
    internal int? WindowBytes { get; init; }

    internal int? SymbolsPerPass { get; init; }

    internal int? ReadBytes { get; init; }
}

/// <summary>What <see cref="SesaBuilder.Build"/> needs besides the two paths.</summary>
internal sealed record SesaBuildOptions
{
    /// <summary>Worker threads of the bucket tasks of PH1 and PH3 and of the B2 tasks of PH4.</summary>
    internal int Threads { get; init; } = 4;

    /// <summary>The codec of the CHILD section (work package 1b.9a): null = by the length of the text (<see cref="SesaBuilder.ChildCodecFor"/>: tagged up to 2^30 - 1 symbols, else packed).</summary>
    internal ChildCodec? ChildCodec { get; init; }

    /// <summary>Write the optional PSI section (kind 12, <c>REQUIRED</c> clear). The cost is one more pass over the suffix array and the text, a scratch BWT of 1 B/char and 4 B/char in the file.</summary>
    internal bool IncludePsi { get; init; } = true;

    /// <summary>X of the difference cover (nominal; a Singer period nearby is used unless <see cref="ExactPeriod"/>: 993 for the default 1024).</summary>
    internal int Period { get; init; } = 1024;

    internal bool ExactPeriod { get; init; }

    /// <summary>The scale-down emulation of ADR-004 section 5 (E3 to E7): 3,000,000,000 for the 100 Mbp runs, 0 = the reference sizes.</summary>
    internal long ScaleN { get; init; }

    /// <summary>E7 series A: the n at which the scale-down rule is evaluated (0 = the text length).</summary>
    internal long RuleN { get; init; }

    /// <summary>Cap of the native buffers of each of PH1, PH2 / PH3 and PH4 (owner decision: 160 MiB).</summary>
    internal long ArenaCapBytes { get; init; } = 160L << 20;

    /// <summary>Cap of the native buffers of PH5 (owner decision: 64 MiB).</summary>
    internal long ChildArenaCapBytes { get; init; } = 64L << 20;

    /// <summary>Cap of the mapped windows and buffers of PH6.</summary>
    internal long Ph6CapBytes { get; init; } = 160L << 20;

    internal SaBuildTuning SaTuning { get; init; } = new();

    internal LcpChildTuning LcpTuning { get; init; } = new();

    internal SesaPh6Tuning Ph6Tuning { get; init; } = new();

    /// <summary>Symbols per block of the two passes over the input of PH0 (a positive multiple of 8).</summary>
    internal int TextBlockSymbols { get; init; } = 1 << 20;

    /// <summary>Window of the mapped writers of PH0.</summary>
    internal int TextWindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    /// <summary>Seed of the sampling; the output does not depend on it.</summary>
    internal ulong Seed { get; init; } = 1;

    /// <summary>The counters of the whole build (a private instance when null).</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>How the container and the BWT scratch are created (reservation of the disk space); null = <see cref="ReservationMode.Auto"/>.</summary>
    internal MappedFileOptions? Files { get; init; }

    /// <summary>Test hook: options of the arenas (guard bands, poison pattern); null = production arenas.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>fsync the file before the magic is written and after it (the atomic publication of ADR-003 section 6; off only in tests that do not care).</summary>
    internal bool Fsync { get; init; } = true;

    /// <summary>Read every slot of the rank file once between PH1 and PH2 / PH3 and check it (an O(S) pass).</summary>
    internal bool VerifyRanks { get; init; }

    /// <summary>
    /// Work package 1b.9e, ADR-004 section 4.3: record the build in <c>&lt;target&gt;.tmp.d/state.json</c> (the fingerprint of the input and what each completed stage wrote, every output forced to storage before its stage is
    /// marked) and, when a valid state of the same input is found, skip the completed stages whose outputs still verify (length and CRC-32C). A build that fails (an exception) keeps its temporary file and state for the next call;
    /// a state that is missing, cut, damaged, of another format, of another input or of another layout is discarded and the build starts clean. Takes the lock <c>&lt;target&gt;.lock</c>, so that two resumable builds of one
    /// target cannot run at once. False (the default): the build is exactly the build of 1b.7e-1b.9d, no state file, no lock, no extra fsync, and a failure removes everything. The output is the same either way.
    /// </summary>
    internal bool Resume { get; init; }

    /// <summary>Sample RssAnon while the build runs (the figure of the gate); unavailable where there is no <c>/proc</c>.</summary>
    internal bool MeasureRss { get; init; } = true;

    /// <summary>Where the kernel's text files come from (null = the real <c>/proc</c>); a test passes a source that fails like a platform without one.</summary>
    internal IProcSource? ProcSource { get; init; }

    /// <summary>The gate of the owner decision of 2026-11-18: RssAnon of the process at most 320 MB (10^6 bytes) at every n with the default parameters and four threads.</summary>
    internal double RssAnonCapMb { get; init; } = 320;

    /// <summary>
    /// Test hook, called with a name of <see cref="SesaHooks"/> at every boundary of the build; an exception thrown here fails the build like a failing step, a call that blocks holds the build at that point. With
    /// <see cref="Resume"/> a boundary is reached whether the build ran the stage before it or found it done (the hooks of a skipped stage are called in order, one after the other), and the state of a stage is
    /// recorded before the hook that follows the stage, so that a failure here is a kill just after the stage.
    /// </summary>
    internal Action<string>? Hook { get; init; }

    /// <summary>Test hook, called between the two passes over the input.</summary>
    internal Action? BetweenPasses { get; init; }

    /// <summary>Test hook passed to PH1 and PH2 / PH3 (<see cref="SaBuildOptions.Trace"/>).</summary>
    internal Action<string>? SaTrace { get; init; }

    /// <summary>Test hook passed to PH4 (<see cref="LcpChildOptions.Trace"/>).</summary>
    internal Action<string>? LcpTrace { get; init; }
}

/// <summary>Wall seconds of the phases of one build.</summary>
internal sealed record SesaPhaseSeconds
{
    internal double Scan { get; init; }

    internal double Create { get; init; }

    internal double Pack { get; init; }

    internal double Ph1 { get; init; }

    internal double Ph23 { get; init; }

    internal double Ph45 { get; init; }

    internal double Bwt { get; init; }

    internal double Psi { get; init; }

    internal double Crc { get; init; }

    internal double Seal { get; init; }

    internal double Total { get; init; }
}

/// <summary>Peak RssAnon of the whole process (MB, 10^6 bytes) during each phase, or null where it is unavailable.</summary>
internal sealed record SesaPhaseRss
{
    internal double? Ph0 { get; init; }

    internal double? Ph1 { get; init; }

    internal double? Ph23 { get; init; }

    internal double? Ph45 { get; init; }

    internal double? Ph6 { get; init; }

    internal double? Seal { get; init; }
}

/// <summary>The checksum pass over one section: bytes, seconds and whether the hardware path ran.</summary>
internal sealed record SesaCrcRow(int Kind, long Bytes, double Seconds, bool Hardware);

/// <summary>What a finished build knows.</summary>
internal sealed class SesaBuildResult
{
    internal required string TargetPath { get; init; }

    internal long N { get; init; }

    internal int Sigma { get; init; }

    internal int Bits { get; init; }

    internal long FileSize { get; init; }

    internal required SesaHeader Header { get; init; }

    /// <summary>The section table as written, in file order.</summary>
    internal required IReadOnlyList<SesaSection> Sections { get; init; }

    internal LcpKind LcpKind { get; init; }

    internal long OverflowCount { get; init; }

    internal long NodeCount { get; init; }

    internal long InternalCount { get; init; }

    internal long LrsDepth { get; init; }

    /// <summary>-1 when the text has no repeated substring.</summary>
    internal long LrsRank { get; init; }

    internal long MaxLcp { get; init; }

    internal bool HasPsi { get; init; }

    /// <summary>The row of the suffix 0 (<c>PARAM0</c> of PSI), -1 without PSI.</summary>
    internal long PrimaryRow { get; init; }

    internal required SaBuildResult Sa { get; init; }

    internal required LcpChildResult Lcp { get; init; }

    internal PsiStats? Psi { get; init; }

    internal BwtResult? Bwt { get; init; }

    /// <summary>The counters of the whole build.</summary>
    internal IoSnapshot Io { get; init; }

    internal required SesaPhaseSeconds Seconds { get; init; }

    internal required SesaPhaseRss Rss { get; init; }

    /// <summary>Peak RssAnon of the whole process during the build, MB, or null.</summary>
    internal double? RssAnonPeakMb { get; init; }

    internal RssVerdict RssVerdict { get; init; }

    /// <summary>The largest peak of the native buffers of any stage (the stages run one after the other).</summary>
    internal long ArenaPeakBytes { get; init; }

    internal long ScratchPeakBytes { get; init; }

    /// <summary>What the creation of the container did about the disk space.</summary>
    internal ReservationOutcome Reservation { get; init; }

    /// <summary>What the growth of the file behind the reserved length did (null when it did not grow): <see cref="ReservationOutcome.Filled"/> when its blocks were reserved by writing zeros before PH4 wrote the LCP (work package 1b.9d), the outcome of the check of the free space when the build was asked for sparse files.</summary>
    internal ReservationOutcome? GrowthReservation { get; init; }

    /// <summary>Bytes the file grew by after PH4-B2 (0 when it did not).</summary>
    internal long GrownBytes { get; init; }

    internal required IReadOnlyList<SesaCrcRow> CrcRows { get; init; }

    /// <summary>What <see cref="SesaBuildOptions.Resume"/> did (null when it was off).</summary>
    internal SesaResumeReport? Resume { get; init; }
}

/// <summary>Cost of recording one completed stage (<see cref="SesaBuildOptions.Resume"/>): the fsync of its outputs, the CRC pass over them, the write of the state file.</summary>
internal sealed record SesaCheckpointRow(string Stage, double FsyncSeconds, double CrcSeconds, double StateSeconds, long CrcBytes);

/// <summary>What a build with <see cref="SesaBuildOptions.Resume"/> found and recorded.</summary>
internal sealed class SesaResumeReport
{
    /// <summary>A state file was read (valid or not does not matter: it existed).</summary>
    internal bool StateFound { get; init; }

    /// <summary>A state of this input and layout was accepted, so that at least the stage <see cref="SesaStages.Packed"/> was skipped.</summary>
    internal bool Accepted { get; init; }

    /// <summary>In words: why the state was not used or why its stages were cut (empty when everything it held was used).</summary>
    internal string Reason { get; init; } = string.Empty;

    /// <summary>The stages that were not run because an earlier run completed them and their outputs verified, in order.</summary>
    internal IReadOnlyList<string> Skipped { get; init; } = [];

    /// <summary>Seconds spent finding, parsing and verifying the state (the CRC pass over the outputs of the skipped stages).</summary>
    internal double VerifySeconds { get; init; }

    /// <summary>The checkpoints this run recorded.</summary>
    internal IReadOnlyList<SesaCheckpointRow> Checkpoints { get; init; } = [];

    /// <summary>State files written and their bytes.</summary>
    internal int StateWrites { get; init; }

    internal long StateBytes { get; init; }
}

/// <summary>
/// The end-to-end builder of the SESA container (work package 1b.7e, ADR-004 sections 3-4, ADR-003 as amended on 2026-11-18): PH0 (the packed text, the alphabet, the symbol counts and the SHA-256) -> PH1-PH3 (the suffix
/// array) -> PH4-PH5 (LCP, child table) -> PH6 (BWT and PSI, optional) -> PH7 (the section checksums, the header, the atomic publication). Every big output is built directly into its section of
/// <c>&lt;target&gt;.tmp</c>, created through <see cref="MappedOutputFile.Create"/>; scratch lives under <c>&lt;target&gt;.tmp.d</c>. The target appears only by the rename at the end; a failure removes the temporary file and
/// the scratch and leaves no target (an existing target stays as it was). The output depends on the input and on <see cref="SesaBuildOptions.IncludePsi"/> only.
/// </summary>
internal static class SesaBuilder
{
    /// <summary>
    /// The codec of the CHILD section for a text of <paramref name="n"/> symbols (work package 1b.9a): <paramref name="requested"/> when it is given, else the tagged codec up to <see cref="TaggedChild.MaxLength"/> symbols (the file form
    /// of the layout that the in-memory tree uses, and the one whose descent touches fewer cache lines) and the packed one above it.
    /// </summary>
    /// <exception cref="NotSupportedException">The tagged codec is requested for a text whose ranks do not fit 30 bits.</exception>
    internal static ChildCodec ChildCodecFor(long n, ChildCodec? requested)
    {
        ChildCodec codec = requested ?? (n <= TaggedChild.MaxLength ? ChildCodec.Tagged : ChildCodec.Packed);
        if (codec == ChildCodec.Tagged && n > TaggedChild.MaxLength) throw new NotSupportedException($"the tagged child table holds ranks in 30 bits (CHILD codec 2); a text of {n} symbols needs the packed codec (CHILD codec 1)");
        return codec;
    }

    /// <summary>Builds the container of the text in the file <paramref name="inputPath"/> (one byte per code unit) at <paramref name="targetPath"/>.</summary>
    /// <exception cref="DirectoryNotFoundException">The directory of the target does not exist.</exception>
    /// <exception cref="NotSupportedException">The text is longer than <see cref="PackedTextFormat.MaxSymbols"/>.</exception>
    /// <exception cref="IOException">A file cannot be created or written (a full disk included); nothing of the build is left.</exception>
    /// <exception cref="InvalidOperationException">The input changed between the two passes.</exception>
    internal static SesaBuildResult Build(string inputPath, string targetPath, SesaBuildOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(inputPath);
        ArgumentNullException.ThrowIfNull(targetPath);
        options ??= new SesaBuildOptions();
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Threads, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ArenaCapBytes);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ChildArenaCapBytes);
        ArgumentOutOfRangeException.ThrowIfNegative(options.Ph6CapBytes);
        ArgumentOutOfRangeException.ThrowIfNegative(options.RssAnonCapMb);
        if (!BitConverter.IsLittleEndian) throw new PlatformNotSupportedException("the file layouts are little endian");
        string target = Path.GetFullPath(targetPath);
        string directory = Path.GetDirectoryName(target) ?? ".";
        if (!Directory.Exists(directory)) throw new DirectoryNotFoundException($"the directory of {targetPath} does not exist");
        string tmp = target + ".tmp";
        string scratch = tmp + ".d";
        using SesaBuildLock? gate = options.Resume ? SesaBuildLock.Acquire(target + ".lock") : null; // held until the build ends, whatever way
        IoCounters counters = options.Counters ?? new IoCounters();
        IReservationSource source = options.Files?.Source ?? SystemReservationSource.Instance;
        long tStart = Stopwatch.GetTimestamp();
        using var whole = options.MeasureRss ? new RssSampler(new RssReader(options.ProcSource), 10) : null;
        whole?.Start();
        void Hook(string point) => options.Hook?.Invoke(point);
        Hook(SesaHooks.Start);

        // ---- PH0: the alphabet, the counts, the hash (pass 1); what the layout needs
        using var ph0 = new PhaseRss(options);
        long t = Stopwatch.GetTimestamp();
        var textOptions = new PackedTextWriterOptions { BlockSymbols = options.TextBlockSymbols, WindowBytes = options.TextWindowBytes, Counters = counters, File = options.Files, Fsync = false };
        PackedTextInfo info = PackedTextWriter.Scan(inputPath, textOptions);
        double secScan = Elapsed(t);
        long n = info.N;
        Hook(SesaHooks.Scanned);
        options.BetweenPasses?.Invoke();
        ChildCodec childCodec = ChildCodecFor(n, options.ChildCodec);
        SesaPlan plan = SesaLayout.Plan(info, options.IncludePsi);
        SesaResumeRun? run = options.Resume ? SesaResumeRun.Load(options, tmp, scratch, info, plan, childCodec, counters) : null;
        bool resumed = run?.Accepted == true;
        long need = plan.ReservedLength + (ScratchBytesPerChar * n);
        if (resumed) need = ResumeNeedBytes(n, options.IncludePsi, run!.Skipped, run.ExistingScratchBytes);
        long free = source.AvailableFreeBytes(directory);
        if (resumed)
        {
            if (free >= 0 && free + run!.ReclaimableBytes < need) throw new IOException($"not enough free space in {directory}: the resumed build needs about {need} bytes for the stages that are not done ({ScratchBytesPerChar} bytes per character of scratch at most), {free} bytes are available and {run.ReclaimableBytes} more are freed first");
        }
        else if (free >= 0 && free < need)
        {
            throw new IOException($"not enough free space in {directory}: the build needs about {need} bytes (the container {plan.ReservedLength} and {ScratchBytesPerChar} bytes per character of scratch), {free} bytes are available");
        }

        if (resumed) run!.Prepare();
        else DeleteStale(tmp, scratch);
        bool published = false;
        MappedOutputFile? container = null;
        MappedTextReader? text = null;
        try
        {
            t = Stopwatch.GetTimestamp();
            Directory.CreateDirectory(scratch);
            if (resumed)
            {
                // the part behind the reserved length that an interrupted growth added is cut off, so that PH4 grows and reserves it again, whole, unless the LCP stage is done and the file is as it left it
                long? cut = run!.Has(SesaStages.Lcp) ? null : plan.ReservedLength;
                container = MappedOutputFile.OpenExisting(tmp, cut);
            }
            else
            {
                container = MappedOutputFile.Create(tmp, plan.ReservedLength, options.Files);
            }
            double secCreate = Elapsed(t);
            ReservationOutcome reservation = container.Outcome;
            Hook(SesaHooks.Created);

            t = Stopwatch.GetTimestamp();
            bool packedKept = run?.Has(SesaStages.Packed) == true;
            if (!packedKept)
            {
                using (MappedWindowWriter w = container.CreateWriter(plan.TextOffset, plan.TextLength, options.TextWindowBytes, counters)) PackedTextWriter.Pack(inputPath, info, w, textOptions);
                WriteSection(container, plan.AlphabetOffset, [.. info.Letters.Select(l => (ushort)l)], options.TextWindowBytes, counters);
                WriteSection(container, plan.CountsOffset, [.. info.Counts.Select(c => (ulong)c)], options.TextWindowBytes, counters);
            }
            double secPack = Elapsed(t);
            if (run is not null && !packedKept) run.CommitPacked(container);
            double? rssPh0 = ph0.Stop();
            Hook(SesaHooks.Packed);
            text = MappedTextReader.OpenRegion(tmp, plan.TextOffset, plan.TextLength, n, info.Bits);

            // ---- PH1-PH3: the suffix array, straight into its section
            var saOptions = new SaBuildOptions
            {
                ScratchDir = scratch,
                Threads = options.Threads,
                Period = options.Period,
                ExactPeriod = options.ExactPeriod,
                ScaleN = options.ScaleN,
                RuleN = options.RuleN,
                ArenaCapBytes = options.ArenaCapBytes,
                Tuning = options.SaTuning,
                Seed = options.Seed,
                Counters = counters,
                Files = options.Files,
                ArenaGuards = options.ArenaGuards,
                VerifyRanks = options.VerifyRanks,
                MeasureRss = options.MeasureRss,
                ProcSource = options.ProcSource,
                RssAnonCapMb = options.RssAnonCapMb,
                Trace = options.SaTrace,
                Checkpoint = run is null ? null : new SaCheckpoint(run, container),
            };
            SaBuildResult sa = run?.Has(SesaStages.Sa) == true ? SkippedSa(n, saOptions, container.Outcome) : SaBuilder.Build(text, container, plan.SaOffset, saOptions);
            Hook(SesaHooks.SaBuilt);

            // ---- PH4-PH5: LCP (placed behind the sections known after PH0, once its width is known), overflow list and child table, straight into their sections
            SesaTail? tail = null;
            ReservationOutcome? growth = null;
            long grown = 0;
            LcpPlacement Place(LcpPlacementRequest request)
            {
                tail = SesaLayout.PlaceTail(plan, request.Kind, request.OverflowCount);
                long growthBytes = tail.FinalLength - container.Length;
                long spillBound = SpillBytesPerEntry * (request.MaxLcp + 1);
                long available = source.AvailableFreeBytes(directory);
                if (available >= 0 && available < growthBytes + spillBound) throw new IOException($"not enough free space in {directory} for the LCP sections and the spill of the scan stack: {growthBytes} + {spillBound} bytes needed, {available} bytes are available");
                if (growthBytes > 0)
                {
                    growth = container.Grow(tail.FinalLength, options.Files, counters);
                    grown = growthBytes;
                }
                Hook(SesaHooks.LcpPlaced);
                return new LcpPlacement(container, plan.LcpOffset, tail.OverflowOffset);
            }
            var lcpOptions = new LcpChildOptions
            {
                ScratchDir = scratch,
                Threads = options.Threads,
                Codec = childCodec,
                ScaleN = options.ScaleN,
                RuleN = options.RuleN,
                ArenaCapBytes = options.ArenaCapBytes,
                ChildArenaCapBytes = options.ChildArenaCapBytes,
                Tuning = options.LcpTuning,
                Counters = counters,
                Files = options.Files,
                ArenaGuards = options.ArenaGuards,
                MeasureRss = options.MeasureRss,
                ProcSource = options.ProcSource,
                RssAnonCapMb = options.RssAnonCapMb,
                Trace = options.LcpTrace,
            };
            LcpChildResult lcp;
            bool lcpKept = run?.Has(SesaStages.Lcp) == true;
            if (lcpKept)
            {
                SesaStageRecord kept = run!.Record(SesaStages.Lcp);
                tail = SesaLayout.PlaceTail(plan, (LcpKind)kept.Value("kind"), kept.Value("overflow"));
                lcp = SkippedLcp(n, lcpOptions, tail, kept, container.Outcome);
                Hook(SesaHooks.LcpPlaced);
            }
            else
            {
                lcp = SectionLcpChildBuilder.Build(tmp, plan.SaOffset, text, container, plan.ChildOffset, Place, lcpOptions);
            }
            SesaTail placed = tail ?? throw new InvalidOperationException("PH4 did not ask for the placement of the LCP");
            if (run is not null && !lcpKept) run.CommitLcp(container, placed, lcp.Scan);
            Hook(SesaHooks.LcpChildBuilt);

            // ---- PH6 (optional): the BWT as scratch, then PSI by the sigma-way distribution, straight into its section
            BwtResult? bwt = null;
            PsiStats? psi = null;
            double secBwt = 0, secPsi = 0;
            double? rssPh6 = null;
            long primaryRow = -1;
            if (options.IncludePsi)
            {
                if (run?.Has(SesaStages.Psi) == true)
                {
                    primaryRow = run.PrimaryRow();
                    Hook(SesaHooks.BwtBuilt);
                    Hook(SesaHooks.PsiBuilt);
                }
                else
                {
                    using var ph6 = new PhaseRss(options);
                    (bwt, psi, secBwt, secPsi) = RunPh6(options, counters, scratch, tmp, plan, info, text, container, Hook, run);
                    rssPh6 = ph6.Stop();
                    primaryRow = bwt.PrimaryRow;
                }
            }

            // ---- PH7: BUILD_INFO, the checksums, the header page, the atomic publication
            using var seal = new PhaseRss(options);
            t = Stopwatch.GetTimestamp();
            string buildInfo = BuildInfoText(info, placed.Kind, placed.OverflowCount, options.IncludePsi, childCodec);
            byte[] infoBytes = Encoding.UTF8.GetBytes(buildInfo);
            CheckBuildInfoFits(infoBytes.Length, plan.BuildInfoCapacity);
            WriteSection(container, plan.BuildInfoOffset, infoBytes, options.TextWindowBytes, counters);
            Hook(SesaHooks.InfoWritten);
            container.Complete(fsync: false);
            text.Dispose();
            text = null;
            container.Dispose();
            CheckContainerLength(new FileInfo(tmp).Length, placed.FinalLength);

            (int Kind, long Offset, long Length)[] ranges = SectionRanges(plan, placed, infoBytes.Length);
            (uint Crc, SesaCrcRow Row)[] sums = SesaChecksums.Compute(tmp, ranges);
            double secCrc = Elapsed(t);
            Hook(SesaHooks.Checksummed);

            t = Stopwatch.GetTimestamp();
            SesaSection[] sections = BuildSections(plan, placed, info, primaryRow, buildInfo, infoBytes.Length, sums.Select(s => s.Crc).ToArray(), ranges, childCodec);
            var header = new SesaHeader
            {
                TextLength = n,
                NodeCount = lcp.Scan.NodeCount,
                FileSize = placed.FinalLength,
                LrsDepth = lcp.Scan.LrsDepth,
                LrsRank = lcp.Scan.LrsDepth == 0 ? -1 : lcp.Scan.LrsRank,
                PosWidth = 4,
                Sigma = info.Sigma,
                TextSha256 = Convert.FromHexString(info.Sha256),
            };
            byte[] page = SesaFormat.EncodeHeaderPage(header, sections, withMagic: false);
            using (var file = new CountedFile(tmp, FileMode.Open, FileAccess.ReadWrite, counters))
            {
                file.WriteAt(page, 0);
                Hook(SesaHooks.HeaderWritten);
                if (options.Fsync) file.Flush();
                Hook(SesaHooks.HeaderSynced);
                file.WriteAt(SesaFormat.Magic, 0); // the magic is written last: a file without it is not an index
                if (options.Fsync) file.Flush();
                Hook(SesaHooks.MagicWritten);
            }
            File.Move(tmp, target, overwrite: true);
            published = true;
            DeleteDirectoryQuietly(scratch);
            double secSeal = Elapsed(t);
            double? rssSeal = seal.Stop();

            RssStats all = whole?.Stop() ?? default;
            double? anonMb = all.Available ? all.MaxRssAnonKb * 1024.0 / 1e6 : null;
            RssVerdict verdict = RssVerdict.Unavailable;
            if (anonMb is { } mb) verdict = mb <= options.RssAnonCapMb ? RssVerdict.Within : RssVerdict.Exceeded;
            return new SesaBuildResult
            {
                TargetPath = target,
                N = n,
                Sigma = info.Sigma,
                Bits = info.Bits,
                FileSize = placed.FinalLength,
                Header = header,
                Sections = sections,
                LcpKind = placed.Kind,
                OverflowCount = placed.OverflowCount,
                NodeCount = lcp.Scan.NodeCount,
                InternalCount = lcp.Scan.InternalCount,
                LrsDepth = lcp.Scan.LrsDepth,
                LrsRank = header.LrsRank,
                MaxLcp = lcp.Scan.MaxLcp,
                HasPsi = options.IncludePsi,
                PrimaryRow = primaryRow,
                Sa = sa,
                Lcp = lcp,
                Psi = psi,
                Bwt = bwt,
                Io = counters.Snapshot(),
                Seconds = new SesaPhaseSeconds
                {
                    Scan = secScan,
                    Create = secCreate,
                    Pack = secPack,
                    Ph1 = sa.Ph1Seconds,
                    Ph23 = sa.Ph23Seconds,
                    Ph45 = lcp.Seconds,
                    Bwt = secBwt,
                    Psi = secPsi,
                    Crc = secCrc,
                    Seal = secSeal,
                    Total = Elapsed(tStart),
                },
                Rss = new SesaPhaseRss { Ph0 = rssPh0, Ph1 = sa.Ph1RssAnonPeakMb, Ph23 = sa.Ph23RssAnonPeakMb, Ph45 = lcp.RssAnonPeakMb, Ph6 = rssPh6, Seal = rssSeal },
                RssAnonPeakMb = anonMb,
                RssVerdict = verdict,
                ArenaPeakBytes = Math.Max(Math.Max(sa.Ph1ArenaPeakBytes, sa.Ph23ArenaPeakBytes), Math.Max(lcp.ArenaPeakBytes + lcp.ChildArenaPeakBytes, psi?.MappedBoundBytes ?? 0)),
                ScratchPeakBytes = counters.ScratchPeakBytes,
                Reservation = reservation,
                GrowthReservation = growth,
                GrownBytes = grown,
                CrcRows = sums.Select(s => s.Row).ToArray(),
                Resume = run?.Report(),
            };
        }
        finally
        {
            text?.Dispose();
            container?.Dispose();
            bool keep = !published && run is { Keeps: true }; // a build that records its stages leaves them for the next call when it fails
            if (!published && !keep) DeleteQuietly(tmp);
            if (!keep) DeleteDirectoryQuietly(scratch);
        }
    }

    /// <summary>Bytes per character of scratch that the pre-flight reserves: the peak of PH4 in base mode (E4: 12.0 B/char) plus the rank file.</summary>
    internal const int ScratchBytesPerChar = 12;

    /// <summary>The bytes of scratch that the stages still to run will need (the pre-flight of a resume), from the stages that are kept and the bytes of scratch files that already exist.</summary>
    internal static long ResumeNeedBytes(long n, bool psi, IReadOnlyCollection<string> kept, long existingScratchBytes)
    {
        long need = ScratchBytesPerChar * n;
        if (kept.Contains(SesaStages.Lcp))
        {
            // only the BWT file is still to be created; PSI is written into the container, which is reserved
            bool bwtToCome = psi && !kept.Contains(SesaStages.Psi);
            need = bwtToCome ? n + 1 : 0;
        }
        return Math.Max(0, need - existingScratchBytes);
    }

    /// <summary>What a build that found PH1-PH3 done reports for them: the sizes that its options would give and nothing that ran.</summary>
    private static SaBuildResult SkippedSa(long n, SaBuildOptions options, ReservationOutcome outcome) => new()
    {
        Parameters = SaBuildRules.Resolve(n, options),
        N = n,
        Dcs = new ExtDcsStats(),
        Sort = new SortStats(),
        Reservation = outcome,
    };

    /// <summary>What a build that found PH4-PH5 done reports for them: the numbers that the stage recorded.</summary>
    private static LcpChildResult SkippedLcp(long n, LcpChildOptions options, SesaTail tail, SesaStageRecord record, ReservationOutcome outcome)
    {
        long depth = record.Value("lrs_depth");
        return new LcpChildResult
        {
            Parameters = LcpChildRules.Resolve(n, options),
            Lcp = new ExtLcpStats { Kind = tail.Kind, N = n, OverflowCount = tail.OverflowCount, Reservation = outcome },
            Scan = new ChildScanResult
            {
                N = n,
                InternalCount = record.Value("internal_count"),
                NodeCount = record.Value("node_count"),
                LrsDepth = (uint)depth,
                LrsRank = depth == 0 ? 0 : (uint)record.Value("lrs_rank"),
                MaxLcp = record.Value("max_lcp"),
            },
            Reservation = outcome,
        };
    }

    /// <summary>The window of <see cref="SaBuilder"/> onto the recording of the build: the rank file after PH1, the suffix array after PH3.</summary>
    private sealed class SaCheckpoint(SesaResumeRun run, MappedOutputFile container) : ISaCheckpoint
    {
        public bool KeepRanks => run.Vouches("ranks.bin");

        public double CheckpointSeconds => run.CheckpointSeconds;

        public bool TryReuseRanks(DifferenceCover cover) => run.Has(SesaStages.Ranks);

        public void RanksBuilt(DifferenceCover cover) => run.CommitRanks(cover);

        public void SaWritten() => run.CommitSa(container);
    }

    /// <summary>Bytes of one entry of the scan stack of PH5 (the upper bound of its spill is this times the depth).</summary>
    internal const int SpillBytesPerEntry = 16;

    private static (BwtResult Bwt, PsiStats Psi, double SecBwt, double SecPsi) RunPh6(SesaBuildOptions options, IoCounters counters, string scratch, string tmp, SesaPlan plan, PackedTextInfo info, MappedTextReader text, MappedOutputFile container, Action<string> hook, SesaResumeRun? run)
    {
        long n = info.N;
        string dir = Path.Combine(scratch, "ph6");
        Directory.CreateDirectory(dir);
        string bwtPath = Path.Combine(dir, "bwt.bin");
        var budget = new ArenaBudget(options.Ph6CapBytes);
        SesaPh6Tuning tune = options.Ph6Tuning;
        int window = tune.WindowBytes ?? MappedWindowWriter.DefaultWindowBytes;
        long t = Stopwatch.GetTimestamp();
        bool bwtKept = run?.Has(SesaStages.Bwt) == true;
        BwtResult bwt = bwtKept
            ? new BwtResult(n, run!.PrimaryRow(), [.. info.Counts]) // an earlier run built it, its counts were those of the text, and the file verified
            : BwtPass.Build(tmp, plan.SaOffset, text, info.Sigma, bwtPath, new BwtOptions
            {
                ReadBytes = tune.ReadBytes ?? 1 << 20,
                WindowBytes = window,
                Budget = budget,
                ArenaGuards = options.ArenaGuards,
                Counters = counters,
                Files = options.Files,
            });
        double secBwt = bwtKept ? 0 : Elapsed(t);
        counters.AddScratch(n + 1); // the BWT file is scratch until PSI is built
        try
        {
            if (!bwt.Counts.AsSpan().SequenceEqual(info.Counts)) throw new InvalidDataException("the symbol counts of the BWT differ from those of the text: the suffix array does not belong to the text");
            if (run is not null && !bwtKept) run.CommitBwt(bwt.PrimaryRow);
            hook(SesaHooks.BwtBuilt);
            t = Stopwatch.GetTimestamp();
            PsiStats psi = PsiStream.Build(bwtPath, n, bwt.PrimaryRow, info.Counts, container, plan.PsiOffset, new PsiOptions
            {
                WindowBytes = window,
                SymbolsPerPass = tune.SymbolsPerPass ?? 64,
                ReadBytes = tune.ReadBytes ?? 1 << 20,
                WindowCapBytes = options.Ph6CapBytes,
                Budget = budget,
                ArenaGuards = options.ArenaGuards,
                Counters = counters,
            });
            double secPsi = Elapsed(t);
            run?.CommitPsi(container, bwt.PrimaryRow);
            File.Delete(bwtPath);
            Directory.Delete(dir);
            hook(SesaHooks.PsiBuilt);
            return (bwt, psi, secBwt, secPsi);
        }
        finally
        {
            counters.RemoveScratch(n + 1);
            if (run is null || !run.Vouches("ph6/bwt.bin")) DeleteQuietly(bwtPath); // the state vouches for it: a failure of PSI must not lose the BWT
        }
    }

    private static void WriteSection<T>(MappedOutputFile container, long offset, T[] values, int window, IoCounters counters)
        where T : unmanaged
        => WriteSection(container, offset, MemoryMarshal.AsBytes(values.AsSpan()), window, counters);

    private static void WriteSection(MappedOutputFile container, long offset, ReadOnlySpan<byte> bytes, int window, IoCounters counters)
    {
        using MappedWindowWriter w = container.CreateWriter(offset, bytes.Length, window, counters);
        w.Append(bytes);
    }

    /// <summary>The payload ranges of the sections in file order (the order of the table).</summary>
    private static (int Kind, long Offset, long Length)[] SectionRanges(SesaPlan plan, SesaTail tail, int infoLength)
    {
        var list = new List<(int, long, long)>
        {
            (SesaFormat.KindText, plan.TextOffset, plan.TextLength),
            (SesaFormat.KindAlphabet, plan.AlphabetOffset, plan.AlphabetLength),
            (SesaFormat.KindCounts, plan.CountsOffset, plan.CountsLength),
            (SesaFormat.KindSa, plan.SaOffset, plan.SaLength),
            (SesaFormat.KindChild, plan.ChildOffset, plan.ChildLength),
        };
        if (plan.Psi) list.Add((SesaFormat.KindPsi, plan.PsiOffset, plan.PsiLength));
        list.Add((SesaFormat.KindBuildInfo, plan.BuildInfoOffset, infoLength));
        list.Add((SesaFormat.KindLcp, plan.LcpOffset, tail.LcpLength));
        if (tail.OverflowCount > 0) list.Add((SesaFormat.KindLcpOverflow, tail.OverflowOffset, tail.OverflowLength));
        return [.. list];
    }

    private static SesaSection[] BuildSections(SesaPlan plan, SesaTail tail, PackedTextInfo info, long primaryRow, string buildInfo, int infoLength, uint[] crcs, (int Kind, long Offset, long Length)[] ranges, ChildCodec childCodec)
    {
        long rows = plan.N + 1;
        const int required = SesaFormat.FlagRequired | SesaFormat.FlagHasCrc;
        const int optional = SesaFormat.FlagHasCrc;
        var sections = new SesaSection[ranges.Length];
        for (int i = 0; i < ranges.Length; i++)
        {
            (int kind, long offset, long length) = ranges[i];
            sections[i] = kind switch
            {
                SesaFormat.KindText => new SesaSection(kind, required, offset, length, plan.N, 0, PackedTextFormat.CodecId, info.Bits, 0, crcs[i]),
                SesaFormat.KindAlphabet => new SesaSection(kind, required, offset, length, info.Sigma, 2, 1, 0, 0, crcs[i]),
                SesaFormat.KindCounts => new SesaSection(kind, optional, offset, length, info.Sigma, 8, 1, 0, 0, crcs[i]),
                SesaFormat.KindSa => new SesaSection(kind, required, offset, length, rows, 4, 1, 0, 0, crcs[i]),
                SesaFormat.KindChild => new SesaSection(kind, required, offset, length, rows, 4, childCodec == ChildCodec.Tagged ? SesaFormat.ChildCodecTagged : SesaFormat.ChildCodecPacked, 0, 0, crcs[i]),
                SesaFormat.KindPsi => new SesaSection(kind, optional, offset, length, rows, 4, 1, primaryRow, 0, crcs[i]),
                SesaFormat.KindBuildInfo => new SesaSection(kind, optional, offset, length, buildInfo.Count(c => c == '\n'), 0, 1, 0, 0, crcs[i]),
                SesaFormat.KindLcp => new SesaSection(kind, required, offset, length, rows, tail.LcpWidth, tail.Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 3 }, tail.OverflowCount, 0, crcs[i]),
                SesaFormat.KindLcpOverflow => new SesaSection(kind, required, offset, length, tail.OverflowCount, 8, 1, 0, 0, crcs[i]),
                _ => throw new InvalidOperationException($"unknown section kind {kind}"),
            };
        }
        _ = infoLength;
        return sections;
    }

    /// <summary>BUILD_INFO is written into a room of fixed size that the layout reserved before it was known (4,096 bytes); a text that would overrun it would overwrite the next section, so it is refused.</summary>
    /// <exception cref="InvalidOperationException">The text of <paramref name="length"/> bytes does not fit <paramref name="capacity"/>.</exception>
    internal static void CheckBuildInfoFits(int length, long capacity)
    {
        if (length > capacity) throw new InvalidOperationException($"BUILD_INFO has {length} bytes, the reserved room is {capacity}");
    }

    /// <summary>The container, closed, has exactly the length that the layout planned (a file that is shorter or longer is not the container that the table describes).</summary>
    internal static void CheckContainerLength(long actual, long planned)
    {
        if (actual != planned) throw new InvalidOperationException($"the container has {actual} bytes, expected {planned}");
    }

    /// <summary>The text of the BUILD_INFO section: everything that is a fact about the output, nothing about the run (no time, no thread count, no sizes), so that the same input gives the same bytes.</summary>
    internal static string BuildInfoText(PackedTextInfo info, LcpKind kind, long overflowCount, bool psi, ChildCodec childCodec)
        => string.Create(CultureInfo.InvariantCulture,
            $"format=SESA 1.0\nwriter=SuffixTree.Esa.External.SesaBuilder 1\nalgorithm=B-DCS (ADR-004)\nn={info.N}\nsigma={info.Sigma}\ntext=codec 4, {info.Bits} bits per symbol\nlcp={kind}, {overflowCount} overflow entries\nchild=layout A, {(childCodec == ChildCodec.Tagged ? "tagged" : "packed")}\npsi={(psi ? "yes" : "no")}\n");

    private static double Elapsed(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;

    /// <summary>What a hard kill of an earlier build of this target left: the temporary file and the scratch directory (a file of that name too) are replaced.</summary>
    private static void DeleteStale(string tmp, string scratch)
    {
        if (File.Exists(tmp)) File.Delete(tmp);
        if (Directory.Exists(scratch)) Directory.Delete(scratch, true);
        else if (File.Exists(scratch)) File.Delete(scratch);
    }

    private static void DeleteQuietly(string path)
    {
        try
        {
            File.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // the failure of the build is what the caller needs
        }
    }

    private static void DeleteDirectoryQuietly(string path)
    {
        try
        {
            if (Directory.Exists(path)) Directory.Delete(path, true);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // reported by the files that are left
        }
    }

    /// <summary>One RssAnon sampler per phase (a no-op where measuring is off or the platform has no <c>/proc</c>).</summary>
    private sealed class PhaseRss : IDisposable
    {
        private readonly RssSampler? _sampler;

        internal PhaseRss(SesaBuildOptions options)
        {
            if (!options.MeasureRss) return;
            _sampler = new RssSampler(new RssReader(options.ProcSource), 10);
            _sampler.Start();
        }

        /// <summary>The peak RssAnon of the process since the start, MB (10^6 bytes), or null.</summary>
        internal double? Stop()
        {
            RssStats stats = _sampler?.Stop() ?? default;
            return stats.Available ? stats.MaxRssAnonKb * 1024.0 / 1e6 : null;
        }

        public void Dispose() => _sampler?.Dispose();
    }
}
