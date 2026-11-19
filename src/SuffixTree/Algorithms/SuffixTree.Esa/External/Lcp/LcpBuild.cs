using System.Diagnostics;
using System.Runtime.ExceptionServices;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>
/// One run of PH4 in its base mode (ADR-004 section 4.2, row PH4). B1 scans the suffix array file in rank order and distributes the tuples <c>(SA[i], SA[i-1], i)</c> by <c>SA[i] / R</c> into P chunk files; B2 takes the
/// chunks (several in parallel), scatters Phi and the rank into two arrays of R entries (and checks that every position of the chunk occurs exactly once: the suffix array is a permutation), runs the PLCP loop
/// (<see cref="PlcpScan"/>), counts the escapes (so that the width is known before any LCP entry is written) and distributes <c>(i, PLCP[j])</c> by <c>i / R</c>; B3 takes the rank chunks in order, scatters them into a dense
/// array, hands it to the chunk consumer (the stream point of PH5) and writes the compact form through mapped windows. The output does not depend on the thread count, R or the buffer sizes.
/// </summary>
internal sealed unsafe class LcpBuild
{
    private readonly string _saPath;
    private readonly MappedTextReader _text;
    private readonly string _basePath;
    private readonly string _overflowPath;
    private readonly ExtLcpOptions _o;
    private readonly ILcpChunkConsumer? _consumer;
    private readonly LcpGeometry _g;
    private readonly ArenaBudget _budget;
    private readonly ExtLcpStats _stats = new();
    private readonly object _gate = new();
    private Exception? _error;
    private long _nextChunk;
    private long _escapes8;
    private long _escapes16;
    private long _steps;
    private long _maxPlcp;
    private long _restart;
    private bool _baseCreated;
    private bool _overflowCreated;

    internal LcpBuild(string saPath, MappedTextReader text, string basePath, string overflowPath, ExtLcpOptions options, ILcpChunkConsumer? consumer)
    {
        ArgumentNullException.ThrowIfNull(saPath);
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(basePath);
        ArgumentNullException.ThrowIfNull(overflowPath);
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ChunkPositions, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Threads, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(options.Workers);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ReadBufferBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.TupleBufferBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(options.SaOffset);
        if (!BitConverter.IsLittleEndian) throw new PlatformNotSupportedException("the file layouts are little endian");
        _saPath = saPath;
        _text = text;
        _basePath = basePath;
        _overflowPath = overflowPath;
        _o = options;
        _consumer = consumer;
        _g = new LcpGeometry(text.Length, options.ChunkPositions);
        _budget = options.Budget ?? new ArenaBudget(long.MaxValue);
    }

    /// <summary>Bytes of native buffers one B2 task holds (see <see cref="ExternalLcp.WorkerBytes"/>).</summary>
    internal static long WorkerBytes(long n, ExtLcpOptions o)
    {
        ArgumentNullException.ThrowIfNull(o);
        var g = new LcpGeometry(n, o.ChunkPositions);
        long chunk = Math.Min(o.ChunkPositions, Math.Max(n, 1));
        long seen = (chunk + 63) / 64 * sizeof(ulong);
        long tuples = Math.Max(1, o.TupleBufferBytes / sizeof(RankTuple)) * sizeof(RankTuple) * g.RankChunks;
        return (8 * chunk) + seen + (Math.Max(1, o.ReadBufferBytes / sizeof(PhiTuple)) * (long)sizeof(PhiTuple)) + tuples;
    }

    internal ExtLcpStats Run()
    {
        long n = _text.Length;
        CheckSuffixArrayFile(n);
        Directory.CreateDirectory(_o.ScratchDir);
        _stats.N = n;
        _stats.TextChunks = _g.TextChunks;
        _stats.RankChunks = _g.RankChunks;
        try
        {
            RunDistributed();
        }
        catch
        {
            if (_baseCreated) DeleteQuietly(_basePath);
            if (_overflowCreated) DeleteQuietly(_overflowPath);
            throw;
        }
        _stats.Escapes8 = _escapes8;
        _stats.Escapes16 = _escapes16;
        _stats.CompareSteps = _steps;
        _stats.MaxPlcp = _maxPlcp;
        _stats.RestartSymbols = _restart;
        return _stats;
    }

    private static void DeleteQuietly(string path)
    {
        try
        {
            File.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // the path is not ours to clean up if it cannot be deleted; the failure of the build is what the caller needs
        }
    }

    private int ReadEntries(int tupleBytes) => Math.Max(1, _o.ReadBufferBytes / tupleBytes);

    private int ChunkEntries(long chunkLength) => checked((int)Math.Min(chunkLength, _g.R));

    /// <summary>The file must hold <c>n + 1</c> entries and start with <c>SA[0] = n</c> (the empty suffix).</summary>
    private void CheckSuffixArrayFile(long n)
    {
        long length = new FileInfo(_saPath).Length;
        long need = (n + 1) * 4;
        bool fits = _o.SaOffset > 0 ? length >= _o.SaOffset + need : length == need; // a section of a larger file, or a file of its own that is exactly the array
        if (!fits) throw new InvalidDataException($"the suffix array needs n + 1 = {n + 1} entries of 4 bytes at offset {_o.SaOffset}, the file has {length} bytes");
        Span<uint> first = stackalloc uint[1];
        using var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read, _o.Counters);
        sa.ReadExactAt(MemoryMarshal.AsBytes(first), _o.SaOffset);
        if (first[0] != n) throw new InvalidDataException($"sa[0] must be n = {n}, the empty suffix, it is {first[0]}");
    }

    private void RunDistributed()
    {
        long n = _text.Length;
        string dir = _o.ScratchDir;
        using var phiSpool = new DistScratch<PhiTuple>(dir, "ph4phi", _g.TextChunks, _o.Counters);
        using var rankSpool = new DistScratch<RankTuple>(dir, "ph4rank", _g.RankChunks, _o.Counters);

        long t0 = Stopwatch.GetTimestamp();
        _o.Trace?.Invoke("b1");
        if (n > 0) DistributePhi(phiSpool);
        long t1 = Stopwatch.GetTimestamp();
        _o.Trace?.Invoke("b2");
        if (n > 0) RunWorkers(WorkerBytes(n, _o), _g.TextChunks, () => PlcpWorker(phiSpool, rankSpool));
        long t2 = Stopwatch.GetTimestamp();
        ChooseKind(n);
        _o.Trace?.Invoke("b3");
        EmitRankChunks(rankSpool);
        long t3 = Stopwatch.GetTimestamp();
        _stats.B1Seconds = (t1 - t0) / (double)Stopwatch.Frequency;
        _stats.B2Seconds = (t2 - t1) / (double)Stopwatch.Frequency;
        _stats.B3Seconds = (t3 - t2) / (double)Stopwatch.Frequency;
    }

    /// <summary>B1: the suffix array in rank order, one tuple per entry into the chunk file of its position.</summary>
    private void DistributePhi(DistScratch<PhiTuple> spool)
    {
        long n = _text.Length;
        using var arena = new BudgetedArena(_budget, _o.ArenaGuards);
        int readEntries = ReadEntries(sizeof(uint));
        Buf<uint> read = arena.Allocate<uint>(readEntries, zero: false, "b1.sa");
        var writer = new DistWriter<PhiTuple>(spool, arena, 0, _g.TextChunks, _o.TupleBufferBytes);
        using var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read, _o.Counters);
        uint previous = (uint)n; // SA[0], the empty suffix
        for (long first = 1; first <= n; first += readEntries)
        {
            int take = (int)Math.Min(readEntries, n + 1 - first);
            sa.ReadExactAt(MemoryMarshal.AsBytes(read.AsSpan(take)), _o.SaOffset + (first * 4));
            for (int i = 0; i < take; i++)
            {
                uint s = read[i];
                if (s >= n) throw new InvalidDataException($"sa[{first + i}] = {s} is outside [0, {n})");
                writer.Add(_g.ChunkOf(s), new PhiTuple(s, previous, (uint)(first + i)));
                previous = s;
            }
        }
        writer.FlushAll();
    }

    /// <summary>B2 of one task: chunks in any order; Phi and rank scattered into arrays of R entries, every position seen exactly once, PLCP in place over Phi, the second distribution by rank.</summary>
    private void PlcpWorker(DistScratch<PhiTuple> phiSpool, DistScratch<RankTuple> rankSpool)
    {
        using var arena = new BudgetedArena(_budget, _o.ArenaGuards);
        int maxChunk = ChunkEntries(_g.TextChunkLength(0));
        Buf<uint> phi = arena.Allocate<uint>(maxChunk, zero: false, "b2.phi");
        Buf<uint> rank = arena.Allocate<uint>(maxChunk, zero: false, "b2.rank");
        Buf<ulong> seen = arena.Allocate<ulong>((maxChunk + 63) / 64, zero: false, "b2.seen");
        int readEntries = ReadEntries(sizeof(PhiTuple));
        Buf<PhiTuple> read = arena.Allocate<PhiTuple>(readEntries, zero: false, "b2.read");
        var writer = new DistWriter<RankTuple>(rankSpool, arena, 0, _g.RankChunks, _o.TupleBufferBytes);
        var total = new PlcpTotals();
        while (Volatile.Read(ref _error) is null)
        {
            long c = Interlocked.Increment(ref _nextChunk) - 1;
            if (c >= _g.TextChunks) break;
            int length = (int)_g.TextChunkLength(c);
            long start = c * _g.R;
            seen.AsSpan((length + 63) / 64).Clear();
            using (DistReader<PhiTuple> reader = phiSpool.OpenReader(c))
            {
                if (reader.Count != length) throw new InvalidDataException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");
                int got;
                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)
                {
                    for (int i = 0; i < got; i++)
                    {
                        PhiTuple t = read[i];
                        long at = t.Position - start;
                        if ((ulong)at >= (ulong)length) throw new InvalidDataException($"position {t.Position} does not belong to text chunk {c}");
                        ref ulong word = ref seen[(int)(at >> 6)];
                        ulong bit = 1UL << (int)(at & 63);
                        if ((word & bit) != 0) throw new InvalidDataException($"position {t.Position} occurs twice in the suffix array");
                        word |= bit;
                        phi[(int)at] = t.Phi;
                        rank[(int)at] = t.Rank;
                    }
                }
            }
            phiSpool.Delete(c); // the chunk is in memory now: its file is scratch no more
            total.Add(PlcpScan.Run(_text, start, phi, length));
            for (int k = 0; k < length; k++)
            {
                uint r = rank[k];
                writer.Add(_g.ChunkOf(r), new RankTuple(r, phi[k]));
            }
        }
        if (Volatile.Read(ref _error) is null) writer.FlushAll();
        Merge(total);
    }

    /// <summary>B3: rank chunks in ascending order, scattered into a dense array, handed to the consumer and written in the chosen width.</summary>
    private void EmitRankChunks(DistScratch<RankTuple> rankSpool)
    {
        long n = _text.Length;
        using var arena = new BudgetedArena(_budget, _o.ArenaGuards);
        int maxChunk = ChunkEntries(_g.RankChunkLength(0));
        Buf<uint> dense = arena.Allocate<uint>(maxChunk, zero: false, "b3.dense");
        int readEntries = ReadEntries(sizeof(RankTuple));
        Buf<RankTuple> read = arena.Allocate<RankTuple>(readEntries, zero: false, "b3.read");
        long baseBytes = (n + 1) * LcpWidth.Bytes(_stats.Kind);
        long overflowBytes = _stats.OverflowCount * 8;
        MappedOutputFile? baseFile = null, overflowFile = null;
        try
        {
            MappedWindowWriter baseWriter, overflowWriter;
            if (_o.Place is { } place)
            {
                // the container: the sections are placed (and the file grown) now that the width is known; the file is the owner's, so nothing is completed, disposed or removed here
                LcpPlacement placement = place(new LcpPlacementRequest(_stats.Kind, _stats.OverflowCount, _maxPlcp));
                baseWriter = placement.File.CreateWriter(placement.BaseStart, baseBytes, _o.WindowBytes, _o.Counters);
                try
                {
                    overflowWriter = placement.File.CreateWriter(placement.OverflowStart, overflowBytes, _o.WindowBytes, _o.Counters);
                }
                catch
                {
                    baseWriter.Dispose();
                    throw;
                }
                _stats.Reservation = placement.File.Outcome;
            }
            else
            {
                baseFile = MappedOutputFile.Create(_basePath, baseBytes, _o.Files);
                _baseCreated = true;
                _stats.Reservation = baseFile.Outcome;
                overflowFile = MappedOutputFile.Create(_overflowPath, overflowBytes, _o.Files);
                _overflowCreated = true;
                baseWriter = baseFile.CreateWriter(_o.WindowBytes, _o.Counters);
                overflowWriter = overflowFile.CreateWriter(_o.WindowBytes, _o.Counters);
            }
            LcpOutputWriter output;
            try
            {
                output = new LcpOutputWriter(arena, _stats.Kind, n, _stats.OverflowCount, baseWriter, overflowWriter, ReadEntries(sizeof(uint)));
            }
            catch
            {
                baseWriter.Dispose();
                overflowWriter.Dispose();
                throw;
            }
            using (output)
            {
                for (long c = 0; c < _g.RankChunks; c++)
                {
                    int length = (int)_g.RankChunkLength(c);
                    long start = c * _g.R;
                    long expected = c == 0 ? length - 1 : length; // rank 0 is the empty suffix: no tuple
                    using (DistReader<RankTuple> reader = rankSpool.OpenReader(c))
                    {
                        if (reader.Count != expected) throw new InvalidDataException($"rank chunk {c} holds {reader.Count} tuples, expected {expected}: the suffix array is not a permutation");
                        int got;
                        while ((got = reader.Read(read.AsSpan(readEntries))) > 0)
                        {
                            for (int i = 0; i < got; i++)
                            {
                                RankTuple t = read[i];
                                long at = t.Rank - start;
                                if ((ulong)at >= (ulong)length) throw new InvalidDataException($"rank {t.Rank} does not belong to rank chunk {c}");
                                dense[(int)at] = t.Value;
                            }
                        }
                    }
                    if (c == 0) dense[0] = 0; // LCP[0] = 0
                    Deliver(output, start, dense.AsSpan(length));
                    rankSpool.Delete(c);
                }
                output.Complete(n);
            }
            baseFile?.Complete(fsync: false);
            overflowFile?.Complete(fsync: false);
        }
        finally
        {
            baseFile?.Dispose();
            overflowFile?.Dispose();
        }
    }

    private void Deliver(LcpOutputWriter output, long firstRank, ReadOnlySpan<uint> dense)
    {
        _consumer?.Consume(firstRank, dense);
        output.Write(firstRank, dense);
    }

    /// <summary>The width rule of ADR-002 section 4 on the counts of the whole text, before the first entry is written.</summary>
    private void ChooseKind(long n)
    {
        _stats.Kind = LcpWidth.Select(n, _escapes8, _escapes16);
        _stats.OverflowCount = LcpWidth.OverflowCount(_stats.Kind, _escapes8, _escapes16);
    }

    /// <summary>Runs <paramref name="work"/> on up to <c>min(Workers or Threads, chunks, what the budget allows)</c> threads (the first on the calling thread) and rethrows the first failure.</summary>
    private void RunWorkers(long perTaskBytes, long chunks, Action work)
    {
        int workers = _o.Workers > 0 ? _o.Workers : _o.Threads;
        if (_o.Workers == 0 && _budget.CapBytes != long.MaxValue)
        {
            long fit = (_budget.CapBytes - _budget.LiveBytes) / Math.Max(1, perTaskBytes);
            workers = (int)Math.Clamp(fit, 1, workers);
        }
        workers = (int)Math.Max(1, Math.Min(workers, chunks));
        _stats.Workers = Math.Max(_stats.Workers, workers);
        Interlocked.Exchange(ref _nextChunk, 0);
        void Guarded()
        {
            try
            {
                work();
            }
#pragma warning disable CA1031 // a failure of any kind in a worker thread is recorded and rethrown, with its own type and stack, on the calling thread (an unhandled exception would end the process)
            catch (Exception e)
            {
                lock (_gate) _error ??= e;
            }
#pragma warning restore CA1031
        }
        var extra = new List<Thread>();
        for (int t = 1; t < workers; t++)
        {
            var th = new Thread(Guarded) { IsBackground = true, Name = $"ph4-b2-{t}" };
            extra.Add(th);
            th.Start();
        }
        Guarded();
        foreach (Thread th in extra) th.Join();
        if (_error is { } error) ExceptionDispatchInfo.Capture(error).Throw();
    }

    private void Merge(PlcpTotals t)
    {
        lock (_gate)
        {
            _escapes8 += t.Escapes8;
            _escapes16 += t.Escapes16;
            _steps += t.Steps;
            _restart += t.Restart;
            _maxPlcp = Math.Max(_maxPlcp, t.Max);
        }
    }

    /// <summary>Sums of the PLCP loops of the chunks of one task.</summary>
    private sealed class PlcpTotals
    {
        internal long Escapes8;
        internal long Escapes16;
        internal long Steps;
        internal long Max;
        internal long Restart;

        internal void Add(PlcpChunkResult r)
        {
            Escapes8 += r.Escapes8;
            Escapes16 += r.Escapes16;
            Steps += r.CompareSteps;
            Restart += r.FirstValue;
            Max = Math.Max(Max, r.MaxPlcp);
        }
    }
}
