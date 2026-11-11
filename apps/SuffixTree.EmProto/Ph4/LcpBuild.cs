using System.Diagnostics;
using System.Runtime.ExceptionServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.X86;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Ph4;

/// <summary>
/// One run of PH4 (ADR-004 section 4.2, row PH4). Base mode: B1 scans the suffix array file in rank order and distributes the tuples <c>(SA[i], SA[i-1], i)</c> by <c>SA[i] / R</c> into P chunk files; B2 takes the
/// chunks (several in parallel), scatters Phi and the rank into two arrays of R entries, runs the PLCP loop (<see cref="PlcpScan"/>), counts the escapes (so that the width is known before any LCP entry is
/// written) and distributes <c>(i, PLCP[j])</c> by <c>i / R</c>; B3 takes the rank chunks in order, scatters them into a dense array, hands it to the chunk consumer (the stream point of PH5) and writes the compact
/// form through mapped windows. Bounded mode: the same B1 and B2 over text-range rounds, the PLCP of every round goes to a text-order file, and B3 gathers <c>PLCP[SA[i]]</c> in rank order.
/// The output does not depend on the thread count, R, the buffer sizes, the mode or the rounds.
/// </summary>
internal sealed unsafe class LcpBuild
{
    private const int GatherPrefetch = 16;

    private readonly string _saPath;
    private readonly PackedText _text;
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
    private long _plcpFileBytes;

    internal LcpBuild(string saPath, PackedText text, string basePath, string overflowPath, ExtLcpOptions options, ILcpChunkConsumer? consumer)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ChunkPositions, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Threads, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(options.Workers);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ReadBufferBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.TupleBufferBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 1);
        ArgumentOutOfRangeException.ThrowIfNegative(options.RoundChunks);
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
        var g = new LcpGeometry(n, o.ChunkPositions);
        long chunk = Math.Min(o.ChunkPositions, Math.Max(n, 1));
        if (o.Mode == LcpMode.Bounded) return (4 * chunk) + (Math.Max(1, o.ReadBufferBytes / sizeof(PhiPair)) * (long)sizeof(PhiPair));
        long tuples = Math.Max(1, o.TupleBufferBytes / sizeof(RankTuple)) * (long)sizeof(RankTuple) * g.RankChunks;
        return (8 * chunk) + (Math.Max(1, o.ReadBufferBytes / sizeof(PhiTuple)) * (long)sizeof(PhiTuple)) + tuples;
    }

    internal ExtLcpStats Run()
    {
        long n = _text.Length;
        CheckSuffixArrayFile(n);
        Directory.CreateDirectory(_o.ScratchDir);
        _stats.N = n;
        _stats.TextChunks = _g.TextChunks;
        _stats.RankChunks = _g.RankChunks;
        if (_o.Mode == LcpMode.Bounded) RunBounded();
        else RunDistributed();
        _stats.Escapes8 = _escapes8;
        _stats.Escapes16 = _escapes16;
        _stats.CompareSteps = _steps;
        _stats.MaxPlcp = _maxPlcp;
        _stats.RestartSymbols = _restart;
        return _stats;
    }

    private int ReadEntries(int tupleBytes) => Math.Max(1, _o.ReadBufferBytes / tupleBytes);

    private int ChunkEntries(long chunkLength) => checked((int)Math.Min(chunkLength, _g.R));

    /// <summary>The file must hold <c>n + 1</c> entries and start with <c>SA[0] = n</c> (the empty suffix).</summary>
    private void CheckSuffixArrayFile(long n)
    {
        long length = new FileInfo(_saPath).Length;
        if (length != (n + 1) * 4) throw new ArgumentException($"the suffix array file needs n + 1 = {n + 1} entries of 4 bytes, it has {length} bytes");
        Span<uint> first = stackalloc uint[1];
        using var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read);
        sa.ReadExactAt(MemoryMarshal.AsBytes(first), 0);
        if (first[0] != n) throw new ArgumentException($"sa[0] must be n = {n}, the empty suffix, it is {first[0]}");
    }

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // base mode

    private void RunDistributed()
    {
        long n = _text.Length;
        _stats.Rounds = 1;
        string dir = _o.ScratchDir;
        using var phiSpool = new Spool<PhiTuple>(dir, "ph4phi", _g.TextChunks);
        using var rankSpool = new Spool<RankTuple>(dir, "ph4rank", _g.RankChunks);

        long t0 = Stopwatch.GetTimestamp();
        if (n > 0) DistributePhi(phiSpool);
        long t1 = Stopwatch.GetTimestamp();
        if (n > 0) RunWorkers(WorkerBytes(n, _o), _g.TextChunks, 0, () => PlcpWorker(phiSpool, rankSpool));
        long t2 = Stopwatch.GetTimestamp();
        ChooseKind(n);
        EmitRankChunks(rankSpool);
        long t3 = Stopwatch.GetTimestamp();
        _stats.B1Seconds = (t1 - t0) / (double)Stopwatch.Frequency;
        _stats.B2Seconds = (t2 - t1) / (double)Stopwatch.Frequency;
        _stats.B3Seconds = (t3 - t2) / (double)Stopwatch.Frequency;
    }

    /// <summary>B1: the suffix array in rank order, one tuple per entry into the chunk file of its position.</summary>
    private void DistributePhi(Spool<PhiTuple> spool)
    {
        long n = _text.Length;
        using var arena = new BudgetedArena(_budget, null, _o.TouchBuffers);
        int readEntries = ReadEntries(sizeof(uint));
        Buf<uint> read = arena.Allocate<uint>(readEntries, zero: false, "b1.sa");
        var writer = new SpoolWriter<PhiTuple>(spool, arena, 0, _g.TextChunks, _o.TupleBufferBytes);
        using var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read);
        uint previous = (uint)n; // SA[0], the empty suffix
        for (long first = 1; first <= n; first += readEntries)
        {
            int take = (int)Math.Min(readEntries, n + 1 - first);
            sa.ReadExactAt(MemoryMarshal.AsBytes(read.AsSpan(take)), first * 4);
            for (int i = 0; i < take; i++)
            {
                uint s = read[i];
                if (s >= n) throw new ArgumentException($"sa[{first + i}] = {s} is outside [0, {n})");
                writer.Add(_g.ChunkOf(s), new PhiTuple(s, previous, (uint)(first + i)));
                previous = s;
            }
        }
        writer.FlushAll();
    }

    /// <summary>B2 of one task: chunks in any order; Phi and rank scattered into arrays of R entries, PLCP in place over Phi, the second distribution by rank.</summary>
    private void PlcpWorker(Spool<PhiTuple> phiSpool, Spool<RankTuple> rankSpool)
    {
        using var arena = new BudgetedArena(_budget, null, _o.TouchBuffers);
        int maxChunk = ChunkEntries(_g.TextChunkLength(0));
        Buf<uint> phi = arena.Allocate<uint>(maxChunk, zero: false, "b2.phi");
        Buf<uint> rank = arena.Allocate<uint>(maxChunk, zero: false, "b2.rank");
        int readEntries = ReadEntries(sizeof(PhiTuple));
        Buf<PhiTuple> read = arena.Allocate<PhiTuple>(readEntries, zero: false, "b2.read");
        var writer = new SpoolWriter<RankTuple>(rankSpool, arena, 0, _g.RankChunks, _o.TupleBufferBytes);
        var total = new PlcpTotals();
        while (Volatile.Read(ref _error) is null)
        {
            long c = Interlocked.Increment(ref _nextChunk) - 1;
            if (c >= _g.TextChunks) break;
            int length = (int)_g.TextChunkLength(c);
            long start = c * _g.R;
            using (SpoolReader<PhiTuple> reader = phiSpool.OpenReader(c))
            {
                if (reader.Count != length) throw new InvalidOperationException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");
                int got;
                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)
                {
                    for (int i = 0; i < got; i++)
                    {
                        PhiTuple t = read[i];
                        long at = t.Position - start;
                        if ((ulong)at >= (ulong)length) throw new InvalidOperationException($"position {t.Position} does not belong to text chunk {c}");
                        phi[(int)at] = t.Phi;
                        rank[(int)at] = t.Rank;
                    }
                }
            }
            total.Add(PlcpScan.Run(_text, start, phi, length));
            for (int k = 0; k < length; k++)
            {
                uint r = rank[k];
                writer.Add(_g.ChunkOf(r), new RankTuple(r, phi[k]));
            }
            phiSpool.Delete(c);
        }
        if (Volatile.Read(ref _error) is null) writer.FlushAll();
        Merge(total);
    }

    /// <summary>B3: rank chunks in ascending order, scattered into a dense array, handed to the consumer and written in the chosen width.</summary>
    private void EmitRankChunks(Spool<RankTuple> rankSpool)
    {
        long n = _text.Length;
        using var arena = new BudgetedArena(_budget, null, _o.TouchBuffers);
        int maxChunk = ChunkEntries(_g.RankChunkLength(0));
        Buf<uint> dense = arena.Allocate<uint>(maxChunk, zero: false, "b3.dense");
        int readEntries = ReadEntries(sizeof(RankTuple));
        Buf<RankTuple> read = arena.Allocate<RankTuple>(readEntries, zero: false, "b3.read");
        using LcpOutputWriter output = OpenOutput(arena, n);
        for (long c = 0; c < _g.RankChunks; c++)
        {
            int length = (int)_g.RankChunkLength(c);
            long start = c * _g.R;
            long expected = c == 0 ? length - 1 : length; // rank 0 is the empty suffix: no tuple
            using (SpoolReader<RankTuple> reader = rankSpool.OpenReader(c))
            {
                if (reader.Count != expected) throw new InvalidOperationException($"rank chunk {c} holds {reader.Count} tuples, expected {expected}: the suffix array is not a permutation");
                int got;
                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)
                {
                    for (int i = 0; i < got; i++)
                    {
                        RankTuple t = read[i];
                        long at = t.Rank - start;
                        if ((ulong)at >= (ulong)length) throw new InvalidOperationException($"rank {t.Rank} does not belong to rank chunk {c}");
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

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // bounded mode

    private void RunBounded()
    {
        long n = _text.Length;
        string dir = _o.ScratchDir;
        string plcpPath = Path.Combine(dir, "ph4plcp.bin");
        long perRound = _o.RoundChunks <= 0 || _o.RoundChunks >= _g.TextChunks ? Math.Max(1, _g.TextChunks) : _o.RoundChunks;
        long rounds = _g.TextChunks == 0 ? 0 : (_g.TextChunks + perRound - 1) / perRound;
        _stats.Rounds = (int)rounds;
        double b1 = 0, b2 = 0;
        try
        {
            if (n > 0)
            {
                using var fs = new FileStream(plcpPath, FileMode.Create, FileAccess.Write, FileShare.ReadWrite);
                fs.SetLength(n * 4);
            }
            for (long round = 0; round < rounds; round++)
            {
                long first = round * perRound;
                long last = Math.Min(_g.TextChunks, first + perRound);
                using var spool = new Spool<PhiPair>(dir, $"ph4phi{round}", _g.TextChunks);
                long t0 = Stopwatch.GetTimestamp();
                DistributePairs(spool, first, last);
                long t1 = Stopwatch.GetTimestamp();
                RunWorkers(WorkerBytes(n, _o), last - first, first, () => PairWorker(spool, plcpPath, last));
                long t2 = Stopwatch.GetTimestamp();
                b1 += (t1 - t0) / (double)Stopwatch.Frequency;
                b2 += (t2 - t1) / (double)Stopwatch.Frequency;
            }
            ChooseKind(n);
            long t3 = Stopwatch.GetTimestamp();
            GatherRankChunks(plcpPath);
            _stats.B3Seconds = (Stopwatch.GetTimestamp() - t3) / (double)Stopwatch.Frequency;
            _stats.B1Seconds = b1;
            _stats.B2Seconds = b2;
        }
        finally
        {
            File.Delete(plcpPath);
            ScratchTracker.Remove(_plcpFileBytes);
        }
    }

    /// <summary>B1 of one round: the whole suffix array is read again, only the tuples of the positions in chunks <c>[first, last)</c> are kept.</summary>
    private void DistributePairs(Spool<PhiPair> spool, long first, long last)
    {
        long n = _text.Length;
        using var arena = new BudgetedArena(_budget, null, _o.TouchBuffers);
        int readEntries = ReadEntries(sizeof(uint));
        Buf<uint> read = arena.Allocate<uint>(readEntries, zero: false, "b1.sa");
        var writer = new SpoolWriter<PhiPair>(spool, arena, first, last - first, _o.TupleBufferBytes);
        using var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read);
        uint previous = (uint)n;
        for (long at = 1; at <= n; at += readEntries)
        {
            int take = (int)Math.Min(readEntries, n + 1 - at);
            sa.ReadExactAt(MemoryMarshal.AsBytes(read.AsSpan(take)), at * 4);
            for (int i = 0; i < take; i++)
            {
                uint s = read[i];
                if (s >= n) throw new ArgumentException($"sa[{at + i}] = {s} is outside [0, {n})");
                uint chunk = _g.ChunkOf(s);
                if (chunk >= first && chunk < last) writer.Add(chunk, new PhiPair(s, previous));
                previous = s;
            }
        }
        writer.FlushAll();
    }

    /// <summary>B2 of one task in the bounded mode: Phi chunks of the round, PLCP written to its place in the text-order file (four bytes per position).</summary>
    private void PairWorker(Spool<PhiPair> spool, string plcpPath, long lastChunk)
    {
        using var arena = new BudgetedArena(_budget, null, _o.TouchBuffers);
        int maxChunk = ChunkEntries(_g.TextChunkLength(0));
        Buf<uint> phi = arena.Allocate<uint>(maxChunk, zero: false, "b2.phi");
        int readEntries = ReadEntries(sizeof(PhiPair));
        Buf<PhiPair> read = arena.Allocate<PhiPair>(readEntries, zero: false, "b2.read");
        using var plcp = new CountedFile(plcpPath, FileMode.Open, FileAccess.Write);
        var total = new PlcpTotals();
        while (Volatile.Read(ref _error) is null)
        {
            long c = Volatile.Read(ref _roundFirst) + Interlocked.Increment(ref _nextChunk) - 1;
            if (c >= lastChunk) break;
            int length = (int)_g.TextChunkLength(c);
            long start = c * _g.R;
            using (SpoolReader<PhiPair> reader = spool.OpenReader(c))
            {
                if (reader.Count != length) throw new InvalidOperationException($"text chunk {c} holds {reader.Count} Phi tuples, expected {length}: the suffix array is not a permutation");
                int got;
                while ((got = reader.Read(read.AsSpan(readEntries))) > 0)
                {
                    for (int i = 0; i < got; i++)
                    {
                        PhiPair t = read[i];
                        long at = t.Position - start;
                        if ((ulong)at >= (ulong)length) throw new InvalidOperationException($"position {t.Position} does not belong to text chunk {c}");
                        phi[(int)at] = t.Phi;
                    }
                }
            }
            total.Add(PlcpScan.Run(_text, start, phi, length));
            plcp.WriteAt(MemoryMarshal.AsBytes(phi.AsSpan(length)), start * 4);
            Interlocked.Add(ref _plcpFileBytes, length * 4L);
            ScratchTracker.Add(length * 4L);
            spool.Delete(c);
        }
        Merge(total);
    }

    private long _roundFirst;

    /// <summary>B3 of the bounded mode: the suffix array in rank order again, <c>LCP[i] = PLCP[SA[i]]</c> through a read-only mapping of the PLCP file.</summary>
    private void GatherRankChunks(string plcpPath)
    {
        long n = _text.Length;
        using var arena = new BudgetedArena(_budget, null, _o.TouchBuffers);
        int maxChunk = ChunkEntries(_g.RankChunkLength(0));
        Buf<uint> dense = arena.Allocate<uint>(maxChunk, zero: false, "b3.dense");
        int readEntries = ReadEntries(sizeof(uint));
        Buf<uint> read = arena.Allocate<uint>(readEntries, zero: false, "b3.sa");
        using LcpOutputWriter output = OpenOutput(arena, n);
        using var sa = new CountedFile(_saPath, FileMode.Open, FileAccess.Read);
        using MappedFile? map = n > 0 ? new MappedFile(plcpPath, writable: false) : null;
        uint* plcp = map is null ? null : (uint*)map.Pointer;
        for (long c = 0; c < _g.RankChunks; c++)
        {
            int length = (int)_g.RankChunkLength(c);
            long start = c * _g.R;
            for (int done = 0; done < length; done += readEntries)
            {
                int take = Math.Min(readEntries, length - done);
                sa.ReadExactAt(MemoryMarshal.AsBytes(read.AsSpan(take)), (start + done) * 4);
                for (int i = 0; i < take; i++)
                {
                    long rank = start + done + i;
                    if (rank == 0)
                    {
                        dense[0] = 0; // LCP[0] = 0
                        continue;
                    }
                    uint s = read[i];
                    if (s >= n) throw new ArgumentException($"sa[{rank}] = {s} is outside [0, {n})");
                    if (Sse.IsSupported && i + GatherPrefetch < take) Sse.Prefetch0(plcp + read[i + GatherPrefetch]);
                    dense[done + i] = plcp[s];
                }
            }
            Deliver(output, start, dense.AsSpan(length));
        }
        output.Complete(n);
    }

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // shared

    private LcpOutputWriter OpenOutput(BudgetedArena arena, long n)
    {
        long overflow = _stats.OverflowCount;
        return new LcpOutputWriter(arena, _stats.Kind, n, overflow, _basePath, _overflowPath, _o.WindowBytes, ReadEntries(sizeof(uint)));
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
    private void RunWorkers(long perTaskBytes, long chunks, long firstChunk, Action work)
    {
        int workers = _o.Workers > 0 ? _o.Workers : _o.Threads;
        if (_o.Workers == 0 && _budget.CapBytes != long.MaxValue)
        {
            long fit = (_budget.CapBytes - _budget.LiveBytes) / Math.Max(1, perTaskBytes);
            workers = (int)Math.Clamp(fit, 1, workers);
        }
        workers = (int)Math.Max(1, Math.Min(workers, chunks));
        _stats.Workers = Math.Max(_stats.Workers, workers);
        Volatile.Write(ref _roundFirst, firstChunk);
        Interlocked.Exchange(ref _nextChunk, 0);
        void Guarded()
        {
            try
            {
                work();
            }
            catch (Exception e)
            {
                lock (_gate) _error ??= e;
            }
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
