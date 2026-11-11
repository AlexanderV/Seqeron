using System.Diagnostics;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Dcs;

/// <summary>
/// Parameters of the external-memory PH1 (ADR-004 section 4.2, prototype "E3b"): everything that sets a buffer size, a bucket capacity or a recursion threshold. The output never depends on
/// any of them (ADR-004 A8, owner answer Q5); the fuzz tests scale them down to K = 3..50, N_base = 4..200, X = 7..64 so that every multi-bucket and multi-level path runs on small texts.
/// </summary>
internal sealed record ExtDcsOptions
{
    /// <summary>K_s: records (sample positions) per sample bucket; a bucket of at most K_s records is sorted in RAM.</summary>
    internal long BucketCapacity { get; init; } = 1 << 20;

    /// <summary>N_base: a reduced string of at most this many symbols is sorted by SA-IS in RAM; longer ones go through the external levels.</summary>
    internal long BaseCase { get; init; } = 4 << 20;

    /// <summary>F_max: largest number of sample buckets one split creates.</summary>
    internal int FanOut { get; init; } = 1024;

    /// <summary>Records kept by the reservoir that chooses the splitters of one split.</summary>
    internal int ReservoirSize { get; init; } = 16384;

    /// <summary>Worker threads of the bucket tasks and of the suffix sort of the external levels.</summary>
    internal int Threads { get; init; } = 1;

    /// <summary>Seed of the reservoir sampling; the output does not depend on it, the shape of the bucket tree does.</summary>
    internal ulong Seed { get; init; } = 1;

    /// <summary>R: slots per merge window and per chunk of the suffix-array inversion (the staging array holds R four-byte entries).</summary>
    internal long ChunkSlots { get; init; } = 4 << 20;

    /// <summary>A bucket run of at most this many bytes stays in memory, larger runs go to a scratch file.</summary>
    internal int RunSpillBytes { get; init; } = 64 * 1024;

    /// <summary>Size of one read or write buffer of the merge, of the tuple files and of the class cursors.</summary>
    internal int IoBufferBytes { get; init; } = 64 * 1024;

    /// <summary>Size of one sliding window of the mapped rank file writer.</summary>
    internal int WindowBytes { get; init; } = 1 << 20;

    /// <summary>K of the suffix sort of a reduced string at the external levels (records per leaf bucket of the position sort).</summary>
    internal long SortLeafCapacity { get; init; } = 1 << 20;

    /// <summary>Directory of the scratch files (created when missing; every file of a build is removed when it ends).</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Shared cap of the native buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }
}

/// <summary>What one level of the recursion did.</summary>
internal sealed class ExtDcsLevelStats
{
    /// <summary>0 = the input text, 1 = the reduced string of level 0, ...</summary>
    internal int Level { get; set; }

    /// <summary>Symbols of the text of this level.</summary>
    internal long TextLength { get; set; }

    /// <summary>Sample positions of this level.</summary>
    internal long SampleSize { get; set; }

    /// <summary>Distinct names (distinct X-prefixes).</summary>
    internal long DistinctNames { get; set; }

    /// <summary>Sample positions whose name is shared with at least one other sample position.</summary>
    internal long NonUniqueRecords { get; set; }

    /// <summary>All names distinct: the names are the ranks and there is no level below.</summary>
    internal bool NamesUnique => DistinctNames == SampleSize;

    /// <summary>The reduced string of this level was sorted by SA-IS in RAM (the recursion ended here).</summary>
    internal bool BaseCase { get; set; }

    /// <summary>Leaf buckets that were sorted.</summary>
    internal long Leaves { get; set; }

    /// <summary>Scans of the sample positions by tasks that split a bucket instead of sorting it (the root and the buckets that overflowed).</summary>
    internal long SplitScans { get; set; }

    /// <summary>Scans of the sample positions by all tasks (one per leaf, one per split, one per overflow).</summary>
    internal long Scans { get; set; }

    /// <summary>Leaves whose run went to a scratch file.</summary>
    internal long FileRuns { get; set; }

    /// <summary>Seconds of the bucket tasks, of the merge into the sink, and of everything below this level.</summary>
    internal double NamingSeconds { get; set; }

    internal double MergeSeconds { get; set; }

    internal double BelowSeconds { get; set; }
}

/// <summary>What <see cref="ExternalDcs.BuildRankFile"/> did, level by level.</summary>
internal sealed class ExtDcsStats
{
    internal List<ExtDcsLevelStats> Levels { get; } = [];

    /// <summary>Index of the deepest level reached.</summary>
    internal int MaxLevel => Levels.Count - 1;
}

/// <summary>
/// PH1 of ADR-004 built in external memory (prototype "E3b"): the rank file <c>rank.L0[slot]</c> of the difference-cover sample, byte for byte the table of <see cref="DcsBuilder"/>.
/// <para>
/// Level l (the text T_l, mapped read-only): the sample positions are named by <see cref="SampleNamer"/> (sample buckets, one scan of the sample per bucket) and the runs are merged into slot order.
/// All names distinct: the names are the ranks and go straight into the rank file. Otherwise the names form the reduced string (class-major, separators 0, see <see cref="SampleGeometry"/>) of
/// <c>Samples + |D|</c> symbols over a 32-bit alphabet, and the ranks of the sample suffixes are the ranks of the suffixes of that string: if it has at most N_base symbols it is sorted by SA-IS in RAM
/// and the inverse is written to the rank file; otherwise it is written to a file (32 bits per symbol, big endian: a <see cref="PackedText"/>), its own rank file is built by this method one level
/// down, its suffix array is built by the bucket sorter of E3 (<see cref="EmBuilder"/>) with that rank file, and <see cref="SaInverter"/> turns the suffix array into the rank file of level l.
/// Every file is deleted as soon as its consumer is done; every big buffer comes from a <see cref="BudgetedArena"/>; the text and the intermediate files are read through read-only mappings or bounded buffers.
/// </para>
/// </summary>
internal static unsafe class ExternalDcs
{
    private const int MaxLevels = 64;

    /// <summary>Builds the rank file of <paramref name="text"/> at <paramref name="rankPath"/> (length <c>Slots * 4</c>), written through bounded mapped windows.</summary>
    internal static ExtDcsStats BuildRankFile(PackedText text, DifferenceCover cover, string rankPath, ExtDcsOptions options)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(cover);
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.BucketCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(options.BucketCapacity, int.MaxValue / 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.BaseCase, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Threads, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ChunkSlots, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.FanOut, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ReservoirSize, 3);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.IoBufferBytes, 16);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.SortLeafCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfNegative(options.RunSpillBytes);
        if (text.Length >= uint.MaxValue - 1) throw new NotSupportedException("positions are 32-bit in the prototype");
        Directory.CreateDirectory(options.ScratchDir);
        string tag = "e3b-" + Guid.NewGuid().ToString("N")[..8];
        var stats = new ExtDcsStats();
        try
        {
            BuildLevel(text, cover, 0, rankPath, false, options, tag, options.Budget ?? new ArenaBudget(long.MaxValue), stats);
        }
        catch
        {
            // a failed build leaves no scratch behind (its files carry this build's tag)
            foreach (string file in Directory.GetFiles(options.ScratchDir, tag + "*"))
            {
                try
                {
                    File.Delete(file);
                }
                catch (IOException)
                {
                    // best effort
                }
            }
            throw;
        }
        return stats;
    }

    private static void BuildLevel(PackedText text, DifferenceCover cover, int level, string rankPath, bool scratchOutput, ExtDcsOptions o, string tag, ArenaBudget budget, ExtDcsStats stats)
    {
        if (level >= MaxLevels) throw new InvalidOperationException("the recursion of the sample ranks did not terminate");
        var ls = new ExtDcsLevelStats { Level = level, TextLength = text.Length };
        stats.Levels.Add(ls);
        var g = new SampleGeometry(text.Length, cover);
        ls.SampleSize = g.Samples;
        long t0 = Stopwatch.GetTimestamp();
        var namer = new SampleNamer(text, g, o, ls, budget, $"{tag}-L{level}");
        using NamingResult naming = namer.Run();
        ls.NamingSeconds = Seconds(t0);

        t0 = Stopwatch.GetTimestamp();
        if (ls.NamesUnique)
        {
            using var sink = new RankFileSink(rankPath, g.Slots, o.WindowBytes, scratchOutput);
            namer.Merge(naming, sink);
            ls.MergeSeconds = Seconds(t0);
            return;
        }

        long reduced = g.RedLength;
        // a reduced string that is not shorter than its text can only occur for tiny texts and tiny periods: sort it in RAM instead of recursing for ever
        if (reduced <= o.BaseCase || reduced >= text.Length)
        {
            ls.BaseCase = true;
            using var arena = new BudgetedArena(budget);
            using var sink = new ReducedRamSink(g, arena);
            namer.Merge(naming, sink);
            ls.MergeSeconds = Seconds(t0);
            t0 = Stopwatch.GetTimestamp();
            naming.Dispose();
            BaseCase(sink.Red, g, ls.DistinctNames + 1, rankPath, scratchOutput, o, arena);
            ls.BelowSeconds = Seconds(t0);
            return;
        }

        string redPath = Path.Combine(o.ScratchDir, $"{tag}-L{level + 1}.red");
        string childRank = Path.Combine(o.ScratchDir, $"{tag}-L{level + 1}.rank");
        string saPath = Path.Combine(o.ScratchDir, $"{tag}-L{level + 1}.sa");
        using (var arena = new BudgetedArena(budget))
        {
            var sink = new ReducedFileSink(redPath, g, o.IoBufferBytes, arena);
            try
            {
                namer.Merge(naming, sink);
                sink.Complete();
            }
            finally
            {
                sink.Dispose();
            }
        }
        ls.MergeSeconds = Seconds(t0);
        naming.Dispose();

        t0 = Stopwatch.GetTimestamp();
        long redBytes = (reduced * 4) + PackedText.PadBytes;
        using (var mappedRed = new MappedFile(redPath, writable: false))
        using (var red = new PackedText(mappedRed.Pointer, reduced, 32, mappedRed, ownsOwner: false))
        {
            BuildLevel(red, cover, level + 1, childRank, true, o, tag, budget, stats);
            long childSlots = new SampleGeometry(reduced, cover).Slots;
            using (var table = new DcsRankTable(cover, new MappedFile(childRank, writable: false), reduced))
            {
                ScratchTracker.Add((reduced + 1) * 4);
                using (var saSink = new MappedSaSink(saPath, reduced + 1))
                {
                    var build = new BuildOptions
                    {
                        LeafCapacity = o.SortLeafCapacity,
                        FanOut = 64,
                        Threads = o.Threads,
                        Seed = o.Seed + (ulong)level + 1,
                        ScratchDir = o.ScratchDir,
                        SpillBytes = Math.Max(64, o.IoBufferBytes),
                        Budget = budget,
                    };
                    _ = EmBuilder.Run(red, build, saSink, table);
                }
            }
            File.Delete(childRank);
            ScratchTracker.Remove(childSlots * 4);
        }
        File.Delete(redPath);
        ScratchTracker.Remove(redBytes);
        SaInverter.Invert(saPath, reduced, g, rankPath, scratchOutput, o, $"{tag}-L{level}", budget);
        File.Delete(saPath);
        ScratchTracker.Remove((reduced + 1) * 4);
        ls.BelowSeconds = Seconds(t0);
    }

    /// <summary>The base case: SA-IS on the reduced string in RAM, the inverse permutation, and the rank file written window by window.</summary>
    private static void BaseCase(Buf<int> red, SampleGeometry g, long sigma, string rankPath, bool scratchOutput, ExtDcsOptions o, BudgetedArena arena)
    {
        int n = checked((int)g.RedLength);
        Buf<int> sa = arena.Allocate<int>(n + 1, zero: false, "base.sa");
        SaisReference.Build(new IntSymbols(red.AsSpan(n)), checked((int)sigma), sa.AsSpan(n + 1));
        Span<int> rank = red.AsSpan(n); // the reduced string is dead now: it receives the inverse permutation
        ReadOnlySpan<int> saSpan = sa.AsSpan(n + 1);
        for (int k = 1; k <= n; k++) rank[saSpan[k]] = k;

        long slots = g.Slots;
        long chunkSlots = Math.Max(1, o.ChunkSlots);
        Buf<uint> chunk = arena.Allocate<uint>((int)Math.Min(chunkSlots, slots), zero: false, "base.chunk");
        using var writer = new MappedWindowWriter(rankPath, slots * 4, o.WindowBytes);
        if (scratchOutput) ScratchTracker.Add(slots * 4);
        for (long start = 0; start < slots; start += chunkSlots)
        {
            int len = (int)Math.Min(chunkSlots, slots - start);
            Span<uint> window = chunk.AsSpan(len);
            for (int i = 0; i < len; i++)
            {
                long slot = start + i;
                window[i] = slot < g.Samples ? (uint)rank[(int)g.RedIndexOfSlot(slot)] : 0;
            }
            writer.Append(MemoryMarshal.AsBytes(window));
        }
    }

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
}

/// <summary>Counters of the writes through mapped windows (the explicit pread / pwrite calls are counted by <see cref="IoCounters"/>).</summary>
internal static class MappedIo
{
    private static long s_writeBytes;
    private static long s_windows;
    private static long s_maxWindow;

    internal static long WriteBytes => Interlocked.Read(ref s_writeBytes);

    /// <summary>Windows (views) created.</summary>
    internal static long Windows => Interlocked.Read(ref s_windows);

    /// <summary>Largest window created, in bytes.</summary>
    internal static long MaxWindowBytes => Interlocked.Read(ref s_maxWindow);

    internal static void Reset()
    {
        Interlocked.Exchange(ref s_writeBytes, 0);
        Interlocked.Exchange(ref s_windows, 0);
        Interlocked.Exchange(ref s_maxWindow, 0);
    }

    internal static void AddWrite(long bytes) => Interlocked.Add(ref s_writeBytes, bytes);

    internal static void AddWindow(long bytes)
    {
        Interlocked.Increment(ref s_windows);
        Atomic.RaiseTo(ref s_maxWindow, bytes);
    }
}
