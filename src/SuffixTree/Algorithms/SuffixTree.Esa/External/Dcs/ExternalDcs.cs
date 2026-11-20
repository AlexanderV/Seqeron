using System.Diagnostics;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>
/// Parameters of the external-memory PH1 (ADR-004 section 4.2, prototype "E3b"): everything that sets a buffer size, a bucket capacity or a recursion threshold. The output never depends on any of them (ADR-004 A8, owner
/// answer Q5); the fuzz tests scale them down to K = 3..50, N_base = 4..200, X = 7..64 so that every multi-bucket and multi-level path runs on small texts.
/// </summary>
internal sealed record ExtDcsOptions
{
    /// <summary>K_s: records (sample positions) per sample bucket; a bucket of at most K_s records is sorted in RAM.</summary>
    internal long BucketCapacity { get; init; } = 1 << 20;

    /// <summary>N_base: a reduced string of at most this many symbols is sorted by SA-IS in RAM; longer ones go through the external levels (the owner's decision of 2026-11-18: 16 Mi symbols).</summary>
    internal long BaseCase { get; init; } = 16 << 20;

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

    /// <summary>Total of the read buffers of the merge (one per run file): each gets <c>min(IoBufferBytes, MergeBufferBudget / runs)</c>, so that this term does not grow with the number of buckets.</summary>
    internal long MergeBufferBudget { get; init; } = 16 << 20;

    /// <summary>Size of one sliding window of the mapped writers (the rank file, the suffix arrays of the levels below; at least 1 MiB in production, the E9 rule).</summary>
    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    /// <summary>K of the suffix sort of a reduced string at the external levels (records per leaf bucket of the position sort).</summary>
    internal long SortLeafCapacity { get; init; } = 1 << 20;

    /// <summary>Directory of the scratch files (created when missing; every file of a build is removed when it ends).</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Shared cap of the native buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Test hook: options of the arenas (guard bands, poison pattern); null = production arenas.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>Receives the explicit I/O, the mapped writes and the scratch bytes (null = not counted).</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>How the files are created (reservation of the disk space); null = <see cref="ReservationMode.Auto"/>.</summary>
    internal MappedFileOptions? Files { get; init; }

    /// <summary>Test hook, called with <c>"naming"</c> when a sample bucket task begins, <c>"merge"</c> before a level's merge, <c>"rank-written"</c> when a level has written its rank file (and, for the sorts of the levels below, as <see cref="SortOptions.Trace"/>); an exception thrown here fails the build like a failing step.</summary>
    internal Action<string>? Trace { get; init; }
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

    /// <summary>What the creation of the rank file did about the disk space; null when the file was created inside DIST (the names were not unique and the reduced string went through the external levels), which does not report it.</summary>
    internal ReservationOutcome? Reservation { get; set; }
}

/// <summary>
/// PH1 of ADR-004 built in external memory (prototype "E3b"): the rank file <c>rank.L0[slot]</c> of the difference-cover sample, byte for byte the table of the in-RAM builder of the prototype.
/// <para>
/// Level l (the text T_l, mapped read-only): the sample positions are named by <see cref="SampleNamer"/> (sample buckets, one scan of the sample per bucket) and the runs are merged into slot order. All names distinct: the names
/// are the ranks and go straight into the rank file. Otherwise the names form the reduced string (class-major, separators 0, see <see cref="SampleGeometry"/>) of <c>Samples + |D|</c> symbols over a 32-bit alphabet, and the ranks
/// of the sample suffixes are the ranks of the suffixes of that string: if it has at most N_base symbols it is sorted by SA-IS in RAM and the inverse is written to the rank file; otherwise it is written to a file (32 bits per
/// symbol, big endian: a reduced text for <see cref="MappedTextReader.OpenReduced"/>), its own rank file is built by this method one level down, its suffix array is built by the bucket sorter of PH2 / PH3
/// (<see cref="BucketSorter"/>) with that rank file, and <see cref="RankInverter"/> turns the suffix array into the rank file of level l. Every file is deleted as soon as its consumer is done; every big buffer comes from a
/// <see cref="BudgetedArena"/>; the text and the intermediate files are read through read-only mappings or bounded buffers. A failed build leaves no rank file, no scratch file (every file carries the tag of the build) and no
/// scratch byte in the counters, and the exception reaches the caller with its own type.
/// </para>
/// </summary>
internal static class ExternalDcs
{
    private const int MaxLevels = 64;

    /// <summary>Builds the rank file of <paramref name="text"/> at <paramref name="rankPath"/> (length <c>Slots * 4</c>, created through <see cref="MappedOutputFile.Create"/>), written through bounded mapped windows.</summary>
    /// <exception cref="IOException">A file cannot be created, written or spilled to (a full disk included); the rank file is removed.</exception>
    internal static ExtDcsStats BuildRankFile(MappedTextReader text, DifferenceCover cover, string rankPath, ExtDcsOptions options)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(cover);
        ArgumentNullException.ThrowIfNull(rankPath);
        ArgumentNullException.ThrowIfNull(options);
        ArgumentNullException.ThrowIfNull(options.ScratchDir);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.BucketCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(options.BucketCapacity, int.MaxValue / 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.BaseCase, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.Threads, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ChunkSlots, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.FanOut, 2);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ReservoirSize, 3);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.IoBufferBytes, 16);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.MergeBufferBudget, 8);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.SortLeafCapacity, 2);
        ArgumentOutOfRangeException.ThrowIfNegative(options.RunSpillBytes);
        if (text.Length > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{text.Length} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");
        string scratch = Path.GetFullPath(options.ScratchDir);
        bool createdRoot = !Directory.Exists(scratch);
        var context = new Context(options, scratch, "e3b-" + Guid.NewGuid().ToString("N")[..8], options.Budget ?? new ArenaBudget(long.MaxValue), new ScratchLedger(options.Counters));
        try
        {
            Directory.CreateDirectory(scratch);
            BuildLevel(text, cover, 0, rankPath, false, context);
        }
        catch
        {
            // a failed build leaves nothing behind: the files carry this build's tag, the rank file is ours from the moment it was asked for, the scratch bytes go back to the counters
            context.Ledger.ReleaseAll();
            if (context.RankTouched) RunStore.DeleteQuietly(rankPath);
            try
            {
                if (Directory.Exists(scratch))
                {
                    foreach (string file in Directory.GetFiles(scratch, context.Tag + "*")) RunStore.DeleteQuietly(file);
                }
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                // reported by the files that are left
            }
            throw;
        }
        finally
        {
            if (createdRoot) DeleteEmptyDirectoryQuietly(scratch);
        }
        return context.Stats;
    }

    private sealed class Context(ExtDcsOptions options, string scratch, string tag, ArenaBudget budget, ScratchLedger ledger)
    {
        internal ExtDcsOptions O { get; } = options;

        internal string Scratch { get; } = scratch;

        internal string Tag { get; } = tag;

        internal ArenaBudget Budget { get; } = budget;

        internal ScratchLedger Ledger { get; } = ledger;

        internal ExtDcsStats Stats { get; } = new();

        /// <summary>The rank file at the caller's path was asked for (so a failure removes it).</summary>
        internal bool RankTouched { get; set; }
    }

    private static void BuildLevel(MappedTextReader text, DifferenceCover cover, int level, string rankPath, bool scratchOutput, Context c)
    {
        if (level >= MaxLevels) throw new InvalidOperationException("the recursion of the sample ranks did not terminate");
        ExtDcsOptions o = c.O;
        var ls = new ExtDcsLevelStats { Level = level, TextLength = text.Length };
        c.Stats.Levels.Add(ls);
        var g = new SampleGeometry(text.Length, cover);
        ls.SampleSize = g.Samples;
        long t0 = Stopwatch.GetTimestamp();
        var runs = new RunStore(c.Scratch, $"{c.Tag}-L{level}", o.RunSpillBytes, c.Ledger, o.Counters);
        var namer = new SampleNamer(text, g, o, ls, c.Budget, runs);
        using NamingResult naming = namer.Run();
        ls.NamingSeconds = Seconds(t0);

        t0 = Stopwatch.GetTimestamp();
        if (ls.NamesUnique)
        {
            if (level == 0) c.RankTouched = true;
            using var sink = new RankFileSink(rankPath, g.Slots, o.WindowBytes, o.Files, o.Counters);
            if (level == 0) c.Stats.Reservation = sink.Outcome;
            if (scratchOutput) c.Ledger.Add(g.Slots * 4);
            namer.Merge(naming, sink);
            ls.MergeSeconds = Seconds(t0);
            o.Trace?.Invoke("rank-written");
            return;
        }

        long reduced = g.RedLength;
        // a reduced string that is not shorter than its text can only occur for tiny texts and tiny periods: sort it in RAM instead of recursing for ever
        if (reduced <= o.BaseCase || reduced >= text.Length)
        {
            ls.BaseCase = true;
            using var arena = new BudgetedArena(c.Budget, o.ArenaGuards);
            using var sink = new ReducedRamSink(g, arena);
            namer.Merge(naming, sink);
            ls.MergeSeconds = Seconds(t0);
            t0 = Stopwatch.GetTimestamp();
            naming.Dispose();
            BaseCase(sink.Red, g, ls.DistinctNames + 1, rankPath, scratchOutput, level, c, arena);
            arena.VerifyCanaries();
            ls.BelowSeconds = Seconds(t0);
            o.Trace?.Invoke("rank-written");
            return;
        }

        string redPath = Path.Combine(c.Scratch, $"{c.Tag}-L{level + 1}.red");
        string childRank = Path.Combine(c.Scratch, $"{c.Tag}-L{level + 1}.rank");
        string saPath = Path.Combine(c.Scratch, $"{c.Tag}-L{level + 1}.sa");
        using (var arena = new BudgetedArena(c.Budget, o.ArenaGuards))
        {
            var sink = new ReducedFileSink(redPath, g, o.IoBufferBytes, arena, c.Ledger, o.Counters);
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
        long redBytes = PackedTextFormat.TotalBytes(reduced, MappedTextReader.ReducedBits);
        using (MappedTextReader red = MappedTextReader.OpenReduced(redPath, reduced))
        {
            BuildLevel(red, cover, level + 1, childRank, true, c);
            long childSlots = new SampleGeometry(reduced, cover).Slots;
            using (RankTable table = RankTable.Open(childRank, cover, reduced))
            {
                using MappedOutputFile saFile = MappedOutputFile.Create(saPath, 4 * (reduced + 1), o.Files);
                c.Ledger.Add(4 * (reduced + 1));
                var sort = new SortOptions
                {
                    LeafCapacity = o.SortLeafCapacity,
                    FanOut = 64,
                    Threads = o.Threads,
                    Seed = o.Seed + (ulong)level + 1,
                    ScratchDir = c.Scratch,
                    SpillBytes = Math.Max(64, o.IoBufferBytes),
                    FilePrefix = $"{c.Tag}-L{level + 1}-pos",
                    Ledger = c.Ledger,
                    Budget = c.Budget,
                    ArenaGuards = o.ArenaGuards,
                    Counters = o.Counters,
                    Trace = o.Trace,
                };
                _ = BucketSorter.Run(red, sort, new MappedSaSliceSink(saFile, o.WindowBytes, o.Counters), table);
                saFile.Complete(fsync: false);
            }
            File.Delete(childRank);
            c.Ledger.Remove(childSlots * 4);
        }
        File.Delete(redPath);
        c.Ledger.Remove(redBytes);
        if (level == 0) c.RankTouched = true;
        if (scratchOutput) c.Ledger.Add(g.Slots * 4);
        var dist = new DistOptions
        {
            ScratchDir = c.Scratch,
            ChunkEntries = o.ChunkSlots,
            IoBufferBytes = o.IoBufferBytes,
            ReadBytes = Math.Max(4, o.IoBufferBytes),
            WindowBytes = o.WindowBytes,
            Budget = c.Budget,
            Counters = o.Counters,
            Arena = o.ArenaGuards,
            Files = o.Files ?? new MappedFileOptions(),
        };
        RankInverter.Invert(saPath, g, rankPath, dist);
        File.Delete(saPath);
        c.Ledger.Remove(4 * (reduced + 1));
        ls.BelowSeconds = Seconds(t0);
        o.Trace?.Invoke("rank-written");
    }

    /// <summary>The base case: SA-IS on the reduced string in RAM, the inverse permutation, and the rank file written window by window.</summary>
    private static void BaseCase(Buf<int> red, SampleGeometry g, long sigma, string rankPath, bool scratchOutput, int level, Context c, BudgetedArena arena)
    {
        ExtDcsOptions o = c.O;
        int n = checked((int)g.RedLength);
        Buf<int> sa = arena.Allocate<int>(n + 1, zero: false, "base.sa");
        SaisReference.Build(new IntSymbols(red.AsSpan(n)), checked((int)sigma), sa.AsSpan(n + 1));
        Span<int> rank = red.AsSpan(n); // the reduced string is dead now: it receives the inverse permutation
        ReadOnlySpan<int> saSpan = sa.AsSpan(n + 1);
        for (int k = 1; k <= n; k++) rank[saSpan[k]] = k;

        long slots = g.Slots;
        long chunkSlots = Math.Max(1, o.ChunkSlots);
        Buf<uint> chunk = arena.Allocate<uint>((int)Math.Min(chunkSlots, slots), zero: false, "base.chunk");
        if (level == 0) c.RankTouched = true;
        using MappedOutputFile file = MappedOutputFile.Create(rankPath, slots * 4, o.Files);
        if (level == 0) c.Stats.Reservation = file.Outcome;
        if (scratchOutput) c.Ledger.Add(slots * 4);
        using MappedWindowWriter writer = file.CreateWriter(o.WindowBytes, o.Counters);
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

    private static void DeleteEmptyDirectoryQuietly(string path)
    {
        try
        {
            if (Directory.Exists(path) && !Directory.EnumerateFileSystemEntries(path).Any()) Directory.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // someone else put a file there between the check and the delete: it is theirs
        }
    }
}
