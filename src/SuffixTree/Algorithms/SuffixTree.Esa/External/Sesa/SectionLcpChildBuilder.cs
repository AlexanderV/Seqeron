using System.Diagnostics;

namespace SuffixTree.Esa.External;

/// <summary>
/// PH4 (base mode) and PH5 (layout A) as <see cref="LcpChildBuilder"/> runs them, with the outputs as sections of one container file instead of three files of their own (work package 1b.7e): the suffix array is read from a
/// byte range of a file (<see cref="ExtLcpOptions.SaOffset"/>), the CHILD section is a range of the container (known up front: <c>4 (n + 1)</c> bytes), and the LCP and its overflow list are placed by a callback once
/// the width is known (after B2, before the first LCP entry is written), so that the container can grow behind the sections that were known after PH0. Nothing is copied. The stage creates and removes its own scratch
/// subdirectories <c>ph4</c> and <c>ph5</c>; the container file is the owner's (it is neither completed nor deleted here, also not on a failure).
/// </summary>
internal static class SectionLcpChildBuilder
{
    /// <summary>
    /// Builds the LCP, its overflow list and the child table of <paramref name="text"/> from the suffix array at byte <paramref name="saOffset"/> of <paramref name="saPath"/>. The child table is written to
    /// <c>[childOffset, childOffset + 4 (n + 1))</c> of <paramref name="file"/>, which must exist at that length when the stage starts; <paramref name="place"/> returns where the LCP and the overflow list go.
    /// </summary>
    /// <exception cref="InvalidDataException">The suffix array is not a permutation of the text positions (wrong range, <c>SA[0] != n</c>, an entry out of range or repeated).</exception>
    internal static LcpChildResult Build(string saPath, long saOffset, MappedTextReader text, MappedOutputFile file, long childOffset, Func<LcpPlacementRequest, LcpPlacement> place, LcpChildOptions options)
    {
        ArgumentNullException.ThrowIfNull(saPath);
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(file);
        ArgumentNullException.ThrowIfNull(place);
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfNegative(saOffset);
        ArgumentOutOfRangeException.ThrowIfNegative(childOffset);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ArenaCapBytes);
        ArgumentOutOfRangeException.ThrowIfNegative(options.ChildArenaCapBytes);
        if (options.Codec is not (ChildCodec.Packed or ChildCodec.Tagged)) throw new ArgumentException("the container holds the child table in CHILD codec 1 (packed) or codec 2 (tagged), nothing else", nameof(options));
        long n = text.Length;
        if (n > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{n} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");
        if (options.Codec == ChildCodec.Tagged && n > TaggedChild.MaxLength) throw new NotSupportedException($"the tagged child table holds ranks in 30 bits (CHILD codec 2); a text of {n} symbols needs the packed codec (CHILD codec 1)");
        if (childOffset > file.Length - (4 * (n + 1))) throw new ArgumentOutOfRangeException(nameof(childOffset), childOffset, $"the child table of {4 * (n + 1)} bytes does not fit the file of {file.Length} bytes at this offset");
        LcpChildParameters p = LcpChildRules.Resolve(n, options);

        string scratch = Path.GetFullPath(options.ScratchDir);
        bool createdRoot = !Directory.Exists(scratch);
        string ph4 = Path.Combine(scratch, "ph4"), ph5 = Path.Combine(scratch, "ph5");
        IoCounters counters = options.Counters ?? new IoCounters();
        var lcpBudget = new ArenaBudget(options.ArenaCapBytes);
        var childBudget = new ArenaBudget(options.ChildArenaCapBytes);
        using var sampler = options.MeasureRss ? new RssSampler(new RssReader(options.ProcSource), 10) : null;
        try
        {
            Directory.CreateDirectory(ph4);
            Directory.CreateDirectory(ph5);
            sampler?.Start();
            var lcpOptions = new ExtLcpOptions
            {
                ScratchDir = ph4,
                ChunkPositions = p.LcpChunk,
                Threads = options.Threads,
                Workers = p.LcpWorkers,
                ReadBufferBytes = p.LcpReadBytes,
                TupleBufferBytes = p.LcpTupleBytes,
                WindowBytes = p.WindowBytes,
                Budget = lcpBudget,
                ArenaGuards = options.ArenaGuards,
                Counters = counters,
                Files = options.Files,
                Trace = options.Trace,
                SaOffset = saOffset,
                Place = place,
            };
            var childOptions = new ChildBuildOptions
            {
                ScratchDir = ph5,
                BlockEntries = p.BlockEntries,
                WindowBytes = p.ChildWindowBytes,
                WindowSlots = p.ChildSlots,
                FlushWindows = true,
                SpillFlushBytes = p.SpillFlushBytes,
                Codec = options.Codec,
                Budget = childBudget,
                ArenaGuards = options.ArenaGuards,
                Counters = counters,
            };
            ExtLcpStats lcp;
            ChildScanResult scan;
            long t0 = Stopwatch.GetTimestamp();
            using (var lazy = new LazySectionScan(n, file, childOffset, childOptions))
            {
                lcp = ExternalLcp.Build(saPath, text, string.Empty, string.Empty, lcpOptions, lazy);
                scan = lazy.Complete();
            }
            double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
            RssStats rss = sampler?.Stop() ?? default;
            double? anonMb = rss.Available ? rss.MaxRssAnonKb * 1024.0 / 1e6 : null;
            RssVerdict verdict = RssVerdict.Unavailable;
            if (anonMb is { } mb) verdict = mb <= options.RssAnonCapMb ? RssVerdict.Within : RssVerdict.Exceeded;
            return new LcpChildResult
            {
                Parameters = p,
                Lcp = lcp,
                Scan = scan,
                Io = counters.Snapshot(),
                Seconds = seconds,
                ArenaPeakBytes = lcpBudget.PeakBytes,
                ChildArenaPeakBytes = childBudget.PeakBytes,
                RssAnonPeakMb = anonMb,
                RssVerdict = verdict,
                Reservation = lcp.Reservation,
            };
        }
        finally
        {
            DeleteDirectoryQuietly(ph4);
            DeleteDirectoryQuietly(ph5);
            if (createdRoot) DeleteEmptyDirectoryQuietly(scratch);
        }
    }

    /// <summary>The child-table scan over a range of the container, created when PH4 delivers its first block (after the width of the LCP is known and the file has grown: the scan holds a writer of the file, which a growth refuses).</summary>
    private sealed class LazySectionScan(long n, MappedOutputFile file, long start, ChildBuildOptions options) : ILcpChunkConsumer, IDisposable
    {
        private ChildScan? _scan;

        public void Consume(long firstRank, ReadOnlySpan<uint> lcp)
        {
            _scan ??= ChildScan.Create(n, file, start, options);
            _scan.Consume(firstRank, lcp);
        }

        internal ChildScanResult Complete() => (_scan ?? throw new InvalidOperationException("PH4 delivered no LCP block")).Complete();

        public void Dispose() => _scan?.Dispose();
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
