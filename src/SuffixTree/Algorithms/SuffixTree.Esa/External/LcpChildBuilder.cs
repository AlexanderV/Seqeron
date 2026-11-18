#pragma warning disable S2325 // stub: the members become instance members of the implementation
namespace SuffixTree.Esa.External;

/// <summary>The files PH4 and PH5 write: the base array of the LCP, its overflow list, the child table (layout A) and, when a path is given, the statistics text.</summary>
/// <param name="Lcp">Base array: <c>n + 1</c> entries of the width the rule of ADR-002 section 4 chooses (little endian).</param>
/// <param name="Overflow">Overflow list: pairs <c>(rank, value)</c> of two 4-byte fields (empty for U32 and when no value escapes).</param>
/// <param name="Child">Child table: <c>n + 1</c> entries of 4 bytes.</param>
/// <param name="Stats">The statistics (<see cref="LcpChildResult.StatsText"/>), or null for none.</param>
internal sealed record LcpChildOutputs(string Lcp, string Overflow, string Child, string? Stats = null);

/// <summary>
/// Explicit values for the sizes that the rule would otherwise derive (null = derive). The tests set every one of them to tiny values (R = 3..50, blocks of 2..5 entries, windows of a few bytes) so that every
/// multi-chunk, spill and window path runs on texts of 10^3-10^5 symbols; the output never depends on any of them.
/// </summary>
internal sealed record LcpChildTuning
{
    /// <summary>PH4: R, positions per Phi chunk and ranks per LCP chunk.</summary>
    internal long? LcpChunk { get; init; }

    /// <summary>PH4: parallel B2 tasks (null = what the arena cap allows at the reference geometry).</summary>
    internal int? LcpWorkers { get; init; }

    /// <summary>PH4: sequential read buffer.</summary>
    internal int? LcpReadBytes { get; init; }

    /// <summary>PH4: write buffer of the tuple distributions.</summary>
    internal int? LcpTupleBytes { get; init; }

    /// <summary>PH4: sliding window of the LCP and overflow writers.</summary>
    internal int? WindowBytes { get; init; }

    /// <summary>PH5: entries per block of the stack.</summary>
    internal int? BlockEntries { get; init; }

    /// <summary>PH5: window of the CHILD writer (a power of two).</summary>
    internal int? ChildWindowBytes { get; init; }

    /// <summary>PH5: windows of the CHILD writer.</summary>
    internal int? ChildSlots { get; init; }

    /// <summary>PH5: fsync of the spill file after this many bytes.</summary>
    internal long? SpillFlushBytes { get; init; }
}

/// <summary>What <see cref="LcpChildBuilder.Build"/> needs besides the files.</summary>
internal sealed record LcpChildOptions
{
    /// <summary>The scratch directory (created when missing; <c>ph4</c> and <c>ph5</c> inside it are created and removed by the build, the directory itself is removed when this build created it and it is empty).</summary>
    internal required string ScratchDir { get; init; }

    /// <summary>Worker threads of B2 (the B1, B3 and PH5 work is sequential).</summary>
    internal int Threads { get; init; } = 4;

    /// <summary>The codec of the child table: the owner's decision of 2026-11-13 is <see cref="ChildCodec.Packed"/>; <see cref="ChildCodec.Tagged"/> holds ranks in 30 bits (texts up to 2^30 - 1 symbols).</summary>
    internal ChildCodec Codec { get; init; } = ChildCodec.Packed;

    /// <summary>
    /// The scale-down emulation of ADR-004 section 5 (E3 to E7): when set (3,000,000,000 for the 100 Mbp runs), R is its reference value for 3x10^9 symbols multiplied by n / ScaleN; 0 = the reference values themselves.
    /// </summary>
    internal long ScaleN { get; init; }

    /// <summary>E7 series A: the n at which the scale-down rule is evaluated (0 = the text length). Ignored when <see cref="ScaleN"/> is 0.</summary>
    internal long RuleN { get; init; }

    /// <summary>Cap of the native buffers of PH4 (owner decision: 160 MiB).</summary>
    internal long ArenaCapBytes { get; init; } = 160L << 20;

    /// <summary>Cap of the native buffers of PH5 (it lives next to PH4 in the fused phase; owner decision: 64 MiB).</summary>
    internal long ChildArenaCapBytes { get; init; } = 64L << 20;

    internal LcpChildTuning Tuning { get; init; } = new();

    /// <summary>The counters of the build (a private instance when null; the result carries a snapshot either way).</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>How the output files are created (reservation of the disk space); null = <see cref="ReservationMode.Auto"/>.</summary>
    internal MappedFileOptions? Files { get; init; }

    /// <summary>Test hook: options of the arenas of PH4 and PH5 (guard bands, poison pattern); null = production arenas.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>The alphabet size, which the statistics text records; required when <see cref="LcpChildOutputs.Stats"/> is set.</summary>
    internal int? Sigma { get; init; }

    /// <summary>fsync every output file when the build is done (not part of the stage time of the prototype).</summary>
    internal bool SyncOutputs { get; init; }

    /// <summary>Sample RssAnon while the build runs (the figure of the gate); unavailable where there is no <c>/proc</c>.</summary>
    internal bool MeasureRss { get; init; } = true;

    /// <summary>Where the kernel's text files come from (null = the real <c>/proc</c>); a test passes a source that fails like a platform without one.</summary>
    internal IProcSource? ProcSource { get; init; }

    /// <summary>The gate of the owner decision of 2026-11-18: RssAnon of the process at most 320 MB (10^6 bytes) at every n with the default parameters and four threads.</summary>
    internal double RssAnonCapMb { get; init; } = 320;

    /// <summary>Test hook, called with <c>"b1"</c>, <c>"b2"</c>, <c>"b3"</c> when a step of PH4 begins; an exception thrown here fails the build like a failing step.</summary>
    internal Action<string>? Trace { get; init; }
}

/// <summary>Every size of one run, resolved from <see cref="LcpChildOptions"/> and the text length (what the report prints, and what the tests pin).</summary>
internal sealed record LcpChildParameters
{
    internal long N { get; init; }

    /// <summary>n / ScaleN (or RuleN / ScaleN), or 1 when unscaled.</summary>
    internal double Factor { get; init; }

    internal bool Scaled { get; init; }

    // PH4
    internal long LcpChunk { get; init; }

    internal int LcpWorkers { get; init; }

    internal int LcpReadBytes { get; init; }

    internal int LcpTupleBytes { get; init; }

    internal int WindowBytes { get; init; }

    // PH5
    internal int BlockEntries { get; init; }

    internal int ChildWindowBytes { get; init; }

    internal int ChildSlots { get; init; }

    internal long SpillFlushBytes { get; init; }
}

/// <summary>The scale-down rule of ADR-004 (E4, E5, E6) for PH4 and PH5 only.</summary>
internal static class LcpChildRules
{
    /// <summary>The reference text length of the sizes: 3x10^9 symbols.</summary>
    internal const long ReferenceN = 3_000_000_000;

    /// <summary>The reference R: 8,000,000 positions per chunk (P = 375 chunks at <see cref="ReferenceN"/>).</summary>
    internal const long ChunkReference = 8_000_000;

    internal static LcpChildParameters Resolve(long n, LcpChildOptions options) => throw new NotImplementedException();
}

/// <summary>The verdict of the RssAnon gate of the builder.</summary>
internal enum RssVerdict
{
    /// <summary>The platform has no <c>/proc</c> (or sampling was off): nothing was measured, the gate was neither passed nor failed.</summary>
    Unavailable,

    /// <summary>The peak RssAnon of the process was at most the cap.</summary>
    Within,

    /// <summary>The peak RssAnon of the process exceeded the cap.</summary>
    Exceeded,
}

/// <summary>What a finished build knows.</summary>
internal sealed class LcpChildResult
{
    internal required LcpChildParameters Parameters { get; init; }

    internal required ExtLcpStats Lcp { get; init; }

    internal required ChildScanResult Scan { get; init; }

    /// <summary>The counters of the build at its end.</summary>
    internal IoSnapshot Io { get; init; }

    /// <summary>Wall seconds of PH4 + PH5 (without the optional fsync).</summary>
    internal double Seconds { get; init; }

    /// <summary>Peak of the native buffers of PH4 and of PH5.</summary>
    internal long ArenaPeakBytes { get; init; }

    internal long ChildArenaPeakBytes { get; init; }

    /// <summary>Peak RssAnon of the whole process in MB (10^6 bytes) while the build ran, or null where it is unavailable.</summary>
    internal double? RssAnonPeakMb { get; init; }

    internal RssVerdict RssVerdict { get; init; }

    /// <summary>What the creation of the LCP base file did about the disk space (the other outputs use the same mode).</summary>
    internal ReservationOutcome Reservation { get; init; }

    /// <summary>The line-oriented statistics text (<c>.stats</c> of the prototype): NODE_COUNT, the internal nodes, LRS_DEPTH and LRS_RANK and the shape of the output.</summary>
    internal string StatsText(int sigma) => throw new NotImplementedException();
}

/// <summary>
/// PH4 (base mode) and PH5 (layout A) of ADR-004 on a given suffix-array file, in one process, fused through the B3 hook exactly as in E6: the compact LCP, its overflow list and the child table are written through
/// bounded mapped windows (every big output is created with <see cref="MappedOutputFile.Create"/>, so a full disk is an <see cref="IOException"/> at creation), nothing is built in RAM and written afterwards.
/// </summary>
internal static class LcpChildBuilder
{
    /// <summary>Builds the outputs of <paramref name="text"/> from the suffix array file <paramref name="saPath"/> (<c>n + 1</c> little-endian <c>uint</c> entries, <c>SA[0] = n</c>).</summary>
    /// <exception cref="InvalidDataException">The suffix array file is not the file of a permutation of the text positions (wrong length, <c>SA[0] != n</c>, an entry out of range or repeated).</exception>
    /// <exception cref="NotSupportedException">The text is longer than <see cref="PackedTextFormat.MaxSymbols"/>, or longer than the tagged codec holds.</exception>
    /// <exception cref="IOException">A file cannot be created, written or spilled to (a full disk included); every output of this build is removed.</exception>
    internal static LcpChildResult Build(string saPath, MappedTextReader text, LcpChildOutputs outputs, LcpChildOptions options) => throw new NotImplementedException();
}
