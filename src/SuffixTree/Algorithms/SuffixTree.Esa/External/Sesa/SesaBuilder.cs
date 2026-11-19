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

    /// <summary>Sample RssAnon while the build runs (the figure of the gate); unavailable where there is no <c>/proc</c>.</summary>
    internal bool MeasureRss { get; init; } = true;

    /// <summary>Where the kernel's text files come from (null = the real <c>/proc</c>); a test passes a source that fails like a platform without one.</summary>
    internal IProcSource? ProcSource { get; init; }

    /// <summary>The gate of the owner decision of 2026-11-18: RssAnon of the process at most 320 MB (10^6 bytes) at every n with the default parameters and four threads.</summary>
    internal double RssAnonCapMb { get; init; } = 320;

    /// <summary>Test hook, called with a name of <see cref="SesaHooks"/> at every boundary of the build; an exception thrown here fails the build like a failing step, a call that blocks holds the build at that point.</summary>
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

    /// <summary>What the growth of the file behind the reserved length did (null when it did not grow).</summary>
    internal ReservationOutcome? GrowthReservation { get; init; }

    /// <summary>Bytes the file grew by after PH4-B2 (0 when it did not).</summary>
    internal long GrownBytes { get; init; }

    internal required IReadOnlyList<SesaCrcRow> CrcRows { get; init; }
}

/// <summary>
/// The end-to-end builder of the SESA container (work package 1b.7e, ADR-004 sections 3-4, ADR-003 as amended on 2026-11-18): PH0 (the packed text, the alphabet, the symbol counts and the SHA-256) -> PH1-PH3 (the suffix
/// array) -> PH4-PH5 (LCP, child table) -> PH6 (BWT and PSI, optional) -> PH7 (the section checksums, the header, the atomic publication). Every big output is built directly into its section of
/// <c>&lt;target&gt;.tmp</c>, created through <see cref="MappedOutputFile.Create"/>; scratch lives under <c>&lt;target&gt;.tmp.d</c>. The target appears only by the rename at the end; a failure removes the temporary file and
/// the scratch and leaves no target (an existing target stays as it was). The output depends on the input and on <see cref="SesaBuildOptions.IncludePsi"/> only.
/// </summary>
internal static class SesaBuilder
{
    /// <summary>Builds the container of the text in the file <paramref name="inputPath"/> (one byte per code unit) at <paramref name="targetPath"/>.</summary>
    /// <exception cref="DirectoryNotFoundException">The directory of the target does not exist.</exception>
    /// <exception cref="NotSupportedException">The text is longer than <see cref="PackedTextFormat.MaxSymbols"/>.</exception>
    /// <exception cref="IOException">A file cannot be created or written (a full disk included); nothing of the build is left.</exception>
    /// <exception cref="InvalidOperationException">The input changed between the two passes.</exception>
    internal static SesaBuildResult Build(string inputPath, string targetPath, SesaBuildOptions? options = null) => throw new NotImplementedException();
}
