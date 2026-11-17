namespace SuffixTree.EmProto.Verify;

/// <summary>Outcome of one check of <see cref="OutputVerifier"/>. <see cref="Skipped"/> means the check could not run (a file it needs is unusable, or its stack bound was exceeded): it is never a pass.</summary>
internal enum CheckStatus
{
    Pass,
    Fail,
    Skipped,
}

internal sealed record CheckResult(string Name, CheckStatus Status, string Detail, double Seconds);

/// <summary>What to verify and how much memory and time the check may use (ADR-004 section 5, row E8).</summary>
internal sealed record VerifyOptions
{
    /// <summary>The output prefix of <c>e6</c>: the files are <c>.pk .pk.meta .sa .lcp .ovf .cld .stats</c> after it.</summary>
    internal required string Prefix { get; init; }

    /// <summary>The ASCII input text (one byte per symbol); when given, the packed text is compared with it symbol by symbol and its SHA-256 with the one in the metadata.</summary>
    internal string? InputPath { get; init; }

    /// <summary>Where the permutation bitmap (n + 1 bits, a temporary file) lives; default: next to the outputs.</summary>
    internal string? ScratchDir { get; init; }

    internal int Threads { get; init; } = 4;

    /// <summary>
    /// Common-prefix lengths of adjacent suffixes are computed directly up to this many symbols (0 = no limit). A pair whose stored LCP is longer ("long pair": the quadratic cost of N runs and of a^n) is checked up to the cap,
    /// at its mismatch position (both symbols, hence the order) and at <see cref="Probes"/> further offsets of its matching part, and counted in the fact <c>long_pairs</c>.
    /// </summary>
    internal long LceCap { get; init; } = 65536;

    /// <summary>Longest stack (16 bytes per entry) the scan of the child table may use; a text with a longer LCP than this is reported as skipped, never as passed.</summary>
    internal long StackCapEntries { get; init; } = 1L << 26;

    /// <summary>SHA-256 of every output file (the identity of two runs).</summary>
    internal bool Hash { get; init; } = true;

    /// <summary>The scan of the LCP stream: statistics and the child table cell by cell.</summary>
    internal bool Child { get; init; } = true;

    /// <summary>Ranks per unit of work of the parallel pass (tests use a few to cross block boundaries).</summary>
    internal int BlockRanks { get; init; } = 1 << 20;

    /// <summary>Probes of the matching part of a long pair.</summary>
    internal int Probes { get; init; } = 8;
}

internal sealed class VerifyResult
{
    internal List<CheckResult> Checks { get; } = [];

    /// <summary>Hex SHA-256 of every output file by extension (<c>pk</c>, <c>pk.meta</c>, <c>sa</c>, <c>lcp</c>, <c>ovf</c>, <c>cld</c>, <c>stats</c>).</summary>
    internal Dictionary<string, string> Sha256 { get; } = [];

    /// <summary>Numbers of the run (n, widths, long pairs, ...) for the rows.</summary>
    internal Dictionary<string, string> Facts { get; } = [];

    /// <summary>Every listed check passed (none failed, none skipped).</summary>
    internal bool AllPassed => Checks.Count > 0 && Checks.All(c => c.Status == CheckStatus.Pass);

    internal CheckResult? Find(string name) => Checks.FirstOrDefault(c => c.Name == name);

    internal bool Passed(string name) => Find(name)?.Status == CheckStatus.Pass;

    internal bool Failed(string name) => Find(name)?.Status == CheckStatus.Fail;
}
