using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Ph6;

internal sealed record PsiOptions
{
    /// <summary>Bytes of one output window (one per symbol of the pass); rounded down to whole 4-byte entries, at least one entry.</summary>
    internal int WindowBytes { get; init; } = 64 * 1024;

    /// <summary>Symbols whose output cursors are open in one pass (the ADR: 64); a larger alphabet takes ceil(sigma / this) passes over the BWT.</summary>
    internal int SymbolsPerPass { get; init; } = 64;

    /// <summary>Bytes of the BWT read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>PSI[0]: false = 0 (the convention of ADR-003 / ADR-004), true = the row of the suffix 0 (the cyclic Psi, a single cycle through all n + 1 rows).</summary>
    internal bool CyclicEntryZero { get; init; }

    /// <summary>Cap on (symbols per pass) x (window bytes), the mapped working set of the pass, charged to the budget like an arena block (ADR-004: arena cap 160 MiB).</summary>
    internal long WindowCapBytes { get; init; } = 160L << 20;

    internal ArenaBudget? Budget { get; init; }
}

internal sealed record PsiStats(int Passes, int SymbolsPerPass, long Windows, long MappedBoundBytes, long BwtBytesRead, long PsiBytesWritten);

/// <summary>PH6 of ADR-004 section 3.3: Psi from the BWT by sigma-way distribution, without the inverse permutation.</summary>
internal static class PsiStream
{
    internal static PsiStats Build(string bwtPath, BwtMeta meta, string psiPath, PsiOptions? options = null) => throw new NotImplementedException();
}
