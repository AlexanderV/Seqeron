namespace SuffixTree.Esa.External;

/// <summary>Options of <see cref="PsiStream"/>.</summary>
internal sealed record PsiOptions
{
    /// <summary>Bytes of one output window (one per symbol of a pass); at least 1 MiB by default (the .NET view accessor flushes every window when it is disposed, ADR-004 section 5a, E9); rounded down to whole 4-byte entries, at least one entry.</summary>
    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    /// <summary>Symbols whose output cursors are open in one pass (ADR-004 section 3.3: 64); a larger alphabet takes <c>ceil(sigma / this)</c> passes over the BWT.</summary>
    internal int SymbolsPerPass { get; init; } = 64;

    /// <summary>Bytes of the BWT read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>Cap on (symbols per pass) x (window bytes), the mapped working set of a pass, charged to the budget like an arena block (ADR-004: arena cap 160 MiB).</summary>
    internal long WindowCapBytes { get; init; } = 160L << 20;

    /// <summary>The shared cap of the native buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Test hook: options of the arena (guard bands, poison pattern); null = production arena.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>The BWT reads and the mapped writes are counted here (null = not counted).</summary>
    internal IoCounters? Counters { get; init; }
}

/// <summary>What one PSI build did: passes over the BWT, symbols per pass, windows mapped, the bound of the mapped windows, BWT bytes read and PSI bytes written.</summary>
internal sealed record PsiStats(int Passes, int SymbolsPerPass, long Windows, long MappedBoundBytes, long BwtBytesRead, long PsiBytesWritten);

/// <summary>
/// PH6 of ADR-004 section 3.3: Psi from the BWT by sigma-way distribution, without the inverse permutation, written into the PSI section of the container. With <c>BWT[j] = T[SA[j] - 1]</c> and <c>C[c]</c> = the number
/// of suffixes that start with a symbol below c (the empty suffix included), <c>Psi[C[c] + k]</c> is the row of the k-th (0-based) occurrence of c in the BWT: one scan of the BWT in row order with one output cursor per
/// symbol, every cursor writing its own contiguous range of the section strictly sequentially through its own mapped window. At most <see cref="PsiOptions.SymbolsPerPass"/> cursors are open at a time. The row of the
/// suffix 0 has no symbol and is skipped; <c>PSI[0]</c> = 0 by convention. The result is <c>PSI[i] = ISA[SA[i] + 1]</c> for <c>SA[i] &lt; n</c> and 0 for <c>SA[i] = n</c>.
/// </summary>
internal static class PsiStream
{
    /// <summary>
    /// Writes the <c>n + 1</c> little-endian <c>uint</c> entries of Psi into the range <c>[psiOffset, psiOffset + 4 (n + 1))</c> of <paramref name="file"/> (every byte of the range is written) from the BWT file
    /// <paramref name="bwtPath"/> (<c>n + 1</c> bytes), the row of the suffix 0 and the symbol counts of the text. The writers are disposed when it returns; the file is not completed.
    /// </summary>
    /// <exception cref="InvalidDataException">The BWT file has the wrong length, or a symbol occurs more or less often than its count.</exception>
    /// <exception cref="ArgumentException">The counts are negative or do not sum to n.</exception>
    internal static PsiStats Build(string bwtPath, long n, long primaryRow, long[] counts, MappedOutputFile file, long psiOffset, PsiOptions? options = null) => throw new NotImplementedException();
}
