namespace SuffixTree.Esa.External;

/// <summary>Options of <see cref="BwtPass"/>.</summary>
internal sealed record BwtOptions
{
    /// <summary>Bytes of the suffix array read at a time (rounded down to whole entries, at least one).</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>Window of the mapped writer of the BWT file (at least 1 MiB by default: the .NET view accessor flushes every window when it is disposed, ADR-004 section 5a, E9).</summary>
    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    /// <summary>The shared cap of the native buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Test hook: options of the arena (guard bands, poison pattern); null = production arena.</summary>
    internal ArenaOptions? ArenaGuards { get; init; }

    /// <summary>The suffix-array reads and the BWT writes are counted here (null = not counted).</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>How the BWT file is created (reservation of the disk space); null = <see cref="ReservationMode.Auto"/>.</summary>
    internal MappedFileOptions? Files { get; init; }
}

/// <summary>What the BWT pass found: the text length, the row of the suffix 0 and the symbol counts of the BWT (those of the text).</summary>
internal sealed record BwtResult(long N, long PrimaryRow, long[] Counts);

/// <summary>
/// The BWT of a text from its suffix array and the packed text (PH6 of ADR-004 section 3.3, route (ii) of E9: a separate pass over the suffix array, the text is read at random through its mapping): <c>BWT[j] = T[SA[j] - 1]</c>,
/// one byte per row (<c>n + 1</c> rows, row 0 is the empty suffix, whose symbol is the last symbol of the text); the row of the suffix 0 has no symbol (the sentinel) and holds the byte 0, its number is the primary row.
/// The suffix array is a range of a file (the SA section of the container); the BWT file is created at its final length through <see cref="MappedOutputFile.Create"/> and written sequentially through a mapped window.
/// </summary>
internal static class BwtPass
{
    /// <summary>
    /// Builds the BWT file <paramref name="bwtPath"/> from the <c>n + 1</c> little-endian <c>uint</c> entries at byte <paramref name="saOffset"/> of <paramref name="saPath"/> and the text <paramref name="text"/> (alphabet size
    /// <paramref name="sigma"/>, at most 256). On any failure the BWT file is removed.
    /// </summary>
    /// <exception cref="InvalidDataException">The suffix array range does not fit the file, <c>SA[0] != n</c>, an entry is above n, the suffix 0 occurs twice or not at all, or a symbol is outside the alphabet.</exception>
    internal static BwtResult Build(string saPath, long saOffset, MappedTextReader text, int sigma, string bwtPath, BwtOptions? options = null) => throw new NotImplementedException();
}
