using System.Runtime.InteropServices;

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
    internal static BwtResult Build(string saPath, long saOffset, MappedTextReader text, int sigma, string bwtPath, BwtOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(saPath);
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(bwtPath);
        options ??= new BwtOptions();
        ArgumentOutOfRangeException.ThrowIfNegative(saOffset);
        ArgumentOutOfRangeException.ThrowIfNegative(sigma);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(sigma, 256); // the BWT is one byte per row
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ReadBytes, sizeof(uint));
        ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 1);
        long n = text.Length;
        if (n > PackedTextFormat.MaxSymbols) throw new NotSupportedException($"{n} symbols: positions are 32-bit, the limit is {PackedTextFormat.MaxSymbols}");
        if (text.Bits > PackedTextFormat.MaxBits) throw new ArgumentException($"a text of {text.Bits} bits per symbol is not a packed text of the file codec", nameof(text));
        long rows = n + 1;
        using var sa = new CountedFile(saPath, FileMode.Open, FileAccess.Read, options.Counters);
        if (saOffset > sa.Length - (sizeof(uint) * rows)) throw new InvalidDataException($"the suffix array needs {sizeof(uint) * rows} bytes at offset {saOffset}, the file has {sa.Length} bytes");

        var counts = new long[sigma];
        long primary = -1;
        int readEntries = (int)Math.Min(options.ReadBytes / sizeof(uint), rows);
        using var arena = new BudgetedArena(options.Budget ?? new ArenaBudget(long.MaxValue), options.ArenaGuards);
        Buf<uint> saBuf = arena.Allocate<uint>(readEntries, zero: false, "bwt.sa");
        Buf<byte> outBuf = arena.Allocate<byte>(readEntries, zero: false, "bwt.out");
        MappedOutputFile file = MappedOutputFile.Create(bwtPath, rows, options.Files);
        try
        {
            using (file)
            {
                using (MappedWindowWriter writer = file.CreateWriter(options.WindowBytes, options.Counters))
                {
                    for (long first = 0; first < rows; first += readEntries)
                    {
                        int take = (int)Math.Min(readEntries, rows - first);
                        Span<uint> entries = saBuf.AsSpan(take);
                        sa.ReadExactAt(MemoryMarshal.AsBytes(entries), saOffset + (first * sizeof(uint)));
                        Span<byte> output = outBuf.AsSpan(take);
                        for (int i = 0; i < take; i++)
                        {
                            uint p = entries[i];
                            if (p > n) throw new InvalidDataException($"SA[{first + i}] = {p} lies beyond n = {n}");
                            if (first + i == 0 && p != n) throw new InvalidDataException($"SA[0] = {p}, expected the empty suffix n = {n}");
                            if (p == 0)
                            {
                                if (primary >= 0) throw new InvalidDataException($"rows {primary} and {first + i} both hold the suffix 0");
                                primary = first + i;
                                output[i] = 0;
                                continue;
                            }
                            uint s = text.Symbol(p - 1);
                            if (s >= sigma) throw new InvalidDataException($"the symbol {s} before suffix {p} is outside the alphabet of {sigma} symbols");
                            output[i] = (byte)s;
                            counts[s]++;
                        }
                        writer.Append(output);
                    }
                }
                if (primary < 0) throw new InvalidDataException("no row holds the suffix 0");
                file.Complete(fsync: false);
            }
        }
        catch
        {
            DeleteQuietly(bwtPath);
            throw;
        }
        return new BwtResult(n, primary, counts);
    }

    private static void DeleteQuietly(string path)
    {
        try
        {
            File.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // the failure of the pass is what the caller needs
        }
    }
}
