using System.Globalization;
using System.Runtime.InteropServices;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.Esa;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Ph6;

/// <summary>
/// The side file of a BWT file: the text length, the row of the suffix 0 (its BWT symbol is the sentinel, which a byte cannot hold for a 256-symbol alphabet, so the BWT file keeps a dummy 0 there and the row is recorded here),
/// and the symbol counts (the counts of the BWT without the sentinel row are those of the text: PH0's C[] in the pipeline).
/// </summary>
internal sealed record BwtMeta(long N, long PrimaryRow, long[] Counts)
{
    internal int Sigma => Counts.Length;

    internal static string MetaPath(string bwtPath) => bwtPath + ".meta";

    internal string Text() => string.Create(CultureInfo.InvariantCulture, $"n={N}\nprimary={PrimaryRow}\nsigma={Sigma}\ncounts={string.Join(',', Counts)}\n");

    internal static BwtMeta Parse(string text)
    {
        var kv = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (string line in text.Split('\n', StringSplitOptions.RemoveEmptyEntries))
        {
            int eq = line.IndexOf('=', StringComparison.Ordinal);
            if (eq > 0) kv[line[..eq]] = line[(eq + 1)..].Trim();
        }
        string Need(string key) => kv.TryGetValue(key, out string? v) ? v : throw new InvalidDataException($"the BWT meta lacks '{key}='");
        long n = long.Parse(Need("n"), CultureInfo.InvariantCulture);
        long primary = long.Parse(Need("primary"), CultureInfo.InvariantCulture);
        int sigma = int.Parse(Need("sigma"), CultureInfo.InvariantCulture);
        string list = Need("counts");
        long[] counts = list.Length == 0 ? [] : [.. list.Split(',').Select(c => long.Parse(c, CultureInfo.InvariantCulture))];
        if (counts.Length != sigma) throw new InvalidDataException($"the BWT meta says sigma={sigma} but lists {counts.Length} counts");
        return new BwtMeta(n, primary, counts);
    }

    internal static BwtMeta Read(string bwtPath) => Parse(File.ReadAllText(MetaPath(bwtPath)));
}

internal sealed record BwtOptions
{
    /// <summary>Bytes of the suffix array read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>Window of the mapped writer of the BWT file.</summary>
    internal int WindowBytes { get; init; } = 1 << 20;

    internal ArenaBudget? Budget { get; init; }
}

/// <summary>
/// The BWT of the text from its suffix array file and the packed text (route (ii) of E9: a separate pass over the SA, the text is read at random through its mapping): <c>BWT[j] = T[SA[j] - 1]</c>, one byte per row
/// (<c>n + 1</c> rows, row 0 is the empty suffix, its symbol is the last symbol of the text); the row of the suffix 0 holds the dummy byte 0 and is recorded in the meta file together with the symbol counts.
/// Reads the SA once (4 B per row, sequentially), writes the BWT once (1 B per row, sequentially through a mapped window); the text accesses are random reads of its mapping.
/// </summary>
internal static class BwtPass
{
    internal static BwtMeta Build(string saPath, PackedText text, int sigma, string bwtPath, BwtOptions? options = null)
    {
        BwtOptions o = options ?? new BwtOptions();
        ArgumentNullException.ThrowIfNull(text);
        ArgumentOutOfRangeException.ThrowIfNegative(sigma);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(sigma, 256); // the BWT is one byte per row
        ArgumentOutOfRangeException.ThrowIfLessThan(o.ReadBytes, 4);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.WindowBytes, 1);
        long n = text.Length;
        if (n >= uint.MaxValue - 2) throw new NotSupportedException("positions are 32-bit in the prototype");
        long rows = n + 1;
        using var sa = new CountedFile(saPath, FileMode.Open, FileAccess.Read);
        if (sa.Length != 4 * rows) throw new InvalidDataException($"the suffix array file has {sa.Length} bytes, expected {4 * rows} for n = {n}");

        var counts = new long[sigma];
        long primary = -1;
        int readEntries = (int)Math.Min(o.ReadBytes / 4, rows);
        using var arena = new BudgetedArena(o.Budget ?? new ArenaBudget(long.MaxValue));
        Buf<uint> saBuf = arena.Allocate<uint>(readEntries, zero: false, "bwt.sa");
        Buf<byte> outBuf = arena.Allocate<byte>(readEntries, zero: false, "bwt.out");
        using (var writer = new MappedWindowWriter(bwtPath, rows, o.WindowBytes))
        {
            for (long first = 0; first < rows; first += readEntries)
            {
                int take = (int)Math.Min(readEntries, rows - first);
                Span<uint> entries = saBuf.AsSpan(take);
                sa.ReadExactAt(MemoryMarshal.AsBytes(entries), first * 4);
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
        var meta = new BwtMeta(n, primary, counts);
        File.WriteAllText(BwtMeta.MetaPath(bwtPath), meta.Text());
        return meta;
    }
}
