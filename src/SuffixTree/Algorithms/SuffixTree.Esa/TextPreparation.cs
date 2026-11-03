using System.Numerics;

namespace SuffixTree.Esa;

/// <summary>Distinct-code-unit scan and order-preserving dense remap of a UTF-16 text (ADR-002 section 2 item 2).</summary>
internal static class TextPreparation
{
    private const int MapWords = 1024; // 65,536 bits

    /// <summary>Number of distinct UTF-16 code units in <paramref name="text"/> (one pass over a 65,536-bit map).</summary>
    internal static int CountDistinct(ReadOnlySpan<char> text)
    {
        Span<ulong> map = stackalloc ulong[MapWords];
        map.Clear();
        MarkUnits(text, map);
        int count = 0;
        foreach (ulong w in map) count += BitOperations.PopCount(w);
        return count;
    }

    /// <summary>
    /// Scans <paramref name="text"/> for its distinct code units and replaces every unit by its rank among them. Lone surrogates,
    /// U+0000 and U+FFFF are ordinary code units; the order is the numeric (ordinal) order, never culture-dependent.
    /// </summary>
    /// <exception cref="NotSupportedException">The suffix array would not fit the 32-bit path.</exception>
    internal static PreparedText Prepare(ReadOnlySpan<char> text)
    {
        EsaLimits.EnsureTextLength(text.Length);
        Span<ulong> map = stackalloc ulong[MapWords];
        map.Clear();
        MarkUnits(text, map);

        // prefix[w] = number of marked units below word w, so rank(c) = prefix[c >> 6] + popcount(map[c >> 6] & below(c))
        Span<int> prefix = stackalloc int[MapWords];
        int sigma = 0;
        for (int w = 0; w < MapWords; w++)
        {
            prefix[w] = sigma;
            sigma += BitOperations.PopCount(map[w]);
        }

        var alphabet = new char[sigma];
        int k = 0;
        for (int w = 0; w < MapWords; w++)
        {
            ulong bits = map[w];
            while (bits != 0)
            {
                alphabet[k++] = (char)((w << 6) + BitOperations.TrailingZeroCount(bits));
                bits &= bits - 1;
            }
        }

        if (sigma <= 256)
        {
            var symbols = new byte[text.Length];
            for (int i = 0; i < text.Length; i++) symbols[i] = (byte)Rank(text[i], map, prefix);
            return new PreparedText(alphabet, symbols, null);
        }
        else
        {
            var symbols = new ushort[text.Length];
            for (int i = 0; i < text.Length; i++) symbols[i] = (ushort)Rank(text[i], map, prefix);
            return new PreparedText(alphabet, null, symbols);
        }
    }

    private static void MarkUnits(ReadOnlySpan<char> text, Span<ulong> map)
    {
        foreach (char c in text) map[c >> 6] |= 1UL << (c & 63);
    }

    private static int Rank(char c, ReadOnlySpan<ulong> map, ReadOnlySpan<int> prefix)
        => prefix[c >> 6] + BitOperations.PopCount(map[c >> 6] & ((1UL << (c & 63)) - 1));
}
