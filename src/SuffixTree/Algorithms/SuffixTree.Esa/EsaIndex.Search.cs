using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

internal sealed partial class EsaIndex
{
    private const int StackPatternLimit = 512;

    /// <summary><c>_symbolOf[c]</c> is the dense symbol of the code unit <c>c</c>, or -1 when <c>c</c> does not occur in the text (the array ends at the largest code unit of the text).</summary>
    private int[]? _symbolOf;

    private int[] SymbolOf => _symbolOf ??= BuildSymbolMap(Prepared.Alphabet);

    private static int[] BuildSymbolMap(char[] alphabet)
    {
        if (alphabet.Length == 0) return [];
        var map = new int[alphabet[^1] + 1];
        Array.Fill(map, -1);
        for (int s = 0; s < alphabet.Length; s++) map[alphabet[s]] = s;
        return map;
    }

    /// <summary>Whether <paramref name="pattern"/> is a substring of the text (see the span overload).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="pattern"/> is null.</exception>
    internal bool Contains(string pattern)
    {
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return Contains(pattern.AsSpan());
    }

    /// <summary>The number of occurrences (see the span overload).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="pattern"/> is null.</exception>
    internal int CountOccurrences(string pattern)
    {
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return CountOccurrences(pattern.AsSpan());
    }

    /// <summary>All start positions (see the span overload).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="pattern"/> is null.</exception>
    internal IReadOnlyList<int> FindAllOccurrences(string pattern)
    {
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return FindAllOccurrences(pattern.AsSpan());
    }

    /// <summary>
    /// Whether <paramref name="pattern"/> is a substring of the text. The empty pattern is a substring of every text, the empty one included
    /// (<c>PersistentSuffixTree.Search.cs:10,18</c>).
    /// </summary>
    internal bool Contains(ReadOnlySpan<char> pattern) => pattern.IsEmpty || TryMatch(pattern, out _, out _);

    /// <summary>The number of occurrences of <paramref name="pattern"/>; the empty pattern counts n (<c>PersistentSuffixTree.Search.cs:35-36</c>).</summary>
    internal int CountOccurrences(ReadOnlySpan<char> pattern)
    {
        if (pattern.IsEmpty) return Length;
        return TryMatch(pattern, out int lo, out int hi) ? hi - lo + 1 : 0;
    }

    /// <summary>
    /// All start positions of <paramref name="pattern"/> in suffix-array order (not ascending, allowed by the contract); the empty pattern returns
    /// <c>[0 .. n-1]</c> (<c>PersistentSuffixTree.Search.cs:58-60</c>).
    /// </summary>
    internal IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
    {
        if (pattern.IsEmpty)
        {
            var all = new int[Length];
            for (int i = 0; i < all.Length; i++) all[i] = i;
            return all;
        }
        if (!TryMatch(pattern, out int lo, out int hi)) return Array.Empty<int>();
        var result = new int[hi - lo + 1];
        Array.Copy(SuffixArray, lo, result, 0, result.Length);
        return result;
    }

    /// <summary>
    /// Lazy locate of the ranks <c>[lo..hi]</c>: writes up to <paramref name="positions"/>.Length start positions after skipping <paramref name="skip"/> of them in
    /// suffix-array order and returns the number written. The interval must be a valid rank interval of non-empty suffixes (<c>1 &lt;= lo &lt;= hi &lt;= n</c>), as returned by
    /// <see cref="TryMatch"/>.
    /// </summary>
    internal int Locate(int lo, int hi, long skip, Span<int> positions)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(skip);
        if ((uint)lo > (uint)hi || hi > Length) throw new ArgumentOutOfRangeException(nameof(lo), "not a rank interval of this index");
        long available = (long)hi - lo + 1 - skip;
        if (available <= 0) return 0;
        int count = (int)Math.Min(available, positions.Length);
        SuffixArray.AsSpan((int)(lo + skip), count).CopyTo(positions);
        return count;
    }

    /// <summary>The same for 64-bit positions (the index contract).</summary>
    internal int Locate(int lo, int hi, long skip, Span<long> positions)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(skip);
        if ((uint)lo > (uint)hi || hi > Length) throw new ArgumentOutOfRangeException(nameof(lo), "not a rank interval of this index");
        long available = (long)hi - lo + 1 - skip;
        if (available <= 0) return 0;
        int count = (int)Math.Min(available, positions.Length);
        int start = (int)(lo + skip);
        for (int i = 0; i < count; i++) positions[i] = SuffixArray[start + i];
        return count;
    }

    /// <summary>
    /// Finds the suffix-array interval of the NON-EMPTY <paramref name="pattern"/>: the ranks <c>[lo..hi]</c> of the suffixes that start with it (never rank 0).
    /// False when the pattern does not occur.
    /// </summary>
    internal bool TryMatch(ReadOnlySpan<char> pattern, out int lo, out int hi)
    {
        lo = 0;
        hi = -1;
        if (pattern.IsEmpty) throw new ArgumentException("the empty pattern has no rank interval; handle it before calling", nameof(pattern));
        if (Length == 0) return false;
        return Prepared.Width == SymbolWidth.Byte
            ? Match<byte>(Prepared.Bytes, pattern, out lo, out hi)
            : Match<ushort>(Prepared.UShorts, pattern, out lo, out hi);
    }

    private bool Match<TSym>(ReadOnlySpan<TSym> text, ReadOnlySpan<char> pattern, out int lo, out int hi)
        where TSym : unmanaged
    {
        lo = 0;
        hi = -1;
        int[] map = SymbolOf;
        Span<TSym> pat = pattern.Length <= StackPatternLimit ? stackalloc TSym[pattern.Length] : new TSym[pattern.Length];
        for (int i = 0; i < pattern.Length; i++)
        {
            char c = pattern[i];
            int s = c < map.Length ? map[c] : -1;
            if (s < 0) return false;
            pat[i] = FromInt<TSym>(s);
        }

        return Child.Kind switch
        {
            ChildLayoutKind.Reference3 => Descend(Child.Reference3, SuffixArray, Lcp, text, pat, out lo, out hi),
            ChildLayoutKind.Packed => Descend(Child.Packed, SuffixArray, Lcp, text, pat, out lo, out hi),
            _ => Descend(Child.Tagged, SuffixArray, Lcp, text, pat, out lo, out hi),
        };
    }

    /// <summary>
    /// Pattern descent through the lcp-interval tree (Abouelhoda, Kurtz, Ohlebusch 2004): at an interval of value <c>ell</c> the child that continues with the next pattern symbol is
    /// found by walking the l-indices (the children are in ascending key order, the terminator first), its edge is compared against the text at the first suffix of the child, and
    /// the walk continues below it. Time O(m + sum of the child counts on the path). The pattern symbols are dense symbols that are present in the text.
    /// </summary>
    internal static bool Descend<TChild, TSym>(TChild child, int[] sa, CompactLcp lcp, ReadOnlySpan<TSym> text, ReadOnlySpan<TSym> pat, out int lo, out int hi)
        where TChild : struct, IChildTable
        where TSym : unmanaged
    {
        int n = text.Length;
        int m = pat.Length;
        int lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        while (c < m)
        {
            // c == ell: choose the child whose key equals pat[c]
            int symbol = ToInt(pat[c]);
            int a = lb;
            int boundary = child.FirstLIndex(lb, rb, isLast);
            int b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                int p = sa[a] + ell;
                int key = p < n ? ToInt(text[p]) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0)
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                a = boundary;
                boundary = child.NextLIndex(boundary, rb);
            }

            bool childIsLast = boundary < 0;
            if (a == b)
            {
                // a leaf: the rest of the pattern must be a prefix of the rest of the suffix
                int start = sa[a] + c;
                if (m - c > n - start || !text.Slice(start, m - c).SequenceEqual(pat.Slice(c)))
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                lo = a;
                hi = a;
                return true;
            }

            int childEll = lcp[child.FirstLIndex(a, b, childIsLast)];
            int limit = Math.Min(childEll, m);
            int from = sa[a] + c;
            if (limit > c + 1 && !text.Slice(from + 1, limit - c - 1).SequenceEqual(pat.Slice(c + 1, limit - c - 1)))
            {
                lo = 0;
                hi = -1;
                return false;
            }
            c = limit;
            lb = a;
            rb = b;
            ell = childEll;
            isLast = childIsLast;
        }
        lo = lb;
        hi = rb;
        return true;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int ToInt<TSym>(TSym value)
        where TSym : unmanaged
        => typeof(TSym) == typeof(byte) ? Unsafe.As<TSym, byte>(ref value) : Unsafe.As<TSym, ushort>(ref value);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static TSym FromInt<TSym>(int value)
        where TSym : unmanaged
    {
        TSym result = default;
        if (typeof(TSym) == typeof(byte)) Unsafe.As<TSym, byte>(ref result) = (byte)value;
        else Unsafe.As<TSym, ushort>(ref result) = (ushort)value;
        return result;
    }
}
