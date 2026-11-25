using System.Runtime.CompilerServices;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa;

internal sealed partial class EsaIndex
{
    private const int StackPatternLimit = 512;

    /// <summary><c>_symbolOf[c]</c> is the dense symbol of the code unit <c>c</c>, or -1 when <c>c</c> does not occur in the text (the array ends at the largest code unit of the text).</summary>
    private int[]? _symbolOf;

    private int[] SymbolOf => _symbolOf ??= BuildSymbolMap(IsMapped ? TextView.Alphabet : Prepared.Alphabet);

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
    internal bool Contains(ReadOnlySpan<char> pattern) => pattern.IsEmpty || TryMatch(pattern, out long _, out long _);

    /// <summary>The number of occurrences of <paramref name="pattern"/>; the empty pattern counts n (<c>PersistentSuffixTree.Search.cs:35-36</c>).</summary>
    /// <exception cref="NotSupportedException">The count does not fit an <see cref="int"/>: <see cref="LongCount(ReadOnlySpan{char})"/> is the 64-bit member.</exception>
    internal int CountOccurrences(ReadOnlySpan<char> pattern)
    {
        if (pattern.IsEmpty) return Length;
        return TryMatch(pattern, out int lo, out int hi) ? hi - lo + 1 : 0;
    }

    /// <summary>The number of occurrences of <paramref name="pattern"/> at 64 bits; the empty pattern counts n.</summary>
    internal long LongCount(ReadOnlySpan<char> pattern)
    {
        if (pattern.IsEmpty) return _length;
        return TryMatch(pattern, out long lo, out long hi) ? hi - lo + 1 : 0;
    }

    /// <summary>The 64-bit count (see the span overload).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="pattern"/> is null.</exception>
    internal long LongCount(string pattern)
    {
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return LongCount(pattern.AsSpan());
    }

    /// <summary>
    /// All start positions of <paramref name="pattern"/> in suffix-array order (not ascending, allowed by the contract); the empty pattern returns
    /// <c>[0 .. n-1]</c> (<c>PersistentSuffixTree.Search.cs:58-60</c>).
    /// </summary>
    /// <exception cref="NotSupportedException">The list does not fit an array of <see cref="int"/> positions: <see cref="Locate(long, long, long, Span{long})"/> pages through it.</exception>
    internal IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
    {
        if (pattern.IsEmpty)
        {
            var all = new int[Length];
            for (int i = 0; i < all.Length; i++) all[i] = i;
            return all;
        }
        if (!TryMatch(pattern, out int lo, out int hi)) return Array.Empty<int>();
        int count = hi - lo + 1;
        return count < SaView.UninitializedFrom ? Sa.ToArray(lo, count) : Sa.ToArrayUninitialized(lo, count);
    }

    /// <summary>
    /// The page of the ranks <c>[lo..hi]</c> that starts <paramref name="skip"/> suffixes in and holds at most <paramref name="capacity"/> of them: the first rank of the page and how many it has (0 when <paramref name="skip"/> is
    /// beyond the interval). The interval must be a valid rank interval of the text of <paramref name="length"/> symbols (<c>0 &lt;= lo &lt;= hi &lt;= n</c>; the callers use <c>1 &lt;= lo</c>: rank 0 is the empty suffix).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="skip"/> is negative, or the interval is not a rank interval of the text.</exception>
    internal static int PageOf(long length, long lo, long hi, long skip, int capacity, out long start)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(skip);
        if ((ulong)lo > (ulong)hi || hi > length) throw new ArgumentOutOfRangeException(nameof(lo), "not a rank interval of this index");
        start = lo;
        long available = hi - lo + 1 - skip;
        if (available <= 0) return 0;
        start = lo + skip;
        return available < capacity ? (int)available : capacity; // (int): below the capacity, which is an int
    }

    /// <summary>
    /// Lazy locate of the ranks <c>[lo..hi]</c>: writes up to <paramref name="positions"/>.Length start positions after skipping <paramref name="skip"/> of them in
    /// suffix-array order and returns the number written. The interval must be a valid rank interval of non-empty suffixes (<c>1 &lt;= lo &lt;= hi &lt;= n</c>), as returned by
    /// <see cref="TryMatch(ReadOnlySpan{char}, out long, out long)"/>.
    /// </summary>
    /// <exception cref="InvalidDataException">A position of the page is not a 32-bit number (a text above 2^31 symbols: use the 64-bit overload).</exception>
    internal int Locate(long lo, long hi, long skip, Span<int> positions)
    {
        int count = PageOf(_length, lo, hi, skip, positions.Length, out long start);
        if (count > 0) Sa.CopyTo(start, positions[..count]);
        return count;
    }

    /// <summary>The same for 64-bit positions (the index contract): the page may start above rank 2^31 and hold positions above 2^31.</summary>
    internal int Locate(long lo, long hi, long skip, Span<long> positions)
    {
        int count = PageOf(_length, lo, hi, skip, positions.Length, out long start);
        if (count > 0) Sa.CopyTo(start, positions[..count]);
        return count;
    }

    /// <summary>
    /// Finds the suffix-array interval of the NON-EMPTY <paramref name="pattern"/>: the ranks <c>[lo..hi]</c> of the suffixes that start with it (never rank 0).
    /// False when the pattern does not occur.
    /// </summary>
    /// <exception cref="NotSupportedException">A rank is above <see cref="int.MaxValue"/> (a text of 2^31 symbols and more): the overload with 64-bit ranks answers it.</exception>
    internal bool TryMatch(ReadOnlySpan<char> pattern, out int lo, out int hi)
    {
        lo = 0;
        hi = -1;
        if (pattern.IsEmpty) throw new ArgumentException("the empty pattern has no rank interval; handle it before calling", nameof(pattern));
        if (_length == 0) return false;
        if (IsMapped) return MatchMapped(pattern, out lo, out hi);
        return Prepared.Width == SymbolWidth.Byte
            ? Match<byte>(Prepared.Bytes, pattern, out lo, out hi)
            : Match<ushort>(Prepared.UShorts, pattern, out lo, out hi);
    }

    /// <summary>
    /// <see cref="TryMatch(ReadOnlySpan{char}, out int, out int)"/> at 64-bit ranks: the ranks of the suffixes that start with the NON-EMPTY <paramref name="pattern"/> (never rank 0), over any index: the managed arrays and the
    /// narrow views answer with the numbers of the 32-bit descent, widened; the other containers with the long descent.
    /// </summary>
    internal bool TryMatch(ReadOnlySpan<char> pattern, out long lo, out long hi)
    {
        lo = 0;
        hi = -1;
        if (pattern.IsEmpty) throw new ArgumentException("the empty pattern has no rank interval; handle it before calling", nameof(pattern));
        if (_length == 0) return false;
        if (IsMapped) return MatchMapped(pattern, out lo, out hi);
        bool found = TryMatch(pattern, out int l, out int h);
        lo = l;
        hi = h;
        return found;
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

    /// <summary>
    /// <see cref="Match{TSym}"/> for the packed text of a container: the pattern is mapped to dense symbols (a code unit that is not in the alphabet is no match) and the descent of <see cref="MappedSearch"/> runs over the mapped sections.
    /// </summary>
    private bool MatchMapped(ReadOnlySpan<char> pattern, out int lo, out int hi)
    {
        lo = 0;
        hi = -1;
        int[] map = SymbolOf;
        Span<byte> pat = pattern.Length <= StackPatternLimit ? stackalloc byte[pattern.Length] : new byte[pattern.Length];
        for (int i = 0; i < pattern.Length; i++)
        {
            char c = pattern[i];
            int s = c < map.Length ? map[c] : -1;
            if (s < 0) return false;
            pat[i] = (byte)s;
        }
        return _mappedSearch!.TryMatch(pat, out lo, out hi);
    }

    /// <summary><see cref="MatchMapped(ReadOnlySpan{char}, out int, out int)"/> at 64-bit ranks.</summary>
    private bool MatchMapped(ReadOnlySpan<char> pattern, out long lo, out long hi)
    {
        lo = 0;
        hi = -1;
        int[] map = SymbolOf;
        Span<byte> pat = pattern.Length <= StackPatternLimit ? stackalloc byte[pattern.Length] : new byte[pattern.Length];
        if (!Dense(map, pattern, pat)) return false;
        return _mappedSearch!.TryMatch(pat, out lo, out hi);
    }

    /// <summary>The pattern as dense symbols; false when a code unit is not in the alphabet of the text.</summary>
    private static bool Dense(int[] map, ReadOnlySpan<char> pattern, Span<byte> dense)
    {
        for (int i = 0; i < pattern.Length; i++)
        {
            char c = pattern[i];
            int s = c < map.Length ? map[c] : -1;
            if (s < 0) return false;
            dense[i] = (byte)s;
        }
        return true;
    }

    /// <summary>
    /// The readable reference of <see cref="MappedSearch"/> (work package 1b.9a replaced it on the hot path): <see cref="Descend{TChild, TSym}"/> over the sections of a container through the typed views (<see cref="SaView"/>, <see cref="CompactLcp"/>,
    /// <see cref="MappedChild"/>, <see cref="MappedTextReader"/>). The tests compare the two on valid and on damaged files. The same walk, with the reads of an untrusted file: the suffix array is a <see cref="SaView"/>, the text is the packed reader (a symbol at or beyond the
    /// end reads as 0, so every comparison is clamped to the suffix as the in-memory one is), and in checked mode the walk is bounded: the l-indices of an interval must increase and stay inside it, and a child's value must be
    /// larger than its parent's, so that every iteration makes progress and a damaged table ends in an <see cref="InvalidDataException"/> instead of a loop.
    /// </summary>
    internal static bool DescendMapped(MappedChild child, SaView sa, CompactLcp lcp, MappedTextReader text, ReadOnlySpan<byte> pat, bool isChecked, out int lo, out int hi)
    {
        int n = IntBridge.Length(text.Length, "EsaIndex.DescendMapped (the 32-bit reference; use DescendMappedWide)");
        int m = pat.Length;
        int lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        Span<byte> edge = stackalloc byte[EdgeChunk];
        while (c < m)
        {
            int symbol = pat[c];
            int a = lb;
            int boundary = child.FirstLIndex(lb, rb, isLast);
            if (isChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamagedChildTable();
            int b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                int p = sa[a] + ell;
                int key = (uint)p < (uint)n ? (int)text.Symbol(p) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0)
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                a = boundary;
                boundary = child.NextLIndex(boundary, rb);
                if (isChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamagedChildTable();
            }

            bool childIsLast = boundary < 0;
            if (a == b)
            {
                // a leaf: the rest of the pattern must be a prefix of the rest of the suffix
                int start = sa[a] + c;
                if ((uint)start > (uint)n || m - c > n - start || !EdgeEquals(text, start, pat.Slice(c), edge))
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
            if (isChecked && (childEll <= ell || childEll > n)) ThrowDamagedChildTable();
            int limit = Math.Min(childEll, m);
            int from = sa[a] + c;
            if (limit > c + 1 && ((uint)from >= (uint)n || !EdgeEquals(text, from + 1, pat.Slice(c + 1, limit - c - 1), edge)))
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

    /// <summary>
    /// The readable reference of the long descent (<see cref="WideSearch"/>): <see cref="DescendMapped"/> at 64-bit ranks, positions and lengths, over the long members of the typed views (<see cref="SaView"/>, <see cref="CompactLcp"/>,
    /// <see cref="MappedChild"/>, <see cref="MappedTextReader"/>), with the same checks in Checked mode. It reads any stored width and any LCP kind; the tests compare it with the long descent and with the 32-bit descent on valid and damaged files.
    /// </summary>
    internal static bool DescendMappedWide(MappedChild child, SaView sa, CompactLcp lcp, MappedTextReader text, ReadOnlySpan<byte> pat, bool isChecked, out long lo, out long hi)
    {
        long n = text.Length;
        long m = pat.Length;
        long lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        Span<byte> edge = stackalloc byte[EdgeChunk];
        while (c < m)
        {
            int symbol = pat[checked((int)c)]; // c < m = pat.Length
            long a = lb;
            long boundary = child.FirstLIndex(lb, rb, isLast);
            if (isChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamagedChildTable();
            long b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                long p = sa[a] + ell;
                int key = (ulong)p < (ulong)n ? (int)text.Symbol(p) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0)
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                a = boundary;
                boundary = child.NextLIndex(boundary, rb);
                if (isChecked && boundary >= 0 && (boundary <= a || boundary > rb)) ThrowDamagedChildTable();
            }

            bool childIsLast = boundary < 0;
            if (a == b)
            {
                long start = sa[a] + c;
                if ((ulong)start > (ulong)n || m - c > n - start || !EdgeEquals(text, start, pat[checked((int)c)..], edge))
                {
                    lo = 0;
                    hi = -1;
                    return false;
                }
                lo = a;
                hi = a;
                return true;
            }

            long childEll = lcp[child.FirstLIndex(a, b, childIsLast)];
            if (isChecked && (childEll <= ell || childEll > n)) ThrowDamagedChildTable();
            long limit = Math.Min(childEll, m);
            long from = sa[a] + c;
            if (limit > c + 1 && ((ulong)from >= (ulong)n || !EdgeEquals(text, from + 1, pat.Slice(checked((int)(c + 1)), checked((int)(limit - c - 1))), edge)))
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

    private const int EdgeChunk = 256;

    /// <summary>Whether the symbols of the text from <paramref name="start"/> equal <paramref name="expected"/> (the caller checked that they lie inside the text), compared in chunks through the packed reader.</summary>
    private static bool EdgeEquals(MappedTextReader text, long start, ReadOnlySpan<byte> expected, Span<byte> buffer)
    {
        if ((ulong)start > (ulong)text.Length || expected.Length > text.Length - start) return false;
        for (int done = 0; done < expected.Length; done += buffer.Length)
        {
            int take = Math.Min(buffer.Length, expected.Length - done);
            text.ReadSymbols(start + done, buffer[..take]);
            if (!buffer[..take].SequenceEqual(expected.Slice(done, take))) return false;
        }
        return true;
    }

    private static void ThrowDamagedChildTable() => throw new InvalidDataException("the child table or the LCP array is inconsistent (the l-indices of an interval must increase and stay inside it, and a child must be deeper than its parent): the file is damaged");

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
