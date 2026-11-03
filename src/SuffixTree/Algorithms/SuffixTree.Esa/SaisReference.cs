namespace SuffixTree.Esa;

/// <summary>
/// Reference suffix-array construction: SA-IS (Nong, Zhang, Chan: "Two efficient algorithms for linear time suffix array
/// construction", IEEE Trans. Computers 60(10), 2011), implemented from the published algorithm in safe code.
/// Clarity and correctness first; speed is the business of the later optimized kernels, which must equal this implementation
/// bit for bit (it is the permanent oracle, ADR-002 section 5, R0).
/// <para>
/// Layout (ADR-002 section 2 item 1): the result has <c>n + 1</c> entries, <c>SA[0] = n</c> is the empty suffix, i.e. a virtual
/// sentinel smaller than every symbol that is never stored in the text, and <c>SA[1..n]</c> are the non-empty suffixes in
/// ascending order. Buckets of the real symbols therefore start at index 1.
/// </para>
/// <para>
/// Memory plan of the recursion (ADR-002 section 2 item 5): the LMS substrings are sorted by one induction round, named
/// (<c>name</c> stored at <c>SA[1 + m + (pos &gt;&gt; 1)]</c>, unique because two LMS positions differ by at least 2) and compacted
/// into the reduced string <c>s1</c> in the tail of <c>SA</c>; the suffix array of <c>s1</c> is built in the head
/// <c>SA[0..m]</c> (head and tail are disjoint because <c>2m &lt;= n</c>), or placed directly when all names are distinct.
/// </para>
/// </summary>
internal static class SaisReference
{
    private const int Empty = -1;

    /// <summary>Builds the suffix array of a string (ordinal UTF-16 order) through the order-preserving dense remap.</summary>
    internal static int[] BuildSuffixArray(ReadOnlySpan<char> text, SaisTrace? trace = null)
    {
        PreparedText prepared = TextPreparation.Prepare(text);
        return prepared.Width == SymbolWidth.Byte
            ? Build(prepared.ByteReader(), prepared.AlphabetSize, trace)
            : Build(prepared.UShortReader(), prepared.AlphabetSize, trace);
    }

    /// <summary>Builds the suffix array of a symbol sequence over <c>[0, alphabetSize)</c> into a new array of <c>n + 1</c> entries.</summary>
    /// <exception cref="NotSupportedException">The suffix array would not fit the 32-bit path (checked before any allocation).</exception>
    internal static int[] Build<TReader>(TReader text, int alphabetSize, SaisTrace? trace = null)
        where TReader : ISymbolReader, allows ref struct
    {
        EsaLimits.EnsureReaderFits(text);
        ArgumentOutOfRangeException.ThrowIfNegative(alphabetSize);
        var sa = new int[text.Length + 1];
        Sais(text, alphabetSize, sa, trace);
        return sa;
    }

    /// <summary>Builds the suffix array into <paramref name="sa"/>, which must have exactly <c>n + 1</c> entries.</summary>
    internal static void Build<TReader>(TReader text, int alphabetSize, Span<int> sa, SaisTrace? trace = null)
        where TReader : ISymbolReader, allows ref struct
    {
        EsaLimits.EnsureReaderFits(text);
        ArgumentOutOfRangeException.ThrowIfNegative(alphabetSize);
        if (sa.Length != text.Length + 1)
            throw new ArgumentException($"the suffix array needs n + 1 = {text.Length + 1} entries, got {sa.Length}", nameof(sa));
        Sais(text, alphabetSize, sa, trace);
    }

    private static void Sais<TReader>(TReader t, int sigma, Span<int> sa, SaisTrace? trace)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        sa[0] = n;
        if (n == 0)
        {
            trace?.Levels.Add(new SaisLevel(0, sigma, 0));
            return;
        }

        int[] counts = CountSymbols(t, sigma);
        bool[] isS = ClassifyTypes(t);
        int m = CountLms(isS);
        trace?.Levels.Add(new SaisLevel(n, sigma, m));
        int[] bkt = new int[sigma];

        // Stage 1: sort the LMS substrings by one induction round, name them, solve the reduced problem.
        SortLmsSubstrings(t, sa, isS, counts, bkt);
        CompactSortedLms(sa, isS, m);
        int names = NameLmsSubstrings(t, sa, isS, m);
        CompactNames(sa, m);
        SolveReduced(sa, isS, m, names, trace);

        // Stage 2: place the LMS suffixes in their final order and induce all others from them.
        PlaceSortedLmsAndInduce(t, sa, isS, counts, bkt, m);
    }

    private static int[] CountSymbols<TReader>(TReader t, int sigma)
        where TReader : ISymbolReader, allows ref struct
    {
        var counts = new int[sigma];
        for (int i = 0; i < t.Length; i++)
        {
            int c = t[i];
            if ((uint)c >= (uint)sigma)
                throw new ArgumentException($"symbol {c} at index {i} is outside [0, {sigma})", nameof(t));
            counts[c]++;
        }
        return counts;
    }

    /// <summary>Type of every position, <c>true</c> = S-type. The virtual sentinel at index n is S-type, so position n - 1 is always L-type.</summary>
    internal static bool[] ClassifyTypes<TReader>(TReader t)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        var isS = new bool[n + 1];
        isS[n] = true;
        for (int i = n - 2; i >= 0; i--)
        {
            int a = t[i], b = t[i + 1];
            isS[i] = a < b || (a == b && isS[i + 1]);
        }
        return isS;
    }

    internal static bool IsLms(bool[] isS, int i) => i > 0 && isS[i] && !isS[i - 1];

    /// <summary>Number of LMS positions in <c>[1, n - 1]</c>; the sentinel position n is the (implicit) smallest LMS position and is not counted.</summary>
    private static int CountLms(bool[] isS)
    {
        int m = 0;
        for (int i = 1; i < isS.Length - 1; i++)
            if (IsLms(isS, i)) m++;
        return m;
    }

    private static void SetBucketHeads(int[] counts, int[] bkt)
    {
        int sum = 1; // index 0 belongs to the sentinel suffix
        for (int c = 0; c < counts.Length; c++)
        {
            bkt[c] = sum;
            sum += counts[c];
        }
    }

    private static void SetBucketTails(int[] counts, int[] bkt)
    {
        int sum = 1;
        for (int c = 0; c < counts.Length; c++)
        {
            sum += counts[c];
            bkt[c] = sum; // exclusive end
        }
    }

    private static void InduceL<TReader>(TReader t, Span<int> sa, bool[] isS, int[] counts, int[] bkt)
        where TReader : ISymbolReader, allows ref struct
    {
        SetBucketHeads(counts, bkt);
        for (int i = 0; i <= t.Length; i++)
        {
            int p = sa[i];
            if (p <= 0) continue; // empty slot, or suffix 0 which has no predecessor
            int j = p - 1;
            if (!isS[j]) sa[bkt[t[j]]++] = j;
        }
    }

    private static void InduceS<TReader>(TReader t, Span<int> sa, bool[] isS, int[] counts, int[] bkt)
        where TReader : ISymbolReader, allows ref struct
    {
        SetBucketTails(counts, bkt);
        for (int i = t.Length; i >= 1; i--)
        {
            int p = sa[i];
            if (p <= 0) continue;
            int j = p - 1;
            if (isS[j]) sa[--bkt[t[j]]] = j;
        }
    }

    /// <summary>Stage 1: LMS positions into the ends of their buckets (any order), then induce L and S; afterwards the LMS substrings are in sorted order.</summary>
    private static void SortLmsSubstrings<TReader>(TReader t, Span<int> sa, bool[] isS, int[] counts, int[] bkt)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        sa.Slice(1).Fill(Empty);
        SetBucketTails(counts, bkt);
        for (int i = n - 1; i >= 1; i--)
            if (IsLms(isS, i)) sa[--bkt[t[i]]] = i;
        InduceL(t, sa, isS, counts, bkt);
        InduceS(t, sa, isS, counts, bkt);
    }

    /// <summary>Moves the sorted LMS positions to <c>SA[1..m]</c> (stable, in place: the write index never passes the read index).</summary>
    private static void CompactSortedLms(Span<int> sa, bool[] isS, int m)
    {
        int k = 0;
        for (int i = 1; i < sa.Length; i++)
        {
            int p = sa[i];
            if (IsLms(isS, p)) sa[1 + k++] = p;
        }
        if (k != m)
            throw new InvalidOperationException($"internal error: {k} LMS positions after induction, expected {m}");
    }

    /// <summary>Names the sorted LMS substrings (equal substrings share a name) and stores the name of the LMS position p at <c>SA[1 + m + (p &gt;&gt; 1)]</c>. Returns the number of distinct names.</summary>
    private static int NameLmsSubstrings<TReader>(TReader t, Span<int> sa, bool[] isS, int m)
        where TReader : ISymbolReader, allows ref struct
    {
        sa.Slice(m + 1).Fill(Empty);
        int name = -1;
        int previous = -1;
        for (int i = 0; i < m; i++)
        {
            int pos = sa[1 + i];
            if (previous < 0 || !LmsSubstringsEqual(t, isS, previous, pos)) name++;
            sa[1 + m + (pos >> 1)] = name;
            previous = pos;
        }
        return name + 1;
    }

    /// <summary>
    /// Equality of the LMS substrings starting at the distinct LMS positions a and b (same length, same symbols, same types).
    /// The LMS substring of a is text[a..next LMS position] inclusive; the one that reaches the sentinel is unique.
    /// Internal so that the tests can check it against the definition on its own.
    /// </summary>
    internal static bool LmsSubstringsEqual<TReader>(TReader t, bool[] isS, int a, int b)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        for (int d = 0; a + d <= n && b + d <= n; d++)
        {
            int x = a + d, y = b + d;
            if (x == n || y == n) return false; // a != b, so at most one of them is the sentinel, which no other substring contains
            if (t[x] != t[y] || isS[x] != isS[y]) return false;
            if (d > 0 && IsLms(isS, x)) return true; // equal so far, so y is LMS as well: both substrings end here
        }
        return false; // not reached: the loop ends at the sentinel of the later substring
    }

    /// <summary>Gathers the names into the tail <c>SA[n + 1 - m .. n]</c> in text order (the reduced string s1).</summary>
    private static void CompactNames(Span<int> sa, int m)
    {
        int n = sa.Length - 1;
        int j = n;
        for (int i = n; i > m; i--)
            if (sa[i] >= 0) sa[j--] = sa[i];
    }

    /// <summary>
    /// Solves the reduced problem and replaces <c>SA[1..m]</c> by the LMS positions in sorted suffix order. The reduced text sits in
    /// the tail, its suffix array (n + 1 layout, so m + 1 entries) in the head; <c>SA[0]</c> is restored to n afterwards.
    /// </summary>
    private static void SolveReduced(Span<int> sa, bool[] isS, int m, int names, SaisTrace? trace)
    {
        int n = sa.Length - 1;
        Span<int> sa1 = sa.Slice(0, m + 1);
        ReadOnlySpan<int> s1 = sa.Slice(n + 1 - m, m);
        if (names < m)
        {
            Sais(new IntSymbols(s1), names, sa1, trace);
        }
        else
        {
            sa1[0] = m;
            for (int i = 0; i < m; i++) sa1[1 + s1[i]] = i; // all names distinct: the name is the rank
        }

        // reuse the tail for the LMS positions in text order, then translate ranks of s1 into text positions
        int tailBase = n + 1 - m;
        int idx = m;
        for (int i = n - 1; i >= 1; i--)
            if (IsLms(isS, i)) sa[tailBase + --idx] = i;
        for (int i = 1; i <= m; i++)
            sa[i] = sa[tailBase + sa[i]];
        sa[0] = n;
    }

    /// <summary>Stage 2: the sorted LMS suffixes go to the ends of their buckets, largest first, then L and S are induced.</summary>
    private static void PlaceSortedLmsAndInduce<TReader>(TReader t, Span<int> sa, bool[] isS, int[] counts, int[] bkt, int m)
        where TReader : ISymbolReader, allows ref struct
    {
        sa.Slice(m + 1).Fill(Empty);
        SetBucketTails(counts, bkt);
        for (int r = m; r >= 1; r--)
        {
            int p = sa[r];
            sa[r] = Empty;
            sa[--bkt[t[p]]] = p; // the target index is >= r, so no unread LMS entry is overwritten
        }
        InduceL(t, sa, isS, counts, bkt);
        InduceS(t, sa, isS, counts, bkt);
    }
}
