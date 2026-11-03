using System.Numerics;

namespace SuffixTree.Esa;

/// <summary>
/// Reference suffix-array construction: SA-IS (Nong, Zhang, Chan: "Two efficient algorithms for linear time suffix array
/// construction", IEEE Trans. Computers 60(10), 2011), implemented from the published algorithm in safe code, generic over the
/// element type of the suffix array and of the buckets. Clarity and correctness first; speed is the business of the later
/// optimized kernels, which must equal this implementation bit for bit (it is the permanent oracle, ADR-002 section 5, R0).
/// <para>
/// Index width (PROGRESS-S1.md 1a.5, review MAJ-5): <c>int</c> is the production width; the <c>long</c> instantiation exists only as an
/// oracle for the index arithmetic and must give the identical array. Spans are <c>int</c>-indexed, so the reference never handles more than
/// <c>Array.MaxLength</c> symbols whatever <typeparamref name="TIdx"/> is. Every conversion is checked (an index type that is too narrow
/// throws <see cref="OverflowException"/> instead of wrapping). The optimized kernels are <c>int</c>-only; the 64-bit kernels of stage S2 are
/// generated from the <c>int</c> source and tested against <c>SaisReference&lt;long&gt;</c>.
/// </para>
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
/// <typeparam name="TIdx">Signed integer type of the suffix array and of the bucket arrays (<c>int</c> or <c>long</c>).</typeparam>
internal static class SaisReference<TIdx>
    where TIdx : unmanaged, IBinaryInteger<TIdx>, ISignedNumber<TIdx>
{
    private static TIdx Empty => TIdx.NegativeOne;

    private static int ToInt(TIdx value) => int.CreateChecked(value);

    private static TIdx FromInt(int value) => TIdx.CreateChecked(value);

    /// <summary>Builds the suffix array of a string (ordinal UTF-16 order) through the order-preserving dense remap.</summary>
    internal static TIdx[] BuildSuffixArray(ReadOnlySpan<char> text, SaisTrace? trace = null)
    {
        PreparedText prepared = TextPreparation.Prepare(text);
        return prepared.Width == SymbolWidth.Byte
            ? Build(prepared.ByteReader(), prepared.AlphabetSize, trace)
            : Build(prepared.UShortReader(), prepared.AlphabetSize, trace);
    }

    /// <summary>Builds the suffix array of a symbol sequence over <c>[0, alphabetSize)</c> into a new array of <c>n + 1</c> entries.</summary>
    /// <exception cref="NotSupportedException">The suffix array would not fit the 32-bit path (checked before any allocation).</exception>
    internal static TIdx[] Build<TReader>(TReader text, int alphabetSize, SaisTrace? trace = null)
        where TReader : ISymbolReader, allows ref struct
    {
        EsaLimits.EnsureReaderFits(text);
        ArgumentOutOfRangeException.ThrowIfNegative(alphabetSize);
        var sa = new TIdx[text.Length + 1];
        Sais(text, alphabetSize, sa, trace);
        return sa;
    }

    /// <summary>Builds the suffix array into <paramref name="sa"/>, which must have exactly <c>n + 1</c> entries.</summary>
    internal static void Build<TReader>(TReader text, int alphabetSize, Span<TIdx> sa, SaisTrace? trace = null)
        where TReader : ISymbolReader, allows ref struct
    {
        EsaLimits.EnsureReaderFits(text);
        ArgumentOutOfRangeException.ThrowIfNegative(alphabetSize);
        if (sa.Length != text.Length + 1)
            throw new ArgumentException($"the suffix array needs n + 1 = {text.Length + 1} entries, got {sa.Length}", nameof(sa));
        Sais(text, alphabetSize, sa, trace);
    }

    private static void Sais<TReader>(TReader t, int sigma, Span<TIdx> sa, SaisTrace? trace)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        sa[0] = FromInt(n);
        if (n == 0)
        {
            trace?.Levels.Add(new SaisLevel(0, sigma, 0));
            return;
        }

        TIdx[] counts = CountSymbols(t, sigma);
        bool[] isS = SaisReference.ClassifyTypes(t);
        int m = SaisReference.CountLms(isS);
        trace?.Levels.Add(new SaisLevel(n, sigma, m));
        var bkt = new TIdx[sigma];

        // Stage 1: sort the LMS substrings by one induction round, name them, solve the reduced problem.
        SortLmsSubstrings(t, sa, isS, counts, bkt);
        CompactSortedLms(sa, isS, m);
        int names = NameLmsSubstrings(t, sa, isS, m);
        CompactNames(sa, m);
        SolveReduced(sa, isS, m, names, trace);

        // Stage 2: place the LMS suffixes in their final order and induce all others from them.
        PlaceSortedLmsAndInduce(t, sa, isS, counts, bkt, m);
    }

    private static TIdx[] CountSymbols<TReader>(TReader t, int sigma)
        where TReader : ISymbolReader, allows ref struct
    {
        var counts = new TIdx[sigma];
        for (int i = 0; i < t.Length; i++)
        {
            int c = t[i];
            if ((uint)c >= (uint)sigma)
                throw new ArgumentException($"symbol {c} at index {i} is outside [0, {sigma})", nameof(t));
            counts[c]++;
        }
        return counts;
    }

    private static void SetBucketHeads(TIdx[] counts, TIdx[] bkt)
    {
        int sum = 1; // index 0 belongs to the sentinel suffix
        for (int c = 0; c < counts.Length; c++)
        {
            bkt[c] = FromInt(sum);
            sum += ToInt(counts[c]);
        }
    }

    private static void SetBucketTails(TIdx[] counts, TIdx[] bkt)
    {
        int sum = 1;
        for (int c = 0; c < counts.Length; c++)
        {
            sum += ToInt(counts[c]);
            bkt[c] = FromInt(sum); // exclusive end
        }
    }

    private static void InduceL<TReader>(TReader t, Span<TIdx> sa, bool[] isS, TIdx[] counts, TIdx[] bkt)
        where TReader : ISymbolReader, allows ref struct
    {
        SetBucketHeads(counts, bkt);
        for (int i = 0; i <= t.Length; i++)
        {
            int p = ToInt(sa[i]);
            if (p <= 0) continue; // empty slot, or suffix 0 which has no predecessor
            int j = p - 1;
            if (isS[j]) continue;
            int c = t[j];
            int slot = ToInt(bkt[c]);
            sa[slot] = FromInt(j);
            bkt[c] = FromInt(slot + 1);
        }
    }

    private static void InduceS<TReader>(TReader t, Span<TIdx> sa, bool[] isS, TIdx[] counts, TIdx[] bkt)
        where TReader : ISymbolReader, allows ref struct
    {
        SetBucketTails(counts, bkt);
        for (int i = t.Length; i >= 1; i--)
        {
            int p = ToInt(sa[i]);
            if (p <= 0) continue;
            int j = p - 1;
            if (!isS[j]) continue;
            int c = t[j];
            int slot = ToInt(bkt[c]) - 1;
            sa[slot] = FromInt(j);
            bkt[c] = FromInt(slot);
        }
    }

    /// <summary>Stage 1: LMS positions into the ends of their buckets (any order), then induce L and S; afterwards the LMS substrings are in sorted order.</summary>
    private static void SortLmsSubstrings<TReader>(TReader t, Span<TIdx> sa, bool[] isS, TIdx[] counts, TIdx[] bkt)
        where TReader : ISymbolReader, allows ref struct
    {
        int n = t.Length;
        sa.Slice(1).Fill(Empty);
        SetBucketTails(counts, bkt);
        for (int i = n - 1; i >= 1; i--)
        {
            if (!SaisReference.IsLms(isS, i)) continue;
            int c = t[i];
            int slot = ToInt(bkt[c]) - 1;
            sa[slot] = FromInt(i);
            bkt[c] = FromInt(slot);
        }
        InduceL(t, sa, isS, counts, bkt);
        InduceS(t, sa, isS, counts, bkt);
    }

    /// <summary>Moves the sorted LMS positions to <c>SA[1..m]</c> (stable, in place: the write index never passes the read index).</summary>
    private static void CompactSortedLms(Span<TIdx> sa, bool[] isS, int m)
    {
        int k = 0;
        for (int i = 1; i < sa.Length; i++)
        {
            if (SaisReference.IsLms(isS, ToInt(sa[i]))) sa[1 + k++] = sa[i];
        }
        if (k != m)
            throw new InvalidOperationException($"internal error: {k} LMS positions after induction, expected {m}");
    }

    /// <summary>Names the sorted LMS substrings (equal substrings share a name) and stores the name of the LMS position p at <c>SA[1 + m + (p &gt;&gt; 1)]</c>. Returns the number of distinct names.</summary>
    private static int NameLmsSubstrings<TReader>(TReader t, Span<TIdx> sa, bool[] isS, int m)
        where TReader : ISymbolReader, allows ref struct
    {
        sa.Slice(m + 1).Fill(Empty);
        int name = -1;
        int previous = -1;
        for (int i = 0; i < m; i++)
        {
            int pos = ToInt(sa[1 + i]);
            if (previous < 0 || !SaisReference.LmsSubstringsEqual(t, isS, previous, pos)) name++;
            sa[1 + m + (pos >> 1)] = FromInt(name);
            previous = pos;
        }
        return name + 1;
    }

    /// <summary>Gathers the names into the tail <c>SA[n + 1 - m .. n]</c> in text order (the reduced string s1).</summary>
    private static void CompactNames(Span<TIdx> sa, int m)
    {
        int n = sa.Length - 1;
        int j = n;
        for (int i = n; i > m; i--)
            if (!TIdx.IsNegative(sa[i])) sa[j--] = sa[i];
    }

    /// <summary>
    /// Solves the reduced problem and replaces <c>SA[1..m]</c> by the LMS positions in sorted suffix order. The reduced text sits in
    /// the tail, its suffix array (n + 1 layout, so m + 1 entries) in the head; <c>SA[0]</c> is restored to n afterwards.
    /// </summary>
    private static void SolveReduced(Span<TIdx> sa, bool[] isS, int m, int names, SaisTrace? trace)
    {
        int n = sa.Length - 1;
        Span<TIdx> sa1 = sa.Slice(0, m + 1);
        ReadOnlySpan<TIdx> s1 = sa.Slice(n + 1 - m, m);
        if (names < m)
        {
            Sais(new IndexSymbols<TIdx>(s1), names, sa1, trace);
        }
        else
        {
            sa1[0] = FromInt(m);
            for (int i = 0; i < m; i++) sa1[1 + ToInt(s1[i])] = FromInt(i); // all names distinct: the name is the rank
        }

        // reuse the tail for the LMS positions in text order, then translate ranks of s1 into text positions
        int tailBase = n + 1 - m;
        int idx = m;
        for (int i = n - 1; i >= 1; i--)
            if (SaisReference.IsLms(isS, i)) sa[tailBase + --idx] = FromInt(i);
        for (int i = 1; i <= m; i++)
            sa[i] = sa[tailBase + ToInt(sa[i])];
        sa[0] = FromInt(n);
    }

    /// <summary>Stage 2: the sorted LMS suffixes go to the ends of their buckets, largest first, then L and S are induced.</summary>
    private static void PlaceSortedLmsAndInduce<TReader>(TReader t, Span<TIdx> sa, bool[] isS, TIdx[] counts, TIdx[] bkt, int m)
        where TReader : ISymbolReader, allows ref struct
    {
        sa.Slice(m + 1).Fill(Empty);
        SetBucketTails(counts, bkt);
        for (int r = m; r >= 1; r--)
        {
            TIdx p = sa[r];
            sa[r] = Empty;
            int c = t[ToInt(p)];
            int slot = ToInt(bkt[c]) - 1;
            sa[slot] = p; // the target index is >= r, so no unread LMS entry is overwritten
            bkt[c] = FromInt(slot);
        }
        InduceL(t, sa, isS, counts, bkt);
        InduceS(t, sa, isS, counts, bkt);
    }
}
