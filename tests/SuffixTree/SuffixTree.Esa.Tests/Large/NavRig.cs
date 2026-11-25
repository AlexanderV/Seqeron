using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>The text of a closed form as the letters that <see cref="EsaNavigator64.TryChild"/> is asked for.</summary>
internal sealed class GenLetters(ClosedFormText form) : ILetters
{
    /// <inheritdoc />
    public int At(long position) => form.CharAt(position);
}

/// <summary>The cells of a form without a child table: any read is a failure of the test.</summary>
internal sealed class NoCells : IWideCells
{
    /// <inheritdoc />
    public long First(long lb, long rb, bool lastChild) => throw new NotSupportedException("this closed form has no child table");

    /// <inheritdoc />
    public long Next(long k, long rb) => throw new NotSupportedException("this closed form has no child table");

    /// <inheritdoc />
    public long Depth(long first) => throw new NotSupportedException("this closed form has no child table");
}

/// <summary>
/// The long navigation engine over the generated sections of a closed-form text at any n (work package S2.5, ADR-005 section 3.7 layer 1): the suffix array, PSI and the RMQ levels as <see cref="PosGen{TWidth}"/> readers of 5 bytes,
/// the LCP array as U64, the packed child cells; no memory at any n. The oracles are the closed form itself: <see cref="ClosedFormText.RangeOfFactor"/> (the interval of a factor), <see cref="ClosedFormText.LcpMin"/>,
/// <see cref="ClosedFormText.Lcp"/>, <see cref="ClosedFormText.SuffixAt"/>. Every read of a generated section is checked.
/// </summary>
internal sealed class NavRig
{
    internal NavRig(ClosedFormText f)
    {
        Form = f;
        long n = f.Length;
        Lcp = GenDescent.Lcp64(f);
        var sa = PosSource.Of(new PosGen<W5>(f, GenSection.Sa));
        var psi = PosSource.Of(new PosGen<W5>(f, GenSection.Psi));
        Ranges = WideLcpPyramid.FromLevels(Lcp, SesaFormat.RmqFanout, PosSource.Of(new PosGen<W5>(f, GenSection.Rmq)));
        Pack = EsaNavigationPack.FromSources(n, sa, psi, Ranges);
        IWideCells cells = f.HasChildTable
            ? new WideCells<CheckedAccess, PosGen<W5>, PosGen<W8>>(new PosGen<W5>(f, GenSection.ChildPacked), new PosGen<W8>(f, GenSection.LcpBase), Lcp, tagged: false, n)
            : new NoCells();
        Nav = new EsaNavigator64(n, sa, cells, new GenLetters(f), () => Pack, mapped: true);
    }

    internal ClosedFormText Form { get; }

    internal CompactLcp Lcp { get; }

    internal WideLcpPyramid Ranges { get; }

    internal EsaNavigationPack Pack { get; }

    internal EsaNavigator64 Nav { get; }

    /// <summary>The node of the rank interval <paramref name="lo"/>..<paramref name="hi"/> at string depth <paramref name="depth"/>, with its parent depth, last-child flag and position as the definitions give them.</summary>
    internal EsaNode64 NodeOf(long lo, long hi, long depth)
    {
        long n = Form.Length;
        long left = Form.Lcp(lo), right = hi < n ? Form.Lcp(hi + 1) : -1;
        return new EsaNode64(lo, hi, depth, Math.Max(left, right), Form.SuffixAt(lo), left > right);
    }

    /// <summary>True when the interval is an lcp-interval of value <paramref name="depth"/> (an internal node of that depth) and not a leaf.</summary>
    internal bool IsInternalNode(long lo, long hi, long depth)
    {
        if (lo >= hi || depth <= 0) return false;
        long n = Form.Length;
        long right = hi < n ? Form.Lcp(hi + 1) : -1;
        return Form.LcpMin(lo + 1, hi) == depth && Form.Lcp(lo) < depth && right < depth;
    }

    /// <summary>The internal node of the factor <c>text[start, start + length)</c>, or null when its interval is not a node (the factor ends inside an edge, or occurs once).</summary>
    internal EsaNode64? NodeOfFactor(long start, long length)
    {
        RankRange r = Form.RangeOfFactor(start, length);
        return !r.IsEmpty && IsInternalNode(r.First, r.Last, length) ? NodeOf(r.First, r.Last, length) : null;
    }

    /// <summary>The factors <c>(start, length)</c> that a test looks at in a text of n symbols: short ones at the start and the end, and the lengths at which a periodic text branches (multiples of the period and their neighbours).</summary>
    internal IEnumerable<(long Start, long Length)> Factors(int startsLimit = 40)
    {
        long n = Form.Length, p = Form.Period;
        var starts = new SortedSet<long>();
        for (long s = 0; s < Math.Min(startsLimit, n); s++) starts.Add(s);
        for (long s = Math.Max(0, n - startsLimit); s < n; s++) starts.Add(s);
        foreach (long s in new[] { n / 2, n / 3, n / 2 + 1, (1L << 31) - 1, 1L << 31, (1L << 31) + 1 }.Where(s => s >= 0 && s < n)) starts.Add(s);
        var lengths = new SortedSet<long>();
        for (long l = 1; l <= 14; l++) lengths.Add(l);
        long q = n / p;
        foreach (long j in new long[] { 1, 2, 3, 4, q / 3, q / 2, q - 3, q - 2, q - 1, q })
            foreach (long d in new long[] { -p - 1, -p, -2, -1, 0, 1, 2, p - 1, p })
                if (j * p + d >= 1) lengths.Add((j * p) + d);
        lengths.Add(n - 1);
        lengths.Add(n - 2);
        foreach (long s in starts)
        {
            foreach (long l in lengths)
            {
                if (l >= 1 && s + l <= n) yield return (s, l);
            }
        }
    }
}
