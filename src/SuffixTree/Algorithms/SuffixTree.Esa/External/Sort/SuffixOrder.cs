using System.Runtime.CompilerServices;

namespace SuffixTree.Esa.External;

/// <summary>Counters of the work done by one task / thread (ADR-004 E1: share of key ties, total LCE work). Summed after the build.</summary>
internal sealed class SortWork
{
    /// <summary>Calls of the suffix comparator.</summary>
    internal long Compares;

    /// <summary>Comparator calls that came from a key tie (the first 64 key bits equal).</summary>
    internal long TiedCompares;

    /// <summary>Longest-common-extension calls (each comparator call has at most one).</summary>
    internal long LceCalls;

    /// <summary>Symbols matched by the LCE calls, beyond the part already known equal from the key.</summary>
    internal long LceSymbols;

    /// <summary>Loop steps of the LCE (vector steps of 256 bits, or scalar words): the cost proxy.</summary>
    internal long LceSteps;

    /// <summary>Sample-rank lookups of the difference-cover comparator.</summary>
    internal long RankLookups;

    /// <summary>Positions of leaf buckets that share their key with another position of the same bucket.</summary>
    internal long TiePositions;

    /// <summary>Positions sorted in leaf buckets.</summary>
    internal long LeafPositions;

    /// <summary>Positions that went through a partition scan (summed over the levels).</summary>
    internal long ScanPositions;

    internal void Add(SortWork o)
    {
        Compares += o.Compares;
        TiedCompares += o.TiedCompares;
        LceCalls += o.LceCalls;
        LceSymbols += o.LceSymbols;
        LceSteps += o.LceSteps;
        RankLookups += o.RankLookups;
        TiePositions += o.TiePositions;
        LeafPositions += o.LeafPositions;
        ScanPositions += o.ScanPositions;
    }
}

/// <summary>
/// The suffix comparator of the bucket partition and the bucket sort: Karkkainen's constant-time comparison with the difference-cover sample ranks. Compare <c>l &lt; X</c> symbols directly, where <c>l</c> is chosen so that both
/// <c>i + l</c> and <c>j + l</c> are sample positions (the cover property), and then compare the two sample ranks. One instance per thread (it owns the <see cref="SortWork"/>); the text and the table are shared.
/// </summary>
internal sealed class SuffixOrder(MappedTextReader text, RankTable ranks, SortWork work)
{
    internal MappedTextReader Text { get; } = text;

    internal RankTable Ranks { get; } = ranks;

    internal SortWork Work { get; } = work;

    /// <summary>Compares the suffixes at <paramref name="i"/> and <paramref name="j"/> from their first symbol (a position at <c>n</c> is the empty suffix).</summary>
    internal int Compare(uint i, uint j) => CompareCore(i, j, tied: false);

    /// <summary>Compares two suffixes whose 64-bit keys are known to be equal (the first <c>KeySymbols</c> symbols are not looked at again when both suffixes are that long).</summary>
    internal int CompareTied(uint i, uint j) => CompareCore(i, j, tied: true);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private int CompareCore(long i, long j, bool tied)
    {
        SortWork c = Work;
        c.Compares++;
        if (tied) c.TiedCompares++;
        if (i == j) return 0;
        MappedTextReader t = Text;
        long maxLen = t.Length - Math.Max(i, j);
        long skip = tied && t.KeySymbols <= maxLen ? t.KeySymbols : 0;
        RankTable d = Ranks;
        FastMod32 mod = d.Mod;
        int x = (int)mod.Divisor;
        uint ri = mod.Mod((uint)i);
        uint rj = mod.Mod((uint)j);
        int delta = (int)rj - (int)ri;
        if (delta < 0) delta += x;
        int a = d.Cover.Witness[delta];
        int l = a - (int)ri;
        if (l < 0) l += x;
        long lim = Math.Min(l, maxLen);
        long s = Math.Min(skip, lim);
        long matched = s;
        if (lim > s)
        {
            c.LceCalls++;
            matched = s + t.Lce(i + s, j + s, lim - s, ref c.LceSteps);
            c.LceSymbols += matched - s;
        }
        if (matched < lim) return t.Symbol(i + matched) < t.Symbol(j + matched) ? -1 : 1;
        if (l > maxLen) return i > j ? -1 : 1; // one suffix has fewer than l symbols and ended inside the equal part
        c.RankLookups += 2;
        return d.Rank(i + l) < d.Rank(j + l) ? -1 : 1;
    }
}

/// <summary>Struct comparer over positions: full comparison (splitter sample).</summary>
internal readonly struct FullComparer(SuffixOrder order) : IComparer<uint>
{
    public int Compare(uint x, uint y) => order.Compare(x, y);
}
