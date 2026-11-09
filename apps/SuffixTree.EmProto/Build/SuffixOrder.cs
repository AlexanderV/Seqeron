using System.Runtime.CompilerServices;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Build;

/// <summary>Counters of the work done by one task / thread (ADR-004 E1: share of key ties, total LCE work). Summed after the build.</summary>
internal sealed class WorkCounters
{
    /// <summary>Calls of the suffix comparator.</summary>
    internal long Compares;

    /// <summary>Comparator calls that came from a key tie (the first 64 key bits equal).</summary>
    internal long TiedCompares;

    /// <summary>Longest-common-extension calls (each comparator call has at most one).</summary>
    internal long LceCalls;

    /// <summary>Symbols matched by the LCE calls, beyond the part already known equal from the key: the "LCE work" of E1.</summary>
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

    internal void Add(WorkCounters o)
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
/// Division and remainder of a 32-bit value by a fixed 32-bit divisor with two multiplications (Lemire, Kaser, Kurz, "Faster remainder by direct computation", 2019): the comparator
/// of the difference-cover sample does five of them per call and a 64-bit hardware division costs 25-90 cycles on the machines of this class.
/// </summary>
internal readonly struct FastMod32
{
    private readonly ulong _m;
    private readonly uint _d;

    internal FastMod32(uint divisor)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(divisor, 2u);
        _d = divisor;
        _m = (ulong.MaxValue / divisor) + 1;
    }

    internal uint Divisor => _d;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Div(uint a) => (uint)Math.BigMul(_m, a, out _);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Mod(uint a) => (uint)Math.BigMul(_m * a, _d, out _);
}

/// <summary>
/// The difference-cover sample ranks: for every position p of <c>[0, n]</c> with <c>p mod X</c> in D the rank of the suffix p among the sample suffixes (a dense
/// <c>uint</c> per slot, slot = <c>(p / X) * |D| + index of (p mod X) in D</c>), in native memory or in a mapped file.
/// </summary>
internal sealed unsafe class DcsRankTable(DifferenceCover cover, IMemoryBlock ranks, long n) : IDisposable
{
    private readonly FastMod32 _mod = new((uint)cover.Period);
    private readonly int[] _indexOf = cover.IndexOf;
    private readonly int _m = cover.Residues.Length;

    internal DifferenceCover Cover { get; } = cover;

    /// <summary>Division and remainder by the period without a hardware division.</summary>
    internal FastMod32 Mod => _mod;

    internal long N { get; } = n;

    internal uint* Ranks { get; } = (uint*)ranks.Pointer;

    internal IMemoryBlock Owner { get; } = ranks;

    /// <summary>Number of slots (including slots of the last block that lie beyond n).</summary>
    internal long Slots => ((N / Cover.Period) + 1) * Cover.Residues.Length;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Rank(long p)
    {
        uint block = _mod.Div((uint)p);
        uint r = (uint)p - (block * _mod.Divisor);
        return Ranks[((long)block * _m) + _indexOf[r]];
    }

    public void Dispose() => Owner.Dispose();
}

/// <summary>
/// The suffix comparator of the bucket partition and the bucket sort. Without a rank table it is the direct comparison (an LCE over the packed text, E1);
/// with a table it is Karkkainen's constant-time comparison: compare <c>l &lt; X</c> symbols directly, where <c>l</c> is chosen so that both <c>i + l</c> and <c>j + l</c> are sample positions
/// (the cover property), and then compare the two sample ranks (E2). One instance per thread (it owns the <see cref="WorkCounters"/>); text and table are shared.
/// </summary>
internal sealed class SuffixOrder(PackedText text, DcsRankTable? dcs, WorkCounters counters)
{
    internal PackedText Text { get; } = text;

    internal DcsRankTable? Dcs { get; } = dcs;

    internal WorkCounters Counters { get; } = counters;

    /// <summary>Compares the suffixes at <paramref name="i"/> and <paramref name="j"/> from their first symbol.</summary>
    internal int Compare(uint i, uint j) => CompareCore(i, j, tied: false);

    /// <summary>Compares two suffixes whose 64-bit keys are known to be equal (the first <c>KeySymbols</c> symbols are not looked at again when both suffixes are that long).</summary>
    internal int CompareTied(uint i, uint j) => CompareCore(i, j, tied: true);

    private int CompareCore(long i, long j, bool tied)
    {
        WorkCounters c = Counters;
        c.Compares++;
        if (tied) c.TiedCompares++;
        if (i == j) return 0;
        PackedText t = Text;
        long maxLen = t.Length - Math.Max(i, j);
        long skip = tied && t.KeySymbols <= maxLen ? t.KeySymbols : 0;
        DcsRankTable? d = Dcs;
        if (d is null)
        {
            c.LceCalls++;
            int r = t.CompareSuffixes(i, j, skip, ref c.LceSteps, out long m);
            c.LceSymbols += m - Math.Min(skip, maxLen);
            return r;
        }

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

/// <summary>Struct comparer over positions for <see cref="Span{T}"/> sorting: key ties of a leaf bucket.</summary>
internal readonly struct TiedComparer(SuffixOrder order) : IComparer<uint>
{
    public int Compare(uint x, uint y) => order.CompareTied(x, y);
}

/// <summary>Struct comparer over positions: full comparison (splitter sample).</summary>
internal readonly struct FullComparer(SuffixOrder order) : IComparer<uint>
{
    public int Compare(uint x, uint y) => order.Compare(x, y);
}
