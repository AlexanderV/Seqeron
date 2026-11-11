namespace SuffixTree.EmProto.Dcs;

/// <summary>
/// The arithmetic of the difference-cover sample of one text (ADR-004 PH1): which positions belong to the sample, where a position lives in the rank file (the <em>slot</em>, block-major:
/// <c>(p / X) * |D| + index of (p mod X) in D</c>, the layout of <see cref="DcsRankTable"/>), and where it lives in the reduced string (class-major: the positions of residue class c in increasing order,
/// then one separator <c>0</c>, class after class). The real slots are exactly <c>0 .. Samples - 1</c> (the unused slots of the last block are the last ones); the reduced string has
/// <c>Samples + |D|</c> symbols, the k-th member of class c at index <c>RedStart[c] + k</c> and the separator of class c at <c>RedStart[c] + ClassCount[c] = RedStart[c + 1] - 1</c>.
/// </summary>
internal sealed class SampleGeometry
{
    internal SampleGeometry(long n, DifferenceCover cover)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        N = n;
        Cover = cover;
        Period = cover.Period;
        Classes = cover.Residues.Length;
        Samples = cover.SampleCount(n);
        Slots = ((n / Period) + 1) * Classes;
        ClassCount = new long[Classes];
        ClassStart = new long[Classes + 1];
        RedStart = new long[Classes + 1];
        for (int c = 0; c < Classes; c++)
        {
            int r = cover.Residues[c];
            ClassCount[c] = r <= n ? ((n - r) / Period) + 1 : 0;
            ClassStart[c + 1] = ClassStart[c] + ClassCount[c];
            RedStart[c + 1] = ClassStart[c + 1] + c + 1;
        }
    }

    /// <summary>The text length n (the suffix at n is the empty suffix and is part of the sample when <c>n mod X</c> is in D).</summary>
    internal long N { get; }

    internal DifferenceCover Cover { get; }

    /// <summary>X.</summary>
    internal int Period { get; }

    /// <summary>|D|.</summary>
    internal int Classes { get; }

    /// <summary>Number of sample positions in <c>[0, n]</c>.</summary>
    internal long Samples { get; }

    /// <summary>Slots of the rank file, including the unused ones of the last block.</summary>
    internal long Slots { get; }

    /// <summary>Members of residue class c.</summary>
    internal long[] ClassCount { get; }

    /// <summary>Number of sample positions in the classes before c (cumulative, <c>Classes + 1</c> entries).</summary>
    internal long[] ClassStart { get; }

    /// <summary>Index in the reduced string of the first member of class c (<c>ClassStart[c] + c</c>; the last entry is the length of the reduced string).</summary>
    internal long[] RedStart { get; }

    /// <summary>Symbols of the reduced string: one per sample position and one separator per class.</summary>
    internal long RedLength => RedStart[Classes];

    /// <summary>Slot of the sample position <paramref name="p"/> (the caller guarantees that <c>p mod X</c> is in D).</summary>
    internal long SlotOf(long p) => (p / Period * Classes) + Cover.IndexOf[(int)(p % Period)];

    /// <summary>Position of the sample record in slot <paramref name="slot"/> (a real slot, below <see cref="Samples"/>).</summary>
    internal long PositionOfSlot(long slot) => (slot / Classes * Period) + Cover.Residues[(int)(slot % Classes)];

    /// <summary>Index in the reduced string of the sample record in slot <paramref name="slot"/>.</summary>
    internal long RedIndexOfSlot(long slot) => RedStart[(int)(slot % Classes)] + (slot / Classes);

    /// <summary>
    /// The slot of the reduced-string index <paramref name="q"/> (<c>0 &lt;= q &lt; RedLength</c>); false when q is the separator of its class.
    /// </summary>
    internal bool TrySlotOfRed(long q, out long slot)
    {
        int lo = 0, hi = Classes - 1; // largest c with RedStart[c] <= q
        while (lo < hi)
        {
            int mid = (lo + hi + 1) >> 1;
            if (RedStart[mid] <= q) lo = mid;
            else hi = mid - 1;
        }
        long k = q - RedStart[lo];
        if (k >= ClassCount[lo])
        {
            slot = -1;
            return false;
        }
        slot = (k * Classes) + lo;
        return true;
    }
}
