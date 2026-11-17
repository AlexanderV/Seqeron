namespace SuffixTree.EmProto.Verify;

/// <summary>An open lcp-interval of the scan: its value, its left boundary, its first and its last l-index (<see cref="IntervalStack.None"/> = none yet). 16 bytes.</summary>
internal struct OpenInterval
{
    internal uint Ell;
    internal uint Lb;
    internal uint First;
    internal uint Last;
}

/// <summary>
/// The stack of the scan of the child-table check: blocks of 2^20 entries allocated on demand, so that a stack of a few entries costs a few kilobytes and a stack of 15 million (an N run of 15 Mbp) 240 MB; a hard bound
/// of <see cref="Capacity"/> entries (a text whose LCP needs more is reported as skipped by the caller).
/// </summary>
internal sealed class IntervalStack(long capacity)
{
    internal const uint None = uint.MaxValue;
    private const int BlockBits = 20;
    private const int BlockSize = 1 << BlockBits;
    private readonly List<OpenInterval[]> _blocks = [];

    internal long Capacity { get; } = capacity;

    internal long Count { get; private set; }

    internal long Peak { get; private set; }

    /// <summary>False when the push would exceed the capacity (nothing is pushed).</summary>
    internal bool TryPush(in OpenInterval e)
    {
        if (Count >= Capacity) return false;
        int block = (int)(Count >> BlockBits);
        if (block == _blocks.Count) _blocks.Add(new OpenInterval[BlockSize]);
        _blocks[block][Count & (BlockSize - 1)] = e;
        Count++;
        if (Count > Peak) Peak = Count;
        return true;
    }

    internal ref OpenInterval Top => ref _blocks[(int)((Count - 1) >> BlockBits)][(Count - 1) & (BlockSize - 1)];

    internal void Pop() => Count--;
}
