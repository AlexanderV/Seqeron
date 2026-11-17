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
internal sealed class IntervalStack(long capacity, int blockBits = 20)
{
    internal const uint None = uint.MaxValue;
    private readonly int _blockBits = blockBits;
    private readonly int _blockSize = 1 << blockBits;
    private readonly List<OpenInterval[]> _blocks = [];

    internal long Capacity { get; } = capacity;

    internal long Count { get; private set; }

    internal long Peak { get; private set; }

    /// <summary>False when the push would exceed the capacity (nothing is pushed).</summary>
    internal bool TryPush(in OpenInterval e)
    {
        if (Count >= Capacity) return false;
        int block = (int)(Count >> _blockBits);
        if (block == _blocks.Count) _blocks.Add(new OpenInterval[_blockSize]);
        _blocks[block][Count & (_blockSize - 1)] = e;
        Count++;
        if (Count > Peak) Peak = Count;
        return true;
    }

    internal ref OpenInterval Top => ref _blocks[(int)((Count - 1) >> _blockBits)][(Count - 1) & (_blockSize - 1)];

    internal void Pop() => Count--;
}
