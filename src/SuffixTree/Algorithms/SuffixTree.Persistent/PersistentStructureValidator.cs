namespace SuffixTree.Persistent;

/// <summary>
/// Opt-in strict structural validation of a persistent suffix tree (ST-S0 WP 0.5).
/// <para>
/// <see cref="PersistentSuffixTree.Load(IStorageProvider)"/> validates only the header. The query loops read
/// nodes and child arrays through raw pointers, so a body that is corrupt behind a valid header can send them
/// outside the mapped file. This validator walks the whole tree once, with an explicit stack and with every
/// offset range-checked <em>before</em> it is dereferenced, and rejects the file with an
/// <see cref="InvalidOperationException"/> ("Invalid storage format: ...") when any of these invariants fails:
/// </para>
/// <list type="bullet">
/// <item>every node and child array lies inside the storage, is 4-byte aligned and is reached exactly once (no cycles, no shared subtrees);</item>
/// <item>child keys are strictly ascending and equal the first symbol of the child's edge (the terminator for a leaf that starts at the text end);</item>
/// <item>edges lie inside the text: internal edges are non-empty and end at or before the text end, leaf edges run to the terminator;</item>
/// <item>leaves are in one-to-one correspondence with the suffixes 0..n of the text, and every leaf count is 1 / the sum of its children's;</item>
/// <item>the header node count equals the number of reachable nodes;</item>
/// <item>suffix links of internal nodes point (directly or through the hybrid jump table) to an internal node or the root;</item>
/// <item>the recorded deepest internal node and LRS depth match the tree.</item>
/// </list>
/// The check guarantees memory-safe traversal. It does not prove that the index matches the text byte for byte
/// (a flipped suffix link to another valid node is not detectable without recomputing the tree); use the SHA-256
/// hash of the serializer export for integrity.
/// Cost: one pass over all nodes with bounds-checked reads, plus a bitmap of one bit per 4 bytes of storage.
/// </summary>
internal sealed class PersistentStructureValidator
{
    private const int HeaderSize = PersistentConstants.HEADER_SIZE_V6;

    private struct Frame
    {
        public long NodeOffset;
        public long ArrayBase;
        public long Depth;          // string depth of the node (length of its path label)
        public long LeafSum;        // leaves counted so far among the children
        public int Count;           // number of children
        public int Index;           // next child to visit
        public int EntrySize;       // bytes per child-array entry
        public bool Entry64;        // entries carry an 8-byte child offset
        public int PreviousKey;
        public uint DeclaredLeaves; // LeafCount stored in the node
    }

    private readonly IStorageProvider _storage;
    private readonly HybridLayout _hybrid;
    private readonly long _rootOffset;
    private readonly ITextSource _text;
    private readonly long _deepestOffset;
    private readonly int _lrsDepth;
    private readonly long _size;
    private readonly int _textLength;
    private readonly ulong[] _nodeSeen;     // one bit per 4 bytes of storage, set at each node start
    private readonly ulong[] _suffixSeen;   // one bit per suffix start 0..n
    private Frame[] _stack = new Frame[64];
    private int _sp;
    private long _nodeCount;
    private long _leafCount;
    private long _maxInternalDepth;
    private long _deepestDepth = -1;

    private PersistentStructureValidator(IStorageProvider storage, in HybridLayout hybrid, long rootOffset,
        ITextSource text, long deepestOffset, int lrsDepth)
    {
        _storage = storage;
        _hybrid = hybrid;
        _rootOffset = rootOffset;
        _text = text;
        _deepestOffset = deepestOffset;
        _lrsDepth = lrsDepth;
        _size = storage.Size;
        _textLength = text.Length;

        long words = ((_size >> 2) + 64) >> 6;
        if (words > Array.MaxLength)
            throw new NotSupportedException(
                $"Storage of {_size} bytes is too large for strict structural validation (limit {(long)Array.MaxLength * 256} bytes).");
        _nodeSeen = new ulong[words];
        _suffixSeen = new ulong[((long)_textLength + 64) >> 6];
    }

    /// <summary>Validates the tree; throws <see cref="InvalidOperationException"/> if any invariant fails.</summary>
    internal static void Validate(IStorageProvider storage, in HybridLayout hybrid, long rootOffset,
        ITextSource text, long deepestOffset, int lrsDepth)
    {
        ArgumentNullException.ThrowIfNull(storage);
        ArgumentNullException.ThrowIfNull(text);
        try
        {
            new PersistentStructureValidator(storage, hybrid, rootOffset, text, deepestOffset, lrsDepth).Run();
        }
        catch (ArgumentOutOfRangeException ex)
        {
            // Every offset is range-checked before use, so this is a safety net, not the expected path.
            throw new InvalidOperationException("Invalid storage format: a read fell outside the storage: " + ex.Message, ex);
        }
    }

    private static InvalidOperationException Invalid(string message) =>
        new("Invalid storage format: " + message);

    private void Run()
    {
        WalkTree();
        CheckGlobalCounts();
        CheckSuffixLinks();
        CheckDeepestNode();
    }

    // ───────────────────────── phase 1: depth-first walk ─────────────────────────

    private void WalkTree()
    {
        if (Enter(_rootOffset, parentDepth: 0, key: 0, isRoot: true, out Frame root))
            Push(root);

        while (_sp > 0)
        {
            ref Frame f = ref _stack[_sp - 1];

            if (f.Index == f.Count)
            {
                if (f.LeafSum != f.DeclaredLeaves)
                    throw Invalid($"node at offset {f.NodeOffset} declares {f.DeclaredLeaves} leaves but its children contain {f.LeafSum}.");
                if (f.NodeOffset != _rootOffset && f.Depth > _maxInternalDepth)
                    _maxInternalDepth = f.Depth;
                if (f.NodeOffset == _deepestOffset)
                    _deepestDepth = f.Depth;
                long finished = f.LeafSum;
                _sp--;
                if (_sp > 0)
                    _stack[_sp - 1].LeafSum += finished;
                continue;
            }

            long entryOffset = f.ArrayBase + (long)f.Index * f.EntrySize;
            uint rawKey = _storage.ReadUInt32(entryOffset + NodeLayout.ChildOffsetKey);
            long child = f.Entry64
                ? _storage.ReadInt64(entryOffset + NodeLayout.ChildOffsetNode)
                : ReadCompactOffset(entryOffset + NodeLayout.ChildOffsetNode);
            int key = unchecked((int)rawKey);
            if (f.Index > 0 && key <= f.PreviousKey)
                throw Invalid($"child keys of the node at offset {f.NodeOffset} are not strictly ascending (entry {f.Index}: {key} after {f.PreviousKey}).");
            f.PreviousKey = key;
            f.Index++;
            long parentDepth = f.Depth;

            if (Enter(child, parentDepth, key, isRoot: false, out Frame next))
                Push(next);
            else
                _stack[_sp - 1].LeafSum += 1;
        }
    }

    private long ReadCompactOffset(long position)
    {
        uint raw = _storage.ReadUInt32(position);
        return raw == uint.MaxValue ? PersistentConstants.NULL_OFFSET : raw;
    }

    private void Push(in Frame frame)
    {
        if (_sp == _stack.Length)
        {
            if (_stack.Length >= Array.MaxLength / 2)
                throw new NotSupportedException("The tree is too deep for strict structural validation.");
            Array.Resize(ref _stack, _stack.Length * 2);
        }
        _stack[_sp++] = frame;
    }

    /// <summary>
    /// Validates the node at <paramref name="offset"/>. Returns <see langword="true"/> with a frame for an
    /// internal node whose children must be visited; <see langword="false"/> for a validated leaf.
    /// </summary>
    private bool Enter(long offset, long parentDepth, int key, bool isRoot, out Frame frame)
    {
        frame = default;
        NodeLayout layout = _hybrid.LayoutForOffset(offset);
        if (offset < HeaderSize || (offset & 3) != 0 || offset > _size - layout.NodeSize)
            throw Invalid($"node offset {offset} is outside the node area [{HeaderSize}, {_size}) or misaligned.");
        if (IsSeen(offset))
            throw Invalid($"node at offset {offset} is reachable more than once (cycle or shared subtree).");
        MarkSeen(offset);
        _nodeCount++;

        uint start = _storage.ReadUInt32(offset + PersistentConstants.OFFSET_START);
        uint end = _storage.ReadUInt32(offset + PersistentConstants.OFFSET_END);
        uint leafCount = _storage.ReadUInt32(offset + layout.OffsetLeafCount);
        int rawCount = _storage.ReadInt32(offset + layout.OffsetChildCount);
        bool isLeaf = end == PersistentConstants.BOUNDLESS;
        long n = _textLength;

        if (isRoot && isLeaf)
            throw Invalid("the root is marked as a leaf.");

        long edgeLength;
        if (isRoot)
        {
            edgeLength = 0;
        }
        else if (isLeaf)
        {
            if (start > n)
                throw Invalid($"leaf at offset {offset} starts at {start}, beyond the text end {n}.");
            edgeLength = n + 1 - start;
        }
        else
        {
            if (start >= end || end > n)
                throw Invalid($"internal node at offset {offset} has the invalid edge [{start}, {end}) for a text of length {n}.");
            edgeLength = (long)end - start;
        }

        long depth = parentDepth + edgeLength;
        if (depth > n + 1)
            throw Invalid($"node at offset {offset} has string depth {depth}, more than the text length plus terminator {n + 1}.");

        if (!isRoot)
        {
            // The first symbol of the edge must be the key under which the parent lists the child.
            int symbol = start < n ? _text[(int)start] : -1;
            if (symbol != key)
                throw Invalid($"child at offset {offset} is listed under key {key} but its edge starts with symbol {symbol}.");
        }

        if (isLeaf)
        {
            if (rawCount != 0)
                throw Invalid($"leaf at offset {offset} has {rawCount} children.");
            if (leafCount != 1)
                throw Invalid($"leaf at offset {offset} declares {leafCount} leaves.");
            long suffixStart = start - parentDepth;
            if (suffixStart < 0 || suffixStart > n)
                throw Invalid($"leaf at offset {offset} (edge start {start}, parent depth {parentDepth}) does not correspond to a suffix of the text.");
            if (!MarkSuffix(suffixStart))
                throw Invalid($"suffix {suffixStart} is represented by more than one leaf (second at offset {offset}).");
            _leafCount++;
            return false;
        }

        if (!isRoot && end < depth)
            throw Invalid($"internal node at offset {offset} has string depth {depth} but its edge ends at text position {end}.");

        bool jumped = (rawCount & int.MinValue) != 0;
        int count = rawCount & int.MaxValue;
        if (count == 0 && !isRoot)
            throw Invalid($"internal node at offset {offset} has no children.");

        long head = layout.ReadOffset(_storage, offset + layout.OffsetChildrenHead);
        long arrayBase;
        int entrySize;
        bool entry64;
        if (jumped)
        {
            if (!_hybrid.IsHybrid || offset >= _hybrid.TransitionOffset)
                throw Invalid($"node at offset {offset} has a hybrid-jump child array outside the compact zone.");
            if (head < _hybrid.JumpTableStart || head > _hybrid.JumpTableEnd - 8 || ((head - _hybrid.JumpTableStart) & 7) != 0)
                throw Invalid($"node at offset {offset} refers to jump-table slot {head} outside [{_hybrid.JumpTableStart}, {_hybrid.JumpTableEnd}).");
            arrayBase = _storage.ReadInt64(head);
            entrySize = NodeLayout.Large.ChildEntrySize;
            entry64 = true;
        }
        else
        {
            arrayBase = head;
            entrySize = layout.ChildEntrySize;
            entry64 = layout.OffsetIs64Bit;
        }

        if (count > 0)
        {
            long arrayBytes = (long)count * entrySize;
            if (arrayBase < HeaderSize || (arrayBase & 3) != 0 || arrayBase > _size - arrayBytes)
                throw Invalid($"child array of the node at offset {offset} ({count} entries at {arrayBase}) lies outside the storage [{HeaderSize}, {_size}) or is misaligned.");
        }

        frame = new Frame
        {
            NodeOffset = offset,
            ArrayBase = arrayBase,
            Depth = depth,
            Count = count,
            EntrySize = entrySize,
            Entry64 = entry64,
            PreviousKey = int.MinValue,
            DeclaredLeaves = leafCount,
        };
        return true;
    }

    // ───────────────────────── phase 2: global checks ─────────────────────────

    private void CheckGlobalCounts()
    {
        int headerNodeCount = _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT);
        if (headerNodeCount != _nodeCount)
            throw Invalid($"the header records {headerNodeCount} nodes but {_nodeCount} are reachable from the root.");

        long expectedLeaves = _textLength == 0 ? 0 : _textLength + 1L;
        if (_leafCount != expectedLeaves)
            throw Invalid($"the tree has {_leafCount} leaves but a text of length {_textLength} needs {expectedLeaves}.");
    }

    private void CheckSuffixLinks()
    {
        for (long w = 0; w < _nodeSeen.LongLength; w++)
        {
            ulong bits = _nodeSeen[w];
            while (bits != 0)
            {
                int b = System.Numerics.BitOperations.TrailingZeroCount(bits);
                bits &= bits - 1;
                long offset = ((w << 6) + b) << 2;
                NodeLayout layout = _hybrid.LayoutForOffset(offset);
                uint end = _storage.ReadUInt32(offset + PersistentConstants.OFFSET_END);
                if (end == PersistentConstants.BOUNDLESS)
                    continue; // leaves carry no suffix link

                long link = layout.ReadOffset(_storage, offset + layout.OffsetSuffixLink);
                if (link == PersistentConstants.NULL_OFFSET)
                    continue;

                if (_hybrid.JumpTableStart >= 0 && link >= _hybrid.JumpTableStart && link < _hybrid.JumpTableEnd)
                {
                    if (link > _hybrid.JumpTableEnd - 8 || ((link - _hybrid.JumpTableStart) & 7) != 0)
                        throw Invalid($"suffix link of the node at offset {offset} points into the middle of a jump-table slot ({link}).");
                    link = _storage.ReadInt64(link);
                }

                if (link < HeaderSize || (link & 3) != 0 || link >= _size || !IsSeen(link))
                    throw Invalid($"suffix link of the node at offset {offset} ({link}) does not point to a node of the tree.");
                if (_storage.ReadUInt32(link + PersistentConstants.OFFSET_END) == PersistentConstants.BOUNDLESS)
                    throw Invalid($"suffix link of the node at offset {offset} points to a leaf ({link}).");
            }
        }
    }

    private void CheckDeepestNode()
    {
        if (_deepestOffset == PersistentConstants.NULL_OFFSET)
            return;

        if (_deepestDepth < 0)
            throw Invalid($"the recorded deepest internal node ({_deepestOffset}) is not an internal node of the tree.");
        if (_deepestDepth != _maxInternalDepth)
            throw Invalid($"the recorded deepest internal node ({_deepestOffset}) has depth {_deepestDepth}, but the deepest internal node has depth {_maxInternalDepth}.");
        if (_lrsDepth >= 0 && _lrsDepth != _deepestDepth)
            throw Invalid($"the recorded longest-repeat depth {_lrsDepth} differs from the depth of the deepest internal node {_deepestDepth}.");
    }

    // ───────────────────────── bitmaps ─────────────────────────

    private bool IsSeen(long offset)
    {
        long bit = offset >> 2;
        return (_nodeSeen[bit >> 6] & (1UL << (int)(bit & 63))) != 0;
    }

    private void MarkSeen(long offset)
    {
        long bit = offset >> 2;
        _nodeSeen[bit >> 6] |= 1UL << (int)(bit & 63);
    }

    /// <summary>Marks suffix <paramref name="start"/> as represented; returns false if it already was.</summary>
    private bool MarkSuffix(long start)
    {
        ulong mask = 1UL << (int)(start & 63);
        ref ulong word = ref _suffixSeen[start >> 6];
        if ((word & mask) != 0)
            return false;
        word |= mask;
        return true;
    }
}
