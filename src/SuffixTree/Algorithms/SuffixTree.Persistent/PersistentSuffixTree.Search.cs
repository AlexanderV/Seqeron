namespace SuffixTree.Persistent;

public sealed partial class PersistentSuffixTree
{
    /// <inheritdoc />
    public bool Contains(string value)
    {
        ThrowIfDisposed();
        SuffixTreeSearchContracts.EnsureNotNull(value, nameof(value));
        if (value.Length == 0) return true;
        return Contains(value.AsSpan());
    }

    /// <inheritdoc />
    public bool Contains(ReadOnlySpan<char> value)
    {
        ThrowIfDisposed();
        if (value.IsEmpty) return true;
        var (_, matched, _) = MatchPatternCore(value);
        return matched;
    }

    /// <inheritdoc />
    public IReadOnlyList<int> FindAllOccurrences(string pattern)
    {
        ThrowIfDisposed();
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return FindAllOccurrences(pattern.AsSpan());
    }

    /// <inheritdoc />
    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
    {
        ThrowIfDisposed();
        if (pattern.IsEmpty)
            return SuffixTreeSearchContracts.BuildAllStartPositions(_textSource.Length);

        var (node, matched, depthFromRoot) = MatchPatternCore(pattern);
        if (!matched) return Array.Empty<int>();

        var results = new List<int>((int)node.LeafCount);
        CollectLeaves(node, depthFromRoot, results);
        return results;
    }

    /// <inheritdoc />
    public int CountOccurrences(string pattern)
    {
        ThrowIfDisposed();
        SuffixTreeSearchContracts.EnsureNotNull(pattern, nameof(pattern));
        return CountOccurrences(pattern.AsSpan());
    }

    /// <inheritdoc />
    public int CountOccurrences(ReadOnlySpan<char> pattern)
    {
        ThrowIfDisposed();
        if (pattern.IsEmpty) return _textSource.Length;
        var (node, matched, _) = MatchPatternCore(pattern);
        return matched ? (int)node.LeafCount : 0;
    }

    // Internal helpers
    private (PersistentSuffixTreeNode node, bool matched, int depthFromRoot) MatchPatternCore(ReadOnlySpan<char> pattern)
    {
        var node = NodeAt(_rootOffset);
        int depthFromRoot = 0; // depth to START of node's edge
        int i = 0;
        while (i < pattern.Length)
        {
            if (!TryGetChildOf(node, (uint)pattern[i], out var child) || child.IsNull)
                return (node, false, depthFromRoot);

            // Read Start/End via fast path (avoids 2 interface dispatches per iteration)
            uint childStart = _hybrid.ReadUInt32Fast(child.Offset);
            uint childEnd = _hybrid.ReadUInt32Fast(child.Offset + 4);
            int edgeStart = (int)childStart;
            int edgeLen = (int)((childEnd == PersistentConstants.BOUNDLESS
                ? (uint)(_textSource.Length + 1) : childEnd) - childStart);
            int remaining = pattern.Length - i;
            int compareLen = edgeLen < remaining ? edgeLen : remaining;

            if (edgeStart + compareLen > _textSource.Length)
                return (node, false, depthFromRoot);

            var patternSlice = pattern.Slice(i, compareLen);
            if (_textSource is ITextPatternMatcher patternMatcher)
            {
                if (!patternMatcher.SequenceEqualAt(edgeStart, patternSlice))
                    return (node, false, depthFromRoot);
            }
            else
            {
                var edgeSpan = _textSource.Slice(edgeStart, compareLen);

                // Use SIMD for longer comparisons (>=8 chars), scalar for short
                if (compareLen >= 8)
                {
                    if (!edgeSpan.SequenceEqual(patternSlice))
                        return (node, false, depthFromRoot);
                }
                else
                {
                    for (int j = 0; j < compareLen; j++)
                    {
                        if (edgeSpan[j] != patternSlice[j])
                            return (node, false, depthFromRoot);
                    }
                }
            }

            i += compareLen;
            int childDFR = depthFromRoot + LengthOf(node);
            node = child;
            depthFromRoot = childDFR;
        }
        return (node, true, depthFromRoot);
    }

    /// <summary>
    /// Zone-aware TryGetChild: reads the child array with the correct entry layout.
    /// Uses direct MMF pointer access when available for zero-overhead binary search.
    /// </summary>
    internal bool TryGetChildOf(PersistentSuffixTreeNode parent, uint key, out PersistentSuffixTreeNode child)
    {
        if (_hybrid.TryGetChildFast(parent.Offset, key, out long childOffset))
        {
            child = NodeAt(childOffset);
            return true;
        }
        child = PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout);
        return false;
    }

    private int LengthOf(PersistentSuffixTreeNode node)
        => (int)((node.End == PersistentConstants.BOUNDLESS ? (uint)(_textSource.Length + 1) : node.End) - node.Start);

    private int GetSymbolAt(int index)
    {
        if (index >= _textSource.Length) return -1;
        return _textSource[index];
    }

    private void CollectLeaves(PersistentSuffixTreeNode node, int depth, List<int> results)
    {
        var nav = CreateNavigator();
        nav.CollectLeaves(node, depth, results);
    }

    private PersistentSuffixTreeNavigator CreateNavigator()
        => new PersistentSuffixTreeNavigator(_storage, _rootOffset, _textSource, _hybrid);
}
