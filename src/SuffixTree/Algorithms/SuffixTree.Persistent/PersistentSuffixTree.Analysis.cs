namespace SuffixTree.Persistent;

public sealed partial class PersistentSuffixTree
{
    /// <inheritdoc />
    public string LongestRepeatedSubstring()
    {
        ThrowIfDisposed();
        if (_cachedLrs != null) return _cachedLrs;

        var (deepest, lrsDepth) = FindDeepestInternalNodeWithDepth(NodeAt(_rootOffset));
        if (deepest.IsNull || deepest.Offset == _rootOffset)
        {
            _cachedLrs = string.Empty;
            return string.Empty;
        }

        int length = lrsDepth;
        int depthFromRoot = length - LengthOf(deepest);
        var nav = CreateNavigator();
        int leafPos = nav.FindAnyLeafPosition(deepest, depthFromRoot);
        string result = leafPos < 0 ? string.Empty : _textSource.Substring(leafPos, length);
        _cachedLrs = result;
        return result;
    }

    /// <inheritdoc />
    public IReadOnlyList<string> GetAllSuffixes()
    {
        ThrowIfDisposed();
        var results = new List<string>();
        foreach (var suffix in EnumerateSuffixes())
            results.Add(suffix);
        return results;
    }

    /// <inheritdoc />
    public IEnumerable<string> EnumerateSuffixes()
    {
        ThrowIfDisposed();
        return EnumerateSuffixesCore(NodeAt(_rootOffset));
    }

    private IEnumerable<string> EnumerateSuffixesCore(PersistentSuffixTreeNode root)
    {
        var stack = new Stack<(long NodeOffset, int Depth)>();
        stack.Push((root.Offset, 0));

        while (stack.Count > 0)
        {
            var (nodeOffset, depth) = stack.Pop();
            var node = NodeAt(nodeOffset);
            int nodeDepth = depth + (node.Offset == _rootOffset ? 0 : LengthOf(node));

            if (node.IsLeaf)
            {
                int suffixStart = (_textSource.Length + 1) - nodeDepth;
                if (suffixStart >= 0 && suffixStart < _textSource.Length)
                {
                    int suffixLen = _textSource.Length - suffixStart;
                    yield return _textSource.Substring(suffixStart, suffixLen);
                }
                continue;
            }

            var (arrayBase, entryLayout, childCount) = ReadChildArrayInfo(node);
            // Push children in reverse order so the stack pops them in ascending (lex) order
            for (int ci = childCount - 1; ci >= 0; ci--)
            {
                long entryOff = arrayBase + (long)ci * entryLayout.ChildEntrySize;
                long childNodeOffset = entryLayout.ReadOffset(_storage, entryOff + NodeLayout.ChildOffsetNode);
                stack.Push((childNodeOffset, nodeDepth));
            }
        }
    }

    /// <inheritdoc />
    public string LongestCommonSubstring(string other)
    {
        ThrowIfDisposed();
        var (substring, _, _) = LongestCommonSubstringInfo(other);
        return substring;
    }

    /// <inheritdoc />
    public string LongestCommonSubstring(ReadOnlySpan<char> other)
    {
        ThrowIfDisposed();
        var (substring, _, _) = LongestCommonSubstringInfo(new string(other));
        return substring;
    }

    /// <inheritdoc />
    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
    {
        ThrowIfDisposed();
        var results = FindAllLcsInternal(other, firstOnly: true);
        if (results.PositionsInText.Count == 0)
            return (string.Empty, -1, -1);
        return (results.Substring, results.PositionsInText[0], results.PositionsInOther[0]);
    }

    /// <inheritdoc />
    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(ReadOnlySpan<char> other)
    {
        ThrowIfDisposed();
        return LongestCommonSubstringInfo(new string(other));
    }

    /// <inheritdoc />
    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
    {
        ThrowIfDisposed();
        var results = FindAllLcsInternal(other, firstOnly: false);
        return (results.Substring, results.PositionsInText, results.PositionsInOther);
    }

    /// <summary>
    /// O(m) LCS using suffix-link-based streaming — delegates to shared <see cref="SuffixTreeAlgorithms"/>.
    /// </summary>
    private (string Substring, List<int> PositionsInText, List<int> PositionsInOther) FindAllLcsInternal(string other, bool firstOnly)
    {
        var nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindAllLcs<PersistentSuffixTreeNode, PersistentSuffixTreeNavigator>(ref nav, other, firstOnly);
    }

    /// <inheritdoc/>
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(
        string query, int minLength)
    {
        ThrowIfDisposed();
        var nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindExactMatchAnchors<PersistentSuffixTreeNode, PersistentSuffixTreeNavigator>(ref nav, query, minLength);
    }

    private (PersistentSuffixTreeNode Node, int Depth) FindDeepestInternalNodeWithDepth(PersistentSuffixTreeNode root)
    {
        // Use pre-computed offset + depth from builder/header if available (O(1) path)
        if (_deepestInternalNodeOffset != PersistentConstants.NULL_OFFSET
            && _deepestInternalNodeOffset != _rootOffset
            && _lrsDepth > 0)
        {
            return (NodeAt(_deepestInternalNodeOffset), _lrsDepth);
        }

        if (root.IsLeaf) return (PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout), 0);

        // Fallback: DFS with on-the-fly depth accumulation (works for all formats including v6 Slim)
        var deepest = root;
        int maxDepth = 0;

        // (nodeOffset, depthFromRoot)
        var stack = new Stack<(long NodeOffset, int Depth)>();
        stack.Push((root.Offset, 0));

        while (stack.Count > 0)
        {
            var (nodeOffset, depthFromRoot) = stack.Pop();
            var node = NodeAt(nodeOffset);

            if (!node.IsLeaf)
            {
                int nodeDepth = depthFromRoot + LengthOf(node);
                if (nodeDepth > maxDepth)
                {
                    maxDepth = nodeDepth;
                    deepest = node;
                }

                var (arrayBase, entryLayout, childCount) = ReadChildArrayInfo(node);
                for (int ci = childCount - 1; ci >= 0; ci--)
                {
                    long entryOffset = arrayBase + (long)ci * entryLayout.ChildEntrySize;
                    long childOffset = entryLayout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                    stack.Push((childOffset, nodeDepth));
                }
            }
        }

        return (deepest, maxDepth);
    }
}
