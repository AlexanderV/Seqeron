using System;
using System.Collections.Generic;
using System.Text;

namespace SuffixTree.Persistent;

/// <summary>
/// A persistent implementation of ISuffixTree that operates on a storage provider.
/// </summary>
public class PersistentSuffixTree : ISuffixTree, IDisposable
{
    private readonly IStorageProvider _storage;
    private readonly long _rootOffset;
    private readonly ITextSource _textSource;
    private readonly bool _ownsTextSource;
    private bool _disposed;

    public PersistentSuffixTree(IStorageProvider storage, long rootOffset, ITextSource? textSource = null)
    {
        _storage = storage ?? throw new ArgumentNullException(nameof(storage));
        _rootOffset = rootOffset;

        if (textSource != null)
        {
            _textSource = textSource;
            _ownsTextSource = false; // Caller retains ownership
        }
        else
        {
            // Try to load text from storage
            long textOff = _storage.ReadInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF);
            int textLen = _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN);

            // If it's a MappedFileStorageProvider, we can use MemoryMappedTextSource for zero-copy
            if (_storage is MappedFileStorageProvider mappedProvider)
            {
                _textSource = new MemoryMappedTextSource(mappedProvider.Accessor, textOff, textLen);
            }
            else
            {
                _textSource = new StringTextSource(LoadStringInternal(textOff, textLen));
            }
            _ownsTextSource = true; // We created it, we dispose it
        }
    }

    private string LoadStringInternal(long textOff, int textLen)
    {
        var sb = new StringBuilder(textLen);
        for (int i = 0; i < textLen; i++)
        {
            sb.Append(_storage.ReadChar(textOff + (i * 2)));
        }
        return sb.ToString();
    }

    public static PersistentSuffixTree Load(IStorageProvider storage)
    {
        long magic = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_MAGIC);
        if (magic != PersistentConstants.MAGIC_NUMBER)
            throw new InvalidOperationException("Invalid storage format: Magic number mismatch.");

        long root = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_ROOT);
        return new PersistentSuffixTree(storage, root);
    }

    public ITextSource Text => _textSource;

    public int NodeCount => _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT);

    public int LeafCount
    {
        get
        {
            uint rawCount = new PersistentSuffixTreeNode(_storage, _rootOffset).LeafCount;
            return rawCount > 0 ? (int)(rawCount - 1) : 0;
        }
    }

    public int MaxDepth => _textSource.Length + 1;

    public bool IsEmpty => _textSource.Length == 0;

    public bool Contains(string value)
    {
        if (string.IsNullOrEmpty(value)) return true;
        return Contains(value.AsSpan());
    }

    public bool Contains(ReadOnlySpan<char> value)
    {
        if (value.IsEmpty) return true;
        var (node, matched) = MatchPatternCore(value);
        return matched;
    }

    public IReadOnlyList<int> FindAllOccurrences(string pattern)
    {
        if (pattern == null) throw new ArgumentNullException(nameof(pattern));
        return FindAllOccurrences(pattern.AsSpan());
    }

    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
    {
        var results = new List<int>();
        if (pattern.IsEmpty) return results;

        var (node, matched) = MatchPatternCore(pattern);
        if (!matched) return results;

        CollectLeaves(node, (int)node.DepthFromRoot, results);
        return results;
    }

    public int CountOccurrences(string pattern)
    {
        if (pattern == null) throw new ArgumentNullException(nameof(pattern));
        return CountOccurrences(pattern.AsSpan());
    }

    public int CountOccurrences(ReadOnlySpan<char> pattern)
    {
        if (pattern.IsEmpty) return 0;
        var (node, matched) = MatchPatternCore(pattern);
        return matched ? (int)node.LeafCount : 0;
    }

    public string LongestRepeatedSubstring()
    {
        // For simplicity, reusing logic or implementing LRS
        // This usually requires a DFS on the tree to find deepest internal node
        var deepest = FindDeepestInternalNode(new PersistentSuffixTreeNode(_storage, _rootOffset));
        if (deepest.IsNull || deepest.Offset == _rootOffset) return string.Empty;

        int length = (int)deepest.DepthFromRoot + LengthOf(deepest);
        // Find one occurrence to get the text
        var occurrences = new List<int>();
        CollectLeaves(deepest, (int)deepest.DepthFromRoot, occurrences);
        if (occurrences.Count == 0) return string.Empty;

        return _textSource.Substring(occurrences[0], length);
    }

    public IReadOnlyList<string> GetAllSuffixes()
    {
        var results = new List<string>();
        foreach (var suffix in EnumerateSuffixes())
            results.Add(suffix);
        return results;
    }

    public IEnumerable<string> EnumerateSuffixes()
    {
        // Lazy enumeration of suffixes (simplified version)
        return EnumerateSuffixesCore(new PersistentSuffixTreeNode(_storage, _rootOffset));
    }

    private IEnumerable<string> EnumerateSuffixesCore(PersistentSuffixTreeNode root)
    {
        var results = new List<string>();
        var stack = new Stack<(PersistentSuffixTreeNode Node, int Depth)>();
        stack.Push((root, 0));

        while (stack.Count > 0)
        {
            var (node, depth) = stack.Pop();
            int nodeDepth = depth + (node.Offset == _rootOffset ? 0 : LengthOf(node));

            if (node.IsLeaf)
            {
                int suffixStart = (_textSource.Length + 1) - nodeDepth;
                if (suffixStart >= 0 && suffixStart < _textSource.Length)
                {
                    int suffixLen = _textSource.Length - suffixStart;
                    results.Add(_textSource.Substring(suffixStart, suffixLen));
                }
                continue;
            }

            long childEntry = node.ChildrenHead;
            while (childEntry != PersistentConstants.NULL_OFFSET)
            {
                var entry = new PersistentChildEntry(_storage, childEntry);
                stack.Push((new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset), nodeDepth));
                childEntry = entry.NextEntryOffset;
            }
        }

        return results;
    }

    public string LongestCommonSubstring(string other) => LongestCommonSubstring(other.AsSpan());

    public string LongestCommonSubstring(ReadOnlySpan<char> other)
    {
        var (substring, _, _) = LongestCommonSubstringInfo(new string(other));
        return substring;
    }

    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
    {
        // LCS implementation for ISuffixTree
        int maxLength = 0;
        int maxPosText = -1;
        int maxPosOther = -1;

        var node = new PersistentSuffixTreeNode(_storage, _rootOffset);

        for (int i = 0; i < other.Length; i++)
        {
            var (matchNode, matchLen) = MatchLongestPrefix(other.AsSpan().Slice(i), node);
            if (matchLen > maxLength)
            {
                maxLength = matchLen;
                maxPosOther = i;

                // Find position in text
                var occurrences = new List<int>();
                CollectLeaves(matchNode, (int)matchNode.DepthFromRoot, occurrences);
                if (occurrences.Count > 0)
                {
                    maxPosText = occurrences[0];
                }
            }
        }

        if (maxLength == 0) return (string.Empty, -1, -1);
        return (other.Substring(maxPosOther, maxLength), maxPosText, maxPosOther);
    }

    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(ReadOnlySpan<char> other)
    {
        // LCS implementation for ISuffixTree
        int maxLength = 0;
        int maxPosText = -1;
        int maxPosOther = -1;

        var node = new PersistentSuffixTreeNode(_storage, _rootOffset);

        for (int i = 0; i < other.Length; i++)
        {
            var (matchNode, matchLen) = MatchLongestPrefix(other.Slice(i), node);
            if (matchLen > maxLength)
            {
                maxLength = matchLen;
                maxPosOther = i;

                // Find position in text
                var occurrences = new List<int>();
                CollectLeaves(matchNode, (int)matchNode.DepthFromRoot, occurrences);
                if (occurrences.Count > 0)
                {
                    maxPosText = occurrences[0];
                }
            }
        }

        if (maxLength == 0) return (string.Empty, -1, -1);
        return (new string(other.Slice(maxPosOther, maxLength)), maxPosText, maxPosOther);
    }

    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
    {
        // Simplified for now
        var info = LongestCommonSubstringInfo(other);
        if (info.PositionInText == -1) return (string.Empty, Array.Empty<int>(), Array.Empty<int>());
        return (info.Substring, new[] { info.PositionInText }, new[] { info.PositionInOther });
    }

    public string PrintTree() => "Persistent Tree Visualization Not Implemented";

    /// <inheritdoc/>
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(
        string query, int minLength)
    {
        ArgumentNullException.ThrowIfNull(query);
        if (query.Length == 0 || _textSource.Length == 0 || minLength <= 0)
            return Array.Empty<(int, int, int)>();

        var results = new List<(int PositionInText, int PositionInQuery, int Length)>();
        var root = new PersistentSuffixTreeNode(_storage, _rootOffset);

        var currentNode = root;
        PersistentSuffixTreeNode currentEdge = PersistentSuffixTreeNode.Null(_storage);
        int edgeOffset = 0;
        int currentMatchLen = 0;

        // Peak tracking
        int peakLen = 0;
        int peakEndInQuery = -1;
        PersistentSuffixTreeNode peakNode = PersistentSuffixTreeNode.Null(_storage);

        for (int i = 0; i < query.Length; i++)
        {
            uint c = (uint)query[i];

            while (true)
            {
                if (!currentEdge.IsNull)
                {
                    if (GetSymbolAt((int)(currentEdge.Start + (uint)edgeOffset)) == (int)c)
                    {
                        edgeOffset++;
                        currentMatchLen++;
                        if (edgeOffset >= LengthOf(currentEdge))
                        {
                            currentNode = currentEdge;
                            currentEdge = PersistentSuffixTreeNode.Null(_storage);
                            edgeOffset = 0;
                        }
                        break;
                    }
                }
                else
                {
                    if (currentNode.TryGetChild(c, out var nextChild) && !nextChild.IsNull)
                    {
                        currentEdge = nextChild;
                        edgeOffset = 1;
                        currentMatchLen++;
                        if (edgeOffset >= LengthOf(currentEdge))
                        {
                            currentNode = currentEdge;
                            currentEdge = PersistentSuffixTreeNode.Null(_storage);
                            edgeOffset = 0;
                        }
                        break;
                    }
                }

                // Cannot extend — follow suffix link
                if (currentMatchLen == 0) break;

                if (currentNode.Offset != _rootOffset)
                {
                    long suffixLink = currentNode.SuffixLink;
                    currentNode = suffixLink != PersistentConstants.NULL_OFFSET
                        ? new PersistentSuffixTreeNode(_storage, suffixLink)
                        : root;
                }
                currentMatchLen--;

                int nodeDepth = GetNodeDepth(currentNode);
                int remaining = currentMatchLen - nodeDepth;

                if (remaining > 0)
                {
                    int pos = i - remaining;
                    currentEdge = PersistentSuffixTreeNode.Null(_storage);
                    edgeOffset = 0;

                    while (remaining > 0)
                    {
                        if (!currentNode.TryGetChild((uint)query[pos], out var nextChild2) || nextChild2.IsNull)
                            break;

                        int edgeLen = LengthOf(nextChild2);
                        if (edgeLen <= remaining)
                        {
                            pos += edgeLen;
                            remaining -= edgeLen;
                            currentNode = nextChild2;
                        }
                        else
                        {
                            currentEdge = nextChild2;
                            edgeOffset = remaining;
                            remaining = 0;
                        }
                    }
                }
                else
                {
                    currentEdge = PersistentSuffixTreeNode.Null(_storage);
                    edgeOffset = 0;
                }
            }

            // Update peak tracking
            if (currentMatchLen >= minLength)
            {
                if (currentMatchLen > peakLen)
                {
                    peakLen = currentMatchLen;
                    peakEndInQuery = i;
                    peakNode = currentEdge.IsNull ? currentNode : currentEdge;
                }
            }
            else if (peakLen >= minLength)
            {
                EmitAnchorFromPeak(results, peakNode, peakEndInQuery, peakLen);
                peakLen = 0;
                peakEndInQuery = -1;
                peakNode = PersistentSuffixTreeNode.Null(_storage);
            }
        }

        // Emit final run
        if (peakLen >= minLength && !peakNode.IsNull)
        {
            EmitAnchorFromPeak(results, peakNode, peakEndInQuery, peakLen);
        }

        return results;
    }

    private void EmitAnchorFromPeak(
        List<(int PositionInText, int PositionInQuery, int Length)> results,
        PersistentSuffixTreeNode node,
        int endInQuery,
        int length)
    {
        int refPos = FindAnyLeafPosition(node);
        if (refPos >= 0)
        {
            results.Add((refPos, endInQuery - length + 1, length));
        }
    }

    private int FindAnyLeafPosition(PersistentSuffixTreeNode node)
    {
        var current = node;
        while (!current.IsLeaf)
        {
            // Walk to any child (first one available)
            long childEntryOffset = current.ChildrenHead;
            PersistentSuffixTreeNode bestChild = PersistentSuffixTreeNode.Null(_storage);
            while (childEntryOffset != PersistentConstants.NULL_OFFSET)
            {
                var entry = new PersistentChildEntry(_storage, childEntryOffset);
                bestChild = new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset);
                // Prefer non-terminator children
                if (entry.Key != PersistentConstants.TERMINATOR_KEY)
                    break;
                childEntryOffset = entry.NextEntryOffset;
            }
            if (bestChild.IsNull) return -1;
            current = bestChild;
        }

        int leafDepth = GetNodeDepth(current);
        int pos = (_textSource.Length + 1) - leafDepth;
        return (pos >= 0 && pos < _textSource.Length) ? pos : -1;
    }

    private int GetNodeDepth(PersistentSuffixTreeNode node)
    {
        if (node.Offset == _rootOffset) return 0;
        return (int)node.DepthFromRoot + LengthOf(node);
    }

    /// <inheritdoc/>
    public void Traverse(ISuffixTreeVisitor visitor)
    {
        ArgumentNullException.ThrowIfNull(visitor);
        TraverseCore(new PersistentSuffixTreeNode(_storage, _rootOffset), 0, visitor);
    }

    private void TraverseCore(PersistentSuffixTreeNode node, int depth, ISuffixTreeVisitor visitor)
    {
        // Frame: node being processed, sorted child keys, current index into keys, depth
        var stack = new Stack<(PersistentSuffixTreeNode Node, List<uint> Keys, int Index, int Depth)>();

        // Visit root node
        visitor.VisitNode((int)node.Start, (int)node.End, (int)node.LeafCount, node.ChildCount, depth);
        if (!node.IsLeaf)
        {
            var rootKeys = CollectSortedKeys(node);
            int rootFullDepth = depth; // root has no edge
            stack.Push((node, rootKeys, 0, rootFullDepth));
        }

        while (stack.Count > 0)
        {
            var (current, keys, index, parentDepth) = stack.Pop();

            if (index >= keys.Count)
            {
                // All children processed — exit branch (unless this is the root frame)
                if (stack.Count > 0)
                    visitor.ExitBranch();
                continue;
            }

            // Push continuation for next sibling
            stack.Push((current, keys, index + 1, parentDepth));

            var key = keys[index];
            if (current.TryGetChild(key, out var child))
            {
                visitor.EnterBranch((int)key);
                visitor.VisitNode((int)child.Start, (int)child.End, (int)child.LeafCount, child.ChildCount, parentDepth + 1);

                if (!child.IsLeaf)
                {
                    var childKeys = CollectSortedKeys(child);
                    int childFullDepth = parentDepth + LengthOf(child);
                    stack.Push((child, childKeys, 0, childFullDepth));
                }
                else
                {
                    visitor.ExitBranch();
                }
            }
        }
    }

    private List<uint> CollectSortedKeys(PersistentSuffixTreeNode node)
    {
        var keys = new List<uint>();
        long offset = node.ChildrenHead;
        while (offset != PersistentConstants.NULL_OFFSET)
        {
            var entry = new PersistentChildEntry(_storage, offset);
            keys.Add(entry.Key);
            offset = entry.NextEntryOffset;
        }
        keys.Sort((a, b) => ((int)a).CompareTo((int)b));
        return keys;
    }

    // Internal helpers
    private (PersistentSuffixTreeNode node, bool matched) MatchPatternCore(ReadOnlySpan<char> pattern)
    {
        var node = new PersistentSuffixTreeNode(_storage, _rootOffset);
        int i = 0;
        while (i < pattern.Length)
        {
            if (!node.TryGetChild((uint)pattern[i], out var child) || child.IsNull)
                return (node, false);

            int edgeLen = LengthOf(child);
            int remaining = pattern.Length - i;
            int compareLen = edgeLen < remaining ? edgeLen : remaining;

            for (int j = 0; j < compareLen; j++)
            {
                if (GetSymbolAt((int)child.Start + j) != pattern[i + j])
                    return (node, false);
            }

            i += compareLen;
            node = child;
        }
        return (node, true);
    }

    private (PersistentSuffixTreeNode node, int length) MatchLongestPrefix(ReadOnlySpan<char> pattern, PersistentSuffixTreeNode startNode)
    {
        var node = startNode;
        int i = 0;
        while (i < pattern.Length)
        {
            if (!node.TryGetChild((uint)pattern[i], out var child) || child.IsNull)
                break;

            int edgeLen = LengthOf(child);
            int remaining = pattern.Length - i;
            int compareLen = edgeLen < remaining ? edgeLen : remaining;

            int j = 0;
            for (; j < compareLen; j++)
            {
                if (GetSymbolAt((int)child.Start + j) != pattern[i + j])
                    break;
            }

            i += j;
            node = child;
            if (j < compareLen) break;
        }
        return (node, i);
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
        var stack = new Stack<(PersistentSuffixTreeNode Node, int Depth)>();
        stack.Push((node, depth));

        while (stack.Count > 0)
        {
            var (current, currentDepth) = stack.Pop();

            if (current.IsLeaf)
            {
                int suffixLength = currentDepth + LengthOf(current);
                int startPosition = (_textSource.Length + 1) - suffixLength;
                if (startPosition < _textSource.Length)
                    results.Add(startPosition);
                continue;
            }

            int childDepth = currentDepth + LengthOf(current);
            long childEntryOffset = current.ChildrenHead;
            while (childEntryOffset != PersistentConstants.NULL_OFFSET)
            {
                var entry = new PersistentChildEntry(_storage, childEntryOffset);
                stack.Push((new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset), childDepth));
                childEntryOffset = entry.NextEntryOffset;
            }
        }
    }

    private PersistentSuffixTreeNode FindDeepestInternalNode(PersistentSuffixTreeNode root)
    {
        if (root.IsLeaf) return PersistentSuffixTreeNode.Null(_storage);

        var deepest = root;
        int maxDepth = (int)root.DepthFromRoot + LengthOf(root);

        var stack = new Stack<PersistentSuffixTreeNode>();
        stack.Push(root);

        while (stack.Count > 0)
        {
            var node = stack.Pop();

            if (!node.IsLeaf)
            {
                int nodeDepth = (int)node.DepthFromRoot + LengthOf(node);
                if (nodeDepth > maxDepth)
                {
                    maxDepth = nodeDepth;
                    deepest = node;
                }

                // Collect non-leaf children, push in reverse order to preserve
                // left-to-right DFS order (matching recursive version tie-breaking)
                var children = new List<PersistentSuffixTreeNode>();
                long childEntryOffset = node.ChildrenHead;
                while (childEntryOffset != PersistentConstants.NULL_OFFSET)
                {
                    var entry = new PersistentChildEntry(_storage, childEntryOffset);
                    var child = new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset);
                    if (!child.IsLeaf)
                        children.Add(child);
                    childEntryOffset = entry.NextEntryOffset;
                }
                for (int i = children.Count - 1; i >= 0; i--)
                    stack.Push(children[i]);
            }
        }

        return deepest;
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;

        // Dispose textSource BEFORE storage: MemoryMappedTextSource must
        // ReleasePointer() while the underlying accessor is still alive.
        if (_ownsTextSource && _textSource is IDisposable disposableText)
            disposableText.Dispose();

        _storage.Dispose();
    }

    private void ThrowIfDisposed()
    {
        if (_disposed)
            throw new ObjectDisposedException(nameof(PersistentSuffixTree));
    }
}
