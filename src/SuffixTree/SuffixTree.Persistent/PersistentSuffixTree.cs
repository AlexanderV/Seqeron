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

    public PersistentSuffixTree(IStorageProvider storage, long rootOffset, string text)
        : this(storage, rootOffset, new StringTextSource(text))
    {
    }

    public PersistentSuffixTree(IStorageProvider storage, long rootOffset, ITextSource? textSource = null)
    {
        _storage = storage ?? throw new ArgumentNullException(nameof(storage));
        _rootOffset = rootOffset;

        if (textSource != null)
        {
            _textSource = textSource;
        }
        else
        {
            // Try to load text from storage
            long textOff = _storage.ReadInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF);
            int textLen = _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN);

            // If it's a MappedFileStorageProvider, we can use MemoryMappedTextSource for zero-copy
            if (_storage is MappedFileStorageProvider mappedProvider)
            {
                // We don't want MemoryMappedTextSource to dispose the accessor it doesn't own
                _textSource = new MemoryMappedTextSource(mappedProvider.Accessor, textOff, textLen);
            }
            else
            {
                _textSource = new StringTextSource(LoadStringInternal(textOff, textLen));
            }
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
        return EnumerateSuffixesCore(new PersistentSuffixTreeNode(_storage, _rootOffset), new StringBuilder());
    }

    private IEnumerable<string> EnumerateSuffixesCore(PersistentSuffixTreeNode node, StringBuilder currentPath)
    {
        if (node.IsLeaf)
        {
            yield return currentPath.ToString();
            yield break;
        }

        long currentChildEntryOffset = node.ChildrenHead;
        while (currentChildEntryOffset != PersistentConstants.NULL_OFFSET)
        {
            var entry = new PersistentChildEntry(_storage, currentChildEntryOffset);
            var child = new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset);

            int charsAdded = 0;
            int edgeLen = LengthOf(child);
            for (int i = 0; i < edgeLen; i++)
            {
                int s = GetSymbolAt((int)child.Start + i);
                if (s == -1) break;
                currentPath.Append((char)s);
                charsAdded++;
            }

            foreach (var suffix in EnumerateSuffixesCore(child, currentPath))
                yield return suffix;

            currentPath.Length -= charsAdded;
            currentChildEntryOffset = entry.NextEntryOffset;
        }
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
    public void Traverse(ISuffixTreeVisitor visitor)
    {
        ArgumentNullException.ThrowIfNull(visitor);
        TraverseCore(new PersistentSuffixTreeNode(_storage, _rootOffset), 0, visitor);
    }

    private void TraverseCore(PersistentSuffixTreeNode node, int depth, ISuffixTreeVisitor visitor)
    {
        visitor.VisitNode((int)node.Start, (int)node.End, (int)node.LeafCount, node.ChildCount, depth);

        if (!node.IsLeaf)
        {
            var keys = new List<uint>();
            long currentChildEntryOffset = node.ChildrenHead;
            while (currentChildEntryOffset != PersistentConstants.NULL_OFFSET)
            {
                var entry = new PersistentChildEntry(_storage, currentChildEntryOffset);
                keys.Add(entry.Key);
                currentChildEntryOffset = entry.NextEntryOffset;
            }

            // Sort as signed int to match reference implementation order (terminator -1 first)
            keys.Sort((a, b) => ((int)a).CompareTo((int)b));

            int nodeFullDepth = depth + (node.Offset == _rootOffset ? 0 : LengthOf(node));

            foreach (var key in keys)
            {
                if (node.TryGetChild(key, out var child))
                {
                    visitor.EnterBranch((int)key);
                    TraverseCore(child, nodeFullDepth, visitor);
                    visitor.ExitBranch();
                }
            }
        }
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
        if (node.IsLeaf)
        {
            int suffixLength = depth + LengthOf(node);
            int startPosition = (_textSource.Length + 1) - suffixLength;
            if (startPosition < _textSource.Length)
                results.Add(startPosition);
            return;
        }

        long currentChildEntryOffset = node.ChildrenHead;
        while (currentChildEntryOffset != PersistentConstants.NULL_OFFSET)
        {
            var entry = new PersistentChildEntry(_storage, currentChildEntryOffset);
            CollectLeaves(new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset), depth + LengthOf(node), results);
            currentChildEntryOffset = entry.NextEntryOffset;
        }
    }

    private PersistentSuffixTreeNode FindDeepestInternalNode(PersistentSuffixTreeNode node)
    {
        if (node.IsLeaf) return PersistentSuffixTreeNode.Null(_storage);

        PersistentSuffixTreeNode deepest = node;
        int maxDepth = (int)node.DepthFromRoot + LengthOf(node);

        long currentChildEntryOffset = node.ChildrenHead;
        while (currentChildEntryOffset != PersistentConstants.NULL_OFFSET)
        {
            var entry = new PersistentChildEntry(_storage, currentChildEntryOffset);
            var child = new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset);
            var deepestInChild = FindDeepestInternalNode(child);

            if (!deepestInChild.IsNull)
            {
                int childDepth = (int)deepestInChild.DepthFromRoot + LengthOf(deepestInChild);
                if (childDepth > maxDepth)
                {
                    maxDepth = childDepth;
                    deepest = deepestInChild;
                }
            }
            currentChildEntryOffset = entry.NextEntryOffset;
        }
        return deepest;
    }

    public void Dispose()
    {
        _storage.Dispose();
    }
}
