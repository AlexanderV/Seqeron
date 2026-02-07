using System;
using System.Buffers;
using System.Collections.Generic;
using System.Text;

namespace SuffixTree.Persistent;

/// <summary>
/// A persistent implementation of ISuffixTree that operates on a storage provider.
/// </summary>
public sealed class PersistentSuffixTree : ISuffixTree, IDisposable
{
    private readonly IStorageProvider _storage;
    private readonly long _rootOffset;
    private readonly ITextSource _textSource;
    private readonly bool _ownsTextSource;
    private volatile bool _disposed;
    private volatile string? _cachedLrs;

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
        long byteLen = (long)textLen * 2;
        if (byteLen > int.MaxValue)
            throw new InvalidOperationException(
                $"Text length {textLen} exceeds maximum loadable size ({int.MaxValue / 2} characters).");
        int byteLenInt = (int)byteLen;
        byte[] bytes = ArrayPool<byte>.Shared.Rent(byteLenInt);
        try
        {
            _storage.ReadBytes(textOff, bytes, 0, byteLenInt);
            return Encoding.Unicode.GetString(bytes, 0, byteLenInt);
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(bytes);
        }
    }

    public static PersistentSuffixTree Load(IStorageProvider storage)
    {
        ArgumentNullException.ThrowIfNull(storage);

        long magic = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_MAGIC);
        if (magic != PersistentConstants.MAGIC_NUMBER)
            throw new InvalidOperationException("Invalid storage format: Magic number mismatch.");

        int version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        if (version != PersistentConstants.CURRENT_VERSION)
            throw new InvalidOperationException($"Unsupported storage format version: {version}. Expected: {PersistentConstants.CURRENT_VERSION}.");

        long root = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_ROOT);
        return new PersistentSuffixTree(storage, root);
    }

    public ITextSource Text { get { ThrowIfDisposed(); return _textSource; } }

    public int NodeCount { get { ThrowIfDisposed(); return _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT); } }

    public int LeafCount
    {
        get
        {
            ThrowIfDisposed();
            uint rawCount = new PersistentSuffixTreeNode(_storage, _rootOffset).LeafCount;
            return rawCount > 0 ? (int)(rawCount - 1) : 0;
        }
    }

    public int MaxDepth { get { ThrowIfDisposed(); return _textSource.Length + 1; } }

    public bool IsEmpty { get { ThrowIfDisposed(); return _textSource.Length == 0; } }

    public bool Contains(string value)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(value);
        if (value.Length == 0) return true;
        return Contains(value.AsSpan());
    }

    public bool Contains(ReadOnlySpan<char> value)
    {
        ThrowIfDisposed();
        if (value.IsEmpty) return true;
        var (node, matched) = MatchPatternCore(value);
        return matched;
    }

    public IReadOnlyList<int> FindAllOccurrences(string pattern)
    {
        ThrowIfDisposed();
        if (pattern == null) throw new ArgumentNullException(nameof(pattern));
        return FindAllOccurrences(pattern.AsSpan());
    }

    public IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern)
    {
        ThrowIfDisposed();
        var results = new List<int>();
        if (pattern.IsEmpty)
        {
            for (int i = 0; i < _textSource.Length; i++) results.Add(i);
            return results;
        }

        var (node, matched) = MatchPatternCore(pattern);
        if (!matched) return results;

        CollectLeaves(node, (int)node.DepthFromRoot, results);
        return results;
    }

    public int CountOccurrences(string pattern)
    {
        ThrowIfDisposed();
        if (pattern == null) throw new ArgumentNullException(nameof(pattern));
        return CountOccurrences(pattern.AsSpan());
    }

    public int CountOccurrences(ReadOnlySpan<char> pattern)
    {
        ThrowIfDisposed();
        if (pattern.IsEmpty) return _textSource.Length;
        var (node, matched) = MatchPatternCore(pattern);
        return matched ? (int)node.LeafCount : 0;
    }

    public string LongestRepeatedSubstring()
    {
        ThrowIfDisposed();
        if (_cachedLrs != null) return _cachedLrs;

        var deepest = FindDeepestInternalNode(new PersistentSuffixTreeNode(_storage, _rootOffset));
        if (deepest.IsNull || deepest.Offset == _rootOffset)
        {
            _cachedLrs = string.Empty;
            return _cachedLrs;
        }

        int length = (int)deepest.DepthFromRoot + LengthOf(deepest);
        var occurrences = new List<int>();
        CollectLeaves(deepest, (int)deepest.DepthFromRoot, occurrences);
        _cachedLrs = occurrences.Count == 0 ? string.Empty : _textSource.Substring(occurrences[0], length);
        return _cachedLrs;
    }

    public IReadOnlyList<string> GetAllSuffixes()
    {
        ThrowIfDisposed();
        var results = new List<string>();
        foreach (var suffix in EnumerateSuffixes())
            results.Add(suffix);
        return results;
    }

    public IEnumerable<string> EnumerateSuffixes()
    {
        ThrowIfDisposed();
        return EnumerateSuffixesCore(new PersistentSuffixTreeNode(_storage, _rootOffset));
    }

    private IEnumerable<string> EnumerateSuffixesCore(PersistentSuffixTreeNode root)
    {
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
                    yield return _textSource.Substring(suffixStart, suffixLen);
                }
                continue;
            }

            long arrayBase = node.ChildrenHead;
            int childCount = node.ChildCount;
            // Push children in reverse order so the stack pops them in ascending (lex) order
            for (int ci = childCount - 1; ci >= 0; ci--)
            {
                var entry = new PersistentChildEntry(_storage, arrayBase + (long)ci * PersistentConstants.CHILD_ENTRY_SIZE);
                stack.Push((new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset), nodeDepth));
            }
        }
    }

    public string LongestCommonSubstring(string other)
    {
        ThrowIfDisposed();
        var (substring, _, _) = LongestCommonSubstringInfo(other);
        return substring;
    }

    public string LongestCommonSubstring(ReadOnlySpan<char> other)
    {
        ThrowIfDisposed();
        var (substring, _, _) = LongestCommonSubstringInfo(new string(other));
        return substring;
    }

    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(string other)
    {
        ThrowIfDisposed();
        var results = FindAllLcsInternal(other, firstOnly: true);
        if (results.PositionsInText.Count == 0)
            return (string.Empty, -1, -1);
        return (results.Substring, results.PositionsInText[0], results.PositionsInOther[0]);
    }

    public (string Substring, int PositionInText, int PositionInOther) LongestCommonSubstringInfo(ReadOnlySpan<char> other)
    {
        ThrowIfDisposed();
        return LongestCommonSubstringInfo(new string(other));
    }

    public (string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther) FindAllLongestCommonSubstrings(string other)
    {
        ThrowIfDisposed();
        var results = FindAllLcsInternal(other, firstOnly: false);
        return (results.Substring, results.PositionsInText, results.PositionsInOther);
    }

    /// <summary>
    /// O(m) LCS using suffix-link-based streaming (same algorithm as FindExactMatchAnchors).
    /// Walks through <paramref name="other"/> maintaining the longest match via suffix links,
    /// instead of restarting from root for each position.
    /// </summary>
    private (string Substring, List<int> PositionsInText, List<int> PositionsInOther) FindAllLcsInternal(string other, bool firstOnly)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (other.Length == 0 || _textSource.Length == 0)
            return (string.Empty, new List<int>(), new List<int>());

        var root = new PersistentSuffixTreeNode(_storage, _rootOffset);
        int maxLen = 0;
        var bestMatches = new List<(PersistentSuffixTreeNode Node, int MatchEndInOther)>();

        var currentNode = root;
        PersistentSuffixTreeNode currentEdge = PersistentSuffixTreeNode.Null(_storage);
        int edgeOffset = 0;
        int currentMatchLen = 0;

        for (int i = 0; i < other.Length; i++)
        {
            uint c = (uint)other[i];

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

                // Rescan from currentNode to restore edge position
                int nodeDepth = GetNodeDepth(currentNode);
                int remaining = currentMatchLen - nodeDepth;

                if (remaining > 0)
                {
                    int pos = i - remaining;
                    currentEdge = PersistentSuffixTreeNode.Null(_storage);
                    edgeOffset = 0;

                    while (remaining > 0)
                    {
                        if (!currentNode.TryGetChild((uint)other[pos], out var nc) || nc.IsNull)
                            break;
                        int edgeLen = LengthOf(nc);
                        if (edgeLen <= remaining)
                        {
                            pos += edgeLen;
                            remaining -= edgeLen;
                            currentNode = nc;
                        }
                        else
                        {
                            currentEdge = nc;
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

            // Track best matches
            if (currentMatchLen > maxLen)
            {
                maxLen = currentMatchLen;
                bestMatches.Clear();
                bestMatches.Add((currentEdge.IsNull ? currentNode : currentEdge, i));
            }
            else if (currentMatchLen == maxLen && maxLen > 0)
            {
                bestMatches.Add((currentEdge.IsNull ? currentNode : currentEdge, i));
            }
        }

        if (maxLen == 0)
            return (string.Empty, new List<int>(), new List<int>());

        var positionsInText = new List<int>();
        var positionsInOther = new List<int>();

        foreach (var match in bestMatches)
        {
            positionsInOther.Add(match.MatchEndInOther - maxLen + 1);

            if (firstOnly)
            {
                // Get one leaf position
                var occurrences = new List<int>();
                CollectLeaves(match.Node, (int)match.Node.DepthFromRoot, occurrences);
                if (occurrences.Count > 0)
                    positionsInText.Add(occurrences[0]);
                break;
            }
            else
            {
                CollectLeaves(match.Node, (int)match.Node.DepthFromRoot, positionsInText);
            }
        }

        string substring = other.Substring(bestMatches[0].MatchEndInOther - maxLen + 1, maxLen);
        return (substring, positionsInText, positionsInOther);
    }

    public string PrintTree()
    {
        ThrowIfDisposed();
        var sb = new StringBuilder(Math.Max(256, _textSource.Length * 100));
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        sb.Append(ci, $"Content length: {_textSource.Length}").AppendLine();
        sb.AppendLine();
        sb.Append(ci, $"0:ROOT").AppendLine();

        var root = new PersistentSuffixTreeNode(_storage, _rootOffset);
        int rootChildCount = root.ChildCount;

        // Frame: (parent node, child count, current index, depth)
        var stack = new Stack<(PersistentSuffixTreeNode Node, int ChildCount, int Index, int Depth)>();
        if (rootChildCount > 0)
            stack.Push((root, rootChildCount, 0, 0));

        while (stack.Count > 0)
        {
            var (node, childCount, index, depth) = stack.Pop();
            if (index >= childCount) continue;

            // Push continuation for next sibling
            stack.Push((node, childCount, index + 1, depth));

            // Read child directly from sorted array (zero allocation)
            long entryOffset = node.ChildrenHead + (long)index * PersistentConstants.CHILD_ENTRY_SIZE;
            long childNodeOffset = _storage.ReadInt64(entryOffset + PersistentConstants.CHILD_OFFSET_NODE);
            var child = new PersistentSuffixTreeNode(_storage, childNodeOffset);

            int childDepth = depth + 1;
            sb.Append(' ', childDepth * 2);
            sb.Append(ci, $"{childDepth}: ");
            AppendLabel(sb, child);

            if (child.IsLeaf)
            {
                sb.AppendLine(" (Leaf)");
            }
            else
            {
                if (child.SuffixLink != PersistentConstants.NULL_OFFSET)
                {
                    var linkNode = new PersistentSuffixTreeNode(_storage, child.SuffixLink);
                    if (linkNode.Offset != _rootOffset)
                    {
                        int firstChar = GetSymbolAt((int)linkNode.Start);
                        if (firstChar >= 0)
                            sb.Append(ci, $" -> [{(char)firstChar}]");
                    }
                }
                sb.AppendLine();

                if (child.ChildCount > 0)
                    stack.Push((child, child.ChildCount, 0, childDepth));
            }
        }

        return sb.ToString();
    }

    private void AppendLabel(StringBuilder sb, PersistentSuffixTreeNode node)
    {
        int len = LengthOf(node);
        for (int i = 0; i < len; i++)
        {
            int s = GetSymbolAt((int)node.Start + i);
            if (s == -1) { sb.Append('#'); break; }
            sb.Append((char)s);
        }
    }

    /// <inheritdoc/>
    public IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> FindExactMatchAnchors(
        string query, int minLength)
    {
        ThrowIfDisposed();
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
            // Walk to any child — prefer non-terminator (sorted: terminator is first if present)
            int childCount = current.ChildCount;
            long arrayBase = current.ChildrenHead;
            PersistentSuffixTreeNode bestChild = PersistentSuffixTreeNode.Null(_storage);
            for (int ci = 0; ci < childCount; ci++)
            {
                var entry = new PersistentChildEntry(_storage, arrayBase + (long)ci * PersistentConstants.CHILD_ENTRY_SIZE);
                bestChild = new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset);
                if (entry.Key != PersistentConstants.TERMINATOR_KEY)
                    break;
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
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(visitor);
        TraverseCore(new PersistentSuffixTreeNode(_storage, _rootOffset), 0, visitor);
    }

    private void TraverseCore(PersistentSuffixTreeNode node, int depth, ISuffixTreeVisitor visitor)
    {
        // Frame: node, child count, current child index, depth (character-based)
        var stack = new Stack<(PersistentSuffixTreeNode Node, int ChildCount, int Index, int Depth)>();

        // Visit root node
        visitor.VisitNode((int)node.Start, (int)node.End, (int)node.LeafCount, node.ChildCount, depth);
        if (!node.IsLeaf)
        {
            stack.Push((node, node.ChildCount, 0, depth));
        }

        while (stack.Count > 0)
        {
            var (current, childCount, index, parentDepth) = stack.Pop();

            if (index >= childCount)
            {
                // All children processed — exit branch (unless this is the root frame)
                if (stack.Count > 0)
                    visitor.ExitBranch();
                continue;
            }

            // Push continuation for next sibling
            stack.Push((current, childCount, index + 1, parentDepth));

            // Read child directly from sorted array (zero allocation)
            long entryOffset = current.ChildrenHead + (long)index * PersistentConstants.CHILD_ENTRY_SIZE;
            uint key = _storage.ReadUInt32(entryOffset + PersistentConstants.CHILD_OFFSET_KEY);
            long childOffset = _storage.ReadInt64(entryOffset + PersistentConstants.CHILD_OFFSET_NODE);
            var child = new PersistentSuffixTreeNode(_storage, childOffset);

            visitor.EnterBranch((int)key);
            visitor.VisitNode((int)child.Start, (int)child.End, (int)child.LeafCount, child.ChildCount, parentDepth);

            int childDepth = parentDepth + LengthOf(child);

            if (!child.IsLeaf)
            {
                stack.Push((child, child.ChildCount, 0, childDepth));
            }
            else
            {
                visitor.ExitBranch();
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
            int childCount = current.ChildCount;
            long arrayBase = current.ChildrenHead;
            for (int ci = 0; ci < childCount; ci++)
            {
                var entry = new PersistentChildEntry(_storage, arrayBase + (long)ci * PersistentConstants.CHILD_ENTRY_SIZE);
                stack.Push((new PersistentSuffixTreeNode(_storage, entry.ChildNodeOffset), childDepth));
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

                // Push non-leaf children in reverse order for left-to-right DFS
                // (no intermediate List — read directly from sorted child array)
                int childCount = node.ChildCount;
                long arrayBase = node.ChildrenHead;
                for (int ci = childCount - 1; ci >= 0; ci--)
                {
                    long entryOffset = arrayBase + (long)ci * PersistentConstants.CHILD_ENTRY_SIZE;
                    long childOffset = _storage.ReadInt64(entryOffset + PersistentConstants.CHILD_OFFSET_NODE);
                    var child = new PersistentSuffixTreeNode(_storage, childOffset);
                    if (!child.IsLeaf)
                        stack.Push(child);
                }
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
