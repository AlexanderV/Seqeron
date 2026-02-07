using System;
using System.Buffers;
using System.Collections.Generic;
using System.Text;
using System.Threading;

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
    private int _disposed;
    private volatile string? _cachedLrs;

    // Pre-computed during build; NULL_OFFSET means not available (loaded trees)
    private readonly long _deepestInternalNodeOffset;

    // Single source of truth for layout + hybrid zone info
    private readonly HybridLayout _hybrid;

    public PersistentSuffixTree(IStorageProvider storage, long rootOffset,
        ITextSource? textSource = null, NodeLayout? layout = null,
        long transitionOffset = -1, long jumpTableStart = -1, long jumpTableEnd = -1,
        long deepestInternalNodeOffset = -1)
        : this(storage, rootOffset, textSource, textSource == null, layout, transitionOffset, jumpTableStart, jumpTableEnd, deepestInternalNodeOffset)
    {
    }

    /// <summary>
    /// Internal constructor allowing explicit control over text source ownership.
    /// </summary>
    internal PersistentSuffixTree(IStorageProvider storage, long rootOffset,
        ITextSource? textSource, bool ownsTextSource, NodeLayout? layout,
        long transitionOffset = -1, long jumpTableStart = -1, long jumpTableEnd = -1,
        long deepestInternalNodeOffset = -1)
    {
        _storage = storage ?? throw new ArgumentNullException(nameof(storage));
        _rootOffset = rootOffset;
        _deepestInternalNodeOffset = deepestInternalNodeOffset;
        _hybrid = new HybridLayout(_storage, layout ?? NodeLayout.Compact, transitionOffset, jumpTableStart, jumpTableEnd);

        if (textSource != null)
        {
            _textSource = textSource;
            _ownsTextSource = ownsTextSource;
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
        var layout = NodeLayout.ForVersion(version);

        long storageSize = storage.Size;
        int headerSize = version == 5 ? PersistentConstants.HEADER_SIZE_V5 : PersistentConstants.HEADER_SIZE;

        // C16: Validate recorded SIZE against actual storage to detect truncated files
        long recordedSize = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_SIZE);
        if (recordedSize > 0 && recordedSize != storageSize)
            throw new InvalidOperationException(
                $"Invalid storage format: header size {recordedSize} does not match actual storage size {storageSize}. File may be truncated or corrupted.");

        long root = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_ROOT);
        if (root < headerSize || root >= storageSize)
            throw new InvalidOperationException(
                $"Invalid storage format: root offset {root} is outside valid range [{headerSize}, {storageSize}).");

        long transitionOffset = -1;
        long jumpTableStart = -1;
        long jumpTableEnd = -1;
        if (version == 5)
        {
            transitionOffset = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_TRANSITION);
            jumpTableStart = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_JUMP_START);
            jumpTableEnd = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_JUMP_END);

            if (transitionOffset < headerSize || transitionOffset > storageSize)
                throw new InvalidOperationException(
                    $"Invalid storage format: transition offset {transitionOffset} is outside valid range [{headerSize}, {storageSize}].");

            if (jumpTableEnd < jumpTableStart)
                throw new InvalidOperationException(
                    $"Invalid storage format: jump table end ({jumpTableEnd}) < start ({jumpTableStart}).");

            if (jumpTableStart < headerSize || jumpTableEnd > storageSize)
                throw new InvalidOperationException(
                    $"Invalid storage format: jump table [{jumpTableStart}, {jumpTableEnd}) is outside valid storage range.");
        }

        // Validate TEXT_OFF and TEXT_LEN
        long textOff = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF);
        int textLen = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN);

        if (textLen < 0)
            throw new InvalidOperationException(
                $"Invalid storage format: text length {textLen} is negative.");

        if (textOff < headerSize || textOff >= storageSize)
            throw new InvalidOperationException(
                $"Invalid storage format: text offset {textOff} is outside valid range [{headerSize}, {storageSize}).");

        long textEnd = textOff + (long)textLen * sizeof(char);
        if (textEnd > storageSize)
            throw new InvalidOperationException(
                $"Invalid storage format: text region [{textOff}, {textEnd}) exceeds storage size {storageSize}.");

        return new PersistentSuffixTree(storage, root, layout: layout,
            transitionOffset: transitionOffset, jumpTableStart: jumpTableStart, jumpTableEnd: jumpTableEnd);
    }

    public ITextSource Text { get { ThrowIfDisposed(); return _textSource; } }

    /// <summary>Whether this tree uses the hybrid v5 format with dual zones.</summary>
    internal bool IsHybrid => _hybrid.IsHybrid;

    /// <summary>Transition offset (compact/large boundary), or -1 for single-format.</summary>
    internal long TransitionOffset => _hybrid.TransitionOffset;

    /// <summary>Start of contiguous jump table, or -1 for single-format.</summary>
    internal long JumpTableStart => _hybrid.JumpTableStart;

    /// <summary>End of jump table, or -1 for single-format.</summary>
    internal long JumpTableEnd => _hybrid.JumpTableEnd;

    /// <summary>
    /// Returns the correct <see cref="NodeLayout"/> for a node at the given offset.
    /// Delegates to <see cref="HybridLayout.LayoutForOffset"/>.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    internal NodeLayout LayoutForOffset(long offset) => _hybrid.LayoutForOffset(offset);

    /// <summary>
    /// Resolves an offset that might be a jump-table entry.
    /// Delegates to <see cref="HybridLayout.ResolveJump"/>.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    internal long ResolveJump(long offset) => _hybrid.ResolveJump(offset);

    /// <summary>
    /// Creates a <see cref="PersistentSuffixTreeNode"/> with the correct layout for its zone.
    /// Delegates to <see cref="HybridLayout.NodeAt"/>.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    internal PersistentSuffixTreeNode NodeAt(long offset) => _hybrid.NodeAt(offset);

    /// <summary>
    /// Reads child information from a parent node, handling hybrid jump entries.
    /// Delegates to <see cref="HybridLayout.ReadChildArrayInfo"/>.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    internal (long ArrayBase, NodeLayout EntryLayout, int Count) ReadChildArrayInfo(PersistentSuffixTreeNode parent)
        => _hybrid.ReadChildArrayInfo(parent);

    /// <summary>
    /// Resolves a suffix link offset, dereferencing through the jump table if needed.
    /// </summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    internal long ResolveSuffixLink(PersistentSuffixTreeNode node)
    {
        long raw = node.SuffixLink;
        if (raw == PersistentConstants.NULL_OFFSET) return raw;
        return ResolveJump(raw);
    }

    public int NodeCount { get { ThrowIfDisposed(); return _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT); } }

    public int LeafCount
    {
        get
        {
            ThrowIfDisposed();
            uint rawCount = NodeAt(_rootOffset).LeafCount;
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

        var deepest = FindDeepestInternalNode(NodeAt(_rootOffset));
        if (deepest.IsNull || deepest.Offset == _rootOffset)
        {
            _cachedLrs = string.Empty;
            return string.Empty;
        }

        int length = (int)deepest.DepthFromRoot + LengthOf(deepest);
        var occurrences = new List<int>();
        CollectLeaves(deepest, (int)deepest.DepthFromRoot, occurrences);
        string result = occurrences.Count == 0 ? string.Empty : _textSource.Substring(occurrences[0], length);
        _cachedLrs = result;
        return result;
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
        return EnumerateSuffixesCore(NodeAt(_rootOffset));
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

            var (arrayBase, entryLayout, childCount) = ReadChildArrayInfo(node);
            // Push children in reverse order so the stack pops them in ascending (lex) order
            for (int ci = childCount - 1; ci >= 0; ci--)
            {
                long entryOff = arrayBase + (long)ci * entryLayout.ChildEntrySize;
                long childNodeOffset = entryLayout.ReadOffset(_storage, entryOff + NodeLayout.ChildOffsetNode);
                stack.Push((NodeAt(childNodeOffset), nodeDepth));
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
    /// O(m) LCS using suffix-link-based streaming — delegates to shared <see cref="SuffixTreeAlgorithms"/>.
    /// </summary>
    private (string Substring, List<int> PositionsInText, List<int> PositionsInOther) FindAllLcsInternal(string other, bool firstOnly)
    {
        var nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindAllLcs<PersistentSuffixTreeNode, PersistentSuffixTreeNavigator>(ref nav, other, firstOnly);
    }

    public string PrintTree()
    {
        ThrowIfDisposed();
        var sb = new StringBuilder(Math.Max(256, Math.Min(_textSource.Length, 100_000) * 100));
        var ci = System.Globalization.CultureInfo.InvariantCulture;
        sb.Append(ci, $"Content length: {_textSource.Length}").AppendLine();
        sb.AppendLine();
        sb.Append(ci, $"0:ROOT").AppendLine();

        var root = NodeAt(_rootOffset);
        int rootChildCount;
        long rootArrayBase;
        NodeLayout rootEntryLayout;
        {
            var (ab, el, cc) = ReadChildArrayInfo(root);
            rootChildCount = cc;
            rootArrayBase = ab;
            rootEntryLayout = el;
        }

        // Frame: (parent arrayBase, entryLayout, child count, current index, depth)
        var stack = new Stack<(long ArrayBase, NodeLayout EntryLayout, int ChildCount, int Index, int Depth)>();
        if (rootChildCount > 0)
            stack.Push((rootArrayBase, rootEntryLayout, rootChildCount, 0, 0));

        while (stack.Count > 0)
        {
            var (arrBase, entryLay, childCount, index, depth) = stack.Pop();
            if (index >= childCount) continue;

            // Push continuation for next sibling
            stack.Push((arrBase, entryLay, childCount, index + 1, depth));

            // Read child directly from sorted array
            long entryOffset = arrBase + (long)index * entryLay.ChildEntrySize;
            long childNodeOffset = entryLay.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
            var child = NodeAt(childNodeOffset);

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
                long suffLink = ResolveSuffixLink(child);
                if (suffLink != PersistentConstants.NULL_OFFSET)
                {
                    var linkNode = NodeAt(suffLink);
                    if (linkNode.Offset != _rootOffset)
                    {
                        int firstChar = GetSymbolAt((int)linkNode.Start);
                        if (firstChar >= 0)
                            sb.Append(ci, $" -> [{(char)firstChar}]");
                    }
                }
                sb.AppendLine();

                var (cBase, cLayout, cCount) = ReadChildArrayInfo(child);
                if (cCount > 0)
                    stack.Push((cBase, cLayout, cCount, 0, childDepth));
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
        var nav = CreateNavigator();
        return SuffixTreeAlgorithms.FindExactMatchAnchors<PersistentSuffixTreeNode, PersistentSuffixTreeNavigator>(ref nav, query, minLength);
    }

    /// <inheritdoc/>
    public void Traverse(ISuffixTreeVisitor visitor)
    {
        ThrowIfDisposed();
        ArgumentNullException.ThrowIfNull(visitor);
        TraverseCore(NodeAt(_rootOffset), 0, visitor);
    }

    private void TraverseCore(PersistentSuffixTreeNode node, int depth, ISuffixTreeVisitor visitor)
    {
        // Frame: arrayBase, entryLayout, child count, current child index, depth (character-based)
        var stack = new Stack<(long ArrayBase, NodeLayout EntryLayout, int ChildCount, int Index, int Depth)>();

        // Visit root node
        var (rootAB, rootEL, rootCC) = ReadChildArrayInfo(node);
        visitor.VisitNode((int)node.Start, (int)node.End, (int)node.LeafCount, rootCC, depth);
        if (!node.IsLeaf)
        {
            stack.Push((rootAB, rootEL, rootCC, 0, depth));
        }

        while (stack.Count > 0)
        {
            var (arrBase, entryLay, childCount, index, parentDepth) = stack.Pop();

            if (index >= childCount)
            {
                // All children processed — exit branch (unless this is the root frame)
                if (stack.Count > 0)
                    visitor.ExitBranch();
                continue;
            }

            // Push continuation for next sibling
            stack.Push((arrBase, entryLay, childCount, index + 1, parentDepth));

            // Read child directly from sorted array
            long entryOffset = arrBase + (long)index * entryLay.ChildEntrySize;
            uint key = _storage.ReadUInt32(entryOffset + NodeLayout.ChildOffsetKey);
            long childOffset = entryLay.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
            var child = NodeAt(childOffset);

            var (cAB, cEL, cCC) = ReadChildArrayInfo(child);
            visitor.EnterBranch((int)key);
            visitor.VisitNode((int)child.Start, (int)child.End, (int)child.LeafCount, cCC, parentDepth);

            int childDepth = parentDepth + LengthOf(child);

            if (!child.IsLeaf)
            {
                stack.Push((cAB, cEL, cCC, 0, childDepth));
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
        var node = NodeAt(_rootOffset);
        int i = 0;
        while (i < pattern.Length)
        {
            if (!TryGetChildOf(node, (uint)pattern[i], out var child) || child.IsNull)
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

    /// <summary>
    /// Zone-aware TryGetChild: reads the child array with the correct entry layout.
    /// </summary>
    internal bool TryGetChildOf(PersistentSuffixTreeNode parent, uint key, out PersistentSuffixTreeNode child)
    {
        var (arrayBase, entryLayout, count) = ReadChildArrayInfo(parent);
        if (count == 0) { child = PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout); return false; }

        int lo = 0, hi = count - 1;
        int signedKey = (int)key;

        while (lo <= hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            long entryOffset = arrayBase + (long)mid * entryLayout.ChildEntrySize;
            int midKey = (int)_storage.ReadUInt32(entryOffset + NodeLayout.ChildOffsetKey);

            if (midKey == signedKey)
            {
                long childOffset = entryLayout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                child = NodeAt(childOffset);
                return true;
            }
            if (midKey < signedKey) lo = mid + 1;
            else hi = mid - 1;
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

    private PersistentSuffixTreeNode FindDeepestInternalNode(PersistentSuffixTreeNode root)
    {
        // Use pre-computed offset from builder if available (O(1) path)
        if (_deepestInternalNodeOffset != PersistentConstants.NULL_OFFSET
            && _deepestInternalNodeOffset != _rootOffset)
        {
            return NodeAt(_deepestInternalNodeOffset);
        }

        if (root.IsLeaf) return PersistentSuffixTreeNode.Null(_storage, _hybrid.Layout);

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

                var (arrayBase, entryLayout, childCount) = ReadChildArrayInfo(node);
                for (int ci = childCount - 1; ci >= 0; ci--)
                {
                    long entryOffset = arrayBase + (long)ci * entryLayout.ChildEntrySize;
                    long childOffset = entryLayout.ReadOffset(_storage, entryOffset + NodeLayout.ChildOffsetNode);
                    var child = NodeAt(childOffset);
                    if (!child.IsLeaf)
                        stack.Push(child);
                }
            }
        }

        return deepest;
    }

    public void Dispose()
    {
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;

        // Dispose textSource BEFORE storage: MemoryMappedTextSource must
        // ReleasePointer() while the underlying accessor is still alive.
        // Use try/finally to guarantee storage is disposed even if textSource throws.
        try
        {
            if (_ownsTextSource && _textSource is IDisposable disposableText)
                disposableText.Dispose();
        }
        finally
        {
            _storage.Dispose();
        }
    }

    private void ThrowIfDisposed()
    {
        if (Volatile.Read(ref _disposed) != 0)
            throw new ObjectDisposedException(nameof(PersistentSuffixTree));
    }
}
