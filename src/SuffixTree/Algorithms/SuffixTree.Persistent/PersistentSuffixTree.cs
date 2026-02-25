using System.Buffers;
using System.Text;

namespace SuffixTree.Persistent;

/// <summary>
/// A persistent implementation of ISuffixTree that operates on a storage provider.
/// </summary>
public sealed partial class PersistentSuffixTree : ISuffixTree, IDisposable
{
    private readonly IStorageProvider _storage;
    private readonly long _rootOffset;
    private readonly ITextSource _textSource;
    private readonly bool _ownsTextSource;
    private int _disposed;
    private volatile string? _cachedLrs;

    // Pre-computed during build; NULL_OFFSET means not available (loaded trees)
    private readonly long _deepestInternalNodeOffset;

    // Pre-computed LRS total depth (v6 Slim): DepthFromRoot + EdgeLen of deepest internal node.
    // -1 means not available (v3/v4/v5 trees compute from node.DepthFromRoot).
    private readonly int _lrsDepth;

    // Single source of truth for layout + hybrid zone info
    private readonly HybridLayout _hybrid;

    /// <summary>Initializes a persistent suffix tree from pre-built storage data.</summary>
    public PersistentSuffixTree(IStorageProvider storage, long rootOffset,
        ITextSource? textSource = null, NodeLayout? layout = null,
        long transitionOffset = -1, long jumpTableStart = -1, long jumpTableEnd = -1,
        long deepestInternalNodeOffset = -1, int lrsDepth = -1)
        : this(storage, rootOffset, textSource, textSource == null, layout, transitionOffset, jumpTableStart, jumpTableEnd, deepestInternalNodeOffset, lrsDepth)
    {
    }

    /// <summary>
    /// Internal constructor allowing explicit control over text source ownership.
    /// </summary>
    internal PersistentSuffixTree(IStorageProvider storage, long rootOffset,
        ITextSource? textSource, bool ownsTextSource, NodeLayout? layout,
        long transitionOffset = -1, long jumpTableStart = -1, long jumpTableEnd = -1,
        long deepestInternalNodeOffset = -1, int lrsDepth = -1)
    {
        _storage = storage ?? throw new ArgumentNullException(nameof(storage));
        _rootOffset = rootOffset;
        _deepestInternalNodeOffset = deepestInternalNodeOffset;
        _lrsDepth = lrsDepth;
        _hybrid = new HybridLayout(_storage, layout ?? NodeLayout.Compact, transitionOffset, jumpTableStart, jumpTableEnd);

        if (textSource != null)
        {
            _textSource = textSource;
            _ownsTextSource = ownsTextSource;
        }
        else
        {
            // Load text from storage — check flags for encoding
            long textOff = _storage.ReadInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF);
            int textLen = _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN);
            int flags = _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_FLAGS);
            bool isAscii = (flags & PersistentConstants.FLAG_TEXT_ASCII) != 0;

            if (_storage is MappedFileStorageProvider mappedProvider)
            {
                _textSource = isAscii
                    ? new AsciiMemoryMappedTextSource(mappedProvider.Accessor, textOff, textLen)
                    : new MemoryMappedTextSource(mappedProvider.Accessor, textOff, textLen);
            }
            else
            {
                _textSource = new StringTextSource(LoadStringInternal(textOff, textLen, isAscii));
            }
            _ownsTextSource = true;
        }
    }

    private string LoadStringInternal(long textOff, int textLen, bool isAscii)
    {
        int bytesPerChar = isAscii ? 1 : 2;
        long byteLen = (long)textLen * bytesPerChar;
        if (byteLen > int.MaxValue)
            throw new InvalidOperationException(
                $"Text length {textLen} exceeds maximum loadable size.");
        int byteLenInt = (int)byteLen;
        byte[] bytes = ArrayPool<byte>.Shared.Rent(byteLenInt);
        try
        {
            _storage.ReadBytes(textOff, bytes, 0, byteLenInt);
            if (isAscii)
                return Encoding.ASCII.GetString(bytes, 0, byteLenInt);
            return Encoding.Unicode.GetString(bytes, 0, byteLenInt);
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(bytes);
        }
    }

    /// <summary>Loads a persistent suffix tree from storage, validating the header and format.</summary>
    public static PersistentSuffixTree Load(IStorageProvider storage)
    {
        ArgumentNullException.ThrowIfNull(storage);

        long magic;
        try
        {
            magic = storage.ReadInt64(PersistentConstants.HEADER_OFFSET_MAGIC);
        }
        catch (ArgumentOutOfRangeException ex)
        {
            throw new InvalidOperationException("Invalid storage format: header is truncated.", ex);
        }

        if (magic != PersistentConstants.MAGIC_NUMBER)
            throw new InvalidOperationException("Invalid storage format: Magic number mismatch.");

        int version;
        try
        {
            version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        }
        catch (ArgumentOutOfRangeException ex)
        {
            throw new InvalidOperationException("Invalid storage format: header is truncated.", ex);
        }

        _ = NodeLayout.ForVersion(version);

        var header = PersistentStorageHeaderReader.Read(storage);
        var metadata = PersistentStorageHeaderValidator.Validate(header, storage.Size);

        return new PersistentSuffixTree(storage, metadata.RootOffset, layout: metadata.Layout,
            transitionOffset: metadata.TransitionOffset,
            jumpTableStart: metadata.JumpTableStart,
            jumpTableEnd: metadata.JumpTableEnd,
            deepestInternalNodeOffset: metadata.DeepestNodeOffset,
            lrsDepth: metadata.LrsDepth);
    }

    /// <inheritdoc />
    public ITextSource Text { get { ThrowIfDisposed(); return _textSource; } }

    /// <summary>Whether this tree uses the hybrid v6 format with dual zones.</summary>
    internal bool IsHybrid => _hybrid.IsHybrid;

    /// <summary>Transition offset (compact/large boundary), or -1 for single-format.</summary>
    internal long TransitionOffset => _hybrid.TransitionOffset;

    /// <summary>Start of contiguous jump table, or -1 for single-format.</summary>
    internal long JumpTableStart => _hybrid.JumpTableStart;

    /// <summary>End of jump table, or -1 for single-format.</summary>
    internal long JumpTableEnd => _hybrid.JumpTableEnd;

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

    /// <inheritdoc />
    public int NodeCount { get { ThrowIfDisposed(); return _storage.ReadInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT); } }

    /// <inheritdoc />
    public int LeafCount
    {
        get
        {
            ThrowIfDisposed();
            uint rawCount = NodeAt(_rootOffset).LeafCount;
            return rawCount > 0 ? (int)(rawCount - 1) : 0;
        }
    }

    /// <inheritdoc />
    public int MaxDepth { get { ThrowIfDisposed(); return _textSource.Length; } }

    /// <inheritdoc />
    public bool IsEmpty { get { ThrowIfDisposed(); return _textSource.Length == 0; } }
}
