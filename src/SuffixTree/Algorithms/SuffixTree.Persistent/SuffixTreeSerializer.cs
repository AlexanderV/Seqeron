using System.Buffers;
using System.Buffers.Binary;
using System.Security.Cryptography;
using System.Text;

namespace SuffixTree.Persistent;

/// <summary>
/// Provides serialization, deserialization, and checksumming for suffix trees.
/// <para>
/// <b>Format v2</b>: stores text + structural hash. Import rebuilds the tree via
/// <see cref="PersistentSuffixTreeBuilder"/> (Ukkonen's algorithm), guaranteeing
/// 100% functionality including suffix links for <c>FindExactMatchAnchors</c>.
/// </para>
/// <para>
/// For direct memory-mapped file persistence, use
/// <see cref="SaveToFile"/> / <see cref="LoadFromFile"/>.
/// </para>
/// </summary>
public static class SuffixTreeSerializer
{
    private const long LOGICAL_MAGIC = 0x53544C4F47494332L; // "STLOGIC2"
    private const int VERSION = 2;

    /// <summary>
    /// Calculates a logical SHA256 hash of the suffix tree.
    /// The hash is identical for trees with the same content regardless of their memory layout.
    /// </summary>
    public static byte[] CalculateLogicalHash(ISuffixTree tree)
    {
        ArgumentNullException.ThrowIfNull(tree);

        using (var sha256 = SHA256.Create())
        {
            var hasher = new HashVisitor(sha256);

            // Hash the text as raw UTF-16LE code units (same bytes as the export payload)
            WriteTextCodeUnits(tree.Text, (buf, count) => sha256.TransformBlock(buf, 0, count, null, 0));

            // Hash the tree structure deterministically
            tree.Traverse(hasher);

            sha256.TransformFinalBlock(Array.Empty<byte>(), 0, 0);
            return sha256.Hash!;
        }
    }

    /// <summary>
    /// Exports the suffix tree to a stream. Stores text and a structural hash
    /// for validation on import.
    /// </summary>
    public static void Export(ISuffixTree tree, Stream stream)
    {
        ArgumentNullException.ThrowIfNull(tree);
        ArgumentNullException.ThrowIfNull(stream);

        var hash = CalculateLogicalHash(tree);

        using (var writer = new BinaryWriter(stream, Encoding.Unicode, leaveOpen: true))
        {
            writer.Write(LOGICAL_MAGIC);
            writer.Write(VERSION);

            // Write text as raw UTF-16LE code units in chunks (no full materialization for
            // large MMF sources). Encoding.Unicode is NOT used: it would replace lone
            // surrogates — and surrogate pairs split by a chunk boundary — with U+FFFD.
            var text = tree.Text;
            writer.Write7BitEncodedInt(text.Length);
            WriteTextCodeUnits(text, (buf, count) => writer.Write(buf, 0, count));

            writer.Write(tree.NodeCount);
            writer.Write(hash.Length);
            writer.Write(hash);
        }
    }

    /// <summary>
    /// Imports a suffix tree from a stream into the specified storage provider.
    /// Rebuilds the tree from the stored text using Ukkonen's algorithm,
    /// guaranteeing full functionality including suffix links.
    /// <para>
    /// Construction starts with Compact (32-bit) layout. If the tree exceeds the
    /// uint32 address space, the builder transparently transitions to Large (64-bit)
    /// mid-build using the v6 hybrid layout jump table.
    /// </para>
    /// </summary>
    public static ISuffixTree Import(Stream stream, IStorageProvider target)
    {
        ArgumentNullException.ThrowIfNull(stream);
        ArgumentNullException.ThrowIfNull(target);

        using (var reader = new BinaryReader(stream, Encoding.Unicode, leaveOpen: true))
        {
            long magic = reader.ReadInt64();
            if (magic != LOGICAL_MAGIC)
                throw new InvalidDataException("Invalid suffix tree format (magic mismatch). " +
                    "This may be a v1 file — only v2 format is supported.");

            int version = reader.ReadInt32();
            if (version != VERSION)
                throw new NotSupportedException($"Format version {version} is not supported (expected {VERSION}).");

            // Read the text payload (Write7BitEncodedInt length + raw UTF-16LE code units)
            int textLen;
            try { textLen = reader.Read7BitEncodedInt(); }
            catch (FormatException ex) { throw new InvalidDataException("Corrupted stream: invalid text length prefix.", ex); }
            if (textLen < 0 || textLen > int.MaxValue / sizeof(char))
                throw new InvalidDataException($"Corrupted stream: invalid text length {textLen}.");

            byte[] textBytes = reader.ReadBytes(textLen * sizeof(char));
            if (textBytes.Length != textLen * sizeof(char))
                throw new InvalidDataException(
                    $"Truncated stream: expected {textLen} characters, got {textBytes.Length / sizeof(char)}.");
            string text = Utf16CodeUnits.Read(textBytes);

            int expectedNodeCount;
            byte[] expectedHash;
            try
            {
                expectedNodeCount = reader.ReadInt32();
                int hashLen = reader.ReadInt32();
                if (hashLen != SHA256.HashSizeInBytes)
                    throw new InvalidDataException(
                        $"Corrupted stream: hash length {hashLen} (expected {SHA256.HashSizeInBytes}).");
                expectedHash = reader.ReadBytes(hashLen);
            }
            catch (EndOfStreamException ex)
            {
                throw new InvalidDataException("Truncated stream: node count / hash trailer is missing.", ex);
            }
            if (expectedHash.Length != SHA256.HashSizeInBytes)
                throw new InvalidDataException(
                    $"Truncated stream: expected {SHA256.HashSizeInBytes} hash bytes, got {expectedHash.Length}.");

            var textSource = new StringTextSource(text);
            return ImportBuild(target, textSource, expectedNodeCount, expectedHash);
        }
    }

    private static ISuffixTree ImportBuild(
        IStorageProvider storage, StringTextSource textSource,
        int expectedNodeCount, byte[] expectedHash)
    {
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        long rootOffset = builder.Build(textSource);

        // Base layout is always Compact for v6 hybrid trees; the large zone is
        // resolved dynamically via TransitionOffset.
        NodeLayout layout = NodeLayout.Compact;
        var tree = new PersistentSuffixTree(storage, rootOffset, textSource, layout,
            builder.TransitionOffset, builder.JumpTableStart, builder.JumpTableEnd,
            builder.DeepestInternalNodeOffset);

        // Validate structural integrity
        if (tree.NodeCount != expectedNodeCount)
            throw new InvalidDataException(
                $"Node count mismatch after rebuild: expected {expectedNodeCount}, got {tree.NodeCount}.");

        var actualHash = CalculateLogicalHash(tree);
        if (!CryptographicOperations.FixedTimeEquals(actualHash, expectedHash))
            throw new InvalidDataException("Structural hash mismatch after rebuild.");

        return tree;
    }

    /// <summary>
    /// Saves a suffix tree to a memory-mapped file. Rebuilds from the tree's text
    /// using Ukkonen's algorithm, producing a native persistent format with full
    /// functionality including suffix links.
    /// </summary>
    /// <param name="tree">The source tree (any <see cref="ISuffixTree"/> implementation).</param>
    /// <param name="filePath">File path for the memory-mapped file.</param>
    /// <returns>A new <see cref="ISuffixTree"/> backed by the MMF (caller must dispose).</returns>
    public static ISuffixTree SaveToFile(ISuffixTree tree, string filePath)
    {
        ArgumentNullException.ThrowIfNull(tree);
        if (string.IsNullOrEmpty(filePath))
            throw new ArgumentException("File path must be provided.", nameof(filePath));

        return PersistentSuffixTreeFactory.Create(tree.Text, filePath);
    }

    /// <summary>
    /// Loads a suffix tree from a memory-mapped file previously created by
    /// <see cref="SaveToFile"/> or <see cref="PersistentSuffixTreeFactory.Create(ITextSource, string?, IProgress{ValueTuple{string, double}}?)"/>.
    /// </summary>
    /// <param name="filePath">Path to the existing tree file.</param>
    /// <returns>A read-only <see cref="ISuffixTree"/> backed by the MMF (caller must dispose).</returns>
    public static ISuffixTree LoadFromFile(string filePath)
    {
        return PersistentSuffixTreeFactory.Load(filePath);
    }

    /// <summary>
    /// Streams <paramref name="text"/> as raw little-endian UTF-16 code units in fixed-size
    /// chunks. Shared by <see cref="CalculateLogicalHash"/> and <see cref="Export"/> so the
    /// hashed bytes and the exported payload are identical and lossless.
    /// </summary>
    private static void WriteTextCodeUnits(ITextSource text, Action<byte[], int> sink)
    {
        const int chunkSize = 4096;
        byte[] byteBuf = ArrayPool<byte>.Shared.Rent(chunkSize * sizeof(char));
        try
        {
            int textLen = text.Length;
            for (int offset = 0; offset < textLen; offset += chunkSize)
            {
                int count = Math.Min(chunkSize, textLen - offset);
                int byteCount = Utf16CodeUnits.Write(text.Slice(offset, count), byteBuf);
                sink(byteBuf, byteCount);
            }
        }
        finally
        {
            ArrayPool<byte>.Shared.Return(byteBuf);
        }
    }

    private sealed class HashVisitor : ISuffixTreeVisitor
    {
        private readonly SHA256 _sha;
        private readonly byte[] _buffer = new byte[4];

        public HashVisitor(SHA256 sha) => _sha = sha;

        private void HashInt(int value)
        {
            BinaryPrimitives.WriteInt32LittleEndian(_buffer, value);
            _sha.TransformBlock(_buffer, 0, 4, null, 0);
        }

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            HashInt(startIndex);
            HashInt(endIndex);
            HashInt(leafCount);
            HashInt(childCount);
        }

        public void EnterBranch(int key) => HashInt(key);
        public void ExitBranch() => HashInt(-999); // Structure sentinel
    }
}
