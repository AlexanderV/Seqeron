using System;
using System.Collections.Generic;
using System.IO;
using System.Security.Cryptography;
using System.Text;

namespace SuffixTree.Persistent
{
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

                // Hash the text first to bind the tree to a specific string
                byte[] textBytes = Encoding.Unicode.GetBytes(tree.Text.ToString() ?? string.Empty);
                sha256.TransformBlock(textBytes, 0, textBytes.Length, null, 0);

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
                writer.Write(tree.Text.ToString() ?? string.Empty);
                writer.Write(tree.NodeCount);
                writer.Write(hash.Length);
                writer.Write(hash);
            }
        }

        /// <summary>
        /// Imports a suffix tree from a stream into the specified storage provider.
        /// Rebuilds the tree from the stored text using Ukkonen's algorithm,
        /// guaranteeing full functionality including suffix links.
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

                string text = reader.ReadString();
                int expectedNodeCount = reader.ReadInt32();
                int hashLen = reader.ReadInt32();
                byte[] expectedHash = reader.ReadBytes(hashLen);

                // Rebuild the tree from text — Ukkonen's creates suffix links natively
                var builder = new PersistentSuffixTreeBuilder(target);
                long rootOffset = builder.Build(text);
                var tree = new PersistentSuffixTree(target, rootOffset, new StringTextSource(text));

                // Validate structural integrity
                if (tree.NodeCount != expectedNodeCount)
                    throw new InvalidDataException(
                        $"Node count mismatch after rebuild: expected {expectedNodeCount}, got {tree.NodeCount}.");

                var actualHash = CalculateLogicalHash(tree);
                if (!CryptographicOperations.FixedTimeEquals(actualHash, expectedHash))
                    throw new InvalidDataException("Structural hash mismatch after rebuild.");

                return tree;
            }
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

            string text = tree.Text.ToString() ?? string.Empty;
            return PersistentSuffixTreeFactory.Create(text, filePath);
        }

        /// <summary>
        /// Loads a suffix tree from a memory-mapped file previously created by
        /// <see cref="SaveToFile"/> or <see cref="PersistentSuffixTreeFactory.Create(string, string?)"/>.
        /// </summary>
        /// <param name="filePath">Path to the existing tree file.</param>
        /// <returns>A read-only <see cref="ISuffixTree"/> backed by the MMF (caller must dispose).</returns>
        public static ISuffixTree LoadFromFile(string filePath)
        {
            return PersistentSuffixTreeFactory.Load(filePath);
        }

        private class HashVisitor : ISuffixTreeVisitor
        {
            private readonly SHA256 _sha;
            private readonly byte[] _buffer = new byte[4];

            public HashVisitor(SHA256 sha) => _sha = sha;

            private void HashInt(int value)
            {
                BitConverter.TryWriteBytes(_buffer, value);
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
}
