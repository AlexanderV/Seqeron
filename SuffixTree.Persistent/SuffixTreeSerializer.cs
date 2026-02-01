using System;
using System.Collections.Generic;
using System.IO;
using System.Security.Cryptography;
using System.Text;

namespace SuffixTree.Persistent
{
    /// <summary>
    /// Provides layout-independent serialization and checksumming for suffix trees.
    /// </summary>
    public static class SuffixTreeSerializer
    {
        private const long LOGICAL_MAGIC = 0x53544C4F47494341L; // "STLOGICA"
        private const int VERSION = 1;

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
                byte[] textBytes = Encoding.Unicode.GetBytes(tree.Text);
                sha256.TransformBlock(textBytes, 0, textBytes.Length, null, 0);
                
                // Hash the tree structure deterministically
                tree.Traverse(hasher);
                
                sha256.TransformFinalBlock(Array.Empty<byte>(), 0, 0);
                return sha256.Hash!;
            }
        }

        /// <summary>
        /// Exports the suffix tree to a stream in a logical, layout-independent format.
        /// </summary>
        public static void Export(ISuffixTree tree, Stream stream)
        {
            ArgumentNullException.ThrowIfNull(tree);
            ArgumentNullException.ThrowIfNull(stream);

            using (var writer = new BinaryWriter(stream, Encoding.Unicode, leaveOpen: true))
            {
                writer.Write(LOGICAL_MAGIC);
                writer.Write(VERSION);
                writer.Write(tree.Text);
                writer.Write(tree.NodeCount);

                var visitor = new ExportVisitor(writer);
                tree.Traverse(visitor);
            }
        }

        /// <summary>
        /// Imports a suffix tree from a logical format stream into the specified storage provider.
        /// </summary>
        public static ISuffixTree Import(Stream stream, IStorageProvider target)
        {
            ArgumentNullException.ThrowIfNull(stream);
            ArgumentNullException.ThrowIfNull(target);

            using (var reader = new BinaryReader(stream, Encoding.Unicode, leaveOpen: true))
            {
                long magic = reader.ReadInt64();
                if (magic != LOGICAL_MAGIC)
                    throw new InvalidDataException("Invalid logical suffix tree format.");
                
                int version = reader.ReadInt32();
                if (version != VERSION)
                    throw new NotSupportedException($"Version {version} is not supported.");

                string text = reader.ReadString();
                int nodeCountCap = reader.ReadInt32();

                // Preliminary allocation for header
                target.Allocate(PersistentConstants.HEADER_SIZE);
                
                var nodeCount = 0;
                long rootOffset = ImportNodeRecursive(reader, target, ref nodeCount);
                
                // Write text
                long textOffset = target.Allocate(text.Length * 2);
                for (int i = 0; i < text.Length; i++)
                    target.WriteChar(textOffset + (i * 2), text[i]);

                // Finalize header
                target.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, PersistentConstants.MAGIC_NUMBER);
                target.WriteInt32(PersistentConstants.HEADER_OFFSET_VERSION, 1);
                target.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, rootOffset);
                target.WriteInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF, textOffset);
                target.WriteInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN, text.Length);
                target.WriteInt32(PersistentConstants.HEADER_OFFSET_NODE_COUNT, nodeCount);
                target.WriteInt64(PersistentConstants.HEADER_OFFSET_SIZE, target.Size);

                return PersistentSuffixTree.Load(target);
            }
        }

        private static long ImportNodeRecursive(BinaryReader reader, IStorageProvider storage, ref int nodeCount)
        {
            int start = reader.ReadInt32();
            int end = reader.ReadInt32();
            int leafCount = reader.ReadInt32();
            int childCount = reader.ReadInt32();
            int depth = reader.ReadInt32();

            nodeCount++;
            long nodeOffset = storage.Allocate(PersistentConstants.NODE_SIZE);
            var node = new PersistentSuffixTreeNode(storage, nodeOffset);
            node.Start = start;
            node.End = end;
            node.LeafCount = leafCount;
            node.DepthFromRoot = depth;
            node.SuffixLink = PersistentConstants.NULL_OFFSET;
            node.ChildrenHead = PersistentConstants.NULL_OFFSET;
            node.ChildCount = 0; // Will be incremented by SetChild

            for (int i = 0; i < childCount; i++)
            {
                int key = reader.ReadInt32();
                long childOffset = ImportNodeRecursive(reader, storage, ref nodeCount);
                node.SetChild(key, new PersistentSuffixTreeNode(storage, childOffset));
            }

            return nodeOffset;
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

        private class ExportVisitor : ISuffixTreeVisitor
        {
            private readonly BinaryWriter _writer;

            public ExportVisitor(BinaryWriter writer) => _writer = writer;

            public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
            {
                _writer.Write(startIndex);
                _writer.Write(endIndex);
                _writer.Write(leafCount);
                _writer.Write(childCount);
                _writer.Write(depth);
            }

            public void EnterBranch(int key) => _writer.Write(key);
            public void ExitBranch() { }
        }
    }
}
