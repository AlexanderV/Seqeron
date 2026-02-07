using System;
using System.IO;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests
{
    /// <summary>
    /// Tests for post-Dispose behavior, read-only mode, TrimToSize,
    /// PrintTree, and FindAllLongestCommonSubstrings.
    /// </summary>
    [TestFixture]
    public class DisposeBehaviorTests
    {
        // ─── Post-Dispose: PersistentSuffixTree ───────────────────────

        [Test]
        public void Dispose_HeapTree_AllMethodsThrowObjectDisposedException()
        {
            var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(
                new StringTextSource("banana"));
            tree.Dispose();

            Assert.Multiple(() =>
            {
                Assert.Throws<ObjectDisposedException>(() => _ = tree.Text);
                Assert.Throws<ObjectDisposedException>(() => _ = tree.NodeCount);
                Assert.Throws<ObjectDisposedException>(() => _ = tree.LeafCount);
                Assert.Throws<ObjectDisposedException>(() => _ = tree.MaxDepth);
                Assert.Throws<ObjectDisposedException>(() => _ = tree.IsEmpty);
                Assert.Throws<ObjectDisposedException>(() => tree.Contains("a"));
                Assert.Throws<ObjectDisposedException>(() => tree.FindAllOccurrences("a"));
                Assert.Throws<ObjectDisposedException>(() => tree.CountOccurrences("a"));
                Assert.Throws<ObjectDisposedException>(() => tree.LongestRepeatedSubstring());
                Assert.Throws<ObjectDisposedException>(() => tree.GetAllSuffixes());
                Assert.Throws<ObjectDisposedException>(() => tree.EnumerateSuffixes());
                Assert.Throws<ObjectDisposedException>(() => tree.LongestCommonSubstring("x"));
                Assert.Throws<ObjectDisposedException>(() => tree.LongestCommonSubstringInfo("x"));
                Assert.Throws<ObjectDisposedException>(() => tree.FindAllLongestCommonSubstrings("x"));
                Assert.Throws<ObjectDisposedException>(() => tree.PrintTree());
                Assert.Throws<ObjectDisposedException>(() => tree.FindExactMatchAnchors("a", 1));
                Assert.Throws<ObjectDisposedException>(() => tree.Traverse(new NullVisitor()));
            });
        }

        [Test]
        public void Dispose_DoubleDispose_DoesNotThrow()
        {
            var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(
                new StringTextSource("test"));
            tree.Dispose();
            Assert.DoesNotThrow(() => tree.Dispose());
        }

        // ─── MappedFileStorageProvider: Read-only write rejection ─────

        [Test]
        public void ReadOnlyProvider_WriteThrowsInvalidOperation()
        {
            string tempFile = Path.Combine(Path.GetTempPath(), $"ST_RO_{Guid.NewGuid():N}.tree");
            try
            {
                // Create a valid tree file first
                using (PersistentSuffixTreeFactory.Create(new StringTextSource("hello"), tempFile) as IDisposable) { }

                // Open read-only
                using var ro = new MappedFileStorageProvider(tempFile, readOnly: true);
                Assert.Multiple(() =>
                {
                    Assert.Throws<InvalidOperationException>(() => ro.WriteInt32(0, 42));
                    Assert.Throws<InvalidOperationException>(() => ro.WriteUInt32(0, 42));
                    Assert.Throws<InvalidOperationException>(() => ro.WriteInt64(0, 42));
                    Assert.Throws<InvalidOperationException>(() => ro.WriteChar(0, 'x'));
                    Assert.Throws<InvalidOperationException>(() => ro.WriteBytes(0, new byte[4], 0, 4));
                    Assert.Throws<InvalidOperationException>(() => ro.EnsureCapacity(long.MaxValue));
                });
            }
            finally
            {
                try { File.Delete(tempFile); } catch { }
            }
        }

        // ─── TrimToSize ──────────────────────────────────────────────

        [Test]
        public void TrimToSize_ReducesFileSize()
        {
            string tempFile = Path.Combine(Path.GetTempPath(), $"ST_Trim_{Guid.NewGuid():N}.tree");
            try
            {
                using (var tree = PersistentSuffixTreeFactory.Create(new StringTextSource("abc"), tempFile) as IDisposable) { }

                var fi = new FileInfo(tempFile);
                // After TrimToSize (called by factory), file should be ≤ initial capacity
                Assert.That(fi.Length, Is.LessThan(65536),
                    "File should be trimmed below initial capacity");
                Assert.That(fi.Length, Is.GreaterThan(0), "File should not be empty");
            }
            finally
            {
                try { File.Delete(tempFile); } catch { }
            }
        }

        // ─── PrintTree ──────────────────────────────────────────────

        [Test]
        public void PrintTree_ContainsRootAndLeaves()
        {
            using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(
                new StringTextSource("ab"));
            string output = tree.PrintTree();

            Assert.Multiple(() =>
            {
                Assert.That(output, Does.Contain("ROOT"));
                Assert.That(output, Does.Contain("(Leaf)"));
                Assert.That(output, Does.Contain("Content length: 2"));
            });
        }

        [Test]
        public void PrintTree_EmptyText_ShowsRootOnly()
        {
            using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(
                new StringTextSource(""));
            string output = tree.PrintTree();

            Assert.That(output, Does.Contain("ROOT"));
            Assert.That(output, Does.Contain("Content length: 0"));
        }

        // ─── FindAllLongestCommonSubstrings ─────────────────────────

        [Test]
        public void FindAllLCS_ReturnsAllPositions()
        {
            using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(
                new StringTextSource("abcabc"));

            var (substring, positionsInText, positionsInOther) =
                tree.FindAllLongestCommonSubstrings("xabcy");

            Assert.Multiple(() =>
            {
                Assert.That(substring, Is.EqualTo("abc"));
                Assert.That(positionsInText, Has.Count.GreaterThanOrEqualTo(1));
                Assert.That(positionsInOther, Has.Count.GreaterThanOrEqualTo(1));
                Assert.That(positionsInOther, Does.Contain(1)); // "abc" at index 1 in "xabcy"
            });
        }

        [Test]
        public void FindAllLCS_NoMatch_ReturnsEmpty()
        {
            using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(
                new StringTextSource("abc"));

            var (substring, positionsInText, positionsInOther) =
                tree.FindAllLongestCommonSubstrings("xyz");

            Assert.That(substring, Is.EqualTo(string.Empty));
            Assert.That(positionsInText, Is.Empty);
            Assert.That(positionsInOther, Is.Empty);
        }

        // ─── S10: Dispose exception safety ────────────────────────

        [Test]
        public void Dispose_TextSourceThrows_StorageStillDisposed()
        {
            // Build a tree with a tracking storage so we can verify storage.Dispose was called
            var storage = new TrackingStorageProvider();
            var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
            long root = builder.Build(new StringTextSource("abc"));

            // Use internal constructor to inject a ThrowingTextSource with ownsTextSource=true
            var tree = new PersistentSuffixTree(storage, root,
                textSource: new ThrowingTextSource("abc"), ownsTextSource: true,
                layout: NodeLayout.Compact);

            // Currently Dispose does NOT use try/finally, so if textSource throws,
            // storage.Dispose() is skipped. After the fix, storage must still be disposed.
            Exception? caughtException = null;
            try { tree.Dispose(); }
            catch (Exception ex) { caughtException = ex; }

            // The throwing text source's exception should propagate (caller should know),
            // but storage must be disposed regardless
            Assert.That(caughtException, Is.Not.Null, "Exception from text source should propagate");
            Assert.That(storage.IsDisposed, Is.True,
                "Storage must be disposed even when textSource.Dispose() throws");
        }

        [Test]
        public void Dispose_StorageThrows_CompletesWithoutCorruption()
        {
            // Even if storage.Dispose() throws, the _disposed flag should be set
            var storage = new ThrowingStorageProvider();
            var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
            long root = builder.Build(new StringTextSource("x"));

            var tree = new PersistentSuffixTree(storage, root, new StringTextSource("x"), NodeLayout.Compact);

            // First Dispose throws from storage
            Assert.Throws<InvalidOperationException>(() => tree.Dispose());

            // Second Dispose must not throw (_disposed is already true)
            Assert.DoesNotThrow(() => tree.Dispose());
        }

        // ─── HeapStorageProvider: Dispose ────────────────────────────

        [Test]
        public void HeapStorageProvider_DoubleDispose_DoesNotThrow()
        {
            var heap = new HeapStorageProvider();
            heap.Dispose();
            Assert.DoesNotThrow(() => heap.Dispose());
        }

        // ─── Helper ─────────────────────────────────────────────────

        private class NullVisitor : ISuffixTreeVisitor
        {
            public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) { }
            public void EnterBranch(int key) { }
            public void ExitBranch() { }
        }

        /// <summary>ITextSource that throws on Dispose for S10 testing.</summary>
        private sealed class ThrowingTextSource : ITextSource, IDisposable
        {
            private readonly string _text;
            public ThrowingTextSource(string text) => _text = text;
            public int Length => _text.Length;
            public char this[int index] => _text[index];
            public string Substring(int start, int length) => _text.Substring(start, length);
            public ReadOnlySpan<char> Slice(int start, int length) => _text.AsSpan(start, length);
            public IEnumerator<char> GetEnumerator() => _text.GetEnumerator();
            System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
            public override string ToString() => _text;
            public void Dispose() => throw new InvalidOperationException("ThrowingTextSource.Dispose boom");
        }

        /// <summary>Storage decorator that tracks Dispose calls for S10 testing.</summary>
        private sealed class TrackingStorageProvider : IStorageProvider
        {
            private readonly HeapStorageProvider _inner = new();
            public bool IsDisposed { get; private set; }
            public long Size => _inner.Size;
            public void EnsureCapacity(long capacity) => _inner.EnsureCapacity(capacity);
            public long Allocate(int size) => _inner.Allocate(size);
            public int ReadInt32(long offset) => _inner.ReadInt32(offset);
            public void WriteInt32(long offset, int value) => _inner.WriteInt32(offset, value);
            public uint ReadUInt32(long offset) => _inner.ReadUInt32(offset);
            public void WriteUInt32(long offset, uint value) => _inner.WriteUInt32(offset, value);
            public long ReadInt64(long offset) => _inner.ReadInt64(offset);
            public void WriteInt64(long offset, long value) => _inner.WriteInt64(offset, value);
            public char ReadChar(long offset) => _inner.ReadChar(offset);
            public void WriteChar(long offset, char value) => _inner.WriteChar(offset, value);
            public void ReadBytes(long offset, byte[] buffer, int bufferOffset, int count) => _inner.ReadBytes(offset, buffer, bufferOffset, count);
            public void WriteBytes(long offset, byte[] data, int dataOffset, int count) => _inner.WriteBytes(offset, data, dataOffset, count);
            public void Dispose() { IsDisposed = true; _inner.Dispose(); }
        }

        /// <summary>Storage decorator that throws on Dispose for S10 testing.</summary>
        private sealed class ThrowingStorageProvider : IStorageProvider
        {
            private readonly HeapStorageProvider _inner = new();
            public long Size => _inner.Size;
            public void EnsureCapacity(long capacity) => _inner.EnsureCapacity(capacity);
            public long Allocate(int size) => _inner.Allocate(size);
            public int ReadInt32(long offset) => _inner.ReadInt32(offset);
            public void WriteInt32(long offset, int value) => _inner.WriteInt32(offset, value);
            public uint ReadUInt32(long offset) => _inner.ReadUInt32(offset);
            public void WriteUInt32(long offset, uint value) => _inner.WriteUInt32(offset, value);
            public long ReadInt64(long offset) => _inner.ReadInt64(offset);
            public void WriteInt64(long offset, long value) => _inner.WriteInt64(offset, value);
            public char ReadChar(long offset) => _inner.ReadChar(offset);
            public void WriteChar(long offset, char value) => _inner.WriteChar(offset, value);
            public void ReadBytes(long offset, byte[] buffer, int bufferOffset, int count) => _inner.ReadBytes(offset, buffer, bufferOffset, count);
            public void WriteBytes(long offset, byte[] data, int dataOffset, int count) => _inner.WriteBytes(offset, data, dataOffset, count);
            public void Dispose() => throw new InvalidOperationException("ThrowingStorageProvider.Dispose boom");
        }
    }
}
