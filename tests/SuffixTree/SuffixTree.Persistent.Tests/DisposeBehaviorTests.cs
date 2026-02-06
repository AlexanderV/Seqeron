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
    }
}
