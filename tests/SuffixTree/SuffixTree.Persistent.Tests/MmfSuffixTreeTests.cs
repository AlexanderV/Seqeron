using System;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class MmfSuffixTreeTests : SuffixTreeTestBase
    {
        private string _tempFile = string.Empty;

        [SetUp]
        public void Setup()
        {
            _tempFile = Path.Combine(Path.GetTempPath(), $"ST_{Guid.NewGuid():N}.tree");
        }

        [TearDown]
        public void Cleanup()
        {
            if (File.Exists(_tempFile))
            {
                try { File.Delete(_tempFile); } catch { /* Ignore */ }
            }
        }

        protected override ISuffixTree CreateTree(string text)
        {
            return PersistentSuffixTreeFactory.Create(text, _tempFile);
        }

        [Test]
        public void Persistence_TreeIsRetrievableAfterReload()
        {
            string text = "abracadabra";
            
            // 1. Build and save
            using (var originalTree = PersistentSuffixTreeFactory.Create(text, _tempFile) as IDisposable)
            {
                // Usage during active mapping
            }

            // 2. Load from file
            using (var loadedTreeDisposable = PersistentSuffixTreeFactory.Load(_tempFile) as IDisposable)
            {
                var loadedTree = (ISuffixTree)loadedTreeDisposable!;
                Assert.Multiple(() =>
                {
                    Assert.That(loadedTree.Text.ToString(), Is.EqualTo(text));
                    Assert.That(loadedTree.Contains("abra"), Is.True);
                    Assert.That(loadedTree.CountOccurrences("a"), Is.EqualTo(5));
                    
                    var results = loadedTree.FindAllOccurrences("ra").ToList();
                    results.Sort();
                    Assert.That(results, Is.EqualTo(new List<int> { 2, 9 }));
                });
            }
        }

        [Test]
        public void LargeTree_PersistenceWorks()
        {
            // Create a larger string (approx 50KB) to test capacity growth and mapping
            string text = string.Concat(Enumerable.Repeat("abcde12345", 5000));
            
            using (var tree = PersistentSuffixTreeFactory.Create(text, _tempFile) as IDisposable)
            {
                Assert.That(tree, Is.Not.Null);
            }

            using (var loadedDisposable = PersistentSuffixTreeFactory.Load(_tempFile) as IDisposable)
            {
                var loaded = (ISuffixTree)loadedDisposable!;
                Assert.That(loaded.Contains("abcde12345abcde12345"), Is.True);
                Assert.That(loaded.CountOccurrences("abcde12345"), Is.EqualTo(5000));
            }
        }

        [Test]
        public void StressTest_VeryLargeData_Works()
        {
            // 1MB string to force multiple re-allocations and cross-page mapping
            int length = 1024 * 1024;
            var sb = new StringBuilder(length);
            for (int i = 0; i < length / 10; i++) sb.Append("ATGC123456");
            string text = sb.ToString();

            using (var tree = PersistentSuffixTreeFactory.Create(text, _tempFile) as IDisposable)
            {
                var st = (ISuffixTree)tree!;
                Assert.That(st.Contains("ATGC123456ATGC"), Is.True);
                Assert.That(st.CountOccurrences("ATGC"), Is.EqualTo(length / 10));
            }
        }

        [Test]
        public void Concurrency_MultipleReaders_Works()
        {
            string text = "concurrent_read_test_data";
            (PersistentSuffixTreeFactory.Create(text, _tempFile) as IDisposable)?.Dispose();

            // Open multiple instances simultaneously (READ ONLY SHARE)
            using (var tree1 = PersistentSuffixTreeFactory.Load(_tempFile) as IDisposable)
            using (var tree2 = PersistentSuffixTreeFactory.Load(_tempFile) as IDisposable)
            {
                var st1 = (ISuffixTree)tree1!;
                var st2 = (ISuffixTree)tree2!;

                Assert.Multiple(() =>
                {
                    Assert.That(st1.Contains("concurrent"), Is.True);
                    Assert.That(st2.Contains("read_test"), Is.True);
                });
            }
        }

        [Test]
        public void Load_WithInvalidFile_ThrowsException()
        {
            File.WriteAllText(_tempFile, "Not a suffix tree file");
            Assert.Throws<InvalidOperationException>(() => PersistentSuffixTreeFactory.Load(_tempFile));
        }
    }
}
