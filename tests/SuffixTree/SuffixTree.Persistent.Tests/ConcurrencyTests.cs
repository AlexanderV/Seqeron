using System;
using System.Collections.Concurrent;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using NUnit.Framework;
using SuffixTree;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class ConcurrencyTests
    {
        private string _tempFile = string.Empty;

        [SetUp]
        public void Setup()
        {
            _tempFile = Path.Combine(Path.GetTempPath(), $"ST_Concurrency_{Guid.NewGuid():N}.tree");
        }

        [TearDown]
        public void Cleanup()
        {
            if (File.Exists(_tempFile))
            {
                try { File.Delete(_tempFile); } catch { }
            }
        }

        [Test]
        public void MmfSuffixTree_ParallelReads_ShouldBeThreadSafe()
        {
            // Setup a decent size tree
            string text = string.Concat(Enumerable.Repeat("abcde12345", 1000)); // 10,000 chars
            using (var builder = PersistentSuffixTreeFactory.Create(new StringTextSource(text), _tempFile) as IDisposable)
            {
                // Ensure it's fully built
            }

            using (var treeDisposable = PersistentSuffixTreeFactory.Load(_tempFile) as IDisposable)
            {
                var tree = (ISuffixTree)treeDisposable!;
                var exceptions = new ConcurrentQueue<Exception>();

                // Hit it with parallel requests
                Parallel.For(0, 1000, new ParallelOptions { MaxDegreeOfParallelism = Environment.ProcessorCount }, i =>
                {
                    try
                    {
                        // 1. Thread-safe method calls
                        Assert.That(tree.Contains("abcde12345"), Is.True);
                        Assert.That(tree.CountOccurrences("123"), Is.EqualTo(1000));
                        Assert.That(tree.Contains("xyz"), Is.False);

                        // 2. Thread-safe property access
                        var len = tree.Text.Length;
                        Assert.That(len, Is.EqualTo(10000));

                        // 3. Thread-safe slice access (via ITextSource)
                        // Random access to mimic real world usage
                        int start = i * 10 % (len - 10);
                        var slice = tree.Text.Slice(start, 5);
                        Assert.That(slice.Length, Is.EqualTo(5));
                    }
                    catch (Exception ex)
                    {
                        exceptions.Enqueue(ex);
                    }
                });

                if (!exceptions.IsEmpty)
                {
                    Assert.Fail($"Encountered {exceptions.Count} exceptions during parallel execution. First: {exceptions.First()}");
                }
            }
        }

        [Test]
        public void InMemorySuffixTree_ParallelReads_ShouldBeThreadSafe()
        {
            // Also test the standard in-memory tree just in case
            string text = string.Concat(Enumerable.Repeat("ATGTGC", 1000));
            var tree = SuffixTree.Build(text);
            var exceptions = new ConcurrentQueue<Exception>();

            Parallel.For(0, 5000, new ParallelOptions { MaxDegreeOfParallelism = Environment.ProcessorCount }, i =>
            {
                try
                {
                    Assert.That(tree.Contains("ATGT"), Is.True);
                    Assert.That(tree.CountOccurrences("GC"), Is.EqualTo(1000));

                    // Stress test LRS (which does traversal)
                    // Note: LRS might cache results, but it should still be thread-safe
                    // Actually LRS does a fresh traversal every time in current implementation unless cached?
                    // Checking implementation: deep traversal.
                    // The _cachedMaxDepth etc are pre-calculated.
                    // But LongestRepeatedSubstring() does logic.
                    // Wait, current implementation of LRS uses `DeepestInternalNode` which is pre-calculated.
                    // So it's O(1) mostly accessing properties.
                    // Let's test a method that does traversal: FindAllOccurrences

                    var occurrences = tree.FindAllOccurrences("ATGT").ToList();
                    Assert.That(occurrences.Count, Is.EqualTo(1000));
                }
                catch (Exception ex)
                {
                    exceptions.Enqueue(ex);
                }
            });

            Assert.That(exceptions, Is.Empty, "Exceptions occurred during parallel execution");
        }
    }
}
