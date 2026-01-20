using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;

namespace SuffixTree.Tests.Robustness
{
    /// <summary>
    /// Tests for thread safety and concurrent access.
    /// </summary>
    [TestFixture]
    public class ThreadSafetyTests
    {
        #region Concurrent Reads

        [Test]
        public void ConcurrentContains_SameTree_ThreadSafe()
        {
            var st = SuffixTree.Build("the quick brown fox jumps over the lazy dog");
            var patterns = new[] { "quick", "brown", "fox", "jumps", "over", "lazy", "dog", "the" };
            var exceptions = new ConcurrentBag<Exception>();

            Parallel.ForEach(Enumerable.Range(0, 1000), i =>
            {
                try
                {
                    var pattern = patterns[i % patterns.Length];
                    var result = st.Contains(pattern);
                    Assert.That(result, Is.True);
                }
                catch (Exception ex)
                {
                    exceptions.Add(ex);
                }
            });

            Assert.That(exceptions, Is.Empty,
                $"Exceptions occurred: {string.Join(", ", exceptions.Select(e => e.Message))}");
        }

        #endregion

        #region Mixed Operations

        [Test]
        public void ConcurrentMixedOperations_ThreadSafe()
        {
            var st = SuffixTree.Build("mississippi");
            var exceptions = new ConcurrentBag<Exception>();

            Parallel.ForEach(Enumerable.Range(0, 200), i =>
            {
                try
                {
                    switch (i % 5)
                    {
                        case 0:
                            _ = st.Contains("issi");
                            break;
                        case 1:
                            _ = st.FindAllOccurrences("i");
                            _ = st.CountOccurrences("ss");
                            break;
                        case 3:
                            _ = st.LongestRepeatedSubstring();
                            break;
                        case 4:
                            _ = st.LongestCommonSubstring("missouri");
                            break;
                    }
                }
                catch (Exception ex)
                {
                    exceptions.Add(ex);
                }
            });

            Assert.That(exceptions, Is.Empty);
        }

        #endregion

        #region Concurrent Construction

        [Test]
        public void ConcurrentBuild_DifferentTrees_ThreadSafe()
        {
            var texts = Enumerable.Range(0, 100)
                .Select(i => $"text_{i}_" + new string((char)('a' + i % 26), 100))
                .ToArray();

            var results = new ConcurrentBag<SuffixTree>();

            Parallel.ForEach(texts, text =>
            {
                var st = SuffixTree.Build(text);
                results.Add(st);
            });

            Assert.That(results, Has.Count.EqualTo(100));

            // Verify each tree is correct
            foreach (var st in results)
            {
                Assert.That(st.Contains(st.Text.Substring(0, 10)), Is.True);
            }
        }

        #endregion

        #region ThreadStatic Buffer Safety

        [Test]
        public void LCS_ConcurrentWithThreadStatic_NoInterference()
        {
            // This tests that ThreadStatic buffer doesn't cause issues
            var st = SuffixTree.Build(new string('x', 500) + "common" + new string('y', 500));

            var tasks = Enumerable.Range(0, 10).Select(i => Task.Run(() =>
            {
                var other = new string((char)('a' + i), 100) + "common" + new string((char)('a' + i), 100);
                var lcs = st.LongestCommonSubstring(other);
                return lcs;
            })).ToArray();

            Task.WaitAll(tasks);

            foreach (var task in tasks)
            {
                Assert.That(task.Result, Is.EqualTo("common"));
            }
        }

        #endregion

        #region Result Consistency

        [Test]
        public void ConcurrentReads_ReturnConsistentResults()
        {
            var st = SuffixTree.Build("abracadabra");
            var expectedPositions = st.FindAllOccurrences("a").OrderBy(x => x).ToList();
            var expectedCount = st.CountOccurrences("a");

            var results = new ConcurrentBag<(List<int> positions, int count)>();

            Parallel.ForEach(Enumerable.Range(0, 100), _ =>
            {
                var positions = st.FindAllOccurrences("a").OrderBy(x => x).ToList();
                var count = st.CountOccurrences("a");
                results.Add((positions, count));
            });

            foreach (var (positions, count) in results)
            {
                Assert.That(positions, Is.EqualTo(expectedPositions));
                Assert.That(count, Is.EqualTo(expectedCount));
            }
        }

        #endregion
    }
}
