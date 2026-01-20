using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Tests for GetAllSuffixes and EnumerateSuffixes.
    /// </summary>
    [TestFixture]
    public class SuffixEnumerationTests
    {
        #region Empty and Single

        [Test]
        public void GetAllSuffixes_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");

            var result = st.GetAllSuffixes().ToList();

            Assert.That(result, Is.Empty);
        }

        [Test]
        public void GetAllSuffixes_SingleChar_ReturnsSingleSuffix()
        {
            var st = SuffixTree.Build("a");

            var result = st.GetAllSuffixes().ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result, Does.Contain("a"));
            });
        }

        [Test]
        public void EnumerateSuffixes_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");

            var result = st.EnumerateSuffixes().ToList();

            Assert.That(result, Is.Empty);
        }

        [Test]
        public void EnumerateSuffixes_SingleChar_ReturnsSingleSuffix()
        {
            var st = SuffixTree.Build("a");

            var result = st.EnumerateSuffixes().ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result, Does.Contain("a"));
            });
        }

        #endregion

        #region All Suffixes Present

        [Test]
        public void GetAllSuffixes_ReturnsAllSuffixes()
        {
            var text = "abc";
            var st = SuffixTree.Build(text);

            var result = st.GetAllSuffixes().OrderBy(s => s).ToList();
            var expected = new[] { "abc", "bc", "c" }.OrderBy(s => s).ToList();

            Assert.That(result, Is.EqualTo(expected));
        }

        [Test]
        public void GetAllSuffixes_Banana_ReturnsAllSuffixes()
        {
            var text = "banana";
            var st = SuffixTree.Build(text);

            var result = st.GetAllSuffixes().OrderBy(s => s).ToList();
            var expected = new[] { "banana", "anana", "nana", "ana", "na", "a" }
                .OrderBy(s => s).ToList();

            Assert.That(result, Is.EqualTo(expected));
        }

        [Test]
        public void GetAllSuffixes_Count_EqualsTextLength()
        {
            var testCases = new[] { "a", "ab", "abc", "banana", "mississippi" };

            foreach (var text in testCases)
            {
                var st = SuffixTree.Build(text);
                var count = st.GetAllSuffixes().Count;

                Assert.That(count, Is.EqualTo(text.Length),
                    $"Suffix count for '{text}' should be {text.Length}");
            }
        }

        #endregion

        #region Uniqueness

        [Test]
        public void GetAllSuffixes_AllUnique()
        {
            var text = "banana";
            var st = SuffixTree.Build(text);

            var result = st.GetAllSuffixes().ToList();
            var uniqueCount = result.Distinct().Count();

            Assert.That(uniqueCount, Is.EqualTo(result.Count),
                "All suffixes should be unique");
        }

        #endregion

        #region Verification Against Manual Calculation

        [Test]
        public void GetAllSuffixes_MatchManualCalculation()
        {
            var testCases = new[] { "hello", "test", "abracadabra" };

            foreach (var text in testCases)
            {
                var st = SuffixTree.Build(text);
                var treeResult = st.GetAllSuffixes().OrderBy(s => s).ToList();

                var expected = Enumerable.Range(0, text.Length)
                    .Select(i => text.Substring(i))
                    .OrderBy(s => s)
                    .ToList();

                Assert.That(treeResult, Is.EqualTo(expected),
                    $"Suffixes mismatch for '{text}'");
            }
        }

        #endregion

        #region EnumerateSuffixes vs GetAllSuffixes Consistency

        [Test]
        public void EnumerateSuffixes_MatchesGetAllSuffixes()
        {
            var testCases = new[] { "a", "abc", "banana", "mississippi" };

            foreach (var text in testCases)
            {
                var st = SuffixTree.Build(text);
                var getAllResult = st.GetAllSuffixes().OrderBy(s => s).ToList();
                var enumerateResult = st.EnumerateSuffixes().OrderBy(s => s).ToList();

                Assert.That(enumerateResult, Is.EqualTo(getAllResult),
                    $"EnumerateSuffixes should match GetAllSuffixes for '{text}'");
            }
        }

        [Test]
        public void EnumerateSuffixes_IsLazyEvaluated()
        {
            // Build a tree with many suffixes
            var st = SuffixTree.Build(new string('a', 1000));

            // Take only first 5 - should not enumerate all
            var first5 = st.EnumerateSuffixes().Take(5).ToList();

            Assert.That(first5, Has.Count.EqualTo(5));
        }

        #endregion

        #region Ordering

        [Test]
        public void GetAllSuffixes_OrderDoesNotMatter_ForCorrectness()
        {
            var text = "abcd";
            var st = SuffixTree.Build(text);

            var result = st.GetAllSuffixes().ToHashSet();
            var expected = new HashSet<string> { "abcd", "bcd", "cd", "d" };

            Assert.That(result.SetEquals(expected), Is.True);
        }

        #endregion
    }
}
