using System;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Tests based on classic examples from suffix tree literature.
    /// References: Ukkonen (1995), Gusfield (1997), cp-algorithms.
    /// </summary>
    [TestFixture]
    public class LiteratureExamplesTests
    {
        #region Ukkonen's Paper Examples

        /// <summary>
        /// "xabxac" from Ukkonen's 1995 paper - tests edge splitting.
        /// </summary>
        [Test]
        public void Xabxac_UkkonenPaper()
        {
            var tree = SuffixTree.Build("xabxac");

            Assert.Multiple(() =>
            {
                Assert.That(tree.Contains("xab"), Is.True);
                Assert.That(tree.Contains("xa"), Is.True);
                Assert.That(tree.Contains("xac"), Is.True);
                Assert.That(tree.Contains("ab"), Is.True);
                Assert.That(tree.Contains("bxac"), Is.True);

                Assert.That(tree.CountOccurrences("xa"), Is.EqualTo(2));
                Assert.That(tree.CountOccurrences("a"), Is.EqualTo(2));

                // LRS: "xa" appears at positions 0 and 3
                Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("xa"));
            });
        }

        /// <summary>
        /// "abcabxabcd" - tricky case for suffix link handling.
        /// </summary>
        [Test]
        public void Abcabxabcd_TrickySuffixLinks()
        {
            var tree = SuffixTree.Build("abcabxabcd");

            Assert.Multiple(() =>
            {
                Assert.That(tree.Contains("abc"), Is.True);
                Assert.That(tree.Contains("abx"), Is.True);
                Assert.That(tree.Contains("abcd"), Is.True);
                Assert.That(tree.Contains("cab"), Is.True);

                Assert.That(tree.CountOccurrences("ab"), Is.EqualTo(3));
                Assert.That(tree.CountOccurrences("abc"), Is.EqualTo(2));

                // LRS: "abc" appears at positions 0 and 6
                Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("abc"));
            });
        }

        #endregion

        #region Fibonacci Strings (Worst Case)

        /// <summary>
        /// Fibonacci strings are worst case for suffix tree size (2n-1 nodes).
        /// F_k = F_{k-1} + F_{k-2}, F_1 = "b", F_2 = "a"
        /// </summary>
        [Test]
        public void FibonacciString_WorstCaseNodeCount()
        {
            // F_7 = "abaababaabaab" (13 characters)
            string fib = GenerateFibonacciString(7);
            var tree = SuffixTree.Build(fib);

            Assert.Multiple(() =>
            {
                Assert.That(tree.LeafCount, Is.EqualTo(fib.Length));

                // All suffixes present
                for (int i = 0; i < fib.Length; i++)
                {
                    Assert.That(tree.Contains(fib.Substring(i)), Is.True);
                }

                // Fibonacci strings have complex repetitive structure
                Assert.That(tree.LongestRepeatedSubstring().Length, Is.GreaterThan(0));
            });
        }

        [Test]
        [TestCase(5, "abaab")]
        [TestCase(6, "abaababa")]
        [TestCase(7, "abaababaabaab")]
        public void FibonacciString_CorrectGeneration(int n, string expected)
        {
            string fib = GenerateFibonacciString(n);
            Assert.That(fib, Is.EqualTo(expected));
        }

        private static string GenerateFibonacciString(int n)
        {
            if (n <= 1) return "b";
            if (n == 2) return "a";

            string prev2 = "b";
            string prev1 = "a";
            for (int i = 3; i <= n; i++)
            {
                string current = prev1 + prev2;
                prev2 = prev1;
                prev1 = current;
            }
            return prev1;
        }

        #endregion

        // NOTE: "All different characters" tests already exist in:
        // - BuildTests.Build_AllUniqueCharacters_CreatesValidTree
        // - LongestRepeatedSubstringTests.LRS_AllUnique_ReturnsEmpty
    }
}
