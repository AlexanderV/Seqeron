using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Property-based tests using brute force verification.
    /// These tests verify correctness by comparing against naive O(n²) algorithms.
    /// </summary>
    [TestFixture]
    public class BruteForceVerificationTests
    {
        #region LCS Brute Force Verification

        /// <summary>
        /// Verify LCS against brute force O(n²m) algorithm.
        /// </summary>
        [Test]
        [TestCase("abcdefgh", "xxcdefxx", 4)]
        [TestCase("banana", "bandana", 3)]
        [TestCase("abc", "def", 0)]
        [TestCase("abcabc", "xabcyabcz", 3)]
        [TestCase("hello", "helloworld", 5)]
        [TestCase("xyz", "abc", 0)]
        public void LCS_MatchesBruteForce(string text, string other, int expectedLength)
        {
            var tree = SuffixTree.Build(text);
            string lcs = tree.LongestCommonSubstring(other);

            // Verify length matches expected
            Assert.That(lcs.Length, Is.EqualTo(expectedLength));

            // Verify LCS exists in both strings
            if (lcs.Length > 0)
            {
                Assert.That(text.Contains(lcs), Is.True, $"LCS '{lcs}' not in text '{text}'");
                Assert.That(other.Contains(lcs), Is.True, $"LCS '{lcs}' not in other '{other}'");
            }

            // Verify matches brute force
            string bruteForceLcs = BruteForceLCS(text, other);
            Assert.That(lcs.Length, Is.EqualTo(bruteForceLcs.Length));
        }

        [Test]
        [Repeat(5)]
        public void LCS_RandomInput_MatchesBruteForce()
        {
            var random = new Random();
            string text = GenerateRandomString(random, random.Next(20, 50), "abc");
            string other = GenerateRandomString(random, random.Next(20, 50), "abc");

            var tree = SuffixTree.Build(text);
            string lcs = tree.LongestCommonSubstring(other);
            string expected = BruteForceLCS(text, other);

            Assert.That(lcs.Length, Is.EqualTo(expected.Length),
                $"Text='{text}', Other='{other}', Got='{lcs}', Expected='{expected}'");
        }

        private static string BruteForceLCS(string s1, string s2)
        {
            string longest = "";
            for (int i = 0; i < s1.Length; i++)
            {
                for (int len = 1; len <= s1.Length - i; len++)
                {
                    string sub = s1.Substring(i, len);
                    if (s2.Contains(sub) && sub.Length > longest.Length)
                    {
                        longest = sub;
                    }
                }
            }
            return longest;
        }

        #endregion

        #region LRS Brute Force Verification

        /// <summary>
        /// Verify LRS against brute force O(n³) algorithm.
        /// NOTE: "banana" and "mississippi" already tested in LongestRepeatedSubstringTests.
        /// These cases test brute force equivalence for cases NOT already covered.
        /// </summary>
        [Test]
        [TestCase("abcdefghij", 0)]  // All different - no LRS
        [TestCase("xyxyxyxy", 6)]    // Highly repetitive
        [TestCase("abacaba", 3)]     // Palindrome-like
        public void LRS_MatchesBruteForce(string text, int expectedLength)
        {
            var tree = SuffixTree.Build(text);
            string lrs = tree.LongestRepeatedSubstring();

            Assert.That(lrs.Length, Is.EqualTo(expectedLength));

            // Verify matches brute force
            string bruteForceLrs = BruteForceLRS(text);
            Assert.That(lrs.Length, Is.EqualTo(bruteForceLrs.Length),
                $"Text='{text}', Got='{lrs}', BruteForce='{bruteForceLrs}'");
        }

        [Test]
        [Repeat(5)]
        public void LRS_RandomInput_MatchesBruteForce()
        {
            var random = new Random();
            string text = GenerateRandomString(random, random.Next(20, 80), "abcd");

            var tree = SuffixTree.Build(text);
            string lrs = tree.LongestRepeatedSubstring();
            string expected = BruteForceLRS(text);

            Assert.That(lrs.Length, Is.EqualTo(expected.Length),
                $"Text='{text}', Got='{lrs}', Expected='{expected}'");
        }

        private static string BruteForceLRS(string text)
        {
            string longest = "";
            for (int i = 0; i < text.Length; i++)
            {
                for (int j = i + 1; j < text.Length; j++)
                {
                    // Find length of common prefix at positions i and j
                    int len = 0;
                    while (j + len < text.Length && text[i + len] == text[j + len])
                    {
                        len++;
                    }
                    if (len > longest.Length)
                    {
                        longest = text.Substring(i, len);
                    }
                }
            }
            return longest;
        }

        #endregion

        // NOTE: FindAllOccurrences brute force tests removed - already exist in
        // FindAllOccurrencesTests.FindAll_AllSubstrings_MatchLinearSearch

        #region Helpers

        private static string GenerateRandomString(Random random, int length, string alphabet)
        {
            return new string(Enumerable.Range(0, length)
                .Select(_ => alphabet[random.Next(alphabet.Length)])
                .ToArray());
        }

        #endregion
    }
}
