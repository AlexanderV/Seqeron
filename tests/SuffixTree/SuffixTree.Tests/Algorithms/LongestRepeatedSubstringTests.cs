using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Tests for LongestRepeatedSubstring (LRS) functionality.
    /// </summary>
    [TestFixture]
    public class LongestRepeatedSubstringTests
    {
        #region Input Validation

        [Test]
        public void LRS_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo(string.Empty));
        }

        #endregion

        #region No Repetitions

        [Test]
        public void LRS_SingleCharacter_ReturnsEmpty()
        {
            var st = SuffixTree.Build("a");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo(string.Empty));
        }

        [Test]
        public void LRS_AllUnique_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abcdef");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo(string.Empty));
        }

        #endregion

        #region Basic Cases

        [Test]
        public void LRS_SimpleRepetition_FindsIt()
        {
            var st = SuffixTree.Build("abcabc");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("abc"));
        }

        [Test]
        public void LRS_SingleCharRepeated_FindsChar()
        {
            var st = SuffixTree.Build("aa");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("a"));
        }

        [Test]
        public void LRS_Banana_FindsAna()
        {
            var st = SuffixTree.Build("banana");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("ana"));
        }

        [Test]
        public void LRS_Mississippi_FindsIssi()
        {
            var st = SuffixTree.Build("mississippi");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("issi"));
        }

        #endregion

        #region Overlapping Patterns

        [Test]
        public void LRS_OverlappingRepetition_FindsLongest()
        {
            var st = SuffixTree.Build("aaaaa");

            var result = st.LongestRepeatedSubstring();

            // "aaaa" repeats at positions 0 and 1 (overlapping)
            Assert.That(result, Is.EqualTo("aaaa"));
        }

        [Test]
        public void LRS_ababab_FindsAbab()
        {
            var st = SuffixTree.Build("ababab");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("abab"));
        }

        #endregion

        #region Multiple Equal-Length LRS

        [Test]
        public void LRS_MultipleEqualLength_ReturnsOne()
        {
            // "ab" repeats twice, "xy" repeats twice, both length 2
            var st = SuffixTree.Build("abxyabxy");

            var result = st.LongestRepeatedSubstring();

            // Should return the longest, which is "abxy" (length 4)
            Assert.That(result, Is.EqualTo("abxy"));
        }

        [Test]
        public void LRS_DisjointRepeats_ReturnsLongest()
        {
            // "ab" and "cd" both repeat, but not as one string
            var st = SuffixTree.Build("abcdabcd");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("abcd"));
        }

        #endregion

        #region All Characters Same

        [Test]
        public void LRS_AllSameCharacter_FindsNMinus1()
        {
            var text = new string('x', 10);
            var st = SuffixTree.Build(text);

            var result = st.LongestRepeatedSubstring();

            // In "xxxxxxxxxx", "xxxxxxxxx" (9 chars) repeats at positions 0 and 1
            Assert.That(result, Is.EqualTo(new string('x', 9)));
        }

        #endregion

        #region Result Validation

        [Test]
        public void LRS_ResultActuallyRepeats()
        {
            var testCases = new[] { "banana", "mississippi", "abcabc", "ababab", "aabbaabb" };

            foreach (var text in testCases)
            {
                var st = SuffixTree.Build(text);
                var lrs = st.LongestRepeatedSubstring();

                if (lrs.Length > 0)
                {
                    var count = st.CountOccurrences(lrs);
                    Assert.That(count, Is.GreaterThanOrEqualTo(2),
                        $"LRS '{lrs}' in '{text}' should occur at least twice");
                }
            }
        }

        [Test]
        public void LRS_NoLongerRepeatingSubstring()
        {
            var testCases = new[] { "banana", "mississippi", "abcabc", "ababab" };

            foreach (var text in testCases)
            {
                var st = SuffixTree.Build(text);
                var lrs = st.LongestRepeatedSubstring();

                // Check that no substring of length (lrs.Length + 1) repeats
                for (int i = 0; i <= text.Length - lrs.Length - 1; i++)
                {
                    var longer = text.Substring(i, lrs.Length + 1);
                    var count = st.CountOccurrences(longer);
                    if (count >= 2)
                    {
                        Assert.Fail($"Found longer repeating substring '{longer}' in '{text}' with count {count}");
                    }
                }
            }
        }

        #endregion

    }
}
