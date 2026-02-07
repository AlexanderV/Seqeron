using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Robustness
{
    /// <summary>
    /// Edge case tests for boundary conditions and unusual inputs.
    /// </summary>
    [TestFixture]
    public class EdgeCaseTests
    {
        #region Very Long Patterns

        [Test]
        public void Contains_PatternLongerThanText_ReturnsFalse()
        {
            var st = SuffixTree.Build("abc");

            Assert.That(st.Contains("abcdefgh"), Is.False);
        }

        [Test]
        public void FindAll_PatternLongerThanText_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");

            var result = st.FindAllOccurrences("abcdefgh").ToList();

            Assert.That(result, Is.Empty);
        }

        [Test]
        public void LCS_OtherMuchLongerThanTree_Works()
        {
            var st = SuffixTree.Build("abc");
            var other = new string('x', 10000) + "abc" + new string('y', 10000);

            var result = st.LongestCommonSubstring(other);

            Assert.That(result, Is.EqualTo("abc"));
        }

        #endregion

        #region Pattern Equals Text

        [Test]
        public void Contains_PatternEqualsText_ReturnsTrue()
        {
            var text = "exactly this text";
            var st = SuffixTree.Build(text);

            Assert.That(st.Contains(text), Is.True);
        }

        [Test]
        public void FindAll_PatternEqualsText_ReturnsZero()
        {
            var text = "exactly";
            var st = SuffixTree.Build(text);

            var result = st.FindAllOccurrences(text).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0], Is.EqualTo(0));
            });
        }

        #endregion

        #region Single Character Strings

        [Test]
        public void SingleCharText_AllOperationsWork()
        {
            var st = SuffixTree.Build("x");

            Assert.Multiple(() =>
            {
                Assert.That(st.Text, Is.EqualTo("x"));
                Assert.That(st.Contains("x"), Is.True);
                Assert.That(st.Contains("y"), Is.False);
                Assert.That(st.FindAllOccurrences("x").ToList(), Is.EqualTo(new[] { 0 }));
                Assert.That(st.CountOccurrences("x"), Is.EqualTo(1));
                Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo(string.Empty));
                Assert.That(st.LongestCommonSubstring("x"), Is.EqualTo("x"));
                Assert.That(st.LongestCommonSubstring("y"), Is.EqualTo(string.Empty));
                Assert.That(st.GetAllSuffixes().ToList(), Is.EqualTo(new[] { "x" }));
            });
        }

        #endregion

        #region Whitespace Only

        [Test]
        public void WhitespaceOnlyText_Works()
        {
            var st = SuffixTree.Build("     ");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains(" "), Is.True);
                Assert.That(st.Contains("  "), Is.True);
                Assert.That(st.Contains("     "), Is.True);
                Assert.That(st.FindAllOccurrences(" ").Count(), Is.EqualTo(5));
            });
        }

        [Test]
        public void MixedWhitespace_Works()
        {
            var st = SuffixTree.Build(" \t\n\r ");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains(" "), Is.True);
                Assert.That(st.Contains("\t"), Is.True);
                Assert.That(st.Contains("\n"), Is.True);
                Assert.That(st.Contains("\r"), Is.True);
                Assert.That(st.Contains(" \t\n"), Is.True);
            });
        }

        #endregion

        #region Alternating Patterns

        [Test]
        public void AlternatingPattern_Works()
        {
            var st = SuffixTree.Build("ababababab");

            Assert.Multiple(() =>
            {
                Assert.That(st.FindAllOccurrences("ab").Count(), Is.EqualTo(5));
                Assert.That(st.FindAllOccurrences("ba").Count(), Is.EqualTo(4));
                Assert.That(st.FindAllOccurrences("aba").Count(), Is.EqualTo(4));
                Assert.That(st.FindAllOccurrences("bab").Count(), Is.EqualTo(4));
            });
        }

        #endregion

        #region Palindromes

        [Test]
        public void Palindrome_Works()
        {
            var st = SuffixTree.Build("abcba");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("abcba"), Is.True);
                Assert.That(st.Contains("bcb"), Is.True);
                Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("a").Or.EqualTo("b"));
            });
        }

        [Test]
        public void LongPalindrome_Works()
        {
            var half = "abcdefg";
            var palindrome = half + new string(half.Reverse().ToArray());
            var st = SuffixTree.Build(palindrome);

            Assert.That(st.Contains(palindrome), Is.True);
        }

        #endregion

        #region Special String Patterns

        [Test]
        public void RunLengthEncoded_Works()
        {
            var st = SuffixTree.Build("aabbccddee");

            Assert.Multiple(() =>
            {
                Assert.That(st.FindAllOccurrences("aa").Count(), Is.EqualTo(1));
                Assert.That(st.FindAllOccurrences("bb").Count(), Is.EqualTo(1));
                Assert.That(st.Contains("aabb"), Is.True);
                Assert.That(st.Contains("bbcc"), Is.True);
            });
        }

        #endregion
    }
}
