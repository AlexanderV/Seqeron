using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Compatibility
{
    /// <summary>
    /// Tests for binary data support (null characters and other control characters).
    /// </summary>
    [TestFixture]
    public class BinaryDataTests
    {
        #region Null Character Construction

        [Test]
        public void Build_WithNullCharacter_Works()
        {
            var text = "a\0b";
            var st = SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(st.Text, Is.EqualTo(text));
                Assert.That(st.Text.Length, Is.EqualTo(3));
            });
        }

        [Test]
        public void Build_StartsWithNull_Works()
        {
            var text = "\0abc";
            var st = SuffixTree.Build(text);

            Assert.That(st.Text, Is.EqualTo(text));
        }

        [Test]
        public void Build_EndsWithNull_Works()
        {
            var text = "abc\0";
            var st = SuffixTree.Build(text);

            Assert.That(st.Text, Is.EqualTo(text));
        }

        [Test]
        public void Build_MultipleNulls_Works()
        {
            var text = "\0\0\0";
            var st = SuffixTree.Build(text);

            Assert.That(st.Text, Is.EqualTo(text));
        }

        [Test]
        public void Build_OnlyNull_Works()
        {
            var text = "\0";
            var st = SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(st.Text, Is.EqualTo(text));
                Assert.That(st.Contains("\0"), Is.True);
            });
        }

        #endregion

        #region Contains with Null

        [Test]
        public void Contains_NullCharacter_FindsIt()
        {
            var st = SuffixTree.Build("ab\0cd");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("\0"), Is.True);
                Assert.That(st.Contains("b\0c"), Is.True);
                Assert.That(st.Contains("\0cd"), Is.True);
                Assert.That(st.Contains("ab\0"), Is.True);
            });
        }

        [Test]
        public void Contains_NullInPattern_NotInTree_ReturnsFalse()
        {
            var st = SuffixTree.Build("abcdef");

            Assert.That(st.Contains("ab\0cd"), Is.False);
        }

        #endregion

        #region FindAllOccurrences with Null

        [Test]
        public void FindAll_NullPattern_Works()
        {
            var st = SuffixTree.Build("a\0b\0c");

            var positions = st.FindAllOccurrences("\0").OrderBy(x => x).ToList();

            Assert.That(positions, Is.EqualTo(new[] { 1, 3 }));
        }

        [Test]
        public void FindAll_PatternWithNull_Works()
        {
            var st = SuffixTree.Build("a\0ba\0b");

            var positions = st.FindAllOccurrences("a\0b").OrderBy(x => x).ToList();

            Assert.That(positions, Is.EqualTo(new[] { 0, 3 }));
        }

        #endregion

        #region LCS with Null

        [Test]
        public void LCS_WithNullCharacters_Works()
        {
            var st = SuffixTree.Build("ab\0cd\0ef");

            var result = st.LongestCommonSubstring("xx\0cd\0yy");

            Assert.That(result, Is.EqualTo("\0cd\0"));
        }

        [Test]
        public void LCS_NullAsOnlyCommon_Works()
        {
            var st = SuffixTree.Build("a\0b");

            var result = st.LongestCommonSubstring("x\0y");

            Assert.That(result, Is.EqualTo("\0"));
        }

        [Test]
        public void LCS_BothContainNull_NoCommonPart()
        {
            var st = SuffixTree.Build("abc\0xyz");

            var result = st.LongestCommonSubstring("def\0uvw");

            Assert.That(result, Is.EqualTo("\0"));
        }

        #endregion

        #region LRS with Null

        [Test]
        public void LRS_WithNullCharacters_Works()
        {
            var st = SuffixTree.Build("a\0ba\0b");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("a\0b"));
        }

        [Test]
        public void LRS_RepeatedNulls_Works()
        {
            var st = SuffixTree.Build("\0\0\0");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("\0\0"));
        }

        #endregion

        #region Control Characters

        [Test]
        public void Build_WithControlCharacters_Works()
        {
            var text = "a\x01b";
            var st = SuffixTree.Build(text);

            Assert.That(st.Contains("a\x01b"), Is.True);
        }

        [Test]
        public void Build_WithTab_Works()
        {
            var st = SuffixTree.Build("hello\tworld");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("\t"), Is.True);
                Assert.That(st.Contains("o\tw"), Is.True);
            });
        }

        [Test]
        public void Build_WithNewline_Works()
        {
            var st = SuffixTree.Build("line1\nline2\r\nline3");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("\n"), Is.True);
                Assert.That(st.Contains("\r\n"), Is.True);
                Assert.That(st.Contains("1\nline2"), Is.True);
            });
        }

        #endregion

        #region Binary-Like Data

        [Test]
        public void Build_AllByteValues_Works()
        {
            // Create a string with all possible char values 0-255
            var chars = Enumerable.Range(0, 256).Select(i => (char)i).ToArray();
            var text = new string(chars);
            var st = SuffixTree.Build(text);

            // Verify all single characters can be found
            for (int i = 0; i < 256; i++)
            {
                Assert.That(st.Contains(((char)i).ToString()), Is.True,
                    $"Character {i} (0x{i:X2}) not found");
            }
        }

        [Test]
        public void Build_HighControlCharacters_Works()
        {
            // DEL character
            var text = "abc\x7Fdef";
            var st = SuffixTree.Build(text);

            Assert.That(st.Contains("abc\x7Fdef"), Is.True);
        }

        #endregion

        #region GetAllSuffixes with Null

        [Test]
        public void GetAllSuffixes_WithNull_ReturnsAllSuffixes()
        {
            var text = "a\0b";
            var st = SuffixTree.Build(text);

            var suffixes = st.GetAllSuffixes().OrderBy(s => s).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(suffixes, Has.Count.EqualTo(3));
                Assert.That(suffixes, Does.Contain("a\0b"));
                Assert.That(suffixes, Does.Contain("\0b"));
                Assert.That(suffixes, Does.Contain("b"));
            });
        }

        #endregion
    }
}
