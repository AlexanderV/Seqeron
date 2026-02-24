using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Compatibility
{
    /// <summary>
    /// Tests for Unicode support across all operations.
    /// </summary>
    [TestFixture]
    [Category("Compatibility")]
    public class UnicodeTests
    {
        #region Construction

        [Test]
        public void Build_WithCyrillicText_Works()
        {
            var st = SuffixTree.Build("привет мир");

            Assert.Multiple(() =>
            {
                Assert.That(st.Text, Is.EqualTo("привет мир"));
                Assert.That(st.Contains("привет"), Is.True);
                Assert.That(st.Contains("мир"), Is.True);
            });
        }

        [Test]
        public void Build_WithChineseText_Works()
        {
            var st = SuffixTree.Build("你好世界");

            Assert.Multiple(() =>
            {
                Assert.That(st.Text, Is.EqualTo("你好世界"));
                Assert.That(st.Contains("你好"), Is.True);
                Assert.That(st.Contains("世界"), Is.True);
            });
        }

        [Test]
        public void Build_WithGreekText_Works()
        {
            var st = SuffixTree.Build("αβγδεζ");

            Assert.Multiple(() =>
            {
                Assert.That(st.Text, Is.EqualTo("αβγδεζ"));
                Assert.That(st.Contains("αβγ"), Is.True);
                Assert.That(st.Contains("δεζ"), Is.True);
            });
        }

        [Test]
        public void Build_WithMixedScripts_Works()
        {
            var text = "Hello世界Привет";
            var st = SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("Hello"), Is.True);
                Assert.That(st.Contains("世界"), Is.True);
                Assert.That(st.Contains("Привет"), Is.True);
                Assert.That(st.Contains("o世界П"), Is.True, "Cross-script substring");
            });
        }

        #endregion

        #region Emoji Support

        [Test]
        public void Build_WithEmoji_Works()
        {
            // Note: Some emoji are surrogate pairs (2 chars in .NET)
            var st = SuffixTree.Build("hello😀world");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("hello"), Is.True);
                Assert.That(st.Contains("world"), Is.True);
                Assert.That(st.Contains("😀"), Is.True);
                Assert.That(st.Contains("o😀w"), Is.True);
            });
        }

        [Test]
        public void Build_WithSurrogatePairs_Works()
        {
            // 𝄞 = musical G clef, is a surrogate pair (U+1D11E)
            var st = SuffixTree.Build("abc𝄞def");

            Assert.That(st.Contains("𝄞"), Is.True);
        }

        #endregion

        #region FindAllOccurrences

        [Test]
        public void FindAll_UnicodePattern_FindsCorrectPositions()
        {
            var st = SuffixTree.Build("αβαβα");

            var positions = st.FindAllOccurrences("αβ").OrderBy(x => x).ToList();

            Assert.That(positions, Is.EqualTo(new[] { 0, 2 }));
        }

        [Test]
        public void FindAll_EmojiPattern_Works()
        {
            var st = SuffixTree.Build("a😀b😀c");

            var positions = st.FindAllOccurrences("😀").OrderBy(x => x).ToList();

            // Position in chars (emoji is 2 chars as surrogate pair)
            Assert.That(positions, Has.Count.EqualTo(2));
        }

        #endregion

        #region LCS with Unicode

        [Test]
        public void LCS_Unicode_FindsCorrectSubstring()
        {
            var st = SuffixTree.Build("αβγδεζ");

            var result = st.LongestCommonSubstring("xxγδεxx");

            Assert.That(result, Is.EqualTo("γδε"));
        }

        [Test]
        public void LCS_MixedScripts_Works()
        {
            var st = SuffixTree.Build("hello世界test");

            var result = st.LongestCommonSubstring("xxx世界yyy");

            Assert.That(result, Is.EqualTo("世界"));
        }

        #endregion

        #region LRS with Unicode

        [Test]
        public void LRS_Unicode_FindsRepeatedSubstring()
        {
            var st = SuffixTree.Build("αβγαβγ");

            var result = st.LongestRepeatedSubstring();

            Assert.That(result, Is.EqualTo("αβγ"));
        }

        #endregion

        #region GetAllSuffixes

        [Test]
        public void GetAllSuffixes_Unicode_ReturnsAllSuffixes()
        {
            var text = "αβγ";
            var st = SuffixTree.Build(text);

            var suffixes = st.GetAllSuffixes().OrderBy(s => s).ToList();

            Assert.That(suffixes, Has.Count.EqualTo(3));
            Assert.That(suffixes, Does.Contain("αβγ"));
            Assert.That(suffixes, Does.Contain("βγ"));
            Assert.That(suffixes, Does.Contain("γ"));
        }

        #endregion

        #region High Unicode Code Points

        [Test]
        public void Build_HighCodePoints_Works()
        {
            // Characters outside BMP (require surrogate pairs)
            var text = "a\U0001F600b"; // 😀
            var st = SuffixTree.Build(text);

            Assert.That(st.Contains("\U0001F600"), Is.True);
        }

        [Test]
        public void Build_MathematicalSymbols_Works()
        {
            var st = SuffixTree.Build("∑∏∫∂∇");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("∑"), Is.True);
                Assert.That(st.Contains("∏∫"), Is.True);
                Assert.That(st.Contains("∂∇"), Is.True);
            });
        }

        #endregion

        #region Zero-Width Characters

        [Test]
        public void Build_WithZeroWidthChars_PreservesThem()
        {
            // Zero-width joiner
            var zwj = "\u200D";
            var text = $"a{zwj}b";
            var st = SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(st.Text, Is.EqualTo(text));
                Assert.That(st.Contains(zwj), Is.True);
                Assert.That(st.Contains($"a{zwj}"), Is.True);
            });
        }

        #endregion

        #region Normalization

        [Test]
        public void Build_DifferentNormalizationForms_AreDistinct()
        {
            // é can be represented as single char (U+00E9) or as e + combining acute (U+0065 U+0301)
            var composed = "é";    // Single char
            var decomposed = "e\u0301"; // Two chars

            // Note: SuffixTree treats these as different (no normalization)
            var st = SuffixTree.Build(composed);

            // They are different at the character level
            Assert.That(st.Contains(decomposed), Is.EqualTo(composed == decomposed));
        }

        #endregion
    }
}

