using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Tests for LongestCommonSubstring (LCS) functionality.
    /// Tests both single result and all results variants.
    /// </summary>
    [TestFixture]
    public class LongestCommonSubstringTests
    {
        #region Input Validation

        [Test]
        public void LCS_NullOther_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.LongestCommonSubstring(null!));
        }

        [Test]
        public void LCS_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");

            var result = st.LongestCommonSubstring("abc");

            Assert.That(result, Is.EqualTo(string.Empty));
        }

        [Test]
        public void LCS_EmptyOther_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");

            var result = st.LongestCommonSubstring("");

            Assert.That(result, Is.EqualTo(string.Empty));
        }

        [Test]
        public void LCS_BothEmpty_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");

            var result = st.LongestCommonSubstring("");

            Assert.That(result, Is.EqualTo(string.Empty));
        }

        #endregion

        #region Basic Cases

        [Test]
        public void LCS_IdenticalStrings_ReturnsFullString()
        {
            var text = "hello";
            var st = SuffixTree.Build(text);

            var result = st.LongestCommonSubstring(text);

            Assert.That(result, Is.EqualTo(text));
        }

        [Test]
        public void LCS_NoCommonCharacters_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");

            var result = st.LongestCommonSubstring("xyz");

            Assert.That(result, Is.EqualTo(string.Empty));
        }

        [Test]
        public void LCS_SingleCommonCharacter_ReturnsThatChar()
        {
            var st = SuffixTree.Build("abc");

            var result = st.LongestCommonSubstring("xay");

            Assert.That(result, Is.EqualTo("a"));
        }

        [Test]
        public void LCS_CommonPrefix_FindsIt()
        {
            var st = SuffixTree.Build("abcdef");

            var result = st.LongestCommonSubstring("abcxyz");

            Assert.That(result, Is.EqualTo("abc"));
        }

        [Test]
        public void LCS_CommonSuffix_FindsIt()
        {
            var st = SuffixTree.Build("xyzdef");

            var result = st.LongestCommonSubstring("abcdef");

            Assert.That(result, Is.EqualTo("def"));
        }

        [Test]
        public void LCS_CommonMiddle_FindsIt()
        {
            var st = SuffixTree.Build("xxxabcyyy");

            var result = st.LongestCommonSubstring("zzzabcwww");

            Assert.That(result, Is.EqualTo("abc"));
        }

        #endregion

        #region Classic Examples

        [Test]
        public void LCS_ClassicExample_OldSiteAndSite()
        {
            var st = SuffixTree.Build("OLDSITE:MREXSITE:MAPS");

            var result = st.LongestCommonSubstring("SITE:MR");

            Assert.That(result, Is.EqualTo("SITE:MR"));
        }

        [Test]
        public void LCS_ABAB_BABA()
        {
            var st = SuffixTree.Build("ABAB");

            var result = st.LongestCommonSubstring("BABA");

            // "ABA" or "BAB" - both length 3
            Assert.That(result.Length, Is.EqualTo(3));
            Assert.That(result, Is.AnyOf("ABA", "BAB"));
        }

        [Test]
        public void LCS_xabxac_abcabxabcd()
        {
            var st = SuffixTree.Build("xabxac");

            var result = st.LongestCommonSubstring("abcabxabcd");

            Assert.That(result, Is.EqualTo("abxa"));
        }

        #endregion

        #region Asymmetric Cases

        [Test]
        public void LCS_TreeShorterThanOther_Works()
        {
            var st = SuffixTree.Build("abc");

            var result = st.LongestCommonSubstring("xxabcxxabcxx");

            Assert.That(result, Is.EqualTo("abc"));
        }

        [Test]
        public void LCS_OtherShorterThanTree_Works()
        {
            var st = SuffixTree.Build("xxabcxxabcxx");

            var result = st.LongestCommonSubstring("abc");

            Assert.That(result, Is.EqualTo("abc"));
        }

        [Test]
        public void LCS_OtherIsSubstringOfTree_ReturnsOther()
        {
            var st = SuffixTree.Build("hello world");

            var result = st.LongestCommonSubstring("world");

            Assert.That(result, Is.EqualTo("world"));
        }

        [Test]
        public void LCS_TreeIsSubstringOfOther_ReturnsTree()
        {
            var st = SuffixTree.Build("world");

            var result = st.LongestCommonSubstring("hello world today");

            Assert.That(result, Is.EqualTo("world"));
        }

        #endregion

        #region Backtracking Edge Cases

        [Test]
        public void LCS_DeepBacktracking_FindsCorrectResult()
        {
            // Force deep matching then backtracking
            var st = SuffixTree.Build("xyzabcdefgh");

            var result = st.LongestCommonSubstring("abcXdefgh");

            // "defgh" is longer than "abc"
            Assert.That(result, Is.EqualTo("defgh"));
        }

        #endregion

        #region FindAllLongestCommonSubstrings

        [Test]
        public void FindAllLCS_NullOther_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.FindAllLongestCommonSubstrings(null!));
        }

        [Test]
        public void FindAllLCS_SingleResult_ReturnsSingle()
        {
            var st = SuffixTree.Build("abcdef");

            var (substring, positionsInText, positionsInOther) = st.FindAllLongestCommonSubstrings("xxabcxx");

            Assert.Multiple(() =>
            {
                Assert.That(substring, Is.EqualTo("abc"));
                Assert.That(positionsInText.Count, Is.GreaterThan(0));
                Assert.That(positionsInOther.Count, Is.GreaterThan(0));
            });
        }

        [Test]
        public void FindAllLCS_ReturnsPositions()
        {
            var st = SuffixTree.Build("abcdef");

            var (substring, positionsInText, positionsInOther) = st.FindAllLongestCommonSubstrings("xxabcxx");

            Assert.Multiple(() =>
            {
                Assert.That(substring, Is.EqualTo("abc"));
                Assert.That(positionsInText, Does.Contain(0)); // "abc" at position 0 in "abcdef"
                Assert.That(positionsInOther, Does.Contain(2)); // "abc" at position 2 in "xxabcxx"
            });
        }

        [Test]
        public void FindAllLCS_NoCommon_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");

            var (substring, positionsInText, positionsInOther) = st.FindAllLongestCommonSubstrings("xyz");

            Assert.Multiple(() =>
            {
                Assert.That(substring, Is.EqualTo(string.Empty));
                Assert.That(positionsInText, Is.Empty);
                Assert.That(positionsInOther, Is.Empty);
            });
        }

        [Test]
        public void FindAllLCS_IdenticalStrings_ReturnsFullString()
        {
            var st = SuffixTree.Build("test");

            var (substring, positionsInText, positionsInOther) = st.FindAllLongestCommonSubstrings("test");

            Assert.Multiple(() =>
            {
                Assert.That(substring, Is.EqualTo("test"));
                Assert.That(positionsInText, Does.Contain(0));
                Assert.That(positionsInOther, Does.Contain(0));
            });
        }

        #endregion

        #region LongestCommonSubstringInfo

        [Test]
        public void LCSInfo_NullOther_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.LongestCommonSubstringInfo(null!));
        }

        [Test]
        public void LCSInfo_ReturnsCorrectPositions()
        {
            var st = SuffixTree.Build("abcdef");

            var (substring, positionInText, positionInOther) = st.LongestCommonSubstringInfo("xxcdexx");

            Assert.Multiple(() =>
            {
                Assert.That(substring, Is.EqualTo("cde"));
                Assert.That(positionInText, Is.EqualTo(2)); // "cde" at position 2 in "abcdef"
                Assert.That(positionInOther, Is.EqualTo(2)); // "cde" at position 2 in "xxcdexx"
            });
        }

        [Test]
        public void LCSInfo_NoCommon_ReturnsMinusOne()
        {
            var st = SuffixTree.Build("abc");

            var (substring, positionInText, positionInOther) = st.LongestCommonSubstringInfo("xyz");

            Assert.Multiple(() =>
            {
                Assert.That(substring, Is.EqualTo(string.Empty));
                Assert.That(positionInText, Is.EqualTo(-1));
                Assert.That(positionInOther, Is.EqualTo(-1));
            });
        }

        [Test]
        public void LCSInfo_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");

            var (substring, positionInText, positionInOther) = st.LongestCommonSubstringInfo("abc");

            Assert.Multiple(() =>
            {
                Assert.That(substring, Is.EqualTo(string.Empty));
                Assert.That(positionInText, Is.EqualTo(-1));
                Assert.That(positionInOther, Is.EqualTo(-1));
            });
        }

        #endregion

        #region Result Validation

        [Test]
        public void LCS_ResultActuallyInBothStrings()
        {
            var testCases = new[]
            {
                ("hello world", "world hello"),
                ("abcdefgh", "xyzabcxyz"),
                ("mississippi", "missouri"),
            };

            foreach (var (text, other) in testCases)
            {
                var st = SuffixTree.Build(text);
                var lcs = st.LongestCommonSubstring(other);

                Assert.Multiple(() =>
                {
                    Assert.That(text.Contains(lcs), Is.True,
                        $"LCS '{lcs}' not in tree text '{text}'");
                    Assert.That(other.Contains(lcs), Is.True,
                        $"LCS '{lcs}' not in other string '{other}'");
                });
            }
        }

        [Test]
        public void LCS_IsActuallyLongest()
        {
            var testCases = new[]
            {
                ("hello world", "world hello"),
                ("abcdefgh", "xyzabcxyz"),
            };

            foreach (var (text, other) in testCases)
            {
                var st = SuffixTree.Build(text);
                var lcs = st.LongestCommonSubstring(other);

                // Verify no longer common substring exists
                for (int i = 0; i <= text.Length - lcs.Length - 1; i++)
                {
                    var longer = text.Substring(i, lcs.Length + 1);
                    if (other.Contains(longer))
                    {
                        Assert.Fail($"Found longer common substring '{longer}' between '{text}' and '{other}'");
                    }
                }
            }
        }

        #endregion

        #region Span Overload

        [Test]
        public void LCS_SpanOverload_MatchesStringOverload()
        {
            var st = SuffixTree.Build("hello world");

            var others = new[] { "world hello", "xyz", "lo wo", "" };

            foreach (var other in others)
            {
                Assert.That(
                    st.LongestCommonSubstring(other.AsSpan()),
                    Is.EqualTo(st.LongestCommonSubstring(other)),
                    $"Span and string overloads should match for '{other}'");
            }
        }

        [Test]
        public void LCS_SpanFromSlice_Works()
        {
            var st = SuffixTree.Build("hello world");
            var text = "xxxworldxxx";

            var result = st.LongestCommonSubstring(text.AsSpan(3, 5));

            Assert.That(result, Is.EqualTo("world"));
        }

        [Test]
        public void LCS_EmptySpan_ReturnsEmpty()
        {
            var st = SuffixTree.Build("hello");

            var result = st.LongestCommonSubstring(ReadOnlySpan<char>.Empty);

            Assert.That(result, Is.EqualTo(string.Empty));
        }

        #endregion

        #region Edge Cases

        [Test]
        public void LCS_SingleCharStrings_Works()
        {
            var st = SuffixTree.Build("a");

            Assert.Multiple(() =>
            {
                Assert.That(st.LongestCommonSubstring("a"), Is.EqualTo("a"));
                Assert.That(st.LongestCommonSubstring("b"), Is.EqualTo(string.Empty));
            });
        }

        [Test]
        public void LCS_RepeatingPatterns_Works()
        {
            var st = SuffixTree.Build("ababab");

            var result = st.LongestCommonSubstring("bababa");

            // "ababa" or "babab" - length 5
            Assert.That(result.Length, Is.EqualTo(5));
        }

        [Test]
        public void LCS_LongStrings_Works()
        {
            var pattern = "commonpart";
            var st = SuffixTree.Build(new string('x', 1000) + pattern + new string('y', 1000));

            var result = st.LongestCommonSubstring(new string('a', 500) + pattern + new string('b', 500));

            Assert.That(result, Is.EqualTo(pattern));
        }

        #endregion
    }
}
