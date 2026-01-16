using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using NUnit.Framework;

namespace SuffixTree.Tests
{
    /// <summary>
    /// Comprehensive unit tests for SuffixTree implementation.
    /// 
    /// Test Hierarchy (115 tests total):
    /// 1. Build (13) - tree creation, validation, edge cases [incl. 6 TestCase variants]
    /// 2. Contains (11) - substring search with string and Span overloads
    /// 3. FindAllOccurrences/CountOccurrences (13) - position finding
    /// 4. LongestRepeatedSubstring (11) - LRS algorithm verification
    /// 5. LongestCommonSubstring (13) - LCS algorithm (all API variants)
    /// 6. GetAllSuffixes/EnumerateSuffixes (9) - suffix enumeration
    /// 7. Statistics (7) - NodeCount, LeafCount, MaxDepth
    /// 8. Unicode and Special Characters (9) - international support
    /// 9. Performance and Large String (5) - stress testing
    /// 10. Thread Safety (3) - concurrent access
    /// 11. Algorithm-Specific (9) - Ukkonen's edge cases
    /// 12. Regression (8) - historical bug prevention
    /// 13. Span API (4) - ReadOnlySpan overloads
    /// </summary>
    [TestFixture]
    public class SuffixTreeTests
    {
        private const int RANDOM_SEED = 42;

        #region 1. Build Tests

        [Test]
        public void Build_NullString_ThrowsArgumentNullException()
        {
            Assert.Throws<ArgumentNullException>(() => SuffixTree.Build(null!));
        }

        [Test]
        public void Build_StringWithTerminator_ThrowsArgumentException()
        {
            Assert.Throws<ArgumentException>(() => SuffixTree.Build("abc\0def"));
        }

        [Test]
        public void Build_EmptyString_CreatesValidTree()
        {
            var st = SuffixTree.Build("");

            Assert.That(st, Is.Not.Null);
            Assert.That(st.Text, Is.EqualTo(""));
            Assert.That(st.Contains(""), Is.True);
            Assert.That(st.Contains("a"), Is.False);
        }

        [Test]
        public void Build_SingleCharacter_CreatesValidTree()
        {
            var st = SuffixTree.Build("a");

            Assert.That(st.Text, Is.EqualTo("a"));
            Assert.That(st.Contains("a"), Is.True);
            Assert.That(st.Contains("b"), Is.False);
        }

        [TestCase("abcdefghij")]
        [TestCase("aaaaaaaaaa")]
        [TestCase("abababab")]
        [TestCase("mississippi")]
        [TestCase("banana")]
        [TestCase("abcabxabcd")]
        public void Build_VariousStrings_DoesNotThrow(string text)
        {
            Assert.DoesNotThrow(() => SuffixTree.Build(text));
        }

        [Test]
        public void Text_ReturnsOriginalString()
        {
            var original = "test string";
            var st = SuffixTree.Build(original);

            Assert.That(st.Text, Is.EqualTo(original));
        }

        [Test]
        public void ToString_ContainsMeaningfulInfo()
        {
            var st = SuffixTree.Build("hello");
            var result = st.ToString();

            Assert.That(result, Does.Contain("SuffixTree"));
            Assert.That(result, Does.Contain("hello"));
        }

        [Test]
        public void PrintTree_DoesNotThrow()
        {
            var st = SuffixTree.Build("banana");
            Assert.DoesNotThrow(() => st.PrintTree());
        }

        #endregion

        #region 2. Contains Tests

        [Test]
        public void Contains_NullPattern_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.Contains(null!));
        }

        [Test]
        public void Contains_PatternWithTerminator_ThrowsArgumentException()
        {
            var st = SuffixTree.Build("hello");
            Assert.Throws<ArgumentException>(() => st.Contains("hel\0lo".AsSpan()));
        }

        [Test]
        public void Contains_EmptyPattern_ReturnsTrue()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.Contains(""), Is.True);
            Assert.That(st.Contains(ReadOnlySpan<char>.Empty), Is.True);
        }

        [Test]
        public void Contains_EmptyTree_NonEmptyPattern_ReturnsFalse()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.Contains("a"), Is.False);
            Assert.That(st.Contains("abc"), Is.False);
        }

        [Test]
        public void Contains_FullString_ReturnsTrue()
        {
            var text = "teststring";
            var st = SuffixTree.Build(text);
            Assert.That(st.Contains(text), Is.True);
        }

        [Test]
        public void Contains_AllSuffixes_ReturnsTrue()
        {
            var text = "banana";
            var st = SuffixTree.Build(text);

            for (int i = 0; i < text.Length; i++)
            {
                var suffix = text.Substring(i);
                Assert.That(st.Contains(suffix), Is.True, $"Suffix '{suffix}' not found");
            }
        }

        [Test]
        public void Contains_AllSubstrings_ReturnsTrue()
        {
            var text = "abracadabra";
            var st = SuffixTree.Build(text);

            for (int i = 0; i < text.Length; i++)
            {
                for (int len = 1; len <= text.Length - i; len++)
                {
                    var substr = text.Substring(i, len);
                    Assert.That(st.Contains(substr), Is.True,
                        $"Substring '{substr}' at [{i},{i + len}) not found");
                }
            }
        }

        [Test]
        public void Contains_NonExistentPatterns_ReturnsFalse()
        {
            var st = SuffixTree.Build("abcdef");

            Assert.That(st.Contains("xyz"), Is.False);
            Assert.That(st.Contains("abd"), Is.False);
            Assert.That(st.Contains("abcdefg"), Is.False);
            Assert.That(st.Contains("aba"), Is.False);
        }

        [Test]
        public void Contains_SpanOverload_MatchesStringOverload()
        {
            var r = new Random(RANDOM_SEED);
            var text = MakeRandomString(r, 100);
            var st = SuffixTree.Build(text);

            for (int i = 0; i < 20; i++)
            {
                int pos = r.Next(0, text.Length - 10);
                int len = r.Next(1, 10);
                var pattern = text.Substring(pos, len);

                Assert.That(st.Contains(pattern.AsSpan()), Is.EqualTo(st.Contains(pattern)));
            }
        }

        [Test]
        public void Contains_OverlappingPatterns_Works()
        {
            var st = SuffixTree.Build("ababababab");

            Assert.That(st.Contains("abab"), Is.True);
            Assert.That(st.Contains("baba"), Is.True);
            Assert.That(st.Contains("ababa"), Is.True);
            Assert.That(st.Contains("ababababab"), Is.True);
            Assert.That(st.Contains("abababababa"), Is.False);
        }

        [Test]
        public void Contains_RepeatingCharacter_Works()
        {
            var st = SuffixTree.Build("aaaaaaaaaa");

            for (int len = 1; len <= 10; len++)
            {
                Assert.That(st.Contains(new string('a', len)), Is.True);
            }
            Assert.That(st.Contains(new string('a', 11)), Is.False);
            Assert.That(st.Contains("b"), Is.False);
        }

        #endregion

        #region 3. FindAllOccurrences / CountOccurrences Tests

        [Test]
        public void FindAllOccurrences_NullPattern_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.FindAllOccurrences(null!));
        }

        [Test]
        public void FindAllOccurrences_EmptyPattern_ReturnsAllPositions()
        {
            var st = SuffixTree.Build("abc");
            var result = st.FindAllOccurrences("");
            Assert.That(result, Is.EquivalentTo(new[] { 0, 1, 2 }));
        }

        [Test]
        public void FindAllOccurrences_NotFound_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abcdef");
            Assert.That(st.FindAllOccurrences("xyz"), Is.Empty);
        }

        [Test]
        public void FindAllOccurrences_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.FindAllOccurrences("a"), Is.Empty);
        }

        [Test]
        public void FindAllOccurrences_SingleOccurrence_ReturnsCorrectPosition()
        {
            var st = SuffixTree.Build("abcdef");
            var result = st.FindAllOccurrences("cde");
            Assert.That(result, Is.EquivalentTo(new[] { 2 }));
        }

        [Test]
        public void FindAllOccurrences_MultipleOccurrences_ReturnsAll()
        {
            var st = SuffixTree.Build("abcabc");
            Assert.That(st.FindAllOccurrences("abc"), Is.EquivalentTo(new[] { 0, 3 }));
        }

        [Test]
        public void FindAllOccurrences_OverlappingOccurrences_ReturnsAll()
        {
            var st = SuffixTree.Build("aaa");
            Assert.That(st.FindAllOccurrences("aa"), Is.EquivalentTo(new[] { 0, 1 }));
        }

        [Test]
        public void FindAllOccurrences_Banana_Ana_ReturnsBothPositions()
        {
            var st = SuffixTree.Build("banana");
            Assert.That(st.FindAllOccurrences("ana"), Is.EquivalentTo(new[] { 1, 3 }));
        }

        [Test]
        public void FindAllOccurrences_MatchesNaiveImplementation()
        {
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < 50; i++)
            {
                var text = MakeRandomString(r, 100);
                var st = SuffixTree.Build(text);

                int pos = r.Next(0, text.Length - 10);
                int len = r.Next(1, 10);
                var pattern = text.Substring(pos, len);

                var stResult = st.FindAllOccurrences(pattern).OrderBy(x => x);
                var naiveResult = NaiveFindAll(text, pattern).OrderBy(x => x);

                Assert.That(stResult, Is.EqualTo(naiveResult));
            }
        }

        [Test]
        public void CountOccurrences_NullPattern_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.CountOccurrences(null!));
        }

        [Test]
        public void CountOccurrences_MatchesFindAllOccurrencesCount()
        {
            var st = SuffixTree.Build("banana");

            Assert.That(st.CountOccurrences("ana"), Is.EqualTo(st.FindAllOccurrences("ana").Count));
            Assert.That(st.CountOccurrences("a"), Is.EqualTo(st.FindAllOccurrences("a").Count));
            Assert.That(st.CountOccurrences("xyz"), Is.EqualTo(st.FindAllOccurrences("xyz").Count));
        }

        [Test]
        public void CountOccurrences_OverlappingPattern_ReturnsCorrectCount()
        {
            var st = SuffixTree.Build("aaaa");
            Assert.That(st.CountOccurrences("aa"), Is.EqualTo(3));
        }

        [Test]
        public void CountOccurrences_EmptyPattern_ReturnsTextLength()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.CountOccurrences(""), Is.EqualTo(3));
        }

        #endregion

        #region 4. LongestRepeatedSubstring Tests

        [Test]
        public void LongestRepeatedSubstring_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo(""));
        }

        [Test]
        public void LongestRepeatedSubstring_SingleChar_ReturnsEmpty()
        {
            var st = SuffixTree.Build("a");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo(""));
        }

        [Test]
        public void LongestRepeatedSubstring_NoRepeats_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abcdef");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo(""));
        }

        [Test]
        public void LongestRepeatedSubstring_SimpleRepeat_ReturnsCorrect()
        {
            var st = SuffixTree.Build("abcabc");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("abc"));
        }

        [Test]
        public void LongestRepeatedSubstring_Banana_ReturnsAna()
        {
            var st = SuffixTree.Build("banana");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("ana"));
        }

        [Test]
        public void LongestRepeatedSubstring_Mississippi_ReturnsIssi()
        {
            var st = SuffixTree.Build("mississippi");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("issi"));
        }

        [Test]
        public void LongestRepeatedSubstring_OverlappingRepeat_ReturnsCorrect()
        {
            var st = SuffixTree.Build("aaa");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("aa"));
        }

        [Test]
        public void LongestRepeatedSubstring_Result_AppearsAtLeastTwice()
        {
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < 30; i++)
            {
                var text = MakeRandomString(r, 100);
                var st = SuffixTree.Build(text);
                var lrs = st.LongestRepeatedSubstring();

                if (lrs.Length > 0)
                {
                    Assert.That(st.CountOccurrences(lrs), Is.GreaterThanOrEqualTo(2),
                        $"LRS '{lrs}' should appear at least twice in '{text}'");
                }
            }
        }

        [Test]
        public void LongestRepeatedSubstring_IsActuallyLongest()
        {
            var r = new Random(RANDOM_SEED);

            for (int i = 0; i < 20; i++)
            {
                var text = MakeRandomString(r, 50);
                var st = SuffixTree.Build(text);
                var lrs = st.LongestRepeatedSubstring();

                // No longer repeated substring should exist
                if (lrs.Length > 0 && lrs.Length < text.Length - 1)
                {
                    for (int j = 0; j < text.Length - lrs.Length; j++)
                    {
                        var candidate = text.Substring(j, lrs.Length + 1);
                        Assert.That(st.CountOccurrences(candidate), Is.LessThan(2),
                            $"Found longer repeated substring '{candidate}' than LRS '{lrs}'");
                    }
                }
            }
        }

        [Test]
        public void LongestRepeatedSubstring_ABAB_Pattern()
        {
            var st = SuffixTree.Build("abab");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("ab"));
        }

        [Test]
        public void LongestRepeatedSubstring_GEEKSFORGEEKS()
        {
            var st = SuffixTree.Build("GEEKSFORGEEKS");
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("GEEKS"));
        }

        #endregion

        #region 5. LongestCommonSubstring Tests

        [Test]
        public void LongestCommonSubstring_NullOther_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.LongestCommonSubstring(null!));
        }

        [Test]
        public void LongestCommonSubstring_OtherWithTerminator_ThrowsArgumentException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentException>(() => st.LongestCommonSubstring("a\0b"));
        }

        [Test]
        public void LongestCommonSubstring_EmptyOther_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.LongestCommonSubstring(""), Is.EqualTo(""));
        }

        [Test]
        public void LongestCommonSubstring_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.LongestCommonSubstring("abc"), Is.EqualTo(""));
        }

        [Test]
        public void LongestCommonSubstring_NoCommon_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");
            Assert.That(st.LongestCommonSubstring("xyz"), Is.EqualTo(""));
        }

        [Test]
        public void LongestCommonSubstring_FullMatch_ReturnsFullString()
        {
            var st = SuffixTree.Build("abcdef");
            Assert.That(st.LongestCommonSubstring("abcdef"), Is.EqualTo("abcdef"));
        }

        [Test]
        public void LongestCommonSubstring_PartialMatch_ReturnsLongestCommon()
        {
            var st = SuffixTree.Build("xyzabcdef");
            Assert.That(st.LongestCommonSubstring("123abc456"), Is.EqualTo("abc"));
        }

        [Test]
        public void LongestCommonSubstring_MultipleMatches_ReturnsLongest()
        {
            var st = SuffixTree.Build("abcdefghijklmnop");
            Assert.That(st.LongestCommonSubstring("xxabcxxdefghixx"), Is.EqualTo("defghi"));
        }

        [Test]
        public void LongestCommonSubstringInfo_ReturnsCorrectPositions()
        {
            var st = SuffixTree.Build("abracadabra");
            var (substring, posInText, posInOther) = st.LongestCommonSubstringInfo("xxcadyy");

            Assert.That(substring, Is.EqualTo("cad"));
            Assert.That(posInText, Is.EqualTo(4));
            Assert.That(posInOther, Is.EqualTo(2));
        }

        [Test]
        public void LongestCommonSubstringInfo_NoMatch_ReturnsNegativePositions()
        {
            var st = SuffixTree.Build("abc");
            var (substring, posInText, posInOther) = st.LongestCommonSubstringInfo("xyz");

            Assert.That(substring, Is.EqualTo(""));
            Assert.That(posInText, Is.EqualTo(-1));
            Assert.That(posInOther, Is.EqualTo(-1));
        }

        [Test]
        public void FindAllLongestCommonSubstrings_MultipleMatches_ReturnsAllPositions()
        {
            var st = SuffixTree.Build("abcabc");
            var result = st.FindAllLongestCommonSubstrings("xabcyabcz");

            Assert.That(result.Substring, Is.EqualTo("abc"));
            Assert.That(result.PositionsInOther, Has.Count.EqualTo(2));
            Assert.That(result.PositionsInOther, Does.Contain(1));
            Assert.That(result.PositionsInOther, Does.Contain(5));
        }

        [Test]
        public void FindAllLongestCommonSubstrings_NoMatch_ReturnsEmpty()
        {
            var st = SuffixTree.Build("abc");
            var result = st.FindAllLongestCommonSubstrings("xyz");

            Assert.That(result.Substring, Is.EqualTo(""));
            Assert.That(result.PositionsInText, Is.Empty);
            Assert.That(result.PositionsInOther, Is.Empty);
        }

        [Test]
        public void FindAllLongestCommonSubstrings_NullOther_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.FindAllLongestCommonSubstrings(null!));
        }

        #endregion

        #region 6. GetAllSuffixes / EnumerateSuffixes Tests

        [Test]
        public void GetAllSuffixes_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.GetAllSuffixes(), Is.Empty);
        }

        [Test]
        public void EnumerateSuffixes_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.EnumerateSuffixes().ToList(), Is.Empty);
        }

        [Test]
        public void GetAllSuffixes_SingleChar_ReturnsSingleSuffix()
        {
            var st = SuffixTree.Build("a");
            var suffixes = st.GetAllSuffixes();

            Assert.That(suffixes, Has.Count.EqualTo(1));
            Assert.That(suffixes[0], Is.EqualTo("a"));
        }

        [Test]
        public void GetAllSuffixes_CountEqualsTextLength()
        {
            var text = "banana";
            var st = SuffixTree.Build(text);

            Assert.That(st.GetAllSuffixes().Count, Is.EqualTo(text.Length));
        }

        [Test]
        public void GetAllSuffixes_Banana_ReturnsCorrectSuffixes()
        {
            var st = SuffixTree.Build("banana");
            var suffixes = st.GetAllSuffixes();

            var expected = new[] { "a", "ana", "anana", "banana", "na", "nana" };
            Assert.That(suffixes, Is.EqualTo(expected));
        }

        [Test]
        public void GetAllSuffixes_AreSortedLexicographically()
        {
            var st = SuffixTree.Build("mississippi");
            var suffixes = st.GetAllSuffixes();

            Assert.That(suffixes, Is.EqualTo(suffixes.OrderBy(s => s).ToList()));
        }

        [Test]
        public void EnumerateSuffixes_MatchesGetAllSuffixes()
        {
            var st = SuffixTree.Build("banana");

            Assert.That(st.EnumerateSuffixes().ToList(), Is.EqualTo(st.GetAllSuffixes()));
        }

        [Test]
        public void EnumerateSuffixes_CanBreakEarly()
        {
            var st = SuffixTree.Build("abcdefghij");
            var firstThree = st.EnumerateSuffixes().Take(3).ToList();

            Assert.That(firstThree.Count, Is.EqualTo(3));
        }

        [Test]
        public void GetAllSuffixes_LargeString_NoStackOverflow()
        {
            var text = new string('a', 10_000);
            var st = SuffixTree.Build(text);

            IReadOnlyList<string> suffixes = null;
            Assert.DoesNotThrow(() => suffixes = st.GetAllSuffixes());
            Assert.That(suffixes.Count, Is.EqualTo(10_000));
        }

        #endregion

        #region 7. Statistics Properties Tests

        [Test]
        public void NodeCount_EmptyTree_IsPositive()
        {
            var st = SuffixTree.Build("");
            // Even empty tree has at least root node
            Assert.That(st.NodeCount, Is.GreaterThan(0));
        }

        [Test]
        public void NodeCount_GrowsWithTextLength()
        {
            var st1 = SuffixTree.Build("a");
            var st2 = SuffixTree.Build("ab");
            var st3 = SuffixTree.Build("abc");

            // More complex text = more nodes
            Assert.That(st2.NodeCount, Is.GreaterThanOrEqualTo(st1.NodeCount));
            Assert.That(st3.NodeCount, Is.GreaterThanOrEqualTo(st2.NodeCount));
        }

        [Test]
        public void LeafCount_EmptyTree_ReturnsZero()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.LeafCount, Is.EqualTo(0));
        }

        [Test]
        public void LeafCount_EqualsTextLengthPlusOne()
        {
            // LeafCount equals text length + 1 (for terminator suffix)
            var texts = new[] { "a", "ab", "banana", "mississippi" };
            foreach (var text in texts)
            {
                var st = SuffixTree.Build(text);
                Assert.That(st.LeafCount, Is.EqualTo(text.Length + 1),
                    $"LeafCount for '{text}' should be {text.Length + 1}");
            }
        }

        [Test]
        public void MaxDepth_EmptyTree_ReturnsZero()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.MaxDepth, Is.EqualTo(0));
        }

        [Test]
        public void MaxDepth_EqualsTextLength()
        {
            var st = SuffixTree.Build("abcdef");
            // Longest suffix is full string
            Assert.That(st.MaxDepth, Is.EqualTo(6));
        }

        [Test]
        public void MaxDepth_RepeatingString_EqualsTextLength()
        {
            var st = SuffixTree.Build("aaa");
            Assert.That(st.MaxDepth, Is.EqualTo(3));
        }

        #endregion

        #region 8. Unicode and Special Characters Tests

        [Test]
        public void Contains_CyrillicText_Works()
        {
            var st = SuffixTree.Build("привет мир");

            Assert.That(st.Contains("привет"), Is.True);
            Assert.That(st.Contains("мир"), Is.True);
            Assert.That(st.Contains("вет м"), Is.True);
        }

        [Test]
        public void Contains_ChineseCharacters_Works()
        {
            var st = SuffixTree.Build("你好世界");

            Assert.That(st.Contains("你好"), Is.True);
            Assert.That(st.Contains("世界"), Is.True);
            Assert.That(st.Contains("你好世界"), Is.True);
        }

        [Test]
        public void Contains_Emoji_Works()
        {
            var st = SuffixTree.Build("hello🌍world");

            Assert.That(st.Contains("hello"), Is.True);
            Assert.That(st.Contains("🌍"), Is.True);
            Assert.That(st.Contains("🌍world"), Is.True);
        }

        [Test]
        public void Contains_MixedCase_IsCaseSensitive()
        {
            var st = SuffixTree.Build("AbCdEf");

            Assert.That(st.Contains("AbCd"), Is.True);
            Assert.That(st.Contains("abcd"), Is.False);
            Assert.That(st.Contains("ABCD"), Is.False);
        }

        [Test]
        public void Contains_Whitespace_Works()
        {
            var st = SuffixTree.Build("a b c");

            Assert.That(st.Contains(" "), Is.True);
            Assert.That(st.Contains("a "), Is.True);
            Assert.That(st.Contains(" b"), Is.True);
        }

        [Test]
        public void Contains_TabsAndNewlines_Works()
        {
            var st = SuffixTree.Build("line1\nline2\tdata");

            Assert.That(st.Contains("\n"), Is.True);
            Assert.That(st.Contains("\t"), Is.True);
            Assert.That(st.Contains("line1\n"), Is.True);
        }

        [Test]
        public void Contains_SpecialCharacters_Works()
        {
            var st = SuffixTree.Build("a!@#$%^&*()b");

            Assert.That(st.Contains("!@#"), Is.True);
            Assert.That(st.Contains("$%^"), Is.True);
            Assert.That(st.Contains("&*()"), Is.True);
        }

        [Test]
        public void Contains_NumericString_Works()
        {
            var st = SuffixTree.Build("123123123");

            Assert.That(st.Contains("123"), Is.True);
            Assert.That(st.Contains("231"), Is.True);
            Assert.That(st.Contains("312"), Is.True);
        }

        [Test]
        public void PrintTree_UnicodeCharacters_DoesNotThrow()
        {
            var st = SuffixTree.Build("привет");
            Assert.DoesNotThrow(() => st.PrintTree());
        }

        #endregion

        #region 9. Performance and Large String Tests

        [Test]
        public void Build_LargeString_Completes()
        {
            var r = new Random(12345);
            var chars = new char[100_000];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = (char)('a' + r.Next(26));
            var text = new string(chars);

            SuffixTree st = null;
            Assert.DoesNotThrow(() => st = SuffixTree.Build(text));
            Assert.That(st, Is.Not.Null);
        }

        [Test]
        public void Contains_LargeString_Works()
        {
            var r = new Random(12345);
            var chars = new char[100_000];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = (char)('a' + r.Next(26));
            var text = new string(chars);
            var st = SuffixTree.Build(text);

            Assert.That(st.Contains(text.Substring(50000, 100)), Is.True);
            Assert.That(st.Contains("xyzxyzxyzxyz"), Is.False);
        }

        [Test]
        public void Build_DeepTree_NoStackOverflow()
        {
            // Worst case: all same characters creates very deep tree
            var text = new string('a', 10_000);
            var st = SuffixTree.Build(text);

            Assert.DoesNotThrow(() => st.Contains(new string('a', 5000)));
            Assert.DoesNotThrow(() => st.LongestRepeatedSubstring());
        }

        [Test]
        public void StressTest_AllSubstrings_SmallAlphabet()
        {
            var r = new Random(555);
            for (int trial = 0; trial < 10; trial++)
            {
                var chars = new char[30];
                for (int i = 0; i < chars.Length; i++)
                    chars[i] = (char)('a' + r.Next(3));

                var text = new string(chars);
                var st = SuffixTree.Build(text);

                // Verify ALL substrings
                for (int i = 0; i < text.Length; i++)
                {
                    for (int len = 1; len <= text.Length - i; len++)
                    {
                        var substr = text.Substring(i, len);
                        Assert.That(st.Contains(substr), Is.True,
                            $"Trial {trial}, string '{text}': Substring '{substr}' not found");
                    }
                }
            }
        }

        [Test]
        public void StressTest_RandomContains_Barrage()
        {
            var r = new Random(RANDOM_SEED);
            var text = MakeRandomString(r, 500);
            var st = SuffixTree.Build(text);

            for (int i = 0; i < 500; i++)
            {
                var pos = r.Next(0, text.Length - 10);
                var len = r.Next(1, 50);
                len = Math.Min(len, text.Length - pos);
                var substr = text.Substring(pos, len);

                Assert.That(st.Contains(substr), Is.True);
            }
        }

        #endregion

        #region 10. Thread Safety Tests

        [Test]
        public void ConcurrentReads_DoNotThrow()
        {
            var st = SuffixTree.Build("abracadabra mississippi banana");
            var patterns = new[] { "abra", "iss", "ana", "ab", "pp", "xyz" };

            var exceptions = new System.Collections.Concurrent.ConcurrentBag<Exception>();

            Parallel.For(0, 1000, i =>
            {
                try
                {
                    var pattern = patterns[i % patterns.Length];
                    _ = st.Contains(pattern);
                    _ = st.CountOccurrences(pattern);
                    _ = st.FindAllOccurrences(pattern);
                }
                catch (Exception ex)
                {
                    exceptions.Add(ex);
                }
            });

            Assert.That(exceptions, Is.Empty);
        }

        [Test]
        public void ConcurrentReads_ReturnConsistentResults()
        {
            var st = SuffixTree.Build("banana");
            var results = new System.Collections.Concurrent.ConcurrentBag<int>();

            Parallel.For(0, 100, _ =>
            {
                results.Add(st.CountOccurrences("ana"));
            });

            Assert.That(results.Distinct().Count(), Is.EqualTo(1));
            Assert.That(results.First(), Is.EqualTo(2));
        }

        [Test]
        public void ConcurrentReads_AllMethods_Stress()
        {
            var st = SuffixTree.Build("the quick brown fox jumps over the lazy dog");

            var tasks = new Task[10];
            for (int i = 0; i < 10; i++)
            {
                tasks[i] = Task.Run(() =>
                {
                    for (int j = 0; j < 50; j++)
                    {
                        st.Contains("the");
                        st.CountOccurrences("o");
                        st.FindAllOccurrences("e");
                        st.LongestRepeatedSubstring();
                        st.LongestCommonSubstring("fox dog");
                    }
                });
            }

            Task.WaitAll(tasks);
        }

        #endregion

        #region 11. Algorithm-Specific Tests (Ukkonen's Edge Cases)

        [Test]
        public void Algorithm_AbcAbxAbcd_TrickyCase()
        {
            // Known tricky case for Ukkonen's - requires correct suffix link handling
            var st = SuffixTree.Build("abcabxabcd");

            Assert.That(st.Contains("abcab"), Is.True);
            Assert.That(st.Contains("bcabx"), Is.True);
            Assert.That(st.Contains("xab"), Is.True);
            Assert.That(st.Contains("abcd"), Is.True);
        }

        [Test]
        public void Algorithm_SuffixLinks_Fibonacci()
        {
            // Fibonacci-like strings stress suffix links
            var st = SuffixTree.Build("abaababaabaab");

            Assert.That(st.Contains("abaab"), Is.True);
            Assert.That(st.Contains("baaba"), Is.True);
            Assert.That(st.Contains("abaababaabaab"), Is.True);
        }

        [Test]
        public void Algorithm_EdgeSplit_AtDifferentPositions()
        {
            var st = SuffixTree.Build("abcabdabeabf");

            // All share "ab" prefix but diverge
            Assert.That(st.Contains("abc"), Is.True);
            Assert.That(st.Contains("abd"), Is.True);
            Assert.That(st.Contains("abe"), Is.True);
            Assert.That(st.Contains("abf"), Is.True);
        }

        [Test]
        public void Algorithm_ImplicitToExplicit_Terminator()
        {
            // Without terminator, "a" would be implicit in "aa"
            var st = SuffixTree.Build("aa");

            Assert.That(st.Contains("aa"), Is.True);
            Assert.That(st.Contains("a"), Is.True);
        }

        [Test]
        public void Algorithm_SuffixLinks_RepeatingPatterns()
        {
            // Repeating patterns heavily exercise suffix links
            var st = SuffixTree.Build("abababababab");

            Assert.That(st.Contains("ababab"), Is.True);
            Assert.That(st.Contains("bababa"), Is.True);
            Assert.That(st.Contains("ababababab"), Is.True);
        }

        [Test]
        public void Algorithm_SuffixLinkChain_DeepNesting()
        {
            // Pattern that creates deep suffix link chains
            var st = SuffixTree.Build("aaabaaabaaab");

            Assert.That(st.Contains("aaabaaab"), Is.True);
            Assert.That(st.Contains("aabaaaba"), Is.True);
            Assert.That(st.Contains("abaaabaa"), Is.True);
            Assert.That(st.Contains("baaabaaab"), Is.True);
        }

        [Test]
        public void Algorithm_Rule1_ActivePointAtRoot()
        {
            // Test Rule 1: when activeNode is root, decrement activeLength
            var st = SuffixTree.Build("aabbaabb");

            Assert.That(st.Contains("aabb"), Is.True);
            Assert.That(st.Contains("abba"), Is.True);
            Assert.That(st.Contains("bbaa"), Is.True);
            Assert.That(st.Contains("baab"), Is.True);
        }

        [Test]
        public void Algorithm_Rule3_FollowSuffixLink()
        {
            // Test Rule 3: follow suffix link when not at root
            var st = SuffixTree.Build("xabxac");

            Assert.That(st.Contains("xab"), Is.True);
            Assert.That(st.Contains("xac"), Is.True);
            Assert.That(st.Contains("abxa"), Is.True);
            Assert.That(st.Contains("bxac"), Is.True);
        }

        [Test]
        public void Algorithm_PartialMatchAtEdgeBoundary()
        {
            // Test case where match ends exactly at edge boundary
            var st = SuffixTree.Build("aaaaab");

            Assert.That(st.Contains("aaaa"), Is.True);
            Assert.That(st.Contains("aaaaa"), Is.True);
            Assert.That(st.Contains("aaaaab"), Is.True);
            Assert.That(st.Contains("ba"), Is.False);
        }

        #endregion

        #region 12. Regression Tests (Specific Edge Cases)

        [Test]
        public void Regression_TerminatorCharacter_SuffixesAreExplicit()
        {
            // Without terminator, "ab" in "aab" would be implicit suffix
            // With terminator, all suffixes are explicit (end at leaf)
            var st = SuffixTree.Build("aab");

            Assert.That(st.Contains("aab"), Is.True);
            Assert.That(st.Contains("ab"), Is.True);
            Assert.That(st.Contains("b"), Is.True);
            Assert.That(st.Contains("a"), Is.True);
            Assert.That(st.Contains("aa"), Is.True);
        }

        [Test]
        public void Regression_MultipleBuildCallsAreIndependent()
        {
            // Verify that multiple trees are truly independent
            var st1 = SuffixTree.Build("hello");
            var st2 = SuffixTree.Build("world");
            var st3 = SuffixTree.Build("test");

            Assert.That(st1.Contains("hello"), Is.True);
            Assert.That(st1.Contains("world"), Is.False);
            Assert.That(st1.Contains("test"), Is.False);

            Assert.That(st2.Contains("world"), Is.True);
            Assert.That(st2.Contains("hello"), Is.False);

            Assert.That(st3.Contains("test"), Is.True);
            Assert.That(st3.Text, Is.EqualTo("test"));
        }

        [Test]
        public void Regression_MississippiClassicCase()
        {
            // Mississippi is the classic suffix tree test case
            var st = SuffixTree.Build("mississippi");

            // All suffixes
            Assert.That(st.Contains("mississippi"), Is.True);
            Assert.That(st.Contains("ississippi"), Is.True);
            Assert.That(st.Contains("i"), Is.True);

            // Key substrings
            Assert.That(st.Contains("issi"), Is.True);
            Assert.That(st.Contains("iss"), Is.True);
            Assert.That(st.Contains("ssi"), Is.True);
            Assert.That(st.Contains("sis"), Is.True);
            Assert.That(st.Contains("pp"), Is.True);
            Assert.That(st.Contains("ississi"), Is.True);
            Assert.That(st.Contains("sissipp"), Is.True);

            // Non-existent
            Assert.That(st.Contains("spa"), Is.False);
        }

        [Test]
        public void Regression_LargeAlphabet_AllCharsFound()
        {
            var alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
            var st = SuffixTree.Build(alphabet);

            foreach (char c in alphabet)
            {
                Assert.That(st.Contains(c.ToString()), Is.True, $"Character '{c}' not found");
            }
        }

        [Test]
        public void Regression_VeryLongRepeatingString()
        {
            var longString = new string('a', 10000) + "b";
            var st = SuffixTree.Build(longString);

            Assert.That(st.Contains(longString), Is.True);
            Assert.That(st.Contains("aaaaaaaaab"), Is.True);
            Assert.That(st.Contains("b"), Is.True);
        }

        [Test]
        public void Regression_EmptyStringToString()
        {
            var st = SuffixTree.Build("");
            var result = st.ToString();
            Assert.That(result, Does.Contain("empty"));
        }

        [Test]
        public void Regression_BinaryAlphabet_Stress()
        {
            var r = new Random(777);
            var chars = new char[500];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = r.Next(2) == 0 ? 'a' : 'b';

            var s = new string(chars);
            var st = SuffixTree.Build(s);

            for (int i = 0; i < s.Length - 20; i += 10)
            {
                for (int len = 1; len <= 20; len++)
                {
                    var substr = s.Substring(i, len);
                    Assert.That(st.Contains(substr), Is.True,
                        $"Substring '{substr}' at pos {i} not found");
                }
            }
        }

        [Test]
        public void Regression_SuffixLinkTraversal_LongRepeat()
        {
            var s = string.Concat(Enumerable.Repeat("abcdefgh", 100));
            var st = SuffixTree.Build(s);

            Assert.That(st.Contains("abcdefghabcdefgh"), Is.True);
            Assert.That(st.Contains("efghabcd"), Is.True);
            Assert.That(st.Contains("ghabcdefghabcdef"), Is.True);
            Assert.That(st.Contains(s), Is.True);
        }

        #endregion

        #region 13. Span API Tests

        [Test]
        public void Contains_Span_FromStringSlice_Works()
        {
            var st = SuffixTree.Build("hello world");
            var text = "xxxhelloxxx";

            Assert.That(st.Contains(text.AsSpan(3, 5)), Is.True);
        }

        [Test]
        public void Contains_Span_FromCharArray_Works()
        {
            var st = SuffixTree.Build("test string");
            var chars = new[] { 't', 'e', 's', 't' };

            Assert.That(st.Contains(chars.AsSpan()), Is.True);
        }

        [Test]
        public void FindAllOccurrences_Span_MatchesStringOverload()
        {
            var st = SuffixTree.Build("banana");
            var pattern = "ana";

            var stringResult = st.FindAllOccurrences(pattern);
            var spanResult = st.FindAllOccurrences(pattern.AsSpan());

            Assert.That(spanResult, Is.EquivalentTo(stringResult));
        }

        [Test]
        public void CountOccurrences_Span_MatchesStringOverload()
        {
            var st = SuffixTree.Build("mississippi");
            var pattern = "issi";

            Assert.That(st.CountOccurrences(pattern.AsSpan()), Is.EqualTo(st.CountOccurrences(pattern)));
        }

        #endregion

        #region Helper Methods

        private static string MakeRandomString(Random r, int len)
        {
            const string SET = "abcdefghijklmnopqrstuvwxyz";
            var res = new char[len];
            for (int i = 0; i < len; i++)
                res[i] = SET[r.Next(SET.Length)];
            return new string(res);
        }

        private static List<int> NaiveFindAll(string text, string pattern)
        {
            var result = new List<int>();
            int idx = 0;
            while ((idx = text.IndexOf(pattern, idx)) != -1)
            {
                result.Add(idx);
                idx++;
            }
            return result;
        }

        #endregion
    }
}
