using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Search
{
    /// <summary>
    /// Tests for Contains method (string and Span overloads).
    /// </summary>
    [TestFixture]
    public class ContainsTests
    {
        #region Input Validation

        [Test]
        public void Contains_NullPattern_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.Contains(null!));
        }

        [Test]
        public void Contains_EmptyPattern_ReturnsTrue()
        {
            // Specification: Empty string is a substring of any string (including empty string).
            // This follows standard string theory definition.
            var st = SuffixTree.Build("abc");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains(""), Is.True);
                Assert.That(st.Contains(ReadOnlySpan<char>.Empty), Is.True);
            });
        }

        [Test]
        public void Contains_EmptyTree_AnyPattern_ReturnsFalse()
        {
            // Specification: Empty tree contains no non-empty substrings.
            // Empty pattern in empty tree returns true (empty is substring of empty).
            var st = SuffixTree.Build("");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("a"), Is.False);
                Assert.That(st.Contains("abc"), Is.False);
                Assert.That(st.Contains(""), Is.True, "Empty pattern in empty tree should return true");
            });
        }

        #endregion

        #region Basic Matching

        [Test]
        public void Contains_FullString_ReturnsTrue()
        {
            var text = "teststring";
            var st = SuffixTree.Build(text);

            Assert.That(st.Contains(text), Is.True);
        }

        [Test]
        public void Contains_NonExistentPatterns_ReturnsFalse()
        {
            var st = SuffixTree.Build("abcdef");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("xyz"), Is.False, "Completely different pattern");
                Assert.That(st.Contains("abd"), Is.False, "Pattern with wrong character");
                Assert.That(st.Contains("abcdefg"), Is.False, "Pattern longer than text");
                Assert.That(st.Contains("aba"), Is.False, "Pattern not in text");
            });
        }

        #endregion

        #region All Substrings (Exhaustive)

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

        #endregion

        #region Overlapping Patterns

        [Test]
        public void Contains_OverlappingPatterns_Works()
        {
            var st = SuffixTree.Build("ababababab");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("abab"), Is.True);
                Assert.That(st.Contains("baba"), Is.True);
                Assert.That(st.Contains("ababa"), Is.True);
                Assert.That(st.Contains("ababababab"), Is.True);
                Assert.That(st.Contains("abababababa"), Is.False, "One char too long");
            });
        }

        [Test]
        public void Contains_RepeatingCharacter_Works()
        {
            var st = SuffixTree.Build("aaaaaaaaaa");

            for (int len = 1; len <= 10; len++)
            {
                Assert.That(st.Contains(new string('a', len)), Is.True, $"Should contain {len} a's");
            }
            Assert.That(st.Contains(new string('a', 11)), Is.False, "Should not contain 11 a's");
            Assert.That(st.Contains("b"), Is.False, "Should not contain 'b'");
        }

        #endregion

        #region Span Overload

        [Test]
        public void Contains_SpanOverload_MatchesStringOverload()
        {
            var st = SuffixTree.Build("hello world");

            var patterns = new[] { "hello", "world", "lo wo", "xyz", "" };

            foreach (var pattern in patterns)
            {
                Assert.That(
                    st.Contains(pattern.AsSpan()),
                    Is.EqualTo(st.Contains(pattern)),
                    $"Span and string overloads should match for '{pattern}'");
            }
        }

        [Test]
        public void Contains_SpanFromSlice_Works()
        {
            var st = SuffixTree.Build("hello world");
            var text = "xxxhelloxxx";

            Assert.That(st.Contains(text.AsSpan(3, 5)), Is.True);
        }

        [Test]
        public void Contains_SpanFromCharArray_Works()
        {
            var st = SuffixTree.Build("test string");
            var chars = new[] { 't', 'e', 's', 't' };

            Assert.That(st.Contains(chars.AsSpan()), Is.True);
        }

        #endregion

        #region Case Sensitivity

        [Test]
        public void Contains_IsCaseSensitive()
        {
            var st = SuffixTree.Build("AbCdEf");

            Assert.Multiple(() =>
            {
                Assert.That(st.Contains("AbCd"), Is.True);
                Assert.That(st.Contains("abcd"), Is.False, "Lowercase should not match");
                Assert.That(st.Contains("ABCD"), Is.False, "Uppercase should not match");
            });
        }

        #endregion
    }
}
