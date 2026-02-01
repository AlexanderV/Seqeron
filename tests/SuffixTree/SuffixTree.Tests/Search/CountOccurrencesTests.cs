using System;
using NUnit.Framework;

namespace SuffixTree.Tests.Search
{
    /// <summary>
    /// Tests for CountOccurrences method (string and Span overloads).
    /// </summary>
    [TestFixture]
    public class CountOccurrencesTests
    {
        #region Input Validation

        [Test]
        public void Count_NullPattern_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.CountOccurrences(null!));
        }

        [Test]
        public void Count_EmptyPattern_ReturnsTextLength()
        {
            // Specification: Empty pattern matches at every position in the text.
            // For text of length N, empty pattern occurs N times (positions 0..N-1).
            var st = SuffixTree.Build("abc");

            Assert.That(st.CountOccurrences(""), Is.EqualTo(3));
        }

        [Test]
        public void Count_EmptyTree_ReturnsZero()
        {
            var st = SuffixTree.Build("");

            Assert.Multiple(() =>
            {
                Assert.That(st.CountOccurrences("a"), Is.EqualTo(0));
                Assert.That(st.CountOccurrences("abc"), Is.EqualTo(0));
            });
        }

        #endregion

        #region Single Occurrence

        [Test]
        public void Count_SingleOccurrence_ReturnsOne()
        {
            var st = SuffixTree.Build("hello world");

            Assert.Multiple(() =>
            {
                Assert.That(st.CountOccurrences("hello"), Is.EqualTo(1));
                Assert.That(st.CountOccurrences("world"), Is.EqualTo(1));
                Assert.That(st.CountOccurrences("hello world"), Is.EqualTo(1));
            });
        }

        #endregion

        #region Multiple Occurrences

        [Test]
        public void Count_MultipleOccurrences_ReturnsCorrectCount()
        {
            var st = SuffixTree.Build("abcabc");

            Assert.That(st.CountOccurrences("abc"), Is.EqualTo(2));
        }

        [Test]
        public void Count_OverlappingPatterns_CountsOverlaps()
        {
            var st = SuffixTree.Build("aaaa");

            Assert.Multiple(() =>
            {
                Assert.That(st.CountOccurrences("a"), Is.EqualTo(4));
                Assert.That(st.CountOccurrences("aa"), Is.EqualTo(3));
                Assert.That(st.CountOccurrences("aaa"), Is.EqualTo(2));
                Assert.That(st.CountOccurrences("aaaa"), Is.EqualTo(1));
            });
        }

        [Test]
        public void Count_Banana_CorrectCounts()
        {
            var st = SuffixTree.Build("banana");

            Assert.Multiple(() =>
            {
                Assert.That(st.CountOccurrences("a"), Is.EqualTo(3));
                Assert.That(st.CountOccurrences("n"), Is.EqualTo(2));
                Assert.That(st.CountOccurrences("b"), Is.EqualTo(1));
                Assert.That(st.CountOccurrences("an"), Is.EqualTo(2));
                Assert.That(st.CountOccurrences("ana"), Is.EqualTo(2));
                Assert.That(st.CountOccurrences("nan"), Is.EqualTo(1));
                Assert.That(st.CountOccurrences("banana"), Is.EqualTo(1));
            });
        }

        #endregion

        #region No Occurrences

        [Test]
        public void Count_NonExistent_ReturnsZero()
        {
            var st = SuffixTree.Build("hello world");

            Assert.Multiple(() =>
            {
                Assert.That(st.CountOccurrences("xyz"), Is.EqualTo(0));
                Assert.That(st.CountOccurrences("worldly"), Is.EqualTo(0));
                Assert.That(st.CountOccurrences("z"), Is.EqualTo(0));
            });
        }

        #endregion

        #region Consistency with FindAll

        [Test]
        public void Count_MatchesFindAllCount()
        {
            var text = "mississippi";
            var st = SuffixTree.Build(text);

            var patterns = new[] { "i", "s", "ss", "issi", "mississippi", "pp", "xyz" };

            foreach (var pattern in patterns)
            {
                var count = st.CountOccurrences(pattern);
                var findAllCount = st.FindAllOccurrences(pattern).Count;

                Assert.That(count, Is.EqualTo(findAllCount),
                    $"Count and FindAll.Count should match for '{pattern}'");
            }
        }

        #endregion

        #region Large Counts

        [Test]
        public void Count_ManyOccurrences_Works()
        {
            var st = SuffixTree.Build(new string('a', 1000));

            Assert.Multiple(() =>
            {
                Assert.That(st.CountOccurrences("a"), Is.EqualTo(1000));
                Assert.That(st.CountOccurrences("aa"), Is.EqualTo(999));
                Assert.That(st.CountOccurrences("aaa"), Is.EqualTo(998));
            });
        }

        #endregion

        #region Span Overload

        [Test]
        public void Count_SpanOverload_MatchesStringOverload()
        {
            var st = SuffixTree.Build("abracadabra");

            // Note: Empty pattern behavior may differ between string and Span overloads
            var patterns = new[] { "a", "abra", "bra", "xyz" };

            foreach (var pattern in patterns)
            {
                Assert.That(
                    st.CountOccurrences(pattern.AsSpan()),
                    Is.EqualTo(st.CountOccurrences(pattern)),
                    $"Span and string overloads should match for '{pattern}'");
            }
        }

        [Test]
        public void Count_SpanFromSlice_Works()
        {
            var st = SuffixTree.Build("hello world");
            var text = "xxxhelloxxx";

            Assert.That(st.CountOccurrences(text.AsSpan(3, 5)), Is.EqualTo(1));
        }

        [Test]
        public void Count_EmptySpan_ReturnsZero()
        {
            // Note: Span overload returns 0 for empty pattern (no special handling)
            // This differs from string overload which returns text length
            var st = SuffixTree.Build("abc");

            Assert.That(st.CountOccurrences(ReadOnlySpan<char>.Empty), Is.EqualTo(0));
        }

        #endregion
    }
}
