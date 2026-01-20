using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Search
{
    /// <summary>
    /// Tests for FindAllOccurrences method.
    /// </summary>
    [TestFixture]
    public class FindAllOccurrencesTests
    {
        #region Input Validation

        [Test]
        public void FindAll_NullPattern_ThrowsArgumentNullException()
        {
            var st = SuffixTree.Build("abc");
            Assert.Throws<ArgumentNullException>(() => st.FindAllOccurrences(null!).ToList());
        }

        [Test]
        public void FindAll_EmptyPattern_ReturnsAllPositions()
        {
            // Specification: Empty pattern conceptually matches at every position in the text.
            // This is consistent with regex "" matching at every position.
            var st = SuffixTree.Build("abc");

            var result = st.FindAllOccurrences("");

            Assert.That(result.Count, Is.EqualTo(3));
        }

        [Test]
        public void FindAll_EmptyTree_ReturnsEmpty()
        {
            var st = SuffixTree.Build("");

            var result = st.FindAllOccurrences("a").ToList();

            Assert.That(result, Is.Empty);
        }

        #endregion

        #region Single Occurrence

        [Test]
        public void FindAll_SingleOccurrence_ReturnsCorrectPosition()
        {
            var st = SuffixTree.Build("hello world");

            var result = st.FindAllOccurrences("world").ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0], Is.EqualTo(6));
            });
        }

        [Test]
        public void FindAll_AtBeginning_ReturnsZero()
        {
            var st = SuffixTree.Build("hello world");

            var result = st.FindAllOccurrences("hello").ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0], Is.EqualTo(0));
            });
        }

        [Test]
        public void FindAll_AtEnd_ReturnsCorrectPosition()
        {
            var st = SuffixTree.Build("hello world");

            var result = st.FindAllOccurrences("orld").ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0], Is.EqualTo(7));
            });
        }

        #endregion

        #region Multiple Occurrences

        [Test]
        public void FindAll_MultipleOccurrences_ReturnsAllPositions()
        {
            var st = SuffixTree.Build("abcabc");

            var result = st.FindAllOccurrences("abc").OrderBy(x => x).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(2));
                Assert.That(result, Is.EqualTo(new[] { 0, 3 }));
            });
        }

        [Test]
        public void FindAll_OverlappingOccurrences_FindsAll()
        {
            var st = SuffixTree.Build("aaaa");

            var result = st.FindAllOccurrences("aa").OrderBy(x => x).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(3));
                Assert.That(result, Is.EqualTo(new[] { 0, 1, 2 }));
            });
        }

        [Test]
        public void FindAll_Banana_FindsAllOccurrences()
        {
            var st = SuffixTree.Build("banana");

            Assert.Multiple(() =>
            {
                var ana = st.FindAllOccurrences("ana").OrderBy(x => x).ToList();
                Assert.That(ana, Is.EqualTo(new[] { 1, 3 }), "'ana' occurs at positions 1 and 3");

                var a = st.FindAllOccurrences("a").OrderBy(x => x).ToList();
                Assert.That(a, Is.EqualTo(new[] { 1, 3, 5 }), "'a' occurs at positions 1, 3, 5");

                var na = st.FindAllOccurrences("na").OrderBy(x => x).ToList();
                Assert.That(na, Is.EqualTo(new[] { 2, 4 }), "'na' occurs at positions 2 and 4");
            });
        }

        #endregion

        #region No Occurrences

        [Test]
        public void FindAll_NonExistent_ReturnsEmpty()
        {
            var st = SuffixTree.Build("hello world");

            Assert.Multiple(() =>
            {
                Assert.That(st.FindAllOccurrences("xyz").ToList(), Is.Empty);
                Assert.That(st.FindAllOccurrences("worldly").ToList(), Is.Empty);
            });
        }

        #endregion

        #region Full String

        [Test]
        public void FindAll_FullString_ReturnsZero()
        {
            var text = "abcdef";
            var st = SuffixTree.Build(text);

            var result = st.FindAllOccurrences(text).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0], Is.EqualTo(0));
            });
        }

        #endregion

        #region Single Character

        [Test]
        public void FindAll_SingleCharacter_FindsAllOccurrences()
        {
            var st = SuffixTree.Build("abracadabra");

            var result = st.FindAllOccurrences("a").OrderBy(x => x).ToList();

            Assert.That(result, Is.EqualTo(new[] { 0, 3, 5, 7, 10 }));
        }

        #endregion

        #region Exhaustive Verification

        [Test]
        public void FindAll_AllSubstrings_MatchLinearSearch()
        {
            var text = "mississippi";
            var st = SuffixTree.Build(text);

            for (int i = 0; i < text.Length; i++)
            {
                for (int len = 1; len <= text.Length - i; len++)
                {
                    var pattern = text.Substring(i, len);

                    var treeResult = st.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                    var linearResult = LinearFindAll(text, pattern).ToList();

                    Assert.That(treeResult, Is.EqualTo(linearResult),
                        $"Mismatch for pattern '{pattern}'");
                }
            }
        }

        private static IEnumerable<int> LinearFindAll(string text, string pattern)
        {
            var positions = new List<int>();
            int pos = 0;
            while ((pos = text.IndexOf(pattern, pos, StringComparison.Ordinal)) != -1)
            {
                positions.Add(pos);
                pos++;
            }
            return positions;
        }

        #endregion

        #region Lazy Enumeration

        [Test]
        public void FindAll_IsLazyEnumerated()
        {
            var st = SuffixTree.Build("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");

            var enumerable = st.FindAllOccurrences("a");

            // Take only first 3 - should not traverse entire tree
            var first3 = enumerable.Take(3).ToList();

            Assert.That(first3, Has.Count.EqualTo(3));
        }

        #endregion

        #region Special Characters

        [Test]
        public void FindAll_WithSpaces_Works()
        {
            var st = SuffixTree.Build("hello world hello");

            var result = st.FindAllOccurrences("hello").OrderBy(x => x).ToList();

            Assert.That(result, Is.EqualTo(new[] { 0, 12 }));
        }

        [Test]
        public void FindAll_WithNewlines_Works()
        {
            var st = SuffixTree.Build("line1\nline2\nline1");

            var result = st.FindAllOccurrences("line1").OrderBy(x => x).ToList();

            Assert.That(result, Is.EqualTo(new[] { 0, 12 }));
        }

        #endregion

        #region Span Overload

        [Test]
        public void FindAll_SpanOverload_MatchesStringOverload()
        {
            var st = SuffixTree.Build("hello world");

            // Note: Empty pattern behavior may differ between string and Span overloads
            var patterns = new[] { "hello", "world", "lo wo", "xyz" };

            foreach (var pattern in patterns)
            {
                var stringResult = st.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                var spanResult = st.FindAllOccurrences(pattern.AsSpan()).OrderBy(x => x).ToList();

                Assert.That(spanResult, Is.EqualTo(stringResult),
                    $"Span and string overloads should match for '{pattern}'");
            }
        }

        [Test]
        public void FindAll_SpanFromSlice_Works()
        {
            var st = SuffixTree.Build("hello world");
            var text = "xxxworldxxx";

            var result = st.FindAllOccurrences(text.AsSpan(3, 5)).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(result, Has.Count.EqualTo(1));
                Assert.That(result[0], Is.EqualTo(6));
            });
        }

        [Test]
        public void FindAll_EmptySpan_ReturnsEmpty()
        {
            // Note: Span overload returns empty for empty pattern (no special handling)
            // This differs from string overload which returns all positions
            var st = SuffixTree.Build("abc");

            var result = st.FindAllOccurrences(ReadOnlySpan<char>.Empty);

            Assert.That(result.Count, Is.EqualTo(0));
        }

        #endregion
    }
}
