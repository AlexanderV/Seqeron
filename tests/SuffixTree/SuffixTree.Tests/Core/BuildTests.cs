#nullable enable
using System;
using System.Collections;
using System.Collections.Generic;
using NUnit.Framework;

namespace SuffixTree.Tests.Core
{
    /// <summary>
    /// Tests for SuffixTree construction and basic properties.
    /// </summary>
    [TestFixture]
    [Category("Core")]
    public class BuildTests
    {
        #region Null and Empty Input

        [Test]
        public void Build_NullString_ThrowsArgumentNullException()
        {
            string nullString = null!;
            Assert.Throws<ArgumentNullException>(() => SuffixTree.Build(nullString));
        }

        [Test]
        public void Build_EmptyString_CreatesValidTree()
        {
            var st = SuffixTree.Build("");

            Assert.Multiple(() =>
            {
                Assert.That(st, Is.Not.Null);
                Assert.That(st.Text.ToString(), Is.EqualTo(""));
                Assert.That(st.NodeCount, Is.EqualTo(1), "Empty tree should have only root node");
                Assert.That(st.LeafCount, Is.EqualTo(0));
                Assert.That(st.MaxDepth, Is.EqualTo(0));
            });
        }

        #endregion

        #region Single Character

        [Test]
        public void Build_SingleCharacter_CreatesValidTree()
        {
            var st = SuffixTree.Build("a");

            Assert.Multiple(() =>
            {
                Assert.That(st.Text.ToString(), Is.EqualTo("a"));
                Assert.That(st.LeafCount, Is.EqualTo(1));
                Assert.That(st.MaxDepth, Is.EqualTo(1));
            });
        }

        #endregion

        #region Classic Test Cases

        [TestCase("banana", Description = "Classic LRS example")]
        [TestCase("mississippi", Description = "Classic suffix tree example")]
        [TestCase("abracadabra", Description = "Repeating patterns")]
        [TestCase("abcabxabcd", Description = "Ukkonen's tricky case")]
        public void Build_ClassicTestCases_CreatesValidTree(string text)
        {
            var st = SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(st.Text.ToString(), Is.EqualTo(text));
                Assert.That(st.LeafCount, Is.EqualTo(text.Length), "LeafCount should equal text length");
                Assert.That(st.MaxDepth, Is.EqualTo(text.Length), "MaxDepth should equal text length");
            });
        }

        #endregion

        #region Special Patterns

        [Test]
        public void Build_AllSameCharacters_CreatesValidTree()
        {
            var st = SuffixTree.Build("aaaaaaaaaa");

            Assert.Multiple(() =>
            {
                Assert.That(st.LeafCount, Is.EqualTo(10));
                Assert.That(st.MaxDepth, Is.EqualTo(10));
            });
        }

        [Test]
        public void Build_AlternatingPattern_CreatesValidTree()
        {
            var st = SuffixTree.Build("abababab");

            Assert.Multiple(() =>
            {
                Assert.That(st.LeafCount, Is.EqualTo(8));
                Assert.That(st.MaxDepth, Is.EqualTo(8));
                // Alternating "ab" pattern: "ab" occurs 4 times
                Assert.That(st.CountOccurrences("ab"), Is.EqualTo(4));
                Assert.That(st.LongestRepeatedSubstring(), Has.Length.EqualTo(6)); // "ababab"
            });
        }

        [Test]
        public void Build_AllUniqueCharacters_CreatesValidTree()
        {
            var st = SuffixTree.Build("abcdefghij");

            Assert.Multiple(() =>
            {
                Assert.That(st.LeafCount, Is.EqualTo(10));
                Assert.That(st.MaxDepth, Is.EqualTo(10));
                // All unique → no repeated substring
                Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo(string.Empty));
                // Each char occurs exactly once
                Assert.That(st.CountOccurrences("a"), Is.EqualTo(1));
            });
        }

        #endregion

        #region Text Property

        [Test]
        public void Text_PreservesOriginalContentExactly()
        {
            var inputs = new[] {
                "test string with spaces",
                "  leading and trailing  ",
                "line1\nline2\rline3",
                "\t\ttabs"
            };

            foreach (var input in inputs)
            {
                var st = SuffixTree.Build(input);
                Assert.That(st.Text.ToString(), Is.EqualTo(input),
                    $"Text property must preserve exact content for: [{input}]");
            }
        }

        #endregion

        #region Multiple Builds Are Independent

        [Test]
        public void Build_MultipleTrees_AreIndependent()
        {
            var st1 = SuffixTree.Build("hello");
            var st2 = SuffixTree.Build("world");
            var st3 = SuffixTree.Build("test");

            Assert.Multiple(() =>
            {
                Assert.That(st1.Text.ToString(), Is.EqualTo("hello"));
                Assert.That(st2.Text.ToString(), Is.EqualTo("world"));
                Assert.That(st3.Text.ToString(), Is.EqualTo("test"));

                Assert.That(st1.Contains("world"), Is.False);
                Assert.That(st2.Contains("hello"), Is.False);
            });
        }

        #endregion

        #region TryBuild API

        [Test]
        public void TryBuild_NullString_ReturnsFalse()
        {
            bool result = SuffixTree.TryBuild((string?)null, out var tree);

            Assert.Multiple(() =>
            {
                Assert.That(result, Is.False);
                Assert.That(tree, Is.Null);
            });
        }

        [Test]
        public void TryBuild_ValidString_ReturnsTrueWithTree()
        {
            bool result = SuffixTree.TryBuild("banana", out var tree);

            Assert.Multiple(() =>
            {
                Assert.That(result, Is.True);
                Assert.That(tree, Is.Not.Null);
                Assert.That(tree!.Text.ToString(), Is.EqualTo("banana"));
                Assert.That(tree.Contains("nan"), Is.True);
            });
        }

        [Test]
        public void TryBuild_EmptyString_ReturnsTrueWithEmptyTree()
        {
            bool result = SuffixTree.TryBuild("", out var tree);

            Assert.Multiple(() =>
            {
                Assert.That(result, Is.True);
                Assert.That(tree, Is.Not.Null);
                Assert.That(tree!.IsEmpty, Is.True);
            });
        }

        #endregion

        #region IsEmpty Property

        [Test]
        public void IsEmpty_CorrelatesWithTextLength()
        {
            Assert.Multiple(() =>
            {
                Assert.That(SuffixTree.Build("").IsEmpty, Is.True);
                Assert.That(SuffixTree.Build("a").IsEmpty, Is.False);
                Assert.That(SuffixTree.Empty.IsEmpty, Is.True);
            });
        }

        [Test]
        public void Empty_StaticProperty_ReturnsEmptyTree()
        {
            var empty = SuffixTree.Empty;

            Assert.Multiple(() =>
            {
                Assert.That(empty, Is.Not.Null);
                Assert.That(empty.IsEmpty, Is.True);
                Assert.That(empty.Text.ToString(), Is.EqualTo(""));
                Assert.That(empty.LeafCount, Is.EqualTo(0));
            });
        }

        #endregion

        #region Memory/Span Build API

        [Test]
        public void Build_ReadOnlyMemory_CreatesValidTree()
        {
            ReadOnlyMemory<char> memory = "banana".AsMemory();
            var tree = SuffixTree.Build(memory);

            Assert.Multiple(() =>
            {
                Assert.That(tree, Is.Not.Null);
                Assert.That(tree.Text.ToString(), Is.EqualTo("banana"));
                Assert.That(tree.Contains("nan"), Is.True);
            });
        }

        [Test]
        public void Build_ReadOnlySpan_CreatesValidTree()
        {
            ReadOnlySpan<char> span = "banana".AsSpan();
            var tree = SuffixTree.Build(span);

            Assert.Multiple(() =>
            {
                Assert.That(tree, Is.Not.Null);
                Assert.That(tree.Text.ToString(), Is.EqualTo("banana"));
                Assert.That(tree.Contains("nan"), Is.True);
            });
        }

        #endregion

        #region Empty Singleton

        [Test]
        public void Empty_ReturnsSameInstance()
        {
            var empty1 = SuffixTree.Empty;
            var empty2 = SuffixTree.Empty;

            Assert.That(empty1, Is.SameAs(empty2));
        }

        #endregion

        #region TryBuild(ITextSource) Overload

        [Test]
        public void TryBuild_NullTextSource_ReturnsFalse()
        {
            bool result = SuffixTree.TryBuild((ITextSource?)null, out var tree);

            Assert.Multiple(() =>
            {
                Assert.That(result, Is.False);
                Assert.That(tree, Is.Null);
            });
        }

        [Test]
        public void TryBuild_ValidTextSource_ReturnsTrue()
        {
            var source = new ArrayTextSource("banana");
            bool result = SuffixTree.TryBuild(source, out var tree);

            Assert.Multiple(() =>
            {
                Assert.That(result, Is.True);
                Assert.That(tree, Is.Not.Null);
                Assert.That(tree!.Contains("nan"), Is.True);
            });
        }

        #endregion

        #region Custom ITextSource (non-StringTextSource) Build Path

        [Test]
        public void Build_CustomTextSource_ConstructsCorrectTree()
        {
            // Using a non-StringTextSource exercises the fallback code paths:
            // - Constructor foreach loop instead of raw string indexing
            // - GetSymbolAt fallback (_text[index] instead of _rawString[index])
            var source = new ArrayTextSource("banana");
            SuffixTree.TryBuild(source, out var tree);

            Assert.Multiple(() =>
            {
                Assert.That(tree!.Text.Length, Is.EqualTo(6));
                Assert.That(tree.LeafCount, Is.EqualTo(6), "Leaf count = text length");
                Assert.That(tree.Contains("banana"), Is.True);
                Assert.That(tree.Contains("ana"), Is.True);
                Assert.That(tree.Contains("xyz"), Is.False);
            });
        }

        [Test]
        public void Build_CustomTextSource_LRSMemoryFallback()
        {
            // LongestRepeatedSubstringMemory() has a fallback path when _rawString is null.
            // A custom ITextSource (not StringTextSource) triggers this path.
            var source = new ArrayTextSource("banana");
            SuffixTree.TryBuild(source, out var tree);

            var lrsMemory = tree!.LongestRepeatedSubstringMemory();
            Assert.That(lrsMemory.ToString(), Is.EqualTo("ana"),
                "LRS memory fallback should return same result as LRS string path");
        }

        [Test]
        public void Build_CustomTextSource_SearchAndCount()
        {
            var source = new ArrayTextSource("abcabc");
            SuffixTree.TryBuild(source, out var tree);

            Assert.Multiple(() =>
            {
                Assert.That(tree!.CountOccurrences("abc"), Is.EqualTo(2));
                Assert.That(tree.FindAllOccurrences("abc"), Is.EquivalentTo(new[] { 0, 3 }));
            });
        }

        /// <summary>
        /// A simple ITextSource backed by a char array (not StringTextSource).
        /// Forces the suffix tree to use fallback ITextSource code paths.
        /// </summary>
        private sealed class ArrayTextSource : ITextSource
        {
            private readonly char[] _data;

            public ArrayTextSource(string text) => _data = text.ToCharArray();

            public int Length => _data.Length;
            public char this[int index] => _data[index];
            public string Substring(int start, int length) => new string(_data, start, length);
            public ReadOnlySpan<char> Slice(int start, int length) => _data.AsSpan(start, length);

            public IEnumerator<char> GetEnumerator()
            {
                foreach (var c in _data) yield return c;
            }

            IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
        }

        #endregion
    }
}

