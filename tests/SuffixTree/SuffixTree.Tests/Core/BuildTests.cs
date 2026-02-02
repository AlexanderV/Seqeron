using System;
using NUnit.Framework;

namespace SuffixTree.Tests.Core
{
    /// <summary>
    /// Tests for SuffixTree construction and basic properties.
    /// </summary>
    [TestFixture]
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

            Assert.That(st.LeafCount, Is.EqualTo(8));
        }

        [Test]
        public void Build_AllUniqueCharacters_CreatesValidTree()
        {
            var st = SuffixTree.Build("abcdefghij");

            Assert.That(st.LeafCount, Is.EqualTo(10));
        }

        #endregion

        #region Text Property

        [Test]
        public void Text_ReturnsOriginalString()
        {
            var original = "test string with spaces";
            var st = SuffixTree.Build(original);

            Assert.That(st.Text.ToString(), Is.EqualTo(original));
        }

        [Test]
        public void Text_PreservesWhitespace()
        {
            var text = "  leading and trailing  ";
            var st = SuffixTree.Build(text);

            Assert.That(st.Text.ToString(), Is.EqualTo(text));
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
        public void IsEmpty_EmptyString_ReturnsTrue()
        {
            var st = SuffixTree.Build("");
            Assert.That(st.IsEmpty, Is.True);
        }

        [Test]
        public void IsEmpty_NonEmptyString_ReturnsFalse()
        {
            var st = SuffixTree.Build("a");
            Assert.That(st.IsEmpty, Is.False);
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
    }
}
