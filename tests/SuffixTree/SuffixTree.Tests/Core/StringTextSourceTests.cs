using System;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Core
{
    /// <summary>
    /// Tests for StringTextSource — the ITextSource implementation that wraps a string.
    /// </summary>
    [TestFixture]
    public class StringTextSourceTests
    {
        #region Constructor Guards

        [Test]
        public void Constructor_NullString_ThrowsArgumentNullException()
        {
            Assert.Throws<ArgumentNullException>(() => new StringTextSource(null!));
        }

        #endregion

        #region Properties

        [Test]
        public void Length_ReturnsStringLength()
        {
            var source = new StringTextSource("hello");
            Assert.That(source.Length, Is.EqualTo(5));
        }

        [Test]
        public void Length_EmptyString_ReturnsZero()
        {
            var source = new StringTextSource("");
            Assert.That(source.Length, Is.EqualTo(0));
        }

        [Test]
        public void Indexer_ReturnsCorrectChar()
        {
            var source = new StringTextSource("abc");
            Assert.Multiple(() =>
            {
                Assert.That(source[0], Is.EqualTo('a'));
                Assert.That(source[1], Is.EqualTo('b'));
                Assert.That(source[2], Is.EqualTo('c'));
            });
        }

        [Test]
        public void Value_ReturnsOriginalString()
        {
            var source = new StringTextSource("test");
            Assert.That(source.Value, Is.EqualTo("test"));
        }

        #endregion

        #region Substring

        [Test]
        public void Substring_ReturnsCorrectSubstring()
        {
            var source = new StringTextSource("hello world");
            Assert.That(source.Substring(6, 5), Is.EqualTo("world"));
        }

        [Test]
        public void Substring_FromStart_Works()
        {
            var source = new StringTextSource("abcdef");
            Assert.That(source.Substring(0, 3), Is.EqualTo("abc"));
        }

        [Test]
        public void Substring_ZeroLength_ReturnsEmpty()
        {
            var source = new StringTextSource("abc");
            Assert.That(source.Substring(1, 0), Is.EqualTo(""));
        }

        #endregion

        #region Slice

        [Test]
        public void Slice_ReturnsCorrectSpan()
        {
            var source = new StringTextSource("hello world");
            var slice = source.Slice(6, 5);
            Assert.That(slice.ToString(), Is.EqualTo("world"));
        }

        [Test]
        public void Slice_ZeroLength_ReturnsEmptySpan()
        {
            var source = new StringTextSource("abc");
            var slice = source.Slice(1, 0);
            Assert.That(slice.Length, Is.EqualTo(0));
        }

        [Test]
        public void Slice_FullString_ReturnsAllChars()
        {
            var source = new StringTextSource("test");
            var slice = source.Slice(0, 4);
            Assert.That(slice.ToString(), Is.EqualTo("test"));
        }

        #endregion

        #region Enumeration

        [Test]
        public void GetEnumerator_YieldsAllChars()
        {
            var source = new StringTextSource("abc");
            var chars = source.ToList();
            Assert.That(chars, Is.EqualTo(new[] { 'a', 'b', 'c' }));
        }

        [Test]
        public void GetEnumerator_EmptyString_YieldsNothing()
        {
            var source = new StringTextSource("");
            var chars = source.ToList();
            Assert.That(chars, Is.Empty);
        }

        #endregion

        #region Implicit Conversion

        [Test]
        public void ImplicitConversion_StringToTextSource_Works()
        {
            StringTextSource source = "hello";
            Assert.Multiple(() =>
            {
                Assert.That(source.Length, Is.EqualTo(5));
                Assert.That(source.Value, Is.EqualTo("hello"));
            });
        }

        [Test]
        public void ImplicitConversion_UsedInBuild_Works()
        {
            // SuffixTree.Build(ITextSource) works with implicit conversion
            var tree = SuffixTree.Build((ITextSource)new StringTextSource("banana"));

            Assert.Multiple(() =>
            {
                Assert.That(tree.LeafCount, Is.EqualTo(6));
                Assert.That(tree.Contains("ana"), Is.True);
            });
        }

        #endregion

        #region ToString

        [Test]
        public void ToString_ReturnsUnderlyingString()
        {
            var source = new StringTextSource("test");
            Assert.That(source.ToString(), Is.EqualTo("test"));
        }

        #endregion
    }
}
