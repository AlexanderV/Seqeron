using System;
using System.Linq;
using NUnit.Framework;
using SuffixTree;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests.Validation
{
    [TestFixture]
    [Category("Validation")]
    public class MathematicalInvariantsTests
    {
        [TestCase("banana")]
        [TestCase("mississippi")]
        [TestCase("aaaaaaaa")]
        [TestCase("abcde")]
        [TestCase("z")]
        [TestCase("")]
        public void LeafCount_IsExactly_TextLength(string text)
        {
            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                int expectedLeaves = text.Length;
                Assert.That(tree.LeafCount, Is.EqualTo(expectedLeaves), $"LeafCount for '{text}' should be {expectedLeaves}");
            }
        }

        [TestCase("abcde", 6, 11)] // 1 root + 5 leaves (at least)
        [TestCase("aaaaa", 6, 11)]
        public void NodeCount_Bounds(string text, int min, int max)
        {
            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                // Node count for string of length N is at most 2N
                Assert.That(tree.NodeCount, Is.GreaterThanOrEqualTo(min));
                Assert.That(tree.NodeCount, Is.LessThanOrEqualTo(max));
            }
        }

        [Test]
        public void SyntheticTerminatorCharacter_IsNotImplicitlyAppended()
        {
            string text = "abc";
            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                string extended = text + "$";
                Assert.That(tree.Contains(extended), Is.False, "Tree must not expose synthetic terminator as user-visible symbol");
                Assert.That(tree.CountOccurrences(extended), Is.EqualTo(0));
            }
        }

        [Test]
        public void SuffixAtZero_DoesNotContainTerminator()
        {
            string text = "test";
            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                var results = tree.FindAllOccurrences("test").ToList();
                Assert.That(results.Count, Is.EqualTo(1));
                Assert.That(results[0], Is.EqualTo(0));

                // Ensure length is exactly what was passed
                Assert.That(tree.Text.Length, Is.EqualTo(text.Length));
            }
        }
    }
}
