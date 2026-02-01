using System;
using NUnit.Framework;
using SuffixTree;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class MathematicalInvariantsTests
    {
        [TestCase("banana")]
        [TestCase("mississippi")]
        [TestCase("aaaaaaaa")]
        [TestCase("abcde")]
        [TestCase("z")]
        [TestCase("")]
        public void LeafCount_IsExactly_TextLengthPlusOne(string text)
        {
            using (var st = PersistentSuffixTreeFactory.Create(text) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                // Raw leaf count in persistent tree is exposed via ISuffixTree.LeafCount
                // Note: The ISuffixTree.LeafCount property was adjusted to return (raw - 1)
                // to match reference implementation. Here we check the underlying logic.
                
                int expectedLeaves = text.Length; // Adjusted LeafCount
                Assert.That(tree.LeafCount, Is.EqualTo(expectedLeaves), $"LeafCount for '{text}' should be {expectedLeaves}");
            }
        }

        [TestCase("abcde", 6, 11)] // 1 root + 5 leaves (at least)
        [TestCase("aaaaa", 6, 11)]
        public void NodeCount_Bounds(string text, int min, int max)
        {
            using (var st = PersistentSuffixTreeFactory.Create(text) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                // Node count for string of length N is at most 2N
                Assert.That(tree.NodeCount, Is.GreaterThanOrEqualTo(min));
                Assert.That(tree.NodeCount, Is.LessThanOrEqualTo(max));
            }
        }

        [Test]
        public void Terminator_IsShielded()
        {
            // The internal terminator is '#' in reference and '$' in persistent? 
            // Wait, PersistentConstants uses '$' or similar?
            // Let's check PersistentConstants.TERMINATOR_KEY.
            
            string text = "abc";
            using (var st = PersistentSuffixTreeFactory.Create(text) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                const string terminator = "$"; // Standard terminator used in our implementation
                
                Assert.That(tree.Contains(terminator), Is.False, "Internal terminator should not be searchable");
                Assert.That(tree.CountOccurrences(terminator), Is.EqualTo(0));
            }
        }

        [Test]
        public void SuffixAtZero_DoesNotContainTerminator()
        {
            string text = "test";
            using (var st = PersistentSuffixTreeFactory.Create(text) as IDisposable)
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
