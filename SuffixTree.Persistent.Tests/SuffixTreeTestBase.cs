using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests
{
    public abstract class SuffixTreeTestBase
    {
        protected abstract ISuffixTree CreateTree(string text);

        [Test]
        public void Contains_BasicPatterns_Works()
        {
            using (var st = CreateTree("banana") as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                Assert.Multiple(() =>
                {
                    Assert.That(tree.Contains("banana"), Is.True);
                    Assert.That(tree.Contains("ana"), Is.True);
                    Assert.That(tree.Contains("nan"), Is.True);
                    Assert.That(tree.Contains("ba"), Is.True);
                    Assert.That(tree.Contains("apple"), Is.False);
                    Assert.That(tree.Contains("bananax"), Is.False);
                });
            }
        }

        [Test]
        public void CountOccurrences_MultiplePatterns_Works()
        {
            using (var st = CreateTree("mississippi") as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                Assert.Multiple(() =>
                {
                    Assert.That(tree.CountOccurrences("i"), Is.EqualTo(4));
                    Assert.That(tree.CountOccurrences("s"), Is.EqualTo(4));
                    Assert.That(tree.CountOccurrences("ssi"), Is.EqualTo(2));
                    Assert.That(tree.CountOccurrences("p"), Is.EqualTo(2));
                    Assert.That(tree.CountOccurrences("mississippi"), Is.EqualTo(1));
                    Assert.That(tree.CountOccurrences("x"), Is.EqualTo(0));
                });
            }
        }

        [Test]
        public void FindAllOccurrences_ReturnsCorrectPositions()
        {
            using (var st = CreateTree("abracadabra") as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                var result = tree.FindAllOccurrences("abra").ToList();
                result.Sort();

                Assert.That(result, Is.EqualTo(new List<int> { 0, 7 }));
            }
        }

        [Test]
        public void LongestRepeatedSubstring_Works()
        {
            using (var st = CreateTree("banana") as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("ana"));
            }
        }

        [Test]
        public void LongestCommonSubstring_Works()
        {
            using (var st = CreateTree("abracadabra") as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                Assert.That(tree.LongestCommonSubstring("cadab"), Is.EqualTo("cadab"));
            }
        }

        [Test]
        public void EmptyTree_HandlesCorrectly()
        {
            using (var st = CreateTree("") as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                Assert.Multiple(() =>
                {
                    Assert.That(tree.Contains(""), Is.True);
                    Assert.That(tree.Contains("a"), Is.False);
                    Assert.That(tree.LeafCount, Is.EqualTo(0));
                });
            }
        }

        [Test]
        public void Unicode_SpecialCharacters_Works()
        {
            // Testing multi-byte character support (C# char is 2 bytes)
            string text = "🧬αβγ🧪$";
            using (var st = CreateTree(text) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                Assert.Multiple(() =>
                {
                    Assert.That(tree.Contains("🧬"), Is.True);
                    Assert.That(tree.Contains("αβγ"), Is.True);
                    Assert.That(tree.Contains("🧪"), Is.True);
                    Assert.That(tree.CountOccurrences("α"), Is.EqualTo(1));
                    Assert.That(tree.LongestCommonSubstring("βγ🧪"), Is.EqualTo("βγ🧪"));
                });
            }
        }

        [Test]
        public void AllSuffixes_AreContained()
        {
            string text = "abracadabra";
            using (var st = CreateTree(text) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                for (int i = 0; i < text.Length; i++)
                {
                    string suffix = text.Substring(i);
                    Assert.That(tree.Contains(suffix), Is.True, $"Suffix '{suffix}' not found");
                }
            }
        }
    }
}
