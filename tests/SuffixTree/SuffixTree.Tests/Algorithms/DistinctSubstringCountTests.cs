using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Distinct-substring counting (per length and total) verified against brute-force enumeration.
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class DistinctSubstringCountTests
    {
        private static long[] BruteForceByLength(string text, int maxLength)
        {
            int m = Math.Min(maxLength, text.Length);
            var counts = new long[m + 1];
            counts[0] = 1;
            for (int len = 1; len <= m; len++)
            {
                var seen = new HashSet<string>(StringComparer.Ordinal);
                for (int i = 0; i + len <= text.Length; i++)
                    seen.Add(text.Substring(i, len));
                counts[len] = seen.Count;
            }
            return counts;
        }

        [Test]
        public void Banana_LiteratureValues()
        {
            var tree = SuffixTree.Build("banana");

            // b,a,n | ba,an,na | ban,ana,nan | bana,anan,nana | banan,anana | banana
            Assert.That(tree.CountDistinctSubstringsByLength(6), Is.EqualTo(new long[] { 1, 3, 3, 3, 3, 2, 1 }));
            Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(15));
        }

        [TestCase("abcd", 10)]      // all n(n+1)/2 substrings distinct
        [TestCase("aaaa", 4)]       // one substring per length
        [TestCase("a", 1)]
        [TestCase("ATTTGGATT", 35)]
        public void CountDistinctSubstrings_KnownTotals(string text, long expected)
        {
            Assert.That(SuffixTree.Build(text).CountDistinctSubstrings(), Is.EqualTo(expected));
        }

        [Test]
        public void EmptyText_HasNoSubstrings()
        {
            var tree = SuffixTree.Build(string.Empty);
            Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(0));
            Assert.That(tree.CountDistinctSubstringsByLength(5), Is.EqualTo(new long[] { 1 }));
        }

        [Test]
        public void MaxLengthBelowOne_Throws()
        {
            var tree = SuffixTree.Build("acgt");
            Assert.Throws<ArgumentOutOfRangeException>(() => tree.CountDistinctSubstringsByLength(0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeAlgorithms.CountDistinctSubstringsByLength(tree, -1));
        }

        [Test]
        public void MaxLengthBeyondText_IsClampedToTextLength()
        {
            var tree = SuffixTree.Build("acgtac");
            Assert.That(tree.CountDistinctSubstringsByLength(int.MaxValue), Is.EqualTo(BruteForceByLength("acgtac", 6)));
        }

        [Test]
        public void RandomTexts_MatchBruteForce_AllPaths([Values(0, 1, 2, 3, 4)] int seed, [Values(2, 4, 26)] int alphabet)
        {
            var rng = new Random(seed * 31 + alphabet);
            for (int round = 0; round < 20; round++)
            {
                int n = rng.Next(1, 120);
                string text = new string(Enumerable.Range(0, n).Select(_ => (char)('a' + rng.Next(alphabet))).ToArray());
                var tree = SuffixTree.Build(text);
                ISuffixTree asInterface = tree;

                foreach (int m in new[] { 1, 3, 7, n, n + 5 })
                {
                    var expected = BruteForceByLength(text, m);
                    Assert.That(tree.CountDistinctSubstringsByLength(m), Is.EqualTo(expected), $"direct walk, text={text}, m={m}");
                    Assert.That(SuffixTreeAlgorithms.CountDistinctSubstringsByLength(asInterface, m), Is.EqualTo(expected),
                        $"visitor path, text={text}, m={m}");
                }

                long total = BruteForceByLength(text, n).Skip(1).Sum();
                Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(total), text);
                Assert.That(SuffixTreeAlgorithms.CountDistinctSubstrings(asInterface), Is.EqualTo(total), text);
            }
        }
    }
}
