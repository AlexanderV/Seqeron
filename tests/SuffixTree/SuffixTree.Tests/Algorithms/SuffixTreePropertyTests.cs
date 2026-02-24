using System;
using System.Collections.Generic;
using System.Linq;
using CsCheck;
using NUnit.Framework;

namespace SuffixTree.Tests.Algorithms
{
    /// <summary>
    /// Focused CsCheck properties that add brute-force oracle coverage.
    /// Broad invariants are covered in Core/InvariantTests and Properties/PropertyBasedTests.
    /// </summary>
    [TestFixture]
    [Category("Algorithms")]
    public class SuffixTreePropertyTests
    {
        private static Gen<string> DnaStringGen() =>
            Gen.Select(Gen.Int[0, 3], i => "ACGT"[i]).Array[1, 300].Select(chars => new string(chars));

        private static Gen<string> AnyPatternGen() =>
            Gen.Select(Gen.Int[0, 3], i => "ACGT"[i]).Array[1, 120].Select(chars => new string(chars));

        private static Gen<(string Text, string Pattern)> ExistingPatternGen() =>
            from text in DnaStringGen()
            from start in Gen.Int[0, text.Length - 1]
            from len in Gen.Int[1, text.Length - start]
            select (text, text.Substring(start, len));

        private static List<int> BruteForcePositions(string text, string pattern)
        {
            var result = new List<int>();
            if (pattern.Length == 0)
            {
                result.AddRange(Enumerable.Range(0, text.Length));
                return result;
            }

            for (int i = 0; i <= text.Length - pattern.Length; i++)
            {
                if (text.AsSpan(i, pattern.Length).SequenceEqual(pattern.AsSpan()))
                    result.Add(i);
            }

            return result;
        }

        [Test]
        public void CountOccurrences_Equals_PositionsCount()
        {
            ExistingPatternGen().Sample(sample =>
            {
                var (text, pattern) = sample;
                var tree = SuffixTree.Build(text);

                int count = tree.CountOccurrences(pattern);
                var positions = tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                var expected = BruteForcePositions(text, pattern);

                Assert.That(positions, Is.EqualTo(expected));
                Assert.That(count, Is.EqualTo(expected.Count));
            });
        }

        [Test]
        public void FindAllOccurrences_Matches_BruteForce()
        {
            Gen.Select(DnaStringGen(), AnyPatternGen()).Sample(sample =>
            {
                var (text, pattern) = sample;
                var tree = SuffixTree.Build(text);

                var expected = BruteForcePositions(text, pattern);
                var actual = tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                Assert.That(actual, Is.EqualTo(expected));
            });
        }

        [Test]
        public void SpanOverloads_MatchStringOverloads_OnRandomDna()
        {
            Gen.Select(DnaStringGen(), AnyPatternGen()).Sample(sample =>
            {
                var (text, pattern) = sample;
                var tree = SuffixTree.Build(text);

                Assert.That(tree.Contains(pattern.AsSpan()), Is.EqualTo(tree.Contains(pattern)));
                Assert.That(tree.CountOccurrences(pattern.AsSpan()), Is.EqualTo(tree.CountOccurrences(pattern)));

                var strPos = tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                var spanPos = tree.FindAllOccurrences(pattern.AsSpan()).OrderBy(x => x).ToList();
                Assert.That(spanPos, Is.EqualTo(strPos));
            });
        }
    }
}

