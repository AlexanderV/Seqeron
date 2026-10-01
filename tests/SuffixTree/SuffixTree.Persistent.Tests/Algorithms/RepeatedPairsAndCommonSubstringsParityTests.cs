using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Persistent.Tests.Algorithms
{
    /// <summary>
    /// FindMaximalRepeatedPairs, FindAllDistinctLongestCommonSubstrings and the k-common substring
    /// (generalized suffix tree) must give identical results from the in-memory tree and the persistent
    /// tree (heap, hybrid, memory-mapped, reloaded).
    /// </summary>
    [TestFixture]
    [Category("Parity")]
    public class RepeatedPairsAndCommonSubstringsParityTests
    {
        private static readonly Func<char, bool>?[] Predicates = { null, c => c == 'N', c => "ACGT".IndexOf(c) < 0 };

        private static string Random(Random rng, int length, string alphabet)
            => new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

        [Test]
        public void KnownCases_PersistentTree()
        {
            using var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("xabcyabcwabcyz"));
            Assert.That(tree.FindMaximalRepeatedPairs(1), Is.EqualTo(new[] { (1, 5, 3), (1, 9, 4), (5, 9, 3) }));

            using var lcsTree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("abcxyzabc"));
            var all = lcsTree.FindAllDistinctLongestCommonSubstrings("xyzqabcpxyz");
            Assert.That(all.Select(a => a.Substring), Is.EqualTo(new[] { "xyz", "abc" }));
            Assert.That(all[1].PositionsInText, Is.EqualTo(new[] { 0, 6 }));
            Assert.That(all[0].PositionsInOther, Is.EqualTo(new[] { 0, 8 }));

            var lcsm = SuffixTreeAlgorithms.FindLongestCommonSubstrings(
                new[] { "GATTACA", "TAGACCA", "ATACA" }, 3, s => PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(s)));
            Assert.That(lcsm, Is.EqualTo(new[] { "AC", "CA", "TA" }));
        }

        [Test]
        public void RandomTexts_PersistentEqualsInMemory([Values(0, 1, 2)] int seed)
        {
            var rng = new Random(75_000 + seed);
            for (int round = 0; round < 20; round++)
            {
                string alphabet = (round % 3) switch { 0 => "AC", 1 => "ACGT", _ => "ACGTN" };
                string text = Random(rng, rng.Next(0, 400), alphabet);
                string other = Random(rng, rng.Next(0, 200), alphabet);
                int minLength = rng.Next(1, 8);
                var unique = Predicates[rng.Next(Predicates.Length)];
                var memory = SuffixTree.Build(text);
                var pairs = memory.FindMaximalRepeatedPairs(minLength, unique);
                var lcs = memory.FindAllDistinctLongestCommonSubstrings(other);
                using (var heap = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text)))
                {
                    Assert.That(heap.FindMaximalRepeatedPairs(minLength, unique), Is.EqualTo(pairs), $"heap {seed}/{round}");
                    AssertSameLcs(lcs, heap.FindAllDistinctLongestCommonSubstrings(other), $"heap {seed}/{round}");
                }
                using (var hybrid = PersistentSuffixTreeFactory.CreateCore(new StringTextSource(text), null, 88 + 24 * 3))
                {
                    Assert.That(hybrid.FindMaximalRepeatedPairs(minLength, unique), Is.EqualTo(pairs), $"hybrid {seed}/{round}");
                    AssertSameLcs(lcs, hybrid.FindAllDistinctLongestCommonSubstrings(other), $"hybrid {seed}/{round}");
                }
            }
        }

        [Test]
        public void KCommonSubstring_PersistentEqualsInMemory([Values(0, 1)] int seed)
        {
            var rng = new Random(76_000 + seed);
            for (int round = 0; round < 15; round++)
            {
                int k = rng.Next(1, 8);
                var texts = Enumerable.Range(0, k).Select(_ => Random(rng, rng.Next(0, 120), "ACG")).ToArray();
                Func<string, ISuffixTree> persistent = s => PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(s));
                Assert.That(SuffixTreeAlgorithms.LongestCommonSubstringLengthsBySupport(texts, persistent),
                    Is.EqualTo(SuffixTree.LongestCommonSubstringLengthsBySupport(texts)), $"{seed}/{round}");
                for (int q = 1; q <= k; q++)
                    Assert.That(SuffixTreeAlgorithms.FindLongestCommonSubstrings(texts, q, persistent),
                        Is.EqualTo(SuffixTree.FindLongestCommonSubstrings(texts, q)), $"{seed}/{round} q={q}");
            }
        }

        [Test]
        public void MemoryMappedTree_EqualsInMemory()
        {
            var rng = new Random(77_000);
            string text = Random(rng, 3000, "ACGTN");
            string other = Random(rng, 1500, "ACGT");
            var memory = SuffixTree.Build(text);
            var pairs = memory.FindMaximalRepeatedPairs(6, c => c == 'N');
            var lcs = memory.FindAllDistinctLongestCommonSubstrings(other);
            var path = Path.Combine(Path.GetTempPath(), $"mrp_{Guid.NewGuid():N}.st");
            try
            {
                using (var mmf = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path))
                {
                    Assert.That(mmf.FindMaximalRepeatedPairs(6, c => c == 'N'), Is.EqualTo(pairs), "mmf");
                    AssertSameLcs(lcs, mmf.FindAllDistinctLongestCommonSubstrings(other), "mmf");
                }
                using var loaded = PersistentSuffixTreeFactory.LoadPersistent(path);
                Assert.That(loaded.FindMaximalRepeatedPairs(6, c => c == 'N'), Is.EqualTo(pairs), "reload");
                AssertSameLcs(lcs, loaded.FindAllDistinctLongestCommonSubstrings(other), "reload");
            }
            finally { File.Delete(path); }
        }

        [Test]
        public void DisposedTree_Throws()
        {
            var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("banana"));
            tree.Dispose();
            Assert.Throws<ObjectDisposedException>(() => tree.FindMaximalRepeatedPairs(1));
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllDistinctLongestCommonSubstrings("ana"));
        }

        private static void AssertSameLcs(
            IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> expected,
            IReadOnlyList<(string Substring, IReadOnlyList<int> PositionsInText, IReadOnlyList<int> PositionsInOther)> actual,
            string label)
        {
            Assert.That(actual.Select(a => a.Substring), Is.EqualTo(expected.Select(e => e.Substring)), label);
            for (int i = 0; i < expected.Count; i++)
            {
                Assert.That(actual[i].PositionsInText, Is.EqualTo(expected[i].PositionsInText), $"{label} #{i} text");
                Assert.That(actual[i].PositionsInOther, Is.EqualTo(expected[i].PositionsInOther), $"{label} #{i} other");
            }
        }
    }
}
