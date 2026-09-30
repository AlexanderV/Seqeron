using System;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Persistent.Tests.Algorithms
{
    /// <summary>
    /// FindAllLongestRepeatedSubstrings must return identical lists (substrings, ascending positions,
    /// first-occurrence order) from the in-memory tree and the persistent tree (heap, hybrid, MMF, reload).
    /// </summary>
    [TestFixture]
    [Category("Parity")]
    public class AllLongestRepeatedSubstringsParityTests
    {
        [Test]
        public void KnownTies_PersistentTree()
        {
            using var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("abcxbcaxcab"));
            var all = tree.FindAllLongestRepeatedSubstrings();
            Assert.That(all.Select(a => a.Substring), Is.EqualTo(new[] { "ab", "bc", "ca" }));
            Assert.That(all.Select(a => a.Positions.ToArray()), Is.EqualTo(new[] { new[] { 0, 9 }, new[] { 1, 4 }, new[] { 5, 8 } }));
            Assert.That(tree.LongestRepeatedSubstring().Length, Is.EqualTo(2));
        }

        [Test]
        public void RandomTexts_PersistentEqualsInMemory([Values(0, 1, 2)] int seed)
        {
            var rng = new Random(73_000 + seed);
            for (int round = 0; round < 20; round++)
            {
                string alphabet = round % 2 == 0 ? "AC" : "ACGT";
                string text = new(Enumerable.Range(0, rng.Next(0, 500)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
                var expected = SuffixTree.Build(text).FindAllLongestRepeatedSubstrings();
                using (var heap = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text)))
                    AssertSame(expected, heap.FindAllLongestRepeatedSubstrings(), $"heap {seed}/{round}");
                using (var hybrid = PersistentSuffixTreeFactory.CreateCore(new StringTextSource(text), null, 88 + 24 * 3))
                    AssertSame(expected, hybrid.FindAllLongestRepeatedSubstrings(), $"hybrid {seed}/{round}");
            }
        }

        [Test]
        public void MemoryMappedTree_EqualsInMemory()
        {
            var rng = new Random(74_000);
            string text = new(Enumerable.Range(0, 3000).Select(_ => "ACG"[rng.Next(3)]).ToArray());
            var expected = SuffixTree.Build(text).FindAllLongestRepeatedSubstrings();
            var path = Path.Combine(Path.GetTempPath(), $"lrs_all_{Guid.NewGuid():N}.st");
            try
            {
                using (var mmf = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path))
                    AssertSame(expected, mmf.FindAllLongestRepeatedSubstrings(), "mmf");
                using var loaded = PersistentSuffixTreeFactory.LoadPersistent(path);
                AssertSame(expected, loaded.FindAllLongestRepeatedSubstrings(), "mmf reload");
            }
            finally { File.Delete(path); }
        }

        [Test]
        public void DisposedTree_Throws()
        {
            var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("banana"));
            tree.Dispose();
            Assert.Throws<ObjectDisposedException>(() => tree.FindAllLongestRepeatedSubstrings());
        }

        private static void AssertSame(
            System.Collections.Generic.IReadOnlyList<(string Substring, System.Collections.Generic.IReadOnlyList<int> Positions)> expected,
            System.Collections.Generic.IReadOnlyList<(string Substring, System.Collections.Generic.IReadOnlyList<int> Positions)> actual,
            string label)
        {
            Assert.That(actual.Select(a => a.Substring), Is.EqualTo(expected.Select(e => e.Substring)), label);
            for (int i = 0; i < expected.Count; i++)
                Assert.That(actual[i].Positions, Is.EqualTo(expected[i].Positions), $"{label} #{i}");
        }
    }
}
