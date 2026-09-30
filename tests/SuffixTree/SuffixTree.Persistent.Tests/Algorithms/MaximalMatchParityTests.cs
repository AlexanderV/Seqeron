using System;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Persistent.Tests.Algorithms
{
    /// <summary>
    /// FindMaximalExactMatches / FindMaximalUniqueMatches must return identical lists from the
    /// in-memory tree and the persistent tree (heap, hybrid compact/large layout, memory-mapped file).
    /// </summary>
    [TestFixture]
    [Category("Parity")]
    public class MaximalMatchParityTests
    {
        // mummer -maxmatch / -mumreference / -mum -l 4 (MUMmer 3.23), 1-based output converted to 0-based.
        [Test]
        public void MummerExample_PersistentTree()
        {
            const string reference = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG";
            const string query = "TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG";
            using var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(reference));
            Assert.That(tree.FindMaximalExactMatches(query, 4), Is.EqualTo(new[]
                { (0, 0, 13), (21, 6, 4), (15, 10, 4), (3, 17, 9), (21, 20, 4), (25, 26, 12), (22, 32, 4) }));
            Assert.That(tree.FindMaximalUniqueMatches(query, 4, MumUniqueness.Reference), Is.EqualTo(new[]
                { (0, 0, 13), (15, 10, 4), (3, 17, 9), (25, 26, 12) }));
            Assert.That(tree.FindMaximalUniqueMatches(query, 4, MumUniqueness.Both), Is.EqualTo(new[]
                { (0, 0, 13), (15, 10, 4), (25, 26, 12) }));
        }

        [Test]
        public void RandomPairs_PersistentEqualsInMemory([Values(0, 1, 2)] int seed)
        {
            var rng = new Random(71_000 + seed);
            for (int round = 0; round < 8; round++)
            {
                string alphabet = round % 2 == 0 ? "AC" : "ACGT";
                string text = RandomString(rng, rng.Next(1, 400), alphabet);
                string query = round % 3 == 0
                    ? RandomString(rng, rng.Next(1, 400), alphabet)
                    : new string(text.Skip(rng.Next(text.Length)).Select(c => rng.Next(12) == 0 ? alphabet[rng.Next(alphabet.Length)] : c).ToArray());
                int minLength = rng.Next(1, 9);
                var memory = SuffixTree.Build(text);

                using (var heap = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text)))
                    AssertSame(memory, heap, query, minLength, $"heap seed {seed} round {round}");
                using (var hybrid = PersistentSuffixTreeFactory.CreateCore(new StringTextSource(text), null, 88 + 24 * 3))
                    AssertSame(memory, hybrid, query, minLength, $"hybrid seed {seed} round {round}");
            }
        }

        [Test]
        public void MemoryMappedTree_EqualsInMemory()
        {
            var rng = new Random(72_000);
            string text = RandomString(rng, 1500, "ACG");
            string query = string.Concat(text.AsSpan(200, 700), RandomString(rng, 300, "ACG"));
            var path = Path.Combine(Path.GetTempPath(), $"mem_parity_{Guid.NewGuid():N}.st");
            try
            {
                using (var mmf = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path))
                    AssertSame(SuffixTree.Build(text), mmf, query, 6, "mmf");
                using var loaded = PersistentSuffixTreeFactory.LoadPersistent(path);
                AssertSame(SuffixTree.Build(text), loaded, query, 6, "mmf reload");
            }
            finally { File.Delete(path); }
        }

        [Test]
        public void DisposedTree_Throws()
        {
            var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("banana"));
            tree.Dispose();
            Assert.Throws<ObjectDisposedException>(() => tree.FindMaximalExactMatches("ana", 2));
            Assert.Throws<ObjectDisposedException>(() => tree.FindMaximalUniqueMatches("ana", 2));
        }

        private static void AssertSame(ISuffixTree expected, ISuffixTree actual, string query, int minLength, string label)
        {
            Assert.That(actual.FindMaximalExactMatches(query, minLength),
                Is.EqualTo(expected.FindMaximalExactMatches(query, minLength)), label + " MEM");
            foreach (var mode in new[] { MumUniqueness.Both, MumUniqueness.Reference })
                Assert.That(actual.FindMaximalUniqueMatches(query, minLength, mode),
                    Is.EqualTo(expected.FindMaximalUniqueMatches(query, minLength, mode)), $"{label} MUM {mode}");
        }

        private static string RandomString(Random rng, int length, string alphabet)
            => new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
    }
}
