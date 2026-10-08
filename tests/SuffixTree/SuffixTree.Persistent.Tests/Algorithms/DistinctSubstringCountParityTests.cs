using System;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Persistent.Tests.Algorithms
{
    /// <summary>
    /// Distinct-substring counts from the persistent tree must equal the in-memory tree's.
    /// </summary>
    [TestFixture]
    public class DistinctSubstringCountParityTests
    {
        [Test]
        public void Banana_LiteratureValues()
        {
            using var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("banana"));
            Assert.That(tree.CountDistinctSubstringsByLength(10), Is.EqualTo(new long[] { 1, 3, 3, 3, 3, 2, 1 }));
            Assert.That(tree.CountDistinctSubstrings(), Is.EqualTo(15));
        }

        [Test]
        public void RandomTexts_MatchInMemoryTree([Values(0, 1, 2)] int seed)
        {
            var rng = new Random(seed);
            for (int round = 0; round < 10; round++)
            {
                int n = rng.Next(1, 300);
                string text = new string(Enumerable.Range(0, n).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
                var memory = SuffixTree.Build(text);
                using var persistent = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));

                foreach (int m in new[] { 1, 5, 13, n })
                    Assert.That(persistent.CountDistinctSubstringsByLength(m),
                        Is.EqualTo(memory.CountDistinctSubstringsByLength(m)), $"text={text}, m={m}");
                Assert.That(persistent.CountDistinctSubstrings(), Is.EqualTo(memory.CountDistinctSubstrings()), text);
            }
        }

        [Test]
        public void DisposedTree_Throws()
        {
            var tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("acgt"));
            tree.Dispose();
            Assert.Throws<ObjectDisposedException>(() => tree.CountDistinctSubstrings());
        }
    }
}
