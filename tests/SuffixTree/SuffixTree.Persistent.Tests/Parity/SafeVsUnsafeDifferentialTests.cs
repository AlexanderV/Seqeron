using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests
{
    /// <summary>
    /// Differential tests: safe path (HeapStorageProvider, _ptr == null)
    /// vs unsafe path (MappedFileStorageProvider, _ptr != null).
    /// 
    /// Both paths must produce identical results for all query operations.
    /// This validates that CollectLeavesUnsafe, TryGetChildFast, ReadUInt32Fast,
    /// ReadInt64Fast produce the same results as their safe counterparts.
    /// </summary>
    [TestFixture]
    public class SafeVsUnsafeDifferentialTests
    {
        private string? _tempDir;

        [SetUp]
        public void SetUp()
        {
            _tempDir = Path.Combine(Path.GetTempPath(), "SuffixTree_Diff_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_tempDir);
        }

        [TearDown]
        public void TearDown()
        {
            if (_tempDir != null && Directory.Exists(_tempDir))
            {
                try { Directory.Delete(_tempDir, true); }
                catch { /* best effort cleanup */ }
            }
        }

        private (ISuffixTree Heap, ISuffixTree Mmf) BuildBoth(string text)
        {
            var textSource = new StringTextSource(text);
            var heap = PersistentSuffixTreeFactory.Create(textSource);

            var textSource2 = new StringTextSource(text);
            var filePath = Path.Combine(_tempDir!, $"tree_{Guid.NewGuid():N}.st");
            var mmf = PersistentSuffixTreeFactory.Create(textSource2, filePath);

            return (heap, mmf);
        }

        #region Fixed inputs

        private static readonly string[] FixedInputs = new[]
        {
            "banana",
            "mississippi",
            "abracadabra",
            "aaaaaaaaaa",
            "abcdefghij",
            "a",
            "xyxyxyxyxy",
            "ACGTACGTACGTACGT",
            "aababcabcdabcde",
            "the quick brown fox jumps over the lazy dog"
        };

        [Test]
        public void Contains_HeapVsMmf_Identical([ValueSource(nameof(FixedInputs))] string text)
        {
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                // Test every substring
                for (int i = 0; i < text.Length; i++)
                {
                    for (int len = 1; len <= Math.Min(10, text.Length - i); len++)
                    {
                        string pattern = text.Substring(i, len);
                        Assert.That(mmf.Contains(pattern), Is.EqualTo(heap.Contains(pattern)),
                            $"Contains mismatch for \"{pattern}\"");
                    }
                }

                // Test non-existing patterns
                Assert.That(mmf.Contains("ZZZZZ"), Is.EqualTo(heap.Contains("ZZZZZ")));
                Assert.That(mmf.Contains(""), Is.EqualTo(heap.Contains("")));
            }
        }

        [Test]
        public void CountOccurrences_HeapVsMmf_Identical([ValueSource(nameof(FixedInputs))] string text)
        {
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                for (int i = 0; i < text.Length; i++)
                {
                    for (int len = 1; len <= Math.Min(8, text.Length - i); len++)
                    {
                        string pattern = text.Substring(i, len);
                        Assert.That(mmf.CountOccurrences(pattern), Is.EqualTo(heap.CountOccurrences(pattern)),
                            $"CountOccurrences mismatch for \"{pattern}\"");
                    }
                }
            }
        }

        [Test]
        public void FindAllOccurrences_HeapVsMmf_Identical([ValueSource(nameof(FixedInputs))] string text)
        {
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                for (int i = 0; i < text.Length; i++)
                {
                    for (int len = 1; len <= Math.Min(8, text.Length - i); len++)
                    {
                        string pattern = text.Substring(i, len);
                        var heapResult = heap.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                        var mmfResult = mmf.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                        Assert.That(mmfResult, Is.EqualTo(heapResult),
                            $"FindAllOccurrences mismatch for \"{pattern}\"");
                    }
                }
            }
        }

        [Test]
        public void LRS_HeapVsMmf_Identical([ValueSource(nameof(FixedInputs))] string text)
        {
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                var heapLrs = heap.LongestRepeatedSubstring();
                var mmfLrs = mmf.LongestRepeatedSubstring();
                Assert.That(mmfLrs.Length, Is.EqualTo(heapLrs.Length),
                    $"LRS length mismatch: heap=\"{heapLrs}\", mmf=\"{mmfLrs}\"");
            }
        }

        [Test]
        public void LCS_HeapVsMmf_Identical([ValueSource(nameof(FixedInputs))] string text)
        {
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                string[] others = { "ana_missi_abra", "ACGTACGT", "quick lazy", "xyzxyz", "" };
                foreach (var other in others)
                {
                    var heapLcs = heap.LongestCommonSubstring(other);
                    var mmfLcs = mmf.LongestCommonSubstring(other);
                    Assert.That(mmfLcs.Length, Is.EqualTo(heapLcs.Length),
                        $"LCS length mismatch for other=\"{other}\": heap=\"{heapLcs}\", mmf=\"{mmfLcs}\"");
                }
            }
        }

        [Test]
        public void FindExactMatchAnchors_HeapVsMmf_Identical([ValueSource(nameof(FixedInputs))] string text)
        {
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                string query = "banana_mississippi_abracadabra";
                var heapAnchors = heap.FindExactMatchAnchors(query, 2);
                var mmfAnchors = mmf.FindExactMatchAnchors(query, 2);

                Assert.That(mmfAnchors.Count, Is.EqualTo(heapAnchors.Count),
                    "Anchor count mismatch");
                for (int i = 0; i < heapAnchors.Count; i++)
                {
                    Assert.That(mmfAnchors[i].PositionInQuery, Is.EqualTo(heapAnchors[i].PositionInQuery),
                        $"Anchor[{i}] PositionInQuery mismatch");
                    Assert.That(mmfAnchors[i].Length, Is.EqualTo(heapAnchors[i].Length),
                        $"Anchor[{i}] Length mismatch");
                }
            }
        }

        [Test]
        public void Statistics_HeapVsMmf_Identical([ValueSource(nameof(FixedInputs))] string text)
        {
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(mmf.LeafCount, Is.EqualTo(heap.LeafCount), "LeafCount");
                    Assert.That(mmf.NodeCount, Is.EqualTo(heap.NodeCount), "NodeCount");
                    Assert.That(mmf.MaxDepth, Is.EqualTo(heap.MaxDepth), "MaxDepth");
                    Assert.That(mmf.IsEmpty, Is.EqualTo(heap.IsEmpty), "IsEmpty");
                });
            }
        }

        #endregion

        #region Random inputs — comprehensive differential

        [Test]
        public void AllQueries_HeapVsMmf_RandomInputs()
        {
            for (int seed = 0; seed < 30; seed++)
            {
                var rng = new Random(seed);
                int len = rng.Next(10, 300);
                string alphabet = seed % 3 == 0 ? "ab" : seed % 3 == 1 ? "ACGT" : "abcdefgh";
                string text = GenerateRandom(rng, len, alphabet);

                var (heap, mmf) = BuildBoth(text);
                using ((IDisposable)heap)
                using ((IDisposable)mmf)
                {
                    // Statistics
                    Assert.That(mmf.LeafCount, Is.EqualTo(heap.LeafCount),
                        $"Seed={seed}: LeafCount mismatch");
                    Assert.That(mmf.NodeCount, Is.EqualTo(heap.NodeCount),
                        $"Seed={seed}: NodeCount mismatch");

                    // Random queries
                    for (int q = 0; q < 20; q++)
                    {
                        int start = rng.Next(text.Length);
                        int plen = rng.Next(1, Math.Min(15, text.Length - start + 1));
                        string pattern = text.Substring(start, plen);

                        Assert.That(mmf.Contains(pattern), Is.EqualTo(heap.Contains(pattern)),
                            $"Seed={seed}: Contains mismatch for \"{pattern}\"");

                        Assert.That(mmf.CountOccurrences(pattern), Is.EqualTo(heap.CountOccurrences(pattern)),
                            $"Seed={seed}: Count mismatch for \"{pattern}\"");

                        var heapPos = heap.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                        var mmfPos = mmf.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                        Assert.That(mmfPos, Is.EqualTo(heapPos),
                            $"Seed={seed}: FindAll mismatch for \"{pattern}\"");
                    }

                    // LRS
                    Assert.That(mmf.LongestRepeatedSubstring().Length,
                        Is.EqualTo(heap.LongestRepeatedSubstring().Length),
                        $"Seed={seed}: LRS length mismatch");

                    // LCS with random other
                    string other = GenerateRandom(rng, rng.Next(10, 100), alphabet);
                    Assert.That(mmf.LongestCommonSubstring(other).Length,
                        Is.EqualTo(heap.LongestCommonSubstring(other).Length),
                        $"Seed={seed}: LCS length mismatch");

                    // Anchors
                    string anchorQuery = text.Substring(0, Math.Min(text.Length, 50));
                    var hAnchors = heap.FindExactMatchAnchors(anchorQuery, 3);
                    var mAnchors = mmf.FindExactMatchAnchors(anchorQuery, 3);
                    Assert.That(mAnchors.Count, Is.EqualTo(hAnchors.Count),
                        $"Seed={seed}: Anchor count mismatch");
                    for (int i = 0; i < hAnchors.Count; i++)
                    {
                        Assert.That(mAnchors[i].PositionInQuery, Is.EqualTo(hAnchors[i].PositionInQuery),
                            $"Seed={seed}: Anchor[{i}] PosInQuery mismatch");
                        Assert.That(mAnchors[i].Length, Is.EqualTo(hAnchors[i].Length),
                            $"Seed={seed}: Anchor[{i}] Length mismatch");
                    }
                }
            }
        }

        #endregion

        #region Topology differential

        [Test]
        public void Topology_HeapVsMmf_Identical([ValueSource(nameof(FixedInputs))] string text)
        {
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                var heapNodes = new List<string>();
                var mmfNodes = new List<string>();

                heap.Traverse(new FlatVisitor(heapNodes));
                mmf.Traverse(new FlatVisitor(mmfNodes));

                Assert.That(mmfNodes.Count, Is.EqualTo(heapNodes.Count),
                    "Traversal node count mismatch");
                for (int i = 0; i < heapNodes.Count; i++)
                {
                    Assert.That(mmfNodes[i], Is.EqualTo(heapNodes[i]),
                        $"Traversal mismatch at index {i}");
                }
            }
        }

        private class FlatVisitor : ISuffixTreeVisitor
        {
            private readonly List<string> _records;
            public FlatVisitor(List<string> records) => _records = records;

            public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
                => _records.Add($"N({startIndex},{endIndex},{leafCount},{childCount},{depth})");

            public void EnterBranch(int key)
                => _records.Add($"E({key})");

            public void ExitBranch()
                => _records.Add("X");
        }

        #endregion

        #region Pathological inputs (stress the unsafe code)

        [Test]
        public void SingleCharRepeat_HeapVsMmf()
        {
            string text = new string('a', 1000);
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                // This creates maximum branching at root with many children
                string pattern = "aaa";
                var heapPos = heap.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                var mmfPos = mmf.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                Assert.That(mmfPos, Is.EqualTo(heapPos));
                Assert.That(mmf.CountOccurrences(pattern), Is.EqualTo(heap.CountOccurrences(pattern)));
            }
        }

        [Test]
        public void LongDna_HeapVsMmf()
        {
            var rng = new Random(42);
            string text = GenerateRandom(rng, 5000, "ACGT");
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                // Verify multiple patterns
                string[] patterns = { "ACGT", "GATTACA", "AAA", "CCCC", text.Substring(100, 20) };
                foreach (var p in patterns)
                {
                    var heapPos = heap.FindAllOccurrences(p).OrderBy(x => x).ToList();
                    var mmfPos = mmf.FindAllOccurrences(p).OrderBy(x => x).ToList();
                    Assert.That(mmfPos, Is.EqualTo(heapPos), $"Pattern \"{p}\"");
                }

                Assert.That(mmf.LongestRepeatedSubstring().Length,
                    Is.EqualTo(heap.LongestRepeatedSubstring().Length));
            }
        }

        [Test]
        public void AlternatingPattern_HeapVsMmf()
        {
            // Alternating pattern creates many internal nodes with 2 children each
            string text = string.Concat(Enumerable.Repeat("ab", 500));
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
                Assert.That(mmf.LeafCount, Is.EqualTo(heap.LeafCount));
                Assert.That(mmf.NodeCount, Is.EqualTo(heap.NodeCount));

                var heapAll = heap.FindAllOccurrences("abab").OrderBy(x => x).ToList();
                var mmfAll = mmf.FindAllOccurrences("abab").OrderBy(x => x).ToList();
                Assert.That(mmfAll, Is.EqualTo(heapAll));
            }
        }

        #endregion

        private static string GenerateRandom(Random rng, int length, string alphabet)
        {
            var chars = new char[length];
            for (int i = 0; i < length; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            return new string(chars);
        }
    }
}
