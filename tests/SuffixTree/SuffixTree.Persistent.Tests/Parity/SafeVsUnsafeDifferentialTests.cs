using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests.Parity
{
    /// <summary>
    /// Differential tests: safe path (HeapStorageProvider, _ptr == null)
    /// vs unsafe path (MappedFileStorageProvider, _ptr != null).
    /// </summary>
    [TestFixture]
    [Category("Parity")]
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
                catch { }
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
                for (int i = 0; i < text.Length; i++)
                {
                    for (int len = 1; len <= Math.Min(10, text.Length - i); len++)
                    {
                        string pattern = text.Substring(i, len);
                        Assert.That(mmf.Contains(pattern), Is.EqualTo(heap.Contains(pattern)),
                            $"Contains mismatch for \"{pattern}\"");
                    }
                }

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

                AssertLrsIsValid(heap, text, heapLrs, "heap/fixed");
                AssertLrsIsValid(mmf, text, mmfLrs, "mmf/fixed");
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

                    AssertLcsIsValid(text, other, heapLcs, "heap/fixed");
                    AssertLcsIsValid(text, other, mmfLcs, "mmf/fixed");
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

                AssertAnchorsAreValid(text, query, heapAnchors, 2, "heap/fixed");
                AssertAnchorsAreValid(text, query, mmfAnchors, 2, "mmf/fixed");

                var heapSignatures = AnchorSignatures(text, heapAnchors);
                var mmfSignatures = AnchorSignatures(text, mmfAnchors);
                Assert.That(mmfSignatures, Is.EqualTo(heapSignatures), "Anchor signature mismatch");
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
                    Assert.That(mmf.LeafCount, Is.EqualTo(heap.LeafCount),
                        $"Seed={seed}: LeafCount mismatch");
                    Assert.That(mmf.NodeCount, Is.EqualTo(heap.NodeCount),
                        $"Seed={seed}: NodeCount mismatch");

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

                    string heapLrs = heap.LongestRepeatedSubstring();
                    string mmfLrs = mmf.LongestRepeatedSubstring();
                    AssertLrsIsValid(heap, text, heapLrs, $"heap/random seed={seed}");
                    AssertLrsIsValid(mmf, text, mmfLrs, $"mmf/random seed={seed}");
                    Assert.That(mmfLrs.Length, Is.EqualTo(heapLrs.Length),
                        $"Seed={seed}: LRS length mismatch");

                    string other = GenerateRandom(rng, rng.Next(10, 100), alphabet);
                    string heapLcs = heap.LongestCommonSubstring(other);
                    string mmfLcs = mmf.LongestCommonSubstring(other);
                    AssertLcsIsValid(text, other, heapLcs, $"heap/random seed={seed}");
                    AssertLcsIsValid(text, other, mmfLcs, $"mmf/random seed={seed}");
                    Assert.That(mmfLcs.Length, Is.EqualTo(heapLcs.Length),
                        $"Seed={seed}: LCS length mismatch");

                    string anchorQuery = text.Substring(0, Math.Min(text.Length, 50));
                    var hAnchors = heap.FindExactMatchAnchors(anchorQuery, 3);
                    var mAnchors = mmf.FindExactMatchAnchors(anchorQuery, 3);
                    AssertAnchorsAreValid(text, anchorQuery, hAnchors, 3, $"heap/random seed={seed}");
                    AssertAnchorsAreValid(text, anchorQuery, mAnchors, 3, $"mmf/random seed={seed}");
                    Assert.That(AnchorSignatures(text, mAnchors), Is.EqualTo(AnchorSignatures(text, hAnchors)),
                        $"Seed={seed}: Anchor signature mismatch");
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

        #region Pathological inputs

        [Test]
        public void SingleCharRepeat_HeapVsMmf()
        {
            string text = new string('a', 1000);
            var (heap, mmf) = BuildBoth(text);
            using ((IDisposable)heap)
            using ((IDisposable)mmf)
            {
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
                string[] patterns = { "ACGT", "GATTACA", "AAA", "CCCC", text.Substring(100, 20) };
                foreach (var pattern in patterns)
                {
                    var heapPos = heap.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                    var mmfPos = mmf.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                    Assert.That(mmfPos, Is.EqualTo(heapPos), $"Pattern \"{pattern}\"");
                }

                string heapLrs = heap.LongestRepeatedSubstring();
                string mmfLrs = mmf.LongestRepeatedSubstring();
                Assert.That(mmfLrs.Length, Is.EqualTo(heapLrs.Length));
                Assert.That(mmf.CountOccurrences(mmfLrs), Is.GreaterThanOrEqualTo(2));
            }
        }

        [Test]
        public void AlternatingPattern_HeapVsMmf()
        {
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

        private static void AssertLcsIsValid(string text, string other, string lcs, string context)
        {
            if (lcs.Length > 0)
            {
                Assert.That(text.Contains(lcs, StringComparison.Ordinal), Is.True,
                    $"{context}: LCS \"{lcs}\" not found in text");
                Assert.That(other.Contains(lcs, StringComparison.Ordinal), Is.True,
                    $"{context}: LCS \"{lcs}\" not found in other");
            }

            int expectedLength = LongestCommonSubstringLength(text, other);
            Assert.That(lcs.Length, Is.EqualTo(expectedLength),
                $"{context}: LCS length {lcs.Length} differs from expected {expectedLength}");
        }

        private static void AssertLrsIsValid(ISuffixTree tree, string text, string lrs, string context)
        {
            if (lrs.Length == 0)
            {
                if (text.Length <= 400)
                {
                    Assert.That(ExpectedLrsLength(text), Is.EqualTo(0),
                        $"{context}: expected non-empty LRS");
                }
                return;
            }

            Assert.That(text.Contains(lrs, StringComparison.Ordinal), Is.True,
                $"{context}: LRS \"{lrs}\" is not a substring of text");
            Assert.That(tree.CountOccurrences(lrs), Is.GreaterThanOrEqualTo(2),
                $"{context}: LRS \"{lrs}\" does not repeat");

            if (text.Length <= 400)
            {
                int expectedLength = ExpectedLrsLength(text);
                Assert.That(lrs.Length, Is.EqualTo(expectedLength),
                    $"{context}: LRS length {lrs.Length} differs from expected {expectedLength}");
            }
        }

        private static void AssertAnchorsAreValid(
            string text,
            string query,
            IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> anchors,
            int minLength,
            string context)
        {
            foreach (var anchor in anchors)
            {
                Assert.That(anchor.PositionInText, Is.GreaterThanOrEqualTo(0),
                    $"{context}: anchor PositionInText < 0");
                Assert.That(anchor.PositionInQuery, Is.GreaterThanOrEqualTo(0),
                    $"{context}: anchor PositionInQuery < 0");
                Assert.That(anchor.Length, Is.GreaterThanOrEqualTo(minLength),
                    $"{context}: anchor length < minLength");
                Assert.That(anchor.PositionInText + anchor.Length, Is.LessThanOrEqualTo(text.Length),
                    $"{context}: anchor exceeds text");
                Assert.That(anchor.PositionInQuery + anchor.Length, Is.LessThanOrEqualTo(query.Length),
                    $"{context}: anchor exceeds query");

                string fromText = text.Substring(anchor.PositionInText, anchor.Length);
                string fromQuery = query.Substring(anchor.PositionInQuery, anchor.Length);
                Assert.That(fromText, Is.EqualTo(fromQuery),
                    $"{context}: anchor mismatch at query={anchor.PositionInQuery}, len={anchor.Length}");
            }
        }

        private static List<(int PositionInQuery, int Length, string Substring)> AnchorSignatures(
            string text,
            IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> anchors)
        {
            var signatures = new List<(int PositionInQuery, int Length, string Substring)>(anchors.Count);
            foreach (var anchor in anchors)
            {
                signatures.Add((
                    anchor.PositionInQuery,
                    anchor.Length,
                    text.Substring(anchor.PositionInText, anchor.Length)));
            }

            return signatures;
        }

        private static int LongestCommonSubstringLength(string a, string b)
        {
            if (a.Length == 0 || b.Length == 0)
                return 0;

            var dp = new int[b.Length + 1];
            int best = 0;

            for (int i = 1; i <= a.Length; i++)
            {
                for (int j = b.Length; j >= 1; j--)
                {
                    if (a[i - 1] == b[j - 1])
                    {
                        dp[j] = dp[j - 1] + 1;
                        if (dp[j] > best)
                            best = dp[j];
                    }
                    else
                    {
                        dp[j] = 0;
                    }
                }
            }

            return best;
        }

        private static int ExpectedLrsLength(string text)
        {
            int best = 0;
            for (int i = 0; i < text.Length; i++)
            {
                for (int j = i + 1; j < text.Length; j++)
                {
                    int len = 0;
                    while (i + len < text.Length && j + len < text.Length && text[i + len] == text[j + len])
                    {
                        len++;
                    }

                    if (len > best)
                        best = len;
                }
            }

            return best;
        }

        private static string GenerateRandom(Random rng, int length, string alphabet)
        {
            var chars = new char[length];
            for (int i = 0; i < length; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            return new string(chars);
        }
    }
}
