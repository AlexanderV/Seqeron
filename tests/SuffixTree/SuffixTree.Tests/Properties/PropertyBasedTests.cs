using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace SuffixTree.Tests.Properties
{
    /// <summary>
    /// Property-based tests for suffix tree.
    /// Each test verifies a universal mathematical property that must hold
    /// for ANY input string, using randomized inputs with fixed seeds for reproducibility.
    /// </summary>
    [TestFixture]
    [Category("Properties")]
    public class PropertyBasedTests
    {
        private const int Iterations = 50;

        #region Generators

        private static string GenerateRandom(Random rng, int minLen, int maxLen, string alphabet)
        {
            int len = rng.Next(minLen, maxLen + 1);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
            return new string(chars);
        }

        private static IEnumerable<(string Text, int Seed)> RandomInputs(
            string alphabet, int minLen = 1, int maxLen = 500)
        {
            for (int seed = 0; seed < Iterations; seed++)
            {
                var rng = new Random(seed);
                yield return (GenerateRandom(rng, minLen, maxLen, alphabet), seed);
            }
        }

        #endregion

        #region P1: LeafCount = Text.Length

        /// <summary>
        /// Gusfield §5.2: A suffix tree for text T of length n with unique
        /// terminator has exactly n leaves.
        /// </summary>
        [Test]
        public void P1_LeafCount_EqualsTextLength_SmallAlphabet()
        {
            foreach (var (text, seed) in RandomInputs("ab", 1, 200))
            {
                var tree = SuffixTree.Build(text);
                Assert.That(tree.LeafCount, Is.EqualTo(text.Length),
                    $"Seed={seed}, Text[0..20]=\"{text[..Math.Min(20, text.Length)]}\"");
            }
        }

        [Test]
        public void P1_LeafCount_EqualsTextLength_DnaAlphabet()
        {
            foreach (var (text, seed) in RandomInputs("ACGT", 1, 500))
            {
                var tree = SuffixTree.Build(text);
                Assert.That(tree.LeafCount, Is.EqualTo(text.Length),
                    $"Seed={seed}");
            }
        }

        #endregion

        #region P2: NodeCount ∈ [n+1, 2n+1]

        /// <summary>
        /// Gusfield bound: n+1 ≤ NodeCount ≤ 2n+1 where n = |T|.
        /// Lower bound: all leaves + root. Upper: at most n-1 internal nodes.
        /// </summary>
        [Test]
        public void P2_NodeCount_WithinGusfieldBounds()
        {
            foreach (var (text, seed) in RandomInputs("abcd"))
            {
                var tree = SuffixTree.Build(text);
                int n = text.Length;
                Assert.That(tree.NodeCount, Is.GreaterThanOrEqualTo(n + 1),
                    $"Seed={seed}: NodeCount below lower bound");
                Assert.That(tree.NodeCount, Is.LessThanOrEqualTo(2 * n + 1),
                    $"Seed={seed}: NodeCount above upper bound");
            }
        }

        #endregion

        #region P3: ∀ suffix s of T: Contains(s) = true

        /// <summary>
        /// Every suffix of T must be findable in the tree.
        /// </summary>
        [Test]
        public void P3_AllSuffixes_AreContained()
        {
            foreach (var (text, seed) in RandomInputs("abc", 1, 100))
            {
                var tree = SuffixTree.Build(text);
                for (int i = 0; i < text.Length; i++)
                {
                    Assert.That(tree.Contains(text.Substring(i)), Is.True,
                        $"Seed={seed}: suffix at {i} not found");
                }
            }
        }

        #endregion

        #region P4: CountOccurrences(p) = |FindAllOccurrences(p)|

        /// <summary>
        /// Count (via LeafCount read) must be consistent with DFS leaf collection.
        /// </summary>
        [Test]
        public void P4_Count_EqualsPositionsCount()
        {
            foreach (var (text, seed) in RandomInputs("abc", 5, 200))
            {
                var tree = SuffixTree.Build(text);
                var rng = new Random(seed + 10000);

                // Test 10 random substrings per input
                for (int q = 0; q < 10; q++)
                {
                    int start = rng.Next(text.Length);
                    int len = rng.Next(1, Math.Min(20, text.Length - start + 1));
                    string pattern = text.Substring(start, len);

                    int count = tree.CountOccurrences(pattern);
                    var positions = tree.FindAllOccurrences(pattern);

                    Assert.That(count, Is.EqualTo(positions.Count),
                        $"Seed={seed}, pattern=\"{pattern}\"");
                }
            }
        }

        #endregion

        #region P5: Contains(p) ⟺ CountOccurrences(p) > 0

        [Test]
        public void P5_Contains_EquivalentToCountPositive()
        {
            foreach (var (text, seed) in RandomInputs("abcd", 5, 200))
            {
                var tree = SuffixTree.Build(text);
                var rng = new Random(seed + 20000);

                for (int q = 0; q < 10; q++)
                {
                    // Mix of existing and non-existing patterns
                    string pattern;
                    if (rng.Next(2) == 0)
                    {
                        int start = rng.Next(text.Length);
                        int len = rng.Next(1, Math.Min(10, text.Length - start + 1));
                        pattern = text.Substring(start, len);
                    }
                    else
                    {
                        pattern = GenerateRandom(rng, 1, 10, "abcde"); // 'e' not always in text
                    }

                    bool contains = tree.Contains(pattern);
                    int count = tree.CountOccurrences(pattern);

                    Assert.That(contains, Is.EqualTo(count > 0),
                        $"Seed={seed}, pattern=\"{pattern}\": Contains={contains}, Count={count}");
                }
            }
        }

        #endregion

        #region P6: FindAllOccurrences positions are valid

        /// <summary>
        /// Every reported position p must satisfy: T[p..p+|pattern|] = pattern.
        /// </summary>
        [Test]
        public void P6_AllPositions_AreCorrect()
        {
            foreach (var (text, seed) in RandomInputs("abcd", 5, 200))
            {
                var tree = SuffixTree.Build(text);
                var rng = new Random(seed + 30000);

                for (int q = 0; q < 10; q++)
                {
                    int start = rng.Next(text.Length);
                    int len = rng.Next(1, Math.Min(15, text.Length - start + 1));
                    string pattern = text.Substring(start, len);

                    var positions = tree.FindAllOccurrences(pattern);

                    foreach (int pos in positions)
                    {
                        Assert.That(pos, Is.GreaterThanOrEqualTo(0),
                            $"Seed={seed}: negative position");
                        Assert.That(pos + pattern.Length, Is.LessThanOrEqualTo(text.Length),
                            $"Seed={seed}: position overflows text");
                        Assert.That(text.Substring(pos, pattern.Length), Is.EqualTo(pattern),
                            $"Seed={seed}: mismatch at pos {pos}");
                    }
                }
            }
        }

        #endregion

        #region P7: FindAllOccurrences is complete (no missed occurrences)

        /// <summary>
        /// Brute-force scan must not find positions that FindAllOccurrences missed.
        /// </summary>
        [Test]
        public void P7_FindAll_IsComplete_VsBruteForce()
        {
            foreach (var (text, seed) in RandomInputs("ab", 1, 100))
            {
                var tree = SuffixTree.Build(text);
                var rng = new Random(seed + 40000);

                for (int q = 0; q < 5; q++)
                {
                    int start = rng.Next(text.Length);
                    int len = rng.Next(1, Math.Min(8, text.Length - start + 1));
                    string pattern = text.Substring(start, len);

                    var treePositions = tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();

                    // Brute force
                    var brutePositions = new List<int>();
                    for (int i = 0; i <= text.Length - pattern.Length; i++)
                    {
                        if (text.Substring(i, pattern.Length) == pattern)
                            brutePositions.Add(i);
                    }

                    Assert.That(treePositions, Is.EqualTo(brutePositions),
                        $"Seed={seed}, pattern=\"{pattern}\"");
                }
            }
        }

        #endregion

        #region P8: LCS(A,B).Length = LCS(B,A).Length (symmetry)

        /// <summary>
        /// LCS length is symmetric: |LCS(A,B)| = |LCS(B,A)|.
        /// The actual substring may differ but length must match.
        /// </summary>
        [Test]
        public void P8_LCS_LengthIsSymmetric()
        {
            foreach (var (text, seed) in RandomInputs("abc", 5, 150))
            {
                var rng = new Random(seed + 50000);
                string other = GenerateRandom(rng, 5, 150, "abc");

                var treeA = SuffixTree.Build(text);
                var treeB = SuffixTree.Build(other);

                string lcsAB = treeA.LongestCommonSubstring(other);
                string lcsBA = treeB.LongestCommonSubstring(text);

                Assert.That(lcsAB.Length, Is.EqualTo(lcsBA.Length),
                    $"Seed={seed}: |LCS(A,B)|={lcsAB.Length} ≠ |LCS(B,A)|={lcsBA.Length}");
            }
        }

        #endregion

        #region P9: LCS(A,B) exists in both A and B

        [Test]
        public void P9_LCS_ExistsInBothStrings()
        {
            foreach (var (text, seed) in RandomInputs("abcd", 5, 200))
            {
                var rng = new Random(seed + 60000);
                string other = GenerateRandom(rng, 5, 200, "abcd");

                var tree = SuffixTree.Build(text);
                string lcs = tree.LongestCommonSubstring(other);

                if (lcs.Length > 0)
                {
                    Assert.That(text.Contains(lcs), Is.True,
                        $"Seed={seed}: LCS \"{lcs}\" not in text");
                    Assert.That(other.Contains(lcs), Is.True,
                        $"Seed={seed}: LCS \"{lcs}\" not in other");
                }
            }
        }

        #endregion

        #region P10: LRS occurs at least twice

        /// <summary>
        /// The longest repeated substring must occur ≥ 2 times (by definition).
        /// </summary>
        [Test]
        public void P10_LRS_OccursAtLeastTwice()
        {
            foreach (var (text, seed) in RandomInputs("abc", 2, 200))
            {
                var tree = SuffixTree.Build(text);
                string lrs = tree.LongestRepeatedSubstring();

                if (lrs.Length > 0)
                {
                    int count = tree.CountOccurrences(lrs);
                    Assert.That(count, Is.GreaterThanOrEqualTo(2),
                        $"Seed={seed}: LRS \"{lrs}\" occurs only {count} time(s)");
                }
            }
        }

        #endregion

        #region P11: No longer repeated substring exists

        /// <summary>
        /// There must be no repeated substring longer than LRS.
        /// Verified by checking that no substring of length LRS+1 occurs ≥ 2 times.
        /// </summary>
        [Test]
        public void P11_LRS_IsMaximal()
        {
            foreach (var (text, seed) in RandomInputs("ab", 2, 80))
            {
                var tree = SuffixTree.Build(text);
                string lrs = tree.LongestRepeatedSubstring();
                int targetLen = lrs.Length + 1;

                if (targetLen > text.Length)
                    continue;

                // Check that no substring of length LRS+1 repeats
                bool foundLonger = false;
                for (int i = 0; i <= text.Length - targetLen; i++)
                {
                    string sub = text.Substring(i, targetLen);
                    if (tree.CountOccurrences(sub) >= 2)
                    {
                        foundLonger = true;
                        break;
                    }
                }

                Assert.That(foundLonger, Is.False,
                    $"Seed={seed}: found repeated substring longer than LRS \"{lrs}\"");
            }
        }

        #endregion

        #region P12: Span overload ≡ string overload

        /// <summary>
        /// ReadOnlySpan overloads must produce identical results to string overloads.
        /// </summary>
        [Test]
        public void P12_SpanOverloads_MatchStringOverloads()
        {
            foreach (var (text, seed) in RandomInputs("abcd", 5, 200))
            {
                var tree = SuffixTree.Build(text);
                var rng = new Random(seed + 70000);

                for (int q = 0; q < 5; q++)
                {
                    int start = rng.Next(text.Length);
                    int len = rng.Next(1, Math.Min(15, text.Length - start + 1));
                    string pattern = text.Substring(start, len);

                    bool containsStr = tree.Contains(pattern);
                    bool containsSpan = tree.Contains(pattern.AsSpan());
                    Assert.That(containsSpan, Is.EqualTo(containsStr),
                        $"Seed={seed}: Contains mismatch for \"{pattern}\"");

                    int countStr = tree.CountOccurrences(pattern);
                    int countSpan = tree.CountOccurrences(pattern.AsSpan());
                    Assert.That(countSpan, Is.EqualTo(countStr),
                        $"Seed={seed}: CountOccurrences mismatch for \"{pattern}\"");

                    var posStr = tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                    var posSpan = tree.FindAllOccurrences(pattern.AsSpan()).OrderBy(x => x).ToList();
                    Assert.That(posSpan, Is.EqualTo(posStr),
                        $"Seed={seed}: FindAllOccurrences mismatch for \"{pattern}\"");
                }
            }
        }

        #endregion

        #region P13: Suffixes in lexicographic order

        /// <summary>
        /// GetAllSuffixes must return suffixes in lexicographic (Ordinal) order.
        /// </summary>
        [Test]
        public void P13_GetAllSuffixes_InLexicographicOrder()
        {
            foreach (var (text, seed) in RandomInputs("abcd", 1, 100))
            {
                var tree = SuffixTree.Build(text);
                var suffixes = tree.GetAllSuffixes();

                var expected = Enumerable.Range(0, text.Length)
                    .Select(i => text.Substring(i))
                    .OrderBy(s => s, StringComparer.Ordinal)
                    .ToList();

                Assert.That(suffixes, Is.EqualTo(expected),
                    $"Seed={seed}: suffix order mismatch");
            }
        }

        #endregion

        #region P14: MaxDepth = Text.Length

        [Test]
        public void P14_MaxDepth_EqualsTextLength()
        {
            foreach (var (text, seed) in RandomInputs("ACGT", 1, 300))
            {
                var tree = SuffixTree.Build(text);
                Assert.That(tree.MaxDepth, Is.EqualTo(text.Length),
                    $"Seed={seed}");
            }
        }

        #endregion

        #region P15: Empty pattern → all positions

        /// <summary>
        /// FindAllOccurrences("") returns all positions [0..n-1].
        /// CountOccurrences("") = n.
        /// </summary>
        [Test]
        public void P15_EmptyPattern_ReturnsAllPositions()
        {
            foreach (var (text, seed) in RandomInputs("ab", 1, 50))
            {
                var tree = SuffixTree.Build(text);

                var positions = tree.FindAllOccurrences("").OrderBy(x => x).ToList();
                var expected = Enumerable.Range(0, text.Length).ToList();

                Assert.That(positions, Is.EqualTo(expected),
                    $"Seed={seed}: empty pattern positions mismatch");
                Assert.That(tree.CountOccurrences(""), Is.EqualTo(text.Length),
                    $"Seed={seed}: empty pattern count mismatch");
            }
        }

        #endregion

        #region P16: LCS is truly the longest

        /// <summary>
        /// No common substring of length |LCS|+1 exists.
        /// Verified by brute force on small inputs.
        /// </summary>
        [Test]
        public void P16_LCS_IsMaximal()
        {
            foreach (var (text, seed) in RandomInputs("ab", 2, 50))
            {
                var rng = new Random(seed + 80000);
                string other = GenerateRandom(rng, 2, 50, "ab");

                var tree = SuffixTree.Build(text);
                string lcs = tree.LongestCommonSubstring(other);
                int targetLen = lcs.Length + 1;

                if (targetLen > text.Length || targetLen > other.Length)
                    continue;

                bool foundLonger = false;
                for (int i = 0; i <= text.Length - targetLen; i++)
                {
                    string sub = text.Substring(i, targetLen);
                    if (other.Contains(sub))
                    {
                        foundLonger = true;
                        break;
                    }
                }

                Assert.That(foundLonger, Is.False,
                    $"Seed={seed}: found common substring longer than LCS \"{lcs}\"");
            }
        }

        #endregion
    }
}

