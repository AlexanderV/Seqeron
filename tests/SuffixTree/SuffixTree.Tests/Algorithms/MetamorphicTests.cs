using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using CsCheck;

namespace SuffixTree.Tests.Algorithms
{
    [TestFixture]
    [Category("Algorithms")]
    public class MetamorphicTests
    {
        private static Gen<string> DnaStringGen() =>
            Gen.Select(Gen.Int[0, 3], i => "ACGT"[i]).Array[1, 120].Select(chars => new string(chars));

        private static Gen<(string Text, string Pattern)> ExistingPatternGen() =>
            from text in DnaStringGen()
            from start in Gen.Int[0, text.Length - 1]
            from len in Gen.Int[1, text.Length - start]
            select (text, text.Substring(start, len));

        // MR1: Build(T) -> Build(T + T) => FindAll(p) in T+T doubles the number of positions compared to T
        [Test]
        public void MR1_FindAll_ContainsMappedPositions_OnDoubledText()
        {
            ExistingPatternGen().Sample(t => {
                var (text, pattern) = t;
                var tree1 = SuffixTree.Build(text);
                var tree2 = SuffixTree.Build(text + text);

                var occ1 = tree1.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                var occ2 = tree2.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                var occ2Set = new HashSet<int>(occ2);

                // Every occurrence in T must reappear in both halves of T+T.
                foreach (int pos in occ1)
                {
                    Assert.That(occ2Set.Contains(pos), Is.True, $"Missing mapped position at {pos}");
                    Assert.That(occ2Set.Contains(pos + text.Length), Is.True, $"Missing shifted mapped position at {pos + text.Length}");
                }

                Assert.That(occ2.Count, Is.GreaterThanOrEqualTo(occ1.Count * 2));
            });
        }

        // MR2: Contains(T, p) -> Build(T + "X" + T) => All positions from T are present. 
        // We test this by checking that if p is in T, it's also in T+"X"+T.
        [Test]
        public void MR2_Contains_PreservedAndAbsentPatternStaysAbsent()
        {
            ExistingPatternGen().Sample(t => {
                var (text, existingPattern) = t;
                var treeT = SuffixTree.Build(text);
                var treeTXT = SuffixTree.Build(text + "X" + text);

                Assert.That(treeT.Contains(existingPattern), Is.True);
                Assert.That(treeTXT.Contains(existingPattern), Is.True);

                // '$' never appears in DNA input and also not in inserted separator 'X'.
                const string missingPattern = "$";
                Assert.That(treeT.Contains(missingPattern), Is.False);
                Assert.That(treeTXT.Contains(missingPattern), Is.False);
            });
        }

        // MR3: LCS(A, B) => LCS(B, A) has the same length
        [Test]
        public void MR3_LCS_Symmetry()
        {
            Gen.Select(DnaStringGen(), DnaStringGen()).Sample(t => {
                var (a, b) = t;
                var treeA = SuffixTree.Build(a);
                var treeB = SuffixTree.Build(b);

                var lcsAB = treeA.LongestCommonSubstring(b);
                var lcsBA = treeB.LongestCommonSubstring(a);

                Assert.That(lcsAB.Length, Is.EqualTo(lcsBA.Length));
            });
        }

        // MR4: LRS(T) -> LRS(T + T) => |LRS(T+T)| >= |T| (because T itself is repeated)
        [Test]
        public void MR4_DoubledText_MakesOriginalTextRepeated()
        {
            DnaStringGen().Sample(text => {
                var doubled = text + text;
                var tree = SuffixTree.Build(doubled);

                Assert.That(tree.CountOccurrences(text), Is.GreaterThanOrEqualTo(2));

                var lrs = tree.LongestRepeatedSubstring();
                Assert.That(lrs.Length, Is.GreaterThanOrEqualTo(text.Length));
            });
        }

        // MR5: FindAll(T, p) -> FindAll(T, p[0..-1]) => |FindAll(p)| <= |FindAll(p[0..-1])|
        [Test]
        public void MR5_FindAll_PatternPositionsAreSubsetOfPrefixPositions()
        {
            ExistingPatternGen().Sample(t => {
                var (text, pattern) = t;
                if (pattern.Length < 2) return;

                var tree = SuffixTree.Build(text);
                var prefix = pattern.Substring(0, pattern.Length - 1);

                var occPattern = tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
                var occPrefix = tree.FindAllOccurrences(prefix).OrderBy(x => x).ToHashSet();

                foreach (int pos in occPattern)
                {
                    Assert.That(occPrefix.Contains(pos), Is.True, $"Prefix misses pattern position {pos}");
                }
            });
        }

        // MR6: NodeCount(T) -> NodeCount(T + "$") => NodeCount changes predictably 
        // Appending a unique terminal symbol exactly forces all implicit suffixes to become explicit leaves,
        // which might increase the node count, but it structurally bounds it.
        [Test]
        public void MR6_NodeCount_PredictableChangeWithTerminal()
        {
            DnaStringGen().Sample(text => {
                var treeT = SuffixTree.Build(text);
                // Assume '$' never appears in DnaStringGen
                var treeTerminal = SuffixTree.Build(text + "$");

                // Adding a unique terminal symbol forces all previously implicit suffixes to become explicit leaves.
                // The new tree must have exactly |text| + 1 leaves.
                Assert.That(treeTerminal.LeafCount, Is.EqualTo(text.Length + 1));
                Assert.That(treeTerminal.CountOccurrences("$"), Is.EqualTo(1));
                Assert.That(treeT.CountOccurrences("$"), Is.EqualTo(0));

                // The total node count bounds still apply strictly.
                Assert.That(treeTerminal.NodeCount, Is.LessThanOrEqualTo(2 * (text.Length + 1) + 1));
            });
        }
    }
}

