using System;
using System.Linq;
using NUnit.Framework;
using CsCheck;

namespace SuffixTree.Tests.Algorithms
{
    [TestFixture]
    public class MetamorphicTests
    {
        private static Gen<string> DnaStringGen() => 
            Gen.Select(Gen.Int[0, 3], i => "ACGT"[i]).Array[1, 100].Select(chars => new string(chars));

        // MR1: Build(T) -> Build(T + T) => FindAll(p) in T+T doubles the number of positions compared to T
        [Test]
        public void MR1_FindAll_DoublesPositions_OnDoubledText()
        {
            Gen.Select(DnaStringGen(), DnaStringGen()).Sample(t => {
                var (text, pattern) = t;
                var tree1 = SuffixTree.Build(text);
                var tree2 = SuffixTree.Build(text + text);

                var occ1 = tree1.FindAllOccurrences(pattern);
                var occ2 = tree2.FindAllOccurrences(pattern);

                int expectedCount = occ1.Count * 2;
                // If the pattern itself overlaps the boundary of T+T, it might appear MORE than 2 times.
                // But it will ALWAYS appear AT LEAST 2 times the original count (once in the first half, once in the second)
                Assert.That(occ2.Count, Is.GreaterThanOrEqualTo(expectedCount));
            });
        }

        // MR2: Contains(T, p) -> Build(T + "X" + T) => All positions from T are present. 
        // We test this by checking that if p is in T, it's also in T+"X"+T.
        [Test]
        public void MR2_Contains_PreservedInConcatenation()
        {
            Gen.Select(DnaStringGen(), DnaStringGen()).Sample(t => {
                var (text, pattern) = t;
                var treeT = SuffixTree.Build(text);
                
                if (treeT.Contains(pattern))
                {
                    var treeTXT = SuffixTree.Build(text + "X" + text);
                    Assert.That(treeTXT.Contains(pattern), Is.True);
                }
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
        public void MR4_LRS_OfDoubledText_IsAtLeastLengthOfText()
        {
            DnaStringGen().Sample(text => {
                var tree = SuffixTree.Build(text + text);
                var lrs = tree.LongestRepeatedSubstring();

                Assert.That(lrs.Length, Is.GreaterThanOrEqualTo(text.Length));
            });
        }

        // MR5: FindAll(T, p) -> FindAll(T, p[0..-1]) => |FindAll(p)| <= |FindAll(p[0..-1])|
        [Test]
        public void MR5_FindAll_SubsetPrefix_HasMoreOrEqualOccurrences()
        {
            Gen.Select(DnaStringGen(), DnaStringGen()).Sample(t => {
                var (text, pattern) = t;
                // Ensure pattern is at least length 2 to take a substring
                if (pattern.Length >= 2)
                {
                    var tree = SuffixTree.Build(text);
                    var prefix = pattern.Substring(0, pattern.Length - 1);

                    var occPattern = tree.FindAllOccurrences(pattern);
                    var occPrefix = tree.FindAllOccurrences(prefix);

                    Assert.That(occPrefix.Count, Is.GreaterThanOrEqualTo(occPattern.Count));
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
                
                // The total node count bounds still apply strictly.
                Assert.That(treeTerminal.NodeCount, Is.LessThanOrEqualTo(2 * (text.Length + 1) + 1));
            });
        }
    }
}
