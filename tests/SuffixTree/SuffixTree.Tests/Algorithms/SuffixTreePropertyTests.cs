using System;
using System.Linq;
using NUnit.Framework;
using CsCheck;

namespace SuffixTree.Tests.Algorithms
{
    [TestFixture]
    public class SuffixTreePropertyTests
    {
        // A simple DNA string generator for CsCheck: length 1 to 500, characters A, C, G, T
        private static Gen<string> DnaStringGen() => 
            Gen.Select(Gen.Int[0, 3], i => "ACGT"[i]).Array[1, 500].Select(chars => new string(chars));

        // 1. Leaf count: LeafCount = |T|
        [Test]
        public void LeafCount_Equals_TextLength()
        {
            DnaStringGen().Sample(text => {
                var tree = SuffixTree.Build(text);
                Assert.That(tree.LeafCount, Is.EqualTo(text.Length));
            });
        }

        // 2. Node bound: NodeCount <= 2|T| + 1
        [Test]
        public void NodeBound_WithinTheoreticalLimit()
        {
            DnaStringGen().Sample(text => {
                var tree = SuffixTree.Build(text);
                Assert.That(tree.NodeCount, Is.LessThanOrEqualTo(2 * text.Length + 1));
            });
        }

        // 3. Suffix completeness: For all i: Contains(T[i..]) = true
        [Test]
        public void SuffixCompleteness_AllSuffixesFound()
        {
            DnaStringGen().Sample(text => {
                var tree = SuffixTree.Build(text);
                
                for (int i = 0; i < text.Length; i++)
                {
                    Assert.That(tree.Contains(text.Substring(i)), Is.True);
                }
            });
        }

        // 4. Count = positions: Count(p) = |FindAll(p)|
        [Test]
        public void CountOccurrences_Equals_PositionsCount()
        {
            Gen.Select(DnaStringGen(), DnaStringGen()).Sample(t => {
                var (text, pattern) = t;
                var tree = SuffixTree.Build(text);
                int count = tree.CountOccurrences(pattern);
                var positions = tree.FindAllOccurrences(pattern);
                
                Assert.That(positions.Count, Is.EqualTo(count));
            });
        }

        // 5. Empty pattern: FindAll("") = [0..|T|-1]
        [Test]
        public void EmptyPattern_FindsAllPositions()
        {
            DnaStringGen().Sample(text => {
                var tree = SuffixTree.Build(text);
                
                var positions = tree.FindAllOccurrences("");
                Assert.That(positions.Count, Is.EqualTo(text.Length));
                Assert.That(positions.SequenceEqual(Enumerable.Range(0, text.Length)), Is.True);
            });
        }

        // 6. LCS symmetry: |LCS(A,B)| == |LCS(B,A)|
        [Test]
        public void LCS_Is_Symmetric()
        {
            Gen.Select(DnaStringGen(), DnaStringGen()).Sample(t => {
                var (a, b) = t;
                var treeA = SuffixTree.Build(a);
                var lcsAB = treeA.LongestCommonSubstring(b);
                
                var treeB = SuffixTree.Build(b);
                var lcsBA = treeB.LongestCommonSubstring(a);
                
                Assert.That(lcsAB.Length, Is.EqualTo(lcsBA.Length));
            });
        }

        // 7. LRS existence: |LRS| == max substring with Count >= 2
        [Test]
        public void LRS_IfExists_OccursAtLeastTwice()
        {
            DnaStringGen().Sample(text => {
                var tree = SuffixTree.Build(text);
                
                var lrs = tree.LongestRepeatedSubstring();
                if (lrs.Length > 0)
                {
                    Assert.That(tree.CountOccurrences(lrs), Is.GreaterThanOrEqualTo(2));
                }
            });
        }

        // 8. Idempotent Contains: Contains(p) = (Count(p) > 0)
        [Test]
        public void Contains_Equals_CountGreaterThanZero()
        {
            Gen.Select(DnaStringGen(), DnaStringGen()).Sample(t => {
                var (text, pattern) = t;
                var tree = SuffixTree.Build(text);
                bool contains = tree.Contains(pattern);
                int count = tree.CountOccurrences(pattern);
                
                Assert.That(contains, Is.EqualTo(count > 0));
            });
        }
    }
}
