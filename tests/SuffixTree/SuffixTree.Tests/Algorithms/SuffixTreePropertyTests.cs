using System;
using System.Linq;
using NUnit.Framework;
using FsCheck;

namespace SuffixTree.Tests.Algorithms
{
    public class NonEmptyDnaString
    {
        public string Value { get; }
        public NonEmptyDnaString(string value) => Value = value;
        public override string ToString() => Value;
    }

    [TestFixture]
    public class SuffixTreePropertyTests
    {
        private Arbitrary<NonEmptyDnaString> _dnaStrings;

        [SetUp]
        public void SetUp()
        {
            var gen = Gen.Elements('A', 'C', 'G', 'T')
                         .ArrayOf()
                         .Where(a => a.Length > 0 && a.Length <= 500)
                         .Select(a => new NonEmptyDnaString(new string(a)));
            _dnaStrings = Arb.From(gen);
        }

        // 1. Leaf count: LeafCount = |T|
        [Test]
        public void LeafCount_Equals_TextLength()
        {
            Prop.ForAll(_dnaStrings, input => {
                var text = input.Value;
                var tree = SuffixTree.Build(text);
                return tree.LeafCount == text.Length;
            }).QuickCheckThrowOnFailure();
        }

        // 2. Node bound: NodeCount <= 2|T| + 1
        [Test]
        public void NodeBound_WithinTheoreticalLimit()
        {
            Prop.ForAll(_dnaStrings, input => {
                var text = input.Value;
                var tree = SuffixTree.Build(text);
                return tree.NodeCount <= 2 * text.Length + 1;
            }).QuickCheckThrowOnFailure();
        }

        // 3. Suffix completeness: For all i: Contains(T[i..]) = true
        [Test]
        public void SuffixCompleteness_AllSuffixesFound()
        {
            Prop.ForAll(_dnaStrings, input => {
                var text = input.Value;
                var tree = SuffixTree.Build(text);
                
                for (int i = 0; i < text.Length; i++)
                {
                    if (!tree.Contains(text.Substring(i)))
                    {
                        return false;
                    }
                }
                return true;
            }).QuickCheckThrowOnFailure();
        }

        // 4. Count = positions: Count(p) = |FindAll(p)|
        [Test]
        public void CountOccurrences_Equals_PositionsCount()
        {
            Prop.ForAll(_dnaStrings, _dnaStrings, (input, patternInput) => {
                var text = input.Value;
                var pattern = patternInput.Value;
                
                var tree = SuffixTree.Build(text);
                int count = tree.CountOccurrences(pattern);
                var positions = tree.FindAllOccurrences(pattern);
                
                return count == positions.Count;
            }).QuickCheckThrowOnFailure();
        }

        // 5. Empty pattern: FindAll("") = [0..|T|-1]
        [Test]
        public void EmptyPattern_FindsAllPositions()
        {
            Prop.ForAll(_dnaStrings, input => {
                var text = input.Value;
                var tree = SuffixTree.Build(text);
                
                var positions = tree.FindAllOccurrences("");
                return positions.Count == text.Length && positions.SequenceEqual(Enumerable.Range(0, text.Length));
            }).QuickCheckThrowOnFailure();
        }

        // 6. LCS symmetry: |LCS(A,B)| == |LCS(B,A)|
        [Test]
        public void LCS_Is_Symmetric()
        {
            Prop.ForAll(_dnaStrings, _dnaStrings, (inputA, inputB) => {
                var a = inputA.Value;
                var b = inputB.Value;
                
                var treeA = SuffixTree.Build(a);
                var lcsAB = treeA.LongestCommonSubstring(b);
                
                var treeB = SuffixTree.Build(b);
                var lcsBA = treeB.LongestCommonSubstring(a);
                
                return lcsAB.Length == lcsBA.Length;
            }).QuickCheckThrowOnFailure();
        }

        // 7. LRS existence: |LRS| == max substring with Count >= 2
        [Test]
        public void LRS_IfExists_OccursAtLeastTwice()
        {
            Prop.ForAll(_dnaStrings, input => {
                var text = input.Value;
                var tree = SuffixTree.Build(text);
                
                var lrs = tree.LongestRepeatedSubstring();
                if (lrs.Length > 0)
                {
                    return tree.CountOccurrences(lrs) >= 2;
                }
                return true;
            }).QuickCheckThrowOnFailure();
        }

        // 8. Idempotent Contains: Contains(p) = (Count(p) > 0)
        [Test]
        public void Contains_Equals_CountGreaterThanZero()
        {
            Prop.ForAll(_dnaStrings, _dnaStrings, (input, patternInput) => {
                var text = input.Value;
                var pattern = patternInput.Value;
                
                var tree = SuffixTree.Build(text);
                bool contains = tree.Contains(pattern);
                int count = tree.CountOccurrences(pattern);
                
                return contains == (count > 0);
            }).QuickCheckThrowOnFailure();
        }
    }
}
