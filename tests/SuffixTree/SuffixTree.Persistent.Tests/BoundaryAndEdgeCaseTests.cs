using System;
using System.Collections.Generic;
using System.IO;
using NUnit.Framework;
using SuffixTree;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class BoundaryAndEdgeCaseTests
    {
        private string _tempFile = "boundary_test.tree";

        [TearDown]
        public void Cleanup()
        {
            if (File.Exists(_tempFile)) File.Delete(_tempFile);
        }

        [Test]
        public void Mmf_FrequentResizing_MaintainsIntegrity()
        {
            // Start with a tiny capacity to force many resizes
            string text = "the quick brown fox jumps over the lazy dog";

            using (var storage = new MappedFileStorageProvider(_tempFile, initialCapacity: 64))
            {
                var builder = new PersistentSuffixTreeBuilder(storage);
                builder.Build(new StringTextSource(text));

                var tree = PersistentSuffixTree.Load(storage);
                Assert.That(tree.Text.ToString(), Is.EqualTo(text));
                Assert.That(tree.Contains("fox jumps"), Is.True);
                Assert.That(tree.Contains("lazy dog"), Is.True);
                Assert.That(tree.CountOccurrences("the"), Is.EqualTo(2));
            }
        }

        [Test]
        public void StructuralInvariant_EdgeCountMatchesNodeCount()
        {
            string text = "abracadabra";
            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                int totalEdges = 0;

                tree.Traverse(new EdgeCounterVisitor((count) => totalEdges += count));

                // In any tree, Number of Edges = Number of Nodes - 1
                Assert.That(totalEdges, Is.EqualTo(tree.NodeCount - 1), "Structural invariant |E| = |V| - 1 failed");
            }
        }

        [Test]
        public void Pathological_RepetitivePattern_Works()
        {
            // abcabcabc...
            string pattern = "abc";
            string text = "";
            for (int i = 0; i < 100; i++) text += pattern;

            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                Assert.That(tree.Contains("abcabc"), Is.True);
                Assert.That(tree.CountOccurrences("abc"), Is.EqualTo(100));

                // For a periodic string of length N and period P, LRS is N-P
                string expectedLrs = text.Substring(0, text.Length - 3);
                Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(expectedLrs));
            }
        }

        [Test]
        public void Encapsulation_UserStringContainsTerminator()
        {
            // If the user text contains the character we use as terminator, we should handle it.
            // Our implementation uses '$' (char 36).
            char terminator = '$';
            string text = $"abc{terminator}def";

            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;

                // The tree should treat the internal terminator as part of the text if it's there.
                // However, Ukkonen ADDS another terminator. 
                // We need to ensure searching for the string works even if it has '$'.
                Assert.That(tree.Contains($"abc{terminator}d"), Is.True);
                Assert.That(tree.Text.ToString(), Is.EqualTo(text));
            }
        }

        [Test]
        public void EdgeCase_OnlyTerminatorString()
        {
            char terminator = '$';
            string text = terminator.ToString();

            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                Assert.That(tree.Text.ToString(), Is.EqualTo(text));
                Assert.That(tree.LeafCount, Is.EqualTo(1)); // The character is there, but ISuffixTree.LeafCount subtracts 1? 
                // Wait, if text is "$", and we add "$", we have "$$". 
                // LeafCount 2, ISuffixTree.LeafCount returns 1. Correct.
            }
        }

        [Test]
        public void Stress_RandomSmallAlphabet_DenseTree()
        {
            var random = new Random(9000);
            // Small alphabet (a, b) creates very deep and dense structures
            string text = new string(Enumerable.Range(0, 1000).Select(_ => (char)random.Next('a', 'c')).ToArray());

            using (var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var tree = (ISuffixTree)st!;
                // Just check that it doesn't crash and satisfies base invariants
                Assert.That(tree.NodeCount, Is.GreaterThan(text.Length));
                Assert.That(tree.LeafCount, Is.EqualTo(text.Length));

                string lrs = tree.LongestRepeatedSubstring();
                Assert.That(tree.CountOccurrences(lrs), Is.GreaterThanOrEqualTo(2));
            }
        }

        private class EdgeCounterVisitor : ISuffixTreeVisitor
        {
            private readonly Action<int> _onNode;
            public EdgeCounterVisitor(Action<int> onNode) => _onNode = onNode;

            public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
            {
                _onNode(childCount);
            }

            public void EnterBranch(int key) { }
            public void ExitBranch() { }
        }
    }
}
