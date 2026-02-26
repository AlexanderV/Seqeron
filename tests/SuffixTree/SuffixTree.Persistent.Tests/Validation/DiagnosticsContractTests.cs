using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests.Validation;

[TestFixture]
[Category("Validation")]
public class DiagnosticsContractTests
{
    [Test]
    public void PrintTree_EmptyTree_HasOnlyRootNodeLine()
    {
        using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(new StringTextSource(string.Empty));
        var parsed = ParsePrintTree(tree.PrintTree());

        Assert.Multiple(() =>
        {
            Assert.That(parsed.ContentLength, Is.EqualTo(0));
            Assert.That(parsed.NodeLines.Count, Is.EqualTo(1));
            Assert.That(parsed.NodeLines[0].Depth, Is.EqualTo(0));
            Assert.That(parsed.NodeLines[0].Label, Is.EqualTo("ROOT"));
        });
    }

    [TestCase("ab")]
    [TestCase("banana")]
    [TestCase("mississippi")]
    public void PrintTree_NodeLinesMatchNodeCount_AndDepthIndentContract(string text)
    {
        using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(new StringTextSource(text));
        var parsed = ParsePrintTree(tree.PrintTree());

        Assert.Multiple(() =>
        {
            Assert.That(parsed.ContentLength, Is.EqualTo(text.Length));
            Assert.That(parsed.NodeLines.Count, Is.EqualTo(tree.NodeCount),
                "Each tree node must be printed exactly once");
            Assert.That(parsed.NodeLines[0].Label, Is.EqualTo("ROOT"));
            Assert.That(parsed.NodeLines[0].Depth, Is.EqualTo(0));
            Assert.That(parsed.NodeLines.All(line => line.IndentSpaces == line.Depth * 2), Is.True,
                "Indentation must be exactly 2 spaces per depth level");
        });
    }

    [TestCase("ab")]
    [TestCase("banana")]
    [TestCase("mississippi")]
    public void PrintTree_NonEmptyTree_LeafLineCountMatchesLeafCountPlusTerminator(string text)
    {
        using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(new StringTextSource(text));
        var parsed = ParsePrintTree(tree.PrintTree());
        int leafLines = parsed.NodeLines.Count(line => line.IsLeaf);

        Assert.That(leafLines, Is.EqualTo(tree.LeafCount + 1),
            "Printed leaves include explicit terminator leaf in addition to API LeafCount");
    }

    [Test]
    public void PrintTree_FirstLevelChildren_AreSortedByLeadingSymbol()
    {
        using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(new StringTextSource("banana"));
        var parsed = ParsePrintTree(tree.PrintTree());
        var depth1Leading = parsed.NodeLines
            .Where(line => line.Depth == 1)
            .Select(line => line.Label[0])
            .ToArray();

        var sorted = depth1Leading.OrderBy(ch => ch).ToArray();
        Assert.That(depth1Leading, Is.EqualTo(sorted),
            "Root children must be printed in deterministic sorted order");
    }

    [Test]
    public void Traverse_EnterAndExitBranches_AreBalanced()
    {
        using var tree = (PersistentSuffixTree)PersistentSuffixTreeFactory.Create(new StringTextSource("banana"));
        var visitor = new CountingVisitor();
        tree.Traverse(visitor);

        Assert.Multiple(() =>
        {
            Assert.That(visitor.NodeCount, Is.EqualTo(tree.NodeCount));
            Assert.That(visitor.EnterCount, Is.EqualTo(visitor.ExitCount));
        });
    }

    private static ParsedPrintTree ParsePrintTree(string output)
    {
        var lines = output
            .Split('\n', StringSplitOptions.RemoveEmptyEntries)
            .Select(static line => line.TrimEnd('\r'))
            .ToArray();
        Assert.That(lines.Length, Is.GreaterThanOrEqualTo(2), "PrintTree output must contain header and root");
        Assert.That(lines[0], Does.StartWith("Content length: "));

        int contentLength = int.Parse(lines[0]["Content length: ".Length..], System.Globalization.CultureInfo.InvariantCulture);
        var nodeLines = new List<PrintedNodeLine>();

        foreach (string line in lines.Skip(1))
        {
            int colon = line.IndexOf(':');
            if (colon <= 0)
                continue;

            string depthPart = line[..colon].Trim();
            if (!int.TryParse(depthPart, out int depth))
                continue;

            int indentSpaces = line.TakeWhile(ch => ch == ' ').Count();
            string content = line[(colon + 1)..].TrimStart();

            bool isLeaf = content.EndsWith(" (Leaf)", StringComparison.Ordinal);
            if (isLeaf)
                content = content[..^" (Leaf)".Length];

            int linkMarker = content.IndexOf(" -> [", StringComparison.Ordinal);
            if (linkMarker >= 0)
                content = content[..linkMarker];

            nodeLines.Add(new PrintedNodeLine(depth, indentSpaces, content, isLeaf));
        }

        return new ParsedPrintTree(contentLength, nodeLines);
    }

    private sealed class CountingVisitor : ISuffixTreeVisitor
    {
        public int NodeCount { get; private set; }
        public int EnterCount { get; private set; }
        public int ExitCount { get; private set; }

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => NodeCount++;
        public void EnterBranch(int key) => EnterCount++;
        public void ExitBranch() => ExitCount++;
    }

    private sealed record PrintedNodeLine(int Depth, int IndentSpaces, string Label, bool IsLeaf);
    private sealed record ParsedPrintTree(int ContentLength, IReadOnlyList<PrintedNodeLine> NodeLines);
}
