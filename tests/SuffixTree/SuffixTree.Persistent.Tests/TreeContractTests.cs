using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for P3 (suffix lex order), P4 (Contains(null) contract),
/// P6 (LCS correctness — sanity), P12 (Traverse depth semantics).
/// Written RED-first: P3, P4, P12 tests expose existing bugs.
/// </summary>
[TestFixture]
public class TreeContractTests
{
    private string _tempFile = string.Empty;

    [SetUp]
    public void SetUp()
    {
        _tempFile = Path.GetTempFileName();
    }

    [TearDown]
    public void TearDown()
    {
        try { if (File.Exists(_tempFile)) File.Delete(_tempFile); } catch { }
    }

    private PersistentSuffixTree BuildTree(string text)
    {
        var storage = new MappedFileStorageProvider(_tempFile);
        var builder = new PersistentSuffixTreeBuilder(storage);
        var textSource = new StringTextSource(text);
        long rootOffset = builder.Build(textSource);
        if (storage is MappedFileStorageProvider mapped) mapped.TrimToSize();
        return new PersistentSuffixTree(storage, rootOffset, textSource);
    }

    // ─── P3: Suffixes must be in lexicographic order (contract) ──────

    [Test]
    public void EnumerateSuffixes_ReturnsSortedOrder()
    {
        using var tree = BuildTree("banana");
        var suffixes = tree.EnumerateSuffixes().ToList();
        var sorted = suffixes.OrderBy(s => s, StringComparer.Ordinal).ToList();
        Assert.That(suffixes, Is.EqualTo(sorted),
            $"Suffixes not sorted.\nActual:   [{string.Join(", ", suffixes)}]\nExpected: [{string.Join(", ", sorted)}]");
    }

    [Test]
    public void EnumerateSuffixes_ShortText_ReturnsSortedOrder()
    {
        using var tree = BuildTree("abcab");
        var suffixes = tree.EnumerateSuffixes().ToList();
        var sorted = suffixes.OrderBy(s => s, StringComparer.Ordinal).ToList();
        Assert.That(suffixes, Is.EqualTo(sorted));
    }

    // ─── P4: Contains(null) must throw ArgumentNullException ─────────

    [Test]
    public void Contains_NullString_ThrowsArgumentNullException()
    {
        using var tree = BuildTree("test");
        Assert.Throws<ArgumentNullException>(() => tree.Contains((string)null!));
    }

    // ─── P6: LCS correctness sanity (should pass — no RED bug) ──────

    [Test]
    public void LCS_StringOverload_ReturnsCorrectResult()
    {
        using var tree = BuildTree("abcdef");
        string lcs = tree.LongestCommonSubstring("xyzcdeww");
        Assert.That(lcs, Is.EqualTo("cde"));
    }

    // ─── P12: Traverse depth follows ISuffixTree convention ────────
    // Non-persistent convention: VisitNode receives depth BEFORE this node's edge
    // (i.e. parent's cumulative character depth), NOT including current node.

    [Test]
    public void Traverse_DepthMatchesNonPersistentConvention()
    {
        // "aab" → "aab$"
        // Root (depth=0):
        //   "$" leaf → VisitNode depth=0 (parent cumDepth=0)
        //   "a" internal → VisitNode depth=0, fullDepth=0+1=1
        //     "ab$" leaf → VisitNode depth=1 (parent cumDepth=1)
        //     "b$"  leaf → VisitNode depth=1
        //   "b$" leaf → VisitNode depth=0
        // Leaf depths: {0, 0, 1, 1}
        using var tree = BuildTree("aab");
        var visitor = new DepthRecordingVisitor();
        tree.Traverse(visitor);

        var leafDepths = visitor.Nodes
            .Where(n => n.ChildCount == 0)
            .Select(n => n.Depth)
            .OrderBy(d => d)
            .ToList();

        Assert.That(leafDepths, Is.EqualTo(new[] { 0, 0, 1, 1 }),
            $"Leaf depths: [{string.Join(", ", leafDepths)}]");
    }

    [Test]
    public void Traverse_RootDepthIsZero()
    {
        using var tree = BuildTree("abc");
        var visitor = new DepthRecordingVisitor();
        tree.Traverse(visitor);

        Assert.That(visitor.Nodes[0].Depth, Is.EqualTo(0), "Root depth must be 0");
    }

    private class DepthRecordingVisitor : ISuffixTreeVisitor
    {
        public List<(int Start, int End, int LeafCount, int ChildCount, int Depth)> Nodes { get; } = new();

        public void VisitNode(int start, int end, int leafCount, int childCount, int depth)
            => Nodes.Add((start, end, leafCount, childCount, depth));

        public void EnterBranch(int key) { }
        public void ExitBranch() { }
    }
}
