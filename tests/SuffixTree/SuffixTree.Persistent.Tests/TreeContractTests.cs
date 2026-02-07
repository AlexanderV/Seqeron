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

    // ─── P12: Traverse depth must be character-based ─────────────────

    [Test]
    public void Traverse_LeafDepthEqualsSuffixLength()
    {
        // "aab" → "aab$" has multi-char edges, so edge-count ≠ character depth.
        // Suffix tree for "aab$":
        //   Root → "$" (leaf, depth 1) 
        //        → "a" (internal, depth 1) → "ab$" (leaf, depth 4)
        //                                  → "b$" (leaf, depth 3)
        //        → "b$" (leaf, depth 2)
        // Leaf character-depths must be {1, 2, 3, 4} = suffix lengths.
        using var tree = BuildTree("aab");
        var visitor = new DepthRecordingVisitor();
        tree.Traverse(visitor);

        var leafDepths = visitor.Nodes
            .Where(n => n.ChildCount == 0)
            .Select(n => n.Depth)
            .OrderBy(d => d)
            .ToList();

        // With the bug, depths are {1, 1, 2, 2} instead of {1, 2, 3, 4}
        Assert.That(leafDepths, Is.EqualTo(new[] { 1, 2, 3, 4 }),
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
