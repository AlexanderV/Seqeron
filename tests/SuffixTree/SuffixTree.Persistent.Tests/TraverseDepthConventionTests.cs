using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Test for Q3: Traverse depth convention must match non-persistent tree.
/// VisitNode depth = cumulative character depth BEFORE this node's edge
/// (i.e., parent's full depth), NOT including current edge length.
/// </summary>
[TestFixture]
public class TraverseDepthConventionTests
{
    private string _tempFile = string.Empty;

    [SetUp]
    public void SetUp() => _tempFile = Path.GetTempFileName();

    [TearDown]
    public void TearDown()
    {
        try { if (File.Exists(_tempFile)) File.Delete(_tempFile); } catch { }
    }

    [Test]
    public void Traverse_ChildDepthIsParentCumulativeDepth()
    {
        // Build tree for "aab" → after terminal: "aab$"
        // Root (depth=0) has children:
        //   '$' leaf edge length=1 → visitor depth must be 0 (parent=root, cumDepth=0)
        //   'a' internal edge length=1 → visitor depth must be 0
        //     'ab$' leaf edge length=3 → visitor depth must be 1 (parent 'a' cumDepth=0+1=1)
        //     'b$' leaf edge length=2 → visitor depth must be 1
        //   'b' leaf edge length=2 → visitor depth must be 0
        //
        // Non-persistent convention: VisitNode gets depth BEFORE the node = parent cumulative.
        // So root=0, all root children=0, children of 'a' node=1.

        var storage = new MappedFileStorageProvider(_tempFile);
        var builder = new PersistentSuffixTreeBuilder(storage);
        builder.Build(new StringTextSource("aab"));
        using var tree = PersistentSuffixTree.Load(storage);

        var visitor = new DepthCollector();
        tree.Traverse(visitor);

        // Root gets depth=0
        Assert.That(visitor.Nodes[0].Depth, Is.EqualTo(0), "Root depth");

        // All direct children of root should get depth=0 (parent=root, cumDepth=0)
        // Nodes at index 1+ that are direct root children were entered via EnterBranch from root
        // We track branch depth to verify
        foreach (var child in visitor.RootChildren)
        {
            Assert.That(child.Depth, Is.EqualTo(0),
                $"Root child (start={child.Start}) should have depth=0 (parent cumulative), got {child.Depth}");
        }
    }

    private class DepthCollector : ISuffixTreeVisitor
    {
        public List<(int Start, int End, int LeafCount, int ChildCount, int Depth)> Nodes { get; } = new();
        public List<(int Start, int End, int Depth)> RootChildren { get; } = new();
        private int _branchDepth = 0;

        public void VisitNode(int start, int end, int leafCount, int childCount, int depth)
        {
            Nodes.Add((start, end, leafCount, childCount, depth));
            if (_branchDepth == 1) // direct child of root
                RootChildren.Add((start, end, depth));
        }

        public void EnterBranch(int key) => _branchDepth++;
        public void ExitBranch() => _branchDepth--;
    }
}
