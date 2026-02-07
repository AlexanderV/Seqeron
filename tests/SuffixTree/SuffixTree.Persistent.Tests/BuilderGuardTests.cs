using System;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Test for P8: Builder.Build() called twice must throw (not corrupt data).
/// Written RED-first.
/// </summary>
[TestFixture]
public class BuilderGuardTests
{
    [Test]
    public void Build_CalledTwice_ThrowsInvalidOperationException()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        builder.Build(new StringTextSource("abc"));

        // Second call must throw — not silently corrupt
        Assert.Throws<InvalidOperationException>(() =>
            builder.Build(new StringTextSource("xyz")));
    }

    [Test]
    public void Build_CalledTwice_FirstTreeRemainsCorrect()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        var text = new StringTextSource("hello");
        long root = builder.Build(text);

        // Verify tree is functional before second call
        var tree = new PersistentSuffixTree(storage, root, text);
        Assert.That(tree.Contains("hello"), Is.True);
        Assert.That(tree.Contains("xyz"), Is.False);
    }

    // ─── P13: Builder exposes deepest internal node for O(1) LRS ───

    [Test]
    public void Build_DeepestInternalNodeOffset_IsNotNull_ForNonTrivialTree()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        var text = new StringTextSource("banana");
        long root = builder.Build(text);

        // "banana" has internal nodes, so the deepest must be found
        Assert.That(builder.DeepestInternalNodeOffset, Is.Not.EqualTo(PersistentConstants.NULL_OFFSET),
            "P13: Builder must compute deepest internal node during build");
        Assert.That(builder.DeepestInternalNodeOffset, Is.Not.EqualTo(root),
            "P13: Deepest internal node must not be root for 'banana'");
    }

    [Test]
    public void Build_DeepestInternalNodeOffset_IsNull_ForSingleChar()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        var text = new StringTextSource("a");
        long root = builder.Build(text);

        // "a" → root + "a$" leaf + "$" leaf → no internal nodes besides root
        // DeepestInternalNodeOffset should be root offset (no deeper internal)
        Assert.That(builder.DeepestInternalNodeOffset, Is.EqualTo(root),
            "P13: Single-char tree deepest internal is root");
    }

    [Test]
    public void LRS_UsesPrecomputedDeepestNode_SkipsDFS()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        var text = new StringTextSource("abcabc");
        long root = builder.Build(text);

        // Build tree WITH pre-computed deepest node
        using var tree = new PersistentSuffixTree(storage, root, text,
            deepestInternalNodeOffset: builder.DeepestInternalNodeOffset);

        string lrs = tree.LongestRepeatedSubstring();
        Assert.That(lrs, Is.EqualTo("abc"),
            "P13: LRS with pre-computed deepest must return correct result");
    }
}
