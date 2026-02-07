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
}
