using System;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for R5 (Factory.Create null text), R6 (Builder.Build null text),
/// R8 (PersistentSuffixTree.Load null storage).
/// All three currently throw NullReferenceException; should throw ArgumentNullException.
/// </summary>
[TestFixture]
public class NullGuardTests
{
    // ─── R5: Factory.Create(null) ───────────────────────────────────

    [Test]
    public void Factory_Create_NullText_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            PersistentSuffixTreeFactory.Create(null!));
    }

    // ─── R6: Builder.Build(null) ────────────────────────────────────

    [Test]
    public void Builder_Build_NullText_ThrowsArgumentNullException()
    {
        using var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        Assert.Throws<ArgumentNullException>(() => builder.Build(null!));
    }

    // ─── R8: PersistentSuffixTree.Load(null) ────────────────────────

    [Test]
    public void Load_NullStorage_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            PersistentSuffixTree.Load(null!));
    }
}
