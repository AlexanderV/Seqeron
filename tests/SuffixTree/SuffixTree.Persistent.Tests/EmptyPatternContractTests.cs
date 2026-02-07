using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for Q1 (FindAllOccurrences("") contract) and Q2 (CountOccurrences("") contract).
/// ISuffixTree contract: empty pattern → all positions / text length.
/// Written RED-first.
/// </summary>
[TestFixture]
public class EmptyPatternContractTests
{
    private string _tempFile = string.Empty;

    [SetUp]
    public void SetUp() => _tempFile = Path.GetTempFileName();

    [TearDown]
    public void TearDown()
    {
        try { if (File.Exists(_tempFile)) File.Delete(_tempFile); } catch { }
    }

    private PersistentSuffixTree BuildTree(string text)
    {
        var storage = new MappedFileStorageProvider(_tempFile);
        var builder = new PersistentSuffixTreeBuilder(storage);
        builder.Build(new StringTextSource(text));
        return PersistentSuffixTree.Load(storage);
    }

    // ─── Q1: FindAllOccurrences("") must return all positions ────────

    [Test]
    public void FindAllOccurrences_EmptyString_ReturnsAllPositions()
    {
        using var tree = BuildTree("abcde");
        var result = tree.FindAllOccurrences("");
        var expected = Enumerable.Range(0, 5).ToList(); // [0,1,2,3,4]
        Assert.That(result.OrderBy(x => x).ToList(), Is.EqualTo(expected));
    }

    [Test]
    public void FindAllOccurrences_EmptySpan_ReturnsAllPositions()
    {
        using var tree = BuildTree("abc");
        var result = tree.FindAllOccurrences(ReadOnlySpan<char>.Empty);
        var expected = Enumerable.Range(0, 3).ToList();
        Assert.That(result.OrderBy(x => x).ToList(), Is.EqualTo(expected));
    }

    // ─── Q2: CountOccurrences("") must return text length ───────────

    [Test]
    public void CountOccurrences_EmptyString_ReturnsTextLength()
    {
        using var tree = BuildTree("hello");
        Assert.That(tree.CountOccurrences(""), Is.EqualTo(5));
    }

    [Test]
    public void CountOccurrences_EmptySpan_ReturnsTextLength()
    {
        using var tree = BuildTree("test");
        Assert.That(tree.CountOccurrences(ReadOnlySpan<char>.Empty), Is.EqualTo(4));
    }
}
