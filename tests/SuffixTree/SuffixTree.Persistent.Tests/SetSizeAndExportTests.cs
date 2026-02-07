using System;
using System.IO;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for Q8 (SetSize validation) and Q10 (Export text serialization).
/// </summary>
[TestFixture]
public class SetSizeAndExportTests
{
    private string _tempFile = string.Empty;

    [SetUp]
    public void SetUp() => _tempFile = Path.GetTempFileName();

    [TearDown]
    public void TearDown()
    {
        try { if (File.Exists(_tempFile)) File.Delete(_tempFile); } catch { }
    }

    // ─── Q8: SetSize must validate input ────────────────────────────

    [Test]
    public void SetSize_AfterDispose_ThrowsObjectDisposedException()
    {
        var p = new MappedFileStorageProvider(_tempFile);
        p.Dispose();
        Assert.Throws<ObjectDisposedException>(() => p.SetSize(100));
    }

    [Test]
    public void SetSize_NegativeSize_ThrowsArgumentOutOfRangeException()
    {
        using var p = new MappedFileStorageProvider(_tempFile);
        Assert.Throws<ArgumentOutOfRangeException>(() => p.SetSize(-1));
    }

    [Test]
    public void SetSize_ValidSize_UpdatesPosition()
    {
        using var p = new MappedFileStorageProvider(_tempFile);
        p.SetSize(100);
        Assert.That(p.Size, Is.EqualTo(100));
    }

    // ─── Q10: Export round-trip correctness (sanity) ────────────────

    [Test]
    public void Export_LargeText_RoundTripsCorrectly()
    {
        // Text > 4096 chars to test chunked export path
        string longText = new string('X', 5000) + "FIND_ME" + new string('Y', 5000);
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        builder.Build(new StringTextSource(longText));
        var tree = PersistentSuffixTree.Load(storage);

        using var ms = new MemoryStream();
        SuffixTreeSerializer.Export(tree, ms);
        ms.Position = 0;

        var target = new HeapStorageProvider();
        var imported = SuffixTreeSerializer.Import(ms, target);

        Assert.That(imported.Text.ToString(), Is.EqualTo(longText));
        Assert.That(imported.Contains("FIND_ME"), Is.True);
    }
}
