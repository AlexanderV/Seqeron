using System;
using System.IO;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for R2: Import must detect truncated streams with a clear message
/// instead of building a wrong tree and failing later on hash mismatch.
/// </summary>
[TestFixture]
public class ImportTruncationTests
{
    [Test]
    public void Import_TruncatedChars_ThrowsInvalidDataException_WithTruncatedMessage()
    {
        // Build a valid export
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        builder.Build(new StringTextSource("hello world"));
        var tree = PersistentSuffixTree.Load(storage);

        using var ms = new MemoryStream();
        SuffixTreeSerializer.Export(tree, ms);

        // magic(8) + version(4) + 7bit-len(1) = 13 bytes of header,
        // then chars follow. Truncate to 15 bytes — only 1 char of data instead of 11.
        byte[] fullBytes = ms.ToArray();
        var truncated = new MemoryStream(fullBytes, 0, 15, writable: false);

        var target = new HeapStorageProvider();
        var ex = Assert.Throws<InvalidDataException>(() =>
            SuffixTreeSerializer.Import(truncated, target));

        Assert.That(ex!.Message, Does.Contain("Truncated"));
    }

    [Test]
    public void Import_StreamEndsAfterLengthPrefix_ThrowsInvalidDataException()
    {
        // Stream has magic + version + 7-bit length but zero char data
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        builder.Build(new StringTextSource("abcdef"));
        var tree = PersistentSuffixTree.Load(storage);

        using var ms = new MemoryStream();
        SuffixTreeSerializer.Export(tree, ms);

        // magic(8) + version(4) + 7bit-len(1) = 13 bytes — exactly at the length prefix boundary
        byte[] fullBytes = ms.ToArray();
        var truncated = new MemoryStream(fullBytes, 0, 13, writable: false);

        var target = new HeapStorageProvider();
        var ex = Assert.Throws<InvalidDataException>(() =>
            SuffixTreeSerializer.Import(truncated, target));

        Assert.That(ex!.Message, Does.Contain("Truncated"));
    }
}
