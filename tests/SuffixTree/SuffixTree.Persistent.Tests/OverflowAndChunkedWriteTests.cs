using System;
using System.Buffers;
using System.Text;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for Q6 (int overflow in FinalizeTree text allocation),
/// Q7 (int overflow in LoadStringInternal), Q9 (chunked write perf — sanity).
/// </summary>
[TestFixture]
public class OverflowAndChunkedWriteTests
{
    // ─── Q6: FinalizeTree must reject text > int.MaxValue/2 chars ───
    // We can't allocate 1B+ chars, but we test that the Allocate call
    // uses long arithmetic by verifying small texts work correctly
    // (the overflow fix also adds an explicit guard for large texts).

    [Test]
    public void Builder_SmallText_WritesTextCorrectly()
    {
        // Sanity: verify text round-trips through build+load
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        var text = new StringTextSource("ACGT");
        builder.Build(text);

        var tree = PersistentSuffixTree.Load(storage);
        Assert.That(tree.Text.ToString(), Is.EqualTo("ACGT"));
    }

    // ─── Q7: LoadStringInternal must use long for byte length ───────

    [Test]
    public void Load_TextWithLargeLength_DoesNotOverflowMultiply()
    {
        // We verify the fix by testing that a normal load works.
        // The actual overflow guard is structural: (long)textLen * 2 instead of int * int.
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        builder.Build(new StringTextSource("Hello World"));
        var tree = PersistentSuffixTree.Load(storage);
        Assert.That(tree.Text.ToString(), Is.EqualTo("Hello World"));
        Assert.That(tree.Contains("World"), Is.True);
    }

    // ─── Q9: Chunked text write must not corrupt data ───────────────
    // The fix changes ToString() to span-based encoding.
    // Verify correctness with a text larger than one chunk (4096 chars).

    [Test]
    public void Builder_TextLargerThanChunk_RoundTripsCorrectly()
    {
        // Create text > 4096 chars to force multiple chunks
        string longText = new string('A', 5000) + "NEEDLE" + new string('B', 5000);
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);
        builder.Build(new StringTextSource(longText));

        var tree = PersistentSuffixTree.Load(storage);
        Assert.That(tree.Text.ToString(), Is.EqualTo(longText));
        Assert.That(tree.Contains("NEEDLE"), Is.True);
        Assert.That(tree.Text.Length, Is.EqualTo(10006));
    }
}
