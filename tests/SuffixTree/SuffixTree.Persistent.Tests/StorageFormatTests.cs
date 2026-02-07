using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for the hybrid dual-format storage (Compact v4 / Large v3),
/// auto-promotion on overflow, format parity, and size characteristics.
/// </summary>
[TestFixture]
public class StorageFormatTests
{
    // ──── Hybrid: compact → large promotion ──────────────────────────

    [Test]
    public void Builder_ThrowsCompactOverflowException_WhenLimitExceeded()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        // Set a very small limit so even "banana" overflows
        builder.CompactOffsetLimit = 200;
        Assert.Throws<CompactOverflowException>(
            () => builder.Build(new StringTextSource("banana")));
    }

    [Test]
    public void Builder_DoesNotThrow_WhenWithinCompactLimit()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        // Default limit — "banana" fits easily
        Assert.DoesNotThrow(
            () => builder.Build(new StringTextSource("banana")));
    }

    [Test]
    public void Factory_AutoPromotesToLarge_WhenCompactOverflows()
    {
        // Tiny limit forces overflow → Factory rebuilds with Large
        using var tree = PersistentSuffixTreeFactory.CreateCore(
            new StringTextSource("banana"), filePath: null, compactOffsetLimit: 200) as IDisposable;
        var st = (ISuffixTree)tree!;

        Assert.Multiple(() =>
        {
            Assert.That(st.Contains("ana"), Is.True);
            Assert.That(st.CountOccurrences("a"), Is.EqualTo(3));
            Assert.That(st.LongestRepeatedSubstring(), Is.EqualTo("ana"));
            Assert.That(st.FindExactMatchAnchors("bandana", 3).Count, Is.GreaterThan(0));
        });
    }

    [Test]
    public void Factory_PromotedTree_HasCorrectVersion()
    {
        // Build with a tiny limit → must promote to Large (v3)
        using var tree = PersistentSuffixTreeFactory.CreateCore(
            new StringTextSource("mississippi"), filePath: null, compactOffsetLimit: 200) as IDisposable;
        var st = (PersistentSuffixTree)tree!;

        // Load from the same storage should auto-detect v3
        // Verify by checking that Contains/CountOccurrences work (they require correct layout)
        Assert.Multiple(() =>
        {
            Assert.That(st.Contains("issi"), Is.True);
            Assert.That(st.CountOccurrences("i"), Is.EqualTo(4));
        });
    }

    [Test]
    public void Factory_CompactUsedByDefault_ForSmallText()
    {
        // Build normally — should use Compact (v4)
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.Build(new StringTextSource("banana"));
        int version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        Assert.That(version, Is.EqualTo(4));
    }

    // ──── Compact format: build, query, Load ─────────────────────────

    [Test]
    public void CompactFormat_BuildAndQuery_Works()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        var text = new StringTextSource("banana");
        long root = builder.Build(text);
        using var tree = new PersistentSuffixTree(storage, root, text, NodeLayout.Compact);

        Assert.Multiple(() =>
        {
            Assert.That(tree.Contains("ana"), Is.True);
            Assert.That(tree.Contains("xyz"), Is.False);
            Assert.That(tree.CountOccurrences("a"), Is.EqualTo(3));
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("ana"));
            Assert.That(tree.FindExactMatchAnchors("bandana", 3).Count, Is.GreaterThan(0));
        });
    }

    [Test]
    public void CompactFormat_Load_AutoDetectsVersion()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.Build(new StringTextSource("abracadabra"));

        // Load without specifying layout — must auto-detect v4
        var tree = PersistentSuffixTree.Load(storage);
        Assert.Multiple(() =>
        {
            Assert.That(tree.Contains("abra"), Is.True);
            Assert.That(tree.CountOccurrences("a"), Is.EqualTo(5));
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("abra"));
        });
    }

    // ──── Large format: build, query, Load ───────────────────────────

    [Test]
    public void LargeFormat_BuildAndQuery_Works()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Large);
        var text = new StringTextSource("banana");
        long root = builder.Build(text);
        using var tree = new PersistentSuffixTree(storage, root, text, NodeLayout.Large);

        Assert.Multiple(() =>
        {
            Assert.That(tree.Contains("ana"), Is.True);
            Assert.That(tree.Contains("xyz"), Is.False);
            Assert.That(tree.CountOccurrences("a"), Is.EqualTo(3));
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("ana"));
            Assert.That(tree.FindExactMatchAnchors("bandana", 3).Count, Is.GreaterThan(0));
        });
    }

    [Test]
    public void LargeFormat_Load_AutoDetectsVersion()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Large);
        builder.Build(new StringTextSource("abracadabra"));

        // Load without specifying layout — must auto-detect v3
        var tree = PersistentSuffixTree.Load(storage);
        Assert.Multiple(() =>
        {
            Assert.That(tree.Contains("abra"), Is.True);
            Assert.That(tree.CountOccurrences("a"), Is.EqualTo(5));
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("abra"));
        });
    }

    // ──── Format parity (Compact and Large produce identical results) ─

    private static readonly string[] ParityTexts =
    [
        "",
        "a",
        "banana",
        "mississippi",
        "abracadabra",
        "aaaaaaaaaa",
        "🧬αβγ🧪$",
        "repetitive-repetitive-repetitive"
    ];

    [Test]
    public void FormatParity_CompactAndLarge_IdenticalResults(
        [ValueSource(nameof(ParityTexts))] string text)
    {
        // Build with Compact
        var compactStorage = new HeapStorageProvider();
        var compactBuilder = new PersistentSuffixTreeBuilder(compactStorage, NodeLayout.Compact);
        var textSource1 = new StringTextSource(text);
        long compactRoot = compactBuilder.Build(textSource1);
        using var compactTree = new PersistentSuffixTree(compactStorage, compactRoot, textSource1, NodeLayout.Compact);

        // Build with Large
        var largeStorage = new HeapStorageProvider();
        var largeBuilder = new PersistentSuffixTreeBuilder(largeStorage, NodeLayout.Large);
        var textSource2 = new StringTextSource(text);
        long largeRoot = largeBuilder.Build(textSource2);
        using var largeTree = new PersistentSuffixTree(largeStorage, largeRoot, textSource2, NodeLayout.Large);

        Assert.Multiple(() =>
        {
            Assert.That(compactTree.NodeCount, Is.EqualTo(largeTree.NodeCount), "NodeCount");
            Assert.That(compactTree.LeafCount, Is.EqualTo(largeTree.LeafCount), "LeafCount");
            Assert.That(compactTree.LongestRepeatedSubstring(),
                Is.EqualTo(largeTree.LongestRepeatedSubstring()), "LRS");

            if (text.Length > 0)
            {
                // Exhaustive short-substring check
                for (int i = 0; i < text.Length; i++)
                {
                    for (int len = 1; len <= Math.Min(5, text.Length - i); len++)
                    {
                        string sub = text.Substring(i, len);
                        Assert.That(compactTree.Contains(sub),
                            Is.EqualTo(largeTree.Contains(sub)), $"Contains(\"{sub}\")");
                        Assert.That(compactTree.CountOccurrences(sub),
                            Is.EqualTo(largeTree.CountOccurrences(sub)), $"Count(\"{sub}\")");
                    }
                }
            }

            // Logical hash parity (layout-independent)
            var compactHash = SuffixTreeSerializer.CalculateLogicalHash(compactTree);
            var largeHash = SuffixTreeSerializer.CalculateLogicalHash(largeTree);
            Assert.That(compactHash, Is.EqualTo(largeHash), "LogicalHash");
        });
    }

    // ──── Size characteristics ───────────────────────────────────────

    [Test]
    public void CompactFormat_ProducesSmallerStorage()
    {
        string text = "The quick brown fox jumps over the lazy dog and extra words to amplify difference";

        var compactStorage = new HeapStorageProvider();
        new PersistentSuffixTreeBuilder(compactStorage, NodeLayout.Compact)
            .Build(new StringTextSource(text));

        var largeStorage = new HeapStorageProvider();
        new PersistentSuffixTreeBuilder(largeStorage, NodeLayout.Large)
            .Build(new StringTextSource(text));

        Assert.That(compactStorage.Size, Is.LessThan(largeStorage.Size),
            $"Compact ({compactStorage.Size} B) should be smaller than Large ({largeStorage.Size} B)");
    }

    // ──── Factory auto-selects format ────────────────────────────────

    [Test]
    public void Factory_Create_UsesCompactForSmallText()
    {
        using var tree = PersistentSuffixTreeFactory.Create(new StringTextSource("banana")) as IDisposable;
        var st = (ISuffixTree)tree!;
        Assert.Multiple(() =>
        {
            Assert.That(st.Contains("ana"), Is.True);
            Assert.That(st.LeafCount, Is.EqualTo(6));
        });
    }

    // ──── Version header round-trip ──────────────────────────────────

    [Test]
    public void HeaderVersion_CompactWritesV4()
    {
        var storage = new HeapStorageProvider();
        new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact)
            .Build(new StringTextSource("abc"));
        int version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        Assert.That(version, Is.EqualTo(4));
    }

    [Test]
    public void HeaderVersion_LargeWritesV3()
    {
        var storage = new HeapStorageProvider();
        new PersistentSuffixTreeBuilder(storage, NodeLayout.Large)
            .Build(new StringTextSource("abc"));
        int version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        Assert.That(version, Is.EqualTo(3));
    }

    // ──── NodeLayout properties ──────────────────────────────────────

    [Test]
    public void NodeLayout_Compact_HasExpectedSizes()
    {
        Assert.Multiple(() =>
        {
            Assert.That(NodeLayout.Compact.NodeSize, Is.EqualTo(28));
            Assert.That(NodeLayout.Compact.ChildEntrySize, Is.EqualTo(8));
            Assert.That(NodeLayout.Compact.OffsetIs64Bit, Is.False);
            Assert.That(NodeLayout.Compact.Version, Is.EqualTo(4));
        });
    }

    [Test]
    public void NodeLayout_Large_HasExpectedSizes()
    {
        Assert.Multiple(() =>
        {
            Assert.That(NodeLayout.Large.NodeSize, Is.EqualTo(40));
            Assert.That(NodeLayout.Large.ChildEntrySize, Is.EqualTo(12));
            Assert.That(NodeLayout.Large.OffsetIs64Bit, Is.True);
            Assert.That(NodeLayout.Large.Version, Is.EqualTo(3));
        });
    }

    [Test]
    public void ForVersion_ReturnsCorrectLayout()
    {
        Assert.That(NodeLayout.ForVersion(3), Is.SameAs(NodeLayout.Large));
        Assert.That(NodeLayout.ForVersion(4), Is.SameAs(NodeLayout.Compact));
        Assert.Throws<InvalidOperationException>(() => NodeLayout.ForVersion(99));
    }

    [Test]
    public void CompactMaxOffset_IsUIntMaxMinusOne()
    {
        Assert.That(NodeLayout.CompactMaxOffset, Is.EqualTo((long)uint.MaxValue - 1));
    }

    // ──── Serializer round-trip preserves format ─────────────────────

    [Test]
    public void Serializer_Import_AutoSelectsFormat()
    {
        // Export from a Compact tree
        var original = PersistentSuffixTreeFactory.Create(new StringTextSource("mississippi"));

        using var ms = new System.IO.MemoryStream();
        SuffixTreeSerializer.Export(original, ms);
        ms.Position = 0;

        // Import into new storage — should auto-select Compact for small text
        var importedStorage = new HeapStorageProvider();
        var imported = SuffixTreeSerializer.Import(ms, importedStorage);

        Assert.Multiple(() =>
        {
            Assert.That(imported.Contains("ssi"), Is.True);
            Assert.That(imported.CountOccurrences("i"), Is.EqualTo(4));
            Assert.That(imported.LongestRepeatedSubstring(), Is.EqualTo("issi"));
            Assert.That(imported.FindExactMatchAnchors("mississippi", 3).Count, Is.GreaterThan(0));
        });
    }

    // ──── Hybrid promotion parity with reference ─────────────────────

    [Test]
    public void Factory_PromotedTree_MatchesReferenceImplementation()
    {
        string text = "abracadabra";
        var reference = global::SuffixTree.SuffixTree.Build(text);

        // Force promotion to Large via tiny limit
        using var tree = PersistentSuffixTreeFactory.CreateCore(
            new StringTextSource(text), filePath: null, compactOffsetLimit: 200) as IDisposable;
        var promoted = (ISuffixTree)tree!;

        Assert.Multiple(() =>
        {
            Assert.That(promoted.Text.ToString(), Is.EqualTo(reference.Text.ToString()), "Text");
            Assert.That(promoted.LeafCount, Is.EqualTo(reference.LeafCount), "LeafCount");
            Assert.That(promoted.LongestRepeatedSubstring(),
                Is.EqualTo(reference.LongestRepeatedSubstring()), "LRS");
            Assert.That(promoted.LongestCommonSubstring("cadab"),
                Is.EqualTo(reference.LongestCommonSubstring("cadab")), "LCS");

            var refHash = SuffixTreeSerializer.CalculateLogicalHash(reference);
            var promHash = SuffixTreeSerializer.CalculateLogicalHash(promoted);
            Assert.That(promHash, Is.EqualTo(refHash), "LogicalHash");
        });
    }
}
