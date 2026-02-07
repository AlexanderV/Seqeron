using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests for the hybrid dual-format storage (Compact v4 / Large v3 / Hybrid v5),
/// hybrid continuation with jump table, format parity, and size characteristics.
/// </summary>
[TestFixture]
public class StorageFormatTests
{
    // ──── Hybrid continuation: compact → large mid-build ─────────────

    [Test]
    public void Builder_TransitionsToLargeZone_WhenLimitExceeded()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        // Set a very small limit so even "banana" triggers transition
        builder.CompactOffsetLimit = 200;
        var text = new StringTextSource("banana");
        long root = builder.Build(text);

        // Builder should have transitioned (hybrid)
        Assert.Multiple(() =>
        {
            Assert.That(builder.IsHybrid, Is.True, "Builder should be hybrid after transition");
            Assert.That(builder.TransitionOffset, Is.GreaterThan(0), "TransitionOffset > 0");
            Assert.That(builder.JumpTableEnd, Is.GreaterThanOrEqualTo(builder.TransitionOffset),
                "JumpTableEnd >= TransitionOffset");
        });
    }

    [Test]
    public void Builder_DoesNotTransition_WhenWithinCompactLimit()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        // Default limit — "banana" fits easily
        var text = new StringTextSource("banana");
        long root = builder.Build(text);
        Assert.That(builder.IsHybrid, Is.False, "Should not transition for small text");
    }

    [Test]
    public void Factory_HybridContinuation_ProducesCorrectTree()
    {
        // Tiny limit forces mid-build transition → Factory produces hybrid tree
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
    public void Factory_HybridTree_HasVersion5Header()
    {
        // Build with a tiny limit → must produce Hybrid v5
        using var tree = PersistentSuffixTreeFactory.CreateCore(
            new StringTextSource("mississippi"), filePath: null, compactOffsetLimit: 200) as IDisposable;
        var pst = (PersistentSuffixTree)tree!;

        Assert.Multiple(() =>
        {
            Assert.That(pst.IsHybrid, Is.True, "Tree should be hybrid");
            Assert.That(pst.TransitionOffset, Is.GreaterThan(0), "TransitionOffset > 0");
            Assert.That(pst.JumpTableEnd, Is.GreaterThanOrEqualTo(pst.TransitionOffset));
            Assert.That(pst.Contains("issi"), Is.True);
            Assert.That(pst.CountOccurrences("i"), Is.EqualTo(4));
        });
    }

    [Test]
    public void Factory_CompactUsedByDefault_ForSmallText()
    {
        // Build normally — should use Compact (v4), no transition
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.Build(new StringTextSource("banana"));
        int version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        Assert.That(version, Is.EqualTo(4));
    }

    [Test]
    public void Builder_HybridHeader_WritesVersion5()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = 200;
        builder.Build(new StringTextSource("banana"));
        int version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        Assert.That(version, Is.EqualTo(5), "Hybrid builds should write version 5");
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
        Assert.That(NodeLayout.ForVersion(5), Is.SameAs(NodeLayout.Compact));
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

    // ──── Hybrid continuation parity with reference ────────────────────

    [Test]
    public void Factory_HybridTree_MatchesReferenceImplementation()
    {
        string text = "abracadabra";
        var reference = global::SuffixTree.SuffixTree.Build(text);

        // Force hybrid continuation via tiny limit
        using var tree = PersistentSuffixTreeFactory.CreateCore(
            new StringTextSource(text), filePath: null, compactOffsetLimit: 200) as IDisposable;
        var hybrid = (ISuffixTree)tree!;

        Assert.Multiple(() =>
        {
            Assert.That(hybrid.Text.ToString(), Is.EqualTo(reference.Text.ToString()), "Text");
            Assert.That(hybrid.LeafCount, Is.EqualTo(reference.LeafCount), "LeafCount");
            Assert.That(hybrid.LongestRepeatedSubstring(),
                Is.EqualTo(reference.LongestRepeatedSubstring()), "LRS");
            Assert.That(hybrid.LongestCommonSubstring("cadab"),
                Is.EqualTo(reference.LongestCommonSubstring("cadab")), "LCS");

            var refHash = SuffixTreeSerializer.CalculateLogicalHash(reference);
            var hybridHash = SuffixTreeSerializer.CalculateLogicalHash(hybrid);
            Assert.That(hybridHash, Is.EqualTo(refHash), "LogicalHash");
        });
    }

    [Test]
    public void Factory_HybridTree_MatchesCompactTree()
    {
        // Build a pure Compact tree and a hybrid tree for the same text
        // They should produce identical query results
        string text = "mississippi";
        var textSource1 = new StringTextSource(text);
        var textSource2 = new StringTextSource(text);

        using var compactTree = PersistentSuffixTreeFactory.Create(textSource1) as IDisposable;
        using var hybridTree = PersistentSuffixTreeFactory.CreateCore(
            textSource2, filePath: null, compactOffsetLimit: 200) as IDisposable;

        var compact = (ISuffixTree)compactTree!;
        var hybrid = (ISuffixTree)hybridTree!;

        Assert.Multiple(() =>
        {
            Assert.That(hybrid.NodeCount, Is.EqualTo(compact.NodeCount), "NodeCount");
            Assert.That(hybrid.LeafCount, Is.EqualTo(compact.LeafCount), "LeafCount");
            Assert.That(hybrid.LongestRepeatedSubstring(),
                Is.EqualTo(compact.LongestRepeatedSubstring()), "LRS");
            Assert.That(hybrid.CountOccurrences("issi"),
                Is.EqualTo(compact.CountOccurrences("issi")), "CountOccurrences");

            var compactHash = SuffixTreeSerializer.CalculateLogicalHash(compact);
            var hybridHash = SuffixTreeSerializer.CalculateLogicalHash(hybrid);
            Assert.That(hybridHash, Is.EqualTo(compactHash), "LogicalHash");
        });
    }
}
