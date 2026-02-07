using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Thorough coverage of the Compact→Large transition zone in the hybrid builder.
/// <para>
/// Strategy: sweep <c>compactOffsetLimit</c> through every possible node allocation
/// boundary for several representative texts.  For each limit value that produces a
/// hybrid tree we verify:
/// <list type="bullet">
///   <item>Full ISuffixTree API parity with a pure-Compact and reference (in-memory) tree.</item>
///   <item>Structural invariants: NodeCount, LeafCount, LogicalHash equality.</item>
///   <item>Jump table consistency: JumpTableStart ≥ TransitionOffset, contiguous.</item>
///   <item>All suffixes enumerable and sorted.</item>
///   <item>Load() round-trip (auto-detect v5).</item>
/// </list>
/// </para>
/// </summary>
[TestFixture]
public class HybridTransitionZoneTests
{
    // ────────── representative texts with different characteristics ──────────

    /// <summary>
    /// Texts that exercise different suffix-tree topologies:
    /// mono-run, repeated pattern, biological-like, unique chars, unicode.
    /// </summary>
    private static readonly string[] TransitionTexts =
    [
        "banana",                                      // classic, short
        "mississippi",                                 // many repeated substrings
        "abracadabra",                                 // overlapping repeats
        "aaaaaaaaaa",                                  // mono-character run
        "abcdefghij",                                  // all unique chars
        "aababcabcdabcde",                             // nested repeats of increasing length
        "ATCGATCGATCG",                                // biological repeat
        "the_quick_brown_fox_jumps_over_the_lazy_dog", // long, many unique
    ];

    // ──────────── Sweep: every possible transition boundary ──────────────

    /// <summary>
    /// For each text, sweep compactOffsetLimit from the header size up to the point
    /// where all nodes fit in Compact.  Every distinct limit that produces a hybrid
    /// tree is checked for full parity with the pure-Compact and reference trees.
    /// </summary>
    [Test]
    public void TransitionSweep_EveryBoundary_MatchesCompactTree(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        // Build reference trees once
        var reference = global::SuffixTree.SuffixTree.Build(text);
        var compactTree = BuildCompact(text);

        // The builder starts at offset 72 (header) + 28 (root) = 100.
        // Each subsequent Compact node adds 28 bytes.
        // Sweep limits from just after the header to beyond where all nodes fit.
        long maxLimit = compactTree.storage.Size;
        int headerPlusRoot = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize;
        int step = NodeLayout.Compact.NodeSize; // 28 bytes — one node boundary

        int hybridCount = 0;
        long? firstHybridLimit = null;
        long? lastHybridLimit = null;

        for (long limit = headerPlusRoot; limit < maxLimit; limit += step)
        {
            var hybrid = BuildHybrid(text, limit);
            if (!hybrid.pst.IsHybrid) continue;

            hybridCount++;
            firstHybridLimit ??= limit;
            lastHybridLimit = limit;

            // ── Structural invariants ──
            Assert.That(hybrid.pst.IsHybrid, Is.True, $"limit={limit}: should be hybrid");
            Assert.That(hybrid.pst.TransitionOffset, Is.GreaterThan(0), $"limit={limit}: TransitionOffset > 0");
            Assert.That(hybrid.pst.JumpTableStart,
                Is.GreaterThanOrEqualTo(hybrid.pst.TransitionOffset).Or.EqualTo(-1),
                $"limit={limit}: JumpTableStart >= TransitionOffset or -1 (no cross-zone refs)");
            if (hybrid.pst.JumpTableStart >= 0)
            {
                Assert.That(hybrid.pst.JumpTableEnd, Is.GreaterThan(hybrid.pst.JumpTableStart),
                    $"limit={limit}: JumpTableEnd > JumpTableStart");
                // Jump table size must be a multiple of 8
                long jtSize = hybrid.pst.JumpTableEnd - hybrid.pst.JumpTableStart;
                Assert.That(jtSize % 8, Is.Zero, $"limit={limit}: Jump table size must be multiple of 8");
            }

            // ── API parity with Compact tree ──
            AssertFullParity(hybrid.tree, compactTree.tree, reference, text, $"limit={limit}");

            hybrid.Dispose();
        }

        // At least some limits should produce hybrid trees
        TestContext.Out.WriteLine(
            $"Text '{text}': {hybridCount} hybrid transitions, first@{firstHybridLimit}, last@{lastHybridLimit}");
        Assert.That(hybridCount, Is.GreaterThan(0), $"Expected at least one hybrid transition for '{text}'");

        compactTree.Dispose();
    }

    // ──────────── Specific interesting transition points ──────────────

    /// <summary>
    /// Transition at the very first possible node after root — absolute minimum compact zone.
    /// </summary>
    [Test]
    public void Transition_AtFirstNodeAfterRoot_ProducesCorrectTree(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        // limit = header + root node → next allocation triggers transition
        long limit = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize;
        var hybrid = BuildHybrid(text, limit);
        var compact = BuildCompact(text);
        var reference = global::SuffixTree.SuffixTree.Build(text);

        Assert.That(hybrid.pst.IsHybrid, Is.True, "Should be hybrid");
        AssertFullParity(hybrid.tree, compact.tree, reference, text, "first-node transition");

        hybrid.Dispose();
        compact.Dispose();
    }

    /// <summary>
    /// Transition at a boundary where only root + 1 compact node fit.
    /// </summary>
    [Test]
    public void Transition_RootPlusOneNode_ProducesCorrectTree(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        long limit = PersistentConstants.HEADER_SIZE_V5 + 2 * NodeLayout.Compact.NodeSize;
        var hybrid = BuildHybrid(text, limit);
        var compact = BuildCompact(text);
        var reference = global::SuffixTree.SuffixTree.Build(text);

        Assert.That(hybrid.pst.IsHybrid, Is.True, "Should be hybrid");
        AssertFullParity(hybrid.tree, compact.tree, reference, text, "root+1 transition");

        hybrid.Dispose();
        compact.Dispose();
    }

    /// <summary>
    /// Transition at roughly the midpoint of the tree (half the nodes in compact zone).
    /// </summary>
    [Test]
    public void Transition_AtMidpoint_ProducesCorrectTree(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var compact = BuildCompact(text);
        long midpoint = PersistentConstants.HEADER_SIZE_V5
            + (long)(compact.tree.NodeCount / 2) * NodeLayout.Compact.NodeSize;
        compact.Dispose();

        var hybrid = BuildHybrid(text, midpoint);
        var compact2 = BuildCompact(text);
        var reference = global::SuffixTree.SuffixTree.Build(text);

        if (hybrid.pst.IsHybrid)
            AssertFullParity(hybrid.tree, compact2.tree, reference, text, "midpoint transition");

        hybrid.Dispose();
        compact2.Dispose();
    }

    /// <summary>
    /// Transition at the very last node — only 1 node ends up in the large zone.
    /// </summary>
    [Test]
    public void Transition_AtLastNode_ProducesCorrectTree(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var compact = BuildCompact(text);
        // Last node boundary: total nodes - 1 → everything except the last node fits in compact
        long limit = PersistentConstants.HEADER_SIZE_V5
            + (long)(compact.tree.NodeCount - 1) * NodeLayout.Compact.NodeSize;
        compact.Dispose();

        var hybrid = BuildHybrid(text, limit);
        var compact2 = BuildCompact(text);
        var reference = global::SuffixTree.SuffixTree.Build(text);

        if (hybrid.pst.IsHybrid)
            AssertFullParity(hybrid.tree, compact2.tree, reference, text, "last-node transition");

        hybrid.Dispose();
        compact2.Dispose();
    }

    // ──────────── v5 header round-trip via Load() ──────────────

    [Test]
    public void HybridTree_LoadAutoDetectsV5_AndQueriesWork(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = PersistentConstants.HEADER_SIZE_V5 + 2 * NodeLayout.Compact.NodeSize;
        builder.Build(new StringTextSource(text));

        if (!builder.IsHybrid)
        {
            Assert.Pass("Text too short to trigger transition at this limit");
            return;
        }

        // Load from storage — must auto-detect v5
        var loaded = PersistentSuffixTree.Load(storage);

        Assert.Multiple(() =>
        {
            Assert.That(loaded.IsHybrid, Is.True, "Loaded tree should be hybrid");
            Assert.That(loaded.TransitionOffset, Is.EqualTo(builder.TransitionOffset), "TransitionOffset persisted");
            Assert.That(loaded.JumpTableStart, Is.EqualTo(builder.JumpTableStart), "JumpTableStart persisted");
            Assert.That(loaded.JumpTableEnd, Is.EqualTo(builder.JumpTableEnd), "JumpTableEnd persisted");

            var reference = global::SuffixTree.SuffixTree.Build(text);
            Assert.That(loaded.NodeCount, Is.EqualTo(reference.NodeCount), "NodeCount after Load");
            Assert.That(loaded.LeafCount, Is.EqualTo(reference.LeafCount), "LeafCount after Load");
            Assert.That(loaded.Contains(text), Is.True, "Contains full text after Load");
            Assert.That(loaded.LongestRepeatedSubstring(),
                Is.EqualTo(reference.LongestRepeatedSubstring()), "LRS after Load");
        });
    }

    // ──────────── Suffix enumeration matches reference ──────────────

    [Test]
    public void HybridTree_EnumerateSuffixes_MatchesReference(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var reference = global::SuffixTree.SuffixTree.Build(text);
        var hybrid = BuildHybrid(text, PersistentConstants.HEADER_SIZE_V5 + 3 * NodeLayout.Compact.NodeSize);
        if (!hybrid.pst.IsHybrid)
        {
            hybrid.Dispose();
            Assert.Pass("Text too short to trigger hybrid");
            return;
        }

        var refSuffixes = reference.GetAllSuffixes().OrderBy(s => s, StringComparer.Ordinal).ToList();
        var hybridSuffixes = hybrid.tree.GetAllSuffixes().OrderBy(s => s, StringComparer.Ordinal).ToList();

        Assert.That(hybridSuffixes, Is.EqualTo(refSuffixes), "All suffixes must match reference");
        hybrid.Dispose();
    }

    // ──────────── FindAllOccurrences parity ──────────────

    [Test]
    public void HybridTree_FindAllOccurrences_MatchesCompact(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var compact = BuildCompact(text);
        var hybrid = BuildHybrid(text, 200);
        if (!hybrid.pst.IsHybrid) { hybrid.Dispose(); compact.Dispose(); Assert.Pass(); return; }

        // Check every 2-char substring
        for (int i = 0; i <= text.Length - 2; i++)
        {
            string pattern = text.Substring(i, 2);
            var compactOcc = compact.tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
            var hybridOcc = hybrid.tree.FindAllOccurrences(pattern).OrderBy(x => x).ToList();
            Assert.That(hybridOcc, Is.EqualTo(compactOcc), $"FindAllOccurrences(\"{pattern}\")");
        }

        hybrid.Dispose();
        compact.Dispose();
    }

    // ──────────── Traverse produces identical structure ──────────────

    [Test]
    public void HybridTree_Traverse_MatchesCompact(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var compact = BuildCompact(text);
        var hybrid = BuildHybrid(text, 200);
        if (!hybrid.pst.IsHybrid) { hybrid.Dispose(); compact.Dispose(); Assert.Pass(); return; }

        var compactVisits = new List<(int Start, int End, int Leaves, int Children, int Depth)>();
        compact.tree.Traverse(new CollectingVisitor(compactVisits));

        var hybridVisits = new List<(int Start, int End, int Leaves, int Children, int Depth)>();
        hybrid.tree.Traverse(new CollectingVisitor(hybridVisits));

        Assert.That(hybridVisits.Count, Is.EqualTo(compactVisits.Count), "Visit count");
        for (int i = 0; i < compactVisits.Count; i++)
        {
            Assert.That(hybridVisits[i], Is.EqualTo(compactVisits[i]),
                $"Visit[{i}]: hybrid={hybridVisits[i]} compact={compactVisits[i]}");
        }

        hybrid.Dispose();
        compact.Dispose();
    }

    // ──────────── Export/Import round-trip for hybrid trees ──────────────

    [Test]
    public void HybridTree_ExportImport_RoundTrip(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var hybrid = BuildHybrid(text, 200);
        if (!hybrid.pst.IsHybrid) { hybrid.Dispose(); Assert.Pass(); return; }

        // Export
        using var ms = new System.IO.MemoryStream();
        SuffixTreeSerializer.Export(hybrid.tree, ms);
        ms.Position = 0;

        // Import
        var importedStorage = new HeapStorageProvider();
        var imported = SuffixTreeSerializer.Import(ms, importedStorage);

        Assert.Multiple(() =>
        {
            Assert.That(imported.NodeCount, Is.EqualTo(hybrid.tree.NodeCount), "NodeCount");
            Assert.That(imported.LeafCount, Is.EqualTo(hybrid.tree.LeafCount), "LeafCount");
            Assert.That(imported.LongestRepeatedSubstring(),
                Is.EqualTo(hybrid.tree.LongestRepeatedSubstring()), "LRS");
            Assert.That(imported.Contains(text), Is.True, "Contains full text");

            var origHash = SuffixTreeSerializer.CalculateLogicalHash(hybrid.tree);
            var importedHash = SuffixTreeSerializer.CalculateLogicalHash(imported);
            Assert.That(importedHash, Is.EqualTo(origHash), "LogicalHash");
        });

        hybrid.Dispose();
    }

    // ──────────── LCS (suffix-link-dependent) works across zones ──────────────

    [Test]
    public void HybridTree_LongestCommonSubstring_MatchesReference(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var reference = global::SuffixTree.SuffixTree.Build(text);
        var hybrid = BuildHybrid(text, 200);
        if (!hybrid.pst.IsHybrid) { hybrid.Dispose(); Assert.Pass(); return; }

        // LCS with itself should be the full text
        Assert.That(hybrid.tree.LongestCommonSubstring(text),
            Is.EqualTo(reference.LongestCommonSubstring(text)), "LCS with self");

        // LCS with a rotated version (exercises suffix link traversal)
        if (text.Length >= 4)
        {
            string rotated = text.Substring(text.Length / 2) + text.Substring(0, text.Length / 2);
            Assert.That(hybrid.tree.LongestCommonSubstring(rotated),
                Is.EqualTo(reference.LongestCommonSubstring(rotated)), "LCS with rotated");
        }

        // LCS with a completely unrelated string
        Assert.That(hybrid.tree.LongestCommonSubstring("zzz"),
            Is.EqualTo(reference.LongestCommonSubstring("zzz")), "LCS with unrelated");

        hybrid.Dispose();
    }

    // ──────────── FindExactMatchAnchors across zones ──────────────

    [Test]
    public void HybridTree_FindExactMatchAnchors_MatchesCompact(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        if (text.Length < 6) { Assert.Pass("Text too short for meaningful anchors test"); return; }

        var compact = BuildCompact(text);
        var hybrid = BuildHybrid(text, 200);
        if (!hybrid.pst.IsHybrid) { hybrid.Dispose(); compact.Dispose(); Assert.Pass(); return; }

        // Use a query that overlaps with the text
        string query = text.Substring(1) + text.Substring(0, 3);

        var compactAnchors = compact.tree.FindExactMatchAnchors(query, 2)
            .OrderBy(a => a.PositionInText).ThenBy(a => a.PositionInQuery).ToList();
        var hybridAnchors = hybrid.tree.FindExactMatchAnchors(query, 2)
            .OrderBy(a => a.PositionInText).ThenBy(a => a.PositionInQuery).ToList();

        Assert.That(hybridAnchors, Is.EqualTo(compactAnchors), "FindExactMatchAnchors");

        hybrid.Dispose();
        compact.Dispose();
    }

    // ══════════════════════════════════════════════════════════════════
    //  Compact upper-bound tests: trees near the uint32 boundary
    // ══════════════════════════════════════════════════════════════════

    /// <summary>
    /// Verify that a tree built with compactOffsetLimit just above the actual tree size
    /// stays pure Compact (no unnecessary transition).
    /// </summary>
    [Test]
    public void UpperBound_LimitJustAboveTreeSize_StaysCompact(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var probe = BuildCompact(text);
        long treeSize = probe.storage.Size;
        probe.Dispose();

        // Build with limit = treeSize + 1 → should NOT transition
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = treeSize + 1;
        builder.Build(new StringTextSource(text));

        Assert.That(builder.IsHybrid, Is.False,
            $"Should stay pure Compact when limit ({treeSize + 1}) > tree size ({treeSize})");

        int version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        Assert.That(version, Is.EqualTo(4), "Version should be 4 (Compact)");
    }

    /// <summary>
    /// Verify that a tree built with compactOffsetLimit == exact tree size triggers hybrid.
    /// (Because the last allocation pushes storage.Size == limit, which still fits,
    /// but subsequent child arrays push past.)
    /// </summary>
    [Test]
    public void UpperBound_LimitAtExactTreeSize_OrJustBelow_IsHybrid(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var probe = BuildCompact(text);
        long treeSize = probe.storage.Size;
        int nodeCount = probe.tree.NodeCount;
        probe.Dispose();

        if (nodeCount <= 2) { Assert.Pass("Too few nodes to test boundary"); return; }

        // Set limit to be just under the last node allocation
        // This should force the last node(s) into large zone
        long limit = treeSize - NodeLayout.Compact.NodeSize;
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = limit;
        long root = builder.Build(new StringTextSource(text));

        if (builder.IsHybrid)
        {
            var tree = new PersistentSuffixTree(storage, root, new StringTextSource(text),
                NodeLayout.Compact, builder.TransitionOffset, builder.JumpTableStart, builder.JumpTableEnd);
            var reference = global::SuffixTree.SuffixTree.Build(text);

            Assert.Multiple(() =>
            {
                Assert.That(tree.Contains(text), Is.True, "Contains full text");
                Assert.That(tree.NodeCount, Is.EqualTo(reference.NodeCount), "NodeCount");
                Assert.That(tree.LeafCount, Is.EqualTo(reference.LeafCount), "LeafCount");
                Assert.That(tree.LongestRepeatedSubstring(),
                    Is.EqualTo(reference.LongestRepeatedSubstring()), "LRS");
            });
        }
    }

    /// <summary>
    /// Verify that the default CompactMaxOffset (~4 GB) does not cause transition for
    /// any of our test texts (they're all small).
    /// </summary>
    [Test]
    public void UpperBound_DefaultLimit_NeverTriggersForSmallTexts(
        [ValueSource(nameof(TransitionTexts))] string text)
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        // No override → default limit = CompactMaxOffset ≈ 4 GB
        builder.Build(new StringTextSource(text));

        Assert.That(builder.IsHybrid, Is.False, "Default limit should not trigger hybrid for small texts");
    }

    /// <summary>
    /// Stress-test: build with a larger text to verify more nodes and deeper transitions.
    /// Uses a 1000-char repetitive string that creates many internal nodes.
    /// </summary>
    [Test]
    public void UpperBound_LargerText_HybridMatchesCompact()
    {
        // Generate a repetitive text with interesting structure
        string baseText = "abcabc_defdef_ghighi_jkljkl_mnomno";
        var sb = new System.Text.StringBuilder(1200);
        while (sb.Length < 1000)
            sb.Append(baseText);
        string text = sb.ToString(0, 1000);

        var compact = BuildCompact(text);
        var reference = global::SuffixTree.SuffixTree.Build(text);

        // Try several transition points
        long maxLimit = compact.storage.Size;
        long[] limits = [
            PersistentConstants.HEADER_SIZE_V5 + 3 * NodeLayout.Compact.NodeSize,   // very early
            maxLimit / 4,                                                              // 25%
            maxLimit / 2,                                                              // 50%
            maxLimit * 3 / 4,                                                          // 75%
            maxLimit - NodeLayout.Compact.NodeSize,                                    // very late
        ];

        foreach (long limit in limits)
        {
            var hybrid = BuildHybrid(text, limit);
            if (!hybrid.pst.IsHybrid) { hybrid.Dispose(); continue; }

            Assert.Multiple(() =>
            {
                Assert.That(hybrid.tree.NodeCount, Is.EqualTo(compact.tree.NodeCount),
                    $"limit={limit}: NodeCount");
                Assert.That(hybrid.tree.LeafCount, Is.EqualTo(compact.tree.LeafCount),
                    $"limit={limit}: LeafCount");
                Assert.That(hybrid.tree.LongestRepeatedSubstring(),
                    Is.EqualTo(compact.tree.LongestRepeatedSubstring()),
                    $"limit={limit}: LRS");

                var compactHash = SuffixTreeSerializer.CalculateLogicalHash(compact.tree);
                var hybridHash = SuffixTreeSerializer.CalculateLogicalHash(hybrid.tree);
                Assert.That(hybridHash, Is.EqualTo(compactHash),
                    $"limit={limit}: LogicalHash");
            });

            hybrid.Dispose();
        }

        compact.Dispose();
    }

    /// <summary>
    /// Verify no off-by-one: transition at limit = headerSize + rootSize + 1 byte.
    /// The +1 makes the root barely fit, but the next node won't.
    /// </summary>
    [Test]
    public void OffByOne_LimitJustPastRoot_TriggersAtSecondNode()
    {
        // Root at 72, size 28 → ends at 100.  Next node starts at 100.
        // If limit = 101, node at 100 with size 28 → 100+28=128 > 101 → transition.
        long limit = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize + 1;
        string text = "banana";

        var hybrid = BuildHybrid(text, limit);
        var compact = BuildCompact(text);
        var reference = global::SuffixTree.SuffixTree.Build(text);

        Assert.That(hybrid.pst.IsHybrid, Is.True, "Should trigger at second node");
        Assert.That(hybrid.pst.TransitionOffset,
            Is.EqualTo(PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize),
            "Transition should be right after root");

        AssertFullParity(hybrid.tree, compact.tree, reference, text, "off-by-one");

        hybrid.Dispose();
        compact.Dispose();
    }

    // ══════════════════════════════════════════════════════════════════
    //  Helpers
    // ══════════════════════════════════════════════════════════════════

    /// <summary>
    /// Full ISuffixTree parity assertion between hybrid, compact and reference trees.
    /// </summary>
    private static void AssertFullParity(
        ISuffixTree hybrid, ISuffixTree compact, ISuffixTree reference,
        string text, string context)
    {
        Assert.Multiple(() =>
        {
            // Structural
            Assert.That(hybrid.NodeCount, Is.EqualTo(compact.NodeCount),
                $"{context}: NodeCount (hybrid vs compact)");
            Assert.That(hybrid.NodeCount, Is.EqualTo(reference.NodeCount),
                $"{context}: NodeCount (hybrid vs reference)");
            Assert.That(hybrid.LeafCount, Is.EqualTo(compact.LeafCount),
                $"{context}: LeafCount");
            Assert.That(hybrid.MaxDepth, Is.EqualTo(compact.MaxDepth),
                $"{context}: MaxDepth");

            // Contains — full text and every single-char
            Assert.That(hybrid.Contains(text), Is.True, $"{context}: Contains(full text)");
            foreach (char c in text.Distinct())
            {
                Assert.That(hybrid.Contains(c.ToString()), Is.EqualTo(compact.Contains(c.ToString())),
                    $"{context}: Contains(\"{c}\")");
            }

            // Contains — every substring up to length 5
            for (int i = 0; i < text.Length; i++)
            {
                for (int len = 1; len <= Math.Min(5, text.Length - i); len++)
                {
                    string sub = text.Substring(i, len);
                    Assert.That(hybrid.Contains(sub), Is.EqualTo(compact.Contains(sub)),
                        $"{context}: Contains(\"{sub}\")");
                    Assert.That(hybrid.CountOccurrences(sub), Is.EqualTo(compact.CountOccurrences(sub)),
                        $"{context}: Count(\"{sub}\")");
                }
            }

            // Absent substrings
            Assert.That(hybrid.Contains("ZZZZZZ"), Is.False, $"{context}: Contains absent");
            Assert.That(hybrid.CountOccurrences("ZZZZZZ"), Is.EqualTo(0), $"{context}: Count absent");

            // LRS
            Assert.That(hybrid.LongestRepeatedSubstring(),
                Is.EqualTo(compact.LongestRepeatedSubstring()),
                $"{context}: LRS");
            Assert.That(hybrid.LongestRepeatedSubstring(),
                Is.EqualTo(reference.LongestRepeatedSubstring()),
                $"{context}: LRS vs reference");

            // LCS with a probe string
            if (text.Length >= 3)
            {
                string probe = text.Substring(text.Length / 3, Math.Min(5, text.Length / 3));
                Assert.That(hybrid.LongestCommonSubstring(probe),
                    Is.EqualTo(compact.LongestCommonSubstring(probe)),
                    $"{context}: LCS(\"{probe}\")");
            }

            // LogicalHash — the definitive structural equality check
            var compactHash = SuffixTreeSerializer.CalculateLogicalHash(compact);
            var hybridHash = SuffixTreeSerializer.CalculateLogicalHash(hybrid);
            Assert.That(hybridHash, Is.EqualTo(compactHash),
                $"{context}: LogicalHash");
        });
    }

    // ──────────── Builder wrappers with lifecycle tracking ──────────────

    private record struct TreeHolder(ISuffixTree tree, PersistentSuffixTree pst, HeapStorageProvider storage) : IDisposable
    {
        public void Dispose() => (tree as IDisposable)?.Dispose();
    }

    private static TreeHolder BuildCompact(string text)
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        var ts = new StringTextSource(text);
        long root = builder.Build(ts);
        var pst = new PersistentSuffixTree(storage, root, ts, NodeLayout.Compact);
        return new TreeHolder(pst, pst, storage);
    }

    private static TreeHolder BuildHybrid(string text, long compactOffsetLimit)
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = compactOffsetLimit;
        var ts = new StringTextSource(text);
        long root = builder.Build(ts);
        var pst = new PersistentSuffixTree(storage, root, ts, NodeLayout.Compact,
            builder.TransitionOffset, builder.JumpTableStart, builder.JumpTableEnd,
            builder.DeepestInternalNodeOffset);
        return new TreeHolder(pst, pst, storage);
    }

    // ──────────── C11: Cross-zone suffix links must be followed correctly ──────────

    [TestCase("abcabcabc", 200)]
    [TestCase("mississippi", 200)]
    [TestCase("banana", 250)]
    [TestCase("aababcabcd", 200)]
    public void HybridTree_CrossZoneSuffixLinks_ProduceCorrectLCS(string text, int limit)
    {
        // LCS relies on suffix links.  If cross-zone suffix links are broken,
        // LCS will produce wrong results compared to a pure compact tree.
        using var compact = BuildCompact(text);
        using var hybrid = BuildHybrid(text, limit);

        Assert.That(((PersistentSuffixTree)hybrid.tree).IsHybrid, Is.True, "Tree must be hybrid for this test");

        string compactLcs = compact.tree.LongestCommonSubstring(text);
        string hybridLcs = hybrid.tree.LongestCommonSubstring(text);

        Assert.That(hybridLcs, Is.EqualTo(compactLcs),
            "C11: Hybrid cross-zone suffix links must produce same LCS as compact");
    }

    // ──────────── Visitor for structural comparison ──────────────

    private class CollectingVisitor : ISuffixTreeVisitor
    {
        private readonly List<(int Start, int End, int Leaves, int Children, int Depth)> _visits;

        public CollectingVisitor(List<(int Start, int End, int Leaves, int Children, int Depth)> visits)
            => _visits = visits;

        public void VisitNode(int start, int end, int leafCount, int childCount, int depth)
            => _visits.Add((start, end, leafCount, childCount, depth));

        public void EnterBranch(int key) { }
        public void ExitBranch() { }
    }
}
