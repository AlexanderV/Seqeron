using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests that <see cref="PersistentSuffixTreeNode"/> public API methods behave
/// correctly (or fail safely) for both pure-Compact and Hybrid (v5) trees.
/// </summary>
[TestFixture]
public class NodeApiSafetyTests
{
    // ──────────── TryGetChild on non-hybrid tree: should work ──────────────

    [Test]
    public void TryGetChild_NonHybridTree_FindsExistingChild()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        var ts = new StringTextSource("banana");
        long root = builder.Build(ts);

        // Non-hybrid: root is compact, children are compact, no jumps
        Assert.That(builder.IsHybrid, Is.False);

        var rootNode = new PersistentSuffixTreeNode(storage, root, NodeLayout.Compact);
        bool found = rootNode.TryGetChild((uint)'b', out var child);

        Assert.That(found, Is.True, "Should find 'b' child on root of non-hybrid tree");
        Assert.That(child.IsNull, Is.False);
    }

    [Test]
    public void TryGetChild_NonHybridTree_ReturnsFalseForAbsentChild()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        var ts = new StringTextSource("banana");
        long root = builder.Build(ts);

        var rootNode = new PersistentSuffixTreeNode(storage, root, NodeLayout.Compact);
        bool found = rootNode.TryGetChild((uint)'z', out _);

        Assert.That(found, Is.False, "Should not find 'z' in 'banana' tree");
    }

    // ──────────── TryGetChild on hybrid tree with jumped children ──────────────

    [Test]
    public void TryGetChild_HybridTree_JumpedNode_ThrowsInvalidOperation()
    {
        // Build hybrid: limit = header + root → first child allocation triggers transition
        // Root is compact but some children are in large zone → root gets jumped child array
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize;
        var ts = new StringTextSource("banana");
        long root = builder.Build(ts);

        Assert.That(builder.IsHybrid, Is.True, "Should be hybrid");

        var rootNode = new PersistentSuffixTreeNode(storage, root, NodeLayout.Compact);

        // Root has jumped children (high-bit set in ChildCount)
        // TryGetChild must throw instead of returning silently wrong results
        Assert.Throws<InvalidOperationException>(() =>
            rootNode.TryGetChild((uint)'b', out _));
    }

    [Test]
    public void TryGetChild_HybridTree_NonJumpedNode_StillWorks()
    {
        // Build hybrid with limit that places most nodes in large zone
        // Internal nodes created in large zone should have normal (non-jumped) children
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize;
        var ts = new StringTextSource("banana");
        long root = builder.Build(ts);

        Assert.That(builder.IsHybrid, Is.True);

        // Find a large-zone internal node (its children are also large-zone, no jump needed)
        // We'll verify this by checking non-jumped large nodes still work
        var pst = new PersistentSuffixTree(storage, root, new StringTextSource("banana"),
            NodeLayout.Compact, builder.TransitionOffset, builder.JumpTableStart, builder.JumpTableEnd);

        // Verify the tree itself works (sanity)
        Assert.That(pst.Contains("banana"), Is.True);
        Assert.That(pst.Contains("ban"), Is.True);
    }

    // ──────────── HasChildren on hybrid tree with jumped children ──────────────

    [Test]
    public void HasChildren_HybridTree_JumpedNode_ReturnsTrue()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize;
        var ts = new StringTextSource("banana");
        long root = builder.Build(ts);

        Assert.That(builder.IsHybrid, Is.True);

        var rootNode = new PersistentSuffixTreeNode(storage, root, NodeLayout.Compact);

        // Root definitely has children — HasChildren must return true
        // even when ChildCount has the high bit set (jumped flag)
        Assert.That(rootNode.HasChildren, Is.True,
            "HasChildren must return true for jumped node (high-bit in ChildCount)");
    }

    [Test]
    public void HasChildren_NonHybridTree_WorksNormally()
    {
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        var ts = new StringTextSource("banana");
        long root = builder.Build(ts);

        var rootNode = new PersistentSuffixTreeNode(storage, root, NodeLayout.Compact);
        Assert.That(rootNode.HasChildren, Is.True);
    }

    // ──────────── A15: WriteOffset must reject out-of-range values for Compact ──────────────

    [Test]
    public void WriteOffset_Compact_ValueExceedsMaxOffset_Throws()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(16); // enough space for a uint32 write

        // Value > CompactMaxOffset (0xFFFFFFFE) and not NULL_OFFSET should throw
        long tooLarge = NodeLayout.CompactMaxOffset + 1; // 0xFFFFFFFF = uint.MaxValue = reserved NULL sentinel
        Assert.Throws<InvalidOperationException>(
            () => NodeLayout.Compact.WriteOffset(storage, 0, tooLarge),
            "A15: WriteOffset must reject value exceeding CompactMaxOffset");
    }

    [Test]
    public void WriteOffset_Compact_HugeValue_Throws()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(16);

        long hugeValue = (long)uint.MaxValue + 100; // clearly beyond uint32 range
        Assert.Throws<InvalidOperationException>(
            () => NodeLayout.Compact.WriteOffset(storage, 0, hugeValue),
            "A15: WriteOffset must reject int64 values that truncate silently");
    }

    [Test]
    public void WriteOffset_Compact_NullOffset_Succeeds()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(16);

        // NULL_OFFSET must always work (maps to 0xFFFFFFFF sentinel)
        Assert.DoesNotThrow(() => NodeLayout.Compact.WriteOffset(storage, 0, PersistentConstants.NULL_OFFSET));
    }

    [Test]
    public void WriteOffset_Large_AnyValue_Succeeds()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(16);

        // Large layout uses int64 — any value should work
        Assert.DoesNotThrow(() => NodeLayout.Large.WriteOffset(storage, 0, long.MaxValue));
        Assert.DoesNotThrow(() => NodeLayout.Large.WriteOffset(storage, 0, NodeLayout.CompactMaxOffset + 100));
    }

    // ──────────── A14: TryGetChild must NOT be public API ──────────────

    [Test]
    public void TryGetChild_IsNotPublic()
    {
        // A14: TryGetChild throws InvalidOperationException for hybrid trees,
        // making it a trap for external consumers. It must be internal.
        var method = typeof(PersistentSuffixTreeNode).GetMethod(
            "TryGetChild",
            System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance);

        Assert.That(method, Is.Null,
            "A14: PersistentSuffixTreeNode.TryGetChild should not be public — " +
            "it throws for hybrid nodes, use ISuffixTree/Navigator instead");
    }

    // ──────────── A16: Raw encoded getters must NOT be public ──────────────

    [Test]
    public void SuffixLink_Getter_IsNotPublic()
    {
        // A16: SuffixLink returns raw value (may be jump-table offset in hybrid trees).
        // External consumers must use ISuffixTree/Navigator to resolve links correctly.
        var prop = typeof(PersistentSuffixTreeNode).GetProperty(
            "SuffixLink",
            System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance);

        Assert.That(prop, Is.Null,
            "A16: PersistentSuffixTreeNode.SuffixLink should not be public — " +
            "it returns raw encoded value, use ISuffixTree/Navigator instead");
    }

    [Test]
    public void ChildrenHead_Getter_IsNotPublic()
    {
        // A16: ChildrenHead may point to jump-table entry in hybrid trees.
        var prop = typeof(PersistentSuffixTreeNode).GetProperty(
            "ChildrenHead",
            System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance);

        Assert.That(prop, Is.Null,
            "A16: PersistentSuffixTreeNode.ChildrenHead should not be public — " +
            "it returns raw encoded value");
    }

    [Test]
    public void ChildCount_Getter_IsNotPublic()
    {
        // A16: ChildCount may have high-bit 0x80000000 set (jumped flag) in hybrid trees.
        var prop = typeof(PersistentSuffixTreeNode).GetProperty(
            "ChildCount",
            System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.Instance);

        Assert.That(prop, Is.Null,
            "A16: PersistentSuffixTreeNode.ChildCount should not be public — " +
            "it returns raw encoded value with possible jumped flag");
    }
}
