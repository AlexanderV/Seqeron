using NUnit.Framework;
using SuffixTree;

namespace SuffixTree.Persistent.Tests;

/// <summary>
/// Tests that <see cref="PersistentSuffixTree.Load"/> validates header fields
/// and rejects corrupted or truncated storage with clear exceptions.
/// </summary>
[TestFixture]
public class LoadValidationTests
{
    // ──────────── Truncated storage ──────────────

    [Test]
    public void Load_StorageSmallerThanHeader_ThrowsInvalidOperation()
    {
        // Storage with only 10 bytes — too small for any header (48 bytes min)
        var storage = new HeapStorageProvider(64);
        storage.Allocate(10); // advance position to 10 bytes

        Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
    }

    // ──────────── Bad MAGIC ──────────────

    [Test]
    public void Load_GarbageMagic_ThrowsInvalidOperation()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE); // 48 bytes
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, unchecked((long)0xDEADBEEFCAFEBABE));

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("Magic"));
    }

    // ──────────── Unsupported version ──────────────

    [Test]
    public void Load_UnsupportedVersion_Throws()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, PersistentConstants.MAGIC_NUMBER);
        storage.WriteInt32(PersistentConstants.HEADER_OFFSET_VERSION, 99);

        Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
    }

    // ──────────── Root offset beyond storage ──────────────

    [Test]
    public void Load_RootOffsetBeyondStorage_ThrowsInvalidOperation()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE);
        WriteValidHeader(storage, version: 4);
        // Root points far beyond actual storage
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, 999_999);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("root"));
    }

    [Test]
    public void Load_NegativeRootOffset_ThrowsInvalidOperation()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE);
        WriteValidHeader(storage, version: 4);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, -42);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("root"));
    }

    // ──────────── S11: TEXT_OFF / TEXT_LEN validation ──────────────

    [Test]
    public void Load_TextOffBeyondStorage_ThrowsInvalidOperation()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE + 100);
        WriteValidHeader(storage, version: 4);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, PersistentConstants.HEADER_SIZE);
        // TEXT_OFF points far beyond storage
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF, 999_999);
        storage.WriteInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN, 5);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("text").IgnoreCase);
    }

    [Test]
    public void Load_TextLenNegative_ThrowsInvalidOperation()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE + 100);
        WriteValidHeader(storage, version: 4);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, PersistentConstants.HEADER_SIZE);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF, PersistentConstants.HEADER_SIZE + 10);
        storage.WriteInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN, -42);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("text").IgnoreCase);
    }

    [Test]
    public void Load_TextRegionExceedsStorage_ThrowsInvalidOperation()
    {
        var storage = new HeapStorageProvider();
        int totalSize = PersistentConstants.HEADER_SIZE + 100;
        storage.Allocate(totalSize);
        WriteValidHeader(storage, version: 4);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, PersistentConstants.HEADER_SIZE);
        // TEXT_OFF is valid, but TEXT_OFF + TEXT_LEN * 2 exceeds storage
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TEXT_OFF, totalSize - 4);
        storage.WriteInt32(PersistentConstants.HEADER_OFFSET_TEXT_LEN, 100); // 100 chars = 200 bytes

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("text").IgnoreCase);
    }

    // ──────────── V5 hybrid header field validation ──────────────

    [Test]
    public void Load_V5_TransitionBeyondStorage_ThrowsInvalidOperation()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE_V5 + 100); // some space
        WriteValidHeader(storage, version: 5);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, PersistentConstants.HEADER_SIZE_V5);
        // Transition points beyond storage
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TRANSITION, 999_999);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_START, PersistentConstants.HEADER_SIZE_V5 + 28);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_END, PersistentConstants.HEADER_SIZE_V5 + 36);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("transition").IgnoreCase);
    }

    [Test]
    public void Load_V5_JumpEndBeforeJumpStart_ThrowsInvalidOperation()
    {
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE_V5 + 200);
        WriteValidHeader(storage, version: 5);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_ROOT, PersistentConstants.HEADER_SIZE_V5);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_TRANSITION, PersistentConstants.HEADER_SIZE_V5 + 50);
        // Jump end < jump start — invalid
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_START, PersistentConstants.HEADER_SIZE_V5 + 100);
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_JUMP_END, PersistentConstants.HEADER_SIZE_V5 + 50);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("jump"));
    }

    // ──────────── C16: Header SIZE must match storage.Size ──────────────

    [Test]
    public void Load_HeaderSizeMismatch_ThrowsInvalidOperation()
    {
        // Build a real tree so header is valid
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.Build(new StringTextSource("banana"));

        // Tamper: write a SIZE that differs from actual storage size
        long actualSize = storage.Size;
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_SIZE, actualSize + 1000);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("size").IgnoreCase);
    }

    [Test]
    public void Load_HeaderSizeSmallerThanStorage_ThrowsInvalidOperation()
    {
        // Build a real tree, then set SIZE to a smaller value (simulates file corruption)
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.Build(new StringTextSource("abc"));

        long actualSize = storage.Size;
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_SIZE, actualSize - 10);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("size").IgnoreCase);
    }

    // ──────────── Valid storage should still load ──────────────

    [Test]
    public void Load_ValidCompactTree_LoadsSuccessfully()
    {
        // Build a real tree → export → Load
        var buildStorage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(buildStorage, NodeLayout.Compact);
        long root = builder.Build(new StringTextSource("banana"));

        using var tree = PersistentSuffixTree.Load(buildStorage);
        Assert.That(tree.Contains("banana"), Is.True);
        Assert.That(tree.Contains("ana"), Is.True);
    }

    [Test]
    public void Load_ValidHybridTree_LoadsSuccessfully()
    {
        var buildStorage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(buildStorage, NodeLayout.Compact);
        builder.CompactOffsetLimit = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize;
        long root = builder.Build(new StringTextSource("banana"));

        Assert.That(builder.IsHybrid, Is.True);

        using var tree = PersistentSuffixTree.Load(buildStorage);
        Assert.That(tree.Contains("banana"), Is.True);
        Assert.That(tree.Contains("ana"), Is.True);
    }

    // ──────────── S19+S21: DEEPEST_NODE must be validated ──────────────

    [Test]
    public void Load_DeepestNodeBeyondStorage_ThrowsInvalidOperation()
    {
        // Build a real v5 (hybrid) tree, then tamper DEEPEST_NODE to point beyond storage
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize;
        builder.Build(new StringTextSource("banana"));

        // Tamper: set deepest node to far beyond storage
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE, storage.Size + 9999);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("deepest").IgnoreCase);
    }

    [Test]
    public void Load_DeepestNodeNegative_ThrowsInvalidOperation()
    {
        // Build a real v5 tree, then tamper DEEPEST_NODE to a negative non-NULL value
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.CompactOffsetLimit = PersistentConstants.HEADER_SIZE_V5 + NodeLayout.Compact.NodeSize;
        builder.Build(new StringTextSource("banana"));

        // -42 is not NULL_OFFSET (-1) and not a valid offset
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE, -42);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("deepest").IgnoreCase);
    }

    [Test]
    public void Load_NonHybrid_ReadsDeepestNodeFromV5Header()
    {
        // Build a non-hybrid tree. Builder always writes v5 header now.
        // DEEPEST_NODE at offset 72 must be read and used for O(1) LRS.
        var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        builder.Build(new StringTextSource("banana"));

        // Verify this is v5 (all trees are v5 now)
        int version = storage.ReadInt32(PersistentConstants.HEADER_OFFSET_VERSION);
        Assert.That(version, Is.EqualTo(5), "All trees must use v5 header");

        using var tree = PersistentSuffixTree.Load(storage);
        Assert.That(tree.Contains("banana"), Is.True);
        Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("ana"));
    }

    // ──────────── Helpers ──────────────

    private static void WriteValidHeader(HeapStorageProvider storage, int version)
    {
        storage.WriteInt64(PersistentConstants.HEADER_OFFSET_MAGIC, PersistentConstants.MAGIC_NUMBER);
        storage.WriteInt32(PersistentConstants.HEADER_OFFSET_VERSION, version);
    }
}
