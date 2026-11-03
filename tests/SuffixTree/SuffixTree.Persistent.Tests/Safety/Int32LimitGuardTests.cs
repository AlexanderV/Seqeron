using System.Buffers.Binary;
using System.Collections;

namespace SuffixTree.Persistent.Tests.Safety;

/// <summary>
/// ST-S0 WP 0.2: every Int32 limit of the persistent builder and its storage providers must be reported
/// with a clear exception instead of a wrapped counter, an <see cref="OutOfMemoryException"/> from the
/// runtime, or silent corruption. None of these tests allocates more than a few kilobytes: huge sizes are
/// reported by a fake <see cref="ITextSource"/> / <see cref="IStorageProvider"/> that never serves content.
/// </summary>
[TestFixture]
[Category("Safety")]
public class Int32LimitGuardTests
{
    // 2n + 1 nodes must fit in Int32: the largest accepted text is int.MaxValue / 2 = 1_073_741_823 characters.
    private const int LargestAcceptedLength = 1_073_741_823;
    private const int SmallestRejectedLength = 1_073_741_824;

    /// <summary>Reports a configurable <see cref="Length"/> and throws if anybody touches the content.</summary>
    private sealed class HugeLengthTextSource(int length) : ITextSource
    {
        public int Length { get; } = length;
        public char this[int index] => throw new NotSupportedException("content must not be read");
        public string Substring(int start, int length) => throw new NotSupportedException("content must not be read");
        public ReadOnlySpan<char> Slice(int start, int length) => throw new NotSupportedException("content must not be read");
        public IEnumerator<char> GetEnumerator() => throw new NotSupportedException("content must not be read");
        IEnumerator IEnumerable.GetEnumerator() => throw new NotSupportedException("content must not be read");
    }

    /// <summary>Storage that serves only a hand-written v6 header and reports a multi-gigabyte <see cref="Size"/>.</summary>
    private sealed class HeaderOnlyStorage : IStorageProvider
    {
        private readonly byte[] _header = new byte[PersistentConstants.HEADER_SIZE_V6];

        public HeaderOnlyStorage(long size, int textLength, bool ascii)
        {
            Size = size;
            BinaryPrimitives.WriteInt64LittleEndian(_header.AsSpan(PersistentConstants.HEADER_OFFSET_MAGIC), PersistentConstants.MAGIC_NUMBER);
            BinaryPrimitives.WriteInt32LittleEndian(_header.AsSpan(PersistentConstants.HEADER_OFFSET_VERSION), 6);
            BinaryPrimitives.WriteInt32LittleEndian(_header.AsSpan(PersistentConstants.HEADER_OFFSET_TEXT_LEN), textLength);
            BinaryPrimitives.WriteInt64LittleEndian(_header.AsSpan(PersistentConstants.HEADER_OFFSET_ROOT), PersistentConstants.HEADER_SIZE_V6);
            BinaryPrimitives.WriteInt64LittleEndian(_header.AsSpan(PersistentConstants.HEADER_OFFSET_TEXT_OFF), PersistentConstants.HEADER_SIZE_V6);
            BinaryPrimitives.WriteInt32LittleEndian(_header.AsSpan(PersistentConstants.HEADER_OFFSET_FLAGS), ascii ? PersistentConstants.FLAG_TEXT_ASCII : 0);
            BinaryPrimitives.WriteInt64LittleEndian(_header.AsSpan(PersistentConstants.HEADER_OFFSET_TRANSITION), -1);
            BinaryPrimitives.WriteInt64LittleEndian(_header.AsSpan(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE), PersistentConstants.NULL_OFFSET);
        }

        public long Size { get; }
        public void EnsureCapacity(long capacity) { }
        public int ReadInt32(long offset) => BinaryPrimitives.ReadInt32LittleEndian(_header.AsSpan((int)offset));
        public void WriteInt32(long offset, int value) => throw new NotSupportedException();
        public uint ReadUInt32(long offset) => BinaryPrimitives.ReadUInt32LittleEndian(_header.AsSpan((int)offset));
        public void WriteUInt32(long offset, uint value) => throw new NotSupportedException();
        public long ReadInt64(long offset) => BinaryPrimitives.ReadInt64LittleEndian(_header.AsSpan((int)offset));
        public void WriteInt64(long offset, long value) => throw new NotSupportedException();
        public char ReadChar(long offset) => throw new NotSupportedException();
        public void WriteChar(long offset, char value) => throw new NotSupportedException();
        public void ReadBytes(long offset, byte[] buffer, int start, int count) => throw new NotSupportedException("content must not be read");
        public void WriteBytes(long offset, byte[] buffer, int start, int count) => throw new NotSupportedException();
        public long Allocate(int size) => throw new NotSupportedException();
        public void Dispose() { }
    }

    // ─── Text length: 2n + 1 nodes must fit in Int32 ─────────────────

    [Test]
    public void Builder_Build_TextOneBeyondLimit_ThrowsArgumentException_BeforeReadingContent()
    {
        using var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);

        var ex = Assert.Throws<ArgumentException>(() => builder.Build(new HugeLengthTextSource(SmallestRejectedLength)));

        Assert.That(ex!.ParamName, Is.EqualTo("text"));
        Assert.That(ex.Message, Does.Contain(SmallestRejectedLength.ToString()).And.Contain("1073741823"));
    }

    [Test]
    public void Builder_Build_TextOfIntMaxValue_ThrowsArgumentException()
    {
        using var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);

        Assert.Throws<ArgumentException>(() => builder.Build(new HugeLengthTextSource(int.MaxValue)));
    }

    [Test]
    public void Factory_Create_HeapMode_TextBeyondLimit_ThrowsArgumentException_WithoutAllocatingStorage()
    {
        var ex = Assert.Throws<ArgumentException>(() =>
            PersistentSuffixTreeFactory.Create(new HugeLengthTextSource(SmallestRejectedLength)));

        Assert.That(ex!.ParamName, Is.EqualTo("text"));
    }

    [Test]
    public void MaxTextLength_IsTheLargestLengthWhoseWorstCaseNodeCountFitsInt32()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SuffixTreeLimits.MaxTextLength, Is.EqualTo(LargestAcceptedLength));
            Assert.That(2L * SuffixTreeLimits.MaxTextLength + 1, Is.LessThanOrEqualTo(int.MaxValue));
            Assert.That(2L * (SuffixTreeLimits.MaxTextLength + 1L) + 1, Is.GreaterThan(int.MaxValue));
        });
    }

    [Test]
    public void Builder_Build_TextExactlyAtLimit_PassesTheLengthGuard()
    {
        // The guard must not reject the largest legal length: the build proceeds and only then trips over
        // the fake source's refusal to serve content.
        using var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);

        Assert.Throws<NotSupportedException>(() => builder.Build(new HugeLengthTextSource(LargestAcceptedLength)));
    }

    [Test]
    public void Builder_Build_RejectedLength_DoesNotConsumeTheBuilder()
    {
        using var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage);

        Assert.Throws<ArgumentException>(() => builder.Build(new HugeLengthTextSource(SmallestRejectedLength)));

        // A rejected argument is not a started build: the same builder can still index a valid text.
        Assert.DoesNotThrow(() => builder.Build(new StringTextSource("banana")));
    }

    // ─── Node counter: _nodeCount++ must never wrap ──────────────────

    [Test]
    public void Builder_NodeCountLimit_ExactlyEnough_Builds_OneLess_ThrowsInvalidOperationException()
    {
        const string text = "mississippi";
        int nodes;
        using (var reference = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text)))
            nodes = reference.NodeCount;

        using (var ok = new HeapStorageProvider())
            Assert.DoesNotThrow(() => new PersistentSuffixTreeBuilder(ok) { NodeCountLimit = nodes }.Build(new StringTextSource(text)));

        using var tooSmall = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(tooSmall) { NodeCountLimit = nodes - 1 };

        var ex = Assert.Throws<InvalidOperationException>(() => builder.Build(new StringTextSource(text)));

        Assert.That(ex!.Message, Does.Contain((nodes - 1).ToString()).And.Contain("nodes"));
    }

    [Test]
    public void Builder_NodeCountLimit_DefaultsToInt32MaxValue()
    {
        using var storage = new HeapStorageProvider();

        Assert.That(new PersistentSuffixTreeBuilder(storage).NodeCountLimit, Is.EqualTo(int.MaxValue));
    }

    // ─── Heap mode: the byte[] behind HeapStorageProvider is capped by the runtime ───

    [Test]
    public void Factory_Create_HeapMode_TextTooLargeForManagedArray_ThrowsClearArgumentException()
    {
        // 40 M characters is below the 2n+1 limit but needs an estimated 72 B/char = 2.9 GB of storage,
        // more than one managed byte[] can hold. The factory used to request a 2 GB array and fail
        // with an OutOfMemoryException from the runtime.
        var ex = Assert.Throws<ArgumentException>(() =>
            PersistentSuffixTreeFactory.Create(new HugeLengthTextSource(40_000_000)));

        Assert.That(ex!.ParamName, Is.EqualTo("text"));
        Assert.That(ex.Message, Does.Contain("filePath"));
    }

    [Test]
    public void HeapProvider_EnsureCapacity_AtIntMaxValue_ThrowsInvalidOperationException()
    {
        // int.MaxValue is larger than the maximum .NET array length (Array.MaxLength = 0x7FFFFFC7):
        // the provider used to hand int.MaxValue to Array.Resize and surface an OutOfMemoryException.
        using var p = new HeapStorageProvider(128);

        var ex = Assert.Throws<InvalidOperationException>(() => p.EnsureCapacity(int.MaxValue));

        Assert.That(ex!.Message, Does.Contain("HeapStorageProvider"));
    }

    [Test]
    public void HeapProvider_EnsureCapacity_JustAboveArrayMaxLength_ThrowsInvalidOperationException()
    {
        using var p = new HeapStorageProvider(128);

        Assert.Throws<InvalidOperationException>(() => p.EnsureCapacity((long)Array.MaxLength + 1));
    }

    [Test]
    public void HeapProvider_Constructor_InitialCapacityAboveArrayMaxLength_ThrowsArgumentOutOfRange()
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => new HeapStorageProvider(int.MaxValue));

        Assert.That(ex!.ParamName, Is.EqualTo("initialCapacity"));
    }

    [Test]
    public void HeapProvider_FailedGrowth_LeavesProviderUsable()
    {
        using var p = new HeapStorageProvider(128);
        p.Allocate(64);
        p.WriteInt32(0, 0x1234);

        Assert.Throws<InvalidOperationException>(() => p.EnsureCapacity(int.MaxValue));

        Assert.That(p.Size, Is.EqualTo(64));
        Assert.That(p.ReadInt32(0), Is.EqualTo(0x1234));
        Assert.That(p.Allocate(16), Is.EqualTo(64));
    }

    // ─── Load into a non-memory-mapped storage: the text must fit in one string / byte[] ───

    [Test]
    public void Load_NonMmfStorage_AsciiTextLongerThanByteArrayLimit_ThrowsInvalidOperationException()
    {
        // 2_147_483_600 ASCII bytes <= int.MaxValue (so the old check passed) but > Array.MaxLength:
        // ArrayPool.Rent used to fail with an OutOfMemoryException.
        using var storage = new HeaderOnlyStorage(size: 4_000_000_000L, textLength: 2_147_483_600, ascii: true);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));

        Assert.That(ex!.Message, Does.Contain("2147483600").And.Contain("memory-mapped"));
    }

    [Test]
    public void Load_NonMmfStorage_Utf16TextLongerThanByteArrayLimit_ThrowsInvalidOperationException()
    {
        // 1_073_741_823 UTF-16 chars = 2_147_483_646 bytes > Array.MaxLength (2_147_483_591).
        using var storage = new HeaderOnlyStorage(size: 4_000_000_000L, textLength: 1_073_741_823, ascii: false);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));

        Assert.That(ex!.Message, Does.Contain("1073741823").And.Contain("memory-mapped"));
    }

    [Test]
    public void Load_NonMmfStorage_Utf16TextBeyondInt32Bytes_ThrowsInvalidOperationException()
    {
        // 1_500_000_000 UTF-16 chars = 3e9 bytes > int.MaxValue (this limit already existed; characterization).
        using var storage = new HeaderOnlyStorage(size: 4_000_000_000L, textLength: 1_500_000_000, ascii: false);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));

        Assert.That(ex!.Message, Does.Contain("1500000000").And.Contain("memory-mapped"));
    }
}
