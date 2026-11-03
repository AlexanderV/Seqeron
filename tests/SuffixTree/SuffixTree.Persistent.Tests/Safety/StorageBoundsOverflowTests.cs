namespace SuffixTree.Persistent.Tests.Safety;

/// <summary>
/// ST-S0 WP 0.5: the bounds checks of the storage providers computed <c>offset + size</c> in 64-bit arithmetic that
/// wraps for offsets within <c>size</c> of <see cref="long.MaxValue"/>. A corrupt offset read from a hostile file
/// then passed the check and was dereferenced (an <see cref="AccessViolationException"/>/NullReference in the
/// memory-mapped provider). The check must be written without overflow.
/// </summary>
[TestFixture]
[Category("Safety")]
public class StorageBoundsOverflowTests
{
    private string _tempFile = string.Empty;

    [SetUp]
    public void SetUp() => _tempFile = Path.Combine(Path.GetTempPath(), $"bounds_{Guid.NewGuid():N}.bin");

    [TearDown]
    public void TearDown()
    {
        try { if (File.Exists(_tempFile)) File.Delete(_tempFile); } catch (IOException) { }
    }

    private static IEnumerable<long> WrappingOffsets()
    {
        yield return long.MaxValue;
        yield return long.MaxValue - 1;
        yield return long.MaxValue - 3;
        yield return long.MaxValue - 7;
        yield return long.MaxValue - 63;
    }

    private static void AssertAllAccessorsRejectOffset(IStorageProvider p, long offset)
    {
        var bytes = new byte[8];
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadInt32(offset))!.Message, Does.Contain("exceeds logical size"), nameof(p.ReadInt32));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadUInt32(offset))!.Message, Does.Contain("exceeds logical size"), nameof(p.ReadUInt32));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadInt64(offset))!.Message, Does.Contain("exceeds logical size"), nameof(p.ReadInt64));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadChar(offset))!.Message, Does.Contain("exceeds logical size"), nameof(p.ReadChar));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.WriteInt32(offset, 1))!.Message, Does.Contain("exceeds logical size"), nameof(p.WriteInt32));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.WriteUInt32(offset, 1))!.Message, Does.Contain("exceeds logical size"), nameof(p.WriteUInt32));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.WriteInt64(offset, 1))!.Message, Does.Contain("exceeds logical size"), nameof(p.WriteInt64));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.WriteChar(offset, 'x'))!.Message, Does.Contain("exceeds logical size"), nameof(p.WriteChar));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadBytes(offset, bytes, 0, 8))!.Message, Does.Contain("exceeds logical size"), nameof(p.ReadBytes));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => p.WriteBytes(offset, bytes, 0, 8))!.Message, Does.Contain("exceeds logical size"), nameof(p.WriteBytes));
        });
    }

    [TestCaseSource(nameof(WrappingOffsets))]
    public void HeapProvider_OffsetNearLongMaxValue_IsRejectedByTheBoundsCheck(long offset)
    {
        using var p = new HeapStorageProvider(128);
        p.Allocate(100);

        AssertAllAccessorsRejectOffset(p, offset);
    }

    [TestCaseSource(nameof(WrappingOffsets))]
    public void MappedProvider_OffsetNearLongMaxValue_IsRejectedByTheBoundsCheck(long offset)
    {
        using var p = new MappedFileStorageProvider(_tempFile, 4096);
        p.Allocate(100);

        AssertAllAccessorsRejectOffset(p, offset);
    }

    [Test]
    public void Providers_NegativeCount_IsRejected()
    {
        using var heap = new HeapStorageProvider(128);
        heap.Allocate(100);
        using var mapped = new MappedFileStorageProvider(_tempFile, 4096);
        mapped.Allocate(100);
        var buffer = new byte[16];

        foreach (IStorageProvider p in new IStorageProvider[] { heap, mapped })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadBytes(50, buffer, 0, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => p.WriteBytes(50, buffer, 0, -1));
        }
    }

    [Test]
    public void Providers_AccessesAtTheLogicalEnd_StillWork()
    {
        // The rewritten check must keep the exact boundary: the last 4 bytes are readable, the next byte is not.
        foreach (IStorageProvider p in new IStorageProvider[] { new HeapStorageProvider(128), new MappedFileStorageProvider(_tempFile, 4096) })
        {
            using (p)
            {
                p.Allocate(100);
                p.WriteInt32(96, 0x1234);
                Assert.That(p.ReadInt32(96), Is.EqualTo(0x1234));
                Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadInt32(97));
                Assert.That(p.ReadInt64(92), Is.EqualTo(0x1234L << 32));
                Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadInt64(93));
                Assert.That(p.ReadChar(98), Is.EqualTo('\0'));
                Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadChar(99));
                Assert.Throws<ArgumentOutOfRangeException>(() => p.ReadInt32(-1));
            }
            TearDown();
        }
    }
}
