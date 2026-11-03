using System.Buffers.Binary;
using System.Text;

namespace SuffixTree.Persistent.Tests.Format;

/// <summary>
/// ST-S0 WP 0.8: pins the on-disk magic number. The constant was documented as the ASCII string "SUFFIXTR" while
/// its value <c>0x5452454558494646</c> spells "TREEXIFF" (most-significant byte first) and is stored as the bytes
/// "FFIXEERT". The format stays unchanged; these tests make the documentation and the code unable to drift apart again.
/// </summary>
[TestFixture]
[Category("Format")]
public class MagicNumberTests
{
    private static readonly byte[] ExpectedFileBytes = { 0x46, 0x46, 0x49, 0x58, 0x45, 0x45, 0x52, 0x54 };

    [Test]
    public void MagicNumber_HasTheDocumentedValue()
    {
        Assert.That(PersistentConstants.MAGIC_NUMBER, Is.EqualTo(0x5452454558494646L));
    }

    [Test]
    public void MagicNumber_SpellsTreexiffMostSignificantByteFirst_AndIsNotSuffixtr()
    {
        var bigEndian = new byte[8];
        BinaryPrimitives.WriteInt64BigEndian(bigEndian, PersistentConstants.MAGIC_NUMBER);

        Assert.Multiple(() =>
        {
            Assert.That(Encoding.ASCII.GetString(bigEndian), Is.EqualTo("TREEXIFF"));
            Assert.That(Encoding.ASCII.GetString(bigEndian), Is.Not.EqualTo("SUFFIXTR"));
        });
    }

    [Test]
    public void MagicNumber_IsStoredAsTheBytes46_46_49_58_45_45_52_54()
    {
        var littleEndian = new byte[8];
        BinaryPrimitives.WriteInt64LittleEndian(littleEndian, PersistentConstants.MAGIC_NUMBER);

        Assert.That(littleEndian, Is.EqualTo(ExpectedFileBytes));
        Assert.That(Encoding.ASCII.GetString(littleEndian), Is.EqualTo("FFIXEERT"));
    }

    [Test]
    public void HeapBuiltTree_StartsWithTheMagicBytes()
    {
        using var storage = new HeapStorageProvider();
        new PersistentSuffixTreeBuilder(storage).Build(new StringTextSource("banana"));

        var first8 = new byte[8];
        storage.ReadBytes(0, first8, 0, 8);

        Assert.That(first8, Is.EqualTo(ExpectedFileBytes));
    }

    [Test]
    public void TreeFile_StartsWithTheMagicBytes_AndLoadsBack()
    {
        string path = Path.Combine(Path.GetTempPath(), $"magic_{Guid.NewGuid():N}.tree");
        try
        {
            PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("banana"), path).Dispose();

            byte[] file = File.ReadAllBytes(path);
            Assert.That(file.AsSpan(0, 8).ToArray(), Is.EqualTo(ExpectedFileBytes));

            // A loader that expected a different magic would refuse this file; conversely the documented bytes are accepted.
            using var tree = PersistentSuffixTreeFactory.LoadPersistent(path);
            Assert.That(tree.Contains("nan"), Is.True);
        }
        finally
        {
            try { File.Delete(path); } catch (IOException) { }
        }
    }

    [Test]
    public void FileWithTheOldDocumentedMagic_IsRejected()
    {
        // The string the documentation used to claim ("SUFFIXTR") is not the format's magic.
        var storage = new HeapStorageProvider();
        storage.Allocate(PersistentConstants.HEADER_SIZE_V6);
        storage.WriteBytes(0, Encoding.ASCII.GetBytes("SUFFIXTR"), 0, 8);

        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage));
        Assert.That(ex!.Message, Does.Contain("Magic"));
    }
}
