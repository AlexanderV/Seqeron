using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The header page of the SESA container (ADR-003 section 3 as amended on 2026-11-18): every field at its offset, the CRC, the alignment rule, the limits; decoded by the independent <see cref="SesaOracle"/>.</summary>
[TestFixture]
public sealed class SesaFormatTests
{
    private static SesaHeader Header() => new()
    {
        TextLength = 0x0102030405060708,
        NodeCount = 0x1112131415161718,
        FileSize = 0x2122232425262728,
        LrsDepth = 0x3132333435363738,
        LrsRank = 0x4142434445464748,
        PosWidth = 4,
        Sigma = 0x51525354,
        TextSha256 = [.. Enumerable.Range(0, 32).Select(i => (byte)(0xA0 + i))],
        BuildFlags = 1,
    };

    private static SesaSection Section(int k) => new(1 + k, k % 4, 4096L * (k + 1), 100 + k, 200 + k, 4 + k, 5 + k, 0x6162636465666768 + k, 0x7172737475767778 + k, 0x81828384u + (uint)k);

    private static SesaOracle Decode(byte[] page) => SesaOracle.FromBytes(page);

    [Test]
    public void EveryHeaderField_IsAtItsOffsetOfADR003()
    {
        SesaOracle o = Decode(SesaFormat.EncodeHeaderPage(Header(), [Section(0), Section(1)], true));
        Assert.That(o.Page, Has.Length.EqualTo(4096));
        Assert.That(o.Magic, Is.EqualTo(new byte[] { 0x89, 0x53, 0x45, 0x53, 0x41, 0x0D, 0x0A, 0x1A }));
        Assert.That((o.VersionMajor, o.VersionMinor), Is.EqualTo((1, 0)));
        Assert.That(o.Page.AsSpan(12, 4).ToArray(), Is.EqualTo(new byte[] { 4, 3, 2, 1 }));
        Assert.That(o.TextLen, Is.EqualTo(0x0102030405060708));
        Assert.That(o.NodeCount, Is.EqualTo(0x1112131415161718));
        Assert.That(o.FileSize, Is.EqualTo(0x2122232425262728));
        Assert.That(o.TableOffset, Is.EqualTo(128));
        Assert.That(o.SectionCount, Is.EqualTo(2));
        Assert.That(o.EntrySize, Is.EqualTo(64));
        Assert.That(o.LrsDepth, Is.EqualTo(0x3132333435363738));
        Assert.That(o.LrsRank, Is.EqualTo(0x4142434445464748));
        Assert.That(o.PosWidth, Is.EqualTo(4));
        Assert.That(o.Sigma, Is.EqualTo(0x51525354));
        Assert.That(o.TextSha256, Is.EqualTo(Enumerable.Range(0, 32).Select(i => (byte)(0xA0 + i)).ToArray()));
        Assert.That(o.BuildFlags, Is.EqualTo(1));
        Assert.That(o.HeaderReserved, Is.Zero);
    }

    [Test]
    public void EverySectionEntryField_IsAtItsOffsetOfADR003_InTheGivenOrder()
    {
        SesaSection[] sections = [Section(0), Section(1), Section(2)];
        SesaOracle o = Decode(SesaFormat.EncodeHeaderPage(Header(), sections, true));
        List<OracleSection> got = o.Sections;
        Assert.That(got, Has.Count.EqualTo(3));
        for (int k = 0; k < 3; k++)
        {
            SesaSection e = sections[k];
            OracleSection g = got[k];
            Assert.That((g.Kind, g.Flags, g.Offset, g.Length, g.Count, g.ElemWidth, g.Codec, g.Param0, g.Param1, g.Crc, g.Reserved), Is.EqualTo((e.Kind, e.Flags, e.Offset, e.Length, e.Count, e.ElemWidth, e.Codec, e.Param0, e.Param1, e.Crc, 0u)), $"entry {k}");
        }
    }

    [Test]
    public void TheHeaderCrc_CoversTheFirst124Bytes_AndTheWholeTable_IncludingTheMagic()
    {
        SesaSection[] sections = [Section(0), Section(1), Section(2), Section(3)];
        byte[] page = SesaFormat.EncodeHeaderPage(Header(), sections, true);
        SesaOracle o = Decode(page);
        Assert.That(o.HeaderCrc, Is.EqualTo(o.ExpectedHeaderCrc()));
        Assert.That(o.HeaderCrc, Is.Not.Zero);
        // every single bit of the header or of the table changes the CRC that the text defines
        for (int i = 0; i < 128 + (64 * 4); i++)
        {
            if (i is >= 124 and < 128) continue;
            byte[] copy = (byte[])page.Clone();
            copy[i] ^= 0x10;
            SesaOracle c = Decode(copy);
            Assert.That(c.ExpectedHeaderCrc(), Is.Not.EqualTo(o.HeaderCrc), $"byte {i}");
        }
    }

    [Test]
    public void WithoutTheMagic_TheFirstEightBytesAreZero_AndTheRestAndTheCrcAreTheSame()
    {
        SesaSection[] sections = [Section(0), Section(1)];
        byte[] withMagic = SesaFormat.EncodeHeaderPage(Header(), sections, true);
        byte[] without = SesaFormat.EncodeHeaderPage(Header(), sections, false);
        Assert.That(without.AsSpan(0, 8).ToArray(), Is.EqualTo(new byte[8]));
        Assert.That(without.AsSpan(8).SequenceEqual(withMagic.AsSpan(8)), Is.True);
        Assert.That(SesaFormat.Magic.ToArray(), Is.EqualTo(withMagic.AsSpan(0, 8).ToArray()));
    }

    [Test]
    public void ThePaddingAfterTheTable_IsZero()
    {
        for (int count = 0; count <= 62; count += 31)
        {
            SesaSection[] sections = [.. Enumerable.Range(0, count).Select(Section)];
            byte[] page = SesaFormat.EncodeHeaderPage(Header(), sections, true);
            Assert.That(page.AsSpan(128 + (64 * count)).IndexOfAnyExcept((byte)0), Is.EqualTo(-1), $"{count} sections");
        }
    }

    [Test]
    public void SixtyTwoSectionsFit_SixtyThreeAreRefused()
    {
        Assert.That(SesaFormat.EncodeHeaderPage(Header(), [.. Enumerable.Range(0, 62).Select(Section)], true), Has.Length.EqualTo(4096));
        Assert.Throws<ArgumentException>(() => SesaFormat.EncodeHeaderPage(Header(), [.. Enumerable.Range(0, 63).Select(Section)], true));
    }

    [Test]
    public void AHashThatIsNot32Bytes_IsRefused()
    {
        Assert.Throws<ArgumentException>(() => SesaFormat.EncodeHeaderPage(Header() with { TextSha256 = new byte[31] }, [], true));
        Assert.Throws<ArgumentException>(() => SesaFormat.EncodeHeaderPage(Header() with { TextSha256 = new byte[33] }, [], true));
    }

    [Test]
    public void TheDefaults_AreTheValuesOfAFileThatHasNoRepeat()
    {
        var h = new SesaHeader();
        Assert.That((h.LrsRank, h.PosWidth, h.BuildFlags), Is.EqualTo((-1L, 4, 1L)));
        Assert.That(h.TextSha256, Has.Length.EqualTo(32));
    }

    [TestCase(4096L, 100L, 4096L)]
    [TestCase(4097L, 100L, 4160L)]
    [TestCase(4160L, 0L, 4160L)]
    [TestCase(4161L, 0L, 4224L)]
    [TestCase(4096L, 1048575L, 4096L)]
    [TestCase(4100L, 1048575L, 4160L)]
    [TestCase(4100L, 1048576L, 8192L)]
    [TestCase(8192L, 1048576L, 8192L)]
    [TestCase(8193L, 5000000000L, 12288L)]
    [TestCase(5000000001L, 5000000000L, 5000003584L)]
    [TestCase(0L, 64L, 0L)]
    [TestCase(1L, 64L, 64L)]
    public void ASectionStartsAt64BytesOrAt4096BytesFromOneMiB(long offset, long length, long expected)
    {
        Assert.That(SesaFormat.AlignedOffset(offset, length), Is.EqualTo(expected));
    }

    [Test]
    public void TheConstants_AreThoseOfTheTextOfADR003()
    {
        Assert.That((SesaFormat.HeaderBytes, SesaFormat.PageBytes, SesaFormat.EntryBytes, SesaFormat.MaxSections), Is.EqualTo((128, 4096, 64, 62)));
        Assert.That((SesaFormat.KindText, SesaFormat.KindAlphabet, SesaFormat.KindSa, SesaFormat.KindLcp, SesaFormat.KindLcpOverflow, SesaFormat.KindChild), Is.EqualTo((1, 2, 3, 4, 5, 6)));
        Assert.That((SesaFormat.KindBuildInfo, SesaFormat.KindPsi, SesaFormat.KindCounts), Is.EqualTo((11, 12, 13)));
        Assert.That((SesaFormat.FlagRequired, SesaFormat.FlagHasCrc), Is.EqualTo((1, 2)));
        Assert.That(SesaFormat.EndianProbe, Is.EqualTo(0x01020304u));
        Assert.That(SesaFormat.CrcOffset, Is.EqualTo(124));
    }
}
