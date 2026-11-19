using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The validation levels of the reader (work package 1b.7f, gate G1; ADR-003 section 6): a flipped byte of the header page is refused at <c>Header</c>, a flipped byte of the body at <c>Sections</c>, and a body that was
/// damaged and sealed again (the section checksums recomputed) at <c>Structure</c>, including the suffix array that is a permutation in the wrong order; a refused file is an <see cref="InvalidDataException"/>
/// (a larger major version a <see cref="NotSupportedException"/>), never a hang, a crash or a read outside the mapping.
/// </summary>
[TestFixture]
public class SesaReaderCorruptionTests
{
    /// <summary>A text with every section: a random part and a copy of 320 of its symbols, so that the LCP has an overflow list.</summary>
    private static readonly Lazy<(string Text, LcpChildCase Case)> WithOverflow = new(() =>
    {
        string random = StringFamilies.Random(1_500, 4, 21);
        string text = random + random[..300];
        return (text, LcpChildCase.Of(text));
    });

    private static void AssertRefused(Action open, string label)
    {
        Exception? caught = Assert.Catch(() => open(), label);
        Assert.That(caught, Is.Not.Null, $"{label}: opened");
        Assert.That(caught is InvalidDataException or NotSupportedException, Is.True, $"{label}: {caught!.GetType().Name}: {caught.Message}");
    }

    private static void OpenAndDispose(string path, EsaValidationLevel level, bool sha = false)
    {
        using EsaSuffixTree tree = SesaTestFile.Open(path, level, EsaAccessMode.Checked, sha);
        Assert.That(tree, Is.Not.Null);
    }

    [Test]
    [CancelAfter(300_000)]
    public void EveryByteOfThePageAndTheTable_Flipped_IsRefusedAtHeader()
    {
        (string text, LcpChildCase c) = WithOverflow.Value;
        byte[] image = SesaTestFile.Image(c);
        Assert.That(SesaOracle.FromBytes(image).SectionCount, Is.EqualTo(9), "every section kind of the container is present");
        using var dir = new TempDir();
        string path = dir.File("flip.sesa");
        foreach (byte mask in new byte[] { 0x01, 0xFF })
        {
            for (int i = 0; i < SesaFormat.PageBytes; i++)
            {
                byte[] copy = (byte[])image.Clone();
                copy[i] ^= mask;
                File.WriteAllBytes(path, copy);
                int at = i;
                AssertRefused(() => OpenAndDispose(path, EsaValidationLevel.Header), $"byte {at} xor {mask:X2} of the header page of '{text.Length}' symbols");
            }
        }
    }

    [Test]
    public void TheFlipsOfTheCrcBytesAndTheMagic_AreRefused_WithTheRightException()
    {
        byte[] image = SesaTestFile.Image(LcpChildCase.Of("banana"));
        using var dir = new TempDir();
        string path = dir.File("m.sesa");
        for (int i = 0; i < 8; i++)
        {
            byte[] copy = (byte[])image.Clone();
            copy[i] ^= 0x80;
            File.WriteAllBytes(path, copy);
            Assert.That(() => OpenAndDispose(path, EsaValidationLevel.Header), Throws.TypeOf<InvalidDataException>(), $"magic byte {i}");
        }
        byte[] crc = (byte[])image.Clone();
        crc[124] ^= 1;
        File.WriteAllBytes(path, crc);
        Assert.That(() => OpenAndDispose(path, EsaValidationLevel.Header), Throws.TypeOf<InvalidDataException>().With.Message.Contains("CRC").IgnoreCase);
    }

    [Test]
    public void EveryTruncation_IsRefused()
    {
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(StringFamilies.Random(200, 3, 4)));
        using var dir = new TempDir();
        string path = dir.File("t.sesa");
        var lengths = new SortedSet<int> { 0, 1, 7, 8, 9, 100, 127, 128, 129, 191, 4095, 4096, 4097, image.Length - 1, image.Length - 64, image.Length / 2 };
        for (int k = 4096; k < image.Length; k += 97) lengths.Add(k);
        foreach (int length in lengths.Where(l => l < image.Length))
        {
            File.WriteAllBytes(path, image.AsSpan(0, length).ToArray());
            int len = length;
            AssertRefused(() => OpenAndDispose(path, EsaValidationLevel.Header), $"truncated to {len} of {image.Length} bytes");
        }
    }

    [Test]
    public void TrailingBytes_AfterTheLastSection_AreRefused()
    {
        byte[] image = SesaTestFile.Image(LcpChildCase.Of("banana"));
        using var dir = new TempDir();
        string path = dir.File("t.sesa");
        File.WriteAllBytes(path, [.. image, .. new byte[64]]);
        Assert.That(() => OpenAndDispose(path, EsaValidationLevel.Header), Throws.TypeOf<InvalidDataException>().With.Message.Contains("size").IgnoreCase);
    }

    [Test]
    [CancelAfter(300_000)]
    public void ThousandRandomBodyFlips_AreCaughtBySections()
    {
        (_, LcpChildCase c) = WithOverflow.Value;
        byte[] image = SesaTestFile.Image(c);
        using var dir = new TempDir();
        string path = dir.File("body.sesa");
        var rng = new SplitMix64(0x1B7F);
        for (int i = 0; i < 1_000; i++)
        {
            byte[] copy = (byte[])image.Clone();
            int at = SesaFormat.PageBytes + rng.NextInt(image.Length - SesaFormat.PageBytes);
            copy[at] ^= (byte)(1 + rng.NextInt(255));
            File.WriteAllBytes(path, copy);
            int where = at;
            AssertRefused(() => OpenAndDispose(path, EsaValidationLevel.Sections), $"flip {i} at byte {where} (file of {image.Length})");
        }
    }

    [Test]
    public void ABodyFlip_IsInvisibleAtHeader_ThatIsTheContractOfTheLevel()
    {
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(StringFamilies.Random(300, 4, 1)));
        OracleSection sa = SesaTestFile.Section(image, SesaFormat.KindSa);
        image[sa.Offset + 40] ^= 0x04;
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Header);
        Assert.That(tree, Is.Not.Null);
        Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Sections, name: "y.sesa"), Throws.TypeOf<InvalidDataException>());
    }

    // ---------------------------------------------------------------- sealed damage: only Structure can see it

    private static (string Name, byte[] Image) Damaged(string name, LcpChildCase c, Action<byte[], SesaOracle> damage, bool psi = true)
    {
        byte[] image = SesaTestFile.Image(c, psi);
        damage(image, SesaOracle.FromBytes(image));
        SesaTestFile.Reseal(image);
        return (name, image);
    }

    private static void AssertOnlyStructureSees(string name, byte[] image, bool headerMaySee = false)
    {
        using var dir = new TempDir();
        string path = dir.File("s.sesa");
        File.WriteAllBytes(path, image);
        if (!headerMaySee)
        {
            OpenAndDispose(path, EsaValidationLevel.Header);
            OpenAndDispose(path, EsaValidationLevel.Sections);
        }
        InvalidDataException? ex = Assert.Throws<InvalidDataException>(() => OpenAndDispose(path, EsaValidationLevel.Structure, sha: true), $"{name}: Structure");
        Assert.That(ex!.Message, Is.Not.Empty, name);
    }

    private static int[] ReadInts(byte[] image, OracleSection s) => System.Runtime.InteropServices.MemoryMarshal.Cast<byte, int>(image.AsSpan((int)s.Offset, (int)s.Length)).ToArray();

    private static void WriteInt(byte[] image, OracleSection s, int index, int value) => BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan((int)s.Offset + (4 * index)), value);

    private static void WriteUInt(byte[] image, OracleSection s, int index, uint value) => BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)s.Offset + (4 * index)), value);

    [TestCase("banana")]
    [TestCase("mississippi")]
    [TestCase("aaaaaaaaaaaa")]
    [TestCase("abababababab")]
    [TestCase("GATTACAGATTACAGGTA")]
    public void TheSuffixArray_WithTwoAdjacentEntriesSwapped_IsCaughtByTheOrderCheck_AtEveryRank(string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        for (int k = 1; k < c.N; k++)
        {
            int rank = k;
            (string name, byte[] image) = Damaged($"swap SA[{rank}] and SA[{rank + 1}] of '{text}'", c, (img, o) =>
            {
                OracleSection sa = o.Get(SesaFormat.KindSa);
                int[] v = ReadInts(img, sa);
                WriteInt(img, sa, rank, v[rank + 1]);
                WriteInt(img, sa, rank + 1, v[rank]);
            });
            AssertOnlyStructureSees(name, image);
        }
    }

    [Test]
    public void TheSuffixArray_WithAWrongOrderPermutation_OfALargerText_IsCaught()
    {
        (_, LcpChildCase c) = WithOverflow.Value;
        var rng = new SplitMix64(5);
        for (int i = 0; i < 25; i++)
        {
            int rank = 1 + rng.NextInt(c.N - 1);
            (string name, byte[] image) = Damaged($"swap at {rank}", c, (img, o) =>
            {
                OracleSection sa = o.Get(SesaFormat.KindSa);
                int[] v = ReadInts(img, sa);
                WriteInt(img, sa, rank, v[rank + 1]);
                WriteInt(img, sa, rank + 1, v[rank]);
            });
            AssertOnlyStructureSees(name, image);
        }
    }

    [Test]
    public void TheSuffixArray_ThatIsNotAPermutation_OrHasNoSentinelFirst_IsCaught()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        var cases = new List<(string, byte[])>
        {
            Damaged("duplicate entry", c, (img, o) => WriteInt(img, o.Get(3), 5, ReadInts(img, o.Get(3))[6])),
            Damaged("entry n + 1", c, (img, o) => WriteInt(img, o.Get(3), 4, c.N + 1)),
            Damaged("entry 0xFFFFFFFF", c, (img, o) => WriteInt(img, o.Get(3), 4, -1)),
            Damaged("entry 0x7FFFFFFF", c, (img, o) => WriteInt(img, o.Get(3), 7, int.MaxValue)),
            Damaged("SA[0] is not n", c, (img, o) =>
            {
                OracleSection sa = o.Get(3);
                int[] v = ReadInts(img, sa);
                int j = Array.IndexOf(v, 0);
                WriteInt(img, sa, 0, 0);
                WriteInt(img, sa, j, c.N);
            }),
        };
        foreach ((string name, byte[] image) in cases) AssertOnlyStructureSees(name, image);
    }

    [Test]
    public void TheLcpArray_WithAWrongValue_IsCaught_EverywhereAndInTheOverflowList()
    {
        (_, LcpChildCase c) = WithOverflow.Value;
        var rng = new SplitMix64(7);
        var cases = new List<(string, byte[])>();
        for (int i = 0; i < 30; i++)
        {
            int rank = 2 + rng.NextInt(c.N - 2);
            int delta = i % 2 == 0 ? 1 : -1;
            cases.Add(Damaged($"LCP[{rank}] {delta:+0;-0}", c, (img, o) =>
            {
                OracleSection lcp = o.Get(SesaFormat.KindLcp);
                int at = (int)lcp.Offset + (rank * lcp.ElemWidth);
                if (img[at] == 255 || img[at] == 0 && delta < 0) img[at] = 7; // an escape slot or a zero: any other value is wrong as well
                else img[at] = (byte)(img[at] + delta);
            }));
        }
        OracleSection ovf = SesaOracle.FromBytes(SesaTestFile.Image(c)).Get(SesaFormat.KindLcpOverflow);
        for (int k = 0; k < Math.Min(8, (int)ovf.Count); k++)
        {
            int pair = k;
            cases.Add(Damaged($"overflow value {pair} + 1", c, (img, o) =>
            {
                OracleSection s = o.Get(SesaFormat.KindLcpOverflow);
                WriteUInt(img, s, (2 * pair) + 1, BinaryPrimitives.ReadUInt32LittleEndian(img.AsSpan((int)s.Offset + (8 * pair) + 4)) + 1);
            }));
            cases.Add(Damaged($"overflow rank {pair} + 1", c, (img, o) =>
            {
                OracleSection s = o.Get(SesaFormat.KindLcpOverflow);
                WriteUInt(img, s, 2 * pair, BinaryPrimitives.ReadUInt32LittleEndian(img.AsSpan((int)s.Offset + (8 * pair))) + 1);
            }));
        }
        foreach ((string name, byte[] image) in cases) AssertOnlyStructureSees(name, image, headerMaySee: true);
    }

    [TestCase("banana")]
    [TestCase("abababababab")]
    [TestCase("GATTACAGATTACAGGTA")]
    public void TheChildTable_WithAnyOneCellChanged_IsCaught(string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        for (int cell = 0; cell <= c.N; cell++)
        {
            int k = cell;
            foreach (uint xor in new uint[] { 1, 4, 0x80000000 })
            {
                (string name, byte[] image) = Damaged($"CHILD[{k}] xor {xor:X} of '{text}'", c, (img, o) =>
                {
                    OracleSection s = o.Get(SesaFormat.KindChild);
                    WriteUInt(img, s, k, BinaryPrimitives.ReadUInt32LittleEndian(img.AsSpan((int)s.Offset + (4 * k))) ^ xor);
                });
                AssertOnlyStructureSees(name, image);
            }
        }
    }

    [Test]
    public void ThePsiArray_WithAnEntryOrThePrimaryRowChanged_IsCaught()
    {
        LcpChildCase c = LcpChildCase.Of("GATTACAGATTACAGGTA");
        var cases = new List<(string, byte[])>();
        for (int i = 0; i <= c.N; i++)
        {
            int at = i;
            cases.Add(Damaged($"PSI[{at}] + 1", c, (img, o) =>
            {
                OracleSection s = o.Get(SesaFormat.KindPsi);
                WriteUInt(img, s, at, (BinaryPrimitives.ReadUInt32LittleEndian(img.AsSpan((int)s.Offset + (4 * at))) + 1) % (uint)(c.N + 1));
            }));
        }
        cases.Add(Damaged("PSI PARAM0 + 1", c, (img, o) => BinaryPrimitives.WriteInt64LittleEndian(img.AsSpan(128 + (64 * o.Get(SesaFormat.KindPsi).Index) + 40), o.Get(SesaFormat.KindPsi).Param0 + 1)));
        foreach ((string name, byte[] image) in cases) AssertOnlyStructureSees(name, image);
    }

    [Test]
    public void TheText_WithASymbolChanged_OrOneAboveSigma_IsCaught()
    {
        LcpChildCase c = LcpChildCase.Of("GATTACAGATTACAGGTA"); // sigma 4: 2 bits
        var cases = new List<(string, byte[])>();
        for (int i = 0; i < c.N; i++)
        {
            int at = i;
            cases.Add(Damaged($"symbol {at} changed", c, (img, o) =>
            {
                OracleSection s = o.Get(SesaFormat.KindText);
                int bit = at * 2;
                img[s.Offset + (bit / 8)] ^= (byte)(0x80 >> (bit % 8)); // flips the high bit of the symbol: another symbol of the same width
            }));
        }
        LcpChildCase three = LcpChildCase.Of("abcabcabcabcab"); // sigma 3: the 2-bit value 3 is not a symbol
        cases.Add(Damaged("symbol 3 of sigma 3", three, (img, o) =>
        {
            OracleSection s = o.Get(SesaFormat.KindText);
            img[s.Offset] |= 0xC0;
        }));
        foreach ((string name, byte[] image) in cases) AssertOnlyStructureSees(name, image);
    }

    [Test]
    public void TheCounts_WithTwoEntriesMovedByOne_AreCaught_TheSumStillMatches()
    {
        LcpChildCase c = LcpChildCase.Of("GATTACAGATTACAGGTA");
        (string name, byte[] image) = Damaged("counts moved", c, (img, o) =>
        {
            OracleSection s = o.Get(SesaFormat.KindCounts);
            ulong a = BinaryPrimitives.ReadUInt64LittleEndian(img.AsSpan((int)s.Offset));
            ulong b = BinaryPrimitives.ReadUInt64LittleEndian(img.AsSpan((int)s.Offset + 8));
            BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)s.Offset), a + 1);
            BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)s.Offset + 8), b - 1);
        });
        AssertOnlyStructureSees(name, image);
    }

    [Test]
    public void TheCounts_ThatDoNotSumToN_AreRefusedAtHeader()
    {
        LcpChildCase c = LcpChildCase.Of("GATTACAGATTACAGGTA");
        (_, byte[] image) = Damaged("counts + 1", c, (img, o) =>
        {
            OracleSection s = o.Get(SesaFormat.KindCounts);
            BinaryPrimitives.WriteUInt64LittleEndian(img.AsSpan((int)s.Offset), BinaryPrimitives.ReadUInt64LittleEndian(img.AsSpan((int)s.Offset)) + 1);
        });
        using var dir = new TempDir();
        Assert.That(() => SesaTestFile.OpenImage(dir, image).Dispose(), Throws.TypeOf<InvalidDataException>().With.Message.Contains("COUNTS"));
    }

    [Test]
    public void TheHeaderStatistics_NodeCount_LrsDepth_LrsRank_AreRecomputedAtStructure()
    {
        LcpChildCase c = LcpChildCase.Of("GATTACAGATTACAGGTAGATTACA");
        var cases = new List<(string, byte[])>();
        foreach ((string label, int offset, long delta) in new[] { ("NODE_COUNT + 1", 24, 1L), ("LRS_DEPTH - 1", 56, -1L), ("LRS_RANK + 1", 64, 1L) })
        {
            byte[] image = SesaTestFile.Image(c);
            BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(offset), BinaryPrimitives.ReadInt64LittleEndian(image.AsSpan(offset)) + delta);
            SesaTestFile.ResealHeader(image);
            cases.Add((label, image));
        }
        foreach ((string name, byte[] image) in cases) AssertOnlyStructureSees(name, image, headerMaySee: true);
    }

    [Test]
    public void TheTextHash_IsCheckedOnlyWhenAsked()
    {
        LcpChildCase c = LcpChildCase.Of("GATTACAGATTACAGGTA");
        byte[] image = SesaTestFile.Image(c);
        image[80] ^= 0xFF;
        SesaTestFile.ResealHeader(image);
        using var dir = new TempDir();
        using (EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, sha: false)) Assert.That(tree, Is.Not.Null);
        Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, sha: true, name: "h.sesa").Dispose(), Throws.TypeOf<InvalidDataException>().With.Message.Contains("SHA-256"));
        using EsaSuffixTree good = SesaTestFile.OpenImage(dir, SesaTestFile.Image(c), EsaValidationLevel.Structure, sha: true, name: "g.sesa");
        Assert.That(good, Is.Not.Null);
    }

    [Test]
    public void AHeaderThatRecordsNoHash_IsNotRefusedWhenTheHashIsAsked()
    {
        LcpChildCase c = LcpChildCase.Of("GATTACA");
        byte[] image = SesaTestFile.Image(c);
        image.AsSpan(80, 32).Clear();
        SesaTestFile.ResealHeader(image);
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, sha: true);
        Assert.That(tree, Is.Not.Null);
    }

    // ---------------------------------------------------------------- versions, kinds, codecs

    private static void Patch(byte[] image, Action<byte[]> edit)
    {
        edit(image);
        SesaTestFile.ResealHeader(image);
    }

    [Test]
    public void ALargerMajorVersion_IsNotSupported_ASmallerOneIsCorrupt_AMinorOneIsOpened()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        using var dir = new TempDir();
        byte[] major2 = SesaTestFile.Image(c);
        Patch(major2, i => BinaryPrimitives.WriteUInt16LittleEndian(i.AsSpan(8), 2));
        Assert.That(() => SesaTestFile.OpenImage(dir, major2).Dispose(), Throws.TypeOf<NotSupportedException>().With.Message.Contains("2").And.Message.Contains("1"));
        byte[] major0 = SesaTestFile.Image(c);
        Patch(major0, i => BinaryPrimitives.WriteUInt16LittleEndian(i.AsSpan(8), 0));
        Assert.That(() => SesaTestFile.OpenImage(dir, major0, name: "z.sesa").Dispose(), Throws.TypeOf<InvalidDataException>());
        byte[] minor7 = SesaTestFile.Image(c);
        Patch(minor7, i => BinaryPrimitives.WriteUInt16LittleEndian(i.AsSpan(10), 7));
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, minor7, EsaValidationLevel.Structure, name: "m.sesa");
        Assert.That(tree.CountOccurrences("an"), Is.EqualTo(2));
    }

    [Test]
    public void ALargerMajor_WithABrokenCrc_IsStillNotSupportedFirst()
    {
        byte[] image = SesaTestFile.Image(LcpChildCase.Of("banana"));
        BinaryPrimitives.WriteUInt16LittleEndian(image.AsSpan(8), 9);
        using var dir = new TempDir();
        Assert.That(() => SesaTestFile.OpenImage(dir, image).Dispose(), Throws.TypeOf<NotSupportedException>(), "the version is read before the checksum so that a file of a later format is not called corrupt");
    }

    [Test]
    public void AnUnknownOptionalSection_IsIgnored_AnUnknownRequiredOne_IsRefused()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        using var dir = new TempDir();
        foreach ((int kind, int flags, bool refused) in new[] { (7, 2, false), (8, 2, false), (9, 2, false), (10, 2, false), (14, 2, false), (0x10000000, 2, false), (99, 3, true), (7, 3, true), (0x10000001, 3, true) })
        {
            byte[] image = WithExtraSection(SesaTestFile.Image(c), kind, flags);
            string name = $"k{kind}f{flags}.sesa";
            if (refused) Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, name: name).Dispose(), Throws.TypeOf<NotSupportedException>(), $"kind {kind} flags {flags}");
            else
            {
                using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, name: name);
                Assert.That(tree.CountOccurrences("an"), Is.EqualTo(2), $"kind {kind} flags {flags}");
            }
        }
    }

    /// <summary>The image with one more section entry whose payload is 100 bytes appended at a 64-byte boundary (FILE_SIZE and the CRCs follow).</summary>
    private static byte[] WithExtraSection(byte[] image, int kind, int flags)
    {
        int count = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(48));
        long offset = (image.Length + 63) / 64 * 64;
        var bigger = new byte[offset + 100];
        image.CopyTo(bigger, 0);
        for (int i = 0; i < 100; i++) bigger[offset + i] = (byte)(i + 1);
        Span<byte> e = bigger.AsSpan(128 + (64 * count), 64);
        BinaryPrimitives.WriteInt32LittleEndian(e, kind);
        BinaryPrimitives.WriteInt32LittleEndian(e[4..], flags);
        BinaryPrimitives.WriteInt64LittleEndian(e[8..], offset);
        BinaryPrimitives.WriteInt64LittleEndian(e[16..], 100);
        BinaryPrimitives.WriteInt64LittleEndian(e[24..], 100);
        BinaryPrimitives.WriteInt32LittleEndian(e[32..], 1);
        BinaryPrimitives.WriteInt32LittleEndian(bigger.AsSpan(48), count + 1);
        BinaryPrimitives.WriteInt64LittleEndian(bigger.AsSpan(32), bigger.Length);
        SesaTestFile.Reseal(bigger);
        return bigger;
    }

    [Test]
    public void AnUnknownCodecOfARequiredSection_AndAPositionWidthOfEight_AreNotSupported()
    {
        LcpChildCase c = LcpChildCase.Of("banana");
        using var dir = new TempDir();
        foreach ((string label, int kind, int codec) in new[] { ("TEXT codec 1", SesaFormat.KindText, 1), ("TEXT codec 3", SesaFormat.KindText, 3), ("CHILD tagged", SesaFormat.KindChild, 2), ("LCP codec 4", SesaFormat.KindLcp, 4), ("SA codec 2", SesaFormat.KindSa, 2) })
        {
            byte[] image = SesaTestFile.Image(c);
            Patch(image, i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(128 + (64 * SesaOracle.FromBytes(i).Get(kind).Index) + 36), codec));
            Assert.That(() => SesaTestFile.OpenImage(dir, image, name: "c.sesa").Dispose(), Throws.TypeOf<NotSupportedException>(), label);
        }
        byte[] wide = SesaTestFile.Image(c);
        Patch(wide, i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(72), 8));
        Assert.That(() => SesaTestFile.OpenImage(dir, wide, name: "w.sesa").Dispose(), Throws.TypeOf<NotSupportedException>().With.Message.Contains("POS_WIDTH"));
    }

    [Test]
    public void TheTableRules_OfHeaderLevel_EachRefusal()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.Random(120, 3, 3));
        using var dir = new TempDir();
        SesaOracle o = SesaOracle.FromBytes(SesaTestFile.Image(c));
        int Entry(int kind) => 128 + (64 * o.Get(kind).Index);
        var edits = new List<(string Label, Action<byte[]> Edit)>
        {
            ("section outside the file", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(3) + 16), 1L << 40)),
            ("negative offset", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(3) + 8), -64)),
            ("offset inside the header page", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(3) + 8), 192)),
            ("misaligned offset", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(3) + 8), o.Get(3).Offset + 4)),
            ("overlap", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(6) + 8), o.Get(3).Offset)),
            ("COUNT * ELEM_WIDTH != LENGTH", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(3) + 24), o.Get(3).Count + 1)),
            ("SA COUNT != n + 1", i => { BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(3) + 24), o.Get(3).Count - 1); BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(3) + 16), o.Get(3).Length - 4); }),
            ("ELEM_WIDTH 0 of the SA", i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(Entry(3) + 32), 0)),
            ("the SA is not required", i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(Entry(3) + 4), 2)),
            ("a required section is missing", i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(Entry(6)), 77)),
            ("a kind twice", i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(Entry(6)), 3)),
            ("TEXT bits too small for sigma", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(1) + 40), 1)),
            ("TEXT LENGTH wrong", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(1) + 16), o.Get(1).Length + 1)),
            ("TEXT_LEN differs from the SA", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(16), o.TextLen + 1)),
            ("SIGMA differs from ALPHABET", i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(76), o.Sigma + 1)),
            ("SECTION_TABLE_OFFSET", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(40), 192)),
            ("SECTION_ENTRY_SIZE", i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(52), 128)),
            ("SECTION_COUNT 63", i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(48), 63)),
            ("SECTION_COUNT 0", i => BinaryPrimitives.WriteInt32LittleEndian(i.AsSpan(48), 0)),
            ("ENDIAN_PROBE", i => BinaryPrimitives.WriteUInt32LittleEndian(i.AsSpan(12), 0x04030201)),
            ("FILE_SIZE", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(32), o.File.Length + 64)),
            ("reserved header word", i => BinaryPrimitives.WriteUInt32LittleEndian(i.AsSpan(120), 1)),
            ("reserved word of an entry", i => BinaryPrimitives.WriteUInt32LittleEndian(i.AsSpan(Entry(3) + 60), 1)),
            ("BUILD_FLAGS without the sentinel", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(112), 0)),
            ("a byte after the table is not zero", i => i[128 + (64 * o.SectionCount) + 10] = 1),
            ("ALPHABET not ascending", i => BinaryPrimitives.WriteUInt16LittleEndian(i.AsSpan((int)o.Get(2).Offset), 0xFF)),
            ("ALPHABET letter above 255", i => BinaryPrimitives.WriteUInt16LittleEndian(i.AsSpan((int)o.Get(2).Offset + 4), 0x100)),
            ("the TEXT padding", i => i[o.Get(1).Offset + o.Get(1).Length - 3] = 1),
            ("NODE_COUNT out of range", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(24), 3L)),
            ("LRS_DEPTH 0 with a rank", i => { BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(56), 0); BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(64), 3); }),
            ("PSI PARAM0 out of range", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(12) + 40), o.TextLen + 1)),
            ("LCP PARAM0 without an overflow list", i => BinaryPrimitives.WriteInt64LittleEndian(i.AsSpan(Entry(4) + 40), 5)),
        };
        foreach ((string label, Action<byte[]> edit) in edits)
        {
            byte[] image = SesaTestFile.Image(c);
            Patch(image, edit);
            string name = "t.sesa";
            AssertRefused(() => SesaTestFile.OpenImage(dir, image, name: name).Dispose(), label);
        }
    }

    [Test]
    public void AnEverythingElseTreeLike_TheOptionalSectionsMayBeAbsent()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        using var dir = new TempDir();
        foreach ((bool psi, bool counts, bool info) in new[] { (false, true, true), (true, false, true), (true, true, false), (false, false, false) })
        {
            using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, SesaTestFile.Image(c, psi, counts, info), EsaValidationLevel.Structure, sha: true, name: $"{psi}{counts}{info}.sesa");
            Assert.That(tree.CountOccurrences("ssi"), Is.EqualTo(2), $"psi {psi} counts {counts} info {info}");
        }
    }

    [Test]
    public void TheGapsBetweenSections_MustBeZero_AtSections()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.Random(150, 4, 2));
        byte[] image = SesaTestFile.Image(c);
        SesaOracle o = SesaOracle.FromBytes(image);
        var covered = new bool[image.Length];
        foreach (OracleSection s in o.Sections)
            for (long b = s.Offset; b < s.Offset + s.Length; b++) covered[b] = true;
        int tested = 0;
        using var dir = new TempDir();
        for (int b = SesaFormat.PageBytes; b < image.Length; b++)
        {
            if (covered[b]) continue;
            byte[] copy = (byte[])image.Clone();
            copy[b] = 0x5A;
            int at = b;
            AssertRefused(() => SesaTestFile.OpenImage(dir, copy, EsaValidationLevel.Sections, name: "g.sesa").Dispose(), $"gap byte {at}");
            tested++;
        }
        Assert.That(tested, Is.GreaterThan(100), "the container has gaps to test");
    }
}
