using System.Runtime.InteropServices;
using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.1, deliverable C (the default set): the 4 -> 5 byte re-encoder <see cref="Sesa5Encoder"/> against the layout of ADR-005 section 6.1. A container of POS_WIDTH 4 is built by the product's test writer
/// (<see cref="SesaTestFile.Image"/>), re-encoded, parsed back by <see cref="Sesa5Image"/> (which knows section 6.1 and nothing of the encoder) and compared with the logical arrays of the in-memory index of S1: suffix array,
/// LCP with its overflow list, the packed child table (a tagged input is converted), PSI, ISA and BUCKET, TEXT, ALPHABET, COUNTS and the header. The S1 reader refuses the output at <c>POS_WIDTH</c> with its documented message.
/// </summary>
[TestFixture]
public class Sesa5EncoderTests
{
    private static readonly string[] SmallTexts =
    [
        string.Empty, "a", "aa", "ab", "ba", "aba", "banana", "mississippi", "abracadabra", "aaaaaaaa", "abababab", "abaababaabaab", "xyzxyz",
    ];

    private static List<(string Name, string Text)> Texts()
    {
        var texts = new List<(string Name, string Text)>();
        foreach (string t in SmallTexts) texts.Add((t.Length == 0 ? "empty" : t, t));
        texts.Add(("dna-300", StringFamilies.RandomDna(300, 11)));
        texts.Add(("dna-3000", StringFamilies.RandomDna(3000, 12)));
        texts.Add(("protein-500", StringFamilies.RandomProtein(500, 13)));
        texts.Add(("fibonacci-400", StringFamilies.Fibonacci(400)));
        texts.Add(("thue-morse-300", StringFamilies.ThueMorse(300)));
        texts.Add(("runs-500", StringFamilies.Runs(500, 7, 14)));
        texts.Add(("periodic-5-400", StringFamilies.Periodic(400, 5, 15)));
        texts.Add(("a^1000", StringFamilies.AllSame(1000)));
        texts.Add(("(a^6 b)* 2000", StringFamilies.AkB(2000, 6)));
        texts.Add(("genome-like-2000", GenomeLike(2000, 16)));
        return texts;
    }

    public static IEnumerable<TestCaseData> Containers()
    {
        foreach ((string name, string text) in Texts())
            foreach (bool psi in new[] { true, false })
                foreach (bool tagged in new[] { false, true })
                    yield return new TestCaseData(text, psi, tagged).SetArgDisplayNames(name, $"psi={psi}", tagged ? "Tagged" : "Packed");
    }

    [TestCaseSource(nameof(Containers))]
    public void Encode_GivesTheSameLogicalIndex_InFiveByteEntries(string text, bool psi, bool tagged)
    {
        ChildCodec codec = tagged ? ChildCodec.Tagged : ChildCodec.Packed;
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] image4 = SesaTestFile.Image(c, psi: psi, codec: codec);
        byte[] image5 = Sesa5Encoder.Encode(image4);
        SameLogicalContent(c, image4, Sesa5Image.Parse(image5), psi);
    }

    public static IEnumerable<TestCaseData> FourAlphabets()
    {
        yield return new TestCaseData("dna", StringFamilies.RandomDna(100_000, 21), true);
        yield return new TestCaseData("genome-like", GenomeLike(100_000, 22), true);
        yield return new TestCaseData("protein", StringFamilies.RandomProtein(100_000, 23), true);
        yield return new TestCaseData("repeat", RepeatText(100_000, 24), false);
        yield return new TestCaseData("dna-1m", StringFamilies.RandomDna(1_000_000, 25), false);
        yield return new TestCaseData("repeat-1m", RepeatText(1_000_000, 26), true);
    }

    [TestCaseSource(nameof(FourAlphabets))]
    public void Encode_RoundTripsTheFourAlphabets_UpToOneMillionSymbols(string name, string text, bool tagged)
    {
        ChildCodec codec = tagged ? ChildCodec.Tagged : ChildCodec.Packed;
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] image4 = SesaTestFile.Image(c, psi: true, codec: codec);
        byte[] image5 = Sesa5Encoder.Encode(image4);
        Assert.That(image5.Length, Is.GreaterThan(image4.Length), $"{name}: 5-byte entries are larger");
        SameLogicalContent(c, image4, Sesa5Image.Parse(image5), psi: true);
    }

    [Test]
    public void TheEntriesAreFiveBytesLittleEndian_AndTheLayoutIsTheOneOfSection61()
    {
        LcpChildCase c = LcpChildCase.Of("mississippi");
        byte[] image4 = SesaTestFile.Image(c, psi: true, codec: ChildCodec.Packed);
        byte[] image5 = Sesa5Encoder.Encode(image4);
        SesaOracle o4 = SesaOracle.FromBytes(image4), o5 = SesaOracle.FromBytes(image5);
        Assert.That((o5.PosWidth, o5.VersionMajor, o5.VersionMinor), Is.EqualTo((5, 1, 1)));
        int[] sa = c.Sa;
        OracleSection s5 = o5.Get(SesaFormat.KindSa);
        Assert.That((s5.Count, s5.ElemWidth, s5.Length, s5.Codec), Is.EqualTo((12L, 5, 12L * 5 + 3, 1)));
        for (int i = 0; i < sa.Length; i++)
            for (int j = 0; j < 5; j++)
                Assert.That(image5[s5.Offset + (5 * i) + j], Is.EqualTo((byte)((long)sa[i] >> (8 * j))), $"SA[{i}] byte {j}");
        for (int j = 0; j < 3; j++) Assert.That(image5[s5.Offset + (5 * 12) + j], Is.Zero, $"padding byte {j}");
        // each section of the container, one entry of the table each, in file order, aligned, with the CRC of its payload
        Assert.That(o5.Sections.Select(s => s.Kind), Is.EqualTo(o4.Sections.OrderBy(s => s.Offset).Select(s => s.Kind)), "the sections keep the file order");
        Assert.That(o5.Sections.Select(s => s.Offset), Is.Ordered.Ascending);
        foreach (OracleSection s in o5.Sections) Assert.That(s.Crc, Is.EqualTo(SesaOracle.Crc32cBitwise(o5.Payload(s))), $"kind {s.Kind}: CRC-32C");
        Assert.That(o5.HeaderCrc, Is.EqualTo(o5.ExpectedHeaderCrc()));
        Assert.That(o5.FileSize, Is.EqualTo(image5.Length));
    }

    [Test]
    public void PackPositions_KnownAnswers_AboveTwoToThe32()
    {
        ulong[] values = [0x12_3456_789A, 0, 0xFF_FFFF_FFFF, 1UL << 32, 0x01_0203_0405, 1, 255, 256];
        byte[] expected =
        [
            0x9A, 0x78, 0x56, 0x34, 0x12,
            0, 0, 0, 0, 0,
            0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
            0, 0, 0, 0, 1,
            0x05, 0x04, 0x03, 0x02, 0x01,
            1, 0, 0, 0, 0,
            255, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 0,
        ];
        byte[] actual = Sesa5Encoder.PackPositions(values);
        Assert.That(actual, Is.EqualTo(expected));
        Assert.That(Sesa5Encoder.PackPositions([]), Is.EqualTo(new byte[3]), "an empty section is its 3 bytes of padding");
        Assert.Throws<ArgumentOutOfRangeException>(() => Sesa5Encoder.PackPositions([1UL << 40]));
        // and the decoder reads them back by one 8-byte load and a mask, whatever the neighbours hold
        var file = new byte[4096 + actual.Length];
        actual.CopyTo(file, 4096);
        for (int i = 0; i < values.Length; i++) Assert.That(Sesa5Image.Pos40At(file, 4096, i), Is.EqualTo(values[i]), $"entry {i}");
    }

    public static IEnumerable<TestCaseData> EveryText() => Texts().Select(t => new TestCaseData(t.Text).SetArgDisplayNames(t.Name));

    [TestCaseSource(nameof(EveryText))]
    public void TheTaggedChildTable_BecomesThePackedTableOfTheInMemoryIndex(string text)
    {
        LcpChildCase c = LcpChildCase.Of(text);
        Assert.That(Sesa5Encoder.TaggedToPacked(c.Tagged), Is.EqualTo(c.Packed));
    }

    [Test]
    public void TheS1Reader_RefusesTheEncodedContainer_WithTheDocumentedDamageMessage()
    {
        // ADR-005 section 6: "SesaHeaderReader.cs:96-97: any other width but 4 is InvalidDataException("POS_WIDTH is 5: 4 or 8")". S2.7 replaces this test by the golden refusal of the frozen 1.0 reader.
        using var dir = new TempDir();
        byte[] image5 = Sesa5Encoder.Encode(SesaTestFile.Image(LcpChildCase.Of("mississippi")));
        string path = dir.File("w5.sesa");
        File.WriteAllBytes(path, image5);
        foreach (EsaValidationLevel level in new[] { EsaValidationLevel.Header, EsaValidationLevel.Sections, EsaValidationLevel.Structure })
        {
            InvalidDataException e = Assert.Throws<InvalidDataException>(() => SesaFile.Open(path, new EsaOpenOptions { Validation = level }).Dispose(), $"{level}")!;
            Assert.That(e.Message, Is.EqualTo("POS_WIDTH is 5: 4 or 8"), $"{level}: the refusal is at the width, so the header CRC, FILE_SIZE and the table of the output are valid");
        }
        Assert.That(Assert.Throws<InvalidDataException>(() => EsaSuffixTreeFactory.Open(path).Dispose())!.Message, Is.EqualTo("POS_WIDTH is 5: 4 or 8"));
        Assert.That(Assert.Throws<InvalidDataException>(() => SuffixIndexFactory.Open(path))!.Message, Is.EqualTo("POS_WIDTH is 5: 4 or 8"));
    }

    [Test]
    public void TheEncoder_RefusesWhatIsNotAContainerOfWidth4()
    {
        byte[] good = SesaTestFile.Image(LcpChildCase.Of("banana"));
        Assert.Throws<InvalidDataException>(() => Sesa5Encoder.Encode(new byte[100]), "shorter than a page");
        Assert.Throws<InvalidDataException>(() => Sesa5Encoder.Encode(new byte[8192]), "no magic");
        byte[] twice = Sesa5Encoder.Encode(good);
        Assert.Throws<InvalidDataException>(() => Sesa5Encoder.Encode(twice), "a container of width 5 is not re-encoded");
        byte[] widthFive = (byte[])good.Clone();
        widthFive[72] = 5;
        SesaTestFile.ResealHeader(widthFive);
        Assert.That(Assert.Throws<InvalidDataException>(() => Sesa5Encoder.Encode(widthFive))!.Message, Does.Contain("POS_WIDTH"), "a version 1.0 container that says POS_WIDTH 5 is refused for its width");
        byte[] badCrc = (byte[])good.Clone();
        badCrc[124] ^= 1;
        Assert.Throws<InvalidDataException>(() => Sesa5Encoder.Encode(badCrc), "the header CRC");
        byte[] badSize = (byte[])good.Clone();
        Array.Resize(ref badSize, good.Length + 64);
        Assert.Throws<InvalidDataException>(() => Sesa5Encoder.Encode(badSize), "FILE_SIZE");
        byte[] minor1 = (byte[])good.Clone();
        minor1[10] = 1;
        SesaTestFile.ResealHeader(minor1);
        Assert.Throws<InvalidDataException>(() => Sesa5Encoder.Encode(minor1), "version 1.1");
    }

    [Test]
    public void AnIsaAndABucketSection_AreEncodedLikeTheOtherPositionSections()
    {
        // the S1 builder writes neither: they are added to a 1.0 container by hand (the S1 reader ignores optional kinds 7 and 8)
        LcpChildCase c = LcpChildCase.Of("abracadabra");
        byte[] image4 = SesaTestFile.Image(c, psi: false);
        int[] sa = c.Sa;
        var isa = new uint[sa.Length];
        for (int r = 0; r < sa.Length; r++) isa[sa[r]] = (uint)r;
        uint[] bucket = [0, 1, 6, 8, 9, 11, 12];
        byte[] withExtras = AddSections(image4,
            (Sesa5Image.KindIsa, MemoryMarshal.AsBytes(isa.AsSpan()).ToArray(), sa.Length, 4, 1, 0L, 0L),
            (Sesa5Image.KindBucket, MemoryMarshal.AsBytes(bucket.AsSpan()).ToArray(), bucket.Length, 4, 1, 1L, 5L));
        using (var dir = new TempDir())
        {
            string path = dir.File("extras.sesa");
            File.WriteAllBytes(path, withExtras);
            using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections });
            Assert.That(file.Sections.Select(s => s.Kind), Does.Contain(Sesa5Image.KindIsa).And.Contain(Sesa5Image.KindBucket), "the S1 reader accepts the two optional kinds");
        }
        Sesa5Image i5 = Sesa5Image.Parse(Sesa5Encoder.Encode(withExtras));
        Assert.That(i5.Positions(Sesa5Image.KindIsa), Is.EqualTo(isa.Select(x => (ulong)x)));
        Assert.That(i5.Positions(Sesa5Image.KindBucket), Is.EqualTo(bucket.Select(x => (ulong)x)));
        OracleSection b = i5.Get(Sesa5Image.KindBucket);
        Assert.That((b.Count, b.ElemWidth, b.Length, b.Param0, b.Param1), Is.EqualTo(((long)bucket.Length, 5, (5L * bucket.Length) + 3, 1L, 5L)));
        OracleSection a = i5.Get(Sesa5Image.KindIsa);
        Assert.That((a.Count, a.ElemWidth, a.Length), Is.EqualTo(((long)sa.Length, 5, (5L * sa.Length) + 3)));
    }

    [Test]
    public void AnOverflowListOfAnLcp_IsEncodedAsPairsOfTenBytes()
    {
        // a text whose LCP has escapes in U8: a 400-symbol block of random DNA that occurs twice in random DNA (about 150 LCP values of 255 or more, the rest of the 10,000 below)
        string block = StringFamilies.RandomDna(400, 6);
        string text = StringFamilies.RandomDna(4000, 7) + block + StringFamilies.RandomDna(4000, 8) + block + StringFamilies.RandomDna(1000, 9);
        LcpChildCase c = LcpChildCase.Of(text);
        Assert.That(c.Lcp.OverflowCount, Is.GreaterThan(0), "the text must have escapes (otherwise the test checks nothing)");
        byte[] image4 = SesaTestFile.Image(c);
        byte[] image5 = Sesa5Encoder.Encode(image4);
        Sesa5Image i5 = Sesa5Image.Parse(image5);
        OracleSection o = i5.Get(SesaFormat.KindLcpOverflow);
        Assert.That((o.Count, o.ElemWidth, o.Length), Is.EqualTo((c.Lcp.OverflowCount, 10, (10 * c.Lcp.OverflowCount) + 3)));
        SameLogicalContent(c, image4, i5, psi: true);
    }

    [Test]
    public void TheSparseWriterAtWidth5_AndTheEncoder_MakeTheSameLogicalContainer()
    {
        // two independent implementations of section 6.1 (the writer of the holes and the re-encoder), one reader that knows only the text of the section
        foreach (ClosedFormText f in new ClosedFormText[] { new PowerOfA(150), new PeriodicAkB(199, 4), PeriodicAkB.AbRepeat(80) })
        {
            using var dir = new TempDir();
            string text = f.Materialize();
            LcpChildCase c = LcpChildCase.Of(text);
            SparseSpec spec = SparseSpec.For(f, 5, psi: true);
            string path = dir.File("w5.sesa");
            SparseContainer.Write(path, spec, SparseContainer.AllPatches(f, SparseContainer.Plan(spec)));
            Sesa5Image written = Sesa5Image.Parse(File.ReadAllBytes(path));
            Sesa5Image encoded = Sesa5Image.Parse(Sesa5Encoder.Encode(SesaTestFile.Image(c, psi: true)));
            foreach (int kind in new[] { SesaFormat.KindSa, SesaFormat.KindChild, SesaFormat.KindPsi })
                Assert.That(written.Positions(kind), Is.EqualTo(encoded.Positions(kind)), $"{f}: kind {kind}");
            Assert.That(written.Lcp(), Is.EqualTo(encoded.Lcp()), $"{f}: LCP");
            Assert.That(written.Oracle.NodeCount, Is.EqualTo(encoded.Oracle.NodeCount), $"{f}: NODE_COUNT");
            Assert.That(written.Oracle.LrsRank, Is.EqualTo(encoded.Oracle.LrsRank), $"{f}: LRS_RANK");
        }
    }

    [Test]
    public void TheDecoder_RefusesWhatSection61Forbids()
    {
        byte[] good = Sesa5Encoder.Encode(SesaTestFile.Image(LcpChildCase.Of("mississippi"), psi: true));
        Assert.DoesNotThrow(() => Sesa5Image.Parse(good));
        SesaOracle o = SesaOracle.FromBytes(good);
        OracleSection sa = o.Get(SesaFormat.KindSa), child = o.Get(SesaFormat.KindChild);

        byte[] Damage(Action<byte[]> action, bool reseal = true)
        {
            byte[] copy = (byte[])good.Clone();
            action(copy);
            if (reseal) SesaTestFile.Reseal(copy);
            return copy;
        }

        for (int j = 1; j <= 3; j++)
        {
            int back = j;
            Assert.Throws<InvalidDataException>(() => Sesa5Image.Parse(Damage(b => b[sa.Offset + sa.Length - back] = 1)), $"the padding is zero (byte {4 - j} of 3)");
        }
        Assert.Throws<InvalidDataException>(() => Sesa5Image.Parse(Damage(b => b[10] = 0)), "POS_WIDTH 5 needs minor 1");
        Assert.Throws<InvalidDataException>(() => Sesa5Image.Parse(Damage(b => b[72] = 4)), "POS_WIDTH 4");
        Assert.Throws<InvalidDataException>(() => Sesa5Image.Parse(Damage(b => b[128 + (64 * child.Index) + 36] = 2)), "CHILD codec 2 does not exist at width 5");
        Assert.Throws<InvalidDataException>(() => Sesa5Image.Parse(Damage(b => b[128 + (64 * sa.Index) + 32] = 4)), "ELEM_WIDTH 4 in a width-5 container");
        Assert.Throws<InvalidDataException>(() => Sesa5Image.Parse(Damage(b => b[child.Offset + 1] ^= 1, reseal: false)), "CRC-32C");
        Assert.Throws<InvalidDataException>(() => Sesa5Image.Parse(good.AsSpan(0, good.Length - 8).ToArray()), "FILE_SIZE");
        Assert.Throws<InvalidDataException>(() => Sesa5Image.Parse(Damage(b => b[124] ^= 1, reseal: false)), "HEADER_CRC32C");
    }

    // ------------------------------------------------------------------------------------------ helpers

    /// <summary>The container of POS_WIDTH 5 holds exactly the logical arrays of the in-memory index and the payloads that section 6.1 leaves unchanged.</summary>
    private static void SameLogicalContent(LcpChildCase c, byte[] image4, Sesa5Image i5, bool psi)
    {
        SesaOracle o4 = SesaOracle.FromBytes(image4);
        SesaOracle o5 = i5.Oracle;
        Assert.That((o5.TextLen, o5.NodeCount, o5.LrsDepth, o5.LrsRank, o5.Sigma, o5.BuildFlags), Is.EqualTo((o4.TextLen, o4.NodeCount, o4.LrsDepth, o4.LrsRank, o4.Sigma, o4.BuildFlags)), "header fields");
        Assert.That(o5.TextSha256, Is.EqualTo(o4.TextSha256), "TEXT_SHA256");
        Assert.That(o5.Sections.Select(s => s.Kind), Is.EqualTo(o4.Sections.OrderBy(s => s.Offset).Select(s => s.Kind)), "the sections");
        foreach (OracleSection s4 in o4.Sections)
        {
            OracleSection s5 = i5.Get(s4.Kind);
            Assert.That((s5.Flags, s5.Count, s5.Param0, s5.Param1), Is.EqualTo((s4.Flags, s4.Count, s4.Param0, s4.Param1)), $"kind {s4.Kind}: flags, COUNT, PARAM0, PARAM1");
        }

        Assert.That(i5.Positions(SesaFormat.KindSa), Is.EqualTo(c.Sa.Select(x => (ulong)x)), "SA");
        OracleSection child = i5.Get(SesaFormat.KindChild);
        Assert.That((child.Codec, child.ElemWidth), Is.EqualTo((SesaFormat.ChildCodecPacked, 5)), "CHILD is packed at width 5");
        Assert.That(i5.Positions(SesaFormat.KindChild), Is.EqualTo(c.Packed.Select(x => (ulong)x)), "CHILD (the packed table of the in-memory index)");
        if (psi)
        {
            Assert.That(i5.Positions(SesaFormat.KindPsi), Is.EqualTo(SesaRun.NaivePsi(c.Sa).Select(x => (ulong)x)), "PSI");
            Assert.That(i5.Get(SesaFormat.KindPsi).Param0, Is.EqualTo(SesaRun.PrimaryRow(c.Sa)), "PSI primary row");
        }
        else
        {
            Assert.That(i5.Find(SesaFormat.KindPsi), Is.Null, "no PSI");
        }

        Assert.That(i5.Lcp(), Is.EqualTo(c.Lcp.ToArray().Select(x => (long)x)), "LCP, escapes resolved");
        Assert.That(i5.Get(SesaFormat.KindLcp).Codec, Is.EqualTo(o4.Get(SesaFormat.KindLcp).Codec), "LCP codec");
        Assert.That(i5.Payload(i5.Get(SesaFormat.KindLcp)).ToArray(), Is.EqualTo(o4.Payload(SesaFormat.KindLcp)), "LCP base array: no position in it, copied");
        if (c.Lcp.OverflowCount > 0)
        {
            ReadOnlySpan<uint> pairs4 = MemoryMarshal.Cast<byte, uint>(o4.Payload(SesaFormat.KindLcpOverflow));
            (ulong Rank, ulong Value)[] pairs5 = i5.Overflow();
            Assert.That(pairs5, Has.Length.EqualTo(pairs4.Length / 2));
            for (int i = 0; i < pairs5.Length; i++) Assert.That(pairs5[i], Is.EqualTo(((ulong)pairs4[2 * i], (ulong)pairs4[(2 * i) + 1])), $"overflow pair {i}");
        }
        foreach (int kind in new[] { SesaFormat.KindText, SesaFormat.KindAlphabet, SesaFormat.KindCounts })
            Assert.That(i5.Payload(i5.Get(kind)).ToArray(), Is.EqualTo(o4.Payload(kind)), $"kind {kind} is copied");
        string info4 = Encoding.UTF8.GetString(o4.Payload(SesaFormat.KindBuildInfo)), info5 = Encoding.UTF8.GetString(i5.Payload(i5.Get(SesaFormat.KindBuildInfo)));
        Assert.That(info5, Is.EqualTo(info4.Replace("format=SESA 1.0", "format=SESA 1.1", StringComparison.Ordinal).Replace("child=layout A, tagged", "child=layout A, packed", StringComparison.Ordinal)), "BUILD_INFO");
        Assert.That(info5.Count(ch => ch == '\n'), Is.EqualTo(i5.Get(SesaFormat.KindBuildInfo).Count), "BUILD_INFO lines");
    }

    /// <summary>Appends sections (after the last one, aligned, with their CRC) to a 1.0 container and reseals it.</summary>
    private static byte[] AddSections(byte[] image, params (int Kind, byte[] Payload, long Count, int ElemWidth, int Codec, long Param0, long Param1)[] extra)
    {
        SesaOracle o = SesaOracle.FromBytes(image);
        var sections = o.Sections.Select(s => new SesaSection(s.Kind, s.Flags, s.Offset, s.Length, s.Count, s.ElemWidth, s.Codec, s.Param0, s.Param1, s.Crc)).ToList();
        long end = sections.Max(s => s.Offset + s.Length);
        var tail = new List<(long Offset, byte[] Payload)>();
        foreach (var e in extra)
        {
            long at = SesaFormat.AlignedOffset(end, e.Payload.Length);
            sections.Add(new SesaSection(e.Kind, SesaFormat.FlagHasCrc, at, e.Payload.Length, e.Count, e.ElemWidth, e.Codec, e.Param0, e.Param1, Crc32C.Compute(e.Payload)));
            tail.Add((at, e.Payload));
            end = at + e.Payload.Length;
        }
        var header = new SesaHeader
        {
            TextLength = o.TextLen,
            NodeCount = o.NodeCount,
            FileSize = end,
            LrsDepth = o.LrsDepth,
            LrsRank = o.LrsRank,
            PosWidth = 4,
            Sigma = o.Sigma,
            TextSha256 = o.TextSha256,
            BuildFlags = o.BuildFlags,
        };
        var file = new byte[end];
        Array.Copy(image, file, image.Length);
        SesaFormat.EncodeHeaderPage(header, sections, withMagic: true).CopyTo(file, 0);
        foreach ((long offset, byte[] payload) in tail) payload.CopyTo(file, offset);
        return file;
    }

    /// <summary>Random ACGT with runs of N, and repeats of earlier blocks with a few substitutions (five letters: 3 bits per symbol).</summary>
    internal static string GenomeLike(int n, ulong seed)
    {
        var rng = new SplitMix64(seed);
        var sb = new StringBuilder(n + 2000);
        const string acgt = "ACGT";
        while (sb.Length < n)
        {
            int kind = rng.NextInt(100);
            if (kind < 3)
            {
                sb.Append('N', 20 + rng.NextInt(400));
            }
            else if (kind < 12 && sb.Length > 2000)
            {
                int length = 50 + rng.NextInt(800), from = rng.NextInt(sb.Length - length);
                for (int i = 0; i < length; i++) sb.Append(rng.NextInt(100) == 0 ? acgt[rng.NextInt(4)] : sb[from + i]);
            }
            else
            {
                int length = 100 + rng.NextInt(900);
                for (int i = 0; i < length; i++) sb.Append(acgt[rng.NextInt(4)]);
            }
        }
        return sb.ToString(0, n);
    }

    /// <summary>A block of random DNA repeated, with 1 % substitutions.</summary>
    internal static string RepeatText(int n, ulong seed)
    {
        var rng = new SplitMix64(seed);
        string block = StringFamilies.RandomDna(Math.Min(5000, n), seed ^ 0xB10C);
        var sb = new StringBuilder(n);
        while (sb.Length < n)
            foreach (char ch in block)
            {
                if (sb.Length >= n) break;
                sb.Append(rng.NextInt(100) == 0 ? "ACGT"[rng.NextInt(4)] : ch);
            }
        return sb.ToString();
    }
}
