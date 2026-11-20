using System.Buffers.Binary;
using System.Runtime.InteropServices;
using System.Text;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.Facade;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// CHILD codec 2 of the container (work package 1b.9a, ADR-003 amendment of 2026-11-19): layout A with the kind of a field in bits 30-31 (the in-memory <c>TaggedChild</c>), written by the builder for a text of at most 2^30 - 1
/// symbols and read in place. The pattern descent of a container touched about 18 % more cache lines and 15 % more pages per query with codec 1 (the LCP array is read to tell the kind of a field), and codec 2 needs no such read.
/// Codec 1 stays readable and writable (<see cref="SesaBuildOptions.ChildCodec"/>).
/// </summary>
[TestFixture]
public class SesaChildCodecTests
{
    private static readonly string[] Named =
    [
        "", "a", "aa", "ab", "ba", "aba", "abab", "banana", "mississippi", "abracadabra", "aaaaaaaaaa", "abababababab", "abcabxabcd", "xabxac", "GATTACAGATTACA", "ACGTACGTTACGT",
        "a$b$", "the quick brown fox jumps over the lazy dog", "aabaabaaabaabaaab", "café ÿres \u0001\u0000x",
    ];

    private static void Compare(string text, ulong seed, string context, ChildCodec codec, EsaAccessMode access, EsaValidationLevel level, int maxPatternLength = 5)
    {
        using EsaSuffixTree tree = SesaTestFile.FileTree(text, psi: seed % 2 == 0, access: access, level: level, codec: codec);
        using var v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        FacadeChecks.CheckEverything(tree, UkkonenTree.Build(text), v6, text, seed, $"[{codec} {access} {level}] {context}", maxPatternLength);
        FacadeChecks.CheckAgainstTree(tree, memory, "esa-memory", text, $"[{codec} {access} {level}] {context}");
    }

    [Test]
    public void NamedTexts_OnATaggedContainer_EveryMember_EqualTheOracles_InBothModesAndAtTheStructureLevel()
    {
        ulong seed = 1;
        foreach (string text in Named)
        {
            Compare(text, seed++, "named", ChildCodec.Tagged, EsaAccessMode.Checked, EsaValidationLevel.Header);
            Compare(text, seed++, "named", ChildCodec.Tagged, EsaAccessMode.Trusted, EsaValidationLevel.Header);
            Compare(text, seed++, "named", ChildCodec.Tagged, EsaAccessMode.Checked, EsaValidationLevel.Structure);
        }
    }

    [Test]
    public void AdversarialFamilies_OnATaggedContainer_EveryMember()
    {
        foreach ((string name, string text) in StringFamilies.Adversarial(150, 3))
            if (SesaTestFile.IsByteText(text)) Compare(text, 7, name, ChildCodec.Tagged, EsaAccessMode.Checked, EsaValidationLevel.Structure, 4);
        foreach (int n in new[] { 1, 2, 3, 4, 5, 8, 13 })
            foreach ((string name, string text) in StringFamilies.Adversarial(n, 11))
                if (SesaTestFile.IsByteText(text)) Compare(text, 5, $"{name} n={n}", ChildCodec.Tagged, EsaAccessMode.Trusted, EsaValidationLevel.Header, 4);
    }

    [Test]
    public void TheTwoCodecs_GiveTheSameTree_AndTheSameTraversal()
    {
        var rng = new SplitMix64(0xC0DEC);
        foreach (int sigma in new[] { 2, 4, 5, 20 })
        {
            string text = StringFamilies.Random(1500, sigma, rng.Next());
            using EsaSuffixTree packed = SesaTestFile.FileTree(text, codec: ChildCodec.Packed);
            using EsaSuffixTree tagged = SesaTestFile.FileTree(text, codec: ChildCodec.Tagged);
            FacadeChecks.CheckAgainstTree(tagged, packed, "packed container", text, $"sigma {sigma}");
            Assert.That(tagged.File!.Parsed.Child.Codec, Is.EqualTo(SesaFormat.ChildCodecTagged));
            Assert.That(packed.File!.Parsed.Child.Codec, Is.EqualTo(SesaFormat.ChildCodecPacked));
            Assert.That(tagged.Index.Child.Kind, Is.EqualTo(ChildLayoutKind.Tagged), "the table of a codec-2 container is the tagged layout");
            Assert.That(packed.Index.Child.Kind, Is.EqualTo(ChildLayoutKind.Packed));
            Assert.That(tagged.Index.Child.Mapped.IsTagged, Is.True);
            Assert.That(packed.Index.Child.Mapped.IsTagged, Is.False);
        }
    }

    [Test]
    public void TheStructureLevel_OfATaggedTable_AcceptsTheValidTable_AndRefusesEveryDamagedCell()
    {
        string text = StringFamilies.Random(120, 4, 5);
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] clean = SesaTestFile.Image(c, codec: ChildCodec.Tagged);
        using var dir = new TempDir();
        using (EsaSuffixTree ok = SesaTestFile.OpenImage(dir, clean, EsaValidationLevel.Structure, sha: true)) Assert.That(ok.File!.Parsed.Child.Codec, Is.EqualTo(SesaFormat.ChildCodecTagged));
        OracleSection child = SesaOracle.FromBytes(clean).Get(SesaFormat.KindChild);
        int cells = (int)child.Count;
        int refused = 0;
        for (int cell = 0; cell < cells; cell++)
        {
            uint original = BinaryPrimitives.ReadUInt32LittleEndian(clean.AsSpan((int)child.Offset + (4 * cell)));
            foreach (uint damaged in new[] { original + 1, original ^ (1u << 30), original ^ (1u << 31), original == 0 ? 5u : 0u })
            {
                byte[] image = (byte[])clean.Clone();
                BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)child.Offset + (4 * cell)), damaged);
                SesaTestFile.Reseal(image);
                Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, name: "d.sesa").Dispose(), Throws.TypeOf<InvalidDataException>(), $"cell {cell}: {original:X8} -> {damaged:X8}");
                refused++;
            }
        }
        Assert.That(refused, Is.EqualTo(4 * cells));
    }

    [TestCase(false)]
    [TestCase(true)]
    public void ACellOutsideTheChildTable_InCheckedMode_IsInvalidData_ForTheNextLIndexOfBothLayouts(bool tagged)
    {
        ChildCodec codec = tagged ? ChildCodec.Tagged : ChildCodec.Packed;
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of("mississippimississippi"), codec: codec);
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, access: EsaAccessMode.Checked);
        MappedChild child = tree.File!.Child.Mapped;
        int cells = (int)(child.PayloadBytes / 4);
        Assert.That(child.IsTagged, Is.EqualTo(tagged));
        foreach (int k in new[] { cells, cells + 1, cells + 1000, int.MaxValue - 1, -1 })
            Assert.That(() => child.NextLIndex(k, int.MaxValue), Throws.TypeOf<InvalidDataException>(), $"{codec} cell {k}");
        Assert.That(() => child.NextLIndex(1, int.MaxValue), Throws.Nothing);
    }

    [Test]
    public void ATaggedTable_ThatTheLayoutOfCodec1Wrote_IsRefusedAtTheStructureLevel_AndTheOtherWayAround()
    {
        string text = "mississippimississippi";
        LcpChildCase c = LcpChildCase.Of(text);
        using var dir = new TempDir();
        foreach ((ChildCodec written, int claimed) in new[] { (ChildCodec.Packed, SesaFormat.ChildCodecTagged), (ChildCodec.Tagged, SesaFormat.ChildCodecPacked) })
        {
            byte[] image = SesaTestFile.Image(c, codec: written);
            OracleSection child = SesaOracle.FromBytes(image).Get(SesaFormat.KindChild);
            BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(128 + (64 * child.Index) + 36), claimed);
            SesaTestFile.ResealHeader(image);
            Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, name: "w.sesa").Dispose(), Throws.TypeOf<InvalidDataException>(), $"{written} cells under codec {claimed}");
        }
    }

    [TestCase(0)]
    [TestCase(3)]
    [TestCase(4)]
    [TestCase(255)]
    public void AnUnknownChildCodec_IsNotSupported_AtTheHeaderLevel(int codec)
    {
        byte[] image = SesaTestFile.Image(LcpChildCase.Of("banana"));
        OracleSection child = SesaOracle.FromBytes(image).Get(SesaFormat.KindChild);
        BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(128 + (64 * child.Index) + 36), codec);
        SesaTestFile.ResealHeader(image);
        using var dir = new TempDir();
        Assert.That(() => SesaTestFile.OpenImage(dir, image).Dispose(), Throws.TypeOf<NotSupportedException>().With.Message.Contains("CHILD"));
    }

    [Test]
    public void TheTaggedCodec_NeedsRanksThatFitThirtyBits_TheHeaderCheckSaysSo()
    {
        Assert.That(() => SesaHeaderReader.CheckChildCodec(SesaFormat.ChildCodecTagged, TaggedChild.MaxLength), Throws.Nothing);
        Assert.That(() => SesaHeaderReader.CheckChildCodec(SesaFormat.ChildCodecTagged, TaggedChild.MaxLength + 1L), Throws.TypeOf<InvalidDataException>().With.Message.Contains("30 bits"));
        Assert.That(() => SesaHeaderReader.CheckChildCodec(SesaFormat.ChildCodecPacked, TaggedChild.MaxLength + 1L), Throws.Nothing);
        Assert.That(() => SesaHeaderReader.CheckChildCodec(SesaFormat.ChildCodecPacked, int.MaxValue - 1L), Throws.Nothing);
        Assert.That(() => SesaHeaderReader.CheckChildCodec(3, 10), Throws.TypeOf<NotSupportedException>());
        Assert.That(() => SesaHeaderReader.CheckChildCodec(0, 10), Throws.TypeOf<NotSupportedException>());
    }

    [Test]
    public void TheBuilder_ChoosesTheCodec_FromTheRequestAndTheLengthOfTheText()
    {
        long max = TaggedChild.MaxLength;
        Assert.That(SesaBuilder.ChildCodecFor(0, null), Is.EqualTo(ChildCodec.Tagged));
        Assert.That(SesaBuilder.ChildCodecFor(1000, null), Is.EqualTo(ChildCodec.Tagged));
        Assert.That(SesaBuilder.ChildCodecFor(max, null), Is.EqualTo(ChildCodec.Tagged), "the largest text whose ranks fit 30 bits");
        Assert.That(SesaBuilder.ChildCodecFor(max + 1, null), Is.EqualTo(ChildCodec.Packed), "above it the default is the packed codec");
        Assert.That(SesaBuilder.ChildCodecFor(3_000_000_000, null), Is.EqualTo(ChildCodec.Packed));
        Assert.That(SesaBuilder.ChildCodecFor(1000, ChildCodec.Packed), Is.EqualTo(ChildCodec.Packed));
        Assert.That(SesaBuilder.ChildCodecFor(1000, ChildCodec.Tagged), Is.EqualTo(ChildCodec.Tagged));
        Assert.That(SesaBuilder.ChildCodecFor(max, ChildCodec.Tagged), Is.EqualTo(ChildCodec.Tagged));
        Assert.That(() => SesaBuilder.ChildCodecFor(max + 1, ChildCodec.Tagged), Throws.TypeOf<NotSupportedException>().With.Message.Contains("30 bits"));
        Assert.That(SesaBuilder.ChildCodecFor(max + 1, ChildCodec.Packed), Is.EqualTo(ChildCodec.Packed));
    }

    [Test]
    public void TheStage_RefusesTheTaggedCodec_ForATextBeyondItsRanks_BeforeItTouchesAnyFile()
    {
        using var dir = new TempDir();
        long n = TaggedChild.MaxLength + 1L; // 2^30 symbols: the ranks need 31 bits
        string sparse = dir.File("long.pk");
        using (var fs = new FileStream(sparse, FileMode.Create)) fs.SetLength(PackedTextFormat.TotalBytes(n, 1));
        using MappedTextReader reader = MappedTextReader.Open(sparse, n, 1);
        using MappedOutputFile file = MappedOutputFile.Create(dir.File("child.bin"), 4096);
        Func<LcpPlacementRequest, LcpPlacement> place = _ => throw new InvalidOperationException("the stage must refuse before it asks for the LCP");
        LcpChildOptions tagged = LcpChildRun.Options(dir.Sub("t"), LcpChildRun.TinyTuning(1), codec: ChildCodec.Tagged);
        Assert.That(() => SectionLcpChildBuilder.Build("no-such-suffix-array", 0, reader, file, 128, place, tagged), Throws.TypeOf<NotSupportedException>().With.Message.Contains("30 bits"));
    }

    [TestCase("mississippi")]
    [TestCase("abracadabra")]
    [TestCase("aaaaaaaaaaaaaaaaaaaa")]
    [TestCase("ACGTACGTTGCAACGT")]
    [TestCase("")]
    [TestCase("x")]
    public void TheBuilder_WritesTheTaggedCodecByDefault_AndThePackedOnRequest_BothReadable(string text)
    {
        using var dir = new TempDir();
        LcpChildCase c = LcpChildCase.Of(text);
        foreach ((ChildCodec? requested, ChildCodec expected) in new (ChildCodec?, ChildCodec)[] { (null, ChildCodec.Tagged), (ChildCodec.Tagged, ChildCodec.Tagged), (ChildCodec.Packed, ChildCodec.Packed) })
        {
            SesaRun.Outcome o = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, ChildCodec = requested }, $"-{requested}");
            byte[] bytes = File.ReadAllBytes(o.Target);
            OracleSection child = SesaOracle.FromBytes(bytes).Get(SesaFormat.KindChild);
            Assert.That(child.Codec, Is.EqualTo(expected == ChildCodec.Tagged ? SesaFormat.ChildCodecTagged : SesaFormat.ChildCodecPacked), $"{requested}: codec");
            Assert.That(MemoryMarshal.Cast<byte, uint>(bytes.AsSpan((int)child.Offset, (int)child.Length)).ToArray(), Is.EqualTo(c.Child(expected)), $"{requested}: cells");
            Assert.That(Encoding.UTF8.GetString(SesaOracle.FromBytes(bytes).Payload(11)), Does.Contain(expected == ChildCodec.Tagged ? "child=layout A, tagged" : "child=layout A, packed"), $"{requested}: BUILD_INFO");
            using EsaSuffixTree tree = SesaTestFile.Open(o.Target, EsaValidationLevel.Structure, EsaAccessMode.Checked, sha: true);
            FacadeChecks.CheckAgainstTree(tree, EsaSuffixTree.Build(text), "esa-memory", text, $"built, {requested}");
        }
    }

    [Test]
    public void TwoBuilds_OfTheSameText_GiveTheSameBytes_ForEveryCodec()
    {
        using var dir = new TempDir();
        string text = StringFamilies.Random(3000, 4, 9);
        foreach (ChildCodec codec in new[] { ChildCodec.Packed, ChildCodec.Tagged })
        {
            SesaRun.Outcome a = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, ChildCodec = codec }, $"-a{codec}");
            SesaRun.Outcome b = SesaRun.Build(dir, text, new SesaBuildOptions { Fsync = false, ChildCodec = codec, Threads = 2 }, $"-b{codec}");
            Assert.That(File.ReadAllBytes(b.Target), Is.EqualTo(File.ReadAllBytes(a.Target)), $"{codec}");
        }
    }
}
