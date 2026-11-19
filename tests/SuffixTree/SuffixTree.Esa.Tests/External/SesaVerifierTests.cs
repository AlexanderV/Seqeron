using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The structure pass of the reader (work package 1b.7f) at its limits: chunks of the text of a few positions, a scan stack that spills to a file after a few frames, the report of what it used, and the statistics it recomputes.
/// </summary>
[TestFixture]
public class SesaVerifierTests
{
    [TearDown]
    public void RestoreTuning() => SesaVerifierTuning.Reset();

    private static readonly string[] Texts =
    [
        "a", "ab", "banana", "mississippi", "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
        "abababababababababababababababababababababababababababab", "GATTACAGATTACAGGTAGATTACAGATTACAGATTACA", StringFamilies.Fibonacci(200), StringFamilies.ThueMorse(150), StringFamilies.Random(300, 3, 8),
    ];

    private static void OpenAtStructure(string path, bool sha = false)
    {
        using EsaSuffixTree tree = SesaTestFile.Open(path, EsaValidationLevel.Structure, EsaAccessMode.Checked, sha);
        Assert.That(tree.File!.Report.StructureSeconds, Is.GreaterThanOrEqualTo(0));
    }

    [TestCase(1, 2)]
    [TestCase(2, 2)]
    [TestCase(3, 3)]
    [TestCase(5, 4)]
    [TestCase(7, 9)]
    [TestCase(64, 64)]
    [TestCase(1 << 23, 1 << 22)]
    public void TinyChunksAndASpillingStack_AcceptEveryValidFile(int chunk, int stack)
    {
        using var dir = new TempDir();
        string spill = dir.Sub("spill");
        Directory.CreateDirectory(spill);
        SesaVerifierTuning.ChunkSymbols = chunk;
        SesaVerifierTuning.StackEntries = stack;
        SesaVerifierTuning.SpillDirectory = spill;
        foreach (string text in Texts)
        {
            string path = SesaTestFile.Write(dir, text, "t.sesa", psi: text.Length % 2 == 0);
            OpenAtStructure(path, sha: true);
            Assert.That(Directory.GetFileSystemEntries(spill), Is.Empty, "the spill file is removed when the pass ends");
        }
    }

    [Test]
    public void TheStackSpills_ToTheDirectoryOfTheTuning_WhenTheTreeIsDeeperThanTheLimit()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, new string('a', 400));
        SesaVerifierTuning.StackEntries = 8;
        SesaVerifierTuning.SpillDirectory = dir.Sub("does-not-exist");
        Assert.That(() => OpenAtStructure(path), Throws.InstanceOf<IOException>(), "a^400 has a chain of 400 intervals: with room for 8 the older ones go to the file, which cannot be created here");
        SesaVerifierTuning.StackEntries = 1000;
        Assert.That(() => OpenAtStructure(path), Throws.Nothing, "with room for all of them nothing is spilled");
    }

    [Test]
    public void TheReport_NamesTheBitmapAndTheBuffers_ThatTheLimitsBound()
    {
        using var dir = new TempDir();
        string text = StringFamilies.Random(5_000, 4, 4);
        string path = SesaTestFile.Write(dir, text);
        SesaVerifierTuning.ChunkSymbols = 100;
        SesaVerifierTuning.StackEntries = 64;
        using EsaSuffixTree tree = SesaTestFile.Open(path, EsaValidationLevel.Structure);
        SesaOpenReport r = tree.File!.Report;
        Assert.That(r.BitmapBytes, Is.EqualTo((5_000 + 1 + 63) / 64 * 8));
        Assert.That(r.BufferBytes, Is.InRange(8 * 100, (8 * 100) + (16 * 64) + (1 << 17)), "two int arrays of the chunk, the frames of the stack, a text block");
    }

    [Test]
    public void TheLongestRepeat_IsTheLeftBoundOfTheFirstDeepestInterval()
    {
        // ADR-003 as amended on 2026-11-19: LRS_RANK is the left bound lb; the first wording ("first rank r with LCP[r] == LRS_DEPTH") is lb + 1
        LcpChildCase c = LcpChildCase.Of("banana");
        Assert.That(c.Sa, Is.EqualTo(new[] { 6, 5, 3, 1, 0, 4, 2 }));
        byte[] image = SesaTestFile.Image(c);
        SesaOracle o = SesaOracle.FromBytes(image);
        Assert.That((o.LrsDepth, o.LrsRank), Is.EqualTo((3L, 2L)));
        Assert.That(c.Lcp.ToArray()[(int)o.LrsRank + 1], Is.EqualTo(o.LrsDepth), "LCP[LRS_RANK + 1] = LRS_DEPTH");
        using var dir = new TempDir();
        using EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure);
        Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("ana"));
        // a header with the other convention (lb + 1) is recomputed and refused at the structure level
        byte[] other = (byte[])image.Clone();
        BinaryPrimitives.WriteInt64LittleEndian(other.AsSpan(64), 3);
        SesaTestFile.ResealHeader(other);
        Assert.That(() => SesaTestFile.OpenImage(dir, other, EsaValidationLevel.Structure, name: "o.sesa").Dispose(), Throws.TypeOf<InvalidDataException>().With.Message.Contains("LRS_RANK"));
    }

    [Test]
    public void TheStatistics_OfAnEmptyAndOfASingleLetterText_AreAccepted()
    {
        using var dir = new TempDir();
        foreach (string text in new[] { string.Empty, "a", "aa", "ab" })
        {
            string path = SesaTestFile.Write(dir, text);
            OpenAtStructure(path, sha: true);
        }
    }

    [Test]
    public void ADeepChain_WithAWrongNodeCount_IsCaught_AtAnyTuning()
    {
        string text = new string('a', 150);
        LcpChildCase c = LcpChildCase.Of(text);
        byte[] image = SesaTestFile.Image(c);
        BinaryPrimitives.WriteInt64LittleEndian(image.AsSpan(24), BinaryPrimitives.ReadInt64LittleEndian(image.AsSpan(24)) - 1);
        SesaTestFile.ResealHeader(image);
        using var dir = new TempDir();
        foreach ((int chunk, int stack) in new[] { (1, 2), (13, 5), (1 << 23, 1 << 22) })
        {
            SesaVerifierTuning.ChunkSymbols = chunk;
            SesaVerifierTuning.StackEntries = stack;
            SesaVerifierTuning.SpillDirectory = dir.Path;
            Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure, name: "n.sesa").Dispose(), Throws.TypeOf<InvalidDataException>().With.Message.Contains("NODE_COUNT").Or.Message.Contains("nodes"), $"chunk {chunk} stack {stack}");
        }
    }

    [Test]
    public void ACarriedPrefix_OfMoreThanSixtyFourSymbols_IsVerifiedAtItsStart()
    {
        // two suffixes that share a long prefix are adjacent; swap the ranks of two LONG equal-prefix suffixes of a^n: the order check at the stored LCP sees it, and so does the Phi pass
        string text = new string('a', 300);
        LcpChildCase c = LcpChildCase.Of(text);
        for (int rank = 1; rank < 300; rank += 37)
        {
            byte[] image = SesaTestFile.Image(c);
            OracleSection sa = SesaOracle.FromBytes(image).Get(SesaFormat.KindSa);
            Span<byte> a = image.AsSpan((int)sa.Offset + (4 * rank), 4), b = image.AsSpan((int)sa.Offset + (4 * rank) + 4, 4);
            byte[] t = a.ToArray();
            b.CopyTo(a);
            t.CopyTo(b);
            SesaTestFile.Reseal(image);
            using var dir = new TempDir();
            Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Structure).Dispose(), Throws.TypeOf<InvalidDataException>(), $"rank {rank}");
        }
    }

    [Test]
    public void ASectionedLevel_WithAGapByte_IsCaughtOnlyFromTheSectionsLevel()
    {
        LcpChildCase c = LcpChildCase.Of(StringFamilies.Random(200, 4, 5));
        byte[] image = SesaTestFile.Image(c);
        SesaOracle o = SesaOracle.FromBytes(image);
        OracleSection text = o.Get(SesaFormat.KindText);
        long gap = text.Offset + text.Length; // the byte after TEXT: padding up to the alignment of the next section
        Assert.That(gap, Is.LessThan(o.Get(SesaFormat.KindAlphabet).Offset), "there is a gap");
        image[gap] = 7;
        using var dir = new TempDir();
        using (EsaSuffixTree tree = SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Header)) Assert.That(tree, Is.Not.Null);
        Assert.That(() => SesaTestFile.OpenImage(dir, image, EsaValidationLevel.Sections, name: "g.sesa").Dispose(), Throws.TypeOf<InvalidDataException>().With.Message.Contains("not zero"));
    }
}
