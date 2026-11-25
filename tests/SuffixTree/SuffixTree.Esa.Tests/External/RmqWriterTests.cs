using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Work package S2.5: the RMQ writer of the container builder (<see cref="SesaBuildOptions.IncludeRmq"/>, <see cref="RmqWriter"/>). The section is the pyramid of the LCP array as ADR-005 section 6.1 says (the oracle of the bytes is
/// <see cref="SesaTestFile.RmqBytes"/>, written from the definition and from nothing of the product), its entry is the one that the reader checks, every other section of the container is the one that the builder writes without the option,
/// and the container with it opens at every validation level and answers like the index in memory.
/// </summary>
[TestFixture]
[HangLimit(60)]
public class RmqWriterTests
{
    private static string RandomDna(int n, ulong seed) => StringFamilies.RandomDna(n, seed);

    public static IEnumerable<TestCaseData> Texts()
    {
        yield return new TestCaseData("banana").SetArgDisplayNames("banana");
        yield return new TestCaseData("mississippi").SetArgDisplayNames("mississippi");
        yield return new TestCaseData("a").SetArgDisplayNames("a");
        yield return new TestCaseData(RandomDna(31, 1)).SetArgDisplayNames("random 31");
        yield return new TestCaseData(RandomDna(32, 2)).SetArgDisplayNames("random 32");
        yield return new TestCaseData(RandomDna(33, 3)).SetArgDisplayNames("random 33");
        yield return new TestCaseData(RandomDna(1_023, 4)).SetArgDisplayNames("random 1023");
        yield return new TestCaseData(RandomDna(1_024, 5)).SetArgDisplayNames("random 1024 (1,025 LCP entries: 33 blocks at the first level)");
        yield return new TestCaseData(RandomDna(40_000, 6)).SetArgDisplayNames("random 40,000 (three levels)");
        yield return new TestCaseData(new string('a', 3_000)).SetArgDisplayNames("a^3000 (escapes of the 8-bit LCP array)");
        string a = RandomDna(5_000, 7);
        yield return new TestCaseData(a + a[1_000..4_000] + RandomDna(100, 8)).SetArgDisplayNames("random with a repeat of 3,000 (escapes, an overflow list)");
        string b = RandomDna(40_000, 9);
        yield return new TestCaseData(b + b[5_000..9_000]).SetArgDisplayNames("random with a repeat of 4,000 in 40,000");
    }

    private static (string Target, SesaBuildResult Result) Build(TempDir dir, string text, string name, bool rmq, bool psi = true)
    {
        string input = SesaRun.WriteInput(dir, text);
        string outDir = dir.Sub("out-" + name);
        Directory.CreateDirectory(outDir);
        string target = Path.Combine(outDir, name + ".sesa");
        SesaBuildResult r = SesaBuilder.Build(input, target, new SesaBuildOptions { Threads = 2, Fsync = false, MeasureRss = false, IncludeRmq = rmq, IncludePsi = psi });
        return (target, r);
    }

    [TestCaseSource(nameof(Texts))]
    public void TheSection_IsThePyramidOfTheLcpArray_WithTheEntryThatTheReaderChecks(string text)
    {
        using var dir = new TempDir();
        (string target, SesaBuildResult r) = Build(dir, text, "t", rmq: true);
        SesaOracle f = SesaOracle.Read(target);
        int[] lcp = EsaIndex.Build(text).Lcp.ToArray();
        byte[] expected = SesaTestFile.RmqBytes(lcp, SesaFormat.RmqFanout);
        if (text.Length == 0 || lcp.Length == 1)
        {
            Assert.That(f.Sections.Any(s => s.Kind == SesaFormat.KindRmq), Is.False, "the LCP array of a single entry has no level");
            return;
        }
        OracleSection s = f.Get(SesaFormat.KindRmq);
        Assert.That((s.Flags, s.Count, s.ElemWidth, s.Codec, s.Param0, s.Param1), Is.EqualTo((SesaFormat.FlagHasCrc, (long)expected.Length / 4, 4, 1, (long)SesaFormat.RmqFanout, 0L)), "optional, with a CRC, 4-byte entries, codec 1, PARAM0 the fanout");
        Assert.That(s.Length, Is.EqualTo(expected.Length));
        Assert.That(f.File.AsSpan((int)s.Offset, (int)s.Length).ToArray(), Is.EqualTo(expected), "the bytes of the section are the levels of the definition");
        Assert.That(s.Offset % (s.Length >= 1 << 20 ? 4096 : 64), Is.Zero);
        Assert.That(s.Offset + s.Length, Is.EqualTo(f.File.Length), "the section is the last: the file ends with it");
        Assert.That(r.FileSize, Is.EqualTo(f.File.Length));
        Assert.That(RmqWriter.CountFor(text.Length + 1), Is.EqualTo(expected.Length / 4));
        Assert.That(RmqWriter.LengthFor(text.Length + 1), Is.EqualTo(expected.Length));
    }

    [TestCaseSource(nameof(Texts))]
    public void EveryOtherSection_IsTheSectionOfTheBuildWithoutTheOption_AndTheContainerOpensAtEveryLevel(string text)
    {
        using var dir = new TempDir();
        (string with, _) = Build(dir, text, "with", rmq: true);
        (string without, _) = Build(dir, text, "without", rmq: false);
        SesaOracle a = SesaOracle.Read(with), b = SesaOracle.Read(without);
        foreach (OracleSection s in b.Sections)
        {
            OracleSection t = a.Get(s.Kind);
            Assert.That((t.Offset, t.Length, t.Count, t.Crc), Is.EqualTo((s.Offset, s.Length, s.Count, s.Crc)), $"section {s.Kind}");
            Assert.That(a.File.AsSpan((int)t.Offset, (int)t.Length).ToArray(), Is.EqualTo(b.File.AsSpan((int)s.Offset, (int)s.Length).ToArray()), $"section {s.Kind}");
        }
        foreach (EsaValidationLevel level in new[] { EsaValidationLevel.Header, EsaValidationLevel.Sections, EsaValidationLevel.Structure })
        {
            using EsaSuffixTree tree = SesaTestFile.Open(with, level, EsaAccessMode.Checked, sha: true);
            using EsaSuffixTree memory = EsaSuffixTree.Build(text);
            Assert.That(tree.NodeCount, Is.EqualTo(memory.NodeCount), $"{level}");
            Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo(memory.LongestRepeatedSubstring()), $"{level}");
        }
    }

    [Test]
    public void TheBuildWithoutTheOption_IsTheBuildOfS1_NoRmqSectionAndTheSameLength()
    {
        using var dir = new TempDir();
        string text = RandomDna(5_000, 11);
        (string target, SesaBuildResult r) = Build(dir, text, "plain", rmq: false);
        SesaOracle f = SesaOracle.Read(target);
        Assert.That(f.Sections.Any(s => s.Kind == SesaFormat.KindRmq), Is.False);
        Assert.That(r.Sections.Any(s => s.Kind == SesaFormat.KindRmq), Is.False);
        Assert.That(new SesaBuildOptions().IncludeRmq, Is.False, "off by default");
        f.AssertConforms("without RMQ", psiExpected: true);
    }

    [Test]
    public void TheOptionWorksWithoutPsi_AndTheContainerWithBothAnswersLinksFromTheSections()
    {
        using var dir = new TempDir();
        string text = RandomDna(2_000, 12);
        (string noPsi, _) = Build(dir, text, "nopsi", rmq: true, psi: false);
        SesaOracle f = SesaOracle.Read(noPsi);
        Assert.That(f.Sections.Any(s => s.Kind == SesaFormat.KindPsi), Is.False);
        Assert.That(f.Sections.Any(s => s.Kind == SesaFormat.KindRmq), Is.True);
        (string both, _) = Build(dir, text, "both", rmq: true, psi: true);
        var m = new EsaSuffixIndex(EsaIndex.Build(text));
        foreach (string path in new[] { noPsi, both })
        {
            using SesaFile file = SesaFile.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Header, Access = EsaAccessMode.Checked });
            var w = new EsaSuffixIndex(EsaIndex.FromFile(file));
            Navigation.WideNavigationTests.CompareAll(Path.GetFileName(path), m, w, text);
        }
    }

    [Test]
    public void TheOptionIsRefusedWithAResumableBuild()
    {
        using var dir = new TempDir();
        string input = SesaRun.WriteInput(dir, "banana");
        Assert.That(() => SesaBuilder.Build(input, dir.File("x.sesa"), new SesaBuildOptions { Threads = 1, IncludeRmq = true, Resume = true }), Throws.TypeOf<ArgumentException>());
    }

    [TestCase(1, 0)]
    [TestCase(2, 1)]
    [TestCase(32, 1)]
    [TestCase(33, 3)]
    [TestCase(64, 3)]
    [TestCase(1_024, 33)]
    [TestCase(1_025, 33 + 2 + 1)]
    public void TheNumberOfEntries_IsTheSumOfTheLevels(long rows, long expected) => Assert.That(RmqWriter.CountFor(rows), Is.EqualTo(expected));

    [Test]
    public void ABadFanoutOrARowCountIsRefused()
    {
        Assert.That(() => RmqWriter.CountFor(0), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => RmqWriter.CountFor(10, 1), Throws.TypeOf<ArgumentOutOfRangeException>());
        Assert.That(() => RmqWriter.CountFor(10, 12), Throws.TypeOf<ArgumentOutOfRangeException>());
    }
}
