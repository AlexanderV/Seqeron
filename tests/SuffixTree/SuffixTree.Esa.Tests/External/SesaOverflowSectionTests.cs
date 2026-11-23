using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// Gap of the audit of 1b.7e (mutant <c>O1</c>): no unit test had an LCP_OVERFLOW list of 131,072 pairs (1 MiB), the size from which ADR-003 section 3 wants the section at a multiple of 4,096, so the alignment of that
/// section was checked only by the 10^8 runs. A repeat of B symbols adds B - 255 ranks whose LCP reaches the escape of a U8 LCP (and B - 65,535 that of a U16 LCP), so a random text of 2.1 million symbols with a repeat of 140,000 has an
/// overflow list of 139,745 pairs next to a U8 LCP (U8 needs the list to be at most n / 16 entries), and one of 1.5 million symbols with a repeat of 200,000 has 134,465 pairs next to a U16 LCP (an a^n text of 300,000 symbols,
/// which the brief of 1b.7f suggested, does not: its LCP is U32, which has no overflow list). The whole file is checked against ADR-003 and read back by the reader at every level.
/// </summary>
[TestFixture]
[HangLimit(30)]
public class SesaOverflowSectionTests
{
    private static string RepeatText(int random, int repeat, ulong seed)
    {
        string a = StringFamilies.RandomDna(random, seed);
        return string.Concat(a, a.AsSpan(random / 3, repeat));
    }

    [TestCase(2_100_000, 140_000, "U8")]
    [TestCase(1_500_000, 200_000, "U16")]
    [CancelAfter(600_000)]
    [HangLimit(60)]
    public void AnOverflowListOfAMebibyte_IsAlignedToTheLargeBoundary_AndTheWholeFileConforms(int random, int repeat, string kind)
    {
        LcpKind expectedKind = Enum.Parse<LcpKind>(kind);
        string text = RepeatText(random, repeat, 31);
        using var dir = new TempDir();
        string input = SesaRun.WriteInput(dir, text);
        string outDir = dir.Sub("out");
        Directory.CreateDirectory(outDir);
        string target = Path.Combine(outDir, "big.sesa");
        SesaBuildResult r = SesaBuilder.Build(input, target, new SesaBuildOptions { Threads = 2, Fsync = false, MeasureRss = false });
        Assert.That(r.LcpKind, Is.EqualTo(expectedKind));
        Assert.That(r.OverflowCount, Is.GreaterThanOrEqualTo(131_072), "the list is at least 1 MiB (131,072 pairs of 8 bytes)");

        SesaOracle f = SesaOracle.Read(target);
        OracleSection overflow = f.Get(SesaFormat.KindLcpOverflow);
        Assert.That(overflow.Length, Is.GreaterThanOrEqualTo(1 << 20));
        Assert.That(overflow.Offset % 4096, Is.Zero, "a section of 1 MiB or more starts at a multiple of 4,096 (ADR-003 section 3)");
        Assert.That(overflow.Offset, Is.GreaterThanOrEqualTo(f.Get(SesaFormat.KindLcp).Offset + f.Get(SesaFormat.KindLcp).Length));
        Assert.That(overflow.Offset + overflow.Length, Is.EqualTo(f.File.Length), "the file ends with the overflow list");
        foreach (OracleSection s in f.Sections) if (s.Length >= (1 << 20)) Assert.That(s.Offset % 4096, Is.Zero, $"section {s.Kind}");
        f.AssertConforms($"overflow {expectedKind}", psiExpected: true);

        // the same file through the reader, at the strongest level, against the in-memory tree on a sample of queries
        using EsaSuffixTree tree = SesaTestFile.Open(target, EsaValidationLevel.Structure, EsaAccessMode.Checked, sha: true);
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        Assert.That(tree.NodeCount, Is.EqualTo(memory.NodeCount));
        Assert.That(tree.LongestRepeatedSubstring().Length, Is.EqualTo(memory.LongestRepeatedSubstring().Length));
        var rng = new SplitMix64(5);
        for (int i = 0; i < 300; i++)
        {
            int len = 1 + rng.NextInt(40);
            string p = text.Substring(rng.NextInt(text.Length - len), len);
            Assert.That(tree.CountOccurrences(p), Is.EqualTo(memory.CountOccurrences(p)), p);
        }
        string longRepeat = text.Substring(random / 3, 2_000);
        Assert.That(tree.CountOccurrences(longRepeat), Is.EqualTo(memory.CountOccurrences(longRepeat)));
        Assert.That(tree.FindAllOccurrences(longRepeat), Is.EquivalentTo(memory.FindAllOccurrences(longRepeat)));
    }
}
