// SEQ-REPLICATION-001 — finisher A1-3(b): SkewIT Skew Index (Lu & Salzberg 2020, PLoS Comput Biol
// 16:e1008439). Expected values are the output of the authors' own src/skewi.py
// (raw.githubusercontent.com/jenniferlu717/SkewIT/master/src/skewi.py, run under Python 3 with
// `--min-len 0`, header containing "complete", output format patched from %0.10f to %r for full precision).
// Evidence: docs/Evidence/SEQ-REPLICATION-001-Evidence.md

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GcSkewCalculator_SkewIndex_Tests
{
    private const string Ba1fExtraResource = "Seqeron.Genomics.Tests.TestData.Rosalind.ba1f_extra_dataset.txt";

    private static string Ba1fExtraGenome()
    {
        var asm = typeof(GcSkewCalculator_SkewIndex_Tests).Assembly;
        using Stream stream = asm.GetManifestResourceStream(Ba1fExtraResource)
            ?? throw new InvalidOperationException($"Embedded resource '{Ba1fExtraResource}' not found.");
        using var reader = new StreamReader(stream);
        return reader.ReadToEnd().Split('\n')[1].Trim();
    }

    private static string SyntheticGenome() =>
        string.Concat(Enumerable.Repeat("ACGTC", 6000))
        + string.Concat(Enumerable.Repeat("AGGTC", 10000))
        + string.Concat(Enumerable.Repeat("ACGTC", 4000));

    // S1 — Small sequences, k = 4 (skewi.py output). 13 windows: h = round(6.5) = 6 (half-to-even),
    // r = round(0.52) = 1; 15 windows with a partial tail; 41 windows.
    [TestCase("GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGG", 0.23076923076923078)]
    [TestCase("GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAAT", 0.6206896551724138)]
    [TestCase("GGGCCAATTGGGTCTCAAGCAAGTCGCGACGGTACAGGGGCCCAGCCTGGCTGCGCGGAGGGGACTGGGAGCTGTTGGTTAGTGGGGAGGACACGCACTAAATCTACTAACTGCCCTCGAAGGGGCACACCGCTACTCCTATCACTCCCCTTTCCTCCGCC", 0.6211180124223602)]
    [TestCase("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", 0.23076923076923078)]
    public void CalculateSkewIndex_SmallSequences_MatchSkewItScript(string seq, double expected)
    {
        Assert.That(GcSkewCalculator.CalculateSkewIndex(seq, 4), Is.EqualTo(expected).Within(1e-15));
        Assert.That(GcSkewCalculator.CalculateSkewIndex(new DnaSequence(seq), 4), Is.EqualTo(expected).Within(1e-15));
    }

    // S2 — skewi.py writes no SkewI when maxDiff ≤ 0: 12 windows (r = round(0.48) = 0) or all windows
    // sign 0 (no G/C) → null.
    [TestCase("GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC")]
    [TestCase("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")]
    public void CalculateSkewIndex_NoValueInSkewIt_ReturnsNull(string seq)
    {
        Assert.That(GcSkewCalculator.CalculateSkewIndex(seq, 4), Is.Null);
    }

    // S3 — BA1F extra dataset (93 523 bp), skewi.py -k 1000 → 0.203158581311549; -k 20 → 0.034430033253851994.
    [TestCase(1000, 0.203158581311549)]
    [TestCase(20, 0.034430033253851994)]
    public void CalculateSkewIndex_Ba1fExtraDataset_MatchesSkewItScript(int k, double expected)
    {
        Assert.That(GcSkewCalculator.CalculateSkewIndex(Ba1fExtraGenome(), k), Is.EqualTo(expected).Within(1e-15));
    }

    // S4 — Ideal two-strand genome: skewi.py caps at 1.0 (-k 1000 and -k 20).
    [TestCase(1000)]
    [TestCase(20)]
    public void CalculateSkewIndex_IdealSkewGenome_CappedAtOne(int k)
    {
        Assert.That(GcSkewCalculator.CalculateSkewIndex(SyntheticGenome(), k), Is.EqualTo(1.0));
    }

    // S5 — Default window is SkewIT's 20 kb; guards; null/empty; case-insensitive (documented deviation:
    // skewi.py counts only upper-case G/C).
    [Test]
    public void CalculateSkewIndex_DefaultsGuardsAndCase()
    {
        string g = SyntheticGenome();
        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.DefaultSkewIndexWindow, Is.EqualTo(20000));
            // 100 kb / 20 kb = 5 windows < 13 → r = 0 → skewi.py reports nothing.
            Assert.That(GcSkewCalculator.CalculateSkewIndex(g), Is.Null);
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateSkewIndex("ACGT", 0));
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.CalculateSkewIndex((DnaSequence)null!));
            Assert.That(GcSkewCalculator.CalculateSkewIndex((string)null!, 4), Is.Null);
            Assert.That(GcSkewCalculator.CalculateSkewIndex(string.Empty, 4), Is.Null);
            Assert.That(GcSkewCalculator.CalculateSkewIndex(new string('g', 52), 4), Is.EqualTo(0.23076923076923078).Within(1e-15));
        });
    }
}
