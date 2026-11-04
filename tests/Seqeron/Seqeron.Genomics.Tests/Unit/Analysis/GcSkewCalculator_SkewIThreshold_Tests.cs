// SEQ-REPLICATION-001 — finisher A2-3 (F31): SkewIT per-genus SkewI thresholds (Lu & Salzberg 2020,
// PLoS Comput Biol 16:e1008439). Table rows below are copied verbatim (CRLF, trailing "STDEV " space,
// empty threshold for genera with < 10 genomes, no final newline) from
// raw.githubusercontent.com/jenniferlu717/SkewIT/master/data/RefSeq97_Bacteria_GenusSkewIThresholds.txt
// (1 147 genus rows, 160 with a threshold). The full table is not bundled: SkewIT is GPL-3.0.
// SkewI values are the output of SkewIT's own src/skewi.py (`-k 1000 --min-len 0`, %r output).
// Evidence: docs/Evidence/SEQ-REPLICATION-001-Evidence.md

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GcSkewCalculator_SkewIThreshold_Tests
{
    private const string TableExcerpt =
        "Genus\tNum_Genomes\tMean\tSTDEV \tThreshold\r\n" +
        "g__Escherichia\t934\t0.8486\t0.0688\t0.7110\r\n" +
        "g__Bordetella\t618\t0.4988\t0.1394\t0.2200\r\n" +
        "g__Mycobacterium\t260\t0.7460\t0.1751\t0.3959\r\n" +
        "g__Streptomyces\t181\t0.3352\t0.1447\t0.046\r\n" +
        "g__Synechococcus\t29\t0.4288\t0.3256\t-0.222\r\n" +
        "g__Synechocystis\t9\t0.1010\t\t\r\n" +
        "g__Thermincola\t1\t1.0000\t\t";

    private const string Ba1fExtraResource = "Seqeron.Genomics.Tests.TestData.Rosalind.ba1f_extra_dataset.txt";

    private static IReadOnlyDictionary<string, double> Table() =>
        GcSkewCalculator.ParseSkewIGenusThresholds(new StringReader(TableExcerpt));

    private static string Ba1fExtraGenome()
    {
        var asm = typeof(GcSkewCalculator_SkewIThreshold_Tests).Assembly;
        using Stream stream = asm.GetManifestResourceStream(Ba1fExtraResource)
            ?? throw new InvalidOperationException($"Embedded resource '{Ba1fExtraResource}' not found.");
        using var reader = new StreamReader(stream);
        return reader.ReadToEnd().Split('\n')[1].Trim();
    }

    private static string SyntheticGenome() =>
        string.Concat(Enumerable.Repeat("ACGTC", 6000))
        + string.Concat(Enumerable.Repeat("AGGTC", 10000))
        + string.Concat(Enumerable.Repeat("ACGTC", 4000));

    // T1 — Parse: one entry per row with a threshold; values exactly as published; rows without a
    // threshold (< 10 genomes) are absent; g__ prefix stripped; case-insensitive.
    [Test]
    public void ParseSkewIGenusThresholds_VerbatimRows_ValuesLockedToFile()
    {
        var t = Table();
        Assert.Multiple(() =>
        {
            Assert.That(t, Has.Count.EqualTo(5));
            Assert.That(t["Escherichia"], Is.EqualTo(0.7110));
            Assert.That(t["Bordetella"], Is.EqualTo(0.2200));
            Assert.That(t["Mycobacterium"], Is.EqualTo(0.3959));
            Assert.That(t["Streptomyces"], Is.EqualTo(0.046));
            Assert.That(t["Synechococcus"], Is.EqualTo(-0.222));
            Assert.That(t.ContainsKey("Synechocystis"), Is.False);
            Assert.That(t.ContainsKey("Thermincola"), Is.False);
            Assert.That(t["escherichia"], Is.EqualTo(0.7110));
        });
    }

    // T2 — LF endings parse identically; malformed threshold / duplicate genus → FormatException.
    [Test]
    public void ParseSkewIGenusThresholds_LineEndingsAndMalformedInput()
    {
        var lf = GcSkewCalculator.ParseSkewIGenusThresholds(new StringReader(TableExcerpt.Replace("\r\n", "\n")));
        Assert.Multiple(() =>
        {
            Assert.That(lf, Is.EquivalentTo(Table()));
            Assert.Throws<FormatException>(() => GcSkewCalculator.ParseSkewIGenusThresholds(
                new StringReader("Genus\tNum_Genomes\tMean\tSTDEV \tThreshold\ng__X\t10\t0.5\t0.1\tabc")));
            Assert.Throws<FormatException>(() => GcSkewCalculator.ParseSkewIGenusThresholds(
                new StringReader("g__X\t10\t0.5\t0.1\t0.3\ng__x\t10\t0.5\t0.1\t0.3")));
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.ParseSkewIGenusThresholds(null!));
        });
    }

    // T3 — Lookup: exact genus, case-insensitive, optional g__ prefix; unknown / blank → false.
    [Test]
    public void TryGetSkewIThreshold_ExactCaseInsensitiveMatch()
    {
        var t = Table();
        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.TryGetSkewIThreshold(t, "ESCHERICHIA", out double e), Is.True);
            Assert.That(e, Is.EqualTo(0.7110));
            Assert.That(GcSkewCalculator.TryGetSkewIThreshold(t, "g__Bordetella", out double b), Is.True);
            Assert.That(b, Is.EqualTo(0.2200));
            Assert.That(GcSkewCalculator.TryGetSkewIThreshold(t, "Escherich", out _), Is.False);
            Assert.That(GcSkewCalculator.TryGetSkewIThreshold(t, "Synechocystis", out _), Is.False);
            Assert.That(GcSkewCalculator.TryGetSkewIThreshold(t, "", out _), Is.False);
            // caller-built, case-sensitive map
            var cs = new Dictionary<string, double> { ["Escherichia"] = 0.7110 };
            Assert.That(GcSkewCalculator.TryGetSkewIThreshold(cs, "escherichia", out double c), Is.True);
            Assert.That(c, Is.EqualTo(0.7110));
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.TryGetSkewIThreshold(null!, "X", out _));
        });
    }

    // T4 — Decision (skewi.py -k 1000): BA1F extra dataset SkewI 0.203158581311549 is below Escherichia
    // 0.7110, Mycobacterium 0.3959 and Bordetella 0.2200, not below Streptomyces 0.046 or Synechococcus
    // −0.222. Ideal synthetic genome SkewI 1.0 → never below.
    [TestCase("Escherichia", true)]
    [TestCase("Mycobacterium", true)]
    [TestCase("Bordetella", true)]
    [TestCase("Streptomyces", false)]
    [TestCase("Synechococcus", false)]
    public void IsSkewIBelowGenusThreshold_Ba1fExtra_MatchesSkewItRule(string genus, bool expected)
    {
        Assert.That(GcSkewCalculator.IsSkewIBelowGenusThreshold(Ba1fExtraGenome(), genus, Table(), 1000),
            Is.EqualTo(expected));
        Assert.That(GcSkewCalculator.IsSkewIBelowGenusThreshold(SyntheticGenome(), genus, Table(), 1000),
            Is.False);
    }

    // T5 — Explicit threshold: strict "<" at the skewi.py value; null SkewI / unknown genus → null; guards.
    [Test]
    public void IsSkewIBelowThreshold_StrictComparisonNullsAndGuards()
    {
        string ba1f = Ba1fExtraGenome();
        const double skewI = 0.203158581311549; // skewi.py -k 1000
        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.IsSkewIBelowThreshold(ba1f, skewI, 1000), Is.False);
            Assert.That(GcSkewCalculator.IsSkewIBelowThreshold(ba1f, 0.2032, 1000), Is.True);
            Assert.That(GcSkewCalculator.IsSkewIBelowThreshold(ba1f, 0.2031, 1000), Is.False);
            // 93 523 bp / 20 kb = 5 windows → skewi.py reports nothing → null.
            Assert.That(GcSkewCalculator.IsSkewIBelowThreshold(ba1f, 0.7110), Is.Null);
            Assert.That(GcSkewCalculator.IsSkewIBelowGenusThreshold(ba1f, "Escherichia", Table()), Is.Null);
            Assert.That(GcSkewCalculator.IsSkewIBelowGenusThreshold(ba1f, "Synechocystis", Table(), 1000), Is.Null);
            Assert.That(GcSkewCalculator.IsSkewIBelowGenusThreshold(ba1f, "Unknownia", Table(), 1000), Is.Null);
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.IsSkewIBelowThreshold(ba1f, double.NaN, 1000));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.IsSkewIBelowThreshold(ba1f, 0.5, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.IsSkewIBelowGenusThreshold(ba1f, "Escherichia", Table(), 0));
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.IsSkewIBelowGenusThreshold(ba1f, "Escherichia", null!, 1000));
        });
    }
}
