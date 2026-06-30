// 08_DIFFERENTIAL_TESTING rows 48, 50 (Chromosome). Independent oracles: a manual arm-ratio + Levan
// classification, and a manual exact-repeat telomere run count. (Rows 49/51/52 — centromere/aneuploidy/
// synteny — are deferred: large-window scans / loose "correlated" comparisons.)

using System;
using NUnit.Framework;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class Chromosome_Differential_Tests
{
    private const double Tol = 1e-12;

    // ---- Row 50: CHROM-KARYO-001 — arm ratio + classification vs manual Levan computation ----

    private static string ClassifyOracle(int centromere, int length)
    {
        double p = centromere, q = length - centromere;
        if (centromere <= 0 || length <= 0 || q <= 0) return "Telocentric";
        double ratio = p / q;
        if (ratio <= 0) return "Telocentric";
        double r = ratio >= 1.0 ? ratio : 1.0 / ratio; // long/short
        if (r <= 1.7) return "Metacentric";
        if (r <= 3.0) return "Submetacentric";
        if (r < 7.0) return "Subtelocentric";
        return "Acrocentric";
    }

    [Test]
    [Category("CHROM-KARYO-001")]
    [TestCase(50, 100, "Metacentric")]      // r = 1.0
    [TestCase(40, 100, "Metacentric")]      // r = 1.5
    [TestCase(30, 100, "Submetacentric")]   // r = 2.33
    [TestCase(20, 100, "Subtelocentric")]   // r = 4.0
    [TestCase(10, 100, "Acrocentric")]      // r = 9.0
    public void ArmRatioClassification_MatchesManualLevan(int centromere, int length, string expected)
    {
        double ratio = ChromosomeAnalyzer.CalculateArmRatio(centromere, length);
        Assert.That(ratio, Is.EqualTo((double)centromere / (length - centromere)).Within(Tol), "arm ratio = p/q");

        string actual = ChromosomeAnalyzer.ClassifyChromosomeByArmRatio(ratio);
        Assert.That(actual, Is.EqualTo(ClassifyOracle(centromere, length)), "matches independent Levan classifier");
        Assert.That(actual, Is.EqualTo(expected), "matches hand-known anchor");
    }

    // ---- Row 48: CHROM-TELO-001 — telomere repeat run vs manual exact-repeat count ----

    [Test]
    [Category("CHROM-TELO-001")]
    public void AnalyzeTelomeres_PerfectRepeats_MatchManualCount()
    {
        // 5' end: 3 perfect CCCTAA (revcomp of TTAGGG); 3' end: 4 perfect TTAGGG; non-telomere filler between.
        const string seq = "CCCTAACCCTAACCCTAA" + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + "TTAGGGTTAGGGTTAGGGTTAGGG";

        // Independent oracle: count leading exact CCCTAA and trailing exact TTAGGG 6-mers.
        int lead = 0;
        while ((lead + 1) * 6 <= seq.Length && seq.Substring(lead * 6, 6) == "CCCTAA") lead++;
        int trail = 0;
        while ((trail + 1) * 6 <= seq.Length && seq.Substring(seq.Length - (trail + 1) * 6, 6) == "TTAGGG") trail++;

        var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr1", seq);
        Assert.That(r.TelomereLength5Prime, Is.EqualTo(lead * 6), "5' telomere length");
        Assert.That(r.TelomereLength3Prime, Is.EqualTo(trail * 6), "3' telomere length");
        Assert.That(lead, Is.EqualTo(3));
        Assert.That(trail, Is.EqualTo(4));
        Assert.That(r.RepeatPurity5Prime, Is.EqualTo(1.0).Within(Tol), "perfect repeats -> purity 1");
        Assert.That(r.RepeatPurity3Prime, Is.EqualTo(1.0).Within(Tol));
    }
}
