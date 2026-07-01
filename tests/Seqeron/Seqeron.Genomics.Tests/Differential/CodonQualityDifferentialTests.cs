// 08_DIFFERENTIAL_TESTING rows 214, 219, 220. Independent oracles: hand-derived RSCU on controlled
// two-codon families, ASCII-offset Phred decoding, and manual quality statistics.

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.IO;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class CodonQualityDifferentialTests
{
    private const double Tol = 1e-9;

    // ---- Row 214: CODON-RSCU-001 — RSCU vs hand-derived synonymous ratios ----

    [Test]
    [Category("CODON-RSCU-001")]
    public void Rscu_MatchesHandDerivedValues()
    {
        // Codons: Phe{TTT,TTC} both once; Lys{AAA,AAG} only AAA; Asp{GAT,GAC} only GAT.
        // RSCU_j = observed / (totalSyn / nSyn).
        var rscu = CodonUsageAnalyzer.CalculateRscu("TTTTTCAAAGAT");

        Assert.That(rscu["TTT"], Is.EqualTo(1.0).Within(Tol)); // 1 / (2/2)
        Assert.That(rscu["TTC"], Is.EqualTo(1.0).Within(Tol));
        Assert.That(rscu["AAA"], Is.EqualTo(2.0).Within(Tol)); // 1 / (1/2)
        Assert.That(rscu["AAG"], Is.EqualTo(0.0).Within(Tol)); // 0 / (1/2)
        Assert.That(rscu["GAT"], Is.EqualTo(2.0).Within(Tol));
        Assert.That(rscu["GAC"], Is.EqualTo(0.0).Within(Tol));
    }

    // ---- Row 219: QUALITY-PHRED-001 — Phred33 = ASCII − 33 ----

    [Test]
    [Category("QUALITY-PHRED-001")]
    [TestCase("IIII")]
    [TestCase("!#5I")]
    [TestCase("")]
    public void ParseQualityString_MatchesAsciiOffset(string qual)
    {
        int[] expected = qual.Select(c => c - 33).ToArray();
        Assert.That(QualityScoreAnalyzer.ParseQualityString(qual), Is.EqualTo(expected));
    }

    // ---- Row 220: QUALITY-STATS-001 — quality statistics vs manual ----

    [Test]
    [Category("QUALITY-STATS-001")]
    [TestCase("!#5I")]              // scores 0,2,20,40
    [TestCase("IIIIIIIIII")]        // all 40
    [TestCase("$$$$5555")]          // 3,3,3,3,20,20,20,20
    public void QualityStatistics_MatchManual(string qual)
    {
        var scores = qual.Select(c => c - 33).ToArray();
        double mean = scores.Average();
        var sorted = scores.OrderBy(x => x).ToArray();
        double median = sorted.Length % 2 == 0
            ? (sorted[sorted.Length / 2 - 1] + sorted[sorted.Length / 2]) / 2.0
            : sorted[sorted.Length / 2];
        double std = Math.Sqrt(scores.Select(x => Math.Pow(x - mean, 2)).Average());

        var st = QualityScoreAnalyzer.CalculateStatistics(qual);
        Assert.That(st.MeanQuality, Is.EqualTo(mean).Within(Tol));
        Assert.That(st.MedianQuality, Is.EqualTo(median).Within(Tol));
        Assert.That(st.MinQuality, Is.EqualTo(scores.Min()));
        Assert.That(st.MaxQuality, Is.EqualTo(scores.Max()));
        Assert.That(st.StandardDeviation, Is.EqualTo(std).Within(Tol));
        Assert.That(st.TotalBases, Is.EqualTo(scores.Length));
        Assert.That(st.BasesAboveQ20, Is.EqualTo(scores.Count(q => q >= 20)));
        Assert.That(st.BasesAboveQ30, Is.EqualTo(scores.Count(q => q >= 30)));
    }
}
