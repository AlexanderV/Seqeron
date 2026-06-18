// QUALITY-STATS-001 — Quality Statistics
// Evidence: docs/Evidence/QUALITY-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/QUALITY-STATS-001.md
// Source: Ewing & Green (1998) Genome Research 8(3):186-194; Illumina "Sequencing Quality Scores";
//         Newcastle Univ. ASK (population stddev); Math is Fun (median); Cock et al. (2010) NAR 38(6).

using NUnit.Framework;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class QualityScoreAnalyzer_CalculateStatistics_Tests
{
    private const QualityScoreAnalyzer.QualityEncoding Phred33 = QualityScoreAnalyzer.QualityEncoding.Phred33;
    private const QualityScoreAnalyzer.QualityEncoding Phred64 = QualityScoreAnalyzer.QualityEncoding.Phred64;

    // Phred+33: '5'(53)->20, '?'(63)->30, 'I'(73)->40 (Cock et al. 2010, ord-33).
    private const string Triplet = "5?I"; // decoded 20, 30, 40

    #region CalculateStatistics(string)

    // M1 — Phred decode 20,30,40 → mean 30, min 20, max 40.
    [Test]
    public void CalculateStatistics_Phred33Triplet_MeanMinMax()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(Triplet, Phred33);

        Assert.Multiple(() =>
        {
            Assert.That(s.MeanQuality, Is.EqualTo(30.0).Within(1e-10),
                "(20+30+40)/3 = 30 (arithmetic mean of decoded Phred scores).");
            Assert.That(s.MinQuality, Is.EqualTo(20), "min decoded score is 20 ('5').");
            Assert.That(s.MaxQuality, Is.EqualTo(40), "max decoded score is 40 ('I').");
        });
    }

    // M2 — Population std dev σ=√((1/N)Σ(x-μ)²)=√(200/3) (Newcastle, ÷N).
    [Test]
    public void CalculateStatistics_Phred33Triplet_PopulationStdDev()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(Triplet, Phred33);

        Assert.That(s.StandardDeviation, Is.EqualTo(8.16496580927726).Within(1e-10),
            "deviations -10,0,10 → variance 200/3 → σ=√(200/3); population divisor ÷N per Newcastle ASK.");
    }

    // M3 — Median odd count (n=3): middle of sorted 20,30,40 = 30 (Math is Fun).
    [Test]
    public void CalculateStatistics_OddCount_MedianMiddle()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(Triplet, Phred33);

        Assert.That(s.MedianQuality, Is.EqualTo(30.0).Within(1e-10),
            "odd-length median is the single middle order statistic (30).");
    }

    // M4 — Median even count (n=4): "5II?" → 20,40,40,30 → sorted 20,30,40,40 → (30+40)/2 (Math is Fun).
    [Test]
    public void CalculateStatistics_EvenCount_MedianAveragesTwoCentral()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics("5II?", Phred33); // 20,40,40,30

        Assert.That(s.MedianQuality, Is.EqualTo(35.0).Within(1e-10),
            "even-length median averages the two central order statistics: (30+40)/2 = 35.");
    }

    // M5 — %≥Q30 inclusive: 2 of 3 bases (30,40) ≥30 (Illumina).
    [Test]
    public void CalculateStatistics_Phred33Triplet_PercentAboveQ30()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(Triplet, Phred33);

        Assert.That(s.PercentAboveQ30, Is.EqualTo(66.66666666666667).Within(1e-10),
            "2 of 3 bases (Q30 and Q40) are ≥ Q30 (inclusive) → 100·2/3.");
    }

    // M6 — %≥Q20 inclusive: all 3 bases ≥20 (Illumina).
    [Test]
    public void CalculateStatistics_Phred33Triplet_PercentAboveQ20()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(Triplet, Phred33);

        Assert.That(s.PercentAboveQ20, Is.EqualTo(100.0).Within(1e-10),
            "all 3 bases (20,30,40) are ≥ Q20 → 100%.");
    }

    // M8 — TotalBases / BasesAboveQ30 / BasesAboveQ20 counts.
    [Test]
    public void CalculateStatistics_Phred33Triplet_Counts()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(Triplet, Phred33);

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalBases, Is.EqualTo(3), "3 quality characters.");
            Assert.That(s.BasesAboveQ30, Is.EqualTo(2), "scores 30 and 40 are ≥ Q30.");
            Assert.That(s.BasesAboveQ20, Is.EqualTo(3), "all scores ≥ Q20.");
        });
    }

    // M9 — INV-05: identical decoded scores under a different encoding yield identical statistics.
    // Phred+64: 't'(116)->52, '~'(126)->62 — different scores; build a P64 string decoding to 20,30,40:
    // '4'? No — use offset: Q+64 → 20→'T'(84),30→'^'(94),40→'h'(104).
    [Test]
    public void CalculateStatistics_EncodingInvariance_SameScoresSameStats()
    {
        var p33 = QualityScoreAnalyzer.CalculateStatistics(Triplet, Phred33);            // 20,30,40
        var p64 = QualityScoreAnalyzer.CalculateStatistics("T^h", Phred64);              // 84-64,94-64,104-64 = 20,30,40

        Assert.Multiple(() =>
        {
            Assert.That(p64.MeanQuality, Is.EqualTo(p33.MeanQuality).Within(1e-10),
                "statistics depend only on decoded scores, not the encoding (INV-05).");
            Assert.That(p64.PercentAboveQ30, Is.EqualTo(p33.PercentAboveQ30).Within(1e-10),
                "same decoded scores → same %≥Q30 regardless of encoding.");
        });
    }

    // S1 — single base "I" (Q40): zero spread identities.
    [Test]
    public void CalculateStatistics_SingleBase_ZeroStdDev()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics("I", Phred33); // 40

        Assert.Multiple(() =>
        {
            Assert.That(s.MeanQuality, Is.EqualTo(40.0).Within(1e-10), "single value → mean 40.");
            Assert.That(s.MedianQuality, Is.EqualTo(40.0).Within(1e-10), "single value → median 40.");
            Assert.That(s.MinQuality, Is.EqualTo(40), "single value → min 40.");
            Assert.That(s.MaxQuality, Is.EqualTo(40), "single value → max 40.");
            Assert.That(s.StandardDeviation, Is.EqualTo(0.0).Within(1e-10), "no spread → σ = 0 (INV-02).");
        });
    }

    // S2 — all bases ≥ Q30: "?I" (30,40) → 100%.
    [Test]
    public void CalculateStatistics_AllAboveQ30_HundredPercent()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics("?I", Phred33); // 30,40

        Assert.That(s.PercentAboveQ30, Is.EqualTo(100.0).Within(1e-10),
            "both bases (30,40) are ≥ Q30 → 100%.");
    }

    // S3 — no bases ≥ Q30 but all ≥ Q20: "5" (20) → Q30%=0, Q20%=100 (INV-03).
    [Test]
    public void CalculateStatistics_NoneAboveQ30_ZeroPercent()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics("5", Phred33); // 20

        Assert.Multiple(() =>
        {
            Assert.That(s.PercentAboveQ30, Is.EqualTo(0.0).Within(1e-10), "Q20 base is below Q30 → 0%.");
            Assert.That(s.PercentAboveQ20, Is.EqualTo(100.0).Within(1e-10), "Q20 base is ≥ Q20 → 100% (INV-03).");
        });
    }

    // S4 — base exactly at Q30 ('?'→30) is counted (inclusive ≥, Illumina).
    [Test]
    public void CalculateStatistics_ExactlyQ30_CountedInclusive()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics("?", Phred33); // 30

        Assert.That(s.PercentAboveQ30, Is.EqualTo(100.0).Within(1e-10),
            "a base at exactly Q30 is counted (threshold is inclusive ≥30).");
    }

    // Edge: empty string → zeroed result (documented contract).
    [Test]
    public void CalculateStatistics_Empty_ReturnsZeroedResult()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics("", Phred33);

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalBases, Is.EqualTo(0), "no bases.");
            Assert.That(s.MeanQuality, Is.EqualTo(0.0).Within(1e-10), "empty → zeroed mean.");
            Assert.That(s.PercentAboveQ30, Is.EqualTo(0.0).Within(1e-10), "empty → 0% Q30.");
        });
    }

    // Edge: null string → zeroed result, no throw (documented contract).
    [Test]
    public void CalculateStatistics_Null_ReturnsZeroedResult()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics((string)null!, Phred33);

        Assert.That(s.TotalBases, Is.EqualTo(0), "null quality string returns a zeroed result without throwing.");
    }

    #endregion

    #region CalculateQ30Percentage

    // M7 — INV-04: CalculateQ30Percentage equals CalculateStatistics(...).PercentAboveQ30.
    [Test]
    public void CalculateQ30Percentage_MatchesStatisticsPercentAboveQ30()
    {
        double q30 = QualityScoreAnalyzer.CalculateQ30Percentage(Triplet, Phred33);
        var s = QualityScoreAnalyzer.CalculateStatistics(Triplet, Phred33);

        Assert.Multiple(() =>
        {
            Assert.That(q30, Is.EqualTo(66.66666666666667).Within(1e-10),
                "2 of 3 bases ≥ Q30 → 100·2/3.");
            Assert.That(q30, Is.EqualTo(s.PercentAboveQ30).Within(1e-10),
                "Q30 percentage == PercentAboveQ30 of the full statistics (INV-04).");
        });
    }

    // Edge: empty/null → 0 for Q30 percentage.
    [Test]
    public void CalculateQ30Percentage_EmptyOrNull_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(QualityScoreAnalyzer.CalculateQ30Percentage("", Phred33), Is.EqualTo(0.0).Within(1e-10),
                "empty input → 0% Q30.");
            Assert.That(QualityScoreAnalyzer.CalculateQ30Percentage(null!, Phred33), Is.EqualTo(0.0).Within(1e-10),
                "null input → 0% Q30 (no throw).");
        });
    }

    #endregion

    #region CalculateStatistics(IEnumerable<string>) — delegate smoke

    // C1 — multi-read overload aggregates scores across reads.
    [Test]
    public void CalculateStatistics_MultiReadDelegate_AggregatesScores()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(new[] { Triplet, Triplet }, Phred33); // 6 scores: 20,30,40 x2

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalBases, Is.EqualTo(6), "two 3-base reads → 6 aggregated scores.");
            Assert.That(s.MeanQuality, Is.EqualTo(30.0).Within(1e-10), "aggregated mean of {20,30,40,20,30,40} = 30.");
            Assert.That(s.PercentAboveQ30, Is.EqualTo(66.66666666666667).Within(1e-10),
                "4 of 6 aggregated bases ≥ Q30 → 100·4/6.");
        });
    }

    #endregion
}
