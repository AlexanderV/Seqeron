// ONCO-MSI-001 — Microsatellite Instability (MSI) Detection
// Evidence: docs/Evidence/ONCO-MSI-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-MSI-001.md
// Source: Niu B et al. (2014). MSIsensor. Bioinformatics 30(7):1015-1016. https://doi.org/10.1093/bioinformatics/btt755
//         niu-lab/msisensor2 README (msi score = msi sites / valid sites; MSI-H: score >= 20%).
//         Boland CR et al. (1998). NCI Workshop on MSI. Cancer Res 58(22):5248-5257.

using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_DetectMSI_Tests
{
    #region CalculateMSIScore

    // M1 — 5 unstable / 25 valid = 0.20 (MSIsensor2: msi sites / valid sites).
    [Test]
    public void CalculateMSIScore_5over25_Returns0_20()
    {
        double score = OncologyAnalyzer.CalculateMSIScore(5, 25);

        Assert.That(score, Is.EqualTo(0.20).Within(1e-10),
            "MSI score = unstable/valid = 5/25 = 0.20 (MSIsensor2 fraction definition)");
    }

    // M2 — 3 / 12 = 0.25 (exact fraction, not the trivial value).
    [Test]
    public void CalculateMSIScore_3over12_Returns0_25()
    {
        double score = OncologyAnalyzer.CalculateMSIScore(3, 12);

        Assert.That(score, Is.EqualTo(0.25).Within(1e-10),
            "MSI score = 3/12 = 0.25 (MSIsensor2 fraction definition)");
    }

    // M3 — 0 / 25 = 0.0 (no unstable loci).
    [Test]
    public void CalculateMSIScore_ZeroUnstable_ReturnsZero()
    {
        double score = OncologyAnalyzer.CalculateMSIScore(0, 25);

        Assert.That(score, Is.EqualTo(0.0).Within(1e-10),
            "0 unstable / 25 valid = 0.0 (INV-01 lower bound)");
    }

    // M4 — 25 / 25 = 1.0 (INV-01 upper bound).
    [Test]
    public void CalculateMSIScore_AllUnstable_ReturnsOne()
    {
        double score = OncologyAnalyzer.CalculateMSIScore(25, 25);

        Assert.That(score, Is.EqualTo(1.0).Within(1e-10),
            "25/25 = 1.0; MSI score is bounded in [0,1] (INV-01)");
    }

    // S1 — zero valid loci: score u/n undefined (division by zero).
    [Test]
    public void CalculateMSIScore_ZeroValid_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateMSIScore(0, 0),
            "MSI score is undefined with no valid evaluable loci (division by zero)");
    }

    // S2 — unstable > valid is invalid (0 <= unstable <= valid).
    [Test]
    public void CalculateMSIScore_UnstableExceedsValid_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateMSIScore(6, 5),
            "unstable loci cannot exceed valid loci");
    }

    // S3 — negative unstable count is invalid.
    [Test]
    public void CalculateMSIScore_NegativeUnstable_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateMSIScore(-1, 5),
            "unstable loci count must be >= 0");
    }

    #endregion

    #region ClassifyMSIStatus

    // M5 — boundary 0.20 is MSI-High (MSIsensor2 cutoff is inclusive: ">= 20%").
    [Test]
    public void ClassifyMSIStatus_At20Percent_ReturnsMsiHigh()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyMSIStatus(0.20);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_High),
            "MSIsensor2: msi score >= 20% is MSI-H; the 20% boundary is inclusive");
    }

    // M6 — 0.16 (4/25) is below 20% -> MSS (not High).
    [Test]
    public void ClassifyMSIStatus_Below20Percent_ReturnsMss()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyMSIStatus(0.16);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSS),
            "0.16 < 0.20 -> not MSI-H (MSIsensor2 cutoff)");
    }

    // M7 — 0.40 is well above the cutoff -> MSI-High.
    [Test]
    public void ClassifyMSIStatus_At40Percent_ReturnsMsiHigh()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyMSIStatus(0.40);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_High),
            "0.40 >= 0.20 -> MSI-H (MSIsensor2 cutoff)");
    }

    // M8 — 0.0 is MSS.
    [Test]
    public void ClassifyMSIStatus_Zero_ReturnsMss()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyMSIStatus(0.0);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSS),
            "0.0 < 0.20 -> MSS");
    }

    // M5b — just below the 20% boundary (0.19999) is MSS: the cutoff is ">= 20%", so a value
    // strictly below 0.20 must NOT be MSI-H. This guards against an off-by-epsilon (> vs >=) regression.
    [Test]
    public void ClassifyMSIStatus_JustBelow20Percent_ReturnsMss()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyMSIStatus(0.19999);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSS),
            "0.19999 < 0.20 -> MSS (MSIsensor2 cutoff 'msi score >= 20%')");
    }

    // M7b — upper bound score 1.0 is a valid input and is MSI-H (>= 0.20).
    [Test]
    public void ClassifyMSIStatus_One_ReturnsMsiHigh()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyMSIStatus(1.0);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_High),
            "1.0 >= 0.20 -> MSI-H; score == 1.0 is a valid in-range input (INV-01 upper bound)");
    }

    // S4 — score outside [0,1] or non-finite throws (each guarded branch: >1, <0, NaN, +Inf, -Inf).
    [Test]
    public void ClassifyMSIStatus_InvalidScore_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyMSIStatus(1.5),
                "score > 1 is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyMSIStatus(-0.1),
                "score < 0 is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyMSIStatus(double.NaN),
                "NaN score is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyMSIStatus(double.PositiveInfinity),
                "+Infinity score is invalid (non-finite)");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyMSIStatus(double.NegativeInfinity),
                "-Infinity score is invalid (non-finite)");
        });
    }

    #endregion

    #region ClassifyBethesdaPanel

    // M9 — 0 of 5 markers unstable -> MSS (Boland 1998).
    [Test]
    public void ClassifyBethesdaPanel_ZeroOf5_ReturnsMss()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyBethesdaPanel(0, 5);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSS),
            "Boland 1998: no marker unstable -> MSS");
    }

    // M10 — exactly 1 of 5 markers unstable -> MSI-L (Boland 1998).
    [Test]
    public void ClassifyBethesdaPanel_OneOf5_ReturnsMsiLow()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyBethesdaPanel(1, 5);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_Low),
            "Boland 1998: only one of five markers unstable -> MSI-L");
    }

    // M11 — 2 of 5 markers unstable -> MSI-H (Boland 1998: >= 2/5).
    [Test]
    public void ClassifyBethesdaPanel_TwoOf5_ReturnsMsiHigh()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyBethesdaPanel(2, 5);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_High),
            "Boland 1998: two or more of five markers unstable -> MSI-H");
    }

    // M12 — 5 of 5 markers unstable -> MSI-H.
    [Test]
    public void ClassifyBethesdaPanel_AllOf5_ReturnsMsiHigh()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyBethesdaPanel(5, 5);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_High),
            "Boland 1998: all markers unstable -> MSI-H");
    }

    // C1 — larger panel, 4 of 10 unstable (>=2) -> MSI-H (Boland count rule).
    [Test]
    public void ClassifyBethesdaPanel_FourOf10_ReturnsMsiHigh()
    {
        OncologyAnalyzer.MsiStatus status = OncologyAnalyzer.ClassifyBethesdaPanel(4, 10);

        Assert.That(status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_High),
            "Boland 1998 count rule: >= 2 unstable markers -> MSI-H");
    }

    // S5 — invalid Bethesda inputs throw.
    [Test]
    public void ClassifyBethesdaPanel_InvalidInputs_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyBethesdaPanel(3, 2),
                "unstable > total is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyBethesdaPanel(-1, 5),
                "negative unstable is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyBethesdaPanel(1, 0),
                "total markers must be > 0");
        });
    }

    #endregion

    #region DetectMSI

    // M13 — 6 unstable / 20 valid = 0.30 >= 0.20 -> MSI-High end-to-end.
    [Test]
    public void DetectMSI_6Unstable20Total_ReturnsScore0_30AndMsiHigh()
    {
        bool[] flags = MakeFlags(unstable: 6, total: 20);

        OncologyAnalyzer.MsiResult result = OncologyAnalyzer.DetectMSI(flags);

        Assert.Multiple(() =>
        {
            Assert.That(result.UnstableLoci, Is.EqualTo(6), "6 unstable flags counted");
            Assert.That(result.TotalLoci, Is.EqualTo(20), "20 valid loci counted");
            Assert.That(result.Score, Is.EqualTo(0.30).Within(1e-10),
                "MSI score = 6/20 = 0.30 (MSIsensor2 fraction)");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSI_High),
                "0.30 >= 0.20 -> MSI-H (MSIsensor2 cutoff)");
        });
    }

    // M14 — 2 unstable / 20 valid = 0.10 < 0.20 -> MSS end-to-end.
    [Test]
    public void DetectMSI_2Unstable20Total_ReturnsScore0_10AndMss()
    {
        bool[] flags = MakeFlags(unstable: 2, total: 20);

        OncologyAnalyzer.MsiResult result = OncologyAnalyzer.DetectMSI(flags);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.EqualTo(0.10).Within(1e-10),
                "MSI score = 2/20 = 0.10 (MSIsensor2 fraction)");
            Assert.That(result.Status, Is.EqualTo(OncologyAnalyzer.MsiStatus.MSS),
                "0.10 < 0.20 -> MSS");
        });
    }

    // S6 — null flags throws.
    [Test]
    public void DetectMSI_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.DetectMSI(null!),
            "null per-locus flags is a guard violation");
    }

    // S7 — empty flags throws (no valid loci -> undefined score).
    [Test]
    public void DetectMSI_Empty_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.DetectMSI(Array.Empty<bool>()),
            "no valid loci -> MSI score undefined");
    }

    #endregion

    #region Sourced constants

    // Lock the source-backed thresholds so a silent constant change is caught.
    // MSIsensor2 README: "msi high: msi score >= 20%". Boland 1998: MSI-H >=2/5, MSI-L exactly 1/5.
    [Test]
    public void Constants_MatchSources()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.MsiHighScoreThreshold, Is.EqualTo(0.20).Within(1e-12),
                "MSIsensor2 recommended cutoff is 20% (>=)");
            Assert.That(OncologyAnalyzer.BethesdaMsiHighMarkerCount, Is.EqualTo(2),
                "Boland 1998: >=2 of 5 markers -> MSI-H");
            Assert.That(OncologyAnalyzer.BethesdaMsiLowMarkerCount, Is.EqualTo(1),
                "Boland 1998: exactly 1 of 5 markers -> MSI-L");
        });
    }

    #endregion

    #region Helpers

    private static bool[] MakeFlags(int unstable, int total)
    {
        var flags = new bool[total];
        for (int i = 0; i < unstable; i++)
        {
            flags[i] = true;
        }

        return flags;
    }

    #endregion
}
