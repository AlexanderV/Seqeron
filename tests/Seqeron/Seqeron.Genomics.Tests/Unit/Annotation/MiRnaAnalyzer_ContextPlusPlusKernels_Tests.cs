// MIRNA-CONTEXT-001 / MIRNA-TARGET-001 (08_DIFFERENTIAL strategy REF) — closed-form differential
// killers for the TargetScan context++ feature/contribution kernels (Agarwal et al. 2015). Each
// per-feature contribution is coeff × min-max-scaled(feature); the kernels are private and only
// observable through the full context++ score (which throws under the strict LimitationPolicy when
// any feature is omitted), so their scaling/indexing mutants survive. They are exposed internal
// (IVT) and asserted against values re-derived BY HAND from the published Agarwal_2015_parameters.txt
// coefficients (independent of the implementation) → mutated arithmetic/index/guard diverges → killed.

using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
[Category("MIRNA-CONTEXT-001")]
public class MiRnaAnalyzer_ContextPlusPlusKernels_Tests
{
    private const double Tol = 1e-9;
    private const TargetSiteType M8 = TargetSiteType.Seed8mer;

    // Published Agarwal_2015_parameters.txt values used as the independent oracle (8mer column).
    private static double Scaled(double raw, double coeff, double min, double max) => coeff * (raw - min) / (max - min);

    // ── ScaledContribution (generic getAgarwalContribution) ──
    [Test]
    public void ScaledContribution_Formula_IsCoeffTimesMinMaxScaled()
    {
        Assert.That(MiRnaAnalyzer.ScaledContribution(2.0, 0.5, 1.0, 3.0), Is.EqualTo(0.25).Within(Tol)); // 0.5*(2-1)/(3-1)
        Assert.That(MiRnaAnalyzer.ScaledContribution(7.0, -0.3, 4.0, 4.0), Is.EqualTo(-2.1).Within(Tol)); // max==min → coeff*raw
    }

    // ── caller-supplied scaled features (pure functions of the input) ──
    [Test]
    public void SpsContribution_8mer_ExactScaled()
        => Assert.That(MiRnaAnalyzer.SpsContribution(-8.0, M8), Is.EqualTo(Scaled(-8.0, 0.210, -11.13, -5.52)).Within(Tol));

    [Test]
    public void TaContribution_8mer_ExactScaled()
        => Assert.That(MiRnaAnalyzer.TaContribution(3.5, M8), Is.EqualTo(Scaled(3.5, 0.222, 3.113, 3.865)).Within(Tol));

    [Test]
    public void LenOrfContribution_8mer_Log10ThenScaled()
        => Assert.That(MiRnaAnalyzer.LenOrfContribution(1000.0, M8), Is.EqualTo(Scaled(3.0, 0.205, 2.788, 3.753)).Within(Tol)); // log10(1000)=3

    [Test]
    public void LenOrfContribution_NonPositive_Log10IsZero()
        => Assert.That(MiRnaAnalyzer.LenOrfContribution(0.0, M8), Is.EqualTo(Scaled(0.0, 0.205, 2.788, 3.753)).Within(Tol));

    [Test]
    public void Orf8mContribution_8mer_RawCount()
        => Assert.That(MiRnaAnalyzer.Orf8mContribution(3, M8), Is.EqualTo(-0.118 * 3).Within(Tol)); // used raw

    [Test]
    public void PctContribution_8mer_ExactScaled()
        => Assert.That(MiRnaAnalyzer.PctContribution(0.5, M8), Is.EqualTo(Scaled(0.5, -0.103, 0.0, 0.816)).Within(Tol));

    // ── features computed from the 3'UTR ──
    [Test]
    public void MinDistContribution_8mer_NearestEndLog10Scaled()
    {
        // length 100, site [10..17]: distTo5 = 10, distTo3 = 100-18 = 82, nearest = 10, log10 = 1.
        string mrna = new string('G', 100);
        Assert.That(MiRnaAnalyzer.MinDistContribution(mrna, 10, 17, M8),
            Is.EqualTo(Scaled(1.0, 0.118, 1.415, 3.113)).Within(Tol));
    }

    [Test]
    public void Len3UtrContribution_8mer_Log10LengthScaled()
    {
        // length 1000 → log10 = 3.
        Assert.That(MiRnaAnalyzer.Len3UtrContribution(new string('G', 1000), M8),
            Is.EqualTo(Scaled(3.0, 0.310, 2.392, 3.637)).Within(Tol));
    }

    [Test]
    public void Off6mContribution_8mer_CountsOffset6merRaw()
    {
        // Offset-6mer pattern = first 6 nt of revcomp(miRNA nt2-8). For miRNA "NAAAAAAAA..." the seed
        // region (nt2-8) = "AAAAAAA" → revcomp = "UUUUUUU" → pattern = "UUUUUU". Count its occurrences.
        var mirna = "CAAAAAAAAUUUUUUUUUUUU"; // nt2-8 = "AAAAAAA"
        string mrna = "UUUUUUUU";             // contains "UUUUUU" at offsets 0,1,2 → 3 occurrences
        int expectedCount = 3;
        Assert.That(MiRnaAnalyzer.Off6mContribution(mrna, mirna, M8),
            Is.EqualTo(-0.020 * expectedCount).Within(Tol));
    }

    // ── sRNA position indicators (binary, raw coeff) ──
    [Test]
    public void SRna1Contribution_8mer_SelectsByNt1()
    {
        Assert.That(MiRnaAnalyzer.SRna1Contribution("ACGUACGU", M8), Is.EqualTo(-0.018).Within(Tol), "nt1=A");
        Assert.That(MiRnaAnalyzer.SRna1Contribution("GCGUACGU", M8), Is.EqualTo(0.060).Within(Tol), "nt1=G");
        Assert.That(MiRnaAnalyzer.SRna1Contribution("UCGUACGU", M8), Is.EqualTo(0.0), "nt1=U → 0");
    }

    [Test]
    public void SRna8Contribution_8mer_SelectsByNt8()
    {
        Assert.That(MiRnaAnalyzer.SRna8Contribution("GGGGGGGG", M8), Is.EqualTo(0.015).Within(Tol), "nt8=G");
        Assert.That(MiRnaAnalyzer.SRna8Contribution("GGGGGGGA", M8), Is.EqualTo(0.022).Within(Tol), "nt8=A");
        Assert.That(MiRnaAnalyzer.SRna8Contribution("GGGGGGG", M8), Is.EqualTo(0.0), "length<8 → 0");
    }

    [Test]
    public void Site8Contribution_OnlyForA1And6mer()
    {
        // 8mer → always 0 (Site8 only defined for 7mer-A1 / 6mer).
        Assert.That(MiRnaAnalyzer.Site8Contribution("CGGGGGGG", 1, M8), Is.EqualTo(0.0), "8mer → 0");
        // 7mer-A1: base opposite pos8 = mrna[siteStart-1]; 'C' → CtxSite8C7merA1 = 0.036.
        Assert.That(MiRnaAnalyzer.Site8Contribution("CGGGGGGG", 1, TargetSiteType.Seed7merA1),
            Is.EqualTo(0.036).Within(Tol), "7mer-A1 pos8=C");
    }

    // ── Local AU positional weighting ──
    [Test]
    public void LocalAuContribution_8mer_PositionalWeighting()
    {
        // 8mer upstream weight = 1/(i+1). Site [2..9]; mrna[0..1] = "GA": upstream i=0→idx1='A'(w=1,AU),
        // i=1→idx0='G'(w=0.5,not AU). No downstream (site ends at last index). fraction = 1/1.5.
        string mrna = "GA" + new string('C', 8); // site occupies [2..9]
        double fraction = 1.0 / 1.5;
        Assert.That(MiRnaAnalyzer.LocalAuContribution(mrna, 2, 9, M8),
            Is.EqualTo(Scaled(fraction, -0.254, 0.308, 0.814)).Within(Tol));
    }

    // ── AlignMiRnaToTarget (public): exact pairing tally ──
    [Test]
    public void AlignMiRnaToTarget_TalliesMatchesWobblesMismatches()
    {
        // miRNA "ACGU" pairs antiparallel with target read 3'->5'. Target "ACGU": target[3..0]=U,G,C,A.
        // miRNA A:U match, C:G match, G:C match, U:A match → 4 matches, 0 mismatch.
        var d = AlignMiRnaToTarget("ACGU", "ACGU");
        Assert.That(d.Matches, Is.EqualTo(4));
        Assert.That(d.Mismatches, Is.EqualTo(0));
        Assert.That(d.AlignmentString, Is.EqualTo("||||"));

        // G:U wobble: miRNA "G" vs target "G" (target[0]=G read opposite) → G:U? target "U": miRNA G vs U = wobble.
        var w = AlignMiRnaToTarget("G", "U");
        Assert.That(w.GUWobbles, Is.EqualTo(1));
        Assert.That(w.Matches, Is.EqualTo(0));

        // Mismatch: miRNA "A" vs target "A" → A:A no pair.
        var mm = AlignMiRnaToTarget("A", "A");
        Assert.That(mm.Mismatches, Is.EqualTo(1));
    }
}
