using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ONCO-* mutation killers (batch 1): exact-value, inclusive-boundary and exception tests for the
/// closed-form clinical-genomics scoring methods of <see cref="OncologyAnalyzer"/> that the existing
/// fixtures left wholly or partly uncovered. Each assertion reproduces the published formula/cutoff
/// so an injected operator/threshold change diverges.
///
/// Evidence: CCF point estimate (McGranahan 2016; Tarabichi 2021); binomial-grid clonality (Landau 2013);
/// TMB ≥10 mut/Mb (Marcus 2021); MSI 20% / Bethesda ≥2-of-5 (msisensor2; Boland 1998); HRD ≥42 (Telli 2016);
/// SBS-96 pyrimidine folding (Alexandrov 2013; COSMIC); MHC IC50 50/500 nM &amp; NetMHCpan %Rank (Reynisson 2020);
/// ctDNA Poisson p = 1−e^(−ndk) (Avanzini 2020).
/// </summary>
[TestFixture]
public class OncologyAnalyzer_MutationKillers_Tests
{
    private const double Tol = 1e-9;

    #region VAF

    [Test]
    public void CalculateVAF_Fraction_AndZeroCoverage()
    {
        Assert.That(CalculateVAF(30, 100), Is.EqualTo(0.3).Within(Tol));
        Assert.That(CalculateVAF(0, 0), Is.EqualTo(0.0).Within(Tol)); // uncovered site ⇒ 0
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateVAF(101, 100));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateVAF(-1, 100));
    }

    #endregion

    #region EstimateCcf

    [Test]
    public void EstimateCcf_ClosedFormAndCap()
    {
        // CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m).
        var a = EstimateCcf(vaf: 0.1, purity: 0.8, tumorCopyNumber: 4, multiplicity: 2);
        Assert.That(a.RawCcf, Is.EqualTo(0.1 * (0.8 * 4 + 2 * 0.2) / (0.8 * 2)).Within(Tol)); // 0.225
        Assert.That(a.Ccf, Is.EqualTo(0.225).Within(Tol));

        // Over-unity raw value is capped to 1 but exposed uncapped.
        var b = EstimateCcf(vaf: 0.5, purity: 0.5, tumorCopyNumber: 2, multiplicity: 1);
        Assert.That(b.RawCcf, Is.EqualTo(2.0).Within(Tol));
        Assert.That(b.Ccf, Is.EqualTo(1.0).Within(Tol)); // kills the Math.Min cap
    }

    [Test]
    public void EstimateCcf_DomainGuards()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => EstimateCcf(1.1, 0.5, 2, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => EstimateCcf(0.5, 0.0, 2, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => EstimateCcf(0.5, 0.5, 0, 1));
        Assert.Throws<ArgumentException>(() => EstimateCcf(0.5, 0.5, 2, 3));
    }

    #endregion

    #region IdentifyClonalMutations & ClassifyClonality

    [Test]
    public void IdentifyClonalMutations_StrictAbove095()
    {
        var idx = IdentifyClonalMutations(new[] { 0.96, 0.5, 1.0, 0.95 });
        Assert.That(idx, Is.EqualTo(new[] { 0, 2 })); // 0.95 is NOT > 0.95
    }

    [Test]
    public void ClassifyClonality_HighVafClonal_LowVafSubclonal()
    {
        var clonal = ClassifyClonality(new[] { new ClonalityVariant(100, 100, 2, 2) }, purity: 1.0);
        Assert.That(clonal.Calls[0].Status, Is.EqualTo(ClonalityStatus.Clonal));
        Assert.That(clonal.ClonalCount, Is.EqualTo(1));
        Assert.That(clonal.SubclonalCount, Is.EqualTo(0));
        Assert.That(clonal.ClonalFraction, Is.EqualTo(1.0).Within(Tol));

        var sub = ClassifyClonality(new[] { new ClonalityVariant(20, 100, 2, 2) }, purity: 1.0);
        Assert.That(sub.Calls[0].Status, Is.EqualTo(ClonalityStatus.Subclonal));
        Assert.That(sub.ClonalCount, Is.EqualTo(0));
        Assert.That(sub.ClonalFraction, Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void ClassifyClonality_EmptyAndGuards()
    {
        var empty = ClassifyClonality(Array.Empty<ClonalityVariant>(), 0.5);
        Assert.That(empty.Calls, Is.Empty);
        Assert.That(empty.ClonalFraction, Is.EqualTo(0.0).Within(Tol));

        Assert.Throws<ArgumentOutOfRangeException>(() => ClassifyClonality(new[] { new ClonalityVariant(10, 100, 2) }, 0.0));
        Assert.Throws<ArgumentException>(() => ClassifyClonality(new[] { new ClonalityVariant(101, 100, 2) }, 0.5));
    }

    #endregion

    #region TMB

    [Test]
    public void CalculateTMB_PerMegabaseAndGuards()
    {
        Assert.That(CalculateTMB(100, 2.0), Is.EqualTo(50.0).Within(Tol));
        Assert.That(CalculateTMB(0, 1.0), Is.EqualTo(0.0).Within(Tol));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateTMB(-1, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateTMB(10, 0.0));
    }

    [Test]
    public void ClassifyTMB_InclusiveHighCutoff()
    {
        Assert.That(ClassifyTMB(10.0), Is.EqualTo(TmbStatus.High)); // ≥10 inclusive
        Assert.That(ClassifyTMB(9.9), Is.EqualTo(TmbStatus.Low));
        Assert.Throws<ArgumentOutOfRangeException>(() => ClassifyTMB(-0.1));
    }

    #endregion

    #region MSI

    [Test]
    public void Msi_ScoreClassificationAndBethesda()
    {
        Assert.That(CalculateMSIScore(1, 5), Is.EqualTo(0.2).Within(Tol));
        Assert.That(ClassifyMSIStatus(0.2), Is.EqualTo(MsiStatus.MSI_High)); // ≥0.20 inclusive
        Assert.That(ClassifyMSIStatus(0.19), Is.EqualTo(MsiStatus.MSS));

        Assert.That(ClassifyBethesdaPanel(2, 5), Is.EqualTo(MsiStatus.MSI_High)); // ≥2
        Assert.That(ClassifyBethesdaPanel(1, 5), Is.EqualTo(MsiStatus.MSI_Low));  // exactly 1
        Assert.That(ClassifyBethesdaPanel(0, 5), Is.EqualTo(MsiStatus.MSS));

        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMSIScore(1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMSIScore(6, 5));
    }

    [Test]
    public void DetectMSI_EndToEnd()
    {
        var r = DetectMSI(new[] { true, false, false, false, false });
        Assert.That(r.UnstableLoci, Is.EqualTo(1));
        Assert.That(r.TotalLoci, Is.EqualTo(5));
        Assert.That(r.Score, Is.EqualTo(0.2).Within(Tol));
        Assert.That(r.Status, Is.EqualTo(MsiStatus.MSI_High));
        Assert.Throws<ArgumentOutOfRangeException>(() => DetectMSI(Array.Empty<bool>()));
    }

    #endregion

    #region HRD

    [Test]
    public void Hrd_ScoreSumAndInclusiveCutoff()
    {
        Assert.That(CalculateHRDScore(20, 12, 10), Is.EqualTo(42));
        Assert.That(CalculateHRDScore(1, 1, 1), Is.EqualTo(3));
        Assert.That(ClassifyHRDStatus(42), Is.EqualTo(HrdStatus.HrdHigh));   // ≥42 inclusive
        Assert.That(ClassifyHRDStatus(41), Is.EqualTo(HrdStatus.HrdNegative));

        var d = DetectHRD(new HrdComponents(20, 12, 10));
        Assert.That(d.Score, Is.EqualTo(42));
        Assert.That(d.Status, Is.EqualTo(HrdStatus.HrdHigh));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateHRDScore(-1, 0, 0));
    }

    #endregion

    #region SBS-96 context

    [Test]
    public void ClassifySbsContext_PyrimidineKept_PurineFolded()
    {
        // Pyrimidine reference (C) is kept on its strand.
        Assert.That(ClassifySbsContext('A', 'C', 'T', 'G'), Is.EqualTo("A[C>T]G"));
        // Purine reference (G) is reverse-complement-folded onto the pyrimidine strand.
        Assert.That(ClassifySbsContext('A', 'G', 'A', 'T'), Is.EqualTo("A[C>T]T"));
        Assert.Throws<ArgumentException>(() => ClassifySbsContext('A', 'C', 'C', 'G'));
        Assert.Throws<ArgumentException>(() => ClassifySbsContext('A', 'X', 'T', 'G'));
    }

    [Test]
    public void Sbs96_ChannelsAndCatalog()
    {
        var channels = EnumerateSbs96Channels();
        Assert.That(channels, Has.Count.EqualTo(96));
        Assert.That(channels[0], Is.EqualTo("A[C>A]A")); // substitution-major, 5'-then-3' order
        Assert.That(channels, Does.Contain("T[T>G]T"));

        var catalog = Build96ContextCatalog(new[] { ('A', 'C', 'T', 'G') });
        Assert.That(catalog.Count, Is.EqualTo(96));
        Assert.That(catalog["A[C>T]G"], Is.EqualTo(1));
        Assert.That(catalog.Values.Sum(), Is.EqualTo(1));
    }

    #endregion

    #region MHC binding

    [Test]
    public void ClassifyBindingAffinity_Ic50Cutoffs()
    {
        Assert.That(ClassifyBindingAffinity(49.0), Is.EqualTo(BindingStrength.Strong));   // <50
        Assert.That(ClassifyBindingAffinity(50.0), Is.EqualTo(BindingStrength.Weak));     // [50,500)
        Assert.That(ClassifyBindingAffinity(499.0), Is.EqualTo(BindingStrength.Weak));
        Assert.That(ClassifyBindingAffinity(500.0), Is.EqualTo(BindingStrength.NonBinder)); // ≥500
        Assert.Throws<ArgumentOutOfRangeException>(() => ClassifyBindingAffinity(0.0));
    }

    [Test]
    public void ClassifyBindingRank_ClassSpecificCutoffs()
    {
        Assert.That(ClassifyBindingRank(0.4, MhcClass.ClassI), Is.EqualTo(BindingStrength.Strong));    // <0.5
        Assert.That(ClassifyBindingRank(0.5, MhcClass.ClassI), Is.EqualTo(BindingStrength.Weak));      // [0.5,2)
        Assert.That(ClassifyBindingRank(2.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.NonBinder)); // ≥2
        Assert.That(ClassifyBindingRank(1.5, MhcClass.ClassII), Is.EqualTo(BindingStrength.Strong));   // <2
        Assert.That(ClassifyBindingRank(9.0, MhcClass.ClassII), Is.EqualTo(BindingStrength.Weak));     // [2,10)
        Assert.That(ClassifyBindingRank(10.0, MhcClass.ClassII), Is.EqualTo(BindingStrength.NonBinder));
    }

    [Test]
    public void IsValidPeptideLength_AndMhcBindingLengthGate()
    {
        Assert.That(IsValidPeptideLength(8, MhcClass.ClassI), Is.True);
        Assert.That(IsValidPeptideLength(11, MhcClass.ClassI), Is.True);
        Assert.That(IsValidPeptideLength(12, MhcClass.ClassI), Is.False);
        Assert.That(IsValidPeptideLength(13, MhcClass.ClassII), Is.True);

        Assert.That(ClassifyMhcBinding(9, 40.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.Strong));
        Assert.That(ClassifyMhcBinding(7, 40.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.NonBinder)); // bad length
    }

    #endregion

    #region ctDNA Poisson

    [Test]
    public void CtDna_DetectionProbabilityAndExpectedMolecules()
    {
        // λ = n·d·k; p = 1 − e^(−λ).
        Assert.That(CtDnaDetectionProbability(15000, 0.001, 1), Is.EqualTo(1.0 - Math.Exp(-15.0)).Within(Tol));
        Assert.That(CtDnaDetectionProbability(0, 0.5), Is.EqualTo(0.0).Within(Tol)); // λ=0 ⇒ p=0
        Assert.That(ExpectedMutantMolecules(15000, 0.001, 1), Is.EqualTo(15.0).Within(Tol));
        Assert.That(ExpectedMutantMolecules(1000, 0.01, 3), Is.EqualTo(30.0).Within(Tol));
        Assert.Throws<ArgumentOutOfRangeException>(() => CtDnaDetectionProbability(-1, 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExpectedMutantMolecules(100, 1.5));
    }

    #endregion
}
