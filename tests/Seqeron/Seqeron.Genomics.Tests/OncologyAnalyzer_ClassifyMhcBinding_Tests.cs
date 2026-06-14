// ONCO-MHC-001 — MHC-Peptide Binding Classification (length filtering + affinity/%rank thresholds)
// Evidence: docs/Evidence/ONCO-MHC-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-MHC-001.md
// Source: Reynisson B et al. (2020). NetMHCpan-4.1. Nucleic Acids Res 48(W1):W449-W454.
//         https://doi.org/10.1093/nar/gkaa379  (class I %Rank SB<0.5/WB<2; class II SB<2/WB<10; len 8-14, default 8-11)
//         Sette A et al. (1994). J Immunol 153(12):5586-92. PMID 7527444  (IC50 ~500 nM, preferably 50 nM)
//         IEDB threshold help: <50 nM high, <500 nM intermediate. IEDB class II tool desc: length 13-25.
//
// Expected categories below are derived from the cited cutoffs (strict "<"), NOT from running the code:
//   IC50 (nM): Strong < 50, Weak < 500, else NonBinder.
//   %Rank class I: Strong < 0.5, Weak < 2.  class II: Strong < 2, Weak < 10.
//   Length: class I 8-11 inclusive; class II 13-25 inclusive.

using System;
using MhcClass = Seqeron.Genomics.Oncology.OncologyAnalyzer.MhcClass;
using BindingStrength = Seqeron.Genomics.Oncology.OncologyAnalyzer.BindingStrength;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_ClassifyMhcBinding_Tests
{
    #region ClassifyBindingAffinity

    // M1 — IC50 10 nM < 50 ⇒ Strong (IEDB <50 high affinity; Sette 1994).
    [Test]
    public void ClassifyBindingAffinity_TenNm_ReturnsStrong()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingAffinity(10.0), Is.EqualTo(BindingStrength.Strong),
            "IC50 10 nM is below the 50 nM high-affinity cutoff ⇒ Strong (IEDB; Sette 1994).");
    }

    // M2 — IC50 exactly 50 nM: strict "<" ⇒ NOT strong ⇒ Weak (50 < 500).
    [Test]
    public void ClassifyBindingAffinity_FiftyNmBoundary_ReturnsWeak()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingAffinity(50.0), Is.EqualTo(BindingStrength.Weak),
            "IEDB states '<50 nM' (strict): 50 nM is not strong; 50 < 500 ⇒ Weak.");
    }

    // M3 — IC50 200 nM < 500 ⇒ Weak (intermediate affinity, IEDB).
    [Test]
    public void ClassifyBindingAffinity_TwoHundredNm_ReturnsWeak()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingAffinity(200.0), Is.EqualTo(BindingStrength.Weak),
            "IC50 200 nM is in [50,500) ⇒ Weak (IEDB '<500 nM intermediate affinity').");
    }

    // M4 — IC50 exactly 500 nM: strict "<" ⇒ NOT weak ⇒ NonBinder (Roomp 2010 demarcation).
    [Test]
    public void ClassifyBindingAffinity_FiveHundredNmBoundary_ReturnsNonBinder()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingAffinity(500.0), Is.EqualTo(BindingStrength.NonBinder),
            "IEDB '<500 nM' (strict): 500 nM is not a binder ⇒ NonBinder (Roomp 2010 binder demarcation).");
    }

    // M5 — IC50 1000 nM ≥ 500 ⇒ NonBinder.
    [Test]
    public void ClassifyBindingAffinity_ThousandNm_ReturnsNonBinder()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingAffinity(1000.0), Is.EqualTo(BindingStrength.NonBinder),
            "IC50 1000 nM is above the 500 nM weak cutoff ⇒ NonBinder.");
    }

    // S1 — IC50 must be > 0 (positive concentration; Registry invariant IC50 > 0).
    [Test]
    public void ClassifyBindingAffinity_ZeroOrNegative_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyBindingAffinity(0.0),
                "IC50 = 0 is not a positive concentration (invariant IC50 > 0).");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyBindingAffinity(-1.0),
                "IC50 < 0 is not a positive concentration.");
        });
    }

    // S2 — IC50 must be finite.
    [Test]
    public void ClassifyBindingAffinity_NonFinite_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyBindingAffinity(double.NaN),
                "NaN IC50 is invalid.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.ClassifyBindingAffinity(double.PositiveInfinity),
                "Infinite IC50 is invalid.");
        });
    }

    #endregion

    #region ClassifyBindingRank

    // M6 — class I %Rank 0.4 < 0.5 ⇒ Strong (Reynisson 2020).
    [Test]
    public void ClassifyBindingRank_ClassIPointFour_ReturnsStrong()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(0.4, MhcClass.ClassI), Is.EqualTo(BindingStrength.Strong),
            "Class I %Rank 0.4 < 0.5 ⇒ Strong (Reynisson 2020 '%Rank < 0.5% ... for SBs').");
    }

    // M7 — class I %Rank exactly 0.5: strict "<" ⇒ Weak (0.5 < 2).
    [Test]
    public void ClassifyBindingRank_ClassIHalfBoundary_ReturnsWeak()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(0.5, MhcClass.ClassI), Is.EqualTo(BindingStrength.Weak),
            "Reynisson '%Rank < 0.5%' (strict): 0.5 is not strong; 0.5 < 2 ⇒ Weak.");
    }

    // M8 — class I %Rank 1.0 < 2 ⇒ Weak.
    [Test]
    public void ClassifyBindingRank_ClassIOne_ReturnsWeak()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(1.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.Weak),
            "Class I %Rank 1.0 is in [0.5,2) ⇒ Weak (Reynisson 2020 '%Rank < 2% ... for WBs').");
    }

    // M9 — class I %Rank exactly 2.0: strict "<" ⇒ NonBinder.
    [Test]
    public void ClassifyBindingRank_ClassITwoBoundary_ReturnsNonBinder()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(2.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.NonBinder),
            "Reynisson '%Rank < 2%' (strict): 2.0 is not a weak binder ⇒ NonBinder.");
    }

    // M10 — class I %Rank 5.0 ⇒ NonBinder.
    [Test]
    public void ClassifyBindingRank_ClassIFive_ReturnsNonBinder()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(5.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.NonBinder),
            "Class I %Rank 5.0 ≥ 2 ⇒ NonBinder.");
    }

    // M11 — class II %Rank 1.5 < 2 ⇒ Strong (Reynisson 2020 class II SB < 2).
    [Test]
    public void ClassifyBindingRank_ClassIIOnePointFive_ReturnsStrong()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(1.5, MhcClass.ClassII), Is.EqualTo(BindingStrength.Strong),
            "Class II %Rank 1.5 < 2 ⇒ Strong (Reynisson 2020 '%Rank < 2% ... for SBs ... class II').");
    }

    // M12 — class II %Rank exactly 10.0: strict "<" ⇒ NonBinder.
    [Test]
    public void ClassifyBindingRank_ClassIITenBoundary_ReturnsNonBinder()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(10.0, MhcClass.ClassII), Is.EqualTo(BindingStrength.NonBinder),
            "Reynisson '%Rank < 10%' (strict): 10.0 is not a class II weak binder ⇒ NonBinder.");
    }

    // M13 — class II %Rank 5.0 in [2,10) ⇒ Weak.
    [Test]
    public void ClassifyBindingRank_ClassIIFive_ReturnsWeak()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(5.0, MhcClass.ClassII), Is.EqualTo(BindingStrength.Weak),
            "Class II %Rank 5.0 is in [2,10) ⇒ Weak (Reynisson 2020 class II WB < 10).");
    }

    // S3 — %Rank must lie in [0,100].
    [Test]
    public void ClassifyBindingRank_OutOfPercentileRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.ClassifyBindingRank(-0.1, MhcClass.ClassI),
                "%Rank < 0 is not a valid percentile.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.ClassifyBindingRank(100.1, MhcClass.ClassI),
                "%Rank > 100 is not a valid percentile.");
        });
    }

    // S4 — %Rank NaN rejected.
    [Test]
    public void ClassifyBindingRank_NaN_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.ClassifyBindingRank(double.NaN, MhcClass.ClassII),
            "NaN %Rank is invalid.");
    }

    // C1 — monotonicity of class I %Rank classification (INV-3): strength is non-increasing as %Rank rises.
    [Test]
    public void ClassifyBindingRank_ClassI_IsMonotoneInRank()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.ClassifyBindingRank(0.1, MhcClass.ClassI), Is.EqualTo(BindingStrength.Strong),
                "0.1 < 0.5 ⇒ Strong.");
            Assert.That(OncologyAnalyzer.ClassifyBindingRank(0.5, MhcClass.ClassI), Is.EqualTo(BindingStrength.Weak),
                "0.5 boundary ⇒ Weak (strict <).");
            Assert.That(OncologyAnalyzer.ClassifyBindingRank(1.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.Weak),
                "1.0 ⇒ Weak.");
            Assert.That(OncologyAnalyzer.ClassifyBindingRank(2.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.NonBinder),
                "2.0 boundary ⇒ NonBinder (strict <).");
            Assert.That(OncologyAnalyzer.ClassifyBindingRank(3.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.NonBinder),
                "3.0 ⇒ NonBinder.");
        });
    }

    #endregion

    #region IsValidPeptideLength

    // M14/M15/M16 — class I length: 9 valid; 7 too short; 12 above the 8-11 default.
    [Test]
    public void IsValidPeptideLength_ClassI_RespectsEightToEleven()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(8, MhcClass.ClassI), Is.True, "8 is the class I min.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(9, MhcClass.ClassI), Is.True, "9 is a valid class I length.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(11, MhcClass.ClassI), Is.True, "11 is the class I max (default).");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(7, MhcClass.ClassI), Is.False, "7 < 8 ⇒ too short.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(12, MhcClass.ClassI), Is.False,
                "12 > 11 ⇒ above the canonical class I default range 8-11.");
        });
    }

    // M17/M18/M19 — class II length: 15 valid; 12 too short; 26 too long (IEDB 13-25).
    [Test]
    public void IsValidPeptideLength_ClassII_RespectsThirteenToTwentyFive()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(13, MhcClass.ClassII), Is.True, "13 is the class II min.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(15, MhcClass.ClassII), Is.True, "15 is a valid class II length.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(25, MhcClass.ClassII), Is.True, "25 is the class II max.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(12, MhcClass.ClassII), Is.False, "12 < 13 ⇒ too short.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(26, MhcClass.ClassII), Is.False, "26 > 25 ⇒ too long.");
        });
    }

    // S5 — non-positive length is never valid.
    [Test]
    public void IsValidPeptideLength_NonPositive_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(0, MhcClass.ClassI), Is.False, "0 is not a valid length.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(-1, MhcClass.ClassII), Is.False, "negative length is invalid.");
        });
    }

    #endregion

    #region ClassifyMhcBinding

    // M20 — invalid length gates out: a 7-mer (class I) is not a candidate even with strong IC50.
    [Test]
    public void ClassifyMhcBinding_InvalidLength_ReturnsNonBinder()
    {
        Assert.That(OncologyAnalyzer.ClassifyMhcBinding(7, 10.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.NonBinder),
            "Length 7 is invalid for class I ⇒ NonBinder regardless of the strong 10 nM affinity.");
    }

    // M21 — valid length + strong affinity ⇒ Strong (delegates to ClassifyBindingAffinity).
    [Test]
    public void ClassifyMhcBinding_ValidLengthStrongAffinity_ReturnsStrong()
    {
        Assert.That(OncologyAnalyzer.ClassifyMhcBinding(9, 10.0, MhcClass.ClassI), Is.EqualTo(BindingStrength.Strong),
            "Length 9 valid (class I 8-11) and IC50 10 nM < 50 ⇒ Strong.");
    }

    #endregion
}
