// ONCO-MHC-001 — MHC-Peptide Binding (classification + matrix-based prediction)
// Evidence: docs/Evidence/ONCO-MHC-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-MHC-001.md
// Source: Reynisson B et al. (2020). NetMHCpan-4.1. Nucleic Acids Res 48(W1):W449-W454.
//         https://doi.org/10.1093/nar/gkaa379  (class I %Rank SB<0.5/WB<2; class II SB<2/WB<10; class I len 8-14)
//         Sette A et al. (1994). J Immunol 153(12):5586-92. PMID 7527444  (IC50 ~500 nM, preferably 50 nM)
//         IEDB threshold help: <50 nM high, <500 nM intermediate. IEDB class II tool desc: length 13-25.
//         Parker KC, Bednarek MA, Coligan JE (1994). J Immunol 152(1):163-175. PMID 8254189 + BIMAS scoring docs
//           (BIMAS: T1/2 = finalConstant * product-of-coefficients; unlisted residue coefficient = 1.0).
//         Peters B, Sette A (2005). BMC Bioinformatics 6:132. https://doi.org/10.1186/1471-2105-6-132 +
//           IEDB log50k = 1 - log(IC50)/log(50000)  =>  IC50 = 50000^(1 - score)  (SMM additive sum).
//
// Expected categories below are derived from the cited cutoffs (strict "<"), NOT from running the code:
//   IC50 (nM): Strong < 50, Weak < 500, else NonBinder.
//   %Rank class I: Strong < 0.5, Weak < 2.  class II: Strong < 2, Weak < 10.
//   Length: class I 8-14 inclusive (NetMHCpan-4.1 window); class II 13-25 inclusive.
// Prediction expected values are computed independently from the published rules:
//   SMM: IC50 = 50000^(1 - score)  =>  score 0 -> 50000, 0.5 -> sqrt(50000) = 223.6067977499790, 1 -> 1.
//   BIMAS: T1/2 = finalConstant * product of per-position coefficients (missing residue = 1.0).

using MhcClass = Seqeron.Genomics.Oncology.OncologyAnalyzer.MhcClass;
using BindingStrength = Seqeron.Genomics.Oncology.OncologyAnalyzer.BindingStrength;
using PmhcScoringMatrix = Seqeron.Genomics.Oncology.OncologyAnalyzer.PmhcScoringMatrix;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

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

    // M11b — class II %Rank exactly 2.0 (strong-cutoff boundary): strict "<" ⇒ NOT strong ⇒ Weak (2.0 < 10).
    [Test]
    public void ClassifyBindingRank_ClassIITwoBoundary_ReturnsWeak()
    {
        Assert.That(OncologyAnalyzer.ClassifyBindingRank(2.0, MhcClass.ClassII), Is.EqualTo(BindingStrength.Weak),
            "Reynisson 2020 class II '%Rank < 2%' (strict): 2.0 is not strong; 2.0 < 10 ⇒ Weak.");
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

    // M14/M15/M16 — class I length: 8/9/14 valid; 7 too short; 15 above the NetMHCpan-4.1 class I window 8-14.
    [Test]
    public void IsValidPeptideLength_ClassI_RespectsEightToFourteen()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(8, MhcClass.ClassI), Is.True, "8 is the class I min.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(9, MhcClass.ClassI), Is.True, "9 is a valid class I length.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(11, MhcClass.ClassI), Is.True, "11 is a valid class I length.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(14, MhcClass.ClassI), Is.True, "14 is the class I max (NetMHCpan-4.1 window).");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(7, MhcClass.ClassI), Is.False, "7 < 8 ⇒ too short.");
            Assert.That(OncologyAnalyzer.IsValidPeptideLength(15, MhcClass.ClassI), Is.False,
                "15 > 14 ⇒ above the NetMHCpan-4.1 class I peptide window 8-14.");
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

    // M21b — valid length but invalid IC50 (≤ 0) propagates ClassifyBindingAffinity's validation (INV-01).
    [Test]
    public void ClassifyMhcBinding_ValidLengthInvalidIc50_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.ClassifyMhcBinding(9, 0.0, MhcClass.ClassI),
            "Valid length passes the length gate, so the IC50 = 0 invariant (IC50 > 0) is enforced ⇒ throws.");
    }

    #endregion

    #region PredictIc50Smm (SMM transform: IC50 = 50000^(1 - score))

    // Builds a single-position SMM matrix whose only residue 'A' contributes `contribution`, with
    // the given intercept (FinalConstant). A 1-mer 'A' then has score = intercept + contribution.
    private static PmhcScoringMatrix SingleSmmPosition(double contribution, double intercept)
        => new(
            new IReadOnlyDictionary<char, double>[]
            {
                new Dictionary<char, double> { ['A'] = contribution },
            },
            intercept);

    // P1 — score 0 (intercept 0, contribution 0) ⇒ IC50 = 50000^1 = 50000 nM (IEDB log50k inverted).
    [Test]
    public void PredictIc50Smm_ScoreZero_ReturnsBase50000()
    {
        double ic50 = OncologyAnalyzer.PredictIc50Smm("A", SingleSmmPosition(0.0, 0.0));
        Assert.That(ic50, Is.EqualTo(50000.0).Within(1e-6),
            "IC50 = 50000^(1 - 0) = 50000 nM — the SMM transform's maximum (IEDB log50k = 1 - log(IC50)/log(50000)).");
    }

    // P2 — score 1 (contribution 1, intercept 0) ⇒ IC50 = 50000^0 = 1 nM.
    [Test]
    public void PredictIc50Smm_ScoreOne_ReturnsOneNm()
    {
        double ic50 = OncologyAnalyzer.PredictIc50Smm("A", SingleSmmPosition(1.0, 0.0));
        Assert.That(ic50, Is.EqualTo(1.0).Within(1e-9),
            "IC50 = 50000^(1 - 1) = 50000^0 = 1 nM (a strong binder by the SMM transform).");
    }

    // P3 — score 0.5 ⇒ IC50 = 50000^0.5 = sqrt(50000) = 223.6067977499790 nM (exact, independent of code).
    [Test]
    public void PredictIc50Smm_ScoreHalf_ReturnsSqrt50000()
    {
        double ic50 = OncologyAnalyzer.PredictIc50Smm("A", SingleSmmPosition(0.5, 0.0));
        Assert.That(ic50, Is.EqualTo(223.6067977499790).Within(1e-9),
            "IC50 = 50000^(1 - 0.5) = sqrt(50000) = 223.6067977499790 nM (IEDB transform).");
    }

    // P4 — the intercept (FinalConstant) is part of the score: intercept 0.3 + contribution 0.2 = 0.5
    // ⇒ same sqrt(50000), proving the intercept is added (not ignored).
    [Test]
    public void PredictIc50Smm_InterceptAddedToScore_ReturnsSqrt50000()
    {
        double ic50 = OncologyAnalyzer.PredictIc50Smm("A", SingleSmmPosition(0.2, 0.3));
        Assert.That(ic50, Is.EqualTo(223.6067977499790).Within(1e-9),
            "score = intercept 0.3 + contribution 0.2 = 0.5 ⇒ IC50 = sqrt(50000); the intercept must be summed in.");
    }

    // P5 — an unlisted residue contributes 0 (additive identity): peptide 'C' (not in the row) ⇒ score 0
    // ⇒ IC50 = 50000 nM. A wrong default (e.g. treating missing as the listed value) would fail.
    [Test]
    public void PredictIc50Smm_UnlistedResidue_ContributesZero()
    {
        double ic50 = OncologyAnalyzer.PredictIc50Smm("C", SingleSmmPosition(1.0, 0.0));
        Assert.That(ic50, Is.EqualTo(50000.0).Within(1e-6),
            "Residue 'C' is absent from the row ⇒ contributes 0 (SMM additive identity) ⇒ score 0 ⇒ IC50 = 50000 nM.");
    }

    // P5b — multi-position additive sum + intercept + a missing residue, hand-computed end-to-end.
    // 3-position matrix: intercept 0.05; pos0 K=0.30; pos1 lists only V (peptide's 'A' is unlisted ⇒ 0);
    // pos2 Y=0.25. Peptide "KAY" ⇒ score = 0.05 + 0.30 + 0 + 0.25 = 0.60 ⇒ IC50 = 50000^(1-0.60) =
    // 75.78582832551992 nM. Exercises the SMM additive position-specific sum (not just a single position),
    // the intercept being summed in, AND the missing-residue additive identity together.
    [Test]
    public void PredictIc50Smm_MultiPositionAdditiveSumWithMissingResidue_HandComputed()
    {
        var matrix = new PmhcScoringMatrix(
            new IReadOnlyDictionary<char, double>[]
            {
                new Dictionary<char, double> { ['K'] = 0.30 },
                new Dictionary<char, double> { ['V'] = 0.20 }, // peptide has 'A' here ⇒ contributes 0
                new Dictionary<char, double> { ['Y'] = 0.25 },
            },
            0.05);
        double ic50 = OncologyAnalyzer.PredictIc50Smm("KAY", matrix);
        Assert.That(ic50, Is.EqualTo(75.78582832551992).Within(1e-9),
            "score = intercept 0.05 + K 0.30 + (A unlisted ⇒ 0) + Y 0.25 = 0.60 ⇒ IC50 = 50000^(1-0.60) = 75.78582832551992 nM.");
    }

    // P5c — a non-amino-acid character in the peptide is simply unlisted at its position ⇒ contributes 0
    // (the SMM additive identity), exactly like any residue absent from the row. Single position lists only
    // 'A'; peptide "1" (a digit) ⇒ score 0 ⇒ IC50 = 50000 nM. No throw — non-AA is treated as neutral.
    [Test]
    public void PredictIc50Smm_NonAminoAcidCharacter_ContributesZero()
    {
        double ic50 = OncologyAnalyzer.PredictIc50Smm("1", SingleSmmPosition(1.0, 0.0));
        Assert.That(ic50, Is.EqualTo(50000.0).Within(1e-6),
            "Non-AA char '1' is absent from the row ⇒ contributes 0 (additive identity) ⇒ score 0 ⇒ IC50 = 50000 nM.");
    }

    // P6 — null peptide ⇒ ArgumentNullException.
    [Test]
    public void PredictIc50Smm_NullPeptide_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.PredictIc50Smm(null!, SingleSmmPosition(1.0, 0.0)),
            "A null peptide is invalid input.");
    }

    // P7 — peptide length ≠ matrix row count ⇒ ArgumentException (one row, 2-residue peptide).
    [Test]
    public void PredictIc50Smm_LengthMismatch_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.PredictIc50Smm("AA", SingleSmmPosition(1.0, 0.0)),
            "Peptide length 2 ≠ matrix position count 1 ⇒ ArgumentException.");
    }

    // P8 — empty matrix (no rows) ⇒ ArgumentException.
    [Test]
    public void PredictIc50Smm_EmptyMatrix_Throws()
    {
        var empty = new PmhcScoringMatrix(Array.Empty<IReadOnlyDictionary<char, double>>(), 0.0);
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.PredictIc50Smm("", empty),
            "A matrix with no position rows cannot score any peptide ⇒ ArgumentException.");
    }

    #endregion

    #region PredictAndClassifySmm (predict → classify chain)

    // The 9-mer SMM matrix for GILGFVFTL (influenza M1 58-66, the paradigm HLA-A*02:01 binder):
    // per-position contributions 0.20,0.10,0.15,0.05,0.10,0.10,0.10,0.05,0.15 sum to exactly 1.0
    // (intercept 0) ⇒ score 1 ⇒ IC50 = 1 nM. Every other residue at each position contributes 0,
    // so a poly-residue peptide that hits none of the listed residues scores 0 ⇒ IC50 = 50000 nM.
    private static PmhcScoringMatrix StrongBinderSmmMatrix()
    {
        const string peptide = "GILGFVFTL";
        double[] contrib = { 0.20, 0.10, 0.15, 0.05, 0.10, 0.10, 0.10, 0.05, 0.15 };
        var rows = new IReadOnlyDictionary<char, double>[peptide.Length];
        for (int i = 0; i < peptide.Length; i++)
        {
            rows[i] = new Dictionary<char, double> { [peptide[i]] = contrib[i] };
        }

        return new PmhcScoringMatrix(rows, 0.0);
    }

    // P9 — known strong binder GILGFVFTL: score 1.0 ⇒ IC50 = 1 nM ⇒ Strong (end-to-end predict→classify).
    [Test]
    public void PredictAndClassifySmm_KnownStrongBinder_ReturnsStrong()
    {
        var (ic50, strength) = OncologyAnalyzer.PredictAndClassifySmm("GILGFVFTL", StrongBinderSmmMatrix());
        Assert.Multiple(() =>
        {
            Assert.That(ic50, Is.EqualTo(1.0).Within(1e-9),
                "GILGFVFTL contributions sum to 1.0 ⇒ IC50 = 50000^0 = 1 nM.");
            Assert.That(strength, Is.EqualTo(BindingStrength.Strong),
                "IC50 = 1 nM < 50 nM ⇒ Strong (chains into ClassifyBindingAffinity).");
        });
    }

    // P10 — non-binder (no listed residues matched): poly-'W' 9-mer ⇒ score 0 ⇒ IC50 = 50000 nM ⇒ NonBinder,
    // and the strong binder's IC50 (1 nM) is far below the non-binder's (50000 nM) — the required ranking.
    [Test]
    public void PredictAndClassifySmm_NonBinderRankedAboveStrong()
    {
        var matrix = StrongBinderSmmMatrix();
        var (strongIc50, strongStrength) = OncologyAnalyzer.PredictAndClassifySmm("GILGFVFTL", matrix);
        var (nonIc50, nonStrength) = OncologyAnalyzer.PredictAndClassifySmm("WWWWWWWWW", matrix);
        Assert.Multiple(() =>
        {
            Assert.That(nonIc50, Is.EqualTo(50000.0).Within(1e-6),
                "Poly-W matches no listed residue ⇒ score 0 ⇒ IC50 = 50000 nM.");
            Assert.That(nonStrength, Is.EqualTo(BindingStrength.NonBinder),
                "IC50 = 50000 nM ≥ 500 nM ⇒ NonBinder.");
            Assert.That(strongIc50, Is.LessThan(nonIc50),
                "The strong binder's predicted IC50 (1 nM) must be far below the non-binder's (50000 nM).");
            Assert.That(strongStrength, Is.EqualTo(BindingStrength.Strong),
                "Sanity: the binder is Strong while the non-binder is NonBinder (ranking holds).");
        });
    }

    // P10b — Weak band of the predict→classify chain: a 2-position matrix with intercept 0.1, pos0 L=0.3,
    // pos1 V=0.2 ⇒ score 0.6 ⇒ IC50 = 50000^(1-0.6) = 75.78582832551992 nM, which is in [50,500) ⇒ Weak.
    // Covers the middle classification band that P9 (Strong) and P10 (NonBinder) leave untested.
    [Test]
    public void PredictAndClassifySmm_WeakBand_ReturnsWeak()
    {
        var matrix = new PmhcScoringMatrix(
            new IReadOnlyDictionary<char, double>[]
            {
                new Dictionary<char, double> { ['L'] = 0.3 },
                new Dictionary<char, double> { ['V'] = 0.2 },
            },
            0.1);
        var (ic50, strength) = OncologyAnalyzer.PredictAndClassifySmm("LV", matrix);
        Assert.Multiple(() =>
        {
            Assert.That(ic50, Is.EqualTo(75.78582832551992).Within(1e-9),
                "score = 0.1 + 0.3 + 0.2 = 0.6 ⇒ IC50 = 50000^(1-0.6) = 75.78582832551992 nM (hand-computed).");
            Assert.That(strength, Is.EqualTo(BindingStrength.Weak),
                "IC50 75.79 nM ∈ [50,500) ⇒ Weak (chains into ClassifyBindingAffinity's middle band).");
        });
    }

    #endregion

    #region PredictBindingHalfLifeBimas (BIMAS product rule)

    // P11 — BIMAS product: T1/2 = finalConstant * (2.0 * 3.0 * 1.5) = 10 * 9 = 90.0 (exact, hand-computed).
    [Test]
    public void PredictBindingHalfLifeBimas_ProductTimesConstant_Returns90()
    {
        var matrix = new PmhcScoringMatrix(
            new IReadOnlyDictionary<char, double>[]
            {
                new Dictionary<char, double> { ['L'] = 2.0 },
                new Dictionary<char, double> { ['M'] = 3.0 },
                new Dictionary<char, double> { ['V'] = 1.5 },
            },
            10.0);
        double t12 = OncologyAnalyzer.PredictBindingHalfLifeBimas("LMV", matrix);
        Assert.That(t12, Is.EqualTo(90.0).Within(1e-9),
            "BIMAS: T1/2 = finalConstant 10 * (2.0*3.0*1.5) = 10*9 = 90.0 (running score starts at 1.0).");
    }

    // P12 — unlisted residues contribute the neutral coefficient 1.0: 'AAA' on the same matrix ⇒
    // T1/2 = 10 * (1*1*1) = 10.0. A wrong default (e.g. 0) would give 0 and fail.
    [Test]
    public void PredictBindingHalfLifeBimas_UnlistedResidues_NeutralOne()
    {
        var matrix = new PmhcScoringMatrix(
            new IReadOnlyDictionary<char, double>[]
            {
                new Dictionary<char, double> { ['L'] = 2.0 },
                new Dictionary<char, double> { ['M'] = 3.0 },
                new Dictionary<char, double> { ['V'] = 1.5 },
            },
            10.0);
        double t12 = OncologyAnalyzer.PredictBindingHalfLifeBimas("AAA", matrix);
        Assert.That(t12, Is.EqualTo(10.0).Within(1e-9),
            "Unlisted residues use coefficient 1.0 (BIMAS) ⇒ T1/2 = 10 * 1 * 1 * 1 = 10.0.");
    }

    // P13 — a strong-anchor peptide outscores a weak one (product rule preserves ranking).
    [Test]
    public void PredictBindingHalfLifeBimas_StrongAnchorOutscoresWeak()
    {
        var matrix = new PmhcScoringMatrix(
            new IReadOnlyDictionary<char, double>[]
            {
                new Dictionary<char, double> { ['L'] = 5.0, ['A'] = 0.1 },
                new Dictionary<char, double> { ['V'] = 4.0, ['A'] = 0.2 },
            },
            1.0);
        double strong = OncologyAnalyzer.PredictBindingHalfLifeBimas("LV", matrix); // 1*5*4 = 20
        double weak = OncologyAnalyzer.PredictBindingHalfLifeBimas("AA", matrix);   // 1*0.1*0.2 = 0.02
        Assert.Multiple(() =>
        {
            Assert.That(strong, Is.EqualTo(20.0).Within(1e-9), "Favorable anchors: 1*5*4 = 20.0.");
            Assert.That(weak, Is.EqualTo(0.02).Within(1e-9), "Unfavorable anchors: 1*0.1*0.2 = 0.02.");
            Assert.That(strong, Is.GreaterThan(weak), "Higher BIMAS T1/2 ⇒ stronger predicted binder.");
        });
    }

    #endregion

    #region LoadScoringMatrix (caller-supplied matrix loader)

    // P14 — loader parses CONST + per-position RESIDUE=VALUE rows, handles comments/blanks, upper-cases
    // residues, and round-trips into a correct BIMAS prediction (10 * 2.0 * 3.0 = 60.0).
    [Test]
    public void LoadScoringMatrix_ParsesConstAndRows_RoundTripsBimas()
    {
        string[] lines =
        {
            "# example BIMAS-style matrix (caller-supplied values)",
            "CONST=10.0",
            "l=2.0 M=1.0",     // lowercase 'l' must upper-case to 'L'
            "",                 // blank line skipped
            "M=3.0, A=0.5",     // comma-separated tokens
        };
        PmhcScoringMatrix matrix = OncologyAnalyzer.LoadScoringMatrix(lines);
        Assert.Multiple(() =>
        {
            Assert.That(matrix.Rows, Has.Count.EqualTo(2), "Two position rows (CONST line is not a position).");
            Assert.That(matrix.FinalConstant, Is.EqualTo(10.0).Within(1e-12), "CONST=10.0 sets the final constant.");
            double t12 = OncologyAnalyzer.PredictBindingHalfLifeBimas("LM", matrix);
            Assert.That(t12, Is.EqualTo(60.0).Within(1e-9),
                "'L' at pos0 = 2.0 (lower-cased input), 'M' at pos1 = 3.0 ⇒ 10*2*3 = 60.0.");
        });
    }

    // P14b — when no CONST line is given, FinalConstant defaults to the multiplicative/additive identity 1.0,
    // so a BIMAS product round-trips unchanged: "LM" on rows L=2.0 / M=3.0 ⇒ 1.0 * 2.0 * 3.0 = 6.0.
    [Test]
    public void LoadScoringMatrix_NoConst_DefaultsToIdentityOne()
    {
        PmhcScoringMatrix matrix = OncologyAnalyzer.LoadScoringMatrix(new[] { "L=2.0", "M=3.0" });
        Assert.Multiple(() =>
        {
            Assert.That(matrix.FinalConstant, Is.EqualTo(1.0).Within(1e-12),
                "No CONST line ⇒ FinalConstant defaults to the identity 1.0.");
            Assert.That(matrix.Rows, Has.Count.EqualTo(2), "Two position rows.");
            double t12 = OncologyAnalyzer.PredictBindingHalfLifeBimas("LM", matrix);
            Assert.That(t12, Is.EqualTo(6.0).Within(1e-9),
                "Identity constant ⇒ T1/2 = 1.0 * 2.0 * 3.0 = 6.0 (constant does not perturb the product).");
        });
    }

    // P15 — malformed token (no '=') ⇒ FormatException.
    [Test]
    public void LoadScoringMatrix_MalformedToken_Throws()
    {
        Assert.Throws<FormatException>(
            () => OncologyAnalyzer.LoadScoringMatrix(new[] { "L 2.0" }),
            "A token without '=' is malformed ⇒ FormatException.");
    }

    // P16 — non-numeric value ⇒ FormatException.
    [Test]
    public void LoadScoringMatrix_NonNumericValue_Throws()
    {
        Assert.Throws<FormatException>(
            () => OncologyAnalyzer.LoadScoringMatrix(new[] { "L=abc" }),
            "A non-numeric value ⇒ FormatException.");
    }

    // P17 — multi-character residue key ⇒ FormatException.
    [Test]
    public void LoadScoringMatrix_MultiCharResidue_Throws()
    {
        Assert.Throws<FormatException>(
            () => OncologyAnalyzer.LoadScoringMatrix(new[] { "LM=2.0" }),
            "A residue key must be a single amino-acid letter ⇒ FormatException.");
    }

    // P18 — null input ⇒ ArgumentNullException.
    [Test]
    public void LoadScoringMatrix_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.LoadScoringMatrix(null!),
            "Null line enumerable is invalid input.");
    }

    #endregion
}
