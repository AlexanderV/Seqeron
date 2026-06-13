// SEQ-SECSTRUCT-001 — Protein Secondary Structure Prediction (Chou-Fasman propensity profile)
// Evidence: docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-SECSTRUCT-001.md
// Source: Chou PY, Fasman GD (1978). Annu Rev Biochem 47:251-276.
//         Pa/Pb/Pt verbatim from Przytycka NCBI lecture + ravihansa3000/ChouFasman ref impl.
//         Window mean = arithmetic mean of member residue propensities (per component).

using System.Linq;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceStatistics_PredictSecondaryStructure_Tests
{
    // Expected values are exact derivations from the Chou-Fasman propensity table.
    private const double Tolerance = 1e-10;

    #region PredictSecondaryStructure

    // M1 — single residue A, window 1: returns A's (Pa, Pb, Pt) = (1.42, 0.83, 0.66).
    // Evidence: A 1.42/0.83/0.66 (Przytycka, ref impl); mean of one value is that value (INV-01).
    [Test]
    public void PredictSecondaryStructure_SingleResidueAla_ReturnsAlaPropensities()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("A", windowSize: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Helix, Is.EqualTo(1.42).Within(Tolerance), "Ala Pa = 1.42 (INV-01)");
            Assert.That(result.Sheet, Is.EqualTo(0.83).Within(Tolerance), "Ala Pb = 0.83 (INV-01)");
            Assert.That(result.Turn, Is.EqualTo(0.66).Within(Tolerance), "Ala Pt = 0.66 (INV-01)");
        });
    }

    // M2 — single residue E, window 1: (1.51, 0.37, 0.74); strongest helix former.
    // Evidence: E 1.51/0.37/0.74 (Przytycka, ref impl).
    [Test]
    public void PredictSecondaryStructure_SingleResidueGlu_ReturnsGluPropensities()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("E", windowSize: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Helix, Is.EqualTo(1.51).Within(Tolerance), "Glu Pa = 1.51");
            Assert.That(result.Sheet, Is.EqualTo(0.37).Within(Tolerance), "Glu Pb = 0.37");
            Assert.That(result.Turn, Is.EqualTo(0.74).Within(Tolerance), "Glu Pt = 0.74");
        });
    }

    // M3 — single residue V, window 1: (1.06, 1.70, 0.50); strongest sheet former.
    // Evidence: V 1.06/1.70/0.50 (CSB|SJU, ref impl).
    [Test]
    public void PredictSecondaryStructure_SingleResidueVal_ReturnsValPropensities()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("V", windowSize: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Helix, Is.EqualTo(1.06).Within(Tolerance), "Val Pa = 1.06");
            Assert.That(result.Sheet, Is.EqualTo(1.70).Within(Tolerance), "Val Pb = 1.70 (highest sheet)");
            Assert.That(result.Turn, Is.EqualTo(0.50).Within(Tolerance), "Val Pt = 0.50");
        });
    }

    // M4 — lysine K, window 1: helix = 1.14 (conflict-resolved value), NOT 1.16.
    // Evidence: K 1.14/0.74/1.01 (Przytycka + ref impl) chosen over CSB|SJU 1.16 (Evidence Assumption 1).
    [Test]
    public void PredictSecondaryStructure_SingleResidueLys_ReturnsResolvedHelixPropensity114()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("K", windowSize: 1).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Helix, Is.EqualTo(1.14).Within(Tolerance),
                "Lys Pa = 1.14 (two-source majority over the disputed 1.16)");
            Assert.That(result.Sheet, Is.EqualTo(0.74).Within(Tolerance), "Lys Pb = 0.74");
            Assert.That(result.Turn, Is.EqualTo(1.01).Within(Tolerance), "Lys Pt = 1.01");
        });
    }

    // M5 — two-residue mean "AE", window 2: per-component arithmetic mean.
    // Evidence: helix (1.42+1.51)/2=1.465; sheet (0.83+0.37)/2=0.60; turn (0.66+0.74)/2=0.70 (INV-02).
    [Test]
    public void PredictSecondaryStructure_PairAE_ReturnsExactComponentMeans()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("AE", windowSize: 2).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Helix, Is.EqualTo(1.465).Within(Tolerance), "(1.42+1.51)/2 = 1.465 (INV-02)");
            Assert.That(result.Sheet, Is.EqualTo(0.60).Within(Tolerance), "(0.83+0.37)/2 = 0.60");
            Assert.That(result.Turn, Is.EqualTo(0.70).Within(Tolerance), "(0.66+0.74)/2 = 0.70");
        });
    }

    // M6 — three-residue mean "AEV", window 3.
    // Evidence: helix (1.42+1.51+1.06)/3; sheet (0.83+0.37+1.70)/3; turn (0.66+0.74+0.50)/3 (INV-02).
    [Test]
    public void PredictSecondaryStructure_TripletAEV_ReturnsExactComponentMeans()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("AEV", windowSize: 3).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Helix, Is.EqualTo((1.42 + 1.51 + 1.06) / 3.0).Within(Tolerance),
                "mean helix over A,E,V (INV-02)");
            Assert.That(result.Sheet, Is.EqualTo((0.83 + 0.37 + 1.70) / 3.0).Within(Tolerance),
                "mean sheet over A,E,V");
            Assert.That(result.Turn, Is.EqualTo((0.66 + 0.74 + 0.50) / 3.0).Within(Tolerance),
                "mean turn over A,E,V");
        });
    }

    // M7 — sliding step + count "AEV" window 2: two windows [A,E] and [E,V], step 1.
    // Evidence: INV-03 count = n-w+1 = 3-2+1 = 2; window means per INV-02.
    [Test]
    public void PredictSecondaryStructure_AEVWindow2_YieldsTwoSlidingWindowsInOrder()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("AEV", windowSize: 2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(2), "n-w+1 = 2 windows (INV-03)");
            // window 0 = [A,E]
            Assert.That(result[0].Helix, Is.EqualTo(1.465).Within(Tolerance), "window0 helix (A,E)");
            Assert.That(result[0].Sheet, Is.EqualTo(0.60).Within(Tolerance), "window0 sheet (A,E)");
            Assert.That(result[0].Turn, Is.EqualTo(0.70).Within(Tolerance), "window0 turn (A,E)");
            // window 1 = [E,V]
            Assert.That(result[1].Helix, Is.EqualTo((1.51 + 1.06) / 2.0).Within(Tolerance), "window1 helix (E,V)");
            Assert.That(result[1].Sheet, Is.EqualTo((0.37 + 1.70) / 2.0).Within(Tolerance), "window1 sheet (E,V)");
            Assert.That(result[1].Turn, Is.EqualTo((0.74 + 0.50) / 2.0).Within(Tolerance), "window1 turn (E,V)");
        });
    }

    // M8 — case-insensitive: "ae" equals "AE".
    // Evidence: implementation uppercases input (INV-04).
    [Test]
    public void PredictSecondaryStructure_LowercaseInput_EqualsUppercase()
    {
        var lower = SequenceStatistics.PredictSecondaryStructure("ae", windowSize: 2).Single();
        var upper = SequenceStatistics.PredictSecondaryStructure("AE", windowSize: 2).Single();

        Assert.Multiple(() =>
        {
            Assert.That(lower.Helix, Is.EqualTo(upper.Helix).Within(Tolerance), "case-insensitive helix (INV-04)");
            Assert.That(lower.Sheet, Is.EqualTo(upper.Sheet).Within(Tolerance), "case-insensitive sheet");
            Assert.That(lower.Turn, Is.EqualTo(upper.Turn).Within(Tolerance), "case-insensitive turn");
        });
    }

    // M9 — unknown residue excluded: "AXE" window 3 averages only A and E.
    // Evidence: X has no propensity; excluded from count and mean (INV-05).
    [Test]
    public void PredictSecondaryStructure_UnknownResidueInWindow_ExcludedFromMean()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("AXE", windowSize: 3).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Helix, Is.EqualTo(1.465).Within(Tolerance), "mean of A,E only (X excluded) (INV-05)");
            Assert.That(result.Sheet, Is.EqualTo(0.60).Within(Tolerance), "sheet mean of A,E only");
            Assert.That(result.Turn, Is.EqualTo(0.70).Within(Tolerance), "turn mean of A,E only");
        });
    }

    // M10 — all-unknown window: "XBZ" window 3 emits nothing.
    // Evidence: count = 0 → no tuple emitted (INV-05).
    [Test]
    public void PredictSecondaryStructure_AllUnknownWindow_EmitsNothing()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("XBZ", windowSize: 3).ToList();

        Assert.That(result, Is.Empty, "window with no known residues yields no tuple (INV-05)");
    }

    #endregion

    #region Edge Cases

    // S1 — null input → empty.
    [Test]
    public void PredictSecondaryStructure_NullInput_ReturnsEmpty()
    {
        var result = SequenceStatistics.PredictSecondaryStructure(null!, windowSize: 7).ToList();

        Assert.That(result, Is.Empty, "null sequence yields empty result (INV-06)");
    }

    // S2 — empty input → empty.
    [Test]
    public void PredictSecondaryStructure_EmptyInput_ReturnsEmpty()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("", windowSize: 7).ToList();

        Assert.That(result, Is.Empty, "empty sequence yields empty result (INV-06)");
    }

    // S3 — window larger than sequence → empty.
    [Test]
    public void PredictSecondaryStructure_WindowLargerThanSequence_ReturnsEmpty()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("AE", windowSize: 7).ToList();

        Assert.That(result, Is.Empty, "window > length yields no scan positions (INV-06)");
    }

    // S4 — non-positive window → empty.
    [Test]
    public void PredictSecondaryStructure_NonPositiveWindow_ReturnsEmpty()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("AEV", windowSize: 0).ToList();

        Assert.That(result, Is.Empty, "window < 1 is invalid → empty result (INV-06)");
    }

    #endregion

    #region Biological Sanity (exact-mean derived)

    // C1 — helix-favouring peptide "AEMLK": mean helix > mean sheet.
    // Evidence: helix (1.42+1.51+1.45+1.21+1.14)/5 = 1.346; sheet (0.83+0.37+1.05+1.30+0.74)/5 = 0.858.
    [Test]
    public void PredictSecondaryStructure_HelixFavouringPeptide_HelixMeanExceedsSheetMean()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("AEMLK", windowSize: 5).Single();
        double expectedHelix = (1.42 + 1.51 + 1.45 + 1.21 + 1.14) / 5.0; // 1.346
        double expectedSheet = (0.83 + 0.37 + 1.05 + 1.30 + 0.74) / 5.0; // 0.858

        Assert.Multiple(() =>
        {
            Assert.That(result.Helix, Is.EqualTo(expectedHelix).Within(Tolerance), "exact helix mean");
            Assert.That(result.Sheet, Is.EqualTo(expectedSheet).Within(Tolerance), "exact sheet mean");
            Assert.That(result.Helix, Is.GreaterThan(result.Sheet), "helix formers dominate → mean helix > mean sheet");
        });
    }

    // C2 — sheet-favouring peptide "VIY": mean sheet > mean helix.
    // Evidence: sheet (1.70+1.60+1.47)/3 = 1.59; helix (1.06+1.08+0.69)/3 = 0.9433...
    [Test]
    public void PredictSecondaryStructure_SheetFavouringPeptide_SheetMeanExceedsHelixMean()
    {
        var result = SequenceStatistics.PredictSecondaryStructure("VIY", windowSize: 3).Single();
        double expectedSheet = (1.70 + 1.60 + 1.47) / 3.0;
        double expectedHelix = (1.06 + 1.08 + 0.69) / 3.0;

        Assert.Multiple(() =>
        {
            Assert.That(result.Sheet, Is.EqualTo(expectedSheet).Within(Tolerance), "exact sheet mean");
            Assert.That(result.Helix, Is.EqualTo(expectedHelix).Within(Tolerance), "exact helix mean");
            Assert.That(result.Sheet, Is.GreaterThan(result.Helix), "sheet formers dominate → mean sheet > mean helix");
        });
    }

    #endregion
}
