// ONCO-ACTION-001 — Clinical Actionability Assessment (OncoKB Therapeutic Levels of Evidence)
// Evidence: docs/Evidence/ONCO-ACTION-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-ACTION-001.md
// Source: Chakravarty D et al. (2017). OncoKB: A Precision Oncology Knowledge Base.
//         JCO Precis Oncol 2017:1-16. https://doi.org/10.1200/PO.17.00011
//         OncoKB Therapeutic Levels of Evidence (V2); oncokb-annotator README (HIGHEST_LEVEL order:
//         LEVEL_R1 > LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4 > LEVEL_R2).

using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Oncology;
using Level = Seqeron.Genomics.Oncology.OncologyAnalyzer.OncoKbLevel;
using Assoc = Seqeron.Genomics.Oncology.OncologyAnalyzer.TherapyAssociation;
using Input = Seqeron.Genomics.Oncology.OncologyAnalyzer.VariantActionabilityInput;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_AssessActionability_Tests
{
    // Helper: build a variant with the given leveled drug associations (one drug per level).
    private static Input Variant(params Level[] levels) =>
        new("GENE", "p.X1Y", levels.Select((lvl, i) => new Assoc($"Drug{i}", lvl)).ToArray());

    #region ClassifyActionabilityLevel — combined order R1 > 1 > 2 > 3A > 3B > 4 > R2

    // M1 — combined order full chain: each adjacent pair ranks per HIGHEST_LEVEL order.
    // Expected order (annotator README): R1 > 1 > 2 > 3A > 3B > 4 > R2.
    [Test]
    public void ClassifyActionabilityLevel_CombinedOrder_AdjacentPairsRankCorrectly()
    {
        // Descending chain copied verbatim from the README HIGHEST_LEVEL order — NOT from code output.
        var descending = new[]
        {
            Level.R1, Level.Level1, Level.Level2, Level.Level3A, Level.Level3B, Level.Level4, Level.R2
        };

        Assert.Multiple(() =>
        {
            for (int i = 0; i < descending.Length - 1; i++)
            {
                var higher = descending[i];
                var lower = descending[i + 1];
                // A variant carrying both the higher and lower level must classify to the higher one.
                var result = OncologyAnalyzer.ClassifyActionabilityLevel(Variant(lower, higher));
                Assert.That(result, Is.EqualTo(higher),
                    $"{higher} must outrank {lower} in the OncoKB combined order R1>1>2>3A>3B>4>R2.");
            }
        });
    }

    // M2 — {1, R1} ⇒ combined R1 (R1 outranks Level 1). README: R1 > 1.
    [Test]
    public void ClassifyActionabilityLevel_R1AndLevel1_ReturnsR1()
    {
        var result = OncologyAnalyzer.ClassifyActionabilityLevel(Variant(Level.Level1, Level.R1));

        Assert.That(result, Is.EqualTo(Level.R1),
            "R1 is the highest combined level, above Level 1 (annotator HIGHEST_LEVEL: R1 > 1).");
    }

    // M3 — {4, R2} ⇒ combined Level 4 (Level 4 outranks R2). README: 4 > R2.
    [Test]
    public void ClassifyActionabilityLevel_Level4AndR2_ReturnsLevel4()
    {
        var result = OncologyAnalyzer.ClassifyActionabilityLevel(Variant(Level.R2, Level.Level4));

        Assert.That(result, Is.EqualTo(Level.Level4),
            "Level 4 outranks R2 in the combined order (annotator HIGHEST_LEVEL: 4 > R2).");
    }

    // M11 — {3B, 2, 4} ⇒ highest combined Level 2. README: 2 > 3A > 3B > 4.
    [Test]
    public void ClassifyActionabilityLevel_MixedSet_ReturnsHighest()
    {
        var result = OncologyAnalyzer.ClassifyActionabilityLevel(Variant(Level.Level3B, Level.Level2, Level.Level4));

        Assert.That(result, Is.EqualTo(Level.Level2),
            "Among {3B,2,4} the highest combined level is 2 (order 2 > 3A > 3B > 4).");
    }

    #endregion

    #region AssessActionability — sensitivity / resistance / combined axes

    // M4 — highest sensitive of {2, 3A} is Level 2. README HIGHEST_SENSITIVE_LEVEL: 1 > 2 > 3A > 3B > 4.
    [Test]
    public void AssessActionability_SensitiveLevel2And3A_ReturnsLevel2()
    {
        var result = OncologyAnalyzer.AssessActionability(new[] { Variant(Level.Level2, Level.Level3A) })[0];

        Assert.That(result.HighestSensitiveLevel, Is.EqualTo(Level.Level2),
            "Highest sensitive level of {2,3A} is Level 2 (1 > 2 > 3A > 3B > 4).");
    }

    // M5 — highest sensitive of {3A, 3B, 4} is Level 3A.
    [Test]
    public void AssessActionability_Sensitive3A3B4_ReturnsLevel3A()
    {
        var result = OncologyAnalyzer.AssessActionability(new[] { Variant(Level.Level3A, Level.Level3B, Level.Level4) })[0];

        Assert.That(result.HighestSensitiveLevel, Is.EqualTo(Level.Level3A),
            "Highest sensitive level of {3A,3B,4} is 3A (3A > 3B > 4).");
    }

    // M6 — highest resistance of {R1, R2} is R1. README HIGHEST_RESISTANCE_LEVEL: R1 > R2.
    [Test]
    public void AssessActionability_ResistanceR1R2_ReturnsR1()
    {
        var result = OncologyAnalyzer.AssessActionability(new[] { Variant(Level.R1, Level.R2) })[0];

        Assert.That(result.HighestResistanceLevel, Is.EqualTo(Level.R1),
            "Highest resistance level of {R1,R2} is R1 (R1 > R2).");
    }

    // M7 — variant with sensitivity Level 1 + resistance R1: both axes reported, combined = R1.
    [Test]
    public void AssessActionability_BothAxes_ReportsEachAndCombined()
    {
        var result = OncologyAnalyzer.AssessActionability(new[] { Variant(Level.Level1, Level.R1) })[0];

        Assert.Multiple(() =>
        {
            Assert.That(result.HighestSensitiveLevel, Is.EqualTo(Level.Level1),
                "Sensitive axis sees only Level 1 (R1 is resistance).");
            Assert.That(result.HighestResistanceLevel, Is.EqualTo(Level.R1),
                "Resistance axis sees only R1 (Level 1 is sensitivity).");
            Assert.That(result.HighestCombinedLevel, Is.EqualTo(Level.R1),
                "Combined axis: R1 outranks Level 1 (R1 > 1).");
            Assert.That(result.IsActionable, Is.True, "Variant has leveled associations ⇒ actionable.");
        });
    }

    // M8 — no associations ⇒ None on all axes, not actionable. Annotator leaves HIGHEST_LEVEL empty.
    [Test]
    public void AssessActionability_NoAssociations_NotActionable()
    {
        var result = OncologyAnalyzer.AssessActionability(new[] { Variant( /* no levels */ ) })[0];

        Assert.Multiple(() =>
        {
            Assert.That(result.HighestSensitiveLevel, Is.EqualTo(Level.None), "No associations ⇒ sensitive None.");
            Assert.That(result.HighestResistanceLevel, Is.EqualTo(Level.None), "No associations ⇒ resistance None.");
            Assert.That(result.HighestCombinedLevel, Is.EqualTo(Level.None), "No associations ⇒ combined None.");
            Assert.That(result.IsActionable, Is.False, "No leveled association ⇒ not actionable (VUS-like).");
        });
    }

    // M9 — single Level 1 association.
    [Test]
    public void AssessActionability_SingleLevel1_ReportsLevel1Sensitivity()
    {
        var result = OncologyAnalyzer.AssessActionability(new[] { Variant(Level.Level1) })[0];

        Assert.Multiple(() =>
        {
            Assert.That(result.HighestSensitiveLevel, Is.EqualTo(Level.Level1), "Single Level 1 ⇒ sensitive Level 1.");
            Assert.That(result.HighestResistanceLevel, Is.EqualTo(Level.None), "No resistance association ⇒ None.");
            Assert.That(result.HighestCombinedLevel, Is.EqualTo(Level.Level1), "Combined = Level 1.");
            Assert.That(result.IsActionable, Is.True, "Level 1 ⇒ actionable.");
        });
    }

    // M10 — AssessActionability preserves input order and count.
    [Test]
    public void AssessActionability_PreservesOrderAndCount()
    {
        var v0 = Variant(Level.Level1);
        var v1 = Variant(); // empty
        var v2 = Variant(Level.R2);

        var results = OncologyAnalyzer.AssessActionability(new[] { v0, v1, v2 });

        Assert.Multiple(() =>
        {
            Assert.That(results, Has.Count.EqualTo(3), "One assessment per input variant.");
            Assert.That(results[0].HighestCombinedLevel, Is.EqualTo(Level.Level1), "Index 0 preserved (Level 1).");
            Assert.That(results[1].HighestCombinedLevel, Is.EqualTo(Level.None), "Index 1 preserved (empty).");
            Assert.That(results[2].HighestCombinedLevel, Is.EqualTo(Level.R2), "Index 2 preserved (R2).");
        });
    }

    #endregion

    #region Validation

    // S1 — null variant list throws.
    [Test]
    public void AssessActionability_NullVariants_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.AssessActionability(null!),
            "Null variant enumerable must be rejected.");
    }

    // S2 — null associations list rejected at construction.
    [Test]
    public void VariantActionabilityInput_NullAssociations_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new Input("GENE", "p.X1Y", null!),
            "Null associations list must be rejected at construction.");
    }

    #endregion

    #region GetTherapyRecommendations (Delegate — ordering smoke)

    // S3 — recommendations ordered by descending combined level: {4,1,3A} ⇒ 1, 3A, 4.
    [Test]
    public void GetTherapyRecommendations_OrdersByDescendingLevel()
    {
        var variant = new Input("GENE", "p.X1Y", new[]
        {
            new Assoc("D4", Level.Level4),
            new Assoc("D1", Level.Level1),
            new Assoc("D3A", Level.Level3A),
        });

        var ordered = OncologyAnalyzer.GetTherapyRecommendations(variant);

        Assert.Multiple(() =>
        {
            Assert.That(ordered.Select(a => a.Level),
                Is.EqualTo(new[] { Level.Level1, Level.Level3A, Level.Level4 }),
                "Therapies ordered most-actionable first (1 > 3A > 4).");
            Assert.That(ordered[0].Drug, Is.EqualTo("D1"), "Top recommendation is the Level 1 drug.");
        });
    }

    // S4 — no associations ⇒ empty list (not null).
    [Test]
    public void GetTherapyRecommendations_NoAssociations_ReturnsEmpty()
    {
        var ordered = OncologyAnalyzer.GetTherapyRecommendations(Variant());

        Assert.That(ordered, Is.Empty, "A variant with no associations yields an empty recommendation list.");
    }

    #endregion

    #region IsStandardCare (SOP grouping)

    // C1 — standard-care levels are 1, 2, R1; investigational/hypothetical are 3A, 3B, 4, R2. SOP v3.
    [Test]
    public void IsStandardCare_LevelGrouping_MatchesSop()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IsStandardCare(Level.Level1), Is.True, "Level 1 is standard care.");
            Assert.That(OncologyAnalyzer.IsStandardCare(Level.Level2), Is.True, "Level 2 is standard care.");
            Assert.That(OncologyAnalyzer.IsStandardCare(Level.R1), Is.True, "R1 is standard care (resistance).");
            Assert.That(OncologyAnalyzer.IsStandardCare(Level.Level3A), Is.False, "3A is investigational.");
            Assert.That(OncologyAnalyzer.IsStandardCare(Level.Level3B), Is.False, "3B is investigational.");
            Assert.That(OncologyAnalyzer.IsStandardCare(Level.Level4), Is.False, "4 is hypothetical.");
            Assert.That(OncologyAnalyzer.IsStandardCare(Level.R2), Is.False, "R2 is investigational resistance.");
            Assert.That(OncologyAnalyzer.IsStandardCare(Level.None), Is.False, "None is not standard care.");
        });
    }

    #endregion
}
