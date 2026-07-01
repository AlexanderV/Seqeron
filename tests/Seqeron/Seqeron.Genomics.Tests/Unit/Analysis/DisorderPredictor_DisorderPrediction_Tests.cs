using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// DISORDER-PRED-001: Disorder Prediction tests for <see cref="DisorderPredictor"/>.
/// Tests the <see cref="DisorderPredictor.PredictDisorder"/> canonical method,
/// <see cref="DisorderPredictor.GetDisorderPropensity"/>,
/// <see cref="DisorderPredictor.IsDisorderPromoting"/>,
/// and the amino acid classification properties.
///
/// Evidence: Campen et al. (2008) TOP-IDP scale (PMC2676888), Dunker et al. (2001)
/// disorder/order classification, Kyte &amp; Doolittle (1982) hydropathy,
/// Uversky et al. (2000) charge-hydropathy model.
/// See docs/Evidence/DISORDER-PRED-001-Evidence.md for full citations.
/// </summary>
[TestFixture]
public class DisorderPredictor_DisorderPrediction_Tests
{
    #region M1: Empty Sequence Returns Empty Result — Standard Edge Case

    [Test]
    public void PredictDisorder_EmptySequence_ReturnsEmptyResult()
    {
        var result = DisorderPredictor.PredictDisorder("");

        Assert.Multiple(() =>
        {
            Assert.That(result.Sequence, Is.Empty);
            Assert.That(result.ResiduePredictions, Is.Empty);
            Assert.That(result.DisorderedRegions, Is.Empty);
            Assert.That(result.OverallDisorderContent, Is.EqualTo(0.0));
            Assert.That(result.MeanDisorderScore, Is.EqualTo(0.0));
        });
    }

    #endregion

    #region M2: Residue Predictions Count Equals Sequence Length — INV-1

    [Test]
    public void PredictDisorder_AllTwentyAminoAcids_ReturnsCorrectLength()
    {
        const string sequence = "ACDEFGHIKLMNPQRSTVWY";

        var result = DisorderPredictor.PredictDisorder(sequence);

        Assert.That(result.ResiduePredictions.Count, Is.EqualTo(20),
            "INV-1: ResiduePredictions.Count must equal sequence length");
    }

    #endregion

    #region M3: All Disorder Scores in [0, 1] — INV-2

    [Test]
    public void PredictDisorder_AllScoresInZeroOneRange()
    {
        // Test with diverse sequences: ordered, disordered, mixed, unknown
        string[] sequences =
        {
            "ACDEFGHIKLMNPQRSTVWY",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPP",
            "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEK",
            "XXXXX",
            "P"
        };

        foreach (string seq in sequences)
        {
            var result = DisorderPredictor.PredictDisorder(seq);

            foreach (var pred in result.ResiduePredictions)
            {
                Assert.That(pred.DisorderScore, Is.InRange(0.0, 1.0),
                    $"INV-2: Score for residue {pred.Residue} at position {pred.Position} in '{seq}' must be in [0,1]");
            }
        }
    }

    #endregion

    #region M4: Hydrophobic Sequence → Low Disorder — Uversky (2000), Kyte-Doolittle

    [Test]
    public void PredictDisorder_HydrophobicSequence_LowDisorderContent()
    {
        // Poly-Ile: highest hydropathy (4.5), TOP-IDP = -0.486
        // Normalized = (-0.486 - (-0.884)) / 1.871 = 0.2127 < 0.542 → all residues ordered
        // Therefore OverallDisorderContent must be exactly 0.0
        string ordered = new string('I', 30);

        var result = DisorderPredictor.PredictDisorder(ordered);

        Assert.That(result.OverallDisorderContent, Is.EqualTo(0.0),
            "Poly-Ile: normalized TOP-IDP = 0.2127 < cutoff 0.542 → all ordered → content = 0.0 — Campen et al. (2008)");
    }

    [Test]
    public void PredictDisorder_MixedHydrophobic_LowDisorderContent()
    {
        // Mixed hydrophobic residues: M(0.260), V(0.408), I(0.213), L(0.298), F(0.100), A(0.505)
        // All normalized TOP-IDP values < 0.542 → every window averages below cutoff → all ordered
        string ordered = "MVILLFFFLLLAAAAIIIIIVVVVVLLLLLL";

        var result = DisorderPredictor.PredictDisorder(ordered);

        Assert.That(result.OverallDisorderContent, Is.EqualTo(0.0),
            "All residues have normalized TOP-IDP < 0.542 → all ordered → content = 0.0 — Campen et al. (2008)");
    }

    #endregion

    #region M5: Charged/Polar Sequence → High Disorder — Uversky (2000), Dunker (2001)

    [Test]
    public void PredictDisorder_ChargedSequence_HighDisorderContent()
    {
        // Poly-Glu: TOP-IDP = 0.736
        // Normalized = (0.736 - (-0.884)) / 1.871 = 0.8660 >= 0.542 → all residues disordered
        // Therefore OverallDisorderContent must be exactly 1.0
        string disordered = new string('E', 30);

        var result = DisorderPredictor.PredictDisorder(disordered);

        Assert.That(result.OverallDisorderContent, Is.EqualTo(1.0),
            "Poly-Glu: normalized TOP-IDP = 0.8660 >= cutoff 0.542 → all disordered → content = 1.0 — Campen et al. (2008)");
    }

    [Test]
    public void PredictDisorder_MixedDisorderPromoting_HighDisorderContent()
    {
        // Mixture of disorder-promoting residues: E(0.866), P(1.000), K(0.786), D(0.575), R(0.569)
        // All normalized TOP-IDP values >= 0.542 → every window averages above cutoff → all disordered
        string disordered = "EPPPPKKKKEEEEDDDDRRRRKKKKEEEEPPPP";

        var result = DisorderPredictor.PredictDisorder(disordered);

        Assert.That(result.OverallDisorderContent, Is.EqualTo(1.0),
            "All residues have normalized TOP-IDP >= 0.542 → all disordered → content = 1.0 — Campen et al. (2008)");
    }

    #endregion

    #region M6: Proline-Rich → High Disorder — Wikipedia IDP, Dunker (2001)

    [Test]
    public void PredictDisorder_ProlineRich_ProducesDisorderedRegions()
    {
        // Proline has highest disorder propensity (TOP-IDP=0.987)
        // Its cyclic side chain disrupts alpha-helices — Wikipedia Amino acid
        string prolineRich = new string('P', 30);

        var result = DisorderPredictor.PredictDisorder(prolineRich, minRegionLength: 5);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverallDisorderContent, Is.EqualTo(1.0),
                "Poly-Pro: normalized TOP-IDP = 1.0 ≥ cutoff 0.542 → all disordered → content = 1.0 — Campen et al. (2008)");
            Assert.That(result.DisorderedRegions, Is.Not.Empty,
                "30× Proline must produce at least one disordered region");
        });
    }

    #endregion

    #region M7: Case Insensitivity — INV-5

    [Test]
    public void PredictDisorder_CaseInsensitive_SameResults()
    {
        var upper = DisorderPredictor.PredictDisorder("PPPPEEEE");
        var lower = DisorderPredictor.PredictDisorder("ppppeeee");
        var mixed = DisorderPredictor.PredictDisorder("PpPpEeEe");

        Assert.Multiple(() =>
        {
            Assert.That(lower.MeanDisorderScore, Is.EqualTo(upper.MeanDisorderScore).Within(0.001),
                "INV-5: Lower case must produce same score as upper case");
            Assert.That(mixed.MeanDisorderScore, Is.EqualTo(upper.MeanDisorderScore).Within(0.001),
                "INV-5: Mixed case must produce same score as upper case");
        });
    }

    #endregion

    // Tests for GetDisorderPropensity / IsDisorderPromoting / classification properties
    // moved to DisorderPredictor_GetDisorderPropensity_Tests.cs (unit DISORDER-PROPENSITY-001),
    // which owns those four methods per the Method Index in ALGORITHMS_CHECKLIST_V2.md.

    #region M8b: Normalized TOP-IDP Score Formula — Campen et al. (2008)

    [Test]
    public void PredictDisorder_HomopolymericSequences_MatchExactNormalizedTopIdpScores()
    {
        // Theory: S = (TOP-IDP(c) - TOP-IDP_min) / (TOP-IDP_max - TOP-IDP_min)
        // where TOP-IDP_min = -0.884 (W), TOP-IDP_max = 0.987 (P), range = 1.871
        // Source: Campen et al. (2008) PMC2676888, Table 2

        // Poly-Trp: lowest propensity → normalized = (-0.884 + 0.884) / 1.871 = 0.0
        var resultW = DisorderPredictor.PredictDisorder(new string('W', 30));
        Assert.That(resultW.ResiduePredictions[15].DisorderScore, Is.EqualTo(0.0).Within(0.0001),
            "Poly-Trp: normalized TOP-IDP = 0.0 (minimum anchor) — Campen et al. (2008)");

        // Poly-Pro: highest propensity → normalized = (0.987 + 0.884) / 1.871 = 1.0
        var resultP = DisorderPredictor.PredictDisorder(new string('P', 30));
        Assert.That(resultP.ResiduePredictions[15].DisorderScore, Is.EqualTo(1.0).Within(0.0001),
            "Poly-Pro: normalized TOP-IDP = 1.0 (maximum anchor) — Campen et al. (2008)");

        // Poly-Glu: normalized = (0.736 + 0.884) / 1.871 = 1.620 / 1.871 ≈ 0.8660
        var resultE = DisorderPredictor.PredictDisorder(new string('E', 30));
        Assert.That(resultE.ResiduePredictions[15].DisorderScore, Is.EqualTo(0.8660).Within(0.0005),
            "Poly-Glu: normalized TOP-IDP = 0.8660 — Campen et al. (2008)");

        // Poly-Ile: normalized = (-0.486 + 0.884) / 1.871 = 0.398 / 1.871 ≈ 0.2127
        var resultI = DisorderPredictor.PredictDisorder(new string('I', 30));
        Assert.That(resultI.ResiduePredictions[15].DisorderScore, Is.EqualTo(0.2127).Within(0.0005),
            "Poly-Ile: normalized TOP-IDP = 0.2127 — Campen et al. (2008)");
    }

    #endregion

    #region M13: Residue Predictions Have Correct Positions — Mathematical Identity

    [Test]
    public void PredictDisorder_ResiduePredictionsHaveCorrectPositionsAndResidues()
    {
        const string sequence = "AAAKKKEEE";

        var result = DisorderPredictor.PredictDisorder(sequence);

        Assert.Multiple(() =>
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                Assert.That(result.ResiduePredictions[i].Position, Is.EqualTo(i),
                    $"Position at index {i} must be {i}");
                Assert.That(result.ResiduePredictions[i].Residue, Is.EqualTo(sequence[i]),
                    $"Residue at index {i} must be '{sequence[i]}'");
            }
        });
    }

    #endregion

    #region S1: Single Residue Handled — Boundary Condition

    [Test]
    public void PredictDisorder_SingleResidue_ReturnsOnePrediction()
    {
        // Single Pro: normalized TOP-IDP = (0.987 + 0.884) / 1.871 = 1.0 — Campen et al. (2008)
        var result = DisorderPredictor.PredictDisorder("P");

        Assert.Multiple(() =>
        {
            Assert.That(result.ResiduePredictions.Count, Is.EqualTo(1));
            Assert.That(result.ResiduePredictions[0].Residue, Is.EqualTo('P'));
            Assert.That(result.ResiduePredictions[0].DisorderScore, Is.EqualTo(1.0).Within(0.0001),
                "Single Pro: normalized TOP-IDP = (0.987+0.884)/1.871 = 1.0 — Campen et al. (2008)");
        });
    }

    #endregion

    #region S2: Short Sequence Handled — Boundary Condition

    [Test]
    public void PredictDisorder_ShortSequence_HandledCorrectly()
    {
        var result = DisorderPredictor.PredictDisorder("PPPP");

        Assert.Multiple(() =>
        {
            Assert.That(result.ResiduePredictions.Count, Is.EqualTo(4));
            Assert.That(result.Sequence, Is.EqualTo("PPPP"));
        });
    }

    #endregion

    #region S3: Unknown Residues Handled — Tolerance

    [Test]
    public void PredictDisorder_UnknownResidues_HandledGracefully()
    {
        // All-unknown: no recognized AA → CalculateDisorderScore returns 0.0
        var result = DisorderPredictor.PredictDisorder("XXXXX");

        Assert.Multiple(() =>
        {
            Assert.That(result.ResiduePredictions.Count, Is.EqualTo(5));
            foreach (var pred in result.ResiduePredictions)
            {
                Assert.That(pred.DisorderScore, Is.EqualTo(0.0).Within(0.0001),
                    "Unknown residues contribute nothing → score = 0.0");
            }
        });
    }

    #endregion

    #region S4: MeanDisorderScore Is Average of Residue Scores — INV-4

    [Test]
    public void PredictDisorder_MeanDisorderScore_IsAverageOfResidueScores()
    {
        const string sequence = "ACDEFGHIKLMNPQRSTVWY";
        var result = DisorderPredictor.PredictDisorder(sequence);

        double expectedMean = result.ResiduePredictions.Average(r => r.DisorderScore);

        Assert.That(result.MeanDisorderScore, Is.EqualTo(expectedMean).Within(0.0001),
            "INV-4: MeanDisorderScore must equal the average of all ResiduePrediction scores");
    }

    #endregion

    #region S5: OverallDisorderContent Is Fraction of Disordered Residues — INV-3

    [Test]
    public void PredictDisorder_OverallDisorderContent_IsFraction()
    {
        const string sequence = "PPPPPPPPPPIIIIIIIIIIIIIIIIIIIII";
        var result = DisorderPredictor.PredictDisorder(sequence);

        int disorderedCount = result.ResiduePredictions.Count(r => r.IsDisordered);
        double expectedContent = disorderedCount / (double)sequence.Length;

        Assert.That(result.OverallDisorderContent, Is.EqualTo(expectedContent).Within(0.0001),
            "INV-3: OverallDisorderContent must equal disorderedCount / sequenceLength");
    }

    #endregion

    // GetDisorderPropensity edge cases (unknown residue, lowercase) moved to
    // DisorderPredictor_GetDisorderPropensity_Tests.cs (DISORDER-PROPENSITY-001).

    #region C1: Mixed Ordered/Disordered Sequence Finds Transition

    [Test]
    public void PredictDisorder_MixedSequence_FindsDisorderedRegions()
    {
        // Ordered flanks (hydrophobic) + disordered middle (charged/polar)
        string sequence = "LLLLLLLLLL" + "PPPPEEEEKKKKDDDD" + "IIIIIIIIII";

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.Multiple(() =>
        {
            Assert.That(result.MeanDisorderScore, Is.GreaterThan(0.0),
                "Mixed sequence must have non-zero mean disorder score");
            Assert.That(result.DisorderedRegions, Is.Not.Empty,
                "Disordered middle (P/E/K/D, all high TOP-IDP) must produce at least one region");
        });
    }

    #endregion

    #region C2: Custom Window Size Respected

    [Test]
    public void PredictDisorder_DifferentWindowSizes_ProduceDifferentScores()
    {
        const string sequence = "PPPPPIIIIIIPPPPPIIIIII";

        var result7 = DisorderPredictor.PredictDisorder(sequence, windowSize: 7);
        var result21 = DisorderPredictor.PredictDisorder(sequence, windowSize: 21);

        // Different window sizes should produce different score distributions
        // because smaller windows capture local composition more sharply
        Assert.Multiple(() =>
        {
            Assert.That(result7.ResiduePredictions.Count, Is.EqualTo(result21.ResiduePredictions.Count),
                "Both window sizes must produce same number of predictions");
            Assert.That(result7.ResiduePredictions.Count, Is.EqualTo(sequence.Length));

            bool anyDifferent = false;
            for (int i = 0; i < result7.ResiduePredictions.Count; i++)
            {
                if (Math.Abs(result7.ResiduePredictions[i].DisorderScore
                    - result21.ResiduePredictions[i].DisorderScore) > 0.0001)
                {
                    anyDifferent = true;
                    break;
                }
            }
            Assert.That(anyDifferent, Is.True,
                "Different window sizes on mixed P/I sequence must produce different per-residue scores");
        });
    }

    #endregion

    // AmbiguousAminoAcids and the classification-set disjointness tests moved to
    // DisorderPredictor_GetDisorderPropensity_Tests.cs (DISORDER-PROPENSITY-001).

    #region C4: CalculateHydropathy Returns Mean Kyte-Doolittle Value

    [Test]
    public void CalculateHydropathy_SingleResidue_ReturnsKyteDoolittleValue()
    {
        // Kyte & Doolittle (1982): I=4.5, W=-0.9, E=-3.5
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.CalculateHydropathy("I"), Is.EqualTo(4.5).Within(0.001));
            Assert.That(DisorderPredictor.CalculateHydropathy("W"), Is.EqualTo(-0.9).Within(0.001));
            Assert.That(DisorderPredictor.CalculateHydropathy("E"), Is.EqualTo(-3.5).Within(0.001));
        });
    }

    [Test]
    public void CalculateHydropathy_EmptySequence_ReturnsZero()
    {
        Assert.That(DisorderPredictor.CalculateHydropathy(""), Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateHydropathy_MixedSequence_ReturnsMean()
    {
        // "AI": A=1.8, I=4.5 → mean = (1.8+4.5)/2 = 3.15
        Assert.That(DisorderPredictor.CalculateHydropathy("AI"), Is.EqualTo(3.15).Within(0.001));
    }

    [Test]
    public void CalculateHydropathy_CaseInsensitive()
    {
        Assert.That(DisorderPredictor.CalculateHydropathy("ai"),
            Is.EqualTo(DisorderPredictor.CalculateHydropathy("AI")).Within(0.001));
    }

    #endregion
}
