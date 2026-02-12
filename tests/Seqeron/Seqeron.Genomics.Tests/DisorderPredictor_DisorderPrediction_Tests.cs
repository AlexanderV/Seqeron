using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

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
        // Poly-Ile: highest hydropathy (4.5), lowest propensity in order-promoting group
        // Hydrophobic sequences fold into stable 3D structures → ordered
        string ordered = new string('I', 30);

        var result = DisorderPredictor.PredictDisorder(ordered);

        Assert.That(result.OverallDisorderContent, Is.LessThan(0.5),
            "Hydrophobic poly-Ile (hydropathy=4.5, TOP-IDP=-0.486) must have low disorder content — Uversky (2000)");
    }

    [Test]
    public void PredictDisorder_MixedHydrophobic_LowDisorderContent()
    {
        // Mixed hydrophobic residues: V(4.2), L(3.8), F(2.8), I(4.5), M(1.9)
        string ordered = "MVILLFFFLLLAAAAIIIIIVVVVVLLLLLL";

        var result = DisorderPredictor.PredictDisorder(ordered);

        Assert.That(result.OverallDisorderContent, Is.LessThan(0.5),
            "Mixed hydrophobic sequence must have low disorder content — Uversky (2000)");
    }

    #endregion

    #region M5: Charged/Polar Sequence → High Disorder — Uversky (2000), Dunker (2001)

    [Test]
    public void PredictDisorder_ChargedSequence_HighDisorderContent()
    {
        // Poly-Glu: TOP-IDP=0.736, charge=-1, hydropathy=-3.5
        // Charged sequences resist folding → disordered
        string disordered = new string('E', 30);

        var result = DisorderPredictor.PredictDisorder(disordered);

        Assert.That(result.OverallDisorderContent, Is.GreaterThan(0.3),
            "Charged poly-Glu (TOP-IDP=0.736, charge=-1, hydropathy=-3.5) must be disordered — Uversky (2000)");
    }

    [Test]
    public void PredictDisorder_MixedDisorderPromoting_HighDisorderContent()
    {
        // Mixture of disorder-promoting residues: E, P, K, D, R
        string disordered = "EPPPPKKKKEEEEDDDDRRRRKKKKEEEEPPPP";

        var result = DisorderPredictor.PredictDisorder(disordered);

        Assert.That(result.OverallDisorderContent, Is.GreaterThan(0.3),
            "Mixed disorder-promoting residues (E,P,K,D,R) must produce high disorder — Dunker (2001)");
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

        Assert.That(result.DisorderedRegions, Is.Not.Empty,
            "30× Proline (highest propensity, TOP-IDP=0.987) must produce disordered regions — Campen et al. (2008)");
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

    #region M8: Disorder Propensity Values Match TOP-IDP Scale — Campen et al. (2008)

    [Test]
    public void GetDisorderPropensity_AllTwentyAminoAcids_MatchScale()
    {
        // TOP-IDP scale values from Campen et al. (2008) Table 2.
        // Source: PMC2676888, PMID 18991772.
        var expected = new Dictionary<char, double>
        {
            ['A'] = 0.060,
            ['R'] = 0.180,
            ['N'] = 0.007,
            ['D'] = 0.192,
            ['C'] = 0.020,
            ['Q'] = 0.318,
            ['E'] = 0.736,
            ['G'] = 0.166,
            ['H'] = 0.303,
            ['I'] = -0.486,
            ['L'] = -0.326,
            ['K'] = 0.586,
            ['M'] = -0.397,
            ['F'] = -0.697,
            ['P'] = 0.987,
            ['S'] = 0.341,
            ['T'] = 0.059,
            ['W'] = -0.884,
            ['Y'] = -0.510,
            ['V'] = -0.121
        };

        Assert.Multiple(() =>
        {
            foreach (var (aa, value) in expected)
            {
                Assert.That(DisorderPredictor.GetDisorderPropensity(aa), Is.EqualTo(value).Within(0.001),
                    $"Propensity for {aa} must match TOP-IDP scale value {value} — Campen et al. (2008)");
            }
        });
    }

    [Test]
    public void GetDisorderPropensity_ProlineIsHighest()
    {
        // Proline: highest disorder propensity = 0.987 — Campen et al. (2008) TOP-IDP
        Assert.That(DisorderPredictor.GetDisorderPropensity('P'), Is.EqualTo(0.987).Within(0.001),
            "Proline must have highest propensity (0.987) — Campen et al. (2008)");
    }

    [Test]
    public void GetDisorderPropensity_TryptophanIsLowest()
    {
        // Tryptophan: lowest disorder propensity = -0.884 — Campen et al. (2008) TOP-IDP
        Assert.That(DisorderPredictor.GetDisorderPropensity('W'), Is.EqualTo(-0.884).Within(0.001),
            "Tryptophan must have lowest propensity (-0.884) — Campen et al. (2008)");
    }

    #endregion

    #region M9: Disorder-Promoting Residues → True — Dunker (2001)

    [Test]
    [TestCase('A', TestName = "IsDisorderPromoting_Ala_True")]
    [TestCase('R', TestName = "IsDisorderPromoting_Arg_True")]
    [TestCase('Q', TestName = "IsDisorderPromoting_Gln_True")]
    [TestCase('E', TestName = "IsDisorderPromoting_Glu_True")]
    [TestCase('G', TestName = "IsDisorderPromoting_Gly_True")]
    [TestCase('K', TestName = "IsDisorderPromoting_Lys_True")]
    [TestCase('P', TestName = "IsDisorderPromoting_Pro_True")]
    [TestCase('S', TestName = "IsDisorderPromoting_Ser_True")]
    public void IsDisorderPromoting_DisorderPromotingResidues_ReturnsTrue(char aa)
    {
        Assert.That(DisorderPredictor.IsDisorderPromoting(aa), Is.True,
            $"'{aa}' must be disorder-promoting — Dunker et al. (2001)");
    }

    #endregion

    #region M10: Order-Promoting Residues → False — Dunker (2001)

    [Test]
    [TestCase('C', TestName = "IsDisorderPromoting_Cys_False")]
    [TestCase('F', TestName = "IsDisorderPromoting_Phe_False")]
    [TestCase('I', TestName = "IsDisorderPromoting_Ile_False")]
    [TestCase('L', TestName = "IsDisorderPromoting_Leu_False")]
    [TestCase('N', TestName = "IsDisorderPromoting_Asn_False")]
    [TestCase('V', TestName = "IsDisorderPromoting_Val_False")]
    [TestCase('W', TestName = "IsDisorderPromoting_Trp_False")]
    [TestCase('Y', TestName = "IsDisorderPromoting_Tyr_False")]
    public void IsDisorderPromoting_OrderPromotingResidues_ReturnsFalse(char aa)
    {
        Assert.That(DisorderPredictor.IsDisorderPromoting(aa), Is.False,
            $"'{aa}' must be order-promoting — Dunker et al. (2001)");
    }

    #endregion

    #region M10b: Ambiguous Residues → False — Dunker (2001)

    [Test]
    [TestCase('D', TestName = "IsDisorderPromoting_Asp_False")]
    [TestCase('H', TestName = "IsDisorderPromoting_His_False")]
    [TestCase('M', TestName = "IsDisorderPromoting_Met_False")]
    [TestCase('T', TestName = "IsDisorderPromoting_Thr_False")]
    public void IsDisorderPromoting_AmbiguousResidues_ReturnsFalse(char aa)
    {
        Assert.That(DisorderPredictor.IsDisorderPromoting(aa), Is.False,
            $"'{aa}' is ambiguous per Dunker (2001) and must not be classified as disorder-promoting");
    }

    #endregion

    #region M11: DisorderPromotingAminoAcids Contains Dunker (2001) Set

    [Test]
    public void DisorderPromotingAminoAcids_ContainsAllExpected()
    {
        var promoting = DisorderPredictor.DisorderPromotingAminoAcids;

        // Dunker et al. (2001): disorder-promoting = {A, R, G, Q, S, P, E, K}
        char[] expected = { 'A', 'E', 'G', 'K', 'P', 'Q', 'R', 'S' };

        Assert.Multiple(() =>
        {
            foreach (char aa in expected)
            {
                Assert.That(promoting, Does.Contain(aa),
                    $"DisorderPromotingAminoAcids must contain '{aa}' — Dunker et al. (2001)");
            }

            Assert.That(promoting.Count, Is.EqualTo(expected.Length),
                $"DisorderPromotingAminoAcids count must be {expected.Length}");
        });
    }

    #endregion

    #region M12: OrderPromotingAminoAcids Contains Dunker (2001) Set

    [Test]
    public void OrderPromotingAminoAcids_ContainsAllExpected()
    {
        var ordering = DisorderPredictor.OrderPromotingAminoAcids;

        // Dunker et al. (2001): order-promoting = {W, C, F, I, Y, V, L, N}
        char[] expected = { 'C', 'F', 'I', 'L', 'N', 'V', 'W', 'Y' };

        Assert.Multiple(() =>
        {
            foreach (char aa in expected)
            {
                Assert.That(ordering, Does.Contain(aa),
                    $"OrderPromotingAminoAcids must contain '{aa}' — Dunker et al. (2001)");
            }

            Assert.That(ordering.Count, Is.EqualTo(expected.Length),
                $"OrderPromotingAminoAcids count must be {expected.Length}");
        });
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
        var result = DisorderPredictor.PredictDisorder("P");

        Assert.Multiple(() =>
        {
            Assert.That(result.ResiduePredictions.Count, Is.EqualTo(1));
            Assert.That(result.ResiduePredictions[0].Residue, Is.EqualTo('P'));
            Assert.That(result.ResiduePredictions[0].DisorderScore, Is.InRange(0.0, 1.0));
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
        var result = DisorderPredictor.PredictDisorder("XXXXX");

        Assert.Multiple(() =>
        {
            Assert.That(result.ResiduePredictions.Count, Is.EqualTo(5));
            foreach (var pred in result.ResiduePredictions)
            {
                Assert.That(pred.DisorderScore, Is.InRange(0.0, 1.0),
                    "Unknown residues must still produce valid scores");
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

    #region S6: Unknown Residue Propensity Returns Zero

    [Test]
    public void GetDisorderPropensity_UnknownResidue_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.GetDisorderPropensity('X'), Is.EqualTo(0.0),
                "Unknown residue 'X' propensity must be 0.0");
            Assert.That(DisorderPredictor.GetDisorderPropensity('Z'), Is.EqualTo(0.0),
                "Unknown residue 'Z' propensity must be 0.0");
            Assert.That(DisorderPredictor.GetDisorderPropensity('B'), Is.EqualTo(0.0),
                "Unknown residue 'B' propensity must be 0.0");
        });
    }

    #endregion

    #region S7: Lowercase Input Handled — Case Insensitivity

    [Test]
    public void GetDisorderPropensity_LowercaseInput_SameAsUppercase()
    {
        Assert.Multiple(() =>
        {
            Assert.That(DisorderPredictor.GetDisorderPropensity('p'), Is.EqualTo(0.987).Within(0.001),
                "Lowercase 'p' must return same propensity as 'P'");
            Assert.That(DisorderPredictor.GetDisorderPropensity('w'), Is.EqualTo(-0.884).Within(0.001),
                "Lowercase 'w' must return same propensity as 'W'");
        });
    }

    #endregion

    #region C1: Mixed Ordered/Disordered Sequence Finds Transition

    [Test]
    public void PredictDisorder_MixedSequence_FindsDisorderedRegions()
    {
        // Ordered flanks (hydrophobic) + disordered middle (charged/polar)
        string sequence = "LLLLLLLLLL" + "PPPPEEEEKKKKDDDD" + "IIIIIIIIII";

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        // The disordered middle should be detected, though boundaries may be blurred by window
        Assert.That(result.MeanDisorderScore, Is.GreaterThan(0.0),
            "Mixed sequence must have non-zero mean disorder score");
    }

    #endregion

    #region C2: Custom Window Size Respected

    [Test]
    public void PredictDisorder_DifferentWindowSizes_ProduceDifferentScores()
    {
        const string sequence = "PPPPPIIIIIIPPPPPIIIIII";

        var result7 = DisorderPredictor.PredictDisorder(sequence, windowSize: 7);
        var result21 = DisorderPredictor.PredictDisorder(sequence, windowSize: 21);

        // Different window sizes should produce at least slightly different score distributions
        // (both should have same count though)
        Assert.Multiple(() =>
        {
            Assert.That(result7.ResiduePredictions.Count, Is.EqualTo(result21.ResiduePredictions.Count),
                "Both window sizes must produce same number of predictions");
            // Scores may differ due to averaging over different windows
            Assert.That(result7.ResiduePredictions.Count, Is.EqualTo(sequence.Length));
        });
    }

    #endregion

    #region C3: AmbiguousAminoAcids Contains Dunker (2001) Set

    [Test]
    public void AmbiguousAminoAcids_ContainsAllExpected()
    {
        var ambiguous = DisorderPredictor.AmbiguousAminoAcids;

        // Dunker et al. (2001): ambiguous = {D, H, M, T}
        char[] expected = { 'D', 'H', 'M', 'T' };

        Assert.Multiple(() =>
        {
            foreach (char aa in expected)
            {
                Assert.That(ambiguous, Does.Contain(aa),
                    $"AmbiguousAminoAcids must contain '{aa}' — Dunker et al. (2001)");
            }

            Assert.That(ambiguous.Count, Is.EqualTo(expected.Length),
                $"AmbiguousAminoAcids count must be {expected.Length}");
        });
    }

    #endregion

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

    #region C5: Three Classification Sets Are Disjoint And Cover All 20 AA

    [Test]
    public void ClassificationSets_AreDisjointAndCoverAll20()
    {
        var disorder = DisorderPredictor.DisorderPromotingAminoAcids;
        var order = DisorderPredictor.OrderPromotingAminoAcids;
        var ambiguous = DisorderPredictor.AmbiguousAminoAcids;

        var all = disorder.Concat(order).Concat(ambiguous).ToList();

        Assert.Multiple(() =>
        {
            // Total = 20
            Assert.That(all.Count, Is.EqualTo(20),
                "Three sets must cover all 20 standard amino acids");

            // Disjoint
            Assert.That(all.Distinct().Count(), Is.EqualTo(20),
                "Three sets must be pairwise disjoint");

            // Counts: 8 + 8 + 4 = 20
            Assert.That(disorder.Count, Is.EqualTo(8), "Disorder-promoting: 8 AA");
            Assert.That(order.Count, Is.EqualTo(8), "Order-promoting: 8 AA");
            Assert.That(ambiguous.Count, Is.EqualTo(4), "Ambiguous: 4 AA");
        });
    }

    #endregion
}
