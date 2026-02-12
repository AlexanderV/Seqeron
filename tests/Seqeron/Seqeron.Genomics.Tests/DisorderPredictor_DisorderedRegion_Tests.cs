using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// DISORDER-REGION-001: Disordered Region Detection
///
/// Tests for <c>IdentifyDisorderedRegions</c> (canonical) and
/// <c>ClassifyDisorderedRegion</c> (classification) from
/// <see cref="DisorderPredictor"/>. Both methods are private static —
/// tested indirectly through the public <c>PredictDisorder()</c> method.
///
/// Sources:
///   Campen et al. (2008) TOP-IDP Scale — PMC2676888, Table 2
///   van der Lee et al. (2014) IDR classification — PMC4095912
///   Ward et al. (2004) Long IDR significance — doi:10.1016/j.jmb.2004.02.002
///   Dunker et al. (2001) disorder/order classification — PMID 11381529
/// </summary>
[TestFixture]
public class DisorderPredictor_DisorderedRegion_Tests
{
    #region DISORDER-REGION-001 — M1–M6: Region Detection

    /// <summary>
    /// M1 — All ordered residues produce no disordered regions.
    /// 30×W has the lowest TOP-IDP propensity (−0.884) → normalized score = 0.0.
    /// Source: Campen et al. (2008) Table 2.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_AllOrdered_NoRegions()
    {
        // W: propensity −0.884, normalized 0.0 — most order-promoting
        string allOrdered = new string('W', 30);

        var result = DisorderPredictor.PredictDisorder(allOrdered, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Is.Empty,
            "30×W (normalized score 0.0) must produce zero disordered regions");
    }

    /// <summary>
    /// M2 — All disordered residues produce exactly one region spanning the full sequence.
    /// 30×P has the highest TOP-IDP propensity (0.987) → normalized score = 1.0.
    /// Source: Campen et al. (2008) Table 2.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_AllDisordered_OneRegion()
    {
        // P: propensity 0.987, normalized 1.0 — most disorder-promoting
        string allDisordered = new string('P', 30);

        var result = DisorderPredictor.PredictDisorder(allDisordered, minRegionLength: 5);

        Assert.Multiple(() =>
        {
            Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
                "30×P (normalized score 1.0) must produce exactly one region");
            Assert.That(result.DisorderedRegions[0].Start, Is.EqualTo(0));
            Assert.That(result.DisorderedRegions[0].End, Is.EqualTo(29));
        });
    }

    /// <summary>
    /// M3 — Region boundaries are exact for a homopolymeric disordered sequence.
    /// 30×E: normalized score ≈ 0.866 → all positions disordered → one region [0, 29].
    /// Source: Campen et al. (2008) Table 2, E propensity = 0.736.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_BoundariesCorrect()
    {
        // E: propensity 0.736, normalized = (0.736+0.884)/1.871 ≈ 0.866 → all disordered
        string sequence = new string('E', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "30×E (normalized ≈ 0.866) must produce exactly one region");
        Assert.Multiple(() =>
        {
            Assert.That(result.DisorderedRegions[0].Start, Is.EqualTo(0),
                "Homopolymeric E: region must start at position 0");
            Assert.That(result.DisorderedRegions[0].End, Is.EqualTo(29),
                "Homopolymeric E: region must end at position 29");
        });
    }

    /// <summary>
    /// M4 — MeanScore for a homopolymeric disordered region equals the normalized TOP-IDP value.
    /// 30×P → every window contains only P → DisorderScore = 1.0 at each position → MeanScore = 1.0.
    /// Source: Campen et al. (2008) Table 2.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_MeanScoreIsAverage()
    {
        string sequence = new string('P', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1));
        // Homopolymeric P: every window is all P → DisorderScore = 1.0 → MeanScore = 1.0
        Assert.That(result.DisorderedRegions[0].MeanScore, Is.EqualTo(1.0).Within(0.001),
            "MeanScore for 30×P must equal 1.0 (all residue scores are 1.0)");
    }

    /// <summary>
    /// M5 — Regions shorter than minRegionLength are excluded.
    /// 30×P with minRegionLength=31 → region length (30) &lt; 31 → no regions.
    /// Source: algorithm definition.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_MinLengthFiltering()
    {
        string sequence = new string('P', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 31);

        Assert.That(result.DisorderedRegions, Is.Empty,
            "Region of 30 residues must be excluded when minRegionLength=31");
    }

    /// <summary>
    /// M6 — Trailing disordered region (at end of sequence) is captured.
    /// W(10)+P(20): P block at positions 10–29. With window=21 (halfWindow=10),
    /// position 11 is the first where ≥12 of 21 window residues are P → score ≈ 0.571 &gt; 0.542.
    /// Region = [11, 29].
    /// Source: algorithm definition — handles "region at end" via end-of-loop check.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_TrailingRegionCaptured()
    {
        // W×10 (ordered) + P×20 (disordered) = 30 residues
        string sequence = new string('W', 10) + new string('P', 20);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "W10+P20 must produce exactly one trailing region");
        Assert.Multiple(() =>
        {
            Assert.That(result.DisorderedRegions[0].Start, Is.EqualTo(11),
                "Window boundary: position 11 has 12P/21 → score 0.571 > 0.542");
            Assert.That(result.DisorderedRegions[0].End, Is.EqualTo(29),
                "Trailing region must extend to the final residue");
        });
    }

    #endregion

    #region DISORDER-REGION-001 — M7–M12: Region Classification

    /// <summary>
    /// M7 — Proline-rich classification.
    /// 30×P → Pro fraction = 1.0 &gt; 0.25 → "Proline-rich".
    /// Source: van der Lee et al. (2014) PMC4095912.
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_ProlineRich()
    {
        string sequence = new string('P', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Is.Not.Empty);
        Assert.That(result.DisorderedRegions[0].RegionType, Is.EqualTo("Proline-rich"),
            "30×P (Pro fraction 1.0) must classify as Proline-rich");
    }

    /// <summary>
    /// M8 — Acidic classification.
    /// 30×E → E/D fraction = 1.0 &gt; 0.25 → "Acidic".
    /// Source: van der Lee et al. (2014) PMC4095912.
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_Acidic()
    {
        // E: propensity 0.736, normalized ≈ 0.866, all disordered
        string sequence = new string('E', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "30×E (normalized ≈ 0.866) must produce exactly one region");
        Assert.That(result.DisorderedRegions[0].RegionType, Is.EqualTo("Acidic"),
            "30×E (E/D fraction 1.0) must classify as Acidic");
    }

    /// <summary>
    /// M9 — Basic classification.
    /// 30×K → K/R fraction = 1.0 &gt; 0.25 → "Basic".
    /// K: propensity 0.586, normalized ≈ 0.786 → all disordered.
    /// Source: van der Lee et al. (2014) PMC4095912.
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_Basic()
    {
        string sequence = new string('K', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "30×K (normalized ≈ 0.786) must produce exactly one region");
        Assert.That(result.DisorderedRegions[0].RegionType, Is.EqualTo("Basic"),
            "30×K (K/R fraction 1.0) must classify as Basic");
    }

    /// <summary>
    /// M10 — Ser/Thr-rich classification.
    /// 30×S → S/T fraction = 1.0 &gt; 0.25 → "Ser/Thr-rich".
    /// S: propensity 0.341, normalized ≈ 0.655 → all disordered.
    /// Source: van der Lee et al. (2014) PMC4095912.
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_SerThrRich()
    {
        string sequence = new string('S', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "30×S (normalized ≈ 0.655) must produce exactly one region");
        Assert.That(result.DisorderedRegions[0].RegionType, Is.EqualTo("Ser/Thr-rich"),
            "30×S (S/T fraction 1.0) must classify as Ser/Thr-rich");
    }

    /// <summary>
    /// M11 — Long IDR classification.
    /// (EKQSP)×8 = 40 residues. Each AA at 20% fraction (none &gt; 0.25), length &gt; 30 → "Long IDR".
    /// All AAs are disorder-promoting with normalized scores &gt; 0.542.
    /// Source: Ward et al. (2004); van der Lee et al. (2014).
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_LongIdr()
    {
        // 5 disorder-promoting AAs: E(0.866), K(0.786), Q(0.642), S(0.655), P(1.0)
        // Each at 8/40 = 0.20 fraction → none > 0.25; length 40 > 30 → "Long IDR"
        string sequence = string.Concat(Enumerable.Repeat("EKQSP", 8));
        Assert.That(sequence.Length, Is.EqualTo(40), "Precondition: sequence length");

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "All residues disordered (avg ≈ 0.79) → exactly one region");
        Assert.Multiple(() =>
        {
            Assert.That(result.DisorderedRegions[0].Start, Is.EqualTo(0));
            Assert.That(result.DisorderedRegions[0].End, Is.EqualTo(39));
            Assert.That(result.DisorderedRegions[0].RegionType, Is.EqualTo("Long IDR"),
                "Region of 40 residues with no dominant AA must classify as Long IDR");
        });
    }

    /// <summary>
    /// M12 — Standard IDR classification.
    /// (EKQSP)×4 = 20 residues. Each AA at 20% fraction (none &gt; 0.25), length ≤ 30 → "Standard IDR".
    /// Source: classification fallback definition.
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_StandardIdr()
    {
        // Same amino acid mix as M11 but shorter: 20 residues ≤ 30
        string sequence = string.Concat(Enumerable.Repeat("EKQSP", 4));
        Assert.That(sequence.Length, Is.EqualTo(20), "Precondition: sequence length");

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "All residues disordered → exactly one region");
        Assert.Multiple(() =>
        {
            Assert.That(result.DisorderedRegions[0].Start, Is.EqualTo(0));
            Assert.That(result.DisorderedRegions[0].End, Is.EqualTo(19));
            Assert.That(result.DisorderedRegions[0].RegionType, Is.EqualTo("Standard IDR"),
                "Region of 20 residues (≤30) with no dominant AA must classify as Standard IDR");
        });
    }

    #endregion

    #region DISORDER-REGION-001 — M13–M14: Structural Invariants

    /// <summary>
    /// M13 — Confidence values are in [0, 1] for all regions.
    /// Confidence = (meanScore − cutoff) / (1.0 − cutoff), clamped to [0, 1].
    /// Source: Campen et al. (2008) TOP-IDP scale; normalized distance from cutoff.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_ConfidenceInRange()
    {
        // Test across multiple sequence types
        string[] sequences =
        {
            new string('P', 30),               // High confidence: score=1, length=30
            new string('P', 5),                 // Lower confidence: length=5
            new string('E', 30),                // Moderate score
            new string('S', 30),                // Lower score, still disordered
        };

        foreach (string seq in sequences)
        {
            var result = DisorderPredictor.PredictDisorder(seq, minRegionLength: 5);
            foreach (var region in result.DisorderedRegions)
            {
                Assert.That(region.Confidence, Is.InRange(0.0, 1.0),
                    $"Confidence must be in [0,1] for sequence starting with '{seq[0]}' (len={seq.Length})");
            }
        }
    }

    /// <summary>
    /// M14 — Regions are non-overlapping and sorted by Start.
    /// W(15)+P(20)+W(15)+P(20)+W(15) → two separated disordered regions.
    /// No region's Start may be ≤ a previous region's End.
    /// Source: single-pass scan guarantees non-overlapping.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_NonOverlapping()
    {
        // Two disordered blocks separated by ordered blocks
        string sequence = new string('W', 15) + new string('P', 20)
                        + new string('W', 15) + new string('P', 20)
                        + new string('W', 15);
        Assert.That(sequence.Length, Is.EqualTo(85), "Precondition: sequence length");

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions.Count, Is.EqualTo(2),
            "Two separated P blocks must produce exactly two regions");

        for (int i = 1; i < result.DisorderedRegions.Count; i++)
        {
            Assert.Multiple(() =>
            {
                Assert.That(result.DisorderedRegions[i].Start,
                    Is.GreaterThan(result.DisorderedRegions[i - 1].End),
                    $"Region {i} Start must be > Region {i - 1} End (non-overlapping)");
                Assert.That(result.DisorderedRegions[i].Start,
                    Is.GreaterThanOrEqualTo(result.DisorderedRegions[i - 1].Start),
                    $"Regions must be sorted by Start");
            });
        }
    }

    #endregion

    #region DISORDER-REGION-001 — S1–S5: SHOULD Tests

    /// <summary>
    /// S2 — Region at the start of the sequence.
    /// P(20)+W(30): disordered segment at start → region.Start = 0.
    /// Source: algorithm definition.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_RegionAtStart()
    {
        string sequence = new string('P', 20) + new string('W', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Is.Not.Empty,
            "Leading P block must produce a region");
        Assert.That(result.DisorderedRegions[0].Start, Is.EqualTo(0),
            "First region must start at position 0");
    }

    /// <summary>
    /// S3 — Region of exactly minRegionLength is included.
    /// 5×P with minRegionLength=5 → region of 5 ≥ 5 → included.
    /// Source: algorithm definition (≥ not >).
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_ExactMinLength_Included()
    {
        string sequence = new string('P', 5);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "Region of exactly minRegionLength=5 must be included");
    }

    /// <summary>
    /// S4 — Region just below minRegionLength is excluded.
    /// 4×P with minRegionLength=5 → region length (4) &lt; 5 → excluded.
    /// Source: algorithm definition.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_BelowMinLength_Excluded()
    {
        string sequence = new string('P', 4);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Is.Empty,
            "Region of 4 residues must be excluded when minRegionLength=5");
    }

    /// <summary>
    /// S5 — Mixed ordered-disordered-ordered sequence: exact central region boundaries.
    /// W(15)+P(20)+W(15) = 50 residues. P block at positions 15–34.
    /// Window=21 (halfWindow=10): position 16 is the first where 12 of 21 residues are P
    /// → score = 12/21 ≈ 0.571 &gt; 0.542. Position 33 is the last disordered.
    /// Region = [16, 33], length = 18.
    /// Source: Campen et al. (2008), window boundary analysis.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_CentralRegion()
    {
        string sequence = new string('W', 15) + new string('P', 20) + new string('W', 15);
        Assert.That(sequence.Length, Is.EqualTo(50), "Precondition: sequence length");

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1),
            "Central P block flanked by W must produce exactly one region");
        Assert.Multiple(() =>
        {
            Assert.That(result.DisorderedRegions[0].Start, Is.EqualTo(16),
                "Position 16: 12P in 21-window → score 0.571 > 0.542");
            Assert.That(result.DisorderedRegions[0].End, Is.EqualTo(33),
                "Position 33: 12P in 21-window → last disordered position");
        });
    }

    #endregion

    #region DISORDER-REGION-001 — C1–C3: COULD Tests

    /// <summary>
    /// C1 — Classification priority: Proline-rich wins over Acidic.
    /// When Pro fraction &gt; 0.25 AND E/D fraction &gt; 0.25, Proline-rich is returned.
    /// Priority: most specific single-AA bias first (Proline has highest TOP-IDP
    /// propensity, 0.987 — Campen et al. 2008) before charge-based classes
    /// (Das &amp; Pappu 2013, f+/f− diagram-of-states).
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_Priority_ProlineOverAcidic()
    {
        // 15P + 15E = 30 chars: P fraction = 0.5, E/D fraction = 0.5 → both > 0.25
        string sequence = new string('P', 15) + new string('E', 15);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1));
        Assert.That(result.DisorderedRegions[0].RegionType, Is.EqualTo("Proline-rich"),
            "When both Pro and Acidic fractions > 0.25, Proline-rich wins (priority order)");
    }

    /// <summary>
    /// C2 — Empty sequence produces no disordered regions.
    /// Source: trivial boundary condition.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_EmptySequence_NoRegions()
    {
        var result = DisorderPredictor.PredictDisorder("");

        Assert.That(result.DisorderedRegions, Is.Empty,
            "Empty sequence must produce zero disordered regions");
    }

    /// <summary>
    /// C3 — Classification priority: Acidic wins over Basic.
    /// When E/D fraction &gt; 0.25 AND K/R fraction &gt; 0.25, Acidic is returned.
    /// Acidic is checked before Basic in the classification chain:
    /// Pro &gt; Acidic &gt; Basic &gt; S/T &gt; Long IDR &gt; Standard IDR.
    /// Source: Das &amp; Pappu (2013) f+/f− diagram-of-states; charge-driven
    /// classification checks net-negative (Acidic) before net-positive (Basic).
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_Priority_AcidicOverBasic()
    {
        // 15E + 15K = 30 chars: E/D fraction = 0.5, K/R fraction = 0.5 → both > 0.25
        string sequence = new string('E', 15) + new string('K', 15);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1));
        Assert.That(result.DisorderedRegions[0].RegionType, Is.EqualTo("Acidic"),
            "When both Acidic and Basic fractions > 0.25, Acidic wins (priority order)");
    }

    #endregion

    #region DISORDER-REGION-001 — Additional Invariant Tests

    /// <summary>
    /// INV-5 — MeanScore for homopolymeric sequences equals the exact normalized TOP-IDP value.
    /// Normalized = (propensity − (−0.884)) / 1.871.
    /// P: 1.0, E: (0.736+0.884)/1.871 ≈ 0.866, K: (0.586+0.884)/1.871 ≈ 0.786, S: (0.341+0.884)/1.871 ≈ 0.655.
    /// Source: Campen et al. (2008) Table 2.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_MeanScoreExactValues()
    {
        // P: propensity 0.987 → normalized 1.0
        var pResult = DisorderPredictor.PredictDisorder(new string('P', 30), minRegionLength: 5);
        Assert.That(pResult.DisorderedRegions[0].MeanScore, Is.EqualTo(1.0).Within(0.001),
            "30×P MeanScore must equal 1.0");

        // E: propensity 0.736 → normalized (0.736+0.884)/1.871 ≈ 0.866
        var eResult = DisorderPredictor.PredictDisorder(new string('E', 30), minRegionLength: 5);
        Assert.That(eResult.DisorderedRegions[0].MeanScore, Is.EqualTo(0.866).Within(0.01),
            "30×E MeanScore must equal (0.736+0.884)/1.871");

        // K: propensity 0.586 → normalized (0.586+0.884)/1.871 ≈ 0.786
        var kResult = DisorderPredictor.PredictDisorder(new string('K', 30), minRegionLength: 5);
        Assert.That(kResult.DisorderedRegions[0].MeanScore, Is.EqualTo(0.786).Within(0.01),
            "30×K MeanScore must equal (0.586+0.884)/1.871");

        // S: propensity 0.341 → normalized (0.341+0.884)/1.871 ≈ 0.655
        var sResult = DisorderPredictor.PredictDisorder(new string('S', 30), minRegionLength: 5);
        Assert.That(sResult.DisorderedRegions[0].MeanScore, Is.EqualTo(0.655).Within(0.01),
            "30×S MeanScore must equal (0.341+0.884)/1.871");
    }

    /// <summary>
    /// INV-7 — RegionType is one of the six valid classification labels.
    /// Source: ClassifyDisorderedRegion definition.
    /// </summary>
    [Test]
    public void ClassifyDisorderedRegion_ValidLabels()
    {
        string[] validTypes =
        {
            "Proline-rich", "Acidic", "Basic", "Ser/Thr-rich", "Long IDR", "Standard IDR"
        };

        // Test multiple sequence types to exercise different classification paths
        string[] sequences =
        {
            new string('P', 30),                                    // Proline-rich
            new string('E', 30),                                    // Acidic
            new string('K', 30),                                    // Basic
            new string('S', 30),                                    // Ser/Thr-rich
            string.Concat(Enumerable.Repeat("EKQSP", 8)),          // Long IDR
            string.Concat(Enumerable.Repeat("EKQSP", 4)),          // Standard IDR
        };

        foreach (string seq in sequences)
        {
            var result = DisorderPredictor.PredictDisorder(seq, minRegionLength: 5);
            foreach (var region in result.DisorderedRegions)
            {
                Assert.That(validTypes, Does.Contain(region.RegionType),
                    $"RegionType '{region.RegionType}' must be one of the six valid labels");
            }
        }
    }

    /// <summary>
    /// INV-3 — Every region length ≥ minRegionLength.
    /// Source: algorithm definition.
    /// </summary>
    [Test]
    public void IdentifyDisorderedRegions_AllRegionsAboveMinLength()
    {
        const int minLen = 5;
        string sequence = new string('W', 15) + new string('P', 20)
                        + new string('W', 15) + new string('P', 20)
                        + new string('W', 15);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: minLen);

        foreach (var region in result.DisorderedRegions)
        {
            int length = region.End - region.Start + 1;
            Assert.That(length, Is.GreaterThanOrEqualTo(minLen),
                $"Region [{region.Start}–{region.End}] length {length} must be ≥ {minLen}");
        }
    }

    /// <summary>
    /// Confidence for a high-scoring region approaches 1.0.
    /// 30×P: meanScore ≈ 1.0, confidence = (1.0−0.542)/(1.0−0.542) = 1.0.
    /// Source: Campen et al. (2008) TOP-IDP scale; normalized distance from cutoff.
    /// </summary>
    [Test]
    public void CalculateConfidence_HighScoreLongRegion_ApproachesOne()
    {
        string sequence = new string('P', 30);

        var result = DisorderPredictor.PredictDisorder(sequence, minRegionLength: 5);

        Assert.That(result.DisorderedRegions, Has.Count.EqualTo(1));
        Assert.That(result.DisorderedRegions[0].Confidence, Is.EqualTo(1.0).Within(0.01),
            "30×P: scoreConfidence=1.0 → confidence=1.0");
    }

    /// <summary>
    /// Confidence for a lower-scoring disordered region is less than for a high-scoring one.
    /// P: confidence = (1.0−0.542)/(1.0−0.542) = 1.0.
    /// S: MeanScore ≈ 0.655, confidence = (0.655−0.542)/(1.0−0.542) ≈ 0.246.
    /// Source: Campen et al. (2008) TOP-IDP scale; normalized distance from cutoff.
    /// </summary>
    [Test]
    public void CalculateConfidence_LowerScoreRegion_LowerConfidence()
    {
        // 30×P has the highest possible score (P = 0.987 on TOP-IDP)
        var highScoreResult = DisorderPredictor.PredictDisorder(new string('P', 30), minRegionLength: 5);
        // 30×S has a moderate disorder score (S = 0.341 propensity, normalized ≈ 0.655)
        var lowerScoreResult = DisorderPredictor.PredictDisorder(new string('S', 30), minRegionLength: 5);

        Assert.That(highScoreResult.DisorderedRegions, Has.Count.EqualTo(1));
        Assert.That(lowerScoreResult.DisorderedRegions, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(highScoreResult.DisorderedRegions[0].Confidence,
                Is.EqualTo(1.0).Within(0.01), "P confidence = 1.0");
            Assert.That(lowerScoreResult.DisorderedRegions[0].Confidence,
                Is.EqualTo(0.246).Within(0.02),
                "S confidence = (0.655−0.542)/(1.0−0.542) ≈ 0.246");
            Assert.That(lowerScoreResult.DisorderedRegions[0].Confidence,
                Is.LessThan(highScoreResult.DisorderedRegions[0].Confidence),
                "Lower mean disorder score produces lower confidence");
        });
    }

    #endregion
}
