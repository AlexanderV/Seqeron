using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// DISORDER-LC-001: Low Complexity Region Detection (SEG Algorithm)
///
/// Tests for <see cref="DisorderPredictor.PredictLowComplexityRegions"/>.
/// Algorithm: Wootton &amp; Federhen (1993) Computers &amp; Chemistry 17(2):149-163;
///            Wootton &amp; Federhen (1996) Methods Enzymol 266:554-571.
/// Default parameters: triggerWindow=12, K1=2.2 bits, K2=2.5 bits.
/// </summary>
[TestFixture]
public class DisorderPredictor_LowComplexity_Tests
{
    #region S — Smoke Tests

    [Test]
    public void S1_Homopolymer_DetectedAsLowComplexity()
    {
        // Poly-Q (26 residues): every 12-window has entropy 0 bits ≤ 2.2.
        // Expected: one LC region spanning entire sequence, typed "Q-rich".
        string polyQ = new string('Q', 26);

        var regions = DisorderPredictor.PredictLowComplexityRegions(polyQ).ToList();

        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(0));
        Assert.That(regions[0].End, Is.EqualTo(25));
        Assert.That(regions[0].Type, Is.EqualTo("Q-rich"));
    }

    [Test]
    public void S2_ComplexSequence_NoLowComplexity()
    {
        // All 20 standard amino acids repeated twice (40 AA).
        // Every 12-window contains ≥12 distinct AAs → entropy ≥ log₂(12) ≈ 3.58 bits >> 2.2.
        string complex = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY";

        var regions = DisorderPredictor.PredictLowComplexityRegions(complex).ToList();

        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void S3_DipeptideBlock_DetectedAsLowComplexity()
    {
        // 12A + 12L = 24 AA. Window at any position has at most 2 types → H ≤ log₂(2) = 1.0 bits ≤ 2.2.
        // Entire sequence triggered → one merged LC region.
        string sequence = new string('A', 12) + new string('L', 12);

        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(0));
        Assert.That(regions[0].End, Is.EqualTo(23));
    }

    #endregion

    #region C — Corner Cases

    [Test]
    public void C1_SequenceShorterThanWindow_Empty()
    {
        // 6 residues < default triggerWindow=12 → yield break.
        string seq = "QQQQQQ";

        var regions = DisorderPredictor.PredictLowComplexityRegions(seq).ToList();

        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void C2_ExactlyWindowLength_HomopolymerDetected()
    {
        // 12 Q's = exactly triggerWindow. One window at position 0.
        // Entropy = 0 ≤ 2.2 → triggers. Region (0, 11).
        string exact = new string('Q', 12);

        var regions = DisorderPredictor.PredictLowComplexityRegions(exact).ToList();

        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(0));
        Assert.That(regions[0].End, Is.EqualTo(11));
        Assert.That(regions[0].Type, Is.EqualTo("Q-rich"));
    }

    [Test]
    public void C3_EmptySequence_Empty()
    {
        var regions = DisorderPredictor.PredictLowComplexityRegions("").ToList();
        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void C4_MinLengthFiltersShortRegions()
    {
        // 12 Q's triggered → region length 12.
        // If minLength=15, the region should be filtered out.
        string seq = new string('Q', 12);

        var regions = DisorderPredictor.PredictLowComplexityRegions(
            seq, minLength: 15).ToList();

        Assert.That(regions, Is.Empty);
    }

    #endregion

    #region M — Method Tests

    [Test]
    public void M1_EntropyBasedTrigger_FourTypesTriggersK1()
    {
        // 4 AA types × 3 each = 12 residues, H = log₂(4) = 2.0 bits ≤ 2.2 → triggers.
        string fourTypes = "AAABBBCCCDDD";

        var regions = DisorderPredictor.PredictLowComplexityRegions(fourTypes).ToList();

        Assert.That(regions, Has.Count.EqualTo(1));
    }

    [Test]
    public void M2_HighEntropy_TwelveDistinctDoesNotTrigger()
    {
        // 12 distinct AAs → H = log₂(12) ≈ 3.58 bits >> 2.2 → no trigger.
        string highEntropy = "ACDEFGHIKLMN";

        var regions = DisorderPredictor.PredictLowComplexityRegions(highEntropy).ToList();

        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void M3_TwoSeparatedLCRegions_BothDetected()
    {
        // Poly-Q + long complex separator (3× all-20 AA = 60 chars) + Poly-A.
        // The separator has entropy log₂(12+) in every 12-window.
        string polyQ = new string('Q', 20);
        string separator = string.Concat(Enumerable.Repeat("ACDEFGHIKLMNPQRSTVWY", 3));
        string polyA = new string('A', 20);
        string sequence = polyQ + separator + polyA;

        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence).ToList();

        // Both Q-rich and A-rich compositions should be detected
        Assert.That(regions.Any(r => r.Type.Contains("Q")), Is.True, "Q-rich region expected");
        Assert.That(regions.Any(r => r.Type.Contains("A")), Is.True, "A-rich region expected");
    }

    [Test]
    public void M4_CustomTriggerThreshold_StricterFiltersMore()
    {
        // With K1=0.5, only homopolymers trigger (4 types × 3: H=2.0 > 0.5).
        string fourTypes = "AAABBBCCCDDD";

        var regions = DisorderPredictor.PredictLowComplexityRegions(
            fourTypes, triggerThreshold: 0.5).ToList();

        Assert.That(regions, Is.Empty);
    }

    [Test]
    public void M5_CaseInsensitive()
    {
        string lower = new string('q', 26);
        string upper = new string('Q', 26);

        var regionsLower = DisorderPredictor.PredictLowComplexityRegions(lower).ToList();
        var regionsUpper = DisorderPredictor.PredictLowComplexityRegions(upper).ToList();

        Assert.That(regionsLower, Has.Count.EqualTo(regionsUpper.Count));
        Assert.That(regionsLower[0].Start, Is.EqualTo(regionsUpper[0].Start));
        Assert.That(regionsLower[0].End, Is.EqualTo(regionsUpper[0].End));
    }

    [Test]
    public void M6_TypeClassification_SingleDominant()
    {
        // Poly-A: single dominant AA > 50% → "A-rich".
        string polyA = new string('A', 20);

        var regions = DisorderPredictor.PredictLowComplexityRegions(polyA).ToList();

        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].Type, Is.EqualTo("A-rich"));
    }

    #endregion

    #region INV — Invariants

    [TestCase("QQQQQQQQQQQQQQQQQQQQQQQQQQ")]
    [TestCase("AAABBBCCCDDDAAABBBCCCDDD")]
    [TestCase("QQQQQQQQQQQQACDEFGHIKLMNQQQQQQQQQQQQQ")]
    public void INV1_ValidBoundaries(string sequence)
    {
        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence).ToList();

        foreach (var region in regions)
        {
            Assert.That(region.Start, Is.GreaterThanOrEqualTo(0),
                "Start must be ≥ 0");
            Assert.That(region.End, Is.LessThan(sequence.Length),
                "End must be < sequence length");
            Assert.That(region.Start, Is.LessThanOrEqualTo(region.End),
                "Start must ≤ End");
            Assert.That(region.Type, Is.Not.Null.And.Not.Empty,
                "Type must be non-empty");
        }
    }

    [TestCase("QQQQQQQQQQQQQQQQQQQQQQQQQQ")]
    [TestCase("AAAAALLLLLAAAAALLLLLAAAAALLLLLL")]
    public void INV2_NoOverlappingRegions(string sequence)
    {
        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence)
            .OrderBy(r => r.Start).ToList();

        for (int i = 1; i < regions.Count; i++)
        {
            Assert.That(regions[i].Start, Is.GreaterThan(regions[i - 1].End),
                $"Region {i} overlaps region {i - 1}");
        }
    }

    #endregion
}
