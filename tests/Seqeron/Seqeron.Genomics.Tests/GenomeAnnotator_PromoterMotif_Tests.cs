using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for GenomeAnnotator.FindPromoterMotifs - bacterial promoter motif detection.
/// Test Unit: ANNOT-PROM-001
///
/// References:
/// - Wikipedia: Promoter (genetics) - bacterial promoter structure
/// - Wikipedia: Pribnow box - -10 element consensus TATAAT
/// - Pribnow (1975) PNAS 72(3):784-788
/// - Harley & Reynolds (1987) NAR 15(5):2343-2361
/// </summary>
[TestFixture]
public class GenomeAnnotator_PromoterMotif_Tests
{
    #region Consensus Motif Detection

    /// <summary>
    /// M01: Full -35 box consensus (TTGACA) should be detected with score 1.0.
    /// Evidence: Wikipedia states -35 box consensus is TTGACA (6 bp).
    /// </summary>
    [Test]
    public void FindPromoterMotifs_FullMinus35Consensus_ReturnsCorrectHit()
    {
        // Arrange - place TTGACA at known position (after 5 G's)
        const string sequence = "GGGGGTTGACAGGGGG";
        const int expectedPosition = 5;

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        var minus35Hits = motifs.Where(m => m.type == "-35 box" && m.sequence == "TTGACA").ToList();
        Assert.That(minus35Hits.Count, Is.GreaterThanOrEqualTo(1), "Should find at least one full -35 box");

        var fullHit = minus35Hits.First(h => h.sequence == "TTGACA");
        Assert.Multiple(() =>
        {
            Assert.That(fullHit.position, Is.EqualTo(expectedPosition), "Position should be 0-based index");
            Assert.That(fullHit.type, Is.EqualTo("-35 box"), "Type should be -35 box");
            Assert.That(fullHit.score, Is.EqualTo(1.0).Within(0.001), "Full 6 bp match should have score 1.0");
        });
    }

    /// <summary>
    /// M02: Full -10 box consensus (TATAAT) should be detected with score 1.0.
    /// Evidence: Wikipedia/Pribnow box consensus is TATAAT (6 bp).
    /// </summary>
    [Test]
    public void FindPromoterMotifs_FullMinus10Consensus_ReturnsCorrectHit()
    {
        // Arrange - place TATAAT at known position (after 5 C's)
        const string sequence = "CCCCCTATAAT CCCCC";
        const int expectedPosition = 5;

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        var minus10Hits = motifs.Where(m => m.type == "-10 box" && m.sequence == "TATAAT").ToList();
        Assert.That(minus10Hits.Count, Is.GreaterThanOrEqualTo(1), "Should find at least one full -10 box");

        var fullHit = minus10Hits.First(h => h.sequence == "TATAAT");
        Assert.Multiple(() =>
        {
            Assert.That(fullHit.position, Is.EqualTo(expectedPosition), "Position should be 0-based index");
            Assert.That(fullHit.type, Is.EqualTo("-10 box"), "Type should be -10 box");
            Assert.That(fullHit.score, Is.EqualTo(1.0).Within(0.001), "Full 6 bp match should have score 1.0");
        });
    }

    /// <summary>
    /// M09: Sequence containing both -35 and -10 boxes should return hits for both types.
    /// Evidence: Both elements are part of bacterial promoter architecture.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_BothMotifTypes_ReturnsBothTypes()
    {
        // Arrange - typical promoter-like sequence with ~17 bp spacing
        const string sequence = "GGGGGTTGACAGGGGGGGGGGGGGGGGGTATAAT GGGGG";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(motifs.Any(m => m.type == "-35 box"), Is.True, "Should find -35 box");
            Assert.That(motifs.Any(m => m.type == "-10 box"), Is.True, "Should find -10 box");
        });
    }

    #endregion

    #region Partial Motif Detection

    /// <summary>
    /// M03: Partial -35 box (5 bp) should return score less than 1.0.
    /// Score = sum of matched position probabilities / total.
    /// TTGAC (positions 1–5): (0.69+0.79+0.61+0.56+0.54) / 3.73 ≈ 0.855.
    /// Source: Wikipedia "Promoter (genetics)" nucleotide occurrence probabilities.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_PartialMinus35_ReturnsLowerScore()
    {
        // Arrange - TTGAC is a 5 bp partial match (positions 1–5 of consensus)
        const string sequence = "GGGGGTTGACGGGGG";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        var partialHit = motifs.FirstOrDefault(m => m.type == "-35 box" && m.sequence == "TTGAC");
        Assert.That(partialHit, Is.Not.EqualTo(default((int, string, string, double))), "Should find partial -35 box TTGAC");
        Assert.That(partialHit.score, Is.LessThan(1.0), "5 bp match should have score < 1.0");
        Assert.That(partialHit.score, Is.EqualTo(0.855).Within(0.001), "Score = (69+79+61+56+54)/373");
    }

    /// <summary>
    /// M04: Partial -10 box (5 bp) should return score less than 1.0.
    /// Score = sum of matched position probabilities / total.
    /// TATAA (positions 1–5): (0.77+0.76+0.60+0.61+0.56) / 4.12 ≈ 0.801.
    /// Source: Wikipedia "Promoter (genetics)" nucleotide occurrence probabilities.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_PartialMinus10_ReturnsLowerScore()
    {
        // Arrange - TATAA is a 5 bp partial match (positions 1–5 of consensus)
        const string sequence = "CCCCCTATAA CCCCC";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        var partialHit = motifs.FirstOrDefault(m => m.type == "-10 box" && m.sequence == "TATAA");
        Assert.That(partialHit, Is.Not.EqualTo(default((int, string, string, double))), "Should find partial -10 box TATAA");
        Assert.That(partialHit.score, Is.LessThan(1.0), "5 bp match should have score < 1.0");
        Assert.That(partialHit.score, Is.EqualTo(0.801).Within(0.001), "Score = (77+76+60+61+56)/412");
    }

    /// <summary>
    /// S03: All -35 box variants detected with probability-weighted scores.
    /// Variants are consensus substrings: full(6bp), prefix-5bp, suffix-5bp, prefix-4bp.
    /// Scores from Wikipedia "Promoter (genetics)" / Harley &amp; Reynolds (1987).
    /// -35 box: T(69%) T(79%) G(61%) A(56%) C(54%) A(54%), total weight 3.73.
    /// </summary>
    [TestCase("TTGACA", 1.000)]
    [TestCase("TTGAC", 0.855)]
    [TestCase("TGACA", 0.815)]
    [TestCase("TTGA", 0.710)]
    public void FindPromoterMotifs_AllMinus35Variants_Detected(string variant, double expectedScore)
    {
        // Arrange
        string sequence = "CCCCC" + variant + "CCCCC";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        var hit = motifs.FirstOrDefault(m => m.type == "-35 box" && m.sequence == variant);
        Assert.That(hit, Is.Not.EqualTo(default((int, string, string, double))), $"Should find -35 variant {variant}");
        Assert.That(hit.score, Is.EqualTo(expectedScore).Within(0.001));
    }

    /// <summary>
    /// S04: All -10 box variants detected with probability-weighted scores.
    /// Variants are consensus substrings: full(6bp), prefix-5bp, suffix-5bp, prefix-4bp.
    /// Scores from Wikipedia "Promoter (genetics)" / Harley &amp; Reynolds (1987).
    /// -10 box: T(77%) A(76%) T(60%) A(61%) A(56%) T(82%), total weight 4.12.
    /// </summary>
    [TestCase("TATAAT", 1.000)]
    [TestCase("TATAA", 0.801)]
    [TestCase("ATAAT", 0.813)]
    [TestCase("TATA", 0.665)]
    public void FindPromoterMotifs_AllMinus10Variants_Detected(string variant, double expectedScore)
    {
        // Arrange
        string sequence = "CCCCC" + variant + "CCCCC";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        var hit = motifs.FirstOrDefault(m => m.type == "-10 box" && m.sequence == variant);
        Assert.That(hit, Is.Not.EqualTo(default((int, string, string, double))), $"Should find -10 variant {variant}");
        Assert.That(hit.score, Is.EqualTo(expectedScore).Within(0.001));
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// M05: Sequence without any promoter motifs should return empty collection.
    /// Using only C nucleotides to avoid any motif matches.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_NoMotifs_ReturnsEmpty()
    {
        // Arrange
        const string sequence = "CCCCCCCCCCCCCCCCCCCC";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        Assert.That(motifs, Is.Empty, "Sequence without T/A should have no promoter motifs");
    }

    /// <summary>
    /// M06: Empty sequence should return empty collection without error.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_EmptySequence_ReturnsEmpty()
    {
        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(string.Empty).ToList();

        // Assert
        Assert.That(motifs, Is.Empty);
    }

    /// <summary>
    /// M07: Mixed case input should be handled correctly (case-insensitive).
    /// Implementation uses ToUpperInvariant() internally.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_MixedCase_HandlesCorrectly()
    {
        // Arrange - lowercase -35 box
        const string sequence = "gggggttgacaggggg";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        Assert.That(motifs.Any(m => m.type == "-35 box"), Is.True,
            "Should find motif regardless of case");
    }

    /// <summary>
    /// M08: Multiple occurrences of the same motif type should all be reported.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_MultipleMotifsOfSameType_ReturnsAll()
    {
        // Arrange - two TTGACA motifs
        const string sequence = "TTGACACCCCCTTGACA";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        var fullMinus35Hits = motifs.Where(m => m.type == "-35 box" && m.sequence == "TTGACA").ToList();
        Assert.That(fullMinus35Hits.Count, Is.EqualTo(2), "Should find both full -35 box occurrences");
        Assert.That(fullMinus35Hits[0].position, Is.EqualTo(0), "First occurrence at position 0");
        Assert.That(fullMinus35Hits[1].position, Is.EqualTo(11), "Second occurrence at position 11");
    }

    #endregion

    #region Position and Score Validation

    /// <summary>
    /// M10: Position should be 0-based index of motif start in the sequence.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_CorrectPositionReporting()
    {
        // Arrange - motif at exact position 10
        const string sequence = "CCCCCCCCCCTTGACACCCC";
        const int expectedPosition = 10;

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert
        var hit = motifs.First(m => m.sequence == "TTGACA");
        Assert.That(hit.position, Is.EqualTo(expectedPosition));
    }

    /// <summary>
    /// S02: Score reflects position-specific probability weights from literature.
    /// -35 box: T(69%) T(79%) G(61%) A(56%) C(54%) A(54%), total 3.73.
    /// Score = sum(matched position probabilities) / 3.73.
    /// Source: Wikipedia "Promoter (genetics)" / Harley &amp; Reynolds (1987).
    /// </summary>
    [Test]
    public void FindPromoterMotifs_Score_ReflectsPositionProbabilityWeights()
    {
        // Arrange - sequence with 4, 5, and 6 bp -35 variants
        const string sequence = "CCCTTGACCCTTGACCCCTTGACACCC";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert - verify probability-weighted scores
        Assert.Multiple(() =>
        {
            var ttga = motifs.FirstOrDefault(m => m.sequence == "TTGA");
            var ttgac = motifs.FirstOrDefault(m => m.sequence == "TTGAC");
            var ttgaca = motifs.FirstOrDefault(m => m.sequence == "TTGACA");

            if (ttga != default)
                Assert.That(ttga.score, Is.EqualTo(0.710).Within(0.001), "TTGA: (69+79+61+56)/373");
            if (ttgac != default)
                Assert.That(ttgac.score, Is.EqualTo(0.855).Within(0.001), "TTGAC: (69+79+61+56+54)/373");
            if (ttgaca != default)
                Assert.That(ttgaca.score, Is.EqualTo(1.0).Within(0.001), "TTGACA: full consensus");
        });
    }

    #endregion

    #region Adjacent and Overlapping Motifs

    /// <summary>
    /// S01: Adjacent/overlapping motifs should all be reported.
    /// When partial motifs overlap, implementation reports all matches.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_AdjacentMotifs_ReportsAllPositions()
    {
        // Arrange - TATA appears at position 0, and again within TATAAT at position 0
        const string sequence = "TATAAT";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(sequence).ToList();

        // Assert - should find multiple overlapping matches
        var tataHits = motifs.Where(m => m.sequence == "TATA").ToList();
        var tataatHits = motifs.Where(m => m.sequence == "TATAAT").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tataatHits.Count, Is.GreaterThanOrEqualTo(1), "Should find full TATAAT");
            Assert.That(tataHits.Count, Is.GreaterThanOrEqualTo(1), "Should find partial TATA");
        });
    }

    #endregion

    #region Real-World Sequence (Could Test)

    /// <summary>
    /// C01: Test with a realistic E. coli promoter-like sequence.
    /// The lac promoter contains recognizable -35 and -10 elements.
    /// Note: This is a simplified test; real promoters may have mismatches.
    /// </summary>
    [Test]
    public void FindPromoterMotifs_RealisticPromoterSequence_FindsMotifs()
    {
        // Arrange - simplified E. coli lac promoter region (stylized)
        // Real lac promoter: -35 at TTTACA (not exact), -10 at TATGTT (not exact)
        // Using consensus sequences for test reliability
        const string lacPromoterStyled = "AAAATTGACACCCCCCCCCCCCCCCCCTATAAT AAA";

        // Act
        var motifs = GenomeAnnotator.FindPromoterMotifs(lacPromoterStyled).ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(motifs.Any(m => m.type == "-35 box"), Is.True, "Should find -35 box region");
            Assert.That(motifs.Any(m => m.type == "-10 box"), Is.True, "Should find -10 box region");
        });
    }

    #endregion
}
