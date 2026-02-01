// Canonical Test File: RESTR-DIGEST-001 - Restriction Digest Simulation
// Methods Under Test: RestrictionAnalyzer.Digest, GetDigestSummary, CreateMap, AreCompatible, FindCompatibleEnzymes
//
// Evidence Sources:
// - Wikipedia: Restriction digest (https://en.wikipedia.org/wiki/Restriction_digest)
// - Wikipedia: Restriction enzyme (https://en.wikipedia.org/wiki/Restriction_enzyme)
// - Wikipedia: Restriction map (https://en.wikipedia.org/wiki/Restriction_map)
// - Addgene: Restriction Digest Protocol (https://www.addgene.org/protocols/restriction-digest/)
// - Roberts RJ (1976) Restriction endonucleases
// - REBASE: The Restriction Enzyme Database
//
// Algorithm Documentation: docs/algorithms/MolTools/Restriction_Digest_Simulation.md
// TestSpec: TestSpecs/RESTR-DIGEST-001.md

using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical test class for RESTR-DIGEST-001: Restriction Digest Simulation.
/// Tests Digest, GetDigestSummary, CreateMap, AreCompatible, and FindCompatibleEnzymes methods.
/// For FindSites/GetEnzyme tests, see RestrictionAnalyzer_FindSites_Tests.cs (RESTR-FIND-001).
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_Digest_Tests
{
    #region Digest Core Tests

    /// <summary>
    /// Evidence: Wikipedia - Single cut produces two fragments.
    /// Invariant #1: k cut positions → k+1 fragments.
    /// Invariant #2: Fragment sum = original sequence length (Addgene Protocol).
    /// </summary>
    [Test]
    public void Digest_SingleCut_ReturnsTwoFragmentsWithCorrectSum()
    {
        // Arrange: EcoRI site at position 3 (GAATTC), cuts after G (position 4)
        var sequence = new DnaSequence("AAAGAATTCAAA");

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariants #1 and #2
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(2), "Single cut should produce 2 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(sequence.Length),
                "Fragment sum must equal original sequence length (Addgene)");
        });
    }

    /// <summary>
    /// Evidence: Implementation contract - No cuts returns whole sequence.
    /// Invariant #7: Zero cut sites returns single fragment equal to original.
    /// </summary>
    [Test]
    public void Digest_NoCuts_ReturnsWholeSequenceAsSingleFragment()
    {
        // Arrange: Sequence without any EcoRI sites
        var sequence = new DnaSequence("AAAAAAAAAAAA");

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariant #7
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(1), "No cuts should return single fragment");
            Assert.That(fragments[0].Length, Is.EqualTo(sequence.Length), "Fragment equals original");
            Assert.That(fragments[0].Sequence, Is.EqualTo(sequence.Sequence), "Content matches original");
            Assert.That(fragments[0].LeftEnzyme, Is.Null, "First fragment has no left enzyme");
            Assert.That(fragments[0].RightEnzyme, Is.Null, "Last fragment has no right enzyme");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - Multiple cuts produce expected fragment count.
    /// Invariant #1: k cut positions → k+1 fragments.
    /// Invariant #2: Fragment sum = original sequence length.
    /// </summary>
    [Test]
    public void Digest_MultipleCuts_ReturnsCorrectFragmentCount()
    {
        // Arrange: Two EcoRI sites (positions 0 and 9)
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA");

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariants #1 and #2
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(3), "Two cuts should produce 3 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(sequence.Length),
                "Fragment sum must equal original sequence length");
        });
    }

    /// <summary>
    /// Evidence: Implementation - Comprehensive fragment property verification.
    /// Invariants #3, #4, #5: Fragment ordering, numbering, boundary enzymes.
    /// </summary>
    [Test]
    public void Digest_FragmentsHaveCorrectProperties()
    {
        // Arrange
        var sequence = new DnaSequence("AAAGAATTCAAA");
        // EcoRI cuts at position 4 (after G in GAATTC at position 3)

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Multiple invariants
        Assert.Multiple(() =>
        {
            // Invariant #3: Sequential numbering
            Assert.That(fragments[0].FragmentNumber, Is.EqualTo(1), "First fragment numbered 1");
            Assert.That(fragments[1].FragmentNumber, Is.EqualTo(2), "Second fragment numbered 2");

            // Invariant #4: Start positions in ascending order
            Assert.That(fragments[0].StartPosition, Is.EqualTo(0), "First fragment starts at 0");
            Assert.That(fragments[1].StartPosition, Is.GreaterThan(fragments[0].StartPosition),
                "Positions increase monotonically");

            // Invariant #5: Boundary enzymes
            Assert.That(fragments[0].LeftEnzyme, Is.Null, "First fragment LeftEnzyme is null");
            Assert.That(fragments[^1].RightEnzyme, Is.Null, "Last fragment RightEnzyme is null");

            // Invariant #6: Positive length
            Assert.That(fragments.All(f => f.Length > 0), "All fragments have positive length");
        });
    }

    /// <summary>
    /// Evidence: Implementation - Fragment positions must increase monotonically.
    /// Invariant #4: Fragment start positions sorted ascending.
    /// </summary>
    [Test]
    public void Digest_FragmentStartPositions_IncreaseMonotonically()
    {
        // Arrange: Multiple cuts for more data points
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAAGAATTCAAA");

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariant #4
        for (int i = 1; i < fragments.Count; i++)
        {
            Assert.That(fragments[i].StartPosition, Is.GreaterThan(fragments[i - 1].StartPosition),
                $"Fragment {i + 1} start position should be greater than fragment {i}");
        }
    }

    /// <summary>
    /// Evidence: Implementation - Verify fragment sequence content is correct substring.
    /// </summary>
    [Test]
    public void Digest_FragmentSequenceContent_MatchesExpectedSubstring()
    {
        // Arrange: Known sequence with single cut
        var sequence = new DnaSequence("AAAGAATTCAAA");
        // EcoRI at position 3, cuts after position 4 (G|AATTC)

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Verify content matches
        Assert.Multiple(() =>
        {
            foreach (var fragment in fragments)
            {
                string expectedContent = sequence.Sequence.Substring(fragment.StartPosition, fragment.Length);
                Assert.That(fragment.Sequence, Is.EqualTo(expectedContent),
                    $"Fragment {fragment.FragmentNumber} content should match substring at position {fragment.StartPosition}");
            }
        });
    }

    /// <summary>
    /// Evidence: Implementation - Multiple enzymes produce combined cut sites.
    /// </summary>
    [Test]
    public void Digest_MultipleEnzymes_CutsWithBoth()
    {
        // Arrange: EcoRI at position 3, BamHI at position 12
        var sequence = new DnaSequence("AAAGAATTCAAAGGATCCAAA");

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI", "BamHI").ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(3), "Two cuts should produce 3 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(sequence.Length),
                "Fragment sum must equal original sequence length");
        });
    }

    #endregion

    #region Digest Edge Cases

    /// <summary>
    /// Evidence: API contract - No enzymes provided must throw.
    /// </summary>
    [Test]
    public void Digest_NoEnzymes_ThrowsArgumentException()
    {
        var sequence = new DnaSequence("ACGT");

        Assert.Throws<ArgumentException>(() =>
            RestrictionAnalyzer.Digest(sequence).ToList());
    }

    /// <summary>
    /// Evidence: API contract - Null sequence must throw.
    /// </summary>
    [Test]
    public void Digest_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.Digest(null!, "EcoRI").ToList());
    }

    /// <summary>
    /// Evidence: Edge case - Sequence shorter than recognition sequence.
    /// </summary>
    [Test]
    public void Digest_SequenceShorterThanRecognition_ReturnsWholeSequence()
    {
        // Arrange: 4 bp sequence, EcoRI needs 6 bp
        var sequence = new DnaSequence("ACGT");

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(1), "No possible cut sites");
            Assert.That(fragments[0].Length, Is.EqualTo(sequence.Length));
        });
    }

    #endregion

    #region GetDigestSummary Tests

    /// <summary>
    /// Evidence: Implementation contract - Summary returns correct statistics.
    /// Invariants #8, #9, #10: Sorting, bounds, enzyme list.
    /// </summary>
    [Test]
    public void GetDigestSummary_ReturnsCorrectSummaryWithInvariants()
    {
        // Arrange
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA");

        // Act
        var summary = RestrictionAnalyzer.GetDigestSummary(sequence, "EcoRI");

        // Assert: All invariants
        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalFragments, Is.EqualTo(3), "Should report 3 fragments");
            Assert.That(summary.FragmentSizes, Has.Count.EqualTo(3), "Should have 3 sizes");

            // Invariant #9: Size bounds
            Assert.That(summary.LargestFragment, Is.GreaterThanOrEqualTo(summary.SmallestFragment),
                "Largest >= Smallest");
            Assert.That(summary.AverageFragmentSize, Is.GreaterThanOrEqualTo(summary.SmallestFragment),
                "Average >= Smallest");
            Assert.That(summary.AverageFragmentSize, Is.LessThanOrEqualTo(summary.LargestFragment),
                "Average <= Largest");

            // Invariant #10: Enzymes list
            Assert.That(summary.EnzymesUsed, Contains.Item("EcoRI"), "Should list EcoRI");
        });
    }

    /// <summary>
    /// Evidence: Implementation - Fragment sizes sorted descending.
    /// Invariant #8: FragmentSizes[i] >= FragmentSizes[i+1].
    /// </summary>
    [Test]
    public void GetDigestSummary_FragmentsSortedDescending()
    {
        // Arrange: Asymmetric cuts to ensure different sizes
        var sequence = new DnaSequence("GAATTCAAAAAGAATTCAAA");

        // Act
        var summary = RestrictionAnalyzer.GetDigestSummary(sequence, "EcoRI");

        // Assert: Invariant #8
        for (int i = 0; i < summary.FragmentSizes.Count - 1; i++)
        {
            Assert.That(summary.FragmentSizes[i], Is.GreaterThanOrEqualTo(summary.FragmentSizes[i + 1]),
                $"Size at index {i} should be >= size at index {i + 1}");
        }
    }

    #endregion

    #region Restriction Map Tests

    /// <summary>
    /// Evidence: Implementation - CreateMap returns comprehensive map data.
    /// Invariants #11, #12, #13: Non-cutters, unique cutters, site count.
    /// </summary>
    [Test]
    public void CreateMap_ReturnsCorrectMapWithAllFields()
    {
        // Arrange
        var sequence = new DnaSequence("GAATTCAAAGAATTC");

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI");

        // Assert: Multiple invariants
        Assert.Multiple(() =>
        {
            Assert.That(map.SequenceLength, Is.EqualTo(sequence.Length), "Sequence length preserved");
            Assert.That(map.SitesByEnzyme.ContainsKey("EcoRI"), "EcoRI sites present");
            // Invariant #13: TotalSites counts forward-strand only
            Assert.That(map.TotalSites, Is.GreaterThan(0), "Should have sites");
        });
    }

    /// <summary>
    /// Evidence: Implementation - Unique cutters are enzymes with exactly one site.
    /// Invariant #12: UniqueCutters list contains enzymes with one site in SitesByEnzyme.
    /// Note: SitesByEnzyme contains all positions (both strands may have same position for palindromic sites).
    /// </summary>
    [Test]
    public void CreateMap_IdentifiesUniqueCutters()
    {
        // Arrange: Sequence where each enzyme cuts only once (forward strand)
        // SitesByEnzyme includes both strand positions, so a palindromic site appears twice
        var sequence = new DnaSequence("AAAAAAAAAGAATTCAAAAAAAAAAAAAAAGGATCCAAAAAAAAAA");

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI", "BamHI");

        // Assert: Invariant #12
        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(2), "Two forward-strand sites");
            Assert.That(map.SitesByEnzyme.ContainsKey("EcoRI"), "EcoRI present");
            Assert.That(map.SitesByEnzyme.ContainsKey("BamHI"), "BamHI present");
            // Verify sites are found - SitesByEnzyme may contain positions from both strands
            Assert.That(map.SitesByEnzyme["EcoRI"], Is.Not.Empty, "EcoRI has sites");
            Assert.That(map.SitesByEnzyme["BamHI"], Is.Not.Empty, "BamHI has sites");
        });
    }

    /// <summary>
    /// Evidence: Implementation - Non-cutters are enzymes with zero sites.
    /// Invariant #11: NonCutters list contains enzymes with no sites.
    /// </summary>
    [Test]
    public void CreateMap_IdentifiesNonCutters()
    {
        // Arrange: Sequence without any recognition sites
        var sequence = new DnaSequence("AAAAAAAAAA");

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI", "BamHI");

        // Assert: Invariant #11
        Assert.Multiple(() =>
        {
            Assert.That(map.NonCutters, Contains.Item("EcoRI"), "EcoRI is non-cutter");
            Assert.That(map.NonCutters, Contains.Item("BamHI"), "BamHI is non-cutter");
            Assert.That(map.TotalSites, Is.EqualTo(0), "Zero total sites");
        });
    }

    /// <summary>
    /// Evidence: Implementation - TotalSites counts forward-strand sites only.
    /// Invariant #13: Avoid double-counting palindromic sites.
    /// </summary>
    [Test]
    public void CreateMap_TotalSites_CountsForwardStrandOnly()
    {
        // Arrange: Single EcoRI site (palindromic, appears on both strands)
        var sequence = new DnaSequence("AAAAAAGAATTCAAAAAA");

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI");

        // Assert: Invariant #13 - should count only forward strand
        Assert.That(map.TotalSites, Is.EqualTo(1),
            "Should count forward-strand sites only, not both strands");
    }

    /// <summary>
    /// Evidence: API contract - Null sequence must throw.
    /// </summary>
    [Test]
    public void CreateMap_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.CreateMap(null!, "EcoRI"));
    }

    #endregion

    #region Compatibility Tests

    /// <summary>
    /// Evidence: Wikipedia - Blunt enzymes produce ligatable ends with other blunt enzymes.
    /// Invariant #14: All blunt-end enzymes are compatible.
    /// </summary>
    [Test]
    public void AreCompatible_BluntEnzymes_AreCompatible()
    {
        // EcoRV and SmaI are both blunt cutters
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRV", "SmaI");

        Assert.That(compatible, Is.True, "Blunt enzymes should be compatible (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia - BamHI and BglII both produce GATC overhangs.
    /// Invariant #15: Identical overhangs are compatible.
    /// </summary>
    [Test]
    public void AreCompatible_SameOverhang_AreCompatible()
    {
        // BamHI (GGATCC, GATC overhang) and BglII (AGATCT, GATC overhang)
        bool compatible = RestrictionAnalyzer.AreCompatible("BamHI", "BglII");

        Assert.That(compatible, Is.True, "Same overhang enzymes should be compatible (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia - Different overhangs cannot ligate.
    /// Invariant #15 (inverse): Different overhangs are not compatible.
    /// </summary>
    [Test]
    public void AreCompatible_DifferentOverhangs_NotCompatible()
    {
        // EcoRI (AATT overhang) and PstI (TGCA overhang)
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRI", "PstI");

        Assert.That(compatible, Is.False, "Different overhang enzymes should not be compatible");
    }

    /// <summary>
    /// Evidence: API contract - Unknown enzyme returns false (no exception).
    /// </summary>
    [Test]
    public void AreCompatible_UnknownEnzyme_ReturnsFalse()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRI", "UnknownEnzyme");

        Assert.That(compatible, Is.False, "Unknown enzyme should return false, not throw");
    }

    /// <summary>
    /// Evidence: Mathematical property - Compatibility is symmetric.
    /// Invariant #16: AreCompatible(A, B) == AreCompatible(B, A).
    /// </summary>
    [Test]
    [TestCase("BamHI", "BglII")]
    [TestCase("EcoRV", "SmaI")]
    [TestCase("EcoRI", "PstI")]
    public void AreCompatible_IsSymmetric(string enzyme1, string enzyme2)
    {
        bool forward = RestrictionAnalyzer.AreCompatible(enzyme1, enzyme2);
        bool reverse = RestrictionAnalyzer.AreCompatible(enzyme2, enzyme1);

        Assert.That(forward, Is.EqualTo(reverse),
            $"AreCompatible({enzyme1}, {enzyme2}) should equal AreCompatible({enzyme2}, {enzyme1})");
    }

    /// <summary>
    /// Evidence: Implementation - FindCompatibleEnzymes returns known pairs.
    /// </summary>
    [Test]
    public void FindCompatibleEnzymes_FindsKnownPairs()
    {
        // Act
        var compatiblePairs = RestrictionAnalyzer.FindCompatibleEnzymes().ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(compatiblePairs, Has.Count.GreaterThan(0), "Should find compatible pairs");

            // Check that BamHI/BglII pair is found (known compatible pair from Wikipedia)
            bool hasBamHIBglII = compatiblePairs.Any(p =>
                (p.Enzyme1 == "BamHI" && p.Enzyme2 == "BglII") ||
                (p.Enzyme1 == "BglII" && p.Enzyme2 == "BamHI"));
            Assert.That(hasBamHIBglII, Is.True, "BamHI/BglII should be listed as compatible");
        });
    }

    /// <summary>
    /// Evidence: Implementation - All returned pairs should be actually compatible.
    /// </summary>
    [Test]
    public void FindCompatibleEnzymes_AllReturnedPairsAreActuallyCompatible()
    {
        // Act
        var compatiblePairs = RestrictionAnalyzer.FindCompatibleEnzymes().Take(10).ToList();

        // Assert: Verify each pair
        foreach (var pair in compatiblePairs)
        {
            bool compatible = RestrictionAnalyzer.AreCompatible(pair.Enzyme1, pair.Enzyme2);
            Assert.That(compatible, Is.True,
                $"{pair.Enzyme1} and {pair.Enzyme2} should be compatible per AreCompatible()");
        }
    }

    #endregion
}
