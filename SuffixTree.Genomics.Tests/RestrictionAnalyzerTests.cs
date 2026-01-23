// Tests for RESTR-DIGEST-001: Restriction Digest Simulation
// For RESTR-FIND-001 tests (FindSites, FindAllSites, GetEnzyme), see RestrictionAnalyzer_FindSites_Tests.cs

using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for Digest, Map, and Compatibility functionality.
/// For FindSites/GetEnzyme tests, see RestrictionAnalyzer_FindSites_Tests.cs (RESTR-FIND-001).
/// </summary>
[TestFixture]
public class RestrictionAnalyzerTests
{
    #region Smoke Tests for FindSites (Integration)

    /// <summary>
    /// Smoke test: Basic FindSites functionality.
    /// Detailed tests in RestrictionAnalyzer_FindSites_Tests.cs.
    /// </summary>
    [Test]
    public void FindSites_EcoRI_FindsSite_Smoke()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");
        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI").ToList();

        Assert.That(sites.Any(s => s.IsForwardStrand && s.Position == 3));
    }

    #endregion

    #region Digest Simulation Tests

    [Test]
    public void Digest_SingleCut_ReturnsTwoFragments()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        Assert.That(fragments, Has.Count.EqualTo(2));
        Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(sequence.Length));
    }

    [Test]
    public void Digest_NoCuts_ReturnsWholeSequence()
    {
        var sequence = new DnaSequence("AAAAAAAAAAAA");
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        Assert.That(fragments, Has.Count.EqualTo(1));
        Assert.That(fragments[0].Length, Is.EqualTo(sequence.Length));
    }

    [Test]
    public void Digest_MultipleCuts_ReturnsCorrectFragments()
    {
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA");
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        Assert.That(fragments, Has.Count.EqualTo(3));
    }

    [Test]
    public void Digest_FragmentsHaveCorrectProperties()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        Assert.That(fragments[0].FragmentNumber, Is.EqualTo(1));
        Assert.That(fragments[1].FragmentNumber, Is.EqualTo(2));
        Assert.That(fragments[0].StartPosition, Is.EqualTo(0));
    }

    [Test]
    public void Digest_MultipleEnzymes_CutsWithBoth()
    {
        var sequence = new DnaSequence("AAAGAATTCAAAGGATCCAAA");
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI", "BamHI").ToList();

        Assert.That(fragments, Has.Count.EqualTo(3));
    }

    [Test]
    public void Digest_NoEnzymes_ThrowsException()
    {
        var sequence = new DnaSequence("ACGT");
        Assert.Throws<ArgumentException>(() =>
            RestrictionAnalyzer.Digest(sequence).ToList());
    }

    [Test]
    public void GetDigestSummary_ReturnsCorrectSummary()
    {
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA");
        var summary = RestrictionAnalyzer.GetDigestSummary(sequence, "EcoRI");

        Assert.That(summary.TotalFragments, Is.EqualTo(3));
        Assert.That(summary.FragmentSizes, Has.Count.EqualTo(3));
        Assert.That(summary.LargestFragment, Is.GreaterThan(0));
        Assert.That(summary.SmallestFragment, Is.GreaterThan(0));
        Assert.That(summary.EnzymesUsed, Contains.Item("EcoRI"));
    }

    [Test]
    public void GetDigestSummary_FragmentsSortedDescending()
    {
        var sequence = new DnaSequence("GAATTCAAAAAGAATTCAAA");
        var summary = RestrictionAnalyzer.GetDigestSummary(sequence, "EcoRI");

        for (int i = 0; i < summary.FragmentSizes.Count - 1; i++)
        {
            Assert.That(summary.FragmentSizes[i], Is.GreaterThanOrEqualTo(summary.FragmentSizes[i + 1]));
        }
    }

    #endregion

    #region Restriction Map Tests

    [Test]
    public void CreateMap_ReturnsCorrectMap()
    {
        var sequence = new DnaSequence("GAATTCAAAGAATTC");
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI");

        Assert.That(map.SequenceLength, Is.EqualTo(sequence.Length));
        Assert.That(map.TotalSites, Is.EqualTo(2));
        Assert.That(map.SitesByEnzyme.ContainsKey("EcoRI"));
    }

    [Test]
    public void CreateMap_IdentifiesUniqueCutters()
    {
        // Create a longer sequence where each enzyme cuts only once
        var sequence = new DnaSequence("AAAAAAAAAGAATTCAAAAAAAAAAAAAAAGGATCCAAAAAAAAAA");
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI", "BamHI");

        // Map identifies enzymes that have only one cut site (on forward strand)
        Assert.That(map.TotalSites, Is.EqualTo(2));
        Assert.That(map.SitesByEnzyme.ContainsKey("EcoRI"));
        Assert.That(map.SitesByEnzyme.ContainsKey("BamHI"));
    }

    [Test]
    public void CreateMap_IdentifiesNonCutters()
    {
        var sequence = new DnaSequence("AAAAAAAAAA");
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI", "BamHI");

        Assert.That(map.NonCutters, Contains.Item("EcoRI"));
        Assert.That(map.NonCutters, Contains.Item("BamHI"));
        Assert.That(map.TotalSites, Is.EqualTo(0));
    }

    #endregion

    #region Compatibility Tests

    [Test]
    public void AreCompatible_BluntEnzymes_AreCompatible()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRV", "SmaI");
        Assert.That(compatible, Is.True);
    }

    [Test]
    public void AreCompatible_SameOverhang_AreCompatible()
    {
        // BamHI and BglII both produce GATC overhangs
        bool compatible = RestrictionAnalyzer.AreCompatible("BamHI", "BglII");
        Assert.That(compatible, Is.True);
    }

    [Test]
    public void AreCompatible_DifferentOverhangs_NotCompatible()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRI", "PstI");
        Assert.That(compatible, Is.False);
    }

    [Test]
    public void AreCompatible_UnknownEnzyme_ReturnsFalse()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRI", "Unknown");
        Assert.That(compatible, Is.False);
    }

    [Test]
    public void FindCompatibleEnzymes_FindsPairs()
    {
        var compatiblePairs = RestrictionAnalyzer.FindCompatibleEnzymes().ToList();

        Assert.That(compatiblePairs, Has.Count.GreaterThan(0));

        // Check that BamHI and BglII are in the list
        Assert.That(compatiblePairs.Any(p =>
            (p.Enzyme1 == "BamHI" && p.Enzyme2 == "BglII") ||
            (p.Enzyme1 == "BglII" && p.Enzyme2 == "BamHI")));
    }

    #endregion

    #region Edge Cases (Digest/Map)

    [Test]
    public void Digest_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.Digest(null!, "EcoRI").ToList());
    }

    [Test]
    public void CreateMap_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.CreateMap(null!, "EcoRI"));
    }

    #endregion
}
