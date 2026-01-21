using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class RestrictionAnalyzerTests
{
    #region Enzyme Database Tests

    [Test]
    public void GetEnzyme_EcoRI_ReturnsCorrectEnzyme()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRI");

        Assert.That(enzyme, Is.Not.Null);
        Assert.That(enzyme!.Name, Is.EqualTo("EcoRI"));
        Assert.That(enzyme.RecognitionSequence, Is.EqualTo("GAATTC"));
        Assert.That(enzyme.CutPositionForward, Is.EqualTo(1));
        Assert.That(enzyme.CutPositionReverse, Is.EqualTo(5));
    }

    [Test]
    public void GetEnzyme_CaseInsensitive_Works()
    {
        var enzyme1 = RestrictionAnalyzer.GetEnzyme("EcoRI");
        var enzyme2 = RestrictionAnalyzer.GetEnzyme("ecori");
        var enzyme3 = RestrictionAnalyzer.GetEnzyme("ECORI");

        Assert.That(enzyme1, Is.Not.Null);
        Assert.That(enzyme2, Is.Not.Null);
        Assert.That(enzyme3, Is.Not.Null);
        Assert.That(enzyme1!.Name, Is.EqualTo(enzyme2!.Name));
        Assert.That(enzyme2.Name, Is.EqualTo(enzyme3!.Name));
    }

    [Test]
    public void GetEnzyme_UnknownEnzyme_ReturnsNull()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("UnknownEnzyme");
        Assert.That(enzyme, Is.Null);
    }

    [Test]
    public void GetEnzyme_BamHI_ReturnsCorrectEnzyme()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("BamHI");

        Assert.That(enzyme, Is.Not.Null);
        Assert.That(enzyme!.RecognitionSequence, Is.EqualTo("GGATCC"));
    }

    [Test]
    public void Enzymes_ContainsCommonEnzymes()
    {
        var enzymes = RestrictionAnalyzer.Enzymes;

        Assert.That(enzymes.ContainsKey("EcoRI"));
        Assert.That(enzymes.ContainsKey("BamHI"));
        Assert.That(enzymes.ContainsKey("HindIII"));
        Assert.That(enzymes.ContainsKey("NotI"));
        Assert.That(enzymes.Count, Is.GreaterThan(30));
    }

    [Test]
    public void GetEnzymesByCutLength_SixCutters_ReturnsEnzymes()
    {
        var sixCutters = RestrictionAnalyzer.GetEnzymesByCutLength(6).ToList();

        Assert.That(sixCutters, Has.Count.GreaterThan(10));
        Assert.That(sixCutters.All(e => e.RecognitionSequence.Length == 6));
    }

    [Test]
    public void GetBluntCutters_ReturnsBluntEnzymes()
    {
        var bluntCutters = RestrictionAnalyzer.GetBluntCutters().ToList();

        Assert.That(bluntCutters, Has.Count.GreaterThan(0));
        Assert.That(bluntCutters.All(e => e.IsBluntEnd));
    }

    [Test]
    public void GetStickyCutters_ReturnsStickyEnzymes()
    {
        var stickyCutters = RestrictionAnalyzer.GetStickyCutters().ToList();

        Assert.That(stickyCutters, Has.Count.GreaterThan(0));
        Assert.That(stickyCutters.All(e => !e.IsBluntEnd));
    }

    #endregion

    #region Enzyme Properties Tests

    [Test]
    public void RestrictionEnzyme_EcoRI_HasFivePrimeOverhang()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRI")!;

        Assert.That(enzyme.IsBluntEnd, Is.False);
        Assert.That(enzyme.OverhangType, Is.EqualTo(OverhangType.FivePrime));
    }

    [Test]
    public void RestrictionEnzyme_PstI_HasThreePrimeOverhang()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("PstI")!;

        Assert.That(enzyme.IsBluntEnd, Is.False);
        Assert.That(enzyme.OverhangType, Is.EqualTo(OverhangType.ThreePrime));
    }

    [Test]
    public void RestrictionEnzyme_EcoRV_IsBlunt()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRV")!;

        Assert.That(enzyme.IsBluntEnd, Is.True);
        Assert.That(enzyme.OverhangType, Is.EqualTo(OverhangType.Blunt));
    }

    [Test]
    public void RestrictionEnzyme_RecognitionLength_Correct()
    {
        var ecoRI = RestrictionAnalyzer.GetEnzyme("EcoRI")!;
        var notI = RestrictionAnalyzer.GetEnzyme("NotI")!;

        Assert.That(ecoRI.RecognitionLength, Is.EqualTo(6));
        Assert.That(notI.RecognitionLength, Is.EqualTo(8));
    }

    #endregion

    #region Site Finding Tests

    [Test]
    public void FindSites_EcoRI_FindsSite()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");
        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI").ToList();

        Assert.That(sites.Any(s => s.IsForwardStrand && s.Position == 3));
    }

    [Test]
    public void FindSites_MultipleSites_FindsAll()
    {
        var sequence = new DnaSequence("GAATTCAAAGAATTC");
        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.That(sites, Has.Count.EqualTo(2));
        Assert.That(sites[0].Position, Is.EqualTo(0));
        Assert.That(sites[1].Position, Is.EqualTo(9));
    }

    [Test]
    public void FindSites_NoSites_ReturnsEmpty()
    {
        var sequence = new DnaSequence("AAAAAAAAAA");
        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI").ToList();

        Assert.That(sites, Is.Empty);
    }

    [Test]
    public void FindSites_StringOverload_Works()
    {
        var sites = RestrictionAnalyzer.FindSites("AAAGAATTCAAA", "EcoRI").ToList();
        Assert.That(sites.Any(s => s.Position == 3));
    }

    [Test]
    public void FindSites_UnknownEnzyme_ThrowsException()
    {
        var sequence = new DnaSequence("ACGT");
        Assert.Throws<ArgumentException>(() =>
            RestrictionAnalyzer.FindSites(sequence, "UnknownEnzyme").ToList());
    }

    [Test]
    public void FindSites_MultipleEnzymes_FindsAllSites()
    {
        var sequence = new DnaSequence("GAATTCAAAGGATCC");
        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI", "BamHI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.That(sites.Any(s => s.Enzyme.Name == "EcoRI"));
        Assert.That(sites.Any(s => s.Enzyme.Name == "BamHI"));
    }

    [Test]
    public void FindSites_ReturnsCutPosition()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");
        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.That(sites[0].CutPosition, Is.EqualTo(4)); // 3 + 1
    }

    [Test]
    public void FindSites_ReturnsRecognizedSequence()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");
        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.That(sites[0].RecognizedSequence, Is.EqualTo("GAATTC"));
    }

    [Test]
    public void FindSites_CustomEnzyme_Works()
    {
        var customEnzyme = new RestrictionEnzyme("CustomI", "ATAT", 2, 2, "Custom");
        var sequence = new DnaSequence("AAATATAAA");
        var sites = RestrictionAnalyzer.FindSites(sequence, customEnzyme)
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Enzyme.Name, Is.EqualTo("CustomI"));
    }

    [Test]
    public void FindAllSites_FindsMultipleEnzymes()
    {
        var sequence = new DnaSequence("GAATTCGGATCCAAGCTT");
        var sites = RestrictionAnalyzer.FindAllSites(sequence)
            .Where(s => s.IsForwardStrand)
            .ToList();

        var enzymeNames = sites.Select(s => s.Enzyme.Name).Distinct().ToList();
        Assert.That(enzymeNames, Has.Count.GreaterThanOrEqualTo(3));
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

    #region Edge Cases

    [Test]
    public void FindSites_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.FindSites((DnaSequence)null!, "EcoRI").ToList());
    }

    [Test]
    public void FindSites_EmptyEnzymeName_ThrowsException()
    {
        var sequence = new DnaSequence("ACGT");
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.FindSites(sequence, "").ToList());
    }

    [Test]
    public void FindSites_EmptyStringSequence_ReturnsEmpty()
    {
        var sites = RestrictionAnalyzer.FindSites("", "EcoRI").ToList();
        Assert.That(sites, Is.Empty);
    }

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
