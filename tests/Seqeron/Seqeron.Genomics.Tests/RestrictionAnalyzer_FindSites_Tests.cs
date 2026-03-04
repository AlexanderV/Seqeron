// Canonical Test File: RESTR-FIND-001 - Restriction Site Detection
// Methods Under Test: RestrictionAnalyzer.FindSites, FindAllSites, GetEnzyme
//
// Evidence Sources:
// - Wikipedia: Restriction enzyme (https://en.wikipedia.org/wiki/Restriction_enzyme)
// - Wikipedia: Restriction site (https://en.wikipedia.org/wiki/Restriction_site)
// - Wikipedia: EcoRI (https://en.wikipedia.org/wiki/EcoRI)
// - Roberts RJ (1976) Restriction endonucleases
// - REBASE: The Restriction Enzyme Database
//
// Algorithm Documentation: docs/algorithms/MolTools/Restriction_Site_Detection.md

using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical test class for RESTR-FIND-001: Restriction Site Detection.
/// Tests FindSites, FindAllSites, and GetEnzyme methods of RestrictionAnalyzer.
/// </summary>
[TestFixture]
public class RestrictionAnalyzer_FindSites_Tests
{
    #region Enzyme Lookup Tests (GetEnzyme)

    /// <summary>
    /// Evidence: Wikipedia EcoRI - "GAATTC" recognition sequence, cuts between G and A.
    /// Invariant: Enzyme lookup returns correct enzyme with expected properties.
    /// </summary>
    [Test]
    public void GetEnzyme_EcoRI_ReturnsEnzymeWithCorrectProperties()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRI");

        Assert.Multiple(() =>
        {
            Assert.That(enzyme, Is.Not.Null, "EcoRI should be in enzyme database");
            Assert.That(enzyme!.Name, Is.EqualTo("EcoRI"), "Name should match");
            Assert.That(enzyme.RecognitionSequence, Is.EqualTo("GAATTC"), "Recognition sequence from Wikipedia");
            Assert.That(enzyme.CutPositionForward, Is.EqualTo(1), "Cuts after G (position 1)");
            Assert.That(enzyme.CutPositionReverse, Is.EqualTo(5), "Reverse cut at position 5");
            Assert.That(enzyme.RecognitionLength, Is.EqualTo(6), "6-cutter enzyme");
            Assert.That(enzyme.OverhangType, Is.EqualTo(OverhangType.FivePrime), "Produces 5' overhang (AATT)");
            Assert.That(enzyme.IsBluntEnd, Is.False, "EcoRI is not a blunt cutter");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - BamHI recognition sequence "GGATCC".
    /// </summary>
    [Test]
    public void GetEnzyme_BamHI_ReturnsEnzymeWithCorrectProperties()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("BamHI");

        Assert.Multiple(() =>
        {
            Assert.That(enzyme, Is.Not.Null, "BamHI should be in enzyme database");
            Assert.That(enzyme!.Name, Is.EqualTo("BamHI"), "Name should match");
            Assert.That(enzyme.RecognitionSequence, Is.EqualTo("GGATCC"), "Recognition sequence from Wikipedia");
            Assert.That(enzyme.CutPositionForward, Is.EqualTo(1), "Cuts after G (position 1)");
            Assert.That(enzyme.CutPositionReverse, Is.EqualTo(5), "Reverse cut at position 5");
            Assert.That(enzyme.RecognitionLength, Is.EqualTo(6), "6-cutter enzyme");
            Assert.That(enzyme.OverhangType, Is.EqualTo(OverhangType.FivePrime), "Produces 5' overhang (GATC)");
            Assert.That(enzyme.IsBluntEnd, Is.False, "BamHI is not a blunt cutter");
        });
    }

    /// <summary>
    /// Invariant #4: Case Insensitivity - GetEnzyme("EcoRI") == GetEnzyme("ecori").
    /// </summary>
    [Test]
    [TestCase("EcoRI")]
    [TestCase("ecori")]
    [TestCase("ECORI")]
    [TestCase("eCoRi")]
    public void GetEnzyme_CaseInsensitive_ReturnsCorrectEnzyme(string enzymeName)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(enzymeName);

        Assert.Multiple(() =>
        {
            Assert.That(enzyme, Is.Not.Null, $"Enzyme should be found for '{enzymeName}'");
            Assert.That(enzyme!.Name, Is.EqualTo("EcoRI"), "Canonical name should be returned");
            Assert.That(enzyme.RecognitionSequence, Is.EqualTo("GAATTC"));
        });
    }

    /// <summary>
    /// Invariant: Unknown enzyme name returns null (does not throw).
    /// </summary>
    [Test]
    [TestCase("UnknownEnzyme")]
    [TestCase("FakeI")]
    [TestCase("")]
    [TestCase("   ")]
    public void GetEnzyme_UnknownOrInvalidName_ReturnsNull(string enzymeName)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(enzymeName);
        Assert.That(enzyme, Is.Null, $"Unknown enzyme '{enzymeName}' should return null");
    }

    #endregion

    #region Enzyme Database Coverage Tests

    /// <summary>
    /// Evidence: REBASE - Comprehensive enzyme database with 30+ common enzymes.
    /// </summary>
    [Test]
    public void Enzymes_ContainsAtLeast30CommonEnzymes()
    {
        var enzymes = RestrictionAnalyzer.Enzymes;

        Assert.Multiple(() =>
        {
            Assert.That(enzymes.Count, Is.GreaterThanOrEqualTo(30), "Database should have 30+ enzymes");

            // Common 6-cutters (Wikipedia)
            Assert.That(enzymes.ContainsKey("EcoRI"), "EcoRI should be present");
            Assert.That(enzymes.ContainsKey("BamHI"), "BamHI should be present");
            Assert.That(enzymes.ContainsKey("HindIII"), "HindIII should be present");
            Assert.That(enzymes.ContainsKey("PstI"), "PstI should be present");
            Assert.That(enzymes.ContainsKey("SalI"), "SalI should be present");
            Assert.That(enzymes.ContainsKey("XhoI"), "XhoI should be present");

            // Rare cutters / 8-cutters
            Assert.That(enzymes.ContainsKey("NotI"), "NotI (8-cutter) should be present");
            Assert.That(enzymes.ContainsKey("SfiI"), "SfiI (8-cutter) should be present");

            // 4-cutters (frequent cutters)
            Assert.That(enzymes.ContainsKey("MspI") || enzymes.ContainsKey("HpaII"),
                "At least one 4-cutter should be present");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - Recognition sequences range from 4-8 bp.
    /// </summary>
    [Test]
    [TestCase(4, Description = "4-cutters (frequent cutters)")]
    [TestCase(6, Description = "6-cutters (most common)")]
    [TestCase(8, Description = "8-cutters (rare cutters)")]
    public void GetEnzymesByCutLength_ReturnsEnzymesOfCorrectLength(int cutLength)
    {
        var enzymes = RestrictionAnalyzer.GetEnzymesByCutLength(cutLength).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(enzymes, Is.Not.Empty, $"Should have enzymes with {cutLength}-bp recognition");
            Assert.That(enzymes.All(e => e.RecognitionLength == cutLength),
                $"All enzymes should have recognition length of {cutLength}");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - Enzymes produce blunt ends, 5' overhangs, or 3' overhangs.
    /// </summary>
    [Test]
    public void GetBluntCutters_ReturnsOnlyBluntEndEnzymes()
    {
        var bluntCutters = RestrictionAnalyzer.GetBluntCutters().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(bluntCutters, Is.Not.Empty, "Should have blunt cutters");
            Assert.That(bluntCutters.All(e => e.IsBluntEnd), "All should be blunt end");
            Assert.That(bluntCutters.All(e => e.OverhangType == OverhangType.Blunt),
                "All should have Blunt overhang type");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - Sticky end enzymes produce overhangs.
    /// </summary>
    [Test]
    public void GetStickyCutters_ReturnsOnlyStickyEndEnzymes()
    {
        var stickyCutters = RestrictionAnalyzer.GetStickyCutters().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(stickyCutters, Is.Not.Empty, "Should have sticky cutters");
            Assert.That(stickyCutters.All(e => !e.IsBluntEnd), "None should be blunt end");
            Assert.That(stickyCutters.All(e =>
                e.OverhangType == OverhangType.FivePrime || e.OverhangType == OverhangType.ThreePrime),
                "All should have 5' or 3' overhang type");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - Enzymes classified by overhang type.
    /// </summary>
    [Test]
    public void EnzymeDatabase_ContainsAllOverhangTypes()
    {
        var enzymes = RestrictionAnalyzer.Enzymes.Values.ToList();

        Assert.Multiple(() =>
        {
            Assert.That(enzymes.Any(e => e.OverhangType == OverhangType.FivePrime),
                "Should have 5' overhang enzymes (e.g., EcoRI)");
            Assert.That(enzymes.Any(e => e.OverhangType == OverhangType.ThreePrime),
                "Should have 3' overhang enzymes (e.g., PstI)");
            Assert.That(enzymes.Any(e => e.OverhangType == OverhangType.Blunt),
                "Should have blunt enzymes (e.g., EcoRV)");
        });
    }

    #endregion

    #region FindSites - Basic Functionality Tests

    /// <summary>
    /// Evidence: Wikipedia EcoRI - Recognizes GAATTC sequence.
    /// Invariant #1: Position Range - 0 ≤ Position ≤ sequence.Length - recognitionLength.
    /// </summary>
    [Test]
    public void FindSites_EcoRI_FindsSiteAtCorrectPosition()
    {
        // Arrange: GAATTC at position 3
        var sequence = new DnaSequence("AAAGAATTCAAA");

        // Act
        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "Should find exactly one EcoRI site");
            Assert.That(sites[0].Position, Is.EqualTo(3), "Site should be at position 3");
            Assert.That(sites[0].Position, Is.GreaterThanOrEqualTo(0), "Invariant: position >= 0");
            Assert.That(sites[0].Position, Is.LessThanOrEqualTo(sequence.Length - 6),
                "Invariant: position <= length - recognition length");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - BamHI recognizes GGATCC.
    /// </summary>
    [Test]
    public void FindSites_BamHI_FindsSiteAtCorrectPosition()
    {
        // Arrange: GGATCC at position 4
        var sequence = new DnaSequence("AAAAGGATCCAAAA");

        // Act
        var sites = RestrictionAnalyzer.FindSites(sequence, "BamHI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "Should find exactly one BamHI site");
            Assert.That(sites[0].Position, Is.EqualTo(4), "Site should be at position 4");
            Assert.That(sites[0].Enzyme.Name, Is.EqualTo("BamHI"));
        });
    }

    /// <summary>
    /// Invariant #2: Cut Position - Position ≤ CutPosition ≤ Position + recognitionLength.
    /// Evidence: Wikipedia EcoRI - Cuts after G (position 1 in recognition sequence).
    /// </summary>
    [Test]
    public void FindSites_EcoRI_CutPositionCalculatedCorrectly()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            var site = sites[0];
            var enzyme = site.Enzyme;

            // Cut position = Position + CutPositionForward = 3 + 1 = 4
            Assert.That(site.CutPosition, Is.EqualTo(4), "Cut position should be 4 (3 + 1)");

            // Invariant check
            Assert.That(site.CutPosition, Is.GreaterThanOrEqualTo(site.Position),
                "Invariant: CutPosition >= Position");
            Assert.That(site.CutPosition, Is.LessThanOrEqualTo(site.Position + enzyme.RecognitionLength),
                "Invariant: CutPosition <= Position + RecognitionLength");
        });
    }

    /// <summary>
    /// Invariant #3: Recognition Sequence Length matches enzyme.
    /// </summary>
    [Test]
    public void FindSites_EcoRI_RecognizedSequenceMatchesEnzyme()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            var site = sites[0];
            Assert.That(site.RecognizedSequence, Is.EqualTo("GAATTC"),
                "Recognized sequence should match the actual bases");
            Assert.That(site.RecognizedSequence.Length, Is.EqualTo(site.Enzyme.RecognitionLength),
                "Invariant: RecognizedSequence.Length == Enzyme.RecognitionLength");
        });
    }

    /// <summary>
    /// Test: Multiple sites in sequence are all found.
    /// </summary>
    [Test]
    public void FindSites_MultipleSites_FindsAllAtCorrectPositions()
    {
        // Two EcoRI sites at positions 0 and 9
        var sequence = new DnaSequence("GAATTCAAAGAATTC");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .OrderBy(s => s.Position)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(2), "Should find both EcoRI sites");
            Assert.That(sites[0].Position, Is.EqualTo(0), "First site at position 0");
            Assert.That(sites[1].Position, Is.EqualTo(9), "Second site at position 9");
        });
    }

    /// <summary>
    /// Test: No false positives when recognition sequence is absent.
    /// </summary>
    [Test]
    public void FindSites_NoSites_ReturnsEmptyCollection()
    {
        var sequence = new DnaSequence("AAAAAAAAAA");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI").ToList();

        Assert.That(sites, Is.Empty, "Should find no EcoRI sites in poly-A sequence");
    }

    /// <summary>
    /// Test: String overload works identically to DnaSequence overload.
    /// </summary>
    [Test]
    public void FindSites_StringOverload_ProducesIdenticalResults()
    {
        var sequenceString = "AAAGAATTCAAA";
        var sequenceDna = new DnaSequence(sequenceString);

        var sitesFromString = RestrictionAnalyzer.FindSites(sequenceString, "EcoRI").ToList();
        var sitesFromDna = RestrictionAnalyzer.FindSites(sequenceDna, "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sitesFromString.Count, Is.EqualTo(sitesFromDna.Count),
                "Both overloads should find same number of sites");

            for (int i = 0; i < sitesFromString.Count; i++)
            {
                Assert.That(sitesFromString[i].Position, Is.EqualTo(sitesFromDna[i].Position),
                    $"Site {i} position should match");
                Assert.That(sitesFromString[i].IsForwardStrand, Is.EqualTo(sitesFromDna[i].IsForwardStrand),
                    $"Site {i} strand should match");
            }
        });
    }

    #endregion

    #region FindSites - Multiple Enzymes Tests

    /// <summary>
    /// Test: FindSites with multiple enzyme names finds sites for all.
    /// </summary>
    [Test]
    public void FindSites_MultipleEnzymes_FindsSitesForAllEnzymes()
    {
        // EcoRI at 0, BamHI at 9
        var sequence = new DnaSequence("GAATTCAAAGGATCC");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI", "BamHI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(2), "Should find one site for each enzyme");
            Assert.That(sites.Select(s => s.Enzyme.Name).Distinct().Count(), Is.EqualTo(2),
                "Should have sites from 2 different enzymes");
            Assert.That(sites.Any(s => s.Enzyme.Name == "EcoRI"), "Should find EcoRI site");
            Assert.That(sites.Any(s => s.Enzyme.Name == "BamHI"), "Should find BamHI site");
        });
    }

    /// <summary>
    /// Test: FindAllSites searches with all enzymes in database.
    /// </summary>
    [Test]
    public void FindAllSites_FindsSitesFromMultipleEnzymes()
    {
        // Contains recognition sites for EcoRI (GAATTC), BamHI (GGATCC), HindIII (AAGCTT)
        var sequence = new DnaSequence("GAATTCGGATCCAAGCTT");

        var sites = RestrictionAnalyzer.FindAllSites(sequence)
            .Where(s => s.IsForwardStrand)
            .ToList();

        var enzymeNames = sites.Select(s => s.Enzyme.Name).Distinct().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.GreaterThanOrEqualTo(3), "Should find at least 3 sites");
            Assert.That(enzymeNames, Has.Count.GreaterThanOrEqualTo(3),
                "Should find sites from at least 3 different enzymes");
            Assert.That(enzymeNames, Does.Contain("EcoRI"));
            Assert.That(enzymeNames, Does.Contain("BamHI"));
            Assert.That(enzymeNames, Does.Contain("HindIII"));
        });
    }

    #endregion

    #region FindSites - Custom Enzyme Tests

    /// <summary>
    /// Test: Custom enzyme can be used for site finding.
    /// </summary>
    [Test]
    public void FindSites_CustomEnzyme_FindsSiteCorrectly()
    {
        // ATAT at position 2 in AAATATAAA
        var customEnzyme = new RestrictionEnzyme("CustomI", "ATAT", 2, 2, "Custom enzyme");
        var sequence = new DnaSequence("AAATATAAA");

        var sites = RestrictionAnalyzer.FindSites(sequence, customEnzyme)
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "Should find custom enzyme site");
            Assert.That(sites[0].Position, Is.EqualTo(2), "ATAT starts at position 2");
            Assert.That(sites[0].Enzyme.Name, Is.EqualTo("CustomI"), "Enzyme name should match");
            Assert.That(sites[0].RecognizedSequence, Is.EqualTo("ATAT"));
            Assert.That(sites[0].CutPosition, Is.EqualTo(4), "Cut at 2 + 2 = 4");
        });
    }

    /// <summary>
    /// Test: Custom enzyme with odd cut positions.
    /// </summary>
    [Test]
    public void FindSites_CustomEnzymeWithAsymmetricCut_CutPositionCorrect()
    {
        // Asymmetric cut: cuts at position 1 forward, position 3 reverse
        var customEnzyme = new RestrictionEnzyme("AsymI", "ACTG", 1, 3, "Asymmetric cutter");
        var sequence = new DnaSequence("AAACTGAAA");

        var sites = RestrictionAnalyzer.FindSites(sequence, customEnzyme)
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Position, Is.EqualTo(2), "Site at position 2");
            Assert.That(sites[0].CutPosition, Is.EqualTo(3), "Cut at 2 + 1 = 3");
        });
    }

    #endregion

    #region Edge Cases and Error Handling

    /// <summary>
    /// Invariant #5: Empty sequence returns no sites.
    /// </summary>
    [Test]
    public void FindSites_EmptyStringSequence_ReturnsEmptyCollection()
    {
        var sites = RestrictionAnalyzer.FindSites("", "EcoRI").ToList();
        Assert.That(sites, Is.Empty, "Empty sequence should have no sites");
    }

    /// <summary>
    /// Test: Null DnaSequence throws ArgumentNullException.
    /// </summary>
    [Test]
    public void FindSites_NullDnaSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.FindSites((DnaSequence)null!, "EcoRI").ToList());
    }

    /// <summary>
    /// Invariant #6: Unknown enzyme name throws ArgumentException.
    /// </summary>
    [Test]
    public void FindSites_UnknownEnzyme_ThrowsArgumentException()
    {
        var sequence = new DnaSequence("ACGT");

        var ex = Assert.Throws<ArgumentException>(() =>
            RestrictionAnalyzer.FindSites(sequence, "UnknownEnzyme").ToList());

        Assert.That(ex!.Message, Does.Contain("UnknownEnzyme").IgnoreCase,
            "Exception message should mention the unknown enzyme name");
    }

    /// <summary>
    /// Test: Empty enzyme name throws ArgumentNullException.
    /// </summary>
    [Test]
    [TestCase("")]
    [TestCase(null)]
    public void FindSites_EmptyOrNullEnzymeName_ThrowsArgumentNullException(string? enzymeName)
    {
        var sequence = new DnaSequence("ACGT");

        Assert.Throws<ArgumentNullException>(() =>
            RestrictionAnalyzer.FindSites(sequence, enzymeName!).ToList());
    }

    /// <summary>
    /// Test: Sequence shorter than recognition sequence returns no sites.
    /// </summary>
    [Test]
    public void FindSites_SequenceShorterThanRecognition_ReturnsEmpty()
    {
        // EcoRI recognition is 6 bp, sequence is 4 bp
        var sequence = new DnaSequence("ACGT");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI").ToList();

        Assert.That(sites, Is.Empty, "Short sequence cannot contain recognition site");
    }

    /// <summary>
    /// Test: Recognition site at very end of sequence is found.
    /// </summary>
    [Test]
    public void FindSites_SiteAtEndOfSequence_IsFound()
    {
        // EcoRI site at positions 6-11 (exactly at end)
        var sequence = new DnaSequence("AAAAAAGAATTC");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "Should find site at end");
            Assert.That(sites[0].Position, Is.EqualTo(6), "Site at position 6");
            Assert.That(sites[0].Position + sites[0].Enzyme.RecognitionLength, Is.EqualTo(sequence.Length),
                "Site should end exactly at sequence end");
        });
    }

    /// <summary>
    /// Test: Recognition site at very beginning of sequence is found.
    /// </summary>
    [Test]
    public void FindSites_SiteAtStartOfSequence_IsFound()
    {
        var sequence = new DnaSequence("GAATTCAAAAAA");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "Should find site at start");
            Assert.That(sites[0].Position, Is.EqualTo(0), "Site at position 0");
        });
    }

    #endregion

    #region Enzyme Overhang Type Tests

    /// <summary>
    /// Evidence: Wikipedia - PstI produces 3' overhang.
    /// </summary>
    [Test]
    public void RestrictionEnzyme_PstI_HasThreePrimeOverhang()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("PstI")!;

        Assert.Multiple(() =>
        {
            Assert.That(enzyme.IsBluntEnd, Is.False, "PstI is not blunt");
            Assert.That(enzyme.OverhangType, Is.EqualTo(OverhangType.ThreePrime),
                "PstI produces 3' overhang");
            // 3' overhang means CutPositionForward > CutPositionReverse
            Assert.That(enzyme.CutPositionForward, Is.GreaterThan(enzyme.CutPositionReverse),
                "3' overhang: forward cut after reverse cut");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - EcoRV is a blunt cutter.
    /// </summary>
    [Test]
    public void RestrictionEnzyme_EcoRV_IsBluntCutter()
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme("EcoRV")!;

        Assert.Multiple(() =>
        {
            Assert.That(enzyme.IsBluntEnd, Is.True, "EcoRV is blunt");
            Assert.That(enzyme.OverhangType, Is.EqualTo(OverhangType.Blunt), "EcoRV produces blunt end");
            // Blunt means CutPositionForward == CutPositionReverse
            Assert.That(enzyme.CutPositionForward, Is.EqualTo(enzyme.CutPositionReverse),
                "Blunt end: forward and reverse cut at same position");
        });
    }

    #endregion

    #region IUPAC Degenerate Recognition Tests

    /// <summary>
    /// Evidence: HincII recognizes GTYRAC (Y = C/T, R = A/G) — a degenerate sequence.
    /// All 4 valid IUPAC combinations must be recognized.
    /// Source: REBASE, Wikipedia: Restriction enzyme
    /// </summary>
    [Test]
    [TestCase("GTCAAC", Description = "Y=C, R=A")]
    [TestCase("GTTAAC", Description = "Y=T, R=A")]
    [TestCase("GTCGAC", Description = "Y=C, R=G")]
    [TestCase("GTTGAC", Description = "Y=T, R=G")]
    public void FindSites_HincII_MatchesAllDegenerateCombinations(string recognition)
    {
        var sequence = new DnaSequence($"AA{recognition}AA");

        var sites = RestrictionAnalyzer.FindSites(sequence, "HincII")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                $"HincII (GTYRAC) should match {recognition}");
            Assert.That(sites[0].Position, Is.EqualTo(2),
                "Site at position 2");
            Assert.That(sites[0].RecognizedSequence, Is.EqualTo(recognition),
                "RecognizedSequence should be the actual bases, not the IUPAC pattern");
        });
    }

    /// <summary>
    /// Evidence: HincII GTYRAC — Y requires C or T only.
    /// A base that does not satisfy the IUPAC code must NOT match.
    /// </summary>
    [Test]
    [TestCase("GTAAAC", Description = "A does not satisfy Y (C/T)")]
    [TestCase("GTGAAC", Description = "G does not satisfy Y (C/T)")]
    [TestCase("GTCTAC", Description = "T does not satisfy R (A/G)")]
    [TestCase("GTCCAC", Description = "C does not satisfy R (A/G)")]
    public void FindSites_HincII_RejectsNonMatchingDegenerateBases(string nonMatch)
    {
        var sequence = new DnaSequence($"AA{nonMatch}AA");

        var sites = RestrictionAnalyzer.FindSites(sequence, "HincII")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.That(sites, Is.Empty,
            $"HincII (GTYRAC) should NOT match {nonMatch}");
    }

    /// <summary>
    /// Evidence: HincII is a blunt cutter (3, 3) even with degenerate recognition.
    /// Verifies that IUPAC matching does not affect cut position calculation.
    /// </summary>
    [Test]
    public void FindSites_HincII_DegenerateEnzymeProducesCorrectCutPosition()
    {
        // GTCAAC at position 3
        var sequence = new DnaSequence("AAAGTCAACAAA");

        var site = RestrictionAnalyzer.FindSites(sequence, "HincII")
            .First(s => s.IsForwardStrand);

        Assert.Multiple(() =>
        {
            Assert.That(site.CutPosition, Is.EqualTo(6), "Cut at 3 + 3 = 6");
            Assert.That(site.Enzyme.IsBluntEnd, Is.True, "HincII is blunt cutter");
            Assert.That(site.Enzyme.OverhangType, Is.EqualTo(OverhangType.Blunt));
        });
    }

    #endregion

    #region External Source Verification Tests

    /// <summary>
    /// Evidence: Wikipedia EcoRI — "5' end overhangs of AATT", "G↓AATTC / CTTAA↓G".
    /// Verifies exact overhang bases match Wikipedia data.
    /// Source: https://en.wikipedia.org/wiki/EcoRI
    /// </summary>
    [Test]
    public void FindSites_EcoRI_OverhangSequenceIsAATT()
    {
        var sequence = new DnaSequence("AAAGAATTCAAA");

        var forwardSite = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .First(s => s.IsForwardStrand);
        var reverseSite = RestrictionAnalyzer.FindSites(sequence, "EcoRI")
            .First(s => !s.IsForwardStrand);

        int fwdCut = forwardSite.CutPosition;
        int revCut = reverseSite.CutPosition;
        string overhang = sequence.Sequence.Substring(fwdCut, revCut - fwdCut);

        Assert.Multiple(() =>
        {
            Assert.That(fwdCut, Is.LessThan(revCut),
                "5' overhang: forward cut before reverse cut");
            Assert.That(overhang, Is.EqualTo("AATT"),
                "EcoRI overhang should be AATT (Wikipedia: EcoRI)");
            Assert.That(overhang.Length, Is.EqualTo(4),
                "EcoRI produces 4-base overhang");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia BamHI — "cleaves just after the 5'-guanine on each strand"
    /// producing 4 bp 5' sticky ends with GATC overhang.
    /// Source: https://en.wikipedia.org/wiki/BamHI
    /// </summary>
    [Test]
    public void FindSites_BamHI_OverhangSequenceIsGATC()
    {
        var sequence = new DnaSequence("AAAAGGATCCAAAA");

        var forwardSite = RestrictionAnalyzer.FindSites(sequence, "BamHI")
            .First(s => s.IsForwardStrand);
        var reverseSite = RestrictionAnalyzer.FindSites(sequence, "BamHI")
            .First(s => !s.IsForwardStrand);

        int fwdCut = forwardSite.CutPosition;
        int revCut = reverseSite.CutPosition;
        string overhang = sequence.Sequence.Substring(fwdCut, revCut - fwdCut);

        Assert.Multiple(() =>
        {
            Assert.That(fwdCut, Is.LessThan(revCut),
                "5' overhang: forward cut before reverse cut");
            Assert.That(overhang, Is.EqualTo("GATC"),
                "BamHI overhang should be GATC (Wikipedia: BamHI)");
            Assert.That(overhang.Length, Is.EqualTo(4),
                "BamHI produces 4-base overhang");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia HindIII — "cleavage of this sequence between the AA's
    /// results in 5' overhangs". Cut pattern: A↓AGCTT / TTCGA↓A.
    /// Source: https://en.wikipedia.org/wiki/HindIII
    /// </summary>
    [Test]
    public void FindSites_HindIII_OverhangSequenceIsAGCT()
    {
        // AAGCTT at position 3
        var sequence = new DnaSequence("AAAAAGCTTAAA");

        var forwardSite = RestrictionAnalyzer.FindSites(sequence, "HindIII")
            .First(s => s.IsForwardStrand);
        var reverseSite = RestrictionAnalyzer.FindSites(sequence, "HindIII")
            .First(s => !s.IsForwardStrand);

        int fwdCut = forwardSite.CutPosition;
        int revCut = reverseSite.CutPosition;
        string overhang = sequence.Sequence.Substring(fwdCut, revCut - fwdCut);

        Assert.Multiple(() =>
        {
            Assert.That(fwdCut, Is.LessThan(revCut),
                "5' overhang: forward cut before reverse cut");
            Assert.That(overhang, Is.EqualTo("AGCT"),
                "HindIII overhang should be AGCT (Wikipedia: HindIII)");
            Assert.That(overhang.Length, Is.EqualTo(4),
                "HindIII produces 4-base overhang");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia PstI — "generates fragments with 3'-cohesive termini",
    /// cut pattern: CTGCA↓G / G↓ACGTC, producing 4 bp 3' overhang.
    /// Source: https://en.wikipedia.org/wiki/PstI
    /// </summary>
    [Test]
    public void FindSites_PstI_ProducesThreePrimeOverhangOfFourBases()
    {
        // CTGCAG at position 2
        var sequence = new DnaSequence("AACTGCAGAA");

        var forwardSite = RestrictionAnalyzer.FindSites(sequence, "PstI")
            .First(s => s.IsForwardStrand);
        var reverseSite = RestrictionAnalyzer.FindSites(sequence, "PstI")
            .First(s => !s.IsForwardStrand);

        Assert.Multiple(() =>
        {
            Assert.That(forwardSite.CutPosition, Is.GreaterThan(reverseSite.CutPosition),
                "3' overhang: forward cut after reverse cut (Wikipedia: PstI)");
            Assert.That(forwardSite.CutPosition - reverseSite.CutPosition, Is.EqualTo(4),
                "PstI produces 4-base 3' overhang (Wikipedia: PstI)");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia EcoRV — "creates blunt ends", "GAT|ATC".
    /// Forward and reverse cuts at identical position.
    /// Source: https://en.wikipedia.org/wiki/EcoRV
    /// </summary>
    [Test]
    public void FindSites_EcoRV_ProducesBluntEndCuts()
    {
        // GATATC at position 3
        var sequence = new DnaSequence("AAAGATATCAAA");

        var forwardSite = RestrictionAnalyzer.FindSites(sequence, "EcoRV")
            .First(s => s.IsForwardStrand);
        var reverseSite = RestrictionAnalyzer.FindSites(sequence, "EcoRV")
            .First(s => !s.IsForwardStrand);

        Assert.That(forwardSite.CutPosition, Is.EqualTo(reverseSite.CutPosition),
            "Blunt end: forward and reverse cut at same position (Wikipedia: EcoRV)");
    }

    /// <summary>
    /// Evidence: Wikipedia Restriction enzyme — "Many of them are palindromic, meaning
    /// the base sequence reads the same backwards and forwards."
    /// Palindromic sites are recognized on both strands at the same position.
    /// Source: https://en.wikipedia.org/wiki/Restriction_enzyme#Recognition_site
    /// </summary>
    [Test]
    public void FindSites_PalindromicSite_FoundOnBothStrandsAtSamePosition()
    {
        // EcoRI recognizes palindromic GAATTC (reverse complement is also GAATTC)
        var sequence = new DnaSequence("AAAGAATTCAAA");

        var sites = RestrictionAnalyzer.FindSites(sequence, "EcoRI").ToList();
        var forwardSites = sites.Where(s => s.IsForwardStrand).ToList();
        var reverseSites = sites.Where(s => !s.IsForwardStrand).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(forwardSites, Has.Count.EqualTo(1),
                "Should find palindromic site on forward strand");
            Assert.That(reverseSites, Has.Count.EqualTo(1),
                "Should find palindromic site on reverse strand");
            Assert.That(forwardSites[0].Position, Is.EqualTo(reverseSites[0].Position),
                "Palindromic site at same forward-strand position on both strands");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia Restriction enzyme examples table — comprehensive verification
    /// of cut positions for all enzymes listed in the Wikipedia table.
    /// Source: https://en.wikipedia.org/wiki/Restriction_enzyme#Examples
    /// </summary>
    [Test]
    [TestCase("EcoRI", "GAATTC", 1, 5, Description = "Wikipedia: G↓AATTC / CTTAA↓G")]
    [TestCase("BamHI", "GGATCC", 1, 5, Description = "Wikipedia: G↓GATCC / CCTAG↓G")]
    [TestCase("HindIII", "AAGCTT", 1, 5, Description = "Wikipedia: A↓AGCTT / TTCGA↓A")]
    [TestCase("TaqI", "TCGA", 1, 3, Description = "Wikipedia: T↓CGA / AGC↓T")]
    [TestCase("NotI", "GCGGCCGC", 2, 6, Description = "Wikipedia: GC↓GGCCGC / CGCCGG↓CG")]
    [TestCase("Sau3AI", "GATC", 0, 4, Description = "Wikipedia: ↓GATC / CTAG↓")]
    [TestCase("AluI", "AGCT", 2, 2, Description = "Wikipedia: AG↓CT / TC↓GA (blunt)")]
    [TestCase("HaeIII", "GGCC", 2, 2, Description = "Wikipedia: GG↓CC / CC↓GG (blunt)")]
    [TestCase("EcoRV", "GATATC", 3, 3, Description = "Wikipedia: GAT↓ATC / CTA↓TAG (blunt)")]
    [TestCase("SmaI", "CCCGGG", 3, 3, Description = "Wikipedia: CCC↓GGG / GGG↓CCC (blunt)")]
    [TestCase("PstI", "CTGCAG", 5, 1, Description = "Wikipedia: CTGCA↓G / G↓ACGTC")]
    [TestCase("KpnI", "GGTACC", 5, 1, Description = "Wikipedia: GGTAC↓C / C↓CATGG")]
    [TestCase("SacI", "GAGCTC", 5, 1, Description = "Wikipedia: GAGCT↓C / C↓TCGAG")]
    [TestCase("SalI", "GTCGAC", 1, 5, Description = "Wikipedia: G↓TCGAC / CAGCT↓G")]
    [TestCase("ScaI", "AGTACT", 3, 3, Description = "Wikipedia: AGT↓ACT / TCA↓TGA (blunt)")]
    [TestCase("SpeI", "ACTAGT", 1, 5, Description = "Wikipedia: A↓CTAGT / TGATC↓A")]
    [TestCase("SphI", "GCATGC", 5, 1, Description = "Wikipedia: GCATG↓C / C↓GTACG")]
    [TestCase("StuI", "AGGCCT", 3, 3, Description = "Wikipedia: AGG↓CCT / TCC↓GGA (blunt)")]
    [TestCase("XbaI", "TCTAGA", 1, 5, Description = "Wikipedia: T↓CTAGA / AGATC↓T")]
    [TestCase("BglII", "AGATCT", 1, 5, Description = "Wikipedia: A↓GATCT / TCTAG↓A")]
    [TestCase("ApaI", "GGGCCC", 5, 1, Description = "Wikipedia: GGGCC↓C / C↓CCGGG")]
    public void EnzymeDatabase_CutPositions_MatchWikipedia(
        string name, string expectedRecognition, int expectedCutFwd, int expectedCutRev)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(name);

        Assert.Multiple(() =>
        {
            Assert.That(enzyme, Is.Not.Null, $"{name} should be in database");
            Assert.That(enzyme!.RecognitionSequence, Is.EqualTo(expectedRecognition),
                $"{name} recognition sequence (Wikipedia)");
            Assert.That(enzyme.CutPositionForward, Is.EqualTo(expectedCutFwd),
                $"{name} forward cut position (Wikipedia)");
            Assert.That(enzyme.CutPositionReverse, Is.EqualTo(expectedCutRev),
                $"{name} reverse cut position (Wikipedia)");
        });
    }

    #endregion
}
