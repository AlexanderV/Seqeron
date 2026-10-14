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

namespace Seqeron.Genomics.Tests.Unit.MolTools;

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

    #region Reference-Locked Site Search (RESTR-FIND-001 review, 2026-10)

    // ─────────────────────────────────────────────────────────────────────────────────────
    // Oracle for this region: Biopython 1.88 `Bio.Restriction`, whose enzyme table is the
    // REBASE EMBOSS files **version 404 (2024)** (stated in the generated
    // Bio/Restriction/Restriction_Dictionary.py header). Every expected number below was
    // produced by running that package this session; `search()` returns 1-BASED top-strand
    // cut positions, so the relation to this library is
    //     biopython value == forward-strand RestrictionSite.CutPosition + 1.
    // Relevant reference code: Restriction.py FormattedSeq.finditer L215-228 (circular search
    // over `data + data[1:size]`), AbstractCut._drop L542-568 (linear cut-feasibility filter
    // `1 < cut_on_watson <= length` and the same for `cut_on_crick = cut_on_watson - ovhg`;
    // circular branch folds out-of-range cuts by ±length), Defined._modify L789-814 /
    // Defined._rev_modify L817-825, Palindromic._search L967-978, NonPalindromic._search
    // L996-1019.
    // ─────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// A cut can only fall on a phosphodiester bond that exists: Sau3AI (^GATC) cannot cut a
    /// GATC sitting at the very 5' end of a linear molecule, and cannot cut one at the very
    /// 3' end either (the bottom-strand cut would fall past the last base).
    /// Reference: Bio.Restriction.Sau3AI.search(Seq("GATCGATCGATC")) == [5] — one cut only,
    /// 1-based, i.e. CutPosition 4 here. (AbstractCut._drop, Restriction.py L542-568.)
    /// </summary>
    [Test]
    public void FindSites_Sau3AI_EdgeSites_AreNotCuttable_MatchesBiopython()
    {
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence("GATCGATCGATC"), "Sau3AI")
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1),
                "GATC occurs at 0, 4 and 8 but only the middle one has both cuts inside the molecule");
            Assert.That(sites[0].Position, Is.EqualTo(4), "the cuttable GATC starts at 4");
            Assert.That(sites[0].CutPosition, Is.EqualTo(4),
                "Biopython Sau3AI.search(Seq(\"GATCGATCGATC\")) == [5] (1-based) == 4 here");
        });
    }

    /// <summary>
    /// Reference: Bio.Restriction.Sau3AI.search(Seq("GATC")) == [] and
    /// Bio.Restriction.MboI.search(Seq("GATC")) == [] — a ^GATC enzyme cannot cut a molecule
    /// that is nothing but its own recognition site.
    /// </summary>
    [Test]
    [TestCase("Sau3AI")]
    [TestCase("MboI")]
    public void FindSites_GatcEdgeCutter_WholeSequenceIsTheSite_ReportsNothing(string enzymeName)
    {
        RestrictionAnalyzer.FindSites(new DnaSequence("GATC"), enzymeName)
            .Should().BeEmpty("neither of this enzyme's two cuts lands on an existing bond");
    }

    /// <summary>
    /// Reference: Bio.Restriction.Sau3AI.search(Seq("AAAAAAAAGATC")) == [] (linear) — the
    /// bottom-strand cut GATC^ would fall past the last base.
    /// </summary>
    [Test]
    public void FindSites_Sau3AI_SiteAtThreePrimeEnd_ReportsNothing()
    {
        RestrictionAnalyzer.FindSites(new DnaSequence("AAAAAAAAGATC"), "Sau3AI")
            .Should().BeEmpty("the GATC^ bottom-strand cut would fall outside the molecule");
    }

    /// <summary>
    /// The feasibility filter must not swallow ordinary sites: an enzyme that cuts strictly
    /// inside its own window is still reported when its site fills the whole molecule.
    /// Reference: Bio.Restriction.EcoRI.search(Seq("GAATTC")) == [2] (1-based) == 1 here.
    /// </summary>
    [Test]
    public void FindSites_EcoRI_SiteFillsWholeSequence_IsStillReported()
    {
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence("GAATTC"), "EcoRI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Count(s => s.IsForwardStrand), Is.EqualTo(1));
            Assert.That(sites.Single(s => s.IsForwardStrand).CutPosition, Is.EqualTo(1),
                "G^AATTC — Biopython reports 2 (1-based)");
            Assert.That(sites.Single(s => !s.IsForwardStrand).CutPosition, Is.EqualTo(5),
                "GAATT^C on the bottom strand");
        });
    }

    /// <summary>
    /// Circular molecules: the recognition window may wrap past the origin.
    /// Reference: Bio.Restriction.EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False) == [12]
    /// while the linear search returns []. (FormattedSeq.finditer, Restriction.py L215-228.)
    /// </summary>
    [Test]
    public void FindSites_EcoRI_Circular_FindsOriginSpanningSite()
    {
        const string plasmid = "TTCAAAAAAAGAA"; // GAATTC spans the origin starting at 10

        var linear = RestrictionAnalyzer.FindSites(new DnaSequence(plasmid), "EcoRI").ToList();
        var circular = RestrictionAnalyzer
            .FindSites(new DnaSequence(plasmid), "EcoRI", MoleculeTopology.Circular)
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(linear, Is.Empty, "no GAATTC lies wholly inside the linear molecule");
            Assert.That(circular, Has.Count.EqualTo(1), "one origin-spanning GAATTC");
            Assert.That(circular[0].Position, Is.EqualTo(10), "the site starts 3 bases before the origin");
            Assert.That(circular[0].CutPosition, Is.EqualTo(11),
                "Biopython EcoRI.search(..., linear=False) == [12] (1-based)");
            Assert.That(circular[0].RecognizedSequence, Is.EqualTo("GAATTC"),
                "the wrapped window reads GAATTC");
        });
    }

    /// <summary>
    /// Reference: Bio.Restriction.BamHI.search(Seq("CCAAAAAAAAGGAT"), linear=False) == [12],
    /// linear == [].
    /// </summary>
    [Test]
    public void FindSites_BamHI_Circular_FindsOriginSpanningSite()
    {
        var circular = RestrictionAnalyzer
            .FindSites(new DnaSequence("CCAAAAAAAAGGAT"), "BamHI", MoleculeTopology.Circular)
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(circular, Has.Count.EqualTo(1));
            Assert.That(circular[0].Position, Is.EqualTo(10));
            Assert.That(circular[0].CutPosition, Is.EqualTo(11),
                "Biopython reports 12 (1-based) for the circular molecule");
        });
    }

    /// <summary>
    /// On a closed circle every bond exists, so the sites the linear filter rejects are real.
    /// Reference: Bio.Restriction.Sau3AI.search(Seq("GATCGATCGATC"), linear=False)
    /// == [1, 5, 9] (1-based) == cut positions 0, 4, 8 here.
    /// </summary>
    [Test]
    public void FindSites_Sau3AI_Circular_ReportsEveryGatc()
    {
        var cuts = RestrictionAnalyzer
            .FindSites(new DnaSequence("GATCGATCGATC"), "Sau3AI", MoleculeTopology.Circular)
            .Where(s => s.IsForwardStrand)
            .Select(s => s.CutPosition)
            .OrderBy(c => c)
            .ToList();

        cuts.Should().Equal([0, 4, 8],
            "Biopython's circular Sau3AI.search returns [1, 5, 9] (1-based)");
    }

    /// <summary>
    /// A cut that would fall after the last base is, on a circle, the cut at the origin, so cut
    /// positions are reduced modulo the length.
    /// Reference: Bio.Restriction.Sau3AI.search(Seq("AAAAAAAAGATC"), linear=False) == [9]
    /// (1-based) — top-strand cut 8 — while the bottom-strand cut GATC^ at 12 == 0 (mod 12).
    /// </summary>
    [Test]
    public void FindSites_Circular_CutPositionIsReducedModuloLength()
    {
        var sites = RestrictionAnalyzer
            .FindSites(new DnaSequence("AAAAAAAAGATC"), "Sau3AI", MoleculeTopology.Circular)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Single(s => s.IsForwardStrand).CutPosition, Is.EqualTo(8),
                "Biopython circular Sau3AI.search == [9] (1-based)");
            Assert.That(sites.Single(s => !s.IsForwardStrand).CutPosition, Is.EqualTo(0),
                "the bottom-strand cut at 12 is the bond at the origin on a 12-bp circle");
        });
    }

    /// <summary>
    /// Documented divergence from the oracle: a recognition site longer than the molecule cannot
    /// be present, so no site is reported. Biopython pads a circular sequence with
    /// `data[1:size]` unconditionally and therefore wraps a short molecule more than once —
    /// Bio.Restriction.ClaI.search(Seq("GATC"), linear=False) == [4] for the 6-bp site ATCGAT
    /// on a 4-bp circle. This library refuses that reading.
    /// </summary>
    [Test]
    public void FindSites_Circular_MoleculeShorterThanSite_ReportsNothing()
    {
        RestrictionAnalyzer.FindSites(new DnaSequence("GATC"), "ClaI", MoleculeTopology.Circular)
            .Should().BeEmpty("a 6-bp recognition site cannot occur on a 4-bp circle");
    }

    /// <summary>
    /// FindAllSites and CreateMap carry the topology through.
    /// </summary>
    [Test]
    public void FindAllSites_Circular_FindsOriginSpanningSite()
    {
        var circular = RestrictionAnalyzer
            .FindAllSites(new DnaSequence("TTCAAAAAAAGAA"), MoleculeTopology.Circular)
            .Where(s => s.IsForwardStrand && s.Enzyme.Name == "EcoRI")
            .ToList();

        circular.Should().HaveCount(1, "the origin-spanning EcoRI site is found by the all-enzyme scan");
        circular[0].CutPosition.Should().Be(11);
    }

    /// <summary>
    /// A reverse-strand record reports the window as read on the strand it was found on. For a
    /// degenerate palindromic enzyme the two strands' windows differ: HincII (GTYRAC) matching
    /// GTCAAC on the forward strand reads GTTGAC on the reverse strand.
    /// Reference: Bio.Restriction.HincII.search(Seq("AAGTCAACAA")) == [6] (1-based) == cut 5.
    /// </summary>
    [Test]
    public void FindSites_HincII_ReverseRecognizedSequenceIsReadOnItsOwnStrand()
    {
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence("AAGTCAACAA"), "HincII").ToList();

        Assert.Multiple(() =>
        {
            var fwd = sites.Single(s => s.IsForwardStrand);
            var rev = sites.Single(s => !s.IsForwardStrand);
            Assert.That(fwd.CutPosition, Is.EqualTo(5), "GTC^AAC — Biopython reports 6 (1-based)");
            Assert.That(fwd.RecognizedSequence, Is.EqualTo("GTCAAC"), "forward-strand window");
            Assert.That(rev.Position, Is.EqualTo(2), "same forward-strand position");
            Assert.That(rev.RecognizedSequence, Is.EqualTo("GTTGAC"),
                "the reverse strand of GTCAAC reads GTTGAC, which also matches GTYRAC");
        });
    }

    /// <summary>
    /// Custom non-palindromic enzyme, reverse strand. CdiI (CATC^G, blunt) is one of only two
    /// REBASE enzymes whose two in-site cuts are asymmetric (fwd + rev != recognition length).
    /// Reference: Bio.Restriction.CdiI — site CATCG, fst5 = 4, fst3 = -1 (so reverse cut 4),
    /// elucidate() == "CATC^_G", ovhg == 0; CdiI.search(Seq("AAACATCGAAACGATGAAA")) == [8, 13]
    /// (1-based) == top-strand cuts 7 and 12.
    /// </summary>
    [Test]
    public void FindSites_CustomNonPalindromicEnzyme_CdiI_ReverseStrandCutMatchesReference()
    {
        var cdiI = new RestrictionEnzyme("CdiI", "CATCG", 4, 4, "REBASE 404");
        var sites = RestrictionAnalyzer
            .FindSites(new DnaSequence("AAACATCGAAACGATGAAA"), cdiI)
            .ToList();

        Assert.Multiple(() =>
        {
            var fwd = sites.Single(s => s.IsForwardStrand);
            var rev = sites.Single(s => !s.IsForwardStrand);
            Assert.That(fwd.Position, Is.EqualTo(3));
            Assert.That(fwd.CutPosition, Is.EqualTo(7), "Biopython reports 8 (1-based)");
            Assert.That(fwd.RecognizedSequence, Is.EqualTo("CATCG"));
            Assert.That(rev.Position, Is.EqualTo(11), "CGATG on the forward strand is CATCG on the reverse");
            Assert.That(rev.CutPosition, Is.EqualTo(12),
                "CdiI is blunt (ovhg 0), so both cuts of the reverse-strand site are at 12 — "
                + "Biopython reports 13 (1-based); the pre-fix formula Position + CutPositionReverse gave 15");
            Assert.That(rev.RecognizedSequence, Is.EqualTo("CATCG"),
                "the reverse-strand window reads the enzyme's own recognition sequence");
        });
    }

    /// <summary>
    /// Custom non-palindromic enzyme with a 3-base 5' overhang. Reference:
    /// Bio.Restriction.SimI — site GGGTC, fst5 = 2, fst3 = 0 (reverse cut 5),
    /// elucidate() == "GG^GTC_N", ovhg == -3; SimI.search(Seq("AAAGGGTCAAAGACCCAAA")) == [6, 12]
    /// (1-based) == top-strand cuts 5 and 11. The reverse-strand record carries the cut on its
    /// own (bottom) strand, which the reference relates to the top-strand cut by
    /// crick = watson - ovhg = 11 + 3 = 14.
    /// </summary>
    [Test]
    public void FindSites_CustomNonPalindromicEnzyme_SimI_ReverseStrandCutMatchesReference()
    {
        var simI = new RestrictionEnzyme("SimI", "GGGTC", 2, 5, "REBASE 404");
        var sites = RestrictionAnalyzer
            .FindSites(new DnaSequence("AAAGGGTCAAAGACCCAAA"), simI)
            .ToList();

        Assert.Multiple(() =>
        {
            var fwd = sites.Single(s => s.IsForwardStrand);
            var rev = sites.Single(s => !s.IsForwardStrand);
            Assert.That(fwd.Position, Is.EqualTo(3));
            Assert.That(fwd.CutPosition, Is.EqualTo(5), "GG^GTC — Biopython reports 6 (1-based)");
            Assert.That(rev.Position, Is.EqualTo(11));
            Assert.That(rev.CutPosition, Is.EqualTo(14),
                "bottom-strand cut of the reverse-strand site = watson 11 - ovhg (-3); "
                + "the pre-fix formula Position + CutPositionReverse gave 16");
            Assert.That(rev.RecognizedSequence, Is.EqualTo("GGGTC"));
        });
    }

    /// <summary>
    /// Regression guard for the built-in table: every entry is palindromic with
    /// CutPositionForward + CutPositionReverse == RecognitionLength, which is exactly the
    /// condition under which the corrected reverse-strand cut formula
    /// (Position + length - CutPositionForward) equals the historical Position +
    /// CutPositionReverse — so no built-in enzyme's output changed.
    /// </summary>
    [Test]
    public void EnzymeDatabase_EveryEntry_HasSymmetricCutGeometry()
    {
        foreach (var enzyme in RestrictionAnalyzer.Enzymes.Values)
        {
            (enzyme.CutPositionForward + enzyme.CutPositionReverse).Should().Be(
                enzyme.RecognitionLength,
                $"{enzyme.Name} is palindromic, so its two cuts are symmetric about the site centre");
            DnaSequence.GetReverseComplementString(enzyme.RecognitionSequence).Should().Be(
                enzyme.RecognitionSequence,
                $"{enzyme.Name}'s recognition sequence is its own reverse complement");
        }
    }

    /// <summary>
    /// A custom enzyme whose recognition sequence is not a usable IUPAC pattern is rejected at
    /// the call with a clear <see cref="ArgumentException"/>, instead of either throwing
    /// <see cref="ArgumentOutOfRangeException"/> from inside the IUPAC matcher (lower case,
    /// non-IUPAC letters) or reporting one "site" per position (empty recognition sequence).
    /// </summary>
    [Test]
    [TestCase("", TestName = "FindSites_CustomEnzyme_EmptyRecognition_Throws")]
    [TestCase("GAXTTC", TestName = "FindSites_CustomEnzyme_NonIupacRecognition_Throws")]
    [TestCase("GAATTC U", TestName = "FindSites_CustomEnzyme_SpaceInRecognition_Throws")]
    public void FindSites_CustomEnzyme_InvalidRecognitionSequence_ThrowsArgumentException(string recognition)
    {
        var enzyme = new RestrictionEnzyme("Custom", recognition, 1, 5, "test");

        Assert.Throws<ArgumentException>(
            () => RestrictionAnalyzer.FindSites(new DnaSequence("AAAGAATTCAAA"), enzyme),
            "an unusable recognition sequence is an argument error, reported eagerly");
    }

    /// <summary>
    /// Recognition sequences are matched case-insensitively, like the sequence itself.
    /// </summary>
    [Test]
    public void FindSites_CustomEnzyme_LowerCaseRecognitionSequence_StillMatches()
    {
        var enzyme = new RestrictionEnzyme("Custom", "gaattc", 1, 5, "test");

        var sites = RestrictionAnalyzer.FindSites(new DnaSequence("AAAGAATTCAAA"), enzyme)
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Position, Is.EqualTo(3));
            Assert.That(sites[0].CutPosition, Is.EqualTo(4));
            Assert.That(sites[0].RecognizedSequence, Is.EqualTo("GAATTC"),
                "the matched window is reported as read, upper case");
        });
    }

    /// <summary>
    /// The raw-string overload validates the enzyme name eagerly, at the call, and does so even
    /// when the sequence is empty — previously an unknown name surfaced only on enumeration and
    /// an empty sequence suppressed the check entirely.
    /// </summary>
    [Test]
    public void FindSites_RawString_UnknownEnzyme_ThrowsAtCallTime()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites("ACGTACGT", "NoSuchEnzyme"));
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites("", "NoSuchEnzyme"));
            Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.FindSites("ACGTACGT", ""));
        });
    }

    #endregion

    #region CreateMap Site Reporting (RESTR-FIND-001 review, 2026-10)

    /// <summary>
    /// NonCutters is taken from the considered enzyme set and compared case-insensitively against
    /// the canonical table names, so an enzyme that cut can never also be listed as a non-cutter.
    /// Before the fix, CreateMap(seq, "ecori") listed EcoRI under SitesByEnzyme *and* under
    /// NonCutters.
    /// </summary>
    [Test]
    public void CreateMap_LowerCaseEnzymeNames_DoesNotListACutterAsNonCutter()
    {
        var map = RestrictionAnalyzer.CreateMap(new DnaSequence("AAAGAATTCAAA"), "ecori", "bamhi");

        Assert.Multiple(() =>
        {
            Assert.That(map.SitesByEnzyme.ContainsKey("EcoRI"), Is.True, "EcoRI cut the sequence");
            Assert.That(map.NonCutters, Does.Not.Contain("EcoRI"), "a cutter is not a non-cutter");
            Assert.That(map.NonCutters, Does.Not.Contain("ecori"), "nor under the caller's spelling");
            Assert.That(map.NonCutters, Is.EquivalentTo(new[] { "BamHI" }),
                "non-cutters are reported under their canonical table name");
            Assert.That(map.UniqueCutters, Is.EquivalentTo(new[] { "EcoRI" }));
        });
    }

    /// <summary>
    /// With no enzyme names the considered set is the whole built-in table (the documented
    /// behaviour of this overload), so NonCutters lists every built-in enzyme without a site —
    /// previously it was always empty. Invariant: SitesByEnzyme ∪ NonCutters covers the
    /// considered set exactly (TestSpec RESTR-DIGEST-001 invariant #11).
    /// </summary>
    [Test]
    public void CreateMap_NoEnzymesSpecified_NonCuttersCoverTheRestOfTheTable()
    {
        var map = RestrictionAnalyzer.CreateMap(new DnaSequence("AAAGAATTCAAA"));
        int tableSize = RestrictionAnalyzer.Enzymes.Values
            .Select(e => e.Name)
            .Distinct(StringComparer.OrdinalIgnoreCase)
            .Count();

        Assert.Multiple(() =>
        {
            Assert.That(map.SitesByEnzyme.Keys, Is.EquivalentTo(new[] { "EcoRI" }));
            Assert.That(map.NonCutters, Does.Not.Contain("EcoRI"));
            Assert.That(map.NonCutters, Has.Count.EqualTo(tableSize - 1),
                "every built-in enzyme other than EcoRI is a non-cutter here");
            Assert.That(map.SitesByEnzyme.Count + map.NonCutters.Count, Is.EqualTo(tableSize),
                "the considered set is partitioned into cutters and non-cutters");
        });
    }

    /// <summary>
    /// A name repeated (in any casing) is considered once, so its sites are not double-counted.
    /// </summary>
    [Test]
    public void CreateMap_RepeatedEnzymeName_IsConsideredOnce()
    {
        var map = RestrictionAnalyzer.CreateMap(new DnaSequence("AAAGAATTCAAA"), "EcoRI", "ecori", "EcoRI");

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(1), "one forward-strand site, counted once");
            Assert.That(map.SitesByEnzyme["EcoRI"], Is.EquivalentTo(new[] { 3, 3 }),
                "SitesByEnzyme holds both strands' records, so the one palindromic site appears twice");
            Assert.That(map.UniqueCutters, Is.EquivalentTo(new[] { "EcoRI" }));
        });
    }

    /// <summary>
    /// An unknown enzyme name is rejected eagerly rather than surfacing from inside the site scan.
    /// </summary>
    [Test]
    public void CreateMap_UnknownEnzymeName_ThrowsArgumentException()
    {
        Assert.Throws<ArgumentException>(
            () => RestrictionAnalyzer.CreateMap(new DnaSequence("AAAGAATTCAAA"), "NoSuchEnzyme"));
    }

    /// <summary>
    /// CreateMap carries the topology through to the site search.
    /// Reference: Bio.Restriction.EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False) == [12].
    /// </summary>
    [Test]
    public void CreateMap_Circular_FindsOriginSpanningSite()
    {
        var linear = RestrictionAnalyzer.CreateMap(new DnaSequence("TTCAAAAAAAGAA"), "EcoRI");
        var circular = RestrictionAnalyzer.CreateMap(
            new DnaSequence("TTCAAAAAAAGAA"), MoleculeTopology.Circular, "EcoRI");

        Assert.Multiple(() =>
        {
            Assert.That(linear.TotalSites, Is.EqualTo(0));
            Assert.That(linear.NonCutters, Is.EquivalentTo(new[] { "EcoRI" }));
            Assert.That(circular.TotalSites, Is.EqualTo(1));
            Assert.That(circular.NonCutters, Is.Empty);
            Assert.That(circular.SitesByEnzyme["EcoRI"], Is.EquivalentTo(new[] { 10, 10 }),
                "both strands' records of the one origin-spanning site");
        });
    }

    #endregion

    #region FindPatternOccurrences (raw IUPAC pattern scan)

    /// <summary>
    /// The raw scan reports 0-based start indices of the leftmost base of each match, ascending.
    /// </summary>
    [Test]
    public void FindPatternOccurrences_ReportsZeroBasedStartIndices()
    {
        RestrictionAnalyzer.FindPatternOccurrences("AAAGAATTCAAAGAATTC", "GAATTC")
            .Should().Equal(3, 12);
    }

    /// <summary>
    /// Overlapping occurrences are all reported, matching the enzyme search (and Biopython,
    /// whose site regexes are lookaheads: Bio.Restriction.PacI.search(Seq("TTAATTAATTAA"))
    /// == [6, 10], i.e. sites starting at 0 and 4).
    /// </summary>
    [Test]
    public void FindPatternOccurrences_ReportsOverlappingOccurrences()
    {
        RestrictionAnalyzer.FindPatternOccurrences("TTAATTAATTAA", "TTAATTAA")
            .Should().Equal([0, 4], "the two TTAATTAA windows overlap by four bases");
    }

    /// <summary>
    /// Ambiguity is interpreted in the pattern, through the canonical
    /// <see cref="IupacHelper.MatchesIupac"/>.
    /// </summary>
    [Test]
    [TestCase("AAGTCGACAA", 2, TestName = "FindPatternOccurrences_GtyracMatchesGtcgac")]
    [TestCase("AAGTTAACAA", 2, TestName = "FindPatternOccurrences_GtyracMatchesGttaac")]
    [TestCase("AAGTCAACAA", 2, TestName = "FindPatternOccurrences_GtyracMatchesGtcaac")]
    [TestCase("AAGTTGACAA", 2, TestName = "FindPatternOccurrences_GtyracMatchesGttgac")]
    public void FindPatternOccurrences_DegeneratePattern_MatchesEveryExpansion(string sequence, int expected)
    {
        RestrictionAnalyzer.FindPatternOccurrences(sequence, "GTYRAC").Should().Equal(expected);
    }

    /// <summary>
    /// GTYRAC must not match an expansion outside its own code set.
    /// </summary>
    [Test]
    public void FindPatternOccurrences_DegeneratePattern_RejectsNonMatchingBases()
    {
        RestrictionAnalyzer.FindPatternOccurrences("AAGTAGACAA", "GTYRAC")
            .Should().BeEmpty("Y is C or T, so GTAGAC does not match");
    }

    /// <summary>
    /// An ambiguity code in the <b>sequence</b> blocks a concrete pattern base rather than
    /// widening it — Biopython 1.88 compiles a site's A/C/G/T to a literal, so
    /// <c>EcoRI.search(Seq("AGAANTCA")) == []</c> and <c>EcoRI.search(Seq("AGAARTCA")) == []</c>.
    /// Only a <b>pattern</b> <c>N</c> accepts a sequence ambiguity code (see the D3 tests below).
    /// </summary>
    [Test]
    [TestCase("AAAGAANTCAAA", TestName = "FindPatternOccurrences_SequenceN_BlocksMatch")]
    [TestCase("AAAGAARTCAAA", TestName = "FindPatternOccurrences_SequenceR_BlocksMatch")]
    public void FindPatternOccurrences_AmbiguousSequenceBase_NeverMatches(string sequence)
    {
        RestrictionAnalyzer.FindPatternOccurrences(sequence, "GAATTC")
            .Should().BeEmpty("a sequence position that is not A, C, G or T matches no pattern code");
    }

    #region D3 — a site N matches a sequence N (Biopython 1.88 Bio.Restriction)

    // Oracle: Biopython 1.88, run 2026-10-14. Bio.Restriction compiles a site's N to the regex "."
    // and every other code to the class of its A/C/G/T expansion (SfiI.compsite ==
    // "(?=(?P<SfiI>GGCC.....GGCC))", HincII.compsite == "(?=(?P<HincII>GT[CT][AG]AC))"), and the
    // searched sequence keeps its IUPAC letters. Search results below are Biopython's 1-based
    // top-strand cuts; Seqeron's CutPosition is the same bond 0-based.

    /// <summary>
    /// <c>SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]</c>: the five site N's accept the five sequence
    /// N's. Before D3 every sequence N blocked the site and this returned nothing.
    /// </summary>
    [TestCase(MoleculeTopology.Linear)]
    [TestCase(MoleculeTopology.Circular)]
    public void FindSites_SfiI_SequenceNInsideSpacer_MatchesBiopython(MoleculeTopology topology)
    {
        var sites = RestrictionAnalyzer.FindSites("GGCCNNNNNGGCC", "SfiI", topology).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Where(s => s.IsForwardStrand).Select(s => (s.Position, s.CutPosition)),
                Is.EqualTo(new[] { (0, 8) }), "Biopython SfiI.search == [9] (1-based)");
            Assert.That(sites.Where(s => !s.IsForwardStrand).Select(s => (s.Position, s.CutPosition)),
                Is.EqualTo(new[] { (0, 5) }), "bottom-strand cut GGCCN_NNN^NGGCC = 13 - 8");
            Assert.That(sites.Select(s => s.RecognizedSequence), Is.All.EqualTo("GGCCNNNNNGGCC"));
        });
    }

    /// <summary>
    /// Flanked, lower-case and origin-spanning variants, each locked to Biopython 1.88:
    /// <c>SfiI.search(Seq("AGGCCNNNNNGGCCA")) == [10]</c>,
    /// <c>SfiI.search(Seq("ggccnnnnnggcc")) == [9]</c>,
    /// <c>SfiI.search(Seq("NNNNNGGCCAGGCC"), linear=False) == [5]</c>.
    /// </summary>
    [TestCase("AGGCCNNNNNGGCCA", MoleculeTopology.Linear, 9)]
    [TestCase("ggccnnnnnggcc", MoleculeTopology.Linear, 8)]
    [TestCase("NNNNNGGCCAGGCC", MoleculeTopology.Circular, 4)]
    public void FindSites_SfiI_SequenceN_TopStrandCutsMatchBiopython(
        string sequence, MoleculeTopology topology, int expectedCut)
    {
        RestrictionAnalyzer.FindSites(sequence, "SfiI", topology)
            .Where(s => s.IsForwardStrand).Select(s => s.CutPosition)
            .Should().Equal(expectedCut);
    }

    /// <summary>
    /// A sequence ambiguity code still blocks a concrete or degenerate site code:
    /// <c>EcoRI.search(Seq("AGAANTCA")) == []</c>, <c>HincII.search(Seq("AGTYRACA")) == []</c>
    /// (Y and R are not in the regex classes [CT] and [AG]), <c>HincII.search(Seq("AGTNNACA")) == []</c>,
    /// and an N in a fixed SfiI half: <c>SfiI.search(Seq("GGCNNNNNNGGCC")) == []</c>.
    /// </summary>
    [TestCase("AGAANTCA", "EcoRI")]
    [TestCase("AGTYRACA", "HincII")]
    [TestCase("AGTNNACA", "HincII")]
    [TestCase("GGCNNNNNNGGCC", "SfiI")]
    [TestCase("NNNNNNNNNNNNNNNN", "EcoRI")]
    public void FindSites_SequenceAmbiguityAtConcreteSiteCode_StillBlocks(string sequence, string enzyme)
    {
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
            RestrictionAnalyzer.FindSites(sequence, enzyme, topology).Should().BeEmpty(
                $"Biopython {enzyme}.search reports no site ({topology})");
    }

    /// <summary>
    /// Any of the 15 IUPAC codes in the spacer is accepted by the site N's
    /// (<c>SfiI.search(Seq("GGCCRYSWKGGCC")) == [9]</c>), and the raw scan follows the same rule.
    /// </summary>
    [Test]
    public void SiteN_AcceptsEveryIupacCode_InFindSitesAndFindPatternOccurrences()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindSites("GGCCRYSWKGGCC", "SfiI")
                .Where(s => s.IsForwardStrand).Select(s => s.CutPosition), Is.EqualTo(new[] { 8 }));
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("GGCCNNNNNGGCC", "GGCCNNNNNGGCC"),
                Is.EqualTo(new[] { 0 }));
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("AATCNGAA", "TCNGA"),
                Is.EqualTo(new[] { 2 }), "Hpy188I compsite TC.GA matches TCNGA");
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("AAGTNNACAA", "GTNNAC"),
                Is.EqualTo(new[] { 2 }), "Hpy8I compsite GT..AC matches GTNNAC");
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("AAGTNNACAA", "GTMKAC"),
                Is.Empty, "AccI compsite GT[AC][GT]AC does not match GTNNAC");
        });
    }

    /// <summary>
    /// Deliberate narrowing (Evidence §3): a character outside the 15 IUPAC codes matches nothing,
    /// not even a site N. Biopython lets "." match any letter, so
    /// <c>SfiI.search(Seq("GGCCXXXXXGGCC")) == [9]</c>; this library reports no site.
    /// </summary>
    [Test]
    public void FindSites_NonIupacCharacterInSpacer_MatchesNothing()
    {
        RestrictionAnalyzer.FindSites("GGCCXXXXXGGCC", "SfiI").Should().BeEmpty();
    }

    #endregion

    /// <summary>
    /// The sequence that <see cref="DnaSequence"/> would reject is exactly the case this overload
    /// exists for: a pattern may be scanned against a string still carrying IUPAC codes.
    /// </summary>
    [Test]
    public void FindPatternOccurrences_AcceptsSequenceThatDnaSequenceRejects()
    {
        const string ambiguous = "GAATTCNNNGAATTC";

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => _ = new DnaSequence(ambiguous));
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences(ambiguous, "GAATTC"),
                Is.EqualTo(new[] { 0, 9 }));
        });
    }

    /// <summary>
    /// Both the sequence and the pattern are matched case-insensitively.
    /// </summary>
    [Test]
    public void FindPatternOccurrences_IsCaseInsensitive()
    {
        RestrictionAnalyzer.FindPatternOccurrences("aaagaattcaaa", "gaattc").Should().Equal(3);
    }

    /// <summary>
    /// startIndex skips occurrences that start before it (the semantics
    /// <c>CodonOptimizer.FindPattern(dna, pattern, from)</c> needs).
    /// </summary>
    [Test]
    [TestCase(0, new[] { 3, 12 })]
    [TestCase(3, new[] { 3, 12 })]
    [TestCase(4, new[] { 12 })]
    [TestCase(13, new int[0])]
    [TestCase(1000, new int[0])]
    public void FindPatternOccurrences_StartIndex_SkipsEarlierOccurrences(int startIndex, int[] expected)
    {
        RestrictionAnalyzer.FindPatternOccurrences("AAAGAATTCAAAGAATTC", "GAATTC", startIndex)
            .Should().Equal(expected);
    }

    /// <summary>
    /// An empty sequence, and a sequence shorter than the pattern, yield no occurrences.
    /// </summary>
    [Test]
    [TestCase("")]
    [TestCase("GAATT")]
    public void FindPatternOccurrences_SequenceTooShort_ReturnsEmpty(string sequence)
    {
        RestrictionAnalyzer.FindPatternOccurrences(sequence, "GAATTC").Should().BeEmpty();
    }

    /// <summary>
    /// Argument validation. The pattern must be a non-empty IUPAC string; the sequence must not
    /// be null; startIndex must not be negative.
    /// </summary>
    [Test]
    public void FindPatternOccurrences_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.FindPatternOccurrences(null!, "GAATTC"));
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindPatternOccurrences("ACGT", null!));
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindPatternOccurrences("ACGT", ""));
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindPatternOccurrences("ACGT", "GAXTTC"),
                "X is not an IUPAC nucleotide code");
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindPatternOccurrences("ACGT", "GAAUUC"),
                "U is not an IUPAC DNA code — callers map U to T themselves");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RestrictionAnalyzer.FindPatternOccurrences("ACGT", "ACGT", -1));
        });
    }

    /// <summary>
    /// The validation happens at the call, not on first enumeration, so a caller that stores the
    /// enumerable cannot be surprised later.
    /// </summary>
    [Test]
    public void FindPatternOccurrences_ValidatesEagerly()
    {
        Assert.Throws<ArgumentException>(
            () => RestrictionAnalyzer.FindPatternOccurrences("ACGTACGT", "GAXTTC"),
            "the exception is thrown by the call itself, before any enumeration");
    }

    /// <summary>
    /// No duplication: the raw scan is the same scan the enzyme search uses, so for every
    /// built-in enzyme it reproduces the forward-strand positions FindSites reports (which, for a
    /// linear molecule, are the feasible subset of the raw occurrences).
    /// </summary>
    [Test]
    public void FindPatternOccurrences_AgreesWithFindSitesForwardPositions()
    {
        const string sequence =
            "GGATCCAAAGAATTCCCCGGGAAGCTTGTCGACGATATCGGTACCAAAGGCCAAAAAGGCCAAATTAATTAAGATCAAA";

        foreach (var enzyme in RestrictionAnalyzer.Enzymes.Values)
        {
            var raw = RestrictionAnalyzer.FindPatternOccurrences(sequence, enzyme.RecognitionSequence).ToList();
            var found = RestrictionAnalyzer.FindSites(new DnaSequence(sequence), enzyme.Name)
                .Where(s => s.IsForwardStrand)
                .Select(s => s.Position)
                .ToList();

            found.Should().BeSubsetOf(raw,
                $"{enzyme.Name}: every reported forward-strand site is a raw pattern occurrence");

            var feasible = raw
                .Where(p => p + enzyme.CutPositionForward >= 1
                            && p + enzyme.CutPositionForward <= sequence.Length - 1
                            && p + enzyme.CutPositionReverse >= 1
                            && p + enzyme.CutPositionReverse <= sequence.Length - 1)
                .ToList();
            found.Should().Equal(feasible,
                $"{enzyme.Name}: the reported sites are exactly the occurrences whose two cuts exist");
        }
    }

    #endregion

    #region Enzyme Table Audit against REBASE 404 (RESTR-FIND-001 review, 2026-10)

    /// <summary>
    /// Every entry of the built-in table, audited against REBASE EMBOSS files <b>version 404
    /// (2024)</b> as shipped in Biopython 1.88 <c>Bio/Restriction/Restriction_Dictionary.py</c>
    /// (the provenance is stated in that generated file's own header; rebase.neb.com itself is
    /// not reachable from this environment). For each enzyme the expected values are
    /// <c>site</c>, <c>fst5</c> (= CutPositionForward) and <c>len(site) + fst3</c>
    /// (= CutPositionReverse), with the overhang class taken from the reference's own
    /// <c>is_blunt()</c> / <c>is_5overhang()</c> / <c>is_3overhang()</c>. All 39 entries matched
    /// on all four fields when this test was written; the Description carries the reference's
    /// <c>elucidate()</c> rendering of the cut.
    /// </summary>
    [Test]
    [TestCase("EcoRI", "GAATTC", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: G^AATT_C")]
    [TestCase("BamHI", "GGATCC", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: G^GATC_C")]
    [TestCase("HindIII", "AAGCTT", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: A^AGCT_T")]
    [TestCase("XhoI", "CTCGAG", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: C^TCGA_G")]
    [TestCase("SalI", "GTCGAC", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: G^TCGA_C")]
    [TestCase("PstI", "CTGCAG", 5, 1, OverhangType.ThreePrime, Description = "REBASE 404: C_TGCA^G")]
    [TestCase("SphI", "GCATGC", 5, 1, OverhangType.ThreePrime, Description = "REBASE 404: G_CATG^C")]
    [TestCase("KpnI", "GGTACC", 5, 1, OverhangType.ThreePrime, Description = "REBASE 404: G_GTAC^C")]
    [TestCase("SacI", "GAGCTC", 5, 1, OverhangType.ThreePrime, Description = "REBASE 404: G_AGCT^C")]
    [TestCase("XbaI", "TCTAGA", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: T^CTAG_A")]
    [TestCase("NcoI", "CCATGG", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: C^CATG_G")]
    [TestCase("NdeI", "CATATG", 2, 4, OverhangType.FivePrime, Description = "REBASE 404: CA^TA_TG")]
    [TestCase("NheI", "GCTAGC", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: G^CTAG_C")]
    [TestCase("SpeI", "ACTAGT", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: A^CTAG_T")]
    [TestCase("AvrII", "CCTAGG", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: C^CTAG_G")]
    [TestCase("ClaI", "ATCGAT", 2, 4, OverhangType.FivePrime, Description = "REBASE 404: AT^CG_AT")]
    [TestCase("AgeI", "ACCGGT", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: A^CCGG_T")]
    [TestCase("NotI", "GCGGCCGC", 2, 6, OverhangType.FivePrime, Description = "REBASE 404: GC^GGCC_GC")]
    [TestCase("SfiI", "GGCCNNNNNGGCC", 8, 5, OverhangType.ThreePrime, Description = "REBASE 404: GGCCN_NNN^NGGCC")]
    [TestCase("PacI", "TTAATTAA", 5, 3, OverhangType.ThreePrime, Description = "REBASE 404: TTA_AT^TAA")]
    [TestCase("AscI", "GGCGCGCC", 2, 6, OverhangType.FivePrime, Description = "REBASE 404: GG^CGCG_CC")]
    [TestCase("FseI", "GGCCGGCC", 6, 2, OverhangType.ThreePrime, Description = "REBASE 404: GG_CCGG^CC")]
    [TestCase("SwaI", "ATTTAAAT", 4, 4, OverhangType.Blunt, Description = "REBASE 404: ATTT^_AAAT")]
    [TestCase("MspI", "CCGG", 1, 3, OverhangType.FivePrime, Description = "REBASE 404: C^CG_G")]
    [TestCase("HpaII", "CCGG", 1, 3, OverhangType.FivePrime, Description = "REBASE 404: C^CG_G")]
    [TestCase("TaqI", "TCGA", 1, 3, OverhangType.FivePrime, Description = "REBASE 404: T^CG_A")]
    [TestCase("AluI", "AGCT", 2, 2, OverhangType.Blunt, Description = "REBASE 404: AG^_CT")]
    [TestCase("RsaI", "GTAC", 2, 2, OverhangType.Blunt, Description = "REBASE 404: GT^_AC")]
    [TestCase("HaeIII", "GGCC", 2, 2, OverhangType.Blunt, Description = "REBASE 404: GG^_CC")]
    [TestCase("DpnI", "GATC", 2, 2, OverhangType.Blunt, Description = "REBASE 404: GA^_TC")]
    [TestCase("Sau3AI", "GATC", 0, 4, OverhangType.FivePrime, Description = "REBASE 404: N^GATC_N")]
    [TestCase("MboI", "GATC", 0, 4, OverhangType.FivePrime, Description = "REBASE 404: N^GATC_N")]
    [TestCase("EcoRV", "GATATC", 3, 3, OverhangType.Blunt, Description = "REBASE 404: GAT^_ATC")]
    [TestCase("SmaI", "CCCGGG", 3, 3, OverhangType.Blunt, Description = "REBASE 404: CCC^_GGG")]
    [TestCase("HincII", "GTYRAC", 3, 3, OverhangType.Blunt, Description = "REBASE 404: GTY^_RAC")]
    [TestCase("ScaI", "AGTACT", 3, 3, OverhangType.Blunt, Description = "REBASE 404: AGT^_ACT")]
    [TestCase("StuI", "AGGCCT", 3, 3, OverhangType.Blunt, Description = "REBASE 404: AGG^_CCT")]
    [TestCase("BglII", "AGATCT", 1, 5, OverhangType.FivePrime, Description = "REBASE 404: A^GATC_T")]
    [TestCase("ApaI", "GGGCCC", 5, 1, OverhangType.ThreePrime, Description = "REBASE 404: G_GGCC^C")]
    public void EnzymeDatabase_AllEntries_MatchRebase404(
        string name, string site, int cutForward, int cutReverse, OverhangType overhang)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(name);

        Assert.Multiple(() =>
        {
            Assert.That(enzyme, Is.Not.Null, $"{name} should be in the built-in table");
            Assert.That(enzyme!.RecognitionSequence, Is.EqualTo(site), $"{name} recognition sequence (REBASE 404)");
            Assert.That(enzyme.CutPositionForward, Is.EqualTo(cutForward), $"{name} fst5 (REBASE 404)");
            Assert.That(enzyme.CutPositionReverse, Is.EqualTo(cutReverse), $"{name} len(site) + fst3 (REBASE 404)");
            Assert.That(enzyme.OverhangType, Is.EqualTo(overhang), $"{name} overhang class (REBASE 404)");
            Assert.That(enzyme.IsBluntEnd, Is.EqualTo(overhang == OverhangType.Blunt),
                $"{name} blunt/sticky classification (REBASE 404)");
        });
    }

    /// <summary>
    /// The audit above covers the whole table — guard against an entry being added without one.
    /// </summary>
    [Test]
    public void EnzymeDatabase_AuditCoversEveryEntry()
    {
        RestrictionAnalyzer.Enzymes.Should().HaveCount(39,
            "the REBASE 404 audit in EnzymeDatabase_AllEntries_MatchRebase404 has one case per entry");
    }

    #endregion
}
