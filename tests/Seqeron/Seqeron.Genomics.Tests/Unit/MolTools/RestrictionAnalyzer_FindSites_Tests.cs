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
    /// Evidence: REBASE. The built-in table holds exactly <b>39</b> common enzymes — the number the
    /// class doc, the three algorithm docs and every MCP tool description state — so the count is
    /// pinned exactly. It used to assert only <c>&gt;= 30</c>, under which nine enzymes could have
    /// vanished silently (review 2026-09 B08 H54).
    /// </summary>
    [Test]
    public void Enzymes_ContainsExactlyThe39CommonEnzymes()
    {
        var enzymes = RestrictionAnalyzer.Enzymes;

        Assert.Multiple(() =>
        {
            Assert.That(enzymes, Has.Count.EqualTo(39), "the documented size of the built-in table");
            Assert.That(enzymes.Keys.Order(StringComparer.Ordinal), Is.EqualTo(new[]
            {
                "AgeI", "AluI", "ApaI", "AscI", "AvrII", "BamHI", "BglII", "ClaI", "DpnI", "EcoRI",
                "EcoRV", "FseI", "HaeIII", "HincII", "HindIII", "HpaII", "KpnI", "MboI", "MspI",
                "NcoI", "NdeI", "NheI", "NotI", "PacI", "PstI", "RsaI", "SacI", "SalI", "Sau3AI",
                "ScaI", "SfiI", "SmaI", "SpeI", "SphI", "StuI", "SwaI", "TaqI", "XbaI", "XhoI",
            }), "the exact membership, so neither a loss nor an unannounced addition passes");

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
    /// Evidence: Wikipedia - Recognition sequences range from 4-8 bp. The invariant the filtering doc
    /// states is "returns <b>exactly</b> the enzymes with that length", so the expected membership is
    /// pinned rather than only soundness plus non-emptiness (review 2026-09 B08 H54). The four
    /// lengths present in the built-in table are 4 (9 enzymes), 6 (24), 8 (5) and 13 (SfiI, the one
    /// interrupted palindrome), which sum to the table's 39.
    /// </summary>
    [TestCase(4, new[] { "AluI", "DpnI", "HaeIII", "HpaII", "MboI", "MspI", "RsaI", "Sau3AI", "TaqI" })]
    [TestCase(6, new[]
    {
        "AgeI", "ApaI", "AvrII", "BamHI", "BglII", "ClaI", "EcoRI", "EcoRV", "HincII", "HindIII",
        "KpnI", "NcoI", "NdeI", "NheI", "PstI", "SacI", "SalI", "ScaI", "SmaI", "SpeI", "SphI",
        "StuI", "XbaI", "XhoI",
    })]
    [TestCase(8, new[] { "AscI", "FseI", "NotI", "PacI", "SwaI" })]
    [TestCase(13, new[] { "SfiI" })]
    [TestCase(5, new string[0])]
    [TestCase(7, new string[0])]
    public void GetEnzymesByCutLength_ReturnsExactlyTheEnzymesOfThatLength(int cutLength, string[] expected)
    {
        var enzymes = RestrictionAnalyzer.GetEnzymesByCutLength(cutLength).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(enzymes.Select(e => e.Name).Order(StringComparer.Ordinal), Is.EqualTo(expected),
                $"the exact set of {cutLength}-bp cutters");
            Assert.That(enzymes.All(e => e.RecognitionLength == cutLength),
                $"All enzymes should have recognition length of {cutLength}");
        });
    }

    /// <summary>
    /// And the per-length sets partition the whole table, so a length added or removed cannot hide
    /// between the cases above (B08 H54).
    /// </summary>
    [Test]
    public void GetEnzymesByCutLength_PerLengthSets_PartitionTheTable()
    {
        var byLength = RestrictionAnalyzer.Enzymes.Values
            .GroupBy(e => e.RecognitionLength)
            .ToDictionary(g => g.Key, g => g.Count());

        Assert.Multiple(() =>
        {
            Assert.That(byLength, Is.EquivalentTo(new Dictionary<int, int>
            {
                [4] = 9, [6] = 24, [8] = 5, [13] = 1,
            }));
            Assert.That(byLength.Values.Sum(), Is.EqualTo(39));
            Assert.That(Enumerable.Range(1, 20).Sum(l => RestrictionAnalyzer.GetEnzymesByCutLength(l).Count()),
                Is.EqualTo(39), "every built-in is returned by exactly one length");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - Enzymes produce blunt ends, 5' overhangs, or 3' overhangs. These two
    /// tests used to assert only <c>Is.Not.Empty</c> plus a soundness predicate, so an
    /// implementation returning just <c>AluI</c> would have passed both (review 2026-09 B08 H54).
    /// The counts are now pinned: 10 blunt and 29 sticky, which partition the table's 39. The exact
    /// blunt membership is locked separately in
    /// <c>RestrictionAnalyzer_Filter_Tests.GetBluntCutters_ReturnsExactlyTheDocumentedBluntSet</c>;
    /// 10 is also the common-set share of the 137 blunt REBASE cutters that Biopython's
    /// <c>EcoRV.compatible_end()</c> returns.
    /// </summary>
    [Test]
    public void GetBluntCutters_ReturnsOnlyBluntEndEnzymes()
    {
        var bluntCutters = RestrictionAnalyzer.GetBluntCutters().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(bluntCutters, Has.Count.EqualTo(10), "the 10 centre-cutting built-ins");
            Assert.That(bluntCutters.All(e => e.IsBluntEnd), "All should be blunt end");
            Assert.That(bluntCutters.All(e => e.OverhangType == OverhangType.Blunt),
                "All should have Blunt overhang type");
            Assert.That(bluntCutters.All(e => e.CutPositionForward == e.CutPositionReverse),
                "a blunt cutter cuts both strands at the same offset");
        });
    }

    /// <summary>
    /// Evidence: Wikipedia - Sticky end enzymes produce overhangs. See the remark above: the count is
    /// pinned so a soundness-only pass is no longer possible (B08 H54).
    /// </summary>
    [Test]
    public void GetStickyCutters_ReturnsOnlyStickyEndEnzymes()
    {
        var stickyCutters = RestrictionAnalyzer.GetStickyCutters().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(stickyCutters, Has.Count.EqualTo(29), "39 built-ins minus the 10 blunt ones");
            Assert.That(stickyCutters.All(e => !e.IsBluntEnd), "None should be blunt end");
            Assert.That(stickyCutters.All(e =>
                e.OverhangType == OverhangType.FivePrime || e.OverhangType == OverhangType.ThreePrime),
                "All should have 5' or 3' overhang type");
            Assert.That(
                stickyCutters.Select(e => e.Name)
                    .Concat(RestrictionAnalyzer.GetBluntCutters().Select(e => e.Name)),
                Is.EquivalentTo(RestrictionAnalyzer.Enzymes.Keys),
                "blunt and sticky partition the table exactly");
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
    /// Test: a null enzyme name throws ArgumentNullException and an empty one ArgumentException, eagerly - before
    /// anything is enumerated (review 2026-09 B08 H24; the null/empty split of H19).
    /// </summary>
    [Test]
    public void FindSites_NullEnzymeName_ThrowsArgumentNullException_EmptyThrowsArgumentException()
    {
        var sequence = new DnaSequence("ACGT");

        Assert.Throws<ArgumentNullException>(() => RestrictionAnalyzer.FindSites(sequence, (string)null!).ToList());
        Assert.Throws<ArgumentNullException>(() => _ = RestrictionAnalyzer.FindSites(sequence, (string)null!),
            "the guard is eager: constructing the query already throws");
        Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites(sequence, "").ToList());
        Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.FindSites(sequence, ""),
            "the guard is eager: constructing the query already throws");
    }

    #region Name resolution, empty lists and eagerness (B08 H13/H24)

    /// <summary>
    /// H13: REBASE spells the shipped enzyme <c>CviKI-1</c> (NEB R0710) with a hyphen, but the frozen
    /// table is keyed by Biopython's class-name spelling <c>CviKI_1</c>, because a Python class name
    /// cannot contain one. Before the fix <c>GetEnzyme("CviKI-1")</c> answered null and every
    /// name-taking member rejected the enzyme's real name. Reference: Bio.Restriction.CviKI_1 —
    /// site RGCY, blunt, <c>ovhgseq == ""</c>, and <c>search(Seq("AGCTAGCT"), linear=True) == [3, 7]</c>
    /// (1-based), i.e. the 0-based bonds 2 and 6. Measured with Biopython 1.88.
    /// The frozen data file is NOT changed; the consuming lookup accepts both spellings.
    /// </summary>
    [TestCase("CviKI-1")]
    [TestCase("CviKI_1")]
    [TestCase("cviki-1")]
    [TestCase("CVIKI_1")]
    public void GetEnzyme_RebaseHyphenSpelling_ResolvesToTheSameEnzyme(string spelling)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(spelling);

        Assert.That(enzyme, Is.Not.Null, $"'{spelling}' must resolve");
        Assert.Multiple(() =>
        {
            Assert.That(enzyme!.Name, Is.EqualTo("CviKI_1"), "the table's own spelling is reported back");
            Assert.That(enzyme.RecognitionSequence, Is.EqualTo("RGCY"));
            Assert.That(RestrictionAnalyzer.GetOverhangSequence(spelling), Is.Empty,
                "Biopython CviKI_1.ovhgseq == '' (blunt)");
            Assert.That(
                RestrictionAnalyzer.Digest(new DnaSequence("AGCTAGCT"), MoleculeTopology.Linear, spelling)
                    .Select(f => f.StartPosition).Where(p => p != 0).Distinct().Order(),
                Is.EqualTo(new[] { 2, 6 }),
                "Biopython CviKI_1.search(Seq(\"AGCTAGCT\"), linear=True) == [3, 7] (1-based)");
            Assert.That(RestrictionAnalyzer.AreCompatible(spelling, "AluI"), Is.True,
                "both are blunt cutters");
        });
    }

    /// <summary>
    /// H13: surrounding whitespace is ignored, as <c>MethylaseSets.FromName</c> already did.
    /// </summary>
    [TestCase(" EcoRI")]
    [TestCase("EcoRI ")]
    [TestCase("  ecori\t")]
    public void GetEnzyme_SurroundingWhitespace_IsIgnored(string spelling)
    {
        Assert.That(RestrictionAnalyzer.GetEnzyme(spelling)?.Name, Is.EqualTo("EcoRI"));
    }

    /// <summary>
    /// H13: a name that is genuinely unknown still answers null, and a hyphen is not a wildcard —
    /// only the underscore position the table actually uses resolves.
    /// </summary>
    [TestCase("NoSuchEnzyme")]
    [TestCase("Eco-RI")]
    [TestCase("Dpi3084I")]
    public void GetEnzyme_StillNullForAnUnresolvableName(string spelling)
    {
        Assert.That(RestrictionAnalyzer.GetEnzyme(spelling), Is.Null);
    }

    /// <summary>
    /// H24: an empty enzyme-name array means the same thing on every member that takes one —
    /// an error — rather than "no sites" on <c>FindSites</c>, a throw on the digest family and the
    /// methylation overload of <c>FindSites</c>, and "every common enzyme" on <c>CreateMap</c>.
    /// <c>CreateMap</c>'s name overload is the one deliberate exception and keeps its documented
    /// meaning, which the second half of this test pins.
    /// </summary>
    [Test]
    public void EmptyEnzymeNameArray_IsAnErrorEverywhereExceptCreateMap()
    {
        var seq = new DnaSequence("GAATTCGGATCC");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.FindSites(seq));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.FindSites(seq, MoleculeTopology.Linear));
            Assert.Throws<ArgumentException>(
                () => _ = RestrictionAnalyzer.FindSites(seq, MoleculeTopology.Linear, MethylaseSets.Dam));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.Digest(seq));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.Digest(seq, MoleculeTopology.Linear));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.PartialDigest(seq));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.GetDigestSummary(seq));

            // The documented exception: no names means the 39 common enzymes.
            var map = RestrictionAnalyzer.CreateMap(seq);
            Assert.That(map.SitesByEnzyme.Keys.Concat(map.NonCutters).Distinct().Count(), Is.EqualTo(39),
                "CreateMap with no names considers every common enzyme");
        });
    }

    /// <summary>
    /// H24: an unknown name in a <c>params string[]</c> is rejected at the call, not on first
    /// enumeration. Both entry points were iterators, so the query could be built and passed around
    /// before failing — although the class's own name-resolution contract and the
    /// <c>AreCompatible</c> remark both state the rejection is eager.
    /// </summary>
    [Test]
    public void UnknownEnzymeNameInAnArray_IsRejectedAtTheCall()
    {
        var seq = new DnaSequence("GAATTCGGATCC");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => _ = RestrictionAnalyzer.FindSites(seq, "NoSuchEnzyme1", "EcoRI"));
            Assert.Throws<ArgumentException>(
                () => _ = RestrictionAnalyzer.FindSites(seq, MoleculeTopology.Circular, "EcoRI", "NoSuchEnzyme1"));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.Digest(seq, "NoSuchEnzyme1"));
            Assert.Throws<ArgumentException>(
                () => _ = RestrictionAnalyzer.Digest(seq, MoleculeTopology.Circular, "EcoRI", "NoSuchEnzyme1"));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.PartialDigest(seq, "NoSuchEnzyme1"));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.GetDigestSummary(seq, "NoSuchEnzyme1"));
        });
    }

    /// <summary>
    /// H25: <c>FindPatternOccurrences</c> promised in its <c>&lt;param&gt;</c> doc that "a null or
    /// empty sequence yields no occurrences" while the code began with
    /// <c>ArgumentNullException.ThrowIfNull</c> — and its own <c>&lt;exception&gt;</c> tag documented
    /// the throw, so two doc statements disagreed with each other and one with the code. The sibling
    /// raw-string <c>FindSites(string, string)</c> makes the same promise and honours it, so the code
    /// now does too. The <i>pattern</i> is still validated eagerly, before the sequence is looked at.
    /// </summary>
    [Test]
    public void FindPatternOccurrences_NullOrEmptySequence_YieldsNoOccurrences()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences(null, "GAATTC"), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("", "GAATTC"), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindSites((string?)null!, "EcoRI"), Is.Empty,
                "the sibling already behaved this way");
            Assert.That(RestrictionAnalyzer.FindSites("", "EcoRI"), Is.Empty);

            // An invalid pattern is still an error, whatever the sequence is.
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.FindPatternOccurrences(null, ""));
            Assert.Throws<ArgumentException>(() => _ = RestrictionAnalyzer.FindPatternOccurrences(null, "GAZTC"));
            Assert.Throws<ArgumentOutOfRangeException>(
                () => _ = RestrictionAnalyzer.FindPatternOccurrences(null, "GAATTC", -1));
        });
    }

    #endregion

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
    /// of cut positions for the enzymes listed in the Wikipedia table.
    /// Source: https://en.wikipedia.org/wiki/Restriction_enzyme#Examples
    /// Attribution (B08 F64, QUEUE R9; page re-fetched 2026-10-15): the 19 cases labelled
    /// "Wikipedia" are on that table. BglII and ApaI are NOT on it (BglII appears only under
    /// "See also"); their site and cut positions come from Biopython 1.88 Bio.Restriction
    /// (site / fst5 / size + fst3, pinned REBASE snapshot — docs/Evidence/B08-Source-Provenance.md),
    /// which also reproduces all 19 Wikipedia rows. The method name is historical.
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
    [TestCase("BglII", "AGATCT", 1, 5, Description = "Biopython 1.88 Bio.Restriction (not on the Wikipedia table): A↓GATCT / TCTAG↓A")]
    [TestCase("ApaI", "GGGCCC", 5, 1, Description = "Biopython 1.88 Bio.Restriction (not on the Wikipedia table): GGGCC↓C / C↓CCGGG")]
    public void EnzymeDatabase_CutPositions_MatchWikipedia(
        string name, string expectedRecognition, int expectedCutFwd, int expectedCutRev)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(name);

        Assert.Multiple(() =>
        {
            Assert.That(enzyme, Is.Not.Null, $"{name} should be in database");
            Assert.That(enzyme!.RecognitionSequence, Is.EqualTo(expectedRecognition),
                $"{name} recognition sequence (Wikipedia table / Biopython 1.88)");
            Assert.That(enzyme.CutPositionForward, Is.EqualTo(expectedCutFwd),
                $"{name} forward cut position (Wikipedia table / Biopython 1.88)");
            Assert.That(enzyme.CutPositionReverse, Is.EqualTo(expectedCutRev),
                $"{name} reverse cut position (Wikipedia table / Biopython 1.88)");
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
    /// A recognition site longer than the molecule still occurs on a circle: the window reads past
    /// the origin and back in. The pre-H6 code applied the linear length guard before the topology
    /// dispatch and answered "no sites", with a comment claiming such a reading wraps the molecule
    /// "more than once" — ATCGAT on the 4-bp circle GATC starts at index 1 and needs only two bases
    /// past the origin, exactly Biopython's own <c>data + data[1:size]</c> extension.
    /// Reference: Bio.Restriction.ClaI.search(Seq("GATC"), linear=False) == [4] (1-based), i.e. the
    /// 0-based top-strand bond 3 (ClaI = AT^CGAT, so the site starts at 1 and cuts at 1 + 2).
    /// Measured with Biopython 1.88 (review 2026-09 B08 H6).
    /// </summary>
    [Test]
    public void FindSites_Circular_MoleculeShorterThanSite_FindsTheWrappedSite()
    {
        var sites = RestrictionAnalyzer
            .FindSites(new DnaSequence("GATC"), "ClaI", MoleculeTopology.Circular)
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1), "ATCGAT occurs once on the 4-bp circle GATC");
            Assert.That(sites[0].Position, Is.EqualTo(1), "the window starts one base past the origin");
            Assert.That(sites[0].RecognizedSequence, Is.EqualTo("ATCGAT"));
            Assert.That(sites[0].CutPosition, Is.EqualTo(3),
                "Biopython ClaI.search(Seq(\"GATC\"), linear=False) == [4] (1-based)");
        });
        Assert.That(RestrictionAnalyzer.FindSites(new DnaSequence("GATC"), "ClaI").ToList(), Is.Empty,
            "the linear molecule is still too short to hold the site");
    }

    /// <summary>
    /// The five enzyme/molecule rows H6 names: the only unexplained rows of a 72384-row Biopython
    /// digest differential before the fix, where Seqeron reported no site at all.
    /// Numbers are Bio.Restriction &lt;enzyme&gt;.search(Seq(molecule), linear=False) under
    /// Biopython 1.88, converted from its 1-based "first base of the 3' fragment" to the 0-based
    /// top-strand bond this library reports, and reduced modulo the molecule length (Biopython's
    /// circular <c>_drop</c> normalises only a prefix/suffix of its own result list, so it can
    /// return 13 for a 6-mer). Measured post-fix over all 754 REBASE enzymes × 84 molecules:
    /// 0 divergences at every length that can hold the site (review 2026-09 B08 H6).
    /// </summary>
    [TestCase("ACGT", "BsiWI", new[] { 2 })]
    [TestCase("ACGT", "Pfl23II", new[] { 2 })]
    [TestCase("ACGT", "PspLI", new[] { 2 })]
    [TestCase("ACGT", "SplI", new[] { 2 })]
    [TestCase("ACATGA", "MslI", new[] { 0 })]
    [TestCase("ACATGA", "DrdV", new[] { 0 })]
    [TestCase("ACATGA", "RseI", new[] { 0 })]
    [TestCase("ACATGA", "SmiMI", new[] { 0 })]
    [TestCase("TATTGCAAAA", "AlfI", new[] { 4, 8 })]
    public void Digest_Circular_MoleculeShorterThanSite_MatchesBiopython(
        string molecule, string enzymeName, int[] expectedCuts)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(enzymeName);
        Assert.That(enzyme, Is.Not.Null, $"{enzymeName} is a REBASE enzyme with a defined cut");
        Assert.That(molecule.Length, Is.LessThan(enzyme!.RecognitionSequence.Length),
            "the row exists precisely because the molecule is shorter than the recognition sequence");

        var sites = RestrictionAnalyzer
            .FindSites(new DnaSequence(molecule), MoleculeTopology.Circular, [enzyme])
            .ToList();
        Assert.That(sites, Is.Not.Empty, "the wrapped window is found");

        var cuts = RestrictionAnalyzer
            .Digest(new DnaSequence(molecule), MoleculeTopology.Circular, [enzyme])
            .Select(f => f.StartPosition)
            .Distinct()
            .OrderBy(c => c)
            .ToList();

        Assert.That(cuts, Is.EqualTo(expectedCuts),
            $"Biopython 1.88 {enzymeName}.search(Seq(\"{molecule}\"), linear=False) reduced mod {molecule.Length}");
    }

    /// <summary>
    /// Where the fix is a strict superset of Biopython: a window that needs more than one extra
    /// turn around the circle. Biopython's extension is <c>data + data[1:size]</c>, whose slice
    /// silently truncates to the molecule, so on the 4-bp circle ACGT it can only reach start
    /// indices 0..2 and never tests index 3 — yet SnaBI's TACGTA is the window at index 3 and is
    /// physically present on the circle. Measured: Bio.Restriction.SnaBI.search(Seq("ACGT"),
    /// linear=False) == [] with Biopython 1.88 (review 2026-09 B08 H6, F100).
    /// </summary>
    [Test]
    public void FindSites_Circular_ShortMolecule_ReportsWindowsBiopythonTruncatesAway()
    {
        var sites = RestrictionAnalyzer
            .FindSites(new DnaSequence("ACGT"), "SnaBI", MoleculeTopology.Circular)
            .Where(s => s.IsForwardStrand)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(1));
            Assert.That(sites[0].Position, Is.EqualTo(3), "the window Biopython's truncated slice cannot reach");
            Assert.That(sites[0].RecognizedSequence, Is.EqualTo("TACGTA"));
            Assert.That(sites[0].CutPosition, Is.EqualTo(2), "TAC^GTA, so 3 + 3 reduced mod 4");
        });
    }

    /// <summary>
    /// The empty molecule has no window on either topology, so the circular branch keeps the
    /// guard the H6 fix removed for every non-empty length.
    /// </summary>
    [Test]
    public void FindSites_Circular_EmptyMolecule_ReportsNothing()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindSites("", "EcoRI", MoleculeTopology.Circular), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindSites("", "EcoRI", MoleculeTopology.Linear), Is.Empty);
        });
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
            Assert.Throws<ArgumentException>(() => RestrictionAnalyzer.FindSites("ACGTACGT", ""));
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

    /// <summary>
    /// Review 2026-09 B08 F49 (D20): U+017F 'ſ' (LATIN SMALL LETTER LONG S) is not an IUPAC code
    /// either, so it matches nothing — even though <see cref="string.ToUpperInvariant"/> folds it
    /// onto <c>'S'</c>, which a site <c>N</c> or <c>S</c> would accept. All three raw-string entry
    /// points fold ASCII-only (<c>SequenceExtensions.ToUpperAscii</c>): the sequence scans see 'ſ'
    /// as a non-IUPAC character (F40 contract), and a pattern containing it is rejected like any
    /// other non-IUPAC pattern. Before F49 each of these three reported a match / accepted the
    /// pattern (measured in F43). Biopython 1.88 cannot even build the input
    /// (<c>Seq("GGCCſſſſſGGCC")</c> raises <c>UnicodeEncodeError</c>).
    /// </summary>
    [Test]
    public void RawStringEntryPoints_LongS_IsNotFoldedOntoIupacS()
    {
        Assert.Multiple(() =>
        {
            Assert.That("ſ".ToUpperInvariant(), Is.EqualTo("S"),
                "precondition: the invariant fold this test guards against really maps U+017F onto 'S'");
            Assert.That(RestrictionAnalyzer.FindSites("GGCCſſſſſGGCC", "SfiI"), Is.Empty,
                "a non-IUPAC spacer matches nothing, exactly like GGCCXXXXXGGCC");
            Assert.That(RestrictionAnalyzer.FindSites("GGCCSSSSSGGCC", "SfiI"), Is.Not.Empty,
                "control: the ASCII IUPAC code S in the spacer is matched by site N");
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("AſC", "ANC"), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("AſC", "ASC"), Is.Empty);
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("ASC", "ANC"), Is.EqualTo(new[] { 0 }),
                "control: the ASCII S is still matched");
        });
    }

    /// <summary>
    /// F49: a recognition pattern containing U+017F 'ſ' is rejected with the same
    /// <see cref="ArgumentException"/> as any other non-IUPAC pattern ('X', F40), on both the raw
    /// pattern scan and a custom enzyme; before F49 it was silently read as <c>ASC</c>.
    /// </summary>
    [Test]
    public void Pattern_WithLongS_IsRejectedLikeAnyNonIupacCharacter()
    {
        Assert.Multiple(() =>
        {
            var ex = Assert.Throws<ArgumentException>(
                () => RestrictionAnalyzer.FindPatternOccurrences("AAASCAAA", "AſC"));
            Assert.That(ex!.Message, Does.Contain("'ſ' at index 1"));
            Assert.Throws<ArgumentException>(
                () => RestrictionAnalyzer.FindPatternOccurrences("AAAXCAAA", "AXC"),
                "control: the same rejection as an ASCII non-IUPAC letter");
            Assert.Throws<ArgumentException>(
                () => RestrictionAnalyzer.FindSites(new DnaSequence("AAAGCCAAA"),
                    new RestrictionEnzyme("Custom", "GſC", 1, 2, "test")));
        });
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
    /// Argument validation. The pattern must be a non-empty IUPAC string and startIndex must not be
    /// negative. The <i>sequence</i> is not validated: a null or empty one yields no occurrences,
    /// which is what the parameter's own doc promises and what the sibling raw-string
    /// <c>FindSites(string, string)</c> has always done (review 2026-09 B08 H25 — this used to throw
    /// <c>ArgumentNullException</c>, contradicting both that doc and the sibling; pinned by
    /// <c>FindPatternOccurrences_NullOrEmptySequence_YieldsNoOccurrences</c>).
    /// </summary>
    [Test]
    public void FindPatternOccurrences_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences(null, "GAATTC"), Is.Empty);
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

    #region D7 — REBASE allenz.104 methylation-site field (RESTR-FIND-001 review, 2026-10)

    /// <summary>
    /// B08 D7 / F52. Every built-in entry's <see cref="RestrictionEnzyme.RebaseMethylationField"/>
    /// equals REBASE <c>allenz.104</c> field <c>&lt;6&gt;</c> ("METHYLATION SITE") for that enzyme.
    /// Expected strings were extracted programmatically from
    /// <c>crisprVerse/crisprBase/inst/enzymes/allenz.txt</c> (sha256
    /// <c>43b12432…69388</c>, see <c>docs/Evidence/B08-Source-Provenance.md</c> §3a) — the
    /// Description carries the record's <c>&lt;5&gt;</c> and <c>&lt;6&gt;</c> lines verbatim. The
    /// field is the site protected by the enzyme's <i>cognate</i> methyltransferase; it says
    /// nothing about Dam/Dcm/CpG sensitivity.
    /// </summary>
    [Test]
    [TestCase("EcoRI", "3(6)", Description = "allenz.104 <5>G^AATTC <6>3(6)")]
    [TestCase("BamHI", "5(4)", Description = "allenz.104 <5>G^GATCC <6>5(4)")]
    [TestCase("HindIII", "1(6)", Description = "allenz.104 <5>A^AGCTT <6>1(6)")]
    [TestCase("XhoI", "5(6)", Description = "allenz.104 <5>C^TCGAG <6>5(6)")]
    [TestCase("SalI", "5(6)", Description = "allenz.104 <5>G^TCGAC <6>5(6)")]
    [TestCase("PstI", "5(6)", Description = "allenz.104 <5>CTGCA^G <6>5(6)")]
    [TestCase("SphI", "", Description = "allenz.104 <5>GCATG^C <6>")]
    [TestCase("KpnI", "4(6)", Description = "allenz.104 <5>GGTAC^C <6>4(6)")]
    [TestCase("SacI", "4(5)", Description = "allenz.104 <5>GAGCT^C <6>4(5)")]
    [TestCase("XbaI", "6(6)", Description = "allenz.104 <5>T^CTAGA <6>6(6)")]
    [TestCase("NcoI", "?(4)", Description = "allenz.104 <5>C^CATGG <6>?(4)")]
    [TestCase("NdeI", "4(6)", Description = "allenz.104 <5>CA^TATG <6>4(6)")]
    [TestCase("NheI", "6(4)", Description = "allenz.104 <5>G^CTAGC <6>6(4)")]
    [TestCase("SpeI", "1(6)", Description = "allenz.104 <5>A^CTAGT <6>1(6)")]
    [TestCase("AvrII", "", Description = "allenz.104 <5>C^CTAGG <6>")]
    [TestCase("ClaI", "5(6)", Description = "allenz.104 <5>AT^CGAT <6>5(6)")]
    [TestCase("AgeI", "3(5)", Description = "allenz.104 <5>A^CCGGT <6>3(5)")]
    [TestCase("NotI", "2(4)", Description = "allenz.104 <5>GC^GGCCGC <6>2(4)")]
    [TestCase("SfiI", "?(4)", Description = "allenz.104 <5>GGCCNNNN^NGGCC <6>?(4)")]
    [TestCase("PacI", "", Description = "allenz.104 <5>TTAAT^TAA <6>")]
    [TestCase("AscI", "?(5)", Description = "allenz.104 <5>GG^CGCGCC <6>?(5)")]
    [TestCase("FseI", "?(5)", Description = "allenz.104 <5>GGCCGG^CC <6>?(5)")]
    [TestCase("SwaI", "?(6)", Description = "allenz.104 <5>ATTT^AAAT <6>?(6)")]
    [TestCase("MspI", "1(5)", Description = "allenz.104 <5>C^CGG <6>1(5)")]
    [TestCase("HpaII", "2(5)", Description = "allenz.104 <5>C^CGG <6>2(5)")]
    [TestCase("TaqI", "4(6)", Description = "allenz.104 <5>T^CGA <6>4(6)")]
    [TestCase("AluI", "3(5)", Description = "allenz.104 <5>AG^CT <6>3(5)")]
    [TestCase("RsaI", "4(4)", Description = "allenz.104 <5>GT^AC <6>4(4)")]
    [TestCase("HaeIII", "3(5)", Description = "allenz.104 <5>GG^CC <6>3(5)")]
    [TestCase("DpnI", "", Description = "allenz.104 <5>GA^TC <6>")]
    [TestCase("Sau3AI", "4(5)", Description = "allenz.104 <5>^GATC <6>4(5)")]
    [TestCase("MboI", "2(6)", Description = "allenz.104 <5>^GATC <6>2(6)")]
    [TestCase("EcoRV", "2(6)", Description = "allenz.104 <5>GAT^ATC <6>2(6)")]
    [TestCase("SmaI", "2(4)", Description = "allenz.104 <5>CCC^GGG <6>2(4)")]
    [TestCase("HincII", "5(6)", Description = "allenz.104 <5>GTY^RAC <6>5(6)")]
    [TestCase("ScaI", "5(4)", Description = "allenz.104 <5>AGT^ACT <6>5(4)")]
    [TestCase("StuI", "", Description = "allenz.104 <5>AGG^CCT <6>")]
    [TestCase("BglII", "5(4)", Description = "allenz.104 <5>A^GATCT <6>5(4)")]
    [TestCase("ApaI", "4(5)", Description = "allenz.104 <5>GGGCC^C <6>4(5)")]
    public void EnzymeDatabase_MethylationField_MatchesRebaseAllenz104(string name, string allenzField6)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(name);

        Assert.That(enzyme, Is.Not.Null, $"{name} should be in the built-in table");
        Assert.That(enzyme!.RebaseMethylationField, Is.EqualTo(allenzField6), $"{name} allenz.104 <6>");

        // Format header: X is "the base within the recognition sequence that is modified" and
        // (6) = N6-methyladenosine, (5) = 5-methylcytosine, (4) = N4-methylcytosine — so a known
        // top-strand position must lie inside the site and name an A for (6), a C for (4)/(5).
        foreach (var site in enzyme.CognateMethylationSites)
        {
            if (site.Position is not int pos || site.IsComplementaryStrand)
                continue;
            Assert.That(pos, Is.InRange(1, enzyme.RecognitionLength), $"{name} position inside the site");
            char expectedBase = site.Type == MethylatedBase.N6Methyladenine ? 'A' : 'C';
            Assert.That(enzyme.RecognitionSequence[pos - 1], Is.EqualTo(expectedBase),
                $"{name} {site.Type} at {pos}");
        }
    }

    /// <summary>The methylation audit above covers the whole table.</summary>
    [Test]
    public void EnzymeDatabase_MethylationAuditCoversEveryEntry()
    {
        var audited = typeof(RestrictionAnalyzer_FindSites_Tests)
            .GetMethod(nameof(EnzymeDatabase_MethylationField_MatchesRebaseAllenz104))!
            .GetCustomAttributes(typeof(TestCaseAttribute), false)
            .Cast<TestCaseAttribute>()
            .Select(a => (string)a.Arguments[0]!)
            .ToList();

        audited.Should().OnlyHaveUniqueItems();
        audited.Should().BeEquivalentTo(RestrictionAnalyzer.Enzymes.Keys,
            "every built-in enzyme carries a transcribed allenz.104 <6> field");
    }

    /// <summary>
    /// Typed parse of single-site entries, one per kind that occurs in the table (values from the
    /// test above): a known position, and <c>?</c> (type known, base not determined).
    /// </summary>
    [Test]
    public void CognateMethylationSites_BuiltIns_ParsedToTypedSites()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetEnzyme("EcoRI")!.CognateMethylationSites,
                Is.EqualTo(new[] { new CognateMethylationSite(3, false, MethylatedBase.N6Methyladenine) }));
            Assert.That(RestrictionAnalyzer.GetEnzyme("HpaII")!.CognateMethylationSites,
                Is.EqualTo(new[] { new CognateMethylationSite(2, false, MethylatedBase.C5Methylcytosine) }));
            Assert.That(RestrictionAnalyzer.GetEnzyme("BamHI")!.CognateMethylationSites,
                Is.EqualTo(new[] { new CognateMethylationSite(5, false, MethylatedBase.N4Methylcytosine) }));
            Assert.That(RestrictionAnalyzer.GetEnzyme("NcoI")!.CognateMethylationSites,
                Is.EqualTo(new[] { new CognateMethylationSite(null, false, MethylatedBase.N4Methylcytosine) }));
            Assert.That(RestrictionAnalyzer.GetEnzyme("DpnI")!.CognateMethylationSites, Is.Empty);
        });
    }

    /// <summary>
    /// The general REBASE notation, on field values that occur in allenz.104 (none in the
    /// built-in table): a negative position is the complementary strand, entries are
    /// comma-separated (with or without a space). A user-defined enzyme defaults to no site.
    /// </summary>
    [Test]
    public void CognateMethylationSites_RebaseNotation_ComplementaryStrandAndLists()
    {
        var both = new RestrictionEnzyme("X", "GATATC", 3, 3, "test", "3(6),-3(6)");
        var spaced = new RestrictionEnzyme("X", "GATATC", 3, 3, "test", "3(5), 2(5)");
        var mixed = new RestrictionEnzyme("X", "GATATC", 3, 3, "test", "2(4),?(6)");
        var plain = new RestrictionEnzyme("X", "GATATC", 3, 3, "test");

        Assert.Multiple(() =>
        {
            Assert.That(both.CognateMethylationSites, Is.EqualTo(new[]
            {
                new CognateMethylationSite(3, false, MethylatedBase.N6Methyladenine),
                new CognateMethylationSite(3, true, MethylatedBase.N6Methyladenine),
            }));
            Assert.That(spaced.CognateMethylationSites, Is.EqualTo(new[]
            {
                new CognateMethylationSite(3, false, MethylatedBase.C5Methylcytosine),
                new CognateMethylationSite(2, false, MethylatedBase.C5Methylcytosine),
            }));
            Assert.That(mixed.CognateMethylationSites, Is.EqualTo(new[]
            {
                new CognateMethylationSite(2, false, MethylatedBase.N4Methylcytosine),
                new CognateMethylationSite(null, false, MethylatedBase.N6Methyladenine),
            }));
            Assert.That(plain.RebaseMethylationField, Is.Empty);
            Assert.That(plain.CognateMethylationSites, Is.Empty);
        });
    }

    /// <summary>
    /// A value outside the documented notation is rejected at construction. <c>(99)</c> occurs
    /// 3× in allenz.104 but is not one of the three types the format header defines.
    /// </summary>
    [TestCase("1(99)")]
    [TestCase("3")]
    [TestCase("(6)")]
    [TestCase("0(6)")]
    [TestCase("-?(6)")]
    [TestCase("3(6),")]
    [TestCase("x(5)")]
    public void RestrictionEnzyme_MalformedMethylationField_Throws(string field)
    {
        Assert.Throws<FormatException>(() => _ = new RestrictionEnzyme("X", "GATATC", 3, 3, "test", field));
    }

    [Test]
    public void RestrictionEnzyme_NullMethylationField_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => _ = new RestrictionEnzyme("X", "GATATC", 3, 3, "test", null!));
    }

    #endregion

    #region Documented order, counts and compatibility ends (B08 H30/H31/H32/H33)

    /// <summary>
    /// H31: the <c>&lt;returns&gt;</c> doc used to promise that both strands' hits come out in
    /// ascending Position order. They do not: the reverse-strand scan runs 5'→3' along the reverse
    /// complement and <c>Position = n - i - len</c> falls as the scan index rises, so the
    /// reverse-strand hits are <b>descending</b> on a linear molecule and not monotone at all on a
    /// circular one, where the position is additionally reduced modulo the length.
    /// </summary>
    [Test]
    public void FindSites_ReverseStrandHits_ComeOutInDescendingPositionOrder()
    {
        var linear = RestrictionAnalyzer
            .FindSites(new DnaSequence("GAATTCAAGAATTCAAGAATTC"), "EcoRI")
            .Select(x => (x.Position, x.IsForwardStrand))
            .ToList();

        var circular = RestrictionAnalyzer
            .FindSites(new DnaSequence("GATCGATCGATC"), "Sau3AI", MoleculeTopology.Circular)
            .Select(x => x.Position)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(linear, Is.EqualTo(new[]
            {
                (0, true), (8, true), (16, true), (16, false), (8, false), (0, false),
            }), "forward ascending, then reverse descending");
            Assert.That(linear.Where(t => t.IsForwardStrand).Select(t => t.Position), Is.Ordered,
                "the forward half is ascending, which is all the doc may promise");
            Assert.That(linear.Where(t => !t.IsForwardStrand).Select(t => t.Position),
                Is.Ordered.Descending);

            Assert.That(circular, Is.EqualTo(new[] { 0, 4, 8, 8, 4, 0 }),
                "circular: the reverse half is not globally monotone against the forward half");
        });
    }

    /// <summary>
    /// H31: the same doc block claimed <c>Sau3AI</c> on <c>GATCGATCGATC</c> "reports the single site
    /// at position 4". One <i>window</i> survives the linear cut-feasibility filter, but it is
    /// reported as two records — the forward one carrying the top-strand cut and the reverse one the
    /// bottom-strand cut — by the two-records-per-palindromic-site rule the same doc states.
    /// Reference: Bio.Restriction.Sau3AI.search(Seq("GATCGATCGATC")) == [5] (1-based), the one
    /// top-strand cut.
    /// </summary>
    [Test]
    public void FindSites_Sau3AI_OneSurvivingWindow_IsTwoRecords()
    {
        var sites = RestrictionAnalyzer
            .FindSites(new DnaSequence("GATCGATCGATC"), "Sau3AI")
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(2), "one window, two strand records");
            Assert.That(sites.Select(x => x.Position), Is.All.EqualTo(4));
            Assert.That(sites.Single(x => x.IsForwardStrand).CutPosition, Is.EqualTo(4),
                "Biopython Sau3AI.search(Seq(\"GATCGATCGATC\")) == [5] (1-based)");
            Assert.That(sites.Single(x => !x.IsForwardStrand).CutPosition, Is.EqualTo(8));
        });
    }

    /// <summary>
    /// H32: <c>TotalSites</c> and <c>UniqueCutters</c> count cleavages, not forward-strand records.
    /// One overload's doc used to say "forward-strand records only", which reads wrong for every
    /// Type IIS enzyme — i.e. every Golden Gate enzyme. A molecule carrying GGTCTC on one strand and
    /// on the other has two BsaI cleavages, so BsaI is not a unique cutter.
    /// </summary>
    [Test]
    public void CreateMap_TypeIISEnzyme_CountsTheReverseStrandSiteAndIsNotAUniqueCutter()
    {
        const string molecule = "AAAAAAAAAAGGTCTCAAAAAAAAAAGAGACCAAAAAAAAAA";
        var map = RestrictionAnalyzer.CreateMap(new DnaSequence(molecule), MoleculeTopology.Linear, "BsaI");
        var sites = RestrictionAnalyzer.FindSites(new DnaSequence(molecule), MoleculeTopology.Linear, "BsaI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites.Count(x => x.IsForwardStrand), Is.EqualTo(1), "one forward-strand record");
            Assert.That(sites.Count(x => !x.IsForwardStrand), Is.EqualTo(1), "and one reverse-strand record");
            Assert.That(map.TotalSites, Is.EqualTo(2), "both are cleavages; 'forward records only' would say 1");
            Assert.That(map.UniqueCutters, Is.Empty, "so BsaI is not a single cutter here");
            Assert.That(map.SitesByEnzyme["BsaI"], Is.EqualTo(new[] { 10, 26 }));
        });
    }

    /// <summary>
    /// H33 / H30: <c>FindCompatibleEnzymes(name)</c>'s end is the <i>queried</i> enzyme's own
    /// overhang window, not a string both enzymes leave. 12 REBASE pairs anneal only through the
    /// reverse complement, and all 12 are <c>GN</c> against <c>CN</c>: BsmI leaves <c>CN</c>, BsrI
    /// leaves <c>GN</c>, and each query reports its own. The MCP description used to say the value
    /// "echoes the shared end".
    /// </summary>
    [Test]
    public void FindCompatibleEnzymes_ReverseComplementOnlyPair_ReportsTheQueriedEnzymesOwnEnd()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetOverhangSequence("BsmI"), Is.EqualTo("CN"));
            Assert.That(RestrictionAnalyzer.GetOverhangSequence("BsrI"), Is.EqualTo("GN"));
            Assert.That(RestrictionAnalyzer.AreCompatible("BsmI", "BsrI"), Is.True,
                "CN anneals to GN through the reverse complement");

            Assert.That(
                RestrictionAnalyzer.FindCompatibleEnzymes("BsmI", EnzymeSet.Rebase)
                    .Single(t => t.Enzyme == "BsrI").CompatibleEnd,
                Is.EqualTo("CN"), "the queried enzyme's own end");
            Assert.That(
                RestrictionAnalyzer.FindCompatibleEnzymes("BsrI", EnzymeSet.Rebase)
                    .Single(t => t.Enzyme == "BsmI").CompatibleEnd,
                Is.EqualTo("GN"), "and the other way round — the value is order-dependent");
        });
    }

    /// <summary>
    /// H33: Biopython's <c>compatible_end()</c> has three implementations, not one 5'-only filter,
    /// so it answers a non-empty list for blunt and 3'-overhang enzymes too — the doc used to claim
    /// the opposite. These counts are this library's side of the same question over the REBASE set;
    /// the reference numbers measured with Biopython 1.88 are
    /// <c>EcoRV.compatible_end()</c> → 137 (every blunt cutter), <c>PstI</c> → 13, <c>SacI</c> → 9.
    /// The blunt row is directly comparable because every blunt end ligates to every other.
    /// </summary>
    [Test]
    public void FindCompatibleEnzymes_BluntAndThreePrimeEnzymes_AreNotEmpty()
    {
        var bluntCount = RestrictionAnalyzer.GetBluntCutters(EnzymeSet.Rebase).Count();
        var ecoRV = RestrictionAnalyzer.FindCompatibleEnzymes("EcoRV", EnzymeSet.Rebase).ToList();
        var pstI = RestrictionAnalyzer.FindCompatibleEnzymes("PstI", EnzymeSet.Rebase).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(bluntCount, Is.EqualTo(137),
                "Biopython EcoRV.compatible_end() returns 137 enzymes, i.e. every blunt cutter");
            Assert.That(ecoRV, Has.Count.EqualTo(137), "EcoRV is blunt, so every blunt cutter qualifies");
            Assert.That(ecoRV.Select(t => t.CompatibleEnd), Is.All.EqualTo("blunt"));
            Assert.That(pstI, Is.Not.Empty, "PstI leaves a 3' overhang and still has partners");
            Assert.That(pstI.Select(t => t.Enzyme), Does.Contain("NsiI"),
                "NsiI is in Biopython's PstI.compatible_end() list");
        });
    }

    /// <summary>
    /// H33: the queried enzyme appears in its own result only when it is a member of the chosen set.
    /// BsaI is a REBASE enzyme, not one of the 39 common ones, so the common-set query cannot list
    /// it although <c>GetEnzyme</c> resolves it — the doc used to say "the enzyme itself included"
    /// without that qualification.
    /// </summary>
    [Test]
    public void FindCompatibleEnzymes_QueryOutsideTheChosenSet_DoesNotListItself()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.GetEnzyme("BsaI"), Is.Not.Null);
            Assert.That(RestrictionAnalyzer.Enzymes.ContainsKey("BsaI"), Is.False, "not a common enzyme");
            Assert.That(
                RestrictionAnalyzer.FindCompatibleEnzymes("BsaI", EnzymeSet.Common).Select(t => t.Enzyme),
                Does.Not.Contain("BsaI"));
            Assert.That(
                RestrictionAnalyzer.FindCompatibleEnzymes("BsaI", EnzymeSet.Rebase).Select(t => t.Enzyme),
                Does.Contain("BsaI"), "it is a member of the REBASE set, so there it does appear");
            Assert.That(
                RestrictionAnalyzer.FindCompatibleEnzymes("EcoRI", EnzymeSet.Common).Select(t => t.Enzyme),
                Does.Contain("EcoRI"), "a common-set enzyme does appear in the common-set query");
        });
    }

    #endregion

    #region Raw-sequence-string entry points (B08 H27)

    /// <summary>
    /// H27: <c>DnaSequence</c> accepts only A/C/G/T, so an N-containing contig or a soft-masked
    /// assembly could not be handed to the digest/map family at all. The raw-string overloads are the
    /// reachable form of the IUPAC-tolerant scan the raw-string <c>FindSites</c> and
    /// <c>FindPatternOccurrences</c> already documented, and each forwards to the same private core —
    /// so on a pure-ACGT molecule they are byte-identical to the <c>DnaSequence</c> overloads.
    /// </summary>
    [TestCase(MoleculeTopology.Linear)]
    [TestCase(MoleculeTopology.Circular)]
    public void RawStringOverloads_OnPureAcgt_EqualTheDnaSequenceOverloads(MoleculeTopology topology)
    {
        const string molecule = "AAAGAATTCAAAGGATCCAAAGAATTCAAA";
        var dna = new DnaSequence(molecule);
        string[] names = ["EcoRI", "BamHI"];

        Assert.Multiple(() =>
        {
            Assert.That(RestrictionAnalyzer.FindSites(molecule, topology, null, names),
                Is.EqualTo(RestrictionAnalyzer.FindSites(dna, topology, null, names)));
            Assert.That(RestrictionAnalyzer.FindAllSites(molecule, topology, EnzymeSet.Common, null),
                Is.EqualTo(RestrictionAnalyzer.FindAllSites(dna, topology, EnzymeSet.Common, null)));
            Assert.That(RestrictionAnalyzer.Digest(molecule, topology, null, names),
                Is.EqualTo(RestrictionAnalyzer.Digest(dna, topology, null, names)));
            AssertSameSummary(
                RestrictionAnalyzer.GetDigestSummary(molecule, topology, null, names),
                RestrictionAnalyzer.GetDigestSummary(dna, topology, null, names));
            Assert.That(RestrictionAnalyzer.PartialDigest(molecule, topology, null, null, names),
                Is.EqualTo(RestrictionAnalyzer.PartialDigest(dna, topology, null, null, names)));
            AssertSameMap(
                RestrictionAnalyzer.CreateMap(molecule, topology, null, names),
                RestrictionAnalyzer.CreateMap(dna, topology, null, names));

            // The enzyme-record and methylation forms too.
            var enzymes = names.Select(n => RestrictionAnalyzer.GetEnzyme(n)!).ToList();
            Assert.That(RestrictionAnalyzer.FindSites(molecule, topology, enzymes, MethylaseSets.Dam),
                Is.EqualTo(RestrictionAnalyzer.FindSites(dna, topology, enzymes, MethylaseSets.Dam)));
            Assert.That(RestrictionAnalyzer.Digest(molecule, topology, enzymes, MethylaseSets.Dam),
                Is.EqualTo(RestrictionAnalyzer.Digest(dna, topology, enzymes, MethylaseSets.Dam)));
            AssertSameMap(
                RestrictionAnalyzer.CreateMap(molecule, topology, enzymes, MethylaseSets.Dam),
                RestrictionAnalyzer.CreateMap(dna, topology, enzymes, MethylaseSets.Dam));
            AssertSameSummary(
                RestrictionAnalyzer.GetDigestSummary(molecule, topology, enzymes, MethylaseSets.Dam),
                RestrictionAnalyzer.GetDigestSummary(dna, topology, enzymes, MethylaseSets.Dam));
            Assert.That(RestrictionAnalyzer.PartialDigest(molecule, topology, null, enzymes, MethylaseSets.Dam),
                Is.EqualTo(RestrictionAnalyzer.PartialDigest(dna, topology, null, enzymes, MethylaseSets.Dam)));
        });
    }

    /// <summary>
    /// H27: and the point of them — a sequence <c>DnaSequence</c> rejects. The masked bases do not
    /// become wildcards: an ambiguity code in the sequence blocks a concrete site code rather than
    /// widening it, which is Biopython's rule, so the <c>nnn</c> runs add no site.
    /// </summary>
    [Test]
    public void RawStringOverloads_AcceptAnIupacSequenceDnaSequenceRejects()
    {
        const string masked = "nnnGAATTCnnnGGATCCNNN";

        Assert.Throws<ArgumentException>(() => _ = new DnaSequence(masked),
            "the defect: DnaSequence accepts only A/C/G/T");

        Assert.Multiple(() =>
        {
            Assert.That(
                RestrictionAnalyzer.FindSites(masked, MoleculeTopology.Linear, null, "EcoRI", "BamHI")
                    .Where(s => s.IsForwardStrand).Select(s => s.Position),
                Is.EquivalentTo(new[] { 3, 12 }));
            Assert.That(
                RestrictionAnalyzer.Digest(masked, MoleculeTopology.Linear, null, "EcoRI").Select(f => f.Length),
                Is.EqualTo(new[] { 4, 17 }), "G^AATTC at 3, so the 0-based top-strand bond is 4");
            Assert.That(
                RestrictionAnalyzer.CreateMap(masked, MoleculeTopology.Linear, null, "EcoRI", "BamHI").TotalSites,
                Is.EqualTo(2));
            // No extra site from the masked runs: a sequence N never matches a concrete site code.
            Assert.That(RestrictionAnalyzer.FindSites(masked, MoleculeTopology.Linear, null, "AluI"), Is.Empty,
                "AluI AGCT has no window here, and the N runs do not create one");
            Assert.That(RestrictionAnalyzer.FindPatternOccurrences("GAANTC", "GAATTC"), Is.Empty);
        });
    }

    /// <summary>
    /// H27: a null or empty raw sequence yields the natural empty result, as it does on the
    /// pre-existing raw-string <c>FindSites(string, string)</c> — while the enzyme arguments are still
    /// validated eagerly, before the sequence is looked at (B08 H24/F104).
    /// </summary>
    [Test]
    public void RawStringOverloads_NullOrEmptySequence_YieldTheEmptyResult()
    {
        foreach (string? sequence in new string?[] { null, "" })
        {
            Assert.Multiple(() =>
            {
                Assert.That(RestrictionAnalyzer.FindSites(sequence, MoleculeTopology.Linear, null, "EcoRI"), Is.Empty);
                Assert.That(RestrictionAnalyzer.FindAllSites(sequence, MoleculeTopology.Linear, EnzymeSet.Common, null), Is.Empty);
                Assert.That(RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Linear, null, "EcoRI").Select(f => f.Length),
                    Is.EqualTo(new[] { 0 }), "one empty fragment, as the DnaSequence overload gives for an empty molecule");
                Assert.That(RestrictionAnalyzer.PartialDigest(sequence, MoleculeTopology.Linear, null, null, "EcoRI"),
                    Is.Not.Null);
                Assert.That(RestrictionAnalyzer.CreateMap(sequence, MoleculeTopology.Linear, null, "EcoRI").SequenceLength,
                    Is.Zero);

                Assert.Throws<ArgumentException>(
                    () => _ = RestrictionAnalyzer.FindSites(sequence, MoleculeTopology.Linear, null, "NoSuchI"));
                Assert.Throws<ArgumentException>(
                    () => _ = RestrictionAnalyzer.Digest(sequence, MoleculeTopology.Linear, null, System.Array.Empty<string>()));
            });
        }
    }

    /// <summary>
    /// <see cref="DigestSummary"/> is a record holding a <c>List&lt;T&gt;</c>, so its generated
    /// equality is reference equality on that list; compare the values.
    /// </summary>
    private static void AssertSameSummary(DigestSummary actual, DigestSummary expected)
    {
        Assert.That(
            (actual.TotalFragments, actual.FragmentSizes, actual.LargestFragment, actual.SmallestFragment,
                actual.AverageFragmentSize, actual.EnzymesUsed, actual.Topology),
            Is.EqualTo((expected.TotalFragments, expected.FragmentSizes, expected.LargestFragment,
                expected.SmallestFragment, expected.AverageFragmentSize, expected.EnzymesUsed, expected.Topology)));
    }

    /// <summary>Same for <see cref="RestrictionMap"/>.</summary>
    private static void AssertSameMap(RestrictionMap actual, RestrictionMap expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(actual.SequenceLength, Is.EqualTo(expected.SequenceLength));
            Assert.That(actual.Sites, Is.EqualTo(expected.Sites));
            Assert.That(actual.SitesByEnzyme, Is.EqualTo(expected.SitesByEnzyme));
            Assert.That(actual.TotalSites, Is.EqualTo(expected.TotalSites));
            Assert.That(actual.UniqueCutters, Is.EqualTo(expected.UniqueCutters));
            Assert.That(actual.NonCutters, Is.EqualTo(expected.NonCutters));
        });
    }

    #endregion
}
