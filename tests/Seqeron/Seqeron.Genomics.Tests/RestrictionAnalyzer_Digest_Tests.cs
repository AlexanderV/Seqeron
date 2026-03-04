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
        // Arrange: EcoRI (G↓AATTC) at position 3, cuts at position 4
        var sequence = new DnaSequence("AAAGAATTCAAA"); // 12 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariants #1 and #2 with exact values
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(2), "Single cut → 2 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(12),
                "Fragment sum must equal original sequence length (Addgene)");
            Assert.That(fragments[0].StartPosition, Is.EqualTo(0));
            Assert.That(fragments[0].Length, Is.EqualTo(4), "First fragment: AAAG");
            Assert.That(fragments[1].StartPosition, Is.EqualTo(4));
            Assert.That(fragments[1].Length, Is.EqualTo(8), "Second fragment: AATTCAAA");
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
        // Arrange: Two EcoRI sites at positions 0 and 9, cuts at 1 and 10
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA"); // 18 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Invariants #1 and #2 with exact values
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(3), "Two cuts → 3 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(18));
            Assert.That(fragments[0].Length, Is.EqualTo(1), "Fragment 1: G");
            Assert.That(fragments[1].Length, Is.EqualTo(9), "Fragment 2: AATTCAAAG");
            Assert.That(fragments[2].Length, Is.EqualTo(8), "Fragment 3: AATTCAAA");
        });
    }

    /// <summary>
    /// Evidence: Comprehensive fragment property verification with exact values.
    /// Invariants #3 (sequential numbering), #4 (monotonic positions),
    /// #5 (boundary enzymes), #6 (positive length).
    /// </summary>
    [Test]
    public void Digest_FragmentsHaveCorrectProperties()
    {
        // Arrange: 3 EcoRI sites at positions 0, 9, 18 → cuts at 1, 10, 19 → 4 fragments
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAAGAATTCAAA"); // 27 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: All fragment properties with exact values
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(4), "3 cuts → 4 fragments");

            // Fragment 1: [0, 1) = "G"
            Assert.That(fragments[0].FragmentNumber, Is.EqualTo(1));
            Assert.That(fragments[0].StartPosition, Is.EqualTo(0));
            Assert.That(fragments[0].Length, Is.EqualTo(1));
            Assert.That(fragments[0].LeftEnzyme, Is.Null, "First fragment: LeftEnzyme null (Invariant #5)");
            Assert.That(fragments[0].RightEnzyme, Is.EqualTo("EcoRI"));

            // Fragment 2: [1, 10) = "AATTCAAAG"
            Assert.That(fragments[1].FragmentNumber, Is.EqualTo(2));
            Assert.That(fragments[1].StartPosition, Is.EqualTo(1));
            Assert.That(fragments[1].Length, Is.EqualTo(9));
            Assert.That(fragments[1].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[1].RightEnzyme, Is.EqualTo("EcoRI"));

            // Fragment 3: [10, 19) = "AATTCAAAG"
            Assert.That(fragments[2].FragmentNumber, Is.EqualTo(3));
            Assert.That(fragments[2].StartPosition, Is.EqualTo(10));
            Assert.That(fragments[2].Length, Is.EqualTo(9));
            Assert.That(fragments[2].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[2].RightEnzyme, Is.EqualTo("EcoRI"));

            // Fragment 4: [19, 27) = "AATTCAAA"
            Assert.That(fragments[3].FragmentNumber, Is.EqualTo(4));
            Assert.That(fragments[3].StartPosition, Is.EqualTo(19));
            Assert.That(fragments[3].Length, Is.EqualTo(8));
            Assert.That(fragments[3].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[3].RightEnzyme, Is.Null, "Last fragment: RightEnzyme null (Invariant #5)");

            // Invariant #4: Monotonic start positions (absorbed from dedicated test)
            for (int i = 1; i < fragments.Count; i++)
                Assert.That(fragments[i].StartPosition, Is.GreaterThan(fragments[i - 1].StartPosition),
                    $"Fragment {i + 1} start > Fragment {i} start");

            // Invariant #6: All positive length
            Assert.That(fragments.All(f => f.Length > 0), "All fragments have positive length");
        });
    }

    /// <summary>
    /// Evidence: Verify fragment sequence content matches exact known substrings.
    /// </summary>
    [Test]
    public void Digest_FragmentSequenceContent_MatchesExpectedSubstring()
    {
        // Arrange: EcoRI (G↓AATTC) at position 3, cuts at position 4
        var sequence = new DnaSequence("AAAGAATTCAAA"); // 12 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI").ToList();

        // Assert: Exact known sequences
        Assert.Multiple(() =>
        {
            Assert.That(fragments[0].Sequence, Is.EqualTo("AAAG"), "Fragment 1: positions [0, 4)");
            Assert.That(fragments[1].Sequence, Is.EqualTo("AATTCAAA"), "Fragment 2: positions [4, 12)");
        });
    }

    /// <summary>
    /// Evidence: Multiple enzymes produce combined cut sites with correct enzyme attribution.
    /// </summary>
    [Test]
    public void Digest_MultipleEnzymes_CutsWithBoth()
    {
        // Arrange: EcoRI at position 3 (cuts at 4), BamHI at position 12 (cuts at 13)
        var sequence = new DnaSequence("AAAGAATTCAAAGGATCCAAA"); // 21 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "EcoRI", "BamHI").ToList();

        // Assert: Exact fragment properties with enzyme attribution
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(3), "Two cuts → 3 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(21));
            Assert.That(fragments[0].Length, Is.EqualTo(4), "Fragment 1: AAAG");
            Assert.That(fragments[0].RightEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[1].Length, Is.EqualTo(9), "Fragment 2: AATTCAAAG");
            Assert.That(fragments[1].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[1].RightEnzyme, Is.EqualTo("BamHI"));
            Assert.That(fragments[2].Length, Is.EqualTo(8), "Fragment 3: GATCCAAA");
            Assert.That(fragments[2].LeftEnzyme, Is.EqualTo("BamHI"));
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

    /// <summary>
    /// Evidence: Edge case — adjacent cut positions produce a 1-bp fragment.
    /// TaqI (T↓CGA) at position 3 cuts at 4; MboI (↓GATC) at position 5 cuts at 5.
    /// All invariants must hold: positive length, sum = original.
    /// </summary>
    [Test]
    public void Digest_AdjacentCutSites_ProducesSmallFragment()
    {
        // Arrange: TaqI cuts at 4, MboI cuts at 5 → adjacent cuts, 1-bp middle fragment
        var sequence = new DnaSequence("AAATCGATCAAA"); // 12 bp

        // Act
        var fragments = RestrictionAnalyzer.Digest(sequence, "TaqI", "MboI").ToList();

        // Assert: Adjacent cuts handled correctly with exact values
        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(3), "Two adjacent cuts → 3 fragments");
            Assert.That(fragments.Sum(f => f.Length), Is.EqualTo(12));
            Assert.That(fragments[0].Sequence, Is.EqualTo("AAAT"), "Fragment 1: [0, 4)");
            Assert.That(fragments[1].Sequence, Is.EqualTo("C"), "Fragment 2: [4, 5) — 1 bp");
            Assert.That(fragments[2].Sequence, Is.EqualTo("GATCAAA"), "Fragment 3: [5, 12)");
            Assert.That(fragments.All(f => f.Length > 0), "All fragments have positive length");
        });
    }

    #endregion

    #region GetDigestSummary Tests

    /// <summary>
    /// Evidence: Summary invariants verified with exact values.
    /// Invariants #8 (sorted descending), #9 (bounds), #10 (enzyme list).
    /// EcoRI sites at positions 0 and 9 → fragments of length 1, 9, 8.
    /// </summary>
    [Test]
    public void GetDigestSummary_ReturnsCorrectSummaryWithInvariants()
    {
        // Arrange: 2 EcoRI sites → fragments [1, 9, 8], sorted desc: [9, 8, 1]
        var sequence = new DnaSequence("GAATTCAAAGAATTCAAA"); // 18 bp

        // Act
        var summary = RestrictionAnalyzer.GetDigestSummary(sequence, "EcoRI");

        // Assert: All invariants with exact values
        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalFragments, Is.EqualTo(3));

            // Invariant #8: Fragment sizes sorted descending — exact values
            Assert.That(summary.FragmentSizes, Is.EqualTo(new[] { 9, 8, 1 }),
                "Fragment sizes sorted descending");

            // Invariant #9: Exact bounds
            Assert.That(summary.LargestFragment, Is.EqualTo(9));
            Assert.That(summary.SmallestFragment, Is.EqualTo(1));
            Assert.That(summary.AverageFragmentSize, Is.EqualTo(6.0));

            // Invariant #10: Exact enzyme list
            Assert.That(summary.EnzymesUsed, Is.EqualTo(new[] { "EcoRI" }));
        });
    }

    #endregion

    #region Restriction Map Tests

    /// <summary>
    /// Evidence: CreateMap returns complete map data with exact values.
    /// Invariants #12, #13: Unique cutters, site count (forward-strand only).
    /// Two EcoRI sites → not a unique cutter.
    /// </summary>
    [Test]
    public void CreateMap_ReturnsCorrectMapWithAllFields()
    {
        // Arrange: Two EcoRI sites at positions 0 and 9 (palindromic → both strands per site)
        var sequence = new DnaSequence("GAATTCAAAGAATTC"); // 15 bp

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI");

        // Assert: Exact values for all fields
        Assert.Multiple(() =>
        {
            Assert.That(map.SequenceLength, Is.EqualTo(15));
            Assert.That(map.TotalSites, Is.EqualTo(2), "Forward-strand sites only (Invariant #13)");
            Assert.That(map.SitesByEnzyme["EcoRI"], Is.EqualTo(new[] { 0, 0, 9, 9 }),
                "Both strands: 2 forward + 2 reverse at positions 0 and 9");
            Assert.That(map.UniqueCutters, Is.Empty, "EcoRI cuts twice — not a unique cutter");
            Assert.That(map.NonCutters, Is.Empty);
        });
    }

    /// <summary>
    /// Evidence: Implementation - Unique cutters are enzymes with exactly one forward-strand site.
    /// Invariant #12: UniqueCutters list contains enzymes with one forward-strand site.
    /// For palindromic enzymes, both strands match at the same position — only forward counts.
    /// </summary>
    [Test]
    public void CreateMap_IdentifiesUniqueCutters()
    {
        // Arrange: Sequence where each enzyme cuts only once (forward strand)
        var sequence = new DnaSequence("AAAAAAAAAGAATTCAAAAAAAAAAAAAAAGGATCCAAAAAAAAAA");

        // Act
        var map = RestrictionAnalyzer.CreateMap(sequence, "EcoRI", "BamHI");

        // Assert: Invariant #12 — unique cutters correctly identified for palindromic sites
        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(2), "Two forward-strand sites");
            Assert.That(map.UniqueCutters, Contains.Item("EcoRI"), "EcoRI is a unique cutter");
            Assert.That(map.UniqueCutters, Contains.Item("BamHI"), "BamHI is a unique cutter");
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

    /// <summary>
    /// Evidence: API contract — When no enzymes specified, CreateMap searches all known enzymes.
    /// Sequence contains EcoRI (pos 9), BamHI (pos 30), and GATC-cutters (Sau3AI, MboI, DpnI at pos 31).
    /// </summary>
    [Test]
    public void CreateMap_NoEnzymesSpecified_SearchesAll()
    {
        // Arrange: Known EcoRI and BamHI sites (BamHI's GGATCC contains GATC for other cutters)
        var sequence = new DnaSequence("AAAAAAAAAGAATTCAAAAAAAAAAAAAAAGGATCCAAAAAAAAAA"); // 45 bp

        // Act: No enzymes specified — should search all
        var map = RestrictionAnalyzer.CreateMap(sequence);

        // Assert: Finds multiple enzyme families
        Assert.Multiple(() =>
        {
            Assert.That(map.SitesByEnzyme.ContainsKey("EcoRI"), Is.True, "Should find EcoRI site");
            Assert.That(map.SitesByEnzyme.ContainsKey("BamHI"), Is.True, "Should find BamHI site");
            Assert.That(map.SitesByEnzyme.ContainsKey("Sau3AI"), Is.True,
                "Should find Sau3AI (GATC within GGATCC)");
            Assert.That(map.TotalSites, Is.EqualTo(5),
                "5 forward-strand sites: EcoRI, BamHI, Sau3AI, MboI, DpnI");
            Assert.That(map.NonCutters, Is.Empty, "No enzyme names specified → no non-cutters");
        });
    }

    #endregion

    #region Compatibility Tests

    /// <summary>
    /// Evidence: Wikipedia (Sticky and blunt ends) — "blunt ends are always compatible with each other."
    /// Invariant #14: All blunt-end enzymes are compatible.
    /// </summary>
    [Test]
    public void AreCompatible_BluntEnzymes_AreCompatible()
    {
        // EcoRV and SmaI are both blunt cutters (Wikipedia examples table: * = blunt ends)
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRV", "SmaI");

        Assert.That(compatible, Is.True, "Blunt enzymes should be compatible (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia (Restriction enzyme) — BamHI and BglII both produce 5' GATC overhangs.
    /// BamHI: 5'---G↓GATCC---3' / 3'---CCTAG↓G---5' → 5' overhang GATC
    /// BglII: 5'---A↓GATCT---3' / 3'---TCTAG↓A---5' → 5' overhang GATC
    /// Invariant #15: Same overhang type + same sequence = compatible.
    /// </summary>
    [Test]
    public void AreCompatible_SameOverhang_AreCompatible_BamHI_BglII()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("BamHI", "BglII");

        Assert.That(compatible, Is.True, "BamHI/BglII both produce 5' GATC overhang (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — SalI and XhoI both produce 5' TCGA overhangs.
    /// SalI:  5'---G↓TCGAC---3' / 3'---CAGCT↓G---5' → 5' overhang TCGA
    /// XhoI:  5'---C↓TCGAG---3' / 3'---GAGCT↓C---5' → 5' overhang TCGA
    /// </summary>
    [Test]
    public void AreCompatible_SameOverhang_AreCompatible_SalI_XhoI()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("SalI", "XhoI");

        Assert.That(compatible, Is.True, "SalI/XhoI both produce 5' TCGA overhang (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — MboI, Sau3AI, BamHI, BglII all produce 5' GATC overhangs.
    /// MboI/Sau3AI: 5'---↓GATC---3' / 3'---CTAG↓---5' → 5' overhang GATC
    /// </summary>
    [Test]
    [TestCase("MboI", "BamHI")]
    [TestCase("Sau3AI", "BglII")]
    [TestCase("MboI", "Sau3AI")]
    public void AreCompatible_GATC_OverhangFamily_AllCompatible(string enzyme1, string enzyme2)
    {
        bool compatible = RestrictionAnalyzer.AreCompatible(enzyme1, enzyme2);

        Assert.That(compatible, Is.True,
            $"{enzyme1}/{enzyme2} both produce 5' GATC overhang (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — NheI, XbaI, SpeI, AvrII all produce 5' CTAG overhangs.
    /// NheI:  5'---G↓CTAGC---3' → 5' overhang CTAG
    /// XbaI:  5'---T↓CTAGA---3' → 5' overhang CTAG
    /// SpeI:  5'---A↓CTAGT---3' → 5' overhang CTAG
    /// AvrII: 5'---C↓CTAGG---3' → 5' overhang CTAG
    /// </summary>
    [Test]
    [TestCase("NheI", "XbaI")]
    [TestCase("NheI", "SpeI")]
    [TestCase("XbaI", "AvrII")]
    public void AreCompatible_CTAG_OverhangFamily_AllCompatible(string enzyme1, string enzyme2)
    {
        bool compatible = RestrictionAnalyzer.AreCompatible(enzyme1, enzyme2);

        Assert.That(compatible, Is.True,
            $"{enzyme1}/{enzyme2} both produce 5' CTAG overhang (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia (Sticky and blunt ends) — "overhangs have to be complementary
    /// in order for the ligase to work." A 5' overhang CANNOT ligate with a 3' overhang,
    /// even when the overhang sequence string is the same palindrome, because both
    /// single-stranded extensions would be on the same strand and cannot base-pair.
    ///
    /// HindIII: 5'---A↓AGCTT---3' / 3'---TTCGA↓A---5' → 5' overhang AGCT
    /// SacI:    5'---GAGCT↓C---3' / 3'---C↓TCGAG---5' → 3' overhang AGCT
    /// </summary>
    [Test]
    public void AreCompatible_CrossTypeOverhang_NotCompatible_HindIII_SacI()
    {
        // HindIII produces 5' overhang AGCT; SacI produces 3' overhang AGCT
        // Same overhang string but different types → NOT compatible
        bool compatible = RestrictionAnalyzer.AreCompatible("HindIII", "SacI");

        Assert.That(compatible, Is.False,
            "HindIII (5' AGCT) and SacI (3' AGCT) are NOT compatible — different overhang types (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — SphI produces 3' overhang CATG, NcoI produces 5' overhang CATG.
    /// SphI:  5'---GCATG↓C---3' / 3'---C↓GTACG---5' → 3' overhang CATG
    /// NcoI:  5'---C↓CATGG---3' / 3'---GGTAC↓C---5' → 5' overhang CATG
    /// Same sequence CATG but different types → NOT compatible.
    /// </summary>
    [Test]
    public void AreCompatible_CrossTypeOverhang_NotCompatible_SphI_NcoI()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("SphI", "NcoI");

        Assert.That(compatible, Is.False,
            "SphI (3' CATG) and NcoI (5' CATG) are NOT compatible — different overhang types (Wikipedia)");
    }

    /// <summary>
    /// Evidence: Wikipedia — Different overhang sequences cannot ligate.
    /// EcoRI (5' AATT overhang) and PstI (3' TGCA overhang): different type AND different sequence.
    /// </summary>
    [Test]
    public void AreCompatible_DifferentOverhangs_NotCompatible()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRI", "PstI");

        Assert.That(compatible, Is.False, "Different overhang enzymes should not be compatible");
    }

    /// <summary>
    /// Evidence: API contract — Unknown enzyme returns false (no exception).
    /// </summary>
    [Test]
    public void AreCompatible_UnknownEnzyme_ReturnsFalse()
    {
        bool compatible = RestrictionAnalyzer.AreCompatible("EcoRI", "UnknownEnzyme");

        Assert.That(compatible, Is.False, "Unknown enzyme should return false, not throw");
    }

    /// <summary>
    /// Evidence: Mathematical property — Compatibility is symmetric.
    /// Invariant #16: AreCompatible(A, B) == AreCompatible(B, A).
    /// </summary>
    [Test]
    [TestCase("BamHI", "BglII")]
    [TestCase("EcoRV", "SmaI")]
    [TestCase("EcoRI", "PstI")]
    [TestCase("HindIII", "SacI")]
    [TestCase("SphI", "NcoI")]
    public void AreCompatible_IsSymmetric(string enzyme1, string enzyme2)
    {
        bool forward = RestrictionAnalyzer.AreCompatible(enzyme1, enzyme2);
        bool reverse = RestrictionAnalyzer.AreCompatible(enzyme2, enzyme1);

        Assert.That(forward, Is.EqualTo(reverse),
            $"AreCompatible({enzyme1}, {enzyme2}) should equal AreCompatible({enzyme2}, {enzyme1})");
    }

    /// <summary>
    /// Evidence: Wikipedia (Restriction enzyme) — Verified enzyme cut positions from Examples table.
    /// EcoRI:  G↓AATTC (forward=1, reverse=5), 5' overhang AATT
    /// BamHI:  G↓GATCC (forward=1, reverse=5), 5' overhang GATC
    /// PstI:   CTGCA↓G (forward=5, reverse=1), 3' overhang TGCA
    /// SmaI:   CCC↓GGG (forward=3, reverse=3), Blunt
    /// NotI:   GC↓GGCCGC (forward=2, reverse=6), 5' overhang GGCC
    /// </summary>
    [Test]
    [TestCase("EcoRI", "GAATTC", 1, 5, false, Description = "Wikipedia: G↓AATTC, 5' overhang")]
    [TestCase("BamHI", "GGATCC", 1, 5, false, Description = "Wikipedia: G↓GATCC, 5' overhang")]
    [TestCase("HindIII", "AAGCTT", 1, 5, false, Description = "Wikipedia: A↓AGCTT, 5' overhang")]
    [TestCase("PstI", "CTGCAG", 5, 1, false, Description = "Wikipedia: CTGCA↓G, 3' overhang")]
    [TestCase("SmaI", "CCCGGG", 3, 3, true, Description = "Wikipedia: CCC↓GGG, Blunt")]
    [TestCase("EcoRV", "GATATC", 3, 3, true, Description = "Wikipedia: GAT↓ATC, Blunt")]
    [TestCase("AluI", "AGCT", 2, 2, true, Description = "Wikipedia: AG↓CT, Blunt")]
    [TestCase("HaeIII", "GGCC", 2, 2, true, Description = "Wikipedia: GG↓CC, Blunt")]
    [TestCase("NotI", "GCGGCCGC", 2, 6, false, Description = "Wikipedia: GC↓GGCCGC, 5' overhang")]
    [TestCase("TaqI", "TCGA", 1, 3, false, Description = "Wikipedia: T↓CGA, 5' overhang")]
    [TestCase("Sau3AI", "GATC", 0, 4, false, Description = "Wikipedia: ↓GATC, 5' overhang")]
    [TestCase("KpnI", "GGTACC", 5, 1, false, Description = "Wikipedia: GGTAC↓C, 3' overhang")]
    [TestCase("SacI", "GAGCTC", 5, 1, false, Description = "Wikipedia: GAGCT↓C, 3' overhang")]
    [TestCase("SphI", "GCATGC", 5, 1, false, Description = "Wikipedia: GCATG↓C, 3' overhang")]
    [TestCase("XbaI", "TCTAGA", 1, 5, false, Description = "Wikipedia: T↓CTAGA, 5' overhang")]
    [TestCase("SalI", "GTCGAC", 1, 5, false, Description = "Wikipedia: G↓TCGAC, 5' overhang")]
    [TestCase("BglII", "AGATCT", 1, 5, false, Description = "Wikipedia: A↓GATCT, 5' overhang")]
    [TestCase("SpeI", "ACTAGT", 1, 5, false, Description = "Wikipedia: A↓CTAGT, 5' overhang")]
    [TestCase("ScaI", "AGTACT", 3, 3, true, Description = "Wikipedia: AGT↓ACT, Blunt")]
    [TestCase("StuI", "AGGCCT", 3, 3, true, Description = "Wikipedia: AGG↓CCT, Blunt")]
    public void EnzymeDatabase_MatchesWikipediaData(
        string name, string recognition, int cutForward, int cutReverse, bool isBlunt)
    {
        var enzyme = RestrictionAnalyzer.GetEnzyme(name);

        Assert.Multiple(() =>
        {
            Assert.That(enzyme, Is.Not.Null, $"Enzyme {name} should exist in database");
            Assert.That(enzyme!.RecognitionSequence, Is.EqualTo(recognition),
                $"{name} recognition sequence (Wikipedia)");
            Assert.That(enzyme.CutPositionForward, Is.EqualTo(cutForward),
                $"{name} forward cut position (Wikipedia)");
            Assert.That(enzyme.CutPositionReverse, Is.EqualTo(cutReverse),
                $"{name} reverse cut position (Wikipedia)");
            Assert.That(enzyme.IsBluntEnd, Is.EqualTo(isBlunt),
                $"{name} blunt end status (Wikipedia)");
        });
    }

    /// <summary>
    /// Evidence: FindCompatibleEnzymes returns known pairs.
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

            // HindIII (5' AGCT) and SacI (3' AGCT) must NOT appear as compatible
            bool hasHindIIISacI = compatiblePairs.Any(p =>
                (p.Enzyme1 == "HindIII" && p.Enzyme2 == "SacI") ||
                (p.Enzyme1 == "SacI" && p.Enzyme2 == "HindIII"));
            Assert.That(hasHindIIISacI, Is.False,
                "HindIII/SacI must NOT be listed — different overhang types (Wikipedia)");
        });
    }

    /// <summary>
    /// Evidence: All returned pairs should be actually compatible.
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
