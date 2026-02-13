// EPIGEN-CPG-001 — CpG Site Detection
// Evidence: docs/Evidence/EPIGEN-CPG-001-Evidence.md
// TestSpec: tests/TestSpecs/EPIGEN-CPG-001.md
// Source: Gardiner-Garden M, Frommer M (1987). J Mol Biol. 196(2):261–282.

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class EpigeneticsAnalyzer_CpGDetection_Tests
{
    #region FindCpGSites Tests

    // M1 — FindCpGSites exact positions
    // Evidence: CpG = C immediately followed by G in 5'→3' direction
    [Test]
    public void FindCpGSites_SimpleCpG_ReturnsExactPositions()
    {
        // "ACGTCGACG" → CpG at positions 1, 4, 7
        //  A C G T C G A C G
        //  0 1 2 3 4 5 6 7 8
        string sequence = "ACGTCGACG";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Has.Count.EqualTo(3), "Expected exactly 3 CpG sites in ACGTCGACG");
        Assert.That(sites[0], Is.EqualTo(1), "First CpG at position 1 (CG at 1-2)");
        Assert.That(sites[1], Is.EqualTo(4), "Second CpG at position 4 (CG at 4-5)");
        Assert.That(sites[2], Is.EqualTo(7), "Third CpG at position 7 (CG at 7-8)");
    }

    // M2 — Adjacent CpG dinucleotides
    // Evidence: "CGCGCG" contains 3 CpG sites at positions 0, 2, 4
    [Test]
    public void FindCpGSites_AdjacentCpG_ReturnsAllPositions()
    {
        // C G C G C G
        // 0 1 2 3 4 5
        string sequence = "CGCGCG";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Has.Count.EqualTo(3), "Expected 3 CpG sites in CGCGCG");
        Assert.That(sites[0], Is.EqualTo(0), "CpG at position 0");
        Assert.That(sites[1], Is.EqualTo(2), "CpG at position 2");
        Assert.That(sites[2], Is.EqualTo(4), "CpG at position 4");
    }

    // M3 — No CpG in sequence
    [Test]
    public void FindCpGSites_NoCpG_ReturnsEmpty()
    {
        string sequence = "AATTAATT";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Is.Empty, "AT-only sequence should have no CpG sites");
    }

    // M4 — Null and empty input
    [Test]
    public void FindCpGSites_NullSequence_ReturnsEmpty()
    {
        var sites = EpigeneticsAnalyzer.FindCpGSites(null!).ToList();

        Assert.That(sites, Is.Empty, "Null sequence should return empty");
    }

    [Test]
    public void FindCpGSites_EmptySequence_ReturnsEmpty()
    {
        var sites = EpigeneticsAnalyzer.FindCpGSites("").ToList();

        Assert.That(sites, Is.Empty, "Empty sequence should return empty");
    }

    // M5 — Case insensitive
    // Evidence: Standard convention — sequences may be lowercase
    [Test]
    public void FindCpGSites_LowercaseSequence_ReturnsExactPositions()
    {
        // "acgtcg" → CpG at positions 1, 4
        string sequence = "acgtcg";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Has.Count.EqualTo(2), "Expected 2 CpG sites in lowercase acgtcg");
        Assert.That(sites[0], Is.EqualTo(1), "First CpG at position 1");
        Assert.That(sites[1], Is.EqualTo(4), "Second CpG at position 4");
    }

    // M6 — CpG at start of sequence
    [Test]
    public void FindCpGSites_CpGAtStart_Detected()
    {
        string sequence = "CGAAA";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Has.Count.EqualTo(1), "Expected 1 CpG site");
        Assert.That(sites[0], Is.EqualTo(0), "CpG at position 0 (start of sequence)");
    }

    // M7 — CpG at end of sequence
    [Test]
    public void FindCpGSites_CpGAtEnd_Detected()
    {
        string sequence = "AACG";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Has.Count.EqualTo(1), "Expected 1 CpG site");
        Assert.That(sites[0], Is.EqualTo(2), "CpG at position 2 (end of sequence)");
    }

    // M8 — Single nucleotide
    [Test]
    public void FindCpGSites_SingleNucleotide_ReturnsEmpty()
    {
        var sites = EpigeneticsAnalyzer.FindCpGSites("C").ToList();

        Assert.That(sites, Is.Empty, "Single nucleotide cannot form a dinucleotide");
    }

    // M18 — GpC is NOT CpG
    // Evidence: CpG is 5'→3' C then G; GC is NOT a CpG site
    [Test]
    public void FindCpGSites_GpCNotCountedAsCpG()
    {
        // G C G C G C → GpC at 0, CpG at 1, GpC at 2, CpG at 3, GpC at 4
        string sequence = "GCGCGC";

        var sites = EpigeneticsAnalyzer.FindCpGSites(sequence).ToList();

        Assert.That(sites, Has.Count.EqualTo(2), "Only CpG (not GpC) should be counted");
        Assert.That(sites[0], Is.EqualTo(1), "CpG at position 1");
        Assert.That(sites[1], Is.EqualTo(3), "CpG at position 3");
    }

    // S1 — Minimal CG dinucleotide
    [Test]
    public void FindCpGSites_MinimalCG_ReturnsOnePosition()
    {
        var sites = EpigeneticsAnalyzer.FindCpGSites("CG").ToList();

        Assert.That(sites, Has.Count.EqualTo(1), "Minimal CG should yield 1 site");
        Assert.That(sites[0], Is.EqualTo(0), "CpG at position 0");
    }

    #endregion

    #region CalculateCpGObservedExpected Tests

    // M9 — O/E for pure CG repeat
    // Evidence: Gardiner-Garden & Frommer (1987) formula:
    //   O/E = CpG_count / (C_count × G_count / Length)
    //   "CGCGCGCGCGCGCGCGCGCG" (20 bp): CpG=10, C=10, G=10
    //   expected = 10*10/20 = 5.0 → O/E = 10/5 = 2.0
    [Test]
    public void CalculateCpGObservedExpected_PureCGRepeat_ReturnsExact2()
    {
        string sequence = "CGCGCGCGCGCGCGCGCGCG";

        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(sequence);

        Assert.That(ratio, Is.EqualTo(2.0).Within(1e-10),
            "O/E for pure CG repeat: CpG=10, exp=10×10/20=5, O/E=10/5=2.0");
    }

    // M10 — O/E for AT-only
    // Evidence: No CpG, no C or G → O/E = 0
    [Test]
    public void CalculateCpGObservedExpected_ATOnly_ReturnsZero()
    {
        string sequence = "AATTAATTAATTAATTAATT";

        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(sequence);

        Assert.That(ratio, Is.EqualTo(0.0), "AT-only sequence: no CpG, O/E = 0");
    }

    // M11 — O/E exact calculation for mixed sequence
    // Evidence: "ACGTCGACG" (9 bp): CpG=3, C=3, G=3
    //   expected = 3*3/9 = 1.0 → O/E = 3/1 = 3.0
    [Test]
    public void CalculateCpGObservedExpected_MixedSequence_ReturnsExact3()
    {
        string sequence = "ACGTCGACG";

        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(sequence);

        Assert.That(ratio, Is.EqualTo(3.0).Within(1e-10),
            "O/E for ACGTCGACG: CpG=3, exp=3×3/9=1.0, O/E=3/1=3.0");
    }

    // M12 — O/E minimal sequence
    // Evidence: "ACGT" (4 bp): CpG=1, C=1, G=1
    //   expected = 1*1/4 = 0.25 → O/E = 1/0.25 = 4.0
    [Test]
    public void CalculateCpGObservedExpected_MinimalSequence_ReturnsExact4()
    {
        string sequence = "ACGT";

        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(sequence);

        Assert.That(ratio, Is.EqualTo(4.0).Within(1e-10),
            "O/E for ACGT: CpG=1, exp=1×1/4=0.25, O/E=1/0.25=4.0");
    }

    // M13 — Null input
    [Test]
    public void CalculateCpGObservedExpected_NullSequence_ReturnsZero()
    {
        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(null!);

        Assert.That(ratio, Is.EqualTo(0.0), "Null sequence should return 0");
    }

    // M14 — Single character
    [Test]
    public void CalculateCpGObservedExpected_SingleChar_ReturnsZero()
    {
        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected("A");

        Assert.That(ratio, Is.EqualTo(0.0), "Single character: length < 2, O/E = 0");
    }

    // S2 — O/E case insensitive
    [Test]
    public void CalculateCpGObservedExpected_Lowercase_ReturnsSameAsUppercase()
    {
        string lower = "cgcgcgcgcgcgcgcgcgcg";
        string upper = "CGCGCGCGCGCGCGCGCGCG";

        double ratioLower = EpigeneticsAnalyzer.CalculateCpGObservedExpected(lower);
        double ratioUpper = EpigeneticsAnalyzer.CalculateCpGObservedExpected(upper);

        Assert.That(ratioLower, Is.EqualTo(ratioUpper).Within(1e-10),
            "Case should not affect O/E calculation");
        Assert.That(ratioLower, Is.EqualTo(2.0).Within(1e-10));
    }

    // C1 — Only C, no G → expected = 0, O/E = 0
    [Test]
    public void CalculateCpGObservedExpected_OnlyCNoG_ReturnsZero()
    {
        double ratio = EpigeneticsAnalyzer.CalculateCpGObservedExpected("CCCC");

        Assert.That(ratio, Is.EqualTo(0.0),
            "No G nucleotides → expected=0 → guarded → O/E=0");
    }

    #endregion

    #region FindCpGIslands Tests

    // M15 — Positive CpG island detection
    // Evidence: Gardiner-Garden & Frommer (1987) criteria: ≥200 bp, GC% > 50%, O/E > 0.6
    [Test]
    public void FindCpGIslands_CpGRichRegion_DetectsIsland()
    {
        // 400 bp CG repeat — GC=100%, O/E=2.0, length=400 ≥ 200
        string cpgIsland = string.Concat(Enumerable.Repeat("CGCG", 100));

        var islands = EpigeneticsAnalyzer.FindCpGIslands(
            cpgIsland,
            minLength: 100,
            minGc: 0.5,
            minCpGRatio: 0.6).ToList();

        Assert.That(islands, Has.Count.GreaterThanOrEqualTo(1),
            "400 bp CGCG repeat should produce at least 1 CpG island");
        Assert.That(islands[0].GcContent, Is.EqualTo(1.0).Within(1e-10),
            "GC content of pure CG sequence should be 1.0");
        Assert.That(islands[0].CpGRatio, Is.GreaterThanOrEqualTo(0.6),
            "CpG O/E ratio must meet minimum threshold");
    }

    // M16 — No CpG island in AT-rich sequence
    [Test]
    public void FindCpGIslands_ATRichSequence_ReturnsEmpty()
    {
        string atRich = string.Concat(Enumerable.Repeat("AATT", 100));

        var islands = EpigeneticsAnalyzer.FindCpGIslands(atRich).ToList();

        Assert.That(islands, Is.Empty, "AT-rich sequence should not contain CpG islands");
    }

    // M17 — Sequence shorter than minimum island length
    [Test]
    public void FindCpGIslands_ShortSequence_ReturnsEmpty()
    {
        string shortSeq = "CGCGCG";

        var islands = EpigeneticsAnalyzer.FindCpGIslands(shortSeq, minLength: 200).ToList();

        Assert.That(islands, Is.Empty,
            "Sequence shorter than minLength should yield no islands");
    }

    // S3 — Island result contains valid GcContent and CpGRatio
    [Test]
    public void FindCpGIslands_ResultContainsValidMetrics()
    {
        string cpgIsland = string.Concat(Enumerable.Repeat("CGCG", 100));

        var islands = EpigeneticsAnalyzer.FindCpGIslands(
            cpgIsland,
            minLength: 100,
            minGc: 0.5,
            minCpGRatio: 0.6).ToList();

        Assert.That(islands, Is.Not.Empty, "Should detect at least 1 island");
        var island = islands[0];
        Assert.That(island.Start, Is.GreaterThanOrEqualTo(0), "Start must be non-negative");
        Assert.That(island.End, Is.GreaterThan(island.Start), "End must be after Start");
        Assert.That(island.GcContent, Is.GreaterThan(0).And.LessThanOrEqualTo(1.0),
            "GC content must be in (0, 1]");
        Assert.That(island.CpGRatio, Is.GreaterThan(0),
            "CpG O/E ratio must be positive for a detected island");
    }

    // S4 — Custom parameters are respected
    [Test]
    public void FindCpGIslands_CustomParameters_Respected()
    {
        // 200 bp CG repeat
        string cpgIsland = string.Concat(Enumerable.Repeat("CGCG", 50));

        // With default minLength=200 and minGc=0.5 → should find island
        var islandsDefault = EpigeneticsAnalyzer.FindCpGIslands(
            cpgIsland, minLength: 200, minGc: 0.5, minCpGRatio: 0.6).ToList();

        // With impossible minGc=1.1 → should find nothing (unreachable threshold)
        var islandsImpossible = EpigeneticsAnalyzer.FindCpGIslands(
            cpgIsland, minLength: 200, minGc: 1.1, minCpGRatio: 0.6).ToList();

        Assert.That(islandsDefault, Is.Not.Empty,
            "Default criteria should detect island in CG repeat");
        Assert.That(islandsImpossible, Is.Empty,
            "Impossible GC threshold should filter out all windows");
    }

    // C2 — Null input
    [Test]
    public void FindCpGIslands_NullSequence_ReturnsEmpty()
    {
        var islands = EpigeneticsAnalyzer.FindCpGIslands(null!).ToList();

        Assert.That(islands, Is.Empty, "Null sequence should return empty");
    }

    #endregion
}
