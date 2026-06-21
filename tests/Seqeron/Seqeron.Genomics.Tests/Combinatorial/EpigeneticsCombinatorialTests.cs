namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Epigenetics area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Epigenetics")]
public class EpigeneticsCombinatorialTests
{
    // Strong CpG island (GC = 1.0, O/E ≈ 2) flanked by AT-rich, CpG-poor filler.
    private static readonly string CpgSeq =
        string.Concat(Enumerable.Repeat("AT", 50)) + string.Concat(Enumerable.Repeat("CG", 200)) + string.Concat(Enumerable.Repeat("AT", 50));

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: EPIGEN-CPG-001 — CpG island detection (Epigenetics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 85.
    // Spec: tests/TestSpecs/EPIGEN-CPG-001.md (canonical FindCpGIslands).
    // Dimensions: windowSize(3) × minOE(3) × minGC(3) × minLen(3). Grid 3⁴ = 81.
    //
    // Model (Gardiner-Garden & Frommer 1987): a CpG island is a region of length ≥ minLength with
    // GC content ≥ minGc and observed/expected CpG ratio ≥ minCpGRatio. The detection window equals
    // minLength internally, so windowSize is realised as the searched-substring length.
    //
    // The combinatorial point: the searched window and the three criteria jointly constrain the
    // output — every reported island simultaneously satisfies the GC, O/E and length thresholds,
    // for all 81 parameter combinations, and lies within the searched window.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void EpigenCpg_ReportedIslandsSatisfyAllCriteria(
        [Values(150, 350, 600)] int windowSize,
        [Values(0.4, 0.6, 1.0)] double minOE,
        [Values(0.4, 0.6, 0.8)] double minGC,
        [Values(100, 200, 300)] int minLen)
    {
        string searched = CpgSeq[..windowSize];
        var islands = EpigeneticsAnalyzer.FindCpGIslands(searched, minLen, minGC, minOE).ToList();

        islands.Should().OnlyContain(isl => isl.GcContent >= minGC, "island GC clears the threshold");
        islands.Should().OnlyContain(isl => isl.CpGRatio >= minOE, "island O/E clears the threshold");
        islands.Should().OnlyContain(isl => isl.End - isl.Start >= minLen, "island length clears the minimum");
        islands.Should().OnlyContain(isl => isl.Start >= 0 && isl.End <= searched.Length, "coordinates within the searched window");
    }

    /// <summary>
    /// Interaction witness: a strong CpG island is detected under permissive criteria, while an
    /// unreachable O/E or length floor removes it.
    /// </summary>
    [Test]
    public void EpigenCpg_StrongIsland_DetectedAndGated()
    {
        EpigeneticsAnalyzer.FindCpGIslands(CpgSeq, 200, 0.4, 0.4).Should().NotBeEmpty("a CG-rich island passes permissive criteria");
        EpigeneticsAnalyzer.FindCpGIslands(CpgSeq, 200, 0.4, 5.0).Should().BeEmpty("no region reaches an O/E of 5");
        EpigeneticsAnalyzer.FindCpGIslands(CpgSeq, 5000, 0.4, 0.4).Should().BeEmpty("no 5000-bp island fits in the sequence");
    }

    /// <summary>
    /// Interaction witness: AT-rich, CpG-poor sequence has no islands regardless of length.
    /// </summary>
    [Test]
    public void EpigenCpg_AtRichSequence_NoIslands()
    {
        EpigeneticsAnalyzer.FindCpGIslands(string.Concat(Enumerable.Repeat("AT", 200)), 200, 0.5, 0.6)
            .Should().BeEmpty("an AT-rich tract is not a CpG island");
    }
}
