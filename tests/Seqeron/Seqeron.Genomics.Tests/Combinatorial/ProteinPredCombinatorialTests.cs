namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the ProteinPred area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("ProteinPred")]
public class ProteinPredCombinatorialTests
{
    private const string AaPool = "SPEKQNDGARWFILVYHTMC";
    private static string Protein(int n) => string.Concat(Enumerable.Range(0, n).Select(i => AaPool[(i * 7 + 3) % AaPool.Length]));

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: DISORDER-PRED-001 — Intrinsic-disorder prediction (ProteinPred)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 80.
    // Spec: tests/TestSpecs/DISORDER-PRED-001.md (canonical PredictDisorder).
    // Dimensions: windowSize(3) × propScale(2) × seqLen(3). Grid 3×2×3 = 18.
    //
    // Model (Campen 2008 TOP-IDP): each residue's disorder score is the window-averaged TOP-IDP
    // propensity (normalised to [0,1]); a residue is disordered when its score ≥ the threshold.
    // A single TOP-IDP scale is implemented, so the propScale axis is realised as the disorder
    // call threshold (the cutoff that defines disorder).
    //
    // The combinatorial point: window size, threshold and length interact — per-residue calls are
    // exactly score ≥ threshold, the disorder content is their fraction, and reported regions meet
    // the minimum length, for every (window, threshold, length) cell.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void DisorderPred_PerResidueCallsAndContent(
        [Values(7, 15, 21)] int windowSize,
        [Values(0.45, 0.55)] double threshold,
        [Values(30, 60, 120)] int seqLen)
    {
        string seq = Protein(seqLen);
        var r = DisorderPredictor.PredictDisorder(seq, windowSize, threshold, minRegionLength: 5);

        r.ResiduePredictions.Should().HaveCount(seqLen, "one prediction per residue");
        r.ResiduePredictions.Should().OnlyContain(rp => rp.DisorderScore >= 0 && rp.DisorderScore <= 1,
            "normalised TOP-IDP scores are in [0,1]");
        r.ResiduePredictions.Should().OnlyContain(rp => rp.IsDisordered == (rp.DisorderScore >= threshold),
            "a residue is disordered iff its score clears the threshold");

        double expectedContent = r.ResiduePredictions.Count(rp => rp.IsDisordered) / (double)seqLen;
        r.OverallDisorderContent.Should().BeApproximately(expectedContent, 1e-9);
        r.MeanDisorderScore.Should().BeApproximately(r.ResiduePredictions.Average(rp => rp.DisorderScore), 1e-9);

        foreach (var region in r.DisorderedRegions)
        {
            (region.End - region.Start + 1).Should().BeGreaterThanOrEqualTo(5, "regions meet minRegionLength");
            for (int p = region.Start; p <= region.End; p++)
                r.ResiduePredictions[p].IsDisordered.Should().BeTrue("region residues are disordered");
        }
    }

    /// <summary>
    /// Interaction witness: raising the disorder threshold can only reduce the disordered-residue
    /// set (monotone), at a fixed window.
    /// </summary>
    [Test]
    public void DisorderPred_Threshold_IsMonotone()
    {
        string seq = Protein(80);
        int low = DisorderPredictor.PredictDisorder(seq, 15, 0.45).ResiduePredictions.Count(rp => rp.IsDisordered);
        int high = DisorderPredictor.PredictDisorder(seq, 15, 0.55).ResiduePredictions.Count(rp => rp.IsDisordered);
        high.Should().BeLessThanOrEqualTo(low, "a higher threshold flags fewer residues");
    }

    /// <summary>
    /// Interaction witness: a disorder-promoting homopolymer (proline) has higher disorder content
    /// than an order-promoting one (tryptophan).
    /// </summary>
    [Test]
    public void DisorderPred_DisorderPromotingResidues_HigherContent()
    {
        double pro = DisorderPredictor.PredictDisorder(new string('P', 40)).OverallDisorderContent;
        double trp = DisorderPredictor.PredictDisorder(new string('W', 40)).OverallDisorderContent;
        pro.Should().BeGreaterThan(trp, "proline is disorder-promoting, tryptophan is order-promoting");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: DISORDER-REGION-001 — Disordered-region segmentation (ProteinPred)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 81.
    // Spec: tests/TestSpecs/DISORDER-REGION-001.md (canonical PredictDisorder regions).
    // Dimensions: threshold(3) × minLen(3) × mergeGap(3). Grid 3×3×3 = 27.
    //
    // Model: a disordered region is a maximal CONTIGUOUS run of disordered residues of length
    // ≥ minRegionLength. Region merging across gaps is NOT a parameter here — regions break at the
    // first ordered residue — so the mergeGap axis is realised as the size of an ordered gap
    // planted between two disordered blocks (a clear gap must NOT be merged away).
    //
    // The combinatorial point: threshold, minimum length and gap size interact — every reported
    // region is a contiguous disordered run meeting the length floor, regions never span an
    // ordered residue, and a clear ordered gap keeps two blocks separate.
    // ═══════════════════════════════════════════════════════════════════════

    private static string GappedDisorder(int gap) => new string('P', 30) + new string('W', gap) + new string('P', 30);

    [Test, Combinatorial]
    public void DisorderRegion_ContiguousRunsHonourBounds(
        [Values(0.45, 0.55, 0.65)] double threshold,
        [Values(3, 5, 10)] int minLen,
        [Values(1, 5, 15)] int gap)
    {
        string seq = GappedDisorder(gap);
        var r = DisorderPredictor.PredictDisorder(seq, windowSize: 7, threshold, minRegionLength: minLen);

        DisorderPredictor.DisorderedRegion? prev = null;
        foreach (var region in r.DisorderedRegions)
        {
            (region.End - region.Start + 1).Should().BeGreaterThanOrEqualTo(minLen, "region meets the length floor");
            for (int p = region.Start; p <= region.End; p++)
                r.ResiduePredictions[p].IsDisordered.Should().BeTrue("a region never spans an ordered residue");
            if (prev.HasValue)
                region.Start.Should().BeGreaterThan(prev.Value.End, "regions are ordered and disjoint");
            prev = region;
        }
    }

    /// <summary>
    /// Interaction witness: a clear ordered gap (15 tryptophans) between two proline blocks is not
    /// merged — at least two separate disordered regions are reported, and none spans the gap.
    /// </summary>
    [Test]
    public void DisorderRegion_ClearGap_NotMerged()
    {
        string seq = GappedDisorder(15);
        var r = DisorderPredictor.PredictDisorder(seq, windowSize: 7, 0.55, minRegionLength: 5);

        r.DisorderedRegions.Should().HaveCountGreaterThanOrEqualTo(2, "a clear ordered gap splits the disorder");
        r.DisorderedRegions.Should().OnlyContain(reg =>
            Enumerable.Range(reg.Start, reg.End - reg.Start + 1).All(p => r.ResiduePredictions[p].IsDisordered));
    }

    /// <summary>
    /// Interaction witness: raising minRegionLength can only drop short regions (monotone in the
    /// reported region set).
    /// </summary>
    [Test]
    public void DisorderRegion_MinLength_IsMonotone()
    {
        string seq = GappedDisorder(8);
        int many = DisorderPredictor.PredictDisorder(seq, 7, 0.55, minRegionLength: 3).DisorderedRegions.Count;
        int few = DisorderPredictor.PredictDisorder(seq, 7, 0.55, minRegionLength: 20).DisorderedRegions.Count;
        few.Should().BeLessThanOrEqualTo(many, "a larger minimum length yields no more regions");
    }
}
