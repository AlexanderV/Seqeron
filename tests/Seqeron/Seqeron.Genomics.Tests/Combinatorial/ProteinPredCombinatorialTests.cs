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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: DISORDER-LC-001 — Low-complexity region detection (SEG) (ProteinPred)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 204.
    // Spec: tests/TestSpecs/DISORDER-LC-001.md (canonical PredictLowComplexityRegions). ADVANCED §10.
    // Dimensions: windowSize(3) × threshold(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Wootton & Federhen 1993, SEG): a window with Shannon entropy ≤ the trigger cutoff K1
    // starts a low-complexity segment, extended over windows with entropy ≤ K2. A homopolymer window
    // has entropy 0.
    //
    // Axis mapping (documented): windowSize → the SEG trigger window; threshold → the trigger cutoff K1
    // (extension cutoff held equal so the grid uses a single complexity cutoff). The combinatorial
    // point: a planted poly-A block (entropy 0) is reported at every window/threshold, with valid
    // coordinates covering the block, while diverse flanks are not low-complexity.
    // ═══════════════════════════════════════════════════════════════════════

    private const string DiverseProteinFlank = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"; // all-20 AA cycle

    [Test, Combinatorial]
    public void DisorderLc_DetectsHomopolymerBlock_AcrossWindowAndThreshold(
        [Values(6, 12, 18)] int windowSize,
        [Values(0.5, 1.0, 2.0)] double threshold)
    {
        string protein = DiverseProteinFlank + new string('A', 30) + DiverseProteinFlank;
        int blockStart = DiverseProteinFlank.Length, blockEnd = blockStart + 29;

        var regions = DisorderPredictor.PredictLowComplexityRegions(
            protein, windowSize, triggerThreshold: threshold, extensionThreshold: threshold, minLength: 1).ToList();

        regions.Should().NotBeEmpty("a 30-residue poly-A block (entropy 0) is low-complexity at any K1 ≥ 0");
        regions.Should().OnlyContain(r => r.Start >= 0 && r.Start <= r.End && r.End < protein.Length, "valid coordinates");
        regions.Should().Contain(r => r.Start <= blockEnd && r.End >= blockStart, "a region covers the poly-A block");
    }

    /// <summary>
    /// Interaction witness — the complexity threshold gates: a two-residue "AB" repeat (entropy 1 bit)
    /// is low-complexity at K1 = 1.5 but not 0.5; a fully diverse window is never low-complexity.
    /// </summary>
    [Test]
    public void DisorderLc_ComplexityThresholdGates()
    {
        string ab = string.Concat(Enumerable.Repeat("AB", 20)); // entropy 1 bit/residue

        DisorderPredictor.PredictLowComplexityRegions(ab, 12, 1.5, 1.5, 1).Should().NotBeEmpty("entropy 1 ≤ 1.5");
        DisorderPredictor.PredictLowComplexityRegions(ab, 12, 0.5, 0.5, 1).Should().BeEmpty("entropy 1 > 0.5");

        DisorderPredictor.PredictLowComplexityRegions(DiverseProteinFlank, 12, 2.0, 2.0, 1)
            .Should().BeEmpty("a fully diverse window is not low-complexity");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: DISORDER-MORF-001 — Molecular-recognition-feature (MoRF) detection (ProteinPred)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 205.
    // Spec: tests/TestSpecs/DISORDER-MORF-001.md (canonical PredictMoRFs). ADVANCED §10.
    // Dimensions: windowSize(3) × threshold(3) × seqLen(3). Grid 3×3×3 = 27 (full, exhaustive).
    //
    // Model (Cheng/Oldfield; Mohan 2006): a MoRF is a short ORDER "dip" (disorder score < 0.5)
    // embedded within a longer disordered region, with length in [minLength, maxLength].
    //
    // Axis mapping (documented — PredictMoRFs takes only minLength/maxLength): windowSize → minLength;
    // threshold → maxLength; seqLen → the flanking disorder length. Engineered construct: a 20-residue
    // order dip (poly-W) flanked by disorder (poly-P). The combinatorial point: the dip is reported as
    // a MoRF exactly when its length 20 lies in [minLength, maxLength], independent of flank length.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void DisorderMorf_DetectsOrderDip_AcrossMinLenMaxLenFlank(
        [Values(10, 18, 25)] int minLength,
        [Values(22, 40, 70)] int maxLength,
        [Values(20, 40, 60)] int flankLen)
    {
        const int dipLen = 20;
        string protein = new string('P', flankLen) + new string('W', dipLen) + new string('P', flankLen);
        int dipStart = flankLen, dipEnd = flankLen + dipLen - 1;

        var morfs = DisorderPredictor.PredictMoRFs(protein, minLength, maxLength).ToList();

        bool expectDetected = minLength <= dipLen && dipLen <= maxLength;
        bool coversDip = morfs.Any(m => m.Start <= dipEnd && m.End >= dipStart);
        coversDip.Should().Be(expectDetected, "the 20-residue order dip is a MoRF iff 20 ∈ [minLength, maxLength]");

        morfs.Should().OnlyContain(m => m.Score >= -1e-9 && m.Score <= 1.0 + 1e-9, "the dip-depth score is in [0,1]");
    }

    /// <summary>
    /// Interaction witness — a pure-disorder protein has no MoRF (no order dip), and the length band
    /// gates a planted dip (too-long minLength removes it).
    /// </summary>
    [Test]
    public void DisorderMorf_NoDipNoMorf_AndLengthBandGates()
    {
        DisorderPredictor.PredictMoRFs(new string('P', 80)).Should().BeEmpty("a fully disordered protein has no order dip");

        string protein = new string('P', 40) + new string('W', 20) + new string('P', 40);
        DisorderPredictor.PredictMoRFs(protein, minLength: 10, maxLength: 70)
            .Should().Contain(m => m.Start <= 59 && m.End >= 40, "a 20-residue dip is a MoRF in the 10–70 band");
        DisorderPredictor.PredictMoRFs(protein, minLength: 30, maxLength: 70)
            .Should().NotContain(m => m.Start <= 59 && m.End >= 40, "a 20-residue dip is too short for minLength 30");
    }
}
