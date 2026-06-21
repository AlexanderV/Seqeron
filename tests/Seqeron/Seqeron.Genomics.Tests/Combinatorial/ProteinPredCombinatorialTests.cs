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
}
