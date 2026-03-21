using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for intrinsically disordered protein prediction.
/// Verifies score range, length preservation, and determinism invariants.
///
/// Test Units: DISORDER-PRED-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("ProteinPred")]
public class DisorderProperties
{
    #region Generators

    /// <summary>
    /// Generates random protein sequences from the 20 standard amino acids.
    /// </summary>
    private static Arbitrary<string> ProteinArbitrary(int minLen = 25) =>
        Gen.Elements('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                     'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// A protein rich in disorder-promoting residues (P, E, S, Q, K, A, G, R).
    /// Evidence: Dunker et al. (2001) classified these as disorder-promoting.
    /// </summary>
    private static string DisorderRichProtein =>
        new(Enumerable.Repeat("PPEESSKK", 10).SelectMany(s => s).ToArray());

    /// <summary>
    /// A protein rich in order-promoting residues (W, C, F, I, Y, V, L, N).
    /// Evidence: Dunker et al. (2001) classified these as order-promoting.
    /// </summary>
    private static string OrderRichProtein =>
        new(Enumerable.Repeat("WWCCFFII", 10).SelectMany(s => s).ToArray());

    #endregion

    #region DISORDER-PRED-001: R: score ∈ [0,1]; P: len(scores) = len(sequence); D: deterministic

    /// <summary>
    /// INV-1: All per-residue disorder scores are in [0, 1].
    /// Evidence: Score = average of normalized TOP-IDP values over sliding window.
    /// TOP-IDP is normalized to [0, 1] before averaging, so the mean is also in [0, 1].
    /// Source: Campen et al. (2008) — TOP-IDP scale; normalization bounds the output.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderPrediction_Scores_InRange()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            return result.ResiduePredictions.All(rp => rp.DisorderScore >= -1e-9 && rp.DisorderScore <= 1.0 + 1e-9)
                .Label("All disorder scores must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-2: Number of per-residue predictions equals sequence length.
    /// Evidence: The sliding window produces one score per residue position.
    /// Source: Campen et al. (2008) — one prediction per residue.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderPrediction_PredictionCount_EqualsSequenceLength()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            return (result.ResiduePredictions.Count == seq.Length)
                .Label($"Predictions count={result.ResiduePredictions.Count}, seqLen={seq.Length}");
        });
    }

    /// <summary>
    /// INV-3: OverallDisorderContent ∈ [0, 1].
    /// Evidence: Content = fraction of residues classified as disordered, bounded by 0 and 1.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderPrediction_OverallContent_InRange()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            return (result.OverallDisorderContent >= -1e-9 && result.OverallDisorderContent <= 1.0 + 1e-9)
                .Label($"OverallDisorderContent={result.OverallDisorderContent:F4} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-4: MeanDisorderScore ∈ [0, 1].
    /// Evidence: Mean of per-residue scores, each in [0, 1], is also in [0, 1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderPrediction_MeanScore_InRange()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            return (result.MeanDisorderScore >= -1e-9 && result.MeanDisorderScore <= 1.0 + 1e-9)
                .Label($"MeanDisorderScore={result.MeanDisorderScore:F4} must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-5: Disordered regions have valid coordinates — start &lt; end ≤ seqLen.
    /// Evidence: Regions are contiguous stretches of disordered residues.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderPrediction_RegionCoordinates_Valid()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            return result.DisorderedRegions.All(r =>
                r.Start >= 0 && r.Start < r.End && r.End <= seq.Length)
                .Label("Disordered region coordinates must satisfy 0 ≤ start < end ≤ seqLen");
        });
    }

    /// <summary>
    /// INV-6: Disorder prediction is deterministic.
    /// Evidence: PredictDisorder uses a deterministic sliding window with fixed TOP-IDP scale.
    /// Source: Campen et al. (2008).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderPrediction_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var r1 = DisorderPredictor.PredictDisorder(seq);
            var r2 = DisorderPredictor.PredictDisorder(seq);
            bool same = r1.ResiduePredictions.Count == r2.ResiduePredictions.Count &&
                        r1.ResiduePredictions.Zip(r2.ResiduePredictions)
                            .All(pair => Math.Abs(pair.First.DisorderScore - pair.Second.DisorderScore) < 1e-10) &&
                        Math.Abs(r1.MeanDisorderScore - r2.MeanDisorderScore) < 1e-10;
            return same.Label("PredictDisorder must be deterministic");
        });
    }

    /// <summary>
    /// INV-7: Per-residue positions are sequential from 0 to N-1.
    /// Evidence: Each residue is assigned a unique sequential position index.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderPrediction_Positions_AreSequential()
    {
        return Prop.ForAll(ProteinArbitrary(), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            bool sequential = result.ResiduePredictions
                .Select((rp, i) => rp.Position == i)
                .All(x => x);
            return sequential.Label("Residue positions must be sequential 0..N-1");
        });
    }

    /// <summary>
    /// INV-8: Disorder-rich protein (mostly P, E, S, K) has higher mean disorder score
    /// than order-rich protein (mostly W, C, F, I).
    /// Evidence: TOP-IDP propensity is &gt; 0 for disorder-promoting and &lt; 0 for
    /// order-promoting amino acids. After normalization and averaging, disorder-rich
    /// sequences score higher.
    /// Source: Dunker et al. (2001) for classification; Campen et al. (2008) for TOP-IDP.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DisorderPrediction_DisorderRich_ScoresHigherThanOrderRich()
    {
        var disorderResult = DisorderPredictor.PredictDisorder(DisorderRichProtein);
        var orderResult = DisorderPredictor.PredictDisorder(OrderRichProtein);

        Assert.That(disorderResult.MeanDisorderScore, Is.GreaterThan(orderResult.MeanDisorderScore),
            $"Disorder-rich mean={disorderResult.MeanDisorderScore:F4} should be > " +
            $"order-rich mean={orderResult.MeanDisorderScore:F4}");
    }

    #endregion
}
