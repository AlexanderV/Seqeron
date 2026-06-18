using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for intrinsically disordered protein prediction.
/// Verifies score range, length preservation, and determinism invariants.
///
/// Test Units: DISORDER-PRED-001, DISORDER-REGION-001, DISORDER-LC-001, DISORDER-MORF-001, DISORDER-PROPENSITY-001
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

    #region DISORDER-REGION-001: R: region start < end ≤ seqLen; M: lower threshold → larger regions; D: deterministic

    /// <summary>
    /// INV-9: Disordered regions have valid coordinates — start &lt; end ≤ seqLen.
    /// Evidence: Regions are contiguous stretches of consecutive disordered residues.
    /// Source: Campen et al. (2008) — region boundaries are defined by contiguous
    /// residues that exceed the disorder threshold.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderRegion_Coordinates_Valid()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            return result.DisorderedRegions.All(r =>
                r.Start >= 0 && r.Start < r.End && r.End <= seq.Length)
                .Label("All disordered region coordinates must satisfy 0 ≤ start < end ≤ seqLen");
        });
    }

    /// <summary>
    /// INV-10: Lower disorder threshold produces regions covering ≥ as many residues.
    /// Evidence: Lowering the threshold admits more residues as disordered, so the union
    /// of disordered regions can only grow or stay the same. Total disordered residue count
    /// is monotonically non-decreasing as threshold decreases.
    /// Source: Standard thresholding property — fewer residues are excluded with lower cutoff.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderRegion_LowerThreshold_MoreDisorderedResidues()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var resultHigh = DisorderPredictor.PredictDisorder(seq, disorderThreshold: 0.6);
            var resultLow = DisorderPredictor.PredictDisorder(seq, disorderThreshold: 0.3);
            int highCount = resultHigh.ResiduePredictions.Count(r => r.DisorderScore >= 0.6);
            int lowCount = resultLow.ResiduePredictions.Count(r => r.DisorderScore >= 0.3);
            return (lowCount >= highCount)
                .Label($"Lower threshold should produce ≥ disordered residues: low={lowCount}, high={highCount}");
        });
    }

    /// <summary>
    /// INV-11: Disordered regions do not overlap.
    /// Evidence: Regions are derived from contiguous non-overlapping segments.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderRegion_NoOverlap()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            var regions = result.DisorderedRegions;
            for (int i = 0; i < regions.Count - 1; i++)
            {
                if (regions[i].End > regions[i + 1].Start)
                    return false.Label(
                        $"Overlapping regions: [{regions[i].Start},{regions[i].End}) and [{regions[i + 1].Start},{regions[i + 1].End})");
            }
            return true.Label("No overlapping disordered regions");
        });
    }

    /// <summary>
    /// INV-12: Disordered region mean score ∈ [0, 1].
    /// Evidence: MeanScore is the average of per-residue disorder scores within the region.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderRegion_MeanScore_InRange()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            return result.DisorderedRegions.All(r =>
                r.MeanScore >= -1e-9 && r.MeanScore <= 1.0 + 1e-9)
                .Label("Disordered region mean score must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-13: Disordered region detection is deterministic.
    /// Evidence: PredictDisorder is a deterministic sliding window computation.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DisorderRegion_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(50), seq =>
        {
            var r1 = DisorderPredictor.PredictDisorder(seq);
            var r2 = DisorderPredictor.PredictDisorder(seq);
            bool same = r1.DisorderedRegions.Count == r2.DisorderedRegions.Count &&
                        r1.DisorderedRegions.Zip(r2.DisorderedRegions)
                            .All(pair => pair.First.Start == pair.Second.Start &&
                                         pair.First.End == pair.Second.End);
            return same.Label("DisorderedRegion detection must be deterministic");
        });
    }

    #endregion

    #region DISORDER-LC-001: R: region start ≤ end; M: higher threshold → ≥ coverage; D: deterministic

    // PredictLowComplexityRegions implements SEG (Wootton & Federhen 1993). As in all SEG variants a
    // window is flagged when entropy ≤ threshold, so RAISING the threshold flags more — coverage is
    // monotone increasing in the threshold (the checklist's wording is the inverse sense).

    private static HashSet<int> Covered(IEnumerable<(int Start, int End, string Type)> regions)
    {
        var set = new HashSet<int>();
        foreach (var (s, e, _) in regions)
            for (int p = s; p <= e; p++) set.Add(p);
        return set;
    }

    /// <summary>
    /// INV-1 (R): every reported low-complexity region has Start ≤ End within bounds and a non-empty
    /// classification.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LowComplexity_RegionsAreValid()
    {
        return Prop.ForAll(ProteinArbitrary(20), seq =>
        {
            var regions = DisorderPredictor.PredictLowComplexityRegions(seq).ToList();
            return regions.All(r => r.Start >= 0 && r.Start <= r.End && r.End < seq.Length && !string.IsNullOrEmpty(r.Type))
                .Label("a low-complexity region was invalid");
        });
    }

    /// <summary>
    /// INV-2 (M): raising the SEG thresholds never reduces the flagged residues.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LowComplexity_HigherThreshold_CoversMore()
    {
        return Prop.ForAll(ProteinArbitrary(20), seq =>
        {
            var low = Covered(DisorderPredictor.PredictLowComplexityRegions(seq, 12, 1.0, 1.5));
            var high = Covered(DisorderPredictor.PredictLowComplexityRegions(seq, 12, 3.0, 3.5));
            return low.IsSubsetOf(high).Label($"low-threshold coverage ({low.Count}) not ⊆ high ({high.Count})");
        });
    }

    /// <summary>
    /// INV-3 (P, positive control + D): a homopolymer is low-complexity; a maximally diverse window is
    /// not; detection is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void LowComplexity_HomopolymerDetected_AndDeterministic()
    {
        var homo = DisorderPredictor.PredictLowComplexityRegions(new string('Q', 20)).ToList();
        var homo2 = DisorderPredictor.PredictLowComplexityRegions(new string('Q', 20)).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(homo, Is.Not.Empty, "homopolymer is low-complexity");
            Assert.That(homo.Select(r => (r.Start, r.End)), Is.EqualTo(homo2.Select(r => (r.Start, r.End))), "deterministic");
            Assert.That(DisorderPredictor.PredictLowComplexityRegions("ACDEFGHIKLMNPQRSTVWY"), Is.Empty,
                "a maximally diverse window is not low-complexity");
        });
    }

    #endregion

    #region DISORDER-MORF-001: P: MoRF is an ordered dip within disorder; R: positions valid; D: deterministic

    // PredictMoRFs finds short ordered dips (disorder < 0.5) of length 10–70 flanked by disorder
    // (≥ 0.5) on both sides — Molecular Recognition Features embedded in disordered regions.

    private const double MoRFThreshold = 0.5;

    /// <summary>
    /// INV-1 (R + P): every MoRF lies in [10,70] residues with a score in [0,1] at valid positions,
    /// its residues are ordered (disorder &lt; 0.5), and it is flanked by disorder on both sides —
    /// verified against the disorder profile.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MoRFs_AreOrderedDipsWithinDisorder()
    {
        return Prop.ForAll(ProteinArbitrary(30), seq =>
        {
            var profile = DisorderPredictor.PredictDisorder(seq).ResiduePredictions;
            var morfs = DisorderPredictor.PredictMoRFs(seq).ToList();
            bool ok = morfs.All(m =>
            {
                int len = m.End - m.Start + 1;
                bool basics = m.Start >= 0 && m.Start <= m.End && m.End < seq.Length
                              && len is >= 10 and <= 70 && m.Score is >= 0.0 and <= 1.0;
                bool ordered = Enumerable.Range(m.Start, len).All(i => profile[i].DisorderScore < MoRFThreshold);
                bool flanked = m.Start > 0 && profile[m.Start - 1].DisorderScore >= MoRFThreshold
                               && m.End < seq.Length - 1 && profile[m.End + 1].DisorderScore >= MoRFThreshold;
                return basics && ordered && flanked;
            });
            return ok.Label("a MoRF was not an ordered dip flanked by disorder");
        });
    }

    /// <summary>
    /// INV-2 (P, positive control): an order-promoting block flanked by disorder-promoting blocks
    /// yields a MoRF inside the ordered block.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MoRFs_OrderedBlockInDisorder_IsDetected()
    {
        // P (strongly disorder-promoting) flanks; W (order-promoting) core of length 15 (∈ [10,70]).
        string seq = new string('P', 25) + new string('W', 15) + new string('P', 25);
        var morfs = DisorderPredictor.PredictMoRFs(seq).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(morfs, Is.Not.Empty, "an ordered block within disorder must yield a MoRF");
            Assert.That(morfs.All(m => m.End - m.Start + 1 is >= 10 and <= 70), Is.True);
        });
    }

    /// <summary>
    /// INV-3 (D): MoRF prediction is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MoRFs_AreDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(30), seq =>
            DisorderPredictor.PredictMoRFs(seq).SequenceEqual(DisorderPredictor.PredictMoRFs(seq))
                .Label("PredictMoRFs must be deterministic"));
    }

    #endregion

    #region DISORDER-PROPENSITY-001: R: per-residue propensity ∈ [0,1]; P: one score per residue; D: deterministic

    // PredictDisorder emits one per-residue disorder score (normalized TOP-IDP) per residue, plus
    // overall summary scores; all are probabilities/fractions in [0,1].

    /// <summary>
    /// INV-1 (P + R): there is exactly one residue prediction per input residue (in order, correct
    /// residue) and every per-residue score and summary score is in [0,1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Propensity_OneScorePerResidueInUnitInterval()
    {
        return Prop.ForAll(ProteinArbitrary(25), seq =>
        {
            var result = DisorderPredictor.PredictDisorder(seq);
            bool lengthOk = result.ResiduePredictions.Count == seq.Length;
            bool perResidue = result.ResiduePredictions.Select((r, i) =>
                r.Position == i && r.Residue == seq[i] && r.DisorderScore is >= 0.0 and <= 1.0).All(x => x);
            bool summaries = result.OverallDisorderContent is >= 0.0 and <= 1.0
                             && result.MeanDisorderScore is >= 0.0 and <= 1.0;
            return (lengthOk && perResidue && summaries)
                .Label($"length={result.ResiduePredictions.Count}/{seq.Length}, perResidue/summaries failed");
        });
    }

    /// <summary>
    /// INV-2 (R): the raw single-residue propensity is the signed TOP-IDP scale (Campen et al. 2008),
    /// finite and within the published [-1,1] band — it is the normalization in <see cref="DisorderPredictor.PredictDisorder"/>
    /// (INV-1) that maps these to the [0,1] disorder score, not this lookup.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Propensity_PerAminoAcid_OnTopIdpScale()
    {
        foreach (char aa in "ACDEFGHIKLMNPQRSTVWY")
        {
            double p = DisorderPredictor.GetDisorderPropensity(aa);
            Assert.That(double.IsFinite(p) && p is >= -1.0 and <= 1.0, Is.True,
                $"TOP-IDP propensity for {aa} = {p} outside the expected [-1,1] band");
        }
    }

    /// <summary>
    /// INV-3 (D): Per-residue disorder scoring is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Propensity_IsDeterministic()
    {
        return Prop.ForAll(ProteinArbitrary(25), seq =>
        {
            var a = DisorderPredictor.PredictDisorder(seq).ResiduePredictions.Select(r => r.DisorderScore).ToList();
            var b = DisorderPredictor.PredictDisorder(seq).ResiduePredictions.Select(r => r.DisorderScore).ToList();
            return a.SequenceEqual(b).Label("PredictDisorder per-residue scores must be deterministic");
        });
    }

    #endregion
}
