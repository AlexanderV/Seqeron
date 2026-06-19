using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the ProteinPred area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: DISORDER-PRED-001 — intrinsic disorder prediction (ProteinPred).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 80.
///
/// API under test (DisorderPredictor.PredictDisorder):
///   Assigns each residue a disorder score = the average normalised TOP-IDP propensity over a
///   sliding window of size w (half-window h = w/2). Disorder-promoting residues (e.g. P, E, K)
///   have high propensity; order-promoting residues (e.g. W, F, I) have low propensity.
///
/// Relations (derived from the windowed propensity average, NOT from output):
///   • INV  (same seq ⇒ same scores): prediction is deterministic and case-insensitive.
///   • MON  (more disorder-promoting residues ⇒ higher disorder): replacing order-promoting
///          residues with disorder-promoting ones (proline) raises the mean disorder score.
///   • INV  (single residue change ⇒ local effect): a point substitution only changes residue
///          scores whose window covers it (|i − p| ≤ h); residues farther away are unchanged.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ProteinPredMetamorphicTests
{
    private static double[] ResidueScores(string sequence) =>
        DisorderPredictor.PredictDisorder(sequence).ResiduePredictions.Select(r => r.DisorderScore).ToArray();

    #region DISORDER-PRED-001 INV — prediction is deterministic and case-insensitive

    [Test]
    [Description("INV: PredictDisorder is a pure, case-insensitive function — repeated calls (and upper/lower case) give identical scores.")]
    public void PredictDisorder_SameSequence_SameScores()
    {
        const string seq = "MKKLLPESTQRGSAADEWFYILVNPPPQEKR";

        var first = DisorderPredictor.PredictDisorder(seq);
        var again = DisorderPredictor.PredictDisorder(seq);
        var lower = DisorderPredictor.PredictDisorder(seq.ToLowerInvariant());

        again.MeanDisorderScore.Should().Be(first.MeanDisorderScore, because: "prediction has no hidden state");
        ResidueScores(seq.ToLowerInvariant()).Should().Equal(ResidueScores(seq),
            because: "the sequence is upper-cased before scoring, so case does not matter");
    }

    #endregion

    #region DISORDER-PRED-001 MON — more disorder-promoting residues raise the disorder score

    [Test]
    [Description("MON: replacing order-promoting residues (W) with the disorder-promoting proline raises the mean disorder score monotonically.")]
    public void PredictDisorder_MoreProline_HigherDisorder()
    {
        const char order = 'W'; // order-promoting (low TOP-IDP propensity)
        const char proline = 'P'; // disorder-promoting

        DisorderPredictor.GetDisorderPropensity(proline)
            .Should().BeGreaterThan(DisorderPredictor.GetDisorderPropensity(order),
                because: "proline has a higher TOP-IDP disorder propensity than tryptophan");

        const int n = 30;
        double previous = double.MinValue;
        foreach (int prolineCount in new[] { 0, 10, 20, 30 })
        {
            string seq = new string(proline, prolineCount) + new string(order, n - prolineCount);
            double mean = DisorderPredictor.PredictDisorder(seq).MeanDisorderScore;

            mean.Should().BeGreaterThan(previous,
                because: $"replacing order-promoting residues with {prolineCount} prolines raises the average disorder propensity");
            previous = mean;
        }
    }

    #endregion

    #region DISORDER-PRED-001 INV — a point substitution has only a local effect

    [Test]
    [Description("INV: changing a single residue only affects scores whose window covers it (|i − p| ≤ h); residues farther than the half-window are unchanged.")]
    public void PredictDisorder_SingleResidueChange_LocalEffectOnly()
    {
        const int windowSize = 21;
        const int halfWindow = windowSize / 2; // 10
        const int p = 25;

        string baseSeq = new string('A', 50);
        char[] mutatedChars = baseSeq.ToCharArray();
        mutatedChars[p] = 'E'; // swap one residue for a disorder-promoting glutamate
        string mutatedSeq = new string(mutatedChars);

        var baseScores = ResidueScores(baseSeq);
        var mutatedScores = ResidueScores(mutatedSeq);

        for (int i = 0; i < baseSeq.Length; i++)
        {
            if (System.Math.Abs(i - p) > halfWindow)
                mutatedScores[i].Should().Be(baseScores[i],
                    because: $"residue {i}'s window does not cover the mutated position {p}, so its score is unchanged");
        }

        // Non-vacuous: at least one residue whose window covers p actually changed.
        Enumerable.Range(0, baseSeq.Length)
            .Any(i => System.Math.Abs(i - p) <= halfWindow && mutatedScores[i] != baseScores[i])
            .Should().BeTrue(because: "residues whose window covers the mutation must feel its local effect");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: DISORDER-REGION-001 — disordered-region calling (ProteinPred).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 81.
    //
    // API under test (DisorderPredictor.PredictDisorder → DisorderedRegions):
    //   A residue is "disordered" when its window score ≥ threshold; disordered regions are the
    //   maximal runs of disordered residues of length ≥ minRegionLength.
    //
    // Relations (derived from the score ≥ threshold predicate, NOT from output):
    //   • MON  (lower threshold ⇒ more/larger regions): the predicate is monotone in the
    //          threshold, so lowering it only adds disordered residues — region coverage grows.
    //   • SUB  (strict ⊆ lenient): the residues covered at a stricter threshold are a subset of
    //          those covered at a more lenient one.
    //   • INV  (distant ordered insert ⇒ unaffected): appending an ordered block at the far 3'
    //          end leaves the scores (and regions) of residues beyond its window reach unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region DISORDER-REGION-001 — Helpers

    // Ordered flanks (W) around a disorder-promoting core (E).
    private const string DisorderTestSeq = "WWWWWWWWWWWWWWW" + "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" + "WWWWWWWWWWWWWWW";

    private static System.Collections.Generic.HashSet<int> CoveredResidues(string sequence, double threshold)
    {
        var covered = new System.Collections.Generic.HashSet<int>();
        foreach (var region in DisorderPredictor.PredictDisorder(sequence, disorderThreshold: threshold).DisorderedRegions)
            for (int i = region.Start; i <= region.End; i++)
                covered.Add(i);
        return covered;
    }

    #endregion

    #region DISORDER-REGION-001 MON / SUB — lowering the threshold grows region coverage

    [Test]
    [Description("MON/SUB: score ≥ threshold is monotone in the threshold, so lowering it only adds covered residues — coverage at a stricter threshold is a subset of coverage at a more lenient one.")]
    public void DisorderedRegions_LowerThreshold_SupersetCoverage()
    {
        double[] descending = { 0.60, 0.50, 0.40, 0.30, 0.20 };

        System.Collections.Generic.HashSet<int>? previous = null;
        foreach (double threshold in descending)
        {
            var covered = CoveredResidues(DisorderTestSeq, threshold);

            if (previous is not null)
                covered.IsSupersetOf(previous).Should().BeTrue(
                    because: $"lowering the threshold to {threshold} can only add disordered residues, never remove them");
            previous = covered;
        }

        CoveredResidues(DisorderTestSeq, 0.20).Count
            .Should().BeGreaterThan(CoveredResidues(DisorderTestSeq, 0.60).Count,
                because: "a lenient threshold calls strictly more disordered residues than a strict one");
    }

    #endregion

    #region DISORDER-REGION-001 INV — a distant ordered insert doesn't affect the disordered region

    [Test]
    [Description("INV: appending an ordered block at the far 3' end leaves the scores (and the disordered region) of residues beyond its window reach unchanged.")]
    public void DisorderedRegions_DistantOrderedInsert_DoesNotAffectDistantDisorder()
    {
        const int halfWindow = 21 / 2; // 10
        int n = DisorderTestSeq.Length; // 60

        // Append a 40-residue ordered (W) block; it can only influence windows of residues ≥ n−h.
        string extended = DisorderTestSeq + new string('W', 40);

        var baseScores = ResidueScores(DisorderTestSeq);
        var extendedScores = ResidueScores(extended);

        for (int i = 0; i < n - halfWindow; i++)
            extendedScores[i].Should().Be(baseScores[i],
                because: $"residue {i} is more than the half-window from the appended block, so its score is unchanged");

        // The disordered E-core (indices 15..44) lies well within that unaffected prefix.
        CoveredResidues(extended, 0.40).Where(i => i < n - halfWindow)
            .Should().BeEquivalentTo(CoveredResidues(DisorderTestSeq, 0.40).Where(i => i < n - halfWindow),
                because: "the disordered region far from the insertion is preserved exactly");
    }

    #endregion
}
