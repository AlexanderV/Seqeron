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
}
