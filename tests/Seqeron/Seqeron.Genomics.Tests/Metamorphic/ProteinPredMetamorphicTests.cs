using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Metamorphic;

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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: DISORDER-LC-001 — low-complexity region prediction (SEG; ProteinPred).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 204.
    //
    // API under test (DisorderPredictor.PredictLowComplexityRegions):
    //   SEG algorithm (Wootton & Federhen 1993; NCBI blast_seg.c). Complexity is the Shannon
    //   entropy H = −Σ pᵢ·log₂(pᵢ) (bits/residue) of a window's residue composition. A window is
    //   "low complexity" when its entropy is LOW. Two-stage scan: stage 1 marks windows with
    //   H ≤ K1 (trigger); stage 2 extends triggered spans while the whole-segment H ≤ K2; a final
    //   length filter keeps segments of length ≥ minLength.
    //
    // Relations (derived from the H ≤ K low-complexity predicate, NOT from output):
    //   • MON  (more permissive complexity cutoff ⇒ superset): the predicate H ≤ K is monotone in
    //          the entropy cutoff K, so raising K1/K2 only adds triggered/extended positions —
    //          low-complexity coverage grows. NOTE the checklist shorthand "lower threshold →
    //          superset" is the score-predicate (score ≥ τ) convention; for SEG the criterion is an
    //          entropy CEILING (H ≤ K), so the permissive direction is a HIGHER cutoff. Same
    //          monotonicity, inverted sign — the test encodes SEG's actual H ≤ K theory.
    //   • MON  (lower minLength ⇒ superset): minLength filters only completed segments by length,
    //          so lowering it admits a superset of the identical region objects (the literal
    //          "lower threshold → superset").
    //   • SHIFT (prepend flank shifts regions): prepending a high-complexity (non-triggering,
    //          extension-stopping) flank translates every reported low-complexity segment's
    //          coordinates by exactly the flank length and preserves its composition type.
    // ───────────────────────────────────────────────────────────────────────────

    #region DISORDER-LC-001 — Helpers

    // A maximally diverse 20-mer (each of the 20 standard amino acids once): every 12-residue
    // window has entropy ≈ log2(12) ≈ 3.585 bits ≫ the default K1/K2 (2.2 / 2.5), so the flank
    // neither triggers a low-complexity call nor permits a triggered neighbour to extend into it.
    private const string DiverseFlank = "ACDEFGHIKLMNPQRSTVWY";

    private static System.Collections.Generic.HashSet<int> LowComplexityCoverage(
        string sequence, double triggerThreshold, double extensionThreshold)
    {
        var covered = new System.Collections.Generic.HashSet<int>();
        foreach (var (start, end, _) in DisorderPredictor.PredictLowComplexityRegions(
                     sequence, triggerThreshold: triggerThreshold, extensionThreshold: extensionThreshold))
            for (int i = start; i <= end; i++)
                covered.Add(i);
        return covered;
    }

    #endregion

    #region DISORDER-LC-001 MON — a more permissive complexity cutoff grows low-complexity coverage

    [Test]
    [Description("MON: H ≤ K is monotone in the entropy cutoff K, so raising the SEG trigger/extension thresholds only adds low-complexity positions — coverage at a stricter cutoff is a subset of coverage at a more permissive one.")]
    public void LowComplexity_HigherEntropyCutoff_SupersetCoverage()
    {
        // A homopolymer run (H = 0, always low-complexity) and a 6-letter periodic block
        // ("ACDEFG"×4, every 12-window H = log2(6) ≈ 2.585) separated by diverse flanks so they
        // form distinct candidates. The periodic block is low-complexity only once K1 ≥ 2.585.
        const string seq =
            DiverseFlank + "AAAAAAAAAAAA" + DiverseFlank +
            "ACDEFGACDEFGACDEFGACDEFG" + DiverseFlank;

        // Ascending entropy cutoffs (K1, K2): each pair is at least as permissive as the previous.
        var cutoffs = new[] { (K1: 2.2, K2: 2.5), (K1: 2.7, K2: 3.0), (K1: 3.2, K2: 3.4) };

        System.Collections.Generic.HashSet<int>? previous = null;
        foreach (var (k1, k2) in cutoffs)
        {
            var covered = LowComplexityCoverage(seq, k1, k2);
            if (previous is not null)
                covered.IsSupersetOf(previous).Should().BeTrue(
                    because: $"raising the entropy cutoff to (K1={k1}, K2={k2}) can only add low-complexity positions, never remove them");
            previous = covered;
        }

        LowComplexityCoverage(seq, 3.2, 3.4).Count
            .Should().BeGreaterThan(LowComplexityCoverage(seq, 2.2, 2.5).Count,
                because: "a permissive cutoff additionally flags the periodic block (H ≈ 2.585) that the strict default rejects");
    }

    #endregion

    #region DISORDER-LC-001 MON — lowering minLength admits a superset of regions

    [Test]
    [Description("MON: minLength filters only completed segments by length, so lowering it yields a superset of the identical low-complexity region objects.")]
    public void LowComplexity_LowerMinLength_SupersetOfRegions()
    {
        // A short and a longer homopolymer run, separated by a wide diverse spacer (3 flanks).
        // A long homopolymer's SEG extension is greedy (the dominant residue keeps whole-segment
        // entropy low), so the spacer must exceed that reach to stop the two runs from merging.
        // After extension the short run stays below 40 residues and the longer run above it.
        const string spacer = DiverseFlank + DiverseFlank + DiverseFlank;
        const string seq =
            DiverseFlank + "AAAAAAAAAAAA" + spacer +
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + DiverseFlank;

        var lenient = DisorderPredictor.PredictLowComplexityRegions(seq, minLength: 1).ToList();
        var strict = DisorderPredictor.PredictLowComplexityRegions(seq, minLength: 40).ToList();

        lenient.Should().HaveCountGreaterThan(1, because: "both homopolymer runs are reported when no length filter applies");
        strict.Should().NotBeEmpty(because: "the long run survives the length filter");

        strict.Should().BeSubsetOf(lenient,
            because: "minLength removes whole segments without altering their boundaries, so the strict set is a subset of the lenient set");
        strict.Should().HaveCountLessThan(lenient.Count,
            because: "the short run is filtered out at minLength = 40 but present at minLength = 1");
    }

    #endregion

    #region DISORDER-LC-001 SHIFT — prepending a diverse flank shifts low-complexity regions

    [Test]
    [Description("SHIFT: prepending a high-complexity (non-triggering, extension-stopping) flank translates every low-complexity segment's coordinates by exactly the flank length and preserves its type.")]
    public void LowComplexity_PrependDiverseFlank_ShiftsRegions()
    {
        // A low-complexity core embedded between diverse flanks; the flanks halt SEG extension well
        // inside themselves, so the reported segment lies strictly interior to the sequence.
        const string core = "AAAAAAAAAAAAAAAA";
        string seq = DiverseFlank + core + DiverseFlank;

        var baseline = DisorderPredictor.PredictLowComplexityRegions(seq).ToList();
        baseline.Should().ContainSingle(because: "the single homopolymer core is the only low-complexity region");
        baseline[0].Start.Should().BeGreaterThan(0, because: "extension stops inside the left flank, not at the sequence start");
        baseline[0].End.Should().BeLessThan(seq.Length - 1, because: "extension stops inside the right flank, not at the sequence end");

        foreach (int flankCount in new[] { 1, 3 })
        {
            string prefix = string.Concat(Enumerable.Repeat(DiverseFlank, flankCount));
            int offset = prefix.Length;
            var shifted = DisorderPredictor.PredictLowComplexityRegions(prefix + seq).ToList();

            shifted.Should().HaveCount(baseline.Count,
                because: "a non-interacting flank changes neither which windows trigger nor how far they extend");

            for (int i = 0; i < baseline.Count; i++)
            {
                shifted[i].Start.Should().Be(baseline[i].Start + offset,
                    because: $"prepending {offset} diverse residues shifts the segment start by {offset}");
                shifted[i].End.Should().Be(baseline[i].End + offset,
                    because: $"prepending {offset} diverse residues shifts the segment end by {offset}");
                shifted[i].Type.Should().Be(baseline[i].Type,
                    because: "the segment's composition (and thus its low-complexity type) is unchanged by translation");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: DISORDER-MORF-001 — Molecular Recognition Feature (MoRF) prediction (ProteinPred).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 205.
    //
    // API under test (DisorderPredictor.PredictMoRFs):
    //   A MoRF is a short region of relative ORDER — a downward "dip" in the per-residue disorder
    //   profile (TOP-IDP score < 0.5) — embedded WITHIN disorder (flanked by predicted-disordered
    //   residues, score ≥ 0.5, on both immediate sides) whose length lies in the 10–70 band
    //   (Cheng/Oldfield et al. PMC2570644; Mohan et al. 2006).
    //
    // Relations (derived from the dip-in-disorder definition, NOT from output):
    //   • INV  (deterministic): PredictMoRFs is a pure, case-insensitive function — repeated calls
    //          and upper/lower case give identical MoRFs.
    //   • SHIFT (prepend flank shifts MoRFs): prepending more of the same disordered residue extends
    //          the disordered context without altering the local profile around the ordered dip, so
    //          every MoRF's coordinates translate by exactly the prepended length.
    // ───────────────────────────────────────────────────────────────────────────

    #region DISORDER-MORF-001 — Helpers

    // Disordered context residue: all-proline scores 1.0 on the normalised TOP-IDP scale (P is the
    // most disorder-promoting residue), comfortably ≥ the 0.5 MoRF order/disorder boundary.
    private const char MorfDisorderedUnit = 'P';
    // Ordered "dip" residue: all-tryptophan scores 0.0 (W is the most order-promoting residue), so a
    // long W block is a deep dip below 0.5 — the MoRF signature.
    private const char MorfOrderedUnit = 'W';

    // A 40-residue ordered dip embedded in 30-residue disordered flanks: PredictDisorder's window-21
    // profile makes the central W block a single MoRF flanked by disorder on both sides.
    private static string MorfSequence(int leftFlank, int dip, int rightFlank) =>
        new string(MorfDisorderedUnit, leftFlank) + new string(MorfOrderedUnit, dip) + new string(MorfDisorderedUnit, rightFlank);

    private static System.Collections.Generic.List<(int Start, int End, double Score)> MoRFs(string sequence) =>
        DisorderPredictor.PredictMoRFs(sequence).ToList();

    #endregion

    #region DISORDER-MORF-001 INV — MoRF prediction is deterministic and case-insensitive

    [Test]
    [Description("INV: PredictMoRFs is a pure, case-insensitive function — repeated calls and upper/lower case give identical MoRFs.")]
    public void PredictMoRFs_SameSequence_SameMoRFs()
    {
        string seq = MorfSequence(30, 40, 30);

        var first = MoRFs(seq);
        first.Should().ContainSingle(because: "the embedded ordered dip is a single MoRF — a non-vacuous fixture");

        MoRFs(seq).Should().Equal(first, because: "MoRF prediction has no hidden state");
        MoRFs(seq.ToLowerInvariant()).Should().Equal(first, because: "the sequence is upper-cased before scoring, so case does not matter");
    }

    #endregion

    #region DISORDER-MORF-001 SHIFT — prepending disordered residues shifts the MoRFs

    [Test]
    [Description("SHIFT: prepending more of the same disordered residue extends the disordered context without changing the profile around the ordered dip, so every MoRF's coordinates shift by exactly the prepended length.")]
    public void PredictMoRFs_PrependDisorderedFlank_ShiftsMoRFs()
    {
        // Left flank ≥ the window half-width (10) so the dip's local profile sees only flank residues
        // near its boundary — prepending further disordered residues cannot reach it.
        string seq = MorfSequence(30, 40, 30);

        var baseline = MoRFs(seq);
        baseline.Should().ContainSingle(because: "the central dip is the only MoRF");
        baseline[0].Start.Should().BeGreaterThan(0, because: "the MoRF is flanked by disorder, so it never touches the N-terminus");

        foreach (int offset in new[] { 15, 25 })
        {
            string prepended = new string(MorfDisorderedUnit, offset) + seq;
            var shifted = MoRFs(prepended);

            shifted.Should().HaveCount(baseline.Count,
                because: "extending the disordered context neither creates nor destroys an ordered dip");

            for (int i = 0; i < baseline.Count; i++)
            {
                shifted[i].Start.Should().Be(baseline[i].Start + offset,
                    because: $"prepending {offset} disordered residues shifts the MoRF start by {offset}");
                shifted[i].End.Should().Be(baseline[i].End + offset,
                    because: $"prepending {offset} disordered residues shifts the MoRF end by {offset}");
                shifted[i].Score.Should().BeApproximately(baseline[i].Score, 1e-9,
                    because: "the dip's residue composition (and thus its depth score) is unchanged by translation");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: DISORDER-PROPENSITY-001 — per-residue disorder-propensity profile (ProteinPred).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 206.
    //
    // API under test (DisorderPredictor.PredictDisorder → ResiduePredictions[].DisorderScore):
    //   The propensity profile is the array of per-residue scores, each the average normalised
    //   TOP-IDP propensity over a sliding window of size w (half-window h = w/2). Position i's score
    //   depends only on the window [i−h, i+h].
    //
    // Relations (derived from the windowed-average definition, NOT from output):
    //   • INV   (deterministic): the profile is a pure, case-insensitive function of the sequence.
    //   • SHIFT (prepend flank shifts profile): because score i sees only [i−h, i+h], prepending a
    //           flank of length L leaves every interior window (original i ≥ h) intact but at
    //           coordinate i+L — so profile_ext[i+L] = profile_orig[i] exactly, independent of the
    //           flank's content. This is the positional-translation covariance of the profile.
    // ───────────────────────────────────────────────────────────────────────────

    #region DISORDER-PROPENSITY-001 — Helpers

    // A compositionally varied protein region (ordered W/F/I/V and disordered P/E/K/S stretches), so
    // the propensity profile is genuinely non-constant — a meaningful, non-vacuous fixture.
    private const string PropensityTestSeq = "MKKLLPESTQRGSAADEWFYILVNPPPQEKRGGSSEEDDKKRRWWFFIIVVNNCCTT";

    private static double[] PropensityProfile(string sequence, int windowSize) =>
        DisorderPredictor.PredictDisorder(sequence, windowSize: windowSize).ResiduePredictions
            .Select(r => r.DisorderScore).ToArray();

    #endregion

    #region DISORDER-PROPENSITY-001 INV — the propensity profile is deterministic and case-insensitive

    [Test]
    [Description("INV: the per-residue propensity profile is a pure, case-insensitive function — repeated calls and upper/lower case give identical profiles.")]
    public void PropensityProfile_SameSequence_SameProfile()
    {
        const int windowSize = 21;

        var profile = PropensityProfile(PropensityTestSeq, windowSize);
        profile.Distinct().Should().HaveCountGreaterThan(1, because: "the test sequence yields a genuinely varying profile — a non-vacuous fixture");

        PropensityProfile(PropensityTestSeq, windowSize).Should().Equal(profile, because: "the profile computation has no hidden state");
        PropensityProfile(PropensityTestSeq.ToLowerInvariant(), windowSize).Should().Equal(profile,
            because: "the sequence is upper-cased before scoring, so case does not matter");
    }

    #endregion

    #region DISORDER-PROPENSITY-001 SHIFT — prepending a flank translates the interior profile

    [Test]
    [Description("SHIFT: score i depends only on the window [i−h, i+h], so prepending a flank of length L translates every interior score (original i ≥ h) to coordinate i+L unchanged, regardless of the flank's content.")]
    public void PropensityProfile_PrependFlank_ShiftsInteriorProfile()
    {
        foreach (int windowSize in new[] { 11, 21 })
        {
            int halfWindow = windowSize / 2;
            var original = PropensityProfile(PropensityTestSeq, windowSize);

            // The flank content is deliberately arbitrary: the relation holds for ANY prepended
            // residues because interior windows never reach into the flank.
            foreach (string flank in new[] { "GGGGGGG", "ACDEFGHIKLMNP" })
            {
                int offset = flank.Length;
                var extended = PropensityProfile(flank + PropensityTestSeq, windowSize);

                extended.Should().HaveCount(original.Length + offset, because: "the profile has one score per residue");

                for (int i = halfWindow; i < PropensityTestSeq.Length; i++)
                    extended[i + offset].Should().Be(original[i],
                        because: $"residue {i}'s window [{i - halfWindow},{i + halfWindow}] is clear of the {offset}-residue flank, so its score only translates by {offset}");
            }
        }
    }

    #endregion
}
