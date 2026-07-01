using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.ImmuneAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology immune-infiltration area — ONCO-IMMUNE-001.
/// The unit under test is the ESTIMATE-style enrichment scorer
/// <see cref="ImmuneAnalyzer.EstimateInfiltration"/> (with supporting coverage of
/// <see cref="ImmuneAnalyzer.DeconvoluteImmuneCells"/> for the same degenerate
/// inputs), implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (KeyNotFoundException / NullReference / DivideByZero / Overflow). Every input
/// must resolve to EITHER a well-defined, theory-correct value OR a *documented,
/// intentional* outcome. For immune-signature scoring the headline hazards are:
///   • a DivideByZeroException when the gene set / overlap is empty (the ssGSEA
///     "1/nMiss" miss-step, or a mean over zero genes);
///   • a silent NaN score that corrupts downstream purity when expression values
///     are NaN — the score must stay finite because ssGSEA weights by integer
///     RANK position, never by the raw expression value;
///   • a KeyNotFoundException / NullReference when a signature gene is absent from
///     the expression profile (an "unknown gene");
///   • a log(0) / -inf or negative-mass crash on all-zero or negative expression;
///   • a TumorPurity escaping the documented [0, 1] clamp.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-IMMUNE-001 — Immune infiltration estimation (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 86.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, empty.
///   • MC = Malformed Content — невалідний контент (NaN, ±Infinity).
///     Targets (checklist row 86): "0 expression, all NaN, negative expression,
///     empty gene set, unknown genes".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// ESTIMATE-style scoring (docs/algorithms/Oncology/Immune_Infiltration_Estimation.md):
///   • EstimateScore = ImmuneScore + StromalScore                      (INV-ESTIMATE-01, §2.A)
///   • TumorPurity = clamp(cos(a + b·EstimateScore), 0, 1),
///       a = 0.6049872018, b = 0.0001467884                            (INV-ESTIMATE-02, §2.A)
///   • An empty effective hit set (no overlap, or an empty supplied gene
///       set) yields a ZERO enrichment contribution from that set       (INV-ESTIMATE-03, §6.1)
///   • Empty expression profile → ImmuneScore=StromalScore=Estimate=0,
///       and TumorPurity ≈ cos(0.6049872018) ≈ 0.8225                   (§3.3, §6.1)
///   • Unknown signature genes (absent from the profile) are simply not
///       in the overlap; they are filtered by ContainsKey               (§4.A step 2)
///   • ssGSEA weights hits by rank^τ (τ = 0.25) where rank = N − i — i.e. by
///       integer RANK POSITION, NOT by the raw expression value. So a NaN /
///       Infinity expression value can only perturb the SORT ORDER, never the
///       arithmetic; the score therefore stays finite.                  (§2.A, source ComputeSsGseaScore)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyImmuneFuzzTests
{
    private const double EstimatePurityCoefA = 0.6049872018;
    private const double EstimatePurityCoefB = 0.0001467884;

    // A small set of real ESTIMATE immune-signature genes (present in
    // DefaultImmuneSignatureGenes) used to build deterministic profiles.
    private static readonly string[] ImmuneMarkers =
        { "CD2", "CD3D", "GZMB", "PRF1", "NKG7", "CD27", "CCL5", "CD48" };

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented numeric contract on EVERY result: all four scores must
    // be FINITE (no NaN / Infinity leaking out), TumorPurity must sit inside the
    // documented [0, 1] clamp, EstimateScore must equal the sum of the two
    // component scores, and the overlap counts must be non-negative. This is what
    // stops a fuzz test from rubber-stamping a NaN-corrupted result green.
    private static void AssertWellFormedInfiltration(InfiltrationResult r)
    {
        double.IsNaN(r.ImmuneScore).Should().BeFalse("immune score must never be NaN");
        double.IsNaN(r.StromalScore).Should().BeFalse("stromal score must never be NaN");
        double.IsNaN(r.EstimateScore).Should().BeFalse("ESTIMATE score must never be NaN");
        double.IsNaN(r.TumorPurity).Should().BeFalse("tumor purity must never be NaN");

        double.IsInfinity(r.ImmuneScore).Should().BeFalse("immune score must be finite");
        double.IsInfinity(r.StromalScore).Should().BeFalse("stromal score must be finite");
        double.IsInfinity(r.EstimateScore).Should().BeFalse("ESTIMATE score must be finite");
        double.IsInfinity(r.TumorPurity).Should().BeFalse("tumor purity must be finite");

        // INV-ESTIMATE-01: EstimateScore = ImmuneScore + StromalScore.
        r.EstimateScore.Should().BeApproximately(r.ImmuneScore + r.StromalScore, 1e-9);

        // INV-ESTIMATE-02: documented [0, 1] purity clamp.
        r.TumorPurity.Should().BeGreaterThanOrEqualTo(0.0);
        r.TumorPurity.Should().BeLessThanOrEqualTo(1.0);

        r.OverlappingImmuneGenes.Should().BeGreaterThanOrEqualTo(0);
        r.OverlappingStromalGenes.Should().BeGreaterThanOrEqualTo(0);
    }

    private static double ExpectedPurity(double estimateScore) =>
        Math.Clamp(Math.Cos(EstimatePurityCoefA + EstimatePurityCoefB * estimateScore), 0.0, 1.0);

    #region ONCO-IMMUNE-001 — Immune infiltration estimation

    // ── POSITIVE sanity: documented score + correct ordering ─────────────────

    [Test]
    public void EstimateInfiltration_EmptyProfile_MatchesDocumentedPurityConstant()
    {
        // Docs §3.3 / §6.1 edge case: empty profile → all scores 0, and tumor
        // purity is the published cosine evaluated at score 0 ≈ 0.8225.
        var result = EstimateInfiltration(new Dictionary<string, double>());

        AssertWellFormedInfiltration(result);
        result.ImmuneScore.Should().Be(0.0);
        result.StromalScore.Should().Be(0.0);
        result.EstimateScore.Should().Be(0.0);
        result.TumorPurity.Should().BeApproximately(Math.Cos(EstimatePurityCoefA), 1e-9);
        result.TumorPurity.Should().BeApproximately(0.8225, 1e-3);
        result.OverlappingImmuneGenes.Should().Be(0);
        result.OverlappingStromalGenes.Should().Be(0);
    }

    [Test]
    public void EstimateInfiltration_HandComputedSingleHit_MatchesSsGseaIntegral()
    {
        // Three genes, one of which (CD2) is the ONLY immune signature hit and is
        // ranked at the TOP (highest expression). With N=3, nHits=1, nMiss=2:
        //   the single hit's normalized weight is rank^τ / totalHitWeight = 1, so
        //   the running sum after the hit is +1, then two miss steps of 1/2 each.
        //   Walk (descending: CD2(hit), then two misses):
        //     i=0 hit:  running = +1            integral = 1
        //     i=1 miss: running = +1 − 0.5=0.5  integral = 1.5
        //     i=2 miss: running = 0.5 − 0.5=0   integral = 1.5
        //   → ImmuneScore = 1.5 (hit at the very top). Stromal overlap empty → 0.
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = 100.0,   // immune hit, top rank
            ["AAA"] = 50.0,    // miss
            ["BBB"] = 10.0     // miss
        };

        var result = EstimateInfiltration(
            profile,
            immuneGenes: new[] { "CD2" },
            stromalGenes: Array.Empty<string>());

        AssertWellFormedInfiltration(result);
        result.ImmuneScore.Should().BeApproximately(1.5, 1e-9);
        result.StromalScore.Should().Be(0.0);
        result.EstimateScore.Should().BeApproximately(1.5, 1e-9);
        result.OverlappingImmuneGenes.Should().Be(1);
        result.OverlappingStromalGenes.Should().Be(0);
    }

    [Test]
    public void EstimateInfiltration_HighImmuneOrdersAboveLowImmune()
    {
        // A profile where immune markers dominate the top ranks must score a
        // strictly higher immune enrichment than one where they sit at the bottom.
        var immune = new[] { "CD2", "CD3D", "GZMB", "PRF1" };

        var highProfile = new Dictionary<string, double>
        {
            ["CD2"] = 100, ["CD3D"] = 99, ["GZMB"] = 98, ["PRF1"] = 97,
            ["b1"] = 5, ["b2"] = 4, ["b3"] = 3, ["b4"] = 2
        };
        var lowProfile = new Dictionary<string, double>
        {
            ["t1"] = 100, ["t2"] = 99, ["t3"] = 98, ["t4"] = 97,
            ["CD2"] = 5, ["CD3D"] = 4, ["GZMB"] = 3, ["PRF1"] = 2
        };

        var high = EstimateInfiltration(highProfile, immuneGenes: immune, stromalGenes: Array.Empty<string>());
        var low = EstimateInfiltration(lowProfile, immuneGenes: immune, stromalGenes: Array.Empty<string>());

        AssertWellFormedInfiltration(high);
        AssertWellFormedInfiltration(low);
        high.ImmuneScore.Should().BeGreaterThan(low.ImmuneScore,
            "immune signature at the top of the ranking must enrich more strongly than at the bottom");
    }

    // ── BE: null profile rejected (documented precondition) ──────────────────

    [Test]
    public void EstimateInfiltration_NullProfile_ThrowsArgumentNullException()
    {
        // Docs §3.1 / §3.3: a null expression profile is rejected.
        Action act = () => EstimateInfiltration(null!);
        act.Should().Throw<ArgumentNullException>();
    }

    // ── BE: all-zero expression (no log(0), no DivideByZero) ─────────────────

    [Test]
    public void EstimateInfiltration_AllZeroExpression_ProducesFiniteScores()
    {
        // Every gene has expression 0. ssGSEA ranks by value (all tie) and weights
        // by rank position, so the score is finite — there must be no log(0)/-inf
        // or DivideByZero leaking out. Purity stays in [0, 1].
        var profile = new Dictionary<string, double>();
        foreach (var g in ImmuneMarkers) profile[g] = 0.0;
        for (int i = 0; i < 20; i++) profile["filler" + i] = 0.0;

        Action act = () => EstimateInfiltration(profile);
        act.Should().NotThrow();

        var result = EstimateInfiltration(profile);
        AssertWellFormedInfiltration(result);
        result.TumorPurity.Should().BeApproximately(ExpectedPurity(result.EstimateScore), 1e-9);
    }

    [Test]
    public void EstimateInfiltration_AllZeroExpression_OnlyImmuneGenesPresent_NoMissDivideByZero()
    {
        // EVERY profile gene is also a signature hit → nMiss = N − nHits = 0.
        // The ssGSEA miss-step is 1/nMiss; the source must guard nMiss==0 and
        // return 0 rather than dividing by zero / emitting Infinity.
        var profile = new Dictionary<string, double>();
        var genes = ImmuneMarkers.Take(4).ToArray();
        foreach (var g in genes) profile[g] = 0.0;

        Action act = () => EstimateInfiltration(profile, immuneGenes: genes, stromalGenes: Array.Empty<string>());
        act.Should().NotThrow();

        var result = EstimateInfiltration(profile, immuneGenes: genes, stromalGenes: Array.Empty<string>());
        AssertWellFormedInfiltration(result);
        result.ImmuneScore.Should().Be(0.0, "with no miss genes the ssGSEA walk has no defined denominator → documented 0");
    }

    // ── BE: negative expression (log-transformed inputs) ─────────────────────

    [Test]
    public void EstimateInfiltration_NegativeExpression_ProducesFiniteScores()
    {
        // Log2-transformed expression legitimately goes negative. Ranking still
        // orders them; the score must stay finite and purity inside [0, 1].
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = -2.0, ["CD3D"] = -1.5, ["GZMB"] = -0.5,
            ["DCN"] = -3.0, ["LUM"] = -1.0, ["other"] = -10.0
        };

        Action act = () => EstimateInfiltration(profile);
        act.Should().NotThrow();

        var result = EstimateInfiltration(profile);
        AssertWellFormedInfiltration(result);
    }

    [Test]
    public void EstimateInfiltration_RandomMixedSignExpression_NeverCorruptsResult()
    {
        // Fuzz: random profiles mixing large positive, negative, and zero values
        // across many seeds. Every result must be well-formed.
        for (int seed = 0; seed < 60; seed++)
        {
            var rng = new Random(seed);
            var profile = new Dictionary<string, double>();
            foreach (var g in ImmuneMarkers)
                profile[g] = (rng.NextDouble() - 0.5) * 2_000.0;
            int filler = rng.Next(0, 30);
            for (int i = 0; i < filler; i++)
                profile["g" + i] = (rng.NextDouble() - 0.5) * 2_000.0;

            var result = EstimateInfiltration(profile);
            AssertWellFormedInfiltration(result);
        }
    }

    // ── MC: all-NaN / ±Infinity expression (no NaN propagation) ──────────────

    [Test]
    public void EstimateInfiltration_AllNaNExpression_DoesNotPropagateNaN()
    {
        // ssGSEA weights hits by INTEGER RANK POSITION, not by the raw expression
        // value (source ComputeSsGseaScore: hitWeight = rank^τ / totalHitWeight).
        // Therefore an all-NaN profile can only perturb the (undefined) sort
        // order — it must NOT leak a NaN into the score or corrupt tumor purity.
        var profile = new Dictionary<string, double>();
        foreach (var g in ImmuneMarkers) profile[g] = double.NaN;
        for (int i = 0; i < 15; i++) profile["x" + i] = double.NaN;

        Action act = () => EstimateInfiltration(profile);
        act.Should().NotThrow();

        var result = EstimateInfiltration(profile);
        AssertWellFormedInfiltration(result);
    }

    [Test]
    public void EstimateInfiltration_InfinityExpression_DoesNotPropagateInfinity()
    {
        // ±Infinity values, like NaN, only affect ordering — never the rank-based
        // arithmetic. Result must remain finite and purity in [0, 1].
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = double.PositiveInfinity,
            ["CD3D"] = double.NegativeInfinity,
            ["GZMB"] = double.PositiveInfinity,
            ["filler1"] = double.NegativeInfinity,
            ["filler2"] = 1.0
        };

        Action act = () => EstimateInfiltration(profile);
        act.Should().NotThrow();

        var result = EstimateInfiltration(profile);
        AssertWellFormedInfiltration(result);
    }

    [Test]
    public void EstimateInfiltration_MixedNaNAndFiniteExpression_NeverThrowsOrCorrupts()
    {
        // Fuzz: each gene is randomly NaN, +Inf, -Inf, or a finite value.
        for (int seed = 0; seed < 50; seed++)
        {
            var rng = new Random(seed);
            var profile = new Dictionary<string, double>();
            var names = ImmuneMarkers.Concat(new[] { "DCN", "LUM", "FAP", "fillerA", "fillerB", "fillerC" });
            foreach (var g in names)
            {
                profile[g] = rng.Next(0, 4) switch
                {
                    0 => double.NaN,
                    1 => double.PositiveInfinity,
                    2 => double.NegativeInfinity,
                    _ => (rng.NextDouble() - 0.5) * 100.0
                };
            }

            var result = EstimateInfiltration(profile);
            AssertWellFormedInfiltration(result);
        }
    }

    // ── BE: empty gene set (no DivideByZero on zero-length walk) ──────────────

    [Test]
    public void EstimateInfiltration_EmptyImmuneAndStromalGeneSets_YieldZeroScores()
    {
        // INV-ESTIMATE-03 / §6.1: an empty (or non-overlapping) gene set yields a
        // zero contribution — the ssGSEA helper returns 0 when the hit set is
        // empty, with no DivideByZero on a zero-length mean.
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = 5.0, ["CD3D"] = 4.0, ["DCN"] = 3.0, ["other"] = 1.0
        };

        var result = EstimateInfiltration(
            profile,
            immuneGenes: Array.Empty<string>(),
            stromalGenes: Array.Empty<string>());

        AssertWellFormedInfiltration(result);
        result.ImmuneScore.Should().Be(0.0);
        result.StromalScore.Should().Be(0.0);
        result.EstimateScore.Should().Be(0.0);
        result.OverlappingImmuneGenes.Should().Be(0);
        result.OverlappingStromalGenes.Should().Be(0);
        // Empty sets ⇒ score 0 ⇒ same purity constant as the empty-profile case.
        result.TumorPurity.Should().BeApproximately(Math.Cos(EstimatePurityCoefA), 1e-9);
    }

    // ── BE/MC: unknown genes (no KeyNotFoundException) ───────────────────────

    [Test]
    public void EstimateInfiltration_AllUnknownSignatureGenes_NoOverlapZeroScore()
    {
        // None of the supplied signature genes exist in the profile → overlap is
        // empty. Must NOT throw KeyNotFoundException; score contribution is 0.
        var profile = new Dictionary<string, double>
        {
            ["realGeneA"] = 9.0, ["realGeneB"] = 8.0, ["realGeneC"] = 7.0
        };
        var unknown = new[] { "NOPE1", "NOPE2", "DOESNOTEXIST", "ZZZ" };

        Action act = () => EstimateInfiltration(profile, immuneGenes: unknown, stromalGenes: unknown);
        act.Should().NotThrow();

        var result = EstimateInfiltration(profile, immuneGenes: unknown, stromalGenes: unknown);
        AssertWellFormedInfiltration(result);
        result.OverlappingImmuneGenes.Should().Be(0);
        result.OverlappingStromalGenes.Should().Be(0);
        result.ImmuneScore.Should().Be(0.0);
        result.StromalScore.Should().Be(0.0);
    }

    [Test]
    public void EstimateInfiltration_UnknownGenesMixedWithKnown_CountsOnlyOverlap()
    {
        // A signature mixing real markers with garbage names must count ONLY the
        // overlapping (real) genes and never dereference a missing key.
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = 10.0, ["CD3D"] = 9.0, ["GZMB"] = 8.0,
            ["bg1"] = 2.0, ["bg2"] = 1.0
        };
        var immune = new[] { "CD2", "CD3D", "GZMB", "GHOST1", "GHOST2", "GHOST3" };

        var result = EstimateInfiltration(profile, immune, stromalGenes: Array.Empty<string>());

        AssertWellFormedInfiltration(result);
        result.OverlappingImmuneGenes.Should().Be(3, "only CD2/CD3D/GZMB are present in the profile");
    }

    // ── BE: single-gene profile (degenerate ranking) ─────────────────────────

    [Test]
    public void EstimateInfiltration_SingleGeneProfile_NoCrash()
    {
        // N = 1: if that gene is a hit then nMiss = 0; if not, nHits = 0. Either
        // way the ssGSEA helper must return a documented 0, not divide by zero.
        var hitOnly = new Dictionary<string, double> { ["CD2"] = 1.0 };
        var missOnly = new Dictionary<string, double> { ["NOTASIG"] = 1.0 };

        var r1 = EstimateInfiltration(hitOnly, immuneGenes: new[] { "CD2" }, stromalGenes: Array.Empty<string>());
        var r2 = EstimateInfiltration(missOnly, immuneGenes: new[] { "CD2" }, stromalGenes: Array.Empty<string>());

        AssertWellFormedInfiltration(r1);
        AssertWellFormedInfiltration(r2);
        r1.ImmuneScore.Should().Be(0.0);
        r2.ImmuneScore.Should().Be(0.0);
    }

    #endregion

    #region ONCO-IMMUNE-001 — Deconvolution under the same degenerate inputs

    // The deconvolution entry point shares the immune-infiltration contract's
    // degenerate-input surface (0 / NaN / negative / empty / unknown). It must be
    // equally robust: no DivideByZero on empty overlap, no NaN fractions, no
    // KeyNotFoundException on unknown genes, fractions non-negative & summing to 1
    // (or all-zero on the no-overlap branch).  — docs §3.3, §6.1, INV-NNLS-01/02.

    private static void AssertWellFormedDeconvolution(DeconvolutionResult r)
    {
        r.CellFractions.Should().NotBeNull();
        foreach (var f in r.CellFractions.Values)
        {
            double.IsNaN(f).Should().BeFalse("cell fractions must never be NaN");
            double.IsInfinity(f).Should().BeFalse("cell fractions must be finite");
            f.Should().BeGreaterThanOrEqualTo(0.0, "INV-NNLS-01: fractions are non-negative");
        }
        double.IsNaN(r.Correlation).Should().BeFalse();
        double.IsNaN(r.Rmse).Should().BeFalse();
        r.Rmse.Should().BeGreaterThanOrEqualTo(0.0);
        r.OverlappingGenes.Should().BeGreaterThanOrEqualTo(0);

        double sum = r.CellFractions.Values.Sum();
        // INV-NNLS-02: fractions sum to 1 when there is positive mass, else 0.
        (Math.Abs(sum - 1.0) < 1e-6 || Math.Abs(sum) < 1e-9).Should().BeTrue(
            $"fraction sum must be ~1 (positive mass) or ~0 (no-overlap branch); was {sum}");
    }

    [Test]
    public void DeconvoluteImmuneCells_NullProfile_Throws()
    {
        Action act = () => DeconvoluteImmuneCells(null!);
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void DeconvoluteImmuneCells_AllUnknownGenes_ZeroFractionsNoThrow()
    {
        // No overlap with the signature matrix → all-zero fractions, zero metrics,
        // zero overlap (no-overlap branch). Must not throw KeyNotFoundException.
        var profile = new Dictionary<string, double> { ["NOPE1"] = 9.0, ["NOPE2"] = 8.0 };

        Action act = () => DeconvoluteImmuneCells(profile);
        act.Should().NotThrow();

        var result = DeconvoluteImmuneCells(profile);
        AssertWellFormedDeconvolution(result);
        result.OverlappingGenes.Should().Be(0);
        result.CellFractions.Values.Should().OnlyContain(f => f == 0.0);
        result.Correlation.Should().Be(0.0);
        result.Rmse.Should().Be(0.0);
    }

    [Test]
    public void DeconvoluteImmuneCells_AllZeroExpression_NoCrash()
    {
        // Every overlapping signature gene has expression 0 → mixture vector is
        // all-zero. NNLS must not divide by zero; fractions stay well-formed.
        var profile = new Dictionary<string, double>
        {
            ["CD8A"] = 0.0, ["CD8B"] = 0.0, ["CD3D"] = 0.0, ["GZMB"] = 0.0, ["PRF1"] = 0.0
        };

        Action act = () => DeconvoluteImmuneCells(profile);
        act.Should().NotThrow();
        AssertWellFormedDeconvolution(DeconvoluteImmuneCells(profile));
    }

    [Test]
    public void DeconvoluteImmuneCells_NegativeExpression_NoCrash()
    {
        var profile = new Dictionary<string, double>
        {
            ["CD8A"] = -1.0, ["CD8B"] = -0.5, ["CD3D"] = -2.0, ["GZMB"] = 1.0, ["PRF1"] = 0.5
        };

        Action act = () => DeconvoluteImmuneCells(profile);
        act.Should().NotThrow();
        AssertWellFormedDeconvolution(DeconvoluteImmuneCells(profile));
    }

    [Test]
    public void DeconvoluteImmuneCells_EmptyProfile_ZeroFractions()
    {
        var result = DeconvoluteImmuneCells(new Dictionary<string, double>());
        AssertWellFormedDeconvolution(result);
        result.OverlappingGenes.Should().Be(0);
        result.CellFractions.Values.Should().OnlyContain(f => f == 0.0);
    }

    #endregion
}
