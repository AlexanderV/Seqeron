using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.ImmuneAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology CIBERSORT-style ν-support-vector-regression immune-cell
/// deconvolution — IMMUNE-NUSVR-001 (the ν-SVR half of ONCO-IMMUNE-001).
/// The unit under test is
/// <see cref="ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(IReadOnlyDictionary{string,double}, IReadOnlyDictionary{string,IReadOnlyDictionary{string,double}}?, IReadOnlyList{double}?)"/>
/// in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/ImmuneAnalyzer.cs (~line 658), with its
/// SMO solver <c>SolveNuSvrLinear</c> (~line 1082), the <c>Standardize</c> / <c>StandardizeColumns</c>
/// z-score helpers (~line 991 / 1023), and the <see cref="ImmuneAnalyzer.NuSvrDeconvolutionResult"/>
/// record (~line 371).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed mixtures and signature matrices and asserts the
/// ν-SVR deconvolution NEVER fails in an undisciplined way: it must not hang (the SMO solver has a
/// <c>maxIterations = 200·ℓ</c> cap and must always terminate — pinned with [CancelAfter]); it must not
/// throw an unhandled runtime exception (DivideByZero in the std-=0 z-score branch, KeyNotFound on a
/// gene present in the signature but not the profile, NullReference); and it must never emit
/// out-of-contract output — a cell fraction &lt; 0 or &gt; 1, fractions not summing to 1 when there is
/// positive post-clip mass, or a NaN / Infinity leaking into a contracted-finite field. Every input
/// resolves to EITHER a well-defined, THEORY-correct result OR a documented validation exception
/// (a null profile → <see cref="ArgumentNullException"/>, §3.1 / §3.3).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzz strategies (docs/checklists/03_FUZZING.md §Description), row 247
/// ───────────────────────────────────────────────────────────────────────────
/// MC = Malformed Content, BE = Boundary Exploitation. Checklist targets for this row:
///   • "gene count mismatch" (BE) — the signature matrix and profile share only a SUBSET of genes, or
///     NO genes; partial overlap fits on the intersection, zero overlap → the documented all-zero
///     branch (§6.1: all fractions 0, BestNu 0, Correlation 0, Rmse 0, OverlappingGenes 0).
///   • "all-zero mixture" (BE) — every overlapping expression value is 0 ⇒ the mixture has zero
///     variance, so <c>Standardize</c> hits its <c>sd &lt; 1e-15</c> guard and returns all-zeros (no
///     DivideByZero); the SMO target is all-zero ⇒ weights 0 ⇒ defined result, no NaN.
///   • "empty matrix" / empty profile (BE) — an empty signature matrix (0 cell types) or empty profile
///     (0 overlap) → the documented all-zero / no-overlap branch, never a crash.
///   • "negative expression" (MC) — log2-transformed expression legitimately goes negative; standardise
///     + regress is defined; fractions stay in [0,1] and sum to 1 (or all-zero).
///   • "NaN" / Infinity (MC) — a NaN / ±Inf expression value must NOT silently propagate into the
///     fractions: the contract fields (CellFractions, BestNu, Correlation, Rmse) must stay finite, or
///     the call raises a documented exception. A NaN fraction leaking out is a REAL bug.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (independently re-derived from §2.C / §4.C, NOT read off code)
/// ───────────────────────────────────────────────────────────────────────────
/// Immune_Infiltration_Estimation.md §2.C / §4.C / §6.1:
///   • The mixture vector m and each column of the signature matrix B are z-score STANDARDISED, then a
///     linear-kernel ν-SVR is solved for each ν ∈ {0.25, 0.5, 0.75}; the ν minimising the reconstruction
///     RMSE is selected; negative weights are zero-clipped and the remainder normalised to sum 1.
///   • INV-NUSVR-01: reported cell fractions are non-negative (negative weights clipped to 0).
///   • INV-NUSVR-02: fractions sum to 1 when total post-clip mass is positive (else all-zero).
///   • INV-NUSVR-03: BestNu ∈ {0.25, 0.5, 0.75} when a fit was performed (0 in the no-overlap branch).
///   • §6.1 no-overlap branch: all fractions 0, BestNu 0, Correlation 0, Rmse 0, OverlappingGenes 0.
///   • Correlation is a Pearson correlation ⇒ ∈ [-1, 1]; Rmse ≥ 0; determinism (same input → identical
///     result — the SMO solver is deterministic, no shared RNG).
///
/// ───────────────────────────────────────────────────────────────────────────
/// CRITICAL standardisation-conditioning note (do NOT over-claim recovery)
/// ───────────────────────────────────────────────────────────────────────────
/// ν-SVR z-score-standardises EACH signature column before regression. That is only well-posed with an
/// adequate marker-gene count (CIBERSORT/LM22 = 547 genes). On a degenerate 2–4-gene toy mixture,
/// standardisation fills disjoint-support zeros with large negatives and DESTROYS the columns'
/// orthogonality, so dominant-component recovery is NOT guaranteed and can fully invert — that is
/// CORRECT behaviour for this conditioning, not a defect. These fuzz tests therefore assert ONLY robust
/// invariants that hold regardless of conditioning (fractions ∈ [0,1], sum to 1 or all-zero, finite
/// metrics, termination, determinism) and DELIBERATELY do not assert that ν-SVR recovers the planted /
/// dominant cell type on a tiny mixture.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng. The assembly bootstrap
/// (LimitationPolicyTestBootstrap, ModuleInitializer) sets DefaultMode = Permissive, so the
/// ONCO-IMMUNE-001 LimitationPolicy guard in DeconvoluteImmuneCellsNuSvr is a no-op here (it only throws
/// under Strict).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyNuSvrFuzzTests
{
    private static readonly IReadOnlyList<double> DefaultNu = new[] { 0.25, 0.5, 0.75 };

    // Genes present in the DefaultSignatureMatrix's representative markers (used to build profiles that
    // genuinely overlap the built-in matrix so the regression path — not just the no-overlap branch — runs).
    private static readonly string[] SignatureGenes =
        { "CD8A", "CD8B", "CD3D", "CD3E", "GZMB", "PRF1", "NKG7", "MS4A1", "CD19", "CD14", "FCGR3A" };

    // ── Robust-invariant assertion helper ────────────────────────────────────
    // Pins the documented ν-SVR contract on EVERY result regardless of conditioning: each fraction is
    // finite and in [0,1]; the fractions sum to 1 (positive post-clip mass) or are all-zero (no
    // overlap / no surviving mass); BestNu is one of the swept ν (or 0 in the no-overlap branch);
    // Correlation ∈ [-1,1]; Rmse ≥ 0 and finite; OverlappingGenes ≥ 0. This is what stops a fuzz test
    // from rubber-stamping a NaN-corrupted or out-of-[0,1] result green.
    private static void AssertWellFormedNuSvr(NuSvrDeconvolutionResult r, IReadOnlyList<double> sweptNu)
    {
        const double Eps = 1e-9;

        r.CellFractions.Should().NotBeNull();
        foreach (var f in r.CellFractions.Values)
        {
            double.IsNaN(f).Should().BeFalse("a ν-SVR cell fraction must never be NaN");
            double.IsInfinity(f).Should().BeFalse("a ν-SVR cell fraction must be finite");
            f.Should().BeGreaterThanOrEqualTo(-Eps, "INV-NUSVR-01: fractions are non-negative (zero-clipped)");
            f.Should().BeLessThanOrEqualTo(1.0 + Eps, "a normalised fraction cannot exceed 1");
        }

        double sum = r.CellFractions.Values.Sum();
        // INV-NUSVR-02: sum to 1 when there is positive mass, else all-zero.
        (Math.Abs(sum - 1.0) < 1e-6 || Math.Abs(sum) < Eps).Should().BeTrue(
            $"fraction sum must be ~1 (positive mass) or ~0 (no-overlap / no surviving mass); was {sum}");

        // INV-NUSVR-03: BestNu ∈ swept set, or 0 in the no-overlap branch.
        (Math.Abs(r.BestNu) < Eps || sweptNu.Any(nu => Math.Abs(nu - r.BestNu) < Eps)).Should().BeTrue(
            $"BestNu must be one of the swept ν values or 0 (no-overlap); was {r.BestNu}");

        double.IsNaN(r.Correlation).Should().BeFalse("Correlation must never be NaN");
        double.IsInfinity(r.Correlation).Should().BeFalse("Correlation must be finite");
        r.Correlation.Should().BeGreaterThanOrEqualTo(-1.0 - Eps, "Pearson correlation ∈ [-1, 1]");
        r.Correlation.Should().BeLessThanOrEqualTo(1.0 + Eps, "Pearson correlation ∈ [-1, 1]");

        double.IsNaN(r.Rmse).Should().BeFalse("Rmse must never be NaN");
        double.IsInfinity(r.Rmse).Should().BeFalse("Rmse must be finite");
        r.Rmse.Should().BeGreaterThanOrEqualTo(0.0, "RMSE is a root-mean-square, never negative");

        r.OverlappingGenes.Should().BeGreaterThanOrEqualTo(0);
    }

    // Builds a profile keyed by the supplied genes with values drawn from a generator.
    private static Dictionary<string, double> Profile(IEnumerable<string> genes, Func<string, double> value)
    {
        var d = new Dictionary<string, double>();
        foreach (var g in genes)
            d[g] = value(g);
        return d;
    }

    #region IMMUNE-NUSVR-001 — null / documented-precondition rejection

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_NullProfile_ThrowsArgumentNullException()
    {
        // Docs §3.1 / §3.3: a null expression profile is rejected before any work.
        Action act = () => DeconvoluteImmuneCellsNuSvr(null!);
        act.Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region IMMUNE-NUSVR-001 — BE: gene-count mismatch (partial / zero overlap)

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_NoOverlappingGenes_AllZeroNoOverlapBranch()
    {
        // §6.1: when no profile gene is in the signature matrix, the method returns the documented
        // no-overlap branch — all fractions 0, BestNu 0, Correlation 0, Rmse 0, OverlappingGenes 0 —
        // and must NOT throw (no KeyNotFound, no DivideByZero on an empty regression).
        var profile = new Dictionary<string, double>
        {
            ["GHOST_GENE_1"] = 9.0, ["GHOST_GENE_2"] = 8.0, ["GHOST_GENE_3"] = 7.0
        };

        Action act = () => DeconvoluteImmuneCellsNuSvr(profile);
        act.Should().NotThrow();

        var result = DeconvoluteImmuneCellsNuSvr(profile);
        AssertWellFormedNuSvr(result, DefaultNu);
        result.OverlappingGenes.Should().Be(0);
        result.BestNu.Should().Be(0.0);
        result.Correlation.Should().Be(0.0);
        result.Rmse.Should().Be(0.0);
        result.CellFractions.Values.Should().OnlyContain(f => f == 0.0,
            "the no-overlap branch returns all-zero fractions");
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_EmptyProfile_AllZeroNoOverlapBranch()
    {
        // An empty profile shares no genes with any signature → no-overlap branch.
        var result = DeconvoluteImmuneCellsNuSvr(new Dictionary<string, double>());
        AssertWellFormedNuSvr(result, DefaultNu);
        result.OverlappingGenes.Should().Be(0);
        result.BestNu.Should().Be(0.0);
        result.CellFractions.Values.Should().OnlyContain(f => f == 0.0);
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_PartialOverlap_StaysInContract()
    {
        // Profile shares only a SUBSET of the signature genes (gene-count mismatch). The regression
        // runs on the intersection only; the result must still satisfy every robust invariant — but we
        // deliberately do NOT assert any particular cell type is recovered (degenerate conditioning).
        var partial = SignatureGenes.Take(3).ToArray();
        var profile = Profile(partial, g => 10.0 + partial.ToList().IndexOf(g));
        // Add some genes the signature has never heard of, to mix the overlap.
        profile["MYSTERY_A"] = 3.0;
        profile["MYSTERY_B"] = 2.0;

        var result = DeconvoluteImmuneCellsNuSvr(profile);
        AssertWellFormedNuSvr(result, DefaultNu);
        result.OverlappingGenes.Should().BeGreaterThan(0, "the profile overlaps a subset of signature genes");
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_EmptySignatureMatrix_NoCellTypesBranch()
    {
        // BE: an explicitly EMPTY signature matrix (0 cell types) → the cellTypes.Count == 0 guard fires,
        // returning the all-zero no-overlap branch with an empty fraction dictionary. Must not crash.
        var emptyMatrix = new Dictionary<string, IReadOnlyDictionary<string, double>>();
        var profile = Profile(SignatureGenes, _ => 5.0);

        Action act = () => DeconvoluteImmuneCellsNuSvr(profile, emptyMatrix);
        act.Should().NotThrow();

        var result = DeconvoluteImmuneCellsNuSvr(profile, emptyMatrix);
        AssertWellFormedNuSvr(result, DefaultNu);
        result.OverlappingGenes.Should().Be(0);
        result.CellFractions.Should().BeEmpty("no cell types ⇒ no fractions");
    }

    #endregion

    #region IMMUNE-NUSVR-001 — BE: all-zero mixture (std=0 guard, no DivideByZero)

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_AllZeroMixture_DefinedNoDivideByZero()
    {
        // BE: every overlapping expression value is 0 ⇒ zero-variance mixture. Standardize must hit its
        // sd < 1e-15 guard and return all-zeros (no DivideByZero); the SMO target is all-zero so the
        // solver finds zero weights. Result must be well-formed and finite throughout.
        var profile = Profile(SignatureGenes, _ => 0.0);

        Action act = () => DeconvoluteImmuneCellsNuSvr(profile);
        act.Should().NotThrow();

        var result = DeconvoluteImmuneCellsNuSvr(profile);
        AssertWellFormedNuSvr(result, DefaultNu);
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_ConstantNonZeroMixture_DefinedNoDivideByZero()
    {
        // A constant (but non-zero) mixture also has zero variance → same sd=0 guard. Must be defined.
        var profile = Profile(SignatureGenes, _ => 42.0);

        Action act = () => DeconvoluteImmuneCellsNuSvr(profile);
        act.Should().NotThrow();
        AssertWellFormedNuSvr(DeconvoluteImmuneCellsNuSvr(profile), DefaultNu);
    }

    #endregion

    #region IMMUNE-NUSVR-001 — MC: negative expression

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_NegativeExpression_StaysInContract()
    {
        // MC: log2-transformed expression legitimately goes negative. Standardise + regress is defined;
        // fractions ∈ [0,1] summing to 1 (or all-zero), metrics finite.
        var profile = Profile(SignatureGenes, g => -(g.Length) - 0.5);

        Action act = () => DeconvoluteImmuneCellsNuSvr(profile);
        act.Should().NotThrow();
        AssertWellFormedNuSvr(DeconvoluteImmuneCellsNuSvr(profile), DefaultNu);
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_RandomMixedSignExpression_NeverCorruptsResult()
    {
        // Fuzz: random profiles mixing large positive, negative, and zero expression across many
        // LOCALLY-seeded runs. Every result must satisfy every robust invariant.
        for (int seed = 0; seed < 80; seed++)
        {
            var rng = new Random(seed);
            var profile = Profile(SignatureGenes, _ => (rng.NextDouble() - 0.5) * 2_000.0);
            int filler = rng.Next(0, 6);
            for (int i = 0; i < filler; i++)
                profile["extra" + i] = (rng.NextDouble() - 0.5) * 2_000.0;

            var result = DeconvoluteImmuneCellsNuSvr(profile);
            AssertWellFormedNuSvr(result, DefaultNu);
        }
    }

    #endregion

    #region IMMUNE-NUSVR-001 — MC: NaN / Infinity expression (no silent NaN propagation)

    // NOTE on the NaN/Infinity contract (source fix made under this row, see class header / git):
    // A NaN or ±Infinity expression value would propagate through z-score standardisation, the SMO
    // solver, and the fit diagnostics, leaking a NaN/Infinity into the contracted-finite output fields
    // (Correlation, Rmse — verified pre-fix). Per §2.C the linear-mixture model is defined only for
    // finite expression, so DeconvoluteImmuneCellsNuSvr now rejects non-finite input up front with a
    // documented ArgumentException rather than emitting a corrupted result. These tests pin that
    // documented validation branch.

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_AllNaNExpression_ThrowsDocumentedArgumentException()
    {
        // MC headline hazard: a NaN expression value must NOT silently propagate into the contract.
        // It is rejected with a documented ArgumentException (never a NaN-corrupted result).
        var profile = Profile(SignatureGenes, _ => double.NaN);

        Action act = () => DeconvoluteImmuneCellsNuSvr(profile);
        act.Should().Throw<ArgumentException>()
            .Which.ParamName.Should().Be("expressionProfile");
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_InfinityExpression_ThrowsDocumentedArgumentException()
    {
        var profile = new Dictionary<string, double>();
        bool flip = false;
        foreach (var g in SignatureGenes)
        {
            profile[g] = flip ? double.PositiveInfinity : double.NegativeInfinity;
            flip = !flip;
        }

        Action act = () => DeconvoluteImmuneCellsNuSvr(profile);
        act.Should().Throw<ArgumentException>()
            .Which.ParamName.Should().Be("expressionProfile");
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_SingleNonFiniteAmongFinite_StillRejected()
    {
        // Even ONE non-finite value anywhere in the profile triggers the documented rejection — the
        // guard scans every value, not just the first.
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(seed);
            var profile = Profile(SignatureGenes, _ => rng.NextDouble() * 100.0);
            // Poison exactly one gene with a randomly chosen non-finite value.
            string poisoned = SignatureGenes[rng.Next(SignatureGenes.Length)];
            profile[poisoned] = rng.Next(0, 3) switch
            {
                0 => double.NaN,
                1 => double.PositiveInfinity,
                _ => double.NegativeInfinity
            };

            Action act = () => DeconvoluteImmuneCellsNuSvr(profile);
            act.Should().Throw<ArgumentException>()
                .Which.ParamName.Should().Be("expressionProfile");
        }
    }

    #endregion

    #region IMMUNE-NUSVR-001 — solver termination (SMO must not hang)

    [Test]
    [CancelAfter(30_000)]
    public void DeconvoluteImmuneCellsNuSvr_DenseManyGeneMixture_SolverTerminates()
    {
        // The SMO solver has a maxIterations = 200·ℓ cap; on a dense, larger, ill-conditioned mixture it
        // must still TERMINATE (not hang) and return an in-contract result. [CancelAfter] fails the test
        // if the solver were to spin past the cap.
        var rng = new Random(1234);
        var genes = SignatureGenes.Concat(Enumerable.Range(0, 40).Select(i => "extra" + i)).ToArray();
        var profile = Profile(genes, _ => rng.NextDouble() * 1_000.0);

        var result = DeconvoluteImmuneCellsNuSvr(profile);
        AssertWellFormedNuSvr(result, DefaultNu);
    }

    [Test]
    [CancelAfter(30_000)]
    public void DeconvoluteImmuneCellsNuSvr_PathologicalNearDuplicateGenes_SolverTerminates()
    {
        // Near-identical rows make the kernel near-singular (tiny curvature) — a classic non-termination
        // trap for SMO. The curvature floor (Eps) + maxIterations cap must keep it terminating.
        var profile = Profile(SignatureGenes, _ => 100.0 + 1e-12);
        profile["CD8A"] = 100.0;  // a hair of variance so std-guard does not trivially zero it

        var result = DeconvoluteImmuneCellsNuSvr(profile);
        AssertWellFormedNuSvr(result, DefaultNu);
    }

    #endregion

    #region IMMUNE-NUSVR-001 — determinism + custom ν sweep boundaries

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_SameInput_IsDeterministic()
    {
        // The SMO solver is deterministic (no shared RNG); identical input → identical result.
        var rng = new Random(7);
        var profile = Profile(SignatureGenes, _ => rng.NextDouble() * 500.0);

        var a = DeconvoluteImmuneCellsNuSvr(profile);
        var b = DeconvoluteImmuneCellsNuSvr(profile);

        a.BestNu.Should().Be(b.BestNu);
        a.OverlappingGenes.Should().Be(b.OverlappingGenes);
        a.Correlation.Should().Be(b.Correlation);
        a.Rmse.Should().Be(b.Rmse);
        a.CellFractions.Keys.Should().BeEquivalentTo(b.CellFractions.Keys);
        foreach (var key in a.CellFractions.Keys)
            a.CellFractions[key].Should().Be(b.CellFractions[key]);
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_CustomSingleNu_BestNuIsThatNu()
    {
        // BE: a single-element ν sweep. When a fit is performed, BestNu must be exactly that ν
        // (INV-NUSVR-03 — the sweep only considers the supplied values).
        var rng = new Random(99);
        var profile = Profile(SignatureGenes, _ => rng.NextDouble() * 300.0 + 1.0);
        var customNu = new[] { 0.4 };

        var result = DeconvoluteImmuneCellsNuSvr(profile, signatureMatrix: null, nuValues: customNu);
        AssertWellFormedNuSvr(result, customNu);
        if (result.OverlappingGenes > 0)
            result.BestNu.Should().Be(0.4, "a single-ν sweep selects that ν when a fit runs");
    }

    [Test]
    public void DeconvoluteImmuneCellsNuSvr_BoundaryNuValues_StayInContract()
    {
        // BE: ν at the edges of (0,1]. The doc constrains ν ∈ (0,1]; sweep near both ends and assert the
        // result stays in contract (no hang, finite, in-[0,1] fractions).
        var rng = new Random(2024);
        var profile = Profile(SignatureGenes, _ => rng.NextDouble() * 200.0 + 5.0);
        var edgeNu = new[] { 0.01, 1.0 };

        var result = DeconvoluteImmuneCellsNuSvr(profile, signatureMatrix: null, nuValues: edgeNu);
        AssertWellFormedNuSvr(result, edgeNu);
    }

    #endregion
}
