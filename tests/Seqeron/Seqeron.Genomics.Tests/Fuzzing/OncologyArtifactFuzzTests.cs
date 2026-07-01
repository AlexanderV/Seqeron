using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology sequencing-artifact-filtering area — ONCO-ARTIFACT-001.
/// The unit under test is the deterministic, rule-based artifact classifier
/// <see cref="OncologyAnalyzer.ClassifyArtifact"/> and its building blocks /
/// collection wrappers
/// <see cref="OncologyAnalyzer.CalculateGivScore"/> (GIV read-orientation imbalance),
/// <see cref="OncologyAnalyzer.CalculateStrandBias"/> (Phred-scaled GATK FisherStrand FS),
/// <see cref="OncologyAnalyzer.DetectOxoGArtifacts"/> and
/// <see cref="OncologyAnalyzer.FilterArtifacts"/>, implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime fault
/// (DivideByZero / Overflow / NaN / unguarded Infinity). Every input must resolve
/// to EITHER a well-defined, theory-correct value OR a *documented, intentional*
/// outcome (here, an <see cref="ArgumentOutOfRangeException"/> for a negative read
/// count or an <see cref="ArgumentNullException"/> for a null variant enumerable).
/// For artifact filtering the headline hazards — exactly the BE targets of row 90
/// ("zero depth, extreme strand bias, all-pass, all-fail") — are:
///   • ZERO DEPTH: a strand table with no reads (or no reads on a strand) and a
///     GIV with a zero read-mate count must NOT throw DivideByZero and must NOT
///     leak a NaN. Documented guards: an all-zero strand table ⇒ Fisher p = 1 ⇒
///     FS = 0 (§3.3, §6.1); GIV(0, 0) ⇒ 1.0 (no imbalance evidence, §6.1).
///   • EXTREME STRAND BIAS: an allele fully segregated to one strand drives the
///     Fisher p toward 0, so FS = −10·log10(max(p, MIN_PVALUE)) SATURATES to a
///     large but FINITE value — MIN_PVALUE = 1E-320 caps it; FS must never be +∞
///     or NaN (§2.2, §3.3). Likewise GIV with a zero read-2 count and a positive
///     read-1 count is the documented maximal one-sided imbalance ⇒ +∞, and that
///     +∞ correctly FLAGS an OxoG substitution as an artifact (GIV > 1.5).
///   • ALL-PASS: a balanced, well-supported NON-artifact substitution set (e.g.
///     A>G transitions, ~50/50 strand split, balanced read mates) must ALL pass
///     the filter — none flagged, FilterArtifacts is the identity on them.
///   • ALL-FAIL: a set of clear artifacts (FFPE C>T/G>A always; OxoG G>T/C>A with
///     GIV > 1.5) must ALL be flagged — FilterArtifacts returns empty.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-ARTIFACT-001 — Sequencing artifact detection / filtering (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 90.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 90): "zero depth, extreme strand bias, all-pass,
///     all-fail".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Sequencing_Artifact_Detection.md (docs/algorithms/Oncology/Sequencing_Artifact_Detection.md):
///   • Substitution classes (disjoint): FFPE deamination {C>T, G>A};
///       OxoG oxidation {G>T, C>A}; everything else None.                  (§2.2, §4.2, INV-04)
///   • GIV = read1Alt / read2Alt; GIV = 1 balanced; GIV > 1.5 = damaged.   (§2.2, §4.2)
///       GIV(0, 0) = 1.0 (no imbalance); GIV(>0, 0) = +∞ (maximal).        (§3.3, §6.1, INV-02)
///   • FS = −10·log10(max(p, MIN_PVALUE)) on the [refFwd, refRev, altFwd,
///       altRev] 2×2 table; FS = 0 for a balanced table (p = 1); FS ≥ 0;
///       MIN_PVALUE = 1E-320 caps FS finite.                               (§2.2, §3.3, INV-03)
///   • Flag rule: FFPE ⇒ always artifact; OxoG ⇒ artifact iff GIV > 1.5;
///       None ⇒ never.                                                     (§4.1, §4.2)
///   • FilterArtifacts output ⊆ input, input order (non-flagged subset).   (§3.2, INV-01)
///   • DetectOxoGArtifacts returns the flagged OxoG calls (GIV > 1.5).     (§3.2)
///   • Bases upper-cased before classification (case-insensitive).         (§3.3)
///   • Negative count ⇒ ArgumentOutOfRangeException; null variants ⇒
///       ArgumentNullException; empty variant list ⇒ empty result.         (§3.3, §6.1)
///   • Worked numbers (independently derived, §7.1 / sibling unit tests):
///       table [20,0,0,20]  ⇒ FS = 108.38365838736458;
///       table [15,5,5,15]  ⇒ FS = 24.148182890180962;
///       table [10,10,10,10]⇒ FS = 0; table [0,0,0,0] ⇒ FS = 0.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyArtifactFuzzTests
{
    private const double GivThreshold = DamagedGivThreshold; // 1.5

    // The four artifact substitutions and a representative non-artifact set.
    private static readonly (char Ref, char Alt)[] FfpeSubs = { ('C', 'T'), ('G', 'A') };
    private static readonly (char Ref, char Alt)[] OxoGSubs = { ('G', 'T'), ('C', 'A') };
    private static readonly (char Ref, char Alt)[] NonArtifactSubs =
    {
        ('A', 'G'), ('T', 'C'), ('A', 'T'), ('T', 'A'), ('A', 'C'), ('T', 'G'), ('G', 'C'), ('C', 'G'),
    };

    private static ArtifactObservation Obs(
        char r, char a,
        int rf = 10, int rr = 10, int af = 10, int ar = 10,
        int r1 = 10, int r2 = 10) =>
        new(r, a, rf, rr, af, ar, r1, r2);

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented numeric contract on EVERY call: FS must be FINITE and
    // ≥ 0 (no unguarded +∞ from a p-value of 0, no NaN from a zero-depth table);
    // GIV must be ≥ 0 and either finite or the documented +∞ sentinel (never NaN);
    // and the IsArtifact flag must EXACTLY follow the substitution-class + GIV rule.
    // This is what stops a fuzz test from rubber-stamping a NaN / unguarded-Inf /
    // mis-flagged result green.
    private static void AssertWellFormedCall(ArtifactCall call)
    {
        // FS: finite, non-negative (§2.2, INV-03; MIN_PVALUE keeps it finite).
        double.IsNaN(call.StrandBiasPhred).Should().BeFalse("FS must never be NaN");
        double.IsInfinity(call.StrandBiasPhred).Should().BeFalse(
            "FS is capped finite by MIN_PVALUE = 1E-320 (§3.3) — never an unguarded +∞");
        call.StrandBiasPhred.Should().BeGreaterThanOrEqualTo(0.0, "FS = −10·log10(p≤1) ≥ 0 (INV-03)");

        // GIV: non-negative, never NaN; +∞ is the documented maximal-imbalance sentinel only.
        double.IsNaN(call.GivScore).Should().BeFalse("GIV must never be NaN (§3.3)");
        call.GivScore.Should().BeGreaterThanOrEqualTo(0.0, "GIV is a ratio of non-negative counts (INV-02)");

        // The flag must be the exact deterministic rule (§4.2).
        bool expected = call.Type switch
        {
            ArtifactType.FfpeDeamination => true,
            ArtifactType.OxoG => call.GivScore > GivThreshold,
            _ => false,
        };
        call.IsArtifact.Should().Be(expected,
            "IsArtifact must follow: FFPE⇒always, OxoG⇒GIV>1.5, None⇒never (§4.2)");
    }

    // Confirms FilterArtifacts output is an order-preserving subsequence of the input (INV-01).
    private static void AssertOrderedSubsetOfInput(
        IReadOnlyList<ArtifactObservation> kept, IReadOnlyList<ArtifactObservation> input)
    {
        kept.Count.Should().BeLessThanOrEqualTo(input.Count, "kept ⊆ input (INV-01)");

        int j = 0;
        foreach (var k in kept)
        {
            while (j < input.Count && !input[j].Equals(k))
            {
                j++;
            }

            j.Should().BeLessThan(input.Count, "each kept variant must appear in input order (INV-01)");
            j++;
        }
    }

    #region ONCO-ARTIFACT-001 — positive sanity (balanced passes; biased/damaged is flagged; FS hand-computed)

    // A balanced A>G transition (non-artifact substitution, 50/50 strand split, balanced read mates)
    // must NOT be flagged: Type = None, FS = 0, GIV = 1, and FilterArtifacts keeps it.
    [Test]
    public void ClassifyArtifact_BalancedNonArtifact_NotFlagged_HandComputed()
    {
        var v = Obs('A', 'G', rf: 20, rr: 20, af: 20, ar: 20, r1: 50, r2: 50);

        var call = ClassifyArtifact(v);

        call.Type.Should().Be(ArtifactType.None, "A>G is not an FFPE/OxoG substitution class (§2.2)");
        call.IsArtifact.Should().BeFalse("a balanced non-artifact substitution is a candidate true variant");
        call.StrandBiasPhred.Should().BeApproximately(0.0, 1e-9,
            "balanced table [20,20,20,20] ⇒ p ≈ 1 ⇒ FS ≈ 0 (INV-03)");
        call.GivScore.Should().Be(1.0, "equal read1/read2 alt counts ⇒ GIV = 1 (INV-02)");
        FilterArtifacts(new[] { v }).Should().ContainSingle().Which.Should().Be(v);
    }

    // A canonical OxoG G>T with a strand-biased table and a damaging GIV (200/100 = 2.0 > 1.5) must be
    // flagged AND the FS must equal the independently-derived documented Phred score for [20,0,0,20].
    [Test]
    public void ClassifyArtifact_StrandBiasedDamagedOxoG_IsFlagged_FsHandComputed()
    {
        // Worked example (§7.1): strand table [20,0,0,20] ⇒ p = 1.4508889e-11 ⇒ FS = 108.38365838736458.
        var v = Obs('G', 'T', rf: 20, rr: 0, af: 0, ar: 20, r1: 200, r2: 100);

        var call = ClassifyArtifact(v);

        call.Type.Should().Be(ArtifactType.OxoG, "G>T is the OxoG substitution class (§2.2)");
        call.GivScore.Should().Be(2.0, "GIV = 200/100 = 2.0 > 1.5 ⇒ damaged (§4.2)");
        call.IsArtifact.Should().BeTrue("OxoG with GIV > 1.5 is flagged (§4.2)");
        call.StrandBiasPhred.Should().BeApproximately(108.38365838736458, 1e-6,
            "table [20,0,0,20] ⇒ two-sided Fisher p = 1.4508889e-11 ⇒ FS = −10·log10(p) (§7.1)");
    }

    // FFPE C>T is flagged by substitution class ALONE, regardless of GIV / strand balance (§4.2).
    [Test]
    public void ClassifyArtifact_FfpeDeamination_AlwaysFlagged_EvenWhenBalanced()
    {
        var v = Obs('C', 'T', rf: 30, rr: 30, af: 25, ar: 25, r1: 40, r2: 40);

        var call = ClassifyArtifact(v);

        call.Type.Should().Be(ArtifactType.FfpeDeamination);
        call.IsArtifact.Should().BeTrue("FFPE deamination is flagged by class alone (§4.2)");
        call.GivScore.Should().Be(1.0, "balanced read mates ⇒ GIV = 1 (does not gate FFPE)");
        call.StrandBiasPhred.Should().BeApproximately(0.0, 1e-9,
            "near-balanced strand table ⇒ p ≈ 1 ⇒ FS ≈ 0 (FFPE flag is class-only)");
    }

    // OxoG substitution but a BALANCED GIV (1.0 ≤ 1.5) ⇒ NOT flagged: GIV gates OxoG (§4.2).
    [Test]
    public void ClassifyArtifact_OxoGSubstitution_BalancedGiv_NotFlagged()
    {
        var v = Obs('C', 'A', rf: 15, rr: 15, af: 12, ar: 13, r1: 30, r2: 30);

        var call = ClassifyArtifact(v);

        call.Type.Should().Be(ArtifactType.OxoG, "C>A is the OxoG class");
        call.GivScore.Should().Be(1.0);
        call.IsArtifact.Should().BeFalse("OxoG with GIV ≤ 1.5 is NOT confirmed as damage (§4.2)");
    }

    // The partially-biased table FS is also pinned to an independently-derived documented value (INV-05).
    [Test]
    public void CalculateStrandBias_PartialBias_FsHandComputed_AndMonotone()
    {
        double balanced = CalculateStrandBias(10, 10, 10, 10);
        double partial = CalculateStrandBias(15, 5, 5, 15);
        double full = CalculateStrandBias(20, 0, 0, 20);

        balanced.Should().BeApproximately(0.0, 1e-9, "p ≈ 1 ⇒ FS ≈ 0 (documented numeric tolerance)");
        partial.Should().BeApproximately(24.148182890180962, 1e-6,
            "table [15,5,5,15] ⇒ p = 0.0038475273 ⇒ FS = 24.1481829");
        full.Should().BeApproximately(108.38365838736458, 1e-6);

        // INV-05: greater strand segregation ⇒ FS non-decreasing.
        partial.Should().BeGreaterThanOrEqualTo(balanced);
        full.Should().BeGreaterThanOrEqualTo(partial);
    }

    // Case-insensitive classification: lowercase ref/alt must classify identically (§3.3).
    [Test]
    public void ClassifyArtifact_LowercaseBases_ClassifiedSameAsUppercase()
    {
        var upper = ClassifyArtifact(Obs('G', 'T', r1: 200, r2: 100));
        var lower = ClassifyArtifact(Obs('g', 't', r1: 200, r2: 100));

        lower.Type.Should().Be(upper.Type, "bases are upper-cased before classification (§3.3)");
        lower.IsArtifact.Should().Be(upper.IsArtifact);
        lower.GivScore.Should().Be(upper.GivScore);
    }

    #endregion

    #region ONCO-ARTIFACT-001 — BE: zero depth (no DivideByZero, no NaN; documented guarded result)

    // An all-zero strand table has no evidence of bias: documented p = 1 ⇒ FS = 0 (§3.3, §6.1).
    // Must NOT throw DivideByZero and must NOT produce NaN.
    [Test]
    public void CalculateStrandBias_AllZeroTable_FsZero_NoThrow_NoNaN()
    {
        double fs = 0;
        FluentActions.Invoking(() => fs = CalculateStrandBias(0, 0, 0, 0))
            .Should().NotThrow("an empty strand table is a documented guarded case (§6.1)");

        double.IsNaN(fs).Should().BeFalse("zero-depth must not leak a NaN FS");
        fs.Should().Be(0.0, "no reads ⇒ Fisher p = 1 ⇒ FS = 0 (§3.3)");
    }

    // Zero reads on one strand for both alleles (a zero strand margin) ⇒ still p = 1, FS = 0.
    [Test]
    public void CalculateStrandBias_ZeroStrandMargin_FsZero_NoNaN()
    {
        double fs = CalculateStrandBias(10, 0, 10, 0); // all reads forward; reverse strand empty

        double.IsNaN(fs).Should().BeFalse();
        fs.Should().BeApproximately(0.0, 1e-9, "a zero-count strand margin gives no strand bias ⇒ p ≈ 1 ⇒ FS ≈ 0");
    }

    // GIV with BOTH read-mate counts zero (zero depth) ⇒ documented 1.0, NOT a DivideByZero / NaN (§6.1).
    [Test]
    public void CalculateGivScore_BothZero_ReturnsOne_NoDivideByZero_NoNaN()
    {
        double giv = 0;
        FluentActions.Invoking(() => giv = CalculateGivScore(0, 0))
            .Should().NotThrow("GIV(0,0) is the documented no-imbalance case (§6.1)");

        double.IsNaN(giv).Should().BeFalse("0/0 must be guarded to 1.0, never NaN");
        giv.Should().Be(1.0, "no read-mate evidence ⇒ no imbalance ⇒ GIV = 1 (§6.1, INV-02)");
    }

    // A zero-depth OxoG observation (no reads at all) must classify without throwing and NOT be flagged:
    // GIV guards to 1.0 (≤ 1.5) and FS guards to 0.
    [Test]
    public void ClassifyArtifact_ZeroDepthOxoG_NotFlagged_NoThrow_WellFormed()
    {
        var v = Obs('G', 'T', rf: 0, rr: 0, af: 0, ar: 0, r1: 0, r2: 0);

        ArtifactCall call = default;
        FluentActions.Invoking(() => call = ClassifyArtifact(v))
            .Should().NotThrow("zero depth must not crash the classifier");

        call.Type.Should().Be(ArtifactType.OxoG);
        call.GivScore.Should().Be(1.0, "GIV(0,0) ⇒ 1.0");
        call.StrandBiasPhred.Should().Be(0.0, "empty table ⇒ FS = 0");
        call.IsArtifact.Should().BeFalse("OxoG with GIV = 1.0 ≤ 1.5 is not confirmed damage");
        AssertWellFormedCall(call);
    }

    // Fuzz: a strand table that is empty on at least one strand or entirely zero, across substitution
    // classes — never DivideByZero, never NaN/Inf FS, always a well-formed call.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifyArtifact_RandomZeroDepthTables_NeverThrow_WellFormed()
    {
        var subs = FfpeSubs.Concat(OxoGSubs).Concat(NonArtifactSubs).ToArray();

        for (int seed = 0; seed < 400; seed++)
        {
            var rng = new Random(seed);
            var (r, a) = subs[rng.Next(subs.Length)];

            // Force at least one zero strand and/or zero read-mate count to probe the divisions.
            int rf = rng.Next(0, 2) == 0 ? 0 : rng.Next(0, 30);
            int rr = rng.Next(0, 2) == 0 ? 0 : rng.Next(0, 30);
            int af = rng.Next(0, 2) == 0 ? 0 : rng.Next(0, 30);
            int ar = rng.Next(0, 2) == 0 ? 0 : rng.Next(0, 30);
            int r1 = rng.Next(0, 2) == 0 ? 0 : rng.Next(0, 100);
            int r2 = rng.Next(0, 2) == 0 ? 0 : rng.Next(0, 100);

            var v = Obs(r, a, rf, rr, af, ar, r1, r2);

            ArtifactCall call = default;
            FluentActions.Invoking(() => call = ClassifyArtifact(v))
                .Should().NotThrow($"seed {seed}: zero-depth table {(r, a, rf, rr, af, ar, r1, r2)} must not crash");

            AssertWellFormedCall(call);
        }
    }

    #endregion

    #region ONCO-ARTIFACT-001 — BE: extreme strand bias / extreme GIV (saturates finite; no Inf leak in FS)

    // Fully-segregated alleles drive Fisher p toward 0, but MIN_PVALUE caps FS finite (§3.3).
    // FS must be large, FINITE, and never +∞ / NaN.
    [Test]
    public void CalculateStrandBias_ExtremeSegregation_SaturatesFinite_NoInfinity()
    {
        // Large, fully strand-segregated table: ref all forward, alt all reverse.
        double fs = CalculateStrandBias(500, 0, 0, 500);

        double.IsInfinity(fs).Should().BeFalse("MIN_PVALUE = 1E-320 caps FS finite (§3.3) — never +∞");
        double.IsNaN(fs).Should().BeFalse();
        fs.Should().BeGreaterThan(108.0, "extreme segregation ⇒ very small p ⇒ very large FS");

        // The absolute ceiling: −10·log10(1E-320) ≈ 3200. FS must stay strictly under it.
        double ceiling = -10.0 * Math.Log10(1E-320);
        fs.Should().BeLessThanOrEqualTo(ceiling + 1e-6, "FS cannot exceed the MIN_PVALUE-floored ceiling");
    }

    // Extreme GIV: positive read-1 count, zero read-2 count ⇒ documented maximal imbalance = +∞ (§6.1),
    // and that +∞ correctly FLAGS an OxoG variant as an artifact (+∞ > 1.5).
    [Test]
    public void CalculateGivScore_ZeroR2_PositiveR1_ReturnsPositiveInfinity_FlagsOxoG()
    {
        double giv = CalculateGivScore(50, 0);

        double.IsPositiveInfinity(giv).Should().BeTrue("read2 = 0, read1 > 0 ⇒ maximal one-sided imbalance (§6.1)");

        // The +∞ GIV must drive the OxoG flag — it is a guarded sentinel, compared as +∞ > 1.5.
        var call = ClassifyArtifact(Obs('G', 'T', r1: 50, r2: 0));
        call.GivScore.Should().Be(double.PositiveInfinity);
        call.IsArtifact.Should().BeTrue("OxoG with GIV = +∞ > 1.5 is flagged (the sentinel must compare, not crash)");
    }

    // Extreme strand bias on a NON-artifact substitution: FS saturates, but flagging is class-driven,
    // so the variant is still NOT flagged (FS is a reported metric, not the flag rule, §4.2).
    [Test]
    public void ClassifyArtifact_ExtremeStrandBias_NonArtifactSubstitution_StillNotFlagged()
    {
        var v = Obs('A', 'G', rf: 200, rr: 0, af: 0, ar: 200, r1: 100, r2: 100);

        var call = ClassifyArtifact(v);

        call.Type.Should().Be(ArtifactType.None);
        call.StrandBiasPhred.Should().BeGreaterThan(100.0, "extreme bias ⇒ large FS");
        double.IsInfinity(call.StrandBiasPhred).Should().BeFalse("FS stays finite even at extreme bias");
        call.IsArtifact.Should().BeFalse("None-class variants are never flagged; FS is reported, not decisive (§4.2)");
    }

    // Fuzz: extreme counts (incl. int.MaxValue strand cells and huge GIV ratios) — FS finite, GIV ≥ 0,
    // no overflow/Inf-leak, flag follows the rule. NB: total reads kept moderate to bound the Fisher sum;
    // int.MaxValue is probed on the GIV ratio path where it is a pure division.
    [Test]
    [CancelAfter(30_000)]
    public void ClassifyArtifact_ExtremeBiasAndGiv_NeverInfFs_WellFormed()
    {
        var subs = FfpeSubs.Concat(OxoGSubs).Concat(NonArtifactSubs).ToArray();

        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            var (r, a) = subs[rng.Next(subs.Length)];

            // One strand heavily loaded, the other near-empty ⇒ extreme bias (moderate magnitudes
            // so the Fisher feasible-cell sum stays fast; the saturation behaviour is identical).
            int heavy = rng.Next(50, 400);
            int light = rng.Next(0, 3);
            bool altReverseBias = rng.Next(0, 2) == 0;
            int rf = altReverseBias ? heavy : light;
            int rr = altReverseBias ? light : heavy;
            int af = altReverseBias ? light : heavy;
            int ar = altReverseBias ? heavy : light;

            // Extreme / boundary GIV: probe int.MaxValue and zero read-2.
            int r1 = rng.Next(0, 4) switch
            {
                0 => int.MaxValue,
                1 => rng.Next(1, 1000),
                2 => 0,
                _ => rng.Next(1, 1_000_000),
            };
            int r2 = rng.Next(0, 4) switch
            {
                0 => 0,            // ⇒ +∞ when r1 > 0
                1 => 1,
                2 => rng.Next(1, 10),
                _ => rng.Next(1, 1000),
            };

            var v = Obs(r, a, rf, rr, af, ar, r1, r2);

            ArtifactCall call = default;
            FluentActions.Invoking(() => call = ClassifyArtifact(v))
                .Should().NotThrow($"seed {seed}: extreme bias/GIV must not crash");

            // FS must remain FINITE (no leaked +∞ from a p-value of exactly 0).
            double.IsInfinity(call.StrandBiasPhred).Should().BeFalse(
                $"seed {seed}: FS must stay finite under extreme bias (MIN_PVALUE cap)");
            AssertWellFormedCall(call);
        }
    }

    #endregion

    #region ONCO-ARTIFACT-001 — BE: all-pass (every variant is a clean non-artifact ⇒ none flagged)

    // A set of balanced, well-supported NON-artifact substitutions ⇒ ALL pass: FilterArtifacts is identity,
    // DetectOxoGArtifacts is empty, no variant flagged.
    [Test]
    public void FilterArtifacts_AllNonArtifact_AllPass_FilterIsIdentity()
    {
        var input = NonArtifactSubs
            .Select(s => Obs(s.Ref, s.Alt, rf: 20, rr: 20, af: 18, ar: 17, r1: 40, r2: 38))
            .ToList();

        var kept = FilterArtifacts(input);

        kept.Should().Equal(input, "no variant is an artifact ⇒ the filter keeps every one, in order");
        DetectOxoGArtifacts(input).Should().BeEmpty("none is an OxoG-class artifact");
        input.Select(ClassifyArtifact).Should().OnlyContain(c => !c.IsArtifact);
    }

    // Even OxoG-class substitutions with BALANCED GIV (≤ 1.5) all pass — GIV gates OxoG flagging.
    [Test]
    public void FilterArtifacts_OxoGClassButBalancedGiv_AllPass()
    {
        var input = OxoGSubs
            .Select(s => Obs(s.Ref, s.Alt, rf: 25, rr: 24, af: 20, ar: 21, r1: 30, r2: 30))
            .ToList();

        FilterArtifacts(input).Should().Equal(input, "OxoG with GIV = 1.0 ≤ 1.5 is not confirmed damage (§4.2)");
        DetectOxoGArtifacts(input).Should().BeEmpty("none exceeds the GIV threshold");
    }

    // Fuzz all-pass: random balanced non-artifact variants — never flagged, filter is total identity.
    [Test]
    [CancelAfter(20_000)]
    public void FilterArtifacts_RandomNonArtifactSets_AlwaysAllPass()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(0, 25);
            var input = new List<ArtifactObservation>(n);
            for (int i = 0; i < n; i++)
            {
                var (r, a) = NonArtifactSubs[rng.Next(NonArtifactSubs.Length)];
                int depth = rng.Next(5, 50);
                // Roughly balanced strands + balanced read mates ⇒ no bias signal anyway, but the
                // None class means it can never be flagged regardless.
                input.Add(Obs(r, a, depth, depth + rng.Next(-2, 3), depth, depth + rng.Next(-2, 3),
                    depth, depth + rng.Next(-2, 3)));
            }

            var kept = FilterArtifacts(input);

            kept.Should().Equal(input, $"seed {seed}: non-artifact substitutions always pass");
            AssertOrderedSubsetOfInput(kept, input);
            input.Select(ClassifyArtifact).Should().OnlyContain(c => !c.IsArtifact,
                $"seed {seed}: no None-class variant may be flagged");
        }
    }

    #endregion

    #region ONCO-ARTIFACT-001 — BE: all-fail (every variant is a clear artifact ⇒ all flagged, filter empty)

    // A set of clear artifacts: FFPE C>T/G>A (always) + OxoG G>T/C>A with damaging GIV (2.0 > 1.5).
    // FilterArtifacts must return EMPTY; every variant flagged.
    [Test]
    public void FilterArtifacts_AllArtifacts_AllFail_FilterReturnsEmpty()
    {
        var input = new List<ArtifactObservation>();
        foreach (var (r, a) in FfpeSubs)
        {
            input.Add(Obs(r, a, r1: 40, r2: 40)); // FFPE flagged by class regardless of GIV
        }

        foreach (var (r, a) in OxoGSubs)
        {
            input.Add(Obs(r, a, r1: 200, r2: 100)); // GIV = 2.0 > 1.5 ⇒ damaged OxoG
        }

        FilterArtifacts(input).Should().BeEmpty("every variant is a flagged artifact ⇒ nothing survives the filter");
        input.Select(ClassifyArtifact).Should().OnlyContain(c => c.IsArtifact);

        // The two OxoG entries must surface from DetectOxoGArtifacts.
        DetectOxoGArtifacts(input).Should().HaveCount(2)
            .And.OnlyContain(c => c.Type == ArtifactType.OxoG && c.IsArtifact);
    }

    // OxoG with a zero read-2 count (GIV = +∞) is also a clear all-fail case.
    [Test]
    public void FilterArtifacts_OxoGInfiniteGiv_AllFail()
    {
        var input = OxoGSubs.Select(s => Obs(s.Ref, s.Alt, r1: 30, r2: 0)).ToList(); // GIV = +∞

        FilterArtifacts(input).Should().BeEmpty("GIV = +∞ > 1.5 ⇒ all flagged");
        DetectOxoGArtifacts(input).Should().HaveCount(input.Count);
    }

    // Fuzz all-fail: random mixes of FFPE (any GIV) and OxoG (forced GIV > 1.5) — always all flagged.
    [Test]
    [CancelAfter(20_000)]
    public void FilterArtifacts_RandomClearArtifactSets_AlwaysAllFail()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(1, 25);
            var input = new List<ArtifactObservation>(n);
            for (int i = 0; i < n; i++)
            {
                if (rng.Next(0, 2) == 0)
                {
                    // FFPE: flagged by class regardless of GIV / strand.
                    var (r, a) = FfpeSubs[rng.Next(FfpeSubs.Length)];
                    input.Add(Obs(r, a, rng.Next(0, 30), rng.Next(0, 30), rng.Next(0, 30), rng.Next(0, 30),
                        rng.Next(0, 100), rng.Next(0, 100)));
                }
                else
                {
                    // OxoG with a GIV strictly above 1.5: r1 = r2*2 + (1..) with r2 ≥ 1, or r2 = 0 (⇒ +∞).
                    var (r, a) = OxoGSubs[rng.Next(OxoGSubs.Length)];
                    int r2 = rng.Next(0, 5);
                    int r1 = r2 == 0 ? rng.Next(1, 100) : r2 * 2 + rng.Next(1, 50);
                    input.Add(Obs(r, a, rng.Next(0, 30), rng.Next(0, 30), rng.Next(0, 30), rng.Next(0, 30),
                        r1, r2));
                }
            }

            var kept = FilterArtifacts(input);

            kept.Should().BeEmpty($"seed {seed}: every variant is a clear artifact ⇒ filter empty");
            input.Select(ClassifyArtifact).Should().OnlyContain(c => c.IsArtifact,
                $"seed {seed}: all flagged");
        }
    }

    #endregion

    #region ONCO-ARTIFACT-001 — BE: empty / null / negative-count guards (documented intentional outcomes)

    [Test]
    public void FilterArtifacts_EmptyList_ReturnsEmpty_NoThrow()
    {
        FilterArtifacts(Array.Empty<ArtifactObservation>()).Should().BeEmpty();
    }

    [Test]
    public void DetectOxoGArtifacts_EmptyList_ReturnsEmpty_NoThrow()
    {
        DetectOxoGArtifacts(Array.Empty<ArtifactObservation>()).Should().BeEmpty();
    }

    [Test]
    public void NullVariants_ThrowArgumentNull()
    {
        FluentActions.Invoking(() => FilterArtifacts(null!))
            .Should().Throw<ArgumentNullException>();
        FluentActions.Invoking(() => DetectOxoGArtifacts(null!))
            .Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void NegativeCounts_ThrowArgumentOutOfRange()
    {
        FluentActions.Invoking(() => CalculateGivScore(-1, 10))
            .Should().Throw<ArgumentOutOfRangeException>();
        FluentActions.Invoking(() => CalculateGivScore(10, -1))
            .Should().Throw<ArgumentOutOfRangeException>();
        FluentActions.Invoking(() => CalculateStrandBias(-1, 0, 0, 0))
            .Should().Throw<ArgumentOutOfRangeException>();
        FluentActions.Invoking(() => CalculateStrandBias(0, -1, 0, 0))
            .Should().Throw<ArgumentOutOfRangeException>();
        FluentActions.Invoking(() => CalculateStrandBias(0, 0, -1, 0))
            .Should().Throw<ArgumentOutOfRangeException>();
        FluentActions.Invoking(() => CalculateStrandBias(0, 0, 0, -1))
            .Should().Throw<ArgumentOutOfRangeException>();
    }

    // Mixed fuzz: arbitrary blends of all classes + boundary depths — DetectOxoGArtifacts and
    // FilterArtifacts are complementary and consistent with per-variant ClassifyArtifact, never crash.
    [Test]
    [CancelAfter(30_000)]
    public void FilterAndDetect_RandomMixedSets_ConsistentWithClassify_NeverThrow()
    {
        var subs = FfpeSubs.Concat(OxoGSubs).Concat(NonArtifactSubs).ToArray();

        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(0, 30);
            var input = new List<ArtifactObservation>(n);
            for (int i = 0; i < n; i++)
            {
                var (r, a) = subs[rng.Next(subs.Length)];
                input.Add(Obs(r, a,
                    rng.Next(0, 40), rng.Next(0, 40), rng.Next(0, 40), rng.Next(0, 40),
                    rng.Next(0, 200), rng.Next(0, 200)));
            }

            IReadOnlyList<ArtifactObservation> kept = null!;
            IReadOnlyList<ArtifactCall> oxoG = null!;
            FluentActions.Invoking(() =>
            {
                kept = FilterArtifacts(input);
                oxoG = DetectOxoGArtifacts(input);
            }).Should().NotThrow($"seed {seed}: mixed set must not crash");

            AssertOrderedSubsetOfInput(kept, input);

            var calls = input.Select(ClassifyArtifact).ToList();
            foreach (var c in calls)
            {
                AssertWellFormedCall(c);
            }

            // kept = exactly the non-flagged variants (INV-01 + flag rule).
            kept.Count.Should().Be(calls.Count(c => !c.IsArtifact),
                $"seed {seed}: FilterArtifacts keeps exactly the non-flagged variants");

            // DetectOxoGArtifacts = exactly the flagged OxoG calls.
            oxoG.Count.Should().Be(calls.Count(c => c.Type == ArtifactType.OxoG && c.IsArtifact),
                $"seed {seed}: DetectOxoGArtifacts returns exactly the flagged OxoG calls");
            oxoG.Should().OnlyContain(c => c.Type == ArtifactType.OxoG && c.IsArtifact);
        }
    }

    #endregion
}
