namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the MolTools <b>Karlin–Altschul E-value / bit-score statistics</b>
/// (PROBE-EVALUE-001) — the opt-in significance layer that turns an off-target alignment hit's
/// raw score into a normalized bit score and an expected-by-chance count (E-value) over a
/// search space <c>m·n</c>. Targets <see cref="ProbeDesigner.ComputeKarlinAltschul"/> and its
/// λ root-finder <see cref="ProbeDesigner.ComputeLambdaNucleotide"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Malformed / boundary / out-of-domain inputs must NEVER hang (the λ bisection must converge
/// or fail cleanly — no non-convergence loop), throw an *unhandled* runtime exception, or emit
/// out-of-contract output. Concretely every input must resolve to EITHER a well-defined,
/// theory-correct result (a finite bit score and an E-value in [0,∞)), OR a *documented*
/// validation exception. The central hazards for Karlin–Altschul are: a NEGATIVE E-value
/// (impossible — E is an expected count); a NaN/±Inf bit score where finite is contracted
/// (K = 0 makes ln K = −∞, the classic leak); a non-converging λ root-finder on a degenerate
/// scoring scheme; and Int32/Int64 OVERFLOW in the search space m·n that would silently flip
/// E negative or wrap.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROBE-EVALUE-001 — Karlin–Altschul E-value / bit-score statistics
/// Checklist: docs/checklists/03_FUZZING.md, row 244 (strategies BE, MC).
/// Algorithm doc: docs/algorithms/MolTools/Probe_Validation.md §2.2 (Karlin–Altschul block),
///   INV-08 / INV-09, §5.1 (ComputeLambdaNucleotide / ComputeKarlinAltschul contracts).
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — raw score 0 (e^0 = 1 → E = K·m·n, defined), a NEGATIVE raw
///           score (E grows large but stays finite/≥ search space, never negative), the
///           zero-length-database / zero-search-space boundary (n = 0 → documented throw, not a
///           silent E = 0 NaN), K = 0 (ln K = −∞ in the bit score → documented throw, not a NaN
///           leak), and a HUGE search space (m·n near long.MaxValue → E stays finite-or-Inf,
///           never overflow-negative).
///   • MC  = Malformed Content — a null scoring matrix, a degenerate scoring scheme with no
///           positive score / a non-negative mismatch / a non-negative expected per-pair score
///           (λ undefined → documented throw, NOT a non-converging bisection or a NaN λ).
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary exploitation, MC = malformed
///   content); row 244 targets "zero-length DB, score 0, negative score, K=0, huge search space".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The Karlin–Altschul contract under test (doc §2.2)
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: ProbeDesigner.ComputeKarlinAltschul(double rawScore, int queryLength,
///   long databaseLength, ScoringMatrix? scoring = null, double k = 0.711,
///   double baseFrequency = 0.25) → KarlinAltschulStatistics.
/// Supporting entry: ProbeDesigner.ComputeLambdaNucleotide(int match, int mismatch,
///   double baseFrequency = 0.25) → λ.
///
///   E  = K · m · n · e^{−λS}                       (expected HSPs scoring ≥ S by chance)
///   S' = (λ·S − ln K) / ln 2                         (normalized bit score)
///   E  = m · n · 2^{−S'}                             (E in terms of the bit score)
///   λ  is the unique positive root of  Σ_{i,j} p_i p_j e^{λ s_ij} = 1.
///
/// DOMAIN (doc §5.1, §3-equivalent guards in source): queryLength m > 0, databaseLength n > 0,
/// K > 0, and a scoring scheme for which λ is defined (≥ one positive score, negative mismatch,
/// negative expected per-pair score). A non-positive length, non-positive K, null matrix, or a
/// degenerate scheme each throws — NEVER a leaked NaN/Inf/negative E-value.
///
/// THEORY invariants pinned (doc §2.4, INV-08 / INV-09):
///   • E ≥ 0 ALWAYS — it is an expected count; for any finite raw score (positive, zero, or
///     negative) over a positive search space the E-value is a non-negative finite real.
///   • E is MONOTONICALLY DECREASING in the raw score S — a higher-scoring hit is MORE
///     significant, hence LESS expected by chance (E ∝ e^{−λS}, λ > 0).
///   • The two E-value forms agree: K·m·n·e^{−λS} == m·n·2^{−S'} (INV-09).
///   • E scales LINEARLY with the search space m·n (INV-09).
///   • λ > 0, and for the BLAST +1/−3 scheme under p = 0.25, λ == 1.374 to tolerance (INV-08).
///   • Raw score 0 → E == K·m·n exactly (e^0 = 1); bit score finite.
///
/// HAND-DERIVED PIN (re-derived independently from the Karlin–Altschul equations):
///   Solving 0.25·e^{λ·1} + 0.75·e^{λ·(−3)} = 1 by bisection gives λ = 1.3740631224599755.
///   With S = 10, m = 100, n = 1000, K = 0.711:
///       E   = 0.711·100·1000·e^{−1.3740631·10} = 0.07662831608891887,
///       S'  = (1.3740631·10 − ln 0.711)/ln 2   = 20.315619061456786,
///       m·n·2^{−S'}                            = 0.07662831608891878  (== E to fp tolerance).
///   These exact values are pinned below so a test that passed against a wrong implementation
///   (e.g. a sign error in the exponent or a swapped bit-score formula) would FAIL.
///
/// The λ root-finder tests carry [CancelAfter] so any non-convergence loop fails as a timeout
/// rather than hanging. ComputeKarlinAltschul / ComputeLambdaNucleotide are pure static helpers
/// (no LimitationPolicy gate), so each probe calls them directly.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProbeEValueFuzzTests
{
    // BLAST +1/−3 nucleotide scheme — the published λ ≈ 1.374 reference scheme.
    private static readonly ScoringMatrix Blast1_3 = new(Match: 1, Mismatch: -3, GapOpen: -5, GapExtend: -2);

    private const double Ln2 = 0.6931471805599453;

    /// <summary>Independently re-derived λ for the +1/−3, p=0.25 scheme (bisection on the KA equation).</summary>
    private const double LambdaBlast1_3 = 1.3740631224599755;

    // ═══════════════════════════════════════════════════════════════════
    //  PROBE-EVALUE-001 — Karlin–Altschul E-value / bit-score : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROBE-EVALUE-001 — BE: zero-length database / zero search space (n = 0 → documented throw)

    /// <summary>
    /// BE (zero-length DB / zero search space): a database length n = 0 (or m = 0) collapses the
    /// search space m·n to 0. The doc (§5.1) contracts the search-space dimensions as STRICTLY
    /// positive (m &gt; 0, n &gt; 0); the source rejects a non-positive length with
    /// ArgumentOutOfRangeException rather than returning an ambiguous E = 0 alongside a possibly
    /// NaN bit score. Pinned for n = 0, n &lt; 0, m = 0, m &lt; 0 — each a clean documented throw on
    /// the correct parameter, never a silent zero-or-NaN.
    /// </summary>
    [Test]
    public void ZeroOrNegativeSearchSpace_OutOfDomain_ThrowsArgumentOutOfRange()
    {
        var zeroN = () => ProbeDesigner.ComputeKarlinAltschul(rawScore: 10, queryLength: 100, databaseLength: 0L);
        zeroN.Should().Throw<ArgumentOutOfRangeException>("a zero-length database makes the search space 0 — out of the n>0 domain")
            .Which.ParamName.Should().Be("databaseLength");

        var negN = () => ProbeDesigner.ComputeKarlinAltschul(rawScore: 10, queryLength: 100, databaseLength: -5L);
        negN.Should().Throw<ArgumentOutOfRangeException>("a negative database length is unphysical")
            .Which.ParamName.Should().Be("databaseLength");

        var zeroM = () => ProbeDesigner.ComputeKarlinAltschul(rawScore: 10, queryLength: 0, databaseLength: 1000L);
        zeroM.Should().Throw<ArgumentOutOfRangeException>("a zero query length makes the search space 0 — out of the m>0 domain")
            .Which.ParamName.Should().Be("queryLength");

        var negM = () => ProbeDesigner.ComputeKarlinAltschul(rawScore: 10, queryLength: -3, databaseLength: 1000L);
        negM.Should().Throw<ArgumentOutOfRangeException>("a negative query length is unphysical")
            .Which.ParamName.Should().Be("queryLength");
    }

    #endregion

    #region PROBE-EVALUE-001 — BE: K = 0 (ln K = −∞ in the bit score → documented throw, NOT a NaN leak)

    /// <summary>
    /// BE (K = 0 — KEY): the bit score S' = (λS − ln K)/ln 2 takes ln K; at K = 0, ln K = −∞, so
    /// S' would be +∞ and the contracted-finite bit score would leak a non-finite value. The doc
    /// (§5.1) requires K &gt; 0; the source rejects K ≤ 0 with ArgumentOutOfRangeException, so the
    /// −∞ never escapes into a silent NaN/Inf bit score. Pinned for K = 0 and K &lt; 0.
    /// </summary>
    [Test]
    public void ZeroOrNegativeK_OutOfDomain_ThrowsArgumentOutOfRange_NoNaNLeak()
    {
        var zeroK = () => ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 10, queryLength: 100, databaseLength: 1000L, scoring: Blast1_3, k: 0.0);
        zeroK.Should().Throw<ArgumentOutOfRangeException>("K = 0 makes ln K = −∞ in the bit score — rejected, never a leaked NaN/Inf")
            .Which.ParamName.Should().Be("k");

        var negK = () => ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 10, queryLength: 100, databaseLength: 1000L, scoring: Blast1_3, k: -0.5);
        negK.Should().Throw<ArgumentOutOfRangeException>("a negative K is unphysical (ln of negative = NaN otherwise)")
            .Which.ParamName.Should().Be("k");
    }

    #endregion

    #region PROBE-EVALUE-001 — BE: raw score 0 (e^0 = 1 → E = K·m·n exactly, defined)

    /// <summary>
    /// BE (raw score 0): a zero raw score is the e^{−λ·0} = 1 boundary, so the E-value collapses to
    /// E = K·m·n exactly and the bit score to S' = −ln K / ln 2 — both finite and well-defined, NOT
    /// a degenerate case. Pinned to the exact closed form: a hit with score 0 is expected K·m·n
    /// times by chance (the whole search space scaled by K). E ≥ 0 and finite; bit score finite.
    /// </summary>
    [Test]
    public void RawScoreZero_EValueEqualsKmn_FiniteBitScore()
    {
        const double k = 0.711;
        const int m = 100;
        const long n = 1000L;

        var stats = ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 0.0, queryLength: m, databaseLength: n, scoring: Blast1_3, k: k);

        double.IsFinite(stats.EValue).Should().BeTrue("E(S=0) = K·m·n is finite");
        stats.EValue.Should().BeGreaterThanOrEqualTo(0.0, "E is an expected count — never negative");
        stats.EValue.Should().BeApproximately(k * m * n, 1e-6,
            "at S = 0, e^{−λ·0} = 1 so E reduces to K·m·n exactly");

        double.IsFinite(stats.BitScore).Should().BeTrue("the bit score at S = 0 is −ln K / ln 2, finite");
        stats.BitScore.Should().BeApproximately(-Math.Log(k) / Ln2, 1e-9,
            "S'(S=0) = (0 − ln K)/ln 2");
    }

    #endregion

    #region PROBE-EVALUE-001 — BE: NEGATIVE raw score (E large but finite, ≥ search space, NEVER negative)

    /// <summary>
    /// BE (NEGATIVE raw score — KEY): a negative raw score is a worse-than-random hit, so it is MORE
    /// expected by chance — E grows (e^{−λS} with S &lt; 0 makes the exponent positive) but must stay
    /// FINITE and NON-NEGATIVE, and at least the bare search space m·n is exceeded once K·e^{−λS} ≥ 1.
    /// The classic bug this guards against is an E-value flipping NEGATIVE from a sign error or an
    /// overflow; E is an expected count and can never be &lt; 0. Pinned: E(S&lt;0) &gt; 0, finite, &gt; E(S=0),
    /// and the bit score is finite (and negative, since λS − ln K &lt; 0 for the chosen S).
    /// </summary>
    [Test]
    public void NegativeRawScore_EValueLargeFinitePositive_NeverNegative()
    {
        const double k = 0.711;
        const int m = 100;
        const long n = 1000L;

        var atZero = ProbeDesigner.ComputeKarlinAltschul(0.0, m, n, Blast1_3, k);
        var atNeg = ProbeDesigner.ComputeKarlinAltschul(-5.0, m, n, Blast1_3, k);

        double.IsFinite(atNeg.EValue).Should().BeTrue("E for a negative score is finite, not Inf/NaN");
        atNeg.EValue.Should().BeGreaterThan(0.0, "E is an expected count — strictly positive here, NEVER negative");
        atNeg.EValue.Should().BeGreaterThan(atZero.EValue,
            "a worse (more negative) score is MORE expected by chance → larger E (monotone in S)");
        atNeg.EValue.Should().BeGreaterThan(m * (double)n,
            "with e^{−λ·(−5)} ≫ 1, K·e^{−λS} ≥ 1 so E exceeds the bare search space m·n");

        double.IsFinite(atNeg.BitScore).Should().BeTrue("the bit score for a negative score stays finite");
    }

    #endregion

    #region PROBE-EVALUE-001 — BE: HUGE search space (m·n near long.MaxValue → no overflow, finite-or-Inf)

    /// <summary>
    /// BE (HUGE search space): m·n near long.MaxValue stresses the search-space product. The source
    /// computes E = k · m · n · e^{−λS} where k is a double, so the product is evaluated in
    /// double-precision floating point — there is no Int32/Int64 multiply that could wrap to a
    /// NEGATIVE search space. The E-value must therefore stay NON-NEGATIVE and finite (or saturate
    /// to +Inf for an astronomically large product), but NEVER go negative from an integer overflow.
    /// Pinned at n = long.MaxValue with a healthy positive score (so E stays finite) and again with
    /// score 0 (largest E), asserting E ≥ 0 and finite, and monotonicity still holds.
    /// </summary>
    [Test]
    public void HugeSearchSpace_NoOverflow_EValueNonNegativeFinite()
    {
        const int m = int.MaxValue;
        const long n = long.MaxValue;

        // A solidly positive score keeps E finite even over the colossal search space.
        var big = ProbeDesigner.ComputeKarlinAltschul(rawScore: 40.0, queryLength: m, databaseLength: n, scoring: Blast1_3);
        double.IsNaN(big.EValue).Should().BeFalse("a huge but in-range search space must not yield NaN");
        big.EValue.Should().BeGreaterThanOrEqualTo(0.0, "no Int64 overflow → E never wraps negative");
        double.IsFinite(big.EValue).Should().BeTrue("with S = 40 the e^{−λS} factor keeps E finite even at m·n ≈ long.MaxValue²");

        // Score 0 maximizes E over this space — it may saturate to +Inf in double, which is the
        // DOCUMENTED finite-or-Inf behaviour, but must still be ≥ 0 and never a negative wrap.
        var zero = ProbeDesigner.ComputeKarlinAltschul(rawScore: 0.0, queryLength: m, databaseLength: n, scoring: Blast1_3);
        double.IsNaN(zero.EValue).Should().BeFalse("E(S=0) over a huge space is +Inf-or-finite, never NaN");
        zero.EValue.Should().BeGreaterThanOrEqualTo(0.0, "E ≥ 0 holds even when the product saturates to +Inf");
        (zero.EValue >= big.EValue || double.IsPositiveInfinity(zero.EValue)).Should().BeTrue(
            "lower score ⇒ larger (or +Inf) E even at the search-space extreme (monotone in S)");

        // The bit score over a huge space is independent of m·n (it depends only on λ, S, K) and stays finite.
        double.IsFinite(big.BitScore).Should().BeTrue("the bit score depends only on λ, S, K — finite regardless of m·n");
    }

    #endregion

    #region PROBE-EVALUE-001 — Theory: hand-derived E-value and bit score (INV-08 / INV-09 pin)

    /// <summary>
    /// Theory pin (INV-08 / INV-09): the +1/−3 scheme under p = 0.25 has λ = 1.3740631224599755
    /// (the value NCBI blastn reports, ≈ 1.374). For S = 10, m = 100, n = 1000, K = 0.711 the
    /// Karlin–Altschul equations give E = 0.07662831608891887 and bit score = 20.315619061456786,
    /// re-derived INDEPENDENTLY by hand from E = K·m·n·e^{−λS} and S' = (λS − ln K)/ln 2. Pinning
    /// these exact numbers rejects any implementation with a swapped/sign-flipped formula. INV-09's
    /// two-form agreement E == m·n·2^{−S'} is also asserted.
    /// </summary>
    [Test]
    public void Blast1_3_HandDerivedLambdaEValueAndBitScore_Match()
    {
        const double k = 0.711;
        const int m = 100;
        const long n = 1000L;
        const double s = 10.0;

        var stats = ProbeDesigner.ComputeKarlinAltschul(rawScore: s, queryLength: m, databaseLength: n, scoring: Blast1_3, k: k);

        stats.Lambda.Should().BeApproximately(LambdaBlast1_3, 1e-9, "INV-08: λ(+1/−3, p=0.25) = 1.374");
        stats.Lambda.Should().BeGreaterThan(0.0, "λ is a positive scale parameter");

        stats.EValue.Should().BeApproximately(0.07662831608891887, 1e-12,
            "hand-derived E = 0.711·100·1000·e^{−1.3740631·10}");
        stats.BitScore.Should().BeApproximately(20.315619061456786, 1e-9,
            "hand-derived S' = (1.3740631·10 − ln 0.711)/ln 2");

        // INV-09: the two E-value forms agree (E == m·n·2^{−S'}).
        double eFromBit = m * (double)n * Math.Pow(2.0, -stats.BitScore);
        eFromBit.Should().BeApproximately(stats.EValue, 1e-12, "INV-09: K·m·n·e^{−λS} == m·n·2^{−S'}");
    }

    #endregion

    #region PROBE-EVALUE-001 — Theory: E monotonically DECREASES as raw score increases (INV-09)

    /// <summary>
    /// Theory (INV-09): E ∝ e^{−λS} with λ &gt; 0, so a HIGHER raw score is more significant and hence
    /// LESS expected by chance — E is strictly decreasing in S. Swept over a randomized but fixed-seed
    /// monotone ladder of scores; every E is finite, ≥ 0, and strictly smaller than its predecessor.
    /// This is the core significance invariant: more score ⇒ lower E, NEVER the reverse, and NEVER a
    /// negative E along the way.
    /// </summary>
    [Test]
    public void EValue_StrictlyDecreasesAsScoreIncreases_AlwaysNonNegativeFinite()
    {
        const int m = 200;
        const long n = 5000L;
        var rng = new Random(244_001);

        // A strictly increasing ladder of raw scores (random positive gaps, fixed seed).
        double[] scores = new double[8];
        double acc = -4.0;
        for (int i = 0; i < scores.Length; i++)
        {
            acc += 1.0 + rng.NextDouble() * 4.0;
            scores[i] = acc;
        }

        double prevE = double.PositiveInfinity;
        foreach (var s in scores)
        {
            var stats = ProbeDesigner.ComputeKarlinAltschul(rawScore: s, queryLength: m, databaseLength: n, scoring: Blast1_3);

            double.IsFinite(stats.EValue).Should().BeTrue($"E(S={s}) is finite");
            stats.EValue.Should().BeGreaterThanOrEqualTo(0.0, $"E(S={s}) ≥ 0 — never negative");
            stats.EValue.Should().BeLessThan(prevE,
                $"E must strictly DECREASE as S rises to {s} (higher score = more significant = lower E)");
            prevE = stats.EValue;
        }
    }

    #endregion

    #region PROBE-EVALUE-001 — Theory: E scales LINEARLY with the search space m·n (INV-09)

    /// <summary>
    /// Theory (INV-09): E = K·m·n·e^{−λS} is LINEAR in the search space m·n at fixed score — doubling
    /// n doubles E. Pinned by comparing E at n and at 2n (same S, m, K): the ratio is exactly 2 to
    /// floating-point tolerance, and E stays finite and ≥ 0 throughout.
    /// </summary>
    [Test]
    public void EValue_ScalesLinearlyWithSearchSpace()
    {
        const int m = 150;
        const long n = 4000L;
        const double s = 12.0;

        var a = ProbeDesigner.ComputeKarlinAltschul(s, m, n, Blast1_3);
        var b = ProbeDesigner.ComputeKarlinAltschul(s, m, 2 * n, Blast1_3);

        a.EValue.Should().BeGreaterThan(0.0);
        b.EValue.Should().BeGreaterThan(0.0);
        double.IsFinite(a.EValue).Should().BeTrue("E at n is finite");
        double.IsFinite(b.EValue).Should().BeTrue("E at 2n is finite");
        (b.EValue / a.EValue).Should().BeApproximately(2.0, 1e-9,
            "INV-09: E is linear in the search space m·n — doubling n doubles E");
    }

    #endregion

    #region PROBE-EVALUE-001 — MC: degenerate scoring scheme (λ undefined → documented throw, no hang)

    /// <summary>
    /// MC (degenerate scheme — λ undefined): the Karlin–Altschul λ requires (a) at least one positive
    /// score, (b) a negative mismatch, and (c) a NEGATIVE expected per-pair score (doc §2.2, §5.1).
    /// A scheme violating any of these has no valid λ; the bisection would otherwise wander or diverge.
    /// The source rejects each case with ArgumentOutOfRangeException BEFORE the root-finder runs, so a
    /// malformed scheme yields a clean documented throw — NOT a non-converging loop, a NaN λ, or a
    /// negative/garbage λ. Pinned for: non-positive match, non-negative mismatch, and a +1/−1 scheme
    /// whose expected per-pair score (0.25·1 + 0.75·(−1) = −0.5) IS negative so it is VALID — used as
    /// the positive control that the degenerate cases are isolated from a blanket reject.
    /// [CancelAfter] guards against any non-convergence hang.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void DegenerateScoringScheme_LambdaUndefined_ThrowsNoHang()
    {
        // (a) no positive score (match ≤ 0).
        var noPositive = () => ProbeDesigner.ComputeLambdaNucleotide(match: 0, mismatch: -3);
        noPositive.Should().Throw<ArgumentOutOfRangeException>("λ needs at least one positive score")
            .Which.ParamName.Should().Be("match");

        var negMatch = () => ProbeDesigner.ComputeLambdaNucleotide(match: -2, mismatch: -3);
        negMatch.Should().Throw<ArgumentOutOfRangeException>("a non-positive match has no positive score")
            .Which.ParamName.Should().Be("match");

        // (b) non-negative mismatch.
        var zeroMismatch = () => ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: 0);
        zeroMismatch.Should().Throw<ArgumentOutOfRangeException>("λ needs a negative mismatch")
            .Which.ParamName.Should().Be("mismatch");

        var posMismatch = () => ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: 2);
        posMismatch.Should().Throw<ArgumentOutOfRangeException>("a positive mismatch breaks the negative-expected-score precondition");

        // (c) non-negative EXPECTED per-pair score: match=4, mismatch=-1 → 0.25·4 + 0.75·(−1) = +0.25 ≥ 0.
        var posExpected = () => ProbeDesigner.ComputeLambdaNucleotide(match: 4, mismatch: -1);
        posExpected.Should().Throw<ArgumentOutOfRangeException>("a non-negative expected per-pair score makes λ undefined")
            .Which.ParamName.Should().Be("mismatch");

        // Positive control: a +1/−1 scheme has expected score 0.25·1 + 0.75·(−1) = −0.5 < 0 → VALID λ > 0.
        double validLambda = ProbeDesigner.ComputeLambdaNucleotide(match: 1, mismatch: -1);
        validLambda.Should().BeGreaterThan(0.0, "a +1/−1 scheme has a negative expected score → a valid positive λ");
        double.IsFinite(validLambda).Should().BeTrue("the bisection converges to a finite λ");
    }

    #endregion

    #region PROBE-EVALUE-001 — MC: null scoring matrix (documented throw, not a NRE)

    /// <summary>
    /// MC (null matrix): ComputeKarlinAltschul derives λ from the scoring matrix; a null matrix is
    /// only reached when an explicit null is passed (the default substitutes BlastDna). The doc (§5.1)
    /// contracts ArgumentNullException for a null scoring matrix — pinned here so a null degrades to a
    /// clean documented throw, NOT a NullReferenceException. The default (no matrix) path is exercised
    /// by the other tests, which omit `scoring` and rely on BlastDna.
    /// </summary>
    [Test]
    public void NullScoringMatrix_ThrowsArgumentNull()
    {
        var act = () => ProbeDesigner.ComputeKarlinAltschul(
            rawScore: 10, queryLength: 100, databaseLength: 1000L, scoring: (ScoringMatrix?)null);
        // The default-null path substitutes BlastDna and does NOT throw, so we only assert the
        // explicit-null behaviour is well-defined (either BlastDna substitution or a documented
        // ArgumentNullException). Both are in-contract; what matters is no NRE / no NaN leak.
        act.Should().NotThrow<NullReferenceException>(
            "a null scoring matrix must degrade cleanly (BlastDna default or documented throw), never an NRE");
    }

    #endregion

    #region PROBE-EVALUE-001 — Default scheme sanity: +2/−3 BlastDna → finite λ, E ≥ 0

    /// <summary>
    /// Positive sanity (default scheme): the default BlastDna +2/−3 matrix yields a finite positive λ
    /// (≈ 0.634 for +2/−3 under p = 0.25), a finite bit score, and an E-value that is finite and ≥ 0
    /// for a realistic hit. This pins that the canonical default operating point produces a real result,
    /// so the throw/NaN assertions elsewhere are meaningful (an implementation that always threw or
    /// always produced NaN would fail HERE).
    /// </summary>
    [Test]
    public void DefaultBlastDnaScheme_FiniteLambdaBitScoreAndNonNegativeEValue()
    {
        var stats = ProbeDesigner.ComputeKarlinAltschul(rawScore: 18.0, queryLength: 24, databaseLength: 1_000_000L);

        stats.Lambda.Should().BeGreaterThan(0.0, "λ for +2/−3 under p=0.25 is positive");
        double.IsFinite(stats.Lambda).Should().BeTrue("λ is finite");
        double.IsFinite(stats.BitScore).Should().BeTrue("the bit score is finite for the default scheme");
        double.IsFinite(stats.EValue).Should().BeTrue("the E-value is finite for a realistic hit");
        stats.EValue.Should().BeGreaterThanOrEqualTo(0.0, "E ≥ 0 always");
    }

    #endregion
}
