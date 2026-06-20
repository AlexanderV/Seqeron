using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Phylogenetics;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Phylogenetic area — phylogenetic distance calculation
/// (PHYLO-DIST-001): the pairwise evolutionary-distance routine and the symmetric
/// distance matrix built from it.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed, out-of-domain and boundary inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang or infinite loop, no
/// nonsense output, and no *unhandled* runtime exception (DivideByZeroException, a NaN
/// or ±Infinity leaking from log-of-non-positive, IndexOutOfRangeException). Every input
/// must resolve to EITHER a well-defined, theory-correct result OR a *documented,
/// intentional* validation exception (ArgumentException / ArgumentNullException). For a
/// distance function the central theory contract is that it behaves like a (pseudo)metric:
/// d(a,a) = 0 (identity of indiscernibles), d(a,b) = d(b,a) (symmetry), and d ≥ 0
/// (non-negativity). A raw runtime exception, a hang, a NaN, or a negative distance on
/// garbage input is a bug, not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8
/// "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-DIST-001 — phylogenetic distance
/// Checklist: docs/checklists/03_FUZZING.md, row 39.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the empty sequence (no comparable sites: a
///          division-by-zero hazard on the p = differences/comparableSites step), the
///          single-sequence matrix (a degenerate 1×1 matrix with only the 0 diagonal),
///          and the saturation boundary p ≥ 0.75 where the Jukes-Cantor / Kimura
///          logarithm argument turns non-positive (log of ≤ 0 → NaN/−Inf hazard).
///   • MC = Malformed Content — non-DNA / gap / ambiguous characters in otherwise
///          aligned sequences, which the implementation skips by pairwise deletion.
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, empty;
///   MC = malformed content); row 39 targets:
///   "Identical seqs, empty seqs, single seq, non-DNA chars".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The distance contract under test (Distance_Matrix.md)
/// ───────────────────────────────────────────────────────────────────────────
/// PhylogeneticAnalyzer exposes ONE shared pairwise scan over aligned sequences behind
/// four distance models selected by the DistanceMethod enum (Distance_Matrix.md §1, §2):
///   • Hamming  — raw mismatch COUNT over comparable sites:  d_H = Σ 1[s1[i] ≠ s2[i]].
///   • PDistance — proportion of differing comparable sites:  p = differences / comparable.
///   • JukesCantor (JC69) — d = −3/4 · ln(1 − 4p/3)            (Distance_Matrix.md §2.C).
///   • Kimura2Parameter (K80) — d = −1/2 · ln((1 − 2S − V)·√(1 − 2V))  (§2.D),
///       S = transitions/comparable, V = transversions/comparable.
///
/// API entry points (Distance_Matrix.md §5.1; PhylogeneticAnalyzer.cs):
///   • CalculatePairwiseDistance(string seq1, string seq2, DistanceMethod = JukesCantor)
///     (lines 223–270) — one scalar distance. null seq → ArgumentNullException;
///     unequal lengths → ArgumentException (lines 226–229).
///   • CalculateDistanceMatrix(IReadOnlyList<string>, DistanceMethod = JukesCantor)
///     (lines 199–218) — a symmetric n×n double[,] with a 0 diagonal, filled only for
///     i &lt; j and mirrored (j,i); the diagonal is left at its default 0.0.
///
/// THE FOUR ROW-39 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Identical seqs (KEY — identity of indiscernibles): two equal sequences have zero
///     differences, so p = 0 and EVERY method returns exactly 0 — Hamming 0, PDistance 0,
///     JC −3/4·ln(1) = 0, K2P −1/2·ln(1·√1) = 0 (Distance_Matrix.md §6.1 "Identical
///     comparable sequences → distance 0 for all methods"). Pinned as the metric's
///     defining property, and as the zero diagonal of the matrix.
///   • Empty seqs (BE — div-by-zero hazard): two empty (or all-gap / all-ambiguous)
///     sequences have comparableSites = 0. The implementation GUARDS this with an explicit
///     `if (comparableSites == 0) return 0;` (PhylogeneticAnalyzer.cs line 256) BEFORE the
///     p = differences/comparableSites division (line 258), so the denominator is never
///     zero — the theory-correct boundary is a finite 0, NOT a DivideByZeroException, NaN,
///     or Infinity (Distance_Matrix.md §3.3, §6.1 "no comparable sites → 0"). This is the
///     central BE probe.
///   • Single seq (BE — degenerate matrix): a 1-element sequence list yields a 1×1 matrix.
///     The fill loop `for j = i+1` never runs, so only the default-0.0 diagonal remains —
///     a well-formed 1×1 matrix whose single entry [0,0] is 0, no throw, no out-of-range
///     indexing (CalculateDistanceMatrix lines 206–215). The pairwise scalar surface
///     instead requires TWO sequences and is exercised via self-distance d(a,a) = 0.
///   • Non-DNA chars (MC — pairwise deletion): gaps ('-') and any non-A/C/G/T symbol are
///     SKIPPED at a site rather than crashing or being miscounted (lines 242–243,
///     IsStandardBase). Junk therefore lowers the comparable-site count but never throws,
///     never invents a difference, and never produces a NaN; an all-junk pair collapses to
///     the comparableSites = 0 → 0 boundary above (Distance_Matrix.md §3.3, §5.2). Case is
///     irrelevant: each site is upper-cased before inspection (line 238–239).
///
/// THE SATURATION BOUNDARY (BE — log of non-positive): for a CORRECTED model the hazard is
/// p → 0.75, where 1 − 4p/3 → 0 and the JC logarithm argument turns non-positive. The
/// source GUARDS this with `if (arg <= 0) return double.PositiveInfinity;` (line 287; the
/// K2P helper guards both of its arguments, lines 296). So a fully-saturated pair returns a
/// DEFINED +Infinity sentinel — NEVER a NaN and never a −Infinity from ln(negative)
/// (Distance_Matrix.md §6.1 "p ≥ 0.75 in JC69 → positive infinity"; INV-JC-02, INV-K2P-02).
/// +Infinity is a finite-arithmetic-safe, non-negative, intentional saturation marker, so it
/// satisfies non-negativity; we pin that it is +Infinity and specifically NOT NaN.
///
/// Documented invariants pinned (Distance_Matrix.md §2): INV-HAMMING-02 / INV-PDIST /
/// identical → 0; INV-PDIST-01 0 ≤ p ≤ 1; INV-JC-01 d_JC ≥ p; INV-JC-02 / INV-K2P-02
/// saturation → +Infinity (not NaN). The pinned exact value: for p = 1/4 (one mismatch in
/// four comparable A/C/G/T sites) JC = −3/4·ln(1 − 1/3) = −3/4·ln(2/3) ≈ 0.3040988.
/// CalculatePairwiseDistance and CalculateDistanceMatrix are pure (no iterators), so every
/// probe calls them directly; deterministic fuzz inputs use a locally fixed-seed Random
/// (no shared static RNG).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PhylogeneticFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    private static readonly PhylogeneticAnalyzer.DistanceMethod[] AllMethods =
    {
        PhylogeneticAnalyzer.DistanceMethod.Hamming,
        PhylogeneticAnalyzer.DistanceMethod.PDistance,
        PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
        PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter,
    };

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PHYLO-DIST-001 — phylogenetic distance : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PHYLO-DIST-001 — phylogenetic distance

    #region Positive sanity — known distance + metric properties

    /// <summary>
    /// Positive-sanity anchor: a hand-checkable pair pins the EXACT documented distance for
    /// every model, and the three metric properties (identical → 0, symmetry, non-negativity).
    /// seq1 = "AAAA", seq2 = "AAAG": one mismatch (A→G, a transition) in four comparable A/C/G/T
    /// sites ⇒ differences = 1, comparable = 4, p = 1/4, transitions = 1, transversions = 0.
    ///   • Hamming  = 1 (raw count).
    ///   • PDistance = 0.25.
    ///   • JC69     = −3/4·ln(1 − 4·0.25/3) = −3/4·ln(2/3) ≈ 0.3040988.
    ///   • K2P      = −1/2·ln((1 − 2·0.25 − 0)·√(1 − 0)) = −1/2·ln(0.5) ≈ 0.3465736.
    /// </summary>
    [Test]
    public void PairwiseDistance_KnownPair_MatchesDocumentedFormula_AndIsAMetric()
    {
        const string a = "AAAA";
        const string b = "AAAG";

        double hamming = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Hamming);
        double pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.PDistance);
        double jc = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
        double k2p = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        hamming.Should().Be(1.0, "one differing site over the comparable A/C/G/T positions (INV-HAMMING-01)");
        pDist.Should().BeApproximately(0.25, 1e-12, "p = differences/comparable = 1/4 (Distance_Matrix.md §2.B)");
        jc.Should().BeApproximately(0.3040988, 1e-6, "JC69 = −3/4·ln(1 − 4p/3) at p = 1/4 (Distance_Matrix.md §2.C)");
        k2p.Should().BeApproximately(0.3465736, 1e-6, "K2P = −1/2·ln((1 − 2S − V)·√(1 − 2V)) at S = 1/4, V = 0 (§2.D)");

        // INV-JC-01: the corrected distance never under-estimates the raw proportion.
        jc.Should().BeGreaterThanOrEqualTo(pDist, "JC69 corrects upward for hidden substitutions (INV-JC-01)");

        foreach (var method in AllMethods)
        {
            double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);

            // Symmetry: d(a,b) = d(b,a).
            double reversed = PhylogeneticAnalyzer.CalculatePairwiseDistance(b, a, method);
            reversed.Should().Be(d, "phylogenetic distance is symmetric: d(a,b) = d(b,a)");

            // Non-negativity and finiteness for this non-saturated pair.
            d.Should().BeGreaterThanOrEqualTo(0.0, "an evolutionary distance is non-negative");
            double.IsNaN(d).Should().BeFalse("a well-defined pair must not produce NaN");

            // Identity of indiscernibles: d(a,a) = 0.
            double self = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, a, method);
            self.Should().Be(0.0, "identical sequences have distance 0 (identity of indiscernibles)");
        }
    }

    #endregion

    #region BE — Identical sequences (identity of indiscernibles)

    /// <summary>
    /// Identical aligned sequences ⇒ zero differences ⇒ p = 0 ⇒ distance 0 for EVERY method,
    /// both as the scalar self-distance and as the zero diagonal of the matrix
    /// (Distance_Matrix.md §6.1; INV-HAMMING-02). Uses a deterministic, locally-seeded random
    /// sequence so the identity holds for arbitrary content, not a hand-picked string.
    /// </summary>
    [Test]
    public void IdenticalSequences_AreDistanceZero_ForAllMethods()
    {
        string s = RandomDna(40, seed: 39_001);

        foreach (var method in AllMethods)
        {
            double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(s, s, method);
            d.Should().Be(0.0, "identical comparable sequences differ at no site → distance 0 ({0})", method);
        }

        // The matrix diagonal is the self-distance: an identical/equal-content matrix has a 0 diagonal.
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(new[] { s, s, s }, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
        for (int i = 0; i < 3; i++)
            matrix[i, i].Should().Be(0.0, "the distance-matrix diagonal is the self-distance, which is 0");
        // All off-diagonal pairs are identical content ⇒ also 0.
        matrix[0, 1].Should().Be(0.0, "two equal sequences have distance 0 off the diagonal too");
        matrix[1, 0].Should().Be(0.0, "the matrix is symmetric");
    }

    #endregion

    #region BE — Empty sequences (no comparable sites → div-by-zero hazard)

    /// <summary>
    /// Two empty sequences (and the all-gap / all-ambiguous degenerate that also leaves zero
    /// comparable sites) hit the `comparableSites == 0` guard and return a finite 0 for every
    /// method — NEVER a DivideByZeroException, NaN, or Infinity from the p = differences/0 step
    /// (Distance_Matrix.md §3.3, §6.1). This is the central boundary probe.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void EmptySequences_NoComparableSites_ReturnZeroNeverNaN_ForAllMethods()
    {
        foreach (var method in AllMethods)
        {
            // Both empty: length 0 ⇒ no sites ⇒ comparableSites == 0 ⇒ guarded 0.
            Action act = () => PhylogeneticAnalyzer.CalculatePairwiseDistance(string.Empty, string.Empty, method);
            act.Should().NotThrow("zero comparable sites are guarded before the division ({0})", method);

            double empty = PhylogeneticAnalyzer.CalculatePairwiseDistance(string.Empty, string.Empty, method);
            empty.Should().Be(0.0, "no comparable sites → distance 0 (no division by zero)");
            double.IsNaN(empty).Should().BeFalse("the zero-site boundary must not produce NaN");
            double.IsInfinity(empty).Should().BeFalse("the zero-site boundary must not produce Infinity");

            // All-gap pair: every site is skipped by pairwise deletion ⇒ also zero comparable sites.
            double allGap = PhylogeneticAnalyzer.CalculatePairwiseDistance("----", "----", method);
            allGap.Should().Be(0.0, "an all-gap pair leaves no comparable site → guarded 0");
        }
    }

    #endregion

    #region BE — Single sequence (degenerate 1×1 matrix)

    /// <summary>
    /// A single-element sequence list yields a well-formed 1×1 matrix whose only entry is the
    /// default-0.0 diagonal: the i&lt;j fill loop never runs, so no division and no out-of-range
    /// indexing occur (CalculateDistanceMatrix). The scalar surface instead requires two
    /// sequences, so a "single sequence" is exercised there as the self-distance d(a,a) = 0.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void SingleSequence_MatrixIsOneByOneZero_AndSelfDistanceIsZero()
    {
        string s = RandomDna(25, seed: 39_002);

        foreach (var method in AllMethods)
        {
            Action act = () => PhylogeneticAnalyzer.CalculateDistanceMatrix(new[] { s }, method);
            act.Should().NotThrow("a 1×1 matrix has no pair to compare; the fill loop never runs ({0})", method);

            var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(new[] { s }, method);
            matrix.GetLength(0).Should().Be(1, "one taxon → a 1×1 matrix");
            matrix.GetLength(1).Should().Be(1);
            matrix[0, 0].Should().Be(0.0, "the lone diagonal entry is the self-distance, 0");

            // The single-sequence notion on the scalar surface: distance to itself is 0.
            double self = PhylogeneticAnalyzer.CalculatePairwiseDistance(s, s, method);
            self.Should().Be(0.0, "a single sequence compared to itself has distance 0");
        }
    }

    #endregion

    #region MC — Non-DNA / gap / ambiguous characters (pairwise deletion)

    /// <summary>
    /// Gaps and non-A/C/G/T symbols are SKIPPED at a site (pairwise deletion), never crash and
    /// never invent a difference: the distance equals that of the cleaned, comparable-only
    /// sub-alignment. Here only positions 0,1,2 are comparable A/C/G/T (one mismatch at index 2),
    /// while positions 3+ carry junk in at least one row and are excluded ⇒ p = 1/3 regardless of
    /// the junk (Distance_Matrix.md §3.3, §5.2). Case-insensitivity is pinned alongside.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void NonDnaCharacters_AreSkipped_NoCrashNoNaN_NoSpuriousDifference()
    {
        //               idx: 0123456
        const string a = "ACG-N*x";
        const string b = "ACT9 Zq";
        // Comparable A/C/G/T-vs-A/C/G/T positions: 0 (A=A), 1 (C=C), 2 (G≠T). Index 3+ excluded.
        // ⇒ differences = 1, comparable = 3, p = 1/3.

        foreach (var method in AllMethods)
        {
            Action act = () => PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);
            act.Should().NotThrow("non-DNA/gap symbols are skipped, never crash ({0})", method);

            double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);
            double.IsNaN(d).Should().BeFalse("pairwise-deleted junk must not produce NaN ({0})", method);
            d.Should().BeGreaterThanOrEqualTo(0.0, "distance stays non-negative under pairwise deletion");
        }

        // Exact pins on the uncorrected models prove junk is excluded, not miscounted.
        PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Hamming)
            .Should().Be(1.0, "exactly one comparable mismatch (index 2); all junk positions are deleted");
        PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.PDistance)
            .Should().BeApproximately(1.0 / 3.0, 1e-12, "p = 1 difference / 3 comparable sites");

        // An all-junk pair collapses to the zero-comparable-sites boundary → 0, not a crash.
        PhylogeneticAnalyzer.CalculatePairwiseDistance("****", "----", PhylogeneticAnalyzer.DistanceMethod.JukesCantor)
            .Should().Be(0.0, "no comparable site survives pairwise deletion → guarded 0");

        // Case-insensitive: lower-case bases compare equal to their upper-case counterparts.
        PhylogeneticAnalyzer.CalculatePairwiseDistance("acgt", "ACGT", PhylogeneticAnalyzer.DistanceMethod.Hamming)
            .Should().Be(0.0, "each site is upper-cased before inspection (case-insensitive)");
    }

    #endregion

    #region BE — Saturation boundary (log of non-positive → +Infinity, not NaN)

    /// <summary>
    /// The corrected-model saturation boundary: when p ≥ 0.75 the JC argument 1 − 4p/3 ≤ 0 and
    /// the K2P argument 1 − 2S − V ≤ 0. The source guards both and returns +Infinity rather than
    /// taking ln of a non-positive value (Distance_Matrix.md §6.1; INV-JC-02, INV-K2P-02). The
    /// theory-correct sentinel is +Infinity — specifically NOT NaN and NOT −Infinity — while the
    /// uncorrected models stay finite (PDistance saturates at 1, Hamming counts).
    /// </summary>
    [Test]
    public void SaturatedPair_CorrectedModelsReturnPositiveInfinity_NeverNaN()
    {
        // All four sites differ as transitions ⇒ p = 1.0 (≥ 0.75), S = 1.0, V = 0.0.
        const string a = "AAAA";
        const string b = "GGGG";

        double jc = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
        double k2p = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        double.IsPositiveInfinity(jc).Should().BeTrue("JC69 with p ≥ 0.75 has a non-positive log argument → +Infinity (INV-JC-02)");
        double.IsNaN(jc).Should().BeFalse("the saturation guard must return +Infinity, never NaN from ln(≤0)");
        double.IsPositiveInfinity(k2p).Should().BeTrue("K2P with a non-positive log argument → +Infinity (INV-K2P-02)");
        double.IsNaN(k2p).Should().BeFalse("the K2P saturation guard must return +Infinity, never NaN");

        // The uncorrected models stay finite at saturation.
        double p = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.PDistance);
        p.Should().Be(1.0, "every site differs → p saturates at 1, finite (INV-PDIST-01 upper bound)");
        PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Hamming)
            .Should().Be(4.0, "all four comparable sites differ");
    }

    #endregion

    #region Validation — null and unequal-length inputs (documented exceptions)

    /// <summary>
    /// The two documented validation throws on the scalar surface: a null sequence →
    /// ArgumentNullException; unequal lengths (an unaligned pair) → ArgumentException
    /// (PhylogeneticAnalyzer.cs lines 226–229; Distance_Matrix.md §6.1). These are
    /// INTENTIONAL, contract-defined rejections — not raw runtime crashes.
    /// </summary>
    [Test]
    public void NullOrUnequalLength_ThrowDocumentedValidationExceptions()
    {
        Action nullFirst = () => PhylogeneticAnalyzer.CalculatePairwiseDistance(null!, "ACGT");
        Action nullSecond = () => PhylogeneticAnalyzer.CalculatePairwiseDistance("ACGT", null!);
        nullFirst.Should().Throw<ArgumentNullException>("a null sequence is an explicit, documented rejection");
        nullSecond.Should().Throw<ArgumentNullException>("a null sequence is an explicit, documented rejection");

        Action unequal = () => PhylogeneticAnalyzer.CalculatePairwiseDistance("ACG", "ACGT");
        unequal.Should().Throw<ArgumentException>("pairwise distance requires aligned (equal-length) sequences");
    }

    #endregion

    #endregion
}
