using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the K-mer area — the alignment-free K-mer EUCLIDEAN DISTANCE
/// between two sequences (KMER-DIST-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain parameter values to a
/// unit and asserts that the code NEVER fails in an undisciplined way: no hang or
/// infinite loop, no state corruption, no nonsense output, and no *unhandled*
/// runtime exception (DivideByZeroException, an out-of-range index, a NaN/Inf
/// distance, an OutOfMemoryException). Every input must resolve to EITHER a
/// well-defined, theory-correct result, OR a *documented, intentional* validation
/// exception (here: ArgumentOutOfRangeException for k ≤ 0). A raw runtime
/// exception, a hang, a NaN distance, or a distance out of the documented range
/// is a bug, not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-DIST-001 — k-mer Euclidean distance
/// Checklist: docs/checklists/03_FUZZING.md, row 158.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate/boundary input shapes called
///          out in the checklist row: identical sequences (the self-distance
///          floor), disjoint sequences (no shared k-mer — the worst-case word
///          overlap), empty/short sequences (the zero-vector / zero-window
///          boundary that could divide by zero), and different k values
///          (including k ≤ 0 and k > the sequence length).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The k-mer-distance contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// K-mer Euclidean distance is an alignment-free dissimilarity between two
/// sequences: each sequence is mapped to its NORMALIZED k-mer frequency vector
/// (count ÷ number of windows L − k + 1), and the distance is the L2 distance
/// between the two vectors taken over the UNION of observed k-mers, with an absent
/// k-mer contributing a 0 component (K-mer_Euclidean_Distance.md §2.2):
///   f_s(w) = count_s(w) / (L_s − k + 1)
///   d(x, y) = sqrt( Σ_{w ∈ W} ( f_x(w) − f_y(w) )² )
/// The API entry under test is
///   KmerAnalyzer.KmerDistance(string seq1, string seq2, int k)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines 212–233).
/// It builds each per-sequence frequency vector via GetKmerFrequencies (count ÷
/// sum of counts = L − k + 1), unions the key sets, sums the squared per-word
/// frequency differences, and returns the square root (KmerAnalyzer.cs lines 217–232).
///
/// THE DOCUMENTED METRIC INVARIANTS (K-mer_Euclidean_Distance.md §2.4) — every
/// positive test below pins these, because they are the load-bearing correctness
/// checks that distinguish a true Euclidean metric from a broken one:
///   • INV-01 IDENTITY:    d(x, x) = 0 — equal frequency vectors give a zero sum
///                         of squares ("identical sequences yield a distance of 0").
///   • INV-02 SYMMETRY:    d(x, y) = d(y, x) — (f_x − f_y)² = (f_y − f_x)².
///   • INV-03 NON-NEGATIVITY/RANGE: d(x, y) ≥ 0 and finite — a square root of a
///                         sum of squares is non-negative; it is also finite
///                         because every frequency lies in [0, 1].
///   • INV-04 DISJOINT:    two sequences each consisting of a SINGLE distinct k-mer
///                         (frequency 1) with DISJOINT word sets give d = √2 — the
///                         frequency vectors are (1,0) and (0,1).
///
/// THE KEY FUZZ CONCERN — the zero-window / zero-vector boundary. A frequency is
/// count ÷ total windows. At zero windows (empty/null sequence, or k > the
/// sequence length, or a sequence shorter than k) GetKmerFrequencies' explicit
/// `if (total == 0) return new Dictionary<string,double>()` guard (KmerAnalyzer.cs
/// lines 182–183) returns the EMPTY map — the ZERO vector — so KmerDistance never
/// divides by zero and never forms a NaN. The distance of an empty (zero) vector
/// against the other vector is then exactly the L2 NORM of that other vector, and
/// 0 when BOTH are empty (K-mer_Euclidean_Distance.md §3.3 ASM-02, §6.1: "Both
/// sequences empty → 0"; "One sequence shorter than k → distance = L2 norm of the
/// other's frequency vector"). These tests pin that boundary so an empty profile
/// can never produce a 0/0 → NaN, a DivideByZeroException, or a crash.
///
/// Documented parameter contract (K-mer_Euclidean_Distance.md §3.1, §3.3, §6.1;
/// KmerAnalyzer.cs lines 206–215):
///   • k ≤ 0 → ArgumentOutOfRangeException (nameof(k), "K must be positive.",
///     KmerAnalyzer.cs lines 214–215; §6.1 "k ≤ 0 → ArgumentOutOfRangeException").
///     This guard runs FIRST, BEFORE the empty/null sequences are inspected, so a
///     degenerate k throws even when both sequences are empty.
///   • null / empty sequence, or a sequence shorter than k → empty (zero) frequency
///     vector (ASM-02); the distance is the L2 norm of the other vector, 0 when both
///     are empty.
///   • k > sequence length → that sequence has zero windows → empty (zero) vector.
///   • inputs are upper-cased before counting, so the distance is case-insensitive
///     (ASM-01).
///   • no alphabet restriction — any character may form a k-mer (§3.3).
///
/// The four checklist BE targets map to these documented behaviours:
///   • identical   → d(x, x) = 0 exactly (INV-01), and case-insensitively (ASM-01).
///   • disjoint    → two single-distinct-k-mer sequences with NO shared word → √2
///                   (INV-04); the maximum-dissimilarity word-overlap case.
///   • empty       → empty/null/short input is the zero vector: distance = L2 norm
///                   of the other vector; both-empty → 0; NO DivideByZero, NO NaN.
///   • different k → k ≤ 0 throws (rejected); k > length / k > one side gives the
///                   zero vector on that side; varying k never crashes and always
///                   yields a finite, non-negative, symmetric distance.
/// A positive-sanity test pins the canonical worked example from the algorithm doc
/// (§7.1: "ATGTGTG" vs "CATGTG", k = 3 → √0.11 ≈ 0.33166247903553997), hand-derived
/// from the union {ATG, CAT, GTG, TGT} with frequencies f_x = (0.2, 0, 0.4, 0.4)
/// and f_y = (0.25, 0.25, 0.25, 0.25).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerDistanceFuzzTests
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

    /// <summary>
    /// Well-formed-result helper: asserts a distance is a properly-formed metric value —
    /// finite (no NaN/±Inf), inside the documented range (≥ 0; INV-03), and SYMMETRIC
    /// (d(x,y) = d(y,x); INV-02). Re-evaluating the swapped order also re-exercises the
    /// code under both argument orders, catching any order-dependent crash or asymmetry.
    /// </summary>
    private static void AssertWellFormed(string a, string b, int k)
    {
        double d = KmerAnalyzer.KmerDistance(a, b, k);
        double dSwapped = KmerAnalyzer.KmerDistance(b, a, k);

        double.IsNaN(d).Should().BeFalse("a k-mer distance must never be NaN (no 0/0 over an empty profile)");
        double.IsInfinity(d).Should().BeFalse("a k-mer distance must be finite — every frequency lies in [0, 1]");
        d.Should().BeGreaterThanOrEqualTo(0.0, "INV-03: a Euclidean distance is non-negative");
        dSwapped.Should().BeApproximately(d, 1e-12, "INV-02: the distance is symmetric, d(x, y) = d(y, x)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-DIST-001 — k-mer Euclidean distance : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-DIST-001 — k-mer Euclidean distance

    #region BE — Target: identical sequences (the self-distance floor, INV-01)

    /// <summary>
    /// BE / KEY: identical inputs are the self-distance floor. Equal frequency vectors
    /// give a zero sum of squares, so d(x, x) = 0 EXACTLY — not merely "small"
    /// (K-mer_Euclidean_Distance.md §2.4 INV-01, §6.1 "Identical sequences → 0"). We pin
    /// the EXACT 0.0 for several sequence shapes (heterogeneous, homopolymer, k = 1) so
    /// the identity edge cannot silently drift to a non-zero floor.
    /// </summary>
    [Test]
    public void KmerDistance_IdenticalSequences_IsExactlyZero()
    {
        KmerAnalyzer.KmerDistance("ACGTACGT", "ACGTACGT", 3).Should().Be(0.0,
            "INV-01: identical sequences have identical frequency vectors, so the distance is exactly 0");
        KmerAnalyzer.KmerDistance("AAAA", "AAAA", 2).Should().Be(0.0,
            "INV-01: a homopolymer is identical to itself, distance 0");
        KmerAnalyzer.KmerDistance("ACGTACGT", "ACGTACGT", 1).Should().Be(0.0,
            "INV-01: identity holds at the smallest legal k = 1 too");

        // The self-distance floor on fixed-seed random content at several k.
        string seq = RandomDna(256, seed: 158_001);
        foreach (int k in new[] { 1, 2, 4, 7 })
            KmerAnalyzer.KmerDistance(seq, seq, k).Should().Be(0.0,
                $"INV-01: a sequence is always at distance 0 from itself (k = {k})");
    }

    /// <summary>
    /// BE: identity is case-INSENSITIVE — the implementation upper-cases before counting
    /// (ASM-01, K-mer_Euclidean_Distance.md §3.3, §6.1 "Lower-case vs upper-case → same
    /// result"). A lower-case copy of a sequence is, after normalization, the SAME
    /// sequence, so its distance to the upper-case form is exactly 0. This pins that the
    /// identity floor is not accidentally broken by letter case.
    /// </summary>
    [Test]
    public void KmerDistance_LowercaseVsUppercase_IsExactlyZero()
    {
        KmerAnalyzer.KmerDistance("acgtacgt", "ACGTACGT", 3).Should().Be(0.0,
            "ASM-01: input is upper-cased before counting, so mixed case is the identical k-mer profile → distance 0");
        KmerAnalyzer.KmerDistance("AcGtAcGt", "ACGTACGT", 2).Should().Be(0.0,
            "ASM-01: case is normalized away before the frequency vectors are compared");
    }

    #endregion

    #region BE — Target: disjoint sequences (no shared k-mer → √2, INV-04)

    /// <summary>
    /// BE / KEY: two sequences sharing NO k-mer is the worst-case word overlap. When each
    /// is a single distinct k-mer (frequency 1) with disjoint word sets, the frequency
    /// vectors over the union are (1, 0) and (0, 1), so d = √((1−0)² + (0−1)²) = √2
    /// (K-mer_Euclidean_Distance.md §2.4 INV-04, §6.1 "Disjoint single-k-mer sequences →
    /// √2"). We pin the EXACT √2 and confirm symmetry — the canonical maximum-dissimilarity
    /// value for normalized single-k-mer profiles.
    /// </summary>
    [Test]
    public void KmerDistance_DisjointSingleKmerSequences_IsSqrtTwo()
    {
        // "AAA" k=3 → {AAA:1.0}; "GGG" k=3 → {GGG:1.0}; union {AAA, GGG}, vectors (1,0) & (0,1).
        KmerAnalyzer.KmerDistance("AAA", "GGG", 3).Should().BeApproximately(Math.Sqrt(2.0), 1e-12,
            "INV-04: disjoint single-k-mer frequency vectors (1,0) and (0,1) give √((1)²+(1)²) = √2");

        // Symmetry and well-formedness of the disjoint case.
        AssertWellFormed("AAA", "GGG", 3);

        // Disjoint with longer single-distinct-k-mer runs: "AAAA"/"CCCC" at k=2 are both
        // single-k-mer profiles ({AA:1.0}, {CC:1.0}) with no shared word → still √2.
        KmerAnalyzer.KmerDistance("AAAA", "CCCC", 2).Should().BeApproximately(Math.Sqrt(2.0), 1e-12,
            "INV-04: {AA:1.0} vs {CC:1.0} are disjoint single-k-mer profiles → √2");
    }

    /// <summary>
    /// BE: general disjoint sequences (no shared k-mer, but each with several distinct
    /// k-mers) are NOT at distance √2 in general — but they ARE bounded and well-formed.
    /// When the word sets are disjoint, the sum of squares is Σf_x² + Σf_y² (no cross
    /// terms cancel), which for two probability-like frequency vectors is ≤ 2, so the
    /// distance is in [0, √2]. We pin that a disjoint multi-k-mer pair stays finite,
    /// non-negative, ≤ √2, and symmetric — never out of range, never a crash.
    /// </summary>
    [Test]
    public void KmerDistance_DisjointMultiKmerSequences_IsInRangeAndSymmetric()
    {
        // "ACACAC" 2-mers → {AC, CA}; "GTGTGT" 2-mers → {GT, TG}; the word sets are disjoint.
        const string x = "ACACAC";
        const string y = "GTGTGT";
        const int k = 2;

        AssertWellFormed(x, y, k);

        double d = KmerAnalyzer.KmerDistance(x, y, k);
        d.Should().BeLessThanOrEqualTo(Math.Sqrt(2.0) + 1e-12,
            "disjoint frequency vectors have no cancelling cross terms; their L2 distance is bounded by √2");
        d.Should().BeGreaterThan(0.0, "the sequences share no k-mer, so they are strictly dissimilar (distance > 0)");
    }

    #endregion

    #region BE — Target: empty / short sequences (zero-vector boundary; no DivideByZero)

    /// <summary>
    /// BE / KEY: both sequences empty is the lower size boundary and the both-zero-vector
    /// case. Each empty/null input yields zero windows → GetKmerFrequencies' `total == 0`
    /// guard returns the EMPTY (zero) vector (KmerAnalyzer.cs lines 182–183), so the union
    /// is empty, the sum of squares is 0, and d = 0 — NOT a 0/0 → NaN and NOT a
    /// DivideByZeroException (K-mer_Euclidean_Distance.md §3.3 ASM-02, §6.1 "Both sequences
    /// empty → 0"). We pin no-throw, the exact 0.0, and the absence of NaN for empty/empty,
    /// null/null, and empty/null mixed.
    /// </summary>
    [Test]
    public void KmerDistance_BothEmptyOrNull_IsZeroNoDivideByZero()
    {
        var emptyEmpty = () => KmerAnalyzer.KmerDistance(string.Empty, string.Empty, 3);
        var nullNull = () => KmerAnalyzer.KmerDistance(null!, null!, 3);
        var emptyNull = () => KmerAnalyzer.KmerDistance(string.Empty, null!, 3);

        emptyEmpty.Should().NotThrow(
            "both empty → both zero vectors → empty union → sum of squares 0; the guard avoids any 0/0 division");
        nullNull.Should().NotThrow("null is treated as empty by the underlying count guard, not as an error");
        emptyNull.Should().NotThrow("an empty and a null input are both the zero vector");

        double d = KmerAnalyzer.KmerDistance(string.Empty, string.Empty, 3);
        d.Should().Be(0.0, "ASM-02: two empty (zero) frequency vectors are at distance 0");
        double.IsNaN(d).Should().BeFalse("the zero-window guard avoids 0/0; no NaN distance is produced");
        KmerAnalyzer.KmerDistance(null!, null!, 3).Should().Be(0.0, "two null (zero) vectors are at distance 0");
        KmerAnalyzer.KmerDistance(string.Empty, null!, 3).Should().Be(0.0, "empty vs null are both the zero vector → 0");
    }

    /// <summary>
    /// BE / KEY: ONE sequence empty (or shorter than k) against a non-empty one is the
    /// asymmetric zero-vector case. The empty side is the zero vector, so the distance
    /// equals the L2 NORM of the other sequence's frequency vector (K-mer_Euclidean_Distance.md
    /// §3.3 ASM-02, §6.1 "One sequence shorter than k → distance = L2 norm of the other's
    /// frequency vector"). We hand-derive that norm and pin it EXACTLY: for "AAAA", k = 2,
    /// the profile is {AA:1.0}, whose L2 norm is √(1²) = 1.0; and we confirm no NaN/throw
    /// and that the empty side is symmetric.
    /// </summary>
    [Test]
    public void KmerDistance_OneEmptyAgainstNonEmpty_IsL2NormOfOther()
    {
        // "AAAA" k=2 → {AA:1.0}; empty → {}. Distance = √((1.0 − 0)²) = 1.0.
        var act = () => KmerAnalyzer.KmerDistance(string.Empty, "AAAA", 2);
        act.Should().NotThrow("the empty side is the zero vector; the distance is the other vector's L2 norm, no 0/0");

        double d = KmerAnalyzer.KmerDistance(string.Empty, "AAAA", 2);
        d.Should().BeApproximately(1.0, 1e-12,
            "ASM-02: empty vs {AA:1.0} = the L2 norm of (1.0) = 1.0");
        double.IsNaN(d).Should().BeFalse("no NaN: the non-empty side is a well-formed frequency vector");

        // Symmetry of the one-empty case.
        KmerAnalyzer.KmerDistance("AAAA", string.Empty, 2).Should().BeApproximately(1.0, 1e-12,
            "INV-02: the empty-vs-non-empty distance is symmetric");

        // A multi-k-mer norm: "ACGTACGT" k=3 → frequencies (2/6,2/6,1/6,1/6); norm = √(Σ f²).
        const double total = 6.0;
        double expectedNorm = Math.Sqrt(
            Math.Pow(2.0 / total, 2) + Math.Pow(2.0 / total, 2) +
            Math.Pow(1.0 / total, 2) + Math.Pow(1.0 / total, 2));
        KmerAnalyzer.KmerDistance(string.Empty, "ACGTACGT", 3).Should().BeApproximately(expectedNorm, 1e-12,
            "ASM-02: empty vs a non-empty profile equals that profile's L2 norm exactly");
    }

    /// <summary>
    /// BE: a sequence SHORTER than k has zero windows (L − k + 1 ≤ 0), so it is also the
    /// zero vector — the same boundary as the empty string, reached via k > length rather
    /// than length 0 (K-mer_Euclidean_Distance.md §3.3 ASM-02). "AC" with k = 5 yields no
    /// window; against another too-short sequence the distance is 0, and against a valid
    /// one it is the valid side's norm. Neither may divide by zero or throw.
    /// </summary>
    [Test]
    public void KmerDistance_SequenceShorterThanK_IsZeroVector()
    {
        // Both sides shorter than k → both zero vectors → 0.
        var act = () => KmerAnalyzer.KmerDistance("AC", "GT", 5);
        act.Should().NotThrow("both sequences shorter than k have zero windows; both are the zero vector, no 0/0");
        KmerAnalyzer.KmerDistance("AC", "GT", 5).Should().Be(0.0,
            "ASM-02: two too-short sequences are both the empty (zero) vector → distance 0");

        // One side too short (k=5 > len 2) is the zero vector; distance = the other's norm.
        // "AAAAA" k=5 → {AAAAA:1.0}, norm 1.0.
        KmerAnalyzer.KmerDistance("AC", "AAAAA", 5).Should().BeApproximately(1.0, 1e-12,
            "ASM-02: the too-short side is the zero vector; distance = L2 norm of {AAAAA:1.0} = 1.0");
        double.IsNaN(KmerAnalyzer.KmerDistance("AC", "AAAAA", 5)).Should().BeFalse("no NaN on the short-input boundary");
    }

    #endregion

    #region BE — Target: different k (k ≤ 0 rejected; varying k never crashes)

    /// <summary>
    /// BE / KEY: k = 0 is the degenerate floor and a meaningless k-mer length — there is no
    /// "length-0 substring" whose frequency could be compared. KmerDistance validates k
    /// FIRST and rejects k ≤ 0 with ArgumentOutOfRangeException(nameof(k), "K must be
    /// positive.") (KmerAnalyzer.cs lines 214–215; K-mer_Euclidean_Distance.md §3.3, §6.1
    /// "k ≤ 0 → ArgumentOutOfRangeException"). We pin that k = 0 throws and carries the
    /// documented "k" parameter name, so a 0-length k-mer can never reach the frequency
    /// vectors.
    /// </summary>
    [Test]
    public void KmerDistance_KZero_ThrowsArgumentOutOfRange()
    {
        var act = () => KmerAnalyzer.KmerDistance("ACGTACGT", "ACGTACGT", 0);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a 0-length k-mer is meaningless; the contract rejects k <= 0 rather than defining a distance over 0-mers")
            .Which.ParamName.Should().Be("k");
    }

    /// <summary>
    /// BE: a negative k is below the floor too and must be rejected the same way — pinning
    /// that the rejection boundary is exactly k ≤ 0, not merely k == 0. The k check runs
    /// BEFORE the sequences are inspected, so a degenerate k throws even when BOTH inputs
    /// are empty (the empty short-circuit does NOT win over k validation here, unlike the
    /// CountKmers string path). We pin both the negative-k throw and the empty-input-with-
    /// degenerate-k throw.
    /// </summary>
    [Test]
    public void KmerDistance_NegativeKOrDegenerateKWithEmpty_Throws()
    {
        var negative = () => KmerAnalyzer.KmerDistance("ACGTACGT", "ACGTACGT", -3);
        negative.Should().Throw<ArgumentOutOfRangeException>(
                "a negative k-mer length is nonsensical; the contract rejects all k <= 0")
            .Which.ParamName.Should().Be("k");

        var emptyDegenerateK = () => KmerAnalyzer.KmerDistance(string.Empty, string.Empty, 0);
        emptyDegenerateK.Should().Throw<ArgumentOutOfRangeException>(
                "KmerDistance validates k FIRST (before inspecting the sequences), so a degenerate k = 0 throws even on empty input")
            .Which.ParamName.Should().Be("k");
    }

    /// <summary>
    /// BE: varying k across the full legal range (k = 1, intermediate, k = L, k > L) on a
    /// fixed pair of sequences must NEVER crash and must always yield a finite,
    /// non-negative, symmetric distance. At k = 1 the comparison is over single-base
    /// composition; as k grows shared windows vanish and the distance grows; at k > BOTH
    /// lengths both sides become the zero vector → distance 0 (the upper boundary). We
    /// pin well-formedness at every k and the exact 0 at the both-too-short ceiling.
    /// </summary>
    [Test]
    public void KmerDistance_VaryingK_AlwaysWellFormed()
    {
        const string x = "ACGTACGTAC"; // length 10
        const string y = "ACGTTTGCAC"; // length 10, shares a prefix then diverges

        foreach (int k in new[] { 1, 2, 3, 5, 10 })
            AssertWellFormed(x, y, k);

        // k = 11 > both lengths → both zero vectors → distance exactly 0 (no windows fit).
        var act = () => KmerAnalyzer.KmerDistance(x, y, 11);
        act.Should().NotThrow("k > both sequence lengths gives two zero vectors; the distance is 0, never a crash");
        KmerAnalyzer.KmerDistance(x, y, 11).Should().Be(0.0,
            "at k > both lengths neither sequence has a window; both are the zero vector → distance 0");
    }

    /// <summary>
    /// BE / RB: a fixed-seed random PAIR of sequences (different lengths) must complete
    /// promptly and satisfy the full metric contract for several k values — finite,
    /// non-negative, in range [0, √2] (normalized frequency vectors), symmetric, and 0
    /// against itself — regardless of the random content. [CancelAfter] guards against any
    /// hang on the largest k scanned.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void KmerDistance_RandomPair_SatisfiesMetricContractForEveryK()
    {
        string a = RandomDna(1500, seed: 158_101);
        string b = RandomDna(2000, seed: 158_202); // deliberately different length

        foreach (int k in new[] { 1, 2, 3, 5, 8, 13 })
        {
            AssertWellFormed(a, b, k);

            double d = KmerAnalyzer.KmerDistance(a, b, k);
            d.Should().BeLessThanOrEqualTo(Math.Sqrt(2.0) + 1e-9,
                $"at k = {k} the distance between two normalized frequency vectors is bounded by √2");
            KmerAnalyzer.KmerDistance(a, a, k).Should().Be(0.0, $"INV-01: a is at distance 0 from itself (k = {k})");
        }
    }

    #endregion

    #region Positive sanity — the canonical worked example (hand-computed)

    /// <summary>
    /// Positive sanity / KEY: the canonical worked example from the algorithm doc
    /// (K-mer_Euclidean_Distance.md §7.1, Zielezinski et al. 2017 Fig. 1). For
    /// x = "ATGTGTG" (L = 7, 5 windows) and y = "CATGTG" (L = 6, 4 windows), k = 3, the
    /// union of words is {ATG, CAT, GTG, TGT}: c_x = (1, 0, 2, 2), c_y = (1, 1, 1, 1),
    /// giving f_x = (0.2, 0, 0.4, 0.4) and f_y = (0.25, 0.25, 0.25, 0.25). The squared
    /// differences (0.0025, 0.0625, 0.0225, 0.0225) sum to 0.11, so
    /// d = √0.11 ≈ 0.33166247903553997. This pins the EXACT distance from a fully
    /// hand-derived frequency computation — the load-bearing check that the boundary
    /// hardening never comes at the cost of the core metric silently breaking. A
    /// count-based (un-normalized) or wrong-divisor implementation would NOT reproduce
    /// √0.11, so this is a genuine discriminating value, not a rubber stamp.
    /// </summary>
    [Test]
    public void KmerDistance_CanonicalWorkedExample_MatchesHandComputedDistance()
    {
        double d = KmerAnalyzer.KmerDistance("ATGTGTG", "CATGTG", 3);

        d.Should().BeApproximately(Math.Sqrt(0.11), 1e-12,
            "the doc's Fig. 1 walk-through: Σ squared frequency diffs = 0.11 → d = √0.11 ≈ 0.33166247903553997");
        d.Should().BeApproximately(0.33166247903553997, 1e-12, "the documented numeric value (§7.1)");

        // Symmetry of the worked example.
        KmerAnalyzer.KmerDistance("CATGTG", "ATGTGTG", 3).Should().BeApproximately(d, 1e-12,
            "INV-02: the worked example is symmetric");
    }

    /// <summary>
    /// Positive sanity: the triangle-inequality direction of the metric on a known triple.
    /// A Euclidean distance is a true metric (K-mer_Euclidean_Distance.md §2.5 "Metric:
    /// yes"), so d(a, c) ≤ d(a, b) + d(b, c). We pin this on three concrete sequences as
    /// an extra structural check that the frequency-vector L2 distance behaves as a metric,
    /// complementing the per-test identity/symmetry/range assertions.
    /// </summary>
    [Test]
    public void KmerDistance_TriangleInequality_HoldsOnKnownTriple()
    {
        const string a = "ACGTACGT";
        const string b = "ACGTTTGC";
        const string c = "TTTTGCGC";
        const int k = 3;

        double dac = KmerAnalyzer.KmerDistance(a, c, k);
        double dab = KmerAnalyzer.KmerDistance(a, b, k);
        double dbc = KmerAnalyzer.KmerDistance(b, c, k);

        dac.Should().BeLessThanOrEqualTo(dab + dbc + 1e-12,
            "the frequency-vector Euclidean distance is a metric: d(a, c) ≤ d(a, b) + d(b, c)");
    }

    #endregion

    #endregion
}
