// RNA-PSEUDOKNOT-001 — RNA pseudoknot (crossing base-pair) detection.
// Fuzz tests (strategy BE = Boundary Exploitation).
// Algorithm doc: docs/algorithms/RnaStructure/Pseudoknot_Detection.md
// Canonical tests: tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_DetectPseudoknots_Tests.cs (RNA-PSEUDOKNOT-001)
// Evidence: docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md
// Source: RnaSecondaryStructure.DetectPseudoknots(IReadOnlyList<BasePair>) — RnaSecondaryStructure.cs.
//         Antczak et al. (2018) Bioinformatics 34(8):1304–1312; Smit et al. (2008) RNA 14(3):410;
//         biotite.structure.pseudoknots; Rivas & Eddy (1999).

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for RNA-PSEUDOKNOT-001 —
/// <see cref="RnaSecondaryStructure.DetectPseudoknots(IReadOnlyList{BasePair})"/>, the exact,
/// deterministic combinatorial test that reports every <em>crossing</em> pair-of-base-pairs
/// (a pseudoknot) in a fixed base-pair set. No thermodynamics, no scoring — the output depends
/// only on nucleotide positions.
/// Lives in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Documented contract (Pseudoknot_Detection.md §2.2, §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   • Two base pairs, written with open &lt; close as (i,j) and (k,l), CROSS (form a pseudoknot)
///       IFF  i &lt; k &lt; j &lt; l  (§2.2 Core Model; verbatim Antczak 2018: ∃ (j,j′) with
///       i &lt; j &lt; i′ &lt; j′). All four inequalities are STRICT.
///   • NESTED  (i &lt; k &lt; l &lt; j, one pair fully inside the other) → NOT a pseudoknot (INV-02,
///       fails j &lt; l).
///   • DISJOINT (j &lt; k, non-overlapping ranges) → NOT a pseudoknot (INV-03, fails k &lt; j).
///   • Each pair is normalized to (open &lt; close) via min/max before comparison, so a pair stored
///       as (close,open) yields the IDENTICAL result (§3.3, §6.1 "pair stored as (close,open)").
///   • The two pairs are reordered by opening position, so the result is independent of input
///       ordering and storage direction (INV-05: deterministic & order-independent).
///   • &lt; 2 pairs, or a `null` set, → EMPTY result, NO exception (INV-04, §3.3, §6.1).
///   • A degenerate pair (i = j) cannot cross anything — silently never part of a pseudoknot
///       (ASM-01): every reported crossing satisfies strict i &lt; k &lt; j &lt; l (INV-01).
///   • Output (§3.2): one <see cref="Pseudoknot"/> per crossing pair-of-pairs, with
///       Start1 &lt; End1, Start2 &lt; End2, Start1 &lt; Start2, and Start1 &lt; Start2 &lt; End1 &lt; End2;
///       `CrossingPairs` holds the two ORIGINAL `BasePair`s.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"
/// ───────────────────────────────────────────────────────────────────────────
/// Feed degenerate / boundary base-pair sets and assert the detector NEVER fails undisciplined:
/// no exception (no IndexOutOfRange, no DivideByZero on empty), no FALSE pseudoknot on a nested
/// or merely-touching structure (off-by-one in the i&lt;k&lt;j&lt;l comparison), no MISSED crossing.
/// Every result is well formed: every reported knot really satisfies strict i&lt;k&lt;j&lt;l, and the
/// reported count matches an independent O(n²) oracle.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy BE = Boundary Exploitation — row 155
/// targets "nested-only, fully crossing, empty"
/// ───────────────────────────────────────────────────────────────────────────
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, -1, MaxInt, empty), row 155.
///
///   • nested-only (BE) — a fully nested structure like (((...))) → NO pseudoknot; the canonical
///       negative. Boundary: pairs that share/touch an endpoint (i&lt;k&lt;j=l, or i=k&lt;j&lt;l) must
///       NOT cross because the inequalities are strict.
///   • fully crossing (BE) — the classic H-type pseudoknot "([)]" / two interleaved stems →
///       pseudoknot DETECTED with the documented crossing pairs; the canonical positive.
///   • empty (BE) — empty set / single pair / `null` → NO pseudoknot, NO crash, NO DivideByZero.
///
/// Watched failure modes: false pseudoknot on a nested or touching structure (≤ vs &lt; off-by-one);
/// missed crossing; crash on empty/null; IndexOutOfRange; non-determinism / order dependence.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaPseudoknotFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static Random Rng(int seed) => new(seed);

    /// <summary>Builds a base pair from two positions (bases are irrelevant to the positional test).</summary>
    private static BasePair Pair(int p1, int p2) =>
        new(p1, p2, 'A', 'U', BasePairType.WatsonCrick);

    /// <summary>
    /// Independent O(n²) oracle for the documented crossing relation: counts the unordered
    /// pairs-of-pairs (a,b) with a&lt;b that satisfy strict i&lt;k&lt;j&lt;l after (open&lt;close) and
    /// open-first normalization — exactly Pseudoknot_Detection.md §2.2 / §4.1.
    /// </summary>
    private static int ExpectedCrossings(IReadOnlyList<BasePair> pairs)
    {
        int count = 0;
        for (int a = 0; a < pairs.Count; a++)
        for (int b = a + 1; b < pairs.Count; b++)
        {
            int i = Math.Min(pairs[a].Position1, pairs[a].Position2);
            int j = Math.Max(pairs[a].Position1, pairs[a].Position2);
            int k = Math.Min(pairs[b].Position1, pairs[b].Position2);
            int l = Math.Max(pairs[b].Position1, pairs[b].Position2);
            if (k < i) (i, j, k, l) = (k, l, i, j);
            if (i < k && k < j && j < l) count++;
        }
        return count;
    }

    /// <summary>
    /// Asserts the FULL documented well-formedness of a detection result for the given input:
    ///   • count equals the independent oracle (no missed / spurious crossing),
    ///   • every reported knot satisfies strict Start1&lt;Start2&lt;End1&lt;End2 (INV-01),
    ///   • Start1&lt;End1 and Start2&lt;End2 (each normalized pair is a proper open&lt;close span),
    ///   • CrossingPairs carries exactly the two original pairs,
    ///   • detection is order-independent (INV-05): shuffling the input does not change the count.
    /// </summary>
    private static IReadOnlyList<Pseudoknot> AssertWellFormed(IReadOnlyList<BasePair> pairs, Random? shuffleRng = null)
    {
        List<Pseudoknot> knots = null!;
        Action act = () => knots = DetectPseudoknots(pairs).ToList();
        act.Should().NotThrow("detection is a pure positional scan and must never throw");

        knots.Count.Should().Be(ExpectedCrossings(pairs), "detected count must equal the oracle");

        foreach (var pk in knots)
        {
            // INV-01: strict crossing i < k < j < l with i=Start1, j=End1, k=Start2, l=End2.
            pk.Start1.Should().BeLessThan(pk.Start2, "Start1 < Start2 (opening order)");
            pk.Start2.Should().BeLessThan(pk.End1, "Start2 < End1 (crossing, not nested)");
            pk.End1.Should().BeLessThan(pk.End2, "End1 < End2 (crossing, not nested)");
            // Each normalized pair is a proper span open < close.
            pk.Start1.Should().BeLessThan(pk.End1, "Start1 < End1 (first pair open < close)");
            pk.Start2.Should().BeLessThan(pk.End2, "Start2 < End2 (second pair open < close)");
            pk.CrossingPairs.Should().HaveCount(2, "a pseudoknot is reported per crossing pair-of-pairs");
        }

        // INV-05: order independence — shuffling input preserves the crossing count.
        if (shuffleRng is not null && pairs.Count > 2)
        {
            var shuffled = pairs.OrderBy(_ => shuffleRng.Next()).ToList();
            DetectPseudoknots(shuffled).Count().Should().Be(knots.Count,
                "detection must be order-independent over the same pair set");
        }

        return knots;
    }

    #endregion

    #region RNA-PSEUDOKNOT-001 — DetectPseudoknots

    // ───────────────────────────────────────────────────────────────────────
    // POSITIVE sanity — the canonical H-type pseudoknot "([)]" and the worked
    // example from the doc (Pseudoknot_Detection.md §7.1).
    // ───────────────────────────────────────────────────────────────────────
    #region Positive sanity — classic crossing

    [Test]
    public void Detect_HTypePseudoknot_BracketSquare_IsDetectedWithDocumentedEndpoints()
    {
        // "([)]" — pairs (0,2) and (1,3): 0 < 1 < 2 < 3 → exactly one crossing (§7.1 worked example).
        var pairs = new[] { Pair(0, 2), Pair(1, 3) };

        var knots = AssertWellFormed(pairs);

        knots.Should().HaveCount(1);
        var pk = knots[0];
        pk.Start1.Should().Be(0);
        pk.End1.Should().Be(2);
        pk.Start2.Should().Be(1);
        pk.End2.Should().Be(3);
        pk.CrossingPairs.Should().BeEquivalentTo(pairs, "both original pairs are carried");
    }

    [Test]
    public void Detect_InterleavedStems_ClassicPseudoknot_IsDetected()
    {
        // Two interleaved stems: (0,4) and (2,6) → 0 < 2 < 4 < 6 → crossing.
        var pairs = new[] { Pair(0, 4), Pair(2, 6) };

        var knots = AssertWellFormed(pairs);

        knots.Should().HaveCount(1);
        knots[0].Start1.Should().Be(0);
        knots[0].Start2.Should().Be(2);
        knots[0].End1.Should().Be(4);
        knots[0].End2.Should().Be(6);
    }

    [Test]
    public void Detect_CrossingIsIndependentOfStorageDirection_AndInputOrder()
    {
        // Same crossing (0,4)+(2,6) but stored as (close,open) and in reverse order: identical result.
        var forward = new[] { Pair(0, 4), Pair(2, 6) };
        var reversedStorage = new[] { Pair(6, 2), Pair(4, 0) };

        DetectPseudoknots(reversedStorage).Should().HaveCount(1,
            "(close,open) storage is normalized to the same crossing");
        DetectPseudoknots(forward).Count().Should().Be(DetectPseudoknots(reversedStorage).Count());
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────
    // BE = Boundary Exploitation — NESTED-ONLY (the canonical negative).
    // A fully nested structure (((...))) has NO crossing pairs. The boundary is
    // the off-by-one: pairs that share/touch an endpoint must NOT cross because
    // the i<k<j<l inequalities are strict (Pseudoknot_Detection.md §2.4, §6.1).
    // ───────────────────────────────────────────────────────────────────────
    #region BE — nested-only / touching boundaries (no pseudoknot)

    [Test]
    public void Detect_FullyNestedStructure_ReportsNoPseudoknot()
    {
        // "(((...)))" — pairs (0,8),(1,7),(2,6): every pair fully inside the previous → all nested.
        var pairs = new[] { Pair(0, 8), Pair(1, 7), Pair(2, 6) };

        var knots = AssertWellFormed(pairs, Rng(1));

        knots.Should().BeEmpty("a fully nested structure has no crossing pairs (INV-02)");
    }

    [Test]
    public void Detect_DisjointPairs_ReportNoPseudoknot()
    {
        // (0,2) then (3,5): j=2 < k=3 → disjoint → no crossing (INV-03).
        var pairs = new[] { Pair(0, 2), Pair(3, 5) };

        AssertWellFormed(pairs).Should().BeEmpty("disjoint ranges do not cross (INV-03)");
    }

    [TestCase(0, 4, 4, 6, TestName = "Detect_ShareInnerOuterEndpoint_j_equals_k_NotCrossing")] // i<k, j=k → fails k<j
    [TestCase(0, 4, 0, 6, TestName = "Detect_ShareOpeningEndpoint_i_equals_k_NotCrossing")]     // i=k → fails i<k
    [TestCase(0, 4, 2, 4, TestName = "Detect_ShareClosingEndpoint_j_equals_l_NotCrossing")]     // j=l → fails j<l (nested-touch)
    [TestCase(0, 6, 2, 4, TestName = "Detect_StrictlyNestedInside_NotCrossing")]                // i<k<l<j → nested
    [TestCase(0, 6, 0, 6, TestName = "Detect_IdenticalPairs_NotCrossing")]                      // duplicate
    public void Detect_TouchingOrNestedBoundary_IsNeverAPseudoknot(int i, int j, int k, int l)
    {
        // The strict-inequality boundary: any arrangement that merely TOUCHES an endpoint or is
        // nested must NOT be flagged. A ≤ off-by-one here would wrongly report a pseudoknot.
        var pairs = new[] { Pair(i, j), Pair(k, l) };

        AssertWellFormed(pairs).Should().BeEmpty(
            $"({i},{j})+({k},{l}) does not satisfy STRICT i<k<j<l");
    }

    [Test]
    public void Detect_DegeneratePair_ZeroWidth_NeverCrosses()
    {
        // ASM-01: a pair with i = j cannot cross anything (it has no interior).
        var pairs = new[] { Pair(3, 3), Pair(0, 6), Pair(2, 4) };

        AssertWellFormed(pairs, Rng(2)).Should().BeEmpty(
            "a zero-width pair (i=j) is never part of a pseudoknot (ASM-01)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────
    // BE = Boundary Exploitation — EMPTY / degenerate cardinality.
    // null, empty, single pair → empty result, NO exception, NO DivideByZero
    // (Pseudoknot_Detection.md INV-04, §3.3, §6.1).
    // ───────────────────────────────────────────────────────────────────────
    #region BE — empty / single / null

    [Test]
    [CancelAfter(5000)]
    public void Detect_NullInput_ReturnsEmptyWithoutException()
    {
        Action act = () => DetectPseudoknots(null!).ToList();
        act.Should().NotThrow("null is documented to yield an empty result (INV-04)");
        DetectPseudoknots(null!).Should().BeEmpty();
    }

    [Test]
    [CancelAfter(5000)]
    public void Detect_EmptySet_ReturnsEmptyWithoutException()
    {
        var pairs = Array.Empty<BasePair>();
        AssertWellFormed(pairs).Should().BeEmpty("< 2 pairs cannot cross (INV-04)");
    }

    [Test]
    [CancelAfter(5000)]
    public void Detect_SinglePair_ReturnsEmpty()
    {
        var pairs = new[] { Pair(0, 9) };
        AssertWellFormed(pairs).Should().BeEmpty("a single pair cannot cross (INV-04)");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────
    // BE = Boundary Exploitation — integer / index boundaries.
    // Extreme positions (0, large, int.MaxValue) must not overflow or throw; the
    // comparison is pure ordering so MaxInt endpoints behave as ordinary positions.
    // ───────────────────────────────────────────────────────────────────────
    #region BE — integer boundaries

    [Test]
    public void Detect_ExtremePositions_CrossingStillDetected_NoOverflow()
    {
        // Crossing using huge positions: 0 < 100 < MaxInt-1 < MaxInt.
        var pairs = new[]
        {
            Pair(0, int.MaxValue - 1),
            Pair(100, int.MaxValue),
        };

        var knots = AssertWellFormed(pairs);
        knots.Should().HaveCount(1, "ordering is unaffected by the magnitude of the positions");
        knots[0].Start1.Should().Be(0);
        knots[0].End2.Should().Be(int.MaxValue);
    }

    [Test]
    public void Detect_NegativeAndZeroPositions_FollowSameOrderingRule()
    {
        // Negative positions are still just integers; crossing -5 < -2 < 0 < 3 holds.
        var pairs = new[] { Pair(-5, 0), Pair(-2, 3) };
        AssertWellFormed(pairs).Should().HaveCount(1,
            "the crossing test is pure integer ordering, independent of sign");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────
    // Randomized fuzz — random base-pair sets (incl. nested, disjoint, crossing,
    // degenerate, reversed-storage) must ALWAYS match the oracle, be well formed,
    // and never throw. Also a guaranteed-nested generator (the canonical negative
    // stress) and a guaranteed-crossing generator (the canonical positive stress).
    // ───────────────────────────────────────────────────────────────────────
    #region Randomized fuzz

    [Test]
    public void Detect_RandomPairSets_AlwaysMatchOracleAndNeverThrow()
    {
        var rng = Rng(20260620);
        for (int iter = 0; iter < 4_000; iter++)
        {
            int n = rng.Next(0, 9);
            var pairs = new List<BasePair>(n);
            for (int p = 0; p < n; p++)
            {
                int x = rng.Next(-3, 15);
                int y = rng.Next(-3, 15);
                // Randomly store as (open,close) or (close,open) to exercise normalization.
                pairs.Add(rng.Next(2) == 0 ? Pair(x, y) : Pair(y, x));
            }
            AssertWellFormed(pairs, rng);
        }
    }

    [Test]
    public void Detect_RandomFullyNestedStructures_NeverReportPseudoknot()
    {
        // Canonical negative stress: build deeply nested pairs (0,2m-1),(1,2m-2),... → 0 crossings.
        var rng = Rng(424242);
        for (int iter = 0; iter < 1_000; iter++)
        {
            int depth = rng.Next(2, 12);
            var pairs = new List<BasePair>(depth);
            for (int d = 0; d < depth; d++)
                pairs.Add(Pair(d, 2 * depth - 1 - d)); // strictly nested onion

            AssertWellFormed(pairs, rng).Should().BeEmpty(
                "a fully nested onion never crosses (INV-02), regardless of depth");
        }
    }

    [Test]
    public void Detect_RandomGuaranteedCrossing_AlwaysReportsAtLeastOne()
    {
        // Canonical positive stress: pick i<k<j<l → always at least one crossing in the set.
        var rng = Rng(13579);
        for (int iter = 0; iter < 1_000; iter++)
        {
            int i = rng.Next(0, 20);
            int k = i + 1 + rng.Next(0, 10);
            int j = k + 1 + rng.Next(0, 10);
            int l = j + 1 + rng.Next(0, 10);
            var pairs = new[] { Pair(i, j), Pair(k, l) };

            var knots = AssertWellFormed(pairs);
            knots.Should().HaveCount(1,
                $"({i},{j})+({k},{l}) satisfies strict i<k<j<l so must cross");
        }
    }

    #endregion

    #endregion
}
