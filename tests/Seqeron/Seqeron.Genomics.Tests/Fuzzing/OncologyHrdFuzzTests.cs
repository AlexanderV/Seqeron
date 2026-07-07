using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology homologous-recombination-deficiency area — ONCO-HRD-001.
/// The units under test are the deterministic genomic-scar entry points
/// <see cref="OncologyAnalyzer.CalculateHRDScore(int, int, int)"/> (the unweighted sum),
/// <see cref="OncologyAnalyzer.ClassifyHRDStatus(int)"/> (the 42 cutoff classifier), and
/// <see cref="OncologyAnalyzer.DetectHRD(OncologyAnalyzer.HrdComponents)"/> (sum + classify),
/// implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, and no *unhandled* runtime fault (a silently overflowed / NEGATIVE
/// HRD score leaking out). Every input must resolve to EITHER a well-defined,
/// theory-correct value OR a *documented, intentional* outcome (here, an
/// <see cref="ArgumentOutOfRangeException"/> for a negative component, a negative
/// score-to-classify, or a sum that genuinely overruns the int range).
/// For HRD the headline hazards (checklist row 94, targets
/// "no events, negative component, extreme counts") are:
///   • no events — all three components 0 ⇒ HRD score 0, classified
///     HRD-negative, no throw (§6.1 "All components 0 ⇒ score 0, HRD-negative");
///   • negative component — any of loh/tai/lst &lt; 0 ⇒ documented
///     ArgumentOutOfRangeException, NEVER a silently negative HRD that slips
///     past the unweighted sum (§3.3, §6.1, INV-04 "each component ≥ 0");
///   • extreme counts — components near Int32.MaxValue ⇒ the sum MUST NOT wrap
///     to a negative score (INV-04 "score ≥ 0"); the implementation sums in a
///     wider width so an in-range total stays exact and non-negative, and a
///     total that truly exceeds Int32.MaxValue is rejected with
///     ArgumentOutOfRangeException rather than silently overflowing.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-HRD-001 — Homologous recombination deficiency score (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 94.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 94): "no events, negative component, extreme counts".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// HRD_Score.md (docs/algorithms/Oncology/HRD_Score.md):
///   • HRD = LOH + TAI + LST, an UNWEIGHTED sum                          (§2.2, INV-01)
///   • Sum is order-independent / commutative                            (INV-02)
///   • status = HrdHigh iff score ≥ 42 (inclusive), else HrdNegative     (§2.2, INV-03)
///   • score ≥ 0 and each component ≥ 0                                  (INV-04)
///   • All components 0 ⇒ score 0, HRD-negative                         (§6.1)
///   • Score exactly 42 ⇒ HrdHigh; score 41 ⇒ HrdNegative                (§6.1)
///   • Any negative component ⇒ ArgumentOutOfRangeException
///       (CalculateHRDScore, DetectHRD)                                 (§3.3)
///   • Negative score-to-classify ⇒ ArgumentOutOfRangeException
///       (ClassifyHRDStatus)                                            (§3.3)
///   • HrdHighScoreThreshold == 42                                       (§4.2)
///
/// SOURCE FIX (this row): CalculateHRDScore previously summed in int
/// (`return loh + tai + lst;`), so extreme counts near Int32.MaxValue wrapped to
/// a NEGATIVE score — violating INV-04 ("score ≥ 0"). The sum is now computed in
/// long; an in-range total is returned exactly, a total that exceeds
/// Int32.MaxValue throws ArgumentOutOfRangeException. No test was weakened.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyHrdFuzzTests
{
    private const int Cutoff = HrdHighScoreThreshold; // 42, Telli 2016, inclusive

    // ── Well-formed-score assertion helper ───────────────────────────────────
    // Pins the documented numeric contract on EVERY returned score: it must be
    // NON-negative (INV-04). This is what stops a fuzz test from rubber-stamping
    // a silently overflowed (wrapped-negative) HRD score green.
    private static void AssertWellFormedScore(int score) =>
        score.Should().BeGreaterThanOrEqualTo(0, "HRD score = LOH + TAI + LST is a count sum (INV-04)");

    #region ONCO-HRD-001 — Homologous recombination deficiency (positive sanity)

    // ── POSITIVE sanity: hand-computed sum and the correct High/Negative class ─
    [Test]
    public void CalculateHRDScore_KnownComponents_HandComputedSumAndHigh()
    {
        // LOH=20, TAI=15, LST=10 ⇒ 45; 45 ≥ 42 ⇒ HRD-high.
        int score = CalculateHRDScore(loh: 20, tai: 15, lst: 10);

        AssertWellFormedScore(score);
        score.Should().Be(45, "20 + 15 + 10 = 45 (INV-01)");
        ClassifyHRDStatus(score).Should().Be(HrdStatus.HrdHigh, "45 ≥ 42 ⇒ HRD-high (INV-03)");
    }

    [Test]
    public void DetectHRD_DocWorkedExample_Score47IsHigh()
    {
        // Docs §7.1: HrdComponents(Loh: 20, Tai: 15, Lst: 12) ⇒ Score 47, HrdHigh.
        HrdResult result = DetectHRD(new HrdComponents(Loh: 20, Tai: 15, Lst: 12));

        AssertWellFormedScore(result.Score);
        result.Score.Should().Be(47, "20 + 15 + 12 = 47 (§7.1)");
        result.Status.Should().Be(HrdStatus.HrdHigh, "47 ≥ 42 ⇒ HRD-high (§7.1)");
        result.Components.Should().Be(new HrdComponents(20, 15, 12), "components echoed back unchanged");
    }

    [Test]
    public void ClassifyHRDStatus_AtAndBelowCutoff_IsInclusiveAt42()
    {
        // §6.1: score exactly 42 ⇒ HrdHigh; score 41 ⇒ HrdNegative.
        HrdHighScoreThreshold.Should().Be(42, "Telli 2016 cutoff (§4.2)");
        ClassifyHRDStatus(Cutoff).Should().Be(HrdStatus.HrdHigh, "cutoff is inclusive ≥ 42 (INV-03)");
        ClassifyHRDStatus(Cutoff - 1).Should().Be(HrdStatus.HrdNegative, "41 < 42 ⇒ HRD-negative (§6.1)");
    }

    #endregion

    #region ONCO-HRD-001 — BE: no events (all components 0 ⇒ score 0, HRD-negative)

    [Test]
    public void CalculateHRDScore_NoEvents_IsZeroAndNegative_NoThrow()
    {
        // §6.1: all components 0 (near-diploid tumour) ⇒ score 0, HRD-negative.
        // The hazard guarded here is that a fully-zero input must NOT trip any
        // guard — it is a perfectly valid, well-defined HRD score of exactly 0.
        int score = CalculateHRDScore(loh: 0, tai: 0, lst: 0);

        AssertWellFormedScore(score);
        score.Should().Be(0, "0 + 0 + 0 = 0 (§6.1)");
        ClassifyHRDStatus(score).Should().Be(HrdStatus.HrdNegative, "0 < 42 ⇒ HRD-negative (§6.1)");
    }

    [Test]
    public void DetectHRD_NoEvents_IsZeroAndNegative()
    {
        HrdResult result = DetectHRD(new HrdComponents(0, 0, 0));

        AssertWellFormedScore(result.Score);
        result.Score.Should().Be(0);
        result.Status.Should().Be(HrdStatus.HrdNegative, "no events ⇒ HRD-negative (§6.1)");
    }

    #endregion

    #region ONCO-HRD-001 — BE: negative component (documented ArgumentOutOfRangeException)

    [Test]
    public void CalculateHRDScore_NegativeLoh_ThrowsArgumentOutOfRange()
    {
        // §3.3 / INV-04: a negative component is rejected; the unweighted sum must
        // NEVER let a negative count leak into a (possibly negative) HRD score.
        Action act = () => CalculateHRDScore(loh: -1, tai: 10, lst: 10);
        act.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("loh");
    }

    [Test]
    public void CalculateHRDScore_NegativeTai_ThrowsArgumentOutOfRange()
    {
        Action act = () => CalculateHRDScore(loh: 10, tai: -1, lst: 10);
        act.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("tai");
    }

    [Test]
    public void CalculateHRDScore_NegativeLst_ThrowsArgumentOutOfRange()
    {
        Action act = () => CalculateHRDScore(loh: 10, tai: 10, lst: -1);
        act.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("lst");
    }

    [Test]
    public void DetectHRD_AnyNegativeComponent_ThrowsArgumentOutOfRange()
    {
        Action act = () => DetectHRD(new HrdComponents(Loh: 30, Tai: -5, Lst: 30));
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void ClassifyHRDStatus_NegativeScore_ThrowsArgumentOutOfRange()
    {
        // §3.3: a negative score-to-classify is itself out of contract.
        Action act = () => ClassifyHRDStatus(-1);
        act.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("score");
    }

    [Test]
    public void CalculateHRDScore_FuzzedNegativeComponents_AlwaysThrow_NeverNegativeScore()
    {
        // BE: scatter a negative into exactly one of the three slots and confirm
        // the documented throw — a negative count must never silently sum.
        var rng = new Random(94_0001);
        for (int i = 0; i < 600; i++)
        {
            int a = rng.Next(0, 50);
            int b = rng.Next(0, 50);
            int neg = -rng.Next(1, int.MaxValue); // strictly negative
            int slot = rng.Next(3);

            (int loh, int tai, int lst) = slot switch
            {
                0 => (neg, a, b),
                1 => (a, neg, b),
                _ => (a, b, neg),
            };

            Action act = () => CalculateHRDScore(loh, tai, lst);
            act.Should().Throw<ArgumentOutOfRangeException>(
                "a negative HRD component is documented out of contract (§3.3)");
        }
    }

    #endregion

    #region ONCO-HRD-001 — BE: extreme counts (no integer overflow / negative wrap)

    [Test]
    public void CalculateHRDScore_ExtremeCountsSummingToIntMax_StaysExactAndNonNegative()
    {
        // The total exactly equals Int32.MaxValue: it is representable, so it must
        // be returned exactly and stay NON-negative (INV-04) — no wrap.
        int score = CalculateHRDScore(loh: int.MaxValue - 30, tai: 20, lst: 10);

        AssertWellFormedScore(score);
        score.Should().Be(int.MaxValue, "the in-range total must be exact (INV-01), not wrapped");
        ClassifyHRDStatus(score).Should().Be(HrdStatus.HrdHigh, "≫ 42 ⇒ HRD-high (INV-03)");
    }

    [Test]
    public void CalculateHRDScore_ExtremeCountsOverflowingIntMax_ThrowsNotWrapNegative()
    {
        // BE headline: components near Int32.MaxValue whose sum overruns the int
        // range MUST NOT wrap to a negative score (the pre-fix int sum did exactly
        // that). The documented out-of-contract outcome is ArgumentOutOfRangeException.
        Action act = () => CalculateHRDScore(loh: int.MaxValue, tai: 10, lst: 5);
        act.Should().Throw<ArgumentOutOfRangeException>(
            "an HRD total exceeding Int32.MaxValue is rejected, never silently wrapped negative (INV-04)");
    }

    [Test]
    public void CalculateHRDScore_MaxComponentAlone_IsExactAndNonNegative()
    {
        // A single component at Int32.MaxValue with the others 0 is in range.
        int score = CalculateHRDScore(loh: int.MaxValue, tai: 0, lst: 0);

        AssertWellFormedScore(score);
        score.Should().Be(int.MaxValue, "MaxValue + 0 + 0 is exact (INV-01)");
    }

    [Test]
    public void CalculateHRDScore_FuzzedExtremeCounts_NeverWrapNegative()
    {
        // BE: hammer large component triples; whatever the verdict (exact score or
        // documented throw), a NEGATIVE score must NEVER leak out.
        var rng = new Random(94_0002);
        for (int i = 0; i < 1000; i++)
        {
            int loh = rng.Next(0, int.MaxValue);
            int tai = rng.Next(0, int.MaxValue);
            int lst = rng.Next(0, int.MaxValue);
            long expected = (long)loh + tai + lst;

            try
            {
                int score = CalculateHRDScore(loh, tai, lst);
                AssertWellFormedScore(score);
                ((long)score).Should().Be(expected,
                    "an in-range total is the exact long sum, never a wrapped value (INV-01/INV-04)");
                expected.Should().BeLessThanOrEqualTo(int.MaxValue,
                    "a value was returned only because the total fit in int");
            }
            catch (ArgumentOutOfRangeException)
            {
                expected.Should().BeGreaterThan(int.MaxValue,
                    "the throw is reserved for totals that genuinely overrun the int range");
            }
        }
    }

    [Test]
    public void DetectHRD_ExtremeInRangeComponents_HighAndNonNegative()
    {
        HrdResult result = DetectHRD(new HrdComponents(int.MaxValue - 100, 60, 40));

        AssertWellFormedScore(result.Score);
        result.Score.Should().Be(int.MaxValue);
        result.Status.Should().Be(HrdStatus.HrdHigh);
    }

    #endregion

    #region ONCO-HRD-001 — INV-02/INV-03: commutativity & cutoff across fuzz

    [Test]
    public void CalculateHRDScore_IsOrderIndependent_AcrossFuzzedTriples()
    {
        // INV-02: the unweighted sum is commutative — any permutation of the three
        // components yields the identical score.
        var rng = new Random(94_0003);
        for (int i = 0; i < 1000; i++)
        {
            int a = rng.Next(0, 1_000_000);
            int b = rng.Next(0, 1_000_000);
            int c = rng.Next(0, 1_000_000);

            int s1 = CalculateHRDScore(a, b, c);
            AssertWellFormedScore(s1);
            CalculateHRDScore(c, a, b).Should().Be(s1, "sum is order-independent (INV-02)");
            CalculateHRDScore(b, c, a).Should().Be(s1, "sum is order-independent (INV-02)");
            CalculateHRDScore(c, b, a).Should().Be(s1, "sum is order-independent (INV-02)");
        }
    }

    [Test]
    public void ClassifyHRDStatus_AgreesWithCutoff_AcrossFuzzedScores()
    {
        // INV-03: status = HrdHigh iff score ≥ 42 (inclusive). Sweep boundary-dense
        // and broad scores and confirm the classifier never disagrees with the cutoff.
        var rng = new Random(94_0004);
        for (int i = 0; i < 2000; i++)
        {
            // Half boundary-dense around 42, half broad.
            int score = i % 2 == 0
                ? Math.Max(0, Cutoff + rng.Next(-5, 6))
                : rng.Next(0, 100_000);

            HrdStatus status = ClassifyHRDStatus(score);
            HrdStatus expected = score >= Cutoff ? HrdStatus.HrdHigh : HrdStatus.HrdNegative;
            status.Should().Be(expected, $"score {score} vs inclusive cutoff 42 (INV-03)");
        }
    }

    #endregion
}
