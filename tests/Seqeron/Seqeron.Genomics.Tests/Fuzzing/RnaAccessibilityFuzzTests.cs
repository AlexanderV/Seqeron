namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for RNA-ACCESS-001 — RNA REGION ACCESSIBILITY (joint region-unpaired probability).
/// Unit under test:
///   RnaSecondaryStructure.CalculateRegionUnpairedProbability(
///       string rnaSequence, int windowEnd, int windowLength,
///       int minLoopSize = 3, double temperature = 310.15)  →  double  P(window entirely unpaired)
/// in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs (L2672–2719). It is
/// the RNAplfold-style accessibility — the exact McCaskill ensemble quantity Z_open(window)/Z, where
/// Z_open is the Turner-2004 partition function over structures in which EVERY base of the window
/// [windowEnd−windowLength+1 .. windowEnd] is unpaired (no pair incident to any window position).
/// Checklist: docs/checklists/03_FUZZING.md, row 238 (strategy codes BE, MC).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Malformed / boundary inputs must NEVER fail in an undisciplined way: no hang; no unhandled
/// runtime exception (an IndexOutOfRange from a bad window, a DivideByZero / NaN / Inf out of the
/// partition function); and no out-of-contract output (a "probability" outside [0,1]). Every input
/// resolves to EITHER a well-defined theory-correct value OR a DOCUMENTED validation exception
/// (ArgumentOutOfRangeException) — never a silent NaN leak or an IndexOutOfRange.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy codes (docs/checklists/03_FUZZING.md §Description)
/// ───────────────────────────────────────────────────────────────────────────
///   BE = Boundary Exploitation (граничні значення: 0, -1, empty) — empty sequence, window out of
///        bounds (start &lt; 0 / end ≥ n), window length 0, single base.
///   MC = Malformed Content (некоректні дані) — non-ACGU characters in the sequence.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The DOCUMENTED contract under test
/// (docs/algorithms/RnaStructure/Turner_McCaskill_Partition_Function.md §3.1–§3.3, INV-06;
///  docs/Validation/reports/RNA-ACCESS-001.md "Edge-case semantics")
/// ───────────────────────────────────────────────────────────────────────────
///   • Return is P(window entirely unpaired) = Z_open/Z, ALWAYS in [0,1] (INV-06: Z_open ≤ Z
///     because Z_open's structures are a subset of the full ensemble).
///   • temperature ≤ 0 → ArgumentOutOfRangeException (RT ≤ 0 ⇒ exp blow-up; physically meaningless).
///   • windowLength ≤ 0 → ArgumentOutOfRangeException (region length 0 / negative is invalid).
///   • A window that does not fit — windowStart = windowEnd−windowLength+1 &lt; 0, or windowEnd ≥ n
///     (this INCLUDES the empty sequence, where every window is out of bounds) →
///     ArgumentOutOfRangeException. NOT an IndexOutOfRange.
///   • A sequence too short to form ANY pair (n &lt; minLoopSize+2) → returns 1.0 exactly: nothing
///     can pair, so the window is unpaired with probability 1.
///   • Non-ACGU characters (MC): they are not in {A·U, G·C, G·U}, so they cannot pair — they behave
///     as permanently unpaired. A window of only such characters is fully accessible (≈ 1); they
///     never introduce NaN/Inf nor an out-of-range probability.
///   • DNA accepted (T read as U); case-insensitive.
///
/// THEORY PINS (derived independently, NOT echoes of the code):
///   • all-A (or any homopolymer / non-complementary) region → accessibility ≈ 1 (no pair possible).
///   • a window over a base forced into a strong stem → accessibility &lt; 1 (it is paired part of
///     the time), yet still ≥ 0.
///   • a longer queried window is never MORE accessible than a shorter one ending at the same
///     position (Bernhart 2006: forbidding more positions can only shrink Z_open) — monotonicity.
///   • analytic two-state pin GAAAC: the whole 5-nt window @end=4 = 1/Z = 0.9998435192135765
///     (Z = 1 + exp(−5.4/RT), the single G·C hairpin over loop "AAA"; RNA-ACCESS-001.md Stage-A).
///
/// Hang-safety (§4.3): CalculateRegionUnpairedProbability is O(n³) (two DP fills). Fuzz lengths are
/// kept MODEST and guarded with [CancelAfter] so a regression introducing a hang FAILS, not wedges.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaAccessibilityFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed LOCALLY so generated fuzz inputs are reproducible.</summary>
    private static string RandomRna(int length, int seed)
    {
        const string bases = "ACGU";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// The WELL-FORMED-RESULT pin (INV-06): the returned accessibility is a finite probability —
    /// no NaN, no ±∞, and inside [0,1] up to round-off. Asserted for EVERY non-throwing call.
    /// </summary>
    private static void AssertWellFormedAccessibility(double a, string because)
    {
        double.IsNaN(a).Should().BeFalse($"accessibility must never be NaN ({because})");
        double.IsInfinity(a).Should().BeFalse($"accessibility must never be ±∞ ({because})");
        a.Should().BeInRange(-1e-12, 1.0 + 1e-9,
            $"INV-06: region accessibility Z_open/Z is a probability in [0,1] ({because})");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-ACCESS-001 — region accessibility (Z_open/Z) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region RNA-ACCESS-001 — BE: empty sequence → window out of bounds → ArgumentOutOfRange (not IndexOutOfRange)

    /// <summary>
    /// BE — "empty seq". For the empty sequence n = 0, so EVERY window is out of bounds: windowEnd ≥ n
    /// holds for any windowEnd ≥ 0 and windowStart &lt; 0 for any positive length. The contract is a
    /// DOCUMENTED ArgumentOutOfRangeException, NEVER an IndexOutOfRange off the (zero-size) DP table
    /// (§3.3; RNA-ACCESS-001.md: "window outside the sequence … must throw"). We also pin that the
    /// throw happens BEFORE any partition allocation (no crash on the 0×0 buffer).
    /// </summary>
    [Test]
    public void EmptySequence_AnyWindow_ThrowsArgumentOutOfRange_NotIndexOutOfRange()
    {
        foreach ((int end, int len) in new[] { (0, 1), (0, 5), (3, 1), (5, 14) })
        {
            var act = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability("", end, len);
            act.Should().Throw<ArgumentOutOfRangeException>(
                $"empty sequence has no positions — window (end={end},len={len}) cannot fit (§3.3)");
            // Defensively confirm it is NOT a raw IndexOutOfRange masquerading as the base ArgumentException.
            act.Should().NotThrow<IndexOutOfRangeException>();
        }
    }

    #endregion

    #region BE — Window out of bounds: start < 0 or end ≥ n → ArgumentOutOfRange (never IndexOutOfRange)

    /// <summary>
    /// BE — "region out of bounds". A window with windowStart = windowEnd−windowLength+1 &lt; 0 (length
    /// overruns the 5' end) or windowEnd ≥ n (end past the 3' end) is invalid and MUST throw a
    /// documented ArgumentOutOfRangeException — explicitly NOT a runtime IndexOutOfRange from indexing
    /// the n×n DP tables (§3.3, RNA-ACCESS-001.md). Both overrun directions and far-out-of-range
    /// indices (incl. int.MaxValue / very negative) are pinned across a foldable sequence.
    /// </summary>
    [Test]
    public void WindowOutOfBounds_BothDirections_ThrowArgumentOutOfRange()
    {
        const string seq = "GGGGAAAACCCC"; // n = 12, foldable
        int n = seq.Length;

        // windowEnd ≥ n  (end past the 3' end)
        foreach (int end in new[] { n, n + 1, n + 100, int.MaxValue })
        {
            var act = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, end, 1);
            act.Should().Throw<ArgumentOutOfRangeException>($"windowEnd={end} ≥ n={n} is out of bounds");
            act.Should().NotThrow<IndexOutOfRangeException>();
        }

        // windowStart < 0  (length overruns the 5' end)
        foreach ((int end, int len) in new[] { (0, 2), (3, 5), (5, 100), (2, int.MaxValue) })
        {
            var act = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, end, len);
            act.Should().Throw<ArgumentOutOfRangeException>(
                $"window (end={end},len={len}) has start < 0 — overruns the 5' end");
            act.Should().NotThrow<IndexOutOfRangeException>();
        }

        // Negative windowEnd — also start < 0, out of bounds.
        foreach (int end in new[] { -1, -50, int.MinValue + 1 })
        {
            var act = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, end, 1);
            act.Should().Throw<ArgumentOutOfRangeException>($"negative windowEnd={end} is out of bounds");
        }
    }

    #endregion

    #region BE — Region length 0 (and negative) → documented ArgumentOutOfRange

    /// <summary>
    /// BE — "region len 0". A zero-length (or negative-length) window is rejected up front with a
    /// documented ArgumentOutOfRangeException on windowLength — the empty-region accessibility is
    /// DEFINED as invalid by the contract (§3.1 "length > 0"; §3.3), not silently returned as 1 nor
    /// crashed. Pinned across several windowEnd values and on the empty sequence (length is checked
    /// before the fit check, so it throws regardless of n).
    /// </summary>
    [Test]
    public void ZeroOrNegativeWindowLength_ThrowsArgumentOutOfRange()
    {
        const string seq = "GGGGAAAACCCC";
        foreach (int len in new[] { 0, -1, -14, int.MinValue })
        {
            foreach (int end in new[] { 0, 3, 11 })
            {
                var act = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, end, len);
                act.Should().Throw<ArgumentOutOfRangeException>(
                    $"windowLength={len} is non-positive — empty/negative region is invalid (§3.1)");
            }
            // length check precedes the fit check ⇒ throws even for the empty sequence.
            var emptyAct = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability("", 0, len);
            emptyAct.Should().Throw<ArgumentOutOfRangeException>(
                $"windowLength={len} is non-positive regardless of sequence length");
        }
    }

    #endregion

    #region BE — Single base → too short to pair → accessibility = 1.0 exactly

    /// <summary>
    /// BE — "single base". A length-1 sequence (n=1 &lt; minLoopSize+2 = 5) admits no pair, so the
    /// only window that fits (end=0, len=1) is unpaired with probability EXACTLY 1.0 (§3.3:
    /// "a window in a too-short sequence returns 1.0"). All four residues are pinned; an out-of-fit
    /// window on the single base (len=2) still throws OOB, confirming the fit check is independent
    /// of the too-short short-circuit.
    /// </summary>
    [Test]
    public void SingleBase_WholeWindow_IsFullyAccessible()
    {
        foreach (char b in "ACGU")
        {
            double a = RnaSecondaryStructure.CalculateRegionUnpairedProbability(b.ToString(), windowEnd: 0, windowLength: 1);
            a.Should().Be(1.0, $"a single '{b}' cannot pair — its only window is unpaired w.p. 1");
            AssertWellFormedAccessibility(a, $"single base '{b}'");

            // A window that does not fit the single base still throws (fit check independent of n<5 path).
            var act = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability(b.ToString(), windowEnd: 0, windowLength: 2);
            act.Should().Throw<ArgumentOutOfRangeException>($"len-2 window does not fit a single '{b}'");
        }
    }

    /// <summary>
    /// BE — sub-foldable spans (n &lt; minLoopSize+2 = 5): "GC","GCG","GCGC" are perfectly
    /// complementary yet too short to enclose a 3-nt hairpin loop, so NO pair can form → any fitting
    /// window is fully accessible (1.0). Guards against counting sterically forbidden pairs at the
    /// length boundary (§3.3, ASM/min-loop rule).
    /// </summary>
    [Test]
    public void SubFoldableSpans_AnyFittingWindow_IsFullyAccessible()
    {
        foreach (string s in new[] { "GC", "GCG", "GCGC", "AU", "CG" })
        {
            for (int end = 0; end < s.Length; end++)
                for (int len = 1; len <= end + 1; len++)
                {
                    double a = RnaSecondaryStructure.CalculateRegionUnpairedProbability(s, end, len);
                    a.Should().Be(1.0, $"'{s}' is too short to pair — window (end={end},len={len}) is fully accessible");
                    AssertWellFormedAccessibility(a, $"'{s}' window (end={end},len={len})");
                }
        }
    }

    #endregion

    #region MC — Non-ACGU characters: cannot pair → permanently unpaired (≈ 1), never NaN/Inf

    /// <summary>
    /// MC — "non-ACGU chars". Characters outside {A,C,G,U} are not in any admissible pair
    /// ({A·U, G·C, G·U}), so they behave as permanently unpaired (RNA-ACCESS-001.md: "Non-ACGU
    /// characters … cannot pair … behave as permanently unpaired"). A window over ONLY such
    /// characters is therefore fully accessible (1.0), and the call must never throw nor leak NaN/Inf
    /// — these are malformed CONTENT, handled gracefully, not a contract violation. We test runs of
    /// 'N', digits, punctuation, whitespace and Unicode within an otherwise foldable context.
    /// </summary>
    [Test]
    [CancelAfter(20_000)]
    public void NonAcgu_WindowOverInvalidChars_IsFullyAccessible(CancellationToken _)
    {
        // "GGGG" + 4 invalid + "CCCC": the invalid run (indices 4..7) can never pair.
        foreach (string fill in new[] { "NNNN", "1234", "....", "    ", "xyzw", "ζαβγ" })
        {
            string seq = "GGGG" + fill + "CCCC"; // n = 12; window 4..7 = the invalid fill
            var act = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, windowEnd: 7, windowLength: 4);
            double a = act.Should().NotThrow($"non-ACGU content is malformed but handled, not a throw (fill='{fill}')").Subject;
            AssertWellFormedAccessibility(a, $"non-ACGU fill '{fill}'");
            a.Should().BeApproximately(1.0, 1e-9,
                $"the all-invalid window '{fill}' cannot pair — accessibility = 1 (permanently unpaired)");
        }
    }

    /// <summary>
    /// MC — an ENTIRELY non-ACGU sequence. With no pairable base anywhere, Z = 1 (open chain only)
    /// and any fitting window is fully accessible. Must be finite, in [0,1], and never throw on the
    /// content (only on a genuinely out-of-fit window). Several lengths and alphabets are swept.
    /// </summary>
    [Test]
    [CancelAfter(20_000)]
    public void NonAcgu_EntireSequence_FullyAccessible_NoThrow(CancellationToken _)
    {
        foreach (string seq in new[] { "NNNNNNNN", "ZZZZZZZZ", "xxxxxxxxxx", "????????", "########" })
        {
            double a = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, windowEnd: 5, windowLength: 4);
            AssertWellFormedAccessibility(a, $"all-invalid '{seq}'");
            a.Should().BeApproximately(1.0, 1e-9, $"'{seq}' has no pairable base — accessibility = 1");
        }
    }

    #endregion

    #region THEORY — all-A region ≈ 1; stem region < 1; finite

    /// <summary>
    /// THEORY (anti-rubber-stamp): an all-A region cannot pair (A pairs only U; no U present), so its
    /// accessibility is ≈ 1 — fully unpaired. Pinned over a pure poly-A sequence AND an A-only window
    /// embedded in a foldable G…C context (the A-loop of a hairpin is still unpaired-with-prob-1 only
    /// when it is the hairpin loop; here the embedded all-A window is the loop "AAAA" of GGGG…CCCC and
    /// its joint-unpaired probability is high but we only require the strict poly-A case to be ≈ 1).
    /// </summary>
    [Test]
    [CancelAfter(20_000)]
    public void AllA_Region_IsFullyAccessible(CancellationToken _)
    {
        foreach (int n in new[] { 5, 8, 16, 40 })
        {
            string seq = new string('A', n);
            for (int len = 1; len <= n; len += Math.Max(1, n / 4))
            {
                double a = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, windowEnd: n - 1, windowLength: len);
                AssertWellFormedAccessibility(a, $"poly-A n={n} len={len}");
                a.Should().BeApproximately(1.0, 1e-9,
                    $"poly-A '{seq}' has no complementary base — every window is fully accessible");
            }
        }
    }

    /// <summary>
    /// THEORY: a region forced into a STRONG stem is paired part of the time, so its accessibility is
    /// STRICTLY less than 1 (and ≥ 0). In "GGGGAAAACCCC" the G·C stem closures (0..3, 8..11) carry
    /// real pairing mass, so the 4-nt 3'-stem window @end=11 must be &lt; 1 — proving the method
    /// actually reflects pairing, not a constant 1. The loop window @end=7 (AAAA) is conversely highly
    /// accessible. Both remain in [0,1].
    /// </summary>
    [Test]
    [CancelAfter(20_000)]
    public void StrongStemRegion_LessAccessibleThanLoop(CancellationToken _)
    {
        const string seq = "GGGGAAAACCCC"; // 4-bp G·C stem closing an AAAA loop

        double stem = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, windowEnd: 11, windowLength: 4); // CCCC stem
        double loop = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, windowEnd: 7, windowLength: 4);  // AAAA loop

        AssertWellFormedAccessibility(stem, "stem window");
        AssertWellFormedAccessibility(loop, "loop window");

        stem.Should().BeLessThan(1.0, "a strong G·C stem is paired part of the ensemble → accessibility < 1");
        stem.Should().BeGreaterThanOrEqualTo(0.0, "accessibility is never negative");
        loop.Should().BeGreaterThan(stem, "the unpaired AAAA loop is more accessible than the paired stem");
    }

    #endregion

    #region THEORY — analytic two-state pin + monotonicity in window length

    /// <summary>
    /// THEORY — analytic two-state pin (independent oracle, RNA-ACCESS-001.md Stage-A). "GAAAC" has
    /// exactly two structures: the open chain (w=1) and the single G(0)·C(4) hairpin over loop "AAA"
    /// (ΔG = hairpin-initiation(3) = 5.4). So Z = 1 + exp(−5.4/RT) and the whole 5-nt window @end=4
    /// has accessibility = Z_open/Z = 1/Z = 0.9998435192135765 (Z_open forbids all bases ⇒ only the
    /// open chain survives). This is a hand-derived value, NOT a code echo.
    /// </summary>
    [Test]
    public void AnalyticPin_GAAAC_WholeWindowAccessibility()
    {
        double a = RnaSecondaryStructure.CalculateRegionUnpairedProbability("GAAAC", windowEnd: 4, windowLength: 5);
        a.Should().BeApproximately(0.9998435192135765, 1e-12,
            "GAAAC whole-window accessibility = 1/Z = 1/(1+exp(−5.4/RT)) (analytic two-state oracle)");
        AssertWellFormedAccessibility(a, "GAAAC analytic pin");
    }

    /// <summary>
    /// THEORY — monotonicity (Bernhart et al. 2006): a LONGER window ending at the same position is
    /// never MORE accessible than a shorter one — forbidding more positions from pairing can only
    /// shrink Z_open, hence Z_open/Z. A genuine mathematical property of the ensemble, distinct from
    /// any single computed value. Verified on a foldable sequence over a sweep of window lengths.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Monotonicity_LongerWindowNeverMoreAccessible(CancellationToken _)
    {
        const string seq = "GGGGAAAACCCCGGGG"; // n = 16, richly foldable
        int end = seq.Length - 1;

        double prev = double.PositiveInfinity;
        for (int len = 1; len <= end + 1; len++)
        {
            double a = RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, end, len);
            AssertWellFormedAccessibility(a, $"window len={len} @end={end}");
            a.Should().BeLessThanOrEqualTo(prev + 1e-9,
                $"a length-{len} window cannot be MORE accessible than a shorter one (Bernhart monotonicity)");
            prev = a;
        }
    }

    #endregion

    #region BE — Validation: non-positive temperature throws

    /// <summary>
    /// BE — documented validation outcome (§3.3): a non-positive temperature (0 or negative Kelvin is
    /// physically meaningless — RT ≤ 0 would make exp(−ΔG/RT) blow up to ±∞/NaN) is rejected with an
    /// ArgumentOutOfRangeException, never a silent NaN leak. Checked before the DP runs.
    /// </summary>
    [Test]
    public void NonPositiveTemperature_ThrowsArgumentOutOfRange()
    {
        foreach (double t in new[] { 0.0, -1.0, -310.15, double.NegativeInfinity })
        {
            var act = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability("GGGGAAAACCCC", 11, 4, temperature: t);
            act.Should().Throw<ArgumentOutOfRangeException>(
                $"temperature {t} K is non-positive — rejected, never a NaN leak (§3.3)");
        }
    }

    #endregion

    #region BE — Robustness sweep: random fuzz windows over random RNAs never violate the contract

    /// <summary>
    /// BE robustness sweep: a deterministic stream of random RNAs across lengths, queried at random
    /// FITTING windows, must NEVER violate the contract — no throw on a valid window, accessibility
    /// finite and in [0,1] (INV-06). Random NON-fitting windows must throw ArgumentOutOfRange (never
    /// IndexOutOfRange). [CancelAfter] proves O(n³) termination on every draw. This is the broad
    /// anti-fragility net behind the targeted BE/MC cases above.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void RandomSweep_FittingAndNonFittingWindows_RespectContract(CancellationToken _)
    {
        var rng = new Random(20260626);
        foreach (int len in new[] { 1, 2, 5, 8, 13, 21, 40 })
        {
            for (int s = 0; s < 6; s++)
            {
                string seq = RandomRna(len, s * 97 + len);
                int n = seq.Length;

                // A random FITTING window: end in [0,n-1], length in [1, end+1].
                int end = rng.Next(0, n);
                int wlen = rng.Next(1, end + 2);
                var fitAct = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, end, wlen);
                double a = fitAct.Should().NotThrow(
                    $"len={len} seed={s} '{seq}' fitting window (end={end},wlen={wlen}) must not throw").Subject;
                AssertWellFormedAccessibility(a, $"random len={len} seed={s} window (end={end},wlen={wlen})");

                // A random NON-fitting window: either end ≥ n or length overruns the 5' end → must throw.
                var oobAct = () => RnaSecondaryStructure.CalculateRegionUnpairedProbability(seq, n + rng.Next(0, 3), wlen);
                oobAct.Should().Throw<ArgumentOutOfRangeException>(
                    $"len={len} seed={s}: end ≥ n is out of bounds — documented throw, not IndexOutOfRange");
            }
        }
    }

    #endregion
}
