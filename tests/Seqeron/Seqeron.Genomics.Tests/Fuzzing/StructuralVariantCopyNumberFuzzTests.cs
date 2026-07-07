using static Seqeron.Genomics.Annotation.StructuralVariantAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for read-depth copy-number variation detection — SV-CNV-001. The unit
/// under test is the windowed read-depth → log2-ratio → integer-copy-number transform in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs:
///   • <see cref="StructuralVariantAnalyzer.DetectCNV"/> — partitions a per-position
///     depth track into non-overlapping fixed-size full windows, takes each window's mean
///     depth, converts it to <c>log2(RD_w / RD_ref)</c> and an integer copy number
///     <c>round(2·2^log2)</c>, emitting one <see cref="StructuralVariantAnalyzer.CopyNumberSegment"/>
///     per non-zero window (canonical);
///   • <see cref="StructuralVariantAnalyzer.SegmentCopyNumber"/> — consumes per-window log2
///     ratios directly and merges maximal runs of equal copy number (segmentation variant).
///
/// SCOPE. This file is scoped strictly to SV-CNV-001 — the deterministic depth→CN numeric
/// transform. It is NOT a statistical segmentation/change-point test (no CBS), and the
/// caller models no tumour purity, GC bias, or allele-specific signal (docs §1, §5.3
/// "Not implemented"). The sibling units split-read breakpoint localization
/// (SV-BREAKPOINT-001) and paired-end SV typing (SV-DETECT-001) are SEPARATE units and
/// are not exercised here.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary inputs to a unit and asserts that the code NEVER
/// fails in an undisciplined way: no hang, no nonsense output, no *unhandled* runtime
/// exception, and no NaN/Infinity leaking into a finite call. Every input must resolve to
/// EITHER a well-defined, theory-correct value OR a documented, intentional outcome
/// (ArgumentNullException for null depth, ArgumentOutOfRangeException for windowSize ≤ 0).
/// The headline hazards for THIS unit, which is built on a division and a log2, are:
///   • zero coverage → a zero-depth window has an undefined ratio (log2 0 = −∞); it MUST be
///     dropped as a no-call, never yield a −∞ / NaN copy number, and never divide by zero
///     (§3.3, §6.1, INV-06);
///   • a single bin / windowSize ≥ length → exactly one (or zero) full window; the trailing
///     partial window is dropped, with no off-by-one or empty-window crash (§6.1, INV-05);
///   • ratio = 1 (RD_w = RD_ref) → log2 = 0 → CN = 2 (neutral); the round-to-even of
///     2·2^0 = 2.0 must land on 2, not 1 or 3 (§6.1 "Window mean = reference", INV-01);
///   • all-zero depth → the median baseline is 0, the reference is unusable ⇒ NO segments,
///     never a divide-by-zero on RD_ref (§3.3, §4.1[2]);
///   • extreme depths (0 / MaxInt) → the per-window sum is accumulated in a <c>long</c>; a
///     MaxInt-filled window must not overflow into a negative mean (Boundary Exploitation
///     MaxInt; §4.1[1]).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SV-CNV-001 — Read-depth copy number variation detection (StructuralVar)
/// Checklist: docs/checklists/03_FUZZING.md, row 202.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 202): "ratio=1, zero coverage, single bin".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (Copy_Number_Variation.md)
/// (docs/algorithms/StructuralVar/Copy_Number_Variation.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • partition depth into floor(n / windowSize) non-overlapping FULL windows; a trailing
///     partial window is dropped                                          (§4.1[1], INV-05)
///   • RD_w = mean of per-position depths inside the window                (§2.2, §4.1[1])
///   • reference depth = supplied value, else overall median of NON-ZERO window means;
///     reference ≤ 0 ⇒ no segments                                  (§3.1, §4.1[2], ASM-01)
///   • log2 = log2(RD_w / RD_ref); CN = max(0, round(2·2^log2))      (§2.2, INV-02, INV-03)
///   • RD_w = RD_ref ⇒ log2 = 0 ⇒ CN = 2 (diploid neutral)                (INV-01, §6.1)
///   • a zero-depth window is a no-call (excluded), not a finite call      (INV-06, §6.1)
///   • CN is monotonically non-decreasing in RD_w                          (INV-04)
///   • defaults: windowSize 100, ploidy 2                                  (§3.1, §4.2)
///   • null depth/logRatios ⇒ ArgumentNullException; windowSize ≤ 0 ⇒
///     ArgumentOutOfRangeException                                         (§3.3, §6.1)
///   • SegmentCopyNumber: NaN log2 ratios are no-calls (dropped); adjacent equal-CN windows
///     merge into one segment (mean log2, window count as probe count)     (§4.1[4], §3.3)
///   • worked example: ref 100, windowSize 4, depth {100×4, 50×4, 150×4} ⇒ window means
///     100/50/150 ⇒ log2 0 / −1 / 0.585 ⇒ CN 2 / 1 / 3                     (§7.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class StructuralVariantCopyNumberFuzzTests
{
    // Documented diploid baseline and literature default window size (§4.2), mirrored
    // locally so the test owns the conversion rule independently of the source constants.
    private const int DiploidPloidy = 2;

    // Independent re-derivation of the documented call rule (§2.2, INV-02/INV-03):
    //   CN = max(0, round_to_even(ploidy · 2^log2)). round-to-even matches CNVkit's NumPy
    // ndarray.round() at half-integers (§5.x note in the doc). This is computed from the
    // spec formula, NOT read off the implementation's output.
    private static int ExpectedCopyNumber(double log2)
    {
        double copies = DiploidPloidy * Math.Pow(2, log2);
        int rounded = (int)Math.Round(copies, MidpointRounding.ToEven);
        return Math.Max(0, rounded);
    }

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural invariants on EVERY accepted segment set, regardless
    // of input. Each emitted segment must carry a finite log2 ratio (no −∞/NaN from a
    // zero-depth window), a non-negative copy number that equals the spec conversion of its
    // own log2 ratio (INV-02/INV-03), and the field conventions of a depth-only caller
    // (BAF = NaN by design, §3.2). This is what stops a fuzz test from rubber-stamping a
    // green call.
    private static void AssertWellFormed(IReadOnlyList<CopyNumberSegment> segments)
    {
        foreach (var seg in segments)
        {
            double.IsNaN(seg.LogRatio).Should().BeFalse("a no-call window is excluded, never emitted with NaN log2 (INV-06)");
            double.IsInfinity(seg.LogRatio).Should().BeFalse("log2 of a positive mean depth is finite (INV-06)");

            seg.CopyNumber.Should().BeGreaterThanOrEqualTo(0, "copy number is clamped non-negative (INV-03)");
            seg.CopyNumber.Should().Be(ExpectedCopyNumber(seg.LogRatio),
                "every segment's CN must equal the spec conversion round(2·2^log2) of its own log2 ratio (INV-02)");

            double.IsNaN(seg.BAlleleFrequency).Should().BeTrue("a depth-only caller has no allele-specific data ⇒ BAF = NaN (§3.2)");
            seg.End.Should().BeGreaterThanOrEqualTo(seg.Start, "window/segment bounds are ordered (§3.2)");
            seg.ProbeCount.Should().BeGreaterThan(0, "a window/merged run covers at least one position/window (§3.2)");
        }
    }

    #region SV-CNV-001 — Positive sanity (documented worked example)

    [Test]
    public void DetectCNV_DocumentedWorkedExample_NeutralLossGain()
    {
        // Docs §7.1 (hand-checked): reference 100, windowSize 4, depth = {100×4, 50×4,
        // 150×4}. Window means 100/50/150 ⇒ log2(1.0)=0, log2(0.5)=−1, log2(1.5)=0.585 ⇒
        // CN round(2·1)=2 (neutral), round(2·0.5)=1 (deletion), round(2·1.5)=3 (duplication).
        int[] depth =
        {
            100, 100, 100, 100,
             50,  50,  50,  50,
            150, 150, 150, 150,
        };

        var segments = DetectCNV(depth, windowSize: 4, referenceDepth: 100).ToList();

        segments.Should().HaveCount(3, "three full windows of 4 positions each (§7.1)");

        segments[0].LogRatio.Should().BeApproximately(0.0, 1e-12, "log2(100/100) = 0 (§7.1)");
        segments[0].CopyNumber.Should().Be(2, "round(2·2^0) = 2 — neutral (INV-01)");

        segments[1].LogRatio.Should().BeApproximately(-1.0, 1e-12, "log2(50/100) = −1 (§7.1)");
        segments[1].CopyNumber.Should().Be(1, "round(2·2^−1) = round(1.0) = 1 — single-copy loss (§7.1)");

        segments[2].LogRatio.Should().BeApproximately(Math.Log2(1.5), 1e-12, "log2(150/100) = log2(1.5) ≈ 0.585 (§7.1)");
        segments[2].CopyNumber.Should().Be(3, "round(2·1.5) = 3 — single-copy gain (§7.1)");

        // Coordinates: inclusive 0-based window bounds; one segment per window (§3.2).
        segments[0].Start.Should().Be(0);
        segments[0].End.Should().Be(3);
        segments[2].Start.Should().Be(8);
        segments[2].End.Should().Be(11);
        segments.Should().OnlyContain(s => s.ProbeCount == 4, "ProbeCount = positions per window (§3.2)");
        AssertWellFormed(segments);
    }

    [Test]
    public void DetectCNV_DefaultMedianReference_NeutralAtOverallMedian()
    {
        // ASM-01: with no explicit reference, the anchor is the overall median of non-zero
        // window means. Window means {40, 80, 80, 160} ⇒ median = (80+80)/2 = 80. The two
        // 80-windows then have RD_w = RD_ref ⇒ log2 0 ⇒ CN 2; 40 ⇒ log2(0.5)=−1 ⇒ CN 1;
        // 160 ⇒ log2(2)=1 ⇒ CN round(2·2)=4. Hand-checked, independent of the code.
        int[] depth =
        {
             40,  40,
             80,  80,
             80,  80,
            160, 160,
        };

        var segments = DetectCNV(depth, windowSize: 2).ToList(); // referenceDepth null ⇒ median

        segments.Should().HaveCount(4);
        segments[0].CopyNumber.Should().Be(1, "40/80 ⇒ log2(0.5) ⇒ CN 1");
        segments[1].CopyNumber.Should().Be(2, "80/80 ⇒ log2 0 ⇒ CN 2 (INV-01)");
        segments[2].CopyNumber.Should().Be(2, "80/80 ⇒ log2 0 ⇒ CN 2 (INV-01)");
        segments[3].CopyNumber.Should().Be(4, "160/80 ⇒ log2 1 ⇒ round(2·2) = 4");
        AssertWellFormed(segments);
    }

    [Test]
    public void SegmentCopyNumber_MergesAdjacentEqualCopyNumberRuns()
    {
        // §4.1[4]: log2 ratios {0,0,0, −1, 1,1} ⇒ CN {2,2,2, 1, 4,4} ⇒ three merged
        // segments (a run of three CN-2, one CN-1, a run of two CN-4). LogRatio of a merged
        // run is the mean of its members; ProbeCount is the window count. Hand-checked.
        var logRatios = new[] { 0.0, 0.0, 0.0, -1.0, 1.0, 1.0 };

        var segments = SegmentCopyNumber(logRatios).ToList();

        segments.Should().HaveCount(3, "three maximal equal-CN runs (§4.1[4])");
        segments[0].CopyNumber.Should().Be(2);
        segments[0].ProbeCount.Should().Be(3, "three CN-2 windows merged");
        segments[0].LogRatio.Should().BeApproximately(0.0, 1e-12, "mean of {0,0,0}");
        segments[1].CopyNumber.Should().Be(1);
        segments[1].ProbeCount.Should().Be(1);
        segments[2].CopyNumber.Should().Be(4);
        segments[2].ProbeCount.Should().Be(2, "two CN-4 windows merged");
        AssertWellFormed(segments);
    }

    #endregion

    #region SV-CNV-001 — BE boundary: ratio = 1 (neutral, no CNV)

    [Test]
    public void DetectCNV_RatioOne_YieldsNeutralCopyNumberTwo()
    {
        // "ratio = 1" boundary (checklist row 202): RD_w exactly equal to RD_ref ⇒ log2 = 0
        // ⇒ CN = round(2·2^0) = round(2.0) = 2, the diploid-neutral call (INV-01, §6.1
        // "Window mean = reference"). Round-to-even of an exact 2.0 must land on 2.
        int[] depth = Enumerable.Repeat(75, 12).ToArray();

        var segments = DetectCNV(depth, windowSize: 4, referenceDepth: 75).ToList();

        segments.Should().HaveCount(3, "12 positions / window 4 = 3 full windows (INV-05)");
        segments.Should().OnlyContain(s => s.LogRatio == 0.0, "RD_w = RD_ref ⇒ log2 0 (INV-01)");
        segments.Should().OnlyContain(s => s.CopyNumber == 2, "log2 0 ⇒ neutral CN 2 (INV-01)");
        AssertWellFormed(segments);

        // SegmentCopyNumber sees the same neutral signal: log2 0 everywhere ⇒ ONE merged
        // CN-2 segment spanning all windows (no spurious gain/loss called at ratio 1).
        var merged = SegmentCopyNumber(Enumerable.Repeat(0.0, 5)).ToList();
        merged.Should().ContainSingle("equal CN merges").Which.CopyNumber.Should().Be(2);
    }

    #endregion

    #region SV-CNV-001 — BE boundary: zero coverage (no divide-by-zero / no −∞)

    [Test]
    public void DetectCNV_ZeroDepthWindow_IsNoCall_Excluded()
    {
        // "zero coverage" boundary (checklist row 202): a window whose mean depth is 0 has
        // an undefined ratio (log2 0 = −∞). It MUST be dropped as a no-call (INV-06), never
        // emitted with a −∞/NaN log2. Here window 1 (positions 4..7) is all-zero; windows 0
        // and 2 are valid neutral/gain calls against ref 50.
        int[] depth =
        {
            50, 50, 50, 50,   // mean 50  → log2 0   → CN 2
             0,  0,  0,  0,   // mean 0   → no-call (excluded, INV-06)
           100,100,100,100,   // mean 100 → log2 1   → CN 4
        };

        var segments = DetectCNV(depth, windowSize: 4, referenceDepth: 50).ToList();

        segments.Should().HaveCount(2, "the zero-depth window is excluded as a no-call (INV-06)");
        segments.Should().OnlyContain(s => s.Start != 4, "no segment is emitted for the zero-depth window");
        segments.Select(s => s.CopyNumber).Should().Equal(2, 4);
        AssertWellFormed(segments);
    }

    [Test]
    public void DetectCNV_AllZeroDepth_YieldsNoSegments_NoDivideByZero()
    {
        // Zero coverage everywhere: every window mean is 0 ⇒ the median of non-zero means is
        // 0 ⇒ reference ≤ 0 ⇒ NO segments (§3.3, §4.1[2]). The division RD_w/RD_ref must
        // never be reached, so there is no divide-by-zero and no NaN segment.
        int[] depth = Enumerable.Repeat(0, 20).ToArray();

        var segments = DetectCNV(depth, windowSize: 5).ToList(); // null reference ⇒ median 0

        segments.Should().BeEmpty("reference depth 0 (all-zero windows) ⇒ no calls (§3.3)");
    }

    [Test]
    public void DetectCNV_NonPositiveExplicitReference_YieldsNoSegments()
    {
        // A supplied reference ≤ 0 is unusable (§3.1: "≤ 0 reference ⇒ no calls"). Both 0
        // and a negative reference must short-circuit to an empty result, never a log2 of a
        // negative ratio (which would be NaN) leaking into a segment.
        int[] depth = Enumerable.Repeat(30, 8).ToArray();

        DetectCNV(depth, windowSize: 4, referenceDepth: 0).Should().BeEmpty("zero reference ⇒ no calls (§3.1)");
        DetectCNV(depth, windowSize: 4, referenceDepth: -10).Should().BeEmpty("negative reference ⇒ no calls (§3.1)");
    }

    [Test]
    public void SegmentCopyNumber_NaNLogRatios_AreNoCalls_Dropped()
    {
        // The SegmentCopyNumber analogue of a zero-coverage no-call: a NaN log2 ratio is
        // unusable and is dropped, breaking the current run (§3.3). Ratios {0,0, NaN, 0}
        // ⇒ a CN-2 run, then the NaN, then a fresh CN-2 run ⇒ TWO segments, none NaN.
        var logRatios = new[] { 0.0, 0.0, double.NaN, 0.0 };

        var segments = SegmentCopyNumber(logRatios).ToList();

        segments.Should().HaveCount(2, "the NaN no-call breaks the run into two CN-2 segments (§3.3)");
        segments.Should().OnlyContain(s => s.CopyNumber == 2);
        AssertWellFormed(segments);
    }

    #endregion

    #region SV-CNV-001 — BE boundary: single bin / windowSize vs length

    [Test]
    public void DetectCNV_SingleBin_WindowSizeEqualsLength_OneSegment()
    {
        // "single bin" boundary (checklist row 202): windowSize == length ⇒ exactly ONE
        // full window spanning the whole track. Mean of {10,20,30,40} = 25; ref 25 ⇒
        // log2 0 ⇒ CN 2. No off-by-one, no second partial window.
        int[] depth = { 10, 20, 30, 40 };

        var segments = DetectCNV(depth, windowSize: 4, referenceDepth: 25).ToList();

        segments.Should().ContainSingle("windowSize = length ⇒ one full window (INV-05)");
        segments[0].Start.Should().Be(0);
        segments[0].End.Should().Be(3, "inclusive bound of the single window");
        segments[0].LogRatio.Should().BeApproximately(0.0, 1e-12, "mean 25 / ref 25 ⇒ log2 0");
        segments[0].CopyNumber.Should().Be(2);
        segments[0].ProbeCount.Should().Be(4);
        AssertWellFormed(segments);
    }

    [Test]
    public void DetectCNV_TrailingPartialWindow_IsDropped()
    {
        // INV-05: a trailing partial window is dropped. 10 positions / windowSize 4 ⇒
        // floor(10/4) = 2 full windows; positions 8..9 form no window. Pins that only full
        // windows are summarised (the partial tail does not appear).
        int[] depth = Enumerable.Repeat(60, 10).ToArray();

        var segments = DetectCNV(depth, windowSize: 4, referenceDepth: 60).ToList();

        segments.Should().HaveCount(2, "floor(10/4) = 2 full windows; partial tail dropped (INV-05)");
        segments.Last().End.Should().Be(7, "the last full window ends at position 7, not 9");
        AssertWellFormed(segments);
    }

    [Test]
    public void DetectCNV_WindowSizeGreaterThanLength_YieldsEmpty()
    {
        // windowSize > length ⇒ floor(n/windowSize) = 0 ⇒ no full window ⇒ empty (§6.1,
        // INV-05). No partial-window crash.
        int[] depth = { 5, 5, 5 };

        DetectCNV(depth, windowSize: 100, referenceDepth: 5)
            .Should().BeEmpty("no full window fits ⇒ empty (INV-05)");
    }

    [Test]
    public void SegmentCopyNumber_SingleLogRatio_OneSegment()
    {
        // Single bin for the segmentation variant: one log2 ratio ⇒ one segment covering
        // window index 0. log2 1 ⇒ CN round(2·2) = 4.
        var segments = SegmentCopyNumber(new[] { 1.0 }).ToList();

        segments.Should().ContainSingle("one log2 ratio ⇒ one segment");
        segments[0].Start.Should().Be(0);
        segments[0].End.Should().Be(0, "single-window index range [0,0]");
        segments[0].CopyNumber.Should().Be(4, "log2 1 ⇒ round(2·2) = 4");
        segments[0].ProbeCount.Should().Be(1);
        AssertWellFormed(segments);
    }

    #endregion

    #region SV-CNV-001 — BE boundary: empty / null / invalid windowSize / extremes

    [Test]
    public void DetectCNV_EmptyDepth_YieldsEmpty()
    {
        // Empty depth ⇒ no window can be formed ⇒ empty (§6.1). No throw on an empty track.
        DetectCNV(Array.Empty<int>(), windowSize: 4, referenceDepth: 10)
            .Should().BeEmpty("no positions ⇒ no windows (§6.1)");

        SegmentCopyNumber(Array.Empty<double>())
            .Should().BeEmpty("no log2 ratios ⇒ no segments (§6.1)");
    }

    [Test]
    public void DetectCNV_NullDepth_ThrowsArgumentNullException()
    {
        // §3.3 / §6.1: null depth ⇒ ArgumentNullException, thrown eagerly (the public method
        // validates before deferring to the iterator).
        Action act = () => DetectCNV(null!, windowSize: 4);
        act.Should().Throw<ArgumentNullException>();

        Action seg = () => SegmentCopyNumber((IEnumerable<double>)null!);
        seg.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void DetectCNV_NonPositiveWindowSize_ThrowsArgumentOutOfRange()
    {
        // §3.3 / §6.1: windowSize ≤ 0 has no meaningful window ⇒ ArgumentOutOfRangeException,
        // thrown eagerly. Covers the 0 and −1 boundaries (BE).
        int[] depth = { 1, 2, 3, 4 };

        ((Action)(() => DetectCNV(depth, windowSize: 0))).Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => DetectCNV(depth, windowSize: -1))).Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void DetectCNV_MaxIntDepth_NoOverflow_FiniteCall()
    {
        // MaxInt boundary: a window filled with int.MaxValue must accumulate in a long
        // without overflowing into a negative mean. With ref = int.MaxValue, RD_w = RD_ref
        // ⇒ log2 0 ⇒ CN 2; the mean must be exactly int.MaxValue, not a wrapped value.
        int[] depth = Enumerable.Repeat(int.MaxValue, 4).ToArray();

        var segments = DetectCNV(depth, windowSize: 4, referenceDepth: int.MaxValue).ToList();

        segments.Should().ContainSingle("one full window");
        segments[0].LogRatio.Should().BeApproximately(0.0, 1e-9,
            "mean of four int.MaxValue = int.MaxValue (long accumulation, no overflow)");
        segments[0].CopyNumber.Should().Be(2, "RD_w = RD_ref ⇒ neutral (INV-01)");
        AssertWellFormed(segments);
    }

    #endregion

    #region SV-CNV-001 — Randomized BE sweep

    [Test]
    [CancelAfter(30_000)]
    public void DetectCNV_RandomizedBoundarySweep_NeverThrows_WellFormed()
    {
        // Randomized BE sweep over degenerate depth tracks: empty / single-window sizes,
        // windowSize ≥ length, zero-coverage windows, MaxInt-adjacent and identical depths,
        // explicit and median references. Every accepted result must satisfy the documented
        // invariants and the call must never throw or hang (CancelAfter 30s). The expected
        // copy number per window is recomputed from the spec formula (§2.2/§4.1), NOT echoed
        // off the implementation.
        var rng = new Random(202_001);
        for (int trial = 0; trial < 4000; trial++)
        {
            int n = rng.Next(0, 40);                 // includes empty
            int windowSize = rng.Next(1, 12);        // positive only (≤0 is the throw test)
            double? reference = rng.Next(0, 3) switch
            {
                0 => null,                            // median baseline
                1 => rng.Next(1, 200),                // explicit positive
                _ => 0.0,                             // unusable reference ⇒ no calls
            };

            var depth = new int[n];
            for (int i = 0; i < n; i++)
            {
                depth[i] = rng.Next(0, 6) switch
                {
                    0 => 0,                           // zero-coverage bait (no-call windows)
                    1 => int.MaxValue - rng.Next(0, 4), // MaxInt neighbourhood (overflow bait)
                    2 => 100,                         // identical-depth bait
                    _ => rng.Next(0, 500),
                };
            }

            List<CopyNumberSegment> segments = null!;
            Action act = () => segments = DetectCNV(depth, windowSize, reference, "chrZ").ToList();
            act.Should().NotThrow($"trial {trial}: degenerate depth track must not throw");

            AssertWellFormed(segments);

            // Independent recomputation of the contract (§4.1): floor(n/windowSize) full
            // windows; mean per window (long sum); reference = explicit, else overall median
            // of non-zero means; emit a segment per positive-mean window when reference > 0.
            int windowCount = n / windowSize;
            var means = new double[windowCount];
            for (int w = 0; w < windowCount; w++)
            {
                long sum = 0;
                for (int j = 0; j < windowSize; j++)
                    sum += depth[w * windowSize + j];
                means[w] = (double)sum / windowSize;
            }

            double refDepth;
            if (reference.HasValue)
            {
                refDepth = reference.Value;
            }
            else
            {
                var pos = means.Where(m => m > 0).OrderBy(m => m).ToList();
                refDepth = pos.Count == 0
                    ? 0
                    : pos.Count % 2 == 1
                        ? pos[pos.Count / 2]
                        : (pos[pos.Count / 2 - 1] + pos[pos.Count / 2]) / 2.0;
            }

            var expectedCns = new List<int>();
            if (refDepth > 0)
            {
                foreach (double m in means)
                {
                    if (m <= 0) continue; // zero-depth no-call (INV-06)
                    expectedCns.Add(ExpectedCopyNumber(Math.Log2(m / refDepth)));
                }
            }

            segments.Select(s => s.CopyNumber).Should().Equal(expectedCns,
                $"trial {trial}: per-window copy numbers must match the spec conversion (§4.1)");
        }
    }

    #endregion
}
