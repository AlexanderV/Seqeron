using Enc = Seqeron.Genomics.IO.QualityScoreAnalyzer.QualityEncoding;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Quality area — QUALITY-STATS-001 (Phred summary statistics).
/// The unit under test is <see cref="QualityScoreAnalyzer.CalculateStatistics(string,
/// QualityScoreAnalyzer.QualityEncoding)"/> — descriptive summary statistics over the
/// per-base Phred scores of a FASTQ quality line (mean, median, min, max, population
/// standard deviation, %≥Q20, %≥Q30) — and its canonical companion
/// <see cref="QualityScoreAnalyzer.CalculateQ30Percentage"/>; implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no NaN/Infinity, no nonsense output, and no *unhandled* runtime
/// exception. Every input must resolve to a well-defined, theory-correct value.
/// For this unit the documented contract is *total* on strings: null / empty
/// yield a zeroed result (TotalBases = 0) rather than throwing or dividing by
/// zero — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing",
/// docs/algorithms/Quality/Quality_Statistics.md §3.3, §6.1.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: QUALITY-STATS-001 — Quality Statistics (Quality)
/// Checklist: docs/checklists/03_FUZZING.md, row 220.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// MAPPING of the checklist's "empty, single base, all-equal quality" targets
/// onto THIS unit's documented contract
/// (docs/algorithms/Quality/Quality_Statistics.md):
///   • "empty" (BE) → null / empty quality string ⇒ the all-zero
///       QualityStatistics record (TotalBases = 0); percentages are 0 by
///       contract, NOT 0/0 = NaN, and no divide-by-zero (§3.3, §6.1, §6.2).
///   • "single base" (BE) → one score ⇒ mean = median = min = max = that score
///       and σ = 0 (zero spread): degenerate but defined, no NaN (§6.1).
///   • "all-equal quality" (BE) → every base the same Q ⇒ mean = median = that
///       value, σ = 0 (INV-02: σ = 0 iff all scores equal), %≥thresholds is 0
///       or 100 according to whether the common value clears the threshold.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (docs/algorithms/Quality/Quality_Statistics.md)
/// Independently re-derived from the cited primary sources, NOT read off the code:
/// ───────────────────────────────────────────────────────────────────────────
///   • Scores qᵢ = ord(char) − offset (offset 33 / 64) before any statistic (§2.2).
///   • Mean μ = (1/N) Σ qᵢ — the ARITHMETIC mean of the decoded Phred scores
///       (the samtools-stats / FastQC "average quality"), NOT the
///       probability-based −10·log₁₀(mean P_error). This distinction is the
///       crux of §6.2 and is pinned explicitly below.
///   • Median = middle of sorted scores (odd N) / mean of the two central order
///       statistics (even N) (§2.2 [5]).
///   • Population σ = √((1/N) Σ (qᵢ − μ)²) — divisor N, not N−1 (§2.2 [6]).
///   • %≥Q20 = 100·|{qᵢ ≥ 20}|/N, %≥Q30 = 100·|{qᵢ ≥ 30}|/N, thresholds
///       INCLUSIVE (≥) (§2.2 [2][3]).
///   • INV-01 min ≤ mean ≤ max (non-empty); INV-02 σ ≥ 0, σ = 0 iff all equal;
///       INV-03 0 ≤ %≥Q30 ≤ %≥Q20 ≤ 100; INV-04 CalculateQ30Percentage ==
///       PercentAboveQ30; INV-05 statistics depend only on the decoded scores.
///   • null / empty ⇒ zeroed record / 0.0 (§3.3, §6.1), no throw.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class QualityStatsFuzzTests
{
    // Q20/Q30 thresholds, independently fixed from Illumina / Ewing & Green (1998) —
    // NOT echoed from the implementation's private constants.
    private const int Q20 = 20;
    private const int Q30 = 30;
    private const int Phred33Offset = 33;
    private const int Phred64Offset = 64;

    // ── Independent oracles, computed directly from §2.2 definitions. Deliberately
    //    written WITHOUT reusing the analyzer, so a wrong implementation cannot make
    //    the test pass by agreeing with itself. ──────────────────────────────────
    private static double ExpectedMean(int[] q) => q.Sum() / (double)q.Length;

    private static double ExpectedMedian(int[] q)
    {
        int[] s = q.OrderBy(x => x).ToArray();
        int n = s.Length;
        return n % 2 == 0 ? (s[n / 2 - 1] + s[n / 2]) / 2.0 : s[n / 2];
    }

    private static double ExpectedPopulationStdDev(int[] q)
    {
        double mean = ExpectedMean(q);
        double var = q.Sum(x => (x - mean) * (x - mean)) / q.Length;   // ÷N (population)
        return Math.Sqrt(var);
    }

    private static double ExpectedPercentAtLeast(int[] q, int threshold)
        => 100.0 * q.Count(x => x >= threshold) / q.Length;

    private static string Encode(int[] q, int offset)
        => new string(q.Select(x => (char)(x + offset)).ToArray());

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-STATS-001 — positive sanity (hand-checkable worked example)
    // ═════════════════════════════════════════════════════════════════════════

    // ── POSITIVE sanity: the doc's worked example (§7.1). "5?I" under Phred+33
    //    decodes to [20, 30, 40] (5=53−33, ?=63−33, I=73−33). Hand-computed:
    //      μ = (20+30+40)/3 = 30
    //      deviations −10,0,10 ⇒ σ² = (100+0+100)/3 = 200/3 ⇒ σ = √(200/3)
    //                                                          = 8.16496580927726
    //      median = 30 (middle of sorted [20,30,40])
    //      %≥Q30 = 100·2/3 = 66.66666666666667 (30 and 40 clear Q30)
    //      %≥Q20 = 100      (all three clear Q20)
    //    Every value is taken from the doc, not from the implementation. ──
    [Test]
    public void CalculateStatistics_DocWorkedExample_HandComputed()
    {
        var s = QualityScoreAnalyzer.CalculateStatistics("5?I", Enc.Phred33);

        s.TotalBases.Should().Be(3);
        s.MeanQuality.Should().BeApproximately(30.0, 1e-12, "μ = (20+30+40)/3 = 30 (§7.1)");
        s.MedianQuality.Should().Be(30.0, "median of sorted [20,30,40] (§2.2)");
        s.MinQuality.Should().Be(20);
        s.MaxQuality.Should().Be(40);
        s.StandardDeviation.Should().BeApproximately(8.16496580927726, 1e-12,
            "population σ = √(200/3) (÷N) per §7.1");
        s.BasesAboveQ30.Should().Be(2);
        s.PercentAboveQ30.Should().BeApproximately(66.66666666666667, 1e-12,
            "100·2/3 — Q30 inclusive (§2.2 [2][3], §7.1)");
        s.PercentAboveQ20.Should().BeApproximately(100.0, 1e-12, "all three ≥ Q20");

        // INV-04: the canonical Q30 helper agrees exactly with the record field.
        QualityScoreAnalyzer.CalculateQ30Percentage("5?I", Enc.Phred33)
            .Should().Be(s.PercentAboveQ30, "CalculateQ30Percentage == PercentAboveQ30 (INV-04)");
    }

    // ── POSITIVE sanity, MEAN-DEFINITION CRUX (§6.2): MeanQuality is the ARITHMETIC
    //    mean of Q, NOT the probability-based −10·log₁₀(mean P_error). These two
    //    differ sharply when the spread is large, so we use a deliberately skewed
    //    quality line to separate them. For scores [2, 40]:
    //      arithmetic mean of Q          = (2 + 40)/2          = 21
    //      probability-based "true" mean = −10·log₁₀((10^−0.2 + 10^−4)/2)
    //                                    ≈ −10·log₁₀(0.31550)  ≈ 5.009  (≈ Q5, far lower)
    //    The documented metric is the arithmetic 21, NOT ≈5. Pinning this prevents a
    //    silent switch to the error-probability convention. ──
    [Test]
    public void CalculateStatistics_MeanIsArithmeticMeanOfQ_NotProbabilityBased()
    {
        int[] q = { 2, 40 };                       // ASCII '#' and 'I' under Phred+33
        string line = Encode(q, Phred33Offset);

        double arithmeticMean = ExpectedMean(q);   // = 21
        double pbar = q.Average(x => Math.Pow(10, -x / 10.0));
        double probabilityBasedMean = -10.0 * Math.Log10(pbar);   // ≈ 5.009

        arithmeticMean.Should().Be(21.0);
        probabilityBasedMean.Should().BeApproximately(5.009, 0.01,
            "fixture sanity: the probability-based mean is far below the arithmetic mean");

        var s = QualityScoreAnalyzer.CalculateStatistics(line, Enc.Phred33);
        s.MeanQuality.Should().BeApproximately(arithmeticMean, 1e-12,
            "MeanQuality is the arithmetic mean of Q (samtools/FastQC), not −10·log₁₀(mean P) (§6.2)");
        s.MeanQuality.Should().NotBeApproximately(probabilityBasedMean, 1.0,
            "the reported mean must NOT be the probability-based convention (§6.2)");
    }

    // ── POSITIVE sanity, EVEN N median rule (§2.2): median of an even-length set is
    //    the mean of the two central order statistics. Scores [10,20,30,40] sorted
    //    ⇒ centre pair (20,30) ⇒ median 25 (a value not present in the data). ──
    [Test]
    public void CalculateStatistics_EvenLength_MedianIsMeanOfCentralPair()
    {
        int[] q = { 40, 10, 30, 20 };              // unsorted on purpose
        var s = QualityScoreAnalyzer.CalculateStatistics(Encode(q, Phred33Offset), Enc.Phred33);

        s.MedianQuality.Should().Be(25.0, "median of even N = mean of central pair (20,30) (§2.2)");
        s.MedianQuality.Should().Be(ExpectedMedian(q));
    }

    // ── POSITIVE sanity, INV-05 (encoding-invariance): the SAME scores written under
    //    Phred+33 and Phred+64 must yield identical statistics, because the statistic
    //    depends only on the decoded score, not the byte. ──
    [Test]
    public void CalculateStatistics_EncodingInvariant_SameScoresSameStats()
    {
        int[] q = { 5, 18, 27, 33, 41 };
        var s33 = QualityScoreAnalyzer.CalculateStatistics(Encode(q, Phred33Offset), Enc.Phred33);
        var s64 = QualityScoreAnalyzer.CalculateStatistics(Encode(q, Phred64Offset), Enc.Phred64);

        s64.MeanQuality.Should().Be(s33.MeanQuality, "statistics are encoding-invariant (INV-05)");
        s64.MedianQuality.Should().Be(s33.MedianQuality);
        s64.StandardDeviation.Should().BeApproximately(s33.StandardDeviation, 1e-12);
        s64.PercentAboveQ20.Should().Be(s33.PercentAboveQ20);
        s64.PercentAboveQ30.Should().Be(s33.PercentAboveQ30);
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-STATS-001 — BE: empty / null input (no NaN, no divide-by-zero)
    // ═════════════════════════════════════════════════════════════════════════

    // ── BE empty/null: the documented contract is a ZEROED record (TotalBases = 0),
    //    NOT an exception and NOT 0/0 = NaN for the percentages (§3.3, §6.1, §6.2).
    //    The crux of the boundary is that the percentages must be a finite 0, never
    //    NaN/Infinity from a 0/0 division. ──
    [TestCase("", Enc.Phred33)]
    [TestCase("", Enc.Phred64)]
    [TestCase("", Enc.Auto)]
    [TestCase(null, Enc.Phred33)]
    public void CalculateStatistics_NullOrEmpty_ReturnsZeroedResult_NoNaN(string? input, Enc enc)
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(input!, enc);

        s.TotalBases.Should().Be(0, "no observations ⇒ zeroed result (§3.3, §6.1)");
        s.MeanQuality.Should().Be(0.0);
        s.MedianQuality.Should().Be(0.0);
        s.MinQuality.Should().Be(0);
        s.MaxQuality.Should().Be(0);
        s.StandardDeviation.Should().Be(0.0);
        s.BasesAboveQ20.Should().Be(0);
        s.BasesAboveQ30.Should().Be(0);

        // The headline boundary: percentages are a defined 0, NOT 0/0 = NaN (§6.2).
        s.PercentAboveQ20.Should().Be(0.0).And.NotBe(double.NaN);
        s.PercentAboveQ30.Should().Be(0.0);
        double.IsNaN(s.PercentAboveQ20).Should().BeFalse("0/0 must not surface as NaN");
        double.IsNaN(s.PercentAboveQ30).Should().BeFalse();
        s.PerPositionMeanQuality.Should().BeEmpty();
    }

    // ── BE empty/null on the canonical Q30 helper: 0.0, never NaN (§3.3, INV-04). ──
    [TestCase("")]
    [TestCase((string?)null)]
    public void CalculateQ30Percentage_NullOrEmpty_ReturnsZero_NoNaN(string? input)
    {
        double q30 = QualityScoreAnalyzer.CalculateQ30Percentage(input!, Enc.Phred33);
        q30.Should().Be(0.0, "empty/null ⇒ 0.0 by contract (§3.3)");
        double.IsNaN(q30).Should().BeFalse();
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-STATS-001 — BE: single base (degenerate, σ = 0, no NaN)
    // ═════════════════════════════════════════════════════════════════════════

    // ── BE single base: one score ⇒ mean = median = min = max = that score, σ = 0
    //    (zero spread), and the percentages are exactly 0 or 100 (§6.1). Swept over
    //    the boundary Q values 0 (worst), 20, 30 (threshold edges), 40, 93 (max). ──
    [TestCase(0)]
    [TestCase(19)]    // just below Q20
    [TestCase(20)]    // exactly Q20 — inclusive ⇒ counts toward %≥Q20
    [TestCase(29)]    // just below Q30
    [TestCase(30)]    // exactly Q30 — inclusive ⇒ counts toward %≥Q30
    [TestCase(40)]
    [TestCase(93)]    // Phred+33 maximum
    public void CalculateStatistics_SingleBase_AllCentralMomentsEqualScore_StdDevZero(int score)
    {
        string line = Encode(new[] { score }, Phred33Offset);
        var s = QualityScoreAnalyzer.CalculateStatistics(line, Enc.Phred33);

        s.TotalBases.Should().Be(1);
        s.MeanQuality.Should().Be(score, "single observation ⇒ mean = the score (§6.1)");
        s.MedianQuality.Should().Be(score, "single observation ⇒ median = the score");
        s.MinQuality.Should().Be(score);
        s.MaxQuality.Should().Be(score);
        s.StandardDeviation.Should().Be(0.0, "zero spread ⇒ σ = 0, not NaN (INV-02, §6.1)");
        double.IsNaN(s.StandardDeviation).Should().BeFalse();

        // Inclusive thresholds: a base exactly at the threshold is counted (§2.2 [2][3]).
        s.PercentAboveQ20.Should().Be(score >= Q20 ? 100.0 : 0.0);
        s.PercentAboveQ30.Should().Be(score >= Q30 ? 100.0 : 0.0);

        // INV-01 holds trivially as equalities; INV-03 holds.
        s.PercentAboveQ30.Should().BeLessThanOrEqualTo(s.PercentAboveQ20, "INV-03");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-STATS-001 — BE: all-equal quality (σ = 0, no NaN)
    // ═════════════════════════════════════════════════════════════════════════

    // ── BE all-equal: every base the same Q over many lengths ⇒ mean = median = that
    //    value and σ = 0 exactly (INV-02: σ = 0 iff all scores equal). No spurious
    //    floating-point NaN/Infinity even at large N. Threshold edges 20 and 30 are
    //    included to pin the inclusive-comparison behaviour for the uniform case. ──
    [TestCase(0, 1)]
    [TestCase(20, 5)]
    [TestCase(30, 50)]
    [TestCase(35, 500)]
    [TestCase(93, 1000)]
    public void CalculateStatistics_AllEqualQuality_StdDevZero_MeanEqualsMedianEqualsValue(int value, int n)
    {
        string line = Encode(Enumerable.Repeat(value, n).ToArray(), Phred33Offset);
        var s = QualityScoreAnalyzer.CalculateStatistics(line, Enc.Phred33);

        s.TotalBases.Should().Be(n);
        s.MeanQuality.Should().Be(value, "uniform input ⇒ mean = the common value");
        s.MedianQuality.Should().Be(value, "uniform input ⇒ median = the common value");
        s.MinQuality.Should().Be(value);
        s.MaxQuality.Should().Be(value);
        s.StandardDeviation.Should().Be(0.0, "σ = 0 iff all scores equal (INV-02)");
        double.IsNaN(s.StandardDeviation).Should().BeFalse();
        double.IsInfinity(s.StandardDeviation).Should().BeFalse();

        s.PercentAboveQ20.Should().Be(value >= Q20 ? 100.0 : 0.0, "uniform ⇒ all-or-nothing, inclusive");
        s.PercentAboveQ30.Should().Be(value >= Q30 ? 100.0 : 0.0);
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-STATS-001 — Robustness: out-of-alphabet / degenerate bytes
    // ═════════════════════════════════════════════════════════════════════════

    // ── Robustness: CalculateStatistics decodes with the lenient codec (no per-char
    //    range validation), so an out-of-alphabet byte yields a possibly-negative or
    //    oversized score rather than throwing. The fuzz bar here is: the call must
    //    not CRASH or hang and must produce a finite, self-consistent record
    //    (no NaN/Infinity in mean/σ/percentages), even on garbage characters. ──
    [TestCase(" QQQ")]      // leading space (ASCII 32) ⇒ a Q −1 under Phred+33
    [TestCase("\0\0\0")]    // NUL bytes ⇒ deeply negative scores
    [TestCase("ÿÿ")] // high bytes ⇒ oversized scores
    [TestCase("AaZz09~!")]  // mixed printable ASCII
    public void CalculateStatistics_OutOfAlphabetBytes_NoCrash_FiniteSelfConsistent(string line)
    {
        var s = QualityScoreAnalyzer.CalculateStatistics(line, Enc.Phred33);

        s.TotalBases.Should().Be(line.Length, "one decoded score per character");
        double.IsNaN(s.MeanQuality).Should().BeFalse();
        double.IsNaN(s.MedianQuality).Should().BeFalse();
        double.IsNaN(s.StandardDeviation).Should().BeFalse();
        double.IsInfinity(s.StandardDeviation).Should().BeFalse();
        s.StandardDeviation.Should().BeGreaterThanOrEqualTo(0.0, "σ ≥ 0 even on garbage (INV-02)");

        // Percentages stay within [0,100] and respect the subset ordering (INV-03).
        s.PercentAboveQ20.Should().BeInRange(0.0, 100.0);
        s.PercentAboveQ30.Should().BeInRange(0.0, 100.0);
        s.PercentAboveQ30.Should().BeLessThanOrEqualTo(s.PercentAboveQ20, "INV-03: {q≥30} ⊆ {q≥20}");

        // INV-01: min ≤ mean ≤ max must still hold on the decoded scores.
        s.MeanQuality.Should().BeGreaterThanOrEqualTo(s.MinQuality)
            .And.BeLessThanOrEqualTo(s.MaxQuality, "INV-01");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region QUALITY-STATS-001 — Randomized boundary sweep
    // ═════════════════════════════════════════════════════════════════════════

    // ── Randomized sweep over VALID quality lines: for many random in-range score
    //    arrays under a random encoding, every documented statistic must match the
    //    independent oracle (mean/median/σ/percentages), all invariants must hold,
    //    and no value may be NaN/Infinity. The sweep deliberately favours short and
    //    uniform lines so it keeps re-hitting the empty/single/all-equal boundaries.
    //    Locally seeded; bounded by [CancelAfter]. ──
    [Test]
    [CancelAfter(20_000)]
    public void RandomizedSweep_ValidLines_MatchOracleAndInvariants()
    {
        var rng = new Random(220_001);
        for (int iter = 0; iter < 5000; iter++)
        {
            bool phred64 = rng.Next(2) == 0;
            Enc enc = phred64 ? Enc.Phred64 : Enc.Phred33;
            int offset = phred64 ? Phred64Offset : Phred33Offset;
            int max = phred64 ? 62 : 93;

            // Bias toward the boundaries: ~⅓ empty/single, occasionally all-equal.
            int len = rng.Next(4) == 0 ? rng.Next(0, 2) : rng.Next(0, 30);
            bool uniform = rng.Next(3) == 0;
            int fixedVal = rng.Next(0, max + 1);

            var q = new int[len];
            for (int i = 0; i < len; i++)
                q[i] = uniform ? fixedVal : rng.Next(0, max + 1);

            string line = Encode(q, offset);
            var s = QualityScoreAnalyzer.CalculateStatistics(line, enc);

            s.TotalBases.Should().Be(len, "one score per character");

            if (len == 0)
            {
                // Empty boundary: zeroed and finite, no NaN.
                s.MeanQuality.Should().Be(0.0);
                s.StandardDeviation.Should().Be(0.0);
                s.PercentAboveQ20.Should().Be(0.0);
                s.PercentAboveQ30.Should().Be(0.0);
                double.IsNaN(s.PercentAboveQ30).Should().BeFalse();
                continue;
            }

            // Match the independent oracle exactly (mean/median) / closely (σ).
            s.MeanQuality.Should().BeApproximately(ExpectedMean(q), 1e-9, "μ = (1/N)Σqᵢ (§2.2)");
            s.MedianQuality.Should().BeApproximately(ExpectedMedian(q), 1e-9, "median rule (§2.2)");
            s.StandardDeviation.Should().BeApproximately(ExpectedPopulationStdDev(q), 1e-9,
                "population σ ÷N (§2.2 [6])");
            s.MinQuality.Should().Be(q.Min());
            s.MaxQuality.Should().Be(q.Max());
            s.PercentAboveQ20.Should().BeApproximately(ExpectedPercentAtLeast(q, Q20), 1e-9);
            s.PercentAboveQ30.Should().BeApproximately(ExpectedPercentAtLeast(q, Q30), 1e-9);

            // Invariants.
            s.MeanQuality.Should().BeGreaterThanOrEqualTo(s.MinQuality)
                .And.BeLessThanOrEqualTo(s.MaxQuality, "INV-01");
            s.StandardDeviation.Should().BeGreaterThanOrEqualTo(0.0, "INV-02");
            if (q.Distinct().Count() == 1)
                s.StandardDeviation.Should().Be(0.0, "σ = 0 iff all scores equal (INV-02)");
            s.PercentAboveQ30.Should().BeLessThanOrEqualTo(s.PercentAboveQ20, "INV-03");
            s.PercentAboveQ20.Should().BeInRange(0.0, 100.0);

            // Finiteness.
            double.IsNaN(s.StandardDeviation).Should().BeFalse();
            double.IsInfinity(s.StandardDeviation).Should().BeFalse();
            double.IsNaN(s.MeanQuality).Should().BeFalse();

            // INV-04: canonical helper agrees with the record.
            QualityScoreAnalyzer.CalculateQ30Percentage(line, enc)
                .Should().Be(s.PercentAboveQ30, "CalculateQ30Percentage == PercentAboveQ30 (INV-04)");
        }
    }

    #endregion
}
