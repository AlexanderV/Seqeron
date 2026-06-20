using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology tumour-ploidy area — ONCO-PLOIDY-001.
/// The unit under test is the deterministic ploidy entry point
/// <see cref="OncologyAnalyzer.EstimatePloidy(IEnumerable{OncologyAnalyzer.AlleleSpecificSegment})"/>,
/// implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output (NaN / Infinity / negative ploidy), and no *unhandled* runtime fault
/// (a DivideByZero on Σ Lᵢ = 0, or an integer-overflow-wrapped Σ CNᵢ·Lᵢ leaking
/// a wrong/negative ploidy). Every input must resolve to EITHER a well-defined,
/// theory-correct value OR a *documented, intentional* outcome (here, an
/// <see cref="ArgumentException"/> for an empty segment set / non-positive length /
/// negative copy number, or <see cref="ArgumentNullException"/> for null input).
/// For ploidy the headline hazards (checklist row 107, targets
/// "flat diploid, fully amplified, empty segments") are:
///   • flat diploid — every segment 1:1 (total CN 2) ⇒ ψ = 2.0 EXACTLY, the
///     canonical neutral 2n genome, irrespective of segment lengths (INV-02, §6.1);
///   • fully amplified — every segment at a high total CN over a huge total
///     length ⇒ ψ = that CN exactly, with NO integer overflow in Σ CNᵢ·Lᵢ even
///     when total length × CN dwarfs Int64.MaxValue (the weighted sum accumulates
///     in double, so an in-range mean stays exact and finite; INV-03);
///   • empty segments — no segments ⇒ documented ArgumentException, NEVER a
///     DivideByZero / NaN from Σ Lᵢ = 0 (§3.3, §6.1 "Empty segment set ⇒
///     ArgumentException — weighted mean undefined").
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-PLOIDY-001 — Tumour ploidy estimation (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 107.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 107): "flat diploid, fully amplified, empty segments".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Tumor_Ploidy_Estimation.md (docs/algorithms/Oncology/Tumor_Ploidy_Estimation.md):
///   • ψ = Σ(CNᵢ·Lᵢ) / Σ(Lᵢ), CNᵢ = Major+Minor, Lᵢ = End−Start    (§2.2)
///   • length-weighted mean ⇒ min CNᵢ ≤ ψ ≤ max CNᵢ                  (INV-03)
///   • all 1:1 segments ⇒ ψ = 2.0 exactly (2n diploid)               (INV-02, §6.1)
///   • ψ > 0 for any non-empty genome with ≥ 1 positive copy number  (INV-01)
///   • empty segment set ⇒ ArgumentException (Σ L = 0 undefined)      (§3.3, §6.1)
///   • segment End ≤ Start ⇒ ArgumentException                       (§6.1)
///   • negative copy number ⇒ ArgumentException                      (§6.1)
///   • null segments ⇒ ArgumentNullException                         (§3.3)
///   • worked example: 1:1@100Mb, 2:2@100Mb, 2:1@50Mb ⇒ ψ = 3.0      (§7.1)
///
/// SOURCE: no bug found. EstimatePloidy already accumulates the weighted copy-number
/// sum in `double` (no integer overflow in Σ CNᵢ·Lᵢ) and guards `totalLength == 0`
/// with an explicit ArgumentException (no DivideByZero / NaN). No test was weakened
/// and no source change was required.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyPloidyFuzzTests
{
    private const long Mb = 1_000_000L;

    // ── Well-formed-ploidy assertion helper ──────────────────────────────────
    // Pins the documented numeric contract on EVERY returned ψ: it must be a
    // FINITE, NON-negative real (INV-01). This is what stops a fuzz test from
    // rubber-stamping a NaN (0/0), an Infinity (x/0), or an overflow-wrapped
    // negative ploidy green.
    private static void AssertWellFormedPloidy(double ploidy)
    {
        double.IsNaN(ploidy).Should().BeFalse("ploidy ψ must never be NaN (no 0/0 on Σ L = 0)");
        double.IsInfinity(ploidy).Should().BeFalse("ploidy ψ must never be ±Infinity (no x/0)");
        ploidy.Should().BeGreaterThanOrEqualTo(0.0, "ψ is a length-weighted mean of non-negative CN (INV-01)");
    }

    private static AlleleSpecificSegment Seg(long start, long end, int major, int minor) =>
        new("1", start, end, major, minor);

    #region ONCO-PLOIDY-001 — Tumour ploidy estimation (positive sanity)

    // ── POSITIVE sanity: the §7.1 worked example, hand-computed exactly ───────
    // ψ = (2·100 + 4·100 + 3·50) Mb / (100+100+50) Mb = 750/250 = 3.0.
    [Test]
    public void EstimatePloidy_WorkedExample_HandComputedWeightedMean()
    {
        var segments = new[]
        {
            Seg(0, 100 * Mb, 1, 1), // total CN 2 over 100 Mb
            Seg(0, 100 * Mb, 2, 2), // total CN 4 over 100 Mb
            Seg(0, 50 * Mb, 2, 1),  // total CN 3 over 50 Mb
        };

        double ploidy = EstimatePloidy(segments);

        ploidy.Should().BeApproximately(3.0, 1e-12, "documented §7.1 walk-through: 750/250 = 3.0");
        AssertWellFormedPloidy(ploidy);
    }

    // ── POSITIVE sanity: UNEVEN lengths weight correctly (not a plain mean) ───
    // CN 2 over 1 Mb and CN 4 over 3 Mb ⇒ (2·1 + 4·3)/(1+3) = 14/4 = 3.5,
    // NOT the unweighted average (3.0). Pins that length actually weights.
    [Test]
    public void EstimatePloidy_UnevenLengths_WeightsByLengthNotCount()
    {
        var segments = new[]
        {
            Seg(0, 1 * Mb, 1, 1), // total CN 2 over 1 Mb
            Seg(0, 3 * Mb, 2, 2), // total CN 4 over 3 Mb
        };

        double ploidy = EstimatePloidy(segments);

        ploidy.Should().BeApproximately(3.5, 1e-12, "length-weighted mean 14/4, not the unweighted 3.0");
        AssertWellFormedPloidy(ploidy);
    }

    #endregion

    #region ONCO-PLOIDY-001 — BE: flat diploid (canonical 2n)

    // ── BE: every segment 1:1 ⇒ ψ = 2.0 EXACTLY, whatever the lengths ────────
    // INV-02 / §6.1 "All 1:1 segments ⇒ ψ = 2.0". Randomised lengths and counts
    // must never perturb the canonical neutral genome away from exactly 2.0.
    [Test]
    [CancelAfter(20_000)]
    public void EstimatePloidy_FlatDiploid_AlwaysExactlyTwo()
    {
        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);
            int n = 1 + rng.Next(50);
            var segments = new List<AlleleSpecificSegment>(n);
            long start = 0;
            for (int i = 0; i < n; i++)
            {
                long len = 1 + rng.Next(5_000_000); // arbitrary positive length
                segments.Add(Seg(start, start + len, 1, 1)); // 1:1 ⇒ total CN 2
                start += len;
            }

            double ploidy = EstimatePloidy(segments);

            ploidy.Should().Be(2.0, "all 1:1 segments ⇒ ψ = 2.0 exactly (INV-02), seed {0}", seed);
            AssertWellFormedPloidy(ploidy);
        }
    }

    // ── BE: a single 1:1 segment is the minimal flat-diploid genome ──────────
    [Test]
    public void EstimatePloidy_SingleDiploidSegment_IsExactlyTwo()
    {
        EstimatePloidy(new[] { Seg(0, 1, 1, 1) }).Should().Be(2.0);
    }

    // ── BE: a uniform genome at ANY single total-CN equals that CN exactly ────
    // INV-03 degenerate case (min = max = CN ⇒ ψ = CN). Sweeps CN 0..6.
    [Test]
    public void EstimatePloidy_UniformCopyNumber_EqualsThatCopyNumber()
    {
        for (int total = 0; total <= 6; total++)
        {
            int major = (total + 1) / 2;
            int minor = total - major; // major ≥ minor ≥ 0, major+minor = total
            var segments = new[]
            {
                Seg(0, 3 * Mb, major, minor),
                Seg(0, 7 * Mb, major, minor),
            };

            double ploidy = EstimatePloidy(segments);

            ploidy.Should().BeApproximately(total, 1e-12,
                "a uniform genome at total CN {0} has ψ = {0} (min = max)", total);
            AssertWellFormedPloidy(ploidy);
        }
    }

    #endregion

    #region ONCO-PLOIDY-001 — BE: fully amplified (huge length × CN, no overflow)

    // ── BE: fully amplified at a high CN over a GENOME-SCALE total length ─────
    // Every segment total CN = HighCn over ~3 Gb; Σ CNᵢ·Lᵢ in *integer* width
    // would be HighCn·3e9 — fine in long here, but the implementation accumulates
    // the weighted sum in DOUBLE, so even a pathologically large total length ×
    // CN cannot wrap to a negative/garbage ψ. ψ must equal HighCn exactly.
    [Test]
    [CancelAfter(20_000)]
    public void EstimatePloidy_FullyAmplified_GenomeScale_NoOverflow_EqualsCn()
    {
        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(seed);
            int major = 50 + rng.Next(200); // high major-allele CN
            int minor = rng.Next(200);
            int total = major + minor;

            int n = 1 + rng.Next(40);
            var segments = new List<AlleleSpecificSegment>(n);
            long start = 0;
            for (int i = 0; i < n; i++)
            {
                long len = 10 * Mb + rng.Next(90_000_000); // tens of Mb each
                segments.Add(Seg(start, start + len, major, minor));
                start += len;
            }

            double ploidy = EstimatePloidy(segments);

            ploidy.Should().BeApproximately(total, 1e-6,
                "fully amplified uniform genome ⇒ ψ = total CN {0} (seed {1})", total, seed);
            AssertWellFormedPloidy(ploidy);
        }
    }

    // ── BE: EXTREME single segment — near-Int64 length × large CN ─────────────
    // A single segment of length ~Int64.MaxValue at total CN = 100. Σ CNᵢ·Lᵢ in
    // long would be 100 · (~9.2e18) ≈ 9.2e20 — FAR beyond Int64.MaxValue (≈9.2e18),
    // which is exactly the overflow the double accumulator must avoid. ψ must be
    // the finite total CN (100), not a wrapped-negative or NaN. No DivideByZero.
    [Test]
    public void EstimatePloidy_ExtremeLengthTimesCn_NoIntegerWrap()
    {
        const int major = 60;
        const int minor = 40; // total CN 100
        long hugeLength = long.MaxValue - 1; // single segment spanning almost the whole long range
        var segments = new[] { Seg(0, hugeLength, major, minor) };

        double ploidy = EstimatePloidy(segments);

        ploidy.Should().BeApproximately(100.0, 1e-6,
            "uniform genome at total CN 100 ⇒ ψ = 100 even when CN·L overruns Int64 (must accumulate wider than long)");
        AssertWellFormedPloidy(ploidy);
    }

    // ── BE: ψ is the length-weighted mean — sum of products overflows long ────
    // Two genome-scale segments at different CNs whose Σ CNᵢ·Lᵢ would overflow a
    // long accumulator; the double-weighted mean must still land between min and
    // max CN (INV-03) and match the hand-computed weighted value.
    [Test]
    public void EstimatePloidy_TwoHugeSegments_WeightedMeanWithinRange()
    {
        long lenLow = 4_000_000_000L;  // 4 Gb
        long lenHigh = 6_000_000_000L; // 6 Gb
        var segments = new[]
        {
            Seg(0, lenLow, 1, 1),  // total CN 2
            Seg(0, lenHigh, 50, 50), // total CN 100
        };

        double ploidy = EstimatePloidy(segments);

        // (2·4e9 + 100·6e9) / (4e9 + 6e9) = (8e9 + 600e9) / 10e9 = 60.8
        ploidy.Should().BeApproximately(60.8, 1e-6, "hand-computed length-weighted mean");
        ploidy.Should().BeInRange(2.0, 100.0, "weighted mean lies within [min CN, max CN] (INV-03)");
        AssertWellFormedPloidy(ploidy);
    }

    #endregion

    #region ONCO-PLOIDY-001 — BE: empty segments and degenerate input (documented guards)

    // ── BE: empty segment set ⇒ ArgumentException, NEVER DivideByZero / NaN ───
    // §3.3 / §6.1 "Empty segment set ⇒ ArgumentException (Σ L = 0 undefined)".
    [Test]
    public void EstimatePloidy_EmptySegments_ThrowsArgumentException()
    {
        Action act = () => EstimatePloidy(Array.Empty<AlleleSpecificSegment>());

        act.Should().Throw<ArgumentException>("the length-weighted mean is undefined for Σ L = 0 (§6.1)");
    }

    [Test]
    public void EstimatePloidy_EmptyEnumerable_ThrowsArgumentException()
    {
        Action act = () => EstimatePloidy(Enumerable.Empty<AlleleSpecificSegment>());

        act.Should().Throw<ArgumentException>();
    }

    // ── BE: null input ⇒ ArgumentNullException (§3.3) ─────────────────────────
    [Test]
    public void EstimatePloidy_NullSegments_ThrowsArgumentNullException()
    {
        Action act = () => EstimatePloidy(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    // ── BE: a zero/negative-length segment ⇒ ArgumentException (§6.1) ─────────
    // End ≤ Start is invalid; it must be rejected, not silently contribute 0/neg
    // length that could zero or wrap the denominator.
    [Test]
    public void EstimatePloidy_ZeroLengthSegment_ThrowsArgumentException()
    {
        Action act = () => EstimatePloidy(new[] { Seg(100, 100, 1, 1) });

        act.Should().Throw<ArgumentException>("End ≤ Start is a non-positive length (§6.1)");
    }

    [Test]
    public void EstimatePloidy_NegativeLengthSegment_ThrowsArgumentException()
    {
        Action act = () => EstimatePloidy(new[] { Seg(200, 100, 1, 1) });

        act.Should().Throw<ArgumentException>();
    }

    // ── BE: a negative copy number ⇒ ArgumentException (§6.1) ─────────────────
    [Test]
    public void EstimatePloidy_NegativeCopyNumber_ThrowsArgumentException()
    {
        Action actMajor = () => EstimatePloidy(new[] { Seg(0, Mb, -1, 0) });
        Action actMinor = () => EstimatePloidy(new[] { Seg(0, Mb, 1, -1) });

        actMajor.Should().Throw<ArgumentException>("negative major CN is invalid input (§6.1)");
        actMinor.Should().Throw<ArgumentException>("negative minor CN is invalid input (§6.1)");
    }

    // ── BE: a genome of pure deletions (total CN 0) ⇒ ψ = 0, finite, no NaN ───
    // Σ CNᵢ·Lᵢ = 0 with positive Σ Lᵢ ⇒ 0/L = 0.0, a well-formed result (NOT a
    // 0/0 NaN). Guards the numerator-zero / denominator-positive corner.
    [Test]
    public void EstimatePloidy_AllZeroCopyNumber_IsExactlyZero_NotNaN()
    {
        var segments = new[]
        {
            Seg(0, 3 * Mb, 0, 0),
            Seg(0, 7 * Mb, 0, 0),
        };

        double ploidy = EstimatePloidy(segments);

        ploidy.Should().Be(0.0, "Σ(0·L) / Σ L = 0 (finite), not 0/0 = NaN");
        AssertWellFormedPloidy(ploidy);
    }

    #endregion

    #region ONCO-PLOIDY-001 — BE/RB: invariants over randomised mixed genomes

    // ── INV-03 over random mixed genomes: min CN ≤ ψ ≤ max CN, always finite ──
    // The bounding invariant is the strongest order-independent check: whatever
    // the random mix of CNs and lengths, the length-weighted mean must lie within
    // the observed copy-number range and never produce a malformed value.
    [Test]
    [CancelAfter(20_000)]
    public void EstimatePloidy_RandomMixedGenome_WithinCopyNumberRange()
    {
        for (int seed = 0; seed < 1_000; seed++)
        {
            var rng = new Random(seed);
            int n = 1 + rng.Next(30);
            var segments = new List<AlleleSpecificSegment>(n);
            int minTotal = int.MaxValue;
            int maxTotal = int.MinValue;
            long start = 0;
            for (int i = 0; i < n; i++)
            {
                int major = rng.Next(0, 12);
                int minor = rng.Next(0, major + 1); // minor ≤ major, both ≥ 0
                int total = major + minor;
                minTotal = Math.Min(minTotal, total);
                maxTotal = Math.Max(maxTotal, total);
                long len = 1 + rng.Next(10_000_000);
                segments.Add(Seg(start, start + len, major, minor));
                start += len;
            }

            double ploidy = EstimatePloidy(segments);

            AssertWellFormedPloidy(ploidy);
            ploidy.Should().BeInRange(minTotal, maxTotal,
                "length-weighted mean lies within [min CN, max CN] (INV-03), seed {0}", seed);
        }
    }

    // ── Order independence: ψ does not depend on the enumeration order ────────
    // A length-weighted mean is a symmetric aggregation; shuffling the segments
    // must yield the identical value (no accumulation-order divergence beyond fp
    // noise).
    [Test]
    [CancelAfter(20_000)]
    public void EstimatePloidy_OrderIndependent_ShuffleYieldsSameValue()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = 2 + rng.Next(20);
            var segments = new List<AlleleSpecificSegment>(n);
            long start = 0;
            for (int i = 0; i < n; i++)
            {
                int major = rng.Next(0, 8);
                int minor = rng.Next(0, major + 1);
                long len = 1 + rng.Next(5_000_000);
                segments.Add(Seg(start, start + len, major, minor));
                start += len;
            }

            double original = EstimatePloidy(segments);

            var shuffled = segments.OrderBy(_ => rng.Next()).ToList();
            double afterShuffle = EstimatePloidy(shuffled);

            afterShuffle.Should().BeApproximately(original, 1e-9,
                "ψ is a symmetric length-weighted mean (order-independent), seed {0}", seed);
        }
    }

    #endregion
}
