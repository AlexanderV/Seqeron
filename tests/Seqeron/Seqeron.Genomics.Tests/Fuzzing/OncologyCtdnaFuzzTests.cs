using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology circulating-tumour-DNA (ctDNA) tumour-fraction area — ONCO-CTDNA-001.
/// The unit under test is ctDNA tumour-fraction estimation / detection from a liquid biopsy,
/// implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs as the
/// public entry points of the <c>ctDNA analysis (ONCO-CTDNA-001)</c> region:
///   • <see cref="OncologyAnalyzer.CalculateTumorFraction(IEnumerable{VariantObservation})"/> —
///       TF = 2 · (mean plasma VAF) over clonal heterozygous copy-neutral diploid SNVs, clamped to [0, 1];
///   • <see cref="OncologyAnalyzer.CalculateMeanVaf(IEnumerable{VariantObservation})"/> —
///       arithmetic mean of per-reporter VAF = alt / total;
///   • <see cref="OncologyAnalyzer.CtDnaDetectionProbability(int, double, int)"/> —
///       Poisson detection probability p = 1 − e^(−n·d·k);
///   • <see cref="OncologyAnalyzer.ExpectedMutantMolecules(int, double, int)"/> — λ = n·d·k;
///   • <see cref="OncologyAnalyzer.IsCtDnaDetected(int, double, int, double)"/> —
///       detected ⇔ λ ≥ 1 AND p ≥ threshold (default 0.95).
/// (Minimal-residual-disease entry points — DetectMRD / TrackVariantsOverTime — are ONCO-MRD-001,
///  row 112, and are intentionally OUT OF SCOPE here.)
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts that the code NEVER
/// fails in an undisciplined way: no hang, no nonsense output, and no *unhandled* runtime exception
/// (DivideByZero / Overflow / NaN). Every input must resolve to EITHER a well-defined, theory-correct
/// tumour fraction / detection call OR a *documented, intentional* outcome (an
/// <see cref="ArgumentNullException"/> for a null collection, <see cref="ArgumentException"/> for an
/// empty collection, <see cref="ArgumentOutOfRangeException"/> for a per-variant VAF > 0.5 or a malformed
/// detection parameter). For ctDNA tumour fraction the headline hazards are:
///   • a DivideByZero on a ZERO-total-reads reporter (depth 0) — the shared VAF helper defines 0/0 ≡ 0,
///     so a tumour fraction is still produced, never an exception (§3.3 read-count delegation);
///   • a tumour fraction that ESCAPES [0, 1] — INV-04 says TF = 2 · mean VAF clamped to [0, 1]; a >1
///     leak must never reach the caller;
///   • a NaN tumour fraction on ULTRA-LOW depth (depth 1, a handful of reads) — mean of finite per-reporter
///     VAFs is finite, never NaN;
///   • a 0/0 mean on an EMPTY collection silently becoming NaN instead of the documented throw (§6.1).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-CTDNA-001 — ctDNA tumour-fraction estimation (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 111.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 111): "zero tumor reads, 100% tumor, ultra-low depth".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Mapping of the BE targets onto the documented contract:
///   • "zero tumor reads" ⇒ no mutant reads at any reporter ⇒ every VAF = 0 ⇒ TF = 2·0 = 0 (INV-04),
///       and the Poisson model gives λ = n·d·k = 0 ⇒ p = 0 ⇒ NOT detected (INV-02). No DivideByZero on
///       reporters whose total reads are also 0 (0/0 ≡ 0 in the VAF helper).
///   • "100% tumor" ⇒ all plasma reads are mutant. The documented diploid model caps a clonal het SNV
///       at VAF = 0.5 (TF = 2·0.5 = 1.0, the upper boundary); a per-reporter VAF > 0.5 (alt > total/2)
///       is a documented throw (§6.1), never a TF > 1 leak. The Poisson model saturates p → 1 (clamped).
///   • "ultra-low depth" ⇒ depth 1 or a handful of total reads. The estimate is still finite and in
///       [0, 1]; the Poisson decision is below the limit of detection when λ < 1 (§2.2 low-burden regime),
///       with no DivideByZero on tiny denominators and no nonsense fraction.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// CtDNA_Analysis.md (docs/algorithms/Oncology/CtDNA_Analysis.md):
///   • Tumour fraction = 2 · (mean clonal het VAF), clamped to [0, 1] (§2.2, INV-04). For a clonal het
///       copy-neutral diploid SNV the expected VAF is v = π/2, so TF = 2·v [4].
///   • Mean VAF = arithmetic mean of per-reporter alt/total (§2.2 "Mean VAF").
///   • Poisson detection: p = 1 − e^(−n·d·k) ∈ [0, 1] (INV-01); p = 0 ⇔ λ = 0 ⇔ n = 0 or d = 0 (INV-02);
///       p is non-decreasing in n, d, k (INV-03).
///   • Detection decision: detected ⇔ λ ≥ 1 AND p ≥ threshold (default 0.95) (§3 step 3).
///   • genomeEquivalents < 0, mutantAlleleFraction ∉ [0,1] or NaN, reporterCount < 1, and
///       minDetectionProbability ∉ (0,1] ⇒ ArgumentOutOfRangeException (§3.3).
///   • Null variant collection ⇒ ArgumentNullException; empty ⇒ ArgumentException (statistic undefined)
///       (§3.3, §6.1).
///   • A per-variant VAF > 0.5 in CalculateTumorFraction ⇒ ArgumentOutOfRangeException (impossible for a
///       diploid het SNV) (§3.3, §6.1). Read-count validation (alt ≤ total, non-negative) is delegated to
///       the shared VAF helper; total = 0 ⇒ VAF = 0 (no DivideByZero) (§9).
///   • Worked example (§6.2): two clonal het SNVs at VAF 0.10 and 0.20 ⇒ TF = 2·0.15 = 0.30.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyCtdnaFuzzTests
{
    // ── Well-formed-tumour-fraction assertion helper ─────────────────────────
    // Pins the documented numeric contract on EVERY accepted tumour fraction:
    // finite, never NaN/Infinity (no DivideByZero leak on depth-0 reporters,
    // no 0/0 mean) and inside [0, 1] (INV-04). This is what stops a fuzz test
    // from rubber-stamping a NaN or an out-of-range fraction such as 1.4.
    private static void AssertWellFormedTumorFraction(double tumorFraction)
    {
        double.IsNaN(tumorFraction).Should().BeFalse(
            "tumour fraction must never be NaN (no 0/0 in the mean VAF, no DivideByZero on depth-0 reporters)");
        double.IsInfinity(tumorFraction).Should().BeFalse("tumour fraction must be finite");
        tumorFraction.Should().BeInRange(0.0, 1.0, "tumour fraction is a fraction ⇒ TF ∈ [0, 1] (INV-04)");
    }

    // A plasma ctDNA reporter (clonal het copy-neutral diploid SNV) expressed as
    // a read-count VariantObservation: alt mutant reads out of `total` covering
    // reads. NormalAltReads/NormalTotalReads are unused by the tumour-fraction
    // path (it reads only TumorAltReads/TumorTotalReads).
    private static VariantObservation Reporter(int alt, int total)
        => new(
            Chromosome: "chr1",
            Position: 1,
            ReferenceAllele: "A",
            AlternateAllele: "T",
            TumorAltReads: alt,
            TumorTotalReads: total,
            NormalAltReads: 0,
            NormalTotalReads: total);

    #region ONCO-CTDNA-001 — Positive sanity (documented formula / detection call on hand-built examples)

    [Test]
    public void CalculateTumorFraction_DocWorkedExample_TwiceMeanVaf()
    {
        // Docs §6.2: two clonal het SNVs at VAF 0.10 (100/1000) and 0.20 (200/1000)
        // ⇒ mean VAF 0.15 ⇒ TF = 2 · 0.15 = 0.30. Pins the headline closed form.
        var variants = new[]
        {
            Reporter(alt: 100, total: 1000), // VAF 0.10
            Reporter(alt: 200, total: 1000), // VAF 0.20
        };

        double tf = CalculateTumorFraction(variants);

        AssertWellFormedTumorFraction(tf);
        tf.Should().BeApproximately(0.30, 1e-12);
        CalculateMeanVaf(variants).Should().BeApproximately(0.15, 1e-12);
    }

    [Test]
    public void CtDnaDetection_DocPoissonExample_DetectedAndProbabilityExact()
    {
        // Docs §6.2 / example: n = 15000 GE, d = 0.001 (0.1% VAF), k = 1 ⇒ λ = 15
        // ⇒ p = 1 − e^(−15) ≈ 0.99999969 ≥ 0.95 ⇒ detected.
        double p = CtDnaDetectionProbability(15_000, 0.001, 1);
        p.Should().BeApproximately(1.0 - Math.Exp(-15.0), 1e-15);
        p.Should().BeInRange(0.0, 1.0);

        ExpectedMutantMolecules(15_000, 0.001, 1).Should().BeApproximately(15.0, 1e-9);
        IsCtDnaDetected(15_000, 0.001, 1).Should().BeTrue("λ = 15 ≥ 1 and p ≈ 1 ≥ 0.95");
    }

    [Test]
    public void CalculateTumorFraction_UpperBoundaryVafHalf_TumorFractionIsOne()
    {
        // The documented upper boundary that PRODUCES TF = 1 is VAF = 0.5
        // (TF = 2 · 0.5 = 1.0): a fully clonal sample at a het diploid locus.
        double tf = CalculateTumorFraction(new[] { Reporter(alt: 500, total: 1000) });

        AssertWellFormedTumorFraction(tf);
        tf.Should().BeApproximately(1.0, 1e-12);
    }

    #endregion

    #region ONCO-CTDNA-001 — BE: zero tumor reads (no mutant reads ⇒ TF 0 / below LoD, no DivideByZero)

    [Test]
    [CancelAfter(10_000)]
    public void CalculateTumorFraction_AllReportersZeroAltReads_TumorFractionIsZero()
    {
        // "zero tumor reads": every reporter has 0 mutant reads out of real depth
        // ⇒ every VAF = 0 ⇒ TF = 2 · 0 = 0. No DivideByZero, no NaN (INV-04).
        var rng = new Random(111_001);
        for (int i = 0; i < 200; i++)
        {
            int n = rng.Next(1, 8);
            var variants = Enumerable.Range(0, n)
                .Select(_ => Reporter(alt: 0, total: rng.Next(1, 5000)))
                .ToArray();

            double tf = CalculateTumorFraction(variants);

            AssertWellFormedTumorFraction(tf);
            tf.Should().Be(0.0, "no mutant reads ⇒ mean VAF 0 ⇒ TF 0");
        }
    }

    [Test]
    [CancelAfter(10_000)]
    public void CalculateTumorFraction_ZeroAltAndZeroTotalReads_NoDivideByZero_TumorFractionZero()
    {
        // BE corner: reporters with BOTH zero alt AND zero total reads (depth-0
        // loci). The shared VAF helper defines 0/0 ≡ 0, so no DivideByZero — the
        // estimate is still a finite 0, not an exception or NaN (§9).
        var variants = new[]
        {
            Reporter(alt: 0, total: 0),
            Reporter(alt: 0, total: 0),
            Reporter(alt: 0, total: 0),
        };

        double tf = CalculateTumorFraction(variants);

        AssertWellFormedTumorFraction(tf);
        tf.Should().Be(0.0);
        CalculateMeanVaf(variants).Should().Be(0.0);
    }

    [Test]
    public void CtDnaDetection_ZeroEvidence_ProbabilityZeroAndNotDetected()
    {
        // λ = n·d·k = 0 when n = 0 (no input molecules) OR d = 0 (no mutant
        // fraction) ⇒ p = 1 − e⁰ = 0 ⇒ never detected (INV-02).
        CtDnaDetectionProbability(0, 0.001, 1).Should().Be(0.0);
        CtDnaDetectionProbability(15_000, 0.0, 1).Should().Be(0.0);
        ExpectedMutantMolecules(0, 0.5, 4).Should().Be(0.0);

        IsCtDnaDetected(0, 0.001, 1).Should().BeFalse();
        IsCtDnaDetected(15_000, 0.0, 1).Should().BeFalse();
    }

    #endregion

    #region ONCO-CTDNA-001 — BE: 100% tumor (upper bound, clamp ≤ 1, VAF>0.5 throws, no >1 leak)

    [Test]
    [CancelAfter(10_000)]
    public void CalculateTumorFraction_MaxValidEvidence_ClampedAtOne_NeverExceedsOne()
    {
        // "100% tumor": maximal per-reporter evidence the diploid model accepts is
        // VAF = 0.5 (alt = total/2). Many such reporters ⇒ TF saturates at 1.0,
        // clamped; it must never leak > 1 (INV-04).
        var rng = new Random(111_002);
        for (int i = 0; i < 200; i++)
        {
            int n = rng.Next(1, 8);
            var variants = Enumerable.Range(0, n)
                .Select(_ =>
                {
                    int total = 2 * rng.Next(1, 2500); // even depth so total/2 is exact
                    return Reporter(alt: total / 2, total: total); // VAF exactly 0.5
                })
                .ToArray();

            double tf = CalculateTumorFraction(variants);

            AssertWellFormedTumorFraction(tf);
            tf.Should().BeApproximately(1.0, 1e-12, "VAF 0.5 ⇒ TF 2·0.5 = 1.0 (upper boundary)");
        }
    }

    [Test]
    [CancelAfter(10_000)]
    public void CalculateTumorFraction_AltExceedsHalf_ThrowsNotLeakingFractionAboveOne()
    {
        // A reporter where every read (or > half) is mutant violates the clonal het
        // diploid model (VAF > 0.5). It is a DOCUMENTED throw (§6.1), NOT a TF > 1.
        var rng = new Random(111_003);
        for (int i = 0; i < 200; i++)
        {
            int total = rng.Next(2, 5000);
            int alt = rng.Next(total / 2 + 1, total + 1); // strictly > total/2 ⇒ VAF > 0.5

            Action act = () => CalculateTumorFraction(new[] { Reporter(alt, total) });

            act.Should().Throw<ArgumentOutOfRangeException>(
                "VAF > 0.5 is impossible for a clonal het diploid SNV; must throw, never leak TF > 1");
        }
    }

    [Test]
    public void CalculateTumorFraction_AllReadsMutant_DepthOne_Throws()
    {
        // Pathological "100% tumor at depth 1": alt = total = 1 ⇒ VAF = 1.0 > 0.5
        // ⇒ documented throw, not a TF of 2.0.
        Action act = () => CalculateTumorFraction(new[] { Reporter(alt: 1, total: 1) });
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void CtDnaDetection_SaturatingEvidence_ProbabilityClampedAtOne()
    {
        // Very large λ = n·d·k ⇒ p → 1 but never exceeds 1 (INV-01); detected.
        double p = CtDnaDetectionProbability(int.MaxValue, 1.0, 1000);
        p.Should().BeInRange(0.0, 1.0);
        p.Should().BeApproximately(1.0, 1e-12);
        IsCtDnaDetected(int.MaxValue, 1.0, 1000).Should().BeTrue();
    }

    #endregion

    #region ONCO-CTDNA-001 — BE: ultra-low depth (depth 1 / handful of reads ⇒ finite TF, below-LoD call)

    [Test]
    [CancelAfter(10_000)]
    public void CalculateTumorFraction_UltraLowDepth_FiniteAndInRange()
    {
        // "ultra-low depth": each reporter has 1..3 total reads. With alt ∈ {0,1}
        // and alt ≤ total/2 enforced, VAF ∈ {0, 0.5}. The estimate must be finite
        // and in [0, 1] regardless — no NaN, no DivideByZero on tiny denominators.
        var rng = new Random(111_004);
        for (int i = 0; i < 300; i++)
        {
            int n = rng.Next(1, 6);
            var variants = Enumerable.Range(0, n)
                .Select(_ =>
                {
                    int total = rng.Next(1, 4); // 1..3 reads
                    int alt = rng.Next(0, total / 2 + 1); // keep VAF ≤ 0.5
                    return Reporter(alt, total);
                })
                .ToArray();

            double tf = CalculateTumorFraction(variants);
            AssertWellFormedTumorFraction(tf);
        }
    }

    [Test]
    public void CalculateTumorFraction_SingleReporterDepthOneZeroAlt_TumorFractionZero()
    {
        // The smallest non-degenerate reporter: depth 1, no mutant read ⇒ VAF 0
        // ⇒ TF 0. Median/mean-of-one must not crash on the tiny denominator.
        double tf = CalculateTumorFraction(new[] { Reporter(alt: 0, total: 1) });
        AssertWellFormedTumorFraction(tf);
        tf.Should().Be(0.0);
    }

    [Test]
    public void CtDnaDetection_UltraLowGenomeEquivalents_BelowLimitOfDetection()
    {
        // Ultra-low input (a single genome equivalent at a sub-LoD allele fraction)
        // gives λ = n·d·k < 1 ⇒ NOT detected (the λ ≥ 1 physical floor), even though
        // the probability is a finite positive number (§3 step 3, low-burden regime).
        double lambda = ExpectedMutantMolecules(1, 0.0001, 1);
        lambda.Should().BeLessThan(1.0);

        double p = CtDnaDetectionProbability(1, 0.0001, 1);
        p.Should().BeInRange(0.0, 1.0);
        double.IsNaN(p).Should().BeFalse();

        IsCtDnaDetected(1, 0.0001, 1).Should().BeFalse("λ < 1 ⇒ below the limit of detection");
    }

    #endregion

    #region ONCO-CTDNA-001 — BE: empty / null collections (documented throws, no 0/0 mean)

    [Test]
    public void CalculateTumorFraction_Null_ThrowsArgumentNullException()
    {
        Action act = () => CalculateTumorFraction(null!);
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void CalculateTumorFraction_Empty_ThrowsArgumentException_NotNaN()
    {
        // Empty reporter set ⇒ 0/0 mean is UNDEFINED ⇒ documented ArgumentException,
        // never a silently leaked NaN tumour fraction (§6.1).
        Action act = () => CalculateTumorFraction(Array.Empty<VariantObservation>());
        act.Should().Throw<ArgumentException>().Which.Should().NotBeOfType<ArgumentNullException>();
    }

    [Test]
    public void CalculateMeanVaf_NullAndEmpty_Throw()
    {
        ((Action)(() => CalculateMeanVaf(null!))).Should().Throw<ArgumentNullException>();
        ((Action)(() => CalculateMeanVaf(Array.Empty<VariantObservation>())))
            .Should().Throw<ArgumentException>();
    }

    #endregion

    #region ONCO-CTDNA-001 — BE: malformed detection scalars (NaN / negative / out-of-range ⇒ throws)

    [Test]
    public void CtDnaDetectionProbability_MalformedScalars_Throw()
    {
        // §3.3 domain limits — each malformed scalar is a documented throw, not a NaN p.
        ((Action)(() => CtDnaDetectionProbability(-1, 0.01, 1)))
            .Should().Throw<ArgumentOutOfRangeException>("genomeEquivalents < 0");
        ((Action)(() => CtDnaDetectionProbability(100, double.NaN, 1)))
            .Should().Throw<ArgumentOutOfRangeException>("mutantAlleleFraction NaN");
        ((Action)(() => CtDnaDetectionProbability(100, -0.01, 1)))
            .Should().Throw<ArgumentOutOfRangeException>("mutantAlleleFraction < 0");
        ((Action)(() => CtDnaDetectionProbability(100, 1.01, 1)))
            .Should().Throw<ArgumentOutOfRangeException>("mutantAlleleFraction > 1");
        ((Action)(() => CtDnaDetectionProbability(100, 0.01, 0)))
            .Should().Throw<ArgumentOutOfRangeException>("reporterCount < 1");
    }

    [Test]
    public void IsCtDnaDetected_MalformedThreshold_Throws()
    {
        // minDetectionProbability must be in (0, 1]; 0, > 1, and NaN are throws.
        ((Action)(() => IsCtDnaDetected(15_000, 0.001, 1, 0.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => IsCtDnaDetected(15_000, 0.001, 1, 1.01)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => IsCtDnaDetected(15_000, 0.001, 1, double.NaN)))
            .Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    [CancelAfter(10_000)]
    public void CalculateTumorFraction_NegativeOrAltAboveTotal_DelegatedReadCountThrows()
    {
        // Read-count validation is delegated to the shared VAF helper: negative
        // counts and alt > total are documented throws (§3.3), never a NaN/negative TF.
        var rng = new Random(111_005);
        for (int i = 0; i < 100; i++)
        {
            int total = rng.Next(1, 1000);
            int altAboveTotal = total + rng.Next(1, 100);
            ((Action)(() => CalculateTumorFraction(new[] { Reporter(altAboveTotal, total) })))
                .Should().Throw<ArgumentOutOfRangeException>("alt > total is invalid");

            ((Action)(() => CalculateTumorFraction(new[] { Reporter(-rng.Next(1, 100), total) })))
                .Should().Throw<ArgumentOutOfRangeException>("negative alt is invalid");
        }
    }

    #endregion

    #region ONCO-CTDNA-001 — Invariant sweep (TF ∈ [0,1] = 2·mean VAF; p ∈ [0,1] monotone in λ)

    [Test]
    [CancelAfter(20_000)]
    public void CalculateTumorFraction_RandomValidReporters_AlwaysTwiceMeanVafClampedToUnitInterval()
    {
        // INV-04 sweep: for any set of valid clonal het reporters (each VAF ≤ 0.5),
        // TF = min(2 · mean VAF, 1) and stays in [0, 1]. Never NaN / > 1 / negative.
        var rng = new Random(111_006);
        for (int i = 0; i < 500; i++)
        {
            int n = rng.Next(1, 12);
            var variants = Enumerable.Range(0, n)
                .Select(_ =>
                {
                    int total = rng.Next(0, 6000);
                    int alt = total == 0 ? 0 : rng.Next(0, total / 2 + 1); // VAF ≤ 0.5
                    return Reporter(alt, total);
                })
                .ToArray();

            double tf = CalculateTumorFraction(variants);
            double meanVaf = CalculateMeanVaf(variants);

            AssertWellFormedTumorFraction(tf);
            tf.Should().BeApproximately(Math.Min(2.0 * meanVaf, 1.0), 1e-12);
        }
    }

    [Test]
    [CancelAfter(20_000)]
    public void CtDnaDetectionProbability_RandomDomain_InUnitIntervalAndMonotoneInLambda()
    {
        // INV-01 / INV-03 sweep: p ∈ [0, 1], finite, and non-decreasing as λ = n·d·k
        // grows (here: increasing d at fixed n, k). Never NaN, never > 1.
        var rng = new Random(111_007);
        for (int i = 0; i < 500; i++)
        {
            int n = rng.Next(0, 200_000);
            int k = rng.Next(1, 64);
            double dLow = rng.NextDouble() * 0.5;
            double dHigh = dLow + rng.NextDouble() * (1.0 - dLow);

            double pLow = CtDnaDetectionProbability(n, dLow, k);
            double pHigh = CtDnaDetectionProbability(n, dHigh, k);

            foreach (double p in new[] { pLow, pHigh })
            {
                double.IsNaN(p).Should().BeFalse();
                p.Should().BeInRange(0.0, 1.0);
            }
            pHigh.Should().BeGreaterThanOrEqualTo(pLow - 1e-12, "p is non-decreasing in d (INV-03)");
        }
    }

    #endregion
}
