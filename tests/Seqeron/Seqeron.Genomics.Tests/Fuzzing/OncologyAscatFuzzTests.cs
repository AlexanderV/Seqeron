using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology allele-specific copy-number derivation area — ONCO-ASCAT-001.
/// The units under test are the ASCAT-style entry points implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
///   • <see cref="OncologyAnalyzer.SegmentAlleleSpecific(IEnumerable{OncologyAnalyzer.AlleleSpecificLocus},double,double,int)"/>
///       — greedy joint logR/mirrored-BAF mean-shift segmentation of per-locus signal;
///   • <see cref="OncologyAnalyzer.FitPurityPloidy(IReadOnlyList{OncologyAnalyzer.AlleleSpecificSegmentSummary},double,double,double,double,double,double,double)"/>
///       — the ASCAT grid fit → ρ, ψ, GoF % and the implied integer allele-specific segments;
///   • <see cref="OncologyAnalyzer.DeriveMultiplicity(double,double,int,int)"/>
///       — McGranahan rounded/clamped mutation multiplicity.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / out-of-domain inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang / infinite grid loop,
/// no nonsense out-of-contract output (a negative or non-integer copy number, a
/// minor > major segment, a GoF > 100 %, a purity escaping its grid, a NaN/Infinity
/// leaking through 2^(r/γ)), and no *unhandled* runtime fault (IndexOutOfRange on an
/// empty locus run, DivideByZero on ρ → 0, an overflow-wrapped length). Every input
/// must resolve to EITHER a well-defined, theory-correct result OR a *documented,
/// intentional* outcome (an <see cref="ArgumentNullException"/> for null loci/segments,
/// an <see cref="ArgumentException"/> for an empty segment set or a null chromosome
/// label, an <see cref="ArgumentOutOfRangeException"/> for an out-of-range threshold,
/// grid bound, or multiplicity argument).
///
/// For ASCAT derivation the headline BE hazards (checklist row 235, targets
/// "empty loci, single locus, all-het, all-hom, extreme logR/BAF") are:
///   • empty loci — SegmentAlleleSpecific over an empty sequence emits ZERO segments
///     (an empty list), NEVER an IndexOutOfRange from BuildSegmentSummary on an empty
///     run; that empty summary list then makes FitPurityPloidy throw the documented
///     ArgumentException (§3.3 "empty segments ⇒ ArgumentException"), not a 0/0 GoF;
///   • single locus — one locus ⇒ exactly one segment, LocusCount = 1, the minimal
///     genome (§6.1 "Single locus per chromosome ⇒ one segment, LocusCount=1");
///   • all-heterozygous (BAF ≈ 0/1 unbalanced) and all-homozygous/balanced (BAF = 0.5)
///     single-segment genomes — the original runASCAT returns rho = NA for all of them
///     (no strict 7×7 local minimum passes its filters; R run 2026-09), so the fit must end
///     in the documented "no optimum" outcome (TryFitPurityPloidy = false,
///     FitPurityPloidy → InvalidOperationException), never a malformed fit;
///   • extreme logR (±large) and extreme BAF (0, 0.5, 1) — 2^(r/γ) can underflow to 0
///     or overflow toward +Infinity; the emitted integer copy numbers must remain
///     finite, non-negative, sorted (major ≥ minor), and a reported ρ lies on its grid
///     (values above 1 reported as 1) — never a NaN/Infinity copy number or a wrapped state.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-ASCAT-001 — Allele-specific copy-number derivation (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 235.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 235): "empty loci, single locus, all-het, all-hom,
///     extreme logR/BAF".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Allele_Specific_Copy_Number_Derivation.md (docs/algorithms/Oncology/):
///   • ASCAT seg_raw rounding (negative allele folded, half-to-even) ⇒ every emitted
///       segment has major ≥ minor ≥ 0 integers (INV-04, §4.1).                  ── FitPurityPloidy
///   • GoF = (1 − d/TheoretMaxdist)·100 with 0 ≤ d ≤ TheoretMaxdist ⇒ GoF ≤ 100 %
///       (INV-03, §2.2).                                                          ── FitPurityPloidy
///   • Reported ρ on the purity grid, ≤ 1; ψ on seq(ploidyMin − 0.5, ploidyMax + 0.5)  (§2.2).
///   • No acceptable ASCAT optimum ⇒ TryFitPurityPloidy false / InvalidOperationException (§3.3).
///   • Every emitted AlleleSpecificSegment.Length > 0 (single-position ⇒ 1 bp)   (§4.2).
///   • Multiplicity m ∈ [1, majorCopyNumber] (INV-02, §2.2 clamp).               ── DeriveMultiplicity
///   • Single locus per chromosome ⇒ one segment, LocusCount = 1 (§6.1).         ── SegmentAlleleSpecific
///   • Single-segment genome ⇒ ASCAT rho = NA (§6.1).
///   • Null loci/segments ⇒ ArgumentNullException; empty segments ⇒
///       ArgumentException; out-of-range thresholds/grid/multiplicity args ⇒
///       ArgumentOutOfRangeException; null chromosome ⇒ ArgumentException (§3.3).
///
/// Tests encode the DOCUMENTED contract, derived independently from the algorithm
/// doc — a test that would still pass against a wrong implementation (e.g. one that
/// emitted minor > major, a negative CN, a GoF of 137 %, or a NaN copy number) would
/// FAIL here. They are not a rubber-stamp of the current return values.
///
/// SOURCE (2026-09, B24 F12): FitPurityPloidy became a runASCAT port; the single-segment
/// expectations were re-derived from the original R code (rho = NA), not bent to the
/// code. SegmentAlleleSpecific returns an empty list for empty loci, DeriveMultiplicity
/// clamps to [1, major], the (ρ, ψ) grid is bounded, so there is no hang.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyAscatFuzzTests
{
    private const double Tol = 1e-12;

    private static AlleleSpecificLocus Locus(string chrom, long pos, double logR, double baf) =>
        new(chrom, pos, logR, baf);

    // ── Well-formed-fit assertion helper ──────────────────────────────────────
    // Pins the documented numeric contract on EVERY fit, whatever the input:
    //   • ρ, ψ, GoF finite (no NaN/Infinity leaking from 2^(r/γ));
    //   • ρ within its grid bounds ⊆ (0, 1], ψ within its grid bounds;
    //   • GoF ≤ 100 % (INV-03);
    //   • every emitted segment: major ≥ minor ≥ 0 INTEGER copy numbers (INV-04)
    //     and a strictly positive length (§4.2).
    // This is what stops a fuzz test from rubber-stamping a malformed fit green.
    private static void AssertWellFormedFit(PurityPloidyFit fit, double purityMin, double purityMax,
        double ploidyMin, double ploidyMax, double ploidyStep = 0.05)
    {
        double.IsNaN(fit.Purity).Should().BeFalse("purity ρ must never be NaN");
        double.IsNaN(fit.Ploidy).Should().BeFalse("ploidy must never be NaN");
        double.IsNaN(fit.GoodnessOfFit).Should().BeFalse("GoF must never be NaN");
        double.IsInfinity(fit.Purity).Should().BeFalse();
        double.IsInfinity(fit.Ploidy).Should().BeFalse();
        double.IsInfinity(fit.GoodnessOfFit).Should().BeFalse("GoF must never be ±Infinity");

        fit.Purity.Should().BeInRange(Math.Min(purityMin, 1.0) - Tol, Math.Min(purityMax, 1.0) + Tol,
            "reported ρ lies on the purity grid, grid points above 1 reported as 1 (runASCAT: rho_opt1 > 1 ⇒ 1)");
        fit.Psi.Should().BeInRange(ploidyMin - 0.5 - Tol, ploidyMax + 0.5 + ploidyStep + Tol,
            "ψ lies on the ASCAT grid seq(min_ploidy − 0.5, max_ploidy + 0.5, step)");
        fit.GoodnessOfFit.Should().BeGreaterThan(80.0, "an accepted ASCAT optimum has GoF > MINGOODNESSOFFIT = 80");
        fit.GoodnessOfFit.Should().BeLessThanOrEqualTo(100.0 + 1e-9, "GoF ≤ 100 % (INV-03)");

        fit.Segments.Should().NotBeNull();
        foreach (AlleleSpecificSegment s in fit.Segments)
        {
            s.MinorCopyNumber.Should().BeGreaterThanOrEqualTo(0, "ASCAT folds a negative allele into the other (INV-04)");
            s.MajorCopyNumber.Should().BeGreaterThanOrEqualTo(s.MinorCopyNumber, "major ≥ minor in every segment (INV-04)");
            s.Length.Should().BeGreaterThan(0, "every emitted segment has a positive span (§4.2)");
        }
    }

    /// <summary>
    /// The disciplined outcomes of the ASCAT fit: either a well-formed optimum, or ASCAT's "no optimum" (rho = NA)
    /// reported as TryFitPurityPloidy = false / FitPurityPloidy → InvalidOperationException.
    /// </summary>
    private static bool AssertWellFormedOrNoOptimum(IReadOnlyList<AlleleSpecificSegmentSummary> segs,
        double purityStep = 0.01, double ploidyStep = 0.05)
    {
        bool found = TryFitPurityPloidy(segs, out PurityPloidyFit fit, purityStep: purityStep, ploidyStep: ploidyStep);
        if (found)
        {
            AssertWellFormedFit(fit, 0.1, 1.05, 1.5, 5.5, ploidyStep);
        }
        else
        {
            ((Action)(() => FitPurityPloidy(segs, purityStep: purityStep, ploidyStep: ploidyStep)))
                .Should().Throw<InvalidOperationException>("ASCAT's rho = NA outcome is a documented exception");
        }

        return found;
    }

    #region ONCO-ASCAT-001 — SegmentAlleleSpecific: positive sanity

    // ── POSITIVE sanity: a copy-neutral-LOH region splits from a balanced region ──
    // §4.1: BAF-aware segmentation must NOT merge a 2:0 LOH run (folded BAF = 1.0,
    // logR ≈ 0) with a balanced 1:1 run (folded BAF = 0.5, logR ≈ 0): identical logR,
    // different BAF ⇒ a BAF mean-shift split. Pins that BAF actually segments.
    [Test]
    public void SegmentAlleleSpecific_LohVsBalancedSameLogR_SplitsOnBaf()
    {
        var loci = new[]
        {
            Locus("1", 0, 0.0, 0.5), Locus("1", 100, 0.0, 0.5),  // balanced 1:1
            Locus("1", 200, 0.0, 1.0), Locus("1", 300, 0.0, 1.0), // copy-neutral LOH 2:0
        };

        var segments = SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, bafChangeThreshold: 0.1);

        segments.Should().HaveCountGreaterThan(1, "a BAF shift at identical logR must split the run (§4.1)");
        segments.Sum(s => s.LocusCount).Should().Be(loci.Length, "every locus is summarised exactly once");
    }

    #endregion

    #region ONCO-ASCAT-001 — BE: empty loci and single locus

    // ── BE: empty loci ⇒ ZERO segments, never an IndexOutOfRange on an empty run ──
    // BuildSegmentSummary is only called for current.Count > 0, so an empty sequence
    // yields an empty list (not a crash). The empty list is the documented input that
    // then makes FitPurityPloidy throw ArgumentException.
    [Test]
    public void SegmentAlleleSpecific_EmptyLoci_ReturnsEmptyListNoThrow()
    {
        var segments = SegmentAlleleSpecific(Array.Empty<AlleleSpecificLocus>(), logRChangeThreshold: 0.2);

        segments.Should().NotBeNull();
        segments.Should().BeEmpty("no loci ⇒ no segments (BuildSegmentSummary never indexes an empty run)");
    }

    [Test]
    public void SegmentAlleleSpecific_EmptyLoci_FeedsFitPurityPloidyDocumentedThrow()
    {
        var segments = SegmentAlleleSpecific(Enumerable.Empty<AlleleSpecificLocus>(), logRChangeThreshold: 0.2);

        Action act = () => FitPurityPloidy(segments);

        act.Should().Throw<ArgumentException>("empty segment set ⇒ documented ArgumentException, not a 0/0 GoF (§3.3)");
    }

    // ── BE: a single locus ⇒ exactly one segment, LocusCount = 1 (§6.1) ──────────
    [Test]
    public void SegmentAlleleSpecific_SingleLocus_OneSegmentLocusCountOne()
    {
        var segments = SegmentAlleleSpecific(new[] { Locus("1", 42, 0.3, 0.5) }, logRChangeThreshold: 0.2);

        segments.Should().HaveCount(1, "a single locus is one segment (§6.1)");
        segments[0].LocusCount.Should().Be(1);
        segments[0].Start.Should().Be(42);
        segments[0].End.Should().Be(42);
    }

    // ── BE: null loci ⇒ ArgumentNullException; null chromosome label ⇒ ArgumentException ──
    [Test]
    public void SegmentAlleleSpecific_NullLoci_ThrowsArgumentNullException()
    {
        Action act = () => SegmentAlleleSpecific(null!, logRChangeThreshold: 0.2);
        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void SegmentAlleleSpecific_NullChromosomeLabel_ThrowsArgumentException()
    {
        var loci = new[] { Locus(null!, 0, 0.0, 0.5) };
        Action act = () => SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2);
        act.Should().Throw<ArgumentException>("a null chromosome label is invalid input (§3.3)");
    }

    // ── BE: non-positive thresholds / minLociPerSegment < 1 ⇒ ArgumentOutOfRangeException ──
    [Test]
    public void SegmentAlleleSpecific_NonPositiveOrNaNThresholds_Throw()
    {
        var loci = new[] { Locus("1", 0, 0.0, 0.5) };
        ((Action)(() => SegmentAlleleSpecific(loci, logRChangeThreshold: 0.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => SegmentAlleleSpecific(loci, logRChangeThreshold: -1.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => SegmentAlleleSpecific(loci, logRChangeThreshold: double.NaN)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, bafChangeThreshold: 0.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, bafChangeThreshold: double.NaN)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 0)))
            .Should().Throw<ArgumentOutOfRangeException>();
    }

    // ── BE: segmentation never invents or drops loci, whatever the random signal ──
    // The structural invariant: Σ LocusCount == #loci and 1 ≤ #segments ≤ #loci, even
    // with extreme/NaN-free logR and out-of-[0,1] BAF, and no segment summary is NaN.
    [Test]
    [CancelAfter(20_000)]
    public void SegmentAlleleSpecific_RandomSignal_PreservesLocusCountNoMalformedSummary()
    {
        for (int seed = 0; seed < 500; seed++)
        {
            var rng = new Random(seed);
            int n = 1 + rng.Next(60);
            var loci = new List<AlleleSpecificLocus>(n);
            for (int i = 0; i < n; i++)
            {
                // extreme logR (±large) and BAF including out-of-domain values 0, 0.5, 1 and beyond.
                double logR = (rng.NextDouble() - 0.5) * 2_000.0;
                double baf = rng.Next(5) switch { 0 => 0.0, 1 => 0.5, 2 => 1.0, 3 => -0.3, _ => rng.NextDouble() };
                loci.Add(Locus(rng.Next(2) == 0 ? "1" : "2", i * 100L, logR, baf));
            }

            var segments = SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2 + rng.NextDouble());

            segments.Sum(s => s.LocusCount).Should().Be(n, "every locus summarised exactly once, seed {0}", seed);
            segments.Count.Should().BeInRange(1, n, "1 ≤ #segments ≤ #loci, seed {0}", seed);
            foreach (var s in segments)
            {
                double.IsNaN(s.MeanLogR).Should().BeFalse("a finite-logR run yields a finite mean, seed {0}", seed);
                double.IsNaN(s.MeanBAF).Should().BeFalse("folded-BAF mean is finite, seed {0}", seed);
            }
        }
    }

    #endregion

    #region ONCO-ASCAT-001 — BE: FitPurityPloidy all-het / all-hom / extreme logR-BAF

    // ── BE: empty / null segments ⇒ documented throws (§3.3) ─────────────────────
    [Test]
    public void FitPurityPloidy_EmptySegments_ThrowsArgumentException()
    {
        ((Action)(() => FitPurityPloidy(Array.Empty<AlleleSpecificSegmentSummary>())))
            .Should().Throw<ArgumentException>("at least one segment is required (§3.3)");
    }

    [Test]
    public void FitPurityPloidy_NullSegments_ThrowsArgumentNullException()
    {
        ((Action)(() => FitPurityPloidy(null!))).Should().Throw<ArgumentNullException>();
    }

    // ── BE: out-of-range grid bounds ⇒ ArgumentOutOfRangeException ────────────────
    [Test]
    public void FitPurityPloidy_OutOfRangeGrid_Throws()
    {
        var seg = new[] { new AlleleSpecificSegmentSummary("1", 0, 1000, 0.0, 0.5, 5) };
        ((Action)(() => FitPurityPloidy(seg, purityMin: 0.0)))
            .Should().Throw<ArgumentOutOfRangeException>("purityMin must be in (0,1]");
        ((Action)(() => FitPurityPloidy(seg, purityMax: double.PositiveInfinity)))
            .Should().Throw<ArgumentOutOfRangeException>("purityMax must be finite");
        ((Action)(() => FitPurityPloidy(seg, purityMin: 0.6, purityMax: 0.5)))
            .Should().Throw<ArgumentOutOfRangeException>("purityMax must be ≥ purityMin");
        ((Action)(() => FitPurityPloidy(seg, purityStep: 0.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => FitPurityPloidy(seg, ploidyMin: 0.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => FitPurityPloidy(seg, ploidyMin: 4.0, ploidyMax: 2.0)))
            .Should().Throw<ArgumentOutOfRangeException>("ploidyMax must be ≥ ploidyMin");
        ((Action)(() => FitPurityPloidy(seg, gamma: 0.0)))
            .Should().Throw<ArgumentOutOfRangeException>();
        ((Action)(() => FitPurityPloidy(seg, purityStep: 1e-9)))
            .Should().Throw<ArgumentOutOfRangeException>("the grid size is bounded");
    }

    // ── BE: a purity grid above 1 is ASCAT's default (max_purity = 1.05) and is accepted ──
    [Test]
    [CancelAfter(30_000)]
    public void FitPurityPloidy_PurityGridAboveOne_AcceptedAndReportedAsOne()
    {
        var segs = new[]
        {
            new AlleleSpecificSegmentSummary("X", 0, 1000, -0.19799990120548305, 0.7800586872805495, 55),
            new AlleleSpecificSegmentSummary("3", 0, 1000, 0.097918462573404322, 0.66269643393554578, 13),
            new AlleleSpecificSegmentSummary("1", 0, 1000, 0.065882791496869667, 0.6385382413628512, 15),
        };

        var fit = FitPurityPloidy(segs, purityMax: 1.5);

        fit.Purity.Should().BeLessThanOrEqualTo(1.0, "runASCAT reports rho_opt1 > 1 as 1");
    }

    // ── BE: single-segment genomes (balanced, all-het unbalanced, zero-span, extreme logR) ──
    // A single segment cannot produce a strict 7 × 7 local minimum that passes the ASCAT filters: the original
    // runASCAT returns rho = NA for every one of these inputs (R 4.3 run; Evidence). They must end in the
    // documented "no optimum" outcome — never NaN, a hang, or an undocumented exception.
    [Test]
    [CancelAfter(60_000)]
    public void FitPurityPloidy_SingleSegmentGenomes_NoAscatOptimumAsRunAscat()
    {
        var inputs = new List<AlleleSpecificSegmentSummary>
        {
            new("1", 0, 1_000_000, 0.0, 0.5, 100),   // balanced
            new("1", 0, 1_000_000, 0.0, 1.0, 100),   // all-het unbalanced (LOH)
            new("1", 500, 500, 0.0, 0.5, 1),         // zero span
        };
        foreach (double r in new[] { -400.0, -50.0, -10.0, 10.0, 50.0, 400.0 })
        {
            foreach (double b in new[] { 0.0, 0.5, 1.0 })
            {
                inputs.Add(new AlleleSpecificSegmentSummary("1", 0, 1_000_000, r, b, 50));
            }
        }

        foreach (AlleleSpecificSegmentSummary seg in inputs)
        {
            AssertWellFormedOrNoOptimum(new[] { seg }).Should().BeFalse(
                $"runASCAT returns rho = NA for the single segment (r = {seg.MeanLogR}, b = {seg.MeanBAF})");
        }
    }

    // ── BE: extreme logR with a fixed (ρ, ψ) (ASCAT rho_manual/psi_manual) stays well-formed ─
    // r = ±400 ⇒ 2^(r/γ) over/underflows towards 1e120 / 0; the integer segments must still be
    // non-negative with major ≥ minor, and GoF ≤ 100 %.
    [Test]
    [CancelAfter(30_000)]
    public void EvaluatePurityPloidy_ExtremeLogR_NoNaNOrMalformedCopyNumber()
    {
        foreach (double r in new[] { -400.0, -50.0, -10.0, 10.0, 50.0, 400.0 })
        {
            foreach (double b in new[] { 0.0, 0.5, 1.0 })
            {
                var seg = new[] { new AlleleSpecificSegmentSummary("1", 0, 1_000_000, r, b, 50) };

                var fit = EvaluatePurityPloidy(seg, 0.7, 2.5);

                double.IsNaN(fit.GoodnessOfFit).Should().BeFalse();
                fit.GoodnessOfFit.Should().BeLessThanOrEqualTo(100.0 + 1e-9);
                fit.Segments[0].MinorCopyNumber.Should().BeGreaterThanOrEqualTo(0);
                fit.Segments[0].MajorCopyNumber.Should().BeGreaterThanOrEqualTo(fit.Segments[0].MinorCopyNumber);
            }
        }
    }

    // ── BE: extreme-logR multi-segment genome over a randomised grid ─────────────
    // WHATEVER the random mix of extreme logR, boundary BAF (0/0.5/1) and probe counts, the fit
    // either returns a well-formed ASCAT optimum (ρ on the grid, GoF in (80, 100], major ≥ minor ≥ 0
    // integers, length > 0, one segment per summary) or ASCAT's documented "no optimum".
    [Test]
    [CancelAfter(60_000)]
    public void FitPurityPloidy_RandomExtremeGenome_AlwaysWellFormedOrNoOptimum()
    {
        for (int seed = 0; seed < 120; seed++)
        {
            var rng = new Random(seed);
            int n = 1 + rng.Next(8);
            var segs = new List<AlleleSpecificSegmentSummary>(n);
            long start = 0;
            for (int i = 0; i < n; i++)
            {
                double r = (rng.NextDouble() - 0.5) * 600.0;          // ±300 logR
                double b = rng.Next(4) switch { 0 => 0.0, 1 => 0.5, 2 => 1.0, _ => 0.5 + rng.NextDouble() * 0.5 };
                long len = 1 + rng.Next(5_000_000);
                segs.Add(new AlleleSpecificSegmentSummary("1", start, start + len, r, b, 1 + rng.Next(200)));
                start += len + 1;
            }

            // Coarsen the grid a little so the random sweep stays fast but still exercises the fit.
            if (AssertWellFormedOrNoOptimum(segs, purityStep: 0.05, ploidyStep: 0.25))
            {
                FitPurityPloidy(segs, purityStep: 0.05, ploidyStep: 0.25).Segments.Should()
                    .HaveCount(n, "one integer segment is emitted per summary, seed {0}", seed);
            }
        }
    }

    // ── BE: zero-span segment (End == Start) ⇒ emitted with a 1 bp span (§4.2) ────
    [Test]
    [CancelAfter(30_000)]
    public void EvaluatePurityPloidy_ZeroSpanSegment_EmittedWithPositiveLength()
    {
        var seg = new[] { new AlleleSpecificSegmentSummary("1", 500, 500, 0.0, 0.5, 1) };

        var fit = EvaluatePurityPloidy(seg, 0.8, 2.0);

        fit.Segments[0].Length.Should().Be(1, "a single-position summary is widened to a 1 bp span (§4.2)");
    }

    #endregion

    #region ONCO-ASCAT-001 — BE: DeriveMultiplicity boundary clamp (INV-02)

    // ── BE: m ∈ [1, major] over a randomised sweep of boundary VAF/purity/CN ─────
    // INV-02 / §2.2: m = clamp(round(n_mut), 1, major). VAF=0 ⇒ rounds to 0 ⇒ clamped
    // up to 1; a huge n_mut ⇒ clamped down to major. Never < 1, never > major, never a
    // crash, whatever the boundary inputs.
    [Test]
    [CancelAfter(20_000)]
    public void DeriveMultiplicity_BoundaryInputs_AlwaysWithinOneToMajor()
    {
        for (int seed = 0; seed < 2_000; seed++)
        {
            var rng = new Random(seed);
            double vaf = rng.Next(3) switch { 0 => 0.0, 1 => 1.0, _ => rng.NextDouble() };
            double purity = rng.Next(3) switch { 0 => 1.0, 1 => 1e-6, _ => rng.NextDouble() }; // (0,1]
            if (purity <= 0.0) purity = 1e-6;
            int total = 1 + rng.Next(12);
            int major = 1 + rng.Next(total); // [1, total]

            int m = DeriveMultiplicity(vaf, purity, total, major);

            m.Should().BeGreaterThanOrEqualTo(1, "a variant sits on ≥ 1 copy (INV-02), seed {0}", seed);
            m.Should().BeLessThanOrEqualTo(major, "m ≤ major-allele copy number (INV-02), seed {0}", seed);
        }
    }

    // ── BE: documented boundary results ──────────────────────────────────────────
    [Test]
    public void DeriveMultiplicity_DocumentedBoundaries_ClampToOneAndMajor()
    {
        DeriveMultiplicity(vaf: 0.0, purity: 0.6, totalCopyNumber: 4, majorCopyNumber: 3)
            .Should().Be(1, "VAF rounds to 0 ⇒ clamp up to 1 (§6.1)");
        DeriveMultiplicity(vaf: 1.0, purity: 0.6, totalCopyNumber: 4, majorCopyNumber: 2)
            .Should().Be(2, "a large n_mut ⇒ clamp down to major CN 2 (§6.1)");
    }

    // ── BE: out-of-range multiplicity arguments ⇒ ArgumentOutOfRangeException ─────
    [Test]
    public void DeriveMultiplicity_OutOfRangeArguments_Throw()
    {
        ((Action)(() => DeriveMultiplicity(-0.1, 0.6, 2, 1))).Should().Throw<ArgumentOutOfRangeException>("vaf < 0");
        ((Action)(() => DeriveMultiplicity(1.1, 0.6, 2, 1))).Should().Throw<ArgumentOutOfRangeException>("vaf > 1");
        ((Action)(() => DeriveMultiplicity(double.NaN, 0.6, 2, 1))).Should().Throw<ArgumentOutOfRangeException>("vaf NaN");
        ((Action)(() => DeriveMultiplicity(0.4, 0.0, 2, 1))).Should().Throw<ArgumentOutOfRangeException>("purity 0");
        ((Action)(() => DeriveMultiplicity(0.4, 1.5, 2, 1))).Should().Throw<ArgumentOutOfRangeException>("purity > 1");
        ((Action)(() => DeriveMultiplicity(0.4, 0.6, 0, 1))).Should().Throw<ArgumentOutOfRangeException>("total < 1");
        ((Action)(() => DeriveMultiplicity(0.4, 0.6, 2, 0))).Should().Throw<ArgumentOutOfRangeException>("major < 1");
        ((Action)(() => DeriveMultiplicity(0.4, 0.6, 2, 3))).Should().Throw<ArgumentOutOfRangeException>("major > total");
    }

    #endregion
}
