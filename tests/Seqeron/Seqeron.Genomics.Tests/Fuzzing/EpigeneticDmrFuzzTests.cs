using static Seqeron.Genomics.Annotation.EpigeneticsAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Epigenetics area — differentially methylated region (DMR)
/// detection (EPIGEN-DMR-001). The public entry point under test lives in
/// <see cref="EpigeneticsAnalyzer"/>:
///   • <see cref="EpigeneticsAnalyzer.FindDMRs(IEnumerable{MethylationSite}, IEnumerable{MethylationSite}, int, double, int)"/>
///     — compares two single-sample methylation profiles using the methylKit
///     tiling-window model: positions are grouped into fixed-width windows
///     [start, start+windowSize); within each window the per-site differences
///     (sample2 − sample1, in fraction units) are averaged, and a region is
///     reported when it has ≥ minCpGCount covered sites and the absolute mean
///     difference STRICTLY exceeds minDifference. Each region carries a two-sided
///     Fisher's exact p-value and a hyper/hypo label by the sign of the mean.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary inputs to a unit and asserts that the
/// code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// nonsense output, and no *unhandled* runtime exception (DivideByZero on a
/// single-site/empty mean, IndexOutOfRange, NullReference). Every input must
/// resolve to EITHER a well-defined, theory-correct value OR a *documented,
/// intentional* outcome. For DMR detection the headline hazards are:
///   • a FALSE DMR on identical methylomes (difference 0 everywhere → no region);
///   • an OFF-BY-ONE at the Δ threshold (the cutoff is STRICT `|meanDiff| > Δ`,
///     so a region exactly at the threshold must NOT be reported);
///   • a DivideByZero / crash on a single-site window or an empty input
///     (the per-window mean divides by the site count);
///   • a NaN mean difference, or a p-value outside [0,1].
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: EPIGEN-DMR-001 — Differentially Methylated Regions (Epigenetics)
/// Checklist: docs/checklists/03_FUZZING.md, row 184.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, empty.
///     Targets (checklist row 184): "identical methylomes, threshold edge,
///     single site".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The DMR contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Region difference (methylKit; Akalin et al. 2012):
///   meanDiff = mean over the window of (level2 − level1), fraction units.
///   A window is a DMR iff it has ≥ minCpGCount covered sites AND
///   |meanDiff| > minDifference (STRICT), and is labelled "Hypermethylated"
///   when meanDiff > 0 (treatment higher than control), "Hypomethylated" when < 0.
///   — docs/algorithms/Epigenetics/Differentially_Methylated_Regions.md §2.2, §4.1.
///
/// Documented invariants exercised here
/// (docs/algorithms/Epigenetics/Differentially_Methylated_Regions.md §2.4):
///   • INV-01: every reported DMR has |MeanDifference| > minDifference (strict).
///   • INV-02: Annotation is "Hypermethylated" iff MeanDifference > 0,
///             "Hypomethylated" iff < 0.
///   • INV-03: 0 ≤ PValue ≤ 1.
///   • INV-04: positions ≥ windowSize apart fall in different windows.
///   • INV-05: a reported DMR contains ≥ minCpGCount covered sites.
///   • INV-06: output is ordered by Start ascending and is deterministic.
///
/// Boundary / degenerate handling (docs §3.3, §6.1):
///   • empty input                          → no DMRs.
///   • window with < minCpGCount sites       → not reported (no DivideByZero).
///   • identical methylomes                  → meanDiff 0 everywhere → no DMRs.
///   • |meanDiff| == minDifference exactly    → NOT reported (strict `>`).
///   • zero coverage in one group in a window → Fisher p = 1.0 (degenerate margin).
///   • a position present in only one sample  → compared against implicit 0.
///   • null sample                           → ArgumentNullException.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class EpigeneticDmrFuzzTests
{
    private const double DefaultMinDifference = 0.25;
    private const int DefaultMinCpGCount = 3;

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins INV-01..INV-06 against an independent re-derivation of the documented
    // rule. For every reported DMR:
    //   • |MeanDifference| > minDifference (INV-01, strict);
    //   • the mean is a finite number in [-1, 1] (no NaN/Inf);
    //   • Annotation is exactly "Hypermethylated" for a positive mean and
    //     "Hypomethylated" for a negative mean (INV-02);
    //   • 0 ≤ PValue ≤ 1 (INV-03);
    //   • CpGCount ≥ minCpGCount (INV-05);
    //   • Start ≤ End and both are non-negative coordinates in the input space.
    // The list as a whole must be ordered by Start ascending (INV-06).
    // This is what stops a test from rubber-stamping a spurious region green.
    private static void AssertWellFormedDmrs(
        IReadOnlyList<DifferentiallyMethylatedRegion> dmrs,
        double minDifference,
        int minCpGCount)
    {
        for (int i = 0; i < dmrs.Count; i++)
        {
            var d = dmrs[i];

            double.IsNaN(d.MeanDifference).Should().BeFalse("a DMR mean difference must never be NaN");
            double.IsInfinity(d.MeanDifference).Should().BeFalse("a DMR mean difference must be finite");
            d.MeanDifference.Should().BeInRange(-1.0, 1.0, "the mean of fraction differences lies in [-1,1]");

            Math.Abs(d.MeanDifference).Should().BeGreaterThan(minDifference,
                "INV-01: the reporting cutoff |meanDiff| > minDifference is strict");

            if (d.MeanDifference > 0)
                d.Annotation.Should().Be("Hypermethylated", "INV-02: positive mean → hyper");
            else
                d.Annotation.Should().Be("Hypomethylated", "INV-02: negative mean → hypo");

            d.PValue.Should().BeInRange(0.0, 1.0, "INV-03: a Fisher p-value lies in [0,1]");

            d.CpGCount.Should().BeGreaterThanOrEqualTo(minCpGCount,
                "INV-05: a DMR is a region of ≥ minCpGCount covered sites");

            d.Start.Should().BeLessThanOrEqualTo(d.End, "a region spans Start..End");
            d.Start.Should().BeGreaterThanOrEqualTo(0, "positions are non-negative 0-based coordinates");

            if (i > 0)
                d.Start.Should().BeGreaterThanOrEqualTo(dmrs[i - 1].Start,
                    "INV-06: DMRs are ordered by Start ascending");
        }
    }

    // Builds a contiguous run of CpG sites (positions start..start+count-1, all in
    // one default window) at a fixed methylation level and coverage.
    private static List<MethylationSite> Run(int start, int count, double level, int coverage) =>
        Enumerable.Range(start, count)
            .Select(p => new MethylationSite(p, MethylationType.CpG, "CG", level, coverage))
            .ToList();

    // ─────────────────────────────────────────────────────────────────────────
    #region EPIGEN-DMR-001 — Differentially Methylated Regions (FindDMRs)

    // ── BE: positive sanity — a clear DMR is detected with the documented sign ──
    // The worked example (docs §7.1): control β=0.0, treatment β=1.0 over 3 CpGs
    // → MeanDifference == 1.0, "Hypermethylated".
    [Test]
    public void FindDMRs_ClearHyperRegion_ReportedWithDocumentedDifferenceAndDirection()
    {
        var control = Run(start: 0, count: 5, level: 0.1, coverage: 20);
        var treatment = Run(start: 0, count: 5, level: 0.9, coverage: 20);

        var dmrs = FindDMRs(control, treatment).ToList();

        AssertWellFormedDmrs(dmrs, DefaultMinDifference, DefaultMinCpGCount);
        dmrs.Should().HaveCount(1, "the five adjacent CpGs form a single window above the cutoff");
        dmrs[0].MeanDifference.Should().BeApproximately(0.8, 1e-9, "mean(0.9 − 0.1) = 0.8");
        dmrs[0].Annotation.Should().Be("Hypermethylated", "treatment β > control β");
        dmrs[0].CpGCount.Should().Be(5);
    }

    [Test]
    public void FindDMRs_ClearHypoRegion_ReportedAsHypomethylated()
    {
        var control = Run(start: 0, count: 4, level: 0.95, coverage: 30);
        var treatment = Run(start: 0, count: 4, level: 0.05, coverage: 30);

        var dmrs = FindDMRs(control, treatment).ToList();

        AssertWellFormedDmrs(dmrs, DefaultMinDifference, DefaultMinCpGCount);
        dmrs.Should().ContainSingle();
        dmrs[0].MeanDifference.Should().BeApproximately(-0.9, 1e-9, "mean(0.05 − 0.95) = −0.9");
        dmrs[0].Annotation.Should().Be("Hypomethylated", "treatment β < control β");
    }

    // ── BE: identical methylomes — the canonical negative (no false positive) ───
    [Test]
    public void FindDMRs_IdenticalMethylomes_ReportsNoDmr()
    {
        var sites = Run(start: 0, count: 6, level: 0.6, coverage: 25);

        // sample1 and sample2 are byte-for-byte identical → diff 0 at every site.
        var dmrs = FindDMRs(sites, sites.ToList()).ToList();

        AssertWellFormedDmrs(dmrs, DefaultMinDifference, DefaultMinCpGCount);
        dmrs.Should().BeEmpty("identical methylomes have zero difference everywhere");
    }

    [Test]
    [CancelAfter(20000)]
    public void FindDMRs_IdenticalRandomMethylomes_NeverReportsDmr([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        int count = rng.Next(1, 40);
        var sites = Enumerable.Range(0, count)
            .Select(p => new MethylationSite(
                p,
                MethylationType.CpG,
                "CG",
                rng.NextDouble(),
                rng.Next(1, 100)))
            .ToList();

        // Pass the SAME sites as both samples (a fresh list, identical values).
        var dmrs = FindDMRs(sites, sites.ToList()).ToList();

        AssertWellFormedDmrs(dmrs, DefaultMinDifference, DefaultMinCpGCount);
        dmrs.Should().BeEmpty("a methylome compared against itself has no differential region");
    }

    // ── BE: threshold edge — strict `>` cutoff, no off-by-one ───────────────────
    [Test]
    public void FindDMRs_MeanDifferenceExactlyAtThreshold_NotReported()
    {
        // Three sites, each with diff exactly 0.25 → meanDiff == 0.25 == minDifference.
        var control = Run(start: 0, count: 3, level: 0.50, coverage: 40);
        var treatment = Run(start: 0, count: 3, level: 0.75, coverage: 40);

        var dmrs = FindDMRs(control, treatment, minDifference: 0.25).ToList();

        AssertWellFormedDmrs(dmrs, 0.25, DefaultMinCpGCount);
        dmrs.Should().BeEmpty("the cutoff is strict |meanDiff| > Δ, so equality is excluded");
    }

    [Test]
    public void FindDMRs_MeanDifferenceJustAboveThreshold_Reported()
    {
        // meanDiff = 0.25 + tiny ε → strictly above the threshold → reported.
        var control = Run(start: 0, count: 3, level: 0.50, coverage: 40);
        var treatment = control
            .Select(s => s with { MethylationLevel = 0.50 + 0.25 + 1e-6 })
            .ToList();

        var dmrs = FindDMRs(control, treatment, minDifference: 0.25).ToList();

        AssertWellFormedDmrs(dmrs, 0.25, DefaultMinCpGCount);
        dmrs.Should().ContainSingle("a difference strictly above the cutoff is a DMR");
        dmrs[0].Annotation.Should().Be("Hypermethylated");
    }

    [Test]
    public void FindDMRs_NegativeMeanExactlyAtThreshold_NotReported()
    {
        // Symmetric negative edge: meanDiff == −0.25, |−0.25| == 0.25 == Δ → excluded.
        var control = Run(start: 0, count: 3, level: 0.75, coverage: 40);
        var treatment = Run(start: 0, count: 3, level: 0.50, coverage: 40);

        var dmrs = FindDMRs(control, treatment, minDifference: 0.25).ToList();

        dmrs.Should().BeEmpty("strict cutoff excludes |meanDiff| == Δ on the negative side too");
    }

    // ── BE: single site — minCpGCount requirement, no DivideByZero ──────────────
    [Test]
    public void FindDMRs_SingleCpG_BelowDefaultMinCount_NotReportedAndNoCrash()
    {
        // A single, maximally-different CpG: the window has 1 site < default 3.
        var control = Run(start: 100, count: 1, level: 0.0, coverage: 50);
        var treatment = Run(start: 100, count: 1, level: 1.0, coverage: 50);

        var act = () => FindDMRs(control, treatment).ToList();

        var dmrs = act.Should().NotThrow("a 1-site window must not crash the per-window mean").Subject;
        dmrs.Should().BeEmpty("a single site is below the default minCpGCount of 3 (INV-05)");
    }

    [Test]
    public void FindDMRs_SingleSite_WithMinCpGCountOne_ReportsAndDoesNotDivideByZero()
    {
        // minCpGCount=1 makes a lone site a valid region; the per-window mean
        // divides by the site count — must be exactly the single difference, no NaN.
        var control = Run(start: 0, count: 1, level: 0.1, coverage: 50);
        var treatment = Run(start: 0, count: 1, level: 0.9, coverage: 50);

        var dmrs = FindDMRs(control, treatment, minCpGCount: 1).ToList();

        AssertWellFormedDmrs(dmrs, DefaultMinDifference, minCpGCount: 1);
        dmrs.Should().ContainSingle();
        dmrs[0].MeanDifference.Should().BeApproximately(0.8, 1e-9, "single-site mean = the lone diff");
        dmrs[0].CpGCount.Should().Be(1);
        dmrs[0].Start.Should().Be(0);
        dmrs[0].End.Should().Be(0);
    }

    [Test]
    public void FindDMRs_SingleSite_PresentInOneSampleOnly_ComparedAgainstImplicitZero()
    {
        // docs §3.3: a position present in only one sample is compared to implicit 0.
        var control = new List<MethylationSite>();
        var treatment = Run(start: 0, count: 1, level: 1.0, coverage: 50);

        var dmrs = FindDMRs(control, treatment, minCpGCount: 1).ToList();

        AssertWellFormedDmrs(dmrs, DefaultMinDifference, minCpGCount: 1);
        dmrs.Should().ContainSingle();
        dmrs[0].MeanDifference.Should().BeApproximately(1.0, 1e-9, "1.0 − implicit 0.0 = 1.0");
        dmrs[0].Annotation.Should().Be("Hypermethylated");
    }

    // ── BE: empty / degenerate inputs ───────────────────────────────────────────
    [Test]
    public void FindDMRs_BothSamplesEmpty_ReturnsNoDmrAndNoCrash()
    {
        var act = () => FindDMRs(
            new List<MethylationSite>(),
            new List<MethylationSite>()).ToList();

        var dmrs = act.Should().NotThrow().Subject;
        dmrs.Should().BeEmpty("no positions → no tiles → no DMRs");
    }

    [Test]
    public void FindDMRs_NullSample_ThrowsArgumentNullException()
    {
        var sites = Run(start: 0, count: 3, level: 0.5, coverage: 10);

        // Validation is eager (not deferred to enumeration) per the source guard.
        var actNull1 = () => FindDMRs(null!, sites);
        var actNull2 = () => FindDMRs(sites, null!);

        actNull1.Should().Throw<ArgumentNullException>();
        actNull2.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void FindDMRs_ZeroCoverageInOneGroup_FisherPValueIsOne()
    {
        // docs §6.1: zero coverage in a group within a window → degenerate margin,
        // Fisher p = 1.0. The level difference still drives the DMR call.
        var control = Run(start: 0, count: 3, level: 0.0, coverage: 0);
        var treatment = Run(start: 0, count: 3, level: 0.9, coverage: 30);

        var dmrs = FindDMRs(control, treatment).ToList();

        AssertWellFormedDmrs(dmrs, DefaultMinDifference, DefaultMinCpGCount);
        dmrs.Should().ContainSingle();
        dmrs[0].PValue.Should().Be(1.0, "a zero row total is a degenerate fixed margin");
        dmrs[0].MeanDifference.Should().BeApproximately(0.9, 1e-9);
    }

    // ── BE: windowSize tiling boundary (INV-04) ─────────────────────────────────
    [Test]
    public void FindDMRs_PositionsAtWindowBoundary_SplitIntoSeparateWindows()
    {
        // windowSize=10: positions 0,1,2 in window [0,10); position 10 opens a new
        // window. The first window has 3 sites (a DMR); the lone site at 10 does not.
        var control = new List<MethylationSite>
        {
            new(0, MethylationType.CpG, "CG", 0.0, 20),
            new(1, MethylationType.CpG, "CG", 0.0, 20),
            new(2, MethylationType.CpG, "CG", 0.0, 20),
            new(10, MethylationType.CpG, "CG", 0.0, 20),
        };
        var treatment = control.Select(s => s with { MethylationLevel = 1.0 }).ToList();

        var dmrs = FindDMRs(control, treatment, windowSize: 10).ToList();

        AssertWellFormedDmrs(dmrs, DefaultMinDifference, DefaultMinCpGCount);
        dmrs.Should().ContainSingle("only the 3-site window qualifies; the boundary site is alone");
        dmrs[0].Start.Should().Be(0);
        dmrs[0].End.Should().Be(2);
    }

    // ── BE: determinism (INV-06) ────────────────────────────────────────────────
    [Test]
    [CancelAfter(20000)]
    public void FindDMRs_RandomMethylomes_AreWellFormedAndDeterministic([Random(1, 1_000_000, 30)] int seed)
    {
        var rng = new Random(seed);
        int count = rng.Next(0, 60);

        var control = Enumerable.Range(0, count)
            .Select(p => new MethylationSite(p, MethylationType.CpG, "CG", rng.NextDouble(), rng.Next(0, 80)))
            .ToList();
        var treatment = Enumerable.Range(0, count)
            .Select(p => new MethylationSite(p, MethylationType.CpG, "CG", rng.NextDouble(), rng.Next(0, 80)))
            .ToList();

        var first = FindDMRs(control, treatment).ToList();
        var second = FindDMRs(control, treatment).ToList();

        AssertWellFormedDmrs(first, DefaultMinDifference, DefaultMinCpGCount);
        first.Should().Equal(second, "INV-06: detection has no randomness; output is deterministic");
    }

    #endregion
}
