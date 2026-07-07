namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the PanGenome area — Heaps' law pan-genome openness estimation
/// (PANGEN-HEAP-001), the Tettelin (2008) / micropan <c>heaps()</c> power-law model
/// that fits <c>n(N) = K · N^(−α)</c> to the new-gene-discovery curve and decides
/// open (α &lt; 1) vs closed (α &gt; 1).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to the fit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop,
/// no NaN / Infinity leaking out of the log-log / power-law regression, and no
/// *unhandled* runtime exception. The headline numeric hazards for a power-law fit
/// are exactly the ones called out by the BE checklist targets:
///   • a regression over a SINGLE data point (1 genome ⇒ the curve N = 2..G is
///     EMPTY ⇒ there is nothing to fit ⇒ a DivideByZero / 0-point regression /
///     NaN must NOT escape — the documented guard is a degenerate fit, not a throw);
///   • a regression over a FLAT (zero-variance) curve (all genomes identical ⇒
///     every genome after the first contributes ZERO new clusters ⇒ y ≡ 0, a
///     constant ZERO curve ⇒ no log(0) blow-up, no NaN — the objective
///     J(K,α) = K·sqrt(Σ x^(−2α))/|x| is uniquely minimized at K = 0 where α is
///     UNIDENTIFIABLE (the model is 0 for every α); micropan's optim, started at
///     (mean y at N=2, 1) = (0, 1), cannot improve, so it leaves α = 1 ⇒ the
///     saturated pan-genome is CLOSED — which is the biologically correct verdict
///     for genomes that add no new genes (the "α = 0, open" §6.1 row is the
///     constant-POSITIVE case, not the zero curve));
///   • log(0) when a new-gene count is 0 (the objective is evaluated directly on
///     y, not on log y, so a 0 count must never feed a log).
/// Every input must resolve to EITHER a well-defined, theory-correct fit (finite
/// K ∈ [0,10000], finite α ∈ [0,2], IsOpen ⇔ α &lt; 1, predictor finite and
/// non-increasing), OR the *documented, intentional* degenerate fit
/// (0, 0, false, predictor→0) for &lt; 2 genomes / null / empty input. A raw
/// exception, a hang, a NaN/Inf α or K, an out-of-bounds parameter, or an
/// IsOpen flag inconsistent with α is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PANGEN-HEAP-001 — Pan-genome growth model (Heaps' law openness)
/// Checklist: docs/checklists/03_FUZZING.md, row 192.
/// Algorithm doc: docs/algorithms/PanGenome/Pan_Genome_Growth_Model.md
///                (Test Unit ID PANGEN-HEAP-001).
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate genome-count / curve boundaries
///          called out in the checklist row:
///            – SINGLE genome (N=1) → the new-gene curve N = 2..G is EMPTY (one data
///              point at most, in fact zero fitted points) → the fit is undefined →
///              documented degenerate result (0, 0, false, predictor→0); NO
///              DivideByZero / NaN from a regression over &lt; 2 points, no log(0)
///              (§3.3, §6.1 "0 or 1 genome").
///            – 2 genomes → the new-gene curve has exactly ONE fitted point (N=2);
///              the MINIMUM data on which the power fit runs → a finite, in-bounds
///              (K, α), IsOpen consistent with α, no crash (§6.1, §2.2).
///            – IDENTICAL genomes → every genome after the first adds zero new
///              clusters → a FLAT (zero-variance) ZERO curve y ≡ 0 → the unique
///              optimum is K = 0 (= mean count) with α UNIDENTIFIABLE; micropan's
///              optim leaves α at its start 1 ⇒ the saturated pan-genome is CLOSED
///              (adding genomes yields no new genes); NO log(0) / NaN from the flat
///              curve. A flat POSITIVE curve (each genome adds the SAME k > 0 new
///              genes) is the genuine §6.1 case: α = 0, K = k, OPEN (tested
///              separately).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes:
///    "BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty)").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The Heaps'-law contract under test (Pan_Genome_Growth_Model.md §2–§3)
/// ───────────────────────────────────────────────────────────────────────────
/// For a binary gene presence/absence matrix (genomes × clusters), order the
/// genomes and record, for the N-th genome (N = 2..G), the number of clusters that
/// appear for the FIRST time at position N (INV-02). Heaps' law fits this new-gene
/// curve as a power law n(N) = K · N^(−α) by minimizing the micropan objective
///   J(K, α) = sqrt( Σ_i (y_i − K · x_i^(−α))² ) / |x|
/// over box constraints K ∈ [0, 10000], α ∈ [0, 2] from start (mean y at N=2, 1)
/// (§2.2, INV-04). The verdict is the strict threshold rule:
///   • INV-01: open ⇔ α &lt; 1; closed ⇔ α &gt; 1 (boundary α = 1 ⇒ not-open).
///   • INV-05: for points exactly on a single power curve within bounds, the
///     recovered (K, α) equal the analytic solution (unique J = 0 minimum).
///   • INV-06: predictor(N) = K · N^(−α) is non-increasing in N for α ≥ 0.
/// The API entries under test are
///   PanGenomeAnalyzer.FitHeapsLaw(IEnumerable&lt;GenePresenceRow&gt; matrix, int permutations = 100)
///   PanGenomeAnalyzer.FitHeapsLaw(IReadOnlyDictionary&lt;...&gt; genomes, double, int)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs,
///    FitHeapsLaw lines 575–599; degenerate guard EmptyHeapsFit lines 708–713;
///    objective + bounded minimizer lines 716–774).
///
/// HAND-CHECKABLE worked example (Pan_Genome_Growth_Model.md §7.1): a fixed-order
/// 3-genome matrix where genome 2 introduces 8 clusters absent in genome 1 and
/// genome 3 introduces 4 clusters absent in genomes 1–2 gives the new-gene curve
/// x = [2, 3], y = [8, 4]. Solving K·2^(−α) = 8, K·3^(−α) = 4:
///   α = ln2 / ln(3/2) = 1.7095113,  K = 8·2^α = 26.1640014.
/// Since α &gt; 1 the pan-genome is CLOSED. Verified below with permutations = 1
/// (natural input order, §5.2) so the fit is exactly reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PanGenomeHeapsFuzzTests
{
    #region Helpers

    private const int Seed = 192_0001; // local-only deterministic seed for the fuzz sweep

    // Heaps' fit bounds (Pan_Genome_Growth_Model.md §2.2 / INV-04).
    private const double KMin = 0.0;
    private const double KMax = 10000.0;
    private const double AlphaMin = 0.0;
    private const double AlphaMax = 2.0;
    private const double OpenThreshold = 1.0;

    /// <summary>
    /// Builds a single presence/absence row from an explicit set of present cluster
    /// IDs over a fixed universe of cluster IDs. A cluster is present iff its ID is in
    /// <paramref name="present"/>. This lets a test hand-craft an exact new-gene curve
    /// (which clusters first appear at which genome) without going through clustering.
    /// </summary>
    private static PanGenomeAnalyzer.GenePresenceRow Row(
        string genomeId, IReadOnlyList<string> universe, ISet<string> present)
    {
        var presence = universe.ToDictionary(c => c, c => present.Contains(c));
        return new PanGenomeAnalyzer.GenePresenceRow(
            GenomeId: genomeId,
            GenePresence: presence,
            TotalGenes: universe.Count,
            PresentGenes: present.Count);
    }

    private static IReadOnlyList<string> Universe(int n) =>
        Enumerable.Range(0, n).Select(i => $"c{i}").ToList();

    private static HashSet<string> Clusters(IEnumerable<int> ids) =>
        new(ids.Select(i => $"c{i}"));

    private static (string GeneId, string Sequence) Gene(string id, string seq) => (id, seq);

    // A pool of distinct 30-bp sequences keyed by an integer "family" so the k-mer
    // (k=7) clusterer assigns distinct families to distinct clusters and identical
    // families to the same cluster (mirrors PanGenomeCoreFuzzTests.Family).
    private static readonly string[] Bases = { "A", "C", "G", "T" };
    private static string Family(int f)
    {
        var rng = new Random(unchecked(1_000_003 * (f + 1)));
        var sb = new System.Text.StringBuilder(30);
        for (int i = 0; i < 30; i++)
            sb.Append(Bases[rng.Next(4)]);
        return sb.ToString();
    }

    /// <summary>
    /// Asserts the result is a WELL-FORMED Heaps' fit per Pan_Genome_Growth_Model.md
    /// §3.2 / §2.4: K and α are FINITE (no NaN/Inf escaping the regression — the
    /// headline fuzz hazard), both lie inside their documented box bounds (INV-04),
    /// the IsOpen flag is consistent with the strict α &lt; 1 rule (INV-01), and the
    /// predictor is finite and non-increasing on N = 1..50 (INV-06). Used by every
    /// test so a NaN α, an out-of-bounds K, or an inconsistent flag fails everywhere.
    /// </summary>
    private static void AssertWellFormedFit(PanGenomeAnalyzer.HeapsLawFit fit)
    {
        // Finiteness: the central hazard — a 1-point / flat-curve / log(0) regression
        // must never leak NaN or Infinity into K or α.
        double.IsNaN(fit.Intercept).Should().BeFalse("Intercept K must be finite, never NaN.");
        double.IsInfinity(fit.Intercept).Should().BeFalse("Intercept K must be finite, never Infinity.");
        double.IsNaN(fit.Alpha).Should().BeFalse("Alpha α must be finite, never NaN.");
        double.IsInfinity(fit.Alpha).Should().BeFalse("Alpha α must be finite, never Infinity.");

        // INV-04: parameters stay inside the micropan box bounds.
        fit.Intercept.Should().BeInRange(KMin, KMax, "INV-04: Intercept K ∈ [0, 10000].");
        fit.Alpha.Should().BeInRange(AlphaMin, AlphaMax, "INV-04: Alpha α ∈ [0, 2].");

        // INV-01: open ⇔ α < 1 (strict; boundary α = 1 ⇒ not-open). The ONE documented
        // exception (§3.3, §6.1) is the degenerate sentinel fit (K=0, α=0, IsOpen=false)
        // returned for < 2 genomes / null / empty input: there α = 0 but IsOpen is the
        // documented `false` sentinel, NOT a genuine "open" verdict. Outside that
        // sentinel the open/closed flag must match α exactly.
        bool isDegenerateSentinel = fit.Intercept == 0.0 && fit.Alpha == 0.0 && !fit.IsOpen;
        if (!isDegenerateSentinel)
        {
            fit.IsOpen.Should().Be(fit.Alpha < OpenThreshold,
                "INV-01: IsOpen must be exactly (α < 1) — the open/closed verdict must match α.");
        }

        // INV-06: predictor(N) = K·N^(−α) is finite and non-increasing in N for α ≥ 0.
        double prev = double.PositiveInfinity;
        for (int n = 1; n <= 50; n++)
        {
            double p = fit.PredictNewGenes(n);
            double.IsNaN(p).Should().BeFalse($"predictor({n}) must be finite, never NaN.");
            double.IsInfinity(p).Should().BeFalse($"predictor({n}) must be finite, never Infinity.");
            p.Should().BeLessThanOrEqualTo(prev + 1e-9,
                $"INV-06: predictor must be non-increasing in N (α ≥ 0); failed at N = {n}.");
            prev = p;
        }
    }

    #endregion

    #region PANGEN-HEAP-001 — Heaps' law pan-genome openness

    #region Positive sanity — documented worked example & known curves

    // Pan_Genome_Growth_Model.md §7.1 worked example, made hand-checkable: a
    // fixed-order 3-genome matrix whose new-gene curve is x=[2,3], y=[8,4].
    //   genome 1 (g1): clusters {0..2}        (3 shared/core clusters, present in all)
    //   genome 2 (g2): {0..2} ∪ {3..10}       (8 NEW clusters at N=2)
    //   genome 3 (g3): {0..2} ∪ {3..10} ∪ {11..14}  (4 NEW clusters at N=3)
    // Analytic fit: α = ln2/ln(3/2) = 1.7095113, K = 8·2^α = 26.1640014 ⇒ CLOSED.
    // permutations = 1 ⇒ natural input order only (§5.2) ⇒ exactly reproducible.
    [Test]
    [CancelAfter(20000)]
    public void Sanity_DocWorkedExample_RecoversAnalyticAlphaAndK_Closed()
    {
        var universe = Universe(15); // c0..c14
        var rows = new[]
        {
            Row("g1", universe, Clusters(Enumerable.Range(0, 3))),                 // {0,1,2}
            Row("g2", universe, Clusters(Enumerable.Range(0, 11))),                // {0..10} → +8 new
            Row("g3", universe, Clusters(Enumerable.Range(0, 15))),                // {0..14} → +4 new
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);
        AssertWellFormedFit(fit);

        double expectedAlpha = Math.Log(2.0) / Math.Log(1.5);   // 1.7095113
        double expectedK = 8.0 * Math.Pow(2.0, expectedAlpha);  // 26.1640014

        fit.Alpha.Should().BeApproximately(expectedAlpha, 1e-4,
            "INV-05: the curve x=[2,3], y=[8,4] lies on an exact power law ⇒ analytic α = ln2/ln(3/2).");
        fit.Intercept.Should().BeApproximately(expectedK, 1e-3,
            "INV-05: analytic K = 8·2^α for the worked example.");
        fit.IsOpen.Should().BeFalse("§7.1: α = 1.71 > 1 ⇒ the pan-genome is CLOSED.");

        // Predictor reproduces the fitted curve at the fitted x-points.
        fit.PredictNewGenes(2).Should().BeApproximately(8.0, 1e-2, "predictor(2) = K·2^(−α) = 8.");
        fit.PredictNewGenes(3).Should().BeApproximately(4.0, 1e-2, "predictor(3) = K·3^(−α) = 4.");
    }

    // A genome set whose new-gene curve lies on an exact OPEN power law (α < 1).
    // Choose α = 0.5: y(2) = K·2^(−0.5), y(3) = K·3^(−0.5). With K = 20:
    //   y(2) = 20/√2 ≈ 14.142  → 14 new,   y(3) = 20/√3 ≈ 11.547 → 12 new.
    // We instead pick exact integer-friendly counts on a shallow-decay curve and
    // assert the openness CLASSIFICATION (the load-bearing verdict), with α < 1.
    [Test]
    [CancelAfter(20000)]
    public void Sanity_ShallowDecayCurve_ClassifiedOpen_AlphaBelowOne()
    {
        // new-gene curve x=[2,3,4], y=[10,9,8] — a slow, near-linear decay ⇒ α < 1.
        var universe = Universe(40);
        var rows = new[]
        {
            Row("g1", universe, Clusters(Enumerable.Range(0, 5))),    // {0..4}
            Row("g2", universe, Clusters(Enumerable.Range(0, 15))),   // +10 new at N=2
            Row("g3", universe, Clusters(Enumerable.Range(0, 24))),   // +9  new at N=3
            Row("g4", universe, Clusters(Enumerable.Range(0, 32))),   // +8  new at N=4
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);
        AssertWellFormedFit(fit);

        fit.Alpha.Should().BeLessThan(OpenThreshold,
            "a slow new-gene decay (10,9,8) fits a shallow power law α < 1.");
        fit.IsOpen.Should().BeTrue("α < 1 ⇒ OPEN pan-genome (INV-01).");
    }

    // A genome set whose new-gene curve decays steeply ⇒ α > 1 ⇒ CLOSED.
    [Test]
    [CancelAfter(20000)]
    public void Sanity_SteepDecayCurve_ClassifiedClosed_AlphaAboveOne()
    {
        // new-gene curve x=[2,3,4], y=[16,4,1] — steep decay ⇒ α > 1 (≈ halving twice).
        var universe = Universe(30);
        var rows = new[]
        {
            Row("g1", universe, Clusters(Enumerable.Range(0, 4))),    // {0..3}
            Row("g2", universe, Clusters(Enumerable.Range(0, 20))),   // +16 new at N=2
            Row("g3", universe, Clusters(Enumerable.Range(0, 24))),   // +4  new at N=3
            Row("g4", universe, Clusters(Enumerable.Range(0, 25))),   // +1  new at N=4
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);
        AssertWellFormedFit(fit);

        fit.Alpha.Should().BeGreaterThan(OpenThreshold,
            "a steep new-gene decay (16,4,1) fits a steep power law α > 1.");
        fit.IsOpen.Should().BeFalse("α > 1 ⇒ CLOSED pan-genome (INV-01).");
    }

    #endregion

    #region BE — Boundary: single genome (N=1 ⇒ empty curve ⇒ degenerate guard)

    // §3.3 / §6.1 "0 or 1 genome": with one genome the curve N = 2..G is EMPTY, so
    // there is nothing to fit. The documented guard is a DEGENERATE fit
    // (0, 0, IsOpen=false, predictor→0) — NOT a DivideByZero, a 1-point regression
    // NaN, or a log(0). This is the headline single-data-point hazard for the row.
    [Test]
    [CancelAfter(20000)]
    public void Single_OneGenome_ReturnsDegenerateFit_NoRegressionOverOnePoint()
    {
        var universe = Universe(5);
        var rows = new[] { Row("only", universe, Clusters(new[] { 0, 1, 2, 3, 4 })) };

        var act = () => PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 50);
        act.Should().NotThrow("§3.3: 1 genome ⇒ empty curve ⇒ degenerate fit, never a throw.");

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 50);
        AssertWellFormedFit(fit);

        fit.Intercept.Should().Be(0.0, "§6.1: 1 genome ⇒ degenerate fit Intercept = 0.");
        fit.Alpha.Should().Be(0.0, "§6.1: 1 genome ⇒ degenerate fit Alpha = 0.");
        fit.IsOpen.Should().BeFalse("§6.1: degenerate fit reports IsOpen = false (closed).");
        fit.PredictNewGenes(10).Should().Be(0.0, "§6.1: degenerate predictor → 0.");
    }

    // Single genome holding zero clusters (empty presence row): still N=1, plus an
    // empty cluster universe. Probes the empty-matrix path crossed with N=1.
    [Test]
    [CancelAfter(20000)]
    public void Single_OneGenomeNoClusters_ReturnsDegenerateFit_NoCrash()
    {
        var universe = Universe(0);
        var rows = new[] { Row("only", universe, new HashSet<string>()) };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 10);
        AssertWellFormedFit(fit);
        fit.Intercept.Should().Be(0.0);
        fit.Alpha.Should().Be(0.0);
        fit.IsOpen.Should().BeFalse();
    }

    // null and empty matrix are the same documented degenerate input (§3.3, §6.1):
    // a degenerate fit, never an exception.
    [Test]
    [CancelAfter(20000)]
    public void Single_NullAndEmptyMatrix_ReturnDegenerateFit_NoThrow()
    {
        var actNull = () => PanGenomeAnalyzer.FitHeapsLaw((IEnumerable<PanGenomeAnalyzer.GenePresenceRow>)null!);
        actNull.Should().NotThrow("§3.3: null matrix ⇒ degenerate fit, never an exception.");

        var fitNull = PanGenomeAnalyzer.FitHeapsLaw((IEnumerable<PanGenomeAnalyzer.GenePresenceRow>)null!);
        AssertWellFormedFit(fitNull);
        fitNull.Intercept.Should().Be(0.0);
        fitNull.Alpha.Should().Be(0.0);

        var fitEmpty = PanGenomeAnalyzer.FitHeapsLaw(Array.Empty<PanGenomeAnalyzer.GenePresenceRow>());
        AssertWellFormedFit(fitEmpty);
        fitEmpty.Intercept.Should().Be(0.0);
        fitEmpty.Alpha.Should().Be(0.0);
    }

    // The dictionary overload with null / empty genomes is the same guard (§3.3).
    [Test]
    [CancelAfter(20000)]
    public void Single_GenomeDictionaryOverload_NullAndEmpty_ReturnDegenerateFit()
    {
        var actNull = () => PanGenomeAnalyzer.FitHeapsLaw(
            (IReadOnlyDictionary<string, IReadOnlyList<(string, string)>>)null!);
        actNull.Should().NotThrow("null genomes ⇒ degenerate fit (§3.3).");

        var fitNull = PanGenomeAnalyzer.FitHeapsLaw(
            (IReadOnlyDictionary<string, IReadOnlyList<(string, string)>>)null!);
        AssertWellFormedFit(fitNull);
        fitNull.Alpha.Should().Be(0.0);

        var empty = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        var fitEmpty = PanGenomeAnalyzer.FitHeapsLaw(empty);
        AssertWellFormedFit(fitEmpty);
        fitEmpty.Alpha.Should().Be(0.0);
    }

    #endregion

    #region BE — Boundary: 2 genomes (exactly one fitted point, the minimum)

    // §6.1 / §2.2: two genomes ⇒ the curve N = 2..G has exactly ONE fitted point
    // (N = 2). This is the MINIMUM data the power fit ever runs on. The fit must be
    // finite and in-bounds, with IsOpen consistent with α — no crash from a
    // single-point pool, no NaN.
    [Test]
    [CancelAfter(20000)]
    public void Two_TwoGenomesOneFittedPoint_FiniteInBoundsFit_NoCrash()
    {
        var universe = Universe(10);
        var rows = new[]
        {
            Row("g1", universe, Clusters(new[] { 0, 1, 2 })),                 // 3 clusters
            Row("g2", universe, Clusters(new[] { 0, 1, 2, 3, 4, 5, 6 })),     // +4 new at N=2
        };

        var act = () => PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);
        act.Should().NotThrow("2 genomes ⇒ a single fitted point, the documented minimum — never a crash.");

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);
        AssertWellFormedFit(fit);

        // With a single point (x=2, y=4) the best fit at α = start (1, then refined)
        // matches the count at N=2 through K; the contract here is finiteness +
        // in-bounds + consistent flag, asserted by AssertWellFormedFit.
        fit.PredictNewGenes(2).Should().BeApproximately(4.0, 1e-2,
            "with a single fitted point the predictor passes through (2, 4).");
    }

    // Two genomes that are DISJOINT (share no clusters): genome 2 contributes ALL of
    // its clusters as new ⇒ a single large positive point. Still finite, in-bounds.
    [Test]
    [CancelAfter(20000)]
    public void Two_TwoDisjointGenomes_AllNewAtSecond_FiniteFit()
    {
        var universe = Universe(12);
        var rows = new[]
        {
            Row("g1", universe, Clusters(new[] { 0, 1, 2, 3, 4, 5 })),
            Row("g2", universe, Clusters(new[] { 6, 7, 8, 9, 10, 11 })),  // +6 new, fully disjoint
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);
        AssertWellFormedFit(fit);
        fit.PredictNewGenes(2).Should().BeApproximately(6.0, 1e-2,
            "the disjoint second genome adds all 6 of its clusters as new at N=2.");
    }

    #endregion

    #region BE — Boundary: identical genomes (flat zero-variance curve ⇒ no log(0)/NaN)

    // IDENTICAL genomes ⇒ every genome after the first adds ZERO new clusters ⇒ y ≡ 0,
    // a flat zero-variance ZERO curve. The fit must NOT blow up: no log(0) (the
    // objective is on y, not log y), no DivideByZero on zero variance, no NaN. The
    // objective J(K,α) = K·sqrt(Σ x^(−2α))/|x| is uniquely minimized at K = 0
    // (= mean count) where α is UNIDENTIFIABLE; micropan's optim, started at (0, 1),
    // cannot improve and leaves α = 1 ⇒ the saturated pan-genome is CLOSED — the
    // biologically correct verdict for genomes that contribute no new genes. (The
    // §6.1 "α = 0, open" row is the constant-POSITIVE curve, tested below.)
    [Test]
    [CancelAfter(20000)]
    public void Identical_AllGenomesIdentical_FlatZeroCurve_KZero_Closed_NoLogZeroNaN()
    {
        var universe = Universe(6);
        var present = Clusters(new[] { 0, 1, 2, 3, 4, 5 });
        var rows = new[]
        {
            Row("g1", universe, new HashSet<string>(present)),
            Row("g2", universe, new HashSet<string>(present)),
            Row("g3", universe, new HashSet<string>(present)),
            Row("g4", universe, new HashSet<string>(present)),
        };

        var act = () => PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 100);
        act.Should().NotThrow("identical genomes ⇒ a flat y≡0 curve; must not throw / log(0) / NaN.");

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 100);
        AssertWellFormedFit(fit);

        fit.Intercept.Should().BeApproximately(0.0, 1e-9,
            "the zero curve is uniquely minimized at K = mean count = 0 (no new genes).");
        fit.IsOpen.Should().BeFalse(
            "identical genomes add zero new genes ⇒ a maximally saturated ⇒ CLOSED pan-genome.");
        fit.PredictNewGenes(7).Should().BeApproximately(0.0, 1e-9,
            "predictor of a zero-K curve is 0 everywhere.");
    }

    // A flat NON-zero curve: each genome after the first adds the SAME positive number
    // of new clusters (a constant > 0). §6.1 "α = 0, K = mean count, open". This
    // exercises the zero-variance-but-positive case (constant y > 0): α = 0, K = the
    // constant, OPEN — and crucially no log(0) (y > 0) and no NaN from flat variance.
    [Test]
    [CancelAfter(20000)]
    public void Identical_FlatPositiveCurve_AlphaZero_KMeanCount_Open()
    {
        // Each genome introduces exactly 3 brand-new clusters ⇒ y = [3,3,3] (flat > 0).
        var universe = Universe(15);
        var rows = new[]
        {
            Row("g1", universe, Clusters(new[] { 0, 1, 2 })),
            Row("g2", universe, Clusters(new[] { 0, 1, 2, 3, 4, 5 })),          // +3
            Row("g3", universe, Clusters(new[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 })), // +3
            Row("g4", universe, Clusters(Enumerable.Range(0, 12))),             // +3
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: 1);
        AssertWellFormedFit(fit);

        fit.Alpha.Should().BeApproximately(0.0, 1e-3,
            "§6.1: a constant positive new-gene curve (3,3,3) ⇒ α = 0.");
        fit.Intercept.Should().BeApproximately(3.0, 1e-2,
            "§6.1: K = the constant new-gene count (mean = 3).");
        fit.IsOpen.Should().BeTrue("α = 0 < 1 ⇒ OPEN (INV-01).");
    }

    #endregion

    #region BE — Boundary: integration via the genome-dictionary overload

    // The clustering overload on identical genomes: real gene sequences, all genomes
    // carrying the SAME families ⇒ flat zero curve after clustering ⇒ K = 0, CLOSED
    // (saturated), no crash from clustering + flat fit composed.
    [Test]
    [CancelAfter(30000)]
    public void Identical_GenomeOverload_SameFamilies_FlatCurve_Closed_NoCrash()
    {
        var families = new[] { Family(1), Family(2), Family(3), Family(4) };
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        for (int g = 1; g <= 4; g++)
        {
            var genes = new List<(string GeneId, string Sequence)>();
            for (int f = 0; f < families.Length; f++)
                genes.Add(Gene($"g{g}_f{f}", families[f])); // same families in every genome
            genomes[$"g{g}"] = genes;
        }

        var act = () => PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 50);
        act.Should().NotThrow("identical genomes through clustering ⇒ flat curve, no crash.");

        var fit = PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 50);
        AssertWellFormedFit(fit);
        fit.Intercept.Should().BeApproximately(0.0, 1e-6,
            "identical genomes ⇒ flat zero new-gene curve ⇒ K = mean count = 0.");
        fit.IsOpen.Should().BeFalse(
            "identical genomes add zero new genes ⇒ saturated ⇒ CLOSED pan-genome.");
    }

    // Single genome through the clustering overload: degenerate guard (§6.1).
    [Test]
    [CancelAfter(20000)]
    public void Single_GenomeOverload_OneGenome_DegenerateFit()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>
        {
            ["only"] = new[] { Gene("only_a", Family(1)), Gene("only_b", Family(2)) },
        };

        var fit = PanGenomeAnalyzer.FitHeapsLaw(genomes, permutations: 10);
        AssertWellFormedFit(fit);
        fit.Intercept.Should().Be(0.0, "§6.1: single genome ⇒ degenerate fit.");
        fit.Alpha.Should().Be(0.0);
        fit.IsOpen.Should().BeFalse();
    }

    #endregion

    #region BE — Boundary: permutations parameter degenerate values

    // permutations ≤ 0 is clamped to ≥ 1 (§3.1 "≥ 1 (clamped)"): must not produce an
    // empty pool / DivideByZero / NaN, and must give the same fit as permutations = 1
    // for a fixed-order (single ordering) curve.
    [Test]
    [CancelAfter(20000)]
    public void Boundary_PermutationsZeroAndNegative_Clamped_NoNaN()
    {
        var universe = Universe(15);
        var rows = new[]
        {
            Row("g1", universe, Clusters(Enumerable.Range(0, 3))),
            Row("g2", universe, Clusters(Enumerable.Range(0, 11))),  // +8
            Row("g3", universe, Clusters(Enumerable.Range(0, 15))),  // +4
        };

        foreach (int perms in new[] { 0, -1, int.MinValue })
        {
            var act = () => PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: perms);
            act.Should().NotThrow($"permutations = {perms} is clamped to ≥ 1 (§3.1), never a crash.");

            var fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: perms);
            AssertWellFormedFit(fit);
            // Clamped to 1 ⇒ natural order only ⇒ matches the §7.1 worked example.
            fit.Alpha.Should().BeApproximately(Math.Log(2.0) / Math.Log(1.5), 1e-4,
                $"permutations = {perms} clamps to 1 ⇒ natural-order fit (the worked example).");
        }
    }

    #endregion

    #region Randomized sweep — never NaN / out-of-bounds / inconsistent flag

    // Broad randomized fuzz over genome counts, cluster universes and random presence
    // matrices including the degenerate boundaries (1 genome, identical rows, all-zero
    // and all-one presence). The fit must ALWAYS be well-formed: finite in-bounds
    // (K, α), IsOpen ⇔ α < 1, non-increasing finite predictor — never a NaN, an
    // out-of-bounds parameter, a hang, or a throw.
    [Test]
    [CancelAfter(60000)]
    public void Fuzz_RandomPresenceMatrices_AlwaysWellFormedFit()
    {
        var rng = new Random(Seed);

        for (int iter = 0; iter < 300; iter++)
        {
            int genomeCount = rng.Next(1, 9);    // include N=1 (degenerate) and small N
            int clusterCount = rng.Next(0, 25);  // include 0 clusters (empty universe)
            var universe = Universe(clusterCount);

            int mode = rng.Next(4);
            var rows = new List<PanGenomeAnalyzer.GenePresenceRow>(genomeCount);
            for (int g = 0; g < genomeCount; g++)
            {
                HashSet<string> present;
                switch (mode)
                {
                    case 0: // fully random presence
                        present = new HashSet<string>(universe.Where(_ => rng.Next(2) == 0));
                        break;
                    case 1: // all absent (zero presence ⇒ flat zero curve)
                        present = new HashSet<string>();
                        break;
                    case 2: // all present (identical full rows ⇒ flat zero curve)
                        present = new HashSet<string>(universe);
                        break;
                    default: // growing presence (each genome a superset of the prior)
                        int upto = clusterCount == 0 ? 0 : Math.Min(clusterCount, (g + 1) * Math.Max(1, clusterCount / genomeCount));
                        present = new HashSet<string>(universe.Take(upto));
                        break;
                }
                rows.Add(Row($"g{g}", universe, present));
            }

            int perms = rng.Next(1, 6);
            PanGenomeAnalyzer.HeapsLawFit fit = default;
            var act = () => fit = PanGenomeAnalyzer.FitHeapsLaw(rows, permutations: perms);
            act.Should().NotThrow(
                $"iter {iter}: genomes={genomeCount}, clusters={clusterCount}, mode={mode}, perms={perms} must never throw.");

            AssertWellFormedFit(fit);

            // Degenerate guard: < 2 genomes ⇒ the documented (0,0,false) fit.
            if (genomeCount < 2)
            {
                fit.Alpha.Should().Be(0.0, $"iter {iter}: < 2 genomes ⇒ degenerate Alpha = 0.");
                fit.Intercept.Should().Be(0.0, $"iter {iter}: < 2 genomes ⇒ degenerate Intercept = 0.");
                fit.IsOpen.Should().BeFalse($"iter {iter}: degenerate fit ⇒ IsOpen false.");
            }
        }
    }

    #endregion

    #endregion
}
