namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Epigenetics area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Epigenetics")]
public class EpigeneticsCombinatorialTests
{
    // Strong CpG island (GC = 1.0, O/E ≈ 2) flanked by AT-rich, CpG-poor filler.
    private static readonly string CpgSeq =
        string.Concat(Enumerable.Repeat("AT", 50)) + string.Concat(Enumerable.Repeat("CG", 200)) + string.Concat(Enumerable.Repeat("AT", 50));

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: EPIGEN-CPG-001 — CpG island detection (Epigenetics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 85.
    // Spec: tests/TestSpecs/EPIGEN-CPG-001.md (canonical FindCpGIslands).
    // Dimensions: windowSize(3) × minOE(3) × minGC(3) × minLen(3). Grid 3⁴ = 81.
    //
    // Model (Gardiner-Garden & Frommer 1987): a CpG island is a region of length ≥ minLength with
    // GC content ≥ minGc and observed/expected CpG ratio ≥ minCpGRatio. The detection window equals
    // minLength internally, so windowSize is realised as the searched-substring length.
    //
    // The combinatorial point: the searched window and the three criteria jointly constrain the
    // output — every reported island simultaneously satisfies the GC, O/E and length thresholds,
    // for all 81 parameter combinations, and lies within the searched window.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void EpigenCpg_ReportedIslandsSatisfyAllCriteria(
        [Values(150, 350, 600)] int windowSize,
        [Values(0.4, 0.6, 1.0)] double minOE,
        [Values(0.4, 0.6, 0.8)] double minGC,
        [Values(100, 200, 300)] int minLen)
    {
        string searched = CpgSeq[..windowSize];
        var islands = EpigeneticsAnalyzer.FindCpGIslands(searched, minLen, minGC, minOE).ToList();

        islands.Should().OnlyContain(isl => isl.GcContent >= minGC, "island GC clears the threshold");
        islands.Should().OnlyContain(isl => isl.CpGRatio >= minOE, "island O/E clears the threshold");
        islands.Should().OnlyContain(isl => isl.End - isl.Start >= minLen, "island length clears the minimum");
        islands.Should().OnlyContain(isl => isl.Start >= 0 && isl.End <= searched.Length, "coordinates within the searched window");
    }

    /// <summary>
    /// Interaction witness: a strong CpG island is detected under permissive criteria, while an
    /// unreachable O/E or length floor removes it.
    /// </summary>
    [Test]
    public void EpigenCpg_StrongIsland_DetectedAndGated()
    {
        EpigeneticsAnalyzer.FindCpGIslands(CpgSeq, 200, 0.4, 0.4).Should().NotBeEmpty("a CG-rich island passes permissive criteria");
        EpigeneticsAnalyzer.FindCpGIslands(CpgSeq, 200, 0.4, 5.0).Should().BeEmpty("no region reaches an O/E of 5");
        EpigeneticsAnalyzer.FindCpGIslands(CpgSeq, 5000, 0.4, 0.4).Should().BeEmpty("no 5000-bp island fits in the sequence");
    }

    /// <summary>
    /// Interaction witness: AT-rich, CpG-poor sequence has no islands regardless of length.
    /// </summary>
    [Test]
    public void EpigenCpg_AtRichSequence_NoIslands()
    {
        EpigeneticsAnalyzer.FindCpGIslands(string.Concat(Enumerable.Repeat("AT", 200)), 200, 0.5, 0.6)
            .Should().BeEmpty("an AT-rich tract is not a CpG island");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: EPIGEN-AGE-001 — Epigenetic-clock age prediction (Epigenetics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 181.
    // Spec: tests/TestSpecs/EPIGEN-AGE-001.md (canonical CalculateEpigeneticAge / HorvathAntiTransform).
    // ADVANCED §10.
    // Dimensions: clock(2) × nSites(3). Grid 2×3 = 6 (full, exhaustive ⊇ pairwise).
    //
    // Model (Horvath 2013): age = anti.trafo(intercept + Σ coefᵢ·βᵢ) over the clock CpGs, where
    // anti.trafo(x) = (1+20)·eˣ − 1 for x < 0, else (1+20)·x + 20 (adult.age = 20). Only CpGs present
    // in BOTH the methylation profile and the coefficient table contribute to the linear predictor.
    //
    // Axis mapping (documented): clock → two coefficient tables — one driving the linear predictor
    // negative (the exponential anti-trafo branch), one non-negative (the linear branch); nSites →
    // the number of clock CpGs (3/6/9). The combinatorial point: across both clocks and every site
    // count the production age equals an independent recomputation of the linear predictor + anti-trafo.
    // ═══════════════════════════════════════════════════════════════════════

    public enum AgeClock { ExpBranch, LinearBranch }

    private static (Dictionary<string, double> Meth, Dictionary<string, double> Coef, double Intercept) BuildClock(AgeClock clock, int nSites)
    {
        var meth = new Dictionary<string, double>();
        var coef = new Dictionary<string, double>();
        double perSite = clock == AgeClock.ExpBranch ? -0.2 : 0.6;
        for (int i = 0; i < nSites; i++)
        {
            string cg = $"cg{i:D3}";
            meth[cg] = (i % 5 + 1) / 10.0; // 0.1..0.5
            coef[cg] = perSite;
        }
        double intercept = clock == AgeClock.ExpBranch ? -0.5 : 1.0;
        return (meth, coef, intercept);
    }

    [Test, Combinatorial]
    public void EpigenAge_MatchesLinearPredictorAndAntiTrafo_AcrossClockAndSites(
        [Values(AgeClock.ExpBranch, AgeClock.LinearBranch)] AgeClock clock,
        [Values(3, 6, 9)] int nSites)
    {
        var (meth, coef, intercept) = BuildClock(clock, nSites);

        // Independent ground truth: linear predictor + Horvath anti-transform.
        double lp = intercept;
        foreach (var (cg, b) in meth)
            if (coef.TryGetValue(cg, out double c)) lp += c * b;
        double expected = lp < 0 ? 21.0 * Math.Exp(lp) - 1.0 : 21.0 * lp + 20.0;

        // The chosen clock exercises the intended anti-trafo branch.
        (lp < 0).Should().Be(clock == AgeClock.ExpBranch, "the clock drives the predictor into the intended branch");

        EpigeneticsAnalyzer.CalculateEpigeneticAge(meth, coef, intercept)
            .Should().BeApproximately(expected, 1e-9, "age = anti.trafo(intercept + Σ coef·β)");
    }

    /// <summary>
    /// Interaction witnesses — the anti-transform is continuous & monotone at the branch boundary
    /// (f(0)=20), and CpGs absent from the coefficient table contribute nothing to the predictor.
    /// </summary>
    [Test]
    public void EpigenAge_AntiTrafoBranchesAndCoefficientLookup()
    {
        EpigeneticsAnalyzer.HorvathAntiTransform(0.0).Should().BeApproximately(20.0, 1e-12, "f(0) = adult.age");
        EpigeneticsAnalyzer.HorvathAntiTransform(-1.0).Should().BeApproximately(21.0 * Math.Exp(-1.0) - 1.0, 1e-12, "exp branch");
        EpigeneticsAnalyzer.HorvathAntiTransform(1.0).Should().BeApproximately(41.0, 1e-12, "linear branch: 21·1+20");
        EpigeneticsAnalyzer.HorvathAntiTransform(-0.5).Should().BeLessThan(EpigeneticsAnalyzer.HorvathAntiTransform(0.5), "monotone increasing");

        var coef = new Dictionary<string, double> { ["cgA"] = 1.0 };
        var withExtra = new Dictionary<string, double> { ["cgA"] = 0.5, ["cgUnknown"] = 0.9 };
        var without = new Dictionary<string, double> { ["cgA"] = 0.5 };
        EpigeneticsAnalyzer.CalculateEpigeneticAge(withExtra, coef, 0.0)
            .Should().Be(EpigeneticsAnalyzer.CalculateEpigeneticAge(without, coef, 0.0),
                "a CpG not in the coefficient table does not change the age");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: EPIGEN-CHROM-001 — Chromatin-state prediction from histone marks (Epigenetics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 183.
    // Spec: tests/TestSpecs/EPIGEN-CHROM-001.md (canonical PredictChromatinState). ADVANCED §10.
    // Dimensions: nMarks(3) × windowSize(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Ernst & Kellis ChromHMM; Roadmap Epigenomics): each mark is binarized at a presence
    // threshold, and the combinatorial signature maps to a chromatin state by a fixed priority
    // (bivalent K4me3+K27me3 → BivalentPromoter; K4me1+K27me3 → BivalentEnhancer; K4me3 → ActivePromoter;
    // K4me1 → Active/WeakEnhancer by K27ac; K36me3 → Transcribed; K27me3 → Repressed; K9me3 → Het; else Low).
    //
    // Axis mapping (documented — PredictChromatinState has no windowSize knob): nMarks → a signature
    // scenario with 1/2/3 active marks; windowSize → the presenceThreshold (a mark held at 0.5 flips
    // present↔absent across thresholds). The combinatorial point: the state equals an independent
    // re-implementation of the Roadmap priority at every (scenario, threshold) cell.
    // ═══════════════════════════════════════════════════════════════════════

    private static EpigeneticsAnalyzer.ChromatinState ExpectedState(
        double k4me3, double k4me1, double k27ac, double k36me3, double k27me3, double k9me3, double thr)
    {
        bool a = k4me3 >= thr, b = k4me1 >= thr, c = k27ac >= thr, d = k36me3 >= thr, e = k27me3 >= thr, f = k9me3 >= thr;
        if (a && e) return EpigeneticsAnalyzer.ChromatinState.BivalentPromoter;
        if (b && e) return EpigeneticsAnalyzer.ChromatinState.BivalentEnhancer;
        if (a) return EpigeneticsAnalyzer.ChromatinState.ActivePromoter;
        if (b) return c ? EpigeneticsAnalyzer.ChromatinState.ActiveEnhancer : EpigeneticsAnalyzer.ChromatinState.WeakEnhancer;
        if (d) return EpigeneticsAnalyzer.ChromatinState.Transcribed;
        if (e) return EpigeneticsAnalyzer.ChromatinState.Repressed;
        if (f) return EpigeneticsAnalyzer.ChromatinState.Heterochromatin;
        return EpigeneticsAnalyzer.ChromatinState.LowSignal;
    }

    [Test, Combinatorial]
    public void EpigenChrom_SignatureMapsToState_AcrossScenarioAndThreshold(
        [Values(1, 2, 3)] int nMarks,
        [Values(0.3, 0.5, 0.7)] double presenceThreshold)
    {
        // (k4me3, k4me1, k27ac, k36me3, k27me3, k9me3); one mark held at 0.5 to be threshold-sensitive.
        double[] s = nMarks switch
        {
            1 => new[] { 1.0, 0, 0, 0, 0, 0 },        // → ActivePromoter
            2 => new[] { 0.0, 1.0, 0.5, 0, 0, 0 },    // K4me1 + (K27ac at 0.5) → Active/WeakEnhancer by threshold
            _ => new[] { 1.0, 0, 0, 1.0, 1.0, 0 },    // K4me3 + K27me3 (+K36me3) → BivalentPromoter
        };

        var state = EpigeneticsAnalyzer.PredictChromatinState(s[0], s[1], s[2], s[3], s[4], s[5], presenceThreshold);
        state.Should().Be(ExpectedState(s[0], s[1], s[2], s[3], s[4], s[5], presenceThreshold),
            "the binarized signature maps to the Roadmap state");
    }

    /// <summary>
    /// Interaction witnesses — bivalent priority (an active mark co-occurring with K27me3 yields a
    /// bivalent state, not the active one) and threshold-driven enhancer flip (K27ac at 0.5).
    /// </summary>
    [Test]
    public void EpigenChrom_BivalentPriority_AndThresholdFlip()
    {
        EpigeneticsAnalyzer.PredictChromatinState(1.0, 0, 0, 0, 1.0, 0, 0.5)
            .Should().Be(EpigeneticsAnalyzer.ChromatinState.BivalentPromoter, "K4me3+K27me3 outranks ActivePromoter");
        EpigeneticsAnalyzer.PredictChromatinState(1.0, 0, 0, 0, 0, 0, 0.5)
            .Should().Be(EpigeneticsAnalyzer.ChromatinState.ActivePromoter, "K4me3 alone is an active promoter");

        EpigeneticsAnalyzer.PredictChromatinState(0, 1.0, 0.5, 0, 0, 0, 0.3)
            .Should().Be(EpigeneticsAnalyzer.ChromatinState.ActiveEnhancer, "K27ac present at threshold 0.3");
        EpigeneticsAnalyzer.PredictChromatinState(0, 1.0, 0.5, 0, 0, 0, 0.7)
            .Should().Be(EpigeneticsAnalyzer.ChromatinState.WeakEnhancer, "K27ac absent at threshold 0.7");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: EPIGEN-DMR-001 — Differentially methylated regions (Epigenetics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 184.
    // Spec: tests/TestSpecs/EPIGEN-DMR-001.md (canonical FindDMRs). ADVANCED §10.
    // Dimensions: threshold(3) × minSites(3) × nSamples(2). Grid 3×3×2 = 18 (full, exhaustive).
    //
    // Model (Akalin 2012 methylKit): tile the genome into windows; a window is a DMR when it has
    // ≥ minCpGCount CpGs and |mean methylation difference| > minDifference (STRICT), with a Fisher's
    // exact p-value over the pooled C/T counts; hyper vs hypo by the sign of the difference.
    //
    // Axis mapping (documented — FindDMRs compares exactly two samples): threshold → minDifference;
    // minSites → minCpGCount; nSamples → the differential DIRECTION (Hyper/Hypo). Engineered construct:
    // one window of 3 CpGs with |diff| = 0.3. The combinatorial point: the window is called exactly
    // when minDifference < 0.3 AND minCpGCount ≤ 3, with the correct hyper/hypo annotation and a valid
    // p-value — verified per cell.
    // ═══════════════════════════════════════════════════════════════════════

    public enum DmrDirection { Hyper, Hypo }

    private static EpigeneticsAnalyzer.MethylationSite Site(int pos, double level) =>
        new(pos, EpigeneticsAnalyzer.MethylationType.CpG, "CpG", level, 100);

    [Test, Combinatorial]
    public void EpigenDmr_WindowCalledExactlyWhenBothFloorsMet_AcrossThresholdSitesDirection(
        [Values(0.1, 0.25, 0.4)] double threshold,
        [Values(2, 3, 4)] int minSites,
        [Values(DmrDirection.Hyper, DmrDirection.Hypo)] DmrDirection direction)
    {
        var positions = new[] { 0, 100, 200 }; // 3 CpGs inside one default 1000-bp window
        double delta = direction == DmrDirection.Hyper ? 0.3 : -0.3;

        var sample1 = positions.Select(p => Site(p, 0.5)).ToList();
        var sample2 = positions.Select(p => Site(p, 0.5 + delta)).ToList();

        var dmrs = EpigeneticsAnalyzer.FindDMRs(sample1, sample2, windowSize: 1000, minDifference: threshold, minCpGCount: minSites).ToList();

        bool expectCalled = threshold < 0.3 - 1e-9 && minSites <= 3;
        (dmrs.Count != 0).Should().Be(expectCalled, "a window is a DMR iff |diff| > minDifference AND CpGs ≥ minCpGCount");

        if (expectCalled)
        {
            var dmr = dmrs.Single();
            dmr.CpGCount.Should().Be(3);
            Math.Abs(dmr.MeanDifference).Should().BeApproximately(0.3, 1e-9);
            dmr.PValue.Should().BeInRange(0.0, 1.0);
            dmr.Annotation.Should().Be(direction == DmrDirection.Hyper ? "Hypermethylated" : "Hypomethylated",
                "the sign of the difference sets hyper/hypo");
        }
    }

    /// <summary>
    /// Interaction witness — each floor independently gates the call: too-high minDifference or
    /// too-high minCpGCount suppresses the DMR.
    /// </summary>
    [Test]
    public void EpigenDmr_EachFloorGatesTheCall()
    {
        var positions = new[] { 0, 100, 200 };
        var s1 = positions.Select(p => Site(p, 0.5)).ToList();
        var s2 = positions.Select(p => Site(p, 0.8)).ToList(); // |diff| = 0.3, 3 CpGs

        EpigeneticsAnalyzer.FindDMRs(s1, s2, 1000, 0.1, 3).Should().ContainSingle("permissive floors call the window");
        EpigeneticsAnalyzer.FindDMRs(s1, s2, 1000, 0.4, 3).Should().BeEmpty("|diff| 0.3 ≤ minDifference 0.4");
        EpigeneticsAnalyzer.FindDMRs(s1, s2, 1000, 0.1, 4).Should().BeEmpty("only 3 CpGs < minCpGCount 4");
    }
}
