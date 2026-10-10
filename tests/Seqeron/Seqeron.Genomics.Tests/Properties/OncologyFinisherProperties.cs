using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// FIN-B24 heavy tier — property-based tests for the behaviour added by the B24 finisher (F44–F64). Oracles are
/// independent restatements of the reference definitions (never routed through production helpers):
/// <list type="bullet">
/// <item>Ckmeans.1d.dp (F46/F47): the returned partition attains the optimal within-cluster sum of squares (an
/// independent O(k·n²) dynamic program over the sorted values), and the BIC-selected k lies in the adjusted bounds and
/// equals the fixed-k clustering for that k; constant input gives one cluster.</item>
/// <item>LICHeE cluster summaries (F42/F43): with zero SDs and robust clusters the per-cluster margins collapse to ε,
/// i.e. the result equals <see cref="OncologyAnalyzer.ReconstructPhylogeny"/>; a bijective id relabelling maps the tree.</item>
/// <item>maftools <c>math.score</c> filters (F56): raising the VAF cutoff can only drop VAFs, so once a sample is skipped
/// (null) every higher cutoff skips it; cutoff 0 with n_min 1 is the unfiltered MATH.</item>
/// <item>Shannon clonal diversity (F57): proportion-only (permutation/scale invariant), 0 ≤ H ≤ ln k.</item>
/// <item>CNAqc peak QC (F34/F62/F64): expected peaks are <c>expected_vaf_fun</c> m·π/(2(1 − π) + π·n) for the QC'd π, so
/// a purity shift moves them exactly along that curve; the verdicts depend on the multiset of VAFs (mutation order is
/// irrelevant); the subclonal R seed only names mutations; BMix auto mode = on iff read counts are supplied.</item>
/// </list>
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Oncology")]
public class OncologyFinisherProperties
{
    // -------------------------------------------------------------------------
    // Ckmeans.1d.dp: optimum (SMAWK fill, F46) and BIC k (F47)
    // -------------------------------------------------------------------------

    private static Arbitrary<(double[] Values, int K, int Seed)> CcfArbitrary() =>
        (from n in Gen.Choose(1, 40)
         from raw in Gen.Choose(0, 1000).ArrayOf(n)
         from k in Gen.Choose(1, Math.Max(1, Math.Min(6, n)))
         from seed in Gen.Choose(0, int.MaxValue)
         select (raw.Select(v => v / 1000.0).ToArray(), k, seed)).ToArbitrary();

    private static double OptimalWcss(double[] values, int k)
    {
        double[] x = values.OrderBy(v => v).ToArray();
        int n = x.Length;
        double Cost(int i, int j) // inclusive [i, j]
        {
            double mean = 0;
            for (int t = i; t <= j; t++) mean += x[t];
            mean /= j - i + 1;
            double s = 0;
            for (int t = i; t <= j; t++) s += (x[t] - mean) * (x[t] - mean);
            return s;
        }

        var d = new double[k + 1, n + 1];
        for (int j = 1; j <= n; j++) d[0, j] = double.PositiveInfinity;
        for (int c = 1; c <= k; c++)
        {
            for (int j = 0; j <= n; j++)
            {
                if (j < c) { d[c, j] = j == 0 ? 0 : double.PositiveInfinity; continue; }
                double best = double.PositiveInfinity;
                for (int i = c; i <= j; i++) // last cluster = x[i-1 .. j-1]
                {
                    best = Math.Min(best, d[c - 1, i - 1] + Cost(i - 1, j - 1));
                }

                d[c, j] = best;
            }
        }

        return d[k, n];
    }

    private static double Wcss(double[] values, OncologyAnalyzer.CcfClustering c) =>
        values.Select((v, i) => (v - c.Centroids[c.Assignments[i]]) * (v - c.Centroids[c.Assignments[i]])).Sum();

    /// <summary>
    /// F46 (SMAWK row fill): the returned clustering attains the brute-force optimal WCSS for k* = min(k, #distinct)
    /// clusters (only the choice among equal-WCSS optima may differ from the quadratic fill), centroids are ascending
    /// cluster means and the clonal index is the last centroid.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClusterCcfValues_AttainsBruteForceOptimalWcss()
    {
        return Prop.ForAll(CcfArbitrary(), t =>
        {
            var c = OncologyAnalyzer.ClusterCcfValues(t.Values, t.K);
            int kEff = Math.Min(t.K, t.Values.Distinct().Count());
            double wcss = Wcss(t.Values, c);
            double optimum = OptimalWcss(t.Values, kEff);
            bool means = Enumerable.Range(0, c.Centroids.Count).All(j =>
                Math.Abs(c.Centroids[j] - t.Values.Where((_, i) => c.Assignments[i] == j).Average()) <= 1e-12);
            bool ascending = c.Centroids.Zip(c.Centroids.Skip(1), (a, b) => a < b).All(ok => ok);
            return (c.Centroids.Count == kEff && Math.Abs(wcss - optimum) <= 1e-9 * (1 + optimum) && means && ascending
                    && c.ClonalClusterIndex == c.Centroids.Count - 1)
                .Label($"k*={kEff}, centroids={c.Centroids.Count}, wcss={wcss}, optimum={optimum}");
        });
    }

    /// <summary>
    /// F47 (Ckmeans.1d.dp <c>select_levels</c> BIC): the selected k lies in [min(k_min, d), min(k_max, d)] (d = #distinct,
    /// R's bound adjustment) and the clustering is the fixed-k clustering for that k; all-equal input gives one cluster.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ClusterCcfValues_BicSelection_WithinBounds_EqualsFixedK()
    {
        return Prop.ForAll(CcfArbitrary(), t =>
        {
            int kMin = Math.Max(1, t.K - 1), kMax = t.K + 2;
            int distinct = t.Values.Distinct().Count();
            var auto = OncologyAnalyzer.ClusterCcfValues(t.Values, kMin, kMax);
            int k = auto.Centroids.Count;
            int lo = distinct < kMin ? distinct : kMin, hi = Math.Min(kMax, distinct);
            var fixedK = OncologyAnalyzer.ClusterCcfValues(t.Values, k);
            var constant = OncologyAnalyzer.ClusterCcfValues(Enumerable.Repeat(t.Values[0], t.Values.Length).ToArray(), kMin, kMax);
            return (k >= lo && k <= hi && auto.Assignments.SequenceEqual(fixedK.Assignments)
                    && auto.Centroids.Zip(fixedK.Centroids, (a, b) => Math.Abs(a - b) <= 1e-12).All(ok => ok)
                    && constant.Centroids.Count == 1)
                .Label($"k={k} ∉ [{lo}, {hi}] or differs from fixed-k");
        });
    }

    // -------------------------------------------------------------------------
    // LICHeE cluster summaries (F42, F43)
    // -------------------------------------------------------------------------

    private static Arbitrary<(double[][] Ccf, int[] IdMap)> PhyloArbitrary() =>
        (from k in Gen.Choose(1, 2)
         from n in Gen.Choose(2, 5)
         from raw in Gen.Choose(0, 10).ArrayOf(k).ArrayOf(n)
         from ids in Gen.Shuffle(Enumerable.Range(10, n).ToArray())
         select (raw.Select((r, i) => r.Select(v => i == 0 ? 1.0 : v / 10.0).ToArray()).ToArray(), ids.ToArray())).ToArbitrary();

    /// <summary>
    /// F42: with SD = 0 the margin max(ε, se_u + se_v) is ε, and with robust clusters F43's <c>fixNetwork</c> drops
    /// nothing, so the summary path equals <see cref="OncologyAnalyzer.TryReconstructPhylogeny"/> (success, edges, root,
    /// error score); a bijective id relabelling maps the summary-path tree edge-for-edge.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ReconstructPhylogenyFromClusterSummaries_ZeroSd_EqualsStaticPath_AndRelabels()
    {
        return Prop.ForAll(PhyloArbitrary(), t =>
        {
            int n = t.Ccf.Length;
            var clusters = Enumerable.Range(0, n).Select(i => new OncologyAnalyzer.CcfCluster(i + 1, t.Ccf[i])).ToArray();
            var summaries = clusters.Select(c => new OncologyAnalyzer.CcfClusterSummary(
                c.Id, c.CcfPerSample, new double[c.CcfPerSample.Count], 10)).ToArray();
            bool okA = OncologyAnalyzer.TryReconstructPhylogeny(clusters, out var a);
            bool okB = OncologyAnalyzer.TryReconstructPhylogenyFromClusterSummaries(summaries, out var b);
            if (okA != okB)
            {
                return false.Label($"static {okA} vs summaries {okB}");
            }

            if (!okA)
            {
                return true.Label("infeasible on both paths");
            }

            var map = Enumerable.Range(0, n).ToDictionary(i => i + 1, i => t.IdMap[i]);
            var relabelled = summaries.Select(s => s with { Id = map[s.Id] }).ToArray();
            bool okC = OncologyAnalyzer.TryReconstructPhylogenyFromClusterSummaries(relabelled, out var c);
            var mappedEdges = b.Edges.Select(e => (e.ParentId == b.RootId ? -1 : map[e.ParentId], map[e.ChildId])).OrderBy(e => e).ToList();
            var relEdges = c.Edges.Select(e => (e.ParentId == c.RootId ? -1 : e.ParentId, e.ChildId)).OrderBy(e => e).ToList();
            return (a.Edges.SequenceEqual(b.Edges) && a.RootId == b.RootId && a.ErrorScore == b.ErrorScore
                    && okC && mappedEdges.SequenceEqual(relEdges))
                .Label($"edges static [{string.Join(",", a.Edges)}] vs summary [{string.Join(",", b.Edges)}]; relabel ok={okC}");
        });
    }

    // -------------------------------------------------------------------------
    // maftools math.score filters (F56), Shannon clonal diversity (F57)
    // -------------------------------------------------------------------------

    /// <summary>
    /// F56: the retained set {VAF ≥ cutoff} shrinks as the cutoff grows, so null (fewer than n_min retained) at one cutoff
    /// implies null at every higher cutoff; a non-null score equals the unfiltered MATH of the retained VAFs; cutoff 0 with
    /// n_min = 1 equals the unfiltered MATH.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateIthWithFilters_CutoffMonotone_AndMatchesUnfilteredOnRetained()
    {
        var arb = (from n in Gen.Choose(1, 30)
                   from raw in Gen.Choose(1, 1000).ArrayOf(n)
                   from c1 in Gen.Choose(0, 600)
                   from dc in Gen.Choose(0, 400)
                   from minM in Gen.Choose(1, 8)
                   select (raw.Select(v => v / 1000.0).ToArray(), c1 / 1000.0, (c1 + dc) / 1000.0, minM)).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            (double[] vafs, double lowCut, double highCut, int minM) = t;
            double? low = OncologyAnalyzer.CalculateITH(vafs, lowCut, minM);
            double? high = OncologyAnalyzer.CalculateITH(vafs, highCut, minM);
            var retained = vafs.Where(v => !(v < highCut)).ToArray();
            bool highMatches = high is null ? retained.Length < minM : high == OncologyAnalyzer.CalculateITH(retained);
            bool unfiltered = OncologyAnalyzer.CalculateITH(vafs, 0.0, 1) == OncologyAnalyzer.CalculateITH(vafs);
            return ((low is not null || high is null) && highMatches && unfiltered)
                .Label($"low={low}, high={high}, retained={retained.Length}, minM={minM}");
        });
    }

    /// <summary>
    /// F57: H over clone frequencies is proportion-only — permutation-invariant and invariant to halving every
    /// frequency (exact) — bounded by 0 ≤ H ≤ ln(#non-zero clones), and equals the canonical Shannon of weights.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateCloneShannonDiversity_ProportionOnly_Bounded()
    {
        var arb = (from n in Gen.Choose(1, 12)
                   from raw in Gen.Choose(0, 1000).ArrayOf(n)
                   from seed in Gen.Choose(0, int.MaxValue)
                   where raw.Any(v => v > 0)
                   select (raw.Select(v => v / 1000.0).ToArray(), seed)).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            var rng = new Random(t.Item2);
            double h = OncologyAnalyzer.CalculateCloneShannonDiversity(t.Item1);
            double perm = OncologyAnalyzer.CalculateCloneShannonDiversity(t.Item1.OrderBy(_ => rng.Next()).ToArray());
            double half = OncologyAnalyzer.CalculateCloneShannonDiversity(t.Item1.Select(f => f / 2).ToArray());
            int k = t.Item1.Count(f => f > 0);
            return (Math.Abs(perm - h) <= 1e-12 && half == h && h >= 0 && h <= Math.Log(k) + 1e-12
                    && h == StatisticsHelper.ShannonIndexOfWeights(t.Item1))
                .Label($"H={h}, perm={perm}, half={half}, ln k={Math.Log(k)}");
        });
    }

    // -------------------------------------------------------------------------
    // CNAqc peak QC (F34, F62, F64)
    // -------------------------------------------------------------------------

    private static double ExpectedVaf(int m, double purity, int copies) => m * purity / (2 * (1 - purity) + purity * copies);

    private static Gen<(double Purity, double QcPurity, (int Major, int Minor)[] Karyotypes, int Seed)> PeakScenarioGen(
        (int Major, int Minor)[] pool) =>
        from purityPct in Gen.Choose(35, 95)
        from shift in Gen.Choose(-15, 15)
        from count in Gen.Choose(1, 2)
        from picks in Gen.Shuffle(pool)
        from seed in Gen.Choose(0, int.MaxValue)
        select (purityPct / 100.0, Math.Clamp(purityPct + shift, 20, 100) / 100.0, picks.Take(count).ToArray(), seed);

    private static List<OncologyAnalyzer.PurityPeakMutation> SimulateMutations(
        Random rng, double purity, (int Major, int Minor)[] karyotypes, int perKaryotype, bool withReads)
    {
        var list = new List<OncologyAnalyzer.PurityPeakMutation>();
        foreach (var (major, minor) in karyotypes)
        {
            for (int i = 0; i < perKaryotype; i++)
            {
                int m = rng.Next(2) == 0 ? 1 : Math.Max(1, major);
                double mean = ExpectedVaf(m, purity, major + minor);
                int depth = 120;
                int alt = 0;
                for (int r = 0; r < depth; r++)
                {
                    if (rng.NextDouble() < mean) alt++;
                }

                alt = Math.Max(1, alt);
                var mutation = new OncologyAnalyzer.PurityPeakMutation((double)alt / depth, major, minor);
                list.Add(withReads ? mutation with { AlternateReads = alt, Depth = depth } : mutation);
            }
        }

        return list;
    }

    private static readonly (int, int)[] SimplePool = { (1, 1), (2, 0), (2, 1), (2, 2), (1, 0) };

    /// <summary>
    /// F34: for any QC'd purity q the expected peaks are CNAqc <c>expected_vaf_fun</c> m·q/(2(1 − q) + q·(Major + minor))
    /// for m ∈ {1, Major} — a purity shift moves every expected peak along that curve — and the analysis depends only on
    /// the multiset of mutations (score, verdict and matched flags are permutation-invariant; KDE peaks to the 1e−12
    /// rounding of the density sums).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 40)]
    public Property AnalyzePurityPeaks_ExpectedPeaksFollowPurity_PermutationInvariant()
    {
        return Prop.ForAll(PeakScenarioGen(SimplePool).ToArbitrary(), s =>
        {
            var rng = new Random(s.Seed);
            var mutations = SimulateMutations(rng, s.Purity, s.Karyotypes, 120, withReads: false);
            var a = OncologyAnalyzer.AnalyzePurityPeaks(mutations, s.QcPurity);
            var b = OncologyAnalyzer.AnalyzePurityPeaks(mutations.OrderBy(_ => rng.Next()).ToList(), s.QcPurity);

            bool expected = a.Matches.All(m =>
                (m.Multiplicity == 1 || m.Multiplicity == m.MajorCopyNumber)
                && Math.Abs(m.ExpectedPeak - ExpectedVaf(m.Multiplicity, s.QcPurity, m.MajorCopyNumber + m.MinorCopyNumber)) <= 1e-15);
            bool sameScore = double.IsNaN(a.Score) ? double.IsNaN(b.Score) : Math.Abs(a.Score - b.Score) <= 1e-12;
            bool invariant = sameScore && a.Pass == b.Pass
                             && a.Matches.Select(m => (m.MajorCopyNumber, m.MinorCopyNumber, m.Multiplicity, m.Matched))
                                 .SequenceEqual(b.Matches.Select(m => (m.MajorCopyNumber, m.MinorCopyNumber, m.Multiplicity, m.Matched)));
            return (a.Matches.Count > 0 && expected && invariant)
                .Label($"π={s.Purity}, q={s.QcPurity}, karyotypes=[{string.Join(",", s.Karyotypes)}], score {a.Score} vs {b.Score}, pass {a.Pass} vs {b.Pass}");
        });
    }

    private static readonly (int, int)[] ComplexPool = { (3, 1), (3, 0), (4, 2), (3, 2), (4, 1) };

    /// <summary>
    /// F62 (<c>analyze_peaks_general</c>): for every analysed complex karyotype the expected peaks are m = 1..Major with
    /// m·q/(2(1 − q) + q·(Major + minor)); matched + mismatched = #expected; the matched proportion and PASS
    /// (prop ≥ 0.5) are invariant to the mutation order.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 40)]
    public Property AnalyzeComplexKaryotypePeaks_ExpectedPeaks_Bookkeeping_PermutationInvariant()
    {
        return Prop.ForAll(PeakScenarioGen(ComplexPool).ToArbitrary(), s =>
        {
            var rng = new Random(s.Seed);
            var mutations = SimulateMutations(rng, s.Purity, s.Karyotypes, 130, withReads: false);
            var a = OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(mutations, s.QcPurity);
            var b = OncologyAnalyzer.AnalyzeComplexKaryotypePeaks(mutations.OrderBy(_ => rng.Next()).ToList(), s.QcPurity);

            bool ok = a.Ran && a.Karyotypes.Count == s.Karyotypes.Length;
            foreach (var k in a.Karyotypes)
            {
                var exp = k.ExpectedPeaks;
                ok &= exp.Select(e => e.Multiplicity).SequenceEqual(Enumerable.Range(1, k.MajorCopyNumber))
                      && exp.All(e => Math.Abs(e.ExpectedPeak - ExpectedVaf(e.Multiplicity, s.QcPurity, k.MajorCopyNumber + k.MinorCopyNumber)) <= 1e-15)
                      && k.MatchedPeaks + k.MismatchedPeaks == exp.Count
                      && k.MatchedPeaks == exp.Count(e => e.Matched)
                      && k.Pass == (k.MatchedProportion >= 0.5);
            }

            ok &= a.Karyotypes.Select(k => (k.MajorCopyNumber, k.MinorCopyNumber, k.MatchedProportion, k.Pass))
                .SequenceEqual(b.Karyotypes.Select(k => (k.MajorCopyNumber, k.MinorCopyNumber, k.MatchedProportion, k.Pass)));
            return ok.Label($"π={s.Purity}, q={s.QcPurity}, karyotypes=[{string.Join(",", s.Karyotypes)}]");
        });
    }

    /// <summary>
    /// F62 (<c>analyze_peaks_subclonal</c>): the R seed only draws the 8-letter mutation identifiers, so expected-peak
    /// values, matches, rankings and the best models do not depend on it; permuting a segment's VAFs leaves them unchanged.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 30)]
    public Property AnalyzeSubclonalPurityPeaks_SeedAndVafOrder_DoNotChangeDecisions()
    {
        var pairs = new[] { ((2, 1), (1, 1)), ((2, 0), (1, 1)), ((1, 0), (1, 1)), ((2, 1), (2, 0)), ((1, 1), (2, 2)) };
        var arb = (from purityPct in Gen.Choose(40, 95)
                   from ccfPct in Gen.Choose(20, 80)
                   from pair in Gen.Elements(pairs)
                   from seed in Gen.Choose(0, int.MaxValue)
                   from rSeed in Gen.Choose(0, 100_000)
                   select (purityPct / 100.0, ccfPct / 100.0, pair, seed, rSeed)).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            (double purity, double ccf, var pair, int seed, int rSeed) = t;
            var rng = new Random(seed);
            ((int ma, int mi), (int ma2, int mi2)) = pair;
            double tumourCopies = ccf * (ma + mi) + (1 - ccf) * (ma2 + mi2);
            var vafs = Enumerable.Range(0, 150).Select(_ =>
            {
                double mean = (rng.Next(2) == 0 ? ccf : 1.0) * purity / (2 * (1 - purity) + purity * tumourCopies);
                return Math.Clamp(mean + (rng.NextDouble() - 0.5) * 0.06, 0.01, 0.99);
            }).ToList();
            var segment = new OncologyAnalyzer.SubclonalPeakSegment("1", 1, 1_000_000, ma, mi, ma2, mi2, ccf, vafs);
            var shuffled = segment with { Vafs = vafs.OrderBy(_ => rng.Next()).ToList() };

            string Run(OncologyAnalyzer.SubclonalPeakSegment seg, int rs)
            {
                try
                {
                    var r = OncologyAnalyzer.AnalyzeSubclonalPurityPeaks(new[] { seg }, purity,
                        OncologyAnalyzer.SubclonalPeakOptions.Default with { Seed = rs });
                    return string.Join(";", r.Select(x =>
                        string.Join(",", x.ExpectedPeaks.Select(e => (e.ModelId, e.ExpectedPeak.ToString("R"), e.Matched)).OrderBy(e => e.ToString(), StringComparer.Ordinal))
                        + "|" + string.Join(",", x.Rankings) + "|" + string.Join(",", x.BestModels)));
                }
                catch (InvalidOperationException ex)
                {
                    return "IOE:" + ex.GetType().Name;
                }
            }

            string baseRun = Run(segment, rSeed);
            return (baseRun == Run(segment, rSeed + 1) && baseRun == Run(shuffled, rSeed))
                .Label($"π={purity}, ccf={ccf}, pair={pair}");
        });
    }

    /// <summary>
    /// F64 (CNAqc <c>analyze_peaks</c> always runs BMix when NV/DP exist): with <see cref="OncologyAnalyzer.PurityPeakOptions.FitMixturePeaks"/>
    /// null, VAF-only input equals FitMixturePeaks = false and read-count input equals FitMixturePeaks = true (same seed).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 8)]
    public Property AnalyzePurityPeaks_MixtureAutoDefault_FollowsReadCountAvailability()
    {
        return Prop.ForAll(PeakScenarioGen(new[] { (1, 1), (2, 1) }).ToArbitrary(), s =>
        {
            var rng = new Random(s.Seed);
            var withReads = SimulateMutations(rng, s.Purity, s.Karyotypes.Take(1).ToArray(), 110, withReads: true);
            var vafOnly = withReads.Select(m => new OncologyAnalyzer.PurityPeakMutation(m.Vaf, m.MajorCopyNumber, m.MinorCopyNumber)).ToList();
            var dflt = OncologyAnalyzer.PurityPeakOptions.Default with { Seed = s.Seed % 1000 };

            static string Key(OncologyAnalyzer.PurityPeakAnalysis r) =>
                $"{r.Score:R}|{r.Pass}|" + string.Join(";", r.Karyotypes.SelectMany(k => k.Peaks).Select(p => $"{p.X:R},{p.Source},{p.Discarded}"));

            bool vafSide = Key(OncologyAnalyzer.AnalyzePurityPeaks(vafOnly, s.QcPurity, dflt))
                           == Key(OncologyAnalyzer.AnalyzePurityPeaks(vafOnly, s.QcPurity, dflt with { FitMixturePeaks = false }));
            bool readSide = Key(OncologyAnalyzer.AnalyzePurityPeaks(withReads, s.QcPurity, dflt))
                            == Key(OncologyAnalyzer.AnalyzePurityPeaks(withReads, s.QcPurity, dflt with { FitMixturePeaks = true }));
            return (vafSide && readSide).Label($"vafSide={vafSide}, readSide={readSide}");
        });
    }
}
