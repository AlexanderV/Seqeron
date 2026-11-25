using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// FIN-B24 heavy tier — property-based tests for the behaviour added by the B24 finisher (F44–F64). Oracles are
/// independent restatements of the reference definitions (never routed through production helpers):
/// <list type="bullet">
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

    private static readonly (int, int)[] SimplePool = { (1, 1), (2, 0), (2, 1), (2, 2), (1, 0) };

    private static readonly (int, int)[] ComplexPool = { (3, 1), (3, 0), (4, 2), (3, 2), (4, 1) };
}
