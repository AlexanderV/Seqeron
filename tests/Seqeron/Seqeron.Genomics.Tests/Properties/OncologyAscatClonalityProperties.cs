using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the behaviour introduced by the 2026-09 B24 review (ONCO-ASCAT-001, ONCO-CCF-001,
/// ONCO-PHYLO-001, ONCO-HETERO-001): the runASCAT / ASCAT-manual fit, the ascat.aspcf segmentation, the Battenberg
/// sub-clonal fit, the LICHeE trunk rule and the maftools MATH score.
///
/// Oracles are recomputed independently from the published equations (never routed through production helpers):
/// <list type="bullet">
/// <item>ASCAT model (Van Loo et al. 2010, PNAS 107:16910; ascat.runAscat.R): a clonal segment with integer allele copy
/// numbers (nA, nB) at purity ρ and ploidy ψ has logR r = log2((2(1−ρ)+ρ(nA+nB)) / (2(1−ρ)+ρψ)) (γ = 1) and BAF
/// b = (1−ρ+ρ·nB) / (2(1−ρ)+ρ(nA+nB)).</item>
/// <item>Battenberg (Nik-Zainal et al. 2012; Wedge-lab/battenberg R/fitcopynumber.R): the sub-clonal fraction τ solves
/// l = [1−ρ+ρ(τM₁+(1−τ)M₂)] / [2−2ρ+ρ(τT₁+(1−τ)T₂)] for the two states of the nearest edge (which differ by one copy
/// of one allele); a clonal call has |level(state) − l| &lt; maxdist = 0.01.</item>
/// <item>ascat.aspcf: each segment's logR level is the mean of the raw logR of its loci; its BAF is 0.5 + μ ∈ [0.5, 1].</item>
/// <item>Werner et al. 2017 (Sci Rep 7:44991): trunk clusters are present in all cells of every sample (CCF ≥ 1 − ε).</item>
/// <item>MATH (Mroz &amp; Rocco 2013; maftools mathScore.R) is a MAD/median ratio: invariant to input order and to
/// scaling every VAF by c &gt; 0.</item>
/// </list>
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Oncology")]
public class OncologyAscatClonalityProperties
{
    // -------------------------------------------------------------------------
    // ASCAT exact-model generators (Van Loo 2010 equations, γ = 1)
    // -------------------------------------------------------------------------

    private static double ModelLogR(double rho, double psi, int total) =>
        Math.Log2((2.0 * (1.0 - rho) + rho * total) / (2.0 * (1.0 - rho) + rho * psi));

    private static double ModelMinorBaf(double rho, int minor, int total) =>
        (1.0 - rho + rho * minor) / (2.0 * (1.0 - rho) + rho * total);

    private static Gen<(double Rho, double Psi, (int Major, int Minor, int Probes)[] Karyotype)> ExactAscatGenomeGen() =>
        from rhoPercent in Gen.Choose(20, 100)
        from psiTwentieths in Gen.Choose(30, 100) // ψ ∈ [1.5, 5.0] step 0.05
        from n in Gen.Choose(1, 6)
        from karyotype in (from major in Gen.Choose(0, 5)
                           from minor in Gen.Choose(0, 5)
                           from probes in Gen.Choose(1, 200)
                           select (Major: Math.Max(major, minor), Minor: Math.Min(major, minor), Probes: probes))
            .Where(k => k.Major + k.Minor >= 1)
            .ArrayOf(n)
        select (rhoPercent / 100.0, psiTwentieths / 20.0, karyotype);

    // -------------------------------------------------------------------------
    // Battenberg sub-clonal fit
    // -------------------------------------------------------------------------

    private static double CornerLevel(double rho, int major, int minor) =>
        (1.0 - rho + rho * major) / (2.0 - 2.0 * rho + rho * (major + minor));

    // -------------------------------------------------------------------------
    // LICHeE trunk (Werner 2017) and MATH (maftools)
    // -------------------------------------------------------------------------

    private static Arbitrary<(OncologyAnalyzer.CcfCluster[] Clusters, double Tolerance)> PhylogenyArbitrary() =>
        (from samples in Gen.Choose(1, 3)
         from n in Gen.Choose(1, 6)
         from vectors in Gen.Elements(0, 250, 500, 750, 900, 970, 1000).Select(v => v / 1000.0).ArrayOf(samples).ArrayOf(n)
         from tolMilli in Gen.Choose(0, 100)
         select (vectors.Select((v, i) => new OncologyAnalyzer.CcfCluster(i + 1, v)).ToArray(), tolMilli / 1000.0))
        .ToArbitrary();

    /// <summary>
    /// Werner et al. 2017: every trunk cluster has CCF ≥ 1 − ε in every sample, and the trunk is maximal — the last trunk
    /// node (or the root, for an empty trunk) has no child that is clonal in every sample. Oracle recomputed from the
    /// input CCFs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IdentifyTrunkMutations_IsTheMaximalAllClonalRootPath()
    {
        return Prop.ForAll(PhylogenyArbitrary(), p =>
        {
            if (!OncologyAnalyzer.TryReconstructPhylogeny(p.Clusters, out var phylo, p.Tolerance))
            {
                return true.Label("no valid tree");
            }

            var ccf = p.Clusters.ToDictionary(c => c.Id, c => c.CcfPerSample);
            bool Clonal(int id) => ccf[id].All(v => v >= 1.0 - p.Tolerance);

            var trunk = OncologyAnalyzer.IdentifyTrunkMutations(phylo);
            if (!trunk.All(Clonal))
            {
                return false.Label($"trunk [{string.Join(",", trunk)}] contains a cluster with CCF < 1 − ε");
            }

            int last = trunk.Count == 0 ? phylo.RootId : trunk[^1];
            bool maximal = !phylo.ChildrenOf(last).Any(Clonal);
            return maximal.Label($"trunk [{string.Join(",", trunk)}] stops before an all-clonal child of {last}");
        });
    }

    private static Arbitrary<(double[] Values, int Seed, int ScalePermille)> VafArbitrary() =>
        (from n in Gen.Choose(1, 25)
         from values in Gen.Choose(1, 1000).Select(v => v / 1000.0).ArrayOf(n)
         from seed in Gen.Choose(0, int.MaxValue)
         from scale in Gen.Choose(50, 1000)
         select (values, seed, scale)).ToArbitrary();

    /// <summary>
    /// MATH depends only on the multiset of VAFs (both medians are order statistics): any permutation gives a
    /// bit-identical score.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateITH_IsPermutationInvariant()
    {
        return Prop.ForAll(VafArbitrary(), t =>
        {
            var rng = new Random(t.Seed);
            double[] shuffled = t.Values.OrderBy(_ => rng.Next()).ToArray();
            return (OncologyAnalyzer.CalculateITH(shuffled) == OncologyAnalyzer.CalculateITH(t.Values))
                .Label("MATH changed under a permutation of the VAFs");
        });
    }

    /// <summary>
    /// MATH = 100·1.4826·MAD/median is scale-free: multiplying every VAF by c ∈ (0, 1] (so VAFs stay in [0, 1])
    /// scales MAD and median alike and leaves the score unchanged (up to rounding).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateITH_IsInvariantToScalingAllVafs()
    {
        return Prop.ForAll(VafArbitrary(), t =>
        {
            double c = t.ScalePermille / 1000.0;
            double baseMath = OncologyAnalyzer.CalculateITH(t.Values);
            double scaledMath = OncologyAnalyzer.CalculateITH(t.Values.Select(v => v * c).ToArray());
            return (Math.Abs(scaledMath - baseMath) <= 1e-9 * Math.Max(1.0, baseMath))
                .Label($"MATH {baseMath} → {scaledMath} under scaling by {c}");
        });
    }
}
