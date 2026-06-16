// ONCO-CCF-001 — Cancer Cell Fraction Estimation and CCF Clustering
// Evidence: docs/Evidence/ONCO-CCF-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CCF-001.md
// Source: McGranahan N et al. (2016). Science 351(6280):1463-1469. https://doi.org/10.1126/science.aaf1490
//         Tarabichi M et al. (2021). Nat. Methods 18:144-155 (Box 1). PMC7867630
//         Zheng J et al. (2022). Bioinformatics 38(15):3677-3683. https://doi.org/10.1093/bioinformatics/btac367
//         Lloyd SP (1982). IEEE Trans. Inf. Theory 28(2):129-137. https://doi.org/10.1109/TIT.1982.1056489
//
// Expected CCF values are computed independently from CCF = VAF·(rho·N_T + 2(1-rho))/(rho·m) — NOT
// copied from the implementation. Cluster centroids/assignments are derived from Lloyd's k-means on
// the sorted values; the clonal cluster is the highest centroid (Tarabichi 2021).

using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_EstimateCcf_Tests
{
    private const double Tolerance = 1e-10;

    #region EstimateCcf

    // M1 — f=0.40, rho=0.80, N_T=2, m=1 -> 0.40·(0.8·2+2·0.2)/(0.8·1) = 0.40·2.0/0.8 = 1.0
    [Test]
    public void EstimateCcf_ClonalDiploid_ReturnsOne()
    {
        OncologyAnalyzer.CcfEstimate result = OncologyAnalyzer.EstimateCcf(0.40, 0.80, 2, 1);

        Assert.Multiple(() =>
        {
            Assert.That(result.Ccf, Is.EqualTo(1.0).Within(Tolerance),
                "CCF = 0.40·2.0/0.8 = 1.0; a clonal diploid heterozygous SNV at purity 0.8.");
            Assert.That(result.RawCcf, Is.EqualTo(1.0).Within(Tolerance),
                "Raw equals capped here because the formula value is exactly 1.0.");
        });
    }

    // M2 — f=0.20, rho=0.80, N_T=2, m=1 -> 0.20·2.0/0.8 = 0.5 (subclonal)
    [Test]
    public void EstimateCcf_Subclonal_ReturnsHalf()
    {
        OncologyAnalyzer.CcfEstimate result = OncologyAnalyzer.EstimateCcf(0.20, 0.80, 2, 1);

        Assert.That(result.Ccf, Is.EqualTo(0.5).Within(Tolerance),
            "CCF = 0.20·(0.8·2+2·0.2)/0.8 = 0.40/0.8 = 0.5; mutation in ~half the cancer cells.");
    }

    // M3 — f=0.50, rho=1.0, N_T=4, m=2 -> 0.50·(1·4+0)/(1·2) = 2.0/2.0 = 1.0 (multi-copy)
    [Test]
    public void EstimateCcf_MultiCopyLocus_ReturnsOne()
    {
        OncologyAnalyzer.CcfEstimate result = OncologyAnalyzer.EstimateCcf(0.50, 1.0, 4, 2);

        Assert.That(result.Ccf, Is.EqualTo(1.0).Within(Tolerance),
            "CCF = 0.50·(1·4+2·0)/(1·2) = 1.0; multiplicity 2 on a 4-copy locus is clonal.");
    }

    // M4 — f=0.25, rho=0.50, N_T=2, m=1 -> 0.25·(0.5·2+2·0.5)/0.5 = 0.5/0.5 = 1.0
    [Test]
    public void EstimateCcf_HalfPurityDiploid_ReturnsOne()
    {
        OncologyAnalyzer.CcfEstimate result = OncologyAnalyzer.EstimateCcf(0.25, 0.50, 2, 1);

        Assert.That(result.Ccf, Is.EqualTo(1.0).Within(Tolerance),
            "CCF = 0.25·(0.5·2+2·0.5)/0.5 = 0.50/0.50 = 1.0; clonal at 50% purity.");
    }

    // M5 — f=0.471, rho=1.0, N_T=2, m=1 -> 0.471·2.0/1.0 = 0.942 (raw < 1, no cap)
    [Test]
    public void EstimateCcf_RawBelowOne_NotCapped()
    {
        OncologyAnalyzer.CcfEstimate result = OncologyAnalyzer.EstimateCcf(0.471, 1.0, 2, 1);

        Assert.Multiple(() =>
        {
            Assert.That(result.Ccf, Is.EqualTo(0.942).Within(Tolerance),
                "CCF = 0.471·(1·2+0)/(1·1) = 0.942; below 1 so reported uncapped.");
            Assert.That(result.RawCcf, Is.EqualTo(0.942).Within(Tolerance),
                "Raw equals reported when below 1.");
        });
    }

    // M6 — f=0.60, rho=0.80, N_T=2, m=1 -> raw 0.60·2.0/0.8 = 1.5; reported capped to 1.0 (INV-CCF-01)
    [Test]
    public void EstimateCcf_RawAboveOne_ReportedCappedRawExposed()
    {
        OncologyAnalyzer.CcfEstimate result = OncologyAnalyzer.EstimateCcf(0.60, 0.80, 2, 1);

        Assert.Multiple(() =>
        {
            Assert.That(result.Ccf, Is.EqualTo(1.0).Within(Tolerance),
                "Reported CCF capped at 1.0 to honour 0<=CCF<=1; raw 1.5 exceeds 1 from over-sampled VAF.");
            Assert.That(result.RawCcf, Is.EqualTo(1.5).Within(Tolerance),
                "Raw = 0.60·(0.8·2+2·0.2)/0.8 = 1.5 exposed uncapped.");
        });
    }

    // C1 — f=0 -> CCF = 0 (INV-CCF-03)
    [Test]
    public void EstimateCcf_ZeroVaf_ReturnsZero()
    {
        OncologyAnalyzer.CcfEstimate result = OncologyAnalyzer.EstimateCcf(0.0, 0.80, 2, 1);

        Assert.That(result.Ccf, Is.EqualTo(0.0).Within(Tolerance),
            "VAF = 0 implies no mutated reads, so CCF = 0.");
    }

    // CNAqc real worked outputs (Caravagna lab, "Computation of Cancer Cell Fractions" vignette,
    // sample purity 89%, diploid N_T=2): VAF=0.08/m=1 -> 0.180, VAF=0.883/m=2 -> 0.993, VAF=0.471/m=1 -> 1.06.
    // These are reference-implementation outputs (not derived from this code), reproduced to 1e-2 by
    // CCF = VAF·(ρ·N_T + 2(1−ρ))/(ρ·m) with ρ=0.89, N_T=2.
    // https://caravagnalab.github.io/CNAqc/articles/a4_ccf_computation.html
    [Test]
    public void EstimateCcf_CnaqcWorkedOutputs_MatchReferenceImplementation()
    {
        const double cnaqcTolerance = 5e-3;

        OncologyAnalyzer.CcfEstimate low = OncologyAnalyzer.EstimateCcf(0.08, 0.89, 2, 1);
        OncologyAnalyzer.CcfEstimate mid = OncologyAnalyzer.EstimateCcf(0.883, 0.89, 2, 2);
        OncologyAnalyzer.CcfEstimate high = OncologyAnalyzer.EstimateCcf(0.471, 0.89, 2, 1);

        Assert.Multiple(() =>
        {
            Assert.That(low.RawCcf, Is.EqualTo(0.180).Within(cnaqcTolerance),
                "CNAqc: VAF=0.08, multiplicity=1, CCF=0.180.");
            Assert.That(mid.RawCcf, Is.EqualTo(0.993).Within(cnaqcTolerance),
                "CNAqc: VAF=0.883, multiplicity=2, CCF=0.993.");
            Assert.That(high.RawCcf, Is.EqualTo(1.06).Within(cnaqcTolerance),
                "CNAqc: VAF=0.471, multiplicity=1, CCF=1.06 (raw exceeds 1 from sampling noise).");
            Assert.That(high.Ccf, Is.EqualTo(1.0).Within(Tolerance),
                "Reported CCF is capped at 1.0 even though the raw CNAqc value is 1.06.");
        });
    }

    // S2 — INV-CCF-02: CCF strictly increases with VAF holding other inputs fixed.
    [Test]
    public void EstimateCcf_IncreasingVaf_IncreasesCcf()
    {
        double low = OncologyAnalyzer.EstimateCcf(0.10, 0.80, 2, 1).RawCcf;
        double high = OncologyAnalyzer.EstimateCcf(0.20, 0.80, 2, 1).RawCcf;

        Assert.That(high, Is.GreaterThan(low),
            "Formula is linear in VAF with positive slope; doubling VAF doubles raw CCF.");
    }

    // M7 — purity outside (0,1] rejected.
    [Test]
    public void EstimateCcf_InvalidPurity_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimateCcf(0.4, 0.0, 2, 1),
                "Purity 0 divides by zero in the formula and must be rejected.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimateCcf(0.4, 1.2, 2, 1),
                "Purity > 1 is not a valid fraction.");
        });
    }

    // M8 — VAF outside [0,1] rejected.
    [Test]
    public void EstimateCcf_InvalidVaf_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimateCcf(-0.1, 0.8, 2, 1),
                "Negative VAF is not a valid fraction.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimateCcf(1.1, 0.8, 2, 1),
                "VAF > 1 is not a valid fraction.");
        });
    }

    // M9 — tumor copy number < 1 rejected.
    [Test]
    public void EstimateCcf_InvalidCopyNumber_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimateCcf(0.4, 0.8, 0, 1),
            "Tumor copy number must be at least 1.");
    }

    // M10 — multiplicity outside [1, copyNumber] rejected.
    [Test]
    public void EstimateCcf_InvalidMultiplicity_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.EstimateCcf(0.4, 0.8, 2, 0),
                "Multiplicity 0 means no mutated copies; not valid.");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.EstimateCcf(0.4, 0.8, 2, 3),
                "Multiplicity cannot exceed the tumor copy number.");
        });
    }

    #endregion

    #region ClusterCcfValues

    // M11/M13 — {1.0,0.98,0.96,0.50,0.48,0.52}, k=2: centroids {0.50,0.98};
    //           low cluster (idx 0) = inputs 3,4,5; high cluster (idx 1) = inputs 0,1,2; clonal = idx 1.
    [Test]
    public void ClusterCcfValues_TwoClones_ReturnsExactCentroidsAndAssignments()
    {
        var values = new[] { 1.0, 0.98, 0.96, 0.50, 0.48, 0.52 };

        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, 2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Centroids[0], Is.EqualTo(0.50).Within(Tolerance),
                "Low centroid = mean(0.48,0.50,0.52) = 0.50.");
            Assert.That(result.Centroids[1], Is.EqualTo(0.98).Within(Tolerance),
                "High centroid = mean(0.96,0.98,1.0) = 0.98.");
            Assert.That(result.Assignments, Is.EqualTo(new[] { 1, 1, 1, 0, 0, 0 }),
                "First three (clonal CCF ~1) -> high cluster 1; last three (~0.5) -> low cluster 0.");
            Assert.That(result.ClonalClusterIndex, Is.EqualTo(1),
                "Clonal cluster is the highest centroid (0.98), index 1.");
        });
    }

    // M12 — determinism: shuffled input groups each value into the same centroid.
    [Test]
    public void ClusterCcfValues_ShuffledInput_ProducesIdenticalCentroids()
    {
        var ordered = new[] { 1.0, 0.98, 0.96, 0.50, 0.48, 0.52 };
        var shuffled = new[] { 0.48, 1.0, 0.52, 0.96, 0.50, 0.98 };

        OncologyAnalyzer.CcfClustering a = OncologyAnalyzer.ClusterCcfValues(ordered, 2);
        OncologyAnalyzer.CcfClustering b = OncologyAnalyzer.ClusterCcfValues(shuffled, 2);

        Assert.Multiple(() =>
        {
            Assert.That(b.Centroids, Is.EqualTo(a.Centroids),
                "Centroids are independent of input order (deterministic quantile seeding).");
            // shuffled[1]=1.0 and shuffled[3]=0.96 are clonal -> cluster 1; shuffled[0]=0.48 -> cluster 0.
            Assert.That(b.Assignments[1], Is.EqualTo(1), "1.0 always lands in the high (clonal) cluster.");
            Assert.That(b.Assignments[0], Is.EqualTo(0), "0.48 always lands in the low cluster.");
        });
    }

    // S1 — k=1: single cluster at the global mean; clonal index 0.
    [Test]
    public void ClusterCcfValues_SingleCluster_ReturnsGlobalMean()
    {
        var values = new[] { 0.3, 0.6, 0.9 };

        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, 1);

        Assert.Multiple(() =>
        {
            Assert.That(result.Centroids, Has.Count.EqualTo(1), "k=1 yields one cluster.");
            Assert.That(result.Centroids[0], Is.EqualTo(0.6).Within(Tolerance),
                "Single centroid = mean(0.3,0.6,0.9) = 0.6.");
            Assert.That(result.Assignments, Is.EqualTo(new[] { 0, 0, 0 }),
                "All values map to the only cluster.");
            Assert.That(result.ClonalClusterIndex, Is.EqualTo(0), "The only cluster is clonal.");
        });
    }

    // M14 — null input rejected.
    [Test]
    public void ClusterCcfValues_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.ClusterCcfValues(null!, 2),
            "Null value list is rejected.");
    }

    // S3 — empty input rejected.
    [Test]
    public void ClusterCcfValues_Empty_Throws()
    {
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClusterCcfValues(Array.Empty<double>(), 1),
            "At least one CCF value is required.");
    }

    // M15 — k outside [1, count] rejected.
    [Test]
    public void ClusterCcfValues_InvalidClusterCount_Throws()
    {
        var values = new[] { 0.3, 0.6 };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClusterCcfValues(values, 0),
                "k must be at least 1.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClusterCcfValues(values, 3),
                "k cannot exceed the number of values.");
        });
    }

    // Validation — non-finite value rejected.
    [Test]
    public void ClusterCcfValues_NonFiniteValue_Throws()
    {
        var values = new[] { 0.3, double.NaN };

        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClusterCcfValues(values, 1),
            "NaN cannot be clustered by squared distance.");
    }

    // Validation — infinite value rejected (separate branch from NaN).
    [Test]
    public void ClusterCcfValues_InfiniteValue_Throws()
    {
        var values = new[] { 0.3, double.PositiveInfinity };

        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClusterCcfValues(values, 1),
            "Infinite CCF cannot be clustered by squared distance.");
    }

    // INV-5 / relabeling — clonal cluster is always the highest centroid and centroids are returned
    // ascending, regardless of input order. Here the high-CCF (clonal) values are listed FIRST, which
    // exercises the ascending-centroid relabeling rather than an identity mapping.
    // Centroids: mean(0.10,0.12,0.14)=0.12 and mean(0.90,0.92,0.94)=0.92; clonal = high = index 1.
    [Test]
    public void ClusterCcfValues_HighCcfValuesFirst_ClonalIsHighestCentroidAscending()
    {
        var values = new[] { 0.94, 0.92, 0.90, 0.14, 0.12, 0.10 };

        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, 2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Centroids[0], Is.EqualTo(0.12).Within(Tolerance),
                "Centroids returned ascending: low = mean(0.10,0.12,0.14) = 0.12.");
            Assert.That(result.Centroids[1], Is.EqualTo(0.92).Within(Tolerance),
                "High centroid = mean(0.90,0.92,0.94) = 0.92.");
            Assert.That(result.Assignments, Is.EqualTo(new[] { 1, 1, 1, 0, 0, 0 }),
                "Input order preserved: the first three (high CCF) map to the high cluster 1.");
            Assert.That(result.ClonalClusterIndex, Is.EqualTo(1),
                "Clonal cluster is the highest centroid (0.92) at index 1 (Tarabichi 2021).");
        });
    }

    #endregion
}
