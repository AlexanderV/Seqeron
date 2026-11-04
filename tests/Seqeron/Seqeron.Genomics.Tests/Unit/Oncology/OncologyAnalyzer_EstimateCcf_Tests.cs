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

namespace Seqeron.Genomics.Tests.Unit.Oncology;

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
                "Centroids are independent of input order (optimal DP on the sorted data, no seeding).");
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

    // B24 review 2026-09 (ONCO-CCF-001, F17): the k-means objective (minimum WCSS, Lloyd 1982) is solved exactly
    // by the Ckmeans.1d.dp dynamic program (Wang & Song 2011, R Journal 3(2):29–33). The former Lloyd iteration
    // with quantile seeding stopped in a local optimum here: centroids {0.4167, 0.97, 1.0}, WCSS 0.071867 — it
    // merged the 0.2 minor subclone with the ~0.52 subclone and split the clonal cluster. Reference (Ckmeans.1d.dp
    // C++ via ckwrap 1.2.3): centers {0.2, 0.525, 0.98}, labels {2,2,2,1,1,0}, WCSS 0.0020500000000000036.
    [Test]
    public void ClusterCcfValues_LloydLocalOptimumCase_ReturnsCkmeansGlobalOptimum()
    {
        var values = new[] { 1.0, 0.98, 0.96, 0.55, 0.50, 0.20 };

        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, 3);

        double wcss = values.Select((v, i) => Math.Pow(v - result.Centroids[result.Assignments[i]], 2)).Sum();
        Assert.Multiple(() =>
        {
            Assert.That(result.Centroids, Is.EqualTo(new[] { 0.2, 0.525, 0.98 }).Within(1e-12),
                "Ckmeans.1d.dp optimum: minor subclone 0.2, subclone mean(0.55,0.50), clonal mean(1.0,0.98,0.96).");
            Assert.That(result.Assignments, Is.EqualTo(new[] { 2, 2, 2, 1, 1, 0 }), "Ckmeans.1d.dp labels.");
            Assert.That(result.ClonalClusterIndex, Is.EqualTo(2));
            Assert.That(wcss, Is.EqualTo(0.00205).Within(1e-12), "Global minimum WCSS (Lloyd local optimum: 0.071867).");
        });
    }

    // F17 second reference: {0.81,0.54,0.82,0.55,0.71,0.31}, k=3. Ckmeans.1d.dp: centers {0.31, 0.545, 0.78},
    // labels {2,1,2,1,2,0}, WCSS 0.00745; former Lloyd: {0.4667, 0.71, 0.815}, WCSS 0.036917.
    [Test]
    public void ClusterCcfValues_UnsortedInput_MatchesCkmeansReference()
    {
        var values = new[] { 0.81, 0.54, 0.82, 0.55, 0.71, 0.31 };

        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, 3);

        Assert.Multiple(() =>
        {
            Assert.That(result.Centroids, Is.EqualTo(new[] { 0.31, 0.545, 0.7799999999999999 }).Within(1e-12));
            Assert.That(result.Assignments, Is.EqualTo(new[] { 2, 1, 2, 1, 2, 0 }));
        });
    }

    // F17: Ckmeans.1d.dp sets Kmax = min(k, number of unique values) (R wrapper cluster.1d.dp and C++ kmeans_1d_dp),
    // so every returned cluster is non-empty. The former code returned centroids {0.5, 0.5, 1.0} with an empty
    // middle cluster.
    [Test]
    public void ClusterCcfValues_FewerDistinctValuesThanK_ReducesToDistinctCount()
    {
        var values = new[] { 0.5, 0.5, 0.5, 0.5, 1.0 };

        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, 3);

        Assert.Multiple(() =>
        {
            Assert.That(result.Centroids, Is.EqualTo(new[] { 0.5, 1.0 }), "k reduced to the 2 distinct values.");
            Assert.That(result.Assignments, Is.EqualTo(new[] { 0, 0, 0, 0, 1 }));
            Assert.That(result.ClonalClusterIndex, Is.EqualTo(1));
            Assert.That(OncologyAnalyzer.InferSubclones(result), Is.EqualTo(result.Centroids.Count),
                "No empty clusters.");
        });
    }

    // F17: the exact DP needs an effective-k × n backtrack matrix; above 10^8 cells the call is rejected up front
    // instead of exhausting memory.
    [Test]
    public void ClusterCcfValues_DpMatrixAboveCellLimit_Throws()
    {
        double[] values = Enumerable.Range(0, 20_001).Select(i => i / 20_000.0).ToArray();

        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClusterCcfValues(values, 5_000));
    }

    // B24 review 2026-09 (ONCO-CCF-001, F46): R's default Ckmeans.1d.dp(x, k) uses method = "linear" (SMAWK row fill,
    // EWL2_fill_SMAWK.cpp); among equal-WCSS optima its tie-breaking picks a different split than the former
    // "loglinear" fill. Here {0.5,0.5,0.5,0.75} and {0.75,1,1,1} both have WCSS 0.046875: R 4.3.6 default puts 0.75
    // with the 0.5s (centers {0, 0.25, 0.5625, 1}); method = "loglinear" (the former port) put it with the 1s
    // (centers {0, 0.25, 0.5, 0.9375}).
    [Test]
    public void ClusterCcfValues_EqualWcssTie_MatchesRDefaultLinearMethod()
    {
        var values = new[] { 0.5, 0.25, 1, 0.5, 0, 0.25, 0.75, 0.25, 0.5, 0, 1, 0, 0, 0, 0, 1 };

        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, 4);

        Assert.Multiple(() =>
        {
            Assert.That(result.Assignments, Is.EqualTo(new[] { 2, 1, 3, 2, 0, 1, 2, 1, 2, 0, 3, 0, 0, 0, 0, 3 }),
                "R Ckmeans.1d.dp(x, 4)$cluster - 1.");
            Assert.That(result.Centroids, Is.EqualTo(new[] { 0.0, 0.25, 0.5625, 1.0 }), "R $centers (exact).");
        });
    }

    // F46: a tie case from the 6000-input R comparison (one of the 21 where R "linear" and "loglinear" disagree);
    // labels and centers are R 4.3.6 Ckmeans.1d.dp(x, 8) output, compared bit for bit.
    [Test]
    public void ClusterCcfValues_EightClusterTieCase_BitIdenticalToR()
    {
        var values = new[]
        {
            0.4, 0.7, 0.3, 0.8, 0.8, 0, 0.5, 0.8, 0.5, 0.6, 0.1, 0.7, 0.2, 0.2, 0.6, 0, 0.5, 0.6, 0.7, 0.4, 0.9, 0.8, 0.9,
        };

        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, 8);

        Assert.Multiple(() =>
        {
            Assert.That(result.Assignments, Is.EqualTo(new[]
            {
                2, 5, 1, 6, 6, 0, 3, 6, 3, 4, 0, 5, 1, 1, 4, 0, 3, 4, 5, 2, 7, 6, 7,
            }));
            Assert.That(result.Centroids, Is.EqualTo(new[]
            {
                0.033333333333333333, 0.23333333333333331, 0.40000000000000002, 0.5, 0.59999999999999998,
                0.69999999999999984, 0.80000000000000004, 0.90000000000000002,
            }));
        });
    }

    // F46: equally spaced inputs have equal-WCSS splits; R 4.3.6 Ckmeans.1d.dp(x, k)$cluster - 1 (linear and
    // loglinear agree here).
    [TestCase(new[] { 0.0, 0.5, 1.0 }, 2, new[] { 0, 0, 1 })]
    [TestCase(new[] { 0.2, 0.4, 0.6, 0.8 }, 3, new[] { 0, 1, 1, 2 })]
    [TestCase(new[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 }, 4, new[] { 0, 1, 2, 2, 3, 3 })]
    [TestCase(new[] { 0.0, 0.25, 0.5, 0.75, 1.0 }, 3, new[] { 0, 0, 1, 1, 2 })]
    public void ClusterCcfValues_EquallySpacedTies_MatchR(double[] values, int k, int[] expected)
    {
        Assert.That(OncologyAnalyzer.ClusterCcfValues(values, k).Assignments, Is.EqualTo(expected));
    }

    // B24 review 2026-10 (ONCO-CCF-001, F47): automatic k = R Ckmeans.1d.dp(x, k = c(kmin, kmax)) (defaults
    // method = "linear", estimate.k = "BIC" → select_levels). Each row: values, kmin, kmax, R $cluster - 1, R $BIC
    // (R 4.3.6 output, set.seed(46) for the random rows). Covers k = 1 selected (rows 2, 6, 9, 19), kmax above the
    // number of distinct values (rows 3, 6, 19: R caps kmax), kmin above it (row 4: both bounds → 2), all values
    // equal (row 5: one cluster, BIC 0), ties (rows 3, 8, 16, 19, 20), n = 2 (row 6).
    private static readonly object[] RAutoKCases =
    {
            new object[] { new[] { 0.97999999999999998, 1.0, 0.96999999999999997, 1.02, 0.5, 0.52000000000000002, 0.47999999999999998, 0.51000000000000001 }, 1, 9, new[] { 1, 1, 1, 1, 0, 0, 0, 0 }, new[] { -4.4652525291599927, 20.839855023899851, 16.591100191352954, 10.365338648141801, 3.1857191312926822, -2.9052107672096952, -9.9231528049821236, -16.654140007666378 } },
            new object[] { new[] { 0.94999999999999996, 0.96999999999999997, 1.0, 0.98999999999999999, 0.95999999999999996, 1.01 }, 1, 5, new[] { 0, 0, 0, 0, 0, 0 }, new[] { 25.31466130849828, 23.032182234478128, 17.480340574682064, 11.726983533042723, 5.881096002825668 } },
            new object[] { new[] { 0.5, 0.5, 0.5, 1.0, 1.0 }, 1, 9, new[] { 0, 0, 0, 1, 1 }, new[] { -3.4569253296857951, 0.8823750858700059 } },
            new object[] { new[] { 0.29999999999999999, 0.29999999999999999, 0.59999999999999998, 0.59999999999999998 }, 3, 5, new[] { 0, 0, 1, 1 }, new[] { 4.1377007945553688 } },
            new object[] { new[] { 0.69999999999999996, 0.69999999999999996, 0.69999999999999996 }, 1, 3, new[] { 0, 0, 0 }, new[] { 0.0 } },
            new object[] { new[] { 0.20000000000000001, 0.90000000000000002 }, 1, 9, new[] { 0, 0 }, new[] { -3.2490543570637609, -6.5910710453828418 } },
            new object[] { new[] { 1.0, 0.97999999999999998, 0.95999999999999996, 0.55000000000000004, 0.5, 0.20000000000000001 }, 1, 9, new[] { 1, 1, 1, 0, 0, 0 }, new[] { -6.350563548530852, 1.1661367085263574, -1.0773909612144852, -7.6557349557707788, -13.33600481869254, -19.18862147506719 } },
            new object[] { new[] { 0.5, 0.25, 1.0, 0.5, 0.0, 0.25, 0.75, 0.25, 0.5, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0 }, 1, 9, new[] { 2, 1, 3, 2, 0, 1, 2, 1, 2, 0, 3, 0, 0, 0, 0, 3 }, new[] { -19.597290748854988, -21.476067854081172, -15.29592944814758, -12.82173398107139, -17.108438369771225 } },
            new object[] { new[] { 0.10000000000000001, 0.20000000000000001, 0.29999999999999999, 0.40000000000000002, 0.5, 0.59999999999999998 }, 1, 6, new[] { 0, 0, 0, 0, 0, 0 }, new[] { 0.50366196804415875, -3.8207489781603146, -8.7207093415845911, -14.825241440685256, -20.877157456368288, -26.019957860022597 } },
            new object[] { new[] { 0.81000000000000005, 0.54000000000000004, 0.81999999999999995, 0.55000000000000004, 0.70999999999999996, 0.31 }, 2, 4, new[] { 2, 1, 2, 1, 2, 0 }, new[] { -3.1187933895979896, -0.79522533232957215, -1.1468908823295827 } },
            new object[] { new[] { 0.97299999999999998, 1.006, 0.97799999999999998, 1.0369999999999999, 1.0349999999999999, 0.98099999999999998, 1.0129999999999999, 1.022, 1.0509999999999999, 1.0069999999999999, 0.47199999999999998, 0.51600000000000001, 0.39200000000000002, 0.439, 0.44400000000000001, 0.42199999999999999, 0.151, 0.161, 0.14199999999999999, 0.158 }, 1, 9, new[] { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 }, new[] { -21.50352366060028, 12.182499216098934, 29.546682492076489, 21.207302182296466, 18.725516424646216, 9.0371364015363227, 2.5050843844952198, -6.4448888561721489, -13.336443358682828 } },
            new object[] { new[] { 0.97999999999999998, 0.96099999999999997, 0.999, 1.081, 0.98299999999999998, 0.95699999999999996, 0.97999999999999998, 1.022, 0.93100000000000005, 1.0369999999999999, 0.96399999999999997, 1.0149999999999999, 0.57599999999999996, 0.57799999999999996, 0.65000000000000002, 0.47399999999999998, 0.69599999999999995 }, 1, 9, new[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 }, new[] { 2.7602154034215047, 20.378213290107006, 12.589844007694822, 5.2886161909918243, 7.4756383493138188, -0.93605417986473327, -9.1164393988463104, -18.807120470346121, -23.832979802189548 } },
            new object[] { new[] { 0.38, 0.72999999999999998, 0.25, 0.93999999999999995, 0.97999999999999998, 0.93000000000000005, 0.77000000000000002, 0.34999999999999998, 0.75, 0.14999999999999999, 0.089999999999999997, 0.59999999999999998, 0.46000000000000002, 0.90000000000000002, 0.93000000000000005 }, 1, 9, new[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, new[] { -11.886181231643665, -13.363128230455306, -20.531310480489388, -20.452413984305657, -24.481357451159937, -32.211039847926713, -39.131042656018238, -47.262358173272929, -55.847551616871172 } },
            new object[] { new[] { 0.95099999999999996, 0.84099999999999997, 0.77900000000000003, 0.85699999999999998, 0.83799999999999997, 0.88300000000000001, 0.91900000000000004, 0.86599999999999999, 0.88200000000000001, 0.95799999999999996, 0.871, 0.82199999999999995, 0.85299999999999998, 0.93400000000000005, 0.86699999999999999, 0.91600000000000004, 0.872, 0.92700000000000005, 0.83199999999999996, 0.89200000000000002 }, 1, 9, new[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, new[] { 61.484317400487441, 53.325983043333473, 46.510702990340079, 39.633146102564652, 34.736022408280924, 26.735201277227787, 22.104340652043142, 15.998133088576353, 7.2720333484649302 } },
            new object[] { new[] { 1.02, 0.99299999999999999, 1.0069999999999999, 1.0069999999999999, 1.0029999999999999, 1.0089999999999999, 0.98399999999999999, 1.004, 0.54800000000000004, 0.51400000000000001, 0.51200000000000001, 0.51200000000000001, 0.47599999999999998, 0.504, 0.46100000000000002, 0.47999999999999998, 0.25900000000000001, 0.218, 0.28899999999999998, 0.26000000000000001, 0.26200000000000001, 0.22800000000000001, 0.27000000000000002, 0.24399999999999999 }, 2, 6, new[] { 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, new[] { 25.188257961064259, 47.077685945827568, 39.064128273882893, 30.471312381728097, 26.359446461233212 } },
            new object[] { new[] { 0.5, 0.25, 0.25, 0.5, 0.5, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25, 0.5 }, 1, 9, new[] { 1, 0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 1 }, new[] { -7.3580474924646992, -3.4094143763100586, 7.2932231929778766 } },
            new object[] { new[] { 1.0249999999999999, 1.03, 0.96199999999999997, 1.0169999999999999, 0.877, 0.29999999999999999 }, 1, 4, new[] { 1, 1, 1, 1, 1, 0 }, new[] { -4.5255125032013979, -0.7315681945079735, -4.3677254085097967, -3.7912682670702278 } },
            new object[] { new[] { 1.1220000000000001, 1.0229999999999999, 1.0780000000000001, 1.1160000000000001, 1.0660000000000001, 0.98599999999999999, 0.97599999999999998, 1.0169999999999999, 0.98899999999999999, 1.0700000000000001, 1.0640000000000001, 1.048, 1.006, 1.0349999999999999, 1.006, 0.999, 0.98999999999999999, 0.91900000000000004, 0.89400000000000002, 1.0289999999999999, 1.0009999999999999, 0.92200000000000004, 1.0069999999999999, 0.98099999999999998, 1.0229999999999999, 0.92100000000000004, 0.93700000000000006, 1.0369999999999999, 1.0589999999999999, 0.93100000000000005, 0.40799999999999997, 0.249, 0.32900000000000001, 0.35699999999999998, 0.32900000000000001, 0.24199999999999999, 0.375, 0.253, 0.307, 0.29799999999999999 }, 3, 9, new[] { 2, 2, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, new[] { 40.785253454003623, 39.835975924481936, 29.545544811843179, 27.178106891454803, 23.787555965670066, 16.137683747365671, 4.6971464277115587 } },
            new object[] { new[] { 0.0, 0.0, 0.0, 1.0 }, 1, 9, new[] { 0, 0, 0, 0 }, new[] { -7.5789195433975998, -7.8322586329281041 } },
            new object[] { new[] { 1.0, 1.0, 1.0, 0.97999999999999998, 0.5 }, 1, 2, new[] { 1, 1, 1, 1, 0 }, new[] { -1.3367428612177954, 13.10861570059506 } },
    };

    [TestCaseSource(nameof(RAutoKCases))]
    public void ClusterCcfValues_AutomaticK_MatchesRCkmeansBic(
        double[] values, int minClusters, int maxClusters, int[] expectedLabels, double[] expectedBic)
    {
        OncologyAnalyzer.CcfClustering result = OncologyAnalyzer.ClusterCcfValues(values, minClusters, maxClusters);
        (OncologyAnalyzer.CcfClustering withBic, double[] bic) =
            OncologyAnalyzer.ClusterCcfValuesSelectingK(values, minClusters, maxClusters);

        int expectedK = expectedLabels.Max() + 1;
        Assert.Multiple(() =>
        {
            Assert.That(result.Assignments, Is.EqualTo(expectedLabels), "R $cluster - 1.");
            Assert.That(result.Centroids, Has.Count.EqualTo(expectedK), "Chosen k = R max($cluster).");
            Assert.That(result.ClonalClusterIndex, Is.EqualTo(expectedK - 1));
            Assert.That(bic, Is.EqualTo(expectedBic), "R $BIC, bit for bit.");
            Assert.That(withBic.Assignments, Is.EqualTo(result.Assignments));
            for (int c = 0; c < expectedK; c++)
            {
                double[] members = values.Where((_, i) => expectedLabels[i] == c).ToArray();
                Assert.That(result.Centroids[c], Is.EqualTo(members.Average()).Within(1e-15), "Centre = block mean.");
            }
        });
    }

    // F47: for minClusters = maxClusters the BIC overload is the fixed-k clustering (same single DP fill).
    [Test]
    public void ClusterCcfValues_AutomaticKWithEqualBounds_EqualsFixedK()
    {
        var values = new[] { 0.5, 0.25, 1, 0.5, 0, 0.25, 0.75, 0.25, 0.5, 0, 1, 0, 0, 0, 0, 1 };

        OncologyAnalyzer.CcfClustering auto = OncologyAnalyzer.ClusterCcfValues(values, 4, 4);
        OncologyAnalyzer.CcfClustering fixedK = OncologyAnalyzer.ClusterCcfValues(values, 4);

        Assert.Multiple(() =>
        {
            Assert.That(auto.Assignments, Is.EqualTo(fixedK.Assignments));
            Assert.That(auto.Centroids, Is.EqualTo(fixedK.Centroids));
        });
    }

    // F47: R stops for k.max <= 0; the port requires 1 <= minClusters <= maxClusters (maxClusters may exceed the
    // number of values — R caps it at the number of distinct values).
    [Test]
    public void ClusterCcfValues_AutomaticK_InvalidArguments_Throw()
    {
        var values = new[] { 0.2, 0.9 };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClusterCcfValues(values, 0, 3));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClusterCcfValues(values, 3, 2));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.ClusterCcfValues(null!, 1, 9));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClusterCcfValues(Array.Empty<double>(), 1, 9));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClusterCcfValues(new[] { 0.5, double.NaN }, 1, 9));
            Assert.DoesNotThrow(() => OncologyAnalyzer.ClusterCcfValues(values, 1, 50));
        });
    }

    #endregion

    // B24 review 2026-09 (ONCO-PURITY-001 dedup): EstimateCcf routes through the canonical CNAqc
    // purity/copy-number correction AdjustVAFForPurity (n_mut = VAF·(ρ·N_T + 2(1−ρ))/ρ, McGranahan 2016),
    // so RawCcf = n_mut / m exactly. Reference (Python, McGranahan 2016 formula):
    // VAF 0.3, ρ 0.7, N_T 3, m 2 ⇒ n_mut = 0.3·2.7/0.7 = 1.157142857…, CCF = 0.578571428…
    [Test]
    public void EstimateCcf_RawCcf_EqualsCanonicalPurityCorrectionOverMultiplicity()
    {
        OncologyAnalyzer.CcfEstimate estimate = OncologyAnalyzer.EstimateCcf(0.3, 0.7, 3, 2);

        Assert.Multiple(() =>
        {
            Assert.That(estimate.RawCcf, Is.EqualTo(OncologyAnalyzer.AdjustVAFForPurity(0.3, 0.7, 3) / 2));
            Assert.That(estimate.RawCcf, Is.EqualTo(0.3 * 2.7 / 0.7 / 2.0).Within(1e-12));
            Assert.That(OncologyAnalyzer.DeriveMultiplicity(0.3, 0.7, 3, 2), Is.EqualTo(1),
                "n_mut = 1.157… rounds to multiplicity 1 (McGranahan 2016)");
        });
    }
}
