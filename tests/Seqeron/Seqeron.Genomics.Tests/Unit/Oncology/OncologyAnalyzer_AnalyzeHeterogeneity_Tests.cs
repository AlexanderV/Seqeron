// ONCO-HETERO-001 — Tumor Heterogeneity Analysis (MATH, Shannon diversity, subclone count, subclonal fraction)
// Evidence: docs/Evidence/ONCO-HETERO-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-HETERO-001.md
// Source: Mroz EA, Rocco JW (2013). Oral Oncology 49(3):211-215. https://pubmed.ncbi.nlm.nih.gov/23079694/
//         Mroz EA et al. (2015). PLOS Medicine 12(2):e1001786. https://doi.org/10.1371/journal.pmed.1001786
//         maftools mathScore.R: pat.math = (median(abs(vaf-median(vaf)))*100)*1.4826/median(vaf)
//         Martinez P et al. (2017). Sci Rep 7:3248 (PMC5468233) — Shannon H = -sum p_i ln(p_i) over clonal frequencies
//         Landau DA et al. (2013). Cell 152(4):714-726 — clonal iff CCF > 0.95, "subclonal otherwise" (CCF <= 0.95)
//
// Expected MATH values are derived independently from MATH = 100*1.4826*median(|f-median(f)|)/median(f),
// and Shannon values from H = -sum p_i ln(p_i) (natural log) over clone fractions — NOT from the implementation.

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_AnalyzeHeterogeneity_Tests
{
    private const double Tolerance = 1e-10;

    #region CalculateITH (MATH score)

    // M1 — VAFs {0.1,0.2,0.3,0.4,0.5}: median=0.30; absdev {0.2,0.1,0,0.1,0.2}; rawMAD=0.10;
    //      MATH = 100*1.4826*0.10/0.30 = 49.42
    [Test]
    public void CalculateITH_OddCountWorkedExample_Returns49Point42()
    {
        double math = OncologyAnalyzer.CalculateITH(new[] { 0.10, 0.20, 0.30, 0.40, 0.50 });

        Assert.That(math, Is.EqualTo(49.42).Within(Tolerance),
            "MATH = 100*1.4826*median(|f-0.30|)/0.30 = 100*1.4826*0.10/0.30 = 49.42 (Mroz & Rocco 2013).");
    }

    // M2 — VAFs {0.2,0.4,0.6,0.8}: median=(0.4+0.6)/2=0.50; absdev {0.3,0.1,0.1,0.3}; rawMAD=(0.1+0.3)/2=0.20;
    //      MATH = 100*1.4826*0.20/0.50 = 59.304
    [Test]
    public void CalculateITH_EvenCountWorkedExample_Returns59Point304()
    {
        double math = OncologyAnalyzer.CalculateITH(new[] { 0.20, 0.40, 0.60, 0.80 });

        Assert.That(math, Is.EqualTo(59.304).Within(Tolerance),
            "Even-count median=0.50, rawMAD=0.20; MATH = 100*1.4826*0.20/0.50 = 59.304 (maftools mathScore.R).");
    }

    // M3 — all identical VAFs: MAD = 0 => MATH = 0 (INV-02)
    [Test]
    public void CalculateITH_AllIdenticalVafs_ReturnsZero()
    {
        double math = OncologyAnalyzer.CalculateITH(new[] { 0.30, 0.30, 0.30 });

        Assert.That(math, Is.EqualTo(0.0).Within(Tolerance),
            "Every VAF equals the median so MAD=0 and MATH=0 (no heterogeneity).");
    }

    // M4 — single VAF: median=value, MAD=0 => MATH=0
    [Test]
    public void CalculateITH_SingleVaf_ReturnsZero()
    {
        double math = OncologyAnalyzer.CalculateITH(new[] { 0.40 });

        Assert.That(math, Is.EqualTo(0.0).Within(Tolerance),
            "A single mutation has median=its VAF and MAD=0, so MATH=0.");
    }

    // S1 — median of 0 => MATH undefined (division by zero) => ArgumentException
    [Test]
    public void CalculateITH_ZeroMedian_Throws()
    {
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateITH(new[] { 0.0, 0.0, 0.40 }),
            "median(VAF)=0 makes MATH=100*MAD/0 undefined; must throw.");
    }

    // S2 — null distribution
    [Test]
    public void CalculateITH_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateITH(null!),
            "Null distribution is invalid input.");
    }

    // S3 — empty distribution
    [Test]
    public void CalculateITH_Empty_Throws()
    {
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateITH(Array.Empty<double>()),
            "An empty distribution has no median; must throw.");
    }

    // S4 — out-of-range VAF
    [Test]
    public void CalculateITH_OutOfRangeVaf_Throws()
    {
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateITH(new[] { 0.30, 1.50 }),
            "VAF must be in [0,1]; 1.50 is invalid.");
    }

    // C1 — INV-01: MATH >= 0 for varied valid inputs
    [Test]
    public void CalculateITH_VariedInputs_AlwaysNonNegative()
    {
        double[][] inputs =
        {
            new[] { 0.05, 0.5, 0.95 },
            new[] { 0.1, 0.1, 0.9, 0.9 },
            new[] { 0.33, 0.34, 0.35 },
        };

        Assert.Multiple(() =>
        {
            foreach (double[] vafs in inputs)
            {
                Assert.That(OncologyAnalyzer.CalculateITH(vafs), Is.GreaterThanOrEqualTo(0.0),
                    "Registry invariant ITH_score >= 0: MAD>=0 and median>0.");
            }
        });
    }

    #endregion

    #region InferSubclones

    // M8 — three well-separated CCFs clustered into k=3 => 3 occupied clusters
    [Test]
    public void InferSubclones_ThreeSeparatedClusters_ReturnsThree()
    {
        OncologyAnalyzer.CcfClustering clustering =
            OncologyAnalyzer.ClusterCcfValues(new[] { 0.20, 0.50, 0.90 }, 3);

        int count = OncologyAnalyzer.InferSubclones(clustering);

        Assert.That(count, Is.EqualTo(3),
            "Three distinct CCF clusters each contain one mutation => richness 3 (Liu & Zhang 2017).");
    }

    // M8b — single cluster => richness 1
    [Test]
    public void InferSubclones_SingleCluster_ReturnsOne()
    {
        OncologyAnalyzer.CcfClustering clustering =
            OncologyAnalyzer.ClusterCcfValues(new[] { 0.20, 0.22, 0.24 }, 1);

        int count = OncologyAnalyzer.InferSubclones(clustering);

        Assert.That(count, Is.EqualTo(1), "One cluster => one clone (monoclonal).");
    }

    // S6 — empty clustering throws
    [Test]
    public void InferSubclones_EmptyClustering_Throws()
    {
        var empty = new OncologyAnalyzer.CcfClustering(Array.Empty<double>(), Array.Empty<int>(), 0);

        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.InferSubclones(empty),
            "A clustering with no centroids/assignments is invalid.");
    }

    #endregion

    #region AnalyzeHeterogeneity (aggregate)

    // M5 — CCFs {0.20,0.22,0.90,0.92}, k=2 => two clusters of size 2 => p={0.5,0.5} => H = -ln 0.5 = 0.6931471805599453
    [Test]
    public void AnalyzeHeterogeneity_TwoEqualClones_ShannonIsLn2()
    {
        var vafs = new[] { 0.10, 0.11, 0.45, 0.46 };
        var ccf = new[] { 0.20, 0.22, 0.90, 0.92 };

        OncologyAnalyzer.HeterogeneityResult result =
            OncologyAnalyzer.AnalyzeHeterogeneity(vafs, ccf, clusterCount: 2);

        Assert.Multiple(() =>
        {
            Assert.That(result.SubcloneCount, Is.EqualTo(2),
                "Two well-separated CCF groups => 2 clones.");
            Assert.That(result.ShannonDiversity, Is.EqualTo(0.6931471805599453).Within(Tolerance),
                "Two equal clones (p=0.5 each): H = -2*(0.5*ln0.5) = -ln0.5 = ln2 (Shannon 1948).");
            Assert.That(result.SubclonalFraction, Is.EqualTo(1.0).Within(Tolerance),
                "All four CCFs < 0.95 => subclonal fraction 1.0 (Landau 2013).");
        });
    }

    // M6 — CCFs {0.10,0.40,0.70,0.95}, k=4 => four clusters size 1 => H = ln 4 = 1.3862943611198906
    [Test]
    public void AnalyzeHeterogeneity_FourEqualClones_ShannonIsLn4()
    {
        var vafs = new[] { 0.05, 0.20, 0.35, 0.48 };
        var ccf = new[] { 0.10, 0.40, 0.70, 0.95 };

        OncologyAnalyzer.HeterogeneityResult result =
            OncologyAnalyzer.AnalyzeHeterogeneity(vafs, ccf, clusterCount: 4);

        Assert.Multiple(() =>
        {
            Assert.That(result.SubcloneCount, Is.EqualTo(4), "Four distinct CCFs => 4 clones.");
            Assert.That(result.ShannonDiversity, Is.EqualTo(1.3862943611198906).Within(Tolerance),
                "Four equal clones (p=0.25 each): H = -4*(0.25*ln0.25) = ln4 (Shannon 1948).");
        });
    }

    // M7 — single clone (k=1) => H = 0
    [Test]
    public void AnalyzeHeterogeneity_SingleClone_ShannonIsZero()
    {
        var vafs = new[] { 0.40, 0.42, 0.44 };
        var ccf = new[] { 0.90, 0.92, 0.94 };

        OncologyAnalyzer.HeterogeneityResult result =
            OncologyAnalyzer.AnalyzeHeterogeneity(vafs, ccf, clusterCount: 1);

        Assert.Multiple(() =>
        {
            Assert.That(result.SubcloneCount, Is.EqualTo(1), "k=1 => one clone.");
            Assert.That(result.ShannonDiversity, Is.EqualTo(0.0).Within(Tolerance),
                "One clone (p=1): H = -1*ln1 = 0 (Shannon 1948).");
        });
    }

    // M9 — subclonal fraction: CCFs {0.40,0.50,0.98,1.0}, two below 0.95 => 0.5
    [Test]
    public void AnalyzeHeterogeneity_SubclonalFraction_IsHalf()
    {
        var vafs = new[] { 0.20, 0.25, 0.49, 0.50 };
        var ccf = new[] { 0.40, 0.50, 0.98, 1.00 };

        OncologyAnalyzer.HeterogeneityResult result =
            OncologyAnalyzer.AnalyzeHeterogeneity(vafs, ccf, clusterCount: 2);

        Assert.That(result.SubclonalFraction, Is.EqualTo(0.5).Within(Tolerance),
            "Exactly 2 of 4 CCFs are < 0.95 (0.40, 0.50) => fraction 0.5 (Landau 2013).");
    }

    // M9b — subclonal threshold boundary (F20): Landau et al. (2013) "classified a mutation as clonal if the CCF
    //       harboring it was >0.95 ... and subclonal otherwise", so CCF exactly 0.95 is SUBCLONAL (same rule as the
    //       canonical IdentifyClonalMutations). CCFs {0.94, 0.95, 0.96, 0.97}: 0.94 and 0.95 => 2/4 = 0.5.
    [Test]
    public void AnalyzeHeterogeneity_SubclonalThresholdBoundary_Exactly0Point95IsSubclonal()
    {
        var vafs = new[] { 0.47, 0.475, 0.48, 0.485 };
        var ccf = new[] { 0.94, 0.95, 0.96, 0.97 };

        OncologyAnalyzer.HeterogeneityResult result =
            OncologyAnalyzer.AnalyzeHeterogeneity(vafs, ccf, clusterCount: 2);

        Assert.That(result.SubclonalFraction, Is.EqualTo(0.5).Within(Tolerance),
            "Clonal iff CCF > 0.95 (Landau 2013); 0.94 and 0.95 are subclonal => 2/4 = 0.5.");
    }

    // F20 — single mutation at CCF 0.95: not clonal (Landau 2013) => subclonal fraction 1; old code reported 0.
    [Test]
    public void AnalyzeHeterogeneity_SingleCcfAtThreshold_AgreesWithIdentifyClonalMutations()
    {
        OncologyAnalyzer.HeterogeneityResult result =
            OncologyAnalyzer.AnalyzeHeterogeneity(new[] { 0.475 }, new[] { 0.95 }, clusterCount: 1);

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IdentifyClonalMutations(new[] { 0.95 }), Is.Empty);
            Assert.That(result.SubclonalFraction, Is.EqualTo(1.0));
        });
    }

    // Shannon for unequal clones: CCFs {0.20, 0.21, 0.22, 0.90}, k = 2 => cluster sizes {3, 1};
    // scipy.stats.entropy([3, 1]) = skbio shannon([3, 1], base=e) = 0.5623351446188083.
    [Test]
    public void AnalyzeHeterogeneity_UnequalClones_ShannonMatchesScipy()
    {
        OncologyAnalyzer.HeterogeneityResult result = OncologyAnalyzer.AnalyzeHeterogeneity(
            new[] { 0.10, 0.11, 0.12, 0.45 }, new[] { 0.20, 0.21, 0.22, 0.90 }, clusterCount: 2);

        Assert.Multiple(() =>
        {
            Assert.That(result.SubcloneCount, Is.EqualTo(2));
            Assert.That(result.ShannonDiversity, Is.EqualTo(0.5623351446188083).Within(1e-15));
        });
    }

    // M10 — aggregate consistency: MATH component equals CalculateITH on the same VAFs
    [Test]
    public void AnalyzeHeterogeneity_MathComponent_MatchesCalculateITH()
    {
        var vafs = new[] { 0.10, 0.20, 0.30, 0.40, 0.50 };
        var ccf = new[] { 0.20, 0.22, 0.90, 0.92, 0.95 };

        OncologyAnalyzer.HeterogeneityResult result =
            OncologyAnalyzer.AnalyzeHeterogeneity(vafs, ccf, clusterCount: 2);

        Assert.That(result.MathScore, Is.EqualTo(49.42).Within(Tolerance),
            "Aggregate MATH over these VAFs equals the independently derived 49.42 (Mroz & Rocco 2013).");
    }

    // S5 — mismatched VAF/CCF lengths
    [Test]
    public void AnalyzeHeterogeneity_MismatchedLengths_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.AnalyzeHeterogeneity(new[] { 0.1, 0.2 }, new[] { 0.5 }, 1),
            "VAF and CCF lists must be aligned (same length).");
    }

    // S5b — null inputs
    [Test]
    public void AnalyzeHeterogeneity_NullInputs_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.AnalyzeHeterogeneity(null!, new[] { 0.5 }, 1), "Null VAFs invalid.");
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.AnalyzeHeterogeneity(new[] { 0.5 }, null!, 1), "Null CCFs invalid.");
        });
    }

    #endregion

    #region F21 — bit-exact maftools conformance, F22 — InferSubclones label validation, canonical helpers

    // F21: R 4.x, maftools mathScore.R arithmetic (median(abs.med.dev)*100)*1.4826/median(vaf), printed %.17g.
    //      {0.16, 0.87} is a case where the former 100*(1.4826*MAD)/median order differed by 1 ulp (…902).
    [TestCase(new[] { 0.16, 0.87 }, 102.19864077669901)]
    [TestCase(new[] { 0.12, 0.31, 0.07, 0.45, 0.26, 0.39 }, 70.22842105263156)]
    public void CalculateITH_MatchesMaftoolsBitExact(double[] vafs, double expected)
    {
        Assert.That(OncologyAnalyzer.CalculateITH(vafs), Is.EqualTo(expected));
    }

    // F22: an assignment label with no centroid is not a valid clustering (count would exceed k).
    [TestCase(1)]
    [TestCase(-1)]
    public void InferSubclones_LabelOutsideCentroids_Throws(int badLabel)
    {
        var clustering = new OncologyAnalyzer.CcfClustering(new[] { 0.3 }, new[] { 0, badLabel }, 0);
        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.InferSubclones(clustering));
    }

    // Canonical helpers introduced for this unit (StatisticsHelper) — numpy.median / scipy.stats.entropy values.
    [Test]
    public void StatisticsHelper_Median_MatchesNumpy()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.Median(new[] { 0.2, 0.4, 0.6, 0.8 }), Is.EqualTo(0.5));
            Assert.That(StatisticsHelper.Median(new[] { 3.0, 1.0, 2.0 }), Is.EqualTo(2.0));
            Assert.That(StatisticsHelper.Median(new[] { 1.0, double.NaN }), Is.NaN);
            Assert.That(StatisticsHelper.Median(new[] { 1e308, 1.7e308 }), Is.EqualTo(1.35e308)); // R median
            Assert.Throws<ArgumentException>(() => StatisticsHelper.Median(Array.Empty<double>()));
        });
    }

    [Test]
    public void StatisticsHelper_ShannonIndex_MatchesScipyEntropy()
    {
        Assert.Multiple(() =>
        {
            Assert.That(StatisticsHelper.ShannonIndex(new[] { 2, 2 }), Is.EqualTo(0.6931471805599453).Within(1e-15));
            Assert.That(StatisticsHelper.ShannonIndex(new[] { 3, 1, 0 }), Is.EqualTo(0.5623351446188083).Within(1e-15));
            Assert.That(StatisticsHelper.ShannonIndex(new[] { 1, 2, 3, 4 }), Is.EqualTo(1.2798542258336676).Within(1e-15));
            Assert.That(StatisticsHelper.ShannonIndex(new[] { 5 }), Is.EqualTo(0.0));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.ShannonIndex(new[] { 0, 0 }));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.ShannonIndex(new[] { 1, -1 }));
        });
    }

    #endregion

    #region F56 — maftools math.score pre-filters, F57 — prevalence-weighted Shannon

    // F56: expected values are maftools R/mathScore.R math.score (unmodified source, R 4.3.3) printed with %.17g.
    // R: dat = dat[!t_vaf < vafCutOff]; if (length(vaf) < 5) skip; pat.math = median(|vaf-median|)*100*1.4826/median(vaf).

    // P1 — no VAF below 0.075: filter is a no-op; R MATH = 49.420000000000016 (= unfiltered overload)
    [Test]
    public void CalculateITH_MaftoolsFilters_NothingFiltered_MatchesR()
    {
        double[] vafs = { 0.10, 0.20, 0.30, 0.40, 0.50 };
        double? math = OncologyAnalyzer.CalculateITH(vafs, OncologyAnalyzer.MaftoolsMathVafCutOff, OncologyAnalyzer.MaftoolsMathMinMutations);

        Assert.Multiple(() =>
        {
            Assert.That(math, Is.EqualTo(49.420000000000016));
            Assert.That(math, Is.EqualTo(OncologyAnalyzer.CalculateITH(vafs)));
        });
    }

    // P2 — VAF exactly 0.075 kept (!t_vaf < 0.075), 0.074 dropped; 6 retained; R MATH = 74.129999999999967
    [Test]
    public void CalculateITH_MaftoolsFilters_VafExactlyAtCutoffKept_MatchesR()
    {
        double? math = OncologyAnalyzer.CalculateITH(
            new[] { 0.075, 0.074, 0.12, 0.25, 0.33, 0.41, 0.48 }, OncologyAnalyzer.MaftoolsMathVafCutOff);

        Assert.That(math, Is.EqualTo(74.129999999999967),
            "maftools keeps VAF = cutoff: {0.075,0.12,0.25,0.33,0.41,0.48} median 0.29, raw MAD 0.145 -> 74.13.");
    }

    // P3 — S3a: 6 VAFs, 2 below cutoff -> 4 < 5 -> maftools skips (no row); S3b: exactly 5 retained -> 60.651818181818179
    [Test]
    public void CalculateITH_MaftoolsFilters_TooFewAfterCutoff_SkippedLikeR()
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CalculateITH(new[] { 0.05, 0.06, 0.20, 0.30, 0.40, 0.45 }, 0.075), Is.Null);
            Assert.That(OncologyAnalyzer.CalculateITH(new[] { 0.08, 0.15, 0.22, 0.31, 0.44 }, 0.075),
                Is.EqualTo(60.651818181818179));
        });
    }

    // P4 — 0.01 and 0.0749999 dropped, 6 retained (even count); R MATH = 85.614929577464792
    [Test]
    public void CalculateITH_MaftoolsFilters_EvenCountAfterCutoff_MatchesR()
    {
        double? math = OncologyAnalyzer.CalculateITH(
            new[] { 0.01, 0.0749999, 0.09, 0.11, 0.35, 0.36, 0.52, 0.61 }, OncologyAnalyzer.MaftoolsMathVafCutOff);

        Assert.That(math, Is.EqualTo(85.614929577464792));
    }

    // P5 — 4 VAFs, none filtered -> maftools skips (< 5)
    [Test]
    public void CalculateITH_MaftoolsFilters_FourMutations_Skipped()
    {
        Assert.That(OncologyAnalyzer.CalculateITH(new[] { 0.2, 0.3, 0.4, 0.5 }, 0.075), Is.Null);
    }

    // P6 — maftools VAF from counts t_alt/(t_ref+t_alt): ref {90,95,70,60,45,30}, alt {10,5,30,40,55,20},
    //      vafCutOff = 0.1 (10/100 exactly at cutoff kept, 5/100 dropped); R MATH = 37.065000000000005
    [Test]
    public void CalculateITH_MaftoolsFilters_VafFromCountsCustomCutoff_MatchesR()
    {
        double[] refCounts = { 90, 95, 70, 60, 45, 30 };
        double[] altCounts = { 10, 5, 30, 40, 55, 20 };
        double[] vafs = refCounts.Select((r, i) => altCounts[i] / (r + altCounts[i])).ToArray();

        Assert.That(OncologyAnalyzer.CalculateITH(vafs, 0.1), Is.EqualTo(37.065000000000005));
    }

    [Test]
    public void CalculateITH_MaftoolsFilters_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateITH(null!, 0.075));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateITH(new[] { 0.2 }, -0.01));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateITH(new[] { 0.2 }, double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.CalculateITH(new[] { 0.2 }, 0.075, 0));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateITH(new[] { 0.2, 1.5 }, 0.075, 1));
            // vafCutOff = 0 keeps zeros: median 0 -> MATH undefined, as in the unfiltered overload.
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateITH(new[] { 0.0, 0.0, 0.3 }, 0.0, 1));
        });
    }

    // F57 — Martinez et al. (2017) Shannon over clonal (cellular) frequencies; expected = scipy.stats.entropy 1.18.1.
    [TestCase(new[] { 0.6, 0.3, 0.1 }, 0.8979457248567798)]
    [TestCase(new[] { 0.45, 0.35, 0.15, 0.05 }, 1.161120818283116)]
    [TestCase(new[] { 0.7, 0.0, 0.3 }, 0.6108643020548935)]
    [TestCase(new[] { 0.2, 0.1, 0.05 }, 0.9556998911125343)] // unnormalised (sum 0.35): scipy normalises
    [TestCase(new[] { 0.25, 0.25, 0.25, 0.25 }, 1.3862943611198906)]
    [TestCase(new[] { 1.0 }, 0.0)]
    public void CalculateCloneShannonDiversity_MatchesScipyEntropy(double[] fractions, double expected)
    {
        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.CalculateCloneShannonDiversity(fractions), Is.EqualTo(expected).Within(1e-15));
            Assert.That(StatisticsHelper.ShannonIndexOfWeights(fractions), Is.EqualTo(expected).Within(1e-15));
        });
    }

    // Prevalence- vs count-weighting differ: 2 clones with 8 and 2 mutations but cellular frequencies 0.3 / 0.7.
    [Test]
    public void CalculateCloneShannonDiversity_DiffersFromCountWeighting()
    {
        double prevalence = OncologyAnalyzer.CalculateCloneShannonDiversity(new[] { 0.3, 0.7 });
        double counts = StatisticsHelper.ShannonIndex(new[] { 8, 2 });

        Assert.Multiple(() =>
        {
            Assert.That(prevalence, Is.EqualTo(0.6108643020548935).Within(1e-15)); // scipy entropy([0.3,0.7])
            Assert.That(counts, Is.EqualTo(0.5004024235381879).Within(1e-15));     // scipy entropy([8,2])
        });
    }

    [Test]
    public void CalculateCloneShannonDiversity_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.CalculateCloneShannonDiversity(null!));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateCloneShannonDiversity(Array.Empty<double>()));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateCloneShannonDiversity(new[] { 0.0, 0.0 }));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateCloneShannonDiversity(new[] { 0.5, 1.2 }));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.CalculateCloneShannonDiversity(new[] { 0.5, double.NaN }));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.ShannonIndexOfWeights(new[] { 0.5, -0.1 }));
            Assert.Throws<ArgumentException>(() => StatisticsHelper.ShannonIndexOfWeights(new[] { double.PositiveInfinity }));
        });
    }

    #endregion
}
