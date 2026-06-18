// ONCO-HETERO-001 — Tumor Heterogeneity Analysis (MATH, Shannon diversity, subclone count, subclonal fraction)
// Evidence: docs/Evidence/ONCO-HETERO-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-HETERO-001.md
// Source: Mroz EA, Rocco JW (2013). Oral Oncology 49(3):211-215. https://pubmed.ncbi.nlm.nih.gov/23079694/
//         Mroz EA et al. (2015). PLOS Medicine 12(2):e1001786. https://doi.org/10.1371/journal.pmed.1001786
//         maftools mathScore.R: pat.math = (median(abs(vaf-median(vaf)))*100)*1.4826/median(vaf)
//         Liu Z, Zhang S (2017). BMC Genomics 18:457 (PMC5468233) — Shannon H = -sum p_i ln(p_i)
//         Landau DA et al. (2013). Cell 152(4):714-726 — subclonal iff CCF < 0.95
//
// Expected MATH values are derived independently from MATH = 100*1.4826*median(|f-median(f)|)/median(f),
// and Shannon values from H = -sum p_i ln(p_i) (natural log) over clone fractions — NOT from the implementation.

using System;
using System.Collections.Generic;

namespace Seqeron.Genomics.Tests;

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

    // M9b — subclonal threshold boundary: CCF exactly 0.95 is clonal (strict CCF < 0.95, Landau 2013).
    //       CCFs {0.94, 0.95, 0.96, 0.97}: only 0.94 is strictly below 0.95 => fraction 1/4 = 0.25.
    [Test]
    public void AnalyzeHeterogeneity_SubclonalThresholdBoundary_ExcludesExactly0Point95()
    {
        var vafs = new[] { 0.47, 0.475, 0.48, 0.485 };
        var ccf = new[] { 0.94, 0.95, 0.96, 0.97 };

        OncologyAnalyzer.HeterogeneityResult result =
            OncologyAnalyzer.AnalyzeHeterogeneity(vafs, ccf, clusterCount: 2);

        Assert.That(result.SubclonalFraction, Is.EqualTo(0.25).Within(Tolerance),
            "Subclonal iff CCF < 0.95 (strict); CCF=0.95 is clonal, so only 0.94 counts => 1/4 = 0.25 (Landau 2013).");
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
}
