// ONCO-EXPR-001 — Tumor Gene Expression Outlier (z-score) and Signature Score
// Evidence: docs/Evidence/ONCO-EXPR-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-EXPR-001.md
// Source: cBioPortal mRNA z-score normalization spec — z = (r-mu)/sigma. https://docs.cbioportal.org/z-score-normalization-script/
//         cBioPortal-core NormalizeExpressionLevels.java — sample SD divisor (n-1); sigma=0 fatal error.
//         cBioPortal FAQ — default threshold ±2, strict ">2 or <-2". https://docs.cbioportal.org/user-guide/faq/
//         Lee E et al. (2008) PLoS Comput Biol 4(11):e1000217 — combined z-score a = (sum z)/sqrt(k). https://doi.org/10.1371/journal.pcbi.1000217
//
// Expected values are derived independently from z = (x - mean)/sd (sd = sample SD, divisor n-1) and
// a = (sum z)/sqrt(k) — NOT from the implementation. Reference cohort {2,2,4,6,6}: mean=4, SS=16,
// var=16/4=4, sd=2; chosen so the sample-SD result (sd=2) differs from the population-SD result (sd=1.789).

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_IdentifyOutlierGenes_Tests
{
    private const double Tolerance = 1e-10;

    // Reference cohort with clean sample-SD statistics: mean=4, sd=2 (var=16/(5-1)=4).
    private static IReadOnlyList<double> Cohort() => new[] { 2.0, 2.0, 4.0, 6.0, 6.0 };

    #region CalculateExpressionZScore

    // M1 — cohort {2,2,4,6,6} (mean=4, sd=2), x=10: z=(10-4)/2 = 3.0
    [Test]
    public void CalculateExpressionZScore_OverexpressedValue_Returns3()
    {
        double z = OncologyAnalyzer.CalculateExpressionZScore(10.0, Cohort());

        Assert.That(z, Is.EqualTo(3.0).Within(Tolerance),
            "z = (10-4)/2 = 3.0 with mean=4 and sample SD=2 (cBioPortal z=(r-mu)/sigma).");
    }

    // M2 — x equal to the mean: z=(4-4)/2 = 0.0 (INV-01)
    [Test]
    public void CalculateExpressionZScore_ValueAtMean_ReturnsZero()
    {
        double z = OncologyAnalyzer.CalculateExpressionZScore(4.0, Cohort());

        Assert.That(z, Is.EqualTo(0.0).Within(Tolerance),
            "z at the cohort mean is 0 (INV-01).");
    }

    // M3 — x=-1: z=(-1-4)/2 = -2.5 (under-expressed, |z|>2)
    [Test]
    public void CalculateExpressionZScore_UnderexpressedValue_ReturnsNegative2Point5()
    {
        double z = OncologyAnalyzer.CalculateExpressionZScore(-1.0, Cohort());

        Assert.That(z, Is.EqualTo(-2.5).Within(Tolerance),
            "z = (-1-4)/2 = -2.5; below -2 is an under-expression outlier.");
    }

    // M4 — sample SD (n-1) must be used: x=6 -> z=(6-4)/2 = 1.0.
    //      A population-SD implementation (divisor n) would give sd=sqrt(16/5)=1.78885 and z=1.11803.
    [Test]
    public void CalculateExpressionZScore_UsesSampleStandardDeviation_NotPopulation()
    {
        double z = OncologyAnalyzer.CalculateExpressionZScore(6.0, Cohort());

        Assert.That(z, Is.EqualTo(1.0).Within(Tolerance),
            "Sample SD (divisor n-1) gives sd=2 -> z=1.0; population SD would give 1.118 (NormalizeExpressionLevels.java std()).");
    }

    // S1 — null reference cohort
    [Test]
    public void CalculateExpressionZScore_NullCohort_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.CalculateExpressionZScore(1.0, null!), "Null reference cohort is invalid.");
    }

    // S2 — cohort of size 1: sample SD (n-1) undefined
    [Test]
    public void CalculateExpressionZScore_CohortOfSizeOne_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.CalculateExpressionZScore(1.0, new[] { 5.0 }),
            "A single reference value has no sample standard deviation (n-1).");
    }

    // M11 — zero-SD (constant) cohort throws, mirroring the reference implementation's fatal error
    [Test]
    public void CalculateExpressionZScore_ZeroStandardDeviationCohort_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.CalculateExpressionZScore(7.0, new[] { 5.0, 5.0, 5.0 }),
            "A constant cohort (sd=0) has no defined z-score (NormalizeExpressionLevels.java fatalError).");
    }

    // C1 — INV-03: reflecting the value about the mean negates z. x=10 -> z=3; 2*mean-x = 8-10 = -2 -> z=-3.
    [Test]
    public void CalculateExpressionZScore_ReflectionAboutMean_NegatesZScore()
    {
        double zHigh = OncologyAnalyzer.CalculateExpressionZScore(10.0, Cohort());
        double zLow = OncologyAnalyzer.CalculateExpressionZScore(-2.0, Cohort()); // 2*4 - 10 = -2

        Assert.That(zLow, Is.EqualTo(-zHigh).Within(Tolerance),
            "z(2*mean - x) = -z(x): z(-2) = -z(10) = -3.0 (INV-03).");
    }

    #endregion

    #region IdentifyOutlierGenes

    private static IReadOnlyDictionary<string, IReadOnlyList<double>> Cohorts() =>
        new Dictionary<string, IReadOnlyList<double>>
        {
            ["GENE_OVER"] = Cohort(),
            ["GENE_UNDER"] = Cohort(),
            ["GENE_BOUND"] = Cohort(),
            ["GENE_NORMAL"] = Cohort(),
        };

    // M5/M6/M8 — over (z=3.0), under (z=-2.5), normal (z=0) classified correctly.
    [Test]
    public void IdentifyOutlierGenes_OverAndUnder_ClassifiedNonOutlierExcluded()
    {
        var sample = new Dictionary<string, double>
        {
            ["GENE_OVER"] = 10.0,    // z = 3.0
            ["GENE_UNDER"] = -1.0,   // z = -2.5
            ["GENE_NORMAL"] = 4.0,   // z = 0.0
        };

        var outliers = OncologyAnalyzer.IdentifyOutlierGenes(sample, Cohorts());

        var map = new Dictionary<string, OncologyAnalyzer.ExpressionOutlier>();
        foreach (var o in outliers)
        {
            map[o.Gene] = o;
        }

        Assert.Multiple(() =>
        {
            Assert.That(outliers, Has.Count.EqualTo(2), "Only GENE_OVER and GENE_UNDER exceed |z|>2.");
            Assert.That(map.ContainsKey("GENE_OVER"), Is.True, "z=3.0 > 2 is an outlier.");
            Assert.That(map["GENE_OVER"].Direction, Is.EqualTo(OncologyAnalyzer.ExpressionDirection.Over),
                "z=3.0 is overexpressed (cBioPortal >2).");
            Assert.That(map["GENE_OVER"].ZScore, Is.EqualTo(3.0).Within(Tolerance), "GENE_OVER z = 3.0.");
            Assert.That(map.ContainsKey("GENE_UNDER"), Is.True, "z=-2.5 < -2 is an outlier.");
            Assert.That(map["GENE_UNDER"].Direction, Is.EqualTo(OncologyAnalyzer.ExpressionDirection.Under),
                "z=-2.5 is underexpressed (cBioPortal <-2).");
            Assert.That(map["GENE_UNDER"].ZScore, Is.EqualTo(-2.5).Within(Tolerance), "GENE_UNDER z = -2.5.");
            Assert.That(map.ContainsKey("GENE_NORMAL"), Is.False, "z=0 is not an outlier.");
        });
    }

    // M7 — boundary: x=8 -> z=(8-4)/2 = 2.0 exactly; strict >2 means NOT an outlier.
    [Test]
    public void IdentifyOutlierGenes_ZScoreExactlyAtThreshold_NotOutlier()
    {
        var sample = new Dictionary<string, double> { ["GENE_BOUND"] = 8.0 }; // z = 2.0

        var outliers = OncologyAnalyzer.IdentifyOutlierGenes(sample, Cohorts());

        Assert.That(outliers, Is.Empty,
            "z=2.0 equals the threshold; the rule is strict (>2 / <-2), so it is not an outlier (cBioPortal FAQ).");
    }

    // S5 — null sample / null cohorts
    [Test]
    public void IdentifyOutlierGenes_NullArguments_Throw()
    {
        var sample = new Dictionary<string, double> { ["G"] = 1.0 };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.IdentifyOutlierGenes(null!, Cohorts()), "Null sample expression is invalid.");
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.IdentifyOutlierGenes(sample, null!), "Null reference cohorts is invalid.");
        });
    }

    // S6 — sampled gene with no reference cohort
    [Test]
    public void IdentifyOutlierGenes_MissingReferenceCohort_Throws()
    {
        var sample = new Dictionary<string, double> { ["UNKNOWN_GENE"] = 10.0 };

        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.IdentifyOutlierGenes(sample, Cohorts()),
            "A gene with no reference cohort cannot be z-scored.");
    }

    // Non-positive threshold guard
    [Test]
    public void IdentifyOutlierGenes_NonPositiveThreshold_Throws()
    {
        var sample = new Dictionary<string, double> { ["GENE_OVER"] = 10.0 };

        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.IdentifyOutlierGenes(sample, Cohorts(), 0.0), "Threshold must be positive.");
    }

    #endregion

    #region CalculateSignatureScore

    // M9 — combined z-score: z={3,1,-1,1}, k=4: a = (3+1-1+1)/sqrt(4) = 4/2 = 2.0
    [Test]
    public void CalculateSignatureScore_FourGenes_Returns2()
    {
        double a = OncologyAnalyzer.CalculateSignatureScore(new[] { 3.0, 1.0, -1.0, 1.0 });

        Assert.That(a, Is.EqualTo(2.0).Within(Tolerance),
            "a = (sum z)/sqrt(k) = 4/sqrt(4) = 2.0 (Lee et al. 2008).");
    }

    // M10 — single-gene signature: a = z/sqrt(1) = z (INV-06)
    [Test]
    public void CalculateSignatureScore_SingleGene_ReturnsThatZScore()
    {
        double a = OncologyAnalyzer.CalculateSignatureScore(new[] { 2.5 });

        Assert.That(a, Is.EqualTo(2.5).Within(Tolerance),
            "a = 2.5/sqrt(1) = 2.5 (INV-06).");
    }

    // C2 — INV-05: k equal z-scores all = c give a = c*sqrt(k). c=1.5, k=4: a = 6/2 = 3.0.
    [Test]
    public void CalculateSignatureScore_EqualZScores_ReturnsCTimesSqrtK()
    {
        double a = OncologyAnalyzer.CalculateSignatureScore(new[] { 1.5, 1.5, 1.5, 1.5 });

        Assert.That(a, Is.EqualTo(3.0).Within(Tolerance),
            "a = (4*1.5)/sqrt(4) = 1.5*2 = 3.0 (INV-05; sqrt(k) denominator, not k).");
    }

    // S3 — empty signature
    [Test]
    public void CalculateSignatureScore_EmptySignature_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.CalculateSignatureScore(Array.Empty<double>()),
            "An empty signature (k=0) has no defined combined z-score.");
    }

    // S4 — null signature
    [Test]
    public void CalculateSignatureScore_NullSignature_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.CalculateSignatureScore(null!), "Null signature is invalid.");
    }

    #endregion
}
