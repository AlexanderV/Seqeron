using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// B09 A2-31: calculate_similarity sketch estimators, distances and ANI. The expected numbers are the exact k-mer values of the Python set oracle
/// (identical rows are locked in GenomicAnalyzer_SimilaritySketch_Oracle_Tests against sourmash 4.9.4) or follow from the definitions.
/// </summary>
[TestFixture]
[Category("McpCore")]
public class CalculateSimilarityEstimatorsTests
{
    [Test]
    public void CalculateSimilarity_Exact_ReportsDistancesContainmentAndAni()
    {
        // k = 2: {AT, TG, GC} vs {AT, TG, GA}: J = 2/4 = 0.5, C(A,B) = C(B,A) = 2/3
        var r = SuffixTreeGenomicsTools.CalculateSimilarity("ATGC", "ATGA", 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.Similarity, Is.EqualTo(50.0).Within(1e-9));
            Assert.That(r.Estimator, Is.EqualTo("Exact"));
            Assert.That(r.JaccardDistance, Is.EqualTo(0.5).Within(1e-12));
            Assert.That(r.MashDistance, Is.EqualTo(-Math.Log(2 * 0.5 / 1.5) / 2).Within(1e-12));
            Assert.That(r.JaccardAni, Is.EqualTo(Math.Pow(2 * 0.5 / 1.5, 0.5)).Within(1e-12));
            Assert.That(r.ContainmentAInB, Is.EqualTo(2.0 / 3).Within(1e-12));
            Assert.That(r.ContainmentBInA, Is.EqualTo(2.0 / 3).Within(1e-12));
            Assert.That(r.MaxContainment, Is.EqualTo(2.0 / 3).Within(1e-12));
            Assert.That(r.ContainmentAniAInB, Is.EqualTo(Math.Pow(2.0 / 3, 0.5)).Within(1e-12));
            Assert.That(r.SharedHashes, Is.Null);
        });
    }

    [Test]
    public void CalculateSimilarity_Sketches_SizeLargerThanTheUnionIsExact()
    {
        const string a = "ACGTACGTTAGCTAGCTAGGATCGATCGATTTAGCGCGATATCGGCTAGCTAGCTAGGGATCTAGCTAGCTTAGCGA";
        const string b = "ACGTACGTTAGCTAGCTAGGATCGATCGATTTAGCGCGATATCGGCTAGGGGATCTAGCTAGCTTAGCGATTTTTT";
        var exact = SuffixTreeGenomicsTools.CalculateSimilarity(a, b, 21, canonical: true);
        var mash = SuffixTreeGenomicsTools.CalculateSimilarity(a, b, 21, canonical: true, sketchSize: 100000);
        var frac = SuffixTreeGenomicsTools.CalculateSimilarity(a, b, 21, canonical: true, scaled: 1);
        Assert.Multiple(() =>
        {
            Assert.That(mash.Similarity, Is.EqualTo(exact.Similarity).Within(1e-9));
            Assert.That(mash.Estimator, Is.EqualTo("MashBottomSketch"));
            Assert.That(mash.MashPValue, Is.Not.Null);
            Assert.That(mash.ContainmentAInB, Is.Null);
            Assert.That(frac.Similarity, Is.EqualTo(exact.Similarity).Within(1e-9));
            Assert.That(frac.Estimator, Is.EqualTo("FracMinHash"));
            Assert.That(frac.ContainmentAInB, Is.EqualTo(exact.ContainmentAInB).Within(1e-12));
            Assert.That(frac.SharedHashes, Is.Not.Null);
        });
    }

    [Test]
    public void CalculateSimilarity_EstimatorArguments_AreValidated()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeGenomicsTools.CalculateSimilarity("ACGTACGT", "ACGTACGA", 3, false, 10, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeGenomicsTools.CalculateSimilarity("ACGTACGT", "ACGTACGA", 3, false, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeGenomicsTools.CalculateSimilarity("ACGTACGT", "ACGTACGA", 3, false, null, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeGenomicsTools.CalculateSimilarity("ACGTACGT", "ACGTACGA", 40, false, 10));
    }
}
