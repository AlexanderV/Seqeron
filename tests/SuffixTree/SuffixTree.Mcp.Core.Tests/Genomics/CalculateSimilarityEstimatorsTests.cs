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

    [Test]
    public void CalculateSimilarity_FracMinHash_ReportsSourmashAniEstimatesAndInterval()
    {
        // library oracle: 1 120 sourmash 4.9.4 rows (GenomicAnalyzer_AniConfidence_Oracle_Tests); here the tool wiring on a 6 kb related pair
        var rnd = new Random(2026);
        string a = new(Enumerable.Range(0, 6000).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
        string b = new(a.Select(c => rnd.NextDouble() < 0.03 ? "ACGT"[rnd.Next(4)] : c).ToArray());
        var plain = SuffixTreeGenomicsTools.CalculateSimilarity(a, b, 21, canonical: true, scaled: 2);
        var r = SuffixTreeGenomicsTools.CalculateSimilarity(a, b, 21, canonical: true, scaled: 2, estimateAniInterval: true, confidence: 0.99);
        Assert.Multiple(() =>
        {
            Assert.That(plain.AniEstimates!.ContainmentAInB.DistanceLow, Is.Null);
            var c = r.AniEstimates!.ContainmentAInB;
            Assert.That(c.DistanceLow, Is.LessThan(c.Distance));
            Assert.That(c.DistanceHigh, Is.GreaterThan(c.Distance));
            Assert.That(c.Ani, Is.EqualTo(r.ContainmentAniAInB).Within(1e-12));
            Assert.That(c.AniLow, Is.EqualTo(1.0 - c.DistanceHigh).Within(1e-15));
            Assert.That(c.SizeIsInaccurate, Is.False);
            Assert.That(r.AniEstimates.Jaccard.JaccardError, Is.Not.Null);
            Assert.That(r.AniEstimates.MaxContainment.Distance, Is.LessThanOrEqualTo(c.Distance + 1e-12));
        });
        // the exact and Mash estimators carry no sourmash ANI results
        Assert.That(SuffixTreeGenomicsTools.CalculateSimilarity("ATGC", "ATGA", 2).AniEstimates, Is.Null);
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeGenomicsTools.CalculateSimilarity(a, b, 21, scaled: 2, estimateAniInterval: true, confidence: 1.0));
    }
}
