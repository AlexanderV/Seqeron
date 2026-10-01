using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_distance</c> MCP tool.
/// Expected values derived from the Euclidean word-composition distance definition
/// (identity => 0; orthogonal monomer vectors => sqrt(2)), NOT the wrapper output.
/// </summary>
[TestFixture]
public class KmerDistanceTests
{
    [Test]
    public void KmerDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmerDistance("AAAA", "TTTT", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("", "TTTT", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("AAAA", "", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance(null!, "TTTT", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("AAAA", null!, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("AAAA", "TTTT", 0));
    }

    [Test]
    public void KmerDistance_Binding_InvokesSuccessfully()
    {
        // Identical sequences => distance 0.
        var same = AnalysisTools.KmerDistance("ACGTACGT", "ACGTACGT", 2).Distance;
        Assert.That(same, Is.EqualTo(0.0).Within(1e-12));

        // {A:1} vs {T:1} => sqrt(1^2 + 1^2) = sqrt(2).
        var orth = AnalysisTools.KmerDistance("AAAA", "TTTT", 1).Distance;
        Assert.That(orth, Is.EqualTo(Math.Sqrt(2)).Within(1e-12));
    }

    [Test]
    public void KmerDistance_Metric_DelegatesToMetricOverload()
    {
        // Zielezinski 2017 Fig. 1 (ATGTGTG / CATGTG, k=3); reference values scipy.spatial.distance on
        // scikit-bio kmer_frequencies (counts / relative), alfpy 1.0.6 word_distance.
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3).Distance, Is.EqualTo(Math.Sqrt(0.11)).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "euclidean").Distance, Is.EqualTo(Math.Sqrt(0.11)).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "squared_euclidean_counts").Distance, Is.EqualTo(3.0));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "manhattan").Distance, Is.EqualTo(0.6).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "chebyshev").Distance, Is.EqualTo(0.25).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "canberra").Distance, Is.EqualTo(1.5726495726495728).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "COSINE").Distance, Is.EqualTo(1.0 / 6).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "d2").Distance, Is.EqualTo(5.0));
            Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("ACGT", "ACGT", 2, "hamming"));
        });
    }
}
