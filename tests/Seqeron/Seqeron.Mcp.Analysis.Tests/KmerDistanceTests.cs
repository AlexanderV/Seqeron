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

    /// <summary>
    /// Audit round 1 WP4: d2star / d2shepherd (+ markovOrder) delegate to KmerAnalyzer.BackgroundAdjustedD2.
    /// Reference values: Python replica of Reinert et al. 2009 / Song et al. 2014 (formula engine checked against the
    /// CAFE binary), pair S1/S2 of KmerAnalyzer_ParallelAndBackgroundD2_Tests.
    /// </summary>
    [Test]
    public void KmerDistance_D2StarMetrics_DelegateToBackgroundAdjustedD2()
    {
        const string s1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
        const string s2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 3, "d2star").Distance, Is.EqualTo(0.44457941706964565).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 3, "d2shepherd", 1).Distance, Is.EqualTo(0.5226793773737308).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 3, "d2s", 1).Distance, Is.EqualTo(0.5226793773737308).Within(1e-12));
        });
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance(s1, s2, 3, "euclidean", 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.KmerDistance(s1, s2, 3, "d2star", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance(s1, s2, 3, "bogus"));
    }

    /// <summary>
    /// Audit round 2 WP6: bothStrands (CAFE -R) and the jensen_shannon / euclidean_counts metrics. Reference values:
    /// Python replica whose formula engine reproduces the CAFE -R binary (20 runs, 6 digits), sequence MLE background;
    /// scipy jensenshannon(base=2)**2.
    /// </summary>
    [Test]
    public void KmerDistance_BothStrandsAndNewMetrics_Delegate()
    {
        const string s1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
        const string s2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 3, "d2star", 0, bothStrands: true).Distance, Is.EqualTo(0.3838876158581438).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 3, "d2shepherd", 1, bothStrands: true).Distance, Is.EqualTo(0.5152519131905425).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 3, "d2star", 0, bothStrands: false).Distance, Is.EqualTo(0.44457941706964565).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "jensen_shannon").Distance, Is.EqualTo(0.1522040934665307).Within(1e-12));
            Assert.That(AnalysisTools.KmerDistance("ATGTGTG", "CATGTG", 3, "euclidean_counts").Distance, Is.EqualTo(Math.Sqrt(3)).Within(1e-12));
        });
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance(s1, s2, 3, "euclidean", 0, bothStrands: true));
    }

    /// <summary>
    /// Audit round 3 WP9: metric ev = spaced -d EV with the pattern 1^k (spaced 1.2.0 -f 1111 on S1/S2 prints 1.2, the
    /// saturation value, with and without -r).
    /// </summary>
    [Test]
    public void KmerDistance_Ev_DelegatesToSpacedEvolutionary()
    {
        const string s1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
        const string s2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 4, "ev").Distance, Is.EqualTo(1.2));
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 4, "ev", 0, bothStrands: true).Distance, Is.EqualTo(1.2));
            Assert.That(AnalysisTools.KmerDistance(s1, s2, 3, "ev").Distance, Is.EqualTo(
                Seqeron.Genomics.Analysis.KmerAnalyzer.SpacedWordDistance(s1, s2, new[] { "111" }, Seqeron.Genomics.Analysis.KmerDistanceMetric.SpacedEvolutionary)));
        });
    }
}
