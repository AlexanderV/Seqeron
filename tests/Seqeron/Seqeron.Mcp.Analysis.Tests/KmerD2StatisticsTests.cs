using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_d2_statistics</c> MCP tool (audit round 2 WP6). Expected values: the Python replica of
/// Reinert et al. 2009 / Song et al. 2014 with a per-sequence maximum-likelihood Markov background (its formula engine
/// reproduces the CAFE binary, single strand and -R, to 6 digits) and the replica BIC (Schwarz 1978), the values locked
/// in KmerAnalyzer_ParallelAndBackgroundD2_Tests / KmerAnalyzer_BothStrandD2AndSpacedWords_Tests. NOT the wrapper output.
/// </summary>
[TestFixture]
public class KmerD2StatisticsTests
{
    private const string S1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
    private const string S2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";

    [Test]
    public void KmerD2Statistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmerD2Statistics(S1, S2, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerD2Statistics("", S2, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerD2Statistics(S1, null!, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.KmerD2Statistics(S1, S2, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.KmerD2Statistics(S1, S2, 13));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.KmerD2Statistics(S1, S2, 3, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerD2Statistics("NNNNNN", S2, 3));
    }

    [Test]
    public void KmerD2Statistics_Binding_ReturnsRawStatisticsOrdersAndBic()
    {
        var r = AnalysisTools.KmerD2Statistics(S1, S2, 3);
        Assert.Multiple(() =>
        {
            Assert.That(r.D2Star, Is.EqualTo(7.136981012975184).Within(1e-12));
            Assert.That(r.D2Shepherd, Is.EqualTo(-0.6884013403313263).Within(1e-12));
            Assert.That(r.D2StarDistance, Is.EqualTo(0.44457941706964565).Within(1e-12));
            Assert.That(r.D2ShepherdDistance, Is.EqualTo(0.5084031847096179).Within(1e-12));
            Assert.That(r.MarkovOrder1, Is.EqualTo(0));
            Assert.That(r.MarkovOrder2, Is.EqualTo(0));
            Assert.That(r.Bic1, Is.EqualTo(new[] { 231.88984329117494, 255.74056400499072, 370.90681674887685 }).Within(1e-9));
            Assert.That(r.Bic2, Is.EqualTo(new[] { 199.825085479026, 217.12856241056033, 322.30742386881116 }).Within(1e-9));
        });

        var both = AnalysisTools.KmerD2Statistics(S1, S2, 3, 1, bothStrands: true);
        Assert.Multiple(() =>
        {
            Assert.That(both.D2Star, Is.EqualTo(-5.014593950349383).Within(1e-12));
            Assert.That(both.D2Shepherd, Is.EqualTo(-1.3145456673597726).Within(1e-12));
            Assert.That(both.D2StarDistance, Is.EqualTo(0.5716659331844273).Within(1e-12));
            Assert.That(both.D2ShepherdDistance, Is.EqualTo(0.5152519131905425).Within(1e-12));
            Assert.That(both.MarkovOrder1, Is.EqualTo(1));
        });

        // markovOrder = -1: the reported orders are the BIC minimisers.
        var auto = AnalysisTools.KmerD2Statistics(S1, S2, 4, -1);
        Assert.Multiple(() =>
        {
            Assert.That(auto.MarkovOrder1, Is.EqualTo(Array.IndexOf(auto.Bic1, auto.Bic1.Min())));
            Assert.That(auto.MarkovOrder2, Is.EqualTo(Array.IndexOf(auto.Bic2, auto.Bic2.Min())));
            Assert.That(auto.Bic1, Has.Length.EqualTo(4));
        });
    }
}
