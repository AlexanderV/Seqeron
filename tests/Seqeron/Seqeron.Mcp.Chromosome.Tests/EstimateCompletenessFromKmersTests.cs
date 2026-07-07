using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>estimate_completeness_from_kmers</c>. GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers
/// derives completeness, error rate and genome size from a k-mer count spectrum. A clean spectrum with
/// all k-mers at the expected coverage has no error k-mers -> completeness 1, error 0.
/// </summary>
[TestFixture]
public class EstimateCompletenessFromKmersTests
{
    private static List<KmerCount> Spectrum(int distinct, int coverage)
    {
        var s = new List<KmerCount>();
        for (int i = 0; i < distinct; i++) s.Add(new KmerCount("k" + i, coverage));
        return s;
    }

    [Test]
    public void EstimateCompletenessFromKmers_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.EstimateCompletenessFromKmers(Spectrum(100, 30), 30));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCompletenessFromKmers(null!, 30));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCompletenessFromKmers(Spectrum(10, 30), expectedCoverage: -1));
    }

    [Test]
    public void EstimateCompletenessFromKmers_CleanSpectrum_FullCompleteness()
    {
        // 100 distinct k-mers all at coverage 30 (== expected): no error k-mers.
        var result = ChromosomeTools.EstimateCompletenessFromKmers(Spectrum(100, 30), 30);

        Assert.Multiple(() =>
        {
            Assert.That(result.Completeness, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(result.ErrorRate, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(result.EstimatedGenomeSize, Is.EqualTo(100));
        });
    }
}
