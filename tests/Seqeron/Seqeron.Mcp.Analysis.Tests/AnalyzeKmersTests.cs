using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>analyze_kmers</c> MCP tool.
/// Expected values taken from KMER-STATS-001 / the algorithm's own unit tests
/// (KmerAnalyzer_AnalyzeKmers_Tests, Wikipedia K-mer worked example GTAGAGCTGT),
/// NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class AnalyzeKmersTests
{
    [Test]
    public void AnalyzeKmers_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.AnalyzeKmers("GTAGAGCTGT", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeKmers("", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeKmers(null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeKmers("GTAG", 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeKmers("GTAG", -1));
    }

    [Test]
    public void AnalyzeKmers_Binding_InvokesSuccessfully()
    {
        // GTAGAGCTGT, k=1: G4 T3 A2 C1 -> total 10, distinct 4, max 4, min 1, avg 2.5.
        // Entropy of {0.4,0.3,0.2,0.1} = 1.846439344671 bits.
        var s1 = AnalysisTools.AnalyzeKmers("GTAGAGCTGT", 1);
        Assert.Multiple(() =>
        {
            Assert.That(s1.TotalKmers, Is.EqualTo(10));
            Assert.That(s1.UniqueKmers, Is.EqualTo(4));
            Assert.That(s1.MaxCount, Is.EqualTo(4));
            Assert.That(s1.MinCount, Is.EqualTo(1));
            Assert.That(s1.AverageCount, Is.EqualTo(2.5).Within(1e-10));
            Assert.That(s1.Entropy, Is.EqualTo(1.846439344671).Within(1e-10));
        });
    }

    [Test]
    public void AnalyzeKmers_AllDistinctTrimers_EntropyIsLog2Eight()
    {
        // GTAGAGCTGT, k=3: 8 windows all distinct -> entropy = log2(8) = 3 exactly.
        var s3 = AnalysisTools.AnalyzeKmers("GTAGAGCTGT", 3);
        Assert.Multiple(() =>
        {
            Assert.That(s3.TotalKmers, Is.EqualTo(8));
            Assert.That(s3.UniqueKmers, Is.EqualTo(8));
            Assert.That(s3.MaxCount, Is.EqualTo(1));
            Assert.That(s3.MinCount, Is.EqualTo(1));
            Assert.That(s3.AverageCount, Is.EqualTo(1.0).Within(1e-10));
            Assert.That(s3.Entropy, Is.EqualTo(3.0).Within(1e-10));
        });
    }

    [Test]
    public void AnalyzeKmers_JellyfishFieldsAndCountFilters()
    {
        // Jellyfish stats compute_stats replica (Python Counter) + scipy entropy(base=2):
        // GTAGAGCTGT k=2: Unique 5, Distinct 7, Total 9, Max 2; with -L 2 -U 2: 0, 2, 4, 2, entropy 1.
        var all = AnalysisTools.AnalyzeKmers("GTAGAGCTGT", 2);
        var filt = AnalysisTools.AnalyzeKmers("GTAGAGCTGT", 2, lowerCount: 2, upperCount: 2);
        Assert.Multiple(() =>
        {
            Assert.That(all.SingletonKmers, Is.EqualTo(5));
            Assert.That(all.DistinctKmers, Is.EqualTo(7));
            Assert.That(all.AverageCount, Is.EqualTo(9.0 / 7.0).Within(1e-15));
            Assert.That(filt.SingletonKmers, Is.EqualTo(0));
            Assert.That(filt.DistinctKmers, Is.EqualTo(2));
            Assert.That(filt.TotalKmers, Is.EqualTo(4));
            Assert.That(filt.MaxCount, Is.EqualTo(2));
            Assert.That(filt.Entropy, Is.EqualTo(1.0).Within(1e-12));
        });
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeKmers("GTAG", 2, lowerCount: -1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AnalyzeKmers("GTAG", 2, upperCount: -1));
    }

    /// <summary>
    /// Audit round 1 WP4: optional canonical / acgtOnly delegate to the option-aware library overload.
    /// Reference: Jellyfish 2.3.1 binary, BA1B sample k=4: count -C + stats → Unique 16, Distinct 20, Total 27, Max 4;
    /// count -C + stats -L 2 → 0, 4, 11, 4.
    /// </summary>
    [Test]
    public void AnalyzeKmers_CanonicalMode_MatchesJellyfishStats()
    {
        const string ba1b = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
        var c = AnalysisTools.AnalyzeKmers(ba1b, 4, canonical: true);
        Assert.Multiple(() =>
        {
            Assert.That(c.SingletonKmers, Is.EqualTo(16));
            Assert.That(c.DistinctKmers, Is.EqualTo(20));
            Assert.That(c.TotalKmers, Is.EqualTo(27));
            Assert.That(c.MaxCount, Is.EqualTo(4));
        });
        var l2 = AnalysisTools.AnalyzeKmers(ba1b, 4, lowerCount: 2, canonical: true);
        Assert.That((l2.SingletonKmers, l2.DistinctKmers, l2.TotalKmers, l2.MaxCount), Is.EqualTo((0, 4, 11, 4)));
        // acgtOnly drops the N windows.
        Assert.That(AnalysisTools.AnalyzeKmers("ACGTNACGT", 4, acgtOnly: true).TotalKmers, Is.EqualTo(2));
    }
}
