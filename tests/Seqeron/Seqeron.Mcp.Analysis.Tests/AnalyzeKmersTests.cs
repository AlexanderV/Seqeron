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
}
