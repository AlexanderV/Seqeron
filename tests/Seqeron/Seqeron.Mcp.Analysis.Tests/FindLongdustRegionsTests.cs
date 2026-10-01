using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_longdust_regions</c> MCP tool. Expected = BED output of compiled lh3/longdust 1.4-r97,
/// locked in SequenceComplexity_CalculateDustScore_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindLongdustRegionsTests
{
    private const string LdVntr =
        "ATCAGTCATTAAACTATAAACCACTTGAACCACAACGATGTCGTTTATAGCGCGCGGGGACGGCAGCTGCGATACCCCCTCGAATCCCCGGCGGCTCTCACCTGCAGGGTGGACGTTTG" +
        "GGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCA" +
        "GGTACCATTGCAGGTACCATTGCAGGTACCATTGCAGGTACCATTGCA" +
        "ACCGAGCCTCAACGGAAAGGCGGCATTGGGCGTAGATCATTGTAAGAATTGAGAGGACTGAGGGATAGGGAAAGGTACGGGCCCCGATTTCCCATGCAGGCATCTCCAAGTGTAAGCACG";

    [Test]
    public void FindLongdustRegions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindLongdustRegions(LdVntr));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLongdustRegions(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLongdustRegions(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindLongdustRegions("XYZ"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindLongdustRegions("ACGT", k: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindLongdustRegions("ACGT", threshold: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindLongdustRegions("ACGT", minStartCount: 1));
    }

    [Test]
    public void FindLongdustRegions_Binding_MatchesCompiledLongdust()
    {
        Assert.That(LdVntr, Has.Length.EqualTo(408));
        // longdust defaults (-k7 -w5000 -t0.6 -e50 -b3, both strands)
        Assert.That(AnalysisTools.FindLongdustRegions(LdVntr).Items.Select(i => (i.Start, i.End, i.Length)),
            Is.EqualTo(new[] { (120, 288, 168) }));
        // longdust -k5 -w100
        Assert.That(AnalysisTools.FindLongdustRegions(LdVntr, k: 5, windowSize: 100).Items.Select(i => (i.Start, i.End)),
            Is.EqualTo(new[] { (117, 291) }));
        Assert.That(AnalysisTools.FindLongdustRegions("NNNNNNNNNN").Items, Is.Empty);
    }
}
