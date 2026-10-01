using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>generate_consensus</c> MCP tool.
/// Expected values from the IUPAC >25% inclusion rule ({A,T} -> W), NOT the wrapper output.
/// </summary>
[TestFixture]
public class GenerateConsensusTests
{
    [Test]
    public void GenerateConsensus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GenerateConsensus(new[] { "ATGC" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateConsensus(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateConsensus(Array.Empty<string>()));
    }

    [Test]
    public void GenerateConsensus_Binding_InvokesSuccessfully()
    {
        // Unanimous columns -> the single base at each position.
        var same = AnalysisTools.GenerateConsensus(new[] { "ATGC", "ATGC", "ATGC" }).Consensus;
        Assert.That(same, Is.EqualTo("ATGC"));

        // A and T each 50% (>25%) -> IUPAC W at every column.
        var amb = AnalysisTools.GenerateConsensus(new[] { "AAAA", "TTTT" }).Consensus;
        Assert.That(amb, Is.EqualTo("WWWW"));
    }

    [Test]
    public void GenerateConsensus_FourWayTie_ReturnsN_UnequalLengthsRejected()
    {
        // Four equally abundant bases → N (Biopython degenerate_consensus / DECIPHER equal-abundance rule).
        Assert.That(AnalysisTools.GenerateConsensus(new[] { "A", "C", "G", "T" }).Consensus, Is.EqualTo("N"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateConsensus(new[] { "ACG", "AC" }));
    }

    [Test]
    public void GenerateConsensus_InclusionThreshold_DelegatesToOverload()
    {
        string[] rows = { "AAAA", "TTTT", "TTTT", "CCCC" };
        // Default 0.25: only T (50 %) is strictly above 25 % → T.
        Assert.That(AnalysisTools.GenerateConsensus(rows).Consensus, Is.EqualTo("TTTT"));
        // 0.2: A, C (25 %) and T pass → H (A/C/T).
        Assert.That(AnalysisTools.GenerateConsensus(rows, 0.2).Consensus, Is.EqualTo("HHHH"));
        Assert.That(AnalysisTools.GenerateConsensus(rows, 0.2).Consensus,
            Is.EqualTo(global::Seqeron.Genomics.Analysis.MotifFinder.GenerateConsensus(rows, 0.2)));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.GenerateConsensus(rows, 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.GenerateConsensus(rows, double.NaN));
    }
}
