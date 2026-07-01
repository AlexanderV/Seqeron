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
}
