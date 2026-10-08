using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// generate_cavener_consensus — delegation to <c>MotifFinder.GenerateCavenerConsensus</c>.
/// Expected values: Biopython-locked cases of MotifFinder_GenerateConsensus_Tests (tutorial WACVC).
/// </summary>
[TestFixture]
public class GenerateCavenerConsensusTests
{
    private static readonly string[] Tutorial = { "TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC" };

    [Test]
    public void GenerateCavenerConsensus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GenerateCavenerConsensus(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateCavenerConsensus(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateCavenerConsensus(Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateCavenerConsensus(new[] { "ACG", "AC" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateCavenerConsensus(new[] { "AC-T" }));
    }

    [Test]
    public void GenerateCavenerConsensus_BiopythonTutorial_ReturnsWacvc()
    {
        Assert.That(AnalysisTools.GenerateCavenerConsensus(Tutorial).Consensus, Is.EqualTo("WACVC"));
        Assert.That(AnalysisTools.GenerateCavenerConsensus(Tutorial).Consensus,
            Is.EqualTo(global::Seqeron.Genomics.Analysis.MotifFinder.GenerateCavenerConsensus(Tutorial)));
    }
}
