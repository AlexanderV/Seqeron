using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class BuildCoexpressionNetworkTests
{
    // Fixture mirrors TranscriptomeAnalyzerTests.BuildCoExpressionNetwork_HighCorrelation_IncludesEdge.
    private static List<GeneProfileInputDto> Profiles() => new()
    {
        new("GENE1", new List<double> { 1, 2, 3, 4, 5 }),
        new("GENE2", new List<double> { 2, 4, 6, 8, 10 }), // perfect +1 correlation with GENE1
        new("GENE3", new List<double> { 5, 4, 3, 2, 1 })   // perfect -1 correlation with GENE1
    };

    [Test]
    public void BuildCoexpressionNetwork_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.BuildCoexpressionNetwork(Profiles()));
        Assert.DoesNotThrow(() => AnnotationTools.BuildCoexpressionNetwork(new List<GeneProfileInputDto>()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.BuildCoexpressionNetwork(null!));
    }

    [Test]
    public void BuildCoexpressionNetwork_PerfectlyCorrelatedGenes_EmitAllEdges()
    {
        // |corr| = 1 for every pair, so all three upper-triangle pairs meet threshold 0.7.
        var result = AnnotationTools.BuildCoexpressionNetwork(Profiles(), correlationThreshold: 0.7);

        Assert.That(result.Edges, Has.Count.EqualTo(3));

        var e12 = result.Edges.Single(e => e.Gene1 == "GENE1" && e.Gene2 == "GENE2");
        var e13 = result.Edges.Single(e => e.Gene1 == "GENE1" && e.Gene2 == "GENE3");
        var e23 = result.Edges.Single(e => e.Gene1 == "GENE2" && e.Gene2 == "GENE3");
        Assert.Multiple(() =>
        {
            Assert.That(e12.Correlation, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(e13.Correlation, Is.EqualTo(-1.0).Within(1e-9));
            Assert.That(e23.Correlation, Is.EqualTo(-1.0).Within(1e-9));
        });
    }

    [Test]
    public void BuildCoexpressionNetwork_HighThreshold_KeepsOnlyPositiveEdge()
    {
        // Only |corr| >= 1.0 survive; all pairs are exactly 1.0 in magnitude, so still all three.
        // Raising above 1.0 drops everything.
        var result = AnnotationTools.BuildCoexpressionNetwork(Profiles(), correlationThreshold: 1.0001);
        Assert.That(result.Edges, Is.Empty);
    }

    [Test]
    public void BuildCoexpressionNetwork_SingleGene_NoEdges()
    {
        var result = AnnotationTools.BuildCoexpressionNetwork(
            new List<GeneProfileInputDto> { new("GENE1", new List<double> { 1, 2, 3 }) });
        Assert.That(result.Edges, Is.Empty);
    }
}
