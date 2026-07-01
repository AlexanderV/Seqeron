using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class OverRepresentationAnalysisTests
{
    private static List<string> DeGenes() => new() { "A", "B", "C", "D", "E" };

    private static List<PathwayInputDto> Pathways() =>
        new() { new PathwayInputDto("P1", "Pathway1", new List<string> { "A", "B", "C", "X", "Y" }) };

    [Test]
    public void OverRepresentationAnalysis_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.OverRepresentationAnalysis(DeGenes(), Pathways(), 1000));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.OverRepresentationAnalysis(null!, Pathways(), 1000));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.OverRepresentationAnalysis(DeGenes(), null!, 1000));
        Assert.Throws<ArgumentException>(() => AnnotationTools.OverRepresentationAnalysis(new List<string>(), Pathways(), 1000));
        Assert.Throws<ArgumentException>(() => AnnotationTools.OverRepresentationAnalysis(DeGenes(), Pathways(), 0));
    }

    [Test]
    public void OverRepresentationAnalysis_Binding_InvokesSuccessfully()
    {
        // Mirrors Seqeron.Genomics.Tests: DE {A..E} vs pathway {A,B,C,X,Y} over 1000 genes ->
        // 3 overlapping, enrichment score > 1.
        var result = AnnotationTools.OverRepresentationAnalysis(DeGenes(), Pathways(), 1000);

        Assert.That(result.Results, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Results[0].PathwayId, Is.EqualTo("P1"));
            Assert.That(result.Results[0].GenesInPathway, Is.EqualTo(5));
            Assert.That(result.Results[0].OverlappingGenes, Is.EqualTo(3));
            Assert.That(result.Results[0].EnrichmentScore, Is.GreaterThan(1.0));
            Assert.That(result.Results[0].Genes, Is.EquivalentTo(new[] { "A", "B", "C" }));
        });
    }

    [Test]
    public void OverRepresentationAnalysis_NoOverlap_NotIncluded()
    {
        var pathways = new List<PathwayInputDto>
        {
            new PathwayInputDto("P1", "Pathway1", new List<string> { "X", "Y", "Z" }),
        };
        var result = AnnotationTools.OverRepresentationAnalysis(new List<string> { "A", "B", "C" }, pathways, 1000);
        Assert.That(result.Results, Is.Empty);
    }
}
