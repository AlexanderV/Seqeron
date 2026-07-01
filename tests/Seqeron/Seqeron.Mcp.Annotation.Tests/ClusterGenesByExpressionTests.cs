using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ClusterGenesByExpressionTests
{
    // G1 & G3 perfectly correlate; G2 & G4 perfectly correlate (and anti-correlate with G1/G3).
    private static List<GeneProfileInputDto> Profiles() =>
        new()
        {
            new GeneProfileInputDto("G1", new List<double> { 1, 2, 3 }),
            new GeneProfileInputDto("G2", new List<double> { 3, 2, 1 }),
            new GeneProfileInputDto("G3", new List<double> { 2, 4, 6 }),
            new GeneProfileInputDto("G4", new List<double> { 6, 4, 2 }),
        };

    [Test]
    public void ClusterGenesByExpression_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ClusterGenesByExpression(Profiles(), numClusters: 2));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.ClusterGenesByExpression(null!));
    }

    [Test]
    public void ClusterGenesByExpression_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.ClusterGenesByExpression: first numClusters genes seed clusters; each
        // remaining gene joins the most-correlated cluster. G1,G2 seed; G3 joins G1's cluster, G4 joins G2's.
        var result = AnnotationTools.ClusterGenesByExpression(Profiles(), numClusters: 2);

        Assert.That(result.Clusters, Has.Count.EqualTo(2));
        var c0 = result.Clusters.Single(c => c.Genes.Contains("G1"));
        var c1 = result.Clusters.Single(c => c.Genes.Contains("G2"));
        Assert.Multiple(() =>
        {
            Assert.That(c0.Genes, Is.EquivalentTo(new[] { "G1", "G3" }));
            Assert.That(c1.Genes, Is.EquivalentTo(new[] { "G2", "G4" }));
        });
    }

    [Test]
    public void ClusterGenesByExpression_EmptyInput_ReturnsEmpty()
    {
        var result = AnnotationTools.ClusterGenesByExpression(new List<GeneProfileInputDto>());
        Assert.That(result.Clusters, Is.Empty);
    }
}
