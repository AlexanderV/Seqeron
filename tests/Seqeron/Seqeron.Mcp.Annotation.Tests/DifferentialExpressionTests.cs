using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class DifferentialExpressionTests
{
    // Mirrors Seqeron.Genomics.Tests TranscriptomeAnalyzerTests.AnalyzeDifferentialExpression_* :
    // GENE1 group1 {10,12,11,13} -> group2 {100,110,105,115} is a clear ~10x upregulation.
    private static List<DiffExprInputDto> Upregulated() =>
        new()
        {
            new DiffExprInputDto("GENE1",
                new List<double> { 10, 12, 11, 13 },
                new List<double> { 100, 110, 105, 115 }),
        };

    [Test]
    public void DifferentialExpression_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.DifferentialExpression(Upregulated()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.DifferentialExpression(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.DifferentialExpression(new List<DiffExprInputDto>()));
    }

    [Test]
    public void DifferentialExpression_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.AnalyzeDifferentialExpression: log2 fold change > 2 (~10x),
        // Regulation "Upregulated", significant (|log2FC| >= 1 and adjusted p < 0.05).
        var result = AnnotationTools.DifferentialExpression(Upregulated());

        Assert.That(result.Results, Has.Count.EqualTo(1));
        var g = result.Results[0];
        Assert.Multiple(() =>
        {
            Assert.That(g.GeneId, Is.EqualTo("GENE1"));
            Assert.That(g.Log2FoldChange, Is.GreaterThan(2.0));
            Assert.That(g.Regulation, Is.EqualTo("Upregulated"));
            Assert.That(g.IsSignificant, Is.True);
            // Single gene: BH-adjusted p equals the raw p (m/rank = 1).
            Assert.That(g.AdjustedPValue, Is.EqualTo(g.PValue).Within(1e-12));
        });
    }

    [Test]
    public void DifferentialExpression_Downregulated_AnnotatesCorrectly()
    {
        // Reversed groups -> negative log2 fold change, "Downregulated".
        var data = new List<DiffExprInputDto>
        {
            new DiffExprInputDto("GENE1",
                new List<double> { 100, 110, 105 },
                new List<double> { 10, 12, 11 }),
        };
        var result = AnnotationTools.DifferentialExpression(data);
        Assert.Multiple(() =>
        {
            Assert.That(result.Results[0].Log2FoldChange, Is.LessThan(0));
            Assert.That(result.Results[0].Regulation, Is.EqualTo("Downregulated"));
        });
    }

    [Test]
    public void DifferentialExpression_NoChange_NotSignificant()
    {
        var data = new List<DiffExprInputDto>
        {
            new DiffExprInputDto("GENE1",
                new List<double> { 100, 102, 98, 101 },
                new List<double> { 101, 99, 100, 102 }),
        };
        var result = AnnotationTools.DifferentialExpression(data);
        Assert.That(result.Results[0].IsSignificant, Is.False);
    }
}
