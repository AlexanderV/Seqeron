using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictGenesTests
{
    // Mirrors GenomeAnnotator_Gene_Tests.CreateMinimalOrf: ATG + (aa-1)*3 A's + TAA.
    private static string CreateMinimalOrf(int aminoAcidCount = 100)
        => "ATG" + new string('A', (aminoAcidCount - 1) * 3) + "TAA";

    [Test]
    public void PredictGenes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictGenes(CreateMinimalOrf(100), minOrfLength: 50));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictGenes(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictGenes(null!));
    }

    [Test]
    public void PredictGenes_AllGenesHaveCdsTypeAndProduct()
    {
        // Mirrors PredictGenes_AllGenesHaveCdsType.
        var result = AnnotationTools.PredictGenes(CreateMinimalOrf(100), minOrfLength: 50);

        Assert.That(result.Genes, Is.Not.Empty);
        Assert.Multiple(() =>
        {
            Assert.That(result.Genes.All(g => g.Type == "CDS"), Is.True);
            Assert.That(result.Genes.All(g => g.Product == "hypothetical protein"), Is.True);
            Assert.That(result.Genes.All(g => g.Strand == '+' || g.Strand == '-'), Is.True);
        });
    }

    [Test]
    public void PredictGenes_AssignsSequentialGeneIds()
    {
        // Mirrors PredictGenes_AssignsSequentialGeneIds: two ORFs -> test_0001, test_0002, ...
        string sequence = CreateMinimalOrf(100) + new string('C', 20) + CreateMinimalOrf(100);

        var result = AnnotationTools.PredictGenes(sequence, minOrfLength: 50, prefix: "test");

        Assert.That(result.Genes.Count, Is.GreaterThanOrEqualTo(2));
        for (int i = 0; i < result.Genes.Count; i++)
            Assert.That(result.Genes[i].GeneId, Is.EqualTo($"test_{(i + 1):D4}"));
    }

    [Test]
    public void PredictGenes_ForwardOrf_HasFrameAndTranslationAttributes()
    {
        var result = AnnotationTools.PredictGenes(CreateMinimalOrf(100), minOrfLength: 50);

        var forward = result.Genes.First(g => g.Strand == '+');
        Assert.Multiple(() =>
        {
            Assert.That(forward.Start, Is.EqualTo(0));
            Assert.That(forward.End, Is.EqualTo(303)); // ATG + 297 A + TAA = 303 nt
            Assert.That(forward.Attributes["translation"], Does.StartWith("M"));
            Assert.That(forward.Attributes["translation"], Does.EndWith("*"));
        });
    }
}
