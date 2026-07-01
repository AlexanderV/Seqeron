using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictImprintedGenesTests
{
    private static List<ImprintedGeneInputDto> Genes() =>
        new()
        {
            new ImprintedGeneInputDto("IGF2", 1000, 2000, 0.9, 0.1),  // diff 0.8 -> imprinted
            new ImprintedGeneInputDto("BALANCED", 3000, 4000, 0.5, 0.5), // diff 0 -> excluded
        };

    [Test]
    public void PredictImprintedGenes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictImprintedGenes(Genes()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.PredictImprintedGenes(null!));
    }

    [Test]
    public void PredictImprintedGenes_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.PredictImprintedGenes: |maternal-paternal| >= minDifference (default 0.4).
        // IGF2: diff=0.8 -> included. origin = maternal>paternal -> "Maternal".
        // score = diff / (maternal+paternal+0.01) = 0.8/1.01 ≈ 0.792079. HasDMR = diff>0.5 -> true.
        var result = AnnotationTools.PredictImprintedGenes(Genes());

        Assert.That(result.Imprinted, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Imprinted[0].GeneId, Is.EqualTo("IGF2"));
            Assert.That(result.Imprinted[0].Start, Is.EqualTo(1000));
            Assert.That(result.Imprinted[0].End, Is.EqualTo(2000));
            Assert.That(result.Imprinted[0].ParentalOrigin, Is.EqualTo("Maternal"));
            Assert.That(result.Imprinted[0].ImprintingScore, Is.EqualTo(0.8 / 1.01).Within(1e-9));
            Assert.That(result.Imprinted[0].HasDMR, Is.True);
        });
    }

    [Test]
    public void PredictImprintedGenes_PaternalOriginAndNoDmr()
    {
        // Paternal higher, diff exactly 0.45 (>=0.4 included, <=0.5 -> HasDMR false).
        var genes = new List<ImprintedGeneInputDto>
        {
            new ImprintedGeneInputDto("PEG", 10, 20, 0.1, 0.55),
        };
        var result = AnnotationTools.PredictImprintedGenes(genes);
        Assert.That(result.Imprinted, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Imprinted[0].ParentalOrigin, Is.EqualTo("Paternal"));
            Assert.That(result.Imprinted[0].HasDMR, Is.False); // diff 0.45 not > 0.5
        });
    }
}
