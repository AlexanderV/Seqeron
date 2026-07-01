using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindDominantIsoformsTests
{
    private static TranscriptIsoformDto Iso(string tx, string gene, double expr) =>
        new(tx, gene, 1000, 5, expr, true, new List<ExonRangeDto>());

    private static List<TranscriptIsoformDto> Isoforms() =>
        new()
        {
            Iso("TX1", "GENE1", 100),
            Iso("TX2", "GENE1", 50),
            Iso("TX3", "GENE1", 25),
        };

    [Test]
    public void FindDominantIsoforms_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindDominantIsoforms(Isoforms()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindDominantIsoforms(null!));
    }

    [Test]
    public void FindDominantIsoforms_Binding_InvokesSuccessfully()
    {
        // Mirrors Seqeron.Genomics.Tests: highest-expression isoform per gene; dominance ratio =
        // its expression / total gene expression = 100 / (100+50+25) = 100/175.
        var result = AnnotationTools.FindDominantIsoforms(Isoforms());

        Assert.That(result.Dominants, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Dominants[0].GeneId, Is.EqualTo("GENE1"));
            Assert.That(result.Dominants[0].DominantIsoform.TranscriptId, Is.EqualTo("TX1"));
            Assert.That(result.Dominants[0].DominanceRatio, Is.EqualTo(100.0 / 175.0).Within(1e-9));
        });
    }

    [Test]
    public void FindDominantIsoforms_MultipleGenes_GroupsPerGene()
    {
        var isoforms = new List<TranscriptIsoformDto>
        {
            Iso("TX1", "GENE1", 100),
            Iso("TX2", "GENE2", 200),
        };
        var result = AnnotationTools.FindDominantIsoforms(isoforms);
        Assert.That(result.Dominants.Select(d => d.GeneId), Is.EquivalentTo(new[] { "GENE1", "GENE2" }));
    }
}
