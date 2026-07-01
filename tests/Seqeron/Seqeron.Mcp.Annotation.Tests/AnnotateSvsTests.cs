using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AnnotateSvsTests
{
    private static StructuralVariantDto Sv(string type, int start, int end) =>
        new("sv1", "chr1", start, end, type, end - start + 1, 60.0, 10, null);

    private static SvGeneInputDto Gene(int start, int end, params (int Start, int End)[] exons) =>
        new("GENE1", "chr1", start, end, exons.Select(e => new SvExonIntervalDto(e.Start, e.End)).ToList());

    [Test]
    public void AnnotateSvs_Schema_ValidatesCorrectly()
    {
        var vs = new List<StructuralVariantDto> { Sv("Deletion", 100, 500) };
        var gs = new List<SvGeneInputDto> { Gene(200, 400, (250, 300)) };
        Assert.DoesNotThrow(() => AnnotationTools.AnnotateSvs(vs, gs));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.AnnotateSvs(null!, gs));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.AnnotateSvs(vs, null!));
    }

    [Test]
    public void AnnotateSvs_Binding_InvokesSuccessfully()
    {
        // Deletion overlapping an exon -> HIGH impact, pathogenic (StructuralVariantAnalyzer.DetermineImpact).
        var vs = new List<StructuralVariantDto> { Sv("Deletion", 100, 500) };
        var gs = new List<SvGeneInputDto> { Gene(200, 400, (250, 300)) };

        var result = AnnotationTools.AnnotateSvs(vs, gs);

        Assert.That(result.Annotations, Has.Count.EqualTo(1));
        var a = result.Annotations[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.SvId, Is.EqualTo("sv1"));
            Assert.That(a.AffectedGenes, Is.EquivalentTo(new[] { "GENE1" }));
            Assert.That(a.AffectedExons, Is.EquivalentTo(new[] { "GENE1:exon1" }));
            Assert.That(a.FunctionalImpact, Is.EqualTo("HIGH"));
            Assert.That(a.IsPathogenic, Is.True);
        });
    }

    [Test]
    public void AnnotateSvs_GeneOverlapButNoExon_IsModifier()
    {
        // Overlaps gene [200,400] but not exon [350,380] -> MODIFIER, not pathogenic.
        var vs = new List<StructuralVariantDto> { Sv("Deletion", 100, 250) };
        var gs = new List<SvGeneInputDto> { Gene(200, 400, (350, 380)) };

        var a = AnnotationTools.AnnotateSvs(vs, gs).Annotations[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.AffectedGenes, Is.EquivalentTo(new[] { "GENE1" }));
            Assert.That(a.AffectedExons, Is.Empty);
            Assert.That(a.FunctionalImpact, Is.EqualTo("MODIFIER"));
            Assert.That(a.IsPathogenic, Is.False);
        });
    }

    [Test]
    public void AnnotateSvs_NoGeneOverlap_IsLow()
    {
        var vs = new List<StructuralVariantDto> { Sv("Deletion", 1000, 2000) };
        var gs = new List<SvGeneInputDto> { Gene(200, 400, (250, 300)) };

        var a = AnnotationTools.AnnotateSvs(vs, gs).Annotations[0];
        Assert.Multiple(() =>
        {
            Assert.That(a.AffectedGenes, Is.Empty);
            Assert.That(a.FunctionalImpact, Is.EqualTo("LOW"));
            Assert.That(a.IsPathogenic, Is.False);
        });
    }
}
