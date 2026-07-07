using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.FindResistanceGenes (motif-containment search).
// A gene hits when its sequence contains the DB motif; Identity = motifLength / geneLength.
// Reference values hand-derived from the algorithm contract.
[TestFixture]
public class FindResistanceGenesTests
{
    [Test]
    public void FindResistanceGenes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.FindResistanceGenes(
            new[] { new GeneInput("g1", "AAACGTACGT") },
            new[] { new ResistanceDatabaseEntry("CGTACGT", "blaX-like", "beta-lactam") }));

        // Empty gene set is defined input (no hits), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.FindResistanceGenes(
            System.Array.Empty<GeneInput>(),
            new[] { new ResistanceDatabaseEntry("CGTACGT", "blaX-like", "beta-lactam") }));
    }

    [Test]
    public void FindResistanceGenes_Binding_InvokesSuccessfully()
    {
        // Gene "AAACGTACGT" (len 10) contains motif "CGTACGT" (len 7) -> hit, identity = 7/10.
        var hits = MetagenomicsTools.FindResistanceGenes(
            new[] { new GeneInput("g1", "AAACGTACGT") },
            new[] { new ResistanceDatabaseEntry("CGTACGT", "blaX-like", "beta-lactam") });

        var hit = hits.Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That(hit.GeneId, Is.EqualTo("g1"));
            Assert.That(hit.ResistanceGene, Is.EqualTo("blaX-like"));
            Assert.That(hit.AntibioticClass, Is.EqualTo("beta-lactam"));
            Assert.That(hit.Identity, Is.EqualTo(7.0 / 10.0).Within(1e-10),
                "Identity = motif length (7) / gene length (10) = 0.7.");
        });

        // A gene not containing the motif yields no hit.
        var noHit = MetagenomicsTools.FindResistanceGenes(
            new[] { new GeneInput("g2", "AAAAAAAAAA") },
            new[] { new ResistanceDatabaseEntry("CGTACGT", "blaX-like", "beta-lactam") });
        Assert.That(noHit.Items, Is.Empty, "Gene without the motif -> no resistance hit.");
    }
}
