using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.SelectPhylogeneticMarkers: single-copy core clusters with >= 1
// parsimony-informative site, ranked by descending PIS, capped at maxMarkers (panX/Roary).
// Reference values from Seqeron.Genomics.Tests PanGenomeAnalyzer_SelectPhylogeneticMarkers_Tests
// (Ding 2018 panX; Page 2015 Roary; Zvelebil 2008).
[TestFixture]
public class SelectPhylogeneticMarkersTests
{
    private static GenomeInput Genome(string id, params (string GeneId, string Seq)[] genes)
        => new(id, genes.Select(g => new GeneInput(g.GeneId, g.Seq)).ToList());

    private static PanGenomeAnalyzer.GeneCluster Cluster(string id, string[] geneIds, string[] genomeIds)
        => new(id, geneIds, genomeIds, genomeIds.Length, 1.0, geneIds.Length > 0 ? geneIds[0] : "");

    // aHi members AC,AC,GG,GG -> 2 PI columns; zLo members AC,AC,GC,GT -> 1 PI column.
    private static GenomeInput[] Genomes() => new[]
    {
        Genome("g1", ("hi1", "AC"), ("lo1", "AC")),
        Genome("g2", ("hi2", "AC"), ("lo2", "AC")),
        Genome("g3", ("hi3", "GG"), ("lo3", "GC")),
        Genome("g4", ("hi4", "GG"), ("lo4", "GT")),
    };

    private static PanGenomeAnalyzer.GeneCluster[] Candidates() => new[]
    {
        Cluster("aHi", new[] { "hi1", "hi2", "hi3", "hi4" }, new[] { "g1", "g2", "g3", "g4" }),
        Cluster("zLo", new[] { "lo1", "lo2", "lo3", "lo4" }, new[] { "g1", "g2", "g3", "g4" }),
    };

    [Test]
    public void SelectPhylogeneticMarkers_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.SelectPhylogeneticMarkers(
            Genomes(), Candidates(), totalGenomes: 4));

        // Empty candidate set is defined (no markers), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.SelectPhylogeneticMarkers(
            Genomes(), System.Array.Empty<PanGenomeAnalyzer.GeneCluster>(), totalGenomes: 4));
    }

    [Test]
    public void SelectPhylogeneticMarkers_Binding_InvokesSuccessfully()
    {
        // Both are single-copy core; ranked by descending PIS -> aHi (2) before zLo (1).
        var all = MetagenomicsTools.SelectPhylogeneticMarkers(Genomes(), Candidates(), totalGenomes: 4);
        Assert.That(all.Items.Select(m => m.ClusterId).ToArray(), Is.EqualTo(new[] { "aHi", "zLo" }),
            "Markers ordered by descending parsimony-informative-site count.");

        // maxMarkers = 1 keeps only the most informative marker.
        var capped = MetagenomicsTools.SelectPhylogeneticMarkers(
            Genomes(), Candidates(), totalGenomes: 4, maxMarkers: 1);
        Assert.That(capped.Items.Select(m => m.ClusterId).ToArray(), Is.EqualTo(new[] { "aHi" }));
    }
}
