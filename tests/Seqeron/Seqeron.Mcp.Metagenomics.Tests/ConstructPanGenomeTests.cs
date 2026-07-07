using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.ConstructPanGenome.
// Reference values from Seqeron.Genomics.Tests PanGenomeAnalyzer_ConstructPanGenome_Tests
// (Tettelin 2005/2008; Kislyuk 2011; Page 2015).
[TestFixture]
public class ConstructPanGenomeTests
{
    private const string SeqCore = "ATGCGATCGATCGATCGATCGATCGATCGA";
    private const string SeqShared2 = "TTACGGCATTACGGCATTACGGCATTACGG";
    private const string SeqU1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    private const string SeqU2 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    private const string SeqU3 = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

    private static GenomeInput Genome(string id, params (string GeneId, string Seq)[] genes)
        => new(id, genes.Select(g => new GeneInput(g.GeneId, g.Seq)).ToList());

    private static GenomeInput[] ThreeGenomes() => new[]
    {
        Genome("g1", ("a", SeqCore), ("b", SeqShared2), ("c", SeqU1)),
        Genome("g2", ("a2", SeqCore), ("b2", SeqShared2), ("d", SeqU2)),
        Genome("g3", ("a3", SeqCore), ("e", SeqU3)),
    };

    [Test]
    public void ConstructPanGenome_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.ConstructPanGenome(ThreeGenomes(), coreFraction: 1.0));

        // Empty genome set is valid (zeroed statistics), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.ConstructPanGenome(System.Array.Empty<GenomeInput>()));
    }

    [Test]
    public void ConstructPanGenome_Binding_InvokesSuccessfully()
    {
        // c1 (SeqCore) in all 3 -> core; SeqShared2 in 2 -> accessory; three singletons -> unique.
        var result = MetagenomicsTools.ConstructPanGenome(ThreeGenomes(), coreFraction: 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(result.Statistics.CoreGeneCount, Is.EqualTo(1));
            Assert.That(result.Statistics.AccessoryGeneCount, Is.EqualTo(1));
            Assert.That(result.Statistics.UniqueGeneCount, Is.EqualTo(3));
            Assert.That(result.Statistics.TotalGenes, Is.EqualTo(5),
                "1 core + 1 accessory + 3 unique = 5 clusters.");
            Assert.That(result.Statistics.CoreFraction, Is.EqualTo(1.0 / 5.0).Within(1e-10),
                "CoreFraction = 1/5 = 0.2.");
            Assert.That(result.Statistics.TotalGenomes, Is.EqualTo(3));
            // Partition lists agree with the statistics.
            Assert.That(result.CoreGenes, Has.Count.EqualTo(1));
            Assert.That(result.AccessoryGenes, Has.Count.EqualTo(1));
            Assert.That(result.UniqueGenes, Has.Count.EqualTo(3));
        });

        // Kislyuk (2011) genome fluidity on the hand-derived example: phi = 10/18.
        var fluid = MetagenomicsTools.ConstructPanGenome(new[]
        {
            Genome("A", ("a_core", SeqCore), ("a_sh", SeqShared2), ("a_u", SeqU1)),
            Genome("B", ("b_core", SeqCore), ("b_sh", SeqShared2), ("b_u", SeqU2)),
            Genome("C", ("c_core", SeqCore), ("c_u1", SeqU3),
                        ("c_u2", "TATATATATATATATATATATATATATATA")),
        }, coreFraction: 1.0);

        Assert.That(fluid.Statistics.GenomeFluidity, Is.EqualTo(10.0 / 18.0).Within(1e-10),
            "Genome fluidity = (1/3)(2/6 + 4/6 + 4/6) = 10/18 (Kislyuk 2011).");
    }
}
