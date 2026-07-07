using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.BuildKmerDatabase (Kraken canonical-k-mer -> taxon LCA DB).
// Reference values from Seqeron.Genomics.Tests MetagenomicsAnalyzer_TaxonomicClassification_Tests
// (Wood & Salzberg 2014).
[TestFixture]
public class BuildKmerDatabaseTests
{
    // Hand-built taxonomy: root(1) -> Escherichia(20) -> {E.coli(100), E.fergusonii(101)}.
    private static TaxonNodeInput[] Taxonomy() => new[]
    {
        new TaxonNodeInput(1, "root", "root", 1),
        new TaxonNodeInput(20, "Escherichia", "genus", 1),
        new TaxonNodeInput(100, "Escherichia coli", "species", 20),
        new TaxonNodeInput(101, "Escherichia fergusonii", "species", 20),
    };

    private static Dictionary<string, int> ToMap(BuildKmerDatabaseResult r)
        => r.Entries.ToDictionary(e => e.Kmer, e => e.TaxonId);

    [Test]
    public void BuildKmerDatabase_Schema_ValidatesCorrectly()
    {
        // Valid reference + taxonomy -> no throw.
        Assert.DoesNotThrow(() => MetagenomicsTools.BuildKmerDatabase(
            new[] { new ReferenceGenomeInput(100, "AAAACAA") }, Taxonomy(), k: 4));

        // A reference taxon absent from the taxonomy tree is rejected by the algorithm.
        Assert.Throws<KeyNotFoundException>(() => MetagenomicsTools.BuildKmerDatabase(
            new[] { new ReferenceGenomeInput(999, "AAAACAA") }, Taxonomy(), k: 4));
    }

    [Test]
    public void BuildKmerDatabase_Binding_InvokesSuccessfully()
    {
        // Single reference "AAAACAA" at k=4 -> 4 distinct self-canonical k-mers, all taxon 100.
        var single = MetagenomicsTools.BuildKmerDatabase(
            new[] { new ReferenceGenomeInput(100, "AAAACAA") }, Taxonomy(), k: 4);
        var singleMap = ToMap(single);

        Assert.Multiple(() =>
        {
            Assert.That(single.Count, Is.EqualTo(4), "7 - 4 + 1 = 4 distinct canonical k-mers.");
            Assert.That(singleMap.Values.All(v => v == 100), Is.True);
            Assert.That(singleMap.ContainsKey("AAAA"), Is.True);
        });

        // Shared k-mer AGCT (self-canonical) appears in both species -> collapses to LCA genus(20).
        var shared = MetagenomicsTools.BuildKmerDatabase(
            new[]
            {
                new ReferenceGenomeInput(100, "AGCTAAAA"),
                new ReferenceGenomeInput(101, "AGCTCCCC"),
            },
            Taxonomy(), k: 4);
        var sharedMap = ToMap(shared);

        Assert.Multiple(() =>
        {
            Assert.That(sharedMap["AGCT"], Is.EqualTo(20),
                "Shared k-mer -> LCA(100,101) = genus Escherichia(20) (Kraken DB build).");
            Assert.That(sharedMap["GCTA"], Is.EqualTo(100), "E.coli-only k-mer stays at species 100.");
            Assert.That(sharedMap["GAGC"], Is.EqualTo(101),
                "E.fergusonii-only k-mer (canonical of GCTC) stays at species 101.");
        });
    }
}
