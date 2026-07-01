using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.ClassifyReads (Kraken k-mer/LCA RTL classification).
// Reference values from Seqeron.Genomics.Tests MetagenomicsAnalyzer_TaxonomicClassification_Tests
// (Wood & Salzberg 2014). Hand-built taxonomy, k = 4, self-canonical k-mers.
[TestFixture]
public class ClassifyReadsTests
{
    // root(1) -> Escherichia(20) -> {E.coli(100), E.fergusonii(101)}.
    private static TaxonNodeInput[] Taxonomy() => new[]
    {
        new TaxonNodeInput(1, "root", "root", 1),
        new TaxonNodeInput(20, "Escherichia", "genus", 1),
        new TaxonNodeInput(100, "Escherichia coli", "species", 20),
        new TaxonNodeInput(101, "Escherichia fergusonii", "species", 20),
    };

    private static KmerDatabaseEntry[] Db(params (string Kmer, int Taxon)[] entries)
        => entries.Select(e => new KmerDatabaseEntry(e.Kmer, e.Taxon)).ToArray();

    [Test]
    public void ClassifyReads_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.ClassifyReads(
            new[] { new ReadInput("r", "AAAACAA") },
            Db(("AAAA", 100), ("AAAC", 100), ("AACA", 100), ("ACAA", 100)),
            Taxonomy(), k: 4));

        // k must be positive (algorithm guard).
        Assert.Throws<System.ArgumentOutOfRangeException>(() => MetagenomicsTools.ClassifyReads(
            new[] { new ReadInput("r", "AAAACAA") },
            Db(("AAAA", 100)), Taxonomy(), k: 0));
    }

    [Test]
    public void ClassifyReads_Binding_InvokesSuccessfully()
    {
        var taxonomy = Taxonomy();

        // Single species: all 4 windows of "AAAACAA" map to E.coli(100).
        var single = MetagenomicsTools.ClassifyReads(
            new[] { new ReadInput("r", "AAAACAA") },
            Db(("AAAA", 100), ("AAAC", 100), ("AACA", 100), ("ACAA", 100)),
            taxonomy, k: 4).Items.Single();

        Assert.Multiple(() =>
        {
            Assert.That(single.TaxonId, Is.EqualTo(100));
            Assert.That(single.TaxonName, Is.EqualTo("Escherichia coli"));
            Assert.That(single.Rank, Is.EqualTo("species"));
            Assert.That(single.RtlScore, Is.EqualTo(4), "RTL score = sum of node weights on root->100.");
            Assert.That(single.TotalKmers, Is.EqualTo(4), "Q = 4 non-ambiguous k-mers.");
            Assert.That(single.MatchedKmers, Is.EqualTo(4));
            Assert.That(single.Confidence, Is.EqualTo(1.0).Within(1e-9));
        });

        // Split within genus: 2 windows -> 100, 2 -> 101; tie -> LCA = genus Escherichia(20).
        var split = MetagenomicsTools.ClassifyReads(
            new[] { new ReadInput("r", "AAAACAA") },
            Db(("AAAA", 100), ("AAAC", 100), ("AACA", 101), ("ACAA", 101)),
            taxonomy, k: 4).Items.Single();

        Assert.Multiple(() =>
        {
            Assert.That(split.TaxonId, Is.EqualTo(20), "Tied species leaves -> LCA genus(20).");
            Assert.That(split.Rank, Is.EqualTo("genus"));
            Assert.That(split.RtlScore, Is.EqualTo(2), "Each tied RTL path scores 2.");
            Assert.That(split.MatchedKmers, Is.EqualTo(4), "C = clade(20) = both species' 4 k-mers.");
        });
    }
}
