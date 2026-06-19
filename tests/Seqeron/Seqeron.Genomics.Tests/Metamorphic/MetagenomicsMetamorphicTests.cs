using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Metagenomics area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-CLASS-001 — Kraken-style taxonomic read classification (Metagenomics).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 53.
///
/// API under test (MetagenomicsAnalyzer.ClassifyReads + BuildKmerDatabase):
///   Each read's canonical k-mers are looked up in a k-mer→taxon database; the best-scoring
///   root-to-leaf taxonomy path assigns the read. Reads with no hits get the root (unclassified).
///
/// Relations (derived from the per-read, hit-driven definition, NOT from output):
///   • INV (duplicate read ⇒ same classification): ClassifyRead is a pure function of the read
///          and the database, so duplicating a read yields identical classifications.
///   • MON (more reference genomes ⇒ ≥ classified reads): adding references with disjoint k-mers
///          only adds database entries (no LCA-collapse), so every read classified under the
///          smaller database stays classified and reads matching the new references become
///          classified — the number of classified reads cannot decrease.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class MetagenomicsMetamorphicTests
{
    #region Helpers

    private const int K = 11;

    // Distinct homopolymer references: canonical(poly-A) = "AAAA…" and canonical(poly-C) = "CCCC…"
    // share no canonical k-mer, so the two reference databases never collapse to an LCA.
    private const int SpeciesA = 11;
    private const int SpeciesB = 21;

    private static TaxonomyTree BuildTaxonomy() => new(new[]
    {
        new TaxonNode(TaxonomyTree.RootId, "root", "no rank", TaxonomyTree.RootId),
        new TaxonNode(SpeciesA, "SpeciesA", "species", TaxonomyTree.RootId),
        new TaxonNode(SpeciesB, "SpeciesB", "species", TaxonomyTree.RootId),
    });

    private static string Poly(char b, int n) => new(b, n);
    private static string Mixed(int n) => string.Concat(Enumerable.Repeat("ACGT", (n + 3) / 4))[..n];

    #endregion

    #region INV — duplicating a read yields identical classifications

    [Test]
    [Description("INV: classifying a read twice (a duplicated input) gives identical taxon assignments and scores, since classification is a pure per-read function.")]
    public void ClassifyReads_DuplicateRead_SameClassification()
    {
        var taxonomy = BuildTaxonomy();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (SpeciesA, Poly('A', 50)), (SpeciesB, Poly('C', 50)) }, taxonomy, K);

        foreach (var read in new[] { ("rA", Poly('A', 40)), ("rB", Poly('C', 40)), ("rC", Mixed(40)) })
        {
            var single = MetagenomicsAnalyzer.ClassifyReads(new[] { read }, db, taxonomy, K).Single();
            var duplicated = MetagenomicsAnalyzer.ClassifyReads(new[] { read, read }, db, taxonomy, K).ToList();

            duplicated.Should().HaveCount(2, because: "two input reads yield two classifications");
            foreach (var c in duplicated)
            {
                c.TaxonId.Should().Be(single.TaxonId, because: "the same read maps to the same taxon regardless of duplication");
                c.RtlScore.Should().Be(single.RtlScore, because: "the score depends only on the read's k-mer hits");
                c.MatchedKmers.Should().Be(single.MatchedKmers);
                c.TotalKmers.Should().Be(single.TotalKmers);
                c.Confidence.Should().Be(single.Confidence);
            }
        }
    }

    #endregion

    #region MON — more reference genomes classify at least as many reads

    [Test]
    [Description("MON: adding a reference genome with disjoint k-mers classifies at least as many reads — reads already classified stay classified and reads matching the new reference become classified.")]
    public void ClassifyReads_MoreReferences_ClassifyAtLeastAsManyReads()
    {
        var taxonomy = BuildTaxonomy();

        var reads = new[]
        {
            ("rA", Poly('A', 40)),   // matches SpeciesA reference
            ("rB", Poly('C', 40)),   // matches SpeciesB reference
            ("rC", Mixed(40)),       // matches no reference (control — always unclassified)
        };

        var smallDb = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (SpeciesA, Poly('A', 50)) }, taxonomy, K);
        var largeDb = MetagenomicsAnalyzer.BuildKmerDatabase(
            new[] { (SpeciesA, Poly('A', 50)), (SpeciesB, Poly('C', 50)) }, taxonomy, K);

        var small = MetagenomicsAnalyzer.ClassifyReads(reads, smallDb, taxonomy, K).ToList();
        var large = MetagenomicsAnalyzer.ClassifyReads(reads, largeDb, taxonomy, K).ToList();

        int ClassifiedCount(IEnumerable<MetagenomicsAnalyzer.TaxonomicClassification> cs) =>
            cs.Count(c => c.TaxonId != taxonomy.Root);

        ClassifiedCount(large).Should().BeGreaterThanOrEqualTo(ClassifiedCount(small),
            because: "references with disjoint k-mers only add database entries, so no previously-classified read becomes unclassified");

        // Every read classified by the smaller database is still classified (to the same taxon) by the larger.
        for (int i = 0; i < reads.Length; i++)
            if (small[i].TaxonId != taxonomy.Root)
                large[i].TaxonId.Should().Be(small[i].TaxonId,
                    because: "adding disjoint references does not change an existing read's hits");

        // The added SpeciesB reference makes the poly-C read classified, strictly increasing the count.
        ClassifiedCount(large).Should().BeGreaterThan(ClassifiedCount(small),
            because: "the poly-C read was unclassified without the SpeciesB reference and is classified once it is added");
    }

    #endregion
}
