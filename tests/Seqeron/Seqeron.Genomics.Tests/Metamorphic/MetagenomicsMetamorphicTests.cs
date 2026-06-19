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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: META-PROF-001 — taxonomic profile / relative abundance (Metagenomics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 54.
    //
    // API under test (MetagenomicsAnalyzer.GenerateTaxonomicProfile):
    //   Counts the classified reads per taxon at each rank and divides by the number of
    //   classified reads to give relative abundances; also reports Shannon/Simpson diversity
    //   over the species abundances.
    //
    // Relations (derived from the count/total definition, NOT from output):
    //   • INV (doubling all reads ⇒ same abundances): duplicating every classification doubles
    //          every count and the total, so the relative abundances — and the diversity indices
    //          computed from them — are unchanged; the profile is also independent of read order.
    //   • COMP (abundances sum to 1): at any rank where every classified read carries a value,
    //          the per-taxon counts partition the classified reads, so the abundances sum to 1.
    // ───────────────────────────────────────────────────────────────────────────

    #region Profile helpers

    private static MetagenomicsAnalyzer.TaxonomicClassification Classified(
        string kingdom, string phylum, string genus, string species) =>
        new(ReadId: "r", TaxonId: 2, TaxonName: species, Rank: "species",
            RtlScore: 1, Confidence: 1.0, MatchedKmers: 1, TotalKmers: 1,
            Kingdom: kingdom, Phylum: phylum, Class: "", Order: "", Family: "",
            Genus: genus, Species: species);

    private static List<MetagenomicsAnalyzer.TaxonomicClassification> SampleClassifications() => new()
    {
        Classified("Bacteria", "Firmicutes", "Bacillus", "B.subtilis"),
        Classified("Bacteria", "Firmicutes", "Bacillus", "B.subtilis"),
        Classified("Bacteria", "Firmicutes", "Bacillus", "B.subtilis"),
        Classified("Bacteria", "Firmicutes", "Bacillus", "B.cereus"),
        Classified("Bacteria", "Firmicutes", "Bacillus", "B.cereus"),
        Classified("Bacteria", "Proteobacteria", "Escherichia", "E.coli"),
        // An unclassified read (excluded from abundances, must not perturb the relative values).
        Classified("Unclassified", "", "", ""),
    };

    private static void AbundancesShouldMatch(
        IReadOnlyDictionary<string, double> actual, IReadOnlyDictionary<string, double> expected, string rank)
    {
        actual.Keys.Should().BeEquivalentTo(expected.Keys, because: $"the {rank} taxa present are unchanged");
        foreach (var (taxon, value) in expected)
            actual[taxon].Should().BeApproximately(value, 1e-12,
                because: $"the relative abundance of {taxon} at {rank} level is unchanged");
    }

    #endregion

    #region INV — doubling reads (and reordering) preserves the relative abundances

    [Test]
    [Description("INV: duplicating every classification doubles all counts and the total, leaving the relative abundances and the diversity indices unchanged; the profile is also order-independent.")]
    public void GenerateTaxonomicProfile_DoublingReads_PreservesAbundances()
    {
        var baseList = SampleClassifications();
        var baseProfile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(baseList);

        var doubled = baseList.Concat(baseList).ToList();
        var doubledProfile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(doubled);

        AbundancesShouldMatch(doubledProfile.SpeciesAbundance, baseProfile.SpeciesAbundance, "species");
        AbundancesShouldMatch(doubledProfile.GenusAbundance, baseProfile.GenusAbundance, "genus");
        AbundancesShouldMatch(doubledProfile.PhylumAbundance, baseProfile.PhylumAbundance, "phylum");
        AbundancesShouldMatch(doubledProfile.KingdomAbundance, baseProfile.KingdomAbundance, "kingdom");
        doubledProfile.ShannonDiversity.Should().BeApproximately(baseProfile.ShannonDiversity, 1e-12,
            because: "Shannon diversity is a function of the (unchanged) relative abundances");
        doubledProfile.SimpsonDiversity.Should().BeApproximately(baseProfile.SimpsonDiversity, 1e-12,
            because: "Simpson diversity is a function of the (unchanged) relative abundances");

        // Order independence.
        var shuffled = baseList.AsEnumerable().Reverse().ToList();
        var shuffledProfile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(shuffled);
        AbundancesShouldMatch(shuffledProfile.SpeciesAbundance, baseProfile.SpeciesAbundance, "species");
    }

    #endregion

    #region COMP — per-rank abundances sum to 1

    [Test]
    [Description("COMP: at every rank where each classified read carries a value, the per-taxon counts partition the classified reads, so the relative abundances sum to 1.")]
    public void GenerateTaxonomicProfile_Abundances_SumToOne()
    {
        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(SampleClassifications());

        profile.SpeciesAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12,
            because: "every classified read has a species, so the species abundances partition probability mass 1");
        profile.GenusAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12,
            because: "every classified read has a genus, so the genus abundances sum to 1");
        profile.PhylumAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12,
            because: "every classified read has a phylum, so the phylum abundances sum to 1");
        profile.KingdomAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12,
            because: "every classified read has a kingdom, so the kingdom abundances sum to 1");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: META-ALPHA-001 — alpha diversity (Shannon / Simpson) (Metagenomics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 55.
    //
    // API under test (MetagenomicsAnalyzer.CalculateAlphaDiversity):
    //   ShannonIndex = −Σ pᵢ·ln pᵢ over the (internally normalised) abundances.
    //
    // Relations (derived from the entropy definition, NOT from output):
    //   • COMP (single species ⇒ Shannon = 0): one species has p = 1 and −1·ln 1 = 0 — a
    //          monoculture carries no diversity.
    //   • MON (equalise abundances ⇒ maximal Shannon): for a fixed richness S the entropy is
    //          maximised by the uniform distribution, where Shannon = ln S; any skew lowers it.
    //   • MON (remove a species ⇒ diversity decreases): a uniform community over S species has
    //          Shannon ln S, strictly increasing in S, so dropping a species lowers diversity.
    // ───────────────────────────────────────────────────────────────────────────

    #region Alpha-diversity helpers

    private static Dictionary<string, double> Uniform(int speciesCount) =>
        Enumerable.Range(0, speciesCount).ToDictionary(i => $"sp{i}", _ => 1.0);

    private static double Shannon(IReadOnlyDictionary<string, double> abundances) =>
        MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances).ShannonIndex;

    #endregion

    #region COMP — a single species has Shannon diversity 0

    [Test]
    [Description("COMP: a community of one species has Shannon diversity 0, since the sole abundance has p = 1 and −ln 1 = 0.")]
    public void AlphaDiversity_SingleSpecies_ShannonIsZero()
    {
        foreach (double abundance in new[] { 1.0, 5.0, 100.0 })   // internally normalised, so the raw value is irrelevant
            Shannon(new Dictionary<string, double> { ["only"] = abundance })
                .Should().BeApproximately(0.0, 1e-12,
                    because: "a monoculture has p = 1 at its single species, contributing −1·ln 1 = 0");
    }

    #endregion

    #region MON — the uniform distribution maximises Shannon (= ln S)

    [Test]
    [Description("MON: for a fixed richness S the uniform distribution attains the maximal Shannon ln S; any skewed distribution over the same species scores strictly less.")]
    public void AlphaDiversity_UniformAbundances_MaximiseShannon()
    {
        foreach (int s in new[] { 2, 3, 5, 8 })
        {
            double uniform = Shannon(Uniform(s));
            uniform.Should().BeApproximately(Math.Log(s), 1e-12,
                because: $"the uniform distribution over {s} species has Shannon = ln {s}");

            var skewed = Uniform(s);
            skewed["sp0"] = 10.0;   // make one species dominate
            Shannon(skewed).Should().BeLessThan(uniform - 1e-9,
                because: "any departure from equal abundances lowers the entropy below its uniform maximum");
        }
    }

    #endregion

    #region MON — removing a species decreases diversity

    [Test]
    [Description("MON: a uniform community over S species has Shannon ln S, strictly increasing in S, so removing a species strictly decreases diversity.")]
    public void AlphaDiversity_RemoveSpecies_DecreasesDiversity()
    {
        double previous = double.NegativeInfinity;
        for (int s = 1; s <= 8; s++)
        {
            double shannon = Shannon(Uniform(s));
            shannon.Should().BeGreaterThan(previous,
                because: "a richer uniform community has strictly higher Shannon diversity, so removing a species lowers it");
            previous = shannon;
        }
    }

    #endregion
}
