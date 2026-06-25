using System.Collections.Generic;
using System.Linq;
using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

using TaxonomicClassification = MetagenomicsAnalyzer.TaxonomicClassification;

/// <summary>
/// Algebraic-law tests for the Metagenomics area (taxonomic profiling, beta
/// diversity).
///
/// Algebraic testing pins the normalization identity of a relative-abundance
/// profile (the abundances form a probability distribution summing to one) and
/// the metric/identity behaviour of the pairwise beta-diversity dissimilarities.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 54, 56.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Metagenomics")]
public class MetagenomicsAlgebraicTests
{
    private static TaxonomicClassification Classify(string readId, string species) =>
        new(
            ReadId: readId,
            TaxonId: 1,
            TaxonName: species,
            Rank: "species",
            RtlScore: 1,
            Confidence: 1.0,
            MatchedKmers: 10,
            TotalKmers: 10,
            Kingdom: "Bacteria",
            Phylum: "Phylum_" + species,
            Class: "Class",
            Order: "Order",
            Family: "Family",
            Genus: "Genus_" + species,
            Species: species);

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-PROF-001 — Taxonomic profile (Metagenomics)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 54.
    //
    // Model: a taxonomic profile normalizes per-taxon read counts by the number of
    //        classified reads, producing relative abundances that form a probability
    //        distribution over the observed taxa.
    //   — docs/algorithms/Metagenomics; MetagenomicsAnalyzer.GenerateTaxonomicProfile.
    //
    // Laws (row 54): DIST — Σ species abundances = 1.0.  ID — a single species
    //                accounts for the entire sample (abundance = 1).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>DIST: species (and kingdom) relative abundances sum to 1.0.</summary>
    [FsCheck.NUnit.Property]
    public Property Profile_Distributive_AbundancesSumToOne()
    {
        var speciesList = Gen.Elements("S1", "S2", "S3", "S4", "S5")
            .NonEmptyListOf().ToArbitrary();
        return Prop.ForAll(speciesList, species =>
        {
            var reads = species.Select((s, i) => Classify($"r{i}", s));
            var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
            double speciesSum = profile.SpeciesAbundance.Values.Sum();
            double kingdomSum = profile.KingdomAbundance.Values.Sum();
            return (Math.Abs(speciesSum - 1.0) < 1e-9 && Math.Abs(kingdomSum - 1.0) < 1e-9)
                .Label($"speciesSum={speciesSum}, kingdomSum={kingdomSum}");
        });
    }

    /// <summary>ID: a sample of one species has that species at abundance 1.0.</summary>
    [Test]
    public void Profile_Identity_SingleSpeciesIsOne()
    {
        var reads = Enumerable.Range(0, 20).Select(i => Classify($"r{i}", "OnlySpecies"));
        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(reads);
        profile.SpeciesAbundance.Should().ContainKey("OnlySpecies");
        profile.SpeciesAbundance["OnlySpecies"].Should().BeApproximately(1.0, 1e-12);
        profile.SpeciesAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-12);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-BETA-001 — Beta diversity (Metagenomics)
    // Checklist: docs/checklists/06_ALGEBRAIC_TESTING.md, row 56.
    //
    // Model: pairwise between-sample dissimilarity. Bray–Curtis = 1 − 2Σmin/Σtotal
    //        and Jaccard distance are symmetric and vanish when a sample is compared
    //        to itself.
    //   — docs/algorithms/Metagenomics; MetagenomicsAnalyzer.CalculateBetaDiversity.
    //
    // Laws (row 56): ID — dist(x, x) = 0.  COMM — dist(a, b) = dist(b, a).
    // ═══════════════════════════════════════════════════════════════════════

    private static Arbitrary<Dictionary<string, double>> SampleArbitrary() =>
        (from keys in Gen.Elements("A", "B", "C", "D", "E").ListOf()
         from vals in Gen.Choose(1, 100).Select(x => x / 100.0).ListOf()
         let pairs = keys.Distinct().Zip(vals, (k, v) => (k, v))
         select pairs.ToDictionary(p => p.k, p => p.v))
        .ToArbitrary();

    /// <summary>ID: beta diversity of a sample against itself is 0 (Bray–Curtis, Jaccard).</summary>
    [FsCheck.NUnit.Property]
    public Property Beta_Identity_SelfDistanceIsZero()
    {
        return Prop.ForAll(SampleArbitrary(), sample =>
        {
            if (sample.Count == 0) return true.ToProperty();
            var bd = MetagenomicsAnalyzer.CalculateBetaDiversity("x", sample, "x", sample);
            return (Math.Abs(bd.BrayCurtis) < 1e-12 && Math.Abs(bd.JaccardDistance) < 1e-12)
                .Label($"bray={bd.BrayCurtis}, jaccard={bd.JaccardDistance}");
        });
    }

    /// <summary>COMM: beta diversity is symmetric — dist(a, b) = dist(b, a).</summary>
    [FsCheck.NUnit.Property]
    public Property Beta_Commutative_Symmetric()
    {
        return Prop.ForAll(SampleArbitrary(), SampleArbitrary(), (a, b) =>
        {
            var ab = MetagenomicsAnalyzer.CalculateBetaDiversity("a", a, "b", b);
            var ba = MetagenomicsAnalyzer.CalculateBetaDiversity("b", b, "a", a);
            return (Math.Abs(ab.BrayCurtis - ba.BrayCurtis) < 1e-12
                    && Math.Abs(ab.JaccardDistance - ba.JaccardDistance) < 1e-12)
                .Label($"bray {ab.BrayCurtis}/{ba.BrayCurtis}, jaccard {ab.JaccardDistance}/{ba.JaccardDistance}");
        });
    }
}
