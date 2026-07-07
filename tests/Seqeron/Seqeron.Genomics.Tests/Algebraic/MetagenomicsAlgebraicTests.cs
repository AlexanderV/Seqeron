using FsCheck;
using FsCheck.Fluent;

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

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-CHECKM-001 — CheckM marker-gene bin quality (Metagenomics), row 248.
    //
    // Model: completeness = 100·(1/|M|)·Σ_s |s∩G_M|/|s|; contamination =
    //        100·(1/|M|)·Σ_s Σ_{g∈s}(N_g−1)/|s| (Parks et al. 2015, CheckM Eqs. 1–2).
    //   — MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts; TestSpec META-CHECKM-001.
    //
    // Laws (row 248): ID — a complete single-copy set (every marker present exactly once)
    //                 → 100% completeness / 0% contamination.  IDEMP — deterministic.
    // ═══════════════════════════════════════════════════════════════════════

    private static MetagenomicsAnalyzer.MarkerSet MarkerSet(params string[] ids) => new(ids);

    [Test]
    public void CheckM_Identity_CompleteSingleCopySetIs100And0()
    {
        // Every marker of every set present exactly once: present_s/|s| = 1 ⇒ comp = 100%;
        // C_g = N−1 = 0 for all ⇒ cont = 0%.
        var sets = new[] { MarkerSet("A", "B"), MarkerSet("C", "D", "E"), MarkerSet("F") };
        var counts = new Dictionary<string, int>
        {
            ["A"] = 1, ["B"] = 1, ["C"] = 1, ["D"] = 1, ["E"] = 1, ["F"] = 1,
        };

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);
        q.Completeness.Should().BeApproximately(100.0, 1e-10);
        q.Contamination.Should().BeApproximately(0.0, 1e-10);
    }

    [Test]
    public void CheckM_Idempotent_Deterministic()
    {
        var sets = new[] { MarkerSet("A", "B"), MarkerSet("C") };
        var counts = new Dictionary<string, int> { ["A"] = 1, ["B"] = 2, ["C"] = 1 };

        var q1 = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);
        var q2 = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);
        q2.Completeness.Should().Be(q1.Completeness);
        q2.Contamination.Should().Be(q1.Contamination);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-TETRA-001 — Tetranucleotide z-scores (Metagenomics), row 249.
    //
    // Model: the z-score of a tetramer is its observed-minus-expected count over the
    //        Markov (n−1/n−2/n−3-mer) standard deviation, the TETRA composition vector
    //        of Teeling et al. (2004) used for binning.
    //   — MetagenomicsAnalyzer.CalculateTetranucleotideZScores; TestSpec META-TETRA-001.
    //
    // Laws (row 249): ID — the analytic strand-symmetric construction gives z(ACGT) = √5
    //                 exactly (a formula pin, not a self-comparison).  IDEMP — deterministic.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Tetra_Identity_AnalyticZScoreIsSqrt5()
    {
        // "ACGTACGTGGCC" (folded with its reverse complement inside the estimator) yields, for the
        // tetramer ACGT: observed 4, expected 3.2, variance 0.128 ⇒ z = 0.8/√0.128 = √5 exactly.
        var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores("ACGTACGTGGCC");
        z["ACGT"].Should().BeApproximately(System.Math.Sqrt(5.0), 1e-10);
    }

    [Test]
    public void Tetra_Idempotent_Deterministic()
    {
        const string seq = "ACGTACGTGGCCATGCATGCTTAA";
        var a = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);
        var b = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);
        a.Should().Equal(b);
    }
}
