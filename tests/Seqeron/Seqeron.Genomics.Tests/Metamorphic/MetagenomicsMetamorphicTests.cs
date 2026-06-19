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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: META-BETA-001 — beta diversity (Bray–Curtis / Jaccard) (Metagenomics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 56.
    //
    // API under test (MetagenomicsAnalyzer.CalculateBetaDiversity):
    //   Bray–Curtis = 1 − 2·Σmin(a,b)/Σ(a+b); Jaccard = 1 − shared/(shared+unique1+unique2).
    //
    // Relations (derived from the formulas, NOT from output):
    //   • SYM (symmetry): both indices are symmetric in the two samples (min, sum, shared and
    //          the union count are order-free), so dist(a,b) = dist(b,a).
    //   • COMP (identical samples ⇒ 0): identical samples have Σmin = ½Σtotal and no unique
    //          species, so Bray–Curtis = Jaccard = 0.
    //   • MON (remove a shared species ⇒ higher distance): on the presence/absence Jaccard,
    //          deleting a shared species from one sample turns it into a species unique to the
    //          other — shared drops while the union is unchanged — so the distance strictly
    //          rises. ── Reconciliation: the abundance-weighted Bray–Curtis is NOT monotone
    //          under this transform (removing a low-abundance shared species shrinks both Σmin
    //          and Σtotal and can move it either way), so the monotone relation is asserted on
    //          Jaccard, the index for which it provably holds.
    // ───────────────────────────────────────────────────────────────────────────

    #region Beta-diversity helpers

    private static Dictionary<string, double> Community(params string[] species) =>
        species.ToDictionary(s => s, _ => 10.0);

    private static MetagenomicsAnalyzer.BetaDiversity Beta(
        IReadOnlyDictionary<string, double> a, IReadOnlyDictionary<string, double> b) =>
        MetagenomicsAnalyzer.CalculateBetaDiversity("a", a, "b", b);

    #endregion

    #region SYM — beta diversity is symmetric

    [Test]
    [Description("SYM: Bray–Curtis and Jaccard distances are symmetric in the two samples.")]
    public void BetaDiversity_IsSymmetric()
    {
        var a = Community("sp1", "sp2", "sp3");
        var b = Community("sp2", "sp3", "sp4", "sp5");

        var ab = Beta(a, b);
        var ba = Beta(b, a);

        ba.BrayCurtis.Should().BeApproximately(ab.BrayCurtis, 1e-12,
            because: "Bray–Curtis uses the symmetric Σmin and Σtotal");
        ba.JaccardDistance.Should().BeApproximately(ab.JaccardDistance, 1e-12,
            because: "Jaccard depends on the shared count and union size, both order-free");
    }

    #endregion

    #region COMP — identical samples have distance 0

    [Test]
    [Description("COMP: a sample compared with itself has Bray–Curtis = Jaccard = 0.")]
    public void BetaDiversity_IdenticalSamples_AreZero()
    {
        var s = Community("sp1", "sp2", "sp3", "sp4");
        var self = Beta(s, s);

        self.BrayCurtis.Should().BeApproximately(0.0, 1e-12,
            because: "identical abundances give Σmin = ½Σtotal, so Bray–Curtis = 0");
        self.JaccardDistance.Should().BeApproximately(0.0, 1e-12,
            because: "identical samples share every species and have none unique, so Jaccard = 0");
    }

    #endregion

    #region MON — removing shared species raises the Jaccard distance

    [Test]
    [Description("MON: deleting a shared species from one sample turns it into a species unique to the other, lowering 'shared' while keeping the union fixed, so the Jaccard distance strictly increases.")]
    public void BetaDiversity_RemoveSharedSpecies_IncreasesJaccard()
    {
        var reference = Community("sp1", "sp2", "sp3", "sp4", "sp5");   // 5 species, the fixed union

        // Start with both samples holding all five species (shared = 5), then strip shared
        // species from one sample one at a time. The union stays {sp1..sp5}.
        var sample = new Dictionary<string, double>(reference);
        double previous = double.NegativeInfinity;

        foreach (var toRemove in new[] { (string?)null, "sp1", "sp2", "sp3" })
        {
            if (toRemove != null) sample.Remove(toRemove);

            double jaccard = Beta(sample, reference).JaccardDistance;
            jaccard.Should().BeGreaterThan(previous,
                because: "removing a shared species lowers the shared count while the union is unchanged, so 1 − shared/union strictly increases");
            previous = jaccard;
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: META-BIN-001 — genome binning (Metagenomics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 57.
    //
    // API under test (MetagenomicsAnalyzer.BinContigs):
    //   Deterministic k-means (GC-sorted centroid seeding, no RNG) over GC / coverage /
    //   tetranucleotide features, with k = min(numBins, contigCount); clusters below
    //   minBinSize are dropped. Tested on STRONGLY-SEPARATED synthetic genomes (GC 0 / 0.5 / 1
    //   and coverages two orders of magnitude apart) with numBins = the genome count, so each
    //   genome forms its own bin — the regime in which the checklist relations provably hold.
    //
    // Relations (derived from the per-genome clustering, NOT from output):
    //   • MON (more contigs/genomes ⇒ ≥ bins): each well-separated genome forms one bin, so
    //          adding a genome adds a bin.
    //   • INV (adding a non-overlapping genome ⇒ existing bins unchanged): the contig partition
    //          of the original genomes is identical whether or not a far-separated genome is added.
    // ───────────────────────────────────────────────────────────────────────────

    #region Binning helpers

    private const double BinMinSize = 2000;

    // A "genome" = several near-identical contigs (same GC/TNF unit) with slightly varied coverage.
    private static List<(string ContigId, string Sequence, double Coverage)> Genome(
        string tag, string unit, double baseCoverage, int contigs = 3, int len = 1000)
    {
        string seq = string.Concat(Enumerable.Repeat(unit, len / unit.Length + 1))[..len];
        return Enumerable.Range(0, contigs)
            .Select(i => ($"{tag}_c{i}", seq, baseCoverage + i * 0.1))
            .ToList();
    }

    // GC 0 (AT), GC 1 (GC), GC 0.5 (ATGC); coverages 5 / 500 / 50 — three well-separated genomes.
    private static List<(string, string, double)> GenomeLowGc() => Genome("low", "AT", 5);
    private static List<(string, string, double)> GenomeHighGc() => Genome("high", "GC", 500);
    private static List<(string, string, double)> GenomeMidGc() => Genome("mid", "ATGC", 50);

    private static List<List<string>> PartitionOf(
        IEnumerable<MetagenomicsAnalyzer.GenomeBin> bins, ISet<string> restrictTo) =>
        bins.Select(b => b.ContigIds.Where(restrictTo.Contains).OrderBy(s => s).ToList())
            .Where(ids => ids.Count > 0)
            .ToList();

    #endregion

    #region MON — more genomes give at least as many bins

    [Test]
    [Description("MON: with each strongly-separated genome forming its own bin, adding a genome (and a matching bin budget) yields at least as many bins.")]
    public void BinContigs_MoreGenomes_GiveAtLeastAsManyBins()
    {
        var two = GenomeLowGc().Concat(GenomeHighGc()).ToList();
        var three = two.Concat(GenomeMidGc()).ToList();

        int binsTwo = MetagenomicsAnalyzer.BinContigs(two, numBins: 2, minBinSize: BinMinSize).Count();
        int binsThree = MetagenomicsAnalyzer.BinContigs(three, numBins: 3, minBinSize: BinMinSize).Count();

        binsTwo.Should().Be(2, because: "two well-separated genomes each form one bin");
        binsThree.Should().BeGreaterThanOrEqualTo(binsTwo,
            because: "adding a third well-separated genome cannot reduce the number of recovered bins");
        binsThree.Should().Be(3, because: "three well-separated genomes each form one bin");
    }

    #endregion

    #region INV — adding a far-separated genome leaves the existing bins' partition unchanged

    [Test]
    [Description("INV: the contig partition of the original genomes is identical whether or not a far-separated extra genome is added.")]
    public void BinContigs_AddNonOverlappingGenome_PreservesExistingBins()
    {
        var two = GenomeLowGc().Concat(GenomeHighGc()).ToList();
        var three = two.Concat(GenomeMidGc()).ToList();

        var originalIds = two.Select(c => c.Item1).ToHashSet();

        var basePartition = PartitionOf(
            MetagenomicsAnalyzer.BinContigs(two, numBins: 2, minBinSize: BinMinSize), originalIds);
        var withExtraPartition = PartitionOf(
            MetagenomicsAnalyzer.BinContigs(three, numBins: 3, minBinSize: BinMinSize), originalIds);

        withExtraPartition.Should().BeEquivalentTo(basePartition,
            because: "a far-separated added genome forms its own bin and does not move the original genomes' contigs between bins");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: META-FUNC-001 — functional annotation (Metagenomics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 194.
    //
    // API under test (MetagenomicsAnalyzer.PredictFunctions):
    //   Transfers the best-hit annotation of a signature database to each gene whose protein
    //   contains a database signature.
    //
    // Relations (derived from per-gene best-hit transfer, NOT from output):
    //   • INV  (read order independent): each gene is annotated independently, so reordering the genes
    //          yields the same annotation set.
    //   • SUB  (larger DB ⇒ ≥ assignments): adding database signatures can only let more genes match,
    //          so the set of annotated genes grows.
    // ───────────────────────────────────────────────────────────────────────────

    private static readonly (string, string)[] FuncGenes =
    {
        ("g1", "MKLVAGWTYSDE"), // contains VAG
        ("g2", "ACDEFGHIK"),    // contains DEF
        ("g3", "PQRSTPQRST"),   // contains QRS
    };

    private static readonly Dictionary<string, (string, string, string)> FuncDbSmall = new()
    {
        ["VAG"] = ("FunctionA", "PathA", "K1"),
        ["DEF"] = ("FunctionB", "PathB", "K2"),
    };

    private static readonly Dictionary<string, (string, string, string)> FuncDbLarge = new(FuncDbSmall)
    {
        ["QRS"] = ("FunctionC", "PathC", "K3"),
    };

    #region META-FUNC-001 INV — annotation is independent of gene order

    [Test]
    [Description("INV: each gene is annotated from its own best database hit, so reordering the genes yields the same set of (gene, function) annotations.")]
    public void Functions_ReadOrder_Invariant()
    {
        var forward = MetagenomicsAnalyzer.PredictFunctions(FuncGenes, FuncDbLarge)
            .Select(a => (a.GeneId, a.Function)).ToHashSet();
        var reversed = MetagenomicsAnalyzer.PredictFunctions(FuncGenes.Reverse(), FuncDbLarge)
            .Select(a => (a.GeneId, a.Function)).ToHashSet();

        reversed.Should().BeEquivalentTo(forward, because: "per-gene annotation does not depend on the order of the genes");
    }

    #endregion

    #region META-FUNC-001 SUB — a larger database yields more assignments

    [Test]
    [Description("SUB: adding signatures to the database can only let more genes match, so the set of annotated genes from the larger DB is a superset of the smaller DB's.")]
    public void Functions_LargerDatabase_MoreAssignments()
    {
        var small = MetagenomicsAnalyzer.PredictFunctions(FuncGenes, FuncDbSmall).Select(a => a.GeneId).ToHashSet();
        var large = MetagenomicsAnalyzer.PredictFunctions(FuncGenes, FuncDbLarge).Select(a => a.GeneId).ToHashSet();

        small.IsSubsetOf(large).Should().BeTrue(because: "the larger DB still contains every signature of the smaller one");
        large.Count.Should().BeGreaterThan(small.Count, because: "the added QRS signature annotates the previously-unmatched g3");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: META-PATHWAY-001 — pathway enrichment (Metagenomics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 195.
    //
    // API under test (MetagenomicsAnalyzer.FindPathwayEnrichment):
    //   Hypergeometric over-representation test; a smaller p-value means stronger enrichment.
    //
    // Relations (derived from the hypergeometric upper tail, NOT from output):
    //   • MON  (more pathway genes ⇒ higher enrichment): increasing the query's overlap with a
    //          pathway (query size fixed) lowers P(X ≥ x), i.e. raises enrichment.
    //   • INV  (gene order independent): the query and background are sets, so reordering the query
    //          genes leaves the enrichment unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    private static readonly string[] EnrichPathway = Enumerable.Range(1, 10).Select(i => $"p{i}").ToArray();
    private static readonly Dictionary<string, IReadOnlyCollection<string>> EnrichDb = new() { ["P"] = EnrichPathway };
    private static readonly List<string> EnrichBackground =
        EnrichPathway.Concat(Enumerable.Range(1, 90).Select(i => $"b{i}")).ToList();

    // Query of fixed size 10 with k pathway genes and (10-k) background non-pathway genes.
    private static IEnumerable<string> EnrichQuery(int pathwayGenes) =>
        EnrichPathway.Take(pathwayGenes).Concat(Enumerable.Range(1, 10 - pathwayGenes).Select(i => $"b{i}"));

    #region META-PATHWAY-001 MON — more pathway genes raise enrichment

    [Test]
    [Description("MON: at fixed query size, increasing the overlap with the pathway lowers the hypergeometric p-value (stronger enrichment).")]
    public void Pathway_MoreOverlap_HigherEnrichment()
    {
        double previous = double.MaxValue;
        foreach (int k in new[] { 2, 4, 6, 8 })
        {
            double p = MetagenomicsAnalyzer.FindPathwayEnrichment(EnrichQuery(k), EnrichDb, EnrichBackground)
                .Single(e => e.Pathway == "P").PValue;
            p.Should().BeLessThan(previous, because: $"a query with {k} of the pathway's genes is more enriched (lower p) than one with fewer");
            previous = p;
        }
    }

    #endregion

    #region META-PATHWAY-001 INV — enrichment is independent of gene order

    [Test]
    [Description("INV: query and background are treated as sets, so reordering the query genes yields the same enrichment p-value.")]
    public void Pathway_GeneOrder_Invariant()
    {
        var query = EnrichQuery(5).ToList();
        double forward = MetagenomicsAnalyzer.FindPathwayEnrichment(query, EnrichDb, EnrichBackground).Single(e => e.Pathway == "P").PValue;
        double reversed = MetagenomicsAnalyzer.FindPathwayEnrichment(((IEnumerable<string>)query).Reverse(), EnrichDb, EnrichBackground).Single(e => e.Pathway == "P").PValue;

        reversed.Should().Be(forward, because: "the hypergeometric test depends on the gene sets, not their order");
    }

    #endregion
}
