// PANGEN-CORE-001 — Core/Accessory Genome (pan-genome construction)
// Evidence: docs/Evidence/PANGEN-CORE-001-Evidence.md
// TestSpec: tests/TestSpecs/PANGEN-CORE-001.md
// Source: Tettelin H et al. (2005) PNAS 102:13950; Tettelin H et al. (2008) Curr Opin Microbiol 11:472;
//         Kislyuk AO et al. (2011) BMC Genomics 12:32; Page AJ et al. (2015) Bioinformatics 31:3691.

using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests.Unit.Metagenomics;

[TestFixture]
public class PanGenomeAnalyzer_ConstructPanGenome_Tests
{
    private static Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string Genome, (string GeneId, string Sequence)[] Genes)[] entries)
    {
        var dict = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        foreach (var (genome, genes) in entries)
            dict[genome] = genes.ToList();
        return dict;
    }

    // Distinct 30-bp sequences (>= k=7) so identical strings cluster together and
    // different strings form separate clusters under the k-mer Jaccard clusterer.
    private const string SeqCore = "ATGCGATCGATCGATCGATCGATCGATCGA"; // shared "core" gene
    private const string SeqShared2 = "TTACGGCATTACGGCATTACGGCATTACGG"; // shared by two genomes
    private const string SeqU1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    private const string SeqU2 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    private const string SeqU3 = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

    #region ConstructPanGenome — Partition (Tettelin 2005/2008; Page 2015)

    // M1 — Core/accessory/unique partition by cluster occupancy. c1 in all 3 -> core;
    // c2 in 2 -> accessory; three singletons -> unique. coreFraction 1.0 => threshold 3.
    [Test]
    public void ConstructPanGenome_ThreeGenomes_PartitionsCoreAccessoryUnique()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqCore), ("b", SeqShared2), ("c", SeqU1) }),
            ("g2", new[] { ("a2", SeqCore), ("b2", SeqShared2), ("d", SeqU2) }),
            ("g3", new[] { ("a3", SeqCore), ("e", SeqU3) }));

        var result = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0);

        Assert.Multiple(() =>
        {
            // c (SeqCore) present in all 3 genomes -> exactly one core cluster.
            Assert.That(result.Statistics.CoreGeneCount, Is.EqualTo(1),
                "SeqCore occurs in all 3 genomes; with coreFraction 1.0 it is the sole core cluster (Tettelin 2008).");
            // SeqShared2 present in 2 of 3 -> accessory/shell.
            Assert.That(result.Statistics.AccessoryGeneCount, Is.EqualTo(1),
                "SeqShared2 occurs in 2 of 3 genomes -> accessory (present in some but not all).");
            // SeqU1, SeqU2, SeqU3 each in exactly one genome -> 3 unique clusters.
            Assert.That(result.Statistics.UniqueGeneCount, Is.EqualTo(3),
                "Three strain-specific sequences each occur in exactly one genome -> unique (Tettelin 2005).");
            Assert.That(result.Statistics.TotalGenes, Is.EqualTo(5),
                "Total gene clusters = 1 core + 1 accessory + 3 unique = 5.");
        });
    }

    // M2 — INV-02: the three partitions sum to TotalGenes (every cluster classified once).
    [Test]
    public void ConstructPanGenome_ThreeGenomes_PartitionCountsSumToTotal()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqCore), ("b", SeqShared2), ("c", SeqU1) }),
            ("g2", new[] { ("a2", SeqCore), ("b2", SeqShared2), ("d", SeqU2) }),
            ("g3", new[] { ("a3", SeqCore), ("e", SeqU3) }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        Assert.That(s.CoreGeneCount + s.AccessoryGeneCount + s.UniqueGeneCount, Is.EqualTo(s.TotalGenes),
            "INV-02: core + accessory + unique must equal the total number of gene clusters.");
    }

    // M10 — INV-06: CoreFraction = CoreGeneCount / TotalGenes.
    [Test]
    public void ConstructPanGenome_ThreeGenomes_CoreFractionEqualsCoreOverTotal()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqCore), ("b", SeqShared2), ("c", SeqU1) }),
            ("g2", new[] { ("a2", SeqCore), ("b2", SeqShared2), ("d", SeqU2) }),
            ("g3", new[] { ("a3", SeqCore), ("e", SeqU3) }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        // 1 core / 5 total = 0.2 exactly.
        Assert.That(s.CoreFraction, Is.EqualTo(1.0 / 5.0).Within(1e-10),
            "INV-06: CoreFraction = CoreGeneCount (1) / TotalGenes (5) = 0.2.");
    }

    // S4 — Core threshold is a FRACTIONAL test: Roary core = "a gene being in at least 99%
    // of samples" (Page et al., 2015). With coreFraction 0.99 and N=3, a cluster is core
    // iff occupancy/3 >= 0.99, i.e. occupancy = 3 only. A 2-of-3 (66.7%) cluster is NOT
    // core (it is accessory/shell) — this guards against the unsourced floor(0.99*3)=2
    // convention that would wrongly include it.
    [Test]
    public void ConstructPanGenome_CoreFraction099_OnlyFullyConservedClusterIsCore()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqCore), ("b", SeqShared2) }),
            ("g2", new[] { ("a2", SeqCore), ("b2", SeqShared2) }),
            ("g3", new[] { ("a3", SeqCore), ("c", SeqU1) }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99).Statistics;

        Assert.Multiple(() =>
        {
            // Only SeqCore (3/3 = 100% >= 99%) is core; SeqShared2 (2/3 = 66.7% < 99%) is NOT.
            Assert.That(s.CoreGeneCount, Is.EqualTo(1),
                "Roary core = present in >= 99% of samples (Page 2015); only the 3/3 cluster qualifies, not 2/3 (66.7%).");
            Assert.That(s.AccessoryGeneCount, Is.EqualTo(1),
                "SeqShared2 in 2 of 3 genomes (66.7% < 99%) is accessory/shell, not core (Page 2015).");
            Assert.That(s.UniqueGeneCount, Is.EqualTo(1),
                "SeqU1 occurs in one genome -> unique.");
        });
    }

    // S4b — Fractional core at the exact 99% boundary with N=100: a 99/100 cluster (99%)
    // IS core; a 98/100 cluster (98% < 99%) is NOT. Guards the float round-off of 0.99*100
    // (= 98.999...) and confirms the boundary matches Roary's "at least 99%" (Page 2015).
    [Test]
    public void GetCoreGeneClusters_CoreFraction099_NinetyNineOfHundredIsCoreNotNinetyEight()
    {
        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("c99", new[] { "g99" }, GenomeNames(99), 99, 1.0, "ATGC"),
            new("c98", new[] { "g98" }, GenomeNames(98), 98, 1.0, "GCTA"),
        };

        var core = PanGenomeAnalyzer.GetCoreGeneClusters(clusters, totalGenomes: 100, threshold: 0.99).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(core, Has.Count.EqualTo(1),
                "At threshold 0.99 over 100 genomes, only occupancy >= 99 (>= 99%) is core (Page 2015).");
            Assert.That(core[0].ClusterId, Is.EqualTo("c99"),
                "99/100 = 99% >= 99% is core; 98/100 = 98% < 99% is not.");
        });
    }

    private static string[] GenomeNames(int count) =>
        Enumerable.Range(1, count).Select(i => $"genome{i}").ToArray();

    #endregion

    #region ConstructPanGenome — Genome Fluidity (Kislyuk 2011)

    // M3 — Exact genome fluidity on the Evidence hand-derived example:
    // A={c1,c2,c3}, B={c1,c2,c4}, C={c1,c5,c6}. phi = (1/3)(2/6+4/6+4/6) = 10/18.
    [Test]
    public void ConstructPanGenome_HandDerivedExample_GenomeFluidityMatchesKislyukFormula()
    {
        var genomes = Genomes(
            ("A", new[] { ("a_core", SeqCore), ("a_sh", SeqShared2), ("a_u", SeqU1) }),
            ("B", new[] { ("b_core", SeqCore), ("b_sh", SeqShared2), ("b_u", SeqU2) }),
            ("C", new[] { ("c_core", SeqCore), ("c_u1", SeqU3),
                          ("c_u2", "TATATATATATATATATATATATATATATA") }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        // phi = 2/(N(N-1)) * sum_{k<l} (U_k+U_l)/(M_k+M_l) = 10/18 (Kislyuk 2011, eq. for phi).
        Assert.That(s.GenomeFluidity, Is.EqualTo(10.0 / 18.0).Within(1e-10),
            "Genome fluidity must equal the Kislyuk (2011) closed form: (1/3)(2/6 + 4/6 + 4/6) = 10/18.");
    }

    // M4 — Identical gene content across all genomes => fluidity 0 (Kislyuk 2011).
    [Test]
    public void ConstructPanGenome_IdenticalGeneContent_GenomeFluidityIsZero()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqCore) }),
            ("g2", new[] { ("a2", SeqCore) }),
            ("g3", new[] { ("a3", SeqCore) }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        Assert.That(s.GenomeFluidity, Is.EqualTo(0.0).Within(1e-10),
            "All genomes share the same single gene cluster (no unique families) -> fluidity 0 (Kislyuk 2011).");
    }

    // M5 — Pairwise-disjoint gene content => fluidity 1 (Kislyuk 2011).
    [Test]
    public void ConstructPanGenome_DisjointGeneContent_GenomeFluidityIsOne()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqU1) }),
            ("g2", new[] { ("b", SeqU2) }),
            ("g3", new[] { ("c", SeqU3) }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        Assert.That(s.GenomeFluidity, Is.EqualTo(1.0).Within(1e-10),
            "Every genome has a distinct cluster only (all families unique per pair) -> fluidity 1 (Kislyuk 2011).");
    }

    // M6 — INV-04: fluidity bounded in [0,1] on the hand-derived example.
    [Test]
    public void ConstructPanGenome_GenomeFluidity_IsWithinUnitInterval()
    {
        var genomes = Genomes(
            ("A", new[] { ("a_core", SeqCore), ("a_sh", SeqShared2), ("a_u", SeqU1) }),
            ("B", new[] { ("b_core", SeqCore), ("b_sh", SeqShared2), ("b_u", SeqU2) }),
            ("C", new[] { ("c_core", SeqCore), ("c_u1", SeqU3) }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        Assert.That(s.GenomeFluidity, Is.InRange(0.0, 1.0),
            "INV-04: genome fluidity must lie in [0,1] (Kislyuk 2011).");
    }

    #endregion

    #region ConstructPanGenome — Open vs Closed (Tettelin 2008; micropan)

    // M7 — Open: each added genome contributes one brand-new unique cluster, so the
    // new-gene curve stays flat (decay exponent alpha ~ 0 < 1) => Open (Tettelin 2008).
    [Test]
    public void ConstructPanGenome_NewGenesPerGenome_ClassifiedOpen()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqCore), ("u1", SeqU1) }),
            ("g2", new[] { ("a2", SeqCore), ("u2", SeqU2) }),
            ("g3", new[] { ("a3", SeqCore), ("u3", SeqU3) }),
            ("g4", new[] { ("a4", SeqCore), ("u4", "TATATATATATATATATATATATATATATA") }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        Assert.That(s.Type, Is.EqualTo(PanGenomeAnalyzer.PanGenomeType.Open),
            "Constant per-genome novelty -> Heaps decay exponent alpha < 1 -> Open (Tettelin 2008; micropan).");
    }

    // M8 — Closed: novelty decays steeply (4, 2, 1 new clusters at k=2,3,4 -> alpha ~ 2 > 1)
    // => Closed (Tettelin 2008; micropan).
    [Test]
    public void ConstructPanGenome_DecayingNovelty_ClassifiedClosed()
    {
        // g1 establishes the shared baseline; each later genome shares it and adds a
        // decaying number of new clusters: g2:+4, g3:+2, g4:+1.
        var baseGenes = new[]
        {
            ("c1", SeqCore), ("c2", SeqShared2), ("c3", SeqU1), ("c4", SeqU2), ("c5", SeqU3),
        };
        var genomes = Genomes(
            ("g1", baseGenes),
            ("g2", baseGenes.Concat(new[]
            {
                ("n1", "ACACACACACACACACACACACACACACAC"),
                ("n2", "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAG"),
                ("n3", "ATATATATATATATATATATATATATATAT"),
                ("n4", "CTCTCTCTCTCTCTCTCTCTCTCTCTCTCT"),
            }).ToArray()),
            ("g3", baseGenes.Concat(new[]
            {
                ("m1", "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"),
                ("m2", "GTGTGTGTGTGTGTGTGTGTGTGTGTGTGT"),
            }).ToArray()),
            ("g4", baseGenes.Concat(new[]
            {
                ("p1", "TCTCTCTCTCTCTCTCTCTCTCTCTCTCTC"),
            }).ToArray()));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        Assert.That(s.Type, Is.EqualTo(PanGenomeAnalyzer.PanGenomeType.Closed),
            "Steeply decaying novelty (4,2,1) -> Heaps decay exponent alpha > 1 -> Closed (Tettelin 2008; micropan).");
    }

    #endregion

    #region Core-gene identification — GetCoreGeneClusters (Page 2015)

    // M9 — Core-gene identification: threshold 1.0 over occupancy {3,2,1}/3 returns only
    // the occupancy-3 cluster (IdentifyCoreGenes referent).
    [Test]
    public void GetCoreGeneClusters_ThresholdOne_ReturnsOnlyFullyConservedCluster()
    {
        var clusters = new List<PanGenomeAnalyzer.GeneCluster>
        {
            new("c1", new[] { "g1" }, new[] { "genome1", "genome2", "genome3" }, 3, 0.95, "ATGC"),
            new("c2", new[] { "g2" }, new[] { "genome1", "genome2" }, 2, 0.95, "GCTA"),
            new("c3", new[] { "g3" }, new[] { "genome1" }, 1, 1.0, "TTTT"),
        };

        var core = PanGenomeAnalyzer.GetCoreGeneClusters(clusters, totalGenomes: 3, threshold: 1.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(core, Has.Count.EqualTo(1),
                "threshold 1.0 over 3 genomes requires occupancy >= 3; only c1 qualifies (Page 2015).");
            Assert.That(core[0].ClusterId, Is.EqualTo("c1"),
                "The fully conserved cluster (occupancy 3) is the sole core gene.");
        });
    }

    // S5 — Empty cluster list => no core genes.
    [Test]
    public void GetCoreGeneClusters_EmptyClusters_ReturnsEmpty()
    {
        var core = PanGenomeAnalyzer.GetCoreGeneClusters(
            new List<PanGenomeAnalyzer.GeneCluster>(), totalGenomes: 5, threshold: 0.99).ToList();

        Assert.That(core, Is.Empty, "No clusters -> no core genes.");
    }

    #endregion

    #region Edge cases — empty / null / single genome

    // S1 — Empty input -> empty result with zeroed statistics.
    [Test]
    public void ConstructPanGenome_EmptyInput_ReturnsEmptyResult()
    {
        var result = PanGenomeAnalyzer.ConstructPanGenome(
            new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>());

        Assert.Multiple(() =>
        {
            Assert.That(result.CoreGenes, Is.Empty, "No genomes -> no core genes.");
            Assert.That(result.AccessoryGenes, Is.Empty, "No genomes -> no accessory genes.");
            Assert.That(result.UniqueGenes, Is.Empty, "No genomes -> no unique genes.");
            Assert.That(result.Statistics.TotalGenomes, Is.EqualTo(0), "TotalGenomes = 0 for empty input.");
            Assert.That(result.Statistics.TotalGenes, Is.EqualTo(0), "TotalGenes = 0 for empty input.");
            Assert.That(result.Statistics.GenomeFluidity, Is.EqualTo(0.0).Within(1e-10),
                "No genome pairs -> fluidity 0.");
            Assert.That(result.Statistics.Type, Is.EqualTo(PanGenomeAnalyzer.PanGenomeType.Closed),
                "Fewer than 3 genomes -> openness undeterminable -> Closed default.");
        });
    }

    // S2 — Null input -> empty result (no throw), matching the documented contract.
    [Test]
    public void ConstructPanGenome_NullInput_ReturnsEmptyResult()
    {
        var result = PanGenomeAnalyzer.ConstructPanGenome(null!);

        Assert.Multiple(() =>
        {
            Assert.That(result.Statistics.TotalGenomes, Is.EqualTo(0), "Null input -> 0 genomes.");
            Assert.That(result.CoreGenes, Is.Empty, "Null input -> empty partitions.");
        });
    }

    // S3 — Single genome: occupancy of every cluster = N = 1, so all clusters are core
    // (intersection over 1 genome), no accessory, and no pairs -> fluidity 0.
    [Test]
    public void ConstructPanGenome_SingleGenome_AllClustersCoreNoFluidity()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqU1), ("b", SeqU2) }));

        var s = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99).Statistics;

        Assert.Multiple(() =>
        {
            Assert.That(s.TotalGenomes, Is.EqualTo(1), "One genome provided.");
            Assert.That(s.TotalGenes, Is.EqualTo(2), "Two distinct sequences -> two clusters.");
            Assert.That(s.CoreGeneCount, Is.EqualTo(2),
                "With N=1 every cluster is present in all (1) genomes -> all core (Tettelin intersection).");
            Assert.That(s.AccessoryGeneCount, Is.EqualTo(0), "No accessory genes possible with one genome.");
            Assert.That(s.GenomeFluidity, Is.EqualTo(0.0).Within(1e-10),
                "Single genome -> no pairs -> fluidity 0 (Kislyuk 2011).");
            Assert.That(s.Type, Is.EqualTo(PanGenomeAnalyzer.PanGenomeType.Closed),
                "N<3 -> openness undeterminable -> Closed default.");
        });
    }

    #endregion

    #region Determinism and property-based (O(g^2 . s) algorithm)

    // C1 — Determinism: identical input yields identical statistics across runs.
    [Test]
    public void ConstructPanGenome_CalledTwice_ProducesIdenticalStatistics()
    {
        var genomes = Genomes(
            ("g1", new[] { ("a", SeqCore), ("b", SeqShared2), ("c", SeqU1) }),
            ("g2", new[] { ("a2", SeqCore), ("b2", SeqShared2), ("d", SeqU2) }),
            ("g3", new[] { ("a3", SeqCore), ("e", SeqU3) }));

        var s1 = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;
        var s2 = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0).Statistics;

        Assert.Multiple(() =>
        {
            Assert.That(s2.CoreGeneCount, Is.EqualTo(s1.CoreGeneCount), "Core count is deterministic.");
            Assert.That(s2.AccessoryGeneCount, Is.EqualTo(s1.AccessoryGeneCount), "Accessory count is deterministic.");
            Assert.That(s2.UniqueGeneCount, Is.EqualTo(s1.UniqueGeneCount), "Unique count is deterministic.");
            Assert.That(s2.GenomeFluidity, Is.EqualTo(s1.GenomeFluidity).Within(1e-12), "Fluidity is deterministic.");
            Assert.That(s2.Type, Is.EqualTo(s1.Type), "Open/closed classification is deterministic.");
        });
    }

    // C2 — Property (INV-04): genome fluidity stays in [0,1] across structured inputs of
    // varying overlap (O(g^2) pairwise computation).
    [Test]
    public void ConstructPanGenome_VaryingOverlap_FluidityAlwaysInUnitInterval()
    {
        string[] pool = { SeqCore, SeqShared2, SeqU1, SeqU2, SeqU3, "TATATATATATATATATATATATATATATA" };

        Assert.Multiple(() =>
        {
            for (int sharedCount = 0; sharedCount <= 3; sharedCount++)
            {
                var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
                for (int g = 0; g < 3; g++)
                {
                    var genes = new List<(string, string)>();
                    for (int sc = 0; sc < sharedCount; sc++)
                        genes.Add(($"shared_{g}_{sc}", pool[sc]));       // shared families
                    genes.Add(($"uniq_{g}", pool[3 + g]));               // one unique family each
                    genomes[$"g{g}"] = genes;
                }

                double fluidity = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0)
                    .Statistics.GenomeFluidity;

                Assert.That(fluidity, Is.InRange(0.0, 1.0),
                    $"INV-04: fluidity must stay in [0,1] (sharedCount={sharedCount}).");
            }
        });
    }

    #endregion
}
