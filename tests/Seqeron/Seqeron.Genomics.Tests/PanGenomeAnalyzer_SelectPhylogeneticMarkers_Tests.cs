// PANGEN-MARKER-001 — Phylogenetic Marker Selection (single-copy core genes ranked by PIS)
// Evidence: docs/Evidence/PANGEN-MARKER-001-Evidence.md
// TestSpec: tests/TestSpecs/PANGEN-MARKER-001.md
// Source: Ding W, Baumdicker F, Neher RA (2018) Nucleic Acids Research 46(1):e5 (panX);
//         Page AJ et al. (2015) Bioinformatics 31(22):3691 (Roary);
//         Zvelebil M, Baum JO (2008) Understanding Bioinformatics (parsimony-informative site).

using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class PanGenomeAnalyzer_SelectPhylogeneticMarkers_Tests
{
    /// <summary>Builds a genome → (geneId, sequence) dictionary.</summary>
    private static Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string Genome, (string GeneId, string Sequence)[] Genes)[] entries)
    {
        var dict = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        foreach (var (genome, genes) in entries)
            dict[genome] = genes.ToList();
        return dict;
    }

    /// <summary>Constructs a GeneCluster with explicit gene ids and genome membership.</summary>
    private static PanGenomeAnalyzer.GeneCluster Cluster(
        string id, string[] geneIds, string[] genomeIds) =>
        new(id, geneIds, genomeIds, genomeIds.Length, 1.0, geneIds.Length > 0 ? geneIds[0] : "");

    #region CountParsimonyInformativeSites (canonical)

    // M1 — Worked alignment (Evidence dataset; Zvelebil 2008). s1=AAAAA, s2=AAACA,
    // s3=AACCG, s4=ACCTG. Only columns 3 (A,A,C,C) and 5 (A,A,G,G) are PI -> 2.
    [Test]
    public void CountParsimonyInformativeSites_WorkedAlignment_ReturnsTwo()
    {
        var aln = new[] { "AAAAA", "AAACA", "AACCG", "ACCTG" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln);

        Assert.That(pis, Is.EqualTo(2),
            "Cols 3 (A,A,C,C) and 5 (A,A,G,G) have two states each in >=2 sequences; "
            + "cols 1 (mono), 2 (singleton C), 4 (four singletons) are not PI (Zvelebil 2008).");
    }

    // M2 — Monomorphic column: all rows identical -> 0 (no signal).
    [Test]
    public void CountParsimonyInformativeSites_MonomorphicColumn_ReturnsZero()
    {
        var aln = new[] { "A", "A", "A", "A" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln);

        Assert.That(pis, Is.EqualTo(0),
            "A single-state column carries no phylogenetic signal (Zvelebil 2008).");
    }

    // M3 — Singleton column: one row differs -> not PI (variant in only one sequence).
    [Test]
    public void CountParsimonyInformativeSites_SingletonColumn_ReturnsZero()
    {
        var aln = new[] { "A", "A", "A", "C" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln);

        Assert.That(pis, Is.EqualTo(0),
            "C occurs in only one sequence -> singleton, not parsimony-informative (Zvelebil 2008).");
    }

    // M4 — Minimal informative column: two states each in two rows -> 1.
    [Test]
    public void CountParsimonyInformativeSites_MinimalInformativeColumn_ReturnsOne()
    {
        var aln = new[] { "A", "A", "C", "C" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln);

        Assert.That(pis, Is.EqualTo(1),
            "Two states (A,C) each in >=2 sequences is the minimal PI pattern (Zvelebil 2008).");
    }

    // M5 — Four-singleton column: four distinct states, none in >=2 rows -> not PI.
    [Test]
    public void CountParsimonyInformativeSites_FourSingletonsColumn_ReturnsZero()
    {
        var aln = new[] { "A", "C", "G", "T" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln);

        Assert.That(pis, Is.EqualTo(0),
            "Four singletons: no state occurs in >=2 sequences, so not PI (Zvelebil 2008).");
    }

    // S2 — Single sequence (fewer than 2 rows) -> 0 (no state can be in >=2 rows).
    [Test]
    public void CountParsimonyInformativeSites_SingleSequence_ReturnsZero()
    {
        var aln = new[] { "ACGT" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln);

        Assert.That(pis, Is.EqualTo(0),
            "With one row no state can occur in >=2 sequences -> 0 PIS.");
    }

    // S3 — Empty / null inputs -> 0, no exception (INV-1 lower bound).
    [Test]
    public void CountParsimonyInformativeSites_EmptyOrNull_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PanGenomeAnalyzer.CountParsimonyInformativeSites(null!), Is.EqualTo(0),
                "Null alignment -> 0, no exception.");
            Assert.That(PanGenomeAnalyzer.CountParsimonyInformativeSites(new string[0]), Is.EqualTo(0),
                "Empty alignment -> 0.");
            Assert.That(PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "", "" }), Is.EqualTo(0),
                "Zero-length rows -> 0 columns -> 0 PIS.");
        });
    }

    // S1 (PIS side) — Unequal-length rows: no common alignment -> 0.
    [Test]
    public void CountParsimonyInformativeSites_UnequalLengths_ReturnsZero()
    {
        var aln = new[] { "AACC", "AAC" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln);

        Assert.That(pis, Is.EqualTo(0),
            "Rows of differing length share no common alignment -> 0 PIS (Assumption 1).");
    }

    // C1 — Row-order invariance: permuting rows of the M1 alignment keeps PIS = 2 (INV-7).
    [Test]
    public void CountParsimonyInformativeSites_RowOrderPermuted_Unchanged()
    {
        var ordered = new[] { "AAAAA", "AAACA", "AACCG", "ACCTG" };
        var permuted = new[] { "ACCTG", "AACCG", "AAAAA", "AAACA" };

        int a = PanGenomeAnalyzer.CountParsimonyInformativeSites(ordered);
        int b = PanGenomeAnalyzer.CountParsimonyInformativeSites(permuted);

        Assert.That(b, Is.EqualTo(a).And.EqualTo(2),
            "PIS is a per-column property independent of row order (INV-7).");
    }

    // C2 — State-relabel invariance: swapping A<->C in M1 alignment keeps PIS = 2 (INV-7).
    [Test]
    public void CountParsimonyInformativeSites_StateRelabeled_Unchanged()
    {
        // Original M1 with A and C swapped (G,T untouched): a bijective relabeling.
        var relabeled = new[] { "CCCCC", "CCCAC", "CCAAG", "CAATG" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(relabeled);

        Assert.That(pis, Is.EqualTo(2),
            "Bijective state relabeling preserves column partitions, so PIS is unchanged (INV-7).");
    }

    #endregion

    #region SelectPhylogeneticMarkers (canonical)

    // M9 / M6 / M7 / M8 — Mixed candidate set: only the single-copy core informative
    // cluster is selected; paralog, non-core, and conserved clusters are excluded.
    [Test]
    public void SelectPhylogeneticMarkers_MixedCandidates_SelectsOnlyInformativeSingleCopyCore()
    {
        // Four single-copy genomes; an informative marker needs 2 states each in >=2 rows.
        var genomes = Genomes(
            ("g1", new[] { ("inf1", "AACC"), ("par1a", "AAAA"), ("par1b", "AAAA"), ("con1", "AAAA"), ("nc1", "AACC") }),
            ("g2", new[] { ("inf2", "AACC"), ("par2", "AAAA"), ("con2", "AAAA"), ("nc2", "AACC") }),
            ("g3", new[] { ("inf3", "GGCC"), ("par3", "AAAA"), ("con3", "AAAA") }),
            ("g4", new[] { ("inf4", "GGCC"), ("par4", "AAAA"), ("con4", "AAAA") }));

        var candidates = new[]
        {
            // Single-copy core, PIS: cols A,A,G,G (col1) and A,A,G,G (col2) -> 2 PI columns.
            Cluster("informative", new[] { "inf1", "inf2", "inf3", "inf4" },
                new[] { "g1", "g2", "g3", "g4" }),
            // Paralog: g1 contributes two genes (par1a, par1b) -> 5 gene ids for 4 genomes.
            Cluster("paralog", new[] { "par1a", "par1b", "par2", "par3", "par4" },
                new[] { "g1", "g2", "g3", "g4" }),
            // Not core: only 2 of 4 genomes (g1,g2).
            Cluster("noncore", new[] { "nc1", "nc2" }, new[] { "g1", "g2" }),
            // Conserved single-copy core: all members identical -> 0 PIS.
            Cluster("conserved", new[] { "con1", "con2", "con3", "con4" },
                new[] { "g1", "g2", "g3", "g4" }),
        };

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, candidates, totalGenomes: 4)
            .Select(m => m.ClusterId)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(markers, Is.EqualTo(new[] { "informative" }),
                "Only the single-copy core cluster with >=1 PIS is selected (panX/Roary; INV-04/05).");
            Assert.That(markers, Does.Not.Contain("paralog"),
                "A genome contributing >=2 genes is not single-copy (Roary paralog filtering).");
            Assert.That(markers, Does.Not.Contain("noncore"),
                "A cluster absent from a genome is not core (panX 'all strains').");
            Assert.That(markers, Does.Not.Contain("conserved"),
                "A fully conserved cluster has 0 variable positions -> excluded (panX).");
        });
    }

    // M10 / M11 — Ranking by descending PIS and maxMarkers cap. Cluster ids are chosen so
    // that an id-only ordering would place the lower-PIS marker first; descending-PIS must
    // override that and put the higher-PIS marker first.
    //   aHi members: AC,AC,GG,GG -> col1 A,A,G,G (PI); col2 C,C,G,G (PI)        => PIS 2
    //   zLo members: AC,AC,GC,GT -> col1 A,A,G,G (PI); col2 C,C,C,T (singleton) => PIS 1
    [Test]
    public void SelectPhylogeneticMarkers_OrdersByDescendingPisAndCaps()
    {
        var genomes = Genomes(
            ("g1", new[] { ("hi1", "AC"), ("lo1", "AC") }),
            ("g2", new[] { ("hi2", "AC"), ("lo2", "AC") }),
            ("g3", new[] { ("hi3", "GG"), ("lo3", "GC") }),
            ("g4", new[] { ("hi4", "GG"), ("lo4", "GT") }));

        var candidates = new[]
        {
            // 'aHi' < 'zLo' ordinally; descending-PIS must still keep 'aHi' first by its higher PIS.
            Cluster("aHi", new[] { "hi1", "hi2", "hi3", "hi4" }, new[] { "g1", "g2", "g3", "g4" }),
            Cluster("zLo", new[] { "lo1", "lo2", "lo3", "lo4" }, new[] { "g1", "g2", "g3", "g4" }),
        };

        int pisHi = PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "AC", "AC", "GG", "GG" });
        int pisLo = PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "AC", "AC", "GC", "GT" });

        var all = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, candidates, totalGenomes: 4)
            .Select(m => m.ClusterId).ToList();
        var capped = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, candidates, totalGenomes: 4, maxMarkers: 1)
            .Select(m => m.ClusterId).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(pisHi, Is.EqualTo(2), "aHi alignment has 2 PI columns.");
            Assert.That(pisLo, Is.EqualTo(1), "zLo alignment has 1 PI column (col2 is a singleton).");
            Assert.That(all, Is.EqualTo(new[] { "aHi", "zLo" }),
                "Markers are ordered by descending PIS (2 before 1), not by cluster id (INV-06).");
            Assert.That(capped, Is.EqualTo(new[] { "aHi" }),
                "maxMarkers=1 returns exactly the most-informative marker (INV-06).");
        });
    }

    // INV-06 tie-break — equal PIS resolves deterministically by ordinal cluster id.
    [Test]
    public void SelectPhylogeneticMarkers_EqualPis_TieBrokenByOrdinalClusterId()
    {
        // Both markers have identical informative alignments (PIS 2); ids 'mA' < 'mB'.
        var genomes = Genomes(
            ("g1", new[] { ("a1", "AC"), ("b1", "AC") }),
            ("g2", new[] { ("a2", "AC"), ("b2", "AC") }),
            ("g3", new[] { ("a3", "GG"), ("b3", "GG") }),
            ("g4", new[] { ("a4", "GG"), ("b4", "GG") }));

        var candidates = new[]
        {
            Cluster("mB", new[] { "b1", "b2", "b3", "b4" }, new[] { "g1", "g2", "g3", "g4" }),
            Cluster("mA", new[] { "a1", "a2", "a3", "a4" }, new[] { "g1", "g2", "g3", "g4" }),
        };

        var ordered = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, candidates, totalGenomes: 4)
            .Select(m => m.ClusterId).ToList();

        Assert.That(ordered, Is.EqualTo(new[] { "mA", "mB" }),
            "Equal PIS ties are broken by ordinal cluster id ('mA' < 'mB') for determinism (INV-06).");
    }

    // S1 — Unequal-length members in a single-copy core cluster -> PIS 0 -> not selected.
    [Test]
    public void SelectPhylogeneticMarkers_UnequalLengthMembers_NotSelected()
    {
        var genomes = Genomes(
            ("g1", new[] { ("u1", "AACC") }),
            ("g2", new[] { ("u2", "AAC") }),   // shorter
            ("g3", new[] { ("u3", "GGCC") }),
            ("g4", new[] { ("u4", "GGCC") }));

        var candidates = new[]
        {
            Cluster("ragged", new[] { "u1", "u2", "u3", "u4" }, new[] { "g1", "g2", "g3", "g4" }),
        };

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, candidates, totalGenomes: 4).ToList();

        Assert.That(markers, Is.Empty,
            "Members of differing length have no common alignment -> PIS 0 -> not selected (Assumption 1).");
    }

    // M12 — Null / empty inputs -> empty result, no exception.
    [Test]
    public void SelectPhylogeneticMarkers_NullOrEmpty_ReturnsEmpty()
    {
        var genomes = Genomes(("g1", new[] { ("x", "AACC") }));
        var clusters = new[] { Cluster("c", new[] { "x" }, new[] { "g1" }) };

        Assert.Multiple(() =>
        {
            Assert.That(PanGenomeAnalyzer.SelectPhylogeneticMarkers(null!, clusters, 1), Is.Empty,
                "Null genomes -> empty.");
            Assert.That(PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, null!, 1), Is.Empty,
                "Null clusters -> empty.");
            Assert.That(PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, clusters, 0), Is.Empty,
                "totalGenomes <= 0 -> empty.");
            Assert.That(PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, clusters, 1, maxMarkers: 0), Is.Empty,
                "maxMarkers <= 0 -> empty.");
            Assert.That(PanGenomeAnalyzer.SelectPhylogeneticMarkers(
                    new Dictionary<string, IReadOnlyList<(string, string)>>(),
                    System.Array.Empty<PanGenomeAnalyzer.GeneCluster>(), 1),
                Is.Empty, "Empty inputs -> empty.");
        });
    }

    #endregion
}
