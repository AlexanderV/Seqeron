// COMPGEN-COMPARE-001 — Comprehensive Genome Comparison (core/dispensable partition + syntenic-gene fraction)
// Evidence: docs/Evidence/COMPGEN-COMPARE-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-COMPARE-001.md
// Source: Tettelin H et al. (2005). PNAS 102(39):13950–13955. https://doi.org/10.1073/pnas.0506758102
//         Moreno-Hagelsieb & Latimer (2008). Bioinformatics 24(3):319–324. (RBH = shared genes)
//         Synteny overview — "fraction of syntenic genes" metric (OverallSynteny).

using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using G = Seqeron.Genomics.Analysis.ComparativeGenomics.Gene;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class ComparativeGenomics_CompareGenomes_Tests
{
    #region Test Data

    // Five DISTINCT >=60-nt sequences, each shared identically by both genomes.
    // Distinct content guarantees an unambiguous reciprocal-best-hit matching
    // (a_i <-> c_i), so the conserved set is exactly these pairs (Source: RBH).
    private static readonly string[] Shared =
    {
        "ATGGCAAAGCTTGATCCGTACGGGTTAACCGGATCAGGTTCAAAGCTTGATCCGTACGGG",
        "TTACCGGATCAGGTTCATGGCAAAGCTTGATCCGTACGGGAATTACCGGATCAGGTTCAT",
        "GGCCAATTGGCCAATTACGTACGTGGCCAATTGGCCAATTACGTACGTGGCCAATTGGCC",
        "CTGACTGACAAATTTGGGCCCCTGACTGACAAATTTGGGCCCCTGACTGACAAATTTGGG",
        "AGAGAGTCTCTCAAAGGGCCCTTTAGAGAGTCTCTCAAAGGGCCCTTTAGAGAGTCTCTC",
    };

    // Genome-specific sequences: mutually dissimilar 60-nt, share no 5-mers with each other
    // or with the Shared set, so they yield no ortholog (dispensable genes). (Source: Tettelin)
    private const string Unique1 = "CCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGGGG";
    private const string Unique2 = "TTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGG";
    private const string Unique3 = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT";
    private const string Unique4 = "ACGTAAACCCGGGTTTACGTAAACCCGGGTTTACGTAAACCCGGGTTTACGTAAACCCGG";

    private static List<G> GenomeOf(string prefix, string genomeId, params string[] seqs)
    {
        var list = new List<G>(seqs.Length);
        for (int i = 0; i < seqs.Length; i++)
            list.Add(new G($"{prefix}{i}", genomeId, i * 100, i * 100 + 60, '+', seqs[i]));
        return list;
    }

    #endregion

    #region CompareGenomes — core/dispensable partition

    // M1 — Tettelin (2005): a gene shared by both genomes is CORE (conserved); a gene unique to
    // one genome is that genome's DISPENSABLE (genome-specific) gene. One shared pair + one unique
    // each => Conserved=1, Specific1=1, Specific2=1. A wrong (e.g. one-directional or no-ortholog)
    // implementation would give Conserved=0 / Specific=2, so this value rejects such defects.
    [Test]
    public void CompareGenomes_OneSharedOneUnique_PartitionsCoreAndSpecific()
    {
        var g1 = GenomeOf("a", "G1", Shared[0], Unique1);
        var g2 = GenomeOf("c", "G2", Shared[0], Unique2);

        var result = ComparativeGenomics.CompareGenomes(g1, g2);

        Assert.Multiple(() =>
        {
            Assert.That(result.ConservedGenes, Is.EqualTo(1),
                "One sequence is shared by both genomes -> exactly 1 core gene (Tettelin 2005 core = present in all).");
            Assert.That(result.Orthologs.Count, Is.EqualTo(1),
                "ConservedGenes must equal the RBH ortholog-pair count (INV-01).");
            Assert.That(result.GenomeSpecificGenes1, Is.EqualTo(1),
                "Unique1 has no ortholog -> 1 dispensable gene in genome 1 (Tettelin: unique to each strain).");
            Assert.That(result.GenomeSpecificGenes2, Is.EqualTo(1),
                "Unique2 has no ortholog -> 1 dispensable gene in genome 2.");
            Assert.That(result.ConservedGenes + result.GenomeSpecificGenes1, Is.EqualTo(g1.Count),
                "core + dispensable must equal genome-1 size (INV-02).");
            Assert.That(result.ConservedGenes + result.GenomeSpecificGenes2, Is.EqualTo(g2.Count),
                "core + dispensable must equal genome-2 size (INV-02).");
        });
    }

    // M2 — Tettelin (2005): with no shared genes, the core is empty and every gene is "unique to
    // each strain" (all dispensable). Conserved=0, Specific1=2, Specific2=2.
    [Test]
    public void CompareGenomes_DisjointContent_NoCoreAllSpecific()
    {
        var g1 = GenomeOf("a", "G1", Unique1, Unique3);
        var g2 = GenomeOf("c", "G2", Unique2, Unique4);

        var result = ComparativeGenomics.CompareGenomes(g1, g2);

        Assert.Multiple(() =>
        {
            Assert.That(result.ConservedGenes, Is.EqualTo(0),
                "No sequence is shared -> empty core genome (Tettelin 2005).");
            Assert.That(result.Orthologs, Is.Empty,
                "No reciprocal best hits exist for disjoint content.");
            Assert.That(result.GenomeSpecificGenes1, Is.EqualTo(2),
                "All of genome 1 is genome-specific (dispensable).");
            Assert.That(result.GenomeSpecificGenes2, Is.EqualTo(2),
                "All of genome 2 is genome-specific (dispensable).");
        });
    }

    // C1 — Tettelin (2005): if every gene is shared, the dispensable genome is empty for both.
    [Test]
    public void CompareGenomes_AllGenesShared_NoGenomeSpecificGenes()
    {
        var g1 = GenomeOf("a", "G1", Shared[0], Shared[1]);
        var g2 = GenomeOf("c", "G2", Shared[0], Shared[1]);

        var result = ComparativeGenomics.CompareGenomes(g1, g2);

        Assert.Multiple(() =>
        {
            Assert.That(result.ConservedGenes, Is.EqualTo(2),
                "Both genes shared -> 2 core genes.");
            Assert.That(result.GenomeSpecificGenes1, Is.EqualTo(0),
                "No gene of genome 1 is unique -> empty dispensable set (Tettelin 2005).");
            Assert.That(result.GenomeSpecificGenes2, Is.EqualTo(0),
                "No gene of genome 2 is unique -> empty dispensable set.");
        });
    }

    #endregion

    #region CompareGenomes — synteny fraction & full result

    // M3 — Core partition + "fraction of syntenic genes" metric. Five collinear shared orthologs
    // form one MCScanX block (score 5*50=250 >= 250); plus one unique gene each.
    // Conserved=5, Specific=1/1. OverallSynteny = genes-in-blocks / min(|g1|,|g2|) = 5/6 = 0.8333...
    // Identity permutation => 0 rearrangements. (Sources: Tettelin core; synteny fraction; MCScanX.)
    [Test]
    public void CompareGenomes_IdenticalCollinearContent_ReportsCoreAndSyntenyFraction()
    {
        var g1 = GenomeOf("a", "G1", Shared[0], Shared[1], Shared[2], Shared[3], Shared[4], Unique1);
        var g2 = GenomeOf("c", "G2", Shared[0], Shared[1], Shared[2], Shared[3], Shared[4], Unique2);

        var result = ComparativeGenomics.CompareGenomes(g1, g2);

        const double expectedSynteny = 5.0 / 6.0; // genes in block (5) / smaller genome (6)

        Assert.Multiple(() =>
        {
            Assert.That(result.ConservedGenes, Is.EqualTo(5),
                "Five sequences are shared -> 5 core genes (Tettelin 2005).");
            Assert.That(result.GenomeSpecificGenes1, Is.EqualTo(1),
                "Unique1 is the single dispensable gene of genome 1.");
            Assert.That(result.GenomeSpecificGenes2, Is.EqualTo(1),
                "Unique2 is the single dispensable gene of genome 2.");
            Assert.That(result.SyntenicBlocks.Count, Is.EqualTo(1),
                "Five collinear orthologs form exactly one MCScanX syntenic block (score 250).");
            Assert.That(result.SyntenicBlocks[0].GeneCount, Is.EqualTo(5),
                "The block contains all 5 collinear conserved genes.");
            Assert.That(result.OverallSynteny, Is.EqualTo(expectedSynteny).Within(1e-10),
                "OverallSynteny = fraction of syntenic genes = 5/6 (Synteny overview metric).");
            Assert.That(result.Rearrangements, Is.Empty,
                "Identity gene order has no breakpoints -> no rearrangements.");
        });
    }

    // S1 — MCScanX threshold boundary (Assumption 2): with only 3 collinear orthologs the block
    // score (3*50=150) is below the 250 report rule, so NO block is reported and OverallSynteny=0,
    // even though there ARE conserved orthologs. Conserved still equals the ortholog count.
    [Test]
    public void CompareGenomes_FewCollinearOrthologs_ConservedButZeroSynteny()
    {
        var g1 = GenomeOf("a", "G1", Shared[0], Shared[1], Shared[2], Unique1);
        var g2 = GenomeOf("c", "G2", Shared[0], Shared[1], Shared[2], Unique2);

        var result = ComparativeGenomics.CompareGenomes(g1, g2);

        Assert.Multiple(() =>
        {
            Assert.That(result.ConservedGenes, Is.EqualTo(3),
                "Three sequences are shared -> 3 core genes regardless of block reporting.");
            Assert.That(result.SyntenicBlocks, Is.Empty,
                "3 collinear anchors score 150 < 250 -> no MCScanX block reported (Assumption 2).");
            Assert.That(result.OverallSynteny, Is.EqualTo(0.0).Within(1e-10),
                "No syntenic block -> fraction of syntenic genes is 0.");
        });
    }

    // S2 — INV-04: the RBH matching is symmetric, so swapping the genomes keeps ConservedGenes and
    // swaps Specific1 <-> Specific2. Use an asymmetric layout (2 unique in g1, 1 unique in g2).
    [Test]
    public void CompareGenomes_SwappedGenomes_SwapsGenomeSpecificCounts()
    {
        var g1 = GenomeOf("a", "G1", Shared[0], Shared[1], Unique1, Unique3);
        var g2 = GenomeOf("c", "G2", Shared[0], Shared[1], Unique2);

        var forward = ComparativeGenomics.CompareGenomes(g1, g2);
        var swapped = ComparativeGenomics.CompareGenomes(g2, g1);

        Assert.Multiple(() =>
        {
            Assert.That(forward.ConservedGenes, Is.EqualTo(2),
                "Two shared genes -> 2 core genes.");
            Assert.That(forward.GenomeSpecificGenes1, Is.EqualTo(2),
                "Genome 1 has 2 unique genes.");
            Assert.That(forward.GenomeSpecificGenes2, Is.EqualTo(1),
                "Genome 2 has 1 unique gene.");
            Assert.That(swapped.ConservedGenes, Is.EqualTo(forward.ConservedGenes),
                "Conserved count is invariant under genome swap (INV-04).");
            Assert.That(swapped.GenomeSpecificGenes1, Is.EqualTo(forward.GenomeSpecificGenes2),
                "Swapping genomes swaps the dispensable counts (INV-04).");
            Assert.That(swapped.GenomeSpecificGenes2, Is.EqualTo(forward.GenomeSpecificGenes1),
                "Swapping genomes swaps the dispensable counts (INV-04).");
        });
    }

    #endregion

    #region CompareGenomes — edge cases

    // M4 — empty genomes: no ortholog pairs possible -> all-zero partition, empty collections.
    [Test]
    public void CompareGenomes_EmptyGenomes_ReturnsEmptyPartition()
    {
        var result = ComparativeGenomics.CompareGenomes(new List<G>(), new List<G>());

        Assert.Multiple(() =>
        {
            Assert.That(result.ConservedGenes, Is.EqualTo(0), "No genes -> no core genes.");
            Assert.That(result.GenomeSpecificGenes1, Is.EqualTo(0), "No genes -> no dispensable genes.");
            Assert.That(result.GenomeSpecificGenes2, Is.EqualTo(0), "No genes -> no dispensable genes.");
            Assert.That(result.OverallSynteny, Is.EqualTo(0.0).Within(1e-10), "No genes -> 0 synteny.");
            Assert.That(result.Orthologs, Is.Empty, "No ortholog pairs.");
            Assert.That(result.SyntenicBlocks, Is.Empty, "No syntenic blocks.");
            Assert.That(result.Rearrangements, Is.Empty, "No rearrangements.");
        });
    }

    [Test]
    public void CompareGenomes_NullGenome_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => ComparativeGenomics.CompareGenomes(null!, new List<G>()),
                NUnit.Framework.Throws.ArgumentNullException, "Null genome 1 must throw ArgumentNullException.");
            Assert.That(() => ComparativeGenomics.CompareGenomes(new List<G>(), null!),
                NUnit.Framework.Throws.ArgumentNullException, "Null genome 2 must throw ArgumentNullException.");
        });
    }

    #endregion
}
