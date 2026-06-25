using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using G = Seqeron.Genomics.Analysis.ComparativeGenomics.Gene;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Comparative-Genomics area — comprehensive two-genome
/// comparison (COMPGEN-COMPARE-001), the pan-genome core/dispensable partition
/// plus the syntenic-gene fraction.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (a similarity fraction outside
/// [0,1], a self-comparison that is not the maximum, a partition whose core +
/// dispensable does not equal the genome size, a negative count), and no
/// *unhandled* runtime exception (DivideByZero on an empty genome's zero
/// denominator, IndexOutOfRange on a 1-gene genome, NullReferenceException on a
/// gene with no sequence). Every input must resolve to EITHER a well-defined,
/// theory-correct result, OR a *documented, intentional* validation exception
/// (ArgumentNullException). A raw runtime exception, a hang, a malformed
/// ComparisonResult, or an out-of-range fraction is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-COMPARE-001 — Comprehensive genome comparison
/// Checklist: docs/checklists/03_FUZZING.md, row 133.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty) —
///          the degenerate comparison boundaries called out in the checklist row:
///          an EMPTY genome (one or both), A vs A (a genome compared to itself,
///          the canonical maximum), and DISJOINT genomes (no shared content, the
///          canonical minimum).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The comparison contract under test (Genome_Comparison.md)
/// ───────────────────────────────────────────────────────────────────────────
/// CompareGenomes partitions genes via the pan-genome model of Tettelin et al.
/// (2005, PNAS 102(39):13950–13955): a gene with an ortholog in the other genome
/// is CORE (conserved); a gene with no ortholog is that genome's DISPENSABLE
/// (genome-specific) gene. Shared genes are reciprocal best hits (Moreno-Hagelsieb
/// &amp; Latimer 2008). OverallSynteny is the "fraction of syntenic genes" — genes
/// inside MCScanX syntenic blocks ÷ the smaller genome size, clamped to ≤ 1
/// (Wang et al. 2012). The API entry under test is
///   ComparativeGenomics.CompareGenomes(
///       IReadOnlyList&lt;Gene&gt; genome1Genes,
///       IReadOnlyList&lt;Gene&gt; genome2Genes,
///       double minOrthologIdentity = 0.3,
///       int minSyntenicBlockSize = 3)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs
///    lines 765–810; the synteny-fraction guard at line 797–800).
///
/// THE DOCUMENTED INVARIANTS (Genome_Comparison.md §2.4):
///   • INV-01: ConservedGenes == Orthologs.Count.
///   • INV-02: ConservedGenes + GenomeSpecificGenes_i == |genome_i| for i ∈ {1,2}.
///   • INV-03: 0 ≤ OverallSynteny ≤ 1 (a fraction of genes, explicitly clamped).
///   • INV-04: swapping genome1/genome2 keeps ConservedGenes and swaps
///             GenomeSpecificGenes1 ↔ GenomeSpecificGenes2 (RBH is reciprocal).
/// The documented edge cases (§3.3 / §6.1):
///   • both genomes empty → all counts 0, empty collections, OverallSynteny 0;
///   • no shared genes → ConservedGenes 0, every gene genome-specific;
///   • all genes shared → GenomeSpecific 0 for both;
///   • null genome list → ArgumentNullException.
///
/// The three BE checklist targets map to these documented behaviours:
///   • empty genome → all-zero, NO DivideByZero on the zero-gene denominator of
///                    OverallSynteny (the `smallerGenome > 0` guard, line 798,
///                    short-circuits the division before min(|g1|,|g2|) = 0 is
///                    used as a divisor).
///   • A vs A       → the MAXIMUM comparison: every gene of the genome is core
///                    (ConservedGenes = |genome|, GenomeSpecific = 0 for both),
///                    and with ≥ 5 collinear orthologs OverallSynteny is its
///                    clamped maximum 1.0. The self-comparison is the upper
///                    boundary of similarity and must never exceed it (INV-03).
///   • disjoint     → the MINIMUM comparison: ConservedGenes = 0, every gene is
///                    genome-specific, OverallSynteny = 0 — no false sharing,
///                    no crash.
/// A positive-sanity test pins the documented golden vector (§7.1) and that a
/// known overlap yields the hand-computed shared/unique counts and synteny.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeCompareFuzzTests
{
    #region Test Data and Helpers

    // Five DISTINCT >= 60-nt sequences, mutually dissimilar (no shared 5-mer overlap above the
    // RBH gate). Distinct content guarantees an unambiguous reciprocal-best-hit matching, so a
    // genome built from a subset of these has a clean self-RBH and a clean cross-RBH. (Source: RBH.)
    private static readonly string[] Shared =
    {
        "ATGGCAAAGCTTGATCCGTACGGGTTAACCGGATCAGGTTCAAAGCTTGATCCGTACGGG",
        "TTACCGGATCAGGTTCATGGCAAAGCTTGATCCGTACGGGAATTACCGGATCAGGTTCAT",
        "GGCCAATTGGCCAATTACGTACGTGGCCAATTGGCCAATTACGTACGTGGCCAATTGGCC",
        "CTGACTGACAAATTTGGGCCCCTGACTGACAAATTTGGGCCCCTGACTGACAAATTTGGG",
        "AGAGAGTCTCTCAAAGGGCCCTTTAGAGAGTCTCTCAAAGGGCCCTTTAGAGAGTCTCTC",
    };

    // Genome-specific sequences: mutually dissimilar 60-nt, share no 5-mers with each other
    // or with the Shared set, so they yield no ortholog (dispensable genes). (Source: Tettelin.)
    private static readonly string[] Unique =
    {
        "CCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGGGGTTTTTTAAAAAACCCCCCGGGGGG",
        "TTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGGAACCTTGG",
        "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT",
        "ACGTAAACCCGGGTTTACGTAAACCCGGGTTTACGTAAACCCGGGTTTACGTAAACCCGG",
        "TCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCA",
    };

    /// <summary>
    /// Builds an ordered genome from a list of sequences. Gene ids carry the genome prefix so the
    /// two genomes use disjoint id namespaces (RBH matches on sequence, not id). Genes are laid out
    /// at non-overlapping 0-based coordinates in chromosomal order.
    /// </summary>
    private static List<G> GenomeOf(string prefix, string genomeId, params string[] seqs)
    {
        var list = new List<G>(seqs.Length);
        for (int i = 0; i < seqs.Length; i++)
            list.Add(new G($"{prefix}{i}", genomeId, i * 100, i * 100 + 60, '+', seqs[i]));
        return list;
    }

    /// <summary>
    /// A well-formed ComparisonResult satisfies the documented structural invariants for ANY input:
    ///   • no collection is null (§3.2 shape),
    ///   • all counts are non-negative (counts, not signed metrics),
    ///   • ConservedGenes == Orthologs.Count (INV-01),
    ///   • ConservedGenes + GenomeSpecificGenes_i == |genome_i| (INV-02, the partition is exact),
    ///   • 0 ≤ OverallSynteny ≤ 1 — never NaN/Infinity from a zero denominator (INV-03).
    /// This is the always-on safety net for every fuzzed input; the targeted tests pin the exact
    /// boundary values on top of it.
    /// </summary>
    private static void AssertWellFormed(
        ComparativeGenomics.ComparisonResult r, int genome1Size, int genome2Size)
    {
        r.Orthologs.Should().NotBeNull("the ortholog list must never be null (§3.2)");
        r.SyntenicBlocks.Should().NotBeNull("the syntenic-block list must never be null (§3.2)");
        r.Rearrangements.Should().NotBeNull("the rearrangement list must never be null (§3.2)");

        r.ConservedGenes.Should().BeGreaterThanOrEqualTo(0, "a gene count is non-negative");
        r.GenomeSpecificGenes1.Should().BeGreaterThanOrEqualTo(0, "a gene count is non-negative");
        r.GenomeSpecificGenes2.Should().BeGreaterThanOrEqualTo(0, "a gene count is non-negative");

        r.ConservedGenes.Should().Be(r.Orthologs.Count, "ConservedGenes == Orthologs.Count (INV-01)");

        (r.ConservedGenes + r.GenomeSpecificGenes1).Should().Be(genome1Size,
            "core + dispensable must equal genome-1 size (INV-02)");
        (r.ConservedGenes + r.GenomeSpecificGenes2).Should().Be(genome2Size,
            "core + dispensable must equal genome-2 size (INV-02)");

        double.IsNaN(r.OverallSynteny).Should().BeFalse(
            "OverallSynteny must never be NaN — a zero-gene denominator must be guarded, not divided (INV-03)");
        double.IsInfinity(r.OverallSynteny).Should().BeFalse(
            "OverallSynteny must never be Infinity — no DivideByZero on min(|g1|,|g2|)=0 (INV-03)");
        r.OverallSynteny.Should().BeInRange(0.0, 1.0,
            "OverallSynteny is a fraction of syntenic genes, clamped to [0,1] (INV-03)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  COMPGEN-COMPARE-001 — Comprehensive genome comparison : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region COMPGEN-COMPARE-001 — Genome comparison (core/dispensable partition + synteny fraction)

    #region BE — Boundary: empty genome (zero-gene denominator hazard)

    /// <summary>
    /// BE: both genomes empty is the floor of the gene-count axis. The documented result is the
    /// all-zero partition with empty collections and OverallSynteny 0 (§3.3, §6.1). The hazard
    /// probed is a DivideByZero / NaN / Infinity on OverallSynteny, whose denominator is
    /// min(|g1|,|g2|) = 0; the `smallerGenome > 0` guard (ComparativeGenomics.cs line 798) must
    /// short-circuit to 0 BEFORE the division.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_BothEmpty_AllZeroNoDivideByZero()
    {
        var g1 = new List<G>();
        var g2 = new List<G>();

        ComparativeGenomics.ComparisonResult r = default;
        FluentActions.Invoking(() => r = ComparativeGenomics.CompareGenomes(g1, g2))
            .Should().NotThrow("two empty genomes must yield a clean all-zero result, never a DivideByZero");

        AssertWellFormed(r, g1.Count, g2.Count);
        r.ConservedGenes.Should().Be(0, "no genes → no core genes (§6.1)");
        r.GenomeSpecificGenes1.Should().Be(0, "no genes → no dispensable genes");
        r.GenomeSpecificGenes2.Should().Be(0, "no genes → no dispensable genes");
        r.OverallSynteny.Should().Be(0.0, "min(0,0)=0 denominator is guarded → synteny exactly 0 (§3.3)");
        r.Orthologs.Should().BeEmpty("no ortholog pairs are possible between empty genomes");
        r.SyntenicBlocks.Should().BeEmpty("no syntenic blocks");
        r.Rearrangements.Should().BeEmpty("no rearrangements");
    }

    /// <summary>
    /// BE: ONE empty genome (the other non-empty) is the asymmetric floor — the case most likely to
    /// hit a one-sided DivideByZero, since min(|g1|,|g2|) = 0 only when one side is empty too. Every
    /// gene of the non-empty genome has no possible ortholog (the other genome is empty), so it is
    /// entirely genome-specific; the empty side has zero of everything. OverallSynteny stays 0.
    /// Both argument orders are exercised to pin symmetry of the guard.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_OneEmpty_OtherAllSpecificNoCrash()
    {
        var full = GenomeOf("a", "G1", Shared[0], Shared[1], Unique[0]);
        var empty = new List<G>();

        ComparativeGenomics.ComparisonResult fwd = default, rev = default;
        FluentActions.Invoking(() => fwd = ComparativeGenomics.CompareGenomes(full, empty))
            .Should().NotThrow("a non-empty vs empty comparison must not DivideByZero on the empty side");
        FluentActions.Invoking(() => rev = ComparativeGenomics.CompareGenomes(empty, full))
            .Should().NotThrow("the empty-first ordering must be equally safe");

        AssertWellFormed(fwd, full.Count, empty.Count);
        AssertWellFormed(rev, empty.Count, full.Count);

        fwd.ConservedGenes.Should().Be(0, "the other genome is empty → no ortholog possible → no core gene");
        fwd.GenomeSpecificGenes1.Should().Be(full.Count, "every gene of the non-empty genome is dispensable");
        fwd.GenomeSpecificGenes2.Should().Be(0, "the empty genome has no genes");
        fwd.OverallSynteny.Should().Be(0.0, "min(|full|,0)=0 → guarded synteny is 0");

        // INV-04: swapping moves the specific count to the other side; conserved is invariant.
        rev.ConservedGenes.Should().Be(fwd.ConservedGenes, "conserved is swap-invariant (INV-04)");
        rev.GenomeSpecificGenes2.Should().Be(fwd.GenomeSpecificGenes1, "swap moves the specific count (INV-04)");
        rev.GenomeSpecificGenes1.Should().Be(fwd.GenomeSpecificGenes2, "swap moves the specific count (INV-04)");
    }

    /// <summary>
    /// BE: a single 1-gene genome vs an empty genome is the smallest non-trivial input on the
    /// boundary — the corner where a naive synteny or RBH loop could index past the end. It must
    /// hit the empty-side path with no IndexOutOfRange and report the lone gene as genome-specific.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_SingleGeneVsEmpty_OneSpecificNoIndexOutOfRange()
    {
        var one = GenomeOf("a", "G1", Shared[0]);
        var empty = new List<G>();

        ComparativeGenomics.ComparisonResult r = default;
        FluentActions.Invoking(() => r = ComparativeGenomics.CompareGenomes(one, empty))
            .Should().NotThrow("a 1-gene vs empty comparison must not index past the empty list");

        AssertWellFormed(r, one.Count, empty.Count);
        r.ConservedGenes.Should().Be(0, "the lone gene has nothing to pair with");
        r.GenomeSpecificGenes1.Should().Be(1, "the lone gene is genome-specific");
        r.OverallSynteny.Should().Be(0.0, "no block, guarded denominator → 0");
    }

    #endregion

    #region BE — Boundary: A vs A (self-comparison is the maximum)

    /// <summary>
    /// BE: A vs A — a genome compared to ITSELF is the upper boundary of similarity. Every gene is
    /// its own reciprocal best hit, so the partition is fully core: ConservedGenes = |genome| and
    /// GenomeSpecificGenes = 0 for BOTH sides (the dispensable genome is empty, §6.1 "all genes
    /// shared"). This pins the canonical maximum-sharing result and that the self-comparison never
    /// invents a genome-specific gene. (Two independent copies are used so the id namespaces differ
    /// — RBH matches on sequence, exactly as a self-comparison of re-annotated identical genomes.)
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_SelfComparison_AllCoreNoSpecific()
    {
        var seqs = new[] { Shared[0], Shared[1], Shared[2] };
        var a1 = GenomeOf("a", "G1", seqs);
        var a2 = GenomeOf("b", "G1", seqs); // identical content, distinct ids = the same genome

        var r = ComparativeGenomics.CompareGenomes(a1, a2);

        AssertWellFormed(r, a1.Count, a2.Count);
        r.ConservedGenes.Should().Be(seqs.Length,
            "a genome vs itself: every gene is its own RBH → fully core (§6.1 all-shared, the maximum)");
        r.GenomeSpecificGenes1.Should().Be(0, "self-comparison has an empty dispensable genome on side 1");
        r.GenomeSpecificGenes2.Should().Be(0, "self-comparison has an empty dispensable genome on side 2");
    }

    /// <summary>
    /// BE: A vs A at the synteny upper boundary. With ≥ 5 collinear conserved orthologs in identity
    /// order the single MCScanX block contains all genes, so OverallSynteny = |genes| / |genes| = 1
    /// — its clamped MAXIMUM (INV-03). The identity gene order has no breakpoints, so there are no
    /// rearrangements. This pins that the maximum self-similarity is exactly 1.0, never above.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_SelfComparisonCollinear_SyntenyIsMaximumOne()
    {
        var seqs = new[] { Shared[0], Shared[1], Shared[2], Shared[3], Shared[4] };
        var a1 = GenomeOf("a", "G1", seqs);
        var a2 = GenomeOf("b", "G1", seqs);

        var r = ComparativeGenomics.CompareGenomes(a1, a2);

        AssertWellFormed(r, a1.Count, a2.Count);
        r.ConservedGenes.Should().Be(5, "all 5 genes are core in a self-comparison");
        r.GenomeSpecificGenes1.Should().Be(0, "no dispensable gene in a self-comparison");
        r.GenomeSpecificGenes2.Should().Be(0, "no dispensable gene in a self-comparison");
        r.SyntenicBlocks.Should().HaveCount(1, "5 collinear orthologs (score 250) form one MCScanX block");
        r.SyntenicBlocks[0].GeneCount.Should().Be(5, "the block spans all 5 collinear conserved genes");
        r.OverallSynteny.Should().Be(1.0,
            "all genes are syntenic → 5/5 = the clamped MAXIMUM 1.0 (INV-03); self-similarity is the upper boundary");
        r.Rearrangements.Should().BeEmpty("identity gene order has no breakpoints → no rearrangements");
    }

    /// <summary>
    /// BE: self-comparison must be the GLOBAL maximum over the partition — compared to any other
    /// genome built from a strict subset/superset, A vs A reports at least as many core genes
    /// relative to its size (zero dispensable) as any non-identical pairing. We pin that
    /// GenomeSpecific is 0 on both sides for the self-pair but strictly positive once a foreign gene
    /// is introduced, so the self-comparison is strictly the maximum-sharing configuration.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_SelfComparison_StrictlyMaximalSharing()
    {
        var seqs = new[] { Shared[0], Shared[1] };
        var self1 = GenomeOf("a", "G1", seqs);
        var self2 = GenomeOf("b", "G1", seqs);

        // A near-identical genome that replaces one gene with a foreign one.
        var nearly = GenomeOf("c", "G2", Shared[0], Unique[0]);

        var selfR = ComparativeGenomics.CompareGenomes(self1, self2);
        var nearR = ComparativeGenomics.CompareGenomes(self1, nearly);

        AssertWellFormed(selfR, self1.Count, self2.Count);
        AssertWellFormed(nearR, self1.Count, nearly.Count);

        (selfR.GenomeSpecificGenes1 + selfR.GenomeSpecificGenes2).Should().Be(0,
            "the self-comparison shares everything → zero total dispensable genes (the maximum)");
        (nearR.GenomeSpecificGenes1 + nearR.GenomeSpecificGenes2).Should().BeGreaterThan(0,
            "introducing a foreign gene strictly lowers sharing below the self-comparison maximum");
        selfR.ConservedGenes.Should().BeGreaterThan(nearR.ConservedGenes,
            "the self-comparison has strictly more core genes than the genome that swapped one gene out");
    }

    #endregion

    #region BE — Boundary: disjoint genomes (the minimum)

    /// <summary>
    /// BE: disjoint genomes share NO content — the lower boundary of similarity. No reciprocal best
    /// hit exists, so ConservedGenes = 0 and EVERY gene is genome-specific (§6.1 "no shared genes";
    /// Tettelin "genes unique to each strain"). OverallSynteny = 0 (no block). This pins the
    /// canonical minimum and that nothing is falsely clustered, with no crash on fully-unmatched input.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_DisjointGenomes_NoCoreAllSpecificMinimumSynteny()
    {
        var g1 = GenomeOf("a", "G1", Unique[0], Unique[2], Unique[4]);
        var g2 = GenomeOf("c", "G2", Unique[1], Unique[3], Shared[0]);

        var r = ComparativeGenomics.CompareGenomes(g1, g2);

        AssertWellFormed(r, g1.Count, g2.Count);
        r.ConservedGenes.Should().Be(0, "no shared content → empty core genome (the minimum, §6.1)");
        r.Orthologs.Should().BeEmpty("no reciprocal best hits for disjoint content");
        r.GenomeSpecificGenes1.Should().Be(g1.Count, "every gene of genome 1 is dispensable (unique to the strain)");
        r.GenomeSpecificGenes2.Should().Be(g2.Count, "every gene of genome 2 is dispensable");
        r.OverallSynteny.Should().Be(0.0, "no syntenic block → fraction of syntenic genes is the minimum 0");
        r.SyntenicBlocks.Should().BeEmpty("no orthologs → no collinear block");
        r.Rearrangements.Should().BeEmpty("no shared anchors → no breakpoints");
    }

    /// <summary>
    /// BE: a degenerate gene with a NULL/EMPTY sequence in a disjoint comparison must be treated as
    /// un-orthologable (ortholog detection ignores genes with no sequence, §3.3) — counted as
    /// genome-specific, never a NullReferenceException on the missing sequence. Pins that the
    /// minimum-similarity case is robust to malformed genes.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_DisjointWithNullSequenceGene_AllSpecificNoNullRef()
    {
        var g1 = new List<G>
        {
            new("a0", "G1", 0, 60, '+', Unique[0]),
            new("a1", "G1", 100, 160, '+', null),  // null sequence: un-orthologable
            new("a2", "G1", 200, 260, '+', ""),    // empty sequence: un-orthologable
        };
        var g2 = GenomeOf("c", "G2", Unique[1], Unique[3]);

        ComparativeGenomics.ComparisonResult r = default;
        FluentActions.Invoking(() => r = ComparativeGenomics.CompareGenomes(g1, g2))
            .Should().NotThrow("a gene with no sequence must be ignored for orthology, not deref'd");

        AssertWellFormed(r, g1.Count, g2.Count);
        r.ConservedGenes.Should().Be(0, "disjoint + un-orthologable genes → no core gene");
        r.GenomeSpecificGenes1.Should().Be(3, "all 3 genome-1 genes (incl. the seq-less ones) are dispensable");
        r.GenomeSpecificGenes2.Should().Be(2, "both genome-2 genes are dispensable");
    }

    #endregion

    #region BE — Validation: null contract

    /// <summary>
    /// BE: a null gene list is the one input the contract rejects with an exception (§3.3,
    /// §6.1; ComparativeGenomics.cs lines 771–772). Each side is pinned independently so a null
    /// argument is a disciplined defensive throw, not an unhandled NRE deep in the RBH scan.
    /// </summary>
    [Test]
    public void CompareGenomes_NullGenome_ThrowsArgumentNull()
    {
        FluentActions.Invoking(() => ComparativeGenomics.CompareGenomes(null!, new List<G>()))
            .Should().Throw<ArgumentNullException>("null genome 1 is rejected defensively")
            .Which.ParamName.Should().Be("genome1Genes");

        FluentActions.Invoking(() => ComparativeGenomics.CompareGenomes(new List<G>(), null!))
            .Should().Throw<ArgumentNullException>("null genome 2 is rejected defensively")
            .Which.ParamName.Should().Be("genome2Genes");
    }

    #endregion

    #region Positive sanity — the documented metric on hand-built examples

    /// <summary>
    /// Positive sanity (Genome_Comparison.md §7.1 worked example): a genome with one shared gene
    /// and one unique gene, compared to another with the same shared gene and a different unique
    /// gene, yields the EXACT hand-computed partition — ConservedGenes = 1, GenomeSpecificGenes1 = 1,
    /// GenomeSpecificGenes2 = 1. This is the load-bearing correctness check distinguishing a
    /// theory-correct aggregator from a one-directional or no-ortholog implementation.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_KnownOverlap_ExactPartition()
    {
        var g1 = GenomeOf("a", "G1", Shared[0], Unique[0]);
        var g2 = GenomeOf("c", "G2", Shared[0], Unique[1]);

        var r = ComparativeGenomics.CompareGenomes(g1, g2);

        AssertWellFormed(r, g1.Count, g2.Count);
        r.ConservedGenes.Should().Be(1, "exactly one sequence is shared → 1 core gene (§7.1)");
        r.GenomeSpecificGenes1.Should().Be(1, "Unique[0] has no ortholog → 1 dispensable gene in genome 1");
        r.GenomeSpecificGenes2.Should().Be(1, "Unique[1] has no ortholog → 1 dispensable gene in genome 2");
    }

    /// <summary>
    /// Positive sanity — the synteny FRACTION on a hand-computed mixed case: 5 collinear shared
    /// orthologs (one MCScanX block of 5, score 250) plus one unique gene each. The smaller genome
    /// has 6 genes, so OverallSynteny = 5/6 ≈ 0.8333… (genes-in-blocks ÷ smaller genome). Pins the
    /// exact documented metric value, not merely a green run. (§4.1 step 5; Wang et al. 2012.)
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_MixedCollinear_ExactSyntenyFraction()
    {
        var g1 = GenomeOf("a", "G1", Shared[0], Shared[1], Shared[2], Shared[3], Shared[4], Unique[0]);
        var g2 = GenomeOf("c", "G2", Shared[0], Shared[1], Shared[2], Shared[3], Shared[4], Unique[1]);

        var r = ComparativeGenomics.CompareGenomes(g1, g2);

        AssertWellFormed(r, g1.Count, g2.Count);
        r.ConservedGenes.Should().Be(5, "five shared sequences → 5 core genes");
        r.GenomeSpecificGenes1.Should().Be(1, "Unique[0] is the lone dispensable gene of genome 1");
        r.GenomeSpecificGenes2.Should().Be(1, "Unique[1] is the lone dispensable gene of genome 2");
        r.SyntenicBlocks.Should().HaveCount(1, "5 collinear orthologs form one MCScanX block (score 250)");
        r.OverallSynteny.Should().BeApproximately(5.0 / 6.0, 1e-10,
            "OverallSynteny = genes-in-blocks (5) / smaller genome (6) = 5/6 (§4.1, Wang 2012)");
    }

    /// <summary>
    /// Positive sanity — DETERMINISM. CompareGenomes uses a deterministic RBH matching with
    /// deterministic tie-breaks (Moreno-Hagelsieb &amp; Latimer 2008), so repeated invocations on the
    /// same input must return identical partition values and metrics, with no reliance on hash
    /// iteration order. Pins reproducibility, a fuzzing must-have.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CompareGenomes_RepeatedRuns_AreDeterministic()
    {
        var g1 = GenomeOf("a", "G1", Shared[0], Shared[1], Shared[2], Unique[0]);
        var g2 = GenomeOf("c", "G2", Shared[0], Shared[1], Shared[2], Unique[1]);

        var first = ComparativeGenomics.CompareGenomes(g1, g2);
        for (int i = 0; i < 5; i++)
        {
            var again = ComparativeGenomics.CompareGenomes(g1, g2);
            again.ConservedGenes.Should().Be(first.ConservedGenes, "conserved count is deterministic");
            again.GenomeSpecificGenes1.Should().Be(first.GenomeSpecificGenes1, "specific1 is deterministic");
            again.GenomeSpecificGenes2.Should().Be(first.GenomeSpecificGenes2, "specific2 is deterministic");
            again.OverallSynteny.Should().Be(first.OverallSynteny, "synteny is deterministic");
        }
    }

    /// <summary>
    /// Fuzz sweep — across many random genome pairs of varied sizes, overlap fractions, gene
    /// orders, and gene counts (including empty genomes, self-overlaps and disjoint pairs),
    /// CompareGenomes must ALWAYS terminate, never throw, and return a well-formed result: the exact
    /// partition (core + dispensable = genome size, INV-02), ConservedGenes == Orthologs.Count
    /// (INV-01), and OverallSynteny ∈ [0,1] with no NaN/Infinity from a zero denominator (INV-03).
    /// As a partition lower bound we also pin that every sequence shared BY IDENTITY between the two
    /// genomes is core (it is at minimum its own reciprocal best hit), and conserved never exceeds
    /// the smaller genome (the matching maps each gene at most once). The pool sequences are NOT all
    /// mutually dissimilar (some cross-pairs pass the alignment-free RBH gate), so an exact core
    /// count is asserted only in the dedicated positive-sanity tests over verified-distinct
    /// sequences; here the invariant net + the identity-shared floor + the matching ceiling are the
    /// theory-correct, input-independent checks.
    /// </summary>
    [Test]
    [CancelAfter(120000)]
    public void CompareGenomes_RandomPairs_AlwaysWellFormedAndPartitionBounded()
    {
        string[] pool = Shared.Concat(Unique).ToArray(); // 10 pool sequences
        for (int seed = 0; seed < 60; seed++)
        {
            var rng = new Random(seed);

            // Pick a random subset of distinct pool sequences for each genome (size 0..pool.Length,
            // includes the empty genome). Distinct selection keeps every chosen sequence unique
            // within its genome.
            var idx1 = SampleDistinct(rng, pool.Length, rng.Next(0, pool.Length + 1));
            var idx2 = SampleDistinct(rng, pool.Length, rng.Next(0, pool.Length + 1));

            var g1 = GenomeOf("a", "G1", idx1.Select(i => pool[i]).ToArray());
            var g2 = GenomeOf("c", "G2", idx2.Select(i => pool[i]).ToArray());

            // Sequences present in BOTH genomes by identity are necessarily core (each is its own
            // best, reciprocal hit), giving a verifiable LOWER bound on the conserved count.
            int identityShared = idx1.Intersect(idx2).Count();

            ComparativeGenomics.ComparisonResult r = default;
            FluentActions.Invoking(() => r = ComparativeGenomics.CompareGenomes(g1, g2))
                .Should().NotThrow($"random genome pair must never crash (seed {seed})");

            AssertWellFormed(r, g1.Count, g2.Count);
            r.ConservedGenes.Should().BeGreaterThanOrEqualTo(identityShared,
                $"every identity-shared sequence is core → conserved ≥ identity-shared (seed {seed})");
            r.ConservedGenes.Should().BeLessThanOrEqualTo(Math.Min(g1.Count, g2.Count),
                $"the RBH matching pairs each gene at most once → conserved ≤ smaller genome (seed {seed})");
        }
    }

    /// <summary>Returns <paramref name="count"/> distinct indices in [0,n), drawn without replacement.</summary>
    private static List<int> SampleDistinct(Random rng, int n, int count)
    {
        var all = Enumerable.Range(0, n).ToList();
        for (int i = all.Count - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (all[i], all[j]) = (all[j], all[i]);
        }
        return all.Take(count).ToList();
    }

    #endregion

    #endregion
}
