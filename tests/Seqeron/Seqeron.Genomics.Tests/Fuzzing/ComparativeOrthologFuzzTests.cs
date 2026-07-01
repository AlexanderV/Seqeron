using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Comparative-Genomics area — Ortholog Identification
/// (COMPGEN-ORTHO-001), the Reciprocal-Best-Hit (RBH/BBH) ortholog detector.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (a false ortholog pair, a gene
/// appearing in two pairs, a non-reciprocal pair), and no *unhandled* runtime
/// exception (NullReferenceException on an empty/null gene set, DivideByZero in the
/// similarity normalization, IndexOutOfRange on a sub-k sequence). Every input must
/// resolve to EITHER a well-defined, theory-correct result, OR a *documented,
/// intentional* validation exception (ArgumentNullException on a null gene list). A
/// raw runtime exception, a hang, or a false/non-reciprocal ortholog pair is a bug,
/// not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-ORTHO-001 — Ortholog Identification (Reciprocal Best Hits)
/// Checklist: docs/checklists/03_FUZZING.md, row 135.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate gene-set boundaries called out in
///          the checklist row: no homologs (no similar pair above threshold),
///          all-identical genes, and the empty gene set (one/both sides empty).
///   • MC = Malformed Content — genes carrying malformed/absent payload: null or
///          empty sequence, sub-k-mer (length &lt; 5) sequence, and non-AA/non-DNA
///          junk characters in the sequence string.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The ortholog contract under test (Ortholog_Identification.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Orthologs are detected by the Reciprocal-Best-Hit criterion (Moreno-Hagelsieb &amp;
/// Latimer 2008; Tatusov et al. 1997): a pair (g1 ∈ G1, g2 ∈ G2) is an ortholog pair
/// iff BH(g1, G2) = g2 AND BH(g2, G1) = g1 (§2.2). A one-directional best hit is NOT
/// an ortholog (§2.2, §6.1). Best hit BH(x, T) = argmax_{y∈T} sim(x, y) among targets
/// whose hit QUALIFIES: sim ≥ minIdentity (0.3) AND coverage ≥ minCoverage (0.5)
/// (§3.1). Ties are broken deterministically: higher score, then higher coverage,
/// then SMALLER ordinal gene id (§4.2) — so the best hit is unique. The score is the
/// alignment-free 5-mer Jaccard similarity; coverage = shared 5-mers / 5-mers of the
/// shorter sequence (§4.2, §5.4 Assumption 1). The API entry under test is
///   ComparativeGenomics.FindOrthologs(IReadOnlyList&lt;Gene&gt; genome1Genes,
///       IReadOnlyList&lt;Gene&gt; genome2Genes,
///       double minIdentity = 0.3, double minCoverage = 0.5)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs
///    lines 334-341), which delegates verbatim to FindReciprocalBestHits
///   (lines 465-516); best hit in the private FindBestHit (lines 410-440);
///   similarity in CalculateSequenceSimilarity (lines 518-549).
///
/// THE DOCUMENTED INVARIANTS (Ortholog_Identification.md §2.4):
///   • INV-01: every ortholog pair is reciprocal (RBH is symmetric by construction).
///   • INV-02: output is a MATCHING — no gene id appears in two pairs (best hit is
///             unique per gene via the deterministic argmax + ordinal tie-break).
///   • INV-04: output is deterministic and order-independent.
///   • INV-05: pairs below minIdentity OR minCoverage are excluded (qualification gate).
/// Every positive-result test pins the documented metric (the reciprocity, the
/// Jaccard/coverage values on a hand-built example, the matching property); this is
/// the load-bearing correctness check that distinguishes a true RBH ortholog from a
/// false or one-directional pairing.
///
/// Documented degenerate / validation contract (Ortholog_Identification.md §3.3, §6.1;
/// ComparativeGenomics.cs lines 471-478, 474-475, 524-525):
///   • Null gene list → ArgumentNullException (ThrowIfNull, lines 471-472).
///   • Genes whose Sequence is null/empty are SKIPPED before pairing (the
///     `Where(g => !string.IsNullOrEmpty(g.Sequence))` filter, lines 474-475).
///   • Empty genome (no sequenced genes on either side) → no orthologs; the
///     `Count == 0` short-circuit (lines 477-478) returns an empty result with no
///     NullReference / DivideByZero.
///   • Sequences shorter than the 5-mer length contribute no k-mers and never qualify
///     (the `seq.Length &lt; k` guard returns (0,0,0), lines 524-525) — no false pair,
///     and no DivideByZero because the `total &gt; 0`/`minKmerCount &gt; 0` guards
///     (lines 540, 545) cover the empty-k-mer-set case.
///   • Sequences are upper-cased for k-mer extraction (case-insensitive, §3.3).
///
/// The BE/MC checklist targets map to these documented behaviours:
///   • no homologs → empty result (INV-05): every cross pair scores below the gate,
///                   so no best hit qualifies and no reciprocal pair forms. No crash,
///                   no false positive.
///   • all-identical → a deterministic MATCHING, NOT N² pairs and NOT a hang: with N
///                   identical genes per genome, every gene's best hit is the
///                   smallest-ordinal-id gene (tie-break, §4.2), so reciprocity holds
///                   only for the two smallest-id genes → exactly ONE pair (INV-02).
///                   [CancelAfter] guards against any quadratic blow-up.
///   • empty gene set → empty result (one/both sides empty), no DivideByZero / NRE.
///   • malformed content → genes with null/empty/sub-k/junk sequences are skipped or
///                   simply never qualify; no crash.
/// A positive-sanity test pins the documented §7.1 worked example (two identical
/// genes across genomes → one ortholog pair, Jaccard 1.0) and pins that a divergent,
/// unrelated gene is NOT paired, and that a one-directional best hit is NOT returned.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeOrthologFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    private static ComparativeGenomics.Gene Gene(string id, string genomeId, string? sequence) =>
        new(id, genomeId, 0, (sequence?.Length ?? 0), '+', sequence);

    /// <summary>
    /// A well-formed ortholog result (Ortholog_Identification.md §2.4): a finite,
    /// non-null set of pairs that is a proper MATCHING — every pair is distinct, no
    /// gene-1 id appears twice (INV-02), no gene-2 id appears twice (INV-02), and
    /// every reported pair's metrics are finite and in [0,1] (INV-05 gate). The pairs
    /// are also validated to only reference genes from the supplied sets.
    /// </summary>
    private static void AssertWellFormedMatching(
        IReadOnlyList<ComparativeGenomics.OrthologPair> pairs,
        IEnumerable<ComparativeGenomics.Gene> genome1,
        IEnumerable<ComparativeGenomics.Gene> genome2)
    {
        pairs.Should().NotBeNull("FindOrthologs must always return a (possibly empty) collection, never null");

        var ids1 = new HashSet<string>(genome1.Select(g => g.Id));
        var ids2 = new HashSet<string>(genome2.Select(g => g.Id));

        pairs.Select(p => p.Gene1Id).Should().OnlyHaveUniqueItems(
            "no genome-1 gene may appear in two ortholog pairs — RBH output is a matching (INV-02)");
        pairs.Select(p => p.Gene2Id).Should().OnlyHaveUniqueItems(
            "no genome-2 gene may appear in two ortholog pairs — RBH output is a matching (INV-02)");

        foreach (var p in pairs)
        {
            ids1.Should().Contain(p.Gene1Id, "an ortholog pair may only reference a genome-1 gene");
            ids2.Should().Contain(p.Gene2Id, "an ortholog pair may only reference a genome-2 gene");

            double.IsNaN(p.Identity).Should().BeFalse("identity must never be NaN (the total>0 Jaccard guard prevents 0/0)");
            double.IsNaN(p.Coverage).Should().BeFalse("coverage must never be NaN (the minKmerCount>0 guard prevents 0/0)");
            p.Identity.Should().BeInRange(0.3, 1.0, "a returned pair must have cleared the minIdentity (0.3) gate and Jaccard ≤ 1 (INV-05)");
            p.Coverage.Should().BeInRange(0.5, 1.0, "a returned pair must have cleared the minCoverage (0.5) gate and coverage ≤ 1 (INV-05)");
        }
    }

    /// <summary>
    /// Asserts the RBH output is genuinely reciprocal (INV-01): for every reported
    /// pair (a, b), running FindOrthologs with the two genomes SWAPPED must report the
    /// symmetric pair (b, a). A one-directional best hit would survive one orientation
    /// but not its swap, so this is the load-bearing reciprocity check.
    /// </summary>
    private static void AssertReciprocal(
        IReadOnlyList<ComparativeGenomics.Gene> genome1,
        IReadOnlyList<ComparativeGenomics.Gene> genome2,
        double minIdentity = 0.3,
        double minCoverage = 0.5)
    {
        var forward = ComparativeGenomics.FindOrthologs(genome1, genome2, minIdentity, minCoverage).ToList();
        var swapped = ComparativeGenomics.FindOrthologs(genome2, genome1, minIdentity, minCoverage).ToList();

        var forwardSet = forward.Select(p => (p.Gene1Id, p.Gene2Id)).ToHashSet();
        var swappedSet = swapped.Select(p => (p.Gene2Id, p.Gene1Id)).ToHashSet();

        forwardSet.Should().BeEquivalentTo(swappedSet,
            "RBH is symmetric (INV-01): swapping the two genomes must yield the mirror-image matching");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  COMPGEN-ORTHO-001 — Ortholog Identification (RBH) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region COMPGEN-ORTHO-001 — Ortholog Identification (Reciprocal Best Hits)

    #region BE — Boundary: no homologs (no similar pair above the gate → no false pairs)

    /// <summary>
    /// BE: two genomes with NO conserved sequence are the lower boundary. Every cross
    /// pair (poly-A genes vs poly-C genes) shares zero 5-mers, so Jaccard = 0 &lt; 0.3
    /// and coverage = 0 &lt; 0.5 — no hit qualifies, no best hit exists, no reciprocal
    /// pair forms. We pin an EMPTY result (INV-05): the unit must never invent a false
    /// ortholog from unrelated genes, and the empty-k-mer-intersection must not
    /// DivideByZero (the total&gt;0 / minKmerCount&gt;0 guards return a clean 0).
    /// </summary>
    [Test]
    public void FindOrthologs_NoHomologs_ReturnsEmptyNoFalsePairs()
    {
        var genome1 = new[]
        {
            Gene("a1", "G1", new string('A', 60)),
            Gene("a2", "G1", new string('A', 80)),
        };
        var genome2 = new[]
        {
            Gene("b1", "G2", new string('C', 60)),
            Gene("b2", "G2", new string('C', 80)),
        };

        var orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2).ToList();

        AssertWellFormedMatching(orthologs, genome1, genome2);
        orthologs.Should().BeEmpty("no cross pair shares any 5-mer, so none clears the 0.3 identity / 0.5 coverage gate → no ortholog (INV-05)");
    }

    /// <summary>
    /// BE: a borderline pair just BELOW the qualification gate must still produce no
    /// ortholog — pinning that INV-05 is a real exclusion, not a near-miss admission.
    /// We raise minIdentity to a value strictly above the only candidate pair's Jaccard
    /// so the otherwise-best, otherwise-reciprocal hit is rejected by the gate alone.
    /// </summary>
    [Test]
    public void FindOrthologs_BestHitBelowThreshold_ExcludedByGate()
    {
        // Two moderately similar genes that WOULD be reciprocal best hits at the default
        // threshold; we push minIdentity above 0.99 so only an exact match could pass.
        string seqA = RandomDna(120, seed: 7001);
        string seqB = seqA[..100] + RandomDna(20, seed: 7002); // ~83% shared prefix → Jaccard < 1
        var genome1 = new[] { Gene("a1", "G1", seqA) };
        var genome2 = new[] { Gene("b1", "G2", seqB) };

        // Sanity: at the default gate they ARE orthologs (the pair is reciprocal & qualifies).
        ComparativeGenomics.FindOrthologs(genome1, genome2).Should().ContainSingle(
            "at the default 0.3 gate the similar pair is a reciprocal best hit");

        // At a near-1.0 identity gate the imperfect pair no longer qualifies.
        var strict = ComparativeGenomics.FindOrthologs(genome1, genome2, minIdentity: 0.999).ToList();
        AssertWellFormedMatching(strict, genome1, genome2);
        strict.Should().BeEmpty("the pair's Jaccard is below 0.999, so the qualification gate excludes it (INV-05)");
    }

    /// <summary>
    /// BE: a strictly ONE-DIRECTIONAL best hit must NOT be returned (INV-01, §6.1).
    /// Construct g2's best match to be a1, but a1's best match to be a different,
    /// more-similar gene b2, so a1→b2 (not b1) while b1→a1. The b1↔a1 relation is not
    /// reciprocal and must be dropped; the only reciprocal pair is a1↔b2.
    /// </summary>
    [Test]
    public void FindOrthologs_OneDirectionalBestHit_NotReturned()
    {
        string anchor = RandomDna(150, seed: 8100);
        // b2 is nearly identical to anchor (a1) → mutual best.
        string b2seq = anchor[..140] + RandomDna(10, seed: 8101);
        // b1 is somewhat similar to anchor but a1's best is b2; b1's best is a1 (no closer gene).
        string b1seq = anchor[..90] + RandomDna(60, seed: 8102);

        var genome1 = new[] { Gene("a1", "G1", anchor) };
        var genome2 = new[]
        {
            Gene("b1", "G2", b1seq),
            Gene("b2", "G2", b2seq),
        };

        var orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2).ToList();

        AssertWellFormedMatching(orthologs, genome1, genome2);
        orthologs.Should().ContainSingle("only the mutually-best a1↔b2 pair is reciprocal");
        orthologs[0].Gene1Id.Should().Be("a1");
        orthologs[0].Gene2Id.Should().Be("b2", "a1's best hit is b2, not the one-directional b1 (INV-01)");
        AssertReciprocal(genome1, genome2);
    }

    #endregion

    #region BE — Boundary: all-identical (deterministic matching, no quadratic hang)

    /// <summary>
    /// BE: every gene identical across both genomes is the degenerate ceiling. With N
    /// identical genes per genome, each gene's best hit is a perfect match (Jaccard 1.0,
    /// coverage 1.0) but the deterministic tie-break (smaller ordinal id, §4.2) makes
    /// EVERY gene's best hit the smallest-ordinal-id gene on the other side. Reciprocity
    /// therefore holds for the two smallest-id genes only → EXACTLY ONE pair (INV-02:
    /// output is a matching, not N² pairs). We pin that single pair, that it is perfect,
    /// and — via [CancelAfter] — that the all-identical case does not hang on a quadratic
    /// blow-up.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindOrthologs_AllIdentical_IsSingleDeterministicPairNoHang()
    {
        string shared = RandomDna(100, seed: 9001);
        var genome1 = Enumerable.Range(0, 30)
            .Select(i => Gene($"a{i:00}", "G1", shared)).ToArray();
        var genome2 = Enumerable.Range(0, 30)
            .Select(i => Gene($"b{i:00}", "G2", shared)).ToArray();

        var orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2).ToList();

        AssertWellFormedMatching(orthologs, genome1, genome2);
        orthologs.Should().ContainSingle(
            "all genes are identical so every best hit collapses (by ordinal tie-break) onto the smallest-id gene; only that one pair is reciprocal → a single pair, not 30 (INV-02)");

        var only = orthologs[0];
        only.Gene1Id.Should().Be("a00", "the smallest-ordinal genome-1 id wins the tie-break");
        only.Gene2Id.Should().Be("b00", "the smallest-ordinal genome-2 id wins the tie-break");
        only.Identity.Should().Be(1.0, "identical sequences have Jaccard 1.0");
        only.Coverage.Should().Be(1.0, "identical sequences share all of the shorter sequence's 5-mers → coverage 1.0");
    }

    /// <summary>
    /// BE: all-identical determinism (INV-04). Shuffling the gene order on both sides
    /// must yield the identical matching — the deterministic best-hit ranking with the
    /// ordinal tie-break makes the output order-independent. This pins that the
    /// all-identical collapse is stable, not order-dependent.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindOrthologs_AllIdentical_IsOrderIndependent()
    {
        string shared = RandomDna(80, seed: 9100);
        var genes1 = Enumerable.Range(0, 12).Select(i => Gene($"a{i:00}", "G1", shared)).ToList();
        var genes2 = Enumerable.Range(0, 12).Select(i => Gene($"b{i:00}", "G2", shared)).ToList();

        var ordered = ComparativeGenomics.FindOrthologs(genes1, genes2)
            .Select(p => (p.Gene1Id, p.Gene2Id)).ToList();

        var rng = new Random(424242);
        var shuffled1 = genes1.OrderBy(_ => rng.Next()).ToList();
        var shuffled2 = genes2.OrderBy(_ => rng.Next()).ToList();
        var fromShuffled = ComparativeGenomics.FindOrthologs(shuffled1, shuffled2)
            .Select(p => (p.Gene1Id, p.Gene2Id)).ToList();

        fromShuffled.Should().BeEquivalentTo(ordered,
            "the deterministic argmax + ordinal tie-break makes the matching independent of input order (INV-04)");
        ordered.Should().BeEquivalentTo(new[] { ("a00", "b00") },
            "the all-identical collapse picks the smallest-ordinal pair regardless of order");
    }

    /// <summary>
    /// BE: a larger all-identical input must still complete promptly (no quadratic hang
    /// in practice for a modest N) and still collapse to one pair. The [CancelAfter]
    /// budget is the hang guard; the single-pair assertion is the correctness guard.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindOrthologs_ManyIdenticalGenes_TerminatesWithSinglePair()
    {
        string shared = RandomDna(60, seed: 9200);
        var genome1 = Enumerable.Range(0, 80).Select(i => Gene($"a{i:000}", "G1", shared)).ToArray();
        var genome2 = Enumerable.Range(0, 80).Select(i => Gene($"b{i:000}", "G2", shared)).ToArray();

        List<ComparativeGenomics.OrthologPair> orthologs = null!;
        FluentActions.Invoking(() => orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2).ToList())
            .Should().NotThrow("a modest all-identical input must complete, not hang or overflow");

        AssertWellFormedMatching(orthologs, genome1, genome2);
        orthologs.Should().ContainSingle("the all-identical collapse always yields exactly one reciprocal pair (INV-02)");
    }

    #endregion

    #region BE — Boundary: empty gene set (zero-side safety, no NRE / DivideByZero)

    /// <summary>
    /// BE: an empty gene set on EITHER side is the degenerate floor of the size axis. A
    /// pair needs one gene from each genome, so an empty side → no orthologs (§6.1). The
    /// `Count == 0` short-circuit (ComparativeGenomics.cs lines 477-478) returns an
    /// empty result before any best-hit search, so there is no NullReference and no
    /// DivideByZero on an empty k-mer set. We probe empty-left, empty-right and
    /// empty-both, all of which must be empty and never throw.
    /// </summary>
    [Test]
    public void FindOrthologs_EmptyGeneSet_ReturnsEmptyNoThrow()
    {
        var real = new[] { Gene("a1", "G1", RandomDna(50, seed: 1500)) };
        var empty = Array.Empty<ComparativeGenomics.Gene>();

        FluentActions.Invoking(() => ComparativeGenomics.FindOrthologs(empty, real).ToList())
            .Should().NotThrow("an empty genome-1 short-circuits to an empty result, never an NRE");
        ComparativeGenomics.FindOrthologs(empty, real).Should().BeEmpty("no genome-1 gene → no pair (§6.1)");

        ComparativeGenomics.FindOrthologs(real, empty).Should().BeEmpty("no genome-2 gene → no pair (§6.1)");
        ComparativeGenomics.FindOrthologs(empty, empty).Should().BeEmpty("both empty → empty result, never a 0/0 (§6.1)");
    }

    /// <summary>
    /// BE: a single gene on each side (the smallest non-empty input) that are clear
    /// mutual best hits must yield exactly one ortholog pair, with no over-indexing on
    /// the single-element best-hit search. The other half of the boundary — a single
    /// gene on each side with NO similarity — must yield no pair.
    /// </summary>
    [Test]
    public void FindOrthologs_SingleGeneEachSide_IsWellFormed()
    {
        string shared = RandomDna(90, seed: 1600);
        var one1 = new[] { Gene("a1", "G1", shared) };
        var one2 = new[] { Gene("b1", "G2", shared) };

        var paired = ComparativeGenomics.FindOrthologs(one1, one2).ToList();
        AssertWellFormedMatching(paired, one1, one2);
        paired.Should().ContainSingle("two identical lone genes are mutual best hits → one ortholog pair");
        paired[0].Gene1Id.Should().Be("a1");
        paired[0].Gene2Id.Should().Be("b1");

        var noMatch1 = new[] { Gene("a1", "G1", new string('A', 60)) };
        var noMatch2 = new[] { Gene("b1", "G2", new string('C', 60)) };
        ComparativeGenomics.FindOrthologs(noMatch1, noMatch2)
            .Should().BeEmpty("two dissimilar lone genes share no 5-mer → no qualifying hit → no pair");
    }

    /// <summary>
    /// BE: a null gene list is the ONE input the contract rejects with an exception
    /// (§3.3, §6.1; ComparativeGenomics.cs lines 471-472 ThrowIfNull). Because the
    /// implementation is an iterator-delegating method, the ThrowIfNull runs eagerly on
    /// the public call, so a null on either side throws ArgumentNullException up-front
    /// (not a deferred NRE on enumeration). We pin the throw for null-left and null-right.
    /// </summary>
    [Test]
    public void FindOrthologs_NullGeneList_ThrowsArgumentNull()
    {
        var real = new[] { Gene("a1", "G1", RandomDna(40, seed: 1700)) };

        FluentActions.Invoking(() => ComparativeGenomics.FindOrthologs(null!, real))
            .Should().Throw<ArgumentNullException>("a null genome-1 list is rejected by ThrowIfNull, not a deferred NRE (§3.3)");
        FluentActions.Invoking(() => ComparativeGenomics.FindOrthologs(real, null!))
            .Should().Throw<ArgumentNullException>("a null genome-2 list is rejected by ThrowIfNull (§3.3)");
    }

    #endregion

    #region MC — Malformed Content (null/empty/sub-k/junk gene sequences)

    /// <summary>
    /// MC: genes whose Sequence is null or empty are MALFORMED payload — they carry no
    /// comparable sequence. The documented behaviour (§3.3, §6.1; the
    /// `Where(g => !string.IsNullOrEmpty(g.Sequence))` filter, lines 474-475) is to SKIP
    /// them, not to crash. A genome that is ALL malformed genes therefore behaves like an
    /// empty genome → no orthologs, no NullReference on the null Sequence. We mix
    /// sequenced and malformed genes and confirm only the real pair survives.
    /// </summary>
    [Test]
    public void FindOrthologs_NullOrEmptySequences_AreSkippedNoCrash()
    {
        string shared = RandomDna(70, seed: 2100);
        var genome1 = new[]
        {
            Gene("a_null", "G1", null),
            Gene("a_empty", "G1", ""),
            Gene("a_real", "G1", shared),
        };
        var genome2 = new[]
        {
            Gene("b_real", "G2", shared),
            Gene("b_null", "G2", null),
            Gene("b_empty", "G2", ""),
        };

        List<ComparativeGenomics.OrthologPair> orthologs = null!;
        FluentActions.Invoking(() => orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2).ToList())
            .Should().NotThrow("null/empty gene sequences must be skipped, never dereferenced (§3.3)");

        AssertWellFormedMatching(orthologs, genome1, genome2);
        orthologs.Should().ContainSingle("only the two real-sequence genes can pair; the malformed ones are skipped");
        orthologs[0].Gene1Id.Should().Be("a_real");
        orthologs[0].Gene2Id.Should().Be("b_real");
    }

    /// <summary>
    /// MC: every gene malformed (all null/empty sequences) on both sides must behave like
    /// the empty-gene-set boundary — an empty result, never a DivideByZero or NRE. This
    /// pins that the skip-filter + Count==0 short-circuit compose correctly for the
    /// all-malformed degenerate case.
    /// </summary>
    [Test]
    public void FindOrthologs_AllMalformedSequences_ReturnsEmptyNoThrow()
    {
        var genome1 = new[] { Gene("a1", "G1", null), Gene("a2", "G1", "") };
        var genome2 = new[] { Gene("b1", "G2", ""), Gene("b2", "G2", null) };

        List<ComparativeGenomics.OrthologPair> orthologs = null!;
        FluentActions.Invoking(() => orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2).ToList())
            .Should().NotThrow("an all-malformed input is filtered to empty, then short-circuits to an empty result");

        orthologs.Should().BeEmpty("after dropping all null/empty-sequence genes there is nothing to pair");
    }

    /// <summary>
    /// MC: a sequence SHORTER than the 5-mer length contributes no k-mers (the
    /// `seq.Length &lt; k` guard returns (0,0,0), lines 524-525). Such a gene therefore
    /// never qualifies (identity 0 &lt; 0.3) and is never falsely paired, and the empty
    /// k-mer sets do NOT DivideByZero (the total&gt;0 / minKmerCount&gt;0 guards). We pin
    /// that a sub-k gene produces no ortholog even against an identical sub-k gene, and
    /// that the call does not throw.
    /// </summary>
    [Test]
    public void FindOrthologs_SubKmerSequences_NeverQualifyNoDivideByZero()
    {
        // Length-3 sequences (< k = 5) on both sides, even though identical.
        var genome1 = new[] { Gene("a1", "G1", "ACG") };
        var genome2 = new[] { Gene("b1", "G2", "ACG") };

        List<ComparativeGenomics.OrthologPair> orthologs = null!;
        FluentActions.Invoking(() => orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2).ToList())
            .Should().NotThrow("a sub-k sequence yields an empty k-mer set; the total>0/minKmerCount>0 guards prevent 0/0");

        AssertWellFormedMatching(orthologs, genome1, genome2);
        orthologs.Should().BeEmpty("a sequence shorter than k=5 contributes no k-mers → identity 0 < 0.3 → never qualifies (§3.3, §6.2)");
    }

    /// <summary>
    /// MC: junk / non-AA / non-DNA characters in a gene sequence must not crash the
    /// k-mer machinery — the similarity is a pure string-k-mer Jaccard over whatever
    /// characters are present (upper-cased, §3.3). Two genes carrying the SAME junk are
    /// still mutual best hits (a perfectly valid string match); two genes carrying
    /// DIFFERENT junk share no k-mer and do not pair. Either way: no exception, a
    /// well-formed matching.
    /// </summary>
    [Test]
    public void FindOrthologs_JunkCharactersInSequence_NoCrashWellFormed()
    {
        string junk = "!!@@##$$%%^^&&**(())___++==[[]]{{}}<<>>??//\\\\||~~``";
        var genome1 = new[] { Gene("a1", "G1", junk) };
        var genome2Same = new[] { Gene("b1", "G2", junk) };

        List<ComparativeGenomics.OrthologPair> same = null!;
        FluentActions.Invoking(() => same = ComparativeGenomics.FindOrthologs(genome1, genome2Same).ToList())
            .Should().NotThrow("junk characters are just k-mer symbols; the Jaccard machinery must not crash on them");
        AssertWellFormedMatching(same, genome1, genome2Same);
        same.Should().ContainSingle("identical junk strings are still mutual best hits (Jaccard 1.0)");

        var genome2Diff = new[] { Gene("b1", "G2", "@@##$$%%^^abcdefghij@@##$$%%^^") };
        var diffGenome1 = new[] { Gene("a1", "G1", "ZZZZZZZZZZZZZZZZZZZZ") };
        List<ComparativeGenomics.OrthologPair> diff = null!;
        FluentActions.Invoking(() => diff = ComparativeGenomics.FindOrthologs(diffGenome1, genome2Diff).ToList())
            .Should().NotThrow("disjoint junk strings must not crash the k-mer comparison");
        diff.Should().BeEmpty("two genes with no shared k-mer (even of junk) do not pair");
    }

    /// <summary>
    /// MC fuzz sweep — across many random genomes seeded with malformed genes (null,
    /// empty, sub-k, junk) interleaved with real sequences, FindOrthologs must ALWAYS
    /// return a well-formed matching and never throw. This is the broad
    /// undisciplined-failure net behind the targeted MC tests.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void FindOrthologs_RandomMalformedMix_AlwaysWellFormedNoThrow()
    {
        for (int seed = 0; seed < 30; seed++)
        {
            var rng = new Random(seed);
            var g1 = new List<ComparativeGenomics.Gene>();
            var g2 = new List<ComparativeGenomics.Gene>();

            int n1 = rng.Next(0, 8);
            int n2 = rng.Next(0, 8);
            for (int i = 0; i < n1; i++) g1.Add(Gene($"a{i}", "G1", RandomSequenceOrJunk(rng, seed * 31 + i)));
            for (int i = 0; i < n2; i++) g2.Add(Gene($"b{i}", "G2", RandomSequenceOrJunk(rng, seed * 53 + i)));

            List<ComparativeGenomics.OrthologPair> orthologs = null!;
            FluentActions.Invoking(() => orthologs = ComparativeGenomics.FindOrthologs(g1, g2).ToList())
                .Should().NotThrow($"a random malformed mix must never crash FindOrthologs (seed {seed})");

            AssertWellFormedMatching(orthologs, g1, g2);
        }
    }

    private static string? RandomSequenceOrJunk(Random rng, int seed)
    {
        return rng.Next(5) switch
        {
            0 => null,
            1 => "",
            2 => RandomDna(rng.Next(1, 4), seed),       // sub-k
            3 => "##" + new string('@', rng.Next(0, 30)), // junk
            _ => RandomDna(rng.Next(20, 120), seed),       // real
        };
    }

    #endregion

    #region Positive sanity — the documented RBH ortholog on hand-built examples

    /// <summary>
    /// Positive sanity (Ortholog_Identification.md §7.1 worked example): two identical
    /// genes across genomes are mutual best hits → exactly one ortholog pair, with the
    /// documented Jaccard 1.0 and coverage 1.0. This pins the load-bearing metric, not
    /// merely a green pass.
    /// </summary>
    [Test]
    public void FindOrthologs_WorkedExample_PinsDocumentedPairAndScore()
    {
        var g1 = new[] { Gene("a1", "G1", "ACGTACGTACGTACGT") };
        var g2 = new[] { Gene("b1", "G2", "ACGTACGTACGTACGT") };

        var orthologs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();

        AssertWellFormedMatching(orthologs, g1, g2);
        orthologs.Should().ContainSingle("two identical genes are reciprocal best hits → one pair (§7.1)");
        orthologs[0].Gene1Id.Should().Be("a1");
        orthologs[0].Gene2Id.Should().Be("b1");
        orthologs[0].Identity.Should().Be(1.0, "identical sequences → 5-mer Jaccard 1.0 (§7.1)");
        orthologs[0].Coverage.Should().Be(1.0, "identical sequences share all 5-mers of the shorter sequence → coverage 1.0");
        AssertReciprocal(g1, g2);
    }

    /// <summary>
    /// Positive sanity — the discriminating case: a clearly orthologous (highly similar)
    /// pair IS reported, while a divergent/unrelated gene in the same genome is NOT
    /// paired. genome1 has the true ortholog (a1↔b1) plus an unrelated gene (a2, poly-G);
    /// genome2 has b1 plus an unrelated gene (b2, poly-T). Only the similar a1↔b1 pair
    /// must appear; the unrelated genes share no 5-mer with anything and stay unpaired.
    /// This is the core business meaning of ortholog detection.
    /// </summary>
    [Test]
    public void FindOrthologs_RelatedPairedUnrelatedUnpaired()
    {
        string orthoSeq = RandomDna(140, seed: 30100);
        // b1 is a1 with a few scattered substitutions → still mutual best, comfortably > gate.
        char[] m = orthoSeq.ToCharArray();
        var rng = new Random(30101);
        for (int i = 0; i < 4; i++)
        {
            int pos = rng.Next(m.Length);
            m[pos] = m[pos] == 'A' ? 'C' : 'A';
        }
        string b1Seq = new string(m);

        var genome1 = new[]
        {
            Gene("a1", "G1", orthoSeq),
            Gene("a2", "G1", new string('G', 100)), // unrelated
        };
        var genome2 = new[]
        {
            Gene("b1", "G2", b1Seq),
            Gene("b2", "G2", new string('T', 100)), // unrelated
        };

        var orthologs = ComparativeGenomics.FindOrthologs(genome1, genome2).ToList();

        AssertWellFormedMatching(orthologs, genome1, genome2);
        orthologs.Should().ContainSingle("only the truly similar a1↔b1 pair is a reciprocal best hit");
        orthologs[0].Gene1Id.Should().Be("a1");
        orthologs[0].Gene2Id.Should().Be("b1", "the divergent/unrelated genes (a2, b2) share no 5-mer and stay unpaired");
        orthologs[0].Identity.Should().BeGreaterThan(0.3, "a few substitutions keep the pair well above the qualification gate");
        AssertReciprocal(genome1, genome2);
    }

    /// <summary>
    /// Positive sanity — case-insensitivity (§3.3). A lowercase gene sequence is
    /// uppercased for k-mer extraction, so a lowercase gene is identical to its uppercase
    /// twin across genomes → one ortholog pair, Jaccard 1.0. Guards the documented
    /// ToUpperInvariant normalization from treating case as divergence.
    /// </summary>
    [Test]
    public void FindOrthologs_LowercaseSequence_IsUppercasedBeforeComparison()
    {
        string seq = RandomDna(120, seed: 30200);
        var g1 = new[] { Gene("a1", "G1", seq.ToLowerInvariant()) };
        var g2 = new[] { Gene("b1", "G2", seq) };

        var orthologs = ComparativeGenomics.FindOrthologs(g1, g2).ToList();

        AssertWellFormedMatching(orthologs, g1, g2);
        orthologs.Should().ContainSingle("case is normalized via ToUpperInvariant, so the lowercase gene is identical to its uppercase twin (§3.3)");
        orthologs[0].Identity.Should().Be(1.0, "after upper-casing the sequences are identical → Jaccard 1.0");
    }

    #endregion

    #endregion
}
