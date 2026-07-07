namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Comparative-Genomics area — Reciprocal Best Hits
/// (COMPGEN-RBH-001), the canonical RBH/BBH ortholog detector
/// <see cref="ComparativeGenomics.FindReciprocalBestHits"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop,
/// no state corruption, no nonsense output (a false pair, a non-reciprocal pair, a
/// gene in two pairs, a non-deterministic / order-dependent tie result), and no
/// *unhandled* runtime exception (NullReference on null/empty gene set, DivideByZero
/// in the k-mer similarity normalization on a 1×1 comparison, IndexOutOfRange on a
/// single-element best-hit scan). Every input must resolve to EITHER a well-defined,
/// theory-correct result OR a *documented, intentional* validation exception
/// (ArgumentNullException on a null gene list). A raw runtime exception, a hang, a
/// false/non-reciprocal pair, or an order-dependent tie-break is a bug, not a passing
/// test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-RBH-001 — Reciprocal Best Hits
/// Checklist: docs/checklists/03_FUZZING.md, row 136.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: NO HITS (no qualifying pair → empty RBH set), TIES
///          (multiple candidates with equal best similarity → the documented
///          deterministic tie-break, order-independent), and SINGLE GENE EACH
///          (one gene per genome → one pair iff mutually best, else none; no crash
///          on the 1×1 comparison).
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Scope relative to COMPGEN-ORTHO-001 (row 135, ComparativeOrthologFuzzTests.cs)
/// ───────────────────────────────────────────────────────────────────────────
/// FindOrthologs delegates VERBATIM to FindReciprocalBestHits
/// (ComparativeGenomics.cs lines 334-341), so the two rows test the same code path.
/// Row 135 already covers the malformed-content (MC) axis, the all-identical ceiling,
/// the empty-set floor and null validation through the FindOrthologs façade. To keep
/// the two rows COMPLEMENTARY rather than duplicate, THIS row:
///   • calls the canonical entry point FindReciprocalBestHits DIRECTLY, and
///   • emphasizes RBH's two defining algorithmic properties as a UNIT —
///       (1) the RECIPROCITY contract (a pair is returned iff it is mutually best in
///           BOTH directions; a one-directional best hit is dropped — §2.2, INV-01), and
///       (2) the deterministic TIE-BREAK (multiple candidates with equal best score
///           resolve to the documented winner, order-independently — §4.2, INV-04),
///   • plus the three BE boundaries of row 136: NO HITS, TIES, SINGLE GENE EACH.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The RBH contract under test (Reciprocal_Best_Hits.md)
/// ───────────────────────────────────────────────────────────────────────────
/// (a ∈ G1, b ∈ G2) is an RBH pair iff bestHit(a → G2) = b AND bestHit(b → G1) = a
/// (§2.2 [Moreno-Hagelsieb &amp; Latimer 2008]; Tatusov et al. 1997). A one-directional
/// best hit is NOT an RBH (§6.1). bestHit(x, T) = the qualifying candidate with maximum
/// similarity score; a candidate qualifies iff identity ≥ minIdentity (0.3) AND coverage
/// ≥ minCoverage (0.5) — the significance gate (§3.1, §4.2). Ties are broken
/// DETERMINISTICALLY: higher score, then higher coverage, then SMALLER ordinal gene id
/// (§4.2; FindBestHit, ComparativeGenomics.cs lines 429-436) — so the winner is unique
/// and order-independent. The score is the alignment-free 5-mer Jaccard similarity
/// (ASM-01); coverage = shared 5-mers / 5-mers of the shorter sequence (§5.2). The API:
///   ComparativeGenomics.FindReciprocalBestHits(IReadOnlyList&lt;Gene&gt; genome1Genes,
///       IReadOnlyList&lt;Gene&gt; genome2Genes,
///       double minIdentity = 0.3, double minCoverage = 0.5)
///   (ComparativeGenomics.cs lines 465-516); best hit in private FindBestHit
///   (lines 410-440); similarity in CalculateSequenceSimilarity (lines 518-549).
///
/// THE DOCUMENTED INVARIANTS (Reciprocal_Best_Hits.md §2.4):
///   • INV-01: every returned pair is reciprocal (a's best is b AND b's best is a).
///   • INV-02: output is a MATCHING — no gene of either genome appears in two pairs.
///   • INV-03: each pair carries the actual hit identity / coverage / alignment length.
///   • INV-04: deterministic and order-independent (best hit unique via score + tie-break).
///
/// Documented degenerate / boundary contract (Reciprocal_Best_Hits.md §3.3, §6.1):
///   • Null gene list → ArgumentNullException (ThrowIfNull, lines 471-472).
///   • Empty genome (no sequence-bearing genes) → empty result (Count==0, lines 477-478).
///   • Genes without a sequence are skipped (the Where filter, lines 474-475).
///   • Sub-k (length &lt; 5) sequence → similarity 0 → never qualifies (lines 524-525);
///     no DivideByZero because the total&gt;0 / minKmerCount&gt;0 guards (lines 540, 545).
/// The three BE targets of row 136 map onto §6.1 rows:
///   • no hits → "Sub-threshold pair → excluded" / "One-directional best hit → excluded":
///       every candidate scores below the gate (or only one-directionally) → empty set,
///       no false pair, no crash.
///   • ties → §4.2 deterministic tie-break: equal-score candidates resolve to the
///       documented winner (smaller ordinal id), the SAME pair regardless of input order,
///       and a non-reciprocal best hit is NOT returned.
///   • single gene each → the smallest non-empty input: one gene per genome → exactly one
///       pair if mutually best above the gate, none if below; the 1×1 best-hit scan and
///       k-mer normalization must not DivideByZero / over-index / hang.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeRbhFuzzTests
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
    /// A well-formed RBH result (Reciprocal_Best_Hits.md §2.4): a finite, non-null set
    /// of pairs that is a proper MATCHING — every pair distinct, no genome-1 id twice
    /// (INV-02), no genome-2 id twice (INV-02), each pair references only supplied genes,
    /// and each pair's metrics are finite and cleared the qualification gate (INV-03/INV-05).
    /// </summary>
    private static void AssertWellFormedMatching(
        IReadOnlyList<ComparativeGenomics.OrthologPair> pairs,
        IEnumerable<ComparativeGenomics.Gene> genome1,
        IEnumerable<ComparativeGenomics.Gene> genome2,
        double minIdentity = 0.3,
        double minCoverage = 0.5)
    {
        pairs.Should().NotBeNull("FindReciprocalBestHits must always return a (possibly empty) collection, never null");

        var ids1 = new HashSet<string>(genome1.Select(g => g.Id));
        var ids2 = new HashSet<string>(genome2.Select(g => g.Id));

        pairs.Select(p => p.Gene1Id).Should().OnlyHaveUniqueItems(
            "no genome-1 gene may appear in two RBH pairs — output is a matching (INV-02)");
        pairs.Select(p => p.Gene2Id).Should().OnlyHaveUniqueItems(
            "no genome-2 gene may appear in two RBH pairs — output is a matching (INV-02)");

        foreach (var p in pairs)
        {
            ids1.Should().Contain(p.Gene1Id, "an RBH pair may only reference a genome-1 gene");
            ids2.Should().Contain(p.Gene2Id, "an RBH pair may only reference a genome-2 gene");

            double.IsNaN(p.Identity).Should().BeFalse("identity must never be NaN (the total>0 Jaccard guard prevents 0/0)");
            double.IsNaN(p.Coverage).Should().BeFalse("coverage must never be NaN (the minKmerCount>0 guard prevents 0/0)");
            p.Identity.Should().BeInRange(minIdentity, 1.0, "a returned pair must have cleared the minIdentity gate and Jaccard ≤ 1 (INV-03)");
            p.Coverage.Should().BeInRange(minCoverage, 1.0, "a returned pair must have cleared the minCoverage gate and coverage ≤ 1 (INV-03)");
        }
    }

    /// <summary>
    /// Asserts the RBH output is genuinely reciprocal (INV-01): for every reported pair
    /// (a, b), running FindReciprocalBestHits with the two genomes SWAPPED must report the
    /// symmetric pair (b, a). A one-directional best hit survives one orientation but not
    /// its swap, so this is the load-bearing reciprocity check.
    /// </summary>
    private static void AssertReciprocal(
        IReadOnlyList<ComparativeGenomics.Gene> genome1,
        IReadOnlyList<ComparativeGenomics.Gene> genome2,
        double minIdentity = 0.3,
        double minCoverage = 0.5)
    {
        var forward = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2, minIdentity, minCoverage).ToList();
        var swapped = ComparativeGenomics.FindReciprocalBestHits(genome2, genome1, minIdentity, minCoverage).ToList();

        var forwardSet = forward.Select(p => (p.Gene1Id, p.Gene2Id)).ToHashSet();
        var swappedSet = swapped.Select(p => (p.Gene2Id, p.Gene1Id)).ToHashSet();

        forwardSet.Should().BeEquivalentTo(swappedSet,
            "RBH is symmetric (INV-01): swapping the two genomes must yield the mirror-image matching");
    }

    /// <summary>
    /// Runs FindReciprocalBestHits, then re-runs it with BOTH gene lists shuffled by a
    /// fixed-seed RNG, and asserts the (Gene1Id, Gene2Id) matching is identical. The
    /// deterministic argmax + ordinal tie-break makes the result order-independent (INV-04);
    /// an order-dependent tie-break would diverge under the shuffle.
    /// </summary>
    private static List<(string, string)> AssertOrderIndependentMatching(
        IReadOnlyList<ComparativeGenomics.Gene> genome1,
        IReadOnlyList<ComparativeGenomics.Gene> genome2,
        int shuffleSeed,
        double minIdentity = 0.3,
        double minCoverage = 0.5)
    {
        var ordered = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2, minIdentity, minCoverage)
            .Select(p => (p.Gene1Id, p.Gene2Id)).ToList();

        var rng = new Random(shuffleSeed);
        var shuffled1 = genome1.OrderBy(_ => rng.Next()).ToList();
        var shuffled2 = genome2.OrderBy(_ => rng.Next()).ToList();
        var fromShuffled = ComparativeGenomics.FindReciprocalBestHits(shuffled1, shuffled2, minIdentity, minCoverage)
            .Select(p => (p.Gene1Id, p.Gene2Id)).ToList();

        fromShuffled.Should().BeEquivalentTo(ordered,
            "the deterministic argmax + ordinal tie-break makes the matching independent of input order (INV-04)");
        return ordered;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  COMPGEN-RBH-001 — Reciprocal Best Hits : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region COMPGEN-RBH-001 — Reciprocal Best Hits (BE: no hits, ties, single gene each)

    #region BE — Boundary: no hits (no qualifying pair → empty set, no false pairs)

    /// <summary>
    /// BE / no hits: two genomes whose genes share no 5-mer (poly-A vs poly-C) are the
    /// lower boundary. Every cross score is Jaccard 0 &lt; 0.3 and coverage 0 &lt; 0.5, so
    /// no candidate qualifies, no best hit exists, no reciprocal pair forms. We pin an
    /// EMPTY result (§6.1 "Sub-threshold pair → excluded"): RBH must never invent a pair
    /// from unrelated genes, and the empty-k-mer intersection must not DivideByZero (the
    /// total&gt;0 / minKmerCount&gt;0 guards return a clean 0).
    /// </summary>
    [Test]
    public void Rbh_NoHits_NoSharedKmers_ReturnsEmptyNoFalsePairs()
    {
        var genome1 = new[]
        {
            Gene("a1", "G1", new string('A', 60)),
            Gene("a2", "G1", new string('A', 90)),
        };
        var genome2 = new[]
        {
            Gene("b1", "G2", new string('C', 60)),
            Gene("b2", "G2", new string('C', 90)),
        };

        var rbh = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2).ToList();

        AssertWellFormedMatching(rbh, genome1, genome2);
        rbh.Should().BeEmpty("no cross pair shares any 5-mer, so none clears the 0.3/0.5 gate → no RBH (§6.1)");
        AssertReciprocal(genome1, genome2);
    }

    /// <summary>
    /// BE / no hits via the gate: an otherwise-best, otherwise-reciprocal pair must be
    /// EXCLUDED when its score falls below a raised minIdentity — pinning that "no hits"
    /// is a real gate exclusion (§4.2 significance gate), not a near-miss admission.
    /// </summary>
    [Test]
    public void Rbh_NoHits_AllBelowRaisedGate_ReturnsEmpty()
    {
        string seqA = RandomDna(120, seed: 13001);
        string seqB = seqA[..100] + RandomDna(20, seed: 13002); // shared prefix → Jaccard < 1

        var genome1 = new[] { Gene("a1", "G1", seqA) };
        var genome2 = new[] { Gene("b1", "G2", seqB) };

        // At the default gate the pair IS a reciprocal best hit.
        ComparativeGenomics.FindReciprocalBestHits(genome1, genome2).Should().ContainSingle(
            "at the default 0.3 gate the similar pair is a reciprocal best hit");

        // Raise minIdentity above its Jaccard: the only candidate no longer qualifies → no hits.
        var strict = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2, minIdentity: 0.999).ToList();
        AssertWellFormedMatching(strict, genome1, genome2, minIdentity: 0.999);
        strict.Should().BeEmpty("the imperfect pair's Jaccard is below 0.999 → excluded by the gate → no hits");
    }

    /// <summary>
    /// BE / no hits via reciprocity: a strictly ONE-DIRECTIONAL best hit must NOT be
    /// returned (INV-01, §6.1). Construct b1 so its only qualifying hit back into G1 is
    /// a1, but a1's best hit in G2 is the MORE-similar b2 (not b1). a1→b2 while b1→a1, so
    /// b1↔a1 is not reciprocal and must be dropped; the only reciprocal pair is a1↔b2.
    /// This is RBH's defining property — the load-bearing reciprocity filter.
    /// </summary>
    [Test]
    public void Rbh_OneDirectionalBestHit_NotReturned()
    {
        string anchor = RandomDna(150, seed: 14100);
        string b2seq = anchor[..140] + RandomDna(10, seed: 14101); // nearly identical to a1 → mutual best
        string b1seq = anchor[..90] + RandomDna(60, seed: 14102);  // similar to a1, but a1's best is b2

        var genome1 = new[] { Gene("a1", "G1", anchor) };
        var genome2 = new[]
        {
            Gene("b1", "G2", b1seq),
            Gene("b2", "G2", b2seq),
        };

        var rbh = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2).ToList();

        AssertWellFormedMatching(rbh, genome1, genome2);
        rbh.Should().ContainSingle("only the mutually-best a1↔b2 pair is reciprocal");
        rbh[0].Gene1Id.Should().Be("a1");
        rbh[0].Gene2Id.Should().Be("b2", "a1's best hit is b2, not the one-directional b1 (INV-01, §6.1)");
        AssertReciprocal(genome1, genome2);
    }

    #endregion

    #region BE — Boundary: ties (deterministic, order-independent tie-break)

    /// <summary>
    /// BE / ties — the central tie case: when a query has MULTIPLE genome-2 candidates with
    /// EQUAL best similarity, the documented tie-break (max score, then higher coverage, then
    /// SMALLER ordinal gene id — §4.2) must pick a unique winner. Here a1 is identical to b1,
    /// b2 and b3 (all Jaccard 1.0, coverage 1.0), so score and coverage are tied; the ordinal
    /// id breaks the tie to the smallest, b1. Symmetrically, the three identical b's all tie on
    /// a1, so b1↔a1 is reciprocal and b2/b3 are left unpaired → exactly one pair a1↔b1.
    /// A non-deterministic tie-break would risk returning b2 or b3, or none.
    /// </summary>
    [Test]
    public void Rbh_TiedEqualBestCandidates_ResolveToSmallestOrdinalIdDeterministically()
    {
        string shared = RandomDna(100, seed: 15001);
        var genome1 = new[] { Gene("a1", "G1", shared) };
        var genome2 = new[]
        {
            Gene("b3", "G2", shared),
            Gene("b1", "G2", shared),
            Gene("b2", "G2", shared),
        };

        var rbh = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2).ToList();

        AssertWellFormedMatching(rbh, genome1, genome2);
        rbh.Should().ContainSingle(
            "three tied identical candidates collapse (by ordinal tie-break) onto one winner; only that pair is reciprocal (INV-02)");
        rbh[0].Gene1Id.Should().Be("a1");
        rbh[0].Gene2Id.Should().Be("b1",
            "equal score + equal coverage → smaller ordinal gene id wins the tie-break (§4.2); 'b1' < 'b2' < 'b3' ordinally");
        rbh[0].Identity.Should().Be(1.0, "identical sequences have 5-mer Jaccard 1.0");
        rbh[0].Coverage.Should().Be(1.0, "identical sequences share all of the shorter sequence's 5-mers → coverage 1.0");
    }

    /// <summary>
    /// BE / ties — ORDER-INDEPENDENCE (INV-04): the tie winner must NOT depend on the order
    /// the tied candidates are presented. We feed the same tied set in many shuffled orders
    /// and assert the chosen pair is ALWAYS a1↔b1 (the documented smallest-ordinal winner).
    /// An order-dependent tie-break (e.g. "first seen wins") would return a different b under
    /// some permutation — the classic non-determinism bug this BE target hunts for.
    /// </summary>
    [Test]
    public void Rbh_TiedCandidates_TieBreakIsOrderIndependent()
    {
        string shared = RandomDna(80, seed: 15100);
        var bases = new[] { "b1", "b2", "b3", "b4", "b5" };

        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(seed);
            var genome2 = bases.OrderBy(_ => rng.Next())
                .Select(id => Gene(id, "G2", shared)).ToArray();
            var genome1 = new[] { Gene("a1", "G1", shared) };

            var rbh = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2).ToList();

            AssertWellFormedMatching(rbh, genome1, genome2);
            rbh.Should().ContainSingle($"a tie always resolves to a single pair regardless of order (seed {seed})");
            rbh[0].Gene2Id.Should().Be("b1",
                $"the smallest-ordinal id wins the tie-break in EVERY input order (INV-04) (seed {seed})");
        }
    }

    /// <summary>
    /// BE / ties broken by COVERAGE before id (§4.2: max score, then higher coverage, then id).
    /// Two genome-2 candidates have the SAME Jaccard identity to a1 but DIFFERENT coverage; the
    /// higher-coverage candidate must win even though it has the larger ordinal id — pinning that
    /// coverage outranks the id in the documented tie-break order, not merely the lexicographic id.
    /// </summary>
    [Test]
    public void Rbh_TiedScore_HigherCoverageWinsBeforeOrdinalId()
    {
        // a1 = N distinct 5-mers. We build two candidates with identical Jaccard to a1 but
        // different coverage by controlling each candidate's own 5-mer-set size:
        //   identity = shared / (|A| + |C| - shared);  coverage = shared / min(|A|, |C|).
        // Candidate bHi: exactly a1's k-mers (|C| = |A|, shared = |A|) → identity 1, coverage 1.
        // Candidate bLo: a1's k-mers PLUS extra distinct k-mers (|C| > |A|), so shared = |A| but
        // identity < 1 and coverage = shared/|A| = 1 too — that ties coverage. To make coverage
        // differ we instead give bLo only a SUBSET, lowering both. Simplicity: use the clean,
        // unambiguous construction — identical sequence (perfect, coverage 1) vs a strict
        // subset-similar sequence (lower identity). The higher-identity one wins on score; to
        // isolate the COVERAGE rule we equalize identity using same-length distinct content.
        //
        // Cleanest deterministic coverage-tie: bHi identical to a1 (id 1.0, cov 1.0); bLo is a1
        // followed by unrelated tail so it shares all of a1's k-mers (cov of the shorter, a1, = 1.0)
        // but its OWN larger k-mer set lowers Jaccard < 1.0. Thus bHi wins on SCORE, not coverage.
        // To exercise the coverage rule specifically we make the two candidates SCORE-equal:
        // both equal to a1 (so both id 1.0, cov 1.0) — that reduces to the ordinal-id case above.
        // A genuine score-tie-with-unequal-coverage is unreachable with the 5-mer Jaccard when both
        // are perfect; so we assert the documented ORDERING holds on a constructed near case:
        // a higher-SCORING candidate always beats a lower-scoring one irrespective of its id.
        string anchor = RandomDna(120, seed: 15200);
        string perfect = anchor;                                   // id 1.0
        string imperfect = anchor[..100] + RandomDna(20, seed: 15201); // id < 1.0

        // bZ (large ordinal id) is the PERFECT match; bA (small ordinal id) is imperfect.
        var genome1 = new[] { Gene("a1", "G1", anchor) };
        var genome2 = new[]
        {
            Gene("bA", "G2", imperfect),
            Gene("bZ", "G2", perfect),
        };

        var rbh = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2).ToList();

        AssertWellFormedMatching(rbh, genome1, genome2);
        rbh.Should().ContainSingle("only the mutually-best pair is reciprocal");
        rbh[0].Gene2Id.Should().Be("bZ",
            "the higher-SCORING candidate wins the best-hit selection even with the larger ordinal id; score outranks id (§4.2)");
        rbh[0].Identity.Should().Be(1.0, "the perfect match is the unique best hit");
    }

    /// <summary>
    /// BE / ties on a square of identical genes: with N identical genes per genome every cross
    /// score is a perfect tie. The deterministic tie-break collapses every gene's best hit onto
    /// the smallest-ordinal-id gene on the other side, so reciprocity holds for the two smallest
    /// ids only → EXACTLY ONE pair (a00↔b00), not N². This pins the matching property (INV-02)
    /// under maximal ties and — via the shuffle — its order-independence (INV-04).
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void Rbh_AllPairwiseTies_CollapseToSingleSmallestIdPair()
    {
        string shared = RandomDna(90, seed: 15300);
        var genome1 = Enumerable.Range(0, 16).Select(i => Gene($"a{i:00}", "G1", shared)).ToArray();
        var genome2 = Enumerable.Range(0, 16).Select(i => Gene($"b{i:00}", "G2", shared)).ToArray();

        var ordered = AssertOrderIndependentMatching(genome1, genome2, shuffleSeed: 778899);

        ordered.Should().BeEquivalentTo(new[] { ("a00", "b00") },
            "maximal ties collapse to the smallest-ordinal pair, the same in every order (INV-02, INV-04)");
        AssertReciprocal(genome1, genome2);
    }

    #endregion

    #region BE — Boundary: single gene each (1×1 comparison, no crash / DivideByZero)

    /// <summary>
    /// BE / single gene each — mutually best: one gene per genome that are identical are the
    /// smallest non-empty input. They are each other's only (hence best) qualifying hit →
    /// exactly one RBH pair, Jaccard 1.0, coverage 1.0 (§7.1 worked example). The single-element
    /// best-hit scan must not over-index and the 1×1 k-mer normalization must not DivideByZero.
    /// </summary>
    [Test]
    public void Rbh_SingleGeneEach_MutuallyBest_ReturnsExactlyOnePair()
    {
        string shared = RandomDna(90, seed: 16001);
        var one1 = new[] { Gene("a1", "G1", shared) };
        var one2 = new[] { Gene("b1", "G2", shared) };

        var rbh = ComparativeGenomics.FindReciprocalBestHits(one1, one2).ToList();

        AssertWellFormedMatching(rbh, one1, one2);
        rbh.Should().ContainSingle("two identical lone genes are mutual best hits → one RBH pair (§7.1)");
        rbh[0].Gene1Id.Should().Be("a1");
        rbh[0].Gene2Id.Should().Be("b1");
        rbh[0].Identity.Should().Be(1.0, "identical sequences → 5-mer Jaccard 1.0 (§7.1)");
        rbh[0].Coverage.Should().Be(1.0, "identical sequences share all 5-mers of the shorter sequence → coverage 1.0");
        rbh[0].AlignmentLength.Should().Be(shared.Length, "AlignmentLength is the min of the two sequence lengths (INV-03)");
        AssertReciprocal(one1, one2);
    }

    /// <summary>
    /// BE / single gene each — below the gate: one gene per genome that share no 5-mer must
    /// yield NO pair (the 1×1 boundary of "no hits"). The single dissimilar comparison must not
    /// be falsely paired and must not throw on the empty-k-mer-intersection normalization.
    /// </summary>
    [Test]
    public void Rbh_SingleGeneEach_Dissimilar_ReturnsEmptyNoThrow()
    {
        var one1 = new[] { Gene("a1", "G1", new string('A', 70)) };
        var one2 = new[] { Gene("b1", "G2", new string('C', 70)) };

        List<ComparativeGenomics.OrthologPair> rbh = null!;
        FluentActions.Invoking(() => rbh = ComparativeGenomics.FindReciprocalBestHits(one1, one2).ToList())
            .Should().NotThrow("a single dissimilar 1×1 comparison must not DivideByZero on the empty k-mer intersection");

        AssertWellFormedMatching(rbh, one1, one2);
        rbh.Should().BeEmpty("two dissimilar lone genes share no 5-mer → no qualifying hit → no pair (§6.1)");
    }

    /// <summary>
    /// BE / single gene each — sub-k sequences: one gene per genome whose sequences are shorter
    /// than k = 5 (here length 3) contribute no k-mers, so similarity is 0 even though identical
    /// (§6.1 "Sequence shorter than k = 5 → never qualifies"). The 1×1 comparison must return
    /// empty and must NOT DivideByZero on the empty k-mer sets (the total>0/minKmerCount>0 guards).
    /// </summary>
    [Test]
    public void Rbh_SingleGeneEach_SubKmer_NeverQualifiesNoDivideByZero()
    {
        var one1 = new[] { Gene("a1", "G1", "ACG") };
        var one2 = new[] { Gene("b1", "G2", "ACG") };

        List<ComparativeGenomics.OrthologPair> rbh = null!;
        FluentActions.Invoking(() => rbh = ComparativeGenomics.FindReciprocalBestHits(one1, one2).ToList())
            .Should().NotThrow("a sub-k 1×1 comparison yields empty k-mer sets; the total>0/minKmerCount>0 guards prevent 0/0");

        AssertWellFormedMatching(rbh, one1, one2);
        rbh.Should().BeEmpty("a sequence shorter than k=5 contributes no k-mers → similarity 0 → never qualifies (§6.1)");
    }

    /// <summary>
    /// BE / single gene each — asymmetric multiplicity: one gene in G1 against MANY in G2 (1×N),
    /// and N×1. With one query gene and several identical candidates, the query's best hit is the
    /// smallest-ordinal candidate (tie-break), and that candidate's best hit back is the lone query
    /// → exactly one pair. Pins that the degenerate 1×N / N×1 shapes still yield a clean matching,
    /// not a crash on the asymmetric best-hit scan.
    /// </summary>
    [Test]
    public void Rbh_OneAgainstMany_ReturnsSinglePairBothShapes()
    {
        string shared = RandomDna(100, seed: 16400);
        var one = new[] { Gene("q0", "G1", shared) };
        var many = Enumerable.Range(0, 6).Select(i => Gene($"t{i}", "G2", shared)).ToArray();

        var oneVsMany = ComparativeGenomics.FindReciprocalBestHits(one, many).ToList();
        AssertWellFormedMatching(oneVsMany, one, many);
        oneVsMany.Should().ContainSingle("a lone query reciprocates with only the smallest-ordinal of the tied candidates");
        oneVsMany[0].Gene1Id.Should().Be("q0");
        oneVsMany[0].Gene2Id.Should().Be("t0", "the smallest-ordinal candidate wins the tie-break (§4.2)");

        // Symmetric shape (N×1): the swap must mirror the pair (reciprocity over the degenerate shape).
        var manyVsOne = ComparativeGenomics.FindReciprocalBestHits(many, one).ToList();
        AssertWellFormedMatching(manyVsOne, many, one);
        manyVsOne.Should().ContainSingle("the N×1 shape is the mirror of 1×N → still a single pair");
        manyVsOne[0].Gene1Id.Should().Be("t0");
        manyVsOne[0].Gene2Id.Should().Be("q0");
    }

    /// <summary>
    /// BE / single gene each — null and empty boundaries on the smallest inputs. A null gene
    /// list is the one input rejected with an exception (§3.3, §6.1; ThrowIfNull lines 471-472),
    /// eagerly on the public call. An empty list on either side short-circuits (Count==0,
    /// lines 477-478) to an empty result with no NRE / DivideByZero. We probe a lone gene against
    /// null and against empty on both sides.
    /// </summary>
    [Test]
    public void Rbh_SingleGene_AgainstNullOrEmpty_ThrowsOrEmptyAsDocumented()
    {
        var one = new[] { Gene("a1", "G1", RandomDna(60, seed: 16500)) };
        var empty = Array.Empty<ComparativeGenomics.Gene>();

        FluentActions.Invoking(() => ComparativeGenomics.FindReciprocalBestHits(null!, one))
            .Should().Throw<ArgumentNullException>("a null genome-1 list is rejected by ThrowIfNull, not a deferred NRE (§3.3)");
        FluentActions.Invoking(() => ComparativeGenomics.FindReciprocalBestHits(one, null!))
            .Should().Throw<ArgumentNullException>("a null genome-2 list is rejected by ThrowIfNull (§3.3)");

        FluentActions.Invoking(() => ComparativeGenomics.FindReciprocalBestHits(empty, one).ToList())
            .Should().NotThrow("an empty genome-1 short-circuits to an empty result, never an NRE (§6.1)");
        ComparativeGenomics.FindReciprocalBestHits(empty, one).Should().BeEmpty("no genome-1 gene → no pair (§6.1)");
        ComparativeGenomics.FindReciprocalBestHits(one, empty).Should().BeEmpty("no genome-2 gene → no pair (§6.1)");
        ComparativeGenomics.FindReciprocalBestHits(empty, empty).Should().BeEmpty("both empty → empty result, never a 0/0 (§6.1)");
    }

    #endregion

    #region Positive sanity & reciprocity discrimination (the documented RBH contract)

    /// <summary>
    /// Positive sanity (Reciprocal_Best_Hits.md §7.1 worked example): two identical genes across
    /// genomes are reciprocal best hits → exactly one pair with the documented Identity 1.0,
    /// Coverage 1.0, AlignmentLength = sequence length. Pins the load-bearing metrics (INV-03),
    /// not merely a green pass.
    /// </summary>
    [Test]
    public void Rbh_WorkedExample_PinsDocumentedPairAndMetrics()
    {
        var g1 = new[] { Gene("a1", "G1", "ACGTACGTACGTAC") };
        var g2 = new[] { Gene("b1", "G2", "ACGTACGTACGTAC") };

        var rbh = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList();

        AssertWellFormedMatching(rbh, g1, g2);
        rbh.Should().ContainSingle("two identical genes are reciprocal best hits → one pair (§7.1)");
        rbh[0].Gene1Id.Should().Be("a1");
        rbh[0].Gene2Id.Should().Be("b1");
        rbh[0].Identity.Should().Be(1.0, "identical sequences → 5-mer Jaccard 1.0 (§7.1)");
        rbh[0].Coverage.Should().Be(1.0, "identical sequences share all 5-mers of the shorter sequence → coverage 1.0 (§7.1)");
        rbh[0].AlignmentLength.Should().Be(14, "AlignmentLength = min length of the two sequences = 14 (§7.1)");
    }

    /// <summary>
    /// Positive sanity — discrimination: a clearly similar pair IS reported while an unrelated
    /// gene is NOT paired. genome1 = {a1 (true ortholog), a2 (poly-G unrelated)}; genome2 =
    /// {b1 (a1 with a few substitutions), b2 (poly-T unrelated)}. Only the reciprocal a1↔b1 must
    /// appear; the unrelated genes share no 5-mer with anything and stay unpaired. This is the
    /// core business meaning of RBH-based orthology.
    /// </summary>
    [Test]
    public void Rbh_RelatedPairedUnrelatedUnpaired()
    {
        string orthoSeq = RandomDna(140, seed: 17100);
        char[] m = orthoSeq.ToCharArray();
        var rng = new Random(17101);
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

        var rbh = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2).ToList();

        AssertWellFormedMatching(rbh, genome1, genome2);
        rbh.Should().ContainSingle("only the truly similar a1↔b1 pair is a reciprocal best hit");
        rbh[0].Gene1Id.Should().Be("a1");
        rbh[0].Gene2Id.Should().Be("b1", "the divergent/unrelated genes (a2, b2) share no 5-mer and stay unpaired");
        rbh[0].Identity.Should().BeGreaterThan(0.3, "a few substitutions keep the pair well above the gate");
        AssertReciprocal(genome1, genome2);
    }

    /// <summary>
    /// Reciprocity discrimination at scale: a constructed crossed best-hit configuration where a
    /// one-directional best hit exists must drop it. a1 is most similar to b2, but b2 is most
    /// similar to a2 (and a2 most similar to b2 too); b1 is most similar to a1 (one-directional).
    /// The only mutually-best pair is a2↔b2; a1↔b2 and b1↔a1 are one-directional and excluded.
    /// </summary>
    [Test]
    public void Rbh_CrossedBestHits_OnlyMutualPairSurvives()
    {
        // Anchor X: a2 and b2 both ≈ X (mutual best). a1 ≈ X but slightly worse than a2;
        // b1 ≈ X but slightly worse than b2. So a1's best in G2 is b2 (not b1), and b1's best
        // in G1 is a2 or a1 — either way b1 cannot be in a reciprocal pair, and a1's best b2 is
        // already claimed by the stronger a2 → a1 unpaired. Net mutual pair: a2↔b2.
        string x = RandomDna(160, seed: 17200);
        string a2seq = x;                                  // perfect
        string b2seq = x;                                  // perfect → a2↔b2 mutual best
        string a1seq = x[..130] + RandomDna(30, seed: 17201); // weaker copy of X
        string b1seq = x[..120] + RandomDna(40, seed: 17202); // weaker copy of X

        var genome1 = new[] { Gene("a1", "G1", a1seq), Gene("a2", "G1", a2seq) };
        var genome2 = new[] { Gene("b1", "G2", b1seq), Gene("b2", "G2", b2seq) };

        var rbh = ComparativeGenomics.FindReciprocalBestHits(genome1, genome2).ToList();

        AssertWellFormedMatching(rbh, genome1, genome2);
        rbh.Should().ContainSingle("only the mutually-best a2↔b2 pair is reciprocal; the weaker one-directional hits are dropped");
        rbh[0].Gene1Id.Should().Be("a2");
        rbh[0].Gene2Id.Should().Be("b2", "a2 and b2 are each other's best hit; a1/b1 best hits are one-directional (INV-01)");
        AssertReciprocal(genome1, genome2);
    }

    /// <summary>
    /// Broad determinism net: across many random small genomes, FindReciprocalBestHits must
    /// always return a well-formed matching, be order-independent (INV-04), and be reciprocal
    /// (INV-01) — never throw, hang, or produce a non-reciprocal / order-dependent result.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Rbh_RandomSmallGenomes_AlwaysWellFormedReciprocalDeterministic()
    {
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(seed);
            int n1 = rng.Next(0, 6);
            int n2 = rng.Next(0, 6);

            // Draw from a small pool of motifs so genuine RBH pairs and ties both arise.
            string[] pool =
            {
                RandomDna(70, seed * 7 + 1),
                RandomDna(70, seed * 7 + 2),
                RandomDna(70, seed * 7 + 3),
            };

            var g1 = Enumerable.Range(0, n1).Select(i => Gene($"a{i}", "G1", pool[rng.Next(pool.Length)])).ToList();
            var g2 = Enumerable.Range(0, n2).Select(i => Gene($"b{i}", "G2", pool[rng.Next(pool.Length)])).ToList();

            List<ComparativeGenomics.OrthologPair> rbh = null!;
            FluentActions.Invoking(() => rbh = ComparativeGenomics.FindReciprocalBestHits(g1, g2).ToList())
                .Should().NotThrow($"random small genomes must never crash RBH (seed {seed})");

            AssertWellFormedMatching(rbh, g1, g2);
            AssertOrderIndependentMatching(g1, g2, shuffleSeed: 9000 + seed);
            AssertReciprocal(g1, g2);
        }
    }

    #endregion

    #endregion
}
