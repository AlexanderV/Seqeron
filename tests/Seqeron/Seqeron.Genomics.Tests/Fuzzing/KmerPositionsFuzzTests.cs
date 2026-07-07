namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the K-mer area — K-MER POSITIONS / pattern location
/// (KMER-POSITIONS-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain values to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception (IndexOutOfRangeException / ArgumentOutOfRangeException leaking from
/// internal indexing, OutOfMemoryException). Every input must resolve to EITHER a
/// well-defined, theory-correct result, OR a *documented, intentional* validation
/// exception. For this unit the documented contract is total — null/empty inputs
/// and an over-long k-mer all yield an EMPTY result, never a throw — so a raw
/// runtime exception, a hang, a false position, a missed overlapping occurrence,
/// or an off-by-one at the last valid start is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-POSITIONS-001 — all start positions of a k-mer in a sequence
/// Checklist: docs/checklists/03_FUZZING.md, row 160.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: an ABSENT k-mer (no match → empty list, no false
///          position), OVERLAPPING occurrences (the headline correctness rule —
///          "AA" in "AAAA" → {0,1,2}, every overlapping start reported), and
///          k &gt; len (the k-mer longer than the sequence → empty list, never an
///          out-of-range Substring). We also pin the exact off-by-one edges
///          k = L (single match at 0) and the last valid start L − k.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The k-mer-positions contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// The API entry is
///   KmerAnalyzer.FindKmerPositions(string sequence, string kmer)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines 432–447).
/// It solves the Pattern Matching Problem (Rosalind BA1D; Compeau &amp; Pevzner,
/// Bioinformatics Algorithms): report every 0-based start position `p` where
/// `sequence[p .. p+k)` equals the k-mer (K-mer_Positions.md §2.2). The documented
/// invariants (K-mer_Positions.md §2.4) are:
///   • INV-01 — every returned position `p` satisfies `T[p .. p+k) = P` (case-folded).
///   • INV-02 — positions are 0-based and strictly ascending.
///   • INV-03 — the count of returned positions equals the OVERLAPPING occurrence
///              count of `P` in `T` (every candidate window is tested; overlaps
///              are all reported — this is the headline BE target).
///   • INV-04 — all positions lie in `[0, L − k]`; the result is EMPTY when k &gt; L.
/// Matching is case-insensitive (both arguments upper-cased via ToUpperInvariant,
/// K-mer_Positions.md §5.2). Null/empty sequence or k-mer, and |kmer| &gt; |sequence|,
/// all yield an EMPTY enumerable — no exception (K-mer_Positions.md §3.3, §6.1).
///
/// The three checklist targets map to these documented behaviours:
///   • absent k-mer  → empty list; only matching starts are reported, no false
///                     position (K-mer_Positions.md §6.1).
///   • overlapping   → every overlapping start reported, e.g. "AA" in "AAAA" →
///                     {0,1,2}, "ATA" in "ATATA" → {0,2} (INV-03; §7.1 walk-through).
///   • k &gt; len      → empty list; the loop bound `i ≤ L − k` is negative so the
///                     scan never runs, no out-of-range Substring (INV-04, §6.1).
/// A positive-sanity test pins the documented Rosalind BA1D worked example
/// ("ATAT" in "GATATATGCATATACTT" → {1,3,9}) and the overlapping example, with the
/// reusable well-formed oracle (every p in [0, L−k], Substring(p,k)==kmer, ascending).
///
/// NOTE ON SOURCE STATE: the unit was independently validated CLEAN
/// (docs/Validation/FINDINGS_REGISTER.md A36); these fuzz tests exercise the BE
/// boundaries and confirm no defect surfaces under degenerate input.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerPositionsFuzzTests
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

    /// <summary>
    /// Asserts the returned positions form a WELL-FORMED occurrence list for the
    /// given sequence and k-mer (the load-bearing structural oracle reused across
    /// the fuzz cases): every position lies in [0, L − k] (INV-04), the substring at
    /// each position equals the k-mer case-insensitively (INV-01), the positions are
    /// strictly ascending and distinct (INV-02), and the COUNT equals the brute-force
    /// overlapping-occurrence count computed by an independent reference scan that
    /// advances by ONE each step (INV-03 — proves no overlapping start was skipped and
    /// no false position invented). — K-mer_Positions.md §2.4.
    /// </summary>
    private static void AssertWellFormedPositions(
        IReadOnlyList<int> positions, string sequence, string kmer)
    {
        string seq = (sequence ?? string.Empty).ToUpperInvariant();
        string km = (kmer ?? string.Empty).ToUpperInvariant();

        // INV-04: every position is a valid start window [0, L − k].
        positions.Should().OnlyContain(p => p >= 0 && p <= seq.Length - km.Length,
            "INV-04: every reported start lies in [0, L − k]; nothing past the last fitting window");

        // INV-01: the substring at each position equals the k-mer (case-folded).
        foreach (int p in positions)
            seq.Substring(p, km.Length).Should().Be(km,
                $"INV-01: sequence[{p} .. {p}+k) must equal the k-mer at every reported position");

        // INV-02: strictly ascending and therefore distinct.
        positions.Should().BeInAscendingOrder("INV-02: positions are produced in ascending order");
        positions.Should().OnlyHaveUniqueItems("INV-02: ascending positions are necessarily distinct");

        // INV-03: count equals an INDEPENDENT brute-force overlapping reference (step +1).
        var reference = new List<int>();
        if (km.Length > 0 && seq.Length >= km.Length)
            for (int i = 0; i <= seq.Length - km.Length; i++)
                if (seq.Substring(i, km.Length) == km)
                    reference.Add(i);
        positions.Should().Equal(reference,
            "INV-03: the positions must match the independent overlapping reference exactly — every overlapping start, no false position");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-POSITIONS-001 — all start positions of a k-mer : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-POSITIONS-001 — k-mer positions

    #region BE — Boundary: overlapping occurrences (the headline rule)

    /// <summary>
    /// BE / KEY: the headline correctness property — OVERLAPPING occurrences are all
    /// reported. "AA" in "AAAA" (L = 4, k = 2) occurs at starts 0, 1, 2 (the windows
    /// AA|AA, A|AA|A, AA|AA overlap by one base). A non-overlapping scanner bug would
    /// advance past each match and report only {0, 2}; the documented contract
    /// (Rosalind BA1D; INV-03; K-mer_Positions.md §6.1) is {0, 1, 2}. We pin the exact
    /// list so a non-overlapping implementation FAILS.
    /// </summary>
    [Test]
    public void FindKmerPositions_SelfOverlappingKmer_ReportsEveryOverlappingStart()
    {
        var positions = KmerAnalyzer.FindKmerPositions("AAAA", "AA").ToList();

        positions.Should().Equal(new[] { 0, 1, 2 },
            "'AA' in 'AAAA' overlaps; every overlapping start is reported (Rosalind BA1D), not just the non-overlapping {0,2}");

        AssertWellFormedPositions(positions, "AAAA", "AA");
    }

    /// <summary>
    /// BE: the documented overlapping walk-through "ATA" in "ATATA" → {0, 2}
    /// (K-mer_Positions task description: ATA occurs at overlapping offsets). A
    /// non-overlapping scan that consumed the first match (positions 0..2) would miss
    /// the start at 2. We pin {0, 2} exactly.
    /// </summary>
    [Test]
    public void FindKmerPositions_OverlappingTriplet_ReportsBothOverlappingStarts()
    {
        var positions = KmerAnalyzer.FindKmerPositions("ATATA", "ATA").ToList();

        positions.Should().Equal(new[] { 0, 2 },
            "'ATA' in 'ATATA' occurs at overlapping starts 0 and 2; both are reported");

        AssertWellFormedPositions(positions, "ATATA", "ATA");
    }

    /// <summary>
    /// BE: a maximally self-overlapping homopolymer. "A" (k = 1) in a run of 50 'A's
    /// occurs at EVERY index 0..49 — the densest possible overlap. The result must be
    /// the full ascending range with no gap and no duplicate; this stresses that the
    /// single-base window never skips a position.
    /// </summary>
    [Test]
    public void FindKmerPositions_SingleBaseInHomopolymer_ReportsEveryIndex()
    {
        const int n = 50;
        string seq = new string('A', n);

        var positions = KmerAnalyzer.FindKmerPositions(seq, "A").ToList();

        positions.Should().Equal(Enumerable.Range(0, n),
            "a length-1 k-mer in a homopolymer of length n occurs at every index 0..n−1");

        AssertWellFormedPositions(positions, seq, "A");
    }

    #endregion

    #region BE — Boundary: absent k-mer (no false position)

    /// <summary>
    /// BE: an ABSENT k-mer over the SAME alphabet must yield the empty list — no false
    /// position, no crash. "GG" never occurs in "ACTACTACT", so the result is empty
    /// even though every base of the k-mer is present individually
    /// (K-mer_Positions.md §6.1: only matching starts are reported).
    /// </summary>
    [Test]
    public void FindKmerPositions_AbsentKmer_IsEmptyAndDoesNotThrow()
    {
        var act = () => KmerAnalyzer.FindKmerPositions("ACTACTACT", "GG").ToList();
        act.Should().NotThrow("an absent k-mer is a normal no-match outcome, not an error");

        KmerAnalyzer.FindKmerPositions("ACTACTACT", "GG").Should().BeEmpty(
            "'GG' does not occur as a substring; only matching starts are reported — no false position");
    }

    /// <summary>
    /// BE: a near-miss k-mer that shares a prefix with real windows but differs in its
    /// last base must NOT be reported. "ACG" appears in "ACGACG" but "ACT" (same first
    /// two bases) does not — a partial-match bug would falsely report it. We pin empty.
    /// </summary>
    [Test]
    public void FindKmerPositions_PrefixSharingButAbsentKmer_ReportsNoFalsePosition()
    {
        KmerAnalyzer.FindKmerPositions("ACGACG", "ACT").Should().BeEmpty(
            "'ACT' shares the prefix 'AC' with real windows but never matches in full; no partial/false position");

        // And the genuinely present prefix-sharing k-mer IS found, proving the test is not vacuous.
        KmerAnalyzer.FindKmerPositions("ACGACG", "ACG").Should().Equal(new[] { 0, 3 },
            "'ACG' genuinely occurs at 0 and 3");
    }

    #endregion

    #region BE — Boundary: k > sequence length

    /// <summary>
    /// BE: the k-mer LONGER than the sequence must yield the empty list, never an
    /// out-of-range Substring. The loop bound `i ≤ L − k` is negative, so the forward
    /// scan never runs (INV-04; K-mer_Positions.md §6.1). We pin no-throw AND emptiness
    /// at a far-oversized k-mer and at the exact off-by-one boundary k = L + 1.
    /// </summary>
    [Test]
    public void FindKmerPositions_KmerLongerThanSequence_IsEmptyAndDoesNotThrow()
    {
        var act = () => KmerAnalyzer.FindKmerPositions("ACGT", new string('A', 1000)).ToList();
        act.Should().NotThrow(
            "k > L makes the window count L − k + 1 ≤ 0; the scan never runs, so nothing is indexed past the end");

        KmerAnalyzer.FindKmerPositions("ACGT", new string('A', 1000)).Should().BeEmpty(
            "no length-1000 window fits a 4-base sequence; the result is empty, not a crash");

        // k = L + 1 is the exact off-by-one boundary one past the longest fitting window.
        KmerAnalyzer.FindKmerPositions("ACGT", "ACGTA").Should().BeEmpty(
            "k = L + 1 is one base too long; still empty, never an out-of-range Substring");
    }

    #endregion

    #region BE — Boundary: off-by-one edges (k = L and last valid start L − k)

    /// <summary>
    /// BE: k = L is the upper edge where exactly ONE window fits. A k-mer equal to the
    /// whole sequence matches at start 0 only (K-mer_Positions.md §6.1: "kmer equals
    /// whole sequence → [0]"). We pin {0} and confirm a same-length non-matching k-mer
    /// yields empty (so the single-window scan still tests the predicate).
    /// </summary>
    [Test]
    public void FindKmerPositions_KmerEqualsWholeSequence_IsSingletonZeroOrEmpty()
    {
        KmerAnalyzer.FindKmerPositions("ACGT", "ACGT").Should().Equal(new[] { 0 },
            "a k-mer equal to the whole sequence matches once, at start 0");

        KmerAnalyzer.FindKmerPositions("ACGT", "TGCA").Should().BeEmpty(
            "a same-length non-matching k-mer fits one window but does not match it → empty");
    }

    /// <summary>
    /// BE: the LAST valid start is L − k — the classic off-by-one boundary. A k-mer that
    /// occurs ONLY as the trailing window must be found at exactly that index and nowhere
    /// else. In "GGGACT" (L = 6), "ACT" (k = 3) occurs only at start 3 = L − k. An
    /// off-by-one that stopped the scan at L − k − 1 would miss it; one that ran to L − k + 1
    /// would index out of range. We pin {3} exactly.
    /// </summary>
    [Test]
    public void FindKmerPositions_MatchOnlyAtLastWindow_FindsStartLMinusK()
    {
        const string seq = "GGGACT";
        const string km = "ACT";

        var positions = KmerAnalyzer.FindKmerPositions(seq, km).ToList();

        positions.Should().Equal(new[] { seq.Length - km.Length },
            "the trailing-only k-mer occurs exactly at the last valid start L − k = 3");

        AssertWellFormedPositions(positions, seq, km);
    }

    #endregion

    #region BE — Boundary: null / empty inputs

    /// <summary>
    /// BE: null/empty sequence or k-mer is the lower size boundary. All four degenerate
    /// combinations yield the EMPTY enumerable with no exception — the guard returns
    /// before any indexing (K-mer_Positions.md §3.3). We pin no-throw and emptiness for
    /// each.
    /// </summary>
    [Test]
    public void FindKmerPositions_NullOrEmptyInputs_AreEmptyAndDoNotThrow()
    {
        var nullSeq = () => KmerAnalyzer.FindKmerPositions(null!, "AC").ToList();
        var emptySeq = () => KmerAnalyzer.FindKmerPositions(string.Empty, "AC").ToList();
        var nullKmer = () => KmerAnalyzer.FindKmerPositions("ACGT", null!).ToList();
        var emptyKmer = () => KmerAnalyzer.FindKmerPositions("ACGT", string.Empty).ToList();

        nullSeq.Should().NotThrow("null sequence is treated as empty, not an error");
        emptySeq.Should().NotThrow("empty sequence has no windows");
        nullKmer.Should().NotThrow("null k-mer is treated as empty, not an error");
        emptyKmer.Should().NotThrow("empty k-mer is treated as empty, not an error");

        KmerAnalyzer.FindKmerPositions(null!, "AC").Should().BeEmpty();
        KmerAnalyzer.FindKmerPositions(string.Empty, "AC").Should().BeEmpty();
        KmerAnalyzer.FindKmerPositions("ACGT", null!).Should().BeEmpty();
        KmerAnalyzer.FindKmerPositions("ACGT", string.Empty).Should().BeEmpty();
    }

    #endregion

    #region BE — Injection-adjacent: non-DNA characters matched literally

    /// <summary>
    /// BE: the alphabet is unrestricted text (no DNA validation; K-mer_Positions.md §6.2).
    /// A k-mer of arbitrary characters must be located literally (after case-folding)
    /// without crashing — digits, punctuation and whitespace are just characters. We pin
    /// an overlapping match of "%%" in "%%%" → {0, 1} to combine the no-validation rule
    /// with the overlapping rule.
    /// </summary>
    [Test]
    public void FindKmerPositions_NonDnaCharacters_AreMatchedLiterallyAndOverlap()
    {
        var positions = KmerAnalyzer.FindKmerPositions("%%%", "%%").ToList();

        positions.Should().Equal(new[] { 0, 1 },
            "no alphabet validation: '%%' is matched literally and its overlapping starts 0,1 are reported");

        AssertWellFormedPositions(positions, "%%%", "%%");
    }

    #endregion

    #region Positive sanity — documented worked examples

    /// <summary>
    /// Positive sanity / KEY: the documented Rosalind BA1D worked example
    /// (K-mer_Positions.md §7.1) — "ATAT" in "GATATATGCATATACTT" → {1, 3, 9}. This pins
    /// overlapping detection (1 and 3 overlap by two bases) on a real, externally-sourced
    /// example. A 1-based or non-overlapping implementation would FAIL ({2,4,10} or
    /// {1,9}). Every reported position is verified to satisfy Substring(p,k)==kmer.
    /// </summary>
    [Test]
    public void FindKmerPositions_RosalindBa1dSample_MatchesDocumentedPositions()
    {
        const string seq = "GATATATGCATATACTT";
        const string km = "ATAT";

        var positions = KmerAnalyzer.FindKmerPositions(seq, km).ToList();

        positions.Should().Equal(new[] { 1, 3, 9 },
            "Rosalind BA1D sample: 'ATAT' in 'GATATATGCATATACTT' → 1 3 9 (0-based, overlapping)");

        AssertWellFormedPositions(positions, seq, km);
    }

    /// <summary>
    /// Positive sanity: case-insensitivity (K-mer_Positions.md §5.2). A lowercase k-mer
    /// and/or lowercase sequence must locate the same positions as the uppercase form,
    /// because both arguments are upper-cased before matching.
    /// </summary>
    [Test]
    public void FindKmerPositions_MixedCase_IsCaseInsensitive()
    {
        var upper = KmerAnalyzer.FindKmerPositions("ACGTACGT", "ACG").ToList();
        var lowerSeq = KmerAnalyzer.FindKmerPositions("acgtacgt", "ACG").ToList();
        var lowerKmer = KmerAnalyzer.FindKmerPositions("ACGTACGT", "acg").ToList();
        var mixed = KmerAnalyzer.FindKmerPositions("AcGtAcGt", "aCg").ToList();

        upper.Should().Equal(new[] { 0, 4 }, "'ACG' occurs at 0 and 4 in 'ACGTACGT'");
        lowerSeq.Should().Equal(upper, "lowercasing the sequence does not change the positions");
        lowerKmer.Should().Equal(upper, "lowercasing the k-mer does not change the positions");
        mixed.Should().Equal(upper, "mixed case on both sides is folded before matching");
    }

    /// <summary>
    /// Positive sanity: lazy enumeration must be stable and repeatable. FindKmerPositions
    /// is implemented as a lazy IEnumerable&lt;int&gt; (yield); re-enumerating the same
    /// query must produce the identical ascending list (deterministic, no shared state).
    /// </summary>
    [Test]
    public void FindKmerPositions_ReEnumeration_IsDeterministic()
    {
        var query = KmerAnalyzer.FindKmerPositions("ATATATAT", "AT");

        var first = query.ToList();
        var second = query.ToList();

        first.Should().Equal(new[] { 0, 2, 4, 6 }, "'AT' in 'ATATATAT' occurs at 0,2,4,6");
        second.Should().Equal(first, "re-enumerating the lazy sequence yields the identical result");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// satisfy the well-formed invariants for several k-mers drawn from inside it AND for
    /// an absent k-mer — the count equals the independent overlapping reference, every
    /// position is in-bounds with a matching substring, and positions are ascending —
    /// regardless of the random content. [CancelAfter] guards against any hang.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindKmerPositions_RandomSequence_SatisfiesPositionInvariants()
    {
        const int length = 4000;
        string seq = RandomDna(length, seed: 160_001);

        foreach (int k in new[] { 1, 2, 3, 5, 8 })
        {
            // A k-mer guaranteed PRESENT (a window taken from the sequence itself).
            string present = seq.Substring(length / 2, k);
            var hits = KmerAnalyzer.FindKmerPositions(seq, present).ToList();
            hits.Should().NotBeEmpty($"a k-mer copied out of the sequence (k = {k}) must occur at least once");
            AssertWellFormedPositions(hits, seq, present);

            // An ABSENT k-mer (a non-DNA filler) → empty list, no false position.
            string absent = new string('Z', k);
            KmerAnalyzer.FindKmerPositions(seq, absent).Should().BeEmpty(
                $"a k-mer of '{absent}' (over an absent symbol) occurs nowhere in a DNA sequence");
        }
    }

    #endregion

    #endregion
}
