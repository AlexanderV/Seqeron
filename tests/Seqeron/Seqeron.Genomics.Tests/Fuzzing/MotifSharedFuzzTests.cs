using System.Text;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Matching area — Shared-Motif finding (MOTIF-SHARED-001),
/// the multi-sequence fixed-length word enumerator with a matching-sequence quorum
/// <see cref="MotifFinder.FindSharedMotifs(IEnumerable{DnaSequence}, int, int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain collections/parameters to a
/// unit and asserts the code NEVER fails in an undisciplined way: no hang/infinite
/// loop (the per-sequence O(Σᵢ nᵢ·k) window scan must always terminate), no state
/// corruption, no nonsense output (a reported word whose length is not k, a word that
/// does NOT actually occur in every index it reports, a word below the quorum, a
/// duplicate index for one word, a Prevalence outside (0,1] or that disagrees with the
/// index count / total), a NON-DETERMINISTIC result set, and no *unhandled* runtime
/// exception — in particular NO NullReference/crash on a one-element or empty collection
/// and NO DivideByZero / Infinity / NaN Prevalence on the empty-collection boundary.
/// Every input must resolve to EITHER a well-defined, theory-correct result OR a
/// *documented, intentional* validation exception (ArgumentNullException for a null
/// collection; ArgumentOutOfRangeException for k &lt; 1 or minSequences &lt; 1 — contract
/// §3.3, §6.1). A raw runtime exception, a hang, a false shared word, a wrong index set,
/// or a non-deterministic output is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MOTIF-SHARED-001 — Shared Motifs (word enumeration, matching-sequence quorum)
/// Checklist: docs/checklists/03_FUZZING.md, row 173.
/// Algorithm doc: docs/algorithms/Motif_Discovery/Shared_Motifs.md
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the row:
///       – single input: a collection of ONE sequence. With the default quorum
///         minSequences=2 nothing can meet the quorum (one sequence cannot reach 2),
///         so the result is EMPTY — never a NullReference/crash (§6.1). With
///         minSequences=1 every distinct length-k word of that one sequence is
///         trivially "shared" (matching count 1 ≥ 1), Prevalence = 1/1 = 1.0,
///         with NO DivideByZero (INV-03, INV-04).
///       – disjoint inputs: sequences that share NO common length-k word → an EMPTY
///         shared set, NO false shared word fabricated (INV-03, INV-02).
///       – identical: all inputs are the SAME sequence. Every distinct length-k word of
///         that sequence is shared by all m inputs; SequenceIndices = {0..m−1} (each
///         index once, no double-count), Prevalence = m/m = 1.0; the output is
///         DETERMINISTIC across repeated calls (INV-02, INV-04).
///       – k &lt; 1 / minSequences &lt; 1 (0, −1, int.MinValue): the documented
///         ArgumentOutOfRangeException, never an empty/garbage result (§3.3, §6.1).
///       – null collection: the documented ArgumentNullException (§3.3, §6.1).
/// — docs/checklists/03_FUZZING.md §Description (BE = граничні значення: 0, -1, MaxInt, empty).
///
/// Note on Malformed Content: each element is a <see cref="DnaSequence"/>, which is
/// uppercased and validated to the {A,C,G,T} alphabet at construction, so out-of-domain
/// residues cannot reach this method; this is therefore a pure boundary (BE) row over
/// the collection shape (single / disjoint / identical / empty) and the integer
/// parameters k / minSequences, exactly as the checklist row specifies.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Shared_Motifs.md §2.2, §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// For input sequences S = {s_0 … s_{m−1}} and word length k, every distinct length-k
/// word w is scored by its matching-sequence set M(w) = { i : s_i contains ≥ 1 exact
/// occurrence of w } (a word repeated within one sequence still contributes 1 — INV-02).
/// A word is reported iff |M(w)| ≥ minSequences (INV-03). Each reported word has length
/// exactly k (INV-01), distinct sequence indices (INV-02), exact matching — a
/// 1-substitution variant is a different word (INV-05), and Prevalence = |M(w)|/m ∈ (0,1]
/// (INV-04). Validation (§3.3): null collection → ArgumentNullException; k &lt; 1 or
/// minSequences &lt; 1 → ArgumentOutOfRangeException; empty collection → no results; a
/// sequence shorter than k contributes no words. Deterministic (§1).
///   MotifFinder.FindSharedMotifs(IEnumerable&lt;DnaSequence&gt;, int k = 6, int minSequences = 2)
///       → IEnumerable&lt;SharedMotif(Sequence, SequenceIndices, Prevalence)&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MotifSharedFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>A random ACGT string of the given length.</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(Alphabet[rng.Next(Alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent oracle implementing the documented decision rule verbatim
    /// (Shared_Motifs.md §2.2, §4.1): for each sequence collect its DISTINCT length-k
    /// words (per-sequence presence), accumulate word → sorted set of containing indices,
    /// and keep every word whose matching-sequence count ≥ minSequences. Returns
    /// word → (sorted distinct indices, prevalence). Built from the spec, not the unit.
    /// </summary>
    private static Dictionary<string, (List<int> Indices, double Prevalence)> Oracle(
        IReadOnlyList<string> seqs, int k, int minSequences)
    {
        var matching = new Dictionary<string, SortedSet<int>>();
        for (int i = 0; i < seqs.Count; i++)
        {
            string s = seqs[i];
            var seen = new HashSet<string>();
            for (int p = 0; p <= s.Length - k; p++)
                seen.Add(s.Substring(p, k));
            foreach (var w in seen)
            {
                if (!matching.TryGetValue(w, out var set))
                    matching[w] = set = new SortedSet<int>();
                set.Add(i);
            }
        }

        var result = new Dictionary<string, (List<int>, double)>();
        foreach (var (w, set) in matching)
            if (set.Count >= minSequences)
                result[w] = (set.ToList(), (double)set.Count / seqs.Count);
        return result;
    }

    /// <summary>
    /// Asserts a shared-motif result is WELL-FORMED per the documented contract,
    /// independent of the (possibly degenerate) input:
    ///   • every reported word has length exactly k (INV-01);
    ///   • word strings are DISTINCT (one record per word);
    ///   • SequenceIndices are distinct, in range [0, m), and at least minSequences of
    ///     them (INV-02, INV-03 — the quorum is honoured, no false shared word);
    ///   • the word GENUINELY occurs in every index it reports — exact match (INV-05);
    ///   • Prevalence is finite, equals |indices|/m, and lies in (0,1] (INV-04 — no
    ///     DivideByZero / NaN / Infinity, no value &gt; 1).
    /// </summary>
    private static void AssertWellFormed(
        IReadOnlyList<SharedMotif> motifs, IReadOnlyList<string> seqs, int k, int minSequences)
    {
        int m = seqs.Count;
        motifs.Select(x => x.Sequence).Should().OnlyHaveUniqueItems("one record per distinct word");

        foreach (var sm in motifs)
        {
            sm.Sequence.Should().HaveLength(k, "every reported word is a length-k oligonucleotide (INV-01)");

            sm.SequenceIndices.Should().OnlyHaveUniqueItems("each sequence counted at most once per word (INV-02)");
            sm.SequenceIndices.Should().OnlyContain(i => i >= 0 && i < m, "indices are valid 0-based positions in the input");
            sm.SequenceIndices.Count.Should().BeGreaterThanOrEqualTo(minSequences,
                "INV-03: only words meeting the matching-sequence quorum are reported");

            foreach (int i in sm.SequenceIndices)
                seqs[i].Should().Contain(sm.Sequence,
                    "INV-05: a reported word genuinely occurs (exact match) in every index it claims");

            double.IsFinite(sm.Prevalence).Should().BeTrue("INV-04: prevalence is finite (no DivideByZero/NaN/Infinity)");
            sm.Prevalence.Should().BeApproximately((double)sm.SequenceIndices.Count / m, 1e-12,
                "INV-04: prevalence == |M(w)| / m");
            sm.Prevalence.Should().BeInRange(double.Epsilon, 1.0, "INV-04: prevalence ∈ (0, 1]");
        }
    }

    /// <summary>
    /// Asserts the unit's reported word SET matches the independent oracle exactly
    /// (same words, same index sets, same prevalence) — the strict theory-correct
    /// cross-check used by the fuzz loops.
    /// </summary>
    private static void AssertMatchesOracle(
        IReadOnlyList<SharedMotif> motifs, IReadOnlyList<string> seqs, int k, int minSequences)
    {
        var oracle = Oracle(seqs, k, minSequences);

        motifs.Select(x => x.Sequence).Should().BeEquivalentTo(oracle.Keys,
            "the reported word set equals the documented |M(w)| ≥ minSequences set (INV-03)");

        foreach (var sm in motifs)
        {
            var (indices, prevalence) = oracle[sm.Sequence];
            sm.SequenceIndices.OrderBy(i => i).Should().Equal(indices,
                "SequenceIndices match the documented matching-sequence set (INV-02)");
            sm.Prevalence.Should().BeApproximately(prevalence, 1e-12, "Prevalence matches |M(w)|/m (INV-04)");
        }
    }

    #endregion

    #region MOTIF-SHARED-001 — Shared Motifs (BE: single input, disjoint inputs, identical)

    #region Positive sanity — hand-computed documented sharing

    // Documented worked example (§7.1): seqs = {ATGATG, ATGCCC, CCCGGG}, k=3, minSequences=2.
    // ATG in seq0 (twice → counts once) and seq1 → indices [0,1], prevalence 2/3.
    // CCC in seq1 and seq2 → indices [1,2], prevalence 2/3. GGG only seq2 → excluded.
    [Test]
    public void FindSharedMotifs_DocumentedWorkedExample_ReportsAtgAndCccOnly()
    {
        var seqs = new[] { new DnaSequence("ATGATG"), new DnaSequence("ATGCCC"), new DnaSequence("CCCGGG") };

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();

        var atg = shared.Single(s => s.Sequence == "ATG");
        atg.SequenceIndices.OrderBy(i => i).Should().Equal(new[] { 0, 1 }, "ATG occurs in seq0 (twice→once) and seq1 (INV-02)");
        atg.Prevalence.Should().BeApproximately(2.0 / 3.0, 1e-12, "matching 2 of 3 sequences (INV-04)");

        var ccc = shared.Single(s => s.Sequence == "CCC");
        ccc.SequenceIndices.OrderBy(i => i).Should().Equal(1, 2);
        ccc.Prevalence.Should().BeApproximately(2.0 / 3.0, 1e-12);

        shared.Should().NotContain(s => s.Sequence == "GGG", "GGG occurs in only 1 sequence < quorum 2 (INV-03)");
        AssertWellFormed(shared, new[] { "ATGATG", "ATGCCC", "CCCGGG" }, 3, 2);
    }

    // A common substring planted in every input ("ACGTACGT") is reported as shared.
    [Test]
    public void FindSharedMotifs_CommonSubstringInAllInputs_IsReportedShared()
    {
        // Every sequence contains the 8-mer ACGTACGT at distinct flanks.
        var raw = new[] { "TTACGTACGTAA", "GGGACGTACGTC", "ACGTACGTGGGG" };
        var seqs = raw.Select(s => new DnaSequence(s)).ToArray();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 8, minSequences: 3).ToList();

        shared.Should().ContainSingle("only ACGTACGT is shared by all three at k=8");
        var w = shared[0];
        w.Sequence.Should().Be("ACGTACGT");
        w.SequenceIndices.OrderBy(i => i).Should().Equal(new[] { 0, 1, 2 }, "present in every input (INV-02)");
        w.Prevalence.Should().Be(1.0, "shared by all m sequences ⇒ prevalence 1 (INV-04)");
        AssertWellFormed(shared, raw, 8, 3);
    }

    #endregion

    #region BE — Boundary: single input (one-element collection → no crash)

    // Single sequence with the default quorum (minSequences=2): one sequence can never
    // reach a quorum of 2 → EMPTY, never a NullReference/crash (§6.1).
    [Test]
    public void FindSharedMotifs_SingleInput_DefaultQuorum_ReturnsEmpty_NoCrash()
    {
        var seqs = new[] { new DnaSequence("ACGTACGTACGT") };

        Action act = () => MotifFinder.FindSharedMotifs(seqs, k: 4).ToList();
        act.Should().NotThrow("a one-element collection is a documented boundary, not an error (§6.1)");

        MotifFinder.FindSharedMotifs(seqs, k: 4)
            .Should().BeEmpty("a single sequence cannot reach the default quorum of 2 (INV-03)");
    }

    // Single sequence with minSequences=1: every DISTINCT length-k word is trivially
    // shared (matching count 1 ≥ 1), each with SequenceIndices=[0], Prevalence=1/1=1.0,
    // and NO DivideByZero (the m=1 denominator).
    [Test]
    public void FindSharedMotifs_SingleInput_Quorum1_AllDistinctWordsSharedWithPrevalenceOne()
    {
        const string seq = "ACGTACG"; // distinct 3-mers: ACG, CGT, GTA, TAC, (ACG again) → {ACG,CGT,GTA,TAC}
        var seqs = new[] { new DnaSequence(seq) };

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 1).ToList();

        shared.Select(s => s.Sequence).Should().BeEquivalentTo(
            new[] { "ACG", "CGT", "GTA", "TAC" }, "every distinct 3-mer of the lone sequence is trivially shared");
        shared.Should().OnlyContain(s => s.SequenceIndices.Count == 1 && s.SequenceIndices[0] == 0,
            "the only containing index is 0 (INV-02)");
        shared.Should().OnlyContain(s => s.Prevalence == 1.0, "1/1 = 1, no DivideByZero (INV-04)");
        AssertWellFormed(shared, new[] { seq }, 3, 1);
        AssertMatchesOracle(shared, new[] { seq }, 3, 1);
    }

    // Fuzz: a one-element collection over random params never throws and is well-formed.
    [Test]
    [CancelAfter(30_000)]
    public void FindSharedMotifs_SingleInput_RandomParams_NeverThrows_MatchesOracle()
    {
        var rng = new Random(173_001);
        for (int trial = 0; trial < 600; trial++)
        {
            int n = rng.Next(0, 40);
            int k = rng.Next(1, 12);
            int minSeq = rng.Next(1, 4);
            string s = RandomDna(rng, n);
            var seqs = new[] { new DnaSequence(s) };

            var shared = MotifFinder.FindSharedMotifs(seqs, k: k, minSequences: minSeq).ToList();

            // Quorum ≥ 2 is unreachable with one sequence ⇒ must be empty; quorum 1 ⇒ all distinct words.
            if (minSeq >= 2)
                shared.Should().BeEmpty("a single sequence cannot reach quorum ≥ 2 (INV-03)");
            AssertWellFormed(shared, new[] { s }, k, minSeq);
            AssertMatchesOracle(shared, new[] { s }, k, minSeq);
        }
    }

    #endregion

    #region BE — Boundary: disjoint inputs (no common word → empty, no false shared word)

    // Two sequences over disjoint k-mer spaces (all-A vs all-C) share NO k-mer → empty.
    [Test]
    public void FindSharedMotifs_DisjointInputs_NoCommonWord_ReturnsEmpty()
    {
        var seqs = new[] { new DnaSequence("AAAAAAAA"), new DnaSequence("CCCCCCCC") };

        MotifFinder.FindSharedMotifs(seqs, k: 4, minSequences: 2)
            .Should().BeEmpty("AAAA-only and CCCC-only sequences share no 4-mer (INV-03)");
    }

    // Three sequences each built from a different residue: still no shared k-mer.
    [Test]
    public void FindSharedMotifs_ThreeDisjointHomopolymers_ReturnsEmpty()
    {
        var raw = new[] { "AAAAAA", "GGGGGG", "TTTTTT" };
        var seqs = raw.Select(s => new DnaSequence(s)).ToArray();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();
        shared.Should().BeEmpty("three single-residue homopolymers share no 3-mer (INV-03)");
        AssertWellFormed(shared, raw, 3, 2); // vacuously well-formed (no false word)
    }

    // Fuzz: sequences drawn so that no word is shared (disjoint residue sub-alphabets per
    // sequence) → never a false shared word, never a crash.
    [Test]
    [CancelAfter(30_000)]
    public void FindSharedMotifs_DisjointResidueSpaces_NeverEmitsFalseSharedWord()
    {
        var rng = new Random(173_002);
        for (int trial = 0; trial < 500; trial++)
        {
            // Seq0 uses {A}, Seq1 uses {C}, Seq2 uses {G}: pairwise-disjoint k-mer spaces.
            var raw = new[]
            {
                new string('A', rng.Next(1, 20)),
                new string('C', rng.Next(1, 20)),
                new string('G', rng.Next(1, 20)),
            };
            var seqs = raw.Select(s => new DnaSequence(s)).ToArray();
            int k = rng.Next(1, 8);
            int minSeq = rng.Next(2, 4);

            var shared = MotifFinder.FindSharedMotifs(seqs, k: k, minSequences: minSeq).ToList();

            shared.Should().BeEmpty("residue-disjoint sequences can share no word (INV-03)");
            AssertMatchesOracle(shared, raw, k, minSeq);
        }
    }

    #endregion

    #region BE — Boundary: identical inputs (every word shared, deterministic, no double-count)

    // m identical copies: every distinct length-k word of that sequence is shared by all
    // m inputs; SequenceIndices = {0..m−1} exactly once each; Prevalence = m/m = 1.0.
    [Test]
    public void FindSharedMotifs_IdenticalInputs_EveryWordSharedByAll_NoDoubleCount()
    {
        const string seq = "ACGTACGTTT";
        var seqs = Enumerable.Range(0, 4).Select(_ => new DnaSequence(seq)).ToArray();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 4, minSequences: 2).ToList();

        // Distinct 4-mers of seq: ACGT, CGTA, GTAC, TACG, ACGT(again), CGTT, GTTT.
        var oracleWords = Oracle(new[] { seq }, 4, 1).Keys; // distinct 4-mers of the one sequence
        shared.Select(s => s.Sequence).Should().BeEquivalentTo(oracleWords,
            "every distinct word of the identical sequence is shared (INV-03)");
        shared.Should().OnlyContain(s => s.SequenceIndices.Count == 4, "shared by all 4 identical inputs");
        shared.Should().OnlyContain(s => s.SequenceIndices.Distinct().Count() == s.SequenceIndices.Count,
            "no index double-counted (INV-02)");
        shared.Should().OnlyContain(s => s.Prevalence == 1.0, "m/m = 1 (INV-04)");
        AssertWellFormed(shared, Enumerable.Repeat(seq, 4).ToList(), 4, 2);
    }

    // Determinism: repeated calls on identical inputs produce the SAME word set, indices,
    // and prevalence — no non-deterministic ordering of results affecting equality.
    [Test]
    public void FindSharedMotifs_IdenticalInputs_DeterministicAcrossCalls()
    {
        var seqs = Enumerable.Range(0, 3).Select(_ => new DnaSequence("ACGTACGTACGT")).ToArray();

        var a = MotifFinder.FindSharedMotifs(seqs, k: 5, minSequences: 2).ToList();
        var b = MotifFinder.FindSharedMotifs(seqs, k: 5, minSequences: 2).ToList();

        a.Select(x => x.Sequence).OrderBy(x => x)
            .Should().Equal(b.Select(x => x.Sequence).OrderBy(x => x), "result set is deterministic (§1)");
        foreach (var w in a)
        {
            var match = b.Single(x => x.Sequence == w.Sequence);
            w.SequenceIndices.OrderBy(i => i).Should().Equal(match.SequenceIndices.OrderBy(i => i),
                "index sets are deterministic");
            w.Prevalence.Should().Be(match.Prevalence, "prevalence is deterministic");
        }
    }

    // Fuzz: m identical copies of a random sequence → all distinct words shared, prevalence 1.
    [Test]
    [CancelAfter(30_000)]
    public void FindSharedMotifs_IdenticalInputs_RandomParams_AllWordsSharedPrevalenceOne()
    {
        var rng = new Random(173_003);
        for (int trial = 0; trial < 500; trial++)
        {
            int n = rng.Next(0, 40);
            int k = rng.Next(1, 12);
            int m = rng.Next(2, 6);
            int minSeq = rng.Next(1, m + 1); // reachable quorum (≤ m)
            string s = RandomDna(rng, n);
            var raw = Enumerable.Repeat(s, m).ToList();
            var seqs = raw.Select(x => new DnaSequence(x)).ToArray();

            var shared = MotifFinder.FindSharedMotifs(seqs, k: k, minSequences: minSeq).ToList();

            // All distinct length-k words of s qualify (each present in all m copies ≥ minSeq).
            var distinctWords = Oracle(new[] { s }, k, 1).Keys;
            shared.Select(x => x.Sequence).Should().BeEquivalentTo(distinctWords,
                "every distinct word of identical inputs is shared by all m copies (INV-03)");
            shared.Should().OnlyContain(x => x.SequenceIndices.Count == m, "present in every identical copy");
            shared.Should().OnlyContain(x => x.Prevalence == 1.0, "m/m = 1 (INV-04)");
            AssertWellFormed(shared, raw, k, minSeq);
        }
    }

    #endregion

    #region BE — Boundary: empty collection (no words → empty, no DivideByZero)

    [Test]
    public void FindSharedMotifs_EmptyCollection_ReturnsEmpty_NoCrash()
    {
        Action act = () => MotifFinder.FindSharedMotifs(Array.Empty<DnaSequence>(), k: 4, minSequences: 1).ToList();
        act.Should().NotThrow("an empty collection has no words to enumerate (§6.1) — no DivideByZero");

        MotifFinder.FindSharedMotifs(Array.Empty<DnaSequence>(), k: 4, minSequences: 1)
            .Should().BeEmpty("no sequences ⇒ no shared words");
    }

    // Sequences all SHORTER than k contribute no words → empty, no crash.
    [Test]
    public void FindSharedMotifs_AllSequencesShorterThanK_ReturnsEmpty()
    {
        var seqs = new[] { new DnaSequence("AC"), new DnaSequence("GT"), new DnaSequence("") };

        MotifFinder.FindSharedMotifs(seqs, k: 6, minSequences: 1)
            .Should().BeEmpty("a sequence shorter than k yields no length-k window (§6.1)");
    }

    #endregion

    #region BE — Boundary: parameter guards (k < 1, minSequences < 1, null → documented exceptions)

    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void FindSharedMotifs_KBelowOne_ThrowsArgumentOutOfRange(int k)
    {
        var seqs = new[] { new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT") };

        Action act = () => MotifFinder.FindSharedMotifs(seqs, k: k, minSequences: 2).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>("k < 1 is the documented validation contract (§3.3, §6.1)")
            .Which.ParamName.Should().Be("k");
    }

    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void FindSharedMotifs_MinSequencesBelowOne_ThrowsArgumentOutOfRange(int minSequences)
    {
        var seqs = new[] { new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT") };

        Action act = () => MotifFinder.FindSharedMotifs(seqs, k: 4, minSequences: minSequences).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
            "minSequences < 1 is the documented validation contract (§3.3, §6.1)")
            .Which.ParamName.Should().Be("minSequences");
    }

    [Test]
    public void FindSharedMotifs_NullCollection_Throws()
    {
        Action act = () => MotifFinder.FindSharedMotifs(null!, k: 4, minSequences: 2).ToList();

        act.Should().Throw<ArgumentNullException>("null collection is the documented validation contract (§3.3)");
    }

    #endregion

    #region BE — Broad fuzz: random collections / k / quorum never crash, match the documented rule

    [Test]
    [CancelAfter(60_000)]
    public void FindSharedMotifs_RandomCollections_NeverThrows_MatchesOracle()
    {
        var rng = new Random(173_004);
        for (int trial = 0; trial < 1500; trial++)
        {
            int m = rng.Next(0, 7);
            int k = rng.Next(1, 10);
            int minSeq = rng.Next(1, 5);
            var raw = new List<string>(m);
            for (int i = 0; i < m; i++)
                raw.Add(RandomDna(rng, rng.Next(0, 30)));
            var seqs = raw.Select(s => new DnaSequence(s)).ToArray();

            var shared = MotifFinder.FindSharedMotifs(seqs, k: k, minSequences: minSeq).ToList();

            AssertWellFormed(shared, raw, k, minSeq);
            AssertMatchesOracle(shared, raw, k, minSeq);
        }
    }

    // Exact matching (INV-05): a 1-substitution variant is a DIFFERENT word and must not
    // be conflated. seq0 has ACGTAC, seq1 has ACATAC (pos2 substitution) — at k=6 they
    // share NO word; the near-miss must not be reported as shared.
    [Test]
    public void FindSharedMotifs_OneSubstitutionVariant_NotTreatedAsShared()
    {
        var raw = new[] { "ACGTAC", "ACATAC" };
        var seqs = raw.Select(s => new DnaSequence(s)).ToArray();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 6, minSequences: 2).ToList();
        shared.Should().BeEmpty("a 1-substitution variant is a different word — matching is exact (INV-05)");
        AssertMatchesOracle(shared, raw, 6, 2);
    }

    // Word repeated within ONE sequence still contributes 1 (matching-sequence semantics,
    // INV-02): seq0 = ATGATGATG (ATG ×3), seq1 = TTTATGTTT (ATG ×1) → ATG indices [0,1].
    [Test]
    public void FindSharedMotifs_RepeatWithinSequence_CountsOncePerSequence()
    {
        var raw = new[] { "ATGATGATG", "TTTATGTTT" };
        var seqs = raw.Select(s => new DnaSequence(s)).ToArray();

        var shared = MotifFinder.FindSharedMotifs(seqs, k: 3, minSequences: 2).ToList();

        var atg = shared.Single(s => s.Sequence == "ATG");
        atg.SequenceIndices.OrderBy(i => i).Should().Equal(new[] { 0, 1 },
            "ATG repeats in seq0 but contributes a single matching-sequence count (INV-02)");
        atg.Prevalence.Should().BeApproximately(1.0, 1e-12, "present in both inputs");
        AssertWellFormed(shared, raw, 3, 2);
    }

    #endregion

    #endregion
}
