using System.Text;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Analysis area — Sequence Similarity (GENOMIC-SIMILARITY-001), the exact
/// k-mer Jaccard similarity measure
/// <see cref="GenomicAnalyzer.CalculateSimilarity(DnaSequence, DnaSequence, int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and extreme sequence pairs / k values to the unit and
/// asserts the code NEVER fails in an undisciplined way: no DivideByZero when the k-mer union is
/// empty (both sequences empty or both shorter than k — the documented empty-union case that must
/// return 0.0, §5.4/§6.1), no IndexOutOfRange when the two operands have very different lengths
/// (the metric is set-based, not positional — it must never index the shorter operand by the
/// longer's positions), no NaN/Infinity, no result outside the documented [0, 100] range
/// (INV-01), and no asymmetry (INV-04). Every input must resolve to a well-defined, theory-correct
/// Jaccard percentage; a raw runtime exception, a NaN, a value &lt; 0 or &gt; 100, or an asymmetric
/// result is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: GENOMIC-SIMILARITY-001 — k-mer Jaccard similarity (alignment-free, exact)
/// Checklist: docs/checklists/03_FUZZING.md, row 179.
/// Algorithm doc: docs/algorithms/Analysis/Sequence_Similarity.md
/// Distinct from row 175 GENOMIC-COMMON-001 (longest common SUBSTRING, positional) — THIS unit is
/// the composition-only, set-resemblance Jaccard index over distinct k-mer SETS.
///
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the row:
///       – identical: sim(A,A) = 100.0 exactly for any non-empty sequence with length ≥ k (A = B
///         ⇒ A∩B = A∪B ⇒ J = 1, INV-02); the canonical self-similarity maximum.
///       – disjoint: sequences whose k-mer SETS share nothing → 0.0 (A∩B = ∅, A∪B ≠ ∅ ⇒ J = 0,
///         INV-03), with NO DivideByZero — and the special empty-union boundary (both empty / both
///         shorter than k) which is the ONLY place the union can be 0; that returns 0.0 (§5.4).
///       – different lengths: sequences of very unequal length → the metric normalizes by the SET
///         UNION (|A|+|B|−|A∩B|), never by a positional alignment, so no IndexOutOfRange and no
///         DivideByZero; the result is still a valid [0,100] percentage and stays symmetric.
///       – k boundaries: k = 1 (every base is a k-mer); k larger than both lengths (both k-mer
///         sets empty → empty union → 0.0); k &lt; 1 → ArgumentOutOfRangeException (§3.3); null
///         operand → ArgumentNullException (§3.3).
/// — docs/checklists/03_FUZZING.md §Description (BE = граничні значення: 0, -1, MaxInt, empty).
///
/// Note on Malformed Content / Injection: each operand is a <see cref="DnaSequence"/>, which is
/// uppercased and validated to the {A,C,G,T} alphabet at construction, so out-of-domain residues,
/// null bytes and unicode cannot reach this method; this is therefore a pure boundary (BE) row over
/// the pair shape (identical / disjoint / different-length / empty) and the integer k, exactly as
/// the checklist row specifies.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Sequence_Similarity.md §2, §3, §5.4, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   result = J(A,B) × 100 where J(A,B) = |A∩B| / |A∪B| over the distinct-k-mer sets A, B;
///   INV-01 0 ≤ result ≤ 100;
///   INV-02 identical non-empty (len ≥ k) sequences → 100.0;
///   INV-03 disjoint k-mer sets → 0.0;
///   INV-04 symmetric: result(a,b,k) = result(b,a,k);
///   INV-05 k-mers compared as SETS — within-sequence repeats counted once;
///   §5.4   empty union (both empty / both shorter than k) → 0.0 (no DivideByZero).
///   GenomicAnalyzer.CalculateSimilarity(DnaSequence, DnaSequence, int kmerSize = 5) → double
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class GenomicSimilarityFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };
    private const double Tolerance = 1e-9;

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
    /// Independent oracle for the k-mer Jaccard percentage, built straight from the spec
    /// (§2.2): the distinct-k-mer SET of each operand, J = |A∩B|/|A∪B|, ×100, with the empty-union
    /// case (both sets empty) returning 0.0 (§5.4). Deliberately a separate implementation from the
    /// unit so the test does not merely re-state the code.
    /// </summary>
    private static double OracleSimilarity(string a, string b, int k)
    {
        HashSet<string> setA = Kmers(a, k);
        HashSet<string> setB = Kmers(b, k);
        int intersection = setA.Count(setB.Contains);
        int union = setA.Count + setB.Count - intersection;
        return union == 0 ? 0.0 : (double)intersection / union * 100.0;
    }

    private static HashSet<string> Kmers(string s, int k)
    {
        var set = new HashSet<string>();
        for (int i = 0; i + k <= s.Length; i++)
            set.Add(s.Substring(i, k));
        return set;
    }

    /// <summary>
    /// Asserts a CalculateSimilarity result is WELL-FORMED per the documented contract:
    ///   • finite (no NaN / Infinity);
    ///   • within the documented [0, 100] range (INV-01);
    ///   • symmetric: result(a,b,k) == result(b,a,k) (INV-04);
    ///   • equal to the independent oracle's Jaccard percentage.
    /// </summary>
    private static void AssertWellFormed(string a, string b, int k)
    {
        double ab = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(b), k);
        double ba = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(b), new DnaSequence(a), k);

        double.IsNaN(ab).Should().BeFalse("a Jaccard percentage is finite (no DivideByZero / 0/0)");
        double.IsInfinity(ab).Should().BeFalse("a Jaccard percentage is finite");
        ab.Should().BeInRange(0.0, 100.0, "0 ≤ J ≤ 1 scaled ×100 (INV-01)");

        ba.Should().BeApproximately(ab, Tolerance, "Jaccard is symmetric: ∩ and ∪ are commutative (INV-04)");
        ab.Should().BeApproximately(OracleSimilarity(a, b, k), Tolerance,
            "result equals the independent Jaccard oracle (J×100)");
    }

    #endregion

    #region GENOMIC-SIMILARITY-001 — k-mer Jaccard similarity (BE: identical, disjoint, different lengths)

    #region Positive sanity — hand-computed documented Jaccard percentage

    // Documented worked example (§7.1): ACGTACGT vs ACGTACGA, k=3 → 80.0.
    //   A = {ACG,CGT,GTA,TAC} (4); B = {ACG,CGT,GTA,TAC,CGA} (5); ∩ = 4; ∪ = 5; J = 4/5 = 0.8 → 80.0.
    [Test]
    public void CalculateSimilarity_DocumentedWorkedExample_Returns80()
    {
        double pct = GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGA"), kmerSize: 3);

        pct.Should().BeApproximately(80.0, Tolerance, "documented walk-through (§7.1): J = 4/5 → 80.0");
        AssertWellFormed("ACGTACGT", "ACGTACGA", 3);
    }

    // A second hand-computed case: ACGT vs ACGT shares one extra path — use k=2 over partly shared.
    //   ACGTT vs ACGAA, k=2: A = {AC,CG,GT,TT} (4); B = {AC,CG,GA,AA} (4); ∩ = {AC,CG} = 2;
    //   ∪ = 4+4-2 = 6; J = 2/6 = 1/3 → 33.333…%.
    [Test]
    public void CalculateSimilarity_HandComputedPartialOverlap_MatchesJaccard()
    {
        double pct = GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence("ACGTT"), new DnaSequence("ACGAA"), kmerSize: 2);

        pct.Should().BeApproximately(100.0 / 3.0, Tolerance, "∩={AC,CG}=2, ∪=6, J=1/3 → 33.33…%");
        AssertWellFormed("ACGTT", "ACGAA", 2);
    }

    // INV-05: within-sequence repeats are counted ONCE (distinct-set semantics). AAAA vs AA, k=2:
    //   both k-mer sets = {AA}; ∩ = ∪ = 1 ⇒ J = 1 ⇒ 100.0, despite very different lengths.
    [Test]
    public void CalculateSimilarity_RepeatedKmersCountedOnce_HomopolymersAreIdentical()
    {
        double pct = GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence("AAAAAAAA"), new DnaSequence("AA"), kmerSize: 2);

        pct.Should().BeApproximately(100.0, Tolerance,
            "both reduce to the distinct set {AA} ⇒ J = 1 (INV-05): repeats counted once");
        AssertWellFormed("AAAAAAAA", "AA", 2);
    }

    #endregion

    #region BE — Boundary: identical (self-similarity = 100.0, the maximum)

    // sim(A,A) = 100.0 exactly for a non-empty sequence with length ≥ k (INV-02).
    [Test]
    public void CalculateSimilarity_IdenticalNonEmpty_Returns100()
    {
        var seq = new DnaSequence("ACGTACGTTACG");
        double pct = GenomicAnalyzer.CalculateSimilarity(seq, seq, kmerSize: 5);

        pct.Should().BeApproximately(100.0, Tolerance, "A = B ⇒ A∩B = A∪B ⇒ J = 1 (INV-02)");
    }

    // Fuzz: identical random pairs (length ≥ k) → exactly 100.0; the canonical maximum self-similarity.
    [Test]
    [CancelAfter(30_000)]
    public void CalculateSimilarity_IdenticalRandom_Always100_WhenLongerThanK()
    {
        var rng = new Random(179_001);
        for (int trial = 0; trial < 800; trial++)
        {
            int k = rng.Next(1, 6);
            string s = RandomDna(rng, rng.Next(k, 40)); // guarantee length ≥ k ⇒ non-empty k-mer set
            var seq = new DnaSequence(s);

            double pct = GenomicAnalyzer.CalculateSimilarity(seq, seq, k);

            pct.Should().BeApproximately(100.0, Tolerance,
                "identical sequences of length ≥ k are maximally self-similar (INV-02)");
            AssertWellFormed(s, s, k);
        }
    }

    // Self-similarity is the maximum: sim(A,A) ≥ sim(A,B) for any B (a corollary of J ≤ 1).
    [Test]
    [CancelAfter(30_000)]
    public void CalculateSimilarity_SelfSimilarityIsMaximum()
    {
        var rng = new Random(179_002);
        for (int trial = 0; trial < 600; trial++)
        {
            int k = rng.Next(1, 5);
            string a = RandomDna(rng, rng.Next(k, 30));
            string b = RandomDna(rng, rng.Next(k, 30));

            double self = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(a), k);
            double cross = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(b), k);

            self.Should().BeApproximately(100.0, Tolerance, "self-similarity is the maximum (INV-02)");
            cross.Should().BeLessThanOrEqualTo(self + Tolerance, "no pair exceeds self-similarity (INV-01)");
        }
    }

    #endregion

    #region BE — Boundary: disjoint (no shared k-mers → 0.0, no DivideByZero)

    // Disjoint k-mer sets (all-A vs all-C, k=1) → 0.0, union non-empty, no DivideByZero (INV-03).
    [Test]
    public void CalculateSimilarity_DisjointHomopolymers_Returns0()
    {
        double pct = GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence("AAAAAAAA"), new DnaSequence("CCCCCCCC"), kmerSize: 1);

        pct.Should().BeApproximately(0.0, Tolerance, "A∩B = ∅, A∪B = {A,C} ≠ ∅ ⇒ J = 0 (INV-03)");
        AssertWellFormed("AAAAAAAA", "CCCCCCCC", 1);
    }

    // Fuzz: operands drawn from disjoint sub-alphabets ({A,C} vs {G,T}) → no shared k-mer at any k
    // ≤ 1 boundary; at k=1 the alphabets are disjoint ⇒ similarity exactly 0.0, never a DivideByZero.
    [Test]
    [CancelAfter(30_000)]
    public void CalculateSimilarity_DisjointAlphabets_Returns0_NeverThrows()
    {
        var rng = new Random(179_003);
        char[] left = { 'A', 'C' };
        char[] right = { 'G', 'T' };
        for (int trial = 0; trial < 600; trial++)
        {
            var sa = new StringBuilder();
            int la = rng.Next(1, 30);
            for (int i = 0; i < la; i++) sa.Append(left[rng.Next(left.Length)]);
            var sb = new StringBuilder();
            int lb = rng.Next(1, 30);
            for (int i = 0; i < lb; i++) sb.Append(right[rng.Next(right.Length)]);
            string a = sa.ToString(), b = sb.ToString();

            double pct = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(b), kmerSize: 1);

            pct.Should().BeApproximately(0.0, Tolerance,
                "disjoint single-base alphabets share no k-mer ⇒ J = 0 (INV-03)");
            AssertWellFormed(a, b, 1);
        }
    }

    #endregion

    #region BE — Boundary: empty union (both empty / both shorter than k → 0.0, the only union==0 path)

    // Both sequences empty → empty union → 0.0, NO DivideByZero (§5.4).
    [Test]
    public void CalculateSimilarity_BothEmpty_Returns0_NoDivideByZero()
    {
        Action act = () => GenomicAnalyzer.CalculateSimilarity(new DnaSequence(""), new DnaSequence(""), kmerSize: 5);
        act.Should().NotThrow("empty union is a documented boundary, not an error (§5.4)");

        GenomicAnalyzer.CalculateSimilarity(new DnaSequence(""), new DnaSequence(""), kmerSize: 5)
            .Should().BeApproximately(0.0, Tolerance, "empty union ⇒ undefined → impl returns 0 (§5.4)");
    }

    // One empty, the other non-empty → empty intersection over non-empty union → 0.0 (§6.1).
    [Test]
    public void CalculateSimilarity_OneEmpty_Returns0()
    {
        double pct = GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence(""), new DnaSequence("ACGTACGT"), kmerSize: 3);

        pct.Should().BeApproximately(0.0, Tolerance, "empty ∩ over non-empty ∪ ⇒ J = 0 (§6.1)");
        AssertWellFormed("", "ACGTACGT", 3);
    }

    // Both shorter than k → both k-mer sets empty → empty union → 0.0, no DivideByZero (§6.1).
    [Test]
    public void CalculateSimilarity_BothShorterThanK_Returns0_NoDivideByZero()
    {
        Action act = () => GenomicAnalyzer.CalculateSimilarity(new DnaSequence("AC"), new DnaSequence("GT"), kmerSize: 5);
        act.Should().NotThrow("both-shorter-than-k yields an empty union, a documented boundary (§6.1)");

        GenomicAnalyzer.CalculateSimilarity(new DnaSequence("AC"), new DnaSequence("GT"), kmerSize: 5)
            .Should().BeApproximately(0.0, Tolerance, "both k-mer sets empty ⇒ empty union → 0 (§6.1)");
    }

    // Fuzz: random pairs both shorter than k → always 0.0, never a DivideByZero.
    [Test]
    [CancelAfter(30_000)]
    public void CalculateSimilarity_RandomBothShorterThanK_Always0()
    {
        var rng = new Random(179_004);
        for (int trial = 0; trial < 600; trial++)
        {
            int k = rng.Next(3, 12);
            string a = RandomDna(rng, rng.Next(0, k));   // length < k
            string b = RandomDna(rng, rng.Next(0, k));   // length < k
            var sa = new DnaSequence(a);
            var sb = new DnaSequence(b);

            double pct = GenomicAnalyzer.CalculateSimilarity(sa, sb, k);

            pct.Should().BeApproximately(0.0, Tolerance,
                "both sequences shorter than k ⇒ empty union → 0 (§5.4), no DivideByZero");
            AssertWellFormed(a, b, k);
        }
    }

    #endregion

    #region BE — Boundary: different lengths (set-normalized, no IndexOutOfRange, stays in range/symmetric)

    // Very unequal lengths: a 1-base sequence vs a long one, k=1. The metric is set-based, so it
    // must NOT index the short operand by the long operand's positions (no IndexOutOfRange).
    //   "A" vs "ACGTACGTACGT", k=1: A = {A}, B = {A,C,G,T}; ∩ = {A} = 1; ∪ = 4; J = 1/4 → 25.0.
    [Test]
    public void CalculateSimilarity_VeryDifferentLengths_SetNormalized_NoIndexOutOfRange()
    {
        double pct = GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence("A"), new DnaSequence("ACGTACGTACGT"), kmerSize: 1);

        pct.Should().BeApproximately(25.0, Tolerance, "{A} ∩ {A,C,G,T} = 1; ∪ = 4; J = 1/4 → 25.0");
        AssertWellFormed("A", "ACGTACGTACGT", 1);
    }

    // Subset relationship across different lengths: the SHORT sequence's k-mers are a subset of the
    // LONG one's. ACG (k-mers {ACG} at k=3) vs ACGACGACG (k-mers {ACG,CGA,GAC} at k=3):
    //   ∩ = {ACG} = 1; ∪ = 3; J = 1/3 → 33.33…%. Normalization is by the union, never positional.
    [Test]
    public void CalculateSimilarity_ShortSubsetOfLong_NormalizesByUnion()
    {
        double pct = GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence("ACG"), new DnaSequence("ACGACGACG"), kmerSize: 3);

        pct.Should().BeApproximately(100.0 / 3.0, Tolerance, "∩={ACG}=1, ∪=3, J=1/3 → 33.33…%");
        AssertWellFormed("ACG", "ACGACGACG", 3);
    }

    // Fuzz: deliberately very unequal lengths (1 vs up to 200) over varied k → never throws, finite,
    // [0,100], symmetric, matches the oracle. Catches IndexOutOfRange / DivideByZero on length skew.
    [Test]
    [CancelAfter(60_000)]
    public void CalculateSimilarity_HighlyUnequalLengths_NeverThrows_WellFormed()
    {
        var rng = new Random(179_005);
        for (int trial = 0; trial < 1200; trial++)
        {
            int k = rng.Next(1, 8);
            string shortS = RandomDna(rng, rng.Next(0, 4));
            string longS = RandomDna(rng, rng.Next(50, 200));

            Action act = () => GenomicAnalyzer.CalculateSimilarity(new DnaSequence(shortS), new DnaSequence(longS), k);
            act.Should().NotThrow("unequal lengths are a set-based boundary, not an error (BE: different lengths)");

            AssertWellFormed(shortS, longS, k);
        }
    }

    #endregion

    #region BE — Boundary: k boundaries (k=1; k > both lengths; k<1 and null → documented exceptions)

    // k larger than both sequence lengths → both k-mer sets empty → 0.0, no crash.
    [Test]
    public void CalculateSimilarity_KLargerThanBoth_Returns0()
    {
        double pct = GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence("ACGT"), new DnaSequence("ACGT"), kmerSize: 99);

        pct.Should().BeApproximately(0.0, Tolerance,
            "k > both lengths ⇒ both k-mer sets empty ⇒ empty union → 0 (§5.4) even though sequences are identical");
    }

    // k = 1: every base is a k-mer; reduces to base-composition Jaccard. ACGT vs ACGT → 100.0.
    [Test]
    public void CalculateSimilarity_KEqualsOne_BaseCompositionJaccard()
    {
        GenomicAnalyzer.CalculateSimilarity(new DnaSequence("ACGT"), new DnaSequence("ACGT"), kmerSize: 1)
            .Should().BeApproximately(100.0, Tolerance, "identical base composition ⇒ J = 1 at k=1");

        // ACGT vs ACG (different length): {A,C,G,T} vs {A,C,G}; ∩=3; ∪=4; J=3/4 → 75.0.
        GenomicAnalyzer.CalculateSimilarity(new DnaSequence("ACGT"), new DnaSequence("ACG"), kmerSize: 1)
            .Should().BeApproximately(75.0, Tolerance, "{A,C,G,T}∩{A,C,G}=3, ∪=4, J=3/4 → 75.0");
    }

    // kmerSize < 1 → ArgumentOutOfRangeException (§3.3). BE: 0, −1, int.MinValue.
    [Test]
    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void CalculateSimilarity_KBelowOne_Throws(int k)
    {
        Action act = () => GenomicAnalyzer.CalculateSimilarity(new DnaSequence("ACGT"), new DnaSequence("ACGT"), k);
        act.Should().Throw<ArgumentOutOfRangeException>("k-mer size must be ≥ 1 (§3.3)");
    }

    // int.MaxValue k → no k-mer fits → empty union → 0.0, no overflow/crash.
    [Test]
    [CancelAfter(15_000)]
    public void CalculateSimilarity_KMaxInt_Returns0_NoOverflow()
    {
        Action act = () => GenomicAnalyzer.CalculateSimilarity(
            new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT"), int.MaxValue);
        act.Should().NotThrow("an unreachable k is a boundary, not an error (BE: MaxInt)");

        GenomicAnalyzer.CalculateSimilarity(new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT"), int.MaxValue)
            .Should().BeApproximately(0.0, Tolerance, "no k-mer of length int.MaxValue exists ⇒ empty union → 0");
    }

    // Null operands → ArgumentNullException (§3.3).
    [Test]
    public void CalculateSimilarity_NullOperands_Throw()
    {
        ((Action)(() => GenomicAnalyzer.CalculateSimilarity(null!, new DnaSequence("ACGT"), 3)))
            .Should().Throw<ArgumentNullException>("sequence1 null ⇒ ArgumentNullException (§3.3)");
        ((Action)(() => GenomicAnalyzer.CalculateSimilarity(new DnaSequence("ACGT"), null!, 3)))
            .Should().Throw<ArgumentNullException>("sequence2 null ⇒ ArgumentNullException (§3.3)");
    }

    #endregion

    #region BE — Broad fuzz: random pairs / k never crash, match the documented contract

    // Random pairs and random k over varied lengths: never throws/hangs, finite, [0,100], symmetric,
    // equal to the independent Jaccard oracle.
    [Test]
    [CancelAfter(60_000)]
    public void CalculateSimilarity_RandomPairsAndK_NeverThrows_WellFormed()
    {
        var rng = new Random(179_006);
        for (int trial = 0; trial < 2000; trial++)
        {
            int k = rng.Next(1, 10);
            string a = RandomDna(rng, rng.Next(0, 60));
            string b = RandomDna(rng, rng.Next(0, 60));

            Action act = () => GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(b), k);
            act.Should().NotThrow("any ACGT pair with k ≥ 1 is a valid input");

            AssertWellFormed(a, b, k);
        }
    }

    // Default k (5): exercise the production default path over random pairs — never throws, well-formed.
    [Test]
    [CancelAfter(30_000)]
    public void CalculateSimilarity_DefaultK_RandomPairs_WellFormed()
    {
        var rng = new Random(179_007);
        for (int trial = 0; trial < 800; trial++)
        {
            string a = RandomDna(rng, rng.Next(0, 50));
            string b = RandomDna(rng, rng.Next(0, 50));

            double ab = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(a), new DnaSequence(b));
            double ba = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(b), new DnaSequence(a));

            double.IsNaN(ab).Should().BeFalse();
            ab.Should().BeInRange(0.0, 100.0, "default-k result in documented range (INV-01)");
            ba.Should().BeApproximately(ab, Tolerance, "symmetric at default k (INV-04)");
            ab.Should().BeApproximately(OracleSimilarity(a, b, 5), Tolerance, "default k = 5 matches oracle");
        }
    }

    #endregion

    #endregion
}
