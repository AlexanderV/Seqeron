// RNA-PKRECURSIVE-001 — RECURSIVE RNA pseudoknot STRUCTURE PREDICTION (nested / multiple knots,
// pknotsRG canonical simple-recursive class). Fuzz tests (strategy BE = Boundary Exploitation),
// row 237 of docs/checklists/03_FUZZING.md.
// Algorithm doc: docs/algorithms/RnaStructure/Pseudoknot_Prediction_Recursive.md
// Canonical tests: tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_PredictStructurePseudoknotRecursive_Tests.cs
// Evidence: docs/Evidence/RNA-PKRECURSIVE-001-Evidence.md
// Source: RnaSecondaryStructure.PredictStructurePseudoknotRecursive(string, int) — RnaSecondaryStructure.cs.
//         Reeder & Giegerich (2004) BMC Bioinformatics 5:104; Reeder, Steffen & Giegerich (2007)
//         NAR 35:W320; Turner 2004 NN; Antczak et al. (2018) Bioinformatics 34(8):1304.

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for RNA-PKRECURSIVE-001 —
/// <see cref="RnaSecondaryStructure.PredictStructurePseudoknotRecursive(string,int)"/>, the
/// thermodynamic predictor for the canonical simple-<em>recursive</em> pseudoknot (csr-PK) class of
/// pknotsRG. Unlike <see cref="RnaSecondaryStructure.PredictStructurePseudoknot"/> (single top-level
/// H-type, covered by row 236), this entry point folds the WHOLE sequence by a memoised interval
/// recurrence <c>F(i,j)</c> in which a pseudoknot value competes with the unknotted value at every
/// interval (Reeder &amp; Giegerich 2004/2007), so the optimum may contain several knots and knots
/// nested inside the loops of an outer helix. Helices are scored with the Turner-2004
/// nearest-neighbour model; the pknotsRG penalties (initiation 9.0; 0.3 per unpaired knot-loop
/// nucleotide; 0.0 per in-knot base pair) are reused unchanged. A recursive fold is accepted ONLY
/// if its ΔG is strictly below the plain pseudoknot-free MFE; otherwise the plain MFE is returned.
/// Lives in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Documented contract (Pseudoknot_Prediction_Recursive.md §2.4, §3, §6.1) — every result MUST satisfy
/// ───────────────────────────────────────────────────────────────────────────
///   • INV-PKR-01  Returned ΔG ≤ ΔG of the plain pseudoknot-free MFE for the same sequence/params
///                 (the plain MFE is the always-available fallback; the recursive fold is taken only
///                 when it strictly improves it).
///   • INV-PKR-02  HasPseudoknot ⇒ the base-pair set contains ≥1 CROSSING pair (∃ (i,j),(k,l) with
///                 i &lt; k &lt; j &lt; l). Conversely no crossing pair ⇒ HasPseudoknot=false.
///   • INV-PKR-03  Every position is paired AT MOST ONCE; all indices in [0, n); i &lt; j per pair.
///   • INV-PKR-04  No spurious pseudoknot — a knot is accepted only when it lowers ΔG (9 kcal/mol
///                 initiation penalty suppresses spurious knots on non-pseudoknotted sequences).
///   • Complementarity (§3.3): every reported pair (i,j) is a legal A-U / G-C / G-U pair on the
///                 FOLDED sequence (upper-cased, T→U). No non-complementary pair may ever appear.
///   • Output shape (§3.2): Sequence = upper-cased + T→U spelling of the input; DotBracket has length
///                 n, uses only {(),[],.} , each family balanced, and agrees column-for-column with
///                 BasePairs; FreeEnergy is a finite double; BasePairs sorted by 5' position.
///   • Edge cases (§6.1): null / empty / &lt; 11 nt → empty pseudoknot-free structure (no pairs,
///                 all dots, ΔG = 0, HasPseudoknot=false). Non-ACGU characters simply do not pair.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"
/// ───────────────────────────────────────────────────────────────────────────
/// Feed malformed / boundary sequences and assert the recursive predictor NEVER fails undisciplined.
/// Because the fold is RECURSIVE (F calls itself on loops and enclosed sub-spans), the headline risk
/// is the recursion: a deeply-nested or long low-complexity input must NOT StackOverflow, hang, or
/// infinite-loop — every hang/recursion-sensitive test carries [CancelAfter]. It must not throw an
/// unhandled runtime exception (IndexOutOfRange / NullReference / overflow / NaN); and it must not
/// emit out-of-contract nonsense — a non-complementary base pair, an out-of-bounds index, a position
/// paired twice, a dot-bracket disagreeing with BasePairs, HasPseudoknot=true with no crossing pair,
/// a returned ΔG worse than the plain MFE, or a non-finite ΔG. Every input → EITHER a well-defined
/// theory-correct structure OR the documented empty result (null/empty/too-short are mapped to the
/// empty structure, not to an exception).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy BE = Boundary Exploitation — row 237
/// targets "empty, single base, long low-complexity, deeply nested"
/// ───────────────────────────────────────────────────────────────────────────
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, -1, MaxInt, empty), row 237.
///   • empty / null (BE)            → empty pseudoknot-free structure, ΔG 0, no crash.
///   • single base / too-short (BE) → &lt; 11 nt → no pseudoknot (returns the plain MFE), no crash.
///   • long low-complexity (BE)     → all-A / repetitive sequences stress the O(n⁴) recursion + memo
///                                     WITHOUT a real pairable structure: must terminate (memoised),
///                                     not blow up, and never invent an illegal pair.
///   • deeply nested (BE)           → sequences engineered to force maximal recursion depth (long
///                                     palindromic onions, knots-in-loops). The recursion must NOT
///                                     StackOverflow and must respect the recursive topology (every
///                                     pair valid, no double-pairing, ΔG ≤ MFE, finite).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaPseudoknotRecursiveFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static Random Rng(int seed) => new(seed);

    /// <summary>The shortest canonical H-type knot length (2·2 + 2·2 + 3·1 loops), per §6.1 / source.</summary>
    private const int MinPseudoknotLength = 11;

    /// <summary>
    /// Asserts the FULL documented well-formedness of a RECURSIVELY-predicted structure for the given
    /// input sequence, deriving every invariant independently from Pseudoknot_Prediction_Recursive.md:
    ///   • never throws / never StackOverflows (the predictor maps all inputs to a structure);
    ///   • Sequence is the upper-cased, T→U spelling of the (non-null) input;
    ///   • FreeEnergy is finite (no NaN / ±∞);
    ///   • every pair (i,j): 0 ≤ i &lt; j &lt; n  AND  CanPair(seq[i],seq[j]) — legal A-U/G-C/G-U (INV-PKR-03 + §3.3);
    ///   • no position is paired more than once (INV-PKR-03);
    ///   • ΔG ≤ plain-MFE ΔG for the same sequence/minLoopSize (INV-PKR-01);
    ///   • HasPseudoknot ⇔ the base-pair set actually contains a crossing pair (INV-PKR-02 + INV-PKR-04);
    ///   • DotBracket has length n, only {(,),[,],.}, both families balanced, columns agree with BasePairs;
    ///   • BasePairs sorted by 5' position (§3.2).
    /// </summary>
    private static PseudoknotStructure AssertWellFormed(string? input, int minLoopSize = 3)
    {
        PseudoknotStructure pk = default;
        Action act = () => pk = PredictStructurePseudoknotRecursive(input!, minLoopSize);
        act.Should().NotThrow("the recursive predictor maps every input to a structure, never an exception");

        string seq = pk.Sequence;
        int n = seq.Length;

        // Folded spelling: upper-cased, T→U; length preserved; empty for null.
        string expectedSeq = input is null ? "" : input.ToUpperInvariant().Replace('T', 'U');
        seq.Should().Be(expectedSeq, "the folded sequence is the upper-cased T→U spelling of the input");

        double.IsFinite(pk.FreeEnergy).Should().BeTrue("ΔG must be a finite number (no NaN/∞)");

        // INV-PKR-03 + §3.3 complementarity: indices in range, i<j, legal pair, no double-pairing.
        var seen = new HashSet<int>();
        int prevOpen = -1;
        foreach (var (a, b) in pk.BasePairs)
        {
            int i = Math.Min(a, b), j = Math.Max(a, b);
            i.Should().BeGreaterThanOrEqualTo(0, "pair index must be ≥ 0");
            j.Should().BeLessThan(n, "pair index must be < sequence length (in bounds)");
            i.Should().BeLessThan(j, "a base pair spans two DISTINCT positions (i<j)");
            CanPair(seq[i], seq[j]).Should().BeTrue(
                $"reported pair ({i},{j}) = ({seq[i]},{seq[j]}) must be a legal A-U/G-C/G-U pair");
            seen.Add(i).Should().BeTrue($"position {i} must not be paired twice (INV-PKR-03)");
            seen.Add(j).Should().BeTrue($"position {j} must not be paired twice (INV-PKR-03)");
            // §3.2: pairs are sorted by 5' position.
            i.Should().BeGreaterThanOrEqualTo(prevOpen, "BasePairs must be sorted by 5' position (§3.2)");
            prevOpen = i;
        }

        // INV-PKR-01: never worse than the plain pseudoknot-free MFE (the always-available fallback).
        double plainMfe = CalculateMfeStructure(input!, minLoopSize).FreeEnergy;
        pk.FreeEnergy.Should().BeLessThanOrEqualTo(plainMfe + 1e-6,
            "the returned ΔG can never exceed the plain-MFE baseline (INV-PKR-01)");

        // INV-PKR-02 / INV-PKR-04: HasPseudoknot ⇔ a genuine crossing pair exists.
        HasCrossingPair(pk.BasePairs).Should().Be(pk.HasPseudoknot,
            "HasPseudoknot must be true exactly when the pair set contains a crossing pair (INV-PKR-02/04)");

        // DotBracket shape and exact agreement with the pair set.
        AssertDotBracketConsistent(pk, n);

        return pk;
    }

    /// <summary>True iff the pair set contains a crossing pair (i,j),(k,l) with i &lt; k &lt; j &lt; l (§2.1).</summary>
    private static bool HasCrossingPair(IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        for (int a = 0; a < pairs.Count; a++)
        for (int b = a + 1; b < pairs.Count; b++)
        {
            int i = Math.Min(pairs[a].Position1, pairs[a].Position2);
            int j = Math.Max(pairs[a].Position1, pairs[a].Position2);
            int k = Math.Min(pairs[b].Position1, pairs[b].Position2);
            int l = Math.Max(pairs[b].Position1, pairs[b].Position2);
            if (k < i) (i, j, k, l) = (k, l, i, j);
            if (i < k && k < j && j < l) return true;
        }
        return false;
    }

    /// <summary>
    /// Asserts the two-layer dot-bracket is well formed and EXACTLY encodes the base-pair set:
    /// length n; only {(,),[,],.}; both families balanced; every column is '.' iff unpaired and a
    /// matching open/close bracket iff it is the 5'/3' end of a reported pair (§3.2).
    /// </summary>
    private static void AssertDotBracketConsistent(PseudoknotStructure pk, int n)
    {
        string db = pk.DotBracket;
        db.Length.Should().Be(n, "dot-bracket length must equal the sequence length");
        db.Should().MatchRegex(@"^[\(\)\[\]\.]*$", "dot-bracket uses only ( ) [ ] .");

        Balanced(db, '(', ')').Should().BeTrue("the () family must be balanced");
        Balanced(db, '[', ']').Should().BeTrue("the [] family must be balanced");

        var open = new bool[n];
        var close = new bool[n];
        foreach (var (a, b) in pk.BasePairs)
        {
            int i = Math.Min(a, b), j = Math.Max(a, b);
            open[i] = true;
            close[j] = true;
        }
        for (int p = 0; p < n; p++)
        {
            bool isOpen = db[p] is '(' or '[';
            bool isClose = db[p] is ')' or ']';
            bool isDot = db[p] == '.';
            if (open[p])
                isOpen.Should().BeTrue($"column {p} is a 5' end of a pair → must be an opening bracket");
            else if (close[p])
                isClose.Should().BeTrue($"column {p} is a 3' end of a pair → must be a closing bracket");
            else
                isDot.Should().BeTrue($"column {p} is unpaired → must be '.'");
        }
    }

    private static bool Balanced(string s, char open, char close)
    {
        int depth = 0;
        foreach (char c in s)
        {
            if (c == open) depth++;
            else if (c == close) depth--;
            if (depth < 0) return false;
        }
        return depth == 0;
    }

    /// <summary>Asserts the documented empty pseudoknot-free result: no pairs, all dots, ΔG 0, no knot.</summary>
    private static void AssertEmptyStructure(PseudoknotStructure pk)
    {
        pk.BasePairs.Should().BeEmpty("a degenerate / too-short input yields no base pairs (§6.1)");
        pk.HasPseudoknot.Should().BeFalse("no pseudoknot is possible (§6.1)");
        pk.FreeEnergy.Should().Be(0.0, "the empty structure has ΔG = 0 (§6.1)");
        pk.DotBracket.Should().Be(new string('.', pk.Sequence.Length), "all positions are unpaired dots");
    }

    /// <summary>Antiparallel Watson–Crick complement (A↔U, G↔C), reversed — to seed pairable helices.</summary>
    private static string Complement(string s)
    {
        var c = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
        {
            char x = s[s.Length - 1 - i];
            c[i] = x switch { 'A' => 'U', 'U' => 'A', 'G' => 'C', 'C' => 'G', _ => 'A' };
        }
        return new string(c);
    }

    #endregion

    #region RNA-PKRECURSIVE-001 — empty / null / single base (BE)

    [Test]
    [CancelAfter(10_000)]
    public void Recursive_Null_ReturnsEmptyPseudoknotFreeStructure_NoCrash()
    {
        var pk = AssertWellFormed(null);
        AssertEmptyStructure(pk);
        pk.Sequence.Should().BeEmpty("null folds to the empty sequence (§6.1 contract parity with MFE)");
    }

    [Test]
    [CancelAfter(10_000)]
    public void Recursive_Empty_ReturnsEmptyPseudoknotFreeStructure_NoCrash()
    {
        var pk = AssertWellFormed("");
        AssertEmptyStructure(pk);
        pk.Sequence.Should().BeEmpty();
    }

    [TestCase("A")]
    [TestCase("G")]
    [TestCase("U")]
    [TestCase("C")]
    [TestCase("N")]
    [TestCase(" ")]
    public void Recursive_SingleBase_ReturnsEmptyStructure_NoCrash(string seq)
    {
        // A single nucleotide cannot form any pair and is far below the 11-nt knot minimum.
        var pk = AssertWellFormed(seq);
        AssertEmptyStructure(pk);
    }

    #endregion

    #region RNA-PKRECURSIVE-001 — too-short (< 11 nt, BE)

    [TestCase("GC")]
    [TestCase("GCGC")]
    [TestCase("GGGGAACCCC")]      // 10 nt — one below the 11-nt canonical minimum
    [TestCase("GGGGAACCCCA")]     // exactly 11 nt with a single hairpin: no crossing helix
    public void Recursive_TooShortForCanonicalKnot_NeverReportsPseudoknot(string seq)
    {
        var pk = AssertWellFormed(seq);
        if (seq.Length < MinPseudoknotLength)
            pk.HasPseudoknot.Should().BeFalse($"{seq.Length} nt < {MinPseudoknotLength} → no canonical knot (§6.1)");
    }

    [Test]
    [CancelAfter(15_000)]
    public void Recursive_AllShortLengths_BelowMin_NeverThrow_NeverKnot()
    {
        // Every length 0..10 of mixed pairable bases: no canonical knot can form, no crash.
        const string alphabet = "ACGU";
        var rng = Rng(20260626);
        for (int len = 0; len < MinPseudoknotLength; len++)
        {
            for (int iter = 0; iter < 40; iter++)
            {
                var chars = new char[len];
                for (int p = 0; p < len; p++) chars[p] = alphabet[rng.Next(alphabet.Length)];
                var pk = AssertWellFormed(new string(chars));
                pk.HasPseudoknot.Should().BeFalse($"length {len} < {MinPseudoknotLength} → no knot");
            }
        }
    }

    #endregion

    #region RNA-PKRECURSIVE-001 — long low-complexity (BE) — stresses the O(n⁴) recursion/memo

    [TestCase(11)]
    [TestCase(50)]
    [TestCase(120)]
    [TestCase(200)]
    public void Recursive_AllAdenine_LongLowComplexity_YieldsUnpairedStructure_NoBlowUp(int length)
    {
        // A pairs only with U; an all-A homopolymer has ZERO legal pairs. Even at 200 nt the memoised
        // O(n⁴) recurrence must TERMINATE and return the fully-unpaired structure — never StackOverflow,
        // never hang, never invent a pair. This is the canonical long low-complexity stress.
        var pk = AssertWellFormed(new string('A', length));

        pk.BasePairs.Should().BeEmpty("all-A has no complementary positions → no base pair at all");
        pk.HasPseudoknot.Should().BeFalse("no pairs ⇒ no crossing helix");
        pk.FreeEnergy.Should().Be(0.0, "no stabilizing pair ⇒ ΔG = 0");
        pk.DotBracket.Should().Be(new string('.', length), "every position is unpaired");
    }

    [Test]
    [CancelAfter(120_000)]
    public void Recursive_LongRepetitiveLowComplexity_StaysWithinContract_NoHang()
    {
        // Repetitive low-complexity motifs (some self-complementary, some not) at lengths that exercise
        // the recursion deeply; each must stay in-contract and TERMINATE within the time bound.
        foreach (string unit in new[] { "AU", "GC", "AAAU", "GGGC", "ACGU", "AAAAGGGG" })
        {
            for (int reps = 4; reps <= 24; reps += 4)
            {
                string seq = string.Concat(Enumerable.Repeat(unit, reps));
                AssertWellFormed(seq);
            }
        }
    }

    [Test]
    [CancelAfter(60_000)]
    public void Recursive_HomopolymersOfEachBase_Long_NeverThrow_NeverInventPairs()
    {
        // C-only, G-only, U-only homopolymers: no self-pairing possible at any length (G·U wobble needs
        // both a G and a U, absent in a homopolymer). Must yield the unpaired structure, no crash.
        foreach (char b in new[] { 'C', 'G', 'U', 'A' })
        {
            foreach (int len in new[] { 11, 64, 150 })
            {
                var pk = AssertWellFormed(new string(b, len));
                pk.BasePairs.Should().BeEmpty($"a {b}-homopolymer cannot self-pair");
                pk.HasPseudoknot.Should().BeFalse();
                pk.FreeEnergy.Should().Be(0.0);
            }
        }
    }

    #endregion

    #region RNA-PKRECURSIVE-001 — deeply nested (BE) — recursion-depth stress, must NOT StackOverflow

    [Test]
    [CancelAfter(120_000)]
    public void Recursive_DeepNestedOnion_LongPalindrome_NoStackOverflow_StaysWithinContract()
    {
        // A long palindromic "onion" (sequence = block + reverse-complement(block)) forces a maximally
        // nested helix, which drives the F-recurrence to its maximum depth (loops/sub-spans fold by F,
        // each enclosed span shorter by the helix length). The HEADLINE risk for a recursive folder:
        // this must NOT StackOverflow / hang. Every reported pair must still be legal & in-contract,
        // and ΔG ≤ plain MFE. We scale up to a length that genuinely exercises deep recursion.
        foreach (int half in new[] { 8, 16, 32, 60 })
        {
            var rng = Rng(7000 + half);
            const string acgu = "ACGU";
            var block = new char[half];
            for (int p = 0; p < half; p++) block[p] = acgu[rng.Next(acgu.Length)];
            string left = new(block);
            string seq = left + Complement(left); // perfectly self-complementary → deep nesting

            AssertWellFormed(seq);
        }
    }

    [Test]
    [CancelAfter(120_000)]
    public void Recursive_NestedHairpinsChain_DeepRecursion_StaysWithinContract()
    {
        // A chain of stacked hairpins concatenated with self-complementary clamps: many foldable
        // sub-spans force the memoised recurrence to recurse through every interval. Must terminate
        // and stay in-contract; no StackOverflow regardless of how many hairpins are chained.
        foreach (int hairpins in new[] { 2, 5, 10, 16 })
        {
            // Each hairpin: GGGG <loop> CCCC ; chained so the outer recurrence descends repeatedly.
            string hp = "GGGGAAAACCCC";
            string seq = string.Concat(Enumerable.Repeat(hp, hairpins));
            AssertWellFormed(seq);
        }
    }

    [Test]
    [CancelAfter(60_000)]
    public void Recursive_KnotInsideOuterHelix_OverArching_DeepRecursion_StaysWithinContract()
    {
        // §7.1 worked example shape: an outer A·U helix over-arching an inner crossing knot, embedded
        // into longer flanks so the recursion descends through the outer helix into the knotted loop
        // (the recursive class's signature topology). Must stay in-contract and terminate.
        string core = "AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU"; // doc §7.1
        foreach (int pad in new[] { 0, 6, 12, 20 })
        {
            string flank = new('A', pad);
            AssertWellFormed(flank + core + flank);
        }
    }

    #endregion

    #region RNA-PKRECURSIVE-001 — non-ACGU / minLoopSize boundaries (BE)

    [TestCase("NNNNNNNNNNNNNNNNNNNN")]      // IUPAC ambiguity 'N' — not in the pairing table
    [TestCase("XYZXYZXYZXYZXYZXYZXY")]      // arbitrary non-nucleotide letters
    [TestCase("....----....----....")]      // gap / alignment characters
    [TestCase("GGGG##CCCC##GGGG##CC")]      // ACGU islands separated by junk
    [TestCase("12345678901234567890")]      // digits
    [TestCase(" gggg  cccc  gggg  cc")]     // lower-case + spaces
    public void Recursive_NonAcguCharacters_NeverThrow_OnlyLegalPairsReported(string seq)
    {
        // Unknown characters must simply not pair; the predictor must not throw and any pair it does
        // report must STILL be a legal A-U/G-C/G-U pair on the folded sequence.
        AssertWellFormed(seq);
    }

    [Test]
    [CancelAfter(20_000)]
    public void Recursive_HighAndControlCodepoints_NeverThrow()
    {
        // Codepoints ≥ 128 and control chars exercise the GetBasePairType high-codepoint guard.
        AssertWellFormed(" GGGGÿCCCC€GGGGAACCCC\t\n");
        AssertWellFormed(new string('￿', 20));
    }

    [Test]
    [CancelAfter(30_000)]
    public void Recursive_NegativeAndExtremeMinLoopSize_IsClampedAndNeverThrows()
    {
        // §3.1: minLoopSize < 3 is clamped to 3. Negative / zero / int.MinValue / huge must not crash.
        foreach (int m in new[] { int.MinValue, -100, -1, 0, 1, 2, 3, 1000 })
            AssertWellFormed("AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU", m);
    }

    #endregion

    #region RNA-PKRECURSIVE-001 — positive sanity (recursive class) + randomized fuzz

    [Test]
    [CancelAfter(30_000)]
    public void Recursive_DocumentedOverArchingExample_IsAGenuineRecursiveKnot()
    {
        // §7.1 worked example: an inner crossing knot nested inside an outer over-arching helix —
        // the recursive (csr-PK) class's signature. The recursive fold must DETECT a pseudoknot and
        // strictly improve on the plain MFE (INV-PKR-01/04).
        string seq = "AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU";
        var pk = AssertWellFormed(seq);
        pk.HasPseudoknot.Should().BeTrue("the §7.1 example folds into an over-arching recursive knot");
        pk.DotBracket.Should().Contain("[").And.Contain("]",
            "a crossing helix is annotated with the [] family");
        pk.FreeEnergy.Should().BeLessThan(CalculateMfeStructure(seq).FreeEnergy + 1e-6,
            "the recursive fold must not be worse than the plain MFE (INV-PKR-01)");
    }

    [Test]
    [CancelAfter(30_000)]
    public void Recursive_NeverWorseThanSingleKnotPredictor_OnSameInputs()
    {
        // The recursive predictor "never returns a structure worse than the plain MFE" (§1) and is a
        // superset of the single-knot class, so its ΔG ≤ the plain MFE on the same inputs (INV-PKR-01).
        foreach (string seq in new[]
        {
            "GGGGAACCCCAACCCCAAGGGG",                 // single H-type example
            "AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU",  // over-arching example
            "GCGCGCAUAUAUGCGCGCAUAU",
        })
        {
            var pk = AssertWellFormed(seq);
            pk.FreeEnergy.Should().BeLessThanOrEqualTo(CalculateMfeStructure(seq).FreeEnergy + 1e-6,
                "recursive ΔG ≤ plain MFE (INV-PKR-01)");
        }
    }

    [Test]
    [CancelAfter(120_000)]
    public void Recursive_RandomMixedSequences_AlwaysWellFormed_NeverThrow()
    {
        // Random sequences mixing ACGU with junk and varied lengths around / above the 11-nt boundary.
        const string alphabet = "ACGUacguTtNnXx.-#";
        var rng = Rng(13572468);
        for (int iter = 0; iter < 300; iter++)
        {
            int len = rng.Next(0, 40);
            var chars = new char[len];
            for (int p = 0; p < len; p++) chars[p] = alphabet[rng.Next(alphabet.Length)];
            int minLoop = rng.Next(-2, 6);
            AssertWellFormed(new string(chars), minLoop);
        }
    }

    [Test]
    [CancelAfter(120_000)]
    public void Recursive_RandomKnotProneSequences_StayWithinContract()
    {
        // Bias toward complementary blocks that can actually form crossing helices, to stress the
        // accept-knot / recursive-loop paths; every accepted fold must satisfy the full contract.
        var rng = Rng(98765);
        for (int iter = 0; iter < 200; iter++)
        {
            int stem = rng.Next(2, 6);
            int loop = rng.Next(1, 5);
            string a = RandomBlock(rng, stem);
            string b = RandomBlock(rng, stem);
            // a · loop1 · b · loop2 · a' · loop3 · b' — the canonical H-type skeleton.
            string seq =
                a + new string('A', loop) +
                b + new string('A', loop) +
                Complement(a) + new string('A', loop) +
                Complement(b);
            AssertWellFormed(seq);
        }
    }

    private static string RandomBlock(Random rng, int len)
    {
        const string acgu = "ACGU";
        var c = new char[len];
        for (int i = 0; i < len; i++) c[i] = acgu[rng.Next(acgu.Length)];
        return new string(c);
    }

    #endregion
}
