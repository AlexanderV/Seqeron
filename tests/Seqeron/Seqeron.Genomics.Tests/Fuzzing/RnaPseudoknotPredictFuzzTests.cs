// RNA-PKPREDICT-001 — RNA pseudoknot STRUCTURE PREDICTION (canonical H-type, pknotsRG class).
// Fuzz tests (strategy BE = Boundary Exploitation), row 236 of docs/checklists/03_FUZZING.md.
// Algorithm doc: docs/algorithms/RnaStructure/Pseudoknot_Prediction.md
// Canonical tests: tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_PredictStructurePseudoknot_Tests.cs
// Evidence: docs/Evidence/RNA-PKPREDICT-001-Evidence.md
// Source: RnaSecondaryStructure.PredictStructurePseudoknot(string, int) — RnaSecondaryStructure.cs.
//         Reeder & Giegerich (2004) BMC Bioinformatics 5:104; Turner 2004 NN; Antczak et al. (2018).

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for RNA-PKPREDICT-001 —
/// <see cref="RnaSecondaryStructure.PredictStructurePseudoknot(string,int)"/>, the thermodynamic
/// predictor of a single canonical H-type pseudoknot (two crossing helices a·a', b·b' with three
/// intervening loops u, v, w) of the pknotsRG canonical simple-recursive class. The two helices are
/// scored with the Turner 2004 nearest-neighbour stacking model, the three loops fold with the
/// pseudoknot-free MFE, and the pknotsRG penalties are added (initiation 9.0; 0.3 per unpaired loop
/// nucleotide). A candidate knot is accepted ONLY if its ΔG is strictly below the plain MFE.
/// Lives in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Documented contract (Pseudoknot_Prediction.md §2.4, §3, §6.1) — what every result MUST satisfy
/// ───────────────────────────────────────────────────────────────────────────
///   • INV-PK-01  Returned ΔG ≤ ΔG of the plain pseudoknot-free MFE for the same sequence/params.
///   • INV-PK-02  HasPseudoknot ⇒ the base-pair set contains a CROSSING pair (∃ (i,j),(k,l) with
///                i &lt; k &lt; j &lt; l). Conversely a single bracket family ⇒ HasPseudoknot=false.
///   • INV-PK-03  Every position is paired AT MOST ONCE; all indices in [0, n); i &lt; j for each pair.
///   • INV-PK-04  No spurious pseudoknot — a knot is accepted only when strictly below the MFE.
///   • Complementarity (§3.3, BuildPairLookup): every reported pair (i,j) is a legal A-U / G-C / G-U
///                wobble pair on the FOLDED sequence (upper-cased, T read as U). No non-complementary
///                pair may ever appear.
///   • Output shape (§3.2): Sequence upper-cased + T→U; DotBracket has length n, only chars
///                {(),[],.} , each family balanced, and brackets correspond exactly to BasePairs;
///                FreeEnergy is a finite double; BasePairs sorted by 5' position.
///   • Edge cases (§6.1): null / empty / &lt; 11 nt → empty pseudoknot-free structure (no pairs,
///                all dots, ΔG = 0, HasPseudoknot=false). Non-ACGU characters simply do not pair.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"
/// ───────────────────────────────────────────────────────────────────────────
/// Feed malformed / boundary sequences and assert the predictor NEVER fails undisciplined: it must
/// not hang or infinite-loop (the H-type enumeration is O(n³)+ — every hang-sensitive test carries
/// [CancelAfter]); must not throw an unhandled runtime exception (IndexOutOfRange / NullReference /
/// overflow / NaN); and must not emit out-of-contract nonsense — a non-complementary base pair, an
/// index out of bounds, a position paired twice, a dot-bracket that disagrees with BasePairs, a
/// HasPseudoknot=true with no crossing pair, or a non-finite ΔG. Every input → EITHER a well-defined
/// theory-correct structure OR (here) the documented empty result; no validation exception is thrown
/// because the contract maps null/empty/too-short to the empty structure, not to an ArgumentException.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy BE = Boundary Exploitation — row 236
/// targets "empty, too-short, all-A (no pairing), non-ACGU"
/// ───────────────────────────────────────────────────────────────────────────
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, -1, MaxInt, empty), row 236.
///   • empty / null (BE)        → empty pseudoknot-free structure, ΔG 0, no crash.
///   • too-short &lt; 11 nt (BE)   → no pseudoknot (returns the plain MFE for that length).
///   • all-A (BE)               → A cannot pair A → NO base pair at all → unpaired structure, ΔG 0,
///                                 HasPseudoknot=false; must NOT crash and must NOT invent a pair.
///   • non-ACGU (BE)            → unknown letters never pair; predictor must not throw and any pair
///                                 it does report must still be a legal A-U/G-C/G-U pair.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaPseudoknotPredictFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static Random Rng(int seed) => new(seed);

    /// <summary>The shortest canonical H-type knot length (2·2 + 2·2 + 3 loops), per §6.1.</summary>
    private const int MinPseudoknotLength = 11;

    /// <summary>
    /// Asserts the FULL documented well-formedness of a predicted structure for the given input
    /// sequence, deriving every invariant independently from Pseudoknot_Prediction.md:
    ///   • never throws (the predictor maps all inputs to a structure, never an exception);
    ///   • Sequence is the upper-cased, T→U spelling of the (non-null) input;
    ///   • FreeEnergy is finite (no NaN / ±∞);
    ///   • every pair (i,j): 0 ≤ i &lt; j &lt; n  AND  CanPair(seq[i],seq[j]) — legal A-U/G-C/G-U (INV-PK-03 + §3.3);
    ///   • no position is paired more than once (INV-PK-03);
    ///   • ΔG ≤ plain-MFE ΔG for the same sequence/minLoopSize (INV-PK-01);
    ///   • HasPseudoknot ⇔ the base-pair set actually contains a crossing pair (INV-PK-02 + INV-PK-04 ⇒
    ///       a [] family appears iff HasPseudoknot);
    ///   • DotBracket has length n, uses only {(,),[,],.}, both bracket families balanced, and the
    ///       paired/unpaired columns agree exactly with BasePairs.
    /// </summary>
    private static PseudoknotStructure AssertWellFormed(string? input, int minLoopSize = 3)
    {
        PseudoknotStructure pk = default;
        Action act = () => pk = PredictStructurePseudoknot(input!, minLoopSize);
        act.Should().NotThrow("the predictor maps every input to a structure, never an exception");

        string seq = pk.Sequence;
        int n = seq.Length;

        // Folded spelling: upper-cased, T→U; length preserved; empty for null.
        string expectedSeq = input is null ? "" : input.ToUpperInvariant().Replace('T', 'U');
        seq.Should().Be(expectedSeq, "the folded sequence is the upper-cased T→U spelling of the input");

        double.IsFinite(pk.FreeEnergy).Should().BeTrue("ΔG must be a finite number (no NaN/∞)");

        // INV-PK-03 + §3.3 complementarity: indices in range, i<j, legal pair, no double-pairing.
        var seen = new HashSet<int>();
        foreach (var (a, b) in pk.BasePairs)
        {
            int i = Math.Min(a, b), j = Math.Max(a, b);
            i.Should().BeGreaterThanOrEqualTo(0, "pair index must be ≥ 0");
            j.Should().BeLessThan(n, "pair index must be < sequence length (in bounds)");
            i.Should().BeLessThan(j, "a base pair spans two DISTINCT positions (i<j)");
            CanPair(seq[i], seq[j]).Should().BeTrue(
                $"reported pair ({i},{j}) = ({seq[i]},{seq[j]}) must be a legal A-U/G-C/G-U pair");
            seen.Add(a).Should().BeTrue($"position {a} must not be paired twice (INV-PK-03)");
            seen.Add(b).Should().BeTrue($"position {b} must not be paired twice (INV-PK-03)");
        }

        // INV-PK-01: never worse than the plain pseudoknot-free MFE.
        double plainMfe = CalculateMfeStructure(input!, minLoopSize).FreeEnergy;
        pk.FreeEnergy.Should().BeLessThanOrEqualTo(plainMfe + 1e-6,
            "the returned ΔG can never exceed the plain-MFE baseline (INV-PK-01)");

        // INV-PK-02: HasPseudoknot ⇔ a genuine crossing pair exists.
        HasCrossingPair(pk.BasePairs).Should().Be(pk.HasPseudoknot,
            "HasPseudoknot must be true exactly when the pair set contains a crossing pair (INV-PK-02/04)");

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
    /// length n; only {(,),[,],.}; both families balanced; and every column is '.' iff unpaired and
    /// a matching open/close bracket iff it is the 5'/3' end of a reported pair (§3.2).
    /// </summary>
    private static void AssertDotBracketConsistent(PseudoknotStructure pk, int n)
    {
        string db = pk.DotBracket;
        db.Length.Should().Be(n, "dot-bracket length must equal the sequence length");
        db.Should().MatchRegex(@"^[\(\)\[\]\.]*$", "dot-bracket uses only ( ) [ ] .");

        // Both bracket families must be balanced.
        Balanced(db, '(', ')').Should().BeTrue("the () family must be balanced");
        Balanced(db, '[', ']').Should().BeTrue("the [] family must be balanced");

        // Column-level agreement with the pair set.
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

    #endregion

    #region RNA-PKPREDICT-001 — empty / null (BE)

    [Test]
    [CancelAfter(10_000)]
    public void Predict_Null_ReturnsEmptyPseudoknotFreeStructure_NoCrash()
    {
        var pk = AssertWellFormed(null);
        AssertEmptyStructure(pk);
        pk.Sequence.Should().BeEmpty("null folds to the empty sequence (§6.1 contract parity with MFE)");
    }

    [Test]
    [CancelAfter(10_000)]
    public void Predict_Empty_ReturnsEmptyPseudoknotFreeStructure_NoCrash()
    {
        var pk = AssertWellFormed("");
        AssertEmptyStructure(pk);
        pk.Sequence.Should().BeEmpty();
    }

    #endregion

    #region RNA-PKPREDICT-001 — too-short (< 11 nt, BE)

    [TestCase("G")]
    [TestCase("GC")]
    [TestCase("GCGC")]
    [TestCase("GGGGAACCCC")]      // 10 nt — one below the 11-nt canonical minimum
    [TestCase("GGGGAACCCCA")]     // exactly 11 nt with a single hairpin: still no crossing helix possible
    public void Predict_TooShortForCanonicalKnot_NeverReportsPseudoknot(string seq)
    {
        var pk = AssertWellFormed(seq);
        // Below MinPseudoknotLength the predictor must return the plain MFE with no knot. At exactly
        // 11 nt a knot is geometrically possible only for very specific sequences; a plain hairpin
        // must never be flagged. AssertWellFormed already ties HasPseudoknot to a real crossing pair.
        if (seq.Length < MinPseudoknotLength)
            pk.HasPseudoknot.Should().BeFalse($"{seq.Length} nt < {MinPseudoknotLength} → no canonical knot (§6.1)");
    }

    [Test]
    [CancelAfter(10_000)]
    public void Predict_AllShortLengths_BelowMin_NeverThrow_NeverKnot()
    {
        // Every length 0..10 of mixed pairable bases: no canonical knot can form, no crash.
        const string alphabet = "ACGU";
        var rng = Rng(20260626);
        for (int len = 0; len < MinPseudoknotLength; len++)
        {
            for (int iter = 0; iter < 50; iter++)
            {
                var chars = new char[len];
                for (int p = 0; p < len; p++) chars[p] = alphabet[rng.Next(alphabet.Length)];
                var pk = AssertWellFormed(new string(chars));
                pk.HasPseudoknot.Should().BeFalse($"length {len} < {MinPseudoknotLength} → no knot");
            }
        }
    }

    #endregion

    #region RNA-PKPREDICT-001 — all-A / no pairing possible (BE)

    [TestCase(11)]
    [TestCase(20)]
    [TestCase(40)]
    [TestCase(64)]
    public void Predict_AllAdenine_NoPairingPossible_YieldsUnpairedStructure(int length)
    {
        // A pairs only with U; an all-A sequence has ZERO legal pairs, so the predictor must return
        // the fully-unpaired structure with ΔG 0 and no pseudoknot — NOT crash, NOT invent a pair.
        var pk = AssertWellFormed(new string('A', length));

        pk.BasePairs.Should().BeEmpty("all-A has no complementary positions → no base pair at all");
        pk.HasPseudoknot.Should().BeFalse("no pairs ⇒ no crossing helix");
        pk.FreeEnergy.Should().Be(0.0, "no stabilizing pair ⇒ ΔG = 0");
        pk.DotBracket.Should().Be(new string('.', length), "every position is unpaired");
    }

    [TestCase("CCCCCCCCCCCCCCCC")] // all-C — C pairs only with G
    [TestCase("UUUUUUUUUUUUUUUU")] // all-U — U pairs A and G(wobble) but a homopolymer has neither
    public void Predict_Homopolymer_NoSelfPairing_YieldsUnpairedStructure(string seq)
    {
        var pk = AssertWellFormed(seq);
        pk.BasePairs.Should().BeEmpty("a single-letter homopolymer cannot self-pair");
        pk.HasPseudoknot.Should().BeFalse();
        pk.FreeEnergy.Should().Be(0.0);
    }

    #endregion

    #region RNA-PKPREDICT-001 — non-ACGU characters (BE)

    [TestCase("NNNNNNNNNNNNNNNN")]            // IUPAC ambiguity 'N' — not in the pairing table
    [TestCase("XYZXYZXYZXYZXYZX")]            // arbitrary non-nucleotide letters
    [TestCase("....-----....----")]           // gap / alignment characters
    [TestCase("GGGG##CCCC##GGGG")]            // ACGU islands separated by junk
    [TestCase("123456789012345678")]          // digits
    [TestCase(" gggg  cccc  gggg ")]          // lower-case + spaces
    public void Predict_NonAcguCharacters_NeverThrow_OnlyLegalPairsReported(string seq)
    {
        // Unknown characters must simply not pair; the predictor must not throw and any pair it does
        // report must STILL be a legal A-U/G-C/G-U pair on the folded sequence (AssertWellFormed
        // checks complementarity, in-bounds indices, no double-pairing, dot-bracket agreement).
        AssertWellFormed(seq);
    }

    [Test]
    [CancelAfter(15_000)]
    public void Predict_HighAndControlCodepoints_NeverThrow()
    {
        // Codepoints ≥ 128 and control chars exercise the (b1|b2) ≥ 128 guard in GetBasePairType.
        AssertWellFormed(" GGGGÿCCCC€GGGG\t\n");
        AssertWellFormed(new string('￿', 16));
    }

    #endregion

    #region RNA-PKPREDICT-001 — positive sanity (the doc worked example) + minLoopSize fuzz

    [Test]
    [CancelAfter(15_000)]
    public void Predict_DocumentedHTypeExample_IsAGenuineCrossingKnot()
    {
        // §7.1 worked example: this sequence DOES fold into a canonical H-type knot.
        var pk = AssertWellFormed("GGGGAACCCCAACCCCAAGGGG");
        pk.HasPseudoknot.Should().BeTrue("the §7.1 example folds into an H-type pseudoknot");
        pk.DotBracket.Should().Contain("[").And.Contain("]",
            "a crossing helix is annotated with the [] family");
        // The knot must strictly improve on the plain MFE (INV-PK-04).
        pk.FreeEnergy.Should().BeLessThan(CalculateMfeStructure("GGGGAACCCCAACCCCAAGGGG").FreeEnergy + 1e-6);
    }

    [Test]
    [CancelAfter(20_000)]
    public void Predict_NegativeMinLoopSize_IsClampedAndNeverThrows()
    {
        // §3.1: minLoopSize < 3 is clamped to 3. Negative / zero / int.MinValue must not crash.
        foreach (int m in new[] { int.MinValue, -100, -1, 0, 1, 2 })
            AssertWellFormed("GGGGAACCCCAACCCCAAGGGG", m);
    }

    [Test]
    [CancelAfter(60_000)]
    public void Predict_RandomMixedSequences_AlwaysWellFormed_NeverThrow()
    {
        // Random sequences mixing ACGU with junk and varied lengths around the 11-nt boundary.
        const string alphabet = "ACGUacguTtNnXx.-#";
        var rng = Rng(13572468);
        for (int iter = 0; iter < 600; iter++)
        {
            int len = rng.Next(0, 36);
            var chars = new char[len];
            for (int p = 0; p < len; p++) chars[p] = alphabet[rng.Next(alphabet.Length)];
            int minLoop = rng.Next(-2, 6);
            AssertWellFormed(new string(chars), minLoop);
        }
    }

    [Test]
    [CancelAfter(60_000)]
    public void Predict_RandomKnotProneSequences_StayWithinContract()
    {
        // Bias toward GC/AU-rich palindromic stretches that can actually form crossing helices, to
        // stress the accept-knot path; every accepted knot must still satisfy the full contract.
        var rng = Rng(98765);
        for (int iter = 0; iter < 300; iter++)
        {
            int stem = rng.Next(2, 6);
            int loop = rng.Next(1, 5);
            // a · loop1 · b · loop2 · a' · loop3 · b'  built from complementary blocks.
            string a = RandomBlock(rng, stem);
            string b = RandomBlock(rng, stem);
            string seq =
                a +
                new string('A', loop) +
                b +
                new string('A', loop) +
                Complement(a) +
                new string('A', loop) +
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

    // Antiparallel Watson–Crick complement (A↔U, G↔C), reversed — to seed pairable helices.
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
}
