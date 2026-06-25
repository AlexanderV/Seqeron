// RNA-PAIR-001 — RNA base-pairing rule / can-pair predicate.
// Fuzz tests (strategy BE = Boundary Exploitation, MC = Malformed Content).
// Algorithm doc: docs/algorithms/RnaStructure/RNA_Base_Pairing.md
// Canonical tests: tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_CanPair_Tests.cs (RNA-PAIR-001)
// Evidence: docs/Evidence/RNA-PAIR-001-Evidence.md
// Source: RnaSecondaryStructure.CanPair / GetBasePairType — RnaSecondaryStructure.cs.
//         Crick (1966) JMB 19(2):548 (wobble); IUPAC-IUB (1970); Biopython complement_rna.

using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for RNA-PAIR-001 — <see cref="RnaSecondaryStructure.CanPair(char,char)"/> and
/// <see cref="RnaSecondaryStructure.GetBasePairType(char,char)"/>, the O(1) RNA base-pairing
/// predicate: whether two ribonucleotides form a hydrogen-bonded pair, and whether that pair is
/// Watson-Crick or wobble.
/// Lives in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Documented contract (RNA_Base_Pairing.md §2.2, §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   • CanPair(b1,b2) == true IFF the unordered pair is one of {A,U}, {G,C}, {G,U} — i.e. the six
///       ordered pairs A·U, U·A, G·C, C·G (Watson-Crick) and G·U, U·G (wobble). Nothing else over
///       the RNA alphabet pairs (§2.2): A pairs only with U, C pairs only with G.
///   • GetBasePairType returns WatsonCrick for {A,U}/{G,C}, Wobble for {G,U}, null otherwise.
///   • INV-01 — symmetry: CanPair(x,y) == CanPair(y,x) (the lookup table is seeded symmetrically).
///   • INV-02 — symmetry: GetBasePairType(x,y) == GetBasePairType(y,x).
///   • INV-03 — CanPair(x,y) == (GetBasePairType(x,y) != null).
///   • INV-04 — G·U and U·G are Wobble, NEVER WatsonCrick (§2.4).
///   • Case-insensitive: inputs are upper-cased before lookup, so lowercase a/c/g/u behave
///       identically to uppercase (§3.3, §6.1).
///   • DNA T is NOT an RNA base for pairing: it does NOT pair (CanPair false / type null) — even
///       T·A, despite GetComplement treating T as U (§3.3, §6.1 "DNA T in CanPair").
///   • Out-of-domain chars (non-RNA letters N/X, digits, punctuation, gap '-'/'.'/' ', control,
///       Unicode, surrogate, char 127, the full 0–65535 range) return false/null with NO exception:
///       the lookup is bounds-checked ((b1|b2) < 128) so no KeyNotFound / IndexOutOfRange (§3.3,
///       §6.1 "Out-of-ASCII char in CanPair").
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"
/// ───────────────────────────────────────────────────────────────────────────
/// Feed degenerate / out-of-domain chars and assert the predicate NEVER fails undisciplined:
/// no exception (no KeyNotFound, no IndexOutOfRange on the 128×128 table), no asymmetry, no
/// case-sensitivity bug, no wobble misclassification. Every input resolves to a documented,
/// deterministic bool / BasePairType?.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy BE = Boundary Exploitation / MC = Malformed Content — row 153
/// targets "non-RNA base, lowercase, gap char"
/// ───────────────────────────────────────────────────────────────────────────
/// — docs/checklists/03_FUZZING.md §Description (BE, MC), row 153.
///
///   • non-RNA base (MC) — a char not in {A,C,G,U}: DNA T, ambiguity N/X, digits, junk letters →
///       does NOT pair with anything, returns false / null, no crash.
///   • lowercase (MC) — a/c/g/u → case-insensitive, pairs identically to uppercase.
///   • gap char (BE) — alignment gap/spacer '-', '.', ' ' → does NOT pair, false / null, no crash.
///   • boundary chars (BE) — '\0', char 127 (ASCII boundary), char 128 (just past the table),
///       char.MaxValue, surrogate halves → bounds-checked false/null, no IndexOutOfRange.
///
/// Watched failure modes: KeyNotFound / IndexOutOfRange on a non-RNA or gap char; case-sensitivity
/// bug (lowercase not pairing); asymmetry; G·U returned as WatsonCrick rather than Wobble; T
/// wrongly accepted as an RNA base in CanPair.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaBasePairFuzzTests
{
    /// <summary>The standard RNA alphabet (uppercase).</summary>
    private const string RnaBases = "ACGU";

    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static Random Rng(int seed) => new(seed);

    /// <summary>
    /// Asserts the full documented contract for one (b1,b2) probe: the two query forms agree
    /// (INV-03), the result is symmetric (INV-01/INV-02), and — when a pair forms — the type is the
    /// expected Watson-Crick / wobble class. Never throws regardless of the input chars.
    /// </summary>
    private static void AssertWellFormed(char b1, char b2)
    {
        bool canPair = CanPair(b1, b2);
        BasePairType? type = GetBasePairType(b1, b2);

        // INV-03: CanPair iff a type is defined.
        canPair.Should().Be(type is not null,
            $"CanPair and GetBasePairType must agree for ('{b1}','{b2}')");

        // INV-01 / INV-02: symmetry.
        CanPair(b2, b1).Should().Be(canPair,
            $"CanPair must be symmetric for ('{b1}','{b2}')");
        GetBasePairType(b2, b1).Should().Be(type,
            $"GetBasePairType must be symmetric for ('{b1}','{b2}')");
    }

    /// <summary>
    /// Independent oracle for the documented pairing rule over the RNA alphabet (T NOT treated as
    /// U for pairing). Returns the expected BasePairType? for two chars after case-folding.
    /// </summary>
    private static BasePairType? Expected(char b1, char b2)
    {
        char u1 = char.ToUpperInvariant(b1);
        char u2 = char.ToUpperInvariant(b2);
        bool Pair(char a, char b) => (u1 == a && u2 == b) || (u1 == b && u2 == a);
        if (Pair('A', 'U') || Pair('G', 'C')) return BasePairType.WatsonCrick;
        if (Pair('G', 'U')) return BasePairType.Wobble;
        return null;
    }

    #endregion

    #region RNA-PAIR-001 — CanPair / GetBasePairType

    // ───────────────────────────────────────────────────────────────────────
    // POSITIVE sanity — the documented valid pairs and the full 16-combination
    // truth table over {A,C,G,U} (RNA_Base_Pairing.md §2.2, §7.1).
    // ───────────────────────────────────────────────────────────────────────
    #region Positive sanity & truth table

    [TestCase('A', 'U')]
    [TestCase('U', 'A')]
    [TestCase('G', 'C')]
    [TestCase('C', 'G')]
    [TestCase('G', 'U')]
    [TestCase('U', 'G')]
    public void CanPair_CanonicalValidPairs_ReturnTrue(char b1, char b2)
    {
        CanPair(b1, b2).Should().BeTrue($"('{b1}','{b2}') is a documented RNA pair");
    }

    [TestCase('A', 'U', BasePairType.WatsonCrick)]
    [TestCase('U', 'A', BasePairType.WatsonCrick)]
    [TestCase('G', 'C', BasePairType.WatsonCrick)]
    [TestCase('C', 'G', BasePairType.WatsonCrick)]
    [TestCase('G', 'U', BasePairType.Wobble)]
    [TestCase('U', 'G', BasePairType.Wobble)]
    public void GetBasePairType_CanonicalValidPairs_ReturnDocumentedType(char b1, char b2, BasePairType expected)
    {
        // INV-04 in particular: G·U / U·G are Wobble, never WatsonCrick.
        GetBasePairType(b1, b2).Should().Be(expected);
    }

    [TestCase('A', 'A')]
    [TestCase('A', 'G')]
    [TestCase('A', 'C')]
    [TestCase('C', 'U')]
    [TestCase('C', 'C')]
    [TestCase('C', 'A')]
    [TestCase('G', 'G')]
    [TestCase('G', 'A')]
    [TestCase('U', 'U')]
    [TestCase('U', 'C')]
    public void CanPair_InvalidCombinationsOverAlphabet_ReturnFalse(char b1, char b2)
    {
        CanPair(b1, b2).Should().BeFalse($"('{b1}','{b2}') is not an RNA pair");
        GetBasePairType(b1, b2).Should().BeNull();
    }

    [Test]
    public void CanPair_FullAlphabetTruthTable_MatchesOracleAndIsSymmetric()
    {
        // All 16 ordered combinations over {A,C,G,U}: exactly 6 must pair, all symmetric.
        int pairing = 0;
        foreach (char b1 in RnaBases)
        foreach (char b2 in RnaBases)
        {
            BasePairType? expected = Expected(b1, b2);
            GetBasePairType(b1, b2).Should().Be(expected, $"truth table for ('{b1}','{b2}')");
            CanPair(b1, b2).Should().Be(expected is not null);
            AssertWellFormed(b1, b2);
            if (expected is not null) pairing++;
        }
        pairing.Should().Be(6, "exactly A·U, U·A, G·C, C·G, G·U, U·G pair over the RNA alphabet");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────
    // MC = Malformed Content — non-RNA bases.
    // A char not in {A,C,G,U} never pairs and never throws (RNA_Base_Pairing.md
    // §3.3, §6.1). T in particular is a DNA base and does NOT pair, even with A.
    // ───────────────────────────────────────────────────────────────────────
    #region MC — non-RNA base

    [TestCase('T')]   // DNA thymine — NOT an RNA base for pairing
    [TestCase('N')]   // IUPAC any
    [TestCase('X')]   // unknown
    [TestCase('R')]   // IUPAC degenerate (purine)
    [TestCase('Y')]   // IUPAC degenerate (pyrimidine)
    [TestCase('I')]   // inosine (out of scope per §5.3 / §6.2)
    [TestCase('Z')]
    [TestCase('B')]
    [TestCase('5')]   // digit
    [TestCase('@')]   // punctuation / junk
    public void CanPair_NonRnaBaseAgainstEveryRnaBase_NeverPairs(char alien)
    {
        foreach (char rna in RnaBases)
        {
            CanPair(alien, rna).Should().BeFalse($"'{alien}' is not an RNA base");
            CanPair(rna, alien).Should().BeFalse($"'{alien}' is not an RNA base (symmetric)");
            GetBasePairType(alien, rna).Should().BeNull();
            GetBasePairType(rna, alien).Should().BeNull();
            AssertWellFormed(alien, rna);
        }
        // Also against itself and another alien.
        CanPair(alien, alien).Should().BeFalse();
        AssertWellFormed(alien, alien);
    }

    [Test]
    public void CanPair_DnaT_DoesNotPairEvenWithA_DespiteComplementTreatingTasU()
    {
        // §6.1: T is NOT an RNA base in CanPair — distinct from GetComplement, where T → A.
        CanPair('T', 'A').Should().BeFalse();
        CanPair('A', 'T').Should().BeFalse();
        CanPair('t', 'a').Should().BeFalse("case-folding must not resurrect T as a pairing base");
        GetBasePairType('T', 'A').Should().BeNull();
        // Sanity: GetComplement DOES treat T as U (separate contract).
        GetComplement('T').Should().Be('A');
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────
    // MC = Malformed Content — lowercase.
    // Inputs are upper-cased before lookup: lowercase behaves identically to
    // uppercase (RNA_Base_Pairing.md §3.3, §6.1 "Lowercase input").
    // ───────────────────────────────────────────────────────────────────────
    #region MC — lowercase / mixed case

    [TestCase('a', 'u', BasePairType.WatsonCrick)]
    [TestCase('g', 'c', BasePairType.WatsonCrick)]
    [TestCase('g', 'u', BasePairType.Wobble)]
    [TestCase('A', 'u', BasePairType.WatsonCrick)]  // mixed case
    [TestCase('g', 'C', BasePairType.WatsonCrick)]
    [TestCase('u', 'G', BasePairType.Wobble)]
    public void CanPair_LowercaseAndMixedCase_PairLikeUppercase(char b1, char b2, BasePairType expected)
    {
        CanPair(b1, b2).Should().BeTrue();
        GetBasePairType(b1, b2).Should().Be(expected);
    }

    [Test]
    public void CanPair_CaseInsensitive_AgreesWithUppercaseForEntireAlphabet()
    {
        foreach (char b1 in RnaBases)
        foreach (char b2 in RnaBases)
        {
            char l1 = char.ToLowerInvariant(b1);
            char l2 = char.ToLowerInvariant(b2);
            CanPair(l1, l2).Should().Be(CanPair(b1, b2), $"lowercase ('{l1}','{l2}') must match upper");
            GetBasePairType(l1, l2).Should().Be(GetBasePairType(b1, b2));
            AssertWellFormed(l1, l2);
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────
    // BE = Boundary Exploitation — gap chars and ASCII/Unicode boundaries.
    // Gap/spacer chars never pair; the bounds-checked lookup ((b1|b2) < 128)
    // returns false/null for any char incl. > 127, with NO exception
    // (RNA_Base_Pairing.md §3.3, §6.1 "Out-of-ASCII char in CanPair").
    // ───────────────────────────────────────────────────────────────────────
    #region BE — gap char & boundary chars

    [TestCase('-')]   // alignment gap
    [TestCase('.')]   // dot-bracket unpaired / gap
    [TestCase(' ')]   // space
    [TestCase('*')]
    [TestCase('~')]
    public void CanPair_GapChars_NeverPair(char gap)
    {
        foreach (char rna in RnaBases)
        {
            CanPair(gap, rna).Should().BeFalse($"'{gap}' is a gap char, not a base");
            CanPair(rna, gap).Should().BeFalse();
            GetBasePairType(gap, rna).Should().BeNull();
            AssertWellFormed(gap, rna);
        }
        CanPair(gap, gap).Should().BeFalse();
    }

    [TestCase('\0')]                 // null char — lower boundary
    [TestCase((char)127)]           // DEL — top of ASCII, in-table but unseeded
    [TestCase((char)128)]           // just past the 128×128 table — bounds guard
    [TestCase((char)255)]           // Latin-1
    [TestCase('￿')]            // char.MaxValue
    [TestCase('\uD800')]            // high surrogate half
    [TestCase('\uDC00')]            // low surrogate half
    [TestCase('λ')]                 // Greek
    [TestCase('日')]                // CJK
    public void CanPair_BoundaryAndOutOfAsciiChars_ReturnFalseWithNoException(char c)
    {
        foreach (char rna in RnaBases)
        {
            Action probe = () => { CanPair(c, rna); GetBasePairType(c, rna); };
            probe.Should().NotThrow($"bounds-checked lookup must not throw for char U+{(int)c:X4}");
            CanPair(c, rna).Should().BeFalse();
            CanPair(rna, c).Should().BeFalse();
            GetBasePairType(c, rna).Should().BeNull();
            AssertWellFormed(c, rna);
        }
        CanPair(c, c).Should().BeFalse();
    }

    [Test]
    public void CanPair_EntireCharRange_NeverThrowsAndOnlyDocumentedPairsHold()
    {
        // Exhaustive boundary sweep of all 65536 chars against each RNA base: the ONLY chars that
        // pair are the documented A/C/G/U (upper and lower), and nothing ever throws.
        foreach (char rna in RnaBases)
        {
            for (int code = 0; code <= char.MaxValue; code++)
            {
                char c = (char)code;
                bool can = CanPair(c, rna);
                can.Should().Be(Expected(c, rna) is not null,
                    $"char U+{code:X4} vs '{rna}' must follow the documented rule");
                CanPair(rna, c).Should().Be(can, "symmetry across the full char range");
            }
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────
    // Randomized fuzz — random char pairs (incl. junk, gaps, Unicode) must
    // always satisfy the well-formedness + oracle contract, never throwing.
    // ───────────────────────────────────────────────────────────────────────
    #region Randomized fuzz

    [Test]
    public void CanPair_RandomCharPairs_AlwaysMatchOracleAndNeverThrow()
    {
        var rng = Rng(20260620);
        for (int i = 0; i < 20_000; i++)
        {
            char b1 = (char)rng.Next(0, char.MaxValue + 1);
            char b2 = (char)rng.Next(0, char.MaxValue + 1);

            Action probe = () => { CanPair(b1, b2); GetBasePairType(b1, b2); };
            probe.Should().NotThrow($"no exception for (U+{(int)b1:X4}, U+{(int)b2:X4})");

            BasePairType? expected = Expected(b1, b2);
            GetBasePairType(b1, b2).Should().Be(expected);
            CanPair(b1, b2).Should().Be(expected is not null);
            AssertWellFormed(b1, b2);
        }
    }

    [Test]
    public void CanPair_RandomDegenerateAlphabet_NeverPairsOutsideDocumentedSet()
    {
        // Bias the alphabet toward the malformed/boundary targets so most probes are junk:
        // lowercase, T, N, gaps, digits — none of which (except a/c/g/u) may pair.
        const string alphabet = "acguACGU TtNn-._@5RYxX";
        var rng = Rng(777);
        for (int i = 0; i < 10_000; i++)
        {
            char b1 = alphabet[rng.Next(alphabet.Length)];
            char b2 = alphabet[rng.Next(alphabet.Length)];
            BasePairType? expected = Expected(b1, b2);
            GetBasePairType(b1, b2).Should().Be(expected, $"('{b1}','{b2}')");
            CanPair(b1, b2).Should().Be(expected is not null);
            AssertWellFormed(b1, b2);
        }
    }

    #endregion

    #endregion
}
