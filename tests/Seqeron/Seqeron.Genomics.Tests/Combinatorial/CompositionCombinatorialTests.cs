namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Composition area.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What combinatorial testing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Combinatorial testing samples the Cartesian product of an algorithm's
/// independent configuration parameters so that every t-tuple of parameter
/// values (here every PAIR — and, where the grid is small, every full
/// combination) is exercised at least once. The motivation is interaction
/// faults: a rule that holds for each parameter in isolation can still break
/// for a specific COMBINATION of values. The grid is built with NUnit's
/// <c>[Combinatorial]</c> generator (the checklist tool of record); a cell is a
/// concrete (parameter₁, …, parameterₙ) assignment and each cell carries a real
/// pass/fail business assertion, not a smoke check.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description;
///   docs/ADVANCED_TESTING_CHECKLIST.md §10 "Combinatorial / Pairwise Testing".
///
/// The doc's pairwise threshold ("apply when full enumeration > 100 combos") is
/// not reached by these small grids, so the exhaustive <c>[Combinatorial]</c>
/// product is used: it is a strict superset of any pairwise sample and therefore
/// the stronger guarantee at no practical cost.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Composition")]
public class CompositionCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-VALID-001 — Sequence validation (Composition)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 4.
    // Dimensions: alphabet(4: DNA/RNA/Protein/Ambiguous) × strict(2) × seqLen(3).
    //             Full grid 4×2×3 = 24 cells (< 100 → exhaustive, not sampled).
    //
    // Model (Sequence_Validation.md §2.2, §3.3, §5.2): validation is a per-symbol
    // set-membership scan over the alphabet selected by (alphabet, strict).
    //   • strict mode is the documented production surface —
    //       DNA   accepts exactly {A,C,G,T}          (IsValidDna, INV-01),
    //       RNA   accepts exactly {A,C,G,U}          (IsValidRna, INV-02),
    //       Protein accepts the 20 standard residues (ProteinSequence),
    //       IUPAC ambiguity codes are REJECTED        (§5.2, §5.3 deviation).
    //   • non-strict mode is the full IUPAC alphabet (NC-IUB 1984), whose code
    //     set is the production primitive IupacHelper.MatchesIupac —
    //       DNA/Ambiguous additionally accept {N,R,Y,S,W,K,M,B,D,H,V},
    //       RNA likewise but U-based (no T),
    //       Protein additionally accepts the X (unknown) placeholder.
    //
    // The combinatorial point: `strict` and `alphabet` INTERACT. The same symbol
    // 'N' is valid for (Ambiguous, non-strict) yet invalid for (Ambiguous,
    // strict); 'T' is valid DNA but invalid RNA; 'X' is valid non-strict protein
    // but invalid strict protein. Each interaction lives in a distinct grid cell.
    // ═══════════════════════════════════════════════════════════════════════

    public enum Alphabet { Dna, Rna, Protein, Ambiguous }

    /// <summary>The 20 standard amino acids — production set ProteinSequence.StandardAminoAcids.</summary>
    private const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>
    /// Characters that VALIDATE for a given (alphabet, strict) cell. Drawn from
    /// the algorithm's alphabet definition (Sequence_Validation.md §2.1–2.2).
    /// </summary>
    private static string AcceptedAlphabet(Alphabet alphabet, bool strict) => (alphabet, strict) switch
    {
        (Alphabet.Dna, true) => "ACGT",
        (Alphabet.Dna, false) => "ACGTNRYSWKMBDHV",         // DNA + IUPAC ambiguity
        (Alphabet.Rna, true) => "ACGU",
        (Alphabet.Rna, false) => "ACGUNRYSWKMBDHV",         // RNA (U-based) + IUPAC ambiguity
        (Alphabet.Protein, true) => StandardAminoAcids,      // 20 residues, no X
        (Alphabet.Protein, false) => StandardAminoAcids + "X", // + unknown placeholder
        (Alphabet.Ambiguous, true) => "ACGT",                // strict collapses ambiguity → canonical bases
        (Alphabet.Ambiguous, false) => "ACGTNRYSWKMBDHV",    // full IUPAC ambiguity DNA
        _ => throw new ArgumentOutOfRangeException(nameof(alphabet)),
    };

    /// <summary>
    /// A symbol that must be REJECTED in a given (alphabet, strict) cell. Each is
    /// chosen to expose the interaction under test: an ambiguity code where strict
    /// mode forbids it, a cross-alphabet base, or an out-of-IUPAC symbol.
    /// </summary>
    private static char ForeignSymbol(Alphabet alphabet, bool strict) => (alphabet, strict) switch
    {
        (Alphabet.Dna, true) => 'N',           // ambiguity code — strict DNA rejects (§5.2)
        (Alphabet.Dna, false) => 'Z',          // not an IUPAC code at all
        (Alphabet.Rna, true) => 'T',           // DNA base — invalid for RNA (INV-02)
        (Alphabet.Rna, false) => 'Z',
        (Alphabet.Protein, true) => 'X',       // unknown placeholder — strict rejects, non-strict accepts
        (Alphabet.Protein, false) => 'J',      // not a defined amino-acid letter
        (Alphabet.Ambiguous, true) => 'N',     // the defining interaction: ambiguity forbidden when strict
        (Alphabet.Ambiguous, false) => 'Z',
        _ => throw new ArgumentOutOfRangeException(nameof(alphabet)),
    };

    /// <summary>
    /// Production validation router. Strict nucleotide/protein paths are the
    /// shipped validators (IsValidDna / IsValidRna / ProteinSequence membership);
    /// the non-strict paths apply the production IUPAC code primitive
    /// (IupacHelper) that defines the NC-IUB 1984 ambiguity alphabet.
    /// </summary>
    private static bool ValidateUnderModel(string sequence, Alphabet alphabet, bool strict) => alphabet switch
    {
        Alphabet.Dna or Alphabet.Ambiguous when strict => sequence.AsSpan().IsValidDna(),
        Alphabet.Dna or Alphabet.Ambiguous => AllSymbols(sequence, IsValidDnaIupacCode),
        Alphabet.Rna when strict => sequence.AsSpan().IsValidRna(),
        Alphabet.Rna => AllSymbols(sequence, IsValidRnaIupacCode),
        Alphabet.Protein when strict => AllSymbols(sequence, c => ProteinSequence.StandardAminoAcids.Contains(char.ToUpperInvariant(c))),
        Alphabet.Protein => AllSymbols(sequence, c => ProteinSequence.ValidCharacters.Contains(char.ToUpperInvariant(c))),
        _ => throw new ArgumentOutOfRangeException(nameof(alphabet)),
    };

    private static bool AllSymbols(string sequence, Func<char, bool> predicate)
    {
        foreach (char c in sequence)
            if (!predicate(c))
                return false;
        return true;
    }

    /// <summary>True iff <paramref name="c"/> is a valid IUPAC DNA code (A,C,G,T + ambiguity).</summary>
    private static bool IsValidDnaIupacCode(char c)
    {
        try
        {
            // IupacHelper.MatchesIupac throws ArgumentOutOfRangeException for non-codes;
            // probing with a fixed nucleotide turns "is a valid code" into a boolean.
            IupacHelper.MatchesIupac('A', char.ToUpperInvariant(c));
            return true;
        }
        catch (ArgumentOutOfRangeException)
        {
            return false;
        }
    }

    /// <summary>True iff <paramref name="c"/> is a valid IUPAC RNA code (U-based: U replaces T).</summary>
    private static bool IsValidRnaIupacCode(char c)
    {
        char u = char.ToUpperInvariant(c);
        if (u == 'U') return true;          // RNA base
        if (u == 'T') return false;         // T is DNA-only, not part of the RNA alphabet
        return IsValidDnaIupacCode(u);      // shared ambiguity codes {N,R,Y,S,W,K,M,B,D,H,V} and A,C,G
    }

    private static string BuildSequence(string pool, int length)
    {
        if (length == 0) return string.Empty;
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = pool[i % pool.Length];
        return new string(chars);
    }

    /// <summary>
    /// Pairwise grid: every (alphabet × strict × seqLen) cell. A sequence built
    /// entirely from the cell's accepted alphabet MUST validate; the same
    /// sequence with one foreign symbol injected MUST fail. This pins the
    /// alphabet × strict interaction across all three length regimes.
    /// </summary>
    [Test, Combinatorial]
    public void SeqValid_AcceptsOwnAlphabet_RejectsForeignSymbol(
        [Values(Alphabet.Dna, Alphabet.Rna, Alphabet.Protein, Alphabet.Ambiguous)] Alphabet alphabet,
        [Values(true, false)] bool strict,
        [Values(0, 8, 120)] int seqLen)
    {
        string pool = AcceptedAlphabet(alphabet, strict);
        string positive = BuildSequence(pool, seqLen);

        ValidateUnderModel(positive, alphabet, strict)
            .Should().BeTrue($"a length-{seqLen} sequence of only {alphabet}/strict={strict} symbols \"{positive}\" is valid");

        if (seqLen == 0)
            return; // empty sequence is vacuously valid (§3.3); no symbol to corrupt

        char foreign = ForeignSymbol(alphabet, strict);
        char[] corrupted = positive.ToCharArray();
        corrupted[seqLen / 2] = foreign;
        string negative = new string(corrupted);

        ValidateUnderModel(negative, alphabet, strict)
            .Should().BeFalse($"foreign symbol '{foreign}' at the middle of \"{negative}\" must fail {alphabet}/strict={strict} validation");
    }

    /// <summary>
    /// Interaction witness: the SAME ambiguity symbol 'N' flips validity purely on
    /// the `strict` axis for the Ambiguous alphabet — the canonical reason this
    /// unit needs a combinatorial (not one-parameter-at-a-time) suite. Asserted
    /// across all three length regimes to rule out a length-coupled shortcut.
    /// </summary>
    [Test, Combinatorial]
    public void SeqValid_AmbiguityCode_ValidIffNonStrict(
        [Values(4, 16, 64)] int seqLen)
    {
        char[] chars = BuildSequence("ACGT", seqLen).ToCharArray();
        chars[seqLen / 2] = 'N';                  // guarantee the ambiguity code is present
        string ambiguous = new string(chars);

        ValidateUnderModel(ambiguous, Alphabet.Ambiguous, strict: false)
            .Should().BeTrue("the full IUPAC alphabet (NC-IUB 1984) accepts 'N'");
        ValidateUnderModel(ambiguous, Alphabet.Ambiguous, strict: true)
            .Should().BeFalse("strict mode rejects ambiguity codes (Sequence_Validation.md §5.2)");
    }

    /// <summary>
    /// Worked truth-table witnesses from Sequence_Validation.md §5.2 so the grid
    /// cannot pass vacuously: production IsValidDna/IsValidRna pinned on the exact
    /// documented examples, including case-insensitivity and N rejection.
    /// </summary>
    [Test]
    public void SeqValid_DocumentedTruthTable()
    {
        "".AsSpan().IsValidDna().Should().BeTrue();
        "".AsSpan().IsValidRna().Should().BeTrue();
        "ACGT".AsSpan().IsValidDna().Should().BeTrue();
        "ACGT".AsSpan().IsValidRna().Should().BeFalse();
        "ACGU".AsSpan().IsValidDna().Should().BeFalse();
        "ACGU".AsSpan().IsValidRna().Should().BeTrue();
        "acgt".AsSpan().IsValidDna().Should().BeTrue();   // case-insensitive (§3.3)
        "ACGN".AsSpan().IsValidDna().Should().BeFalse();  // ambiguity code rejected
        "ACGN".AsSpan().IsValidRna().Should().BeFalse();
        "AC GT".AsSpan().IsValidDna().Should().BeFalse(); // whitespace rejected
    }
}
