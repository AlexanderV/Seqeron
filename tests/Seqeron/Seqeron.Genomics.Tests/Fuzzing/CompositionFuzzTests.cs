using System;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Composition area — GC content.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, OverflowException, …). Every input must result in
/// EITHER a well-defined, theory-correct value, OR a *documented, intentional*
/// validation exception (ArgumentException / ArgumentNullException). A raw
/// runtime exception or a hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-GC-001 — GC content (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 1.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single char, extremely long.
///   • INJ = Injection — non-ACGT characters, null, unicode (combining marks,
///           astral/surrogate-pair code points, null byte '\0').
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The GC% contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// GC% = (G + C) / (A + T + G + C) × 100, case-insensitive.
///   — docs/algorithms/Statistics/GC_Content_Profile.md §2.2 [Wikipedia GC-content];
///     docs/algorithms/Sequence_Composition/Sequence_Composition.md §2.2.
///
/// API entry: DnaSequence.GcContent()
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs), backed by
///   SequenceExtensions.CalculateGcContent
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs).
///
/// Documented input handling along the public DnaSequence path
/// (DnaSequence.cs lines 22–33, 112–124):
///   • null or empty string  → an empty sequence; GcContent() == 0
///     (string.IsNullOrEmpty short-circuit). This is a *defined result*, NOT an
///     exception: the public surface treats null as "no sequence".
///   • input is case-folded with ToUpperInvariant, then validated; so lowercase
///     a/c/g/t round-trips to the same GC% as uppercase.
///   • ANY character that is not A/C/G/T after upper-casing (digits, whitespace,
///     N/IUPAC ambiguity codes, U, '\0', unicode letters, combining marks,
///     astral code points) → a *documented, intentional* ArgumentException from
///     ValidateSequence. This is the contract's validation gate, not a crash.
///   • a valid, extremely long sequence → computed without overflow/hang;
///     result stays in [0, 100].
///
/// The backing span primitive SequenceExtensions.CalculateGcContent is *lenient*
/// by separate documented design: it excludes non-A/T/G/C/U symbols from BOTH
/// numerator and denominator and returns 0 when no valid base is present
/// (SequenceExtensions.cs lines 17–58, matching Biopython gc_fraction "remove").
/// We pin that contract too, so the boundary between the strict public API and
/// the lenient primitive is explicit and cannot silently drift.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-COMP-001 — DNA complement (Composition)
/// Checklist: docs/checklists/03_FUZZING.md, row 2.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — empty string, single char, null.
///   • INJ = Injection — non-DNA characters, mixed case, unicode (accented
///           Latin, Greek, combining marks, full-width look-alikes, and
///           astral/surrogate-pair code points).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The DNA-complement contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// DNA complement maps A↔T and C↔G (Watson–Crick base pairing).
///   — docs/algorithms/Sequence_Composition/RNA_Complement.md §2.1–2.2 (the DNA
///     sibling GetComplementBase, explicitly cited there as SEQ-COMP-001);
///     docs/mcp/tools/sequence/complement_base.md.
///
/// SEQ-COMP-001 has TWO documented surfaces with DIFFERENT, intentional contracts.
/// Fuzzing must pin both, and pin the boundary between them so neither can drift:
///
/// (1) The STRICT public path — DnaSequence.Complement()
///     (src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs lines 54–62),
///     which validates its input at construction (DnaSequence ctor +
///     ValidateSequence, lines 22–33 / 112–124):
///       • null or empty string  → an empty sequence; Complement() is the empty
///         sequence; no exception (string.IsNullOrEmpty short-circuit). A defined
///         result, NOT an error.
///       • input is case-folded with ToUpperInvariant before validation, so
///         lowercase / mixed case a-c-g-t is accepted and complements identically
///         to uppercase.
///       • ANY character that is not A/C/G/T after upper-casing (digits,
///         whitespace, N/IUPAC ambiguity codes, U, unicode letters, combining
///         marks, astral code points) → a *documented, intentional*
///         ArgumentException from ValidateSequence. The validation gate, not a
///         crash and not a silent mis-complement.
///
/// (2) The LENIENT char/span primitive — SequenceExtensions.GetComplementBase
///     (SequenceExtensions.cs lines 137–157) and the span overload
///     ReadOnlySpan&lt;char&gt;.TryGetComplement (lines 204–215). By separate
///     documented design (RNA_Complement.md §3.3, the DNA sibling) this surface:
///       • is IUPAC-complete: A↔T, C↔G, U→A, and the eleven ambiguity codes
///         (R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N);
///       • is case-insensitive and always emits UPPERCASE for recognized symbols;
///       • passes ANY non-IUPAC character (gaps, digits, whitespace, '\0',
///         unicode letters, surrogate halves) THROUGH UNCHANGED and NEVER throws;
///       • on an empty span succeeds and writes nothing — no exception, no hang.
///     We pin this so the lenient/strict boundary is explicit: the same garbage
///     that the public DnaSequence path REJECTS, the primitive must carry through
///     without crashing.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class CompositionFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>Deterministic RNG — seed fixed so generated fuzz inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260620);

    /// <summary>Generates a random valid DNA string of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(bases.Length)];
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-GC-001 — GC content : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region BE — Boundary: empty string

    /// <summary>
    /// BE: the empty string is the lower size boundary. The public DnaSequence
    /// path defines this as an empty sequence whose GC% is 0 (no division by
    /// zero, no exception) — DnaSequence.cs lines 24–28; GC% definition has an
    /// empty denominator, resolved to 0 by the repository zero-division
    /// convention (GC_Content_Profile.md INV-05).
    /// </summary>
    [Test]
    public void GcContent_EmptyString_IsZeroAndDoesNotThrow()
    {
        var act = () =>
        {
            var seq = new DnaSequence(string.Empty);
            seq.Length.Should().Be(0);
            seq.GcContent().Should().Be(0.0,
                because: "an empty sequence has no bases; GC% is defined as 0 by the zero-division convention");
        };

        act.Should().NotThrow("the empty string is a defined boundary input, not an error");
    }

    #endregion

    #region BE — Boundary: single character

    /// <summary>
    /// BE: a one-base sequence is the minimal non-empty input. GC% is the binary
    /// extreme — exactly 100 for a single G/C, exactly 0 for a single A/T — with
    /// no rounding drift. Lowercase is accepted (case-folded) and yields the same
    /// value. Verified over all four bases in both cases.
    /// </summary>
    [TestCase('G', 100.0)]
    [TestCase('C', 100.0)]
    [TestCase('A', 0.0)]
    [TestCase('T', 0.0)]
    [TestCase('g', 100.0)]
    [TestCase('c', 100.0)]
    [TestCase('a', 0.0)]
    [TestCase('t', 0.0)]
    public void GcContent_SingleCharacter_IsBinaryExtreme(char baseChar, double expectedGc)
    {
        var seq = new DnaSequence(baseChar.ToString());

        seq.GcContent().Should().BeApproximately(expectedGc, Tolerance,
            because: $"a single '{baseChar}' is {(expectedGc == 100.0 ? "a GC base → 100%" : "an AT base → 0%")}");
    }

    #endregion

    #region INJ — Injection: non-ACGT characters

    /// <summary>
    /// INJ: characters that are not A/C/G/T after upper-casing must be rejected
    /// by the public DnaSequence path with the *documented, intentional*
    /// ArgumentException (DnaSequence.ValidateSequence, lines 112–124). This is
    /// the validation contract — NOT a silent miscount and NOT a raw runtime
    /// exception. Covers digits, whitespace, punctuation, the ambiguity code N,
    /// the RNA base U (DNA does not accept U), and an embedded null byte.
    /// </summary>
    [TestCase("N", TestName = "GcContent_NonAcgt_AmbiguityCodeN_Throws")]
    [TestCase("ACGTN", TestName = "GcContent_NonAcgt_TrailingN_Throws")]
    [TestCase("ACGU", TestName = "GcContent_NonAcgt_RnaBaseU_Throws")]
    [TestCase("ACGT123", TestName = "GcContent_NonAcgt_Digits_Throws")]
    [TestCase("AC GT", TestName = "GcContent_NonAcgt_Whitespace_Throws")]
    [TestCase("ACGT-", TestName = "GcContent_NonAcgt_Punctuation_Throws")]
    [TestCase("ACG\0T", TestName = "GcContent_NonAcgt_EmbeddedNullByte_Throws")]
    public void GcContent_NonAcgtCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input);

        act.Should().Throw<ArgumentException>(
            "non-ACGT input is rejected at construction by the documented validation gate, " +
            "not miscounted and not crashed");
    }

    #endregion

    #region INJ / BE — Injection: null

    /// <summary>
    /// INJ/BE: a null reference is the boundary of "no input". The public
    /// DnaSequence path defines null as an empty sequence (string.IsNullOrEmpty
    /// short-circuit, DnaSequence.cs lines 24–28) — so it must NOT throw
    /// NullReferenceException and GcContent() must be the defined 0. This pins
    /// that null is handled gracefully rather than crashing.
    /// </summary>
    [Test]
    public void GcContent_NullSequence_IsTreatedAsEmptyAndDoesNotThrow()
    {
        var act = () =>
        {
            var seq = new DnaSequence(null!);
            seq.Length.Should().Be(0);
            seq.GcContent().Should().Be(0.0);
        };

        act.Should().NotThrow<NullReferenceException>(
            "null must be handled by the documented IsNullOrEmpty gate, never dereferenced");
        act.Should().NotThrow(
            "null is a defined 'empty sequence' input on the public path, not an error");
    }

    #endregion

    #region INJ — Injection: unicode

    /// <summary>
    /// INJ: unicode injection — non-ASCII letters, combining diacritics,
    /// full-width look-alikes, and astral/surrogate-pair code points. None of
    /// these are A/C/G/T, so the public DnaSequence path must reject every one
    /// with the documented ArgumentException — never an
    /// IndexOutOfRange/encoding surprise from surrogate handling. The astral case
    /// (😀, a surrogate pair) specifically guards char-by-char validation against
    /// crashing on the high/low surrogate halves.
    /// </summary>
    [TestCase("ÀCGT", TestName = "GcContent_Unicode_AccentedLatin_Throws")]
    [TestCase("ACGTα", TestName = "GcContent_Unicode_GreekLetter_Throws")]
    [TestCase("ÁCGT", TestName = "GcContent_Unicode_CombiningAcute_Throws")]
    [TestCase("ＡＣＧＴ", TestName = "GcContent_Unicode_FullWidthLatin_Throws")]
    [TestCase("ACG😀T", TestName = "GcContent_Unicode_AstralSurrogatePair_Throws")]
    public void GcContent_UnicodeCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input);

        act.Should().Throw<ArgumentException>(
            "unicode characters are not valid nucleotides; the validation gate must reject them " +
            "via ArgumentException, including surrogate-pair (astral) code points");
    }

    #endregion

    #region BE — Boundary: extremely long

    /// <summary>
    /// BE/OVF: an extremely long valid sequence (1,000,000 bases) must compute
    /// without overflow, hang, or precision blow-up, and the result must stay in
    /// the closed range [0, 100]. The denominator is an int count; this guards
    /// that the (double)gc / valid * 100 arithmetic does not overflow or drift
    /// out of range at scale. A known-composition long input pins the exact
    /// value too.
    /// </summary>
    [Test]
    public void GcContent_ExtremelyLong_StaysInRangeAndDoesNotHang()
    {
        const int length = 1_000_000;

        // Known composition: exactly half GC — "AG" repeated is 50% A (AT base)
        // and 50% G (GC base), so GC% is precisely 50.
        var halfGc = new DnaSequence(string.Concat(Enumerable.Repeat("AG", length / 2)));
        halfGc.Length.Should().Be(length);
        halfGc.GcContent().Should().BeApproximately(50.0, Tolerance,
            because: "a sequence that is exactly half G/C has GC% = 50, at any length");

        // Fixed-seed random long sequence: result must be a valid percentage.
        var random = new DnaSequence(RandomDna(length));
        random.GcContent().Should().BeInRange(0.0, 100.0,
            because: "GC% is a proportion times 100; it can never escape [0, 100], even at scale");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-GC-001 — backing primitive SequenceExtensions.CalculateGcContent
    //  (documented as *lenient*: excludes non-A/T/G/C/U from num & denom)
    // ═══════════════════════════════════════════════════════════════════

    #region BE/INJ — backing span primitive contract

    /// <summary>
    /// The lenient backing primitive must, by its documented contract
    /// (SequenceExtensions.cs lines 17–58), return 0 — not NaN, not a crash — for
    /// the empty span and for input that contains no valid nucleotide at all
    /// (avoiding 0/0). Pins the zero-division convention at the primitive layer.
    /// </summary>
    [TestCase("", TestName = "CalculateGcContent_Empty_IsZero")]
    [TestCase("N", TestName = "CalculateGcContent_OnlyAmbiguity_IsZero")]
    [TestCase("12345", TestName = "CalculateGcContent_OnlyDigits_IsZero")]
    [TestCase("----", TestName = "CalculateGcContent_OnlyPunctuation_IsZero")]
    public void CalculateGcContent_NoValidNucleotide_IsZero(string input)
    {
        double result = input.AsSpan().CalculateGcContent();

        result.Should().Be(0.0,
            because: "the lenient primitive excludes invalid symbols and returns 0 when the denominator is empty");
    }

    /// <summary>
    /// The lenient primitive must exclude non-A/T/G/C/U symbols from BOTH the
    /// numerator and the denominator (Biopython "remove" mode). So "GC" with any
    /// amount of injected garbage (N, digits, whitespace, null byte) interspersed
    /// must still read as 100% GC, and the canonical Biopython case
    /// gc_fraction("ACTGN") = 0.50 → 50% must hold. This pins that injection does
    /// not silently corrupt the count at the primitive layer.
    /// </summary>
    [TestCase("ACTGN", 50.0, TestName = "CalculateGcContent_BiopythonAcgtnRemoveMode_Is50")]
    [TestCase("G N C ", 100.0, TestName = "CalculateGcContent_GarbageAroundGc_Is100")]
    [TestCase("AT12", 0.0, TestName = "CalculateGcContent_DigitsAroundAt_Is0")]
    [TestCase("g\0c", 100.0, TestName = "CalculateGcContent_NullByteBetweenGc_Is100")]
    [TestCase("acgu", 50.0, TestName = "CalculateGcContent_LowercaseWithU_Is50")]
    public void CalculateGcContent_ExcludesInvalidSymbols(string input, double expected)
    {
        double result = input.AsSpan().CalculateGcContent();

        result.Should().BeApproximately(expected, Tolerance,
            because: "non-A/T/G/C/U symbols are excluded from both numerator and denominator (remove mode)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-COMP-001 — DNA complement
    //  Strict public path: DnaSequence.Complement() (validates at construction)
    //  Lenient primitive:  SequenceExtensions.GetComplementBase / TryGetComplement
    // ═══════════════════════════════════════════════════════════════════

    #region SEQ-COMP-001 — DNA complement

    #region INJ — Injection: non-DNA characters (strict path rejects)

    /// <summary>
    /// INJ: characters that are not A/C/G/T after upper-casing must be rejected by
    /// the strict public DnaSequence path with the *documented, intentional*
    /// ArgumentException (DnaSequence ctor → ValidateSequence, lines 22–33 /
    /// 112–124) BEFORE any complement is taken — never a silent mis-complement and
    /// never a raw runtime exception. Covers digits, whitespace, punctuation/gap,
    /// the ambiguity code N, the RNA base U (DNA does not accept U), and an
    /// embedded null byte.
    /// </summary>
    [TestCase("N", TestName = "Complement_NonDna_AmbiguityCodeN_Throws")]
    [TestCase("ACGTN", TestName = "Complement_NonDna_TrailingN_Throws")]
    [TestCase("ACGU", TestName = "Complement_NonDna_RnaBaseU_Throws")]
    [TestCase("ACGT123", TestName = "Complement_NonDna_Digits_Throws")]
    [TestCase("AC GT", TestName = "Complement_NonDna_Whitespace_Throws")]
    [TestCase("ACGT-", TestName = "Complement_NonDna_GapDash_Throws")]
    [TestCase("ACG\0T", TestName = "Complement_NonDna_EmbeddedNullByte_Throws")]
    public void Complement_NonDnaCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input).Complement();

        act.Should().Throw<ArgumentException>(
            "non-A/C/G/T input is rejected at construction by the documented validation gate, " +
            "so the complement is never computed on garbage and never crashes");
    }

    #endregion

    #region BE — Boundary: empty string (strict path)

    /// <summary>
    /// BE: the empty string is the lower size boundary. The strict DnaSequence
    /// path defines it as an empty sequence (DnaSequence.cs lines 24–28); its
    /// complement is therefore the empty sequence — no division, no indexing, no
    /// exception. Complement of nothing is nothing.
    /// </summary>
    [Test]
    public void Complement_EmptyString_IsEmptyAndDoesNotThrow()
    {
        var act = () =>
        {
            var complement = new DnaSequence(string.Empty).Complement();
            complement.Length.Should().Be(0);
            complement.Sequence.Should().BeEmpty(
                because: "the complement of an empty sequence is the empty sequence");
        };

        act.Should().NotThrow("the empty string is a defined boundary input, not an error");
    }

    #endregion

    #region INJ / BE — Injection: null (strict path treats as empty)

    /// <summary>
    /// INJ/BE: a null reference is the boundary of "no input". The strict
    /// DnaSequence path defines null as an empty sequence (string.IsNullOrEmpty
    /// short-circuit, DnaSequence.cs lines 24–28), so Complement() must NOT throw
    /// NullReferenceException and must yield the empty sequence. Pins that null is
    /// handled gracefully rather than dereferenced.
    /// </summary>
    [Test]
    public void Complement_NullSequence_IsTreatedAsEmptyAndDoesNotThrow()
    {
        var act = () =>
        {
            var complement = new DnaSequence(null!).Complement();
            complement.Length.Should().Be(0);
            complement.Sequence.Should().BeEmpty();
        };

        act.Should().NotThrow<NullReferenceException>(
            "null must be handled by the documented IsNullOrEmpty gate, never dereferenced");
        act.Should().NotThrow(
            "null is a defined 'empty sequence' input on the public path, not an error");
    }

    #endregion

    #region INJ — Injection: mixed case (strict path accepts, case-folded)

    /// <summary>
    /// INJ: mixed and lower case input is accepted by the strict path — it is
    /// upper-cased before validation (ToUpperInvariant, DnaSequence.cs line 30) —
    /// and must complement IDENTICALLY to the uppercase form, always emitting
    /// uppercase A/T/G/C. This guards that case-folding neither rejects valid DNA
    /// nor corrupts the A↔T / C↔G mapping.
    /// </summary>
    [TestCase("acgt", "ACGT", TestName = "Complement_MixedCase_AllLower_FoldsAndComplements")]
    [TestCase("AcGt", "ACGT", TestName = "Complement_MixedCase_Alternating_FoldsAndComplements")]
    [TestCase("aCgT", "ACGT", TestName = "Complement_MixedCase_AlternatingInverse_FoldsAndComplements")]
    public void Complement_MixedCase_FoldsToUppercaseAndComplements(string input, string upper)
    {
        // The complement of A/C/G/T is T/G/C/A; computed once on the canonical
        // uppercase form, it must equal the complement of the mixed-case input.
        string expected = new DnaSequence(upper).Complement().Sequence;

        var complement = new DnaSequence(input).Complement();

        complement.Sequence.Should().Be(expected,
            because: "input is case-folded before complementing; case must not change the result");
        complement.Sequence.Should().MatchRegex("^[ACGT]*$",
            because: "the complement always emits uppercase canonical bases");
    }

    #endregion

    #region INJ — Injection: unicode (strict path rejects)

    /// <summary>
    /// INJ: unicode injection — accented Latin, Greek letters, combining
    /// diacritics, full-width look-alikes, and astral/surrogate-pair code points.
    /// None are A/C/G/T, so the strict DnaSequence path must reject every one with
    /// the documented ArgumentException — never an IndexOutOfRange/encoding
    /// surprise from surrogate handling. The astral case (😀, a surrogate pair)
    /// specifically guards char-by-char validation against crashing on the
    /// high/low surrogate halves before the complement is ever taken.
    /// </summary>
    [TestCase("ÀCGT", TestName = "Complement_Unicode_AccentedLatin_Throws")]
    [TestCase("ACGTα", TestName = "Complement_Unicode_GreekLetter_Throws")]
    [TestCase("ÁCGT", TestName = "Complement_Unicode_CombiningAcute_Throws")]
    [TestCase("ＡＣＧＴ", TestName = "Complement_Unicode_FullWidthLatin_Throws")]
    [TestCase("ACG😀T", TestName = "Complement_Unicode_AstralSurrogatePair_Throws")]
    public void Complement_UnicodeCharacters_ThrowDocumentedArgumentException(string input)
    {
        var act = () => _ = new DnaSequence(input).Complement();

        act.Should().Throw<ArgumentException>(
            "unicode characters are not valid nucleotides; the validation gate must reject them " +
            "via ArgumentException, including surrogate-pair (astral) code points");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────
    //  Lenient primitive: GetComplementBase / span TryGetComplement
    //  (IUPAC-complete; non-IUPAC passes through unchanged; never throws)
    // ───────────────────────────────────────────────────────────────────

    #region INJ — lenient primitive: non-DNA / unicode pass through, never throw

    /// <summary>
    /// The lenient char primitive must NEVER throw on any char and, by its
    /// documented contract (RNA_Complement.md §3.3, the DNA sibling), must pass
    /// any non-IUPAC character through UNCHANGED. Recognized IUPAC symbols are
    /// complemented and emitted uppercase (A↔T, C↔G, U→A, N↔N, R↔Y). This pins
    /// that injection (digits, gap, whitespace, null byte, unicode letters,
    /// surrogate halves) cannot crash or silently corrupt the recognized mapping.
    /// </summary>
    [TestCase('A', 'T', TestName = "GetComplementBase_A_IsT")]
    [TestCase('C', 'G', TestName = "GetComplementBase_C_IsG")]
    [TestCase('U', 'A', TestName = "GetComplementBase_U_IsA")]
    [TestCase('N', 'N', TestName = "GetComplementBase_AmbiguityN_IsN")]
    [TestCase('R', 'Y', TestName = "GetComplementBase_AmbiguityR_IsY")]
    [TestCase('1', '1', TestName = "GetComplementBase_Digit_PassesThrough")]
    [TestCase('-', '-', TestName = "GetComplementBase_GapDash_PassesThrough")]
    [TestCase(' ', ' ', TestName = "GetComplementBase_Whitespace_PassesThrough")]
    [TestCase('\0', '\0', TestName = "GetComplementBase_NullByte_PassesThrough")]
    [TestCase('α', 'α', TestName = "GetComplementBase_GreekLetter_PassesThrough")]
    [TestCase('Z', 'Z', TestName = "GetComplementBase_NonIupacLetter_PassesThrough")]
    [TestCase('\uD83D', '\uD83D', TestName = "GetComplementBase_HighSurrogateHalf_PassesThrough")]
    [TestCase('\uDE00', '\uDE00', TestName = "GetComplementBase_LowSurrogateHalf_PassesThrough")]
    public void GetComplementBase_AnyChar_NeverThrowsAndPassesNonIupacThrough(char input, char expected)
    {
        char result = '￿';
        var act = () => result = SequenceExtensions.GetComplementBase(input);

        act.Should().NotThrow("the lenient primitive is total over char and never throws");
        result.Should().Be(expected,
            because: "recognized IUPAC symbols are complemented (uppercase); everything else passes through unchanged");
    }

    /// <summary>
    /// The lenient span complement must carry a mix of recognized and non-IUPAC
    /// garbage through char-by-char without throwing: recognized bases are
    /// complemented and uppercased, every other character (gap, digit, null byte)
    /// is preserved verbatim. This pins that injected garbage interspersed with
    /// real bases neither crashes nor shifts the complement of the valid bases.
    /// </summary>
    [Test]
    public void TryGetComplement_GarbageInterspersed_ComplementsBasesAndPreservesGarbage()
    {
        const string input = "aC-G1N\0T";       // mixed case, gap, digit, N, null byte
        Span<char> destination = new char[input.Length];

        bool ok = input.AsSpan().TryGetComplement(destination);

        ok.Should().BeTrue("the destination is exactly the source length");
        new string(destination).Should().Be("TG-C1N\0A",
            because: "A→T C→G G→C T→A N→N (uppercase), and non-IUPAC '-','1','\\0' pass through unchanged");
    }

    #endregion

    #region BE — lenient primitive: empty span

    /// <summary>
    /// BE: the empty span is the lower size boundary for the lenient primitive.
    /// TryGetComplement must succeed (destination length ≥ source length holds
    /// trivially) and write nothing — no exception, no hang, no out-of-range.
    /// </summary>
    [Test]
    public void TryGetComplement_EmptySpan_SucceedsAndWritesNothing()
    {
        var act = () =>
        {
            bool ok = ReadOnlySpan<char>.Empty.TryGetComplement(Span<char>.Empty);
            ok.Should().BeTrue("an empty complement always fits an empty destination");
        };

        act.Should().NotThrow("the empty span is a defined boundary, not an error");
    }

    #endregion

    #endregion
}
