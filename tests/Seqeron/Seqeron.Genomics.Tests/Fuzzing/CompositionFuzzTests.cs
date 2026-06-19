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
}
