using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Quality area (Phred char ↔ score encoding).
///
/// Algebraic testing pins the encode∘decode round-trip isomorphism between Phred
/// scores and their ASCII quality characters, and the per-character definition.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, row 219.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Quality")]
public class QualityAlgebraicTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: QUALITY-PHRED-001 — Phred quality encoding (Quality), row 219.
    //
    // Model: Phred+33 maps a score q to the character (char)(q+33) and back via
    //        c−33 — a bijection between scores [0,93] and printable ASCII [33,126].
    //   — docs/algorithms/Quality; QualityScoreAnalyzer.CharToPhred / PhredToChar.
    //
    // Laws (row 219): RT — PhredToChar(CharToPhred(c)) = c and
    //                 CharToPhred(PhredToChar(q)) = q over the valid domains.
    //                 ID — defined per char: CharToPhred(c) = c − 33.
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property Phred_RoundTrip_ScoreCharScoreIsIdentity()
    {
        return Prop.ForAll(Gen.Choose(0, 93).ToArbitrary(), q =>
        {
            char c = QualityScoreAnalyzer.PhredToChar(q);
            int back = QualityScoreAnalyzer.CharToPhred(c);
            return (back == q).Label($"score {q} -> '{c}' -> {back}");
        });
    }

    [FsCheck.NUnit.Property]
    public Property Phred_RoundTrip_CharScoreCharIsIdentity()
    {
        return Prop.ForAll(Gen.Choose(33, 126).Select(i => (char)i).ToArbitrary(), c =>
        {
            int q = QualityScoreAnalyzer.CharToPhred(c);
            char back = QualityScoreAnalyzer.PhredToChar(q);
            return (back == c).Label($"'{c}' -> {q} -> '{back}'");
        });
    }

    [FsCheck.NUnit.Property]
    public Property Phred_Identity_DefinedPerChar()
    {
        return Prop.ForAll(Gen.Choose(33, 126).Select(i => (char)i).ToArbitrary(), c =>
            (QualityScoreAnalyzer.CharToPhred(c) == c - 33).Label($"CharToPhred('{c}')"));
    }
}
