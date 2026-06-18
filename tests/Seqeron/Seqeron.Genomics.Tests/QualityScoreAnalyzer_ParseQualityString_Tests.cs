// QUALITY-PHRED-001 — Phred Score Handling
// Evidence: docs/Evidence/QUALITY-PHRED-001-Evidence.md
// TestSpec: tests/TestSpecs/QUALITY-PHRED-001.md
// Source: Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM (2010). Nucleic Acids Research 38(6):1767-1771.

using NUnit.Framework;
using Seqeron.Genomics.IO;
using System;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class QualityScoreAnalyzer_ParseQualityString_Tests
{
    private const QualityScoreAnalyzer.QualityEncoding Phred33 = QualityScoreAnalyzer.QualityEncoding.Phred33;
    private const QualityScoreAnalyzer.QualityEncoding Phred64 = QualityScoreAnalyzer.QualityEncoding.Phred64;

    #region ParseQualityString

    // M1 — Cock et al. (2010): Phred+33 offset 33; '!'(33)->0, '~'(126)->93.
    [Test]
    public void ParseQualityString_Phred33Boundaries_DecodesZeroAnd93()
    {
        var scores = QualityScoreAnalyzer.ParseQualityString("!~", Phred33);

        Assert.That(scores, Is.EqualTo(new[] { 0, 93 }),
            "Phred+33 boundary chars '!'(ASCII 33) and '~'(ASCII 126) decode to Phred 0 and 93 (offset 33).");
    }

    // M2 — Cock et al. (2010): Phred+33 interior '5'(53)->20, '?'(63)->30, 'I'(73)->40.
    [Test]
    public void ParseQualityString_Phred33Interior_DecodesExactScores()
    {
        var scores = QualityScoreAnalyzer.ParseQualityString("5?I", Phred33);

        Assert.That(scores, Is.EqualTo(new[] { 20, 30, 40 }),
            "Phred+33 chars '5','?','I' decode to Phred 20,30,40 (ASCII - 33).");
    }

    // M3 — Cock et al. (2010): Phred+64 offset 64; '@'(64)->0, 'h'(104)->40, '~'(126)->62.
    [Test]
    public void ParseQualityString_Phred64_DecodesBoundaryAndInterior()
    {
        var scores = QualityScoreAnalyzer.ParseQualityString("@h~", Phred64);

        Assert.That(scores, Is.EqualTo(new[] { 0, 40, 62 }),
            "Phred+64 chars '@','h','~' decode to Phred 0,40,62 (offset 64; max valid 62).");
    }

    // M8 — Cock et al. (2010) corner case: a char below the offset decodes to negative Phred (invalid).
    [Test]
    public void ParseQualityString_CharBelowOffset_ThrowsOutOfRange()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => QualityScoreAnalyzer.ParseQualityString(" ", Phred33),
            "Space (ASCII 32) decodes to Phred -1 under Phred+33, outside valid [0,93]; must throw.");
    }

    // S1 — empty input is the identity boundary.
    [Test]
    public void ParseQualityString_Empty_ReturnsEmptyArray()
    {
        var scores = QualityScoreAnalyzer.ParseQualityString("", Phred33);

        Assert.That(scores, Is.Empty, "Empty quality string parses to an empty score array.");
    }

    // S3 — public-API failure mode.
    [Test]
    public void ParseQualityString_Null_ThrowsArgumentNull()
    {
        Assert.Throws<ArgumentNullException>(() => QualityScoreAnalyzer.ParseQualityString(null!, Phred33),
            "Null quality string must throw ArgumentNullException.");
    }

    #endregion

    #region ToQualityString

    // M4 — Cock et al. (2010): encode char = chr(Q + 33). [0,20,30,40,93] -> "!5?I~".
    [Test]
    public void ToQualityString_Phred33_EncodesExactCharacters()
    {
        string s = QualityScoreAnalyzer.ToQualityString(new[] { 0, 20, 30, 40, 93 }, Phred33);

        Assert.That(s, Is.EqualTo("!5?I~"),
            "Phred 0,20,30,40,93 encode to '!','5','?','I','~' under Phred+33 (Q + 33).");
    }

    // M5 — Cock et al. (2010): Phred+64 [0,40,62] -> "@h~".
    [Test]
    public void ToQualityString_Phred64_EncodesExactCharacters()
    {
        string s = QualityScoreAnalyzer.ToQualityString(new[] { 0, 40, 62 }, Phred64);

        Assert.That(s, Is.EqualTo("@h~"),
            "Phred 0,40,62 encode to '@','h','~' under Phred+64 (Q + 64).");
    }

    // M10 — Cock et al. (2010): Phred+33 max valid score is 93; 94 is out of range.
    [Test]
    public void ToQualityString_ScoreAbovePhred33Max_ThrowsOutOfRange()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => QualityScoreAnalyzer.ToQualityString(new[] { 94 }, Phred33),
            "Phred 94 exceeds the Phred+33 maximum of 93; must throw.");
    }

    // S2 — empty input is the identity boundary.
    [Test]
    public void ToQualityString_Empty_ReturnsEmptyString()
    {
        string s = QualityScoreAnalyzer.ToQualityString(Array.Empty<int>(), Phred33);

        Assert.That(s, Is.Empty, "Empty score array encodes to an empty string.");
    }

    // S4 — public-API failure mode.
    [Test]
    public void ToQualityString_Null_ThrowsArgumentNull()
    {
        Assert.Throws<ArgumentNullException>(() => QualityScoreAnalyzer.ToQualityString(null!, Phred33),
            "Null scores must throw ArgumentNullException.");
    }

    #endregion

    #region ConvertEncoding

    // M6 — Cock et al. (2010): Phred score invariant; Phred+64 "@h~"(Q 0,40,62) -> Phred+33 "!I_".
    [Test]
    public void ConvertEncoding_Phred64ToPhred33_PreservesScores()
    {
        string converted = QualityScoreAnalyzer.ConvertEncoding("@h~", Phred64, Phred33);

        Assert.Multiple(() =>
        {
            Assert.That(converted, Is.EqualTo("!I_"),
                "Phred+64 '@h~' (Q 0,40,62) re-offsets to Phred+33 '!I_' (same scores).");
            Assert.That(QualityScoreAnalyzer.ParseQualityString(converted, Phred33),
                Is.EqualTo(new[] { 0, 40, 62 }),
                "Decoding the converted Phred+33 string yields the original Phred scores 0,40,62.");
        });
    }

    // M7 — Cock et al. (2010): Phred+33 "!I"(Q 0,40) -> Phred+64 "@h".
    [Test]
    public void ConvertEncoding_Phred33ToPhred64_PreservesScores()
    {
        string converted = QualityScoreAnalyzer.ConvertEncoding("!I", Phred33, Phred64);

        Assert.That(converted, Is.EqualTo("@h"),
            "Phred+33 '!I' (Q 0,40) re-offsets to Phred+64 '@h' (same scores).");
    }

    // M9 — Cock et al. (2010) corner case: Phred+33 Q>62 has no Phred+64 representation; '~'=93.
    [Test]
    public void ConvertEncoding_Phred33ToPhred64_OverflowThrows()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => QualityScoreAnalyzer.ConvertEncoding("~", Phred33, Phred64),
            "Phred+33 '~' (Q 93) exceeds the Phred+64 maximum of 62 and cannot be converted.");
    }

    // S5 — public-API failure mode.
    [Test]
    public void ConvertEncoding_Null_ThrowsArgumentNull()
    {
        Assert.Throws<ArgumentNullException>(() => QualityScoreAnalyzer.ConvertEncoding(null!, Phred33, Phred64),
            "Null quality string must throw ArgumentNullException.");
    }

    // S6 — converting to the same encoding is a no-op.
    [Test]
    public void ConvertEncoding_SameEncoding_ReturnsIdentity()
    {
        string converted = QualityScoreAnalyzer.ConvertEncoding("5?I", Phred33, Phred33);

        Assert.That(converted, Is.EqualTo("5?I"),
            "Phred+33 -> Phred+33 conversion is an identity re-offset.");
    }

    #endregion

    #region Properties (INV-03)

    // C1 — INV-03: ToQualityString(ParseQualityString(s)) == s for any valid Phred+33 string (fixed seed).
    [Test]
    public void ParseThenEncode_Phred33_RoundTripsToIdentity()
    {
        // Build a valid Phred+33 string from a fixed seed over the full valid ASCII range [33,126].
        var rng = new Random(20260613);
        var chars = new char[256];
        for (int i = 0; i < chars.Length; i++)
            chars[i] = (char)rng.Next(33, 127); // inclusive 33..126 = valid Phred+33 ASCII
        string original = new string(chars);

        int[] scores = QualityScoreAnalyzer.ParseQualityString(original, Phred33);
        string roundTrip = QualityScoreAnalyzer.ToQualityString(scores, Phred33);

        Assert.That(roundTrip, Is.EqualTo(original),
            "Encode after decode reconstructs the original Phred+33 quality string (INV-03).");
    }

    #endregion
}
