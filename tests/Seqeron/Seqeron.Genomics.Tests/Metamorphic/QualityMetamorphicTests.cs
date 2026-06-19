using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Quality area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: QUALITY-PHRED-001 — FASTQ Phred quality encode/decode (Quality).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 219.
///
/// API under test (QualityScoreAnalyzer.ToQualityString / ParseQualityString):
///   Encode: char = chr(Q + offset). Decode: Q = ord(char) − offset. Offset = 33 (Phred+33,
///   Sanger/Illumina 1.8+) or 64 (Phred+64, Illumina 1.3–1.7) — Cock et al. (2010).
///
/// Relations (derived from the offset bijection, NOT from output):
///   • RT  (encode∘decode = identity): for valid scores, decode(encode(s)) = s; and for a valid
///         quality string, encode(decode(q)) = q — the offset map is a bijection on the valid range.
///   • INV (offset consistency): the only thing distinguishing the two encodings is the additive
///         offset, so the SAME Phred score encodes to ASCII characters exactly (64 − 33) = 31 apart,
///         and the SAME character decodes to scores 31 apart.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class QualityMetamorphicTests
{
    private const int OffsetDifference = 64 - 33; // Phred+64 vs Phred+33

    // Scores valid in BOTH encodings (Phred+64 caps at 62), so every relation can be checked on both.
    private static readonly int[] SharedScores = { 0, 1, 2, 10, 20, 30, 40, 50, 62 };

    #region QUALITY-PHRED-001 RT — encode then decode is the identity

    [Test]
    [Description("RT: decode(encode(scores)) = scores under each encoding — the offset map is a bijection on the valid range.")]
    public void Phred_EncodeThenDecode_RecoversScores()
    {
        foreach (var encoding in new[] { QualityScoreAnalyzer.QualityEncoding.Phred33, QualityScoreAnalyzer.QualityEncoding.Phred64 })
        {
            string encoded = QualityScoreAnalyzer.ToQualityString(SharedScores, encoding);
            int[] decoded = QualityScoreAnalyzer.ParseQualityString(encoded, encoding);

            decoded.Should().Equal(SharedScores,
                because: $"chr(Q + offset) then ord(c) − offset recovers Q exactly under {encoding}");
        }
    }

    [Test]
    [Description("RT: encode(decode(q)) = q for a valid quality string — the round-trip is lossless on the character side too.")]
    public void Phred_DecodeThenEncode_RecoversQualityString()
    {
        // A varied Phred+33 quality line (ASCII 33..73 ⇒ Q0..Q40), all within range.
        const string quality = "!\"#$%5?@ABHI";

        int[] scores = QualityScoreAnalyzer.ParseQualityString(quality, QualityScoreAnalyzer.QualityEncoding.Phred33);
        string roundTripped = QualityScoreAnalyzer.ToQualityString(scores, QualityScoreAnalyzer.QualityEncoding.Phred33);

        roundTripped.Should().Be(quality, because: "ord(c) − offset then chr(Q + offset) restores the original characters");
    }

    #endregion

    #region QUALITY-PHRED-001 INV — the two encodings differ only by a constant offset

    [Test]
    [Description("INV: the same Phred score encodes to characters exactly 31 apart between Phred+64 and Phred+33.")]
    public void Phred_SameScore_EncodingsDifferByOffset()
    {
        string p33 = QualityScoreAnalyzer.ToQualityString(SharedScores, QualityScoreAnalyzer.QualityEncoding.Phred33);
        string p64 = QualityScoreAnalyzer.ToQualityString(SharedScores, QualityScoreAnalyzer.QualityEncoding.Phred64);

        p64.Length.Should().Be(p33.Length);
        for (int i = 0; i < p33.Length; i++)
            (p64[i] - p33[i]).Should().Be(OffsetDifference,
                because: $"score {SharedScores[i]} encodes 31 ASCII positions higher under Phred+64 than Phred+33");
    }

    [Test]
    [Description("INV: the same character decodes to scores exactly 31 apart between Phred+33 and Phred+64.")]
    public void Phred_SameCharacter_DecodingsDifferByOffset()
    {
        // Characters in ASCII 64..126 are valid under both encodings (decode ≥ 0 in each).
        string shared = new string(Enumerable.Range(64, 30).Select(a => (char)a).ToArray());

        int[] asP33 = QualityScoreAnalyzer.ParseQualityString(shared, QualityScoreAnalyzer.QualityEncoding.Phred33);
        int[] asP64 = QualityScoreAnalyzer.ParseQualityString(shared, QualityScoreAnalyzer.QualityEncoding.Phred64);

        for (int i = 0; i < shared.Length; i++)
            (asP33[i] - asP64[i]).Should().Be(OffsetDifference,
                because: "decoding subtracts a 31-larger offset under Phred+64, so its scores are 31 lower for the same character");
    }

    #endregion
}
