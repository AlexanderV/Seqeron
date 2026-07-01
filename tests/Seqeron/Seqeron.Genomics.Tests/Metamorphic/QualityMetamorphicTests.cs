using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests.Metamorphic;

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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: QUALITY-STATS-001 — quality-score statistics (Quality).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 220.
    //
    // API under test (QualityScoreAnalyzer.CalculateStatistics(string)):
    //   Decodes the quality line to Phred scores and reports MeanQuality, Min/Max, StandardDeviation,
    //   TotalBases, and BasesAboveQ20/Q30 (inclusive thresholds).
    //
    // Relations (derived from the aggregate definitions, NOT from output):
    //   • INV (order independent for mean): the mean, extrema, spread and threshold COUNTS depend only
    //         on the multiset of scores, so permuting the quality string leaves them unchanged.
    //   • ADD (counts additive): TotalBases and BasesAboveQ20/Q30 — and the total score sum
    //         (mean·total) — are additive over concatenation of quality strings.
    // ───────────────────────────────────────────────────────────────────────────

    // A varied Phred+33 line: '!'=Q0, '+'=Q10, '5'=Q20, '?'=Q30, 'I'=Q40 — a mix above and below
    // the Q20/Q30 thresholds so the counts are non-trivial.
    private const string QualityA = "!+5?I?5+!I?5";
    private const string QualityB = "I?5+!!+5?III";

    private static QualityScoreAnalyzer.QualityStatistics Stats(string quality) =>
        QualityScoreAnalyzer.CalculateStatistics(quality);

    #region QUALITY-STATS-001 INV — mean/extrema/counts are independent of order

    [Test]
    [Description("INV: the mean, min, max, spread and threshold counts depend only on the score multiset, so permuting the quality string leaves them unchanged.")]
    public void Stats_Permutation_PreservesAggregateStatistics()
    {
        var original = Stats(QualityA);

        var chars = QualityA.ToCharArray();
        var rng = new Random(20260620);
        for (int i = chars.Length - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (chars[i], chars[j]) = (chars[j], chars[i]);
        }
        var permuted = Stats(new string(chars));

        permuted.MeanQuality.Should().BeApproximately(original.MeanQuality, 1e-9, because: "the mean depends only on the score multiset");
        permuted.MedianQuality.Should().Be(original.MedianQuality, because: "the median depends only on the sorted scores");
        permuted.MinQuality.Should().Be(original.MinQuality);
        permuted.MaxQuality.Should().Be(original.MaxQuality);
        permuted.StandardDeviation.Should().BeApproximately(original.StandardDeviation, 1e-9, because: "spread is order-independent");
        permuted.TotalBases.Should().Be(original.TotalBases);
        permuted.BasesAboveQ20.Should().Be(original.BasesAboveQ20, because: "the count of Q≥20 bases ignores their order");
        permuted.BasesAboveQ30.Should().Be(original.BasesAboveQ30, because: "the count of Q≥30 bases ignores their order");
    }

    #endregion

    #region QUALITY-STATS-001 ADD — base counts are additive over concatenation

    [Test]
    [Description("ADD: TotalBases and BasesAboveQ20/Q30 — and the total score sum — are additive over a concatenation of quality strings.")]
    public void Stats_Concatenation_CountsAreAdditive()
    {
        var a = Stats(QualityA);
        var b = Stats(QualityB);
        var ab = Stats(QualityA + QualityB);

        ab.TotalBases.Should().Be(a.TotalBases + b.TotalBases, because: "every base of both reads is counted once");
        ab.BasesAboveQ20.Should().Be(a.BasesAboveQ20 + b.BasesAboveQ20, because: "the Q≥20 count is a per-base tally, hence additive");
        ab.BasesAboveQ30.Should().Be(a.BasesAboveQ30 + b.BasesAboveQ30, because: "the Q≥30 count is a per-base tally, hence additive");

        // The score SUM (mean·total) is additive; the pooled mean is their length-weighted average.
        double sumA = a.MeanQuality * a.TotalBases;
        double sumB = b.MeanQuality * b.TotalBases;
        (ab.MeanQuality * ab.TotalBases).Should().BeApproximately(sumA + sumB, 1e-9,
            because: "concatenation pools the scores, so the total sum is additive");

        ab.MaxQuality.Should().Be(Math.Max(a.MaxQuality, b.MaxQuality), because: "the pooled max is the max of the parts");
        ab.MinQuality.Should().Be(Math.Min(a.MinQuality, b.MinQuality), because: "the pooled min is the min of the parts");
    }

    #endregion
}
