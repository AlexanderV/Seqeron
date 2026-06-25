namespace Seqeron.Genomics.Tests.Combinatorial;

using Seqeron.Genomics.IO;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Quality area (QualityScoreAnalyzer,
/// Seqeron.Genomics.IO).
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of combinatorial testing.
/// Each grid cell carries a real business assertion; small grids use the exhaustive
/// <c>[Combinatorial]</c> product (a strict superset of pairwise).
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Quality")]
public class QualityCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: QUALITY-PHRED-001 — Phred quality encode/decode (Quality)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 219.
    // Spec: tests/TestSpecs/QUALITY-PHRED-001.md (canonical ParseQualityString / ToQualityString /
    //       ConvertEncoding). ADVANCED §10.
    // Dimensions: encoding(3) × seqLen(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Cock 2010): Q = ord(char) − offset (33 for Phred+33, 64 for Phred+64); encode is the
    // inverse; converting between encodings re-offsets the characters while preserving the Phred score.
    //
    // Axis mapping (documented): encoding → the round-trip mode {Phred33, Phred64, Cross-convert};
    // seqLen → number of scores. The combinatorial point: encode∘decode is the identity for each
    // encoding, and a Phred33→Phred64 conversion preserves the decoded Phred scores, at every length.
    // ═══════════════════════════════════════════════════════════════════════

    public enum QualMode { Phred33RoundTrip, Phred64RoundTrip, CrossConvert }

    [Test, Combinatorial]
    public void QualityPhred_RoundTripAndConversion_AcrossModeAndLength(
        [Values(QualMode.Phred33RoundTrip, QualMode.Phred64RoundTrip, QualMode.CrossConvert)] QualMode mode,
        [Values(4, 16, 40)] int seqLen)
    {
        // Deterministic scores in [0,40] (valid under both Phred+33 and Phred+64).
        int[] scores = Enumerable.Range(0, seqLen).Select(i => (i * 7) % 41).ToArray();

        switch (mode)
        {
            case QualMode.Phred33RoundTrip:
                QualityScoreAnalyzer.ParseQualityString(
                    QualityScoreAnalyzer.ToQualityString(scores, QualityScoreAnalyzer.QualityEncoding.Phred33), QualityScoreAnalyzer.QualityEncoding.Phred33)
                    .Should().Equal(scores, "encode∘decode is the identity under Phred+33");
                break;
            case QualMode.Phred64RoundTrip:
                QualityScoreAnalyzer.ParseQualityString(
                    QualityScoreAnalyzer.ToQualityString(scores, QualityScoreAnalyzer.QualityEncoding.Phred64), QualityScoreAnalyzer.QualityEncoding.Phred64)
                    .Should().Equal(scores, "encode∘decode is the identity under Phred+64");
                break;
            default:
                string p33 = QualityScoreAnalyzer.ToQualityString(scores, QualityScoreAnalyzer.QualityEncoding.Phred33);
                string p64 = QualityScoreAnalyzer.ConvertEncoding(p33, QualityScoreAnalyzer.QualityEncoding.Phred33, QualityScoreAnalyzer.QualityEncoding.Phred64);
                QualityScoreAnalyzer.ParseQualityString(p64, QualityScoreAnalyzer.QualityEncoding.Phred64)
                    .Should().Equal(scores, "conversion preserves the decoded Phred scores");
                break;
        }
    }

    /// <summary>
    /// Interaction witness — the offset convention: '!' (ASCII 33) is Q0 in Phred+33, '@' (ASCII 64)
    /// is Q0 in Phred+64, and conversion shifts the character by the offset difference (31).
    /// </summary>
    [Test]
    public void QualityPhred_OffsetConvention()
    {
        QualityScoreAnalyzer.ParseQualityString("!", QualityScoreAnalyzer.QualityEncoding.Phred33).Should().Equal(new[] { 0 });
        QualityScoreAnalyzer.ParseQualityString("@", QualityScoreAnalyzer.QualityEncoding.Phred64).Should().Equal(new[] { 0 });

        // Q20: char '5' (53) in Phred+33, char 'T' (84) in Phred+64 — same score, char shifted by 31.
        QualityScoreAnalyzer.ToQualityString(new[] { 20 }, QualityScoreAnalyzer.QualityEncoding.Phred33).Should().Be("5");
        QualityScoreAnalyzer.ToQualityString(new[] { 20 }, QualityScoreAnalyzer.QualityEncoding.Phred64).Should().Be("T");
        QualityScoreAnalyzer.ConvertEncoding("5", QualityScoreAnalyzer.QualityEncoding.Phred33, QualityScoreAnalyzer.QualityEncoding.Phred64).Should().Be("T");
    }
}
