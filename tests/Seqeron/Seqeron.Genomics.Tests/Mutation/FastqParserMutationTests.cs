using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.IO;
using static Seqeron.Genomics.IO.FastqParser;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for FastqParser.cs (checklist 04 row 65, PARSE-FASTQ-001).
/// The canonical suite left the quality-encoding boundaries (Phred+33/+64 detection per the
/// Cock et al. 2010 FASTQ specification), the length/quality filters, end-trimming, adapter
/// trimming and the summary statistics under-pinned. These pin the exact published behaviour
/// so the boundary/arithmetic mutants diverge.
/// </summary>
[TestFixture]
public class FastqParserMutationTests
{
    // ── DetectEncoding: Phred+33 uses '!'(33)..'I'(73); Phred+64 uses '@'(64)..'h'(104) ──

    [Test]
    public void DetectEncoding_CharBelowAt_IsPhred33()
    {
        // '5' (53) < '@' (64) ⇒ Phred33 (kills `c < '@'` → `c <= '@'`/`>`).
        DetectEncoding("5555").Should().Be(QualityEncoding.Phred33);
    }

    [Test]
    public void DetectEncoding_AtThenAboveI_IsPhred64()
    {
        // '@'(64) is in the ambiguous overlap so parsing continues; 'J'(74) > 'I'(73) ⇒ Phred64.
        // Kills `c < '@'` → `c <= '@'` (which would early-return Phred33 on the leading '@').
        DetectEncoding("@J").Should().Be(QualityEncoding.Phred64);
    }

    [Test]
    public void DetectEncoding_BoundaryI_StaysAmbiguousThenDefaultsPhred33()
    {
        // 'I'(73) is the top of the Sanger range and not > 'I' ⇒ no Phred64 trigger ⇒ default Phred33.
        DetectEncoding("III").Should().Be(QualityEncoding.Phred33);
    }

    // ── DecodeQualityScores: score = max(0, char - offset) ────────────────────────────

    [Test]
    public void DecodeQualityScores_Phred33_SubtractsThirtyThree()
    {
        // '!'=33→0, '+'=43→10, 'I'=73→40.
        DecodeQualityScores("!+I", QualityEncoding.Phred33).Should().Equal(0, 10, 40);
    }

    [Test]
    public void DecodeQualityScores_Phred64_SubtractsSixtyFour()
    {
        // '@'=64→0, 'h'=104→40.
        DecodeQualityScores("@h", QualityEncoding.Phred64).Should().Equal(0, 40);
    }

    // ── ParseHeader: split only when the space is past the first column ─────────────────

    [Test]
    public void Parse_HeaderWithDescription_SplitsIdAndDescription()
    {
        var rec = Parse("@read1 left mate\nACGT\n+\nIIII\n").Single();
        rec.Id.Should().Be("read1");
        rec.Description.Should().Be("left mate");
        rec.Sequence.Should().Be("ACGT");
        rec.QualityString.Should().Be("IIII");
        rec.QualityScores.Should().Equal(40, 40, 40, 40);
    }

    [Test]
    public void Parse_HeaderWithLeadingSpace_KeepsWholeHeaderAsId()
    {
        // spaceIndex == 0 ⇒ guard `spaceIndex > 0` is false ⇒ no split (kills `> 0` → `>= 0`,
        // which would yield an empty Id and shift the description).
        var rec = Parse("@ odd\nAC\n+\nII\n").Single();
        rec.Id.Should().Be(" odd");
        rec.Description.Should().Be("");
    }

    // ── CalculateStatistics: exact EMBOSS-style aggregate values ───────────────────────

    [Test]
    public void CalculateStatistics_TwoReads_ComputesExactAggregates()
    {
        var reads = new[]
        {
            new FastqRecord("r1", "", "ACGT", "", new[] { 20, 30, 40, 10 }), // 4 bp, GC=2
            new FastqRecord("r2", "", "GGC",  "", new[] { 25, 35, 15 })       // 3 bp, GC=3
        };

        var stats = CalculateStatistics(reads);

        stats.TotalReads.Should().Be(2);
        stats.TotalBases.Should().Be(7);
        stats.MinReadLength.Should().Be(3);          // kills `< minLength` boundary
        stats.MaxReadLength.Should().Be(4);          // kills `> maxLength` boundary
        stats.MeanReadLength.Should().Be(3.5);
        stats.MeanQuality.Should().BeApproximately(175.0 / 7.0, 1e-9);
        // Q20 counts scores ≥ 20 (the 20 itself must be included → kills `>= 20` → `> 20`):
        stats.Q20Percentage.Should().BeApproximately(100.0 * 5 / 7, 1e-9);
        // Q30 counts scores ≥ 30 (the 30 itself must be included → kills `>= 30` → `> 30`):
        stats.Q30Percentage.Should().BeApproximately(100.0 * 3 / 7, 1e-9);
        stats.GcContent.Should().BeApproximately(5.0 / 7, 1e-9);
    }

    [Test]
    public void CalculateStatistics_EmptySequenceRead_GuardsDivisionByZero()
    {
        // totalReads = 1 but totalBases = 0 ⇒ the `totalBases > 0 ? … : 0` guards must short-circuit
        // to 0 (kills `> 0` → `>= 0`, which would divide 0/0 → NaN).
        var stats = CalculateStatistics(new[] { new FastqRecord("e", "", "", "", System.Array.Empty<int>()) });

        stats.TotalReads.Should().Be(1);
        stats.TotalBases.Should().Be(0);
        stats.MeanQuality.Should().Be(0);
        stats.Q20Percentage.Should().Be(0);
        stats.Q30Percentage.Should().Be(0);
        stats.GcContent.Should().Be(0);
    }

    [Test]
    public void CalculatePositionQuality_ComputesMeanAndPopulationStdDev()
    {
        var reads = new[]
        {
            new FastqRecord("a", "", "A", "", new[] { 10 }),
            new FastqRecord("b", "", "A", "", new[] { 20 })
        };

        var pos = CalculatePositionQuality(reads);

        pos.Should().HaveCount(1);
        pos[0].Position.Should().Be(1);                 // 1-based reporting
        pos[0].MeanQuality.Should().Be(15);
        // variance = ((10-15)² + (20-15)²)/2 = 25 ⇒ stddev = 5 (kills `(s-mean)*(s-mean)` → `(s+mean)…`).
        pos[0].StdDev.Should().BeApproximately(5.0, 1e-9);
    }

    // ── FilterByLength / FilterByQuality: inclusive bounds, empty-quality exclusion ────

    [Test]
    public void FilterByLength_BoundsAreInclusive()
    {
        var reads = new[]
        {
            new FastqRecord("short", "", "AC", "", new[] { 40, 40 }),    // len 2
            new FastqRecord("min",   "", "ACG", "", new[] { 40, 40, 40 }), // len 3 (== min, kept)
            new FastqRecord("max",   "", "ACGTA", "", new int[5]),        // len 5 (== max, kept)
            new FastqRecord("long",  "", "ACGTAC", "", new int[6])        // len 6 (> max, dropped)
        };

        FilterByLength(reads, 3, 5).Select(r => r.Id)
            .Should().Equal("min", "max");
    }

    [Test]
    public void FilterByQuality_RequiresNonEmptyScores()
    {
        var reads = new[]
        {
            new FastqRecord("good",  "", "AC", "", new[] { 30, 30 }),
            new FastqRecord("empty", "", "",   "", System.Array.Empty<int>()) // Count == 0 ⇒ excluded
        };

        FilterByQuality(reads, 25).Select(r => r.Id).Should().Equal("good");
    }

    // ── TrimByQuality: trim from both ends until score ≥ minQuality ────────────────────

    [Test]
    public void TrimByQuality_TrimsLowQualityEnds_KeepsInteriorInclusiveOfThreshold()
    {
        // scores below 20 at both ends trimmed; the score == 20 is kept (kills `< minQuality`
        // edge and the `start < end` / `end > start` boundaries).
        var rec = new FastqRecord("r", "", "AACGTT", "??????",
            new[] { 5, 19, 20, 40, 10, 3 });

        var trimmed = TrimByQuality(rec, 20);

        // start advances past {5,19} to index 2 (score 20, kept); end retreats past {10,3} to
        // index 4 ⇒ kept window is indices [2,4) = "CG" with scores {20,40}.
        trimmed.Sequence.Should().Be("CG");
        trimmed.QualityScores.Should().Equal(20, 40);
    }

    [Test]
    public void TrimByQuality_AllBelowThreshold_ReturnsEmpty()
    {
        var rec = new FastqRecord("r", "", "ACG", "###", new[] { 2, 3, 4 });
        var trimmed = TrimByQuality(rec, 20);
        trimmed.Sequence.Should().Be("");
        trimmed.QualityScores.Should().BeEmpty();
    }

    // ── TrimAdapter: 3' overlap and internal full-adapter removal ──────────────────────

    [Test]
    public void TrimAdapter_RemovesThreePrimeOverlap()
    {
        // sequence ends with the adapter prefix "ADAP" (overlap 4 ≥ minOverlap) ⇒ trimmed off.
        var rec = new FastqRecord("r", "", "GGGGADAP", "IIIIIIII",
            new[] { 40, 40, 40, 40, 30, 30, 30, 30 });

        var trimmed = TrimAdapter(rec, "ADAPTER", minOverlap: 4);

        trimmed.Sequence.Should().Be("GGGG");
        trimmed.QualityString.Should().Be("IIII");
        trimmed.QualityScores.Should().Equal(40, 40, 40, 40);
    }

    [Test]
    public void TrimAdapter_TooShortToOverlap_ReturnsUnchanged()
    {
        var rec = new FastqRecord("r", "", "ACGT", "IIII", new[] { 40, 40, 40, 40 });
        // adapter shorter than minOverlap ⇒ returned unchanged (kills `adapter.Length < minOverlap`).
        TrimAdapter(rec, "AC", minOverlap: 5).Sequence.Should().Be("ACGT");
    }
}
