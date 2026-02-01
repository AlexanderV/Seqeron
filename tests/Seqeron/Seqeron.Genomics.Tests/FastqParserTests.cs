using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the FASTQ format parser.
/// </summary>
[TestFixture]
public class FastqParserTests
{
    #region Test Data

    private const string SimpleFastq = @"@SEQ_ID_1 description
GATCGATCGATCGATC
+
IIIIIIIIIIIIIIII
@SEQ_ID_2
ACGTACGTACGTACGT
+
HHHHHHHHHHHHHHHH";

    private const string FastqWithVariousQuality = @"@read1
ACGTACGTACGT
+
!!!!!!!!!!!!
@read2
ACGTACGTACGT
+
IIIIIIIIIIII
@read3
ACGTACGTACGT
+
~~~~~~~~~~~~";

    private const string Phred64Fastq = @"@read1
ACGT
+
hhhh";

    #endregion

    #region Basic Parsing Tests

    [Test]
    public void Parse_SimpleFastq_ReturnsCorrectRecords()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();

        Assert.That(records, Has.Count.EqualTo(2));
        Assert.That(records[0].Id, Is.EqualTo("SEQ_ID_1"));
        Assert.That(records[0].Description, Is.EqualTo("description"));
        Assert.That(records[0].Sequence, Is.EqualTo("GATCGATCGATCGATC"));
        Assert.That(records[0].QualityString, Is.EqualTo("IIIIIIIIIIIIIIII"));
    }

    [Test]
    public void Parse_EmptyContent_ReturnsEmpty()
    {
        var records = FastqParser.Parse("").ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_NullContent_ReturnsEmpty()
    {
        var records = FastqParser.Parse((string)null!).ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void Parse_WithNoDescription_ParsesCorrectly()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();

        Assert.That(records[1].Id, Is.EqualTo("SEQ_ID_2"));
        Assert.That(records[1].Description, Is.Null.Or.Empty);
    }

    [Test]
    public void Parse_RecordSequenceLength_MatchesQualityLength()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();

        foreach (var record in records)
        {
            Assert.That(record.Sequence.Length, Is.EqualTo(record.QualityString.Length),
                $"Sequence and quality length mismatch for {record.Id}");
        }
    }

    #endregion

    #region Quality Encoding Tests

    [Test]
    public void DetectEncoding_Phred33_ReturnsPhred33()
    {
        // Quality string with chars < '@' indicates Phred33
        var encoding = FastqParser.DetectEncoding("!!!!!IIIII");
        Assert.That(encoding, Is.EqualTo(FastqParser.QualityEncoding.Phred33));
    }

    [Test]
    public void DetectEncoding_Phred64_ReturnsPhred64()
    {
        // Quality string with chars > 'I' indicates Phred64
        var encoding = FastqParser.DetectEncoding("hhhhhhhh");
        Assert.That(encoding, Is.EqualTo(FastqParser.QualityEncoding.Phred64));
    }

    [Test]
    public void DecodeQualityScores_Phred33_ReturnsCorrectScores()
    {
        // '!' = 33, so Phred33 score = 0
        // 'I' = 73, so Phred33 score = 40
        var scores = FastqParser.DecodeQualityScores("!I", FastqParser.QualityEncoding.Phred33);

        Assert.That(scores[0], Is.EqualTo(0));
        Assert.That(scores[1], Is.EqualTo(40));
    }

    [Test]
    public void DecodeQualityScores_Phred64_ReturnsCorrectScores()
    {
        // '@' = 64, so Phred64 score = 0
        // 'h' = 104, so Phred64 score = 40
        var scores = FastqParser.DecodeQualityScores("@h", FastqParser.QualityEncoding.Phred64);

        Assert.That(scores[0], Is.EqualTo(0));
        Assert.That(scores[1], Is.EqualTo(40));
    }

    #endregion

    #region Filtering Tests

    [Test]
    public void FilterByQuality_FiltersLowQuality()
    {
        var records = FastqParser.Parse(FastqWithVariousQuality).ToList();
        var filtered = FastqParser.FilterByQuality(records, 30).ToList();

        // Only high quality reads should pass
        Assert.That(filtered.Count, Is.LessThan(records.Count));
    }

    [Test]
    public void FilterByLength_FiltersShortReads()
    {
        const string fastq = @"@short
ACGT
+
IIII
@long
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();
        var filtered = FastqParser.FilterByLength(records, minLength: 10).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
        Assert.That(filtered[0].Id, Is.EqualTo("long"));
    }

    [Test]
    public void FilterByLength_WithMaxLength_FiltersBoth()
    {
        const string fastq = @"@short
ACGT
+
IIII
@medium
ACGTACGT
+
IIIIIIII
@long
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();
        var filtered = FastqParser.FilterByLength(records, minLength: 5, maxLength: 10).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
        Assert.That(filtered[0].Id, Is.EqualTo("medium"));
    }

    #endregion

    #region Trimming Tests

    [Test]
    public void TrimByQuality_TrimsLowQualityEnds()
    {
        const string fastq = @"@read1
ACGTACGTACGT
+
!!IIIIIIII!!";

        var records = FastqParser.Parse(fastq).ToList();
        var trimmed = records.Select(r => FastqParser.TrimByQuality(r, minQuality: 30)).ToList();

        Assert.That(trimmed, Has.Count.EqualTo(1));
        Assert.That(trimmed[0].Sequence.Length, Is.LessThan(12));
    }

    [Test]
    public void TrimAdapter_RemovesAdapter()
    {
        const string adapter = "AGATCGGAAGAG";
        const string fastq = @"@read1
ACGTACGTACGTAAAAGATCGGAAGAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();
        var trimmed = records.Select(r => FastqParser.TrimAdapter(r, adapter)).ToList();

        Assert.That(trimmed, Has.Count.EqualTo(1));
        Assert.That(trimmed[0].Sequence, Does.Not.Contain(adapter));
    }

    [Test]
    public void TrimAdapter_NoAdapter_ReturnsUnchanged()
    {
        const string adapter = "AGATCGGAAGAG";
        const string fastq = @"@read1
ACGTACGTACGT
+
IIIIIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();
        var trimmed = records.Select(r => FastqParser.TrimAdapter(r, adapter)).ToList();

        Assert.That(trimmed[0].Sequence, Is.EqualTo("ACGTACGTACGT"));
    }

    #endregion

    #region Statistics Tests

    [Test]
    public void CalculateStatistics_ReturnsCorrectStats()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.TotalReads, Is.EqualTo(2));
        Assert.That(stats.TotalBases, Is.EqualTo(32)); // 16 + 16
        Assert.That(stats.MeanReadLength, Is.EqualTo(16));
        Assert.That(stats.MinReadLength, Is.EqualTo(16));
        Assert.That(stats.MaxReadLength, Is.EqualTo(16));
    }

    [Test]
    public void CalculateStatistics_VariousLengths_CorrectMinMax()
    {
        const string fastq = @"@short
ACGT
+
IIII
@long
ACGTACGTACGTACGT
+
IIIIIIIIIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.MinReadLength, Is.EqualTo(4));
        Assert.That(stats.MaxReadLength, Is.EqualTo(16));
    }

    [Test]
    public void CalculatePositionQuality_ReturnsQualityPerPosition()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var positionQuality = FastqParser.CalculatePositionQuality(records);

        Assert.That(positionQuality.Count, Is.EqualTo(16));
        Assert.That(positionQuality.All(q => q.MeanQuality > 0), Is.True);
    }

    #endregion

    #region Paired-End Tests

    [Test]
    public void InterleavePairedReads_CombinesReads()
    {
        const string r1 = @"@read1/1
ACGTACGT
+
IIIIIIII";

        const string r2 = @"@read1/2
TGCATGCA
+
HHHHHHHH";

        var reads1 = FastqParser.Parse(r1).ToList();
        var reads2 = FastqParser.Parse(r2).ToList();

        var interleaved = FastqParser.InterleavePairedReads(reads1, reads2).ToList();

        Assert.That(interleaved, Has.Count.EqualTo(2));
        Assert.That(interleaved[0].Sequence, Is.EqualTo("ACGTACGT"));
        Assert.That(interleaved[1].Sequence, Is.EqualTo("TGCATGCA"));
    }

    [Test]
    public void SplitInterleavedReads_SeparatesReads()
    {
        const string interleaved = @"@read1/1
ACGTACGT
+
IIIIIIII
@read1/2
TGCATGCA
+
HHHHHHHH
@read2/1
AAAAAAAA
+
IIIIIIII
@read2/2
TTTTTTTT
+
HHHHHHHH";

        var records = FastqParser.Parse(interleaved).ToList();
        var (r1, r2) = FastqParser.SplitInterleavedReads(records);

        var reads1 = r1.ToList();
        var reads2 = r2.ToList();

        Assert.That(reads1, Has.Count.EqualTo(2));
        Assert.That(reads2, Has.Count.EqualTo(2));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void Parse_MultiplePlusLines_ParsesCorrectly()
    {
        const string fastq = @"@read1
ACGT+ACGT
+
IIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        // Sequence may contain + character
    }

    [Test]
    public void Parse_EmptyRecords_Skipped()
    {
        const string fastq = @"@read1
ACGT
+
IIII

@read2
TGCA
+
HHHH";

        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records, Has.Count.EqualTo(2));
    }

    [Test]
    public void FilterByQuality_EmptyInput_ReturnsEmpty()
    {
        var filtered = FastqParser.FilterByQuality(Array.Empty<FastqParser.FastqRecord>(), 30).ToList();
        Assert.That(filtered, Is.Empty);
    }

    #endregion

    #region File I/O Tests

    [Test]
    public void ParseFile_NonexistentFile_ReturnsEmpty()
    {
        var records = FastqParser.ParseFile("nonexistent.fastq").ToList();
        Assert.That(records, Is.Empty);
    }

    [Test]
    public void ParseFile_ValidFile_ParsesRecords()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            File.WriteAllText(tempFile, SimpleFastq);
            var records = FastqParser.ParseFile(tempFile).ToList();

            Assert.That(records, Has.Count.EqualTo(2));
        }
        finally
        {
            File.Delete(tempFile);
        }
    }

    #endregion
}
