using NUnit.Framework;
using Seqeron.Genomics.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the FASTQ format parser.
/// Evidence: Wikipedia (FASTQ format), Cock et al. (2009), NCBI SRA File Format Guide
/// Test Unit: PARSE-FASTQ-001
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

    private const string HighQualityFastq = @"@high_quality
ACGTACGT
+
IIIIIIII";

    private const string LowQualityFastq = @"@low_quality
ACGTACGT
+
!!!!!!!!";

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

    [Test]
    public void Parse_HeaderWithSpace_SeparatesIdAndDescription()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();

        Assert.That(records[0].Id, Is.EqualTo("SEQ_ID_1"));
        Assert.That(records[0].Description, Is.EqualTo("description"));
    }

    #endregion

    #region Quality Encoding Tests

    [Test]
    public void DetectEncoding_Phred33_ReturnsPhred33()
    {
        // Quality string with chars < '@' indicates Phred33 (Evidence: Wikipedia)
        var encoding = FastqParser.DetectEncoding("!!!!!IIIII");
        Assert.That(encoding, Is.EqualTo(FastqParser.QualityEncoding.Phred33));
    }

    [Test]
    public void DetectEncoding_Phred64_ReturnsPhred64()
    {
        // Quality string with chars > 'I' indicates Phred64 (Evidence: Wikipedia)
        var encoding = FastqParser.DetectEncoding("hhhhhhhh");
        Assert.That(encoding, Is.EqualTo(FastqParser.QualityEncoding.Phred64));
    }

    [Test]
    public void DetectEncoding_AmbiguousRange_DefaultsToPhred33()
    {
        // '@' to 'I' range is ambiguous; should default to Phred33 (modern standard)
        var encoding = FastqParser.DetectEncoding("ABCDEFGHI");
        Assert.That(encoding, Is.EqualTo(FastqParser.QualityEncoding.Phred33));
    }

    [Test]
    public void DetectEncoding_EmptyString_ReturnsPhred33()
    {
        var encoding = FastqParser.DetectEncoding("");
        Assert.That(encoding, Is.EqualTo(FastqParser.QualityEncoding.Phred33));
    }

    [Test]
    public void DecodeQualityScores_Phred33_ReturnsCorrectScores()
    {
        // '!' = ASCII 33, Phred33 score = 0
        // 'I' = ASCII 73, Phred33 score = 40
        // Evidence: Wikipedia FASTQ encoding table
        var scores = FastqParser.DecodeQualityScores("!I", FastqParser.QualityEncoding.Phred33);

        Assert.That(scores[0], Is.EqualTo(0));
        Assert.That(scores[1], Is.EqualTo(40));
    }

    [Test]
    public void DecodeQualityScores_Phred64_ReturnsCorrectScores()
    {
        // '@' = ASCII 64, Phred64 score = 0
        // 'h' = ASCII 104, Phred64 score = 40
        // Evidence: Wikipedia FASTQ encoding table
        var scores = FastqParser.DecodeQualityScores("@h", FastqParser.QualityEncoding.Phred64);

        Assert.That(scores[0], Is.EqualTo(0));
        Assert.That(scores[1], Is.EqualTo(40));
    }

    [Test]
    public void DecodeQualityScores_EmptyString_ReturnsEmptyArray()
    {
        var scores = FastqParser.DecodeQualityScores("", FastqParser.QualityEncoding.Phred33);
        Assert.That(scores, Is.Empty);
    }

    [Test]
    public void DecodeQualityScores_NullString_ReturnsEmptyArray()
    {
        var scores = FastqParser.DecodeQualityScores(null!, FastqParser.QualityEncoding.Phred33);
        Assert.That(scores, Is.Empty);
    }

    #endregion

    #region Phred Mathematics Tests

    [Test]
    public void PhredToErrorProbability_Q10_Returns0Point1()
    {
        // Q = -10 × log₁₀(p), so Q10 → p = 0.1 (Evidence: Wikipedia Phred formula)
        var probability = FastqParser.PhredToErrorProbability(10);
        Assert.That(probability, Is.EqualTo(0.1).Within(0.0001));
    }

    [Test]
    public void PhredToErrorProbability_Q20_Returns0Point01()
    {
        // Q20 → p = 0.01 (1% error rate)
        var probability = FastqParser.PhredToErrorProbability(20);
        Assert.That(probability, Is.EqualTo(0.01).Within(0.0001));
    }

    [Test]
    public void PhredToErrorProbability_Q30_Returns0Point001()
    {
        // Q30 → p = 0.001 (0.1% error rate)
        var probability = FastqParser.PhredToErrorProbability(30);
        Assert.That(probability, Is.EqualTo(0.001).Within(0.0001));
    }

    [Test]
    public void PhredToErrorProbability_Q40_Returns0Point0001()
    {
        // Q40 → p = 0.0001 (0.01% error rate)
        var probability = FastqParser.PhredToErrorProbability(40);
        Assert.That(probability, Is.EqualTo(0.0001).Within(0.00001));
    }

    [Test]
    public void ErrorProbabilityToPhred_0Point1_ReturnsQ10()
    {
        // p = 0.1 → Q = 10 (Evidence: Wikipedia inverse formula)
        var phred = FastqParser.ErrorProbabilityToPhred(0.1);
        Assert.That(phred, Is.EqualTo(10));
    }

    [Test]
    public void ErrorProbabilityToPhred_0Point01_ReturnsQ20()
    {
        var phred = FastqParser.ErrorProbabilityToPhred(0.01);
        Assert.That(phred, Is.EqualTo(20));
    }

    [Test]
    public void ErrorProbabilityToPhred_0Point001_ReturnsQ30()
    {
        var phred = FastqParser.ErrorProbabilityToPhred(0.001);
        Assert.That(phred, Is.EqualTo(30));
    }

    [Test]
    public void ErrorProbabilityToPhred_ZeroOrNegative_ReturnsMaxQuality()
    {
        // Zero probability → max quality (Q40 typical max)
        var phred = FastqParser.ErrorProbabilityToPhred(0);
        Assert.That(phred, Is.EqualTo(40));
    }

    #endregion

    #region Quality Encoding Round-Trip Tests

    [Test]
    public void EncodeQualityScores_Phred33_EncodesCorrectly()
    {
        // Q0 → '!', Q40 → 'I'
        var encoded = FastqParser.EncodeQualityScores(new[] { 0, 40 }, FastqParser.QualityEncoding.Phred33);
        Assert.That(encoded, Is.EqualTo("!I"));
    }

    [Test]
    public void EncodeQualityScores_Phred64_EncodesCorrectly()
    {
        // Q0 → '@', Q40 → 'h'
        var encoded = FastqParser.EncodeQualityScores(new[] { 0, 40 }, FastqParser.QualityEncoding.Phred64);
        Assert.That(encoded, Is.EqualTo("@h"));
    }

    [Test]
    public void EncodeDecodeRoundTrip_Phred33_PreservesScores()
    {
        var originalScores = new[] { 0, 10, 20, 30, 40 };
        var encoded = FastqParser.EncodeQualityScores(originalScores, FastqParser.QualityEncoding.Phred33);
        var decoded = FastqParser.DecodeQualityScores(encoded, FastqParser.QualityEncoding.Phred33);

        Assert.That(decoded, Is.EqualTo(originalScores));
    }

    [Test]
    public void EncodeDecodeRoundTrip_Phred64_PreservesScores()
    {
        var originalScores = new[] { 0, 10, 20, 30, 40 };
        var encoded = FastqParser.EncodeQualityScores(originalScores, FastqParser.QualityEncoding.Phred64);
        var decoded = FastqParser.DecodeQualityScores(encoded, FastqParser.QualityEncoding.Phred64);

        Assert.That(decoded, Is.EqualTo(originalScores));
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

    [Test]
    public void WriteToFile_CreatesValidFastq()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            var records = FastqParser.Parse(SimpleFastq).ToList();
            FastqParser.WriteToFile(tempFile, records);

            Assert.That(File.Exists(tempFile), Is.True);
            var content = File.ReadAllText(tempFile);
            Assert.That(content, Does.Contain("@SEQ_ID_1"));
            Assert.That(content, Does.Contain("GATCGATCGATCGATC"));
        }
        finally
        {
            File.Delete(tempFile);
        }
    }

    [Test]
    public void WriteAndParseRoundTrip_PreservesRecords()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            var originalRecords = FastqParser.Parse(SimpleFastq).ToList();
            FastqParser.WriteToFile(tempFile, originalRecords);
            var parsedRecords = FastqParser.ParseFile(tempFile).ToList();

            Assert.That(parsedRecords, Has.Count.EqualTo(originalRecords.Count));
            for (int i = 0; i < originalRecords.Count; i++)
            {
                Assert.That(parsedRecords[i].Id, Is.EqualTo(originalRecords[i].Id));
                Assert.That(parsedRecords[i].Sequence, Is.EqualTo(originalRecords[i].Sequence));
                Assert.That(parsedRecords[i].QualityString, Is.EqualTo(originalRecords[i].QualityString));
            }
        }
        finally
        {
            File.Delete(tempFile);
        }
    }

    [Test]
    public void ToFastqString_FormatsCorrectly()
    {
        var record = new FastqParser.FastqRecord("test_id", "description", "ACGT", "IIII", new[] { 40, 40, 40, 40 });
        var fastqString = FastqParser.ToFastqString(record);

        Assert.That(fastqString, Does.StartWith("@test_id description"));
        Assert.That(fastqString, Does.Contain("ACGT"));
        Assert.That(fastqString, Does.Contain("+"));
        Assert.That(fastqString, Does.Contain("IIII"));
    }

    #endregion

    #region Additional Filtering Tests

    [Test]
    public void FilterByQuality_KeepsRecordsAtThreshold()
    {
        // Q40 records should pass threshold of 40
        var records = FastqParser.Parse(HighQualityFastq).ToList();
        var filtered = FastqParser.FilterByQuality(records, 40).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
    }

    [Test]
    public void FilterByQuality_RemovesRecordsBelowThreshold()
    {
        // Q0 records should not pass threshold of 30
        var records = FastqParser.Parse(LowQualityFastq).ToList();
        var filtered = FastqParser.FilterByQuality(records, 30).ToList();

        Assert.That(filtered, Is.Empty);
    }

    #endregion

    #region Additional Trimming Tests

    [Test]
    public void TrimByQuality_AllHighQuality_ReturnsUnchanged()
    {
        var records = FastqParser.Parse(HighQualityFastq).ToList();
        var trimmed = FastqParser.TrimByQuality(records[0], minQuality: 30);

        Assert.That(trimmed.Sequence, Is.EqualTo("ACGTACGT"));
        Assert.That(trimmed.Sequence.Length, Is.EqualTo(8));
    }

    [Test]
    public void TrimByQuality_AllLowQuality_ReturnsEmptySequence()
    {
        var records = FastqParser.Parse(LowQualityFastq).ToList();
        var trimmed = FastqParser.TrimByQuality(records[0], minQuality: 30);

        Assert.That(trimmed.Sequence, Is.Empty);
    }

    #endregion

    #region Additional Statistics Tests

    [Test]
    public void CalculateStatistics_EmptyInput_ReturnsZeros()
    {
        var stats = FastqParser.CalculateStatistics(Array.Empty<FastqParser.FastqRecord>());

        Assert.Multiple(() =>
        {
            Assert.That(stats.TotalReads, Is.EqualTo(0));
            Assert.That(stats.TotalBases, Is.EqualTo(0));
            Assert.That(stats.MeanReadLength, Is.EqualTo(0));
            Assert.That(stats.MeanQuality, Is.EqualTo(0));
        });
    }

    [Test]
    public void CalculateStatistics_Q20Percentage_InValidRange()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.Q20Percentage, Is.InRange(0, 100));
    }

    [Test]
    public void CalculateStatistics_Q30Percentage_InValidRange()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.Q30Percentage, Is.InRange(0, 100));
    }

    [Test]
    public void CalculateStatistics_GcContent_InValidRange()
    {
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.GcContent, Is.InRange(0, 1));
    }

    [Test]
    public void CalculateStatistics_HighQualityReads_HasHighQ30()
    {
        var records = FastqParser.Parse(HighQualityFastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        // All Q40 should mean 100% Q30
        Assert.That(stats.Q30Percentage, Is.EqualTo(100));
    }

    #endregion

    #region Additional Paired-End Tests

    [Test]
    public void InterleavePairedReads_UnequalLengths_StopsAtShorter()
    {
        const string r1 = @"@read1/1
ACGT
+
IIII
@read2/1
TGCA
+
HHHH";

        const string r2 = @"@read1/2
GGGG
+
IIII";

        var reads1 = FastqParser.Parse(r1).ToList();
        var reads2 = FastqParser.Parse(r2).ToList();

        var interleaved = FastqParser.InterleavePairedReads(reads1, reads2).ToList();

        // Should only interleave up to the shorter list
        Assert.That(interleaved, Has.Count.EqualTo(2)); // 1 pair = 2 records
    }

    #endregion
}
