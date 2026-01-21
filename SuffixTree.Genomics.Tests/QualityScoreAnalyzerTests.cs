using NUnit.Framework;
using SuffixTree.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class QualityScoreAnalyzerTests
{
    #region Phred Conversion Tests

    [Test]
    public void CharToPhred_Phred33_ConvertsCorrectly()
    {
        // '!' = ASCII 33 = Phred 0
        Assert.That(QualityScoreAnalyzer.CharToPhred('!'), Is.EqualTo(0));
        // 'I' = ASCII 73 = Phred 40
        Assert.That(QualityScoreAnalyzer.CharToPhred('I'), Is.EqualTo(40));
        // '5' = ASCII 53 = Phred 20
        Assert.That(QualityScoreAnalyzer.CharToPhred('5'), Is.EqualTo(20));
    }

    [Test]
    public void CharToPhred_Phred64_ConvertsCorrectly()
    {
        // '@' = ASCII 64 = Phred 0 in Phred+64
        Assert.That(QualityScoreAnalyzer.CharToPhred('@', QualityScoreAnalyzer.QualityEncoding.Phred64), Is.EqualTo(0));
        // 'h' = ASCII 104 = Phred 40 in Phred+64
        Assert.That(QualityScoreAnalyzer.CharToPhred('h', QualityScoreAnalyzer.QualityEncoding.Phred64), Is.EqualTo(40));
    }

    [Test]
    public void PhredToChar_Phred33_ConvertsCorrectly()
    {
        Assert.That(QualityScoreAnalyzer.PhredToChar(0), Is.EqualTo('!'));
        Assert.That(QualityScoreAnalyzer.PhredToChar(40), Is.EqualTo('I'));
        Assert.That(QualityScoreAnalyzer.PhredToChar(20), Is.EqualTo('5'));
    }

    [Test]
    public void PhredToChar_Phred64_ConvertsCorrectly()
    {
        Assert.That(QualityScoreAnalyzer.PhredToChar(0, QualityScoreAnalyzer.QualityEncoding.Phred64), Is.EqualTo('@'));
        Assert.That(QualityScoreAnalyzer.PhredToChar(40, QualityScoreAnalyzer.QualityEncoding.Phred64), Is.EqualTo('h'));
    }

    [Test]
    public void QualityStringToPhred_ConvertsString()
    {
        string quality = "!5I";
        var phred = QualityScoreAnalyzer.QualityStringToPhred(quality);

        Assert.That(phred, Has.Length.EqualTo(3));
        Assert.That(phred[0], Is.EqualTo(0));
        Assert.That(phred[1], Is.EqualTo(20));
        Assert.That(phred[2], Is.EqualTo(40));
    }

    [Test]
    public void PhredToQualityString_ConvertsArray()
    {
        var phred = new[] { 0, 20, 40 };
        string quality = QualityScoreAnalyzer.PhredToQualityString(phred);

        Assert.That(quality, Is.EqualTo("!5I"));
    }

    [Test]
    public void RoundTrip_PreservesQuality()
    {
        string original = "IIIIIIIIIII55555!!!!!";
        var phred = QualityScoreAnalyzer.QualityStringToPhred(original);
        string result = QualityScoreAnalyzer.PhredToQualityString(phred);

        Assert.That(result, Is.EqualTo(original));
    }

    #endregion

    #region Encoding Detection Tests

    [Test]
    public void DetectEncoding_LowAscii_ReturnsPhred33()
    {
        string quality = "!!!!!55555IIIII"; // Contains ASCII 33 ('!')

        var encoding = QualityScoreAnalyzer.DetectEncoding(quality);

        Assert.That(encoding, Is.EqualTo(QualityScoreAnalyzer.QualityEncoding.Phred33));
    }

    [Test]
    public void DetectEncoding_HighAscii_ReturnsPhred64()
    {
        string quality = "hhhhhhhhhhhhhhhhhhhh"; // High ASCII, no low chars

        var encoding = QualityScoreAnalyzer.DetectEncoding(quality);

        Assert.That(encoding, Is.EqualTo(QualityScoreAnalyzer.QualityEncoding.Phred64));
    }

    [Test]
    public void DetectEncoding_Empty_ReturnsPhred33()
    {
        var encoding = QualityScoreAnalyzer.DetectEncoding("");
        Assert.That(encoding, Is.EqualTo(QualityScoreAnalyzer.QualityEncoding.Phred33));
    }

    #endregion

    #region Error Probability Tests

    [Test]
    public void PhredToErrorProbability_CorrectConversion()
    {
        // Q10 = 10% error rate
        Assert.That(QualityScoreAnalyzer.PhredToErrorProbability(10), Is.EqualTo(0.1).Within(0.0001));
        // Q20 = 1% error rate
        Assert.That(QualityScoreAnalyzer.PhredToErrorProbability(20), Is.EqualTo(0.01).Within(0.0001));
        // Q30 = 0.1% error rate
        Assert.That(QualityScoreAnalyzer.PhredToErrorProbability(30), Is.EqualTo(0.001).Within(0.0001));
    }

    [Test]
    public void ErrorProbabilityToPhred_CorrectConversion()
    {
        Assert.That(QualityScoreAnalyzer.ErrorProbabilityToPhred(0.1), Is.EqualTo(10));
        Assert.That(QualityScoreAnalyzer.ErrorProbabilityToPhred(0.01), Is.EqualTo(20));
        Assert.That(QualityScoreAnalyzer.ErrorProbabilityToPhred(0.001), Is.EqualTo(30));
    }

    [Test]
    public void ErrorProbabilityToPhred_EdgeCases()
    {
        Assert.That(QualityScoreAnalyzer.ErrorProbabilityToPhred(0), Is.EqualTo(60));
        Assert.That(QualityScoreAnalyzer.ErrorProbabilityToPhred(1), Is.EqualTo(0));
    }

    #endregion

    #region Quality Statistics Tests

    [Test]
    public void CalculateStatistics_SingleString_CalculatesCorrectly()
    {
        string quality = "IIIII"; // All Q40

        var stats = QualityScoreAnalyzer.CalculateStatistics(quality);

        Assert.That(stats.MeanQuality, Is.EqualTo(40));
        Assert.That(stats.MedianQuality, Is.EqualTo(40));
        Assert.That(stats.MinQuality, Is.EqualTo(40));
        Assert.That(stats.MaxQuality, Is.EqualTo(40));
        Assert.That(stats.TotalBases, Is.EqualTo(5));
        Assert.That(stats.BasesAboveQ30, Is.EqualTo(5));
        Assert.That(stats.PercentAboveQ30, Is.EqualTo(100));
    }

    [Test]
    public void CalculateStatistics_MixedQuality_CalculatesCorrectly()
    {
        string quality = "!5I"; // Q0, Q20, Q40

        var stats = QualityScoreAnalyzer.CalculateStatistics(quality);

        Assert.That(stats.MeanQuality, Is.EqualTo(20));
        Assert.That(stats.MedianQuality, Is.EqualTo(20));
        Assert.That(stats.MinQuality, Is.EqualTo(0));
        Assert.That(stats.MaxQuality, Is.EqualTo(40));
        Assert.That(stats.BasesAboveQ20, Is.EqualTo(2)); // Q20 and Q40
        Assert.That(stats.BasesAboveQ30, Is.EqualTo(1)); // Only Q40
    }

    [Test]
    public void CalculateStatistics_MultipleStrings_CombinesCorrectly()
    {
        var qualities = new[] { "IIIII", "55555" }; // 5x Q40, 5x Q20

        var stats = QualityScoreAnalyzer.CalculateStatistics(qualities);

        Assert.That(stats.MeanQuality, Is.EqualTo(30));
        Assert.That(stats.TotalBases, Is.EqualTo(10));
        Assert.That(stats.PerPositionMeanQuality, Has.Count.EqualTo(5));
    }

    [Test]
    public void CalculateStatistics_EmptyString_ReturnsZeros()
    {
        var stats = QualityScoreAnalyzer.CalculateStatistics("");

        Assert.That(stats.TotalBases, Is.EqualTo(0));
        Assert.That(stats.MeanQuality, Is.EqualTo(0));
    }

    #endregion

    #region Quality Trimming Tests

    [Test]
    public void QualityTrim_LowQualityEnds_Trims()
    {
        string sequence = "ACGTACGT";
        string quality = "!!IIII!!"; // Low at ends

        var result = QualityScoreAnalyzer.QualityTrim(sequence, quality, minQuality: 20);

        Assert.That(result.Sequence, Is.EqualTo("GTAC"));
        Assert.That(result.TrimmedFromStart, Is.EqualTo(2));
        Assert.That(result.TrimmedFromEnd, Is.EqualTo(2));
    }

    [Test]
    public void QualityTrim_AllHighQuality_NoTrimming()
    {
        string sequence = "ACGTACGT";
        string quality = "IIIIIIII";

        var result = QualityScoreAnalyzer.QualityTrim(sequence, quality, minQuality: 20);

        Assert.That(result.Sequence, Is.EqualTo("ACGTACGT"));
        Assert.That(result.TrimmedFromStart, Is.EqualTo(0));
        Assert.That(result.TrimmedFromEnd, Is.EqualTo(0));
    }

    [Test]
    public void QualityTrim_AllLowQuality_ReturnsEmpty()
    {
        string sequence = "ACGTACGT";
        string quality = "!!!!!!!!";

        var result = QualityScoreAnalyzer.QualityTrim(sequence, quality, minQuality: 20);

        Assert.That(result.Sequence, Is.EqualTo(""));
        Assert.That(result.FinalLength, Is.EqualTo(0));
    }

    [Test]
    public void SlidingWindowTrim_LowQualityEnd_Trims()
    {
        string sequence = "ACGTACGTACGT";
        string quality = "IIIIIIII!!!!"; // Low quality at end

        var result = QualityScoreAnalyzer.SlidingWindowTrim(
            sequence, quality, windowSize: 4, minAverageQuality: 20);

        Assert.That(result.FinalLength, Is.LessThan(result.OriginalLength));
        Assert.That(result.TrimmedFromEnd, Is.GreaterThan(0));
    }

    [Test]
    public void SlidingWindowTrim_AllHighQuality_NoTrimming()
    {
        string sequence = "ACGTACGT";
        string quality = "IIIIIIII";

        var result = QualityScoreAnalyzer.SlidingWindowTrim(
            sequence, quality, windowSize: 4, minAverageQuality: 20);

        Assert.That(result.Sequence, Is.EqualTo(sequence));
    }

    #endregion

    #region Read Filtering Tests

    [Test]
    public void FilterReads_ByLength_FiltersCorrectly()
    {
        var reads = new List<QualityScoreAnalyzer.FastqRecord>
        {
            new("read1", "ACGT", "IIII"),
            new("read2", "ACGTACGT", "IIIIIIII"),
            new("read3", "AC", "II")
        };

        var filtered = QualityScoreAnalyzer.FilterReads(reads, minLength: 4).ToList();

        Assert.That(filtered, Has.Count.EqualTo(2));
        Assert.That(filtered.All(r => r.Sequence.Length >= 4), Is.True);
    }

    [Test]
    public void FilterReads_ByMeanQuality_FiltersCorrectly()
    {
        var reads = new List<QualityScoreAnalyzer.FastqRecord>
        {
            new("good", "ACGT", "IIII"),  // Q40
            new("bad", "ACGT", "!!!!"),   // Q0
            new("medium", "ACGT", "5555") // Q20
        };

        var filtered = QualityScoreAnalyzer.FilterReads(reads, minMeanQuality: 25).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
        Assert.That(filtered[0].Id, Is.EqualTo("good"));
    }

    [Test]
    public void FilterReads_ByExpectedErrors_FiltersCorrectly()
    {
        var reads = new List<QualityScoreAnalyzer.FastqRecord>
        {
            new("good", "ACGT", "IIII"),  // Very low expected errors
            new("bad", "ACGT", "!!!!"),   // ~4 expected errors
        };

        var filtered = QualityScoreAnalyzer.FilterReads(reads, maxExpectedErrors: 1.0).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
        Assert.That(filtered[0].Id, Is.EqualTo("good"));
    }

    #endregion

    #region Expected Errors Tests

    [Test]
    public void CalculateExpectedErrors_HighQuality_LowErrors()
    {
        string quality = "IIIIIIIIII"; // Q40 x 10

        double errors = QualityScoreAnalyzer.CalculateExpectedErrors(quality);

        Assert.That(errors, Is.LessThan(0.01));
    }

    [Test]
    public void CalculateExpectedErrors_LowQuality_HighErrors()
    {
        string quality = "!!!!!!!!!!"; // Q0 x 10

        double errors = QualityScoreAnalyzer.CalculateExpectedErrors(quality);

        Assert.That(errors, Is.EqualTo(10).Within(0.1));
    }

    [Test]
    public void CalculateExpectedErrors_Q20_OnePercentEach()
    {
        string quality = "5555555555"; // Q20 x 10, each base 1% error

        double errors = QualityScoreAnalyzer.CalculateExpectedErrors(quality);

        Assert.That(errors, Is.EqualTo(0.1).Within(0.01));
    }

    #endregion

    #region Base Masking Tests

    [Test]
    public void MaskLowQualityBases_MasksLowQuality()
    {
        string sequence = "ACGTACGT";
        string quality = "II!!II!!";

        string masked = QualityScoreAnalyzer.MaskLowQualityBases(sequence, quality, minQuality: 20);

        Assert.That(masked, Is.EqualTo("ACNNACNN"));
    }

    [Test]
    public void MaskLowQualityBases_AllHighQuality_NoMasking()
    {
        string sequence = "ACGTACGT";
        string quality = "IIIIIIII";

        string masked = QualityScoreAnalyzer.MaskLowQualityBases(sequence, quality, minQuality: 20);

        Assert.That(masked, Is.EqualTo("ACGTACGT"));
    }

    [Test]
    public void MaskLowQualityBases_AllLowQuality_AllMasked()
    {
        string sequence = "ACGTACGT";
        string quality = "!!!!!!!!";

        string masked = QualityScoreAnalyzer.MaskLowQualityBases(sequence, quality, minQuality: 20);

        Assert.That(masked, Is.EqualTo("NNNNNNNN"));
    }

    #endregion

    #region FASTQ Parsing Tests

    [Test]
    public void ParseFastq_ValidRecord_ParsesCorrectly()
    {
        var lines = new[]
        {
            "@read1 description",
            "ACGTACGT",
            "+",
            "IIIIIIII"
        };

        var records = QualityScoreAnalyzer.ParseFastq(lines).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].Id, Is.EqualTo("read1"));
        Assert.That(records[0].Sequence, Is.EqualTo("ACGTACGT"));
        Assert.That(records[0].QualityString, Is.EqualTo("IIIIIIII"));
        Assert.That(records[0].Description, Is.EqualTo("description"));
    }

    [Test]
    public void ParseFastq_MultipleRecords_ParsesAll()
    {
        var lines = new[]
        {
            "@read1", "ACGT", "+", "IIII",
            "@read2", "TGCA", "+", "5555"
        };

        var records = QualityScoreAnalyzer.ParseFastq(lines).ToList();

        Assert.That(records, Has.Count.EqualTo(2));
    }

    [Test]
    public void ParseFastq_NoDescription_DescriptionIsNull()
    {
        var lines = new[] { "@read1", "ACGT", "+", "IIII" };

        var records = QualityScoreAnalyzer.ParseFastq(lines).ToList();

        Assert.That(records[0].Description, Is.Null);
    }

    [Test]
    public void ToFastq_FormatsCorrectly()
    {
        var record = new QualityScoreAnalyzer.FastqRecord("read1", "ACGT", "IIII", "desc");

        var lines = QualityScoreAnalyzer.ToFastq(record).ToList();

        Assert.That(lines, Has.Count.EqualTo(4));
        Assert.That(lines[0], Is.EqualTo("@read1 desc"));
        Assert.That(lines[1], Is.EqualTo("ACGT"));
        Assert.That(lines[2], Is.EqualTo("+"));
        Assert.That(lines[3], Is.EqualTo("IIII"));
    }

    [Test]
    public void ToFastq_NoDescription_OmitsDescription()
    {
        var record = new QualityScoreAnalyzer.FastqRecord("read1", "ACGT", "IIII");

        var lines = QualityScoreAnalyzer.ToFastq(record).ToList();

        Assert.That(lines[0], Is.EqualTo("@read1"));
    }

    [Test]
    public void RoundTrip_FastqParsing_PreservesData()
    {
        var original = new QualityScoreAnalyzer.FastqRecord("read1", "ACGTACGT", "IIIIIIII", "test");

        var lines = QualityScoreAnalyzer.ToFastq(original).ToList();
        var parsed = QualityScoreAnalyzer.ParseFastq(lines).First();

        Assert.That(parsed.Id, Is.EqualTo(original.Id));
        Assert.That(parsed.Sequence, Is.EqualTo(original.Sequence));
        Assert.That(parsed.QualityString, Is.EqualTo(original.QualityString));
        Assert.That(parsed.Description, Is.EqualTo(original.Description));
    }

    #endregion

    #region Quality Distribution Tests

    [Test]
    public void GetQualityDistribution_CountsCorrectly()
    {
        string quality = "IIIII55!!!"; // 5x Q40, 2x Q20, 3x Q0

        var dist = QualityScoreAnalyzer.GetQualityDistribution(quality);

        Assert.That(dist[40], Is.EqualTo(5));
        Assert.That(dist[20], Is.EqualTo(2));
        Assert.That(dist[0], Is.EqualTo(3));
    }

    #endregion

    #region Low Quality Regions Tests

    [Test]
    public void FindLowQualityRegions_IdentifiesLowRegions()
    {
        string quality = "IIIIIIII!!!!!!!!!!!!IIIIIIII"; // Low in middle

        var regions = QualityScoreAnalyzer.FindLowQualityRegions(
            quality, windowSize: 5, maxQuality: 10).ToList();

        Assert.That(regions.Count, Is.GreaterThan(0));
    }

    [Test]
    public void FindLowQualityRegions_AllHighQuality_NoRegions()
    {
        string quality = "IIIIIIIIIIIIIIIIIIII";

        var regions = QualityScoreAnalyzer.FindLowQualityRegions(
            quality, windowSize: 5, maxQuality: 10).ToList();

        Assert.That(regions.Count, Is.EqualTo(0));
    }

    #endregion

    #region Consensus Quality Tests

    [Test]
    public void CalculateConsensusQuality_TakesMaximum()
    {
        var qualities = new[] { "!5I", "I5!", "555" };

        string consensus = QualityScoreAnalyzer.CalculateConsensusQuality(qualities);

        var phred = QualityScoreAnalyzer.QualityStringToPhred(consensus);
        Assert.That(phred[0], Is.EqualTo(40)); // Max of 0, 40, 20
        Assert.That(phred[1], Is.EqualTo(20)); // All 20
        Assert.That(phred[2], Is.EqualTo(40)); // Max of 40, 0, 20
    }

    [Test]
    public void CalculateConsensusQuality_Empty_ReturnsEmpty()
    {
        string consensus = QualityScoreAnalyzer.CalculateConsensusQuality(Array.Empty<string>());
        Assert.That(consensus, Is.EqualTo(""));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void QualityStringToPhred_Empty_ReturnsEmpty()
    {
        var phred = QualityScoreAnalyzer.QualityStringToPhred("");
        Assert.That(phred, Is.Empty);
    }

    [Test]
    public void QualityTrim_Empty_ReturnsEmpty()
    {
        var result = QualityScoreAnalyzer.QualityTrim("", "");
        Assert.That(result.Sequence, Is.EqualTo(""));
    }

    [Test]
    public void MaskLowQualityBases_Empty_ReturnsEmpty()
    {
        string masked = QualityScoreAnalyzer.MaskLowQualityBases("", "");
        Assert.That(masked, Is.EqualTo(""));
    }

    #endregion
}
