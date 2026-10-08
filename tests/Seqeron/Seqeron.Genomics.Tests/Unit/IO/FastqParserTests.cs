namespace Seqeron.Genomics.Tests.Unit.IO;

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
        // '!' = ASCII 33 → Q0, 'I' = ASCII 73 → Q40, '~' = ASCII 126 → Q93
        // Evidence: Wikipedia FASTQ encoding table, Sanger format ASCII 33-126
        var scores = FastqParser.DecodeQualityScores("!I~", FastqParser.QualityEncoding.Phred33);

        Assert.That(scores[0], Is.EqualTo(0));
        Assert.That(scores[1], Is.EqualTo(40));
        Assert.That(scores[2], Is.EqualTo(93));
    }

    [Test]
    public void DecodeQualityScores_Phred64_ReturnsCorrectScores()
    {
        // '@' = ASCII 64 → Q0, 'h' = ASCII 104 → Q40, '~' = ASCII 126 → Q62
        // Evidence: Wikipedia FASTQ encoding table, Phred+64 ASCII 64-126
        var scores = FastqParser.DecodeQualityScores("@h~", FastqParser.QualityEncoding.Phred64);

        Assert.That(scores[0], Is.EqualTo(0));
        Assert.That(scores[1], Is.EqualTo(40));
        Assert.That(scores[2], Is.EqualTo(62));
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
    public void PhredToErrorProbability_Q0_Returns1()
    {
        // Q0 → p = 1.0 (100% error, lowest quality)
        // Evidence: Wikipedia Phred quality score Symbols table: '!' Q0 → P = 1.000
        var probability = FastqParser.PhredToErrorProbability(0);
        Assert.That(probability, Is.EqualTo(1.0).Within(0.0001));
    }

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
    public void ErrorProbabilityToPhred_Zero_ReturnsMaxQuality()
    {
        // Zero probability → max representable quality (Q93 per Sanger/Phred+33 range)
        // Evidence: Cock et al. 2010 — Sanger encodes Q 0-93 (ASCII 33-126)
        var phred = FastqParser.ErrorProbabilityToPhred(0);
        Assert.That(phred, Is.EqualTo(93));
    }

    [Test]
    public void ErrorProbabilityToPhred_TinyProbability_CappedAtQ93_Monotone()
    {
        // Review 2026-09 F8: -10·log10(1e-10) = 100 exceeds the Phred+33 maximum Q93 (Cock et al. 2010;
        // Biopython truncates written Sanger qualities at 93). Previously 1e-10 → 100 while 0 → 93,
        // i.e. a SMALLER error probability gave a LOWER score. Now capped: monotone non-increasing.
        Assert.That(FastqParser.ErrorProbabilityToPhred(1e-10), Is.EqualTo(93));
        Assert.That(FastqParser.ErrorProbabilityToPhred(1e-9), Is.EqualTo(90));
        Assert.That(FastqParser.ErrorProbabilityToPhred(1.0), Is.EqualTo(0));
        Assert.That(FastqParser.ErrorProbabilityToPhred(0.5), Is.EqualTo(3)); // round(3.0103)
    }

    [TestCase(-0.1)]
    [TestCase(1.5)]
    [TestCase(double.NaN)]
    public void ErrorProbabilityToPhred_OutsideProbabilityDomain_Throws(double p)
    {
        // A probability lies in [0, 1]; previously 1.5 → -2 (negative Phred) and NaN → garbage.
        Assert.Throws<ArgumentOutOfRangeException>(() => FastqParser.ErrorProbabilityToPhred(p));
    }

    #endregion

    #region Quality Encoding Round-Trip Tests

    [Test]
    public void EncodeQualityScores_Phred33_EncodesCorrectly()
    {
        // Q0 → '!' (ASCII 33), Q40 → 'I' (ASCII 73), Q93 → '~' (ASCII 126)
        // Evidence: Wikipedia FASTQ - Sanger encodes Q 0-93 using ASCII 33-126
        var encoded = FastqParser.EncodeQualityScores(new[] { 0, 40, 93 }, FastqParser.QualityEncoding.Phred33);
        Assert.That(encoded, Is.EqualTo("!I~"));
    }

    [Test]
    public void EncodeQualityScores_Phred64_EncodesCorrectly()
    {
        // Q0 → '@' (ASCII 64), Q40 → 'h' (ASCII 104), Q62 → '~' (ASCII 126)
        // Evidence: Wikipedia FASTQ - Phred+64 encodes Q 0-62 using ASCII 64-126
        var encoded = FastqParser.EncodeQualityScores(new[] { 0, 40, 62 }, FastqParser.QualityEncoding.Phred64);
        Assert.That(encoded, Is.EqualTo("@h~"));
    }

    [Test]
    public void EncodeDecodeRoundTrip_Phred33_PreservesScores()
    {
        // Round-trip must work across full Sanger range Q 0-93
        // Evidence: Wikipedia FASTQ - Sanger encodes Q 0-93 (ASCII 33-126)
        var originalScores = new[] { 0, 10, 20, 30, 40, 50, 60, 93 };
        var encoded = FastqParser.EncodeQualityScores(originalScores, FastqParser.QualityEncoding.Phred33);
        var decoded = FastqParser.DecodeQualityScores(encoded, FastqParser.QualityEncoding.Phred33);

        Assert.That(decoded, Is.EqualTo(originalScores));
    }

    [Test]
    public void EncodeDecodeRoundTrip_Phred64_PreservesScores()
    {
        // Round-trip must work across full Phred+64 range Q 0-62
        // Evidence: Wikipedia FASTQ - Phred+64 encodes Q 0-62 (ASCII 64-126)
        var originalScores = new[] { 0, 10, 20, 30, 40, 50, 62 };
        var encoded = FastqParser.EncodeQualityScores(originalScores, FastqParser.QualityEncoding.Phred64);
        var decoded = FastqParser.DecodeQualityScores(encoded, FastqParser.QualityEncoding.Phred64);

        Assert.That(decoded, Is.EqualTo(originalScores));
    }

    #endregion

    #region Filtering Tests

    [Test]
    public void FilterByQuality_FiltersLowQuality()
    {
        // read1: Q0 (Phred33 '!'), read2: Q40 (Phred33 'I'), read3: Q62 (Phred64 '~')
        var records = FastqParser.Parse(FastqWithVariousQuality).ToList();
        var filtered = FastqParser.FilterByQuality(records, 30).ToList();

        Assert.That(filtered, Has.Count.EqualTo(2));
        Assert.That(filtered[0].Id, Is.EqualTo("read2"));
        Assert.That(filtered[1].Id, Is.EqualTo("read3"));
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
        // Quality: Q0,Q0,Q40×8,Q0,Q0 → trim positions 0-1 and 10-11
        const string fastq = @"@read1
ACGTACGTACGT
+
!!IIIIIIII!!";

        var records = FastqParser.Parse(fastq).ToList();
        var trimmed = FastqParser.TrimByQuality(records[0], minQuality: 30);

        Assert.That(trimmed.Sequence, Is.EqualTo("GTACGTAC"));
        Assert.That(trimmed.QualityString, Is.EqualTo("IIIIIIII"));
    }

    [Test]
    public void TrimAdapter_RemovesAdapter()
    {
        // Adapter "AGATCGGAAGAG" starts at position 15 → keep first 15 bases
        const string adapter = "AGATCGGAAGAG";
        const string fastq = @"@read1
ACGTACGTACGTAAAAGATCGGAAGAG
+
IIIIIIIIIIIIIIIIIIIIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();
        var trimmed = FastqParser.TrimAdapter(records[0], adapter);

        Assert.That(trimmed.Sequence, Is.EqualTo("ACGTACGTACGTAAA"));
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
        // SimpleFastq: Q40 ('I') + Q39 ('H') at each of 16 positions
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var positionQuality = FastqParser.CalculatePositionQuality(records);

        Assert.That(positionQuality.Count, Is.EqualTo(16));

        // S3.2: Position numbering is 1-based
        Assert.That(positionQuality[0].Position, Is.EqualTo(1));
        Assert.That(positionQuality[15].Position, Is.EqualTo(16));

        // S3.1: Mean of [Q40, Q39] = 39.5
        Assert.That(positionQuality[0].MeanQuality, Is.EqualTo(39.5).Within(0.01));

        // S3.3: StdDev of [40, 39] = sqrt(((40-39.5)²+(39-39.5)²)/2) = 0.5
        Assert.That(positionQuality[0].StdDev, Is.EqualTo(0.5).Within(0.01));
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
        Assert.That(reads1[0].Sequence, Is.EqualTo("ACGTACGT"));
        Assert.That(reads1[1].Sequence, Is.EqualTo("AAAAAAAA"));
        Assert.That(reads2[0].Sequence, Is.EqualTo("TGCATGCA"));
        Assert.That(reads2[1].Sequence, Is.EqualTo("TTTTTTTT"));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void Parse_MultiplePlusLines_ParsesCorrectly()
    {
        // '+' in sequence data is allowed per Wikipedia FASTQ spec
        const string fastq = @"@read1
ACGT+ACGT
+
IIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].Sequence, Is.EqualTo("ACGT+ACGT"));
        Assert.That(records[0].QualityString, Is.EqualTo("IIIIIIIII"));
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
            Assert.That(records[0].Id, Is.EqualTo("SEQ_ID_1"));
            Assert.That(records[0].Sequence, Is.EqualTo("GATCGATCGATCGATC"));
            Assert.That(records[1].Id, Is.EqualTo("SEQ_ID_2"));
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

            var lines = File.ReadAllLines(tempFile);
            Assert.That(lines[0], Is.EqualTo("@SEQ_ID_1 description"));
            Assert.That(lines[1], Is.EqualTo("GATCGATCGATCGATC"));
            Assert.That(lines[2], Is.EqualTo("+"));
            Assert.That(lines[3], Is.EqualTo("IIIIIIIIIIIIIIII"));
            Assert.That(lines[4], Is.EqualTo("@SEQ_ID_2"));
            Assert.That(lines[5], Is.EqualTo("ACGTACGTACGTACGT"));
            Assert.That(lines[6], Is.EqualTo("+"));
            Assert.That(lines[7], Is.EqualTo("HHHHHHHHHHHHHHHH"));
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

        var lines = fastqString.TrimEnd().Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
        Assert.That(lines[0], Is.EqualTo("@test_id description"));
        Assert.That(lines[1], Is.EqualTo("ACGT"));
        Assert.That(lines[2], Is.EqualTo("+"));
        Assert.That(lines[3], Is.EqualTo("IIII"));
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
        // SimpleFastq: Q40 ('I') + Q39 ('H') — all bases ≥ Q20 → 100%
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.Q20Percentage, Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateStatistics_Q30Percentage_InValidRange()
    {
        // SimpleFastq: Q40 ('I') + Q39 ('H') — all bases ≥ Q30 → 100%
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.Q30Percentage, Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateStatistics_GcContent_InValidRange()
    {
        // SimpleFastq: "GATCGATCGATCGATC" (8 GC/16) + "ACGTACGTACGTACGT" (8 GC/16) = 50%
        var records = FastqParser.Parse(SimpleFastq).ToList();
        var stats = FastqParser.CalculateStatistics(records);

        Assert.That(stats.GcContent, Is.EqualTo(0.5));
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

    #region Missing Tests (S1.3, S2.1, S2.2, C1.2, C1.3)

    [Test]
    public void Parse_DescriptionWithSpecialCharacters_ParsedCorrectly()
    {
        // S1.3: Description can contain special characters
        const string fastq = "@read1 sample=A;lane=3;barcode=ACGT\nACGT\n+\nIIII";

        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records[0].Id, Is.EqualTo("read1"));
        Assert.That(records[0].Description, Is.EqualTo("sample=A;lane=3;barcode=ACGT"));
    }

    [Test]
    public void Parse_MultiLineSequence_AssembledCorrectly()
    {
        // S2.1: Multi-line sequence lines assembled into single sequence
        // Evidence: Wikipedia — Legacy Sanger files may split sequences across lines
        const string fastq = "@read1\nACGT\nTGCA\n+\nIIIIIIII";

        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].Sequence, Is.EqualTo("ACGTTGCA"));
    }

    [Test]
    public void Parse_MultiLineQuality_AssembledCorrectly()
    {
        // S2.2: Multi-line quality assembled to match sequence length
        const string fastq = "@read1\nACGTACGT\n+\nIIII\nHHHH";

        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].QualityString, Is.EqualTo("IIIIHHHH"));
    }

    [Test]
    public void Parse_VeryLongSequence_HandledCorrectly()
    {
        // C1.2: Very long sequences (10kb+) handled
        var sequence = new string('A', 10000) + new string('C', 10000);
        var quality = new string('I', 20000);
        var fastq = $"@long_read\n{sequence}\n+\n{quality}";

        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].Sequence, Has.Length.EqualTo(20000));
        Assert.That(records[0].Sequence, Is.EqualTo(sequence));
    }

    [Test]
    public void Parse_UnicodeInHeader_HandledGracefully()
    {
        // C1.3: Unicode characters in header handled gracefully
        const string fastq = "@read_g\u00E8ne_\u03B1 description_\u03B2\nACGT\n+\nIIII";

        var records = FastqParser.Parse(fastq).ToList();

        Assert.That(records, Has.Count.EqualTo(1));
        Assert.That(records[0].Id, Is.EqualTo("read_g\u00E8ne_\u03B1"));
        Assert.That(records[0].Description, Is.EqualTo("description_\u03B2"));
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

    #region Review 2026-09 (PARSE-FASTQ-001) — sourced regression tests

    // Reference values: Biopython 1.88 (Bio.SeqIO 'fastq' / 'fastq-illumina', FastqGeneralIterator),
    // cutadapt 5.2, FastQC PhredEncoding (lowest quality char over the whole file decides the offset).

    [Test]
    public void Parse_Auto_Phred64FileWithIllumina15BRead_DetectedFileLevel()
    {
        // F5: Illumina 1.5+ marks trimmed tails with 'B' (Q2). A read made only of 'B' lies in the
        // Phred+33/+64 overlap; per-record detection decoded it as Phred+33 Q33. The encoding is a file
        // property (FastQC): 'h' (104) in r1 ⇒ Phred+64 for every record.
        // Biopython 'fastq-illumina': r1 [40,40,40,40], r2 [2,2,2,2].
        const string fastq = "@r1\nACGT\n+\nhhhh\n@r2\nACGT\n+\nBBBB\n";
        var records = FastqParser.Parse(fastq).ToList();
        Assert.That(records[0].QualityScores, Is.EqualTo(new[] { 40, 40, 40, 40 }));
        Assert.That(records[1].QualityScores, Is.EqualTo(new[] { 2, 2, 2, 2 }));
    }

    [Test]
    public void Parse_Auto_Phred33FileWithQ41Read_DetectedFileLevel()
    {
        // F5: 'J' = Q41 (Illumina 1.8+ Phred+33 ceiling). Old per-record rule `c > 'I'` decoded "JJJJ"
        // as Phred+64 Q10. The '!' in r2 proves Phred+33 for the file.
        // Biopython 'fastq': r1 [41,41,41,41], r2 [0,0,0,0].
        const string fastq = "@r1\nACGT\n+\nJJJJ\n@r2\nACGT\n+\n!!!!\n";
        var records = FastqParser.Parse(fastq).ToList();
        Assert.That(records[0].QualityScores, Is.EqualTo(new[] { 41, 41, 41, 41 }));
        Assert.That(records[1].QualityScores, Is.EqualTo(new[] { 0, 0, 0, 0 }));
    }

    [Test]
    public void Parse_Auto_TextReaderAndFile_UseFileLevelEncoding()
    {
        const string fastq = "@r1\nACGT\n+\nhhhh\n@r2\nACGT\n+\nBBBB\n";
        using var reader = new StringReader(fastq);
        Assert.That(FastqParser.Parse(reader).Last().QualityScores, Is.EqualTo(new[] { 2, 2, 2, 2 }));

        string path = Path.GetTempFileName();
        try
        {
            File.WriteAllText(path, fastq);
            Assert.That(FastqParser.ParseFile(path).Last().QualityScores, Is.EqualTo(new[] { 2, 2, 2, 2 }));
        }
        finally { File.Delete(path); }
    }

    [Test]
    public void DetectEncoding_LowCharAfterHighChar_IsPhred33()
    {
        // F5: the old scan returned Phred+64 on the first char > 'I' ('J') before seeing '5' (53),
        // which cannot occur in Phred+64 (ASCII 64-126). Biopython 'fastq' "J5" ⇒ [41, 20].
        Assert.That(FastqParser.DetectEncoding("J5"), Is.EqualTo(FastqParser.QualityEncoding.Phred33));
        Assert.That(FastqParser.Parse("@x\nAC\n+\nJ5\n").Single().QualityScores, Is.EqualTo(new[] { 41, 20 }));
    }

    [Test]
    public void DecodeQualityScores_CharOutsideEncodingRange_Throws()
    {
        // F7: '5' (53) - 64 < 0 is not a Phred+64 symbol. Biopython 'fastq-illumina' raises
        // InvalidCharError for "J5hh"; the old code clamped it silently to Q0.
        Assert.Throws<ArgumentOutOfRangeException>(
            () => FastqParser.DecodeQualityScores("J5hh", FastqParser.QualityEncoding.Phred64));
        // An explicit Phred+64 parse of Phred+33 data is rejected rather than silently corrupted.
        Assert.Throws<FormatException>(
            () => FastqParser.Parse("@x\nACGT\n+\nJ5hh\n", FastqParser.QualityEncoding.Phred64).ToList());
        // Valid Phred+64 (Biopython 'fastq-illumina' "@h~~" ⇒ [0, 40, 62, 62]).
        Assert.That(FastqParser.DecodeQualityScores("@h~~", FastqParser.QualityEncoding.Phred64),
            Is.EqualTo(new[] { 0, 40, 62, 62 }));
    }

    [Test]
    public void Parse_BiopythonTrickyExample_ParsedByQualityLength()
    {
        // Biopython FastqGeneralIterator docstring ("Quality/tricky.fastq"): quality lines may start with
        // '@' or contain '+', the '+' line may repeat the title, and a record may be split over lines.
        const string tricky =
            "@071113_EAS56_0053:1:1:998:236\nTTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA\n+071113_EAS56_0053:1:1:998:236\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III\n" +
            "@071113_EAS56_0053:1:1:182:712\nACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG\n+\n@IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+\n" +
            "@071113_EAS56_0053:1:1:153:10\nTGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT\n+\nIIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6\n" +
            "@071113_EAS56_0053:1:3:990:501\nTGGGAGGTTTTATGTGGA\nAAGCAGCAATGTACAAGA\n+\nIIIIIII.IIIIII1@44\n@-7.%<&+/$/%4(++(%\n";
        var records = FastqParser.Parse(tricky).ToList();

        Assert.That(records, Has.Count.EqualTo(4));
        Assert.That(records[0].QualityString, Is.EqualTo("IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III"));
        Assert.That(records[1].QualityString, Is.EqualTo("@IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+"));
        Assert.That(records[1].QualityScores.Take(3), Is.EqualTo(new[] { 31, 40, 40 }));
        Assert.That(records[3].Id, Is.EqualTo("071113_EAS56_0053:1:3:990:501"));
        Assert.That(records[3].Sequence, Is.EqualTo("TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA"));
        Assert.That(records[3].QualityString, Is.EqualTo("IIIIIII.IIIIII1@44@-7.%<&+/$/%4(++(%"));
        Assert.That(records[3].QualityScores.TakeLast(5), Is.EqualTo(new[] { 7, 10, 10, 7, 4 }));
    }

    [Test]
    public void Parse_PlusCaptionDiffers_ThrowsFormatException()
    {
        // F6: Biopython: "Sequence and quality captions differ." (Cock et al. 2010: optional repeat must match).
        Assert.Throws<FormatException>(() => FastqParser.Parse("@r1\nACGT\n+r2\nIIII\n").ToList());
        // A matching repeated caption is accepted.
        Assert.That(FastqParser.Parse("@r1 d\nACGT\n+r1 d\nIIII\n").Single().Sequence, Is.EqualTo("ACGT"));
    }

    [Test]
    public void Parse_WhitespaceInSequence_ThrowsFormatException()
    {
        // F6: Biopython: "Whitespace is not allowed in the sequence."
        Assert.Throws<FormatException>(() => FastqParser.Parse("@r1\nAC GT\n+\nIIII\n").ToList());
    }

    [Test]
    public void Parse_QualityLengthMismatch_ThrowsFormatException()
    {
        // F6 / INV-01: previously yielded records whose quality was shorter than the sequence, which then
        // crashed TrimByQuality with IndexOutOfRange. Biopython: "... differs for rec (8 and 2)."
        var ex = Assert.Throws<FormatException>(() => FastqParser.Parse("@rec\nACGTACGT\n+\nII\n").ToList());
        Assert.That(ex!.Message, Does.Contain("(8 and 2)"));
    }

    [Test]
    public void Parse_ZeroLengthRecordAndBlankSeparators_Accepted()
    {
        // Biopython: "@r1\n\n+\n\n@r2\nA\n+\nI\n" ⇒ [('r1','',''), ('r2','A','I')];
        // blank lines between records are tolerated.
        var records = FastqParser.Parse("@r1\n\n+\n\n@r2\nA\n+\nI\n\n\n").ToList();
        Assert.That(records.Select(r => r.Id), Is.EqualTo(new[] { "r1", "r2" }));
        Assert.That(records[0].Sequence, Is.Empty);
        Assert.That(records[1].QualityScores, Is.EqualTo(new[] { 40 }));
    }

    [Test]
    public void Parse_TabInTitle_SplitsIdOnAnyWhitespace()
    {
        // F11: Biopython SeqIO 'fastq' "@r1\tdesc here" ⇒ id 'r1' (title.split(None, 1)[0]).
        var rec = FastqParser.Parse("@r1\tdesc here\nACGT\n+\nIIII\n").Single();
        Assert.That(rec.Id, Is.EqualTo("r1"));
        Assert.That(rec.Description, Is.EqualTo("desc here"));
    }

    private static FastqParser.FastqRecord Rec(string seq) =>
        new("r", "", seq, new string('I', seq.Length), Enumerable.Repeat(40, seq.Length).ToList());

    [TestCase("AGATCGGAAGAGCACACGTC", "")]                    // full adapter at position 0 ⇒ empty read
    [TestCase("ACGTACGTAGATCGGAAGAGCTTTTTAGATC", "ACGTACGT")] // internal full beats 3' partial
    [TestCase("CCCCAGATCGGAAGAGCGGGGAGATCGGAAGAGCTT", "CCCC")] // leftmost of two full occurrences
    [TestCase("ACGTACGTACGTAAAAGATCG", "ACGTACGTACGTAAA")]     // 3' partial (6 >= minOverlap 5)
    [TestCase("ACGTACGTACGTAGATcggaa", "ACGTACGTACGT")]        // case-insensitive partial
    [TestCase("ACGTACGTACGTAGAT", "ACGTACGTACGTAGAT")]         // partial 4 < 5 ⇒ unchanged
    [TestCase("AGATCGGAAG", "")]                              // read shorter than adapter, all adapter
    public void TrimAdapter_MatchesCutadapt(string read, string expected)
    {
        // F9: cutadapt 5.2 `-a AGATCGGAAGAGC -e 0 -O 5` output for the same reads.
        var trimmed = FastqParser.TrimAdapter(Rec(read), "AGATCGGAAGAGC", 5);
        Assert.That(trimmed.Sequence, Is.EqualTo(expected));
        Assert.That(trimmed.QualityString, Has.Length.EqualTo(expected.Length));
        Assert.That(trimmed.QualityScores, Has.Count.EqualTo(expected.Length));
    }

    [Test]
    public void WriteToFile_NoByteOrderMark_FirstByteIsAt()
    {
        // F10: Encoding.UTF8 wrote EF BB BF before '@'; Biopython ("Records in Fastq files should start
        // with '@' character") and cutadapt ("Input file format not recognized") both reject such a file.
        string path = Path.GetTempFileName();
        try
        {
            FastqParser.WriteToFile(path, FastqParser.Parse("@r1\nACGT\n+\nIIII\n"));
            byte[] bytes = File.ReadAllBytes(path);
            Assert.That(bytes[0], Is.EqualTo((byte)'@'));
        }
        finally { File.Delete(path); }
    }

    #endregion
}
