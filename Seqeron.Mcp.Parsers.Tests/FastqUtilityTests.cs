using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class FastqUtilityTests
{
    // ========================
    // fastq_detect_encoding Tests
    // ========================

    [Test]
    public void FastqDetectEncoding_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.FastqDetectEncoding("IIII"));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqDetectEncoding(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqDetectEncoding(null!));
    }

    [Test]
    public void FastqDetectEncoding_Binding_DetectsPhred33()
    {
        // Characters below '@' (64) indicate Phred+33
        var result = ParsersTools.FastqDetectEncoding("!\"#$%&'()*+,-./0123456789:");

        Assert.That(result.Encoding, Is.EqualTo("Phred33"));
        Assert.That(result.Offset, Is.EqualTo(33));
    }

    [Test]
    public void FastqDetectEncoding_Binding_DetectsPhred64()
    {
        // Characters above 'I' (73) indicate Phred+64
        var result = ParsersTools.FastqDetectEncoding("efghijklmnopqrstuvwxyz");

        Assert.That(result.Encoding, Is.EqualTo("Phred64"));
        Assert.That(result.Offset, Is.EqualTo(64));
    }

    // ========================
    // fastq_encode_quality Tests
    // ========================

    [Test]
    public void FastqEncodeQuality_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.FastqEncodeQuality(new List<int> { 20, 30, 40 }));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqEncodeQuality(new List<int>()));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqEncodeQuality(null!));
    }

    [Test]
    public void FastqEncodeQuality_Binding_EncodesPhred33()
    {
        var scores = new List<int> { 0, 10, 20, 30, 40 };
        var result = ParsersTools.FastqEncodeQuality(scores, "phred33");

        Assert.That(result.Length, Is.EqualTo(5));
        Assert.That(result.QualityString, Is.EqualTo("!+5?I")); // 0+33='!', 10+33='+', etc.
    }

    [Test]
    public void FastqEncodeQuality_Binding_EncodesPhred64()
    {
        var scores = new List<int> { 0, 10, 20, 30, 40 };
        var result = ParsersTools.FastqEncodeQuality(scores, "phred64");

        Assert.That(result.Length, Is.EqualTo(5));
        Assert.That(result.QualityString, Is.EqualTo("@JT^h")); // 0+64='@', 10+64='J', etc.
    }

    [Test]
    public void FastqEncodeQuality_Binding_InvalidEncodingThrows()
    {
        var scores = new List<int> { 20, 30 };
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqEncodeQuality(scores, "invalid"));
    }

    // ========================
    // fastq_phred_to_error Tests
    // ========================

    [Test]
    public void FastqPhredToError_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.FastqPhredToError(20));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqPhredToError(-1));
    }

    [Test]
    public void FastqPhredToError_Binding_ConvertsCorrectly()
    {
        // Q=10 -> P=0.1
        var result10 = ParsersTools.FastqPhredToError(10);
        Assert.That(result10.PhredScore, Is.EqualTo(10));
        Assert.That(result10.ErrorProbability, Is.EqualTo(0.1).Within(0.0001));

        // Q=20 -> P=0.01
        var result20 = ParsersTools.FastqPhredToError(20);
        Assert.That(result20.ErrorProbability, Is.EqualTo(0.01).Within(0.0001));

        // Q=30 -> P=0.001
        var result30 = ParsersTools.FastqPhredToError(30);
        Assert.That(result30.ErrorProbability, Is.EqualTo(0.001).Within(0.0001));
    }

    // ========================
    // fastq_error_to_phred Tests
    // ========================

    [Test]
    public void FastqErrorToPhred_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.FastqErrorToPhred(0.01));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqErrorToPhred(-0.1));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqErrorToPhred(1.5));
    }

    [Test]
    public void FastqErrorToPhred_Binding_ConvertsCorrectly()
    {
        // P=0.1 -> Q=10
        var result01 = ParsersTools.FastqErrorToPhred(0.1);
        Assert.That(result01.ErrorProbability, Is.EqualTo(0.1));
        Assert.That(result01.PhredScore, Is.EqualTo(10));

        // P=0.01 -> Q=20
        var result001 = ParsersTools.FastqErrorToPhred(0.01);
        Assert.That(result001.PhredScore, Is.EqualTo(20));

        // P=0.001 -> Q=30
        var result0001 = ParsersTools.FastqErrorToPhred(0.001);
        Assert.That(result0001.PhredScore, Is.EqualTo(30));
    }

    // ========================
    // fastq_trim_quality Tests
    // ========================

    [Test]
    public void FastqTrimQuality_Schema_ValidatesCorrectly()
    {
        var fastq = "@seq1\nATGCATGC\n+\nIIIIIIII";
        Assert.DoesNotThrow(() => ParsersTools.FastqTrimQuality(fastq));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqTrimQuality(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqTrimQuality(null!));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqTrimQuality(fastq, -1));
    }

    [Test]
    public void FastqTrimQuality_Binding_TrimsLowQualityEnds()
    {
        // '!' = Q0, 'I' = Q40
        var fastq = "@seq1\nAAAATGCAAAA\n+\n!!!!!IIII!!";
        var result = ParsersTools.FastqTrimQuality(fastq, 20);

        Assert.That(result.Count, Is.EqualTo(1));
        Assert.That(result.Entries[0].Sequence.Length, Is.LessThan(11));
        Assert.That(result.OriginalBases, Is.EqualTo(11));
        Assert.That(result.TrimmedBases, Is.LessThan(11));
        Assert.That(result.TrimmedPercentage, Is.GreaterThan(0));
    }

    [Test]
    public void FastqTrimQuality_Binding_HighQualityUnchanged()
    {
        var fastq = "@seq1\nATGCATGC\n+\nIIIIIIII";
        var result = ParsersTools.FastqTrimQuality(fastq, 20);

        Assert.That(result.Count, Is.EqualTo(1));
        Assert.That(result.Entries[0].Sequence, Is.EqualTo("ATGCATGC"));
        Assert.That(result.TrimmedPercentage, Is.EqualTo(0));
    }

    // ========================
    // fastq_trim_adapter Tests
    // ========================

    [Test]
    public void FastqTrimAdapter_Schema_ValidatesCorrectly()
    {
        var fastq = "@seq1\nATGCATGC\n+\nIIIIIIII";
        Assert.DoesNotThrow(() => ParsersTools.FastqTrimAdapter(fastq, "ATGC"));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqTrimAdapter("", "ATGC"));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqTrimAdapter(fastq, ""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqTrimAdapter(fastq, "ATGC", 0));
    }

    [Test]
    public void FastqTrimAdapter_Binding_TrimsAdapter()
    {
        var fastq = "@seq1\nATGCATGCAAAAAAAA\n+\nIIIIIIIIIIIIIIII";
        var result = ParsersTools.FastqTrimAdapter(fastq, "AAAAAAAA", 5);

        Assert.That(result.Count, Is.EqualTo(1));
        Assert.That(result.Entries[0].Sequence, Is.EqualTo("ATGCATGC"));
        Assert.That(result.ReadsWithAdapter, Is.EqualTo(1));
        Assert.That(result.OriginalBases, Is.EqualTo(16));
        Assert.That(result.TrimmedBases, Is.EqualTo(8));
    }

    [Test]
    public void FastqTrimAdapter_Binding_NoAdapterUnchanged()
    {
        var fastq = "@seq1\nATGCATGC\n+\nIIIIIIII";
        var result = ParsersTools.FastqTrimAdapter(fastq, "GGGGGGGG", 5);

        Assert.That(result.Count, Is.EqualTo(1));
        Assert.That(result.Entries[0].Sequence, Is.EqualTo("ATGCATGC"));
        Assert.That(result.ReadsWithAdapter, Is.EqualTo(0));
    }

    // ========================
    // fastq_format Tests
    // ========================

    [Test]
    public void FastqFormat_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.FastqFormat("seq1", "ATGC", "IIII"));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqFormat("", "ATGC", "IIII"));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqFormat("seq1", "", "IIII"));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqFormat("seq1", "ATGC", ""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqFormat("seq1", "ATGC", "II")); // length mismatch
    }

    [Test]
    public void FastqFormat_Binding_FormatsCorrectly()
    {
        var result = ParsersTools.FastqFormat("seq1", "ATGCATGC", "IIIIIIII", "Test sequence");

        Assert.That(result.Fastq, Does.StartWith("@seq1 Test sequence"));
        Assert.That(result.Fastq, Does.Contain("ATGCATGC"));
        Assert.That(result.Fastq, Does.Contain("+"));
        Assert.That(result.Fastq, Does.Contain("IIIIIIII"));
    }

    [Test]
    public void FastqFormat_Binding_WithoutDescription()
    {
        var result = ParsersTools.FastqFormat("seq1", "ATGC", "IIII");

        Assert.That(result.Fastq, Does.StartWith("@seq1"));
        Assert.That(result.Fastq, Does.Not.Contain("seq1 ")); // No space after ID without description
    }

    // ========================
    // fastq_write Tests
    // ========================

    [Test]
    public void FastqWrite_Schema_ValidatesCorrectly()
    {
        var fastq = "@seq1\nATGC\n+\nIIII";
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqWrite("", fastq));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqWrite("test.fq", ""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqWrite("test.fq", fastq, "invalid"));
    }

    [Test]
    public void FastqWrite_Binding_WritesToFile()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            var fastq = "@seq1\nATGCATGC\n+\nIIIIIIII\n@seq2\nGGGCCC\n+\nHHHHHH";
            var result = ParsersTools.FastqWrite(tempFile, fastq);

            Assert.That(result.FilePath, Is.EqualTo(tempFile));
            Assert.That(result.RecordsWritten, Is.EqualTo(2));
            Assert.That(result.TotalBases, Is.EqualTo(14));
            Assert.That(File.Exists(tempFile), Is.True);

            var content = File.ReadAllText(tempFile);
            Assert.That(content, Does.Contain("@seq1"));
            Assert.That(content, Does.Contain("ATGCATGC"));
        }
        finally
        {
            if (File.Exists(tempFile))
                File.Delete(tempFile);
        }
    }
}
