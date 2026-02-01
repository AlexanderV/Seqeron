using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class FastqStatisticsTests
{
    [Test]
    public void FastqStatistics_Schema_ValidatesCorrectly()
    {
        var validFastq = "@seq1\nATGC\n+\nIIII";
        Assert.DoesNotThrow(() => ParsersTools.FastqStatistics(validFastq));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqStatistics(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqStatistics(null!));
    }

    [Test]
    public void FastqStatistics_Binding_InvokesSuccessfully()
    {
        var fastq = "@seq1\nATGCATGC\n+\nIIIIIIII\n@seq2\nGGGCCC\n+\nHHHHHH";
        var result = ParsersTools.FastqStatistics(fastq);

        Assert.That(result.TotalReads, Is.EqualTo(2));
        Assert.That(result.TotalBases, Is.EqualTo(14)); // 8 + 6
        Assert.That(result.MeanReadLength, Is.EqualTo(7.0));
        Assert.That(result.MinReadLength, Is.EqualTo(6));
        Assert.That(result.MaxReadLength, Is.EqualTo(8));
        Assert.That(result.MeanQuality, Is.GreaterThan(30)); // High quality scores
        Assert.That(result.Q20Percentage, Is.EqualTo(100.0));
        Assert.That(result.Q30Percentage, Is.EqualTo(100.0));
        Assert.That(result.GcContent, Is.GreaterThan(0).And.LessThan(1));
    }

    [Test]
    public void FastqStatistics_Encoding_ParsesPhred33()
    {
        var fastq = "@test\nATGC\n+\nIIII";
        var result = ParsersTools.FastqStatistics(fastq, "phred33");

        Assert.That(result.TotalReads, Is.EqualTo(1));
        Assert.That(result.MeanQuality, Is.EqualTo(40)); // 'I' = Q40
    }

    [Test]
    public void FastqStatistics_Encoding_InvalidEncodingThrows()
    {
        var fastq = "@test\nATGC\n+\nIIII";
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqStatistics(fastq, "invalid"));
    }

    [Test]
    public void FastqStatistics_GcContent_CalculatesCorrectly()
    {
        // Sequence with 50% GC (GCGC out of ATGCGCAT)
        var fastq = "@test\nGCGCATAT\n+\nIIIIIIII";
        var result = ParsersTools.FastqStatistics(fastq);

        Assert.That(result.GcContent, Is.EqualTo(0.5).Within(0.01)); // 50% GC
    }
}
