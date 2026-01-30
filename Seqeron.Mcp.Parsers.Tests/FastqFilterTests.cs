using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class FastqFilterTests
{
    [Test]
    public void FastqFilter_Schema_ValidatesCorrectly()
    {
        var validFastq = "@seq1\nATGC\n+\nIIII";
        Assert.DoesNotThrow(() => ParsersTools.FastqFilter(validFastq, 20));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqFilter("", 20));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqFilter(null!, 20));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqFilter(validFastq, -1));
    }

    [Test]
    public void FastqFilter_Binding_FiltersCorrectly()
    {
        // High quality read (Q40) and low quality read (Q0)
        var fastq = "@high_quality\nATGC\n+\nIIII\n@low_quality\nATGC\n+\n!!!!";
        var result = ParsersTools.FastqFilter(fastq, 30, "phred33");

        Assert.That(result.TotalCount, Is.EqualTo(2));
        Assert.That(result.PassedCount, Is.EqualTo(1));
        Assert.That(result.Entries[0].Id, Is.EqualTo("high_quality"));
        Assert.That(result.PassedPercentage, Is.EqualTo(50.0));
    }

    [Test]
    public void FastqFilter_Binding_AllPass()
    {
        var fastq = "@seq1\nATGC\n+\nIIII\n@seq2\nGGGG\n+\nHHHH";
        var result = ParsersTools.FastqFilter(fastq, 30);

        Assert.That(result.PassedCount, Is.EqualTo(2));
        Assert.That(result.TotalCount, Is.EqualTo(2));
        Assert.That(result.PassedPercentage, Is.EqualTo(100.0));
    }

    [Test]
    public void FastqFilter_Binding_NonePass()
    {
        // Very low quality reads
        var fastq = "@seq1\nATGC\n+\n!!!!\n@seq2\nGGGG\n+\n!!!!";
        var result = ParsersTools.FastqFilter(fastq, 30, "phred33");

        Assert.That(result.PassedCount, Is.EqualTo(0));
        Assert.That(result.TotalCount, Is.EqualTo(2));
        Assert.That(result.PassedPercentage, Is.EqualTo(0.0));
    }

    [Test]
    public void FastqFilter_Encoding_InvalidEncodingThrows()
    {
        var fastq = "@test\nATGC\n+\nIIII";
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqFilter(fastq, 20, "invalid"));
    }
}
