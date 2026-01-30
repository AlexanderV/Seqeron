using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class FastqParseTests
{
    [Test]
    public void FastqParse_Schema_ValidatesCorrectly()
    {
        var validFastq = "@seq1\nATGC\n+\nIIII";
        Assert.DoesNotThrow(() => ParsersTools.FastqParse(validFastq));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqParse(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqParse(null!));
    }

    [Test]
    public void FastqParse_Binding_InvokesSuccessfully()
    {
        var fastq = "@seq1 Human gene\nATGCATGC\n+\nIIIIIIII\n@seq2\nGGGCCC\n+\nHHHHHH";
        var result = ParsersTools.FastqParse(fastq);

        Assert.That(result.Count, Is.EqualTo(2));
        Assert.That(result.Entries[0].Id, Is.EqualTo("seq1"));
        Assert.That(result.Entries[0].Description, Is.EqualTo("Human gene"));
        Assert.That(result.Entries[0].Sequence, Is.EqualTo("ATGCATGC"));
        Assert.That(result.Entries[0].Length, Is.EqualTo(8));
        Assert.That(result.Entries[0].QualityString, Is.EqualTo("IIIIIIII"));
        Assert.That(result.Entries[1].Id, Is.EqualTo("seq2"));
    }

    [Test]
    public void FastqParse_Encoding_ParsesPhred33()
    {
        var fastq = "@test\nATGC\n+\nIIII";
        var result = ParsersTools.FastqParse(fastq, "phred33");

        Assert.That(result.Count, Is.EqualTo(1));
        Assert.That(result.Entries[0].QualityScores, Has.All.EqualTo(40)); // 'I' = 73, 73-33 = 40
    }

    [Test]
    public void FastqParse_Encoding_InvalidEncodingThrows()
    {
        var fastq = "@test\nATGC\n+\nIIII";
        Assert.Throws<ArgumentException>(() => ParsersTools.FastqParse(fastq, "invalid"));
    }
}
