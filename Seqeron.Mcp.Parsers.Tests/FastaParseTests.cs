using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class FastaParseTests
{
    [Test]
    public void FastaParse_Schema_ValidatesCorrectly()
    {
        var validFasta = ">seq1 test\nATGCATGC\n>seq2\nGGGCCC";
        Assert.DoesNotThrow(() => ParsersTools.FastaParse(validFasta));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastaParse(""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastaParse(null!));
    }

    [Test]
    public void FastaParse_Binding_InvokesSuccessfully()
    {
        var fasta = ">seq1 Human gene\nATGCATGCATGC\n>seq2 Mouse gene\nGGGCCCAAATTT";
        var result = ParsersTools.FastaParse(fasta);

        Assert.That(result.Count, Is.EqualTo(2));
        Assert.That(result.Entries[0].Id, Is.EqualTo("seq1"));
        Assert.That(result.Entries[0].Description, Is.EqualTo("Human gene"));
        Assert.That(result.Entries[0].Sequence, Is.EqualTo("ATGCATGCATGC"));
        Assert.That(result.Entries[0].Length, Is.EqualTo(12));
        Assert.That(result.Entries[1].Id, Is.EqualTo("seq2"));
    }
}
