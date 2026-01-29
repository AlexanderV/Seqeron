using NUnit.Framework;
using SuffixTree.Mcp.Parsers.Tools;

namespace SuffixTree.Mcp.Parsers.Tests;

[TestFixture]
public class FastaFormatTests
{
    [Test]
    public void FastaFormat_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.FastaFormat("seq1", "ATGCATGC"));
        Assert.DoesNotThrow(() => ParsersTools.FastaFormat("seq1", "ATGCATGC", "description", 60));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastaFormat("", "ATGC"));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastaFormat(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastaFormat("seq1", ""));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastaFormat("seq1", null!));
        Assert.Throws<ArgumentException>(() => ParsersTools.FastaFormat("seq1", "ATGC", null, 0));
    }

    [Test]
    public void FastaFormat_Binding_InvokesSuccessfully()
    {
        var result = ParsersTools.FastaFormat("seq1", "ATGCATGCATGC", "Human gene");

        Assert.That(result.Fasta, Does.StartWith(">seq1 Human gene"));
        Assert.That(result.Fasta, Does.Contain("ATGCATGCATGC"));

        // Test without description
        var result2 = ParsersTools.FastaFormat("NM_001", "GGGGCCCC");
        Assert.That(result2.Fasta, Does.StartWith(">NM_001"));
        Assert.That(result2.Fasta, Does.Contain("GGGGCCCC"));

        // Test line wrapping
        var longSeq = new string('A', 100);
        var wrapped = ParsersTools.FastaFormat("test", longSeq, null, 50);
        var lines = wrapped.Fasta.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        Assert.That(lines.Length, Is.EqualTo(3)); // header + 2 sequence lines
    }
}
