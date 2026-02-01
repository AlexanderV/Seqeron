using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class RnaFromDnaTests
{
    [Test]
    public void RnaFromDna_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.RnaFromDna("ATGC"));
        Assert.Throws<ArgumentException>(() => SequenceTools.RnaFromDna(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.RnaFromDna(null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.RnaFromDna("ATGX")); // Invalid DNA
    }

    [Test]
    public void RnaFromDna_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.RnaFromDna("ATGC");
        Assert.That(result.Rna, Is.EqualTo("AUGC"));

        var allT = SequenceTools.RnaFromDna("TTTT");
        Assert.That(allT.Rna, Is.EqualTo("UUUU"));

        var noT = SequenceTools.RnaFromDna("ACGACG");
        Assert.That(noT.Rna, Is.EqualTo("ACGACG"));
    }
}
