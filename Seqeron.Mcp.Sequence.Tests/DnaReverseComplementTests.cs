using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class DnaReverseComplementTests
{
    [Test]
    public void DnaReverseComplement_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.DnaReverseComplement("ATGC"));
        Assert.Throws<ArgumentException>(() => SequenceTools.DnaReverseComplement(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.DnaReverseComplement(null!));
        Assert.Throws<ArgumentException>(() => SequenceTools.DnaReverseComplement("ATGX")); // Invalid DNA
    }

    [Test]
    public void DnaReverseComplement_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.DnaReverseComplement("ATGC");
        Assert.That(result.ReverseComplement, Is.EqualTo("GCAT"));

        var selfComplement = SequenceTools.DnaReverseComplement("AATT");
        Assert.That(selfComplement.ReverseComplement, Is.EqualTo("AATT"));

        var single = SequenceTools.DnaReverseComplement("A");
        Assert.That(single.ReverseComplement, Is.EqualTo("T"));
    }
}
