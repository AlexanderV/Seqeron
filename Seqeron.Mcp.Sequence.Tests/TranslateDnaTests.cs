using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class TranslateDnaTests
{
    [Test]
    public void TranslateDna_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.TranslateDna("ATGCGA", 0, false));
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("", 0, false));
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna(null!, 0, false));
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ATGCGA", -1, false)); // invalid frame
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ATGCGA", 3, false)); // invalid frame
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateDna("ZZZCGA", 0, false)); // invalid DNA
    }

    [Test]
    public void TranslateDna_Binding_InvokesSuccessfully()
    {
        // ATG = Met (M), CGA = Arg (R)
        var result = SequenceTools.TranslateDna("ATGCGA", 0, false);
        Assert.That(result.Protein, Is.EqualTo("MR"));
        Assert.That(result.Frame, Is.EqualTo(0));
        Assert.That(result.DnaLength, Is.EqualTo(6));

        // Test frame 1
        var result2 = SequenceTools.TranslateDna("AATGCGA", 1, false);
        Assert.That(result2.Frame, Is.EqualTo(1));
    }
}
