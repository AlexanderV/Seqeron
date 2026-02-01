using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class TranslateRnaTests
{
    [Test]
    public void TranslateRna_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.TranslateRna("AUGCGA", 0, false));
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateRna("", 0, false));
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateRna(null!, 0, false));
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateRna("AUGCGA", -1, false)); // invalid frame
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateRna("AUGCGA", 3, false)); // invalid frame
        Assert.Throws<ArgumentException>(() => SequenceTools.TranslateRna("ATGCGA", 0, false)); // T is not valid in RNA
    }

    [Test]
    public void TranslateRna_Binding_InvokesSuccessfully()
    {
        // AUG = Met (M), CGA = Arg (R)
        var result = SequenceTools.TranslateRna("AUGCGA", 0, false);
        Assert.That(result.Protein, Is.EqualTo("MR"));
        Assert.That(result.Frame, Is.EqualTo(0));
        Assert.That(result.RnaLength, Is.EqualTo(6));

        // Test frame 1
        var result2 = SequenceTools.TranslateRna("AAUGCGA", 1, false);
        Assert.That(result2.Frame, Is.EqualTo(1));
    }
}
