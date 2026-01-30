using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplexityMaskLowTests
{
    [Test]
    public void ComplexityMaskLow_Schema_ValidatesCorrectly()
    {
        // Need a sequence longer than window size
        var longSeq = "ATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCG";
        Assert.DoesNotThrow(() => SequenceTools.ComplexityMaskLow(longSeq, 64, 2.0, 'N'));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityMaskLow("", 64, 2.0, 'N'));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityMaskLow(null!, 64, 2.0, 'N'));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityMaskLow(longSeq, 0, 2.0, 'N')); // windowSize < 1
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityMaskLow("ZZZZATGC", 4, 2.0, 'N')); // invalid DNA
    }

    [Test]
    public void ComplexityMaskLow_Binding_InvokesSuccessfully()
    {
        var seq = "ATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCGATGCGATCGATCG";
        var result = SequenceTools.ComplexityMaskLow(seq, 64, 2.0, 'N');
        Assert.That(result.MaskedSequence, Is.Not.Null.And.Not.Empty);
        Assert.That(result.OriginalLength, Is.EqualTo(seq.Length));
        Assert.That(result.MaskChar, Is.EqualTo('N'));
    }
}
