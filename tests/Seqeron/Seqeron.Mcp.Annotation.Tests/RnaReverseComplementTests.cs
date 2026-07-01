using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class RnaReverseComplementTests
{
    [Test]
    public void RnaReverseComplement_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.RnaReverseComplement("AUGC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.RnaReverseComplement(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.RnaReverseComplement(null!));
    }

    [Test]
    public void RnaReverseComplement_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.GetReverseComplement: A<->U, G<->C, reversed.
        // AUGC reversed-complemented -> GCAU.
        Assert.That(AnnotationTools.RnaReverseComplement("AUGC").ReverseComplement, Is.EqualTo("GCAU"));
        // let-7a seed GAGGUAG -> CUACCUC (the seed reverse complement used in target scanning).
        Assert.That(AnnotationTools.RnaReverseComplement("GAGGUAG").ReverseComplement, Is.EqualTo("CUACCUC"));
    }

    [Test]
    public void RnaReverseComplement_DnaTIsComplementedAsU_UnknownBecomesN()
    {
        // DNA T is treated like U (complement A); unknown bases map to N.
        Assert.That(AnnotationTools.RnaReverseComplement("ATGC").ReverseComplement, Is.EqualTo("GCAU"));
        Assert.That(AnnotationTools.RnaReverseComplement("AXGC").ReverseComplement, Is.EqualTo("GCNU"));
    }
}
