using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CreateMiRnaTests
{
    [Test]
    public void CreateMiRna_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CreateMiRna("let-7a", "UGAGGUAGUAGGUUGUAUAGUU"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CreateMiRna("", "UGAGGUAG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CreateMiRna("x", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CreateMiRna(null!, "UGAGGUAG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CreateMiRna("x", null!));
    }

    [Test]
    public void CreateMiRna_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.CreateMiRna upper-cases and T->U normalises the sequence, then extracts the
        // positions-2-8 seed with SeedStart=1, SeedEnd=7.
        // DNA "TGAGGTAGTAGGTTGTATAGTT" -> RNA "UGAGGUAGUAGGUUGUAUAGUU", seed "GAGGUAG".
        var result = AnnotationTools.CreateMiRna("let-7a", "TGAGGTAGTAGGTTGTATAGTT");

        Assert.Multiple(() =>
        {
            Assert.That(result.MiRna.Name, Is.EqualTo("let-7a"));
            Assert.That(result.MiRna.Sequence, Is.EqualTo("UGAGGUAGUAGGUUGUAUAGUU"));
            Assert.That(result.MiRna.SeedSequence, Is.EqualTo("GAGGUAG"));
            Assert.That(result.MiRna.SeedStart, Is.EqualTo(1));
            Assert.That(result.MiRna.SeedEnd, Is.EqualTo(7));
        });
    }
}
