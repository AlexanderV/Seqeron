using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class MiRnaSeedSequenceTests
{
    // let-7a-5p mature sequence.
    private const string Let7a = "UGAGGUAGUAGGUUGUAUAGUU";

    [Test]
    public void MiRnaSeedSequence_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.MiRnaSeedSequence(Let7a));
        Assert.Throws<ArgumentException>(() => AnnotationTools.MiRnaSeedSequence(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.MiRnaSeedSequence(null!));
        // Fewer than 8 nt cannot yield a positions-2-8 seed.
        Assert.Throws<ArgumentException>(() => AnnotationTools.MiRnaSeedSequence("UGAGGU"));
    }

    [Test]
    public void MiRnaSeedSequence_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.GetSeedSequence returns positions 2-8 (0-based Substring(1,7)), upper-cased.
        // let-7a: positions 2-8 of UGAGGUAG... = GAGGUAG.
        var result = AnnotationTools.MiRnaSeedSequence(Let7a);
        Assert.That(result.Seed, Is.EqualTo("GAGGUAG"));
    }

    [Test]
    public void MiRnaSeedSequence_UpperCasesInput()
    {
        // Lower-case input is upper-cased; bases are not T->U converted here.
        var result = AnnotationTools.MiRnaSeedSequence("ugagguaguag");
        Assert.That(result.Seed, Is.EqualTo("GAGGUAG"));
    }
}
