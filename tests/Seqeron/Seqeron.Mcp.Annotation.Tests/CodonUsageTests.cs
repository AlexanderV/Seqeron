using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CodonUsageTests
{
    [Test]
    public void CodonUsage_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CodonUsage("ATGAAATAA"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CodonUsage(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CodonUsage(null!));
    }

    [Test]
    public void CodonUsage_CountsInFrameCodons()
    {
        // GetCodonUsage(string) counts in-frame unambiguous triplets: ATG, ATG, AAA, TAA.
        var result = AnnotationTools.CodonUsage("ATGATGAAATAA");

        Assert.Multiple(() =>
        {
            Assert.That(result.Usage["ATG"], Is.EqualTo(2));
            Assert.That(result.Usage["AAA"], Is.EqualTo(1));
            Assert.That(result.Usage["TAA"], Is.EqualTo(1));
            Assert.That(result.Usage.Values.Sum(), Is.EqualTo(4));
        });
    }

    [Test]
    public void CodonUsage_IgnoresPartialTrailingCodon()
    {
        // "ATGAT" has one full codon ATG; the trailing "AT" is dropped.
        var result = AnnotationTools.CodonUsage("ATGAT");

        Assert.Multiple(() =>
        {
            Assert.That(result.Usage["ATG"], Is.EqualTo(1));
            Assert.That(result.Usage.Values.Sum(), Is.EqualTo(1));
        });
    }

    [Test]
    public void CodonUsage_RnaInput_UReadAsT_AmbiguousTripletSkipped()
    {
        // Delegates to Translator.CountCodons: AUG → ATG (CodonW ident_codon), NNN skipped, frame kept.
        var result = AnnotationTools.CodonUsage("AUGNNNAAA");

        Assert.That(result.Usage, Is.EquivalentTo(new Dictionary<string, int> { ["ATG"] = 1, ["AAA"] = 1 }));
    }

    [Test]
    public void CodonUsage_IsCaseInsensitive()
    {
        var result = AnnotationTools.CodonUsage("atgatg");
        Assert.That(result.Usage["ATG"], Is.EqualTo(2));
    }
}
