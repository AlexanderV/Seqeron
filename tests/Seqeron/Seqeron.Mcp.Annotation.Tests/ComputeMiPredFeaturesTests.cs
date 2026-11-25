using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ComputeMiPredFeaturesTests
{
    [Test]
    public void ComputeMiPredFeatures_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeMiPredFeatures("", "(((...)))", -3, 0.5));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeMiPredFeatures("GGGAAACCC", "", -3, 0.5));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeMiPredFeatures("GGGAAACCC", "(((...))", -3, 0.5)); // length mismatch (library)
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ComputeMiPredFeatures("GGGAAACCC", "(((...)))", double.NaN, 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ComputeMiPredFeatures("GGGAAACCC", "(((...)))", -3, 1.5));
    }

    // Library test MP1: GAUCC / "((.((" -> V2 (A,"((."), V11 (U,"(.("), V29 (C,".((") = 1/5 each; V33 = MFE, V34 = p.
    [Test]
    public void ComputeMiPredFeatures_Binding_HandComputedTriplets()
    {
        var r = AnnotationTools.ComputeMiPredFeatures("GAUCC", "((.((", -12.5, 0.25);
        Assert.That(r.Features, Has.Count.EqualTo(34));
        for (int i = 0; i < 32; i++)
            Assert.That(r.Features[i], Is.EqualTo(i is 1 or 10 or 28 ? 0.2 : 0.0), $"V{i + 1}");
        Assert.That(r.Features[32], Is.EqualTo(-12.5));
        Assert.That(r.Features[33], Is.EqualTo(0.25));
    }

    [Test]
    public void ComputeMiPredFeatures_Binding_MatchesLibrary_DnaEqualsRna()
    {
        string seq = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";
        string db = "((((((((((((((((.(((((.(((((.((((.((...))))))))))).)))))))))))))))))))))";
        var r = AnnotationTools.ComputeMiPredFeatures(seq.Replace('U', 'T'), db, -34.6, 0.001);
        Assert.That(r.Features, Is.EqualTo(MiRnaAnalyzer.ComputeMiPredFeatures(seq, db, -34.6, 0.001)));
    }
}
