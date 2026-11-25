using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ClassifyPreMiRnaTests
{
    private const string Mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";

    [Test]
    public void ClassifyPreMiRna_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyPreMiRna(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyPreMiRna(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRna(Mir21, 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRna(Mir21, double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRna(Mir21, 0.5, -1));
    }

    [Test]
    public void ClassifyPreMiRna_Binding_MatchesLibrary_HsaMir21IsNatural()
    {
        var r = AnnotationTools.ClassifyPreMiRna(Mir21);
        var lib = MiRnaAnalyzer.ClassifyPreMiRna(Mir21)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(r.Found, Is.True);
            Assert.That(r.IsNatural, Is.True);
            Assert.That(r.NaturalProbability, Is.GreaterThan(0.95));   // library test ClassifyPreMiRna_HsaMir21_IsNatural
            Assert.That(r.NaturalProbability, Is.EqualTo(lib.NaturalProbability).Within(1e-15));
            Assert.That(r.Features!.FreeEnergy, Is.EqualTo(lib.Features.FreeEnergy).Within(1e-12));
            Assert.That(r.Features.Length, Is.EqualTo(Mir21.Length));
            Assert.That(r.Features.StemBasePairs, Is.EqualTo(lib.Features.StemBasePairs));
        });
    }

    [Test]
    public void ClassifyPreMiRna_Threshold_IsHonoured()
    {
        Assert.That(AnnotationTools.ClassifyPreMiRna(Mir21, 1.0).IsNatural, Is.False); // P < 1 (library ..._Threshold_IsHonoured)
    }

    [Test]
    public void ClassifyPreMiRna_TooShortToFold_StillReportsOrNotFound_NoThrow()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ClassifyPreMiRna("ACGU"));
    }
}
