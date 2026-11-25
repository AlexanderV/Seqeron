using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ScorePreMiRnaFeaturesTests
{
    private const string Mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";

    [Test]
    public void ScorePreMiRnaFeatures_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScorePreMiRnaFeatures(double.NaN, 48, 1, 0.5, 0.9));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScorePreMiRnaFeatures(-35, -1, 1, 0.5, 0.9));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScorePreMiRnaFeatures(-35, 48, double.PositiveInfinity, 0.5, 0.9));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScorePreMiRnaFeatures(-35, 48, 1, 1.1, 0.9));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScorePreMiRnaFeatures(-35, 48, 1, 0.5, -0.1));
    }

    [Test]
    public void ScorePreMiRnaFeatures_EqualsClassifyProbability_ExtractFoldedIntoClassify()
    {
        // classify_pre_mirna returns the library's extracted features; scoring them reproduces its probability,
        // so a separate extract tool would be redundant. hsa-mir-21 (library CL1): dG -35.13, AMFE 48.79166666666667, MFEI 1.0037142857142858.
        var c = AnnotationTools.ClassifyPreMiRna(Mir21);
        var f = c.Features!;
        Assert.That(f.FreeEnergy, Is.EqualTo(-35.13).Within(1e-9));
        Assert.That(f.Amfe, Is.EqualTo(48.79166666666667).Within(1e-12));
        Assert.That(f.Mfei, Is.EqualTo(1.0037142857142858).Within(1e-12));
        var lib = MiRnaAnalyzer.ExtractPreMiRnaFeatures(Mir21)!.Value;
        Assert.That(f.Mfei, Is.EqualTo(lib.Mfei));
        Assert.That(f.PairedFraction, Is.EqualTo(lib.PairedFraction));

        var s = AnnotationTools.ScorePreMiRnaFeatures(f.FreeEnergy, f.Amfe, f.Mfei, f.GcContent, f.PairedFraction);
        Assert.That(s.NaturalProbability, Is.EqualTo(c.NaturalProbability));
        Assert.That(s.NaturalProbability, Is.EqualTo(MiRnaAnalyzer.ScorePreMiRnaFeatures(lib)));
        Assert.That(s.NaturalProbability, Is.InRange(0.95, 1.0));
    }

    [Test]
    public void ScorePreMiRnaFeatures_BackgroundLikeFeatures_ScoreLow()
    {
        // weak, unpaired, low-MFEI feature vector -> background
        var s = AnnotationTools.ScorePreMiRnaFeatures(-5.0, 7.0, 0.15, 0.45, 0.2);
        Assert.That(s.NaturalProbability, Is.LessThan(0.5));
    }
}
