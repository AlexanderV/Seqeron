using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CalculateRandfoldPValueTests
{
    private const string Mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";
    private const string Unstructured = "GAUCAGUCAUGACUAGCUAGAUUGACAAGCUUACGAUCAGGACUUAGCAUCGUAGCAUGACUAG";

    [Test]
    public void CalculateRandfoldPValue_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.CalculateRandfoldPValue(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CalculateRandfoldPValue(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateRandfoldPValue(Mir21, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateRandfoldPValue(Mir21, 10001));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.CalculateRandfoldPValue(Mir21, 1, 1, -1));
    }

    [Test]
    public void CalculateRandfoldPValue_Binding_MatchesLibraryLockedValues()
    {
        // Library test CalculateRandfoldPValue_CountsAsRandfold (randfold.c counting; same System.Random(seed) draws).
        var r = AnnotationTools.CalculateRandfoldPValue(Unstructured, 199, 11);
        var m = AnnotationTools.CalculateRandfoldPValue(Mir21, 99, 7);
        Assert.Multiple(() =>
        {
            Assert.That(r.Found, Is.True);
            Assert.That(r.LowerOrEqualCount, Is.EqualTo(20));
            Assert.That(r.PValue, Is.EqualTo(21.0 / 200.0).Within(1e-15));
            Assert.That(r.FreeEnergy, Is.EqualTo(-14.45).Within(1e-9));
            Assert.That(r.Iterations, Is.EqualTo(199));
            Assert.That(r.Seed, Is.EqualTo(11));
            Assert.That(m.PValue, Is.EqualTo(0.01).Within(1e-15));
            Assert.That(m.FreeEnergy, Is.EqualTo(-35.13).Within(1e-9));
        });
    }

    [Test]
    public void CalculateRandfoldPValue_SameSeed_IsDeterministic()
    {
        var a = AnnotationTools.CalculateRandfoldPValue(Unstructured, 30, 5);
        var b = AnnotationTools.CalculateRandfoldPValue(Unstructured, 30, 5);
        Assert.That(b, Is.EqualTo(a));
    }

    [Test]
    public void CalculateRandfoldPValue_ZeroIterations_PIsOne()
    {
        var r = AnnotationTools.CalculateRandfoldPValue(Mir21, 0);
        Assert.That(r.PValue, Is.EqualTo(1.0));
        Assert.That(r.LowerOrEqualCount, Is.EqualTo(0));
    }
}
