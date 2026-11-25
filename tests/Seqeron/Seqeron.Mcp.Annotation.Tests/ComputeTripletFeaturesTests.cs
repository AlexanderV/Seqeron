using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ComputeTripletFeaturesTests
{
    [Test]
    public void ComputeTripletFeatures_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeTripletFeatures("", "((((....))))"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeTripletFeatures("GGGGAAAACCCC", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ComputeTripletFeatures(null!, "((((....))))"));
    }

    [Test]
    public void ComputeTripletFeatures_Binding_HandComputedAndLibraryEqual()
    {
        // Library TS6: stem GGGG/CCCC, loop AAAA. G: ".((" = 011 (3), "(((" x2 = 7, "((." = 110 (6); C likewise; index = nt(A0,G1,C2,U3)*8 + bits; 8 nt -> count/8.
        var r = AnnotationTools.ComputeTripletFeatures("GGGGAAAACCCC", "((((....))))");
        Assert.That(r.Found, Is.True);
        var expected = new double[32];
        expected[1 * 8 + 3] = 1 / 8.0; expected[1 * 8 + 7] = 2 / 8.0; expected[1 * 8 + 6] = 1 / 8.0;
        expected[2 * 8 + 3] = 1 / 8.0; expected[2 * 8 + 7] = 2 / 8.0; expected[2 * 8 + 6] = 1 / 8.0;
        Assert.That(r.Features, Is.EqualTo(expected).Within(1e-15));
        Assert.That(r.Features, Is.EqualTo(MiRnaAnalyzer.ComputeTripletFeatures("GGGGAAAACCCC", "((((....))))")));
        Assert.That(r.Features!.Sum(), Is.EqualTo(1.0).Within(1e-12));
    }

    [Test]
    public void ComputeTripletFeatures_NotSingleArmPair_NotFound()
    {
        Assert.That(AnnotationTools.ComputeTripletFeatures("GGGGAAAACCCCGGGGAAAACCCC", "((((....))))((((....))))").Found, Is.False); // two arm pairs (TS6)
        Assert.That(AnnotationTools.ComputeTripletFeatures("GGGGAAAACCCC", "((((...))))").Found, Is.False);   // length mismatch
        Assert.That(AnnotationTools.ComputeTripletFeatures("GGGGNAAACCCC", "((((....))))").Found, Is.False);  // non-ACGU
    }
}
