using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class IsConservedByBranchLengthTests
{
    [Test]
    public void IsConservedByBranchLength_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.IsConservedByBranchLength(1.0, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsConservedByBranchLength(1.0, "Bogus"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsConservedByBranchLength(1.0, null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsConservedByBranchLength(1.0, "3"), "numeric enum");
        Assert.Throws<ArgumentException>(() => AnnotationTools.IsConservedByBranchLength(1.0, "Seed8mer", "Bogus"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.IsConservedByBranchLength(-0.1, "Seed8mer"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.IsConservedByBranchLength(double.PositiveInfinity, "Seed8mer"));
    }

    // TargetScanHuman7: 8mer 1.8, 7mer-m8 2.8, 7mer-1a 3.6; Mouse7: 0.6, 1.8, 2.5 (targetscan_70_BL_PCT.pl).
    [TestCase("Seed8mer", 1.8, "TargetScanHuman7", true, 1.8)]
    [TestCase("Seed8mer", 1.7999, "TargetScanHuman7", false, 1.8)]
    [TestCase("Seed7merM8", 2.8, "TargetScanHuman7", true, 2.8)]
    [TestCase("Seed7merM8", 2.7, "TargetScanHuman7", false, 2.8)]
    [TestCase("Seed7merA1", 3.6, "TargetScanHuman7", true, 3.6)]
    [TestCase("Seed7merA1", 3.5, "TargetScanHuman7", false, 3.6)]
    [TestCase("Seed8mer", 0.6, "TargetScanMouse7", true, 0.6)]
    [TestCase("Seed7merM8", 1.8, "TargetScanMouse7", true, 1.8)]
    [TestCase("Seed7merA1", 2.4, "TargetScanMouse7", false, 2.5)]
    public void IsConservedByBranchLength_Binding_Thresholds(string type, double bls, string th, bool expected, double threshold)
    {
        var r = AnnotationTools.IsConservedByBranchLength(bls, type, th);
        Assert.That(r.Conserved, Is.EqualTo(expected));
        Assert.That(r.Threshold, Is.EqualTo(threshold));
        var lib = MiRnaAnalyzer.IsConservedByBranchLength(bls, Enum.Parse<MiRnaAnalyzer.TargetSiteType>(type),
            th == "TargetScanMouse7" ? MiRnaAnalyzer.BranchLengthConservationThresholds.TargetScanMouse7 : MiRnaAnalyzer.BranchLengthConservationThresholds.TargetScanHuman7);
        Assert.That(r.Conserved, Is.EqualTo(lib));
    }

    [TestCase("Seed6mer")]
    [TestCase("Offset6mer")]
    public void IsConservedByBranchLength_6mer_NeverConserved(string type)
    {
        var r = AnnotationTools.IsConservedByBranchLength(100.0, type);
        Assert.That(r.Conserved, Is.False);
        Assert.That(r.Threshold, Is.Null);
    }
}
