using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PctFromBranchLengthTests
{
    // let-7 (GAGGUAG) 7mer-m8 row of 7mer-m8_PCT_parameters.txt; kit: Bls 1.6379 -> PCT 0.0058.
    private const double B0 = -0.221367319203409, B1 = 1.24081845475817, B2 = 1.17370692612852, B3 = 3.41806979609985;

    [Test]
    public void PctFromBranchLength_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PctFromBranchLength(1.0, "Seed7merM8", "GAGGUAG"));
        Assert.DoesNotThrow(() => AnnotationTools.PctFromBranchLength(1.0, null, null, B0, B1, B2, B3));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PctFromBranchLength(1.0), "no mode");
        Assert.Throws<ArgumentException>(() => AnnotationTools.PctFromBranchLength(1.0, "Seed7merM8", "GAGGUAG", B0, B1, B2, B3), "both modes");
        Assert.Throws<ArgumentException>(() => AnnotationTools.PctFromBranchLength(1.0, null, null, B0, B1), "partial sigmoid");
        Assert.Throws<ArgumentException>(() => AnnotationTools.PctFromBranchLength(1.0, null, "GAGGUAG"), "lookup needs siteType");
        Assert.Throws<ArgumentException>(() => AnnotationTools.PctFromBranchLength(1.0, "Seed6mer", "GAGGUAG"), "no 6mer row");
        Assert.Throws<ArgumentException>(() => AnnotationTools.PctFromBranchLength(1.0, "Bogus", "GAGGUAG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PctFromBranchLength(1.0, "7", "GAGGUAG"), "numeric enum");
        Assert.Throws<ArgumentException>(() => AnnotationTools.PctFromBranchLength(1.0, "Seed8mer", "GAGGUA"), "6 nt seed");
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.PctFromBranchLength(-1.0, "Seed8mer", "GAGGUAG"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.PctFromBranchLength(double.NaN, null, null, B0, B1, B2, B3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.PctFromBranchLength(-1.0, "Seed8mer", "ACACACU"), "Bls validated even for NA");
    }

    [Test]
    public void PctFromBranchLength_Binding_KitRow_Bls1_6379_Pct0_0058()
    {
        var byLookup = AnnotationTools.PctFromBranchLength(1.6379, "Seed7merM8", "GAGGUAG");
        var bySigmoid = AnnotationTools.PctFromBranchLength(1.6379, "Seed7merM8", null, B0, B1, B2, B3);
        Assert.Multiple(() =>
        {
            Assert.That(byLookup.Pct, Is.EqualTo(0.0058));
            Assert.That(byLookup.IsNa, Is.False);
            Assert.That(byLookup.Parameters, Is.EqualTo(new PctSigmoidParametersDto(B0, B1, B2, B3)));
            Assert.That(bySigmoid.Pct, Is.EqualTo(0.0058));
            Assert.That(bySigmoid.Pct, Is.EqualTo(MiRnaAnalyzer.PctFromBranchLength(1.6379, new MiRnaAnalyzer.PctSigmoidParameters(B0, B1, B2, B3), MiRnaAnalyzer.TargetSiteType.Seed7merM8)));
        });
    }

    [Test]
    public void PctFromBranchLength_Edges()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.PctFromBranchLength(0.0, "Seed7merM8", "GAGGUAG").Pct, Is.EqualTo(0.0), "Bls 0");
            Assert.That(AnnotationTools.PctFromBranchLength(1.6379, null, null, 0, 1, 1, 0).Pct, Is.EqualTo(0.8372), "1/(1+e^-1.6379)");
            Assert.That(AnnotationTools.PctFromBranchLength(1.6379, "Seed6mer", null, 0, 1, 1, 0).Pct, Is.EqualTo(0.0), "6mer with sigmoid -> 0");
            var na = AnnotationTools.PctFromBranchLength(1.6379, "Seed8mer", "ACACACU");
            Assert.That(na.Pct, Is.Null);
            Assert.That(na.IsNa, Is.True);
            Assert.That(na.Parameters, Is.Null);
            Assert.That(AnnotationTools.PctFromBranchLength(1.6379, "seed7merm8", "gaggtag").Pct, Is.EqualTo(0.0058), "case / T as U");
        });
    }
}
