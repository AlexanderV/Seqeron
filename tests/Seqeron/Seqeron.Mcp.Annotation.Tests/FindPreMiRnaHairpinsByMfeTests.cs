using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindPreMiRnaHairpinsByMfeTests
{
    // Library MF7: 57 nt designed hairpin (25 bp stems + 7 nt loop); full window dG -48.48, 27 bp.
    private const string ValidHairpin57 =
        "GCAUAGCUAGCUAGCUAGCUAGCUA" + "GAAAUUU" + "UAGCUAGCUAGCUAGCUAGCUAUGC";

    [Test]
    public void FindPreMiRnaHairpinsByMfe_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindPreMiRnaHairpinsByMfe(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindPreMiRnaHairpinsByMfe(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FindPreMiRnaHairpinsByMfe(ValidHairpin57, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FindPreMiRnaHairpinsByMfe(ValidHairpin57, 55, 54));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FindPreMiRnaHairpinsByMfe(ValidHairpin57, 55, 120, double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FindPreMiRnaHairpinsByMfe(ValidHairpin57, 55, 120, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FindPreMiRnaHairpinsByMfe(ValidHairpin57, 55, 120, 0.85, -1));
    }

    [Test]
    public void FindPreMiRnaHairpinsByMfe_Binding_EqualsLibraryAndCarriesMfeFeatures()
    {
        var r = AnnotationTools.FindPreMiRnaHairpinsByMfe(ValidHairpin57);
        var lib = MiRnaAnalyzer.FindPreMiRnaHairpinsByMfe(ValidHairpin57).ToList();
        Assert.That(lib, Is.Not.Empty);
        Assert.That(r.Hairpins, Has.Count.EqualTo(lib.Count));
        for (int i = 0; i < lib.Count; i++)
        {
            Assert.That(r.Hairpins[i].Start, Is.EqualTo(lib[i].Start));
            Assert.That(r.Hairpins[i].End, Is.EqualTo(lib[i].End));
            Assert.That(r.Hairpins[i].DotBracket, Is.EqualTo(lib[i].DotBracket));
            Assert.That(r.Hairpins[i].Mfei, Is.EqualTo(lib[i].Mfei));
            Assert.That(r.Hairpins[i].Amfe, Is.EqualTo(lib[i].Amfe));
            Assert.That(r.Hairpins[i].StemBasePairs, Is.GreaterThanOrEqualTo(16));
        }
        // Full 57-nt window: dG -48.48 (library MF7 / doc), 27 bp.
        var full = r.Hairpins.Single(h => h.Start == 0 && h.End == 56);
        Assert.That(full.FreeEnergy, Is.EqualTo(-48.48).Within(1e-9));
        Assert.That(full.StemBasePairs, Is.EqualTo(27));
    }

    [Test]
    public void FindPreMiRnaHairpinsByMfe_NoHairpin_ReturnsEmpty()
    {
        var r = AnnotationTools.FindPreMiRnaHairpinsByMfe(new string('A', 60));
        Assert.That(r.Hairpins, Is.Empty);
        // too short for the window
        Assert.That(AnnotationTools.FindPreMiRnaHairpinsByMfe("ACGU").Hairpins, Is.Empty);
    }
}
