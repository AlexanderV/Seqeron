using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindPreMiRnaHairpinsTests
{
    // Mirrors Seqeron.Genomics.Tests MiRnaAnalyzer_PreMiRna_Tests: 57 nt hairpin =
    // 25 nt 5' stem + 7 nt loop + 25 nt 3' stem (reverse complement).
    private const string ValidHairpin57 =
        "GCAUAGCUAGCUAGCUAGCUAGCUA" + "GAAAUUU" + "UAGCUAGCUAGCUAGCUAGCUAUGC";

    // 47 nt hairpin (20 nt stems) — below the default/used min length of 55.
    private const string ShortHairpin47 =
        "GCAUAGCUAGCUAGCUAGCU" + "GAAAUUU" + "AGCUAGCUAGCUAGCUAUGC";

    [Test]
    public void FindPreMiRnaHairpins_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindPreMiRnaHairpins(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindPreMiRnaHairpins(null!));
    }

    [Test]
    public void FindPreMiRnaHairpins_Binding_InvokesSuccessfully()
    {
        // 57 nt input with minHairpinLength=55 yields exactly 2 candidates:
        // the full 57-nt window (Start=0, End=56, stem=23) and a 55-nt sub-window at offset 1.
        var result = AnnotationTools.FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55);

        Assert.That(result.Hairpins, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(result.Hairpins[0].Start, Is.EqualTo(0));
            Assert.That(result.Hairpins[0].End, Is.EqualTo(56));
            Assert.That(result.Hairpins[0].Sequence, Is.EqualTo(ValidHairpin57));
        });
    }

    [Test]
    public void FindPreMiRnaHairpins_ShortSequence_ReturnsEmpty()
    {
        // 47 nt is below the 55 nt minimum hairpin length -> no candidates.
        var result = AnnotationTools.FindPreMiRnaHairpins(ShortHairpin47, minHairpinLength: 55);
        Assert.That(result.Hairpins, Is.Empty);
    }
}
