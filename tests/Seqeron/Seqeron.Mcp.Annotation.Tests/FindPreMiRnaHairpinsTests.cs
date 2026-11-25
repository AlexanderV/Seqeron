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
        // Default = MFE fold per window (F40). The 57-nt hairpin (27 bp, 3-nt loop) yields exactly 6 accepted
        // windows (starts 0..2, ends 54..56); the full window folds to ΔG° −48.48 (ViennaRNA 2.8 −48.40).
        var result = AnnotationTools.FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55);

        Assert.That(result.Hairpins, Has.Count.EqualTo(6));
        var full = result.Hairpins.Single(h => h.Start == 0 && h.End == 56);
        Assert.Multiple(() =>
        {
            Assert.That(full.Sequence, Is.EqualTo(ValidHairpin57));
            Assert.That(full.Structure, Is.EqualTo(new string('(', 27) + "..." + new string(')', 27)));
            Assert.That(full.FreeEnergy, Is.EqualTo(-48.48).Within(0.01));
            Assert.That(full.MatureSequence, Is.EqualTo("GCAUAGCUAGCUAGCUAGCUAG"));
            Assert.That(full.StarSequence, Is.EqualTo("AGCUAGCUAGCUAGCUAUGC"), "find_star: 2-nt 3' overhang, clipped at the window end");
        });
    }

    [Test]
    public void FindPreMiRnaHairpins_BulgedRealPreMiRna_DetectedByDefault()
    {
        // hsa-mir-21 (MI0000077): bulged stem, missed by the consecutive-pairing screen; ViennaRNA 2.8 −34.70.
        const string mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";
        var h = AnnotationTools.FindPreMiRnaHairpins(mir21).Hairpins.Single(x => x.Start == 1 && x.End == 71);
        Assert.That(h.FreeEnergy, Is.EqualTo(-34.70).Within(0.5));
    }

    [Test]
    public void FindPreMiRnaHairpins_ShortSequence_ReturnsEmpty()
    {
        // 47 nt is below the 55 nt minimum hairpin length -> no candidates.
        var result = AnnotationTools.FindPreMiRnaHairpins(ShortHairpin47, minHairpinLength: 55);
        Assert.That(result.Hairpins, Is.Empty);
    }

    [Test]
    public void FindPreMiRnaHairpins_ConsecutiveScreen_Option_MatchesLibraryScreen()
    {
        // useConsecutiveScreen=true delegates to the fast screen (no MFE fold); parity with the library.
        var tool = AnnotationTools.FindPreMiRnaHairpins(ValidHairpin57, 55, 120, 22, useConsecutiveScreen: true).Hairpins;
        var lib = Seqeron.Genomics.Annotation.MiRnaAnalyzer
            .FindPreMiRnaHairpins(ValidHairpin57, 55, 120, 22, useConsecutiveScreen: true).ToList();
        Assert.That(tool, Has.Count.EqualTo(lib.Count));
        Assert.That(tool.Select(h => (h.Start, h.End)), Is.EqualTo(lib.Select(h => (h.Start, h.End))));
        // blunt labels: mature = first 22 nt of the window
        Assert.That(tool[0].MatureSequence, Is.EqualTo(tool[0].Sequence[..22]));
    }

    [Test]
    public void FindPreMiRnaHairpins_MinMfei_FiltersAndIsValidated()
    {
        Assert.That(AnnotationTools.FindPreMiRnaHairpins(ValidHairpin57, minMfei: 1000).Hairpins, Is.Empty);
        Assert.That(AnnotationTools.FindPreMiRnaHairpins(ValidHairpin57, minMfei: 0).Hairpins, Has.Count.EqualTo(6));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FindPreMiRnaHairpins(ValidHairpin57, minMfei: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FindPreMiRnaHairpins(ValidHairpin57, minMfei: double.NaN));
    }
}
