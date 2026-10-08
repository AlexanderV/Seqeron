using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// B07 audit round 3, A3-28 (F56): the probe tools expose the F55 THAL_MAX_ALIGN opt-in
/// (thermodynamic_screen_max_length / max_align_length; default 60 = Primer3). Reference values: thal.c + thal_parameters.c
/// of primer3-py 2.3.1 compiled with -DTHAL_MAX_ALIGN=10000 (the 60-nt build and primer3-py raise "At least one sequence
/// must be equal to or shorter than 60bp"), 50 mM / 0 / 0 / 50 nM, 37 °C, max loop 30.
/// </summary>
[TestFixture]
public class ProbeThalMaxAlignTests
{
    private const double TmTol = 1e-9;

    // 12-bp stem + TTTT loop with random flanks (82 nt): hairpin 68.102992186803021 °C (F55).
    private const string StemLoop82 =
        "CTTAGGTGGATAGGGAGTGAGCAACAAACGACAGAAGGTATGTTTTCATACCTTCTGTGATCGTTTCTCCCATGCCAAGTTG";
    private const string Random120 =
        "GTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAATCTACCCCCTGTTATGCGCG";
    private static readonly string Probe75 = Random120.Substring(0, 75);
    // Probe75 with substitutions at 10, 35, 60 inside flanks: aligned site [14, 88], 75 nt.
    private const string NonTarget =
        "GATTACAGATTACAGTTTCGGAACATGCGTTTTAGGTATGTCTTAGTGAATCTAAATACCAAGGCAGTCCTCGAACCGTTCCTAATAAGCCTTGGAACC";
    private static readonly ProbeDesigner.ProbeParameters Window82 = new(82, 82, -1000, 1000, 0, 1, 100, true, 1.0);
    private const string HairpinWarning = "Hairpin Tm 68.1°C exceeds 47°C (ntthal)";

    [Test]
    public void ValidateProbe_ThermodynamicScreenMaxLength_ScreensLongProbeAndDuplex()
    {
        var byDefault = MolToolsTools.validate_probe(Probe75, Array.Empty<string>(), non_target_sequences: new[] { NonTarget });
        var optIn = MolToolsTools.validate_probe(Probe75, Array.Empty<string>(), non_target_sequences: new[] { NonTarget },
            max_duplex_tm: 60, thermodynamic_screen_max_length: 120);
        Assert.Multiple(() =>
        {
            Assert.That(byDefault.ThermodynamicScreen, Is.False);
            Assert.That(byDefault.CrossHybridization[0].DuplexTm, Is.Null, "probe and site both > 60 nt");
            Assert.That(optIn.ThermodynamicScreen, Is.True);
            Assert.That(optIn.HairpinTm, Is.EqualTo(36.318725375690178).Within(TmTol));
            Assert.That(optIn.SelfDimerTm, Is.EqualTo(23.325175622075108).Within(TmTol));
            Assert.That(optIn.SelfEndDimerTm, Is.EqualTo(17.579447109879538).Within(TmTol));
            Assert.That(optIn.CrossHybridization[0].DuplexTm, Is.EqualTo(66.457703046655695).Within(TmTol));
            Assert.That(optIn.CrossHybridization[0].ExceedsDuplexTmThreshold, Is.True);
            Assert.That(() => MolToolsTools.validate_probe(Probe75, Array.Empty<string>(), thermodynamic_screen_max_length: 59),
                Throws.InstanceOf<ArgumentException>());
        });
    }

    [Test]
    public void DesignProbes_ThermodynamicScreenMaxLength_ScreensLongProbe()
    {
        var byDefault = MolToolsTools.design_probes(StemLoop82, Window82).Probes.Single();
        var optIn = MolToolsTools.design_probes(StemLoop82, Window82, thermodynamic_screen_max_length: 82).Probes.Single();
        var viaParameters = MolToolsTools.design_probes(StemLoop82, Window82 with { ThermodynamicScreenMaxLength = 82 }).Probes.Single();
        Assert.Multiple(() =>
        {
            Assert.That(byDefault.Warnings, Has.None.Contain("(ntthal)"));
            Assert.That(optIn.Warnings, Has.Some.EqualTo(HairpinWarning));
            Assert.That(viaParameters.Warnings, Is.EqualTo(optIn.Warnings), "null keeps the parameters' value");
            Assert.That(() => MolToolsTools.design_probes(StemLoop82, Window82, thermodynamic_screen_max_length: 10001),
                Throws.InstanceOf<ArgumentException>());
        });
    }

    [Test]
    public void DesignTilingAndAntisenseProbes_ThermodynamicScreenMaxLength_ScreensLongProbe()
    {
        var tilingDefault = MolToolsTools.design_tiling_probes(StemLoop82, 82, 20, Window82).Probes.Single();
        var tiling = MolToolsTools.design_tiling_probes(StemLoop82, 82, 20, Window82, thermodynamic_screen_max_length: 82).Probes.Single();
        string mrna = DnaSequence.GetReverseComplementString(StemLoop82);
        var antisense = MolToolsTools.design_antisense_probes(mrna, Window82, thermodynamic_screen_max_length: 82).Probes.Single();
        Assert.Multiple(() =>
        {
            Assert.That(tilingDefault.Warnings, Has.None.Contain("(ntthal)"));
            Assert.That(tiling.Warnings, Has.Some.EqualTo(HairpinWarning));
            Assert.That(antisense.Sequence, Is.EqualTo(StemLoop82));
            Assert.That(antisense.Warnings, Has.Some.EqualTo(HairpinWarning));
        });
    }

    [Test]
    public void DesignMolecularBeacon_MaxAlignLength_ReportsStemLoopTmOfLongBeacon()
    {
        // 74-nt beacon GGGCCCC + Random120[0..60) + GGGGCCC: hairpin 45.449128344157543 °C.
        var byDefault = MolToolsTools.design_molecular_beacon(Random120, 60, 7).Probe!.Value;
        var optIn = MolToolsTools.design_molecular_beacon(Random120, 60, 7, max_align_length: 100).Probe!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(optIn.Sequence.Length, Is.EqualTo(74));
            Assert.That(byDefault.Warnings, Has.None.Contain("(ntthal)"));
            Assert.That(optIn.Warnings, Has.Some.EqualTo("Stem-loop (hairpin) Tm 45.4°C (ntthal)"));
            Assert.That(() => MolToolsTools.design_molecular_beacon(Random120, 60, 7, max_align_length: 59),
                Throws.InstanceOf<ArgumentException>());
        });
    }
}
