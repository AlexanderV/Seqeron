using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// MCP wrapper of <c>ProbeDesigner.DesignProbesPrimer3</c> (Primer3 <c>pick_hyb_probe_only</c>).
/// Reference: primer3-py 2.3.1 <c>design_primers</c> (values from PROBE-DESIGN-001, B07 F19).
/// </summary>
[TestFixture]
public class DesignProbesPrimer3Tests
{
    // random.seed(1); 120 nt over ACGT (ProbeDesigner_Primer3Probe_Tests).
    private const string Template =
        "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACC";

    [Test]
    public void DesignProbesPrimer3_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.design_probes_primer3(Template));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(Template, num_return: -1));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(Template, min_size: 0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(Template, min_size: 20, max_size: 19));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(Template, max_size: 37));
    }

    [Test]
    public void DesignProbesPrimer3_Defaults_MatchPrimer3PickHybProbeOnly()
    {
        // primer3-py design_primers({SEQUENCE_TEMPLATE}, {PRIMER_TASK: pick_hyb_probe_only, PRIMER_NUM_RETURN: 5}):
        // PRIMER_INTERNAL_0 = (27, 22), TM 57.57366111818777, PENALTY 4.426338881812228,
        // SELF_ANY_TH 15.431688390643387, SELF_END_TH 1.5231521226383506, HAIRPIN_TH 37.472510801120904.
        var probes = MolToolsTools.design_probes_primer3(Template).Probes;

        Assert.Multiple(() =>
        {
            Assert.That(probes.Select(p => (p.Start, p.Length)),
                Is.EqualTo(new[] { (27, 22), (27, 23), (67, 24), (68, 23), (67, 23) }));
            Assert.That(probes[0].Tm, Is.EqualTo(57.57366111818777).Within(1e-9));
            Assert.That(probes[0].Penalty, Is.EqualTo(4.426338881812228).Within(1e-9));
            Assert.That(probes[0].SelfAnyTh, Is.EqualTo(15.431688390643387).Within(1e-9));
            Assert.That(probes[0].SelfEndTh, Is.EqualTo(1.5231521226383506).Within(1e-9));
            Assert.That(probes[0].HairpinTh, Is.EqualTo(37.472510801120904).Within(1e-9));
        });
    }

    [Test]
    public void DesignProbesPrimer3_CustomSettings_MatchPrimer3()
    {
        // PRIMER_INTERNAL_OPT_TM 68, MIN_TM 65, MAX_TM 72, MAX_SIZE 30, OPT_SIZE 24, SALT_DIVALENT 3.0,
        // DNTP_CONC 0.8, DNA_CONC 250, PRIMER_NUM_RETURN 3 (primer3-py 2.3.1).
        var probes = MolToolsTools.design_probes_primer3(Template.ToLowerInvariant(), num_return: 3,
            opt_size: 24, max_size: 30, min_tm: 65, opt_tm: 68, max_tm: 72,
            divalent_mm: 3.0, dntp_mm: 0.8, dna_conc_nm: 250).Probes;

        Assert.Multiple(() =>
        {
            Assert.That(probes.Select(p => p.Sequence), Is.EqualTo(new[]
            {
                "CCCACCTGGTGATCCTATGCTTGT", "GTCCCACCTGGTGATCCTATGCTT", "TCCCACCTGGTGATCCTATGCTTG",
            }));
            Assert.That(probes[0].Tm, Is.EqualTo(68.03185451461161).Within(1e-9));
            Assert.That(probes[0].Penalty, Is.EqualTo(0.031854514611609375).Within(1e-9));
        });
    }
}
