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
        // Primer3 _pr_data_control: "PRIMER_NUM_RETURN < 1".
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(Template, num_return: 0));
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

    [Test]
    public void DesignProbesPrimer3_AlignmentMode_MatchesPrimer3AndSerializes()
    {
        // primer3-py 2.3.1 pick_hyb_probe_only, PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_INTERNAL_MAX_SELF_ANY 6:
        // PRIMER_INTERNAL_0 (8, 23), PENALTY 4.850867227162098, SELF_ANY 6, SELF_END 2.
        const string t = "ATTCTCAGAGGCTCGTACAAACGTATGCCCTAGCTTTTACCACTTAACGCCGTCAAAATGTGCCTATTTTGGAACGAAGGATTCTGTCCTTCGTTCCTTCTTAGTAT";
        var r = MolToolsTools.design_probes_primer3(t, thermodynamic_oligo_alignment: false, max_self_any: 6);
        var p0 = r.Probes[0];

        Assert.Multiple(() =>
        {
            Assert.That((p0.Start, p0.Length), Is.EqualTo((8, 23)));
            Assert.That(p0.Penalty, Is.EqualTo(4.850867227162098).Within(1e-9));
            Assert.That(p0.SelfAny, Is.EqualTo(6.0));
            Assert.That(p0.SelfEnd, Is.EqualTo(2.0));
            Assert.That(p0.HairpinTh, Is.NaN);
            Assert.That(System.Text.Json.JsonSerializer.Serialize(r), Does.Contain("\"NaN\""), "NaN Th fields serialize");
        });
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.design_probes_primer3(t, max_self_end: -1));
    }

    [Test]
    public void DesignProbesPrimer3_MishybLibrary_MatchesPrimer3()
    {
        // primer3-py 2.3.1 design_primers({SEQUENCE_TEMPLATE}, {PRIMER_TASK: pick_hyb_probe_only}, mishyb_lib = lib):
        // PRIMER_INTERNAL_k = (27,22) (27,23) (26,23) (25,24) (24,25), PRIMER_INTERNAL_k_LIBRARY_MISPRIMING (8.0, site*2);
        // with PRIMER_INTERNAL_WT_LIBRARY_MISHYB 0.5 and PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 30: PRIMER_INTERNAL_0_PENALTY
        // 8.426338881812228.
        var lib = new Dictionary<string, string>
        {
            ["site*2"] = "CCCACCTGGTGATCCTATGCTTGTG",
            ["mut"] = "CCCACCTGGTGTTCCTATGC",
            ["amb"] = "GTCCCRCCTGGNGATCCTAYG",
            ["tiny"] = "CA",
        };
        var probes = MolToolsTools.design_probes_primer3(Template, mishyb_library: lib).Probes;
        var weighted = MolToolsTools.design_probes_primer3(Template, mishyb_library: lib, wt_library_mishyb: 0.5,
            max_library_mishyb: 30).Probes;

        Assert.Multiple(() =>
        {
            Assert.That(probes.Select(p => (p.Start, p.Length)),
                Is.EqualTo(new[] { (27, 22), (27, 23), (26, 23), (25, 24), (24, 25) }));
            Assert.That(probes.Select(p => p.LibraryMishyb), Is.All.EqualTo(8.0));
            Assert.That(probes.Select(p => p.LibraryMishybName), Is.All.EqualTo("site*2"));
            Assert.That(weighted[0].Penalty, Is.EqualTo(8.426338881812228).Within(1e-9));
        });
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(Template, wt_library_mishyb: 1));
    }

    [Test]
    public void DesignProbesPrimer3_SequenceQuality_MatchesPrimer3()
    {
        // primer3-py 2.3.1 pick_hyb_probe_only on random.Random(11) 120 nt with SEQUENCE_QUALITY q[i] = 12 if i % 37 == 4,
        // 26 if i % 13 == 2, else 34 + i % 7; PRIMER_INTERNAL_MIN_QUALITY 20, PRIMER_INTERNAL_WT_SEQ_QUAL 0.3 →
        // PRIMER_INTERNAL_0 [85,25], PENALTY 30.044504099519497, MIN_SEQ_QUALITY 26.
        const string p = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC";
        int[] q = Enumerable.Range(0, p.Length).Select(i => i % 37 == 4 ? 12 : i % 13 == 2 ? 26 : 34 + i % 7).ToArray();
        var r = MolToolsTools.design_probes_primer3(p, sequence_quality: q, min_quality: 20, wt_seq_qual: 0.3, wt_end_qual: 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.Probes[0].Start, Is.EqualTo(85));
            Assert.That(r.Probes[0].Length, Is.EqualTo(25));
            Assert.That(r.Probes[0].Penalty, Is.EqualTo(30.044504099519497).Within(1e-9));
            Assert.That(r.Probes[0].MinSequenceQuality, Is.EqualTo(26));
            Assert.Throws<ArgumentException>(() => MolToolsTools.design_probes_primer3(p, min_quality: 20));
        });
    }
}
