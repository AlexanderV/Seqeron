// PROBE-DESIGN-001 — Primer3 hybridization-probe picker, thermodynamic probe screens, probe Tm,
// molecular weight and nearest-neighbour extinction coefficient.
// Reference: primer3-py 2.3.1 design_primers(PRIMER_TASK = pick_hyb_probe_only), calc_tm;
//            Biopython 1.88 Bio.SeqUtils.molecular_weight; Cantor, Warshaw & Shapiro (1970) ε260 table.
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class ProbeDesigner_Primer3Probe_Tests
{
    // random.seed(1); 120 nt over ACGT (see B07.md F17).
    private const string Template =
        "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACC";

    [Test]
    public void DesignProbesPrimer3_Defaults_MatchPrimer3PickHybProbeOnly()
    {
        // primer3-py design_primers({SEQUENCE_TEMPLATE}, {PRIMER_TASK: pick_hyb_probe_only, PRIMER_NUM_RETURN: 5})
        var expected = new (int Start, int Len, double Tm, double Penalty, double Any, double End, double Hp)[]
        {
            (27, 22, 57.57366111818777, 4.426338881812228, 15.431688390643387, 1.5231521226383506, 37.472510801120904),
            (27, 23, 57.926326193810496, 5.073673806189504, 15.431688390643387, 0.0, 37.472510801120904),
            (67, 24, 58.65366346952396, 5.346336530476037, 24.172680221686505, 6.868622949063649, 30.923663896648236),
            (68, 23, 57.61036380090496, 5.389636199095037, 24.172680221686505, 0.0, 30.923663896648236),
            (67, 23, 57.61036380090496, 5.389636199095037, 24.172680221686505, 5.339097117014546, 30.923663896648236),
        };

        var probes = ProbeDesigner.DesignProbesPrimer3(Template);

        Assert.That(probes, Has.Count.EqualTo(5));
        Assert.Multiple(() =>
        {
            for (int i = 0; i < expected.Length; i++)
            {
                var (p, e) = (probes[i], expected[i]);
                Assert.That((p.Start, p.Length), Is.EqualTo((e.Start, e.Len)), $"rank {i}");
                Assert.That(p.Sequence, Is.EqualTo(Template.Substring(e.Start, e.Len)));
                Assert.That(p.Tm, Is.EqualTo(e.Tm).Within(1e-9));
                Assert.That(p.Penalty, Is.EqualTo(e.Penalty).Within(1e-9));
                Assert.That(p.SelfAnyTh, Is.EqualTo(e.Any).Within(1e-9));
                Assert.That(p.SelfEndTh, Is.EqualTo(e.End).Within(1e-9));
                Assert.That(p.HairpinTh, Is.EqualTo(e.Hp).Within(1e-9));
            }
        });
    }

    [Test]
    public void DesignProbesPrimer3_TaqManLikeSettingsAndPcrBuffer_MatchPrimer3()
    {
        // PRIMER_INTERNAL_OPT_TM 68, MIN_TM 65, MAX_TM 72, MAX_SIZE 30, OPT_SIZE 24,
        // SALT_DIVALENT 3.0, DNTP_CONC 0.8, DNA_CONC 250, PRIMER_NUM_RETURN 3.
        var s = new ProbeDesigner.Primer3ProbeSettings(OptSize: 24, MaxSize: 30, MinTm: 65, OptTm: 68, MaxTm: 72,
            DivalentMillimolar: 3.0, DntpMillimolar: 0.8, DnaConcentrationNanomolar: 250);
        var probes = ProbeDesigner.DesignProbesPrimer3(Template.ToLowerInvariant(), s, 3);

        Assert.Multiple(() =>
        {
            Assert.That(probes.Select(p => p.Sequence), Is.EqualTo(new[]
            {
                "CCCACCTGGTGATCCTATGCTTGT", "GTCCCACCTGGTGATCCTATGCTT", "TCCCACCTGGTGATCCTATGCTTG",
            }));
            Assert.That(probes.Select(p => p.Start), Is.EqualTo(new[] { 70, 68, 69 }));
            Assert.That(probes[0].Tm, Is.EqualTo(68.03185451461161).Within(1e-9));
            Assert.That(probes[0].Penalty, Is.EqualTo(0.031854514611609375).Within(1e-9));
            Assert.That(probes[2].SelfEndTh, Is.EqualTo(11.082413688123154).Within(1e-9));
            Assert.That(probes[0].HairpinTh, Is.EqualTo(36.45687764130065).Within(1e-9));
        });
    }

    [Test]
    public void DesignProbesPrimer3_ThermodynamicLimits_RejectSelfComplementaryWindows()
    {
        // primer3-py explain for this template with MIN_GC 0, MAX_POLY_X 20, MIN_TM 0, MAX_TM 100:
        // "considered 99, high any compl 15, high hairpin stability 26, ok 58" (a self-any failure stops the
        // 5' extension of that 3' end); with the ntthal limits lifted: "considered 225, ok 225".
        const string t = "AAAAAAAAAAGCGCGCGCAAAATTTTGCGCGCGCAAAAAAAAAA";
        var s = new ProbeDesigner.Primer3ProbeSettings(MinGcPercent: 0, MaxPolyX: 20, MinTm: 0, MaxTm: 100);
        var all = ProbeDesigner.DesignProbesPrimer3(t, s, 1000);
        var unscreened = ProbeDesigner.DesignProbesPrimer3(t,
            s with { MaxSelfAnyTh = 1000, MaxSelfEndTh = 1000, MaxHairpinTh = 1000 }, 1000);
        Assert.Multiple(() =>
        {
            Assert.That(all, Has.Count.EqualTo(58));
            Assert.That(unscreened, Has.Count.EqualTo(225));
            Assert.That(all.All(p => p.SelfAnyTh <= 47 && p.SelfEndTh <= 47 && p.HairpinTh <= 47), Is.True);
        });
    }

    [Test]
    public void DesignProbesPrimer3_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => ProbeDesigner.DesignProbesPrimer3(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            ProbeDesigner.DesignProbesPrimer3(Template, new ProbeDesigner.Primer3ProbeSettings(MaxSize: 37)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbesPrimer3(Template, null, -1));
        Assert.That(ProbeDesigner.DesignProbesPrimer3("ACGTNNNNNNNNNNNNNNNNNNNNACGT"), Is.Empty,
            "PRIMER_INTERNAL_MAX_NS_ACCEPTED = 0");
    }

    [Test]
    public void DesignProbes_ProbeTm_IsPrimer3SeqtmAtParameterConditions()
    {
        // qPCR preset, Primer3 probe defaults (50 nM, 50 mM, no Mg/dNTP) and a PCR buffer (1.5 mM Mg, 0.6 mM dNTP).
        var p = ProbeDesigner.Defaults.qPCR with { MinTm = 0, MaxTm = 100 };
        var probe = ProbeDesigner.DesignProbes(Template, p, 1).Single();
        var pcr = ProbeDesigner.DesignProbes(Template, p with { DivalentMillimolar = 1.5, DntpMillimolar = 0.6 }, 50)
            .First(x => x.Start == probe.Start && x.Sequence.Length == probe.Sequence.Length);
        Assert.Multiple(() =>
        {
            Assert.That(probe.Tm, Is.EqualTo(
                PrimerDesigner.CalculateMeltingTemperaturePrimer3(probe.Sequence, 50, 50, 0, 0)).Within(1e-12));
            Assert.That(pcr.Tm, Is.EqualTo(
                PrimerDesigner.CalculateMeltingTemperaturePrimer3(probe.Sequence, 50, 50, 1.5, 0.6)).Within(1e-12));
        });
    }

    [Test]
    public void DesignProbes_ThermodynamicScreen_FlagsSelfDimerByNtthalTm()
    {
        // GGGGCCCC-rich palindrome: ntthal self-dimer Tm far above 47 °C; random-like probe below.
        const string palindrome = "ACGCGCGCGCGCGCGCGCGT"; // self-complementary 20-mer
        var st = PrimerDesigner.CalculatePrimer3OligoStructure(palindrome, 50, 0, 0, 50)!.Value;
        var p = new ProbeDesigner.ProbeParameters(20, 20, -1000, 1000, 0, 1, 100, true, 1.0);
        var probe = ProbeDesigner.DesignProbes(palindrome, p).Single();
        var heuristic = ProbeDesigner.DesignProbes(palindrome, p with { StructureScreen = ProbeDesigner.ProbeStructureScreen.Heuristic }).Single();
        Assert.Multiple(() =>
        {
            Assert.That(st.SelfAnyTh, Is.GreaterThan(47));
            Assert.That(probe.Warnings, Has.Some.Contains("Self-dimer Tm"));
            // Fallback screen = Primer3 alignment-mode self_any (dpal.c: 20.00 for this palindrome) > 12.00.
            Assert.That(heuristic.Warnings, Has.Some.EqualTo("Self-complementarity: Primer3 self_any exceeds 12.00"));
        });
    }

    [Test]
    public void DesignProbes_LazyStructureScreen_EqualsExhaustiveRanking()
    {
        var rng = new Random(17);
        foreach (var (preset, length) in new[] { (ProbeDesigner.Defaults.qPCR, 300), (ProbeDesigner.Defaults.Microarray, 70) })
        {
            string seq = new(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
            var top = ProbeDesigner.DesignProbes(seq, preset, 10).ToList();
            var all = ProbeDesigner.DesignProbes(seq, preset, int.MaxValue).Take(10).ToList();
            Assert.That(top.Select(p => (p.Start, p.Sequence.Length, p.Score)),
                Is.EqualTo(all.Select(p => (p.Start, p.Sequence.Length, p.Score))));
        }
    }

    [Test]
    public void DesignMolecularBeacon_DetectionTemperature_AppliesTyagiKramerRules()
    {
        const string target = "ATGCGTACGTTAGCCTAGGCATCGATCGGATCCAGTTACGCAATTGCGATCGTA";
        var b = ProbeDesigner.DesignMolecularBeacon(target, 25, 5, detectionTemperatureCelsius: 50)!.Value;
        string beacon = b.Sequence;
        var hp = PrimerDesigner.CalculateHairpinThermodynamicsNtthal(beacon, 0.05, 0, 0)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(b.Tm, Is.EqualTo(PrimerDesigner.CalculateMeltingTemperaturePrimer3(
                beacon.Substring(5, 25), 50, 50, 0, 0)).Within(1e-12));
            Assert.That(b.Warnings, Has.Some.Contains($"Stem-loop (hairpin) Tm {hp.TmCelsius:F1}"));
            Assert.That(b.Warnings.Any(w => w.StartsWith("Stem-loop Tm")), Is.EqualTo(hp.TmCelsius < 57));
            Assert.That(b.Warnings.Any(w => w.StartsWith("Probe Tm")), Is.EqualTo(b.Tm < 57));
        });
    }

    [Test]
    public void ExtinctionCoefficientNearestNeighbor_MatchesCantorTable()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("ACGT"), Is.EqualTo(40300.0));
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("A"), Is.EqualTo(15400.0));
            // AA AA: 2·27400 − 15400 = 39400
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("aaa"), Is.EqualTo(39400.0));
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("ACGU", isDna: false), Is.EqualTo(41300.0));
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor("ACGU"), Is.NaN);
            Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor(""), Is.EqualTo(0.0));
        });
    }
}
