// PROBE-DESIGN-001 — Primer3 hybridization-probe picker, thermodynamic probe screens, probe Tm,
// molecular weight and nearest-neighbour extinction coefficient.
// Reference: primer3-py 2.3.1 design_primers(PRIMER_TASK = pick_hyb_probe_only), calc_tm;
//            Biopython 1.88 Bio.SeqUtils.molecular_weight; Cantor, Warshaw & Shapiro (1970) ε260 table.
using Seqeron.Genomics.Core;
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

    // ---- Assay presets: sourced Tm windows reachable on their own Tm scale (audit round 3, A3-10, B07.md F50) ----

    [Test]
    public void Defaults_Microarray_UsesOligoArrayConditionsAndTmWindow()
    {
        // OligoArray 2.0 (Rouillard, Zuker & Gulari 2003): NN Tm over the whole oligo at [Na+] 1 M, 1 µM oligo,
        // Tm window 82–90 °C; length 50–60 nt (Kane 2000 50-mers, Agilent 60-mers), G+C 40–60 %.
        var p = ProbeDesigner.Defaults.Microarray;
        Assert.Multiple(() =>
        {
            Assert.That((p.MinLength, p.MaxLength), Is.EqualTo((50, 60)));
            Assert.That((p.MinTm, p.MaxTm), Is.EqualTo((82.0, 90.0)));
            Assert.That((p.MinGc, p.MaxGc), Is.EqualTo((0.40, 0.60)));
            Assert.That(p.MonovalentMillimolar, Is.EqualTo(1000.0));
            Assert.That(p.DnaConcentrationNanomolar, Is.EqualTo(1000.0));
            Assert.That((p.DivalentMillimolar, p.DntpMillimolar), Is.EqualTo((0.0, 0.0)));
            Assert.That(p.MaxNearestNeighborLength, Is.EqualTo(60));
            // Other presets keep the Primer3 probe conditions and Primer3's MAX_NN_TM_LENGTH.
            foreach (var q in new[] { ProbeDesigner.Defaults.qPCR, ProbeDesigner.Defaults.FISH,
                         ProbeDesigner.Defaults.NorthernBlot, ProbeDesigner.Defaults.SouthernBlot })
            {
                Assert.That((q.DnaConcentrationNanomolar, q.MonovalentMillimolar, q.DivalentMillimolar, q.DntpMillimolar),
                    Is.EqualTo((50.0, 50.0, 0.0, 0.0)));
                Assert.That(q.MaxNearestNeighborLength, Is.EqualTo(36));
            }
        });
    }

    [Test]
    public void Defaults_OldMicroarrayWindow_WasUnreachableAtPrimer3Conditions()
    {
        // The defect: 75–85 °C at 50 nM / 50 mM (long_seq_tm for 50–60 nt) — primer3-py calc_tm(max_nn_length 36):
        // 60-mer with 60 % G+C (the window maximum) = 74.5029020719779 °C < 75.
        double max = PrimerDesigner.CalculateMeltingTemperaturePrimer3(new string('G', 36) + new string('A', 24), 50, 50, 0, 0);
        Assert.That(max, Is.EqualTo(74.5029020719779).Within(1e-9));
        Assert.That(max, Is.LessThan(75.0));
    }

    // Witness probes (length and G+C inside the preset window) whose Tm on the preset's own scale lies inside its Tm
    // window; expected Tm = primer3-py 2.3.1 calc_tm(seq, mv, dv 0, dntp 0, dna, max_nn_length).
    private static IEnumerable<TestCaseData> PresetWitnesses()
    {
        yield return new TestCaseData("Microarray", "AGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAAT", 86.25468810488348);
        yield return new TestCaseData("qPCR", "AGGAGCTTCACATCTGGCGCCGTGTGCCT", 68.58490905493761);
        yield return new TestCaseData("FISH", string.Concat(Enumerable.Repeat("ACGT", 50)), 77.4029020719779);
        yield return new TestCaseData("NorthernBlot", string.Concat(Enumerable.Repeat("ACGT", 25)), 74.4029020719779);
        yield return new TestCaseData("SouthernBlot", string.Concat(Enumerable.Repeat("AACGT", 40)), 73.30290207197791);
    }

    private static ProbeDesigner.ProbeParameters Preset(string name) => name switch
    {
        "Microarray" => ProbeDesigner.Defaults.Microarray,
        "qPCR" => ProbeDesigner.Defaults.qPCR,
        "FISH" => ProbeDesigner.Defaults.FISH,
        "NorthernBlot" => ProbeDesigner.Defaults.NorthernBlot,
        "SouthernBlot" => ProbeDesigner.Defaults.SouthernBlot,
        _ => throw new ArgumentException(name),
    };

    [TestCaseSource(nameof(PresetWitnesses))]
    public void Defaults_TmWindow_IsReachableForLengthAndGcWindow(string preset, string witness, double expectedTm)
    {
        var p = Preset(preset);
        double gc = witness.CalculateGcFractionFast();
        double tm = PrimerDesigner.CalculateMeltingTemperaturePrimer3(witness, p.DnaConcentrationNanomolar,
            p.MonovalentMillimolar, p.DivalentMillimolar, p.DntpMillimolar, p.MaxNearestNeighborLength);

        Assert.Multiple(() =>
        {
            Assert.That(witness.Length, Is.InRange(p.MinLength, p.MaxLength), "length window");
            Assert.That(gc, Is.InRange(p.MinGc, p.MaxGc), "G+C window");
            Assert.That(tm, Is.EqualTo(expectedTm).Within(1e-9), "primer3-py calc_tm at the preset conditions");
            Assert.That(tm, Is.InRange(p.MinTm, p.MaxTm), "Tm window reachable");
            // DesignProbes reports the same Tm and no Tm penalty for the witness window.
            var probe = ProbeDesigner.DesignProbes(witness, p with { MinLength = witness.Length, MaxLength = witness.Length }, 1).Single();
            Assert.That(probe.Tm, Is.EqualTo(expectedTm).Within(1e-9));
        });
    }

    // Exact attainable long_seq_tm range over the length × G+C window of the long-probe presets (the Tm depends only
    // on length and G+C count above MAX_NN_TM_LENGTH): primer3-py calc_tm(max_nn_length 36, mv 50, dv 0, dntp 0, dna 50)
    // minimised / maximised over every (N, #GC) in the window.
    [TestCase("FISH", 71.2529020719779, 85.35290207197791)]
    [TestCase("NorthernBlot", 70.30290207197791, 82.5029020719779)]
    [TestCase("SouthernBlot", 70.32012061502427, 85.35290207197791)]
    public void Defaults_LongProbePresets_AttainableTmRangeOverlapsWindow(string preset, double expectedMin, double expectedMax)
    {
        var p = Preset(preset);
        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        for (int n = p.MinLength; n <= p.MaxLength; n++)
        {
            for (int g = (int)Math.Ceiling(p.MinGc * n - 1e-9); g <= (int)Math.Floor(p.MaxGc * n + 1e-9); g++)
            {
                double tm = PrimerDesigner.CalculateMeltingTemperaturePrimer3(new string('G', g) + new string('A', n - g),
                    p.DnaConcentrationNanomolar, p.MonovalentMillimolar, p.DivalentMillimolar, p.DntpMillimolar,
                    p.MaxNearestNeighborLength);
                min = Math.Min(min, tm);
                max = Math.Max(max, tm);
            }
        }

        Assert.Multiple(() =>
        {
            Assert.That(min, Is.EqualTo(expectedMin).Within(1e-9));
            Assert.That(max, Is.EqualTo(expectedMax).Within(1e-9));
            Assert.That(Math.Max(min, p.MinTm), Is.LessThanOrEqualTo(Math.Min(max, p.MaxTm)), "window ∩ attainable range ≠ ∅");
        });
    }

    [Test]
    public void CalculateMeltingTemperaturePrimer3_MaxNearestNeighborLength_MatchesPrimer3SeqTm()
    {
        // primer3-py calc_tm(seq, mv 1000, dv 0, dntp 0, dna 1000, max_nn_length 60 / 36) — oligotm.c seqtm:
        // len > nn_max_len → long_seq_tm. Biopython Tm_NN(DNA_NN3, Na 1000, dnac1 = dnac2 = 500, saltcorr 0) = 86.2547.
        const string s = "AGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAAT";
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3(s, 1000, 1000, 0, 0, 60), Is.EqualTo(86.25468810488348).Within(1e-9));
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3(s, 1000, 1000, 0, 0, 36), Is.EqualTo(87.53999999999999).Within(1e-9));
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3(s, 50, 50, 0, 0, 36),
                Is.EqualTo(PrimerDesigner.CalculateMeltingTemperaturePrimer3(s, 50, 50, 0, 0)));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateMeltingTemperaturePrimer3(s, 50, 50, 0, 0, -1));
        });
    }

    [Test]
    public void DefaultConditions_OfValidateAnalyzeBeacon_StayPrimer3ProbeConditions()
    {
        // AnalyzeOligo keeps the Primer3 probe conditions (50 nM / 50 mM, NN ≤ 36 nt): the 50-mer is long_seq_tm.
        const string s = "AGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAAT";
        Assert.That(ProbeDesigner.AnalyzeOligo(s).Tm, Is.EqualTo(65.9429020719779).Within(1e-9));
    }
}
