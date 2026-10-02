// PRIMER-DESIGN-001 — Primer3 reaction conditions and GC optimum in EvaluatePrimer / DesignPrimers (audit round 3, A3-1):
// PRIMER_SALT_MONOVALENT / _SALT_DIVALENT / _DNTP_CONC / _DNA_CONC / _OPT_GC_PERCENT and the internal-oligo
// PRIMER_INTERNAL_OPT_GC_PERCENT / _WT_GC_PERCENT_GT / _LT.
// Source: primer3 libprimer3.c — calc_and_check_oligo_features (seqtm with po_args conditions), create_thal_arg_holder(p_args)
//         (ntthal self-any / self-end / hairpin and the characterize_pair compl_any / compl_end use the primer conditions),
//         characterize_pair long_seq_tm(p_args salt/divalent/dNTP) for PRIMER_PAIR_k_PRODUCT_TM, p_obj_fn GC terms around
//         p_args/o_args.opt_gc_content, _pr_data_control ("Illegal value for primer salt or dna concentration", "… divalent salt
//         or dNTP concentration").
// Expected values: primer3-py 2.3.1 design_primers (PRIMER_LEFT/RIGHT/PAIR/INTERNAL_k_*).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_ReactionConditions_Tests
{
    private const string T = "GGTCCTAATTGGAGCGCCCAGTTACCGGCCGAGTGCTACGGGCACTCGTTGGTAGTGGGCTCCCTAAGTCGGCGCATCCGTTCCTAGCTTTAAAATATCCGTTGAAAGAATGTTCTGAGTCTCGCCTAGTGAAAGCCAACTCCTTTGGATTTTGTCATA";
    private const string A1 = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";

    // PRIMER_SALT_MONOVALENT 20, PRIMER_SALT_DIVALENT 3.0, PRIMER_DNTP_CONC 0.8, PRIMER_DNA_CONC 250,
    // PRIMER_OPT_GC_PERCENT 45, PRIMER_WT_GC_PERCENT_GT 0.5, PRIMER_WT_GC_PERCENT_LT 1.0.
    private static PrimerParameters Conditions(PrimerStructureScreen screen = PrimerStructureScreen.Primer3Thermodynamic) =>
        PrimerDesigner.Primer3DefaultParameters with
        {
            StructureScreen = screen,
            MonovalentMillimolar = 20.0,
            DivalentMillimolar = 3.0,
            DntpMillimolar = 0.8,
            DnaConcentrationNanomolar = 250.0,
            OptimalGcPercent = 45.0,
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 0.5, GcLt = 1.0 },
        };

    [Test]
    public void EvaluatePrimer_Conditions_MatchPrimer3TmStructureAndPenalty()
    {
        // design_primers(T, target [54,21], conditions above): PRIMER_LEFT_0 TACGGGCACTCGTTGGTA.
        var c = PrimerDesigner.EvaluatePrimer("TACGGGCACTCGTTGGTA", 36, true, Conditions());
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3("TACGGGCACTCGTTGGTA", 250.0, 20.0, 3.0, 0.8),
                Is.EqualTo(61.39018622811568).Within(1e-9), "PRIMER_LEFT_0_TM");
            Assert.That(c.MeltingTemperature, Is.EqualTo(61.4));
            Assert.That(c.Penalty, Is.EqualTo(8.66796400589346).Within(1e-9), "PRIMER_LEFT_0_PENALTY (Tm 1.39 + size 2 + GC 0.5·10.56)");
            Assert.That(c.SelfAnyTh!.Value, Is.EqualTo(4.700639676931019).Within(1e-9), "PRIMER_LEFT_0_SELF_ANY_TH");
            Assert.That(c.SelfEndTh!.Value, Is.EqualTo(0.0).Within(1e-12), "PRIMER_LEFT_0_SELF_END_TH");
            Assert.That(c.HairpinTh!.Value, Is.EqualTo(44.063749383794345).Within(1e-9), "PRIMER_LEFT_0_HAIRPIN_TH");
            Assert.That(c.IsValid, Is.True);
        });
    }

    private static readonly (int L, int LL, int R, int RL, double Pen, double LPen, double RPen, double ProductTm)[] ThermoConditions =
    {
        (36, 18, 146, 20, 9.27131815258577, 8.66796400589346, 0.6033541466923111, 84.7341429364544),
        (2, 18, 146, 20, 9.703063401178618, 9.099709254486307, 0.6033541466923111, 87.17595716447552),
        (2, 18, 127, 20, 9.853866427986254, 9.099709254486307, 0.7541571734999479, 86.86391556622704),
        (2, 18, 124, 20, 10.404866002928095, 9.099709254486307, 1.3051567484417887, 86.94618424605669),
        (36, 18, 144, 20, 10.443121582149622, 8.66796400589346, 1.7751575762561629, 84.63157247608576),
    };

    private static readonly (int L, int LL, int R, int RL, double Pen, double LPen, double RPen, double ProductTm)[] AlignmentConditions =
    {
        (74, 20, 267, 20, 1.3224169404189752, 0.46963265989268166, 0.8527842805262935, 87.44278972757542),
        (13, 20, 267, 20, 1.515246045265883, 0.6624617647395894, 0.8527842805262935, 87.73011537948567),
        (14, 20, 267, 20, 1.8050239469779967, 0.9522396664517032, 0.8527842805262935, 87.80061097930039),
        (74, 20, 246, 20, 1.8293329050669627, 0.46963265989268166, 1.359700245174281, 87.44851048319188),
        (13, 20, 281, 20, 1.8323347965294374, 0.6624617647395894, 1.169873031789848, 87.71270886297572),
    };

    private static void AssertPairs(IReadOnlyList<PrimerPairResult> pairs,
        (int L, int LL, int R, int RL, double Pen, double LPen, double RPen, double ProductTm)[] expected)
    {
        Assert.That(pairs, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var e = expected[k];
                Assert.That(pairs[k].Forward!.Position, Is.EqualTo(e.L), $"rank {k} left start");
                Assert.That(pairs[k].Forward!.Length, Is.EqualTo(e.LL), $"rank {k} left length");
                Assert.That(pairs[k].Reverse!.Position + pairs[k].Reverse!.Length - 1, Is.EqualTo(e.R), $"rank {k} right 3′ end");
                Assert.That(pairs[k].Reverse!.Length, Is.EqualTo(e.RL), $"rank {k} right length");
                Assert.That(pairs[k].PairPenalty!.Value, Is.EqualTo(e.Pen).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(pairs[k].Forward!.Penalty, Is.EqualTo(e.LPen).Within(1e-9), $"rank {k} PRIMER_LEFT_PENALTY");
                Assert.That(pairs[k].Reverse!.Penalty, Is.EqualTo(e.RPen).Within(1e-9), $"rank {k} PRIMER_RIGHT_PENALTY");
                Assert.That(pairs[k].ProductTm!.Value, Is.EqualTo(e.ProductTm).Within(1e-9), $"rank {k} PRIMER_PAIR_PRODUCT_TM");
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_ThermodynamicMode_Conditions_MatchPrimer3()
    {
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, Conditions(), PrimerPairOptions.Primer3Defaults);
        AssertPairs(pairs, ThermoConditions);
        Assert.Multiple(() =>
        {
            // PRIMER_RIGHT_0 CAAAGGAGTTGGCTTTCACT and PRIMER_PAIR_0_COMPL_ANY_TH / _COMPL_END_TH at the same conditions.
            Assert.That(pairs[0].Reverse!.Sequence, Is.EqualTo("CAAAGGAGTTGGCTTTCACT"));
            Assert.That(pairs[0].Reverse!.HairpinTh!.Value, Is.EqualTo(33.768163014786126).Within(1e-9), "PRIMER_RIGHT_0_HAIRPIN_TH");
            Assert.That(pairs[0].ComplAnyTh!.Value, Is.EqualTo(9.358916559593922).Within(1e-9), "PRIMER_PAIR_0_COMPL_ANY_TH");
            Assert.That(pairs[0].ComplEndTh!.Value, Is.EqualTo(0.0).Within(1e-12), "PRIMER_PAIR_0_COMPL_END_TH");
            // Product Tm = long_seq_tm with the primer conditions (PRIMER_SALT_* / PRIMER_DNTP_CONC).
            string product = T.Substring(36, 146 - 36 + 1);
            Assert.That(PrimerDesigner.CalculateProductMeltingTemperaturePrimer3(product, 20.0, 3.0, 0.8),
                Is.EqualTo(84.7341429364544).Within(1e-9));
        });
    }

    [Test]
    public void DesignPrimerPairs_AlignmentMode_Conditions_MatchPrimer3()
    {
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(A1), 100, 130,
            Conditions(PrimerStructureScreen.Primer3Alignment), PrimerPairOptions.Primer3Defaults);
        AssertPairs(pairs, AlignmentConditions);
        Assert.Multiple(() =>
        {
            Assert.That(pairs[0].Forward!.Sequence, Is.EqualTo("TGGTGATCCTATGCTTGTGA"));
            Assert.That(pairs[0].Forward!.SelfAny, Is.EqualTo(4.0));
            Assert.That(pairs[0].Reverse!.SelfAny, Is.EqualTo(6.0));
            Assert.That(pairs[0].ComplAny, Is.EqualTo(6.0));
            Assert.That(pairs[0].ComplEnd, Is.EqualTo(2.0));
        });
    }

    [Test]
    public void DesignPrimerPairs_InternalOligoConditionsAndGcOptimum_MatchPrimer3()
    {
        // + PRIMER_PICK_INTERNAL_OLIGO 1, PRIMER_INTERNAL_SALT_MONOVALENT 120, _SALT_DIVALENT 2, _DNTP_CONC 0.2, _DNA_CONC 100,
        //   PRIMER_INTERNAL_OPT_GC_PERCENT 40, PRIMER_INTERNAL_WT_GC_PERCENT_GT 0.5, _LT 0.5.
        var intl = new ProbeDesigner.Primer3ProbeSettings(MonovalentMillimolar: 120.0, DivalentMillimolar: 2.0,
            DntpMillimolar: 0.2, DnaConcentrationNanomolar: 100.0)
        { OptGcPercent = 40.0, WeightGcPercentGt = 0.5, WeightGcPercentLt = 0.5 };
        var opts = PrimerPairOptions.Primer3Defaults with { PickInternalOligo = true, InternalOligo = intl };
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(A1), 100, 130, Conditions(), opts);

        var expected = new (int L, int R, double Pen, double LPen, double RPen, double ProductTm)[]
        {
            (74, 278, 0.7052466341623926, 0.46963265989268166, 0.23561397426971098, 87.49740375825182),
            (74, 272, 0.8842756639249956, 0.46963265989268166, 0.4146430040323139, 87.41217272503695),
            (13, 278, 0.8980757390093004, 0.6624617647395894, 0.23561397426971098, 87.76032325210839),
            (13, 272, 1.0771047687719033, 0.6624617647395894, 0.4146430040323139, 87.70115610346757),
            (14, 278, 1.1878536407214142, 0.9522396664517032, 0.23561397426971098, 87.82800661145015),
        };
        Assert.That(pairs, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var e = expected[k];
                Assert.That(pairs[k].Forward!.Position, Is.EqualTo(e.L), $"rank {k} left start");
                Assert.That(pairs[k].Reverse!.Position + pairs[k].Reverse!.Length - 1, Is.EqualTo(e.R), $"rank {k} right 3′ end");
                Assert.That(pairs[k].PairPenalty!.Value, Is.EqualTo(e.Pen).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(pairs[k].ProductTm!.Value, Is.EqualTo(e.ProductTm).Within(1e-9), $"rank {k} PRIMER_PAIR_PRODUCT_TM");
                // PRIMER_INTERNAL_k = 122,20 (GC 40 % = the optimum) for every rank; with optimum 55 % Primer3 picks 178.
                var io = pairs[k].InternalOligo!.Value;
                Assert.That(io.Start, Is.EqualTo(122), $"rank {k} PRIMER_INTERNAL start");
                Assert.That(io.Length, Is.EqualTo(20), $"rank {k} PRIMER_INTERNAL length");
                Assert.That(io.Penalty, Is.EqualTo(2.1790501883629076).Within(1e-9), $"rank {k} PRIMER_INTERNAL_PENALTY");
                Assert.That(io.Tm, Is.EqualTo(57.82094981163709).Within(1e-9), $"rank {k} PRIMER_INTERNAL_TM");
            }
            // PRIMER_RIGHT_0 / PRIMER_PAIR_0 structure values at the primer conditions.
            Assert.That(pairs[0].Forward!.HairpinTh!.Value, Is.EqualTo(34.537937050227754).Within(1e-9), "PRIMER_LEFT_0_HAIRPIN_TH");
            Assert.That(pairs[0].ComplAnyTh!.Value, Is.EqualTo(3.421368805685802).Within(1e-9), "PRIMER_PAIR_0_COMPL_ANY_TH");
            Assert.That(pairs[0].ComplEndTh!.Value, Is.EqualTo(17.360051093275388).Within(1e-9), "PRIMER_PAIR_0_COMPL_END_TH");
        });
    }

    [Test]
    public void DesignProbesPrimer3_InternalGcOptimumAndWeights_MatchPrimer3()
    {
        // PRIMER_TASK=pick_hyb_probe_only with the internal conditions / GC optimum of the previous test.
        var s = new ProbeDesigner.Primer3ProbeSettings(MonovalentMillimolar: 120.0, DivalentMillimolar: 2.0,
            DntpMillimolar: 0.2, DnaConcentrationNanomolar: 100.0)
        { OptGcPercent = 40.0, WeightGcPercentGt = 0.5, WeightGcPercentLt = 0.5 };
        var probes = ProbeDesigner.DesignProbesPrimer3(A1, s);
        var expected = new (int Start, double Penalty, double Tm)[]
        {
            (261, 1.1685803999403106, 58.83141960005969),
            (260, 1.3931870962722428, 58.60681290372776),
            (267, 1.974124850932185, 58.025875149067815),
            (122, 2.1790501883629076, 57.82094981163709),
            (16, 2.2008397115195066, 57.79916028848049),
        };
        Assert.That(probes, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                Assert.That(probes[k].Start, Is.EqualTo(expected[k].Start), $"rank {k} start");
                Assert.That(probes[k].Length, Is.EqualTo(20), $"rank {k} length");
                Assert.That(probes[k].Penalty, Is.EqualTo(expected[k].Penalty).Within(1e-9), $"rank {k} PRIMER_INTERNAL_PENALTY");
                Assert.That(probes[k].Tm, Is.EqualTo(expected[k].Tm).Within(1e-9), $"rank {k} PRIMER_INTERNAL_TM");
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_NullConditions_EqualPrimer3Defaults()
    {
        // Default conditions (50 mM / 1.5 mM / 0.6 mM / 50 nM): unchanged primer3-py default ranking on T
        // (PRIMER_PAIR_0_PENALTY 1.7645480870700112, PRODUCT_TM 84.91708592812545).
        var explicitDefaults = PrimerDesigner.Primer3DefaultParameters with
        {
            MonovalentMillimolar = PrimerDesigner.Primer3MonovalentMillimolar,
            DivalentMillimolar = PrimerDesigner.Primer3DivalentMillimolar,
            DntpMillimolar = PrimerDesigner.Primer3DntpMillimolar,
            DnaConcentrationNanomolar = PrimerDesigner.Primer3DnaConcentrationNanomolar,
            OptimalGcPercent = 50.0,
        };
        var implicitPairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75,
            PrimerDesigner.Primer3DefaultParameters, PrimerPairOptions.Primer3Defaults);
        var explicitPairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, explicitDefaults, PrimerPairOptions.Primer3Defaults);
        Assert.Multiple(() =>
        {
            Assert.That(implicitPairs[0].PairPenalty!.Value, Is.EqualTo(1.7645480870700112).Within(1e-9));
            Assert.That(implicitPairs[0].ProductTm!.Value, Is.EqualTo(84.91708592812545).Within(1e-9));
            Assert.That(implicitPairs[0].Forward!.Sequence, Is.EqualTo("GTTACCGGCCGAGTGCTAC"));
            Assert.That(explicitPairs.Select(p => (p.Forward!.Position, p.Reverse!.Position, p.PairPenalty)),
                Is.EqualTo(implicitPairs.Select(p => (p.Forward!.Position, p.Reverse!.Position, p.PairPenalty))));
            Assert.That(PrimerDesigner.Primer3DefaultParameters.EffectiveMonovalentMillimolar, Is.EqualTo(50.0));
            Assert.That(PrimerDesigner.Primer3DefaultParameters.EffectiveDivalentMillimolar, Is.EqualTo(1.5));
            Assert.That(PrimerDesigner.Primer3DefaultParameters.EffectiveDntpMillimolar, Is.EqualTo(0.6));
            Assert.That(PrimerDesigner.Primer3DefaultParameters.EffectiveDnaConcentrationNanomolar, Is.EqualTo(50.0));
            Assert.That(PrimerDesigner.Primer3DefaultParameters.OptimalGcPercent, Is.Null,
                "PRIMER_OPT_GC_PERCENT undefined by default (libprimer3.c DEFAULT_OPT_GC_PERCENT = PR_UNDEFINED_INT_OPT).");
        });
    }

    [Test]
    public void OptimalGcPercent_InertWithoutGcWeights()
    {
        // p_obj_fn adds the GC terms only for a non-zero PRIMER_WT_GC_PERCENT_GT/_LT (default 0).
        var a = PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 20, true, PrimerDesigner.Primer3DefaultParameters);
        var b = PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 20, true,
            PrimerDesigner.Primer3DefaultParameters with { OptimalGcPercent = 30.0 });
        Assert.That(b.Penalty, Is.EqualTo(a.Penalty).And.EqualTo(1.5179132593960958).Within(1e-9));
    }

    private static IEnumerable<TestCaseData> IllegalConditions()
    {
        // _pr_data_control: salt ≤ 0 or DNA ≤ 0 → "Illegal value for primer salt or dna concentration";
        // divalent < 0 → "… divalent salt or dNTP concentration"; NaN / ∞ / negative dNTP also rejected.
        yield return new TestCaseData(0.0, 1.5, 0.6, 50.0).SetName("Monovalent0");
        yield return new TestCaseData(-5.0, 1.5, 0.6, 50.0).SetName("MonovalentNegative");
        yield return new TestCaseData(50.0, 1.5, 0.6, 0.0).SetName("Dna0");
        yield return new TestCaseData(50.0, -1.0, 0.6, 50.0).SetName("DivalentNegative");
        yield return new TestCaseData(50.0, 1.5, -0.1, 50.0).SetName("DntpNegative");
        yield return new TestCaseData(double.NaN, 1.5, 0.6, 50.0).SetName("MonovalentNaN");
        yield return new TestCaseData(50.0, 1.5, 0.6, double.PositiveInfinity).SetName("DnaInfinity");
    }

    [TestCaseSource(nameof(IllegalConditions))]
    public void IllegalConditions_Throw(double mv, double dv, double dntp, double dna)
    {
        var p = PrimerDesigner.Primer3DefaultParameters with
        {
            MonovalentMillimolar = mv, DivalentMillimolar = dv, DntpMillimolar = dntp, DnaConcentrationNanomolar = dna,
        };
        var intl = new ProbeDesigner.Primer3ProbeSettings(MonovalentMillimolar: mv, DivalentMillimolar: dv,
            DntpMillimolar: dntp, DnaConcentrationNanomolar: dna);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 0, true, p));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, p, PrimerPairOptions.Primer3Defaults));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, PrimerDesigner.Primer3DefaultParameters,
                    PrimerPairOptions.Primer3Defaults with { PickInternalOligo = true, InternalOligo = intl }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbesPrimer3(A1, intl));
        });
    }

    [Test]
    public void DivalentWithoutDntp_Accepted()
    {
        // Primer3 only warns for PRIMER_SALT_DIVALENT > 0 with PRIMER_DNTP_CONC = 0 (free Mg²⁺ = all of it).
        var p = PrimerDesigner.Primer3DefaultParameters with { DntpMillimolar = 0.0 };
        Assert.DoesNotThrow(() => PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 0, true, p));
        Assert.That(PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 0, true, p).MeltingTemperature,
            Is.EqualTo(Math.Round(PrimerDesigner.CalculateMeltingTemperaturePrimer3("GTTACCGGCCGAGTGCTAC", 50, 50, 1.5, 0), 1)));
    }
}
