// PRIMER-DESIGN-001 — Primer3 template mispriming (audit round 3, A3-4): PRIMER_MAX_TEMPLATE_MISPRIMING[_TH],
// PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING[_TH], PRIMER_WT_TEMPLATE_MISPRIMING[_TH], PRIMER_PAIR_WT_TEMPLATE_MISPRIMING[_TH],
// PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT.
// Source: primer3 libprimer3.c oligo_template_mispriming, primer_mispriming_to_template (dpal DPAL_LOCAL_END),
//         primer_mispriming_to_template_thermod (thal THAL_END1, use_end_for_th_template_mispriming = 1),
//         _pr_need_[pair_]template_mispriming[_thermod], calc_and_check_oligo_features (pick-time scoring when weighted),
//         characterize_pair (pair score max(left.T + right.T_r, left.T_r + right.T), PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH
//         "non-zero and exceeded" rule), p_obj_fn / obj_fn terms (linear; _TH: temp_cutoff rule), _pr_data_control.
// Expected values: primer3-py 2.3.1 design_primers (PRIMER_TASK check_primers for the alignment-mode single-primer maxima;
//         primer3.calc_end_stability = thal END1 for the thermodynamic per-strand values).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_TemplateMispriming_Tests
{
    // 240-nt random template made repetitive: rc(T[40:70]) copied to [170,200), T[150:175] (1 mutation) copied to [10,35).
    private const string T = "AGACTTTCAATGTCAAACCAATCTACCCCCTTCCGTATTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAATCTACCCCCTTCCGAAACACAATGAGAATTGGTAACAAACAGCGCAGCGGCAGATCAAGCAGGAGGCGGAATGTAAACA";
    private const int TargetStart = 100, TargetEnd = 120; // SEQUENCE_TARGET [100, 20]

    private static IReadOnlyList<PrimerPairResult> Design(PrimerParameters p, PrimerPairOptions o) =>
        PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), TargetStart, TargetEnd, p, o);

    // (left start, left len, right start (Primer3: 3'-end coordinate), right len, PAIR_PENALTY, LEFT_PENALTY, RIGHT_PENALTY,
    //  LEFT_TEMPLATE_MISPRIMING, RIGHT_TEMPLATE_MISPRIMING, PAIR_TEMPLATE_MISPRIMING[_TH]).
    private static void AssertPairs(IReadOnlyList<PrimerPairResult> pairs,
        (int L, int LLen, int R, int RLen, double Pen, double LPen, double RPen, double? LT, double? RT, double? PT)[] expected,
        bool comparePrimerValues = true)
    {
        Assert.That(pairs, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var e = expected[k];
                var p = pairs[k];
                Assert.That(p.Forward!.Position, Is.EqualTo(e.L), $"rank {k} left");
                Assert.That(p.Forward.Length, Is.EqualTo(e.LLen), $"rank {k} left length");
                Assert.That(p.Reverse!.Position + p.Reverse.Length - 1, Is.EqualTo(e.R), $"rank {k} right");
                Assert.That(p.Reverse.Length, Is.EqualTo(e.RLen), $"rank {k} right length");
                Assert.That(p.PairPenalty!.Value, Is.EqualTo(e.Pen).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(p.Forward.Penalty, Is.EqualTo(e.LPen).Within(1e-9), $"rank {k} PRIMER_LEFT_PENALTY");
                Assert.That(p.Reverse.Penalty, Is.EqualTo(e.RPen).Within(1e-9), $"rank {k} PRIMER_RIGHT_PENALTY");
                if (comparePrimerValues)
                {
                    Assert.That(p.Forward.TemplateMispriming, Is.EqualTo(e.LT), $"rank {k} PRIMER_LEFT_TEMPLATE_MISPRIMING");
                    Assert.That(p.Reverse.TemplateMispriming, Is.EqualTo(e.RT), $"rank {k} PRIMER_RIGHT_TEMPLATE_MISPRIMING");
                }
                if (e.PT is null)
                    Assert.That(p.TemplateMispriming, Is.Null, $"rank {k} PRIMER_PAIR_TEMPLATE_MISPRIMING");
                else
                    Assert.That(p.TemplateMispriming!.Value, Is.EqualTo(e.PT.Value).Within(1e-9), $"rank {k} PRIMER_PAIR_TEMPLATE_MISPRIMING");
            }
        });
    }

    [TestCase(58, 20, true, 9.0)]
    [TestCase(16, 20, true, 18.0)]   // the mutated copy at [150,175) primes again
    [TestCase(189, 20, false, 11.0)] // PRIMER_RIGHT [208,20]
    [TestCase(216, 20, false, 6.0)]  // PRIMER_RIGHT [235,20]
    public void CalculateTemplateMispriming_AlignmentMode_MatchesPrimer3CheckPrimers(int position, int length, bool forward, double expected)
    {
        // check_primers (SEQUENCE_PRIMER / SEQUENCE_PRIMER_REVCOMP, PRIMER_MAX_TEMPLATE_MISPRIMING 100):
        // PRIMER_LEFT/RIGHT_0_TEMPLATE_MISPRIMING = max(template_mispriming, template_mispriming_r).
        var s = PrimerDesigner.CalculateTemplateMispriming(new DnaSequence(T), position, length, forward);
        Assert.That(s.Max, Is.EqualTo(expected));
    }

    [TestCase(58, 20, true, 0.0, 21.936861234153127)]
    [TestCase(16, 20, true, 49.59451282825728, 0.0)]
    [TestCase(189, 20, false, 1.5815135678902834, 29.881724966610534)]
    [TestCase(216, 20, false, 12.246642332362342, 5.604136631924121)]
    public void CalculateTemplateMispriming_Thermodynamic_MatchesNtthalEnd1(int position, int length, bool forward,
        double sameStrand, double otherStrand)
    {
        // primer_mispriming_to_template_thermod: thal END1 of the primer against the complementary strand 5′ and 3′ of
        // its own site (max) and against the strand it lies on; reference primer3.calc_end_stability (thal END1, Primer3
        // default conditions) on the same segments.
        var s = PrimerDesigner.CalculateTemplateMispriming(new DnaSequence(T), position, length, forward, thermodynamic: true);
        Assert.Multiple(() =>
        {
            Assert.That(s.SameStrand, Is.EqualTo(sameStrand).Within(1e-9));
            Assert.That(s.OtherStrand, Is.EqualTo(otherStrand).Within(1e-9));
            Assert.That(s.Max, Is.EqualTo(Math.Max(sameStrand, otherStrand)).Within(1e-9));
        });
    }

    [Test]
    public void Defaults_DoNotComputeTemplateMispriming()
    {
        // design_primers with only SEQUENCE_TEMPLATE / SEQUENCE_TARGET: no PRIMER_*_TEMPLATE_MISPRIMING keys; rank 0 is
        // the left primer [16,20] that the alignment limit 12 rejects below (score 18).
        var pairs = Design(PrimerDesigner.Primer3DefaultParameters, PrimerPairOptions.Primer3Defaults);
        AssertPairs(pairs, new (int, int, int, int, double, double, double, double?, double?, double?)[]
        {
            (16, 20, 208, 20, 0.28312152398996204, 0.2530203683126615, 0.03010115567730054, null, null, null),
            (16, 20, 235, 20, 0.3607233439241213, 0.2530203683126615, 0.1077029756114598, null, null, null),
            (16, 20, 232, 20, 0.43205473658974825, 0.2530203683126615, 0.17903436827708674, null, null, null),
            (16, 20, 233, 20, 0.7181173090994548, 0.2530203683126615, 0.4650969407867933, null, null, null),
            (16, 20, 183, 20, 0.7341482757721565, 0.2530203683126615, 0.48112790745949496, null, null, null),
        });
    }

    [Test]
    public void AlignmentLimits_RejectMisprimingPrimersAndMatchPrimer3()
    {
        // PRIMER_MAX_TEMPLATE_MISPRIMING 12, PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING 24 (left "high template mispriming score 4").
        var p = PrimerDesigner.Primer3DefaultParameters with { MaxTemplateMispriming = 12 };
        var o = PrimerPairOptions.Primer3Defaults with { MaxTemplateMispriming = 24 };
        AssertPairs(Design(p, o), new (int, int, int, int, double, double, double, double?, double?, double?)[]
        {
            (58, 20, 208, 20, 1.1184070025607298, 1.0883058468834292, 0.03010115567730054, 9.0, 11.0, 15.0),
            (57, 20, 208, 20, 1.1184070025607298, 1.0883058468834292, 0.03010115567730054, 11.0, 11.0, 17.0),
            (58, 20, 235, 20, 1.196008822494889, 1.0883058468834292, 0.1077029756114598, 9.0, 6.0, 15.0),
            (57, 20, 235, 20, 1.196008822494889, 1.0883058468834292, 0.1077029756114598, 11.0, 6.0, 17.0),
            (58, 20, 232, 20, 1.267340215160516, 1.0883058468834292, 0.17903436827708674, 9.0, 5.0, 14.0),
        });
    }

    [Test]
    public void AlignmentWeights_EnterPenaltiesAndMatchPrimer3()
    {
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_MAX_TEMPLATE_MISPRIMING 9, PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING 16,
        // PRIMER_WT_TEMPLATE_MISPRIMING 0.5 (p_obj_fn), PRIMER_PAIR_WT_TEMPLATE_MISPRIMING 0.2 (obj_fn).
        var p = PrimerDesigner.Primer3DefaultParameters with
        {
            StructureScreen = PrimerStructureScreen.Primer3Alignment,
            MaxTemplateMispriming = 9,
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { TemplateMispriming = 0.5 },
        };
        var o = PrimerPairOptions.Primer3Defaults with
        {
            MaxTemplateMispriming = 16,
            Weights = new Primer3PairWeights { TemplateMispriming = 0.2 },
        };
        var expected = new (int, int, int, int, double, double, double, double?, double?, double?)[]
        {
            (59, 21, 237, 20, 10.14045393659626, 5.521996326128544, 2.618457610467715, 6.0, 4.0, 10.0),
            (58, 22, 237, 20, 10.413618704434054, 5.59516109396634, 2.618457610467715, 7.0, 4.0, 11.0),
            (58, 21, 237, 20, 10.487765937678057, 5.469308327210342, 2.618457610467715, 8.0, 4.0, 12.0),
            (59, 20, 237, 20, 10.50136176843365, 5.682904157965936, 2.618457610467715, 7.0, 4.0, 11.0),
            (60, 20, 237, 20, 10.542062364391768, 5.923604753924053, 2.618457610467715, 6.0, 4.0, 10.0),
        };
        AssertPairs(Design(p, o), expected);
        Assert.That(PrimerDesigner.DesignPrimers(new DnaSequence(T), TargetStart, TargetEnd, p, o).PairPenalty!.Value,
            Is.EqualTo(10.14045393659626).Within(1e-9));
    }

    [Test]
    public void ThermodynamicLimits_MatchPrimer3()
    {
        // PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT 1, PRIMER_MAX_TEMPLATE_MISPRIMING_TH 40,
        // PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH 70 → PRIMER_PAIR_k_TEMPLATE_MISPRIMING_TH (primer3-py does not report the
        // per-primer _TH values in this mode, thermoanalysis.pyx).
        var p = PrimerDesigner.Primer3DefaultParameters with { ThermodynamicTemplateAlignment = true, MaxTemplateMisprimingTh = 40 };
        var o = PrimerPairOptions.Primer3Defaults with { MaxTemplateMisprimingTh = 70 };
        var pairs = Design(p, o);
        AssertPairs(pairs, new (int, int, int, int, double, double, double, double?, double?, double?)[]
        {
            (58, 20, 208, 20, 1.1184070025607298, 1.0883058468834292, 0.03010115567730054, null, null, 29.881724966610534),
            (57, 20, 208, 20, 1.1184070025607298, 1.0883058468834292, 0.03010115567730054, null, null, 37.53906279399433),
            (58, 20, 235, 20, 1.196008822494889, 1.0883058468834292, 0.1077029756114598, null, null, 34.18350356651547),
            (57, 20, 235, 20, 1.196008822494889, 1.0883058468834292, 0.1077029756114598, null, null, 43.14456200410228),
            (58, 20, 232, 20, 1.267340215160516, 1.0883058468834292, 0.17903436827708674, null, null, 28.450035905362824),
        }, comparePrimerValues: false);
        // The pair value is max(left same + right other, left other + right same) of the per-primer END1 Tm values.
        Assert.That(pairs[0].Forward!.TemplateMispriming!.Value, Is.EqualTo(21.936861234153127).Within(1e-9));
        Assert.That(pairs[0].Reverse!.TemplateMispriming!.Value, Is.EqualTo(29.881724966610534).Within(1e-9));
    }

    [Test]
    public void ThermodynamicWeights_TempCutoffRule_MatchPrimer3()
    {
        // PRIMER_MAX_TEMPLATE_MISPRIMING_TH 45, PRIMER_WT_TEMPLATE_MISPRIMING_TH 0.5, PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH 0
        // (no pair limit), PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH 0.2: the per-primer and pair terms use the 5 °C temp_cutoff
        // rule (w / (Tm − 4 − s) below Tm − 5).
        var p = PrimerDesigner.Primer3DefaultParameters with
        {
            ThermodynamicTemplateAlignment = true,
            MaxTemplateMisprimingTh = 45,
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { TemplateMisprimingTh = 0.5 },
        };
        var o = PrimerPairOptions.Primer3Defaults with
        {
            MaxTemplateMisprimingTh = 0,
            Weights = new Primer3PairWeights { TemplateMisprimingTh = 0.2 },
        };
        AssertPairs(Design(p, o), new (int, int, int, int, double, double, double, double?, double?, double?)[]
        {
            (58, 20, 208, 20, 1.1607262752978076, 1.1034689259876893, 0.049266927967217836, null, null, 29.881724966610534),
            (57, 20, 208, 20, 1.169906518813047, 1.1091272300603123, 0.049266927967217836, null, null, 37.53906279399433),
            (58, 20, 235, 20, 1.2322202302654168, 1.1034689259876893, 0.11910260874087336, null, null, 34.18350356651547),
            (57, 20, 235, 20, 1.2452263335687535, 1.1091272300603123, 0.11910260874087336, null, null, 43.14456200410228),
            (58, 20, 232, 20, 1.3001286772497616, 1.1034689259876893, 0.18910164594125498, null, null, 28.450035905362824),
        }, comparePrimerValues: false);
    }

    [Test]
    public void ThermodynamicPairWeightWithDefaultPairLimit_FailsEveryPair_AsPrimer3()
    {
        // characterize_pair: `if (pair_max_template_mispriming_th && score > pair_max_template_mispriming_th)` — with the
        // default −100 and PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH 0.2 every pair fails (PRIMER_PRODUCT_SIZE_RANGE 150-152:
        // "high template mispriming score 3104, ok 0"); PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH 0 lifts the limit.
        var p = PrimerDesigner.Primer3DefaultParameters with { ThermodynamicTemplateAlignment = true };
        var o = PrimerPairOptions.Primer3Defaults with
        {
            ProductSizeRanges = [new ProductSizeRange(150, 152)],
            Weights = new Primer3PairWeights { TemplateMisprimingTh = 0.2 },
        };
        var none = Design(p, o);
        var failure = PrimerDesigner.DesignPrimers(new DnaSequence(T), TargetStart, TargetEnd, p, o);
        Assert.Multiple(() =>
        {
            Assert.That(none, Is.Empty);
            Assert.That(failure.IsValid, Is.False);
            Assert.That(failure.Message, Does.Contain("template-mispriming"));
        });
        AssertPairs(Design(p, o with { MaxTemplateMisprimingTh = 0 }), new (int, int, int, int, double, double, double, double?, double?, double?)[]
        {
            (58, 20, 208, 20, 1.1263974239036303, 1.0883058468834292, 0.03010115567730054, null, null, 29.881724966610534),
            (57, 20, 208, 20, 1.1299193633462468, 1.0883058468834292, 0.03010115567730054, null, null, 37.53906279399433),
            (57, 21, 208, 20, 1.5072070679635288, 1.4693083272103422, 0.03010115567730054, null, null, 29.881724966610534),
            (58, 21, 208, 20, 1.509101954044297, 1.4693083272103422, 0.03010115567730054, null, null, 34.896119047552304),
            (57, 22, 208, 20, 2.1347527160871898, 2.0951610939663396, 0.03010115567730054, null, null, 34.896119047552304),
        }, comparePrimerValues: false);
    }

    [Test]
    public void EvaluatePrimer_WithoutTemplate_HasNoTemplateTerm()
    {
        // No template → no template mispriming score (Primer3 never computes it for a primer alone); the weights are inert.
        var plain = PrimerDesigner.EvaluatePrimer("TGTGTTTCGGAACTTGCGTT", 58, true, PrimerDesigner.Primer3DefaultParameters);
        var weighted = PrimerDesigner.EvaluatePrimer("TGTGTTTCGGAACTTGCGTT", 58, true, PrimerDesigner.Primer3DefaultParameters with
        {
            ThermodynamicTemplateAlignment = true,
            MaxTemplateMisprimingTh = 10,
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { TemplateMisprimingTh = 1.0 },
        });
        Assert.Multiple(() =>
        {
            Assert.That(weighted.Penalty, Is.EqualTo(plain.Penalty));
            Assert.That(weighted.TemplateMispriming, Is.Null);
            Assert.That(weighted.IsValid, Is.True);
        });
    }

    [Test]
    public void CalculatePrimer3Penalty_TemplateTerms()
    {
        // p_obj_fn: alignment mode + w·s; thermodynamic mode s ≥ Tm − 5 → w·(s − (Tm − 6)), else w / (Tm − 4 − s).
        var inputs = new Primer3PenaltyInputs(60.0, 20, 50.0, 0, 0, 0, 0, 0) { TemplateMispriming = 10.0 };
        var w = new Primer3PenaltyWeights { TemplateMispriming = 0.5, TemplateMisprimingTh = 2.0 };
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs, w), Is.EqualTo(5.0).Within(1e-12));
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs, w with { ThermodynamicTemplateAlignment = true }),
                Is.EqualTo(2.0 / (60.0 - 4.0 - 10.0)).Within(1e-12));
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs with { TemplateMispriming = 57.0 },
                w with { ThermodynamicTemplateAlignment = true }), Is.EqualTo(2.0 * (57.0 - 54.0)).Within(1e-12));
        });
    }

    [Test]
    public void DataControl_IllegalSettings_Throw()
    {
        var dna = new DnaSequence(T);
        var p = PrimerDesigner.Primer3DefaultParameters;
        var o = PrimerPairOptions.Primer3Defaults;
        Assert.Multiple(() =>
        {
            // _pr_data_control: "Value too large at tag PRIMER_MAX_TEMPLATE_MISPRIMING" / PRIMER_PAIR_… (alignment mode).
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd,
                p with { MaxTemplateMispriming = 40000 }, o));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd,
                p, o with { MaxTemplateMispriming = 40000 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd,
                p with { PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { TemplateMispriming = -1 } }, o));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd,
                p, o with { Weights = new Primer3PairWeights { TemplateMisprimingTh = -1 } }));
            // thal: template > THAL_MAX_SEQ (10000) in thermodynamic template alignment.
            Assert.Throws<ArgumentException>(() => PrimerDesigner.CalculateTemplateMispriming(
                new DnaSequence(string.Concat(Enumerable.Repeat("ACGTTGCA", 1251))), 0, 20, true, thermodynamic: true));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateTemplateMispriming(dna, 230, 20, true));
            // Library scored at pick time (PRIMER_WT_LIBRARY_MISPRIMING ≠ 0) + a pair template limit without a template
            // weight: Primer3 never scores the primers' template mispriming and aborts (PR_ASSERT) — reported as an error.
            var lib = new PrimerMisprimingLibrary(new[] { KeyValuePair.Create("x", "GGGGGGGGGGCCCCCCCCCC") });
            Assert.Throws<InvalidOperationException>(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd,
                p with
                {
                    MisprimingLibrary = lib,
                    PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { LibraryMispriming = 0.5 },
                }, o with { MaxTemplateMispriming = 20 }));
        });
    }
}
