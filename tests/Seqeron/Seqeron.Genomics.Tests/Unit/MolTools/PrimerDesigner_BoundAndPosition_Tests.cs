// PRIMER-DESIGN-001 / PROBE-DESIGN-001 — Primer3 p_obj_fn terms bound and pos_penalty (audit round 3, A3-5 part 1):
// PRIMER_ANNEALING_TEMP, PRIMER_MIN/MAX/OPT_BOUND, PRIMER_WT_BOUND_GT/LT, PRIMER_INTERNAL_MIN/MAX/OPT_BOUND,
// PRIMER_INTERNAL_WT_BOUND_GT/LT, PRIMER_INSIDE_PENALTY / PRIMER_OUTSIDE_PENALTY, PRIMER_WT_POS_PENALTY.
// Source: primer3-py 2.3.1 libprimer3.c calc_and_check_oligo_features (bound check after the Tm checks when
//         annealing_temp > 0; position penalty / overlap-target logic), compute_position_penalty, pair_spans_target,
//         make_detection_primer_lists (search region with non-default position penalties), p_obj_fn (primer bound terms
//         gated by annealing_temp > 0, internal-oligo bound terms ungated; pos_penalty term), obj_fn PR_ASSERT(sum >= 0),
//         _pr_data_control; oligotm.c oligotm (bound, santalucia salt correction), seqtm (MAX_NN_TM_LENGTH 36).
// Expected values: primer3-py 2.3.1 design_primers (PRIMER_TASK check_primers for single primers, pick_hyb_probe_only
//         for probes, pair designs on the 300-nt template below with SEQUENCE_TARGET [140, 20]).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_BoundAndPosition_Tests
{
    // random.Random(5), 300 nt.
    private const string T = "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG";
    private const int TargetStart = 140, TargetEnd = 160; // SEQUENCE_TARGET [140, 20]

    // random.Random(11), 120 nt (pick_hyb_probe_only).
    private const string P = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC";

    private static readonly PrimerParameters P3 = PrimerDesigner.Primer3DefaultParameters;
    private static readonly PrimerPairOptions O3 = PrimerPairOptions.Primer3Defaults;

    private static IReadOnlyList<PrimerPairResult> Design(PrimerParameters p, PrimerPairOptions o) =>
        PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), TargetStart, TargetEnd, p, o);

    private static PrimerParameters WithWeights(PrimerParameters p, Func<Primer3PenaltyWeights, Primer3PenaltyWeights> f) =>
        p with { PenaltyWeights = f(p.PenaltyWeights ?? PrimerDesigner.DefaultPrimer3Weights) };

    // (left start, left len, right 3'-end coordinate (Primer3 PRIMER_RIGHT start), right len, PAIR_PENALTY, LEFT_PENALTY,
    //  RIGHT_PENALTY, LEFT_BOUND, RIGHT_BOUND, LEFT_POSITION_PENALTY, RIGHT_POSITION_PENALTY).
    private static void AssertPairs(IReadOnlyList<PrimerPairResult> pairs,
        (int L, int LLen, int R, int RLen, double Pen, double LPen, double RPen, double? LB, double? RB, double? LP, double? RP)[] expected)
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
                AssertNullable(p.Forward.Bound, e.LB, $"rank {k} PRIMER_LEFT_BOUND");
                AssertNullable(p.Reverse.Bound, e.RB, $"rank {k} PRIMER_RIGHT_BOUND");
                AssertNullable(p.Forward.PositionPenalty, e.LP, $"rank {k} PRIMER_LEFT_POSITION_PENALTY");
                AssertNullable(p.Reverse.PositionPenalty, e.RP, $"rank {k} PRIMER_RIGHT_POSITION_PENALTY");
            }
        });
    }

    private static void AssertNullable(double? actual, double? expected, string what)
    {
        if (expected is null)
            Assert.That(actual, Is.Null, what);
        else
            Assert.That(actual!.Value, Is.EqualTo(expected.Value).Within(1e-9), what);
    }

    // The default design (no annealing temperature, default position penalties) of T.
    private static readonly (int, int, int, int, double, double, double, double?, double?, double?, double?)[] DefaultRanks =
    [
        (46, 20, 233, 20, 0.2229611791929642, 0.1161425661219937, 0.10681861307097051, null, null, null, null),
        (116, 20, 233, 20, 0.5002207880304468, 0.39340217495947627, 0.10681861307097051, null, null, null, null),
        (117, 20, 233, 20, 0.5003559657240544, 0.39353735265308387, 0.10681861307097051, null, null, null, null),
        (46, 20, 240, 20, 0.5299746133399026, 0.1161425661219937, 0.4138320472179089, null, null, null, null),
        (46, 20, 236, 20, 0.5739934657308936, 0.1161425661219937, 0.4578508996088999, null, null, null, null),
    ];

    #region Fraction bound of one oligo (oligotm)

    [TestCase("AGCTAGCTAGCTAGCTAGCT", 60.0, 50.0, 50.0, 1.5, 0.6, 34.19287577114419)]
    [TestCase("AGCTAGCTAGCTAGCTAGCT", 50.0, 50.0, 50.0, 1.5, 0.6, 94.6818573933472)]
    [TestCase("ACGTACGTACGTACGTACGT", 55.0, 50.0, 50.0, 1.5, 0.6, 84.65679613064577)]   // self-complementary: C/1e9 (Equation A)
    [TestCase("GGCGCGCCATGGCGCGCC", 72.0, 50.0, 50.0, 1.5, 0.6, 68.75377597587205)]     // self-complementary
    [TestCase("TTGCCAGTACCATTGACGTA", 58.0, 250.0, 100.0, 3.0, 0.8, 79.18943387635454)]
    [TestCase("ATATATTTAAATATATTA", 30.0, 50.0, 50.0, 1.5, 0.6, 47.84518704158839)]
    [TestCase("CCGTAGCTAGGCTTAGCAGGCATCGATCGATGCAT", 65.0, 50.0, 50.0, 1.5, 0.6, 99.27112164830403)] // 35 nt
    public void CalculateFractionBoundPrimer3_MatchesPrimer3(string primer, double ta, double dna, double mono, double div,
        double dntp, double expected)
    {
        // primer3-py check_primers (SEQUENCE_PRIMER, PRIMER_ANNEALING_TEMP = ta) → PRIMER_LEFT_0_BOUND.
        Assert.That(PrimerDesigner.CalculateFractionBoundPrimer3(primer, ta, dna, mono, div, dntp),
            Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void CalculateFractionBoundPrimer3_UndefinedAndIllegal()
    {
        Assert.Multiple(() =>
        {
            // seqtm: > 36 nt → long_seq_tm, bound stays OLIGOTM_ERROR; non-ACGT / < 2 nt → OLIGOTM_ERROR.
            Assert.That(PrimerDesigner.CalculateFractionBoundPrimer3(new string('A', 20) + new string('G', 17), 60), Is.NaN);
            Assert.That(PrimerDesigner.CalculateFractionBoundPrimer3("ACGTNACGTACGTACGTACG", 60), Is.NaN);
            Assert.That(PrimerDesigner.CalculateFractionBoundPrimer3("A", 60), Is.NaN);
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateFractionBoundPrimer3("AGCTAGCTAGCTAGCTAGCT", 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateFractionBoundPrimer3("AGCTAGCTAGCTAGCTAGCT", 100.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateFractionBoundPrimer3("AGCTAGCTAGCTAGCTAGCT", 60, dnaConcentrationNanomolar: 0));
            // The Tm is unchanged by the refactoring onto the shared seqtm core (primer3-py calc_tm 58.101 °C).
            Assert.That(PrimerDesigner.CalculateMeltingTemperaturePrimer3("AGCTAGCTAGCTAGCTAGCT"), Is.EqualTo(58.10101033826538).Within(1e-9));
        });
    }

    [Test]
    public void EvaluatePrimer_Bound_MatchesPrimer3CheckPrimers()
    {
        // check_primers AGCTAGCTAGCTAGCTAGCT (Tm 58.101): PRIMER_LEFT_0_PENALTY 1.8989896617346176 without an annealing
        // temperature (bound weights inert, no PRIMER_LEFT_0_BOUND), 14.460414507505782 with Ta 60 + WT_BOUND_LT 0.2 /
        // _GT 0.3 (opt 97), 6.303546879738778 with Ta 50 + opt 80; Ta 60 + MIN_BOUND 40 → "low fraction bound";
        // Ta 50 + MAX_BOUND 80 (opt 70) → "high fraction bound".
        const string s = "AGCTAGCTAGCTAGCTAGCT";
        var loose = P3 with { MinTm = 0, MaxTm = 100 };
        var weighted = WithWeights(loose, w => w with { BoundLt = 0.2, BoundGt = 0.3 });
        var none = PrimerDesigner.EvaluatePrimer(s, 0, true, weighted);
        var ta60 = PrimerDesigner.EvaluatePrimer(s, 0, true, weighted with { AnnealingTemperature = 60 });
        var ta50 = PrimerDesigner.EvaluatePrimer(s, 0, true, weighted with { AnnealingTemperature = 50, OptBound = 80 });
        var low = PrimerDesigner.EvaluatePrimer(s, 0, true, loose with { AnnealingTemperature = 60, MinBound = 40 });
        var high = PrimerDesigner.EvaluatePrimer(s, 0, true, loose with { AnnealingTemperature = 50, MaxBound = 80, OptBound = 70 });
        Assert.Multiple(() =>
        {
            Assert.That(none.Bound, Is.Null);
            Assert.That(none.Penalty, Is.EqualTo(1.8989896617346176).Within(1e-9));
            Assert.That(none.Issues, Has.None.Contains("BOUND"));
            Assert.That(ta60.Bound!.Value, Is.EqualTo(34.19287577114419).Within(1e-9));
            Assert.That(ta60.Penalty, Is.EqualTo(14.460414507505782).Within(1e-9));
            Assert.That(ta60.Issues, Has.None.Contains("BOUND"));
            Assert.That(ta50.Bound!.Value, Is.EqualTo(94.6818573933472).Within(1e-9));
            Assert.That(ta50.Penalty, Is.EqualTo(6.303546879738778).Within(1e-9));
            Assert.That(low.IsValid, Is.False);
            Assert.That(low.Issues, Has.Some.Contains("PRIMER_MIN_BOUND"));
            Assert.That(high.IsValid, Is.False);
            Assert.That(high.Issues, Has.Some.Contains("PRIMER_MAX_BOUND"));
            Assert.That(high.PositionPenalty, Is.Null); // no target in EvaluatePrimer
        });
    }

    [Test]
    public void CalculatePrimer3Penalty_BoundAndPositionTerms()
    {
        // p_obj_fn: bound_gt · (bound − opt) / bound_lt · (opt − bound) only with a bound value; pos_penalty · position.
        var w = PrimerDesigner.DefaultPrimer3Weights with { BoundGt = 0.3, BoundLt = 0.2 };
        var o = PrimerDesigner.DefaultPrimer3Optima with { OptBound = 90 };
        var inputs = new Primer3PenaltyInputs(60.0, 20, 50.0);
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.DefaultPrimer3Weights.PositionPenalty, Is.EqualTo(1.0)); // PRIMER_WT_POS_PENALTY default
            Assert.That(PrimerDesigner.DefaultPrimer3Optima.OptBound, Is.EqualTo(97.0));     // PRIMER_OPT_BOUND default
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs, w, o), Is.EqualTo(0.0));
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs with { Bound = 95 }, w, o), Is.EqualTo(0.3 * 5).Within(1e-12));
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs with { Bound = 40 }, w, o), Is.EqualTo(0.2 * 50).Within(1e-12));
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs with { PositionPenalty = 1.5 }, w with { PositionPenalty = 2 }, o),
                Is.EqualTo(3.0).Within(1e-12));
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs with { PositionPenalty = -20 }, w, o), Is.EqualTo(-20.0));
        });
    }

    #endregion

    #region Position penalty (compute_position_penalty)

    // Target [140, 160): left 3' end 139 → outside distance 0; 3' end 159 → inside 20; 3' end 160 → not allowed.
    [TestCase(108, 20, true, 2.0, 0.1, 1.2000000000000002)]   // PRIMER_LEFT [108,20]: 3' end 127, 12 bases before → 0.1 × 12
    [TestCase(117, 20, true, 2.0, 0.1, 0.30000000000000004)] // 3' end 136
    [TestCase(140, 20, true, -1.0, 0.5, -20.0)]               // 3' end 159 inside (20 × default −1)
    [TestCase(120, 20, true, 2.0, 0.1, 0.0)]                  // 3' end 139, adjacent to the target
    [TestCase(192, 20, false, 2.0, 0.1, 3.2)]                 // PRIMER_RIGHT [211,20]: 3' end 192, 32 bases after the target
    [TestCase(220, 20, false, -1.0, 0.5, 30.0)]               // PRIMER_RIGHT [239,20]
    [TestCase(150, 20, false, 2.0, 0.1, 20.0)]                // right 3' end 150 inside: (159 − 150 + 1) × 2
    public void CalculatePositionPenaltyPrimer3_MatchesPrimer3(int position, int length, bool forward, double inside,
        double outside, double expected)
    {
        Assert.That(PrimerDesigner.CalculatePositionPenaltyPrimer3(position, length, forward, TargetStart, TargetEnd, inside, outside)!.Value,
            Is.EqualTo(expected).Within(1e-12));
    }

    [Test]
    public void CalculatePositionPenaltyPrimer3_ThreePrimeEndBeyondTarget_IsInfinite()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePositionPenaltyPrimer3(141, 20, true, TargetStart, TargetEnd, 2, 0.1), Is.Null);
            Assert.That(PrimerDesigner.CalculatePositionPenaltyPrimer3(139, 20, false, TargetStart, TargetEnd, 2, 0.1), Is.Null);
            Assert.That(PrimerDesigner.CalculatePositionPenaltyPrimer3(140, 20, false, TargetStart, TargetEnd, 2, 0.1), Is.EqualTo(40.0));
        });
    }

    #endregion

    #region Pair designs

    [Test]
    public void DesignPrimerPairs_Defaults_NoBoundNoPositionValues()
    {
        AssertPairs(Design(P3, O3), DefaultRanks);
        // Bound weights without an annealing temperature are inert for primers (p_obj_fn gate).
        AssertPairs(Design(WithWeights(P3, w => w with { BoundLt = 0.1, BoundGt = 0.1 }), O3), DefaultRanks);
    }

    [Test]
    public void DesignPrimerPairs_AnnealingTemperature_ReportsBound()
    {
        // PRIMER_ANNEALING_TEMP 60: same pairs as the default, with PRIMER_LEFT/RIGHT_k_BOUND.
        AssertPairs(Design(P3 with { AnnealingTemperature = 60 }, O3),
        [
            (46, 20, 233, 20, 0.2229611791929642, 0.1161425661219937, 0.10681861307097051, 49.00904260079596, 50.93138362747912, null, null),
            (116, 20, 233, 20, 0.5002207880304468, 0.39340217495947627, 0.10681861307097051, 53.45996470116724, 50.93138362747912, null, null),
            (117, 20, 233, 20, 0.5003559657240544, 0.39353735265308387, 0.10681861307097051, 53.458931104716726, 50.93138362747912, null, null),
            (46, 20, 240, 20, 0.5299746133399026, 0.1161425661219937, 0.4138320472179089, 49.00904260079596, 46.49003072307379, null, null),
            (46, 20, 236, 20, 0.5739934657308936, 0.1161425661219937, 0.4578508996088999, 49.00904260079596, 45.93683824553874, null, null),
        ]);
    }

    [Test]
    public void DesignPrimerPairs_BoundLimitsAndWeights_MatchPrimer3()
    {
        // PRIMER_ANNEALING_TEMP 60, PRIMER_MIN_BOUND 50, PRIMER_OPT_BOUND 90, PRIMER_WT_BOUND_LT/GT 0.1.
        var p = WithWeights(P3 with { AnnealingTemperature = 60, MinBound = 50, OptBound = 90 }, w => w with { BoundLt = 0.1, BoundGt = 0.1 });
        AssertPairs(Design(p, O3),
        [
            (116, 20, 233, 20, 8.061085955165812, 4.0474057048427525, 4.013680250323059, 53.45996470116724, 50.93138362747912, null, null),
            (117, 20, 233, 20, 8.06132449250447, 4.047644242181411, 4.013680250323059, 53.458931104716726, 50.93138362747912, null, null),
            (108, 20, 233, 20, 8.093774633151943, 4.080094382828884, 4.013680250323059, 55.96372073905036, 50.93138362747912, null, null),
            (107, 20, 233, 20, 8.09419258021352, 4.08051232989046, 4.013680250323059, 55.96269694339334, 50.93138362747912, null, null),
            (116, 20, 247, 20, 8.133951250029043, 4.0474057048427525, 4.08654554518629, 53.45996470116724, 54.64585050032973, null, null),
        ]);
    }

    [Test]
    public void DesignPrimerPairs_MaxBoundAlignmentMode_MatchPrimer3()
    {
        // PRIMER_ANNEALING_TEMP 55, PRIMER_MAX_BOUND 80, PRIMER_OPT_BOUND 70, PRIMER_WT_BOUND_GT 0.2,
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0.
        var p = WithWeights(P3 with
        {
            AnnealingTemperature = 55, MaxBound = 80, OptBound = 70, StructureScreen = PrimerStructureScreen.Primer3Alignment,
        }, w => w with { BoundGt = 0.2 });
        AssertPairs(Design(p, O3),
        [
            (6, 20, 249, 20, 5.295322285016098, 2.6659860693652746, 2.6293362156508238, 70.1012795440591, 69.8784460302355, null, null),
            (71, 20, 249, 20, 5.312752404248988, 2.6834161885981644, 2.6293362156508238, 71.56716515091478, 69.8784460302355, null, null),
            (106, 20, 249, 20, 5.31333751684889, 2.684001301198066, 2.6293362156508238, 69.57381081947702, 69.8784460302355, null, null),
            (7, 20, 249, 20, 5.348663685959934, 2.7193274703091106, 2.6293362156508238, 72.12439002312152, 69.8784460302355, null, null),
            (6, 20, 231, 20, 5.369616776421802, 2.6659860693652746, 2.703630707056527, 70.1012795440591, 71.44318573021405, null, null),
        ]);
    }

    [Test]
    public void DesignPrimerPairs_InsideOutsidePenalty_MatchPrimer3()
    {
        // PRIMER_INSIDE_PENALTY 2, PRIMER_OUTSIDE_PENALTY 0.1: primers may now reach into the target.
        AssertPairs(Design(P3, O3 with { InsidePenalty = 2, OutsidePenalty = 0.1 }),
        [
            (108, 20, 211, 20, 5.6992044449417225, 1.8764664567339198, 3.8227379882078028, null, null, 1.2000000000000002, 3.2),
            (108, 20, 212, 20, 5.726416644075584, 1.8764664567339198, 3.849950187341665, null, null, 1.2000000000000002, 3.3000000000000003),
            (107, 20, 211, 20, 5.799520012437597, 1.9767820242297944, 3.8227379882078028, null, null, 1.3, 3.2),
            (107, 20, 212, 20, 5.82673221157146, 1.9767820242297944, 3.849950187341665, null, null, 1.3, 3.3000000000000003),
            (117, 20, 221, 21, 5.860395628099548, 0.6935373526530839, 5.166858275446464, null, null, 0.30000000000000004, 4.1000000000000005),
        ]);
    }

    [Test]
    public void DesignPrimerPairs_InsidePenaltyOnly_SamePairsWithZeroPositionPenalty()
    {
        // PRIMER_INSIDE_PENALTY 0.5 (outside stays 0): the default pairs, now with POSITION_PENALTY 0.0.
        var expected = DefaultRanks.Select(e => (e.Item1, e.Item2, e.Item3, e.Item4, e.Item5, e.Item6, e.Item7,
            (double?)null, (double?)null, (double?)0.0, (double?)0.0)).ToArray();
        AssertPairs(Design(P3, O3 with { InsidePenalty = 0.5 }), expected);
    }

    [Test]
    public void DesignPrimerPairs_OutsidePenaltyWithDefaultInside_NegativeInsidePenalty()
    {
        // PRIMER_OUTSIDE_PENALTY 0.5 only: the default inside penalty −1 multiplies inside distances as Primer3 does
        // (PRIMER_LEFT_0_POSITION_PENALTY −20, PRIMER_LEFT_0_PENALTY −19.16…).
        AssertPairs(Design(P3, O3 with { OutsidePenalty = 0.5 }),
        [
            (140, 20, 239, 20, 11.372928793929077, -19.16472395712634, 30.53765275105542, null, null, -20.0, 30.0),
            (139, 21, 238, 21, 11.69012839391496, -18.788787973004958, 30.47891636691992, null, null, -20.0, 29.0),
            (139, 21, 239, 20, 11.74886477805046, -18.788787973004958, 30.53765275105542, null, null, -20.0, 30.0),
            (140, 20, 240, 20, 11.749108090091568, -19.16472395712634, 30.91383204721791, null, null, -20.0, 30.5),
            (139, 21, 240, 20, 12.125044074212951, -18.788787973004958, 30.91383204721791, null, null, -20.0, 30.5),
        ]);
    }

    [Test]
    public void DesignPrimerPairs_NegativePairPenalty_ThrowsLikePrimer3Assert()
    {
        // PRIMER_OUTSIDE_PENALTY 0.05 only: a pair penalty < 0 — primer3-py aborts (obj_fn PR_ASSERT(sum >= 0.0)).
        Assert.Throws<InvalidOperationException>(() => Design(P3, O3 with { OutsidePenalty = 0.05 }));
    }

    [Test]
    public void DesignPrimerPairs_PositionWeightAndProductRange_MatchPrimer3()
    {
        // PRIMER_INSIDE_PENALTY 1, PRIMER_OUTSIDE_PENALTY 0.2, PRIMER_WT_POS_PENALTY 0.5, PRIMER_PRODUCT_SIZE_RANGE 60-150.
        var p = WithWeights(P3, w => w with { PositionPenalty = 0.5 });
        var o = O3 with { InsidePenalty = 1, OutsidePenalty = 0.2, ProductSizeRanges = [new ProductSizeRange(60, 150)] };
        AssertPairs(Design(p, o),
        [
            (117, 20, 199, 19, 4.479302227330118, 0.6935373526530839, 3.785764874677034, null, null, 0.6000000000000001, 4.2),
            (117, 20, 193, 18, 4.481483102278185, 0.6935373526530839, 3.7879457496251008, null, null, 0.6000000000000001, 3.2),
            (117, 20, 211, 20, 4.5162753408608864, 0.6935373526530839, 3.8227379882078028, null, null, 0.6000000000000001, 6.4),
            (117, 20, 212, 20, 4.543487539994749, 0.6935373526530839, 3.849950187341665, null, null, 0.6000000000000001, 6.6000000000000005),
            (116, 20, 199, 19, 4.57916704963651, 0.7934021749594763, 3.785764874677034, null, null, 0.8, 4.2),
        ]);
        // DesignPrimers returns rank 0.
        var best = PrimerDesigner.DesignPrimers(new DnaSequence(T), TargetStart, TargetEnd, p, o);
        Assert.That(best.PairPenalty!.Value, Is.EqualTo(4.479302227330118).Within(1e-9));
    }

    [Test]
    public void DesignPrimerPairs_InternalOligoBound_MatchPrimer3()
    {
        // PRIMER_PICK_INTERNAL_OLIGO 1, PRIMER_ANNEALING_TEMP 60 (global), PRIMER_INTERNAL_WT_BOUND_LT 0.1,
        // PRIMER_PAIR_WT_IO_PENALTY 1 → PRIMER_INTERNAL_k [152,20], PENALTY 4.7238266716664326, BOUND 51.561222614043864.
        var o = O3 with
        {
            PickInternalOligo = true,
            InternalOligo = new ProbeDesigner.Primer3ProbeSettings { WeightBoundLt = 0.1 },
            Weights = new Primer3PairWeights(InternalOligoPenalty: 1),
        };
        var pairs = Design(P3 with { AnnealingTemperature = 60 }, o);
        double[] pen = [4.946787850859397, 5.224047459696879, 5.224182637390487, 5.253801285006335, 5.297820137397326];
        Assert.That(pairs, Has.Count.EqualTo(5));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < 5; k++)
            {
                Assert.That(pairs[k].PairPenalty!.Value, Is.EqualTo(pen[k]).Within(1e-9), $"rank {k}");
                var io = pairs[k].InternalOligo!.Value;
                Assert.That((io.Start, io.Length), Is.EqualTo((152, 20)));
                Assert.That(io.Penalty, Is.EqualTo(4.7238266716664326).Within(1e-9));
                Assert.That(io.Bound!.Value, Is.EqualTo(51.561222614043864).Within(1e-9));
            }
            Assert.That(pairs[0].Forward!.Bound!.Value, Is.EqualTo(49.00904260079596).Within(1e-9));
        });
    }

    [Test]
    public void DesignPrimerPairs_InternalOligoBoundWithoutAnnealingTemperature_UsesOligoTmError()
    {
        // PRIMER_INTERNAL_WT_BOUND_LT 1e-6 without PRIMER_ANNEALING_TEMP: the internal-oligo p_obj_fn has no gate, so the
        // oligo's bound OLIGOTM_ERROR −999999.9999 adds 1e-6 × (97 + 999999.9999): PENALTY 1.1800459329708184, no BOUND.
        var o = O3 with { PickInternalOligo = true, InternalOligo = new ProbeDesigner.Primer3ProbeSettings { WeightBoundLt = 1e-6 } };
        var pairs = Design(P3, o);
        Assert.That(pairs, Has.Count.EqualTo(5));
        var io = pairs[0].InternalOligo!.Value;
        Assert.Multiple(() =>
        {
            Assert.That((io.Start, io.Length), Is.EqualTo((152, 20)));
            Assert.That(io.Penalty, Is.EqualTo(1.1800459329708184).Within(1e-9));
            Assert.That(io.Bound, Is.Null);
            Assert.That(pairs[0].PairPenalty!.Value, Is.EqualTo(0.2229611791929642).Within(1e-9));
        });
    }

    [Test]
    public void DesignPrimerPairs_InternalOligoMinBound_MatchPrimer3()
    {
        // PRIMER_ANNEALING_TEMP 60, PRIMER_INTERNAL_MIN_BOUND 60: internal oligos below 60 % bound are rejected.
        var o = O3 with { PickInternalOligo = true, InternalOligo = new ProbeDesigner.Primer3ProbeSettings { MinBound = 60 } };
        var pairs = Design(P3 with { AnnealingTemperature = 60 }, o);
        (int Start, double Pen, double Bound)[] expected =
        [
            (87, 1.1533159279721872, 60.29348152377712),
            (153, 1.1827027295453263, 60.24672823198022),
            (153, 1.1827027295453263, 60.24672823198022),
            (87, 1.1533159279721872, 60.29348152377712),
            (87, 1.1533159279721872, 60.29348152377712),
        ];
        Assert.That(pairs, Has.Count.EqualTo(5));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < 5; k++)
            {
                var io = pairs[k].InternalOligo!.Value;
                Assert.That(io.Start, Is.EqualTo(expected[k].Start), $"rank {k}");
                Assert.That(io.Length, Is.EqualTo(20), $"rank {k}");
                Assert.That(io.Penalty, Is.EqualTo(expected[k].Pen).Within(1e-9), $"rank {k}");
                Assert.That(io.Bound!.Value, Is.EqualTo(expected[k].Bound).Within(1e-9), $"rank {k}");
                Assert.That(pairs[k].PairPenalty!.Value, Is.EqualTo(DefaultRanks[k].Item5).Within(1e-9), $"rank {k}");
            }
        });
    }

    #endregion

    #region pick_hyb_probe_only

    [Test]
    public void DesignProbesPrimer3_AnnealingTemperature_ReportsBound()
    {
        var probes = ProbeDesigner.DesignProbesPrimer3(P, new ProbeDesigner.Primer3ProbeSettings { AnnealingTemperature = 60 });
        AssertProbes(probes,
        [
            (98, 22, 4.280897854548641, 29.15677603957223),
            (96, 24, 4.2993382742763515, 46.82527326747804),
            (97, 23, 4.534902697528139, 34.695810436472364),
            (96, 23, 5.163588669432841, 29.218841012366177),
            (95, 25, 5.286580402581137, 53.19609678990071),
        ]);
    }

    [Test]
    public void DesignProbesPrimer3_BoundLimitsAndWeights_MatchPrimer3()
    {
        var s = new ProbeDesigner.Primer3ProbeSettings
        {
            AnnealingTemperature = 60, MinBound = 50, MaxBound = 90, OptBound = 70, WeightBoundGt = 0.1, WeightBoundLt = 0.05,
        };
        AssertProbes(ProbeDesigner.DesignProbesPrimer3(P, s),
        [
            (95, 25, 6.126775563086102, 53.19609678990071),
            (94, 25, 6.126775563086102, 53.19609678990071),
            (94, 26, 7.994476279626942, 70.79890640386272),
        ]);
    }

    [Test]
    public void DesignProbesPrimer3_BoundWeightsWithoutAnnealingTemperature_MatchPrimer3()
    {
        // No annealing temperature: WT_BOUND_LT 1e-6 × (97 + 999999.9999) is added to every probe, WT_BOUND_GT is inert.
        var s = new ProbeDesigner.Primer3ProbeSettings { WeightBoundGt = 0.5, WeightBoundLt = 1e-6 };
        AssertProbes(ProbeDesigner.DesignProbesPrimer3(P, s),
        [
            (98, 22, 5.280994854448641, null),
            (96, 24, 5.299435274176352, null),
            (97, 23, 5.534999697428139, null),
            (96, 23, 6.163685669332841, null),
            (95, 25, 6.286677402481137, null),
        ]);
    }

    private static void AssertProbes(IReadOnlyList<ProbeDesigner.Primer3Probe> probes, (int Start, int Length, double Pen, double? Bound)[] expected)
    {
        Assert.That(probes, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                Assert.That((probes[k].Start, probes[k].Length), Is.EqualTo((expected[k].Start, expected[k].Length)), $"probe {k}");
                Assert.That(probes[k].Penalty, Is.EqualTo(expected[k].Pen).Within(1e-9), $"probe {k} PRIMER_INTERNAL_PENALTY");
                AssertNullable(probes[k].Bound, expected[k].Bound, $"probe {k} PRIMER_INTERNAL_BOUND");
            }
        });
    }

    #endregion

    #region Data control

    [Test]
    public void DataControl_MatchesPrimer3()
    {
        var dna = new DnaSequence(T);
        Assert.Multiple(() =>
        {
            // "Optimum primer fraction binding lower than minimum or higher than maximum".
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.DesignPrimerPairs(dna, 140, 160, P3 with { OptBound = 120 }, O3));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.DesignPrimerPairs(dna, 140, 160, P3 with { AnnealingTemperature = 60, MinBound = 99 }, O3));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.EvaluatePrimer("AGCTAGCTAGCTAGCTAGCT", 0, true, P3 with { OptBound = 120 }));
            // "Annealing temperature higher than 100 C".
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.DesignPrimerPairs(dna, 140, 160, P3 with { AnnealingTemperature = 101 }, O3));
            // "Optimum internal oligo fraction binding lower than minimum or higher than maximum".
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbesPrimer3(P, new ProbeDesigner.Primer3ProbeSettings { OptBound = -20 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.DesignPrimerPairs(dna, 140, 160, P3,
                O3 with { PickInternalOligo = true, InternalOligo = new ProbeDesigner.Primer3ProbeSettings { OptBound = -20 } }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbesPrimer3(P, new ProbeDesigner.Primer3ProbeSettings { AnnealingTemperature = 101 }));
        });
    }

    #endregion
}
