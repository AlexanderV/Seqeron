// PRIMER-DESIGN-001 / PROBE-DESIGN-001 — Primer3 sequence quality (audit round 3, A3-5 part 2a) and the
// PRIMER_PAIR_WT_IO_PENALTY consistency check (A3-23): SEQUENCE_QUALITY, PRIMER_MIN_QUALITY, PRIMER_MIN_END_QUALITY,
// PRIMER_QUALITY_RANGE_MIN/MAX, PRIMER_WT_SEQ_QUAL, PRIMER_WT_END_QUAL, PRIMER_INTERNAL_MIN_QUALITY,
// PRIMER_INTERNAL_WT_SEQ_QUAL, PRIMER_INTERNAL_WT_END_QUAL.
// Source: primer3-py 2.3.1 libprimer3.c sequence_quality_is_ok (seq_quality = min over the oligo, seq_end_quality = min
//         over the five 3'-most bases, both starting from PRIMER_QUALITY_RANGE_MAX; internal oligos: no end check),
//         calc_and_check_oligo_features (quality check after the end-GC check; OP_LOW_[END_]SEQUENCE_QUALITY are
//         five-prime problems), p_obj_fn (weights.seq_quality · (quality_range_max − seq_quality); weights.end_quality is
//         never read), choose_pair_or_triple (an empty internal-oligo list ends the search), _pr_data_control.
// Expected values: primer3-py 2.3.1 design_primers on the templates below with the quality arrays Q / Q2.
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_SequenceQuality_Tests
{
    // random.Random(5), 300 nt (same template as PrimerDesigner_BoundAndPosition_Tests).
    private const string T = "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG";
    private const int TargetStart = 140, TargetEnd = 160; // SEQUENCE_TARGET [140, 20]

    // random.Random(11), 120 nt (pick_hyb_probe_only).
    private const string P = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC";

    // SEQUENCE_QUALITY of T: q[i] = 12 if i % 61 == 5, 28 if i % 23 == 7, otherwise 35 + i % 6.
    private static readonly int[] Q = Enumerable.Range(0, T.Length)
        .Select(i => i % 61 == 5 ? 12 : i % 23 == 7 ? 28 : 35 + i % 6).ToArray();

    // SEQUENCE_QUALITY of P: q[i] = 12 if i % 37 == 4, 26 if i % 13 == 2, otherwise 34 + i % 7.
    private static readonly int[] Q2 = Enumerable.Range(0, P.Length)
        .Select(i => i % 37 == 4 ? 12 : i % 13 == 2 ? 26 : 34 + i % 7).ToArray();

    private static readonly PrimerParameters P3 = PrimerDesigner.Primer3DefaultParameters;
    private static readonly PrimerPairOptions O3 = PrimerPairOptions.Primer3Defaults with { SequenceQuality = Q };

    private static IReadOnlyList<PrimerPairResult> Design(PrimerParameters p, PrimerPairOptions o) =>
        PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), TargetStart, TargetEnd, p, o);

    private static PrimerParameters WithWeights(PrimerParameters p, Func<Primer3PenaltyWeights, Primer3PenaltyWeights> f) =>
        p with { PenaltyWeights = f(p.PenaltyWeights ?? PrimerDesigner.DefaultPrimer3Weights) };

    // (left start, left len, right 3'-end coordinate (PRIMER_RIGHT start), right len, PAIR_PENALTY, LEFT_PENALTY,
    //  RIGHT_PENALTY, LEFT_MIN_SEQ_QUALITY, RIGHT_MIN_SEQ_QUALITY).
    private static void AssertPairs(IReadOnlyList<PrimerPairResult> pairs,
        (int L, int LLen, int R, int RLen, double Pen, double LPen, double RPen, int LQ, int RQ)[] expected)
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
                Assert.That(p.Forward.MinSequenceQuality, Is.EqualTo(e.LQ), $"rank {k} PRIMER_LEFT_MIN_SEQ_QUALITY");
                Assert.That(p.Reverse.MinSequenceQuality, Is.EqualTo(e.RQ), $"rank {k} PRIMER_RIGHT_MIN_SEQ_QUALITY");
            }
        });
    }

    #region sequence_quality_is_ok

    [TestCase(46, 20, true, 28, 36)]   // left [46,20]: Q[53] = 28; 3' window 61..65 → 36
    [TestCase(217, 20, false, 35, 36)] // right PRIMER_RIGHT [236,20] (top strand 217..236): 3' window 217..221 → 36
    [TestCase(64, 3, true, 12, 12)]    // shorter than five bases: the whole oligo (Q[66] = 12)
    [TestCase(48, 10, false, 28, 35)]  // right 48..57: Q[53] = 28 lies outside the 3' window 48..52 (min 35)
    [TestCase(48, 10, true, 28, 28)]   // left 48..57: the 3' window 53..57 holds Q[53] = 28
    public void CalculateSequenceQualityPrimer3_MinAndEndMin(int position, int length, bool isForward, int min, int minEnd)
    {
        var (qMin, qEnd) = PrimerDesigner.CalculateSequenceQualityPrimer3(Q, position, length, isForward);
        Assert.Multiple(() =>
        {
            Assert.That(qMin, Is.EqualTo(min), "seq_quality");
            Assert.That(qEnd, Is.EqualTo(minEnd), "seq_end_quality");
        });
    }

    [Test]
    public void CalculateSequenceQualityPrimer3_StartsFromQualityRangeMax()
    {
        // q starts at PRIMER_QUALITY_RANGE_MAX: with range max 30 every minimum is capped at 30.
        var (qMin, qEnd) = PrimerDesigner.CalculateSequenceQualityPrimer3(Q, 217, 20, false, qualityRangeMax: 30);
        Assert.Multiple(() =>
        {
            Assert.That(qMin, Is.EqualTo(30));
            Assert.That(qEnd, Is.EqualTo(30));
        });
        Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateSequenceQualityPrimer3(Q, 290, 20, true));
    }

    #endregion

    #region Pair designs

    [Test]
    public void DesignPrimerPairs_MinQualityMinEndQualityAndWeight_MatchPrimer3()
    {
        // PRIMER_MIN_QUALITY 20, PRIMER_MIN_END_QUALITY 30, PRIMER_WT_SEQ_QUAL 0.1 (PRIMER_LEFT_EXPLAIN: low sequence quality 73).
        var p = WithWeights(P3 with { MinQuality = 20, MinEndQuality = 30 }, w => w with { SequenceQuality = 0.1 });
        AssertPairs(Design(p, O3),
        [
            (46, 20, 236, 20, 14.273993465730893, 7.316142566121994, 6.9578508996089, 28, 35),
            (46, 20, 212, 20, 14.366092753463658, 7.316142566121994, 7.049950187341665, 28, 35),
            (46, 20, 211, 20, 14.438880554329796, 7.316142566121994, 7.122737988207803, 28, 35),
            (46, 20, 235, 20, 14.500947172201233, 7.316142566121994, 7.18480460607924, 28, 35),
            (100, 21, 236, 20, 14.521319470098717, 7.5634685704898175, 6.9578508996089, 35, 35),
        ]);
        var single = PrimerDesigner.DesignPrimers(new DnaSequence(T), TargetStart, TargetEnd, p, O3);
        Assert.That(single.PairPenalty!.Value, Is.EqualTo(14.273993465730893).Within(1e-9), "DesignPrimers rank 0");
    }

    [Test]
    public void DesignPrimerPairs_EndQualityWeight_HasNoEffect()
    {
        // PRIMER_WT_END_QUAL 5: Primer3's p_obj_fn never reads weights.end_quality — the default ranking, with the
        // MIN_SEQ_QUALITY values reported because SEQUENCE_QUALITY is given.
        var p = WithWeights(P3, w => w with { EndQuality = 5.0 });
        AssertPairs(Design(p, O3),
        [
            (46, 20, 233, 20, 0.2229611791929642, 0.1161425661219937, 0.10681861307097051, 28, 28),
            (116, 20, 233, 20, 0.5002207880304468, 0.39340217495947627, 0.10681861307097051, 12, 28),
            (117, 20, 233, 20, 0.5003559657240544, 0.39353735265308387, 0.10681861307097051, 12, 28),
            (46, 20, 240, 20, 0.5299746133399026, 0.1161425661219937, 0.4138320472179089, 28, 28),
            (46, 20, 236, 20, 0.5739934657308936, 0.1161425661219937, 0.4578508996088999, 28, 35),
        ]);
    }

    [Test]
    public void DesignPrimerPairs_WithoutQuality_NoMinSequenceQuality()
    {
        var pairs = Design(P3, PrimerPairOptions.Primer3Defaults);
        Assert.That(pairs[0].Forward!.MinSequenceQuality, Is.Null);
    }

    [Test]
    public void DesignPrimerPairs_InternalOligoQuality_MatchPrimer3()
    {
        // PRIMER_PICK_INTERNAL_OLIGO 1, PRIMER_INTERNAL_MIN_QUALITY 25, PRIMER_INTERNAL_WT_SEQ_QUAL 0.2, PRIMER_PAIR_WT_IO_PENALTY 1,
        // PRIMER_QUALITY_RANGE_MAX 60, PRIMER_WT_SEQ_QUAL 0.05 → PRIMER_INTERNAL_k [31,20], PENALTY 5.3889380835424845,
        // MIN_SEQ_QUALITY 35 for every rank.
        var p = WithWeights(P3 with { QualityRangeMax = 60 }, w => w with { SequenceQuality = 0.05 });
        var o = O3 with
        {
            PickInternalOligo = true,
            InternalOligo = new ProbeDesigner.Primer3ProbeSettings { MinQuality = 25, WeightSequenceQuality = 0.2 },
            Weights = new Primer3PairWeights(InternalOligoPenalty: 1.0),
        };
        var pairs = Design(p, o);
        AssertPairs(pairs,
        [
            (11, 20, 233, 20, 9.380853713532826, 2.2850970169193716, 1.7068186130709706, 28, 28),
            (11, 20, 236, 20, 9.381886000070756, 2.2850970169193716, 1.7078508996089, 28, 35),
            (11, 20, 212, 20, 9.47398528780352, 2.2850970169193716, 1.7999501873416648, 28, 35),
            (11, 20, 211, 20, 9.546773088669658, 2.2850970169193716, 1.8727379882078026, 28, 35),
            (11, 20, 235, 20, 9.608839706541096, 2.2850970169193716, 1.9348046060792399, 28, 35),
        ]);
        Assert.Multiple(() =>
        {
            foreach (var pair in pairs)
            {
                var io = pair.InternalOligo!.Value;
                Assert.That((io.Start, io.Length), Is.EqualTo((31, 20)), "PRIMER_INTERNAL_k");
                Assert.That(io.Penalty, Is.EqualTo(5.3889380835424845).Within(1e-9), "PRIMER_INTERNAL_k_PENALTY");
                Assert.That(io.MinSequenceQuality, Is.EqualTo(35), "PRIMER_INTERNAL_k_MIN_SEQ_QUALITY");
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_NoAcceptableInternalOligo_NoPairs()
    {
        // PRIMER_INTERNAL_MIN_QUALITY 41: every internal oligo has low sequence quality (PRIMER_INTERNAL_EXPLAIN
        // "considered 287, GC content failed 4, low sequence quality 283, ok 0"), so Primer3 returns no pair.
        var o = O3 with { PickInternalOligo = true, InternalOligo = new ProbeDesigner.Primer3ProbeSettings { MinQuality = 41 } };
        Assert.That(Design(P3, o), Is.Empty);
        var single = PrimerDesigner.DesignPrimers(new DnaSequence(T), TargetStart, TargetEnd, P3, o);
        Assert.Multiple(() =>
        {
            Assert.That(single.IsValid, Is.False);
            Assert.That(single.Message, Does.Contain("internal oligo"));
        });
    }

    #endregion

    #region Hybridization probes

    [Test]
    public void DesignProbesPrimer3_Quality_MatchesPrimer3()
    {
        // pick_hyb_probe_only, PRIMER_INTERNAL_MIN_QUALITY 20, PRIMER_INTERNAL_WT_SEQ_QUAL 0.3 (low sequence quality 73).
        var s = new ProbeDesigner.Primer3ProbeSettings { MinQuality = 20, WeightSequenceQuality = 0.3 };
        var probes = ProbeDesigner.DesignProbesPrimer3(P, s, 5, Q2);
        (int Start, int Len, double Pen, int Q)[] expected =
        [
            (85, 25, 30.044504099519497, 26),
            (15, 26, 30.450610118501526, 26),
            (14, 26, 30.450610118501526, 26),
            (14, 27, 30.54957664505208, 26),
            (13, 27, 31.155945034107287, 26),
        ];
        Assert.That(probes, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                Assert.That((probes[k].Start, probes[k].Length), Is.EqualTo((expected[k].Start, expected[k].Len)), $"probe {k}");
                Assert.That(probes[k].Penalty, Is.EqualTo(expected[k].Pen).Within(1e-9), $"probe {k} PENALTY");
                Assert.That(probes[k].MinSequenceQuality, Is.EqualTo(expected[k].Q), $"probe {k} MIN_SEQ_QUALITY");
            }
        });
    }

    #endregion

    #region Penalty term

    [Test]
    public void CalculatePrimer3Penalty_SequenceQualityTerm()
    {
        var w = new Primer3PenaltyWeights(0, 0, 0, 0, 0, 0, 0, 0, 0) { SequenceQuality = 0.05, EndQuality = 3.0, PositionPenalty = 0 };
        var inputs = new Primer3PenaltyInputs(60, 20, 50) { SequenceQuality = 28, QualityRangeMax = 60 };
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs, w), Is.EqualTo(0.05 * 32).Within(1e-12));
            // No quality data: seq_quality = range max → no term.
            Assert.That(PrimerDesigner.CalculatePrimer3Penalty(inputs with { SequenceQuality = null }, w), Is.EqualTo(0.0));
        });
    }

    #endregion

    #region _pr_data_control

    [Test]
    public void DataControl_MatchesPrimer3Messages()
    {
        var dna = new DnaSequence(T);
        var none = PrimerPairOptions.Primer3Defaults;
        Assert.Multiple(() =>
        {
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, P3 with { MinQuality = 20 }, none),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Sequence quality data missing"));
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, P3, none with { SequenceQuality = Q[..^1] }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Error in sequence quality data"));
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, P3 with { QualityRangeMax = 39 }, O3),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Sequence quality score out of range"));
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, P3 with { MinQuality = 70, QualityRangeMax = 60 }, O3),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("PRIMER_MIN_QUALITY > PRIMER_QUALITY_RANGE_MAX"));
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, P3 with { QualityRangeMin = 12 },
                    O3 with { SequenceQuality = Q.Select(v => Math.Max(v, 12)).ToArray(), InternalOligo = new ProbeDesigner.Primer3ProbeSettings { MinQuality = 10 } }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("PRIMER_INTERNAL_MIN_QUALITY < PRIMER_QUALITY_RANGE_MIN"));
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, WithWeights(P3, w => w with { SequenceQuality = 1 }), none),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Sequence quality is part of objective function but sequence quality is not defined"));
            // The internal-oligo weight counts even without PRIMER_PICK_INTERNAL_OLIGO.
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, P3,
                    none with { InternalOligo = new ProbeDesigner.Primer3ProbeSettings { WeightSequenceQuality = 1 } }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Sequence quality is part of objective function but sequence quality is not defined"));
            // A3-23: PRIMER_PAIR_WT_IO_PENALTY without PRIMER_PICK_INTERNAL_OLIGO.
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, TargetStart, TargetEnd, P3,
                    none with { Weights = new Primer3PairWeights(InternalOligoPenalty: 0.5) }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Internal oligo quality is part of objective function while internal oligo choice is not required"));
            Assert.That(() => PrimerDesigner.DesignPrimers(dna, TargetStart, TargetEnd, P3,
                    none with { Weights = new Primer3PairWeights(InternalOligoPenalty: 0.5) }),
                NUnit.Framework.Throws.ArgumentException);
            // Probe picker and single-primer evaluation (no template quality).
            Assert.That(() => ProbeDesigner.DesignProbesPrimer3(P, new ProbeDesigner.Primer3ProbeSettings { MinQuality = 20 }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Sequence quality data missing"));
            Assert.That(() => ProbeDesigner.DesignProbesPrimer3(P, new ProbeDesigner.Primer3ProbeSettings { QualityRangeMax = 39 }, 5, Q2),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Sequence quality score out of range"));
            Assert.That(() => PrimerDesigner.EvaluatePrimer("TGCTTCAGAGTATGTATACC", 56, true, P3 with { MinQuality = 10 }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Sequence quality data missing"));
            Assert.That(() => PrimerDesigner.EvaluatePrimer("TGCTTCAGAGTATGTATACC", 56, true, WithWeights(P3, w => w with { SequenceQuality = 1 })),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith("Sequence quality is part of objective function"));
        });
    }

    #endregion
}
