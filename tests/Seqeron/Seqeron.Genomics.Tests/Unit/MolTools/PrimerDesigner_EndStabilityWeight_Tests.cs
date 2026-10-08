// PRIMER-DESIGN-001 / PRIMER-TM-001 — PRIMER_WT_END_STABILITY in EvaluatePrimer / DesignPrimers (audit round 3, A3-2).
// Source: primer3 libprimer3.c calc_and_check_oligo_features (h->end_stability = end_oligodg(seq, 5, santalucia), left/right
//         primers only) and p_obj_fn (sum += weights.end_stability · h->end_stability; the OT_INTL branch has no such term);
//         oligotm.c end_oligodg / oligodg.
// Expected values: primer3-py 2.3.1 design_primers (PRIMER_LEFT/RIGHT/PAIR/INTERNAL_k_*).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_EndStabilityWeight_Tests
{
    private const string T = "GGTCCTAATTGGAGCGCCCAGTTACCGGCCGAGTGCTACGGGCACTCGTTGGTAGTGGGCTCCCTAAGTCGGCGCATCCGTTCCTAGCTTTAAAATATCCGTTGAAAGAATGTTCTGAGTCTCGCCTAGTGAAAGCCAACTCCTTTGGATTTTGTCATA";
    private const string A1 = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";

    private static PrimerParameters WithEndStabilityWeight(double weight, PrimerStructureScreen screen = PrimerStructureScreen.Primer3Thermodynamic) =>
        PrimerDesigner.Primer3DefaultParameters with
        {
            StructureScreen = screen,
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { EndStability = weight },
        };

    [Test]
    public void EvaluatePrimer_EndStabilityWeight_AddsWeightTimesPrimer3EndStability()
    {
        // GTTACCGGCCGAGTGCTAC: PRIMER_LEFT_0_PENALTY 1.5179132593960958 (defaults), END_STABILITY 3.58 (TGCTAC),
        // 3.307913259396096 with PRIMER_WT_END_STABILITY 0.5.
        const string primer = "GTTACCGGCCGAGTGCTAC";
        var unweighted = PrimerDesigner.EvaluatePrimer(primer, 20, true, PrimerDesigner.Primer3DefaultParameters);
        var weighted = PrimerDesigner.EvaluatePrimer(primer, 20, true, WithEndStabilityWeight(0.5));
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.Calculate3PrimeStability(primer), Is.EqualTo(-3.58).Within(1e-12));
            Assert.That(unweighted.Penalty, Is.EqualTo(1.5179132593960958).Within(1e-9));
            Assert.That(weighted.Penalty, Is.EqualTo(3.307913259396096).Within(1e-9));
            Assert.That(weighted.IsValid, Is.EqualTo(unweighted.IsValid), "the weight changes the penalty only");
        });
    }

    private static readonly (int L, int LL, int R, int RL, double Pen, double LPen, double RPen)[] ThermoWt05 =
    {
        (20, 19, 137, 20, 5.234548087070011, 3.307913259396096, 1.9266348276739154),
        (19, 19, 137, 20, 5.351612109514906, 3.4249772818409907, 1.9266348276739154),
        (35, 19, 137, 20, 5.436791535391983, 3.5101567077180675, 1.9266348276739154), // rank 3 without the weight
        (20, 19, 136, 20, 5.460945351389999, 3.307913259396096, 2.1530320919939028),
        (19, 19, 136, 20, 5.5780093738348935, 3.4249772818409907, 2.1530320919939028),
    };

    private static readonly (int L, int LL, int R, int RL, double Pen, double LPen, double RPen)[] AlignmentWt1 =
    {
        (26, 20, 201, 20, 7.645601268359103, 5.107102963975174, 2.5384983043839293),
        (26, 20, 198, 20, 7.741262536698873, 5.107102963975174, 2.6341595727236995),
        (27, 20, 201, 20, 8.256828032536855, 5.718329728152926, 2.5384983043839293),
        (26, 20, 167, 20, 8.284611185810249, 5.107102963975174, 3.1775082218350756),
        (27, 20, 198, 20, 8.352489300876625, 5.718329728152926, 2.6341595727236995),
    };

    private static void AssertPairs(IReadOnlyList<PrimerPairResult> pairs,
        (int L, int LL, int R, int RL, double Pen, double LPen, double RPen)[] expected)
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
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_ThermodynamicMode_EndStabilityWeight_MatchesPrimer3()
    {
        // PRIMER_WT_END_STABILITY 0.5, SEQUENCE_TARGET [54, 21]; the weight re-ranks [35,19] from rank 3 to rank 2.
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, WithEndStabilityWeight(0.5), PrimerPairOptions.Primer3Defaults);
        AssertPairs(pairs, ThermoWt05);
    }

    [Test]
    public void DesignPrimerPairs_AlignmentMode_EndStabilityWeight_MatchesPrimer3()
    {
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_WT_END_STABILITY 1, SEQUENCE_TARGET [47, 17].
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(A1), 47, 64,
            WithEndStabilityWeight(1.0, PrimerStructureScreen.Primer3Alignment), PrimerPairOptions.Primer3Defaults);
        AssertPairs(pairs, AlignmentWt1);
    }

    [Test]
    public void DesignPrimers_EndStabilityWeight_ReturnsPrimer3RankZero()
    {
        var r = PrimerDesigner.DesignPrimers(new DnaSequence(T), 54, 75, WithEndStabilityWeight(0.5), PrimerPairOptions.Primer3Defaults);
        Assert.Multiple(() =>
        {
            Assert.That(r.IsValid, Is.True);
            Assert.That(r.PairPenalty!.Value, Is.EqualTo(5.234548087070011).Within(1e-9));
            Assert.That(r.Forward!.Penalty, Is.EqualTo(3.307913259396096).Within(1e-9));
        });
    }

    [Test]
    public void DesignPrimerPairs_PickInternalOligo_EndStabilityAppliesToPrimersOnly()
    {
        // PRIMER_WT_END_STABILITY 1 and PRIMER_INTERNAL_WT_END_STABILITY 1, PRIMER_PICK_INTERNAL_OLIGO 1: primer penalties
        // grow, PRIMER_INTERNAL_k_PENALTY stays 0.7844583918849821 (p_obj_fn OT_INTL has no end-stability term).
        var o = PrimerPairOptions.Primer3Defaults with { PickInternalOligo = true };
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, WithEndStabilityWeight(1.0), o);
        AssertPairs(pairs, new (int, int, int, int, double, double, double)[]
        {
            (20, 19, 137, 20, 8.70454808707001, 5.097913259396096, 3.6066348276739153),
            (35, 19, 137, 20, 8.741791535391982, 5.1351567077180675, 3.6066348276739153),
            (19, 19, 137, 20, 8.776612109514906, 5.169977281840991, 3.6066348276739153),
            (20, 19, 136, 20, 8.870945351389999, 5.097913259396096, 3.773032091993903),
            (35, 19, 136, 20, 8.90818879971197, 5.1351567077180675, 3.773032091993903),
        });
        Assert.Multiple(() =>
        {
            foreach (var p in pairs)
            {
                Assert.That(p.InternalOligo!.Value.Start, Is.EqualTo(54));
                Assert.That(p.InternalOligo!.Value.Length, Is.EqualTo(20));
                Assert.That(p.InternalOligo!.Value.Penalty, Is.EqualTo(0.7844583918849821).Within(1e-9));
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_ZeroEndStabilityWeight_UnchangedPrimer3Defaults()
    {
        // Primer3 default weight 0: ranks of primer3-py with only SEQUENCE_TEMPLATE / SEQUENCE_TARGET.
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, WithEndStabilityWeight(0.0), PrimerPairOptions.Primer3Defaults);
        AssertPairs(pairs, new (int, int, int, int, double, double, double)[]
        {
            (20, 19, 137, 20, 1.7645480870700112, 1.5179132593960958, 0.24663482767391542),
            (19, 19, 137, 20, 1.926612109514906, 1.6799772818409906, 0.24663482767391542),
            (20, 19, 136, 20, 2.0509453513899985, 1.5179132593960958, 0.5330320919939027),
            (35, 19, 137, 20, 2.131791535391983, 1.8851567077180675, 0.24663482767391542),
            (18, 20, 137, 20, 2.1875181726898063, 1.9408833450158909, 0.24663482767391542),
        });
    }
}
