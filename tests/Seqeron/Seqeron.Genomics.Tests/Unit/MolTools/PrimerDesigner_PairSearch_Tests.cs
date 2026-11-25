// PRIMER-DESIGN-001 — Primer3 pair search: PRIMER_PRODUCT_SIZE_RANGE, PRIMER_PAIR_MAX_DIFF_TM,
// PRIMER_NUM_RETURN, SEQUENCE_INCLUDED_REGION, pair objective PRIMER_PAIR_WT_*, product Tm,
// PRIMER_PICK_INTERNAL_OLIGO (audit round 1, L1/L2/L17).
// Source: primer3 libprimer3.cc (make_detection_primer_lists, pick_primer_range, choose_pair_or_triple,
//         characterize_pair, obj_fn, choose_internal_oligo, _pr_data_control, pr_set_default_global_args_1);
//         oligotm.c long_seq_tm.
// Expected values: primer3-py 2.3.1 design_primers (PRIMER_LEFT/RIGHT/INTERNAL_k, PRIMER_PAIR_k_PENALTY,
//         _PRODUCT_TM, _COMPL_ANY_TH, _COMPL_END_TH, _PRODUCT_SIZE), run with SEQUENCE_TARGET = [start, end − start].
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_PairSearch_Tests
{
    // Random templates (Python random.seed(20261001)).
    private const string T1 = "CGGTACCACACTTCGCGGTCGACCGTTCGTATTTAATTCTCGATTACTTCTGCCACCCAGCTAGGACTTTGTGCCGTATCCATTCACCTGCCCGGCCGCATTCCCCGTGACCATACCTCTAGGTGGTCGCGAGCTATTCGGCCGGCATCACGGTCCATACCGGTAGCGTACCCGGTCATCACAACATTCGCCAGATACAGCTTTCGGTTGAAGTAGTCAGGATAACTTCAGGCTAAGAGCAGTAGGTCGCATGCTCCAGCGGGAGAGGGTGGTACGGTCAGCTAGATGCGATGCAATCCTAAAATCGGCCAGCCTGTTGAGCTTGTACCAGCGACTATTCGCTCACCGGCGTGTAGTGATGCACTTCCAGGAGTGAGAGAGAAACCAGGCGCTCAATGGCACAGCTTGAGCCGCCGCTGT";
    private const string T2 = "GGTTATTCCGACGGTAAAACAAGAGTATGGACTGGAGACCTGAATATCTTTAGACTGACGAGCGTAGTCAGTACATGGTAACTTTATTGTTTGACGAAGAGTAAGCCTAATGTCGGTAGTTCTTCAGGTCGAGTCTAAGCTAACAGTTGTGTTGTCACACGAAGCTTAAATCTAAGATTGATTCACGTAAATTGAGCAGCTGGATTTATTAATTGACAGCTGACCTGTGAGGTACTCTGAGACCCAGAGGCTGTTGAATTCGCGCAAATCCCCTCCGAGTTCGTACTCTAAGTCCGTGGTAAGTCGGAAGTCGCAGTACTTGTAGTCAACAGAGATTCAGTATAGTTCACCTTCTCGCCT";
    private const string T3 = "GGCTTGCTAAGGCTACTTAGGTTAACTATTCTGGTGGCGCAGTGGCTCATCCTGAACACCGCCTGCATCCTGCCGGTACTATAATTCTCGTTTAACCCCAATACGTAAACAAGGCCCACGACGGATGATGCTAATGCTGGATATCGTGGACCTAAGATAAGATCAACCACGTGGGTCGCTTCCCACCTCTGCCAGCAGAAATACGCTCAGTCGGGCTTGGACGCACTCCCAGCCCGAGTCCCGGCGAAGGCGACTGGCAACCAGTCCGACTGTTAGTACAGGTGTGATAAGTCCGAAGCGGTCAAGCTGTAGTTGACTAGCAAAGGAAGGGGTCCTCGGTAACCTCGATGTCTGACGTATCGGCCACGCCTACCGACGGATGAGAAGTAGTGGTTAGGGCCAAGCTGCCCGCGAAGGTCTCATGGTTTTGTATTGATCCCGGTAGTCTATCCTGGCGAGTGCTGAATGTCTCGTAACCTCCACCTCTTGTGGACGGAATA";

    // One primer3-py pair row: PRIMER_LEFT_k (start, length), PRIMER_RIGHT_k (3′ position, length),
    // PRIMER_PAIR_k_PENALTY / _PRODUCT_TM / _COMPL_ANY_TH / _COMPL_END_TH / _PRODUCT_SIZE and, when picked,
    // PRIMER_INTERNAL_k (start, length) / _PENALTY / _TM / _SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH.
    private sealed record P3Pair(
        int LeftStart, int LeftLength, int Right3Prime, int RightLength, double Penalty, double ProductTm,
        double ComplAnyTh, double ComplEndTh, int ProductSize,
        int IntStart = -1, int IntLength = 0, double IntPenalty = 0, double IntTm = 0,
        double IntSelfAnyTh = 0, double IntSelfEndTh = 0, double IntHairpinTh = 0);

    private static void AssertPairs(IReadOnlyList<PrimerPairResult> actual, P3Pair[] expected)
    {
        Assert.That(actual, Has.Count.EqualTo(expected.Length), "PRIMER_PAIR_NUM_RETURNED");
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var a = actual[k];
                var e = expected[k];
                Assert.That(a.IsValid, Is.True, $"rank {k}");
                Assert.That(a.Forward!.Position, Is.EqualTo(e.LeftStart), $"rank {k} PRIMER_LEFT start");
                Assert.That(a.Forward.Length, Is.EqualTo(e.LeftLength), $"rank {k} PRIMER_LEFT length");
                Assert.That(a.Reverse!.Position + a.Reverse.Length - 1, Is.EqualTo(e.Right3Prime), $"rank {k} PRIMER_RIGHT position");
                Assert.That(a.Reverse.Length, Is.EqualTo(e.RightLength), $"rank {k} PRIMER_RIGHT length");
                Assert.That(a.PairPenalty!.Value, Is.EqualTo(e.Penalty).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(a.ProductTm!.Value, Is.EqualTo(e.ProductTm).Within(1e-9), $"rank {k} PRIMER_PAIR_PRODUCT_TM");
                Assert.That(a.ComplAnyTh!.Value, Is.EqualTo(e.ComplAnyTh).Within(1e-9), $"rank {k} PRIMER_PAIR_COMPL_ANY_TH");
                Assert.That(a.ComplEndTh!.Value, Is.EqualTo(e.ComplEndTh).Within(1e-9), $"rank {k} PRIMER_PAIR_COMPL_END_TH");
                Assert.That(a.ProductSize, Is.EqualTo(e.ProductSize), $"rank {k} PRIMER_PAIR_PRODUCT_SIZE");
                if (e.IntStart < 0)
                {
                    Assert.That(a.InternalOligo, Is.Null);
                    continue;
                }
                var io = a.InternalOligo!.Value;
                Assert.That(io.Start, Is.EqualTo(e.IntStart), $"rank {k} PRIMER_INTERNAL start");
                Assert.That(io.Length, Is.EqualTo(e.IntLength), $"rank {k} PRIMER_INTERNAL length");
                Assert.That(io.Penalty, Is.EqualTo(e.IntPenalty).Within(1e-9), $"rank {k} PRIMER_INTERNAL_PENALTY");
                Assert.That(io.Tm, Is.EqualTo(e.IntTm).Within(1e-9), $"rank {k} PRIMER_INTERNAL_TM");
                Assert.That(io.SelfAnyTh, Is.EqualTo(e.IntSelfAnyTh).Within(1e-9), $"rank {k} PRIMER_INTERNAL_SELF_ANY_TH");
                Assert.That(io.SelfEndTh, Is.EqualTo(e.IntSelfEndTh).Within(1e-9), $"rank {k} PRIMER_INTERNAL_SELF_END_TH");
                Assert.That(io.HairpinTh, Is.EqualTo(e.IntHairpinTh).Within(1e-9), $"rank {k} PRIMER_INTERNAL_HAIRPIN_TH");
            }
        });
    }

    [Test]
    public void CalculateProductMeltingTemperaturePrimer3_MatchesPrimer3ProductTm()
    {
        // PRIMER_PAIR_0_PRODUCT_TM of the Primer3-defaults case: product T1[54..292], 239 bp (long_seq_tm at
        // 50 mM monovalent + 120·√(1.5 − 0.6) mM, oligotm.c).
        Assert.That(PrimerDesigner.CalculateProductMeltingTemperaturePrimer3(T1.Substring(54, 239)),
            Is.EqualTo(89.27954249148974).Within(1e-9));
        Assert.That(PrimerDesigner.CalculateProductMeltingTemperaturePrimer3(T1.Substring(54, 239).ToLowerInvariant()),
            Is.EqualTo(89.27954249148974).Within(1e-9));
        Assert.Throws<ArgumentException>(() => PrimerDesigner.CalculateProductMeltingTemperaturePrimer3(""));
    }
}
