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
    public void DesignPrimerPairs_Primer3Defaults_MatchesPrimer3Ranks0To4()
    {
        // primer3.design_primers(SEQUENCE_TARGET=[200,30]) with every global setting at its default
        // (PRIMER_PRODUCT_SIZE_RANGE 100-300, PRIMER_NUM_RETURN 5, PRIMER_PAIR_MAX_DIFF_TM 100, sizes 18/20/27,
        // GC 20-80, poly-X 5). PRIMER_PAIR_EXPLAIN "considered 5, ok 5".
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T1), 200, 230,
            PrimerDesigner.Primer3DefaultParameters, PrimerPairOptions.Primer3Defaults);

        AssertPairs(pairs, new P3Pair[]
        {
            new(54, 20, 292, 20, 0.06818840804209003, 89.27954249148974, 19.17438697053791, 14.577563563433046, 239),
            new(166, 20, 292, 20, 0.06938078863225883, 85.68780526617428, 10.567995496614856, 0.0, 127),
            new(100, 20, 292, 20, 0.07344529090630658, 88.29365570067309, 0.0, 0.0, 193),
            new(172, 20, 292, 20, 0.12828968702643806, 85.1867315367541, 5.773113916838213, 19.035341488499114, 121),
            new(171, 20, 292, 20, 0.12828968702643806, 85.3856883176781, 10.567995496614856, 5.257952270280612, 122),
        });
    }

    [Test]
    public void DesignPrimers_ReturnsRankZeroOfDesignPrimerPairs()
    {
        var single = PrimerDesigner.DesignPrimers(new DnaSequence(T1), 200, 230,
            PrimerDesigner.Primer3DefaultParameters, PrimerPairOptions.Primer3Defaults);
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T1), 200, 230,
            PrimerDesigner.Primer3DefaultParameters, PrimerPairOptions.Primer3Defaults with { NumReturn = 1 });

        Assert.That(pairs, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(single.Forward!.Sequence, Is.EqualTo(pairs[0].Forward!.Sequence));
            Assert.That(single.Reverse!.Sequence, Is.EqualTo(pairs[0].Reverse!.Sequence));
            Assert.That(single.PairPenalty, Is.EqualTo(pairs[0].PairPenalty));
            Assert.That(single.ProductSize, Is.EqualTo(239));
            Assert.That(single.Message, Is.EqualTo("Valid primer pair found."));
        });
    }

    [Test]
    public void DesignPrimerPairs_NonDefaultPairWeights_MatchesPrimer3PairPenalty()
    {
        // PRIMER_PAIR_WT_PR_PENALTY 0.5, _DIFF_TM 1, _COMPL_ANY_TH 0.2, _COMPL_END_TH 1, PRIMER_PRODUCT_OPT_TM 82,
        // _PRODUCT_TM_LT 0.5, _PRODUCT_TM_GT 1, PRIMER_PRODUCT_OPT_SIZE 150, _PRODUCT_SIZE_LT 0.05, _GT 0.1
        // (obj_fn). PRIMER_PAIR_EXPLAIN "considered 854, ok 854".
        var options = PrimerPairOptions.Primer3Defaults with
        {
            ProductOptTm = 82.0,
            ProductOptSize = 150,
            Weights = new Primer3PairWeights(PrimerPenalty: 0.5, DiffTm: 1.0, ComplAnyTh: 0.2, ComplEndTh: 1.0,
                ProductTmLt: 0.5, ProductTmGt: 1.0, ProductSizeLt: 0.05, ProductSizeGt: 0.1),
        };
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T2), 170, 190,
            PrimerDesigner.Primer3DefaultParameters, options);

        AssertPairs(pairs, new P3Pair[]
        {
            new(110, 21, 257, 20, 1.4544162076632987, 81.85810745813002, 0.0, 0.0, 148),
            new(109, 22, 255, 20, 1.4974518405367778, 81.94925472784873, 0.0, 0.0, 147),
            new(111, 21, 257, 20, 1.5668555470359617, 81.94925472784873, 0.0, 0.0, 147),
            new(110, 21, 255, 20, 1.6259206260798702, 82.0416505903033, 0.0, 0.0, 146),
            new(111, 22, 255, 20, 1.6277635610157002, 82.13532087844689, 0.0, 0.0, 145),
        });
    }

    [Test]
    public void DesignPrimerPairs_SizeRangesInOrderAndProductTmLimits_MatchesPrimer3()
    {
        // PRIMER_PRODUCT_SIZE_RANGE "60-90 150-250" (ranges tried in order), PRIMER_PRODUCT_MIN_TM 80,
        // PRIMER_PRODUCT_MAX_TM 86, PRIMER_PAIR_MAX_DIFF_TM 2. PRIMER_PAIR_EXPLAIN "considered 297822,
        // unacceptable product size 297158, high product Tm 621, tm diff too large 38, ok 5": all five pairs
        // come from the first range.
        var options = PrimerPairOptions.Primer3Defaults with
        {
            ProductSizeRanges = [new ProductSizeRange(60, 90), new ProductSizeRange(150, 250)],
            ProductMinTm = 80.0,
            ProductMaxTm = 86.0,
            MaxTmDifference = 2.0,
        };
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T3), 240, 260,
            PrimerDesigner.Primer3DefaultParameters, options);

        AssertPairs(pairs, new P3Pair[]
        {
            new(216, 18, 298, 22, 6.683348115717195, 85.92933832755367, 0.0, 0.0, 83),
            new(215, 18, 296, 23, 6.828214973893068, 85.64238563874967, 0.0, 0.0, 82),
            new(214, 18, 296, 23, 6.828214973893068, 85.92933832755367, 0.0, 0.0, 83),
            new(215, 18, 297, 24, 6.901952550723195, 85.92933832755367, 0.0, 0.0, 83),
            new(222, 18, 299, 22, 7.760865031832907, 85.99792034794292, 0.0, 0.0, 78),
        });
    }

    [Test]
    public void DesignPrimerPairs_PickInternalOligo_MatchesPrimer3Triples()
    {
        // PRIMER_PICK_INTERNAL_OLIGO 1, PRIMER_PAIR_WT_IO_PENALTY 1 (PRIMER_INTERNAL_* defaults): the lowest-penalty
        // internal oligo strictly between the primers; pair penalty = primer penalties + internal penalty.
        // PRIMER_PAIR_EXPLAIN "considered 13, ok 13".
        var options = PrimerPairOptions.Primer3Defaults with
        {
            PickInternalOligo = true,
            Weights = new Primer3PairWeights(InternalOligoPenalty: 1.0),
        };
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T1), 200, 230,
            PrimerDesigner.Primer3DefaultParameters, options);

        AssertPairs(pairs, new P3Pair[]
        {
            new(54, 20, 292, 20, 0.09291595413179721, 89.27954249148974, 19.17438697053791, 14.577563563433046, 239, 128, 20, 0.024727546089707175, 59.97527245391029, 27.802919764544356, 16.35283407747835, 45.88049029595203 /* GCGAGCTATTCGGCCGGCAT */),
            new(100, 20, 292, 20, 0.09817283699601376, 88.29365570067309, 0.0, 0.0, 193, 128, 20, 0.024727546089707175, 59.97527245391029, 27.802919764544356, 16.35283407747835, 45.88049029595203 /* GCGAGCTATTCGGCCGGCAT */),
            new(54, 20, 315, 20, 0.16881632105111066, 89.32968781711497, 5.009602532675672, 0.0, 262, 128, 20, 0.024727546089707175, 59.97527245391029, 27.802919764544356, 16.35283407747835, 45.88049029595203 /* GCGAGCTATTCGGCCGGCAT */),
            new(100, 20, 315, 20, 0.1740732039153272, 88.45945880948139, 0.0, 0.0, 216, 128, 20, 0.024727546089707175, 59.97527245391029, 27.802919764544356, 16.35283407747835, 45.88049029595203 /* GCGAGCTATTCGGCCGGCAT */),
            new(166, 20, 315, 20, 0.17619779345295683, 86.32612547614804, 0.0, 0.0, 150, 258, 20, 0.030916637901384547, 59.969083362098615, 0.0, 0.0, 0.0 /* GCGGGAGAGGGTGGTACGGT */),
        });
        foreach (var p in pairs)
        {
            var io = p.InternalOligo!.Value;
            Assert.That(io.Start, Is.GreaterThan(p.Forward!.Position + p.Forward.Length - 1), "internal oligo after the left primer");
            Assert.That(io.Start + io.Length - 1, Is.LessThan(p.Reverse!.Position), "internal oligo before the right primer");
            Assert.That(io.Sequence, Is.EqualTo(T1.Substring(io.Start, io.Length)));
        }
    }

    [Test]
    public void DesignPrimerPairs_IncludedRegion_MatchesPrimer3()
    {
        // SEQUENCE_INCLUDED_REGION [150, 200]: PRIMER_PAIR_EXPLAIN "considered 7, ok 7".
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T1), 200, 230,
            PrimerDesigner.Primer3DefaultParameters, PrimerPairOptions.Primer3Defaults with { IncludedRegion = (150, 200) });

        var expected = new (int L, int R3, double Penalty)[]
        {
            (166, 292, 0.06938078863225883), (172, 292, 0.12828968702643806), (171, 292, 0.12828968702643806),
            (166, 315, 0.14528115555157228), (172, 315, 0.20419005394575152),
        };
        Assert.That(pairs, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                Assert.That(pairs[k].Forward!.Position, Is.EqualTo(expected[k].L));
                Assert.That(pairs[k].Reverse!.Position + pairs[k].Reverse!.Length - 1, Is.EqualTo(expected[k].R3));
                Assert.That(pairs[k].PairPenalty!.Value, Is.EqualTo(expected[k].Penalty).Within(1e-9));
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

    [Test]
    public void DesignPrimers_DefaultProductSizeRange_Is100To300()
    {
        // Library default options = Primer3 PRIMER_PRODUCT_SIZE_RANGE 100-300 (pr_set_default_global_args_1).
        Assert.That(PrimerPairOptions.Default.ProductSizeRanges, Is.EqualTo(new[] { new ProductSizeRange(100, 300) }));
        Assert.That(PrimerPairOptions.Default.MaxTmDifference, Is.EqualTo(PrimerDesigner.MaxPairTmDifference));
        Assert.That(PrimerPairOptions.Default.NumReturn, Is.EqualTo(5));
        Assert.That(PrimerPairOptions.Primer3Defaults.MaxTmDifference, Is.EqualTo(100.0));
        Assert.That(PrimerPairOptions.Default.Weights, Is.EqualTo(new Primer3PairWeights(1, 0, 0, 0, 0, 0, 0, 0, 0)));

        // A 48-bp template cannot hold a 100-bp product.
        var result = PrimerDesigner.DesignPrimers(new DnaSequence("TTGACCACAGCCAGGTTTAATTTTTTTTCAAATACGGTCACGCGCGGA"), 20, 28);
        Assert.That(result.IsValid, Is.False);
    }

    [Test]
    public void DesignPrimers_InvalidOptions_ThrowAsPrimer3DataControl()
    {
        var dna = new DnaSequence(T1);
        var p3 = PrimerDesigner.Primer3DefaultParameters;
        // primer3-py: "Product size is part of objective function while optimum size is not defined".
        Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimers(dna, 200, 230, p3,
            PrimerPairOptions.Primer3Defaults with { Weights = new Primer3PairWeights(ProductSizeLt: 1.0) }));
        Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimers(dna, 200, 230, p3,
            PrimerPairOptions.Primer3Defaults with { Weights = new Primer3PairWeights(ProductTmGt: 1.0) }));
        // "PRIMER_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE" (27 > 20).
        Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimers(dna, 200, 230, p3,
            PrimerPairOptions.Primer3Defaults with { ProductSizeRanges = [new ProductSizeRange(20, 300)] }));
        // "PRIMER_NUM_RETURN < 1".
        Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, 200, 230, p3,
            PrimerPairOptions.Primer3Defaults with { NumReturn = 0 }));
        // "TARGET outside of INCLUDED_REGION".
        Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimers(dna, 200, 230, p3,
            PrimerPairOptions.Primer3Defaults with { IncludedRegion = (0, 210) }));
        Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimers(dna, 200, 230, p3,
            PrimerPairOptions.Primer3Defaults with { ProductSizeRanges = [new ProductSizeRange(300, 100)] }));
        Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimers(dna, 200, 230, p3,
            PrimerPairOptions.Primer3Defaults with { MaxTmDifference = -1 }));
    }
}
