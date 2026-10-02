// PRIMER-STRUCT-001 / PRIMER-DESIGN-001 — Primer3 alignment mode (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0):
// dpal LOCAL self_any / compl_any, the Primer3Alignment structure screen (PRIMER_MAX_SELF_ANY/_END,
// PRIMER_PAIR_MAX_COMPL_ANY/_END, PRIMER_INTERNAL_MAX_SELF_ANY/_END), PRIMER_WT_SELF_ANY/_END,
// PRIMER_PAIR_WT_COMPL_ANY/_END, and per-primer PRIMER_WT_*_TH weights (audit round 2, A1).
// Source: primer3 dpal.c (_dpal_long_nopath_maxgap1_local, set_dpal_args), libprimer3.cc (align, oligo_compl,
//         characterize_pair, choose_internal_oligo, calc_and_check_oligo_features, p_obj_fn, obj_fn,
//         pr_set_default_global_args_1, _pr_data_control), raw.githubusercontent.com/primer3-org/primer3.
// Expected values: dpal.c compiled with Primer3's align() wrapper; primer3-py 2.3.1 design_primers.
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_AlignmentMode_Tests
{
    // (sequence 1, sequence 2, self_any(1), self_end(1), compl_any(1, 2)) from dpal.c + align(). self_end of |X| ≤ 3
    // is dpal.c _dpal_generic (force_generic; the fast GLOBAL_END routine reads past the string there).
    private static readonly object[] DpalCases =
    {
        new object[] { "AC", "GT", 2.0, 0.0, 2.0 },                     // second sequence < 3 nt: align() returns its length
        new object[] { "A", "T", 1.0, 0.0, 1.0 },
        new object[] { "ACG", "CGT", 2.0, 2.0, 3.0 },
        new object[] { "ACGTACGTACGTACGTACGT", "AAAAAAAAAAAAAAAAAAAA", 20.0, 20.0, 1.0 },
        new object[] { "GAATTCGAATTCGAATTC", "CCGGAATTCCAGTCAGT", 18.0, 18.0, 6.0 },
        new object[] { "AGCTAGCTNNAGCTAGCTAG", "CTAGCTNGCTAGCT", 15.5, 9.0, 8.75 }, // N scores −0.25
        new object[] { "TTGACAGCTAGCTCAGTCCTAGG", "CCTAGGACTGAGCTAGCTGTCAA", 8.0, 6.0, 23.0 },
        new object[] { "CGCGCGATATCGCGCG", "ATGCATGCATGCAT", 16.0, 16.0, 2.0 },
        new object[]
        {
            "ACGTTGCAAGGCTTAACCGGTACGATCGATGCATGCTAGCTAGCTGATCGATCGTAGCTAGCTAGCTGATCGATGCATGCATCGATGCTAGCTAGCTAGCTAGCTGATCGATCGATGCATGCATGCAT",
            "GGGCCCAAATTTGGGCCCAAATTTACGTACGTACGT", 58.0, 58.0, 5.0,
        },
    };

    [TestCaseSource(nameof(DpalCases))]
    public void DpalScores_MatchCompiledDpal(string s1, string s2, double selfAny, double selfEnd, double complAny)
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerSelfAnyComplementarity(s1), Is.EqualTo(selfAny), "self_any");
            Assert.That(PrimerDesigner.CalculatePrimerSelfEndComplementarity(s1), Is.EqualTo(selfEnd), "self_end");
            Assert.That(PrimerDesigner.CalculatePrimerDimerAnyComplementarity(s1, s2), Is.EqualTo(complAny), "compl_any");
            Assert.That(PrimerDesigner.CalculatePrimerSelfAnyComplementarity(s1.ToLowerInvariant()), Is.EqualTo(selfAny), "case-insensitive");
        });
    }

    [Test]
    public void DpalScores_NullOrEmpty_ReturnZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerSelfAnyComplementarity(""), Is.EqualTo(0.0));
            Assert.That(PrimerDesigner.CalculatePrimerSelfAnyComplementarity(null!), Is.EqualTo(0.0));
            Assert.That(PrimerDesigner.CalculatePrimerDimerAnyComplementarity("ACGT", ""), Is.EqualTo(0.0));
        });
    }

    [Test]
    public void EvaluatePrimer_Primer3AlignmentScreen_AppliesSelfAnyAndSelfEndLimits()
    {
        // GAATTC×3: self_any = self_end = 18.00 > PRIMER_MAX_SELF_ANY 8 / PRIMER_MAX_SELF_END 3.
        var p = PrimerDesigner.Primer3DefaultParameters with { StructureScreen = PrimerStructureScreen.Primer3Alignment, MinTm = 0 };
        var c = PrimerDesigner.EvaluatePrimer("GAATTCGAATTCGAATTC", 0, true, p);
        var relaxed = PrimerDesigner.EvaluatePrimer("GAATTCGAATTCGAATTC", 0, true, p with { MaxSelfAny = 18, MaxSelfEnd = 18 });

        Assert.Multiple(() =>
        {
            Assert.That(c.SelfAny, Is.EqualTo(18.0));
            Assert.That(c.SelfEnd, Is.EqualTo(18.0));
            Assert.That(c.SelfAnyTh, Is.Null);
            Assert.That(c.HairpinTh, Is.Null);
            Assert.That(c.HasHairpin, Is.False, "alignment mode has no hairpin value");
            Assert.That(c.Issues, Has.Some.EqualTo("Self-complementarity 18.00 exceeds 8.00 (Primer3 PRIMER_MAX_SELF_ANY)"));
            Assert.That(c.Issues, Has.Some.EqualTo("3' self-complementarity 18.00 exceeds 3.00 (Primer3 PRIMER_MAX_SELF_END)"));
            Assert.That(relaxed.Issues, Has.None.Contain("self-complementarity").IgnoreCase);
        });
    }

    // One primer3-py pair row (alignment mode).
    private sealed record P3AlnPair(
        int LeftStart, int LeftLength, int Right3Prime, int RightLength, double Penalty, double ComplAny, double ComplEnd,
        double LeftSelfAny, double LeftSelfEnd, double RightSelfAny, double RightSelfEnd, double LeftPenalty, double RightPenalty,
        double ProductTm, int IntStart = -1, int IntLength = 0, double IntPenalty = 0, double IntSelfAny = 0, double IntSelfEnd = 0);

    private static void AssertPairs(IReadOnlyList<PrimerPairResult> actual, P3AlnPair[] expected)
    {
        Assert.That(actual, Has.Count.EqualTo(expected.Length), "PRIMER_PAIR_NUM_RETURNED");
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var a = actual[k];
                var e = expected[k];
                Assert.That(a.Forward!.Position, Is.EqualTo(e.LeftStart), $"rank {k} PRIMER_LEFT start");
                Assert.That(a.Forward.Length, Is.EqualTo(e.LeftLength), $"rank {k} PRIMER_LEFT length");
                Assert.That(a.Reverse!.Position + a.Reverse.Length - 1, Is.EqualTo(e.Right3Prime), $"rank {k} PRIMER_RIGHT");
                Assert.That(a.Reverse.Length, Is.EqualTo(e.RightLength), $"rank {k} PRIMER_RIGHT length");
                Assert.That(a.PairPenalty!.Value, Is.EqualTo(e.Penalty).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(a.ComplAny, Is.EqualTo(e.ComplAny), $"rank {k} PRIMER_PAIR_COMPL_ANY");
                Assert.That(a.ComplEnd, Is.EqualTo(e.ComplEnd), $"rank {k} PRIMER_PAIR_COMPL_END");
                Assert.That(a.ComplAnyTh, Is.Null, $"rank {k} no _TH value in alignment mode");
                Assert.That(a.Forward.SelfAny, Is.EqualTo(e.LeftSelfAny), $"rank {k} PRIMER_LEFT_SELF_ANY");
                Assert.That(a.Forward.SelfEnd, Is.EqualTo(e.LeftSelfEnd), $"rank {k} PRIMER_LEFT_SELF_END");
                Assert.That(a.Reverse.SelfAny, Is.EqualTo(e.RightSelfAny), $"rank {k} PRIMER_RIGHT_SELF_ANY");
                Assert.That(a.Reverse.SelfEnd, Is.EqualTo(e.RightSelfEnd), $"rank {k} PRIMER_RIGHT_SELF_END");
                Assert.That(a.Forward.Penalty, Is.EqualTo(e.LeftPenalty).Within(1e-9), $"rank {k} PRIMER_LEFT_PENALTY");
                Assert.That(a.Reverse.Penalty, Is.EqualTo(e.RightPenalty).Within(1e-9), $"rank {k} PRIMER_RIGHT_PENALTY");
                Assert.That(a.ProductTm!.Value, Is.EqualTo(e.ProductTm).Within(1e-9), $"rank {k} PRIMER_PAIR_PRODUCT_TM");
                if (e.IntStart < 0)
                {
                    Assert.That(a.InternalOligo, Is.Null);
                    continue;
                }
                var io = a.InternalOligo!.Value;
                Assert.That(io.Start, Is.EqualTo(e.IntStart), $"rank {k} PRIMER_INTERNAL start");
                Assert.That(io.Length, Is.EqualTo(e.IntLength), $"rank {k} PRIMER_INTERNAL length");
                Assert.That(io.Penalty, Is.EqualTo(e.IntPenalty).Within(1e-9), $"rank {k} PRIMER_INTERNAL_PENALTY");
                Assert.That(io.SelfAny, Is.EqualTo(e.IntSelfAny), $"rank {k} PRIMER_INTERNAL_SELF_ANY");
                Assert.That(io.SelfEnd, Is.EqualTo(e.IntSelfEnd), $"rank {k} PRIMER_INTERNAL_SELF_END");
                Assert.That(io.HairpinTh, Is.NaN, $"rank {k} no hairpin in alignment mode");
            }
        });
    }

    private const string A1 = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";
    private const string A2 = "GTAGAATTTCCCAAAAATTCTAAGAACCCAATAGCATTCCTCTGACTTTCTCGCAGCCTGCGAGAAACGATATGATGGCTTGTCCTGGTACTATTTATTGGCCCCTCCAATAAATGATACTAAAGGGTCGATTCTAAGAGTCAAGTTATCCGCGGTTTGACGCGGCCCCTCTGCCA";

    [Test]
    public void DesignPrimerPairs_Primer3AlignmentDefaults_MatchesPrimer3()
    {
        // primer3.design_primers(SEQUENCE_TARGET = [47, 17], PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0), all else default.
        var p = PrimerDesigner.Primer3DefaultParameters with { StructureScreen = PrimerStructureScreen.Primer3Alignment };
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(A1), 47, 64, p, PrimerPairOptions.Primer3Defaults);

        AssertPairs(pairs, new P3AlnPair[]
        {
            new(26, 20, 172, 20, 0.353112963359024, 4, 2, 5, 3, 4, 2, 0.3171029639751737, 0.0360099993838503, 85.01728193873308),
            new(26, 20, 167, 20, 0.4246111858102495, 4, 1, 5, 3, 4, 2, 0.3171029639751737, 0.1075082218350758, 85.0228390911715),
            new(26, 20, 180, 20, 0.42552744112833807, 5, 2, 5, 3, 4, 3, 0.3171029639751737, 0.10842447715316439, 85.22074913206201),
            new(26, 20, 181, 20, 0.42552744112833807, 5, 2, 5, 3, 4, 2, 0.3171029639751737, 0.10842447715316439, 85.11330496332754),
            new(26, 20, 201, 20, 0.565601268359103, 5, 0, 5, 3, 4, 2, 0.3171029639751737, 0.24849830438392928, 85.78332244584502),
        });
    }

    [Test]
    public void DesignPrimerPairs_Primer3AlignmentNonDefaultLimitsWeightsAndInternalOligo_MatchesPrimer3()
    {
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_MAX_SELF_END 1, PRIMER_PAIR_MAX_COMPL_ANY 5,
        // PRIMER_WT_SELF_ANY 1, PRIMER_PAIR_WT_COMPL_ANY 0.1, PRIMER_PAIR_WT_COMPL_END 1,
        // PRIMER_PICK_INTERNAL_OLIGO 1, PRIMER_INTERNAL_MAX_SELF_ANY 4; SEQUENCE_TARGET [54, 9].
        var p = PrimerDesigner.Primer3DefaultParameters with
        {
            StructureScreen = PrimerStructureScreen.Primer3Alignment,
            MaxSelfEnd = 1.0,
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { SelfAny = 1.0 },
        };
        var o = PrimerPairOptions.Primer3Defaults with
        {
            MaxComplAny = 5.0,
            Weights = new Primer3PairWeights { ComplAny = 0.1, ComplEnd = 1.0 },
            PickInternalOligo = true,
            InternalOligo = new ProbeDesigner.Primer3ProbeSettings { MaxSelfAny = 4.0 },
        };
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(A2), 54, 63, p, o);

        AssertPairs(pairs, new P3AlnPair[]
        {
            new(25, 22, 161, 21, 14.72040871267958, 4, 0, 4, 1, 6, 1, 7.118531412319669, 7.201877300359911, 82.63464129123321, 55, 24, 5.426594526493034, 4, 2),
            new(25, 22, 155, 21, 14.967700954757323, 3, 0, 4, 1, 6, 0, 7.118531412319669, 7.549169542437653, 82.34495499268748, 55, 24, 5.426594526493034, 4, 2),
            new(25, 23, 155, 20, 15.008747465067689, 3, 0, 4, 0, 6, 1, 7.525198807594052, 7.1835486574736365, 82.34495499268748, 55, 24, 5.426594526493034, 4, 2),
            new(25, 21, 155, 20, 15.051026694640246, 3, 0, 4, 1, 6, 1, 7.567478037166609, 7.1835486574736365, 82.34495499268748, 55, 24, 5.426594526493034, 4, 2),
            new(25, 23, 161, 21, 15.127076107953963, 4, 0, 4, 0, 6, 1, 7.525198807594052, 7.201877300359911, 82.63464129123321, 55, 24, 5.426594526493034, 4, 2),
        });
    }

    [Test]
    public void DesignPrimerPairs_ThermodynamicPerPrimerWeights_MatchesPrimer3()
    {
        // Default thermodynamic mode with PRIMER_WT_SELF_ANY_TH 0.5, PRIMER_WT_SELF_END_TH 0.5, PRIMER_WT_HAIRPIN_TH 0.1
        // (p_obj_fn temp_cutoff terms over the ntthal values); SEQUENCE_TARGET [54, 21].
        const string t = "GGTCCTAATTGGAGCGCCCAGTTACCGGCCGAGTGCTACGGGCACTCGTTGGTAGTGGGCTCCCTAAGTCGGCGCATCCGTTCCTAGCTTTAAAATATCCGTTGAAAGAATGTTCTGAGTCTCGCCTAGTGAAAGCCAACTCCTTTGGATTTTGTCATA";
        var p = PrimerDesigner.Primer3DefaultParameters with
        {
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { SelfAnyTh = 0.5, SelfEndTh = 0.5, HairpinTh = 0.1 },
        };
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(t), 54, 75, p, PrimerPairOptions.Primer3Defaults);

        var expected = new (int L, int LL, int R, int RL, double Pen, double LPen, double RPen)[]
        {
            (20, 19, 137, 20, 1.8289866358036504, 1.556416885012145, 0.2725697507915054),
            (19, 19, 137, 20, 1.9888587152777533, 1.7162889644862478, 0.2725697507915054),
            (20, 19, 136, 20, 2.1093257363245557, 1.556416885012145, 0.5529088513124109),
            (35, 19, 137, 20, 2.1845454860826283, 1.911975735291123, 0.2725697507915054),
            (18, 20, 137, 20, 2.248120771188721, 1.9755510203972155, 0.2725697507915054),
        };
        Assert.That(pairs, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var e = expected[k];
                Assert.That(pairs[k].Forward!.Position, Is.EqualTo(e.L), $"rank {k}");
                Assert.That(pairs[k].Forward!.Length, Is.EqualTo(e.LL), $"rank {k}");
                Assert.That(pairs[k].Reverse!.Position + pairs[k].Reverse!.Length - 1, Is.EqualTo(e.R), $"rank {k}");
                Assert.That(pairs[k].Reverse!.Length, Is.EqualTo(e.RL), $"rank {k}");
                Assert.That(pairs[k].PairPenalty!.Value, Is.EqualTo(e.Pen).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(pairs[k].Forward!.Penalty, Is.EqualTo(e.LPen).Within(1e-9), $"rank {k} PRIMER_LEFT_PENALTY");
                Assert.That(pairs[k].Reverse!.Penalty, Is.EqualTo(e.RPen).Within(1e-9), $"rank {k} PRIMER_RIGHT_PENALTY");
            }
        });
    }

    [Test]
    public void DesignPrimers_Primer3AlignmentScreen_ReturnsRankZero()
    {
        var p = PrimerDesigner.Primer3DefaultParameters with { StructureScreen = PrimerStructureScreen.Primer3Alignment };
        var r = PrimerDesigner.DesignPrimers(new DnaSequence(A1), 47, 64, p, PrimerPairOptions.Primer3Defaults);
        Assert.Multiple(() =>
        {
            Assert.That(r.IsValid, Is.True);
            Assert.That(r.PairPenalty!.Value, Is.EqualTo(0.353112963359024).Within(1e-9));
            Assert.That(r.ComplAny, Is.EqualTo(4.0));
            Assert.That(r.ComplEnd, Is.EqualTo(2.0));
        });
    }

    [Test]
    public void DesignPrimerPairs_IllegalComplementarityLimits_Throw()
    {
        // _pr_data_control: complementarity limits must lie in [0, SHRT_MAX].
        var dna = new DnaSequence(A1);
        var p = PrimerDesigner.Primer3DefaultParameters with { StructureScreen = PrimerStructureScreen.Primer3Alignment };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, 47, 64, p,
                PrimerPairOptions.Primer3Defaults with { MaxComplAny = -1 }));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, 47, 64, p,
                PrimerPairOptions.Primer3Defaults with { MaxComplEnd = 40000 }));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, 47, 64, p with { MaxSelfAny = -0.5 },
                PrimerPairOptions.Primer3Defaults));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, 47, 64, p,
                PrimerPairOptions.Primer3Defaults with { Weights = new Primer3PairWeights { ComplEnd = -1 } }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.DesignPrimerPairs(dna, 47, 64, p,
                PrimerPairOptions.Primer3Defaults with
                {
                    PickInternalOligo = true,
                    InternalOligo = new ProbeDesigner.Primer3ProbeSettings { MaxSelfEnd = -1 },
                }));
        });
    }

    [Test]
    public void DesignProbesPrimer3_AlignmentMode_MatchesPrimer3PickHybProbeOnly()
    {
        // PRIMER_TASK pick_hyb_probe_only, PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_INTERNAL_MAX_SELF_ANY 6.
        const string t = "ATTCTCAGAGGCTCGTACAAACGTATGCCCTAGCTTTTACCACTTAACGCCGTCAAAATGTGCCTATTTTGGAACGAAGGATTCTGTCCTTCGTTCCTTCTTAGTAT";
        var s = new ProbeDesigner.Primer3ProbeSettings { ThermodynamicOligoAlignment = false, MaxSelfAny = 6.0 };
        var probes = ProbeDesigner.DesignProbesPrimer3(t, s);

        var expected = new (int Start, int Length, double Penalty, double SelfAny, double SelfEnd, double Tm)[]
        {
            (8, 23, 4.850867227162098, 6, 2, 58.1491327728379),
            (9, 22, 4.918996372128504, 6, 2, 57.081003627871496),
            (8, 22, 4.918996372128504, 6, 4, 57.081003627871496),
            (7, 23, 5.154764230098806, 6, 4, 57.845235769901194),
            (7, 24, 5.163222934230475, 6, 2, 58.836777065769525),
        };
        Assert.That(probes, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var e = expected[k];
                Assert.That(probes[k].Start, Is.EqualTo(e.Start), $"rank {k}");
                Assert.That(probes[k].Length, Is.EqualTo(e.Length), $"rank {k}");
                Assert.That(probes[k].Penalty, Is.EqualTo(e.Penalty).Within(1e-9), $"rank {k} PRIMER_INTERNAL_PENALTY");
                Assert.That(probes[k].SelfAny, Is.EqualTo(e.SelfAny), $"rank {k} PRIMER_INTERNAL_SELF_ANY");
                Assert.That(probes[k].SelfEnd, Is.EqualTo(e.SelfEnd), $"rank {k} PRIMER_INTERNAL_SELF_END");
                Assert.That(probes[k].Tm, Is.EqualTo(e.Tm).Within(1e-9), $"rank {k} PRIMER_INTERNAL_TM");
                Assert.That(probes[k].SelfAnyTh, Is.NaN);
            }
        });
    }

    [Test]
    public void DesignProbes_LongProbeFallback_UsesPrimer3AlignmentSelfAny()
    {
        // 80-nt stem-loop probe (35-bp inverted repeat around T10): beyond ntthal's 60 nt; dpal.c self_any = 60.00
        // > PRIMER_INTERNAL_MAX_SELF_ANY 12.00. A random 80-mer: self_any 7.00, self_end 1.00 (no self-dimer warning).
        const string hairpin = "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCTTTTTTTTTTGCCGGGGTTGGAGTGAGCAGTGTAGACTGTGATCC";
        const string random = "CCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCA";
        var p = new ProbeDesigner.ProbeParameters(80, 80, -1000, 1000, 0, 1, 100, true, 1.0);
        var flagged = ProbeDesigner.DesignProbes(hairpin, p).Single();
        var clean = ProbeDesigner.DesignProbes(random, p).Single();
        var validation = ProbeDesigner.ValidateProbe(random, Enumerable.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculatePrimerSelfAnyComplementarity(hairpin), Is.EqualTo(60.0));
            Assert.That(flagged.Warnings, Has.Some.EqualTo("Self-complementarity: Primer3 self_any exceeds 12.00"));
            Assert.That(clean.Warnings, Has.None.StartWith("Self-complementarity"));
            Assert.That(validation.SelfAny, Is.EqualTo(7.0));
            Assert.That(validation.SelfEnd, Is.EqualTo(1.0));
            Assert.That(validation.ThermodynamicScreen, Is.False);
        });
    }
}
