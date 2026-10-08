// PRIMER-DESIGN-001 / PROBE-DESIGN-001 — audit round 3, A3-25 + A3-26:
//  * PRIMER_[INTERNAL_]OPT_GC_PERCENT undefined by default; a GC weight without it is Primer3's _pr_data_control error.
//  * PRIMER_LOWERCASE_MASKING on a case-preserving template (DesignPrimers / DesignPrimerPairs string overloads,
//    Primer3ProbeSettings.LowercaseMasking).
// Source: primer3-py 2.3.1 libprimer3.c pr_set_default_global_args_1 (p_args / o_args.opt_gc_content =
//         DEFAULT_OPT_GC_PERCENT = PR_UNDEFINED_INT_OPT), _pr_data_control ("Primer GC content is part of objective
//         function while optimum gc_content is not defined", "Hyb probe GC content …" — o_args checked whether or not an
//         internal oligo is picked), calc_and_check_oligo_features → is_lowercase_masked (3'-terminal base of a left primer /
//         internal oligo = its rightmost template base, of a right primer its leftmost; on trimmed_orig_seq, the mixed-case
//         template; a/c/g/t only), thermoanalysis.pyx (lowercase_masking = mask_template).
// Expected values: primer3-py 2.3.1 design_primers on the templates below (scratch h10/ref.py).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_LowercaseMasking_Tests
{
    // random.Random(2026), 240 nt; lower case at [10,30), [70,74), [150,175), 205 and 40, 47, 52, 120, 121, 128, 131, 190, 196, 199.
    private const string T = "AGACTTTCAAagatatgctgggtagaggtcGAGGTTATTAtTTGTTAcCAATtCTCATTGTGTTTCGGAActtgCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGatCCGTTCcTAaTAAGGAATGGTGATTCCCtgtcataccaatctaccccctgttaTGCGCGTTTGTCGTTaGACCAaTGtCAGCGcAGCGGCAGATCAAGCAGGAGGCGGAATGTAAACA";
    // T with 210 and 222 also lower case.
    private const string T2 = "AGACTTTCAAagatatgctgggtagaggtcGAGGTTATTAtTTGTTAcCAATtCTCATTGTGTTTCGGAActtgCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGatCCGTTCcTAaTAAGGAATGGTGATTCCCtgtcataccaatctaccccctgttaTGCGCGTTTGTCGTTaGACCAaTGtCAGCGcAGCGgCAGATCAAGCAgGAGGCGGAATGTAAACA";
    private const double Tol = 1e-9;

    private static readonly PrimerPairOptions Pick = PrimerPairOptions.Primer3Defaults with { PickInternalOligo = true };

    // (left start, len, right start (3' end, Primer3 coordinates), len, PRIMER_PAIR_k_PENALTY, internal start, len, penalty)
    private static void AssertPairs(IReadOnlyList<PrimerPairResult> got,
        (int L, int Ll, int R, int Rl, double Pen, int I, int Il, double IPen)[] exp)
    {
        Assert.That(got, Has.Count.EqualTo(exp.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < exp.Length; k++)
            {
                var g = got[k];
                Assert.That((g.Forward!.Position, g.Forward.Length, g.Reverse!.Position + g.Reverse.Length - 1, g.Reverse.Length),
                    Is.EqualTo((exp[k].L, exp[k].Ll, exp[k].R, exp[k].Rl)), $"rank {k}");
                Assert.That(g.PairPenalty!.Value, Is.EqualTo(exp[k].Pen).Within(Tol), $"rank {k} PRIMER_PAIR_{k}_PENALTY");
                Assert.That((g.InternalOligo!.Value.Start, g.InternalOligo.Value.Length), Is.EqualTo((exp[k].I, exp[k].Il)), $"rank {k} internal");
                Assert.That(g.InternalOligo.Value.Penalty, Is.EqualTo(exp[k].IPen).Within(Tol), $"rank {k} PRIMER_INTERNAL_{k}_PENALTY");
            }
        });
    }

    #region A3-26 PRIMER_LOWERCASE_MASKING

    [Test]
    public void DesignPrimerPairs_LowercaseMasking_MatchesPrimer3()
    {
        // primer3-py: SEQUENCE_TARGET [100,20], PRIMER_LOWERCASE_MASKING 1, PRIMER_PICK_INTERNAL_OLIGO 1 (explain: left / right /
        // internal "lowercase masking of 3' end" 155 / 330 / 485).
        var got = PrimerDesigner.DesignPrimerPairs(T, 100, 120, PrimerDesigner.Primer3DefaultParameters, Pick with { LowercaseMasking = true });
        AssertPairs(got, new[]
        {
            (14, 20, 235, 20, 0.8648905951553161, 190, 21, 1.04795185124226),
            (14, 20, 232, 20, 0.936221987820943, 190, 21, 1.04795185124226),
            (16, 20, 235, 20, 1.0791502276797473, 190, 21, 1.04795185124226),
            (14, 20, 205, 20, 1.1274602336220596, 104, 22, 2.339100021051138),
            (16, 20, 232, 20, 1.1504816203453743, 190, 21, 1.04795185124226),
        });
        Assert.Multiple(() =>
        {
            foreach (var p in got)
            {
                Assert.That(T[p.Forward!.Position + p.Forward.Length - 1], Is.Not.AnyOf('a', 'c', 'g', 't'), "left 3' base");
                Assert.That(T[p.Reverse!.Position], Is.Not.AnyOf('a', 'c', 'g', 't'), "right 3' base");
                Assert.That(T[p.InternalOligo!.Value.Start + p.InternalOligo.Value.Length - 1], Is.Not.AnyOf('a', 'c', 'g', 't'), "internal 3' base");
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_LowercaseMaskingOff_IgnoresCase()
    {
        // primer3-py PRIMER_LOWERCASE_MASKING 0 (default): case is ignored.
        var expected = new[]
        {
            (14, 20, 235, 20, 0.8648905951553161, 190, 21, 1.04795185124226),
            (14, 20, 188, 20, 0.9201203739298762, 104, 22, 2.339100021051138),
            (14, 20, 232, 20, 0.936221987820943, 190, 21, 1.04795185124226),
            (14, 20, 186, 20, 0.9363871256945799, 104, 22, 2.339100021051138),
            (16, 20, 235, 20, 1.0791502276797473, 190, 21, 1.04795185124226),
        };
        AssertPairs(PrimerDesigner.DesignPrimerPairs(T, 100, 120, PrimerDesigner.Primer3DefaultParameters, Pick), expected);
        // A DnaSequence template is upper case: LowercaseMasking has nothing to mask.
        AssertPairs(PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 100, 120, PrimerDesigner.Primer3DefaultParameters,
            Pick with { LowercaseMasking = true }), expected);
    }

    [Test]
    public void DesignPrimers_StringOverload_ReturnsRankZero()
    {
        var r = PrimerDesigner.DesignPrimers(T, 100, 120, PrimerDesigner.Primer3DefaultParameters,
            Pick with { LowercaseMasking = true, NumReturn = 1 });
        Assert.Multiple(() =>
        {
            Assert.That(r.IsValid, Is.True);
            Assert.That((r.Forward!.Position, r.Reverse!.Position + r.Reverse.Length - 1), Is.EqualTo((14, 235)));
            Assert.That(r.PairPenalty!.Value, Is.EqualTo(0.8648905951553161).Within(Tol));
        });
    }

    [Test]
    public void DesignProbesPrimer3_LowercaseMasking_MatchesPrimer3()
    {
        // primer3-py pick_hyb_probe_only (PRIMER_PICK_LEFT/RIGHT_PRIMER 0, PRIMER_PICK_INTERNAL_OLIGO 1) on T2.
        var on = ProbeDesigner.DesignProbesPrimer3(T2, new ProbeDesigner.Primer3ProbeSettings { LowercaseMasking = true });
        var off = ProbeDesigner.DesignProbesPrimer3(T2);
        Assert.Multiple(() =>
        {
            Assert.That(on.Select(p => (p.Start, p.Length)),
                Is.EqualTo(new[] { (204, 21), (204, 20), (207, 21), (190, 20), (203, 21) }));
            Assert.That(on.Select(p => p.Penalty), Is.EqualTo(new[]
                { 1.216988073303753, 1.319185186492689, 1.4403313670503053, 2.001596012010225, 2.0162713139086463 }).Within(Tol));
            Assert.That(off.Select(p => (p.Start, p.Length)),
                Is.EqualTo(new[] { (203, 20), (190, 21), (204, 21), (204, 20), (207, 21) }));
            Assert.That(off.Select(p => p.Penalty), Is.EqualTo(new[]
                { 0.9613419511969141, 1.04795185124226, 1.216988073303753, 1.319185186492689, 1.4403313670503053 }).Within(Tol));
        });
    }

    [Test]
    public void DesignPrimers_StringOverload_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => PrimerDesigner.DesignPrimers((string)null!, 0, 10), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => PrimerDesigner.DesignPrimerPairs((string)null!, 0, 10), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(T[..100] + "nnnn" + T[104..], 100, 120), NUnit.Framework.Throws.ArgumentException);
        });
    }

    #endregion

    #region A3-25 undefined GC optimum

    private const string U = "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAATCTACCCCCTGTTATGCGCGTTTGTCGTTAGACCAATGTCAGCGCAGCGGCAGATCAAGCAGGAGGCGGAATGTAAACA";
    private const string PrimerMsg = "Primer GC content is part of objective function while optimum gc_content is not defined";
    private const string ProbeMsg = "Hyb probe GC content is part of objective function while optimum gc_content is not defined";

    [Test]
    public void GcWeightWithoutOptimum_ThrowsPrimer3DataControlErrors()
    {
        var dna = new DnaSequence(U);
        var p = PrimerDesigner.Primer3DefaultParameters;
        Assert.Multiple(() =>
        {
            // primer3-py: PRIMER_WT_GC_PERCENT_GT 0.5 / _LT 1 without PRIMER_OPT_GC_PERCENT.
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, 100, 120,
                    p with { PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 0.5 } }, PrimerPairOptions.Primer3Defaults),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith(PrimerMsg));
            Assert.That(() => PrimerDesigner.DesignPrimers(dna, 100, 120,
                    p with { PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { GcLt = 1.0 } }, PrimerPairOptions.Primer3Defaults),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith(PrimerMsg));
            Assert.That(() => PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 0, true,
                    p with { PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 0.5 } }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith(PrimerMsg));
            // PRIMER_INTERNAL_WT_GC_PERCENT_GT 0.5 without the optimum: rejected even without PRIMER_PICK_INTERNAL_OLIGO.
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, 100, 120, p, PrimerPairOptions.Primer3Defaults with
                    { InternalOligo = new ProbeDesigner.Primer3ProbeSettings { WeightGcPercentGt = 0.5 } }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith(ProbeMsg));
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, 100, 120, p, Pick with
                    { InternalOligo = new ProbeDesigner.Primer3ProbeSettings { WeightGcPercentLt = 0.5 } }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith(ProbeMsg));
            Assert.That(() => ProbeDesigner.DesignProbesPrimer3(U, new ProbeDesigner.Primer3ProbeSettings { WeightGcPercentGt = 0.5 }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith(ProbeMsg));
            // The primer check precedes the internal-oligo one (_pr_data_control order).
            Assert.That(() => PrimerDesigner.DesignPrimerPairs(dna, 100, 120,
                    p with { PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 0.5 } }, PrimerPairOptions.Primer3Defaults with
                    { InternalOligo = new ProbeDesigner.Primer3ProbeSettings { WeightGcPercentGt = 0.5 } }),
                NUnit.Framework.Throws.ArgumentException.With.Message.StartsWith(PrimerMsg));
        });
    }

    [Test]
    public void GcWeightWithOptimum_MatchesPrimer3()
    {
        // primer3-py: PRIMER_WT_GC_PERCENT_GT 0.5, PRIMER_OPT_GC_PERCENT 45, PRIMER_INTERNAL_WT_GC_PERCENT_LT 1,
        // PRIMER_INTERNAL_OPT_GC_PERCENT 55, PRIMER_PICK_INTERNAL_OLIGO 1.
        var got = PrimerDesigner.DesignPrimerPairs(new DnaSequence(U), 100, 120,
            PrimerDesigner.Primer3DefaultParameters with
            {
                OptimalGcPercent = 45, PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 0.5 },
            },
            Pick with { InternalOligo = new ProbeDesigner.Primer3ProbeSettings { OptGcPercent = 55, WeightGcPercentLt = 1 } });
        AssertPairs(got.Take(3).ToList(), new[]
        {
            (58, 20, 189, 20, 1.670788242645017, 104, 22, 2.339100021051138),
            (57, 20, 189, 20, 1.670788242645017, 104, 22, 2.339100021051138),
            (58, 21, 189, 20, 2.05179072297193, 104, 22, 2.339100021051138),
        });
    }

    #endregion
}
