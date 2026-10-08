// PROBE-DESIGN-001 / PRIMER-DESIGN-001 — Primer3 internal-oligo mishybridization library (audit round 3, A3-3 part 2):
// PRIMER_INTERNAL_MISHYB_LIBRARY, PRIMER_INTERNAL_MAX_LIBRARY_MISHYB, PRIMER_INTERNAL_WT_LIBRARY_MISHYB,
// PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS for the hybridization-probe picker (PRIMER_TASK=pick_hyb_probe_only) and the
// internal oligo of a primer pair (PRIMER_PICK_INTERNAL_OLIGO); PRIMER_NUM_RETURN < 1 (A3-15).
// Source: primer3 libprimer3.c oligo_repeat_library_mispriming (OT_INTL: align(s, seqs[i], consensus ? local_ambig :
//         local)), calc_and_check_oligo_features (scored while picking for list output or when weighted),
//         choose_internal_oligo (otherwise after the postponed structure checks), p_obj_fn OT_INTL repeat_sim term,
//         five_prime_problem (OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ), _pr_data_control.
// Expected values: primer3-py 2.3.1 design_primers(..., mishyb_lib=...) (PRIMER_TASK check_primers for a single oligo;
//         primer3-py reports the score as PRIMER_INTERNAL_n_LIBRARY_MISPRIMING).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class ProbeDesigner_MishybLibrary_Tests
{
    // Probe template (ProbeDesigner_Primer3Probe_Tests) and pair template (PrimerDesigner_MisprimingLibrary_Tests).
    private const string P =
        "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACC";
    private const string T = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";
    private const string Oligo = "GTCCCACCTGGTGATCCTATGC"; // P[68:90]

    private static PrimerMisprimingLibrary LibA() => new(new[]
    {
        KeyValuePair.Create("site*2", "CCCACCTGGTGATCCTATGCTTGTG"), // P[70:95], weight 2
        KeyValuePair.Create("mut", "CCCACCTGGTGTTCCTATGC"),
        KeyValuePair.Create("amb", "GTCCCRCCTGGNGATCCTAYG"),
        KeyValuePair.Create("tiny", "CA"),
    });

    private static PrimerMisprimingLibrary LibC() => new(new[]
    {
        KeyValuePair.Create("five", "GTCCCACCTGGTGAAAAAAAA"),
        KeyValuePair.Create("amb", "GTCCCRCCTGGNGATCCTAYG"),
        KeyValuePair.Create("mut", "CCCACCTGGTGTTCCTATGC"),
    });

    private static PrimerMisprimingLibrary LibB() => new(new[]
    {
        KeyValuePair.Create("s1", "GCGGTGTTAAGTGTCGAGCTACATCACTTCTC"),                     // T[118:150]
        KeyValuePair.Create("s2*0.5", "atgtagccagaaggctgcaactcatcgactctatg"),             // T[150:185], weight 0.5
        KeyValuePair.Create("iu", "GGTGTTAAGTGTCRAGCTACAYC"),
    });

    [Test]
    public void CalculateLibraryMishyb_MatchesPrimer3CheckPrimers()
    {
        // check_primers (SEQUENCE_INTERNAL_OLIGO, PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 100): libA → (40.0, site*2) for both
        // consensus settings; libC consensus 0 → (18.0, mut) (IUPAC codes never align), consensus 1 → (21.0, amb).
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateLibraryMishyb(Oligo, LibA()), Is.EqualTo(new LibraryMisprimingScore(40.0, "site*2")));
            Assert.That(PrimerDesigner.CalculateLibraryMishyb(Oligo, LibA(), true), Is.EqualTo(new LibraryMisprimingScore(40.0, "site*2")));
            Assert.That(PrimerDesigner.CalculateLibraryMishyb(Oligo.ToLowerInvariant(), LibC()), Is.EqualTo(new LibraryMisprimingScore(18.0, "mut")));
            Assert.That(PrimerDesigner.CalculateLibraryMishyb(Oligo, LibC(), true), Is.EqualTo(new LibraryMisprimingScore(21.0, "amb")));
        });
    }

    [Test]
    public void CalculateLibraryMishyb_ShortEntryScoresItsLength_AndNullsThrow()
    {
        // align(): an entry shorter than 3 scores strlen; check_primers (SEQUENCE_INTERNAL_OLIGO P[0:20], x*3 = "GA") → (6.0, x*3).
        var lib = new PrimerMisprimingLibrary(new[] { KeyValuePair.Create("x*3", "GA") });
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateLibraryMishyb(P[..20], lib), Is.EqualTo(new LibraryMisprimingScore(6.0, "x*3")));
            Assert.Throws<ArgumentNullException>(() => PrimerDesigner.CalculateLibraryMishyb(null!, lib));
            Assert.Throws<ArgumentNullException>(() => PrimerDesigner.CalculateLibraryMishyb(Oligo, null!));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.CalculateLibraryMishyb("", lib));
        });
    }

    private static void AssertProbes(IReadOnlyList<ProbeDesigner.Primer3Probe> probes,
        (int Start, int Len, double Penalty, double Score, string Name)[] expected)
    {
        Assert.That(probes, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int i = 0; i < expected.Length; i++)
            {
                var (p, e) = (probes[i], expected[i]);
                Assert.That((p.Start, p.Length), Is.EqualTo((e.Start, e.Len)), $"rank {i}");
                Assert.That(p.Penalty, Is.EqualTo(e.Penalty).Within(1e-9), $"rank {i} penalty");
                Assert.That(p.LibraryMishyb, Is.EqualTo(e.Score), $"rank {i} score");
                Assert.That(p.LibraryMishybName, Is.EqualTo(e.Name), $"rank {i} name");
            }
        });
    }

    [Test]
    public void DesignProbesPrimer3_MishybLibrary_MatchesPrimer3PickHybProbeOnly()
    {
        // pick_hyb_probe_only, defaults + mishyb_lib libA: "high repeat similarity 27, ok 62"; ranks 2–4 change
        // (without the library: [67,24] [68,23] [67,23], which align with site*2 above 12).
        AssertProbes(ProbeDesigner.DesignProbesPrimer3(P, new ProbeDesigner.Primer3ProbeSettings { MishybLibrary = LibA() }), new[]
        {
            (27, 22, 4.426338881812228, 8.0, "site*2"),
            (27, 23, 5.073673806189504, 8.0, "site*2"),
            (26, 23, 5.808229853859814, 8.0, "site*2"),
            (25, 24, 6.080699673260824, 8.0, "site*2"),
            (24, 25, 6.119601327140799, 8.0, "site*2"),
        });
    }

    [Test]
    public void DesignProbesPrimer3_WeightedMishyb_AddsPrimer3PenaltyTerm()
    {
        // PRIMER_INTERNAL_WT_LIBRARY_MISHYB 0.5, PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 30: penalty + 0.5 × score.
        var s = new ProbeDesigner.Primer3ProbeSettings { MishybLibrary = LibA(), WeightLibraryMishyb = 0.5, MaxLibraryMishyb = 30 };
        AssertProbes(ProbeDesigner.DesignProbesPrimer3(P, s), new[]
        {
            (27, 22, 8.426338881812228, 8.0, "site*2"),
            (27, 23, 9.073673806189504, 8.0, "site*2"),
            (26, 23, 9.808229853859814, 8.0, "site*2"),
            (25, 24, 10.080699673260824, 8.0, "site*2"),
            (24, 25, 10.1196013271408, 8.0, "site*2"),
        });
    }

    [Test]
    public void DesignProbesPrimer3_AlignmentModeConsensusMishyb_MatchesPrimer3()
    {
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS 1, MAX_LIBRARY_MISHYB 8:
        // "high repeat similarity 38, ok 37".
        var s = new ProbeDesigner.Primer3ProbeSettings
        {
            ThermodynamicOligoAlignment = false, MishybLibrary = LibA(), LibraryAmbiguityCodesConsensus = true, MaxLibraryMishyb = 8,
        };
        AssertProbes(ProbeDesigner.DesignProbesPrimer3(P, s), new[]
        {
            (27, 22, 4.426338881812228, 8.0, "site*2"),
            (27, 23, 5.073673806189504, 8.0, "site*2"),
            (26, 23, 5.808229853859814, 8.0, "site*2"),
            (25, 24, 6.080699673260824, 8.0, "site*2"),
            (24, 25, 6.119601327140799, 8.0, "site*2"),
        });
    }

    [Test]
    public void DesignProbesPrimer3_WithoutLibrary_ReportsNoMishybScore()
    {
        var probes = ProbeDesigner.DesignProbesPrimer3(P);
        Assert.That(probes.Select(p => p.LibraryMishyb), Is.All.Null);
        Assert.That(probes.Select(p => (p.Start, p.Length)).Skip(2), Is.EqualTo(new[] { (67, 24), (68, 23), (67, 23) }));
    }

    [Test]
    public void DesignProbesPrimer3_Primer3DataControl()
    {
        // primer3-py: PRIMER_NUM_RETURN 0 → "PRIMER_NUM_RETURN < 1"; PRIMER_INTERNAL_WT_LIBRARY_MISHYB 1 without a library →
        // "Internal oligo mispriming score is part of objective function while mishyb library is not defined";
        // alignment mode with PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 40000 → "Value too large at tag PRIMER_INTERNAL_MAX_LIBRARY_MISHYB".
        Assert.Multiple(() =>
        {
            var e0 = Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbesPrimer3(P, numReturn: 0));
            Assert.That(e0!.Message, Does.Contain("PRIMER_NUM_RETURN < 1"));
            Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbesPrimer3(P, numReturn: -1));
            var e1 = Assert.Throws<ArgumentException>(() => ProbeDesigner.DesignProbesPrimer3(P,
                new ProbeDesigner.Primer3ProbeSettings { WeightLibraryMishyb = 1.0 }));
            Assert.That(e1!.Message, Does.Contain("mishyb library is not defined"));
            var e2 = Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbesPrimer3(P,
                new ProbeDesigner.Primer3ProbeSettings { ThermodynamicOligoAlignment = false, MaxLibraryMishyb = 40000 }));
            Assert.That(e2!.Message, Does.Contain("PRIMER_INTERNAL_MAX_LIBRARY_MISHYB"));
            // Thermodynamic mode: no SHRT_MAX check, and Primer3 compares with (short) 40000 = −25536, so every oligo is
            // rejected (primer3-py: "high repeat similarity 49, ok 0"); without a library the limit is unused.
            Assert.That(ProbeDesigner.DesignProbesPrimer3(P,
                new ProbeDesigner.Primer3ProbeSettings { MaxLibraryMishyb = 40000, MishybLibrary = LibA() }), Is.Empty);
            Assert.That(ProbeDesigner.DesignProbesPrimer3(P,
                new ProbeDesigner.Primer3ProbeSettings { MaxLibraryMishyb = 40000 }), Has.Count.EqualTo(5));
        });
    }

    private static IReadOnlyList<PrimerPairResult> Pairs(PrimerParameters p, ProbeDesigner.Primer3ProbeSettings s, double ioWeight = 0) =>
        PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 100, 130, p, PrimerPairOptions.Primer3Defaults with
        {
            PickInternalOligo = true,
            NumReturn = 5,
            InternalOligo = s,
            Weights = new Primer3PairWeights { InternalOligoPenalty = ioWeight },
        });

    private static void AssertPairs(IReadOnlyList<PrimerPairResult> pairs,
        (int L, int R, double PairPenalty, int IStart, int ILen, double IPenalty, double Score)[] expected, string name)
    {
        Assert.That(pairs, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int i = 0; i < expected.Length; i++)
            {
                var (r, e) = (pairs[i], expected[i]);
                var io = r.InternalOligo!.Value;
                Assert.That((r.Forward!.Position, r.Reverse!.Position + r.Reverse.Length - 1), Is.EqualTo((e.L, e.R)), $"rank {i}");
                Assert.That(r.PairPenalty, Is.EqualTo(e.PairPenalty).Within(1e-9), $"rank {i} pair penalty");
                Assert.That((io.Start, io.Length), Is.EqualTo((e.IStart, e.ILen)), $"rank {i} internal oligo");
                Assert.That(io.Penalty, Is.EqualTo(e.IPenalty).Within(1e-9), $"rank {i} internal penalty");
                Assert.That(io.LibraryMishyb, Is.EqualTo(e.Score), $"rank {i} score");
                Assert.That(io.LibraryMishybName, Is.EqualTo(name), $"rank {i} name");
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_InternalOligoMishyb_PostponedCheck_MatchesPrimer3()
    {
        // PRIMER_PICK_INTERNAL_OLIGO 1, PRIMER_PAIR_MAX_DIFF_TM 100, target [100,30], mishyb_lib libB.
        // Default limit 12: the best internal oligo [108,20] scores 10 (s1) and is kept.
        AssertPairs(Pairs(PrimerDesigner.Primer3DefaultParameters, new ProbeDesigner.Primer3ProbeSettings { MishybLibrary = LibB() }), new[]
        {
            (68, 180, 0.21776674892095116, 108, 20, 0.04565117322755441, 10.0),
            (68, 181, 0.21776674892095116, 108, 20, 0.04565117322755441, 10.0),
            (68, 201, 0.35784057615171605, 108, 20, 0.04565117322755441, 10.0),
            (68, 231, 0.3611742891033032, 108, 20, 0.04565117322755441, 10.0),
            (68, 269, 0.36238613341350856, 108, 20, 0.04565117322755441, 10.0),
        }, "s1");

        // PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 9: [108,20] is rejected in choose_internal_oligo ("high repeat similarity 80").
        AssertPairs(Pairs(PrimerDesigner.Primer3DefaultParameters,
            new ProbeDesigner.Primer3ProbeSettings { MishybLibrary = LibB(), MaxLibraryMishyb = 9 }), new[]
        {
            (68, 180, 0.21776674892095116, 108, 19, 1.2417931164713423, 9.0),
            (68, 181, 0.21776674892095116, 108, 19, 1.2417931164713423, 9.0),
            (68, 201, 0.35784057615171605, 108, 19, 1.2417931164713423, 9.0),
            (68, 231, 0.3611742891033032, 108, 19, 1.2417931164713423, 9.0),
            (68, 269, 0.36238613341350856, 194, 21, 1.2143916852965049, 5.0),
        }, "s1");
    }

    [Test]
    public void DesignPrimerPairs_WeightedInternalOligoMishyb_MatchesPrimer3()
    {
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_INTERNAL_WT_LIBRARY_MISHYB 1, PRIMER_PAIR_WT_IO_PENALTY 1,
        // PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 25: scored while picking ("high repeat similarity 4"); the internal-oligo
        // penalty (+ score) enters the pair objective. (Alignment mode keeps the ~10⁴ characterized pairs fast; the
        // thermodynamic-mode equivalent is part of the B07 F42 primer3-py cross-check.)
        var p = PrimerDesigner.Primer3DefaultParameters with { StructureScreen = PrimerStructureScreen.Primer3Alignment };
        var s = new ProbeDesigner.Primer3ProbeSettings { MishybLibrary = LibB(), WeightLibraryMishyb = 1.0, MaxLibraryMishyb = 25 };
        AssertPairs(Pairs(p, s, ioWeight: 1.0), new[]
        {
            (68, 172, 6.27132541268503, 103, 20, 6.125973141533393, 5.0),
            (68, 167, 6.342823635136256, 103, 20, 6.125973141533393, 5.0),
            (68, 180, 6.343739890454344, 103, 20, 6.125973141533393, 5.0),
            (68, 181, 6.343739890454344, 103, 20, 6.125973141533393, 5.0),
            (26, 172, 6.479086104892417, 103, 20, 6.125973141533393, 5.0),
        }, "s1");
    }

    [Test]
    public void DesignPrimerPairs_AlignmentModeConsensusMishyb_UsesGlobalConsensusSetting()
    {
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS 1 (global: PrimerParameters),
        // PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 10.
        var p = PrimerDesigner.Primer3DefaultParameters with
        {
            StructureScreen = PrimerStructureScreen.Primer3Alignment,
            LibraryAmbiguityCodesConsensus = true,
        };
        AssertPairs(Pairs(p, new ProbeDesigner.Primer3ProbeSettings { MishybLibrary = LibB(), MaxLibraryMishyb = 10 }), new[]
        {
            (68, 172, 0.14535227115163707, 108, 20, 0.04565117322755441, 10.0),
            (68, 167, 0.21685049360286257, 108, 20, 0.04565117322755441, 10.0),
            (68, 180, 0.21776674892095116, 108, 20, 0.04565117322755441, 10.0),
            (68, 181, 0.21776674892095116, 108, 20, 0.04565117322755441, 10.0),
            (26, 172, 0.353112963359024, 108, 20, 0.04565117322755441, 10.0),
        }, "s1");
    }

    [Test]
    public void DesignPrimerPairs_InternalOligoMishybDataControl()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => Pairs(PrimerDesigner.Primer3DefaultParameters,
                new ProbeDesigner.Primer3ProbeSettings { WeightLibraryMishyb = 1.0 }));
            // The SHRT_MAX limit follows the (global) alignment mode of the primer screen.
            Assert.Throws<ArgumentOutOfRangeException>(() => Pairs(
                PrimerDesigner.Primer3DefaultParameters with { StructureScreen = PrimerStructureScreen.Primer3Alignment },
                new ProbeDesigner.Primer3ProbeSettings { MaxLibraryMishyb = 40000 }));
        });
    }
}
