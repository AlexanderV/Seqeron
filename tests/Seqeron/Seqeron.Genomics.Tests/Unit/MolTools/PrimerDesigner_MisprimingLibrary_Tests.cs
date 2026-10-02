// PRIMER-DESIGN-001 — Primer3 mispriming library (audit round 3, A3-3 part 1): PRIMER_MISPRIMING_LIBRARY,
// PRIMER_MAX_LIBRARY_MISPRIMING, PRIMER_PAIR_MAX_LIBRARY_MISPRIMING, PRIMER_WT_LIBRARY_MISPRIMING,
// PRIMER_PAIR_WT_LIBRARY_MISPRIMING, PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS.
// Source: primer3 libprimer3.c oligo_repeat_library_mispriming (left: align(s, seqs[i], local_end[_ambig]); right:
//         align(s_r, rev_compl_seqs[i], consensus ? local_end_ambig : local); repeat_sim.max quirk), pair_repeat_sim,
//         p_obj_fn / obj_fn repeat_sim terms, _pr_data_control; p3_seq_lib.c (parse_seq_name, upcase_and_check_char,
//         reverse_complement_seq_lib); dpal.c (set_dpal_args, dpal_set_ambiguity_code_matrix,
//         _dpal_long_nopath_maxgap1_local_end).
// Expected values: primer3-py 2.3.1 design_primers(..., misprime_lib=...) (PRIMER_TASK check_primers for single primers).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_MisprimingLibrary_Tests
{
    private const string T = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";
    private const string Left = "CCAGAAAATAGCGACGGACC";
    private const string Right = "AGCTACTTACAACCTTTCG";

    // lib2 / lib3 of the cross-check script (insertion order = Primer3 scan order).
    private static PrimerMisprimingLibrary Lib2() => new(new[]
    {
        KeyValuePair.Create("frag", "CCAGAAAATAGCGACGGACCGCGG"),
        KeyValuePair.Create("amb*0.5", "CCRGAAAATNGCGACGGWCC"),
        KeyValuePair.Create("rc", "AGGTTGTAAGTAGCTNNNN"),
        KeyValuePair.Create("tiny", "GA"),
    });

    private static PrimerMisprimingLibrary Lib3() => new(new[]
    {
        KeyValuePair.Create("site1", "TGTATAGTCCCACCTGGTGATCCTATGCTTGTGAG"),            // T[60:95]
        KeyValuePair.Create("site2*0.5", "gccagaaggctgcaactcatcgactctatg"),             // T[155:185], weight 0.5
        KeyValuePair.Create("iupac", "NNRYCCAGAAAATAGCGWSKMBDHV"),
    });

    [Test]
    public void CalculateLibraryMispriming_MatchesPrimer3CheckPrimers()
    {
        // check_primers (SEQUENCE_PRIMER / SEQUENCE_PRIMER_REVCOMP): PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS 1 → LEFT (20.0, frag),
        // RIGHT (6.0, reverse rc); 0 (Primer3's default) → LEFT (20.0, frag), RIGHT (15.0, rc) — the right primer is then
        // aligned with dpal LOCAL (not 3′-anchored).
        // lib3: consensus 1 → LEFT (18.0, iupac), RIGHT (5.0, iupac); consensus 0 → LEFT (3.0, site1), RIGHT (5.0, site1).
        var lib2 = Lib2();
        var lib3 = Lib3();
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Left, true, lib2, true), Is.EqualTo(new LibraryMisprimingScore(20.0, "frag")));
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Right, false, lib2, true), Is.EqualTo(new LibraryMisprimingScore(6.0, "reverse rc")));
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Left, true, lib2, false), Is.EqualTo(new LibraryMisprimingScore(20.0, "frag")));
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Right, false, lib2, false), Is.EqualTo(new LibraryMisprimingScore(15.0, "rc")));
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Left, true, lib3, true), Is.EqualTo(new LibraryMisprimingScore(18.0, "iupac")));
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Right, false, lib3, true), Is.EqualTo(new LibraryMisprimingScore(5.0, "iupac")));
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Left, true, lib3, false), Is.EqualTo(new LibraryMisprimingScore(3.0, "site1")));
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Right, false, lib3, false), Is.EqualTo(new LibraryMisprimingScore(5.0, "site1")));
            // Default = Primer3's PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS 0 (pr_set_default_global_args_2).
            Assert.That(PrimerDesigner.CalculateLibraryMispriming(Left, true, lib3), Is.EqualTo(new LibraryMisprimingScore(3.0, "site1")));
        });
    }

    [Test]
    public void CalculateLibraryMispriming_ShortEntryScoresItsLength()
    {
        // libprimer3.c align(): for DPAL_LOCAL / DPAL_LOCAL_END a second sequence shorter than 3 scores strlen(s2);
        // check_primers (template (ACGT)5 T10 (ACGT)5, SEQUENCE_PRIMER T10) → PRIMER_LEFT_0_LIBRARY_MISPRIMING (6.0, x*3).
        var lib = new PrimerMisprimingLibrary(new[] { KeyValuePair.Create("x*3", "GA") });
        Assert.That(PrimerDesigner.CalculateLibraryMispriming("TTTTTTTTTT", true, lib), Is.EqualTo(new LibraryMisprimingScore(6.0, "x*3")));
    }

    [Test]
    public void Library_BuildsPrimer3EntriesWeightsAndReverseComplements()
    {
        // p3_seq_lib.c: weight after '*' (strtod), whitespace removed, upper-cased, IUPAC kept, other characters → N
        // (warning), then "reverse <name>" entries with the IUPAC-aware reverse complement and the same weight.
        var lib = new PrimerMisprimingLibrary(new[]
        {
            KeyValuePair.Create("a*2.5x", "acg tR\n"),
            KeyValuePair.Create("b", "AXZ"),
        });
        Assert.Multiple(() =>
        {
            Assert.That(lib.Count, Is.EqualTo(2));
            Assert.That(lib.Names, Is.EqualTo(new[] { "a*2.5x", "b", "reverse a*2.5x", "reverse b" }));
            Assert.That(lib.SequenceAt(0), Is.EqualTo("ACGTR"));
            Assert.That(lib.SequenceAt(1), Is.EqualTo("ANN"));
            Assert.That(lib.SequenceAt(2), Is.EqualTo("YACGT"));
            Assert.That(lib.SequenceAt(3), Is.EqualTo("NNT"));
            Assert.That(lib.WeightAt(0), Is.EqualTo(2.5));
            Assert.That(lib.WeightAt(1), Is.EqualTo(1.0));
            Assert.That(lib.WeightAt(2), Is.EqualTo(2.5));
            Assert.That(lib.Warnings, Has.Count.EqualTo(1));
        });
    }

    [TestCase("e*")]       // no number after '*'
    [TestCase("e*-1")]     // negative
    [TestCase("e*100.5")]  // > PR_MAX_LIBRARY_WT
    [TestCase("e*inf")]
    public void Library_IllegalWeight_Throws(string name) =>
        Assert.Throws<ArgumentException>(() => _ = new PrimerMisprimingLibrary(new[] { KeyValuePair.Create(name, "ACGT") }));

    [Test]
    public void Library_EmptySequence_Throws() =>
        Assert.Throws<ArgumentException>(() => _ = new PrimerMisprimingLibrary(new[] { KeyValuePair.Create("e", "") }));

    private static readonly (int L, int R, double Pen, double LPen, double RPen, double LLib, string LName,
        double RLib, string RName, double PLib, string PName)[] Lib3Defaults =
    {
        (26, 180, 0.42552744112833807, 0.3171029639751737, 0.10842447715316439, 4.0, "reverse site1", 10.0, "site2*0.5", 11.0, "site2*0.5"),
        (26, 181, 0.42552744112833807, 0.3171029639751737, 0.10842447715316439, 4.0, "reverse site1", 10.0, "site2*0.5", 11.0, "site2*0.5"),
        (25, 180, 0.48597824727789884, 0.37755377012473446, 0.10842447715316439, 3.5, "reverse site2*0.5", 10.0, "site2*0.5", 12.0, "site2*0.5"),
        (25, 181, 0.48597824727789884, 0.37755377012473446, 0.10842447715316439, 3.5, "reverse site2*0.5", 10.0, "site2*0.5", 12.0, "site2*0.5"),
        (26, 201, 0.565601268359103, 0.3171029639751737, 0.24849830438392928, 4.0, "reverse site1", 5.0, "site1", 7.0, "site1"),
    };

    private static readonly (int L, int R, double Pen, double LPen, double RPen, double LLib, string LName,
        double RLib, string RName, double PLib, string PName)[] Lib3AlignmentWeighted =
    {
        (29, 231, 6.061258774732972, 2.409426757397455, 2.2518320173355164, 3.5, "reverse site2*0.5", 4.0, "site1", 7.0, "site1"),
        (26, 231, 6.16893498131069, 2.3171029639751737, 2.2518320173355164, 4.0, "reverse site1", 4.0, "site1", 8.0, "reverse site1"),
        (28, 231, 6.1801617454884425, 2.528329728152926, 2.2518320173355164, 3.0, "site1", 4.0, "site1", 7.0, "site1"),
        (27, 231, 6.4301617454884425, 2.778329728152926, 2.2518320173355164, 3.5, "reverse site2*0.5", 4.0, "site1", 7.0, "site1"),
        (26, 201, 6.465601268359103, 2.3171029639751737, 2.7484983043839293, 4.0, "reverse site1", 5.0, "site1", 7.0, "site1"),
    };

    private static readonly (int L, int R, double Pen, double LPen, double RPen, double LLib, string LName,
        double RLib, string RName, double PLib, string PName)[] Lib3NoConsensusTight =
    {
        (26, 201, 0.565601268359103, 0.3171029639751737, 0.24849830438392928, 4.0, "reverse site1", 5.0, "site1", 7.0, "site1"),
        (26, 231, 0.5689349813106901, 0.3171029639751737, 0.2518320173355164, 4.0, "reverse site1", 4.0, "site1", 8.0, "reverse site1"),
        (26, 269, 0.5701468256208955, 0.3171029639751737, 0.2530438616457218, 4.0, "reverse site1", 8.0, "site1", 10.0, "site1"),
        (25, 201, 0.6260520745086637, 0.37755377012473446, 0.24849830438392928, 3.5, "reverse site2*0.5", 5.0, "site1", 8.0, "site1"),
        (25, 231, 0.6293857874602509, 0.37755377012473446, 0.2518320173355164, 3.5, "reverse site2*0.5", 4.0, "site1", 7.0, "site1"),
    };

    private static readonly (int L, int R, double Pen, double LPen, double RPen, double LLib, string LName,
        double RLib, string RName, double PLib, string PName)[] Lib3Consensus =
    {
        (26, 180, 0.42552744112833807, 0.3171029639751737, 0.10842447715316439, 9.0, "iupac", 10.0, "site2*0.5", 17.0, "iupac"),
        (26, 181, 0.42552744112833807, 0.3171029639751737, 0.10842447715316439, 9.0, "iupac", 10.0, "site2*0.5", 16.0, "iupac"),
        (25, 180, 0.48597824727789884, 0.37755377012473446, 0.10842447715316439, 9.0, "reverse iupac", 10.0, "site2*0.5", 16.0, "iupac"),
        (25, 181, 0.48597824727789884, 0.37755377012473446, 0.10842447715316439, 9.0, "reverse iupac", 10.0, "site2*0.5", 15.0, "iupac"),
        (26, 201, 0.565601268359103, 0.3171029639751737, 0.24849830438392928, 9.0, "iupac", 8.0, "iupac", 17.0, "iupac"),
    };

    private static void AssertPairs(IReadOnlyList<PrimerPairResult> pairs,
        (int L, int R, double Pen, double LPen, double RPen, double LLib, string LName,
            double RLib, string RName, double PLib, string PName)[] expected)
    {
        Assert.That(pairs, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var e = expected[k];
                var p = pairs[k];
                Assert.That(p.Forward!.Position, Is.EqualTo(e.L), $"rank {k} left");
                Assert.That(p.Reverse!.Position + p.Reverse.Length - 1, Is.EqualTo(e.R), $"rank {k} right");
                Assert.That(p.PairPenalty!.Value, Is.EqualTo(e.Pen).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(p.Forward.Penalty, Is.EqualTo(e.LPen).Within(1e-9), $"rank {k} PRIMER_LEFT_PENALTY");
                Assert.That(p.Reverse.Penalty, Is.EqualTo(e.RPen).Within(1e-9), $"rank {k} PRIMER_RIGHT_PENALTY");
                Assert.That(p.Forward.LibraryMispriming, Is.EqualTo(e.LLib), $"rank {k} PRIMER_LEFT_LIBRARY_MISPRIMING");
                Assert.That(p.Forward.LibraryMisprimingName, Is.EqualTo(e.LName));
                Assert.That(p.Reverse.LibraryMispriming, Is.EqualTo(e.RLib), $"rank {k} PRIMER_RIGHT_LIBRARY_MISPRIMING");
                Assert.That(p.Reverse.LibraryMisprimingName, Is.EqualTo(e.RName));
                Assert.That(p.LibraryMispriming, Is.EqualTo(e.PLib), $"rank {k} PRIMER_PAIR_LIBRARY_MISPRIMING");
                Assert.That(p.LibraryMisprimingName, Is.EqualTo(e.PName));
            }
        });
    }

    [Test]
    public void DesignPrimerPairs_Library_RejectsMisprimingPrimerAndMatchesPrimer3()
    {
        // design_primers(SEQUENCE_TARGET [100, 30], PRIMER_PAIR_MAX_DIFF_TM 100, misprime_lib = lib3): the library-free
        // rank 0 left primer [68,20] scores 20 > 12 against site1 and is rejected ("high repeat similarity 1").
        var noLib = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 100, 130, PrimerDesigner.Primer3DefaultParameters,
            PrimerPairOptions.Primer3Defaults);
        var pairs = PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 100, 130,
            PrimerDesigner.Primer3DefaultParameters with { MisprimingLibrary = Lib3() }, PrimerPairOptions.Primer3Defaults);
        Assert.Multiple(() =>
        {
            Assert.That(noLib[0].Forward!.Position, Is.EqualTo(68));
            Assert.That(noLib[0].LibraryMispriming, Is.Null);
            Assert.That(noLib[0].Forward!.LibraryMispriming, Is.Null);
        });
        AssertPairs(pairs, Lib3Defaults);
    }

    [Test]
    public void DesignPrimerPairs_AlignmentModeWithLibraryWeights_MatchesPrimer3()
    {
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_WT_LIBRARY_MISPRIMING 0.5, PRIMER_PAIR_WT_LIBRARY_MISPRIMING 0.2.
        var p = PrimerDesigner.Primer3DefaultParameters with
        {
            StructureScreen = PrimerStructureScreen.Primer3Alignment,
            MisprimingLibrary = Lib3(),
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { LibraryMispriming = 0.5 },
        };
        var o = PrimerPairOptions.Primer3Defaults with { Weights = new Primer3PairWeights { LibraryMispriming = 0.2 } };
        AssertPairs(PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 100, 130, p, o), Lib3AlignmentWeighted);
        Assert.That(PrimerDesigner.DesignPrimers(new DnaSequence(T), 100, 130, p, o).PairPenalty!.Value,
            Is.EqualTo(6.061258774732972).Within(1e-9));
    }

    [Test]
    public void DesignPrimerPairs_NoAmbiguityConsensusAndTightLimits_MatchesPrimer3()
    {
        // PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS 0 (explicit; also the default), PRIMER_MAX_LIBRARY_MISPRIMING 9.9 (compared as (short) 9),
        // PRIMER_PAIR_MAX_LIBRARY_MISPRIMING 10.
        var p = PrimerDesigner.Primer3DefaultParameters with
        {
            MisprimingLibrary = Lib3(),
            LibraryAmbiguityCodesConsensus = false,
            MaxLibraryMispriming = 9.9,
        };
        var o = PrimerPairOptions.Primer3Defaults with { MaxLibraryMispriming = 10 };
        AssertPairs(PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 100, 130, p, o), Lib3NoConsensusTight);
    }

    [Test]
    public void DesignPrimerPairs_AmbiguityCodesConsensus_MatchesPrimer3()
    {
        // PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS 1: the IUPAC entry now aligns (N/R/Y/W/S/K/M/B/D/H/V match their bases).
        var p = PrimerDesigner.Primer3DefaultParameters with { MisprimingLibrary = Lib3(), LibraryAmbiguityCodesConsensus = true };
        AssertPairs(PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 100, 130, p, PrimerPairOptions.Primer3Defaults), Lib3Consensus);
    }

    [Test]
    public void DesignPrimers_EveryPrimerOverLibraryLimit_NoPair()
    {
        // PRIMER_MAX_LIBRARY_MISPRIMING 2: PRIMER_LEFT_EXPLAIN "high repeat similarity 328, ok 0" → no pair.
        var p = PrimerDesigner.Primer3DefaultParameters with { MisprimingLibrary = Lib3(), MaxLibraryMispriming = 2 };
        var r = PrimerDesigner.DesignPrimers(new DnaSequence(T), 100, 130, p, PrimerPairOptions.Primer3Defaults);
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 100, 130, p, PrimerPairOptions.Primer3Defaults), Is.Empty);
            Assert.That(r.IsValid, Is.False);
        });
    }

    [Test]
    public void EvaluatePrimer_Library_ReportsScoreIssueAndWeightedPenalty()
    {
        // Left primer T[68..87] (the library-free rank-0 left primer, PRIMER_LEFT_0_PENALTY 0.10934227176778677):
        // check_primers with lib3 → "high repeat similarity 1" (score 20.0 against site1 > 12); with
        // PRIMER_MAX_LIBRARY_MISPRIMING 20 and PRIMER_WT_LIBRARY_MISPRIMING 0.5 → LIBRARY_MISPRIMING (20.0, site1),
        // PRIMER_LEFT_0_PENALTY 10.109342271767787.
        string primer = T.Substring(68, 20);
        var p = PrimerDesigner.Primer3DefaultParameters with { MisprimingLibrary = Lib3() };
        var plain = PrimerDesigner.EvaluatePrimer(primer, 68, true, PrimerDesigner.Primer3DefaultParameters);
        var c = PrimerDesigner.EvaluatePrimer(primer, 68, true, p);
        var weighted = PrimerDesigner.EvaluatePrimer(primer, 68, true,
            p with { MaxLibraryMispriming = 20, PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { LibraryMispriming = 0.5 } });
        Assert.Multiple(() =>
        {
            Assert.That(plain.IsValid, Is.True);
            Assert.That(plain.Penalty, Is.EqualTo(0.10934227176778677).Within(1e-9));
            Assert.That(plain.LibraryMispriming, Is.Null);
            Assert.That(c.LibraryMispriming, Is.EqualTo(20.0));
            Assert.That(c.LibraryMisprimingName, Is.EqualTo("site1"));
            Assert.That(c.IsValid, Is.False);
            Assert.That(c.Issues, Has.Some.Contains("PRIMER_MAX_LIBRARY_MISPRIMING"));
            Assert.That(weighted.IsValid, Is.True);
            Assert.That(weighted.LibraryMispriming, Is.EqualTo(20.0));
            Assert.That(weighted.Penalty, Is.EqualTo(10.109342271767787).Within(1e-9));
        });
    }

    [Test]
    public void LibraryWeightsWithoutLibrary_AndOversizedLimit_Throw()
    {
        // _pr_data_control: "Mispriming score is part of objective function, but mispriming library is not defined";
        // PRIMER_MAX_LIBRARY_MISPRIMING / PRIMER_PAIR_MAX_LIBRARY_MISPRIMING > SHRT_MAX in alignment mode.
        var dna = new DnaSequence(T);
        var weighted = PrimerDesigner.Primer3DefaultParameters with
        {
            PenaltyWeights = PrimerDesigner.DefaultPrimer3Weights with { LibraryMispriming = 1 },
        };
        var pairWeighted = PrimerPairOptions.Primer3Defaults with { Weights = new Primer3PairWeights { LibraryMispriming = 1 } };
        var aln = PrimerDesigner.Primer3DefaultParameters with { StructureScreen = PrimerStructureScreen.Primer3Alignment };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => PrimerDesigner.EvaluatePrimer(Left, 68, true, weighted));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, 100, 130, weighted, PrimerPairOptions.Primer3Defaults));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, 100, 130,
                PrimerDesigner.Primer3DefaultParameters, pairWeighted));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.EvaluatePrimer(Left, 68, true,
                aln with { MaxLibraryMispriming = 40000 }));
            Assert.Throws<ArgumentException>(() => PrimerDesigner.DesignPrimerPairs(dna, 100, 130, aln,
                PrimerPairOptions.Primer3Defaults with { MaxLibraryMispriming = 40000 }));
        });
    }
}
