// PRIMER-DESIGN-001 — Primer3 3′-end checks and 3′-distance pair constraint (audit round 3, A3-6 + A3-7).
// Source: primer3 libprimer3.c calc_and_check_oligo_features (PRIMER_GC_CLAMP over the gc_clamp 3′-most bases,
//         PRIMER_MAX_END_GC over the five 3′-most bases when < 5, PRIMER_MAX_END_STABILITY: end_oligodg(seq, 5) > max fails;
//         left/right primers only), _pr_data_control (GC_CLAMP > MIN_SIZE, MAX_END_GC ∉ [0, 5], MAX_END_STABILITY < 0,
//         min_*_three_prime_distance < −1), choose_pair_or_triple + left/right_oligo_in_pair_overlaps_used_oligo;
//         read_boulder.c PRIMER_MIN_THREE_PRIME_DISTANCE / _LEFT_ / _RIGHT_ tags.
// Expected values: primer3-py 2.3.1 design_primers (design and check_primers tasks).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class PrimerDesigner_ThreePrimeEnd_Tests
{
    private const string T = "GGTCCTAATTGGAGCGCCCAGTTACCGGCCGAGTGCTACGGGCACTCGTTGGTAGTGGGCTCCCTAAGTCGGCGCATCCGTTCCTAGCTTTAAAATATCCGTTGAAAGAATGTTCTGAGTCTCGCCTAGTGAAAGCCAACTCCTTTGGATTTTGTCATA";

    private static readonly PrimerParameters P3 = PrimerDesigner.Primer3DefaultParameters;

    private static bool HasIssue(PrimerCandidate c, string tag) => c.Issues.Any(i => i.Contains(tag, StringComparison.Ordinal));

    #region EvaluatePrimer — per-primer 3′-end checks (primer3-py check_primers)

    [TestCase("GTTACCGGCCGAGTGCTAC", 1, false)]
    [TestCase("GTTACCGGCCGAGTGCTAC", 2, true)]   // "GC clamp failed 1"
    [TestCase("CAGTTACCGGCCGAGTGCGC", 3, false)]
    [TestCase("AGTTACCGGCCGAGTGCTA", 1, true)]
    [TestCase("TCCTAATTGGAGCGCCCA", 1, true)]
    public void EvaluatePrimer_GcClamp_MatchesPrimer3(string primer, int clamp, bool fails)
    {
        var c = PrimerDesigner.EvaluatePrimer(primer, 0, true, P3 with { GcClamp = clamp });
        Assert.That(HasIssue(c, "PRIMER_GC_CLAMP"), Is.EqualTo(fails));
    }

    [TestCase("GTTACCGGCCGAGTGCTAC", 2, true)]   // TGCTAC: 3 G/C
    [TestCase("GTTACCGGCCGAGTGCTAC", 3, false)]
    [TestCase("CAGTTACCGGCCGAGTGCGC", 3, true)]  // TGCGC: 4 G/C
    [TestCase("CAGTTACCGGCCGAGTGCGC", 4, false)]
    [TestCase("AGTTACCGGCCGAGTGCTA", 2, false)]
    [TestCase("TCCTAATTGGAGCGCCCA", 3, true)]    // GCCCA: 4 G/C
    [TestCase("TCCTAATTGGAGCGCCCA", 4, false)]
    public void EvaluatePrimer_MaxEndGc_MatchesPrimer3(string primer, int maxEndGc, bool fails)
    {
        var c = PrimerDesigner.EvaluatePrimer(primer, 0, true, P3 with { MaxEndGc = maxEndGc });
        Assert.That(HasIssue(c, "PRIMER_MAX_END_GC"), Is.EqualTo(fails));
    }

    [TestCase("GTTACCGGCCGAGTGCTAC", 3.5, true)]   // END_STABILITY 3.58
    [TestCase("GTTACCGGCCGAGTGCTAC", 6.85, false)]
    [TestCase("CAGTTACCGGCCGAGTGCGC", 3.5, true)]  // 6.09
    [TestCase("AGTTACCGGCCGAGTGCTA", 3.5, false)]  // 3.49
    [TestCase("TCCTAATTGGAGCGCCCA", 3.5, true)]    // 5.36
    [TestCase("TCCTAATTGGAGCGCCCA", 5.36, false)]  // strictly greater fails
    public void EvaluatePrimer_MaxEndStability_MatchesPrimer3(string primer, double max, bool fails)
    {
        var c = PrimerDesigner.EvaluatePrimer(primer, 0, true, P3 with { MaxEndStability = max });
        Assert.That(HasIssue(c, "PRIMER_MAX_END_STABILITY"), Is.EqualTo(fails));
    }

    [Test]
    public void Primer3EndChecks_Defaults_AreInactive_AndMatchPrimer3()
    {
        Assert.Multiple(() =>
        {
            Assert.That(P3.GcClamp, Is.EqualTo(PrimerDesigner.Primer3GcClamp).And.EqualTo(0));
            Assert.That(P3.EffectiveMaxEndGc, Is.EqualTo(5));
            Assert.That(P3.EffectiveMaxEndStability, Is.EqualTo(100.0));
            // Largest Primer3 end_stability of an ACGT pentamer is 6.86 (GCGCG / CGCGC): the default limit 100 and the
            // former −9 kcal/mol gate (Check3PrimeStability) can never reject a primer.
            Assert.That(-PrimerDesigner.Calculate3PrimeStability("GCGCG"), Is.EqualTo(6.86).Within(1e-12));
            var strong = PrimerDesigner.EvaluatePrimer("CAGTTACCGGCCGAGTGCGCG", 0, true,
                PrimerDesigner.DefaultParameters with { Check3PrimeStability = true });
            Assert.That(strong.Issues.Any(i => i.Contains("3' end", StringComparison.Ordinal) && i.Contains("stab", StringComparison.Ordinal)), Is.False);
        });
    }

    [Test]
    public void EvaluatePrimer_Avoid3PrimeGC_LibraryRuleKept()
    {
        // Deprecated library rule (not Primer3): requires ≥ 1 G/C among the last two bases.
        var p = PrimerDesigner.DefaultParameters with { Avoid3PrimeGC = true };
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.EvaluatePrimer("ACGTACGTACGTACGTACAT", 0, true, p).Issues, Has.Some.EqualTo("No GC clamp at 3' end"));
            Assert.That(PrimerDesigner.EvaluatePrimer("ACGTACGTACGTACGTACGA", 0, true, p).Issues, Has.None.EqualTo("No GC clamp at 3' end"));
        });
    }

    [Test]
    public void EvaluatePrimer_IllegalEndChecks_Throw()
    {
        Assert.Multiple(() =>
        {
            // "PRIMER_MAX_END_GC must be between 0 to 5", "PRIMER_MAX_END_STABILITY must be non-negative",
            // "PRIMER_GC_CLAMP > PRIMER_MIN_SIZE" (primer3-py ValueError).
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 0, true, P3 with { MaxEndGc = 6 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 0, true, P3 with { MaxEndGc = -1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 0, true, P3 with { MaxEndStability = -1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.EvaluatePrimer("GTTACCGGCCGAGTGCTAC", 0, true, P3 with { MaxEndStability = double.NaN }));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, P3 with { GcClamp = 19 }, PrimerPairOptions.Primer3Defaults));
            // A negative clamp requires nothing (Primer3 accepts it).
            Assert.That(PrimerDesigner.EvaluatePrimer("AGTTACCGGCCGAGTGCTA", 0, true, P3 with { GcClamp = -1 }).Issues,
                Has.None.Contains("PRIMER_GC_CLAMP"));
        });
    }

    #endregion

    #region DesignPrimerPairs — ranked pairs (primer3-py design_primers, SEQUENCE_TARGET [54, 21], PRIMER_PAIR_MAX_DIFF_TM 100)

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
                Assert.That(pairs[k].Reverse!.Position + pairs[k].Reverse!.Length - 1, Is.EqualTo(e.R), $"rank {k} right 5′ end");
                Assert.That(pairs[k].Reverse!.Length, Is.EqualTo(e.RL), $"rank {k} right length");
                Assert.That(pairs[k].PairPenalty!.Value, Is.EqualTo(e.Pen).Within(1e-9), $"rank {k} PRIMER_PAIR_PENALTY");
                Assert.That(pairs[k].Forward!.Penalty, Is.EqualTo(e.LPen).Within(1e-9), $"rank {k} PRIMER_LEFT_PENALTY");
                Assert.That(pairs[k].Reverse!.Penalty, Is.EqualTo(e.RPen).Within(1e-9), $"rank {k} PRIMER_RIGHT_PENALTY");
            }
        });
    }

    private static IReadOnlyList<PrimerPairResult> Design(PrimerParameters p, PrimerPairOptions o) =>
        PrimerDesigner.DesignPrimerPairs(new DnaSequence(T), 54, 75, p, o);

    [Test]
    public void DesignPrimerPairs_GcClamp2_MatchesPrimer3()
    {
        AssertPairs(Design(P3 with { GcClamp = 2 }, PrimerPairOptions.Primer3Defaults), new[]
        {
            (34, 18, 143, 21, 3.8176034305392363, 2.813907584633739, 1.0036958459054972),
            (18, 18, 143, 21, 3.8176034305392363, 2.813907584633739, 1.0036958459054972),
            (34, 18, 141, 20, 3.893410502139716, 2.813907584633739, 1.0795029175059767),
            (18, 18, 141, 20, 3.893410502139716, 2.813907584633739, 1.0795029175059767),
            (21, 19, 143, 21, 4.058775840031501, 3.055079994126004, 1.0036958459054972),
        });
    }

    [Test]
    public void DesignPrimerPairs_MaxEndGc2_MatchesPrimer3()
    {
        AssertPairs(Design(P3 with { MaxEndGc = 2 }, PrimerPairOptions.Primer3Defaults), new[]
        {
            (19, 19, 136, 20, 2.2130093738348933, 1.6799772818409906, 0.5330320919939027),
            (35, 19, 136, 20, 2.41818879971197, 1.8851567077180675, 0.5330320919939027),
            (18, 20, 136, 20, 2.4739154370097935, 1.9408833450158909, 0.5330320919939027),
            (34, 20, 136, 20, 2.749782019584643, 2.2167499275907403, 0.5330320919939027),
            (19, 19, 135, 21, 2.8825596319761644, 1.6799772818409906, 1.2025823501351738),
        });
    }

    [Test]
    public void DesignPrimerPairs_MaxEndStability_MatchesPrimer3()
    {
        AssertPairs(Design(P3 with { MaxEndStability = 3.5 }, PrimerPairOptions.Primer3Defaults), new[]
        {
            (19, 19, 137, 20, 1.926612109514906, 1.6799772818409906, 0.24663482767391542),
            (35, 19, 137, 20, 2.131791535391983, 1.8851567077180675, 0.24663482767391542),
            (18, 20, 137, 20, 2.1875181726898063, 1.9408833450158909, 0.24663482767391542),
            (19, 19, 136, 20, 2.2130093738348933, 1.6799772818409906, 0.5330320919939027),
            (35, 19, 136, 20, 2.41818879971197, 1.8851567077180675, 0.5330320919939027),
        });
        // PRIMER_MAX_END_STABILITY 3.0: every left primer fails ("high 3' stability"), no pair.
        Assert.That(Design(P3 with { MaxEndStability = 3.0 }, PrimerPairOptions.Primer3Defaults), Is.Empty);
    }

    [Test]
    public void DesignPrimerPairs_MinThreePrimeDistance3_MatchesPrimer3()
    {
        AssertPairs(Design(P3, PrimerPairOptions.Primer3Defaults with { MinThreePrimeDistance = 3 }), new[]
        {
            (20, 19, 137, 20, 1.7645480870700112, 1.5179132593960958, 0.24663482767391542),
            (35, 19, 140, 20, 2.570253724637439, 1.8851567077180675, 0.6850970169193715),
            (17, 18, 135, 21, 3.328255021347104, 2.1256726712119303, 1.2025823501351738),
            (2, 18, 130, 21, 5.547508345073766, 3.734330471554017, 1.8131778735197486),
            (9, 18, 126, 21, 6.252386006043082, 4.077184009505231, 2.175201996537851),
        });
    }

    [Test]
    public void DesignPrimerPairs_MinLeft5Right0_MatchesPrimer3()
    {
        var o = PrimerPairOptions.Primer3Defaults with { MinLeftThreePrimeDistance = 5, MinRightThreePrimeDistance = 0 };
        AssertPairs(Design(P3, o), new[]
        {
            (20, 19, 137, 20, 1.7645480870700112, 1.5179132593960958, 0.24663482767391542),
            (35, 19, 136, 20, 2.41818879971197, 1.8851567077180675, 0.5330320919939027),
            (16, 18, 140, 20, 4.073043969173341, 3.3879469522539694, 0.6850970169193715),
            (2, 18, 138, 20, 4.707962795985452, 3.734330471554017, 0.9736323244314349),
            (9, 18, 143, 21, 5.080879855410728, 4.077184009505231, 1.0036958459054972),
        });
        Assert.That(o.MinThreePrimeDistance, Is.Null, "left ≠ right");
    }

    [Test]
    public void DesignPrimerPairs_MinThreePrimeDistance0_NoPrimerReused_MatchesPrimer3()
    {
        AssertPairs(Design(P3, PrimerPairOptions.Primer3Defaults with { MinThreePrimeDistance = 0 }), new[]
        {
            (20, 19, 137, 20, 1.7645480870700112, 1.5179132593960958, 0.24663482767391542),
            (19, 19, 136, 20, 2.2130093738348933, 1.6799772818409906, 0.5330320919939027),
            (35, 19, 140, 20, 2.570253724637439, 1.8851567077180675, 0.6850970169193715),
            (18, 20, 138, 20, 2.9145156694473258, 1.9408833450158909, 0.9736323244314349),
            (19, 20, 143, 21, 2.945694609380382, 1.941998763474885, 1.0036958459054972),
        });
    }

    [Test]
    public void DesignPrimerPairs_DefaultDistance_ReusesPrimers_MatchesPrimer3()
    {
        // PRIMER_MIN_*_THREE_PRIME_DISTANCE = −1 (default): ranks 0, 1 and 3 share the right primer.
        var o = PrimerPairOptions.Primer3Defaults;
        Assert.Multiple(() =>
        {
            Assert.That(o.MinLeftThreePrimeDistance, Is.EqualTo(-1));
            Assert.That(o.MinRightThreePrimeDistance, Is.EqualTo(-1));
            Assert.That(o.MinThreePrimeDistance, Is.EqualTo(-1));
        });
        AssertPairs(Design(P3, o), new[]
        {
            (20, 19, 137, 20, 1.7645480870700112, 1.5179132593960958, 0.24663482767391542),
            (19, 19, 137, 20, 1.926612109514906, 1.6799772818409906, 0.24663482767391542),
            (20, 19, 136, 20, 2.0509453513899985, 1.5179132593960958, 0.5330320919939027),
            (35, 19, 137, 20, 2.131791535391983, 1.8851567077180675, 0.24663482767391542),
            (18, 20, 137, 20, 2.1875181726898063, 1.9408833450158909, 0.24663482767391542),
        });
    }

    [Test]
    public void DesignPrimerPairs_MinThreePrimeDistanceBelowMinusOne_Throws()
    {
        // "Minimum 3' distance must be >= -1 (min_*_three_prime_distance)".
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => Design(P3, PrimerPairOptions.Primer3Defaults with { MinThreePrimeDistance = -2 }));
            Assert.Throws<ArgumentException>(() => Design(P3, PrimerPairOptions.Primer3Defaults with { MinLeftThreePrimeDistance = -2 }));
            Assert.Throws<ArgumentException>(() => Design(P3, PrimerPairOptions.Primer3Defaults with { MinRightThreePrimeDistance = -5 }));
        });
    }

    #endregion
}
