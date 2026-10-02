using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class EvaluatePrimerTests
{
    // 20-mer, 50% GC. Marmur-Doty Tm (GC=10, N=20): 64.9 + 41*(10-16.4)/20 = 51.78 -> 51.8.
    private const string Primer = "ATCGATCGATCGATCGATCG";

    [Test]
    public void EvaluatePrimer_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.evaluate_primer(Primer, 0, true));
        Assert.Throws<ArgumentException>(() => MolToolsTools.evaluate_primer("", 0, true));
        Assert.Throws<ArgumentException>(() => MolToolsTools.evaluate_primer(null!, 0, true));
    }

    [Test]
    public void EvaluatePrimer_Binding_InvokesSuccessfully()
    {
        var c = MolToolsTools.evaluate_primer(Primer, 7, is_forward: false);

        Assert.Multiple(() =>
        {
            Assert.That(c.Sequence, Is.EqualTo(Primer));
            Assert.That(c.Position, Is.EqualTo(7));
            Assert.That(c.IsForward, Is.False);
            Assert.That(c.Length, Is.EqualTo(20));
            // GC = 50% (deterministic).
            Assert.That(c.GcContent, Is.EqualTo(50.0).Within(1e-9));
            // Primer3-default Tm (primer3-py 2.3.1 calc_tm = 57.363116…) rounded to 1 dp.
            Assert.That(c.MeltingTemperature, Is.EqualTo(57.4).Within(1e-9));
            // Primer3 per-primer penalty = |Tm − 60| + |20 − 20| = 2.636883760360035.
            Assert.That(c.Penalty, Is.EqualTo(60.0 - 57.363116239639965).Within(1e-9));
            // No homopolymer run beyond 1 in this alternating primer.
            Assert.That(c.HomopolymerLength, Is.EqualTo(1));
            // 3' stability must equal the standalone three_prime_stability of the primer (rounded 1 dp).
            Assert.That(c.Stability3Prime,
                Is.EqualTo(Math.Round(MolToolsTools.three_prime_stability(Primer).DeltaG, 1)).Within(1e-9));
        });
    }

    [Test]
    public void EvaluatePrimer_ThermodynamicStructure_MatchesPrimer3()
    {
        // primer3-py 2.3.1 design_primers(check_primers) / calc_hairpin / calc_homodimer /
        // calc_end_stability at the Primer3 defaults (mv 50, dv 1.5, dNTP 0.6 mM, 50 nM):
        // PRIMER_LEFT_0_HAIRPIN_TH = 70.32453138616256, _SELF_ANY_TH = _SELF_END_TH = 56.93752320069052.
        var c = MolToolsTools.evaluate_primer(Primer, 0, true);
        Assert.Multiple(() =>
        {
            Assert.That(c.HairpinTh, Is.EqualTo(70.32453138616256).Within(1e-9));
            Assert.That(c.SelfAnyTh, Is.EqualTo(56.93752320069052).Within(1e-9));
            Assert.That(c.SelfEndTh, Is.EqualTo(56.93752320069052).Within(1e-9));
            Assert.That(c.HasHairpin, Is.True, "Hairpin Tm exceeds Primer3's PRIMER_MAX_HAIRPIN_TH = 47 °C.");
            Assert.That(c.IsValid, Is.False);
        });
    }
    [Test]
    public void EvaluatePrimer_ReactionConditions_MatchPrimer3()
    {
        // primer3-py 2.3.1 design_primers (PRIMER_SALT_MONOVALENT 20, _SALT_DIVALENT 3.0, PRIMER_DNTP_CONC 0.8,
        // PRIMER_DNA_CONC 250, PRIMER_OPT_GC_PERCENT 45, PRIMER_WT_GC_PERCENT_GT 0.5, _LT 1.0): PRIMER_LEFT_0
        // TACGGGCACTCGTTGGTA TM 61.39018622811568, PENALTY 8.66796400589346, SELF_ANY_TH 4.700639676931019,
        // HAIRPIN_TH 44.063749383794345.
        var c = MolToolsTools.evaluate_primer("TACGGGCACTCGTTGGTA", 36, true,
            Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters,
            salt_monovalent: 20, salt_divalent: 3.0, dntp_conc: 0.8, dna_conc: 250,
            opt_gc_percent: 45, wt_gc_percent_gt: 0.5, wt_gc_percent_lt: 1.0);
        Assert.Multiple(() =>
        {
            Assert.That(c.MeltingTemperature, Is.EqualTo(61.4));
            Assert.That(c.Penalty, Is.EqualTo(8.66796400589346).Within(1e-9));
            Assert.That(c.SelfAnyTh!.Value, Is.EqualTo(4.700639676931019).Within(1e-9));
            Assert.That(c.HairpinTh!.Value, Is.EqualTo(44.063749383794345).Within(1e-9));
        });
        // No condition arguments = Primer3 defaults (penalty |57.363 − 60| on the library defaults, as above).
        Assert.That(MolToolsTools.evaluate_primer(Primer, 0, true).Penalty, Is.EqualTo(60.0 - 57.363116239639965).Within(1e-9));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.evaluate_primer(Primer, 0, true, salt_monovalent: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.evaluate_primer(Primer, 0, true, dntp_conc: -0.5));
    }

    [Test]
    public void EvaluatePrimer_ThreePrimeEndChecks_MatchPrimer3()
    {
        // primer3-py 2.3.1 check_primers GTTACCGGCCGAGTGCTAC: PRIMER_GC_CLAMP 2 → "GC clamp failed 1", PRIMER_MAX_END_GC 2
        // → rejected (TGCTAC: 3 G/C in the last 5), PRIMER_MAX_END_STABILITY 3.5 → "high 3' stability 1" (END_STABILITY 3.58).
        const string primer = "GTTACCGGCCGAGTGCTAC";
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.evaluate_primer(primer, 20, true, p3, gc_clamp: 2).Issues, Has.Some.Contains("PRIMER_GC_CLAMP"));
            Assert.That(MolToolsTools.evaluate_primer(primer, 20, true, p3, gc_clamp: 1).Issues, Has.None.Contains("PRIMER_GC_CLAMP"));
            Assert.That(MolToolsTools.evaluate_primer(primer, 20, true, p3, max_end_gc: 2).Issues, Has.Some.Contains("PRIMER_MAX_END_GC"));
            Assert.That(MolToolsTools.evaluate_primer(primer, 20, true, p3, max_end_gc: 3).Issues, Has.None.Contains("PRIMER_MAX_END_GC"));
            Assert.That(MolToolsTools.evaluate_primer(primer, 20, true, p3, max_end_stability: 3.5).Issues, Has.Some.Contains("PRIMER_MAX_END_STABILITY"));
            Assert.That(MolToolsTools.evaluate_primer(primer, 20, true, p3, max_end_stability: 3.58).Issues, Has.None.Contains("PRIMER_MAX_END_STABILITY"));
        });
    }

    [Test]
    public void EvaluatePrimer_FractionBound_MatchesPrimer3()
    {
        // primer3-py 2.3.1 check_primers AGCTAGCTAGCTAGCTAGCT: PRIMER_ANNEALING_TEMP 60 → PRIMER_LEFT_0_BOUND
        // 34.19287577114419; + PRIMER_WT_BOUND_LT 0.2, _GT 0.3 → PRIMER_LEFT_0_PENALTY 14.460414507505782;
        // PRIMER_MIN_BOUND 40 → "low fraction bound".
        const string primer = "AGCTAGCTAGCTAGCTAGCT";
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters with { MinTm = 0, MaxTm = 100 };
        var c = MolToolsTools.evaluate_primer(primer, 0, true, p3, annealing_temp: 60, wt_bound_lt: 0.2, wt_bound_gt: 0.3);
        Assert.Multiple(() =>
        {
            Assert.That(c.Bound!.Value, Is.EqualTo(34.19287577114419).Within(1e-9));
            Assert.That(c.Penalty, Is.EqualTo(14.460414507505782).Within(1e-9));
            Assert.That(MolToolsTools.evaluate_primer(primer, 0, true, p3).Bound, Is.Null);
            Assert.That(MolToolsTools.evaluate_primer(primer, 0, true, p3, annealing_temp: 60, min_bound: 40).Issues,
                Has.Some.Contains("PRIMER_MIN_BOUND"));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.evaluate_primer(primer, 0, true, p3, opt_bound: 120));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.evaluate_primer(primer, 0, true, p3, annealing_temp: 101));
        });
    }
}
