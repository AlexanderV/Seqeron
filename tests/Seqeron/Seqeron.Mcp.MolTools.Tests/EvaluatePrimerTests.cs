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
}
