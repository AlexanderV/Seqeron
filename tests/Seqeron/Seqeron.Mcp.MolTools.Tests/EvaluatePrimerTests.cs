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
            // Marmur-Doty Tm rounded to 1 dp.
            Assert.That(c.MeltingTemperature, Is.EqualTo(Math.Round(64.9 + 41.0 * (10 - 16.4) / 20.0, 1)).Within(1e-9));
            // No homopolymer run beyond 1 in this alternating primer.
            Assert.That(c.HomopolymerLength, Is.EqualTo(1));
            // 3' stability must equal the standalone three_prime_stability of the primer (rounded 1 dp).
            Assert.That(c.Stability3Prime,
                Is.EqualTo(Math.Round(MolToolsTools.three_prime_stability(Primer).DeltaG, 1)).Within(1e-9));
        });
    }
}
