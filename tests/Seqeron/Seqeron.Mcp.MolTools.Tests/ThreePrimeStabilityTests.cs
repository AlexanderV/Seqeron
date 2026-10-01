using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class ThreePrimeStabilityTests
{
    [Test]
    public void ThreePrimeStability_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.three_prime_stability("GCGCG"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.three_prime_stability(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.three_prime_stability(null!));
    }

    [Test]
    public void ThreePrimeStability_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Documented reference values (Primer3 PRIMER_MAX_END_STABILITY):
            //   GCGCG = GC+CG+GC+CG (-2.24-2.17-2.24-2.17) + 0.98 + 0.98 = -6.86.
            Assert.That(MolToolsTools.three_prime_stability("GCGCG").DeltaG, Is.EqualTo(-6.86).Within(1e-9));
            //   TATAT = TA+AT+TA+AT (-0.58-0.88-0.58-0.88) + 1.03 + 1.03 = -0.86.
            Assert.That(MolToolsTools.three_prime_stability("TATAT").DeltaG, Is.EqualTo(-0.86).Within(1e-9));

            // Only the last 5 bases matter: a longer sequence ending GCGCG gives the same value.
            Assert.That(MolToolsTools.three_prime_stability("AAAAAGCGCG").DeltaG, Is.EqualTo(-6.86).Within(1e-9));

            // Shorter primers are scored whole, as Primer3 end_oligodg (oligotm.c compiled from
            // source: end_oligodg("ACGT", 5) = 2.56, i.e. ΔG = −2.56).
            Assert.That(MolToolsTools.three_prime_stability("ACGT").DeltaG, Is.EqualTo(-2.56).Within(1e-9));
            // A character other than A/C/G/T/N in the 3' window is rejected (Primer3 OLIGOTM_ERROR).
            Assert.Throws<ArgumentException>(() => MolToolsTools.three_prime_stability("ACGTU"));
        });
    }
}
