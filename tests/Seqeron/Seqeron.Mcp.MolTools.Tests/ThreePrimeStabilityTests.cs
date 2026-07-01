using System;
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

            // Sequences shorter than 5 bases return 0.
            Assert.That(MolToolsTools.three_prime_stability("ACGT").DeltaG, Is.EqualTo(0.0).Within(1e-9));
        });
    }
}
