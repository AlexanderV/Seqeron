using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class InbreedingFromRohTests
{
    private static RohSegmentItem S(int start, int end) => new(start, end);

    [Test]
    public void InbreedingFromRoh_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.InbreedingFromRoh(
            new[] { S(0, 10_000_000) }, 100_000_000));

        // Non-positive genome length has no denominator → returns 0 (not an error).
        Assert.DoesNotThrow(() => PopulationTools.InbreedingFromRoh(new[] { S(0, 1_000_000) }, 0));
    }

    [Test]
    public void InbreedingFromRoh_Binding_InvokesSuccessfully()
    {
        // F_ROH = (10M + 10M) / 100M = 0.20 (McQuillan et al. 2008).
        var two = PopulationTools.InbreedingFromRoh(
            new[] { S(0, 10_000_000), S(50_000_000, 60_000_000) },
            100_000_000);
        Assert.That(two.InbreedingCoefficient, Is.EqualTo(0.20).Within(1e-10));

        // Whole-genome ROH → F_ROH = 1.
        var whole = PopulationTools.InbreedingFromRoh(new[] { S(0, 2_673_768) }, 2_673_768);
        Assert.That(whole.InbreedingCoefficient, Is.EqualTo(1.0).Within(1e-10));
    }

    [Test]
    public void InbreedingFromRoh_Binding_EdgeCases()
    {
        Assert.Multiple(() =>
        {
            // No segments → 0.
            var none = PopulationTools.InbreedingFromRoh(Array.Empty<RohSegmentItem>(), 300_000_000);
            Assert.That(none.InbreedingCoefficient, Is.EqualTo(0.0).Within(1e-10));

            // Non-positive genome length → 0.
            var badLen = PopulationTools.InbreedingFromRoh(new[] { S(0, 1_000_000) }, 0);
            Assert.That(badLen.InbreedingCoefficient, Is.EqualTo(0.0).Within(1e-10));
        });
    }
}
