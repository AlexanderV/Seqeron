using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ThermodynamicsTests
{
    [Test]
    public void Thermodynamics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.Thermodynamics("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.Thermodynamics(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.Thermodynamics(null!));
    }

    [Test]
    public void Thermodynamics_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.Thermodynamics("ATGCGATCGATCG");
        Assert.That(result.DeltaH, Is.TypeOf<double>());
        Assert.That(result.DeltaS, Is.TypeOf<double>());
        Assert.That(result.DeltaG, Is.TypeOf<double>());
        Assert.That(result.MeltingTemperature, Is.TypeOf<double>());

        // Longer sequence = higher Tm
        var shortResult = SequenceTools.Thermodynamics("ATGC");
        var longResult = SequenceTools.Thermodynamics("ATGCGATCGATCGATCGATCG");
        Assert.That(longResult.MeltingTemperature, Is.GreaterThan(shortResult.MeltingTemperature));
    }
}
