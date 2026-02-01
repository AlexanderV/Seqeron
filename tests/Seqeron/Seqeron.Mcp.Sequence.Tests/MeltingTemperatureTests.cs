using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class MeltingTemperatureTests
{
    [Test]
    public void MeltingTemperature_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.MeltingTemperature("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.MeltingTemperature(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.MeltingTemperature(null!));
    }

    [Test]
    public void MeltingTemperature_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.MeltingTemperature("ATGCGATCGATCG");
        Assert.That(result.Tm, Is.TypeOf<double>());
        Assert.That(result.Unit, Is.EqualTo("°C"));

        // GC-rich sequence has higher Tm
        var atRichResult = SequenceTools.MeltingTemperature("ATATATATAT");
        var gcRichResult = SequenceTools.MeltingTemperature("GCGCGCGCGC");
        Assert.That(gcRichResult.Tm, Is.GreaterThan(atRichResult.Tm));
    }
}
