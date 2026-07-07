using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class OligoExtinctionCoefficientTests
{
    [Test]
    public void OligoExtinctionCoefficient_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.oligo_extinction_coefficient("ACGT"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.oligo_extinction_coefficient(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.oligo_extinction_coefficient(null!));
    }

    [Test]
    public void OligoExtinctionCoefficient_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // A+C+G+T = 15400 + 7400 + 11500 + 8700 = 43000.
            Assert.That(MolToolsTools.oligo_extinction_coefficient("ACGT").ExtinctionCoefficient,
                Is.EqualTo(43000).Within(1e-9));

            // Single U = 9900.
            Assert.That(MolToolsTools.oligo_extinction_coefficient("U").ExtinctionCoefficient,
                Is.EqualTo(9900).Within(1e-9));

            // Unknown base falls back to 10000.
            Assert.That(MolToolsTools.oligo_extinction_coefficient("N").ExtinctionCoefficient,
                Is.EqualTo(10000).Within(1e-9));
        });
    }
}
