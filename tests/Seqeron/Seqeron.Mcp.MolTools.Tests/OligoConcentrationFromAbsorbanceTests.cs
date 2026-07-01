using System;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class OligoConcentrationFromAbsorbanceTests
{
    [Test]
    public void OligoConcentrationFromAbsorbance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.oligo_concentration_from_absorbance(1.0, 10000));
        Assert.Throws<ArgumentException>(() => MolToolsTools.oligo_concentration_from_absorbance(1.0, 0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.oligo_concentration_from_absorbance(1.0, -100));
        Assert.Throws<ArgumentException>(() => MolToolsTools.oligo_concentration_from_absorbance(1.0, 10000, 0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.oligo_concentration_from_absorbance(1.0, 10000, -1));
    }

    [Test]
    public void OligoConcentrationFromAbsorbance_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // c = A / (eps * l) * 1e6. A=1, eps=10000, l=1 -> 1/10000*1e6 = 100 uM.
            Assert.That(MolToolsTools.oligo_concentration_from_absorbance(1.0, 10000, 1.0).ConcentrationMicromolar,
                Is.EqualTo(100.0).Within(1e-9));

            // A=0.5, eps=43000, l=1 -> 0.5/43000*1e6 = 11.62790697... uM.
            Assert.That(MolToolsTools.oligo_concentration_from_absorbance(0.5, 43000).ConcentrationMicromolar,
                Is.EqualTo(0.5 / 43000.0 * 1e6).Within(1e-9));

            // Path length halves concentration: A=1, eps=10000, l=2 -> 50 uM.
            Assert.That(MolToolsTools.oligo_concentration_from_absorbance(1.0, 10000, 2.0).ConcentrationMicromolar,
                Is.EqualTo(50.0).Within(1e-9));
        });
    }
}
