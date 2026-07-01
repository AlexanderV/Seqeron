using System;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class PrimerMeltingTemperatureSaltTests
{
    [Test]
    public void PrimerMeltingTemperatureSalt_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.primer_melting_temperature_salt("ACGT"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.primer_melting_temperature_salt(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.primer_melting_temperature_salt(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.primer_melting_temperature_salt("ACGT", 0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.primer_melting_temperature_salt("ACGT", -5));
    }

    [Test]
    public void PrimerMeltingTemperatureSalt_Binding_InvokesSuccessfully()
    {
        // base Tm(ACGT) = 12. Salt correction at 50 mM = 16.6*log10(0.05) = -21.597...
        // 12 + (-21.597) = -9.597 -> round(1) = -9.6.
        var expected50 = Math.Round(12.0 + 16.6 * Math.Log10(50.0 / 1000.0), 1);
        Assert.That(MolToolsTools.primer_melting_temperature_salt("ACGT", 50).Tm,
            Is.EqualTo(expected50).Within(1e-9));
        Assert.That(expected50, Is.EqualTo(-9.6).Within(1e-9));

        // At [Na+] = 1000 mM the correction is 16.6*log10(1) = 0, so salt Tm = round(baseTm) = 12.0.
        Assert.That(MolToolsTools.primer_melting_temperature_salt("ACGT", 1000).Tm,
            Is.EqualTo(12.0).Within(1e-9));
    }
}
