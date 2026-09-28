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
        // OligoCalc salt-adjusted Tm (< 14 nt): Wallace(ACGT) = 12 at the 50 mM reference, shifted by
        // 16.6*log10([Na+]/0.050). 50 mM -> 12.0; 1000 mM -> 12 + 16.6*log10(20) = 33.597 -> 33.6.
        Assert.That(MolToolsTools.primer_melting_temperature_salt("ACGT", 50).Tm, Is.EqualTo(12.0).Within(1e-9));
        Assert.That(MolToolsTools.primer_melting_temperature_salt("ACGT", 1000).Tm, Is.EqualTo(33.6).Within(1e-9));

        // >= 14 nt: 100.5 + 41*GC/N - 820/N + 16.6*log10([Na+]); 20-mer 50 %GC @ 50 mM = 58.403 -> 58.4.
        Assert.That(MolToolsTools.primer_melting_temperature_salt("ACGTACGTACGTACGTACGT", 50).Tm, Is.EqualTo(58.4).Within(1e-9));
    }
}
