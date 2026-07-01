using System;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class PrimerMeltingTemperatureTests
{
    [Test]
    public void PrimerMeltingTemperature_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.primer_melting_temperature("ACGT"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.primer_melting_temperature(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.primer_melting_temperature(null!));
    }

    [Test]
    public void PrimerMeltingTemperature_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Wallace (len 4 < 14): AT=2, GC=2 -> 2*2 + 4*2 = 12.
            Assert.That(MolToolsTools.primer_melting_temperature("ACGT").Tm, Is.EqualTo(12.0).Within(1e-9));

            // Wallace all-GC 6-mer: 4*6 = 24.
            Assert.That(MolToolsTools.primer_melting_temperature("GCGCGC").Tm, Is.EqualTo(24.0).Within(1e-9));

            // Marmur-Doty (len 20 >= 14), GC=20: 64.9 + 41*(20-16.4)/20 = 72.28.
            Assert.That(MolToolsTools.primer_melting_temperature("GCGCGCGCGCGCGCGCGCGC").Tm,
                Is.EqualTo(64.9 + 41.0 * (20 - 16.4) / 20.0).Within(1e-9));
        });
    }
}
