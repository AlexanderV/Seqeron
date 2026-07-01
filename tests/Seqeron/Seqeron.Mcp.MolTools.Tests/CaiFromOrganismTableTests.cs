using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CaiFromOrganismTableTests
{
    private static CodonUsageTableInput PheTable() =>
        new(CodonFrequencies: new Dictionary<string, double> { ["UUU"] = 0.8, ["UUC"] = 0.2 });

    [Test]
    public void CaiFromOrganismTable_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.cai_from_organism_table("TTC", PheTable()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.cai_from_organism_table("", PheTable()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.cai_from_organism_table(null!, PheTable()));
        // Null table -> resolver throws.
        Assert.Throws<ArgumentNullException>(() => MolToolsTools.cai_from_organism_table("TTC", null!));
    }

    [Test]
    public void CaiFromOrganismTable_Binding_InvokesSuccessfully()
    {
        var table = PheTable();

        // Single codon UUC: w = f(UUC)/max(f(UUU),f(UUC)) = 0.2/0.8 = 0.25.
        // CAI = geometric mean of {0.25} = 0.25.
        var single = MolToolsTools.cai_from_organism_table("TTC", table);

        // Two codons UUU (w=1) and UUC (w=0.25): CAI = exp((ln1 + ln0.25)/2) = 0.5.
        var pair = MolToolsTools.cai_from_organism_table("TTTTTC", table);

        Assert.Multiple(() =>
        {
            Assert.That(single.Cai, Is.EqualTo(0.25).Within(1e-9));
            Assert.That(pair.Cai, Is.EqualTo(0.5).Within(1e-9));
        });
    }
}
