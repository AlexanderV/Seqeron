using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class GetEnzymeTests
{
    [Test]
    public void GetEnzyme_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.get_enzyme("EcoRI"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.get_enzyme(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.get_enzyme("   "));
        Assert.Throws<ArgumentException>(() => MolToolsTools.get_enzyme(null!));
    }

    [Test]
    public void GetEnzyme_Binding_InvokesSuccessfully()
    {
        // EcoRI: GAATTC, cut forward 1 / reverse 5, Escherichia coli.
        var ecoRI = MolToolsTools.get_enzyme("EcoRI").Enzyme;
        Assert.Multiple(() =>
        {
            Assert.That(ecoRI, Is.Not.Null);
            Assert.That(ecoRI!.Name, Is.EqualTo("EcoRI"));
            Assert.That(ecoRI.RecognitionSequence, Is.EqualTo("GAATTC"));
            Assert.That(ecoRI.CutPositionForward, Is.EqualTo(1));
            Assert.That(ecoRI.CutPositionReverse, Is.EqualTo(5));
            Assert.That(ecoRI.Organism, Is.EqualTo("Escherichia coli"));
        });

        // Lookup is case-insensitive.
        Assert.That(MolToolsTools.get_enzyme("ecori").Enzyme, Is.Not.Null);

        // Unknown enzyme -> null (documented, not an error).
        Assert.That(MolToolsTools.get_enzyme("NotAnEnzyme").Enzyme, Is.Null);
    }

    /// <summary>
    /// RESTR-FILTER-001 / F28, F33: the lookup also reports the enzyme-level overhang sequence
    /// (REBASE / Biopython <c>ovhgseq</c>) and the same-site enzymes split into equischizomers
    /// (same cut) and neoschizomers (different cut). Values read off Biopython 1.88 this session.
    /// </summary>
    [Test]
    public void GetEnzyme_ReportsOverhangSequenceAndIsoschizomers()
    {
        var ecoRI = MolToolsTools.get_enzyme("EcoRI");
        var pstI = MolToolsTools.get_enzyme("PstI");
        var ecoRV = MolToolsTools.get_enzyme("EcoRV");
        var sfiI = MolToolsTools.get_enzyme("SfiI");
        var mboI = MolToolsTools.get_enzyme("MboI");
        var dpnI = MolToolsTools.get_enzyme("DpnI");
        var unknown = MolToolsTools.get_enzyme("NotAnEnzyme");

        Assert.Multiple(() =>
        {
            // ovhgseq / ovhg from Biopython 1.88.
            Assert.That(ecoRI.OverhangSequence, Is.EqualTo("AATT"));      // ovhg -4, G^AATT_C
            Assert.That(ecoRI.Enzyme!.OverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That(pstI.OverhangSequence, Is.EqualTo("TGCA"));       // ovhg +4, C_TGCA^G
            Assert.That(pstI.Enzyme!.OverhangType, Is.EqualTo(OverhangType.ThreePrime));
            Assert.That(ecoRV.OverhangSequence, Is.Empty);                // ovhg 0, GAT^_ATC
            Assert.That(sfiI.OverhangSequence, Is.EqualTo("NNN"));        // ovhg +3, GGCCN_NNN^NGGCC
            Assert.That(sfiI.Enzyme!.RecognitionLength, Is.EqualTo(13));  // Biopython SfiI.size

            // Same site GATC: MboI/Sau3AI cut 0/4 (equischizomers), DpnI cuts 2/2 (neoschizomer).
            Assert.That(mboI.Equischizomers, Is.EquivalentTo(new[] { "Sau3AI" }));
            Assert.That(mboI.Neoschizomers, Is.EquivalentTo(new[] { "DpnI" }));
            Assert.That(dpnI.Equischizomers, Is.Empty);
            Assert.That(dpnI.Neoschizomers, Is.EquivalentTo(new[] { "MboI", "Sau3AI" }));
            Assert.That(ecoRI.Equischizomers, Is.Empty);
            Assert.That(ecoRI.Neoschizomers, Is.Empty);

            // Unknown name: every added field is null alongside enzyme = null.
            Assert.That(unknown.OverhangSequence, Is.Null);
            Assert.That(unknown.Equischizomers, Is.Null);
            Assert.That(unknown.Neoschizomers, Is.Null);
        });
    }

    /// <summary>
    /// B08 D7 / F52: the enzyme record carries REBASE allenz.104 field &lt;6&gt; verbatim plus its
    /// typed parse — the site the enzyme's COGNATE methyltransferase modifies (EcoRI "3(6)" = m6A at
    /// base 3 of GAATTC; MboI "2(6)"; DpnI has none in REBASE).
    /// </summary>
    [Test]
    public void GetEnzyme_ReportsRebaseCognateMethylationSite()
    {
        var ecoRI = MolToolsTools.get_enzyme("EcoRI").Enzyme!;
        var mboI = MolToolsTools.get_enzyme("MboI").Enzyme!;
        var dpnI = MolToolsTools.get_enzyme("DpnI").Enzyme!;
        var json = System.Text.Json.JsonSerializer.Serialize(
            ecoRI, Seqeron.Mcp.MolTools.AppJsonContext.Default.RestrictionEnzyme);

        Assert.Multiple(() =>
        {
            Assert.That(ecoRI.RebaseMethylationField, Is.EqualTo("3(6)"));
            Assert.That(ecoRI.CognateMethylationSites, Is.EqualTo(new[]
                { new CognateMethylationSite(3, false, MethylatedBase.N6Methyladenine) }));
            Assert.That(mboI.RebaseMethylationField, Is.EqualTo("2(6)"));
            Assert.That(dpnI.RebaseMethylationField, Is.Empty);
            Assert.That(dpnI.CognateMethylationSites, Is.Empty);
            Assert.That(json, Does.Contain("\"RebaseMethylationField\":\"3(6)\"").IgnoreCase);
            Assert.That(json, Does.Contain("CognateMethylationSites").IgnoreCase);
        });
    }
}
