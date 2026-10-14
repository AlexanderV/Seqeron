using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Validation 2026-09, B02 (audit A3-13, FIN-B02 WP23): <c>codon_usage_statistics</c> optionally reports
/// CodonW 1.4.4's Fop and CBI (<c>include_fop_cbi</c>, <c>fop_cbi_organism</c> = CodonW <c>-f_type</c>,
/// <c>optimal_codons</c> = CodonW <c>-fop_file</c>/<c>-cbi_file</c>). Values locked against the CodonW
/// binary in Seqeron.Genomics.Tests (CodonUsageAnalyzer_FopCbi_Tests).
/// </summary>
[TestFixture]
public class CodonUsageStatisticsFopCbiTests
{
    // CodonW 1.4.4 -fop -cbi (default E. coli): Fop 0.727 (8/11), CBI 0.561 (23/41).
    private const string Gene = "ATGCTGCTGCTGCTGAAAAAGTTTTTCGGCGGTGGATAA";

    [Test]
    public void Default_FopCbiNull_OutputUnchanged()
    {
        var s = MolToolsTools.codon_usage_statistics(Gene);
        var api = CodonUsageAnalyzer.GetStatistics(Gene, GeneticCode.Standard);
        Assert.Multiple(() =>
        {
            Assert.That(s.Fop, Is.Null);
            Assert.That(s.Cbi, Is.Null);
            Assert.That(s.Enc, Is.EqualTo(api.Enc));
            Assert.That(s.Gc3s, Is.EqualTo(api.Gc3s));
        });
    }

    [Test]
    public void IncludeFopCbi_DefaultEColi_CodonWValues()
    {
        var s = MolToolsTools.codon_usage_statistics(Gene, include_fop_cbi: true);
        Assert.That(s.Fop, Is.EqualTo(8.0 / 11.0).Within(1e-12));
        Assert.That(s.Cbi, Is.EqualTo(23.0 / 41.0).Within(1e-12));
        Assert.That(s.TotalCodons, Is.EqualTo(13));
    }

    [Test]
    public void IncludeFopCbi_OrganismAndCode_CodonWValues()
    {
        // CodonW -f_type 1 -code 1: Fop 0.273, CBI -0.055.
        var s = MolToolsTools.codon_usage_statistics(Gene, 2, include_fop_cbi: true, fop_cbi_organism: "Bacillus subtilis");
        Assert.That(s.Fop, Is.EqualTo(3.0 / 11.0).Within(1e-12));
        Assert.That(s.Cbi, Is.EqualTo(-0.055).Within(0.0005));
    }

    [Test]
    public void IncludeFopCbi_UserOptimalCodons_EqualsApi()
    {
        string[] user = ["cug", "GGC", "AAG"];
        var s = MolToolsTools.codon_usage_statistics(Gene, include_fop_cbi: true, optimal_codons: user);
        var set = new HashSet<string>(user);
        Assert.That(s.Fop, Is.EqualTo(CodonUsageAnalyzer.CalculateFop(Gene, set)));
        Assert.That(s.Cbi, Is.EqualTo(CodonUsageAnalyzer.CalculateCbi(Gene, set)));
        Assert.That(s.Fop, Is.EqualTo(6.0 / 9.0).Within(1e-12), "CTG×4 + AAG + GGC of Leu/Lys/Gly (9 codons)");
    }

    [Test]
    public void IncludeFopCbi_NotCalculable_Null()
    {
        var s = MolToolsTools.codon_usage_statistics("ATGTGGTAA", include_fop_cbi: true);
        Assert.That(s.Fop, Is.Null);
        Assert.That(s.Cbi, Is.Null);
    }

    [Test]
    public void IncludeFopCbi_Serialization_HasFields()
    {
        var web = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        string json = JsonSerializer.Serialize(MolToolsTools.codon_usage_statistics(Gene, include_fop_cbi: true), web);
        Assert.That(json, Does.Contain("\"fop\":0.72"));
        Assert.That(json, Does.Contain("\"cbi\":0.56"));
    }

    [Test]
    public void UnknownOrganism_Throws()
    {
        var ex = Assert.Throws<ArgumentException>(() =>
            MolToolsTools.codon_usage_statistics(Gene, include_fop_cbi: true, fop_cbi_organism: "Homo sapiens"));
        Assert.That(ex!.ParamName, Is.EqualTo("fop_cbi_organism"));
        Assert.That(ex.Message, Does.Contain("EscherichiaColi"));
    }

    [Test]
    public void InvalidOptimalCodon_Throws()
    {
        var ex = Assert.Throws<ArgumentException>(() =>
            MolToolsTools.codon_usage_statistics(Gene, include_fop_cbi: true, optimal_codons: ["CTGA"]));
        Assert.That(ex!.ParamName, Is.EqualTo("optimal_codons"));
    }
}
