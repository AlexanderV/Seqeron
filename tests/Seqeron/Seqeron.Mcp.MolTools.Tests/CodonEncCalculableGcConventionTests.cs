using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Validation 2026-09, B02 F36/F37 (audit A1-11, A1-14):
/// <c>effective_number_of_codons</c> reports <c>calculable</c> (CodonW 1.4.4 prints <c>*****</c>
/// when Nc is not calculated) and <c>codon_usage_statistics</c> reports <c>isEncCalculable</c> and
/// accepts <c>gc_convention</c> (EmbossCusp default = EMBOSS cusp, stops included; CodonW =
/// CodonW <c>gc_out</c>, stop codons of the genetic code excluded).
/// </summary>
[TestFixture]
public class CodonEncCalculableGcConventionTests
{
    private static string Repeat(string s, int n) => string.Concat(Enumerable.Repeat(s, n));

    // Gene M3 of the ENC tests (CodonW -enc: 41.29 under NCBI 1, 45.00 under NCBI 2).
    private static readonly string BiasedGene =
        Repeat("TTT", 4) + "TTC" +
        Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
        Repeat("ATT", 3) + Repeat("ATC", 2) + "ATA" +
        Repeat("GTG", 4) + "GTC" +
        Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
        Repeat("CGC", 4) + Repeat("CGT", 2) +
        Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

    // CodonW -code 0: Nc *****, Len_aa 10, GC 0.600, GC1 0.600, GC2 0.500, GC3 0.700, GC3s 0.625;
    // -code 1 (NCBI 2): Len_aa 11, GC 0.576, GC1 0.545, GC2 0.545, GC3 0.636.
    private const string StopsAmbiguous = "ATGGCTTAAGCCNNNGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA";

    [Test]
    public void EffectiveNumberOfCodons_CodonWStars_CalculableFalse()
    {
        var r = MolToolsTools.effective_number_of_codons("ATGAAAGAGCTGTTCGCCAAA");
        Assert.That(r.Enc, Is.EqualTo(0.0));
        Assert.That(r.Calculable, Is.False);

        var r2 = MolToolsTools.effective_number_of_codons(StopsAmbiguous, 2);
        Assert.That(r2.Calculable, Is.False);
    }

    [TestCase(1, 41.288461538461526)]
    [TestCase(2, 45.0)]
    public void EffectiveNumberOfCodons_CalculableGene_TrueAndCodonWValue(int table, double expected)
    {
        var r = MolToolsTools.effective_number_of_codons(BiasedGene, table);
        Assert.That(r.Calculable, Is.True);
        Assert.That(r.Enc, Is.EqualTo(expected).Within(1e-9));
        Assert.That(new EncResult(r.Enc), Is.EqualTo(r), "back-compatible EncResult(double) ctor derives Calculable");
    }

    [Test]
    public void EffectiveNumberOfCodons_Serialization_HasCalculableField()
    {
        var web = new JsonSerializerOptions(JsonSerializerDefaults.Web);
        string json = JsonSerializer.Serialize(MolToolsTools.effective_number_of_codons("TTTTTTTTTTTC"), web);
        Assert.That(json, Does.Contain("\"calculable\":false"));
        Assert.That(json, Does.Contain("\"enc\":0"));
        string stats = JsonSerializer.Serialize(MolToolsTools.codon_usage_statistics("TTTTTTTTTTTC"), web);
        Assert.That(stats, Does.Contain("\"isEncCalculable\":false"));
    }

    [Test]
    public void CodonUsageStatistics_IsEncCalculable()
    {
        Assert.That(MolToolsTools.codon_usage_statistics(BiasedGene).IsEncCalculable, Is.True);
        Assert.That(MolToolsTools.codon_usage_statistics("ATGTTTAAA").IsEncCalculable, Is.False);
    }

    [Test]
    public void CodonUsageStatistics_DefaultConvention_UnchangedCusp()
    {
        var omitted = MolToolsTools.codon_usage_statistics(StopsAmbiguous);
        var cusp = MolToolsTools.codon_usage_statistics(StopsAmbiguous, 1, CodonGcConvention.EmbossCusp);
        var api = CodonUsageAnalyzer.GetStatistics(StopsAmbiguous, GeneticCode.Standard);

        Assert.Multiple(() =>
        {
            Assert.That(omitted.TotalCodons, Is.EqualTo(14));
            Assert.That(omitted.Gc1, Is.EqualTo(100.0 * 6 / 14).Within(1e-9));
            Assert.That(omitted.OverallGc, Is.EqualTo(100.0 * 20 / 42).Within(1e-9));
            Assert.That(cusp.Gc1, Is.EqualTo(omitted.Gc1));
            Assert.That(cusp.Gc3, Is.EqualTo(omitted.Gc3));
            Assert.That(api.Gc2, Is.EqualTo(omitted.Gc2));
        });
    }

    [TestCase(1, 10, 6, 5, 7, 0.600)]
    [TestCase(2, 11, 6, 6, 7, 0.576)]
    public void CodonUsageStatistics_CodonWConvention_MatchesCodonW(int table, int lenAa, int g1, int g2, int g3, double codonWGc)
    {
        var s = MolToolsTools.codon_usage_statistics(StopsAmbiguous, table, CodonGcConvention.CodonW);
        var api = CodonUsageAnalyzer.GetStatistics(StopsAmbiguous, GeneticCode.GetByTableNumber(table), CodonGcConvention.CodonW);

        Assert.Multiple(() =>
        {
            Assert.That(s.Gc1, Is.EqualTo(100.0 * g1 / lenAa).Within(1e-9));
            Assert.That(s.Gc2, Is.EqualTo(100.0 * g2 / lenAa).Within(1e-9));
            Assert.That(s.Gc3, Is.EqualTo(100.0 * g3 / lenAa).Within(1e-9));
            Assert.That(Math.Round(s.OverallGc / 100, 3), Is.EqualTo(codonWGc).Within(1e-12));
            Assert.That(s.Gc1, Is.EqualTo(api.Gc1));
            Assert.That(s.OverallGc, Is.EqualTo(api.OverallGc));
        });
    }

    [Test]
    public void CodonUsageStatistics_UnknownConvention_Throws()
    {
        var ex = Assert.Throws<ArgumentException>(() =>
            MolToolsTools.codon_usage_statistics("ATG", 1, (CodonGcConvention)9));
        Assert.That(ex!.ParamName, Is.EqualTo("gc_convention"));
    }
}
