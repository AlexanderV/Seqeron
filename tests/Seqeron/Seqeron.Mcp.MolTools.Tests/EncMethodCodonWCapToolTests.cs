using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Validation 2026-09, B02 F50/F51 (audit A2-7, A2-8): <c>effective_number_of_codons</c> optional
/// <c>method</c> (Novembre2002 / SunYangXia2013 as in codonbias 0.5.0), <c>pseudocount</c>,
/// <c>background_correction</c>, <c>background</c> and <c>codonw_cap</c> (CodonW 1.4.4 fixed 61).
/// Defaults unchanged (Wright 1990, sense-codon cap).
/// </summary>
[TestFixture]
public class EncMethodCodonWCapToolTests
{
    private static string Repeat(string s, int n) => string.Concat(Enumerable.Repeat(s, n));

    // Gene M3 of the ENC tests (CodonW -enc 41.29 under NCBI 1).
    private static readonly string M3 =
        Repeat("TTT", 4) + "TTC" +
        Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
        Repeat("ATT", 3) + Repeat("ATC", 2) + "ATA" +
        Repeat("GTG", 4) + "GTC" +
        Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
        Repeat("CGC", 4) + Repeat("CGT", 2) +
        Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

    // CODON-ENC-001 M12: codon_i × ((7i mod 5) + 1) — Wright Nc overshoots every cap; CodonW prints 61.00.
    private static readonly string Over = string.Concat(
        (from a in "ACGT" from b in "ACGT" from c in "ACGT" select $"{a}{b}{c}")
        .Select((c, i) => Repeat(c, (7 * i) % 5 + 1)));

    [Test]
    public void Default_IsWright()
    {
        var r = MolToolsTools.effective_number_of_codons(M3);
        Assert.Multiple(() =>
        {
            Assert.That(r, Is.EqualTo(MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.Wright1990, null, null, null, false)));
            Assert.That(r.Enc, Is.EqualTo(41.288461538461526).Within(1e-9));
            Assert.That(r.Calculable, Is.True);
        });
    }

    // codonbias 0.5.0 EffectiveNumberOfCodons(bg_correction=True, robust=False, pseudocount=0, mean='unweighted')
    // and EffectiveNumberOfCodons() (Sun, Yang & Xia defaults).
    [TestCase(EncMethod.Novembre2002, 1, 45.49951134064658)]
    [TestCase(EncMethod.Novembre2002, 2, 46.946215313358586)]
    [TestCase(EncMethod.SunYangXia2013, 1, 50.957809330628805)]
    [TestCase(EncMethod.SunYangXia2013, 2, 50.70834399908201)]
    public void Method_MatchesCodonbiasAndLibrary(EncMethod method, int table, double expected)
    {
        var r = MolToolsTools.effective_number_of_codons(M3, table, method);
        double lib = CodonUsageAnalyzer.CalculateEnc(M3, GeneticCode.GetByTableNumber(table), new EncOptions { Method = method });
        Assert.Multiple(() =>
        {
            Assert.That(r.Enc, Is.EqualTo(expected).Within(1e-9));
            Assert.That(r.Enc, Is.EqualTo(lib));
            Assert.That(r.Calculable, Is.True);
        });
    }

    [Test]
    public void Options_Forwarded()
    {
        Assert.Multiple(() =>
        {
            // codonbias (bg_correction=True) Sun–Yang–Xia: 47.56338039024118.
            Assert.That(MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.SunYangXia2013, background_correction: true).Enc,
                Is.EqualTo(47.56338039024118).Within(1e-9));
            // pseudocount forwarded (mean is not exposed over MCP): equals the library call.
            Assert.That(MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.SunYangXia2013, pseudocount: 0).Enc,
                Is.EqualTo(CodonUsageAnalyzer.CalculateEnc(M3, GeneticCode.Standard,
                    new EncOptions { Method = EncMethod.SunYangXia2013, Pseudocount = 0 })));
            // codonbias get_score(seq, background="GGGCCCGCGATATGCGCGGCCGCAAT"): 38.480446976398305.
            Assert.That(MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.Novembre2002, background: "GGGCCCGCGATATGCGCGGCCGCAAT").Enc,
                Is.EqualTo(38.480446976398305).Within(1e-9));
        });
    }

    [Test]
    public void NoSenseCodon_NotCalculable()
    {
        var r = MolToolsTools.effective_number_of_codons("TAATAG", 1, EncMethod.SunYangXia2013);
        Assert.That(r.Calculable, Is.False);
        Assert.That(r.Enc, Is.EqualTo(0.0));
    }

    // CodonW 1.4.4 -enc -code c prints 61.00 for M12 under NCBI 1, 2, 3, 9.
    [TestCase(1, 61.0, 61.0)]
    [TestCase(2, 60.0, 61.0)]
    [TestCase(3, 62.0, 61.0)]
    [TestCase(9, 62.0, 61.0)]
    public void CodonWCap(int table, double defaultCap, double codonW)
    {
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.effective_number_of_codons(Over, table).Enc, Is.EqualTo(defaultCap));
            Assert.That(MolToolsTools.effective_number_of_codons(Over, table, codonw_cap: true).Enc, Is.EqualTo(codonW));
        });
    }

    [Test]
    public void InvalidCombinations_Rejected()
    {
        void Rejected(TestDelegate call, string param) =>
            Assert.That(call, Throws.ArgumentException.With.Property("ParamName").EqualTo(param));

        Assert.Multiple(() =>
        {
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, (EncMethod)7), "method");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.Wright1990, pseudocount: 1), "method");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.Wright1990, background_correction: true), "method");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.Wright1990, background: "ACGT"), "method");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.Novembre2002, codonw_cap: true), "codonw_cap");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.SunYangXia2013, pseudocount: -0.5), "pseudocount");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.SunYangXia2013, pseudocount: double.NaN), "pseudocount");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.SunYangXia2013, background: "ACGT"), "background");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.Novembre2002, background_correction: false, background: "ACGT"), "background");
            Rejected(() => MolToolsTools.effective_number_of_codons(M3, 1, EncMethod.Novembre2002, background: ""), "background");
        });
    }
}
