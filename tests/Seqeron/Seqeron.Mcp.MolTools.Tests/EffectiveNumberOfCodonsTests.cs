using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class EffectiveNumberOfCodonsTests
{
    private static string Repeat(string s, int n) => string.Concat(Enumerable.Repeat(s, n));

    [Test]
    public void EffectiveNumberOfCodons_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.effective_number_of_codons("ATGAAAGAGCTGTTCGCCAAA"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.effective_number_of_codons(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.effective_number_of_codons(null!));
    }

    [Test]
    public void EffectiveNumberOfCodons_Binding_InvokesSuccessfully()
    {
        // Fully-populated biased gene (Genomics test M3, independently verified reference):
        //   F2=0.6, F3=0.2666..., F4=0.4333..., F6=0.3333...
        //   Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6 = 41.288461538461526.
        string biased =
            Repeat("TTT", 4) + "TTC" +
            Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
            Repeat("ATT", 3) + Repeat("ATC", 2) + "ATA" +
            Repeat("GTG", 4) + "GTC" +
            Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
            Repeat("CGC", 4) + Repeat("CGT", 2) +
            Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

        Assert.That(MolToolsTools.effective_number_of_codons(biased).Enc,
            Is.EqualTo(41.288461538461526).Within(1e-9));

        // Short gene whose 3-, 4- and 6-fold classes have no estimable amino acid: CodonW 1.4.4
        // (-enc) prints "*****" (Nc not calculated); the tool returns 0. (Validation 2026-09, F15:
        // previously asserted a value in [20, 61] produced by the non-sourced full-count fallback.)
        var enc = MolToolsTools.effective_number_of_codons("ATGAAAGAGCTGTTCGCCAAA").Enc;
        Assert.That(enc, Is.EqualTo(0.0));
    }
}
