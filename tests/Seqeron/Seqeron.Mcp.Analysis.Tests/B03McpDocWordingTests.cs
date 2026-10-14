using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// FIN-B03 WP16 (B03.md F42, audit row A30): locks the behaviour stated by the reworded
/// dinucleotide_frequencies / codon_frequencies / entropy_profile tool descriptions.
/// </summary>
[TestFixture]
public class B03McpDocWordingTests
{
    // Denominator = in-alphabet dinucleotides (seqinr count(freq = TRUE)), not N − 1: ACNGT has 4 pairs, 2 counted.
    [Test]
    public void DinucleotideFrequencies_DenominatorIsInAlphabetPairs()
    {
        var f = AnalysisTools.DinucleotideFrequencies("ACNGT").Frequencies;
        Assert.That(f, Is.EquivalentTo(new Dictionary<string, double> { ["AC"] = 0.5, ["GT"] = 0.5 }));
    }

    // RNA input: U read as T, codons reported in DNA spelling (EMBOSS cusp).
    [Test]
    public void CodonFrequencies_RnaInput_ReadAsDna()
    {
        var rna = AnalysisTools.CodonFrequencies("AUGAUGUUU").Frequencies;
        var dna = AnalysisTools.CodonFrequencies("ATGATGTTT").Frequencies;
        Assert.Multiple(() =>
        {
            Assert.That(rna, Is.EquivalentTo(dna));
            Assert.That(rna["ATG"], Is.EqualTo(2.0 / 3.0).Within(1e-15));
        });
    }

    // Letter alphabet: N is its own symbol, '-' ignored, T and U distinct → "ACGTN" and "ACTGU" = log2 5, "AC-GT" = 2.
    [Test]
    public void EntropyProfile_LetterAlphabet()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.EntropyProfile("ACGTN", 5, 1).Values[0], Is.EqualTo(Math.Log2(5)).Within(1e-12));
            Assert.That(AnalysisTools.EntropyProfile("ACTGU", 5, 1).Values[0], Is.EqualTo(Math.Log2(5)).Within(1e-12));
            Assert.That(AnalysisTools.EntropyProfile("AC-GT", 5, 1).Values[0], Is.EqualTo(2.0).Within(1e-12));
        });
    }
}
