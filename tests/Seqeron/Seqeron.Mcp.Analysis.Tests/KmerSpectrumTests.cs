using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_spectrum</c> MCP tool.
/// Expected spectra derived by hand from KmerAnalyzer.GetKmerSpectrum's documented
/// frequency-of-frequencies definition, NOT the wrapper's output.
/// </summary>
[TestFixture]
public class KmerSpectrumTests
{
    [Test]
    public void KmerSpectrum_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmerSpectrum("AAAA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerSpectrum("", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerSpectrum(null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerSpectrum("AAAA", 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerSpectrum("AAAA", -1));
    }

    [Test]
    public void KmerSpectrum_Binding_InvokesSuccessfully()
    {
        // "GTAGAGCTGT" k=1 -> monomer counts G=4,T=3,A=2,C=1 -> {4:1,3:1,2:1,1:1}.
        var mono = AnalysisTools.KmerSpectrum("GTAGAGCTGT", 1).Spectrum;
        Assert.Multiple(() =>
        {
            Assert.That(mono[4], Is.EqualTo(1));
            Assert.That(mono[3], Is.EqualTo(1));
            Assert.That(mono[2], Is.EqualTo(1));
            Assert.That(mono[1], Is.EqualTo(1));
            Assert.That(mono, Has.Count.EqualTo(4));
        });

        // "ATGATG" k=3 -> ATG:2, TGA:1, GAT:1 -> {2:1, 1:2}.
        var tri = AnalysisTools.KmerSpectrum("ATGATG", 3).Spectrum;
        Assert.Multiple(() =>
        {
            Assert.That(tri[2], Is.EqualTo(1));
            Assert.That(tri[1], Is.EqualTo(2));
            Assert.That(tri, Has.Count.EqualTo(2));
        });

        // "AAAA" k=2 -> AA appears 3 times -> {3:1}.
        var homo = AnalysisTools.KmerSpectrum("AAAA", 2).Spectrum;
        Assert.Multiple(() =>
        {
            Assert.That(homo[3], Is.EqualTo(1));
            Assert.That(homo, Has.Count.EqualTo(1));
        });
    }

    [Test]
    public void KmerSpectrum_HistoOptions_MatchJellyfishHisto()
    {
        const string ba1b = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
        // No histo option: histogram stays null (backward-compatible result).
        Assert.That(AnalysisTools.KmerSpectrum(ba1b, 4).Histogram, Is.Null);

        // Jellyfish 2.3.1: count -C -m 4 + histo -> 1 16 2 2 3 1 4 1; histo -f -l 2 -h 6 -> 1 16 2 2 3 1 4 1 5 0 6 0 7 0.
        var canonical = AnalysisTools.KmerSpectrum(ba1b, 4, canonical: true);
        Assert.That(canonical.Spectrum, Is.EquivalentTo(new Dictionary<int, int> { [1] = 16, [2] = 2, [3] = 1, [4] = 1 }));
        var full = AnalysisTools.KmerSpectrum(ba1b, 4, canonical: true, low: 2, high: 6, full: true).Histogram!;
        Assert.That(string.Join(" ", full.Select(r => $"{r.Bin} {r.Frequency}")), Is.EqualTo("1 16 2 2 3 1 4 1 5 0 6 0 7 0"));

        // Homopolymer-rich input, plain count + histo -h 5 -> 1 6 3 2 4 2 6 1 (count 29 pooled in cap bin 6).
        var capped = AnalysisTools.KmerSpectrum("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTACGTACGTACGTTTGCA", 3, acgtOnly: true, high: 5).Histogram!;
        Assert.That(string.Join(" ", capped.Select(r => $"{r.Bin} {r.Frequency}")), Is.EqualTo("1 6 3 2 4 2 6 1"));
    }

    [Test]
    public void KmerSpectrum_InvalidHistoOptions_ThrowArgumentException()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerSpectrum("ACGTACGT", 2, low: 5, high: 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerSpectrum("ACGTACGT", 2, increment: 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerSpectrum("ACGTACGT", 2, low: -1));
    }
}
