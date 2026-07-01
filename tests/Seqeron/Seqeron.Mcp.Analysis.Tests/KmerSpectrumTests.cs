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
}
