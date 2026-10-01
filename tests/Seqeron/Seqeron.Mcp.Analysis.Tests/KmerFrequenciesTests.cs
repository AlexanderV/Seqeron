using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_frequencies</c> MCP tool.
/// Expected frequencies computed by hand from count/(L-k+1), NOT the wrapper output.
/// </summary>
[TestFixture]
public class KmerFrequenciesTests
{
    [Test]
    public void KmerFrequencies_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmerFrequencies("AAAA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerFrequencies("", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerFrequencies(null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerFrequencies("AAAA", 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerFrequencies("AAAA", -1));
    }

    [Test]
    public void KmerFrequencies_Binding_InvokesSuccessfully()
    {
        // "ATGC" k=1 -> each base 1/4 = 0.25, sum = 1.
        var mono = AnalysisTools.KmerFrequencies("ATGC", 1).Frequencies;
        Assert.Multiple(() =>
        {
            Assert.That(mono["A"], Is.EqualTo(0.25).Within(1e-12));
            Assert.That(mono["T"], Is.EqualTo(0.25).Within(1e-12));
            Assert.That(mono["G"], Is.EqualTo(0.25).Within(1e-12));
            Assert.That(mono["C"], Is.EqualTo(0.25).Within(1e-12));
            Assert.That(mono.Values.Sum(), Is.EqualTo(1.0).Within(1e-12));
        });

        // "AAAA" k=2 -> AA = 3/3 = 1.0, single entry.
        var homo = AnalysisTools.KmerFrequencies("AAAA", 2).Frequencies;
        Assert.Multiple(() =>
        {
            Assert.That(homo["AA"], Is.EqualTo(1.0).Within(1e-12));
            Assert.That(homo, Has.Count.EqualTo(1));
        });
    }

    [Test]
    public void KmerFrequencies_CanonicalAndAcgtOnly_EqualJellyfishAndKpalProfiles()
    {
        // jellyfish count -m 3 -C + dump -c on ACGTNACGTAAcgtRTT: ACG 6, AAC 1, GTA 1, TAA 1 (sum 9).
        var canonical = AnalysisTools.KmerFrequencies("ACGTNACGTAAcgtRTT", 3, canonical: true).Frequencies;
        Assert.That(canonical, Has.Count.EqualTo(4));
        Assert.That(canonical["ACG"], Is.EqualTo(6.0 / 9).Within(1e-15));
        Assert.That(canonical["TAA"], Is.EqualTo(1.0 / 9).Within(1e-15));
        // kPAL profile (= jellyfish without -C) of AAAANTTTTGGGGuCCCC k=3: AAA CCC GGG TTT 2, TGG TTG 1 (sum 10).
        var acgt = AnalysisTools.KmerFrequencies("AAAANTTTTGGGGuCCCC", 3, acgtOnly: true).Frequencies;
        Assert.That(acgt, Has.Count.EqualTo(6));
        Assert.That(acgt["AAA"], Is.EqualTo(0.2).Within(1e-15));
        Assert.That(acgt["TGG"], Is.EqualTo(0.1).Within(1e-15));
    }
}
