using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for the <c>analyze_centromere</c> MCP tool.
/// Expected values come from ChromosomeAnalyzer.AnalyzeCentromere behavior documented in
/// Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Centromere_Tests.cs (CHROM-CENT-001,
/// Levan et al. 1964) — NOT from the wrapper's own output.
/// </summary>
[TestFixture]
public class AnalyzeCentromereTests
{
    [Test]
    public void AnalyzeCentromere_Schema_ValidatesCorrectly()
    {
        // Valid call must not throw.
        Assert.DoesNotThrow(() =>
            ChromosomeTools.AnalyzeCentromere("chr1", new string('A', 300000), 10000, 0.3));

        // Guards.
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeCentromere("chr1", ""));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeCentromere("chr1", null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeCentromere("", "ACGTACGT"));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeCentromere(null!, "ACGTACGT"));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeCentromere("chr1", "ACGTACGT", windowSize: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeCentromere("chr1", "ACGTACGT", minAlphaSatelliteContent: 1.5));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeCentromere("chr1", "ACGTACGT", minAlphaSatelliteContent: -0.1));
    }

    [Test]
    public void AnalyzeCentromere_Binding_InvokesSuccessfully()
    {
        // Uniform 300 kb 'A' run is maximally repetitive: the algorithm anchors the first
        // window (start=0), extends right by windowSize/2 until end-5000, and classifies by
        // arm ratio. These exact values are deterministic for this input.
        var result = ChromosomeTools.AnalyzeCentromere(
            "chr1", new string('A', 300000), windowSize: 10000, minAlphaSatelliteContent: 0.3);

        Assert.Multiple(() =>
        {
            Assert.That(result.Chromosome, Is.EqualTo("chr1"));
            Assert.That(result.Start, Is.EqualTo(0));
            Assert.That(result.End, Is.EqualTo(295000));
            Assert.That(result.Length, Is.EqualTo(295000));
            Assert.That(result.Length, Is.EqualTo(result.End!.Value - result.Start!.Value));
            Assert.That(result.CentromereType, Is.EqualTo("Metacentric"));
            Assert.That(result.AlphaSatelliteContent, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(result.IsAcrocentric, Is.False);
        });
    }

    [Test]
    public void AnalyzeCentromere_NonRepetitiveSequence_ReturnsUnknown()
    {
        // Random, non-repetitive input: no centromere detected (mirrors
        // AnalyzeCentromere_WithNoRepetitiveRegions_ReturnsUnknown).
        var random = new Random(12345);
        var bases = new[] { 'A', 'C', 'G', 'T' };
        var buf = new char[500000];
        for (int i = 0; i < buf.Length; i++) buf[i] = bases[random.Next(4)];
        var sequence = new string(buf);

        var result = ChromosomeTools.AnalyzeCentromere(
            "chr1", sequence, windowSize: 50000, minAlphaSatelliteContent: 0.3);

        Assert.Multiple(() =>
        {
            Assert.That(result.CentromereType, Is.EqualTo("Unknown"));
            Assert.That(result.Start, Is.Null);
            Assert.That(result.End, Is.Null);
            Assert.That(result.Length, Is.EqualTo(0));
            Assert.That(result.AlphaSatelliteContent, Is.LessThan(0.3));
        });
    }
}
