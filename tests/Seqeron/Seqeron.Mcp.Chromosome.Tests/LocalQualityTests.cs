using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>local_quality</c>. GenomeAssemblyAnalyzer.CalculateLocalQuality slides a window
/// (step = stepSize, default windowSize/2) and reports GC fraction, N count and linguistic complexity per window.
/// </summary>
[TestFixture]
public class LocalQualityTests
{
    [Test]
    public void LocalQuality_Schema_ValidatesCorrectly()
    {
        var seqs = new List<NamedSequence> { new("s", "ACGT") };
        Assert.DoesNotThrow(() => ChromosomeTools.LocalQuality(seqs));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.LocalQuality(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.LocalQuality(seqs, windowSize: 0));
    }

    [Test]
    public void LocalQuality_Binding_InvokesSuccessfully()
    {
        // "GCGCGCGCGC" (10 bp) with windowSize 1000 -> single window, GC = 1.0, N = 0.
        var result = ChromosomeTools.LocalQuality(new List<NamedSequence> { new("s", "GCGCGCGCGC") }, 1000);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var w = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(w.SequenceId, Is.EqualTo("s"));
            Assert.That(w.Position, Is.EqualTo(0));
            Assert.That(w.WindowSize, Is.EqualTo(10));
            Assert.That(w.GcContent, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(w.NCount, Is.EqualTo(0));
        });
    }

    [Test]
    public void LocalQuality_StepSizeAndTruncatedLastWindow()
    {
        var q = ChromosomeTools.LocalQuality(new List<NamedSequence> { new("s", "GGGGAAAACC") }, 4, 3).Items;
        Assert.That(q.Select(w => (w.Position, w.WindowSize)), Is.EqualTo(new[] { (0, 4), (3, 4), (6, 4), (9, 1) }));
        Assert.That(q.Select(w => w.GcContent), Is.EqualTo(new[] { 1.0, 0.25, 0.5, 1.0 }).Within(1e-9));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.LocalQuality(new List<NamedSequence> { new("s", "ACGT") }, 4, 0));
    }

    [Test]
    public void LocalQuality_WindowOfOne_Terminates()
        => Assert.That(ChromosomeTools.LocalQuality(new List<NamedSequence> { new("s", "ACGN") }, 1).Items, Has.Count.EqualTo(4));
}
