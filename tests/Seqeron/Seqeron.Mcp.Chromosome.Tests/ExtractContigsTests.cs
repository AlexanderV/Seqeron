using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>extract_contigs</c>. GenomeAssemblyAnalyzer.ExtractContigs (QUAST broken scaffolds) splits
/// a scaffold at N-runs of length >= minGapLength (default 10) and emits each piece of length
/// >= minContigLength as "{id}_contig{n}".
/// </summary>
[TestFixture]
public class ExtractContigsTests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    [Test]
    public void ExtractContigs_Schema_ValidatesCorrectly()
    {
        var scaffolds = new List<NamedSequence> { new("s", Rep("A", 300)) };
        Assert.DoesNotThrow(() => ChromosomeTools.ExtractContigs(scaffolds));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ExtractContigs(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ExtractContigs(scaffolds, minContigLength: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.ExtractContigs(scaffolds, minGapLength: 0));
    }

    [Test]
    public void ExtractContigs_Binding_InvokesSuccessfully()
    {
        // 300 A + 20 N + 250 C -> two contigs (both >= 200 bp).
        var scaffolds = new List<NamedSequence> { new("s", Rep("A", 300) + Rep("N", 20) + Rep("C", 250)) };
        var result = ChromosomeTools.ExtractContigs(scaffolds, 200);

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(2));
            Assert.That(result.Items[0].Id, Is.EqualTo("s_contig1"));
            Assert.That(result.Items[0].Sequence.Length, Is.EqualTo(300));
            Assert.That(result.Items[1].Id, Is.EqualTo("s_contig2"));
            Assert.That(result.Items[1].Sequence.Length, Is.EqualTo(250));
        });
    }

    [Test]
    public void ExtractContigs_MinLengthFiltersShortRuns()
    {
        // Second run (100 C) is below the 200 threshold and is dropped.
        var scaffolds = new List<NamedSequence> { new("s", Rep("A", 300) + Rep("N", 20) + Rep("C", 100)) };
        var result = ChromosomeTools.ExtractContigs(scaffolds, 200);
        Assert.That(result.Items, Has.Count.EqualTo(1));
    }

    [Test]
    public void ExtractContigs_MinGapLength_DelegatesQuastBreakThreshold()
    {
        // 9 Ns < default threshold 10: one contig with the Ns inside; minGapLength 1 cuts at every N.
        string seq = Rep("A", 300) + Rep("N", 9) + Rep("C", 300);
        var scaffolds = new List<NamedSequence> { new("s", seq) };

        Assert.That(ChromosomeTools.ExtractContigs(scaffolds, 200).Items.Select(i => i.Sequence),
            Is.EqualTo(new[] { seq }));
        Assert.That(ChromosomeTools.ExtractContigs(scaffolds, 200, minGapLength: 1).Items.Select(i => i.Sequence.Length),
            Is.EqualTo(new[] { 300, 300 }));
    }
}
