using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>extract_contigs</c>. GenomeAssemblyAnalyzer.ExtractContigs splits a scaffold at N-runs
/// and emits each gap-free run of length >= minContigLength as "{id}_contig{n}".
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
}
