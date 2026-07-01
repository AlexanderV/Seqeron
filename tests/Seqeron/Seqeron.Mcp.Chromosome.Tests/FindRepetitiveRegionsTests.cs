using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_repetitive_regions</c>. GenomeAssemblyAnalyzer.FindRepetitiveRegions flags spans
/// where k-mers recur at least minCopies times. A pure 8-bp periodic sequence collapses to one region.
/// </summary>
[TestFixture]
public class FindRepetitiveRegionsTests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    [Test]
    public void FindRepetitiveRegions_Schema_ValidatesCorrectly()
    {
        var seqs = new List<NamedSequence> { new("s", Rep("ACGTACGT", 100)) };
        Assert.DoesNotThrow(() => ChromosomeTools.FindRepetitiveRegions(seqs));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindRepetitiveRegions(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindRepetitiveRegions(seqs, kmerSize: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindRepetitiveRegions(seqs, minCopies: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindRepetitiveRegions(seqs, windowSize: 0));
    }

    [Test]
    public void FindRepetitiveRegions_Binding_InvokesSuccessfully()
    {
        // "ACGTACGT" x 100 = 800 bp: an 8-mer recurs ~199 times; one repetitive region [0, 789].
        var result = ChromosomeTools.FindRepetitiveRegions(
            new List<NamedSequence> { new("s", Rep("ACGTACGT", 100)) }, 8, 3, 50);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var r = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceId, Is.EqualTo("s"));
            Assert.That(r.Start, Is.EqualTo(0));
            Assert.That(r.End, Is.EqualTo(789));
            Assert.That(r.Copies, Is.EqualTo(199));
        });
    }
}
