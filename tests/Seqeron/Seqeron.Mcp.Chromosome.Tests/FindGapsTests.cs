using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_gaps</c>. GenomeAssemblyAnalyzer.FindGaps reports each N-run of length
/// >= minGapLength as [start, end] inclusive, with a length-class type (5 -> "Short", per ClassifyGap).
/// </summary>
[TestFixture]
public class FindGapsTests
{
    [Test]
    public void FindGaps_Schema_ValidatesCorrectly()
    {
        var seqs = new List<NamedSequence> { new("s", "AAANNNGGG") };
        Assert.DoesNotThrow(() => ChromosomeTools.FindGaps(seqs));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindGaps(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindGaps(seqs, minGapLength: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindGaps(seqs, minGapLength: -1));
    }

    [Test]
    public void FindGaps_Binding_InvokesSuccessfully()
    {
        // "AAANNNNNGGG": one 5-N gap at [3, 7].
        var result = ChromosomeTools.FindGaps(new List<NamedSequence> { new("s", "AAANNNNNGGG") }, 1);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var g = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(g.SequenceId, Is.EqualTo("s"));
            Assert.That(g.Start, Is.EqualTo(3));
            Assert.That(g.End, Is.EqualTo(7));
            Assert.That(g.Length, Is.EqualTo(5));
            Assert.That(g.GapType, Is.EqualTo("Short"));
        });
    }

    [Test]
    public void FindGaps_MinLengthFiltersShortRuns()
    {
        // A 5-N run with minGapLength=10 is filtered out.
        var result = ChromosomeTools.FindGaps(new List<NamedSequence> { new("s", "AAANNNNNGGG") }, 10);
        Assert.That(result.Items, Is.Empty);
    }
}
