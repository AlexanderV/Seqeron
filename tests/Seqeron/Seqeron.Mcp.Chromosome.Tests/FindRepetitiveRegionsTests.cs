using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_repetitive_regions</c>. GenomeAssemblyAnalyzer.FindRepetitiveRegions flags spans
/// runs the NCBI WindowMasker k-mer-count model. A pure 8-bp periodic sequence collapses to one region.
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
        // "ACGTACGT" x 100 = 800 bp: one WindowMasker region [0, 799]. ACGTACGT is its own reverse complement and
        // occurs 199 times per strand, so its both-strand unit count (win_mask_gen_counts.cpp) is 398.
        var result = ChromosomeTools.FindRepetitiveRegions(
            new List<NamedSequence> { new("s", Rep("ACGTACGT", 100)) }, 8, 3, 50);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var r = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceId, Is.EqualTo("s"));
            Assert.That(r.Start, Is.EqualTo(0));
            Assert.That(r.End, Is.EqualTo(799));
            Assert.That(r.Copies, Is.EqualTo(398));
        });
    }

    [Test]
    public void FindRepetitiveRegions_StrandAndExtendOptions_AreForwarded()
    {
        // k = W = 5: ACGTG twice forward + once as its reverse complement CACGT.
        var seqs = new List<NamedSequence> { new("a", "ACGTG"), new("b", "ACGTG"), new("c", "CACGT") };

        Assert.That(ChromosomeTools.FindRepetitiveRegions(seqs, 5, 3, 5, bothStrands: false).Items, Is.Empty);
        var both = ChromosomeTools.FindRepetitiveRegions(seqs, 5, 3, 5, bothStrands: true).Items;
        Assert.That(both.Select(r => (r.SequenceId, r.Start, r.End, r.Copies)),
            Is.EqualTo(new[] { ("a", 0, 4, 3), ("b", 0, 4, 3), ("c", 0, 4, 3) }));

        Assert.Throws<ArgumentOutOfRangeException>(() => ChromosomeTools.FindRepetitiveRegions(seqs, 5, 3, 5, true, extendCopies: 4));
    }
}
