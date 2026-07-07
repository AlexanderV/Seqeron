using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_syntenic_blocks_assemblies</c>. GenomeAssemblyAnalyzer.FindSyntenicBlocks anchors
/// shared k-mers between two assemblies and clusters them into blocks, flagging inversions. Identical
/// assemblies yield a forward (non-inverted) block spanning the shared sequence.
/// </summary>
[TestFixture]
public class FindSyntenicBlocksAssembliesTests
{
    private static string Rand(int seed, int len)
    {
        var r = new Random(seed);
        var b = new[] { 'A', 'C', 'G', 'T' };
        var a = new char[len];
        for (int i = 0; i < len; i++) a[i] = b[r.Next(4)];
        return new string(a);
    }

    [Test]
    public void FindSyntenicBlocksAssemblies_Schema_ValidatesCorrectly()
    {
        var a = new List<NamedSequence> { new("s1", Rand(7, 5000)) };
        Assert.DoesNotThrow(() => ChromosomeTools.FindSyntenicBlocksAssemblies(a, a));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindSyntenicBlocksAssemblies(null!, a));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindSyntenicBlocksAssemblies(a, null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindSyntenicBlocksAssemblies(a, a, minBlockSize: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindSyntenicBlocksAssemblies(a, a, kmerSize: 0));
    }

    [Test]
    public void FindSyntenicBlocksAssemblies_IdenticalAssemblies_ForwardBlock()
    {
        var a = new List<NamedSequence> { new("s1", Rand(7, 5000)) };
        var result = ChromosomeTools.FindSyntenicBlocksAssemblies(a, a, 1000, 21);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var b = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(b.Seq1, Is.EqualTo("s1"));
            Assert.That(b.Seq2, Is.EqualTo("s1"));
            Assert.That(b.Start1, Is.EqualTo(0));
            Assert.That(b.IsInverted, Is.False);
        });
    }
}
