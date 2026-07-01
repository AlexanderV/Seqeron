using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>compare_assemblies</c>. GenomeAssemblyAnalyzer.CompareAssemblies compares k-mer sets:
/// alignedFraction = shared / |kmers|; identity = mean of the two aligned fractions.
/// </summary>
[TestFixture]
public class CompareAssembliesTests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    [Test]
    public void CompareAssemblies_Schema_ValidatesCorrectly()
    {
        var a = new List<NamedSequence> { new("s", Rep("ACGT", 30)) };
        Assert.DoesNotThrow(() => ChromosomeTools.CompareAssemblies(a, a));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.CompareAssemblies(null!, a));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.CompareAssemblies(a, null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.CompareAssemblies(a, a, kmerSize: 0));
    }

    [Test]
    public void CompareAssemblies_IdenticalAssemblies_FullOverlap()
    {
        var a = new List<NamedSequence> { new("s", Rep("ACGT", 30)) };
        var cmp = ChromosomeTools.CompareAssemblies(a, a, "A", "B", 21);

        Assert.Multiple(() =>
        {
            Assert.That(cmp.Assembly1Name, Is.EqualTo("A"));
            Assert.That(cmp.Assembly2Name, Is.EqualTo("B"));
            Assert.That(cmp.AlignedFraction1, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(cmp.AlignedFraction2, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(cmp.SequenceIdentity, Is.EqualTo(1.0).Within(1e-9));
        });
    }
}
