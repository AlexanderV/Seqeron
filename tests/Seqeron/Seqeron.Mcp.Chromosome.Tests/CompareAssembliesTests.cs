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
            Assert.That(cmp.Containment1, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(cmp.Containment2, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(cmp.JaccardIndex, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(cmp.MashDistance, Is.EqualTo(0.0).Within(1e-9));
        });
    }

    [Test]
    public void CompareAssemblies_PartialOverlap_ReportsMashValues()
    {
        // canonical 4-mers {ACGT,CGTA,GTAC} vs {ACGT,CGTA,GTAA,TAAA,AAAA}: J = 2/6, D = ln2/4 (Ondov 2016 eq. 4).
        var cmp = ChromosomeTools.CompareAssemblies(
            new List<NamedSequence> { new("s", "ACGTACGT") }, new List<NamedSequence> { new("s", "ACGTAAAA") }, "x", "y", 4);
        Assert.Multiple(() =>
        {
            Assert.That(cmp.SharedKmers, Is.EqualTo(2));
            Assert.That(cmp.JaccardIndex, Is.EqualTo(1.0 / 3).Within(1e-12));
            Assert.That(cmp.MashDistance, Is.EqualTo(Math.Log(2) / 4).Within(1e-12));
        });
    }
}
