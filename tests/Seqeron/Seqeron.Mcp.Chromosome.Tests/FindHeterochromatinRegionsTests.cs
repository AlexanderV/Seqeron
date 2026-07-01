using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_heterochromatin_regions</c>. ChromosomeAnalyzer.FindHeterochromatinRegions marks
/// windows whose k-mer repeat content >= minRepeatContent and classifies the merged region by position
/// (Telomeric near ends, Centromeric near the middle, else Constitutive).
/// </summary>
[TestFixture]
public class FindHeterochromatinRegionsTests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    [Test]
    public void FindHeterochromatinRegions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.FindHeterochromatinRegions("ACGTACGT"));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindHeterochromatinRegions(""));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindHeterochromatinRegions(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindHeterochromatinRegions("ACGT", windowSize: 0));
    }

    [Test]
    public void FindHeterochromatinRegions_Binding_InvokesSuccessfully()
    {
        // 300 kb uniform 'A' is maximally repetitive: one region spanning [0, len-1], Constitutive.
        var result = ChromosomeTools.FindHeterochromatinRegions(Rep("A", 300000), 100000, 0.5);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var r = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(r.Start, Is.EqualTo(0));
            Assert.That(r.End, Is.EqualTo(299999));
            Assert.That(r.Type, Is.EqualTo("Constitutive"));
        });
    }

    [Test]
    public void FindHeterochromatinRegions_NonRepetitive_ReturnsNone()
    {
        // A short non-repetitive sequence shorter than a window yields no regions.
        var result = ChromosomeTools.FindHeterochromatinRegions("ACGTACGTACGT", 100000, 0.5);
        Assert.That(result.Items, Is.Empty);
    }
}
