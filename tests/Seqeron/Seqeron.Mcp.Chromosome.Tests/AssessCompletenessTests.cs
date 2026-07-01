using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>assess_completeness</c>. GenomeAssemblyAnalyzer.AssessCompleteness aligns marker genes
/// to the assembly (BUSCO-like). A single marker fully contained in the assembly is complete/single-copy.
/// </summary>
[TestFixture]
public class AssessCompletenessTests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    [Test]
    public void AssessCompleteness_Schema_ValidatesCorrectly()
    {
        var asm = new List<NamedSequence> { new("scaf", Rep("ACGTTGCA", 50)) };
        var markers = new List<NamedSequence> { new("m1", Rep("ACGTTGCA", 50)) };
        Assert.DoesNotThrow(() => ChromosomeTools.AssessCompleteness(asm, markers));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AssessCompleteness(null!, markers));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AssessCompleteness(asm, null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AssessCompleteness(asm, markers, identityThreshold: 1.5));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AssessCompleteness(asm, markers, coverageThreshold: -0.1));
    }

    [Test]
    public void AssessCompleteness_MarkerPresent_ReportsComplete()
    {
        var asm = new List<NamedSequence> { new("scaf", Rep("ACGTTGCA", 50)) };
        var markers = new List<NamedSequence> { new("m1", Rep("ACGTTGCA", 50)) };

        var result = ChromosomeTools.AssessCompleteness(asm, markers, 0.9, 0.9);

        Assert.Multiple(() =>
        {
            Assert.That(result.TotalGenes, Is.EqualTo(1));
            Assert.That(result.Complete, Is.EqualTo(1));
            Assert.That(result.CompleteSingleCopy, Is.EqualTo(1));
            Assert.That(result.Missing, Is.EqualTo(0));
            Assert.That(result.CompletenessPercent, Is.EqualTo(100.0).Within(1e-9));
        });
    }
}
