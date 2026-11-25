using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>assess_completeness</c>. GenomeAssemblyAnalyzer.AssessCompleteness aligns nucleotide marker
/// genes to both strands of the assembly and classifies them as compleasm does (ASSEMBLY-COMPLETE-001);
/// percentages are BUSCO's (over n).
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
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AssessCompleteness(asm, markers, minLengthFraction: 2));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AssessCompleteness(asm, new List<NamedSequence> { new("m", "") }));
    }

    [Test]
    public void AssessCompleteness_MarkerPresent_ReportsComplete()
    {
        // The periodic marker also matches its own reverse complement on the same locus: both hits overlap,
        // so compleasm's region clustering keeps it single-copy.
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
            Assert.That(result.BuscoNotation, Is.EqualTo("C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%,n:1"));
        });
    }

    [Test]
    public void AssessCompleteness_DuplicatedAndMissing_BuscoPercentagesOverN()
    {
        // m1 on two scaffolds (duplicated), m2 absent on both strands ⇒ n = 2, D% = 1/2 = 50, M% = 50.
        string m1 = Rep("A", 40) + Rep("C", 40);
        var asm = new List<NamedSequence> { new("s1", m1), new("s2", m1) };
        var markers = new List<NamedSequence> { new("m1", m1), new("m2", Rep("GT", 40)) };

        var result = ChromosomeTools.AssessCompleteness(asm, markers);

        Assert.Multiple(() =>
        {
            Assert.That((result.CompleteSingleCopy, result.CompleteDuplicated, result.Fragmented, result.Missing),
                Is.EqualTo((0, 1, 0, 1)));
            Assert.That(result.DuplicationPercent, Is.EqualTo(50.0).Within(1e-9));
            Assert.That(result.BuscoNotation, Is.EqualTo("C:50.0%[S:0.0%,D:50.0%],F:0.0%,M:50.0%,n:2"));
        });
    }
}
