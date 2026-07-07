using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_suspicious_regions</c>. GenomeAssemblyAnalyzer.FindSuspiciousRegions flags windows
/// by GC deviation from the global GC and by low linguistic complexity. A long homopolymer run embedded
/// in normal sequence is flagged for low complexity.
/// </summary>
[TestFixture]
public class FindSuspiciousRegionsTests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    private static string Rand(int seed, int len)
    {
        var r = new Random(seed);
        var b = new[] { 'A', 'C', 'G', 'T' };
        var a = new char[len];
        for (int i = 0; i < len; i++) a[i] = b[r.Next(4)];
        return new string(a);
    }

    [Test]
    public void FindSuspiciousRegions_Schema_ValidatesCorrectly()
    {
        var seqs = new List<NamedSequence> { new("s", Rand(1, 1000)) };
        Assert.DoesNotThrow(() => ChromosomeTools.FindSuspiciousRegions(seqs));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindSuspiciousRegions(null!));
    }

    [Test]
    public void FindSuspiciousRegions_LowComplexityRun_IsFlagged()
    {
        // 2 kb random | 2 kb poly-A | 2 kb random: the middle is a low-complexity misassembly signal.
        var seq = Rand(1, 2000) + Rep("A", 2000) + Rand(2, 2000);
        var result = ChromosomeTools.FindSuspiciousRegions(new List<NamedSequence> { new("s", seq) }, 0.15, 0.3);

        Assert.That(result.Items, Is.Not.Empty);
        Assert.That(result.Items.Any(r => r.Reason.Contains("Low complexity")), Is.True,
            "The poly-A run must be flagged for low complexity");
    }
}
