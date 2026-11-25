using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_suspicious_regions</c>. GenomeAssemblyAnalyzer.FindSuspiciousRegions flags windows
/// by N content, GC deviation from the pooled GC and SDUST low-complexity intervals (no composite score). A long homopolymer run embedded
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
        var result = ChromosomeTools.FindSuspiciousRegions(new List<NamedSequence> { new("s", seq) }, 0.15, 2.0);

        Assert.That(result.Items, Is.Not.Empty);
        Assert.That(result.Items.Any(r => r.Reason.Contains("Low complexity")), Is.True,
            "The poly-A run must be flagged for low complexity");
    }

    [Test]
    public void FindSuspiciousRegions_PolyA_ReportsMeasuredDustScoreAndNoScore()
    {
        // 600-bp poly(A): one SDUST interval [0,599], DUST = 598*597/2/597 = 299 (Morgulis 2006).
        var r = ChromosomeTools.FindSuspiciousRegions(new List<NamedSequence> { new("s", Rep("A", 600)) }).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That((r.Start, r.End, r.Reason), Is.EqualTo((0, 599, "Low complexity")));
            Assert.That(r.MaxDustScore, Is.EqualTo(299.0).Within(1e-9));
            Assert.That(r.MaxGcDeviation, Is.EqualTo(0.0));
            Assert.That(r.MaxNFraction, Is.EqualTo(0.0));
        });
    }

    [Test]
    public void FindSuspiciousRegions_NRun_IsFlaggedAndParametersAreForwarded()
    {
        var seqs = new List<NamedSequence> { new("s", Rand(6, 1500) + Rep("N", 400) + Rand(7, 1500)) };
        var r = ChromosomeTools.FindSuspiciousRegions(seqs, 0.15, 1e6, 0.1, 500).Items.Single();
        Assert.That((r.Start, r.End, r.Reason), Is.EqualTo((1250, 2249, "High N content")));
        Assert.That(ChromosomeTools.FindSuspiciousRegions(seqs, 0.15, 1e6, 0.9, 500).Items, Is.Empty);
    }
}
