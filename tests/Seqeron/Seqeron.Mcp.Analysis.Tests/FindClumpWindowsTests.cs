using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_clump_windows</c> MCP tool. Expected values: a Python brute force that counts every
/// window Genome[i..i+L-1] (occurrence at p counts when i &lt;= p &lt;= i+L-k; Compeau &amp; Pevzner ch. 1,
/// Rosalind BA1E) and compresses the qualifying i of each k-mer into maximal runs (K-mer_Search.md §7.3).
/// NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindClumpWindowsTests
{
    private const string Ba1e =
        "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC";

    [Test]
    public void FindClumpWindows_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindClumpWindows("AAAA", 2, 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumpWindows("", 2, 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumpWindows(null!, 2, 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumpWindows("AAAA", 0, 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumpWindows("AAAA", 3, 2, 3)); // window < k
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumpWindows("AAAA", 2, 4, 0)); // t <= 0
    }

    [Test]
    public void FindClumpWindows_Binding_MatchesBruteForce()
    {
        // Rosalind BA1E sample (k=5, L=75, t=4): brute force CGACA [0,6], GAAGA [0,16], AATGT [16,21].
        var ba1e = AnalysisTools.FindClumpWindows(Ba1e, 5, 75, 4).Clumps;
        Assert.That(ba1e.Select(c => c.Kmer), Is.EqualTo(new[] { "CGACA", "GAAGA", "AATGT" }));
        Assert.That(ba1e.Select(c => c.FirstWindowStart), Is.EqualTo(new[] { 0, 0, 16 }));
        Assert.That(ba1e[1].WindowRuns, Is.EqualTo(new[] { new ClumpWindowRunItem(0, 16) }));
        Assert.That(ba1e[2].WindowRuns, Is.EqualTo(new[] { new ClumpWindowRunItem(16, 21) }));

        // BA1B sample, k=4, L=12, t=2: split runs GCAT [4,5],[11,12]; CATG [5,6],[12,13]; ATGA [13,14].
        var split = AnalysisTools.FindClumpWindows("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 12, 2).Clumps;
        Assert.That(split.Select(c => c.Kmer), Is.EqualTo(new[] { "GCAT", "CATG", "ATGA" }));
        Assert.That(split[0].WindowRuns, Is.EqualTo(new[] { new ClumpWindowRunItem(4, 5), new ClumpWindowRunItem(11, 12) }));
        Assert.That(split[1].WindowRuns, Is.EqualTo(new[] { new ClumpWindowRunItem(5, 6), new ClumpWindowRunItem(12, 13) }));

        // Same k-mer set as find_clumps.
        Assert.That(ba1e.Select(c => c.Kmer), Is.EquivalentTo(AnalysisTools.FindClumps(Ba1e, 5, 75, 4).Kmers));

        // No clump: empty list.
        Assert.That(AnalysisTools.FindClumpWindows("AAAA", 2, 4, 4).Clumps, Is.Empty);
    }
}
