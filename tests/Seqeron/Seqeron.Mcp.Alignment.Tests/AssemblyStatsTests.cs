using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class AssemblyStatsTests
{
    [Test]
    public void AssemblyStats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.AssemblyStats(new[] { "AAAA" }, 1));
        var empty = AlignmentTools.AssemblyStats(Array.Empty<string>(), 3);
        Assert.That(empty.Contigs, Is.Empty);
        Assert.That(empty.TotalReads, Is.EqualTo(3));
        Assert.That(empty.TotalLength, Is.EqualTo(0));
    }

    [Test]
    public void AssemblyStats_Binding_InvokesSuccessfully()
    {
        // Contig lengths 10, 20, 30 (total 60, half 30). Sorted desc: 30, 20, 10.
        // Cumulative reaches 30 (>= half) at the first contig -> N50 = 30, longest = 30.
        var contigs = new[] { new string('A', 10), new string('A', 20), new string('A', 30) };
        var r = AlignmentTools.AssemblyStats(contigs, 5);
        Assert.Multiple(() =>
        {
            Assert.That(r.TotalReads, Is.EqualTo(5));
            Assert.That(r.LongestContig, Is.EqualTo(30));
            Assert.That(r.TotalLength, Is.EqualTo(60));
            Assert.That(r.N50, Is.EqualTo(30));
        });
    }

    [Test]
    public void AssemblyStats_OddTotal_N50MatchesQuast_AndAssembledReadsIsOptional()
    {
        // QUAST N50([4,3,2]) = 3 (4 covers 4/9 < 50 %); the old total/2 = 4 cut-off gave 4.
        var contigs = new[] { "AAAA", "CCC", "GG" };
        var byDefault = AlignmentTools.AssemblyStats(contigs, 6);
        var explicitCount = AlignmentTools.AssemblyStats(contigs, 6, 5);
        Assert.Multiple(() =>
        {
            Assert.That(byDefault.N50, Is.EqualTo(3));
            Assert.That(byDefault.AssembledReads, Is.EqualTo(6));
            Assert.That(explicitCount.AssembledReads, Is.EqualTo(5));
            Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.AssemblyStats(contigs, 6, 7));
        });
    }
}
