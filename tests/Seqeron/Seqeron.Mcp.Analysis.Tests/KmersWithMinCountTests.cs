using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmers_with_min_count</c> MCP tool.
/// Expected values from the recurrent-k-mer definition (count >= t, sorted desc),
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class KmersWithMinCountTests
{
    [Test]
    public void KmersWithMinCount_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmersWithMinCount("ATGATG", 3, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmersWithMinCount("", 3, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmersWithMinCount(null!, 3, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmersWithMinCount("ATGATG", 0, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmersWithMinCount("ATGATG", -1, 2));
    }

    [Test]
    public void KmersWithMinCount_Binding_InvokesSuccessfully()
    {
        // minCount=2 -> only ATG (count 2) qualifies.
        var recurrent = AnalysisTools.KmersWithMinCount("ATGATG", 3, 2).Items;
        Assert.Multiple(() =>
        {
            Assert.That(recurrent, Has.Length.EqualTo(1));
            Assert.That(recurrent[0].Kmer, Is.EqualTo("ATG"));
            Assert.That(recurrent[0].Count, Is.EqualTo(2));
        });

        // minCount=1 -> all 3 distinct trimers, ATG (count 2) first (sorted desc).
        var all = AnalysisTools.KmersWithMinCount("ATGATG", 3, 1).Items;
        Assert.Multiple(() =>
        {
            Assert.That(all, Has.Length.EqualTo(3));
            Assert.That(all[0].Kmer, Is.EqualTo("ATG"));
            Assert.That(all[0].Count, Is.EqualTo(2));
            Assert.That(all.Select(i => i.Kmer), Is.EquivalentTo(new[] { "ATG", "TGA", "GAT" }));
            Assert.That(all.Where(i => i.Count == 1).Select(i => i.Kmer),
                Is.EquivalentTo(new[] { "TGA", "GAT" }));
        });
    }
}
