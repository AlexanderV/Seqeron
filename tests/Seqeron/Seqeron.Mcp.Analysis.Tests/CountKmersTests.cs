using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>count_kmers</c> MCP tool.
/// Expected values derived from KMER-COUNT-001 / the algorithm's own unit tests
/// (KmerAnalyzer_CountKmers_Tests, Wikipedia L-k+1 invariant), NOT the wrapper's output.
/// </summary>
[TestFixture]
public class CountKmersTests
{
    [Test]
    public void CountKmers_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CountKmers("AAAA", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CountKmers("", 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CountKmers(null!, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CountKmers("AAAA", 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CountKmers("AAAA", -3));
    }

    [Test]
    public void CountKmers_Binding_InvokesSuccessfully()
    {
        // "AAAA" k=2 -> AA appears 3 times (= L-k+1).
        var homo = AnalysisTools.CountKmers("AAAA", 2).Counts;
        Assert.Multiple(() =>
        {
            Assert.That(homo["AA"], Is.EqualTo(3));
            Assert.That(homo, Has.Count.EqualTo(1));
        });

        // "ATGATG" k=3 -> ATG:2, TGA:1, GAT:1 (sum = 4 = L-k+1).
        var mixed = AnalysisTools.CountKmers("ATGATG", 3).Counts;
        Assert.Multiple(() =>
        {
            Assert.That(mixed["ATG"], Is.EqualTo(2));
            Assert.That(mixed["TGA"], Is.EqualTo(1));
            Assert.That(mixed["GAT"], Is.EqualTo(1));
            Assert.That(mixed.Values.Sum(), Is.EqualTo(4));
        });
    }
}
