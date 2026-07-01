using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_syntenic_blocks</c> MCP tool.
/// Expected values from ComparativeGenomics' own unit test
/// (ComparativeGenomics_FindSyntenicBlocks_Tests.FiveAdjacentForwardAnchors: 5 collinear
/// gi->hi anchors form one forward block of 5), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindSyntenicBlocksTests
{
    private static GeneInput[] Genome(string prefix, string genomeId) =>
        Enumerable.Range(0, 5)
            .Select(i => new GeneInput($"{prefix}{i}", genomeId, i * 100, i * 100 + 50, '+'))
            .ToArray();

    [Test]
    public void FindSyntenicBlocks_Schema_ValidatesCorrectly()
    {
        var g1 = Genome("g", "genome1");
        var g2 = Genome("h", "genome2");
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");
        Assert.DoesNotThrow(() => AnalysisTools.FindSyntenicBlocks(g1, g2, map));
    }

    [Test]
    public void FindSyntenicBlocks_Binding_InvokesSuccessfully()
    {
        var genome1 = Genome("g", "genome1");
        var genome2 = Genome("h", "genome2");
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");

        // gi -> hi (identical order) -> one forward block of 5.
        var blocks = AnalysisTools.FindSyntenicBlocks(genome1, genome2, map).Items;
        Assert.Multiple(() =>
        {
            Assert.That(blocks, Has.Length.EqualTo(1));
            Assert.That(blocks[0].GeneCount, Is.EqualTo(5));
            Assert.That(blocks[0].IsInverted, Is.False);
        });
    }
}
