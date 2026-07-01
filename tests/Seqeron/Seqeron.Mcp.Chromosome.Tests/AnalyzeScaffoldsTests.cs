using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for the <c>analyze_scaffolds</c> MCP tool.
/// Expected values derive from GenomeAssemblyAnalyzer.AnalyzeScaffolds: contigs are gap-free runs
/// named "{id}_contig{n}" over inclusive [start, end]; an N-run of length >= minGapLength becomes a
/// gap classified by length (10 -> "Medium", per ClassifyGap).
/// </summary>
[TestFixture]
public class AnalyzeScaffoldsTests
{
    [Test]
    public void AnalyzeScaffolds_Schema_ValidatesCorrectly()
    {
        var valid = new List<NamedSequence> { new("scaf1", "AAAAANNNNNNNNNNTTTTT") };

        Assert.DoesNotThrow(() => ChromosomeTools.AnalyzeScaffolds(valid));
        Assert.DoesNotThrow(() => ChromosomeTools.AnalyzeScaffolds(new List<NamedSequence>()));

        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeScaffolds(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeScaffolds(valid, minGapLength: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeScaffolds(valid, minGapLength: -5));
    }

    [Test]
    public void AnalyzeScaffolds_Binding_InvokesSuccessfully()
    {
        // "AAAAA" (5) + "N"*10 + "TTTTT" (5): two contigs split by one 10-bp gap.
        var scaffolds = new List<NamedSequence>
        {
            new("scaf1", "AAAAA" + new string('N', 10) + "TTTTT"),
        };

        var result = ChromosomeTools.AnalyzeScaffolds(scaffolds, minGapLength: 10);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var s = result.Items[0];

        Assert.Multiple(() =>
        {
            Assert.That(s.ScaffoldId, Is.EqualTo("scaf1"));
            Assert.That(s.TotalLength, Is.EqualTo(20));
            Assert.That(s.ContigLength, Is.EqualTo(10));
            Assert.That(s.GapLength, Is.EqualTo(10));

            Assert.That(s.Contigs, Has.Count.EqualTo(2));
            Assert.That(s.Contigs[0].ContigId, Is.EqualTo("scaf1_contig1"));
            Assert.That(s.Contigs[0].Start, Is.EqualTo(0));
            Assert.That(s.Contigs[0].End, Is.EqualTo(4));
            Assert.That(s.Contigs[1].ContigId, Is.EqualTo("scaf1_contig2"));
            Assert.That(s.Contigs[1].Start, Is.EqualTo(15));
            Assert.That(s.Contigs[1].End, Is.EqualTo(19));

            Assert.That(s.Gaps, Has.Count.EqualTo(1));
            Assert.That(s.Gaps[0].Start, Is.EqualTo(5));
            Assert.That(s.Gaps[0].End, Is.EqualTo(14));
            Assert.That(s.Gaps[0].Length, Is.EqualTo(10));
            Assert.That(s.Gaps[0].GapType, Is.EqualTo("Medium"));
        });
    }

    [Test]
    public void AnalyzeScaffolds_ShortGapBelowThreshold_NotSplit()
    {
        // A 5-N run with minGapLength=10 is not recorded as a gap; the whole scaffold is one contig.
        var scaffolds = new List<NamedSequence> { new("s", "AAAAANNNNNTTTTT") };

        var result = ChromosomeTools.AnalyzeScaffolds(scaffolds, minGapLength: 10);
        var s = result.Items[0];

        Assert.Multiple(() =>
        {
            Assert.That(s.Gaps, Is.Empty);
            Assert.That(s.GapLength, Is.EqualTo(0));
        });
    }
}
