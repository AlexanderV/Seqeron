using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// The MCP-visible side of audit round 7 (B08 A7-1, A7-3): what a caller of the tools sees for the defects the confirmation
/// audit of round 6 found. Each test failed on the code the auditor read (cbcecfd3).
/// </summary>
[TestFixture]
public class RestrictionCrisprRound7WrapperTests
{
    private const string Guide = "GACGCATAAAGATGAGACGC";

    private static readonly string RepeatGenome = string.Join("TTTTT", Enumerable.Repeat(Guide + "TGG", 3));

    // A7-1: an enum arrives over MCP as any number; whether it is refused must not depend on the genome (a repeat guide returned 0).
    [TestCase(99, false)]
    [TestCase(99, true)]
    [TestCase(-1, false)]
    [TestCase(2, true)]
    public void CrisprCfdSpecificityScore_UndefinedAggregation_IsRefused_ForARepeatGuideToo(int aggregation, bool alternativePams)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.crispr_cfd_specificity_score(
            Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, (CfdAggregation)aggregation, alternativePams, max_occurrences: 2));
    }

    // A7-3: the recognition sequence and the name of a custom enzyme are bounded, like the stagger
    [Test]
    public void CustomEnzymeTools_RecognitionSequenceBeyondTheLimit_IsRefusedAtOnce()
    {
        string hugePattern = new('N', 200_000);

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_sites_custom_enzyme(
                new string('A', 2000), hugePattern, 1, 5, topology: MoleculeTopology.Circular));
            Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest_custom_enzyme(
                new string('A', 2000), hugePattern, 1, 5, topology: MoleculeTopology.Circular));
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_sites_custom_enzyme(
                new string('A', 200), new string('N', 65), 1, 5));
        });
    }

    [Test]
    public void CustomEnzymeTools_NameBeyondTheLimit_IsRefusedAtOnce()
    {
        string longName = new('X', 65);

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => MolToolsTools.find_sites_custom_enzyme(new string('A', 200), "A", 0, 0, name: longName));
            Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest_custom_enzyme(new string('A', 200), "A", 0, 0, name: longName));
        });
    }

    [Test]
    public void CustomEnzymeTools_RecognitionSequenceAndNameAtTheLimit_Work()
    {
        var sites = MolToolsTools.find_sites_custom_enzyme(
            new string('A', 100), new string('N', 64), 32, 32, name: new string('X', 64));

        Assert.That(sites.Sites, Has.Count.EqualTo(2 * (100 - 64 + 1)));
    }
}
