using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DigestSummaryTests
{
    [Test]
    public void DigestSummary_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.digest_summary("AAAGAATTCAAA", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.digest_summary("", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.digest_summary("AAAGAATTCAAA", Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.digest_summary("AAAGAATTCAAA", null!));
    }

    [Test]
    public void DigestSummary_Binding_InvokesSuccessfully()
    {
        // EcoRI on "AAAGAATTCAAA" -> fragments of length 4 and 8.
        var summary = MolToolsTools.digest_summary("AAAGAATTCAAA", new[] { "EcoRI" });

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalFragments, Is.EqualTo(2));
            // Sizes are sorted descending.
            Assert.That(summary.FragmentSizes, Is.EqualTo(new[] { 8, 4 }));
            Assert.That(summary.LargestFragment, Is.EqualTo(8));
            Assert.That(summary.SmallestFragment, Is.EqualTo(4));
            Assert.That(summary.AverageFragmentSize, Is.EqualTo(6.0).Within(1e-9));
            Assert.That(summary.EnzymesUsed, Is.EqualTo(new[] { "EcoRI" }));
            Assert.That(summary.Topology, Is.EqualTo(MoleculeTopology.Linear),
                "the tool defaults to Linear");
        });
    }

    /// <summary>
    /// The tool exposes <c>topology</c>, so a circular (plasmid) summary is reachable over MCP.
    /// On the 19 bp circle "TTCGGATCCAAAAAAAGAA" with EcoRI + BamHI, Biopython 1.88
    /// <c>RestrictionBatch([EcoRI, BamHI]).search(Seq(...), linear=False)
    /// == {EcoRI: [18], BamHI: [5]}</c> — 0-based cuts {4, 17} -> 2 circular fragments of 13 and
    /// 6 bp, whereas the linear digest sees only BamHI (cut 4) -> 15 and 4 bp.
    /// </summary>
    [Test]
    public void DigestSummary_Topology_CircularSeesOriginSpanningSite()
    {
        var linear = MolToolsTools.digest_summary("TTCGGATCCAAAAAAAGAA", new[] { "EcoRI", "BamHI" });
        var circular = MolToolsTools.digest_summary(
            "TTCGGATCCAAAAAAAGAA", new[] { "EcoRI", "BamHI" }, MoleculeTopology.Circular);

        Assert.Multiple(() =>
        {
            Assert.That(linear.FragmentSizes, Is.EqualTo(new[] { 15, 4 }));
            Assert.That(linear.Topology, Is.EqualTo(MoleculeTopology.Linear));

            Assert.That(circular.TotalFragments, Is.EqualTo(2));
            Assert.That(circular.FragmentSizes, Is.EqualTo(new[] { 13, 6 }));
            Assert.That(circular.FragmentSizes.Sum(), Is.EqualTo(19));
            Assert.That(circular.Topology, Is.EqualTo(MoleculeTopology.Circular));
        });
    }
}
