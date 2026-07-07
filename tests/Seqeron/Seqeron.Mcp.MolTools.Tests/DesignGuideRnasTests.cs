using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DesignGuideRnasTests
{
    // From Seqeron.Genomics.Tests CrisprDesigner_GuideRNA_Tests (S-003): a single SpCas9
    // guide with an NGG PAM ("AGG") in the region [20,45] at position 24, score 100.
    private const string Sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGG";

    [Test]
    public void DesignGuideRnas_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.design_guide_rnas(Sequence, 20, 45));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_guide_rnas("", 0, 10));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_guide_rnas(null!, 0, 10));
        // region_start out of range.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_guide_rnas(Sequence, -1, 10));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_guide_rnas(Sequence, Sequence.Length, Sequence.Length));
        // region_end out of range / before start.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_guide_rnas(Sequence, 0, Sequence.Length));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_guide_rnas(Sequence, 30, 10));
    }

    [Test]
    public void DesignGuideRnas_Binding_InvokesSuccessfully()
    {
        var guides = MolToolsTools.design_guide_rnas(Sequence, 20, 45, CrisprSystemType.SpCas9).Guides;

        Assert.Multiple(() =>
        {
            Assert.That(guides.Count, Is.EqualTo(1));
            var g = guides[0];
            Assert.That(g.Position, Is.EqualTo(24));
            Assert.That(g.Score, Is.EqualTo(100.0).Within(1e-9));
            Assert.That(g.IsForwardStrand, Is.True);
            Assert.That(g.GcContent, Is.EqualTo(50.0).Within(1e-9));
            Assert.That(g.SeedGcContent, Is.EqualTo(50.0).Within(1e-9));
            Assert.That(g.HasPolyT, Is.False);
            Assert.That(g.Sequence, Is.EqualTo("ACGTACGTACGTACGTACGT"));
        });
    }
}
