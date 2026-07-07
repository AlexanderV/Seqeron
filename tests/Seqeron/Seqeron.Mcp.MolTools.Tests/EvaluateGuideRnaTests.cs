using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class EvaluateGuideRnaTests
{
    // 20-nt guide, 50% GC, no TTTT. SpCas9 seed = last 10 = "GCATGCATGC" -> 6/10 GC = 60%.
    private const string Guide = "ATGCATGCATGCATGCATGC";

    [Test]
    public void EvaluateGuideRna_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.evaluate_guide_rna(Guide));
        Assert.Throws<ArgumentException>(() => MolToolsTools.evaluate_guide_rna(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.evaluate_guide_rna(null!));
    }

    [Test]
    public void EvaluateGuideRna_Binding_InvokesSuccessfully()
    {
        var c = MolToolsTools.evaluate_guide_rna(Guide, CrisprSystemType.SpCas9);

        Assert.Multiple(() =>
        {
            Assert.That(c.Sequence, Is.EqualTo(Guide));
            Assert.That(c.Position, Is.EqualTo(-1));           // ad-hoc evaluation
            Assert.That(c.IsForwardStrand, Is.True);
            Assert.That(c.GcContent, Is.EqualTo(50.0).Within(1e-9));
            Assert.That(c.SeedGcContent, Is.EqualTo(60.0).Within(1e-9));
            Assert.That(c.HasPolyT, Is.False);
            Assert.That(c.Score, Is.InRange(0.0, 100.0));
        });

        // A guide containing TTTT is flagged as having a polyT terminator.
        var polyT = MolToolsTools.evaluate_guide_rna("ATGCATGCATGCATTTTGCA");
        Assert.That(polyT.HasPolyT, Is.True);
    }
}
