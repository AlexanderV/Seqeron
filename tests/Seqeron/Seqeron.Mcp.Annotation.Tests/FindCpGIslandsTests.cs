using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindCpGIslandsTests
{
    // 200 bp of perfect CpG dinucleotides: every window has GC=1.0 and CpG O/E=2.0.
    private static string CpGIsland() => string.Concat(Enumerable.Repeat("CG", 100));

    [Test]
    public void FindCpGIslands_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindCpGIslands(CpGIsland()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindCpGIslands(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindCpGIslands(null!));
    }

    [Test]
    public void FindCpGIslands_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.FindCpGIslands (Gardiner-Garden & Frommer): default minLength=200,
        // minGc=0.5, minCpGRatio=0.6. A 200 bp perfect CpG run is a single island 0..200,
        // GC content 1.0, CpG O/E 2.0.
        var result = AnnotationTools.FindCpGIslands(CpGIsland());

        Assert.That(result.Islands, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Islands[0].Start, Is.EqualTo(0));
            Assert.That(result.Islands[0].End, Is.EqualTo(200));
            Assert.That(result.Islands[0].GcContent, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(result.Islands[0].CpGRatio, Is.EqualTo(2.0).Within(1e-9));
        });
    }

    [Test]
    public void FindCpGIslands_ShorterThanMinLength_ReturnsEmpty()
    {
        // Below the default 200 bp minimum -> no islands.
        var shortSeq = string.Concat(Enumerable.Repeat("CG", 50)); // 100 bp
        Assert.That(AnnotationTools.FindCpGIslands(shortSeq).Islands, Is.Empty);
    }

    [Test]
    public void FindCpGIslands_LowGc_ReturnsEmpty()
    {
        // 300 bp AT-only: GC=0 < minGc -> no islands.
        var atRich = new string('A', 150) + new string('T', 150);
        Assert.That(AnnotationTools.FindCpGIslands(atRich).Islands, Is.Empty);
    }
}
