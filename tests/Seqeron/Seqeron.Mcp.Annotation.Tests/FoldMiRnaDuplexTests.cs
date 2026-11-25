using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FoldMiRnaDuplexTests
{
    private const string Let7a = "UGAGGUAGUAGGUUGUAUAGUU";

    [Test]
    public void FoldMiRnaDuplex_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FoldMiRnaDuplex("GUGUGAAUCG", "CGAUUCAUAC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FoldMiRnaDuplex("", "UUUU"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FoldMiRnaDuplex(null!, "UUUU"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FoldMiRnaDuplex("AAAA", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FoldMiRnaDuplex("AAAA", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.FoldMiRnaDuplex("AAAA", "UUUU", -1));
    }

    [Test]
    public void FoldMiRnaDuplex_Binding_PerfectLet7aDuplex()
    {
        // Mirrors Seqeron.Genomics.Tests FoldMiRnaDuplex_PerfectComplement (ViennaRNA 2.8 duplexfold: same structure, -35.6;
        // NNDB Turner 2004: -35.90).
        string target = MiRnaAnalyzer.GetReverseComplement(Let7a);
        var r = AnnotationTools.FoldMiRnaDuplex(Let7a, target);

        Assert.That(r.Found, Is.True);
        Assert.Multiple(() =>
        {
            Assert.That(r.Duplex!.Structure, Is.EqualTo(new string('(', 22) + "&" + new string(')', 22)));
            Assert.That(r.Duplex.BasePairs, Has.Count.EqualTo(22));
            Assert.That(r.Duplex.BasePairs[0], Is.EqualTo(new DuplexBasePairDto(0, 21)));
            Assert.That(r.Duplex.FreeEnergy, Is.EqualTo(-35.90).Within(1e-9));
        });
    }

    [Test]
    public void FoldMiRnaDuplex_Binding_OneNtBulge_MatchesSourcedStructure()
    {
        // Library test Fold_OneNtBulge_A12: NNDB -32.10 (ViennaRNA 2.8 duplexfold -31.8).
        var r = AnnotationTools.FoldMiRnaDuplex(Let7a, "AACUAUACAACACUACUACCUCA");
        Assert.That(r.Duplex!.Structure, Is.EqualTo("((((((((((((((((((((((&))))))))))).)))))))))))"));
        Assert.That(r.Duplex.FreeEnergy, Is.EqualTo(-32.10).Within(1e-9));
    }

    [Test]
    public void FoldMiRnaDuplex_NoPairPossible_NotFound()
    {
        var r = AnnotationTools.FoldMiRnaDuplex("AAAA", "AAAA");
        Assert.That(r.Found, Is.False);
        Assert.That(r.Duplex, Is.Null);
    }

    [Test]
    public void FoldMiRnaDuplex_DelegatesToLibrary_IncludingMaxLoopSize()
    {
        const string target = "AACUAUACAACACUACUACCUCA";
        foreach (int loop in new[] { 0, 1, 30 })
        {
            var direct = MiRnaAnalyzer.FoldMiRnaDuplex(Let7a, target, loop)!.Value;
            var d = AnnotationTools.FoldMiRnaDuplex(Let7a, target, loop).Duplex!;
            Assert.Multiple(() =>
            {
                Assert.That(d.Structure, Is.EqualTo(direct.Structure), $"maxLoopSize {loop}");
                Assert.That(d.FreeEnergy, Is.EqualTo(direct.FreeEnergy));
                Assert.That(d.MiRnaStart, Is.EqualTo(direct.MiRnaStart));
                Assert.That(d.MiRnaEnd, Is.EqualTo(direct.MiRnaEnd));
                Assert.That(d.TargetStart, Is.EqualTo(direct.TargetStart));
                Assert.That(d.TargetEnd, Is.EqualTo(direct.TargetEnd));
                Assert.That(d.BasePairs.Select(p => (p.MiRnaIndex, p.TargetIndex)), Is.EqualTo(direct.BasePairs));
            });
        }
    }
}
