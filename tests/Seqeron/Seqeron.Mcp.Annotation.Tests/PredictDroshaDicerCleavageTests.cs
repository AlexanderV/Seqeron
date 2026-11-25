using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictDroshaDicerCleavageTests
{
    // Mirrors Seqeron.Genomics.Tests MiRnaAnalyzer_StructureCleavageKernels_Tests: 5-nt flank, 35-bp stem,
    // 8-nt loop, 35-bp stem, 25-nt flank (108 nt); pair k = (5+k, 82-k).
    private static readonly string Pri =
        new string('A', 5) + new string('G', 35) + new string('A', 8) + new string('C', 35) + new string('A', 25);
    private static readonly string Fold =
        new string('.', 5) + new string('(', 35) + new string('.', 8) + new string(')', 35) + new string('.', 25);

    [Test]
    public void PredictDroshaDicerCleavage_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictDroshaDicerCleavage(Pri, 5, Fold));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictDroshaDicerCleavage("", 0));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictDroshaDicerCleavage(null!, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.PredictDroshaDicerCleavage(Pri, -1, Fold));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.PredictDroshaDicerCleavage(Pri, Pri.Length, Fold));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictDroshaDicerCleavage(Pri, 5, Fold, dicerRule: "Bogus"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictDroshaDicerCleavage(Pri, 5, Fold, dicerRule: "1"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictDroshaDicerCleavage(Pri, 5, Fold, droshaRule: "Bogus"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictDroshaDicerCleavage(Pri, 5, "(((", "FivePrimeCounting", "BasalJunction"),
            "unbalanced / wrong-length dot-bracket is rejected by the library");
    }

    [Test]
    public void PredictDroshaDicerCleavage_Binding_ExactFoldGeometry()
    {
        // Library test: junction 5 -> 12th pair (16, 71): Drosha 5' cut 16, 3' cut 71+2 = 73; 5p [16..37]; 3p [52..73].
        var r = AnnotationTools.PredictDroshaDicerCleavage(Pri, 5, Fold);
        Assert.That(r.Predicted, Is.True);
        var c = r.Cleavage!;
        Assert.Multiple(() =>
        {
            Assert.That(c.DroshaCut5Prime, Is.EqualTo(16));
            Assert.That(c.MatureStart, Is.EqualTo(16));
            Assert.That(c.MatureEnd, Is.EqualTo(37));
            Assert.That(c.StarStart, Is.EqualTo(52));
            Assert.That(c.StarEnd, Is.EqualTo(73));
            Assert.That(c.DroshaCut3Prime, Is.EqualTo(73));
            Assert.That(c.ThreePrimeOverhang, Is.EqualTo(2));
            Assert.That(c.DicerRule, Is.EqualTo("FivePrimeCounting"));
            Assert.That(c.DroshaRule, Is.EqualTo("BasalJunction"));
        });
    }

    [Test]
    public void PredictDroshaDicerCleavage_ShortStem_NotPredicted()
    {
        // A 10-bp stem has no 12th pair (library returns null).
        string seq = "AAAAA" + new string('G', 10) + "AAAA" + new string('C', 10) + "AAAAA";
        string db = "....." + new string('(', 10) + "...." + new string(')', 10) + ".....";
        var r = AnnotationTools.PredictDroshaDicerCleavage(seq, 5, db);
        Assert.That(r.Predicted, Is.False);
        Assert.That(r.Cleavage, Is.Null);
    }

    [Test]
    public void PredictDroshaDicerCleavage_DelegatesToLibrary_AllRuleCombinations()
    {
        foreach (var dicer in Enum.GetValues<MiRnaAnalyzer.DicerCountingRule>())
        foreach (var drosha in Enum.GetValues<MiRnaAnalyzer.DroshaCountingRule>())
        {
            var direct = MiRnaAnalyzer.PredictDroshaDicerCleavage(Pri, 5, Fold, dicer, drosha);
            var wrapped = AnnotationTools.PredictDroshaDicerCleavage(Pri, 5, Fold, dicer.ToString().ToLowerInvariant(), drosha.ToString());
            Assert.That(wrapped.Predicted, Is.EqualTo(direct.HasValue), $"{dicer}/{drosha}");
            if (!direct.HasValue) continue;
            var d = direct.Value;
            var c = wrapped.Cleavage!;
            Assert.Multiple(() =>
            {
                Assert.That((c.DroshaCut5Prime, c.DroshaCut3Prime, c.MatureStart, c.MatureEnd, c.StarStart, c.StarEnd),
                    Is.EqualTo((d.DroshaCut5Prime, d.DroshaCut3Prime, d.MatureStart, d.MatureEnd, d.StarStart, d.StarEnd)), $"{dicer}/{drosha}");
                Assert.That(c.MatureSequence, Is.EqualTo(d.MatureSequence));
                Assert.That(c.StarSequence, Is.EqualTo(d.StarSequence));
                Assert.That(c.ApicalJunction, Is.EqualTo(d.ApicalJunction));
                Assert.That(c.HasCnncMotif, Is.EqualTo(d.HasCnncMotif));
                Assert.That(c.HasBasalUgMotif, Is.EqualTo(d.HasBasalUgMotif));
                Assert.That(c.DicerRule, Is.EqualTo(dicer.ToString()));
                Assert.That(c.DroshaRule, Is.EqualTo(drosha.ToString()));
            });
        }
    }

    [Test]
    public void PredictDroshaDicerCleavage_NullDotBracket_UsesMfeFold_LikeLibrary()
    {
        var direct = MiRnaAnalyzer.PredictDroshaDicerCleavage(Pri, 5);
        var wrapped = AnnotationTools.PredictDroshaDicerCleavage(Pri, 5);
        Assert.That(wrapped.Predicted, Is.EqualTo(direct.HasValue));
        if (direct.HasValue)
        {
            Assert.That(wrapped.Cleavage!.Structure, Is.EqualTo(direct.Value.Structure));
            Assert.That(wrapped.Cleavage.MatureSequence, Is.EqualTo(direct.Value.MatureSequence));
        }
    }
}
