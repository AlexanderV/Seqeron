using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AlignMirnaToTargetTests
{
    [Test]
    public void AlignMirnaToTarget_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.AlignMiRnaToTarget("AAAA", "UUUU"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AlignMiRnaToTarget("", "UUUU"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AlignMiRnaToTarget(null!, "UUUU"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AlignMiRnaToTarget("AAAA", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AlignMiRnaToTarget("AAAA", null!));
    }

    [Test]
    public void AlignMirnaToTarget_Binding_InvokesSuccessfully()
    {
        // Expected values mirror Seqeron.Genomics.Tests M-011/M-012/M-013
        // (MiRnaAnalyzer_TargetPrediction_Tests). A pairs with U antiparallel → all Watson-Crick.
        var perfect = AnnotationTools.AlignMiRnaToTarget("AAAA", "UUUU");
        Assert.Multiple(() =>
        {
            Assert.That(perfect.Matches, Is.EqualTo(4));
            Assert.That(perfect.Mismatches, Is.EqualTo(0));
            Assert.That(perfect.GuWobbles, Is.EqualTo(0));
            Assert.That(perfect.Gaps, Is.EqualTo(0));
            Assert.That(perfect.MiRnaSequence, Is.EqualTo("AAAA"));
            Assert.That(perfect.TargetSequence, Is.EqualTo("UUUU"));
            Assert.That(perfect.AlignmentString, Is.EqualTo("||||"));
            // NNDB Turner 2004: init 4.09 + 3 × 5'AA/3'UU (−0.93) + 2 × AU end 0.45 = +2.20
            // (ViennaRNA 2.8 eval_structure AAAA&UUUU = +2.40, 0.1-rounded parameters).
            Assert.That(perfect.FreeEnergy, Is.EqualTo(2.20).Within(1e-9));
        });

        // let-7a vs its perfect complement: Turner 2004 −35.90 (ViennaRNA 2.8 eval_structure −35.6).
        var let7a = AnnotationTools.AlignMiRnaToTarget("UGAGGUAGUAGGUUGUAUAGUU", "AACUAUACAACCUACUACCUCA");
        Assert.That(let7a.FreeEnergy, Is.EqualTo(-35.90).Within(1e-9));

        // G:U wobble pairs (Crick 1966) — 4 wobbles, 0 Watson-Crick matches.
        var wobble = AnnotationTools.AlignMiRnaToTarget("GGGG", "UUUU");
        Assert.Multiple(() =>
        {
            Assert.That(wobble.GuWobbles, Is.EqualTo(4));
            Assert.That(wobble.Matches, Is.EqualTo(0));
            Assert.That(wobble.AlignmentString, Is.EqualTo("::::"));
        });

        // A does not pair with A → all mismatches.
        var mismatch = AnnotationTools.AlignMiRnaToTarget("AAAA", "AAAA");
        Assert.That(mismatch.Mismatches, Is.EqualTo(4));

        // DNA T is normalised to RNA U in the reported sequences.
        var dna = AnnotationTools.AlignMiRnaToTarget("AAAA", "TTTT");
        Assert.Multiple(() =>
        {
            Assert.That(dna.TargetSequence, Is.EqualTo("UUUU"));
            Assert.That(dna.Matches, Is.EqualTo(4));
        });
    }
}
