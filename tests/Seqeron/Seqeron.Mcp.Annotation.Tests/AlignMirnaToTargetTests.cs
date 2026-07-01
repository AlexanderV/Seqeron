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
        });

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
