using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ClassifyPreMiRnaMiPredTests
{
    // hsa-mir-21 hairpin part (MiRnaAnalyzer_MiPred_Tests MP8: 781 votes with the default seed and 999 shuffles).
    private const string Mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";

    [Test]
    public void ClassifyPreMiRnaMiPred_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyPreMiRnaMiPred(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyPreMiRnaMiPred(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRnaMiPred(Mir21, 1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRnaMiPred(Mir21, 1, 10001));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRnaMiPred(Mir21, 1, 99, -1));
    }

    [Test]
    public void ClassifyPreMiRnaMiPred_Binding_Mir21_Locked781Votes()
    {
        var r = AnnotationTools.ClassifyPreMiRnaMiPred(Mir21);
        Assert.Multiple(() =>
        {
            Assert.That(r.Found, Is.True);
            Assert.That(r.IsHairpinLike, Is.True);
            Assert.That(r.RealVotes, Is.EqualTo(781));          // library test MP8 (hsa-mir-21, default seed, 999 shuffles)
            Assert.That(r.TreeCount, Is.EqualTo(1000));
            Assert.That(r.RealProbability, Is.EqualTo(0.781));
            Assert.That(r.IsPreMiRna, Is.True);
            Assert.That(r.PValue, Is.LessThanOrEqualTo(0.01));
            Assert.That(r.Features, Has.Count.EqualTo(34));
            Assert.That(r.Features![32], Is.EqualTo(r.FreeEnergy));
            Assert.That(r.Features[33], Is.EqualTo(r.PValue));
            Assert.That(r.Length, Is.EqualTo(Mir21.Length));
            Assert.That(r.Seed, Is.EqualTo(MiRnaAnalyzer.MiPredDefaultSeed));
        });
    }

    [Test]
    public void ClassifyPreMiRnaMiPred_Binding_MatchesLibrary_SmallShuffleCount()
    {
        var r = AnnotationTools.ClassifyPreMiRnaMiPred(Mir21, 7, 49);
        var lib = MiRnaAnalyzer.ClassifyPreMiRnaMiPred(Mir21, new Random(7), 49)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(r.DotBracket, Is.EqualTo(lib.DotBracket));
            Assert.That(r.FreeEnergy, Is.EqualTo(lib.FreeEnergy));
            Assert.That(r.PValue, Is.EqualTo(lib.PValue));
            Assert.That(r.Features, Is.EqualTo(lib.Features));
            Assert.That(r.RealVotes, Is.EqualTo(lib.RealVotes));
            Assert.That(r.ShuffleCount, Is.EqualTo(49));
        });
    }

    [Test]
    public void ClassifyPreMiRnaMiPred_NotHairpinLike_NotClassified()
    {
        // 12 nt cannot satisfy the package pre-judgement (51..137 nt): no p-value, no features, no call.
        var r = AnnotationTools.ClassifyPreMiRnaMiPred("GGGGAAAACCCC");
        Assert.Multiple(() =>
        {
            Assert.That(r.Found, Is.True);
            Assert.That(r.IsHairpinLike, Is.False);
            Assert.That(r.PValue, Is.Null);
            Assert.That(r.Features, Is.Null);
            Assert.That(r.RealProbability, Is.Null);
            Assert.That(r.IsPreMiRna, Is.False);
        });
    }

    [Test]
    public void ClassifyPreMiRnaMiPred_NothingLeftAfterCleanUp_NotFound()
    {
        var r = AnnotationTools.ClassifyPreMiRnaMiPred("NNNN--");
        Assert.That(r.Found, Is.False);
        Assert.That(r.IsPreMiRna, Is.False);
    }
}
