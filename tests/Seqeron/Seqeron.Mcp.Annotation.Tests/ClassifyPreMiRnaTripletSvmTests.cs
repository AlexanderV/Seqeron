using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ClassifyPreMiRnaTripletSvmTests
{
    private const string Mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";

    [Test]
    public void ClassifyPreMiRnaTripletSvm_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyPreMiRnaTripletSvm(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyPreMiRnaTripletSvm(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRnaTripletSvm(Mir21, null, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRnaTripletSvm(Mir21, null, -3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ClassifyPreMiRnaTripletSvm(Mir21, null, 22, -1));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "((((...))))", 3)); // length mismatch
    }

    [Test]
    public void ClassifyPreMiRnaTripletSvm_Binding_MatchesLibrary()
    {
        // Parity with the library (itself cross-checked against Xue's Perl scripts + libsvm in MiRnaAnalyzer_TripletSvm_Tests).
        var r = AnnotationTools.ClassifyPreMiRnaTripletSvm(Mir21, null, 16);
        var lib = MiRnaAnalyzer.ClassifyPreMiRnaTripletSvm(Mir21, null, 16);
        Assert.That(lib, Is.Not.Empty);
        Assert.That(r.Found, Is.True);
        Assert.That(r.Candidates, Has.Count.EqualTo(lib.Count));
        for (int i = 0; i < lib.Count; i++)
        {
            Assert.That(r.Candidates[i].DecisionValue, Is.EqualTo(lib[i].DecisionValue).Within(1e-15));
            Assert.That(r.Candidates[i].IsPreMiRna, Is.EqualTo(lib[i].IsPreMiRna));
            Assert.That(r.Candidates[i].Start, Is.EqualTo(lib[i].Start));
            Assert.That(r.Candidates[i].Features, Has.Count.EqualTo(MiRnaAnalyzer.TripletSvmFeatureCount));
        }
    }

    [Test]
    public void ClassifyPreMiRnaTripletSvm_SuppliedDotBracket_IsUsed()
    {
        var r = AnnotationTools.ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "((((....))))", 3);
        Assert.That(r.Found, Is.True);
        Assert.That(r.Candidates.Single().DotBracket, Is.EqualTo("((((....))))"));
        Assert.That((r.Candidates[0].Start, r.Candidates[0].End), Is.EqualTo((1, 12)));
    }

    [Test]
    public void ClassifyPreMiRnaTripletSvm_NoStemLoop_FoundFalse()
    {
        // Default minBasePairs = 22 > the 4 bp of this stem -> nothing retrieved (null/empty result -> flag).
        var r = AnnotationTools.ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "((((....))))");
        Assert.That(r.Found, Is.False);
        Assert.That(r.Candidates, Is.Empty);
    }
}
