using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CompareSeedRegionsTests
{
    private static MiRnaDto Make(string name, string seq) =>
        AnnotationTools.CreateMiRna(name, seq).MiRna;

    // Seeds: A = GAGGUAG, C = GAGGUCG (differ at one position).
    private static MiRnaDto MirA() => Make("a", "UGAGGUAGUAGGUUGUAUAGUU");
    private static MiRnaDto MirC() => Make("c", "UGAGGUCGUAGGUUGUAUAGUU");

    [Test]
    public void CompareSeedRegions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CompareSeedRegions(MirA(), MirC()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.CompareSeedRegions(null!, MirC()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.CompareSeedRegions(MirA(), null!));
        var emptySeq = new MiRnaDto("x", "", "", 1, 7);
        Assert.Throws<ArgumentException>(() => AnnotationTools.CompareSeedRegions(emptySeq, MirC()));
    }

    [Test]
    public void CompareSeedRegions_SameFamily()
    {
        // Identical seeds (GAGGUAG) -> 7 matches, 0 mismatches, same family.
        var result = AnnotationTools.CompareSeedRegions(MirA(), MirA());
        Assert.Multiple(() =>
        {
            Assert.That(result.Matches, Is.EqualTo(7));
            Assert.That(result.Mismatches, Is.EqualTo(0));
            Assert.That(result.IsSameFamily, Is.True);
        });
    }

    [Test]
    public void CompareSeedRegions_OneMismatch_DifferentFamily()
    {
        // GAGGUAG vs GAGGUCG differ at one position -> 6 matches, 1 mismatch, not same family.
        var result = AnnotationTools.CompareSeedRegions(MirA(), MirC());
        Assert.Multiple(() =>
        {
            Assert.That(result.Matches, Is.EqualTo(6));
            Assert.That(result.Mismatches, Is.EqualTo(1));
            Assert.That(result.IsSameFamily, Is.False);
        });
    }
}
