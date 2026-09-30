// PAT-APPROX-002 — Linear-space (Hirschberg 1975) edit alignment.
// TestSpec: tests/TestSpecs/PAT-APPROX-002.md (reference cross-check table; no separate Evidence file for this unit)
// Sources: Hirschberg D.S. (1975) Commun. ACM 18(6):341–343; Myers & Miller (1988) CABIOS 4(1):11–17.
// Reference: edlib 1.3 edlib.align(q, t, mode='NW', task='path') editDistance (values 2026-09-30).
// The CIGAR of a co-optimal path is implementation-defined; each expected CIGAR below was checked to
// replay to both strings with cost = the edlib distance.
namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class ApproximateMatcher_EditAlignmentLinearSpace_Tests
{
    [TestCase("ACGT", "ACGT", 0, "4=")]
    [TestCase("ACGT", "AGT", 1, "1=1I2=")]
    [TestCase("kitten", "sitting", 3, "1X3=1X1=1D")]
    [TestCase("survey", "surgery", 2, "3=1X1=1D1=")]
    [TestCase("GATTACA", "GCATGCU", 4, "1=2X1=1X1=1X")]
    [TestCase("AAAAAAAAAACCCCC", "CCCCCAAAAAAAAAA", 10, "5X5=5X")]
    [TestCase("CGATACGT", "CGATACGTT", 1, "7=1D1=")]
    [TestCase("TGCATAT", "ATCCGAT", 4, "2X1=2X2=")]
    [TestCase("ACGTACGTACGTACGT", "TTACGTAACGTACGAT", 5, "2I1X6=1D6=1D1=")]
    [TestCase("A", "AA", 1, "1D1=")]
    [TestCase("A", "CC", 2, "1D1X")]
    [TestCase("", "ACG", 3, "3D")]
    [TestCase("ACG", "", 3, "3I")]
    [TestCase("", "", 0, "")]
    public void GetEditAlignmentLinearSpace_EdlibDistance_ValidCigar(string query, string target, int edlibDistance, string cigar)
    {
        var a = ApproximateMatcher.GetEditAlignmentLinearSpace(query, target);
        Assert.Multiple(() =>
        {
            Assert.That(a.Distance, Is.EqualTo(edlibDistance));
            Assert.That(a.Cigar, Is.EqualTo(cigar));
            Assert.That(ApproximateMatcher_EditAlignment_Tests.ReplayCost(a.Cigar, query, target), Is.EqualTo(edlibDistance));
            Assert.That(a.Distance, Is.EqualTo(ApproximateMatcher.GetEditAlignment(query, target).Distance));
        });
    }

    [Test]
    public void GetEditAlignmentLinearSpace_CoOptimalPathMayDifferFromFullTraceback()
    {
        // Same distance 5, different co-optimal paths (the documented reason for a separate method).
        var linear = ApproximateMatcher.GetEditAlignmentLinearSpace("ACGTACGTACGTACGT", "TTACGTAACGTACGAT");
        var full = ApproximateMatcher.GetEditAlignment("ACGTACGTACGTACGT", "TTACGTAACGTACGAT");
        Assert.That(linear.Distance, Is.EqualTo(full.Distance));
        Assert.That(linear.Operations, Is.EqualTo("IIX======D======D="));
        Assert.That(full.Operations, Is.EqualTo("IIX=====D=======D="));
    }

    [Test]
    public void GetEditAlignmentLinearSpace_AlignedStringsAndMetadata()
    {
        var a = ApproximateMatcher.GetEditAlignmentLinearSpace("survey", "surgery");
        Assert.Multiple(() =>
        {
            Assert.That(a.AlignedQuery, Is.EqualTo("surve-y"));
            Assert.That(a.AlignedTarget, Is.EqualTo("surgery"));
            Assert.That(a.SubstitutionPositions, Is.EqualTo(new[] { 3 }));
            Assert.That(a.HasIndels, Is.True);
            Assert.That(a.StandardCigar, Is.EqualTo("5M1D1M"));
        });
    }

    [Test]
    public void GetEditAlignmentLinearSpace_LongInputs_DistanceEqualsMyers()
    {
        // 3000 × 3300: the full traceback would hold ~10^7 cells; the linear-space form holds 4 columns.
        string q = new string('A', 2000) + new string('C', 1000);
        string t = new string('C', 1500) + new string('A', 1800);
        var a = ApproximateMatcher.GetEditAlignmentLinearSpace(q, t);
        Assert.That(a.Distance, Is.EqualTo(ApproximateMatcher.EditDistance(q, t)));
        Assert.That(ApproximateMatcher_EditAlignment_Tests.ReplayCost(a.Cigar, q, t), Is.EqualTo(a.Distance));
    }

    [Test]
    public void GetEditAlignmentLinearSpace_NullArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.GetEditAlignmentLinearSpace(null!, "A"));
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.GetEditAlignmentLinearSpace("A", null!));
    }
}
