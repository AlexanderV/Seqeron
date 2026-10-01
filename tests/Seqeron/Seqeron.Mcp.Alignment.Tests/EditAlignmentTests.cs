using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// edit_alignment — delegation to <c>ApproximateMatcher.GetEditAlignment</c> /
/// <c>GetEditAlignmentLinearSpace</c>. Expected CIGARs are the edlib-locked cases of
/// ApproximateMatcher_EditAlignment_Tests and ApproximateMatcher_EditAlignmentLinearSpace_Tests.
/// </summary>
[TestFixture]
public class EditAlignmentTests
{
    [Test]
    public void EditAlignment_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.EditAlignment("ACGT", "AGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.EditAlignment("", "ACG"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.EditAlignment("ACG", null!));
        Assert.Throws<ArgumentException>(() => AlignmentTools.EditAlignment(null!, "ACG", linearSpace: true));
    }

    [Test]
    public void EditAlignment_SurveySurgery_EqualsEdlib()
    {
        var r = AlignmentTools.EditAlignment("survey", "surgery");
        Assert.Multiple(() =>
        {
            Assert.That(r.Distance, Is.EqualTo(2));
            Assert.That(r.Cigar, Is.EqualTo("3=1X1=1D1="));
            Assert.That(r.StandardCigar, Is.EqualTo("5M1D1M"));
            Assert.That(r.Operations, Is.EqualTo("===X=D="));
            Assert.That(r.AlignedQuery, Is.EqualTo("surve-y"));
            Assert.That(r.AlignedTarget, Is.EqualTo("surgery"));
            Assert.That(r.SubstitutionPositions, Is.EqualTo(new[] { 3 }));
            Assert.That(r.HasIndels, Is.True);
        });
    }

    [TestCase("kitten", "sitting", 3, "1X3=1X1=1D")]
    [TestCase("ACGT", "AGT", 1, "1=1I2=")]
    [TestCase("AAAAAAAAAACCCCC", "CCCCCAAAAAAAAAA", 10, "5X5=5X")]
    public void EditAlignment_LinearSpace_EqualsLockedHirschberg(string query, string target, int distance, string cigar)
    {
        var r = AlignmentTools.EditAlignment(query, target, linearSpace: true);
        var lib = global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetEditAlignmentLinearSpace(query, target);
        Assert.Multiple(() =>
        {
            Assert.That(r.Distance, Is.EqualTo(distance));
            Assert.That(r.Cigar, Is.EqualTo(cigar));
            Assert.That(r.Operations, Is.EqualTo(lib.Operations));
            Assert.That(r.SubstitutionPositions, Is.EqualTo(lib.SubstitutionPositions));
        });
    }

    [Test]
    public void EditAlignment_FullMatrix_EqualsLibrary_DiagonalFirstTieBreak()
    {
        // GetEditAlignment_ACGT_ACGGT_DiagonalFirstTieBreak: library "2=1D2=" (edlib "3=1D1=").
        var r = AlignmentTools.EditAlignment("ACGT", "ACGGT");
        Assert.That(r.Cigar, Is.EqualTo("2=1D2="));
        var lib = global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetEditAlignment("ACGT", "ACGGT");
        Assert.That((r.Distance, r.Operations, r.AlignedQuery, r.AlignedTarget, r.HasIndels),
            Is.EqualTo((lib.Distance, lib.Operations, lib.AlignedQuery, lib.AlignedTarget, lib.HasIndels)));
    }
}
