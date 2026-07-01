using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class MergeOverlappingSvsTests
{
    private static StructuralVariantDto Sv(string id, int start, int end, string type, double q, int support) =>
        new(id, "chr1", start, end, type, end - start, q, support, null);

    // Two deletions overlapping by 0.5 of the min length (adjacent after sorting) -> merge to one.
    private static List<StructuralVariantDto> TwoDeletions() =>
        new()
        {
            Sv("del1", 1000, 2000, "Deletion", 30, 5),
            Sv("del2", 1500, 2500, "Deletion", 40, 3),
        };

    [Test]
    public void MergeOverlappingSvs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.MergeOverlappingSvs(TwoDeletions()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.MergeOverlappingSvs(null!));
    }

    [Test]
    public void MergeOverlappingSvs_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.MergeOverlappingSVs merges same-type SVs whose reciprocal overlap
        // fraction >= threshold. del1/del2 overlap 500/1000 = 0.5 -> merged span 1000..2500,
        // length 1500, support summed (8), quality max (40).
        var result = AnnotationTools.MergeOverlappingSvs(TwoDeletions());

        Assert.That(result.Merged, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Merged[0].Start, Is.EqualTo(1000));
            Assert.That(result.Merged[0].End, Is.EqualTo(2500));
            Assert.That(result.Merged[0].Length, Is.EqualTo(1500));
            Assert.That(result.Merged[0].SupportingReads, Is.EqualTo(8));
            Assert.That(result.Merged[0].Quality, Is.EqualTo(40.0).Within(1e-9));
        });
    }

    [Test]
    public void MergeOverlappingSvs_DifferentType_NotMerged()
    {
        // Same span but different type -> not merged (adjacent but type mismatch).
        var variants = new List<StructuralVariantDto>
        {
            Sv("del", 1000, 2000, "Deletion", 30, 5),
            Sv("dup", 1000, 2000, "Duplication", 25, 4),
        };
        var result = AnnotationTools.MergeOverlappingSvs(variants);
        Assert.That(result.Merged, Has.Count.EqualTo(2));
    }

    [Test]
    public void MergeOverlappingSvs_NonOverlapping_NotMerged()
    {
        var variants = new List<StructuralVariantDto>
        {
            Sv("d1", 1000, 2000, "Deletion", 30, 5),
            Sv("d2", 5000, 6000, "Deletion", 30, 5),
        };
        var result = AnnotationTools.MergeOverlappingSvs(variants);
        Assert.That(result.Merged, Has.Count.EqualTo(2));
    }
}
