using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindOrfsTests
{
    [Test]
    public void FindOrfs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindOrfs("ATGAAAAAAAAATAA", minLength: 1));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindOrfs(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindOrfs(null!));
    }

    [Test]
    public void FindOrfs_SimpleAtgTaaOrf_DetectsOrf()
    {
        // ATG + 3x AAA (Lys) + TAA. Mirrors GenomeAnnotator_ORF_Tests.CreateOrf("ATG", n, "TAA").
        var result = AnnotationTools.FindOrfs("ATGAAAAAAAAATAA", minLength: 1, searchBothStrands: false);

        Assert.That(result.Orfs, Has.Count.EqualTo(1));
        var orf = result.Orfs[0];
        Assert.Multiple(() =>
        {
            Assert.That(orf.Start, Is.EqualTo(0));
            Assert.That(orf.End, Is.EqualTo(15));
            Assert.That(orf.Frame, Is.EqualTo(1));
            Assert.That(orf.IsReverseComplement, Is.False);
            Assert.That(orf.Sequence, Is.EqualTo("ATGAAAAAAAAATAA"));
            Assert.That(orf.ProteinSequence, Is.EqualTo("MKKK*"));
        });
    }

    [Test]
    public void FindOrfs_BelowMinLength_Excluded()
    {
        // The MKKK ORF is 4 aa (incl. stop); minLength=10 should exclude it.
        var result = AnnotationTools.FindOrfs("ATGAAAAAAAAATAA", minLength: 10, searchBothStrands: false);
        Assert.That(result.Orfs, Is.Empty);
    }

    [Test]
    public void FindOrfs_NoStartCodon_RequireStart_ReturnsEmpty()
    {
        // Mirrors GenomeAnnotator_ORF_Tests.FindOrfs_NoStartCodon_RequireStart_ReturnsEmpty.
        var result = AnnotationTools.FindOrfs("GGGGGGGGGGTAAGGGGGGGGG", minLength: 1, searchBothStrands: false, requireStartCodon: true);
        Assert.That(result.Orfs, Is.Empty);
    }
}
