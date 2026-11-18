using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class AlignmentStatisticsTests
{
    [Test]
    public void AlignmentStatistics_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.AlignmentStatistics("ACGT", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentStatistics("", "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentStatistics("ACGT", null!));
        // Unequal aligned lengths are invalid.
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentStatistics("ACGT", "ACG"));
    }

    [Test]
    public void AlignmentStatistics_Binding_InvokesSuccessfully()
    {
        // EMBOSS needle convention: denominator is the full alignment length (incl. gap columns).
        // "ACGT-A" vs "ACGTGA": 5 identical columns, 0 mismatch, 1 gap column, length 6.
        var r = AlignmentTools.AlignmentStatistics("ACGT-A", "ACGTGA");
        Assert.Multiple(() =>
        {
            Assert.That(r.Matches, Is.EqualTo(5));
            Assert.That(r.Mismatches, Is.EqualTo(0));
            Assert.That(r.Gaps, Is.EqualTo(1));
            Assert.That(r.AlignmentLength, Is.EqualTo(6));
            Assert.That(r.Identity, Is.EqualTo(500.0 / 6).Within(1e-9));
            Assert.That(r.Similarity, Is.EqualTo(500.0 / 6).Within(1e-9));
            Assert.That(r.GapPercent, Is.EqualTo(100.0 / 6).Within(1e-9));
        });

        // Perfect identity.
        var id = AlignmentTools.AlignmentStatistics("ACGT", "ACGT");
        Assert.Multiple(() =>
        {
            Assert.That(id.Matches, Is.EqualTo(4));
            Assert.That(id.Gaps, Is.EqualTo(0));
            Assert.That(id.Identity, Is.EqualTo(100.0));
            Assert.That(id.GapPercent, Is.EqualTo(0.0));
        });
    }

    // ALIGN-STATS-001 (B09 F11-F14): the published EMBOSS needle example (HBA_HUMAN vs HBB_HUMAN, EBLOSUM62):
    // Length 149, Identity 65, Similarity 90, Gaps 9 — reproduced with the real BLOSUM62 table through the tool.
    private const string Hba =
        "MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
    private const string Hbb =
        "MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";

    [Test]
    public void AlignmentStatistics_Blosum62_ReproducesPublishedEmbossExample()
    {
        var r = AlignmentTools.AlignmentStatistics(Hba, Hbb, "blosum62");
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignmentLength, Is.EqualTo(149));
            Assert.That(r.Matches, Is.EqualTo(65));
            Assert.That(r.Positives, Is.EqualTo(90));
            Assert.That(r.Gaps, Is.EqualTo(9));
            Assert.That(r.GapOpens, Is.EqualTo(4));
            Assert.That(r.Identity, Is.EqualTo(43.6).Within(0.05));
            Assert.That(r.Similarity, Is.EqualTo(60.4).Within(0.05));
            Assert.That(r.GapPercent, Is.EqualTo(6.0).Within(0.05));
        });
    }

    [Test]
    public void AlignmentStatistics_ScalarModel_PositiveMismatchIsSimilar_AndGapOpensCounted()
    {
        // Default scalar model (+1/-1): similarity == identity. mismatch = +1: the C/G column is similar.
        var d = AlignmentTools.AlignmentStatistics("AC--GT", "AGTTGT");
        Assert.That(d.GapOpens, Is.EqualTo(1));
        Assert.That(d.Similarity, Is.EqualTo(d.Identity));
        var p = AlignmentTools.AlignmentStatistics("AC--GT", "AGTTGT", mismatch: 1);
        Assert.That(p.Positives, Is.EqualTo(p.Matches + p.Mismatches));
        Assert.That(p.Similarity, Is.EqualTo(100.0 * 4 / 6).Within(1e-9));
    }

    [Test]
    public void AlignmentStatistics_LowerCase_IsCaseInsensitive_AndUnknownMatrixRejected()
    {
        Assert.That(AlignmentTools.AlignmentStatistics("acgt", "ACGT").Matches, Is.EqualTo(4));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentStatistics("ACGT", "ACGT", "BLOSUM99"));
    }
}
