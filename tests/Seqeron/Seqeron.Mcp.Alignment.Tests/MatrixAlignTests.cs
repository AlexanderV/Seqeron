using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// MCP pairwise tools with the additive <c>matrix</c> parameter (B09 F41, audit A1-2): BLOSUM62 / PAM250 / EDNAFULL residue-pair
/// scoring and fractional gap scores delegate to <c>SequenceAligner.*(…, MatrixScoring)</c>. Expected values: EMBOSS needle/water
/// 6.6.0 (-gapopen 10 -gapextend 0.5, EBLOSUM62; needle with -endweight Y) = Biopython 1.88 PairwiseAligner on HBA_HUMAN/HBB_HUMAN.
/// </summary>
[TestFixture]
public class MatrixAlignTests
{
    private const string Hba = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
    private const string Hbb = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";

    [Test]
    public void GlobalAlignAffine_Blosum62_EmbossDefaults_Score292_5()
    {
        var r = AlignmentTools.GlobalAlignAffine(Hba, Hbb, gapOpen: -9.5, gapExtend: -0.5, matrix: "BLOSUM62");
        Assert.Multiple(() =>
        {
            Assert.That(r.ExactScore, Is.EqualTo(292.5));
            Assert.That(r.Score, Is.EqualTo(293));
            Assert.That(r.AlignmentType, Is.EqualTo("Global"));
        });
    }

    [Test]
    public void LocalAlignAffine_Blosum62_EmbossWaterDefaults_Score293_5()
    {
        var r = AlignmentTools.LocalAlignAffine(Hba, Hbb, gapOpen: -9.5, gapExtend: -0.5, matrix: "blosum62");
        Assert.That(r.ExactScore, Is.EqualTo(293.5));
    }

    [Test]
    public void LinearTools_Matrix_ScoreFractionalGapPerPosition()
    {
        // linear 0.5 per gap position: CW / CAW = 9 + 11 - 0.5
        var global = AlignmentTools.GlobalAlign("CW", "CAW", gapExtend: -0.5, matrix: "BLOSUM62");
        var local = AlignmentTools.LocalAlign("CW", "CAW", gapExtend: -0.5, matrix: "BLOSUM62");
        var semi = AlignmentTools.SemiGlobalAlign("CW", "CAW", gapExtend: -0.5, matrix: "BLOSUM62");
        Assert.Multiple(() =>
        {
            Assert.That(global.ExactScore, Is.EqualTo(19.5));
            Assert.That(local.ExactScore, Is.EqualTo(19.5));
            Assert.That(semi.ExactScore, Is.EqualTo(19.5));
            Assert.That(semi.AlignmentType, Is.EqualTo("SemiGlobal"));
        });
    }

    [Test]
    public void SemiGlobal_Matrix_AcceptsProtein_WhileScalarModeStillRequiresDna()
    {
        Assert.DoesNotThrow(() => AlignmentTools.SemiGlobalAlignAffine("MVHL", Hbb, gapOpen: -9.5, gapExtend: -0.5, matrix: "BLOSUM62"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SemiGlobalAlignAffine("MVHL", Hbb));
    }

    [Test]
    public void EdnaFull_ScoresIupacSymbols()
    {
        // EDNAFULL: A/A 5, N/A -2 (EMBOSS data file)
        var r = AlignmentTools.GlobalAlign("ACGN", "ACGA", matrix: "EDNAFULL", gapExtend: -2);
        Assert.That(r.ExactScore, Is.EqualTo(5 + 5 + 5 - 2));
    }

    [Test]
    public void ScalarMode_UnchangedByAdditiveParameters_ExactScoreEqualsScore()
    {
        var r = AlignmentTools.GlobalAlign("TTTTTTTTTTACGT", "ACGT", 2, -3, -5, -2);
        Assert.Multiple(() =>
        {
            Assert.That(r.Score, Is.EqualTo(-12));
            Assert.That(r.ExactScore, Is.EqualTo(-12.0));
        });
    }

    [Test]
    public void FractionalGap_WithoutMatrix_IsRejectedWithGuidance()
    {
        var ex = Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlignAffine("ACGT", "ACGA", gapExtend: -0.5));
        Assert.That(ex!.Message, Does.Contain("matrix"));
    }

    [Test]
    public void UnknownMatrix_Rejected_UnknownResidue_Rejected()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("ACGT", "ACGT", matrix: "BLOSUM45"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("AC1T", "ACGT", matrix: "BLOSUM62"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlignAffine("", "ACGT", matrix: "BLOSUM62"));
    }

    /// <summary>BLAST+ report values of alignment_statistics (B09 F42): integer percentages (NCBI GetPercentMatch), positives, pident.</summary>
    [Test]
    public void AlignmentStatistics_BlastReportFields()
    {
        // 199 identities of 200 columns -> 99.5 % would round to 100 but BLAST shows 99
        string a = new string('A', 199) + "C";
        string b = new string('A', 200);
        var r = AlignmentTools.AlignmentStatistics(a, b);
        var protein = AlignmentTools.AlignmentStatistics("AXL", "AXI", "BLOSUM62");
        Assert.Multiple(() =>
        {
            Assert.That(r.BlastIdentityPercent, Is.EqualTo(99));
            Assert.That(r.BlastPident, Is.EqualTo("99.500"));
            Assert.That(r.BlastReport, Is.EqualTo(" Identities = 199/200 (99%), Gaps = 0/200 (0%)"));
            Assert.That(protein.BlastPositives, Is.EqualTo(3));
            Assert.That(protein.BlastReport, Is.EqualTo(" Identities = 2/3 (67%), Positives = 3/3 (100%), Gaps = 0/3 (0%)"));
        });
    }
}
