using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// needle_align / water_align / clustalw_guide_tree / tcoffee_blast_msa_library / tcoffee_tblastx_library (audit A3-17, B09 F141): the tools delegate to
/// <see cref="SequenceAligner"/>; the library tests lock the rows against the installed EMBOSS 6.6.0 / ClustalW 2.1 / t_coffee 13.41 on their committed oracles.
/// </summary>
[TestFixture]
public class NeedleWaterGuideTreeBlastToolTests
{
    private const string A = "ACGTACGTAGCTAGCTAGCATCGATCGATCGAT";
    private const string B = "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT";

    [Test]
    public void NeedleAlign_EqualsTheLibrary_AndTheDocumentedEmbossScore()
    {
        var tool = AlignmentTools.NeedleAlign(A, B);
        var lib = SequenceAligner.NeedleAlign(A, B);
        Assert.Multiple(() =>
        {
            Assert.That(tool.AlignedSequence1, Is.EqualTo(lib.AlignedSequence1));
            Assert.That(tool.AlignedSequence2, Is.EqualTo(lib.AlignedSequence2));
            Assert.That(tool.ExactScore, Is.EqualTo(lib.ExactScore));
            Assert.That(tool.ExactScore, Is.EqualTo(146.0), "needle a.fa b.fa (EDNAFULL, 10 / 0.5): Score: 146.0");
            Assert.That(tool.AlignmentType, Is.EqualTo("Global"));
            Assert.That(tool.Formatted, Is.Null);
        });
    }

    [Test]
    public void NeedleAlign_EndWeight_ChargesTheEndGaps()
    {
        // end gaps free (default): the overhang costs nothing; -endweight with -endopen / -endextend charges it
        var free = AlignmentTools.NeedleAlign("TTTTACGTACGT", "ACGTACGT");
        var weighted = AlignmentTools.NeedleAlign("TTTTACGTACGT", "ACGTACGT", endWeight: true, endOpen: 10, endExtend: 0.5);
        var lib = SequenceAligner.NeedleAlign("TTTTACGTACGT", "ACGTACGT", endWeight: true, endOpen: 10, endExtend: 0.5);
        Assert.Multiple(() =>
        {
            Assert.That(free.ExactScore, Is.GreaterThan(weighted.ExactScore));
            Assert.That(weighted.ExactScore, Is.EqualTo(lib.ExactScore));
            Assert.That(weighted.AlignedSequence1, Is.EqualTo(lib.AlignedSequence1));
        });
    }

    [Test]
    public void NeedleAlign_FormatsTheAlignment_AndValidates()
    {
        var formatted = AlignmentTools.NeedleAlign(A, B, format: "fasta", names: ["a", "b"]);
        Assert.That(formatted.Formatted, Does.StartWith(">a\n"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.NeedleAlign(null!, A));
        Assert.Throws<ArgumentException>(() => AlignmentTools.NeedleAlign(A, B, matrix: "NOSUCH"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.NeedleAlign(A, B, gapOpen: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.NeedleAlign(A, B, gapExtend: 11));
    }

    [Test]
    public void WaterAlign_EqualsTheLibrary_ReportsStartCoordinates()
    {
        const string s1 = "GGGGACGTACGTAGCTAGCTAAAA";
        const string s2 = "TTTTACGTACGTAGCTAGCTCCCC";
        var tool = AlignmentTools.WaterAlign(s1, s2);
        var lib = SequenceAligner.WaterAlign(s1, s2);
        Assert.Multiple(() =>
        {
            Assert.That(tool.AlignedSequence1, Is.EqualTo(lib.AlignedSequence1));
            Assert.That(tool.AlignedSequence2, Is.EqualTo(lib.AlignedSequence2));
            Assert.That(tool.StartPosition1, Is.EqualTo(lib.StartPosition1));
            Assert.That(tool.StartPosition2, Is.EqualTo(lib.StartPosition2));
            Assert.That(tool.ExactScore, Is.EqualTo(lib.ExactScore));
            Assert.That(tool.AlignmentType, Is.EqualTo("Local"));
            Assert.That(tool.AlignedSequence1, Does.Contain("ACGTACGTAGCTAGCT"));
            Assert.That(tool.StartPosition1, Is.EqualTo(4));
        });
    }

    [Test]
    public void WaterAlign_FormatsAndValidates()
    {
        var formatted = AlignmentTools.WaterAlign("ACGTACGTAGCT", "CGTACGTAG", format: "fasta");
        Assert.That(formatted.Formatted, Does.StartWith(">seq1"));
        Assert.That(AlignmentTools.WaterAlign(A, "").AlignedSequence1, Is.Empty, "water: a sequence shorter than two residues gives the empty alignment");
        Assert.Throws<ArgumentException>(() => AlignmentTools.WaterAlign(null!, A));
        Assert.Throws<ArgumentException>(() => AlignmentTools.WaterAlign(A, B, matrix: "NOSUCH"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.WaterAlign(A, B, gapOpen: 101));
    }

    private const string BlastpXml =
        "<BlastOutput><BlastOutput_iterations>" +
        "<Iteration><Iteration_iter-num>1</Iteration_iter-num><Iteration_hits>" +
        "<Hit><Hit_num>1</Hit_num><Hit_def>b</Hit_def><Hit_hsps><Hsp><Hsp_query-from>1</Hsp_query-from><Hsp_hit-from>1</Hsp_hit-from><Hsp_identity>4</Hsp_identity><Hsp_align-len>5</Hsp_align-len>" +
        "<Hsp_qseq>ACDEF</Hsp_qseq><Hsp_hseq>ACD-F</Hsp_hseq></Hsp></Hit_hsps></Hit></Iteration_hits></Iteration>" +
        "<Iteration><Iteration_iter-num>2</Iteration_iter-num><Iteration_hits></Iteration_hits></Iteration>" +
        "</BlastOutput_iterations></BlastOutput>";
}
