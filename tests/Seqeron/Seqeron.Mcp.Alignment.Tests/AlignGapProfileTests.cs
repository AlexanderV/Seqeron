using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>align_gap_profile (audit A3-10, B09 F137): Biopython 1.88 per-direction / per-end gap scores. The tool delegates to SequenceAligner.*Align*(…, GapScoreProfile); the 1000-row oracle comparison lives in the library tests.</summary>
[TestFixture]
public class AlignGapProfileTests
{
    private static Dictionary<string, double> LeftInsertionFree() => new()
    {
        ["open_gap_score"] = -5,
        ["extend_gap_score"] = -5,
        ["open_left_insertion_score"] = 0,
        ["extend_left_insertion_score"] = 0,
    };

    [Test]
    public void AlignGapProfile_BiopythonWorkedValues()
    {
        // PairwiseAligner(mode="global", match 1, mismatch -1, gap -5, left_insertion_score = 0): score("AC", "GGAC") = 2.0, score("GGAC", "AC") = -8.0 (run on Biopython 1.88)
        Assert.That(AlignmentTools.AlignGapProfile("AC", "GGAC", gapScores: LeftInsertionFree()).Score, Is.EqualTo(2));
        Assert.That(AlignmentTools.AlignGapProfile("GGAC", "AC", gapScores: LeftInsertionFree()).Score, Is.EqualTo(-8));
        // the deprecated names resolve to the same slot (target = sequence 1 ... a gap in sequence 1 = insertion)
        var old = new Dictionary<string, double> { ["gap_score"] = -5, ["target_left_gap_score"] = 0 };
        Assert.That(AlignmentTools.AlignGapProfile("AC", "GGAC", gapScores: old).Score, Is.EqualTo(2));
    }

    [Test]
    public void AlignGapProfile_ModesEndGapsAndMatrix()
    {
        // free overhangs via endGaps = Fitting reproduce semi_global_align; the result is tagged SemiGlobal and the fitted interval is excluded from the coordinates
        var fit = AlignmentTools.AlignGapProfile("GATTACA", "CCGATTTACAGG", "semi_global", affine: false, match: 2, mismatch: -3, gapOpen: -2, gapExtend: -2, endGaps: "Fitting");
        var legacy = AlignmentTools.SemiGlobalAlign("GATTACA", "CCGATTTACAGG", 2, -3, -2, -2);
        Assert.That((fit.Score, fit.AlignmentType, fit.StartPosition2, fit.EndPosition2), Is.EqualTo((legacy.Score, "SemiGlobal", legacy.StartPosition2, legacy.EndPosition2)));
        // default profile = the ordinary affine global alignment
        var plain = AlignmentTools.GlobalAlignAffine("ACGTTGCA", "ACGGTGCA", 1, -1, -2, -1);
        var viaProfile = AlignmentTools.AlignGapProfile("ACGTTGCA", "ACGGTGCA", gapOpen: -2, gapExtend: -1);
        Assert.That((viaProfile.AlignedSequence1, viaProfile.AlignedSequence2, viaProfile.Score), Is.EqualTo((plain.AlignedSequence1, plain.AlignedSequence2, plain.Score)));
        // matrix + fractional scores: exactScore is a multiple of 1/4 here
        var m = AlignmentTools.AlignGapProfile("HEAGAWGHEE", "PAWHEAEHEE", gapOpen: -10, gapExtend: -0.5, matrix: "BLOSUM62", gapScores: new() { ["end_extend_gap_score"] = -0.25, ["end_open_gap_score"] = -1.5 });
        Assert.That(Math.Abs(m.ExactScore * 4 - Math.Round(m.ExactScore * 4)), Is.LessThan(1e-9));
        // local: only the internal scores count
        var loc = AlignmentTools.AlignGapProfile("TTACGTACGTTT", "GGACGTACGGG", "local", gapScores: new() { ["left_gap_score"] = 0, ["right_gap_score"] = 0 });
        Assert.That(loc.AlignmentType, Is.EqualTo("Local"));
    }

    [Test]
    public void AlignGapProfile_RejectsBadArguments()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignGapProfile("ACGT", "ACG", mode: "fitting"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignGapProfile("ACGT", "ACG", gapScores: new() { ["left_gap_scorex"] = 1 }));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignGapProfile("ACGT", "ACG", gapScores: new() { ["extend_gap_score"] = 1 }));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignGapProfile("ACGT", "ACG", endGaps: "FreeMiddle"));
    }
}
