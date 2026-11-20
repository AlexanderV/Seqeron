using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// B09 A3-8 / A3-12 / A3-13: MCP tools of the EMBOSS program ports. Expected values are those printed by the installed EMBOSS 6.6.0 programs (also locked by the oracle tests of
/// the library).
/// </summary>
[TestFixture]
public class EmbossProgramToolsTests
{
    private const string TandemSeq = "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACACGTTAACGTTAACGTTAACGTTAACGTTAacgtNN";

    [Test]
    public void FindTandemRepeatsEmboss_Etandem_ReportsScoreIdentityAndConsensus()
    {
        // etandem -minrepeat 2 -maxrepeat 8 -threshold 10: "33 61 96 3 12 100.0 cag" and "26 137 172 6 6 94.4 acgtta"
        var r = AnalysisTools.FindTandemRepeatsEmboss(TandemSeq, "etandem", 2, 8, 10);
        Assert.That(r.Program, Is.EqualTo("etandem"));
        Assert.That(r.Repeats.Select(h => (h.Position + 1, h.Position + h.Length, h.Score, h.RepeatSize, h.RepeatCount, h.Identity, h.Consensus)),
            Is.EqualTo(new[] { (61, 96, 33, 3, 12, (double?)100.0, "cag"), (137, 172, 26, 6, 6, (double?)94.4, "acgtta") }));
    }

    [Test]
    public void FindTandemRepeatsEmboss_Equicktandem_ReportsDiagonalScores()
    {
        // equicktandem -maxrepeat 6 -threshold 10 on (ACG)x8: "21 1 24 3 8"
        var r = AnalysisTools.FindTandemRepeatsEmboss("ACGACGACGACGACGACGACGACG", "EquickTandem", maxRepeat: 6, threshold: 10);
        var h = r.Repeats.Single();
        Assert.That((h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount, h.Identity, h.Consensus), Is.EqualTo((0, 24, 21, 3, 8, (double?)null, (string?)null)));
    }

    [Test]
    public void FindTandemRepeatsEmboss_RejectsBadArguments()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindTandemRepeatsEmboss("ACGT", "tandem"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindTandemRepeatsEmboss("", "etandem"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindTandemRepeatsEmboss("ACGT", "etandem", 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindTandemRepeatsEmboss("ACGT", "equicktandem", threshold: 0));
    }

    [Test]
    public void SeqMatchAll_ReportsMaximalExactMatchesOfEveryPair()
    {
        // seqmatchall -wordsize 6 on these three sequences: "10 a 1..10 b 3..12", "6 a 1..6 b 7..12", "6 a 5..10 b 3..8" (1-based); the third is shorter than the word size
        var r = AnalysisTools.SeqMatchAll(["ACGTACGTAC", "GGACGTACGTACTT", "TTTT"], 6);
        Assert.That(r.Matches.Select(m => (m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length)).ToArray(),
            Is.EqualTo(new[] { (0, 1, 0, 2, 10), (0, 1, 0, 6, 6), (0, 1, 4, 2, 6) }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SeqMatchAll([], 4));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SeqMatchAll(["ACGT", ""], 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.SeqMatchAll(["ACGT", "ACGT"], 0));
    }
}
