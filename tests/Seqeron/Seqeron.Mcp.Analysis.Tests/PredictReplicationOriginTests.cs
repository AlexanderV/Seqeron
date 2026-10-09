using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_replication_origin</c> MCP tool.
/// Expected values from GcSkewCalculator's own unit test
/// (GcSkewCalculator_PredictReplicationOrigin_Tests: "CCGGGG" origin 2/-2 terminus 6/+2;
/// "GGGCCC" terminus 3/+3 origin 0/0), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictReplicationOriginTests
{
    [Test]
    public void PredictReplicationOrigin_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictReplicationOrigin("CCGGGG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictReplicationOrigin(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictReplicationOrigin(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictReplicationOrigin("XYZ"));
    }

    [Test]
    public void PredictReplicationOrigin_Binding_InvokesSuccessfully()
    {
        // "CCGGGG": cumulative skew min -2 at prefix 2, max +2 at prefix 6.
        var r1 = AnalysisTools.PredictReplicationOrigin("CCGGGG");
        Assert.Multiple(() =>
        {
            Assert.That(r1.PredictedOrigin, Is.EqualTo(2));
            Assert.That(r1.OriginSkew, Is.EqualTo(-2.0).Within(1e-10));
            Assert.That(r1.PredictedTerminus, Is.EqualTo(6));
            Assert.That(r1.TerminusSkew, Is.EqualTo(2.0).Within(1e-10));
        });

        // "GGGCCC": max +3 at prefix 3, min 0 at prefix 0.
        var r2 = AnalysisTools.PredictReplicationOrigin("GGGCCC");
        Assert.Multiple(() =>
        {
            Assert.That(r2.PredictedTerminus, Is.EqualTo(3));
            Assert.That(r2.TerminusSkew, Is.EqualTo(3.0).Within(1e-10));
            Assert.That(r2.PredictedOrigin, Is.EqualTo(0));
            Assert.That(r2.OriginSkew, Is.EqualTo(0.0).Within(1e-10));
        });
    }

    // Finisher A1-2/A1-3: all BA1F minimizers (Rosalind sample answer "53 97"; maximizers [16, 20, 21]
    // by python brute force) and the circular flag ("GGGCCC": linear minimizers [0, 6] → circular [0]).
    private const string Ba1f =
        "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG";

    [Test]
    public void PredictReplicationOrigin_AllPositionsAndCircular()
    {
        var r = AnalysisTools.PredictReplicationOrigin(Ba1f);
        var lin = AnalysisTools.PredictReplicationOrigin("GGGCCC");
        var circ = AnalysisTools.PredictReplicationOrigin("GGGCCC", circular: true);
        var circG = AnalysisTools.PredictReplicationOrigin("G", circular: true);
        Assert.Multiple(() =>
        {
            Assert.That(r.OriginPositions, Is.EqualTo(new[] { 53, 97 }));
            Assert.That(r.TerminusPositions, Is.EqualTo(new[] { 16, 20, 21 }));
            Assert.That(r.PredictedOrigin, Is.EqualTo(53));
            Assert.That(lin.OriginPositions, Is.EqualTo(new[] { 0, 6 }));
            Assert.That(circ.OriginPositions, Is.EqualTo(new[] { 0 }));
            Assert.That(circ.TerminusPositions, Is.EqualTo(new[] { 3 }));
            Assert.That(circG.PredictedTerminus, Is.EqualTo(0));
            Assert.That(AnalysisTools.PredictReplicationOrigin("G").PredictedTerminus, Is.EqualTo(1));
        });
    }

    // Finisher A2-1: windowSize = Grigoriev windowed prediction. Reference: Biopython 1.88
    // numpy.cumsum(GC_skew(ba1f, w)[:n//w]) first argmin/argmax → window centre i·w + w/2:
    // w10 → 45 / −0.3095238095238095, 15 / 0.5; w25 → 37 / −0.15384615384615385, 62 / 0.11888111888111885.
    [TestCase(10, 45, -0.3095238095238095, 15, 0.5)]
    [TestCase(25, 37, -0.15384615384615385, 62, 0.11888111888111885)]
    public void PredictReplicationOrigin_WindowSize_MatchesBiopythonCumsum(
        int w, int origin, double originSkew, int terminus, double terminusSkew)
    {
        var r = AnalysisTools.PredictReplicationOrigin(Ba1f, windowSize: w);
        Assert.Multiple(() =>
        {
            Assert.That(r.PredictedOrigin, Is.EqualTo(origin));
            Assert.That(r.OriginSkew, Is.EqualTo(originSkew).Within(1e-12));
            Assert.That(r.PredictedTerminus, Is.EqualTo(terminus));
            Assert.That(r.TerminusSkew, Is.EqualTo(terminusSkew).Within(1e-12));
            Assert.That(r.IsSignificant, Is.True);
            Assert.That(r.OriginPositions, Is.Empty);
            Assert.That(r.TerminusPositions, Is.Empty);
            Assert.That(r.SkewIndex, Is.Null);
        });
    }

    [Test]
    public void PredictReplicationOrigin_WindowSize_Guards()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => AnalysisTools.PredictReplicationOrigin(Ba1f, circular: true, windowSize: 10));
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PredictReplicationOrigin(Ba1f, windowSize: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PredictReplicationOrigin(Ba1f, skewIndexWindow: 0));
            // Shorter than one window: no point → zero prediction, not significant (Core W3).
            var r = AnalysisTools.PredictReplicationOrigin("ACG", windowSize: 10);
            Assert.That((r.PredictedOrigin, r.PredictedTerminus, r.IsSignificant), Is.EqualTo((0, 0, false)));
        });
    }

    // Finisher A2-1: skewIndexWindow = SkewIT SkewI. Reference: the authors' skewi.py run with
    // `-k 4 --min-len 0` (headers containing "complete", %r output): 0.23076923076923078,
    // 0.6206896551724138, BA1F sample 0.16; 12 windows (r = round(0.48) = 0) → no output → null.
    [TestCase("GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGG", 0.23076923076923078)]
    [TestCase("GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAAT", 0.6206896551724138)]
    [TestCase(Ba1f, 0.16)]
    [TestCase("GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCC", null)]
    public void PredictReplicationOrigin_SkewIndex_MatchesSkewItScript(string seq, double? expected)
    {
        var r = AnalysisTools.PredictReplicationOrigin(seq, skewIndexWindow: 4);
        var plain = AnalysisTools.PredictReplicationOrigin(seq);
        Assert.Multiple(() =>
        {
            if (expected is null)
                Assert.That(r.SkewIndex, Is.Null);
            else
                Assert.That(r.SkewIndex, Is.EqualTo(expected.Value).Within(1e-15));
            Assert.That(plain.SkewIndex, Is.Null);
            Assert.That(r.OriginPositions, Is.EqualTo(plain.OriginPositions));
            Assert.That(AnalysisTools.PredictReplicationOrigin(seq, windowSize: 4, skewIndexWindow: 4).SkewIndex,
                Is.EqualTo(r.SkewIndex));
        });
    }
}
