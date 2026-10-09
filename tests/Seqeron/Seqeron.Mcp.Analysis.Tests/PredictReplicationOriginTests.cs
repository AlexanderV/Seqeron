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
    [Test]
    public void PredictReplicationOrigin_AllPositionsAndCircular()
    {
        const string ba1f =
            "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG";
        var r = AnalysisTools.PredictReplicationOrigin(ba1f);
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
}
