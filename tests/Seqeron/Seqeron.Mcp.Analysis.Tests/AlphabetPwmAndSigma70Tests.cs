using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// create_alphabet_pwm / scan_with_alphabet_pwm / find_sigma70_promoters / predict_sigma70_promoters — delegation to
/// <c>MotifFinder.CreateAlphabetPwm</c>, <c>CalculateAlphabetPwmScores</c>, <c>ScanWithAlphabetPwm</c>,
/// <c>FindSigma70Promoters</c>, <c>PredictSigma70Promoters</c>. Expected values: Biopython 1.88 (protein PWM) and the
/// Promoter Calculator v1.0 reference Python, locked in MotifFinder_AlphabetPwm_Tests / MotifFinder_Sigma70Promoters_Tests.
/// </summary>
[TestFixture]
public class AlphabetPwmAndSigma70Tests
{
    private const string Protein = "ACDEFGHIKLMNPQRSTVWY";
    private static readonly string[] Instances = { "MKVLAT", "MRVLGT", "MKILAS", "LKVLAT", "MKVIAT", "MEVLAT" };

    private const string LacUv5 =
        "TCAGCATTCGAGCTTACGGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCT";

    [Test]
    public void CreateAlphabetPwm_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CreateAlphabetPwm(Instances, Protein, 0.5));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreateAlphabetPwm(Array.Empty<string>(), Protein));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreateAlphabetPwm(Instances, ""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreateAlphabetPwm(Instances, Protein, pseudocounts: new double[3]));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CreateAlphabetPwm(Instances, Protein, background: new double[4]));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.CreateAlphabetPwm(Instances, Protein, -1));
    }

    [Test]
    public void CreateAlphabetPwm_Binding_EqualsBiopython()
    {
        var r = AnalysisTools.CreateAlphabetPwm(Instances, Protein, 0.5);
        Assert.Multiple(() =>
        {
            Assert.That(r.Matrix, Has.Length.EqualTo(20));
            Assert.That(r.Matrix[Protein.IndexOf('M')][0], Is.EqualTo(2.7813597135246595).Within(1e-12));
            Assert.That(r.Consensus, Is.EqualTo("MKVLAT"));
            Assert.That(r.Anticonsensus, Is.EqualTo("AAAACA"));
            Assert.That(r.MaxScore, Is.EqualTo(16.398651663952972).Within(1e-12));
            Assert.That(r.MinScore, Is.EqualTo(-4.068431430675826).Within(1e-12));
            Assert.That(r.Mean, Is.EqualTo(3.809139711610947).Within(1e-12));
            Assert.That(r.Std, Is.EqualTo(3.8318293687371425).Within(1e-12));
        });

        var zero = AnalysisTools.CreateAlphabetPwm(Instances, Protein);
        Assert.That(zero.MinScore, Is.Null, "−∞ has no JSON form");
        Assert.That(zero.Matrix[Protein.IndexOf('W')][0], Is.Null);
    }

    [Test]
    public void ScanWithAlphabetPwm_Binding_ScoresNaNWindowsAndHits()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            AnalysisTools.ScanWithAlphabetPwm("MKVLAT", Instances, Protein, double.NegativeInfinity));

        var r = AnalysisTools.ScanWithAlphabetPwm("GGMKVLATxxMRvLGTPPMKXLAS", Instances, Protein, 5.0, 0.5);
        Assert.Multiple(() =>
        {
            Assert.That(r.Scores, Has.Length.EqualTo(19));
            Assert.That(r.Scores[2], Is.EqualTo(16.398651663952972).Within(1e-12));
            Assert.That(r.Scores[3], Is.Null);
            Assert.That(r.InvalidWindows, Is.EqualTo(new[] { 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 18 }));
            Assert.That(r.Hits.Select(h => h.Position), Is.EqualTo(new[] { 2, 10 }));
            Assert.That(r.Hits[1].Score, Is.EqualTo(12.939220045315675).Within(1e-12));
        });
    }

    [Test]
    public void FindSigma70Promoters_SchemaAndBinding()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindSigma70Promoters(""));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindSigma70Promoters("ACGT", maxMismatches35: 7));

        var r = AnalysisTools.FindSigma70Promoters(LacUv5, 1, 0);
        var h = r.Items.Single();
        Assert.That((h.Strand, h.Minus35Start, h.Minus35, h.Minus10Start, h.Minus10, h.Spacer, h.TotalMismatches),
            Is.EqualTo(("+", 68, "TTTACA", 92, "TATAAT", 18, 1)));
    }

    [Test]
    public void PredictSigma70Promoters_Binding_EqualsReference()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictSigma70Promoters("ACGN"));

        var r = AnalysisTools.PredictSigma70Promoters(LacUv5, bothStrands: false);
        Assert.That(r.Items, Has.Length.EqualTo(65));
        var best = r.Best!;
        Assert.Multiple(() =>
        {
            Assert.That(best.Tss, Is.EqualTo(108));
            Assert.That((best.Minus35, best.Minus10, best.Spacer.Length), Is.EqualTo(("TTTACA", "TATAAT", 18)));
            Assert.That(best.DeltaGTotal, Is.EqualTo(-3.044998584365458));
            Assert.That(best.TranscriptionRate, Is.EqualTo(6123.861140216731).Within(1e-9));
        });

        Assert.That(AnalysisTools.PredictSigma70Promoters(new string('A', 77)).Best, Is.Null);
    }
}
