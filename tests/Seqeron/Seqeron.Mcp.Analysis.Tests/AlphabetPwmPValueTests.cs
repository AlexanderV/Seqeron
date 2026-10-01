using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// alphabet_pwm_score_pvalue / alphabet_pwm_score_thresholds — delegation to <c>MotifFinder.AlphabetPwmScorePValue</c>,
/// <c>AlphabetPwmScoreThresholdForPValue</c> and <c>AlphabetPositionWeightMatrix.ScoreDistribution</c>. Expected values:
/// exhaustive enumeration of all 20^6 words and Biopython 1.88 <c>pssm.distribution()</c>, locked in
/// MotifFinder_PwmPValueOptions_Tests.
/// </summary>
[TestFixture]
public class AlphabetPwmPValueTests
{
    private const string Protein = "ACDEFGHIKLMNPQRSTVWY";
    private static readonly string[] Instances = { "MKVLAT", "MRVLGT", "MKILAS", "LKVLAT", "MKVIAT", "MEVLAT" };

    [Test]
    public void AlphabetPwmScorePValue_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.AlphabetPwmScorePValue(Instances, Protein, score: 0.0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlphabetPwmScorePValue(Instances, Protein));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlphabetPwmScorePValue(Instances, Protein, 0.0, 0.01));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlphabetPwmScorePValue(Instances, Protein, score: 0.0, pseudocount: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.AlphabetPwmScorePValue(Instances, Protein, pValue: 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlphabetPwmScorePValue(Instances, Protein, score: 0, background: new[] { 0.5, 0.5 }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.AlphabetPwmScoreThresholds(Instances, Protein, pseudocount: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.AlphabetPwmScoreThresholds(Instances, Protein, fpr: 2));
    }

    [Test]
    public void AlphabetPwmScorePValue_Binding_EqualsExhaustiveEnumeration()
    {
        var consensus = AnalysisTools.AlphabetPwmScorePValue(Instances, Protein, score: 16.398651663952972);
        Assert.That(consensus.PValue, Is.EqualTo(1.5625000000000006e-08).Within(1e-21));
        Assert.That(consensus.IsExact, Is.True);

        var t = AnalysisTools.AlphabetPwmScorePValue(Instances, Protein, pValue: 1e-5);
        Assert.That(t.Score, Is.EqualTo(11.064750927399535).Within(1e-12));
        Assert.That(t.PValue, Is.EqualTo(8.906250000000003e-06).Within(1e-18));

        var none = AnalysisTools.AlphabetPwmScorePValue(Instances, Protein, pValue: 1e-9);
        Assert.That(none.ScoreAboveMaximum, Is.True);
        Assert.That(none.Score, Is.Null);
    }

    [Test]
    public void AlphabetPwmScoreThresholds_Binding_EqualsBiopython()
    {
        var r = AnalysisTools.AlphabetPwmScoreThresholds(Instances, Protein);
        Assert.Multiple(() =>
        {
            Assert.That(r.MinScore, Is.EqualTo(-4.068431430675826).Within(1e-15));
            Assert.That(r.Step, Is.EqualTo(0.0034117491406282377).Within(1e-17));
            Assert.That(r.PointCount, Is.EqualTo(6000));
            Assert.That(r.ThresholdFpr, Is.EqualTo(2.843772328236984).Within(1e-12));
            Assert.That(r.ThresholdFnr, Is.EqualTo(-0.6157413003600496).Within(1e-12));
            Assert.That(r.ThresholdBalanced, Is.EqualTo(-0.6157413003600496).Within(1e-12));
            Assert.That(r.ThresholdPatser, Is.EqualTo(0.970722050032081).Within(1e-12));
        });
    }
}
