using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class TripletSvmDecisionValueTests
{
    [Test]
    public void TripletSvmDecisionValue_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.TripletSvmDecisionValue(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.TripletSvmDecisionValue(Array.Empty<double>()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.TripletSvmDecisionValue(new double[5]));
        var nan = new double[32]; nan[3] = double.NaN;
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.TripletSvmDecisionValue(nan));
    }

    [Test]
    public void TripletSvmDecisionValue_Binding_EqualsLibsvmLockedValues()
    {
        // Library TS4: libsvm svm_predict_values on the package model at vectors outside the training data.
        var uniform = Enumerable.Repeat(1.0 / 32, 32).ToArray();
        var unit0 = new double[32]; unit0[0] = 1.0;
        var two = new double[32]; two[7] = 0.5; two[23] = 0.5;
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.TripletSvmDecisionValue(uniform).DecisionValue, Is.EqualTo(-3.1988414255639706).Within(1e-12));
            Assert.That(AnnotationTools.TripletSvmDecisionValue(unit0).DecisionValue, Is.EqualTo(-2.7426894505258965).Within(1e-12));
            Assert.That(AnnotationTools.TripletSvmDecisionValue(two).DecisionValue, Is.EqualTo(-2.6955952771281178).Within(1e-12));
            Assert.That(AnnotationTools.TripletSvmDecisionValue(uniform).IsPreMiRna, Is.False);
            Assert.That(AnnotationTools.TripletSvmDecisionValue(two).DecisionValue, Is.EqualTo(MiRnaAnalyzer.TripletSvmDecisionValue(two)));
        });
    }

    [Test]
    public void TripletSvmDecisionValue_ChainedFromFeatures_MatchesClassifyTool()
    {
        // compute_triplet_features -> triplet_svm_decision_value reproduces classify_pre_mirna_triplet_svm's candidate decision value.
        const string mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";
        var cls = AnnotationTools.ClassifyPreMiRnaTripletSvm(mir21, null, 16);
        Assert.That(cls.Found, Is.True);
        var cand = cls.Candidates[0];
        var f = AnnotationTools.ComputeTripletFeatures(cand.Sequence, cand.DotBracket);
        Assert.That(f.Found, Is.True);
        var dv = AnnotationTools.TripletSvmDecisionValue(f.Features!.ToArray());
        Assert.That(dv.DecisionValue, Is.EqualTo(cand.DecisionValue).Within(1e-15));
        Assert.That(dv.IsPreMiRna, Is.EqualTo(cand.IsPreMiRna));
    }
}
