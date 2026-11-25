using System.Globalization;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;
using static Seqeron.Genomics.Tests.Unit.Annotation.MiRnaTripletSvmReferenceData;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// MIRNA-CLASSIFY-001: published Triplet-SVM pre-miRNA classifier (Xue et al. 2005, BMC Bioinformatics 6:310) with the
/// authors' libsvm model. Oracles: the package's own example run (40 queries: 1.txt → predict_format.txt →
/// predict_result.txt), the package Perl scripts + libsvm 3.36 svm_predict_values on this library's folds, and a
/// 981-candidate random cross-check (300 ViennaRNA 2.8 folds × min arm 22/12/6, 0 differences; not committed).
/// </summary>
[TestFixture]
public class MiRnaAnalyzer_TripletSvm_Tests
{
    private static double[] ParseSparse(string line)
    {
        var v = new double[TripletSvmFeatureCount];
        foreach (string tok in line.Split(' ', StringSplitOptions.RemoveEmptyEntries))
        {
            var p = tok.Split(':');
            v[int.Parse(p[0], CultureInfo.InvariantCulture) - 1] = double.Parse(p[1], CultureInfo.InvariantCulture);
        }
        return v;
    }

    // TS1 — the 40 package examples: stem-loop span, 32 features (= shipped predict_format.txt), decision value (libsvm), label (= shipped predict_result.txt).
    [Test]
    public void TS1_ShippedPackageExamples_ReproduceFeaturesAndLabels()
    {
        Assert.That(Shipped, Has.Length.EqualTo(40));
        int positives = 0;
        foreach (var row in Shipped)
        {
            var res = ClassifyPreMiRnaTripletSvm(row.Sequence, row.Structure);
            Assert.That(res, Has.Count.EqualTo(1), row.Name);
            var r = res[0];
            Assert.That((r.Start, r.End), Is.EqualTo((row.Start, row.End)), row.Name);
            double[] expected = ParseSparse(row.FeatureLine);
            for (int k = 0; k < expected.Length; k++)
                Assert.That(r.Features[k], Is.EqualTo(expected[k]).Within(1e-12), $"{row.Name} feature {k + 1}");
            Assert.That(r.DecisionValue, Is.EqualTo(row.DecisionValue).Within(1e-9), row.Name);
            Assert.That(r.IsPreMiRna, Is.EqualTo(row.Label == 1), row.Name);
            if (r.IsPreMiRna) positives++;
        }
        Assert.That(positives, Is.EqualTo(37)); // hsa-mir-193b, -509, -510 are called negative by the package itself
    }

    // TS2 — Seqeron MFE folds of the 5 pri-miRNAs + 5 dinucleotide shuffles (Random(2024)): Perl + libsvm reference.
    [Test]
    public void TS2_SeqeronFolds_EqualPerlAndLibsvm([Values(22, 12)] int minBasePairs)
    {
        foreach (var row in SeqeronFolds)
        {
            var exp = minBasePairs == 22 ? row.Min22 : row.Min12;
            var got = ClassifyPreMiRnaTripletSvm(row.Sequence, row.Structure, minBasePairs);
            Assert.That(got, Has.Count.EqualTo(exp.Length), row.Name);
            for (int k = 0; k < exp.Length; k++)
            {
                Assert.That((got[k].Start, got[k].End), Is.EqualTo((exp[k].Start, exp[k].End)), row.Name);
                Assert.That(got[k].DecisionValue, Is.EqualTo(exp[k].DecisionValue).Within(1e-9), row.Name);
                Assert.That(got[k].IsPreMiRna, Is.EqualTo(exp[k].Label == 1), row.Name);
            }
        }
    }

    // TS3 — default folder = canonical MFE engine, so it equals the explicit-structure call on the same fold.
    [Test]
    public void TS3_DefaultFold_EqualsExplicitMfeStructure()
    {
        foreach (var row in SeqeronFolds.Take(5))
        {
            var viaFold = ClassifyPreMiRnaTripletSvm(row.Sequence);
            var viaStructure = ClassifyPreMiRnaTripletSvm(row.Sequence, row.Structure);
            Assert.That(viaFold, Has.Count.EqualTo(viaStructure.Count));
            for (int k = 0; k < viaFold.Count; k++)
                Assert.That(viaFold[k].DecisionValue, Is.EqualTo(viaStructure[k].DecisionValue));
            Assert.That(viaFold, Is.Not.Empty, row.Name);
            Assert.That(viaFold[0].IsPreMiRna, Is.True, row.Name); // all five miRBase precursors are called real
        }
    }

    // TS4 — model evaluation (libsvm svm_predict_values on the package model) at vectors outside the training data.
    [Test]
    public void TS4_DecisionValue_EqualsLibsvm()
    {
        var uniform = Enumerable.Repeat(1.0 / 32, 32).ToArray();
        var unit0 = new double[32]; unit0[0] = 1.0;
        var two = new double[32]; two[7] = 0.5; two[23] = 0.5;
        Assert.That(TripletSvmDecisionValue(uniform), Is.EqualTo(-3.1988414255639706).Within(1e-12));
        Assert.That(TripletSvmDecisionValue(unit0), Is.EqualTo(-2.7426894505258965).Within(1e-12));
        Assert.That(TripletSvmDecisionValue(two), Is.EqualTo(-2.6955952771281178).Within(1e-12));
    }

    // TS5 — input handling: DNA / lower case normalised, other characters dropped (package step 1), short stems give no candidate.
    [Test]
    public void TS5_Normalisation_AndGuards()
    {
        var row = Shipped[0];
        var dna = ClassifyPreMiRnaTripletSvm(row.Sequence.Replace('U', 'T').ToLowerInvariant(), row.Structure);
        Assert.That(dna, Has.Count.EqualTo(1));
        Assert.That(dna[0].DecisionValue, Is.EqualTo(row.DecisionValue).Within(1e-9));
        Assert.That(ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "((((....))))"), Is.Empty);          // 4 bp < 22
        Assert.That(ClassifyPreMiRnaTripletSvm("GGGGNAAACCCC", "((((....))))", 3), Is.Empty);       // non-ACGU query is deleted
        Assert.That(ClassifyPreMiRnaTripletSvm("", null), Is.Empty);
        Assert.That(ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "((((....))))", 3), Has.Count.EqualTo(1));
        Assert.Throws<ArgumentNullException>(() => ClassifyPreMiRnaTripletSvm(null!));
        Assert.Throws<ArgumentException>(() => ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "((((...))))"));
        Assert.Throws<ArgumentException>(() => ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "[[[[....]]]]"));
        Assert.Throws<ArgumentException>(() => TripletSvmDecisionValue(new double[5]));
    }

    // TS5b — argument validation (WP-R2-2c, exposed by the MCP wrapper): a stem arm needs >= 1 paired base, so
    // minBasePairs < 1 is rejected (with 0 the package's `< min` test would admit every unpaired piece); a negative
    // minLoopSize is meaningless (0..2 are raised to the NNDB minimum 3 by the folder, as documented).
    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void TS5b_MinBasePairs_BelowOne_Throws(int minBasePairs)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "((((....))))", minBasePairs));
        Assert.That(ex!.ParamName, Is.EqualTo("minBasePairs"));
        Assert.That(ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", "((((....))))", 1), Has.Count.EqualTo(1));
    }

    [Test]
    public void TS5c_MinLoopSize_Negative_Throws_ZeroToTwoAccepted()
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", null, 4, -1));
        Assert.That(ex!.ParamName, Is.EqualTo("minLoopSize"));
        Assert.DoesNotThrow(() => ClassifyPreMiRnaTripletSvm("GGGGAAAACCCC", null, 4, 0));
    }

    // TS6 — hand-computed triplet coding: stem GGGG/CCCC, loop AAAA: left arm "((((" → G at i=0 ".((" = 011 (3), inner "(((" = 7, last "((." = 110 (6);
    // right arm "))))" likewise for C. Index = nt(A0,G1,C2,U3)·8 + bits; 8 nt → each count/8.
    [Test]
    public void TS6_TripletCoding_HandComputed()
    {
        var f = ComputeTripletFeatures("GGGGAAAACCCC", "((((....))))");
        Assert.That(f, Is.Not.Null);
        var expected = new double[32];
        expected[1 * 8 + 3] = 1 / 8.0; expected[1 * 8 + 7] = 2 / 8.0; expected[1 * 8 + 6] = 1 / 8.0;
        expected[2 * 8 + 3] = 1 / 8.0; expected[2 * 8 + 7] = 2 / 8.0; expected[2 * 8 + 6] = 1 / 8.0;
        Assert.That(f, Is.EqualTo(expected).Within(1e-15));
        Assert.That(ComputeTripletFeatures("GGGGAAAACCCCGGGGAAAACCCC", "((((....))))((((....))))"), Is.Null); // two arm pairs: package skips
    }
}
