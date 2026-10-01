using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// PAT-PWM-001 — exact PWM p-values beyond the i.i.d. DNA case (B05 F32): search options (TFM-Pvalue
/// initial / maximal granularity, decrease factor, budgets, exhaustive mode), order-m Markov backgrounds
/// (RSAT <c>MarkovFromOligoFrequencies</c>) and generic-alphabet (protein) PWMs, plus the K-row
/// Biopython <c>pssm.distribution</c>.
/// Oracles: exhaustive enumeration of all words (Python / C, left-to-right double window sums, word probability
/// = RSAT segment_proba P(prefix)·∏P(b | context)); a meet-in-the-middle C enumeration of all 4^20 words of the
/// 20-column matrix (exact recomputation of every pair within 1e-9 of the threshold); MACRO-APE 3.0.6
/// <c>ru.autosome.ape.di.FindPvalue --from-mono -d 16</c> with a 16-value dinucleotide background (dyadic matrix,
/// so its ceil discretisation is exact); Biopython 1.88 <c>pssm.distribution(background, precision)</c> thresholds.
/// </summary>
[TestFixture]
public class MotifFinder_PwmPValueOptions_Tests
{
    private static readonly string[] WikipediaSequences =
    {
        "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
        "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
    };

    private static PositionWeightMatrix WikipediaPwm() => MotifFinder.CreatePwm(WikipediaSequences, 0.25);

    private static void AssertExact(PwmPValueResult r, double expected, double relTol = 1e-12)
    {
        Assert.That(r.IsExact, Is.True, "exact");
        Assert.That(r.PValue, Is.EqualTo(expected).Within(Math.Max(expected * relTol, 1e-300)), "p-value");
        Assert.That(r.PValueLowerBound, Is.EqualTo(r.PValue));
        Assert.That(r.PValueUpperBound, Is.EqualTo(r.PValue));
    }

    /// <summary>The 20-column matrix of MotifFinder_PwmPValue_Tests.PwmScorePValue_BudgetExhausted_ReturnsCertifiedBounds.</summary>
    private static PositionWeightMatrix Random20Pwm() => new PositionWeightMatrix(new[,]
    {
        { 0.3125902303080262, -0.04222823538929744, 0.06342538014163379, -6.022367813028454, 1.036220187912753, 1.2352164616940315, -0.02345897282398891, 0.5849625007211562, -5.426264754702098, -1.084888897586513, 0.5555187228286743, 0.9857861407802992, 0.874469117916141, -5.087462841250339, 0.5555187228286743, 1.0677446066358343, 0.5613112326565108, 0.6814704815745026, 0.44312980630819676, 0.2606517545227997 },
        { -0.9098021910284219, -0.22470628717469449, -0.13588342808177303, 0.6918777046376682, -1.024662054234269, -1.1926450779423958, 0.8139880143900512, -1.8231222379159209, 1.1128940564059333, -1.7369655941662063, -0.4296842752432452, -2.502500340529183, 0.4405725913859814, 1.6807214835265871, -1.4005379295837288, -1.308122295362332, -1.713695814843359, 0.6119295483214257, 0.44312980630819676, 0.2606517545227997 },
        { 0.5849625007211562, 0.21842351913350247, 0.39463128861700686, -0.18947779886371285, 0.45720695352278334, -1.8845227825800641, -1.0962153152593033, 0.6662626028230043, -0.7824085649273732, 0.7410817026384381, -1.6520766965796931, -1.5849625007211563, -1.304854581528421, -5.087462841250339, 0.874469117916141, -0.49098635251214234, -0.24902754783991468, -1.0840642647884744, 0.0, 0.16196747966339312 },
        { -0.46234321405720047, 0.013805799525030428, -0.4533656179379432, 0.5775450291586732, -3.560714954474479, -0.09310940439148147, -0.33948646627166706, -0.6655809609294407, 0.3016556998611012, 0.63890130783196, 0.5555187228286743, 0.5943611987234058, -1.304854581528421, -0.4436066514756146, -1.6520766965796931, -0.3428877135230085, 0.4626269577971039, -1.3356030317844387, -1.831877241191673, -1.0435016386365865 },
    }, 20);

    #region Options

    [Test]
    [Description("Default options: the 20-column case is not resolved (bounds 0.058402542608746444 … 0.05840680768687889). " +
                 "With a 4× state budget it is exact and equals the meet-in-the-middle enumeration of all 4^20 words: " +
                 "0.058404839602189895 (no word within 1e-9 of the threshold).")]
    public void PwmScorePValue_LargerBudget_ResolvesBudgetExhaustedCase()
    {
        var pwm = Random20Pwm();
        Assert.That(MotifFinder.PwmScorePValue(pwm, 0.0, null, null).IsExact, Is.False, "default budget");
        var r = MotifFinder.PwmScorePValue(pwm, 0.0, null, new PwmPValueOptions { MaxStates = 1 << 23, MaxSuffixSet = 1 << 22 });
        AssertExact(r, 0.058404839602189895, 1e-13);
        Assert.That(r.PValue, Is.GreaterThanOrEqualTo(0.058402542608746444).And.LessThanOrEqualTo(0.05840680768687889),
            "inside the default certified bounds");
    }

    [Test]
    [Description("Inverse on the 20-column matrix at p = 1e-3 (default options): 7.189559567375462 with P = 9.99999972009391e-4; " +
                 "the next lower word score 7.1895595327135027 has P = 1.0000000256695785e-3 (meet-in-the-middle enumeration).")]
    public void PwmScoreThresholdForPValue_Random20_EqualsMeetInTheMiddleEnumeration()
    {
        var r = MotifFinder.PwmScoreThresholdForPValue(Random20Pwm(), 1e-3, null, null);
        Assert.That(r.Score, Is.EqualTo(7.189559567375462));
        AssertExact(r, 0.000999999972009391, 1e-12);
    }

    [Test]
    [Description("A state budget of 4 cannot resolve any band, so the result is a certified bound; the exhaustive mode resolves " +
                 "the same case exactly (values = exhaustive enumeration of the 4^9 words, MotifFinder_PwmPValue_Tests).")]
    public void Exhaustive_ResolvesWhatTheBudgetCannot()
    {
        var pwm = WikipediaPwm();
        var tiny = new PwmPValueOptions { MaxStates = 4, MaxSuffixSet = 4 };
        var bounded = MotifFinder.PwmScorePValue(pwm, 4.77915994208994, null, tiny);
        Assert.That(bounded.IsExact, Is.False);
        Assert.That(bounded.PValueLowerBound, Is.LessThanOrEqualTo(0.006160736083984375));
        Assert.That(bounded.PValueUpperBound, Is.GreaterThanOrEqualTo(0.006160736083984375));

        var exhaustive = tiny with { Exhaustive = true };
        AssertExact(MotifFinder.PwmScorePValue(pwm, 4.77915994208994, null, exhaustive), 0.006160736083984375);
        AssertExact(MotifFinder.PwmScorePValue(pwm, 9.20753026510891, null, exhaustive), 0.000102996826171875);

        Assert.That(MotifFinder.PwmScoreThresholdForPValue(pwm, 1e-3, null, tiny).IsExact, Is.False);
        foreach (var (p, score, pv) in new[]
                 {
                     (1e-2, 4.028050165603465, 0.009918212890625),
                     (1e-3, 7.150252343998389, 0.000980377197265625),
                     (1e-4, 9.31606123696847, 9.918212890625e-05),
                 })
        {
            var t = MotifFinder.PwmScoreThresholdForPValue(pwm, p, null, exhaustive);
            Assert.That(t.Score, Is.EqualTo(score).Within(1e-12), $"p = {p}");
            AssertExact(t, pv);
        }
    }

    [Test]
    public void DefaultOptions_AreBitIdenticalToTheThreeArgumentOverloads()
    {
        var pwm = WikipediaPwm();
        double[] bg = { 0.3, 0.2, 0.2, 0.3 };
        foreach (double s in new[] { 11.70753026510891, 4.77915994208994, 0.0, -3.0 })
        {
            Assert.That(MotifFinder.PwmScorePValue(pwm, s, bg, PwmPValueOptions.Default), Is.EqualTo(MotifFinder.PwmScorePValue(pwm, s, bg)));
            Assert.That(MotifFinder.PwmScorePValue(pwm, s, bg, null), Is.EqualTo(MotifFinder.PwmScorePValue(pwm, s, bg)));
            Assert.That(MotifFinder.PwmScorePValue(pwm, s, bg, PwmPValueOptions.Exact), Is.EqualTo(MotifFinder.PwmScorePValue(pwm, s, bg)));
        }
        foreach (double p in new[] { 1e-2, 1e-4 })
            Assert.That(MotifFinder.PwmScoreThresholdForPValue(pwm, p, bg, PwmPValueOptions.Default),
                Is.EqualTo(MotifFinder.PwmScoreThresholdForPValue(pwm, p, bg)));
    }

    [Test]
    [Description("TFM-Pvalue granularity schedule: initial granularity 0.01 with decrease factor 100 tries g = 100, 10^4, …; " +
                 "a maximal granularity of 0.1 allows only g = 10.")]
    public void GranularitySchedule_FollowsTheOptions()
    {
        var pwm = WikipediaPwm();
        var coarse = new PwmPValueOptions { InitialGranularity = 0.01, DecreaseFactor = 100 };
        var r = MotifFinder.PwmScorePValue(pwm, 9.20753026510891, null, coarse);
        AssertExact(r, 0.000102996826171875);
        Assert.That(r.Granularity, Is.EqualTo(100).Or.EqualTo(1e4).Or.EqualTo(1e6));

        var onlyTen = MotifFinder.PwmScorePValue(Random20Pwm(), 0.0, null, new PwmPValueOptions { MaxGranularity = 0.1 });
        Assert.That(onlyTen.IsExact, Is.False);
        Assert.That(onlyTen.Granularity, Is.EqualTo(10));
    }

    [Test]
    public void InvalidOptions_Throw()
    {
        var pwm = WikipediaPwm();
        foreach (var o in new[]
                 {
                     new PwmPValueOptions { InitialGranularity = 0 },
                     new PwmPValueOptions { InitialGranularity = double.NaN },
                     new PwmPValueOptions { MaxGranularity = 0.5 },
                     new PwmPValueOptions { MaxGranularity = -1 },
                     new PwmPValueOptions { DecreaseFactor = 1 },
                     new PwmPValueOptions { MaxStates = 0 },
                     new PwmPValueOptions { MaxSuffixSet = 0 },
                 })
        {
            Assert.Throws<ArgumentException>(() => MotifFinder.PwmScorePValue(pwm, 1.0, null, o));
        }
    }

    #endregion

    #region Markov background

    private const string MarkovText = "ATGCGTACGTTAGCCGATAGGCTTACGATCGGATCCATGGCATTAGCAAGT";

    private static Dictionary<string, double> Count(string text, int k)
    {
        var d = new Dictionary<string, double>();
        for (int i = 0; i + k <= text.Length; i++)
        {
            string w = text.Substring(i, k);
            d[w] = d.GetValueOrDefault(w) + 1;
        }
        return d;
    }

    [TestCase(11.70753026510891, 3.033617247605823e-06)]
    [TestCase(4.77915994208994, 0.0056047467361406205)]
    [TestCase(0.0, 0.06962309429876835)]
    [TestCase(-100.0, 1.0)]
    public void Markov1_Wikipedia_EqualsExhaustiveEnumeration(double score, double expected)
    {
        var bg = OligoBackgroundModel.MarkovFromOligoFrequencies(Count(MarkovText, 2));
        AssertExact(MotifFinder.PwmMarkovScorePValue(WikipediaPwm(), score, bg), expected, 1e-13);
    }

    [TestCase(1e-2, 3.8616221022819137, 0.009851341698078267)]   // next lower word 3.8616221022819133: P = 0.01003988831187636
    [TestCase(1e-3, 6.898713577002425, 0.0009711428119647605)]   // next 6.849549269981338: P = 0.0011260438177962255
    [TestCase(1e-4, 9.1714773648687, 9.836356325619014e-05)]     // next 9.01120665544005: P = 1.163138875158734e-4
    public void Markov1_Threshold_EqualsExhaustiveEnumeration(double p, double score, double pv)
    {
        var bg = OligoBackgroundModel.MarkovFromOligoFrequencies(Count(MarkovText, 2));
        var r = MotifFinder.PwmMarkovScoreThresholdForPValue(WikipediaPwm(), p, bg);
        Assert.That(r.Score, Is.EqualTo(score));
        AssertExact(r, pv, 1e-13);
    }

    [Test]
    [Description("Order 2, strand-insensitive table (RSAT 2str), ψ = 0.01 — exhaustive enumeration.")]
    public void Markov2StrandInsensitive_EqualsExhaustiveEnumeration()
    {
        var bg = OligoBackgroundModel.MarkovFromOligoFrequencies(Count(MarkovText, 3), 0.01, strandInsensitive: true);
        var pwm = WikipediaPwm();
        AssertExact(MotifFinder.PwmMarkovScorePValue(pwm, 11.70753026510891, bg), 4.124380700286126e-07, 1e-13);
        AssertExact(MotifFinder.PwmMarkovScorePValue(pwm, 4.77915994208994, bg), 0.008893066420964513, 1e-13);
        var t = MotifFinder.PwmMarkovScoreThresholdForPValue(pwm, 1e-4, bg);
        Assert.That(t.Score, Is.EqualTo(9.250117814108265)); // next lower word 9.1714773648687: P = 1.5779878428662038e-4
        AssertExact(t, 4.7495286790552905e-06, 1e-13);
    }

    [Test]
    [Description("ψ = 0 and no successor of C in the table: words through C have probability 0 (RSAT segment_proba), the word " +
                 "measure has total mass 0.6289127081283757 and P(S ≥ min) equals it — exhaustive enumeration.")]
    public void SubStochasticMarkovTable_UsesTheWordMeasure()
    {
        var sub = new Dictionary<string, double>
        {
            ["AA"] = 3, ["AC"] = 2, ["AG"] = 4, ["AT"] = 1, ["GA"] = 2, ["GG"] = 5, ["GT"] = 1, ["TA"] = 2, ["TT"] = 3, ["TG"] = 1,
        };
        var bg = OligoBackgroundModel.MarkovFromOligoFrequencies(sub, 0.0);
        var pwm = WikipediaPwm();
        AssertExact(MotifFinder.PwmMarkovScorePValue(pwm, -100, bg), 0.6289127081283757, 1e-13);
        AssertExact(MotifFinder.PwmMarkovScorePValue(pwm, 11.70753026510891, bg), 1.3020833333333334e-05, 1e-13);
        AssertExact(MotifFinder.PwmMarkovScorePValue(pwm, 0.0, bg), 0.07390705340952931, 1e-13);
        var t = MotifFinder.PwmMarkovScoreThresholdForPValue(pwm, 1e-4, bg);
        Assert.That(t.Score, Is.EqualTo(9.941995518745932)); // next lower possible word 9.385602170221548: P = 2.3423936631944444e-4
        AssertExact(t, 9.765625e-05, 1e-13);
    }

    [TestCase(2.0, 0.2654045414462081)]
    [TestCase(3.0, 0.12487874779541447)]
    [TestCase(4.5, 0.01626984126984127)]
    [Description("MACRO-APE 3.0.6 ru.autosome.ape.di.FindPvalue t.pat T --from-mono -d 16 -b <16 circular dinucleotide " +
                 "frequencies of ACGTTGCAAGCTTAGGCATCGATCCGATGGCTAAT> (first-letter = second-letter marginals, so MACRO-APE's " +
                 "initial distribution equals RSAT's).")]
    public void Markov1_EqualsMacroApeDinucleotideBackground(double threshold, double ape)
    {
        var pwm = new PositionWeightMatrix(new[,]
        {
            { 1.5, -1, 0.75, -0.5 },
            { -0.5, 2, -1.5, 0.5 },
            { 0.25, 0.5, 1, -0.75 },
            { -1, -0.25, 0.125, 1.25 },
        }, 4);
        var table = new Dictionary<string, double>
        {
            ["AA"] = 2, ["AC"] = 1, ["AG"] = 2, ["AT"] = 4, ["CA"] = 2, ["CC"] = 1, ["CG"] = 3, ["CT"] = 2,
            ["GA"] = 2, ["GC"] = 4, ["GG"] = 2, ["GT"] = 1, ["TA"] = 3, ["TC"] = 2, ["TG"] = 2, ["TT"] = 2,
        };
        var bg = OligoBackgroundModel.MarkovFromOligoFrequencies(table, 0.0);
        AssertExact(MotifFinder.PwmMarkovScorePValue(pwm, threshold, bg), ape, 1e-14);
    }

    [Test]
    public void BernoulliModels_EqualTheIidOverload()
    {
        var pwm = WikipediaPwm();
        double[] q = { 0.3, 0.2, 0.2, 0.3 };
        Assert.That(MotifFinder.PwmMarkovScorePValue(pwm, 4.77915994208994, OligoBackgroundModel.Bernoulli(q)),
            Is.EqualTo(MotifFinder.PwmScorePValue(pwm, 4.77915994208994, q)));
        Assert.That(MotifFinder.PwmMarkovScorePValue(pwm, 4.77915994208994, OligoBackgroundModel.Equiprobable),
            Is.EqualTo(MotifFinder.PwmScorePValue(pwm, 4.77915994208994)));
        Assert.That(MotifFinder.PwmMarkovScoreThresholdForPValue(pwm, 1e-3, OligoBackgroundModel.Bernoulli(q)),
            Is.EqualTo(MotifFinder.PwmScoreThresholdForPValue(pwm, 1e-3, q)));
    }

    [Test]
    public void Markov_UnsupportedModelsAndArguments_Throw()
    {
        var pwm = WikipediaPwm();
        Assert.Throws<ArgumentException>(() => MotifFinder.PwmMarkovScorePValue(pwm, 0, OligoBackgroundModel.MarkovFromInput(1)));
        Assert.Throws<ArgumentException>(() => MotifFinder.PwmMarkovScorePValue(pwm, 0, OligoBackgroundModel.BernoulliFromInput));
        Assert.Throws<ArgumentException>(() => MotifFinder.PwmMarkovScorePValue(pwm, 0, OligoBackgroundModel.Lexicon));
        Assert.Throws<ArgumentNullException>(() => MotifFinder.PwmMarkovScorePValue(pwm, 0, null!));
        Assert.Throws<ArgumentNullException>(() => MotifFinder.PwmMarkovScorePValue(null!, 0, OligoBackgroundModel.Equiprobable));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmMarkovScorePValue(pwm, double.NaN, OligoBackgroundModel.Equiprobable));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmMarkovScoreThresholdForPValue(pwm, 2, OligoBackgroundModel.Equiprobable));
    }

    #endregion

    #region Generic alphabet

    private const string Protein = "ACDEFGHIKLMNPQRSTVWY";
    private static readonly string[] ProteinInstances = { "MKVLAT", "MRVLGT", "MKILAS", "LKVLAT", "MKVIAT", "MEVLAT" };
    private static AlphabetPositionWeightMatrix ProteinPwm() => MotifFinder.CreateAlphabetPwm(ProteinInstances, Protein, 0.5);
    private static double[] ProteinBackground() => Enumerable.Range(0, 20).Select(i => 1.0 + (i % 5)).ToArray();

    [TestCase(16.398651663952972, 1.5625000000000006e-08, 2.5720164609053499e-09)]   // consensus MKVLAT
    [TestCase(12.939220045315675, 1.9531250000000007e-06, 5.4070216049382716e-07)]   // MRVLGT
    [TestCase(8.0, 0.00010960937500000005, 4.7479938271604937e-05)]
    [TestCase(0.0, 0.098773750000002737, 0.087475587191358164)]
    [Description("Exhaustive enumeration of all 20^6 words (C, long double accumulation), uniform and 1,2,3,4,5,… background.")]
    public void AlphabetPwmScorePValue_Protein_EqualsExhaustiveEnumeration(double score, double uniform, double weighted)
    {
        AssertExact(MotifFinder.AlphabetPwmScorePValue(ProteinPwm(), score), uniform, 1e-12);
        AssertExact(MotifFinder.AlphabetPwmScorePValue(ProteinPwm(), score, ProteinBackground()), weighted, 1e-12);
    }

    [TestCase(1e-3, 6.3098634252360659, 0.00091400000000000053, 6.0203568080410816, 0.00063985127314814807)]
    [TestCase(1e-5, 11.064750927399535, 8.906250000000003e-06, 9.7692950438733632, 9.6516203703703694e-06)]
    [TestCase(1e-7, 14.524182546036833, 6.2500000000000024e-08, 13.228726662510661, 5.9413580246913582e-08)]
    [Description("Exhaustive enumeration of the 20^6 words: smallest word score with P ≤ p; e.g. p = 1e-3 uniform: next lower " +
                 "word 6.309863425236065 has P = 0.0013443124999999995.")]
    public void AlphabetPwmScoreThresholdForPValue_Protein_EqualsExhaustiveEnumeration(
        double p, double score, double pv, double weightedScore, double weightedPv)
    {
        var r = MotifFinder.AlphabetPwmScoreThresholdForPValue(ProteinPwm(), p);
        Assert.That(r.Score, Is.EqualTo(score));
        AssertExact(r, pv, 1e-12);
        var w = MotifFinder.AlphabetPwmScoreThresholdForPValue(ProteinPwm(), p, ProteinBackground());
        Assert.That(w.Score, Is.EqualTo(weightedScore));
        AssertExact(w, weightedPv, 1e-12);
    }

    [Test]
    [Description("Biopython 1.88: motifs.create(instances, alphabet=protein).counts.normalize(pseudocounts=0.5).log_odds()" +
                 ".distribution(background, precision=1000) — threshold_fpr(0.01), threshold_fnr(0.1), threshold_balanced(), threshold_patser().")]
    public void AlphabetScoreDistribution_Protein_EqualsBiopython()
    {
        var pwm = ProteinPwm();
        var d = pwm.ScoreDistribution();
        Assert.Multiple(() =>
        {
            Assert.That(d.MinScore, Is.EqualTo(-4.068431430675826));
            Assert.That(d.Step, Is.EqualTo(0.0034117491406282377));
            Assert.That(d.PointCount, Is.EqualTo(6000));
            Assert.That(d.ThresholdFpr(0.01), Is.EqualTo(2.843772328236984));
            Assert.That(d.ThresholdFnr(0.1), Is.EqualTo(-0.6157413003600496));
            Assert.That(d.ThresholdBalanced(), Is.EqualTo(-0.6157413003600496));
            Assert.That(d.ThresholdPatser(), Is.EqualTo(0.970722050032081));
        });
        var w = pwm.ScoreDistribution(ProteinBackground());
        Assert.Multiple(() =>
        {
            Assert.That(w.ThresholdFpr(0.01), Is.EqualTo(2.843772328236984));
            Assert.That(w.ThresholdFnr(0.1), Is.EqualTo(-0.9023282281728213));
            Assert.That(w.ThresholdBalanced(), Is.EqualTo(-2.488791578564952));
            Assert.That(w.ThresholdPatser(), Is.EqualTo(-0.6157413003600496));
        });
    }

    [Test]
    [Description("For the alphabet ACGT the K-row engine and distribution equal the DNA PositionWeightMatrix ones.")]
    public void AlphabetAcgt_EqualsDnaEngine()
    {
        var apwm = MotifFinder.CreateAlphabetPwm(WikipediaSequences, "ACGT", 0.25);
        var dna = WikipediaPwm();
        double[] bg = { 0.3, 0.2, 0.2, 0.3 };
        foreach (double s in new[] { 11.70753026510891, 4.77915994208994, 0.0 })
            Assert.That(MotifFinder.AlphabetPwmScorePValue(apwm, s, bg).PValue, Is.EqualTo(MotifFinder.PwmScorePValue(dna, s, bg).PValue).Within(1e-15));
        var t1 = MotifFinder.AlphabetPwmScoreThresholdForPValue(apwm, 1e-3);
        var t2 = MotifFinder.PwmScoreThresholdForPValue(dna, 1e-3);
        Assert.That(t1.Score, Is.EqualTo(t2.Score).Within(1e-12));
        Assert.That(t1.PValue, Is.EqualTo(t2.PValue).Within(1e-15));
        Assert.That(apwm.ScoreDistribution(bg, 100).ThresholdFpr(0.01), Is.EqualTo(dna.ScoreDistribution(bg, 100).ThresholdFpr(0.01)).Within(1e-12));
    }

    [Test]
    public void AlphabetPValue_InvalidArguments_Throw()
    {
        var pwm = ProteinPwm();
        var withInf = MotifFinder.CreateAlphabetPwm(ProteinInstances, Protein); // zero pseudocount → −∞ cells
        Assert.Throws<ArgumentNullException>(() => MotifFinder.AlphabetPwmScorePValue(null!, 0));
        Assert.Throws<ArgumentException>(() => MotifFinder.AlphabetPwmScorePValue(withInf, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.AlphabetPwmScorePValue(pwm, double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.AlphabetPwmScoreThresholdForPValue(pwm, -0.5));
        Assert.Throws<ArgumentException>(() => MotifFinder.AlphabetPwmScorePValue(pwm, 0, new double[] { 0.5, 0.5 }));
        Assert.Throws<InvalidOperationException>(() => withInf.ScoreDistribution());
    }

    #endregion
}
