// SPDX-License-Identifier: AGPL-3.0-only
// Tests of the port of RSAT oligo-analysis / dyad-analysis; part of Seqeron.Genomics.OverrepresentedWords.Tests, licensed under the GNU Affero General Public License, version 3 (see LICENSE in the production project).
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

    [TestCase(11.70753026510891, 3.033617247605823e-06)]
    [TestCase(4.77915994208994, 0.0056047467361406205)]
    [TestCase(0.0, 0.06962309429876835)]
    [TestCase(-100.0, 1.0)]
    public void Markov1_Wikipedia_EqualsExhaustiveEnumeration(double score, double expected)
    {
        var bg = OligoBackgroundModel.MarkovFromOligoFrequencies(Count(MarkovText, 2));
        AssertExact(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(WikipediaPwm(), score, bg), expected, 1e-13);
    }

    [TestCase(1e-2, 3.8616221022819137, 0.009851341698078267)]   // next lower word 3.8616221022819133: P = 0.01003988831187636
    [TestCase(1e-3, 6.898713577002425, 0.0009711428119647605)]   // next 6.849549269981338: P = 0.0011260438177962255
    [TestCase(1e-4, 9.1714773648687, 9.836356325619014e-05)]     // next 9.01120665544005: P = 1.163138875158734e-4
    public void Markov1_Threshold_EqualsExhaustiveEnumeration(double p, double score, double pv)
    {
        var bg = OligoBackgroundModel.MarkovFromOligoFrequencies(Count(MarkovText, 2));
        var r = OverrepresentedWordAnalyzer.PwmMarkovScoreThresholdForPValue(WikipediaPwm(), p, bg);
        Assert.That(r.Score, Is.EqualTo(score));
        AssertExact(r, pv, 1e-13);
    }

    [Test]
    [Description("Order 2, strand-insensitive table (RSAT 2str), ψ = 0.01 — exhaustive enumeration.")]
    public void Markov2StrandInsensitive_EqualsExhaustiveEnumeration()
    {
        var bg = OligoBackgroundModel.MarkovFromOligoFrequencies(Count(MarkovText, 3), 0.01, strandInsensitive: true);
        var pwm = WikipediaPwm();
        AssertExact(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 11.70753026510891, bg), 4.124380700286126e-07, 1e-13);
        AssertExact(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 4.77915994208994, bg), 0.008893066420964513, 1e-13);
        var t = OverrepresentedWordAnalyzer.PwmMarkovScoreThresholdForPValue(pwm, 1e-4, bg);
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
        AssertExact(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, -100, bg), 0.6289127081283757, 1e-13);
        AssertExact(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 11.70753026510891, bg), 1.3020833333333334e-05, 1e-13);
        AssertExact(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 0.0, bg), 0.07390705340952931, 1e-13);
        var t = OverrepresentedWordAnalyzer.PwmMarkovScoreThresholdForPValue(pwm, 1e-4, bg);
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
        AssertExact(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, threshold, bg), ape, 1e-14);
    }

    [Test]
    public void BernoulliModels_EqualTheIidOverload()
    {
        var pwm = WikipediaPwm();
        double[] q = { 0.3, 0.2, 0.2, 0.3 };
        Assert.That(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 4.77915994208994, OligoBackgroundModel.Bernoulli(q)),
            Is.EqualTo(MotifFinder.PwmScorePValue(pwm, 4.77915994208994, q)));
        Assert.That(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 4.77915994208994, OligoBackgroundModel.Equiprobable),
            Is.EqualTo(MotifFinder.PwmScorePValue(pwm, 4.77915994208994)));
        Assert.That(OverrepresentedWordAnalyzer.PwmMarkovScoreThresholdForPValue(pwm, 1e-3, OligoBackgroundModel.Bernoulli(q)),
            Is.EqualTo(MotifFinder.PwmScoreThresholdForPValue(pwm, 1e-3, q)));
    }

    [Test]
    public void Markov_UnsupportedModelsAndArguments_Throw()
    {
        var pwm = WikipediaPwm();
        Assert.Throws<ArgumentException>(() => OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 0, OligoBackgroundModel.MarkovFromInput(1)));
        Assert.Throws<ArgumentException>(() => OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 0, OligoBackgroundModel.BernoulliFromInput));
        Assert.Throws<ArgumentException>(() => OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 0, OligoBackgroundModel.Lexicon));
        Assert.Throws<ArgumentNullException>(() => OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, 0, null!));
        Assert.Throws<ArgumentNullException>(() => OverrepresentedWordAnalyzer.PwmMarkovScorePValue(null!, 0, OligoBackgroundModel.Equiprobable));
        Assert.Throws<ArgumentOutOfRangeException>(() => OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, double.NaN, OligoBackgroundModel.Equiprobable));
        Assert.Throws<ArgumentOutOfRangeException>(() => OverrepresentedWordAnalyzer.PwmMarkovScoreThresholdForPValue(pwm, 2, OligoBackgroundModel.Equiprobable));
    }


    private static PositionWeightMatrix WikipediaPwm() => MotifFinder.CreatePwm(WikipediaSequences, 0.25);

    private static void AssertExact(PwmPValueResult r, double expected, double relTol = 1e-12)
    {
        Assert.That(r.IsExact, Is.True, "exact");
        Assert.That(r.PValue, Is.EqualTo(expected).Within(Math.Max(expected * relTol, 1e-300)), "p-value");
        Assert.That(r.PValueLowerBound, Is.EqualTo(r.PValue));
        Assert.That(r.PValueUpperBound, Is.EqualTo(r.PValue));
    }

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

    private static readonly string[] WikipediaSequences =
    {
        "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
        "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
    };
}
