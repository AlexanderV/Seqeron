using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// PAT-PWM-001 — exact PWM score p-values (Touzet &amp; Varré 2007, TFM-Pvalue) and the per-base / JASPAR
/// pseudocount <c>CreatePwm</c> overloads.
/// Oracles: exhaustive enumeration of all 4^L words (C program, left-to-right double window sums as
/// ScorePwmWindow, probabilities accumulated in long double) for the Wikipedia PWM (L = 9) and the Bucher
/// TATA box POL012.1 (L = 15, 4^15 words); TFM-Pvalue C++ (CRAN TFMPvalue 1.0.0 src/Matrix.cpp, driver
/// replicating testScoreToPvalue); Biopython 1.88 for the pseudocount PWMs.
/// </summary>
[TestFixture]
public class MotifFinder_PwmPValue_Tests
{
    private static readonly string[] WikipediaSequences =
    {
        "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
        "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
    };

    private const string WikipediaTarget = "CCTAGGTAAGTAACAGGTCAGTGG";
    private static readonly double[] NonUniform = { 0.3, 0.2, 0.2, 0.3 };

    private static PositionWeightMatrix WikipediaPwm() => MotifFinder.CreatePwm(WikipediaSequences, 0.25);

    /// <summary>3 columns, +2 for the column's own base (A, C, G), −1 otherwise: S = 3·matches − 3.</summary>
    private static PositionWeightMatrix ToyPwm() => new(new double[,]
    {
        { 2, -1, -1 },
        { -1, 2, -1 },
        { -1, -1, 2 },
        { -1, -1, -1 },
    }, 3);

    private static void AssertExact(PwmPValueResult r, double expected, double relTol = 1e-12)
    {
        Assert.That(r.IsExact, Is.True, "exact");
        Assert.That(r.PValue, Is.EqualTo(expected).Within(Math.Max(expected * relTol, 1e-300)), "p-value");
        Assert.That(r.PValueLowerBound, Is.EqualTo(r.PValue));
        Assert.That(r.PValueUpperBound, Is.EqualTo(r.PValue));
    }

    #region PwmScorePValue

    [TestCase(6.0, 1.0 / 64)]
    [TestCase(3.0, 10.0 / 64)]   // ≥ 2 matches: 3·(1/4)²(3/4) + 1/64
    [TestCase(2.5, 10.0 / 64)]
    [TestCase(0.0, 37.0 / 64)]   // ≥ 1 match: 1 − 27/64
    [TestCase(-3.0, 1.0)]
    [TestCase(6.5, 0.0)]
    public void PwmScorePValue_ToyMatrix_EqualsBinomialTail(double score, double expected)
    {
        AssertExact(MotifFinder.PwmScorePValue(ToyPwm(), score), expected);
    }

    [TestCase(11.70753026510891, 3.814697265625e-06)]   // consensus = observed site at 2 (4^-9)
    [TestCase(9.20753026510891, 0.000102996826171875)]
    [TestCase(5.853765132554455, 0.002780914306640625)]
    [TestCase(4.77915994208994, 0.006160736083984375)]  // observed site at 6
    [TestCase(0.0, 0.06875991821289062)]
    public void PwmScorePValue_Wikipedia_EqualsExhaustiveEnumeration(double score, double expected)
    {
        AssertExact(MotifFinder.PwmScorePValue(WikipediaPwm(), score), expected);
    }

    [TestCase(11.70753026510891, 5.832000000000001e-06)]
    [TestCase(4.77915994208994, 0.007509384)]
    public void PwmScorePValue_WikipediaBackground_EqualsExhaustiveEnumeration(double score, double expected)
    {
        AssertExact(MotifFinder.PwmScorePValue(WikipediaPwm(), score, NonUniform), expected, 1e-11);
    }

    [Test]
    [Description("The p-value of an observed window's score counts that window (TFM-Pvalue's rounding gives 0 for the consensus score).")]
    public void PwmScorePValue_ObservedWindowScores_AreCounted()
    {
        var pwm = WikipediaPwm();
        double[] scores = MotifFinder.CalculatePwmScores(WikipediaTarget, pwm);
        Assert.That(scores[2], Is.EqualTo(pwm.MaxScore));
        AssertExact(MotifFinder.PwmScorePValue(pwm, scores[2]), Math.Pow(0.25, 9));
        AssertExact(MotifFinder.PwmScorePValue(pwm, scores[6]), 0.006160736083984375);
        Assert.That(MotifFinder.PwmScorePValue(pwm, Math.BitIncrement(scores[2])).PValue, Is.Zero);
    }

    [TestCase(15.782082872405402, 9.313225746154785e-10)]   // consensus, 4^-15
    [TestCase(13.282082872405402, 6.019137799739838e-06)]
    [TestCase(7.891041436202701, 0.000597665086388588)]
    [TestCase(0.0, 0.02533565554767847)]
    public void PwmScorePValue_BucherTataBox_EqualsExhaustiveEnumeration(double score, double expected)
    {
        var pwm = MotifFinder.BucherPromoterMatrices.TataBox.Pwm;
        AssertExact(MotifFinder.PwmScorePValue(pwm, score), expected, 1e-10);
    }

    [Test]
    [Description("Non-tied thresholds: identical to TFM-Pvalue sc2pv (CRAN TFMPvalue src, converged ppv == pv).")]
    public void PwmScorePValue_Wikipedia_EqualsTfmPvalue()
    {
        var pwm = WikipediaPwm();
        // TFM-Pvalue sc2pv (granularity 0.1 → 1e-9) on the same log2 matrix: converged at g = 100 / 1000.
        AssertExact(MotifFinder.PwmScorePValue(pwm, 9.20753026510891), 0.000102996826171875);
        AssertExact(MotifFinder.PwmScorePValue(pwm, 5.853765132554455), 0.002780914306640625);
        AssertExact(MotifFinder.PwmScorePValue(pwm, 0.0), 0.06875991821289062);
        var tata = MotifFinder.BucherPromoterMatrices.TataBox.Pwm;
        AssertExact(MotifFinder.PwmScorePValue(tata, 13.282082872405402), 6.019137799739838e-06, 1e-10);
        AssertExact(MotifFinder.PwmScorePValue(tata, 7.891041436202701), 0.000597665086388588, 1e-10);
        // At a score attained by a word TFM-Pvalue does not converge (pv 0 / 0.006023406982421875 at g = 1e9 for
        // the Wikipedia sites at 2 and 6); enumeration gives 4^-9 and 0.006160736083984375 (tests above).
    }

    [Test]
    [Description("TFM-Pvalue pv2sc returns the same exact p-value but a rounded score (7.11, 9.27, 10.5073242, 14.32642).")]
    public void PwmScoreThresholdForPValue_PValuesEqualTfmPvalue()
    {
        Assert.That(MotifFinder.PwmScoreThresholdForPValue(WikipediaPwm(), 1e-3).PValue, Is.EqualTo(0.000980377197265625).Within(1e-18));
        Assert.That(MotifFinder.PwmScoreThresholdForPValue(WikipediaPwm(), 1e-4).PValue, Is.EqualTo(9.918212890625e-05).Within(1e-19));
        var tata = MotifFinder.BucherPromoterMatrices.TataBox.Pwm;
        Assert.That(MotifFinder.PwmScoreThresholdForPValue(tata, 1e-4).PValue, Is.EqualTo(9.999983012676239e-05).Within(1e-15));
        Assert.That(MotifFinder.PwmScoreThresholdForPValue(tata, 1e-6).PValue, Is.EqualTo(9.993091225624084e-07).Within(1e-17));
    }

    [Test]
    [Description("Biopython's ScoreDistribution is a grid approximation; the exact p-value of its FPR threshold differs.")]
    public void PwmScorePValue_ContrastsWithGridDistribution()
    {
        var pwm = WikipediaPwm();
        double gridThreshold = pwm.ScoreDistribution().ThresholdFpr(0.01); // Biopython 4.028388324862519
        var exact = MotifFinder.PwmScoreThresholdForPValue(pwm, 0.01);
        Assert.That(gridThreshold, Is.EqualTo(4.028388324862519).Within(1e-9));
        Assert.That(exact.Score, Is.EqualTo(4.028050165603465).Within(1e-12));
        Assert.That(MotifFinder.PwmScorePValue(pwm, gridThreshold).PValue, Is.LessThanOrEqualTo(0.01));
    }

    #endregion

    #region PwmScoreThresholdForPValue

    [TestCase(0.2, 3.0, 10.0 / 64)]
    [TestCase(0.6, 0.0, 37.0 / 64)]
    [TestCase(0.5, 3.0, 10.0 / 64)]
    [TestCase(1.0 / 64, 6.0, 1.0 / 64)]
    [TestCase(1.0, -3.0, 1.0)]
    public void PwmScoreThresholdForPValue_ToyMatrix(double p, double score, double pv)
    {
        var r = MotifFinder.PwmScoreThresholdForPValue(ToyPwm(), p);
        Assert.That(r.Score, Is.EqualTo(score));
        AssertExact(r, pv);
    }

    [Test]
    public void PwmScoreThresholdForPValue_BelowConsensusProbability_IsPositiveInfinity()
    {
        var r = MotifFinder.PwmScoreThresholdForPValue(ToyPwm(), 0.01);
        Assert.That(r.Score, Is.EqualTo(double.PositiveInfinity));
        AssertExact(r, 0.0);
        Assert.That(MotifFinder.PwmScoreThresholdForPValue(WikipediaPwm(), 1e-6).Score, Is.EqualTo(double.PositiveInfinity));
    }

    /// <summary>
    /// Expected: smallest word score t with P(S ≥ t) ≤ p by exhaustive enumeration; the next lower word
    /// score has P &gt; p (e.g. p = 0.01: next word score 4.028050165603464, one ulp below, P = 0.01003265380859375).
    /// </summary>
    [TestCase(1e-2, 4.028050165603465, 0.009918212890625)]
    [TestCase(1e-3, 7.150252343998389, 0.000980377197265625)]
    [TestCase(1e-4, 9.31606123696847, 9.918212890625e-05)]
    [TestCase(1e-5, 11.320507141999663, 7.62939453125e-06)]
    public void PwmScoreThresholdForPValue_Wikipedia_EqualsExhaustiveEnumeration(double p, double score, double pv)
    {
        var r = MotifFinder.PwmScoreThresholdForPValue(WikipediaPwm(), p);
        Assert.That(r.Score, Is.EqualTo(score).Within(1e-12));
        AssertExact(r, pv);
    }

    [TestCase(1e-2, 4.279588932599427, 0.009998928)]
    [TestCase(1e-4, 9.567600003964435, 9.5256e-05)]
    [TestCase(1e-5, 11.70753026510891, 5.832000000000001e-06)]
    public void PwmScoreThresholdForPValue_WikipediaBackground_EqualsExhaustiveEnumeration(double p, double score, double pv)
    {
        var r = MotifFinder.PwmScoreThresholdForPValue(WikipediaPwm(), p, NonUniform);
        Assert.That(r.Score, Is.EqualTo(score).Within(1e-12));
        AssertExact(r, pv, 1e-11);
    }

    [TestCase(1e-3, 7.041061417317205, 0.0009999992325901985)]
    [TestCase(1e-4, 10.50732497620838, 9.999983012676239e-05)]
    [TestCase(1e-6, 14.326493167908382, 9.993091225624084e-07)]
    public void PwmScoreThresholdForPValue_BucherTataBox_EqualsExhaustiveEnumeration(double p, double score, double pv)
    {
        var r = MotifFinder.PwmScoreThresholdForPValue(MotifFinder.BucherPromoterMatrices.TataBox.Pwm, p);
        Assert.That(r.Score, Is.EqualTo(score).Within(1e-12));
        AssertExact(r, pv, 1e-10);
    }

    /// <summary>8-column random log-odds matrix (seed 11, case 14 of the TFM-Pvalue cross-check) and its background.</summary>
    private static readonly double[] Random8Background = { 0.2795086578305266, 0.05418444709870725, 0.4216071278300304, 0.24469976724073575 };

    private static PositionWeightMatrix Random8Pwm() => new(new[,]
    {
        { -0.7360795913742273, -1.8867899131704242, -1.5235349559941267, -0.7710183582934051, 0.26842791564649543, -0.6953013042606067, -1.088815090853583, 0.036974501002490306 },
        { 2.1246772824688533, 2.71044996743149, -0.3789851161504548, 1.4909543432851726, 2.8475234136582253, 2.7443972993548273, 1.6484956202716525, 2.3537310050143936 },
        { -0.28558423910334624, -3.0103109162349093, 0.36150605444481204, -0.32654993924044745, -1.3680800093557965, -0.21555124819099777, -0.2191944194221844, -0.556031785363215 },
        { 0.03091536576591566, 1.0803840420302158, 0.38455232038061604, 0.5068351071917819, -1.1535092053715117, -1.2739392157625051, 0.4597586645697902, -0.7107726386426937 },
    }, 8);

    [Test]
    [Description("TFM-Pvalue pv2sc returns 8.6866 / P = 9.923731321091039e-05; the largest achievable p-value ≤ 1e-4 is " +
                 "9.930678854038951e-05 at 8.685556713812995 (next lower word score 8.685479274754845 has P = 1.00020651560101e-4) — exhaustive enumeration.")]
    public void PwmScoreThresholdForPValue_LargestAchievablePValue_WhereTfmIsNotMaximal()
    {
        var r = MotifFinder.PwmScoreThresholdForPValue(Random8Pwm(), 1e-4, Random8Background);
        Assert.That(r.Score, Is.EqualTo(8.685556713812995).Within(1e-12));
        AssertExact(r, 9.930678854038951e-05, 1e-10);
        Assert.That(r.PValue, Is.GreaterThan(9.923731321091039e-05));
    }

    [Test]
    [Description("20-column random matrix at score 0: the band of undecided words exceeds the budget at every scale, " +
                 "so certified bounds are returned (conservative upper bound as the p-value).")]
    public void PwmScorePValue_BudgetExhausted_ReturnsCertifiedBounds()
    {
        var pwm = new PositionWeightMatrix(new[,]
        {
            { 0.3125902303080262, -0.04222823538929744, 0.06342538014163379, -6.022367813028454, 1.036220187912753, 1.2352164616940315, -0.02345897282398891, 0.5849625007211562, -5.426264754702098, -1.084888897586513, 0.5555187228286743, 0.9857861407802992, 0.874469117916141, -5.087462841250339, 0.5555187228286743, 1.0677446066358343, 0.5613112326565108, 0.6814704815745026, 0.44312980630819676, 0.2606517545227997 },
            { -0.9098021910284219, -0.22470628717469449, -0.13588342808177303, 0.6918777046376682, -1.024662054234269, -1.1926450779423958, 0.8139880143900512, -1.8231222379159209, 1.1128940564059333, -1.7369655941662063, -0.4296842752432452, -2.502500340529183, 0.4405725913859814, 1.6807214835265871, -1.4005379295837288, -1.308122295362332, -1.713695814843359, 0.6119295483214257, 0.44312980630819676, 0.2606517545227997 },
            { 0.5849625007211562, 0.21842351913350247, 0.39463128861700686, -0.18947779886371285, 0.45720695352278334, -1.8845227825800641, -1.0962153152593033, 0.6662626028230043, -0.7824085649273732, 0.7410817026384381, -1.6520766965796931, -1.5849625007211563, -1.304854581528421, -5.087462841250339, 0.874469117916141, -0.49098635251214234, -0.24902754783991468, -1.0840642647884744, 0.0, 0.16196747966339312 },
            { -0.46234321405720047, 0.013805799525030428, -0.4533656179379432, 0.5775450291586732, -3.560714954474479, -0.09310940439148147, -0.33948646627166706, -0.6655809609294407, 0.3016556998611012, 0.63890130783196, 0.5555187228286743, 0.5943611987234058, -1.304854581528421, -0.4436066514756146, -1.6520766965796931, -0.3428877135230085, 0.4626269577971039, -1.3356030317844387, -1.831877241191673, -1.0435016386365865 },
        }, 20);

        var r = MotifFinder.PwmScorePValue(pwm, 0.0);
        Assert.That(r.IsExact, Is.False);
        Assert.That(r.PValueLowerBound, Is.LessThan(r.PValueUpperBound));
        Assert.That(r.PValue, Is.EqualTo(r.PValueUpperBound));
        Assert.That(r.PValueUpperBound - r.PValueLowerBound, Is.LessThan(1e-4 * r.PValue));
        Assert.That(r.PValue, Is.EqualTo(0.05840680768687889).Within(1e-12));
        Assert.That(r.PValueLowerBound, Is.EqualTo(0.058402542608746444).Within(1e-12));
    }

    [Test]
    public void PwmScoreThresholdForPValue_RoundTripsThroughScorePValue()
    {
        var pwm = WikipediaPwm();
        foreach (double p in new[] { 0.3, 0.05, 1e-2, 1e-3, 1e-4, 1e-5 })
        {
            var t = MotifFinder.PwmScoreThresholdForPValue(pwm, p, NonUniform);
            var back = MotifFinder.PwmScorePValue(pwm, t.Score, NonUniform);
            Assert.That(back.PValue, Is.EqualTo(t.PValue).Within(1e-12 * t.PValue), $"p = {p}"); // same word set; summation order differs
            Assert.That(t.PValue, Is.LessThanOrEqualTo(p));
        }
    }

    #endregion

    #region Validation

    [Test]
    public void PwmPValue_InvalidArguments_Throw()
    {
        var pwm = WikipediaPwm();
        Assert.Throws<ArgumentNullException>(() => MotifFinder.PwmScorePValue(null!, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmScorePValue(pwm, double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmScorePValue(pwm, double.PositiveInfinity));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmScoreThresholdForPValue(pwm, -0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmScoreThresholdForPValue(pwm, 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmScoreThresholdForPValue(pwm, double.NaN));
        Assert.Throws<ArgumentException>(() => MotifFinder.PwmScorePValue(pwm, 0, new[] { 0.5, 0.5 }));

        var infinite = MotifFinder.CreatePwm(new[] { "ACGT" }, 0.0);
        Assert.Throws<ArgumentException>(() => MotifFinder.PwmScorePValue(infinite, 0));
        Assert.Throws<ArgumentException>(() => MotifFinder.PwmScoreThresholdForPValue(infinite, 0.01));
    }

    [Test]
    public void PwmPValue_EmptyMatrix_ScoresZero()
    {
        var empty = new PositionWeightMatrix(new double[4, 0], 0);
        AssertExact(MotifFinder.PwmScorePValue(empty, 0.0), 1.0);
        AssertExact(MotifFinder.PwmScorePValue(empty, 0.5), 0.0);
        Assert.That(MotifFinder.PwmScoreThresholdForPValue(empty, 1.0).Score, Is.EqualTo(0.0));
    }

    #endregion

    #region CreatePwm per-base / JASPAR pseudocounts (Biopython 1.88)

    private static readonly string[] BioSeqs = { "TACAA", "TACGC", "TACAC", "TACCC", "AACCC", "AATGC", "AATGC" };

    private static void AssertMatrix(PositionWeightMatrix pwm, double[,] expected)
    {
        Assert.That(pwm.Length, Is.EqualTo(expected.GetLength(1)));
        for (int b = 0; b < 4; b++)
            for (int j = 0; j < pwm.Length; j++)
                Assert.That(pwm.Matrix[b, j], Is.EqualTo(expected[b, j]).Within(1e-12), $"[{b},{j}]");
    }

    [Test]
    [Description("motifs.create(seqs).counts.normalize(pseudocounts={A:.1,C:.4,G:.2,T:.3}).log_odds({A:.3,C:.2,G:.2,T:.3})")]
    public void CreatePwm_PerBasePseudocounts_EqualsBiopython()
    {
        var pwm = MotifFinder.CreatePwm(BioSeqs, new[] { 0.1, 0.4, 0.2, 0.3 }, NonUniform);
        AssertMatrix(pwm, new[,]
        {
            { 0.3692338096657191, 1.5647846187835261, -4.584962500721156, -0.19264507794239571, -1.1255308820838588 },
            { -2.0, -2.0, 1.7548875021634687, 0.584962500721156, 2.0 },
            { -3.0, -3.0, -3.0, 1.0, -3.0 },
            { 0.8413022539809418, -3.0, -0.06140054466414342, -3.0, -3.0 },
        });
    }

    [Test]
    [Description("m.background = {A:.3,C:.2,G:.2,T:.3}; normalize(pseudocounts=jaspar.calculate_pseudocounts(m)).log_odds(bg); pseudocounts √7·q.")]
    public void CreatePwmWithJasparPseudocounts_Background_EqualsBiopython()
    {
        var pwm = MotifFinder.CreatePwmWithJasparPseudocounts(BioSeqs, NonUniform);
        AssertMatrix(pwm, new[,]
        {
            { 0.39068723333471306, 1.4293850794872591, -1.866216153661624, -0.05073780119647756, -0.6899689797767886 },
            { -1.866216153661624, -1.8662161536616237, 1.5190922596283898, 0.39068723333471306, 1.758929724302279 },
            { -1.866216153661624, -1.8662161536616237, -1.866216153661624, 0.8713553378958487, -1.866216153661624 },
            { 0.728219246609793, -1.8662161536616237, -0.05073780119647772, -1.8662161536616237, -1.866216153661624 },
        });
        double[] pc = MotifFinder.JasparPseudocounts(new double[,] { { 4 }, { 0 }, { 0 }, { 3 } }, NonUniform);
        Assert.That(pc, Is.EqualTo(new[] { 0.7937253933193772, 0.5291502622129182, 0.5291502622129182, 0.7937253933193772 }).Within(1e-15));
    }

    [Test]
    public void CreatePwmWithJasparPseudocounts_Uniform_EqualsBiopython()
    {
        var pwm = MotifFinder.CreatePwmWithJasparPseudocounts(BioSeqs);
        AssertMatrix(pwm, new[,]
        {
            { 0.6025166839942642, 1.6677215545261992, -1.8662161536616237, 0.14231225004334136, -0.5374613073671782 },
            { -1.8662161536616237, -1.866216153661624, 1.2312748842246843, 0.14231225004334136, 1.465939992501705 },
            { -1.8662161536616237, -1.866216153661624, -1.8662161536616237, 0.602516683994264, -1.866216153661624 },
            { 0.950881410368785, -1.866216153661624, 0.14231225004334164, -1.866216153661624, -1.866216153661624 },
        });
    }

    [Test]
    public void CreatePwm_PerBase_EqualsScalarAndFromCounts()
    {
        var scalar = MotifFinder.CreatePwm(WikipediaSequences, 0.25, NonUniform);
        var perBase = MotifFinder.CreatePwm(WikipediaSequences, new[] { 0.25, 0.25, 0.25, 0.25 }, NonUniform);
        for (int b = 0; b < 4; b++)
            for (int j = 0; j < 9; j++)
                Assert.That(perBase.Matrix[b, j], Is.EqualTo(scalar.Matrix[b, j]).Within(1e-15));
    }

    [Test]
    public void CreatePwm_PerBase_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => MotifFinder.CreatePwm(null!, new[] { 1.0, 1, 1, 1 }));
        Assert.Throws<ArgumentNullException>(() => MotifFinder.CreatePwm(BioSeqs, (IReadOnlyList<double>)null!));
        Assert.Throws<ArgumentException>(() => MotifFinder.CreatePwm(BioSeqs, new[] { 1.0, 1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.CreatePwm(BioSeqs, new[] { 1.0, -1, 1, 1 }));
        Assert.Throws<ArgumentException>(() => MotifFinder.CreatePwm(Array.Empty<string>(), new[] { 1.0, 1, 1, 1 }));
        Assert.Throws<ArgumentException>(() => MotifFinder.CreatePwmWithJasparPseudocounts(new[] { "ACGT", "AC" }));
        Assert.Throws<ArgumentException>(() => MotifFinder.CreatePwmWithJasparPseudocounts(Array.Empty<string>()));
    }

    #endregion
}
