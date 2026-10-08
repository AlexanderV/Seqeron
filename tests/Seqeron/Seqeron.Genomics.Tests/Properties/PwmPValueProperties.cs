namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the review 2026-09 B05 audit additions on PAT-PWM-001
/// (docs/Validation/review-2026-09/B05.md F29 exact PWM p-values (Touzet &amp; Varré 2007), F30 generic-alphabet PWM,
/// F32 p-value options / Markov backgrounds / alphabet p-values, F35 Max/MinScore with −∞ columns).
///
/// Oracle: exhaustive enumeration of all K^L words for small L. A word's score is taken from
/// <c>CalculatePwmScores</c> / <c>CalculateAlphabetPwmScores</c> (the same summation the engine is defined on), its
/// probability is Π q(w_i) (i.i.d.) or the documented RSAT segment_proba of the Markov model
/// (P(x) = (1 − ψ)·S(x)/F + ψ/4^m, P(b | x) = (1 − ψ)·f(xb)/S(x) + ψ/4).
///   P1 (exact)        PwmScorePValue(α) == Σ_{w: S(w) ≥ α} P(w), IsExact, bounds bracket the value.
///   P2 (monotone)     p(α) is non-increasing in α; p(MinScore) = 1; p(&gt; MaxScore) = 0.
///   P3 (inverse)      PwmScoreThresholdForPValue(p).PValue == max{ P(S ≥ s) ≤ p } over word scores s, ≤ p,
///                      and PwmScorePValue(threshold.Score) gives the same p-value.
///   P4 (options)      PwmPValueOptions.Exact == Default when the default result is exact.
///   P5 (Markov)       Bernoulli(q) / Equiprobable Markov p-values == the i.i.d. engine; an order-1 table == the
///                      brute-force Markov word probabilities.
///   P6 (alphabet)     CreateAlphabetPwm over "ACGT" == CreatePwm (matrix, scores, Max/Min, p-values, thresholds);
///                      alphabet p-values for K = 3 == brute force.
///   P7 (extrema)      Max/MinScore == max/min word score by enumeration, including −∞ cells (F35).
/// Fixed seeds, L ≤ 6 (4^6 = 4096 words), deterministic.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Matching")]
public class PwmPValueProperties
{
    private const string Acgt = "ACGT";

    private static double[] RandomBackground(Random rng, int k)
    {
        var q = Enumerable.Range(0, k).Select(_ => 0.05 + rng.NextDouble()).ToArray();
        double s = q.Sum();
        return q.Select(x => x / s).ToArray();
    }

    private static PositionWeightMatrix RandomPwm(Random rng, int length)
    {
        var seqs = Enumerable.Range(0, rng.Next(1, 12))
            .Select(_ => new string(Enumerable.Range(0, length).Select(_ => Acgt[rng.Next(4)]).ToArray())).ToList();
        double pseudo = rng.Next(3) switch { 0 => 0.25, 1 => 0.5 + rng.NextDouble(), _ => 0.01 };
        return MotifFinder.CreatePwm(seqs, pseudo, RandomBackground(rng, 4));
    }

    private static IEnumerable<string> AllWords(string alphabet, int length)
    {
        if (length == 0) { yield return ""; yield break; }
        foreach (var prefix in AllWords(alphabet, length - 1))
            foreach (char c in alphabet)
                yield return prefix + c;
    }

    private static List<(double Score, double Prob)> Enumerate(PositionWeightMatrix pwm, Func<string, double> prob) =>
        AllWords(Acgt, pwm.Length).Select(w => (MotifFinder.CalculatePwmScores(w, pwm)[0], prob(w))).ToList();

    private static Func<string, double> Iid(IReadOnlyList<double> q, string alphabet) =>
        w => w.Aggregate(1.0, (p, c) => p * q[alphabet.IndexOf(c)]);

    private static double Tail(List<(double Score, double Prob)> words, double alpha) =>
        words.Where(x => x.Score >= alpha).Sum(x => x.Prob);

    private static void AssertClose(double actual, double expected, string msg) =>
        Assert.That(actual, Is.EqualTo(expected).Within(1e-12 + 1e-10 * Math.Abs(expected)), msg);

    private static IEnumerable<double> Probes(Random rng, List<(double Score, double Prob)> words)
    {
        var scores = words.Select(x => x.Score).Distinct().OrderBy(s => s).ToList();
        for (int i = 0; i < 8; i++)
            yield return scores[rng.Next(scores.Count)];
        yield return scores[0];
        yield return scores[^1];
        yield return Math.BitIncrement(scores[^1]);
        yield return scores[0] - 1;
        yield return scores[0] + rng.NextDouble() * (scores[^1] - scores[0]);
    }

    [Test]
    public void P1_P2_P4_ExactPValue_EqualsEnumeration_MonotoneInScore()
    {
        var rng = new Random(290_001);
        for (int trial = 0; trial < 60; trial++)
        {
            var pwm = RandomPwm(rng, rng.Next(1, 7));
            var q = trial % 3 == 0 ? new[] { 0.25, 0.25, 0.25, 0.25 } : RandomBackground(rng, 4);
            var words = Enumerate(pwm, Iid(q, Acgt));
            var probes = Probes(rng, words).OrderBy(x => x).ToList();
            double previous = double.PositiveInfinity;
            foreach (double alpha in probes)
            {
                var r = MotifFinder.PwmScorePValue(pwm, alpha, q);
                string msg = $"trial {trial} L={pwm.Length} alpha={alpha:R}";
                Assert.That(r.IsExact, Is.True, msg);
                AssertClose(r.PValue, Tail(words, alpha), msg);
                Assert.That(r.PValueLowerBound, Is.LessThanOrEqualTo(r.PValue), msg);
                Assert.That(r.PValueUpperBound, Is.GreaterThanOrEqualTo(r.PValue), msg);
                Assert.That(r.PValue, Is.LessThanOrEqualTo(previous * (1 + 1e-12)), msg);
                previous = r.PValue;
                Assert.That(MotifFinder.PwmScorePValue(pwm, alpha, q, PwmPValueOptions.Exact).PValue,
                    Is.EqualTo(r.PValue).Within(1e-12 + 1e-10 * r.PValue), msg);
            }
            AssertClose(MotifFinder.PwmScorePValue(pwm, pwm.MinScore, q).PValue, 1.0, $"trial {trial} min");
            Assert.That(MotifFinder.PwmScorePValue(pwm, pwm.MaxScore + 1e-6, q).PValue, Is.Zero, $"trial {trial} max");
        }
    }

    [Test]
    public void P3_ThresholdForPValue_IsInverseOfEnumeration()
    {
        var rng = new Random(290_003);
        for (int trial = 0; trial < 50; trial++)
        {
            var pwm = RandomPwm(rng, rng.Next(1, 7));
            var q = RandomBackground(rng, 4);
            var words = Enumerate(pwm, Iid(q, Acgt));
            var tails = words.Select(x => x.Score).Distinct().Select(s => (Score: s, P: Tail(words, s))).ToList();
            foreach (double p in new[] { 1.0, 0.5, 0.1, 0.02, 1e-3, 1e-4, rng.NextDouble(), Math.Pow(10, -rng.Next(1, 6)) })
            {
                string msg = $"trial {trial} L={pwm.Length} p={p:R}";
                var t = MotifFinder.PwmScoreThresholdForPValue(pwm, p, q);
                var admissible = tails.Where(x => x.P <= p * (1 + 1e-12)).ToList();
                if (admissible.Count == 0)
                {
                    Assert.That(t.Score, Is.EqualTo(double.PositiveInfinity), msg);
                    continue;
                }
                double best = admissible.Max(x => x.P);
                AssertClose(t.PValue, best, msg);
                Assert.That(t.PValue, Is.LessThanOrEqualTo(p * (1 + 1e-12)), msg);
                AssertClose(MotifFinder.PwmScorePValue(pwm, t.Score, q).PValue, t.PValue, msg + " round trip");
            }
        }
    }

    [Test]
    public void P5_MarkovOrder0_EqualsIid_AndOrder1Table_EqualsBruteForce()
    {
        var rng = new Random(290_005);
        for (int trial = 0; trial < 25; trial++)
        {
            var pwm = RandomPwm(rng, rng.Next(1, 6));
            var q = RandomBackground(rng, 4);
            double alpha = pwm.MinScore + rng.NextDouble() * (pwm.MaxScore - pwm.MinScore);
            string msg = $"trial {trial} alpha={alpha:R}";
            AssertClose(MotifFinder.PwmMarkovScorePValue(pwm, alpha, OligoBackgroundModel.Bernoulli(q)).PValue,
                MotifFinder.PwmScorePValue(pwm, alpha, q).PValue, msg + " bernoulli");
            AssertClose(MotifFinder.PwmMarkovScorePValue(pwm, alpha, OligoBackgroundModel.Equiprobable).PValue,
                MotifFinder.PwmScorePValue(pwm, alpha).PValue, msg + " equiprobable");

            // Order-1 table of dinucleotide frequencies (some zero) with pseudo-frequency ψ.
            var table = new Dictionary<string, double>();
            foreach (var w in AllWords(Acgt, 2))
                table[w] = rng.Next(4) == 0 ? 0 : rng.Next(1, 50);
            if (table.Values.Sum() == 0) table["AC"] = 1;
            double psi = rng.Next(3) switch { 0 => 0.0, 1 => 0.01, _ => rng.NextDouble() * 0.5 };
            double total = table.Values.Sum();
            double S(char x) => Acgt.Sum(b => table[$"{x}{b}"]);
            double Prefix(char x) => (1 - psi) * S(x) / total + psi / 4;
            double Trans(char x, char b) => S(x) == 0 ? 0.25 : (1 - psi) * table[$"{x}{b}"] / S(x) + psi / 4;
            double Prob(string w)
            {
                double p = Prefix(w[0]);
                for (int i = 1; i < w.Length; i++) p *= Trans(w[i - 1], w[i]);
                return p;
            }
            var model = OligoBackgroundModel.MarkovFromOligoFrequencies(table, psi);
            var words = Enumerate(pwm, Prob);
            AssertClose(words.Sum(x => x.Prob), 1.0, msg + " normalisation");
            foreach (double a in Probes(rng, words))
            {
                var r = MotifFinder.PwmMarkovScorePValue(pwm, a, model);
                AssertClose(r.PValue, Tail(words, a), $"{msg} markov a={a:R} psi={psi}");
            }
            var t = MotifFinder.PwmMarkovScoreThresholdForPValue(pwm, 0.05, model);
            if (!double.IsPositiveInfinity(t.Score))
            {
                Assert.That(t.PValue, Is.LessThanOrEqualTo(0.05 * (1 + 1e-12)), msg);
                AssertClose(MotifFinder.PwmMarkovScorePValue(pwm, t.Score, model).PValue, t.PValue, msg + " markov round trip");
            }
        }
    }

    [Test]
    public void P6_AlphabetPwmOverAcgt_EqualsDnaPwm()
    {
        var rng = new Random(290_006);
        for (int trial = 0; trial < 60; trial++)
        {
            int length = rng.Next(1, 9);
            var seqs = Enumerable.Range(0, rng.Next(1, 15))
                .Select(_ => new string(Enumerable.Range(0, length).Select(_ => Acgt[rng.Next(4)]).ToArray())).ToList();
            double pseudo = rng.Next(2) == 0 ? 0.25 : rng.NextDouble() * 2 + 0.01;
            var q = RandomBackground(rng, 4);
            var dna = MotifFinder.CreatePwm(seqs, pseudo, q);
            var alpha = MotifFinder.CreateAlphabetPwm(seqs, Acgt, pseudo, q);
            string msg = $"trial {trial}";
            Assert.That(alpha.GetMatrix(), Is.EqualTo(dna.Matrix), msg);
            Assert.That(alpha.MaxScore, Is.EqualTo(dna.MaxScore), msg);
            Assert.That(alpha.MinScore, Is.EqualTo(dna.MinScore), msg);
            Assert.That(alpha.Consensus, Is.EqualTo(dna.Consensus), msg);

            string text = new(Enumerable.Range(0, rng.Next(length, 60)).Select(_ => Acgt[rng.Next(4)]).ToArray());
            Assert.That(MotifFinder.CalculateAlphabetPwmScores(text, alpha), Is.EqualTo(MotifFinder.CalculatePwmScores(text, dna)), msg);

            if (length <= 6)
            {
                double score = dna.MinScore + rng.NextDouble() * (dna.MaxScore - dna.MinScore);
                var a = MotifFinder.AlphabetPwmScorePValue(alpha, score, q);
                var d = MotifFinder.PwmScorePValue(dna, score, q);
                AssertClose(a.PValue, d.PValue, msg + " p-value");
                var ta = MotifFinder.AlphabetPwmScoreThresholdForPValue(alpha, 0.01, q);
                var td = MotifFinder.PwmScoreThresholdForPValue(dna, 0.01, q);
                Assert.That(ta.Score, Is.EqualTo(td.Score).Within(1e-9), msg + " threshold");
                AssertClose(ta.PValue, td.PValue, msg + " threshold p");
            }
        }
    }

    [Test]
    public void P6b_ThreeLetterAlphabet_PValues_EqualEnumeration()
    {
        var rng = new Random(290_016);
        const string abc = "XYZ";
        for (int trial = 0; trial < 40; trial++)
        {
            int length = rng.Next(1, 8);
            var seqs = Enumerable.Range(0, rng.Next(1, 10))
                .Select(_ => new string(Enumerable.Range(0, length).Select(_ => abc[rng.Next(3)]).ToArray())).ToList();
            var q = RandomBackground(rng, 3);
            var pwm = MotifFinder.CreateAlphabetPwm(seqs, abc, 0.1 + rng.NextDouble(), q);
            var words = AllWords(abc, length)
                .Select(w => (Score: MotifFinder.CalculateAlphabetPwmScores(w, pwm)[0], Prob: Iid(q, abc)(w))).ToList();
            foreach (double a in Probes(rng, words))
            {
                var r = MotifFinder.AlphabetPwmScorePValue(pwm, a, q);
                Assert.That(r.IsExact, Is.True);
                AssertClose(r.PValue, Tail(words, a), $"trial {trial} a={a:R}");
            }
            Assert.That(pwm.MaxScore, Is.EqualTo(words.Max(x => x.Score)).Within(1e-12));
            Assert.That(pwm.MinScore, Is.EqualTo(words.Min(x => x.Score)).Within(1e-12));
        }
    }

    [Test]
    public void P7_MaxMinScore_EqualEnumeration_IncludingNegativeInfinityCells()
    {
        var rng = new Random(290_007);
        for (int trial = 0; trial < 80; trial++)
        {
            int length = rng.Next(1, 6);
            var m = new double[4, length];
            for (int r = 0; r < 4; r++)
                for (int c = 0; c < length; c++)
                    m[r, c] = rng.Next(5) == 0 ? double.NegativeInfinity : (rng.NextDouble() - 0.5) * 8;
            if (trial % 7 == 0)
                for (int r = 0; r < 4; r++) m[r, 0] = double.NegativeInfinity; // all −∞ column
            var pwm = new PositionWeightMatrix(m, length);
            var scores = AllWords(Acgt, length).Select(w => MotifFinder.CalculatePwmScores(w, pwm)[0]).ToList();
            string msg = $"trial {trial}";
            Assert.That(pwm.MaxScore, double.IsFinite(scores.Max()) ? Is.EqualTo(scores.Max()).Within(1e-12) : Is.EqualTo(scores.Max()), msg);
            Assert.That(pwm.MinScore, double.IsFinite(scores.Min()) ? Is.EqualTo(scores.Min()).Within(1e-12) : Is.EqualTo(scores.Min()), msg);
        }
    }
}
