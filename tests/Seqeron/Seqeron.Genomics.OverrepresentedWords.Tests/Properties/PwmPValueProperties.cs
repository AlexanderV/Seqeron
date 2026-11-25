// SPDX-License-Identifier: AGPL-3.0-only
// Tests of the port of RSAT oligo-analysis / dyad-analysis; part of Seqeron.Genomics.OverrepresentedWords.Tests, licensed under the GNU Affero General Public License, version 3 (see LICENSE in the production project).
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
            AssertClose(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, alpha, OligoBackgroundModel.Bernoulli(q)).PValue,
                MotifFinder.PwmScorePValue(pwm, alpha, q).PValue, msg + " bernoulli");
            AssertClose(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, alpha, OligoBackgroundModel.Equiprobable).PValue,
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
                var r = OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, a, model);
                AssertClose(r.PValue, Tail(words, a), $"{msg} markov a={a:R} psi={psi}");
            }
            var t = OverrepresentedWordAnalyzer.PwmMarkovScoreThresholdForPValue(pwm, 0.05, model);
            if (!double.IsPositiveInfinity(t.Score))
            {
                Assert.That(t.PValue, Is.LessThanOrEqualTo(0.05 * (1 + 1e-12)), msg);
                AssertClose(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, t.Score, model).PValue, t.PValue, msg + " markov round trip");
            }
        }
    }

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
}
