// PwmBackgroundChain.FromWordFrequencies + PwmMarkovScorePValue / PwmMarkovScoreThresholdForPValue (chain overloads): the Markov
// background of the exact PWM p-value from a table of (m+1)-mer frequencies. Oracles are hand-computed from the definition
// P(b | x) = (1 − ψ)·f(xb)/S(x) + ψ/4, P(x) = (1 − ψ)·S(x)/F + ψ/4^m and a brute-force enumeration of all 4^L words.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_PwmBackgroundChain_Tests
{
    private const double Tol = 1e-12;

    [Test]
    public void FromWordFrequencies_Order1_PseudoFrequencyZero_GivesRelativeFrequencies()
    {
        // f(AA) = f(AC) = 1: S(A) = 2, F = 2, every other prefix empty.
        var chain = PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AA"] = 1, ["AC"] = 1 }, 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(chain.Order, Is.EqualTo(1));
            Assert.That(chain.K, Is.EqualTo(4));
            Assert.That(chain.PrefixProbabilities[1], Is.EqualTo(new[] { 1.0, 0.0, 0.0, 0.0 }).Within(Tol));
            Assert.That(chain.Transitions[..4], Is.EqualTo(new[] { 0.5, 0.5, 0.0, 0.0 }).Within(Tol));
            Assert.That(chain.Transitions[4..8], Is.EqualTo(new[] { 0.25, 0.25, 0.25, 0.25 }).Within(Tol), "an empty prefix row is uniform");
        });
    }

    [Test]
    public void FromWordFrequencies_Order1_PseudoFrequency_SmoothsPrefixAndTransitions()
    {
        const double psi = 0.2;
        var chain = PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AA"] = 1, ["AC"] = 1 }, psi);

        Assert.Multiple(() =>
        {
            Assert.That(chain.PrefixProbabilities[1][0], Is.EqualTo((1 - psi) * 2.0 / 2.0 + psi / 4).Within(Tol));
            Assert.That(chain.PrefixProbabilities[1][1], Is.EqualTo(psi / 4).Within(Tol));
            Assert.That(chain.Transitions[0], Is.EqualTo((1 - psi) * 0.5 + psi / 4).Within(Tol));
            Assert.That(chain.Transitions[2], Is.EqualTo(psi / 4).Within(Tol));
            Assert.That(chain.PrefixProbabilities[1].Sum(), Is.EqualTo(1.0).Within(Tol));
        });
    }

    [Test]
    public void FromWordFrequencies_Order0_IsResidueComposition()
    {
        var chain = PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["A"] = 3, ["T"] = 1 }, 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(chain.Order, Is.EqualTo(0));
            Assert.That(chain.Residues, Is.EqualTo(new[] { 0.75, 0.0, 0.0, 0.25 }).Within(Tol));
        });
    }

    [Test]
    public void FromWordFrequencies_StrandInsensitive_AddsReverseComplementAndHalves()
    {
        // AAC (f = 1) without GTT: GTT takes 1, then both halve to 0.5 (neither is a palindrome).
        var chain = PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AAC"] = 1 }, 0.0, strandInsensitive: true);

        int aa = 0, gt = (2 * 4) + 3;
        Assert.Multiple(() =>
        {
            Assert.That(chain.PrefixProbabilities[2][aa], Is.EqualTo(0.5).Within(Tol));
            Assert.That(chain.PrefixProbabilities[2][gt], Is.EqualTo(0.5).Within(Tol));
            Assert.That(chain.Transitions[(aa * 4) + 1], Is.EqualTo(1.0).Within(Tol));
            Assert.That(chain.Transitions[(gt * 4) + 3], Is.EqualTo(1.0).Within(Tol));
        });
    }

    [Test]
    public void FromWordFrequencies_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => PwmBackgroundChain.FromWordFrequencies(null!));
            Assert.Throws<ArgumentException>(() => PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double>()));
            Assert.Throws<ArgumentException>(() => PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AA"] = 1, ["AAC"] = 1 }), "mixed lengths");
            Assert.Throws<ArgumentException>(() => PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AN"] = 1 }), "non-ACGT");
            Assert.Throws<ArgumentException>(() => PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AA"] = 0 }), "zero total");
            Assert.Throws<ArgumentOutOfRangeException>(() => PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AA"] = -1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AA"] = double.NaN }));
            Assert.Throws<ArgumentOutOfRangeException>(() => PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AA"] = 1 }, 1.5));
            Assert.Throws<ArgumentException>(() => PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { [new string('A', 12)] = 1 }), "order above 10");
        });
    }

    [TestCase(0.0)]
    [TestCase(0.05)]
    [TestCase(0.5)]
    public void PwmMarkovScorePValue_MatchesBruteForceEnumeration(double psi)
    {
        var table = new Dictionary<string, double>();
        int i = 1;
        foreach (char a in "ACGT")
            foreach (char b in "ACGT")
                table[$"{a}{b}"] = (i++ * 7 % 11) + 1;
        var chain = PwmBackgroundChain.FromWordFrequencies(table, psi);
        var pwm = new PositionWeightMatrix(new double[,]
        {
            { 1.2, -0.5, 0.3, 0.9 },
            { -0.7, 1.1, 0.2, -0.4 },
            { 0.4, 0.6, -1.0, 0.1 },
            { -0.3, -0.2, 0.8, -0.6 },
        }, 4);

        // enumerate all 4^4 words with their chain probability and window score
        var words = new List<(double Score, double Prob)>();
        for (int code = 0; code < 256; code++)
        {
            int[] w = { code >> 6 & 3, code >> 4 & 3, code >> 2 & 3, code & 3 };
            double p = chain.PrefixProbabilities[1][w[0]];
            double score = pwm.Matrix[w[0], 0];
            for (int c = 1; c < 4; c++)
            {
                p *= chain.Transitions[(w[c - 1] * 4) + w[c]];
                score += pwm.Matrix[w[c], c];
            }

            words.Add((score, p));
        }

        foreach (double threshold in new[] { -1.0, 0.0, 1.0, 2.0, 3.0 })
        {
            double expected = words.Where(x => x.Score >= threshold).Sum(x => x.Prob);
            var r = MotifFinder.PwmMarkovScorePValue(pwm, threshold, chain, PwmPValueOptions.Exact);
            Assert.That(r.PValue, Is.EqualTo(expected).Within(1e-9), $"threshold {threshold}");
        }
    }

    [Test]
    public void PwmMarkovScore_NullAndInvalidArguments_Throw()
    {
        var pwm = MotifFinder.CreatePwm(new[] { "ACGT", "ACGA" });
        var chain = PwmBackgroundChain.FromWordFrequencies(new Dictionary<string, double> { ["AC"] = 1 });
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MotifFinder.PwmMarkovScorePValue(null!, 0, chain));
            Assert.Throws<ArgumentNullException>(() => MotifFinder.PwmMarkovScorePValue(pwm, 0, (PwmBackgroundChain)null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmMarkovScorePValue(pwm, double.NaN, chain));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.PwmMarkovScoreThresholdForPValue(pwm, 2.0, chain));
            Assert.Throws<ArgumentException>(() => MotifFinder.PwmMarkovScorePValue(pwm, 0, new PwmBackgroundChain(5, new[] { 0.2, 0.2, 0.2, 0.2, 0.2 })));
        });
    }
}
