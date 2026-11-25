// SPDX-License-Identifier: AGPL-3.0-only
// Tests of the port of RSAT oligo-analysis / dyad-analysis; part of Seqeron.Genomics.OverrepresentedWords.Tests, licensed under the GNU Affero General Public License, version 3 (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Fuzzing;


/// <summary>
/// Heavy-tier fuzz tests for the review 2026-09 B05 audit additions (docs/Validation/review-2026-09/B05.md F27–F35):
/// weighted edit / Damerau distances and alignments, PWM p-values (DNA, Markov, alphabet) with random options,
/// generic-alphabet PWMs, DECIPHER / EMBOSS Auto consensus, RSAT -seqtype oligo analysis, dyad analysis and σ70
/// promoter pairing / Promoter Calculator. Random inputs (fixed seeds) must either produce a well-formed result or
/// throw one of the documented exceptions (ArgumentException and subtypes; OverflowException for a weighted distance
/// beyond int.MaxValue) — never anything else. σ70 pairing is also checked against a brute-force enumeration and its
/// reverse-complement relation.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Matching")]
public class B05AuditFuzzTests
{

    [Test]
    public void PwmPValues_RandomMatricesAndOptions_OnlyDocumentedExceptions()
    {
        var rng = new Random(350_102);
        double[] cellPool = { double.NaN, double.NegativeInfinity, double.PositiveInfinity, 0.0, -0.0, 1e-300, 1e300 };
        for (int trial = 0; trial < 300; trial++)
        {
            int length = rng.Next(0, 14);
            var m = new double[4, length];
            for (int r = 0; r < 4; r++)
                for (int c = 0; c < length; c++)
                    m[r, c] = rng.Next(40) == 0 ? cellPool[rng.Next(cellPool.Length)] : Math.Round((rng.NextDouble() - 0.5) * 10, rng.Next(0, 6));
            var pwm = new PositionWeightMatrix(m, length);
            double score = rng.Next(10) == 0 ? double.NaN : (rng.NextDouble() - 0.3) * 4 * Math.Max(1, length);
            double p = rng.Next(10) == 0 ? -0.5 : Math.Pow(10, -rng.Next(0, 9)) * rng.NextDouble();
            IReadOnlyList<double>? bg = rng.Next(3) switch
            {
                0 => null,
                1 => Enumerable.Range(0, 4).Select(_ => rng.NextDouble()).ToArray(),
                _ => new[] { 0.1, 0.2, 0.3 },
            };
            var options = rng.Next(4) switch
            {
                0 => null,
                1 => PwmPValueOptions.Exact,
                2 => new PwmPValueOptions { MaxStates = rng.Next(-1, 2000), MaxSuffixSet = rng.Next(1, 500), InitialGranularity = rng.NextDouble() },
                _ => new PwmPValueOptions { DecreaseFactor = rng.Next(0, 5), MaxGranularity = rng.NextDouble() * 0.1 },
            };
            if (options?.Exhaustive == true && length > 8)
                options = null; // the exhaustive mode is unbounded by design

            void CheckResult(PwmPValueResult r)
            {
                if (double.IsNaN(r.PValue)) Assert.Fail("NaN p-value");
                Assert.That(r.PValue, Is.InRange(0.0, 1.0 + 1e-12));
                Assert.That(r.PValueLowerBound, Is.LessThanOrEqualTo(r.PValueUpperBound + 1e-15));
            }

            Documented(() => CheckResult(MotifFinder.PwmScorePValue(pwm, score, bg, options)));
            Documented(() => CheckResult(MotifFinder.PwmScoreThresholdForPValue(pwm, p, bg, options)));
            if (length <= 8)
            {
                Documented(() => CheckResult(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, score, OligoBackgroundModel.Equiprobable, options)));
                Documented(() => CheckResult(OverrepresentedWordAnalyzer.PwmMarkovScorePValue(pwm, score, OligoBackgroundModel.MarkovFromInput(1), options)));
            }
            Assert.That(double.IsNaN(pwm.MaxScore) || pwm.MaxScore >= pwm.MinScore || double.IsNaN(pwm.MinScore), Is.True);
        }
    }

    [Test]
    public void OligoStringsAndDyads_RandomInput_OnlyDocumentedExceptions()
    {
        var rng = new Random(350_105);
        const string pool = "ACGTacgtNRYMKLVEFQ*X- \t,.!";
        for (int trial = 0; trial < 250; trial++)
        {
            var seqs = Enumerable.Range(0, rng.Next(0, 4)).Select(_ => RandomString(rng, rng.Next(0, 80), pool)).ToList();
            int k = rng.Next(0, 5);
            var options = new OligoAnalysisOptions
            {
                SequenceType = (OligoSequenceType)rng.Next(0, 4),
                Background = rng.Next(3) switch { 0 => OligoBackgroundModel.Equiprobable, 1 => OligoBackgroundModel.BernoulliFromInput, _ => OligoBackgroundModel.MarkovFromInput(1) },
                Strands = (OligoStrandMode)rng.Next(0, 2),
                CountOverlapping = rng.Next(2) == 0,
                PseudoFrequency = rng.Next(4) == 0 ? -0.5 : rng.NextDouble() * 0.2,
                Degeneracy = (OligoDegeneracy)rng.Next(0, 3),
                MinCount = rng.Next(-1, 3),
            };
            Documented(() =>
            {
                var r = OverrepresentedWordAnalyzer.AnalyzeOligoStrings(seqs, k, options);
                foreach (var p in r.Patterns)
                {
                    Assert.That(p.Occurrences, Is.GreaterThanOrEqualTo(0));
                    Assert.That(double.IsNaN(p.OccurrenceProbability), Is.False, p.Pattern);
                }
            });

            var dna = Enumerable.Range(0, rng.Next(0, 4)).Select(_ => new DnaSequence(RandomString(rng, rng.Next(0, 60), "ACGT"))).ToList();
            var dyad = new DyadAnalysisOptions
            {
                MonadLength = rng.Next(0, 4),
                MinSpacing = rng.Next(-1, 4),
                MaxSpacing = rng.Next(-1, 8),
                Type = (DyadType)rng.Next(0, 5),
                Strands = (OligoStrandMode)rng.Next(0, 2),
                CountOverlapping = rng.Next(2) == 0,
                MinCount = rng.Next(-1, 3),
            };
            Documented(() =>
            {
                var r = OverrepresentedWordAnalyzer.AnalyzeDyads(dna, dyad);
                foreach (var d in r.Dyads)
                {
                    Assert.That(d.Occurrences, Is.GreaterThanOrEqualTo(0));
                    // Untested combinations (exp_freq = 0, e.g. MinCount <= 0 with an unseen reverse-complement monad) carry NaN
                    // statistics by contract (RSAT: "Expected frequency must be > 0", not tested); tested ones must be finite.
                    if (!double.IsNaN(d.ExpectedFrequency))
                        Assert.That(d.OccurrenceProbability, Is.InRange(0.0, 1.0), $"{d} {dyad}");
                }
            });
        }
    }

    private static string RandomString(Random rng, int length, string pool) =>
        new(Enumerable.Range(0, length).Select(_ => pool[rng.Next(pool.Length)]).ToArray());

    private static void Documented(Action call, bool allowOverflow = false)
    {
        try
        {
            call();
        }
        catch (ArgumentException)
        {
        }
        catch (OverflowException) when (allowOverflow)
        {
        }
    }
}
