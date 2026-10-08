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

    [Test]
    public void WeightedEditAndDamerau_RandomInputs_OnlyDocumentedExceptions()
    {
        var rng = new Random(350_101);
        const string pool = "ACGTacgtNn-αβ中😀 ";
        int[] costPool = { -1, 0, 1, 2, 3, 7, 1000, int.MaxValue / 3, int.MaxValue };
        for (int trial = 0; trial < 1500; trial++)
        {
            string a = rng.Next(20) == 0 ? null! : RandomString(rng, rng.Next(0, 40), pool);
            string b = RandomString(rng, rng.Next(0, 40), pool);
            int I() => costPool[rng.Next(costPool.Length)];
            var e = new EditCosts(I(), I(), I());
            var d = new DamerauCosts(I(), I(), I(), I());

            Documented(() =>
            {
                int v = ApproximateMatcher.EditDistance(a, b, e);
                Assert.That(v, Is.GreaterThanOrEqualTo(0));
            }, allowOverflow: true);
            Documented(() =>
            {
                var al = ApproximateMatcher.GetEditAlignment(a, b, e);
                Assert.That(al.Operations.Count(c => c != 'D'), Is.EqualTo(a.Length));
            }, allowOverflow: true);
            Documented(() => ApproximateMatcher.GetEditAlignmentLinearSpace(a, b, e), allowOverflow: true);
            Documented(() => ApproximateMatcher.FindEditEndPositions(b, a ?? "", rng.Next(-1, 10), e).ToList(), allowOverflow: true);
            Documented(() => ApproximateMatcher.OptimalStringAlignmentDistance(a, b, d), allowOverflow: true);
            Documented(() => ApproximateMatcher.DamerauLevenshteinDistance(a, b, d), allowOverflow: true);
            Documented(() => ApproximateMatcher.GetOptimalStringAlignment(a, b, d), allowOverflow: true);
            Documented(() => ApproximateMatcher.GetDamerauLevenshteinAlignment(a, b, d), allowOverflow: true);
        }
    }

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
                Documented(() => CheckResult(MotifFinder.PwmMarkovScorePValue(pwm, score, OligoBackgroundModel.Equiprobable, options)));
                Documented(() => CheckResult(MotifFinder.PwmMarkovScorePValue(pwm, score, OligoBackgroundModel.MarkovFromInput(1), options)));
            }
            Assert.That(double.IsNaN(pwm.MaxScore) || pwm.MaxScore >= pwm.MinScore || double.IsNaN(pwm.MinScore), Is.True);
        }
    }

    [Test]
    public void AlphabetPwm_RandomAlphabetsAndSequences_OnlyDocumentedExceptions()
    {
        var rng = new Random(350_103);
        const string pool = "ACDEFGHIKLMNPQRSTVWYacgtX-*?ä";
        for (int trial = 0; trial < 400; trial++)
        {
            string alphabet = RandomString(rng, rng.Next(0, 8), pool);
            int length = rng.Next(0, 8);
            var seqs = Enumerable.Range(0, rng.Next(0, 6)).Select(_ => RandomString(rng, rng.Next(2) == 0 ? length : rng.Next(0, 8), pool)).ToList();
            double pseudo = rng.Next(6) switch { 0 => -1, 1 => double.NaN, 2 => 0, _ => rng.NextDouble() };
            Documented(() =>
            {
                var pwm = MotifFinder.CreateAlphabetPwm(seqs, alphabet, pseudo, null, rng.Next(2) == 0);
                string text = RandomString(rng, rng.Next(0, 30), pool);
                var scores = MotifFinder.CalculateAlphabetPwmScores(text, pwm);
                Assert.That(scores.Length, Is.EqualTo(Math.Max(0, text.Length - pwm.Length + 1)));
                _ = MotifFinder.ScanWithAlphabetPwm(text, pwm, pwm.MinScore).ToList();
                if (pwm.Length <= 6)
                {
                    _ = MotifFinder.AlphabetPwmScorePValue(pwm, (pwm.MaxScore + pwm.MinScore) / 2);
                    _ = MotifFinder.AlphabetPwmScoreThresholdForPValue(pwm, 0.01);
                }
            });
        }
    }

    [Test]
    public void DecipherAndEmbossAuto_RandomInput_OnlyDocumentedExceptions()
    {
        var rng = new Random(350_104);
        const string pool = "ACGTUNRYKMSWBDHVacgtu-.+~*?XEFILPQZJO 1é";
        for (int trial = 0; trial < 1000; trial++)
        {
            var rows = Enumerable.Range(0, rng.Next(0, 8)).Select(_ => RandomString(rng, rng.Next(0, 20), pool)).ToArray();
            int maxLen = rows.Length == 0 ? 0 : rows.Max(r => r.Length);
            Documented(() =>
            {
                string c = MotifFinder.GenerateDecipherConsensus(rows, (DecipherSequenceType)rng.Next(0, 4), rng.NextDouble() * 1.2 - 0.1,
                    rng.Next(2) == 0, "+-NX#"[rng.Next(5)], rng.Next(3) == 0 ? null : rng.NextDouble() * 1.2, rng.Next(2) == 0, rng.Next(2) == 0);
                Assert.That(c.Length, Is.EqualTo(maxLen));
            });
            Documented(() =>
            {
                string c = MotifFinder.GenerateEmbossConsensus(rows, true, ConsensusResidueType.Auto, rng.Next(3) == 0 ? null : rng.Next(-1, 5));
                Assert.That(c.Length, Is.EqualTo(maxLen));
            });
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
                var r = MotifFinder.AnalyzeOligoStrings(seqs, k, options);
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
                var r = MotifFinder.AnalyzeDyads(dna, dyad);
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

    private static int Hamming(string a, string b) => a.Zip(b).Count(p => p.First != p.Second);

    [Test]
    public void Sigma70Pairing_EqualsBruteForce_AndReverseComplementRelation()
    {
        var rng = new Random(350_106);
        const string c35 = MotifFinder.Sigma70Minus35Consensus, c10 = MotifFinder.Sigma70Minus10Consensus;
        for (int trial = 0; trial < 120; trial++)
        {
            // Plant consensus-like boxes in random DNA so that candidates exist.
            var chars = RandomString(rng, rng.Next(0, 160), "ACGT").ToCharArray();
            for (int plant = 0; plant < 3 && chars.Length >= 40; plant++)
            {
                int at = rng.Next(0, chars.Length - 40);
                int sp = rng.Next(14, 23);
                for (int i = 0; i < 6; i++) { chars[at + i] = c35[i]; chars[at + 6 + sp + i] = c10[i]; }
                chars[at + rng.Next(6)] = "ACGT"[rng.Next(4)];
            }
            string seq = new(chars);
            int max35 = rng.Next(0, 4), max10 = rng.Next(0, 4), minSp = rng.Next(12, 18), maxSp = minSp + rng.Next(0, 6);

            var expected = new List<(int, int, int)>();
            for (int i = 0; i + 6 <= seq.Length; i++)
                for (int sp = minSp; sp <= maxSp; sp++)
                {
                    int j = i + 6 + sp;
                    if (j + 6 > seq.Length) break;
                    if (Hamming(seq.Substring(i, 6), c35) <= max35 && Hamming(seq.Substring(j, 6), c10) <= max10)
                        expected.Add((i, j, sp));
                }
            var plus = MotifFinder.FindSigma70Promoters(new DnaSequence(seq), max35, max10, minSp, maxSp).ToList();
            string msg = $"trial {trial} {seq}";
            Assert.That(plus.Select(c => (c.Minus35Start, c.Minus10Start, c.Spacer)), Is.EquivalentTo(expected), msg);
            foreach (var c in plus)
            {
                Assert.That(c.Mismatches35, Is.EqualTo(Hamming(c.Minus35, c35)), msg);
                Assert.That(c.Mismatches10, Is.EqualTo(Hamming(c.Minus10, c10)), msg);
                Assert.That(c.TotalMismatches, Is.EqualTo(c.Mismatches35 + c.Mismatches10), msg);
                Assert.That(c.SpacerDeviation, Is.EqualTo(Math.Abs(c.Spacer - 17)), msg);
            }

            string rc = DnaSequence.GetReverseComplementString(seq);
            var both = MotifFinder.FindSigma70Promoters(new DnaSequence(seq), max35, max10, minSp, maxSp, bothStrands: true).ToList();
            var rcPlus = MotifFinder.FindSigma70Promoters(new DnaSequence(rc), max35, max10, minSp, maxSp).ToList();
            Assert.That(both.Where(c => c.Strand == '+'), Is.EqualTo(plus), msg);
            Assert.That(both.Where(c => c.Strand == '-').Select(c => (c.Minus35, c.Minus10, c.Spacer, seq.Length - c.Minus35Start - 6, seq.Length - c.Minus10Start - 6)),
                Is.EquivalentTo(rcPlus.Select(c => (c.Minus35, c.Minus10, c.Spacer, c.Minus35Start, c.Minus10Start))), msg);

            Documented(() =>
            {
                foreach (var p in MotifFinder.PredictSigma70Promoters(new DnaSequence(seq), rng.Next(2) == 0, rng.Next(2) == 0))
                    Assert.That(double.IsFinite(p.DeltaGTotal), Is.True, msg);
            });
        }
    }
}
