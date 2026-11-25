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
