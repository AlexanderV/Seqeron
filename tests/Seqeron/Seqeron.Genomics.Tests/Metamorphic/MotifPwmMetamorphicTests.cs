namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic tests for behaviour added by review 2026-09 batch B05
/// (docs/Validation/review-2026-09/B05.md).
///
/// Relations:
///   • EDIT-SELLERS (F2): the set of window end positions of
///     <see cref="ApproximateMatcher.FindWithEdits(string, string, int)"/> equals the Sellers end set of
///     <see cref="ApproximateMatcher.FindEditEndPositions(string, string, int)"/>, and the minimum window
///     distance per end equals the Sellers distance C[m, j]; every window's distance is its true
///     Levenshtein distance.
///   • CAVENER (F13c, Cavener 1987): consensus depends only on per-column counts → invariant under
///     row permutation and case; a single sequence is its own consensus; output ⊂ IUPAC; every
///     maximal-count base of a column is contained in the emitted code.
///   • PWM-BG (F6, Wasserman &amp; Sandelin 2004 / Biopython log_odds): W_q[b,j] = W_¼[b,j] + log2(¼ / q̂[b])
///     with q̂ the normalised background; a constant background equals the default overload.
///   • PWM-RC (F7, Biopython reverse_complement): RC is an involution and
///     score_RC(revcomp(s)) = score(s).
///   • PWM-BOTH (B05 follow-up, Biopython search both=True): the both-strand hits of s and of revcomp(s)
///     are mirror images — (i, strand, score) ↔ (n − m − i, opposite strand, same score); the plus subset
///     equals ScanWithPwm; CalculatePwmScores equals the all-window ScanWithPwm scores.
///   • PWM-DIST (Biopython thresholds.py): background density sums to 1; ThresholdFpr is non-increasing in
///     the FPR; ThresholdFnr is non-decreasing in the FNR; thresholds lie on the grid.
///   • REG-BOTH (B05 follow-up): FindRegulatoryElements(s, true) mirrors FindRegulatoryElements(revcomp(s), true)
///     for the both-strand elements (CAAT, GC box, NF-κB) and the palindromic ones (AP-1, E-box, CREB).
///   • DISCOVER-BG (F9/F10, RSAT oligo-analysis): a constant background reproduces the default
///     overload; log2(Enrichment) = log2(count) − Σ log2 q̂[w_i] − log2(N − k + 1), finite for long k.
///
/// All inputs seeded-random with fixed seeds; bounded sizes.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class MotifPwmMetamorphicTests
{
    private const string Acgt = "ACGT";

    private static string RandomString(Random rng, int length, string alphabet)
    {
        var c = new char[length];
        for (int i = 0; i < length; i++)
            c[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(c);
    }

    private static List<string> RandomAlignment(Random rng, int rows, int width)
        => Enumerable.Range(0, rows).Select(_ => RandomString(rng, width, Acgt)).ToList();

    private static string RevComp(string s)
        => new(s.Reverse().Select(c => c switch { 'A' => 'T', 'C' => 'G', 'G' => 'C', _ => 'A' }).ToArray());

    private static int Levenshtein(string a, string b)
    {
        var prev = new int[b.Length + 1];
        var curr = new int[b.Length + 1];
        for (int j = 0; j <= b.Length; j++) prev[j] = j;
        for (int i = 1; i <= a.Length; i++)
        {
            curr[0] = i;
            for (int j = 1; j <= b.Length; j++)
                curr[j] = Math.Min(prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1), Math.Min(prev[j] + 1, curr[j - 1] + 1));
            (prev, curr) = (curr, prev);
        }
        return prev[b.Length];
    }

    // NC-IUB 1984 code → base set.
    private static readonly Dictionary<char, string> Iupac = new()
    {
        ['A'] = "A", ['C'] = "C", ['G'] = "G", ['T'] = "T",
        ['R'] = "AG", ['Y'] = "CT", ['S'] = "CG", ['W'] = "AT", ['K'] = "GT", ['M'] = "AC",
        ['B'] = "CGT", ['D'] = "AGT", ['H'] = "ACT", ['V'] = "ACG", ['N'] = "ACGT",
    };

    #region EDIT-SELLERS (F2)

    [Test]
    public void FindWithEdits_EndSetAndMinDistance_EqualSellers()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(7_000 + seed);
            string alphabet = seed % 2 == 0 ? "AC" : Acgt;
            string text = RandomString(rng, rng.Next(1, 30), alphabet);
            string pattern = RandomString(rng, rng.Next(1, 7), alphabet);
            int k = rng.Next(0, 4);

            var windows = ApproximateMatcher.FindWithEdits(text, pattern, k).ToList();
            foreach (var w in windows)
            {
                Assert.That(w.MatchedSequence, Is.EqualTo(text.Substring(w.Position, w.MatchedSequence.Length)), $"seed={seed}");
                Assert.That(w.Distance, Is.EqualTo(Levenshtein(pattern, w.MatchedSequence)), $"seed={seed} window {w.MatchedSequence}");
                Assert.That(w.Distance, Is.LessThanOrEqualTo(k), $"seed={seed}");
            }

            var fromWindows = windows
                .GroupBy(w => w.Position + w.MatchedSequence.Length - 1)
                .Select(g => (EndPosition: g.Key, Distance: g.Min(w => w.Distance)))
                .OrderBy(e => e.EndPosition)
                .ToList();
            var sellers = ApproximateMatcher.FindEditEndPositions(text, pattern, k).ToList();

            Assert.That(fromWindows, Is.EqualTo(sellers), $"seed={seed} text={text} pattern={pattern} k={k}");
        }
    }

    #endregion

    #region CAVENER (F13c)

    [Test]
    public void CavenerConsensus_RowPermutationAndCaseInvariant_AlphabetIupac_ContainsMaxBases()
    {
        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(13_000 + seed);
            var aln = RandomAlignment(rng, rng.Next(1, 12), rng.Next(1, 15));
            string consensus = MotifFinder.GenerateCavenerConsensus(aln);

            Assert.That(consensus.Length, Is.EqualTo(aln[0].Length), $"seed={seed}");
            Assert.That(consensus.All(Iupac.ContainsKey), Is.True, $"seed={seed}: non-IUPAC symbol in {consensus}");

            var shuffled = aln.OrderBy(_ => rng.Next()).ToList();
            Assert.That(MotifFinder.GenerateCavenerConsensus(shuffled), Is.EqualTo(consensus), $"seed={seed}: row permutation");

            var mixedCase = aln.Select(s => new string(s.Select(c => rng.Next(2) == 0 ? char.ToLowerInvariant(c) : c).ToArray())).ToList();
            Assert.That(MotifFinder.GenerateCavenerConsensus(mixedCase), Is.EqualTo(consensus), $"seed={seed}: case");

            for (int j = 0; j < consensus.Length; j++)
            {
                var counts = Acgt.Select(b => aln.Count(s => s[j] == b)).ToArray();
                int max = counts.Max();
                for (int b = 0; b < 4; b++)
                    if (counts[b] == max)
                        Assert.That(Iupac[consensus[j]], Does.Contain(Acgt[b].ToString()),
                            $"seed={seed} col {j}: max-count base {Acgt[b]} not in {consensus[j]}");
            }
        }
    }

    [Test]
    public void CavenerConsensus_SingleSequence_IsTheSequence()
    {
        var rng = new Random(20260930);
        for (int i = 0; i < 100; i++)
        {
            string s = RandomString(rng, rng.Next(1, 40), "ACGTacgt");
            Assert.That(MotifFinder.GenerateCavenerConsensus(new[] { s }), Is.EqualTo(s.ToUpperInvariant()));
            // Duplicated rows do not change per-column proportions.
            Assert.That(MotifFinder.GenerateCavenerConsensus(new[] { s, s, s }), Is.EqualTo(s.ToUpperInvariant()));
        }
    }

    #endregion

    #region PWM-BG / PWM-RC (F6, F7)

    [Test]
    public void CreatePwm_BackgroundShiftsEachRowByLogOddsOfBackground()
    {
        for (int seed = 0; seed < 100; seed++)
        {
            var rng = new Random(17_000 + seed);
            var aln = RandomAlignment(rng, rng.Next(1, 10), rng.Next(1, 12));
            double p = seed % 5 == 0 ? 1.0 : 0.25 + rng.NextDouble();
            var reference = MotifFinder.CreatePwm(aln, p);

            double c = 0.1 + 5 * rng.NextDouble();
            var constant = MotifFinder.CreatePwm(aln, p, new[] { c, c, c, c });
            var unit = MotifFinder.CreatePwm(aln, p, new[] { 1.0, 1.0, 1.0, 1.0 });

            var raw = Enumerable.Range(0, 4).Select(_ => 0.05 + rng.NextDouble()).ToArray();
            double sum = raw.Sum();
            var bg = MotifFinder.CreatePwm(aln, p, raw);

            for (int j = 0; j < reference.Length; j++)
                for (int b = 0; b < 4; b++)
                {
                    Assert.That(unit.Matrix[b, j], Is.EqualTo(reference.Matrix[b, j]), $"seed={seed} unit background");
                    Assert.That(constant.Matrix[b, j], Is.EqualTo(reference.Matrix[b, j]).Within(1e-12), $"seed={seed} constant background");
                    double shift = Math.Log2(0.25 / (raw[b] / sum));
                    Assert.That(bg.Matrix[b, j], Is.EqualTo(reference.Matrix[b, j] + shift).Within(1e-9), $"seed={seed} b={b} j={j}");
                }
            Assert.That(unit.Consensus, Is.EqualTo(reference.Consensus), $"seed={seed}");
        }
    }

    [Test]
    public void PwmReverseComplement_IsInvolution_AndScoresReverseComplementEqually()
    {
        for (int seed = 0; seed < 100; seed++)
        {
            var rng = new Random(19_000 + seed);
            var aln = RandomAlignment(rng, rng.Next(1, 10), rng.Next(1, 12));
            var pwm = seed % 2 == 0
                ? MotifFinder.CreatePwm(aln, 0.5)
                : MotifFinder.CreatePwm(aln, 0.3, new[] { 0.3, 0.2, 0.2, 0.3 });
            var rc = pwm.ReverseComplement();
            var rcrc = rc.ReverseComplement();

            Assert.That(rc.Length, Is.EqualTo(pwm.Length));
            for (int j = 0; j < pwm.Length; j++)
                for (int b = 0; b < 4; b++)
                {
                    Assert.That(rcrc.Matrix[b, j], Is.EqualTo(pwm.Matrix[b, j]), $"seed={seed}: involution");
                    Assert.That(rc.Matrix[3 - b, pwm.Length - 1 - j], Is.EqualTo(pwm.Matrix[b, j]), $"seed={seed}: RC layout");
                }
            Assert.That(rc.MaxScore, Is.EqualTo(pwm.MaxScore).Within(1e-9));
            Assert.That(rc.MinScore, Is.EqualTo(pwm.MinScore).Within(1e-9));

            for (int t = 0; t < 5; t++)
            {
                string s = RandomString(rng, pwm.Length, Acgt);
                double forward = MotifFinder.ScanWithPwm(new DnaSequence(s), pwm, double.NegativeInfinity).Single().Score;
                double minus = MotifFinder.ScanWithPwm(new DnaSequence(RevComp(s)), rc, double.NegativeInfinity).Single().Score;
                Assert.That(minus, Is.EqualTo(forward).Within(1e-9), $"seed={seed} s={s}");
            }
        }
    }

    #endregion

    #region DISCOVER-BG (F9, F10)

    [Test]
    public void DiscoverMotifs_ConstantBackground_EqualsDefault()
    {
        for (int seed = 0; seed < 60; seed++)
        {
            var rng = new Random(23_000 + seed);
            var seq = new DnaSequence(RandomString(rng, rng.Next(10, 200), seed % 3 == 0 ? "AC" : Acgt));
            int k = rng.Next(1, 7);
            var reference = MotifFinder.DiscoverMotifs(seq, k, 2).ToList();
            var unit = MotifFinder.DiscoverMotifs(seq, k, 2, new[] { 2.0, 2.0, 2.0, 2.0 }).ToList();
            double c = 0.1 + rng.NextDouble();
            var constant = MotifFinder.DiscoverMotifs(seq, k, 2, new[] { c, c, c, c }).ToList();

            Assert.That(unit.Select(m => (m.Sequence, m.Count, m.Enrichment)),
                Is.EqualTo(reference.Select(m => (m.Sequence, m.Count, m.Enrichment))), $"seed={seed} (exact)");
            Assert.That(constant.Select(m => m.Sequence), Is.EqualTo(reference.Select(m => m.Sequence)), $"seed={seed}");
            for (int i = 0; i < reference.Count; i++)
            {
                Assert.That(constant[i].Positions, Is.EqualTo(reference[i].Positions), $"seed={seed}");
                Assert.That(constant[i].Enrichment, Is.EqualTo(reference[i].Enrichment).Within(1e-10).Percent, $"seed={seed}");
            }
        }
    }

    [Test]
    public void DiscoverMotifs_LongK_EnrichmentFinite_MatchesLog2Oracle()
    {
        var rng = new Random(29_000);
        var q = new[] { 0.3, 0.2, 0.2, 0.3 };
        foreach (int k in new[] { 64, 200, 333, 500 })
        {
            string x = RandomString(rng, k + rng.Next(0, 8), Acgt);
            var seq = new DnaSequence(x + x);
            double windows = seq.Length - k + 1.0;

            var uniform = MotifFinder.DiscoverMotifs(seq, k, 2).ToList();
            Assert.That(uniform, Is.Not.Empty, $"k={k}");
            foreach (var m in uniform)
            {
                Assert.That(double.IsFinite(m.Enrichment) && m.Enrichment > 0, Is.True, $"k={k}: {m.Enrichment}");
                double expectedLog2 = Math.Log2(m.Count) + 2.0 * k - Math.Log2(windows);
                Assert.That(Math.Log2(m.Enrichment), Is.EqualTo(expectedLog2).Within(1e-9 * expectedLog2), $"k={k}");
            }

            if (k > 300) continue; // Π q̂ ≥ 0.2^k: ratio ≤ 2^(2.33k) must stay representable
            foreach (var m in MotifFinder.DiscoverMotifs(seq, k, 2, q))
            {
                Assert.That(double.IsFinite(m.Enrichment) && m.Enrichment > 0, Is.True, $"k={k} bg: {m.Enrichment}");
                double sumLog = m.Sequence.Sum(ch => Math.Log2(q[Acgt.IndexOf(ch)]));
                double expectedLog2 = Math.Log2(m.Count) - sumLog - Math.Log2(windows);
                Assert.That(Math.Log2(m.Enrichment), Is.EqualTo(expectedLog2).Within(1e-9 * expectedLog2), $"k={k} bg");
            }
        }
    }

    #endregion

    #region PWM-BOTH / PWM-DIST / REG-BOTH (B05 follow-up)

    [Test]
    [Description("PWM-BOTH: both-strand hits of s mirror those of revcomp(s); plus subset = ScanWithPwm; calculate = all-window scan")]
    public void PwmBothStrands_ReverseComplementMirror([Values(31_001, 31_002, 31_003, 31_004, 31_005)] int seed)
    {
        var rng = new Random(seed);
        int width = rng.Next(3, 10);
        var pwm = MotifFinder.CreatePwm(RandomAlignment(rng, rng.Next(2, 12), width), 0.25 + rng.NextDouble());
        string s = RandomString(rng, rng.Next(width, 120), Acgt);
        var seq = new DnaSequence(s);
        var rcSeq = seq.ReverseComplement();
        int n = s.Length;
        double threshold = rng.NextDouble() * 4 - 2;

        var a = MotifFinder.ScanWithPwmBothStrands(seq, pwm, threshold).ToList();
        var b = MotifFinder.ScanWithPwmBothStrands(rcSeq, pwm, threshold).ToList();
        var mirrored = b.Select(h => (Pos: n - width - h.Position, Strand: h.Strand == '+' ? '-' : '+', h.MatchedSequence, h.Score))
            .OrderBy(h => h.Pos).ThenBy(h => h.Strand == '+' ? 0 : 1).ToList();

        Assert.That(a.Select(h => (h.Position, h.Strand, h.MatchedSequence)),
            Is.EqualTo(mirrored.Select(h => (h.Pos, h.Strand, h.MatchedSequence))), $"seed {seed}");
        Assert.That(a.Select(h => h.Score), Is.EqualTo(mirrored.Select(h => h.Score)).Within(1e-9));
        Assert.That(a.Where(h => h.Strand == '+').Select(h => (h.Position, h.Score)),
            Is.EqualTo(MotifFinder.ScanWithPwm(seq, pwm, threshold).Select(m => (m.Position, m.Score))));
        Assert.That(MotifFinder.CalculatePwmScores(seq, pwm),
            Is.EqualTo(MotifFinder.ScanWithPwm(seq, pwm, double.NegativeInfinity).Select(m => m.Score)));
    }

    [Test]
    [Description("PWM-DIST: densities normalised, ThresholdFpr non-increasing, ThresholdFnr non-decreasing, on-grid")]
    public void PwmScoreDistribution_Monotone([Values(32_001, 32_002, 32_003)] int seed)
    {
        var rng = new Random(seed);
        // Motif density q·2^W is a distribution when the PWM is built against the same background q.
        var q = new[] { 0.2 + rng.NextDouble(), 0.2 + rng.NextDouble(), 0.2 + rng.NextDouble(), 0.2 + rng.NextDouble() };
        var pwm = MotifFinder.CreatePwm(RandomAlignment(rng, rng.Next(3, 10), rng.Next(3, 9)), 0.5, q);
        var d = pwm.ScoreDistribution(q, 200);
        Assert.That(d.BackgroundDensity.Sum(), Is.EqualTo(1.0).Within(1e-9));
        Assert.That(d.MotifDensity.Sum(), Is.EqualTo(1.0).Within(1e-6));

        double[] rates = { 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5 };
        var fpr = rates.Select(d.ThresholdFpr).ToArray();
        var fnr = rates.Select(d.ThresholdFnr).ToArray();
        for (int i = 1; i < rates.Length; i++)
        {
            Assert.That(fpr[i], Is.LessThanOrEqualTo(fpr[i - 1]));
            Assert.That(fnr[i], Is.GreaterThanOrEqualTo(fnr[i - 1]));
        }
        foreach (double t in fpr.Concat(fnr))
        {
            double k = (t - d.MinScore) / d.Step;
            Assert.That(k, Is.EqualTo(Math.Round(k)).Within(1e-6), "threshold is a grid point");
        }
    }

    [Test]
    [Description("REG-BOTH: FindRegulatoryElements(s, true) mirrors FindRegulatoryElements(revcomp(s), true) for both-strand/palindromic elements")]
    public void RegulatoryBothStrands_ReverseComplementMirror([Values(33_001, 33_002, 33_003, 33_004)] int seed)
    {
        var rng = new Random(seed);
        string[] plants = { "CCAAT", "ATTGG", "GGGCGG", "CCGCCC", "GGGACTTTCC", "GGAAAGTCCC", "TGAGTCA", "CACGTG", "TGACGTCA" };
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i < 12; i++)
            sb.Append(RandomString(rng, rng.Next(0, 8), Acgt)).Append(plants[rng.Next(plants.Length)]);
        string s = sb.ToString();
        int n = s.Length;
        var mirrorable = new HashSet<string> { "CAAT Box", "GC Box", "NF-κB", "AP-1", "E-box", "CREB" };

        var a = MotifFinder.FindRegulatoryElements(new DnaSequence(s), true).Where(e => mirrorable.Contains(e.Name)).ToList();
        var b = MotifFinder.FindRegulatoryElements(new DnaSequence(s).ReverseComplement(), true).Where(e => mirrorable.Contains(e.Name)).ToList();
        bool palindromic(string name) => name is "AP-1" or "E-box" or "CREB";
        var mirrored = b.Select(e => (e.Name, Pos: n - e.Sequence.Length - e.Position,
                Strand: palindromic(e.Name) ? '+' : (e.Strand == '+' ? '-' : '+'),
                Seq: palindromic(e.Name) ? DnaSequence.GetReverseComplementString(e.Sequence) : e.Sequence))
            .OrderBy(e => Array.IndexOf(MotifFinder.OrientationIndependentRegulatoryElements.ToArray(), e.Name))
            .ThenBy(e => e.Pos).ThenBy(e => e.Strand == '+' ? 0 : 1);
        var actual = a.Select(e => (e.Name, Pos: e.Position, e.Strand, Seq: e.Sequence))
            .OrderBy(e => Array.IndexOf(MotifFinder.OrientationIndependentRegulatoryElements.ToArray(), e.Name))
            .ThenBy(e => e.Pos).ThenBy(e => e.Strand == '+' ? 0 : 1);
        Assert.That(actual, Is.EqualTo(mirrored), $"seed {seed}: {s}");
        Assert.That(a, Is.Not.Empty);
    }

    #endregion
}
