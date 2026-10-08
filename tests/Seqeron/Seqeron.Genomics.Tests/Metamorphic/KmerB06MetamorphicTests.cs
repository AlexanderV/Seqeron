namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// B06 heavy tier — metamorphic relations for the k-mer distance and sketch behaviour added by batch B06
/// (docs/Validation/review-2026-09/B06.md): the word-vector metrics (F11, F15, F19 D2*/D2S incl. CAFE <c>-R</c>, the
/// d2* clamp of F31), canonical Jaccard / Mash / containment (F11, F24), Mash and FracMinHash sketches (F23, F30),
/// multiple-pattern spaced words incl. <c>spaced</c>'s reverse-complement mode and <c>-d EV</c> (F22, F27, F28), and
/// kPAL ACGT-only both-strand counting (F25).
///
/// Relations (derived from the definitions in the XML documentation, not from observed output):
///   • SYM  — every metric whose formula is symmetric in (X, Y) gives d(a, b) = d(b, a): Euclidean, squared Euclidean
///            on counts, Manhattan, Chebyshev, Canberra, cosine, D2, Jensen–Shannon, Euclidean on counts, D2*, D2S
///            (single- and both-strand, any Markov order incl. BIC), canonical Jaccard and Mash distance.
///   • ID   — identical inputs: every dissimilarity is 0 (d2*/d2S ≈ 0; F31 clamps the −1e-16 rounding to 0).
///   • RC   — canonical k-mer sets, Mash/FracMinHash canonical sketches, the CAFE <c>-R</c> counts and order-0 background,
///            the kPAL both-strand table and <c>spaced</c>'s both-strand seq1 vector are unchanged when the relevant
///            sequence is replaced by its reverse complement.
///   • EQV  — the single all-'1' spaced pattern of length k is the contiguous k-mer distance (Leimeister et al. 2014).
///   • MON  — the Mash distance −ln(2J/(1+J))/k is non-increasing in J.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class KmerB06MetamorphicTests
{
    #region Helpers

    private static readonly KmerDistanceMetric[] SymmetricWordVectorMetrics =
    [
        KmerDistanceMetric.Euclidean, KmerDistanceMetric.SquaredEuclideanCounts, KmerDistanceMetric.Manhattan,
        KmerDistanceMetric.Chebyshev, KmerDistanceMetric.Canberra, KmerDistanceMetric.Cosine, KmerDistanceMetric.D2,
        KmerDistanceMetric.JensenShannon, KmerDistanceMetric.EuclideanCounts,
    ];

    private static string RandomSeq(Random rng, int length, string alphabet)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    private static string Related(Random rng, string s, int every, string alphabet)
    {
        var chars = s.ToCharArray();
        for (int i = 0; i < chars.Length; i++)
            if (rng.Next(every) == 0)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    private static string Rc(string s) => DnaSequence.GetReverseComplementString(s);

    private static void ShouldBeClose(double actual, double expected, double tol, string because)
    {
        if (double.IsNaN(expected))
        {
            double.IsNaN(actual).Should().BeTrue(because);
            return;
        }
        if (double.IsInfinity(expected))
        {
            actual.Should().Be(expected, because);
            return;
        }
        actual.Should().BeApproximately(expected, tol * Math.Max(1.0, Math.Abs(expected)), because);
    }

    #endregion

    #region SYM / ID — word-vector metrics (F11)

    /// <summary>
    /// SYM + ID: the nine word-vector metrics are symmetric functions of the two vectors over the union of keys, and every
    /// dissimilarity is 0 for identical inputs (cosine: 0 for a non-zero vector; D2 is a similarity, so only SYM applies).
    /// </summary>
    [Test]
    public void WordVectorMetrics_AreSymmetric_AndZeroOnIdenticalInputs()
    {
        var rng = new Random(606001);
        for (int iter = 0; iter < 120; iter++)
        {
            string a = RandomSeq(rng, rng.Next(0, 160), "ACGTacgtN");
            string b = rng.Next(3) == 0 ? RandomSeq(rng, rng.Next(0, 160), "ACGTN") : Related(rng, a, 5, "ACGT");
            int k = rng.Next(1, 6);
            foreach (var metric in SymmetricWordVectorMetrics)
            {
                double d = KmerAnalyzer.KmerDistance(a, b, k, metric);
                double dSwap = KmerAnalyzer.KmerDistance(b, a, k, metric);
                ShouldBeClose(dSwap, d, 1e-12, $"SYM: {metric} k={k}");
                if (metric == KmerDistanceMetric.D2)
                    continue;
                double self = KmerAnalyzer.KmerDistance(a, a, k, metric);
                bool zeroVector = a.Length < k;
                if (metric == KmerDistanceMetric.Cosine && zeroVector)
                    self.Should().Be(1.0, "a zero vector has cosine similarity 0 (documented convention)");
                else
                    self.Should().BeApproximately(0.0, 1e-12, $"ID: {metric} k={k} d(a, a) = 0");
            }
        }
    }

    #endregion

    #region SYM / ID / RC — D2* and D2S (F15, F19, F31)

    /// <summary>
    /// SYM + ID: d2* and d2S (Reinert et al. 2009; Song et al. 2014) are symmetric in (X̃, Ỹ) for any Markov order
    /// (including the per-sequence BIC choice −1, which swaps with the sequences) and in both strand modes; identical
    /// sequences give ≈ 0 (Cauchy–Schwarz equality), never below 0 (F31 clamp); raw D2*/D2S are symmetric too.
    /// </summary>
    [Test]
    public void D2StarAndShepherd_AreSymmetric_AndZeroOnIdenticalSequences()
    {
        var rng = new Random(606002);
        for (int iter = 0; iter < 60; iter++)
        {
            string a = RandomSeq(rng, rng.Next(40, 220), "ACGTACGTacgtN");
            string b = rng.Next(3) == 0 ? RandomSeq(rng, rng.Next(40, 220), "ACGT") : Related(rng, a, 4, "ACGT");
            int k = rng.Next(1, 6);
            int order = rng.Next(4) == 0 ? KmerAnalyzer.AutoMarkovOrder : rng.Next(0, Math.Min(k, 3));
            bool both = rng.Next(2) == 0;

            var ab = KmerAnalyzer.BackgroundAdjustedD2(a, b, k, order, both);
            var ba = KmerAnalyzer.BackgroundAdjustedD2(b, a, k, order, both);
            ShouldBeClose(ba.D2StarDistance, ab.D2StarDistance, 1e-9, $"SYM d2* k={k} r={order} both={both}");
            ShouldBeClose(ba.D2ShepherdDistance, ab.D2ShepherdDistance, 1e-9, $"SYM d2S k={k} r={order} both={both}");
            ShouldBeClose(ba.D2Star, ab.D2Star, 1e-9, "SYM raw D2*");
            ShouldBeClose(ba.D2Shepherd, ab.D2Shepherd, 1e-9, "SYM raw D2S");
            ba.MarkovOrder1.Should().Be(ab.MarkovOrder2, "the per-sequence orders swap with the sequences");

            KmerAnalyzer.KmerDistance(a, b, k, KmerDistanceMetric.D2Star, order, both).Should().Be(ab.D2StarDistance,
                "the metric overload delegates to BackgroundAdjustedD2");

            var self = KmerAnalyzer.BackgroundAdjustedD2(a, a, k, order, both);
            if (!double.IsNaN(self.D2StarDistance))
                self.D2StarDistance.Should().BeInRange(0.0, 1e-9, "ID: d2*(a, a) ≈ 0 and clamped at 0 (F31)");
            if (!double.IsNaN(self.D2ShepherdDistance))
                self.D2ShepherdDistance.Should().BeInRange(0.0, 1e-9, "ID: d2S(a, a) ≈ 0 and clamped at 0 (F31)");
            foreach (var d in new[] { ab.D2StarDistance, ab.D2ShepherdDistance })
                if (!double.IsNaN(d))
                    d.Should().BeInRange(0.0, 1.0, "d2*/d2S are dissimilarities in [0, 1]");
        }
    }

    /// <summary>
    /// RC: in CAFE's both-strand mode (F19) the counts are X(w) + X(RC w) and, for order 0, the background
    /// p^R(w) = ½(p̂(w) + p̂(RC w)) with p̂(RC w) computed on RC(s) equal to p̂ on s; so replacing either sequence by its
    /// reverse complement leaves d2*, d2S, D2* and D2S unchanged, and d^R(a, RC a) ≈ 0.
    /// </summary>
    [Test]
    public void BothStrandD2_IsInvariantUnderReverseComplementOfEitherSequence()
    {
        var rng = new Random(606003);
        for (int iter = 0; iter < 60; iter++)
        {
            string a = RandomSeq(rng, rng.Next(40, 220), "ACGTACGTacgtN");
            string b = Related(rng, a, 3, "ACGT");
            int k = rng.Next(1, 6);
            var baseline = KmerAnalyzer.BackgroundAdjustedD2(a, b, k, 0, bothStrands: true);
            foreach (var (x, y, label) in new[] { (Rc(a), b, "RC(a)"), (a, Rc(b), "RC(b)"), (Rc(a), Rc(b), "RC both") })
            {
                var v = KmerAnalyzer.BackgroundAdjustedD2(x, y, k, 0, bothStrands: true);
                ShouldBeClose(v.D2StarDistance, baseline.D2StarDistance, 1e-9, $"RC d2* {label} k={k}");
                ShouldBeClose(v.D2ShepherdDistance, baseline.D2ShepherdDistance, 1e-9, $"RC d2S {label} k={k}");
                ShouldBeClose(v.D2Star, baseline.D2Star, 1e-9, $"RC D2* {label}");
                ShouldBeClose(v.D2Shepherd, baseline.D2Shepherd, 1e-9, $"RC D2S {label}");
            }
            var selfRc = KmerAnalyzer.BackgroundAdjustedD2(a, Rc(a), k, 0, bothStrands: true);
            if (!double.IsNaN(selfRc.D2StarDistance))
                selfRc.D2StarDistance.Should().BeInRange(0.0, 1e-9, "a sequence and its reverse complement have equal both-strand profiles");
        }
    }

    #endregion

    #region RC / SYM / MON — canonical Jaccard, Mash, sketches (F11, F23, F24, F30)

    /// <summary>
    /// RC: the canonical k-mer set of RC(s) is that of s, so the canonical Jaccard index, containment and Mash distance,
    /// the Mash bottom-s sketch and the FracMinHash sketch (with abundances) are unchanged under reverse complement of
    /// either input; SYM for Jaccard and the Mash distance.
    /// </summary>
    [Test]
    public void CanonicalJaccardMashAndSketches_AreInvariantUnderReverseComplement()
    {
        var rng = new Random(606004);
        var canonical = new KmerCountingOptions(Canonical: true);
        for (int iter = 0; iter < 80; iter++)
        {
            string a = RandomSeq(rng, rng.Next(0, 250), "ACGTACGTacgtNRY");
            string b = Related(rng, a, 6, "ACGT");
            int k = rng.Next(1, 32);
            double j = KmerAnalyzer.JaccardSimilarity(a, b, k, canonical);
            double m = KmerAnalyzer.MashDistance(a, b, k, canonical);
            double c = KmerAnalyzer.ContainmentIndex(a, b, k, canonical);
            KmerAnalyzer.JaccardSimilarity(Rc(a), b, k, canonical).Should().Be(j, "RC(a) has the same canonical set");
            KmerAnalyzer.JaccardSimilarity(a, Rc(b), k, canonical).Should().Be(j, "RC(b) has the same canonical set");
            KmerAnalyzer.JaccardSimilarity(b, a, k, canonical).Should().Be(j, "SYM Jaccard");
            KmerAnalyzer.MashDistance(Rc(a), Rc(b), k, canonical).Should().Be(m, "RC Mash distance");
            KmerAnalyzer.MashDistance(b, a, k, canonical).Should().Be(m, "SYM Mash distance");
            KmerAnalyzer.ContainmentIndex(Rc(a), b, k, canonical).Should().Be(c, "RC containment");

            if (k <= KmerAnalyzer.MaxMashKmerSize)
            {
                var sa = KmerAnalyzer.CreateMinHashSketch(a, k, 200);
                var sRc = KmerAnalyzer.CreateMinHashSketch(Rc(a), k, 200);
                sRc.Hashes.Should().Equal(sa.Hashes, "the canonical Mash sketch depends only on the canonical set");
                var cmp = KmerAnalyzer.CompareMinHashSketches(sa, KmerAnalyzer.CreateMinHashSketch(b, k, 200));
                var cmpRc = KmerAnalyzer.CompareMinHashSketches(sRc, KmerAnalyzer.CreateMinHashSketch(Rc(b), k, 200));
                cmpRc.Should().Be(cmp, "RC of both inputs leaves the mash dist line unchanged");
            }

            var fa = KmerAnalyzer.CreateFracMinHashSketch(a, k, 3, trackAbundance: true);
            var fRc = KmerAnalyzer.CreateFracMinHashSketch(Rc(a), k, 3, trackAbundance: true);
            fRc.Hashes.Should().Equal(fa.Hashes, "FracMinHash keeps the hashes of the canonical k-mers");
            fRc.Abundances.Should().Equal(fa.Abundances, "the canonical counts of RC(s) equal those of s");
        }
    }

    /// <summary>
    /// MON: D(J) = −ln(2J/(1+J))/k is strictly decreasing in J on (0, 1] (derivative −1/(kJ(1+J)) &lt; 0), with D(0) = 1,
    /// D(1) = 0 and the cap at 1, so J₁ ≤ J₂ ⇒ D(J₁) ≥ D(J₂); on sequences, a pair with the larger exact canonical Jaccard
    /// has the smaller (or equal) Mash distance.
    /// </summary>
    [Test]
    public void MashDistance_IsMonotoneNonIncreasingInJaccard()
    {
        var rng = new Random(606005);
        for (int iter = 0; iter < 2000; iter++)
        {
            int k = rng.Next(1, 33);
            double j1 = rng.Next(5) == 0 ? 0.0 : rng.NextDouble();
            double j2 = rng.Next(5) == 0 ? 1.0 : j1 + (1.0 - j1) * rng.NextDouble();
            double d1 = KmerAnalyzer.MashDistanceFromJaccard(j1, k);
            double d2 = KmerAnalyzer.MashDistanceFromJaccard(j2, k);
            d1.Should().BeGreaterThanOrEqualTo(d2, $"MON: J {j1} ≤ {j2} at k={k}");
            d1.Should().BeInRange(0.0, 1.0);
        }

        var canonical = new KmerCountingOptions(Canonical: true);
        for (int iter = 0; iter < 60; iter++)
        {
            string a = RandomSeq(rng, rng.Next(20, 200), "ACGT");
            string b = Related(rng, a, rng.Next(2, 30), "ACGT");
            string c = Related(rng, a, rng.Next(2, 30), "ACGT");
            int k = rng.Next(3, 22);
            double jb = KmerAnalyzer.JaccardSimilarity(a, b, k, canonical);
            double jc = KmerAnalyzer.JaccardSimilarity(a, c, k, canonical);
            double db = KmerAnalyzer.MashDistance(a, b, k, canonical);
            double dc = KmerAnalyzer.MashDistance(a, c, k, canonical);
            if (jb <= jc)
                db.Should().BeGreaterThanOrEqualTo(dc, $"MON on sequences: J {jb} ≤ {jc}");
            else
                db.Should().BeLessThanOrEqualTo(dc, $"MON on sequences: J {jb} > {jc}");
        }
    }

    #endregion

    #region EQV / RC / SYM — spaced words (F22, F27, F28)

    /// <summary>
    /// EQV: with the single all-'1' pattern of length k the spaced word at window i is the k-mer at i, so the
    /// multiple-pattern spaced distance equals <c>KmerDistance</c> at k for every word-vector metric (literal words; also
    /// ACGT-only on gap-free ACGT input, where no word is dropped); EV through <c>KmerDistance</c> equals the spaced
    /// EV distance with that pattern in both strand modes.
    /// </summary>
    [Test]
    public void SpacedAllOnesPattern_EqualsContiguousKmerDistance()
    {
        var rng = new Random(606006);
        for (int iter = 0; iter < 80; iter++)
        {
            string a = RandomSeq(rng, rng.Next(0, 150), "ACGTacgtN");
            string b = Related(rng, a, 4, "ACGTN");
            int k = rng.Next(1, 6);
            string[] ones = [new string('1', k)];
            foreach (var metric in SymmetricWordVectorMetrics)
            {
                ShouldBeClose(KmerAnalyzer.SpacedWordDistance(a, b, ones, metric), KmerAnalyzer.KmerDistance(a, b, k, metric),
                    1e-12, $"EQV literal {metric} k={k}");
            }

            string pa = RandomSeq(rng, rng.Next(0, 150), "ACGTacgt");
            string pb = Related(rng, pa, 4, "ACGT");
            foreach (var metric in SymmetricWordVectorMetrics)
            {
                ShouldBeClose(
                    KmerAnalyzer.SpacedWordDistance(pa, pb, ones, metric, new KmerCountingOptions(AcgtOnly: true)),
                    KmerAnalyzer.KmerDistance(pa, pb, k, metric), 1e-12, $"EQV ACGT-only {metric} k={k}");
            }

            if (pa.Length >= k && pb.Length >= k)
            {
                foreach (bool both in new[] { false, true })
                {
                    ShouldBeClose(
                        KmerAnalyzer.KmerDistance(pa, pb, k, KmerDistanceMetric.SpacedEvolutionary, 0, both),
                        KmerAnalyzer.SpacedWordDistance(pa, pb, ones, KmerDistanceMetric.SpacedEvolutionary,
                            new KmerCountingOptions(AcgtOnly: true), both),
                        0, $"EQV EV both={both}");
                }
            }
        }
    }

    /// <summary>
    /// F27 argument-order contract: in <c>spaced</c>'s reverse-complement mode seq1 is read on both strands (vector
    /// F₁ + R₁ over 2·W₁ windows) and seq2 on its forward strand only, so the value is unchanged when seq1 is replaced
    /// by RC(seq1) (F and R swap) — for every word-vector metric and for EV; with bothStrands = false the symmetric
    /// metrics are symmetric under swapping the arguments (multi-pattern averages included).
    /// </summary>
    [Test]
    public void SpacedBothStrands_DependsOnSeq1OnlyThroughBothStrands_SingleStrandIsSymmetric()
    {
        var rng = new Random(606007);
        string[][] patternSets = [["1101"], ["1101", "1011"], ["11011", "10111", "11101"], ["1"], ["101"]];
        var acgt = new KmerCountingOptions(AcgtOnly: true);
        for (int iter = 0; iter < 80; iter++)
        {
            string a = RandomSeq(rng, rng.Next(5, 160), "ACGTACGTacgtN");
            string b = Related(rng, rng.Next(2) == 0 ? a : Rc(a), 5, "ACGTN");
            var patterns = patternSets[rng.Next(patternSets.Length)];
            foreach (var metric in SymmetricWordVectorMetrics.Append(KmerDistanceMetric.SpacedEvolutionary))
            {
                double d = KmerAnalyzer.SpacedWordDistance(a, b, patterns, metric, acgt, bothStrands: true);
                double dRc = KmerAnalyzer.SpacedWordDistance(Rc(a), b, patterns, metric, acgt, bothStrands: true);
                ShouldBeClose(dRc, d, 1e-12, $"RC(seq1) in both-strand mode, {metric} [{string.Join(',', patterns)}]");

                if (metric == KmerDistanceMetric.SpacedEvolutionary)
                    continue;
                foreach (var options in new[] { KmerCountingOptions.Default, acgt })
                {
                    double s = KmerAnalyzer.SpacedWordDistance(a, b, patterns, metric, options);
                    double sSwap = KmerAnalyzer.SpacedWordDistance(b, a, patterns, metric, options);
                    ShouldBeClose(sSwap, s, 1e-12, $"SYM single-strand {metric} {options}");
                }
            }
        }
    }

    #endregion

    #region RC — kPAL ACGT-only both-strand counting (F25)

    /// <summary>
    /// F25: the kPAL both-strand table count[w] = f[w] + f[RC w] is unchanged under reverse complement of the input
    /// (literal and ACGT-only), and on canonical keys it relates to Jellyfish <c>-C</c> by the documented palindrome rule:
    /// equal for w ≠ RC(w), doubled for a palindrome.
    /// </summary>
    [Test]
    public void BothStrandCounts_InvariantUnderReverseComplement_AndMatchCanonicalUpToPalindromes()
    {
        var rng = new Random(606008);
        var acgt = new KmerCountingOptions(AcgtOnly: true);
        for (int iter = 0; iter < 150; iter++)
        {
            string s = RandomSeq(rng, rng.Next(0, 200), "ACGTACGTacgtNRY");
            int k = rng.Next(1, 7);
            foreach (var options in new[] { KmerCountingOptions.Default, acgt })
            {
                KmerAnalyzer.CountKmersBothStrands(Rc(s), k, options)
                    .Should().BeEquivalentTo(KmerAnalyzer.CountKmersBothStrands(s, k, options), $"RC invariance {options} k={k}");
            }
            var both = KmerAnalyzer.CountKmersBothStrands(s, k, acgt);
            foreach (var (w, n) in KmerAnalyzer.CountKmers(s, k, new KmerCountingOptions(Canonical: true)))
            {
                int expected = w == Rc(w) ? 2 * n : n;
                both[w].Should().Be(expected, $"canonical {w}: kPAL balance vs Jellyfish -C (palindromes doubled)");
            }
        }
    }

    #endregion
}
