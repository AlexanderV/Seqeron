// KMER-COUNT-001 / KMER-ASYNC-001 / KMER-DIST-001 (audit round 1, WP4) — single span-lookup counting loop,
// DnaSequence null arguments, opt-in parallel counting, background-adjusted D2* / D2S.
// Evidence: docs/Evidence/KMER-DIST-001-Evidence.md; TestSpecs: KMER-COUNT-001.md, KMER-ASYNC-001.md, KMER-DIST-001.md
// Sources: .NET Framework Design Guidelines (argument validation: ArgumentNullException); Marçais & Kingsford 2011
//          (Jellyfish multi-threaded counting); Reinert, Chew, Sun & Waterman 2009 (J Comput Biol 16:1615) D2*/D2S;
//          Song et al. 2014 (Brief Bioinform 15:343) d2*/d2S dissimilarities; CAFE (Lu et al. 2017) dist_model.cpp.
// Reference values: independent Python replica of the published formulas (order-r Markov background fitted per
//          sequence), whose formula engine reproduces the CAFE binary (built from github.com/younglululu/CAFE,
//          with Jellyfish 2.3.1) to 6 digits when given CAFE's prefix-marginal estimator.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_ParallelAndBackgroundD2_Tests
{
    private const double Tol = 1e-12;

    // Python random.seed(2026): S1 uniform 80 nt, S2 A/T-rich 70 nt, S3 mixed case with N/R (ACGT windows only).
    private const string S1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
    private const string S2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
    private const string S3 = "acgtNNacgtacgRtTTGCAnAGGATTACAggattacaTTTAAACCCGGGTTANcaataccgagcccatggagc";

    private static string RandomSequence(Random rnd, int length, string alphabet)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rnd.Next(alphabet.Length)];
        return new string(chars);
    }

    /// <summary>Independent naive counter: substring per window, upper-cased, optional Jellyfish ACGT rule.</summary>
    private static Dictionary<string, int> NaiveCount(string s, int k, bool acgtOnly)
    {
        var counts = new Dictionary<string, int>();
        var u = s.ToUpperInvariant();
        for (int i = 0; i + k <= u.Length; i++)
        {
            var w = u.Substring(i, k);
            if (acgtOnly && w.Any(c => c is not ('A' or 'C' or 'G' or 'T')))
                continue;
            counts[w] = counts.TryGetValue(w, out var c) ? c + 1 : 1;
        }
        return counts;
    }

    /// <summary>O(n) table equality (NUnit's Is.EquivalentTo is quadratic on large dictionaries).</summary>
    private static void AssertSameCounts(IReadOnlyDictionary<string, int> actual, IReadOnlyDictionary<string, int> expected, string message)
    {
        Assert.That(actual.Count, Is.EqualTo(expected.Count), message);
        foreach (var (kmer, count) in expected)
        {
            if (!actual.TryGetValue(kmer, out var c) || c != count)
                Assert.Fail($"{message}: {kmer} expected {count}, got {(actual.ContainsKey(kmer) ? c : "missing")}");
        }
    }

    #region Single counting loop (span lookup) — unchanged results, CountKmersSpan contract

    [Test]
    public void CountKmers_SpanLookupLoop_EqualsNaiveCounter_OnRandomInputs()
    {
        var rnd = new Random(4);
        for (int t = 0; t < 300; t++)
        {
            var s = RandomSequence(rnd, rnd.Next(0, 3000), "ACGTacgtNR");
            int k = rnd.Next(1, 9);
            AssertSameCounts(KmerAnalyzer.CountKmers(s, k), NaiveCount(s, k, acgtOnly: false), $"literal t={t}");
            AssertSameCounts(KmerAnalyzer.CountKmers(s, k, new KmerCountingOptions(AcgtOnly: true)),
                NaiveCount(s, k, acgtOnly: true), $"acgt t={t}");
        }
    }

    [Test]
    public void CountKmersSpan_EqualsCountKmers_IncludingLowerCase()
    {
        const string s = "acgtNNacgtacgRtTTGCAnA";
        for (int k = 1; k <= s.Length + 1; k++)
            Assert.That(KmerAnalyzer.CountKmersSpan(s.AsSpan(), k), Is.EquivalentTo(KmerAnalyzer.CountKmers(s, k)), $"k={k}");
    }

    [TestCase(0)]
    [TestCase(-3)]
    public void CountKmersSpan_EmptySpan_ReturnsEmptyForAnyK_LikeCountKmers(int k)
    {
        Assert.That(KmerAnalyzer.CountKmersSpan(ReadOnlySpan<char>.Empty, k), Is.Empty);
        Assert.That(KmerAnalyzer.CountKmers(string.Empty, k), Is.Empty);
    }

    [TestCase(0)]
    [TestCase(-1)]
    public void CountKmersSpan_NonEmpty_NonPositiveK_Throws(int k)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.CountKmersSpan("ACGT".AsSpan(), k));
        Assert.That(ex!.ParamName, Is.EqualTo("k"));
    }

    #endregion

    #region Null DnaSequence arguments

    [Test]
    public void DnaSequenceOverloads_Null_ThrowArgumentNullException()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.CountKmers((DnaSequence)null!, 3))!.ParamName,
                Is.EqualTo("dna"));
            Assert.That(Assert.Throws<ArgumentNullException>(
                    () => KmerAnalyzer.CountKmers((DnaSequence)null!, 3, CancellationToken.None))!.ParamName,
                Is.EqualTo("dna"));
            Assert.That(Assert.Throws<ArgumentNullException>(
                    () => KmerAnalyzer.CountKmersBothStrands((DnaSequence)null!, 3))!.ParamName,
                Is.EqualTo("dna"));
        });
    }

    #endregion

    #region CountKmersParallel — exact equality with the serial count, cancellation, progress

    [Test]
    public void CountKmersParallel_EqualsSerialCount_OnRandomInputs_AllOptions()
    {
        var rnd = new Random(11);
        var optionSets = new[] { KmerCountingOptions.Default, new KmerCountingOptions(AcgtOnly: true), new KmerCountingOptions(Canonical: true) };
        for (int t = 0; t < 8; t++)
        {
            // 140k–420k windows: 2–6 ranges of >= 65,536 windows each (the partitioned path).
            var s = RandomSequence(rnd, rnd.Next(140_000, 420_000), t % 2 == 0 ? "ACGT" : "ACGTacgtN");
            int k = rnd.Next(1, 14);
            foreach (var options in optionSets)
            {
                var serial = KmerAnalyzer.CountKmers(s, k, options);
                foreach (int degree in new[] { 2, 3, 4, 8, -1 })
                    AssertSameCounts(KmerAnalyzer.CountKmersParallel(s, k, options, degree), serial,
                        $"t={t} k={k} options={options} degree={degree}");
            }
        }
    }

    [Test]
    public void CountKmersParallel_SmallInputsAndDegreeOne_UseSerialPath_SameResult()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.CountKmersParallel("GTAGAGCTGT", 2), Is.EquivalentTo(KmerAnalyzer.CountKmers("GTAGAGCTGT", 2)));
            Assert.That(KmerAnalyzer.CountKmersParallel(null!, 3), Is.Empty);
            Assert.That(KmerAnalyzer.CountKmersParallel(string.Empty, 0), Is.Empty);
            Assert.That(KmerAnalyzer.CountKmersParallel("ACG", 5), Is.Empty);
            var big = RandomSequence(new Random(3), 200_000, "ACGT");
            AssertSameCounts(KmerAnalyzer.CountKmersParallel(big, 6, default, 1), KmerAnalyzer.CountKmers(big, 6), "degree 1");
        });
    }

    [TestCase(0)]
    [TestCase(-2)]
    public void CountKmersParallel_InvalidDegree_Throws(int degree)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.CountKmersParallel("ACGT", 2, default, degree));
        Assert.That(ex!.ParamName, Is.EqualTo("maxDegreeOfParallelism"));
    }

    [Test]
    public void CountKmersParallel_NonPositiveK_Throws()
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.CountKmersParallel("ACGT", 0));
        Assert.That(ex!.ParamName, Is.EqualTo("k"));
    }

    [Test]
    public void CountKmersParallel_PreCanceled_ThrowsWithToken_NoProgress()
    {
        var big = RandomSequence(new Random(5), 300_000, "ACGT");
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        var reports = new List<double>();
        var ex = Assert.Catch<OperationCanceledException>(
            () => KmerAnalyzer.CountKmersParallel(big, 5, default, 4, cts.Token, new SyncProgress(reports.Add)));
        Assert.That(ex!.CancellationToken, Is.EqualTo(cts.Token));
        Assert.That(reports, Is.Empty);
    }

    [Test]
    public void CountKmersParallel_CanceledMidRun_ThrowsWithToken()
    {
        var big = RandomSequence(new Random(6), 400_000, "ACGT");
        using var cts = new CancellationTokenSource();
        var ex = Assert.Catch<OperationCanceledException>(() => KmerAnalyzer.CountKmersParallel(
            big, 7, default, 4, cts.Token, new SyncProgress(_ => cts.Cancel())));
        Assert.That(ex!.CancellationToken, Is.EqualTo(cts.Token));
    }

    [Test]
    public void CountKmersParallel_Progress_NonDecreasing_BelowOne_ThenSingleFinalOne()
    {
        var big = RandomSequence(new Random(8), 400_000, "ACGT");
        var reports = new List<double>();
        var gate = new object();
        KmerAnalyzer.CountKmersParallel(big, 6, default, 4, default, new SyncProgress(p => { lock (gate) reports.Add(p); }));

        Assert.That(reports, Is.Not.Empty);
        Assert.That(reports[^1], Is.EqualTo(1.0));
        Assert.That(reports.Count(p => p == 1.0), Is.EqualTo(1));
        Assert.That(reports.Take(reports.Count - 1), Is.All.InRange(0.0, 0.9999999));
        Assert.That(reports, Is.Ordered);
    }

    private sealed class SyncProgress(Action<double> onReport) : IProgress<double>
    {
        public void Report(double value) => onReport(value);
    }

    #endregion

    #region D2* / D2S — reference values (Python replica of Reinert 2009 / Song 2014 formulas)

    // Columns: seqA, seqB, k, Markov order r, D2*, D2S, d2*, d2S.
    private static IEnumerable<TestCaseData> D2StarCases()
    {
        yield return new TestCaseData(S1, S2, 2, 0, 1.8888612188959537, -0.06301942672096741, 0.43174180477939367, 0.5016118932249416).SetName("D2Star_S1S2_k2_r0");
        yield return new TestCaseData(S1, S2, 3, 0, 7.136981012975184, -0.6884013403313263, 0.44457941706964565, 0.5084031847096179).SetName("D2Star_S1S2_k3_r0");
        yield return new TestCaseData(S1, S2, 3, 1, -2.9337165550068374, -1.3932480549083965, 0.5426786050393649, 0.5226793773737308).SetName("D2Star_S1S2_k3_r1");
        yield return new TestCaseData(S1, S2, 4, 0, 7.194688991447289, 8.933007449531877, 0.4851244755594203, 0.44374461061164805).SetName("D2Star_S1S2_k4_r0");
        yield return new TestCaseData(S1, S2, 4, 2, 3.1149079908124557, 2.0271033958606988, 0.3759436961250997, 0.47599391302619776).SetName("D2Star_S1S2_k4_r2");
        yield return new TestCaseData(S1, S2, 5, 1, -45.67867820189822, 9.680437548004168, 0.5307449284905376, 0.45501071948433036).SetName("D2Star_S1S2_k5_r1");
        yield return new TestCaseData(S1, S3, 2, 0, -6.209385055236531, -8.829850477014112, 0.7186758460218481, 0.7130516334970357).SetName("D2Star_S1S3_k2_r0");
        yield return new TestCaseData(S1, S3, 3, 0, -8.11343648465283, -2.0818763897269443, 0.5672461669553764, 0.5253895649261623).SetName("D2Star_S1S3_k3_r0");
        yield return new TestCaseData(S1, S3, 3, 1, 2.2822994597132555, 2.038220365866943, 0.46132223213115664, 0.4643608357055673).SetName("D2Star_S1S3_k3_r1");
        yield return new TestCaseData(S1, S3, 4, 0, -9.96036383220566, 8.153419580674003, 0.5205984051117221, 0.44561507973368714).SetName("D2Star_S1S3_k4_r0");
        yield return new TestCaseData(S1, S3, 4, 2, 0.7666959645775582, 1.4472648349778436, 0.4855506868801013, 0.4836048522010693).SetName("D2Star_S1S3_k4_r2");
        yield return new TestCaseData(S1, S3, 5, 1, 27.25084993741092, 15.703707350850248, 0.47723568878974626, 0.41877711459909545).SetName("D2Star_S1S3_k5_r1");
        yield return new TestCaseData(S2, S3, 2, 0, -1.1112383756551796, -0.08368275725436958, 0.5337000828939916, 0.5020458542355434).SetName("D2Star_S2S3_k2_r0");
        yield return new TestCaseData(S2, S3, 3, 0, -10.762534477434523, -2.413979659201138, 0.5685895244149088, 0.5290950989157802).SetName("D2Star_S2S3_k3_r0");
        yield return new TestCaseData(S2, S3, 3, 1, -0.9958183635156561, -1.966785597693506, 0.5164776105446381, 0.5323555318676962).SetName("D2Star_S2S3_k3_r1");
        yield return new TestCaseData(S2, S3, 4, 0, -23.714332789584113, 7.29572604255774, 0.5434633294126657, 0.4491759905434676).SetName("D2Star_S2S3_k4_r0");
        yield return new TestCaseData(S2, S3, 4, 2, -5.390055862578997, -1.8730868082828067, 0.6341957867903316, 0.5252615021631372).SetName("D2Star_S2S3_k4_r2");
        yield return new TestCaseData(S2, S3, 5, 1, -27.466205254757217, 9.863830096809652, 0.5229072948662741, 0.4444945828242881).SetName("D2Star_S2S3_k5_r1");
    }

    [TestCaseSource(nameof(D2StarCases))]
    public void BackgroundAdjustedD2_MatchesReferenceReplica(string a, string b, int k, int r,
        double d2Star, double d2S, double d2StarDistance, double d2SDistance)
    {
        var stats = KmerAnalyzer.BackgroundAdjustedD2(a, b, k, r);
        Assert.Multiple(() =>
        {
            Assert.That(stats.D2Star, Is.EqualTo(d2Star).Within(Tol));
            Assert.That(stats.D2Shepherd, Is.EqualTo(d2S).Within(Tol));
            Assert.That(stats.D2StarDistance, Is.EqualTo(d2StarDistance).Within(Tol));
            Assert.That(stats.D2ShepherdDistance, Is.EqualTo(d2SDistance).Within(Tol));
        });
    }

    [Test]
    public void KmerDistance_D2StarMetrics_EqualOrderZeroStatistics_AndAreSymmetric()
    {
        var stats = KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 3);
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.KmerDistance(S1, S2, 3, KmerDistanceMetric.D2Star), Is.EqualTo(stats.D2StarDistance));
            Assert.That(KmerAnalyzer.KmerDistance(S1, S2, 3, KmerDistanceMetric.D2Shepherd), Is.EqualTo(stats.D2ShepherdDistance));
            Assert.That(KmerAnalyzer.KmerDistance(S2, S1, 3, KmerDistanceMetric.D2Star), Is.EqualTo(stats.D2StarDistance).Within(Tol));
            Assert.That(KmerAnalyzer.KmerDistance(S2, S1, 3, KmerDistanceMetric.D2Shepherd), Is.EqualTo(stats.D2ShepherdDistance).Within(Tol));
        });
    }

    [Test]
    public void BackgroundAdjustedD2_IdenticalSequences_DistanceZero_AndInUnitInterval()
    {
        var self = KmerAnalyzer.BackgroundAdjustedD2(S1, S1, 4, 1);
        Assert.That(self.D2StarDistance, Is.EqualTo(0.0).Within(Tol));
        Assert.That(self.D2ShepherdDistance, Is.EqualTo(0.0).Within(Tol));
        var other = KmerAnalyzer.BackgroundAdjustedD2(S1, S3, 4, 1);
        Assert.That(other.D2StarDistance, Is.InRange(0.0, 1.0));
        Assert.That(other.D2ShepherdDistance, Is.InRange(0.0, 1.0));
    }

    [Test]
    public void BackgroundAdjustedD2_CaseInsensitive_NonAcgtWindowsSkipped()
    {
        var mixed = KmerAnalyzer.BackgroundAdjustedD2(S3, S2, 3);
        var upper = KmerAnalyzer.BackgroundAdjustedD2(S3.ToUpperInvariant(), S2.ToLowerInvariant(), 3);
        Assert.That(upper, Is.EqualTo(mixed));
    }

    [Test]
    public void BackgroundAdjustedD2_SequenceMatchingItsBackground_DistanceIsNaN()
    {
        // Homopolymer, order 0: p(A) = 1 so E_X(AAA) = X(AAA) and every centred count is 0 -> 0/0.
        var stats = KmerAnalyzer.BackgroundAdjustedD2("AAAAAAAAAA", S1, 3);
        Assert.That(double.IsNaN(stats.D2StarDistance), Is.True);
        Assert.That(stats.D2Star, Is.EqualTo(0.0));
    }

    [Test]
    public void BackgroundAdjustedD2_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.BackgroundAdjustedD2(null!, S1, 3))!.ParamName, Is.EqualTo("seq1"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.BackgroundAdjustedD2(S1, null!, 3))!.ParamName, Is.EqualTo("seq2"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 0))!.ParamName, Is.EqualTo("k"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(
                () => KmerAnalyzer.BackgroundAdjustedD2(S1, S2, KmerAnalyzer.MaxBackgroundAdjustedK + 1))!.ParamName, Is.EqualTo("k"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 3, 3))!.ParamName, Is.EqualTo("markovOrder"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 3, -2))!.ParamName, Is.EqualTo("markovOrder"));
            Assert.That(Assert.Throws<ArgumentException>(() => KmerAnalyzer.BackgroundAdjustedD2("NNNNNN", S1, 3))!.ParamName, Is.EqualTo("seq1"));
            Assert.That(Assert.Throws<ArgumentException>(() => KmerAnalyzer.BackgroundAdjustedD2(S1, "AC", 3))!.ParamName, Is.EqualTo("seq2"));
            Assert.Throws<ArgumentException>(() => KmerAnalyzer.KmerDistance(
                KmerAnalyzer.CountKmers(S1, 3), KmerAnalyzer.CountKmers(S2, 3), KmerDistanceMetric.D2Star));
        });
    }

    // BIC(r) = -2 ln L_r + 3·4^r·ln N_r (Schwarz 1978; Katz 1981; CAFE -M -1), Python replica values.
    private const string Periodic = "ACGGTCATGCGATACCGTAGACGGTCATGCGATACCGTAGACGGTCATGCGATACCGTAGACGGTCATGCGATACCGTAGACGGTCATGCGTTACCGTAG";

    [TestCase(S1, new[] { 231.88984329117494, 255.74056400499072, 370.90681674887685, 914.1030179340269 })]
    [TestCase(S2, new[] { 199.825085479026, 217.12856241056033, 322.30742386881116, 879.0755936944657 })]
    [TestCase(S3, new[] { 202.97055742881076, 204.33298318566528, 301.40145051437, 802.2954904698054 })]
    [TestCase(Periodic, new[] { 289.4173873713581, 280.22454714161717, 324.38058560820787, 888.755266614427 })]
    public void MarkovOrderBic_MatchesReplica(string sequence, double[] expected)
    {
        for (int r = 0; r < expected.Length; r++)
            Assert.That(KmerAnalyzer.MarkovOrderBic(sequence, r), Is.EqualTo(expected[r]).Within(1e-9), $"r={r}");
    }

    [Test]
    public void BackgroundAdjustedD2_AutoOrder_UsesBicChoicePerSequence()
    {
        Assert.That(KmerAnalyzer.SelectMarkovOrder(Periodic, 3), Is.EqualTo(1));
        Assert.That(KmerAnalyzer.SelectMarkovOrder(S1, 3), Is.EqualTo(0));
        var auto = KmerAnalyzer.BackgroundAdjustedD2(Periodic, S1, 4, KmerAnalyzer.AutoMarkovOrder);
        Assert.Multiple(() =>
        {
            Assert.That((auto.MarkovOrder1, auto.MarkovOrder2), Is.EqualTo((1, 0)));
            Assert.That(auto.D2Star, Is.EqualTo(-2.128065438849073).Within(Tol));
            Assert.That(auto.D2Shepherd, Is.EqualTo(8.798084968717967).Within(Tol));
            Assert.That(auto.D2StarDistance, Is.EqualTo(0.5040674994481708).Within(Tol));
            Assert.That(auto.D2ShepherdDistance, Is.EqualTo(0.4586027673502613).Within(Tol));
        });
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 3, -2));
        Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.MarkovOrderBic(null!, 0));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.MarkovOrderBic("NNNN", 0));
    }

    [TestCase("euclidean", KmerDistanceMetric.Euclidean)]
    [TestCase(null, KmerDistanceMetric.Euclidean)]
    [TestCase(" D2STAR ", KmerDistanceMetric.D2Star)]
    [TestCase("d2shepherd", KmerDistanceMetric.D2Shepherd)]
    [TestCase("d2s", KmerDistanceMetric.D2Shepherd)]
    [TestCase("squared_euclidean_counts", KmerDistanceMetric.SquaredEuclideanCounts)]
    [TestCase("d2", KmerDistanceMetric.D2)]
    public void ParseDistanceMetric_KnownNames(string? name, KmerDistanceMetric expected)
        => Assert.That(KmerAnalyzer.ParseDistanceMetric(name), Is.EqualTo(expected));

    [Test]
    public void ParseDistanceMetric_Unknown_Throws()
        => Assert.Throws<ArgumentException>(() => KmerAnalyzer.ParseDistanceMetric("hamming"));

    #endregion
}
