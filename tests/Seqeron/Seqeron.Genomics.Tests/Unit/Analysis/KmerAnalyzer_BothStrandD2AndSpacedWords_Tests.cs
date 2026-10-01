// KMER-DIST-001 (audit round 2, WP6) — CAFE both-strand (-R) D2*/D2S, sparse Markov tables for high orders,
// null-as-empty for every KmerDistance metric, Jensen–Shannon / count-Euclidean metrics, multiple-pattern
// spaced-word distance.
// Evidence: docs/Evidence/KMER-DIST-001-Evidence.md; TestSpecs: KMER-DIST-001.md §4.6; algorithm doc §7.5.
// Sources: CAFE (Lu et al. 2017) kmer.cpp KmerModel::load (count of w += count of RC(w)) and
//          KmerProbEnsembDelegate::getKmerlogProb (p = ½(p(w) + p(RC(w)))); Leimeister, Boden, Horwege, Lindner &
//          Morgenstern 2014 (Bioinformatics 30:1991) multiple-pattern distance = average of per-pattern distances
//          (Euclidean / Jensen–Shannon of relative spaced-word frequencies); Lin 1991 (JS divergence).
// Reference values:
//   - D2*/D2S both strands: Python replica (d2rref.py), whose formula engine with CAFE's estimator reproduces the
//     CAFE -R binary (Jellyfish 2.3.1) on 20 runs to 6 digits; values below use the sequence MLE background.
//   - Spaced words: the spaced 1.2.0 binary (Ubuntu package, source sort.h spacedDNA) run with -r -f <patterns>
//     -d JS / -d EU (printed to 12 significant digits), the Python replica, and scipy 1.17.1
//     (jensenshannon(p, q, base=2)**2, euclidean on relative frequencies).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_BothStrandD2AndSpacedWords_Tests
{
    private const double Tol = 1e-12;

    private const string S1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
    private const string S2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
    private const string S3 = "acgtNNacgtacgRtTTGCAnAGGATTACAggattacaTTTAAACCCGGGTTANcaataccgagcccatggagc";

    // The 300/250-nt pair of the WP4 CAFE runs (A.fa / B.fa).
    private const string A = "CAGAGCAGACAACTAAGTGCTATCAACTAGGCGAAAGCCGCCTGAGGTGCTACTACAGTGTCGGGCTTCGAGGTTCCGACAAATAACTACGTTTCCCTTAAAACGCACCGTGGGGATTTTCCAGAAAACAAAACATGACCCATTCCAACCTAATGAGAGTTGCCATGTCATTTTTGAGCAACCGTCTTTCAAAAGTTCGTAGTTTCATCTTCCTGAAATTATTGCGAATGGTCTTACCAGCCATCCGTCTGGGACAATACGGCGGTAGACTGGTTCGGGGCGCTGTTCGTTAACAAAGTT";
    private const string B = "ATTATTATCGTTAGGCACTAGGACTGATTTAAATAACTTAATGTAATCATTTATATGCGTAAGAAAAACCTCGACAAATGAGTATGGTCGTTCTTTCCAAATAAATTTAACGAGATTGATCCACGGAGAAACAACATGTTTTCCTATTATTTTTAGGTTTTTTTAAACATGTTTGATTGGTATAAATATTGCGTTTTAAATCGAAATTTCGGGATATTAGTTGAATATAGTATGTTATGAAGGCACCTAT";

    // Python random.seed(7): concatenations of 12-mers drawn from a pool of 6, so 9-/10-mers repeat (order 8 is
    // not degenerate) — exercises the sparse Markov tables (r + 1 > 8).
    private const string R1 = "ACGTCAGCACGAAACTTGTTGGCCCTTAAGGGTTAACTTAAGGGTTAAGCTAAAGACAATGCTAAAGACAATCTTAAGGGTTAACTTAAGGGTTAAACGTCAGCACGACTTAAGGGTTAACAGTGTGAATCGCTTAAGGGTTAAAACTTGTTGGCCACGTCAGCACGACTTAAGGGTTAAAACTTGTTGGCCCTTAAGGGTTAAACGTCAGCACGAGCTAAAGACAATAACTTGTTGGCCACGTCAGCACGATACATAACATACCAGTGTGAATCGGCTAAAGACAATAACTTGTTGGCCGCTAAAGACAATTACATAACATACACGTCAGCACGATACATAACATACCTTAAGGGTTAATACATAACATACAACTTGTTGGCCAACTTGTTGGCCAACTTGTTGGCCGCTAAAGACAATTACATAACATACAACTTGTTGGCCAACTTGTTGGCCCAGTGTGAATCGACGTCAGCACGA";
    private const string R2 = "TACATAACATACAACTTGTTGGCCCAGTGTGAATCGACGTCAGCACGACTTAAGGGTTAAAACTTGTTGGCCACGTCAGCACGACTTAAGGGTTAAAACTTGTTGGCCTACATAACATACTACATAACATACGCTAAAGACAATTACATAACATACTACATAACATACTACATAACATACCTTAAGGGTTAATACATAACATACGCTAAAGACAATAACTTGTTGGCCCAGTGTGAATCGTACATAACATACACGTCAGCACGAACGTCAGCACGAGCTAAAGACAATTACATAACATACAACTTGTTGGCCCAGTGTGAATCGACGTCAGCACGACAGTGTGAATCGCAGTGTGAATCGACGTCAGCACGATACATAACATACCTTAAGGGTTAACAGTGTGAATCGCAGTGTGAATCG";

    #region D2* / D2S both strands (CAFE -R)

    private static IEnumerable<TestCaseData> BothStrandCases()
    {
        yield return new TestCaseData(S1, S2, 2, 0, 8.301008844517101, 8.708952662259415, 0.23440687751350275, 0.3448243392474428).SetName("D2StarBoth_S1S2_k2_r0");
        yield return new TestCaseData(S1, S2, 3, 0, 15.34415418761105, 8.071125717307769, 0.3838876158581438, 0.4345248995381073).SetName("D2StarBoth_S1S2_k3_r0");
        yield return new TestCaseData(S1, S2, 3, 1, -5.014593950349383, -1.3145456673597726, 0.5716659331844273, 0.5152519131905425).SetName("D2StarBoth_S1S2_k3_r1");
        yield return new TestCaseData(S1, S2, 4, 2, 8.443961462020846, 4.871990517853244, 0.39926744985723295, 0.4625333848223152).SetName("D2StarBoth_S1S2_k4_r2");
        yield return new TestCaseData(S1, S2, 5, 1, -47.72646950483297, 17.04308080418049, 0.5332230620383042, 0.4541988628386579).SetName("D2StarBoth_S1S2_k5_r1");
        yield return new TestCaseData(S1, S3, 2, 0, -9.90968949306743, -15.47971096418623, 0.8207440289946288, 0.7613563478346214).SetName("D2StarBoth_S1S3_k2_r0");
        yield return new TestCaseData(S1, S3, 3, 0, -7.245296905311575, -1.6687474857188895, 0.5556090643283262, 0.5133693804472232).SetName("D2StarBoth_S1S3_k3_r0");
        yield return new TestCaseData(S1, S3, 3, 1, 11.251668423001119, 8.969021843202096, 0.34979466258044856, 0.3952242630040085).SetName("D2StarBoth_S1S3_k3_r1");
        yield return new TestCaseData(S1, S3, 4, 2, 12.003892875760247, 5.4142613404933, 0.3774242096096003, 0.46371667790637694).SetName("D2StarBoth_S1S3_k4_r2");
        yield return new TestCaseData(S1, S3, 5, 1, 71.47472465481589, 34.40727629863311, 0.4604083804025727, 0.40200466168738014).SetName("D2StarBoth_S1S3_k5_r1");
        yield return new TestCaseData(S2, S3, 2, 0, -2.2340882668278494, 2.7926077715892728, 0.5752990255796153, 0.45016226041946134).SetName("D2StarBoth_S2S3_k2_r0");
        yield return new TestCaseData(S2, S3, 3, 0, -26.764124787433442, -17.36707960546063, 0.6634709176745822, 0.6330691167598625).SetName("D2StarBoth_S2S3_k3_r0");
        yield return new TestCaseData(S2, S3, 3, 1, -14.13032895548443, -9.637486419279977, 0.6562875903831356, 0.5994003559542167).SetName("D2StarBoth_S2S3_k3_r1");
        yield return new TestCaseData(S2, S3, 4, 2, 1.7539153086494077, 2.356110650153855, 0.47862802572142515, 0.48036792060212286).SetName("D2StarBoth_S2S3_k4_r2");
        yield return new TestCaseData(S2, S3, 5, 1, -24.375446739666646, 22.353229824080042, 0.5168739781345639, 0.43031257820937474).SetName("D2StarBoth_S2S3_k5_r1");
        yield return new TestCaseData(A, B, 3, 0, 18.903110077470167, 22.49160643946854, 0.35164895823974857, 0.39369557028358376).SetName("D2StarBoth_AB_k3_r0");
        yield return new TestCaseData(A, B, 4, 2, 0.0831585146242248, -9.136104139262747, 0.49967091530489954, 0.5285428590309681).SetName("D2StarBoth_AB_k4_r2");
        yield return new TestCaseData(A, B, 5, 1, 65.83233571592478, 31.770607381468807, 0.4656307641321891, 0.46046441660595405).SetName("D2StarBoth_AB_k5_r1");
    }

    [TestCaseSource(nameof(BothStrandCases))]
    public void BackgroundAdjustedD2_BothStrands_MatchesReplica(string a, string b, int k, int r,
        double d2Star, double d2S, double d2StarDistance, double d2SDistance)
    {
        var stats = KmerAnalyzer.BackgroundAdjustedD2(a, b, k, r, bothStrands: true);
        Assert.Multiple(() =>
        {
            Assert.That(stats.D2Star, Is.EqualTo(d2Star).Within(Tol));
            Assert.That(stats.D2Shepherd, Is.EqualTo(d2S).Within(Tol));
            Assert.That(stats.D2StarDistance, Is.EqualTo(d2StarDistance).Within(Tol));
            Assert.That(stats.D2ShepherdDistance, Is.EqualTo(d2SDistance).Within(Tol));
        });
    }

    [Test]
    public void BackgroundAdjustedD2_BothStrandsFalse_IsTheSingleStrandStatistic()
    {
        // WP4 locked value (S1/S2 k=3 r=0) through the new overload and the 4-argument one.
        var single = KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 3, 0, bothStrands: false);
        Assert.Multiple(() =>
        {
            Assert.That(single.D2StarDistance, Is.EqualTo(0.44457941706964565).Within(Tol));
            Assert.That(single.D2ShepherdDistance, Is.EqualTo(0.5084031847096179).Within(Tol));
            Assert.That(KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 3, 0), Is.EqualTo(single));
        });
    }

    [Test]
    public void BackgroundAdjustedD2_BothStrands_InvariantUnderReverseComplementOfEitherSequence()
    {
        // X(w) + X(RC w) is the same table for a sequence and its reverse complement; with r = 0 the background of
        // RC(s) is the complemented letter distribution, so p(w) + p(RC w) is unchanged too.
        string rc1 = new DnaSequence(S1).ReverseComplement().Sequence;
        var direct = KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 4, 0, bothStrands: true);
        var flipped = KmerAnalyzer.BackgroundAdjustedD2(rc1, S2, 4, 0, bothStrands: true);
        Assert.Multiple(() =>
        {
            Assert.That(flipped.D2Star, Is.EqualTo(direct.D2Star).Within(1e-9));
            Assert.That(flipped.D2Shepherd, Is.EqualTo(direct.D2Shepherd).Within(1e-9));
            Assert.That(flipped.D2StarDistance, Is.EqualTo(direct.D2StarDistance).Within(Tol));
        });
    }

    [Test]
    public void KmerDistance_BothStrands_EqualsStatistics_AndAutoOrderReportsBicChoice()
    {
        var stats = KmerAnalyzer.BackgroundAdjustedD2(S1, S3, 4, 2, bothStrands: true);
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.KmerDistance(S1, S3, 4, KmerDistanceMetric.D2Star, 2, bothStrands: true), Is.EqualTo(stats.D2StarDistance));
            Assert.That(KmerAnalyzer.KmerDistance(S1, S3, 4, KmerDistanceMetric.D2Shepherd, 2, bothStrands: true), Is.EqualTo(stats.D2ShepherdDistance));
            Assert.That(KmerAnalyzer.KmerDistance(S1, S3, 4, KmerDistanceMetric.D2Star, 2), Is.EqualTo(KmerAnalyzer.BackgroundAdjustedD2(S1, S3, 4, 2).D2StarDistance));
        });

        var auto = KmerAnalyzer.BackgroundAdjustedD2(S1, S2, 4, KmerAnalyzer.AutoMarkovOrder, bothStrands: true);
        Assert.Multiple(() =>
        {
            Assert.That(auto.MarkovOrder1, Is.EqualTo(KmerAnalyzer.SelectMarkovOrder(S1, 3)));
            Assert.That(auto.MarkovOrder2, Is.EqualTo(KmerAnalyzer.SelectMarkovOrder(S2, 3)));
        });
    }

    [Test]
    public void KmerDistance_BothStrands_OnPlainMetric_Throws()
    {
        var ex = Assert.Throws<ArgumentException>(() => KmerAnalyzer.KmerDistance(S1, S2, 3, KmerDistanceMetric.Euclidean, 0, bothStrands: true));
        Assert.That(ex!.ParamName, Is.EqualTo("bothStrands"));
        Assert.That(KmerAnalyzer.KmerDistance(S1, S2, 3, KmerDistanceMetric.Euclidean, 0, bothStrands: false),
            Is.EqualTo(KmerAnalyzer.KmerDistance(S1, S2, 3)));
    }

    #endregion

    #region High Markov orders (sparse tables) and null-as-empty

    [TestCase(false, -0.0044688633703816516, 0.8149552746135663, 0.510290283666086, 0.3826552653739205)]
    [TestCase(true, -0.008937726740763303, 1.6299105492271277, 0.510290283666086, 0.3826552653739207)]
    public void BackgroundAdjustedD2_Order8_SparseTables_MatchReplica(bool bothStrands,
        double d2Star, double d2S, double d2StarDistance, double d2SDistance)
    {
        var stats = KmerAnalyzer.BackgroundAdjustedD2(R1, R2, 10, 8, bothStrands);
        Assert.Multiple(() =>
        {
            Assert.That(stats.D2Star, Is.EqualTo(d2Star).Within(1e-9));
            Assert.That(stats.D2Shepherd, Is.EqualTo(d2S).Within(1e-9));
            Assert.That(stats.D2StarDistance, Is.EqualTo(d2StarDistance).Within(1e-9));
            Assert.That(stats.D2ShepherdDistance, Is.EqualTo(d2SDistance).Within(1e-9));
        });
    }

    [Test]
    public void BackgroundAdjustedD2_K12Order11_NoDenseFourPowerKTable()
    {
        // Before: a 4^12-double transition array (134 MB) per sequence. Now the Markov tables hold only the r-mers and
        // (r+1)-mers present; the allocation is dominated by the k-mer count tables of the two short sequences.
        long before = GC.GetAllocatedBytesForCurrentThread();
        var stats = KmerAnalyzer.BackgroundAdjustedD2(R1, R2, 12, 11);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.Multiple(() =>
        {
            Assert.That(allocated, Is.LessThan(16L * 1024 * 1024));
            Assert.That(stats.MarkovOrder1, Is.EqualTo(11));
        });
    }

    [Test]
    public void KmerDistance_NullSequence_TreatedAsEmpty_ForEveryMetric()
    {
        // Word-vector metrics: null = "" = the zero vector.
        Assert.That(KmerAnalyzer.KmerDistance(null!, S1, 3, KmerDistanceMetric.Manhattan, 0),
            Is.EqualTo(KmerAnalyzer.KmerDistance("", S1, 3, KmerDistanceMetric.Manhattan, 0)));

        // D2*/D2S: null is rejected exactly like the empty sequence (no ACGT k-mer => no background), not as a null argument.
        foreach (var metric in new[] { KmerDistanceMetric.D2Star, KmerDistanceMetric.D2Shepherd })
        {
            var nullEx = Assert.Throws<ArgumentException>(() => KmerAnalyzer.KmerDistance(null!, S1, 3, metric));
            var emptyEx = Assert.Throws<ArgumentException>(() => KmerAnalyzer.KmerDistance("", S1, 3, metric));
            var nullEx2 = Assert.Throws<ArgumentException>(() => KmerAnalyzer.KmerDistance(S1, null!, 3, metric, 0, true));
            Assert.Multiple(() =>
            {
                Assert.That(nullEx, Is.Not.InstanceOf<ArgumentNullException>());
                Assert.That(nullEx!.ParamName, Is.EqualTo(emptyEx!.ParamName).And.EqualTo("seq1"));
                Assert.That(nullEx2!.ParamName, Is.EqualTo("seq2"));
            });
        }
    }

    #endregion

    #region Jensen–Shannon and count-Euclidean metrics

    [Test]
    public void JensenShannon_MatchesScipy_AndConventions()
    {
        // ATGTGTG / CATGTG k=3: ATG 1, TGT 2, GTG 2 vs CAT 1, ATG 1, TGT 1, GTG 1.
        // scipy.spatial.distance.jensenshannon([0,1,2,2],[1,1,1,1], base=2)**2 = 0.1522040934665307.
        var a = KmerAnalyzer.CountKmers("ATGTGTG", 3);
        var b = KmerAnalyzer.CountKmers("CATGTG", 3);
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.KmerDistance(a, b, KmerDistanceMetric.JensenShannon), Is.EqualTo(0.1522040934665307).Within(Tol));
            Assert.That(KmerAnalyzer.KmerDistance(b, a, KmerDistanceMetric.JensenShannon), Is.EqualTo(0.1522040934665307).Within(Tol));
            Assert.That(KmerAnalyzer.KmerDistance(a, a, KmerDistanceMetric.JensenShannon), Is.EqualTo(0.0));
            // Disjoint supports: JS = 1 (maximum, base 2).
            Assert.That(KmerAnalyzer.KmerDistance(KmerAnalyzer.CountKmers("AAAA", 2), KmerAnalyzer.CountKmers("CCCC", 2), KmerDistanceMetric.JensenShannon),
                Is.EqualTo(1.0).Within(Tol));
            // Zero vector against a distribution: ½ (documented convention); both empty: 0.
            Assert.That(KmerAnalyzer.KmerDistance(new Dictionary<string, int>(), a, KmerDistanceMetric.JensenShannon), Is.EqualTo(0.5).Within(Tol));
            Assert.That(KmerAnalyzer.KmerDistance(new Dictionary<string, int>(), new Dictionary<string, int>(), KmerDistanceMetric.JensenShannon), Is.EqualTo(0.0));
            // Count Euclidean = sqrt of the squared count Euclidean (WP2 value 3 => sqrt 3).
            Assert.That(KmerAnalyzer.KmerDistance(a, b, KmerDistanceMetric.EuclideanCounts), Is.EqualTo(Math.Sqrt(3)).Within(Tol));
        });
    }

    [TestCase("jensen_shannon", KmerDistanceMetric.JensenShannon)]
    [TestCase(" JS ", KmerDistanceMetric.JensenShannon)]
    [TestCase("euclidean_counts", KmerDistanceMetric.EuclideanCounts)]
    public void ParseDistanceMetric_NewNames(string name, KmerDistanceMetric expected)
        => Assert.That(KmerAnalyzer.ParseDistanceMetric(name), Is.EqualTo(expected));

    #endregion

    #region Multiple-pattern spaced-word distance (Leimeister et al. 2014; spaced 1.2.0)

    private static readonly string[] PatternsW4 = { "11011", "10111", "11101" };
    private static readonly string[] PatternsW5 = { "1101011", "1011101", "1110011" };
    private static readonly string[] Contiguous4 = { "1111" };

    // spaced -r -f <patterns> -d JS / -d EU (12 significant digits) = Python replica (shown) = scipy.
    private static IEnumerable<TestCaseData> SpacedCases()
    {
        yield return new TestCaseData(S1, S2, PatternsW4, 0.8163228541607376, 12.40897581662776, 0.17567404832368613).SetName("Spaced_S1S2_w4");
        yield return new TestCaseData(S1, S2, PatternsW5, 0.9467819156191134, 11.660417980553781, 0.16946832495600062).SetName("Spaced_S1S2_w5");
        yield return new TestCaseData(S1, S2, Contiguous4, 0.8003989533664195, 12.489995996796797, 0.17553282636413806).SetName("Spaced_S1S2_contiguous4");
        yield return new TestCaseData(A, B, PatternsW4, 0.41451725541237366, 25.704410975393987, 0.09756721735575496).SetName("Spaced_AB_w4");
        yield return new TestCaseData(A, B, PatternsW5, 0.7579794284965499, 23.845599914433468, 0.0898354125393724).SetName("Spaced_AB_w5");
        yield return new TestCaseData(A, B, Contiguous4, 0.4189421165247839, 25.80697580112788, 0.09773005238795092).SetName("Spaced_AB_contiguous4");
    }

    [TestCaseSource(nameof(SpacedCases))]
    public void SpacedWordDistance_MatchesSpacedProgram_AndPaperEuclidean(string a, string b, string[] patterns,
        double spacedJs, double spacedEu, double frequencyEuclidean)
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.SpacedWordDistance(a, b, patterns, KmerDistanceMetric.JensenShannon), Is.EqualTo(spacedJs).Within(Tol));
            Assert.That(KmerAnalyzer.SpacedWordDistance(a, b, patterns, KmerDistanceMetric.EuclideanCounts), Is.EqualTo(spacedEu).Within(1e-10));
            Assert.That(KmerAnalyzer.SpacedWordDistance(a, b, patterns), Is.EqualTo(frequencyEuclidean).Within(Tol));
        });
    }

    [Test]
    public void SpacedWordDistance_ContiguousPattern_EqualsKmerDistance_AndIsAverageOfPatterns()
    {
        Assert.That(KmerAnalyzer.SpacedWordDistance(A, B, Contiguous4), Is.EqualTo(KmerAnalyzer.KmerDistance(A, B, 4)));
        double mean = PatternsW4.Average(p => KmerAnalyzer.KmerDistance(
            KmerAnalyzer.CountSpacedWords(S1, p), KmerAnalyzer.CountSpacedWords(S3, p), KmerDistanceMetric.Manhattan));
        Assert.That(KmerAnalyzer.SpacedWordDistance(S1, S3, PatternsW4, KmerDistanceMetric.Manhattan), Is.EqualTo(mean).Within(Tol));
        Assert.That(KmerAnalyzer.SpacedWordDistance(S1, S1, PatternsW5), Is.EqualTo(0.0));
    }

    [Test]
    public void SpacedWordDistance_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, null!))!.ParamName, Is.EqualTo("patterns"));
            Assert.That(Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, Array.Empty<string>()))!.ParamName, Is.EqualTo("patterns"));
            Assert.That(Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, new[] { "11", null! }))!.ParamName, Is.EqualTo("patterns"));
            Assert.That(Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, new[] { "1101", "11" }))!.ParamName, Is.EqualTo("patterns"));
            Assert.That(Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, new[] { "0110" }))!.ParamName, Is.EqualTo("pattern"));
            Assert.That(Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance(S1, S2, PatternsW4, KmerDistanceMetric.D2Star))!.ParamName, Is.EqualTo("metric"));
        });
    }

    #endregion
}
