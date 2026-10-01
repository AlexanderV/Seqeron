// KMER-DIST-001 (audit round 1, WP2) — k-mer distance metrics, exact k-mer Jaccard / Mash distance, spaced words
// Evidence: docs/Evidence/KMER-DIST-001-Evidence.md
// TestSpec: tests/TestSpecs/KMER-DIST-001.md
// Sources: Blaisdell 1986 (PNAS 83:5155) d_E; Vinga & Almeida 2003 (Bioinformatics 19:513); Torney et al. 1990 /
//          Reinert et al. 2009 (D2); Jaccard 1901/1912; Ondov et al. 2016 (Genome Biol 17:132, Mash eq. 1/4;
//          Mash 2.x CommandDistance.cpp / Sketch.cpp); Leimeister et al. 2014 (Bioinformatics 30:1991, spaced words).
// Reference values: scipy 1.17.1 spatial.distance on scikit-bio 0.7.4 kmer_frequencies; alfpy 1.0.6 word_distance;
//          Python set Jaccard; sourmash 4.9.4 MinHash(scaled=1); the real Mash 2.3 binary (mash dist -s 100000 [-n]).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_DistanceMetrics_Tests
{
    private const double Tol = 1e-12;

    // Deterministic 120-nt pair (Python random.seed(2026); every 9th base of R1 substituted A>C>G>T>A in R2).
    private const string R1 = "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCG";
    private const string R2 = "CGACTTTCACAGATATGCAGGGTAGAGTTCGAGGTTCTTATTTGTAACCAATTCACATTGTGTATCGGAACTAGCGTTTTATGTATGTCTAAGTGACTCAAAATACCACGGCAGTCCACG";

    private static readonly KmerCountingOptions Canonical = new(Canonical: true);
    private static readonly KmerCountingOptions AcgtOnly = new(AcgtOnly: true);

    #region Metrics — reference values (scipy / alfpy)

    // Columns: Euclidean (freq), SquaredEuclideanCounts, Manhattan (freq), Chebyshev (freq), Canberra (freq), Cosine (counts), D2 (counts).
    private static IEnumerable<TestCaseData> MetricCases()
    {
        yield return new TestCaseData("ATGTGTG", "CATGTG", 3,
            new[] { 0.33166247903554, 3.0, 0.6, 0.25, 1.5726495726495728, 0.16666666666666663, 5.0 })
            .SetName("Metrics_Zielezinski2017Fig1_k3");
        yield return new TestCaseData("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2,
            new[] { 0.26600633232367216, 11.0, 0.7692307692307692, 0.15384615384615385, 7.568421052631578, 0.29704050843336227, 13.0 })
            .SetName("Metrics_ACGTTGCAACGGT_k2");
        yield return new TestCaseData("GATTACAGATTACA", "GATTACCGATTTCA", 3,
            new[] { 0.31180478223116176, 14.0, 1.0, 0.16666666666666666, 9.666666666666666, 0.3603978509331687, 12.0 })
            .SetName("Metrics_GATTACA_k3");
        yield return new TestCaseData("acgtNNacgtacgRtTTGCAnA", "ACGTACGTTTGCAAA", 3,
            new[] { 0.2179788817739193, 13.0, 0.9, 0.07692307692307693, 13.49750671269659, 0.26664120237743094, 16.0 })
            .SetName("Metrics_LowerCaseAndIupac_Literal_k3");
        yield return new TestCaseData(R1, R2, 5,
            new[] { 0.09829098492233948, 130.0, 1.0517241379310347, 0.017241379310344827, 114.66666666666666, 0.5199846392626807, 60.0 })
            .SetName("Metrics_R1R2_k5");
        yield return new TestCaseData(R1, R2, 11,
            new[] { 0.13483997249264842, 220.0, 2.0, 0.00909090909090909, 220.0, 1.0, 0.0 })
            .SetName("Metrics_R1R2_k11_NoSharedKmer");
    }

    private static readonly KmerDistanceMetric[] MetricOrder =
    {
        KmerDistanceMetric.Euclidean, KmerDistanceMetric.SquaredEuclideanCounts, KmerDistanceMetric.Manhattan,
        KmerDistanceMetric.Chebyshev, KmerDistanceMetric.Canberra, KmerDistanceMetric.Cosine, KmerDistanceMetric.D2,
    };

    [TestCaseSource(nameof(MetricCases))]
    public void KmerDistance_Metric_MatchesScipyAndAlfpy(string a, string b, int k, double[] expected)
    {
        Assert.Multiple(() =>
        {
            for (int i = 0; i < MetricOrder.Length; i++)
            {
                Assert.That(KmerAnalyzer.KmerDistance(a, b, k, MetricOrder[i]), Is.EqualTo(expected[i]).Within(Tol), MetricOrder[i].ToString());
                Assert.That(KmerAnalyzer.KmerDistance(b, a, k, MetricOrder[i]), Is.EqualTo(expected[i]).Within(Tol), "symmetry " + MetricOrder[i]);
            }
        });
    }

    [Test]
    public void SquaredEuclideanCounts_Fig1_IsBlaisdellValue3()
    {
        // Vinga & Almeida 2003 / Blaisdell 1986 d_E on the Zielezinski 2017 Fig. 1 counts (1,0,2,2) vs (1,1,1,1): 0+1+1+1 = 3.
        Assert.That(KmerAnalyzer.KmerDistance("ATGTGTG", "CATGTG", 3, KmerDistanceMetric.SquaredEuclideanCounts), Is.EqualTo(3.0));
    }

    [TestCase("ATGTGTG", "CATGTG", 3)]
    [TestCase("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2)]
    [TestCase(R1, R2, 5)]
    [TestCase("", "ACGT", 2)]
    public void EuclideanMetric_IsBitIdenticalToLegacyKmerDistance(string a, string b, int k)
    {
        Assert.That(KmerAnalyzer.KmerDistance(a, b, k, KmerDistanceMetric.Euclidean), Is.EqualTo(KmerAnalyzer.KmerDistance(a, b, k)));
    }

    [Test]
    public void Metrics_IdenticalSequences_DistanceZero_D2IsSelfInnerProduct()
    {
        foreach (var m in MetricOrder.Where(m => m != KmerDistanceMetric.D2))
            Assert.That(KmerAnalyzer.KmerDistance(R1, R1, 4, m), Is.EqualTo(0.0).Within(Tol), m.ToString());

        // D2(x, x) = Σ c(w)² : ATGTGTG k=3 counts ATG 1, TGT 2, GTG 2 → 1 + 4 + 4 = 9.
        Assert.That(KmerAnalyzer.KmerDistance("ATGTGTG", "ATGTGTG", 3, KmerDistanceMetric.D2), Is.EqualTo(9.0));
    }

    [Test]
    public void Metrics_EmptyVectors_Conventions()
    {
        Assert.Multiple(() =>
        {
            // Both zero vectors: every additive metric is 0; cosine similarity of a zero vector is 0 → distance 1.
            foreach (var m in MetricOrder.Where(m => m != KmerDistanceMetric.Cosine))
                Assert.That(KmerAnalyzer.KmerDistance("", "AC", 3, m), Is.EqualTo(0.0), m.ToString());
            Assert.That(KmerAnalyzer.KmerDistance("", "AC", 3, KmerDistanceMetric.Cosine), Is.EqualTo(1.0));
            Assert.That(KmerAnalyzer.KmerDistance(null!, "ACGT", 2, KmerDistanceMetric.Cosine), Is.EqualTo(1.0));
            // One zero vector: Manhattan = Σ f = 1, Canberra = number of words, SquaredEuclideanCounts = Σ c².
            Assert.That(KmerAnalyzer.KmerDistance("", "ACGTA", 2, KmerDistanceMetric.Manhattan), Is.EqualTo(1.0).Within(Tol));
            Assert.That(KmerAnalyzer.KmerDistance("", "ACGTA", 2, KmerDistanceMetric.Canberra), Is.EqualTo(4.0));
            Assert.That(KmerAnalyzer.KmerDistance("", "AAAA", 2, KmerDistanceMetric.SquaredEuclideanCounts), Is.EqualTo(9.0));
        });
    }

    [Test]
    public void KmerDistance_Contracts()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.KmerDistance("ACGT", "ACGT", 0, KmerDistanceMetric.Manhattan));
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.KmerDistance("ACGT", "ACGT", 2, (KmerDistanceMetric)99));
        Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.KmerDistance(null!, new Dictionary<string, int>(), KmerDistanceMetric.D2));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.KmerDistance(
            new Dictionary<string, int> { ["A"] = -1 }, new Dictionary<string, int>(), KmerDistanceMetric.D2));
    }

    [Test]
    public void KmerDistance_OnCanonicalCountTables_UsesTheGivenTables()
    {
        // ATGTGTG vs CATGTG k=3 canonical (-C) tables: {ACA:2, ATG:1, CAC:2} vs {ACA:1, ATG:2, CAC:1} (CAT→ATG).
        var c1 = KmerAnalyzer.CountKmers("ATGTGTG", 3, Canonical);
        var c2 = KmerAnalyzer.CountKmers("CATGTG", 3, Canonical);
        Assert.That(KmerAnalyzer.KmerDistance(c1, c2, KmerDistanceMetric.SquaredEuclideanCounts), Is.EqualTo(3.0));
        Assert.That(KmerAnalyzer.KmerDistance(c1, c2, KmerDistanceMetric.D2), Is.EqualTo(2 + 2 + 2));
    }

    #endregion

    #region Jaccard / Mash — reference values (Python sets, sourmash scaled=1, Mash 2.3 binary)

    // Columns: literal J, ACGT-only J, canonical J, Mash D canonical, Mash D ACGT-only (non-canonical, mash -n).
    private static IEnumerable<TestCaseData> JaccardCases()
    {
        yield return new TestCaseData("ATGTGTG", "CATGTG", 3,
            new[] { 0.75, 0.75, 1.0, 0.0, 0.05138355994241945 }).SetName("Jaccard_Fig1_k3");
        yield return new TestCaseData("ACGTTGCAACGGT", "ACGTAGCATCGGTA", 2,
            new[] { 0.46153846153846156, 0.46153846153846156, 0.5, 0.20273255405408222, 0.2297661646892201 }).SetName("Jaccard_ACGTTGCAACGGT_k2");
        yield return new TestCaseData("GATTACAGATTACA", "GATTACCGATTTCA", 3,
            new[] { 0.3076923076923077, 0.3076923076923077, 0.3076923076923077, 0.2512572674587934, 0.2512572674587934 }).SetName("Jaccard_GATTACA_k3");
        yield return new TestCaseData("acgtNNacgtacgRtTTGCAnA", "ACGTACGTTTGCAAA", 3,
            new[] { 0.4, 0.7272727272727273, 0.8333333333333334, 0.03177005993477496, 0.05728341897555305 }).SetName("Jaccard_LowerCaseAndIupac_k3");
        yield return new TestCaseData(R1, R2, 5,
            new[] { 0.3273809523809524, 0.3273809523809524, 0.3670886075949367, 0.12433764331556005, 0.1413382811335405 }).SetName("Jaccard_R1R2_k5");
        yield return new TestCaseData(R1, R2, 11,
            new[] { 0.0, 0.0, 0.0, 1.0, 1.0 }).SetName("Jaccard_R1R2_k11_NoSharedKmer");
    }

    [TestCaseSource(nameof(JaccardCases))]
    public void JaccardAndMash_MatchReferences(string a, string b, int k, double[] e)
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.JaccardSimilarity(a, b, k), Is.EqualTo(e[0]).Within(Tol), "literal");
            Assert.That(KmerAnalyzer.JaccardSimilarity(a, b, k, AcgtOnly), Is.EqualTo(e[1]).Within(Tol), "acgt-only");
            Assert.That(KmerAnalyzer.JaccardSimilarity(a, b, k, Canonical), Is.EqualTo(e[2]).Within(Tol), "canonical (Mash/sourmash)");
            Assert.That(KmerAnalyzer.JaccardSimilarity(b, a, k, Canonical), Is.EqualTo(e[2]).Within(Tol), "symmetry");
            Assert.That(KmerAnalyzer.MashDistance(a, b, k, Canonical), Is.EqualTo(e[3]).Within(Tol), "mash dist");
            Assert.That(KmerAnalyzer.MashDistance(a, b, k, AcgtOnly), Is.EqualTo(e[4]).Within(Tol), "mash dist -n");
        });
    }

    [Test]
    public void MashBinaryReportedSharedHashes_EqualExactSetCounts()
    {
        // mash dist -k 5 -s 100000 R1 R2 → 0.124338, shared-hashes 58/158; -n → 0.141338, 55/168.
        Assert.That(KmerAnalyzer.JaccardSimilarity(R1, R2, 5, Canonical), Is.EqualTo(58.0 / 158));
        Assert.That(KmerAnalyzer.JaccardSimilarity(R1, R2, 5, AcgtOnly), Is.EqualTo(55.0 / 168));
        Assert.That(KmerAnalyzer.MashDistance(R1, R2, 5, Canonical), Is.EqualTo(0.124338).Within(5e-7));
        Assert.That(KmerAnalyzer.MashDistance(R1, R2, 5, AcgtOnly), Is.EqualTo(0.141338).Within(5e-7));
    }

    [Test]
    public void Jaccard_Conventions()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.JaccardSimilarity("", "", 3), Is.EqualTo(0.0), "both empty → 0 (GenomicAnalyzer convention)");
            Assert.That(KmerAnalyzer.JaccardSimilarity(null!, "AC", 3), Is.EqualTo(0.0));
            Assert.That(KmerAnalyzer.JaccardSimilarity("", "ACGT", 2), Is.EqualTo(0.0));
            Assert.That(KmerAnalyzer.JaccardSimilarity("acgtac", "ACGTAC", 3), Is.EqualTo(1.0), "case-insensitive");
            // Mash CommandDistance: common == denom (incl. 0 == 0) → 0.
            Assert.That(KmerAnalyzer.MashDistance("", "", 3, Canonical), Is.EqualTo(0.0));
            Assert.That(KmerAnalyzer.MashDistance("NNNN", "ACGT", 3, Canonical), Is.EqualTo(1.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.JaccardSimilarity("A", "A", 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.JaccardSimilarity("", "", 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.MashDistance("A", "A", -1, Canonical));
        });
    }

    [Test]
    public void Jaccard_EqualsGenomicAnalyzerFractionOnUpperCaseDna()
    {
        // DUP_MAP §16: the canonical must reproduce GenomicAnalyzer.CalculateSimilarity (percent) / 100.
        foreach (var k in new[] { 2, 3, 5, 8 })
        {
            double expected = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(R1), new DnaSequence(R2), k) / 100.0;
            Assert.That(KmerAnalyzer.JaccardSimilarity(R1, R2, k), Is.EqualTo(expected).Within(Tol), $"k={k}");
        }
    }

    [TestCase(1.0, 21, 0.0)]
    [TestCase(0.0, 21, 1.0)]
    [TestCase(0.5, 2, 0.20273255405408222)]     // −ln(2/3)/2
    [TestCase(0.75, 3, 0.05138355994241945)]    // −ln(6/7)/3
    [TestCase(1e-6, 1, 1.0)]                    // −ln(2e-6/(1+1e-6)) = 13.1 → capped at 1 (Mash)
    public void MashDistanceFromJaccard_OndovEq4WithMashBoundaries(double j, int k, double expected)
    {
        Assert.That(KmerAnalyzer.MashDistanceFromJaccard(j, k), Is.EqualTo(expected).Within(Tol));
    }

    [TestCase(-0.1)]
    [TestCase(1.1)]
    [TestCase(double.NaN)]
    public void MashDistanceFromJaccard_OutOfRange_Throws(double j)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.MashDistanceFromJaccard(j, 3));
    }

    #endregion

    #region Spaced words (Leimeister et al. 2014) — Python replica values

    [Test]
    public void CountSpacedWords_ReferenceTables()
    {
        Assert.That(KmerAnalyzer.CountSpacedWords("ATGTGTG", "101"),
            Is.EquivalentTo(new Dictionary<string, int> { ["AG"] = 1, ["GG"] = 2, ["TT"] = 2 }));
        Assert.That(KmerAnalyzer.CountSpacedWords("ACGTTGCAACGGT", "1101"),
            Is.EquivalentTo(new Dictionary<string, int>
            {
                ["AAG"] = 1, ["ACG"] = 1, ["ACT"] = 1, ["CAC"] = 1, ["CGT"] = 2, ["GCA"] = 1, ["GTG"] = 1, ["TGA"] = 1, ["TTC"] = 1,
            }));
        Assert.That(KmerAnalyzer.CountSpacedWords("gattacaGATTACA", "11011"),
            Is.EquivalentTo(new Dictionary<string, int>
            {
                ["ACGA"] = 1, ["AGTT"] = 1, ["ATAC"] = 2, ["CAAT"] = 1, ["GATA"] = 2, ["TAAG"] = 1, ["TTCA"] = 2,
            }));
    }

    [Test]
    public void CountSpacedWords_R1_Pattern1100111_DistinctAndTotal()
    {
        var t = KmerAnalyzer.CountSpacedWords(R1, "1100111");
        Assert.That(t.Count, Is.EqualTo(108));
        Assert.That(t.Values.Sum(), Is.EqualTo(R1.Length - 7 + 1)); // 114 windows
    }

    [Test]
    public void CountSpacedWords_AllOnesPattern_EqualsContiguousKmers()
    {
        Assert.That(KmerAnalyzer.CountSpacedWords(R1, "11111"), Is.EquivalentTo(KmerAnalyzer.CountKmers(R1, 5)));
        Assert.That(KmerAnalyzer.CountSpacedWords("ATGTGTG", "111"),
            Is.EquivalentTo(new Dictionary<string, int> { ["ATG"] = 1, ["GTG"] = 2, ["TGT"] = 2 }));
    }

    [Test]
    public void SpacedWordDistance_ViaCountTableOverload()
    {
        // Pattern 1101: scipy euclidean on frequencies 0.3149183286488868; sqeuclidean on counts 12.
        var a = KmerAnalyzer.CountSpacedWords("GATTACAGATTACA", "1101");
        var b = KmerAnalyzer.CountSpacedWords("GATTACCGATTTCA", "1101");
        Assert.That(KmerAnalyzer.KmerDistance(a, b, KmerDistanceMetric.Euclidean), Is.EqualTo(0.3149183286488868).Within(Tol));
        Assert.That(KmerAnalyzer.KmerDistance(a, b, KmerDistanceMetric.SquaredEuclideanCounts), Is.EqualTo(12.0));
    }

    [TestCase("")]
    [TestCase("0110")]
    [TestCase("1100")]
    [TestCase("1021")]
    public void CountSpacedWords_InvalidPattern_Throws(string pattern)
    {
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.CountSpacedWords("ACGTACGT", pattern));
    }

    [Test]
    public void CountSpacedWords_EdgeCases()
    {
        Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.CountSpacedWords("ACGT", null!));
        Assert.That(KmerAnalyzer.CountSpacedWords(null!, "101"), Is.Empty);
        Assert.That(KmerAnalyzer.CountSpacedWords("AC", "101"), Is.Empty);
        Assert.That(KmerAnalyzer.CountSpacedWords("ACG", "101"), Is.EquivalentTo(new Dictionary<string, int> { ["AG"] = 1 }));
    }

    #endregion
}
