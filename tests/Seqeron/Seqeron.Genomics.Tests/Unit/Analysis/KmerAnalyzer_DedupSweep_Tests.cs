// KMER-DIST-001 / KMER-COUNT-001 (B06 duplication sweep) — the helpers factored out of KmerAnalyzer.cs.
// Evidence: docs/Evidence/KMER-DIST-001-Evidence.md; TestSpec: KMER-DIST-001.md (rows X1–X4).
// Sources: Schwarz 1978 BIC / Katz 1981 Markov-order estimation (MarkovOrderBic); CAFE (Lu et al. 2017) MAX_ORDER = 10
//          and -M -1 order search over [0, min(k − 1, 10)].
// Reference values: the replica BIC values of S1/S2 already locked in KmerAnalyzer_ParallelAndBackgroundD2_Tests and
//          Seqeron.Mcp.Analysis.Tests/KmerD2StatisticsTests (independent Python replica, Schwarz 1978); NOT the code output.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class KmerAnalyzer_DedupSweep_Tests
{
    private const string S1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
    private const string S2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";

    // X1: MarkovOrderBics = MarkovOrderBic per order (replica values, Schwarz 1978).
    [Test]
    public void MarkovOrderBics_EqualsReplicaBicPerOrder()
    {
        Assert.Multiple(() =>
        {
            Assert.That(KmerAnalyzer.MarkovOrderBics(S1, 2),
                Is.EqualTo(new[] { 231.88984329117494, 255.74056400499072, 370.90681674887685 }).Within(1e-9));
            Assert.That(KmerAnalyzer.MarkovOrderBics(S2, 2),
                Is.EqualTo(new[] { 199.825085479026, 217.12856241056033, 322.30742386881116 }).Within(1e-9));
            for (int r = 0; r <= 4; r++)
                Assert.That(KmerAnalyzer.MarkovOrderBics(S1, 4)[r], Is.EqualTo(KmerAnalyzer.MarkovOrderBic(S1, r)), $"r={r}");
        });
    }

    // X2: SelectMarkovOrder is the first arg-min of MarkovOrderBics (ties → smaller order).
    [TestCase(S1, 3)]
    [TestCase(S2, 3)]
    [TestCase("ACACACACACACACACACACACACACACACACACACACAC", 3)]
    public void SelectMarkovOrder_IsFirstArgMinOfBics(string sequence, int maxOrder)
    {
        var bics = KmerAnalyzer.MarkovOrderBics(sequence, maxOrder);
        Assert.That(KmerAnalyzer.SelectMarkovOrder(sequence, maxOrder), Is.EqualTo(Array.IndexOf(bics, bics.Min())));
    }

    // X3: contract of MarkovOrderBics (same exceptions as SelectMarkovOrder / MarkovOrderBic).
    [Test]
    public void MarkovOrderBics_Validation()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.MarkovOrderBics(S1, -1));
            Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.MarkovOrderBics(null!, 1));
            Assert.Throws<ArgumentException>(() => KmerAnalyzer.MarkovOrderBics("NNNN", 0));
        });
    }

    // X4: AutoMarkovOrderLimit = min(k − 1, 10) (CAFE MAX_ORDER = 10); every KmerDistance overload rejects k ≤ 0
    // through the single check with the same exception and parameter name.
    [TestCase(1, 0)]
    [TestCase(5, 4)]
    [TestCase(11, 10)]
    [TestCase(12, 10)]
    public void AutoMarkovOrderLimit_IsMinOfKMinusOneAndTen(int k, int expected)
        => Assert.That(KmerAnalyzer.AutoMarkovOrderLimit(k), Is.EqualTo(expected));

    [TestCase(0)]
    [TestCase(-3)]
    public void KmerDistanceOverloads_NonPositiveK_SameException(int k)
    {
        var thrown = new[]
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.KmerDistance(S1, S2, k)),
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.KmerDistance(S1, S2, k, KmerDistanceMetric.Cosine)),
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.KmerDistance(S1, S2, k, KmerDistanceMetric.D2Star, 0)),
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.KmerDistance(S1, S2, k, KmerDistanceMetric.D2Star, 0, true)),
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.MashDistanceFromJaccard(0.5, k)),
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.JaccardSimilarity(S1, S2, k)),
            Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.GenerateAllKmers(k)),
        };
        Assert.That(thrown.Select(e => (e!.ParamName, e.Message)).Distinct().Count(), Is.EqualTo(1));
        Assert.That(thrown[0]!.ParamName, Is.EqualTo("k"));
    }
}
