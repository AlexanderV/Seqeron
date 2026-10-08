// REP-APPROX-001 — TRF detection pipeline parity (B04 audit WP7, L7): apparent-size criterion, active distance
// ranges, best-period list for d > 250, narrow-band wraparound DP for patterns > 20.
// Evidence: docs/Evidence/REP-APPROX-001-Evidence.md (§WP7); TestSpec: tests/TestSpecs/REP-APPROX-001.md (D1..D9).
// Sources: Benson G (1999) Nucleic Acids Res 27(2):573-580; TRF 4.10.0 README ("Apparent Size Distribution",
//          "Narrow Band Alignment", "Multiple Reporting of Repeat at Different Pattern Sizes", What's New 4.04/4.07b).
//
// Every expected row was produced by the compiled TRF 4.10.0 binary (github.com/Benson-Genomics-Lab/TRF, commit
// 355c1f9): `trf case.fa <Match> <Mismatch> <Delta> <PM> <PI> <Minscore> <MaxPeriod> -h -d -ngs`. Each sequence was
// selected (seeded random background + mutated tandem arrays, Evidence §WP7) because the previous implementation
// (commit ffe2267) disagreed with TRF on it and the disagreement disappears only with the component the test is
// named after (ablation builds, Evidence §WP7). TRF prints 1-based indices, copy numbers to one decimal and
// truncated percentages; the asserts compare against those printed values.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RepeatFinder_TrfDetection_Tests
{
    private const string Band31 =
        "CTAGGGTTGGGCTGGTGATCGGAGGACATTCCCCACTTCGATGTACTTTTGAGCTGGCCAGCAACCGCACTGGAGGTGTAATAGTCGTCTCTATTGCGGA"
        + "GGGTTAAGATGGCCCGCTGCCGTAAGGCCTCCGGTTAAGATGGCCCGATGCCGTAAAGCATCGATTTAAGGTGGCCCGCTGCGGGAACGCCTCGGGTTAC"
        + "TTTTCTCAGCGATCCACAAGTCGGGCTTTGCCGGATTAAAATATTTTGGCTAAACCAACCTCAGTG";

    private const string Band36 =
        "GTTTATTTTCCTCAGGGACTTAAGCAAGTCGCGGGTATAAATATTCAACCATTTATCCGGACCTTATCTAGCAACACTTCTAAACTTATTTCACGCCCCA"
        + "ACTCTAGCAACACTTCTAAACCTATTCACGCCCCAACTTCTAGCAACACTTCTAAACCTTTTCACAAGCTATCAGAGGTGGGGGGCTTTCTAGGAACGTT"
        + "TGTGGGAGGGATGCTTCCCGTATTGCCTATTTTACTTCGCTATGTACCAACTACATGTATTAGTCTAAACTGATAAATC";

    private const string Band264 =
        "ACTACAACGCTGGATGTATAGGGTGCTGGAGATACGCTCAGATGTGCGCATACTCTACACACGTCACCTGCTGCCTCACGCAGGTGTAATATGACAGAAT"
        + "TGAGTTGAACCAGATACCGGCGTCTTACTGAGCACTTGGCCCACCCCCGAAAGGCAGATCTCCGTCGATCCGGCAACGAGCCGTACAATTGAGGGAGAAG"
        + "TGCGGAATTTCCGCCATAGTGTTTGTAGATCCTCCGCCTACTTTAAACAATGGGAAGTGGGCGGCCCCTATAGAATTGTACCACCGTGCACCCTTGCCTC"
        + "GAAGTCTGTTGTCGAGAGAATAAATTAGATTAGATTTGGATGAAGAGGTCATAAATGTAAGCAATTGAGTTGACCAGTACCAGCGTCTTACTGGACTTTG"
        + "CCCACGCCCCGAAAGGCAGACCGCCGTCGGTCCGGCAAACAAAGATGTACAATTGGGGAGATGGTAAATGTCCGCAGCGAGTGTTTGTAGCTCCTTCCGC"
        + "CTATCTTCATACAATATGAAGTGGGCCGGCCCCTTAGAATTGGTACCACCAGTGCACCCCTCCTCGAAGTCTGTTGTCGGAAATGAAATTGGATTAGATT"
        + "TGGATGAAGAGGCCATAAATGTAAGAAGTGAGTTGACCAGATGCCGGCGTCTTACGGCACATGGCCCACGCCCCGACAAGGCAATCGCCGTAGCTCATAG"
        + "TCGCAGATCCGGGGGGTAAATCCTTACCCATGCTCAACGAGCCAGTCGAAGCACCCGAAA";

    private const string Spread28 =
        "CAGTCACGGGCTCTGGATCCAGCAGCAGTGCAGCATGTTGGTACCCTATCCCCATACGACACTGTTTGGCGCTGTTGGTTTATGCACGAGTCGTTACTAT"
        + "ATAAAGACCTCGAAGTGCCAGAATTCATCTTTGACCTCAGCGCGTTCGTACTCCGATCGGAACCGCCCGTTCACTGTACTCCGATCGGAACCGCCCCGAT"
        + "ATGTACTCCATTAATCGTCCCTTTGAATTCGGAGATACGCGTGACGGACGTATCGCGTCTCCATTCTTAGCCGACTCCACGACCTCCTTAATGGTTAATC"
        + "AACATAAGAATATTCCCAGGAG";

    private const string Spread22 =
        "CTTATTGGGGTACCTCCCGCCCAAATCGTCACTTTGTACCAGTCTACTCAGTGGTCCGGACGTACAGGGGTAATAATGAGGGATCCTAGAACTAACCTTC"
        + "TTGAACCCGCAGAACGTGCATTTGGGCCTATTCAGGAGTCCCACCTACGATATTCGGATTCCGCAGCTACGATATTTCAGGATCTTGCCACCAACGATAC"
        + "GTGGCCTGGCTGGGAAGACTCTCCACGCTGGTAGCGATTTGACGTACACGGCCGTCTCAGTTTATCAGTTTGTAACGTAAAACGGGTGTGTGCTCTAAAT"
        + "TACACGATCGTCGGGGGATCCTGGGAATTAGCATTCGACTGATGCCTCGCT";

    private const string Spread267 =
        "CTCCAGCTTGGAATTCACAACCCCCATCCGTGCGACACCTTCATCGAGCTGAAATGCGCGATCCTAGTGCCTATTCTTACCTAGGCTATTAGGTGACCCT"
        + "GCGTACTTCATCGGCTACAGCTAAGAGCCGGAGGTAGTCGGTATTGAGTTCCATGCGCGCTAATCGATGGGCCAGAAGTTCTAATTAGCGGCAGGACTGA"
        + "TATTTTTGAACGCAATACTAGAGGGTAGGTGTTTAAACCCTGTTTGATCTTAGTAAGGTGATCCCTAACTCTCGGGAGAGTGGTACCCTACCACATGCGT"
        + "CGTCAGCTTGCACTACTAAAAGGGAAGGATAGTGCCTATTCTTAACCCTAGGCTATTAAGTGACACTGCGTACTTCATCGGCTACAGCTAGAGCCGGAGT"
        + "AGTCGGTATTCGACTTCAATGCGGCAATCGATGGCAGAAGTTCTAATTAACGGCAAGGGGACTGATATTTTGAACGCGAATAATAAGGGTAGGTGTTAAA"
        + "CCCTGTTTCGATCTTCAAGTAAAGATGATCCCCTACTCTCGGTGAGAGTGGTACCCCTCACATGCAGTTCTGCAGCTGCACTACTAAAAGGGAAGACAGT"
        + "GCCTATTCTTATCCTAGGACTATTGAGTCGAACCTGCGGTATTCATCGGCTACAGCTATGACCAAAAGTGGGTGCGCCATCTCCGTCGGCGGCTAACCGG"
        + "TCCCATTACGTCCCCGGAGAGTC";

    private const string Active43 =
        "ACTACGTCCGAATTGTTCGCACTCCAGGACGATTTAGGAGTTCACTACAGTCCGAAGTTGTTCGCACTCCAGACGAATAGTAGTACACTACGGCCCGGAA"
        + "TTGTTCGCACTCCAGGACTGATTAGGAGTTCACTACGTGCGAATTGTTCGCCTCCAGGACGATTAGGAGTTACAC";

    private const string BestList278 =
        "TCCCGACGATAATTTTGCAGAGGACACTTTCAGCGCACGAGATCCTCGAGCTGGCTGCGTATCGTAAGAAGCCAGTTAACATTAAGAGTGAGCTCGGACT"
        + "GCCCTAGCGTGCGACAATTAAGAGTTGACATAACGAGTTATGAAGTGGTCCCTCATTCTCCCACGAGGGGCGGTCGCCTTGACGCTCTCGATGTAAAGAG"
        + "TGTCTTGGCTCTGGGTGTACCCGACCTTCTTTAATCGGACTCTGAGGCATGCTGTCCACCCATGTCCAAGAGGAGCCCTAAAGCATGCAGCAATAACGGC"
        + "GTAGATGTGCCAGAGATCCCCCAGCTGGGTCCGTATCTTAAGAAGCCAGTTACGATAAGAGTGAGCGCGACTGGCCTAGCGTGCAACAATTAAGAGTTGA"
        + "ATTCGCGAACTATGGAGTGTCCCTCACTCTACCAGGAGGGGCGGACGCCTTGTCGGCTCTCGATGTAATGAGTGTCACTGGCTCTGCGTGGTACCCGCCT"
        + "TCTGTATTCGGACTCTGAGGCATGTTCAACACTCTCCGACTGAGGAGCCCGATGCATGCAGACAATAACCGCGTAGATGTGGTACGAGATCCCGGGTTGC"
        + "GAGTTCGTATCTTAAGAAGCCAGTACGTTAAGAGTCGAGCGCGGACTCCCCTAGCGTGAAACAATTAAGAGTTGAATCGGCGAATTCTGGAGTGGTCCCT"
        + "CATTCTCTCAAGAGGGGCGGTCGCCTTGACCTCTCGATGTGATGAGTGGTCTTGGCTCTGGGTTGTACCCGTCCATTCTTTAATCGGACTCTGAGGCATG"
        + "CTGTCACCGATGTCCAACTGAGGAGCCCGCTGCTGCAGCAATAACTGGCGTAAGATGTTGCCACGAGATCCCGAGCTTGAGTTGGGAGCGCAAGAACTAG"
        + "TAACGTTTAGGGTGAGCGCGGACTGCCCTAGCGTGCAACAATTAGGAGCTCGAAATCGCGAATATGGAGTGGTCCCTTATTCTCCCAAGAAGGGGGATCG"
        + "CCTTGACGCTCTTGATGTACTGAGTGTGCTGTTGCTCTGGGTTGTACCCGCGCTTCTTTTATCGGACTGTGAGGCAATGCTTGTGACCCCAGTCACACTG"
        + "AAGAGCCCGATGCATTCAGCAAAACGGCGTAGATGTGCCCGGAGATCCCGAGCTGTGAGTTAGTATCGTAAGAAAGCCAGTTACTCTAAGAGTGAGAGCG"
        + "GGCTATGATGAGTCTCGAGTCGGTATAACCCCCATGAAACCTGAGGTCATTTCAAGCGGACCACAACACCTACTTACAAGCAGTAAGTATTAAGCAGTTG"
        + "AGATACATAATGTGTCTCAAGTCGGTATTAGCCCCACGAAACCTGAGGTCATCAACATGCGGACAACAACACCT";

    private static TandemRepeatsFinderParameters Set(int match, int mismatch, int delta, int pm, int pi, int minScore, int maxPeriod) =>
        new()
        {
            MatchWeight = match, MismatchPenalty = mismatch, IndelPenalty = delta,
            MatchProbability = pm, IndelProbability = pi, MinScore = minScore, MaxPeriod = maxPeriod,
        };

    /// <summary>Asserts one row against TRF's printed .dat values.</summary>
    private static void AssertTrfRow(
        ApproximateTandemRepeatResult r, int first1, int last1, int period, double copies, int size,
        int matches, int indels, int score, int a, int c, int g, int t, double entropy, string consensus)
    {
        Assert.Multiple(() =>
        {
            Assert.That(r.Start + 1, Is.EqualTo(first1), "TRF first index (1-based)");
            Assert.That(r.Start + r.SpanLength, Is.EqualTo(last1), "TRF last index");
            Assert.That(r.Period, Is.EqualTo(period), "period");
            Assert.That(r.CopyNumber, Is.EqualTo(copies).Within(0.05 + 1e-9), "copy number (TRF prints %.1f)");
            Assert.That(r.ConsensusSize, Is.EqualTo(size), "consensus size");
            Assert.That((int)(r.PercentMatches + 1e-9), Is.EqualTo(matches), "% matches (TRF truncates)");
            Assert.That((int)(r.PercentIndels + 1e-9), Is.EqualTo(indels), "% indels (TRF truncates)");
            Assert.That(r.AlignmentScore, Is.EqualTo(score), "score");
            Assert.That((int)(r.PercentA + 1e-9), Is.EqualTo(a), "A");
            Assert.That((int)(r.PercentC + 1e-9), Is.EqualTo(c), "C");
            Assert.That((int)(r.PercentG + 1e-9), Is.EqualTo(g), "G");
            Assert.That((int)(r.PercentT + 1e-9), Is.EqualTo(t), "T");
            Assert.That(Math.Round(r.EntropyTrf, 2), Is.EqualTo(entropy), "TRF entropy (%.2f)");
            Assert.That(r.Consensus, Is.EqualTo(consensus), "consensus");
        });
    }

    #region D1-D3 — apparent-size criterion (exact derivation of the distribution TRF simulates)

    // TRF README, "Apparent Size Distribution": "if PM = .75, k = 5 and d = 100, then the criterion is 56".
    [Test]
    public void ApparentSize_ReadmeExample_Pm75_K5_D100_Is56()
    {
        Assert.That(RepeatFinder.TrfApparentSize(100, 75), Is.EqualTo(56));
    }

    // Offsets (max(d,20) - y - 1) equal to TRF 4.10.0's simulated waitdata tables at these d (the derived table
    // agrees with TRF's entry exactly at 825/2000 (PM 80) and 713/2000 (PM 75) distances and within its simulation
    // noise elsewhere: |difference| <= 3 / 4, mean +0.47; Evidence §WP7).
    [TestCase(1, 80, 18)]
    [TestCase(20, 80, 18)]
    [TestCase(21, 80, 19)]
    [TestCase(30, 80, 28)]
    [TestCase(44, 80, 29)]
    [TestCase(60, 80, 31)]
    [TestCase(250, 80, 65)]
    [TestCase(1000, 80, 67)]
    [TestCase(2000, 80, 68)]
    [TestCase(1, 75, 15)]
    [TestCase(25, 75, 15)]
    [TestCase(30, 75, 25)]
    [TestCase(50, 75, 39)]
    [TestCase(150, 75, 44)]
    public void ApparentSizeOffset_EqualsTrfWaitTableEntry(int d, int pm, int trfEntry)
    {
        Assert.That(RepeatFinder.TrfApparentSizeOffset(d, pm), Is.EqualTo(trfEntry));
    }

    [Test]
    public void ApparentSizeOffset_LiesInsideTheDistanceWindow_ForEveryDistance()
    {
        Assert.Multiple(() =>
        {
            foreach (int pm in new[] { 80, 75 })
            {
                for (int d = 1; d <= RepeatFinder.MaxApproximatePeriod; d++)
                {
                    int offset = RepeatFinder.TrfApparentSizeOffset(d, pm);
                    Assert.That(offset, Is.InRange(1, Math.Max(d, 20) - 1), $"d={d}, PM={pm}");
                }
            }
        });
    }


    // TRF reports nothing; without the apparent-size test a period-28 candidate (148..209, score 92) whose tuple matches all
    // lie at the right end of the distance window is analysed and reported (previous implementation).
    [Test]
    public void ApparentSize_RejectsClusteredMatches_Period28()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(Spread28, Set(2, 7, 7, 80, 10, 50, 500)).ToList();
        Assert.That(rows, Is.Empty, "TRF reports no repeat");
    }


    // As above for a period-22 candidate (129..182, score 67).
    [Test]
    public void ApparentSize_RejectsClusteredMatches_Period22()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(Spread22, Set(2, 7, 7, 80, 10, 50, 500)).ToList();
        Assert.That(rows, Is.Empty, "TRF reports no repeat");
    }


    // The apparent-size test changes which candidate of this period-267 array is examined first; without it the reported
    // row has a 265-bp consensus and score 790 instead of TRF's 264 / 788.
    [Test]
    public void ApparentSize_DelaysCandidate_Period267()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(Spread267, Set(2, 7, 7, 80, 10, 50, 500)).ToList();
        Assert.That(rows, Has.Count.EqualTo(1));
        var r = rows[0];
        AssertTrfRow(r, 65, 658, 267, 2.2, 264, 85, 11, 788, 25, 22, 24, 27, 2.00, "TAGTGCCTATTCTTACCTAGGCTATTAAGTGACACTGCGTACTTCATCGGCTACAGCTAGAGCCGGAGGTAGTCGGTATTGACTTCAATGCGCGCTAATCGATGGGCCAGAAGTTCTAATTAACGGCAGGACTGATATTTTTGAACGCAATAATAGAGGGTAGGTGTTTAAACCCTGTTTGATCTTAGTAAGATGATCCCTAACTCTCGGTGAGAGTGGTACCCCTCACATGCAGTTCTGCAGCTGCACTACTAAAAGGGAAGA");
    }


    #endregion

    #region D4-D6 — narrow-band WDP for patterns > 20


    // Full (optimal) WDP: 102..162, consensus 31, score 95; TRF's band (radius max(6, floor(2.3*sqrt(.1*31))) = 6
    // backward, min(12, 10) = 10 forward) recentres on the matching diagonal and reaches 102..182, score 99.
    [Test]
    public void NarrowBand_Period31_ReproducesTrfExtentAndScore()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(Band31, Set(2, 7, 7, 80, 10, 50, 500)).ToList();
        Assert.That(rows, Has.Count.EqualTo(1));
        var r = rows[0];
        AssertTrfRow(r, 102, 182, 31, 2.6, 30, 84, 3, 99, 19, 27, 32, 20, 1.97, "GGTTAAGATGGCCCGCTGCCGTAAAGCATC");
    }


    // Band alignment changes the majority-rule consensus (37 bp vs 36 bp with the full WDP) and the score (168 vs 164).
    [Test]
    public void NarrowBand_Period36_ReproducesTrfConsensus()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(Band36, Set(2, 7, 7, 80, 10, 50, 500)).ToList();
        Assert.That(rows, Has.Count.EqualTo(1));
        var r = rows[0];
        AssertTrfRow(r, 67, 165, 36, 2.8, 37, 92, 6, 168, 30, 36, 5, 28, 1.79, "TCTAGCAACACTTCTAAACCTATTTCACGCCCCAACT");
    }


    // Large pattern (radius floor(2.3*sqrt(.1*264)) = 11 backward, 22 forward): TRF ends the repeat at 691 (score 761); the full WDP stops at 683 (752).
    [Test]
    public void NarrowBand_Period264_ReproducesTrfExtent()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(Band264, Set(2, 7, 7, 80, 10, 50, 500)).ToList();
        Assert.That(rows, Has.Count.EqualTo(1));
        var r = rows[0];
        AssertTrfRow(r, 96, 691, 264, 2.3, 263, 83, 9, 761, 26, 24, 26, 22, 2.00, "AGAATTGAGTTGACCAGATACCGGCGTCTTACTGGCACTTGGCCCACGCCCCGAAAGGCAGATCGCCGTCGATCCGGCAACAAGACGTACAATTGAGGGAGAAGTGCGAAATGTCCGCCAGAGTGTTTGTAGATCCTCCGCCTACTTCAAACAATAGGAAGTGGGCGGCCCCTATAGAATTGTACCACCGTGCACCCCTGCCTCGAAGTCTGTTGTCGAGAGAATAAATTAGATTAGATTTGGATGAAGAGGCCATAAATGTA");
    }


    [Test]
    public void NarrowBand_LegacyOverload_UsesTheSameAnalysis()
    {
        var legacy = RepeatFinder.FindApproximateTandemRepeats(Band31, 1, 500, 50).ToList();
        Assert.That(legacy, Has.Count.EqualTo(1));
        var r = legacy[0];
        AssertTrfRow(r, 102, 182, 31, 2.6, 30, 84, 3, 99, 19, 27, 32, 20, 1.97, "GGTTAAGATGGCCCGCTGCCGTAAAGCATC");
    }

    #endregion

    #region D7-D8 — active distances and the best-period list


    // A distance joins the random-walk range sums of its neighbours only after it has itself reached a criteria test
    // (TRF links a distance list when it is tested). Summing every distance with matches lets d = 23 pass and adds a
    // spurious 20..175 period-23 row (and, previously, a period-88 row).
    [Test]
    public void RangeSums_UseOnlyActiveDistances()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(Active43, Set(2, 5, 5, 75, 10, 30, 100)).ToList();
        Assert.That(rows, Has.Count.EqualTo(1));
        var r = rows[0];
        AssertTrfRow(r, 1, 175, 43, 4.0, 42, 83, 9, 256, 25, 25, 24, 25, 2.00, "ACTACGTCCGAATTGTTCGCACTCCAGGACGATTAGGAGTTC");
    }


    // Distances > 250 inside a region analysed earlier must be among that region's five best periods (TRF best-period
    // list). Without it a period-278 candidate is examined earlier and the reported repeat starts at 36 (score 1185).
    [Test]
    public void BestPeriodList_BlocksMultipleInsideAnalysedRegion()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(BestList278, Set(2, 7, 7, 80, 10, 50, 2000)).ToList();
        Assert.That(rows, Has.Count.EqualTo(1));
        var r = rows[0];
        AssertTrfRow(r, 58, 1201, 278, 4.1, 275, 79, 10, 1140, 22, 24, 28, 23, 1.99, "CGTATCTTAAGAAGCCAGTACGTTAAGAGTGAGCGCGGACTGCCCTAGCGTGCAACAATTAAGAGTTGAAATCGCGAACTATGGAGTGGTCCCTCATTCTCCCAAGAGGGGCGATCGCCTTGACGCTCTCGATGTAATGAGTGTGCTTGGCTCTGGGTTGTACCCGACCTTCTTTAATCGGACTCTGAGGCATGCTGTACACCCATGTCCAACTGAGGAGCCCGATGCATGCAGCAATAACGGCGTAGATGTGCCCGAGATCCCGAGCTTGAGTT");
    }


    #endregion
}
