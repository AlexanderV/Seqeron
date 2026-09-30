// REP-APPROX-001 — Approximate (Imperfect/Interrupted) Tandem-Repeat Detection (Tandem Repeats Finder model)
// Evidence: docs/Evidence/REP-APPROX-001-Evidence.md
// TestSpec: tests/TestSpecs/REP-APPROX-001.md
// Sources: Benson G (1999). Tandem repeats finder: a program to analyze DNA sequences.
//          Nucleic Acids Research 27(2):573-580. https://doi.org/10.1093/nar/27.2.573
//          TRF 4.10.0 README (github.com/Benson-Genomics-Lab/TRF): parameters, "Table Explanation",
//          "Alignment Explanation", "How does Tandem Repeats Finder work?", test_seqs expected tables.
//
// Every expected row below was produced by the compiled TRF 4.10.0 binary
// (`trf seq.fa 2 7 7 80 10 <minscore> <maxperiod> -h -d`, .dat columns: start end period copies
// consensus-size %match %indel score A C G T entropy consensus). TRF prints 1-based indices and truncates
// percentages to integers; the library reports 0-based Start and exact percentages, so exact values are
// asserted and the TRF integer is noted alongside. Counts behind the exact percentages (matches /
// mismatches / indels between adjacent copies) come from an instrumented TRF build (same algorithm,
// extra stderr print in get_statistics).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RepeatFinder_ApproximateTandemRepeats_Tests
{
    private const int LowMinScore = 10;
    private const double Tol = 1e-9;

    private static void AssertRow(
        ApproximateTandemRepeatResult r, int trfStart1, int trfEnd1, int period, double copies, int consensusSize,
        double percentMatches, double percentIndels, int score, string consensus)
    {
        Assert.Multiple(() =>
        {
            Assert.That(r.Start, Is.EqualTo(trfStart1 - 1), "0-based start = TRF 1-based start - 1");
            Assert.That(r.Start + r.SpanLength, Is.EqualTo(trfEnd1), "exclusive end = TRF 1-based inclusive end");
            Assert.That(r.Period, Is.EqualTo(period), "period");
            Assert.That(r.CopyNumber, Is.EqualTo(copies).Within(Tol), "copy number");
            Assert.That(r.ConsensusSize, Is.EqualTo(consensusSize), "consensus size");
            Assert.That(r.PercentMatches, Is.EqualTo(percentMatches).Within(Tol), "% matches between adjacent copies");
            Assert.That(r.PercentIndels, Is.EqualTo(percentIndels).Within(Tol), "% indels between adjacent copies");
            Assert.That(r.AlignmentScore, Is.EqualTo(score), "WDP alignment score");
            Assert.That(r.Consensus, Is.EqualTo(consensus), "consensus (TRF phase)");
        });
    }

    #region TRF README test sequences (published expected tables)

    // TRF README "Testing the Installation": test_seqs s1..s3 (`trf test_seqs.fasta 2 5 7 80 10 50 2000`) —
    // s1 1--35 period 7 copies 5.0 size 7 100 0 score 70 A14 C28 G28 T28 entropy 1.95;
    // s2 1--84 12 7.0 12 100 0 168 16 41 25 16 1.89; s3 1--1225 35 35.0 35 100 0 2450 28 22 20 28 1.98.
    // (Perfect arrays: identical with the recommended 2 7 7 parameters, confirmed with the binary.)
    [Test]
    public void FindApproximateTandemRepeats_TrfReadmeTestSequences_ReproducePublishedTables()
    {
        const string s1 = "TCATCGGTCATCGGTCATCGGTCATCGGTCATCGG";
        string s2 = string.Concat(Enumerable.Repeat("ACCCCTCAGGGT", 7));
        string s3 = string.Concat(Enumerable.Repeat("TGACTATATCCGCAAATGAAGGCTGTTCTCTGACA", 35));

        var r1 = RepeatFinder.FindApproximateTandemRepeats(s1, 1, 2000).Single();
        var r2 = RepeatFinder.FindApproximateTandemRepeats(s2, 1, 2000).Single();
        var r3 = RepeatFinder.FindApproximateTandemRepeats(s3, 1, 2000).Single();

        AssertRow(r1, 1, 35, 7, 5.0, 7, 100, 0, 70, "TCATCGG");
        AssertRow(r2, 1, 84, 12, 7.0, 12, 100, 0, 168, "ACCCCTCAGGGT");
        AssertRow(r3, 1, 1225, 35, 35.0, 35, 100, 0, 2450, "TGACTATATCCGCAAATGAAGGCTGTTCTCTGACA");
        Assert.Multiple(() =>
        {
            Assert.That(r1.PercentA, Is.EqualTo(100.0 * 5 / 35).Within(Tol), "TRF A 14");
            Assert.That(r1.PercentC, Is.EqualTo(100.0 * 10 / 35).Within(Tol), "TRF C 28");
            Assert.That(r1.PercentG, Is.EqualTo(100.0 * 10 / 35).Within(Tol), "TRF G 28");
            Assert.That(r1.PercentT, Is.EqualTo(100.0 * 10 / 35).Within(Tol), "TRF T 28");
            Assert.That(Math.Round(r1.Entropy, 2), Is.EqualTo(1.95), "TRF entropy 1.95");
            Assert.That(Math.Round(r2.Entropy, 2), Is.EqualTo(1.89), "TRF entropy 1.89");
            Assert.That(Math.Round(r3.Entropy, 2), Is.EqualTo(1.98), "TRF entropy 1.98");
            Assert.That((int)r2.PercentC, Is.EqualTo(41), "TRF C 41");
            Assert.That((int)r3.PercentG, Is.EqualTo(20), "TRF G 20");
        });
    }

    #endregion

    #region Perfect and substituted tracts

    // A1 — perfect (CA)x5, minscore 10. TRF: 1 10 2 5.0 2 100 0 20 50 50 0 0 1.00 CA.
    [Test]
    public void FindApproximateTandemRepeats_PerfectDinucleotide_MatchesTrf()
    {
        var r = RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 1, 6, LowMinScore).Single();

        AssertRow(r, 1, 10, 2, 5.0, 2, 100, 0, 20, "CA");
        Assert.That(r.Entropy, Is.EqualTo(1.0).Within(Tol), "50/50 composition = 1 bit (TRF 1.00)");
    }

    // A2 — (CAG)x6 with copy 4 = TAG. TRF: 1 18 3 6.0 3 86 0 27 33 27 33 5 1.80 CAG.
    // Adjacent-copy comparison (instrumented TRF): 13 matches, 2 mismatches, 0 indels -> 86.67 % (TRF prints 86).
    // The OLD implementation reported 94.44 % = 17/18 (sequence vs consensus), contrary to the TRF definition
    // "Percent of matches between adjacent copies overall".
    [Test]
    public void FindApproximateTandemRepeats_OneSubstitution_PercentMatchesBetweenAdjacentCopies()
    {
        var r = RepeatFinder.FindApproximateTandemRepeats("CAGCAGCAGTAGCAGCAG", 3, 3, LowMinScore).Single();

        AssertRow(r, 1, 18, 3, 6.0, 3, 1300.0 / 15, 0, 27, "CAG");
        Assert.That(r.PercentT, Is.EqualTo(100.0 / 18).Within(Tol), "TRF T 5");
    }

    // A2b — the same interrupted tract under the PERFECT detector fragments into one short CAG x3 run.
    [Test]
    public void FindMicrosatellites_OnInterruptedTract_FragmentsWherePerfectDetectorBreaks()
    {
        var perfect = RepeatFinder.FindMicrosatellites("CAGCAGCAGTAGCAGCAG", 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(perfect, Has.Count.EqualTo(1), "perfect detector reports only the leading run");
            Assert.That(perfect[0].RepeatUnit, Is.EqualTo("CAG"));
            Assert.That(perfect[0].RepeatCount, Is.EqualTo(3), "CAG x3 before the TAG interruption");
            Assert.That(perfect[0].Position, Is.EqualTo(0));
        });
    }

    // A3 — (CA)x6 with index 6 = T. TRF reports nothing even at minscore 10: the only run of >= 4 matches at
    // distance 2 has 4 heads, below the sum-of-heads criterion 5 for d = 2 (k = 4, floor k+1).
    [Test]
    public void FindApproximateTandemRepeats_ShortInterruptedTract_BelowSumOfHeadsCriterion_NotReported()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindApproximateTandemRepeats("CACACATACACA", 2, 2, LowMinScore), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats("CACACATACACA", 1, 500, LowMinScore), Is.Empty);
        });
    }

    #endregion

    #region Indels and overlapping periods

    // A4 — (CAG)x10 with one base deleted (29 bp). TRF: 1 29 3 10.0 3 92 7 51 34 31 34 0 1.58 CAG.
    // Adjacent copies: 25 matches, 0 mismatches, 2 indels (27 trials) -> 92.59 % / 7.41 %. Copy number 10.0
    // (30 aligned consensus columns / 3). OLD implementation: 9.667 copies, 96.67 % / 3.33 % (vs consensus).
    [Test]
    public void FindApproximateTandemRepeats_OneDeletion_MatchesTrf()
    {
        var r = RepeatFinder.FindApproximateTandemRepeats("CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", 3, 3, LowMinScore).Single();

        AssertRow(r, 1, 29, 3, 10.0, 3, 2500.0 / 27, 200.0 / 27, 51, "CAG");
    }

    // A4b — same tract, maxPeriod 500, minscore 10: TRF reports three overlapping repeats of different
    // periods (redundancy only removes multiples of a period / same-period duplicates):
    //   1 29 3 10.0 3 92 7 51 CAG;  5 26 11 2.0 11 100 0 44 AGCAGCAGCAG;  2 29 14 2.0 14 100 0 56 AGCAGCAGCAGCAG.
    // At minscore 50 only the 51 and 56 rows remain. Output is ordered by start.
    [Test]
    public void FindApproximateTandemRepeats_DeletionTract_ReportsOverlappingPeriodsLikeTrf()
    {
        const string seq = "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG";
        var low = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 500, LowMinScore).ToList();
        var def = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 500).ToList();

        Assert.That(low, Has.Count.EqualTo(3));
        AssertRow(low[0], 1, 29, 3, 10.0, 3, 2500.0 / 27, 200.0 / 27, 51, "CAG");
        AssertRow(low[1], 2, 29, 14, 2.0, 14, 100, 0, 56, "AGCAGCAGCAGCAG");
        AssertRow(low[2], 5, 26, 11, 2.0, 11, 100, 0, 44, "AGCAGCAGCAG");
        Assert.That(def.Select(r => (r.Period, r.AlignmentScore)), Is.EqualTo(new[] { (3, 51), (14, 56) }));
    }

    // A4c — flanked, and lowercase: TRF uppercases input; indices shift by the 2-bp flank
    // (TRF: 3 31 3 ...; 7 28 11 ...; 4 31 14 ...).
    [Test]
    public void FindApproximateTandemRepeats_FlankedLowercase_SameRowsShiftedIndices()
    {
        const string upper = "TTCAGCAGCAGCAGCAGAGCAGCAGCAGCAGTTGACCA";
        var a = RepeatFinder.FindApproximateTandemRepeats(upper, 1, 500, 20).ToList();
        var b = RepeatFinder.FindApproximateTandemRepeats(upper.ToLowerInvariant(), 1, 500, 20).ToList();

        Assert.That(b, Is.EqualTo(a), "case-insensitive");
        Assert.That(a.Select(r => (r.Start + 1, r.Start + r.SpanLength, r.Period, r.AlignmentScore)),
            Is.EqualTo(new[] { (3, 31, 3, 51), (4, 31, 14, 56), (7, 28, 11, 44) }));
    }

    #endregion

    #region N, period 1, redundancy, random snapshots

    // A14 — TRF's scoring matrix matches only identical A/C/G/T ("changed to use Similarity Matrix to avoid N
    // matching itself", TRF source 2006); an all-N run is not a repeat. OLD implementation: N matched N and
    // "NNNNNNNN" was a perfect period-1 repeat.
    [Test]
    public void FindApproximateTandemRepeats_AllN_NotReported()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindApproximateTandemRepeats("NNNNNNNN", 1, 6, LowMinScore), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(new string('N', 200), 1, 6, LowMinScore), Is.Empty);
        });
    }

    // N inside a repeat is a mismatch. TRF: 5 34 3 10.0 3 92 0 51 33 30 33 0 1.58 CAG
    // (the N column is a mismatch between adjacent copies; composition denominators include the N).
    [Test]
    public void FindApproximateTandemRepeats_NInsideRepeat_IsMismatch()
    {
        var r = RepeatFinder.FindApproximateTandemRepeats("GATCCAGCAGCAGCAGNAGCAGCAGCAGCAGCAGGTTACG", 1, 500).Single();

        AssertRow(r, 5, 34, 3, 10.0, 3, 2500.0 / 27, 0, 51, "CAG");
        Assert.That(r.PercentC, Is.EqualTo(30.0).Within(Tol), "9 C of 30 region bases (TRF C 30)");
    }

    // Period 1 needs >= 80 % of one base (TRF). A29 G A26: TRF 1 57 1 57.0 1 96 0 105 98 0 1 0 0.13 A.
    [Test]
    public void FindApproximateTandemRepeats_Homopolymer_MatchesTrf()
    {
        string seq = new string('A', 30) + "G" + new string('A', 26);
        var r = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 500).Single();

        AssertRow(r, 1, 57, 1, 57.0, 1, 2700.0 / 28, 0, 105, "A");
        Assert.That(Math.Round(r.Entropy, 2), Is.EqualTo(0.13), "TRF entropy 0.13");
    }

    // (CA)x30: TRF reports only period 2 (1 60 2 30.0 2 100 0 120); periods 4, 6 of the same array are removed
    // as redundant multiples (score <= 1.1 x).
    [Test]
    public void FindApproximateTandemRepeats_LongMicrosatellite_MultiplesRemovedAsRedundant()
    {
        string seq = string.Concat(Enumerable.Repeat("CA", 30));
        var r = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 500).Single();

        AssertRow(r, 1, 60, 2, 30.0, 2, 100, 0, 120, "CA");
    }

    // Random sequences with embedded approximate repeats (generated, then run through TRF 4.10.0 with
    // `2 7 7 80 10 50 500 -h -d`); every TRF row is reproduced (indices, period, copies, size, %match,
    // %indel, score, consensus). Exact % from adjacent-copy counts; TRF integers in comments.
    [Test]
    public void FindApproximateTandemRepeats_RandomEmbeddedRepeats_ReproduceTrfRows()
    {
        const string seqA = "TCGATCGTTGATTAACGCTTCACTACTTTCGCCTAACTTATCGCACTAACTTCGCACTAATTCCGCATGAACTTCTCGTACTAACTTTGCGCACTAACTTTCGCCTAACCTTTCGCACTAAATTTCTGAATGAATCAGCTATTATTAGGGAGTATTT";
        const string seqB = "CCATGCCTCGTGCCGCATCAACCGGTGCCCGGTACCCGGGTCCCCGGTAACCCGGTACCCTGGGTACTCGGTGCCCGGTACCCGGTACCCCGGTACCCGGTACCCGGTACCGCGGTACCCGTACAGCGTCCCACTTACACAATATACATGCCCGGGTTGGCCT";
        const string seqC = "AAGCACAGAAAAGTCTGCGACCTTGTATCGTCCGTGCGCGCTCCCTTTCTGCGTATACGTAGCATGTCAATACTCGACTTACTCGGCTTACTCGGCTCACTCGGCCTACTAGGCTAACTCGTTTTACTCGGCTTACTCGGCTGATCCCGGTATCTTTAGCAAGGTTGCATACAGAAATTCGTTCCACTATTGCGCTCGTTCCACTATTGCGTCGTTCACTATTGCGGCGGCACC";

        var a = RepeatFinder.FindApproximateTandemRepeats(seqA, 1, 500).ToList();
        var b = RepeatFinder.FindApproximateTandemRepeats(seqB, 1, 500).ToList();
        var c = RepeatFinder.FindApproximateTandemRepeats(seqC, 1, 500).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(a, Has.Count.EqualTo(2));
            Assert.That(b, Has.Count.EqualTo(3));
            Assert.That(c, Has.Count.EqualTo(2));
        });
        // TRF: 21 126 13 8.9 12 75 20 98 CACTAACTTTCG / 21 126 24 4.5 23 73 20 94 CACTAACTTTCGCCTAACTTTCG
        AssertRow(a[0], 21, 126, 13, 107.0 / 12, 12, 75.23809523809524, 20.0, 98, "CACTAACTTTCG");
        AssertRow(a[1], 21, 126, 24, 103.0 / 23, 23, 73.91304347826087, 20.652173913043477, 94, "CACTAACTTTCGCCTAACTTTCG");
        // TRF: 22 121 7 13.4 7 80 12 110 / 22 121 15 6.7 14 79 12 92 / 22 121 22 4.4 23 82 9 118
        AssertRow(b[0], 22, 121, 7, 94.0 / 7, 7, 80.8080808080808, 12.121212121212121, 110, "CCGGTAC");
        AssertRow(b[1], 22, 121, 15, 94.0 / 14, 14, 79.12087912087912, 12.087912087912088, 92, "CCGGTCCCCGGTAC");
        AssertRow(b[2], 22, 121, 22, 102.0 / 23, 23, 82.92682926829268, 9.75609756097561, 118, "CCGGTACCCGGTACCCGGGTACC");
        // TRF: 71 142 9 8.0 9 79 0 81 TACTCGGCT / 179 226 17 2.9 17 93 6 82 TCGTTCCACTATTGCGC
        AssertRow(c[0], 71, 142, 9, 8.0, 9, 79.36507936507937, 0, 81, "TACTCGGCT");
        AssertRow(c[1], 179, 226, 17, 50.0 / 17, 17, 93.93939393939394, 6.0606060606060606, 82, "TCGTTCCACTATTGCGC");
    }

    #endregion

    #region Minscore threshold

    // A5 — Benson (1999): "Only those repeats scoring at least 50 ... are reported." CA x5 scores 20.
    [Test]
    public void FindApproximateTandemRepeats_DefaultMinScore_SuppressesLowScoringTract()
    {
        Assert.That(RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 1, 6), Is.Empty);
    }

    // A6 — the deletion tract scores 51 >= 50.
    [Test]
    public void FindApproximateTandemRepeats_DefaultMinScore_ReportsTractAtOrAboveThreshold()
    {
        var results = RepeatFinder.FindApproximateTandemRepeats("CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", 3, 3).ToList();

        Assert.That(results.Single().AlignmentScore, Is.EqualTo(51));
        Assert.That(RepeatFinder.DefaultApproximateMinScore, Is.EqualTo(50));
        Assert.That(RepeatFinder.FindApproximateTandemRepeats("CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", 3, 3, 52), Is.Empty,
            "score 51 < 52");
    }

    #endregion

    #region Sum-of-heads criterion

    // Benson (1999): R(d,k,PM) = heads in runs >= k of a Bernoulli(PM) sequence of length d, normal
    // approximation, 95 % cut-off; tuple sizes k = 4 / 5 / 7 for d <= 29 / 30..159 / >= 160; minimum k+1.
    // Values = TRF 4.10.0 `sumdata80` table (tr30dat.c), which the exact-moment computation reproduces for all
    // d = 1..2000.
    [TestCase(1, 5)]
    [TestCase(22, 5)]
    [TestCase(23, 6)]
    [TestCase(29, 9)]
    [TestCase(30, 6)]
    [TestCase(50, 15)]
    [TestCase(100, 39)]
    [TestCase(159, 69)]
    [TestCase(160, 43)]
    [TestCase(500, 177)]
    [TestCase(2000, 818)]
    public void TrfSumOfHeadsCriterion_MatchesTrfTable(int d, int expected)
    {
        Assert.That(RepeatFinder.TrfSumOfHeadsCriterion(d), Is.EqualTo(expected));
    }

    #endregion

    #region Edge cases and validation

    [Test]
    public void FindApproximateTandemRepeats_EmptyOrNoRepeat_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindApproximateTandemRepeats("", 1, 6, LowMinScore), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats((string)null!, 1, 6, LowMinScore), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats("ACGTGCAT", 1, 6, LowMinScore), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats("A", 1, 6, LowMinScore), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats("ACG", 6, 6, LowMinScore), Is.Empty);
        });
    }

    // Validation is eager on both overloads (also for empty input): minPeriod >= 1, maxPeriod >= minPeriod,
    // maxPeriod <= 2000 (TRF 4.10.0 README: MaxPeriod over 2000 is an error), minScore >= 1 (TRF: "all weights,
    // penalties, and scores are positive").
    [Test]
    public void FindApproximateTandemRepeats_InvalidParameters_ThrowEagerly()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 0, 6, LowMinScore));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 4, 2, LowMinScore));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 1, 2001, LowMinScore));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 1, 6, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats("", 0, 6));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats(new DnaSequence("ACGT"), 1, 2001));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindApproximateTandemRepeats((DnaSequence)null!, 1, 6));
        });
    }

    [Test]
    public void FindApproximateTandemRepeats_RunTwiceAndDnaOverload_Agree()
    {
        const string seq = "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG";
        var first = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 20, LowMinScore).ToList();
        var second = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 20, LowMinScore).ToList();
        var viaDna = RepeatFinder.FindApproximateTandemRepeats(new DnaSequence(seq), 1, 20, LowMinScore).ToList();

        Assert.That(second, Is.EqualTo(first));
        Assert.That(viaDna, Is.EqualTo(first));
    }

    // minPeriod filters the reported period after redundancy elimination: it never resurrects a redundant
    // multiple. (CA)x30 with minPeriod 3 -> nothing (TRF's only row has period 2).
    [Test]
    public void FindApproximateTandemRepeats_MinPeriod_DoesNotResurrectRedundantMultiples()
    {
        string seq = string.Concat(Enumerable.Repeat("CA", 30));
        Assert.That(RepeatFinder.FindApproximateTandemRepeats(seq, 3, 500), Is.Empty);
    }

    #endregion

    #region ComputeBernoulliStatistics — TRF Bernoulli model

    // Benson (1999): alignment of two adjacent copies = Bernoulli trials, heads = match, PM = P(heads),
    // PI = indel probability; TRF statistics are "between adjacent copies ... not between the sequence and
    // the consensus pattern". Expected counts = TRF 4.10.0 (instrumented get_statistics) on the same tract.

    // B1 — perfect CA x5: TRF m=8 mm=0 ind=0.
    [Test]
    public void ComputeBernoulliStatistics_PerfectTract_MatchProbabilityIsOne()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CACACACACA", period: 2);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Period, Is.EqualTo(2));
            Assert.That(stats.AdjacentCopyPairs, Is.EqualTo(4), "5 copies -> 4 adjacent pairs");
            Assert.That(stats.BernoulliTrials, Is.EqualTo(8));
            Assert.That(stats.Matches, Is.EqualTo(8));
            Assert.That(stats.Mismatches, Is.EqualTo(0));
            Assert.That(stats.Indels, Is.EqualTo(0));
            Assert.That(stats.MatchProbability, Is.EqualTo(1.0).Within(Tol));
            Assert.That(stats.IndelProbability, Is.EqualTo(0.0).Within(Tol));
            Assert.That(stats.PercentMatches, Is.EqualTo(100.0).Within(Tol));
            Assert.That(stats.ExpectedMatches, Is.EqualTo(8.0).Within(Tol));
            Assert.That(stats.MeetsExpectedMatchProbability, Is.True);
        });
    }

    // B2 — CAG x6 with copy 4 = TAG: TRF m=13 mm=2 ind=0 (86 %).
    [Test]
    public void ComputeBernoulliStatistics_OneSubstitution_EstimatesAdjacentCopyPm()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGTAGCAGCAG", period: 3);

        Assert.Multiple(() =>
        {
            Assert.That(stats.AdjacentCopyPairs, Is.EqualTo(5));
            Assert.That(stats.BernoulliTrials, Is.EqualTo(15));
            Assert.That(stats.Matches, Is.EqualTo(13));
            Assert.That(stats.Mismatches, Is.EqualTo(2));
            Assert.That(stats.Indels, Is.EqualTo(0));
            Assert.That(stats.MatchProbability, Is.EqualTo(13.0 / 15.0).Within(Tol));
            Assert.That(stats.ExpectedMatches, Is.EqualTo(13.0).Within(Tol));
        });
    }

    // B3 — CA x6 with index 6 = T: copies CA|CA|CA|TA|CA|CA -> 8 matches, 2 mismatches; PM = 0.80.
    [Test]
    public void ComputeBernoulliStatistics_PmEqualToDefault_MeetsThreshold()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CACACATACACA", period: 2);

        Assert.Multiple(() =>
        {
            Assert.That(stats.BernoulliTrials, Is.EqualTo(10));
            Assert.That(stats.Matches, Is.EqualTo(8));
            Assert.That(stats.Mismatches, Is.EqualTo(2));
            Assert.That(stats.MatchProbability, Is.EqualTo(0.80).Within(Tol));
            Assert.That(stats.MeetsExpectedMatchProbability, Is.True, "0.80 >= 0.80");
        });
    }

    // B4 — ACACTGTG, period 4: the WDP aligns only one copy (TGTG), so no two copies can be compared.
    [Test]
    public void ComputeBernoulliStatistics_NoTwoAlignedCopies_ZeroTrials()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("ACACTGTG", period: 4);

        Assert.Multiple(() =>
        {
            Assert.That(stats.AdjacentCopyPairs, Is.EqualTo(0));
            Assert.That(stats.BernoulliTrials, Is.EqualTo(0));
            Assert.That(stats.MatchProbability, Is.EqualTo(0.0));
            Assert.That(stats.MeetsExpectedMatchProbability, Is.False);
        });
    }

    // B5 — caller-supplied PM threshold on B2 (PM 0.8667).
    [Test]
    public void ComputeBernoulliStatistics_CustomExpectedPm_FlagsAgainstThatThreshold()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGTAGCAGCAG", 3, 0.80).MeetsExpectedMatchProbability, Is.True);
            Assert.That(RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGTAGCAGCAG", 3, 0.90).MeetsExpectedMatchProbability, Is.False);
        });
    }

    // B6 — deletion tract: TRF m=25 mm=0 ind=2 (92 % / 7 %), 10.0 copies. The OLD implementation cut the
    // tract at fixed period boundaries, so a single deletion shifted the frame of every later copy.
    [Test]
    public void ComputeBernoulliStatistics_DeletionTract_MatchesTrfAdjacentCopyCounts()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", period: 3);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Matches, Is.EqualTo(25));
            Assert.That(stats.Mismatches, Is.EqualTo(0));
            Assert.That(stats.Indels, Is.EqualTo(2));
            Assert.That(stats.BernoulliTrials, Is.EqualTo(27));
            Assert.That(stats.AdjacentCopyPairs, Is.EqualTo(9));
            Assert.That(stats.MatchProbability, Is.EqualTo(25.0 / 27).Within(Tol));
            Assert.That(stats.IndelProbability, Is.EqualTo(2.0 / 27).Within(Tol));
        });
    }

    // B7 — PM + mismatch fraction + PI = 1.
    [Test]
    public void ComputeBernoulliStatistics_Probabilities_PartitionTheTrials()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", period: 3);

        double mismatchFraction = (double)stats.Mismatches / stats.BernoulliTrials;
        Assert.That(stats.MatchProbability + mismatchFraction + stats.IndelProbability, Is.EqualTo(1.0).Within(Tol));
    }

    // B8/B9 — validation.
    [Test]
    public void ComputeBernoulliStatistics_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RepeatFinder.ComputeBernoulliStatistics("CAG", period: 3));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.ComputeBernoulliStatistics(null!, 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.ComputeBernoulliStatistics("CACA", 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.ComputeBernoulliStatistics("CACA", 2001));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.ComputeBernoulliStatistics("CACA", 2, 1.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.ComputeBernoulliStatistics("CACA", 2, double.NaN));
        });
    }

    // B10 — Benson (1999) defaults PM = .80, PI = .10.
    [Test]
    public void ComputeBernoulliStatistics_DocumentedDefaults_MatchBenson1999()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.TrfDefaultMatchProbability, Is.EqualTo(0.80).Within(Tol));
            Assert.That(RepeatFinder.TrfDefaultIndelProbability, Is.EqualTo(0.10).Within(Tol));
        });
    }

    // The Bernoulli statistics of a detected repeat equal the detector's reported percentages.
    [Test]
    public void ComputeBernoulliStatistics_OnDetectedRegion_EqualsReportedPercentages()
    {
        const string seq = "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG";
        var r = RepeatFinder.FindApproximateTandemRepeats(seq, 3, 3).Single();
        var stats = RepeatFinder.ComputeBernoulliStatistics(seq.Substring(r.Start, r.SpanLength), r.Period);

        Assert.Multiple(() =>
        {
            Assert.That(stats.PercentMatches, Is.EqualTo(r.PercentMatches).Within(Tol));
            Assert.That(stats.PercentIndels, Is.EqualTo(r.PercentIndels).Within(Tol));
        });
    }

    #endregion
}
