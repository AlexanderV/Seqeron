// REP-APPROX-001 — Tandem Repeats Finder parameter sets and outputs (B04 audit WP6: L8 entropy parity,
// L9 TRF parameters, L10 masked sequence / flanks / alignment rows).
// Evidence: docs/Evidence/REP-APPROX-001-Evidence.md
// TestSpec: tests/TestSpecs/REP-APPROX-001.md (section "TRF parameters and outputs", P1..P14)
// Sources: Benson G (1999) Nucleic Acids Res 27(2):573-580; TRF 4.10.0 README (parameters, -m, -f, -r, -l,
//          "Table Explanation", "Alignment Explanation").
//
// Every expected row was produced by the compiled TRF 4.10.0 binary (github.com/Benson-Genomics-Lab/TRF,
// commit 355c1f9): `trf craft.fa <Match> <Mismatch> <Delta> <PM> <PI> <Minscore> <MaxPeriod> -h -d -ngs`
// (.dat columns: start end period copies consensus-size %match %indel score A C G T entropy consensus
// repeat left-50-flank right-50-flank), `-m` for the masked file, `-f -d` (HTML alignment file) for the
// alignment rows, and `-r` for redundancy off. `-l` rows come from a TRF build whose only change is that the
// `-l` value is taken in bp instead of millions (so small caps can be exercised). TRF prints 1-based indices,
// copy numbers to one decimal (float) and truncates percentages; the asserts compare against those printed values.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RepeatFinder_TrfParameters_Tests
{
    // Crafted sequences (seeded random background + mutated tandem arrays), see Evidence §WP6.
    private const string U1 =
        "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA";
    private const string U2 =
        "TAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTATCGTCCCTGCTACCATAGTTCGCATCCGTCCCTTCTAACGTAGTTCGCAATCGTCCCTTCTAACCTAGTTCGCATCGTCCGTTCTAACTAAGTTCGCAATCGTCCCTTCTAACTAGCCCGCACCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTT";
    private const string U3 =
        "CACACACACACACACACGCACACACACACCGACACACACACACACACACACACACACACAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC";
    private const string U4 =
        "AAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTTCACTGCCCTCTTCGTGCCCTCTTCACTGCCCTCGTCACTGCCCTGTTCACTGCCCTCTTCCCTGCCCTCTTACACTGCCCTCTTCACTGCCCGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAG";
    private const string U5 =
        "ATCGCGTCGAATGAGGGAGTTAGTCCTCGT" + /* (CA) x 100 */ "CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA" + "TCCAGCTGGTAATTGTTTTACCGCTTGGGA";

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

    #region L8 — entropy parity with TRF for regions containing N

    // TRF 4.10.0 get_statistics: p_b = count_b / (all non-gap symbols of the repeat, N included), entropy =
    // |sum p_b log2 p_b| over A, C, G, T. U1 repeat 61..124 = 64 symbols: A 9, C 9, G 17, T 27, N 2.
    // TRF row: 61 124 7 9.1 7 92 0 110 14 14 26 42 1.83 TCATTGG.
    [Test]
    public void EntropyTrf_RegionWithN_MatchesTrfColumn_LegacyEntropyDoesNot()
    {
        var r = RepeatFinder.FindApproximateTandemRepeats(U1, TandemRepeatsFinderParameters.Recommended).Single();
        AssertTrfRow(r, 61, 124, 7, 9.1, 7, 92, 0, 110, 14, 14, 26, 42, 1.83, "TCATTGG");

        double p9 = 9.0 / 64, p17 = 17.0 / 64, p27 = 27.0 / 64;
        double expectedTrf = -(p27 * Math.Log2(p27) + p17 * Math.Log2(p17) + 2 * p9 * Math.Log2(p9));
        double q9 = 9.0 / 62, q17 = 17.0 / 62, q27 = 27.0 / 62;
        double acgtOnly = -(q27 * Math.Log2(q27) + q17 * Math.Log2(q17) + 2 * q9 * Math.Log2(q9));
        Assert.Multiple(() =>
        {
            Assert.That(r.EntropyTrf, Is.EqualTo(expectedTrf).Within(1e-12), "TRF definition (denominator 64)");
            Assert.That(r.EntropyTrf, Is.EqualTo(1.829258111162015).Within(1e-12));
            Assert.That(r.Entropy, Is.EqualTo(acgtOnly).Within(1e-12), "legacy value normalised over A/C/G/T (62)");
            Assert.That(Math.Round(r.Entropy, 2), Is.EqualTo(1.84), "legacy Entropy would print 1.84, TRF prints 1.83");
        });
    }

    [Test]
    public void EntropyTrf_PureAcgtRegion_EqualsLegacyEntropy()
    {
        // U4 row 41..136 (TRF 1.66) and U3 row 1..60 (TRF 1.18) contain only A/C/G/T.
        var u4 = RepeatFinder.FindApproximateTandemRepeats(U4, TandemRepeatsFinderParameters.Recommended).First();
        var u3 = RepeatFinder.FindApproximateTandemRepeats(U3, TandemRepeatsFinderParameters.Recommended).Single();
        Assert.Multiple(() =>
        {
            Assert.That(u4.EntropyTrf, Is.EqualTo(u4.Entropy).Within(1e-12));
            Assert.That(Math.Round(u4.EntropyTrf, 2), Is.EqualTo(1.66));
            Assert.That(u3.EntropyTrf, Is.EqualTo(u3.Entropy).Within(1e-12));
            Assert.That(Math.Round(u3.EntropyTrf, 2), Is.EqualTo(1.18));
        });
    }

    #endregion

    #region L9 — TRF parameters

    [Test]
    public void Recommended_IsTrfRecommendedCommandLine()
    {
        var p = TandemRepeatsFinderParameters.Recommended;
        Assert.Multiple(() =>
        {
            Assert.That((p.MatchWeight, p.MismatchPenalty, p.IndelPenalty), Is.EqualTo((2, 7, 7)));
            Assert.That((p.MatchProbability, p.IndelProbability, p.MinScore, p.MaxPeriod), Is.EqualTo((80, 10, 50, 500)));
            Assert.That(p.MaxRepeatLength, Is.EqualTo(2_000_000), "TRF -l default 2 (million)");
            Assert.That(p.EliminateRedundancy, Is.True);
            Assert.That(p.FlankLength, Is.EqualTo(0));
        });
    }

    // `trf craft.fa 2 7 7 80 10 50 500`: U4 → 41 136 12 8.0 12 86 4 140 7 48 11 32 1.66 TCTTCACTGCCC and
    // 43 136 49 2.0 47 91 4 152 (consensus of 47 > 20: TRF uses its narrow-band WDP there; row locked only
    // for indices/period/score).
    [Test]
    public void RecommendedParameters_U4_ReproduceTrfRows()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(U4, TandemRepeatsFinderParameters.Recommended).ToList();
        Assert.That(rows, Has.Count.EqualTo(2));
        AssertTrfRow(rows[0], 41, 136, 12, 8.0, 12, 86, 4, 140, 7, 48, 11, 32, 1.66, "TCTTCACTGCCC");
        Assert.That((rows[1].Start + 1, rows[1].Start + rows[1].SpanLength, rows[1].Period, rows[1].AlignmentScore),
            Is.EqualTo((43, 136, 49, 152)));
    }

    // Weights change the alignment: `trf 2 3 5 80 10 40 200` → U3 1 60 2 30.0 2 89 0 105 (recommended 2 7 7:
    // 90 6 95); U4 → 41 136 12 8.0 12 86 4 160; U2 → 51 173 25 5.0 24 82 6 178 ATCGTCCCTTCTAACTAGTTCGCA and
    // 51 173 49 2.5 50 86 2 196; U1 → 61 124 7 9.1 7 92 0 118.
    [Test]
    public void MismatchAndIndelWeights_2_3_5_ReproduceTrfRows()
    {
        var p = Set(2, 3, 5, 80, 10, 40, 200);
        var u3 = RepeatFinder.FindApproximateTandemRepeats(U3, p).Single();
        var u4 = RepeatFinder.FindApproximateTandemRepeats(U4, p).First();
        var u2 = RepeatFinder.FindApproximateTandemRepeats(U2, p).ToList();
        var u1 = RepeatFinder.FindApproximateTandemRepeats(U1, p).Single();
        AssertTrfRow(u3, 1, 60, 2, 30.0, 2, 89, 0, 105, 46, 50, 3, 0, 1.18, "CA");
        AssertTrfRow(u4, 41, 136, 12, 8.0, 12, 86, 4, 160, 7, 48, 11, 32, 1.66, "TCTTCACTGCCC");
        AssertTrfRow(u1, 61, 124, 7, 9.1, 7, 92, 0, 118, 14, 14, 26, 42, 1.83, "TCATTGG");
        Assert.That(u2, Has.Count.EqualTo(2));
        AssertTrfRow(u2[0], 51, 173, 25, 5.0, 24, 82, 6, 178, 19, 35, 14, 30, 1.92, "ATCGTCCCTTCTAACTAGTTCGCA");
        AssertTrfRow(u2[1], 51, 173, 49, 2.5, 50, 86, 2, 196, 19, 35, 14, 30, 1.92,
            "ATCGTCCCTTCTAACCTAGTTCGCATCCGTCCCTTCTAACGAAGTTCGCA");
    }

    // PM = 75 (tuple sizes 3/4/5/7, sumdata75): `trf 2 7 7 75 20 50 500` → U2 51 167 25 4.8 24 84 10 148
    // ATCGTCCCTTCTAACTAGTTCGCA; `trf 2 5 5 75 10 30 100` → U2 51 167 25 4.8 24 84 10 168, U3 1 60 2 30.0 2 90 6 101,
    // U1 61 124 7 9.1 7 92 0 114.
    [Test]
    public void MatchProbability75_ReproducesTrfRows()
    {
        var p1 = Set(2, 7, 7, 75, 20, 50, 500);
        var p2 = Set(2, 5, 5, 75, 10, 30, 100);
        var a = RepeatFinder.FindApproximateTandemRepeats(U2, p1).Single(r => r.Period == 25);
        var b = RepeatFinder.FindApproximateTandemRepeats(U2, p2).Single(r => r.Period == 25);
        AssertTrfRow(a, 51, 167, 25, 4.8, 24, 84, 10, 148, 19, 34, 14, 31, 1.92, "ATCGTCCCTTCTAACTAGTTCGCA");
        AssertTrfRow(b, 51, 167, 25, 4.8, 24, 84, 10, 168, 19, 34, 14, 31, 1.92, "ATCGTCCCTTCTAACTAGTTCGCA");
        AssertTrfRow(RepeatFinder.FindApproximateTandemRepeats(U3, p2).Single(), 1, 60, 2, 30.0, 2, 90, 6, 101, 46, 50, 3, 0, 1.18, "CA");
        AssertTrfRow(RepeatFinder.FindApproximateTandemRepeats(U1, p2).Single(), 61, 124, 7, 9.1, 7, 92, 0, 114, 14, 14, 26, 42, 1.83, "TCATTGG");
    }

    // TRF 4.10.0 `sumdata75` table (PM = 75): the exact-moment derivation with tuple sizes 3/4/5/7 and floor
    // max(k + 1, 5) reproduces all d = 1..2000. (README text example "PM = .75, k = 5, d = 100 → 26"; the table
    // TRF actually uses, and this derivation, give 27.)
    [TestCase(1, 5)]
    [TestCase(18, 5)]
    [TestCase(20, 5)]
    [TestCase(21, 6)]
    [TestCase(29, 10)]
    [TestCase(30, 6)]
    [TestCase(43, 11)]
    [TestCase(44, 7)]
    [TestCase(100, 27)]
    [TestCase(159, 50)]
    [TestCase(160, 24)]
    [TestCase(500, 116)]
    [TestCase(2000, 567)]
    public void TrfSumOfHeadsCriterion_Pm75_MatchesTrfTable(int d, int expected)
    {
        Assert.That(RepeatFinder.TrfSumOfHeadsCriterion(d, 75), Is.EqualTo(expected));
    }

    [TestCase(1, 5)]
    [TestCase(30, 6)]
    [TestCase(160, 43)]
    [TestCase(2000, 818)]
    public void TrfSumOfHeadsCriterion_Pm80Overload_EqualsDefault(int d, int expected)
    {
        Assert.That(RepeatFinder.TrfSumOfHeadsCriterion(d, 80), Is.EqualTo(expected));
        Assert.That(RepeatFinder.TrfSumOfHeadsCriterion(d), Is.EqualTo(expected));
    }

    // `-r` (no redundancy elimination), recommended weights: U1 reports periods 7, 14, 21 over 61..124 (all score
    // 110); U3 periods 2, 4, 6 over 1..60; U5 periods 2, 4, 6 over 31..230 (copies 100.0 / 50.0 / 33.3, score 400).
    [Test]
    public void EliminateRedundancyOff_ReportsTrfMultiples()
    {
        var p = TandemRepeatsFinderParameters.Recommended with { EliminateRedundancy = false };
        var u1 = RepeatFinder.FindApproximateTandemRepeats(U1, p).ToList();
        var u3 = RepeatFinder.FindApproximateTandemRepeats(U3, p).ToList();
        var u5 = RepeatFinder.FindApproximateTandemRepeats(U5, p).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(u1.Select(r => (r.Start + 1, r.Start + r.SpanLength, r.Period, r.AlignmentScore)),
                Is.EqualTo(new[] { (61, 124, 7, 110), (61, 124, 14, 110), (61, 124, 21, 110) }));
            Assert.That(u3.Select(r => (r.Start + 1, r.Start + r.SpanLength, r.Period, r.AlignmentScore)),
                Is.EqualTo(new[] { (1, 60, 2, 95), (1, 60, 4, 95), (1, 60, 6, 95) }));
            Assert.That(u5.Select(r => (r.Period, Math.Round(r.CopyNumber, 1), r.AlignmentScore)),
                Is.EqualTo(new[] { (2, 100.0, 400), (4, 50.0, 400), (6, 33.3, 400) }));
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(U5, TandemRepeatsFinderParameters.Recommended).Single().Period,
                Is.EqualTo(2), "with redundancy elimination only the period-2 row remains (TRF)");
        });
    }

    // `-l` caps the wraparound alignment at that many rows on each side of the candidate (TRF newwrap: backward
    // and forward loops stop at maxwraplength). Instrumented TRF (-l in bp): -l 150 → U5 79 228 2 75.0 2 100 0 300;
    // -l 60 → U5 169 228 2 30.0 2 100 0 120 and U1 61 117 7 8.1 7 92 0 96; default → U5 31 230 2 100.0 score 400.
    [Test]
    public void MaxRepeatLength_CapsAlignmentRows_LikeTrf()
    {
        var p150 = TandemRepeatsFinderParameters.Recommended with { MaxRepeatLength = 150 };
        var p60 = TandemRepeatsFinderParameters.Recommended with { MaxRepeatLength = 60 };
        AssertTrfRow(RepeatFinder.FindApproximateTandemRepeats(U5, p150).Single(), 79, 228, 2, 75.0, 2, 100, 0, 300, 50, 50, 0, 0, 1.00, "CA");
        AssertTrfRow(RepeatFinder.FindApproximateTandemRepeats(U5, p60).Single(), 169, 228, 2, 30.0, 2, 100, 0, 120, 50, 50, 0, 0, 1.00, "CA");
        AssertTrfRow(RepeatFinder.FindApproximateTandemRepeats(U1, p60).Single(), 61, 117, 7, 8.1, 7, 92, 0, 96, 14, 14, 26, 42, 1.83, "TCATTGG");
        AssertTrfRow(RepeatFinder.FindApproximateTandemRepeats(U5, TandemRepeatsFinderParameters.Recommended).Single(),
            31, 230, 2, 100.0, 2, 100, 0, 400, 50, 50, 0, 0, 1.00, "CA");
    }

    // Validation mirrors TRF: weights / Minscore / PI positive; PM only 80 or 75 ("No sum table file for PM=..."
    // otherwise); MaxPeriod 1..2000; -l >= 1; flank length >= 0; minPeriod >= 1.
    [Test]
    public void InvalidParameters_Throw()
    {
        var ok = TandemRepeatsFinderParameters.Recommended;
        Assert.Multiple(() =>
        {
            foreach (var bad in new[]
            {
                ok with { MatchWeight = 0 }, ok with { MismatchPenalty = 0 }, ok with { IndelPenalty = -1 },
                ok with { MatchProbability = 70 }, ok with { MatchProbability = 85 }, ok with { IndelProbability = 0 },
                ok with { IndelProbability = 101 }, ok with { MinScore = 0 }, ok with { MaxPeriod = 0 },
                ok with { MaxPeriod = 2001 }, ok with { MaxRepeatLength = 0 }, ok with { FlankLength = -1 },
            })
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats("CACACACA", bad), bad.ToString());
            }
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats("CACACACA", ok, 0));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindApproximateTandemRepeats("CACACACA", (TandemRepeatsFinderParameters)null!));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindApproximateTandemRepeats((DnaSequence)null!, ok));
            Assert.That(RepeatFinder.FindApproximateTandemRepeats("", ok), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats((string)null!, ok), Is.Empty);
        });
    }

    [Test]
    public void DnaSequenceOverload_EqualsStringOverload()
    {
        var p = TandemRepeatsFinderParameters.Recommended with { FlankLength = 20 };
        var fromString = RepeatFinder.FindApproximateTandemRepeats(U4, p).ToList();
        var fromDna = RepeatFinder.FindApproximateTandemRepeats(new DnaSequence(U4), p).ToList();
        Assert.That(fromDna, Is.EqualTo(fromString));
    }

    #endregion

    #region L10 — masked sequence, flanks, alignment rows

    // `trf craft.fa 2 7 7 80 10 50 500 -m`: every position of a reported repeat → N.
    [TestCase(U1, "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGC" + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" + "AGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA")]
    [TestCase(U3, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" + "AGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC")]
    [TestCase(U4, "AAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATA" + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" + "GAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAG")]
    public void MaskApproximateTandemRepeats_EqualsTrfMaskedFile(string sequence, string trfMask)
    {
        Assert.That(RepeatFinder.MaskApproximateTandemRepeats(sequence), Is.EqualTo(trfMask));
    }

    [Test]
    public void MaskApproximateTandemRepeats_SoftMask_LowercasesExactlyTheHardMaskedPositions()
    {
        string hard = RepeatFinder.MaskApproximateTandemRepeats(U2);
        string soft = RepeatFinder.MaskApproximateTandemRepeats(U2, softMask: true);
        Assert.Multiple(() =>
        {
            Assert.That(hard, Is.EqualTo(U2[..50] + new string('N', 117) + U2[167..]), "TRF -m: 51..167 masked (union of both rows)");
            Assert.That(soft, Is.EqualTo(U2[..50] + U2[50..167].ToLowerInvariant() + U2[167..]));
        });
    }

    [Test]
    public void MaskApproximateTandemRepeats_ResultsOverload_KeepsInputCase_AndValidates()
    {
        string input = "acgt" + string.Concat(Enumerable.Repeat("CA", 13)) + "tt";
        var repeats = new[] { new ApproximateTandemRepeatResult(4, 26, 2, 2, "CA", 13, 100, 0, 52) };
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.MaskApproximateTandemRepeats(input, repeats), Is.EqualTo("acgt" + new string('N', 26) + "tt"));
            Assert.That(RepeatFinder.MaskApproximateTandemRepeats(input, repeats, softMask: true), Is.EqualTo("acgt" + string.Concat(Enumerable.Repeat("ca", 13)) + "tt"));
            Assert.That(RepeatFinder.MaskApproximateTandemRepeats(input, Array.Empty<ApproximateTandemRepeatResult>()), Is.EqualTo(input));
            Assert.That(RepeatFinder.MaskApproximateTandemRepeats(""), Is.EqualTo(""));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.MaskApproximateTandemRepeats("ACGT", repeats));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.MaskApproximateTandemRepeats(null!));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.MaskApproximateTandemRepeats("ACGT", (IEnumerable<ApproximateTandemRepeatResult>)null!));
        });
    }

    // `-ngs` appends 50-bp flanks ('.' when the repeat touches the sequence end); `-f` prints 500-bp flanks in the
    // alignment file ("Left flanking sequence: None" at the start).
    [Test]
    public void FlankLength_ReproducesTrfFlanks()
    {
        var p50 = TandemRepeatsFinderParameters.Recommended with { FlankLength = 50 };
        var p500 = TandemRepeatsFinderParameters.Recommended with { FlankLength = 500 };
        var u1 = RepeatFinder.FindApproximateTandemRepeats(U1, p50).Single();
        var u3 = RepeatFinder.FindApproximateTandemRepeats(U3, p50).Single();
        var u3Long = RepeatFinder.FindApproximateTandemRepeats(U3, p500).Single();
        Assert.Multiple(() =>
        {
            Assert.That(u1.LeftFlank, Is.EqualTo("CCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGC"));
            Assert.That(u1.RightFlank, Is.EqualTo("AGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCG"));
            Assert.That(u3.LeftFlank, Is.Empty, "TRF '.' / None: the repeat starts at position 1");
            Assert.That(u3.RightFlank, Is.EqualTo("AGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTA"));
            Assert.That(u3Long.RightFlank, Is.EqualTo(U3[60..]), "500-bp flank truncated at the sequence end (Indices 61 -- 140)");
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(U3, TandemRepeatsFinderParameters.Recommended).Single().LeftFlank,
                Is.Null, "no flanks unless requested");
        });
    }

    // TRF alignment file (`-f -d`), U4 repeat 41--136: pairs (sequence / consensus) concatenated.
    [Test]
    public void AlignmentRows_ReproduceTrfAlignmentFile()
    {
        var r = RepeatFinder.FindApproximateTandemRepeats(U4, TandemRepeatsFinderParameters.Recommended).First();
        const string trfTop =
            "TCTTCACTGCCC" + "TCTTC-GTGCCC" + "TCTTCACTGCCC" + "TCGTCACTGCCC" + "TGTTCACTGCCC" + "TCTTCCCTGCCC" + "TCTTACACTGCCC" + "TCTTCACTGCCC";
        const string trfBottom =
            "TCTTCACTGCCC" + "TCTTCACTGCCC" + "TCTTCACTGCCC" + "TCTTCACTGCCC" + "TCTTCACTGCCC" + "TCTTCACTGCCC" + "TCTT-CACTGCCC" + "TCTTCACTGCCC";
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequence, Is.EqualTo(trfTop));
            Assert.That(r.AlignedConsensus, Is.EqualTo(trfBottom));
            Assert.That(r.AlignedSequence!.Replace("-", ""), Is.EqualTo(U4.Substring(r.Start, r.SpanLength)));
        });
    }

    [Test]
    public void AlignmentRows_AlsoOnLegacyOverload_AndConsistentWithStatistics()
    {
        var r = RepeatFinder.FindApproximateTandemRepeats(U1, 1, 500).Single();
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequence!.Length, Is.EqualTo(r.AlignedConsensus!.Length));
            Assert.That(r.AlignedSequence.Replace("-", ""), Is.EqualTo(U1.Substring(r.Start, r.SpanLength)));
            Assert.That(r.AlignedConsensus.Replace("-", ""), Does.StartWith(r.Consensus), "consensus copies start at the reported phase");
            Assert.That(r.LeftFlank, Is.Null);
        });
    }

    #endregion
}
