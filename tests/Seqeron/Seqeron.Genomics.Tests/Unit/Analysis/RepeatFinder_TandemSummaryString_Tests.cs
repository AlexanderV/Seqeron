namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// REP-TANDEM-001 raw-string (N/IUPAC-tolerant) <c>GetTandemRepeatSummary</c> overloads (B04 completeness audit WP11).
/// Expected counts are copied from misa.pl v1.0 <c>.statistics</c> runs (perl 5.38; "Total size of examined sequences",
/// "Total number of identified SSRs", "Distribution to different repeat type classes"); misa.pl scans the raw sequence
/// with <c>[acgt]</c> (case-insensitive), so N/IUPAC symbols never belong to an SSR and <c>length $seq</c> counts them.
/// </summary>
[TestFixture]
public class RepeatFinder_TandemSummaryString_Tests
{
    /// <summary>misa.pl (default 1-10 2-6 3-5 4-5 5-5 6-5): size 93, 6 SSRs, classes 1:3 2:1 3:1 4:1
    /// (.misa row "c* (T)11n(A)12(CG)6nnnnnn(CAG)5ry(GATA)5(A)15* 93 1 93"). Uniform 3 copies: identical counts.</summary>
    private const string S1 =
        "TTTTTTTTTTTNAAAAAAAAAAAACGCGCGCGCGCGNNNNNNCAGCAGCAGCAGCAGryGATAGATAGATAGATAGATAaaaaaAAaaaAAAA";

    /// <summary>misa.pl default: size 87, 3 SSRs, classes 1:1 2:2 ("c (CA)8gtn…n(AT)9ryswkm(G)15 77 11 87");
    /// uniform 3 copies (1-3 … 6-3): 4 SSRs, classes 1:1 2:3 ("c (AC)4an(CA)8gtn…").</summary>
    private const string S2 =
        "ACACACACANCACACACACACACACAGTNNNNNNNNNNNNNNNNNNNNATATATATATATATATATRYSWKMggggggggggggggg";

    private static readonly IReadOnlyDictionary<int, int> Uniform3 =
        Enumerable.Range(1, 6).ToDictionary(p => p, _ => 3);

    private static int[] Classes(TandemRepeatSummary s) =>
    [
        s.MononucleotideRepeats, s.DinucleotideRepeats, s.TrinucleotideRepeats,
        s.TetranucleotideRepeats, s.PentanucleotideRepeats, s.HexanucleotideRepeats,
    ];

    [Test]
    public void MisaRegex_StringWithNAndIupac_EqualsMisaPlStatistics()
    {
        var d1 = RepeatFinder.GetTandemRepeatSummary(S1, RepeatFinder.MisaDefaultMinRepeats, MicrosatelliteScanMode.MisaRegex);
        var d2 = RepeatFinder.GetTandemRepeatSummary(S2, RepeatFinder.MisaDefaultMinRepeats, MicrosatelliteScanMode.MisaRegex);
        var u1 = RepeatFinder.GetTandemRepeatSummary(S1, Uniform3, MicrosatelliteScanMode.MisaRegex);
        var u2 = RepeatFinder.GetTandemRepeatSummary(S2, Uniform3, MicrosatelliteScanMode.MisaRegex);
        Assert.Multiple(() =>
        {
            Assert.That(S1, Has.Length.EqualTo(93));
            Assert.That(S2, Has.Length.EqualTo(87));
            Assert.That(d1.TotalRepeats, Is.EqualTo(6));
            Assert.That(Classes(d1), Is.EqualTo(new[] { 3, 1, 1, 1, 0, 0 }));
            Assert.That(d2.TotalRepeats, Is.EqualTo(3));
            Assert.That(Classes(d2), Is.EqualTo(new[] { 1, 2, 0, 0, 0, 0 }));
            Assert.That(u1.TotalRepeats, Is.EqualTo(6));
            Assert.That(Classes(u1), Is.EqualTo(new[] { 3, 1, 1, 1, 0, 0 }));
            Assert.That(u2.TotalRepeats, Is.EqualTo(4));
            Assert.That(Classes(u2), Is.EqualTo(new[] { 1, 3, 0, 0, 0, 0 }));
        });
    }

    [Test]
    public void Percentage_UsesFullLengthIncludingN()
    {
        // misa.pl: (A)10 1-10, "Total size of examined sequences (bp): 20" → 10 of 20 bases covered.
        var s = RepeatFinder.GetTandemRepeatSummary("AAAAAAAAAANNNNNNNNNN", RepeatFinder.MisaDefaultMinRepeats);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(1));
            Assert.That(s.TotalRepeatBases, Is.EqualTo(10));
            Assert.That(s.PercentageOfSequence, Is.EqualTo(50.0));
        });
    }

    [Test]
    public void NInterruptsRun_BothConventions()
    {
        // CACACA | N | CACACA: N is not [acgt], so two (CA)3 runs (misa.pl 2-3: (CA)3 1-6, (CA)3 8-13); 12 of 13 bases.
        const string seq = "CACACANCACACA";
        var max = RepeatFinder.GetTandemRepeatSummary(seq, 3);
        var misa = RepeatFinder.GetTandemRepeatSummary(seq, Uniform3, MicrosatelliteScanMode.MisaRegex);
        Assert.Multiple(() =>
        {
            Assert.That(max.TotalRepeats, Is.EqualTo(2));
            Assert.That(max.DinucleotideRepeats, Is.EqualTo(2));
            Assert.That(max.PercentageOfSequence, Is.EqualTo(12.0 / 13 * 100).Within(1e-12));
            Assert.That(max.MostFrequentUnit, Is.EqualTo("CA"));
            Assert.That(misa, Is.EqualTo(max));
        });
    }

    [Test]
    public void AcgtOnlyInput_StringEqualsDnaSequenceOverloads()
    {
        const string seq = "aaaagaaaagaaaagaaaagaCCTTAGGGTTAGGGTTAGGGTTAGGGcagcagcagTTTTTTTTTTTT";
        var dna = new DnaSequence(seq);
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.GetTandemRepeatSummary(seq, 3), Is.EqualTo(RepeatFinder.GetTandemRepeatSummary(dna, 3)));
            Assert.That(RepeatFinder.GetTandemRepeatSummary(seq, RepeatFinder.MisaDefaultMinRepeats),
                Is.EqualTo(RepeatFinder.GetTandemRepeatSummary(dna, RepeatFinder.MisaDefaultMinRepeats)));
            Assert.That(RepeatFinder.GetTandemRepeatSummary(seq, Uniform3, MicrosatelliteScanMode.MisaRegex),
                Is.EqualTo(RepeatFinder.GetTandemRepeatSummary(dna, Uniform3, MicrosatelliteScanMode.MisaRegex)));
        });
    }

    [Test]
    public void MaximalRuns_ClassCounts_EqualSumOverAcgtPieces()
    {
        // Runs never cross a non-ACGT symbol, so the per-class counts of the whole string equal the sums over its
        // maximal ACGT pieces summarised as DnaSequence.
        foreach (var seq in new[] { S1, S2 })
        {
            var whole = RepeatFinder.GetTandemRepeatSummary(seq, 3);
            var pieces = System.Text.RegularExpressions.Regex.Split(seq.ToUpperInvariant(), "[^ACGT]+")
                .Where(p => p.Length > 0)
                .Select(p => RepeatFinder.GetTandemRepeatSummary(new DnaSequence(p), 3))
                .ToList();
            var sum = new int[6];
            foreach (var p in pieces)
            {
                var c = Classes(p);
                for (int k = 0; k < 6; k++) sum[k] += c[k];
            }

            Assert.Multiple(() =>
            {
                Assert.That(Classes(whole), Is.EqualTo(sum));
                Assert.That(whole.TotalRepeats, Is.EqualTo(pieces.Sum(p => p.TotalRepeats)));
                Assert.That(whole.TotalRepeatBases, Is.EqualTo(pieces.Sum(p => p.TotalRepeatBases)));
            });
        }
    }

    [Test]
    public void NullOrEmpty_YieldsEmptySummary_InvalidArgumentsThrow()
    {
        Assert.Multiple(() =>
        {
            foreach (var s in new[]
                     {
                         RepeatFinder.GetTandemRepeatSummary((string)null!, 3),
                         RepeatFinder.GetTandemRepeatSummary("", RepeatFinder.MisaDefaultMinRepeats),
                         RepeatFinder.GetTandemRepeatSummary("", Uniform3, MicrosatelliteScanMode.MisaRegex),
                         RepeatFinder.GetTandemRepeatSummary("NNNNNNNNNNNN", 2),
                     })
            {
                Assert.That(s.TotalRepeats, Is.Zero);
                Assert.That(s.PercentageOfSequence, Is.Zero);
                Assert.That(s.LongestRepeat, Is.Null);
                Assert.That(s.MostFrequentUnit, Is.Null);
            }

            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.GetTandemRepeatSummary("ACGT", 1));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.GetTandemRepeatSummary("ACGT", (IReadOnlyDictionary<int, int>)null!));
            Assert.Throws<ArgumentException>(() => RepeatFinder.GetTandemRepeatSummary("ACGT", new Dictionary<int, int>()));
            // Unit sizes above 6 are valid (misa.pl accepts any size in its def line; B04 F61); size 0 is not.
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.GetTandemRepeatSummary("ACGT", new Dictionary<int, int> { [0] = 3 }));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RepeatFinder.GetTandemRepeatSummary("ACGT", Uniform3, (MicrosatelliteScanMode)99));
        });
    }
}
