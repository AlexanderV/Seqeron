// REP-APPROX-001 — caller-supplied apparent-size table and TRF 4.10.0 output formatters (B04 audit WP14).
// Evidence: docs/Evidence/REP-APPROX-001-Evidence.md (§WP14); TestSpec: tests/TestSpecs/REP-APPROX-001.md (O1..O12).
// Sources: Benson G (1999) Nucleic Acids Res 27(2):573-580; TRF 4.10.0 README ("Apparent Size Distribution", "Data
//          file", "-ngs", "Table Explanation"); TRF 4.10.0 source read for the output layout only (trfrun.h .dat / -ngs /
//          summary writers, trfclean.h OutputHTML, tr30dat.c get_statistics: IL field types and truncations).
//
// Every expected text below was written by the compiled TRF 4.10.0 binary (github.com/Benson-Genomics-Lab/TRF, commit
// 355c1f9): `trf <file> <params> -d -h` (.dat), `-ngs -h` (stdout) and without -h (.1.html / .summary.html). TRF's own
// apparent-size (waitdata) tables are AGPL source data and are NOT embedded: the tests use the library's exact table,
// synthetic tables, and single entries whose TRF value is documented by the WP7 bisection (B04 F46).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RepeatFinder_TrfOutput_Tests
{
    // `trf one.fa 2 7 7 80 10 50 500`, FASTA header ">seq1 test".
    private const string Seq1 =
        "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTAC"
        + "CCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTC"
        + "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAGATG"
        + "TCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCTGGCCGCCGAGACACA"
        + "CACACACACACACACACACACACACACACACACACACACACACACACACACACACATAGCTGAGCGGCGAACCACTAGAAAAGGTTCAGACCCCGGAGCC"
        + "CAGCCGTCACGATTGT";

    // set700 sequence s342 (WP7 bisection: TRF's simulated table differs from the exact one only at d = 43 for this case).
    private const string S342 =
        "CGACGGCTGTGGAACGTCAAACAACACATACCTGTACTACAGCAGGGGAGGTTGGTTACACTCGATTGACGTATCCGGGGGTTTCGGCCTCAGACGCTAC"
        + "CCACATATACGTACGGGGTACGGGAGCCGGGACGCAATGCCGGCCGTCCCGCTTACGGACGTAGCCAACAACCCGCTCAAGATTCGGTCTTTATTCTTTC"
        + "TCTCGGTAACCGAATGGACTTAATGCGGTTCTACCCGGAACGCTTGTCCATGCAGGGCCAACCCTGTTACACTATTCCGGGGTACGATTGAACGACTAAT"
        + "CGTTGATCTACTTTACTTAGTAAGCAACCCGCCCGGTAAGCTTAGCACTTAGTATGAGTTGATGCCAAATGAGCGGGTGAAGTTAGAGATCATCTGGCCA"
        + "TTGCACTCCTATATATTGCACACAATTCCTAAAAATAGCTCACTTAGAAAACCACGACTGTGAGGAAATGTTACTAGAAAGTCATCTCCTAATTAAAATA"
        + "AGGGCCATGTTGCCACGCNNNNCAGGTAGGTGTCGGACGCTAGTTGTTGGATCATGCAGACTGCTAGCTGAAGTACCAGAAGCCTGAGCAGAGGGTGATA"
        + "CATGTCCCCCTGAGATCGCGTTGGGGCCATTAGCTACACATATGGAACCGGTCCCGGCAGAATAAGGCAACTTACTTTCCCGCCAGCGGCTTCGGGAATC"
        + "ATTCAGGTTCTGGTGCGCCCGCCAGCGGCTCCGGGAAGTCATTCAAAGGTTCTGNTCGCCCAGCCAGCGGCCTCCGGAATGTCATTCAAGTTCTGGTGCG"
        + "CCCCAGCGGCTCGGGAAGTCATNCAAGTTCTGAGCATTTGACTCGATCACAAAAAAAAAATCGGACCGGTTCTGGGAAGGCACATGCCTTGATCCCCCGA"
        + "ACACTATCAGTCAGGAAGATTGCATATTA";

    // set700 sequence s350 (WP7 bisection: single differing entry d = 24, PM 75).
    private const string S350 =
        "AGCGGTTTGACGTTCGCCTCGGCCTCCGCGGGTTGCTGCTGAAGACATTTGTCGTATGCAGTGAGCTGTCACGGATCACAATCAGTTACTCTGGCAGNNC"
        + "TTATGCCTCTGCACTGACTCCTCTATGCAAATGTTTGAAGCGTTGCTTGCAAGCAATGGAAGTCTTCGCGGCACATGAAACGTCGCACGTGCCTTCTCGT"
        + "GTGGTTTTCCGAGCATTGATTTGTATCAAGGTAGCAGAAGCCCTGGGCCGNGCCCNAGGCTGGGCCTAACGCCCTCAGCCGCCCGTGTGTCTCATGCACT"
        + "CTGAGGTNGCAGAAGCCCNGGGCCGTNCCCAAGGCTGGGCCTAACGCCCTCAGCCGCCCGTGTGTCTCATGCACTCTGAGGTNGCAGAAGCCCTGGGCCG"
        + "TGCCCAAGGCTGGGCCTAACGCCCTCAGCCGCCCGTGTGTCTCATGCACTCTGAGGTAGCAGAAGCCCTGGGCCGTGCCCAAGGCTGGGCCTAACGCCCT"
        + "CAGCCGCCCNTGTGTCTCATGCACTCTGAGGTAGCAGAAGCCCTGGGCCGTGCCCAAGGCTGGGCCTAACGCCCTCTGAGTTTGGTGGTATGATCATTGC"
        + "ATCGGCTGCTCTGATTAACGGCCCTCCTTGGTGTACTGCGGTTGAGCCACATCCAATTAGCTGCATAAGAGCTTAATCGCACTTTCTTCCGTATCTGTTC"
        + "TCAGGTCTACTTNCGTCATATCGATTTATACTGCCTACGTATCTGTTCTCAGGTCTACTTTCGTCATANCGACTTATACTGCCTACGTATCTGTTCTCAG"
        + "GTGTACTGCCATTGGTCTGGATCCCTAGGGACTCCTTGAATCCCGATGATGTATAGCGTAACCATTTTCGGATGATTCCCCCAACACCCAGTGTGCTNAG"
        + "GCAACACGATTCGCCCCATACCGGGCCGATTGGTACCGAGTCGGCTAAGGTAACACGATTCGCCCCATTCCTGGCACAATGATACCCAGTCTGCTCCGGT"
        + "AACACNATTNGCCCCATACTGGGCNGGTTGGTGCCGAGTCTGCTCATGTAACACGANTCGCCCCANACCGGGCCGATTGGTNCGGTCGCCCGGCCATACG"
        + "ATCTTCTCTTAAA";

    // RepeatFinder_TrfDetection_Tests D2/D3 sequences: TRF reports nothing; the apparent-size test rejects a period-28
    // (148..209, score 92) and a period-22 (129..182, score 67) candidate whose tuple matches cluster at the window's right end.
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

    private const string U1 =
        "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA";

    private static TandemRepeatsFinderParameters Set(int match, int mismatch, int delta, int pm, int pi, int minScore, int maxPeriod) =>
        new()
        {
            MatchWeight = match, MismatchPenalty = mismatch, IndelPenalty = delta,
            MatchProbability = pm, IndelProbability = pi, MinScore = minScore, MaxPeriod = maxPeriod,
        };

    private static int[] ExactWith(int pm, int d, int y)
    {
        var table = TandemRepeatsFinderParameters.ExactApparentSizeTable(pm);
        table[d] = y;
        return table;
    }

    #region O1-O6 — ApparentSizeTable

    [Test]
    public void ExactTable_HasTrfLayoutAndReadmeValue()
    {
        var t75 = TandemRepeatsFinderParameters.ExactApparentSizeTable(75);
        var t80 = TandemRepeatsFinderParameters.ExactApparentSizeTable(80);
        Assert.Multiple(() =>
        {
            Assert.That(t75, Has.Length.EqualTo(2001));
            Assert.That(t80, Has.Length.EqualTo(TandemRepeatsFinderParameters.ApparentSizeTableLength));
            Assert.That(t75[100], Is.EqualTo(56), "TRF README: PM .75, k 5, d 100 -> 56");
            Assert.That(t75[0], Is.EqualTo(0));
            Assert.Throws<ArgumentOutOfRangeException>(() => TandemRepeatsFinderParameters.ExactApparentSizeTable(70));
        });
        t80[43] = -5; // a fresh copy each call
        Assert.That(TandemRepeatsFinderParameters.ExactApparentSizeTable(80)[43], Is.Not.EqualTo(-5));
    }

    // Supplying the exact table is the same as supplying none (null = exact computed table).
    [TestCase(80)]
    [TestCase(75)]
    public void SuppliedExactTable_EqualsDefault(int pm)
    {
        var p = Set(2, 5, 5, pm, 10, 30, 100);
        var withTable = p with { ApparentSizeTable = TandemRepeatsFinderParameters.ExactApparentSizeTable(pm) };
        foreach (var seq in new[] { S342, S350, Seq1, U1 })
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(seq, withTable), Is.EqualTo(RepeatFinder.FindApproximateTandemRepeats(seq, p)));
    }

    // WP7 bisection (B04 F46): for s342 TRF's simulated entry at d = 43 is y = 13 (waitdata80[43] = 29) where the
    // exact value is 12. With that one entry the library prints TRF's -ngs block byte for byte; with the exact table
    // it additionally reports a period-42 row TRF rejects.
    [Test]
    public void SingleTrfEntry_Pm80_D43_ReproducesTrfNgsBlock()
    {
        const string trf = "@s342\n"
        + "679 832 39 3.9 40 81 14 185 16 31 29 20 1.94 CCCGCCAGCGGCTCCGGGAAGTCATTCAAGTTCTGGTGCG CCCGCCAGCGGCTTCGGGAATCATTCAGGTTCTGGTGCGCCCGCCAGCGGCTCCGGGAAGTCATTCAAAGGTTCTGNTCGCCCAGCCAGCGGCCTCCGGAATGTCATTCAAGTTCTGGTGCGCCCCAGCGGCTCGGGAAGTCATNCAAGTTCTG ATTAGCTACACATATGGAACCGGTCCCGGCAGAATAAGGCAACTTACTTT AGCATTTGACTCGATCACAAAAAAAAAATCGGACCGGTTCTGGGAAGGCA\n";
        var p = Set(2, 7, 7, 80, 10, 50, 500);
        Assert.That(TandemRepeatsFinderParameters.ExactApparentSizeTable(80)[43], Is.EqualTo(12));
        var tuned = p with { ApparentSizeTable = ExactWith(80, 43, 13) };
        Assert.That(RepeatFinder.FormatTrfDatLines(S342, RepeatFinder.FindApproximateTandemRepeats(S342, tuned), "s342", tuned, TrfDatLayout.Ngs),
            Is.EqualTo(trf));

        var exact = RepeatFinder.FindApproximateTandemRepeats(S342, p).ToList();
        Assert.That(exact.Select(r => (r.Start + 1, r.Period)), Is.EqualTo(new[] { (679, 42), (679, 39) }));
    }

    // As above for PM 75 (s350, d = 24: TRF y = 8 = 24 - waitdata75[24] - 1 with waitdata75[24] = 15; exact y = 7).
    [Test]
    public void SingleTrfEntry_Pm75_D24_ReproducesTrfNgsBlock()
    {
        const string trf = "@s350\n"
        + "229 576 75 4.6 75 96 0 647 15 35 30 16 1.89 AGGTAGCAGAAGCCCTGGGCCGTGCCCAAGGCTGGGCCTAACGCCCTCAGCCGCCCGTGTGTCTCATGCACTCTG AGGTAGCAGAAGCCCTGGGCCGNGCCCNAGGCTGGGCCTAACGCCCTCAGCCGCCCGTGTGTCTCATGCACTCTGAGGTNGCAGAAGCCCNGGGCCGTNCCCAAGGCTGGGCCTAACGCCCTCAGCCGCCCGTGTGTCTCATGCACTCTGAGGTNGCAGAAGCCCTGGGCCGTGCCCAAGGCTGGGCCTAACGCCCTCAGCCGCCCGTGTGTCTCATGCACTCTGAGGTAGCAGAAGCCCTGGGCCGTGCCCAAGGCTGGGCCTAACGCCCTCAGCCGCCCNTGTGTCTCATGCACTCTGAGGTAGCAGAAGCCCTGGGCCGTGCCCAAGGCTGGGCCTAACGCCCTC AACGTCGCACGTGCCTTCTCGTGTGGTTTTCCGAGCATTGATTTGTATCA TGAGTTTGGTGGTATGATCATTGCATCGGCTGCTCTGATTAACGGCCCTC\n"
        + "719 791 24 3.2 21 56 32 52 20 26 13 38 1.90 TATCGATTTATACTGCCTACG TATCGATTTATACTGCCTACGTATCTGTTCTCAGGTCTACTTTCGTCATANCGACTTATACTGCCTACGTATC GAGCTTAATCGCACTTTCTTCCGTATCTGTTCTCAGGTCTACTTNCGTCA TGTTCTCAGGTGTACTGCCATTGGTCTGGATCCCTAGGGACTCCTTGAAT\n"
        + "690 807 48 2.5 48 94 0 208 17 25 16 38 1.90 CGTATCTGTTCTCAGGTCTACTTTCGTCATATCGACTTATACTGCCTA CGTATCTGTTCTCAGGTCTACTTNCGTCATATCGATTTATACTGCCTACGTATCTGTTCTCAGGTCTACTTTCGTCATANCGACTTATACTGCCTACGTATCTGTTCTCAGGTGTACT GGTTGAGCCACATCCAATTAGCTGCATAAGAGCTTAATCGCACTTTCTTC GCCATTGGTCTGGATCCCTAGGGACTCCTTGAATCCCGATGATGTATAGC\n"
        + "727 810 27 3.3 27 55 24 77 17 27 16 36 1.92 TATACTGCCTACGTATCTGTTCTCAGG TATACTGCCTACGTATCTGTTCTCAGGTCTACTTTCGTCATANCGACTTATACTGCCTACGTATCTGTTCTCAGGTGTACTGCC TCGCACTTTCTTCCGTATCTGTTCTCAGGTCTACTTNCGTCATATCGATT ATTGGTCTGGATCCCTAGGGACTCCTTGAATCCCGATGATGTATAGCGTA\n"
        + "886 1081 49 4.0 49 75 2 233 20 31 25 19 1.96 ACCGAGTCTGCTCAGGTAACACGATTCGCCCCATACCGGGCCGATTGGT ACCCAGTGTGCTNAGGCAACACGATTCGCCCCATACCGGGCCGATTGGTACCGAGTCGGCTAAGGTAACACGATTCGCCCCATTCCTGGCACAATGATACCCAGTCTGCTCCGGTAACACNATTNGCCCCATACTGGGCNGGTTGGTGCCGAGTCTGCTCATGTAACACGANTCGCCCCANACCGGGCCGATTGGT TTGAATCCCGATGATGTATAGCGTAACCATTTTCGGATGATTCCCCCAAC NCGGTCGCCCGGCCATACGATCTTCTCTTAAA\n"
        + "886 1065 98 1.8 98 82 0 262 21 32 23 20 1.96 ACCCAGTCTGCTCAGGCAACACGATTCGCCCCATACCGGGCCGATTGGTACCGAGTCGGCTAAGGTAACACGATTCGCCCCATTCCTGGCACAATGAT ACCCAGTGTGCTNAGGCAACACGATTCGCCCCATACCGGGCCGATTGGTACCGAGTCGGCTAAGGTAACACGATTCGCCCCATTCCTGGCACAATGATACCCAGTCTGCTCCGGTAACACNATTNGCCCCATACTGGGCNGGTTGGTGCCGAGTCTGCTCATGTAACACGANTCGCCCCA TTGAATCCCGATGATGTATAGCGTAACCATTTTCGGATGATTCCCCCAAC NACCGGGCCGATTGGTNCGGTCGCCCGGCCATACGATCTTCTCTTAAA\n";
        var p = Set(2, 5, 5, 75, 10, 30, 100);
        Assert.That(TandemRepeatsFinderParameters.ExactApparentSizeTable(75)[24], Is.EqualTo(7));
        var tuned = p with { ApparentSizeTable = TandemRepeatsFinderParameters.ApparentSizeTableFromWaitingTimes(
            WaitingTimesFrom(TandemRepeatsFinderParameters.ExactApparentSizeTable(75), 24, 15)) };
        Assert.That(RepeatFinder.FormatTrfDatLines(S350, RepeatFinder.FindApproximateTandemRepeats(S350, tuned), "s350", tuned, TrfDatLayout.Ngs),
            Is.EqualTo(trf));

        // Exact table: every TRF row but the period-24 one is reproduced; the period-24 row differs.
        string exact = RepeatFinder.FormatTrfDatLines(S350, RepeatFinder.FindApproximateTandemRepeats(S350, p), "s350", p, TrfDatLayout.Ngs);
        var trfRows = trf.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        var csRows = exact.Split('\n', StringSplitOptions.RemoveEmptyEntries);
        Assert.That(csRows.Except(trfRows).Select(l => l.Split(' ')[2]), Is.EqualTo(new[] { "24" }));
        Assert.That(trfRows.Except(csRows).Select(l => l.Split(' ')[2]), Is.EqualTo(new[] { "24" }));
    }

    /// <summary>Waiting times w(d) = max(d,20) − y(d) − 1 of an apparent-size table, with w(d) replaced at one d.</summary>
    private static int[] WaitingTimesFrom(int[] y, int d, int w)
    {
        var waits = new int[y.Length];
        for (int t = 1; t < y.Length; t++)
            waits[t] = Math.Max(t, 20) - y[t] - 1;
        waits[d] = w;
        return waits;
    }

    // Synthetic tables. y = 0 everywhere (most permissive: offset max(d,20) − 1, any earlier tuple match qualifies)
    // disables the apparent-size rejection, so the clustered period-28 / period-22 candidates TRF rejects are reported
    // (the D2/D3 ablation rows). y = max(d,20) − 1 (strictest: offset 0) still detects a perfect run, because the oldest
    // tuple-match run in the window starts before the window (TRF counts runs, not single tuples).
    [Test]
    public void SyntheticTables_ChangeDetectionPredictably()
    {
        var p = TandemRepeatsFinderParameters.Recommended;
        var zero = new int[TandemRepeatsFinderParameters.ApparentSizeTableLength];
        var strict = new int[TandemRepeatsFinderParameters.ApparentSizeTableLength];
        for (int d = 1; d < strict.Length; d++)
            strict[d] = Math.Max(d, 20) - 1;
        string micro = "GATTACAGGCTTAGCCATGCAAGT" + string.Concat(Enumerable.Repeat("AC", 30)) + "TTGCAGGCATTCGGATCAGTCGA";

        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(Spread28, p), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(Spread28, p with { ApparentSizeTable = zero })
                .Select(r => (r.Start + 1, r.Start + r.SpanLength, r.Period, r.AlignmentScore)), Is.EqualTo(new[] { (148, 209, 28, 92) }));
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(Spread22, p), Is.Empty);
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(Spread22, p with { ApparentSizeTable = zero })
                .Select(r => (r.Start + 1, r.Start + r.SpanLength, r.Period, r.AlignmentScore)), Is.EqualTo(new[] { (129, 182, 22, 67) }));
            // Same repeat; only the detection point (WP17 DetectionPosition, TRF "Found at i:") moves later.
            Assert.That(RepeatFinder.FindApproximateTandemRepeats(micro, p with { ApparentSizeTable = strict }).Select(r => r with { DetectionPosition = 0 }),
                Is.EqualTo(RepeatFinder.FindApproximateTandemRepeats(micro, p).Select(r => r with { DetectionPosition = 0 })));
        });
    }

    [Test]
    public void ApparentSizeTable_IsValidatedAndCopied()
    {
        var ok = TandemRepeatsFinderParameters.ExactApparentSizeTable(80);
        var p = TandemRepeatsFinderParameters.Recommended with { ApparentSizeTable = ok };
        ok[10] = 999; // the parameter set keeps its own copy
        Assert.Multiple(() =>
        {
            Assert.That(p.ApparentSizeTable![10], Is.EqualTo(TandemRepeatsFinderParameters.ExactApparentSizeTable(80)[10]));
            Assert.DoesNotThrow(p.Validate);
            Assert.Throws<ArgumentException>(() => (TandemRepeatsFinderParameters.Recommended with { ApparentSizeTable = new int[2000] }).Validate());
            var negative = new int[2001]; negative[5] = -1;
            Assert.Throws<ArgumentOutOfRangeException>(() => (TandemRepeatsFinderParameters.Recommended with { ApparentSizeTable = negative }).Validate());
            var tooLarge = new int[2001]; tooLarge[30] = 30; // must be <= max(d,20) - 1 = 29
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateTandemRepeats("ACACACAC", TandemRepeatsFinderParameters.Recommended with { ApparentSizeTable = tooLarge }));
            var edge = new int[2001]; edge[30] = 29; edge[5] = 19; edge[0] = -7; // entry 0 is ignored
            Assert.DoesNotThrow(() => (TandemRepeatsFinderParameters.Recommended with { ApparentSizeTable = edge }).Validate());
        });
    }

    [Test]
    public void WaitingTimes_ConvertToApparentSize()
    {
        var w = new int[2001];
        w[1] = 18; w[21] = 19; w[100] = 43; w[2000] = 68;
        var y = TandemRepeatsFinderParameters.ApparentSizeTableFromWaitingTimes(w);
        Assert.Multiple(() =>
        {
            Assert.That(y[0], Is.EqualTo(0));
            Assert.That(y[1], Is.EqualTo(20 - 18 - 1));
            Assert.That(y[21], Is.EqualTo(21 - 19 - 1));
            Assert.That(y[100], Is.EqualTo(56), "README example: waiting time 43 at d = 100 <-> apparent size 56");
            Assert.That(y[2000], Is.EqualTo(2000 - 68 - 1));
            Assert.That(y[50], Is.EqualTo(49));
            Assert.Throws<ArgumentException>(() => TandemRepeatsFinderParameters.ApparentSizeTableFromWaitingTimes(new int[10]));
            var bad = new int[2001]; bad[3] = 20;
            Assert.Throws<ArgumentOutOfRangeException>(() => TandemRepeatsFinderParameters.ApparentSizeTableFromWaitingTimes(bad));
            Assert.Throws<ArgumentNullException>(() => TandemRepeatsFinderParameters.ApparentSizeTableFromWaitingTimes(null!));
        });
    }

    #endregion

    #region O7-O12 — .dat / -ngs / HTML formatters

    // `trf one.fa 2 7 7 80 10 50 500 -d -h` → one.fa.2.7.7.80.10.50.500.dat, byte for byte.
    [Test]
    public void DatFile_MatchesTrfByteForByte()
    {
        const string trf = "Tandem Repeats Finder Program written by:\n"
        + "\n"
        + "Gary Benson\n"
        + "Program in Bioinformatics\n"
        + "Boston University\n"
        + "Version 4.10.0\n"
        + "\n"
        + "\n"
        + "Sequence: seq1 test\n"
        + "\n"
        + "\n"
        + "\n"
        + "Parameters: 2 7 7 80 10 50 500\n"
        + "\n"
        + "\n"
        + "201 296 8 12.0 8 100 0 192 25 25 25 25 2.00 ACGTTGCA ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCA\n"
        + "396 456 2 30.5 2 100 0 122 50 49 0 0 1.00 AC ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA\n";
        var p = TandemRepeatsFinderParameters.Recommended;
        var rows = RepeatFinder.FindApproximateTandemRepeats(Seq1, p);
        Assert.That(RepeatFinder.FormatTrfDatFileHeader() + RepeatFinder.FormatTrfDatLines(Seq1, rows, "seq1 test", p), Is.EqualTo(trf));
    }

    // `trf U1.fa 2 7 7 80 10 50 500 -ngs -h` (stdout): '@' line, row + 50-bp flanks; a soft-masked (lower-case)
    // input prints upper case as TRF does.
    [Test]
    public void NgsBlock_MatchesTrf_AndUpperCases()
    {
        const string trf = "@U1\n61 124 7 9.1 7 92 0 110 14 14 26 42 1.83 TCATTGG TCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGT CCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGC AGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCG\n";
        var p = TandemRepeatsFinderParameters.Recommended;
        Assert.That(RepeatFinder.FormatTrfDatLines(U1, RepeatFinder.FindApproximateTandemRepeats(U1, p), "U1", p, TrfDatLayout.Ngs), Is.EqualTo(trf));
        string lower = U1.ToLowerInvariant();
        Assert.That(RepeatFinder.FormatTrfDatLines(lower, RepeatFinder.FindApproximateTandemRepeats(lower, p), "U1", p, TrfDatLayout.Ngs), Is.EqualTo(trf));
    }

    // TRF -ngs flanks: '.' when the repeat touches a sequence end; shorter flanks near the ends. No repeats → empty
    // -ngs block (TRF prints no '@' line) and a bare Sequence/Parameters .dat block.
    [Test]
    public void NgsBlock_SequenceEndsAndEmpty()
    {
        string seq = string.Concat(Enumerable.Repeat("AC", 30)) + "GATTACAGGT";
        var p = TandemRepeatsFinderParameters.Recommended;
        var rows = RepeatFinder.FindApproximateTandemRepeats(seq, p).ToList();
        Assert.That(rows.Select(r => (r.Start, r.SpanLength)), Is.EqualTo(new[] { (0, 60) }));
        string line = RepeatFinder.FormatTrfDatLines(seq, rows, "x", p, TrfDatLayout.Ngs).Split('\n')[1];
        Assert.That(line, Does.EndWith(" . GATTACAGGT"));
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FormatTrfDatLines("GATTACA", [], "empty", p, TrfDatLayout.Ngs), Is.Empty);
            Assert.That(RepeatFinder.FormatTrfDatLines("GATTACA", [], "empty", p),
                Is.EqualTo("\n\nSequence: empty\n\n\n\nParameters: 2 7 7 80 10 50 500\n\n\n"));
        });
    }

    // TRF writes rows in report order (trfclean.h SortByCount), not by start: with TRF's d = 24 entry, s350's period-24 row
    // (found at a later position) precedes the period-48 row that starts earlier. OutputIndex carries TRF's OUTPUTcount: one.fa's anchors
    // are "201--296,8,12.0,8,1" and "396--456,2,30.5,2,4".
    [Test]
    public void Rows_FollowTrfReportOrder_AndOutputIndex()
    {
        var rows = RepeatFinder.FindApproximateTandemRepeats(Seq1, TandemRepeatsFinderParameters.Recommended).ToList();
        Assert.That(rows.Select(r => r.OutputIndex), Is.EqualTo(new[] { 1, 4 }));

        var p = Set(2, 5, 5, 75, 10, 30, 100) with { ApparentSizeTable = ExactWith(75, 24, 8) };
        var found = RepeatFinder.FindApproximateTandemRepeats(S350, p).ToList();
        Assert.That(found.Select(r => r.Start + 1).ToList(), Is.Ordered, "the API keeps start order");
        string text = RepeatFinder.FormatTrfDatLines(S350, found, "s350", p, TrfDatLayout.Ngs);
        var starts = text.Split('\n', StringSplitOptions.RemoveEmptyEntries).Skip(1).Select(l => int.Parse(l.Split(' ')[0])).ToList();
        Assert.That(starts, Is.EqualTo(new[] { 229, 719, 690, 727, 886, 886 }), "TRF's -ngs order");
        Assert.That(starts, Is.Not.Ordered);
    }

    // `trf U1.fa 2 7 7 80 10 50 500` → U1.fa.2.7.7.80.10.50.500.1.html, byte for byte.
    [Test]
    public void HtmlTable_MatchesTrfByteForByte()
    {
        const string trf = "<HTML><HEAD><TITLE>U1.fa.2.7.7.80.10.50.500.1.html</TITLE><BASE TARGET=\"U1.fa.2.7.7.80.10.50.500.txt.html\"></HEAD><BODY bgcolor=\"#FBF8BC\"><BR><PRE>Tandem Repeats Finder Program written by:</PRE><PRE><CENTER>Gary Benson<BR>Program in Bioinformatics<BR>Boston University<BR>Version 4.10.0<BR></CENTER>\n"
        + "Please cite:\n"
        + "G. Benson,\n"
        + "\"Tandem repeats finder: a program to analyze DNA sequences\"\n"
        + "Nucleic Acid Research(1999)\n"
        + "Vol. 27, No. 2, pp. 573-580.\n"
        + "\n"
        + "Sequence: U1\n"
        + "Parameters: 2 7 7 80 10 50 500\n"
        + "Length:  183</PRE>\n"
        + "\n"
        + "<P><PRE>Tables:   1   \n"
        + "\n"
        + "This is table  1  of  1  ( 1 repeats found )\n"
        + "</PRE><PRE>\n"
        + "Click on indices to view alignment\n"
        + "</PRE><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#table\" target = \"explanation\">Table Explanation</A><BR><BR>\n"
        + "<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n"
        + "<TR><TD WIDTH=140><CENTER>Indices</CENTER></TD><TD WIDTH=80><CENTER>Period<BR>Size </CENTER></TD><TD WIDTH=70><CENTER>Copy<BR>Number</CENTER></TD><TD WIDTH=70><CENTER>Consensus<BR>Size</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Matches</CENTER></TD><TD WIDTH=70><CENTER>Percent<BR>Indels</CENTER></TD><TD WIDTH=60><CENTER>Score</CENTER></TD><TD WIDTH=40><CENTER>A</CENTER></TD><TD WIDTH=40><CENTER>C</CENTER></TD><TD WIDTH=40><CENTER>G</CENTER></TD><TD WIDTH=40><CENTER>T</CENTER></TD><TD WIDTH=70><CENTER>Entropy<BR>(0-2)</CENTER></TD></TR>\n"
        + "<TR><TD><CENTER><A HREF=\"U1.fa.2.7.7.80.10.50.500.1.txt.html#61--124,7,9.1,7,1\">61--124</A></CENTER></TD><TD><CENTER>7</CENTER></TD><TD><CENTER>9.1</CENTER></TD><TD><CENTER>7</CENTER></TD><TD><CENTER>92</CENTER></TD><TD><CENTER>0</CENTER></TD><TD><CENTER>110</CENTER></TD><TD><CENTER>14</CENTER></TD><TD><CENTER>14</CENTER></TD><TD><CENTER>26</CENTER></TD><TD><CENTER>42</CENTER></TD><TD><CENTER>1.83</CENTER></TD></TR>\n"
        + "\n"
        + "</TABLE>\n"
        + "\n"
        + "<P><PRE>Tables:   1   \n"
        + "</PRE><P>The End!\n"
        + "\n"
        + "</BODY></HTML>\n";
        var p = TandemRepeatsFinderParameters.Recommended;
        var pages = RepeatFinder.FormatTrfHtmlTables(U1, RepeatFinder.FindApproximateTandemRepeats(U1, p), "U1", p, "U1.fa");
        Assert.That(pages, Has.Count.EqualTo(1));
        Assert.That(pages[0].FileName, Is.EqualTo("U1.fa.2.7.7.80.10.50.500.1.html"));
        Assert.That(pages[0].Html, Is.EqualTo(trf));
    }

    // TRF EO_MAX_TBL = 120 rows per page, heading row every 22 rows, cross-links between pages; no repeats → one page
    // with a heading row and "No Repeats Found!".
    [Test]
    public void HtmlTable_PagingAndEmpty()
    {
        var p = TandemRepeatsFinderParameters.Recommended;
        string seq = new('A', 2000);
        var fake = Enumerable.Range(0, 121).Select(k => new ApproximateTandemRepeatResult(k * 10, 8, 2, 2, "AC", 4.0, 100, 0, 16)
        { CopyMatches = 6, OutputIndex = k + 1, PercentA = 50, PercentC = 50, EntropyTrf = 1 }).ToList();
        var pages = RepeatFinder.FormatTrfHtmlTables(seq, fake, "big", p, "big.fa");
        Assert.Multiple(() =>
        {
            Assert.That(pages.Select(x => x.FileName), Is.EqualTo(new[] { "big.fa.2.7.7.80.10.50.500.1.html", "big.fa.2.7.7.80.10.50.500.2.html" }));
            Assert.That(pages[0].Html, Does.Contain("This is table  1  of  2  ( 121 repeats found )"));
            Assert.That(pages[0].Html, Does.Contain("Tables:   1   <A HREF=\"big.fa.2.7.7.80.10.50.500.2.html\" target=\"_self\">2</A>   "));
            Assert.That(CountOf(pages[0].Html, "<TD WIDTH=140>"), Is.EqualTo(6), "heading every 22 rows: rows 0,22,..,110");
            Assert.That(CountOf(pages[0].Html, "<TR><TD><CENTER><A HREF"), Is.EqualTo(120));
            Assert.That(CountOf(pages[1].Html, "<TR><TD><CENTER><A HREF"), Is.EqualTo(1));
            Assert.That(pages[0].Html, Does.Not.Contain("The End!"));
            Assert.That(pages[1].Html, Does.Contain("<P>The End!"));
            Assert.That(pages[1].Html, Does.Contain("big.fa.2.7.7.80.10.50.500.2.txt.html#1201--1208,2,4.0,2,121"));
        });

        var empty = RepeatFinder.FormatTrfHtmlTables("GATTACA", [], "none", p, "none.fa");
        Assert.That(empty, Has.Count.EqualTo(1));
        Assert.That(empty[0].Html, Does.Contain("( 0 repeats found )"));
        Assert.That(empty[0].Html, Does.Contain("<TD WIDTH=140>"));
        Assert.That(empty[0].Html, Does.Contain("\n</TABLE>\n\nNo Repeats Found!<BR>"));
    }

    private static int CountOf(string text, string part)
    {
        int count = 0;
        for (int at = text.IndexOf(part, StringComparison.Ordinal); at >= 0; at = text.IndexOf(part, at + 1, StringComparison.Ordinal))
            count++;
        return count;
    }

    // `trf two.fa 2 7 7 80 10 50 500` (seq1 test + a 300-bp sequence without repeats) → two.fa.2.7.7.80.10.50.500.summary.html.
    [Test]
    public void HtmlSummary_MatchesTrfByteForByte()
    {
        const string trf = "<HTML><HEAD><TITLE>Output Summary</TITLE></HEAD><BODY bgcolor=\"#FBF8BC\"><PRE>\n"
        + "Tandem Repeats Finder Program written by:<CENTER>\n"
        + "Gary Benson\n"
        + "Program in Bioinformatics\n"
        + "Boston University\n"
        + "Version 4.10.0</CENTER>\n"
        + "\n"
        + "Please cite:\n"
        + "G. Benson,\n"
        + "\"Tandem repeats finder: a program to analyze DNA sequences\"\n"
        + "Nucleic Acid Research(1999)\n"
        + "Vol. 27, No. 2, pp. 573-580.\n"
        + "\n"
        + "\n"
        + "<B>Multiple Sequence Summary</B>\n"
        + "\n"
        + "Only sequences containing repeats are shown!\n"
        + "\n"
        + "Click on sequence description to view repeat table.\n"
        + "\n"
        + "<TABLE BORDER=1 CELLSPACING=0 CELLPADDING=0>\n"
        + "<TR><TD WIDTH=80><CENTER>Sequence\n"
        + "Index</CENTER></TD><TD WIDTH=400><CENTER>Sequence\n"
        + "Description</CENTER></TD><TD WIDTH=80><CENTER>Number of\n"
        + "Repeats</CENTER></TD></TR>\n"
        + "<TR><TD><CENTER>1</CENTER></TD><TD><CENTER><A TARGET=\"two.fa.s1.2.7.7.80.10.50.500.1.html\" HREF=\"two.fa.s1.2.7.7.80.10.50.500.1.html\">seq1 test</A></CENTER></TD><TD><CENTER>2</CENTER></TD></TR>\n"
        + "</TABLE>\n"
        + "\n"
        + "</BODY></HTML>\n";
        var page = RepeatFinder.FormatTrfHtmlSummary([("seq1 test", 2), ("none", 0)], TandemRepeatsFinderParameters.Recommended, "two.fa");
        Assert.That(page.FileName, Is.EqualTo("two.fa.2.7.7.80.10.50.500.summary.html"));
        Assert.That(page.Html, Is.EqualTo(trf));
        Assert.That(RepeatFinder.FormatTrfHtmlSummary([("a", 0)], TandemRepeatsFinderParameters.Recommended, "x.fa").Html,
            Does.EndWith("\n</TABLE>\n\nNo Repeats Found!<BR>\n</BODY></HTML>\n"));
    }

    // C printf("%.Nf") rounds the exact binary value, ties to even; .NET "F" rounds ties away from zero.
    [TestCase(2.25, 1, "2.2")]
    [TestCase(2.75, 1, "2.8")]
    [TestCase(0.125, 2, "0.12")]
    [TestCase(0.375, 2, "0.38")]
    [TestCase(1.0, 2, "1.00")]
    [TestCase(15.65, 1, "15.7")] // the double nearest 15.65 is above it
    [TestCase(0.0, 1, "0.0")]
    [TestCase(1234.5678, 2, "1234.57")]
    public void FormatCFixed_RoundsLikeC(double value, int decimals, string expected)
    {
        Assert.That(RepeatFinder.FormatCFixed(value, decimals), Is.EqualTo(expected));
    }

    [Test]
    public void Formatters_RejectInvalidInput()
    {
        var p = TandemRepeatsFinderParameters.Recommended;
        var outside = new[] { new ApproximateTandemRepeatResult(5, 10, 2, 2, "AC", 5, 100, 0, 20) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FormatTrfDatLines("ACACACAC", outside, "x", p));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FormatTrfDatLines(null!, [], "x", p));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FormatTrfDatLines("AC", [], null!, p));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FormatTrfDatLines("AC", [], "x", null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FormatTrfDatLines("AC", [], "x", p, (TrfDatLayout)7));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FormatTrfHtmlTables("AC", [], "x", p, null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FormatTrfHtmlSummary([("a", -1)], p, "x"));
        });
    }

    #endregion
}
