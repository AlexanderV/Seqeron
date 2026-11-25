using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// ASSEMBLY-REPEAT-001 (B10 sweep U3): <see cref="FindRepetitiveRegions"/> (NCBI WindowMasker model),
/// <see cref="FindTandemRepeats"/> (Tandem Repeats Finder, Benson 1999) and <see cref="CalculateRepeatContent"/>
/// (RepeatMasker "bases masked"), plus the canonical linguistic complexity behind <see cref="CalculateLocalQuality"/>.
/// Oracles: a literal Python port of WindowMasker (3000/3000 random cases identical); an independent cyclic
/// wraparound-DP scorer for TRF (346/348 scores identical; the two others are patterns &gt; 20, TRF's narrow band);
/// a Python linguistic-complexity implementation (300/300).
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_RepeatAnalysis_Tests
{
    private const double Tol = 1e-9;

    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    #region FindRepetitiveRegions (WindowMasker)

    // seqs, k, T_threshold, window, T_extend, bothStrands, expected (id, start, end, copies) — from the Python port of
    // CSeqMasker::DoMask (seq_masker.cpp) with CSeqMaskerWindow / CSeqMaskerScoreMean and win_mask_gen_counts counting.
    private static readonly object[][] WindowMaskerOracle =
    {
        new object[] { new[] { "CARTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGTCGAGCTACTTCGAGCTACTCGAGCTACTCGAGCTACTCGAGCTACTCGAGCTACTCGAGCTACCGAGCTACTCGAGCTACTCGAGCTACGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATT" }, 4, 1, 5, 1, true, new (string, int, int, int)[] { ("s0", 3, 210, 22) } },
        new object[] { new[] { "ACTAACGCCGCCGGGTCGTTACTCGAAAAGCAGGTGGAATTGGTGTATTCATCTTTCGTCTCNTTAnGATCTTTCGTCTCATTANGAGCTTGCTCGATTTGATCGATCTGCAAGGTGCTGTCTAGATAGATACCATGGCCCGGAAGTACGGG" }, 6, 1, 18, 1, false, new (string, int, int, int)[] { ("s0", 0, 61, 2), ("s0", 85, 151, 2) } },
        new object[] { new[] { "CATGCATGAAACAGATACATCGCTTGGGCCACGTAACAAAATGTGCTCCAACAAAATGTGCTCCAACAAAATGTGCTCCAGTCTAGAGCACACTAAATGGAGCACATTTTGTTGAGACATCTTAGAGGAGATAGGCGTAGATCCGGTRTGGAGCACATTTTGTTGGAGCACATTTTGTAnTAGCCGTGATGRAAGGTGGGNGAACGGGA" }, 6, 2, 19, 2, false, new (string, int, int, int)[] { ("s0", 27, 86, 3), ("s0", 89, 119, 4), ("s0", 148, 178, 4) } },
        new object[] { new[] { "TATCAATACTGAATAGGCTTTGGGGTTAAATTGGGGTTAAATTGGGGTTAAAACYTATGTTATACTCCYTGGGGTTAAAGTGTCGTAAGGATGAnGGCTCCTTGGGGTTAAATTGGGGTTAAATTGGGGTTAAATTGGGGTTAAATTGGGGTTAAAGC", "CCTCAGCCGTTGACCGCAACACCGTGAAGCACGGGTAAGTTGGGGTTAAATTGGGGTTAAATTGGGGTTAAATTGGGGTTAAAGCAGCAGTTGGGGTTAAATTGGGGTTAAAAAAGGCGAGAACTGCAGGAGAGCGTATTTGCGCAACCCTGAGGGTCTAGAGAGTCCACCT" }, 8, 3, 20, 2, false, new (string, int, int, int)[] { ("s0", 8, 53, 15), ("s0", 69, 89, 15), ("s0", 95, 157, 15), ("s1", 28, 122, 15) } },
        new object[] { new[] { "TCTCGGCAGCCCTTAACTCCGCGGATTTGGTTGGCCCCCTTGGCCCCCCCCCCTTGGCCCCCTTGGCCCCCTTGGCCCCCTTGGCCCCCTTGGCCCCCTTGGCCCCCTATCCCAGAGCAAATGATTGCTGG" }, 7, 5, 13, 2, true, new (string, int, int, int)[] { ("s0", 28, 109, 8) } },
        new object[] { new[] { "acgtcatatgcttggtaraacgtccctattatagttattagagggacaattatagttattaragggacaattatagttattattatagttattagagggacaattatagttattagagggacaattatagttattagagggacaagagggacaattatagttattaragggacaattatygttattagagggacaaaggggggttattatagttattagagggacattggtcctggttagtgtctctccgagcttggcatgagtttatgtcgcctaagcttctcactgg" }, 15, 5, 16, 1, false, new (string, int, int, int)[] { ("s0", 27, 60, 6), ("s0", 62, 165, 6), ("s0", 180, 225, 6) } },
        new object[] { new[] { "atgactactatccaaaagagttattgtttctttaggcgaacaaggacttattataaccttgcgccccccacttgttatctgagactgctggaagttgttttaatgcaagactacctacgtgccagttgcagtccccgagctgcttaggcactcgtcggg", "AAATGCAACCCATCCNGATGGCACATTCGAGnGTGARAGCAGCAAAGCAGTTGACCGAGCGCTTTGACCACAGG" }, 3, 3, 3, 2, true, new (string, int, int, int)[] { ("s0", 0, 158, 15), ("s1", 0, 14, 15), ("s1", 16, 30, 15), ("s1", 32, 35, 6), ("s1", 37, 73, 15) } },
        new object[] { new[] { "CATATCTGTCACCTTTGGAGATTCCGATATTATAACGTGGGCTCCTACCCGCACTAGGGTCGTACTCGGATTTGATTCGAGTCGTGTACCACGGCCTGGACTGGTGGTAAAGGCTCCGATTGGAGCGGCATAGCGGCATAGCGGCATTATC" }, 4, 1, 10, 1, false, new (string, int, int, int)[] { ("s0", 0, 150, 4) } },
        new object[] { new[] { "CCTATTAAGTGTATGTACGGCACGTGGGTCCTCAACAAATACGCCTATAACTTCCTGTATGGTCTTCCTGTATGGTCTTCCTGTATGGTCTTCCTGTATGGTCTTCCTGTATGGTTGTCGCCTG", "GTCTCGACCTCATGTTCACCATACCATACAGGAAGACCATACAGGAAGACCATACAGGAAGACCATACAGGAAGACAGGAAGACCATACCTTCCTGTATGGTCTTCCTGTATGGTCTTCCTGTATGGTAGGAAGACCATACAGGAAGACCATACAGGAAGCAACTCTGTAAAGCCTGTGCTTAACGTG" }, 6, 4, 7, 4, true, new (string, int, int, int)[] { ("s0", 8, 15, 15), ("s0", 49, 115, 16), ("s1", 16, 160, 16) } },
        new object[] { new[] { "AGATCCCGTGGAGCCCTGTTCCCAAGGACGTCAAGGACGTCAAGGACGTCAAGGACGTTCGCCCGAACAGACTTAAACTTGCCTCCGTTGCCACCAGCAGTCCGCCCTCCCAGCTTGCYAAAGTAAGGGCCGCCGGGGAA", "GTAGNATTGTCGCAGATATATCTGACCNGCGGATGATATAACCATTCACCTGGACCACGCAAGGACGTCAAGGACGTCAAGGACGTCAACGTCCTTGAGGYCGTCAAGRACGTCAAGGACGTGGTGTGCATCGAGCGGGCGGGTATCTCCGTTAAGCTAGCGGTTCGCCTGA" }, 6, 1, 7, 1, true, new (string, int, int, int)[] { ("s0", 0, 117, 12), ("s0", 119, 139, 1), ("s1", 5, 26, 3), ("s1", 28, 99, 12), ("s1", 101, 107, 8), ("s1", 109, 171, 9) } },
        new object[] { new[] { "ccctcaattgtgaaaacctcaagaccgcgccattcatgaacaaccttcgctgaacaaccttcgctgaacgcgaaggttgttcagcgaaggttgttcaaaccttcgctctcgcggatgcgaaggttgttcagcgaaggttgttcagcgaaggttgttcagcgaaggttgttcagcgaaggttgttcagcgaaggttgttcac", "tttntctccatgctagaatactttgactagtggttcggcgagatggaagttcaggaggcattgggtcgattgtaanttagccataacaaaatttgt" }, 8, 4, 9, 1, true, new (string, int, int, int)[] { ("s0", 35, 200, 11) } },
        new object[] { new[] { "tgaagtacgcaaatactcgtaaaagccttcgatacagcttagacagcaaacgattcaagagactggggtgaggtaacaacttgcaccagctatgagaccatgtccgaggatcgtgacgttgaccgcttgtgtgacgttgaccgcttgtgtgacgttgaccgcttgtgtgacgttgaccgcttgttgggtcagtggcgataatgtgacgttgaccgcttgtgtgacgttgaccgcttgtgtgacgttgaccgcttgtgtgacgttgaccgcttgtgtgacgttgaccgcttgttc" }, 15, 4, 15, 4, true, new (string, int, int, int)[] { ("s0", 112, 183, 9), ("s0", 201, 291, 9) } },
        new object[] { new[] { "TTCTGGGGACGAGGCCACTTTGGATCAAACAGAGTNCTTCAAACAGAGTACTTCAAACAGAGTACTTCAAACAGAGTACTTCAAACAGAGTACTTNAAACAGAGTACTTGAAGTACTCTGTTTGGATCTCGCTGCATGGGTCACT" }, 3, 6, 4, 3, true, new (string, int, int, int)[] { ("s0", 0, 34, 15), ("s0", 36, 94, 15), ("s0", 96, 144, 15) } },
        new object[] { new[] { "ACCTGCCCGCACTATGATGTGCCAAATCTAGGAAATTGTTCGCTCTGCCTACAAATGGACGACAAATTCGACGACAAATTCCCCGTGACCAGTGACGACAAACTGGTCACGGGGAATTTGTCGTCACTGGTCACGGGGAATTTGTCGTCATTCCCCGTGACCAGTGACGACAAATTCCCCGTGACCAGTCCCGTGACCAGTGACGACAAATTCCCCGTGACCAGTGACGACAAATTCCCCGTGACCAGTCGGAACTGCACCCTTTGGACAATTGGAAGCCGGAGAGACnCCGAGTTGTGnCGCT" }, 15, 5, 53, 4, true, new (string, int, int, int)[] { ("s0", 97, 153, 7), ("s0", 188, 252, 7) } },
    };

    [TestCaseSource(nameof(WindowMaskerOracle))]
    public void FindRepetitiveRegions_MatchesWindowMaskerOracle(
        string[] seqs, int k, int t, int w, int e, bool both, (string, int, int, int)[] expected)
    {
        var named = seqs.Select((s, i) => ($"s{i}", s)).ToList();

        var got = FindRepetitiveRegions(named, k, t, w, both, e)
            .Select(r => (r.SequenceId, r.Start, r.End, r.Copies)).ToArray();

        Assert.That(got, Is.EqualTo(expected));
    }

    [Test]
    public void FindRepetitiveRegions_TandemArray_OneRegionWithMaxUnitCount()
    {
        // (CAG)×200: 6-mers CAGCAG/AGCAGC/GCAGCA occur 199/198/198 times ⇒ [0, 599], Copies 199.
        var r = FindRepetitiveRegions(new[] { ("s", Rep("CAG", 200)) }, 6, 3, 60).Single();

        Assert.That((r.Start, r.End, r.Copies), Is.EqualTo((0, 599, 199)));
    }

    [Test]
    public void FindRepetitiveRegions_UniqueSequence_NoRegion()
    {
        // de Bruijn-like unique 8-mers: every unit count is 1 < T_threshold 3.
        const string seq = "ACGTTGCAAGCTTCGATCCGGATAGCTAGGCTTAACCGGTTAGCCATGCAAT";
        Assert.That(FindRepetitiveRegions(new[] { ("s", seq) }, 8, 3, 20), Is.Empty);
    }

    [Test]
    public void FindRepetitiveRegions_CountsArePooledOverSequencesAndBothStrands()
    {
        // k = W = 5: a window is one unit and its score is the unit count. ACGTG ×2 on the forward strand and its
        // reverse complement CACGT once: forward-only counts are 2 (< 3, nothing); both strands count 3 each.
        var seqs = new[] { ("a", "ACGTG"), ("b", "ACGTG"), ("c", "CACGT") };

        Assert.That(FindRepetitiveRegions(seqs, 5, 3, 5, bothStrands: false), Is.Empty);

        var both = FindRepetitiveRegions(seqs, 5, 3, 5, bothStrands: true).ToList();
        Assert.That(both.Select(r => (r.SequenceId, r.Start, r.End, r.Copies)),
            Is.EqualTo(new[] { ("a", 0, 4, 3), ("b", 0, 4, 3), ("c", 0, 4, 3) }));
    }

    [Test]
    public void FindRepetitiveRegions_KmersWithAnyNonAcgtSymbolAreNotCounted()
    {
        // Three identical 5-mer units, each in its own sequence, but the third carries an IUPAC R: its unit is
        // skipped (Jellyfish / WindowMasker ambig rule) so the count of ACGTG stays 2 (< 3).
        var seqs = new[] { ("a", "ACGTG"), ("b", "ACGTG"), ("c", "ACRTG") };
        Assert.That(FindRepetitiveRegions(seqs, 5, 3, 5, bothStrands: false), Is.Empty);
    }

    [Test]
    public void FindRepetitiveRegions_IsCaseInsensitive()
    {
        string seq = Rep("CAG", 200);
        var upper = FindRepetitiveRegions(new[] { ("s", seq) }, 6, 3, 60).ToList();
        var lower = FindRepetitiveRegions(new[] { ("s", seq.ToLowerInvariant()) }, 6, 3, 60).ToList();
        Assert.That(lower, Is.EqualTo(upper));
    }

    [Test]
    public void FindRepetitiveRegions_InvalidArguments_ThrowEagerly()
    {
        var seqs = new[] { ("s", "ACGT") };
        Assert.Throws<ArgumentNullException>(() => FindRepetitiveRegions(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindRepetitiveRegions(seqs, kmerSize: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindRepetitiveRegions(seqs, minCopies: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindRepetitiveRegions(seqs, kmerSize: 15, windowSize: 14));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindRepetitiveRegions(seqs, minCopies: 3, extendCopies: 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindRepetitiveRegions(seqs, extendCopies: 0));
    }

    [Test]
    public void FindRepetitiveRegions_ShortOrEmptySequences_YieldNothing()
    {
        Assert.That(FindRepetitiveRegions(new[] { ("s", "ACGT"), ("e", "") }, 5, 1, 10), Is.Empty);
    }

    #endregion

    #region FindTandemRepeats (TRF)

    [Test]
    public void FindTandemRepeats_PerfectArray_TrfValues()
    {
        // (CAG)×30: TRF period 3, 30 copies, 100 % matches, score 2·90 = 180.
        var r = FindTandemRepeats(new[] { ("s", Rep("CAG", 30)) }).Single();

        Assert.That((r.Start, r.End, r.Unit), Is.EqualTo((0, 89, "CAG")));
        Assert.That(r.Copies, Is.EqualTo(30.0).Within(Tol));
        Assert.That(r.Purity, Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void FindTandemRepeats_FlankedArray_FractionalCopyNumberAndPhase()
    {
        // ACGTGCATTG · (CAG)×9 · TGCATC: the G before the array extends it by one base, so TRF reports
        // [9, 36] (28 bp, score 56), consensus GCA (phase of the first repeat base), 28/3 copies.
        var r = FindTandemRepeats(
            new[] { ("s", "ACGTGCATTG" + Rep("CAG", 9) + "TGCATC") }, 3, 5, 5).Single();

        Assert.That((r.Start, r.End, r.Unit), Is.EqualTo((9, 36, "GCA")));
        Assert.That(r.Copies, Is.EqualTo(28.0 / 3.0).Within(Tol));
        Assert.That(r.Purity, Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void FindTandemRepeats_MatchesRepeatFinderTrfRows()
    {
        // The method is the TRF finder (RepeatFinder.FindApproximateTandemRepeats, recommended parameters
        // 2 7 7 80 10 50) with period / copy / score filters, not a second implementation.
        string seq = "GATTACAGATTACA" + Rep("ACGTTGA", 14) + "TTGCAGGCATC";
        var parameters = new TandemRepeatsFinderParameters { MaxPeriod = 50, MinScore = 50 };
        var expected = RepeatFinder.FindApproximateTandemRepeats(seq, parameters, 2)
            .Where(r => r.CopyNumber >= 3)
            .Select(r => ("s", r.Start, r.Start + r.SpanLength - 1, r.Consensus, r.CopyNumber, r.PercentMatches / 100.0))
            .ToList();

        var got = FindTandemRepeats(new[] { ("s", seq) }).ToList();

        Assert.That(got, Is.Not.Empty);
        Assert.That(got, Is.EqualTo(expected));
    }

    [Test]
    public void FindTandemRepeats_BelowMinScore_NotReported()
    {
        // (AT)×4 scores 2·8 = 16 < Minscore 50; lowering Minscore to 16 reports it.
        var seq = new[] { ("s", "GCGCTAGGAC" + "ATATATAT" + "CCTAGGGTCA") };
        Assert.That(FindTandemRepeats(seq, 2, 50, 3), Is.Empty);
        Assert.That(FindTandemRepeats(seq, 2, 50, 3, minScore: 16).Any(r => r.Unit.Length == 2), Is.True);
    }

    [Test]
    public void FindTandemRepeats_MinCopiesAndPeriodRange_Filter()
    {
        var seq = new[] { ("s", Rep("ACGTT", 20)) };
        Assert.That(FindTandemRepeats(seq, 2, 50, 21), Is.Empty);               // 20 copies < 21
        Assert.That(FindTandemRepeats(seq, 6, 50, 3), Is.Empty);                // period 5 < minUnitLength 6
        Assert.That(FindTandemRepeats(seq, 2, 4, 3), Is.Empty);                 // period 5 > MaxPeriod 4
        Assert.That(FindTandemRepeats(seq, 2, 50, 3).Single().Unit.Length, Is.EqualTo(5));
    }

    [Test]
    public void FindTandemRepeats_NNeverMatches()
    {
        // N never matches in TRF (every N column costs −7): a 12-N run (−84) is bridged by the 2·45 = 90
        // gain of the second array (one repeat spanning it), a 40-N run (−280) is not: two separate arrays.
        string bridged = Rep("CAG", 15) + new string('N', 12) + Rep("CAG", 15);
        var one = FindTandemRepeats(new[] { ("s", bridged) }).ToList();
        Assert.That(one, Has.Count.EqualTo(1));
        Assert.That((one[0].Start, one[0].End), Is.EqualTo((0, bridged.Length - 1)));

        string split = Rep("CAG", 15) + new string('N', 40) + Rep("CAG", 15);
        var two = FindTandemRepeats(new[] { ("s", split) }).ToList();
        Assert.That(two.Select(x => (x.Start, x.End)), Is.EqualTo(new[] { (0, 44), (85, 129) }));
    }

    [Test]
    public void FindTandemRepeats_NoRepeat_Empty()
    {
        Assert.That(FindTandemRepeats(new[] { ("s", "ACGTGCATACGTGCATACGT") }), Is.Empty);
        Assert.That(FindTandemRepeats(new[] { ("s", "") }), Is.Empty);
    }

    [Test]
    public void FindTandemRepeats_InvalidArguments_ThrowEagerly()
    {
        var seqs = new[] { ("s", "ACGT") };
        Assert.Throws<ArgumentNullException>(() => FindTandemRepeats(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindTandemRepeats(seqs, minUnitLength: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindTandemRepeats(seqs, minUnitLength: 5, maxUnitLength: 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindTandemRepeats(seqs, maxUnitLength: 2001));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindTandemRepeats(seqs, minCopies: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindTandemRepeats(seqs, minScore: 0));
    }

    #endregion

    #region CalculateRepeatContent (RepeatMasker bases masked)

    private static RepeatAnnotation Ann(string seq, int s, int e, string cls) => new(seq, s, e, cls, cls, 0, '+');

    [Test]
    public void CalculateRepeatContent_OverlappingAnnotations_AreCountedOnce()
    {
        // LINE 0–99, SINE 50–149 (overlap 50), SINE 120–130 (contained), same sequence: union 0–149 = 150 bases.
        // The overlap goes to the first annotation (sweep in start order): LINE 100, SINE 50.
        var (total, pct, classes) = CalculateRepeatContent(
            new[] { Ann("s", 0, 99, "LINE"), Ann("s", 50, 149, "SINE"), Ann("s", 120, 130, "SINE") }, 1000);

        Assert.That(total, Is.EqualTo(150));
        Assert.That(pct, Is.EqualTo(15.0).Within(Tol));
        Assert.That(classes["LINE"], Is.EqualTo(100));
        Assert.That(classes["SINE"], Is.EqualTo(50));
        Assert.That(classes.Values.Sum(), Is.EqualTo(total));
    }

    [Test]
    public void CalculateRepeatContent_SameCoordinatesOnDifferentSequences_AreSeparate()
    {
        var (total, _, _) = CalculateRepeatContent(new[] { Ann("a", 0, 99, "LINE"), Ann("b", 0, 99, "LINE") }, 1000);
        Assert.That(total, Is.EqualTo(200));
    }

    [Test]
    public void CalculateRepeatContent_AdjacentAnnotations_AddUp_AndInputOrderIsIrrelevantForTotal()
    {
        var (total, _, classes) = CalculateRepeatContent(new[] { Ann("s", 100, 149, "SINE"), Ann("s", 0, 99, "LINE") }, 1000);
        Assert.That(total, Is.EqualTo(150));
        Assert.That(classes["LINE"], Is.EqualTo(100));
        Assert.That(classes["SINE"], Is.EqualTo(50));
    }

    [Test]
    public void CalculateRepeatContent_RandomAnnotations_MatchPerBaseOracle()
    {
        var rng = new Random(20261125);
        for (int trial = 0; trial < 2000; trial++)
        {
            int nSeq = rng.Next(1, 4);
            var anns = new List<RepeatAnnotation>();
            int count = rng.Next(0, 12);
            for (int i = 0; i < count; i++)
            {
                int s = rng.Next(0, 80), e = s + rng.Next(0, 30);
                anns.Add(Ann("q" + rng.Next(nSeq), s, e, "C" + rng.Next(3)));
            }

            // Naive oracle: mark every base once, attribute it to the first annotation in (start, input order).
            var owner = new Dictionary<(string, int), string>();
            foreach (var (a, _) in anns.Select((a, i) => (a, i)).OrderBy(x => x.a.Start).ThenBy(x => x.i))
                for (int p = a.Start; p <= a.End; p++)
                    owner.TryAdd((a.SequenceId, p), a.RepeatClass);

            var (total, pct, classes) = CalculateRepeatContent(anns, 100_000);

            Assert.That(total, Is.EqualTo(owner.Count));
            Assert.That(pct, Is.EqualTo(owner.Count * 100.0 / 100_000).Within(Tol));
            foreach (var g in owner.GroupBy(kv => kv.Value))
                Assert.That(classes[g.Key], Is.EqualTo(g.Count()));
            Assert.That(classes.Values.Sum(), Is.EqualTo(total));
        }
    }

    [Test]
    public void CalculateRepeatContent_InvalidInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => CalculateRepeatContent(null!, 100));
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRepeatContent(Array.Empty<RepeatAnnotation>(), -1));
        Assert.Throws<ArgumentException>(() => CalculateRepeatContent(new[] { Ann("s", 10, 5, "X") }, 100));
        Assert.Throws<ArgumentException>(() => CalculateRepeatContent(new[] { Ann("s", -1, 5, "X") }, 100));
        // more masked bases than genome bases is impossible
        Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRepeatContent(new[] { Ann("s", 0, 199, "X") }, 100));
    }

    #endregion

    #region Linguistic complexity (canonical) in CalculateLocalQuality

    [Test]
    public void CalculateLocalQuality_Complexity_IsCanonicalLinguisticComplexity()
    {
        // Window values: Python Troyanskaya LC oracle (m = 10, alphabet 4 extended by N) 300/300 identical.
        string[] windows = { "ACGTACGT", Rep("A", 50), "ACGTNNACGTTGCA", "ACGATCGATTAGCCGATAGGCTA" };
        foreach (var w in windows)
        {
            var q = CalculateLocalQuality(new[] { ("s", w) }).Single();
            Assert.That(q.Complexity, Is.EqualTo(SequenceComplexity.CalculateLinguisticComplexity(w)).Within(Tol), w);
        }
    }

    [Test]
    public void CalculateLocalQuality_Complexity_CountsWordsContainingN()
    {
        // Canonical LC keeps N as a fifth symbol (Orlov & Potapov 2004): "ACNGT" has V = 5,4,3,2,1 and
        // a = 5 ⇒ Σmin(5^i, 6−i) = 5+4+3+2+1: LC = 15/15 = 1 (the old private routine skipped N-words, ≠ 1 here).
        Assert.That(CalculateLocalQuality(new[] { ("s", "ACNGT") }).Single().Complexity, Is.EqualTo(1.0).Within(Tol));
    }

    #endregion
}
