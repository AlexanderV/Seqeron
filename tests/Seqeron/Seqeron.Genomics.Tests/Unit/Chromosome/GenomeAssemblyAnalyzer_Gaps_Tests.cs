// ASSEMBLY-GAPS-001 — Assembly gaps, scaffold structure and contig extraction
// Evidence: docs/Evidence/ASSEMBLY-GAPS-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-GAPS-001.md
// Sources: QUAST quast_libs/qutils.py split_by_ns / broke_scaffolds (Ns_break_threshold = 10);
//          seqtk.c stk_gap (maximal N-runs, BED [i-l, i)); NCBI "Gapped Format for Genome Submissions"
//          (every run of >= 10 Ns is an assembly gap); AGP v2.1 (short N runs inside contigs, 1-based
//          inclusive coordinates); numpy.median. Oracle: tools/review/oracles/B10/gaps_oracle.py.

using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

[TestFixture]
public class GenomeAssemblyAnalyzer_Gaps_Tests
{
    private static string Rep(char c, int n) => new(c, n);

    #region FindGaps

    // Expected = seqtk stk_gap transcription (gaps_oracle.py), End made inclusive. Default nonAcgtAsGap = true
    // (QUAST correct_seq; seqtk gap: every non-ACGT symbol is a gap symbol), so R/Y/K/M join the N runs.
    [TestCase("NNAANNNN", 1, new[] { 0, 1, 2, 4, 7, 4 })]
    [TestCase("AAANNNNNGGG", 1, new[] { 3, 7, 5 })]
    [TestCase("AAANNNNNGGG", 10, new int[0])]
    [TestCase("nnACnNnGT", 1, new[] { 0, 1, 2, 4, 6, 3 })]
    [TestCase("ACRYNNKMN", 1, new[] { 2, 8, 7 })]
    [TestCase("ACRYNNKMN", 8, new int[0])]
    [TestCase("ACRYACKMN", 1, new[] { 2, 3, 2, 6, 8, 3 })]
    [TestCase("ACGU-*T", 1, new[] { 3, 5, 3 })]
    [TestCase("NNNNNNNNNN", 1, new[] { 0, 9, 10 })]
    [TestCase("", 1, new int[0])]
    [TestCase("ACGT", 1, new int[0])]
    [TestCase("ANNNNNNNNNNCNNNNNNNNNG", 10, new[] { 1, 10, 10 })]
    public void FindGaps_MatchesSeqtkGap(string sequence, int minGapLength, int[] expectedTriples)
    {
        var gaps = FindGaps(new[] { ("s", sequence) }, minGapLength)
            .SelectMany(g => new[] { g.Start, g.End, g.Length }).ToArray();
        Assert.That(gaps, Is.EqualTo(expectedTriples));
    }

    [Test]
    public void FindGaps_NonAcgtAsGapFalse_IupacCodesAreBases_LowercaseNIsAGap()
    {
        // AGP v2.1 / NCBI: a gap is a run of N; R/Y/K/M are ambiguous bases (QUAST/seqtk call them N: the default).
        var gaps = FindGaps(new[] { ("s", "ARYnNKMA") }, nonAcgtAsGap: false).ToList();
        Assert.That(gaps, Has.Count.EqualTo(1));
        Assert.That((gaps[0].Start, gaps[0].End, gaps[0].Length), Is.EqualTo((3, 4, 2)));
    }

    [Test]
    public void FindGaps_Default_IupacCodesAndOtherNonAcgtSymbolsAreGapSymbols()
    {
        // QUAST correct_seq: M K R Y W S V B H D -> N (upper-cased first); seqtk gap: every non-ACGT is a gap.
        var gaps = FindGaps(new[] { ("s", "ARYnNKMA") }).ToList();
        var iupacOnly = FindGaps(new[] { ("s", "AmkryWsvbhdT") }).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(gaps.Select(g => (g.Start, g.End, g.Length)), Is.EqualTo(new[] { (1, 6, 6) }));
            Assert.That(iupacOnly.Select(g => (g.Start, g.End, g.Length)), Is.EqualTo(new[] { (1, 10, 10) }));
        });
    }

    [Test]
    public void ScaffoldAndContigApis_NonAcgtAsGap_AppliesConsistently()
    {
        // 6 contig bases, a 10-symbol IUPAC/U run, 6 contig bases (QUAST split_by_ns on the corrected sequence).
        var scaffold = new[] { ("sc", "ACGTAC" + "RYKMSWBDHV" + "GTACGT") };
        var asGap = AnalyzeScaffolds(scaffold).Single();
        var asBase = AnalyzeScaffolds(scaffold, nonAcgtAsGap: false).Single();
        Assert.Multiple(() =>
        {
            Assert.That(asGap.Contigs.Select(c => (c.Start, c.End)), Is.EqualTo(new[] { (0, 5), (16, 21) }));
            Assert.That(asGap.GapLength, Is.EqualTo(10));
            Assert.That(ExtractContigs(scaffold, minContigLength: 1).Select(c => c.Sequence), Is.EqualTo(new[] { "ACGTAC", "GTACGT" }));
            Assert.That(asBase.Contigs, Has.Count.EqualTo(1), "IUPAC codes are bases when nonAcgtAsGap is false");
            Assert.That(ExtractContigs(scaffold, minContigLength: 1, nonAcgtAsGap: false).Single().Sequence, Is.EqualTo(scaffold[0].Item2));
        });
    }

    [Test]
    public void FindGaps_DefaultThreshold_CountsTheGapsOfCalculateStatistics()
    {
        var seqs = new[] { ("a", "NNACGTNAANNNNNNNNNNNNC"), ("b", "ACGT"), ("c", "nnnn") };
        var gaps = FindGaps(seqs).ToList();
        var stats = CalculateStatistics(seqs);
        Assert.That(gaps, Has.Count.EqualTo(stats.TotalGaps));
        Assert.That(gaps.Sum(g => (long)g.Length), Is.EqualTo(stats.TotalGapLength));
    }

    [Test]
    public void FindGaps_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => FindGaps(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindGaps(new[] { ("s", "ANA") }, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => FindGaps(new[] { ("s", "ANA") }, -3));
        Assert.Throws<ArgumentException>(() => FindGaps(new[] { ("s", (string)null!) }).ToList());
    }

    #endregion

    #region AnalyzeGapDistribution

    // numpy.median / numpy.mean of the gap lengths.
    [TestCase(new[] { 2, 4 }, 3.0, 3.0)]
    [TestCase(new[] { 1, 5, 100 }, 5.0, 35.333333333333336)]
    [TestCase(new[] { 10, 100, 1000, 5 }, 55.0, 278.75)]
    [TestCase(new[] { 7 }, 7.0, 7.0)]
    public void AnalyzeGapDistribution_MedianAndMean_MatchNumpy(int[] lengths, double median, double mean)
    {
        var gaps = lengths.Select((l, i) => new GapInfo("s", i * 2000, i * 2000 + l - 1, l, "x")).ToList();
        var d = AnalyzeGapDistribution(gaps);
        Assert.That(d.MedianLength, Is.EqualTo(median).Within(1e-12));
        Assert.That(d.MeanLength, Is.EqualTo(mean).Within(1e-9));
        Assert.That(d.MaxLength, Is.EqualTo(lengths.Max()));
        Assert.That(d.Count, Is.EqualTo(lengths.Length));
    }

    [Test]
    public void AnalyzeGapDistribution_EmptyAndNull()
    {
        var empty = AnalyzeGapDistribution(Array.Empty<GapInfo>());
        Assert.That((empty.Count, empty.MeanLength, empty.MedianLength, empty.MaxLength), Is.EqualTo((0, 0.0, 0.0, 0)));
        Assert.That(empty.TypeCounts, Is.Empty);
        Assert.Throws<ArgumentNullException>(() => AnalyzeGapDistribution(null!));
    }

    #endregion

    #region AnalyzeScaffolds

    // Expected contigs (0-based, end exclusive) = QUAST split_by_ns pieces (gaps_oracle.py).
    [TestCase("AAAAANNNNNTTTTT", 10, new[] { 0, 15 })]
    [TestCase("AAANNNNNNNNNNAAA", 10, new[] { 0, 3, 13, 16 })]
    [TestCase("NNNNNNNNNNAAA", 10, new[] { 10, 13 })]
    [TestCase("AAAANNNNNNNNNN", 10, new[] { 0, 4 })]
    [TestCase("ACnnGTNNNNNNNNNNNNacNNNNNNNNNt", 10, new[] { 0, 6, 18, 30 })]
    [TestCase("NNN", 10, new[] { 0, 3 })]
    [TestCase("NNNNNNNNNN", 10, new int[0])]
    [TestCase("AAANNAAAA", 1, new[] { 0, 3, 5, 9 })]
    [TestCase("AAANNAAAA", 3, new[] { 0, 9 })]
    public void AnalyzeScaffolds_ContigsMatchQuastSplitByNs(string sequence, int minGapLength, int[] expectedPairs)
    {
        var s = AnalyzeScaffolds(new[] { ("sc", sequence) }, minGapLength).Single();
        Assert.That(s.Contigs.SelectMany(c => new[] { c.Start, c.End + 1 }).ToArray(), Is.EqualTo(expectedPairs));
        Assert.That(s.Contigs.Select(c => c.ContigId),
            Is.EqualTo(Enumerable.Range(1, s.Contigs.Count).Select(i => $"sc_contig{i}")));
        Assert.That(s.ContigLength + s.GapLength, Is.EqualTo(s.TotalLength));
    }

    [Test]
    public void AnalyzeScaffolds_NRunBelowMinGap_StaysInsideTheContig()
    {
        // Old code: two contigs [0,4],[10,14], no gap, ContigLength 10 of TotalLength 15.
        // QUAST split_by_ns(threshold 10) / NCBI (< 10 Ns = ambiguous bases): one contig [0,14].
        var s = AnalyzeScaffolds(new[] { ("s", "AAAAANNNNNTTTTT") }).Single();
        Assert.That(s.Contigs, Is.EqualTo(new[] { ("s_contig1", 0, 14) }));
        Assert.That(s.Gaps, Is.Empty);
        Assert.That((s.TotalLength, s.ContigLength, s.GapLength), Is.EqualTo((15, 15, 0)));
    }

    [Test]
    public void AnalyzeScaffolds_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => AnalyzeScaffolds(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalyzeScaffolds(new[] { ("s", "ANA") }, 0));
        Assert.Throws<ArgumentException>(() => AnalyzeScaffolds(new[] { ("s", (string)null!) }).ToList());
    }

    #endregion

    #region ExtractContigs

    [Test]
    public void ExtractContigs_DefaultBreaksOnlyAtTenOrMoreNs()
    {
        // QUAST broken scaffolds (Ns_break_threshold = 10) / NCBI assembly-gap convention.
        string nine = Rep('A', 300) + Rep('N', 9) + Rep('C', 300);
        string ten = Rep('A', 300) + Rep('N', 10) + Rep('C', 300);

        var fromNine = ExtractContigs(new[] { ("s", nine) }).ToList();
        Assert.That(fromNine, Is.EqualTo(new[] { ("s_contig1", nine) }));

        var fromTen = ExtractContigs(new[] { ("s", ten) }).ToList();
        Assert.That(fromTen, Is.EqualTo(new[] { ("s_contig1", Rep('A', 300)), ("s_contig2", Rep('C', 300)) }));

        // minGapLength 1 restores a cut at every N.
        Assert.That(ExtractContigs(new[] { ("s", nine) }, 200, minGapLength: 1).Select(c => c.Sequence.Length),
            Is.EqualTo(new[] { 300, 300 }));
    }

    [Test]
    public void ExtractContigs_EqualsAnalyzeScaffoldsContigsFilteredByLength()
    {
        var scaffolds = new[]
        {
            ("x", "ACnnGTNNNNNNNNNNNNacNNNNNNNNNt"),
            ("y", Rep('N', 12) + "ACGTAC" + Rep('n', 3) + "GG" + Rep('N', 15)),
        };
        foreach (int minGap in new[] { 1, 3, 10 })
        {
            var expected = AnalyzeScaffolds(scaffolds, minGap)
                .SelectMany(st => st.Contigs.Select(c => scaffolds.First(s => s.Item1 == st.ScaffoldId).Item2[c.Start..(c.End + 1)]))
                .Where(c => c.Length >= 4)
                .ToList();
            var extracted = ExtractContigs(scaffolds, 4, minGap).Select(c => c.Sequence).ToList();
            Assert.That(extracted, Is.EqualTo(expected), $"minGapLength {minGap}");
        }
    }

    // A10-2 — QUAST split_by_ns numbers every piece of the scaffold and --min-contig filters afterwards, so
    // a dropped short contig leaves a hole in the numbering; AnalyzeScaffolds names the same contigs identically.
    [Test]
    public void ExtractContigs_NumbersFromTheFullContigList_AgreeingWithAnalyzeScaffolds()
    {
        // 5-bp, 300-bp, 2-bp, 250-bp pieces separated by 10-N runs.
        string scaffold = Rep('A', 5) + Rep('N', 10) + Rep('C', 300) + Rep('N', 10) + "GG" + Rep('N', 10) + Rep('T', 250);
        var scaffolds = new[] { ("s", scaffold) };

        var extracted = ExtractContigs(scaffolds, 200).ToList();

        Assert.That(extracted.Select(c => c.Id), Is.EqualTo(new[] { "s_contig2", "s_contig4" }),
            "contigs 1 and 3 are dropped (< 200) but keep their numbers");
        Assert.That(extracted.Select(c => c.Sequence.Length), Is.EqualTo(new[] { 300, 250 }));

        var analyzed = AnalyzeScaffolds(scaffolds).Single().Contigs;
        Assert.That(analyzed.Select(c => c.ContigId), Is.EqualTo(new[] { "s_contig1", "s_contig2", "s_contig3", "s_contig4" }));
        Assert.That(extracted.Select(c => c.Id), Is.SubsetOf(analyzed.Select(c => c.ContigId)),
            "ExtractContigs and AnalyzeScaffolds name a contig identically");
    }

    [Test]
    public void ExtractContigs_LowercaseCasePreserved_MinContigZeroKeepsAll()
    {
        var contigs = ExtractContigs(new[] { ("s", "acgt" + Rep('n', 10) + "g") }, 0).ToList();
        Assert.That(contigs, Is.EqualTo(new[] { ("s_contig1", "acgt"), ("s_contig2", "g") }));
    }

    [Test]
    public void ExtractContigs_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => ExtractContigs(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExtractContigs(new[] { ("s", "A") }, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => ExtractContigs(new[] { ("s", "A") }, 1, 0));
        Assert.Throws<ArgumentException>(() => ExtractContigs(new[] { ("s", (string)null!) }).ToList());
    }

    #endregion

    #region FilterByLength / SortByLength

    [Test]
    public void SortByLength_IsStableOnTies()
    {
        var seqs = new[] { ("a", "AA"), ("b", "CCC"), ("c", "GG"), ("d", "T"), ("e", "AAA") };
        Assert.That(SortByLength(seqs).Select(s => s.Item1), Is.EqualTo(new[] { "b", "e", "a", "c", "d" }));
        Assert.That(SortByLength(seqs, descending: false).Select(s => s.Item1), Is.EqualTo(new[] { "d", "a", "c", "b", "e" }));
    }

    [Test]
    public void FilterByLength_MinAboveMax_IsEmpty_AndNullsThrow()
    {
        var seqs = new[] { ("a", "AAA") };
        Assert.That(FilterByLength(seqs, 4, 2), Is.Empty);
        Assert.Throws<ArgumentNullException>(() => FilterByLength(null!));
        Assert.Throws<ArgumentNullException>(() => SortByLength(null!));
        Assert.Throws<ArgumentException>(() => FilterByLength(new[] { ("s", (string)null!) }).ToList());
        Assert.Throws<ArgumentException>(() => SortByLength(new[] { ("s", (string)null!) }).ToList());
    }

    #endregion
}
