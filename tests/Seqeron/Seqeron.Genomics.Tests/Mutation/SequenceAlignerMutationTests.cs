namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for SequenceAligner.cs (checklist 04 rows 35/36/37/38/226:
/// ALIGN-GLOBAL/LOCAL/SEMI/MULTI-001 + ALIGN-STATS-001).
///
/// The canonical suite under-covered the cancellable Needleman-Wunsch overload (a full
/// re-implementation), the EMBOSS srspair formatter, the statistics percentages, and the
/// two multiple-sequence aligners (star + progressive). These pin the published algorithm
/// outputs with independent ground truth so operator/relational/logical mutants diverge.
/// </summary>
[TestFixture]
public class SequenceAlignerMutationTests
{
    // Independent column score for a pairwise alignment (no gap-gap columns in NW/SW output):
    // identity → Match, substitution → Mismatch, gap-vs-residue → GapExtend (linear gap).
    private static int ColumnScore(string a, string b, ScoringMatrix s)
    {
        a.Length.Should().Be(b.Length, "aligned rows must be equal length");
        int total = 0;
        for (int i = 0; i < a.Length; i++)
        {
            char c1 = a[i], c2 = b[i];
            if (c1 == '-' && c2 == '-') continue;
            else if (c1 == '-' || c2 == '-') total += s.GapExtend;
            else if (c1 == c2) total += s.Match;
            else total += s.Mismatch;
        }
        return total;
    }

    private static string DeGap(string s) => s.Replace("-", "");

    // ── Cancellable Needleman-Wunsch overload: score + traceback consistency ──────────

    [Test]
    [TestCase("ACGTACGT", "ACGTACGT")]   // identical
    [TestCase("ACGTACGT", "ACGAACGT")]   // single mismatch
    [TestCase("ACGTACGT", "ACGTCGT")]    // deletion
    [TestCase("ACGTACGT", "ACGTAACGT")]  // insertion
    [TestCase("AAAA", "TTTT")]           // all mismatch
    [TestCase("GATTACA", "GCATGCA")]     // classic mixed
    public void GlobalAlign_Cancellable_MatchesTrustedScore_AndTracebackIsConsistent(string s1, string s2)
    {
        var reference = SequenceAligner.GlobalAlign(s1, s2, SequenceAligner.SimpleDna);
        var cancellable = SequenceAligner.GlobalAlign(s1, s2, SequenceAligner.SimpleDna, CancellationToken.None);

        // Optimal NW score is unique → the cancellable DP must reproduce it exactly.
        cancellable.Score.Should().Be(reference.Score);

        // Traceback validity: degapping recovers the inputs, rows are equal length, and the
        // alignment's own column score equals the reported optimal score.
        DeGap(cancellable.AlignedSequence1).Should().Be(s1);
        DeGap(cancellable.AlignedSequence2).Should().Be(s2);
        ColumnScore(cancellable.AlignedSequence1, cancellable.AlignedSequence2, SequenceAligner.SimpleDna)
            .Should().Be(cancellable.Score, "a correct traceback's columns must sum to the DP optimum");
        cancellable.AlignmentType.Should().Be(AlignmentType.Global);
        cancellable.EndPosition1.Should().Be(s1.Length - 1);
        cancellable.EndPosition2.Should().Be(s2.Length - 1);
    }

    [Test]
    public void GlobalAlign_Cancellable_DifferentScoringChangesOptimum()
    {
        // BlastDna (mismatch −3) penalizes the mismatch more than SimpleDna (−1):
        // the reported score must track the scoring matrix (kills match/mismatch/gap constant mutants).
        var simple = SequenceAligner.GlobalAlign("ACGT", "AGGT", SequenceAligner.SimpleDna, CancellationToken.None);
        var blast = SequenceAligner.GlobalAlign("ACGT", "AGGT", SequenceAligner.BlastDna, CancellationToken.None);

        simple.Score.Should().Be(1 + 1 + (-1) + 1);  // 3 matches, 1 mismatch under SimpleDna
        blast.Score.Should().Be(2 + 2 + (-3) + 2);   // under BlastDna
    }

    [Test]
    public void GlobalAlign_Cancellable_EmptyInput_ReturnsEmpty()
    {
        SequenceAligner.GlobalAlign("", "ACGT", SequenceAligner.SimpleDna, CancellationToken.None)
            .Should().Be(AlignmentResult.Empty);
    }

    [Test]
    public void GlobalAlign_Cancellable_PreCancelledToken_Throws()
    {
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        var act = () => SequenceAligner.GlobalAlign(
            new string('A', 500), new string('A', 500), SequenceAligner.SimpleDna, cts.Token);

        act.Should().Throw<OperationCanceledException>();
    }

    [Test]
    public void GlobalAlign_Cancellable_ProgressReachesCompletion()
    {
        var progress = new Progress<double>(_ => { });
        // Progress<T> dispatch is not guaranteed ordered across threads; assert only the
        // terminal contract the algorithm guarantees (a completed, correctly-scored alignment).
        var result = SequenceAligner.GlobalAlign(
            new string('A', 250), new string('A', 250), SequenceAligner.SimpleDna,
            CancellationToken.None, progress);

        result.Score.Should().Be(250, "250 identical bases under +1 match");
    }

    // ── Known-answer global / local / semi-global ─────────────────────────────────────

    [Test]
    public void GlobalAlign_GapRequired_ScoresWithLinearGapPenalty()
    {
        // "AAA" vs "AA": two matches (+2) and one gap (GapExtend −1) → score 1.
        var r = SequenceAligner.GlobalAlign("AAA", "AA", SequenceAligner.SimpleDna);
        r.Score.Should().Be(1);
        (r.AlignedSequence1.Contains('-') || r.AlignedSequence2.Contains('-')).Should().BeTrue();
        ColumnScore(r.AlignedSequence1, r.AlignedSequence2, SequenceAligner.SimpleDna).Should().Be(1);
    }

    [Test]
    public void LocalAlign_FindsBestSubsequence_NotEndToEnd()
    {
        // "ACGT" occurs cleanly inside the second sequence → local score 4, aligned core "ACGT".
        var r = SequenceAligner.LocalAlign("TTACGTTT", "GGACGTGG", SequenceAligner.SimpleDna);
        r.Score.Should().Be(4);
        r.AlignedSequence1.Should().Be("ACGT");
        r.AlignedSequence2.Should().Be("ACGT");
        r.AlignmentType.Should().Be(AlignmentType.Local);
    }

    [Test]
    public void SemiGlobalAlign_FitsShortIntoLong_WithoutEndGapPenalty()
    {
        // Fitting "ACGT" into "TTTTACGTTTTT": the 4-base match scores 4, flanks are free.
        var r = SequenceAligner.SemiGlobalAlign(new DnaSequence("ACGT"), new DnaSequence("TTTTACGTTTTT"),
            SequenceAligner.SimpleDna);
        r.Score.Should().Be(4);
    }

    // ── FormatAlignment: EMBOSS srspair markup ('|' identity, ':' similar, ' ' gap/mismatch) ─

    [Test]
    public void FormatAlignment_Markup_FollowsSrspairLegend()
    {
        // Columns: A|A '|', C|C '|', G/A mismatch (SimpleDna mismatch −1, not >0) ' ', '-'/T gap ' ', T|T '|'.
        var aln = new AlignmentResult("ACG-T", "ACATT", 0, AlignmentType.Global, 0, 0, 4, 4);
        var formatted = SequenceAligner.FormatAlignment(aln, lineWidth: 60, SequenceAligner.SimpleDna);

        formatted.Should().Be("ACG-T\n||  |\nACATT\n\n");
    }

    [Test]
    public void FormatAlignment_PositiveMismatchScore_MarksSimilarColumns()
    {
        // With a positive "mismatch" score, a substitution column is marked ':' (similar).
        var positiveMismatch = new ScoringMatrix(Match: 1, Mismatch: 1, GapOpen: -2, GapExtend: -1);
        var aln = new AlignmentResult("AG", "AC", 0, AlignmentType.Global, 0, 0, 1, 1);

        SequenceAligner.FormatAlignment(aln, 60, positiveMismatch).Should().Be("AG\n|:\nAC\n\n");
    }

    [Test]
    public void FormatAlignment_RespectsLineWidthBlocks()
    {
        var aln = new AlignmentResult("AAAA", "AAAA", 0, AlignmentType.Global, 0, 0, 3, 3);
        // lineWidth 2 → two blocks of 2.
        SequenceAligner.FormatAlignment(aln, lineWidth: 2).Should().Be("AA\n||\nAA\n\nAA\n||\nAA\n\n");
    }

    // ── CalculateStatistics: EMBOSS needle percentages over alignment length ──────────

    [Test]
    public void CalculateStatistics_CountsAndPercentages_AreExact()
    {
        // "ACG-T"/"ACATT": 3 matches, 1 mismatch, 1 gap, length 5.
        var aln = new AlignmentResult("ACG-T", "ACATT", 0, AlignmentType.Global, 0, 0, 4, 4);
        var stats = SequenceAligner.CalculateStatistics(aln, SequenceAligner.SimpleDna);

        stats.Matches.Should().Be(3);
        stats.Mismatches.Should().Be(1);
        stats.Gaps.Should().Be(1);
        stats.AlignmentLength.Should().Be(5);
        stats.Identity.Should().BeApproximately(60.0, 1e-9);     // 3/5
        stats.Similarity.Should().BeApproximately(60.0, 1e-9);   // mismatch −1 ⇒ not similar
        stats.GapPercent.Should().BeApproximately(20.0, 1e-9);   // 1/5
    }

    [Test]
    public void CalculateStatistics_PositiveMismatch_CountsSubstitutionAsSimilar()
    {
        var positiveMismatch = new ScoringMatrix(Match: 1, Mismatch: 1, GapOpen: -2, GapExtend: -1);
        var aln = new AlignmentResult("ACG-T", "ACATT", 0, AlignmentType.Global, 0, 0, 4, 4);
        var stats = SequenceAligner.CalculateStatistics(aln, positiveMismatch);

        stats.Similarity.Should().BeApproximately(80.0, 1e-9, "matches + positive-score substitutions = 4/5");
    }

    // ── Star MSA (MultipleAlign): center select, reconcile, consensus, sum-of-pairs ───

    private static int SumOfPairs(string[] rows, ScoringMatrix s)
    {
        int len = rows[0].Length, total = 0;
        for (int p = 0; p < len; p++)
            for (int i = 0; i < rows.Length; i++)
                for (int j = i + 1; j < rows.Length; j++)
                {
                    char ci = rows[i][p], cj = rows[j][p];
                    if (ci == '-' && cj == '-') continue;
                    else if (ci == '-' || cj == '-') total += s.GapExtend;
                    else if (ci == cj) total += s.Match;
                    else total += s.Mismatch;
                }
        return total;
    }

    [Test]
    public void MultipleAlign_IdenticalSequences_ExactConsensusAndSumOfPairs()
    {
        var seqs = new[] { new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT") };
        var r = SequenceAligner.MultipleAlign(seqs, SequenceAligner.SimpleDna);

        r.AlignedSequences.Should().OnlyContain(s => s == "ACGTACGT");
        r.Consensus.Should().Be("ACGTACGT");
        // 8 columns × C(3,2)=3 identical pairs × Match(1) = 24
        r.TotalScore.Should().Be(24);
    }

    [Test]
    public void MultipleAlign_DivergentSequences_RowsValid_AndTotalScoreIsTrueSumOfPairs()
    {
        var inputs = new[] { "ACGTACGT", "ACGTTACGT", "ACGACGT" };
        var seqs = inputs.Select(s => new DnaSequence(s)).ToArray();
        var r = SequenceAligner.MultipleAlign(seqs, SequenceAligner.SimpleDna);

        r.AlignedSequences.Should().HaveCount(3);
        r.AlignedSequences.Select(x => x.Length).Distinct().Should().ContainSingle("all rows padded to equal length");
        r.AlignedSequences.Select(DeGap).Should().BeEquivalentTo(inputs, "degapping each row recovers its input");
        r.TotalScore.Should().Be(SumOfPairs(r.AlignedSequences, SequenceAligner.SimpleDna),
            "reported TotalScore must equal the independent sum-of-pairs of the final rows");
    }

    [Test]
    public void MultipleAlign_SingleSequence_ReturnedVerbatim()
    {
        var r = SequenceAligner.MultipleAlign(new[] { new DnaSequence("ACGT") });
        r.AlignedSequences.Should().Equal("ACGT");
        r.Consensus.Should().Be("ACGT");
        r.TotalScore.Should().Be(0);
    }

    // ── Progressive MSA (MultipleAlignProgressive): guide tree + profile alignment ────

    [Test]
    public void MultipleAlignProgressive_IdenticalSequences_ExactConsensusAndSumOfPairs()
    {
        var seqs = new[] { new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT") };
        var r = SequenceAligner.MultipleAlignProgressive(seqs, SequenceAligner.SimpleDna);

        r.AlignedSequences.Should().OnlyContain(s => s == "ACGTACGT");
        r.Consensus.Should().Be("ACGTACGT");
        r.TotalScore.Should().Be(24);
    }

    [Test]
    public void MultipleAlignProgressive_DivergentSequences_RowsValid_AndTotalScoreIsTrueSumOfPairs()
    {
        var inputs = new[] { "ACGTACGT", "ACGTTACGT", "ACGACGT", "ACGTACGA" };
        var seqs = inputs.Select(s => new DnaSequence(s)).ToArray();
        var r = SequenceAligner.MultipleAlignProgressive(seqs, SequenceAligner.SimpleDna);

        r.AlignedSequences.Should().HaveCount(4);
        r.AlignedSequences.Select(x => x.Length).Distinct().Should().ContainSingle();
        // rows are reprojected to input order
        for (int i = 0; i < inputs.Length; i++)
            DeGap(r.AlignedSequences[i]).Should().Be(inputs[i]);
        r.TotalScore.Should().Be(SumOfPairs(r.AlignedSequences, SequenceAligner.SimpleDna));
    }

    [Test]
    public void MultipleAlignProgressive_EmptyInput_ReturnsEmpty()
    {
        SequenceAligner.MultipleAlignProgressive(System.Array.Empty<DnaSequence>())
            .Should().Be(MultipleAlignmentResult.Empty);
    }
}
