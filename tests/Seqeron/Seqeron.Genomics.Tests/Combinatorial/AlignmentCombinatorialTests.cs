namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Alignment area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Alignment")]
public class AlignmentCombinatorialTests
{
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    /// <summary>Independent linear-gap Needleman-Wunsch optimal score (ground truth).</summary>
    private static int NwOptimalScore(string a, string b, int match, int mismatch, int gap)
    {
        int m = a.Length, n = b.Length;
        var dp = new int[m + 1, n + 1];
        for (int i = 0; i <= m; i++) dp[i, 0] = i * gap;
        for (int j = 0; j <= n; j++) dp[0, j] = j * gap;
        for (int i = 1; i <= m; i++)
            for (int j = 1; j <= n; j++)
            {
                int s = a[i - 1] == b[j - 1] ? match : mismatch;
                dp[i, j] = Math.Max(dp[i - 1, j - 1] + s, Math.Max(dp[i - 1, j] + gap, dp[i, j - 1] + gap));
            }
        return dp[m, n];
    }

    /// <summary>Independent linear-gap Smith-Waterman optimal local score (ground truth).</summary>
    private static int SwOptimalScore(string a, string b, int match, int mismatch, int gap)
    {
        int m = a.Length, n = b.Length, best = 0;
        var dp = new int[m + 1, n + 1];
        for (int i = 1; i <= m; i++)
            for (int j = 1; j <= n; j++)
            {
                int s = a[i - 1] == b[j - 1] ? match : mismatch;
                dp[i, j] = Math.Max(0, Math.Max(dp[i - 1, j - 1] + s, Math.Max(dp[i - 1, j] + gap, dp[i, j - 1] + gap)));
                best = Math.Max(best, dp[i, j]);
            }
        return best;
    }

    /// <summary>Independent semi-global (fitting) optimum: query q fully aligned, free leading/trailing reference gaps.</summary>
    private static int SemiGlobalOptimalScore(string q, string r, int match, int mismatch, int gap)
    {
        int m = q.Length, n = r.Length;
        var dp = new int[m + 1, n + 1];
        for (int j = 0; j <= n; j++) dp[0, j] = 0;     // free leading reference gaps
        for (int i = 1; i <= m; i++) dp[i, 0] = i * gap; // query leading deletions are penalised
        for (int i = 1; i <= m; i++)
            for (int j = 1; j <= n; j++)
            {
                int s = q[i - 1] == r[j - 1] ? match : mismatch;
                dp[i, j] = Math.Max(dp[i - 1, j - 1] + s, Math.Max(dp[i - 1, j] + gap, dp[i, j - 1] + gap));
            }
        int best = int.MinValue;
        for (int j = 0; j <= n; j++) best = Math.Max(best, dp[m, j]); // free trailing reference gaps
        return best;
    }

    /// <summary>Score of a concrete column-wise alignment under the linear-gap model.</summary>
    private static int ScoreColumns(string a1, string a2, int match, int mismatch, int gap)
    {
        int s = 0;
        for (int i = 0; i < a1.Length; i++)
        {
            if (a1[i] == '-' || a2[i] == '-') s += gap;
            else s += a1[i] == a2[i] ? match : mismatch;
        }
        return s;
    }

    /// <summary>
    /// A set of related-but-distinct DNA sequences: a shared base, each variant carrying one
    /// substitution and (for index &gt; 0) one deletion — so anchors exist, indels appear in the
    /// MSA, and no two inputs are identical.
    /// </summary>
    private static string[] MutatedDnaSet(int n, int baseLen, uint seed)
    {
        string baseSeq = DiverseDna(baseLen, seed);
        var set = new string[n];
        set[0] = baseSeq;
        for (int s = 1; s < n; s++)
        {
            var chars = baseSeq.ToCharArray().ToList();
            int subPos = (s * 3) % baseLen;
            chars[subPos] = NextBase(chars[subPos], s); // always a different base
            int delPos = (s * 5 + 1) % chars.Count;
            chars.RemoveAt(delPos);                      // one deletion -> forces an indel
            set[s] = new string(chars.ToArray());
        }
        return set;
    }

    private static char NextBase(char c, int salt)
    {
        const string b = "ACGT";
        int idx = b.IndexOf(c);
        return b[(idx + 1 + (salt & 1)) % 4]; // rotates by 1 or 2 -> never the same base
    }

    /// <summary>
    /// Independent sum-of-pairs score over an MSA per the Wikipedia definition: for every column
    /// and every C(k,2) row pair, gap-gap is neutral, gap-nucleotide costs <paramref name="gap"/>,
    /// equal residues score <paramref name="match"/>, otherwise <paramref name="mismatch"/>.
    /// </summary>
    private static int SumOfPairsScore(IReadOnlyList<string> rows, int match, int mismatch, int gap)
    {
        if (rows.Count < 2) return 0;
        int len = rows[0].Length;
        int total = 0;
        for (int c = 0; c < len; c++)
            for (int i = 0; i < rows.Count; i++)
                for (int j = i + 1; j < rows.Count; j++)
                {
                    char a = rows[i][c], b = rows[j][c];
                    if (a == '-' && b == '-') continue;
                    if (a == '-' || b == '-') total += gap;
                    else if (a == b) total += match;
                    else total += mismatch;
                }
        return total;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ALIGN-GLOBAL-001 — Global (Needleman-Wunsch) alignment (Alignment)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 35.
    // Spec: tests/TestSpecs/ALIGN-GLOBAL-001.md (canonical SequenceAligner.GlobalAlign).
    // Dimensions: matchScore(3) × mismatchPen(3) × gapPen(3) × seqLen(3). Grid 3⁴ = 81.
    //
    // Model (Needleman & Wunsch 1970): the optimal global alignment maximises the total
    // column score over end-to-end alignments, F(i,j)=max(F(i−1,j−1)+s, F(i−1,j)+d,
    // F(i,j−1)+d) with linear gap d (this implementation uses GapExtend as d). The returned
    // alignment must (a) reconstruct both inputs when gaps are removed, (b) be equal length,
    // (c) have a column score equal to the reported Score, and (d) that Score must equal the
    // independent NW optimum for the given (match, mismatch, gap).
    //
    // The combinatorial point: the three scoring weights and the length jointly determine the
    // optimum and the traceback; correctness is asserted against an independent DP for every
    // weight combination, so no single weight can silently corrupt the result.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void AlignGlobal_OptimalAndConsistent_AcrossScoringAndLength(
        [Values(1, 2, 5)] int matchScore,
        [Values(-1, -3, 0)] int mismatchPen,
        [Values(-1, -2, -5)] int gapPen,
        [Values(4, 8, 12)] int seqLen)
    {
        string s1 = DiverseDna(seqLen, 0x1234u);
        string s2 = DiverseDna(seqLen, 0x9ABCu);
        var scoring = new ScoringMatrix(matchScore, mismatchPen, GapOpen: gapPen, GapExtend: gapPen);

        var r = SequenceAligner.GlobalAlign(s1, s2, scoring);

        // Valid end-to-end alignment of both inputs.
        r.AlignedSequence1.Replace("-", "").Should().Be(s1, "ungapped row 1 reconstructs sequence 1");
        r.AlignedSequence2.Replace("-", "").Should().Be(s2, "ungapped row 2 reconstructs sequence 2");
        r.AlignedSequence1.Length.Should().Be(r.AlignedSequence2.Length, "alignment rows share a column count");

        // Score is self-consistent with the traceback and equals the NW optimum.
        ScoreColumns(r.AlignedSequence1, r.AlignedSequence2, matchScore, mismatchPen, gapPen)
            .Should().Be(r.Score, "the reported score is the column score of the returned alignment");
        r.Score.Should().Be(NwOptimalScore(s1, s2, matchScore, mismatchPen, gapPen),
            "the score is the Needleman-Wunsch optimum for these weights");
        r.AlignmentType.Should().Be(AlignmentType.Global);
    }

    /// <summary>
    /// Worked example: identical sequences align gaplessly along the diagonal, scoring
    /// length × matchScore — independent of the gap and mismatch weights.
    /// </summary>
    [Test, Combinatorial]
    public void AlignGlobal_IdenticalSequences_ScoreIsLengthTimesMatch(
        [Values(1, 2, 5)] int matchScore,
        [Values(-1, -5)] int gapPen)
    {
        string s = DiverseDna(10, 0x55AAu);
        var r = SequenceAligner.GlobalAlign(s, s, new ScoringMatrix(matchScore, -3, gapPen, gapPen));

        r.Score.Should().Be(s.Length * matchScore);
        r.AlignedSequence1.Should().Be(s, "no gaps are introduced for identical inputs");
        r.AlignedSequence2.Should().Be(s);
    }

    /// <summary>
    /// Interaction witness: a single insertion is recovered by exactly one gap — the optimum
    /// aligns the shared bases and pays one gap penalty rather than cascading mismatches.
    /// </summary>
    [Test]
    public void AlignGlobal_SingleInsertion_RecoveredByOneGap()
    {
        string a = "ACGTACGT";          // length 8
        string b = "ACGTAACGT";         // one 'A' inserted after position 5 (length 9)
        var r = SequenceAligner.GlobalAlign(a, b, new ScoringMatrix(2, -1, -2, -2));

        (r.AlignedSequence1 + r.AlignedSequence2).Count(c => c == '-').Should().Be(1, "exactly one gap column");
        r.Score.Should().Be(8 * 2 + 1 * (-2), "8 matches and one gap");
    }

    /// <summary>
    /// Interaction witness: a punishing gap penalty drives equal-length sequences to a
    /// gapless diagonal — any indel pair would cost more than a mismatch.
    /// </summary>
    [Test]
    public void AlignGlobal_HarshGap_ForcesGaplessDiagonal()
    {
        string a = "ACGTACGT";
        string b = "TGCATGCA"; // equal length, fully divergent
        var r = SequenceAligner.GlobalAlign(a, b, new ScoringMatrix(1, -1, -100, -100));

        r.AlignedSequence1.Should().NotContain("-");
        r.AlignedSequence2.Should().NotContain("-");
        r.AlignedSequence1.Length.Should().Be(a.Length, "a gapless diagonal has one column per position");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ALIGN-LOCAL-001 — Local (Smith-Waterman) alignment (Alignment)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 36.
    // Spec: tests/TestSpecs/ALIGN-LOCAL-001.md (canonical SequenceAligner.LocalAlign).
    // Dimensions: matchScore(3) × mismatchPen(3) × gapPen(3) × seqLen(3). Grid 3⁴ = 81.
    //
    // Model (Smith & Waterman 1981): like Needleman-Wunsch but with a zero floor in the
    // recurrence, H(i,j)=max(0, H(i−1,j−1)+s, H(i−1,j)+d, H(i,j−1)+d); the score is the
    // matrix maximum and the alignment is traced back from it to the first zero. The optimal
    // local score is therefore always ≥ 0, the returned rows are substrings of the inputs in
    // [StartPosition, EndPosition], and their column score equals the reported Score.
    //
    // The combinatorial point: the three weights and the length jointly determine the best
    // local segment; the reported score is checked against an independent SW DP for every
    // weight combination, and the traceback's column score must reproduce it.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void AlignLocal_OptimalNonNegativeAndConsistent_AcrossScoringAndLength(
        [Values(1, 2, 5)] int matchScore,
        [Values(-1, -3, 0)] int mismatchPen,
        [Values(-1, -2, -5)] int gapPen,
        [Values(4, 8, 12)] int seqLen)
    {
        string s1 = DiverseDna(seqLen, 0x1234u);
        string s2 = DiverseDna(seqLen, 0x9ABCu);
        var scoring = new ScoringMatrix(matchScore, mismatchPen, GapOpen: gapPen, GapExtend: gapPen);

        var r = SequenceAligner.LocalAlign(s1, s2, scoring);

        r.Score.Should().BeGreaterThanOrEqualTo(0, "local alignment has a zero floor");
        r.Score.Should().Be(SwOptimalScore(s1, s2, matchScore, mismatchPen, gapPen),
            "the score is the Smith-Waterman optimum for these weights");
        r.AlignmentType.Should().Be(AlignmentType.Local);

        ScoreColumns(r.AlignedSequence1, r.AlignedSequence2, matchScore, mismatchPen, gapPen)
            .Should().Be(r.Score, "the reported score is the column score of the returned local alignment");

        if (r.Score > 0)
        {
            // The gap-free rows are the reported substrings of each input.
            r.AlignedSequence1.Replace("-", "")
                .Should().Be(s1.Substring(r.StartPosition1, r.EndPosition1 - r.StartPosition1 + 1));
            r.AlignedSequence2.Replace("-", "")
                .Should().Be(s2.Substring(r.StartPosition2, r.EndPosition2 - r.StartPosition2 + 1));
        }
    }

    /// <summary>
    /// Worked example: an embedded common substring is recovered exactly — the best local
    /// alignment is that substring with score length × matchScore and no gaps.
    /// </summary>
    [Test]
    public void AlignLocal_EmbeddedCommonSubstring_RecoveredExactly()
    {
        string a = "TTTT" + "ACGTACGT" + "TTTT";
        string b = "GGGG" + "ACGTACGT" + "GGGG";
        var r = SequenceAligner.LocalAlign(a, b, new ScoringMatrix(2, -1, -2, -2));

        r.Score.Should().Be(8 * 2, "the 8-nt common core scores 8×match with no gaps");
        r.AlignedSequence1.Should().Be("ACGTACGT");
        r.AlignedSequence2.Should().Be("ACGTACGT");
    }

    /// <summary>
    /// Interaction witness: with no shared bases the zero floor yields score 0 (an empty
    /// local alignment) rather than a negative score.
    /// </summary>
    [Test]
    public void AlignLocal_NoSharedBases_ScoresZero()
    {
        var r = SequenceAligner.LocalAlign("AAAA", "CCCC", new ScoringMatrix(1, -1, -2, -2));
        r.Score.Should().Be(0);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ALIGN-SEMI-001 — Semi-global (fitting) alignment (Alignment)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 37.
    // Spec: tests/TestSpecs/ALIGN-SEMI-001.md (canonical SequenceAligner.SemiGlobalAlign).
    // Dimensions: matchScore(3) × mismatchPen(3) × gapPen(3). Grid 3³ = 27.
    //
    // Model (fitting alignment): the query (seq1) is aligned end-to-end while leading and
    // trailing gaps in the reference (seq2) are FREE — DP first row is 0 and the optimum is
    // the maximum of the last row. So the reported Score excludes the free reference flanks
    // (it is NOT the raw column score of the whole alignment), and equals the fitting optimum.
    //
    // The combinatorial point: the three weights jointly determine where the query fits and
    // the optimal score; correctness is checked against an independent fitting DP for every
    // weight combination, with both inputs reconstructed and the trimmed core score matching.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void AlignSemi_FittingOptimumAndReconstruction_AcrossScoring(
        [Values(1, 2, 5)] int matchScore,
        [Values(-1, -3, 0)] int mismatchPen,
        [Values(-1, -2, -5)] int gapPen)
    {
        string query = DiverseDna(8, 0x1234u);
        string reference = DiverseDna(20, 0x9ABCu);
        var scoring = new ScoringMatrix(matchScore, mismatchPen, GapOpen: gapPen, GapExtend: gapPen);

        var r = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), scoring);

        r.AlignmentType.Should().Be(AlignmentType.SemiGlobal);
        r.AlignedSequence1.Replace("-", "").Should().Be(query, "the query is aligned end-to-end");
        r.AlignedSequence2.Replace("-", "").Should().Be(reference, "the reference is fully represented");
        r.AlignedSequence1.Length.Should().Be(r.AlignedSequence2.Length);

        r.Score.Should().Be(SemiGlobalOptimalScore(query, reference, matchScore, mismatchPen, gapPen),
            "the score is the fitting optimum (free reference end-gaps)");

        // Trimming the free leading/trailing reference flanks (query-side gaps at the ends)
        // leaves a core whose column score equals the reported Score.
        string a1 = r.AlignedSequence1, a2 = r.AlignedSequence2;
        int lo = 0; while (lo < a1.Length && a1[lo] == '-') lo++;
        int hi = a1.Length - 1; while (hi >= 0 && a1[hi] == '-') hi--;
        ScoreColumns(a1[lo..(hi + 1)], a2[lo..(hi + 1)], matchScore, mismatchPen, gapPen)
            .Should().Be(r.Score, "the core (non-free) columns reproduce the score");
    }

    /// <summary>
    /// Worked example: a query embedded in a longer reference fits exactly — leading and
    /// trailing reference flanks are free, so the score is queryLength × matchScore.
    /// </summary>
    [Test]
    public void AlignSemi_QueryEmbeddedInReference_FreeEndGaps()
    {
        string query = "ACGTACGT";
        string reference = "TTTT" + query + "TTTT";
        var r = SequenceAligner.SemiGlobalAlign(new DnaSequence(query), new DnaSequence(reference), new ScoringMatrix(2, -1, -2, -2));

        r.Score.Should().Be(query.Length * 2, "the free flanks are not penalised");
        r.AlignedSequence1.Replace("-", "").Should().Be(query);
        r.AlignedSequence1.Should().StartWith("-", "the leading reference flank is a free gap in the query row");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ALIGN-MULTI-001 — Multiple sequence alignment (Alignment)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 38.
    // Spec: tests/TestSpecs/ALIGN-MULTI-001.md (canonical MultipleAlign / MultipleAlignProgressive).
    // Dimensions: nSeqs(3) × seqLen(3) × gapPen(3) × guideTree(2). Grid 3×3×3×2 = 54.
    //
    // Model (Wikipedia MSA): an MSA arranges k sequences into a rectangular block by inserting
    // gaps only (never editing residues); the score is the sum-of-pairs (SP) over all C(k,2)
    // row pairs, with gap-gap columns neutral. guideTree(2) selects the star anchor method
    // (MultipleAlign) vs the progressive guide-tree method (MultipleAlignProgressive).
    //
    // The combinatorial point: regardless of nSeqs, length, gap weight or assembly strategy,
    // every result must (a) keep all k rows, (b) be rectangular, (c) reproduce each input when
    // gaps are stripped, and (d) carry a TotalScore equal to the SP score of the returned block.
    // ═══════════════════════════════════════════════════════════════════════

    public enum MsaMethod { StarAnchor, ProgressiveTree }

    private static char NextBase(char c) => "ACGT"["ACGT".IndexOf(c) is var k && k >= 0 ? (k + 1) % 4 : 0];

    /// <summary>A deterministic family of related sequences; one member carries a deletion (forces gaps).</summary>
    private static List<DnaSequence> MakeFamily(int nSeqs, int seqLen, uint seed)
    {
        string baseSeq = DiverseDna(seqLen, seed);
        var family = new List<DnaSequence>();
        for (int i = 0; i < nSeqs; i++)
        {
            var chars = baseSeq.ToCharArray();
            for (int k = 0; k < i; k++) { int p = (k * 5 + i) % chars.Length; chars[p] = NextBase(chars[p]); }
            string s = new string(chars);
            if (i == 1) s = s.Remove(seqLen / 2, 1); // a single-base deletion in one member
            family.Add(new DnaSequence(s));
        }
        return family;
    }

    /// <summary>Independent sum-of-pairs score mirroring the documented MSA convention (gap-gap neutral).</summary>
    private static int SumOfPairs(string[] rows, int match, int mismatch, int gap)
    {
        int len = rows.Max(r => r.Length), total = 0;
        for (int pos = 0; pos < len; pos++)
            for (int i = 0; i < rows.Length; i++)
            {
                char ci = pos < rows[i].Length ? rows[i][pos] : '-';
                for (int j = i + 1; j < rows.Length; j++)
                {
                    char cj = pos < rows[j].Length ? rows[j][pos] : '-';
                    if (ci == '-' && cj == '-') continue;
                    else if (ci == '-' || cj == '-') total += gap;
                    else total += ci == cj ? match : mismatch;
                }
            }
        return total;
    }

    [Test, Combinatorial]
    public void AlignMulti_RectangularReconstructsInputs_AndScoreIsSumOfPairs(
        [Values(3, 4, 5)] int nSeqs,
        [Values(8, 12, 20)] int seqLen,
        [Values(-1, -2, -5)] int gapPen,
        [Values(MsaMethod.StarAnchor, MsaMethod.ProgressiveTree)] MsaMethod method)
    {
        var family = MakeFamily(nSeqs, seqLen, 0x77u);
        var scoring = new ScoringMatrix(1, -1, GapOpen: gapPen, GapExtend: gapPen);

        var result = method == MsaMethod.StarAnchor
            ? SequenceAligner.MultipleAlign(family, scoring)
            : SequenceAligner.MultipleAlignProgressive(family, scoring);

        result.AlignedSequences.Should().HaveCount(nSeqs, "every input keeps a row");
        int width = result.AlignedSequences[0].Length;
        result.AlignedSequences.Should().OnlyContain(r => r.Length == width, "the block is rectangular");
        result.Consensus.Length.Should().Be(width, "the consensus spans the alignment width");

        result.AlignedSequences.Select(r => r.Replace("-", ""))
            .Should().BeEquivalentTo(family.Select(f => f.Sequence), "gaps-removed rows reproduce the inputs");

        result.TotalScore.Should().Be(SumOfPairs(result.AlignedSequences, 1, -1, gapPen),
            "TotalScore is the sum-of-pairs score of the returned block");
    }

    /// <summary>
    /// Interaction witness: identical sequences need no gaps under either strategy — the block
    /// is the sequence repeated, the consensus is the sequence, and SP = C(k,2)·L·match.
    /// </summary>
    [Test, Combinatorial]
    public void AlignMulti_IdenticalSequences_TrivialBlock(
        [Values(3, 4)] int nSeqs,
        [Values(MsaMethod.StarAnchor, MsaMethod.ProgressiveTree)] MsaMethod method)
    {
        string s = DiverseDna(12, 0x2BADu);
        var family = Enumerable.Repeat(s, nSeqs).Select(x => new DnaSequence(x)).ToList();
        var scoring = new ScoringMatrix(1, -1, -2, -2);

        var result = method == MsaMethod.StarAnchor
            ? SequenceAligner.MultipleAlign(family, scoring)
            : SequenceAligner.MultipleAlignProgressive(family, scoring);

        result.AlignedSequences.Should().OnlyContain(r => r == s, "no gaps for identical inputs");
        result.Consensus.Should().Be(s);
        result.TotalScore.Should().Be(nSeqs * (nSeqs - 1) / 2 * s.Length, "C(k,2)·L matches at +1 each");
    }

    /// <summary>
    /// Worked example: a single divergent column resolves to the majority base in the consensus.
    /// </summary>
    [Test]
    public void AlignMulti_Consensus_IsColumnMajority()
    {
        var family = new[] { "ACGT", "ACGT", "ACTT" }.Select(s => new DnaSequence(s)).ToList();
        var result = SequenceAligner.MultipleAlign(family, SequenceAligner.SimpleDna);

        result.AlignedSequences.Should().OnlyContain(r => r.Length == 4, "no gaps needed");
        result.Consensus.Should().Be("ACGT", "column 3 majority is G (2 of 3)");
    }
}
