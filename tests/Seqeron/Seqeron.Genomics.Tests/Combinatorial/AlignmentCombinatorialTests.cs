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
}
