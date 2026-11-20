namespace Seqeron.Genomics.Alignment;

// B09 A3-10: per-direction and per-end gap scores (Biopython PairwiseAligner's twelve {target,query}_{internal,left,right}_{open,extend}_gap_score properties) on the
// existing linear / affine dynamic programs. There is no second program: the gap scores of the one Needleman-Wunsch / Gotoh / Smith-Waterman core become functions
// of the matrix position (PairScoring.XExt / XOpen / YExt / YOpen): a vertical run (gap in sequence 2) stays in one column and a horizontal run (gap in sequence 1) in
// one row, so a run in column / row 0 is a start run, in the last column / row an end run, anywhere else interior.
public static partial class SequenceAligner
{
    #region Gap profiles (per-direction / per-end gap scores)

    /// <summary>
    /// Linear-gap global alignment with the direction- and position-specific gap scores of a <see cref="GapScoreProfile"/> (Biopython
    /// <c>left_</c> / <c>right_</c> / <c>internal_</c> × <c>target_</c> / <c>query_</c> gap scores; <see cref="GapScore.Extend"/> is the score per gap position, <see cref="GapScore.Open"/> is
    /// ignored as <see cref="ScoringMatrix.GapOpen"/> is). Members left null use <see cref="ScoringMatrix.GapExtend"/>; <see cref="GapScore.Free"/> members are free overhangs (the
    /// <see cref="EndGaps"/> policies, and any mixture of free and priced ends). Ties are broken diagonal &gt; up &gt; left, as in <see cref="GlobalAlign(string, string, ScoringMatrix?)"/>.
    /// </summary>
    /// <remarks>
    /// Coordinates: a free (zero-score) start / end run is an overhang and is excluded from <see cref="AlignmentResult.StartPosition1"/> … <see cref="AlignmentResult.EndPosition2"/> (the
    /// semi-global convention); the rows still show it against gaps. A priced end run belongs to the alignment (coordinates span the full sequence).
    /// Verified against Biopython 1.88 <c>PairwiseAligner(mode = "global")</c> on 600 random schemes (docs/Evidence/ALIGN-SEMI-001-Evidence.md). Sequences longer than the
    /// full-matrix limit (<c>LinearSpaceCellThreshold</c> cells) are refused: the linear-space program has one gap model.
    /// </remarks>
    /// <param name="sequence1">First sequence (Biopython target).</param>
    /// <param name="sequence2">Second sequence (Biopython query).</param>
    /// <param name="scoring">Match / mismatch scoring (default <see cref="SimpleDna"/>).</param>
    /// <param name="gapScores">The gap scores.</param>
    /// <exception cref="ArgumentNullException">A sequence or <paramref name="gapScores"/> is null.</exception>
    /// <exception cref="ArgumentException">A gap score is positive / not finite / finer than 1e-6, or the input is beyond the full-matrix limit.</exception>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult GlobalAlign(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, scoring, null, gapScores, affine: false, AlignmentType.Global);

    /// <summary>Affine-gap <see cref="GlobalAlign(string, string, ScoringMatrix?, GapScoreProfile)"/> (a gap of length k scores <see cref="GapScore.Open"/> + k·<see cref="GapScore.Extend"/>; ties M &gt; X &gt; Y).</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult GlobalAlignAffine(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, scoring, null, gapScores, affine: true, AlignmentType.Global);

    /// <summary><see cref="GlobalAlign(string, string, ScoringMatrix?, GapScoreProfile)"/> with a substitution matrix and fractional gap scores.</summary>
    public static AlignmentResult GlobalAlign(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, null, scoring, gapScores, affine: false, AlignmentType.Global);

    /// <summary><see cref="GlobalAlignAffine(string, string, ScoringMatrix?, GapScoreProfile)"/> with a substitution matrix and fractional gap scores.</summary>
    public static AlignmentResult GlobalAlignAffine(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, null, scoring, gapScores, affine: true, AlignmentType.Global);

    /// <summary>
    /// <see cref="GlobalAlign(string, string, ScoringMatrix?, GapScoreProfile)"/> for the semi-global family: the same dynamic program, the result tagged
    /// <see cref="AlignmentType.SemiGlobal"/> when a start / end member is <see cref="GapScore.Free"/> (<see cref="GapScoreProfile.FromEndGaps"/> reproduces every
    /// <see cref="EndGaps"/> policy; the rows are one co-optimal alignment, not necessarily the tie-break of <see cref="SemiGlobalAlign(string, string, ScoringMatrix?, EndGaps)"/>).
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlign(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, scoring, null, gapScores, affine: false, AlignmentType.SemiGlobal);

    /// <summary>Affine-gap <see cref="SemiGlobalAlign(string, string, ScoringMatrix?, GapScoreProfile)"/>.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult SemiGlobalAlignAffine(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, scoring, null, gapScores, affine: true, AlignmentType.SemiGlobal);

    /// <summary><see cref="SemiGlobalAlign(string, string, ScoringMatrix?, GapScoreProfile)"/> with a substitution matrix and fractional gap scores.</summary>
    public static AlignmentResult SemiGlobalAlign(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, null, scoring, gapScores, affine: false, AlignmentType.SemiGlobal);

    /// <summary><see cref="SemiGlobalAlignAffine(string, string, ScoringMatrix?, GapScoreProfile)"/> with a substitution matrix and fractional gap scores.</summary>
    public static AlignmentResult SemiGlobalAlignAffine(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, null, scoring, gapScores, affine: true, AlignmentType.SemiGlobal);

    /// <summary>
    /// Local (Smith-Waterman) alignment with direction-specific <b>interior</b> gap scores (<see cref="GapScoreProfile.Interior1"/> / <see cref="GapScoreProfile.Interior2"/>: Biopython
    /// <c>query_internal_*</c> / <c>target_internal_*</c>); a local alignment never starts or ends in a gap, so the start / end members have no effect (Biopython's local mode ignores
    /// the left / right scores the same way).
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult LocalAlign(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, scoring, null, gapScores, affine: false, AlignmentType.Local);

    /// <summary>Affine-gap <see cref="LocalAlign(string, string, ScoringMatrix?, GapScoreProfile)"/>.</summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static AlignmentResult LocalAlignAffine(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, scoring, null, gapScores, affine: true, AlignmentType.Local);

    /// <summary><see cref="LocalAlign(string, string, ScoringMatrix?, GapScoreProfile)"/> with a substitution matrix and fractional gap scores.</summary>
    public static AlignmentResult LocalAlign(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, null, scoring, gapScores, affine: false, AlignmentType.Local);

    /// <summary><see cref="LocalAlignAffine(string, string, ScoringMatrix?, GapScoreProfile)"/> with a substitution matrix and fractional gap scores.</summary>
    public static AlignmentResult LocalAlignAffine(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores) =>
        ProfileAlign(sequence1, sequence2, null, scoring, gapScores, affine: true, AlignmentType.Local);

    private static AlignmentResult ProfileAlign(
        string sequence1, string sequence2, ScoringMatrix? scalar, MatrixScoring? matrix, GapScoreProfile profile, bool affine, AlignmentType type)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(profile);

        string a = sequence1.ToUpperInvariant();
        string b = sequence2.ToUpperInvariant();
        PairScoring ps = matrix is not null
            ? PairScoring.Create(matrix, a, b, profile)
            : PairScoring.Create(scalar ?? SimpleDna, profile);

        bool local = type == AlignmentType.Local;
        AlignmentResult r;
        if (affine)
            r = AffineAlignCore(a, b, ps, local ? AlignmentType.Local : AlignmentType.Global);
        else if (local)
            r = LocalAlignCore(a, b, ps);
        else
            r = LinearAlignCore(a, b, ps, EndGaps.None);

        if (!local)
            r = FreeOverhangCoordinates(r, profile, profile.HasFreeEnd ? AlignmentType.SemiGlobal : AlignmentType.Global);

        double exact = (double)r.Score / ps.Scale;
        return ps.Scale == 1 ? r : r with { Score = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactScore = exact };
    }

    /// <summary>Tags the result and excludes a FREE (zero-score) start / end run from the coordinates (the semi-global overhang convention); a priced run stays inside.</summary>
    private static AlignmentResult FreeOverhangCoordinates(AlignmentResult r, GapScoreProfile profile, AlignmentType tag)
    {
        string row1 = r.AlignedSequence1, row2 = r.AlignedSequence2;
        int cols = row1.Length;
        int lead1 = 0, lead2 = 0, trail1 = 0, trail2 = 0;
        // leading run of columns with a gap in sequence 2 (sequence-1 residues) / in sequence 1
        while (lead1 < cols && row2[lead1] == '-') lead1++;
        while (lead2 < cols && row1[lead2] == '-') lead2++;
        while (trail1 < cols && row2[cols - 1 - trail1] == '-') trail1++;
        while (trail2 < cols && row1[cols - 1 - trail2] == '-') trail2++;
        bool IsFree(GapScore? g) => g is { Open: 0, Extend: 0 };
        int s1 = IsFree(profile.Start1) ? lead1 : 0;
        int s2 = IsFree(profile.Start2) ? lead2 : 0;
        int e1 = IsFree(profile.End1) ? trail1 : 0;
        int e2 = IsFree(profile.End2) ? trail2 : 0;
        int len1 = row1.Count(c => c != '-'), len2 = row2.Count(c => c != '-');
        return r with
        {
            AlignmentType = tag,
            StartPosition1 = s1,
            StartPosition2 = s2,
            EndPosition1 = len1 - e1 - 1,
            EndPosition2 = len2 - e2 - 1,
        };
    }

    #endregion
}
