using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// R14-10 (B09 F303): every co-optimal alignment (the set Biopython's PairwiseAligner.align iterates) under the per-direction / per-end gap scores of a GapScoreProfile. The enumerators of the global and local problems
// take the PairScoring that carries the profile (the very scoring of GlobalAlign / LocalAlign with a profile), so the set is the one of the dynamic program that finds the single alignment.
public static partial class SequenceAligner
{
    private static IReadOnlyList<AlignmentResult> CoOptimalProfile(string sequence1, string sequence2, ScoringMatrix? scalar, MatrixScoring? matrix, GapScoreProfile profile, bool affine, bool local, int maxAlignments)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(profile);
        string a = sequence1.ToUpperInvariant(), b = sequence2.ToUpperInvariant();
        PairScoring ps = matrix is not null ? PairScoring.Create(matrix, a, b, profile) : PairScoring.Create(scalar ?? SimpleDna, profile);
        IReadOnlyList<AlignmentResult> raw = local
            ? LocalCoOptimalCore(a, b, ps, affine, maxAlignments)
            : GlobalCoOptimalCore(a, b, ps, affine, FreeEnds(profile), maxAlignments);
        var list = new List<AlignmentResult>(raw.Count);
        foreach (var r in raw)
        {
            var x = r;
            if (!local)
                x = FreeOverhangCoordinates(r, profile, profile.HasFreeEnd ? AlignmentType.SemiGlobal : AlignmentType.Global);
            if (ps.Scale != 1)
                x = WithExactScore(x, ps);
            list.Add(x);
        }

        return list;
    }

    // the free (zero-score) start / end members of a profile as the end-gap flags of the enumerator
    private static EndGaps FreeEnds(GapScoreProfile profile)
    {
        static bool Free(GapScore? g) => g is { Open: 0, Extend: 0 };
        var ends = EndGaps.None;
        if (Free(profile.Start1))
            ends |= EndGaps.FreeStart1;
        if (Free(profile.End1))
            ends |= EndGaps.FreeEnd1;
        if (Free(profile.Start2))
            ends |= EndGaps.FreeStart2;
        if (Free(profile.End2))
            ends |= EndGaps.FreeEnd2;
        return ends;
    }

    /// <summary>Every co-optimal linear-gap global alignment under the gap scores of a <see cref="GapScoreProfile"/> (see <see cref="GlobalAlign(string, string, ScoringMatrix?, GapScoreProfile)"/>); capped at <paramref name="maxAlignments"/>.</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignCoOptimal(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores, int maxAlignments = 1000) =>
        CoOptimalProfile(sequence1, sequence2, scoring, null, gapScores, affine: false, local: false, maxAlignments);

    /// <summary>Every co-optimal affine-gap global alignment under the gap scores of a <see cref="GapScoreProfile"/>.</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignAffineCoOptimal(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores, int maxAlignments = 1000) =>
        CoOptimalProfile(sequence1, sequence2, scoring, null, gapScores, affine: true, local: false, maxAlignments);

    /// <summary><see cref="GlobalAlignCoOptimal(string, string, ScoringMatrix?, GapScoreProfile, int)"/> with a substitution matrix.</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignCoOptimal(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores, int maxAlignments = 1000) =>
        CoOptimalProfile(sequence1, sequence2, null, scoring, gapScores, affine: false, local: false, maxAlignments);

    /// <summary><see cref="GlobalAlignAffineCoOptimal(string, string, ScoringMatrix?, GapScoreProfile, int)"/> with a substitution matrix.</summary>
    public static IReadOnlyList<AlignmentResult> GlobalAlignAffineCoOptimal(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores, int maxAlignments = 1000) =>
        CoOptimalProfile(sequence1, sequence2, null, scoring, gapScores, affine: true, local: false, maxAlignments);

    /// <summary>Every co-optimal linear-gap local alignment under the interior gap scores of a <see cref="GapScoreProfile"/> (the start / end members have no effect on local alignments).</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignCoOptimal(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores, int maxAlignments = 1000) =>
        CoOptimalProfile(sequence1, sequence2, scoring, null, gapScores, affine: false, local: true, maxAlignments);

    /// <summary>Every co-optimal affine-gap local alignment under the interior gap scores of a <see cref="GapScoreProfile"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineCoOptimal(string sequence1, string sequence2, ScoringMatrix? scoring, GapScoreProfile gapScores, int maxAlignments = 1000) =>
        CoOptimalProfile(sequence1, sequence2, scoring, null, gapScores, affine: true, local: true, maxAlignments);

    /// <summary><see cref="LocalAlignCoOptimal(string, string, ScoringMatrix?, GapScoreProfile, int)"/> with a substitution matrix.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignCoOptimal(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores, int maxAlignments = 1000) =>
        CoOptimalProfile(sequence1, sequence2, null, scoring, gapScores, affine: false, local: true, maxAlignments);

    /// <summary><see cref="LocalAlignAffineCoOptimal(string, string, ScoringMatrix?, GapScoreProfile, int)"/> with a substitution matrix.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineCoOptimal(string sequence1, string sequence2, MatrixScoring scoring, GapScoreProfile gapScores, int maxAlignments = 1000) =>
        CoOptimalProfile(sequence1, sequence2, null, scoring, gapScores, affine: true, local: true, maxAlignments);
}
