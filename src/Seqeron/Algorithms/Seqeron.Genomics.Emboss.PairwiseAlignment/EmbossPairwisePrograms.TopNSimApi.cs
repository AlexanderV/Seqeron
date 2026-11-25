// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.PairwiseAlignment;

public static partial class EmbossPairwisePrograms
{
    #region top-N local alignments in the order of EMBOSS matcher (Huang-Miller SIM)

    // the sequences of matcher are read with the EMBOSS reader: the gap characters and every other non-letter vanish; with the Waterman-Eggert order the MIT aligner is called as it is.

    /// <summary>
    /// <see cref="SequenceAligner.LocalAlignTopN(string, string, ScoringMatrix?, int, int)"/> with the choice of the order of co-optimal alignments (<paramref name="tieOrder"/>; the default is the order of EMBOSS <c>matcher</c>, <see cref="TopNTieOrder.Sim"/>).
    /// </summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0, TopNTieOrder tieOrder = TopNTieOrder.Sim)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return tieOrder switch
        {
            TopNTieOrder.WatermanEggert => SequenceAligner.LocalAlignTopN(sequence1, sequence2, scoring, maxAlignments, minScore),
            TopNTieOrder.Sim => SimTopN(sequence1, sequence2, scoring ?? SequenceAligner.SimpleDna, affine: false, maxAlignments, minScore),
            _ => throw new ArgumentOutOfRangeException(nameof(tieOrder), tieOrder, "Unknown tie order."),
        };
    }

    /// <summary><see cref="LocalAlignTopN(string, string, ScoringMatrix?, int, int, TopNTieOrder)"/> on <see cref="DnaSequence"/>s.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0, TopNTieOrder tieOrder = TopNTieOrder.Sim)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalAlignTopN(sequence1.Sequence, sequence2.Sequence, scoring, maxAlignments, minScore, tieOrder);
    }

    /// <summary><see cref="LocalAlignTopN(string, string, ScoringMatrix?, int, int, TopNTieOrder)"/> that observes <paramref name="cancellationToken"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(string sequence1, string sequence2, ScoringMatrix? scoring, int maxAlignments, int minScore, TopNTieOrder tieOrder, CancellationToken cancellationToken)
    {
        using var cancelScope = AlignmentCancellation.Enter(cancellationToken);
        return LocalAlignTopN(sequence1, sequence2, scoring, maxAlignments, minScore, tieOrder);
    }

    /// <summary>Affine-gap twin of <see cref="LocalAlignTopN(string, string, ScoringMatrix?, int, int, TopNTieOrder)"/> (a gap of length k scores GapOpen + k * GapExtend).</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(string sequence1, string sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0, TopNTieOrder tieOrder = TopNTieOrder.Sim)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return tieOrder switch
        {
            TopNTieOrder.WatermanEggert => SequenceAligner.LocalAlignAffineTopN(sequence1, sequence2, scoring, maxAlignments, minScore),
            TopNTieOrder.Sim => SimTopN(sequence1, sequence2, scoring ?? SequenceAligner.SimpleDna, affine: true, maxAlignments, minScore),
            _ => throw new ArgumentOutOfRangeException(nameof(tieOrder), tieOrder, "Unknown tie order."),
        };
    }

    /// <summary><see cref="LocalAlignAffineTopN(string, string, ScoringMatrix?, int, int, TopNTieOrder)"/> on <see cref="DnaSequence"/>s.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(DnaSequence sequence1, DnaSequence sequence2, ScoringMatrix? scoring = null, int maxAlignments = 10, int minScore = 0, TopNTieOrder tieOrder = TopNTieOrder.Sim)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LocalAlignAffineTopN(sequence1.Sequence, sequence2.Sequence, scoring, maxAlignments, minScore, tieOrder);
    }

    /// <summary><see cref="LocalAlignAffineTopN(string, string, ScoringMatrix?, int, int, TopNTieOrder)"/> that observes <paramref name="cancellationToken"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(string sequence1, string sequence2, ScoringMatrix? scoring, int maxAlignments, int minScore, TopNTieOrder tieOrder, CancellationToken cancellationToken)
    {
        using var cancelScope = AlignmentCancellation.Enter(cancellationToken);
        return LocalAlignAffineTopN(sequence1, sequence2, scoring, maxAlignments, minScore, tieOrder);
    }

    /// <summary>Top-N local alignments with a substitution matrix and linear gap scores in the chosen order (<see cref="AlignmentResult.ExactScore"/> = the real score).</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0, TopNTieOrder tieOrder = TopNTieOrder.Sim) =>
        TopNWithMatrix(sequence1, sequence2, scoring, affine: false, maxAlignments, minScore, tieOrder);

    /// <summary>Top-N local alignments with a substitution matrix and affine gap scores in the chosen order.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments = 10, double minScore = 0, TopNTieOrder tieOrder = TopNTieOrder.Sim) =>
        TopNWithMatrix(sequence1, sequence2, scoring, affine: true, maxAlignments, minScore, tieOrder);

    /// <summary><see cref="LocalAlignAffineTopN(string, string, MatrixScoring, int, double, TopNTieOrder)"/> that observes <paramref name="cancellationToken"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignAffineTopN(string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments, double minScore, TopNTieOrder tieOrder, CancellationToken cancellationToken)
    {
        using var cancelScope = AlignmentCancellation.Enter(cancellationToken);
        return LocalAlignAffineTopN(sequence1, sequence2, scoring, maxAlignments, minScore, tieOrder);
    }

    /// <summary><see cref="LocalAlignTopN(string, string, MatrixScoring, int, double, TopNTieOrder)"/> that observes <paramref name="cancellationToken"/>.</summary>
    public static IReadOnlyList<AlignmentResult> LocalAlignTopN(string sequence1, string sequence2, MatrixScoring scoring, int maxAlignments, double minScore, TopNTieOrder tieOrder, CancellationToken cancellationToken)
    {
        using var cancelScope = AlignmentCancellation.Enter(cancellationToken);
        return LocalAlignTopN(sequence1, sequence2, scoring, maxAlignments, minScore, tieOrder);
    }

    private static IReadOnlyList<AlignmentResult> SimTopN(string sequence1, string sequence2, PairScoring scoring, bool affine, int maxAlignments, long minScaled) =>
        LocalTopNSimCore(EmbossSequenceReader.Clean(sequence1).ToUpperInvariant(), EmbossSequenceReader.Clean(sequence2).ToUpperInvariant(), scoring, affine, maxAlignments, minScaled);

    private static IReadOnlyList<AlignmentResult> TopNWithMatrix(string sequence1, string sequence2, MatrixScoring scoring, bool affine, int maxAlignments, double minScore, TopNTieOrder tieOrder)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentNullException.ThrowIfNull(scoring);
        if (double.IsNaN(minScore))
            throw new ArgumentOutOfRangeException(nameof(minScore), "minScore must not be NaN.");
        switch (tieOrder)
        {
            case TopNTieOrder.WatermanEggert:
                return affine ? SequenceAligner.LocalAlignAffineTopN(sequence1, sequence2, scoring, maxAlignments, minScore) : SequenceAligner.LocalAlignTopN(sequence1, sequence2, scoring, maxAlignments, minScore);
            case TopNTieOrder.Sim:
                break;
            default:
                throw new ArgumentOutOfRangeException(nameof(tieOrder), tieOrder, "Unknown tie order.");
        }

        // matcher reads its sequences with the EMBOSS reader (A4-3): the gap characters and every other non-letter vanish
        string a = EmbossSequenceReader.Clean(sequence1).ToUpperInvariant();
        string b = EmbossSequenceReader.Clean(sequence2).ToUpperInvariant();
        // SIM: q = gapopen - gapextend may be negative (matcher -gapopen below -gapextend)
        var ps = PairScoring.Create(scoring, a, b, null, allowPositiveGapOpen: affine);
        long minScaled = (long)Math.Ceiling((decimal)Math.Clamp(minScore, -1e12, 1e12) * ps.Scale);
        var raw = LocalTopNSimCore(a, b, ps, affine, maxAlignments, minScaled);
        var list = new List<AlignmentResult>(raw.Count);
        foreach (var r in raw)
        {
            double exact = (double)r.Score / ps.Scale;
            list.Add(r with { Score = (int)Math.Round(exact, MidpointRounding.AwayFromZero), ExactScore = exact });
        }

        return list;
    }

    #endregion
}
