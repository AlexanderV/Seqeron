#pragma warning disable S3776
namespace Seqeron.Genomics.Analysis;

/// <summary>The conservation score of a column (Capra and Singh 2007, <see cref="ConservationScores"/>).</summary>
public enum ConservationMethod
{
    /// <summary>Jensen-Shannon divergence from a background distribution (the default method of the script).</summary>
    JsDivergence,

    /// <summary>The Shannon information score <c>1 - h</c>.</summary>
    ShannonEntropy,

    /// <summary>The information score over the groups of Mirny and Shakhnovich.</summary>
    PropertyEntropy,

    /// <summary>Relative entropy (Kullback-Leibler divergence from the background).</summary>
    RelativeEntropy,

    /// <summary>The sum-of-pairs score with a similarity matrix (BLOSUM62 by default).</summary>
    SumOfPairs,
}

public static partial class MultipleAlignmentStatistics
{
    /// <summary>
    /// The conservation score of every column of aligned rows as the main program of <c>score_conservation.py</c> (Capra and Singh 2007) computes it: the rows are upper-cased and read as the script's alignment readers do (B is D, Z is Q,
    /// X and '.' are gaps), each column is scored with <paramref name="method"/> (sequence weights of Henikoff and Henikoff when <paramref name="sequenceWeights"/>), and a window of <paramref name="windowSize"/> columns on either
    /// side smooths the profile (<see cref="ConservationScores.WindowScore"/>; 0 = none).
    /// </summary>
    /// <param name="alignedRows">Aligned rows of equal length.</param>
    /// <param name="method">The score of a column (default Jensen-Shannon divergence).</param>
    /// <param name="background">Background distribution of 20 or 21 values in the order <see cref="ConservationScores.AminoAcids"/> (default BLOSUM62; JSD and relative entropy).</param>
    /// <param name="similarity">20 x 20 similarity matrix in the order of <see cref="ConservationScores.AminoAcids"/> (default BLOSUM62; sum of pairs).</param>
    /// <param name="sequenceWeights">Weight the rows (Henikoff and Henikoff 1994).</param>
    /// <param name="gapPenalty">Multiply each score by the weighted fraction of non-gap rows (default true).</param>
    /// <param name="windowSize">Window length on either side (default 0).</param>
    public static double[] ConservationProfile(IReadOnlyList<string> alignedRows, ConservationMethod method = ConservationMethod.JsDivergence, IReadOnlyList<double>? background = null, double[,]? similarity = null,
        bool sequenceWeights = false, bool gapPenalty = true, int windowSize = 0)
    {
        ArgumentNullException.ThrowIfNull(alignedRows);
        if (alignedRows.Count == 0)
            throw new ArgumentException("At least one aligned row is required.", nameof(alignedRows));
        int length = alignedRows[0].Length;
        if (alignedRows.Any(r => r.Length != length))
            throw new ArgumentException("Rows of unequal length.", nameof(alignedRows));
        if (windowSize < 0)
            throw new ArgumentOutOfRangeException(nameof(windowSize), "The window size is not negative.");
        var rows = alignedRows.Select(r => r.ToUpperInvariant().Replace('B', 'D').Replace('Z', 'Q').Replace('X', '-').Replace('.', '-')).ToList();
        var weights = sequenceWeights ? ConservationScores.SequenceWeights(rows) : Enumerable.Repeat(1.0, rows.Count).ToArray();
        var bg = background ?? ConservationScores.BlosumBackground;
        var scores = new double[length];
        for (int i = 0; i < length; i++)
        {
            string column = new([.. rows.Select(r => r[i])]);
            scores[i] = method switch
            {
                ConservationMethod.ShannonEntropy => ConservationScores.ShannonEntropy(column, weights, gapPenalty),
                ConservationMethod.PropertyEntropy => ConservationScores.PropertyEntropy(column, weights, gapPenalty),
                ConservationMethod.RelativeEntropy => ConservationScores.RelativeEntropy(column, bg, weights, gapPenalty),
                ConservationMethod.SumOfPairs => ConservationScores.SumOfPairs(column, similarity, weights, gapPenalty),
                _ => ConservationScores.JsDivergence(column, bg, weights, gapPenalty),
            };
        }

        return windowSize > 0 ? ConservationScores.WindowScore(scores, windowSize) : scores;
    }
}
