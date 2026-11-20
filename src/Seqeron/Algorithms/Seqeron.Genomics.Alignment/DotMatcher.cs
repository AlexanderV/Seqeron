namespace Seqeron.Genomics.Alignment;

/// <summary>A segment of the plot of <see cref="DotMatcher.Compute"/>: the diagonal run of windows scoring at least the threshold, 1-based, inclusive, on both sequences.</summary>
/// <param name="Start1">First position on the first sequence.</param>
/// <param name="Start2">First position on the second sequence.</param>
/// <param name="End1">Last position on the first sequence.</param>
/// <param name="End2">Last position on the second sequence.</param>
public readonly record struct DotMatcherSegment(int Start1, int Start2, int End1, int End2);

/// <summary>
/// EMBOSS <c>dotmatcher</c> (6.6.0, <c>emboss/dotmatcher.c</c>): the dot plot of two sequences by windowed substitution-matrix scores. On every diagonal a window of <c>windowSize</c> residue pairs slides along; a
/// maximal run of windows whose summed matrix score reaches the threshold is drawn as one segment from the first residue of its first window to the last residue of its last window.
/// </summary>
public static class DotMatcher
{
    /// <summary>
    /// The segments of the plot in the order the program draws them. Measured on the binary: the program's scan runs over the SECOND sequence as its first one (the first diagonals start on the first residue of
    /// <paramref name="sequence2"/> and at residue j of <paramref name="sequence1"/>, j = 0 .. length1 - window - 1; then those that start at residue i of <paramref name="sequence2"/>, i = 0 .. length2 - window - 1)
    /// and the plot's x axis is <paramref name="sequence2"/>; the segments are reported here with the positions on <paramref name="sequence1"/> in <c>Start1</c> / <c>End1</c> (the plot's y values).
    /// </summary>
    /// <param name="sequence1">First sequence (x axis).</param>
    /// <param name="sequence2">Second sequence (y axis).</param>
    /// <param name="windowSize">Window size (<c>-windowsize</c>, default 10; at least 1).</param>
    /// <param name="threshold">Score a window must reach (<c>-threshold</c>, default 23).</param>
    /// <param name="matrix">Scoring matrix (<c>-matrixfile</c>); null = EDNAFULL when the first sequence consists of nucleotide codes only, EBLOSUM62 otherwise. A character that is not a letter of the matrix scores 0.</param>
    /// <returns>The segments (empty when a sequence is not longer than the window, where the program reads beyond the sequence).</returns>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> &lt; 1.</exception>
    public static IReadOnlyList<DotMatcherSegment> Compute(
        string sequence1, string sequence2, int windowSize = 10, int threshold = 23, SubstitutionMatrix? matrix = null)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        matrix ??= IsNucleotide(sequence1) ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;
        int len1 = sequence1.Length, len2 = sequence2.Length;
        var segments = new List<DotMatcherSegment>();
        if (len1 <= windowSize || len2 <= windowSize)
            return segments;

        // the program's s1 is the second sequence (index i, the plot's x) and its s2 the first (index k, the plot's y)
        int Score(int i, int k) => matrix.TryGetScore(sequence2[i], sequence1[k], out int s) ? s : 0;

        for (int j = 0; j < len1 - windowSize; j++)
            Diagonal(0, j);
        for (int i = 0; i < len2 - windowSize; i++)
            Diagonal(i, 0);
        return segments;

        // one diagonal from (i0, k0): the window sum is updated by one pair per step
        void Diagonal(int i0, int k0)
        {
            int i = i0, k = k0, total = 0, starti = 0, startj = 0;
            bool above = false;
            for (int l = 0; l < windowSize; l++)
                total += Score(i++, k++);
            if (total >= threshold)
            {
                above = true;
                starti = i - windowSize;
                startj = k - windowSize;
            }

            while (i < len2 && k < len1)
            {
                total = total - Score(i - windowSize, k - windowSize) + Score(i, k);
                if (above)
                {
                    if (total < threshold)
                    {
                        above = false;
                        segments.Add(new DotMatcherSegment(startj + 1, starti + 1, k, i));
                    }
                }
                else if (total >= threshold)
                {
                    starti = i - windowSize;
                    startj = k - windowSize;
                    above = true;
                }

                i++;
                k++;
            }

            if (above)
                segments.Add(new DotMatcherSegment(startj + 1, starti + 1, k, i));
        }
    }

    // ajSeqIsNuc of an untyped sequence: every character is a nucleotide code or a gap
    private static bool IsNucleotide(string sequence) =>
        sequence.All(c => "ACGTUBDHKMNRSVWXY?.~-".Contains(char.ToUpperInvariant(c), StringComparison.Ordinal));
}
