using Seqeron.Genomics.Infrastructure;

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
    /// <param name="matrix">Scoring matrix (<c>-matrixfile</c>); null = the program's default: EDNAFULL when the first sequence consists of nucleotide codes only, EBLOSUM62 otherwise (measured: the type of the first sequence decides, also when the second is of the other type). Letters are read in either case; a character that is not a letter of the matrix (<c>X</c> / <c>O</c> / <c>J</c> in a nucleotide matrix, ...) scores 0. Both sequences are read like EMBOSS reads them (<see cref="EmbossSequenceReader.Clean(string)"/>): gap characters, digits, blanks and every other non-letter except <c>*</c> are removed first and the segment positions count the cleaned sequences.</param>
    /// <param name="cancellationToken">Cancellation token, polled once per diagonal (A4-1: 100 kb x 100 kb is about 200 000 diagonals and minutes of work).</param>
    /// <returns>The segments. A window may be as long as a sequence (measured: the program then scores the one diagonal set that fits, e.g. a window equal to the length of the second sequence gives the segments of the diagonals that start in the first). A
    /// window LONGER than a sequence makes the program read beyond the end of its sequence buffer: the terminator and the zero-filled allocation slack score 0 (this port scores every position beyond the end as 0 and reports segments that end beyond the
    /// sequence: the window is counted whole), but beyond that slack - from about 8 to 24 residues past the end, depending on the allocator's chunk size of the sequence string - the program reads uninitialised memory and usually crashes: only overhangs up to
    /// 8 residues are reproducible (measured on 600 rows); a window longer than the longer sequence gives no segment.</returns>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">The first sequence is nucleic and the second is not (the program: "Sequence is not nucleic").</exception>
    public static IReadOnlyList<DotMatcherSegment> Compute(
        string sequence1, string sequence2, int windowSize = 10, int threshold = 23, SubstitutionMatrix? matrix = null, CancellationToken cancellationToken = default)
    {
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        sequence1 = EmbossSequenceReader.Clean(sequence1);
        sequence2 = EmbossSequenceReader.Clean(sequence2);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        // ACD: the type of bsequence follows asequence: a nucleic first sequence needs a nucleic second one ("Error: Sequence is not nucleic"); a protein first sequence accepts any letters in the second (measured on 37 rows)
        if (EmbossSequenceReader.IsNucleotide(sequence1) && !EmbossSequenceReader.IsNucleotide(sequence2))
            throw new ArgumentException("EMBOSS dotmatcher: the first sequence is nucleic, so the second sequence must be nucleic too (the program stops with \"Sequence is not nucleic\").", nameof(sequence2));
        matrix ??= EmbossSequenceReader.IsNucleotide(sequence1) ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;
        int len1 = sequence1.Length, len2 = sequence2.Length;
        var segments = new List<DotMatcherSegment>();

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
            AlignmentCancellation.ThrowIfRequested();
            int i = i0, k = k0, total = 0, starti = 0, startj = 0;
            bool above = false;
            // the first window: positions beyond the end of a sequence score 0 (the program reads the zero terminator / slack there), so only the part inside both sequences is summed
            for (int l = 0, inside = Math.Min(windowSize, Math.Min(len2 - i0, len1 - k0)); l < inside; l++)
                total += Score(i + l, k + l);
            i += windowSize;
            k += windowSize;
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
                        Add(new DotMatcherSegment(startj + 1, starti + 1, k, i));
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
                Add(new DotMatcherSegment(startj + 1, starti + 1, k, i));
        }

        // A4-2: the segment list is bounded (a threshold of 0 or less over an alternating window score gives about one segment per four cells)
        void Add(DotMatcherSegment segment)
        {
            if (segments.Count >= SequenceAligner.DotMatcherMaxSegments)
                throw SequenceAligner.ResourceLimit("EMBOSS dotmatcher", "keeps the segments of the plot in memory (16 bytes each)", $"more than {segments.Count} segments for {len1} x {len2} residues and a threshold of {threshold}",
                    SequenceAligner.DotMatcherMaxSegments, nameof(SequenceAligner) + "." + nameof(SequenceAligner.DotMatcherMaxSegments), "The program streams the segments to the graph file and has no such limit", "raise the threshold or the window size, or plot shorter sequences");
            segments.Add(segment);
        }
    }
}
