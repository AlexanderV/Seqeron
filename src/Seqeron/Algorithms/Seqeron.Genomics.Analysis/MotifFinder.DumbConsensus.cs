using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>Biopython <c>SummaryInfo.dumb_consensus</c> (majority threshold) and the configurable-threshold form of <c>GenerateConsensus</c>.</summary>
public static partial class MotifFinder
{

    #region Biopython dumb_consensus

    /// <summary>
    /// Majority-threshold consensus with the semantics of Biopython
    /// <c>Bio.Align.AlignInfo.SummaryInfo.dumb_consensus(threshold, ambiguous, require_multiple)</c>
    /// (Biopython ≤ 1.85; deprecated in 1.82 and removed in later releases). Per column the
    /// residues other than the gaps <c>-</c> and <c>.</c> are counted (case-sensitive, any
    /// alphabet — DNA, RNA or protein); if exactly one residue has the maximum count and that count
    /// divided by the number of non-gap residues is ≥ <paramref name="threshold"/>, it is emitted;
    /// otherwise (ties, below threshold, all-gap column, or a single non-gap residue when
    /// <paramref name="requireMultiple"/> is set) <paramref name="ambiguous"/> is emitted.
    /// </summary>
    /// <param name="alignedSequences">Aligned sequences of equal length.</param>
    /// <param name="threshold">Required fraction of the non-gap residues (Biopython default 0.7).</param>
    /// <param name="ambiguous">Symbol for columns without consensus (Biopython default <c>X</c>).</param>
    /// <param name="requireMultiple">
    /// When true, a column with exactly one non-gap residue gives <paramref name="ambiguous"/>.
    /// </param>
    /// <returns>The consensus; "" for an empty collection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alignedSequences"/> is null.</exception>
    /// <exception cref="ArgumentException">A null element or sequences of unequal length.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="threshold"/> is NaN.</exception>
    public static string GenerateDumbConsensus(
        IEnumerable<string> alignedSequences,
        double threshold = 0.7,
        char ambiguous = 'X',
        bool requireMultiple = false)
    {
        ArgumentNullException.ThrowIfNull(alignedSequences);
        if (double.IsNaN(threshold))
            throw new ArgumentOutOfRangeException(nameof(threshold), "Threshold cannot be NaN.");

        List<string> rows = MaterializeAligned(alignedSequences, nameof(alignedSequences));
        if (rows.Count == 0) return "";

        int length = rows[0].Length;
        var consensus = new StringBuilder(length);
        var counts = new Dictionary<char, int>();

        for (int n = 0; n < length; n++)
        {
            counts.Clear();
            int numAtoms = 0;
            foreach (string row in rows)
            {
                char c = row[n];
                if (c != '-' && c != '.')
                {
                    counts[c] = counts.GetValueOrDefault(c) + 1;
                    numAtoms++;
                }
            }

            char best = '\0';
            int maxSize = 0;
            int maxAtoms = 0;
            foreach (var (atom, count) in counts)
            {
                if (count > maxSize)
                {
                    best = atom;
                    maxSize = count;
                    maxAtoms = 1;
                }
                else if (count == maxSize)
                {
                    maxAtoms++;
                }
            }

            if (requireMultiple && numAtoms == 1)
                consensus.Append(ambiguous);
            else if (maxAtoms == 1 && (double)maxSize / numAtoms >= threshold)
                consensus.Append(best);
            else
                consensus.Append(ambiguous);
        }

        return consensus.ToString();
    }

    #endregion

    #region Configurable IUPAC threshold

    /// <summary>
    /// <see cref="GenerateConsensus(IEnumerable{string})"/> with a configurable per-base inclusion
    /// threshold: at each column the bases whose count is strictly greater than
    /// <paramref name="inclusionThreshold"/> × n (n = number of sequences) form the NC-IUB symbol;
    /// when no base passes, the bases tied at the maximum count are encoded (a column without
    /// A/C/G/T gives <c>N</c>). <paramref name="inclusionThreshold"/> = 0.25 is bit-identical to the
    /// parameterless overload; 0 includes every base present; values ≥ the maximum column
    /// frequency reduce every column to its (tied) most frequent bases.
    /// </summary>
    /// <param name="sequences">Aligned sequences of equal length.</param>
    /// <param name="inclusionThreshold">Per-base frequency cut, in [0, 1].</param>
    /// <returns>Consensus over the 15 IUPAC symbols; "" for an empty collection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">A null element, or sequences of unequal length.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="inclusionThreshold"/> is NaN or outside [0, 1].</exception>
    public static string GenerateConsensus(IEnumerable<string> sequences, double inclusionThreshold)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        if (double.IsNaN(inclusionThreshold) || inclusionThreshold < 0 || inclusionThreshold > 1)
            throw new ArgumentOutOfRangeException(nameof(inclusionThreshold), "Threshold must be in [0, 1].");

        return GenerateConsensusCore(sequences, inclusionThreshold);
    }

    #endregion
}
