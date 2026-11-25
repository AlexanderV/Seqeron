using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;


/// <summary>Per-column statistics of a multiple alignment.</summary>
/// <param name="Index">0-based column.</param>
/// <param name="Counts">Count of each upper-cased symbol in the column, gap '-' included (Biopython <c>Alignment.frequencies</c>).</param>
/// <param name="GapCount">Rows with '-' in the column.</param>
/// <param name="GapFraction">GapCount / number of rows.</param>
/// <param name="MostCommonResidue">Most frequent non-gap symbol (smallest character on ties); '-' for an all-gap column.</param>
/// <param name="MostCommonCount">Its count.</param>
/// <param name="Conservation">ClustalW conservation symbol: '*' all rows carry the same residue (no gap), ':' all residues in one
/// "strongly similar" group, '.' all in one "weakly similar" group, ' ' otherwise (any column with a gap is ' ').</param>
public sealed record MsaColumnStatistics(
    int Index,
    IReadOnlyDictionary<char, int> Counts,
    int GapCount,
    double GapFraction,
    char MostCommonResidue,
    int MostCommonCount,
    char Conservation)
{
    /// <summary>The letters of the column (every symbol of <see cref="Counts"/> repeated by its count; the order of the rows is not kept, which the unweighted scores do not need).</summary>
    public string Column => string.Concat(Counts.Select(kv => new string(kv.Key, kv.Value)));

    /// <summary>The number of rows of the column (gaps included).</summary>
    public int RowCount => Counts.Values.Sum();

    /// <summary>
    /// Shannon entropy of the column in bits, <c>-sum p log2 p</c> over the distinct non-gap symbols with p their share of the non-gap rows (0 for a column that is all gaps or has one symbol; at most log2 of the number of different symbols).
    /// </summary>
    public double ShannonEntropy
    {
        get
        {
            int n = Counts.Where(kv => kv.Key != '-').Sum(kv => kv.Value);
            if (n == 0)
                return 0.0;
            double h = 0;
            foreach (var (symbol, count) in Counts)
            {
                if (symbol == '-' || count == 0)
                    continue;
                double p = (double)count / n;
                h -= p * Math.Log2(p);
            }

            return h == 0 ? 0.0 : h;
        }
    }

    /// <summary>
    /// Information content of the column in <paramref name="logBase"/> units as <c>Bio.Align.AlignInfo.SummaryInfo.information_content</c> counts it for a gap-ignoring column: <c>sum p log(p / e) / log(logBase)</c> over the non-gap symbols
    /// with p their share of the non-gap rows and e their expected frequency (<paramref name="expected"/>; null = uniform over the 20 amino acids, or 4 for a column of A, C, G, T / U only).
    /// </summary>
    /// <exception cref="ArgumentException">A symbol of the column has no expected frequency.</exception>
    public double InformationContent(IReadOnlyDictionary<char, double>? expected = null, double logBase = 2)
    {
        int n = Counts.Where(kv => kv.Key != '-').Sum(kv => kv.Value);
        if (n == 0)
            return 0.0;
        if (expected is null)
        {
            bool nucleotide = Counts.Keys.All(c => c is 'A' or 'C' or 'G' or 'T' or 'U' or '-');
            double uniform = nucleotide ? 0.25 : 0.05;
            expected = Counts.Keys.Where(c => c != '-').ToDictionary(c => c, _ => uniform);
        }

        double total = 0;
        foreach (var (symbol, count) in Counts)
        {
            if (symbol == '-' || count == 0)
                continue;
            if (!expected.TryGetValue(symbol, out double e))
                throw new ArgumentException($"Frequency table provided does not contain observed letter {symbol}");
            double p = (double)count / n;
            total += p * Math.Log(p / e) / Math.Log(logBase);
        }

        return total;
    }

    /// <summary>The Jensen-Shannon divergence of the column from a background distribution (Capra and Singh 2007, <see cref="ConservationScores.JsDivergence"/>; default background BLOSUM62, gap penalty on).</summary>
    public double JensenShannonDivergence(IReadOnlyList<double>? background = null, bool gapPenalty = true) => ConservationScores.JsDivergence(Column, background ?? ConservationScores.BlosumBackground, null, gapPenalty);

    /// <summary>The sum-of-pairs column score (Capra and Singh 2007, <see cref="ConservationScores.SumOfPairs"/>; default matrix BLOSUM62, gap penalty on).</summary>
    public double SumOfPairs(double[,]? similarity = null, bool gapPenalty = true) => ConservationScores.SumOfPairs(Column, similarity, null, gapPenalty);
}

/// <summary>Summary of a multiple alignment.</summary>
/// <param name="SequenceCount">Rows.</param>
/// <param name="AlignmentLength">Columns.</param>
/// <param name="Columns">Per-column statistics.</param>
/// <param name="FullyConservedColumns">Columns with Conservation '*'.</param>
/// <param name="GapFreeColumns">Columns without any gap.</param>
/// <param name="GapPercent">100·gap symbols / (rows·columns).</param>
/// <param name="PairwiseIdentity">Symmetric matrix of 100·matches/(columns where neither row has a gap) (percent identity, diagonal 100;
/// NaN when two rows share no gap-free column) — EMBOSS distmat uncorrected with gap weight 0 = 100·(1 − D).</param>
/// <param name="MeanPairwiseIdentity">Mean of the off-diagonal upper triangle of <paramref name="PairwiseIdentity"/> (NaN entries skipped; NaN when none; 100 for one row).</param>
public sealed record MsaStatistics(
    int SequenceCount,
    int AlignmentLength,
    IReadOnlyList<MsaColumnStatistics> Columns,
    int FullyConservedColumns,
    int GapFreeColumns,
    double GapPercent,
    double[,] PairwiseIdentity,
    double MeanPairwiseIdentity);

/// <summary>Distance model of <see cref="MultipleAlignmentStatistics.CalculateDistanceMatrix(IReadOnlyList{string}, MsaDistanceMethod, ConsensusResidueType, double)"/>.</summary>
public enum MsaDistanceMethod
{
    /// <summary>EMBOSS distmat <c>-nucmethod 0 / -protmethod 0</c>: D = 1 − matches/(L − g + g·gapweight) (p-distance).</summary>
    Uncorrected,

    /// <summary>EMBOSS distmat Jukes-Cantor: −b·ln(1 − D/b), b = 3/4 (nucleotides) or 19/20 (proteins).</summary>
    JukesCantor,

    /// <summary>EMBOSS distmat <c>-nucmethod 2</c> Kimura 2-parameter (M Kimura 1980, J Mol Evol 16:111) for nucleotides, <c>-protmethod 2</c> "Kimura Protein" (D' = -ln(1 - D - 0.2 D²)) for proteins.</summary>
    Kimura,

    /// <summary>EMBOSS distmat <c>-nucmethod 3</c> Tamura (K Tamura 1992, Mol Biol Evol 9:678); nucleotides only.</summary>
    Tamura,

    /// <summary>EMBOSS distmat <c>-nucmethod 4</c> Tajima-Nei (Tajima &amp; Nei 1984, Mol Biol Evol 1:269); nucleotides only.</summary>
    TajimaNei,

    /// <summary>EMBOSS distmat <c>-nucmethod 5</c> Jin-Nei Gamma distance (Jin &amp; Nei 1990, Mol Biol Evol 7:82); nucleotides only.</summary>
    JinNei,
}

/// <summary>
/// Statistics of a multiple sequence alignment: per-sequence comparison with a reference or consensus
/// (EMBOSS <c>infoalign</c>), pairwise distance / identity matrix (EMBOSS <c>distmat</c>; Biopython
/// <c>DistanceCalculator('identity')</c> as cross-check), per-column gap / count / ClustalW conservation statistics
/// (Biopython <c>Alignment.frequencies</c>, ClustalW conservation line).
/// </summary>
/// <remarks>
/// Rows are upper-cased, '.' and '~' read as the gap '-' (EMBOSS reader). Rows of unequal length are rejected. The consensus
/// reference of <c>infoalign</c> is the EMBOSS <c>cons</c> consensus (plurality 50 % of the total weight, identity 0), i.e.
/// <see cref="MotifFinder.GenerateEmbossConsensus(IEnumerable{string}, ConsensusResidueType, float?, int, float?, IReadOnlyList{float}?)"/>
/// — called, not re-implemented.
/// </remarks>
public static partial class MultipleAlignmentStatistics
{
    public const char Gap = '-';

    // ClustalW conservation groups (ClustalW 2.1, clustalw.h: strongGroups / weakGroups).
    private static readonly string[] StrongGroups = { "STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW" };
    private static readonly string[] WeakGroups = { "CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY" };

    /// <summary>Full statistics of a multiple alignment result.</summary>
    public static MsaStatistics Calculate(MultipleAlignmentResult alignment, ConsensusResidueType residueType = ConsensusResidueType.Auto)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        return Calculate(alignment.AlignedSequences, residueType);
    }

    /// <summary>
    /// Per-column statistics, ClustalW conservation line, gap statistics and the percent-identity matrix of aligned rows.
    /// </summary>
    /// <param name="alignedRows">At least one aligned row (equal lengths).</param>
    /// <param name="residueType">Nucleotide: the conservation line only marks identical columns ('*'), as ClustalW does for DNA;
    /// Protein: also the strong ':' / weak '.' similarity groups; Auto: nucleotide when the first row consists of nucleotide
    /// (IUPAC) symbols and gaps only (the EMBOSS rule of <c>cons</c>/<c>infoalign</c>; ClustalW itself needs <c>-type</c> or uses its own guess).</param>
    public static MsaStatistics Calculate(IReadOnlyList<string> alignedRows, ConsensusResidueType residueType = ConsensusResidueType.Auto)
    {
        var rows = Normalize(alignedRows);
        bool nucleotide = IsNucleotide(rows, residueType);
        int n = rows.Length;
        int len = n == 0 ? 0 : rows[0].Length;

        var columns = new List<MsaColumnStatistics>(len);
        int conserved = 0, gapFree = 0;
        long gapSymbols = 0;
        for (int k = 0; k < len; k++)
        {
            var counts = new SortedDictionary<char, int>();
            for (int i = 0; i < n; i++)
            {
                counts.TryGetValue(rows[i][k], out int c);
                counts[rows[i][k]] = c + 1;
            }
            counts.TryGetValue(Gap, out int gaps);
            char best = Gap;
            int bestCount = 0;
            foreach (var (ch, c) in counts)
            {
                if (ch != Gap && c > bestCount)
                {
                    best = ch;
                    bestCount = c;
                }
            }

            char conservation = ClustalConservationSymbol(counts.Keys, gaps > 0, nucleotide);
            if (conservation == '*') conserved++;
            if (gaps == 0) gapFree++;
            gapSymbols += gaps;
            columns.Add(new MsaColumnStatistics(k, counts, gaps, n == 0 ? 0 : (double)gaps / n, best, bestCount, conservation));
        }

        var identity = PairwiseIdentity(rows);
        double sum = 0;
        int cnt = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (double.IsNaN(identity[i, j])) continue;
                sum += identity[i, j];
                cnt++;
            }
        }
        double mean = cnt == 0 ? double.NaN : sum / cnt;
        if (n == 1)
            mean = 100.0;

        return new MsaStatistics(n, len, columns, conserved, gapFree,
            n * len == 0 ? 0.0 : 100.0 * gapSymbols / ((double)n * len), identity, mean);
    }

    /// <summary>The ClustalW conservation line ('*', ':', '.', ' ' per column) of aligned rows.</summary>
    public static string ClustalConservationLine(IReadOnlyList<string> alignedRows, ConsensusResidueType residueType = ConsensusResidueType.Auto)
    {
        var rows = Normalize(alignedRows);
        bool nucleotide = IsNucleotide(rows, residueType);
        int len = rows.Length == 0 ? 0 : rows[0].Length;
        var sb = new StringBuilder(len);
        for (int k = 0; k < len; k++)
        {
            var symbols = new SortedSet<char>();
            bool gap = false;
            foreach (var r in rows)
            {
                symbols.Add(r[k]);
                gap |= r[k] == Gap;
            }
            sb.Append(ClustalConservationSymbol(symbols, gap, nucleotide));
        }
        return sb.ToString();
    }

    // the residue type "decided from the first row": every character is an IUPAC nucleotide code (or ? / a gap . ~ -)
    private static bool IsNucleotideRow(string row)
    {
        foreach (char ch in row)
        {
            switch (char.ToUpperInvariant(ch))
            {
                case 'A' or 'C' or 'G' or 'T' or 'U' or 'B' or 'D' or 'H' or 'K' or 'M' or 'N' or 'R' or 'S' or 'V' or 'W' or 'X' or 'Y' or '?' or '.' or '~' or '-':
                    break;
                default:
                    return false;
            }
        }

        return true;
    }

    public static bool IsNucleotide(string[] rows, ConsensusResidueType residueType) =>
        residueType == ConsensusResidueType.Auto
            ? rows.Length > 0 && IsNucleotideRow(rows[0])
            : residueType == ConsensusResidueType.Nucleotide;

    private static char ClustalConservationSymbol(IEnumerable<char> symbols, bool hasGap, bool nucleotide)
    {
        if (hasGap)
            return ' ';
        var set = new HashSet<char>(symbols);
        if (set.Count == 1)
            return '*';
        if (nucleotide)
            return ' ';
        if (StrongGroups.Any(g => set.All(g.Contains)))
            return ':';
        if (WeakGroups.Any(g => set.All(g.Contains)))
            return '.';
        return ' ';
    }

    /// <summary>
    /// Percent-identity matrix 100·matches/(L − g + g·<paramref name="gapWeight"/>) = 100·(1 − distmat uncorrected D) with L columns and
    /// g columns where either row has a gap (diagonal 100; NaN for a zero denominator).
    /// <paramref name="gapWeight"/> = 0 (default): identity over the columns where neither row has a gap (EMBOSS distmat default);
    /// 1: identity over the whole alignment length, gap columns counting as differences (as a needle-style identity). Biopython's
    /// <c>DistanceCalculator('identity')</c> (<c>score / len(seq)</c>, <c>skip_letters = ()</c>) equals this at gap weight 1 plus
    /// the columns where both rows have '-' (it compares characters, so '-' = '-' scores as an identity).
    /// </summary>
    public static double[,] CalculatePairwiseIdentityMatrix(IReadOnlyList<string> alignedRows, double gapWeight = 0)
    {
        if (double.IsNaN(gapWeight) || double.IsInfinity(gapWeight))
            throw new ArgumentOutOfRangeException(nameof(gapWeight));
        return PairwiseIdentity(Normalize(alignedRows), gapWeight);
    }

    private static double[,] PairwiseIdentity(string[] rows, double gapWeight = 0)
    {
        int n = rows.Length;
        var identity = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            identity[i, i] = 100.0;
            for (int j = i + 1; j < n; j++)
            {
                int matches = 0, gaps = 0;
                for (int k = 0; k < rows[i].Length; k++)
                {
                    char x = rows[i][k], y = rows[j][k];
                    if (x == Gap || y == Gap)
                    {
                        gaps++;
                        continue;
                    }
                    if (x == y) matches++;
                }
                double denominator = rows[i].Length - gaps + gaps * gapWeight;
                double v = denominator == 0 ? double.NaN : 100.0 * matches / denominator;
                identity[i, j] = v;
                identity[j, i] = v;
            }
        }
        return identity;
    }

    public static string[] Normalize(IReadOnlyList<string> alignedRows)
    {
        ArgumentNullException.ThrowIfNull(alignedRows);
        var rows = new string[alignedRows.Count];
        for (int i = 0; i < rows.Length; i++)
        {
            if (alignedRows[i] is null)
                throw new ArgumentException($"Row {i} is null.", nameof(alignedRows));
            rows[i] = alignedRows[i].ToUpperInvariant().Replace('.', Gap).Replace('~', Gap);
        }
        if (rows.Length > 0 && rows.Any(r => r.Length != rows[0].Length))
            throw new ArgumentException("Aligned rows must have equal length.", nameof(alignedRows));
        return rows;
    }
}
