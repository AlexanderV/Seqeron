using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>Per-sequence statistics of one row of a multiple alignment against a reference row or the consensus
/// (EMBOSS <c>infoalign</c> columns).</summary>
/// <param name="Index">0-based row index.</param>
/// <param name="SequenceLength">Residues of the row: <c>Identities + Similarities + Differences</c> (infoalign <c>SeqLen</c>).</param>
/// <param name="AlignmentLength">Columns from the row's first to its last residue (leading/trailing gap columns excluded; infoalign <c>AlignLen</c>).</param>
/// <param name="GapRuns">Internal gap runs of the row (infoalign <c>Gaps</c>).</param>
/// <param name="GapLength">Internal gap positions of the row (infoalign <c>GapLen</c>).</param>
/// <param name="Identities">Residue = reference residue (case-insensitive; infoalign <c>Ident</c>).</param>
/// <param name="Similarities">Not identical, substitution score &gt; 0 (infoalign <c>Similar</c>).</param>
/// <param name="Differences">Every other residue, including residues facing a reference gap (infoalign <c>Differ</c>).</param>
/// <param name="PercentChange">100·(AlignmentLength − Identities)/AlignmentLength (infoalign <c>% Change</c>; similarities and gap positions count as change).</param>
public sealed record MsaSequenceStatistics(
    int Index,
    int SequenceLength,
    int AlignmentLength,
    int GapRuns,
    int GapLength,
    int Identities,
    int Similarities,
    int Differences,
    double PercentChange);

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
    char Conservation);

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
    private const char Gap = '-';

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

    private static bool IsNucleotide(string[] rows, ConsensusResidueType residueType) =>
        residueType == ConsensusResidueType.Auto
            ? rows.Length > 0 && EmbossSequenceReader.IsNucleotide(rows[0])
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
    /// Pairwise distance matrix of aligned rows as EMBOSS <c>distmat</c> (uncorrected and Jukes-Cantor; no ambiguity option, all codon positions - the full option set is <see cref="CalculateDistanceMatrix(IReadOnlyList{string}, MsaDistanceOptions)"/>): matches =
    /// columns where both rows carry the same non-gap symbol; g = columns where either row has a gap;
    /// <see cref="MsaDistanceMethod.Uncorrected"/> D = 1 − matches/(L − g + g·<paramref name="gapWeight"/>);
    /// <see cref="MsaDistanceMethod.JukesCantor"/> D' = −b·ln(1 − D/b). A pair with a zero denominator is NaN
    /// (distmat divides by zero). The matrix is symmetric with a zero diagonal.
    /// </summary>
    /// <param name="alignedRows">At least two aligned rows.</param>
    /// <param name="method">Distance model.</param>
    /// <param name="residueType">Selects b = 3/4 (nucleotide) or 19/20 (protein) for Jukes-Cantor; Auto = EMBOSS type of the first row.</param>
    /// <param name="gapWeight">distmat <c>-gapweight</c> (default 0).</param>
    public static double[,] CalculateDistanceMatrix(
        IReadOnlyList<string> alignedRows, MsaDistanceMethod method = MsaDistanceMethod.Uncorrected,
        ConsensusResidueType residueType = ConsensusResidueType.Auto, double gapWeight = 0)
    {
        return CalculateDistanceMatrix(alignedRows, new MsaDistanceOptions { Method = method, ResidueType = residueType, GapWeight = gapWeight });
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

    /// <summary>
    /// EMBOSS <c>infoalign</c> statistics of every row against a reference row or against the EMBOSS <c>cons</c> consensus.
    /// </summary>
    /// <param name="alignedRows">Aligned rows (at least 2 when the consensus is the reference).</param>
    /// <param name="referenceIndex">0-based reference row, or null for the consensus (infoalign <c>-refseq 0</c>).</param>
    /// <param name="residueType">Matrix choice (EDNAFULL / EBLOSUM62) and consensus type; Auto = type of the first row as EMBOSS.</param>
    /// <remarks>
    /// infoalign_Compare (EMBOSS 6.6.0 <c>emboss/infoalign.c</c>): leading and trailing gaps of the row are ignored; inside, each gap run is
    /// one "Gap" and every '-' a gap position; each residue is Identical (case-insensitive equal to the reference), else Similar when the
    /// substitution score with the reference symbol is &gt; 0, else Different (also a reference gap, or a row longer than the reference);
    /// %Change = 100·(AlignLen − Ident)/AlignLen. A row consisting only of gaps yields zeros (EMBOSS reads past the buffer).
    /// </remarks>
    public static IReadOnlyList<MsaSequenceStatistics> CalculateSequenceStatistics(
        IReadOnlyList<string> alignedRows, int? referenceIndex = null, ConsensusResidueType residueType = ConsensusResidueType.Auto)
    {
        var rows = Normalize(alignedRows);
        bool nucleotide = IsNucleotide(rows, residueType);
        var matrix = nucleotide ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;

        string reference;
        if (referenceIndex is int r)
        {
            if (r < 0 || r >= rows.Length)
                throw new ArgumentOutOfRangeException(nameof(referenceIndex), "Reference sequence number is outside the alignment.");
            reference = rows[r];
        }
        else
        {
            // infoalign: embConsCalc(plurality 50 % of the total weight, setcase 0, identity 0) == cons defaults except setcase.
            reference = MotifFinder.GenerateEmbossConsensus(rows, residueType, plurality: null, identity: 0, setcase: 0f);
        }

        var result = new List<MsaSequenceStatistics>(rows.Length);
        for (int i = 0; i < rows.Length; i++)
            result.Add(Compare(i, reference, rows[i], matrix));
        return result;
    }

    private static MsaSequenceStatistics Compare(int index, string reference, string seq, SubstitutionMatrix matrix)
    {
        int begin = 0;
        while (begin < seq.Length && seq[begin] == Gap) begin++;
        int end = seq.Length - 1;
        while (end >= 0 && seq[end] == Gap) end--;
        if (begin > end)
            return new MsaSequenceStatistics(index, 0, 0, 0, 0, 0, 0, 0, 0);

        int gaps = 0, gapLength = 0, ident = 0, sim = 0, diff = 0;
        bool inGap = false;
        for (int i = begin; i <= end; i++)
        {
            char s = seq[i];
            if (s == Gap)
            {
                if (!inGap)
                {
                    inGap = true;
                    gaps++;
                }
                gapLength++;
                continue;
            }
            inGap = false;
            if (i >= reference.Length)
            {
                diff++;
                continue;
            }
            char r = char.ToUpperInvariant(reference[i]);
            if (r == s)
                ident++;
            else if (matrix.IsPositive(r, s))
                sim++;
            else
                diff++;
        }

        int alignLength = end - begin + 1;
        return new MsaSequenceStatistics(index, ident + sim + diff, alignLength, gaps, gapLength, ident, sim, diff,
            (float)(alignLength - ident) * 100f / alignLength);
    }

    private static string[] Normalize(IReadOnlyList<string> alignedRows)
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
