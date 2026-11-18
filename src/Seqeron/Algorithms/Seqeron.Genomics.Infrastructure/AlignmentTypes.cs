namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// Scalar scoring parameters for sequence alignment (one score for every identical pair, one for every
/// non-identical pair, and two gap parameters).
/// </summary>
/// <remarks>
/// <para>
/// Gap convention. The affine aligners (<c>GlobalAlignAffine</c>, <c>LocalAlignAffine</c>,
/// <c>SemiGlobalAlignAffine</c>) charge a gap of length k as <see cref="GapOpen"/> + k·<see cref="GapExtend"/>
/// (Gotoh 1982 w(k) = v + u·k; the NCBI BLAST+ convention, "existence 5, extension 2" = 5 + 2k; Biopython
/// <c>open_gap_score</c> = <see cref="GapOpen"/> + <see cref="GapExtend"/>). The linear aligners
/// (<c>GlobalAlign</c>, <c>LocalAlign</c>, <c>SemiGlobalAlign</c>, the multiple-alignment methods) use only
/// <see cref="GapExtend"/> per gap position. EMBOSS needle/water instead charge gapopen + (k−1)·gapextend for
/// the same two numbers, so EMBOSS "-gapopen 10 -gapextend 0.5" corresponds to GapOpen = −9.5 here (not
/// representable with integers).
/// </para>
/// <para>
/// The presets live on <c>SequenceAligner</c> (<c>SimpleDna</c>, <c>BlastDna</c>, <c>HighIdentityDna</c>). For
/// residue-pair specific scores (IUPAC DNA, protein) use <see cref="SubstitutionMatrix"/>.
/// </para>
/// </remarks>
public sealed record ScoringMatrix(
    int Match,
    int Mismatch,
    int GapOpen,
    int GapExtend);

/// <summary>
/// Residue-pair scoring for pairwise alignment: a <see cref="SubstitutionMatrix"/> (BLOSUM62, PAM250, EDNAFULL, ...)
/// plus <b>real-valued</b> gap scores. Accepted by the <c>GlobalAlign</c>/<c>LocalAlign</c>/<c>SemiGlobalAlign</c>
/// (and <c>*Affine</c>) overloads of <c>SequenceAligner</c>; protein, IUPAC DNA and any alphabet of the matrix work.
/// </summary>
/// <remarks>
/// <para>
/// Gap convention = <see cref="ScoringMatrix"/>: a gap is a <b>score</b> (negative = penalty). The affine aligners
/// charge a gap of length k as <see cref="GapOpen"/> + k·<see cref="GapExtend"/>; the linear aligners charge
/// <see cref="GapExtend"/> per gap position (<see cref="GapOpen"/> unused). EMBOSS needle/water and Biopython
/// differ only in how the same numbers are named: EMBOSS <c>-gapopen G -gapextend E</c> (positive penalties; a gap of
/// length k costs G + (k-1)·E) is <c>GapOpen = -(G - E)</c>, <c>GapExtend = -E</c> — use
/// <see cref="Emboss(SubstitutionMatrix, double, double)"/>. Biopython <c>open_gap_score</c> = <c>GapOpen + GapExtend</c>,
/// <c>extend_gap_score</c> = <c>GapExtend</c>.
/// </para>
/// <para>
/// Gap scores are resolved exactly to 1e-6 (more decimal places are rejected); internally every score is an
/// exact multiple of 1/scale (scale = smallest integer making the gaps integral: 2 for the EMBOSS default 0.5, 10 for 0.1),
/// so ties and tracebacks are exact (no floating-point drift) and the integer dynamic programs are reused unchanged.
/// </para>
/// </remarks>
/// <param name="Matrix">Residue-pair scores; every sequence symbol must be in its alphabet (case-insensitive).</param>
/// <param name="GapOpen">Score added once per gap run (affine aligners only); at most 0.</param>
/// <param name="GapExtend">Score per gap position; at most 0.</param>
public sealed record MatrixScoring(SubstitutionMatrix Matrix, double GapOpen, double GapExtend)
{
    /// <summary>
    /// Scoring from EMBOSS needle/water penalties (positive numbers): <paramref name="gapOpenPenalty"/> is the cost of a gap
    /// of length 1, every further position costs <paramref name="gapExtendPenalty"/> (needle/water defaults: 10 and 0.5).
    /// </summary>
    public static MatrixScoring Emboss(SubstitutionMatrix matrix, double gapOpenPenalty = 10, double gapExtendPenalty = 0.5) =>
        new(matrix, -(gapOpenPenalty - gapExtendPenalty), -gapExtendPenalty);
}

/// <summary>
/// Type of alignment performed.
/// </summary>
public enum AlignmentType
{
    Global,
    Local,
    SemiGlobal
}

/// <summary>
/// End-gap policy of an <b>ends-free</b> (semi-global) alignment: which of the four alignment ends may overhang
/// without penalty. "Start"/"End" are the left/right ends of the alignment, "1"/"2" the sequence whose residues
/// overhang: <see cref="FreeStart1"/> = the leading residues of sequence 1 that precede the alignment are free (they are
/// shown against '-' in sequence 2), and so on. Flags not set are charged like any other gap
/// (linear d per position; affine GapOpen + k·GapExtend). No flag = global (Needleman-Wunsch) alignment.
/// </summary>
/// <remarks>
/// Equivalent settings in the reference tools (sequence 1 = Biopython target = parasail query/s1, sequence 2 = Biopython
/// query = parasail database/s2; a free overhang of sequence 1 is a gap in sequence 2, i.e. a Biopython <i>deletion</i>):
/// <see cref="FreeStart1"/> = Biopython <c>left_deletion_score = 0</c> = parasail <c>_qb</c>; <see cref="FreeEnd1"/> =
/// <c>right_deletion_score = 0</c> = <c>_qe</c>; <see cref="FreeStart2"/> = <c>left_insertion_score = 0</c> = <c>_db</c>;
/// <see cref="FreeEnd2"/> = <c>right_insertion_score = 0</c> = <c>_de</c>; <see cref="Fitting"/> = <c>sg_dx</c>
/// (<c>end_insertion_score = 0</c>); <see cref="All"/> = <c>sg</c> = EMBOSS needle's default (<c>-endweight N</c>).
/// </remarks>
[Flags]
public enum EndGaps
{
    /// <summary>All end gaps are penalised: global (Needleman-Wunsch) alignment.</summary>
    None = 0,

    /// <summary>Leading residues of sequence 1 before the alignment are free (Biopython <c>left_deletion_score = 0</c>, parasail <c>_qb</c>).</summary>
    FreeStart1 = 1,

    /// <summary>Trailing residues of sequence 1 after the alignment are free (Biopython <c>right_deletion_score = 0</c>, parasail <c>_qe</c>).</summary>
    FreeEnd1 = 2,

    /// <summary>Leading residues of sequence 2 before the alignment are free (Biopython <c>left_insertion_score = 0</c>, parasail <c>_db</c>).</summary>
    FreeStart2 = 4,

    /// <summary>Trailing residues of sequence 2 after the alignment are free (Biopython <c>right_insertion_score = 0</c>, parasail <c>_de</c>).</summary>
    FreeEnd2 = 8,

    /// <summary>Both flanks of sequence 2 are free: sequence 1 is fitted end to end into sequence 2 ("fitting"/"glocal"; the classic
    /// <c>SemiGlobalAlign</c>; Biopython <c>end_insertion_score = 0</c>, parasail <c>sg_dx</c>, Rosalind SIMS).</summary>
    Fitting = FreeStart2 | FreeEnd2,

    /// <summary>Both flanks of sequence 1 are free: sequence 2 is fitted end to end into sequence 1 (parasail <c>sg_qx</c>).</summary>
    FittingSequence2 = FreeStart1 | FreeEnd1,

    /// <summary>Overlap (suffix-prefix): a suffix of sequence 1 overlaps a prefix of sequence 2 — the leading part of sequence 1
    /// and the trailing part of sequence 2 are free (parasail <c>sg_qb_de</c>; sequence assembly "dovetail" overlap).</summary>
    OverlapSuffix1Prefix2 = FreeStart1 | FreeEnd2,

    /// <summary>Overlap (suffix-prefix), other orientation: a suffix of sequence 2 overlaps a prefix of sequence 1
    /// (parasail <c>sg_qe_db</c>).</summary>
    OverlapSuffix2Prefix1 = FreeEnd1 | FreeStart2,

    /// <summary>All four end gaps are free (parasail <c>sg</c>, Rosalind SMGB-like, EMBOSS needle's default <c>-endweight N</c>).</summary>
    All = FreeStart1 | FreeEnd1 | FreeStart2 | FreeEnd2,
}

/// <summary>
/// Result of pairwise sequence alignment.
/// </summary>
/// <param name="AlignedSequence1">Aligned (gapped, '-') segment of sequence 1.</param>
/// <param name="AlignedSequence2">Aligned (gapped, '-') segment of sequence 2.</param>
/// <param name="Score">Alignment score under the scoring model used.</param>
/// <param name="AlignmentType">Kind of alignment that produced the result.</param>
/// <param name="StartPosition1">0-based inclusive start of the aligned segment in sequence 1
/// (0 for global and for semi-global results in which sequence 1 is aligned in full; for an ends-free result
/// (<see cref="EndGaps"/>) the first index of sequence 1 consumed by the alignment path, i.e. the number of free leading
/// sequence-1 residues; −1 for a local result with no positive-scoring region).</param>
/// <param name="StartPosition2">0-based inclusive start of the aligned segment in sequence 2: 0 for global; for a
/// <b>semi-global</b> result the first index of sequence 2 that the alignment consumes, i.e. the number of
/// free leading sequence-2 residues (fitting: Biopython <c>alignment.coordinates</c>/<c>aligned</c> start of the fitted
/// interval); −1 for a local result with no positive-scoring region.</param>
/// <param name="EndPosition1">0-based inclusive end of the aligned segment in sequence 1
/// (length − 1 for global and fitting, so −1 for an empty sequence; for an ends-free result the last index of sequence 1
/// consumed by the path, the free trailing sequence-1 residues lie after it; −1 for an empty local result).</param>
/// <param name="EndPosition2">0-based inclusive end of the aligned segment in sequence 2: length − 1 for global; for a
/// <b>semi-global</b> result the last index of sequence 2 that the alignment consumes (the free trailing
/// residues lie after it; parasail <c>sg_dx</c> <c>end_ref</c>). A consumed interval that is empty
/// (a sequence aligned entirely against gaps, or empty) has <c>End = Start − 1</c>.
/// −1 for an empty local result.</param>
/// <remarks>
/// Empty inputs: a global alignment of an empty side is the all-gap alignment (Needleman-Wunsch optimum, Biopython
/// <c>PairwiseAligner.score</c>); two empty inputs give score 0, empty rows and coordinates 0/0/−1/−1. A local
/// alignment with no positive-scoring pair is <c>("", "", 0, Local, −1, −1, −1, −1)</c>.
/// </remarks>
public sealed record AlignmentResult(
    string AlignedSequence1,
    string AlignedSequence2,
    int Score,
    AlignmentType AlignmentType,
    int StartPosition1,
    int StartPosition2,
    int EndPosition1,
    int EndPosition2)
{
    private readonly double? _exactScore;

    /// <summary>
    /// Alignment score as a real number. Equal to <see cref="Score"/> for every integer-scored aligner; the
    /// <see cref="MatrixScoring"/> aligners (substitution matrix, possibly fractional gap penalties such as
    /// EMBOSS <c>-gapextend 0.5</c>) set it to the exact score (e.g. 49.5) while <see cref="Score"/> then holds
    /// that value rounded half away from zero.
    /// </summary>
    public double ExactScore
    {
        get => _exactScore ?? Score;
        init => _exactScore = value;
    }

    public static AlignmentResult Empty => new("", "", 0, AlignmentType.Global, 0, 0, 0, 0);
}

/// <summary>
/// Statistics calculated from a pairwise alignment (EMBOSS needle/water + BLAST conventions; see
/// <c>SequenceAligner.CalculateStatistics</c>).
/// </summary>
/// <param name="Matches">Identical aligned residue pairs (EMBOSS/BLAST "identities"; Biopython <c>identities</c>).</param>
/// <param name="Mismatches">Aligned residue pairs that are not identical (Biopython <c>mismatches</c>).</param>
/// <param name="Gaps">Columns in which at least one row has a gap (EMBOSS "Gaps" count; equals the total gap length for a
/// valid pairwise alignment, Biopython <c>gaps</c>).</param>
/// <param name="AlignmentLength">Number of alignment columns, gap columns included (EMBOSS "Length"); the denominator of all percentages.</param>
/// <param name="Identity">100·Matches / AlignmentLength.</param>
/// <param name="Similarity">100·Positives / AlignmentLength (EMBOSS "Similarity", BLAST "Positives").</param>
/// <param name="GapPercent">100·Gaps / AlignmentLength.</param>
/// <param name="Positives">Aligned residue pairs (identical or not) whose substitution score is &gt; 0 (EMBOSS similarity count,
/// BLAST positives, Biopython <c>positives</c>).</param>
/// <param name="GapOpens">Number of gap runs: maximal runs of '-' in either row (BLAST tabular "gapopen", Biopython
/// <c>open_gaps</c>). A gap run in one row directly followed by a gap run in the other row counts as two.</param>
/// <param name="BlastPositives">BLAST "Positives" count: identical columns plus non-identical columns with a positive
/// substitution score (BLAST+ <c>x_FillIdentityInfo</c>: match + '+'). Differs from <paramref name="Positives"/> only for identical
/// pairs whose own score is not positive (X/X in BLOSUM62), which BLAST counts as identities.</param>
public readonly record struct AlignmentStatistics(
    int Matches,
    int Mismatches,
    int Gaps,
    int AlignmentLength,
    double Identity,
    double Similarity,
    double GapPercent,
    int Positives = 0,
    int GapOpens = 0,
    int BlastPositives = 0)
{
    public static AlignmentStatistics Empty => new(0, 0, 0, 0, 0, 0, 0);

    /// <summary>
    /// BLAST pairwise-report percentage (NCBI BLAST+ <c>CAlignFormatUtil::GetPercentMatch</c>, used for the
    /// <c>Identities = a/L (p%)</c>, <c>Positives</c> and <c>Gaps</c> lines): an <b>integer</b>, rounded half up
    /// (<c>(int)(0.5 + 100·n/d)</c>) but never 100 unless n = d (an alignment with mismatches is shown as at most 99%).
    /// A zero-length alignment gives 0.
    /// </summary>
    public static int BlastPercent(int numerator, int denominator)
    {
        if (denominator <= 0)
            return 0;
        if (numerator == denominator)
            return 100;
        return Math.Min(99, (int)(0.5 + 100.0 * numerator / denominator));
    }

    /// <summary>BLAST "Identities" percentage, <see cref="BlastPercent"/>(<see cref="Matches"/>, <see cref="AlignmentLength"/>).</summary>
    public int BlastPercentIdentity => BlastPercent(Matches, AlignmentLength);

    /// <summary>BLAST "Positives" percentage, <see cref="BlastPercent"/>(<see cref="BlastPositives"/>, <see cref="AlignmentLength"/>).</summary>
    public int BlastPercentPositives => BlastPercent(BlastPositives, AlignmentLength);

    /// <summary>BLAST "Gaps" percentage, <see cref="BlastPercent"/>(<see cref="Gaps"/>, <see cref="AlignmentLength"/>).</summary>
    public int BlastPercentGaps => BlastPercent(Gaps, AlignmentLength);

    /// <summary>
    /// BLAST tabular (<c>-outfmt 6</c>) <c>pident</c> text: 100·Matches/AlignmentLength with three decimals (BLAST+ 2.12;
    /// "100.000" when every column matches), invariant culture.
    /// </summary>
    public string BlastTabularPidentText
    {
        get
        {
            double pident = 0.0;
            if (AlignmentLength > 0)
                pident = Matches == AlignmentLength ? 100.0 : 100.0 * Matches / AlignmentLength;
            return pident.ToString("F3", System.Globalization.CultureInfo.InvariantCulture);
        }
    }

    /// <summary>
    /// The BLAST pairwise-report statistics line: <c> Identities = 382/403 (95%), Gaps = 9/403 (2%)</c> for nucleotides, with
    /// <c>, Positives = n/L (p%)</c> after the identities for proteins (<paramref name="protein"/> = true).
    /// </summary>
    public string ToBlastString(bool protein = false)
    {
        var sb = new System.Text.StringBuilder();
        sb.Append(" Identities = ").Append(Matches).Append('/').Append(AlignmentLength).Append(" (").Append(BlastPercentIdentity).Append("%)");
        if (protein)
            sb.Append(", Positives = ").Append(BlastPositives).Append('/').Append(AlignmentLength).Append(" (").Append(BlastPercentPositives).Append("%)");
        sb.Append(", Gaps = ").Append(Gaps).Append('/').Append(AlignmentLength).Append(" (").Append(BlastPercentGaps).Append("%)");
        return sb.ToString();
    }
}

/// <summary>
/// Residue-pair substitution matrix (symbol × symbol → integer score), e.g. BLOSUM62 or the IUPAC nucleotide
/// matrix EDNAFULL. Used to decide which aligned pairs are <i>similar</i> (positive score) in
/// <c>SequenceAligner.CalculateStatistics</c> and <c>FormatAlignment</c>.
/// </summary>
/// <remarks>
/// Lookups are case-insensitive. The built-in tables are the NCBI/EMBOSS matrices
/// (<c>EBLOSUM62</c>, <c>EPAM250</c>, <c>EDNAFULL</c> from the EMBOSS 6.6 data directory; BLOSUM62 = Henikoff &amp;
/// Henikoff 1992 PNAS 89:10915, PAM250 = Dayhoff, Schwartz &amp; Orcutt 1978, EDNAFULL = Todd Lowe 1992 = NCBI NUC.4.4 +
/// U). Every cell equals Biopython 1.88 <c>substitution_matrices.load('BLOSUM62' / 'PAM250' / 'NUC.4.4')</c>
/// (EDNAFULL additionally defines U).
/// </remarks>
public sealed class SubstitutionMatrix
{
    private readonly int[,] _scores;
    private readonly sbyte[] _index = new sbyte[128];

    /// <summary>Matrix name (e.g. "BLOSUM62").</summary>
    public string Name { get; }

    /// <summary>Upper-case symbols in row/column order.</summary>
    public string Alphabet { get; }

    /// <summary>Creates a matrix from a square score table in <paramref name="alphabet"/> order.</summary>
    /// <exception cref="ArgumentException">Duplicate/non-ASCII/lower-case symbols or a non-square table.</exception>
    public SubstitutionMatrix(string name, string alphabet, int[,] scores)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        ArgumentException.ThrowIfNullOrEmpty(alphabet);
        ArgumentNullException.ThrowIfNull(scores);
        if (alphabet.Length > sbyte.MaxValue || scores.GetLength(0) != alphabet.Length || scores.GetLength(1) != alphabet.Length)
            throw new ArgumentException("Score table must be square with one row/column per alphabet symbol.", nameof(scores));

        Array.Fill(_index, (sbyte)-1);
        for (int i = 0; i < alphabet.Length; i++)
        {
            char c = alphabet[i];
            if (c >= 128 || c != char.ToUpperInvariant(c) || _index[c] >= 0)
                throw new ArgumentException($"Alphabet symbol '{c}' is duplicated, lower-case or not ASCII.", nameof(alphabet));
            _index[c] = (sbyte)i;
        }

        Name = name;
        Alphabet = alphabet;
        _scores = (int[,])scores.Clone();
    }

    /// <summary>Looks up the score of an aligned pair; false when either symbol is not in <see cref="Alphabet"/>.</summary>
    public bool TryGetScore(char a, char b, out int score)
    {
        score = 0;
        a = char.ToUpperInvariant(a);
        b = char.ToUpperInvariant(b);
        if (a >= 128 || b >= 128 || _index[a] < 0 || _index[b] < 0)
            return false;
        score = _scores[_index[a], _index[b]];
        return true;
    }

    /// <summary>True when <paramref name="symbol"/> (case-insensitive) is in <see cref="Alphabet"/>.</summary>
    public bool Contains(char symbol)
    {
        symbol = char.ToUpperInvariant(symbol);
        return symbol < 128 && _index[symbol] >= 0;
    }

    /// <summary>Largest absolute cell value.</summary>
    public int MaxAbsScore
    {
        get
        {
            int max = 0;
            foreach (int v in _scores)
                max = Math.Max(max, Math.Abs(v));
            return max;
        }
    }

    /// <summary>Row/column index of an upper-case symbol, -1 when absent (alignment dynamic programs).</summary>
    public int IndexOf(char upperSymbol) => upperSymbol < 128 ? _index[upperSymbol] : -1;

    /// <summary>Raw cell by row/column index (see <see cref="IndexOf"/>).</summary>
    public int ScoreAt(int row, int column) => _scores[row, column];

    /// <summary>Score of an aligned pair.</summary>
    /// <exception cref="ArgumentException">A symbol is not in <see cref="Alphabet"/>.</exception>
    public int GetScore(char a, char b) =>
        TryGetScore(a, b, out int s) ? s : throw new ArgumentException($"Symbol pair '{a}','{b}' is not in matrix {Name} ({Alphabet}).");

    /// <summary>BLOSUM62 (EMBOSS <c>EBLOSUM62</c> = Biopython <c>BLOSUM62</c>).</summary>
    public static SubstitutionMatrix Blosum62 => LazyBlosum62.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyBlosum62 = new(() => Build("BLOSUM62", "ARNDCQEGHILKMFPSTWYVBZX*", Blosum62Rows()));

    /// <summary>PAM250 (EMBOSS <c>EPAM250</c> = Biopython <c>PAM250</c>).</summary>
    public static SubstitutionMatrix Pam250 => LazyPam250.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyPam250 = new(() => Build("PAM250", "ARNDCQEGHILKMFPSTWYVBZX*", Pam250Rows()));

    /// <summary>IUPAC nucleotide matrix EDNAFULL (= NCBI NUC.4.4, match 5 / mismatch −4, plus U).</summary>
    public static SubstitutionMatrix EdnaFull => LazyEdnaFull.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyEdnaFull = new(() => Build("EDNAFULL", "ATGCSWRYKMBVHDNU", EdnaFullRows()));

    /// <summary>
    /// ClustalW 2.1's GONNET250 matrix (<c>gon250mt</c> of <c>matrices.h</c>; the default PROTEIN pairwise matrix <c>-PWMATRIX=GONNET</c>): Gonnet, Cohen &amp;
    /// Benner 1992 (Science 256:1443) log-odds at 250 PAM in tenth-bits, i.e. every cell equals <c>round(10 × Biopython GONNET1992)</c> on the 20 amino acids.
    /// Alphabet <c>ABCDEFGHIKLMNPQRSTVWXYZ</c> (ClustalW's packing order; B = Asx, Z = Glx, X = unknown).
    /// </summary>
    public static SubstitutionMatrix ClustalWGonnet250 => LazyClustalWGonnet250.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWGonnet250 = new(() => BuildTriangle("CLUSTALW-GONNET250", ClustalWGon250Triangle()));

    /// <summary>
    /// ClustalW 2.1's BLOSUM30 matrix (<c>blosum30mt</c>; the PROTEIN pairwise matrix <c>-PWMATRIX=BLOSUM</c>): Henikoff &amp; Henikoff 1992, equal on all 23 × 23 cells to
    /// the EMBOSS data file <c>EBLOSUM30</c>.
    /// </summary>
    public static SubstitutionMatrix ClustalWBlosum30 => LazyClustalWBlosum30.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWBlosum30 = new(() => BuildTriangle("CLUSTALW-BLOSUM30", ClustalWBlosum30Triangle()));

    /// <summary>
    /// ClustalW 2.1's PAM350 matrix (<c>pam350mt</c>; the PROTEIN pairwise matrix <c>-PWMATRIX=PAM</c>): the Dayhoff mutation-data log-odds at 350 PAM as distributed
    /// with ClustalW (the table of the original program; no other local reference carries PAM350 — its own PAM250 sibling <c>pam250mt</c> equals Biopython <c>PAM250</c> on the 20 amino acids).
    /// </summary>
    public static SubstitutionMatrix ClustalWPam350 => LazyClustalWPam350.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWPam350 = new(() => BuildTriangle("CLUSTALW-PAM350", ClustalWPam350Triangle()));

    /// <summary>ClustalW 2.1's identity matrix (<c>idmat</c>; <c>-PWMATRIX=ID</c>): 10 for every identical pair of residues (incl. B, X, Z), 0 otherwise.</summary>
    public static SubstitutionMatrix ClustalWIdentity => LazyClustalWIdentity.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyClustalWIdentity = new(() => BuildTriangle("CLUSTALW-ID", ClustalWIdTriangle()));

    /// <summary>ClustalW's amino-acid packing order of the lower-triangular tables of <c>matrices.h</c> (23 residues).</summary>
    public const string ClustalWAminoAcidOrder = "ABCDEFGHIKLMNPQRSTVWXYZ";

    private static SubstitutionMatrix BuildTriangle(string name, short[] triangle)
    {
        const string alphabet = ClustalWAminoAcidOrder;
        int n = alphabet.Length;
        var m = new int[n, n];
        int ix = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                m[i, j] = triangle[ix];
                m[j, i] = triangle[ix];
                ix++;
            }

        return new SubstitutionMatrix(name, alphabet, m);
    }

    private static SubstitutionMatrix Build(string name, string alphabet, string[] rows)
    {
        int n = alphabet.Length;
        var m = new int[n, n];
        for (int i = 0; i < n; i++)
        {
            var parts = rows[i].Split(' ');
            for (int j = 0; j < n; j++)
                m[i, j] = int.Parse(parts[j], System.Globalization.CultureInfo.InvariantCulture);
        }
        return new SubstitutionMatrix(name, alphabet, m);
    }

    private static string[] Blosum62Rows() => new[]
    {
        "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0 -2 -1 0 -4",
        "-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3 -1 0 -1 -4",
        "-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3 3 0 -1 -4",
        "-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3 4 1 -1 -4",
        "0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4",
        "-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2 0 3 -1 -4",
        "-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2 1 4 -1 -4",
        "0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3 -1 -2 -1 -4",
        "-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3 0 0 -1 -4",
        "-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3 -3 -3 -1 -4",
        "-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1 -4 -3 -1 -4",
        "-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2 0 1 -1 -4",
        "-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1 -3 -1 -1 -4",
        "-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1 -3 -3 -1 -4",
        "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2 -2 -1 -2 -4",
        "1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2 0 0 0 -4",
        "0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0 -1 -1 0 -4",
        "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3 -4 -3 -2 -4",
        "-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1 -3 -2 -1 -4",
        "0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4 -3 -2 -1 -4",
        "-2 -1 3 4 -3 0 1 -1 0 -3 -4 0 -3 -3 -2 0 -1 -4 -3 -3 4 1 -1 -4",
        "-1 0 0 1 -3 3 4 -2 0 -3 -3 1 -1 -3 -1 0 -1 -3 -2 -2 1 4 -1 -4",
        "0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 0 0 -2 -1 -1 -1 -1 -1 -4",
        "-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 1"
    };

    private static string[] Pam250Rows() => new[]
    {
        "2 -2 0 0 -2 0 0 1 -1 -1 -2 -1 -1 -3 1 1 1 -6 -3 0 0 0 0 -8",
        "-2 6 0 -1 -4 1 -1 -3 2 -2 -3 3 0 -4 0 0 -1 2 -4 -2 -1 0 -1 -8",
        "0 0 2 2 -4 1 1 0 2 -2 -3 1 -2 -3 0 1 0 -4 -2 -2 2 1 0 -8",
        "0 -1 2 4 -5 2 3 1 1 -2 -4 0 -3 -6 -1 0 0 -7 -4 -2 3 3 -1 -8",
        "-2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3 0 -2 -8 0 -2 -4 -5 -3 -8",
        "0 1 1 2 -5 4 2 -1 3 -2 -2 1 -1 -5 0 -1 -1 -5 -4 -2 1 3 -1 -8",
        "0 -1 1 3 -5 2 4 0 1 -2 -3 0 -2 -5 -1 0 0 -7 -4 -2 3 3 -1 -8",
        "1 -3 0 1 -3 -1 0 5 -2 -3 -4 -2 -3 -5 0 1 0 -7 -5 -1 0 0 -1 -8",
        "-1 2 2 1 -3 3 1 -2 6 -2 -2 0 -2 -2 0 -1 -1 -3 0 -2 1 2 -1 -8",
        "-1 -2 -2 -2 -2 -2 -2 -3 -2 5 2 -2 2 1 -2 -1 0 -5 -1 4 -2 -2 -1 -8",
        "-2 -3 -3 -4 -6 -2 -3 -4 -2 2 6 -3 4 2 -3 -3 -2 -2 -1 2 -3 -3 -1 -8",
        "-1 3 1 0 -5 1 0 -2 0 -2 -3 5 0 -5 -1 0 0 -3 -4 -2 1 0 -1 -8",
        "-1 0 -2 -3 -5 -1 -2 -3 -2 2 4 0 6 0 -2 -2 -1 -4 -2 2 -2 -2 -1 -8",
        "-3 -4 -3 -6 -4 -5 -5 -5 -2 1 2 -5 0 9 -5 -3 -3 0 7 -1 -4 -5 -2 -8",
        "1 0 0 -1 -3 0 -1 0 0 -2 -3 -1 -2 -5 6 1 0 -6 -5 -1 -1 0 -1 -8",
        "1 0 1 0 0 -1 0 1 -1 -1 -3 0 -2 -3 1 2 1 -2 -3 -1 0 0 0 -8",
        "1 -1 0 0 -2 -1 0 0 -1 0 -2 0 -1 -3 0 1 3 -5 -3 0 0 -1 0 -8",
        "-6 2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4 0 -6 -2 -5 17 0 -6 -5 -6 -4 -8",
        "-3 -4 -2 -4 0 -4 -4 -5 0 -1 -1 -4 -2 7 -5 -3 -3 0 10 -2 -3 -4 -2 -8",
        "0 -2 -2 -2 -2 -2 -2 -1 -2 4 2 -2 2 -1 -1 -1 0 -6 -2 4 -2 -2 -1 -8",
        "0 -1 2 3 -4 1 3 0 1 -2 -3 1 -2 -4 -1 0 0 -5 -3 -2 3 2 -1 -8",
        "0 0 1 3 -5 3 3 0 2 -2 -3 0 -2 -5 0 0 -1 -6 -4 -2 2 3 -1 -8",
        "0 -1 0 -1 -3 -1 -1 -1 -1 -1 -1 -1 -1 -2 -1 0 0 -4 -2 -1 -1 -1 -1 -8",
        "-8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 1"
    };

    private static short[] ClustalWGon250Triangle() => new short[]
    {
        24,
        0, 0,
        5, 0, 115,
        -3, 0, -32, 47,
        0, 0, -30, 27, 36,
        -23, 0, -8, -45, -39, 70,
        5, 0, -20, 1, -8, -52, 66,
        -8, 0, -13, 4, 4, -1, -14, 60,
        -8, 0, -11, -38, -27, 10, -45, -22, 40,
        -4, 0, -28, 5, 12, -33, -11, 6, -21, 32,
        -12, 0, -15, -40, -28, 20, -44, -19, 28, -21, 40,
        -7, 0, -9, -30, -20, 16, -35, -13, 25, -14, 28, 43,
        -3, 0, -18, 22, 9, -31, 4, 12, -28, 8, -30, -22, 38,
        3, 0, -31, -7, -5, -38, -16, -11, -26, -6, -23, -24, -9, 76,
        -2, 0, -24, 9, 17, -26, -10, 12, -19, 15, -16, -10, 7, -2, 27,
        -6, 0, -22, -3, 4, -32, -10, 6, -24, 27, -22, -17, 3, -9, 15, 47,
        11, 0, 1, 5, 2, -28, 4, -2, -18, 1, -21, -14, 9, 4, 2, -2, 22,
        6, 0, -5, 0, -1, -22, -11, -3, -6, 1, -13, -6, 5, 1, 0, -2, 15, 25,
        1, 0, 0, -29, -19, 1, -33, -20, 31, -17, 18, 16, -22, -18, -15, -20, -10, 0, 34,
        -36, 0, -10, -52, -43, 36, -40, -8, -18, -35, -7, -10, -36, -50, -27, -16, -33, -35, -26, 142,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        -22, 0, -5, -28, -27, 51, -40, 22, -7, -21, 0, -2, -14, -31, -17, -18, -19, -19, -11, 41, 0, 78,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    private static short[] ClustalWBlosum30Triangle() => new short[]
    {
        4,
        0, 5,
        -3, -2, 17,
        0, 5, -3, 9,
        0, 0, 1, 1, 6,
        -2, -3, -3, -5, -4, 10,
        0, 0, -4, -1, -2, -3, 8,
        -2, -2, -5, -2, 0, -3, -3, 14,
        0, -2, -2, -4, -3, 0, -1, -2, 6,
        0, 0, -3, 0, 2, -1, -1, -2, -2, 4,
        -1, -1, 0, -1, -1, 2, -2, -1, 2, -2, 4,
        1, -2, -2, -3, -1, -2, -2, 2, 1, 2, 2, 6,
        0, 4, -1, 1, -1, -1, 0, -1, 0, 0, -2, 0, 8,
        -1, -2, -3, -1, 1, -4, -1, 1, -3, 1, -3, -4, -3, 11,
        1, -1, -2, -1, 2, -3, -2, 0, -2, 0, -2, -1, -1, 0, 8,
        -1, -2, -2, -1, -1, -1, -2, -1, -3, 1, -2, 0, -2, -1, 3, 8,
        1, 0, -2, 0, 0, -1, 0, -1, -1, 0, -2, -2, 0, -1, -1, -1, 4,
        1, 0, -2, -1, -2, -2, -2, -2, 0, -1, 0, 0, 1, 0, 0, -3, 2, 5,
        1, -2, -2, -2, -3, 1, -3, -3, 4, -2, 1, 0, -2, -4, -3, -1, -1, 1, 5,
        -5, -5, -2, -4, -1, 1, 1, -5, -3, -2, -2, -3, -7, -3, -1, 0, -3, -5, -3, 20,
        0, -1, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, -2, -1,
        -4, -3, -6, -1, -2, 3, -3, 0, -1, -1, 3, -1, -4, -2, -1, 0, -2, -1, 1, 5, -1, 9,
        0, 0, 0, 0, 5, -4, -2, 0, -3, 1, -1, -1, -1, 0, 4, 0, -1, -1, -3, -1, 0, -2, 4,
    };

    private static short[] ClustalWPam350Triangle() => new short[]
    {
        2,
        1, 3,
        -2, -5, 18,
        1, 3, -6, 4,
        1, 3, -6, 4, 4,
        -4, -5, -5, -6, -6, 13,
        2, 1, -4, 1, 1, -6, 5,
        -1, 1, -4, 1, 1, -2, -2, 7,
        0, -2, -3, -2, -2, 2, -2, -2, 5,
        -1, 1, -6, 1, 0, -6, -1, 1, -2, 5,
        -2, -4, -7, -4, -4, 3, -4, -2, 4, -3, 8,
        -1, -2, -6, -3, -2, 1, -3, -2, 3, 0, 5, 6,
        0, 2, -4, 2, 2, -4, 1, 2, -2, 1, -3, -2, 2,
        1, 0, -3, 0, 0, -5, 0, 0, -2, -1, -3, -2, 0, 6,
        0, 2, -6, 2, 3, -5, -1, 3, -2, 1, -2, -1, 1, 1, 4,
        -1, 0, -4, -1, 0, -5, -2, 2, -2, 4, -3, 0, 1, 0, 2, 7,
        1, 1, 0, 1, 0, -4, 1, -1, -1, 0, -3, -2, 1, 1, 0, 0, 1,
        1, 0, -2, 0, 0, -3, 1, -1, 0, 0, -2, -1, 1, 1, 0, -1, 1, 2,
        0, -2, -2, -2, -2, -1, -1, -2, 4, -2, 3, 2, -2, -1, -2, -3, -1, 0, 5,
        -7, -6, -10, -8, -8, 1, -8, -3, -6, -4, -2, -5, -5, -7, -5, 4, -3, -6, -7, 27,
        0, 0, -3, -1, 0, -2, -1, 0, 0, -1, -1, 0, 0, 0, 0, -1, 0, 0, 0, -5, -1,
        -4, -4, 1, -5, -5, 11, -6, 0, 0, -5, 0, -2, -3, -6, -5, -5, -3, -3, -2, 1, -2, 14,
        0, 2, -6, 3, 3, -6, 0, 2, -2, 1, -3, -2, 2, 0, 3, 1, 0, 0, -2, -7, 0, -5, 3,
    };

    private static short[] ClustalWIdTriangle() => new short[]
    {
        10,
        0, 10,
        0, 0, 10,
        0, 0, 0, 10,
        0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,
    };

    private static string[] EdnaFullRows() => new[]
    {
        "5 -4 -4 -4 -4 1 1 -4 -4 1 -4 -1 -1 -1 -2 -4",
        "-4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2 5",
        "-4 -4 5 -4 1 -4 1 -4 1 -4 -1 -1 -4 -1 -2 -4",
        "-4 -4 -4 5 1 -4 -4 1 -4 1 -1 -1 -1 -4 -2 -4",
        "-4 -4 1 1 -1 -4 -2 -2 -2 -2 -1 -1 -3 -3 -1 -4",
        "1 1 -4 -4 -4 -1 -2 -2 -2 -2 -3 -3 -1 -1 -1 1",
        "1 -4 1 -4 -2 -2 -1 -4 -2 -2 -3 -1 -3 -1 -1 -4",
        "-4 1 -4 1 -2 -2 -4 -1 -2 -2 -1 -3 -1 -3 -1 1",
        "-4 1 1 -4 -2 -2 -2 -2 -1 -4 -1 -3 -3 -1 -1 1",
        "1 -4 -4 1 -2 -2 -2 -2 -4 -1 -3 -1 -1 -3 -1 -4",
        "-4 -1 -1 -1 -1 -3 -3 -1 -1 -3 -1 -2 -2 -2 -1 -1",
        "-1 -4 -1 -1 -1 -3 -1 -3 -3 -1 -2 -1 -2 -2 -1 -4",
        "-1 -1 -4 -1 -3 -1 -3 -1 -3 -1 -2 -2 -1 -2 -1 -1",
        "-1 -1 -1 -4 -3 -1 -1 -3 -1 -3 -2 -2 -2 -1 -1 -1",
        "-2 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2",
        "-4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2 5"
    };
}

/// <summary>
/// Result of multiple sequence alignment.
/// </summary>
public sealed record MultipleAlignmentResult(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore)
{
    private readonly double? _exactTotalScore;

    /// <summary>
    /// Sum-of-pairs score as a real number. Equal to <see cref="TotalScore"/> for every integer-scored (<see cref="ScoringMatrix"/>)
    /// aligner; the <see cref="MatrixScoring"/> aligners (substitution matrix, possibly fractional gap penalties such as 0.5) set it to
    /// the exact score while <see cref="TotalScore"/> then holds that value rounded half away from zero (as
    /// <see cref="AlignmentResult.ExactScore"/> for pairwise alignments).
    /// </summary>
    public double ExactTotalScore
    {
        get => _exactTotalScore ?? TotalScore;
        init => _exactTotalScore = value;
    }

    public static MultipleAlignmentResult Empty => new(Array.Empty<string>(), "", 0);
}

/// <summary>
/// Source of the pairwise distances that build the guide tree of the progressive, iterative and consistency
/// multiple aligners.
/// </summary>
public enum MsaDistanceMode
{
    /// <summary>
    /// Default (historical behaviour): d = 1 − identities / min(len1, len2) of the repository's own global
    /// alignment under the caller's <see cref="ScoringMatrix"/> (linear-gap Needleman–Wunsch; the affine Gotoh alignment when
    /// <see cref="MsaOptions.GapModel"/> is <see cref="MsaGapModel.Affine"/>).
    /// </summary>
    Linear = 0,

    /// <summary>
    /// ClustalW 2.1 "slow" pairwise stage (Thompson, Higgins &amp; Gibson 1994): affine Smith–Waterman–Gotoh/Myers–Miller
    /// alignment with ClustalW's own DNA matrix (IUB, ×100) and pairwise gap penalties (open 15, extend 6.66), d = 1 −
    /// identities / min(len1, len2) of that alignment (<c>SequenceAligner.ClustalWPairwiseAlign</c>). The caller's scoring matrix is
    /// then used only for the profile alignment, not for the distances.
    /// </summary>
    ClustalWAffine = 1,

    /// <summary>
    /// ClustalW 2.1 "fast" pairwise stage (<c>-QUICKTREE</c>, <c>FastPairwiseAlign.cpp</c>; Wilbur &amp; Lipman 1983 k-tuple diagonal
    /// scoring): d = (100 − 100·fragmentScore / min(len1, len2)) / 100 in double precision, with the k-tuple / window / top-diagonals /
    /// gap parameters of <see cref="MsaOptions.QuickPairwise"/> (ClustalW's DNA defaults 2 / 4 / 4 / 5). The fragment score is not an
    /// alignment identity, so the distances (and therefore the guide tree) differ from <see cref="ClustalWAffine"/>
    /// (<c>SequenceAligner.ClustalWQuickPairwiseAlign</c>). The caller's scoring matrix is used only for the profile alignment.
    /// </summary>
    ClustalWQuick = 2,
}

/// <summary>Tree-building method of the guide tree of the progressive, iterative and consistency multiple aligners.</summary>
public enum MsaGuideTree
{
    /// <summary>
    /// Default (historical behaviour): UPGMA (Sokal &amp; Michener 1958; the Feng–Doolittle 1987 / CLUSTAL-1 choice) from
    /// <c>PhylogeneticAnalyzer</c>; ties between distances equal within 1e-9 go to the lowest (i, j) index pair.
    /// </summary>
    Upgma = 0,

    /// <summary>
    /// Neighbour-joining (Saitou &amp; Nei 1987) as ClustalW 2.1 uses it by default (<c>-CLUSTERING=NJ</c>; also the default
    /// tree of T-Coffee / Clustal Omega-style pipelines): the unrooted NJ tree is built with ClustalW's own neighbour-joining
    /// arithmetic (<c>NJTree::generateTree</c>; same topology as Saitou–Nei, exact ties resolved as clustalw does), re-rooted
    /// with ClustalW's <c>Tree::reRoot</c> rule and its merges are taken in ClustalW's order (left subtree first). With
    /// <see cref="MsaDistanceMode.ClustalWAffine"/> the distances are ClustalW's single-precision values, so the guide tree
    /// equals the one <c>clustalw -ALIGN</c> writes to its <c>.dnd</c> file.
    /// </summary>
    NeighborJoining = 1,
}

/// <summary>Gap model of the multiple aligners and of their sum-of-pairs objective (<see cref="MultipleAlignmentResult.TotalScore"/>).</summary>
public enum MsaGapModel
{
    /// <summary>
    /// Default (historical behaviour): every gap position scores <see cref="ScoringMatrix.GapExtend"/> and <see cref="ScoringMatrix.GapOpen"/>
    /// is unused; <c>TotalScore</c> is the column sum of pairs (gap–gap columns neutral).
    /// </summary>
    Linear = 0,

    /// <summary>
    /// Affine gaps (Gotoh 1982): a gap of length k costs <see cref="ScoringMatrix.GapOpen"/> + k·<see cref="ScoringMatrix.GapExtend"/>.
    /// <c>TotalScore</c> is the sum-of-pairs score with <b>natural</b> affine gap costs (Altschul 1989: every pair of rows is projected
    /// out of the alignment, gap–gap columns are skipped, and each maximal run of gap columns of the projection pays the opening once);
    /// the star methods align every sequence to the center with the exact affine pairwise aligner, the profile methods use an affine
    /// (Gotoh three-state) profile–profile dynamic program with occupancy-scaled opening/extension.
    /// </summary>
    Affine = 1,
}

/// <summary>
/// Additive options of the multiple aligners (<c>MultipleAlign</c>, <c>MultipleAlignAnchored</c>, <c>MultipleAlignProgressive</c>,
/// <c>MultipleAlignIterative</c>, <c>MultipleAlignConsistency</c>). The default instance reproduces the parameterless behaviour.
/// </summary>
public sealed record MsaOptions
{
    /// <summary>Pairwise distance used for the guide tree of the progressive family (default <see cref="MsaDistanceMode.Linear"/>).</summary>
    public MsaDistanceMode Distance { get; init; } = MsaDistanceMode.Linear;

    /// <summary>Gap model (default <see cref="MsaGapModel.Linear"/>).</summary>
    public MsaGapModel GapModel { get; init; } = MsaGapModel.Linear;

    /// <summary>Guide-tree method of the progressive family (default <see cref="MsaGuideTree.Upgma"/>).</summary>
    public MsaGuideTree GuideTree { get; init; } = MsaGuideTree.Upgma;

    /// <summary>
    /// Parameters of the fast pairwise stage, used only with <see cref="MsaDistanceMode.ClustalWQuick"/> (default: ClustalW's DNA defaults).
    /// </summary>
    public ClustalWQuickOptions QuickPairwise { get; init; } = new();

    /// <summary>
    /// ClustalW PROTEIN pairwise parameters (matrix, <c>-PWGAPOPEN</c> / <c>-PWGAPEXT</c>, fast-stage options), used by
    /// <see cref="MsaDistanceMode.ClustalWAffine"/> / <see cref="MsaDistanceMode.ClustalWQuick"/> when the sequences are scored with a protein
    /// <see cref="MatrixScoring"/> (the DNA stage is used for nucleotide matrices and scalar scoring).
    /// </summary>
    public ClustalWProteinPairwiseOptions ClustalWProtein { get; init; } = new();

    /// <summary>
    /// Progressive stage of <c>MultipleAlignConsistency</c> (default <see cref="MsaConsistencyEngine.PaperExtension"/>, the published formula).
    /// </summary>
    public MsaConsistencyEngine ConsistencyEngine { get; init; } = MsaConsistencyEngine.PaperExtension;

    /// <summary>
    /// <c>MultipleAlignExact</c> only: complexity guard — the largest number of dynamic-programming values (surviving lattice cells ×
    /// gap-history states) the exact sum-of-pairs program may hold (default 50 000 000 ≈ 200 MB of 4-byte scores); above it the
    /// call throws a clear <see cref="ArgumentException"/> instead of exhausting memory.
    /// </summary>
    public long ExactMaxStateCells { get; init; } = 50_000_000;

    /// <summary>
    /// <c>MultipleAlignExact</c> only: prune the k-dimensional lattice with the Carrillo–Lipman (1988) bound (default <c>true</c>;
    /// the result is identical, only the number of cells changes). <c>false</c> evaluates the whole lattice.
    /// </summary>
    public bool ExactCarrilloLipman { get; init; } = true;
}

/// <summary>Progressive stage of the consistency-based (T-Coffee) multiple aligner.</summary>
public enum MsaConsistencyEngine
{
    /// <summary>The paper's objective: extended library weight (direct + sum of triplet minima) summed over every cross pair of two columns, zero gap penalty, full DP (default).</summary>
    PaperExtension = 0,

    /// <summary>
    /// The T-Coffee 13.41 program's progressive stage (<c>SequenceAligner.MultipleAlignTCoffee</c>: library relaxation, sparse linked-pair-wise DP).
    /// Identical to the program for the same primary library and guide tree; the library is this repository's (global + local alignments, weight 10 x percent identity) and the tree the configured guide tree.
    /// </summary>
    TCoffeeLinked = 1,
}

/// <summary>
/// Parameters of ClustalW's fast (k-tuple) pairwise stage for DNA: <c>-KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -SCORE</c>
/// (<c>FastPairwiseAlign.cpp</c>, ClustalW 2.1); defaults are the DNA defaults of <c>UserParameters.cpp</c>.
/// </summary>
public sealed record ClustalWQuickOptions
{
    /// <summary>Word size (<c>-KTUPLE</c>), 1–4 for DNA; default 2.</summary>
    public int KTuple { get; init; } = 2;

    /// <summary>Number of diagonals scanned on each side of a top diagonal (<c>-WINDOW</c>); default 4.</summary>
    public int Window { get; init; } = 4;

    /// <summary>Number of best diagonals (<c>-TOPDIAGS</c>); default 4.</summary>
    public int TopDiagonals { get; init; } = 4;

    /// <summary>
    /// Gap penalty of the fragment chain (<c>-PAIRGAP</c>); default 5 (ClustalW's DNA default). Note: ClustalW's command line sets
    /// the pair gap to <see cref="KTuple"/> + 4 whenever <c>-KTUPLE</c> is given (and a later <c>-PAIRGAP</c> overrides it), so
    /// <c>-KTUPLE=2</c> alone means <c>PairGap = 6</c>, not 5; callers that want the command-line behaviour pass that value.
    /// </summary>
    public int PairGap { get; init; } = 5;

    /// <summary><c>-SCORE=PERCENT</c> (default): score = 100·fragmentScore / min(len1, len2); <c>false</c> = <c>-SCORE=ABSOLUTE</c> (the raw fragment score).</summary>
    public bool Percent { get; init; } = true;
}

/// <summary>Residue type of ClustalW's pairwise stages (<c>-TYPE=DNA | PROTEIN</c>).</summary>
public enum ClustalWSequenceType
{
    /// <summary>Nucleotides (IUB matrix, gap 15 / 6.66, k-tuple 2).</summary>
    Dna = 0,

    /// <summary>Amino acids (GONNET / BLOSUM30 / PAM350 / identity matrix, gap 10 / 0.1, k-tuple 1).</summary>
    Protein = 1,
}

/// <summary>The pairwise substitution matrix of ClustalW's PROTEIN full pairwise stage (<c>-PWMATRIX</c>).</summary>
public enum ClustalWPairwiseMatrix
{
    /// <summary><c>GONNET</c> (default): <c>gon250mt</c> × intScale 10 (<see cref="SubstitutionMatrix.ClustalWGonnet250"/>).</summary>
    Gonnet = 0,

    /// <summary><c>BLOSUM</c>: <c>blosum30mt</c> × 100 (<see cref="SubstitutionMatrix.ClustalWBlosum30"/>).</summary>
    Blosum = 1,

    /// <summary><c>PAM</c>: <c>pam350mt</c> × 100 (<see cref="SubstitutionMatrix.ClustalWPam350"/>).</summary>
    Pam = 2,

    /// <summary><c>ID</c>: identity matrix (<see cref="SubstitutionMatrix.ClustalWIdentity"/>) × 100.</summary>
    Identity = 3,
}

/// <summary>
/// Parameters of ClustalW 2.1's PROTEIN pairwise stages (<c>clustalw -TYPE=PROTEIN</c>): the full dynamic-programming stage (<c>-PWMATRIX</c>,
/// <c>-PWGAPOPEN</c>, <c>-PWGAPEXT</c>) and the fast k-tuple stage (<c>-QUICKTREE</c>: <c>-KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -SCORE</c>); defaults are
/// the PROTEIN defaults of <c>UserParameters.cpp</c>.
/// </summary>
public sealed record ClustalWProteinPairwiseOptions
{
    /// <summary>Pairwise matrix (<c>-PWMATRIX</c>, default <see cref="ClustalWPairwiseMatrix.Gonnet"/>).</summary>
    public ClustalWPairwiseMatrix Matrix { get; init; } = ClustalWPairwiseMatrix.Gonnet;

    /// <summary>Pairwise gap-opening penalty (<c>-PWGAPOPEN</c>, PROTEIN default 10.0; a positive penalty, stored in single precision like ClustalW).</summary>
    public double GapOpen { get; init; } = 10.0;

    /// <summary>Pairwise gap-extension penalty (<c>-PWGAPEXT</c>, PROTEIN default 0.1).</summary>
    public double GapExtend { get; init; } = 0.1;

    /// <summary>
    /// Parameters of the fast stage; PROTEIN defaults <c>KTuple</c> 1 (1-2), <c>Window</c> 5, <c>TopDiagonals</c> 5, <c>PairGap</c> 3 (ClustalW's command line sets the pair
    /// gap to <c>ktuple + 3</c> when <c>-KTUPLE</c> is given).
    /// </summary>
    public ClustalWQuickOptions Quick { get; init; } = new() { KTuple = 1, Window = 5, TopDiagonals = 5, PairGap = 3 };
}

/// <summary>
/// Result of ClustalW's fast k-tuple pairwise stage (<c>SequenceAligner.ClustalWQuickPairwiseAlign</c>, port of
/// <c>FastPairwiseAlign::pairAlign</c>, ClustalW 2.1).
/// </summary>
/// <param name="FragmentScore">Score of the best k-tuple fragment chain (<c>accum[0][maxSoFar]</c>; 0 when no fragment is found).</param>
/// <param name="Score">The number ClustalW prints as <c>Score:</c> (<c>%lg</c>): the percent 100·FragmentScore / min(len1, len2), or the raw fragment score with <c>-SCORE=ABSOLUTE</c>.</param>
/// <param name="Distance">(100 − Score) / 100 (double precision): ClustalW's guide-tree distance (1 when not aligned; with the absolute score it can be negative, as in ClustalW).</param>
/// <param name="Aligned"><c>Score &gt; 0.1</c>: ClustalW prints <c>Aligned. Score: …</c>, otherwise <c>Not Aligned</c>.</param>
/// <param name="PartialAlignment">The fragment list exceeded ClustalW's bound (2 · 2 · longest sequence) and the scan stopped early (<c>(Partial alignment)</c>).</param>
public sealed record ClustalWQuickPairwiseResult(
    int FragmentScore,
    double Score,
    double Distance,
    bool Aligned,
    bool PartialAlignment);

/// <summary>Pairwise stage that produces the distance matrix of the ClustalW multiple alignment (<see cref="ClustalWOptions.PairwiseMode"/>).</summary>
public enum ClustalWPairwiseMode
{
    /// <summary>ClustalW's default "slow" full dynamic-programming stage (<c>FullPairwiseAlign</c>; <c>SequenceAligner.ClustalWPairwiseAlign</c>).</summary>
    Full = 0,

    /// <summary>ClustalW's fast k-tuple stage (<c>-QUICKTREE</c>, <c>FastPairwiseAlign</c>; <c>SequenceAligner.ClustalWQuickPairwiseAlign</c>).</summary>
    Quick = 1,
}

/// <summary>
/// Parameters of <c>SequenceAligner.MultipleAlignClustalW</c>, the port of ClustalW 2.1's multiple-alignment stage for DNA
/// (<c>clustalw -TYPE=DNA -ALIGN</c>). Defaults are ClustalW's DNA defaults.
/// </summary>
public sealed record ClustalWOptions
{
    /// <summary>Gap-opening penalty of the profile alignment (<c>-GAPOPEN</c>, DNA default 15.0; stored in single precision like ClustalW).</summary>
    public double GapOpen { get; init; } = 15.0;

    /// <summary>Gap-extension penalty of the profile alignment (<c>-GAPEXT</c>, DNA default 6.66).</summary>
    public double GapExtend { get; init; } = 6.66;

    /// <summary>Pairwise stage that produces the distances of the guide tree (default <see cref="ClustalWPairwiseMode.Full"/>).</summary>
    public ClustalWPairwiseMode PairwiseMode { get; init; } = ClustalWPairwiseMode.Full;

    /// <summary>Pairwise gap-opening penalty of the full pairwise stage (<c>-PWGAPOPEN</c>, default 15).</summary>
    public double PairwiseGapOpen { get; init; } = 15.0;

    /// <summary>Pairwise gap-extension penalty of the full pairwise stage (<c>-PWGAPEXT</c>, default 6.66).</summary>
    public double PairwiseGapExtend { get; init; } = 6.66;

    /// <summary>Parameters of the fast pairwise stage (<see cref="ClustalWPairwiseMode.Quick"/> only).</summary>
    public ClustalWQuickOptions Quick { get; init; } = new();

    /// <summary>
    /// Weight of transitions (A–G, C–T/U) in the profile substitution matrix (<c>-TRANSWEIGHT</c>, default 0.5; 0 = transitions
    /// score as mismatches, 1 = as matches).
    /// </summary>
    public float TransitionWeight { get; init; } = 0.5f;

    /// <summary>
    /// <c>-MAXDIV</c> (default 30): a sequence whose highest percent identity to any other sequence (as derived from the guide
    /// tree) is at most this value is delayed and added after the rest has been aligned (delayed divergent sequences).
    /// </summary>
    public int DivergenceCutoff { get; init; } = 30;

    /// <summary><c>-NOWEIGHTS</c>: no sequence weights from the guide tree (all sequences weigh the same).</summary>
    public bool NoWeights { get; init; }

    /// <summary><c>-NEGATIVE</c>: use the substitution matrix as it is (negative mismatch scores) instead of shifting it to non-negative scores.</summary>
    public bool NegativeMatrix { get; init; }
}

/// <summary>
/// Result of ClustalW's "slow" pairwise alignment (<c>SequenceAligner.ClustalWPairwiseAlign</c>, port of
/// <c>FullPairwiseAlign::pairwiseAlign</c>, ClustalW 2.1): a local affine alignment whose identities define ClustalW's percent
/// identity and distance.
/// </summary>
/// <param name="Score">Alignment score in ClustalW's integer units (matrix ×100, gap penalties ×100).</param>
/// <param name="Identities">Aligned residue pairs with identical residue codes in the local alignment.</param>
/// <param name="PercentIdentity">100 · identities / min(len1, len2) (0 when a sequence is empty).</param>
/// <param name="Distance">1 − identities / min(len1, len2) (1 when a sequence is empty): ClustalW's guide-tree distance.</param>
/// <param name="ClustalWScore">The integer "Score:" ClustalW prints per pair: (int)(100·identities / min(len1, len2)) in single precision.</param>
/// <param name="AlignedSequence1">Aligned region of sequence 1 (gap '-').</param>
/// <param name="AlignedSequence2">Aligned region of sequence 2.</param>
/// <param name="Start1">0-based first aligned residue of sequence 1 (-1 when the alignment is empty).</param>
/// <param name="End1">0-based last aligned residue of sequence 1 (-1 when empty).</param>
/// <param name="Start2">0-based first aligned residue of sequence 2 (-1 when empty).</param>
/// <param name="End2">0-based last aligned residue of sequence 2 (-1 when empty).</param>
public sealed record ClustalWPairwiseResult(
    int Score,
    int Identities,
    double PercentIdentity,
    double Distance,
    int ClustalWScore,
    string AlignedSequence1,
    string AlignedSequence2,
    int Start1,
    int End1,
    int Start2,
    int End2);
