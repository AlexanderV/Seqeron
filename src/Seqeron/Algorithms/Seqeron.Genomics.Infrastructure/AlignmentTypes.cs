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
/// Type of alignment performed.
/// </summary>
public enum AlignmentType
{
    Global,
    Local,
    SemiGlobal
}

/// <summary>
/// Result of pairwise sequence alignment.
/// </summary>
/// <param name="AlignedSequence1">Aligned (gapped, '-') segment of sequence 1.</param>
/// <param name="AlignedSequence2">Aligned (gapped, '-') segment of sequence 2.</param>
/// <param name="Score">Alignment score under the scoring model used.</param>
/// <param name="AlignmentType">Kind of alignment that produced the result.</param>
/// <param name="StartPosition1">0-based inclusive start of the aligned segment in sequence 1
/// (0 for global/semi-global; −1 for a local result with no positive-scoring region).</param>
/// <param name="StartPosition2">0-based inclusive start of the aligned segment in sequence 2 (see <paramref name="StartPosition1"/>).</param>
/// <param name="EndPosition1">0-based inclusive end of the aligned segment in sequence 1
/// (length − 1 for global/semi-global; −1 for an empty local result).</param>
/// <param name="EndPosition2">0-based inclusive end of the aligned segment in sequence 2 (see <paramref name="EndPosition1"/>).</param>
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
public readonly record struct AlignmentStatistics(
    int Matches,
    int Mismatches,
    int Gaps,
    int AlignmentLength,
    double Identity,
    double Similarity,
    double GapPercent,
    int Positives = 0,
    int GapOpens = 0)
{
    public static AlignmentStatistics Empty => new(0, 0, 0, 0, 0, 0, 0);
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
    public static MultipleAlignmentResult Empty => new(Array.Empty<string>(), "", 0);
}
