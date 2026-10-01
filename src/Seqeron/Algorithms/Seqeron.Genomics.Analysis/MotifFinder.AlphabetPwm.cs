using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Position weight matrices over an arbitrary alphabet (protein, RNA, gapped DNA, …) — the generic part of
/// Biopython <c>Bio.motifs</c>: <c>motifs.create(instances, alphabet)</c>, <c>counts.normalize(pseudocounts)</c>,
/// <c>log_odds(background)</c>, <c>consensus</c> / <c>anticonsensus</c>, <c>max</c> / <c>min</c>,
/// <c>mean</c> / <c>std</c>, and the window scoring of <c>calculate</c> / <c>search(both=False)</c>.
/// </summary>
/// <remarks>
/// The log-odds, mean/std and window-scoring arithmetic is the same kernel the DNA
/// <see cref="PositionWeightMatrix"/> uses (<see cref="LogOddsFromCounts"/>, <see cref="PwmMean"/>,
/// <see cref="PwmStd"/>, <c>ScorePwmWindow&lt;TRowIndex&gt;</c>), so for the alphabet "ACGT" both types give the
/// same matrices and scores.
/// </remarks>
public static partial class MotifFinder
{
    /// <summary>
    /// Shared PWM column-extremum kernel (Biopython 1.88 <c>Bio/motifs/matrix.py</c>
    /// <c>PositionSpecificScoringMatrix.max</c>/<c>min</c>): Σ over columns of Python's builtin
    /// <c>max</c>/<c>min</c> of the column — seeded with the first row and replaced only on a strict
    /// <c>&gt;</c>/<c>&lt;</c>, so an all −∞ column contributes −∞ (not <c>double.MinValue</c>), a NaN in
    /// the first row propagates and a NaN in a later row is ignored. Used by both
    /// <see cref="PositionWeightMatrix"/> and <see cref="AlphabetPositionWeightMatrix"/>.
    /// </summary>
    internal static double ColumnExtremumSum(double[,] matrix, bool maximum)
    {
        double score = 0.0;
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        for (int j = 0; j < cols; j++)
        {
            double best = matrix[0, j];
            for (int a = 1; a < rows; a++)
            {
                double w = matrix[a, j];
                if (maximum ? w > best : w < best)
                    best = w;
            }

            score += best;
        }

        return score;
    }

    #region Generic-alphabet PWM

    /// <summary>
    /// Creates a log-odds PWM over an arbitrary alphabet from aligned instances — Biopython
    /// <c>motifs.create(instances, alphabet).counts.normalize(pseudocounts=p).log_odds(background)</c>:
    /// W[a,j] = log2( ((c[a,j] + p) / (Σ_a c[a,j] + K·p)) / q[a] ), K = |alphabet|.
    /// </summary>
    /// <param name="sequences">Aligned instances of equal length.</param>
    /// <param name="alphabet">
    /// The alphabet, one row per symbol in this order (e.g. "ACDEFGHIKLMNPQRSTVWY" for protein). Symbols must be
    /// distinct ignoring case; instances and scanned sequences are matched case-insensitively (as the DNA
    /// <see cref="CreatePwm(IEnumerable{string}, double)"/> and Biopython's DNA <c>calculate</c>).
    /// </param>
    /// <param name="pseudocount">Pseudocount per cell (finite, ≥ 0; default 0 = Biopython <c>pseudocounts=None</c>).</param>
    /// <param name="background">Background in alphabet order (finite, &gt; 0, normalised to sum 1); uniform when null.</param>
    /// <param name="ignoreUnknownSymbols">
    /// false (default): a symbol outside the alphabet throws. true: such symbols (gaps, X, …) are not counted, exactly as
    /// Biopython (<c>alignment.frequencies</c> keeps only alphabet letters), so that column's total is smaller.
    /// </param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">Empty collection, null element, unequal lengths, invalid alphabet,
    /// unknown symbol (when not ignored), or a column with no counts and no pseudocounts.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Negative / non-finite pseudocount or invalid background value.</exception>
    public static AlphabetPositionWeightMatrix CreateAlphabetPwm(
        IEnumerable<string> sequences,
        string alphabet,
        double pseudocount = 0.0,
        IReadOnlyList<double>? background = null,
        bool ignoreUnknownSymbols = false)
    {
        ArgumentNullException.ThrowIfNull(alphabet);
        if (!double.IsFinite(pseudocount) || pseudocount < 0)
            throw new ArgumentOutOfRangeException(nameof(pseudocount), pseudocount,
                "Pseudocount must be finite and non-negative.");
        var pseudo = new double[alphabet.Length];
        Array.Fill(pseudo, pseudocount);
        return CreateAlphabetPwm(sequences, alphabet, pseudo, background, ignoreUnknownSymbols);
    }

    /// <summary>
    /// Creates a log-odds PWM over an arbitrary alphabet with per-symbol pseudocounts — Biopython
    /// <c>counts.normalize(pseudocounts={letter: p}).log_odds(background)</c>:
    /// W[a,j] = log2( ((c[a,j] + p[a]) / (Σ_a c[a,j] + Σ_a p[a])) / q[a] ).
    /// </summary>
    /// <param name="sequences">Aligned instances of equal length.</param>
    /// <param name="alphabet">The alphabet (row order); see <see cref="CreateAlphabetPwm(IEnumerable{string}, string, double, IReadOnlyList{double}?, bool)"/>.</param>
    /// <param name="pseudocounts">One pseudocount per alphabet symbol (finite, ≥ 0).</param>
    /// <param name="background">Background in alphabet order; uniform when null.</param>
    /// <param name="ignoreUnknownSymbols">Skip (true) or reject (false) symbols outside the alphabet.</param>
    public static AlphabetPositionWeightMatrix CreateAlphabetPwm(
        IEnumerable<string> sequences,
        string alphabet,
        IReadOnlyList<double> pseudocounts,
        IReadOnlyList<double>? background = null,
        bool ignoreUnknownSymbols = false)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(alphabet);
        ArgumentNullException.ThrowIfNull(pseudocounts);
        var rowIndex = AlphabetRowIndex.Create(alphabet, nameof(alphabet));

        List<string> seqList = MaterializeAligned(sequences, nameof(sequences));
        if (seqList.Count == 0)
            throw new ArgumentException("At least one sequence is required.", nameof(sequences));

        int k = alphabet.Length;
        int length = seqList[0].Length;
        var counts = new double[k, length];
        for (int s = 0; s < seqList.Count; s++)
        {
            string seq = seqList[s];
            for (int i = 0; i < length; i++)
            {
                int row = rowIndex.RowOf(seq[i]);
                if (row < 0)
                {
                    if (ignoreUnknownSymbols)
                        continue;
                    throw new ArgumentException(
                        $"Invalid character '{seq[i]}' at position {i} in sequence {s}: not in the alphabet \"{alphabet}\".",
                        nameof(sequences));
                }

                counts[row, i]++;
            }
        }

        return AlphabetPositionWeightMatrix.FromCounts(alphabet, counts, pseudocounts, background);
    }

    /// <summary>
    /// Score of every window — the generic-alphabet form of Biopython <c>PositionSpecificScoringMatrix.calculate</c>:
    /// element i = Σ_j W[s[i+j], j] for i = 0 … n − m; a window containing a symbol outside the alphabet scores NaN
    /// (Biopython <c>_pwm.c</c> semantics). Case-insensitive.
    /// </summary>
    /// <remarks>Biopython 1.88 <c>calculate</c> itself accepts only the DNA alphabet (it raises for any other); this
    /// method applies the same per-window sum to any alphabet. It returns doubles (Biopython: float32) and an empty array
    /// when the sequence is shorter than the motif.</remarks>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static double[] CalculateAlphabetPwmScores(string sequence, AlphabetPositionWeightMatrix pwm)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(pwm);

        int windows = sequence.Length - pwm.Length + 1;
        if (windows <= 0)
            return Array.Empty<double>();

        var scores = new double[windows];
        for (int i = 0; i < windows; i++)
            scores[i] = ScorePwmWindow(sequence, i, pwm.RawMatrix, pwm.Length, pwm.RowIndex);
        return scores;
    }

    /// <summary>
    /// Forward scan with a generic-alphabet PWM — Biopython <c>pssm.search(sequence, threshold, both=False)</c>:
    /// every window with score ≥ <paramref name="threshold"/> in ascending position order (NaN windows never match).
    /// </summary>
    /// <param name="sequence">Sequence to scan (case-insensitive).</param>
    /// <param name="pwm">The matrix.</param>
    /// <param name="threshold">Minimum score (inclusive).</param>
    /// <returns>Hits; <see cref="MotifMatch.MatchedSequence"/> is the window as given, <see cref="MotifMatch.Pattern"/> the consensus.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static IEnumerable<MotifMatch> ScanWithAlphabetPwm(
        string sequence, AlphabetPositionWeightMatrix pwm, double threshold = 0.0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(pwm);
        return ScanWithAlphabetPwmCore(sequence, pwm, threshold);
    }

    private static IEnumerable<MotifMatch> ScanWithAlphabetPwmCore(
        string sequence, AlphabetPositionWeightMatrix pwm, double threshold)
    {
        int m = pwm.Length;
        for (int i = 0; i <= sequence.Length - m; i++)
        {
            double score = ScorePwmWindow(sequence, i, pwm.RawMatrix, m, pwm.RowIndex);
            if (score >= threshold)
                yield return new MotifMatch(i, sequence.Substring(i, m), pwm.Consensus, score);
        }
    }

    /// <summary>
    /// Case-insensitive symbol → row map of an alphabet: ASCII through a 128-entry table, other symbols through a
    /// dictionary on the invariant upper case.
    /// </summary>
    internal readonly struct AlphabetRowIndex : IPwmRowIndex
    {
        private readonly int[] _ascii;
        private readonly Dictionary<char, int>? _other;

        private AlphabetRowIndex(int[] ascii, Dictionary<char, int>? other)
        {
            _ascii = ascii;
            _other = other;
        }

        internal static AlphabetRowIndex Create(string alphabet, string paramName)
        {
            if (alphabet.Length == 0)
                throw new ArgumentException("Alphabet must contain at least one symbol.", paramName);

            var ascii = new int[128];
            Array.Fill(ascii, -1);
            Dictionary<char, int>? other = null;
            var seen = new HashSet<char>();
            for (int r = 0; r < alphabet.Length; r++)
            {
                char c = alphabet[r];
                if (char.IsSurrogate(c))
                    throw new ArgumentException("Alphabet symbols must be single UTF-16 code units.", paramName);
                char folded = char.ToUpperInvariant(c);
                if (!seen.Add(folded))
                    throw new ArgumentException(
                        $"Alphabet symbol '{c}' occurs more than once (symbols are compared ignoring case).", paramName);

                if (folded < 128)
                {
                    ascii[folded] = r;
                    char lower = char.ToLowerInvariant(folded);
                    if (lower < 128)
                        ascii[lower] = r;
                }
                else
                {
                    (other ??= new Dictionary<char, int>())[folded] = r;
                }
            }

            return new AlphabetRowIndex(ascii, other);
        }

        public int RowOf(char symbol)
        {
            if (symbol < 128)
                return _ascii[symbol];
            if (_other is null)
                return -1;
            return _other.TryGetValue(char.ToUpperInvariant(symbol), out int row) ? row : -1;
        }
    }

    #endregion
}

/// <summary>
/// A log-odds position weight matrix over an arbitrary alphabet (Biopython
/// <c>Bio.motifs.matrix.PositionSpecificScoringMatrix</c> with any <c>alphabet</c>), K = |alphabet| rows × L columns.
/// Build it with <see cref="MotifFinder.CreateAlphabetPwm(IEnumerable{string}, string, double, IReadOnlyList{double}?, bool)"/>
/// or <see cref="FromCounts(string, double[,], double, IReadOnlyList{double}?)"/>; score with
/// <see cref="MotifFinder.CalculateAlphabetPwmScores"/> / <see cref="MotifFinder.ScanWithAlphabetPwm"/>.
/// </summary>
public sealed class AlphabetPositionWeightMatrix
{
    private readonly double[,] _matrix;

    /// <summary>Creates a matrix from K × L log-odds (rows in alphabet order).</summary>
    /// <param name="alphabet">Distinct symbols (case-insensitive), one per row.</param>
    /// <param name="matrix">Log-odds, |alphabet| × L; ±∞ allowed (unseen symbol / zero pseudocount), NaN not.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">Invalid alphabet, wrong row count or a NaN cell.</exception>
    public AlphabetPositionWeightMatrix(string alphabet, double[,] matrix)
    {
        ArgumentNullException.ThrowIfNull(alphabet);
        ArgumentNullException.ThrowIfNull(matrix);
        RowIndex = MotifFinder.AlphabetRowIndex.Create(alphabet, nameof(alphabet));
        if (matrix.GetLength(0) != alphabet.Length)
            throw new ArgumentException(
                $"Matrix must have {alphabet.Length} rows (one per alphabet symbol); got {matrix.GetLength(0)}.", nameof(matrix));
        foreach (double w in matrix)
        {
            if (double.IsNaN(w))
                throw new ArgumentException("Matrix cells must not be NaN.", nameof(matrix));
        }

        Alphabet = alphabet;
        _matrix = (double[,])matrix.Clone();
        Length = matrix.GetLength(1);
        Consensus = Extremum(maximum: true);
        Anticonsensus = Extremum(maximum: false);
        MaxScore = ColumnExtremumSum(maximum: true);
        MinScore = ColumnExtremumSum(maximum: false);
    }

    /// <summary>
    /// Builds a log-odds matrix from a K × L count matrix with a scalar pseudocount — Biopython
    /// <c>FrequencyPositionMatrix(alphabet, counts).normalize(pseudocounts=p).log_odds(background)</c>.
    /// </summary>
    /// <param name="alphabet">Row symbols.</param>
    /// <param name="counts">Counts (finite, ≥ 0), |alphabet| × L.</param>
    /// <param name="pseudocount">Pseudocount per cell (finite, ≥ 0).</param>
    /// <param name="background">Background in alphabet order; uniform when null.</param>
    public static AlphabetPositionWeightMatrix FromCounts(
        string alphabet, double[,] counts, double pseudocount = 0.0, IReadOnlyList<double>? background = null)
    {
        ArgumentNullException.ThrowIfNull(alphabet);
        var pseudo = new double[alphabet.Length];
        Array.Fill(pseudo, pseudocount);
        return FromCounts(alphabet, counts, pseudo, background);
    }

    /// <summary>
    /// Builds a log-odds matrix from a count matrix with per-symbol pseudocounts — Biopython
    /// <c>normalize(pseudocounts={letter: p}).log_odds(background)</c>; the log-odds kernel is the DNA PWM's.
    /// </summary>
    /// <param name="alphabet">Row symbols.</param>
    /// <param name="counts">Counts (finite, ≥ 0), |alphabet| × L.</param>
    /// <param name="pseudocounts">One pseudocount per symbol (finite, ≥ 0).</param>
    /// <param name="background">Background in alphabet order (finite, &gt; 0, normalised); uniform when null.</param>
    /// <exception cref="ArgumentException">Shape mismatch, or a column whose counts and pseudocounts are all 0
    /// (Biopython divides by a zero total).</exception>
    /// <exception cref="ArgumentOutOfRangeException">Negative / non-finite count, pseudocount or background value.</exception>
    public static AlphabetPositionWeightMatrix FromCounts(
        string alphabet, double[,] counts, IReadOnlyList<double> pseudocounts, IReadOnlyList<double>? background = null)
    {
        ArgumentNullException.ThrowIfNull(alphabet);
        ArgumentNullException.ThrowIfNull(counts);
        ArgumentNullException.ThrowIfNull(pseudocounts);
        MotifFinder.AlphabetRowIndex.Create(alphabet, nameof(alphabet));
        int k = alphabet.Length;
        MotifFinder.ValidateCountMatrix(counts, k, $"Count matrix must have {k} rows (one per alphabet symbol).");
        double[,] logOdds = MotifFinder.CountsToLogOdds(counts, pseudocounts,
            $"Exactly {k} pseudocounts (one per alphabet symbol) are required.", () => ResolveBackground(background, k));
        return new AlphabetPositionWeightMatrix(alphabet, logOdds);
    }

    /// <summary>The alphabet (row order).</summary>
    public string Alphabet { get; }

    /// <summary>Copy of the K × L log-odds matrix (rows in <see cref="Alphabet"/> order).</summary>
    public double[,] GetMatrix() => (double[,])_matrix.Clone();

    /// <summary>Motif length L.</summary>
    public int Length { get; }

    /// <summary>Biopython <c>consensus</c>: per column the first symbol (alphabet order) with the largest log-odds.</summary>
    public string Consensus { get; }

    /// <summary>Biopython <c>anticonsensus</c>: per column the first symbol (alphabet order) with the smallest log-odds.</summary>
    public string Anticonsensus { get; }

    /// <summary>Biopython <c>max</c>: Σ_j max_a W[a,j] (the consensus score).</summary>
    public double MaxScore { get; }

    /// <summary>Biopython <c>min</c>: Σ_j min_a W[a,j] (the anticonsensus score).</summary>
    public double MinScore { get; }

    internal MotifFinder.AlphabetRowIndex RowIndex { get; }

    /// <summary>The matrix without a defensive copy (scoring kernels only).</summary>
    internal double[,] RawMatrix => _matrix;

    /// <summary>Log-odds of <paramref name="symbol"/> (case-insensitive) at column <paramref name="position"/>.</summary>
    /// <exception cref="ArgumentException">The symbol is not in the alphabet.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The position is outside 0 … L − 1.</exception>
    public double this[char symbol, int position]
    {
        get
        {
            int row = RowIndex.RowOf(symbol);
            if (row < 0)
                throw new ArgumentException($"Symbol '{symbol}' is not in the alphabet \"{Alphabet}\".", nameof(symbol));
            if ((uint)position >= (uint)Length)
                throw new ArgumentOutOfRangeException(nameof(position));
            return _matrix[row, position];
        }
    }

    /// <summary>Expected score of a random background window — Biopython <c>pssm.mean(background)</c>.</summary>
    /// <param name="background">Background in alphabet order; uniform when null.</param>
    public double Mean(IReadOnlyList<double>? background = null)
        => MotifFinder.PwmMean(_matrix, Length, ResolveBackground(background, Alphabet.Length));

    /// <summary>
    /// Discretised score distribution under the motif and background models — Biopython
    /// <c>pssm.distribution(background, precision)</c> (<c>Bio.motifs.thresholds.ScoreDistribution</c>, which loops over
    /// the PSSM's own alphabet, so it applies unchanged to protein / any-alphabet PSSMs); the same kernel as
    /// <see cref="PositionWeightMatrix.ScoreDistribution"/> with K rows.
    /// </summary>
    /// <param name="background">Background in alphabet order; uniform when null.</param>
    /// <param name="precision">Grid points per motif position (Biopython default 10³).</param>
    /// <exception cref="InvalidOperationException">The matrix has a non-finite cell (the grid is undefined).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> × length &lt; 2.</exception>
    public PwmScoreDistribution ScoreDistribution(IReadOnlyList<double>? background = null, int precision = 1000)
    {
        double[] bg = ResolveBackground(background, Alphabet.Length);
        return new PwmScoreDistribution(_matrix, Alphabet.Length, Length, MinScore, MaxScore,
            () => MotifFinder.PwmMean(_matrix, Length, bg), bg, precision);
    }

    /// <summary>Standard deviation of the score of a random background window — Biopython <c>pssm.std(background)</c>.</summary>
    /// <param name="background">Background in alphabet order; uniform when null.</param>
    public double Std(IReadOnlyList<double>? background = null)
        => MotifFinder.PwmStd(_matrix, Length, ResolveBackground(background, Alphabet.Length));

    internal static double[] ResolveBackground(IReadOnlyList<double>? background, int k)
    {
        if (background is null)
        {
            var uniform = new double[k];
            Array.Fill(uniform, 1.0 / k);
            return uniform;
        }

        return MotifFinder.NormalizeBackground(background, k,
            $"Background must have exactly {k} probabilities (one per alphabet symbol).");
    }

    /// <summary>
    /// Biopython <c>consensus</c> / <c>anticonsensus</c>: start at ∓∞ and take a symbol only on a strict improvement, so
    /// ties go to the earlier symbol. A column that never improves (all −∞ for the consensus) yields the first symbol.
    /// </summary>
    private string Extremum(bool maximum)
    {
        var sb = new StringBuilder(Length);
        int k = Alphabet.Length;
        for (int j = 0; j < Length; j++)
        {
            double best = maximum ? double.NegativeInfinity : double.PositiveInfinity;
            int bestRow = 0;
            for (int a = 0; a < k; a++)
            {
                double w = _matrix[a, j];
                if (maximum ? w > best : w < best)
                {
                    best = w;
                    bestRow = a;
                }
            }

            sb.Append(Alphabet[bestRow]);
        }

        return sb.ToString();
    }

    /// <summary>Biopython <c>max</c> / <c>min</c>: Σ over columns of Python's <c>max</c>/<c>min</c> of the column.</summary>
    private double ColumnExtremumSum(bool maximum) => MotifFinder.ColumnExtremumSum(_matrix, maximum);
}
