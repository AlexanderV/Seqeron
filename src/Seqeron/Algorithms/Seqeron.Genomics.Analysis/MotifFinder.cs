using System.Runtime.InteropServices;
using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Finds conserved motifs and patterns in DNA sequences.
/// Supports exact and degenerate motif searching, position weight matrices, and consensus sequences.
/// </summary>
public static partial class MotifFinder
{
    #region Exact Motif Finding

    /// <summary>
    /// Finds all occurrences of an exact motif in a sequence.
    /// Uses SuffixTree for O(m+k) pattern matching where m=motif length, k=occurrences.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="motif">Motif pattern to find.</param>
    /// <returns>Positions where the motif occurs.</returns>
    public static IEnumerable<int> FindExactMotif(DnaSequence sequence, string motif)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindExactMotifCore(sequence, motif);
    }

    private static IEnumerable<int> FindExactMotifCore(DnaSequence sequence, string motif)
    {
        if (string.IsNullOrEmpty(motif)) yield break;

        string motifUpper = motif.ToUpperInvariant();

        // Use SuffixTree for efficient pattern matching
        var positions = sequence.SuffixTree.FindAllOccurrences(motifUpper);
        foreach (int pos in positions.OrderBy(p => p))
        {
            yield return pos;
        }
    }

    #endregion

    #region Degenerate Motif Finding

    /// <summary>
    /// Validates that all characters in a motif pattern are valid IUPAC codes.
    /// Membership is decided by the canonical <see cref="IupacHelper.IsNucleotideCode"/>
    /// (the 15 IUPAC-IUB / NC-IUB 1984 DNA codes, = Biopython <c>IUPACData.ambiguous_dna_letters</c>).
    /// </summary>
    /// <param name="motif">The motif pattern to validate (expected already uppercased).</param>
    /// <exception cref="ArgumentException">Thrown when the pattern contains non-IUPAC characters.</exception>
    private static void ValidateIupacPattern(string motif)
    {
        for (int i = 0; i < motif.Length; i++)
        {
            if (!IupacHelper.IsNucleotideCode(motif[i]))
            {
                throw new ArgumentException(
                    $"Invalid IUPAC code '{motif[i]}' at position {i} in motif pattern. " +
                    "Valid codes: A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V.",
                    nameof(motif));
            }
        }
    }

    /// <summary>
    /// Finds all occurrences of a degenerate motif using IUPAC codes.
    /// </summary>
    /// <remarks>
    /// Pattern-degenerate matching (IUPAC-IUB 1970; NC-IUB 1984, Cornish-Bowden 1985): a window
    /// <c>S[i..i+m-1]</c> matches when every sequence base belongs to the base set of the motif code at
    /// that position (<see cref="IupacHelper.MatchesIupac"/>). The sequence side is literal — an ambiguity
    /// symbol in the sequence (e.g. <c>N</c>) matches no motif code — exactly as Biopython
    /// <c>Bio.SeqUtils.nt_search</c>. Positions are 0-based, ascending and overlapping.
    /// </remarks>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="motif">Motif pattern with IUPAC ambiguity codes.</param>
    /// <returns>Motif matches with positions.</returns>
    public static IEnumerable<MotifMatch> FindDegenerateMotif(DnaSequence sequence, string motif)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindDegenerateMotifCore(sequence, motif);
    }

    private static IEnumerable<MotifMatch> FindDegenerateMotifCore(DnaSequence sequence, string motif)
    {
        // Iterator: motif validation stays deferred to enumeration (unchanged contract).
        if (string.IsNullOrEmpty(motif)) yield break;

        // DnaSequence is already upper-case ACGT; the motif is validated even for an empty sequence.
        string motifUpper = motif.ToUpperInvariant();
        ValidateIupacPattern(motifUpper);
        foreach (var match in ScanDegenerate(sequence.Sequence, motifUpper, CancellationToken.None))
            yield return match;
    }

    /// <summary>
    /// Finds all occurrences of a degenerate motif with cancellation support.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="motif">Motif pattern with IUPAC ambiguity codes.</param>
    /// <param name="cancellationToken">Cancellation token for long-running operations.</param>
    /// <returns>Motif matches with positions.</returns>
    public static IEnumerable<MotifMatch> FindDegenerateMotif(
        DnaSequence sequence,
        string motif,
        CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindDegenerateMotifCore(sequence.Sequence, motif, cancellationToken);
    }

    /// <summary>
    /// Finds degenerate motif in a raw string with cancellation support.
    /// The sequence is upper-cased; any non-ACGT sequence character (including IUPAC ambiguity
    /// symbols, <c>U</c> and gaps) matches no motif code.
    /// </summary>
    public static IEnumerable<MotifMatch> FindDegenerateMotif(
        string sequence,
        string motif,
        CancellationToken cancellationToken)
    {
        return FindDegenerateMotifCore(sequence, motif, cancellationToken);
    }

    private static IEnumerable<MotifMatch> FindDegenerateMotifCore(
        string sequence,
        string motif,
        CancellationToken cancellationToken)
    {
        if (string.IsNullOrEmpty(sequence) || string.IsNullOrEmpty(motif))
            yield break;

        var seq = sequence.ToUpperInvariant();
        var motifUpper = motif.ToUpperInvariant();
        ValidateIupacPattern(motifUpper);
        foreach (var match in ScanDegenerate(seq, motifUpper, cancellationToken))
            yield return match;
    }

    /// <summary>
    /// Single O(n·m) window scan shared by all <c>FindDegenerateMotif</c> overloads; per-position
    /// membership is the canonical <see cref="IupacHelper.MatchesIupac"/> (no private code table).
    /// </summary>
    private static IEnumerable<MotifMatch> ScanDegenerate(
        string seq,
        string motifUpper,
        CancellationToken cancellationToken)
    {
        const int checkInterval = 1000;

        for (int i = 0; i <= seq.Length - motifUpper.Length; i++)
        {
            if (i % checkInterval == 0)
                cancellationToken.ThrowIfCancellationRequested();

            bool matches = true;
            for (int j = 0; j < motifUpper.Length && matches; j++)
                matches = IupacHelper.MatchesIupac(seq[i + j], motifUpper[j]);

            if (matches)
            {
                yield return new MotifMatch(
                    Position: i,
                    MatchedSequence: seq.Substring(i, motifUpper.Length),
                    Pattern: motifUpper,
                    Score: 1.0);
            }
        }
    }

    #endregion

    #region Position Weight Matrix

    /// <summary>Number of rows of a DNA PWM (A, C, G, T).</summary>
    private const int PwmAlphabetSize = 4;

    /// <summary>
    /// Creates a log-odds Position Weight Matrix (PWM) from aligned sequences against a
    /// uniform background (b = 0.25 for every base).
    /// <para>
    /// W[b,j] = log2( ((c[b,j] + p) / (N + 4p)) / 0.25 ), where c is the position frequency
    /// (count) matrix, N the number of sequences and p the pseudocount added to every cell.
    /// Identical to Biopython <c>motif.counts.normalize(pseudocounts=p).log_odds()</c>.
    /// </para>
    /// </summary>
    /// <param name="sequences">Aligned sequences of equal length over A/C/G/T (case-insensitive).</param>
    /// <param name="pseudocount">
    /// Pseudocount added to each of the four cells of every column (default: 0.25, i.e. one
    /// pseudo-observation per column spread uniformly). Must be finite and ≥ 0; 0 gives −∞ for
    /// unseen bases.
    /// </param>
    /// <returns>Position Weight Matrix.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">Empty collection, null element, unequal lengths or non-ACGT character.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pseudocount"/> is negative, NaN or infinite.</exception>
    public static PositionWeightMatrix CreatePwm(IEnumerable<string> sequences, double pseudocount = 0.25)
        => CreatePwm(sequences, pseudocount, UniformBackground);

    private static readonly double[] UniformBackground = { 0.25, 0.25, 0.25, 0.25 };

    /// <summary>
    /// Creates a log-odds Position Weight Matrix against an arbitrary background distribution
    /// (Wasserman &amp; Sandelin 2004, Nat Rev Genet 5:276; Biopython
    /// <c>counts.normalize(pseudocounts=p).log_odds(background=...)</c>):
    /// W[b,j] = log2( ((c[b,j] + p) / (N + 4p)) / q[b] ), with q the background normalised to sum 1.
    /// </summary>
    /// <param name="sequences">Aligned sequences of equal length over A/C/G/T (case-insensitive).</param>
    /// <param name="pseudocount">Pseudocount added to each cell (finite, ≥ 0).</param>
    /// <param name="background">
    /// Background probabilities in the order A, C, G, T (4 finite, strictly positive values;
    /// normalised to sum 1 as in Biopython).
    /// </param>
    /// <returns>Position Weight Matrix.</returns>
    public static PositionWeightMatrix CreatePwm(
        IEnumerable<string> sequences,
        double pseudocount,
        IReadOnlyList<double> background)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(background);
        if (!double.IsFinite(pseudocount) || pseudocount < 0)
            throw new ArgumentOutOfRangeException(nameof(pseudocount), pseudocount,
                "Pseudocount must be finite and non-negative.");
        double[] bg = NormalizeBackground(background);

        double[,] countMatrix = SequenceCountMatrix(sequences);
        int length = countMatrix.GetLength(1);

        var pseudocounts = new[] { pseudocount, pseudocount, pseudocount, pseudocount };
        return new PositionWeightMatrix(
            LogOddsFromCounts(countMatrix, pseudocounts, PwmAlphabetSize * pseudocount, bg), length);
    }

    /// <summary>
    /// Creates a log-odds PWM from aligned sequences with per-base pseudocounts — Biopython
    /// <c>motifs.create(seqs).counts.normalize(pseudocounts={'A':pA,'C':pC,'G':pG,'T':pT}).log_odds(background)</c>:
    /// W[b,j] = log2( ((c[b,j] + p[b]) / (N + Σ p)) / q[b] ). The count matrix is the one of
    /// <see cref="CreatePwm(IEnumerable{string}, double, IReadOnlyList{double})"/>; the log-odds are those of
    /// <see cref="PositionWeightMatrix.FromCounts(double[,], IReadOnlyList{double}, IReadOnlyList{double}?)"/>.
    /// </summary>
    /// <remarks>
    /// JASPAR / Wasserman &amp; Sandelin (2004) pseudocounts √N · q[b] for any background are obtained with
    /// <see cref="CreatePwmWithJasparPseudocounts"/> (Biopython <c>Bio.motifs.jaspar.calculate_pseudocounts</c>).
    /// </remarks>
    /// <param name="sequences">Aligned sequences of equal length over A/C/G/T (case-insensitive).</param>
    /// <param name="pseudocounts">Pseudocounts for A, C, G, T (finite, ≥ 0).</param>
    /// <param name="background">Background (A, C, G, T), normalised to sum 1; uniform when null.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> or <paramref name="pseudocounts"/> is null.</exception>
    /// <exception cref="ArgumentException">Empty collection, null element, unequal lengths, non-ACGT character, not 4 pseudocounts.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Negative / non-finite pseudocount or invalid background value.</exception>
    public static PositionWeightMatrix CreatePwm(
        IEnumerable<string> sequences,
        IReadOnlyList<double> pseudocounts,
        IReadOnlyList<double>? background = null)
    {
        ArgumentNullException.ThrowIfNull(pseudocounts);
        return PwmFromCounts(SequenceCountMatrix(sequences), pseudocounts, background);
    }

    /// <summary>
    /// Creates a log-odds PWM from aligned sequences with the JASPAR pseudocounts p[b] = √N · q[b]
    /// (<see cref="JasparPseudocounts"/>, Biopython <c>Bio.motifs.jaspar.calculate_pseudocounts</c>; Wasserman &amp;
    /// Sandelin 2004) against the same background q — Biopython
    /// <c>m.counts.normalize(pseudocounts=jaspar.calculate_pseudocounts(m)).log_odds(background)</c> with
    /// <c>m.background = background</c>.
    /// </summary>
    /// <param name="sequences">Aligned sequences of equal length over A/C/G/T (case-insensitive).</param>
    /// <param name="background">Background (A, C, G, T), normalised to sum 1; uniform when null.</param>
    public static PositionWeightMatrix CreatePwmWithJasparPseudocounts(
        IEnumerable<string> sequences,
        IReadOnlyList<double>? background = null)
    {
        double[,] counts = SequenceCountMatrix(sequences);
        return PwmFromCounts(counts, JasparPseudocounts(counts, background), background);
    }

    /// <summary>The <see cref="BuildCountMatrix"/> of a non-empty alignment as doubles (shared by the per-base pseudocount overloads).</summary>
    private static double[,] SequenceCountMatrix(IEnumerable<string> sequences)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        int[,] counts = BuildCountMatrix(sequences, nameof(sequences), out int count);
        if (count == 0)
            throw new ArgumentException("At least one sequence is required.", nameof(sequences));

        int length = counts.GetLength(1);
        var countMatrix = new double[PwmAlphabetSize, length];
        for (int i = 0; i < length; i++)
            for (int b = 0; b < PwmAlphabetSize; b++)
                countMatrix[b, i] = counts[b, i];
        return countMatrix;
    }

    /// <summary>
    /// Shared log-odds kernel of <see cref="CreatePwm(IEnumerable{string}, double, IReadOnlyList{double})"/> and
    /// <see cref="PositionWeightMatrix.FromCounts(double[,], IReadOnlyList{double}, IReadOnlyList{double}?)"/>
    /// (Biopython <c>counts.normalize(pseudocounts).log_odds(background)</c>):
    /// W[b,j] = log2( ((c[b,j] + p[b]) / (Σ_b c[b,j] + Σ_b p[b])) / q[b] ).
    /// </summary>
    /// <param name="counts">K × L count matrix (K = 4, rows A, C, G, T, for DNA; K = |alphabet| for
    /// <see cref="AlphabetPositionWeightMatrix"/>), validated by the caller.</param>
    /// <param name="pseudocounts">Per-row pseudocounts.</param>
    /// <param name="pseudocountSum">Σ p[b] (passed so the scalar path keeps its exact N + 4p total).</param>
    /// <param name="background">Normalised background (one value per row).</param>
    internal static double[,] LogOddsFromCounts(double[,] counts, double[] pseudocounts, double pseudocountSum, double[] background)
    {
        // Row count taken from the matrix: 4 for DNA, |alphabet| for AlphabetPositionWeightMatrix (one shared kernel).
        int rows = counts.GetLength(0);
        int length = counts.GetLength(1);
        var matrix = new double[rows, length];
        for (int i = 0; i < length; i++)
        {
            double columnCount = 0;
            for (int b = 0; b < rows; b++)
                columnCount += counts[b, i];
            double total = columnCount + pseudocountSum;

            // Pseudocount-smoothed probabilities → log2 odds against the background.
            for (int b = 0; b < rows; b++)
            {
                double freq = (counts[b, i] + pseudocounts[b]) / total;
                matrix[b, i] = Math.Log2(freq / background[b]);
            }
        }

        return matrix;
    }

    internal static double[] NormalizeBackground(IReadOnlyList<double> background)
        => NormalizeBackground(background, PwmAlphabetSize, "Background must have exactly 4 probabilities (A, C, G, T).");

    /// <summary>
    /// Background of <paramref name="size"/> finite, strictly positive values normalised to sum 1 (Biopython
    /// <c>log_odds</c> / <c>mean</c> / <c>std</c> divide by the total); shared by the DNA and the generic-alphabet PWM.
    /// </summary>
    internal static double[] NormalizeBackground(IReadOnlyList<double> background, int size, string countMessage)
    {
        if (background.Count != size)
            throw new ArgumentException(countMessage, nameof(background));

        double sum = 0;
        for (int b = 0; b < size; b++)
        {
            double v = background[b];
            if (!double.IsFinite(v) || v <= 0)
                throw new ArgumentOutOfRangeException(nameof(background), v,
                    "Background probabilities must be finite and strictly positive.");
            sum += v;
        }

        var result = new double[size];
        for (int b = 0; b < size; b++)
            result[b] = background[b] / sum;
        return result;
    }

    /// <summary>Row index of a PWM for an upper-case base (A=0, C=1, G=2, T=3), or −1.</summary>
    internal static int AcgtIndex(char c) => c switch
    {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => -1
    };

    /// <summary>Upper-case bases in PWM row order (inverse of <see cref="AcgtIndex"/>); shared with <see cref="PositionWeightMatrix"/>.</summary>
    internal static readonly char[] AcgtBases = { 'A', 'C', 'G', 'T' };

    /// <summary>
    /// Materialises an aligned-sequence collection and enforces the shared alignment contract of
    /// <see cref="BuildCountMatrix"/> and <see cref="GenerateConsensus(IEnumerable{string})"/>:
    /// no null element (checked while enumerating) and all rows of equal length.
    /// </summary>
    /// <exception cref="ArgumentException">Null element or unequal lengths.</exception>
    private static List<string> MaterializeAligned(IEnumerable<string> sequences, string paramName)
    {
        List<string> seqList = MaterializeRows(sequences, paramName);

        for (int s = 1; s < seqList.Count; s++)
        {
            if (seqList[s].Length != seqList[0].Length)
                throw new ArgumentException("All sequences must have the same length.", paramName);
        }

        return seqList;
    }

    /// <summary>Materialises the rows of an alignment, rejecting a null element (checked while enumerating).</summary>
    /// <exception cref="ArgumentException">Null element.</exception>
    private static List<string> MaterializeRows(IEnumerable<string> sequences, string paramName)
    {
        var rows = new List<string>();
        foreach (var s in sequences)
        {
            if (s is null)
                throw new ArgumentException("Sequences cannot contain null elements.", paramName);
            rows.Add(s);
        }

        return rows;
    }

    /// <summary>
    /// Builds the 4 × L position frequency (count) matrix — Rosalind CONS "profile matrix",
    /// Biopython <c>motif.counts</c> — of equal-length A/C/G/T sequences (case-insensitive),
    /// rows A, C, G, T. Shared by <see cref="CreatePwm(IEnumerable{string}, double, IReadOnlyList{double})"/>
    /// and <see cref="CreateConsensusFromAlignment"/>. An empty collection yields a 4 × 0 matrix
    /// with <paramref name="sequenceCount"/> = 0 (callers decide whether that is an error).
    /// </summary>
    /// <exception cref="ArgumentException">Null element, unequal lengths or non-ACGT character.</exception>
    private static int[,] BuildCountMatrix(IEnumerable<string> sequences, string paramName, out int sequenceCount)
    {
        List<string> seqList = MaterializeAligned(sequences, paramName);

        sequenceCount = seqList.Count;
        if (sequenceCount == 0)
            return new int[PwmAlphabetSize, 0];

        int length = seqList[0].Length;
        var counts = new int[PwmAlphabetSize, length];
        for (int s = 0; s < seqList.Count; s++)
        {
            var seq = seqList[s];
            for (int i = 0; i < length; i++)
            {
                char c = char.ToUpperInvariant(seq[i]);
                int baseIndex = AcgtIndex(c);
                if (baseIndex < 0)
                    throw new ArgumentException(
                        $"Invalid character '{c}' at position {i} in sequence {s}. " +
                        "Only A, C, G, T are valid nucleotide characters.",
                        paramName);

                counts[baseIndex, i]++;
            }
        }

        return counts;
    }

    /// <summary>
    /// Scans the forward strand of a sequence with a PWM and returns every window whose score
    /// (sum of per-position log-odds) is ≥ <paramref name="threshold"/>, in ascending position
    /// order (Biopython <c>pssm.search(seq, threshold, both=False)</c>). Windows containing a
    /// non-ACGT symbol are skipped. For both strands use
    /// <see cref="ScanWithPwmBothStrands(DnaSequence, PositionWeightMatrix, double)"/> (Biopython
    /// <c>both=True</c>); for every window score use <see cref="CalculatePwmScores(string, PositionWeightMatrix)"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence to scan.</param>
    /// <param name="pwm">Position Weight Matrix.</param>
    /// <param name="threshold">Minimum score threshold (inclusive).</param>
    /// <returns>Matches with scores.</returns>
    public static IEnumerable<MotifMatch> ScanWithPwm(
        DnaSequence sequence,
        PositionWeightMatrix pwm,
        double threshold = 0.0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(pwm);
        return ScanWithPwmCore(sequence, pwm, threshold);
    }

    private static IEnumerable<MotifMatch> ScanWithPwmCore(DnaSequence sequence, PositionWeightMatrix pwm, double threshold)
    {
        string seq = sequence.Sequence;
        int motifLen = pwm.Length;

        for (int i = 0; i <= seq.Length - motifLen; i++)
        {
            // NaN (window with a non-ACGT symbol) never satisfies score >= threshold.
            double score = ScorePwmWindow(seq, i, pwm);
            if (score >= threshold)
            {
                yield return new MotifMatch(
                    Position: i,
                    MatchedSequence: seq.Substring(i, motifLen),
                    Pattern: pwm.Consensus,
                    Score: score);
            }
        }
    }

    #endregion

    #region Consensus Sequence

    /// <summary>
    /// Generates an IUPAC-degenerate consensus from aligned sequences with a per-base frequency
    /// threshold: at each column the bases occurring in strictly more than 25 % of the sequences
    /// are combined into their NC-IUB 1984 symbol (Cornish-Bowden, NAR 13(9):3021). When no base
    /// passes the threshold, the column emits the IUPAC symbol for the set of bases sharing the
    /// maximum count (DECIPHER <c>ConsensusSequence</c>: "degeneracy codes are always used in
    /// cases where multiple characters are equally abundant"), so four equally frequent bases
    /// give <c>N</c>; a column with no A/C/G/T at all gives <c>N</c> (any base).
    /// </summary>
    /// <remarks>
    /// The 25 % per-base cut is this overload's fixed default (see <c>IupacInclusionThreshold</c>);
    /// any other cut is available through <see cref="GenerateConsensus(IEnumerable{string}, double)"/>.
    /// It is not the Cavener (1987) / TRANSFAC / Biopython <c>degenerate_consensus</c> rule (use
    /// <see cref="GenerateCavenerConsensus"/>) nor DECIPHER's <c>ConsensusSequence</c> (use
    /// <see cref="GenerateDecipherConsensus"/>).
    /// Characters other than A/C/G/T (gaps, N, IUPAC codes) are not counted but still count
    /// towards the number of sequences <c>n</c>. Case-insensitive.
    /// </remarks>
    /// <param name="sequences">Aligned sequences of equal length.</param>
    /// <returns>Consensus sequence over the 15 IUPAC symbols; "" for an empty collection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentException">A null element, or sequences of unequal length.</exception>
    public static string GenerateConsensus(IEnumerable<string> sequences)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return GenerateConsensusCore(sequences, IupacInclusionThreshold);
    }

    /// <summary>Shared body of both <c>GenerateConsensus</c> overloads (threshold already validated).</summary>
    private static string GenerateConsensusCore(IEnumerable<string> sequences, double inclusionThreshold)
    {
        List<string> seqList = MaterializeAligned(sequences, nameof(sequences));
        if (seqList.Count == 0) return "";

        int length = seqList[0].Length;
        var consensus = new StringBuilder(length);
        var counts = new int[PwmAlphabetSize];

        for (int i = 0; i < length; i++)
        {
            Array.Clear(counts);
            foreach (var seq in seqList)
            {
                int baseIndex = AcgtIndex(char.ToUpperInvariant(seq[i]));
                if (baseIndex >= 0)
                    counts[baseIndex]++;
            }

            consensus.Append(GetIupacCode(counts, seqList.Count, inclusionThreshold));
        }

        return consensus.ToString();
    }

    /// <summary>
    /// Generates the degenerate consensus by the Cavener (1987) rules (Nucleic Acids Res.
    /// 15(4):1353–1361), as used by TRANSFAC and Biopython <c>Bio.motifs</c>
    /// <c>degenerate_consensus</c>. Per column, with base counts sorted in decreasing order
    /// c1 ≥ c2 ≥ c3 ≥ c4 (ties in A, C, G, T order):
    /// <list type="number">
    /// <item>a single base if c1 &gt; c2 + c3 + c4 (more than 50 %) and c1 &gt; 2·c2;</item>
    /// <item>otherwise the two-base IUPAC code of the top two bases if c1 + c2 &gt; 75 % of the column;</item>
    /// <item>otherwise the three-base code of the top three bases if the fourth base is absent (c4 = 0);</item>
    /// <item>otherwise <c>N</c>.</item>
    /// </list>
    /// The set → symbol mapping is NC-IUB 1984.
    /// </summary>
    /// <param name="alignedSequences">Aligned DNA sequences of equal length over {A, C, G, T} (case-insensitive).</param>
    /// <returns>The degenerate consensus; "" for an empty collection.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="alignedSequences"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// A null element, sequences of unequal length, or a non-ACGT character (gaps are not accepted).
    /// </exception>
    public static string GenerateCavenerConsensus(IEnumerable<string> alignedSequences)
    {
        ArgumentNullException.ThrowIfNull(alignedSequences);

        int[,] counts = BuildCountMatrix(alignedSequences, nameof(alignedSequences), out _);
        int length = counts.GetLength(1);
        var consensus = new StringBuilder(length);
        var order = new int[PwmAlphabetSize];
        var c = new int[PwmAlphabetSize];

        for (int col = 0; col < length; col++)
        {
            // Stable sort of base indices by decreasing count (ties keep A, C, G, T order),
            // identical to Biopython's sorted(..., key=count, reverse=True).
            for (int b = 0; b < PwmAlphabetSize; b++)
                order[b] = b;
            for (int x = 1; x < PwmAlphabetSize; x++)
            {
                int cur = order[x];
                int y = x - 1;
                while (y >= 0 && counts[order[y], col] < counts[cur, col])
                {
                    order[y + 1] = order[y];
                    y--;
                }
                order[y + 1] = cur;
            }

            for (int k = 0; k < PwmAlphabetSize; k++)
                c[k] = counts[order[k], col];
            int total = c[0] + c[1] + c[2] + c[3];

            int take;
            if (c[0] > c[1] + c[2] + c[3] && c[0] > 2 * c[1])
                take = 1;
            else if (4 * (c[0] + c[1]) > 3 * total)
                take = 2;
            else if (c[3] == 0)
                take = 3;
            else
                take = 4;

            var set = new char[take];
            for (int k = 0; k < take; k++)
                set[k] = AcgtBases[order[k]];
            consensus.Append(Core.IupacDnaSequence.GetIupacCode(set));
        }

        return consensus.ToString();
    }

    /// <summary>
    /// Creates a consensus sequence from a multiple alignment by selecting, at each column,
    /// the most frequent nucleotide (the symbol with the maximum count in that column of the
    /// profile matrix). This is the classical "most common symbol per position" consensus
    /// (Rosalind CONS; Wikipedia "Consensus sequence"), not the IUPAC-degenerate variant
    /// produced by <see cref="GenerateConsensus(IEnumerable{string})"/>.
    /// </summary>
    /// <param name="alignedSequences">
    /// Aligned DNA sequences of equal length over the alphabet {A, C, G, T} (case-insensitive).
    /// </param>
    /// <returns>
    /// The consensus string. Each position is the most frequent base in that column; ties are
    /// broken alphabetically (A &lt; C &lt; G &lt; T) for determinism. Returns an empty string for
    /// an empty collection.
    /// </returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="alignedSequences"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the collection contains a null element, the sequences are not all of equal
    /// length, or a sequence contains a non-ACGT character (gaps are not accepted).
    /// </exception>
    public static string CreateConsensusFromAlignment(IEnumerable<string> alignedSequences)
    {
        ArgumentNullException.ThrowIfNull(alignedSequences);

        int[,] counts = BuildCountMatrix(alignedSequences, nameof(alignedSequences), out _);
        int length = counts.GetLength(1);
        var consensus = new StringBuilder(length);

        for (int col = 0; col < length; col++)
        {
            // Profile-column maximum. Rows are in alphabetical order (A, C, G, T) and only a
            // strictly greater count displaces the incumbent, so a tie resolves to the
            // alphabetically-earliest base — the same scan as Biopython
            // Bio.motifs GenericPositionMatrix.consensus ("if count > maximum"), and the
            // LANL/Geneious "residue letter occurring earlier in the alphabet" tie-break.
            int bestIndex = 0;
            for (int b = 1; b < PwmAlphabetSize; b++)
            {
                if (counts[b, col] > counts[bestIndex, col])
                    bestIndex = b;
            }

            consensus.Append(AcgtBases[bestIndex]);
        }

        return consensus.ToString();
    }

    /// <summary>
    /// Minimum column frequency (as a fraction of the number of aligned sequences) a base must
    /// strictly exceed to be included in the position's IUPAC degeneracy code in
    /// <see cref="GenerateConsensus(IEnumerable{string})"/>. The set→symbol mapping follows
    /// NC-IUB 1984 (Cornish-Bowden, NAR 13(9):3021). The per-base 0.25 cut with a strict '&gt;'
    /// boundary is the default of the parameterless overload: it belongs to the threshold-consensus
    /// family but is NOT Bioconductor DECIPHER's rule (implemented as
    /// <see cref="GenerateDecipherConsensus"/>) nor the Cavener 1987 rule (implemented as
    /// <see cref="GenerateCavenerConsensus"/>). Any other cut is available through
    /// <see cref="GenerateConsensus(IEnumerable{string}, double)"/> (MCP <c>generate_consensus</c>
    /// parameter <c>inclusionThreshold</c>). See docs/Evidence/MOTIF-GENERATE-001-Evidence.md.
    /// </summary>
    private const double IupacInclusionThreshold = 0.25;

    private static char GetIupacCode(int[] counts, int total, double inclusionThreshold)
    {
        double threshold = total * inclusionThreshold; // base count must be strictly > threshold

        var present = new List<char>(PwmAlphabetSize);
        for (int b = 0; b < PwmAlphabetSize; b++)
        {
            if (counts[b] > threshold)
                present.Add(AcgtBases[b]);
        }

        if (present.Count == 0)
        {
            // No base passes: encode every base sharing the maximum count (DECIPHER
            // ConsensusSequence — "degeneracy codes are always used in cases where multiple
            // characters are equally abundant"); a column without any A/C/G/T is unknown → N.
            int max = counts.Max();
            if (max == 0)
                return 'N';
            for (int b = 0; b < PwmAlphabetSize; b++)
            {
                if (counts[b] == max)
                    present.Add(AcgtBases[b]);
            }
        }

        // Base set → NC-IUB symbol via the canonical inverse map; present is a non-empty
        // subset of {A,C,G,T}, and the full set {A,C,G,T} maps to N.
        return Core.IupacDnaSequence.GetIupacCode(present);
    }

    #endregion

    #region Motif Discovery

    /// <summary>
    /// Discovers overrepresented k-mers that may represent motifs.
    /// </summary>
    /// <remarks>
    /// Overrepresentation is measured by the observed-over-expected (O/E) ratio. Under the
    /// zero-order (i.i.d. uniform) background model in which each of the four nucleotides is
    /// equally likely (probability 1/4), the expected number of occurrences of any specific
    /// k-mer in a string of length N is E = (N − k + 1) / 4^k (Compeau &amp; Pevzner,
    /// <i>Bioinformatics Algorithms: An Active Learning Approach</i>; the (N − k + 1) factor
    /// is the number of length-k windows and 4^k is the number of distinct k-mers). The
    /// <see cref="DiscoveredMotif.Enrichment"/> field is the observed count divided by E, so a
    /// value &gt; 1 means the k-mer occurs more often than chance predicts. This is the
    /// "equiprobable" model of RSAT <c>oligo-analysis</c> (van Helden et al. 1998); use
    /// <see cref="DiscoverMotifs(DnaSequence, int, int, IReadOnlyList{double})"/> for a
    /// non-uniform (Bernoulli) background. Occurrences overlap (every window is counted).
    /// Results are yielded in order of each k-mer's first occurrence. For RSAT binomial significance
    /// (occ_P / occ_E / occ_sig), Markov backgrounds and both-strand counting use
    /// <see cref="DiscoverMotifs(DnaSequence, int, int, OligoBackgroundModel, OligoStrandMode, bool)"/>.
    /// </remarks>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="k">K-mer length (default: 6).</param>
    /// <param name="minCount">Minimum occurrence count (default: 2).</param>
    /// <returns>Overrepresented k-mers with their counts, positions, and O/E enrichment.</returns>
    public static IEnumerable<DiscoveredMotif> DiscoverMotifs(
        DnaSequence sequence,
        int k = 6,
        int minCount = 2)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        return DiscoverMotifsCore(sequence, k, minCount, UniformBackground);
    }

    /// <summary>
    /// Discovers overrepresented k-mers against a Bernoulli (independent, non-uniform nucleotide)
    /// background model.
    /// </summary>
    /// <remarks>
    /// RSAT <c>oligo-analysis</c> Bernoulli model (van Helden, André &amp; Collado-Vides 1998,
    /// J Mol Biol 281:827): the expected frequency of a word w = w₁…w_k is the product of its
    /// residue probabilities, p(w) = ∏ q[w_i], and its expected number of occurrences is
    /// E = p(w) · (N − k + 1) (the number of overlapping length-k windows). Enrichment = Count / E
    /// (RSAT "ratio" column). With q = (¼, ¼, ¼, ¼) this equals
    /// <see cref="DiscoverMotifs(DnaSequence, int, int)"/> exactly. The ratio is evaluated with
    /// exact power-of-two rescaling, so it stays finite whenever the true value is representable
    /// (the product p(w) itself underflows for long k).
    /// </remarks>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="k">K-mer length (≥ 1).</param>
    /// <param name="minCount">Minimum occurrence count for a k-mer to be returned.</param>
    /// <param name="background">
    /// Background residue probabilities in the order A, C, G, T (4 finite, strictly positive
    /// values; normalised to sum 1, as for <see cref="CreatePwm(IEnumerable{string}, double, IReadOnlyList{double})"/>).
    /// </param>
    /// <returns>Overrepresented k-mers with their counts, positions, and O/E enrichment.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="background"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> &lt; 1, or a background value is not finite and positive.</exception>
    /// <exception cref="ArgumentException"><paramref name="background"/> does not have exactly 4 values.</exception>
    public static IEnumerable<DiscoveredMotif> DiscoverMotifs(
        DnaSequence sequence,
        int k,
        int minCount,
        IReadOnlyList<double> background)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(background);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        return DiscoverMotifsCore(sequence, k, minCount, NormalizeBackground(background));
    }

    private static IEnumerable<DiscoveredMotif> DiscoverMotifsCore(
        DnaSequence sequence, int k, int minCount, double[] background)
    {
        string seq = sequence.Sequence;
        Dictionary<string, List<int>> kmerPositions = CollectKmerPositions(seq, k, minCount);

        // Number of overlapping length-k windows, N − k + 1 (≥ 1 whenever a k-mer was counted).
        double windowCount = seq.Length - k + 1.0;

        foreach (var (kmer, positions) in kmerPositions)
        {
            yield return new DiscoveredMotif(
                Sequence: kmer,
                Count: positions.Count,
                Positions: positions.AsReadOnly(),
                Enrichment: ObservedOverExpected(positions.Count, windowCount, kmer, background));
        }
    }

    /// <summary>
    /// Overlapping k-mer counts via the canonical <see cref="SequenceExtensions.CountKmersSpan"/>,
    /// then 0-based start positions (ascending) collected only for k-mers with count ≥
    /// <paramref name="minCount"/>. Insertion order = order of first occurrence.
    /// </summary>
    private static Dictionary<string, List<int>> CollectKmerPositions(string seq, int k, int minCount)
    {
        Dictionary<string, int> counts = seq.AsSpan().CountKmersSpan(k);

        var positions = new Dictionary<string, List<int>>();
        foreach (var (kmer, count) in counts)
        {
            if (count >= minCount)
                positions.Add(kmer, new List<int>(count));
        }

        if (positions.Count == 0)
            return positions;

        var lookup = positions.GetAlternateLookup<ReadOnlySpan<char>>();
        for (int i = 0; i <= seq.Length - k; i++)
        {
            if (lookup.TryGetValue(seq.AsSpan(i, k), out var list))
                list.Add(i);
        }

        return positions;
    }

    // 2^-500: threshold/scale for renormalising the running product of residue probabilities.
    private const int ProbabilityRescaleExponent = 500;
    private static readonly double ProbabilityRescaleThreshold = Math.ScaleB(1.0, -ProbabilityRescaleExponent);

    /// <summary>
    /// O/E ratio Count / (W · ∏ q[w_i]) for an upper-case ACGT word. The product is kept as
    /// mantissa · 2^exponent (exact power-of-two rescaling), so a product that would underflow
    /// (e.g. 4^-k for k ≥ 512) does not turn a finite ratio into +∞. For q = ¼ every step is
    /// exact, giving the correctly rounded Count · 4^k / W.
    /// </summary>
    private static double ObservedOverExpected(int count, double windowCount, string kmer, double[] background)
    {
        double mantissa = 1.0;
        int exponent = 0;
        foreach (char c in kmer)
        {
            mantissa *= background[AcgtIndex(c)];
            if (mantissa < ProbabilityRescaleThreshold)
            {
                mantissa = Math.ScaleB(mantissa, ProbabilityRescaleExponent);
                exponent -= ProbabilityRescaleExponent;
            }
        }

        return Math.ScaleB(count / (windowCount * mantissa), -exponent);
    }

    // Default oligonucleotide (word) length for shared-motif enumeration. RSAT oligo-analysis
    // permits any oligo length in [1,8]; 6 is a common default and sits inside that range.
    // Source: RSAT oligo-analysis manual (https://rsat.eead.csic.es/plants/help.oligo-analysis.html).
    private const int DefaultSharedMotifLength = 6;

    // Default quorum: a word must occur in at least this many input sequences ("matching
    // sequences") to be reported. 2 is the minimum meaningful "shared" threshold.
    // Source: word-enumeration quorum (Das & Dai 2007, BMC Bioinformatics 8(S7):S21).
    private const int DefaultMinMatchingSequences = 2;

    /// <summary>
    /// Finds fixed-length words (oligonucleotides) shared across multiple DNA sequences,
    /// using the "matching sequences" quorum of the van Helden / RSAT oligo-analysis method:
    /// each length-<paramref name="k"/> word is scored by the number of input sequences that
    /// contain at least one (exact) occurrence of it, and words whose matching-sequence count
    /// is at least <paramref name="minSequences"/> are reported.
    /// A word repeated several times within one sequence still contributes 1 to its
    /// matching-sequence count. Matching is exact (no degenerate/substituted matches).
    /// </summary>
    /// <remarks>
    /// Equals the RSAT <c>oligo-analysis</c> <c>mseq</c> column for <c>-1str -ovlp</c> (single strand;
    /// the reverse complement is not merged; RSAT ms_P / ms_E / ms_sig, background models and <c>-2str</c> are provided by
    /// <see cref="FindSharedMotifs(IEnumerable{DnaSequence}, int, int, OligoBackgroundModel, OligoStrandMode)"/>). Words are enumerated with the canonical
    /// <see cref="SequenceExtensions.CountKmersSpan"/>. Results are yielded in order of each word's
    /// first occurrence (lowest sequence index, then lowest position); <see cref="SharedMotif.SequenceIndices"/>
    /// is strictly ascending. This is not the longest-common-substring (Rosalind LCSM) problem:
    /// k is fixed and membership is a quorum. Enumeration is deferred: a null element throws
    /// <see cref="ArgumentException"/> when the result is enumerated.
    /// </remarks>
    /// <param name="sequences">Collection of DNA sequences (each scanned for its distinct words).</param>
    /// <param name="k">Word (oligonucleotide) length; must be ≥ 1. Default: 6.</param>
    /// <param name="minSequences">Quorum: minimum number of distinct sequences a word must occur in. Must be ≥ 1. Default: 2.</param>
    /// <returns>Shared words with their distinct sequence indices and prevalence (matching sequences / total sequences).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">When <paramref name="k"/> &lt; 1 or <paramref name="minSequences"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">When an element of <paramref name="sequences"/> is null (thrown on enumeration).</exception>
    public static IEnumerable<SharedMotif> FindSharedMotifs(
        IEnumerable<DnaSequence> sequences,
        int k = DefaultSharedMotifLength,
        int minSequences = DefaultMinMatchingSequences)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(minSequences, 1);
        return FindSharedMotifsCore(sequences, k, minSequences);
    }

    private static IEnumerable<SharedMotif> FindSharedMotifsCore(IEnumerable<DnaSequence> sequences, int k, int minSequences)
    {
        var seqList = sequences.ToList();

        // word → ascending indices of the input sequences containing it (RSAT "matching sequences").
        var matchingSequences = new Dictionary<string, List<int>>();

        for (int seqIdx = 0; seqIdx < seqList.Count; seqIdx++)
        {
            DnaSequence dna = seqList[seqIdx]
                ?? throw new ArgumentException($"Sequence at index {seqIdx} is null.", nameof(sequences));

            // Distinct words of this sequence (canonical overlapping k-mer counter); each word
            // contributes this sequence once regardless of its multiplicity. Dictionary keys keep
            // first-occurrence order, so words are registered in (sequence, position) order.
            foreach (string kmer in dna.Sequence.AsSpan().CountKmersSpan(k).Keys)
            {
                ref List<int>? indices = ref CollectionsMarshal.GetValueRefOrAddDefault(matchingSequences, kmer, out _);
                (indices ??= new List<int>()).Add(seqIdx);
            }
        }

        foreach (var (kmer, seqIndices) in matchingSequences)
        {
            if (seqIndices.Count >= minSequences)
            {
                yield return new SharedMotif(
                    Sequence: kmer,
                    SequenceIndices: seqIndices.AsReadOnly(),
                    Prevalence: (double)seqIndices.Count / seqList.Count);
            }
        }
    }

    #endregion

    #region Regulatory Motif Patterns

    /// <summary>
    /// Known regulatory motif patterns: published 5'→3' consensus strings, written in IUPAC
    /// nucleotide code (NC-IUB 1985) exactly as the cited source states them. Degenerate positions
    /// (R, S, W, Y, N) are kept — they are matched by the canonical IUPAC scan
    /// (<see cref="FindDegenerateMotif(DnaSequence, string)"/>), never expanded or collapsed to a
    /// single representative string.
    /// </summary>
    public static class KnownMotifs
    {
        // Eukaryotic core-promoter elements (Bucher 1990 weight-matrix analysis of 502 Pol II promoters).
        /// <summary>TATA box core consensus: TATAAA (eukaryotic RNA Pol II core promoter; Bucher 1990).</summary>
        public const string TataBox = "TATAAA";

        /// <summary>CCAAT box core pentanucleotide: CCAAT (Bucher 1990; ~30 % of promoters).</summary>
        public const string CaatBox = "CCAAT";

        /// <summary>
        /// GC box (Sp1 site) core hexanucleotide: GGGCGG (Dynan &amp; Tjian 1983; Gidoni, Dynan &amp;
        /// Tjian 1984 — the six GGGCGG copies of the SV40 21-bp repeats).
        /// </summary>
        public const string GcBox = "GGGCGG";

        // Prokaryotic sigma-70 promoter hexamers per Harley &amp; Reynolds (1987) compilation.
        /// <summary>-10 (Pribnow) box consensus hexamer: TATAAT (Pribnow 1975; Harley &amp; Reynolds 1987).</summary>
        public const string MinusTenBox = "TATAAT";

        /// <summary>-35 box consensus hexamer: TTGACA (Harley &amp; Reynolds 1987).</summary>
        public const string MinusThirtyFiveBox = "TTGACA";

        /// <summary>
        /// Kozak vertebrate initiation consensus GCCGCC(A/G)CCATGG = <c>GCCGCCRCCATGG</c>
        /// (Kozak 1987, 699 vertebrate mRNAs; positions −9..+4, ATG = +1..+3, purine at −3).
        /// </summary>
        public const string Kozak = "GCCGCCRCCATGG";

        /// <summary>
        /// Shine-Dalgarno (bacterial RBS) consensus: AGGAGG, the complement of the 3'-terminal
        /// CCUCCU of E. coli 16S rRNA (Shine &amp; Dalgarno 1974).
        /// </summary>
        public const string ShineDalgarno = "AGGAGG";

        /// <summary>Poly(A) signal hexamer: AATAAA (Proudfoot &amp; Brownlee 1976).</summary>
        public const string PolyASignal = "AATAAA";

        /// <summary>E-box consensus (IUPAC): CANNTG (Massari &amp; Murre 2000).</summary>
        public const string EBox = "CANNTG";

        /// <summary>
        /// AP-1 site / TPA-response element (TRE) consensus TGA(C/G)TCA = <c>TGASTCA</c>
        /// (Lee, Mitchell &amp; Tjian 1987; Angel et al. 1987 — the collagenase TRE is TGAGTCA).
        /// The pattern is its own reverse complement.
        /// </summary>
        public const string Ap1 = "TGASTCA";

        /// <summary>
        /// NF-κB κB-site consensus <c>GGGRNWYYCC</c> (Gilmore 2006); includes the Ig κ enhancer
        /// site GGGACTTTCC of Sen &amp; Baltimore (1986).
        /// </summary>
        public const string NfKb = "GGGRNWYYCC";

        /// <summary>CREB CRE palindrome: TGACGTCA (Montminy et al. 1986).</summary>
        public const string Creb = "TGACGTCA";
    }

    /// <summary>
    /// The fixed element library scanned by <see cref="FindRegulatoryElements(DnaSequence)"/>, in report order.
    /// <c>OrientationIndependent</c> marks elements documented to act in either orientation:
    /// CCAAT box — "found in the forward or reverse orientation" (Mantovani 1998, NAR 26:1135, survey of 178
    /// NF-Y sites); GC box — Sp1 binds the SV40 GC boxes in both orientations, driving bidirectional
    /// transcription (Gidoni et al. 1985, Science 230:511); AP-1, NF-κB, E-box and CREB sites are enhancer
    /// elements, which act "in either orientation" (Banerji, Rusconi &amp; Schaffner 1981, Cell 27:299).
    /// TATA, −10/−35 boxes, Kozak, Shine–Dalgarno and the poly(A) signal act on their own strand.
    /// </summary>
    private static readonly (string Name, string Pattern, string Description, bool OrientationIndependent)[] RegulatoryLibrary =
    {
        ("TATA Box", KnownMotifs.TataBox, "Eukaryotic core promoter element", false),
        ("CAAT Box", KnownMotifs.CaatBox, "Promoter element", true),
        ("GC Box", KnownMotifs.GcBox, "Sp1 binding site", true),
        ("-10 Box", KnownMotifs.MinusTenBox, "Prokaryotic Pribnow box", false),
        ("-35 Box", KnownMotifs.MinusThirtyFiveBox, "Prokaryotic -35 promoter element", false),
        ("Kozak", KnownMotifs.Kozak, "Translation initiation", false),
        ("Shine-Dalgarno", KnownMotifs.ShineDalgarno, "Bacterial ribosome binding", false),
        ("Poly(A) Signal", KnownMotifs.PolyASignal, "Polyadenylation signal", false),
        ("E-box", KnownMotifs.EBox, "bHLH transcription factor binding", true),
        ("AP-1", KnownMotifs.Ap1, "AP-1 transcription factor binding", true),
        ("NF-κB", KnownMotifs.NfKb, "NF-κB binding site", true),
        ("CREB", KnownMotifs.Creb, "CREB transcription factor binding", true)
    };

    /// <summary>
    /// Names of the <see cref="FindRegulatoryElements(DnaSequence)"/> library elements that act in either
    /// orientation (sources on the library): CAAT Box, GC Box, E-box, AP-1, NF-κB, CREB.
    /// </summary>
    public static IReadOnlyList<string> OrientationIndependentRegulatoryElements { get; } =
        RegulatoryLibrary.Where(e => e.OrientationIndependent).Select(e => e.Name).ToArray();

    /// <summary>
    /// Scans the given strand of a DNA sequence for the <see cref="KnownMotifs"/> consensus library.
    /// </summary>
    /// <remarks>
    /// Each library pattern is matched with the canonical IUPAC degenerate scan
    /// (<see cref="FindDegenerateMotif(DnaSequence, string)"/> → <see cref="IupacHelper.MatchesIupac"/>):
    /// every 0-based start <c>0 ≤ i ≤ n−m</c> whose window lies in the IUPAC sets of the pattern is
    /// reported, overlapping occurrences included. Results are grouped by library entry (library order,
    /// see <see cref="KnownMotifs"/>) and ascending by position within an entry. Only the given strand is
    /// scanned; the AP-1, E-box and CREB patterns are their own reverse complements, so their hits cover
    /// both orientations. For reverse-orientation hits of the other orientation-independent elements use
    /// <see cref="FindRegulatoryElements(DnaSequence, bool)"/>.
    /// </remarks>
    /// <param name="sequence">DNA sequence to scan.</param>
    /// <returns>Found regulatory elements.</returns>
    public static IEnumerable<RegulatoryElement> FindRegulatoryElements(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindRegulatoryElementsCore(sequence);
    }

    private static IEnumerable<RegulatoryElement> FindRegulatoryElementsCore(DnaSequence sequence)
    {
        foreach (var (name, pattern, description, _) in RegulatoryLibrary)
        {
            foreach (var match in FindDegenerateMotif(sequence, pattern))
            {
                yield return new RegulatoryElement(
                    Name: name,
                    Position: match.Position,
                    Sequence: match.MatchedSequence,
                    Pattern: pattern,
                    Description: description);
            }
        }
    }

    /// <summary>
    /// Scans for the <see cref="KnownMotifs"/> library with strand annotation. With
    /// <paramref name="bothStrands"/> = false the result is exactly
    /// <see cref="FindRegulatoryElements(DnaSequence)"/> (all hits on '+'). With
    /// <paramref name="bothStrands"/> = true, the orientation-independent elements
    /// (<see cref="OrientationIndependentRegulatoryElements"/>) whose IUPAC pattern differs from its own
    /// reverse complement — CAAT Box, GC Box, NF-κB — are additionally matched on the minus strand.
    /// </summary>
    /// <remarks>
    /// The minus strand is scanned through the same canonical IUPAC path by matching the reverse
    /// complement of the pattern (canonical IUPAC complement,
    /// <see cref="DnaSequence.GetReverseComplementString"/>) on the forward sequence; this is identical to
    /// matching the pattern on the reverse-complement sequence (Biopython <c>nt_search</c> on
    /// <c>seq.reverse_complement()</c>) with coordinates mapped back. Self-reverse-complementary patterns
    /// (AP-1 TGASTCA, E-box CANNTG, CREB TGACGTCA) are not rescanned: every minus-strand occurrence is
    /// the same window as a plus-strand one. Strand-specific elements (TATA, −10/−35, Kozak,
    /// Shine–Dalgarno, poly(A)) are reported on the given strand only; for a minus-strand gene scan
    /// <see cref="DnaSequence.ReverseComplement"/>.
    /// Order: library order; within an entry ascending forward position, '+' before '-' at equal positions.
    /// <see cref="StrandedRegulatoryElement.Position"/> is the 0-based forward start of the window and
    /// <see cref="StrandedRegulatoryElement.Sequence"/> is the site read 5'→3' on its own strand.
    /// </remarks>
    /// <param name="sequence">DNA sequence to scan.</param>
    /// <param name="bothStrands">Also report minus-strand hits of orientation-independent elements.</param>
    /// <returns>Strand-annotated regulatory elements.</returns>
    public static IEnumerable<StrandedRegulatoryElement> FindRegulatoryElements(DnaSequence sequence, bool bothStrands)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindStrandedRegulatoryElementsCore(sequence, bothStrands);
    }

    private static IEnumerable<StrandedRegulatoryElement> FindStrandedRegulatoryElementsCore(DnaSequence sequence, bool bothStrands)
    {
        foreach (var (name, pattern, description, orientationIndependent) in RegulatoryLibrary)
        {
            var hits = FindDegenerateMotif(sequence, pattern)
                .Select(m => new StrandedRegulatoryElement(name, m.Position, m.MatchedSequence, pattern, description, '+'));

            string rcPattern = DnaSequence.GetReverseComplementString(pattern);
            if (bothStrands && orientationIndependent && rcPattern != pattern)
            {
                var minus = FindDegenerateMotif(sequence, rcPattern)
                    .Select(m => new StrandedRegulatoryElement(name, m.Position,
                        DnaSequence.GetReverseComplementString(m.MatchedSequence), pattern, description, '-'));
                hits = hits.Concat(minus).OrderBy(e => e.Position).ThenBy(e => e.Strand == '+' ? 0 : 1);
            }

            foreach (var element in hits)
                yield return element;
        }
    }

    #endregion
}

/// <summary>
/// A motif match in a sequence.
/// </summary>
public readonly record struct MotifMatch(
    int Position,
    string MatchedSequence,
    string Pattern,
    double Score);

/// <summary>
/// A discovered motif from de novo discovery.
/// </summary>
public readonly record struct DiscoveredMotif(
    string Sequence,
    int Count,
    IReadOnlyList<int> Positions,
    double Enrichment);

/// <summary>
/// A motif shared between multiple sequences.
/// </summary>
public readonly record struct SharedMotif(
    string Sequence,
    IReadOnlyList<int> SequenceIndices,
    double Prevalence);

/// <summary>
/// A regulatory element found in a sequence.
/// </summary>
public readonly record struct RegulatoryElement(
    string Name,
    int Position,
    string Sequence,
    string Pattern,
    string Description);

/// <summary>
/// A strand-annotated regulatory element (<see cref="MotifFinder.FindRegulatoryElements(DnaSequence, bool)"/>).
/// </summary>
/// <param name="Name">Library element name.</param>
/// <param name="Position">0-based forward-strand start of the window.</param>
/// <param name="Sequence">The site read 5'→3' on its own strand (matches <paramref name="Pattern"/>).</param>
/// <param name="Pattern">Library IUPAC pattern.</param>
/// <param name="Description">Element description.</param>
/// <param name="Strand">'+' or '-'.</param>
public readonly record struct StrandedRegulatoryElement(
    string Name,
    int Position,
    string Sequence,
    string Pattern,
    string Description,
    char Strand);

/// <summary>
/// Position Weight Matrix for motif scoring.
/// </summary>
public sealed partial class PositionWeightMatrix
{
    public double[,] Matrix { get; }
    public int Length { get; }
    public string Consensus { get; }

    /// <summary>Creates a PWM from a 4 × <paramref name="length"/> log-odds matrix (rows A, C, G, T).</summary>
    /// <exception cref="ArgumentNullException"><paramref name="matrix"/> is null.</exception>
    /// <exception cref="ArgumentException">The matrix is not 4 × <paramref name="length"/>.</exception>
    public PositionWeightMatrix(double[,] matrix, int length)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        if (length < 0 || matrix.GetLength(0) != 4 || matrix.GetLength(1) != length)
            throw new ArgumentException(
                $"Matrix must be 4 × {length} (rows A, C, G, T); got {matrix.GetLength(0)} × {matrix.GetLength(1)}.",
                nameof(matrix));

        Matrix = matrix;
        Length = length;
        Consensus = GenerateConsensus();
    }

    private string GenerateConsensus()
    {
        var sb = new StringBuilder(Length);

        for (int i = 0; i < Length; i++)
        {
            int maxIdx = 0;
            double maxVal = Matrix[0, i];

            for (int b = 1; b < 4; b++)
            {
                if (Matrix[b, i] > maxVal)
                {
                    maxVal = Matrix[b, i];
                    maxIdx = b;
                }
            }

            sb.Append(MotifFinder.AcgtBases[maxIdx]);
        }

        return sb.ToString();
    }

    /// <summary>
    /// Returns the PWM of the reverse-complement strand: column order reversed and the
    /// A↔T, C↔G rows swapped (Biopython <c>PositionSpecificScoringMatrix.reverse_complement</c>).
    /// Scanning a sequence with it scores the minus strand at the same forward window start.
    /// </summary>
    public PositionWeightMatrix ReverseComplement()
    {
        var rc = new double[4, Length];
        for (int j = 0; j < Length; j++)
        {
            int src = Length - 1 - j;
            for (int b = 0; b < 4; b++)
                rc[3 - b, j] = Matrix[b, src]; // A(0)↔T(3), C(1)↔G(2)
        }
        return new PositionWeightMatrix(rc, Length);
    }

    /// <summary>
    /// Gets the maximum possible score for this PWM.
    /// </summary>
    public double MaxScore
    {
        get
        {
            double max = 0;
            for (int i = 0; i < Length; i++)
            {
                double posMax = double.MinValue;
                for (int b = 0; b < 4; b++)
                    posMax = Math.Max(posMax, Matrix[b, i]);
                max += posMax;
            }
            return max;
        }
    }

    /// <summary>
    /// Gets the minimum possible score for this PWM.
    /// </summary>
    public double MinScore
    {
        get
        {
            double min = 0;
            for (int i = 0; i < Length; i++)
            {
                double posMin = double.MaxValue;
                for (int b = 0; b < 4; b++)
                    posMin = Math.Min(posMin, Matrix[b, i]);
                min += posMin;
            }
            return min;
        }
    }
}
