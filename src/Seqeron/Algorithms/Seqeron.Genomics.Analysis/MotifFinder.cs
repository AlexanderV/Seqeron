using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Finds conserved motifs and patterns in DNA sequences.
/// Supports exact and degenerate motif searching, position weight matrices, and consensus sequences.
/// </summary>
public static class MotifFinder
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

        int[,] counts = BuildCountMatrix(sequences, nameof(sequences), out int count);
        if (count == 0)
            throw new ArgumentException("At least one sequence is required.", nameof(sequences));

        int length = counts.GetLength(1);
        var matrix = new double[PwmAlphabetSize, length];

        // Pseudocount-smoothed probabilities → log2 odds against the background.
        double total = count + PwmAlphabetSize * pseudocount;
        for (int i = 0; i < length; i++)
        {
            for (int b = 0; b < PwmAlphabetSize; b++)
            {
                double freq = (counts[b, i] + pseudocount) / total;
                matrix[b, i] = Math.Log2(freq / bg[b]);
            }
        }

        return new PositionWeightMatrix(matrix, length);
    }

    private static double[] NormalizeBackground(IReadOnlyList<double> background)
    {
        if (background.Count != PwmAlphabetSize)
            throw new ArgumentException(
                "Background must have exactly 4 probabilities (A, C, G, T).", nameof(background));

        double sum = 0;
        for (int b = 0; b < PwmAlphabetSize; b++)
        {
            double v = background[b];
            if (!double.IsFinite(v) || v <= 0)
                throw new ArgumentOutOfRangeException(nameof(background), v,
                    "Background probabilities must be finite and strictly positive.");
            sum += v;
        }

        var result = new double[PwmAlphabetSize];
        for (int b = 0; b < PwmAlphabetSize; b++)
            result[b] = background[b] / sum;
        return result;
    }

    /// <summary>Row index of a PWM for an upper-case base (A=0, C=1, G=2, T=3), or −1.</summary>
    private static int AcgtIndex(char c) => c switch
    {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => -1
    };

    /// <summary>Upper-case bases in PWM row order (inverse of <see cref="AcgtIndex"/>).</summary>
    private static readonly char[] AcgtBases = { 'A', 'C', 'G', 'T' };

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
        var seqList = new List<string>();
        foreach (var s in sequences)
        {
            if (s is null)
                throw new ArgumentException("Sequences cannot contain null elements.", paramName);
            seqList.Add(s);
        }

        sequenceCount = seqList.Count;
        if (sequenceCount == 0)
            return new int[PwmAlphabetSize, 0];

        int length = seqList[0].Length;
        for (int s = 1; s < seqList.Count; s++)
        {
            if (seqList[s].Length != length)
                throw new ArgumentException("All sequences must have the same length.", paramName);
        }

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
    /// non-ACGT symbol are skipped. To scan the reverse strand, scan with
    /// <see cref="PositionWeightMatrix.ReverseComplement"/>; positions are then forward-strand
    /// window starts.
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
            double score = 0;
            bool valid = true;

            for (int j = 0; j < motifLen; j++)
            {
                int baseIndex = AcgtIndex(seq[i + j]);
                if (baseIndex < 0)
                {
                    valid = false;
                    break;
                }

                score += pwm.Matrix[baseIndex, j];
            }

            if (valid && score >= threshold)
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
    /// Generates a consensus sequence from aligned sequences.
    /// </summary>
    /// <param name="sequences">Aligned sequences of equal length.</param>
    /// <returns>Consensus sequence using IUPAC codes.</returns>
    public static string GenerateConsensus(IEnumerable<string> sequences)
    {
        ArgumentNullException.ThrowIfNull(sequences);

        var seqList = sequences.Select(s => s.ToUpperInvariant()).ToList();
        if (seqList.Count == 0) return "";

        int length = seqList[0].Length;
        var consensus = new StringBuilder(length);

        for (int i = 0; i < length; i++)
        {
            var counts = new Dictionary<char, int> { ['A'] = 0, ['C'] = 0, ['G'] = 0, ['T'] = 0 };

            foreach (var seq in seqList)
            {
                if (i < seq.Length && counts.TryGetValue(seq[i], out int value))
                    counts[seq[i]] = ++value;
            }

            consensus.Append(GetIupacCode(counts, seqList.Count));
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
    /// exceed to be included in the position's IUPAC degeneracy code. The "combine the bases
    /// that pass a frequency threshold into the IUPAC symbol for that set" rule is the standard
    /// threshold-consensus mechanism (Bioconductor DECIPHER <c>ConsensusSequence</c>: "removes
    /// the least frequent characters … so long as they represent less than <c>threshold</c> …");
    /// the set→symbol mapping itself follows NC-IUB 1984 (Cornish-Bowden, NAR 13(9):3021). The
    /// 0.25 cut and strict '&gt;' boundary are this implementation's documented design constant
    /// (a base must be present in more than a quarter of the sequences). See
    /// docs/Evidence/MOTIF-GENERATE-001-Evidence.md.
    /// </summary>
    private const double IupacInclusionThreshold = 0.25;

    private static char GetIupacCode(Dictionary<char, int> counts, int total)
    {
        double threshold = total * IupacInclusionThreshold; // base count must be strictly > threshold

        var present = counts.Where(kv => kv.Value > threshold)
                           .Select(kv => kv.Key)
                           .OrderBy(c => c)
                           .ToList();

        if (present.Count == 0)
        {
            return counts.MaxBy(kv => kv.Value).Key;
        }

        // Base set → NC-IUB symbol via the canonical inverse map; present is a non-empty
        // subset of {A,C,G,T}, and the full set {A,C,G,T} maps to N.
        return Core.IupacDnaSequence.GetIupacCode(present);
    }

    #endregion

    #region Motif Discovery

    /// <summary>
    /// Size of the DNA nucleotide alphabet {A, C, G, T}; the base of the
    /// 4^k count of distinct k-mers used in the expected-occurrence formula.
    /// </summary>
    private const int DnaAlphabetSize = 4;

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
    /// value &gt; 1 means the k-mer occurs more often than chance predicts.
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
        return DiscoverMotifsCore(sequence, k, minCount);
    }

    private static IEnumerable<DiscoveredMotif> DiscoverMotifsCore(DnaSequence sequence, int k, int minCount)
    {

        string seq = sequence.Sequence;
        var kmerPositions = new Dictionary<string, List<int>>();

        // Count k-mers at every length-k window (0-based start positions).
        for (int i = 0; i <= seq.Length - k; i++)
        {
            string kmer = seq.Substring(i, k);

            if (!kmerPositions.ContainsKey(kmer))
                kmerPositions[kmer] = new List<int>();

            kmerPositions[kmer].Add(i);
        }

        // Expected occurrences of a specific k-mer under the i.i.d. uniform background:
        // E = (N - k + 1) / 4^k  (Compeau & Pevzner, Bioinformatics Algorithms).
        // windowCount = N - k + 1 is the number of length-k windows; it is >= 1 whenever
        // at least one k-mer was counted, so E is strictly positive here.
        double windowCount = seq.Length - k + 1.0;
        double expectedCount = windowCount / Math.Pow(DnaAlphabetSize, k);

        // Return overrepresented k-mers with their observed/expected (O/E) ratio.
        foreach (var (kmer, positions) in kmerPositions)
        {
            if (positions.Count >= minCount)
            {
                double enrichment = positions.Count / expectedCount;

                yield return new DiscoveredMotif(
                    Sequence: kmer,
                    Count: positions.Count,
                    Positions: positions.AsReadOnly(),
                    Enrichment: enrichment);
            }
        }
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
    /// <param name="sequences">Collection of DNA sequences (each scanned for its distinct words).</param>
    /// <param name="k">Word (oligonucleotide) length; must be ≥ 1. Default: 6.</param>
    /// <param name="minSequences">Quorum: minimum number of distinct sequences a word must occur in. Must be ≥ 1. Default: 2.</param>
    /// <returns>Shared words with their distinct sequence indices and prevalence (matching sequences / total sequences).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">When <paramref name="k"/> &lt; 1 or <paramref name="minSequences"/> &lt; 1.</exception>
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
        var kmerOccurrences = new Dictionary<string, List<int>>();

        // Find k-mers in each sequence
        for (int seqIdx = 0; seqIdx < seqList.Count; seqIdx++)
        {
            var seq = seqList[seqIdx].Sequence;
            var seenInSeq = new HashSet<string>();

            for (int i = 0; i <= seq.Length - k; i++)
            {
                string kmer = seq.Substring(i, k);

                if (seenInSeq.Add(kmer))
                {

                    if (!kmerOccurrences.ContainsKey(kmer))
                        kmerOccurrences[kmer] = new List<int>();

                    kmerOccurrences[kmer].Add(seqIdx);
                }
            }
        }

        // Return shared motifs
        foreach (var (kmer, seqIndices) in kmerOccurrences)
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
    /// Known regulatory motif patterns.
    /// </summary>
    public static class KnownMotifs
    {
        // Eukaryotic core-promoter consensus per Bucher (1990) weight-matrix analysis of 502 promoters.
        /// <summary>TATA box consensus: TATAAA (eukaryotic RNA Pol II core promoter, Bucher 1990).</summary>
        public const string TataBox = "TATAAA";

        /// <summary>CCAAT box consensus pentanucleotide: CCAAT (Bucher 1990).</summary>
        public const string CaatBox = "CCAAT";

        /// <summary>GC box (Sp1) consensus: GGGCGG (Lundin, Nehlin &amp; Ronne 1994).</summary>
        public const string GcBox = "GGGCGG";

        // Prokaryotic sigma-70 promoter hexamers per Harley &amp; Reynolds (1987) compilation.
        /// <summary>-10 (Pribnow) box consensus hexamer: TATAAT (Pribnow 1975; Harley &amp; Reynolds 1987).</summary>
        public const string MinusTenBox = "TATAAT";

        /// <summary>-35 box consensus hexamer: TTGACA (Harley &amp; Reynolds 1987).</summary>
        public const string MinusThirtyFiveBox = "TTGACA";

        /// <summary>Kozak optimal-context sequence: GCCGCCACCATGG (Kozak 1987, most-preferred bases -9..+4).</summary>
        public const string Kozak = "GCCGCCACCATGG";

        /// <summary>Shine-Dalgarno (bacterial RBS) consensus: AGGAGG (complementary to 3' end of 16S rRNA).</summary>
        public const string ShineDalgarno = "AGGAGG";

        /// <summary>Poly(A) signal hexamer: AATAAA (Proudfoot &amp; Brownlee 1976).</summary>
        public const string PolyASignal = "AATAAA";

        /// <summary>E-box consensus (IUPAC): CANNTG (Massari &amp; Murre 2000).</summary>
        public const string EBox = "CANNTG";

        /// <summary>AP-1 (TRE) recognition motif: TGACTCA (Lee, Mitchell &amp; Tjian 1987).</summary>
        public const string Ap1 = "TGACTCA";

        /// <summary>NF-κB κB site: GGGACTTTCC (consensus GGGRNWYYCC; Sen &amp; Baltimore 1986).</summary>
        public const string NfKb = "GGGACTTTCC";

        /// <summary>CREB CRE palindrome: TGACGTCA (Montminy et al. 1986).</summary>
        public const string Creb = "TGACGTCA";
    }

    /// <summary>
    /// Scans for known regulatory motifs.
    /// </summary>
    /// <param name="sequence">DNA sequence to scan.</param>
    /// <returns>Found regulatory elements.</returns>
    public static IEnumerable<RegulatoryElement> FindRegulatoryElements(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindRegulatoryElementsCore(sequence);
    }

    private static IEnumerable<RegulatoryElement> FindRegulatoryElementsCore(DnaSequence sequence)
    {
        var patterns = new (string Name, string Pattern, string Description)[]
        {
            ("TATA Box", KnownMotifs.TataBox, "Eukaryotic core promoter element"),
            ("CAAT Box", KnownMotifs.CaatBox, "Promoter element"),
            ("GC Box", KnownMotifs.GcBox, "Sp1 binding site"),
            ("-10 Box", KnownMotifs.MinusTenBox, "Prokaryotic Pribnow box"),
            ("-35 Box", KnownMotifs.MinusThirtyFiveBox, "Prokaryotic -35 promoter element"),
            ("Kozak", KnownMotifs.Kozak, "Translation initiation"),
            ("Shine-Dalgarno", KnownMotifs.ShineDalgarno, "Bacterial ribosome binding"),
            ("Poly(A) Signal", KnownMotifs.PolyASignal, "Polyadenylation signal"),
            ("E-box", KnownMotifs.EBox, "bHLH transcription factor binding"),
            ("AP-1", KnownMotifs.Ap1, "AP-1 transcription factor binding"),
            ("NF-κB", KnownMotifs.NfKb, "NF-κB binding site"),
            ("CREB", KnownMotifs.Creb, "CREB transcription factor binding")
        };

        foreach (var (name, pattern, description) in patterns)
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
/// Position Weight Matrix for motif scoring.
/// </summary>
public sealed class PositionWeightMatrix
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
        char[] bases = { 'A', 'C', 'G', 'T' };

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

            sb.Append(bases[maxIdx]);
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
