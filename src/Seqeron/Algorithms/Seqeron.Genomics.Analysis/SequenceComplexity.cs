namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Finite-sample bias correction applied by
/// <see cref="SequenceComplexity.CalculateKmerEntropy(string, int, KmerEntropyCorrection, bool)"/>.
/// </summary>
public enum KmerEntropyCorrection
{
    /// <summary>Plug-in (maximum-likelihood) estimate −Σ p ln p, no correction (Shannon 1948).</summary>
    None = 0,

    /// <summary>Miller (1955): plug-in + (D − 1)/(2N) nats, D = number of observed k-mers, N = number of k-mers.</summary>
    MillerMadow = 1,

    /// <summary>Grassberger (2003, arXiv:physics/0307138, eq. 35): ln N − (1/N) Σ n_i G(n_i).</summary>
    Grassberger = 2,
}

/// <summary>
/// Calculates various sequence complexity metrics for detecting low-complexity regions,
/// repetitive sequences, and information content.
/// </summary>
public static partial class SequenceComplexity
{
    #region Linguistic Complexity

    /// <summary>
    /// Calculates linguistic complexity (LC), the summation form
    /// LC = Σ_{i=1..m} V_i / Σ_{i=1..m} V_max,i with V_max,i = min(a^i, N − i + 1),
    /// where V_i is the number of distinct subwords of length i and a is the alphabet size — 4 for DNA/RNA,
    /// extended by any other symbol present (e.g. N), so LC ≤ 1 always (Orlov &amp; Potapov 2004, NAR 32:W628,
    /// word length limited by m ≤ N). With m ≥ N this is exactly the Troyanskaya et al. (2002,
    /// Bioinformatics 18:679) definition LC = A(s)/M(s) over all lengths 1..N (Rosalind LING).
    /// This is not Trifonov's (1990) product form C = Π U_i.
    /// LC = 1.0 for maximum complexity, lower values indicate repeats/low complexity.
    /// </summary>
    /// <remarks>
    /// Small m uses direct hash enumeration; larger m counts V_i from the sequence's suffix tree
    /// (Troyanskaya et al. 2002) via the shared <c>ISuffixTree.CountDistinctSubstringsByLength</c>
    /// (V_i = number of suffix-tree edges spanning depth i), so the full-length LC is linear-time. Both paths return identical values.
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="maxWordLength">Maximum word length m to consider (values ≥ N give Troyanskaya's all-length LC).</param>
    /// <returns>Linguistic complexity (0 to 1).</returns>
    public static double CalculateLinguisticComplexity(DnaSequence sequence, int maxWordLength = 10)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxWordLength, 1);

        return CalculateLinguisticComplexityDna(sequence, maxWordLength, LcAlphabetSize(sequence.Sequence));
    }

    /// <summary>
    /// Calculates linguistic complexity over a caller-supplied alphabet of size <paramref name="alphabetSize"/>:
    /// LC = Σ_{i=1..m} V_i / Σ_{i=1..m} min(a^i, N − i + 1) with a fixed (Troyanskaya et al. 2002, Bioinformatics 18:679,
    /// "a text of size n over an alphabet of cardinality a"; Rosalind LING: m(a, n) with a = 4 for DNA). Unlike the
    /// two-argument overload the alphabet is not inferred from the sequence, so e.g. a DNA window that happens to lack
    /// T is still scored against a = 4, and a protein can be scored with a = 20.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="maxWordLength">Maximum word length m (≥ 1; values ≥ N give the all-length LC).</param>
    /// <param name="alphabetSize">Alphabet size a (≥ 1 and ≥ the number of distinct symbols in the sequence).</param>
    /// <returns>Linguistic complexity in [0, 1]; 0 for an empty sequence.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxWordLength"/> or <paramref name="alphabetSize"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">The sequence contains more distinct symbols than <paramref name="alphabetSize"/>.</exception>
    public static double CalculateLinguisticComplexity(DnaSequence sequence, int maxWordLength, int alphabetSize)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxWordLength, 1);
        ValidateLcAlphabetSize(sequence.Sequence, alphabetSize);

        return CalculateLinguisticComplexityDna(sequence, maxWordLength, alphabetSize);
    }

    /// <summary>
    /// Fixed-alphabet linguistic complexity of a raw string (see
    /// <see cref="CalculateLinguisticComplexity(DnaSequence, int, int)"/>); the input is upper-cased, and the number of
    /// distinct (upper-cased) symbols must not exceed <paramref name="alphabetSize"/>. Null/empty input returns 0.
    /// Rosalind LING sample: <c>ATTTGGATT</c>, a = 4, m ≥ 9 → 35/40 = 0.875.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxWordLength"/> or <paramref name="alphabetSize"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">The sequence contains more distinct symbols than <paramref name="alphabetSize"/>.</exception>
    public static double CalculateLinguisticComplexity(string sequence, int maxWordLength, int alphabetSize)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(maxWordLength, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(alphabetSize, 1);
        if (string.IsNullOrEmpty(sequence)) return 0;

        string seq = sequence.ToUpperInvariant();
        ValidateLcAlphabetSize(seq, alphabetSize);
        return CalculateLinguisticComplexityCore(seq, maxWordLength, alphabetSize);
    }

    private static double CalculateLinguisticComplexityDna(DnaSequence sequence, int maxWordLength, int alphabetSize)
    {
        string seq = sequence.Sequence;
        if (seq.Length == 0) return 0;
        int m = Math.Min(maxWordLength, seq.Length);
        return m > LcHashEnumerationMaxWordLength
            ? LinguisticComplexityFromCounts(
                sequence.SuffixTree.CountDistinctSubstringsByLength(m), seq.Length, alphabetSize)
            : CalculateLinguisticComplexityCore(seq, maxWordLength, alphabetSize);
    }

    private static void ValidateLcAlphabetSize(string seq, int alphabetSize)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(alphabetSize, 1);
        int distinct = new HashSet<char>(seq).Count;
        if (distinct > alphabetSize)
            throw new ArgumentException(
                $"The sequence contains {distinct} distinct symbols, more than alphabetSize = {alphabetSize}.",
                nameof(alphabetSize));
    }

    /// <summary>
    /// Calculates linguistic complexity from a raw sequence string (same definition as the
    /// <see cref="CalculateLinguisticComplexity(DnaSequence, int)"/> overload; input is upper-cased;
    /// null/empty input or <paramref name="maxWordLength"/> &lt; 1 returns 0).
    /// </summary>
    public static double CalculateLinguisticComplexity(string sequence, int maxWordLength = 10)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        string seq = sequence.ToUpperInvariant();
        return CalculateLinguisticComplexityCore(seq, maxWordLength, LcAlphabetSize(seq));
    }

    /// <summary>Largest word length for which direct hash enumeration is used instead of the suffix tree.</summary>
    private const int LcHashEnumerationMaxWordLength = 12;

    private static double CalculateLinguisticComplexityCore(string seq, int maxWordLength, int alphabetSize)
    {
        if (seq.Length == 0) return 0;

        int m = Math.Min(maxWordLength, seq.Length);
        if (m < 1) return 0;

        if (m > LcHashEnumerationMaxWordLength)
            return LinguisticComplexityFromCounts(
                global::SuffixTree.SuffixTree.Build(seq).CountDistinctSubstringsByLength(m), seq.Length, alphabetSize);

        // V_i = number of distinct overlapping i-words = key count of the canonical k-mer tally (KMER-COUNT-001).
        var counts = new long[m + 1];
        for (int wordLen = 1; wordLen <= m; wordLen++)
            counts[wordLen] = KmerAnalyzer.CountKmers(seq, wordLen).Count;

        return LinguisticComplexityFromCounts(counts, seq.Length, alphabetSize);
    }

    /// <summary>
    /// Alphabet size <c>a</c> of the LC denominator (Troyanskaya et al. 2002; Rosalind LING: M = Σ min(a^i, N − i + 1)
    /// "for an alphabet of size a"): the nucleotide alphabet {A, C, G, T} (U in place of T for RNA, i.e. when U occurs
    /// and T does not) extended by every other symbol that occurs in the (upper-cased) sequence. Pure DNA/RNA gives
    /// a = 4; e.g. <c>ACGTN</c> gives a = 5. Because every observed symbol is in the alphabet, V_i ≤ min(a^i, N − i + 1)
    /// and LC ≤ 1 for any input.
    /// </summary>
    private static int LcAlphabetSize(string seq)
    {
        var symbols = new HashSet<char>(seq) { 'A', 'C', 'G' };
        if (!(symbols.Contains('U') && !symbols.Contains('T')))
            symbols.Add('T');
        return symbols.Count;
    }

    /// <summary>
    /// LC = Σ V_i / Σ min(a^i, N − i + 1) for i = 1..counts.Length−1 (a = alphabet size, see <see cref="LcAlphabetSize"/>).
    /// </summary>
    private static double LinguisticComplexityFromCounts(long[] counts, int n, int alphabetSize)
    {
        long wordsOverAlphabet = 1; // a^i, saturated once it exceeds N (then min(a^i, N − i + 1) = N − i + 1)
        long observedTotal = 0;
        long possibleTotal = 0;

        for (int wordLen = 1; wordLen < counts.Length; wordLen++)
        {
            observedTotal += counts[wordLen];

            // V_max,i = min(a^i, N − i + 1); a ≤ int.MaxValue and a^(i−1) ≤ N ≤ int.MaxValue before the multiply,
            // so the product stays below 2^62 (no overflow) and saturation keeps it there.
            if (wordsOverAlphabet <= n)
                wordsOverAlphabet *= alphabetSize;
            long positions = n - wordLen + 1;
            possibleTotal += Math.Min(wordsOverAlphabet, positions);
        }

        return possibleTotal > 0 ? (double)observedTotal / possibleTotal : 0;
    }

    #endregion

    #region Shannon Entropy

    /// <summary>
    /// Calculates Shannon entropy for the sequence (bits per base).
    /// Maximum entropy for DNA is 2 bits (log2(4)).
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <returns>Shannon entropy (0 to 2 for DNA).</returns>
    public static double CalculateShannonEntropy(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateShannonEntropyCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates Shannon entropy H = −Σ p_i·log₂(p_i) (bits per base, Shannon 1948) from a raw
    /// nucleotide string over the 4-letter nucleotide alphabet {A, C, G, T/U}.
    /// Case-insensitive; RNA uracil (U) is counted as the fourth nucleotide, the same class as
    /// DNA thymine (IUPAC-IUB 1970 nucleotide notation; consistent with the canonical GC counter
    /// <c>SequenceExtensions.CountGcAndValidNucleotides</c>), so the maximum stays log₂4 = 2 bits.
    /// Other symbols (N, IUPAC ambiguity codes, gaps) are excluded from numerator and denominator.
    /// Null, empty, or input with no A/C/G/T/U returns 0 (empty-sum convention).
    /// </summary>
    public static double CalculateShannonEntropy(string sequence)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateShannonEntropyCore(sequence.ToUpperInvariant());
    }

    private static double CalculateShannonEntropyCore(string seq)
    {
        if (seq.Length == 0) return 0;

        // Counts over {A, C, G, T/U}: RNA U is the same nucleotide class as DNA T (IUPAC-IUB 1970).
        int a = 0, c = 0, g = 0, t = 0;
        foreach (char ch in seq)
        {
            switch (ch)
            {
                case 'A': a++; break;
                case 'C': c++; break;
                case 'G': g++; break;
                case 'T':
                case 'U': t++; break;
            }
        }

        return ShannonEntropyBits([a, c, g, t]);
    }

    /// <summary>
    /// Shannon entropy H = −Σ p_i·log₂(p_i), p_i = n_i / Σ n, of a frequency table (Shannon 1948), in bits.
    /// Zero counts contribute nothing (0·log 0 := 0); an empty or all-zero table has entropy 0.
    /// Delegates to the canonical <see cref="StatisticsHelper.ShannonIndex"/> (natural log) and converts
    /// to bits by dividing by ln 2 — the same computation as <c>scipy.stats.entropy(counts, base=2)</c>
    /// (shared with <c>SequenceStatistics.CalculateShannonEntropy</c>). <c>ShannonIndex</c> rejects a zero
    /// total, so the empty-sum convention (entropy 0) is handled here.
    /// </summary>
    private static double ShannonEntropyBits(IReadOnlyList<int> counts)
    {
        long total = 0;
        foreach (int count in counts)
            total += count;

        return total == 0 ? 0 : StatisticsHelper.ShannonIndex(counts) / Ln2;
    }

    // Base conversion ln → log₂ (bits): H₂ = H_e / ln 2 (scipy.stats.entropy divides by log(base)).
    private static readonly double Ln2 = Math.Log(2.0);

    /// <summary>
    /// Calculates the Shannon entropy (in bits) of the overlapping k-mer frequency
    /// distribution of the sequence.
    /// </summary>
    /// <remarks>
    /// The sequence is decomposed into its L-k+1 overlapping k-mers (sliding window,
    /// one base step). With n_i the count of distinct k-mer i and N = L-k+1 the total
    /// number of k-mers, p_i = n_i / N and the entropy is H = -Σ p_i · log₂(p_i)
    /// (Shannon 1948). Entropy is reported in bits (log base 2): it is 0 when only one
    /// distinct k-mer occurs (deterministic distribution) and reaches log₂(N) when every k-mer
    /// is distinct (uniform distribution). This is the order-k "block entropy" H_k of the overlapping k-word distribution used in DNA
    /// entropy analysis (Herzel, Ebeling &amp; Schmitt 1994, Phys. Rev. E 50:5061; Schmitt &amp; Herzel 1997,
    /// J. Theor. Biol. 188:369), the same raw quantity BBDuk's EntropyTracker computes before its
    /// 1/ln(N) normalisation. It is the plug-in (maximum-likelihood) estimate: no finite-sample bias
    /// correction is applied, so it underestimates the true block entropy when N ≪ 4^k (Schmitt &amp;
    /// Herzel 1997). Every length-k substring is a symbol (no IUPAC/N filtering).
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="k">K-mer size (default: 2 for dinucleotides). Must be ≥ 1.</param>
    /// <returns>Shannon entropy (bits) of the k-mer frequency distribution; 0 when L &lt; k.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> &lt; 1.</exception>
    public static double CalculateKmerEntropy(DnaSequence sequence, int k = 2)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);

        return CalculateKmerEntropyCore(sequence.Sequence, k);
    }

    /// <summary>
    /// Calculates the Shannon entropy (in bits) of the overlapping k-mer frequency
    /// distribution from a raw sequence string. The string is upper-cased to match the
    /// normalization applied by <see cref="DnaSequence"/>.
    /// </summary>
    /// <param name="sequence">Raw sequence string; null or empty yields 0.</param>
    /// <param name="k">K-mer size (default: 2 for dinucleotides). Must be ≥ 1.</param>
    /// <returns>Shannon entropy (bits) of the k-mer frequency distribution; 0 when L &lt; k.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> &lt; 1.</exception>
    public static double CalculateKmerEntropy(string sequence, int k = 2)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateKmerEntropyCore(sequence.ToUpperInvariant(), k);
    }

    private static double CalculateKmerEntropyCore(string seq, int k)
        => CalculateKmerEntropyCore(seq, k, KmerEntropyCorrection.None, normalize: false);

    /// <summary>
    /// K-mer (block) entropy in bits with an optional finite-sample bias correction and optional normalisation.
    /// </summary>
    /// <remarks>
    /// With n_i the counts of the D distinct overlapping k-mers and N = L − k + 1 = Σ n_i (all in nats, converted to bits
    /// by ÷ ln 2):
    /// <list type="bullet">
    /// <item><see cref="KmerEntropyCorrection.None"/>: plug-in H_ML = −Σ (n_i/N) ln(n_i/N) — identical to
    /// <see cref="CalculateKmerEntropy(DnaSequence, int)"/>.</item>
    /// <item><see cref="KmerEntropyCorrection.MillerMadow"/>: H_MM = H_ML + (D − 1)/(2N) (Miller 1955; D = number of
    /// k-mers actually observed — R <c>entropy::entropy.MillerMadow</c>, <c>m = sum(y &gt; 0)</c>).</item>
    /// <item><see cref="KmerEntropyCorrection.Grassberger"/>: H_G = ln N − (1/N) Σ n_i G(n_i) with
    /// G(n) = ψ(n) + ½(−1)ⁿ[ψ((n + 1)/2) − ψ(n/2)] (Grassberger 2003, arXiv:physics/0307138, eq. 35; the estimator of
    /// Python <c>ndd.estimators.Grassberger</c>). For N = 1 it is γ + ln 2 nats (1.8327 bits), not 0.</item>
    /// </list>
    /// The bias of the plug-in block entropy H_n for N ≪ 4^k is the problem these estimators address (Herzel, Schmitt
    /// &amp; Ebeling 1994; Schmitt &amp; Herzel 1997). With <paramref name="normalize"/> the (possibly corrected) value is
    /// divided by log₂ N, the maximum of the plug-in estimate (all N k-mers distinct) — the 0–1 scale of BBTools
    /// <c>EntropyTracker.calcEntropy</c> (multiplier 1/ln(windowKmers)). A corrected value can exceed 1 after
    /// normalisation. When N ≤ 1 the normaliser log₂ N is 0 and the normalised result is defined as 0.
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="k">K-mer size (≥ 1).</param>
    /// <param name="correction">Bias correction (<see cref="KmerEntropyCorrection.None"/> = plug-in).</param>
    /// <param name="normalize">Divide by log₂ N (N = L − k + 1).</param>
    /// <returns>Entropy in bits (or normalised); 0 when L &lt; k.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> &lt; 1 or <paramref name="correction"/> undefined.</exception>
    public static double CalculateKmerEntropy(
        DnaSequence sequence, int k, KmerEntropyCorrection correction, bool normalize = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateKmerEntropyOptions(k, correction);
        return CalculateKmerEntropyCore(sequence.Sequence, k, correction, normalize);
    }

    /// <summary>
    /// String form of <see cref="CalculateKmerEntropy(DnaSequence, int, KmerEntropyCorrection, bool)"/>; the input is
    /// upper-cased, every length-k substring is a symbol; null/empty input returns 0.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> &lt; 1 or <paramref name="correction"/> undefined.</exception>
    public static double CalculateKmerEntropy(
        string sequence, int k, KmerEntropyCorrection correction, bool normalize = false)
    {
        ValidateKmerEntropyOptions(k, correction);
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateKmerEntropyCore(sequence.ToUpperInvariant(), k, correction, normalize);
    }

    /// <summary>
    /// Parses an MCP / text correction name (case-insensitive): <c>none</c>, <c>millermadow</c> (or <c>miller-madow</c>,
    /// <c>mm</c>), <c>grassberger</c>.
    /// </summary>
    /// <exception cref="ArgumentException">Unknown name.</exception>
    public static KmerEntropyCorrection ParseKmerEntropyCorrection(string? name)
    {
        string key = (name ?? "none").Trim().Replace("-", "").Replace("_", "").ToLowerInvariant();
        return key switch
        {
            "" or "none" or "plugin" or "ml" => KmerEntropyCorrection.None,
            "millermadow" or "mm" => KmerEntropyCorrection.MillerMadow,
            "grassberger" => KmerEntropyCorrection.Grassberger,
            _ => throw new ArgumentException(
                $"Unknown k-mer entropy correction '{name}'; expected none, millerMadow or grassberger.", nameof(name)),
        };
    }

    private static void ValidateKmerEntropyOptions(int k, KmerEntropyCorrection correction)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        if (!Enum.IsDefined(correction))
            throw new ArgumentOutOfRangeException(nameof(correction), correction, "Undefined k-mer entropy correction.");
    }

    private static double CalculateKmerEntropyCore(string seq, int k, KmerEntropyCorrection correction, bool normalize)
    {
        if (seq.Length < k) return 0;

        // Overlapping k-mer tally (N = L − k + 1 windows) via the canonical counter (KMER-COUNT-001);
        // p_i = n_i / N because Σ n_i = N.
        int[] counts = KmerAnalyzer.CountKmers(seq, k).Values.ToArray();
        int n = seq.Length - k + 1;

        double bits = correction switch
        {
            KmerEntropyCorrection.MillerMadow => ShannonEntropyBits(counts) + (counts.Length - 1) / (2.0 * n) / Ln2,
            KmerEntropyCorrection.Grassberger => GrassbergerEntropyNats(counts, n) / Ln2,
            _ => ShannonEntropyBits(counts),
        };

        if (!normalize) return bits;
        return n > 1 ? bits / Math.Log2(n) : 0;
    }

    /// <summary>Grassberger (2003) eq. 35: H = ln N − (1/N) Σ n_i G(n_i), nats.</summary>
    private static double GrassbergerEntropyNats(IReadOnlyList<int> counts, int total)
    {
        double sum = 0;
        foreach (int c in counts)
            sum += c * GrassbergerG(c);
        return Math.Log(total) - sum / total;
    }

    /// <summary>G(n) = ψ(n) + ½(−1)ⁿ[ψ((n + 1)/2) − ψ(n/2)] (Grassberger 2003, eq. 35), n ≥ 1.</summary>
    internal static double GrassbergerG(int n)
    {
        double half = 0.5 * (StatisticsHelper.Digamma((n + 1) / 2.0) - StatisticsHelper.Digamma(n / 2.0));
        return StatisticsHelper.Digamma(n) + ((n & 1) == 0 ? half : -half);
    }

    #endregion

    #region Sliding Window Complexity

    // Default per-window linguistic-complexity word-length cap m. Gabrielian & Bolshoy (1999, Comput. Chem.
    // 23:263-274, doi:10.1016/S0097-8485(99)00007-8) bound the word lengths of the windowed LC "not in the range
    // of 2 to N-1 but only up to W" for efficiency; the value W is a free parameter there (universalmotif
    // sequence_complexity uses trifonov.max.word.size = 7). 6 is this library's historical default, kept for
    // backward compatibility and exposed as the lcMaxWordLength parameter.
    private const int WindowLcMaxWordLength = 6;

    /// <summary>
    /// Calculates complexity across the sequence using a sliding window (a complexity
    /// profile, in the sense of Troyanskaya et al. (2002)). For each window fully contained
    /// in the sequence the per-window Shannon entropy (bits, Shannon 1948) and linguistic
    /// complexity (summation form, word lengths 1..min(<paramref name="lcMaxWordLength"/>, w),
    /// windowed LC with a bounded word length as in Gabrielian &amp; Bolshoy 1999) are reported with the window's coordinates.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the sliding window (default: 64).</param>
    /// <param name="stepSize">Step size for window movement (default: 10).</param>
    /// <param name="lcMaxWordLength">Per-window LC word-length cap m (default 6, the library's historical value; ≥ 1).
    /// Values ≥ <paramref name="windowSize"/> give the all-length LC of each window (Troyanskaya et al. 2002).</param>
    /// <returns>Complexity values with positions.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/>, <paramref name="stepSize"/> or
    /// <paramref name="lcMaxWordLength"/> &lt; 1.</exception>
    public static IEnumerable<ComplexityPoint> CalculateWindowedComplexity(
        DnaSequence sequence,
        int windowSize = 64,
        int stepSize = 10,
        int lcMaxWordLength = WindowLcMaxWordLength)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateWindowedParameters(windowSize, stepSize, lcMaxWordLength);

        return CalculateWindowedComplexityCore(sequence.Sequence, windowSize, stepSize, lcMaxWordLength);
    }

    /// <summary>
    /// Sliding-window complexity profile of a raw nucleotide string (same metrics as the
    /// <see cref="CalculateWindowedComplexity(DnaSequence, int, int, int)"/> overload; input is upper-cased).
    /// Windows that contain an undefined symbol (anything other than A/C/G/T/U — N, IUPAC ambiguity codes, gaps)
    /// are not evaluated and produce no point, following BBTools BBDuk, which scores only windows with
    /// <c>EntropyTracker.ns() &lt; 1</c> (<c>maskLowEntropy</c>/<c>markLowEntropy</c>). For A/C/G/T input the result
    /// is identical to the <see cref="DnaSequence"/> overload.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/>, <paramref name="stepSize"/> or
    /// <paramref name="lcMaxWordLength"/> &lt; 1.</exception>
    public static IEnumerable<ComplexityPoint> CalculateWindowedComplexity(
        string sequence,
        int windowSize = 64,
        int stepSize = 10,
        int lcMaxWordLength = WindowLcMaxWordLength)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateWindowedParameters(windowSize, stepSize, lcMaxWordLength);

        return CalculateWindowedComplexityCore(sequence.ToUpperInvariant(), windowSize, stepSize, lcMaxWordLength);
    }

    private static void ValidateWindowedParameters(int windowSize, int stepSize, int lcMaxWordLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(lcMaxWordLength, 1);
    }

    private static IEnumerable<ComplexityPoint> CalculateWindowedComplexityCore(
        string seq,
        int windowSize,
        int stepSize,
        int lcMaxWordLength)
    {
        int[] undefinedPrefix = UndefinedBasePrefixCounts(seq);
        int m = Math.Min(lcMaxWordLength, windowSize);

        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            if (undefinedPrefix[i + windowSize] != undefinedPrefix[i]) continue; // BBDuk: only windows with ns() < 1

            string window = seq.Substring(i, windowSize);
            double entropy = CalculateShannonEntropyCore(window);
            double lc = WindowLinguisticComplexity(window, m);

            yield return new ComplexityPoint(
                Position: i + windowSize / 2,
                ShannonEntropy: entropy,
                LinguisticComplexity: lc,
                WindowStart: i,
                WindowEnd: i + windowSize - 1);
        }
    }

    /// <summary>
    /// Smallest per-window word-length cap m for which the window's V_1..V_m are read from the window's suffix tree
    /// instead of m hash enumerations. Measured on 1 Mb random DNA (w = 64, s = 10): m = 1/2/3 hash 320/665/1 067 ms vs
    /// tree 1 104/977/1 204 ms; m = 4/6 hash 1 485/2 155 ms vs tree 1 192/1 090 ms (tree cost is independent of m).
    /// </summary>
    private const int WindowLcSuffixTreeMinWordLength = 4;

    /// <summary>
    /// LC of one window. For m ≥ <see cref="WindowLcSuffixTreeMinWordLength"/> the subword counts V_i come from the
    /// window's suffix tree (Troyanskaya et al. 2002: V_i = number of edges spanning string depth i, via the shared
    /// <c>ISuffixTree.CountDistinctSubstringsByLength</c>), O(w) per window whatever m is; otherwise from
    /// <see cref="KmerAnalyzer.CountKmers"/>. Both produce the same integer counts, so the value is bit-identical.
    /// </summary>
    private static double WindowLinguisticComplexity(string window, int m)
    {
        int alphabetSize = LcAlphabetSize(window);
        return m >= WindowLcSuffixTreeMinWordLength
            ? LinguisticComplexityFromCounts(
                global::SuffixTree.SuffixTree.Build(window).CountDistinctSubstringsByLength(m), window.Length, alphabetSize)
            : CalculateLinguisticComplexityCore(window, m, alphabetSize);
    }

    /// <summary>
    /// BBTools nucleotide code (<c>AminoAcid.baseToNumber</c>): A/C/G/T/U in either case → 0..3 (U = T),
    /// anything else → −1 ("not fully defined"; counted by <c>EntropyTracker.ns()</c>).
    /// </summary>
    private static int BbtoolsBaseCode(char c) => c is 'U' or 'u' ? 3 : AcgtCode(c);

    /// <summary>prefix[i] = number of undefined symbols (see <see cref="BbtoolsBaseCode"/>) in seq[0..i).</summary>
    private static int[] UndefinedBasePrefixCounts(string seq)
    {
        var prefix = new int[seq.Length + 1];
        for (int i = 0; i < seq.Length; i++)
            prefix[i + 1] = prefix[i] + (BbtoolsBaseCode(seq[i]) < 0 ? 1 : 0);
        return prefix;
    }

    #endregion

    #region Low Complexity Regions

    /// <summary>
    /// Finds low-complexity regions by a per-base Shannon-entropy window scan: every window of length
    /// <paramref name="windowSize"/> (step 1) whose per-base Shannon entropy (bits, Shannon 1948; the canonical
    /// <see cref="CalculateShannonEntropy(DnaSequence)"/> kernel, i.e. 1-mer entropy, not normalised) is strictly below
    /// <paramref name="entropyThreshold"/> is flagged, and a region is a maximal run of positions covered by flagged
    /// windows (the union of the low-entropy windows; overlapping or abutting flagged windows merge).
    /// </summary>
    /// <remarks>
    /// The window-union reporting rule is the one of BBTools BBDuk <c>maskLowEntropy</c> (each failing window sets bits
    /// [left, right] of a bit mask). The window statistic is <b>not</b> BBDuk's default: BBDuk scores k-mer entropy
    /// (entropyk = 5, entropywindow = 50) normalised by ln(window k-mers) — use
    /// <see cref="FindLowEntropyRegionsBbduk"/> for that. This method equals BBDuk with <c>entropyk=1</c> and
    /// <c>entropy = entropyThreshold / log₂(windowSize)</c> (BBDuk's 1-mer value is H_bits / log₂ w), up to BBDuk's
    /// single-precision comparison.
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Window size for analysis (default: 64). Must be ≥ 1.</param>
    /// <param name="entropyThreshold">Entropy threshold (bits, finite, ≥ 0); windows with entropy strictly below it are low complexity (default: 1.0).</param>
    /// <returns>Low-complexity regions (0-based, inclusive <c>End</c>), in ascending order and pairwise disjoint;
    /// <c>MinEntropy</c> is the lowest window entropy among the windows forming the region.
    /// Empty when the sequence is shorter than the window.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="windowSize"/> &lt; 1 or
    /// <paramref name="entropyThreshold"/> is NaN, infinite or negative.</exception>
    public static IEnumerable<LowComplexityRegion> FindLowComplexityRegions(
        DnaSequence sequence,
        int windowSize = 64,
        double entropyThreshold = 1.0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ValidateEntropyThreshold(entropyThreshold);

        return FindLowComplexityRegionsCore(sequence.Sequence, windowSize, entropyThreshold);
    }

    /// <summary>
    /// Low-complexity regions of a raw nucleotide string (same rule as the
    /// <see cref="FindLowComplexityRegions(DnaSequence, int, double)"/> overload; input is upper-cased). A window that
    /// contains an undefined symbol (anything other than A/C/G/T/U) is never flagged, exactly as BBDuk
    /// <c>maskLowEntropy</c> only tests windows with <c>ns() &lt; 1</c>; a region may still span such a symbol when
    /// flagged windows on both sides overlap it. For A/C/G/T input the result equals the <see cref="DnaSequence"/> overload.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> &lt; 1 or
    /// <paramref name="entropyThreshold"/> is NaN, infinite or negative.</exception>
    public static IEnumerable<LowComplexityRegion> FindLowComplexityRegions(
        string sequence,
        int windowSize = 64,
        double entropyThreshold = 1.0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ValidateEntropyThreshold(entropyThreshold);

        return FindLowComplexityRegionsCore(sequence.ToUpperInvariant(), windowSize, entropyThreshold);
    }

    private static void ValidateEntropyThreshold(double entropyThreshold)
    {
        if (!double.IsFinite(entropyThreshold) || entropyThreshold < 0)
            throw new ArgumentOutOfRangeException(nameof(entropyThreshold), entropyThreshold,
                "Entropy threshold must be a finite, non-negative number of bits.");
    }

    private static IEnumerable<LowComplexityRegion> FindLowComplexityRegionsCore(
        string seq,
        int windowSize,
        double entropyThreshold)
    {
        if (seq.Length < windowSize) yield break;

        int[] undefinedPrefix = UndefinedBasePrefixCounts(seq);

        // Current region = union of flagged windows [regionStart, regionEnd] (inclusive).
        int regionStart = -1;
        int regionEnd = -1;
        double minEntropy = double.MaxValue;

        for (int i = 0; i + windowSize <= seq.Length; i++)
        {
            if (undefinedPrefix[i + windowSize] != undefinedPrefix[i]) continue; // BBDuk: only windows with ns() < 1

            double entropy = CalculateShannonEntropyCore(seq.Substring(i, windowSize));
            if (!(entropy < entropyThreshold)) continue;

            int windowEnd = i + windowSize - 1;
            if (regionStart >= 0 && i <= regionEnd + 1)
            {
                // Overlaps or abuts the current region: extend the union.
                regionEnd = windowEnd;
                minEntropy = Math.Min(minEntropy, entropy);
                continue;
            }

            if (regionStart >= 0)
                yield return MakeLowComplexityRegion(seq, regionStart, regionEnd, minEntropy);

            regionStart = i;
            regionEnd = windowEnd;
            minEntropy = entropy;
        }

        if (regionStart >= 0)
            yield return MakeLowComplexityRegion(seq, regionStart, regionEnd, minEntropy);
    }

    private static LowComplexityRegion MakeLowComplexityRegion(string seq, int start, int end, double minEntropy) =>
        new(
            Start: start,
            End: end,
            Length: end - start + 1,
            MinEntropy: minEntropy,
            Sequence: seq.Substring(start, end - start + 1));

    /// <summary>BBDuk's default entropy window, <c>EntropyTracker.defaultWindowBases</c> (<c>entropywindow=50</c>).</summary>
    public const int BbdukDefaultEntropyWindow = 50;

    /// <summary>BBDuk's default entropy k-mer length, <c>EntropyTracker.defaultK</c> (<c>entropyk=5</c>).</summary>
    public const int BbdukDefaultEntropyK = 5;

    /// <summary>
    /// Returns the regions that BBTools BBDuk masks with
    /// <c>bbduk.sh entropy=<paramref name="entropyCutoff"/> entropymask=t entropywindow=<paramref name="windowSize"/>
    /// entropyk=<paramref name="k"/></c> (port of <c>BBDuk.maskLowEntropy</c> + <c>tracker/EntropyTracker</c>, BBMap 40.02).
    /// </summary>
    /// <remarks>
    /// For every window of <paramref name="windowSize"/> bases (step 1) that contains no undefined base (A/C/G/T/U in
    /// either case are defined; <c>EntropyTracker.ns() &lt; 1</c>), the normalised k-mer entropy
    /// e = −Σ p_j ln p_j / ln(W_k), p_j = c_j / W_k, W_k = windowSize − k + 1 (overlapping k-mers in the window), is
    /// computed; the window fails when e &lt; cutoff (single precision, as BBDuk's <c>float</c> comparison) and all its
    /// bases are masked. Regions are the maximal masked runs (the BBDuk bit set). A sequence shorter than the window is
    /// never masked. Note that the normaliser is ln(W_k) even when 4^k &lt; W_k (so for tiny k the maximum is below 1).
    /// The window entropy is maintained incrementally exactly as BBDuk's default <c>FAST</c> mode (running sum of the
    /// precomputed p·ln p table, k-mers containing an undefined base encoded as A), so values agree with BBDuk to the bit.
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence (any symbols; see remarks).</param>
    /// <param name="entropyCutoff">BBDuk <c>entropy=</c> value in [0, 1]; windows with entropy strictly below it are masked.</param>
    /// <param name="windowSize">BBDuk <c>entropywindow</c> (default 50); must exceed <paramref name="k"/>.</param>
    /// <param name="k">BBDuk <c>entropyk</c> (default 5), 1..15.</param>
    /// <returns>Masked regions (0-based, inclusive <c>End</c>), ascending and disjoint; <c>MinEntropy</c> is the lowest
    /// normalised entropy of the failing windows forming the region; <c>Sequence</c> is the input substring.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> outside 1..15, <paramref name="windowSize"/> ≤ k,
    /// or <paramref name="entropyCutoff"/> outside [0, 1] / NaN.</exception>
    public static IReadOnlyList<LowComplexityRegion> FindLowEntropyRegionsBbduk(
        string sequence,
        double entropyCutoff,
        int windowSize = BbdukDefaultEntropyWindow,
        int k = BbdukDefaultEntropyK)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        // EntropyTracker constructor assertions: k > 0 && k <= 15 && k < windowBases; 0 <= cutoff <= 1.
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(k, 15);
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(windowSize, k);
        if (!(entropyCutoff >= 0 && entropyCutoff <= 1))
            throw new ArgumentOutOfRangeException(nameof(entropyCutoff), entropyCutoff,
                "BBDuk entropy cutoff must lie in [0, 1].");

        return FindLowEntropyRegionsBbdukCore(sequence, (float)entropyCutoff, windowSize, k);
    }

    private static List<LowComplexityRegion> FindLowEntropyRegionsBbdukCore(string seq, float cutoff, int windowBases, int k)
    {
        var regions = new List<LowComplexityRegion>();
        if (seq.Length < windowBases) return regions; // maskLowEntropy: r.length() < window → nothing masked

        // EntropyTracker constructor.
        int windowKmers = windowBases - k + 1;
        var entropy = new double[windowKmers + 2]; // makeEntropyArray: entropy[c] = pk·ln pk, pk = c·(1/W_k)
        double mult = 1d / windowKmers;
        for (int c = 1; c < entropy.Length; c++)
        {
            double pk = c * mult;
            entropy[c] = pk * Math.Log(pk);
        }
        double entropyMult = -1 / Math.Log(windowKmers);
        int mask = ~(-1 << (2 * k));
        var counts = new Dictionary<int, int>();

        // EntropyTracker.clear().
        var ring = new char[windowBases];
        int pos = 0, pos2 = -windowBases + k - 1, len = 0, kmer = 0, kmer2 = 0, ns = 0;
        double currentEsum = 0;

        int regionStart = -1, regionEnd = -1; // inclusive
        double minEntropy = double.MaxValue;

        for (int i = 0; i < seq.Length; i++)
        {
            // EntropyTracker.add(b).
            char b = seq[i];
            char oldBase = ring[pos];
            len++;
            ring[pos] = b;
            int code = BbtoolsBaseCode(b);
            kmer = ((kmer << 2) | Math.Max(code, 0)) & mask; // symbolToNumber0: undefined → 0
            if (code < 0) ns++;
            if (len >= k)
            {
                counts.TryGetValue(kmer, out int oldCount);
                counts[kmer] = oldCount + 1;
                currentEsum = currentEsum + entropy[oldCount + 1] - entropy[oldCount];
            }
            if (pos2 >= 0)
            {
                char b2 = k > 1 ? ring[pos2] : oldBase;
                kmer2 = ((kmer2 << 2) | Math.Max(BbtoolsBaseCode(b2), 0)) & mask;
                if (len > windowBases)
                {
                    if (BbtoolsBaseCode(oldBase) < 0) ns--;
                    int oldCount = counts[kmer2];
                    counts[kmer2] = oldCount - 1;
                    currentEsum = currentEsum + entropy[oldCount - 1] - entropy[oldCount];
                }
            }
            if (++pos >= windowBases) pos = 0;
            if (++pos2 >= windowBases) pos2 = 0;

            // BBDuk.maskLowEntropy: if (i >= window−1 && ns() < 1 && !passes()) mask [leftPos, rightPos].
            if (i < windowBases - 1 || ns >= 1) continue;
            float e = (float)(currentEsum * entropyMult); // calcEntropyFast
            if (!(e > 0)) e = 0;
            if (!(e < cutoff)) continue;

            int left = len - windowBases, right = len - 1;
            if (regionStart >= 0 && left <= regionEnd + 1)
            {
                regionEnd = right;
                minEntropy = Math.Min(minEntropy, e);
                continue;
            }
            if (regionStart >= 0)
                regions.Add(MakeLowComplexityRegion(seq, regionStart, regionEnd, minEntropy));
            regionStart = left;
            regionEnd = right;
            minEntropy = e;
        }

        if (regionStart >= 0)
            regions.Add(MakeLowComplexityRegion(seq, regionStart, regionEnd, minEntropy));
        return regions;
    }

    #endregion

    #region Dust Score

    // DUST/SDUST uses overlapping nucleotide triplets (3-mers) as the word.
    // k = 3 is hardcoded in the reference implementations; Morgulis et al. (2006),
    // J Comput Biol 13(5):1028-1040, doi:10.1089/cmb.2006.13.1028; NCBI dustmasker
    // (symdust.cpp); lh3/sdust (SD_WLEN = 3).
    private const int DustWordSize = 3;

    // Mask threshold for the DUST score: 2.0, i.e. the reference default level T = 20
    // (score(x) > T/10). NCBI symdust.hpp DEFAULT_LEVEL = 20; lh3/sdust "int T = 20".
    private const double DustMaskThreshold = 2.0;

    // Default SDUST window length (bases): NCBI symdust DEFAULT_WINDOW = 64; lh3/sdust "int W = 64".
    private const int DustWindowSize = 64;

    // Number of distinct triplet codes (4^3) and the 2-bit rolling mask (lh3/sdust SD_WTOT / SD_WMSK).
    private const int DustTripletCodes = 1 << (DustWordSize << 1);
    private const int DustTripletMask = DustTripletCodes - 1;

    /// <summary>
    /// Calculates the DUST low-complexity score of a sequence (Morgulis et al. 2006).
    /// For a sequence with ℓ overlapping triplets, where triplet t occurs c_t times,
    /// score = Σ_t c_t·(c_t−1)/2 / (ℓ − 1). This is the score thresholded by NCBI
    /// <c>dustmasker</c> (symdust: <c>10·r &gt; level·(ℓ−1)</c>) and by lh3/sdust
    /// (<c>new_r·10 &gt; T·new_l</c> with <c>new_l</c> = ℓ − 1). A HIGHER score indicates
    /// LOWER complexity; all-distinct words give 0; a homopolymer of length L scores (L−2)/2.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="wordSize">Word size; must be 3. The DUST score is defined for triplets only
    /// (Morgulis et al. 2006; NCBI symdust <c>triplet_type</c>; lh3/sdust <c>SD_WLEN = 3</c>). The
    /// parameter is kept for source compatibility. For a sourced k-mer generalisation use
    /// <see cref="CalculateLongdustScore(string, int, double?)"/> / <see cref="FindLongdustRegions(string, int, int, double, int, int, bool, bool, double?)"/>
    /// (Li &amp; Li 2025, longdust).</param>
    /// <returns>DUST score (≥ 0); 0 when fewer than two words exist (ℓ − 1 ≤ 0).</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="wordSize"/> ≠ 3.</exception>
    public static double CalculateDustScore(DnaSequence sequence, int wordSize = DustWordSize)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateDustWordSize(wordSize);
        return CalculateDustScoreCore(sequence.Sequence, wordSize);
    }

    // DUST is defined for triplets only; any other word size was an unsourced extrapolation (B04 F34).
    private static void ValidateDustWordSize(int wordSize)
    {
        if (wordSize != DustWordSize)
            throw new ArgumentOutOfRangeException(nameof(wordSize), wordSize,
                "The DUST score is defined for triplets only (word size 3; Morgulis et al. 2006, NCBI symdust, lh3/sdust). " +
                "Use CalculateLongdustScore / FindLongdustRegions (longdust, Li & Li 2025) for other k-mer lengths.");
    }

    /// <summary>
    /// Calculates the DUST low-complexity score from a raw sequence string
    /// (score = Σ_t c_t·(c_t−1)/2 / (ℓ − 1), see <see cref="CalculateDustScore(DnaSequence, int)"/>).
    /// The string is upper-cased to match the normalization applied by <see cref="DnaSequence"/>.
    /// </summary>
    /// <param name="sequence">Raw sequence string; null or empty yields 0.</param>
    /// <param name="wordSize">Word size; must be 3 (DUST is defined for triplets only).</param>
    /// <returns>DUST score (≥ 0); 0 when null/empty or fewer than two words exist.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="wordSize"/> ≠ 3.</exception>
    public static double CalculateDustScore(string sequence, int wordSize = DustWordSize)
    {
        ValidateDustWordSize(wordSize);
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateDustScoreCore(sequence.ToUpperInvariant(), wordSize);
    }

    private static double CalculateDustScoreCore(string seq, int wordSize)
    {
        // ℓ = number of overlapping words (L − 2 for triplets). The DUST normalisation is
        // ℓ − 1 (Morgulis et al. 2006; NCBI symdust thresholds_[ℓ−1] = (ℓ−1)·level;
        // lh3/sdust new_l = kdq_size − i − 1). With ℓ ≤ 1 no pair of words exists: score 0.
        int wordCount = seq.Length - wordSize + 1;
        if (wordCount < 2) return 0;

        // Word tally via the canonical k-mer counter (KMER-COUNT-001).
        var wordCounts = KmerAnalyzer.CountKmers(seq, wordSize);

        // Numerator Σ_t c_t·(c_t−1)/2 (= lh3/sdust's running "rw += cw[t]++").
        // Promote to double before multiplying: count·(count−1) overflows Int32 for L ≳ 4.6·10⁴.
        double sum = 0;
        foreach (int count in wordCounts.Values)
            sum += (double)count * (count - 1) / 2.0;

        return sum / (wordCount - 1);
    }

    /// <summary>
    /// Masks low-complexity regions with the symmetric DUST (SDUST) algorithm of
    /// Morgulis et al. (2006), as implemented by NCBI <c>dustmasker</c> and lh3/sdust.
    /// Every <em>perfect interval</em> — a subsequence x of at most <paramref name="windowSize"/>
    /// bases whose DUST score exceeds <paramref name="threshold"/> and is not exceeded by the
    /// score of any of its sub-intervals — found inside any window is masked; overlapping or
    /// adjacent masked intervals are merged. The result is symmetric and context-insensitive.
    /// This is a line-by-line port of lh3/sdust <c>sdust_core</c> (which reproduces dustmasker
    /// output); non-ACGT characters break the input into independently scanned pieces.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">SDUST window length W in bases (default: 64; must be ≥ 3).</param>
    /// <param name="threshold">DUST score threshold (default: 2.0 = level 20); an interval is
    /// low-complexity when its score is strictly greater than this value. Must be ≥ 0.</param>
    /// <param name="maskChar">Character to use for masking (default: 'N').</param>
    /// <returns>Masked sequence (same length as the input).</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="windowSize"/> &lt; 3
    /// or <paramref name="threshold"/> is negative, NaN or infinite.</exception>
    public static string MaskLowComplexity(
        DnaSequence sequence,
        int windowSize = DustWindowSize,
        double threshold = DustMaskThreshold,
        char maskChar = 'N')
        => MaskLowComplexity(sequence, windowSize, threshold, maskChar, DustDefaultLinker, softMask: false);

    /// <summary>
    /// SDUST masking (see <see cref="MaskLowComplexity(DnaSequence, int, double, char)"/>) with the
    /// NCBI <c>dustmasker</c> interval <paramref name="linker"/> and optional soft masking.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">SDUST window length W in bases (must be ≥ 3).</param>
    /// <param name="threshold">DUST score threshold (strictly greater ⇒ low complexity; ≥ 0).</param>
    /// <param name="maskChar">Hard-mask character (ignored when <paramref name="softMask"/> is true).</param>
    /// <param name="linker">dustmasker <c>-linker</c>: consecutive masked intervals are merged when
    /// the number of unmasked bases between them is &lt; <paramref name="linker"/> (NCBI symdust
    /// <c>save_masked_regions</c>: <c>prev.last + linker ≥ next.first</c>, closed coordinates).
    /// 1 (dustmasker's default) merges only overlapping/adjacent intervals, which is exactly the
    /// lh3/sdust behaviour. Must be 1–32, the range symdust accepts (dustmasker silently substitutes the
    /// default 1 for any other value; here that is rejected instead).</param>
    /// <param name="softMask">When true, masked bases are written in lower case and all other bases in
    /// upper case (dustmasker <c>-outfmt fasta</c>); when false, masked bases become <paramref name="maskChar"/>.</param>
    /// <param name="engine">DUST implementation (default <see cref="DustEngine.Sdust"/>, the lh3/sdust port).
    /// <see cref="DustEngine.Dustmasker"/> reproduces NCBI dustmasker 2.12.0 exactly (symdust core +
    /// <c>GetDustMasks_SkipNs</c>): IUPAC codes are scanned as bases, only N runs longer than the window
    /// (and leading/trailing N runs) cut the scan and are reported as masked; window 8–64, 10·threshold an
    /// integer 2–64, IUPAC DNA input only (B04 F53).</param>
    /// <returns>Masked sequence (same length as the input).</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown on an invalid window, threshold or linker (or one outside the dustmasker ranges with <see cref="DustEngine.Dustmasker"/>).</exception>
    /// <exception cref="ArgumentException">Thrown with <see cref="DustEngine.Dustmasker"/> when the input is not IUPAC DNA.</exception>
    public static string MaskLowComplexity(
        DnaSequence sequence,
        int windowSize,
        double threshold,
        char maskChar,
        int linker,
        bool softMask = false,
        DustEngine engine = DustEngine.Sdust)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return MaskLowComplexity(sequence.Sequence, windowSize, threshold, maskChar, linker, softMask, engine);
    }

    /// <summary>
    /// SDUST masking of a raw nucleotide string that may contain N, IUPAC codes or other
    /// non-ACGT symbols. As lh3/sdust specifies ("N effectively breaks input into pieces of
    /// independent sequences"), every non-ACGT character ends the current window: each maximal
    /// ACGT run is scanned exactly as a separate sdust input (output = sdust run per piece, shifted
    /// to input coordinates) and a non-ACGT symbol is never part of a perfect interval (only a <paramref name="linker"/> &gt; 1 can join two
    /// intervals across one). The output is upper-cased (the <see cref="DnaSequence"/> normalisation);
    /// with <paramref name="softMask"/> masked bases are lower case instead (dustmasker <c>-outfmt fasta</c>).
    /// </summary>
    /// <param name="sequence">Nucleotide string (any symbols; ACGT case-insensitive).</param>
    /// <param name="windowSize">SDUST window length W (default 64; ≥ 3).</param>
    /// <param name="threshold">DUST score threshold (default 2.0 = level 20; ≥ 0).</param>
    /// <param name="maskChar">Hard-mask character (default 'N'; ignored when soft-masking).</param>
    /// <param name="linker">dustmasker linker (default 1 = sdust/dustmasker default; 1–32).</param>
    /// <param name="softMask">Lower-case masking instead of <paramref name="maskChar"/> (default false).</param>
    /// <param name="engine">DUST implementation (default <see cref="DustEngine.Sdust"/>, the lh3/sdust port).
    /// <see cref="DustEngine.Dustmasker"/> reproduces NCBI dustmasker 2.12.0 exactly (symdust core +
    /// <c>GetDustMasks_SkipNs</c>): IUPAC codes are scanned as bases, only N runs longer than the window
    /// (and leading/trailing N runs) cut the scan and are reported as masked; window 8–64, 10·threshold an
    /// integer 2–64, IUPAC DNA input only (B04 F53).</param>
    /// <returns>Masked sequence (same length as the input).</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown on an invalid window, threshold or linker (or one outside the dustmasker ranges with <see cref="DustEngine.Dustmasker"/>).</exception>
    /// <exception cref="ArgumentException">Thrown with <see cref="DustEngine.Dustmasker"/> when the input is not IUPAC DNA.</exception>
    public static string MaskLowComplexity(
        string sequence,
        int windowSize = DustWindowSize,
        double threshold = DustMaskThreshold,
        char maskChar = 'N',
        int linker = DustDefaultLinker,
        bool softMask = false,
        DustEngine engine = DustEngine.Sdust)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        string upper = sequence.ToUpperInvariant();
        return ApplySdustMask(upper, FindDustIntervals(upper, windowSize, threshold, linker, engine), maskChar, softMask);
    }

    /// <summary>
    /// Returns the SDUST low-complexity intervals as 0-based half-open [Start, End) pairs in
    /// ascending order — the native output of lh3/sdust (<c>name start end</c>) and, with
    /// End − 1, of dustmasker <c>-outfmt interval</c> (closed coordinates). Non-ACGT symbols
    /// split the input into independently scanned ACGT runs (sdust's contract).
    /// </summary>
    /// <param name="sequence">Nucleotide string (any symbols; ACGT case-insensitive).</param>
    /// <param name="windowSize">SDUST window length W (default 64; ≥ 3).</param>
    /// <param name="threshold">DUST score threshold (default 2.0; ≥ 0).</param>
    /// <param name="linker">dustmasker linker (default 1; 1–32), see
    /// <see cref="MaskLowComplexity(string, int, double, char, int, bool, DustEngine)"/>.</param>
    /// <param name="engine">DUST implementation (default <see cref="DustEngine.Sdust"/>, the lh3/sdust port).
    /// <see cref="DustEngine.Dustmasker"/> reproduces NCBI dustmasker 2.12.0 exactly (symdust core +
    /// <c>GetDustMasks_SkipNs</c>): IUPAC codes are scanned as bases, only N runs longer than the window
    /// (and leading/trailing N runs) cut the scan and are reported as masked; window 8–64, 10·threshold an
    /// integer 2–64, IUPAC DNA input only (B04 F53).</param>
    /// <returns>Merged masked intervals.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown on an invalid window, threshold or linker (or one outside the dustmasker ranges with <see cref="DustEngine.Dustmasker"/>).</exception>
    /// <exception cref="ArgumentException">Thrown with <see cref="DustEngine.Dustmasker"/> when the input is not IUPAC DNA.</exception>
    public static IReadOnlyList<(int Start, int End)> FindLowComplexityIntervals(
        string sequence,
        int windowSize = DustWindowSize,
        double threshold = DustMaskThreshold,
        int linker = DustDefaultLinker,
        DustEngine engine = DustEngine.Sdust)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindDustIntervals(sequence, windowSize, threshold, linker, engine);
    }

    /// <summary>
    /// <see cref="FindLowComplexityIntervals(string, int, double, int, DustEngine)"/> for a <see cref="DnaSequence"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Window length W.</param>
    /// <param name="threshold">DUST score threshold.</param>
    /// <param name="linker">dustmasker linker (1–32).</param>
    /// <param name="engine">DUST implementation (default sdust).</param>
    /// <returns>Merged masked intervals, half-open.</returns>
    public static IReadOnlyList<(int Start, int End)> FindLowComplexityIntervals(
        DnaSequence sequence,
        int windowSize = DustWindowSize,
        double threshold = DustMaskThreshold,
        int linker = DustDefaultLinker,
        DustEngine engine = DustEngine.Sdust)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindLowComplexityIntervals(sequence.Sequence, windowSize, threshold, linker, engine);
    }

    // Engine dispatch shared by the mask and interval entry points (validation per engine).
    private static List<(int Start, int End)> FindDustIntervals(string seq, int windowSize, double threshold, int linker, DustEngine engine)
    {
        ValidateSdustParameters(windowSize, threshold, linker);
        switch (engine)
        {
            case DustEngine.Sdust:
                return FindSdustIntervals(seq, windowSize, threshold, linker);
            case DustEngine.Dustmasker:
                string upper = seq.ToUpperInvariant();
                int level = ValidateDustmaskerParameters(upper, windowSize, threshold);
                return FindDustmaskerIntervals(upper, windowSize, level, linker);
            default:
                throw new ArgumentOutOfRangeException(nameof(engine), engine, "Unknown DUST engine.");
        }
    }

    // dustmasker DEFAULT_LINKER = 1 (symdust.hpp); identical to lh3/sdust's adjacency merge.
    private const int DustDefaultLinker = 1;

    // symdust constructor: linker_( (linker >= 1 && linker <= 32) ? linker : DEFAULT_LINKER ).
    private const int DustMaxLinker = 32;

    private static void ValidateSdustParameters(int windowSize, double threshold, int linker)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, DustWordSize);
        if (double.IsNaN(threshold) || double.IsInfinity(threshold) || threshold < 0)
            throw new ArgumentOutOfRangeException(nameof(threshold), threshold, "Threshold must be a finite value ≥ 0.");
        ArgumentOutOfRangeException.ThrowIfLessThan(linker, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(linker, DustMaxLinker);
    }

    private static string ApplySdustMask(string seq, List<(int Start, int End)> intervals, char maskChar, bool softMask)
    {
        if (intervals.Count == 0) return seq;

        var masked = seq.ToCharArray();
        foreach (var (start, end) in intervals)
        {
            if (softMask)
            {
                for (int i = start; i < end; i++) masked[i] = char.ToLowerInvariant(masked[i]);
            }
            else
            {
                Array.Fill(masked, maskChar, start, end - start);
            }
        }

        return new string(masked);
    }

    // Perfect interval [Start, Finish) with raw score R = Σ c(c−1)/2 and normaliser L = ℓ − 1.
    private struct SdustPerfectInterval
    {
        public int Start, Finish, R, L;
    }

    /// <summary>
    /// SDUST core: returns the merged, 0-based half-open masked intervals [start, end).
    /// Port of lh3/sdust <c>sdust_core</c>/<c>shift_window</c>/<c>find_perfect</c>/
    /// <c>save_masked_regions</c>; integer tests <c>x·10 &gt; T·y</c> become <c>x &gt; threshold·y</c>.
    /// The merge of consecutive intervals uses the NCBI symdust linker rule (linker 1 = sdust).
    /// A non-ACGT symbol resets the window, so every maximal ACGT run is scanned as an independent
    /// sequence (sdust's documented intent; upstream sdust_core leaks the window across the break).
    /// </summary>
    private static List<(int Start, int End)> FindSdustIntervals(string seq, int windowSize, double threshold, int linker)
    {
        var res = new List<(int Start, int End)>();
        var perfect = new List<SdustPerfectInterval>(); // descending start, then ascending finish
        var window = new SdustWindow(windowSize - DustWordSize + 1);
        var cw = new int[DustTripletCodes];
        var cv = new int[DustTripletCodes];
        var scratch = new int[DustTripletCodes];
        int rw = 0, rv = 0, suffixLen = 0;
        int l = 0, t = 0; // l = length of the current contiguous ACGT run; t = current triplet code

        for (int i = 0; i <= seq.Length; i++)
        {
            int b = i < seq.Length ? AcgtCode(seq[i]) : -1;
            if (b >= 0)
            {
                ++l;
                t = ((t << 2) | b) & DustTripletMask;
                if (l >= DustWordSize)
                {
                    int start = Math.Max(l - windowSize, 0) + (i + 1 - l);
                    SdustSaveMaskedRegions(res, perfect, start, linker);
                    SdustShiftWindow(t, window, threshold, ref suffixLen, ref rw, ref rv, cw, cv);
                    if (rw > threshold * suffixLen)
                        SdustFindPerfect(perfect, window, threshold, start, suffixLen, rv, cv, scratch);
                }
            }
            else
            {
                // Non-ACGT or end of input: flush all pending perfect intervals.
                int start = Math.Max(l - windowSize + 1, 0) + (i + 1 - l);
                while (perfect.Count > 0) SdustSaveMaskedRegions(res, perfect, start++, linker);
                l = t = 0;

                // sdust's stated contract is that "N effectively breaks input into pieces of independent
                // sequences", but sdust_core resets only l and t: the triplet window, its counts and the
                // suffix state leak across the N, so later intervals get shifted coordinates (they can
                // even end past the sequence end, e.g. 35–72 on a 53-bp input). Resetting the window
                // state makes each ACGT run exactly an independent sdust run (B04 F36).
                window.Clear();
                Array.Clear(cw);
                Array.Clear(cv);
                rw = rv = suffixLen = 0;
            }
        }

        return res;
    }

    /// <summary>
    /// 2-bit nucleotide code shared by SDUST (this class) and <see cref="RepeatFinder"/> (direct repeats, TRF):
    /// A/C/G/T (either case) → 0..3, any other symbol (N, IUPAC, U, gap, …) → −1. Same mapping as sdust's
    /// <c>seq_nt4_table</c> (non-ACGT = 4 there, i.e. "not a nucleotide").
    /// </summary>
    internal static int AcgtCode(char c) => c switch
    {
        'A' or 'a' => 0,
        'C' or 'c' => 1,
        'G' or 'g' => 2,
        'T' or 't' => 3,
        _ => -1,
    };

    private static void SdustShiftWindow(
        int t, SdustWindow w, double threshold, ref int suffixLen, ref int rw, ref int rv, int[] cw, int[] cv)
    {
        if (w.Count >= w.Capacity)
        {
            int s = w.Shift();
            rw -= --cw[s];
            if (suffixLen > w.Count)
            {
                --suffixLen;
                rv -= --cv[s];
            }
        }

        w.Push(t);
        ++suffixLen;
        rw += cw[t]++;
        rv += cv[t]++;
        if (cv[t] > 2 * threshold)
        {
            int s;
            do
            {
                s = w[w.Count - suffixLen];
                rv -= --cv[s];
                --suffixLen;
            } while (s != t);
        }
    }

    private static void SdustSaveMaskedRegions(
        List<(int Start, int End)> res, List<SdustPerfectInterval> perfect, int start, int linker)
    {
        if (perfect.Count == 0 || perfect[^1].Start >= start) return;

        var p = perfect[^1];
        bool saved = false;
        if (res.Count > 0)
        {
            var (s, f) = res[^1];
            // NCBI symdust: merge when prev.last + linker >= next.first (closed coordinates), i.e.
            // p.Start <= f + linker − 1 with f half-open; linker = 1 is lh3/sdust's "p->start <= f".
            if (p.Start <= f + linker - 1)
            {
                saved = true;
                res[^1] = (s, Math.Max(f, p.Finish));
            }
        }
        if (!saved) res.Add((p.Start, p.Finish));

        int k = perfect.Count - 1;
        while (k >= 0 && perfect[k].Start < start) --k; // drop intervals that fell out of the window
        perfect.RemoveRange(k + 1, perfect.Count - (k + 1));
    }

    private static void SdustFindPerfect(
        List<SdustPerfectInterval> perfect, SdustWindow w, double threshold, int start, int suffixLen, int rv, int[] cv, int[] c)
    {
        Array.Copy(cv, c, cv.Length);
        int r = rv, maxR = 0, maxL = 0;
        for (int i = w.Count - suffixLen - 1; i >= 0; --i)
        {
            int t = w[i];
            r += c[t]++;
            int newR = r, newL = w.Count - i - 1; // newL = ℓ − 1 for the interval's ℓ triplets
            if (newR > threshold * newL)
            {
                int j;
                for (j = 0; j < perfect.Count && perfect[j].Start >= i + start; ++j)
                {
                    var p = perfect[j];
                    if (maxR == 0 || (long)p.R * maxL > (long)maxR * p.L)
                    {
                        maxR = p.R;
                        maxL = p.L;
                    }
                }
                if (maxR == 0 || (long)newR * maxL >= (long)maxR * newL)
                {
                    maxR = newR;
                    maxL = newL;
                    perfect.Insert(j, new SdustPerfectInterval
                    {
                        Start = i + start,
                        Finish = w.Count + (DustWordSize - 1) + start,
                        R = newR,
                        L = newL,
                    });
                }
            }
        }
    }

    /// <summary>Fixed-capacity FIFO of triplet codes with random access (lh3 kdq_t(int)).</summary>
    private sealed class SdustWindow
    {
        private readonly int[] _buf;
        private int _front;

        public SdustWindow(int capacity)
        {
            _buf = new int[capacity];
        }

        public int Capacity => _buf.Length;
        public int Count { get; private set; }

        public int this[int index] => _buf[(_front + index) % _buf.Length];

        public void Push(int value)
        {
            _buf[(_front + Count) % _buf.Length] = value;
            Count++;
        }

        public int Shift()
        {
            int value = _buf[_front];
            _front = (_front + 1) % _buf.Length;
            Count--;
            return value;
        }

        public void Clear()
        {
            _front = 0;
            Count = 0;
        }
    }

    #region Longdust (k-mer generalisation of DUST; Li & Li 2025)

    // Defaults of longdust 1.4-r97 ld_opt_init (github.com/lh3/longdust, MIT):
    // kmer = 7, ws = 5000, thres = 0.6, xdrop_len = 50, min_start_cnt = 3, approx = 0, gc disabled.
    private const int LongdustDefaultK = 7;
    private const int LongdustDefaultWindow = 5000;
    private const double LongdustDefaultThreshold = 0.6;
    private const int LongdustDefaultXdrop = 50;
    private const int LongdustDefaultMinStartCount = 3;
    private const int LongdustMaxK = 14;          // longdust: assert(k < LD_MAX_K = 15)
    private const int LongdustMaxWindow = 0xfffe; // longdust: assert(ws < 0xffff) (16-bit counts)

    /// <summary>
    /// Longdust complexity score of a whole sequence x (Li &amp; Li 2025, arXiv:2509.07357; lh3/longdust):
    /// S_L(x) = Σ_t log c_x(t)! − f(ℓ(x)/4^k), with c_x(t) the count of k-mer t in x,
    /// ℓ(x) = |x| − k + 1 and f(λ) = 4^k e^{−λ} Σ_n log(n!) λ^n/n! (the expected Σ log c! of a
    /// random sequence). This is the sourced k-mer generalisation of the DUST score: DUST's
    /// Σ c(c−1)/2 over triplets is replaced by the composite-likelihood term Σ log c!.
    /// Higher = lower complexity; longdust calls x low-complexity when S_L(x) − T·ℓ(x) &gt; 0 (T = 0.6).
    /// </summary>
    /// <param name="sequence">Nucleotide string; k-mers containing a non-ACGT symbol are not counted
    /// but still count in ℓ, as in longdust's window.</param>
    /// <param name="k">k-mer length (default 7; 1–14).</param>
    /// <param name="gcContent">Optional genome GC fraction in (0, 1) for longdust's GC correction
    /// (<c>-g</c>); null = uniform base composition (longdust default).</param>
    /// <returns>S_L(x); 0 when the sequence holds no k-mer position.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> or <paramref name="gcContent"/> is out of range.</exception>
    public static double CalculateLongdustScore(string sequence, int k = LongdustDefaultK, double? gcContent = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateLongdustK(k);
        ValidateLongdustGc(gcContent);
        int positions = sequence.Length - k + 1;
        if (positions <= 0) return 0;

        var f = LongdustF(k, positions, gcContent);
        var counts = new Dictionary<int, int>();
        int mask = (1 << (2 * k)) - 1, x = 0, run = 0;
        double s = 0;
        foreach (char ch in sequence)
        {
            int b = AcgtCode(ch);
            if (b < 0) { run = 0; continue; }
            x = ((x << 2) | b) & mask;
            if (++run < k) continue;
            counts.TryGetValue(x, out int c);
            counts[x] = ++c;
            if (c >= 2) s += Math.Log(c); // Σ_t log c_t! accumulated as log 2 + … + log c_t
        }
        return s - f[positions];
    }

    /// <summary>
    /// Finds low-complexity regions with longdust (Li &amp; Li 2025, arXiv:2509.07357), the k-mer
    /// generalisation of SDUST for long windows (STRs, VNTRs, satellites). A port of lh3/longdust
    /// 1.4-r97 (<c>ld_dust1</c>/<c>ld_dust2</c>, MIT licence): at each position a backward then forward
    /// scan over the last <paramref name="windowSize"/> k-mers finds a good interval with
    /// S_L(x) − T·ℓ(x) &gt; 0 (see <see cref="CalculateLongdustScore"/>), with X-drop and the reference's
    /// speed-ups; overlapping hits are merged; by default the union over both strands is returned.
    /// Non-ACGT symbols make the overlapping k-mers ambiguous (they score −T).
    /// </summary>
    /// <param name="sequence">Nucleotide string (any symbols; ACGT case-insensitive).</param>
    /// <param name="k">k-mer length (<c>-k</c>, default 7; 1–14).</param>
    /// <param name="windowSize">Window size in k-mers (<c>-w</c>, default 5000; 1–65534).</param>
    /// <param name="threshold">Score threshold T per k-mer (<c>-t</c>, default 0.6; finite, &gt; 0).</param>
    /// <param name="xdropLength">X-drop length (<c>-e</c>, default 50; 0 disables X-drop, i.e. uses the window).</param>
    /// <param name="minStartCount">Minimum count of the current k-mer in the window before a search starts
    /// (<c>-b</c>, default 3; ≥ 2).</param>
    /// <param name="forwardOnly">Scan the forward strand only (<c>-f</c>); default false = union of both strands.</param>
    /// <param name="approximate">Guaranteed O(Lw) mode with one forward pass (<c>-a</c>).</param>
    /// <param name="gcContent">Genome GC fraction in (0, 1) for GC correction (<c>-g</c>); null = off.</param>
    /// <returns>0-based half-open [Start, End) intervals in ascending order (longdust BED output).</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when a parameter is out of range.</exception>
    public static IReadOnlyList<(int Start, int End)> FindLongdustRegions(
        string sequence,
        int k = LongdustDefaultK,
        int windowSize = LongdustDefaultWindow,
        double threshold = LongdustDefaultThreshold,
        int xdropLength = LongdustDefaultXdrop,
        int minStartCount = LongdustDefaultMinStartCount,
        bool forwardOnly = false,
        bool approximate = false,
        double? gcContent = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateLongdustK(k);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(windowSize, LongdustMaxWindow);
        if (double.IsNaN(threshold) || double.IsInfinity(threshold) || threshold <= 0)
            throw new ArgumentOutOfRangeException(nameof(threshold), threshold, "Threshold must be a finite value > 0.");
        ArgumentOutOfRangeException.ThrowIfNegative(xdropLength);
        ArgumentOutOfRangeException.ThrowIfLessThan(minStartCount, 2);
        ValidateLongdustGc(gcContent);

        var ld = new LongdustScanner(k, windowSize, threshold, xdropLength, minStartCount, approximate, gcContent);
        int n = sequence.Length;
        var fwd = new int[n];
        for (int i = 0; i < n; i++) fwd[i] = AcgtCode(sequence[i]);
        var forward = ld.Dust1(fwd);
        if (forwardOnly) return forward;

        // ld_dust2: reverse complement, map back, merge the two sorted interval lists.
        var rev = new int[n];
        for (int i = 0; i < n; i++) rev[n - 1 - i] = fwd[i] < 0 ? -1 : 3 - fwd[i];
        var revRaw = ld.Dust1(rev);
        var reverse = new List<(int Start, int End)>(revRaw.Count);
        for (int i = revRaw.Count - 1; i >= 0; i--) reverse.Add((n - revRaw[i].End, n - revRaw[i].Start));

        var merged = new List<(int Start, int End)>();
        int st = 0, en = 0, j0 = 0, j1 = 0;
        while (j0 < forward.Count || j1 < reverse.Count)
        {
            bool takeReverse = j0 >= forward.Count || (j1 < reverse.Count && forward[j0].Start >= reverse[j1].Start);
            var p = takeReverse ? reverse[j1++] : forward[j0++];
            if (p.Start <= en)
            {
                en = Math.Max(en, p.End);
            }
            else
            {
                if (en > st) merged.Add((st, en));
                st = p.Start;
                en = p.End;
            }
        }
        if (en > st) merged.Add((st, en));
        return merged;
    }

    private static void ValidateLongdustK(int k)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(k, LongdustMaxK);
    }

    private static void ValidateLongdustGc(double? gc)
    {
        if (gc is { } g && !(g > 0.0 && g < 1.0))
            throw new ArgumentOutOfRangeException(nameof(gc), g, "GC content must lie in (0, 1).");
    }

    // f[l] for l = 0..maxL (longdust ld_cal_f / ld_cal_f2): f(λ) = 4^k e^{−λ} Σ_{n≥2} log(n!) λ^n/n!,
    // λ = l/4^k (per GC class when GC correction is on); Stirling form for λ ≥ 30.
    private static double[] LongdustF(int k, int maxL, double? gc)
    {
        double[] dr;
        int[] nDr;
        uint nKmer = 1U << (2 * k);
        if (gc is { } g)
        {
            dr = new double[k + 1];
            nDr = new int[k + 1];
            for (int i = 0; i <= k; ++i)
                dr[i] = Math.Pow(g / 0.5, i) * Math.Pow((1.0 - g) / 0.5, k - i);
            for (uint x = 0; x < nKmer; ++x)
            {
                int nGc = 0;
                for (int i = 0; i < k; ++i)
                {
                    uint b = (x >> (2 * i)) & 3;
                    if (b == 1 || b == 2) ++nGc;
                }
                nDr[nGc]++;
            }
        }
        else
        {
            dr = new[] { 1.0 };
            nDr = new[] { (int)nKmer };
        }

        const double eps = 1e-9;
        const int maxN = 10000;
        var f = new double[maxL + 1];
        for (int l = 1; l <= maxL; ++l)
        {
            for (int i = 0; i < dr.Length; ++i)
            {
                double lambda = (double)l / nKmer * dr[i];
                double fli;
                if (lambda < 30.0)
                {
                    double x = 0.0, sn = 0.0, y = lambda;
                    for (int n = 2; n <= maxN; ++n)
                    {
                        sn += Math.Log(n);
                        y *= lambda / n;
                        double z = y * sn;
                        if (z < x * eps) break;
                        x += z;
                    }
                    fli = x * Math.Exp(-lambda);
                }
                else
                {
                    // Stirling series (longdust ld_f_large).
                    double xl = 0.5 * Math.Log(2.0 * Math.PI * Math.E * lambda)
                        - 1.0 / 12.0 / lambda * (1.0 + 0.5 / lambda + 19.0 / 30.0 / lambda / lambda);
                    fli = xl + lambda * (Math.Log(lambda) - 1.0);
                }
                f[l] += fli * nDr[i];
            }
        }
        return f;
    }

    /// <summary>Port of longdust's <c>ld_data_t</c> and the forward-strand scan <c>ld_dust1</c>.</summary>
    private sealed class LongdustScanner
    {
        private readonly int _k, _ws, _xdropLen, _minStartCnt, _maxTest;
        private readonly double _thres;
        private readonly bool _approx;
        private readonly double[] _f, _c;
        private readonly ushort[] _ht;          // ld->ht: counts of the current backward/forward scan
        private readonly SdustWindow _q;        // entries x<<1 | ambiguous
        private readonly int[] _forPos;
        private readonly double[] _forMax;

        public LongdustScanner(int k, int ws, double thres, int xdropLen, int minStartCnt, bool approx, double? gc)
        {
            _k = k; _ws = ws; _thres = thres; _xdropLen = xdropLen; _minStartCnt = minStartCnt; _approx = approx;
            _ht = new ushort[1 << (2 * k)];
            _f = LongdustF(k, ws + 1, gc);
            // c[i] = log i (c[0] = c[1] = 0). One extra slot: ld_extend may read c[count + 1] with count = ws.
            _c = new double[ws + 2];
            for (int i = 2; i <= ws + 1; ++i) _c[i] = Math.Log(i);
            _q = new SdustWindow(ws);
            _forPos = new int[ws + 1];
            _forMax = new double[ws + 1];

            // max_test: the maximum step used by IfBackward (i − 1 + k = minimum detectable homopolymer).
            int m;
            double s = 0.0;
            for (m = 1; m < ws; ++m)
            {
                s += _c[m] - thres;
                double sl = s - _f[m];
                if (sl > 0.0) break;
            }
            _maxTest = (int)(m * Math.Log(m) / thres);
        }

        private int Forward(int i0, double maxBack)
        {
            Array.Clear(_ht);
            int maxI = -1, qn = _q.Count;
            double s = 0.0, maxSf = 0.0;
            for (int i = i0, l = 1; i < qn; ++i, ++l)
            {
                int x = _q[i];
                s += ((x & 1) != 0 ? 0 : _c[++_ht[x >> 1]]) - _thres;
                double sl = s - _f[l];
                if (sl >= maxSf) { maxSf = sl; maxI = i; }
                if (sl > maxBack + 1e-6) break;
            }
            return maxI;
        }

        private int Backward(ushort[] winHt)
        {
            double xdrop = _thres * (_xdropLen > 0 ? _xdropLen : _ws);
            int maxI = -1, qn = _q.Count;
            double s = 0.0, sw = 0.0, maxSb = 0.0, lastSl = -1.0;

            Array.Clear(_ht);
            int nForPos = 0;
            for (int i = qn - 1, l = 1; i >= 0; --i, ++l)
            {
                int x = _q[i];
                s += ((x & 1) != 0 ? 0 : _c[++_ht[x >> 1]]) - _thres;
                double sl = s - _f[l];
                sw += ((x & 1) != 0 ? 0 : _c[winHt[x >> 1] + 1 - _ht[x >> 1]]) - _thres;
                if (sw - _f[l] < 0.0) break; // the forward pass cannot reach the current position
                if (sl < lastSl && lastSl > 0.0 && lastSl == maxSb)
                {
                    _forPos[nForPos] = i + 1;
                    _forMax[nForPos++] = maxSb;
                }
                if (sl >= maxSb)
                {
                    maxSb = sl;
                    maxI = i;
                }
                else if (maxI >= 0 && maxSb - sl > xdrop)
                {
                    break; // X-drop
                }
                lastSl = sl;
            }
            if (maxI < 0) return -1;
            if (nForPos == 0 || maxI < _forPos[nForPos - 1])
            {
                _forPos[nForPos] = maxI;
                _forMax[nForPos++] = maxSb;
            }
            for (int i = nForPos - 1, maxEnd = -1; i >= 0; --i)
            {
                if (_forPos[i] < maxEnd) continue;
                int e = Forward(_forPos[i], _forMax[i]);
                if (e == qn - 1) return _forPos[i];
                if (_approx) break; // approximate mode: one forward pass only
                maxEnd = Math.Max(maxEnd, e);
            }
            return -1;
        }

        private int Extend()
        {
            int x = _q[_q.Count - 1];
            int l = _q.Count - 1;
            if ((x & 1) != 0) return -1;
            double diff = _c[_ht[x >> 1] + 1] - (_f[l + 1] - _f[l]);
            if (diff < _thres) return -1; // extending would not increase the score
            ++_ht[x >> 1];
            return 0;
        }

        private bool IfBackward(ushort[] winHt)
        {
            double s = 0.0;
            for (int i = _q.Count - 1, j = 0; i >= 0 && j < _maxTest; --i, ++j)
            {
                int x = _q[i];
                s += ((x & 1) != 0 ? 0 : _c[winHt[x >> 1]]) - _thres;
                if (s < 0.0) return false;
            }
            return true;
        }

        private static void SaveInterval(List<(int Start, int End)> intv, int st, int en)
        {
            int k;
            for (k = intv.Count - 1; k >= 0; --k) // sorted by end
                if (st > intv[k].End) break;
            ++k;
            if (k < intv.Count)
            {
                // overlaps one or more saved intervals: widen the leftmost and drop the rest
                intv[k] = (Math.Min(intv[k].Start, st), Math.Max(intv[k].End, en));
                intv.RemoveRange(k + 1, intv.Count - (k + 1));
            }
            else
            {
                intv.Add((st, en));
            }
        }

        /// <summary>ld_dust1 over 2-bit codes (−1 = non-ACGT).</summary>
        public List<(int Start, int End)> Dust1(int[] seq)
        {
            var intv = new List<(int Start, int End)>();
            int mask = (1 << (2 * _k)) - 1;
            int st = -1, en = -1, lastQ = -1, x = 0, l = 0;
            var ht = new ushort[mask + 1];
            double htSum = 0.0;
            _q.Clear();

            for (int i = 0; i <= seq.Length; ++i)
            {
                int b = i < seq.Length ? seq[i] : -1;
                int ambi;
                if (b >= 0)
                {
                    x = ((x << 2) | b) & mask;
                    ++l;
                    ambi = l < _k ? 1 : 0;
                }
                else
                {
                    l = 0;
                    ambi = 1;
                }
                if (_q.Count >= _ws)
                {
                    int p = _q.Shift();
                    if ((p & 1) == 0) htSum -= _c[ht[p >> 1]--];
                    if (lastQ == 0)
                    {
                        if ((p & 1) == 0 && _ht[p >> 1] > 0) --_ht[p >> 1];
                    }
                    else
                    {
                        --lastQ;
                    }
                }
                _q.Push((x << 1) | ambi);
                if (ambi != 0) continue;
                htSum += _c[++ht[x]];

                int j = -1;
                if (ht[x] >= _minStartCnt)
                {
                    int qn = _q.Count;
                    double swin = htSum - _f[qn] - qn * _thres; // full-window score
                    if (i == en && (lastQ == 0 || i - st >= qn) && swin > 0.0)
                        j = Extend();
                    if (j < 0 && IfBackward(ht))
                        j = Backward(ht);
                }
                if (j >= 0)
                {
                    int st2 = i - (_q.Count - 1 - j) - (_k - 1); // LCR start
                    if (st2 < en)
                    {
                        if (st < 0 || st2 < st) st = st2;
                    }
                    else
                    {
                        if (st >= 0) SaveInterval(intv, st, en);
                        st = st2;
                    }
                    en = i + 1;
                    lastQ = j;
                }
            }
            if (st >= 0) SaveInterval(intv, st, en);
            return intv;
        }
    }

    #endregion

    #endregion

    #region Lempel-Ziv Complexity (compression-based)

    // Lempel-Ziv (1976) complexity c(S): the number of components of the exhaustive
    // history of S. A component (factor) starting at position p is extended while it is
    // still reproducible from the text preceding its last symbol (i.e. it occurs starting
    // at some position < p, overlap allowed); the first non-reproducible extension closes
    // the component (a trailing reproducible remainder is also a component).
    // Ref: Lempel A, Ziv J (1976) "On the Complexity of Finite Sequences",
    //      IEEE Trans. Inf. Theory 22(1):75-81, doi:10.1109/TIT.1976.1055501
    //      (example 0001101001000101 = 0.001.10.100.1000.101, c = 6).
    // Reference scan: Kaspar F, Schuster HG (1987) Phys. Rev. A 36(2):842-848, as implemented
    //      in antropy.entropy._lz_complexity; computed here via the Longest-Previous-Factor
    //      array (Crochemore & Ilie 2008) in O(n log² n), value-identical to the scan.
    // NOTE: this is NOT the LZ78 incremental ("set of seen phrases") parse, which gives
    //      8 for 1001111011000010 instead of the LZ76 value 6.

    /// <summary>
    /// Calculates the raw Lempel–Ziv (1976) complexity c(S) of a DNA sequence: the number
    /// of components in the exhaustive history of S (Lempel &amp; Ziv 1976); values are identical
    /// to the Kaspar–Schuster (1987) scan (antropy), computed via the longest-previous-factor
    /// array in O(n log² n). Higher values indicate more complex (less compressible)
    /// sequences; a homopolymer of length ≥ 2 has c = 2.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <returns>Number of Lempel–Ziv (LZ76) components (≥ 0).</returns>
    public static int CalculateLempelZivComplexity(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateLempelZivComplexityCore(sequence.Sequence.ToUpperInvariant());
    }

    /// <summary>
    /// Calculates the raw Lempel–Ziv (1976) complexity from a raw sequence string
    /// (upper-cased; any symbol alphabet). Null/empty → 0.
    /// </summary>
    public static int CalculateLempelZivComplexity(string sequence)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateLempelZivComplexityCore(sequence.ToUpperInvariant());
    }

    /// <summary>
    /// Calculates the normalized Lempel–Ziv complexity (Zhang et al. 2009):
    /// c / (n / log_b(n)), where c is the raw LZ76 complexity, n the sequence length and
    /// b the alphabet size (number of distinct symbols present). Random sequences give
    /// values near 1; repetitive sequences give smaller values.
    /// Following the reference implementation (antropy <c>lziv_complexity</c>), when fewer
    /// than two distinct symbols are present the base is clamped to 2 (b := max(b, 2)).
    /// For a length-1 input log_b(1) = 0 makes the formula undefined (the reference raises
    /// a division-by-zero); this implementation returns the raw count (1) instead.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <returns>Normalized Lempel–Ziv complexity.</returns>
    public static double CalculateNormalizedLempelZivComplexity(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateNormalizedLempelZivComplexityCore(sequence.Sequence.ToUpperInvariant());
    }

    /// <summary>
    /// Calculates the normalized Lempel–Ziv complexity from a raw sequence string.
    /// </summary>
    public static double CalculateNormalizedLempelZivComplexity(string sequence)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateNormalizedLempelZivComplexityCore(sequence.ToUpperInvariant());
    }

    /// <summary>
    /// Estimates sequence complexity using a compression-based measure.
    /// Returns the normalized Lempel–Ziv (1976) complexity c / (n / log_b(n)) (Zhang et al.
    /// 2009); lower values indicate more repetitive/compressible sequences, ≈ 1 for random.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <returns>Normalized Lempel–Ziv complexity.</returns>
    public static double EstimateCompressionRatio(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateNormalizedLempelZivComplexity(sequence);
    }

    /// <summary>
    /// Estimates compression-based complexity (normalized Lempel–Ziv) from a raw
    /// sequence string.
    /// </summary>
    public static double EstimateCompressionRatio(string sequence)
    {
        return CalculateNormalizedLempelZivComplexity(sequence);
    }

    private static int CalculateLempelZivComplexityCore(string s)
    {
        // LZ76 exhaustive history from the Longest-Previous-Factor array:
        // LPF[q] = max_{j<q} lcp(S[q..], S[j..]) (overlap allowed). The component starting
        // at q is S[q .. q+LPF[q]] (the longest copyable prefix plus one new symbol); if the
        // copy reaches the end of S, the reproducible remainder is the last component.
        // Equivalent to the Kaspar–Schuster (1987) scan (antropy _lz_complexity) but
        // O(n log n) instead of O(n²/log n) on random input.
        int n = s.Length;
        if (n == 0) return 0;

        int[] lpf = ComputeLongestPreviousFactor(s);

        int complexity = 0;
        int q = 0;
        while (q < n)
        {
            complexity++;
            q += lpf[q] + 1; // q + lpf[q] >= n ⇒ remainder reproducible ⇒ last component
        }

        return complexity;
    }

    /// <summary>
    /// Longest Previous Factor array (Crochemore &amp; Ilie 2008): for every position q, the
    /// length of the longest prefix of S[q..] that also starts at some j &lt; q.
    /// Built from the suffix array (prefix doubling) and the Kasai et al. (2001) LCP array;
    /// the best earlier suffix is the nearest suffix with a smaller text position on either
    /// side in suffix-array order.
    /// </summary>
    private static int[] ComputeLongestPreviousFactor(string s)
    {
        int n = s.Length;
        var symbols = new int[n];
        for (int i = 0; i < n; i++) symbols[i] = s[i];
        int[] sa = BuildSuffixArray(symbols);
        int[] lcp = BuildLcpArray(symbols, sa);

        var lpf = new int[n];
        var stackPos = new int[n];
        var stackMin = new int[n];

        // Forward pass: nearest previous rank with a smaller text position.
        int top = -1;
        for (int r = 0; r < n; r++)
        {
            int running = r > 0 ? lcp[r] : 0;
            while (top >= 0 && stackPos[top] > sa[r])
            {
                top--;
                if (top >= 0) running = Math.Min(running, stackMin[top]);
            }
            if (top >= 0)
            {
                lpf[sa[r]] = running;
                stackMin[top] = running; // min LCP over ranks (top, r]
            }
            stackPos[++top] = sa[r];
            stackMin[top] = int.MaxValue;
        }

        // Backward pass: nearest following rank with a smaller text position.
        top = -1;
        for (int r = n - 1; r >= 0; r--)
        {
            int running = r + 1 < n ? lcp[r + 1] : 0;
            while (top >= 0 && stackPos[top] > sa[r])
            {
                top--;
                if (top >= 0) running = Math.Min(running, stackMin[top]);
            }
            if (top >= 0)
            {
                if (running > lpf[sa[r]]) lpf[sa[r]] = running;
                stackMin[top] = running;
            }
            stackPos[++top] = sa[r];
            stackMin[top] = int.MaxValue;
        }

        return lpf;
    }

    /// <summary>
    /// Suffix array by prefix doubling (Manber &amp; Myers 1993), O(n log² n) with a comparison
    /// sort on (rank[i], rank[i+k]) keys. Symbols are arbitrary integers (compared by value);
    /// shared by the LZ76 factorization here and by the maximal-repeat enumeration in
    /// <see cref="RepeatFinder"/> (REP-DIRECT-001), which encodes non-ACGT symbols as unique values.
    /// </summary>
    internal static int[] BuildSuffixArray(int[] s)
    {
        int n = s.Length;
        var sa = new int[n];
        if (n == 0) return sa;
        var rank = new int[n];
        var tmp = new int[n];
        var keys = new long[n];

        // Initial ranks: dense 0..d-1 by symbol value (d ≤ n), so every key fits in rank·(n+2)+second.
        for (int i = 0; i < n; i++) sa[i] = i;
        var initial = (int[])s.Clone();
        Array.Sort(initial, sa);
        rank[sa[0]] = 0;
        for (int r = 1; r < n; r++)
            rank[sa[r]] = rank[sa[r - 1]] + (initial[r] != initial[r - 1] ? 1 : 0);

        long radix = n + 2L;
        int k = 1;
        while (rank[sa[n - 1]] != n - 1 && k < n)
        {
            for (int r = 0; r < n; r++)
            {
                int i = sa[r];
                long second = i + k < n ? rank[i + k] + 1L : 0L;
                keys[r] = rank[i] * radix + second;
            }
            Array.Sort(keys, sa);

            tmp[sa[0]] = 0;
            for (int r = 1; r < n; r++)
                tmp[sa[r]] = tmp[sa[r - 1]] + (keys[r] != keys[r - 1] ? 1 : 0);
            Array.Copy(tmp, rank, n);
            k <<= 1;
        }

        return sa;
    }

    /// <summary>
    /// LCP array (Kasai et al. 2001): <c>lcp[r]</c> = length of the longest common prefix of the
    /// suffixes at <c>sa[r-1]</c> and <c>sa[r]</c>; <c>lcp[0] = 0</c>. O(n).
    /// </summary>
    internal static int[] BuildLcpArray(int[] s, int[] sa)
    {
        int n = s.Length;
        var rank = new int[n];
        for (int r = 0; r < n; r++) rank[sa[r]] = r;

        var lcp = new int[n];
        int h = 0;
        for (int i = 0; i < n; i++)
        {
            if (rank[i] > 0)
            {
                int j = sa[rank[i] - 1];
                while (i + h < n && j + h < n && s[i + h] == s[j + h]) h++;
                lcp[rank[i]] = h;
                if (h > 0) h--;
            }
            else
            {
                h = 0;
            }
        }

        return lcp;
    }

    private static double CalculateNormalizedLempelZivComplexityCore(string seq)
    {
        int n = seq.Length;
        if (n == 0) return 0;

        int c = CalculateLempelZivComplexityCore(seq);

        // Alphabet size b = number of distinct symbols actually present.
        var alphabet = new HashSet<char>();
        foreach (char ch in seq) alphabet.Add(ch);
        int b = alphabet.Count;

        // antropy reference: `base = 2 if base < 2 else base`.
        if (b < MinAlphabetForNormalization) b = MinAlphabetForNormalization;

        // b(n) = n / log_b(n); normalized complexity = c / b(n).
        double logBaseN = Math.Log(n) / Math.Log(b);
        if (logBaseN <= 0) return c; // n == 1 ⇒ log_b(1) = 0: undefined, return raw count (1)

        double upperBound = n / logBaseN;
        return c / upperBound;
    }

    // Reference (antropy lziv_complexity) clamps the log base to 2 when fewer than 2
    // distinct symbols are present, so log_b(n) stays defined.
    // Ref: Zhang et al. (2009) normalized LZ; antropy lziv_complexity.
    private const int MinAlphabetForNormalization = 2;

    #endregion
}

/// <summary>
/// A point in complexity analysis.
/// </summary>
public readonly record struct ComplexityPoint(
    int Position,
    double ShannonEntropy,
    double LinguisticComplexity,
    int WindowStart,
    int WindowEnd);

/// <summary>
/// A low-complexity region detected in the sequence.
/// </summary>
public readonly record struct LowComplexityRegion(
    int Start,
    int End,
    int Length,
    double MinEntropy,
    string Sequence);
