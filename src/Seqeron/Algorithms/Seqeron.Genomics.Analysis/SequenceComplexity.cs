namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Calculates various sequence complexity metrics for detecting low-complexity regions,
/// repetitive sequences, and information content.
/// </summary>
public static class SequenceComplexity
{
    #region Linguistic Complexity

    /// <summary>
    /// Calculates linguistic complexity (LC) as the ratio of observed to possible subwords.
    /// LC = 1.0 for maximum complexity, lower values indicate repeats/low complexity.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="maxWordLength">Maximum word length to consider.</param>
    /// <returns>Linguistic complexity (0 to 1).</returns>
    public static double CalculateLinguisticComplexity(DnaSequence sequence, int maxWordLength = 10)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxWordLength, 1);

        return CalculateLinguisticComplexityCore(sequence.Sequence, maxWordLength);
    }

    /// <summary>
    /// Calculates linguistic complexity from a raw sequence string.
    /// </summary>
    public static double CalculateLinguisticComplexity(string sequence, int maxWordLength = 10)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateLinguisticComplexityCore(sequence.ToUpperInvariant(), maxWordLength);
    }

    private static double CalculateLinguisticComplexityCore(string seq, int maxWordLength)
    {
        if (seq.Length == 0) return 0;

        long observedTotal = 0;
        long possibleTotal = 0;

        for (int wordLen = 1; wordLen <= Math.Min(maxWordLength, seq.Length); wordLen++)
        {
            var observedWords = new HashSet<string>();

            for (int i = 0; i <= seq.Length - wordLen; i++)
            {
                observedWords.Add(seq.Substring(i, wordLen));
            }

            observedTotal += observedWords.Count;

            // Maximum possible words of length wordLen
            long maxPossible = Math.Min(
                (long)Math.Pow(4, wordLen),       // 4^wordLen possible DNA words
                seq.Length - wordLen + 1);         // Can't observe more than available positions

            possibleTotal += maxPossible;
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
    /// Calculates Shannon entropy from a raw sequence string.
    /// </summary>
    public static double CalculateShannonEntropy(string sequence)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateShannonEntropyCore(sequence.ToUpperInvariant());
    }

    private static double CalculateShannonEntropyCore(string seq)
    {
        if (seq.Length == 0) return 0;

        var frequencies = new Dictionary<char, int> { ['A'] = 0, ['T'] = 0, ['G'] = 0, ['C'] = 0 };

        foreach (char c in seq)
        {
            if (frequencies.TryGetValue(c, out int value))
                frequencies[c] = ++value;
        }

        double entropy = 0;
        int total = frequencies.Values.Sum();

        if (total == 0) return 0;

        foreach (int count in frequencies.Values)
        {
            if (count > 0)
            {
                double p = (double)count / total;
                entropy -= p * Math.Log2(p);
            }
        }

        return entropy;
    }

    /// <summary>
    /// Calculates the Shannon entropy (in bits) of the overlapping k-mer frequency
    /// distribution of the sequence.
    /// </summary>
    /// <remarks>
    /// The sequence is decomposed into its L-k+1 overlapping k-mers (sliding window,
    /// one base step). With n_i the count of distinct k-mer i and N = L-k+1 the total
    /// number of k-mers, p_i = n_i / N and the entropy is H = -Σ p_i · log₂(p_i)
    /// (Shannon 1948). Entropy is reported in bits (log base 2): it is 0 when a single
    /// k-mer dominates (deterministic distribution) and reaches log₂(N) when every k-mer
    /// is distinct (uniform distribution). See longdust (Li 2025) for the k-mer-frequency
    /// formulation used to detect low-complexity DNA.
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
    {
        if (seq.Length < k) return 0;

        var kmerCounts = new Dictionary<string, int>();
        int total = 0;

        for (int i = 0; i <= seq.Length - k; i++)
        {
            string kmer = seq.Substring(i, k);
            if (kmerCounts.TryGetValue(kmer, out int value))
                kmerCounts[kmer] = ++value;
            else
                kmerCounts[kmer] = 1;
            total++;
        }

        double entropy = 0;
        foreach (int count in kmerCounts.Values)
        {
            double p = (double)count / total;
            entropy -= p * Math.Log2(p);
        }

        return entropy;
    }

    #endregion

    #region Sliding Window Complexity

    // Per-window linguistic-complexity vocabulary cap. Following Gabrielian & Bolshoy (1999),
    // the linguistic-complexity assessment limits vocabulary evaluation to a bounded set of
    // word lengths (W) rather than all N-1 lengths, for computational efficiency. Sequence
    // complexity and DNA curvature, Comput. Chem. 23(3-4):263-274. doi:10.1016/S0097-8485(99)00007-8
    private const int WindowLcMaxWordLength = 6;

    /// <summary>
    /// Calculates complexity across the sequence using a sliding window (a complexity
    /// profile, in the sense of Troyanskaya et al. (2002)). For each window fully contained
    /// in the sequence the per-window Shannon entropy (bits, Shannon 1948) and linguistic
    /// complexity (summation form) are reported with the window's coordinates.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Size of the sliding window (default: 64).</param>
    /// <param name="stepSize">Step size for window movement (default: 10).</param>
    /// <returns>Complexity values with positions.</returns>
    public static IEnumerable<ComplexityPoint> CalculateWindowedComplexity(
        DnaSequence sequence,
        int windowSize = 64,
        int stepSize = 10)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(stepSize, 1);

        return CalculateWindowedComplexityCore(sequence.Sequence, windowSize, stepSize);
    }

    private static IEnumerable<ComplexityPoint> CalculateWindowedComplexityCore(
        string seq,
        int windowSize,
        int stepSize)
    {
        for (int i = 0; i + windowSize <= seq.Length; i += stepSize)
        {
            string window = seq.Substring(i, windowSize);
            double entropy = CalculateShannonEntropyCore(window);
            double lc = CalculateLinguisticComplexityCore(window, Math.Min(WindowLcMaxWordLength, windowSize));

            yield return new ComplexityPoint(
                Position: i + windowSize / 2,
                ShannonEntropy: entropy,
                LinguisticComplexity: lc,
                WindowStart: i,
                WindowEnd: i + windowSize - 1);
        }
    }

    #endregion

    #region Low Complexity Regions

    /// <summary>
    /// Finds low-complexity regions in the sequence.
    /// Uses a combination of entropy and linguistic complexity.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Window size for analysis (default: 64).</param>
    /// <param name="entropyThreshold">Entropy threshold below which regions are considered low complexity (default: 1.0).</param>
    /// <returns>Low-complexity regions.</returns>
    public static IEnumerable<LowComplexityRegion> FindLowComplexityRegions(
        DnaSequence sequence,
        int windowSize = 64,
        double entropyThreshold = 1.0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, 1);

        return FindLowComplexityRegionsCore(sequence.Sequence, windowSize, entropyThreshold);
    }

    private static IEnumerable<LowComplexityRegion> FindLowComplexityRegionsCore(
        string seq,
        int windowSize,
        double entropyThreshold)
    {
        if (seq.Length < windowSize) yield break;

        int? regionStart = null;
        double minEntropy = double.MaxValue;

        for (int i = 0; i + windowSize <= seq.Length; i++)
        {
            string window = seq.Substring(i, windowSize);
            double entropy = CalculateShannonEntropyCore(window);

            if (entropy < entropyThreshold)
            {
                if (regionStart == null)
                {
                    regionStart = i;
                    minEntropy = entropy;
                }
                else
                {
                    minEntropy = Math.Min(minEntropy, entropy);
                }
            }
            else if (regionStart != null)
            {
                // End of low-complexity region
                int end = i + windowSize - 1;
                yield return new LowComplexityRegion(
                    Start: regionStart.Value,
                    End: end,
                    Length: end - regionStart.Value + 1,
                    MinEntropy: minEntropy,
                    Sequence: seq.Substring(regionStart.Value, end - regionStart.Value + 1));

                regionStart = null;
                minEntropy = double.MaxValue;
            }
        }

        // Handle region at end of sequence
        if (regionStart != null)
        {
            int end = seq.Length - 1;
            yield return new LowComplexityRegion(
                Start: regionStart.Value,
                End: end,
                Length: end - regionStart.Value + 1,
                MinEntropy: minEntropy,
                Sequence: seq.Substring(regionStart.Value));
        }
    }

    #endregion

    #region Dust Score

    // DUST/SDUST uses overlapping nucleotide triplets (3-mers) as the default word.
    // k = 3 is hardcoded in the reference implementation; Morgulis et al. (2006),
    // J Comput Biol 13(5):1028-1040, doi:10.1089/cmb.2006.13.1028, and lh3/sdust.
    private const int DustWordSize = 3;

    // Mask/low-complexity threshold for the DUST score: 2.0, corresponding to the
    // reference default level T = 20 (lh3/sdust: "rw*10 > L*T" ⇔ score > T/10 = 2.0).
    private const double DustMaskThreshold = 2.0;

    /// <summary>
    /// Calculates the DUST low-complexity score of a sequence (Morgulis et al. 2006).
    /// The score is Σ_t c_t·(c_t−1)/2 over all overlapping words t, divided by the number
    /// of words (L − wordSize + 1, equal to L − 2 for triplets). A HIGHER score indicates
    /// LOWER complexity (more repeated words); fully distinct words give 0.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="wordSize">Word size (default: 3, as defined by DUST/SDUST).</param>
    /// <returns>DUST score (≥ 0); 0 when the sequence is shorter than one word.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="wordSize"/> &lt; 1.</exception>
    public static double CalculateDustScore(DnaSequence sequence, int wordSize = DustWordSize)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(wordSize, 1);
        return CalculateDustScoreCore(sequence.Sequence, wordSize);
    }

    /// <summary>
    /// Calculates the DUST low-complexity score from a raw sequence string. The string is
    /// upper-cased to match the normalization applied by <see cref="DnaSequence"/>.
    /// </summary>
    /// <param name="sequence">Raw sequence string; null or empty yields 0.</param>
    /// <param name="wordSize">Word size (default: 3, as defined by DUST/SDUST).</param>
    /// <returns>DUST score (≥ 0); 0 when null/empty or shorter than one word.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="wordSize"/> &lt; 1.</exception>
    public static double CalculateDustScore(string sequence, int wordSize = DustWordSize)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(wordSize, 1);
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateDustScoreCore(sequence.ToUpperInvariant(), wordSize);
    }

    private static double CalculateDustScoreCore(string seq, int wordSize)
    {
        if (seq.Length < wordSize) return 0;

        var wordCounts = new Dictionary<string, int>();

        // Number of overlapping words; equals L − 2 for the default triplet word.
        int wordCount = seq.Length - wordSize + 1;

        for (int i = 0; i < wordCount; i++)
        {
            string word = seq.Substring(i, wordSize);
            if (wordCounts.TryGetValue(word, out int value))
                wordCounts[word] = ++value;
            else
                wordCounts[word] = 1;
        }

        // DUST score numerator: Σ_t c_t·(c_t−1)/2 over all distinct words
        // (Morgulis et al. 2006; longdust restatement S = Σ c(c−1)/2 / (L−2), Li 2025).
        double sum = 0;
        foreach (int count in wordCounts.Values)
        {
            // Promote to double before multiplying: for a highly repetitive sequence a
            // single word's count can approach L, and count·(count−1) would overflow a
            // 32-bit int (e.g. L ≈ 2·10⁵ ⇒ count·(count−1) ≈ 4·10¹⁰ > int.MaxValue),
            // silently corrupting the Σ c(c−1)/2 numerator (Morgulis 2006; Li 2025).
            sum += (double)count * (count - 1) / 2.0;
        }

        // Normalize by the number of words (L − wordSize + 1 = L − 2 for triplets), per
        // the 1/(L−2) factor in Li (2025) and lh3/sdust's per-triplet running length.
        return sum / wordCount;
    }

    /// <summary>
    /// Masks low-complexity regions using DUST algorithm.
    /// </summary>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="windowSize">Window size for masking (default: 64).</param>
    /// <param name="threshold">DUST threshold above which to mask (default: 2.0).</param>
    /// <param name="maskChar">Character to use for masking (default: 'N').</param>
    /// <returns>Masked sequence.</returns>
    public static string MaskLowComplexity(
        DnaSequence sequence,
        int windowSize = 64,
        double threshold = DustMaskThreshold,
        char maskChar = 'N')
    {
        ArgumentNullException.ThrowIfNull(sequence);

        return MaskLowComplexityCore(sequence.Sequence, windowSize, threshold, maskChar);
    }

    private static string MaskLowComplexityCore(string seq, int windowSize, double threshold, char maskChar)
    {
        if (seq.Length < windowSize) return seq;

        var masked = new char[seq.Length];
        seq.CopyTo(0, masked, 0, seq.Length);

        for (int i = 0; i + windowSize <= seq.Length; i++)
        {
            string window = seq.Substring(i, windowSize);
            double dustScore = CalculateDustScoreCore(window, DustWordSize);

            if (dustScore > threshold)
            {
                for (int j = i; j < i + windowSize; j++)
                {
                    masked[j] = maskChar;
                }
            }
        }

        return new string(masked);
    }

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
        int[] sa = BuildSuffixArray(s);

        var rank = new int[n];
        for (int r = 0; r < n; r++) rank[sa[r]] = r;

        // Kasai LCP: lcp[r] = lcp(S[sa[r-1]..], S[sa[r]..]), lcp[0] = 0.
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
    /// sort on (rank[i], rank[i+k]) keys.
    /// </summary>
    private static int[] BuildSuffixArray(string s)
    {
        int n = s.Length;
        var sa = new int[n];
        var rank = new int[n];
        var tmp = new int[n];
        var keys = new long[n];

        for (int i = 0; i < n; i++)
        {
            sa[i] = i;
            rank[i] = s[i];
        }

        int k = 1;
        while (true)
        {
            for (int r = 0; r < n; r++)
            {
                int i = sa[r];
                long second = i + k < n ? rank[i + k] + 1L : 0L;
                keys[r] = rank[i] * (long)(Math.Max(n, char.MaxValue) + 2) + second;
            }
            Array.Sort(keys, sa);

            tmp[sa[0]] = 0;
            for (int r = 1; r < n; r++)
                tmp[sa[r]] = tmp[sa[r - 1]] + (keys[r] != keys[r - 1] ? 1 : 0);
            Array.Copy(tmp, rank, n);

            if (rank[sa[n - 1]] == n - 1 || k >= n) break;
            k <<= 1;
        }

        return sa;
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
