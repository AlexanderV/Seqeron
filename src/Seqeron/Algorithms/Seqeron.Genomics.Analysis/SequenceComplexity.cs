namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Calculates various sequence complexity metrics for detecting low-complexity regions,
/// repetitive sequences, and information content.
/// </summary>
public static class SequenceComplexity
{
    #region Linguistic Complexity

    /// <summary>
    /// Calculates linguistic complexity (LC), the summation form
    /// LC = Σ_{i=1..m} V_i / Σ_{i=1..m} V_max,i with V_max,i = min(4^i, N − i + 1),
    /// where V_i is the number of distinct subwords of length i (Orlov &amp; Potapov 2004, NAR 32:W628,
    /// word length limited by m ≤ N). With m ≥ N this is exactly the Troyanskaya et al. (2002,
    /// Bioinformatics 18:679) definition LC = A(s)/M(s) over all lengths 1..N (Rosalind LING).
    /// This is not Trifonov's (1990) product form C = Π U_i.
    /// LC = 1.0 for maximum complexity, lower values indicate repeats/low complexity.
    /// </summary>
    /// <remarks>
    /// Small m uses direct hash enumeration; larger m counts V_i from the sequence's suffix tree
    /// (Troyanskaya et al. 2002): V_i equals the number of suffix-tree edges spanning depth i, so
    /// the full-length LC is computed in linear time. Both paths return identical values.
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="maxWordLength">Maximum word length m to consider (values ≥ N give Troyanskaya's all-length LC).</param>
    /// <returns>Linguistic complexity (0 to 1).</returns>
    public static double CalculateLinguisticComplexity(DnaSequence sequence, int maxWordLength = 10)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxWordLength, 1);

        string seq = sequence.Sequence;
        if (seq.Length == 0) return 0;
        int m = Math.Min(maxWordLength, seq.Length);
        return m > LcHashEnumerationMaxWordLength
            ? LinguisticComplexityFromCounts(DistinctSubwordCountsFromSuffixTree(sequence.SuffixTree, seq.Length, m), seq.Length)
            : CalculateLinguisticComplexityCore(seq, maxWordLength);
    }

    /// <summary>
    /// Calculates linguistic complexity from a raw sequence string (same definition as the
    /// <see cref="CalculateLinguisticComplexity(DnaSequence, int)"/> overload; input is upper-cased;
    /// null/empty input or <paramref name="maxWordLength"/> &lt; 1 returns 0).
    /// </summary>
    public static double CalculateLinguisticComplexity(string sequence, int maxWordLength = 10)
    {
        if (string.IsNullOrEmpty(sequence)) return 0;
        return CalculateLinguisticComplexityCore(sequence.ToUpperInvariant(), maxWordLength);
    }

    /// <summary>Largest word length for which direct hash enumeration is used instead of the suffix tree.</summary>
    private const int LcHashEnumerationMaxWordLength = 12;

    private static double CalculateLinguisticComplexityCore(string seq, int maxWordLength)
    {
        if (seq.Length == 0) return 0;

        int m = Math.Min(maxWordLength, seq.Length);
        if (m < 1) return 0;

        if (m > LcHashEnumerationMaxWordLength)
            return LinguisticComplexityFromCounts(
                DistinctSubwordCountsFromSuffixTree(global::SuffixTree.SuffixTree.Build(seq), seq.Length, m), seq.Length);

        var counts = new long[m + 1];
        for (int wordLen = 1; wordLen <= m; wordLen++)
        {
            var observedWords = new HashSet<string>();

            for (int i = 0; i <= seq.Length - wordLen; i++)
            {
                observedWords.Add(seq.Substring(i, wordLen));
            }

            counts[wordLen] = observedWords.Count;
        }

        return LinguisticComplexityFromCounts(counts, seq.Length);
    }

    /// <summary>
    /// LC = Σ V_i / Σ min(4^i, N − i + 1) for i = 1..counts.Length−1.
    /// </summary>
    private static double LinguisticComplexityFromCounts(long[] counts, int n)
    {
        long observedTotal = 0;
        long possibleTotal = 0;

        for (int wordLen = 1; wordLen < counts.Length; wordLen++)
        {
            observedTotal += counts[wordLen];

            // V_max,i = min(4^i, N − i + 1); 4^i > N for i ≥ 16 (N ≤ int.MaxValue), so avoid overflow.
            long positions = n - wordLen + 1;
            long maxPossible = wordLen < 16 ? Math.Min(1L << (2 * wordLen), positions) : positions;
            possibleTotal += maxPossible;
        }

        return possibleTotal > 0 ? (double)observedTotal / possibleTotal : 0;
    }

    /// <summary>
    /// Distinct-subword counts V_1..V_m from a suffix tree (Troyanskaya et al. 2002): each substring
    /// of length i is a unique point at depth i on exactly one edge, so V_i = number of edges whose
    /// depth range covers i. Leaf edges end with the terminator, which is excluded.
    /// </summary>
    private static long[] DistinctSubwordCountsFromSuffixTree(global::SuffixTree.ISuffixTree tree, int n, int m)
    {
        var visitor = new SubwordDepthVisitor(n, m);
        tree.Traverse(visitor);

        var counts = new long[m + 1];
        long running = 0;
        for (int i = 1; i <= m; i++)
        {
            running += visitor.Diff[i];
            counts[i] = running;
        }
        return counts;
    }

    private sealed class SubwordDepthVisitor : global::SuffixTree.ISuffixTreeVisitor
    {
        private readonly int _n;
        private readonly int _m;
        private bool _rootSeen;

        public SubwordDepthVisitor(int n, int m)
        {
            _n = n;
            _m = m;
            Diff = new long[m + 2];
        }

        public long[] Diff { get; }

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (!_rootSeen)
            {
                _rootSeen = true; // root: no edge
                return;
            }

            // Edge covers text[startIndex, endIndex); a leaf edge (endIndex < 0) runs to the text end
            // (the terminator position _n is not a character of the sequence).
            int end = endIndex < 0 || endIndex > _n ? _n : endIndex;
            int edgeLength = end - startIndex;
            if (edgeLength <= 0) return;

            int from = depth + 1;
            if (from > _m) return;
            int to = Math.Min(depth + edgeLength, _m);

            Diff[from]++;
            Diff[to + 1]--;
        }

        public void EnterBranch(int key) { }

        public void ExitBranch() { }
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
    /// Single entropy kernel shared by the per-base and k-mer entropy methods of this class.
    /// </summary>
    private static double ShannonEntropyBits(IEnumerable<int> counts)
    {
        long total = 0;
        foreach (int count in counts)
            total += count;

        if (total == 0) return 0;

        double entropy = 0;
        foreach (int count in counts)
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
    {
        if (seq.Length < k) return 0;

        // Overlapping k-mer tally (N = L − k + 1 windows) via the canonical counter (KMER-COUNT-001);
        // p_i = n_i / N because Σ n_i = N.
        var kmerCounts = KmerAnalyzer.CountKmers(seq, k);
        return ShannonEntropyBits(kmerCounts.Values);
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
    /// <param name="wordSize">Word size (default: 3, as defined by DUST/SDUST). Values other
    /// than 3 are an extrapolation (divisor = number of words − 1); only k = 3 is source-defined.</param>
    /// <returns>DUST score (≥ 0); 0 when fewer than two words exist (ℓ − 1 ≤ 0).</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="wordSize"/> &lt; 1.</exception>
    public static double CalculateDustScore(DnaSequence sequence, int wordSize = DustWordSize)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(wordSize, 1);
        return CalculateDustScoreCore(sequence.Sequence, wordSize);
    }

    /// <summary>
    /// Calculates the DUST low-complexity score from a raw sequence string
    /// (score = Σ_t c_t·(c_t−1)/2 / (ℓ − 1), see <see cref="CalculateDustScore(DnaSequence, int)"/>).
    /// The string is upper-cased to match the normalization applied by <see cref="DnaSequence"/>.
    /// </summary>
    /// <param name="sequence">Raw sequence string; null or empty yields 0.</param>
    /// <param name="wordSize">Word size (default: 3, as defined by DUST/SDUST).</param>
    /// <returns>DUST score (≥ 0); 0 when null/empty or fewer than two words exist.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="wordSize"/> &lt; 1.</exception>
    public static double CalculateDustScore(string sequence, int wordSize = DustWordSize)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(wordSize, 1);
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
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowSize, DustWordSize);
        if (double.IsNaN(threshold) || double.IsInfinity(threshold) || threshold < 0)
            throw new ArgumentOutOfRangeException(nameof(threshold), threshold, "Threshold must be a finite value ≥ 0.");

        return MaskLowComplexityCore(sequence.Sequence, windowSize, threshold, maskChar);
    }

    private static string MaskLowComplexityCore(string seq, int windowSize, double threshold, char maskChar)
    {
        var intervals = FindSdustIntervals(seq, windowSize, threshold);
        if (intervals.Count == 0) return seq;

        var masked = seq.ToCharArray();
        foreach (var (start, end) in intervals)
            Array.Fill(masked, maskChar, start, end - start);

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
    /// </summary>
    private static List<(int Start, int End)> FindSdustIntervals(string seq, int windowSize, double threshold)
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
            int b = i < seq.Length ? NucleotideCode(seq[i]) : 4;
            if (b < 4)
            {
                ++l;
                t = ((t << 2) | b) & DustTripletMask;
                if (l >= DustWordSize)
                {
                    int start = Math.Max(l - windowSize, 0) + (i + 1 - l);
                    SdustSaveMaskedRegions(res, perfect, start);
                    SdustShiftWindow(t, window, threshold, ref suffixLen, ref rw, ref rv, cw, cv);
                    if (rw > threshold * suffixLen)
                        SdustFindPerfect(perfect, window, threshold, start, suffixLen, rv, cv, scratch);
                }
            }
            else
            {
                // Non-ACGT or end of input: flush all pending perfect intervals.
                int start = Math.Max(l - windowSize + 1, 0) + (i + 1 - l);
                while (perfect.Count > 0) SdustSaveMaskedRegions(res, perfect, start++);
                l = t = 0;
            }
        }

        return res;
    }

    private static int NucleotideCode(char c) => c switch
    {
        'A' or 'a' => 0,
        'C' or 'c' => 1,
        'G' or 'g' => 2,
        'T' or 't' => 3,
        _ => 4,
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

    private static void SdustSaveMaskedRegions(List<(int Start, int End)> res, List<SdustPerfectInterval> perfect, int start)
    {
        if (perfect.Count == 0 || perfect[^1].Start >= start) return;

        var p = perfect[^1];
        bool saved = false;
        if (res.Count > 0)
        {
            var (s, f) = res[^1];
            if (p.Start <= f) // overlapping with or adjacent to the previous interval
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
