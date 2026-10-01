using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// RSAT <c>oligo-analysis</c> options beyond the basic significance run (<see cref="MotifFinder.DiscoverMotifs(DnaSequence, int, int, OligoBackgroundModel, OligoStrandMode, bool)"/>):
/// multi-sequence input, z-scores (<c>-return zscore</c>), the pseudo-frequency on expected frequencies (<c>-pseudo</c>),
/// degenerate words with one ambiguous position (<c>-oneN</c> / <c>-onedeg</c>), lexicon backgrounds (<c>-lexicon</c>)
/// and calibration tables (<c>-calibN</c> / <c>-calib1</c>, negative binomial / Poisson P-values).
/// Ported from rsa-tools/rsat-code master 10043f2: <c>perl-scripts/oligo-analysis</c> v1.169 (<c>CountOligos</c>,
/// <c>CalibrateSetFromSingleSequence</c>, <c>CalcSubWordFrequencies</c>, <c>Degenerate</c>, <c>SumReverseComplements</c>,
/// <c>CalcExpected</c>, <c>CalcOverlapCoefficient</c>, <c>CalcZscore</c>, <c>CalcProba</c>), <c>perl-scripts/lib/RSA.lib</c>
/// (<c>ReadCalibration</c>), <c>perl-scripts/lib/RSA.seq.lib</c> (<c>OverlapCoeff</c>), <c>perl-scripts/lib/RSA.disco.lib</c>
/// (<c>NbPossibleOligos</c>, <c>MultiTestCorrections</c>), <c>perl-scripts/lib/RSAT/stats.pm</c> (<c>sum_of_poisson</c>,
/// <c>sum_of_negbin2</c>).
/// </summary>
public static partial class MotifFinder
{
    #region RSAT oligo-analysis options

    /// <summary>
    /// RSAT <c>oligo-analysis</c> occurrence statistics of the length-<paramref name="k"/> words of a sequence set, with
    /// the extended RSAT options in <paramref name="options"/>: per pattern the counts, expected frequency/occurrences,
    /// expected variance and z-score, observed/expected ratio, and occ_P / occ_E / occ_sig.
    /// </summary>
    /// <remarks>
    /// <para>Steps in RSAT order, with S sequences, nb_pos = Σ (Lᵢ − k + 1) over sequences with Lᵢ ≥ k and n = nb_pos
    /// (RSAT <c>sum_occurrences</c>: every window is either counted or an overlap):</para>
    /// <list type="number">
    /// <item>Counting: overlapping windows, or <c>-noov</c> (3′ → 5′ per sequence; a window is discarded when the last
    /// counted occurrence of the same word — or, with both strands, of its reverse complement — starts fewer than k
    /// positions downstream; the discarded window is credited as an overlap to that word, RSAT
    /// <c>$patterns{$rc}->{overlaps}++</c>).</item>
    /// <item><c>-oneN</c> / <c>-onedeg</c>: every word is replaced by the words with one position replaced by each
    /// IUPAC code that contains its residue there (N, or R Y W S M K H B V D N); occ(D) = Σ occ(w) over the words w
    /// matching D (RSAT <c>Degenerate</c>).</item>
    /// <item>Both strands: occ(W|W') = occ(W) + occ(W'), the pair reported under the lexicographically smaller member.</item>
    /// <item>Expected frequency: the background model's exp_freq(w) (for a degenerate word, the sum over its matching
    /// words); with a pseudo-frequency ψ &gt; 0, exp_freq ← (1 − ψ)·exp_freq + ψ / NPO for each strand
    /// (<c>CalcExpected</c>, before the reverse-complement sum), NPO = <see cref="OligoAnalysisReport.PossibleOligos"/>;
    /// exp_freq(W|W') = exp_freq(W) + exp_freq(W') (palindromes once); exp_occ = exp_freq · n.
    /// With a calibration: exp_occ and exp_var come from the table (× S for <see cref="OligoCalibrationMode.PerSequence"/>),
    /// exp_freq = exp_occ / (nb_pos − forbidden positions) — the RSAT display value, doubled for a non-palindromic pair
    /// because RSAT copies the kept member's value to its partner before summing.</item>
    /// <item>Overlap coefficient (Pevzner, Borodovsky &amp; Mironov 1989; RSAT <c>OverlapCoeff</c>):
    /// ovlp(w) = 1 + Σ_{i: w has period i} ∏_{j&lt;i} q(w_j), with q the model's residue probabilities (Bernoulli
    /// models; ¼ for Markov models); a degenerate code counts as the sum of its residues.</item>
    /// <item>Variance and z-score (<c>CalcZscore</c>): calibration → exp_var from the table; <c>-noov</c> →
    /// exp_var = exp_occ; otherwise exp_var = nb_pos · p · (2·ovlp − 1 − (2k + 1)·p) with p = exp_freq;
    /// z = (occ − exp_occ) / √exp_var when exp_var &gt; 0 (NaN otherwise).</item>
    /// <item>P-value over the patterns with occ ≥ <see cref="OligoAnalysisOptions.MinCount"/> and exp_freq &gt; 0:
    /// binomial P(X ≥ occ), X ~ Bin(n, exp_freq); with a calibration P(occ ≤ X ≤ n) under a negative binomial
    /// (exp_occ &lt; exp_var: p = var/mean − 1, size = mean/p) or a Poisson (mean exp_occ); occ_E = occ_P · T with T the
    /// number of patterns tested; occ_sig = −log₁₀ occ_E. Patterns with exp_freq ≤ 0 are reported with NaN statistics
    /// and not tested (RSAT "Cannot calculate probability").</item>
    /// </list>
    /// <para>With default options this returns exactly the values of
    /// <see cref="DiscoverMotifs(DnaSequence, int, int, OligoBackgroundModel, OligoStrandMode, bool)"/> for one sequence.</para>
    /// <para>Deliberate differences from the RSAT code (see the algorithm document): the RSAT 1.169 <c>-oneN</c>/<c>-onedeg</c>
    /// code path produces no output (it reads the undefined global <c>%IUPAC</c> and stores scalars in <c>%patterns</c>),
    /// so the documented behaviour is implemented with exp_freq summed over the matching words (as RSAT does for frequency
    /// files) and n = the number of windows; RSAT negative-binomial terms are rounded to 5 significant digits
    /// (<c>LogToEng</c>), the exact sum is used here; z-scores are not rounded to 2 decimals.</para>
    /// Results are in order of first occurrence (sequence, position), then pattern; calibration-only patterns last.
    /// </remarks>
    /// <param name="sequences">Input DNA sequences (no null elements; sequences shorter than k contribute no windows).</param>
    /// <param name="k">Word length (≥ 1).</param>
    /// <param name="options">Options (null = defaults: Bernoulli from input, single strand, overlapping, no threshold beyond 1).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> &lt; 1, an invalid option value, or a Markov/lexicon order not allowed for k.</exception>
    /// <exception cref="ArgumentException">A null sequence, a calibration word of another length, or degenerate words combined with a calibration.</exception>
    public static OligoAnalysisReport AnalyzeOligos(IEnumerable<DnaSequence> sequences, int k, OligoAnalysisOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        options ??= new OligoAnalysisOptions();
        ValidateOligoOptions(k, options);

        List<string> seqs = SequenceStrings(sequences);

        // RSAT -seqtype prot / other: the residues of the DNA sequences are analysed as letters of that alphabet.
        if (options.SequenceType != OligoSequenceType.Dna)
            return AnalyzeOligosCore(seqs, k, options, OligoResidueAlphabet.For(options.SequenceType, seqs, k));
        return AnalyzeOligosCore(seqs, k, options, null);
    }

    /// <summary>Option checks shared by the <see cref="AnalyzeOligos"/> entry points (thrown before the sequences are read).</summary>
    private static void ValidateOligoOptions(int k, OligoAnalysisOptions options)
    {
        options.Validate(nameof(options));
        var calibration = options.Calibration;
        bool degenerate = options.Degeneracy != OligoDegeneracy.None;
        if (calibration is null)
            options.Background.ValidateFor(k, nameof(options));
        if (calibration is not null && degenerate)
            throw new ArgumentException("Degenerate words cannot be combined with a calibration table.", nameof(options));
        if (calibration is not null && calibration.WordLength != k)
            throw new ArgumentException(
                $"The calibration table holds {calibration.WordLength}-mers, the analysis uses k = {k}.", nameof(options));
        if (options.SequenceType != OligoSequenceType.Dna)
            ValidateResidueAlphabetOptions(options);
    }

    /// <summary>
    /// The RSAT <c>oligo-analysis</c> pipeline over prepared sequences. <paramref name="alphabet"/> null = upper-case ACGT DNA
    /// (no residue filtering, the original code path); otherwise windows with residues outside the alphabet are discarded,
    /// residues outside it are not counted, and NPO / equiprobable / residue probabilities use the alphabet (RSAT <c>-seqtype</c>).
    /// </summary>
    private static OligoAnalysisReport AnalyzeOligosCore(List<string> seqs, int k, OligoAnalysisOptions options, OligoResidueAlphabet? alphabet)
    {
        var background = options.Background;
        var calibration = options.Calibration;
        bool degenerate = options.Degeneracy != OligoDegeneracy.None;
        bool both = options.Strands == OligoStrandMode.Both;
        bool residueAlphabet = alphabet is { IsDna: false };

        // 1. Counting (RSAT CountOligos; windows with residues outside the alphabet are discarded).
        var wordTallies = CountOligoTallies(seqs, k, both, options.CountOverlapping, out long possiblePositions, alphabet);
        var tallies = wordTallies;

        // 2. Degenerate words (RSAT Degenerate).
        char[] codes = options.Degeneracy switch
        {
            OligoDegeneracy.OneN => OneNCodes,
            OligoDegeneracy.OneDegenerate => OneDegenerateCodes,
            _ => Array.Empty<char>(),
        };
        if (degenerate)
            tallies = DegenerateTallies(tallies, k, codes);

        // 3. Calibration entries (RSAT ReadCalibration / CalibrateSetFromSingleSequence).
        Dictionary<string, (double Mean, double Variance)>? calibrated = null;
        if (calibration is not null)
        {
            double factor = calibration.Mode == OligoCalibrationMode.PerSequence ? seqs.Count : 1.0;
            calibrated = calibration.Entries.ToDictionary(
                e => e.Key, e => (e.Value.Mean * factor, e.Value.Variance * factor), StringComparer.Ordinal);
        }

        // 4. Reverse-complement grouping.
        var patterns = GroupTallies(tallies, both, options.MinCount <= 0 ? calibrated?.Keys : null);

        double npo = residueAlphabet ? alphabet!.PossibleOligos(k) : RsatPossibleOligos(k, both, codes.Length);
        double logNpo = residueAlphabet ? alphabet!.LogPossibleOligos(k) : Math.Log(npo);
        long n = possiblePositions;
        double logN = Math.Log(n);
        double psi = options.PseudoFrequency;

        Func<string, double>? logWord = null;
        OligoLexicon? lexicon = null;
        Func<string, double> overlapCoefficient;
        if (residueAlphabet)
        {
            var q = alphabet!.OverlapResidueProbabilities(background, seqs, k);
            overlapCoefficient = w => OverlapCoefficient(w, c => q.GetValueOrDefault(c));
        }
        else
        {
            double[] residues = calibration is null
                ? background.OverlapResidueProbabilities(seqs, k, both)
                : OligoBackgroundModel.InputResidueProbabilities(seqs, k, both);
            overlapCoefficient = w => OverlapCoefficient(w, residues);
        }
        if (calibration is null && patterns.Count > 0)
        {
            if (background.IsLexicon)
            {
                // RSAT CalcSubWordFrequencies: prefixes of occ + overlaps per (plain, single-strand) word.
                lexicon = OligoLexicon.FromWordCounts(
                    wordTallies.Select(t => new KeyValuePair<string, long>(t.Key, t.Value.Occurrences + t.Value.Overlaps)), k);
                var lex = lexicon;
                logWord = w => Math.Log(lex.ExpectedFrequency(w).Frequency);
            }
            else
            {
                logWord = alphabet is null
                    ? background.CreateLogProbability(seqs, k, both)
                    : background.CreateLogProbability(seqs, k, both, alphabet);
            }
        }

        var rows = new List<RowDraft>(patterns.Count);
        foreach (var pattern in patterns)
        {
            var draft = new RowDraft(pattern);
            if (pattern.Occurrences < options.MinCount)
                continue;

            if (calibrated is not null)
            {
                var (mean, variance) = calibrated.TryGetValue(pattern.Word, out var cv) ? cv : (0.0, 0.0);
                long notForbidden = n - pattern.Forbidden;
                double expFreq = notForbidden > 0 ? mean / notForbidden : double.NaN;
                if (psi > 0 && !double.IsNaN(expFreq))
                    expFreq = expFreq * (1 - psi) + psi / npo;
                if (pattern.ReverseComplement != pattern.Word)
                    expFreq += expFreq;
                draft.ExpectedFrequency = expFreq;
                draft.ExpectedOccurrences = mean;
                draft.ExpectedVariance = variance;
                draft.Overlap = overlapCoefficient(pattern.Word);
                // RSAT tests exp_freq > 0; a zero mean (possible only with -pseudo) crashes RSAT's sum_of_poisson.
                draft.Testable = expFreq > 0 && mean > 0;
            }
            else
            {
                double logP = PatternLogFrequency(pattern, logWord!, codes, psi, logNpo);
                draft.LogP = logP;
                draft.ExpectedFrequency = Math.Exp(logP);
                draft.ExpectedOccurrences = Math.Exp(logP + logN);
                draft.Testable = !double.IsNegativeInfinity(logP);
                draft.Overlap = overlapCoefficient(pattern.Word);
                if (!options.CountOverlapping)
                    draft.ExpectedVariance = draft.ExpectedOccurrences;
                else
                {
                    double p = draft.ExpectedFrequency;
                    draft.ExpectedVariance = n * p * (2 * draft.Overlap - 1 - (2 * k + 1) * p);
                }
                if (lexicon is not null && !degenerate)
                    draft.Segmentation = lexicon.ExpectedFrequency(pattern.Word).Segmentation;
            }

            rows.Add(draft);
        }

        int tested = rows.Count(r => r.Testable);
        double logTested = Math.Log(tested);
        var result = new List<OligoStatistics>(rows.Count);
        foreach (var d in rows)
        {
            var pattern = d.Pattern;
            int count = pattern.Occurrences;
            double z = d.ExpectedVariance > 0 ? (count - d.ExpectedOccurrences) / Math.Sqrt(d.ExpectedVariance) : double.NaN;
            double ratio, logOccP;
            var fitted = OligoFittedDistribution.None;
            if (calibrated is not null)
            {
                ratio = d.ExpectedOccurrences == 0 ? 0.0 : count / d.ExpectedOccurrences;
                logOccP = double.NaN;
                if (d.Testable)
                {
                    if (d.ExpectedOccurrences < d.ExpectedVariance)
                    {
                        fitted = OligoFittedDistribution.NegativeBinomial;
                        logOccP = StatisticsHelper.LogNegativeBinomialRangeProbability(
                            count, n, d.ExpectedOccurrences, d.ExpectedVariance);
                    }
                    else
                    {
                        fitted = OligoFittedDistribution.Poisson;
                        logOccP = StatisticsHelper.LogPoissonRangeProbability(count, n, d.ExpectedOccurrences);
                    }
                }
            }
            else
            {
                ratio = d.Testable ? Math.Exp(Math.Log(count) - logN - d.LogP) : 0.0;
                logOccP = double.NaN;
                if (d.Testable)
                {
                    fitted = OligoFittedDistribution.Binomial;
                    logOccP = StatisticsHelper.LogBinomialUpperTail(count, n, d.LogP);
                }
            }

            double logOccE = logOccP + logTested;
            double expectedOccurrences = d.Testable ? d.ExpectedOccurrences : double.NaN; // RSAT: exp_occ "NA"
            result.Add(new OligoStatistics(
                Pattern: pattern.Word,
                ReverseComplement: both ? pattern.ReverseComplement : null,
                Occurrences: count,
                Overlaps: pattern.Overlaps,
                Positions: pattern.Positions.AsReadOnly(),
                ObservedFrequency: (double)count / n,
                ExpectedFrequency: d.ExpectedFrequency,
                ExpectedOccurrences: expectedOccurrences,
                ExpectedVariance: d.ExpectedVariance,
                OverlapCoefficient: d.Overlap,
                ZScore: z,
                Ratio: ratio,
                OccurrenceProbability: Math.Exp(logOccP),
                OccurrenceEValue: Math.Exp(logOccE),
                OccurrenceSignificance: -logOccE / Ln10,
                FittedDistribution: fitted,
                LexiconSegmentation: d.Segmentation));
        }

        return new OligoAnalysisReport(
            Patterns: result.AsReadOnly(),
            OligoLength: k,
            Strands: options.Strands,
            CountOverlapping: options.CountOverlapping,
            Degeneracy: options.Degeneracy,
            SequenceCount: seqs.Count,
            PossiblePositions: possiblePositions,
            TotalOccurrences: n,
            TestedPatterns: tested,
            PossibleOligos: npo)
        {
            SequenceType = alphabet?.Type ?? OligoSequenceType.Dna,
            AlphabetSize = alphabet?.Size ?? 4,
        };
    }

    /// <summary>
    /// <see cref="FindSharedMotifs(IEnumerable{DnaSequence}, int, int, OligoBackgroundModel, OligoStrandMode)"/> with the RSAT
    /// <c>-pseudo</c> correction of the expected frequencies: exp_freq ← (1 − ψ)·exp_freq + ψ / NPO for each strand
    /// before the reverse-complement sum (RSAT <c>CalcExpected</c>); ψ = 0 gives the plain overload exactly.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="pseudoFrequency"/> outside [0, 1] (or the plain overload's conditions).</exception>
    public static SharedMotifAnalysisResult FindSharedMotifs(
        IEnumerable<DnaSequence> sequences,
        int k,
        int minSequences,
        OligoBackgroundModel background,
        OligoStrandMode strands,
        double pseudoFrequency)
    {
        ValidatePseudoFrequency(pseudoFrequency, nameof(pseudoFrequency));
        if (pseudoFrequency == 0.0)
            return FindSharedMotifs(sequences, k, minSequences, background, strands);
        ArgumentNullException.ThrowIfNull(background);
        double logNpo = LogPossibleOligos(k, strands == OligoStrandMode.Both);
        return FindSharedMotifsCore(sequences, k, minSequences, background, strands,
            logWord => w => StatisticsHelper.LogAddExp(Math.Log(1 - pseudoFrequency) + logWord(w), Math.Log(pseudoFrequency) - logNpo));
    }

    /// <summary>RSAT <c>-pseudo</c> / Markov-table pseudo-frequency guard shared by every entry point: ψ ∈ [0, 1].</summary>
    internal static void ValidatePseudoFrequency(double pseudoFrequency, string paramName)
    {
        if (!(pseudoFrequency >= 0.0 && pseudoFrequency <= 1.0))
            throw new ArgumentOutOfRangeException(paramName, pseudoFrequency, "Pseudo-frequency must be in [0, 1].");
    }

    private static readonly char[] OneNCodes = { 'N' };
    private static readonly char[] OneDegenerateCodes = { 'R', 'Y', 'W', 'S', 'M', 'K', 'H', 'B', 'V', 'D', 'N' };

    /// <summary>Per-word counts (RSAT <c>occ</c>, <c>overlaps</c>, <c>forbocc</c>) with the counted positions.</summary>
    private sealed class OligoTally
    {
        public int Occurrences;
        public int Overlaps;
        public long Forbidden;
        public readonly List<OligoOccurrence> Positions = new();
    }

    /// <summary>A pattern (word or pair) after strand grouping.</summary>
    private sealed record GroupedOligo(string Word, string ReverseComplement, int Occurrences, int Overlaps, long Forbidden,
        List<OligoOccurrence> Positions);

    private sealed class RowDraft(GroupedOligo pattern)
    {
        public GroupedOligo Pattern { get; } = pattern;
        public double LogP = double.NegativeInfinity;
        public double ExpectedFrequency;
        public double ExpectedOccurrences;
        public double ExpectedVariance;
        public double Overlap = double.NaN;
        public bool Testable;
        public LexiconSegmentation? Segmentation;
    }

    /// <summary>
    /// RSAT <c>CountOligos</c> over all sequences of length ≥ k. Overlapping mode counts every window; <c>-noov</c> scans each
    /// sequence 3′ → 5′, discards a window overlapping the last counted occurrence of the same word (or of its reverse
    /// complement with both strands) and credits it as an overlap to that word; every counted occurrence adds
    /// min(k − 1, (L − k) − pos) forbidden positions to the word (and to its reverse complement with both strands).
    /// </summary>
    private static Dictionary<string, OligoTally> CountOligoTallies(
        List<string> seqs, int k, bool both, bool overlapping, out long possiblePositions, OligoResidueAlphabet? alphabet = null)
    {
        var tallies = new Dictionary<string, OligoTally>(StringComparer.Ordinal);
        var rcCache = new Dictionary<string, string>(StringComparer.Ordinal);
        possiblePositions = 0;
        OligoTally Tally(string word)
        {
            if (!tallies.TryGetValue(word, out var t))
            {
                t = new OligoTally();
                tallies.Add(word, t);
            }
            return t;
        }

        for (int si = 0; si < seqs.Count; si++)
        {
            string s = seqs[si];
            if (s.Length < k) continue;
            // RSAT: windows with a residue outside the alphabet are counted, then deleted together with their
            // occurrences and overlaps (nb_possible_pos -= discarded); they never interact with valid words, so
            // skipping them is equivalent.
            int[]? invalidPrefix = alphabet?.InvalidPrefixCounts(s);
            bool Valid(int pos) => invalidPrefix is null || invalidPrefix[pos + k] == invalidPrefix[pos];
            possiblePositions += s.Length - k + 1;
            if (invalidPrefix is not null)
            {
                for (int pos = 0; pos + k <= s.Length; pos++)
                {
                    if (!Valid(pos)) possiblePositions--;
                }
            }
            if (overlapping)
            {
                for (int pos = 0; pos + k <= s.Length; pos++)
                {
                    if (!Valid(pos)) continue;
                    var t = Tally(s.Substring(pos, k));
                    t.Occurrences++;
                    t.Positions.Add(new OligoOccurrence(si, pos));
                }
                continue;
            }

            int lastPos = s.Length - k;
            var lastCounted = new Dictionary<string, int>(StringComparer.Ordinal);
            for (int pos = lastPos; pos >= 0; pos--)
            {
                if (!Valid(pos)) continue;
                string word = s.Substring(pos, k);
                if (lastCounted.TryGetValue(word, out int last) && last - pos < k)
                {
                    Tally(word).Overlaps++;
                    continue;
                }
                string rc = both ? CachedReverseComplement(word, rcCache) : word;
                if (both && lastCounted.TryGetValue(rc, out int lastRc) && lastRc - pos < k)
                {
                    Tally(rc).Overlaps++;
                    continue;
                }

                long forbidden = Math.Min(k - 1, lastPos - pos);
                var t = Tally(word);
                t.Forbidden += forbidden;
                if (both) Tally(rc).Forbidden += forbidden;
                t.Occurrences++;
                t.Positions.Add(new OligoOccurrence(si, pos));
                lastCounted[word] = pos;
            }
        }

        foreach (var t in tallies.Values)
            t.Positions.Sort(OligoOccurrence.Compare);
        return tallies;
    }

    /// <summary>
    /// RSAT <c>Degenerate</c>: for every counted word, each position, and each code containing its residue there, the
    /// degenerate word receives the word's occurrences (and overlaps / positions).
    /// </summary>
    private static Dictionary<string, OligoTally> DegenerateTallies(Dictionary<string, OligoTally> tallies, int k, char[] codes)
    {
        var result = new Dictionary<string, OligoTally>(StringComparer.Ordinal);
        var buffer = new char[k];
        foreach (var (word, tally) in tallies)
        {
            if (tally.Occurrences == 0 && tally.Overlaps == 0) continue;
            for (int l = 0; l < k; l++)
            {
                foreach (char code in codes)
                {
                    if (!IupacHelper.MatchesIupac(word[l], code)) continue;
                    word.CopyTo(0, buffer, 0, k);
                    buffer[l] = code;
                    string degenerate = new(buffer);
                    if (!result.TryGetValue(degenerate, out var d))
                    {
                        d = new OligoTally();
                        result.Add(degenerate, d);
                    }
                    d.Occurrences += tally.Occurrences;
                    d.Overlaps += tally.Overlaps;
                    d.Positions.AddRange(tally.Positions);
                }
            }
        }

        foreach (var d in result.Values)
            d.Positions.Sort(OligoOccurrence.Compare);
        return result;
    }

    /// <summary>
    /// Strand grouping (RSAT <c>SumReverseComplements</c> + <c>GroupRC</c>): with both strands occ / overlaps / positions of
    /// W and W' are summed under the lexicographically smaller member, whose own forbidden-position count is kept.
    /// <paramref name="extraWords"/> (calibration-only words, tested with occ 0 when no occurrence threshold is set) are
    /// appended. Order: first occurrence, then word.
    /// </summary>
    private static List<GroupedOligo> GroupTallies(Dictionary<string, OligoTally> tallies, bool both, IEnumerable<string>? extraWords)
    {
        var words = new HashSet<string>(tallies.Keys, StringComparer.Ordinal);
        if (extraWords is not null)
            words.UnionWith(extraWords);

        var grouped = new List<GroupedOligo>();
        var seen = new HashSet<string>(StringComparer.Ordinal);
        foreach (string word in words)
        {
            string rc = both ? DnaSequence.GetReverseComplementString(word) : word;
            string key = string.CompareOrdinal(rc, word) < 0 ? rc : word;
            string partner = key == word ? rc : word;
            if (!seen.Add(key)) continue;

            tallies.TryGetValue(key, out var a);
            OligoTally? b = null;
            if (partner != key) tallies.TryGetValue(partner, out b);
            var positions = new List<OligoOccurrence>();
            if (a is not null) positions.AddRange(a.Positions);
            if (b is not null) positions.AddRange(b.Positions);
            if (a is not null && b is not null) positions.Sort(OligoOccurrence.Compare);
            grouped.Add(new GroupedOligo(key, partner,
                (a?.Occurrences ?? 0) + (b?.Occurrences ?? 0),
                (a?.Overlaps ?? 0) + (b?.Overlaps ?? 0),
                a?.Forbidden ?? 0,
                positions));
        }

        grouped.Sort((x, y) =>
        {
            bool xe = x.Positions.Count == 0, ye = y.Positions.Count == 0;
            if (xe != ye) return xe ? 1 : -1;
            if (!xe)
            {
                int c = OligoOccurrence.Compare(x.Positions[0], y.Positions[0]);
                if (c != 0) return c;
            }
            return string.CompareOrdinal(x.Word, y.Word);
        });
        return grouped;
    }

    /// <summary>
    /// ln exp_freq of a grouped pattern: per strand the word's background frequency (a degenerate word: the sum over its
    /// matching words), optionally corrected by the pseudo-frequency, then summed over W and W'.
    /// </summary>
    private static double PatternLogFrequency(GroupedOligo pattern, Func<string, double> logWord, char[] codes, double psi, double logNpo)
    {
        double Strand(string word)
        {
            double logQ = codes.Length == 0 ? logWord(word) : DegenerateLogFrequency(word, logWord);
            if (psi > 0)
                logQ = StatisticsHelper.LogAddExp(Math.Log(1 - psi) + logQ, Math.Log(psi) - logNpo);
            return logQ;
        }

        double logP = Strand(pattern.Word);
        if (pattern.ReverseComplement != pattern.Word)
            logP = StatisticsHelper.LogAddExp(logP, Strand(pattern.ReverseComplement));
        return Math.Min(logP, 0.0);
    }

    // Σ over the words matching a word with (at most) one IUPAC code (A, C, G, T order at that position).
    private static double DegenerateLogFrequency(string word, Func<string, double> logWord)
    {
        int at = -1;
        for (int i = 0; i < word.Length; i++)
        {
            if (AcgtIndex(word[i]) < 0) { at = i; break; }
        }
        if (at < 0) return logWord(word);

        double sum = double.NegativeInfinity;
        var buffer = word.ToCharArray();
        foreach (char b in AcgtBases)
        {
            if (!IupacHelper.MatchesIupac(b, word[at])) continue;
            buffer[at] = b;
            sum = StatisticsHelper.LogAddExp(sum, logWord(new string(buffer)));
        }
        return sum;
    }

    /// <summary>
    /// RSAT <c>OverlapCoeff</c> (Pevzner, Borodovsky &amp; Mironov 1989): 1 + Σ over the periods i of the word (w_j = w_{j+i})
    /// of ∏_{j&lt;i} q(w_j); q of an IUPAC code is the sum of its residues' probabilities (RSAT <c>CalcExpected</c>).
    /// Reverse-complement overlaps are not added (RSAT's <c>$sum_strands</c> is never set by oligo-analysis).
    /// </summary>
    internal static double OverlapCoefficient(string word, double[] acgt)
        => OverlapCoefficient(word, c => ResidueProbability(c, acgt));

    /// <summary><see cref="OverlapCoefficient(string, double[])"/> with any residue probability function (RSAT <c>%residue_proba</c>).</summary>
    private static double OverlapCoefficient(string word, Func<char, double> q)
    {
        double coeff = 1;
        for (int i = 1; i < word.Length; i++)
        {
            bool period = true;
            for (int j = i; j < word.Length && period; j++)
                period = word[j] == word[j - i];
            if (!period) continue;
            double add = 1;
            for (int j = 0; j < i; j++)
                add *= q(word[j]);
            coeff += add;
        }
        return coeff;
    }

    // Σ q[b] over the bases of an IUPAC code (canonical IupacHelper sets, summed in A, C, G, T order); N is exactly 1
    // (RSAT), any other symbol 0.
    private static double ResidueProbability(char c, double[] q)
    {
        if (c == 'N')
            return 1.0;
        if (!IupacHelper.IsNucleotideCode(c))
            return 0.0;
        double p = 0.0;
        for (int b = 0; b < PwmAlphabetSize; b++)
        {
            if (IupacHelper.MatchesIupac(AcgtBases[b], c))
                p += q[b];
        }
        return p;
    }

    /// <summary>
    /// RSAT <c>NbPossibleOligos</c>: 4^k, or k·|codes|·4^(k−1) with one degenerate position; with both strands
    /// NPO − (NPO − P)/2, P = 4^(k/2) for even k and 0 for odd k (RSAT's palindrome count, also used for degenerate words).
    /// </summary>
    private static double RsatPossibleOligos(int k, bool both, int codes)
    {
        double npo = codes == 0 ? Math.ScaleB(1.0, 2 * k) : k * (double)codes * Math.ScaleB(1.0, 2 * (k - 1));
        if (!both) return npo;
        double palindromes = k % 2 == 0 ? Math.ScaleB(1.0, k) : 0.0;
        return npo - (npo - palindromes) / 2;
    }

    #endregion
}

/// <summary>Degenerate-word mode of RSAT <c>oligo-analysis</c>.</summary>
public enum OligoDegeneracy
{
    /// <summary>Plain words.</summary>
    None,

    /// <summary>One N at any position (RSAT <c>-oneN</c>).</summary>
    OneN,

    /// <summary>One ambiguous IUPAC code (R Y W S M K H B V D N) at any position (RSAT <c>-onedeg</c>).</summary>
    OneDegenerate,
}

/// <summary>How a calibration table scales with the number of input sequences.</summary>
public enum OligoCalibrationMode
{
    /// <summary>Mean and variance per sequence set, used as given (RSAT <c>-calibN</c>).</summary>
    PerSet,

    /// <summary>Single-sequence mean and variance, multiplied by the number of input sequences (RSAT <c>-calib1</c>).</summary>
    PerSequence,
}

/// <summary>Distribution used for occ_P.</summary>
public enum OligoFittedDistribution
{
    /// <summary>Not computed (exp_freq ≤ 0).</summary>
    None,

    /// <summary>Binomial(n, exp_freq) upper tail.</summary>
    Binomial,

    /// <summary>Poisson(exp_occ), calibration with exp_occ ≥ exp_var.</summary>
    Poisson,

    /// <summary>Negative binomial with mean exp_occ and variance exp_var, calibration with exp_occ &lt; exp_var.</summary>
    NegativeBinomial,
}

/// <summary>
/// RSAT oligo calibration table (<c>-calibN</c> / <c>-calib1</c>): the mean and variance of the occurrences of each word in
/// random sequence sets, as produced by RSAT <c>calibrate-oligos</c> → <c>fit-distribution</c>.
/// </summary>
public sealed class OligoCalibration
{
    private OligoCalibration(Dictionary<string, (double Mean, double Variance)> entries, int wordLength, OligoCalibrationMode mode)
    {
        Entries = entries;
        WordLength = wordLength;
        Mode = mode;
    }

    /// <summary>Upper-case word → (mean occurrences, variance), after RSAT's reverse-complement inference.</summary>
    public IReadOnlyDictionary<string, (double Mean, double Variance)> Entries { get; }

    /// <summary>Length of the calibrated words.</summary>
    public int WordLength { get; }

    /// <summary>Per-set or per-sequence calibration.</summary>
    public OligoCalibrationMode Mode { get; }

    /// <summary>
    /// Builds a calibration from word → (mean, variance). As RSAT <c>ReadCalibration</c>, a word whose reverse complement
    /// has no positive mean passes its values to it (words visited in ordinal order).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="entries"/> is null.</exception>
    /// <exception cref="ArgumentException">Empty table, non-ACGT word, mixed lengths, or a negative / non-finite mean, or a non-finite variance.</exception>
    public static OligoCalibration FromEntries(
        IEnumerable<KeyValuePair<string, (double Mean, double Variance)>> entries,
        OligoCalibrationMode mode = OligoCalibrationMode.PerSet)
    {
        ArgumentNullException.ThrowIfNull(entries);
        var table = new Dictionary<string, (double Mean, double Variance)>(StringComparer.Ordinal);
        int length = -1;
        foreach (var (raw, value) in entries)
        {
            string word = ValidateWord(raw, ref length);
            if (!double.IsFinite(value.Mean) || value.Mean < 0)
                throw new ArgumentException($"Invalid expected occurrences {value.Mean} for '{raw}' (must be a real number >= 0).", nameof(entries));
            if (!double.IsFinite(value.Variance))
                throw new ArgumentException($"Invalid variance {value.Variance} for '{raw}'.", nameof(entries));
            table[word] = value;
        }
        if (table.Count == 0)
            throw new ArgumentException("The calibration table is empty.", nameof(entries));

        foreach (string word in table.Keys.Order(StringComparer.Ordinal).ToList())
        {
            string rc = DnaSequence.GetReverseComplementString(word);
            var own = table[word];
            if (own.Mean != 0 && !(table.TryGetValue(rc, out var other) && other.Mean > 0))
                table[rc] = own;
        }

        return new OligoCalibration(table, length, mode);
    }

    /// <summary>
    /// Parses an RSAT calibration file (the format read by <c>oligo-analysis -calibN/-calib1</c>, i.e. RSAT
    /// <c>fit-distribution</c> output): lines starting with <c>;</c>, <c>#</c> or <c>--</c> and blank lines are skipped;
    /// fields are whitespace-separated; field 1 is the word (or a <c>word|rc</c> pair, the first member is used), an
    /// optional <c>word|rc</c> identifier may follow; then the mean (avg), the standard deviation and the variance
    /// (RSAT reads the 1st and 3rd value after the pattern / identifier). Later lines override earlier ones.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    /// <exception cref="FormatException">A data line without a real mean ≥ 0 and a real variance.</exception>
    /// <exception cref="ArgumentException">Non-ACGT word, mixed lengths, or no data line.</exception>
    public static OligoCalibration Parse(string text, OligoCalibrationMode mode = OligoCalibrationMode.PerSet)
    {
        ArgumentNullException.ThrowIfNull(text);
        var entries = new List<KeyValuePair<string, (double Mean, double Variance)>>();
        int lineNumber = 0;
        foreach (string rawLine in text.Split('\n'))
        {
            lineNumber++;
            string line = rawLine.TrimEnd('\r');
            if (line.StartsWith("--", StringComparison.Ordinal) || line.StartsWith(';') || line.StartsWith('#')
                || string.IsNullOrWhiteSpace(line))
                continue;
            var fields = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries).ToList();
            string pattern = fields[0];
            fields.RemoveAt(0);
            int bar = pattern.IndexOf('|');
            if (bar > 0) pattern = pattern[..bar];
            if (fields.Count > 0 && fields[0].Length > 1 && fields[0].AsSpan(1).Contains('|'))
                fields.RemoveAt(0);
            if (fields.Count < 3
                || !double.TryParse(fields[0], NumberStyles.Float, CultureInfo.InvariantCulture, out double mean)
                || !(mean >= 0) || double.IsInfinity(mean)
                || !double.TryParse(fields[2], NumberStyles.Float, CultureInfo.InvariantCulture, out double variance)
                || !double.IsFinite(variance))
                throw new FormatException($"Calibration line {lineNumber}: expected <pattern> [id] <mean> <sd> <variance>: '{line}'.");
            entries.Add(new(pattern, (mean, variance)));
        }
        return FromEntries(entries, mode);
    }

    private static string ValidateWord(string raw, ref int length)
    {
        if (string.IsNullOrEmpty(raw))
            throw new ArgumentException("Empty word in the calibration table.");
        string word = raw.ToUpperInvariant();
        foreach (char c in word)
        {
            if (MotifFinder.AcgtIndex(c) < 0)
                throw new ArgumentException($"Calibration word '{raw}' contains a non-ACGT character.");
        }
        if (length < 0) length = word.Length;
        else if (word.Length != length)
            throw new ArgumentException("All calibration words must have the same length.");
        return word;
    }
}

/// <summary>Options of <see cref="MotifFinder.AnalyzeOligos"/> (RSAT <c>oligo-analysis</c> command-line options).</summary>
public sealed record OligoAnalysisOptions
{
    /// <summary>Occurrence threshold (RSAT <c>-lth occ</c>, default 1); ≤ 0 also tests calibration words without occurrences.</summary>
    public int MinCount { get; init; } = 1;

    /// <summary>Background model (default: Bernoulli from the input, the RSAT default). Ignored for P-values with a calibration.</summary>
    public OligoBackgroundModel Background { get; init; } = OligoBackgroundModel.BernoulliFromInput;

    /// <summary>Single strand (<c>-1str</c>, default) or reverse-complement pairs (<c>-2str</c>).</summary>
    public OligoStrandMode Strands { get; init; } = OligoStrandMode.Single;

    /// <summary>true = <c>-ovlp</c> (default), false = <c>-noov</c>.</summary>
    public bool CountOverlapping { get; init; } = true;

    /// <summary>Pseudo-frequency ψ ∈ [0, 1] on the expected frequencies (RSAT <c>-pseudo</c>, default 0).</summary>
    public double PseudoFrequency { get; init; }

    /// <summary>Degenerate words (RSAT <c>-oneN</c> / <c>-onedeg</c>).</summary>
    public OligoDegeneracy Degeneracy { get; init; } = OligoDegeneracy.None;

    /// <summary>Calibration table (RSAT <c>-calibN</c> / <c>-calib1</c>); null = background model.</summary>
    public OligoCalibration? Calibration { get; init; }

    /// <summary>
    /// Sequence type (RSAT <c>-seqtype dna|prot|other</c>, default <see cref="OligoSequenceType.Dna"/>). Protein and other
    /// sequence types analyse single-strand words only, without degenerate words or calibration tables, with the
    /// equiprobable, input Bernoulli, input Markov or lexicon background (see <see cref="OligoSequenceType"/>).
    /// </summary>
    public OligoSequenceType SequenceType { get; init; } = OligoSequenceType.Dna;

    internal void Validate(string paramName)
    {
        if (Background is null)
            throw new ArgumentException("A background model is required.", paramName);
        MotifFinder.ValidatePseudoFrequency(PseudoFrequency, paramName);
        if (!Enum.IsDefined(Strands))
            throw new ArgumentOutOfRangeException(paramName, Strands, "Unknown strand mode.");
        if (!Enum.IsDefined(Degeneracy))
            throw new ArgumentOutOfRangeException(paramName, Degeneracy, "Unknown degeneracy mode.");
        if (!Enum.IsDefined(SequenceType))
            throw new ArgumentOutOfRangeException(paramName, SequenceType, "Unknown sequence type.");
    }
}

/// <summary>A counted occurrence: sequence index and 0-based window start.</summary>
public readonly record struct OligoOccurrence(int SequenceIndex, int Position)
{
    /// <summary>Order by sequence index, then position.</summary>
    internal static int Compare(OligoOccurrence a, OligoOccurrence b)
    {
        int c = a.SequenceIndex.CompareTo(b.SequenceIndex);
        return c != 0 ? c : a.Position.CompareTo(b.Position);
    }
}

/// <summary>Best lexicon segmentation of a word (RSAT <c>segments</c> column): prefix, suffix and their maximal frequencies.</summary>
public readonly record struct LexiconSegmentation(string Prefix, string Suffix, double PrefixFrequency, double SuffixFrequency);

/// <summary>A word (or pair) scored by <see cref="MotifFinder.AnalyzeOligos"/>; NaN = RSAT "NA".</summary>
/// <param name="Pattern">Word (smaller pair member with both strands; may contain one IUPAC code with degeneracy).</param>
/// <param name="ReverseComplement">Pair partner with both strands; null for one strand.</param>
/// <param name="Occurrences">occ (pair sum with both strands).</param>
/// <param name="Overlaps">Discarded overlapping windows credited to the pattern (<c>-noov</c>; RSAT <c>ovl_occ</c>).</param>
/// <param name="Positions">Counted occurrences (sequence, window start), ascending.</param>
/// <param name="ObservedFrequency">occ / n (RSAT <c>obs_freq</c>).</param>
/// <param name="ExpectedFrequency">exp_freq.</param>
/// <param name="ExpectedOccurrences">exp_occ (NaN when the pattern is not tested, RSAT "NA").</param>
/// <param name="ExpectedVariance">exp_var.</param>
/// <param name="OverlapCoefficient">ovlp (RSAT <c>ovlp</c> column; with a calibration computed from the input composition, unused by the statistics).</param>
/// <param name="ZScore">(occ − exp_occ) / √exp_var (NaN when exp_var ≤ 0).</param>
/// <param name="Ratio">occ / exp_occ (0 when exp_occ = 0, as RSAT).</param>
/// <param name="OccurrenceProbability">occ_P (NaN when not tested).</param>
/// <param name="OccurrenceEValue">occ_E = occ_P · tested patterns.</param>
/// <param name="OccurrenceSignificance">occ_sig = −log₁₀ occ_E.</param>
/// <param name="FittedDistribution">Distribution behind occ_P.</param>
/// <param name="LexiconSegmentation">Best segmentation with the lexicon background (plain words), otherwise null.</param>
public sealed record OligoStatistics(
    string Pattern,
    string? ReverseComplement,
    int Occurrences,
    int Overlaps,
    IReadOnlyList<OligoOccurrence> Positions,
    double ObservedFrequency,
    double ExpectedFrequency,
    double ExpectedOccurrences,
    double ExpectedVariance,
    double OverlapCoefficient,
    double ZScore,
    double Ratio,
    double OccurrenceProbability,
    double OccurrenceEValue,
    double OccurrenceSignificance,
    OligoFittedDistribution FittedDistribution,
    LexiconSegmentation? LexiconSegmentation);

/// <summary>Result of <see cref="MotifFinder.AnalyzeOligos"/>.</summary>
/// <param name="Patterns">Patterns with occ ≥ the threshold (tested or not), first-occurrence order.</param>
/// <param name="OligoLength">k.</param>
/// <param name="Strands">Strand mode.</param>
/// <param name="CountOverlapping">-ovlp (true) or -noov (false).</param>
/// <param name="Degeneracy">Degenerate-word mode.</param>
/// <param name="SequenceCount">Number of input sequences (RSAT <c>sequence_number</c>).</param>
/// <param name="PossiblePositions">nb_pos = Σ (Lᵢ − k + 1) over sequences with Lᵢ ≥ k.</param>
/// <param name="TotalOccurrences">Binomial trials n (= nb_pos).</param>
/// <param name="TestedPatterns">Patterns with a P-value (the occ_E multiplier).</param>
/// <param name="PossibleOligos">RSAT <c>nb_possible_oligos</c> (the <c>-pseudo</c> denominator).</param>
public sealed record OligoAnalysisReport(
    IReadOnlyList<OligoStatistics> Patterns,
    int OligoLength,
    OligoStrandMode Strands,
    bool CountOverlapping,
    OligoDegeneracy Degeneracy,
    int SequenceCount,
    long PossiblePositions,
    long TotalOccurrences,
    int TestedPatterns,
    double PossibleOligos)
{
    /// <summary>Sequence type of the analysis (RSAT <c>-seqtype</c>).</summary>
    public OligoSequenceType SequenceType { get; init; } = OligoSequenceType.Dna;

    /// <summary>RSAT <c>alphabet_size</c>: 4 (DNA), 20 (protein) or the number of distinct residues of the sequences of length ≥ k (other).</summary>
    public int AlphabetSize { get; init; } = 4;
}

/// <summary>
/// RSAT <c>-lexicon</c> sub-word tables: relative frequencies of the w-mer prefixes (1 ≤ w &lt; k) of the counted k-mer
/// windows and the maximal segmentation frequencies M (see <see cref="OligoBackgroundModel.Lexicon"/>).
/// </summary>
internal sealed class OligoLexicon
{
    private readonly Dictionary<string, double>[] _max; // _max[w][x] = M(x) for observed sub-words of length w

    private OligoLexicon(Dictionary<string, double>[] max, int k)
    {
        _max = max;
        K = k;
    }

    public int K { get; }

    /// <summary>Lexicon of the overlapping k-mer windows of the sequences of length ≥ k.</summary>
    public static OligoLexicon FromWindows(IReadOnlyList<string> sequences, int k)
    {
        var counts = new Dictionary<string, long>(StringComparer.Ordinal);
        foreach (string s in sequences)
        {
            if (s.Length < k) continue;
            foreach (var (word, c) in s.AsSpan().CountKmersSpan(k))
                counts[word] = counts.GetValueOrDefault(word) + c;
        }
        return FromWordCounts(counts, k);
    }

    /// <summary>Lexicon from k-mer counts (RSAT occ + overlaps per word).</summary>
    public static OligoLexicon FromWordCounts(IEnumerable<KeyValuePair<string, long>> wordCounts, int k)
    {
        // RSAT CalcSubWordFrequencies: sub_word[w]{prefix} occurrences, relative to their (common) total.
        var occ = new Dictionary<string, long>[k];
        for (int w = 1; w < k; w++) occ[w] = new Dictionary<string, long>(StringComparer.Ordinal);
        foreach (var (word, c) in wordCounts)
        {
            if (c == 0) continue;
            for (int w = 1; w < k; w++)
            {
                string prefix = word[..w];
                occ[w][prefix] = occ[w].GetValueOrDefault(prefix) + c;
            }
        }

        var max = new Dictionary<string, double>[k];
        max[0] = new Dictionary<string, double>(StringComparer.Ordinal);
        for (int w = 1; w < k; w++)
        {
            double total = occ[w].Values.Sum();
            max[w] = occ[w].ToDictionary(e => e.Key, e => e.Value / total, StringComparer.Ordinal);
        }

        var lexicon = new OligoLexicon(max, k);
        // Levels 2 … k − 1 in increasing length: M(x) = max(f(x), best segmentation).
        for (int l = 2; l < k; l++)
        {
            var level = max[l];
            foreach (string x in level.Keys.ToList())
            {
                double best = lexicon.BestSegmentation(x).Frequency;
                level[x] = Math.Max(level[x], best);
            }
        }
        return lexicon;
    }

    private double MaxFrequency(string x) => _max[x.Length].TryGetValue(x, out double v) ? v : 0.0;

    // max over s of M(x[..s]) · M(x[s..]) with the first strictly larger split kept (RSAT '>' from 0).
    private (double Frequency, int Split) BestSegmentation(string x)
    {
        double best = 0;
        int split = 0;
        for (int s = 1; s < x.Length; s++)
        {
            double seg = MaxFrequency(x[..s]) * MaxFrequency(x[s..]);
            if (seg > best)
            {
                best = seg;
                split = s;
            }
        }
        return (best, split);
    }

    /// <summary>exp_freq of a k-mer and its best segmentation (null when every segmentation has frequency 0).</summary>
    public (double Frequency, LexiconSegmentation? Segmentation) ExpectedFrequency(string word)
    {
        var (best, split) = BestSegmentation(word);
        if (split == 0) return (0.0, null);
        string prefix = word[..split], suffix = word[split..];
        return (best, new LexiconSegmentation(prefix, suffix, MaxFrequency(prefix), MaxFrequency(suffix)));
    }
}
