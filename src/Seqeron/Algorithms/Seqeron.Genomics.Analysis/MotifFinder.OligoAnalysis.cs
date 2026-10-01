using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// RSAT <c>oligo-analysis</c> word statistics (van Helden, André &amp; Collado-Vides 1998, J Mol Biol 281:827–842):
/// binomial occurrence significance (occ_P / occ_E / occ_sig), matching-sequence significance (ms_P / ms_E / ms_sig),
/// Bernoulli and Markov background models, and single- or both-strand (<c>-1str</c> / <c>-2str</c>) counting.
/// Ported from rsa-tools/rsat-code (<c>perl-scripts/oligo-analysis</c>, <c>perl-scripts/lib/RSA.disco.lib</c>,
/// <c>perl-scripts/lib/RSA.lib</c>, <c>perl-scripts/lib/RSAT/stats.pm</c>, <c>perl-scripts/lib/RSAT/MarkovModel.pm</c>).
/// </summary>
public static partial class MotifFinder
{
    #region RSAT oligo-analysis significance

    private static readonly double Ln4 = Math.Log(4.0);
    private static readonly double Ln10 = Math.Log(10.0);

    /// <summary>
    /// RSAT <c>oligo-analysis</c> over-representation analysis of the length-<paramref name="k"/> words of one sequence:
    /// counts, expected frequency under <paramref name="background"/>, observed/expected ratio and the binomial
    /// right-tail significance occ_P / occ_E / occ_sig.
    /// </summary>
    /// <remarks>
    /// <para>For every word w with occ(w) ≥ <paramref name="minCount"/> (RSAT <c>-lth occ</c>), with n = the total number
    /// of word occurrences (RSAT <c>sum_occurrences</c>: the N − k + 1 windows, overlapping or not) and p = exp_freq(w):</para>
    /// <list type="bullet">
    /// <item>exp_occ = p · n; ratio = occ / exp_occ;</item>
    /// <item>occ_P = P(X ≥ occ), X ~ Binomial(n, p) (RSAT <c>sum_of_binomials</c>, right tail);</item>
    /// <item>occ_E = occ_P · T, where T = <see cref="OligoAnalysisResult.TestedPatterns"/> is the number of words (pairs
    /// with <c>-2str</c>) that passed the occurrence threshold — RSAT <c>MultiTestCorrections($nb_tested_patterns)</c>;
    /// <see cref="OligoAnalysisResult.PossibleOligos"/> gives the number of possible oligos used by the manual's
    /// definition (equal to T with RSAT <c>-zeroocc</c>);</item>
    /// <item>occ_sig = −log₁₀ occ_E.</item>
    /// </list>
    /// <para>With <see cref="OligoStrandMode.Both"/> (<c>-2str</c>) a word and its reverse complement form one pattern:
    /// occ(W|W') = occ(W) + occ(W'), exp_freq(W|W') = exp_freq(W) + exp_freq(W') (RSAT manual, "Specific treatment for
    /// double strand counts"); a reverse palindrome is counted once with its single-strand expectation. The pattern is
    /// reported under the lexicographically smaller member (RSAT <c>GroupRC</c>). Note: current RSAT code copies the kept
    /// member's expected frequency to its partner (2·exp_freq(min(W, W'))), which equals the documented sum only for a
    /// strand-symmetric background (equiprobable, input Bernoulli); the documented sum is implemented here.</para>
    /// <para><paramref name="countOverlapping"/> = false is RSAT <c>-noov</c>: scanning from the 3′ end, an occurrence
    /// overlapping a previously counted occurrence of the same word (or, with <c>-2str</c>, of its reverse complement) is
    /// discarded; n still counts all windows.</para>
    /// <para>Probabilities are computed in log space (<see cref="StatisticsHelper.LogBinomialUpperTail"/>), so occ_sig stays
    /// exact where RSAT's occ_P underflows (RSAT then prints its cap of 350).</para>
    /// Results are in order of each pattern's first occurrence; positions are 0-based window starts on the given strand.
    /// </remarks>
    /// <param name="sequence">DNA sequence.</param>
    /// <param name="k">Oligonucleotide length (≥ 1).</param>
    /// <param name="minCount">Occurrence threshold (RSAT <c>-lth occ</c>); values ≤ 1 test every observed word.</param>
    /// <param name="background">Background model (RSAT <c>-bg</c> / <c>-markov</c> / <c>-bgfile</c>).</param>
    /// <param name="strands">Single strand (<c>-1str</c>) or both strands grouped by reverse-complement pairs (<c>-2str</c>).</param>
    /// <param name="countOverlapping">true = RSAT <c>-ovlp</c> (default), false = <c>-noov</c>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="background"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> &lt; 1, or a Markov order not allowed for this k.</exception>
    /// <exception cref="ArgumentException">The background assigns probability 0 to an observed word.</exception>
    public static OligoAnalysisResult DiscoverMotifs(
        DnaSequence sequence,
        int k,
        int minCount,
        OligoBackgroundModel background,
        OligoStrandMode strands = OligoStrandMode.Single,
        bool countOverlapping = true)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(background);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        background.ValidateFor(k, nameof(background));

        bool both = strands == OligoStrandMode.Both;
        string seq = sequence.Sequence;
        long windows = Math.Max(0, seq.Length - k + 1);

        Dictionary<string, List<int>> wordPositions = countOverlapping
            ? CollectKmerPositions(seq, k, 1)
            : CollectNonOverlappingPositions(seq, k, both);

        var patterns = GroupPatterns(wordPositions, both);
        var tested = patterns.Where(p => p.Positions.Count >= minCount).ToList();

        var motifs = new List<SignificantMotif>(tested.Count);
        if (tested.Count > 0)
        {
            Func<string, double> logWordProbability = background.CreateLogProbability(new[] { seq }, k, both);
            double logN = Math.Log(windows);
            double logTested = Math.Log(tested.Count);
            foreach (var pattern in tested)
            {
                int count = pattern.Positions.Count;
                double logP = PatternLogProbability(pattern, logWordProbability);
                double logOccP = StatisticsHelper.LogBinomialUpperTail(count, windows, logP);
                double logOccE = logOccP + logTested;
                motifs.Add(new SignificantMotif(
                    Sequence: pattern.Word,
                    ReverseComplement: both ? pattern.ReverseComplement : null,
                    Count: count,
                    Positions: pattern.Positions.AsReadOnly(),
                    ExpectedFrequency: Math.Exp(logP),
                    ExpectedOccurrences: Math.Exp(logP + logN),
                    Ratio: Math.Exp(Math.Log(count) - logN - logP),
                    OccurrenceProbability: Math.Exp(logOccP),
                    OccurrenceEValue: Math.Exp(logOccE),
                    OccurrenceSignificance: -logOccE / Ln10));
            }
        }

        return new OligoAnalysisResult(
            Motifs: motifs.AsReadOnly(),
            OligoLength: k,
            Strands: strands,
            CountOverlapping: countOverlapping,
            TotalOccurrences: windows,
            TestedPatterns: tested.Count,
            PossibleOligos: PossibleOligos(k, both));
    }

    /// <summary>
    /// Shared words with RSAT <c>oligo-analysis</c> matching-sequence significance (<c>-return mseq,proba</c>): for each
    /// word (or reverse-complement pair with <c>-2str</c>) present in at least <paramref name="minSequences"/> sequences
    /// (RSAT <c>-lth mseq</c>), the binomial probability ms_P of observing that many matching sequences.
    /// </summary>
    /// <remarks>
    /// RSAT <c>CalcExpected</c> / <c>CalcProba</c>: with S input sequences, nb_pos = Σ (Lᵢ − k + 1) over sequences with
    /// Lᵢ ≥ k, positions per sequence π = nb_pos / S and p = exp_freq(w), the probability that one sequence contains
    /// the word is P₁ = 1 − (1 − p)^π (evaluated as −expm1(π·log1p(−p))); exp_ms = S·P₁;
    /// ms_P = P(X ≥ mseq), X ~ Binomial(S, P₁); ms_E = ms_P · NPO with NPO = <see cref="SharedMotifAnalysisResult.PossibleOligos"/>
    /// (4^k, or (4^k + 4^(k/2))/2 pairs with <c>-2str</c>; RSAT <c>NbPossibleOligos</c>); ms_sig = −log₁₀ ms_E.
    /// With <c>-2str</c> a sequence matches the pair W|W' when it contains W or W'. Input-estimated backgrounds are
    /// estimated from all sequences of length ≥ k. Results are in order of first occurrence (sequence, position).
    /// </remarks>
    /// <param name="sequences">Input sequences (no null elements).</param>
    /// <param name="k">Word length (≥ 1).</param>
    /// <param name="minSequences">Matching-sequence quorum (≥ 1).</param>
    /// <param name="background">Background model.</param>
    /// <param name="strands">Single strand or reverse-complement pairs.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> or <paramref name="background"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> &lt; 1, <paramref name="minSequences"/> &lt; 1, or an invalid Markov order.</exception>
    /// <exception cref="ArgumentException">A sequence is null, or the background assigns probability 0 to an observed word.</exception>
    public static SharedMotifAnalysisResult FindSharedMotifs(
        IEnumerable<DnaSequence> sequences,
        int k,
        int minSequences,
        OligoBackgroundModel background,
        OligoStrandMode strands = OligoStrandMode.Single)
        => FindSharedMotifsCore(sequences, k, minSequences, background, strands, null);

    // Shared body; adjustWordProbability (RSAT -pseudo) maps the background's ln exp_freq(word) before the pair sum.
    private static SharedMotifAnalysisResult FindSharedMotifsCore(
        IEnumerable<DnaSequence> sequences,
        int k,
        int minSequences,
        OligoBackgroundModel background,
        OligoStrandMode strands,
        Func<Func<string, double>, Func<string, double>>? adjustWordProbability)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(background);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(minSequences, 1);
        background.ValidateFor(k, nameof(background));

        bool both = strands == OligoStrandMode.Both;
        var seqs = new List<string>();
        foreach (var dna in sequences)
        {
            if (dna is null)
                throw new ArgumentException($"Sequence at index {seqs.Count} is null.", nameof(sequences));
            seqs.Add(dna.Sequence);
        }

        long possiblePositions = 0;
        var rcCache = new Dictionary<string, string>();
        var matching = new Dictionary<string, (string Rc, List<int> Indices)>();
        for (int i = 0; i < seqs.Count; i++)
        {
            string s = seqs[i];
            if (s.Length < k) continue;
            possiblePositions += s.Length - k + 1;
            // Distinct words of this sequence via the canonical counter (keys in first-occurrence order).
            foreach (string word in s.AsSpan().CountKmersSpan(k).Keys)
            {
                string key = word, rc = word;
                if (both)
                {
                    rc = CachedReverseComplement(word, rcCache);
                    if (string.CompareOrdinal(rc, word) < 0) (key, rc) = (rc, word);
                }

                if (!matching.TryGetValue(key, out var entry))
                {
                    entry = (rc, new List<int>());
                    matching.Add(key, entry);
                }

                if (entry.Indices.Count == 0 || entry.Indices[^1] != i)
                    entry.Indices.Add(i);
            }
        }

        double logNpo = LogPossibleOligos(k, both);
        var motifs = new List<SignificantSharedMotif>();
        var shared = matching.Where(e => e.Value.Indices.Count >= minSequences).ToList();
        if (shared.Count > 0)
        {
            Func<string, double> logWordProbability = background.CreateLogProbability(seqs, k, both);
            if (adjustWordProbability is not null)
                logWordProbability = adjustWordProbability(logWordProbability);
            int sequenceCount = seqs.Count;
            double positionsPerSequence = (double)possiblePositions / sequenceCount;
            foreach (var (word, (rc, indices)) in shared)
            {
                var pattern = new OligoPattern(word, rc, new List<int>());
                double logP = PatternLogProbability(pattern, logWordProbability);
                motifs.Add(SharedMotifStatistics(word, both ? rc : null, indices, sequenceCount, positionsPerSequence, logP, logNpo));
            }
        }

        return new SharedMotifAnalysisResult(
            Motifs: motifs.AsReadOnly(),
            OligoLength: k,
            Strands: strands,
            SequenceCount: seqs.Count,
            PossiblePositions: possiblePositions,
            PossibleOligos: PossibleOligos(k, both));
    }

    /// <summary>
    /// RSAT <c>CalcExpected</c> / <c>CalcProba</c> matching-sequence statistics of one pattern with ln exp_freq
    /// <paramref name="logP"/>: P₁ = 1 − (1 − p)^π, exp_ms = S·P₁, ms_P = P(X ≥ mseq) with X ~ Bin(S, P₁), ms_E = ms_P·NPO.
    /// </summary>
    private static SignificantSharedMotif SharedMotifStatistics(string word, string? reverseComplement, List<int> indices,
        int sequenceCount, double positionsPerSequence, double logP, double logNpo)
    {
        double p = Math.Exp(logP);
        double oneSequence = -StatisticsHelper.ExpM1(positionsPerSequence * StatisticsHelper.Log1P(-p));
        double logOneSequence = oneSequence > 0
            ? Math.Log(oneSequence)
            : Math.Log(positionsPerSequence) + logP; // P₁ ≈ π·p when p underflows
        double logMsP = StatisticsHelper.LogBinomialUpperTail(indices.Count, sequenceCount, logOneSequence);
        double logMsE = logMsP + logNpo;
        return new SignificantSharedMotif(
            Sequence: word,
            ReverseComplement: reverseComplement,
            SequenceIndices: indices.AsReadOnly(),
            Prevalence: (double)indices.Count / sequenceCount,
            ExpectedFrequency: p,
            ExpectedMatchingSequences: sequenceCount * oneSequence,
            MatchingSequenceProbability: Math.Exp(logMsP),
            MatchingSequenceEValue: Math.Exp(logMsE),
            MatchingSequenceSignificance: -logMsE / Ln10);
    }

    /// <summary>A word (or reverse-complement pair) with its 0-based positions on the given strand.</summary>
    private sealed record OligoPattern(string Word, string ReverseComplement, List<int> Positions);

    /// <summary>ln exp_freq of a pattern: single word, palindrome, or W|W' pair (ln(p(W) + p(W'))).</summary>
    private static double PatternLogProbability(OligoPattern pattern, Func<string, double> logWordProbability)
    {
        double logP = logWordProbability(pattern.Word);
        if (pattern.ReverseComplement != pattern.Word)
            logP = LogAddExp(logP, logWordProbability(pattern.ReverseComplement));
        if (double.IsNegativeInfinity(logP))
            throw new ArgumentException(
                $"The background model assigns probability 0 to the observed word {pattern.Word}.");
        return Math.Min(logP, 0.0);
    }

    private static double LogAddExp(double a, double b)
    {
        if (double.IsNegativeInfinity(a)) return b;
        if (double.IsNegativeInfinity(b)) return a;
        double max = Math.Max(a, b);
        return max + Math.Log(Math.Exp(a - max) + Math.Exp(b - max));
    }

    /// <summary>
    /// RSAT <c>NbPossibleOligos</c>: 4^k (single strand) or 4^k − (4^k − P)/2 = (4^k + P)/2 with P = 4^(k/2) reverse
    /// palindromes for even k, 0 for odd k (both strands). Exact powers of two; +∞ once 4^k exceeds the double range.
    /// </summary>
    private static double PossibleOligos(int k, bool both)
    {
        double all = Math.ScaleB(1.0, 2 * k);
        if (!both) return all;
        double palindromes = k % 2 == 0 ? Math.ScaleB(1.0, k) : 0.0;
        return (all + palindromes) / 2;
    }

    /// <summary>ln <see cref="PossibleOligos"/>, finite for every k.</summary>
    private static double LogPossibleOligos(int k, bool both)
    {
        double npo = PossibleOligos(k, both);
        if (double.IsFinite(npo)) return Math.Log(npo);
        double logHalf = k * Ln4 - Math.Log(2.0);
        return both ? logHalf : k * Ln4; // 4^(k/2) is negligible against 4^k here (k ≥ 512)
    }

    /// <summary>
    /// Groups single-strand words (first-occurrence order, ascending positions) into patterns: identity with one strand;
    /// reverse-complement pairs keyed by the lexicographically smaller member with merged positions with both strands.
    /// </summary>
    private static List<OligoPattern> GroupPatterns(Dictionary<string, List<int>> wordPositions, bool both)
    {
        var result = new List<OligoPattern>(wordPositions.Count);
        if (!both)
        {
            foreach (var (word, positions) in wordPositions)
                result.Add(new OligoPattern(word, word, positions));
            return result;
        }

        var byKey = new Dictionary<string, OligoPattern>();
        foreach (var (word, positions) in wordPositions)
        {
            string rc = DnaSequence.GetReverseComplementString(word);
            string key = string.CompareOrdinal(rc, word) < 0 ? rc : word;
            if (byKey.TryGetValue(key, out var existing))
            {
                existing.Positions.AddRange(positions);
                existing.Positions.Sort();
                continue;
            }

            var pattern = new OligoPattern(key, key == word ? rc : word, new List<int>(positions));
            byKey.Add(key, pattern);
            result.Add(pattern);
        }

        return result;
    }

    /// <summary>
    /// RSAT <c>-noov</c> counting (oligo-analysis <c>CountOligos</c>): windows are read from the 3′ end; a window is
    /// discarded when the last counted occurrence of the same word — or, with both strands, of its reverse complement —
    /// starts fewer than k positions downstream. Returns counted positions per word, first-occurrence order, ascending.
    /// </summary>
    private static Dictionary<string, List<int>> CollectNonOverlappingPositions(string seq, int k, bool both)
    {
        var lastCounted = new Dictionary<string, int>();
        var descending = new Dictionary<string, List<int>>();
        var rcCache = new Dictionary<string, string>();
        for (int pos = seq.Length - k; pos >= 0; pos--)
        {
            string word = seq.Substring(pos, k);
            if (lastCounted.TryGetValue(word, out int last) && last - pos < k)
                continue;
            if (both && lastCounted.TryGetValue(CachedReverseComplement(word, rcCache), out int lastRc) && lastRc - pos < k)
                continue;

            lastCounted[word] = pos;
            if (!descending.TryGetValue(word, out var list))
            {
                list = new List<int>();
                descending.Add(word, list);
            }
            list.Add(pos);
        }

        var result = new Dictionary<string, List<int>>(descending.Count);
        foreach (var (word, list) in descending.OrderBy(e => e.Value[^1]))
        {
            list.Reverse();
            result.Add(word, list);
        }

        return result;
    }

    private static string CachedReverseComplement(string word, Dictionary<string, string> cache)
    {
        if (!cache.TryGetValue(word, out var rc))
        {
            rc = DnaSequence.GetReverseComplementString(word);
            cache.Add(word, rc);
        }
        return rc;
    }

    #endregion
}

/// <summary>Strand handling of RSAT <c>oligo-analysis</c>.</summary>
public enum OligoStrandMode
{
    /// <summary>Given strand only (RSAT <c>-1str</c>).</summary>
    Single,

    /// <summary>Both strands; each word grouped with its reverse complement (RSAT <c>-2str</c>, <c>-grouprc</c>).</summary>
    Both,
}

/// <summary>
/// Background model giving the expected frequency of a word for RSAT <c>oligo-analysis</c> statistics.
/// </summary>
public sealed class OligoBackgroundModel
{
    private enum ModelKind { Equiprobable, Bernoulli, BernoulliFromInput, MarkovFromInput, MarkovTable, Lexicon }

    private readonly ModelKind _kind;
    private readonly double[]? _residues;                 // Bernoulli: normalised A, C, G, T
    private readonly int _order;                          // Markov order (0 for Bernoulli kinds)
    private readonly Dictionary<string, double>? _logPrefix;      // Markov table: ln P(prefix) for observed prefixes
    private readonly Dictionary<string, double[]>? _logTransition; // Markov table: ln P(suffix | prefix), A, C, G, T
    private readonly double _logAbsentPrefix;             // Markov table: ln P(prefix) for prefixes absent from the table
    private readonly double _logAbsentTransition;         // Markov table: ln P(suffix | absent prefix) = ln ¼

    private OligoBackgroundModel(ModelKind kind, int order = 0, double[]? residues = null,
        Dictionary<string, double>? logPrefix = null, Dictionary<string, double[]>? logTransition = null,
        double logAbsentPrefix = double.NegativeInfinity, double logAbsentTransition = double.NegativeInfinity)
    {
        _kind = kind;
        _order = order;
        _residues = residues;
        _logPrefix = logPrefix;
        _logTransition = logTransition;
        _logAbsentPrefix = logAbsentPrefix;
        _logAbsentTransition = logAbsentTransition;
    }

    /// <summary>Markov order of the model (0 for the Bernoulli / equiprobable models).</summary>
    public int MarkovOrder => _order;

    /// <summary>Equiprobable residues, exp_freq = 4^−k (RSAT <c>-bg equi</c>).</summary>
    public static OligoBackgroundModel Equiprobable { get; } = new(ModelKind.Equiprobable);

    /// <summary>
    /// Bernoulli model with residue frequencies estimated from the input sequences (RSAT default, <c>-bg input</c>,
    /// <c>CalcAlphabet</c>): q(b) = count(b) / total; with both strands q(b) = (count(b) + count(b̄)) / (2·total).
    /// </summary>
    public static OligoBackgroundModel BernoulliFromInput { get; } = new(ModelKind.BernoulliFromInput);

    /// <summary>
    /// Bernoulli model with given residue probabilities (A, C, G, T; finite, strictly positive, normalised to sum 1):
    /// exp_freq(w) = ∏ q(wᵢ).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="acgt"/> is null.</exception>
    /// <exception cref="ArgumentException">Not exactly 4 values.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A value is not finite and positive.</exception>
    public static OligoBackgroundModel Bernoulli(IReadOnlyList<double> acgt)
    {
        ArgumentNullException.ThrowIfNull(acgt);
        return new OligoBackgroundModel(ModelKind.Bernoulli, residues: MotifFinder.NormalizeBackground(acgt));
    }

    /// <summary>
    /// Markov chain of order m estimated from the input sequences (RSAT <c>-markov m</c>, slow-count mode):
    /// exp_freq(w) = f(w[0..m]) · ∏_{o=1}^{k−m−1} f(w[o..o+m]) / f(w[o..o+m−1]), where f are the relative frequencies of
    /// the overlapping (m+1)-mers and m-mers of the input (single strand, sequences of length ≥ k; oligo-analysis
    /// <c>CalcSubWordFrequencies</c> / <c>CalcExpected</c>). Order 0 is the single-strand residue composition.
    /// RSAT requires m ≤ k − 2 for m &gt; 0; this is checked when the model is used.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="order"/> &lt; 0.</exception>
    public static OligoBackgroundModel MarkovFromInput(int order)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(order);
        return new OligoBackgroundModel(ModelKind.MarkovFromInput, order);
    }

    /// <summary>
    /// Markov chain loaded from a table of (m+1)-mer frequencies (RSAT <c>-bgfile</c>, <c>oligos</c> format;
    /// <c>RSAT::MarkovModel</c> <c>load_from_file_oligos</c> → <c>add_pseudo_freq</c> → <c>normalize_transition_frequencies</c>
    /// → <c>segment_proba</c>). With pseudo-frequency ψ, S(x) the table sum of prefix x and F the table total:
    /// P(b | x) = (1 − ψ)·f(xb)/S(x) + ψ/4 (¼ when S(x) = 0), P(x) = (1 − ψ)·S(x)/F + ψ/4^m, and
    /// exp_freq(w) = P(w[0..m−1]) · ∏_{c=m}^{k−1} P(w_c | w[c−m..c−1]). The table scale is irrelevant.
    /// RSAT's <c>segment_proba</c> returns the value printed with <c>%5g</c> (6 significant digits); the unrounded value
    /// is used here.
    /// </summary>
    /// <param name="frequencies">(m+1)-mer → frequency (A/C/G/T words of one length, finite, ≥ 0, positive total; missing = 0).</param>
    /// <param name="pseudoFrequency">ψ in [0, 1] (RSAT <c>bg_pseudo</c>, default 0.01).</param>
    /// <param name="strandInsensitive">
    /// The table holds strand-insensitive pair frequencies (RSAT "2str" file: a missing reverse complement takes its
    /// partner's value, then every non-palindromic frequency is halved).
    /// </param>
    /// <exception cref="ArgumentNullException"><paramref name="frequencies"/> is null.</exception>
    /// <exception cref="ArgumentException">Empty table, mixed lengths, non-ACGT word, or zero total.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A negative / non-finite frequency, or ψ outside [0, 1].</exception>
    public static OligoBackgroundModel MarkovFromOligoFrequencies(
        IReadOnlyDictionary<string, double> frequencies,
        double pseudoFrequency = 0.01,
        bool strandInsensitive = false)
    {
        ArgumentNullException.ThrowIfNull(frequencies);
        if (!(pseudoFrequency >= 0.0 && pseudoFrequency <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(pseudoFrequency), pseudoFrequency, "Pseudo-frequency must be in [0, 1].");
        if (frequencies.Count == 0)
            throw new ArgumentException("The frequency table is empty.", nameof(frequencies));

        int length = -1;
        var table = new Dictionary<string, double>(StringComparer.Ordinal);
        foreach (var (rawWord, value) in frequencies)
        {
            if (string.IsNullOrEmpty(rawWord))
                throw new ArgumentException("Empty word in the frequency table.", nameof(frequencies));
            string word = rawWord.ToUpperInvariant();
            if (length < 0) length = word.Length;
            if (word.Length != length)
                throw new ArgumentException("All words of a Markov frequency table must have the same length.", nameof(frequencies));
            foreach (char c in word)
            {
                if (MotifFinder.AcgtIndex(c) < 0)
                    throw new ArgumentException($"Word '{rawWord}' contains a non-ACGT character.", nameof(frequencies));
            }
            if (!double.IsFinite(value) || value < 0)
                throw new ArgumentOutOfRangeException(nameof(frequencies), value, "Frequencies must be finite and non-negative.");
            if (!table.TryAdd(word, value))
                throw new ArgumentException($"Duplicate word '{rawWord}' in the frequency table.", nameof(frequencies));
        }

        if (strandInsensitive)
        {
            foreach (string word in table.Keys.Order(StringComparer.Ordinal).ToList())
            {
                string rc = DnaSequence.GetReverseComplementString(word);
                if (!table.ContainsKey(rc) && table[word] != 0)
                    table[rc] = table[word];
            }
            foreach (string word in table.Keys.ToList())
            {
                if (DnaSequence.GetReverseComplementString(word) != word)
                    table[word] /= 2;
            }
        }

        int order = length - 1;
        var prefixSum = new Dictionary<string, double>(StringComparer.Ordinal);
        double total = 0;
        foreach (var (word, value) in table)
        {
            string prefix = word[..order];
            prefixSum[prefix] = prefixSum.GetValueOrDefault(prefix) + value;
            total += value;
        }
        if (!(total > 0))
            throw new ArgumentException("The frequency table has a zero total.", nameof(frequencies));

        double psi = pseudoFrequency;
        double uniformPrefix = Math.Pow(4.0, -order);
        var logPrefix = new Dictionary<string, double>(StringComparer.Ordinal);
        var logTransition = new Dictionary<string, double[]>(StringComparer.Ordinal);
        foreach (var (prefix, sum) in prefixSum)
        {
            logPrefix[prefix] = Math.Log((1 - psi) * sum / total + psi * uniformPrefix);
            var row = new double[4];
            for (int b = 0; b < 4; b++)
            {
                double f = table.GetValueOrDefault(prefix + MotifFinder.AcgtBases[b]);
                // RSAT: an all-zero prefix row gets transitions ¼ from the pseudo-frequency (0 without one).
                double transition = psi > 0 ? 0.25 : 0.0;
                if (sum > 0) transition = (1 - psi) * f / sum + psi / 4;
                row[b] = Math.Log(transition);
            }
            logTransition[prefix] = row;
        }

        // The same probabilities in linear form, indexed by base-4 prefix code (A = 0 … T = 3, first letter most
        // significant), for the exact PWM p-value DP (MotifFinder.PwmMarkovScorePValue).
        int contexts = 1 << (2 * order);
        var initial = new double[contexts];
        var transitions = new double[contexts * 4];
        double absentTransition = psi > 0 ? 0.25 : 0.0;
        var prefixChars = new char[order];
        for (int x = 0; x < contexts; x++)
        {
            for (int i = 0, code = x; i < order; i++, code >>= 2)
                prefixChars[order - 1 - i] = MotifFinder.AcgtBases[code & 3];
            string prefix = new(prefixChars);
            if (prefixSum.TryGetValue(prefix, out double sum))
            {
                initial[x] = (1 - psi) * sum / total + psi * uniformPrefix;
                for (int b = 0; b < 4; b++)
                {
                    double f = table.GetValueOrDefault(prefix + MotifFinder.AcgtBases[b]);
                    double transition = psi > 0 ? 0.25 : 0.0;
                    if (sum > 0) transition = (1 - psi) * f / sum + psi / 4;
                    transitions[x * 4 + b] = transition;
                }
            }
            else
            {
                initial[x] = psi * uniformPrefix;
                for (int b = 0; b < 4; b++)
                    transitions[x * 4 + b] = absentTransition;
            }
        }

        return new OligoBackgroundModel(ModelKind.MarkovTable, order,
            logPrefix: logPrefix, logTransition: logTransition,
            logAbsentPrefix: Math.Log(psi * uniformPrefix),
            logAbsentTransition: psi > 0 ? Math.Log(0.25) : double.NegativeInfinity)
        {
            _chainInitial = initial,
            _chainTransitions = transitions,
        };
    }

    /// <summary>Markov table: P(prefix) and P(b | prefix) in linear form (base-4 prefix codes).</summary>
    private double[]? _chainInitial;
    private double[]? _chainTransitions;

    /// <summary>Largest Markov order accepted by the PWM p-value DP (4^10 contexts).</summary>
    internal const int MaxPValueMarkovOrder = 10;

    /// <summary>
    /// The model as an explicit word distribution for the exact PWM p-value DP: equiprobable / Bernoulli → i.i.d.;
    /// Markov table → order-m chain with RSAT's <c>segment_proba</c> prefix and transition probabilities.
    /// </summary>
    /// <exception cref="ArgumentException">An input-estimated or lexicon model (not an explicit probability distribution), or order &gt; 10.</exception>
    internal PwmBackgroundChain ToPValueChain(string paramName)
    {
        switch (_kind)
        {
            case ModelKind.Equiprobable:
                return new PwmBackgroundChain(4, new[] { 0.25, 0.25, 0.25, 0.25 });
            case ModelKind.Bernoulli:
                return new PwmBackgroundChain(4, (double[])_residues!.Clone());
            case ModelKind.MarkovTable:
                if (_order > MaxPValueMarkovOrder)
                    throw new ArgumentException(
                        $"Markov order {_order} exceeds the maximum ({MaxPValueMarkovOrder}) supported by the p-value DP.", paramName);
                return _order == 0
                    ? new PwmBackgroundChain(4, (double[])_chainTransitions!.Clone())
                    : new PwmBackgroundChain(4, _order, _chainInitial!, _chainTransitions!);
            default:
                throw new ArgumentException(
                    "PWM p-values need an explicit word distribution: use Equiprobable, Bernoulli(acgt) or " +
                    "MarkovFromOligoFrequencies (input-estimated models need sequences and RSAT's input Markov estimate is " +
                    "not a normalised distribution; the lexicon model is a segmentation frequency, not a probability).",
                    paramName);
        }
    }

    /// <summary>
    /// RSAT <c>oligo-analysis -lexicon</c>: expected word frequencies from sub-word frequencies of the input, in the
    /// spirit of Bussemaker's dictionary segmentation. With f_w(x) the relative frequency of the w-mer x among the
    /// w-mer prefixes of the k-mer windows of the input (single strand, sequences of length ≥ k; RSAT
    /// <c>CalcSubWordFrequencies</c>, no trailing sub-words), the maximal segmentation frequency is
    /// M(x) = f₁(x) for a residue and, for an observed sub-word of length 2 ≤ l ≤ k − 1,
    /// M(x) = max(f_l(x), max_s M(x[0..s)) · M(x[s..l))) (0 for an unobserved sub-word); the expected frequency of a
    /// k-mer is exp_freq(w) = max_{1 ≤ s &lt; k} M(w[0..s)) · M(w[s..k)) (RSAT <c>CalcExpected</c>, bg_method "lexicon").
    /// Requires k ≥ 2. A word whose every segmentation contains an unobserved sub-word gets probability 0.
    /// </summary>
    public static OligoBackgroundModel Lexicon { get; } = new(ModelKind.Lexicon);

    /// <summary>True for <see cref="Lexicon"/>.</summary>
    internal bool IsLexicon => _kind == ModelKind.Lexicon;

    /// <summary>True for the models estimated from the input or alphabet-independent (equiprobable, input Bernoulli, input Markov, lexicon).</summary>
    internal bool IsAlphabetGeneric => _kind is ModelKind.Equiprobable or ModelKind.BernoulliFromInput
        or ModelKind.MarkovFromInput or ModelKind.Lexicon;

    /// <summary>True for the Markov models (RSAT leaves <c>%residue_proba</c> empty for them).</summary>
    internal bool IsMarkov => _kind is ModelKind.MarkovFromInput or ModelKind.MarkovTable;

    /// <summary>True for <see cref="Equiprobable"/>.</summary>
    internal bool IsEquiprobable => _kind == ModelKind.Equiprobable;

    /// <summary>
    /// Residue probabilities (A, C, G, T) RSAT <c>OverlapCoeff</c> uses for this model: the RSAT Bernoulli residue
    /// probabilities when the model is a Bernoulli one (equiprobable, given, or estimated from the input — the RSAT
    /// default, also behind <c>-lexicon</c>; pooled over complementary residues with both strands), otherwise (Markov
    /// models, for which RSAT leaves <c>%residue_proba</c> empty) ¼ each.
    /// </summary>
    internal double[] OverlapResidueProbabilities(IReadOnlyList<string> sequences, int k, bool bothStrands)
    {
        switch (_kind)
        {
            case ModelKind.Bernoulli:
                return (double[])_residues!.Clone();
            case ModelKind.BernoulliFromInput:
            case ModelKind.Lexicon:
                return InputResidueProbabilities(sequences, k, bothStrands);
            default:
                return new[] { 0.25, 0.25, 0.25, 0.25 };
        }
    }

    /// <summary>RSAT <c>CalcAlphabet</c>: residue frequencies of the sequences of length ≥ k (pooled with both strands).</summary>
    internal static double[] InputResidueProbabilities(IReadOnlyList<string> sequences, int k, bool bothStrands)
    {
        var counts = new double[4];
        foreach (string s in sequences)
        {
            if (s.Length < k) continue;
            foreach (var (residue, n) in s.AsSpan().CountKmersSpan(1))
            {
                int b = MotifFinder.AcgtIndex(residue[0]);
                if (b >= 0) counts[b] += n; // RSAT: residues outside the DNA alphabet are discarded
            }
        }
        double total = counts.Sum();
        var q = new double[4];
        for (int b = 0; b < 4; b++)
            q[b] = bothStrands ? (counts[b] + counts[3 - b]) / (2 * total) : counts[b] / total;
        return q;
    }

    /// <summary>Checks that the model can score words of length <paramref name="k"/> (RSAT order constraints).</summary>
    internal void ValidateFor(int k, string paramName)
    {
        if (_kind == ModelKind.Lexicon && k < 2)
            throw new ArgumentOutOfRangeException(paramName, k,
                "The lexicon background needs a word length of at least 2 (RSAT oligo-analysis -lexicon).");
        if (_kind == ModelKind.MarkovFromInput && _order > 0 && _order > k - 2)
            throw new ArgumentOutOfRangeException(paramName, _order,
                $"Markov order ({_order}) cannot be higher than word length - 2 ({k - 2}) (RSAT oligo-analysis).");
        if (_kind == ModelKind.MarkovTable && _order > k - 1)
            throw new ArgumentOutOfRangeException(paramName, _order,
                $"The word length ({k}) must be larger than the Markov order ({_order}) (RSAT segment_proba).");
    }

    /// <summary>ln exp_freq(word) for upper-case ACGT words of length k, estimated (where applicable) from <paramref name="sequences"/>.</summary>
    internal Func<string, double> CreateLogProbability(IReadOnlyList<string> sequences, int k, bool bothStrands)
    {
        switch (_kind)
        {
            case ModelKind.Equiprobable:
            {
                double logP = -k * Math.Log(4.0);
                return _ => logP;
            }
            case ModelKind.Bernoulli:
                return ResidueProduct(_residues!);
            case ModelKind.BernoulliFromInput:
                return ResidueProduct(InputResidueProbabilities(sequences, k, bothStrands));
            case ModelKind.MarkovFromInput:
                return MarkovFromInputLogProbability(sequences, k);
            case ModelKind.Lexicon:
            {
                var lexicon = OligoLexicon.FromWindows(sequences, k);
                return word => Math.Log(lexicon.ExpectedFrequency(word).Frequency);
            }
            default:
                return MarkovTableLogProbability;
        }
    }

    /// <summary>
    /// ln exp_freq(word) for words of an RSAT <c>-seqtype</c> alphabet (sequences prepared by
    /// <see cref="MotifFinder.AnalyzeOligoStrings"/>): equiprobable 1 / |A|^k; input Bernoulli over the residues of the
    /// alphabet; input Markov chains from RSAT's sub-word counts (<see cref="RsatSubWordLogFrequencies"/>); DNA-only
    /// models (given Bernoulli, Markov table) as for DNA.
    /// </summary>
    internal Func<string, double> CreateLogProbability(IReadOnlyList<string> sequences, int k, bool bothStrands, OligoResidueAlphabet alphabet)
    {
        switch (_kind)
        {
            case ModelKind.Equiprobable:
            {
                double logP = alphabet.IsDna ? -k * Math.Log(4.0) : -alphabet.LogPossibleOligos(k);
                return _ => logP;
            }
            case ModelKind.BernoulliFromInput:
            {
                if (alphabet.IsDna)
                    return ResidueProduct(InputResidueProbabilities(sequences, k, bothStrands));
                var logQ = alphabet.InputResidueProbabilities(sequences, k)
                    .ToDictionary(e => e.Key, e => Math.Log(e.Value));
                return word =>
                {
                    double sum = 0;
                    foreach (char c in word) sum += logQ.TryGetValue(c, out double l) ? l : double.NegativeInfinity;
                    return sum;
                };
            }
            case ModelKind.MarkovFromInput:
            {
                int m = _order;
                var longer = RsatSubWordLogFrequencies(sequences, k, m + 1, alphabet);
                var shorter = m > 0 ? RsatSubWordLogFrequencies(sequences, k, m, alphabet) : null;
                return word =>
                {
                    double sum = 0;
                    for (int o = 0; o + m < word.Length; o++)
                    {
                        if (!longer.TryGetValue(word.Substring(o, m + 1), out double lf))
                            return double.NegativeInfinity;
                        sum += lf;
                        if (m > 0 && o > 0) sum -= shorter![word.Substring(o, m)];
                    }
                    return sum;
                };
            }
            default:
                return CreateLogProbability(sequences, k, bothStrands);
        }
    }

    /// <summary>
    /// RSAT <c>CountOligos</c> + <c>CalcSubWordFrequencies</c> sub-word frequencies for <c>-markov</c>: ln(count / total) of the
    /// w-mers (w &lt; k) that are prefixes of the valid k-mer windows of the sequences of length ≥ k, plus the w-mers starting
    /// at L − k + 1 … L − w (RSAT's trailing sub-words, counted unfiltered; their total includes residues outside the
    /// alphabet). Without such residues these are all overlapping w-mers.
    /// </summary>
    private static Dictionary<string, double> RsatSubWordLogFrequencies(IReadOnlyList<string> sequences, int k, int w, OligoResidueAlphabet alphabet)
    {
        var counts = new Dictionary<string, long>(StringComparer.Ordinal);
        long total = 0;
        foreach (string s in sequences)
        {
            if (s.Length < k) continue;
            int[] invalid = alphabet.InvalidPrefixCounts(s);
            for (int pos = 0; pos + w <= s.Length; pos++)
            {
                bool windowPrefix = pos + k <= s.Length;
                if (windowPrefix && invalid[pos + k] != invalid[pos]) continue; // prefix of a discarded window
                string x = s.Substring(pos, w);
                counts[x] = counts.GetValueOrDefault(x) + 1;
                total++;
            }
        }
        double logTotal = Math.Log(total);
        return counts.ToDictionary(e => e.Key, e => Math.Log(e.Value) - logTotal, StringComparer.Ordinal);
    }

    private static Func<string, double> ResidueProduct(double[] q)
    {
        var logQ = q.Select(x => Math.Log(x)).ToArray();
        return word =>
        {
            double sum = 0;
            foreach (char c in word) sum += logQ[MotifFinder.AcgtIndex(c)];
            return sum;
        };
    }

    private Func<string, double> MarkovFromInputLogProbability(IReadOnlyList<string> sequences, int k)
    {
        int m = _order;
        var longer = RelativeWordLogFrequencies(sequences, k, m + 1);
        var shorter = m > 0 ? RelativeWordLogFrequencies(sequences, k, m) : null;
        return word =>
        {
            double sum = 0;
            for (int o = 0; o + m < word.Length; o++)
            {
                if (!longer.TryGetValue(word.Substring(o, m + 1), out double lf))
                    return double.NegativeInfinity; // unobserved (m+1)-mer: probability 0
                sum += lf;
                if (m > 0 && o > 0) sum -= shorter![word.Substring(o, m)];
            }
            return sum;
        };
    }

    // ln(count / total) of the overlapping w-mers of the sequences of length ≥ k (canonical k-mer counter).
    private static Dictionary<string, double> RelativeWordLogFrequencies(IReadOnlyList<string> sequences, int k, int w)
    {
        var counts = new Dictionary<string, long>(StringComparer.Ordinal);
        long total = 0;
        foreach (string s in sequences)
        {
            if (s.Length < k) continue;
            foreach (var (word, n) in s.AsSpan().CountKmersSpan(w))
            {
                counts[word] = counts.GetValueOrDefault(word) + n;
                total += n;
            }
        }
        double logTotal = Math.Log(total);
        return counts.ToDictionary(e => e.Key, e => Math.Log(e.Value) - logTotal, StringComparer.Ordinal);
    }

    private double MarkovTableLogProbability(string word)
    {
        int m = _order;
        double sum = 0;
        if (m > 0)
            sum = _logPrefix!.TryGetValue(word[..m], out double lp) ? lp : _logAbsentPrefix;
        for (int c = m; c < word.Length; c++)
        {
            int b = MotifFinder.AcgtIndex(word[c]);
            sum += _logTransition!.TryGetValue(word.Substring(c - m, m), out var row) ? row[b] : _logAbsentTransition;
        }
        return sum;
    }
}

/// <summary>
/// A word (or reverse-complement pair) scored by RSAT <c>oligo-analysis</c> occurrence statistics.
/// </summary>
/// <param name="Sequence">The word; with both strands the lexicographically smaller member of the pair.</param>
/// <param name="ReverseComplement">Pair partner with both strands (equal to <paramref name="Sequence"/> for a reverse palindrome); null for a single strand.</param>
/// <param name="Count">Occurrences (occ; pair sum with both strands).</param>
/// <param name="Positions">0-based start positions of the counted occurrences on the given strand, ascending.</param>
/// <param name="ExpectedFrequency">exp_freq (0 when it underflows the double range).</param>
/// <param name="ExpectedOccurrences">exp_occ = exp_freq · total occurrences.</param>
/// <param name="Ratio">occ / exp_occ.</param>
/// <param name="OccurrenceProbability">occ_P = P(X ≥ occ), X ~ Binomial(total occurrences, exp_freq).</param>
/// <param name="OccurrenceEValue">occ_E = occ_P · tested patterns.</param>
/// <param name="OccurrenceSignificance">occ_sig = −log₁₀ occ_E.</param>
public readonly record struct SignificantMotif(
    string Sequence,
    string? ReverseComplement,
    int Count,
    IReadOnlyList<int> Positions,
    double ExpectedFrequency,
    double ExpectedOccurrences,
    double Ratio,
    double OccurrenceProbability,
    double OccurrenceEValue,
    double OccurrenceSignificance);

/// <summary>Result of an RSAT <c>oligo-analysis</c> occurrence run.</summary>
/// <param name="Motifs">Tested patterns, in order of first occurrence.</param>
/// <param name="OligoLength">k.</param>
/// <param name="Strands">Strand mode.</param>
/// <param name="CountOverlapping">true for <c>-ovlp</c>, false for <c>-noov</c>.</param>
/// <param name="TotalOccurrences">Binomial trials n (RSAT <c>sum_occurrences</c> = N − k + 1 windows).</param>
/// <param name="TestedPatterns">Patterns tested (occ ≥ threshold), the occ_E multiplier.</param>
/// <param name="PossibleOligos">RSAT <c>nb_possible_oligos</c> (4^k, or pairs with both strands; +∞ beyond the double range).</param>
public sealed record OligoAnalysisResult(
    IReadOnlyList<SignificantMotif> Motifs,
    int OligoLength,
    OligoStrandMode Strands,
    bool CountOverlapping,
    long TotalOccurrences,
    int TestedPatterns,
    double PossibleOligos);

/// <summary>A word (or pair) scored by RSAT <c>oligo-analysis</c> matching-sequence statistics.</summary>
/// <param name="Sequence">The word (smaller pair member with both strands).</param>
/// <param name="ReverseComplement">Pair partner with both strands; null for a single strand.</param>
/// <param name="SequenceIndices">Ascending indices of the matching sequences (mseq = count).</param>
/// <param name="Prevalence">mseq / number of sequences.</param>
/// <param name="ExpectedFrequency">exp_freq.</param>
/// <param name="ExpectedMatchingSequences">exp_ms = S · (1 − (1 − exp_freq)^(nb_pos / S)).</param>
/// <param name="MatchingSequenceProbability">ms_P = P(X ≥ mseq), X ~ Binomial(S, exp_ms / S).</param>
/// <param name="MatchingSequenceEValue">ms_E = ms_P · nb_possible_oligos.</param>
/// <param name="MatchingSequenceSignificance">ms_sig = −log₁₀ ms_E.</param>
public readonly record struct SignificantSharedMotif(
    string Sequence,
    string? ReverseComplement,
    IReadOnlyList<int> SequenceIndices,
    double Prevalence,
    double ExpectedFrequency,
    double ExpectedMatchingSequences,
    double MatchingSequenceProbability,
    double MatchingSequenceEValue,
    double MatchingSequenceSignificance);

/// <summary>Result of an RSAT <c>oligo-analysis</c> matching-sequence run.</summary>
/// <param name="Motifs">Words meeting the quorum, in order of first occurrence.</param>
/// <param name="OligoLength">k.</param>
/// <param name="Strands">Strand mode.</param>
/// <param name="SequenceCount">S, all input sequences (RSAT <c>sequence_number</c>).</param>
/// <param name="PossiblePositions">nb_pos = Σ (Lᵢ − k + 1) over sequences with Lᵢ ≥ k.</param>
/// <param name="PossibleOligos">RSAT <c>nb_possible_oligos</c>, the ms_E multiplier.</param>
public sealed record SharedMotifAnalysisResult(
    IReadOnlyList<SignificantSharedMotif> Motifs,
    int OligoLength,
    OligoStrandMode Strands,
    int SequenceCount,
    long PossiblePositions,
    double PossibleOligos)
{
    /// <summary>Degenerate-word mode (RSAT <c>-oneN</c> / <c>-onedeg</c>); <see cref="OligoDegeneracy.None"/> for plain words.</summary>
    public OligoDegeneracy Degeneracy { get; init; } = OligoDegeneracy.None;
}
