using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// RSAT <c>dyad-analysis</c> (van Helden, Rios &amp; Collado-Vides 2000, Nucleic Acids Res 28:1808–1818): over-represented
/// spaced dyads — pairs of monads (oligonucleotides of length m) separated by a fixed spacing — with expected
/// frequencies from the monad frequencies of the input (or a background table) and binomial significance.
/// Ported from rsa-tools/rsat-code master 10043f2, <c>perl-scripts/dyad-analysis</c> v1.78 (<c>CountDyads</c>,
/// <c>CalcPossibleDyads</c>, <c>SecondElement</c>, <c>CalcOccSum</c>, <c>SumRCDyads</c>, <c>CalcMonadFrequencies</c>,
/// <c>CalcExpFreqFromMonads</c>, <c>CalcDyadFrequencies</c>, <c>CalcExpectedOcc</c>, <c>CalcProba</c>) with
/// <c>RSA.seq.lib</c> <c>OverlapCoeff</c> and <c>RSA.disco.lib</c> <c>MultiTestCorrections</c>.
/// </summary>
public static partial class MotifFinder
{
    #region RSAT dyad-analysis

    /// <summary>
    /// RSAT <c>dyad-analysis</c>: counts every dyad M₁ n{s} M₂ (monads of length m, s ∈ [min, max] unspecified residues
    /// between them) of the requested type and scores its over-representation.
    /// </summary>
    /// <remarks>
    /// <para>Per spacing s, with Lᵢ the sequence lengths: the possible positions are T_s = Σ max(0, Lᵢ − 2m − s + 1);
    /// a dyad occurs at p when M₁ = S[p, p + m) and M₂ = S[p + m + s, p + 2m + s). <see cref="DyadAnalysisOptions.CountOverlapping"/>
    /// = false (RSAT default <c>-noov</c>): scanning 5′ → 3′, an occurrence starting fewer than 2m + s positions after the
    /// last counted occurrence of the same dyad (or, with both strands, of its reverse complement) is discarded as an
    /// overlap. occ_sum_s / ovl_sum_s = the counted / discarded single-strand occurrences at spacing s.</para>
    /// <para>Both strands (RSAT default <c>-2str</c>): a dyad D and D' = rc(M₂) n{s} rc(M₁) form one pattern with
    /// occ = occ(D) + occ(D') (a reverse palindrome once), reported under the lexicographically smaller member.</para>
    /// <para>Expected frequency from monads (RSAT default, <c>-bg monads</c>): f(M) = occ(M) / Σ occ over all monad windows
    /// (single strand, overlapping); exp_freq(D) = f(M₁)·f(M₂), plus f(rc M₁)·f(rc M₂) with both strands unless D is a
    /// reverse palindrome (<c>CalcExpFreqFromMonads</c>). With <see cref="DyadBackgroundModel.DyadFrequencies"/> the table
    /// value (summed over D and D' with both strands); a missing value falls back to the monad formula when the
    /// maximal spacing exceeds 20 or the dyad has no occurrence, otherwise the dyad is not tested (RSAT).</para>
    /// <para>exp_occ = exp_freq · occ_sum_s (<c>CalcExpectedOcc</c>); obs_freq = occ / (occ_sum_s + ovl_sum_s);
    /// ratio = occ / exp_occ; ovlp = RSAT <c>OverlapCoeff</c> of the concatenated monads M₁M₂ with equiprobable residues;
    /// exp_var = T_s · p · (2·ovlp − 1 − (4m + 1)·p) and z = (occ − exp_occ)/√exp_var (NaN when exp_var ≤ 0);
    /// occ_P = P(X ≥ occ), X ~ Binomial(T_s, exp_freq) (<c>binomial_boe</c>); occ_E = occ_P · tested dyads;
    /// occ_sig = −log₁₀ occ_E.</para>
    /// <para>Tested dyads: those with occ ≥ <see cref="DyadAnalysisOptions.MinCount"/> and exp_freq &gt; 0. MinCount ≤ 0
    /// reproduces RSAT without <c>-lth occ</c>: every combination of an observed first monad, a spacing, and a second monad
    /// allowed by the type (any: observed monads; dr: M₁; ir: rc(M₁); rep: both) is tested, with occ 0 when unseen.</para>
    /// <para>Deliberate differences from the RSAT 1.78 code (see the algorithm document, with the runs that show them):
    /// with both strands RSAT's <c>SumRCDyads</c> only sums a pair when its lexicographically smaller member was observed,
    /// so pairs seen only as the larger member are lost (19 of 139 pairs on RSAT's t2 test input); here every observed pair
    /// is summed. z-scores are not rounded to 2 decimals (RSAT <c>%7.2f</c>). Background-table monad frequencies are used
    /// for every monad (RSAT only assigns them to monads seen in the input).</para>
    /// Results are ordered by first occurrence (sequence, position), then spacing, then pattern; unseen combinations last.
    /// </remarks>
    /// <param name="sequences">Input DNA sequences (no null elements).</param>
    /// <param name="options">Options (null = RSAT defaults: m = 3, spacing 0–20, any dyad, both strands, -noov, occ ≥ 1, monad background).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Monad length &lt; 1, negative or inverted spacing range, or an unknown enum value.</exception>
    /// <exception cref="ArgumentException">A null sequence, or a background table with words of another length.</exception>
    public static DyadAnalysisReport AnalyzeDyads(IEnumerable<DnaSequence> sequences, DyadAnalysisOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        options ??= new DyadAnalysisOptions();
        options.Validate(nameof(options));
        int m = options.MonadLength, minSp = options.MinSpacing, maxSp = options.MaxSpacing;
        bool both = options.Strands == OligoStrandMode.Both;
        var type = options.Type;
        options.Background.ValidateFor(m, nameof(options));

        List<string> seqs = SequenceStrings(sequences);

        int spacings = maxSp - minSp + 1;
        var valid = new long[spacings];
        var occSum = new long[spacings];
        var ovlSum = new long[spacings];
        var monadOcc = new Dictionary<string, long>(StringComparer.Ordinal);
        long monadTotal = 0;
        var tallies = new Dictionary<string, DyadTally>(StringComparer.Ordinal);
        var rcCache = new Dictionary<string, string>(StringComparer.Ordinal);

        DyadTally Tally(string pattern, string m1, int s, string m2)
        {
            if (!tallies.TryGetValue(pattern, out var t))
            {
                t = new DyadTally(m1, s, m2);
                tallies.Add(pattern, t);
            }
            return t;
        }

        // RSAT CountDyads.
        for (int si = 0; si < seqs.Count; si++)
        {
            string seq = seqs[si];
            for (int s = minSp; s <= maxSp; s++)
                valid[s - minSp] += Math.Max(0, seq.Length - 2 * m - s + 1);

            int lastPos = seq.Length - m;
            var lastCounted = new Dictionary<string, int>(StringComparer.Ordinal);
            for (int pos1 = 0; pos1 <= lastPos; pos1++)
            {
                string m1 = seq.Substring(pos1, m);
                monadOcc[m1] = monadOcc.GetValueOrDefault(m1) + 1;
                monadTotal++;
                for (int s = minSp; s <= maxSp; s++)
                {
                    int pos2 = pos1 + m + s;
                    if (pos2 > lastPos) break;
                    string m2 = seq.Substring(pos2, m);
                    if (!DyadTypeAccepts(type, m1, m2, rcCache)) continue;

                    string pattern = DyadPattern(m1, s, m2);
                    if (!options.CountOverlapping)
                    {
                        if (lastCounted.TryGetValue(pattern, out int last) && pos1 - last < 2 * m + s)
                        {
                            Tally(pattern, m1, s, m2).Overlaps++;
                            ovlSum[s - minSp]++;
                            continue;
                        }
                        lastCounted[pattern] = pos1;
                        if (both)
                            lastCounted[DyadPattern(CachedReverseComplement(m2, rcCache), s, CachedReverseComplement(m1, rcCache))] = pos1;
                    }

                    var t = Tally(pattern, m1, s, m2);
                    t.Occurrences++;
                    t.Positions.Add(new OligoOccurrence(si, pos1));
                    occSum[s - minSp]++;
                }
            }
        }

        // Candidate patterns, grouped by reverse-complement pairs.
        var groups = new Dictionary<string, DyadGroup>(StringComparer.Ordinal);
        void AddCandidate(string m1, int s, string m2)
        {
            string pattern = DyadPattern(m1, s, m2);
            string rm1 = CachedReverseComplement(m2, rcCache), rm2 = CachedReverseComplement(m1, rcCache);
            string rcPattern = both ? DyadPattern(rm1, s, rm2) : pattern;
            bool swap = string.CompareOrdinal(rcPattern, pattern) < 0;
            string key = swap ? rcPattern : pattern;
            if (groups.ContainsKey(key)) return;
            var g = swap ? new DyadGroup(rcPattern, rm1, s, rm2, pattern) : new DyadGroup(pattern, m1, s, m2, rcPattern);
            foreach (string member in g.Partner == g.Pattern ? new[] { g.Pattern } : new[] { g.Pattern, g.Partner })
            {
                if (!tallies.TryGetValue(member, out var t)) continue;
                g.Occurrences += t.Occurrences;
                g.Overlaps += t.Overlaps;
                g.Positions.AddRange(t.Positions);
            }
            g.Positions.Sort(OligoOccurrence.Compare);
            groups.Add(key, g);
        }

        foreach (var t in tallies.Values)
            AddCandidate(t.First, t.Spacing, t.Second);
        if (options.MinCount <= 0)
        {
            var observedMonads = monadOcc.Keys.Order(StringComparer.Ordinal).ToList();
            foreach (string m1 in observedMonads)
            {
                for (int s = minSp; s <= maxSp; s++)
                {
                    foreach (string m2 in SecondElements(type, m1, observedMonads, rcCache))
                        AddCandidate(m1, s, m2);
                }
            }
        }

        // Expected frequencies, P-values.
        var background = options.Background;
        double MonadFrequency(string monad) => background.MonadFrequency(monad, monadOcc, monadTotal);
        double FromMonads(DyadGroup g)
        {
            double e = MonadFrequency(g.First) * MonadFrequency(g.Second);
            if (both)
            {
                string rc1 = CachedReverseComplement(g.First, rcCache), rc2 = CachedReverseComplement(g.Second, rcCache);
                if (rc1 != g.Second)
                    e += MonadFrequency(rc1) * MonadFrequency(rc2);
            }
            return e;
        }

        var drafts = new List<(DyadGroup G, double ExpFreq, bool FromMonads)>();
        foreach (var g in groups.Values)
        {
            if (g.Occurrences < options.MinCount) continue;
            bool fromMonads = true;
            double expFreq;
            if (background.IsDyadTable)
            {
                expFreq = background.DyadFrequency(g.Pattern);
                if (both && g.Partner != g.Pattern)
                    expFreq += background.DyadFrequency(g.Partner);
                fromMonads = false;
                if (!(expFreq > 0))
                {
                    if (maxSp > 20 || g.Occurrences <= 0)
                    {
                        expFreq = FromMonads(g);
                        fromMonads = true;
                    }
                    else
                    {
                        expFreq = double.NaN; // RSAT: "Expected frequency must be > 0", not tested
                    }
                }
            }
            else
            {
                expFreq = FromMonads(g);
            }
            drafts.Add((g, expFreq, fromMonads));
        }

        int tested = drafts.Count(d => d.ExpFreq > 0);
        double logTested = Math.Log(tested);
        var equiprobable = new[] { 0.25, 0.25, 0.25, 0.25 };
        var result = new List<DyadStatistics>(drafts.Count);
        foreach (var (g, expFreq, fromMonads) in drafts)
        {
            int si = g.Spacing - minSp;
            int occ = g.Occurrences;
            bool testable = expFreq > 0;
            double expOcc = testable ? expFreq * occSum[si] : double.NaN;
            long total = occSum[si] + ovlSum[si];
            double obsFreq = occSum[si] == 0 ? 0.0 : (double)occ / total;
            double ratio = !testable || expOcc == 0 ? 0.0 : occ / expOcc;
            double ovlp = OverlapCoefficient(g.First + g.Second, equiprobable);
            double variance = testable ? valid[si] * expFreq * (2 * ovlp - 1 - (4 * m + 1) * expFreq) : double.NaN;
            double z = variance > 0 ? (occ - expOcc) / Math.Sqrt(variance) : double.NaN;
            double logOccP = testable
                ? StatisticsHelper.LogBinomialUpperTail(occ, valid[si], Math.Min(0.0, Math.Log(expFreq)))
                : double.NaN;
            double logOccE = logOccP + logTested;
            result.Add(new DyadStatistics(
                Pattern: g.Pattern,
                FirstMonad: g.First,
                Spacing: g.Spacing,
                SecondMonad: g.Second,
                ReverseComplement: both ? g.Partner : null,
                Occurrences: occ,
                Overlaps: g.Overlaps,
                Positions: g.Positions.AsReadOnly(),
                ObservedFrequency: obsFreq,
                ExpectedFrequency: testable ? expFreq : double.NaN,
                ExpectedOccurrences: expOcc,
                ExpectedVariance: variance,
                OverlapCoefficient: ovlp,
                ZScore: z,
                Ratio: ratio,
                OccurrenceProbability: Math.Exp(logOccP),
                OccurrenceEValue: Math.Exp(logOccE),
                OccurrenceSignificance: -logOccE / Ln10,
                IsDirectRepeat: g.First == g.Second,
                IsReversePalindrome: CachedReverseComplement(g.First, rcCache) == g.Second,
                ExpectedFromMonads: fromMonads));
        }

        result.Sort((x, y) =>
        {
            bool xe = x.Positions.Count == 0, ye = y.Positions.Count == 0;
            if (xe != ye) return xe ? 1 : -1;
            int c = xe ? 0 : OligoOccurrence.Compare(x.Positions[0], y.Positions[0]);
            if (c == 0) c = x.Spacing.CompareTo(y.Spacing);
            return c != 0 ? c : string.CompareOrdinal(x.Pattern, y.Pattern);
        });

        var spacingStats = new List<DyadSpacingStatistics>(spacings);
        for (int s = minSp; s <= maxSp; s++)
            spacingStats.Add(new DyadSpacingStatistics(s, valid[s - minSp], occSum[s - minSp], ovlSum[s - minSp]));

        return new DyadAnalysisReport(
            Dyads: result.AsReadOnly(),
            MonadLength: m,
            MinSpacing: minSp,
            MaxSpacing: maxSp,
            Type: type,
            Strands: options.Strands,
            CountOverlapping: options.CountOverlapping,
            SequenceCount: seqs.Count,
            MonadOccurrences: monadTotal,
            Spacings: spacingStats.AsReadOnly(),
            TestedPatterns: tested,
            PossibleDyads: PossibleDyads(m, minSp, maxSp, type, both));
    }

    /// <summary>RSAT dyad notation: first monad, "n{s}", second monad (e.g. <c>ACGn{4}TTA</c>).</summary>
    internal static string DyadPattern(string first, int spacing, string second)
        => string.Concat(first, "n{", spacing.ToString(CultureInfo.InvariantCulture), "}", second);


    private static bool DyadTypeAccepts(DyadType type, string m1, string m2, Dictionary<string, string> rcCache) => type switch
    {
        DyadType.DirectRepeat => m1 == m2,
        DyadType.InvertedRepeat => m2 == CachedReverseComplement(m1, rcCache),
        DyadType.Repeat => m1 == m2 || m2 == CachedReverseComplement(m1, rcCache),
        _ => true,
    };

    // RSAT SecondElement: the second monads combined with m1 for the dyad type.
    private static IEnumerable<string> SecondElements(DyadType type, string m1, List<string> observedMonads, Dictionary<string, string> rcCache)
        => type switch
        {
            DyadType.DirectRepeat => new[] { m1 },
            DyadType.InvertedRepeat => new[] { CachedReverseComplement(m1, rcCache) },
            DyadType.Repeat => new[] { m1, CachedReverseComplement(m1, rcCache) },
            _ => observedMonads,
        };

    /// <summary>
    /// RSAT <c>CalcPossibleDyads</c>: N = 4^m; dr: N (N/2 with both strands); ir: N; rep: dr + ir; any: N² (both strands:
    /// N² − (N² − N)/2); × (max − min + 1) when max &gt; min.
    /// </summary>
    private static double PossibleDyads(int m, int minSp, int maxSp, DyadType type, bool both)
    {
        double oligos = Math.ScaleB(1.0, 2 * m);
        double dr = both ? oligos / 2 : oligos;
        double ir = oligos;
        double n = type switch
        {
            DyadType.DirectRepeat => dr,
            DyadType.InvertedRepeat => ir,
            DyadType.Repeat => dr + ir,
            _ => both ? oligos * oligos - (oligos * oligos - oligos) / 2 : oligos * oligos,
        };
        if (maxSp > minSp) n *= maxSp - minSp + 1;
        return n;
    }

    private sealed class DyadTally(string first, int spacing, string second)
    {
        public string First { get; } = first;
        public int Spacing { get; } = spacing;
        public string Second { get; } = second;
        public int Occurrences;
        public int Overlaps;
        public readonly List<OligoOccurrence> Positions = new();
    }

    private sealed class DyadGroup(string pattern, string first, int spacing, string second, string partner)
    {
        public string Pattern { get; } = pattern;
        public string First { get; } = first;
        public int Spacing { get; } = spacing;
        public string Second { get; } = second;
        public string Partner { get; } = partner;
        public int Occurrences;
        public int Overlaps;
        public readonly List<OligoOccurrence> Positions = new();
    }

    #endregion
}

/// <summary>Dyad types of RSAT <c>dyad-analysis -type</c>.</summary>
public enum DyadType
{
    /// <summary>Any pair of monads (<c>any</c>, default).</summary>
    Any,

    /// <summary>Direct repeats, M₂ = M₁ (<c>dr</c>).</summary>
    DirectRepeat,

    /// <summary>Inverted repeats, M₂ = reverse complement of M₁ (<c>ir</c>).</summary>
    InvertedRepeat,

    /// <summary>Direct or inverted repeats (<c>rep</c>).</summary>
    Repeat,
}

/// <summary>Background for RSAT <c>dyad-analysis</c> expected dyad frequencies.</summary>
public sealed class DyadBackgroundModel
{
    private readonly Dictionary<string, double>? _monads;
    private readonly Dictionary<string, double>? _dyads;
    private readonly int _wordLength;

    private DyadBackgroundModel(Dictionary<string, double>? monads, Dictionary<string, double>? dyads, int wordLength)
    {
        _monads = monads;
        _dyads = dyads;
        _wordLength = wordLength;
    }

    /// <summary>Monad frequencies observed in the input (RSAT default, <c>-bg monads</c>).</summary>
    public static DyadBackgroundModel Monads { get; } = new(null, null, 0);

    /// <summary>
    /// Expected monad frequencies from a table (RSAT <c>-mncf</c> / organism monad files): exp_freq(D) is computed from these
    /// instead of the input monad frequencies. Missing monads have frequency 0.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="frequencies"/> is null.</exception>
    /// <exception cref="ArgumentException">Empty table, non-ACGT or mixed-length words, or a negative / non-finite value.</exception>
    public static DyadBackgroundModel MonadFrequencies(IReadOnlyDictionary<string, double> frequencies)
    {
        ArgumentNullException.ThrowIfNull(frequencies);
        var table = new Dictionary<string, double>(StringComparer.Ordinal);
        int length = -1;
        foreach (var (raw, value) in frequencies)
        {
            string word = (raw ?? string.Empty).ToUpperInvariant();
            if (word.Length == 0 || word.Any(c => MotifFinder.AcgtIndex(c) < 0))
                throw new ArgumentException($"Monad '{raw}' is not an A/C/G/T word.", nameof(frequencies));
            if (length < 0) length = word.Length;
            else if (word.Length != length)
                throw new ArgumentException("All monads must have the same length.", nameof(frequencies));
            if (!double.IsFinite(value) || value < 0)
                throw new ArgumentException($"Invalid frequency {value} for '{raw}'.", nameof(frequencies));
            table[word] = value;
        }
        if (table.Count == 0)
            throw new ArgumentException("The monad frequency table is empty.", nameof(frequencies));
        return new DyadBackgroundModel(table, null, length);
    }

    /// <summary>
    /// Expected single-strand dyad frequencies (RSAT <c>-expfreq</c> / organism dyad files), keyed in RSAT notation
    /// <c>acgn{4}tta</c> (case-insensitive). As RSAT <c>ReadExpectedFrequencies</c>, a dyad whose reverse complement is
    /// missing passes its value to it; with both strands exp_freq(D|D') = exp_freq(D) + exp_freq(D').
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="frequencies"/> is null.</exception>
    /// <exception cref="ArgumentException">Empty table, a key not in <c>MONADn{s}MONAD</c> form, mixed monad lengths, or a negative / non-finite value.</exception>
    public static DyadBackgroundModel DyadFrequencies(IReadOnlyDictionary<string, double> frequencies)
    {
        ArgumentNullException.ThrowIfNull(frequencies);
        var table = new Dictionary<string, double>(StringComparer.Ordinal);
        int length = -1;
        foreach (var (raw, value) in frequencies)
        {
            var (m1, s, m2) = ParseDyad(raw, nameof(frequencies));
            if (length < 0) length = m1.Length;
            else if (m1.Length != length)
                throw new ArgumentException("All dyads must have the same monad length.", nameof(frequencies));
            if (!double.IsFinite(value) || value < 0)
                throw new ArgumentException($"Invalid frequency {value} for '{raw}'.", nameof(frequencies));
            table[MotifFinder.DyadPattern(m1, s, m2)] = value;
        }
        if (table.Count == 0)
            throw new ArgumentException("The dyad frequency table is empty.", nameof(frequencies));
        foreach (string key in table.Keys.Order(StringComparer.Ordinal).ToList())
        {
            var (m1, s, m2) = ParseDyad(key, nameof(frequencies));
            string rc = MotifFinder.DyadPattern(DnaSequence.GetReverseComplementString(m2), s, DnaSequence.GetReverseComplementString(m1));
            if (!table.ContainsKey(rc) && table[key] != 0)
                table[rc] = table[key];
        }
        return new DyadBackgroundModel(null, table, length);
    }

    internal bool IsDyadTable => _dyads is not null;

    internal void ValidateFor(int monadLength, string paramName)
    {
        if (_wordLength > 0 && _wordLength != monadLength)
            throw new ArgumentException(
                $"The background table uses monads of length {_wordLength}, the analysis uses {monadLength}.", paramName);
    }

    internal double MonadFrequency(string monad, Dictionary<string, long> observed, long total)
    {
        if (_monads is not null)
            return _monads.GetValueOrDefault(monad);
        return observed.TryGetValue(monad, out long c) ? (double)c / total : 0.0;
    }

    internal double DyadFrequency(string pattern) => _dyads!.GetValueOrDefault(pattern);

    private static (string First, int Spacing, string Second) ParseDyad(string raw, string paramName)
    {
        string text = raw ?? string.Empty;
        int open = text.IndexOf("n{", StringComparison.OrdinalIgnoreCase);
        int close = open < 0 ? -1 : text.IndexOf('}', open);
        if (open <= 0 || close < 0
            || !int.TryParse(text.AsSpan(open + 2, close - open - 2), NumberStyles.None, CultureInfo.InvariantCulture, out int s))
            throw new ArgumentException($"Dyad '{raw}' is not in MONADn{{s}}MONAD notation.", paramName);
        string m1 = text[..open].ToUpperInvariant(), m2 = text[(close + 1)..].ToUpperInvariant();
        if (m1.Length != m2.Length || m1.Any(c => MotifFinder.AcgtIndex(c) < 0) || m2.Any(c => MotifFinder.AcgtIndex(c) < 0))
            throw new ArgumentException($"Dyad '{raw}' must join two A/C/G/T monads of equal length.", paramName);
        return (m1, s, m2);
    }
}

/// <summary>Options of <see cref="MotifFinder.AnalyzeDyads"/> (RSAT <c>dyad-analysis</c> command-line options; defaults = RSAT's).</summary>
public sealed record DyadAnalysisOptions
{
    /// <summary>Monad length m (RSAT <c>-l</c>, default 3).</summary>
    public int MonadLength { get; init; } = 3;

    /// <summary>Minimal spacing (RSAT <c>-sp min-max</c>, default 0).</summary>
    public int MinSpacing { get; init; }

    /// <summary>Maximal spacing (default 20).</summary>
    public int MaxSpacing { get; init; } = 20;

    /// <summary>Dyad type (RSAT <c>-type</c>, default any).</summary>
    public DyadType Type { get; init; } = DyadType.Any;

    /// <summary>Both strands (RSAT default <c>-2str</c>) or a single strand (<c>-1str</c>).</summary>
    public OligoStrandMode Strands { get; init; } = OligoStrandMode.Both;

    /// <summary>false = RSAT default <c>-noov</c>; true = <c>-ovlp</c>.</summary>
    public bool CountOverlapping { get; init; }

    /// <summary>Occurrence threshold (RSAT <c>-lth occ</c>, default 1); ≤ 0 tests every monad combination (RSAT without a threshold).</summary>
    public int MinCount { get; init; } = 1;

    /// <summary>Expected-frequency model (default: input monad frequencies).</summary>
    public DyadBackgroundModel Background { get; init; } = DyadBackgroundModel.Monads;

    internal void Validate(string paramName)
    {
        if (Background is null)
            throw new ArgumentException("A background model is required.", paramName);
        if (MonadLength < 1)
            throw new ArgumentOutOfRangeException(paramName, MonadLength, "The monad length must be >= 1.");
        if (MinSpacing < 0 || MaxSpacing < MinSpacing)
            throw new ArgumentOutOfRangeException(paramName, $"{MinSpacing}-{MaxSpacing}", "Spacing range must satisfy 0 <= min <= max.");
        if (!Enum.IsDefined(Type))
            throw new ArgumentOutOfRangeException(paramName, Type, "Unknown dyad type.");
        if (!Enum.IsDefined(Strands))
            throw new ArgumentOutOfRangeException(paramName, Strands, "Unknown strand mode.");
    }
}

/// <summary>A dyad (or reverse-complement pair) scored by <see cref="MotifFinder.AnalyzeDyads"/>; NaN = RSAT "NA".</summary>
/// <param name="Pattern">RSAT notation <c>M1n{s}M2</c> (smaller pair member with both strands).</param>
/// <param name="FirstMonad">M₁.</param>
/// <param name="Spacing">s.</param>
/// <param name="SecondMonad">M₂.</param>
/// <param name="ReverseComplement">Pair partner with both strands (equal to Pattern for a reverse palindrome); null for one strand.</param>
/// <param name="Occurrences">occ (pair sum with both strands).</param>
/// <param name="Overlaps">Discarded overlapping occurrences (<c>-noov</c>; RSAT <c>ovl_occ</c>).</param>
/// <param name="Positions">Counted occurrences (sequence, start of M₁), ascending.</param>
/// <param name="ObservedFrequency">obs_freq = occ / (occ_sum_s + ovl_sum_s).</param>
/// <param name="ExpectedFrequency">exp_freq (NaN when not tested).</param>
/// <param name="ExpectedOccurrences">exp_occ = exp_freq · occ_sum_s.</param>
/// <param name="ExpectedVariance">RSAT <c>occ_var</c>.</param>
/// <param name="OverlapCoefficient">RSAT <c>ov_coef</c> of M₁M₂.</param>
/// <param name="ZScore">z-score (NaN when the variance is ≤ 0).</param>
/// <param name="Ratio">occ / exp_occ.</param>
/// <param name="OccurrenceProbability">occ_P.</param>
/// <param name="OccurrenceEValue">occ_E.</param>
/// <param name="OccurrenceSignificance">occ_sig.</param>
/// <param name="IsDirectRepeat">M₁ = M₂ (RSAT remark <c>dir_rep</c>).</param>
/// <param name="IsReversePalindrome">The dyad equals its reverse complement (RSAT remark <c>rc_pal</c>).</param>
/// <param name="ExpectedFromMonads">exp_freq from monad frequencies (always true for the monad backgrounds).</param>
public sealed record DyadStatistics(
    string Pattern,
    string FirstMonad,
    int Spacing,
    string SecondMonad,
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
    bool IsDirectRepeat,
    bool IsReversePalindrome,
    bool ExpectedFromMonads);

/// <summary>Per-spacing counts of <see cref="MotifFinder.AnalyzeDyads"/>.</summary>
/// <param name="Spacing">s.</param>
/// <param name="PossiblePositions">T_s = Σ max(0, L − 2m − s + 1), the binomial trials.</param>
/// <param name="Occurrences">Counted single-strand occurrences (occ_sum_s).</param>
/// <param name="Overlaps">Discarded single-strand overlaps (ovl_sum_s).</param>
public sealed record DyadSpacingStatistics(int Spacing, long PossiblePositions, long Occurrences, long Overlaps);

/// <summary>Result of <see cref="MotifFinder.AnalyzeDyads"/>.</summary>
/// <param name="Dyads">Dyads with occ ≥ the threshold.</param>
/// <param name="MonadLength">m.</param>
/// <param name="MinSpacing">Minimal spacing.</param>
/// <param name="MaxSpacing">Maximal spacing.</param>
/// <param name="Type">Dyad type.</param>
/// <param name="Strands">Strand mode.</param>
/// <param name="CountOverlapping">-ovlp (true) or -noov (false).</param>
/// <param name="SequenceCount">Number of input sequences.</param>
/// <param name="MonadOccurrences">Monad windows counted (RSAT <c>sum_oligo_count</c>).</param>
/// <param name="Spacings">Per-spacing positions and counts.</param>
/// <param name="TestedPatterns">Dyads with a P-value (occ_E multiplier).</param>
/// <param name="PossibleDyads">RSAT <c>nb_possible_dyads</c>.</param>
public sealed record DyadAnalysisReport(
    IReadOnlyList<DyadStatistics> Dyads,
    int MonadLength,
    int MinSpacing,
    int MaxSpacing,
    DyadType Type,
    OligoStrandMode Strands,
    bool CountOverlapping,
    int SequenceCount,
    long MonadOccurrences,
    IReadOnlyList<DyadSpacingStatistics> Spacings,
    int TestedPatterns,
    double PossibleDyads);
