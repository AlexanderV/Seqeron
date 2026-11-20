using Seqeron.Genomics.Infrastructure;
namespace Seqeron.Genomics.Analysis;

// B09 A3-8 (GENOMIC-TANDEM-001): the two EMBOSS 6.6.0 tandem-repeat programs, ported statement by statement from etandem.c / equicktandem.c
// (raw.githubusercontent.com/kimrutherford/EMBOSS, emboss/etandem.c and emboss/equicktandem.c) and checked on the installed binaries
// (scripts/analysis/generate_emboss_tandem_oracle.py). They are diagonal-score models, not exact-run finders: GenomicAnalyzer.FindTandemRepeats stays the exact definition.

/// <summary>A tandem repeat reported by EMBOSS <c>etandem</c> or <c>equicktandem</c>.</summary>
/// <param name="Position">0-based start (EMBOSS <c>Start</c> − 1).</param>
/// <param name="Length">Length in bases (EMBOSS <c>End</c> − <c>Start</c> + 1).</param>
/// <param name="Score">The program's score (<c>etandem</c>: the best running consensus score; <c>equicktandem</c>: the best diagonal score).</param>
/// <param name="RepeatSize">Repeat unit size (<c>Size</c> column).</param>
/// <param name="RepeatCount">Number of copies (<c>Count</c> column; integer division as in the programs).</param>
/// <param name="Identity"><c>etandem</c> only: percent identity of the repeat with its consensus, as the program's 32-bit float, otherwise NaN.</param>
/// <param name="Consensus"><c>etandem</c> only: the consensus unit (lower case, <c>acgtn</c>), phased as printed by the program; otherwise null.</param>
public readonly record struct EmbossTandemRepeat(int Position, int Length, int Score, int RepeatSize, int RepeatCount, float Identity, string? Consensus)
{
    /// <summary>EMBOSS 1-based start.</summary>
    public int Start => Position + 1;

    /// <summary>EMBOSS 1-based inclusive end.</summary>
    public int End => Position + Length;
}

public static partial class GenomicAnalyzer
{
    #region EMBOSS etandem / equicktandem

    private const int EmbossTandemMaxRepeatCells = 200_000_000;

    /// <summary>
    /// EMBOSS <c>etandem</c> (Sanger Centre <c>tandem</c> by J. Thierry-Mieg and R. Durbin, EMBOSS 6.6.0): tandem repeats of every unit size
    /// <paramref name="minRepeat"/>..<paramref name="maxRepeat"/> found by a running consensus score, with the program's reporting rules.
    /// </summary>
    /// <remarks>
    /// <para>
    /// For a unit size n, a candidate starts at every position and slides in steps of n: each step adds one n-mer to a per-column base-count table; a column scores +1 when the
    /// new base is (or becomes) the column's most frequent base and −1 otherwise (ties move the consensus to the new base and score −1); an <c>N</c> (or any non-ACGT symbol)
    /// counts as a fifth base, or as a plain −1 with <paramref name="mismatch"/>. A candidate starts at −n, remembers its best score and where it was reached, and dies when its
    /// score goes below 0 (a candidate that is not above its best score and below 0). Candidates of one size with identical consensus tables are collapsed (the one with the
    /// lower best score goes). Reported candidates need a best score of at least <paramref name="threshold"/> and, unless <paramref name="uniform"/>, a consensus with at least two
    /// different bases; the final list takes the highest score first (ties: the smaller unit size, then the later-found one) and drops every overlapping candidate.
    /// </para>
    /// <para>
    /// <c>Identity = 100 · (score + n·(copies + 1)) / (2·n·copies)</c> in 32-bit float arithmetic, <c>copies = (end − start + 1) / n</c> (integer division). The consensus is
    /// printed from the column after the last base's phase and wraps around, in lower case. Verified against etandem 6.6.0 on 1 600 random and planted sequences
    /// (docs/Evidence/GENOMIC-TANDEM-001-Evidence.md). Hits come in the program's report order (highest score first), not position order.
    /// </para>
    /// </remarks>
    /// <param name="sequence">Nucleotide string; A/C/G/T in either case, everything else counts as N.</param>
    /// <param name="minRepeat">Smallest unit size (<c>-minrepeat</c>, ≥ 2).</param>
    /// <param name="maxRepeat">Largest unit size (<c>-maxrepeat</c>, ≥ <paramref name="minRepeat"/>).</param>
    /// <param name="threshold">Minimum score (<c>-threshold</c>, default 20).</param>
    /// <param name="mismatch">Count N as a mismatch (<c>-mismatch</c>) instead of a fifth base.</param>
    /// <param name="uniform">Also report single-base consensus repeats (<c>-uniform</c>).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minRepeat"/> &lt; 2 or <paramref name="maxRepeat"/> &lt; <paramref name="minRepeat"/>, or the
    /// work (sequence length × unit sizes) is beyond the supported limit.</exception>
    /// <exception cref="InvalidOperationException"><paramref name="threshold"/> is <see cref="int.MinValue"/> and a repeat qualifies: the program divides by zero (SIGFPE) there (see the remarks).</exception>
    public static IReadOnlyList<EmbossTandemRepeat> FindTandemRepeatsEtandem(
        string sequence, int minRepeat, int maxRepeat, int threshold = 20, bool mismatch = false, bool uniform = false, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minRepeat < 2) throw new ArgumentOutOfRangeException(nameof(minRepeat), minRepeat, "minRepeat must be 2 or higher (etandem -minrepeat).");
        if (maxRepeat < minRepeat) throw new ArgumentOutOfRangeException(nameof(maxRepeat), maxRepeat, "maxRepeat must be at least minRepeat (etandem -maxrepeat).");
        sequence = EmbossSequenceReader.Clean(sequence);       // the reader drops the gap characters and every other non-letter but * (A4-3)
        int nbase = sequence.Length;
        int effectiveMax = Math.Min(maxRepeat, nbase);
        if ((long)nbase * Math.Max(1, effectiveMax - minRepeat + 1) > EmbossTandemMaxRepeatCells)
            throw new ArgumentOutOfRangeException(nameof(maxRepeat), maxRepeat, $"sequence length x number of unit sizes exceeds {EmbossTandemMaxRepeatCells}.");

        // etandem: ajSeqcvtNewC("ACGTN") -> 1..5, anything else 0; base = code - 1, negative -> 4 (N).
        var codes = new int[nbase];
        for (int i = 0; i < nbase; i++) codes[i] = EtandemBaseCode(sequence[i]);

        var reportRoot = new EtandemCons(0) { BestScore = unchecked(threshold - 1) };       // etandem.c: reportRoot->bestScore = thresh - 1 (a 32-bit wrap for int.MinValue)
        var ring = new int[Math.Max(nbase, 1)];

        for (int n = minRepeat; n <= maxRepeat; n++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (n > nbase) break; // ibase < n - 1 for every position: no candidate ever starts (same for every larger size)
            var root = new EtandemCons(0);
            for (int ibase = 0; ibase < nbase; ibase++)
            {
                if ((ibase & 1023) == 0) cancellationToken.ThrowIfCancellationRequested();
                int phase = ibase % n;
                ring[phase] = codes[ibase];
                if (ibase < n - 1) continue;

                root.Next = new EtandemCons(n)
                {
                    Start = ibase - n + 1,
                    Score = -n,
                    BestScore = -n,
                    Phase = phase, // phase of the last base of the pattern
                    Repeat = n,
                    Next = root.Next,
                };

                // add the last n-mer to the active candidates
                var olda = root;
                var a = olda.Next;
                while (a is not null)
                {
                    if (a.Phase == phase)
                    {
                        EtandemAddNmer(a, ring, n, mismatch);
                        if (a.Score > a.BestScore)
                        {
                            a.BestScore = a.Score;
                            a.IBest = ibase;
                            Array.Copy(a.Max, a.BestMax, n);
                        }
                        else if (a.Score < 0)
                        {
                            olda.Next = a.Next;
                            EtandemReport(a, threshold, uniform, reportRoot);
                            a = olda;
                        }
                    }
                    olda = a;
                    a = a.Next;
                }

                // remove duplicate max tables
                olda = root;
                a = olda.Next;
                while (a is not null)
                {
                    if (a.Phase == phase) a = EtandemDeduplicate(olda, a, n);
                    olda = a;
                    a = a.Next;
                }
            }

            EtandemCons? c;
            while ((c = root.Next) is not null)
            {
                root.Next = c.Next;
                EtandemReport(c, threshold, uniform, reportRoot);
            }
        }

        // etandem_finalReport: highest score first (ties: smaller unit size), drop every overlapping entry including the reported one.
        // The program starts its search for the top entry at the report root, whose best score is threshold - 1; for threshold = int.MinValue that wraps to int.MaxValue, nothing beats it, the root is reported
        // as the top entry and its repeat size 0 divides (measured on 6.6.0: "Floating point exception", SIGFPE, exit status 136, whenever one repeat qualifies; a sequence without a qualifying repeat prints nothing).
        if (reportRoot.Next is not null && reportRoot.BestScore == int.MaxValue)
            throw new InvalidOperationException("etandem -threshold -2147483648 wraps the report root's best score threshold - 1 to 2147483647: the program reports its empty root and divides by the repeat size 0 (SIGFPE, \"Floating point exception\") as soon as one repeat qualifies; use a threshold above -2147483648.");
        var result = new List<EmbossTandemRepeat>();
        while (reportRoot.Next is not null)
        {
            var top = reportRoot;
            for (var a = reportRoot.Next; a is not null; a = a.Next)
                if (a.BestScore > top.BestScore || (a.BestScore == top.BestScore && a.Repeat < top.Repeat)) top = a;

            result.Add(EtandemBasicReport(top));
            int start = top.Start, end = top.IBest;
            var olda = reportRoot;
            var cur = olda.Next;
            while (cur is not null)
            {
                if (cur.IBest >= start && cur.Start <= end)
                {
                    olda.Next = cur.Next;
                    cur = olda;
                }
                olda = cur;
                cur = cur.Next;
            }
        }
        return result;
    }

    /// <summary>
    /// EMBOSS <c>equicktandem</c> (the quick pre-screen of the same program family, EMBOSS 6.6.0): for every unit size (gap) 1..<paramref name="maxRepeat"/> the sequence is
    /// compared with itself shifted by that gap; a +1 / -1 diagonal score runs along the shift and a stretch whose best score reaches <paramref name="threshold"/> is reported
    /// and masked out of every later (larger) gap.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Exactly the program's scan: <c>back</c> / <c>front = back + gap</c> walk the numeric sequence; equal symbols raise the score (a new run starts at score 1 after the score has
    /// fallen to 0 or below), different symbols lower it by 1, a masked symbol restarts the scan after it; the best score and its end are remembered; when a run restarts, or at a
    /// masked symbol, or at the end, a best score at least the threshold is reported as <c>Position = (start of the best run)</c>, <c>Length</c> up to the best end, and the scan resumes at the
    /// best end. <c>Count = (length) / gap</c>. Symbols are A, C, G, T (case-insensitive), N, other (all other characters share one class); the end of the sequence acts as one
    /// more "other" symbol, as in the C code (<c>front - sq &lt;= len</c> reads the terminator). Results are in report order (by gap, then position). No consensus and no
    /// identity are computed (<see cref="EmbossTandemRepeat.Consensus"/> null, <see cref="EmbossTandemRepeat.Identity"/> NaN).
    /// </para>
    /// <para>Verified against equicktandem 6.6.0 on 1 300 random and planted sequences (docs/Evidence/GENOMIC-TANDEM-001-Evidence.md).</para>
    /// </remarks>
    /// <param name="sequence">Nucleotide string.</param>
    /// <param name="maxRepeat">Largest unit size (<c>-maxrepeat</c>, at least 1; limited to the sequence length as in the program).</param>
    /// <param name="threshold">Minimum score (<c>-threshold</c>, at least 1: the program reads uninitialised memory for a threshold of 0 or less).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxRepeat"/> &lt; 1 or <paramref name="threshold"/> &lt; 1.</exception>
    public static IReadOnlyList<EmbossTandemRepeat> FindTandemRepeatsEquicktandem(
        string sequence, int maxRepeat, int threshold = 20, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (maxRepeat < 1) throw new ArgumentOutOfRangeException(nameof(maxRepeat), maxRepeat, "maxRepeat must be at least 1 (equicktandem -maxrepeat).");
        if (threshold < 1) throw new ArgumentOutOfRangeException(nameof(threshold), threshold, "threshold must be at least 1 (equicktandem reads uninitialised memory below that).");
        sequence = EmbossSequenceReader.Clean(sequence);       // the reader drops the gap characters and every other non-letter but * (A4-3)
        int len = sequence.Length;
        if (maxRepeat > len) maxRepeat = len;
        if ((long)len * maxRepeat > EmbossTandemMaxRepeatCells)
            throw new ArgumentOutOfRangeException(nameof(maxRepeat), maxRepeat, $"sequence length x unit sizes exceeds {EmbossTandemMaxRepeatCells}.");

        // ajSeqcvtNewNumberC("ACGTN"): 1-4 = ACGT, 5 = N, 0 = other; the terminator (0) follows the last base; 'Z' (90) marks reported repeats.
        var sq = new byte[len + 1];
        for (int i = 0; i < len; i++)
        {
            sq[i] = sequence[i] switch
            {
                'A' or 'a' => 1,
                'C' or 'c' => 2,
                'G' or 'g' => 3,
                'T' or 't' => 4,
                'N' or 'n' => 5,
                _ => 0,
            };
        }
        const byte Z = (byte)'Z';

        var result = new List<EmbossTandemRepeat>();
        int maxback = 0, maxfront = 0;

        void Report(int gap, int best)
        {
            result.Add(new EmbossTandemRepeat(maxback, maxfront - maxback + 1, best, gap, (maxfront - maxback + 1) / gap, float.NaN, null));
            for (int cp = maxback; cp <= maxfront; cp++) sq[cp] = Z;
        }

        for (int gap = 1; gap <= maxRepeat; gap++)
        {
            cancellationToken.ThrowIfCancellationRequested();
            int back = 0, front = gap, score = 0, max = 0;
            while (front <= len)
            {
                if (sq[front] == Z)
                {
                    if (max >= threshold)
                    {
                        Report(gap, max);
                        back = maxfront;
                    }
                    else back = front;
                    front = back + gap;
                    score = max = 0;
                }
                else if (sq[front] != sq[back]) score--;
                else if (score <= 0)
                {
                    if (max >= threshold)
                    {
                        Report(gap, max);
                        back = maxfront;
                        front = back + gap;
                        score = max = 0;
                    }
                    else
                    {
                        maxback = back;
                        score = 1;
                    }
                }
                else if (++score > max)
                {
                    max = score;
                    maxfront = front;
                }
                back++;
                front++;
            }

            if (max >= threshold) Report(gap, max);
        }
        return result;
    }

    /// <summary>Adds the ring's n-mer to a candidate's column tables (etandem main loop).</summary>
    private static void EtandemAddNmer(EtandemCons a, int[] ring, int n, bool mismatch)
    {
        for (int i = 0; i < n; i++)
        {
            int x = ring[i];
            if (x == 4 && mismatch) { a.Score--; continue; }
            a.Tab[x + 5 * i]++;
            if (x == a.Max[i]) a.Score++;
            else if (a.Tab[x + 5 * i] > a.Tab[a.Max[i] + 5 * i]) { a.Max[i] = x; a.Score++; }
            else
            {
                a.Score--;
                if (a.Tab[x + 5 * i] == a.Tab[a.Max[i] + 5 * i]) a.Max[i] = x;
            }
        }
    }

    /// <summary>
    /// "Remove duplicate max tables" for the candidate <paramref name="a"/> (preceded by <paramref name="olda"/>): against every later candidate with the same consensus table
    /// the one with the lower best score goes (ties: <paramref name="a"/>). Returns the node the outer scan continues from (<paramref name="olda"/> when <paramref name="a"/> was removed).
    /// </summary>
    private static EtandemCons EtandemDeduplicate(EtandemCons olda, EtandemCons a, int n)
    {
        var oldb = a;
        var b = a.Next;
        while (b is not null)
        {
            bool same = true;
            for (int j = 0; j < n && same; j++) same = a.Max[j] == b.Max[j];
            if (same)
            {
                if (a.BestScore > b.BestScore)
                {
                    oldb.Next = b.Next; // remove b
                    b = oldb;
                }
                else
                {
                    olda.Next = a.Next; // remove a
                    return olda;
                }
            }
            oldb = b;
            b = b.Next;
        }
        return a;
    }

    private static int EtandemBaseCode(char c) => c switch
    {
        'A' or 'a' => 0,
        'C' or 'c' => 1,
        'G' or 'g' => 2,
        'T' or 't' => 3,
        _ => 4,
    };

    private static void EtandemReport(EtandemCons a, int threshold, bool uniform, EtandemCons reportRoot)
    {
        if (a.BestScore >= threshold)
        {
            bool good = uniform;
            if (!good)
            {
                int first = a.BestMax[0];
                for (int j = 1; j < a.Repeat; j++)
                    if (a.BestMax[j] != first) { good = true; break; }
            }
            if (good)
            {
                a.Next = reportRoot.Next;
                reportRoot.Next = a;
            }
        }
    }

    private static EmbossTandemRepeat EtandemBasicReport(EtandemCons a)
    {
        int n = a.Repeat;
        int copies = (a.IBest - a.Start + 1) / n;
        float perc = 100.0f * (a.BestScore + n * (copies + 1)) / (2.0f * n * copies);
        const string letters = "acgtn";
        var consensus = new System.Text.StringBuilder(n);
        for (int j = (a.Phase + 1) % n; j < n; j++) consensus.Append(letters[a.BestMax[j]]);
        if ((a.Phase + 1) % n != 0)
            for (int j = 0; j <= a.Phase; j++) consensus.Append(letters[a.BestMax[j]]);
        return new EmbossTandemRepeat(a.Start, a.IBest - a.Start + 1, a.BestScore, n, copies, perc, consensus.ToString());
    }

    private sealed class EtandemCons
    {
        public EtandemCons(int n)
        {
            Tab = new int[5 * n + 5];
            Max = new int[n + 1];
            BestMax = new int[n + 1];
        }

        public int[] Tab { get; }
        public int[] Max { get; }
        public int[] BestMax { get; }
        public int Start, Score, BestScore, IBest, Phase, Repeat;
        public EtandemCons? Next;
    }

    #endregion
}
