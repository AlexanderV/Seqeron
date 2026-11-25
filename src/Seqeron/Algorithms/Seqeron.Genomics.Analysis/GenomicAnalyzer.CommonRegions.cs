using Seqeron.Genomics.Alignment;
namespace Seqeron.Genomics.Analysis;

// B09 A2-30 / A2-35: strand, mode, k-mismatch and multi-sequence variants of the common-region search, and the string overloads.
// The exact engines are the canonical SuffixTree ones (matching statistics, FindMaximalExactMatches, FindMaximalUniqueMatches,
// FindAllLongestCommonSubstrings, FindLongestCommonSubstrings); nothing here re-implements them.
public static partial class GenomicAnalyzer
{
    // A non-matching symbol of a query. It is outside every text's separator pool (0x100..0xD7FF) and outside A-Z, so it equals
    // nothing in the tree and no match can contain it.
    public const char QueryBreak = '';

    // Upper-case query with every symbol outside A/C/G/T replaced by QueryBreak.
    public static string QueryText(string upper)
    {
        var chars = upper.ToCharArray();
        for (int i = 0; i < chars.Length; i++)
            if (chars[i] is not ('A' or 'C' or 'G' or 'T'))
                chars[i] = QueryBreak;
        return new string(chars);
    }

    #region Longest common region on one or both strands

    /// <summary>
    /// <see cref="FindLongestCommonRegion(DnaSequence, DnaSequence)"/> on the selected strands: the longest substring shared by the first
    /// sequence and the second sequence (forward) and/or the reverse complement of the second sequence (<c>mummer -r</c> / <c>-b</c>).
    /// </summary>
    /// <remarks>
    /// With <see cref="StrandSelection.ReverseComplement"/> the region is a substring of the first sequence whose reverse complement occurs in the
    /// second; <see cref="CommonRegion.PositionInSecond"/> is then the 0-based forward-strand start of that window and
    /// <see cref="CommonRegion.IsReverseComplement"/> is true. With <see cref="StrandSelection.Both"/> the longer of the two strands wins, a tie goes to the
    /// forward strand (same tie rule inside a strand as the 2-argument method).
    /// </remarks>
    /// <param name="sequence1">First (indexed) sequence.</param>
    /// <param name="sequence2">Second (query) sequence.</param>
    /// <param name="strands">Strands to search (default: forward only, the 2-argument behaviour).</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    public static CommonRegion FindLongestCommonRegion(DnaSequence sequence1, DnaSequence sequence2, StrandSelection strands)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return LongestCommonRegionCore(sequence1.SuffixTree, sequence2.Sequence, strands);
    }

    /// <summary>
    /// <see cref="FindLongestCommonRegion(DnaSequence, DnaSequence, StrandSelection)"/> on nucleotide <b>strings</b> (N runs, IUPAC codes, lower case):
    /// only A/C/G/T match, any other symbol matches nothing and ends a region.
    /// </summary>
    /// <param name="sequence1">First (indexed) sequence.</param>
    /// <param name="sequence2">Second (query) sequence.</param>
    /// <param name="strands">Strands to search (default forward).</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    public static CommonRegion FindLongestCommonRegion(string sequence1, string sequence2, StrandSelection strands = StrandSelection.Forward)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        if (sequence1.Length == 0 || sequence2.Length == 0)
            return CommonRegion.None;
        return LongestCommonRegionCore(BuildMaskedTree(sequence1, BuildLetterMask(NucleotideAlphabet)), sequence2.ToUpperInvariant(), strands);
    }

    private static CommonRegion LongestCommonRegionCore(global::SuffixTree.SuffixTree tree, string queryUpper, StrandSelection strands)
    {
        var best = CommonRegion.None;
        string forwardQuery = QueryText(queryUpper);
        if (strands.HasFlag(StrandSelection.Forward))
            best = LongestOnStrand(tree, forwardQuery, reverse: false);
        if (strands.HasFlag(StrandSelection.ReverseComplement))
        {
            var reverse = LongestOnStrand(tree, QueryText(ReverseComplementAcgt(queryUpper)), reverse: true);
            if (reverse.Length > best.Length)
                best = reverse;
        }

        return best;

        static CommonRegion LongestOnStrand(global::SuffixTree.SuffixTree tree, string query, bool reverse)
        {
            var (lcs, positionsInText, positionsInQuery) = tree.FindAllLongestCommonSubstrings(query);
            if (string.IsNullOrEmpty(lcs))
                return CommonRegion.None;
            int start = positionsInQuery[0];
            return new CommonRegion(lcs, positionsInText[0], reverse ? query.Length - start - lcs.Length : start, reverse);
        }
    }

    #endregion

    #region Common regions by mode and strand

    /// <summary>
    /// Common regions of two sequences under an explicit <see cref="CommonRegionMode"/> and strand selection — the MUMmer 3 modes
    /// (<c>-maxmatch</c>, <c>-mumreference</c>, <c>-mum</c>) with <c>-r</c> / <c>-b</c>, plus the right-maximal set of
    /// <see cref="FindCommonRegions(DnaSequence, DnaSequence, int)"/>.
    /// </summary>
    /// <remarks>
    /// Regions are listed forward strand first (ordered by their start in the second sequence), then the reverse-complement ones in the order of the
    /// reverse-complemented query (what <c>mummer</c> prints); a reverse region has <see cref="CommonRegion.IsReverseComplement"/> = true and
    /// <see cref="CommonRegion.PositionInSecond"/> = the forward-strand start of the window that spells the reverse complement of
    /// <see cref="CommonRegion.Sequence"/>. The exact engines are the canonical suffix-tree ones (see <see cref="CommonRegionMode"/>); uniqueness on the reverse strand
    /// is measured against the reverse-complemented query, as in MUMmer. Only A/C/G/T match. Verified against <c>mummer -maxmatch / -mumreference / -mum -b -l L</c>
    /// (3.23): see docs/Evidence/GENOMIC-COMMON-001-Evidence.md.
    /// </remarks>
    /// <param name="sequence1">Reference (indexed) sequence.</param>
    /// <param name="sequence2">Query sequence.</param>
    /// <param name="minLength">Minimum region length (<c>mummer -l</c>); values below 1 are treated as 1.</param>
    /// <param name="mode">Which common regions to report.</param>
    /// <param name="strands">Strands to search (default forward).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    public static IReadOnlyList<CommonRegion> FindCommonRegions(
        DnaSequence sequence1, DnaSequence sequence2, int minLength, CommonRegionMode mode,
        StrandSelection strands = StrandSelection.Forward, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return CommonRegionsCore(() => sequence1.SuffixTree, new Lazy<MatchingStatisticsAutomaton>(() => AutomatonOf(sequence1, cancellationToken)), sequence2.Sequence, minLength, mode, strands, cancellationToken);
    }

    /// <summary>
    /// <see cref="FindCommonRegions(DnaSequence, DnaSequence, int, CommonRegionMode, StrandSelection, CancellationToken)"/> on nucleotide <b>strings</b>
    /// (N runs, IUPAC codes, lower case; a genome assembly with gaps): only A/C/G/T match, any other symbol matches nothing and ends a region (MUMmer <c>-n</c>);
    /// an ordinary match next to an N run is reported exactly as without it. With the default mode this is the right-maximal set of
    /// <see cref="FindCommonRegions(DnaSequence, DnaSequence, int)"/>.
    /// </summary>
    /// <param name="sequence1">Reference (indexed) sequence.</param>
    /// <param name="sequence2">Query sequence.</param>
    /// <param name="minLength">Minimum region length; values below 1 are treated as 1.</param>
    /// <param name="mode">Which common regions to report (default right-maximal).</param>
    /// <param name="strands">Strands to search (default forward).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    public static IReadOnlyList<CommonRegion> FindCommonRegions(
        string sequence1, string sequence2, int minLength, CommonRegionMode mode = CommonRegionMode.RightMaximal,
        StrandSelection strands = StrandSelection.Forward, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        if (sequence1.Length == 0 || sequence2.Length == 0)
            return Array.Empty<CommonRegion>();
        return CommonRegionsCore(() => BuildMaskedTree(sequence1, BuildLetterMask(NucleotideAlphabet)), new Lazy<MatchingStatisticsAutomaton>(() => new MatchingStatisticsAutomaton(AutomatonText(sequence1.ToUpperInvariant()), AutomatonAlphabet, cancellationToken)), sequence2.ToUpperInvariant(), minLength, mode, strands, cancellationToken);
    }

    private static List<CommonRegion> CommonRegionsCore(
        Func<global::SuffixTree.SuffixTree> treeFactory, Lazy<MatchingStatisticsAutomaton> automaton, string queryUpper, int minLength, CommonRegionMode mode,
        StrandSelection strands, CancellationToken cancellationToken)
    {
        AlignmentCancellation.Poll(cancellationToken);
        if (!Enum.IsDefined(mode))
            throw new ArgumentOutOfRangeException(nameof(mode), mode, "Unknown common-region mode.");
        int effectiveMinLength = Math.Max(1, minLength);
        var tree = new Lazy<global::SuffixTree.SuffixTree>(treeFactory);
        var result = new List<CommonRegion>();
        if (strands.HasFlag(StrandSelection.Forward))
            result.AddRange(RegionsOnStrand(queryUpper, reverse: false));
        if (strands.HasFlag(StrandSelection.ReverseComplement))
            result.AddRange(RegionsOnStrand(ReverseComplementAcgt(queryUpper), reverse: true));
        return result;

        IEnumerable<CommonRegion> RegionsOnStrand(string strandQuery, bool reverse)
        {
            IEnumerable<CommonRegion> regions;
            switch (mode)
            {
                case CommonRegionMode.RightMaximal:
                    // R14-4a: the suffix automaton of the text (no suffix tree); a region of the reverse strand is converted below like the others
                    regions = AutomatonSpans(automaton.Value, strandQuery, effectiveMinLength, reverse: false, cancellationToken).Select(span => SpanToRegion(span, strandQuery));
                    break;
                case CommonRegionMode.MaximalExactMatches:
                    {
                        string query = QueryText(strandQuery);
                        regions = FromMatches(tree.Value.FindMaximalExactMatches(query, effectiveMinLength), query);
                        break;
                    }
                case CommonRegionMode.MaximalUniqueReference:
                    {
                        string query = QueryText(strandQuery);
                        regions = FromMatches(tree.Value.FindMaximalUniqueMatches(query, effectiveMinLength, global::SuffixTree.MumUniqueness.Reference), query);
                        break;
                    }
                case CommonRegionMode.MaximalUnique:
                    {
                        string query = QueryText(strandQuery);
                        regions = FromMatches(tree.Value.FindMaximalUniqueMatches(query, effectiveMinLength, global::SuffixTree.MumUniqueness.Both), query);
                        break;
                    }
                default:
                    throw new InvalidOperationException("Unknown common-region mode.");
            }

            foreach (var region in regions)
            {
                AlignmentCancellation.Poll(cancellationToken);
                yield return reverse
                    ? region.WithPositions(region.PositionInFirst, strandQuery.Length - region.PositionInSecond - region.Length, isReverseComplement: true)
                    : region;
            }
        }

        static IEnumerable<CommonRegion> FromMatches(IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> matches, string query)
            => matches.Select(m => new CommonRegion(query, m.PositionInQuery, m.Length, m.PositionInText, m.PositionInQuery, isReverseComplement: false));
    }

    #endregion

    #region Right-maximal common regions by suffix automaton (R14-4a)

    // Symbols of the automaton: A C G T = 0..3 and one separator code (4) for every other symbol of the text; a query symbol outside A/C/G/T is code -1,
    // which no text symbol equals, so a match never contains it or a text separator (MUMmer -n: any other symbol matches nothing).
    private const int AutomatonAlphabet = 5;

    private static byte[] AutomatonText(string upper)
    {
        var codes = new byte[upper.Length];
        for (int i = 0; i < codes.Length; i++)
        {
            codes[i] = upper[i] switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, _ => 4 };
        }

        return codes;
    }

    private static sbyte[] AutomatonQuery(string upper)
    {
        var codes = new sbyte[upper.Length];
        for (int i = 0; i < codes.Length; i++)
        {
            codes[i] = upper[i] switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, _ => -1 };
        }

        return codes;
    }

    // The automaton of a DnaSequence is cached with the sequence (a reference searched with many queries pays the O(n) build once; the cache lives as long as the sequence).
    private static readonly System.Runtime.CompilerServices.ConditionalWeakTable<DnaSequence, MatchingStatisticsAutomaton> AutomatonCache = new();

    private static MatchingStatisticsAutomaton AutomatonOf(DnaSequence sequence, CancellationToken cancellationToken)
        => AutomatonCache.GetValue(sequence, s => new MatchingStatisticsAutomaton(AutomatonText(s.Sequence), AutomatonAlphabet, cancellationToken));

    // Eager matching statistics (cancellation is polled there), lazy enumeration of the distinct regions in order of their first start in the query.
    private static IEnumerable<CommonRegionSpan> AutomatonSpans(MatchingStatisticsAutomaton automaton, string queryUpper, int minLength, bool reverse, CancellationToken cancellationToken)
    {
        var (lengths, states) = automaton.Match(AutomatonQuery(queryUpper), cancellationToken);
        int queryLength = queryUpper.Length;
        return Enumerate();

        IEnumerable<CommonRegionSpan> Enumerate()
        {
            var seen = new HashSet<(int State, int Length)>();   // not a packed long: state ^ length is 0 for every region of a poly-A pair
            for (int i = 0; i < lengths.Length; i++)
            {
                if ((i & 0x3FF) == 0)
                    AlignmentCancellation.Poll(cancellationToken);
                int length = lengths[i];
                if (length < minLength)
                    continue;
                // a (state, length) pair is one distinct string
                if (!seen.Add((states[i], length)))
                    continue;
                yield return new CommonRegionSpan(automaton.FirstOccurrence(states[i]), reverse ? queryLength - i - length : i, length, reverse);
            }
        }
    }

    private static CommonRegion SpanToRegion(CommonRegionSpan span, string strandQuery)
        => new(strandQuery.Substring(span.IsReverseComplement ? strandQuery.Length - span.PositionInSecond - span.Length : span.PositionInSecond, span.Length), span.PositionInFirst, span.PositionInSecond, span.IsReverseComplement);

    /// <summary>
    /// The right-maximal common regions of <see cref="FindCommonRegions(DnaSequence, DnaSequence, int)"/> as <b>coordinates only</b> (R14-4a): no string is built, so the cost is
    /// O(n + m + R) however large the regions are (a poly-A pair of length n has n regions of average length n / 2, 200 M characters at n = 20 000, which the string result must
    /// materialise). Each region is the substring of the first sequence at <see cref="CommonRegionSpan.PositionInFirst"/> (its lowest occurrence) of <see cref="CommonRegionSpan.Length"/>
    /// characters; the order, the positions and the reverse-strand convention are those of <see cref="FindCommonRegions(DnaSequence, DnaSequence, int, CommonRegionMode, StrandSelection, CancellationToken)"/>
    /// with the default mode (forward regions by their start in the second sequence, then the reverse-complement ones).
    /// </summary>
    /// <param name="sequence1">Reference (indexed) sequence.</param>
    /// <param name="sequence2">Query sequence.</param>
    /// <param name="minLength">Minimum region length; values below 1 are treated as 1.</param>
    /// <param name="strands">Strands to search (default forward).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    public static IReadOnlyList<CommonRegionSpan> FindCommonRegionSpans(
        DnaSequence sequence1, DnaSequence sequence2, int minLength, StrandSelection strands = StrandSelection.Forward, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return SpansCore(() => AutomatonOf(sequence1, cancellationToken), sequence2.Sequence, minLength, strands, cancellationToken);
    }

    /// <summary>
    /// <see cref="FindCommonRegionSpans(DnaSequence, DnaSequence, int, StrandSelection, CancellationToken)"/> on nucleotide <b>strings</b> (N runs, IUPAC codes, lower case): only A/C/G/T match, any other
    /// symbol matches nothing and ends a region.
    /// </summary>
    /// <param name="sequence1">Reference (indexed) sequence.</param>
    /// <param name="sequence2">Query sequence.</param>
    /// <param name="minLength">Minimum region length; values below 1 are treated as 1.</param>
    /// <param name="strands">Strands to search (default forward).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    public static IReadOnlyList<CommonRegionSpan> FindCommonRegionSpans(
        string sequence1, string sequence2, int minLength, StrandSelection strands = StrandSelection.Forward, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        if (sequence1.Length == 0 || sequence2.Length == 0)
            return Array.Empty<CommonRegionSpan>();
        return SpansCore(() => new MatchingStatisticsAutomaton(AutomatonText(sequence1.ToUpperInvariant()), AutomatonAlphabet, cancellationToken), sequence2.ToUpperInvariant(), minLength, strands, cancellationToken);
    }

    private static List<CommonRegionSpan> SpansCore(Func<MatchingStatisticsAutomaton> automaton, string queryUpper, int minLength, StrandSelection strands, CancellationToken cancellationToken)
    {
        AlignmentCancellation.Poll(cancellationToken);
        int effectiveMinLength = Math.Max(1, minLength);
        var built = automaton();
        var result = new List<CommonRegionSpan>();
        if (strands.HasFlag(StrandSelection.Forward))
            result.AddRange(AutomatonSpans(built, queryUpper, effectiveMinLength, reverse: false, cancellationToken));
        if (strands.HasFlag(StrandSelection.ReverseComplement))
            result.AddRange(AutomatonSpans(built, ReverseComplementAcgt(queryUpper), effectiveMinLength, reverse: true, cancellationToken));
        return result;
    }

    #endregion

    #region Longest common substring with k mismatches

    /// <summary>
    /// Longest common substring with at most <paramref name="maxMismatches"/> mismatches (the k-mismatch LCS of Flouri, Giaquinta, Kobert &amp; Ukkonen 2015 and
    /// Thankachan, Apostolico &amp; Aluru 2016): the longest pair of equal-length substrings, one of each sequence, that differ in at most k positions (Hamming
    /// distance, no gaps).
    /// </summary>
    /// <remarks>
    /// <para>
    /// With k = 0 this is <see cref="FindLongestCommonRegion(DnaSequence, DnaSequence)"/>. For every diagonal (shift between the two sequences) the longest
    /// window with ≤ k mismatches is found with two pointers; time O(n·m), space O(1) beyond the result. Ties: the window whose start in the second sequence is
    /// smallest, then whose start in the first sequence is smallest (for k = 0 this is the documented tie rule of the exact method). The window can start or end on a
    /// mismatching pair (it is the longest, not the longest exactly-matching-at-the-ends).
    /// </para>
    /// <para>
    /// A pair of symbols agrees only when both are the same A/C/G/T letter: with the string overloads an N (or any other symbol) is a mismatch against everything,
    /// itself included. With <see cref="StrandSelection.ReverseComplement"/> the second sequence is reverse-complemented first and
    /// <see cref="ApproximateCommonRegion.PositionInSecond"/> is the forward-strand start of the window (whose reverse complement is
    /// <see cref="ApproximateCommonRegion.Sequence2"/>); with <see cref="StrandSelection.Both"/> the longer strand wins, a tie goes to the forward strand.
    /// Oracle: brute force over every (i, j, length) — docs/Evidence/GENOMIC-COMMON-001-Evidence.md.
    /// </para>
    /// </remarks>
    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="maxMismatches">Largest number of mismatching positions, ≥ 0.</param>
    /// <param name="strands">Strands of the second sequence to search (default forward).</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static ApproximateCommonRegion FindLongestCommonRegionWithMismatches(
        DnaSequence sequence1, DnaSequence sequence2, int maxMismatches, StrandSelection strands = StrandSelection.Forward)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return KMismatchCore(sequence1.Sequence, sequence2.Sequence, maxMismatches, strands);
    }

    /// <summary>
    /// <see cref="FindLongestCommonRegionWithMismatches(DnaSequence, DnaSequence, int, StrandSelection)"/> on nucleotide <b>strings</b>
    /// (case-insensitive; every symbol other than A/C/G/T is a mismatch with everything).
    /// </summary>
    /// <param name="sequence1">First sequence.</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="maxMismatches">Largest number of mismatching positions, ≥ 0.</param>
    /// <param name="strands">Strands of the second sequence to search (default forward).</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static ApproximateCommonRegion FindLongestCommonRegionWithMismatches(
        string sequence1, string sequence2, int maxMismatches, StrandSelection strands = StrandSelection.Forward)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return KMismatchCore(sequence1.ToUpperInvariant(), sequence2.ToUpperInvariant(), maxMismatches, strands);
    }

    private static ApproximateCommonRegion KMismatchCore(string s1, string s2, int maxMismatches, StrandSelection strands)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        var best = ApproximateCommonRegion.None;
        if (strands.HasFlag(StrandSelection.Forward))
            best = BestWindow(s1, s2, maxMismatches, reverse: false);
        if (strands.HasFlag(StrandSelection.ReverseComplement))
        {
            var reverse = BestWindow(s1, ReverseComplementAcgt(s2), maxMismatches, reverse: true);
            if (reverse.Length > best.Length)
                best = reverse;
        }

        return best;

        static ApproximateCommonRegion BestWindow(string a, string b, int k, bool reverse)
        {
            int n = a.Length, m = b.Length;
            int bestLength = 0, bestI = 0, bestJ = 0, bestMismatches = 0;
            for (int d = -(n - 1); d <= m - 1; d++)
            {
                // diagonal d: pairs (i, i + d); two pointers over the window [lo, hi] of i with at most k disagreements
                int i0 = Math.Max(0, -d), i1 = Math.Min(n - 1, m - 1 - d);
                int lo = i0, mismatches = 0;
                for (int hi = i0; hi <= i1; hi++)
                {
                    if (!Agree(a[hi], b[hi + d]))
                        mismatches++;
                    while (mismatches > k)
                    {
                        if (!Agree(a[lo], b[lo + d]))
                            mismatches--;
                        lo++;
                    }

                    int length = hi - lo + 1;
                    int j = lo + d;
                    if (length > bestLength || (length == bestLength && length > 0 && (j < bestJ || (j == bestJ && lo < bestI))))
                    {
                        bestLength = length;
                        bestI = lo;
                        bestJ = j;
                        bestMismatches = mismatches;
                    }
                }
            }

            if (bestLength == 0)
                return ApproximateCommonRegion.None;
            string w1 = a.Substring(bestI, bestLength);
            string w2 = b.Substring(bestJ, bestLength);
            if (!reverse)
                return new ApproximateCommonRegion(w1, w2, bestI, bestJ, bestLength, bestMismatches);
            // b is the reverse complement of the original second sequence: report the window in forward-strand coordinates
            return new ApproximateCommonRegion(w1, w2, bestI, m - bestJ - bestLength, bestLength, bestMismatches, IsReverseComplement: true);
        }

        static bool Agree(char x, char y) => x == y && x is 'A' or 'C' or 'G' or 'T';
    }

    #endregion

    #region Longest common substring of several sequences

    /// <summary>
    /// Longest substring shared by at least <paramref name="minSupport"/> of the sequences (the k-common substring problem, Gusfield 1997 §7.6; Rosalind LCSM for
    /// <c>minSupport = N</c>), via the canonical generalized suffix tree <see cref="SuffixTree.SuffixTree.FindLongestCommonSubstrings(IReadOnlyList{string}, int?)"/>.
    /// Every distinct longest string is returned (ordinal order) with its first occurrence in each sequence.
    /// </summary>
    /// <param name="sequences">At least two sequences.</param>
    /// <param name="minSupport">Number of sequences that must contain the substring, 2 … N (default all N).</param>
    /// <returns>The longest common strings; empty when no non-empty string is shared by <paramref name="minSupport"/> sequences.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> or an element is null.</exception>
    /// <exception cref="ArgumentException">Fewer than two sequences.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minSupport"/> outside [2, N].</exception>
    public static IReadOnlyList<MultiCommonRegion> FindLongestCommonRegion(IReadOnlyList<DnaSequence> sequences, int? minSupport = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        RequireNoNull(sequences);
        return MultiCommonRegionsCore(sequences.Select(s => s.Sequence).ToList(), minSupport);
    }

    /// <summary>
    /// <see cref="FindLongestCommonRegion(IReadOnlyList{DnaSequence}, int?)"/> on nucleotide <b>strings</b> (case-insensitive): only A/C/G/T match, N and every other
    /// symbol ends a region (the runs of non-matching symbols of different sequences are unique separators of the generalized suffix tree).
    /// </summary>
    /// <param name="sequences">At least two sequences.</param>
    /// <param name="minSupport">Number of sequences that must contain the substring, 2 … N (default all N).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> or an element is null.</exception>
    /// <exception cref="ArgumentException">Fewer than two sequences.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minSupport"/> outside [2, N].</exception>
    public static IReadOnlyList<MultiCommonRegion> FindLongestCommonRegion(IReadOnlyList<string> sequences, int? minSupport = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        RequireNoNull(sequences);
        var letters = BuildLetterMask(NucleotideAlphabet);
        int next = FirstSeparator;
        return MultiCommonRegionsCore(sequences.Select(s => MaskText(s, letters, ref next)).ToList(), minSupport);
    }

    private static void RequireNoNull<T>(IReadOnlyList<T> sequences) where T : class
    {
        foreach (var s in sequences)
            if (s is null)
                throw new ArgumentNullException(nameof(sequences), "A sequence in the list is null.");
    }

    private static IReadOnlyList<MultiCommonRegion> MultiCommonRegionsCore(List<string> texts, int? minSupport)
    {
        if (texts.Count < 2)
            throw new ArgumentException("At least two sequences are required.");
        int support = minSupport ?? texts.Count;
        if (support < 2 || support > texts.Count)
            throw new ArgumentOutOfRangeException(nameof(minSupport), minSupport, $"minSupport must be between 2 and the number of sequences ({texts.Count}).");

        return global::SuffixTree.SuffixTree.FindLongestCommonSubstrings(texts, support)
            .Select(w => new MultiCommonRegion(w, texts.Select(t => t.IndexOf(w, StringComparison.Ordinal)).ToArray()))
            .ToList();
    }

    #endregion
}

/// <summary>Which common regions <see cref="GenomicAnalyzer.FindCommonRegions(DnaSequence, DnaSequence, int, CommonRegionMode, StrandSelection, CancellationToken)"/> reports.</summary>
public enum CommonRegionMode
{
    /// <summary>
    /// For each start of the second sequence the single longest substring that also occurs in the first (matching statistics, Chang &amp; Lawler 1994),
    /// distinct strings once, with the first occurrence in the first sequence — the set of
    /// <see cref="GenomicAnalyzer.FindCommonRegions(DnaSequence, DnaSequence, int)"/>.
    /// </summary>
    RightMaximal = 0,

    /// <summary>
    /// Every maximal exact match (MEM): one entry per pair of occurrences (position in the first, position in the second) that cannot be extended
    /// on either side — <c>mummer -maxmatch</c>, canonical <c>SuffixTree.FindMaximalExactMatches</c>.
    /// </summary>
    MaximalExactMatches = 1,

    /// <summary>
    /// Maximal matches whose string occurs exactly once in the first (reference) sequence and may repeat in the second — <c>mummer -mumreference</c>
    /// (MUM candidates / MAMs), canonical <c>SuffixTree.FindMaximalUniqueMatches(Reference)</c>.
    /// </summary>
    MaximalUniqueReference = 2,

    /// <summary>
    /// Maximal unique matches (MUMs): the string occurs exactly once in each sequence — <c>mummer -mum</c> (Delcher et al. 1999),
    /// canonical <c>SuffixTree.FindMaximalUniqueMatches(Both)</c>.
    /// </summary>
    MaximalUnique = 3,
}

/// <summary>
/// A longest common substring with mismatches (<see cref="GenomicAnalyzer.FindLongestCommonRegionWithMismatches(DnaSequence, DnaSequence, int, StrandSelection)"/>):
/// the window <see cref="Sequence1"/> of the first sequence and the equal-length window <see cref="Sequence2"/> of the second (or of its reverse complement), which differ in
/// <see cref="Mismatches"/> positions.
/// </summary>
public readonly record struct ApproximateCommonRegion(
    string Sequence1, string Sequence2, int PositionInFirst, int PositionInSecond, int Length, int Mismatches, bool IsReverseComplement = false)
{
    /// <summary>The empty result (no common window; also when one sequence is empty).</summary>
    public static readonly ApproximateCommonRegion None = new(string.Empty, string.Empty, -1, -1, 0, 0);

    /// <summary>True for <see cref="None"/>.</summary>
    public bool IsEmpty => Length == 0;
}

/// <summary>
/// A longest common substring of several sequences (<see cref="GenomicAnalyzer.FindLongestCommonRegion(IReadOnlyList{DnaSequence}, int?)"/>): the string and its first
/// 0-based occurrence in every input sequence (−1 where the sequence does not contain it, possible when fewer than all sequences are required).
/// </summary>
public readonly record struct MultiCommonRegion(string Sequence, IReadOnlyList<int> Positions)
{
    /// <summary>Length of the common string.</summary>
    public int Length => Sequence.Length;

    /// <summary>Number of sequences that contain the string.</summary>
    public int Support => Positions.Count(p => p >= 0);
}
