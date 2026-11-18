namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Genomic analysis algorithms based on suffix trees.
/// Provides efficient solutions for common bioinformatics problems.
/// </summary>
public static class GenomicAnalyzer
{
    #region Repeat Finding

    /// <summary>
    /// Finds the longest repeated substring (LRS) of a DNA sequence: the longest
    /// substring that occurs at least twice (occurrences may overlap). This is the
    /// classical suffix-tree application — the deepest internal node, whose path label
    /// spells a longest substring occurring ≥ 2 times (CMU 15-451 Lecture #10 §2.1;
    /// Wikipedia, Longest repeated substring problem). Returns <see cref="RepeatInfo.None"/>
    /// when no substring repeats (including the empty sequence).
    /// Time complexity: O(n) for the deepest-node query after O(n) suffix-tree construction
    /// (Ukkonen); occurrence enumeration adds O(occ). See docs/algorithms/Repeat_Analysis/Repeat_Detection.md.
    /// </summary>
    public static RepeatInfo FindLongestRepeat(DnaSequence sequence)
    {
        var tree = sequence.SuffixTree;
        string lrs = tree.LongestRepeatedSubstring();

        if (string.IsNullOrEmpty(lrs))
        {
            return RepeatInfo.None;
        }

        var positions = tree.FindAllOccurrences(lrs).OrderBy(p => p).ToList();
        return new RepeatInfo(lrs, positions);
    }

    /// <summary>
    /// Finds every distinct substring that occurs at least twice and has length ≥ <paramref name="minLength"/>.
    /// A substring occurs ≥ 2 times iff it ends on, or at, an edge into an internal node of the suffix tree:
    /// each internal node with path-label length d whose parent has path-label length p carries the repeats
    /// of lengths p + 1 … d (all with the node's occurrence set), so the method enumerates exactly those
    /// loci with length ≥ <paramref name="minLength"/> (Gusfield 1997 §7.1; CMU 15-451 Lecture #10 §2.1;
    /// GeeksforGeeks, Suffix Tree Application 3). Results are ordered as the sorted suffixes would
    /// produce them: by the lexicographically smallest suffix that starts with the repeat, then by
    /// increasing length; positions are ascending. Time O(n) for the traversal plus
    /// O(Σ (label length + occurrences)) for the output, which is the size of the result itself
    /// (the earlier implementation materialised and sorted all n suffixes, O(n²) time and memory).
    /// See docs/algorithms/Repeat_Analysis/Repeat_Detection.md.
    /// </summary>
    public static IEnumerable<RepeatInfo> FindRepeats(DnaSequence sequence, int minLength)
    {
        var tree = sequence.SuffixTree;

        // A repeat is a NON-EMPTY substring occurring >=2 times (the path label of an internal
        // node, never the root). The empty string is not a repeat, so the effective minimum is
        // at least 1 even when the caller passes minLength <= 0. (CMU 15-451 Lecture #10 §2.1.)
        int effectiveMinLength = Math.Max(1, minLength);

        // Pre-order traversal with children in ascending key order (the terminator, a leaf, sorts first):
        // internal nodes come out in the order of their lexicographically smallest suffix, and a node
        // always precedes the internal nodes below it.
        var loci = new RepeatLocusCollector(effectiveMinLength);
        tree.Traverse(loci);

        foreach (var (edgeEnd, parentDepth, depth) in loci.Nodes)
        {
            // Path label of the node: the `depth` characters of the text ending at the end of its incoming edge
            // (an internal node's edge never contains the terminator).
            string label = tree.Text.Substring(edgeEnd - depth, depth);
            var positions = new List<int>(tree.FindAllOccurrences(label));
            positions.Sort();

            for (int length = Math.Max(parentDepth + 1, effectiveMinLength); length < depth; length++)
                yield return new RepeatInfo(label.Substring(0, length), new List<int>(positions));

            yield return new RepeatInfo(label, positions);
        }
    }

    /// <summary>
    /// Collects the internal (non-root) nodes of a suffix tree whose path-label length is at least
    /// <c>minLength</c>, in traversal order, as (end of incoming edge in the text, parent path-label length, own path-label length).
    /// </summary>
    private sealed class RepeatLocusCollector(int minLength) : global::SuffixTree.ISuffixTreeVisitor
    {
        private bool _rootSeen;

        public List<(int EdgeEnd, int ParentDepth, int Depth)> Nodes { get; } = new();

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (!_rootSeen)
            {
                _rootSeen = true; // the root has no incoming edge and is not a repeat
                return;
            }
            if (childCount == 0)
                return; // a leaf is a single suffix, not a repeat

            int stringDepth = depth + (endIndex - startIndex);
            if (stringDepth >= minLength)
                Nodes.Add((endIndex, depth, stringDepth));
        }

        public void EnterBranch(int key) { }

        public void ExitBranch() { }
    }

    /// <summary>
    /// Finds tandem repeats (consecutive repeating units like "ATGATGATG").
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minUnitLength">
    /// Minimum candidate repeat-unit length. Must be &gt;= 1: a zero-length unit is degenerate (the empty
    /// string "tiles" every position infinitely and never advances the scan), so it is rejected to avoid a
    /// non-terminating loop rather than silently hanging.
    /// </param>
    /// <param name="minRepetitions">
    /// Minimum number of consecutive unit copies. Must be &gt;= 2: a tandem repeat requires at least two
    /// copies (<c>k &gt;= 2</c>), and a value of <c>0</c> would divide by zero in the unit-length bound while
    /// a value of <c>1</c> would report every substring as a trivial single-copy "repeat".
    /// </param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="minUnitLength"/> is &lt; 1, or <paramref name="minRepetitions"/> is &lt; 2.
    /// </exception>
    public static IEnumerable<TandemRepeat> FindTandemRepeats(DnaSequence sequence, int minUnitLength = 2, int minRepetitions = 2)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minUnitLength < 1)
            throw new ArgumentOutOfRangeException(nameof(minUnitLength), minUnitLength,
                "Minimum repeat-unit length must be at least 1; a zero-length unit never terminates the scan.");
        if (minRepetitions < 2)
            throw new ArgumentOutOfRangeException(nameof(minRepetitions), minRepetitions,
                "A tandem repeat requires at least two consecutive copies; minRepetitions must be >= 2.");

        return FindTandemRepeatsCore(sequence, minUnitLength, minRepetitions);
    }

    private static IEnumerable<TandemRepeat> FindTandemRepeatsCore(DnaSequence sequence, int minUnitLength, int minRepetitions)
    {
        string seq = sequence.Sequence;

        for (int unitLen = minUnitLength; unitLen <= seq.Length / minRepetitions; unitLen++)
        {
            // `start` advances by 1 normally, but jumps past a found tandem — a while loop
            // makes the in-body skip explicit (a for loop would flag S127).
            int start = 0;
            while (start <= seq.Length - unitLen * minRepetitions)
            {
                string unit = seq.Substring(start, unitLen);
                int repetitions = 1;
                int pos = start + unitLen;

                while (pos + unitLen <= seq.Length &&
                       seq.Substring(pos, unitLen) == unit)
                {
                    repetitions++;
                    pos += unitLen;
                }

                if (repetitions >= minRepetitions)
                {
                    yield return new TandemRepeat(unit, start, repetitions);
                    start = pos - unitLen; // Skip to end of this tandem
                }

                start++;
            }
        }
    }

    #endregion

    #region Motif Finding

    /// <summary>
    /// Finds all occurrences of an exact motif (pattern) in a sequence: the 0-based start of every
    /// window <c>S[i..i+m-1]</c> equal to the motif, overlapping occurrences included, in
    /// <b>ascending</b> order. Matching is case-insensitive and exact: an IUPAC ambiguity code in the
    /// motif is a literal symbol (a DNA sequence never contains one, so it matches nothing) and only
    /// the forward strand is searched. Use <see cref="FindKnownMotifHits"/> for degenerate (IUPAC)
    /// motifs and both-strand searches. A null or empty motif yields an empty list (the empty string
    /// is not a motif); a null sequence throws <see cref="ArgumentNullException"/>.
    /// <para>
    /// Delegates to the canonical <see cref="MotifFinder.FindExactMotif"/> (same engine, same order;
    /// the two used to be copies that differed only in sorting and null handling). That engine is the
    /// sequence's cached suffix tree: an O(n) build with a large constant (measured 0.9 s / 170 MB at
    /// 1 Mbp, 23 s / 1.7 GB at 10 Mbp) paid once per <see cref="DnaSequence"/>, then O(m + occ) per
    /// query. For a single one-off query over a very large sequence the linear scan
    /// <see cref="KmerAnalyzer.FindKmerPositions"/> (KMP, 71 ms at 10 Mbp, no index) is the cheaper
    /// call and returns the same positions; <see cref="FindKnownMotifs(DnaSequence, IEnumerable{string})"/>
    /// chooses between the two itself.
    /// </para>
    /// </summary>
    public static IReadOnlyList<int> FindMotif(DnaSequence sequence, string motif)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return MotifFinder.FindExactMotif(sequence, motif).ToList();
    }

    /// <summary>
    /// Finds palindromic sequences (restriction enzyme recognition sites).
    /// A DNA palindrome reads the same 5'→3' on both strands.
    /// Example: GAATTC (EcoRI site) - complement is CTTAAG, reversed is GAATTC.
    /// </summary>
    public static IEnumerable<PalindromeInfo> FindPalindromes(DnaSequence sequence, int minLength = 4, int maxLength = 12)
    {
        string seq = sequence.Sequence;

        for (int len = minLength; len <= Math.Min(maxLength, seq.Length); len += 2) // Palindromes must be even length
        {
            for (int start = 0; start <= seq.Length - len; start++)
            {
                string subseq = seq.Substring(start, len);
                string revComp = DnaSequence.GetReverseComplementString(subseq);

                if (subseq == revComp)
                {
                    yield return new PalindromeInfo(subseq, start);
                }
            }
        }
    }

    /// <summary>
    /// Number of distinct (strand-expanded) queries from which a sequence's suffix tree is built for
    /// <see cref="FindKnownMotifs(DnaSequence, IEnumerable{string})"/> instead of scanning the text once
    /// per query. Measured (Release, random DNA): tree build = 73 (100 kb), 135 (1 Mbp) and 197 (4 Mbp)
    /// KMP scans, so below 128 queries the linear scans are cheaper and need no index memory.
    /// </summary>
    private const int IndexedSearchMinQueries = 128;

    /// <summary>
    /// Searches a sequence for a set of known motifs, returning, for each motif that occurs,
    /// the 0-based start positions of <b>all</b> its occurrences. This is the classical exact
    /// set-matching problem ("find all occurrences of each pattern P in text T") (Gusfield 1997,
    /// ISBN 0-521-58519-8; exact-matching definition, Tufts COMP 150GEN exact.html).
    /// <para>
    /// <b>Overlapping occurrences are all reported</b> — e.g. motif "AAA" in "AAAAA" yields
    /// {0, 1, 2}, mirroring Biopython's overlap-aware semantics (<c>Seq("AAAA").count_overlap("AA") == 3</c>).
    /// </para>
    /// Motifs are upper-cased before searching (DNA is processed upper-cased, per Biopython
    /// <c>Bio.Seq</c>); the result is keyed by the upper-cased motif. Empty or whitespace-only
    /// motifs are skipped (the empty string is not a motif; mirrors <see cref="FindMotif"/>), and a
    /// motif with no occurrence is omitted from the result. Positions for each motif are returned
    /// <b>sorted ascending</b>. Duplicate motifs that normalize to the same upper-cased key collapse to
    /// a single entry. Forward strand only, IUPAC symbols literal (see <see cref="FindKnownMotifHits"/>
    /// for degenerate motifs and both strands).
    /// <para>
    /// Complexity: with fewer than 128 distinct motifs each one is a linear KMP scan of the text
    /// (<see cref="KmerAnalyzer.FindKmerPositions"/>, O(n + |m| + occ), no index); from 128 motifs on the
    /// sequence's cached suffix tree is built once (O(n), Ukkonen) and queried in O(|m| + occ) with an
    /// O(occ·log occ) sort. Both routes return identical positions. See
    /// docs/algorithms/Motif_Analysis/Known_Motif_Search.md.
    /// </para>
    /// </summary>
    public static Dictionary<string, IReadOnlyList<int>> FindKnownMotifs(
        DnaSequence sequence,
        IEnumerable<string> motifs)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(motifs);

        var found = SearchKnownMotifs(sequence, motifs, bothStrands: false, allowIupac: false);
        var result = new Dictionary<string, IReadOnlyList<int>>(found.Count);
        foreach (var (key, (forward, _)) in found)
            result[key] = forward;
        return result;
    }

    /// <summary>
    /// Strand-aware known-motif search with optional IUPAC degeneracy, the generalisation of
    /// <see cref="FindKnownMotifs(DnaSequence, IEnumerable{string})"/> that EMBOSS <c>fuzznuc</c> offers
    /// (<c>-pattern</c> with IUPAC codes, <c>-complement</c>). For every motif that occurs the result
    /// holds its hits sorted by position (a forward hit before a reverse hit at the same position).
    /// <list type="bullet">
    /// <item>Forward hit: <c>S[p..p+m-1]</c> matches the motif.</item>
    /// <item>Reverse hit (<paramref name="bothStrands"/>): the window <c>S[p..p+m-1]</c> matches the
    /// <i>reverse complement</i> of the motif, i.e. the motif occurs on the opposite strand reading 5'→3'.
    /// <see cref="KnownMotifHit.Position"/> is always the 0-based start of the window on the <b>forward</b>
    /// strand (fuzznuc/GFF convention: reverse-strand coordinates are forward-strand coordinates), so a
    /// palindromic site such as EcoRI GAATTC is reported on both strands at the same position.</item>
    /// <item><paramref name="allowIupac"/>: each motif position is an IUPAC NC-IUB 1984 code matching its
    /// base set (R = A|G, W = A|T, …; <see cref="MotifFinder.FindDegenerateMotif(DnaSequence, string)"/>);
    /// the reverse-strand pattern is the IUPAC reverse complement (R↔Y, K↔M, B↔V, D↔H, S/W/N self).
    /// Characters that are not IUPAC DNA codes throw <see cref="ArgumentException"/>. Without it
    /// the motif is exact and ambiguity codes are literal.</item>
    /// </list>
    /// Case-insensitive, overlapping hits included, null/empty/whitespace motifs skipped, results keyed by
    /// the upper-cased motif, duplicates collapsed, motifs without a hit omitted.
    /// Verified against EMBOSS fuzznuc 6.6.0 (<c>-pmismatch 0</c>, with and without <c>-complement</c>).
    /// </summary>
    /// <param name="sequence">Subject sequence.</param>
    /// <param name="motifs">Motifs to search.</param>
    /// <param name="bothStrands">Also report occurrences on the reverse strand.</param>
    /// <param name="allowIupac">Interpret IUPAC ambiguity codes in the motifs.</param>
    public static Dictionary<string, IReadOnlyList<KnownMotifHit>> FindKnownMotifHits(
        DnaSequence sequence,
        IEnumerable<string> motifs,
        bool bothStrands = false,
        bool allowIupac = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(motifs);

        var found = SearchKnownMotifs(sequence, motifs, bothStrands, allowIupac);
        var result = new Dictionary<string, IReadOnlyList<KnownMotifHit>>(found.Count);
        foreach (var (key, (forward, reverse)) in found)
            result[key] = MergeStrandHits(forward, reverse);
        return result;
    }

    private static Dictionary<string, (List<int> Forward, List<int> Reverse)> SearchKnownMotifs(
        DnaSequence sequence, IEnumerable<string> motifs, bool bothStrands, bool allowIupac)
    {
        // The empty string is not a motif (an exact search for "" would match every position).
        // Skip empty/whitespace motifs, consistent with FindMotif's empty-pattern guard.
        var keys = new List<string>();
        var seen = new HashSet<string>(StringComparer.Ordinal);
        foreach (var motif in motifs)
        {
            if (string.IsNullOrWhiteSpace(motif))
                continue;
            string normalized = motif.ToUpperInvariant();
            if (seen.Add(normalized))
                keys.Add(normalized);
        }

        bool useIndex = !allowIupac && keys.Count * (bothStrands ? 2L : 1L) >= IndexedSearchMinQueries;
        var result = new Dictionary<string, (List<int>, List<int>)>(keys.Count);
        foreach (string key in keys)
        {
            var forward = FindMotifPositions(sequence, key, allowIupac, useIndex);
            var reverse = bothStrands
                ? FindMotifPositions(sequence, DnaSequence.GetReverseComplementString(key), allowIupac, useIndex)
                : new List<int>();
            if (forward.Count > 0 || reverse.Count > 0)
                result[key] = (forward, reverse);
        }

        return result;
    }

    private static List<int> FindMotifPositions(DnaSequence sequence, string motif, bool allowIupac, bool useIndex)
    {
        if (allowIupac)
        {
            // MotifFinder scans left to right: ascending, overlapping.
            return MotifFinder.FindDegenerateMotif(sequence, motif).Select(m => m.Position).ToList();
        }

        if (useIndex)
        {
            // Suffix-tree occurrence enumeration is DFS-order (not sorted).
            var positions = sequence.SuffixTree.FindAllOccurrences(motif);
            var sorted = new List<int>(positions);
            sorted.Sort();
            return sorted;
        }

        return KmerAnalyzer.FindKmerPositions(sequence.Sequence, motif).ToList();
    }

    private static List<KnownMotifHit> MergeStrandHits(List<int> forward, List<int> reverse)
    {
        var merged = new List<KnownMotifHit>(forward.Count + reverse.Count);
        int f = 0, r = 0;
        while (f < forward.Count || r < reverse.Count)
        {
            if (r >= reverse.Count || (f < forward.Count && forward[f] <= reverse[r]))
                merged.Add(new KnownMotifHit(forward[f++], IsReverseStrand: false));
            else
                merged.Add(new KnownMotifHit(reverse[r++], IsReverseStrand: true));
        }

        return merged;
    }

    #endregion

    #region Sequence Comparison

    /// <summary>
    /// Finds the longest common <b>substring</b> (a longest <i>contiguous</i> string that is a
    /// substring of both sequences), not a gapped subsequence — Wikipedia "Longest common substring";
    /// the contiguous distinction is explicit there. Computed via the generalized suffix tree of
    /// <paramref name="sequence1"/>: the deepest node whose subtree has leaves from both strings spells
    /// the answer (Gusfield 1997, ISBN 0-521-58519-8; GeeksforGeeks "Suffix Tree Application 5"),
    /// found by streaming <paramref name="sequence2"/> through the tree with suffix links
    /// (matching statistics, Chang &amp; Lawler 1994).
    /// <para>
    /// Ties (Wikipedia: "BADANAT"/"CANADAS" share "ADA" and "ANA"): the substring whose first match
    /// in <paramref name="sequence2"/> ends earliest is returned (documented tie-break of
    /// <see cref="SuffixTree.SuffixTree.FindAllLongestCommonSubstrings(string)"/>). Both positions are the
    /// <b>first (lowest) 0-based occurrence</b> of that substring in the respective sequence, so the
    /// result does not depend on the internal suffix-tree layout.
    /// </para>
    /// Returns <see cref="CommonRegion.None"/> when there is no shared substring (including any empty input).
    /// Sequences are the validated upper-case A/C/G/T strings of <see cref="DnaSequence"/> (lower case is
    /// folded at construction, IUPAC/other symbols are rejected there), so matching is exact and case-insensitive.
    /// Time complexity: O(n + m) with the generalized suffix tree (Gusfield 1997).
    /// See docs/algorithms/Sequence_Comparison/Common_Region_Detection.md.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence1"/> or <paramref name="sequence2"/> is null.</exception>
    public static CommonRegion FindLongestCommonRegion(DnaSequence sequence1, DnaSequence sequence2)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        var (lcs, positionsInFirst, positionsInSecond) =
            sequence1.SuffixTree.FindAllLongestCommonSubstrings(sequence2.Sequence);

        if (string.IsNullOrEmpty(lcs))
        {
            return CommonRegion.None;
        }

        // Position lists are ascending and duplicate-free: [0] is the first occurrence in each sequence.
        return new CommonRegion(lcs, positionsInFirst[0], positionsInSecond[0]);
    }

    /// <summary>
    /// For each start position in <paramref name="sequence2"/>, finds the <b>single longest</b> contiguous
    /// substring (per Wikipedia "Longest common substring") of length ≥ <paramref name="minLength"/> that
    /// also occurs in <paramref name="sequence1"/> — i.e. the right-maximal common substring at that start,
    /// whose length is the matching statistic ms(i) of Chang &amp; Lawler (1994; Gusfield 1997 ch. 7).
    /// Distinct such substrings are reported once, in order of their first start in <paramref name="sequence2"/>,
    /// with the <b>first (lowest) 0-based occurrence</b> in <paramref name="sequence1"/> and that start position
    /// in <paramref name="sequence2"/>.
    /// <para>
    /// This is the set of right-maximal common substrings keyed by start position in
    /// <paramref name="sequence2"/>, <b>not</b> the full set of every common substring: shorter prefixes
    /// that share a start position with a longer match are not reported (e.g. for <c>ACGTACGT</c> vs
    /// <c>TTACGTGG</c> with minLength 3 the result is <c>{TACGT, ACGT, CGT}</c> — the prefixes
    /// <c>TAC</c>, <c>TACG</c>, <c>ACG</c> of those matches are omitted).
    /// </para>
    /// A common substring is non-empty, so values of
    /// <paramref name="minLength"/> below 1 are treated as 1 (the empty string is not a region).
    /// <para>
    /// Algorithm: ms(i) values below K = ⌈log₄(n·m)⌉ + 2 are found by extending the match from ms(i−1) − 1
    /// (ms(i) ≥ ms(i−1) − 1) with <see cref="SuffixTree.SuffixTree.Contains(ReadOnlySpan{char})"/> probes of at most K
    /// characters; every ms(i) ≥ max(K, minLength) is recovered from the canonical suffix-link
    /// maximal-exact-match stream <see cref="SuffixTree.SuffixTree.FindMaximalExactMatches(string, int)"/>
    /// (MUMmer <c>-maxmatch</c>): the longest match at i lies on a diagonal whose maximal run is a MEM, so
    /// ms(i) = max over MEMs starting at or before i of (end − i). The expected number of MEMs of length ≥ K between
    /// unrelated random sequences is O(1), so no start position ever re-scans a long match.
    /// Time O(n + m·K + R) plus the cost of materialising the reported regions (the output itself is
    /// Σ|region|); R is the number of maximal matches of length ≥ max(K, minLength). No O(m²) substring is built
    /// (the former per-start binary search copied up to m characters per probe).
    /// See docs/algorithms/Sequence_Comparison/Common_Region_Detection.md.
    /// </para>
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence1"/> or <paramref name="sequence2"/> is null.</exception>
    public static IEnumerable<CommonRegion> FindCommonRegions(
        DnaSequence sequence1,
        DnaSequence sequence2,
        int minLength)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);

        // A common region is a NON-EMPTY contiguous substring; the empty string is not a region,
        // so the effective minimum length is at least 1 (Wikipedia "Longest common substring").
        int effectiveMinLength = Math.Max(1, minLength);
        int[] matchLengths = MatchingStatistics(sequence1.SuffixTree, sequence2.Sequence, sequence1.Length, effectiveMinLength);
        return EnumerateCommonRegions(sequence1.SuffixTree, sequence2.Sequence, matchLengths, effectiveMinLength);
    }

    private static IEnumerable<CommonRegion> EnumerateCommonRegions(
        global::SuffixTree.SuffixTree tree, string seq2, int[] matchLengths, int effectiveMinLength)
    {
        var found = new HashSet<string>(StringComparer.Ordinal);
        var foundLookup = found.GetAlternateLookup<ReadOnlySpan<char>>();

        for (int i = 0; i < matchLengths.Length; i++)
        {
            int length = matchLengths[i];
            if (length < effectiveMinLength)
            {
                continue;
            }

            ReadOnlySpan<char> match = seq2.AsSpan(i, length);
            if (!foundLookup.Add(match))
            {
                continue;
            }

            // First (lowest) occurrence in sequence1: FindAllOccurrences is in tree order, not sorted.
            var occurrences = tree.FindAllOccurrences(match);
            int first = occurrences[0];
            for (int k = 1; k < occurrences.Count; k++)
            {
                if (occurrences[k] < first)
                {
                    first = occurrences[k];
                }
            }

            yield return new CommonRegion(match.ToString(), first, i);
        }
    }

    /// <summary>
    /// Matching statistics of <paramref name="query"/> against the text of <paramref name="tree"/>:
    /// entry i is ms(i) = the length of the longest prefix of query[i..] occurring in the text — exact
    /// wherever it is ≥ <paramref name="minLength"/>, and a value &lt; <paramref name="minLength"/> elsewhere.
    /// </summary>
    private static int[] MatchingStatistics(global::SuffixTree.SuffixTree tree, string query, int textLength, int minLength)
    {
        int m = query.Length;
        var ms = new int[m];
        if (m == 0 || textLength == 0)
        {
            return ms;
        }

        // K ≈ log4(n·m) + 2: between random sequences the expected number of right-maximal matches of
        // length ≥ K is n·m/4^K ≤ 1/16, so the maximal-exact-match stream below is output-sized.
        int shortLimit = (int)Math.Ceiling(Math.Log((double)textLength + 1, 4) + Math.Log((double)m + 1, 4)) + 2;
        int memLength = Math.Max(minLength, shortLimit);

        bool needMems = true;
        if (minLength < shortLimit)
        {
            // Phase 1: exact ms(i) up to shortLimit by extension from ms(i-1)-1 (ms(i) >= ms(i-1)-1).
            needMems = false;
            int length = 0;
            for (int i = 0; i < m; i++)
            {
                if (length > 0)
                {
                    length--;
                }

                while (length < shortLimit && i + length < m && tree.Contains(query.AsSpan(i, length + 1)))
                {
                    length++;
                }

                ms[i] = length;
                needMems |= length >= shortLimit;
            }
        }

        if (!needMems)
        {
            return ms;
        }

        // Phase 2: ms(i) >= memLength from the maximal exact matches (sorted by query position).
        var mems = tree.FindMaximalExactMatches(query, memLength);
        int next = 0;
        int maxEnd = 0;
        for (int i = 0; i < m; i++)
        {
            while (next < mems.Count && mems[next].PositionInQuery <= i)
            {
                maxEnd = Math.Max(maxEnd, mems[next].PositionInQuery + mems[next].Length);
                next++;
            }

            if (minLength >= shortLimit || ms[i] >= shortLimit)
            {
                ms[i] = maxEnd - i;
            }
        }

        return ms;
    }

    // Default k-mer length. Mash uses k=21 for whole genomes; for short DNA
    // sequences a smaller default is appropriate. The value only sets the
    // resolution of the comparison and does not change the Jaccard formula.
    private const int DefaultKmerSize = 5;

    // Jaccard index is defined on [0, 1] (Jaccard 1901). This factor only
    // re-expresses that ratio as a percentage for reporting; it is not part
    // of the formal coefficient.
    private const double PercentScale = 100.0;

    /// <summary>
    /// Computes the k-mer Jaccard similarity index between two DNA sequences,
    /// reported as a percentage in [0, 100].
    /// </summary>
    /// <remarks>
    /// <para>
    /// Each sequence is decomposed into its <b>set</b> of distinct length-<paramref name="kmerSize"/>
    /// substrings (k-mers). The Jaccard index is the size of the intersection of the two
    /// k-mer sets divided by the size of their union:
    /// <c>J(A,B) = |A ∩ B| / |A ∪ B|</c> (Jaccard, 1901). Applied to k-mer sets this is the
    /// fraction of shared k-mers out of all distinct k-mers in the two sequences (Ondov et al.,
    /// Mash, 2016). The result is multiplied by 100 to report a percentage.
    /// </para>
    /// <para>
    /// The Jaccard index is undefined when the union is empty (Jaccard 1901 is stated for
    /// non-empty sets). This occurs only when both k-mer sets are empty (e.g. both sequences
    /// are empty or shorter than <paramref name="kmerSize"/>); this method returns 0 in that case.
    /// </para>
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="kmerSize">k-mer length; must be ≥ 1.</param>
    /// <returns>Jaccard similarity as a percentage in [0, 100].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence1"/> or <paramref name="sequence2"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> is less than 1.</exception>
    public static double CalculateSimilarity(DnaSequence sequence1, DnaSequence sequence2, int kmerSize = DefaultKmerSize)
    {
        ArgumentNullException.ThrowIfNull(sequence1);

        ArgumentNullException.ThrowIfNull(sequence2);

        if (kmerSize < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(kmerSize), kmerSize, "k-mer size must be at least 1.");
        }

        HashSet<string> kmers1 = GetKmers(sequence1.Sequence, kmerSize);
        HashSet<string> kmers2 = GetKmers(sequence2.Sequence, kmerSize);

        // |A ∪ B| = |A| + |B| − |A ∩ B|, computed on the distinct-k-mer sets.
        int intersection = kmers1.Count(kmers2.Contains);
        int union = kmers1.Count + kmers2.Count - intersection;

        // Union is empty only when both k-mer sets are empty; Jaccard is undefined there.
        return union == 0 ? 0.0 : (double)intersection / union * PercentScale;
    }

    #endregion

    #region Open Reading Frames

    // Standard genetic code (NCBI transl_table=1): start codon ATG, stop codons TAA/TAG/TGA. The stop codons
    // come from the canonical GeneticCode table (no private copy); ATG is the explicit Rosalind / ORFfinder
    // "ATG only" start policy.
    // NCBI Genetic Codes, https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (accessed 2026-06-14).
    private const string StartCodon = "ATG";

    // A codon is three nucleotides.
    private const int CodonLength = 3;

    /// <summary>
    /// Finds potential Open Reading Frames (ORFs) in all six reading frames (three on the
    /// forward strand and three on the reverse complement).
    /// </summary>
    /// <remarks>
    /// <para>
    /// An ORF starts at a start codon (ATG) and ends at the first in-frame stop codon
    /// (TAA, TAG, or TGA; the canonical <see cref="GeneticCode.Standard"/> stop set). Following the
    /// Rosalind "Open Reading Frames" problem (and NCBI <c>COrf::FindOrfs</c> with
    /// <c>longest_orfs = false</c>), <b>every</b> in-frame ATG that is terminated by a downstream in-frame
    /// stop is reported — including nested/overlapping ORFs that share the same stop codon. A reading begun at
    /// an ATG that has no downstream in-frame stop is not a complete ORF and is not reported. The reported
    /// <see cref="OrfInfo.Sequence"/> spans the start codon through the stop codon inclusive, so its length is
    /// divisible by 3; the encoded protein candidate is the translation up to (excluding) the stop.
    /// </para>
    /// <para>
    /// Coordinates: <see cref="OrfInfo.Position"/> is the 0-based start in the <i>scanned strand</i> — the input
    /// for forward ORFs, its reverse complement for <see cref="OrfInfo.IsReverseComplement"/> ORFs.
    /// <see cref="OrfInfo.GetForwardStrandRange(int)"/> converts to input (forward-strand) coordinates.
    /// </para>
    /// <para>
    /// <paramref name="minLength"/> is the minimum length in nucleotides <i>including the stop codon</i> — the
    /// <see cref="OrfInfo.Length"/> of the reported ORF (inclusive lower bound; values ≤ 0 admit every ORF). NCBI
    /// ORFfinder (<c>COrf</c>: "length in base pairs not including the stop") and EMBOSS getorf <c>-minsize</c>
    /// do not count the stop: use <see cref="OrfSearchOptions.MinLengthExcludesStop"/> for that convention.
    /// </para>
    /// <para>
    /// Complexity: O(n + output) — one pass per frame; pending starts are released at each stop. Use
    /// <see cref="FindOpenReadingFrames(DnaSequence, int, OrfSearchOptions)"/> for other genetic codes, start-codon
    /// sets, NCBI ORFfinder's longest-ORF-per-stop reporting and open-ended (partial) ORFs. The input is a
    /// <see cref="DnaSequence"/> (A/C/G/T); for IUPAC-ambiguous input use
    /// <see cref="Translator.FindOrfs(string, GeneticCode?, int, bool)"/> / <see cref="OrfScanner"/>.
    /// </para>
    /// </remarks>
    /// <param name="sequence">The DNA sequence to scan.</param>
    /// <param name="minLength">Minimum ORF length in nucleotides (stop codon included, inclusive lower bound).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static IEnumerable<OrfInfo> FindOpenReadingFrames(DnaSequence sequence, int minLength = 100)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindOpenReadingFramesCore(sequence, minLength, OrfSearchOptions.Default);
    }

    /// <summary>
    /// Finds ORFs with an explicit genetic code, start-codon policy and reporting convention
    /// (see <see cref="OrfSearchOptions"/>); with <see cref="OrfSearchOptions.Default"/> the result equals
    /// <see cref="FindOpenReadingFrames(DnaSequence, int)"/>.
    /// </summary>
    /// <remarks>
    /// With a non-default <see cref="OrfSearchOptions"/> this reproduces NCBI ORFfinder (<c>COrf::FindOrfs</c>:
    /// <c>-g</c> genetic code, <c>-s</c> start-codon mode, longest ORF per stop, partial ORFs, minimum length
    /// without the stop) and, for <see cref="OrfSearchOptions.NestedOrfs"/> = false, EMBOSS getorf <c>-find 3</c>.
    /// </remarks>
    /// <param name="sequence">The DNA sequence to scan.</param>
    /// <param name="minLength">Minimum ORF length in nucleotides (inclusive lower bound; the stop codon counts unless
    /// <see cref="OrfSearchOptions.MinLengthExcludesStop"/>).</param>
    /// <param name="options">The search options.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="options"/> is null.</exception>
    /// <exception cref="ArgumentException">The genetic code has dual-coding stop codons (NCBI tables 27, 28, 31).</exception>
    public static IEnumerable<OrfInfo> FindOpenReadingFrames(DnaSequence sequence, int minLength, OrfSearchOptions options)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(options);
        return FindOpenReadingFramesCore(sequence, minLength, options);
    }

    private static IEnumerable<OrfInfo> FindOpenReadingFramesCore(DnaSequence sequence, int minLength, OrfSearchOptions options)
    {
        var code = options.GeneticCode ?? GeneticCode.Standard;
        Func<string, bool>? isStart = options.StartCodons switch
        {
            OrfStartCodons.AtgOnly => static codon => codon == StartCodon,
            OrfStartCodons.GeneticCodeStarts => code.IsStartCodon,
            OrfStartCodons.AnySenseCodon => null,
            _ => throw new ArgumentOutOfRangeException(nameof(options), options.StartCodons, "Unknown start-codon mode."),
        };

        // Eager validation (genetic code with dual-coding stops) before the lazy enumeration starts.
        string forward = sequence.Sequence;
        var forwardRegions = OrfScanner.FindRegions(forward, code, isStart, options.NestedOrfs, options.IncludeOpenEnded);
        string? reverse = options.SearchBothStrands ? DnaSequence.GetReverseComplementString(forward) : null;
        var reverseRegions = reverse is null
            ? null
            : OrfScanner.FindRegions(reverse, code, isStart, options.NestedOrfs, options.IncludeOpenEnded);

        return EnumerateOrfs(forward, forwardRegions, reverse, reverseRegions, minLength, options.MinLengthExcludesStop);
    }

    private static IEnumerable<OrfInfo> EnumerateOrfs(string forward, IEnumerable<OrfRegion> forwardRegions,
        string? reverse, IEnumerable<OrfRegion>? reverseRegions, int minLength, bool excludeStop)
    {
        foreach (var orf in ToOrfInfos(forward, forwardRegions, isReverseComplement: false, minLength, excludeStop))
        {
            yield return orf;
        }

        if (reverse is null || reverseRegions is null)
        {
            yield break;
        }

        foreach (var orf in ToOrfInfos(reverse, reverseRegions, isReverseComplement: true, minLength, excludeStop))
        {
            yield return orf;
        }
    }

    private static IEnumerable<OrfInfo> ToOrfInfos(string strand, IEnumerable<OrfRegion> regions,
        bool isReverseComplement, int minLength, bool excludeStop)
    {
        foreach (var region in regions)
        {
            int length = region.Length;
            int measured = excludeStop && region.IsTerminated ? length - CodonLength : length;
            if (measured < minLength)
            {
                continue;
            }

            yield return new OrfInfo(
                strand.Substring(region.Start, length),
                region.Start,
                region.FrameOffset + 1,
                isReverseComplement,
                region.IsTerminated);
        }
    }

    #endregion

    #region Helper Methods

    private static HashSet<string> GetKmers(string sequence, int k)
    {
        var kmers = new HashSet<string>();
        for (int i = 0; i <= sequence.Length - k; i++)
        {
            kmers.Add(sequence.Substring(i, k));
        }
        return kmers;
    }

    #endregion
}

#region Result Types

/// <summary>
/// A known-motif hit (<see cref="GenomicAnalyzer.FindKnownMotifHits"/>): the 0-based start of the matching
/// window on the forward strand and whether the motif was found on the reverse strand (the window matches
/// the motif's reverse complement).
/// </summary>
public readonly record struct KnownMotifHit(int Position, bool IsReverseStrand);

/// <summary>
/// Information about a repeated region.
/// </summary>
public readonly struct RepeatInfo
{
    public static readonly RepeatInfo None = new(string.Empty, Array.Empty<int>());

    public RepeatInfo(string sequence, IReadOnlyList<int> positions)
    {
        Sequence = sequence;
        Positions = positions;
    }

    public string Sequence { get; }
    public IReadOnlyList<int> Positions { get; }
    public int Length => Sequence.Length;
    public int Count => Positions.Count;
    public bool IsEmpty => string.IsNullOrEmpty(Sequence);
}

/// <summary>
/// Information about a tandem repeat (consecutive repeating unit).
/// </summary>
public readonly struct TandemRepeat
{
    public TandemRepeat(string unit, int position, int repetitions)
    {
        Unit = unit;
        Position = position;
        Repetitions = repetitions;
    }

    public string Unit { get; }
    public int Position { get; }
    public int Repetitions { get; }
    public int TotalLength => Unit.Length * Repetitions;
    public string FullSequence => string.Concat(Enumerable.Repeat(Unit, Repetitions));
}

/// <summary>
/// Information about a palindromic sequence.
/// </summary>
public readonly struct PalindromeInfo
{
    public PalindromeInfo(string sequence, int position)
    {
        Sequence = sequence;
        Position = position;
    }

    public string Sequence { get; }
    public int Position { get; }
    public int Length => Sequence.Length;
}

/// <summary>
/// Information about a common region between two sequences.
/// </summary>
public readonly struct CommonRegion
{
    public static readonly CommonRegion None = new(string.Empty, -1, -1);

    public CommonRegion(string sequence, int positionInFirst, int positionInSecond)
    {
        Sequence = sequence;
        PositionInFirst = positionInFirst;
        PositionInSecond = positionInSecond;
    }

    public string Sequence { get; }
    public int PositionInFirst { get; }
    public int PositionInSecond { get; }
    public int Length => Sequence.Length;
    public bool IsEmpty => string.IsNullOrEmpty(Sequence);
}

/// <summary>
/// Information about an Open Reading Frame found by <see cref="GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int)"/>.
/// </summary>
public readonly struct OrfInfo
{
    public OrfInfo(string sequence, int position, int frame, bool isReverseComplement, bool isTerminated = true)
    {
        Sequence = sequence;
        Position = position;
        Frame = frame;
        IsReverseComplement = isReverseComplement;
        IsTerminated = isTerminated;
    }

    /// <summary>The ORF's nucleotides, start codon through the stop codon inclusive (through the last complete codon for an open-ended ORF), 5'→3' on the scanned strand.</summary>
    public string Sequence { get; }

    /// <summary>0-based start of the ORF in the <i>scanned strand</i>: the input sequence, or its reverse complement when <see cref="IsReverseComplement"/>. See <see cref="GetForwardStrandRange(int)"/> for input coordinates.</summary>
    public int Position { get; }

    /// <summary>Reading frame 1–3 within the scanned strand (offset <c>Position % 3</c> + 1).</summary>
    public int Frame { get; }

    /// <summary>True when the ORF was found on the reverse complement of the input.</summary>
    public bool IsReverseComplement { get; }

    /// <summary>True when the ORF ends with a stop codon; false for an open-ended (3'-partial) ORF that reaches the strand end (only reported with <see cref="OrfSearchOptions.IncludeOpenEnded"/>).</summary>
    public bool IsTerminated { get; }

    public int Length => Sequence.Length;
    public int CodonCount => Sequence.Length / 3;

    /// <summary>
    /// The ORF's location in input (forward-strand) coordinates: 0-based inclusive <c>(Start, End)</c> with
    /// <c>Start ≤ End</c>, the stop codon included. For a reverse-complement ORF the span is the one whose reverse
    /// complement is <see cref="Sequence"/>: <c>Start = L − Position − Length</c>, <c>End = L − 1 − Position</c>.
    /// </summary>
    /// <param name="sequenceLength">Length L of the searched input sequence.</param>
    public (int Start, int End) GetForwardStrandRange(int sequenceLength)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(sequenceLength, Position + Length);
        return IsReverseComplement
            ? (sequenceLength - Position - Length, sequenceLength - 1 - Position)
            : (Position, Position + Length - 1);
    }
}

/// <summary>Which codons open an ORF in <see cref="GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int, OrfSearchOptions)"/> (NCBI ORFfinder start-codon modes).</summary>
public enum OrfStartCodons
{
    /// <summary>ATG only (Rosalind; ORFfinder default).</summary>
    AtgOnly,

    /// <summary>ATG and the alternative initiation codons of the genetic code: the codons its NCBI <c>sncbieaa</c> string marks as starts (<see cref="GeneticCode.StartCodons"/>, e.g. TTG/CTG/ATG in table 1).</summary>
    GeneticCodeStarts,

    /// <summary>Any sense (non-stop) codon: stop-to-stop regions (ORFfinder "any sense codon").</summary>
    AnySenseCodon,
}

/// <summary>
/// Options of <see cref="GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int, OrfSearchOptions)"/>. The defaults
/// reproduce <see cref="GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int)"/> (Rosalind semantics).
/// </summary>
public sealed record OrfSearchOptions
{
    /// <summary>The default options (Standard code, ATG only, every in-frame start, terminated ORFs only, both strands).</summary>
    public static OrfSearchOptions Default { get; } = new();

    /// <summary>Genetic code whose stop codons end ORFs (and whose starts are used by <see cref="OrfStartCodons.GeneticCodeStarts"/>); <see langword="null"/> = <see cref="GeneticCode.Standard"/>. Tables with dual-coding stop codons (27, 28, 31) are rejected.</summary>
    public GeneticCode? GeneticCode { get; init; }

    /// <summary>Start-codon policy (default <see cref="OrfStartCodons.AtgOnly"/>).</summary>
    public OrfStartCodons StartCodons { get; init; } = OrfStartCodons.AtgOnly;

    /// <summary>
    /// true (default): every in-frame start that reaches a stop is reported, nested ORFs included (Rosalind; NCBI
    /// <c>COrf</c> <c>longest_orfs = false</c>). false: only the most upstream start per stop (NCBI ORFfinder default,
    /// EMBOSS getorf <c>-find 1/3</c>).
    /// </summary>
    public bool NestedOrfs { get; init; } = true;

    /// <summary>true: also report ORFs that run to the sequence end without a stop (3'-partial, <see cref="OrfInfo.IsTerminated"/> = false), ending at the last complete codon. Default false.</summary>
    public bool IncludeOpenEnded { get; init; }

    /// <summary>true: <c>minLength</c> is compared with the ORF length <i>without</i> its stop codon (NCBI ORFfinder <c>-ml</c>, EMBOSS getorf <c>-minsize</c>); default false (<see cref="OrfInfo.Length"/> ≥ minLength).</summary>
    public bool MinLengthExcludesStop { get; init; }

    /// <summary>Scan the reverse complement as well (default true).</summary>
    public bool SearchBothStrands { get; init; } = true;
}

#endregion
