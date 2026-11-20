using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Genomic analysis algorithms based on suffix trees.
/// Provides efficient solutions for common bioinformatics problems.
/// </summary>
public static partial class GenomicAnalyzer
{
    #region Repeat Finding

    /// <summary>
    /// Finds the longest repeated substring (LRS) of a DNA sequence: the longest
    /// substring that occurs at least twice (occurrences may overlap). This is the
    /// classical suffix-tree application — the deepest internal node, whose path label
    /// spells a longest substring occurring ≥ 2 times (CMU 15-451 Lecture #10 §2.1;
    /// Wikipedia, Longest repeated substring problem). Returns <see cref="RepeatInfo.None"/>
    /// when no substring repeats (including the empty sequence).
    /// <para>
    /// Several different substrings can tie for the maximum length; the sources require only "a" longest
    /// one. This method is deterministic: it returns the tied substring whose first occurrence is
    /// leftmost (the first entry of <see cref="FindAllLongestRepeats(DnaSequence)"/>), independent of the tree's
    /// internal child order.
    /// </para>
    /// Time complexity: O(n) suffix-tree construction (Ukkonen, cached on the sequence) plus one O(n)
    /// traversal for the deepest nodes and O(occ) for their occurrences.
    /// See docs/algorithms/Repeat_Analysis/Repeat_Detection.md.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static RepeatInfo FindLongestRepeat(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var all = FindAllLongestRepeats(sequence);
        return all.Count == 0 ? RepeatInfo.None : all[0];
    }

    /// <summary>
    /// Finds every longest repeated substring (all length ties) with all their start positions — the path
    /// labels of all internal suffix-tree nodes of maximal string depth (Gusfield 1997 §7.1). Overlapping
    /// occurrences count. Results are ordered by first occurrence (ascending), positions ascending. Empty
    /// when nothing repeats (including the empty sequence). Delegates to the canonical
    /// <c>SuffixTree.FindAllLongestRepeatedSubstrings</c>.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IReadOnlyList<RepeatInfo> FindAllLongestRepeats(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return sequence.SuffixTree.FindAllLongestRepeatedSubstrings()
            .Select(r => new RepeatInfo(r.Substring, r.Positions))
            .ToList();
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
        => EnumerateRepeats(sequence, minLength, maximalOnly: false, default);

    /// <summary>
    /// Finds every <b>maximal repeat</b> of length ≥ <paramref name="minLength"/>: a substring α that occurs at
    /// least twice and is both <i>left-diverse</i> (its occurrences are not all preceded by the same character;
    /// the start of the text counts as its own symbol) and <i>right-diverse</i> (not all followed by the same
    /// character; the end of the text likewise) — Gusfield 1997 §7.12 ("maximal pair" / "maximal repeat"), the
    /// strings behind REPuter / MUMmer <c>-maxmatch</c> self matches. Equivalently: the internal nodes of the
    /// suffix tree whose occurrence set is left-diverse (every internal node is right-diverse because of the
    /// terminator), so there are at most n of them. Each result carries every occurrence (overlaps allowed,
    /// positions ascending); results are in the same order as <see cref="FindRepeats(DnaSequence, int)"/> (smallest suffix that
    /// starts with the repeat). This is the subset of <see cref="FindRepeats(DnaSequence, int)"/> that cannot be extended on
    /// either side without losing an occurrence; <c>ACGTACGT</c> has the single maximal repeat <c>ACGT</c> @ {0,4}
    /// whereas <see cref="FindRepeats(DnaSequence, int)"/> also lists its prefixes and suffixes. Left-diversity is computed
    /// bottom-up in one linear traversal (a bit set of preceding characters per node), so the cost is
    /// O(n) plus the size of the output. For maximal repeated <i>pairs</i> (each pair of occurrences, REPuter /
    /// <c>repeat-match</c>) use <see cref="RepeatFinder.FindDirectRepeats(DnaSequence,int,int,int)"/>; for
    /// supermaximal repeats <see cref="RepeatFinder.FindSupermaximalRepeats(DnaSequence,int)"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1 (the empty string is not a repeat).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IEnumerable<RepeatInfo> FindMaximalRepeats(DnaSequence sequence, int minLength)
        => EnumerateRepeats(sequence, minLength, maximalOnly: true, default);

    /// <summary>
    /// <see cref="FindMaximalRepeats(DnaSequence, int)"/> that observes <paramref name="cancellationToken"/>: it is polled while the suffix tree is
    /// traversed (every 4096 nodes) and between the reported repeats; a cancelled token throws <see cref="OperationCanceledException"/> from the enumeration.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IEnumerable<RepeatInfo> FindMaximalRepeats(DnaSequence sequence, int minLength, CancellationToken cancellationToken)
        => EnumerateRepeats(sequence, minLength, maximalOnly: true, cancellationToken);

    /// <summary>
    /// <see cref="FindRepeats(DnaSequence, int)"/> that observes <paramref name="cancellationToken"/> (see
    /// <see cref="FindMaximalRepeats(DnaSequence, int, CancellationToken)"/>).
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IEnumerable<RepeatInfo> FindRepeats(DnaSequence sequence, int minLength, CancellationToken cancellationToken)
        => EnumerateRepeats(sequence, minLength, maximalOnly: false, cancellationToken);

    private static IEnumerable<RepeatInfo> EnumerateRepeats(DnaSequence sequence, int minLength, bool maximalOnly, CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return EnumerateRepeatsCore(sequence.SuffixTree, minLength, maximalOnly, cancellationToken);
    }

    // The engine behind every repeat enumeration: it works on any suffix tree whose text holds upper-case letters
    // (A-Z) and, optionally, separator characters above 'Z' (see MaskedText). A separator is never part of a repeat and
    // always ends one: a node whose path label starts with a separator is skipped (it can only be a run of the same
    // separator), and a separator before a suffix counts as a left context different from every other.
    private static IEnumerable<RepeatInfo> EnumerateRepeatsCore(global::SuffixTree.SuffixTree tree, int minLength, bool maximalOnly, CancellationToken cancellationToken)
    {
        AlignmentCancellation.Poll(cancellationToken);

        // A repeat is a NON-EMPTY substring occurring >=2 times (the path label of an internal
        // node, never the root). The empty string is not a repeat, so the effective minimum is
        // at least 1 even when the caller passes minLength <= 0. (CMU 15-451 Lecture #10 §2.1.)
        int effectiveMinLength = Math.Max(1, minLength);

        // Pre-order traversal with children in ascending key order (the terminator, a leaf, sorts first):
        // internal nodes come out in the order of their lexicographically smallest suffix, and a node
        // always precedes the internal nodes below it.
        var loci = new RepeatLocusCollector(effectiveMinLength, maximalOnly, tree.Text, cancellationToken);
        tree.Traverse(loci);

        foreach (var (edgeEnd, parentDepth, depth) in loci.Nodes)
        {
            AlignmentCancellation.Poll(cancellationToken);
            // Path label of the node: the `depth` characters of the text ending at the end of its incoming edge
            // (an internal node's edge never contains the terminator).
            string label = tree.Text.Substring(edgeEnd - depth, depth);
            var positions = new List<int>(tree.FindAllOccurrences(label));
            positions.Sort();

            // A maximal repeat is the full node label only; every shorter length on the edge is extendable to
            // the right without losing an occurrence.
            if (maximalOnly)
            {
                yield return new RepeatInfo(label, positions);
                continue;
            }

            for (int length = Math.Max(parentDepth + 1, effectiveMinLength); length < depth; length++)
                yield return new RepeatInfo(label.Substring(0, length), new List<int>(positions));

            yield return new RepeatInfo(label, positions);
        }
    }

    /// <summary>
    /// Collects the internal (non-root) nodes of a suffix tree whose path-label length is at least
    /// <c>minLength</c>, in traversal order, as (end of incoming edge in the text, parent path-label length, own path-label length).
    /// With <c>maximalOnly</c> only nodes whose occurrence set is left-diverse are kept: each leaf contributes the
    /// character before its suffix (a sentinel bit for suffix 0), a node's set is the union over its subtree, and
    /// the node qualifies when that set has ≥ 2 members (Gusfield 1997 §7.12.1).
    /// </summary>
    private sealed class RepeatLocusCollector(int minLength, bool maximalOnly, global::SuffixTree.ITextSource text, CancellationToken cancellationToken = default) : global::SuffixTree.ISuffixTreeVisitor
    {
        private bool _rootSeen;
        private int _visited;
        private readonly List<(int EdgeEnd, int ParentDepth, int Depth)> _nodes = new();
        private readonly List<bool> _leftDiverse = new();
        private readonly Stack<(int Mask, int Index)> _frames = new();

        public IEnumerable<(int EdgeEnd, int ParentDepth, int Depth)> Nodes
            => maximalOnly ? _nodes.Where((_, i) => _leftDiverse[i]) : _nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if ((++_visited & 0xFFF) == 0)
                AlignmentCancellation.Poll(cancellationToken);
            if (!_rootSeen)
            {
                _rootSeen = true; // the root has no incoming edge and is not a repeat
                if (maximalOnly)
                    _frames.Push((0, -1));
                return;
            }

            if (childCount == 0)
            {
                // A leaf is a single suffix, not a repeat. Its suffix starts at (edge start - parent depth); the character
                // before it is the left context (bit 0 = start of text; 'A','C','G','T' & 31 = 1,3,7,20, all distinct).
                if (maximalOnly)
                {
                    int p = startIndex - depth;
                    _frames.Push((p == 0 ? 1 : LeftContextBit(text[p - 1]), -1));
                }
                return;
            }

            int stringDepth = depth + (endIndex - startIndex);
            int index = -1;
            if (stringDepth >= minLength && text[endIndex - stringDepth] <= 'Z')
            {
                index = _nodes.Count;
                _nodes.Add((endIndex, depth, stringDepth));
                _leftDiverse.Add(!maximalOnly);
            }
            if (maximalOnly)
                _frames.Push((0, index));
        }

        // A letter A-Z & 31 is 1..26 (all distinct; bit 0 = start of text); a separator is its own left context, so it
        // sets two bits at once (any set containing it is already left-diverse).
        private static int LeftContextBit(char c) => c is >= 'A' and <= 'Z' ? 1 << (c & 31) : unchecked((int)0xC0000000);

        public void EnterBranch(int key) { }

        public void ExitBranch()
        {
            if (!maximalOnly)
                return;
            var (mask, index) = _frames.Pop();
            if (index >= 0 && (mask & (mask - 1)) != 0)
                _leftDiverse[index] = true;
            var (parentMask, parentIndex) = _frames.Pop();
            _frames.Push((parentMask | mask, parentIndex));
        }
    }

    /// <summary>
    /// Finds exact tandem repeats (consecutive copies of a unit such as "ATGATGATG"): for every unit
    /// length <c>p ≥ minUnitLength</c> each <b>maximal run</b> of period <c>p</c> with at least
    /// <paramref name="minRepetitions"/> complete copies and a <b>primitive</b> unit, reported once.
    /// <para>
    /// <b>Definition.</b> A tandem repeat is "two or more contiguous … copies of a pattern of nucleotides"
    /// (Benson 1999, Tandem Repeats Finder; this detector accepts exact copies only — approximate copies are
    /// <see cref="RepeatFinder.FindApproximateTandemRepeats(DnaSequence,int,int,int)"/>). For unit length <c>p</c>
    /// the repeat is a maximal interval <c>S[a..e)</c> with <c>S[x] = S[x+p]</c> for all <c>a ≤ x &lt; e−p</c> that
    /// cannot be extended left or right (a maximal repetition, Kolpakov &amp; Kucherov 1999; Gusfield 1997 §9),
    /// reported at its left end <c>a</c> with <see cref="TandemRepeat.Repetitions"/> = <c>⌊(e−a)/p⌋</c> complete
    /// copies (a trailing partial copy is not counted). The unit <c>S[a..a+p)</c> must be <b>primitive</b> (not a
    /// power of a shorter word): <c>ATATATAT</c> is AT×4 at 0 and not additionally ATAT×2 (microsatellite = a primitive
    /// unit repeated; MISA "reject false type motifs", Thiel et al. 2003). Rotations of the same run (TA inside
    /// <c>ATATATA</c>) are not reported again; runs of different unit lengths may overlap (e.g. in <c>AAAGAAAGAAAG</c>, with
    /// minUnitLength 1 and minRepetitions 2, AAAG×3 at 0 and A×3 at 0 are both runs) and are reported independently.
    /// </para>
    /// <para>
    /// <b>Order:</b> unit length ascending, then start ascending (deterministic). Case is irrelevant
    /// (a <see cref="DnaSequence"/> is upper-case A/C/G/T). This method used to be a separate O(n²·m) substring
    /// scan that also reported non-primitive units and had no upper bound on the unit length; it is now a facade
    /// over the canonical <see cref="RepeatFinder.FindMicrosatellites(DnaSequence,int,int,int)"/> (maximal runs,
    /// O(n) per unit length) and agrees with an independent brute-force oracle (Kolpakov–Kucherov runs; see
    /// docs/Evidence/GENOMIC-TANDEM-001-Evidence.md).
    /// </para>
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minUnitLength">
    /// Minimum repeat-unit length (period). Must be &gt;= 1: a zero-length unit is degenerate (the empty
    /// string "tiles" every position infinitely), so it is rejected rather than silently hanging.
    /// </param>
    /// <param name="minRepetitions">
    /// Minimum number of complete consecutive unit copies. Must be &gt;= 2: a tandem repeat requires at least two
    /// copies (<c>k &gt;= 2</c>).
    /// </param>
    /// <returns>The repeats, ordered by unit length then start; see the remarks above.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="minUnitLength"/> is &lt; 1, or <paramref name="minRepetitions"/> is &lt; 2.
    /// </exception>
    public static IEnumerable<TandemRepeat> FindTandemRepeats(DnaSequence sequence, int minUnitLength = 2, int minRepetitions = 2)
        => FindTandemRepeats(sequence, minUnitLength, minRepetitions, int.MaxValue);

    /// <summary>
    /// <see cref="FindTandemRepeats(DnaSequence,int,int)"/> restricted to unit lengths
    /// <c>[minUnitLength, maxUnitLength]</c> (e.g. 6 for microsatellites, 60 for minisatellites; TRF's period limit
    /// is a similar bound). A unit longer than <c>n / minRepetitions</c> cannot repeat and is never reported, so
    /// <see cref="int.MaxValue"/> (the 3-argument overload) means "no upper bound".
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minUnitLength">Minimum unit length, &gt;= 1.</param>
    /// <param name="minRepetitions">Minimum number of complete copies, &gt;= 2.</param>
    /// <param name="maxUnitLength">Maximum unit length, &gt;= <paramref name="minUnitLength"/>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="minUnitLength"/> &lt; 1, <paramref name="minRepetitions"/> &lt; 2 or
    /// <paramref name="maxUnitLength"/> &lt; <paramref name="minUnitLength"/>.
    /// </exception>
    public static IEnumerable<TandemRepeat> FindTandemRepeats(
        DnaSequence sequence, int minUnitLength, int minRepetitions, int maxUnitLength)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateTandemParameters(minUnitLength, minRepetitions, maxUnitLength);

        // The engine (one linear scan per unit length, or the run enumeration for a wide range such as the default
        // "no upper bound" on a megabase sequence) is chosen by the canonical RepeatFinder.FindMicrosatellites.
        return RepeatFinder.FindMicrosatellites(sequence, minUnitLength, maxUnitLength, minRepetitions)
            .Select(m => new TandemRepeat(m.RepeatUnit, m.Position, m.RepeatCount));
    }

    /// <summary>
    /// <see cref="FindTandemRepeats(DnaSequence,int,int,int)"/> with an optional <b>canonical rotation</b> of the reported unit.
    /// </summary>
    /// <remarks>
    /// A run of period p has p equivalent spellings of its unit (<c>ATG</c>, <c>TGA</c>, <c>GAT</c> for the same run);
    /// the detector reports the spelling that starts at the run's left end. With <paramref name="canonicalRotation"/> the
    /// reported <see cref="TandemRepeat.Unit"/> is instead the <b>lexicographically least rotation</b> of that unit (ordinal
    /// A &lt; C &lt; G &lt; T; Booth 1980, "Lexicographically least circular substrings"; the representative used by MISA's
    /// repeat-type classes, <see cref="RepeatFinder.GetCanonicalMotifClass(string)"/>), so that rotations of one motif group
    /// together. The run itself is unchanged: <see cref="TandemRepeat.Position"/> is still the left end of the run and
    /// <see cref="TandemRepeat.Repetitions"/> the number of complete copies, <see cref="TandemRepeat.Rotation"/> says how far
    /// the unit was rotated and <see cref="TandemRepeat.OriginalUnit"/> / <see cref="TandemRepeat.FullSequence"/> still spell the
    /// actual bases of the sequence. Units are primitive, so the least rotation is unique. With
    /// <c>canonicalRotation = false</c> the result is that of the 4-argument overload.
    /// </remarks>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minUnitLength">Minimum unit length, &gt;= 1.</param>
    /// <param name="minRepetitions">Minimum number of complete copies, &gt;= 2.</param>
    /// <param name="maxUnitLength">Maximum unit length, &gt;= <paramref name="minUnitLength"/>.</param>
    /// <param name="canonicalRotation">Report the least rotation of each unit.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid unit / copy bounds (see the 4-argument overload).</exception>
    public static IEnumerable<TandemRepeat> FindTandemRepeats(
        DnaSequence sequence, int minUnitLength, int minRepetitions, int maxUnitLength, bool canonicalRotation)
    {
        var repeats = FindTandemRepeats(sequence, minUnitLength, minRepetitions, maxUnitLength);
        return canonicalRotation ? repeats.Select(CanonicalizeRotation) : repeats;
    }

    /// <summary>
    /// <see cref="FindTandemRepeats(string,int,int,int)"/> with the canonical-rotation option of
    /// <see cref="FindTandemRepeats(DnaSequence,int,int,int,bool)"/>.
    /// </summary>
    /// <param name="sequence">Nucleotide string to analyze. Must not be <c>null</c>.</param>
    /// <param name="minUnitLength">Minimum unit length, &gt;= 1.</param>
    /// <param name="minRepetitions">Minimum number of complete copies, &gt;= 2.</param>
    /// <param name="maxUnitLength">Maximum unit length, &gt;= <paramref name="minUnitLength"/>.</param>
    /// <param name="canonicalRotation">Report the least rotation of each unit.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid unit / copy bounds.</exception>
    public static IEnumerable<TandemRepeat> FindTandemRepeats(
        string sequence, int minUnitLength, int minRepetitions, int maxUnitLength, bool canonicalRotation)
    {
        var repeats = FindTandemRepeats(sequence, minUnitLength, minRepetitions, maxUnitLength);
        return canonicalRotation ? repeats.Select(CanonicalizeRotation) : repeats;
    }

    // The unit of a run is primitive, so its least rotation is unique and sits at a single offset of unit+unit.
    private static TandemRepeat CanonicalizeRotation(TandemRepeat repeat)
    {
        string canonical = RepeatFinder.GetMinimalRotation(repeat.Unit);
        int rotation = (repeat.Unit + repeat.Unit).IndexOf(canonical, StringComparison.Ordinal);
        return new TandemRepeat(canonical, repeat.Position, repeat.Repetitions, rotation);
    }

    /// <summary>
    /// <see cref="FindTandemRepeats(DnaSequence,int,int)"/> on a nucleotide <b>string</b>, for sequences that
    /// <see cref="DnaSequence"/> (A/C/G/T only) cannot hold: a genome with N runs, IUPAC ambiguity codes, RNA or lower case.
    /// </summary>
    /// <remarks>
    /// The contract is that of the canonical <see cref="RepeatFinder.FindMicrosatellites(string,int,int,int)"/>: the input is
    /// case-insensitive and only units spelled with the unambiguous bases A, C, G, T can repeat (MISA <c>[acgt]</c>; pytrf skips
    /// <c>N</c>), so a window that contains <c>N</c>, another IUPAC code, <c>U</c>, a gap or any other symbol never forms a
    /// repeat — <c>NNNN</c> is not a mononucleotide repeat and <c>ACNACN</c> is not an ACN repeat — while an ordinary repeat
    /// next to an N run is reported exactly as without it (<c>ACACNACAC</c> → AC×2 at 0 and AC×2 at 5). Everything else
    /// (primitive units, one maximal run per unit length, ordering, coordinates, parameter checks) is as documented on the
    /// <see cref="DnaSequence"/> overload; for an A/C/G/T string the result is identical. A <c>null</c> string throws
    /// <see cref="ArgumentNullException"/>, an empty one yields no repeats.
    /// </remarks>
    /// <param name="sequence">Nucleotide string to analyze. Must not be <c>null</c>.</param>
    /// <param name="minUnitLength">Minimum unit length, &gt;= 1.</param>
    /// <param name="minRepetitions">Minimum number of complete copies, &gt;= 2.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minUnitLength"/> &lt; 1 or <paramref name="minRepetitions"/> &lt; 2.</exception>
    public static IEnumerable<TandemRepeat> FindTandemRepeats(string sequence, int minUnitLength = 2, int minRepetitions = 2)
        => FindTandemRepeats(sequence, minUnitLength, minRepetitions, int.MaxValue);

    /// <summary>
    /// <see cref="FindTandemRepeats(string,int,int)"/> restricted to unit lengths <c>[minUnitLength, maxUnitLength]</c>.
    /// </summary>
    /// <param name="sequence">Nucleotide string to analyze. Must not be <c>null</c>.</param>
    /// <param name="minUnitLength">Minimum unit length, &gt;= 1.</param>
    /// <param name="minRepetitions">Minimum number of complete copies, &gt;= 2.</param>
    /// <param name="maxUnitLength">Maximum unit length, &gt;= <paramref name="minUnitLength"/>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="minUnitLength"/> &lt; 1, <paramref name="minRepetitions"/> &lt; 2 or
    /// <paramref name="maxUnitLength"/> &lt; <paramref name="minUnitLength"/>.
    /// </exception>
    public static IEnumerable<TandemRepeat> FindTandemRepeats(
        string sequence, int minUnitLength, int minRepetitions, int maxUnitLength)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ValidateTandemParameters(minUnitLength, minRepetitions, maxUnitLength);

        return RepeatFinder.FindMicrosatellites(sequence, minUnitLength, maxUnitLength, minRepetitions)
            .Select(m => new TandemRepeat(m.RepeatUnit, m.Position, m.RepeatCount));
    }

    private static void ValidateTandemParameters(int minUnitLength, int minRepetitions, int maxUnitLength)
    {
        if (minUnitLength < 1)
            throw new ArgumentOutOfRangeException(nameof(minUnitLength), minUnitLength,
                "Minimum repeat-unit length must be at least 1; a zero-length unit never terminates the scan.");
        if (minRepetitions < 2)
            throw new ArgumentOutOfRangeException(nameof(minRepetitions), minRepetitions,
                "A tandem repeat requires at least two consecutive copies; minRepetitions must be >= 2.");
        if (maxUnitLength < minUnitLength)
            throw new ArgumentOutOfRangeException(nameof(maxUnitLength), maxUnitLength,
                "Maximum repeat-unit length must be at least minUnitLength.");
    }

    #endregion

    #region Motif Finding

    /// <summary>
    /// Finds all occurrences of an exact motif (pattern) in a sequence: the 0-based start of every
    /// window <c>S[i..i+m-1]</c> equal to the motif, overlapping occurrences included, in
    /// <b>ascending</b> order. Matching is case-insensitive and exact: an IUPAC ambiguity code in the
    /// motif is a literal symbol (a DNA sequence never contains one, so it matches nothing) and only
    /// the forward strand is searched. Use <see cref="FindKnownMotifHits(DnaSequence, IEnumerable{string}, bool, bool, int)"/> for degenerate (IUPAC)
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
    /// Finds palindromic sequences (restriction enzyme recognition sites): every window of length
    /// <paramref name="minLength"/>…<paramref name="maxLength"/> that equals its own reverse complement, i.e. reads
    /// the same 5'→3' on both strands. Example: GAATTC (EcoRI site) — complement CTTAAG, reversed GAATTC.
    /// Nested and overlapping windows are all reported (Rosalind REVP convention), ordered by position, then length;
    /// positions are 0-based.
    /// <para>
    /// A thin facade over the canonical <see cref="RepeatFinder.FindPalindromes(DnaSequence,int,int)"/> (REP-PALIN-001;
    /// this method used to be a verbatim copy of its scan without the validation), mapping each
    /// <see cref="PalindromeResult"/> to a <see cref="PalindromeInfo"/>. Contract: <paramref name="minLength"/> must be even and
    /// ≥ 4 and <paramref name="maxLength"/> ≥ <paramref name="minLength"/> (otherwise
    /// <see cref="ArgumentOutOfRangeException"/>, thrown at the call, not on enumeration); a <c>null</c> sequence throws
    /// <see cref="ArgumentNullException"/>. A reverse-complement palindrome always has even length, so an odd
    /// <paramref name="maxLength"/> is effectively rounded down. Cost O(n·maxLength) plus output.
    /// </para>
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="minLength"/> is odd or &lt; 4, or <paramref name="maxLength"/> &lt; <paramref name="minLength"/>.
    /// </exception>
    public static IEnumerable<PalindromeInfo> FindPalindromes(DnaSequence sequence, int minLength = 4, int maxLength = 12)
        => RepeatFinder.FindPalindromes(sequence, minLength, maxLength)
            .Select(p => new PalindromeInfo(p.Sequence, p.Position));

    /// <summary>
    /// <see cref="FindPalindromes(DnaSequence,int,int)"/> on a nucleotide <b>string</b> that may hold N, other IUPAC codes, RNA
    /// or lower case (<see cref="DnaSequence"/> is A/C/G/T only). The contract is that of the canonical
    /// <see cref="RepeatFinder.FindPalindromes(string,int,int)"/>: case-insensitive, and only the unambiguous bases A, C, G, T pair
    /// (canonical complement), so a window containing N, another IUPAC code, U, a gap or any other symbol is never a
    /// palindrome (<c>NNNN</c> and <c>SSSS</c> equal their symbolic reverse complement but need not be palindromic once
    /// resolved); an A/C/G/T string gives the same result as the <see cref="DnaSequence"/> overload. A <c>null</c> string throws
    /// <see cref="ArgumentNullException"/>, an empty one yields no palindromes; the length checks are those of the
    /// <see cref="DnaSequence"/> overload.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="minLength"/> is odd or &lt; 4, or <paramref name="maxLength"/> &lt; <paramref name="minLength"/>.
    /// </exception>
    public static IEnumerable<PalindromeInfo> FindPalindromes(string sequence, int minLength = 4, int maxLength = 12)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return RepeatFinder.FindPalindromes(sequence, minLength, maxLength)
            .Select(p => new PalindromeInfo(p.Sequence, p.Position));
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
    /// motifs are skipped (the empty string is not a motif; mirrors <see cref="FindMotif(DnaSequence, string)"/>), and a
    /// motif with no occurrence is omitted from the result. Positions for each motif are returned
    /// <b>sorted ascending</b>. Duplicate motifs that normalize to the same upper-cased key collapse to
    /// a single entry. Forward strand only, IUPAC symbols literal (see <see cref="FindKnownMotifHits(DnaSequence, IEnumerable{string}, bool, bool, int)"/>
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
    /// <para>
    /// <paramref name="maxMismatches"/> &gt; 0 is the mismatch-tolerant search of fuzznuc <c>-pmismatch k</c>: a window
    /// is a hit when at most k of its m positions disagree with the motif (Hamming distance; substitutions only, no
    /// gaps). Without <paramref name="allowIupac"/> every motif symbol is literal, so a degenerate symbol such as
    /// <c>N</c> costs one mismatch at each position; with it a position agrees when the base is in the code's base set
    /// (the distance is then the Hamming distance to the nearest concrete expansion of the motif). Each hit carries its
    /// own <see cref="KnownMotifHit.Mismatches"/> count; on the reverse strand the reverse-complement motif is
    /// matched with the same bound (a palindromic site can be one hit per strand with different counts). A bound of m
    /// or more makes every window a hit. The scan is the canonical <see cref="ApproximateMatcher"/> Hamming matcher
    /// (O(n·m) per query, no index); exact search (k = 0) keeps the KMP / suffix-tree routes.
    /// </para>
    /// Verified against EMBOSS fuzznuc 6.6.0 (<c>-pmismatch 0..3</c>, with and without <c>-complement</c>).
    /// </summary>
    /// <param name="sequence">Subject sequence.</param>
    /// <param name="motifs">Motifs to search.</param>
    /// <param name="bothStrands">Also report occurrences on the reverse strand.</param>
    /// <param name="allowIupac">Interpret IUPAC ambiguity codes in the motifs.</param>
    /// <param name="maxMismatches">Maximum number of mismatching positions per hit (default 0 = exact).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static Dictionary<string, IReadOnlyList<KnownMotifHit>> FindKnownMotifHits(
        DnaSequence sequence,
        IEnumerable<string> motifs,
        bool bothStrands = false,
        bool allowIupac = false,
        int maxMismatches = 0)
        => FindKnownMotifHits(sequence, motifs, bothStrands, allowIupac, maxMismatches, CancellationToken.None);

    /// <summary>
    /// <see cref="FindKnownMotifHits(DnaSequence, IEnumerable{string}, bool, bool, int)"/> that observes <paramref name="cancellationToken"/>: it is polled
    /// before every (motif, strand) search and inside the mismatch-tolerant Hamming scan; a cancelled token throws <see cref="OperationCanceledException"/>
    /// and no partial dictionary is returned.
    /// </summary>
    /// <param name="sequence">Subject sequence.</param>
    /// <param name="motifs">Motifs to search.</param>
    /// <param name="bothStrands">Also report occurrences on the reverse strand.</param>
    /// <param name="allowIupac">Interpret IUPAC ambiguity codes in the motifs.</param>
    /// <param name="maxMismatches">Maximum number of mismatching positions per hit.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxMismatches"/> is negative.</exception>
    public static Dictionary<string, IReadOnlyList<KnownMotifHit>> FindKnownMotifHits(
        DnaSequence sequence,
        IEnumerable<string> motifs,
        bool bothStrands,
        bool allowIupac,
        int maxMismatches,
        CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(motifs);
        ArgumentOutOfRangeException.ThrowIfNegative(maxMismatches);
        AlignmentCancellation.Poll(cancellationToken);

        if (maxMismatches > 0)
            return FindKnownMotifHitsApproximate(sequence, motifs, bothStrands, allowIupac, maxMismatches, cancellationToken);

        var found = SearchKnownMotifs(sequence, motifs, bothStrands, allowIupac, cancellationToken);
        var result = new Dictionary<string, IReadOnlyList<KnownMotifHit>>(found.Count);
        foreach (var (key, (forward, reverse)) in found)
            result[key] = MergeStrandHits(forward, reverse);
        return result;
    }

    // Mismatch-tolerant route: one canonical ApproximateMatcher Hamming scan per (motif, strand).
    private static Dictionary<string, IReadOnlyList<KnownMotifHit>> FindKnownMotifHitsApproximate(
        DnaSequence sequence, IEnumerable<string> motifs, bool bothStrands, bool allowIupac, int maxMismatches, CancellationToken cancellationToken)
        => FindKnownMotifHitsApproximate(sequence.Sequence, motifs, bothStrands, allowIupac, maxMismatches, cancellationToken);

    private static Dictionary<string, IReadOnlyList<KnownMotifHit>> FindKnownMotifHitsApproximate(
        string text, IEnumerable<string> motifs, bool bothStrands, bool allowIupac, int maxMismatches, CancellationToken cancellationToken)
    {
        var result = new Dictionary<string, IReadOnlyList<KnownMotifHit>>();
        foreach (string key in NormalizeMotifs(motifs))
        {
            AlignmentCancellation.Poll(cancellationToken);
            var forward = ApproximateMatcher.FindWithMismatches(text, key, maxMismatches, allowIupac, cancellationToken)
                .Select(m => new KnownMotifHit(m.Position, IsReverseStrand: false, m.Distance)).ToList();
            var reverse = bothStrands
                ? ApproximateMatcher.FindWithMismatches(text, DnaSequence.GetReverseComplementString(key), maxMismatches, allowIupac, cancellationToken)
                    .Select(m => new KnownMotifHit(m.Position, IsReverseStrand: true, m.Distance)).ToList()
                : new List<KnownMotifHit>();
            if (forward.Count > 0 || reverse.Count > 0)
                result[key] = MergeStrandHits(forward, reverse);
        }

        return result;
    }

    // Distinct upper-cased non-blank motifs in input order (the empty string is not a motif).
    private static List<string> NormalizeMotifs(IEnumerable<string> motifs)
    {
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

        return keys;
    }

    private static Dictionary<string, (List<int> Forward, List<int> Reverse)> SearchKnownMotifs(
        DnaSequence sequence, IEnumerable<string> motifs, bool bothStrands, bool allowIupac, CancellationToken cancellationToken = default)
        => SearchKnownMotifs(sequence.Sequence, sequence, motifs, bothStrands, allowIupac, cancellationToken);

    // `indexed` is the sequence whose cached suffix tree serves the many-motif route; null for a plain string (linear scans only).
    private static Dictionary<string, (List<int> Forward, List<int> Reverse)> SearchKnownMotifs(
        string text, DnaSequence? indexed, IEnumerable<string> motifs, bool bothStrands, bool allowIupac, CancellationToken cancellationToken = default)
    {
        // The empty string is not a motif (an exact search for "" would match every position).
        // Skip empty/whitespace motifs, consistent with FindMotif's empty-pattern guard.
        var keys = NormalizeMotifs(motifs);

        bool useIndex = indexed is not null && !allowIupac && keys.Count * (bothStrands ? 2L : 1L) >= IndexedSearchMinQueries;
        var result = new Dictionary<string, (List<int>, List<int>)>(keys.Count);
        foreach (string key in keys)
        {
            AlignmentCancellation.Poll(cancellationToken);
            var forward = FindMotifPositions(text, indexed, key, allowIupac, useIndex, cancellationToken);
            var reverse = bothStrands
                ? FindMotifPositions(text, indexed, DnaSequence.GetReverseComplementString(key), allowIupac, useIndex, cancellationToken)
                : new List<int>();
            if (forward.Count > 0 || reverse.Count > 0)
                result[key] = (forward, reverse);
        }

        return result;
    }

    private static List<int> FindMotifPositions(string text, DnaSequence? indexed, string motif, bool allowIupac, bool useIndex, CancellationToken cancellationToken)
    {
        if (allowIupac)
        {
            // MotifFinder scans left to right: ascending, overlapping.
            var matches = indexed is not null
                ? MotifFinder.FindDegenerateMotif(indexed, motif) // validates the motif even for an empty sequence
                : MotifFinder.FindDegenerateMotif(text, motif, cancellationToken);
            return matches.Select(m => m.Position).ToList();
        }

        if (useIndex)
        {
            // Suffix-tree occurrence enumeration is DFS-order (not sorted).
            var positions = indexed!.SuffixTree.FindAllOccurrences(motif);
            var sorted = new List<int>(positions);
            sorted.Sort();
            return sorted;
        }

        return KmerAnalyzer.FindKmerPositions(text, motif).ToList();
    }

    private static List<KnownMotifHit> MergeStrandHits(List<int> forward, List<int> reverse)
        => MergeStrandHits(
            forward.ConvertAll(p => new KnownMotifHit(p, IsReverseStrand: false)),
            reverse.ConvertAll(p => new KnownMotifHit(p, IsReverseStrand: true)));

    // Both lists ascending by position; a forward hit precedes a reverse hit at the same position.
    private static List<KnownMotifHit> MergeStrandHits(List<KnownMotifHit> forward, List<KnownMotifHit> reverse)
    {
        var merged = new List<KnownMotifHit>(forward.Count + reverse.Count);
        int f = 0, r = 0;
        while (f < forward.Count || r < reverse.Count)
        {
            if (r >= reverse.Count || (f < forward.Count && forward[f].Position <= reverse[r].Position))
                merged.Add(forward[f++]);
            else
                merged.Add(reverse[r++]);
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
        => FindCommonRegions(sequence1, sequence2, minLength, CancellationToken.None);

    /// <summary>
    /// <see cref="FindCommonRegions(DnaSequence, DnaSequence, int)"/> that observes <paramref name="cancellationToken"/>: it is polled every 1024 positions of the
    /// matching-statistics scan (which runs eagerly, before the first region is returned) and between the reported regions; a cancelled token throws
    /// <see cref="OperationCanceledException"/>.
    /// </summary>
    /// <param name="sequence1">Text sequence (indexed by its suffix tree).</param>
    /// <param name="sequence2">Query sequence.</param>
    /// <param name="minLength">Minimum region length; values below 1 are treated as 1.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence1"/> or <paramref name="sequence2"/> is null.</exception>
    public static IEnumerable<CommonRegion> FindCommonRegions(
        DnaSequence sequence1,
        DnaSequence sequence2,
        int minLength,
        CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        AlignmentCancellation.Poll(cancellationToken);

        // A common region is a NON-EMPTY contiguous substring; the empty string is not a region,
        // so the effective minimum length is at least 1 (Wikipedia "Longest common substring").
        int effectiveMinLength = Math.Max(1, minLength);
        int[] matchLengths = MatchingStatistics(sequence1.SuffixTree, sequence2.Sequence, sequence1.Length, effectiveMinLength, cancellationToken);
        return EnumerateCommonRegions(sequence1.SuffixTree, sequence2.Sequence, matchLengths, effectiveMinLength, cancellationToken);
    }

    private static IEnumerable<CommonRegion> EnumerateCommonRegions(
        global::SuffixTree.SuffixTree tree, string seq2, int[] matchLengths, int effectiveMinLength, CancellationToken cancellationToken)
    {
        var found = new HashSet<string>(StringComparer.Ordinal);
        var foundLookup = found.GetAlternateLookup<ReadOnlySpan<char>>();

        for (int i = 0; i < matchLengths.Length; i++)
        {
            if ((i & 0x3FF) == 0)
                AlignmentCancellation.Poll(cancellationToken);
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
    private static int[] MatchingStatistics(global::SuffixTree.SuffixTree tree, string query, int textLength, int minLength, CancellationToken cancellationToken = default)
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
                if ((i & 0x3FF) == 0)
                    AlignmentCancellation.Poll(cancellationToken);
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
            if ((i & 0x3FF) == 0)
                AlignmentCancellation.Poll(cancellationToken);
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
    /// This overload compares the <b>literal</b> (strand-specific) k-mers, i.e. a Python
    /// <c>set</c> of the windows / scipy <c>jaccard</c> on presence vectors: a sequence and its reverse
    /// complement are generally <i>not</i> similar. Use
    /// <see cref="CalculateSimilarity(DnaSequence, DnaSequence, int, bool)"/> with <c>canonical: true</c> for
    /// the strand-agnostic index that Mash and sourmash estimate.
    /// </para>
    /// <para>
    /// The Jaccard index is undefined when the union is empty (Jaccard 1901 is stated for
    /// non-empty sets). This occurs only when both k-mer sets are empty (e.g. both sequences
    /// are empty or shorter than <paramref name="kmerSize"/>); this method returns 0 in that case.
    /// The computation is delegated to <see cref="KmerAnalyzer.JaccardSimilarity(string, string, int)"/>.
    /// </para>
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="kmerSize">k-mer length; must be ≥ 1.</param>
    /// <returns>Jaccard similarity as a percentage in [0, 100].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence1"/> or <paramref name="sequence2"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> is less than 1.</exception>
    public static double CalculateSimilarity(DnaSequence sequence1, DnaSequence sequence2, int kmerSize = DefaultKmerSize)
        => CalculateSimilarity(sequence1, sequence2, kmerSize, canonical: false);

    /// <summary>
    /// Computes the k-mer Jaccard similarity index (percentage in [0, 100]) with an explicit strand mode.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <paramref name="canonical"/> = <c>false</c>: literal k-mers (see
    /// <see cref="CalculateSimilarity(DnaSequence, DnaSequence, int)"/>).
    /// <paramref name="canonical"/> = <c>true</c>: every k-mer is replaced by min(w, reverse-complement(w)),
    /// the strand-agnostic k-mer set of Mash (<c>mash sketch</c> default, Ondov et al. 2016) and sourmash;
    /// the result is the exact (unsketched) Jaccard index those tools estimate, so a sequence and its reverse
    /// complement score 100. <see cref="DnaSequence"/> is upper-case ACGT, so no ambiguous-base windows exist
    /// and the Mash <c>-n</c> / Jellyfish ACGT-only rule is satisfied by construction.
    /// </para>
    /// </remarks>
    /// <param name="sequence1">First DNA sequence.</param>
    /// <param name="sequence2">Second DNA sequence.</param>
    /// <param name="kmerSize">k-mer length; must be ≥ 1.</param>
    /// <param name="canonical">Use canonical (strand-agnostic) k-mers.</param>
    /// <returns>Jaccard similarity as a percentage in [0, 100]; 0 when both k-mer sets are empty.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence1"/> or <paramref name="sequence2"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="kmerSize"/> is less than 1.</exception>
    public static double CalculateSimilarity(DnaSequence sequence1, DnaSequence sequence2, int kmerSize, bool canonical)
    {
        ArgumentNullException.ThrowIfNull(sequence1);

        ArgumentNullException.ThrowIfNull(sequence2);

        if (kmerSize < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(kmerSize), kmerSize, "k-mer size must be at least 1.");
        }

        // Canonical k-mer Jaccard (DUP_MAP §16); the fraction is only re-expressed as a percentage here.
        return KmerAnalyzer.JaccardSimilarity(
            sequence1.Sequence, sequence2.Sequence, kmerSize, new KmerCountingOptions(Canonical: canonical)) * PercentScale;
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
        return FindOpenReadingFramesCore(sequence, minLength, OrfSearchOptions.Default, CancellationToken.None);
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
        return FindOpenReadingFramesCore(sequence, minLength, options, CancellationToken.None);
    }

    /// <summary>
    /// <see cref="FindOpenReadingFrames(DnaSequence, int, OrfSearchOptions)"/> that observes <paramref name="cancellationToken"/>: it is polled while the
    /// (lazy) frame scan advances (every 65 536 codons) and between the reported ORFs; a cancelled token throws <see cref="OperationCanceledException"/>
    /// from the enumeration.
    /// </summary>
    /// <param name="sequence">The DNA sequence to scan.</param>
    /// <param name="minLength">Minimum ORF length in nucleotides.</param>
    /// <param name="options">The search options.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="options"/> is null.</exception>
    public static IEnumerable<OrfInfo> FindOpenReadingFrames(DnaSequence sequence, int minLength, OrfSearchOptions options, CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(options);
        AlignmentCancellation.Poll(cancellationToken);
        return FindOpenReadingFramesCore(sequence, minLength, options, cancellationToken);
    }

    /// <summary>
    /// <see cref="FindOpenReadingFrames(DnaSequence, int)"/> on a nucleotide <b>string</b> that may hold N, other IUPAC codes, RNA
    /// or lower case — e.g. a genome with N runs, which <see cref="DnaSequence"/> (A/C/G/T only) cannot hold.
    /// </summary>
    /// <remarks>See <see cref="FindOpenReadingFrames(string, int, OrfSearchOptions)"/> for the ambiguity rules.</remarks>
    /// <param name="sequence">IUPAC nucleotide string (DNA or RNA, case-insensitive).</param>
    /// <param name="minLength">Minimum ORF length in nucleotides (stop codon included, inclusive lower bound).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="sequence"/> contains a symbol that is not an IUPAC nucleotide code.</exception>
    public static IEnumerable<OrfInfo> FindOpenReadingFrames(string sequence, int minLength = 100)
        => FindOpenReadingFrames(sequence, minLength, OrfSearchOptions.Default);

    /// <summary>
    /// <see cref="FindOpenReadingFrames(DnaSequence, int, OrfSearchOptions)"/> on a nucleotide <b>string</b> that may hold N, other
    /// IUPAC codes, RNA or lower case.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The scan is the shared <see cref="OrfScanner"/> core, whose ambiguity rules are those of EMBOSS getorf (verified with the
    /// 6.6.0 binary) and NCBI <c>COrf</c>: a codon containing an IUPAC ambiguity symbol is a <b>stop</b> only when <i>every</i>
    /// concrete codon it stands for is a stop of the genetic code (<c>TAR</c> and <c>TRA</c> in table 1), it is <b>never a start</b>
    /// (only unambiguous A/C/G/T codons open an ORF, also in <see cref="OrfStartCodons.AnySenseCodon"/> mode) and otherwise it is
    /// an ordinary sense codon the ORF reads through. So <c>NNN</c> neither opens nor ends an ORF, while <c>ATG</c> … <c>NNN</c> …
    /// <c>TAA</c> is one ORF. Input is case-insensitive; <c>U</c> is read as <c>T</c>; the reported
    /// <see cref="OrfInfo.Sequence"/> is the upper-case DNA spelling, the reverse strand is the canonical IUPAC reverse complement
    /// (<see cref="DnaSequence.GetReverseComplementString(string)"/>: R↔Y, K↔M, N↔N, …), and
    /// <see cref="OrfInfo.GetForwardStrandRange(int)"/> takes the length of the input string. For an A/C/G/T string the result
    /// equals the <see cref="DnaSequence"/> overload. An empty string yields no ORFs.
    /// </para>
    /// </remarks>
    /// <param name="sequence">IUPAC nucleotide string (DNA or RNA, case-insensitive).</param>
    /// <param name="minLength">Minimum ORF length in nucleotides (inclusive lower bound; the stop codon counts unless
    /// <see cref="OrfSearchOptions.MinLengthExcludesStop"/>).</param>
    /// <param name="options">The search options.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="options"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="sequence"/> contains a symbol that is not an IUPAC nucleotide code, or the
    /// genetic code has dual-coding stop codons (NCBI tables 27, 28, 31).</exception>
    public static IEnumerable<OrfInfo> FindOpenReadingFrames(string sequence, int minLength, OrfSearchOptions options)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(options);
        return FindOpenReadingFramesCore(SequenceExtensions.NormalizeIupacDna(sequence), minLength, options, CancellationToken.None);
    }

    /// <summary>
    /// <see cref="FindOpenReadingFrames(string, int, OrfSearchOptions)"/> that observes <paramref name="cancellationToken"/> (see
    /// <see cref="FindOpenReadingFrames(DnaSequence, int, OrfSearchOptions, CancellationToken)"/>).
    /// </summary>
    /// <param name="sequence">IUPAC nucleotide string (DNA or RNA, case-insensitive).</param>
    /// <param name="minLength">Minimum ORF length in nucleotides.</param>
    /// <param name="options">The search options.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="options"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="sequence"/> contains a symbol that is not an IUPAC nucleotide code.</exception>
    public static IEnumerable<OrfInfo> FindOpenReadingFrames(string sequence, int minLength, OrfSearchOptions options, CancellationToken cancellationToken)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(options);
        AlignmentCancellation.Poll(cancellationToken);
        return FindOpenReadingFramesCore(SequenceExtensions.NormalizeIupacDna(sequence), minLength, options, cancellationToken);
    }

    private static IEnumerable<OrfInfo> FindOpenReadingFramesCore(DnaSequence sequence, int minLength, OrfSearchOptions options, CancellationToken cancellationToken)
        => FindOpenReadingFramesCore(sequence.Sequence, minLength, options, cancellationToken);

    private static IEnumerable<OrfInfo> FindOpenReadingFramesCore(string forward, int minLength, OrfSearchOptions options, CancellationToken cancellationToken)
    {
        var code = options.GeneticCode ?? GeneticCode.Standard;
        Func<string, bool>? isStart = options.StartCodons switch
        {
            OrfStartCodons.AtgOnly => static codon => codon == StartCodon,
            OrfStartCodons.GeneticCodeStarts => code.IsStartCodon,
            OrfStartCodons.AnySenseCodon => null,
            _ => throw new ArgumentOutOfRangeException(nameof(options), options.StartCodons, "Unknown start-codon mode."),
        };
        if (!Enum.IsDefined(options.Strands))
            throw new ArgumentOutOfRangeException(nameof(options), options.Strands, "Unknown strand choice.");

        // NCBI ORFfinder -b/-e: the region is cut out first and scanned as a sequence of its own (ORFs that run off a region end are
        // open there; reverse frames are anchored at the region end), coordinates are mapped back by regionStart.
        var (regionStart, regionLength) = ResolveRegion(forward.Length, options);
        string scanned = regionStart == 0 && regionLength == forward.Length ? forward : forward.Substring(regionStart, regionLength);

        // Eager validation (genetic code with dual-coding stops) before the lazy enumeration starts.
        bool plus = options.Strands != OrfStrands.Minus;
        bool minus = options.Strands != OrfStrands.Plus;
        var forwardRegions = OrfScanner.FindRegions(scanned, code, isStart, options.NestedOrfs, options.IncludeOpenEnded,
            options.DualCodingStopsTerminate, cancellationToken);
        string? reverse = minus ? DnaSequence.GetReverseComplementString(scanned) : null;
        var reverseRegions = reverse is null
            ? null
            : OrfScanner.FindRegions(reverse, code, isStart, options.NestedOrfs, options.IncludeOpenEnded,
                options.DualCodingStopsTerminate, cancellationToken);

        int effectiveMinLength = EffectiveMinimumLength(minLength, options);
        var orfs = EnumerateOrfs(scanned, plus ? forwardRegions : null, reverse, reverseRegions, effectiveMinLength,
            options.MinLengthExcludesStop, regionStart, regionLength);
        return options.IgnoreNestedOrfs ? RemoveNestedOrfs(orfs, regionLength) : orfs;
    }

    // NCBI ORFfinder -ml: "Value less than 30 is automatically changed by 30", then the nucleotide length is used in whole codons
    // (measured on the binary, A1-18: with -ml 31 an ORF of 30 coding nt is reported, -ml 74 admits 72 nt, -ml 76 admits 75 nt).
    private const int OrfFinderMinimumLength = 30;

    private static int EffectiveMinimumLength(int minLength, OrfSearchOptions options)
    {
        if (!options.OrfFinderMinimumLength)
            return minLength;

        int clamped = Math.Max(minLength, OrfFinderMinimumLength);
        return clamped - clamped % CodonLength;
    }

    // NCBI ORFfinder -b/-e, 1-based inclusive: -b 0 is -b 1, -e 0 and -e beyond the sequence mean "to the end"; a negative start,
    // an end below the start (measured: "Interval Start value is expected to be smaller than Stop", "Requested interval start -5
    // is out of range") are errors. Returns the 0-based offset and the length of the region.
    private static (int Start, int Length) ResolveRegion(int sequenceLength, OrfSearchOptions options)
    {
        int begin = options.RegionStart;
        int end = options.RegionEnd;
        if (begin < 0)
            throw new ArgumentOutOfRangeException(nameof(options), begin, "RegionStart must not be negative (1-based; 0 is read as 1).");
        if (end < 0)
            throw new ArgumentOutOfRangeException(nameof(options), end, "RegionEnd must not be negative (1-based inclusive; 0 = to the end).");
        if (sequenceLength == 0)
            return (0, 0);

        begin = Math.Max(begin, 1);
        if (end == 0 || end > sequenceLength)
            end = sequenceLength;
        if (begin > end)
            throw new ArgumentOutOfRangeException(nameof(options), begin,
                $"RegionStart {begin} lies beyond RegionEnd {end} (sequence length {sequenceLength}).");

        return (begin - 1, end - begin + 1);
    }

    private static IEnumerable<OrfInfo> EnumerateOrfs(string forward, IEnumerable<OrfRegion>? forwardRegions,
        string? reverse, IEnumerable<OrfRegion>? reverseRegions, int minLength, bool excludeStop, int regionStart, int regionLength)
    {
        if (forwardRegions is not null)
        {
            foreach (var orf in ToOrfInfos(forward, forwardRegions, isReverseComplement: false, minLength, excludeStop, regionStart, regionLength))
            {
                yield return orf;
            }
        }

        if (reverse is null || reverseRegions is null)
        {
            yield break;
        }

        foreach (var orf in ToOrfInfos(reverse, reverseRegions, isReverseComplement: true, minLength, excludeStop, regionStart, regionLength))
        {
            yield return orf;
        }
    }

    private static IEnumerable<OrfInfo> ToOrfInfos(string strand, IEnumerable<OrfRegion> regions,
        bool isReverseComplement, int minLength, bool excludeStop, int regionStart, int regionLength)
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
                region.IsTerminated,
                regionStart,
                regionLength);
        }
    }

    // ORFfinder -n: drop every ORF whose forward span lies completely within the span of another ORF of the list (OrfScanner.FindNestedSpans);
    // the surviving ORFs keep their enumeration order.
    private static IEnumerable<OrfInfo> RemoveNestedOrfs(IEnumerable<OrfInfo> orfs, int regionLength)
    {
        var list = orfs.ToList();
        var spans = new OrfScanner.OrfSpan[list.Count];
        for (int i = 0; i < spans.Length; i++)
        {
            var orf = list[i];
            int start = orf.IsReverseComplement ? regionLength - orf.Position - orf.Length : orf.Position;
            spans[i] = new OrfScanner.OrfSpan(start, start + orf.Length - 1, orf.IsReverseComplement);
        }

        var nested = OrfScanner.FindNestedSpans(spans);
        for (int i = 0; i < list.Count; i++)
        {
            if (!nested[i])
            {
                yield return list[i];
            }
        }
    }

    #endregion
}

#region Result Types

/// <summary>
/// A known-motif hit (<see cref="GenomicAnalyzer.FindKnownMotifHits(DnaSequence, IEnumerable{string}, bool, bool, int)"/>): the 0-based start of the matching
/// window on the forward strand and whether the motif was found on the reverse strand (the window matches
/// the motif's reverse complement). <c>Mismatches</c> is the Hamming distance of the window to the motif (or to its
/// reverse complement): always 0 for an exact search, ≤ <c>maxMismatches</c> otherwise (EMBOSS fuzznuc <c>*mismatch</c> note).
/// </summary>
public readonly record struct KnownMotifHit(int Position, bool IsReverseStrand, int Mismatches = 0);

/// <summary>
/// Information about a repeated region.
/// </summary>
public readonly struct RepeatInfo
{
    /// <summary>The empty result: no repeated substring.</summary>
    public static readonly RepeatInfo None = new(string.Empty, Array.Empty<int>());

    /// <summary>Creates a repeat record.</summary>
    /// <param name="sequence">The repeated substring.</param>
    /// <param name="positions">The 0-based start positions of its occurrences.</param>
    public RepeatInfo(string sequence, IReadOnlyList<int> positions)
    {
        Sequence = sequence;
        Positions = positions;
    }

    /// <summary>The repeated substring.</summary>
    public string Sequence { get; }
    /// <summary>The 0-based start positions of its occurrences.</summary>
    public IReadOnlyList<int> Positions { get; }
    /// <summary>Length of the repeated substring.</summary>
    public int Length => Sequence.Length;
    /// <summary>Number of occurrences.</summary>
    public int Count => Positions.Count;
    /// <summary>True for <see cref="None"/> (no repeated substring).</summary>
    public bool IsEmpty => string.IsNullOrEmpty(Sequence);
}

/// <summary>
/// Information about a tandem repeat (consecutive repeating unit).
/// </summary>
public readonly struct TandemRepeat
{
    /// <summary>Creates a tandem-repeat record.</summary>
    /// <param name="unit">The repeat unit (the least rotation when canonicalised).</param>
    /// <param name="position">0-based start of the array.</param>
    /// <param name="repetitions">Number of complete copies of the unit.</param>
    /// <param name="rotation">Left rotation (in bases) from the spelled unit to <paramref name="unit"/>; 0 when not canonicalised.</param>
    public TandemRepeat(string unit, int position, int repetitions, int rotation = 0)
    {
        Unit = unit;
        Position = position;
        Repetitions = repetitions;
        Rotation = rotation;
    }

    /// <summary>The unit: the bases at <see cref="Position"/> (or, with canonical rotation, their least rotation).</summary>
    public string Unit { get; }
    /// <summary>0-based start of the array in the sequence.</summary>
    public int Position { get; }
    /// <summary>Number of complete copies of the unit.</summary>
    public int Repetitions { get; }

    /// <summary>
    /// 0 unless the unit was canonicalised: <c>Unit</c> is then the actual unit rotated left by this many bases
    /// (<c>Unit = OriginalUnit[Rotation..] + OriginalUnit[..Rotation]</c>).
    /// </summary>
    public int Rotation { get; }

    /// <summary>The unit exactly as it is spelled in the sequence starting at <see cref="Position"/>.</summary>
    public string OriginalUnit => Rotation == 0 ? Unit : string.Concat(Unit.AsSpan(Unit.Length - Rotation), Unit.AsSpan(0, Unit.Length - Rotation));

    /// <summary>Length of the array in bases (unit length times repetitions).</summary>
    public int TotalLength => Unit.Length * Repetitions;

    /// <summary>The bases of the repeat as they occur in the sequence (complete copies only).</summary>
    public string FullSequence => string.Concat(Enumerable.Repeat(OriginalUnit, Repetitions));
}

/// <summary>
/// Information about a palindromic sequence.
/// </summary>
public readonly struct PalindromeInfo
{
    /// <summary>Creates a palindrome record.</summary>
    /// <param name="sequence">The palindromic substring.</param>
    /// <param name="position">Its 0-based start position.</param>
    public PalindromeInfo(string sequence, int position)
    {
        Sequence = sequence;
        Position = position;
    }

    /// <summary>The palindromic substring.</summary>
    public string Sequence { get; }
    /// <summary>0-based start position.</summary>
    public int Position { get; }
    /// <summary>Length of the palindrome.</summary>
    public int Length => Sequence.Length;
}

/// <summary>
/// Information about a common region between two sequences.
/// </summary>
public readonly struct CommonRegion
{
    /// <summary>The empty result: no common region (positions -1).</summary>
    public static readonly CommonRegion None = new(string.Empty, -1, -1);

    /// <summary>Creates a common-region record.</summary>
    /// <param name="sequence">The common string as it reads in the first sequence.</param>
    /// <param name="positionInFirst">0-based start in the first sequence.</param>
    /// <param name="positionInSecond">0-based forward-strand start in the second sequence.</param>
    /// <param name="isReverseComplement">True when the region is shared with the reverse complement of the second sequence.</param>
    public CommonRegion(string sequence, int positionInFirst, int positionInSecond, bool isReverseComplement = false)
    {
        Sequence = sequence;
        PositionInFirst = positionInFirst;
        PositionInSecond = positionInSecond;
        IsReverseComplement = isReverseComplement;
    }

    /// <summary>The common string as it reads in the first sequence.</summary>
    public string Sequence { get; }

    /// <summary>0-based start in the first sequence.</summary>
    public int PositionInFirst { get; }

    /// <summary>
    /// 0-based start in the second sequence on the <b>forward</b> strand: for a reverse-complement region
    /// (<see cref="IsReverseComplement"/>) the window of the second sequence that spells the reverse complement of <see cref="Sequence"/>
    /// (MUMmer <c>-c</c> convention: positions of reverse matches refer to the original query strand).
    /// </summary>
    public int PositionInSecond { get; }

    /// <summary>True when the region is shared with the reverse complement of the second sequence (<c>mummer -r / -b</c>).</summary>
    public bool IsReverseComplement { get; }
    /// <summary>Length of the common region.</summary>
    public int Length => Sequence.Length;
    /// <summary>True for <see cref="None"/> (no common region).</summary>
    public bool IsEmpty => string.IsNullOrEmpty(Sequence);
}

/// <summary>
/// Information about an Open Reading Frame found by <see cref="GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int)"/>.
/// </summary>
public readonly struct OrfInfo
{
    /// <summary>Creates an ORF record.</summary>
    /// <param name="sequence">The ORF's nucleotides on the scanned strand.</param>
    /// <param name="position">0-based start in the scanned strand.</param>
    /// <param name="frame">Reading frame 1-3.</param>
    /// <param name="isReverseComplement">True when found on the reverse complement.</param>
    /// <param name="isTerminated">True when the ORF ends with a stop codon.</param>
    /// <param name="regionStart">0-based offset of the searched region in the input.</param>
    /// <param name="regionLength">Length of the searched region; -1 when unknown.</param>
    public OrfInfo(string sequence, int position, int frame, bool isReverseComplement, bool isTerminated = true,
        int regionStart = 0, int regionLength = -1)
    {
        Sequence = sequence;
        Position = position;
        Frame = frame;
        IsReverseComplement = isReverseComplement;
        IsTerminated = isTerminated;
        RegionStart = regionStart;
        RegionLength = regionLength;
    }

    /// <summary>The ORF's nucleotides, start codon through the stop codon inclusive (through the last complete codon for an open-ended ORF), 5'→3' on the scanned strand.</summary>
    public string Sequence { get; }

    /// <summary>0-based start of the ORF in the <i>scanned strand</i>: the searched region of the input sequence (the whole input by default), or its reverse complement when <see cref="IsReverseComplement"/>. See <see cref="GetForwardStrandRange()"/> for input coordinates.</summary>
    public int Position { get; }

    /// <summary>Reading frame 1–3 within the scanned strand (offset <c>Position % 3</c> + 1).</summary>
    public int Frame { get; }

    /// <summary>True when the ORF was found on the reverse complement of the input.</summary>
    public bool IsReverseComplement { get; }

    /// <summary>True when the ORF ends with a stop codon; false for an open-ended (3'-partial) ORF that reaches the strand end (only reported with <see cref="OrfSearchOptions.IncludeOpenEnded"/>).</summary>
    public bool IsTerminated { get; }

    /// <summary>
    /// 0-based offset of the searched region in the input sequence (<see cref="OrfSearchOptions.RegionStart"/> − 1; 0 for a whole-sequence
    /// search): <see cref="Position"/> is relative to it.
    /// </summary>
    public int RegionStart { get; }

    /// <summary>Length of the searched region (the scanned strand); −1 for an <see cref="OrfInfo"/> constructed without region information (then the input length passed to <see cref="GetForwardStrandRange(int)"/> is the strand length).</summary>
    public int RegionLength { get; }

    /// <summary>
    /// True when the ORF begins within the first three nucleotides of the scanned strand (the searched region), so that no complete
    /// codon precedes it and an upstream in-frame stop was not observed — the ORF may extend further upstream (5'-partial).
    /// This is NCBI ORFfinder's rule (measured on the 0.4.3 binary: <c>fuzz-from</c> / "no-left" exactly for starts at strand
    /// offset 0, 1, 2, also for an ORF that begins with ATG). The complement for the other end is <c>!<see cref="IsTerminated"/></c>
    /// (3'-partial).
    /// </summary>
    public bool IsStartPartial => Position < CodonNucleotides;

    private const int CodonNucleotides = 3;

    /// <summary>Length of the ORF in nucleotides (stop codon included).</summary>
    public int Length => Sequence.Length;
    /// <summary>Number of complete codons (stop codon included).</summary>
    public int CodonCount => Sequence.Length / 3;

    /// <summary>
    /// The ORF's location in input (forward-strand) coordinates: 0-based inclusive <c>(Start, End)</c> with
    /// <c>Start ≤ End</c>, the stop codon included. For a reverse-complement ORF the span is the one whose reverse
    /// complement is <see cref="Sequence"/>: <c>Start = RegionStart + L − Position − Length</c>, <c>End = RegionStart + L − 1 − Position</c>
    /// with L the length of the scanned strand (the input length, or the region length when a region was searched).
    /// </summary>
    /// <param name="sequenceLength">Length of the searched input sequence.</param>
    public (int Start, int End) GetForwardStrandRange(int sequenceLength)
    {
        int strandLength = RegionLength >= 0 ? RegionLength : sequenceLength;
        ArgumentOutOfRangeException.ThrowIfLessThan(sequenceLength, RegionStart + Math.Max(strandLength, Position + Length));
        return IsReverseComplement
            ? (RegionStart + strandLength - Position - Length, RegionStart + strandLength - 1 - Position)
            : (RegionStart + Position, RegionStart + Position + Length - 1);
    }

    /// <summary>
    /// <see cref="GetForwardStrandRange(int)"/> for an ORF reported by <c>FindOpenReadingFrames</c>, which knows the region it was found in
    /// (no input length needed).
    /// </summary>
    /// <exception cref="InvalidOperationException">The <see cref="OrfInfo"/> was constructed without <see cref="RegionLength"/>.</exception>
    public (int Start, int End) GetForwardStrandRange()
    {
        if (RegionLength < 0)
            throw new InvalidOperationException("This OrfInfo carries no region length; use GetForwardStrandRange(int sequenceLength).");

        return GetForwardStrandRange(RegionStart + RegionLength);
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

    /// <summary>Scan the reverse complement as well (default true). A compatible view of <see cref="Strands"/>: <c>true</c> = <see cref="OrfStrands.Both"/>, <c>false</c> = <see cref="OrfStrands.Plus"/>; reading it gives <c>Strands == Both</c>.</summary>
    public bool SearchBothStrands
    {
        get => Strands == OrfStrands.Both;
        init => Strands = value ? OrfStrands.Both : OrfStrands.Plus;
    }

    /// <summary>Which strands are scanned (NCBI ORFfinder <c>-strand both|plus|minus</c>; default <see cref="OrfStrands.Both"/>).</summary>
    public OrfStrands Strands { get; init; } = OrfStrands.Both;

    /// <summary>
    /// NCBI ORFfinder <c>-n true</c> ("Ignore nested ORFs (completely placed within another)"): drop every ORF whose forward-strand span
    /// lies completely inside the span of another reported ORF of <i>any</i> reading frame on <i>either</i> strand (see
    /// <see cref="OrfScanner.FindNestedSpans"/>; a shared start or end counts, identical spans: the plus-strand ORF is dropped). Not to be
    /// confused with <see cref="NestedOrfs"/>, which is about in-frame starts sharing one stop. The reference ORFs are the ORFs that
    /// pass <c>minLength</c> and <see cref="Strands"/>. Default false.
    /// </summary>
    public bool IgnoreNestedOrfs { get; init; }

    /// <summary>
    /// NCBI ORFfinder <c>-b</c>: 1-based first nucleotide of the searched region (default 1; 0 is read as 1). The region is cut out and
    /// scanned as a sequence of its own: ORFs that run off its ends are open there (<see cref="OrfInfo.IsStartPartial"/> /
    /// <see cref="OrfInfo.IsTerminated"/>), reverse-strand frames are anchored at the region end. <see cref="OrfInfo.Position"/> is relative to
    /// the region; <see cref="OrfInfo.GetForwardStrandRange()"/> gives input coordinates.
    /// </summary>
    public int RegionStart { get; init; } = 1;

    /// <summary>NCBI ORFfinder <c>-e</c>: 1-based last nucleotide (inclusive) of the searched region; 0 (default) and values beyond the sequence mean "to the end of the sequence".</summary>
    public int RegionEnd { get; init; }

    /// <summary>
    /// NCBI ORFfinder <c>-ml</c> rule for the <c>minLength</c> argument: values below 30 are changed to 30 (USAGE: "Value less than 30 is
    /// automatically changed by 30"), and the length is used in whole codons — rounded down to a multiple of three (measured:
    /// <c>-ml 31</c> reports ORFs of 30 coding nucleotides, <c>-ml 74</c> 72, <c>-ml 76</c> 75). Combine with
    /// <see cref="MinLengthExcludesStop"/> for the complete ORFfinder convention. Default false.
    /// </summary>
    public bool OrfFinderMinimumLength { get; init; }

    /// <summary>
    /// true: for the genetic codes whose stop codons also code for an amino acid (NCBI tables 27, 28, 31) the codons the table lists as
    /// stops (<see cref="GeneticCode.StopCodons"/>) end ORFs instead of the code being rejected — what NCBI ORFfinder 0.4.3 does for
    /// <c>-g 27/28/31</c> (measured). Default false (reject).
    /// </summary>
    public bool DualCodingStopsTerminate { get; init; }

    /// <summary>
    /// The NCBI ORFfinder 0.4.3 convention, measured against the binary (A1-18, 1750 runs): the start codons of the genetic code
    /// (<c>-s 1</c>; use <c>StartCodons = AtgOnly</c> for <c>-s 0</c>, <c>AnySenseCodon</c> for <c>-s 2</c>), the most upstream start per
    /// stop, 3'-partial ORFs, <c>minLength</c> without the stop codon with the <c>-ml</c> floor of 30 / whole codons, and the dual-coding
    /// tables 27/28/31. Pass ORFfinder's <c>-ml</c> (default 75) as <c>minLength</c>; add <see cref="GeneticCode"/> (<c>-g</c>),
    /// <see cref="Strands"/>, <see cref="IgnoreNestedOrfs"/> (<c>-n</c>), <see cref="RegionStart"/>/<see cref="RegionEnd"/> (<c>-b/-e</c>) with
    /// <c>with { ... }</c>. <c>-c</c> (circular) has no counterpart: the binary ignores it.
    /// </summary>
    public static OrfSearchOptions OrfFinderCompatible { get; } = new()
    {
        StartCodons = OrfStartCodons.GeneticCodeStarts,
        NestedOrfs = false,
        IncludeOpenEnded = true,
        MinLengthExcludesStop = true,
        OrfFinderMinimumLength = true,
        DualCodingStopsTerminate = true,
    };
}

/// <summary>Strands scanned by <see cref="GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int, OrfSearchOptions)"/> (NCBI ORFfinder <c>-strand</c>).</summary>
public enum OrfStrands
{
    /// <summary>Forward strand and its reverse complement (default).</summary>
    Both,

    /// <summary>Forward (plus) strand only.</summary>
    Plus,

    /// <summary>Reverse complement (minus) strand only.</summary>
    Minus,
}

#endregion
