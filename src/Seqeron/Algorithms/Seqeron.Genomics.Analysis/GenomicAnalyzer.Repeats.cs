using Seqeron.Genomics.Alignment;
using System.Text;

namespace Seqeron.Genomics.Analysis;

// B09 A2-33 / A2-35: repeat search on both strands, repeat pairs (MUMmer repeat-match) and the string overloads of the
// suffix-tree repeat finders. Every method here reuses the one repeat engine (EnumerateRepeatsCore: suffix-tree
// traversal for repeats / Gusfield maximal repeats, in GenomicAnalyzer.cs) or the canonical RepeatFinder pair engines.
public static partial class GenomicAnalyzer
{
    #region Text preparation (non-matching symbols, reverse strand)

    /// <summary>The nucleotide alphabet of the string overloads: only A, C, G and T can match (case-insensitive).</summary>
    private const string NucleotideAlphabet = "ACGT";

    // Separator characters are the code points above 'Z'; each maximal run of non-matching symbols gets its own.
    private const int FirstSeparator = 0x100;
    private const int LastSeparator = 0xD7FF;

    public static bool[] BuildLetterMask(string alphabet)
    {
        ArgumentNullException.ThrowIfNull(alphabet);
        var letters = new bool[128];
        foreach (char c in alphabet)
        {
            char u = char.ToUpperInvariant(c);
            if (u is < 'A' or > 'Z')
                throw new ArgumentException($"The alphabet may contain only the letters A-Z; found '{c}'.", nameof(alphabet));
            letters[u] = true;
        }

        return letters;
    }

    private static char NextSeparator(ref int next)
    {
        if (next > LastSeparator)
            throw new ArgumentException(
                $"The sequence has more than {LastSeparator - FirstSeparator + 1} separate runs of non-matching symbols; split it at them and search the pieces.");
        return (char)next++;
    }

    // Upper-cases the input and replaces every maximal run of symbols outside the alphabet by a run of one separator
    // character of the same length, unique to that run: positions are preserved, a separator never equals a letter or
    // another run's separator, so no repeat can contain it and none can be extended over it.
    private static string MaskText(string sequence, bool[] letters, ref int nextSeparator)
    {
        var sb = new StringBuilder(sequence.Length);
        char separator = '\0';
        bool inRun = false;
        foreach (char raw in sequence)
        {
            char c = char.ToUpperInvariant(raw);
            if (c < 128 && letters[c])
            {
                sb.Append(c);
                inRun = false;
            }
            else
            {
                if (!inRun)
                {
                    separator = NextSeparator(ref nextSeparator);
                    inRun = true;
                }

                sb.Append(separator);
            }
        }

        return sb.ToString();
    }

    public static global::SuffixTree.SuffixTree BuildMaskedTree(string sequence, bool[] letters)
    {
        int next = FirstSeparator;
        return global::SuffixTree.SuffixTree.Build(MaskText(sequence, letters, ref next));
    }

    // Reverse complement of the A/C/G/T symbols; every other symbol becomes N (it never matches anyway).
    private static string ReverseComplementAcgt(string upper)
    {
        var chars = new char[upper.Length];
        for (int i = 0; i < upper.Length; i++)
        {
            char c = upper[upper.Length - 1 - i];
            chars[i] = c is 'A' or 'C' or 'G' or 'T' ? SequenceExtensions.GetComplementBase(c) : 'N';
        }

        return new string(chars);
    }

    #endregion

    #region String overloads of the repeat finders (A2-35)

    /// <summary>
    /// <see cref="FindLongestRepeat(DnaSequence)"/> on a nucleotide <b>string</b> (or any letter alphabet): a genome with N runs,
    /// IUPAC codes or lower case, which <see cref="DnaSequence"/> (A/C/G/T only) cannot hold.
    /// </summary>
    /// <remarks>
    /// Contract (the B04 / B05 ambiguity rule, MUMmer <c>-n</c>): the input is case-insensitive and only the symbols of
    /// <paramref name="alphabet"/> (default A, C, G, T) can repeat; every other symbol (N, another IUPAC code, U, a gap, ...)
    /// matches nothing and ends a repeat, so <c>NNNN…</c> is never a repeat and an ordinary repeat next to an N run is reported
    /// exactly as without it. Positions are 0-based in the input string. For an A/C/G/T string the result equals the
    /// <see cref="DnaSequence"/> overload. A protein is searched with <c>alphabet: "ACDEFGHIKLMNPQRSTVWY"</c>. An empty string has no repeat.
    /// </remarks>
    /// <param name="sequence">Sequence to search. Must not be <c>null</c>.</param>
    /// <param name="alphabet">The matching symbols, letters A-Z (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="alphabet"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentException">The alphabet holds a non-letter, or the sequence has more than 55 000 separate runs of non-matching symbols.</exception>
    public static RepeatInfo FindLongestRepeat(string sequence, string alphabet = NucleotideAlphabet)
    {
        var all = FindAllLongestRepeats(sequence, alphabet);
        return all.Count == 0 ? RepeatInfo.None : all[0];
    }

    /// <summary>
    /// <see cref="FindAllLongestRepeats(DnaSequence)"/> on a string; the matching rule is that of
    /// <see cref="FindLongestRepeat(string, string)"/>.
    /// </summary>
    /// <param name="sequence">Sequence to search. Must not be <c>null</c>.</param>
    /// <param name="alphabet">The matching symbols, letters A-Z (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="alphabet"/> is <c>null</c>.</exception>
    public static IReadOnlyList<RepeatInfo> FindAllLongestRepeats(string sequence, string alphabet = NucleotideAlphabet)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var letters = BuildLetterMask(alphabet);
        if (sequence.Length == 0)
            return Array.Empty<RepeatInfo>();
        return LongestRepeatsCore(BuildMaskedTree(sequence, letters), default);
    }

    /// <summary>
    /// <see cref="FindRepeats(DnaSequence, int)"/> on a string; the matching rule is that of
    /// <see cref="FindLongestRepeat(string, string)"/> (non-matching symbols end every repeat).
    /// </summary>
    /// <param name="sequence">Sequence to search. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IEnumerable<RepeatInfo> FindRepeats(string sequence, int minLength)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: false, NucleotideAlphabet, default);

    /// <summary><see cref="FindRepeats(string, int)"/> that observes <paramref name="cancellationToken"/>.</summary>
    /// <param name="sequence">Sequence to search. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    public static IEnumerable<RepeatInfo> FindRepeats(string sequence, int minLength, CancellationToken cancellationToken)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: false, NucleotideAlphabet, cancellationToken);

    /// <summary><see cref="FindRepeats(string, int)"/> over an explicit matching alphabet (e.g. the 20 amino acids for a protein).</summary>
    /// <param name="sequence">Sequence to search. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="alphabet">The matching symbols, letters A-Z (case-insensitive).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    public static IEnumerable<RepeatInfo> FindRepeats(string sequence, int minLength, string alphabet, CancellationToken cancellationToken = default)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: false, alphabet, cancellationToken);

    /// <summary>
    /// <see cref="FindMaximalRepeats(DnaSequence, int)"/> on a string. A non-matching symbol is its own left and right context, so an
    /// occurrence next to an N is not extendable (<c>ACGTNxxx</c> / <c>ACGTNyyy</c> share the maximal repeat <c>ACGT</c>).
    /// </summary>
    /// <param name="sequence">Sequence to search. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IEnumerable<RepeatInfo> FindMaximalRepeats(string sequence, int minLength)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: true, NucleotideAlphabet, default);

    /// <summary><see cref="FindMaximalRepeats(string, int)"/> that observes <paramref name="cancellationToken"/>.</summary>
    /// <param name="sequence">Sequence to search. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    public static IEnumerable<RepeatInfo> FindMaximalRepeats(string sequence, int minLength, CancellationToken cancellationToken)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: true, NucleotideAlphabet, cancellationToken);

    /// <summary><see cref="FindMaximalRepeats(string, int)"/> over an explicit matching alphabet.</summary>
    /// <param name="sequence">Sequence to search. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="alphabet">The matching symbols, letters A-Z (case-insensitive).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    public static IEnumerable<RepeatInfo> FindMaximalRepeats(string sequence, int minLength, string alphabet, CancellationToken cancellationToken = default)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: true, alphabet, cancellationToken);

    /// <summary><see cref="FindRepeats(string, int)"/> with explicit <see cref="RepeatOutputLimits"/> (R14-4b; see <see cref="FindRepeats(DnaSequence, int, RepeatOutputLimits, CancellationToken)"/>).</summary>
    /// <param name="sequence">Nucleotide string. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="limits">The limit of the output.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="limits"/> is <c>null</c>.</exception>
    /// <exception cref="RepeatOutputLimitException">(on enumeration) the output exceeds <paramref name="limits"/>.</exception>
    public static IEnumerable<RepeatInfo> FindRepeats(string sequence, int minLength, RepeatOutputLimits limits, CancellationToken cancellationToken = default)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: false, NucleotideAlphabet, cancellationToken, limits ?? throw new ArgumentNullException(nameof(limits)));

    /// <summary><see cref="FindMaximalRepeats(string, int)"/> with explicit <see cref="RepeatOutputLimits"/> (R14-4b).</summary>
    /// <param name="sequence">Nucleotide string. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="limits">The limit of the output.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="limits"/> is <c>null</c>.</exception>
    /// <exception cref="RepeatOutputLimitException">(on enumeration) the output exceeds <paramref name="limits"/>.</exception>
    public static IEnumerable<RepeatInfo> FindMaximalRepeats(string sequence, int minLength, RepeatOutputLimits limits, CancellationToken cancellationToken = default)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: true, NucleotideAlphabet, cancellationToken, limits ?? throw new ArgumentNullException(nameof(limits)));

    /// <summary><see cref="FindRepeats(string, int, string, CancellationToken)"/> with explicit <see cref="RepeatOutputLimits"/> (R14-4b).</summary>
    /// <param name="sequence">Text. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="alphabet">The matching alphabet (letters A-Z).</param>
    /// <param name="limits">The limit of the output.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="limits"/> is <c>null</c>.</exception>
    /// <exception cref="RepeatOutputLimitException">(on enumeration) the output exceeds <paramref name="limits"/>.</exception>
    public static IEnumerable<RepeatInfo> FindRepeats(string sequence, int minLength, string alphabet, RepeatOutputLimits limits, CancellationToken cancellationToken = default)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: false, alphabet, cancellationToken, limits ?? throw new ArgumentNullException(nameof(limits)));

    /// <summary><see cref="FindMaximalRepeats(string, int, string, CancellationToken)"/> with explicit <see cref="RepeatOutputLimits"/> (R14-4b).</summary>
    /// <param name="sequence">Text. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="alphabet">The matching alphabet (letters A-Z).</param>
    /// <param name="limits">The limit of the output.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="limits"/> is <c>null</c>.</exception>
    /// <exception cref="RepeatOutputLimitException">(on enumeration) the output exceeds <paramref name="limits"/>.</exception>
    public static IEnumerable<RepeatInfo> FindMaximalRepeats(string sequence, int minLength, string alphabet, RepeatOutputLimits limits, CancellationToken cancellationToken = default)
        => EnumerateTextRepeats(sequence, minLength, maximalOnly: true, alphabet, cancellationToken, limits ?? throw new ArgumentNullException(nameof(limits)));

    private static IEnumerable<RepeatInfo> EnumerateTextRepeats(string sequence, int minLength, bool maximalOnly, string alphabet, CancellationToken cancellationToken, RepeatOutputLimits? limits = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var letters = BuildLetterMask(alphabet);
        return Iterate();

        IEnumerable<RepeatInfo> Iterate()
        {
            AlignmentCancellation.Poll(cancellationToken);
            if (sequence.Length == 0)
                yield break;
            foreach (var repeat in EnumerateRepeatsCore(BuildMaskedTree(sequence, letters), minLength, maximalOnly, cancellationToken, limits))
                yield return repeat;
        }
    }

    // Every longest repeat of a (possibly separator-masked) tree: the deepest repeat nodes whose label is separator-free,
    // ordered by first occurrence (ascending) like SuffixTree.FindAllLongestRepeatedSubstrings; positions ascending.
    private static List<RepeatInfo> LongestRepeatsCore(global::SuffixTree.SuffixTree tree, CancellationToken cancellationToken)
    {
        var loci = new RepeatLocusCollector(1, maximalOnly: false, tree.Text, cancellationToken);
        tree.Traverse(loci);

        int maxDepth = 0;
        foreach (var (_, _, depth, _) in loci.Nodes)
            maxDepth = Math.Max(maxDepth, depth);

        var result = new List<RepeatInfo>();
        if (maxDepth == 0)
            return result;
        foreach (var (edgeEnd, _, depth, _) in loci.Nodes)
        {
            if (depth != maxDepth)
                continue;
            string label = tree.Text.Substring(edgeEnd - depth, depth);
            var positions = new List<int>(tree.FindAllOccurrences(label));
            positions.Sort();
            result.Add(new RepeatInfo(label, positions));
        }

        result.Sort((a, b) => a.Positions[0].CompareTo(b.Positions[0]));
        return result;
    }

    #endregion

    #region Both-strand repeats (A2-33)

    /// <summary>
    /// Repeated substrings on <b>both strands</b> (MUMmer <c>repeat-match</c> default, Vmatch <c>-p</c>, REPuter forward + palindromic
    /// repeats): the classes of substrings that occur at least twice in <c>S # revcomp(S)</c>, i.e. at least twice counting
    /// occurrences of the string itself on the forward strand and of its reverse complement on the other.
    /// </summary>
    /// <remarks>
    /// <para>
    /// A class is a string α and its reverse complement α′, reported once as the lexicographically (ordinal) smaller of the two;
    /// <see cref="StrandedRepeatInfo.Positions"/> are the 0-based starts of α in the sequence and
    /// <see cref="StrandedRepeatInfo.ReversePositions"/> the 0-based forward-strand starts of the windows that spell α′ (the occurrences of α on the
    /// reverse strand — fuzznuc / GFF coordinates). A palindrome (α = α′, GAATTC) occurs on both strands at each of its sites, so one site is
    /// already a repeat (<c>repeat-match</c> prints it as the pair <c>i = k</c>); both position lists are then equal.
    /// With <paramref name="maximalOnly"/> only maximal repeats are kept: α is left- and right-diverse over all its occurrences in
    /// <c>S # revcomp(S)</c>, equivalently α is the string of at least one maximal repeated pair of
    /// <see cref="FindRepeatPairs(DnaSequence, int, StrandSelection)"/> (Gusfield 1997 §7.12 on the doubled text).
    /// </para>
    /// <para>
    /// Only A/C/G/T match (case-insensitive); N, other IUPAC codes, U and gaps never match and end a repeat (MUMmer <c>-n</c>, B04 convention).
    /// The classes are ordered by their string (ordinal). One suffix tree of the doubled text (O(n) construction, Ukkonen) is built per call
    /// and traversed once; the result is materialised because the order is by string.
    /// Verified against <c>repeat-match</c> 3.23 (the class set of its output equals the maximal classes on every ACGT case) and a brute-force
    /// class enumeration, see docs/Evidence/GENOMIC-REPEAT-001-Evidence.md.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="maximalOnly">Report only maximal repeats (default: every repeated substring).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="limits">The limit of the output (R14-4b; <c>null</c> = <see cref="RepeatOutputLimits.Default"/>): a larger one throws <see cref="RepeatOutputLimitException"/> before anything is materialised.</param>
    /// <exception cref="RepeatOutputLimitException">The output exceeds <paramref name="limits"/>.</exception>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IReadOnlyList<StrandedRepeatInfo> FindRepeatsBothStrands(
        DnaSequence sequence, int minLength, bool maximalOnly = false, CancellationToken cancellationToken = default, RepeatOutputLimits? limits = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return StrandedRepeatsCore(sequence.Sequence, minLength, maximalOnly, longestOnly: false, cancellationToken, limits);
    }

    /// <summary>
    /// <see cref="FindRepeatsBothStrands(DnaSequence, int, bool, CancellationToken)"/> on a nucleotide <b>string</b> (N runs, IUPAC codes, lower
    /// case): non-A/C/G/T symbols never match and end a repeat; positions are 0-based in the input string.
    /// </summary>
    /// <param name="sequence">Nucleotide string. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length; values ≤ 0 are treated as 1.</param>
    /// <param name="maximalOnly">Report only maximal repeats.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <param name="limits">The limit of the output (R14-4b; <c>null</c> = <see cref="RepeatOutputLimits.Default"/>): a larger one throws <see cref="RepeatOutputLimitException"/> before anything is materialised.</param>
    /// <exception cref="RepeatOutputLimitException">The output exceeds <paramref name="limits"/>.</exception>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IReadOnlyList<StrandedRepeatInfo> FindRepeatsBothStrands(
        string sequence, int minLength, bool maximalOnly = false, CancellationToken cancellationToken = default, RepeatOutputLimits? limits = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return StrandedRepeatsCore(sequence, minLength, maximalOnly, longestOnly: false, cancellationToken, limits);
    }

    /// <summary>
    /// Every longest repeat class over both strands (the class of the longest substring that occurs at least twice in
    /// <c>S # revcomp(S)</c>), ordered by string; empty when nothing repeats. Longest repeats are maximal.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IReadOnlyList<StrandedRepeatInfo> FindAllLongestRepeatsBothStrands(DnaSequence sequence, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return StrandedRepeatsCore(sequence.Sequence, 1, maximalOnly: false, longestOnly: true, cancellationToken);
    }

    /// <summary><see cref="FindAllLongestRepeatsBothStrands(DnaSequence, CancellationToken)"/> on a nucleotide string.</summary>
    /// <param name="sequence">Nucleotide string. Must not be <c>null</c>.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static IReadOnlyList<StrandedRepeatInfo> FindAllLongestRepeatsBothStrands(string sequence, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return StrandedRepeatsCore(sequence, 1, maximalOnly: false, longestOnly: true, cancellationToken);
    }

    /// <summary>
    /// The first of <see cref="FindAllLongestRepeatsBothStrands(DnaSequence, CancellationToken)"/> (smallest string), or
    /// <see cref="StrandedRepeatInfo.None"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    public static StrandedRepeatInfo FindLongestRepeatBothStrands(DnaSequence sequence)
    {
        var all = FindAllLongestRepeatsBothStrands(sequence);
        return all.Count == 0 ? StrandedRepeatInfo.None : all[0];
    }

    private static List<StrandedRepeatInfo> StrandedRepeatsCore(string sequence, int minLength, bool maximalOnly, bool longestOnly, CancellationToken cancellationToken, RepeatOutputLimits? limits = null)
    {
        AlignmentCancellation.Poll(cancellationToken);
        var result = new List<StrandedRepeatInfo>();
        int n = sequence.Length;
        if (n == 0)
            return result;

        // T = S # revcomp(S): every symbol outside A/C/G/T is a run separator, '#' is one more; the two halves use disjoint
        // separators (an N before the same site read from the other strand is a different context).
        var letters = BuildLetterMask(NucleotideAlphabet);
        int next = FirstSeparator;
        string forward = MaskText(sequence, letters, ref next);
        char middle = NextSeparator(ref next);
        string reverse = MaskText(ReverseComplementAcgt(sequence.ToUpperInvariant()), letters, ref next);
        var tree = global::SuffixTree.SuffixTree.Build(forward + middle + reverse);

        IEnumerable<RepeatInfo> repeats = longestOnly
            ? LongestRepeatsCore(tree, cancellationToken)
            : EnumerateRepeatsCore(tree, minLength, maximalOnly, cancellationToken, limits);
        foreach (var repeat in repeats)
        {
            AlignmentCancellation.Poll(cancellationToken);
            string alpha = repeat.Sequence;
            int length = alpha.Length;
            if (string.CompareOrdinal(alpha, DnaSequence.GetReverseComplementString(alpha)) > 0)
                continue; // the class is reported through its smaller member

            var forwardStarts = new List<int>();
            var reverseStarts = new List<int>();
            foreach (int p in repeat.Positions)
            {
                if (p < n)
                    forwardStarts.Add(p);
                else
                    reverseStarts.Add(n - (p - n - 1) - length); // an occurrence of alpha in revcomp(S) = a window of S spelling revcomp(alpha)
            }

            forwardStarts.Sort();
            reverseStarts.Sort();
            result.Add(new StrandedRepeatInfo(alpha, forwardStarts, reverseStarts));
        }

        result.Sort((a, b) => string.CompareOrdinal(a.Sequence, b.Sequence));
        return result;
    }

    #endregion

    #region Repeat pairs (MUMmer repeat-match)

    /// <summary>
    /// All maximal exact repeated <b>pairs</b> of a sequence with itself, forward and/or reverse-complement — the output of MUMmer
    /// <c>repeat-match</c> (default both strands; <c>-f</c> forward only): <c>S[i..i+L) = S[j..j+L)</c> (forward) or
    /// <c>S[i..i+L) = revcomp(S[k..k+L))</c> (reverse complement, <c>i ≤ k</c>, <c>i = k</c> is a palindromic site), maximal in both directions.
    /// </summary>
    /// <remarks>
    /// This is the union of the canonical <see cref="RepeatFinder.FindDirectRepeats(DnaSequence,int,int,int)"/> and
    /// <see cref="RepeatFinder.FindReverseComplementRepeats(DnaSequence,int,int,int)"/> with every length / spacing filter lifted
    /// (<c>maxLength = int.MaxValue</c>, <c>minSpacing = int.MinValue</c>), so it is exactly the unfiltered <c>repeat-match</c> list
    /// (0-based forward-strand starts; see those methods for the coordinate conventions and the MUMmer <c>Start2</c> conversion).
    /// Only A/C/G/T match. Pairs are ordered by (first start, second start, length), a forward pair before a reverse pair on a tie.
    /// </remarks>
    /// <param name="sequence">DNA sequence to analyze. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length, ≥ 2 (<c>repeat-match -n</c>).</param>
    /// <param name="strands">Which strands to report (default both).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> &lt; 2.</exception>
    public static IReadOnlyList<RepeatPair> FindRepeatPairs(DnaSequence sequence, int minLength = 5, StrandSelection strands = StrandSelection.Both)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return MergePairs(
            strands.HasFlag(StrandSelection.Forward) ? RepeatFinder.FindDirectRepeats(sequence, minLength, int.MaxValue, int.MinValue) : [],
            strands.HasFlag(StrandSelection.ReverseComplement) ? RepeatFinder.FindReverseComplementRepeats(sequence, minLength, int.MaxValue, int.MinValue) : []);
    }

    /// <summary>
    /// <see cref="FindRepeatPairs(DnaSequence, int, StrandSelection)"/> on a nucleotide <b>string</b> (case-insensitive; N, IUPAC codes, U and gaps never
    /// match), the contract of the canonical string overloads of <see cref="RepeatFinder"/>. An empty string has no pairs.
    /// </summary>
    /// <param name="sequence">Nucleotide string. Must not be <c>null</c>.</param>
    /// <param name="minLength">Minimum repeat length, ≥ 2.</param>
    /// <param name="strands">Which strands to report (default both).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minLength"/> &lt; 2.</exception>
    public static IReadOnlyList<RepeatPair> FindRepeatPairs(string sequence, int minLength = 5, StrandSelection strands = StrandSelection.Both)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return MergePairs(
            strands.HasFlag(StrandSelection.Forward) ? RepeatFinder.FindDirectRepeats(sequence, minLength, int.MaxValue, int.MinValue) : [],
            strands.HasFlag(StrandSelection.ReverseComplement) ? RepeatFinder.FindReverseComplementRepeats(sequence, minLength, int.MaxValue, int.MinValue) : []);
    }

    private static List<RepeatPair> MergePairs(IEnumerable<DirectRepeatResult> direct, IEnumerable<ReverseComplementRepeatResult> reverseComplement)
    {
        var pairs = new List<RepeatPair>();
        foreach (var d in direct)
            pairs.Add(new RepeatPair(d.FirstPosition, d.SecondPosition, d.Length, IsReverseComplement: false));
        foreach (var r in reverseComplement)
            pairs.Add(new RepeatPair(r.FirstPosition, r.SecondPosition, r.Length, IsReverseComplement: true));
        pairs.Sort((a, b) =>
        {
            int c = a.FirstPosition.CompareTo(b.FirstPosition);
            if (c != 0) return c;
            c = a.SecondPosition.CompareTo(b.SecondPosition);
            if (c != 0) return c;
            c = a.Length.CompareTo(b.Length);
            return c != 0 ? c : a.IsReverseComplement.CompareTo(b.IsReverseComplement);
        });
        return pairs;
    }

    #endregion
}

/// <summary>
/// A repeat class over both strands (<see cref="GenomicAnalyzer.FindRepeatsBothStrands(DnaSequence, int, bool, CancellationToken)"/>): the string
/// <see cref="Sequence"/> (the smaller of a substring and its reverse complement), its forward-strand starts and the forward-strand starts of
/// the windows that spell its reverse complement (its occurrences on the reverse strand). For a palindrome both lists are equal.
/// </summary>
public readonly struct StrandedRepeatInfo
{
    /// <summary>The empty result.</summary>
    public static readonly StrandedRepeatInfo None = new(string.Empty, Array.Empty<int>(), Array.Empty<int>());

    /// <summary>Creates a stranded repeat.</summary>
    public StrandedRepeatInfo(string sequence, IReadOnlyList<int> positions, IReadOnlyList<int> reversePositions)
    {
        Sequence = sequence;
        Positions = positions;
        ReversePositions = reversePositions;
    }

    /// <summary>The class representative: the ordinal-smaller of the substring and its reverse complement.</summary>
    public string Sequence { get; }

    /// <summary>0-based starts of <see cref="Sequence"/> on the forward strand, ascending.</summary>
    public IReadOnlyList<int> Positions { get; }

    /// <summary>0-based forward-strand starts of the windows equal to the reverse complement of <see cref="Sequence"/>, ascending.</summary>
    public IReadOnlyList<int> ReversePositions { get; }

    /// <summary>Length of the repeat.</summary>
    public int Length => Sequence.Length;

    /// <summary>Number of occurrences over both strands (a palindromic site counts once per strand).</summary>
    public int Count => Positions.Count + ReversePositions.Count;

    /// <summary>True for the empty result.</summary>
    public bool IsEmpty => string.IsNullOrEmpty(Sequence);
}

/// <summary>Which strands a search reports (repeat pairs, common regions, pattern hits).</summary>
[Flags]
public enum StrandSelection
{
    /// <summary>The forward strand only (direct repeats, <c>repeat-match -f</c>, <c>mummer -n</c> default).</summary>
    Forward = 1,

    /// <summary>The reverse-complement strand only (<c>mummer -r</c>).</summary>
    ReverseComplement = 2,

    /// <summary>Both strands (the <c>repeat-match</c> default, <c>mummer -b</c>).</summary>
    Both = Forward | ReverseComplement,
}

/// <summary>
/// A maximal exact repeated pair (<see cref="GenomicAnalyzer.FindRepeatPairs(DnaSequence, int, StrandSelection)"/>): 0-based forward-strand starts of
/// both copies; for a reverse-complement pair the window at <see cref="SecondPosition"/> spells the reverse complement of the window at
/// <see cref="FirstPosition"/> (<c>FirstPosition ≤ SecondPosition</c>; equal for a palindromic site).
/// </summary>
public readonly record struct RepeatPair(int FirstPosition, int SecondPosition, int Length, bool IsReverseComplement)
{
    /// <summary>Bases between the copies (negative when they overlap): <c>SecondPosition − FirstPosition − Length</c>.</summary>
    public int Spacing => SecondPosition - FirstPosition - Length;
}
