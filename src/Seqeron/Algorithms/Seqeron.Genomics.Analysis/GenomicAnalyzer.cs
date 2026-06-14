using System;
using System.Collections.Generic;
using System.Linq;

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
    /// Each such repeated substring is the longest common prefix of two adjacent suffixes in
    /// sorted order; only those occurring ≥ 2 times are returned (a repeated substring maps to an
    /// internal suffix-tree node, never a leaf — CMU 15-451 Lecture #10 §2.1; GeeksforGeeks,
    /// Suffix Tree Application 3). Worst-case time O(n²) because adjacent-suffix prefix comparison
    /// and occurrence enumeration are over O(n) suffixes of up to O(n) length.
    /// See docs/algorithms/Repeat_Analysis/Repeat_Detection.md.
    /// </summary>
    public static IEnumerable<RepeatInfo> FindRepeats(DnaSequence sequence, int minLength)
    {
        var tree = sequence.SuffixTree;
        var found = new HashSet<string>();

        // A repeat is a NON-EMPTY substring occurring >=2 times (the path label of an internal
        // node, never the root). The empty string is not a repeat, so the effective minimum is
        // at least 1 even when the caller passes minLength <= 0. (CMU 15-451 Lecture #10 §2.1.)
        int effectiveMinLength = Math.Max(1, minLength);

        // Get all suffixes and find common prefixes
        var suffixes = tree.GetAllSuffixes().OrderBy(s => s).ToList();

        for (int i = 0; i < suffixes.Count - 1; i++)
        {
            string s1 = suffixes[i];
            string s2 = suffixes[i + 1];

            // Find longest common prefix
            int lcpLen = 0;
            int maxLen = Math.Min(s1.Length, s2.Length);
            while (lcpLen < maxLen && s1[lcpLen] == s2[lcpLen])
            {
                lcpLen++;
            }

            if (lcpLen >= effectiveMinLength)
            {
                string repeat = s1.Substring(0, lcpLen);
                if (!found.Contains(repeat))
                {
                    found.Add(repeat);
                    var positions = tree.FindAllOccurrences(repeat).OrderBy(p => p).ToList();
                    if (positions.Count >= 2)
                    {
                        yield return new RepeatInfo(repeat, positions);
                    }
                }
            }
        }
    }

    /// <summary>
    /// Finds tandem repeats (consecutive repeating units like "ATGATGATG").
    /// </summary>
    public static IEnumerable<TandemRepeat> FindTandemRepeats(DnaSequence sequence, int minUnitLength = 2, int minRepetitions = 2)
    {
        string seq = sequence.Sequence;

        for (int unitLen = minUnitLength; unitLen <= seq.Length / minRepetitions; unitLen++)
        {
            for (int start = 0; start <= seq.Length - unitLen * minRepetitions; start++)
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
            }
        }
    }

    #endregion

    #region Motif Finding

    /// <summary>
    /// Finds all occurrences of a motif (pattern) in a sequence.
    /// Time complexity: O(m) where m is motif length.
    /// </summary>
    public static IReadOnlyList<int> FindMotif(DnaSequence sequence, string motif)
    {
        if (string.IsNullOrEmpty(motif))
            return Array.Empty<int>();

        string normalizedMotif = motif.ToUpperInvariant();
        return sequence.SuffixTree.FindAllOccurrences(normalizedMotif);
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
    /// Searches a sequence for a set of known motifs, returning, for each motif that occurs,
    /// the 0-based start positions of <b>all</b> its occurrences. This is the classical exact
    /// set-matching problem ("find all occurrences of each pattern P in text T") solved over the
    /// generalized index of one text: each query is matched against this sequence's suffix tree
    /// (Gusfield 1997, ISBN 0-521-58519-8; exact-matching definition, Tufts COMP 150GEN exact.html).
    /// <para>
    /// <b>Overlapping occurrences are all reported</b> — e.g. motif "AAA" in "AAAAA" yields
    /// {0, 1, 2}, mirroring Biopython's overlap-aware semantics (<c>Seq("AAAA").count_overlap("AA") == 3</c>).
    /// </para>
    /// Motifs are upper-cased before searching (DNA is processed upper-cased, per Biopython
    /// <c>Bio.Seq</c>); the result is keyed by the upper-cased motif. Empty or whitespace-only
    /// motifs are skipped (the empty string is not a motif; mirrors <see cref="FindMotif"/>), and a
    /// motif with no occurrence is omitted from the result. Positions for each motif are returned
    /// <b>sorted ascending</b> for a deterministic, stable contract (the suffix tree enumerates
    /// occurrences in DFS order, which is not inherently sorted). Duplicate motifs that normalize to
    /// the same upper-cased key collapse to a single entry.
    /// Time complexity: O(n) suffix-tree construction (Ukkonen) plus O(|m| + occ) per motif query
    /// and O(occ·log occ) to sort each motif's positions. See
    /// docs/algorithms/Motif_Analysis/Known_Motif_Search.md.
    /// </summary>
    public static Dictionary<string, IReadOnlyList<int>> FindKnownMotifs(
        DnaSequence sequence,
        IEnumerable<string> motifs)
    {
        if (motifs is null)
        {
            throw new ArgumentNullException(nameof(motifs));
        }

        var result = new Dictionary<string, IReadOnlyList<int>>();
        var tree = sequence.SuffixTree;

        foreach (var motif in motifs)
        {
            // The empty string is not a motif (suffix-tree FindAllOccurrences("") would return
            // every position, which is not a meaningful match). Skip empty/whitespace motifs,
            // consistent with FindMotif's empty-pattern guard.
            if (string.IsNullOrWhiteSpace(motif))
            {
                continue;
            }

            string normalized = motif.ToUpperInvariant();
            if (result.ContainsKey(normalized))
            {
                continue; // duplicate motif key — already searched
            }

            var positions = tree.FindAllOccurrences(normalized);
            if (positions.Count > 0)
            {
                // Suffix-tree occurrence enumeration is DFS-order (not sorted); sort ascending
                // so each motif's positions form the deterministic exact-matching set.
                result[normalized] = positions.OrderBy(p => p).ToList();
            }
        }

        return result;
    }

    #endregion

    #region Sequence Comparison

    /// <summary>
    /// Finds the longest common <b>substring</b> (a longest <i>contiguous</i> string that is a
    /// substring of both sequences), not a gapped subsequence — Wikipedia "Longest common substring";
    /// the contiguous distinction is explicit there. Computed via the generalized suffix tree of
    /// <paramref name="sequence1"/>: the deepest node whose subtree has leaves from both strings spells
    /// the answer (Gusfield 1997, ISBN 0-521-58519-8; GeeksforGeeks "Suffix Tree Application 5").
    /// On a length tie, the substring first found in <paramref name="sequence2"/> is returned
    /// (documented tie-break of <see cref="SuffixTree.SuffixTree.LongestCommonSubstringInfo(string)"/>).
    /// Returns <see cref="CommonRegion.None"/> when there is no shared substring (including any empty input).
    /// Time complexity: O(n + m) with the generalized suffix tree (Gusfield 1997).
    /// See docs/algorithms/Sequence_Comparison/Common_Region_Detection.md.
    /// </summary>
    public static CommonRegion FindLongestCommonRegion(DnaSequence sequence1, DnaSequence sequence2)
    {
        var tree = sequence1.SuffixTree;
        var (lcs, pos1, pos2) = tree.LongestCommonSubstringInfo(sequence2.Sequence);

        if (string.IsNullOrEmpty(lcs))
        {
            return CommonRegion.None;
        }

        return new CommonRegion(lcs, pos1, pos2);
    }

    /// <summary>
    /// Finds every distinct common <b>substring</b> (contiguous, per Wikipedia "Longest common substring")
    /// of length ≥ <paramref name="minLength"/> that occurs in both sequences. For each start position in
    /// <paramref name="sequence2"/> the longest substring also present in <paramref name="sequence1"/>
    /// (located via the generalized suffix tree) is taken; distinct substrings are reported once, with the
    /// first occurrence position in <paramref name="sequence1"/> and the start position in
    /// <paramref name="sequence2"/>. A common substring is non-empty, so values of
    /// <paramref name="minLength"/> below 1 are treated as 1 (the empty string is not a region).
    /// Time complexity: O(n + m·log m) — O(n + m) suffix-tree construction plus a per-position
    /// binary search using O(m) <see cref="SuffixTree.SuffixTree.Contains(string)"/> lookups
    /// (Gusfield 1997, ISBN 0-521-58519-8). See docs/algorithms/Sequence_Comparison/Common_Region_Detection.md.
    /// </summary>
    public static IEnumerable<CommonRegion> FindCommonRegions(
        DnaSequence sequence1,
        DnaSequence sequence2,
        int minLength)
    {
        var tree = sequence1.SuffixTree;
        string seq2 = sequence2.Sequence;
        var found = new HashSet<string>();

        // A common region is a NON-EMPTY contiguous substring; the empty string is not a region,
        // so the effective minimum length is at least 1 (Wikipedia "Longest common substring").
        int effectiveMinLength = Math.Max(1, minLength);

        // Slide through sequence2 and find matches
        for (int i = 0; i < seq2.Length - effectiveMinLength + 1; i++)
        {
            // Binary search for longest match at this position
            int lo = effectiveMinLength, hi = seq2.Length - i;
            string? bestMatch = null;

            while (lo <= hi)
            {
                int mid = (lo + hi) / 2;
                string candidate = seq2.Substring(i, mid);

                if (tree.Contains(candidate))
                {
                    bestMatch = candidate;
                    lo = mid + 1;
                }
                else
                {
                    hi = mid - 1;
                }
            }

            if (bestMatch != null && !found.Contains(bestMatch))
            {
                found.Add(bestMatch);
                var positions1 = tree.FindAllOccurrences(bestMatch);
                yield return new CommonRegion(bestMatch, positions1[0], i);
            }
        }
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
        if (sequence1 is null)
        {
            throw new ArgumentNullException(nameof(sequence1));
        }

        if (sequence2 is null)
        {
            throw new ArgumentNullException(nameof(sequence2));
        }

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

    // Standard genetic code (NCBI transl_table=1): start codon ATG, stop codons TAA/TAG/TGA.
    // NCBI Genetic Codes, https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (accessed 2026-06-14).
    private const string StartCodon = "ATG";
    private static readonly string[] StopCodons = { "TAA", "TAG", "TGA" };

    // A codon is three nucleotides.
    private const int CodonLength = 3;

    // Number of reading frames per strand (offsets 0,1,2). Six-frame search uses both strands.
    // Open reading frame, https://en.wikipedia.org/wiki/Open_reading_frame (accessed 2026-06-14).
    private const int FramesPerStrand = 3;

    /// <summary>
    /// Finds potential Open Reading Frames (ORFs) in all six reading frames (three on the
    /// forward strand and three on the reverse complement).
    /// </summary>
    /// <remarks>
    /// An ORF starts at a start codon (ATG) and ends at the first in-frame stop codon
    /// (TAA, TAG, or TGA). Following the canonical ORF definition (Rosalind "Open Reading
    /// Frames"; NCBI ORFfinder with "ATG only"), <b>every</b> in-frame ATG that is terminated
    /// by a downstream in-frame stop is reported — including nested/overlapping ORFs that share
    /// the same stop codon. A reading begun at an ATG that has no downstream in-frame stop is
    /// not a complete ORF and is not reported. The reported <see cref="OrfInfo.Sequence"/>
    /// spans the start codon through the stop codon inclusive, so its length is divisible by 3;
    /// the encoded protein candidate is the translation up to (excluding) the stop.
    /// </remarks>
    /// <param name="sequence">The DNA sequence to scan.</param>
    /// <param name="minLength">
    /// Minimum ORF length in nucleotides (inclusive lower bound), matching NCBI ORFfinder's
    /// nucleotide length filter. ORFs with length &lt; <paramref name="minLength"/> are excluded.
    /// </param>
    public static IEnumerable<OrfInfo> FindOpenReadingFrames(DnaSequence sequence, int minLength = 100)
    {
        if (sequence is null)
        {
            throw new ArgumentNullException(nameof(sequence));
        }

        string seq = sequence.Sequence;
        for (int frame = 0; frame < FramesPerStrand; frame++)
        {
            foreach (var orf in FindOrfsInFrame(seq, frame, minLength, isReverseComplement: false))
            {
                yield return orf;
            }
        }

        string revComp = sequence.ReverseComplement().Sequence;
        for (int frame = 0; frame < FramesPerStrand; frame++)
        {
            foreach (var orf in FindOrfsInFrame(revComp, frame, minLength, isReverseComplement: true))
            {
                yield return orf;
            }
        }
    }

    private static IEnumerable<OrfInfo> FindOrfsInFrame(
        string seq, int frame, int minLength, bool isReverseComplement)
    {
        // Scan codon positions in this frame. At each ATG, find the first in-frame stop
        // downstream; report that ATG→stop span. Every ATG is considered independently, so
        // nested ORFs sharing a stop are all reported (canonical Rosalind ORF semantics).
        for (int start = frame; start <= seq.Length - CodonLength; start += CodonLength)
        {
            if (!IsCodon(seq, start, StartCodon))
            {
                continue;
            }

            for (int i = start; i <= seq.Length - CodonLength; i += CodonLength)
            {
                if (!IsStopCodon(seq, i))
                {
                    continue;
                }

                int length = i + CodonLength - start;
                if (length >= minLength)
                {
                    yield return new OrfInfo(
                        seq.Substring(start, length),
                        start,
                        frame + 1,
                        isReverseComplement);
                }

                break; // ORF ends at the first in-frame stop.
            }
        }
    }

    private static bool IsCodon(string seq, int index, string codon) =>
        string.CompareOrdinal(seq, index, codon, 0, CodonLength) == 0;

    private static bool IsStopCodon(string seq, int index)
    {
        foreach (string stop in StopCodons)
        {
            if (IsCodon(seq, index, stop))
            {
                return true;
            }
        }

        return false;
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
/// Information about an Open Reading Frame.
/// </summary>
public readonly struct OrfInfo
{
    public OrfInfo(string sequence, int position, int frame, bool isReverseComplement)
    {
        Sequence = sequence;
        Position = position;
        Frame = frame;
        IsReverseComplement = isReverseComplement;
    }

    public string Sequence { get; }
    public int Position { get; }
    public int Frame { get; }
    public bool IsReverseComplement { get; }
    public int Length => Sequence.Length;
    public int CodonCount => Sequence.Length / 3;
}

#endregion
