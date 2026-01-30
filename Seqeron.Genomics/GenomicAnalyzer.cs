using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics
{
    /// <summary>
    /// Genomic analysis algorithms based on suffix trees.
    /// Provides efficient solutions for common bioinformatics problems.
    /// </summary>
    public static class GenomicAnalyzer
    {
        #region Repeat Finding

        /// <summary>
        /// Finds the longest repeated region in a DNA sequence.
        /// Useful for detecting tandem repeats, transposable elements, etc.
        /// Time complexity: O(n) using suffix tree.
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
        /// Finds all repeated substrings of at least the specified length.
        /// </summary>
        public static IEnumerable<RepeatInfo> FindRepeats(DnaSequence sequence, int minLength)
        {
            var tree = sequence.SuffixTree;
            var found = new HashSet<string>();

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

                if (lcpLen >= minLength)
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
        /// Searches for a set of known motifs in a sequence.
        /// </summary>
        public static Dictionary<string, IReadOnlyList<int>> FindKnownMotifs(
            DnaSequence sequence,
            IEnumerable<string> motifs)
        {
            var result = new Dictionary<string, IReadOnlyList<int>>();
            var tree = sequence.SuffixTree;

            foreach (var motif in motifs)
            {
                string normalized = motif.ToUpperInvariant();
                var positions = tree.FindAllOccurrences(normalized);
                if (positions.Count > 0)
                {
                    result[normalized] = positions;
                }
            }

            return result;
        }

        #endregion

        #region Sequence Comparison

        /// <summary>
        /// Finds the longest common subsequence between two DNA sequences.
        /// Useful for identifying conserved regions, gene homology, etc.
        /// Time complexity: O(n + m) using suffix tree.
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
        /// Finds all common regions of at least the specified length.
        /// </summary>
        public static IEnumerable<CommonRegion> FindCommonRegions(
            DnaSequence sequence1,
            DnaSequence sequence2,
            int minLength)
        {
            var tree = sequence1.SuffixTree;
            string seq2 = sequence2.Sequence;
            var found = new HashSet<string>();

            // Slide through sequence2 and find matches
            for (int i = 0; i < seq2.Length - minLength + 1; i++)
            {
                // Binary search for longest match at this position
                int lo = minLength, hi = seq2.Length - i;
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

        /// <summary>
        /// Calculates sequence similarity as percentage of matching k-mers.
        /// </summary>
        public static double CalculateSimilarity(DnaSequence sequence1, DnaSequence sequence2, int kmerSize = 5)
        {
            var kmers1 = GetKmers(sequence1.Sequence, kmerSize);
            var kmers2 = GetKmers(sequence2.Sequence, kmerSize);

            int intersection = kmers1.Intersect(kmers2).Count();
            int union = kmers1.Union(kmers2).Count();

            return union == 0 ? 0 : (double)intersection / union * 100;
        }

        #endregion

        #region Open Reading Frames

        /// <summary>
        /// Finds potential Open Reading Frames (ORFs) - regions that could encode proteins.
        /// An ORF starts with ATG (start codon) and ends with TAA, TAG, or TGA (stop codons).
        /// </summary>
        public static IEnumerable<OrfInfo> FindOpenReadingFrames(DnaSequence sequence, int minLength = 100)
        {
            string seq = sequence.Sequence;
            var startCodon = "ATG";
            var stopCodons = new[] { "TAA", "TAG", "TGA" };

            // Check all 3 reading frames on forward strand
            for (int frame = 0; frame < 3; frame++)
            {
                foreach (var orf in FindOrfsInFrame(seq, frame, startCodon, stopCodons, minLength, false))
                {
                    yield return orf;
                }
            }

            // Check all 3 reading frames on reverse complement
            string revComp = sequence.ReverseComplement().Sequence;
            for (int frame = 0; frame < 3; frame++)
            {
                foreach (var orf in FindOrfsInFrame(revComp, frame, startCodon, stopCodons, minLength, true))
                {
                    yield return orf;
                }
            }
        }

        private static IEnumerable<OrfInfo> FindOrfsInFrame(
            string seq, int frame, string startCodon, string[] stopCodons, int minLength, bool isReverseComplement)
        {
            int? orfStart = null;

            for (int i = frame; i <= seq.Length - 3; i += 3)
            {
                string codon = seq.Substring(i, 3);

                if (orfStart == null && codon == startCodon)
                {
                    orfStart = i;
                }
                else if (orfStart != null && stopCodons.Contains(codon))
                {
                    int length = i + 3 - orfStart.Value;
                    if (length >= minLength)
                    {
                        yield return new OrfInfo(
                            seq.Substring(orfStart.Value, length),
                            orfStart.Value,
                            frame + 1,
                            isReverseComplement);
                    }
                    orfStart = null;
                }
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
}
