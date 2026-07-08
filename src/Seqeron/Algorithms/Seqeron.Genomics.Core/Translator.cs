using System.Text;

namespace Seqeron.Genomics.Core
{
    /// <summary>
    /// Translates DNA or RNA sequences to protein sequences.
    /// </summary>
    public static class Translator
    {
        // A codon is exactly three nucleotides (NCBI genetic code tables).
        private const int CodonLength = 3;

        // A double-stranded sequence has three forward reading frames at
        // offsets 0, 1, 2 (Biopython six_frame_translations; EMBOSS transeq).
        private const int ReadingFramesPerStrand = 3;
        /// <summary>
        /// Translates a DNA sequence to protein using the specified genetic code.
        /// </summary>
        /// <param name="dna">The DNA sequence to translate.</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <param name="frame">Reading frame (0, 1, or 2).</param>
        /// <param name="toFirstStop">Stop translation at first stop codon.</param>
        /// <returns>The translated protein sequence.</returns>
        public static ProteinSequence Translate(DnaSequence dna, GeneticCode? geneticCode = null,
            int frame = 0, bool toFirstStop = false)
        {
            ArgumentNullException.ThrowIfNull(dna);

            return TranslateSequence(dna.Sequence, geneticCode ?? GeneticCode.Standard, frame, toFirstStop);
        }

        /// <summary>
        /// Translates an RNA sequence to protein using the specified genetic code.
        /// </summary>
        /// <param name="rna">The RNA sequence to translate.</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <param name="frame">Reading frame (0, 1, or 2).</param>
        /// <param name="toFirstStop">Stop translation at first stop codon.</param>
        /// <returns>The translated protein sequence.</returns>
        public static ProteinSequence Translate(RnaSequence rna, GeneticCode? geneticCode = null,
            int frame = 0, bool toFirstStop = false)
        {
            ArgumentNullException.ThrowIfNull(rna);

            return TranslateSequence(rna.Sequence, geneticCode ?? GeneticCode.Standard, frame, toFirstStop);
        }

        /// <summary>
        /// Translates a sequence string to protein.
        /// </summary>
        /// <param name="sequence">The DNA or RNA sequence string.</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <param name="frame">Reading frame (0, 1, or 2).</param>
        /// <param name="toFirstStop">Stop translation at first stop codon.</param>
        /// <returns>The translated protein sequence.</returns>
        public static ProteinSequence Translate(string sequence, GeneticCode? geneticCode = null,
            int frame = 0, bool toFirstStop = false)
        {
            if (string.IsNullOrEmpty(sequence))
                return new ProteinSequence("");

            return TranslateSequence(sequence.ToUpperInvariant(), geneticCode ?? GeneticCode.Standard, frame, toFirstStop);
        }

        /// <summary>
        /// Finds all Open Reading Frames (ORFs) in a DNA sequence.
        /// An ORF starts with a start codon and ends with a stop codon.
        /// </summary>
        /// <param name="dna">The DNA sequence to search.</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <param name="minLength">Minimum ORF length in amino acids (default: 100).</param>
        /// <param name="searchBothStrands">Search both forward and reverse complement strands.</param>
        /// <returns>Enumerable of ORF results.</returns>
        public static IEnumerable<OrfResult> FindOrfs(DnaSequence dna, GeneticCode? geneticCode = null,
            int minLength = 100, bool searchBothStrands = true)
        {
            ArgumentNullException.ThrowIfNull(dna);
            return FindOrfsCore(dna, geneticCode, minLength, searchBothStrands);
        }

        private static IEnumerable<OrfResult> FindOrfsCore(DnaSequence dna, GeneticCode? geneticCode, int minLength, bool searchBothStrands)
        {
            var code = geneticCode ?? GeneticCode.Standard;

            // Search forward strand in all three frames
            foreach (var orf in FindOrfsInSequence(dna.Sequence, code, minLength, false))
                yield return orf;

            // Search reverse complement strand
            if (searchBothStrands)
            {
                var revComp = dna.ReverseComplement();
                foreach (var orf in FindOrfsInSequence(revComp.Sequence, code, minLength, true))
                    yield return orf;
            }
        }

        /// <summary>
        /// Translates all six reading frames of a DNA sequence.
        /// </summary>
        /// <remarks>
        /// Reverse-frame numbering follows the Biopython
        /// <c>SeqUtils.six_frame_translations</c> convention: frame -k is the
        /// translation of the reverse complement read at offset (k-1), i.e.
        /// <c>frames[-(i+1)] = translate(reverse_complement(seq)[i:])</c>. This
        /// is the "alternative" convention explicitly documented by EMBOSS transeq
        /// (frame -1 = frame 1 of the reverse complement), as opposed to the
        /// EMBOSS phase-locked default. Frames render internal stop codons as '*'
        /// (translation is not terminated early).
        /// </remarks>
        /// <param name="dna">The DNA sequence to translate.</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <returns>Dictionary with frame keys (-3 to +3, excluding 0) and protein values.</returns>
        public static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFrames(DnaSequence dna,
            GeneticCode? geneticCode = null)
        {
            ArgumentNullException.ThrowIfNull(dna);

            var code = geneticCode ?? GeneticCode.Standard;
            var result = new Dictionary<int, ProteinSequence>();

            // Forward strand: frames +1, +2, +3 at offsets 0, 1, 2.
            var revComp = dna.ReverseComplement();
            for (int offset = 0; offset < ReadingFramesPerStrand; offset++)
            {
                // Forward frame numbers are 1-based: offset 0 -> +1, etc.
                result[offset + 1] = TranslateSequence(dna.Sequence, code, offset, false);

                // Reverse complement: frame -k = reverse complement at offset (k-1).
                result[-(offset + 1)] = TranslateSequence(revComp.Sequence, code, offset, false);
            }

            return result;
        }

        private static ProteinSequence TranslateSequence(string sequence, GeneticCode geneticCode,
            int frame, bool toFirstStop)
        {
            if (frame < 0 || frame > 2)
                throw new ArgumentOutOfRangeException(nameof(frame), "Frame must be 0, 1, or 2.");

            // Convert T to U for translation
            var rnaSequence = sequence.Replace('T', 'U');
            var sb = new StringBuilder();

            // Trailing nucleotides that cannot form a full codon are ignored
            // (Biopython six_frame_translations: fragment_length = 3*((len-i)//3)).
            for (int i = frame; i + CodonLength <= rnaSequence.Length; i += CodonLength)
            {
                string codon = rnaSequence.Substring(i, CodonLength);
                char aa = geneticCode.Translate(codon);

                if (toFirstStop && aa == '*')
                    break;

                sb.Append(aa);
            }

            return new ProteinSequence(sb.ToString());
        }

        private static IEnumerable<OrfResult> FindOrfsInSequence(string sequence, GeneticCode geneticCode,
            int minLength, bool isReverseComplement)
        {
            var rnaSequence = sequence.Replace('T', 'U');

            // ORF = region from a START codon to a STOP codon
            // (EMBOSS getorf -find 1: "a region that begins with a START codon
            // and ends with a STOP codon"). Scanned in all three frames.
            for (int frame = 0; frame < ReadingFramesPerStrand; frame++)
            {
                int? currentOrfStart = null;
                var currentProtein = new StringBuilder();

                for (int i = frame; i + CodonLength <= rnaSequence.Length; i += CodonLength)
                {
                    string codon = rnaSequence.Substring(i, CodonLength);
                    char aa = geneticCode.Translate(codon);

                    if (currentOrfStart == null)
                    {
                        // Looking for start codon
                        if (geneticCode.IsStartCodon(codon))
                        {
                            currentOrfStart = i;
                            currentProtein.Clear();
                            currentProtein.Append(aa);
                        }
                    }
                    else
                    {
                        // In an ORF, looking for stop
                        if (aa == '*')
                        {
                            // Found stop codon
                            if (currentProtein.Length >= minLength)
                            {
                                yield return new OrfResult(
                                    currentOrfStart.Value,
                                    // Inclusive end = last base of the stop codon
                                    // (EMBOSS getorf positions include the STOP).
                                    i + (CodonLength - 1),
                                    isReverseComplement ? -(frame + 1) : frame + 1,
                                    new ProteinSequence(currentProtein.ToString())
                                );
                            }
                            currentOrfStart = null;
                        }
                        else
                        {
                            currentProtein.Append(aa);
                        }
                    }
                }

                // Handle ORF that extends to end of sequence
                if (currentOrfStart != null && currentProtein.Length >= minLength)
                {
                    yield return new OrfResult(
                        currentOrfStart.Value,
                        rnaSequence.Length - 1,
                        isReverseComplement ? -(frame + 1) : frame + 1,
                        new ProteinSequence(currentProtein.ToString())
                    );
                }
            }
        }
    }

    /// <summary>
    /// Represents an Open Reading Frame (ORF) found in a sequence.
    /// </summary>
    public readonly record struct OrfResult(
        int StartPosition,
        int EndPosition,
        int Frame,
        ProteinSequence Protein)
    {
        /// <summary>
        /// Gets the length of the ORF in nucleotides.
        /// </summary>
        public int NucleotideLength => EndPosition - StartPosition + 1;

        /// <summary>
        /// Gets the length of the ORF in amino acids.
        /// </summary>
        public int AminoAcidLength => Protein.Length;
    }
}
