using System;
using System.Collections.Generic;
using System.Text;

namespace Seqeron.Genomics
{
    /// <summary>
    /// Translates DNA or RNA sequences to protein sequences.
    /// </summary>
    public static class Translator
    {
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
            if (dna == null)
                throw new ArgumentNullException(nameof(dna));

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
            if (rna == null)
                throw new ArgumentNullException(nameof(rna));

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
            if (dna == null)
                throw new ArgumentNullException(nameof(dna));

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
        /// <param name="dna">The DNA sequence to translate.</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <returns>Dictionary with frame keys (-3 to +3, excluding 0) and protein values.</returns>
        public static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFrames(DnaSequence dna,
            GeneticCode? geneticCode = null)
        {
            if (dna == null)
                throw new ArgumentNullException(nameof(dna));

            var code = geneticCode ?? GeneticCode.Standard;
            var result = new Dictionary<int, ProteinSequence>();

            // Forward strand: frames +1, +2, +3
            result[1] = TranslateSequence(dna.Sequence, code, 0, false);
            result[2] = TranslateSequence(dna.Sequence, code, 1, false);
            result[3] = TranslateSequence(dna.Sequence, code, 2, false);

            // Reverse complement: frames -1, -2, -3
            var revComp = dna.ReverseComplement();
            result[-1] = TranslateSequence(revComp.Sequence, code, 0, false);
            result[-2] = TranslateSequence(revComp.Sequence, code, 1, false);
            result[-3] = TranslateSequence(revComp.Sequence, code, 2, false);

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

            for (int i = frame; i + 3 <= rnaSequence.Length; i += 3)
            {
                string codon = rnaSequence.Substring(i, 3);
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

            for (int frame = 0; frame < 3; frame++)
            {
                int? currentOrfStart = null;
                var currentProtein = new StringBuilder();

                for (int i = frame; i + 3 <= rnaSequence.Length; i += 3)
                {
                    string codon = rnaSequence.Substring(i, 3);
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
                                    i + 2,
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
