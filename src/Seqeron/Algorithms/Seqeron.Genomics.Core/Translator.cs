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

        // The initiator residue of an ORF (NCBI The Genetic Codes: the initiator codon is by
        // default translated as methionine).
        private const char InitiatorMethionine = 'M';
        /// <summary>
        /// Translates a DNA sequence to protein using the specified genetic code.
        /// </summary>
        /// <param name="dna">The DNA sequence to translate.</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <param name="frame">Reading frame (0, 1, or 2).</param>
        /// <param name="toFirstStop">Stop translation at first stop codon.</param>
        /// <returns>The translated protein sequence.</returns>
        /// <exception cref="ArgumentException"><paramref name="toFirstStop"/> is true and the genetic code has
        /// dual-coding stop codons (tables 27, 28, 31), or the sequence contains a non-IUPAC codon.</exception>
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
        /// <exception cref="ArgumentException"><paramref name="toFirstStop"/> is true and the genetic code has
        /// dual-coding stop codons (tables 27, 28, 31), or the sequence contains a non-IUPAC codon.</exception>
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
        /// <exception cref="ArgumentException"><paramref name="toFirstStop"/> is true and the genetic code has
        /// dual-coding stop codons (tables 27, 28, 31), or the sequence contains a non-IUPAC codon.</exception>
        public static ProteinSequence Translate(string sequence, GeneticCode? geneticCode = null,
            int frame = 0, bool toFirstStop = false)
        {
            // null/empty yield an empty protein, but the frame / toFirstStop arguments are still
            // validated exactly as for the DnaSequence / RnaSequence overloads.
            return TranslateSequence((sequence ?? string.Empty).ToUpperInvariant(),
                geneticCode ?? GeneticCode.Standard, frame, toFirstStop);
        }

        /// <summary>
        /// Finds all Open Reading Frames (ORFs) in a DNA sequence.
        /// An ORF starts with a start codon and ends with a stop codon.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Model: EMBOSS getorf <c>-find 1</c> (START→STOP, linear sequence): in each of the three
        /// frames of a strand, scanning opens an ORF at the first START codon (only when no ORF is
        /// open, so nested in-frame starts are not reported separately) and closes it at the next
        /// in-frame STOP codon. An ORF that reaches the end of the strand without a STOP is reported
        /// as an open (3'-incomplete) ORF.
        /// </para>
        /// <para>
        /// The initiator residue is always reported as Met ('M'), whatever the start codon
        /// (e.g. TTG, CTG, GTG) — EMBOSS getorf <c>-methionine</c> (default Y: "Change initial
        /// START codons to Methionine"); Biopython <c>translate(cds=True)</c>; NCBI The Genetic
        /// Codes ("The initiator codon ... is by default translated as methionine").
        /// </para>
        /// <para>
        /// Coordinates are 0-based and inclusive, in the scanned strand's coordinates (reverse-strand
        /// ORFs: coordinates of the reverse complement, frames −1..−3). For a terminated ORF,
        /// <see cref="OrfResult.EndPosition"/> is the last base of the STOP codon — the INSDC
        /// feature-table CDS convention ("location includes stop codon"; also what Biopython
        /// <c>translate(cds=True)</c> expects); getorf itself prints the range without the STOP,
        /// i.e. its end = <c>EndPosition − 3</c> (+1 for 1-based). For an open ORF, EndPosition is the
        /// last base of the last complete codon (getorf <c>WriteORF(start, pos+2)</c>), so
        /// <see cref="OrfResult.NucleotideLength"/> is always a multiple of three.
        /// </para>
        /// <para>
        /// Stop codons are the codons the table translates as '*'. Tables with dual-coding
        /// (context-dependent) stop codons — NCBI 27, 28, 31 — are rejected: in each of them every
        /// stop codon also codes for an amino acid, so no codon unambiguously ends an ORF. This
        /// mirrors Biopython, which refuses "translate to the first stop" (<c>to_stop=True</c>) for
        /// such tables, and <see cref="Translate(DnaSequence, GeneticCode?, int, bool)"/> with
        /// <c>toFirstStop</c>.
        /// </para>
        /// </remarks>
        /// <param name="dna">The DNA sequence to search.</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <param name="minLength">Minimum ORF length in amino acids (default: 100).</param>
        /// <param name="searchBothStrands">Search both forward and reverse complement strands.</param>
        /// <returns>Enumerable of ORF results.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="dna"/> is null.</exception>
        /// <exception cref="ArgumentException">The genetic code has dual-coding stop codons (tables 27, 28, 31).</exception>
        public static IEnumerable<OrfResult> FindOrfs(DnaSequence dna, GeneticCode? geneticCode = null,
            int minLength = 100, bool searchBothStrands = true)
        {
            ArgumentNullException.ThrowIfNull(dna);
            var code = geneticCode ?? GeneticCode.Standard;
            if (HasDualCodingStopCodons(code))
                throw new ArgumentException(
                    $"ORF finding cannot be used with genetic code table {code.TableNumber} " +
                    "because its stop codons also code for an amino acid (context-dependent termination).",
                    nameof(geneticCode));
            return FindOrfsCore(dna, code, minLength, searchBothStrands);
        }

        private static IEnumerable<OrfResult> FindOrfsCore(DnaSequence dna, GeneticCode code, int minLength, bool searchBothStrands)
        {

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

            // "Translate to the first stop" is undefined when the table has dual-coding
            // (context-dependent) stop codons, e.g. NCBI tables 27, 28, 31: such codons are
            // translated as their amino acid, so a genuine terminator would be read through.
            // Biopython Bio.Seq._translate_str raises ValueError for to_stop=True here.
            if (toFirstStop && HasDualCodingStopCodons(geneticCode))
                throw new ArgumentException(
                    $"toFirstStop cannot be used with genetic code table {geneticCode.TableNumber} " +
                    "because it contains codons that code for both STOP and an amino acid.",
                    nameof(toFirstStop));

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

        // A stop codon that the table also translates as an amino acid (NCBI gc.prt tables
        // 27, 28, 31; Biopython "dual_coding" check in Bio.Seq._translate_str).
        private static bool HasDualCodingStopCodons(GeneticCode geneticCode)
        {
            foreach (var stop in geneticCode.StopCodons)
            {
                if (geneticCode.CodonTable.TryGetValue(stop, out char aa) && aa != '*')
                    return true;
            }
            return false;
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
                            // Initiator is read as Met whatever the start codon (EMBOSS getorf
                            // -methionine default Y; Biopython translate(cds=True); NCBI gc).
                            currentProtein.Append(InitiatorMethionine);
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
                                    // Inclusive end = last base of the stop codon (INSDC
                                    // feature table: CDS "location includes stop codon").
                                    // getorf prints the range without the STOP (WriteORF(start, pos-1)).
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

                // ORF that runs off the end of the strand: it ends at the last base of the last
                // complete codon (EMBOSS getorf WriteORF(start, pos+2)); trailing partial-codon
                // bases are not part of the ORF.
                if (currentOrfStart != null && currentProtein.Length >= minLength)
                {
                    yield return new OrfResult(
                        currentOrfStart.Value,
                        currentOrfStart.Value + currentProtein.Length * CodonLength - 1,
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
