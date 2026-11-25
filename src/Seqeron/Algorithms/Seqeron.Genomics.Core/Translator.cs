using System.Text;
using Seqeron.Genomics.Infrastructure;

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

        // Signed 1-based frame label of a 0-based strand offset: +1..+3 on the forward strand,
        // −1..−3 on the reverse complement (Biopython / EMBOSS getorf, transeq numbering).
        private static int FrameLabel(int offset, bool isReverseComplement) =>
            isReverseComplement ? -(offset + 1) : offset + 1;

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
            return TranslateSequence(SequenceExtensions.ToUpperAscii(sequence ?? string.Empty),
                geneticCode ?? GeneticCode.Standard, frame, toFirstStop);
        }

        /// <summary>
        /// Translates a sequence string to protein, optionally accepting the nucleotide symbol
        /// <c>X</c> with Biopython 1.88 semantics.
        /// </summary>
        /// <remarks>
        /// With <paramref name="allowNucleotideX"/> = <see langword="false"/> this is identical to
        /// <see cref="Translate(string, GeneticCode?, int, bool)"/>. With <see langword="true"/>, every
        /// codon is resolved by <see cref="GeneticCode.Translate(string, bool)"/>: <c>X</c>/<c>x</c>
        /// stands for any base (Biopython <c>ambiguous_dna_values['X'] == 'GATC'</c>), so e.g.
        /// "GCXATG" → "AM", while a possible-stop codon containing X (e.g. TAX, XXX) throws
        /// <see cref="ArgumentException"/> like Biopython's <c>TranslationError</c> — unless translation
        /// already ended at an earlier stop with <paramref name="toFirstStop"/>.
        /// </remarks>
        /// <param name="sequence">The DNA or RNA sequence string.</param>
        /// <param name="geneticCode">The genetic code to use (null: Standard).</param>
        /// <param name="frame">Reading frame (0, 1, or 2).</param>
        /// <param name="toFirstStop">Stop translation at first stop codon.</param>
        /// <param name="allowNucleotideX">Accept nucleotide X as in Biopython (opt-in; see remarks).</param>
        /// <returns>The translated protein sequence.</returns>
        /// <exception cref="ArgumentException"><paramref name="toFirstStop"/> is true and the genetic code has
        /// dual-coding stop codons (tables 27, 28, 31), the sequence contains a non-IUPAC codon (other than
        /// X when allowed), or an X-containing possible-stop codon is reached.</exception>
        public static ProteinSequence Translate(string sequence, GeneticCode? geneticCode, int frame,
            bool toFirstStop, bool allowNucleotideX)
        {
            return TranslateSequence(SequenceExtensions.ToUpperAscii(sequence ?? string.Empty),
                geneticCode ?? GeneticCode.Standard, frame, toFirstStop, allowNucleotideX);
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
            var code = RequireOrfGeneticCode(geneticCode);
            return FindOrfsCore(dna.Sequence, code, minLength, searchBothStrands);
        }

        /// <summary>
        /// Finds START→STOP ORFs (getorf <c>-find 1</c> model) in a nucleotide string that may contain
        /// IUPAC ambiguity codes — e.g. a genome with N runs — exactly like
        /// <see cref="FindOrfs(DnaSequence, GeneticCode?, int, bool)"/> on an A/C/G/T sequence.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Input: case-insensitive IUPAC DNA or RNA (A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N; U is
        /// read as T); coordinates and the reverse strand are those of the upper-case DNA spelling
        /// (reverse complement by the canonical IUPAC complement,
        /// <see cref="DnaSequence.GetReverseComplementString(string)"/>).
        /// </para>
        /// <para>
        /// Ambiguous codons follow the getorf port rules (see
        /// <see cref="FindOrfs(string, OrfFindMode, GeneticCode?, int, int, int, bool, bool, bool)"/>): a
        /// codon is a STOP only when every concrete codon it stands for is a stop (e.g. TAR, TRA in table 1);
        /// it is a START only when it is an unambiguous start codon (EMBOSS <c>ajTrnCodonstrTypeC</c>: the
        /// Starts table holds only A/C/G/T codons, so YTG never opens an ORF); residues are resolved by
        /// <see cref="GeneticCode.Translate(string)"/> (Biopython ambiguous tables: GCN → A, RAY → B,
        /// NNN → X).
        /// </para>
        /// </remarks>
        /// <param name="sequence">IUPAC nucleotide sequence (DNA or RNA, case-insensitive).</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <param name="minLength">Minimum ORF length in amino acids (default: 100).</param>
        /// <param name="searchBothStrands">Search both forward and reverse complement strands.</param>
        /// <returns>Enumerable of ORF results.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="sequence"/> contains a non-IUPAC symbol, or the
        /// genetic code has dual-coding stop codons (tables 27, 28, 31).</exception>
        public static IEnumerable<OrfResult> FindOrfs(string sequence, GeneticCode? geneticCode = null,
            int minLength = 100, bool searchBothStrands = true)
        {
            string normalized = SequenceExtensions.NormalizeIupacDna(sequence);
            var code = RequireOrfGeneticCode(geneticCode);
            return FindOrfsCore(normalized, code, minLength, searchBothStrands);
        }

        private static IEnumerable<OrfResult> FindOrfsCore(string sequence, GeneticCode code, int minLength, bool searchBothStrands)
        {
            // Search forward strand in all three frames
            foreach (var orf in FindOrfsInSequence(sequence, code, minLength, false))
                yield return orf;

            // Search reverse complement strand
            if (searchBothStrands)
            {
                foreach (var orf in FindOrfsInSequence(DnaSequence.GetReverseComplementString(sequence), code, minLength, true))
                    yield return orf;
            }
        }

        // ORF finding needs codons that unambiguously terminate: tables with dual-coding stops are rejected.
        /// <summary>The genetic code of an ORF search (<c>null</c> = Standard); the tables with dual-coding stop codons are refused.</summary>
        public static GeneticCode RequireOrfGeneticCode(GeneticCode? geneticCode)
        {
            var code = geneticCode ?? GeneticCode.Standard;
            if (HasDualCodingStopCodons(code))
                throw new ArgumentException(
                    $"ORF finding cannot be used with genetic code table {code.TableNumber} " +
                    "because its stop codons also code for an amino acid (context-dependent termination).",
                    nameof(geneticCode));
            return code;
        }

        // getorf codon classification (EMBOSS ajTrnCodonstrTypeC): START only when the Starts table marks
        // the codon 'M' — EMBOSS fills it for the 64 A/C/G/T codons only, so an ambiguous codon is never a
        // START (verified with the getorf 6.6.0 binary: YTG does not open an ORF in table 1).
        /// <summary>Whether <paramref name="codon"/> is an unambiguous start codon of <paramref name="code"/>.</summary>
        public static bool IsOrfStartCodon(GeneticCode code, string codon) =>
            IsUnambiguousDnaCodon(codon) && code.IsStartCodon(codon);

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

            return TranslateSixFramesCore(dna.Sequence, dna.ReverseComplement().Sequence,
                geneticCode ?? GeneticCode.Standard);
        }

        /// <summary>
        /// Translates all six reading frames of a nucleotide string that may contain IUPAC ambiguity
        /// codes (e.g. a genome with N runs).
        /// </summary>
        /// <remarks>
        /// Biopython semantics: frame +k is <c>Seq(s)[k-1:].translate(table)</c> and frame −k is
        /// <c>Seq(s).reverse_complement()[k-1:].translate(table)</c> (the
        /// <see cref="TranslateSixFrames(DnaSequence, GeneticCode?)"/> numbering; trailing bases that do not
        /// form a codon are dropped). Input is case-insensitive IUPAC DNA or RNA (A, C, G, T, U, R, Y, S,
        /// W, K, M, B, D, H, V, N; U is read as T); the reverse strand is the canonical IUPAC reverse
        /// complement (<see cref="DnaSequence.GetReverseComplementString(string)"/>). Each codon is
        /// resolved by <see cref="GeneticCode.Translate(string)"/> — Biopython's ambiguous codon tables:
        /// GCN → A, RAY → B, TAR → '*', a possible stop such as TAN → X, NNN → X; dual-coding codons of
        /// tables 27, 28, 31 read as their amino acid, as in Biopython <c>translate()</c>. Internal stops
        /// are rendered as '*'.
        /// </remarks>
        /// <param name="sequence">IUPAC nucleotide sequence (DNA or RNA, case-insensitive).</param>
        /// <param name="geneticCode">The genetic code to use (default: Standard).</param>
        /// <returns>Dictionary with frame keys (-3 to +3, excluding 0) and protein values.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="sequence"/> contains a non-IUPAC symbol.</exception>
        public static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFrames(string sequence,
            GeneticCode? geneticCode = null)
        {
            string forward = SequenceExtensions.NormalizeIupacDna(sequence);
            return TranslateSixFramesCore(forward, DnaSequence.GetReverseComplementString(forward),
                geneticCode ?? GeneticCode.Standard);
        }

        /// <summary>The six frames of an already normalised strand pair under the Biopython convention.</summary>
        public static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFramesCore(string forward,
            string reverseComplement, GeneticCode code)
        {
            var result = new Dictionary<int, ProteinSequence>();
            for (int offset = 0; offset < ReadingFramesPerStrand; offset++)
            {
                // Forward frame numbers are 1-based: offset 0 -> +1, etc.
                result[FrameLabel(offset, isReverseComplement: false)] = TranslateSequence(forward, code, offset, false);

                // Reverse complement: frame -k = reverse complement at offset (k-1).
                result[FrameLabel(offset, isReverseComplement: true)] = TranslateSequence(reverseComplement, code, offset, false);
            }

            return result;
        }

        /// <summary>
        /// Splits a coding sequence into its complete, non-overlapping in-frame triplets, starting at
        /// offset <paramref name="frame"/>: element <c>k</c> is the triplet at nucleotide position
        /// <c>frame + 3k</c>, in upper-case DNA spelling. Input is case-insensitive DNA or RNA (U is
        /// read as T; CodonW <c>ident_codon</c>, EMBOSS <c>ajBaseAlphaToBin</c>). A triplet containing
        /// any symbol other than A/C/G/T (IUPAC ambiguity codes, N, gaps, …) is returned as
        /// <see langword="null"/> so that callers skip it without shifting the frame (EMBOSS
        /// <c>ajCodSetTripletsS</c>: "Skips triplets with ambiguity codes and any incomplete triplet
        /// at the end"); a trailing partial triplet is dropped.
        /// </summary>
        /// <remarks>
        /// This is the single codon-splitting core of the library, placed in Core so that every layer
        /// can call it (<c>CodonUsageAnalyzer</c>, <c>CodonOptimizer</c>, and — cross-batch —
        /// <c>SequenceStatistics.CalculateCodonFrequencies</c>, <c>GenomeAnnotator.GetCodonUsage</c>).
        /// </remarks>
        /// <param name="sequence">Coding sequence; null/empty yields an empty array.</param>
        /// <param name="frame">0-based offset of the first codon (0, 1, 2; larger values are plain
        /// offsets, as EMBOSS compseq <c>-frame</c>).</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="frame"/> is negative.</exception>
        public static string?[] SplitInFrameCodons(string? sequence, int frame = 0)
        {
            ArgumentOutOfRangeException.ThrowIfNegative(frame);
            if (string.IsNullOrEmpty(sequence) || sequence.Length - frame < CodonLength)
                return Array.Empty<string?>();

            string normalized = Infrastructure.SequenceCase.ToUpperDna(sequence);
            var codons = new string?[(normalized.Length - frame) / CodonLength];
            for (int k = 0; k < codons.Length; k++)
            {
                string codon = normalized.Substring(frame + CodonLength * k, CodonLength);
                codons[k] = IsUnambiguousDnaCodon(codon) ? codon : null;
            }

            return codons;
        }

        /// <summary>
        /// Counts the complete, unambiguous in-frame codons of a coding sequence (upper-case DNA keys;
        /// codons that do not occur are absent). Same codon set as <see cref="SplitInFrameCodons"/>:
        /// case-insensitive DNA or RNA, triplets with other symbols skipped without shifting the frame,
        /// trailing partial triplet ignored (EMBOSS <c>cusp</c> / <c>ajCodSetTripletsS</c>).
        /// </summary>
        /// <param name="sequence">Coding sequence; null/empty yields an empty dictionary.</param>
        /// <param name="frame">0-based offset of the first codon (see <see cref="SplitInFrameCodons"/>).</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="frame"/> is negative.</exception>
        public static Dictionary<string, int> CountCodons(string? sequence, int frame = 0)
        {
            var counts = new Dictionary<string, int>();
            foreach (string? codon in SplitInFrameCodons(sequence, frame))
            {
                if (codon is null)
                    continue;
                counts[codon] = counts.GetValueOrDefault(codon) + 1;
            }

            return counts;
        }

        private static bool IsUnambiguousDnaCodon(string codon)
        {
            foreach (char c in codon)
            {
                if (c is not ('A' or 'C' or 'G' or 'T'))
                    return false;
            }
            return true;
        }

        /// <summary>Translates <paramref name="sequence"/> from offset <paramref name="frame"/> (the shared core of the translation entry points; used by the EMBOSS ports).</summary>
        public static ProteinSequence TranslateSequence(string sequence, GeneticCode geneticCode,
            int frame, bool toFirstStop, bool allowNucleotideX = false)
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

            // GeneticCode.Translate normalises case and T/U itself (single normalisation point).
            var sb = new StringBuilder();

            // Trailing nucleotides that cannot form a full codon are ignored
            // (Biopython six_frame_translations: fragment_length = 3*((len-i)//3)).
            for (int i = frame; i + CodonLength <= sequence.Length; i += CodonLength)
            {
                string codon = sequence.Substring(i, CodonLength);
                char aa = geneticCode.Translate(codon, allowNucleotideX);

                if (toFirstStop && aa == '*')
                    break;

                sb.Append(aa);
            }

            return new ProteinSequence(sb.ToString());
        }

        // A stop codon that the table also translates as an amino acid (NCBI gc.prt tables
        // 27, 28, 31; Biopython "dual_coding" check in Bio.Seq._translate_str).
        internal static bool HasDualCodingStopCodons(GeneticCode geneticCode)
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
            // ORF = region from a START codon to a STOP codon (EMBOSS getorf -find 1: "a region that begins with a
            // START codon and ends with a STOP codon"; nested in-frame STARTs do not open a new ORF), scanned in all
            // three frames by the shared OrfScanner core (also used by GenomicAnalyzer.FindOpenReadingFrames).
            // An ORF that runs off the end of the strand is reported up to the last complete codon (getorf
            // WriteORF(start, pos+2)); trailing partial-codon bases are not part of the ORF.
            // GeneticCode.Translate / IsStartCodon normalise T/U themselves.
            foreach (var region in OrfScanner.FindRegions(sequence, geneticCode,
                         codon => IsOrfStartCodon(geneticCode, codon), nestedStarts: false, includeOpenEnded: true))
            {
                int codingEnd = region.IsTerminated ? region.End - CodonLength : region.End;
                var protein = new StringBuilder();
                // Initiator is read as Met whatever the start codon (EMBOSS getorf
                // -methionine default Y; Biopython translate(cds=True); NCBI gc).
                protein.Append(InitiatorMethionine);
                for (int i = region.Start + CodonLength; i < codingEnd; i += CodonLength)
                    protein.Append(geneticCode.Translate(sequence.Substring(i, CodonLength)));

                if (protein.Length >= minLength)
                {
                    yield return new OrfResult(
                        region.Start,
                        // Inclusive end = last base of the stop codon (INSDC feature table: CDS "location
                        // includes stop codon"; getorf prints the range without the STOP: WriteORF(start, pos-1)).
                        region.End - 1,
                        FrameLabel(region.FrameOffset, isReverseComplement),
                        new ProteinSequence(protein.ToString()));
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

        /// <summary>
        /// Gets whether the ORF ends with a STOP codon (its stop-inclusive range is one codon longer
        /// than the protein); false for an open ORF that runs to the end of the strand.
        /// </summary>
        public bool IsTerminated => NucleotideLength > Protein.Length * 3;

    }

}
