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
            string normalized = NormalizeIupacNucleotides(sequence);
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
        private static GeneticCode RequireOrfGeneticCode(GeneticCode? geneticCode)
        {
            var code = geneticCode ?? GeneticCode.Standard;
            if (HasDualCodingStopCodons(code))
                throw new ArgumentException(
                    $"ORF finding cannot be used with genetic code table {code.TableNumber} " +
                    "because its stop codons also code for an amino acid (context-dependent termination).",
                    nameof(geneticCode));
            return code;
        }

        // Validates an IUPAC nucleotide string (A C G T U R Y S W K M B D H V N, any case) and returns
        // its upper-case DNA spelling (U → T), the alphabet of the IUPAC-aware string overloads.
        private static string NormalizeIupacNucleotides(string sequence)
        {
            ArgumentNullException.ThrowIfNull(sequence);
            string upper = Infrastructure.SequenceCase.ToUpperDna(sequence);
            for (int i = 0; i < upper.Length; i++)
            {
                if (!IupacHelper.IsNucleotideCode(upper[i]))
                    throw new ArgumentException(
                        $"Invalid nucleotide '{sequence[i]}' at position {i}: expected an IUPAC nucleotide code " +
                        "(A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N).", nameof(sequence));
            }
            return upper;
        }

        // getorf codon classification (EMBOSS ajTrnCodonstrTypeC): START only when the Starts table marks
        // the codon 'M' — EMBOSS fills it for the 64 A/C/G/T codons only, so an ambiguous codon is never a
        // START (verified with the getorf 6.6.0 binary: YTG does not open an ORF in table 1).
        private static bool IsOrfStartCodon(GeneticCode code, string codon) =>
            IsUnambiguousDnaCodon(codon) && code.IsStartCodon(codon);

        // EMBOSS getorf.acd defaults: -minsize 30, -maxsize 1000000, -flanking 100.
        private const int GetOrfDefaultMinSize = 30;
        private const int GetOrfDefaultMaxSize = 1_000_000;
        private const int GetOrfDefaultFlanking = 100;

        /// <summary>
        /// Finds ORFs (or the nucleotides flanking their START / STOP codons) with any of the seven
        /// EMBOSS getorf <c>-find</c> modes, on a linear or (getorf <c>-circular</c>) circular sequence.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Port of EMBOSS getorf <c>getorf_FindORFs</c> / <c>getorf_WriteORF</c> / <c>getorf_AppORF</c>
        /// (emboss/getorf.c). Two ORF definitions (getorf.acd <c>-find</c> help): a region free of
        /// STOP codons (modes 0, 2, 5, 6; on a linear sequence every frame is treated as already open
        /// at the sequence start, so the leading region before the first STOP and the trailing region
        /// after the last STOP are reported too), or a region that begins with a START codon and ends
        /// with a STOP codon (modes 1, 3, 4; nested in-frame STARTs do not open a new ORF). An ORF
        /// that reaches the strand end without a STOP ends at the last base of the last complete
        /// codon (getorf <c>WriteORF(start, pos+2)</c>).
        /// </para>
        /// <para>
        /// Size filter (getorf <c>-minsize</c> / <c>-maxsize</c>, in nucleotides): for the
        /// protein-storing modes 0, 1 and 4 getorf divides both bounds by 3 (integer division) and
        /// compares them with the translated length in residues; for modes 2, 3, 5 and 6 they are
        /// compared with the ORF's nucleotide length (STOP excluded). For the flanking modes the
        /// filter applies to the underlying ORF, and a window that would cross either sequence end
        /// is not reported (getorf.acd <c>-flanking</c>: "no output is given for this codon").
        /// </para>
        /// <para>
        /// Coordinates are 0-based and inclusive, in the scanned strand's coordinates (reverse-strand
        /// results: coordinates of the reverse complement, frames −1..−3), as in
        /// <see cref="FindOrfs(DnaSequence, GeneticCode?, int, bool)"/>: for modes 0–3 the region is
        /// the ORF with a stop-inclusive end (see <see cref="OrfFeature.IsTerminated"/>); for modes
        /// 4–6 it is the flank window <c>[codon − flanking, codon + 2 + flanking]</c>.
        /// <see cref="OrfFeature.Frame"/> is ±(offset of the ORF's first base mod 3, plus 1).
        /// <see cref="OrfFeature.ToGetOrfRange(int)"/> gives getorf's 1-based <c>[s - e]</c> header.
        /// Results are yielded in getorf output order (forward strand, then reverse strand; within a
        /// strand in the order the ORFs close), so the k-th result is getorf's <c>name_k</c>.
        /// </para>
        /// <para>
        /// <paramref name="circular"/> (getorf <c>-circular</c>): as getorf, the sequence of length L is
        /// scanned three times end to end ("append it to itself to triple its length so can deal easily
        /// with wrapped ORFs") and its reverse complement likewise; every frame starts closed (no
        /// leading region), a START at or past position L is ignored, and a frame stops at its first
        /// STOP / end codon at or past L. A reported region that starts at or past L is shifted back by
        /// L (getorf <c>while(s &gt; len)</c>), so <see cref="OrfFeature.StartPosition"/> lies in
        /// [0, L) while <see cref="OrfFeature.EndPosition"/> may reach past L − 1: such an ORF crosses
        /// the origin (getorf "(ORF crosses the breakpoint)"; base i ≥ L is base i mod L), and its
        /// <see cref="OrfFeature.Sequence"/> runs across the origin. Flank windows (modes 4–6) are
        /// taken from the tripled sequence and dropped only when they leave [0, 3L). In modes 0 and 2,
        /// a frame of the tripled sequence that contains no STOP at all is reported, after the
        /// scan, as one open record from offset f to the last base (getorf "full length of a circular
        /// genome" block: start f, end 3L − 1); in mode 0 its 1–2 trailing bases that do not form a
        /// codon are translated as EMBOSS <c>ajTrnSeqOrig</c> does (two bases b1 b2 as b1 b2 N, one base
        /// as X), in mode 2 they are kept.
        /// </para>
        /// <para>
        /// Documented divergences from the getorf 6.6.0 binary: (1) in modes 1, 3, 4 on a linear
        /// sequence, a START that is the last complete codon of its frame forms a one-codon open ORF
        /// that is reported when it passes the size filter; getorf never writes it (its START branch
        /// precedes the end-of-sequence branch, so the ORF is left open when the loop ends). (2) The
        /// circular full-length record of mode 2 holds the frame's nucleotides and is size-filtered by
        /// nucleotide length (getorf writes the protein and compares its residue count with the
        /// nucleotide bounds), and the frame-3 full-length record uses the same <c>≤ maxsize</c> test
        /// as frames 1–2 (getorf.c tests <c>&gt;= maxsize</c> for frame 3, a typo).
        /// </para>
        /// </remarks>
        /// <param name="dna">The DNA sequence to search.</param>
        /// <param name="mode">The getorf <c>-find</c> mode.</param>
        /// <param name="geneticCode">The genetic code (default: Standard, NCBI table 1).</param>
        /// <param name="minSize">getorf <c>-minsize</c>: minimum ORF size in nucleotides (default 30).</param>
        /// <param name="maxSize">getorf <c>-maxsize</c>: maximum ORF size in nucleotides (default 1 000 000).</param>
        /// <param name="flanking">getorf <c>-flanking</c>: bases reported on each side of the codon in modes 4–6 (default 100).</param>
        /// <param name="searchBothStrands">getorf <c>-reverse</c>: also search the reverse complement (default true).</param>
        /// <param name="methionine">getorf <c>-methionine</c>: report the initiator residue as Met in mode 1 (default true).</param>
        /// <param name="circular">getorf <c>-circular</c>: treat the sequence as circular (default false; see remarks).</param>
        /// <returns>The reported ORFs / flank windows in getorf output order.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="dna"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="mode"/> is undefined, or a size / flanking value is negative.</exception>
        /// <exception cref="ArgumentException">The genetic code has dual-coding stop codons (tables 27, 28, 31).</exception>
        public static IEnumerable<OrfFeature> FindOrfs(DnaSequence dna, OrfFindMode mode, GeneticCode? geneticCode = null,
            int minSize = GetOrfDefaultMinSize, int maxSize = GetOrfDefaultMaxSize, int flanking = GetOrfDefaultFlanking,
            bool searchBothStrands = true, bool methionine = true, bool circular = false)
        {
            ArgumentNullException.ThrowIfNull(dna);
            var settings = CreateGetOrfSettings(mode, geneticCode, minSize, maxSize, flanking, methionine);
            return FindOrfsGetOrfCore(dna.Sequence, settings, searchBothStrands, circular);
        }

        /// <summary>
        /// Finds ORFs / flank windows with the EMBOSS getorf <c>-find</c> modes in a nucleotide string
        /// that may contain IUPAC ambiguity codes (e.g. a genome with N runs), as getorf scans it.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Same algorithm, defaults, coordinates and output order as
        /// <see cref="FindOrfs(DnaSequence, OrfFindMode, GeneticCode?, int, int, int, bool, bool, bool)"/>.
        /// Input: case-insensitive IUPAC DNA or RNA (A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N; U is
        /// read as T, as getorf does); nucleotide output (modes 2–6) is the upper-case DNA spelling, and
        /// the reverse strand is the canonical IUPAC reverse complement
        /// (<see cref="DnaSequence.GetReverseComplementString(string)"/>; R↔Y, K↔M, B↔V, D↔H).
        /// </para>
        /// <para>
        /// Ambiguous codons. START: only an unambiguous start codon opens an ORF — EMBOSS
        /// <c>ajTrnCodonstrTypeC</c> reads the Starts table, which EMBOSS fills for the 64 A/C/G/T codons
        /// only (an ambiguous codon such as YTG in table 1 is never a START). STOP: a codon ends an ORF
        /// only when every concrete codon it stands for is a stop (<see cref="GeneticCode.Translate(string)"/>
        /// returns '*': TAR, TRA in table 1); a possible stop such as TAN or NNN reads through as 'X'.
        /// Residues (modes 0, 1): <see cref="GeneticCode.Translate(string)"/> (Biopython ambiguous tables).
        /// EMBOSS instead resolves an ambiguous codon from a precomputed table (<c>trnComplete</c>: third-base
        /// wobble, plus one "most ambiguous codon" per amino acid built by OR-ing the codons' bases); that
        /// table is right for the wobble codons and translates every other ambiguous codon as 'X', but its
        /// OR-ed entries can be wrong — table 1: TRR → '*' although TGG is Trp, TRA (TAA/TGA, all stops) → X,
        /// YTN → L although TTT/TTC are Phe. Those codons are the documented differences from the getorf
        /// binary (see the TRANS-SIXFRAME-001 Evidence); every other ambiguous codon gives getorf's ORFs,
        /// and the residue differs from getorf only where getorf prints 'X' and the genetic code still
        /// determines the residue (e.g. RAY → B in table 1).
        /// </para>
        /// </remarks>
        /// <param name="sequence">IUPAC nucleotide sequence (DNA or RNA, case-insensitive).</param>
        /// <param name="mode">The getorf <c>-find</c> mode.</param>
        /// <param name="geneticCode">The genetic code (default: Standard, NCBI table 1).</param>
        /// <param name="minSize">getorf <c>-minsize</c>: minimum ORF size in nucleotides (default 30).</param>
        /// <param name="maxSize">getorf <c>-maxsize</c>: maximum ORF size in nucleotides (default 1 000 000).</param>
        /// <param name="flanking">getorf <c>-flanking</c>: bases reported on each side of the codon in modes 4–6 (default 100).</param>
        /// <param name="searchBothStrands">getorf <c>-reverse</c>: also search the reverse complement (default true).</param>
        /// <param name="methionine">getorf <c>-methionine</c>: report the initiator residue as Met in mode 1 (default true).</param>
        /// <param name="circular">getorf <c>-circular</c>: treat the sequence as circular (default false).</param>
        /// <returns>The reported ORFs / flank windows in getorf output order.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="mode"/> is undefined, or a size / flanking value is negative.</exception>
        /// <exception cref="ArgumentException"><paramref name="sequence"/> contains a non-IUPAC symbol, or the
        /// genetic code has dual-coding stop codons (tables 27, 28, 31).</exception>
        public static IEnumerable<OrfFeature> FindOrfs(string sequence, OrfFindMode mode, GeneticCode? geneticCode = null,
            int minSize = GetOrfDefaultMinSize, int maxSize = GetOrfDefaultMaxSize, int flanking = GetOrfDefaultFlanking,
            bool searchBothStrands = true, bool methionine = true, bool circular = false)
        {
            string normalized = NormalizeIupacNucleotides(sequence);
            var settings = CreateGetOrfSettings(mode, geneticCode, minSize, maxSize, flanking, methionine);
            return FindOrfsGetOrfCore(normalized, settings, searchBothStrands, circular);
        }

        private static GetOrfSettings CreateGetOrfSettings(OrfFindMode mode, GeneticCode? geneticCode,
            int minSize, int maxSize, int flanking, bool methionine)
        {
            if (!Enum.IsDefined(mode))
                throw new ArgumentOutOfRangeException(nameof(mode), mode, "Unknown getorf -find mode (0–6).");
            ArgumentOutOfRangeException.ThrowIfNegative(minSize);
            ArgumentOutOfRangeException.ThrowIfNegative(maxSize);
            ArgumentOutOfRangeException.ThrowIfNegative(flanking);
            var code = RequireOrfGeneticCode(geneticCode);
            return new GetOrfSettings(mode, code, minSize, maxSize, flanking, methionine);
        }

        private readonly record struct GetOrfSettings(OrfFindMode Mode, GeneticCode Code,
            int MinSize, int MaxSize, int Flanking, bool Methionine);

        private static IEnumerable<OrfFeature> FindOrfsGetOrfCore(string sequence, GetOrfSettings settings,
            bool searchBothStrands, bool circular)
        {
            int length = sequence.Length;
            if (length == 0)
                yield break;

            // getorf main(): "if the sequence is circular, append it to itself to triple its length so can
            // deal easily with wrapped ORFs, but don't update len"; the reverse pass reverse-complements
            // the tripled sequence (ajSeqReverseForce).
            string scanned = circular ? string.Concat(sequence, sequence, sequence) : sequence;

            foreach (var feature in ScanGetOrf(scanned, length, settings, isReverseComplement: false, circular))
                yield return feature;

            if (searchBothStrands)
            {
                foreach (var feature in ScanGetOrf(DnaSequence.GetReverseComplementString(scanned), length, settings,
                    isReverseComplement: true, circular))
                    yield return feature;
            }
        }

        // getorf_FindORFs: one pass over every codon position of the (possibly tripled) sequence;
        // frame = pos % 3; len = length of the original sequence.
        private static IEnumerable<OrfFeature> ScanGetOrf(string sequence, int len, GetOrfSettings s,
            bool isReverseComplement, bool circular)
        {
            var mode = s.Mode;
            // getorf main(): "get the minimum size converted to protein length if storing protein sequences".
            bool storesProtein = mode is OrfFindMode.TranslationStopToStop or OrfFindMode.TranslationStartToStop
                or OrfFindMode.FlankingStart;
            int min = storesProtein ? s.MinSize / CodonLength : s.MinSize;
            int max = storesProtein ? s.MaxSize / CodonLength : s.MaxSize;
            bool startToStop = mode is OrfFindMode.TranslationStartToStop or OrfFindMode.NucleotidesStartToStop
                or OrfFindMode.FlankingStart;

            int n = sequence.Length;
            var open = new bool[ReadingFramesPerStrand];
            var start = new int[ReadingFramesPerStrand];
            var buffer = new StringBuilder[ReadingFramesPerStrand];
            // LASTORF: the frame hit its last STOP / end codon past the end of a circular genome.
            var finished = new bool[ReadingFramesPerStrand];
            // GOTSTOP: the frame contains at least one STOP (circular STOP→STOP modes).
            var gotStop = new bool[ReadingFramesPerStrand];
            for (int f = 0; f < ReadingFramesPerStrand; f++)
            {
                buffer[f] = new StringBuilder();
                // Linear STOP→STOP modes: "assume already in a ORF so we get ORFs at the start of the sequence".
                open[f] = !circular && !startToStop;
                start[f] = f;
            }

            for (int pos = 0; pos + CodonLength <= n; pos++)
            {
                int frame = pos % ReadingFramesPerStrand;
                // "don't want to find extra ORFs when already been round circ"
                if (finished[frame])
                    continue;

                string codon = sequence.Substring(pos, CodonLength);
                char aa = s.Code.Translate(codon);
                bool isStop = aa == '*';
                // getorf: pos >= seqlen-5 is the last complete codon of this frame.
                bool isLastCodon = pos + 2 * CodonLength > n;
                bool pastGenomeEnd = circular && pos >= len;

                if (!startToStop)
                {
                    if (isStop)
                        gotStop[frame] = true;
                    if (isStop || isLastCodon)
                    {
                        if (!isStop)
                            AppendGetOrf(buffer[frame], sequence, pos, aa, storesProtein);
                        if (open[frame] && SizeOk(buffer[frame].Length, min, max)
                            && TryBuildGetOrfFeature(sequence, len, s, isReverseComplement, start[frame],
                                pos + CodonLength - 1, isStop, buffer[frame], out var feature))
                            yield return feature;
                        buffer[frame].Clear();
                        if (pastGenomeEnd)
                        {
                            open[frame] = false;
                            finished[frame] = true;
                            if (finished[0] && finished[1] && finished[2])
                                break;
                        }
                        else
                        {
                            open[frame] = true;
                            start[frame] = pos + CodonLength;
                        }
                    }
                    else if (open[frame])
                    {
                        AppendGetOrf(buffer[frame], sequence, pos, aa, storesProtein);
                    }
                }
                else if (!open[frame] && IsOrfStartCodon(s.Code, codon))
                {
                    // A START past the end of a circular genome opens nothing (getorf: if(pos < len)).
                    if (pos < len)
                    {
                        buffer[frame].Clear();
                        open[frame] = true;
                        start[frame] = pos;
                        AppendGetOrf(buffer[frame], sequence, pos, s.Methionine ? InitiatorMethionine : aa, storesProtein);
                    }
                }
                else if (isStop || isLastCodon)
                {
                    if (open[frame])
                    {
                        open[frame] = false;
                        if (!isStop)
                            AppendGetOrf(buffer[frame], sequence, pos, aa, storesProtein);
                        if (SizeOk(buffer[frame].Length, min, max)
                            && TryBuildGetOrfFeature(sequence, len, s, isReverseComplement, start[frame],
                                pos + CodonLength - 1, isStop, buffer[frame], out var feature))
                            yield return feature;
                    }
                    if (pastGenomeEnd)
                    {
                        finished[frame] = true;
                        if (finished[0] && finished[1] && finished[2])
                            break;
                    }
                    buffer[frame].Clear();
                }
                else if (open[frame])
                {
                    AppendGetOrf(buffer[frame], sequence, pos, aa, storesProtein);
                }
            }

            // Documented divergence from getorf: on a linear sequence a START that is the last complete
            // codon of its frame is still open here (getorf drops it); report it as a one-codon open ORF.
            if (startToStop && !circular)
            {
                for (int f = 0; f < ReadingFramesPerStrand; f++)
                {
                    if (open[f] && SizeOk(buffer[f].Length, min, max)
                        && TryBuildGetOrfFeature(sequence, len, s, isReverseComplement, start[f],
                            start[f] + CodonLength - 1, false, buffer[f], out var feature))
                        yield return feature;
                }
            }

            // getorf: "Currently miss reporting a STOP-to-STOP ORF that is the full length of a circular
            // genome when there are no STOP codons in that frame" — such a frame is reported here as one
            // open record over its complete codons from offset f (getorf WriteORF(f, seqlen-1)).
            if (circular && mode is OrfFindMode.TranslationStopToStop or OrfFindMode.NucleotidesStopToStop)
            {
                for (int f = 0; f < ReadingFramesPerStrand; f++)
                {
                    if (gotStop[f])
                        continue;
                    var whole = new StringBuilder();
                    int pos = f;
                    for (; pos + CodonLength <= n; pos += CodonLength)
                        AppendGetOrf(whole, sequence, pos, s.Code.Translate(sequence.Substring(pos, CodonLength)), storesProtein);
                    AppendDanglingBases(whole, sequence, pos, s.Code, storesProtein);
                    if (SizeOk(whole.Length, min, max)
                        && TryBuildGetOrfFeature(sequence, len, s, isReverseComplement, f, n - 1, false, whole, out var feature))
                        yield return feature;
                }
            }
        }

        private static bool SizeOk(int length, int min, int max) => length >= min && length <= max;

        // The circular full-length record covers the frame up to the last base (getorf WriteORF(f, seqlen-1)).
        // Its translation is EMBOSS ajTrnSeqOrig, which translates dangling bases (ajTrnSeqDangleC): two
        // bases b1 b2 as the codon b1 b2 N (an amino acid when the third base does not matter, else X), one
        // base as X. The nucleotide form keeps the dangling bases.
        private static void AppendDanglingBases(StringBuilder buffer, string sequence, int pos, GeneticCode code,
            bool storesProtein)
        {
            int dangle = sequence.Length - pos;
            if (dangle <= 0)
                return;
            if (!storesProtein)
                buffer.Append(sequence, pos, dangle);
            else if (dangle == 2)
                buffer.Append(code.Translate(string.Concat(sequence.AsSpan(pos, 2), "N")));
            else
                buffer.Append('X');
        }

        // getorf_AppORF: protein modes store the residue, nucleotide modes the codon.
        private static void AppendGetOrf(StringBuilder buffer, string sequence, int pos, char aa, bool storesProtein)
        {
            if (storesProtein)
                buffer.Append(aa);
            else
                buffer.Append(sequence, pos, CodonLength);
        }

        // getorf_WriteORF: a region that starts past the end of a circular genome is shifted back by len;
        // ORF modes report the stored sequence; flanking modes the window [codon − flanking,
        // codon + 2 + flanking] of the scanned sequence, skipped when it crosses a sequence end.
        private static bool TryBuildGetOrfFeature(string sequence, int len, GetOrfSettings s, bool isReverseComplement,
            int orfStart, int orfEnd, bool isTerminated, StringBuilder stored, out OrfFeature feature)
        {
            // "Move the reported positions back to start in the range 1..len for readability."
            while (orfStart >= len)
            {
                orfStart -= len;
                orfEnd -= len;
            }

            int frameOffset = orfStart % ReadingFramesPerStrand;
            int frameLabel = FrameLabel(frameOffset, isReverseComplement);
            int codon;
            switch (s.Mode)
            {
                case OrfFindMode.FlankingStart:
                    codon = orfStart;
                    break;
                case OrfFindMode.FlankingInitialStop:
                    codon = orfStart - CodonLength; // < 0 for the leading region: no initial STOP
                    break;
                case OrfFindMode.FlankingEndingStop:
                    // getorf: codonpos = e + 1 (the STOP after the ORF; for an open ORF the codon after
                    // its last codon, whose window passes the end of a linear sequence).
                    codon = isTerminated ? orfEnd - (CodonLength - 1) : orfEnd + 1;
                    break;
                default:
                    feature = new OrfFeature(s.Mode, orfStart, orfEnd, frameLabel, isTerminated, -1, stored.ToString());
                    return true;
            }

            int windowStart = codon - s.Flanking;
            int windowEnd = codon + (CodonLength - 1) + s.Flanking;
            if (windowStart < 0 || windowEnd >= sequence.Length)
            {
                feature = default;
                return false;
            }

            feature = new OrfFeature(s.Mode, windowStart, windowEnd, frameLabel, isTerminated, codon,
                sequence.Substring(windowStart, windowEnd - windowStart + 1));
            return true;
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
            string forward = NormalizeIupacNucleotides(sequence);
            return TranslateSixFramesCore(forward, DnaSequence.GetReverseComplementString(forward),
                geneticCode ?? GeneticCode.Standard);
        }

        /// <summary>
        /// Translates all six reading frames of a DNA sequence under the chosen frame convention.
        /// </summary>
        /// <remarks>
        /// <see cref="SixFrameConvention.Biopython"/> is identical to
        /// <see cref="TranslateSixFrames(DnaSequence, GeneticCode?)"/>. The EMBOSS conventions follow
        /// EMBOSS 6.6.0 <c>transeq</c> (<c>ajTrnSeqOrig</c>): see <see cref="SixFrameConvention"/>.
        /// Frame keys are always −3..−1, +1..+3.
        /// </remarks>
        /// <param name="dna">The DNA sequence to translate.</param>
        /// <param name="geneticCode">The genetic code to use (null = Standard).</param>
        /// <param name="convention">Reverse-frame phasing and trailing-base convention.</param>
        /// <exception cref="ArgumentNullException"><paramref name="dna"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="convention"/> is not a defined value.</exception>
        public static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFrames(DnaSequence dna,
            GeneticCode? geneticCode, SixFrameConvention convention)
        {
            ArgumentNullException.ThrowIfNull(dna);
            return TranslateSixFramesCore(dna.Sequence, dna.ReverseComplement().Sequence,
                geneticCode ?? GeneticCode.Standard, convention);
        }

        /// <summary>
        /// Translates all six reading frames of an IUPAC nucleotide string under the chosen frame
        /// convention (input handling as <see cref="TranslateSixFrames(string, GeneticCode?)"/>).
        /// </summary>
        /// <remarks>
        /// <see cref="SixFrameConvention.Biopython"/> is identical to
        /// <see cref="TranslateSixFrames(string, GeneticCode?)"/>. The EMBOSS conventions follow EMBOSS
        /// 6.6.0 <c>transeq</c>; see <see cref="SixFrameConvention"/>. Ambiguous codons keep this
        /// library's IUPAC resolution (<see cref="GeneticCode.Translate(string)"/>, Biopython tables);
        /// transeq's own ambiguity table differs for a few dozen IUPAC codons per table (it never
        /// emits B/Z/J and resolves some S/W, R/Y codons differently), so only those residues can
        /// differ from the transeq binary.
        /// </remarks>
        /// <param name="sequence">IUPAC nucleotide sequence (DNA or RNA, case-insensitive).</param>
        /// <param name="geneticCode">The genetic code to use (null = Standard).</param>
        /// <param name="convention">Reverse-frame phasing and trailing-base convention.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="sequence"/> contains a non-IUPAC symbol.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="convention"/> is not a defined value.</exception>
        public static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFrames(string sequence,
            GeneticCode? geneticCode, SixFrameConvention convention)
        {
            string forward = NormalizeIupacNucleotides(sequence);
            return TranslateSixFramesCore(forward, DnaSequence.GetReverseComplementString(forward),
                geneticCode ?? GeneticCode.Standard, convention);
        }

        /// <summary>
        /// Translates all six reading frames of a DNA sequence as EMBOSS 6.6.0 <c>transeq</c> does with its
        /// post-processing qualifiers <c>-regions</c>, <c>-methionine</c>, <c>-trim</c> and <c>-clean</c>.
        /// See <see cref="TranslateSixFrames(string, GeneticCode?, SixFrameConvention, TranseqOptions)"/>.
        /// </summary>
        /// <exception cref="ArgumentNullException"><paramref name="dna"/> or <paramref name="options"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="convention"/> is <see cref="SixFrameConvention.Biopython"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A region is invalid, or <paramref name="convention"/> is undefined.</exception>
        public static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFrames(DnaSequence dna,
            GeneticCode? geneticCode, SixFrameConvention convention, TranseqOptions options)
        {
            ArgumentNullException.ThrowIfNull(dna);
            return TranslateSixFrames(dna.Sequence, geneticCode, convention, options);
        }

        /// <summary>
        /// Translates all six reading frames of an IUPAC nucleotide string as EMBOSS 6.6.0 <c>transeq</c>
        /// does with its post-processing qualifiers <c>-regions</c>, <c>-methionine</c>, <c>-trim</c> and
        /// <c>-clean</c> (emboss/transeq.c, acd/transeq.acd).
        /// </summary>
        /// <remarks>
        /// Steps, in transeq's order:
        /// <list type="number">
        /// <item><see cref="TranseqOptions.Regions"/> (<c>-regions</c>): the 1-based inclusive regions are
        /// cut out and concatenated in the given order (<c>ajRangeSeqExtract</c>); the concatenation is
        /// then translated as one sequence in all six frames (its reverse complement for the reverse
        /// frames). An end past the sequence end is clamped to it, as transeq does.</item>
        /// <item><see cref="TranseqOptions.Methionine"/> (<c>-methionine</c>, default on): only when
        /// regions are given (transeq: <c>methionine &amp;&amp; !wholeseq</c>), the first codon of each frame
        /// is read as M when the genetic code lists it as a START codon (<c>ajTrnSeqInit</c>, Starts
        /// table: unambiguous codons only). transeq's phase-locked reverse frames
        /// (<c>ajTrnSeqRevinitC</c>, <c>iend &gt; 2</c>) need at least two whole codons for this;
        /// forward and <c>-alternative</c> reverse frames need one.</item>
        /// <item><see cref="TranseqOptions.Trim"/> (<c>-trim</c>): trailing 'X' and '*' are removed
        /// (<c>transeq_Trim</c>).</item>
        /// <item><see cref="TranseqOptions.Clean"/> (<c>-clean</c>, after trimming): every '*' becomes 'X'
        /// (<c>transeq_Clean</c>).</item>
        /// </list>
        /// With <see cref="TranseqOptions.Default"/> the result equals
        /// <see cref="TranslateSixFrames(string, GeneticCode?, SixFrameConvention)"/>. Only the EMBOSS
        /// conventions are accepted. Codon resolution of ambiguous bases is this library's (see the
        /// convention overload).
        /// </remarks>
        /// <param name="sequence">IUPAC nucleotide sequence (DNA or RNA, case-insensitive).</param>
        /// <param name="geneticCode">The genetic code to use (null = Standard).</param>
        /// <param name="convention"><see cref="SixFrameConvention.EmbossTranseq"/> or
        /// <see cref="SixFrameConvention.EmbossTranseqAlternative"/>.</param>
        /// <param name="options">transeq post-processing qualifiers.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="options"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="sequence"/> contains a non-IUPAC symbol, or
        /// <paramref name="convention"/> is <see cref="SixFrameConvention.Biopython"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A region has Start &lt; 1, End &lt; Start or Start past the
        /// sequence end, or <paramref name="convention"/> is undefined.</exception>
        public static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFrames(string sequence,
            GeneticCode? geneticCode, SixFrameConvention convention, TranseqOptions options)
        {
            ArgumentNullException.ThrowIfNull(options);
            if (convention == SixFrameConvention.Biopython)
                throw new ArgumentException("transeq options require an EMBOSS six-frame convention.",
                    nameof(convention));

            string forward = NormalizeIupacNucleotides(sequence);
            bool wholeSequence = options.Regions is null || options.Regions.Count == 0;
            if (!wholeSequence)
                forward = ExtractTranseqRegions(forward, options.Regions!);

            var code = geneticCode ?? GeneticCode.Standard;
            var frames = TranslateSixFramesCore(forward, DnaSequence.GetReverseComplementString(forward), code,
                convention);

            var result = new Dictionary<int, ProteinSequence>(frames.Count);
            foreach (var (frame, protein) in frames)
            {
                var residues = protein.Sequence.ToCharArray();
                int length = residues.Length;
                if (options.Methionine && !wholeSequence && length > 0
                    && FirstTranseqCodonIsStart(forward, frame, convention, code))
                    residues[0] = 'M';
                if (options.Trim)
                    while (length > 0 && (residues[length - 1] == 'X' || residues[length - 1] == '*'))
                        length--;
                if (options.Clean)
                    for (int i = 0; i < length; i++)
                        if (residues[i] == '*') residues[i] = 'X';
                result[frame] = new ProteinSequence(new string(residues, 0, length));
            }

            return result;
        }

        // transeq -regions (ajRangeSeqExtract): concatenation of seq[start-1 .. end-1] in the given order;
        // ajStrAppendSubS clamps an end past the sequence end.
        private static string ExtractTranseqRegions(string forward, IReadOnlyList<(int Start, int End)> regions)
        {
            var sb = new StringBuilder();
            foreach (var (start, end) in regions)
            {
                if (start < 1 || end < start || start > forward.Length)
                    throw new ArgumentOutOfRangeException(nameof(regions), $"{start}-{end}",
                        $"Region must satisfy 1 <= Start <= End and Start <= sequence length ({forward.Length}).");
                sb.Append(forward, start - 1, Math.Min(end, forward.Length) - start + 1);
            }
            return sb.ToString();
        }

        // ajTrnSeqInit: the first codon read in the frame becomes M when the Starts table marks it.
        //  forward k (ajTrnSeqInitC): forward[k-1 .. k+1], needs one whole codon (lenmod3 > 2);
        //  -alternative −k (ajTrnSeqAltRevinitC): revcomp[k-1 .. k+1], needs one whole codon (iend > 1);
        //  phase-locked −k (ajTrnSeqRevinitC): the last whole frame-k codon, reverse-complemented, which is
        //  revcomp[(len−k+1) mod 3 ..]; needs two whole codons (iend = 3·⌊(len−k+1)/3⌋ − 1 > 2).
        private static bool FirstTranseqCodonIsStart(string forward, int frame, SixFrameConvention convention,
            GeneticCode code)
        {
            int offset = Math.Abs(frame) - 1;
            int available = forward.Length - offset;
            string strand;
            int start;
            if (frame > 0)
            {
                if (available < CodonLength) return false;
                strand = forward;
                start = offset;
            }
            else if (convention == SixFrameConvention.EmbossTranseqAlternative)
            {
                if (available < CodonLength) return false;
                strand = DnaSequence.GetReverseComplementString(forward);
                start = offset;
            }
            else
            {
                if (available < 2 * CodonLength) return false;
                strand = DnaSequence.GetReverseComplementString(forward);
                start = available % CodonLength;
            }

            return IsOrfStartCodon(code, strand.Substring(start, CodonLength));
        }

        private static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFramesCore(string forward,
            string reverseComplement, GeneticCode code, SixFrameConvention convention)
        {
            switch (convention)
            {
                case SixFrameConvention.Biopython:
                    return TranslateSixFramesCore(forward, reverseComplement, code);
                case SixFrameConvention.EmbossTranseq:
                case SixFrameConvention.EmbossTranseqAlternative:
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(convention), convention,
                        "Unknown six-frame convention.");
            }

            var result = new Dictionary<int, ProteinSequence>();
            int len = forward.Length;
            for (int offset = 0; offset < ReadingFramesPerStrand; offset++)
            {
                // transeq frame k = ajTrnSeqOrig(seq, k): codons from base k, then the dangling bases.
                result[FrameLabel(offset, isReverseComplement: false)] = TranslateEmbossStrand(forward, offset, code);

                if (convention == SixFrameConvention.EmbossTranseqAlternative)
                {
                    // -alternative Y: ajTrnSeqOrig(seq, -4-offset) — frame -k starts at the k-th base of the
                    // reverse complement (= Biopython), dangling bases translated as for forward frames.
                    result[FrameLabel(offset, isReverseComplement: true)] = TranslateEmbossStrand(reverseComplement, offset, code);
                }
                else
                {
                    result[FrameLabel(offset, isReverseComplement: true)] = TranslateEmbossPhaseLocked(reverseComplement, len, offset, code);
                }
            }

            return result;
        }

        // ajTrnSeqFrameC (frames 1..3, -4..-6) + ajTrnSeqDangleS: full codons of strand[offset..], then
        // (strand.Length − offset) mod 3 dangling bases — two bases b1 b2 read as the codon b1 b2 N
        // (an amino acid when the third base does not matter, else X), one base as X.
        private static ProteinSequence TranslateEmbossStrand(string strand, int offset, GeneticCode code)
        {
            if (strand.Length <= offset)
                return new ProteinSequence(string.Empty);

            var sb = new StringBuilder(TranslateSequence(strand, code, offset, false).Sequence);
            int remaining = (strand.Length - offset) % CodonLength;
            AppendDanglingBases(sb, strand, strand.Length - remaining, code, storesProtein: true);
            return new ProteinSequence(sb.ToString());
        }

        // transeq default ("frame -1 uses frame 1 codons", ajTrnSeqRevC on seq[offset..] + ajTrnSeqDangleS):
        // the reverse complement of the forward frame-(offset+1) codons, then the offset leading forward
        // bases as the dangling end (frame −2: X; frame −3: codon comp(b2) comp(b1) N). The forward bases
        // after the last full frame codon are not translated. On the reverse complement the codons start
        // at (len − offset) mod 3 and the last offset bases dangle.
        private static ProteinSequence TranslateEmbossPhaseLocked(string reverseComplement, int len, int offset,
            GeneticCode code)
        {
            if (len == 0)
                return new ProteinSequence(string.Empty);
            if (len < offset)
            {
                // Only len = 1, frame −3: transeq still emits its dangling-pair residue, and with one
                // base the pair is never decidable (ajTrnSeqDangleS reads past the base) → X.
                return new ProteinSequence("X");
            }

            int start = (len - offset) % CodonLength;
            var sb = new StringBuilder();
            if (len - start >= CodonLength)
                sb.Append(TranslateSequence(reverseComplement, code, start, false).Sequence);
            AppendDanglingBases(sb, reverseComplement, len - offset, code, storesProtein: true);
            return new ProteinSequence(sb.ToString());
        }

        private static IReadOnlyDictionary<int, ProteinSequence> TranslateSixFramesCore(string forward,
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

        private static ProteinSequence TranslateSequence(string sequence, GeneticCode geneticCode,
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
            // GeneticCode.Translate / IsStartCodon normalise T/U themselves.
            // ORF = region from a START codon to a STOP codon
            // (EMBOSS getorf -find 1: "a region that begins with a START codon
            // and ends with a STOP codon"). Scanned in all three frames.
            for (int frame = 0; frame < ReadingFramesPerStrand; frame++)
            {
                int? currentOrfStart = null;
                var currentProtein = new StringBuilder();

                for (int i = frame; i + CodonLength <= sequence.Length; i += CodonLength)
                {
                    string codon = sequence.Substring(i, CodonLength);
                    char aa = geneticCode.Translate(codon);

                    if (currentOrfStart == null)
                    {
                        // Looking for start codon
                        if (IsOrfStartCodon(geneticCode, codon))
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
                                    FrameLabel(frame, isReverseComplement),
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
                        FrameLabel(frame, isReverseComplement),
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

        /// <summary>
        /// Gets whether the ORF ends with a STOP codon (its stop-inclusive range is one codon longer
        /// than the protein); false for an open ORF that runs to the end of the strand.
        /// </summary>
        public bool IsTerminated => NucleotideLength > Protein.Length * 3;

        /// <summary>
        /// Gets EMBOSS getorf's 1-based <c>[s - e]</c> header range for this ORF, in forward-strand
        /// coordinates (getorf_WriteORF) — the same mapping as <see cref="OrfFeature.ToGetOrfRange(int)"/>
        /// for a <see cref="OrfFindMode.TranslationStartToStop"/> record: the terminating STOP is excluded
        /// and a reverse-strand ORF is mapped with <c>s = L − s + 1</c>, <c>e = L − e + 1</c> (so s &gt; e).
        /// </summary>
        /// <param name="sequenceLength">Length L of the searched sequence.</param>
        /// <returns>The (start, end) pair getorf prints.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="sequenceLength"/> ≤ 0.</exception>
        public (int Start, int End) ToGetOrfRange(int sequenceLength) =>
            OrfFeature.MapToGetOrfRange(StartPosition, IsTerminated ? EndPosition - 3 : EndPosition,
                Frame < 0, sequenceLength);
    }

    /// <summary>
    /// Reverse-frame phasing and trailing-base convention of
    /// <see cref="Translator.TranslateSixFrames(string, GeneticCode?, SixFrameConvention)"/>.
    /// </summary>
    public enum SixFrameConvention
    {
        /// <summary>
        /// Biopython <c>SeqUtils.six_frame_translations</c> (default): frame +k translates
        /// <c>seq[k-1:]</c>, frame −k translates <c>reverse_complement(seq)[k-1:]</c>; trailing bases
        /// that do not form a codon are dropped.
        /// </summary>
        Biopython = 0,

        /// <summary>
        /// EMBOSS <c>transeq</c> default (<c>-alternative N</c>, "phase-locked"): frame −k is the reverse
        /// complement of the codons of frame +k (the Staden convention: "frame -1 is the
        /// reverse-complement of the set of codons used in frame 1"), so forward bases after the last
        /// full frame +k codon are not translated, and the k − 1 leading forward bases form a dangling
        /// end (frame −2: X; frame −3: the codon comp(b2) comp(b1) N). In every frame dangling bases
        /// are translated as transeq does (ajTrnSeqDangleS): two bases b1 b2 as the codon b1 b2 N
        /// (an amino acid when the third base does not matter, else X), one base as X.
        /// </summary>
        EmbossTranseq = 1,

        /// <summary>
        /// EMBOSS <c>transeq -alternative Y</c>: frame −k starts at the k-th base of the reverse
        /// complement ("frame -1 as using the set of codons starting with the last codon of the
        /// sequence") — the Biopython phasing — with dangling bases translated as in
        /// <see cref="EmbossTranseq"/>.
        /// </summary>
        EmbossTranseqAlternative = 2,
    }

    /// <summary>
    /// EMBOSS 6.6.0 <c>transeq</c> post-processing qualifiers (emboss/acd/transeq.acd) for
    /// <see cref="Translator.TranslateSixFrames(string, GeneticCode?, SixFrameConvention, TranseqOptions)"/>.
    /// Defaults equal transeq's: no regions, no trimming, no cleaning, <c>-methionine Y</c>.
    /// </summary>
    public sealed record TranseqOptions
    {
        /// <summary>transeq defaults (whole sequence; the output equals the plain convention overload).</summary>
        public static TranseqOptions Default { get; } = new();

        /// <summary>
        /// <c>-trim</c>: "This removes all 'X' and '*' characters from the right end of the translation"
        /// (default false).
        /// </summary>
        public bool Trim { get; init; }

        /// <summary>
        /// <c>-clean</c>: "This changes all STOP codon positions from the '*' character to 'X'" — applied
        /// after <see cref="Trim"/> (default false).
        /// </summary>
        public bool Clean { get; init; }

        /// <summary>
        /// <c>-regions</c>: 1-based inclusive (Start, End) regions, concatenated in the given order and
        /// translated as one sequence (null or empty = whole sequence). An End past the sequence end is
        /// clamped to it, as transeq does.
        /// </summary>
        public IReadOnlyList<(int Start, int End)>? Regions { get; init; }

        /// <summary>
        /// <c>-methionine</c> (transeq default true): "Change initial START codons to Methionine". transeq
        /// applies it only when <see cref="Regions"/> are given (<c>methionine &amp;&amp; !wholeseq</c>).
        /// </summary>
        public bool Methionine { get; init; } = true;
    }

    /// <summary>
    /// EMBOSS getorf <c>-find</c> modes (emboss/acd/getorf.acd; numeric values equal getorf's codes).
    /// </summary>
    public enum OrfFindMode
    {
        /// <summary>0: translation of regions between STOP codons.</summary>
        TranslationStopToStop = 0,

        /// <summary>1: translation of regions between START and STOP codons.</summary>
        TranslationStartToStop = 1,

        /// <summary>2: nucleic sequences between STOP codons.</summary>
        NucleotidesStopToStop = 2,

        /// <summary>3: nucleic sequences between START and STOP codons.</summary>
        NucleotidesStartToStop = 3,

        /// <summary>4: nucleotides flanking START codons (of START→STOP ORFs).</summary>
        FlankingStart = 4,

        /// <summary>5: nucleotides flanking initial STOP codons (of STOP→STOP ORFs).</summary>
        FlankingInitialStop = 5,

        /// <summary>6: nucleotides flanking ending STOP codons (of STOP→STOP ORFs).</summary>
        FlankingEndingStop = 6,
    }

    /// <summary>
    /// One record reported by
    /// <see cref="Translator.FindOrfs(DnaSequence, OrfFindMode, GeneticCode?, int, int, int, bool, bool, bool)"/>
    /// (one getorf output entry).
    /// </summary>
    /// <param name="Mode">The getorf <c>-find</c> mode that produced the record.</param>
    /// <param name="StartPosition">0-based first base of the region, in the scanned strand's coordinates
    /// (reverse-complement coordinates when <see cref="Frame"/> is negative).</param>
    /// <param name="EndPosition">0-based inclusive last base of the region. Modes 0–3: the last base of
    /// the terminating STOP codon when <see cref="IsTerminated"/> (INSDC CDS convention), otherwise of
    /// the last complete codon. Modes 4–6: the last base of the flank window.</param>
    /// <param name="Frame">Frame label: +1..+3 forward (offset 0..2), −1..−3 reverse complement (offset 0..2).</param>
    /// <param name="IsTerminated">Whether the underlying ORF ends with a STOP codon (false: it runs to the strand end).</param>
    /// <param name="CodonPosition">Modes 4–6: 0-based first base of the flanked START / STOP codon; −1 otherwise.</param>
    /// <param name="Sequence">The getorf output sequence: protein (modes 0, 1; STOP excluded), the ORF's
    /// nucleotides (modes 2, 3; STOP excluded) or the flank window's nucleotides (modes 4–6).</param>
    public readonly record struct OrfFeature(
        OrfFindMode Mode,
        int StartPosition,
        int EndPosition,
        int Frame,
        bool IsTerminated,
        int CodonPosition,
        string Sequence)
    {
        /// <summary>Gets whether the record lies on the reverse-complement strand.</summary>
        public bool IsReverseStrand => Frame < 0;

        /// <summary>
        /// Gets getorf's 1-based <c>[s - e]</c> header range for this record, in forward-strand
        /// coordinates (getorf_WriteORF). ORF modes exclude the terminating STOP; a reverse-strand
        /// record is mapped with <c>s = L − s + 1</c>, <c>e = L − e + 1</c> (so s &gt; e), including
        /// getorf's shift of out-of-range values for empty ORFs.
        /// </summary>
        /// <param name="sequenceLength">Length L of the searched sequence.</param>
        /// <returns>The (start, end) pair getorf prints.</returns>
        public (int Start, int End) ToGetOrfRange(int sequenceLength)
        {
            bool orfMode = Mode is OrfFindMode.TranslationStopToStop or OrfFindMode.TranslationStartToStop
                or OrfFindMode.NucleotidesStopToStop or OrfFindMode.NucleotidesStartToStop;
            return MapToGetOrfRange(StartPosition, orfMode && IsTerminated ? EndPosition - 3 : EndPosition,
                IsReverseStrand, sequenceLength);
        }

        // getorf_WriteORF header mapping, shared with OrfResult.ToGetOrfRange: 0-based inclusive
        // [startPosition, lastBase] (STOP already excluded) in the scanned strand's coordinates →
        // getorf's 1-based forward-strand (s, e), including its range shifts for empty ORFs.
        internal static (int Start, int End) MapToGetOrfRange(int startPosition, int lastBase,
            bool reverseStrand, int sequenceLength)
        {
            ArgumentOutOfRangeException.ThrowIfNegativeOrZero(sequenceLength);
            int s = startPosition + 1;
            int e = lastBase + 1;
            if (!reverseStrand)
                return (s, e);

            s = sequenceLength - s + 1;
            e = sequenceLength - e + 1;
            while (e > sequenceLength)
            {
                s -= sequenceLength;
                e -= sequenceLength;
            }
            while (e < 0 || s < 0)
            {
                s += sequenceLength;
                e += sequenceLength;
            }
            return (s, e);
        }
    }
}
