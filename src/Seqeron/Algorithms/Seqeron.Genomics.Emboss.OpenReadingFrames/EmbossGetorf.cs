// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.OpenReadingFrames, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Emboss.SequenceFormats;

namespace Seqeron.Genomics.Emboss.OpenReadingFrames;

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
/// <see cref="EmbossGetorf.FindOrfs(DnaSequence, OrfFindMode, GeneticCode?, int, int, int, bool, bool, bool)"/>
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

/// <summary>The model of the EMBOSS 6.6.0 <c>getorf</c> program (<c>getorf_FindORFs</c>, modes 0-6), moved out of the MIT <c>Translator</c> (B09 WP80).</summary>
public static class EmbossGetorf
{
    private const int CodonLength = 3;
    private const int ReadingFramesPerStrand = 3;

    // Signed 1-based frame label of a 0-based strand offset: +1..+3 on the forward strand, -1..-3 on the reverse complement.
    private static int FrameLabel(int offset, bool isReverseComplement) =>
        isReverseComplement ? -(offset + 1) : offset + 1;

    private const char InitiatorMethionine = 'M';

    private static GeneticCode RequireOrfGeneticCode(GeneticCode? code) => Translator.RequireOrfGeneticCode(code);

    private static void AppendDanglingBases(StringBuilder buffer, string sequence, int pos, GeneticCode code, bool storesProtein) => EmbossTranslator.AppendDanglingBases(buffer, sequence, pos, code, storesProtein);

    private static bool IsOrfStartCodon(GeneticCode code, string codon) => Translator.IsOrfStartCodon(code, codon);

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
        string normalized = SequenceExtensions.NormalizeIupacDna(sequence);
        var settings = CreateGetOrfSettings(mode, geneticCode, minSize, maxSize, flanking, methionine);
        return FindOrfsGetOrfCore(normalized, settings, searchBothStrands, circular);
    }
    /// <summary>
    /// The model of the EMBOSS 6.6.0 <c>getorf</c> PROGRAM (B09 WP61, F280): <c>getorf_FindORFs</c> of <c>getorf.c</c> statement by statement on the residues as the program holds them (any case, U kept, the
    /// reverse strand by <c>ajSeqstrReverse</c>), with the translation tables of the EMBOSS genetic code files <c>EGC.n</c> (<c>ajtrn.c</c>: table 0 = standard with AUG start only; the ambiguity codes by
    /// <c>trnComplete</c>). Differences from <see cref="FindOrfs(string, OrfFindMode, GeneticCode?, int, int, int, bool, bool, bool)"/> (which classifies codons by the NCBI tables and corrects four quirks of the
    /// program): a START codon at the last complete codon of its frame opens an ORF that is never written, the full-length record of a circular frame without STOP is always the translation and its frame-3 size
    /// test is <c>&gt;= maxsize</c>, ambiguity codes follow <c>trnComplete</c>.
    /// </summary>
    /// <param name="sequence">The residues of the sequence (nucleotide letters of any case).</param>
    /// <param name="mode">The getorf <c>-find</c> mode.</param>
    /// <param name="emboss">The number of the EMBOSS genetic code file (0..6, 9..16, 21..23).</param>
    /// <param name="minSize">getorf <c>-minsize</c> (nucleotides).</param>
    /// <param name="maxSize">getorf <c>-maxsize</c> (nucleotides).</param>
    /// <param name="flanking">getorf <c>-flanking</c>.</param>
    /// <param name="searchBothStrands">getorf <c>-reverse</c>.</param>
    /// <param name="methionine">getorf <c>-methionine</c>.</param>
    /// <param name="circular">getorf <c>-circular</c>.</param>
    /// <returns>The records in getorf output order; the k-th is <c>name_k</c>.</returns>
    /// <exception cref="ArgumentException">There is no genetic code file <c>EGC.n</c>.</exception>
    public static IEnumerable<OrfFeature> FindOrfsGetorfProgram(string sequence, OrfFindMode mode, int emboss = 0,
        int minSize = GetOrfDefaultMinSize, int maxSize = GetOrfDefaultMaxSize, int flanking = GetOrfDefaultFlanking,
        bool searchBothStrands = true, bool methionine = true, bool circular = false)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (!Enum.IsDefined(mode))
            throw new ArgumentOutOfRangeException(nameof(mode), mode, "Unknown getorf -find mode (0–6).");
        EmbossTrn trn;
        try
        {
            trn = EmbossTrn.New(emboss);
        }
        catch (EmbossFatalException e)
        {
            throw new ArgumentException(e.Message.TrimEnd('\n'), nameof(emboss), e);
        }

        var settings = new GetOrfSettings(mode, GeneticCode.Standard, minSize, maxSize, flanking, methionine, trn);
        return FindOrfsGetOrfCore(sequence, settings, searchBothStrands, circular);
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
    // Trn != null: the model of the getorf PROGRAM (B09 WP61, F280): the codons are classified by the EMBOSS translation table of that genetic code file (ajtrn.c: wobble and
    // "most ambiguous codon" resolution of ambiguity codes, start codons of the Starts line), the residues are scanned as they are (case, U) and the control flow is that of getorf_FindORFs
    // statement by statement (a START at the last codon of its frame opens an ORF that is never written; the full-length record of a circular frame that has no STOP is the translation, in every
    // mode, and its frame-3 size test is the one of getorf.c).
    private readonly record struct GetOrfSettings(OrfFindMode Mode, GeneticCode Code,
        int MinSize, int MaxSize, int Flanking, bool Methionine, EmbossTrn? Trn = null);
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
            string reverse = settings.Trn is not null ? EmbossNuc.ReverseComplement(scanned) : DnaSequence.GetReverseComplementString(scanned);
            foreach (var feature in ScanGetOrf(reverse, length, settings, isReverseComplement: true, circular))
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

            char aa;
            bool isStop, isStartCodon;
            if (s.Trn is { } trn)
            {
                int type = trn.CodonType(sequence, pos, out aa);
                isStop = type == -1;
                isStartCodon = type == 1;
            }
            else
            {
                string codon = sequence.Substring(pos, CodonLength);
                aa = s.Code.Translate(codon);
                isStop = aa == '*';
                isStartCodon = startToStop && !open[frame] && IsOrfStartCodon(s.Code, codon);
            }

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
            else if (!open[frame] && isStartCodon)
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
        if (startToStop && !circular && s.Trn is null)
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
                if (s.Trn is { } emboss)
                {
                    // getorf.c: ajTrnSeqOrig of the frame, whatever the mode; the size test of frame 3 compares with the maximum the wrong way round (>= maxsize)
                    string pep = emboss.SeqOrig(sequence, f + 1, false);
                    if (pep.Length >= min && (f == 2 ? pep.Length >= max : pep.Length <= max)
                        && TryBuildGetOrfFeature(sequence, len, s, isReverseComplement, f, n - 1, false, new StringBuilder(pep), out var embossFeature))
                        yield return embossFeature;
                    continue;
                }

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
}

/// <summary>The getorf header range of an <see cref="OrfResult"/> (<c>getorf_WriteORF</c>).</summary>
public static class OrfResultGetorfExtensions
{
    /// <summary>
    /// Gets EMBOSS getorf's 1-based <c>[s - e]</c> header range for this ORF, in forward-strand coordinates (getorf_WriteORF) - the same mapping as <see cref="OrfFeature.ToGetOrfRange(int)"/>
    /// for a <see cref="OrfFindMode.TranslationStartToStop"/> record: the terminating STOP is excluded and a reverse-strand ORF is mapped with <c>s = L - s + 1</c>, <c>e = L - e + 1</c> (so s &gt; e).
    /// </summary>
    /// <param name="orf">The ORF.</param>
    /// <param name="sequenceLength">Length L of the searched sequence.</param>
    /// <returns>The (start, end) pair getorf prints.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="sequenceLength"/> &lt;= 0.</exception>
    public static (int Start, int End) ToGetOrfRange(this OrfResult orf, int sequenceLength) =>
        OrfFeature.MapToGetOrfRange(orf.StartPosition, orf.IsTerminated ? orf.EndPosition - 3 : orf.EndPosition, orf.Frame < 0, sequenceLength);
}
