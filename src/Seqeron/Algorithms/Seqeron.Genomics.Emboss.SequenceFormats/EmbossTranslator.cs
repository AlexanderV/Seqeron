// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

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

/// <summary>The translation conventions of the EMBOSS 6.6.0 <c>transeq</c> program (<c>ajTrnSeqOrig</c>, dangling bases, <c>-regions</c> / <c>-trim</c> / <c>-clean</c>), moved out of the MIT <c>Translator</c> (B09 WP80).</summary>
public static class EmbossTranslator
{
    private const int CodonLength = 3;
    private const int ReadingFramesPerStrand = 3;

    // Signed 1-based frame label of a 0-based strand offset: +1..+3 on the forward strand, -1..-3 on the reverse complement.
    private static int FrameLabel(int offset, bool isReverseComplement) =>
        isReverseComplement ? -(offset + 1) : offset + 1;

    private static bool IsOrfStartCodon(GeneticCode code, string codon) => Translator.IsOrfStartCodon(code, codon);

    private static ProteinSequence TranslateSequence(string sequence, GeneticCode geneticCode, int frame, bool toFirstStop, bool allowNucleotideX = false) =>
        Translator.TranslateSequence(sequence, geneticCode, frame, toFirstStop, allowNucleotideX);

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
        string forward = SequenceExtensions.NormalizeIupacDna(sequence);
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

        string forward = SequenceExtensions.NormalizeIupacDna(sequence);
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
                return Translator.TranslateSixFramesCore(forward, reverseComplement, code);
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
    // The circular full-length record covers the frame up to the last base (getorf WriteORF(f, seqlen-1)).
    // Its translation is EMBOSS ajTrnSeqOrig, which translates dangling bases (ajTrnSeqDangleC): two
    // bases b1 b2 as the codon b1 b2 N (an amino acid when the third base does not matter, else X), one
    // base as X. The nucleotide form keeps the dangling bases.
    public static void AppendDanglingBases(StringBuilder buffer, string sequence, int pos, GeneticCode code,
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
}
