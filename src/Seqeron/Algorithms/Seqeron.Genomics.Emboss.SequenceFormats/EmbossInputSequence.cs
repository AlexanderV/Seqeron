// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Alignment;

// A5-7 (B09 F172): the sequence-input qualifiers -sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein, applied once for every EMBOSS port
// (the reader's cleaning stays in Infrastructure/EmbossSequenceReader; this class needs the canonical nucleotide complement of Core).
public static class EmbossInputSequence
{
    /// <summary>
    /// A plain input sequence as an EMBOSS program receives it with the qualifiers of <paramref name="qualifiers"/>: the reader's cleaning (<see cref="EmbossSequenceReader.Clean(string)"/>), the repair the forced type imposes under the ACD type
    /// of the input (<paramref name="type"/>), <c>-slower</c> / <c>-supper</c>, <c>-sreverse</c> and the range of <c>-sbegin</c> / <c>-send</c> (<c>ajSeqTrim</c>).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Sources (EMBOSS 6.6, raw.githubusercontent.com/kimrutherford/EMBOSS, read in full for these functions): <c>ajacd.c</c> <c>acdSetSeq</c> / <c>acdSetSeqall</c> (the qualifiers; <c>-sreverse</c> takes effect only when <c>ajSeqIsNuc</c> of the
    /// sequence is true), <c>ajseqtype.c</c> <c>ajSeqTypeCheckIn</c> (a forced type: <c>-snucleotide</c> / <c>-sprotein</c> set the sequence type, a nucleotide ACD type refuses a protein, a protein type refuses a nucleotide, and a letter outside the nucleotide
    /// codes of a sequence forced to nucleotide becomes <c>N</c>: <c>seqTypeFix</c>), <c>ajseqread.c</c> (<c>-supper</c>, <c>-slower</c> at read time), <c>ajseq.c</c> <c>ajSeqReverseDo</c> (reverse complement, the range
    /// becomes <c>Begin = -End</c>, <c>End = -Begin</c>), <c>ajSeqTrim</c> (<c>ajCvtSposToPosStart</c>: both ends clipped into the sequence, the end not before the start) and <c>ajmath.c</c> <c>ajCvtSposToPosStart</c>.
    /// Measured on the installed 6.6.0 on random sequences and qualifiers (scripts/alignment/generate_emboss_seqopts_oracle.py, 180 rows per program).
    /// </para>
    /// <para>
    /// The range of a reversed sequence refers to the original: <c>-sbegin 3 -send 20 -sreverse</c> uses the reverse complement of original positions 3..20, and the program reports positions of the reversed strand
    /// (<see cref="EmbossSequenceView.Offset"/> counts from the start of the reverse complement).
    /// </para>
    /// </remarks>
    /// <param name="sequence">The sequence as given (gap characters and other non-letters are dropped first).</param>
    /// <param name="qualifiers">The qualifiers; null = none.</param>
    /// <param name="type">The ACD type of the program's input (decides the type check of a forced type).</param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentException">The sequence is empty after cleaning ("Unable to read sequence"), or the forced type contradicts the ACD type ("Sequence is not nucleic" / "Sequence is not a protein").</exception>
    public static EmbossSequenceView Read(string sequence, EmbossSequenceQualifiers? qualifiers = null, EmbossSequenceType type = EmbossSequenceType.Any)
    {
        string text = EmbossSequenceReader.Clean(sequence);
        return ReadCleaned(text, qualifiers ?? EmbossSequenceQualifiers.None, type);
    }

    /// <summary><see cref="Read(string, EmbossSequenceQualifiers?, EmbossSequenceType)"/> of a sequence that is already clean (the reader's cleaning was applied).</summary>
    public static EmbossSequenceView ReadCleaned(string text, EmbossSequenceQualifiers qualifiers, EmbossSequenceType type)
        => ReadCleaned(text, qualifiers, type, EmbossInputKind.Sequence, null, 0);

    /// <summary>
    /// <see cref="ReadCleaned(string, EmbossSequenceQualifiers, EmbossSequenceType)"/> of a member of a <c>seqall</c> stream or a <c>seqset</c> (<paramref name="kind"/>): the reversal is decided once for the whole input
    /// (<paramref name="inputIsNucleotide"/>: <c>acdSetSeqall</c> asks <c>ajSeqIsNuc</c> of the first sequence, <c>acdSetSeqset</c> <c>ajSeqsetIsNuc</c> of the set), and the quirks of the program's reading
    /// of later stream members and of set ranges are those measured on 6.6.0 (see the remarks of <see cref="EmbossInputKind"/>).
    /// </summary>
    internal static EmbossSequenceView ReadCleaned(string text, EmbossSequenceQualifiers qualifiers, EmbossSequenceType type, EmbossInputKind kind, bool? inputIsNucleotide, int index)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(qualifiers);
        if (text.Length == 0)
            throw new ArgumentException("Unable to read sequence: it is empty (EMBOSS reads no residues from it).", nameof(text));

        bool forcedProtein = qualifiers.Protein;
        bool forcedNucleotide = qualifiers.Nucleotide && !qualifiers.Protein;
        bool nucleotideType = type is EmbossSequenceType.Dna or EmbossSequenceType.Nucleotide;
        bool proteinType = type is EmbossSequenceType.Protein or EmbossSequenceType.StopProtein;
        if (nucleotideType && forcedProtein)
            throw new ArgumentException("Sequence is not nucleic (-sprotein on a nucleotide input).", nameof(qualifiers));
        if (proteinType && forcedNucleotide)
            throw new ArgumentException("Sequence is not a protein (-snucleotide on a protein input).", nameof(qualifiers));
        if (nucleotideType && forcedNucleotide)
            text = FixNucleotide(text, type == EmbossSequenceType.Dna);
        else if (nucleotideType && !EmbossSequenceReader.IsNucleotide(text))
            throw new ArgumentException("Sequence is not nucleic (a letter outside the nucleotide codes in a nucleotide input).", nameof(text));

        if (qualifiers.Lower)
            text = text.ToLowerInvariant();
        else if (qualifiers.Upper)
            text = SequenceCase.ToUpperAscii(text);

        // ajSeqTypeCheckIn: passing the check of a nucleotide / protein ACD type sets the sequence type (ajSeqSetNuc / ajSeqSetProt), so it decides ajSeqIsNuc whatever the letters say
        bool isNuc = nucleotideType || (!proteinType && !forcedProtein && (forcedNucleotide || EmbossSequenceReader.IsNucleotide(text)));
        int begin = qualifiers.Begin, end = qualifiers.End;
        bool reversed = false;
        if (qualifiers.Reverse && (inputIsNucleotide ?? isNuc))
        {
            text = ReverseComplement(text);
            (begin, end) = (-end, -begin);      // ajSeqReverseDo, also for the members of a set (ajSeqsetReverse calls it for each; the set's own begin / end are not used by the programs)

            reversed = true;
        }

        int n = text.Length;
        int jbegin = begin > 0 ? begin - 1 : begin;
        int jend = end > 0 ? end - 1 : end;
        int b = 1 + PositionStart(n, 0, jbegin);
        int e = 1 + PositionStart(n, b - 1, jend);
        string residues = text;
        int offset = 0, offend = 0;
        bool trimmed = false;
        if (end != 0)
        {
            if (e >= b)
            {
                offend = n - e;
                residues = residues[..e];
                trimmed = true;
            }
            else
            {
                begin = 0;                      // ajSeqTrim returns before cutting anything
            }
        }

        if (begin != 0 && !(end != 0 && !trimmed))
        {
            residues = residues[(b - 1)..];
            offset = b - 1;
            trimmed = true;
        }

        return new EmbossSequenceView(residues, offset, offend, reversed, text, begin, end, isNuc, trimmed, reversed && kind == EmbossInputKind.Stream && index > 0 && qualifiers.Begin == 0 && qualifiers.End != 0);
    }

    /// <summary><c>ajCvtSposToPosStart</c> (ajmath.c): the 0-based position of <paramref name="position"/> in a sequence of <paramref name="length"/> residues (negative = from the end), clipped into <c>[minimum, length - 1]</c>.</summary>
    public static int PositionStart(int length, int minimum, long position)
    {
        long j = position < 0 ? length + position : position;
        if (j >= length)
            j = length - 1;
        if (j < minimum)
            j = minimum;
        return (int)j;
    }

    /// <summary>
    /// The reverse complement as <c>ajSeqstrReverse</c> makes it: the IUPAC nucleotide complement of <c>ajBaseComp</c> with the case of every letter kept (<c>U</c> and <c>T</c> become <c>A</c>, <c>X</c> becomes <c>N</c>); every other
    /// character is kept. (Measured on 6.6.0 for all 52 letters.)
    /// </summary>
    public static string ReverseComplement(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return string.Create(sequence.Length, sequence, static (dest, src) =>
        {
            for (int i = 0; i < src.Length; i++)
            {
                char c = src[src.Length - 1 - i];
                bool lower = char.IsAsciiLetterLower(c);
                char u = SequenceCase.ToUpperAscii(c);
                char comp = u == 'X' ? 'N' : SequenceExtensions.GetComplementBase(u);
                dest[i] = lower && char.IsAsciiLetterUpper(comp) ? (char)(comp + 32) : comp;
            }
        });
    }

    // seqTypeFix of a nucleotide ACD type under -snucleotide: every character that is not a nucleotide code becomes N, then the type's conversion (dna: U -> T, X -> N; nucleotide: X -> N), case kept for the converted letters.
    private static string FixNucleotide(string text, bool dna)
    {
        return string.Create(text.Length, (text, dna), static (dest, s) =>
        {
            for (int i = 0; i < s.text.Length; i++)
            {
                char c = s.text[i];
                char u = SequenceCase.ToUpperAscii(c);
                bool lower = c != u;
                if (!(u is 'A' or 'C' or 'G' or 'T' or 'U' or 'B' or 'D' or 'H' or 'K' or 'M' or 'N' or 'R' or 'S' or 'V' or 'W' or 'X' or 'Y'))
                {
                    dest[i] = 'N';
                    continue;
                }

                if (u == 'X')
                    u = 'N';
                else if (u == 'U' && s.dna)
                    u = 'T';
                dest[i] = lower ? (char)(u + 32) : u;
            }
        });
    }
}

/// <summary>How an input is read: a single <c>sequence</c>, a <c>seqall</c> stream (needle's <c>bsequence</c>), or a <c>seqset</c> (needleall's <c>asequence</c>).</summary>
/// <remarks>
/// Measured on 6.6.0 (needle / needleall with several sequences and <c>-sreverse2</c> plus a range): the FIRST sequence of a stream is numbered correctly on the reverse strand (<c>10 ... 1</c> for <c>-send2 10</c>), every LATER sequence of the same stream, when only an END is given (<c>-sbegin</c> absent),
/// has its residues cut correctly but is numbered as if it started at the other end (<c>20 ... 11</c>) (the AJAX sequence object is reused and the offsets are not swapped the second time). Sets are reversed by <c>ajSeqsetReverse</c>
/// and decide the reversal once for the whole set from the first member's letters (every member then reverses through <c>ajSeqReverseDo</c>, a protein-typed one included).
/// </remarks>
public enum EmbossInputKind
{
    Sequence,
    Stream,
    Set,
}

/// <summary>The input handling shared by the pairwise EMBOSS programs (needle, water, stretcher, matcher, merger, dotmatcher ...): the first sequence of type <c>any</c>, the second of the type the first one decides.</summary>
public static class EmbossPairInput
{
    /// <summary>
    /// Reads both sequences of a pairwise program as the ACD declares them (<c>asequence</c>: <c>any</c>; <c>bsequence</c>: <c>@($(acdprotein) ? stopprotein : nucleotide)</c>, <c>acdprotein</c> = the first sequence is not nucleotide, forced type included).
    /// </summary>
    /// <returns>The two views and whether the program takes the protein defaults (matrix, penalties).</returns>
    public static (EmbossSequenceView First, EmbossSequenceView Second, bool Protein) Read(EmbossSequence first, EmbossSequence second)
    {
        ArgumentNullException.ThrowIfNull(first);
        ArgumentNullException.ThrowIfNull(second);
        EmbossInputContract.RequireSequence(first, nameof(first));
        EmbossInputContract.RequireSequence(second, nameof(second));
        var v1 = EmbossInputSequence.Read(first.Residues, first.Qualifiers, EmbossSequenceType.Any);
        bool protein = !v1.IsNucleotide;
        var v2 = EmbossInputSequence.Read(second.Residues, second.Qualifiers, protein ? EmbossSequenceType.StopProtein : EmbossSequenceType.Nucleotide);
        return (v1, v2, protein);
    }
}

public static class EmbossSetInput
{
    // ajSeqsetIsNuc: the set's forced nucleotide type or the letters of the FIRST member (ajSeqTypeGapnucS), whatever the ACD type made of each member
    private static bool SetIsNucleotide(IReadOnlyList<EmbossSequence> set, EmbossSequenceType type)
        => set.Count > 0 && (set[0].Qualifiers is { Nucleotide: true, Protein: false } || EmbossSequenceReader.IsNucleotide(EmbossSequenceReader.Clean(set[0].Residues)));

    /// <summary>The program's default matrix: by the type of the first sequence when the qualifier-aware reading decided it (<paramref name="proteinByFirst"/>), else EDNAFULL when every residue string is nucleotide.</summary>
    public static SubstitutionMatrix DefaultMatrix(bool? proteinByFirst, IEnumerable<EmbossSequence> all)
    {
        if (proteinByFirst is { } protein)
            return protein ? SubstitutionMatrix.Blosum62 : SubstitutionMatrix.EdnaFull;
        return all.All(s => EmbossSequenceReader.IsNucleotide(s.Residues)) ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;
    }

    /// <summary>True when any member of the sets carries qualifiers: the qualifier-aware reading is then used (otherwise the programs keep the plain reader).</summary>
    public static bool Any(params IReadOnlyList<EmbossSequence>[] sets) => sets.Any(set => set.Any(s => s.Qualifiers is not null));

    /// <summary>Reads every member of a set (<see cref="EmbossInputSequence.Read"/>): the trimmed sequences (name and description kept, qualifiers dropped) and their views.</summary>
    public static (IReadOnlyList<EmbossSequence> Sequences, IReadOnlyList<EmbossSequenceView> Views) ReadAll(IReadOnlyList<EmbossSequence> set, EmbossSequenceType type, EmbossInputKind kind = EmbossInputKind.Stream)
    {
        bool? inputIsNucleotide = null;
        var views = new List<EmbossSequenceView>(set.Count);
        var seqs = new List<EmbossSequence>(set.Count);
        foreach (var s in set)
        {
            string clean = EmbossSequenceReader.Clean(s.Residues);
            var q = s.Qualifiers ?? EmbossSequenceQualifiers.None;
            if (inputIsNucleotide is null && views.Count == 0)
                inputIsNucleotide = EmbossInputSequence.ReadCleaned(clean, q, type).IsNucleotide;
            var v = EmbossInputSequence.ReadCleaned(clean, q, type, kind, kind == EmbossInputKind.Set ? SetIsNucleotide(set, type) : inputIsNucleotide, views.Count);
            if (v.OffsetsNotSwapped)
                v = v with { Offset = views[0].Offset, Offend = views[0].Offend };        // the stream's sequence object is reused: Offset / Offend keep the values of the first sequence (not recomputed)
            views.Add(v);
            seqs.Add(s with { Residues = v.Residues, Qualifiers = null });
        }

        return (seqs, views);
    }
}
