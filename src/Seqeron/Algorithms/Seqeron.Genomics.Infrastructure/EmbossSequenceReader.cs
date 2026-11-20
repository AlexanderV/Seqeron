namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// The one canonical cleaner of the EMBOSS ports (audit A4-3, B09 F153): how EMBOSS 6.6 hands a PLAIN input sequence to a program, shared by every port of this repository (needle, water, stretcher, matcher, supermatcher,
/// needleall, merger, megamerger, seqmatchall, fuzznuc, fuzzpro, fuzztran, etandem, equicktandem, dotmatcher, palindrome ...; the alignment readers of consambig / showalign / cons use <see cref="CleanGapped(string)"/>).
/// </summary>
/// <remarks>
/// <para>
/// Source (EMBOSS 6.6 <c>ajax</c>, read at raw.githubusercontent.com/kimrutherford/EMBOSS): the file reader <c>seqAppend</c> (<c>ajseqread.c</c>) keeps the letters and <c>* . ~ ? # + -</c> of every line
/// (<c>ajStrKeepSetAlphaC(.., "*.~?#+-")</c>: digits, blanks, quotes, <c>_</c>, control characters and every non-ASCII byte vanish silently), and the sequence-type check
/// <c>ajSeqTypeCheckIn</c> (<c>ajseqtype.c</c>) of every UNGAPPED ACD type (<c>any</c>, <c>dna</c>, <c>nucleotide</c>, <c>protein</c>, <c>stopprotein</c>, ...; all the programs above) calls <c>ajStrRemoveGap</c>
/// (<c>ajstr.c</c>), which keeps exactly <c>A-Z a-z *</c>. So a plain sequence reaches the program as its letters and <c>*</c>: <c>-</c> <c>.</c> <c>~</c> <c>?</c> <c>#</c> <c>+</c>, digits, blanks and
/// <c>_</c> are gone BEFORE any length check, <c>-sbegin</c> offset, coordinate or score is computed, and every position a port reports counts the cleaned sequence. The case of the letters is untouched
/// (the programs fold it later themselves). Measured on the installed 6.6.0 with a character inserted into a sequence of every program above (all 252 ASCII / high bytes other than newline, carriage return and
/// <c>&gt;</c>, both sequences, nucleotide and protein; scripts/analysis/generate_emboss_gap_oracle.py): the output is identical to the one without the character for exactly the 199 characters this class removes.
/// A sequence that is empty afterwards is "Unable to read sequence" (zero length) in every program.
/// </para>
/// <para>
/// The raw route stays available: the non-EMBOSS-named entry points (<c>GlobalAlign</c>, <c>LocalAlign</c>, <c>FindMotif</c>, ...) keep every character literally.
/// </para>
/// </remarks>
public static class EmbossSequenceReader
{
    /// <summary>
    /// <c>ajStrRemoveGap</c> as applied by <c>ajSeqTypeCheckIn</c> to an ungapped type: keeps <c>A-Z</c>, <c>a-z</c> and <c>*</c> and drops every other character (<c>- . ~ ? # +</c>, digits, blanks, <c>_</c>,
    /// control characters, non-ASCII letters). Returns the argument itself when nothing is removed.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static string Clean(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        int i = 0;
        while (i < sequence.Length && IsPlainSequenceCharacter(sequence[i]))
            i++;
        if (i == sequence.Length)
            return sequence;

        var chars = new char[sequence.Length];
        sequence.CopyTo(0, chars, 0, i);
        int n = i;
        for (i++; i < sequence.Length; i++)
            if (IsPlainSequenceCharacter(sequence[i]))
                chars[n++] = sequence[i];
        return new string(chars, 0, n);
    }

    /// <summary>
    /// The reading of an alignment row (<c>ajseqread.c seqAppend</c> then <c>ajSeqTypeCheckIn</c> of the gapped type <c>gapany</c> of consambig / showalign / cons / distmat): the letters and the gap characters <c>- . ~</c> and <c>*</c> are
    /// kept, every other character (digits, blanks, <c>_</c>, punctuation, control characters, non-ASCII) is dropped, and the characters <c>#</c> and <c>+</c>, which the reader keeps, read exactly like <c>?</c> (the unknown
    /// residue: measured on consambig 6.6.0, the rows <c>ACGT#CGT</c>, <c>ACGT+CGT</c> and <c>ACGT?CGT</c> give the same consensus and the same set type), so they are returned as <c>?</c>. The later conversion of
    /// <c>. ~</c> to <c>-</c> and of <c>?</c> to <c>N</c> / <c>X</c> is the alignment normaliser's. Non-ASCII bytes are dropped here as in the C locale; the binary's <c>isalpha</c> keeps some of them under other locales (not reproducible).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="row"/> is null.</exception>
    public static string CleanGapped(string row)
    {
        ArgumentNullException.ThrowIfNull(row);
        int i = 0;
        while (i < row.Length && IsKeptByGappedReader(row[i]) && row[i] is not ('#' or '+'))
            i++;
        if (i == row.Length)
            return row;

        var chars = new char[row.Length];
        row.CopyTo(0, chars, 0, i);
        int n = i;
        for (; i < row.Length; i++)
        {
            char c = row[i];
            if (IsKeptByGappedReader(c))
                chars[n++] = c is '#' or '+' ? '?' : c;
        }

        return new string(chars, 0, n);
    }

    /// <summary><see cref="Clean(string)"/> of every residue string of a set (<see cref="EmbossSequence"/>: the reader's sequences of needleall / supermatcher); names and descriptions are kept, the list itself when nothing changes.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequences"/> is null.</exception>
    public static IReadOnlyList<EmbossSequence> Clean(IReadOnlyList<EmbossSequence> sequences)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        List<EmbossSequence>? cleaned = null;
        for (int i = 0; i < sequences.Count; i++)
        {
            string residues = Clean(sequences[i].Residues);
            if (cleaned is null && ReferenceEquals(residues, sequences[i].Residues))
                continue;
            cleaned ??= [.. sequences.Take(i)];
            cleaned.Add(sequences[i] with { Residues = residues });
        }

        if (cleaned is null)
            return sequences;
        for (int i = cleaned.Count; i < sequences.Count; i++)
            cleaned.Add(sequences[i]);
        return cleaned;
    }

    /// <summary>A character <see cref="Clean(string)"/> keeps: an ASCII letter or <c>*</c> (<c>ajStrRemoveGap</c>).</summary>
    public static bool IsPlainSequenceCharacter(char c) => char.IsAsciiLetter(c) || c == '*';

    /// <summary>A character the EMBOSS file reader keeps (<c>seqAppend</c>): an ASCII letter or one of <c>* . ~ ? # + -</c>.</summary>
    public static bool IsKeptByGappedReader(char c) => char.IsAsciiLetter(c) || c is '*' or '.' or '~' or '?' or '#' or '+' or '-';

    /// <summary>
    /// <c>ajSeqIsNuc</c> of an untyped sequence (<c>ajSeqTypeGapnucS</c>, <c>ajseqtype.c</c>: <c>ajStrIsCharsetCaseS</c> over <c>seqCharNuc + seqCharGap</c>): every character, in either case, is one of the
    /// nucleotide codes <c>ACGTUBDHKMNRSVWXY</c>, <c>?</c> or a gap <c>. ~ -</c> (an empty sequence is nucleotide). The single definition for the programs that choose their default matrix or alphabet by it
    /// (needle, water, stretcher, dotmatcher, merger, consambig, showalign, distmat ...).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static bool IsNucleotide(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        foreach (char ch in sequence)
        {
            switch (SequenceCase.ToUpperAscii(ch))
            {
                case 'A' or 'C' or 'G' or 'T' or 'U' or 'B' or 'D' or 'H' or 'K' or 'M' or 'N' or 'R' or 'S' or 'V' or 'W' or 'X' or 'Y' or '?' or '.' or '~' or '-':
                    break;
                default:
                    return false;
            }
        }

        return true;
    }
}
