namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// The EMBOSS 6.6 sequence-input qualifiers every program accepts for each of its sequence inputs (<c>-sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein</c>;
/// <c>-sformat</c> selects the file format and has no meaning for the residue strings of this API, <c>-sask</c> prompts). Applied by <c>EmbossInputSequence.Read</c> (Alignment assembly)
/// exactly as the AJAX library does (audit A5-7, B09 F172): see that method for the order of the steps and the sources.
/// </summary>
public sealed record EmbossSequenceQualifiers
{
    /// <summary>No qualifier given (every field at its default): the sequence is used whole, as read.</summary>
    public static EmbossSequenceQualifiers None { get; } = new();

    /// <summary><c>-sbegin</c>: first position to use, 1-based; a negative value counts from the end (-1 = the last residue); 0 = the start.</summary>
    public int Begin { get; init; }

    /// <summary><c>-send</c>: last position to use, 1-based; a negative value counts from the end (-1 = the last residue); 0 = the end.</summary>
    public int End { get; init; }

    /// <summary><c>-sreverse</c>: use the reverse complement (only for a nucleotide sequence; a protein is left alone). The range then refers to the ORIGINAL sequence.</summary>
    public bool Reverse { get; init; }

    /// <summary><c>-slower</c>: make the residues lower case.</summary>
    public bool Lower { get; init; }

    /// <summary><c>-supper</c>: make the residues upper case (<see cref="Lower"/> wins when both are set, as in the program).</summary>
    public bool Upper { get; init; }

    /// <summary><c>-scircular</c>: the sequence is circular. A flag of the sequence object that none of the sixteen alignment / pattern / repeat programs reads (measured on 6.6.0: the output is identical with and without it); kept for the record.</summary>
    public bool Circular { get; init; }

    /// <summary><c>-snucleotide</c>: the sequence is nucleotide whatever its letters say.</summary>
    public bool Nucleotide { get; init; }

    /// <summary><c>-sprotein</c>: the sequence is protein whatever its letters say (wins over <see cref="Nucleotide"/>, as in <c>ajSeqTypeCheckIn</c>).</summary>
    public bool Protein { get; init; }

    /// <summary>
    /// The qualifiers of one input from a command-line fragment such as <c>-sbegin1 3 -send1 20 -sreverse1 -slower1</c>: <paramref name="suffix"/> is the input's number (<c>1</c>, <c>2</c>; empty for a program with one sequence input).
    /// Qualifiers with another suffix are ignored (one text can describe all the inputs of a program); a boolean qualifier may be written <c>-sreverse1</c> or <c>-nosreverse1</c>; <c>-sbegin</c> / <c>-send</c> take an integer
    /// (negative counts from the end). <c>-sformat</c> / <c>-sask</c> / <c>-sdbname</c> ... are accepted and ignored (they select the file format and prompting, which residue strings do not have).
    /// </summary>
    /// <exception cref="ArgumentException">Unknown token, missing or non-integer value.</exception>
    public static EmbossSequenceQualifiers Parse(string? commandLine, string suffix = "")
    {
        var q = None;
        if (string.IsNullOrWhiteSpace(commandLine))
            return q;
        var tokens = commandLine.Split([' ', '\t', '\r', '\n'], StringSplitOptions.RemoveEmptyEntries);
        int next = 0;
        while (next < tokens.Length)
        {
            string t = tokens[next++];
            if (!t.StartsWith('-'))
                throw new ArgumentException($"Unexpected token '{t}': a sequence qualifier starts with '-' (for example -sbegin1 3 -sreverse1).", nameof(commandLine));
            string name = t.TrimStart('-').ToLowerInvariant();
            bool negated = name.StartsWith("nos", StringComparison.Ordinal) && name.Length > 3;
            if (negated)
                name = name[2..];
            string stem = name.TrimEnd('0', '1', '2', '3', '4', '5', '6', '7', '8', '9');
            bool mine = name[stem.Length..] == suffix;
            switch (stem)
            {
                case "sbegin" or "send":
                    {
                        if (negated || next >= tokens.Length || !int.TryParse(tokens[next], System.Globalization.NumberStyles.AllowLeadingSign, System.Globalization.CultureInfo.InvariantCulture, out int value))
                            throw new ArgumentException($"{t} needs an integer value.", nameof(commandLine));
                        next++;
                        if (mine)
                            q = stem == "sbegin" ? q with { Begin = value } : q with { End = value };
                        break;
                    }
                case "sreverse":
                    q = mine ? q with { Reverse = !negated } : q;
                    break;
                case "slower":
                    q = mine ? q with { Lower = !negated } : q;
                    break;
                case "supper":
                    q = mine ? q with { Upper = !negated } : q;
                    break;
                case "scircular":
                    q = mine ? q with { Circular = !negated } : q;
                    break;
                case "snucleotide":
                    q = mine ? q with { Nucleotide = !negated } : q;
                    break;
                case "sprotein":
                    q = mine ? q with { Protein = !negated } : q;
                    break;
                case "sformat" or "sdbname" or "sid" or "sopenfile" or "ssubdb" or "iquery" or "ioffset":
                    next++;      // takes a value that residue strings have no use for
                    break;
                case "squick" or "sask":
                    break;
                default:
                    throw new ArgumentException($"Unknown sequence qualifier '{t}'. Supported: -sbegin, -send, -sreverse, -slower, -supper, -scircular, -snucleotide, -sprotein (with the input number, e.g. -sbegin1).", nameof(commandLine));
            }
        }

        return q;
    }

    /// <summary>True when no qualifier differs from the default.</summary>
    public bool IsDefault => Begin == 0 && End == 0 && !Reverse && !Lower && !Upper && !Circular && !Nucleotide && !Protein;
}

/// <summary>The ACD sequence types of the programs (<c>acdSetSeq</c> type attribute): the check and the repair <c>ajSeqTypeCheckIn</c> applies to an input sequence.</summary>
public enum EmbossSequenceType
{
    /// <summary><c>any</c>: no check (the first sequence of the alignment programs).</summary>
    Any,
    /// <summary><c>gapany</c>: any sequence, gap characters kept (<c>distmat</c>, <c>consambig</c>, <c>showalign</c>).</summary>
    GapAny,
    /// <summary><c>dna</c>: a nucleotide sequence; <c>U</c> reads as <c>T</c>, <c>X</c> and <c>?</c> as <c>N</c> (<c>fuzznuc</c>, <c>fuzztran</c>, <c>etandem</c>, <c>equicktandem</c>, <c>megamerger</c>).</summary>
    Dna,
    /// <summary><c>nucleotide</c>: a nucleotide sequence (the second sequence of an alignment whose first sequence is not protein).</summary>
    Nucleotide,
    /// <summary><c>protein</c>: a protein sequence (<c>fuzzpro</c>).</summary>
    Protein,
    /// <summary><c>stopprotein</c>: a protein sequence with possible stops (the second sequence of an alignment whose first sequence is protein).</summary>
    StopProtein,
}

/// <summary>
/// A sequence after the EMBOSS reader and the range / orientation of its qualifiers (<c>EmbossInputSequence.Read</c> (Alignment assembly)): <c>ajSeqTrim</c> state.
/// </summary>
/// <param name="Residues">The residues the program works on: the sequence as read, case-changed, reverse-complemented and cut to the range (<c>ajSeqGetSeqC</c> after <c>ajSeqTrim</c>).</param>
/// <param name="Offset">Residues cut from the start of the (reversed) sequence (<c>ajSeqGetOffset</c>: the program adds it to every position it reports).</param>
/// <param name="Offend">Residues cut from the end (<c>ajSeqGetOffend</c>).</param>
/// <param name="Reversed">The sequence was reverse-complemented (<c>-sreverse</c> on a nucleotide sequence): reported positions then count down from the original end.</param>
/// <param name="Full">The sequence before the range was applied (as read, case-changed, reverse-complemented): <c>Offset + Residues.Length + Offend</c> residues.</param>
/// <param name="Begin">The range start as <c>ajSeqTrim</c> saw it after the reversal (<c>seq-&gt;Begin</c>: 1-based, negative = from the end, 0 = unset).</param>
/// <param name="End">The range end as <c>ajSeqTrim</c> saw it after the reversal.</param>
/// <param name="IsNucleotide"><c>ajSeqIsNuc</c> of the sequence as the ACD processing saw it (before the range): the forced type, else the letters.</param>
/// <param name="OffsetsNotSwapped">A reversed LATER member of a <c>seqall</c> stream: the program numbers it from the other end (the offsets are not swapped when the report counts the original strand), see <c>EmbossInputKind</c>.</param>
/// <param name="Trimmed">The range was applied (a begin or an end was given and the range was not empty after clipping, <c>ajSeqTrim</c> does nothing otherwise).</param>
public sealed record EmbossSequenceView(string Residues, int Offset, int Offend, bool Reversed, string Full, int Begin, int End, bool IsNucleotide, bool Trimmed, bool OffsetsNotSwapped = false)
{
    /// <summary>
    /// The 1-based inclusive range a program reports for the window <c>Residues[position .. position + length)</c> (0-based <paramref name="position"/> in <see cref="Residues"/>): <c>Offset + position + 1 .. Offset + position + length</c>
    /// on the sequence as used, and on the ORIGINAL strand of a reversed sequence <c>start &gt; end</c> is avoided by listing the lower position first (<c>Full.Length - Offset - position - length + 1 .. Full.Length - Offset - position</c>).
    /// </summary>
    /// <param name="position">0-based start of the window in <see cref="Residues"/>.</param>
    /// <param name="length">Length of the window.</param>
    public (int Start, int End) ReportedRange(int position, int length)
    {
        if (Reversed)
        {
            int end = Full.Length - Offset - position;
            return (end - length + 1, end);
        }

        return (Offset + position + 1, Offset + position + length);
    }
}
