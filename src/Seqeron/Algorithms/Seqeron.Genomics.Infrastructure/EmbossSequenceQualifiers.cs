namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// The EMBOSS 6.6 sequence-input qualifiers every program accepts for each of its sequence inputs (<c>-sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein</c>;
/// <c>-sformat</c> selects the file format and has no meaning for the residue strings of this API, <c>-sask</c> prompts). Applied by <c>EmbossInputSequence.Read</c> (Alignment assembly)
/// exactly as the AJAX library does (audit A5-7, B09 F172): see that method for the order of the steps and the sources.
/// </summary>
public sealed record EmbossSequenceQualifiers
{
    /// <summary>
    /// The rules by which the qualifier text of the MCP tools is resolved (<see cref="ParseInputs"/>), for the tool descriptions.
    /// </summary>
    public const string CommandLineRules = "The text is resolved exactly as the EMBOSS program resolves it on its command line after the sequence arguments (ACD, 6.6.0): -sbegin1 3, -sbegin1=3, -nosreverse1, -sreverse1 N, any case, -/--/ prefix, "
        + "and every UNIQUE abbreviation (-sbeg1, -sre2; uniqueness is decided among all qualifiers of the program, an ambiguous or unknown one is refused). The number is the sequence (1, and 2 for the second sequence of a two-sequence program; 1 for a program with one sequence). "
        + "An unnumbered qualifier belongs to the LAST sequence given (-sbegin 3 is -sbegin2 3 in a two-sequence program) until a qualifier of another sequence ends that tie; after that it applies to all sequences. "
        + "-sformat -sask -squick -sdbname -sid -iquery -ufo -fformat -fopenfile are accepted without effect (residue strings have no file); any other qualifier or a stray word is refused.";

    /// <summary>No qualifier given (every field at its default): the sequence is used whole, as read.</summary>
    public static EmbossSequenceQualifiers None { get; } = new();

    /// <summary><c>-sbegin</c>: first position to use, 1-based; a negative value counts from the end (-1 = the last residue); 0 = the start.</summary>
    public int Begin { get; init; }

    /// <summary><c>-send</c>: last position to use, 1-based; a negative value counts from the end (-1 = the last residue); 0 = the end.</summary>
    public int End { get; init; }

    /// <summary><c>-sformat</c>: the input sequence format as typed (read by <see cref="EmbossEntryReader"/>: fasta, genbank, embl); null = the program finds the format itself.</summary>
    public string? Format { get; init; }

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
    /// The qualifiers of one input from a command-line fragment such as <c>-sbegin1 3 -send1 20 -sreverse1 -slower1</c>, resolved as the ACD library of the EMBOSS programs does it (see <see cref="ParseInputs"/>) for a program
    /// that has no other qualifiers: <paramref name="suffix"/> is the number of the input (<c>1</c>, <c>2</c>; empty for a program with one sequence input). Use <see cref="ParseInputs"/> with the program name when the program is known
    /// (an abbreviation such as <c>-se</c> is only valid when it is unique among ALL qualifiers of the program).
    /// </summary>
    /// <exception cref="ArgumentException">Unknown / ambiguous qualifier, missing or invalid value, a token that is not a qualifier, a qualifier that is not a sequence qualifier, an input number that does not exist.</exception>
    public static EmbossSequenceQualifiers Parse(string? commandLine, string suffix = "")
    {
        int index = 0, inputs = 1;
        if (suffix.Length > 0)
        {
            if (!int.TryParse(suffix, System.Globalization.NumberStyles.None, System.Globalization.CultureInfo.InvariantCulture, out int n) || n is < 1 or > 2)
                throw new ArgumentException($"Input number '{suffix}' must be 1 or 2.", nameof(suffix));
            index = n - 1;
            inputs = 2;
        }

        return Resolve(EmbossAcdProfiles.SequenceOnly(inputs), commandLine)[index];
    }

    /// <summary>
    /// The sequence qualifiers of EVERY sequence input of an EMBOSS program from the command-line text that follows its sequence arguments (<c>needle a.fa b.fa -sbegin 3 -send2 5</c>: the text is <c>-sbegin 3 -send2 5</c>), as the ACD library
    /// resolves them (audit A6-2, B09 F177; <see cref="AcdCommandLine"/> for the rules): both <c>-sbegin1 5</c> (also for a program with one input) and <c>-sbegin=5</c>, <c>-nosreverse2</c>, <c>-sreverse2 N</c>, any unique abbreviation of a
    /// qualifier name (<c>-sbeg1</c>, <c>-sre2</c>; ambiguity is decided over all qualifiers of the program), and an unnumbered qualifier, which belongs to the sequence argument that precedes it (the LAST one, <c>-sbegin 3</c> after two files is
    /// <c>-sbegin2 3</c>) until a qualifier of another input or of the program ends that tie, and then applies to all inputs. <c>-sformat</c> / <c>-sask</c> / <c>-squick</c> / <c>-sdbname</c> / <c>-sid</c> / <c>-iquery</c> / <c>-ufo</c> / <c>-fformat</c> /
    /// <c>-fopenfile</c> are accepted (their values are not interpreted: residue strings have no file format, prompt or database; they are kept as the program accepts them); <c>-ioffset</c> must be 0.
    /// </summary>
    /// <param name="program">The EMBOSS program: needle, water, stretcher, matcher, supermatcher, needleall, merger, megamerger, seqmatchall, fuzznuc, fuzzpro, fuzztran, etandem, equicktandem, dotmatcher, distmat.</param>
    /// <param name="commandLine">The qualifier text; null or blank = none.</param>
    /// <returns>One entry per sequence input of the program (1 or 2), in parameter order.</returns>
    /// <exception cref="ArgumentException">Unknown program; a token the program would refuse (unknown or ambiguous qualifier, missing or invalid value, <c>no</c> prefix with a value) or that is not a qualifier of a sequence input.</exception>
    public static IReadOnlyList<EmbossSequenceQualifiers> ParseInputs(string program, string? commandLine)
    {
        ArgumentNullException.ThrowIfNull(program);
        return Resolve(EmbossAcdProfiles.Get(program), commandLine);
    }

    /// <summary>
    /// <see cref="ParseInputs"/> for a text that may also carry the qualifiers of the program's REPORT output (<c>-rformat -rname -rextension -rdirectory -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>, with the
    /// output's number or without: <c>-rformat2 table</c>, <c>-rformat table</c>; audit A5-8b, B09 F184); every other qualifier is refused as in <see cref="ParseInputs"/>.
    /// </summary>
    /// <param name="program">The EMBOSS program.</param>
    /// <param name="commandLine">The qualifier text; null or blank = none.</param>
    /// <returns>The sequence qualifiers per input and the report qualifiers the text set.</returns>
    public static (IReadOnlyList<EmbossSequenceQualifiers> Inputs, EmbossReportQualifiers Report) ParseInputsAndReport(string program, string? commandLine)
    {
        ArgumentNullException.ThrowIfNull(program);
        var profile = EmbossAcdProfiles.Get(program);
        var inputs = Resolve(profile, commandLine, out var report);
        return (inputs, report);
    }

    private static EmbossSequenceQualifiers[] Resolve(AcdProfile profile, string? commandLine) => Resolve(profile, commandLine, out _, allowReport: false);

    private static EmbossSequenceQualifiers[] Resolve(AcdProfile profile, string? commandLine, out EmbossReportQualifiers report, bool allowReport = true)
    {
        report = new EmbossReportQualifiers();
        var result = new EmbossSequenceQualifiers[profile.Inputs.Count];
        Array.Fill(result, None);
        if (string.IsNullOrWhiteSpace(commandLine))
            return result;
        var tokens = commandLine.Split([' ', '\t', '\r', '\n'], StringSplitOptions.RemoveEmptyEntries);
        var line = new AcdCommandLine(profile);
        line.Run(tokens);

        var sequenceItems = new HashSet<AcdItem>(profile.Inputs.SelectMany(i => i.Assoc));
        foreach (var (item, raw) in line.Assignments)
        {
            if (allowReport && item.Master is { Name: "outfile" } && ReportItems.Contains(item.Name))
            {
                report = ApplyReport(report, item, raw, commandLine);
                continue;
            }

            if (!sequenceItems.Contains(item))
                throw new ArgumentException($"-{item.Name} is not a sequence qualifier: this text takes the qualifiers of the sequence inputs only (-sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein, with the input number, e.g. -sbegin1).", nameof(commandLine));
        }

        for (int k = 0; k < result.Length; k++)
        {
            var q = None;
            foreach (var item in profile.Inputs[k].Assoc)
            {
                string? raw = line.Value(item);
                if (raw is null)
                    continue;
                q = Apply(q, item, raw, commandLine);
            }

            result[k] = q;
        }

        return result;
    }

    private static readonly HashSet<string> ReportItems = ["rformat", "rname", "rextension", "rdirectory", "raccshow", "rdesshow", "rscoreshow", "rstrandshow", "rusashow", "rmaxall", "rmaxseq"];

    private static EmbossReportQualifiers ApplyReport(EmbossReportQualifiers q, AcdItem item, string raw, string commandLine)
    {
        if (item.Type == "integer")
        {
            if (!AcdCommandLine.TryToInt(raw, out int number))
                throw new ArgumentException($"Invalid integer value '{raw}' for -{item.Name}.", nameof(commandLine));
            if (number < 0)
                throw new ArgumentException($"-{item.Name} must be 0 or more (0 = no limit).", nameof(commandLine));
            return item.Name == "rmaxall" ? q with { MaxAll = number } : q with { MaxSeq = number };
        }

        if (item.Type == "boolean")
        {
            if (!AcdCommandLine.TryToBool(raw, out bool flag))
                throw new ArgumentException($"Invalid Y/N value '{raw}' for -{item.Name}.", nameof(commandLine));
            return item.Name switch
            {
                "raccshow" => q with { ShowAccession = flag },
                "rdesshow" => q with { ShowDescription = flag },
                "rscoreshow" => q with { ShowScore = flag },
                "rstrandshow" => q with { ShowStrand = flag },
                _ => q with { ShowUsa = flag },
            };
        }

        return item.Name == "rformat" ? q with { Format = raw } : q;            // rname / rextension / rdirectory: the file name is the caller's business
    }

    private static EmbossSequenceQualifiers Apply(EmbossSequenceQualifiers q, AcdItem item, string raw, string commandLine)
    {
        if (item.Type == "integer")
        {
            if (!AcdCommandLine.TryToInt(raw, out int number))
                throw new ArgumentException($"Invalid integer value '{raw}' for -{item.Name}{(item.PNum > 0 ? item.PNum.ToString(System.Globalization.CultureInfo.InvariantCulture) : "")}.", nameof(commandLine));
            switch (item.Name)
            {
                case "sbegin":
                    return q with { Begin = number };
                case "send":
                    return q with { End = number };
                case "ioffset" when number != 0:
                    throw new ArgumentException("-ioffset must be 0: a residue string has no file to start reading at an offset.", nameof(commandLine));
                default:
                    return q;
            }
        }

        if (item.Type == "boolean")
        {
            if (!AcdCommandLine.TryToBool(raw, out bool flag))
                throw new ArgumentException($"Invalid Y/N value '{raw}' for -{item.Name}{(item.PNum > 0 ? item.PNum.ToString(System.Globalization.CultureInfo.InvariantCulture) : "")}.", nameof(commandLine));
            return item.Name switch
            {
                "sreverse" => q with { Reverse = flag },
                "slower" => q with { Lower = flag },
                "supper" => q with { Upper = flag },
                "scircular" => q with { Circular = flag },
                "snucleotide" => q with { Nucleotide = flag },
                "sprotein" => q with { Protein = flag },
                _ => q,                                    // sask, squick: accepted, no effect on a residue string
            };
        }

        return item.Name == "sformat" ? q with { Format = raw } : q;      // iquery, sdbname, sid, ufo, fformat, fopenfile: accepted, not interpreted
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

/// <summary>The qualifiers of the report output of an EMBOSS program (<c>-rformat -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>) as resolved from a command-line text; null / default = the program's own default.</summary>
/// <param name="Format"><c>-rformat</c> as typed (null = the program's default format).</param>
/// <param name="ShowAccession"><c>-raccshow</c>.</param>
/// <param name="ShowDescription"><c>-rdesshow</c>.</param>
/// <param name="ShowScore"><c>-rscoreshow</c> (null = the default, Y).</param>
/// <param name="ShowStrand"><c>-rstrandshow</c> (null = the default, Y).</param>
/// <param name="ShowUsa"><c>-rusashow</c>.</param>
/// <param name="MaxAll"><c>-rmaxall</c> (0 = no limit).</param>
/// <param name="MaxSeq"><c>-rmaxseq</c> (0 = no limit).</param>
public sealed record EmbossReportQualifiers(string? Format = null, bool ShowAccession = false, bool ShowDescription = false, bool? ShowScore = null, bool? ShowStrand = null, bool ShowUsa = false, int MaxAll = 0, int MaxSeq = 0);
