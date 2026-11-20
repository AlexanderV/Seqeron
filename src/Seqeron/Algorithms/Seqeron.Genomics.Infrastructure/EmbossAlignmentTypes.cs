namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// The EMBOSS 6.6 alignment output formats (<c>-aformat</c>, <c>ajalign.c</c> <c>alignFormat[]</c>) written by <c>EmbossAlignmentWriter</c>. The names are the EMBOSS names
/// (<see cref="Paup"/> / <see cref="Paupnon"/> are the aliases of <see cref="Nexus"/> / <see cref="Nexusnon"/>, <see cref="A2m"/> the alias of <see cref="Fasta"/>);
/// <see cref="Unknown"/>, <see cref="Multiple"/>, <see cref="Pair"/> and <see cref="Simple"/> share one layout.
/// </summary>
public enum EmbossAlignFormat
{
    /// <summary>EMBOSS <c>unknown</c> (written as <c>simple</c>).</summary>
    Unknown,
    /// <summary>FASTA, 60 columns (<c>fasta</c>).</summary>
    Fasta,
    /// <summary>A2M (<c>a2m</c>; the FASTA writer).</summary>
    A2m,
    /// <summary>GCG MSF (<c>msf</c>).</summary>
    Msf,
    /// <summary>ClustalW 1.83 style (<c>clustal</c>).</summary>
    Clustal,
    /// <summary>MEGA interleaved (<c>mega</c>).</summary>
    Mega,
    /// <summary>MEGA non-interleaved (<c>meganon</c>).</summary>
    Meganon,
    /// <summary>NEXUS interleaved (<c>nexus</c>).</summary>
    Nexus,
    /// <summary>Alias of <see cref="Nexus"/> (<c>paup</c>).</summary>
    Paup,
    /// <summary>NEXUS non-interleaved (<c>nexusnon</c>).</summary>
    Nexusnon,
    /// <summary>Alias of <see cref="Nexusnon"/> (<c>paupnon</c>).</summary>
    Paupnon,
    /// <summary>PHYLIP interleaved (<c>phylip</c>).</summary>
    Phylip,
    /// <summary>PHYLIP sequential (<c>phylipnon</c>).</summary>
    Phylipnon,
    /// <summary>SELEX (<c>selex</c>).</summary>
    Selex,
    /// <summary>Treecon (<c>treecon</c>).</summary>
    Treecon,
    /// <summary>Pearson MARKX0 (<c>markx0</c>).</summary>
    Markx0,
    /// <summary>Pearson MARKX1.</summary>
    Markx1,
    /// <summary>Pearson MARKX2.</summary>
    Markx2,
    /// <summary>Pearson MARKX3.</summary>
    Markx3,
    /// <summary>Pearson MARKX10.</summary>
    Markx10,
    /// <summary>Start and end of the match between a pair (<c>match</c>).</summary>
    Match,
    /// <summary>Simple multiple alignment (<c>multiple</c>; the <c>simple</c> layout).</summary>
    Multiple,
    /// <summary>Simple pairwise alignment (<c>pair</c>).</summary>
    Pair,
    /// <summary>Simple multiple alignment (<c>simple</c>).</summary>
    Simple,
    /// <summary>SAM, one record per pair (<c>sam</c>).</summary>
    Sam,
    /// <summary>Score values for pairs (<c>score</c>).</summary>
    Score,
    /// <summary>Simple multiple sequence format for SRS (<c>srs</c>).</summary>
    Srs,
    /// <summary>Simple pairwise sequence format for SRS (<c>srspair</c>; the default of needle and water).</summary>
    SrsPair,
    /// <summary>T-Coffee library (<c>tcoffee</c>).</summary>
    TCoffee,
}

/// <summary>
/// One aligned sequence of an EMBOSS alignment (EMBOSS <c>AlignData</c>): the gapped row and where it sits in the original sequence.
/// </summary>
/// <param name="Name">The sequence name (the first word of the FASTA header).</param>
/// <param name="Row">The aligned row, <c>-</c> (also <c>.</c> <c>~</c> <c>?</c> and space) for a gap.</param>
/// <param name="Offset">Residues of the original sequence before the first residue of <paramref name="Row"/> (water: the local start − 1; <c>-sbegin</c> − 1).</param>
/// <param name="Offend">Residues of the original sequence after the last residue of <paramref name="Row"/>.</param>
/// <param name="Reversed">The sequence was reverse-complemented (<c>-sreverse</c>).</param>
/// <param name="Accession">Accession number (shown by <c>-aaccshow</c>; written after the name by the FASTA writer).</param>
/// <param name="Description">Description (the rest of the FASTA header; shown by <c>-adesshow</c>).</param>
/// <param name="Usa">The full USA (<c>fasta::file:name</c>, shown instead of the name by <c>-ausashow</c>).</param>
/// <param name="Weight">Sequence weight (MSF <c>Weight:</c>, SELEX <c>#=SQ</c>), 1.0 for EMBOSS alignments.</param>
public sealed record EmbossAlignedSequence(
    string Name,
    string Row,
    int Offset = 0,
    int Offend = 0,
    bool Reversed = false,
    string Accession = "",
    string Description = "",
    string Usa = "",
    double Weight = 1.0)
{
    /// <summary>The sequence as aligned after the range / orientation of <paramref name="view"/> (A5-7): the cut residues join the offsets, the reverse-complement flag is set.</summary>
    public EmbossAlignedSequence WithView(EmbossSequenceView view)
    {
        ArgumentNullException.ThrowIfNull(view);
        // the report counts residues in the ORIGINAL strand (ajalign.c alignSeqBegin: a reversed sequence starts at Len - Offend - Start + 1): when the sequence is reversed, the residues before the aligned part in the
        // reverse-complement are the residues AFTER it in the original, so the two offsets swap sides (ajSeqReverseDo swaps seq->Offset and seq->Offend).
        return view.Reversed && !view.OffsetsNotSwapped
            ? this with { Offset = Offend + view.Offend, Offend = Offset + view.Offset, Reversed = true }
            : this with { Offset = Offset + view.Offset, Offend = Offend + view.Offend, Reversed = Reversed || view.Reversed };
    }
}

/// <summary>One pairwise alignment as reported by needle / water / stretcher / matcher / supermatcher.</summary>
/// <param name="First">The first sequence (<c>-asequence</c>).</param>
/// <param name="Second">The second sequence (<c>-bsequence</c>).</param>
/// <param name="Score">The alignment score (printed with up to three decimals, at least one).</param>
/// <param name="ExtraSubHeader">Extra sub-header text after the standard Length / Identity / Similarity / Gaps / Score block (needle <c>-nobrief</c>: the <c>Longest_Identity</c> lines).</param>
public sealed record EmbossPairAlignment(
    EmbossAlignedSequence First,
    EmbossAlignedSequence Second,
    double Score,
    string ExtraSubHeader = "")
{
    /// <summary>
    /// The EMBOSS view of an <see cref="AlignmentResult"/> of two named sequences: the rows, the residues before (<see cref="AlignmentResult.StartPosition1"/>) and after (length − <see cref="AlignmentResult.EndPosition1"/>) the aligned part
    /// of each sequence, and <see cref="AlignmentResult.ExactScore"/>.
    /// </summary>
    public static EmbossPairAlignment From(AlignmentResult alignment, string name1, int length1, string name2, int length2)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        return new EmbossPairAlignment(
            new EmbossAlignedSequence(name1, alignment.AlignedSequence1, alignment.StartPosition1, Math.Max(0, length1 - alignment.EndPosition1)),
            new EmbossAlignedSequence(name2, alignment.AlignedSequence2, alignment.StartPosition2, Math.Max(0, length2 - alignment.EndPosition2)),
            alignment.ExactScore);
    }
    /// <summary>
    /// <see cref="From(AlignmentResult, string, int, string, int)"/> for sequences read with the sequence qualifiers (<c>-sbegin -send -sreverse</c>, A5-7): the alignment was computed on <see cref="EmbossSequenceView.Residues"/>, so the residues
    /// the range cut from each side (<see cref="EmbossSequenceView.Offset"/> / <see cref="EmbossSequenceView.Offend"/>) are added to the offsets and the sequence is marked <see cref="EmbossAlignedSequence.Reversed"/> when it was reverse-complemented.
    /// </summary>
    public static EmbossPairAlignment From(AlignmentResult alignment, string name1, EmbossSequenceView view1, string name2, EmbossSequenceView view2)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        ArgumentNullException.ThrowIfNull(view1);
        ArgumentNullException.ThrowIfNull(view2);
        var pair = From(alignment, name1, view1.Residues.Length, name2, view2.Residues.Length);
        return pair with { First = pair.First.WithView(view1), Second = pair.Second.WithView(view2) };
    }
}

/// <summary>The program context an EMBOSS alignment file reports: command line, date, matrix, gap penalties and the <c>-a*</c> qualifiers.</summary>
public sealed record EmbossAlignOptions
{
    /// <summary>The program name (<c>needle</c>, <c>water</c>, <c>stretcher</c>, ...), printed in the header and the SAM <c>@PG</c> line.</summary>
    public string Program { get; init; } = "needle";

    /// <summary>The EMBOSS version (the 6.6.0 binaries print <c>6.6.0.0</c>).</summary>
    public string ProgramVersion { get; init; } = string.Join('.', 6, 6, 0, 0);

    /// <summary>The command-line arguments after the program name, one per header line (<c>#    -asequence a.fa</c>: a qualifier with its value is one entry); the SAM <c>CL:</c> field joins them by blanks.</summary>
    public IReadOnlyList<string> CommandLine { get; init; } = [];

    /// <summary>The run date (header <c>Rundate</c>, MSF / MEGA / NEXUS dates). EMBOSS uses the local time of the machine.</summary>
    public DateTime RunDate { get; init; } = new(2000, 1, 1, 0, 0, 0, DateTimeKind.Unspecified);

    /// <summary>The report file name (<c>Report_file:</c> and the MSF header): <c>stdout</c> for <c>-stdout</c>.</summary>
    public string ReportFile { get; init; } = "stdout";

    /// <summary>The sequence type of the output set when the program has declared one (<c>seqret</c> sets <c>N</c> or <c>P</c> from the input; the alignment programs leave it unset, so MSF guesses it from the first row and NEXUS always writes <c>datatype=DNA</c>). null = unset.</summary>
    public char? SeqoutType { get; init; }

    /// <summary><c>-awidth</c>: residues per line; 0 = 50.</summary>
    public int Width { get; init; }

    /// <summary><c>-aglobal</c>: "show the full sequence in the alignment" (needle and stretcher default Y, water N). Changes the <c>match</c> and <c>sam</c> output.</summary>
    public bool Global { get; init; }

    /// <summary><c>-aaccshow</c>: <c>(accession)</c> after the name in the header.</summary>
    public bool ShowAccession { get; init; }

    /// <summary><c>-adesshow</c>: the description after the name in the header.</summary>
    public bool ShowDescription { get; init; }

    /// <summary><c>-ausashow</c>: the USA instead of the name.</summary>
    public bool ShowUsa { get; init; }

    /// <summary>The matrix deciding the <c>Similarity</c> count and the similarity marks (null: <see cref="SubstitutionMatrix.EdnaFull"/> for nucleotide rows, <see cref="SubstitutionMatrix.Blosum62"/> otherwise).</summary>
    public SubstitutionMatrix? Matrix { get; init; }

    /// <summary>The matrix name printed (<c>Matrix: EBLOSUM62</c>); null = the EMBOSS file name of <see cref="Matrix"/> (<c>E</c> + name), empty = no line.</summary>
    public string? MatrixName { get; init; }

    /// <summary>The program works with integer scores (stretcher: <c>ajAlignSetGapI</c> / <c>ajAlignSetScoreI</c>): gap penalties and the score are printed as integers (<c>Score: 140</c>) instead of with at least one decimal (<c>Score: 146.0</c>).</summary>
    public bool IntegerScores { get; init; }

    /// <summary>The program stores the alignment statistics itself before writing (stretcher: <c>ajAlignSetStats(align, -1, length, identities, -1, -1)</c>): only the <c>tcoffee</c> format sees the difference (its percent identity is 0 when no program set it: needle, water).</summary>
    public bool StatisticsSetByProgram { get; init; }

    /// <summary>The gap-open penalty printed (<c>Gap_penalty</c>, SAM-less formats); null = no line.</summary>
    public double? GapPenalty { get; init; }

    /// <summary>The gap-extension penalty printed (<c>Extend_penalty</c>); null = no line.</summary>
    public double? ExtendPenalty { get; init; }
}

/// <summary>An input sequence of an EMBOSS alignment program (name = first word of the FASTA header).</summary>
/// <param name="Name">The sequence name.</param>
/// <param name="Residues">The residues, case as read.</param>
/// <param name="Description">The rest of the FASTA header.</param>
/// <param name="Qualifiers">The <c>-sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein</c> of this input (A5-7); null = none.</param>
public sealed record EmbossSequence(string Name, string Residues, string Description = "", EmbossSequenceQualifiers? Qualifiers = null);

/// <summary>The output of a multi-alignment EMBOSS program (<c>supermatcher</c>, <c>needleall</c>, <c>wordmatch</c> ...): the alignments in the order the program writes them and the text of its <c>-errfile</c>.</summary>
/// <param name="Alignments">The pairwise alignments, one <c>ajAlignWrite</c> each.</param>
/// <param name="ErrorText">What the program prints to its error file (one line per failed pair), empty when nothing failed.</param>
/// <param name="MatrixName">The EMBOSS name of the scoring matrix the program used (<c>EDNAFULL</c>, <c>EBLOSUM62</c>, ...; the <c>Matrix:</c> line of its report); empty when the program did not record it.</param>
public sealed record EmbossProgramResult(IReadOnlyList<EmbossPairAlignment> Alignments, string ErrorText, string MatrixName = "");
