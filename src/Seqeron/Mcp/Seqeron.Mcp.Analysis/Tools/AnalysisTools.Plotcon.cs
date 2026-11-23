using System.ComponentModel;
using ModelContextProtocol.Server;
using global::Seqeron.Genomics.Analysis;

namespace Seqeron.Mcp.Analysis.Tools;

// B09 WP67 (F308, R14-13): EMBOSS plotcon as an MCP tool; it builds the command line and runs the ported program (EmbossProgramRunner).
public partial class AnalysisTools
{
    [McpServerTool(Name = "plotcon", Title = "Alignment — EMBOSS plotcon (similarity plot of a multiple alignment)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 plotcon: plots the similarity of a multiple alignment along its length (the data of the plot: x = the column, y = the similarity; -graph data). Every column scores the sum over all pairs of sequences of the integer matrix score of the two residues (scorefile: EBLOSUM62 for protein, EDNAFULL for nucleic sets, or a matrix file in files; a residue that the matrix does not contain, such as a gap, scores nothing) times the weight of the sequence; the score of the window of winsize columns that ends at column k is the sum of the column scores of every row over the window, put at the middle of the window, divided by numseq (numseq - 1) winsize^2 in single precision as the program computes it (a set of one sequence or winsize 0 gives NaN: the data file prints -000nan as the program does); the first winsize / 2 columns and the last ones have no window and stay 0. sequences is the text of the alignment (any input format of seqret) or the sequences; shorter sequences are padded with gaps. The result text is the data file plotcon1.dat (the other files and messages in files / messages); graph qualifiers (-gtitle, -gxtitle, -goutfile, -gdirectory ...) are accepted as typed on the command line, the title and the y label are set by the program. VERIFICATION: standard output, error stream, exit status and written files equal to the installed plotcon 6.6.0 on 400 recorded runs plus a fresh check of 1 500 runs (every input format, ranges, matrices, window sizes including 0 and negative, the graph devices that write text, the messages of the library).")]
    public static EmbossProgramToolResult Plotcon(
        [Description("The text of the alignment (any input format of seqret, e.g. aligned FASTA, Clustal, MSF), or sequences as text; at least one sequence.")] string sequences,
        [Description("Number of columns to average the similarity over (-winsize, default 4).")] int? winsize = null,
        [Description("The scoring matrix (-scorefile): an EMBOSS matrix name (EBLOSUM62 default for protein, EDNAFULL default for nucleic acids, EBLOSUM30 ... EBLOSUM90, EPAM10 ... EPAM500, EDNAMAT, EDNASIMPLE, ENUC.4.4) or the name of a matrix file given in files.")] string? scorefile = null,
        [Description("The input format (-sformat); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin -send -sreverse, -slower -supper, -snucleotide -sprotein, -sid name, -sdbname db, -squick, -iquery query, -ioffset bytes.")] string? qualifiers = null,
        [Description("Graph qualifiers as typed on the command line: -gtitle T -gsubtitle T -gxtitle T -gytitle T -gdesc T -goutfile NAME -gdirectory DIR.")] string? graphQualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a matrix file, a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequences))
            throw new ArgumentException("The sequences cannot be empty", nameof(sequences));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequences.dat"] = sequences;
        args.AddRange(["-sequences", SequenceSpec(null, sequences, "sequences.dat")]);
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (winsize.HasValue)
            args.AddRange(["-winsize", winsize.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(scorefile))
            args.AddRange(["-scorefile", scorefile]);
        args.AddRange(["-graph", "data"]);
        AddQualifierTokens(args, graphQualifiers);
        return RunEmbossMultiProgram("plotcon", args, all, nameof(sequences), "plotcon1.dat");
    }
}
