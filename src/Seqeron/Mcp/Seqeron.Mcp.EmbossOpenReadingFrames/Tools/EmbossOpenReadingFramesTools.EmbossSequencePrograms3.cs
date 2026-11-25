// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossOpenReadingFrames, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossOpenReadingFrames.Tools;

public partial class EmbossOpenReadingFramesTools
{
    [McpServerTool(Name = "tcode", Title = "Genome — EMBOSS tcode (Fickett TESTCODE coding statistic)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 tcode: the Fickett TESTCODE statistic of the windows of window bases (at least 200) every step bases, from the base compositions and the position asymmetries of the three codon positions with the probability tables of the data file Etcode.dat (or the datafile given); every window is reported with its score and the estimation Non-coding (below 0.74), No opinion or Coding (0.95 and above) as the report out.txt in the report format reportFormat (default table: every report format of EMBOSS, the report qualifiers of reportQualifiers); with plot the values are plotted against the window mid positions (graph: data writes the data file tcode1.dat, bedgraph the track file tcode.bedgraph, none draws nothing; the devices that draw with plplot (ps, png, svg ...) are refused). Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.")]
    public static EmbossProgramToolResult Tcode(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The length of the sliding window (-window, default 200: at least 200 and at most the length of the first sequence).")] int? window = null,
        [Description("The step of the window (-step, default 3).")] int? step = null,
        [Description("The text of a testcode data file in the format of Etcode.dat (-datafile; default the file of EMBOSS).")] string? dataFile = null,
        [Description("Plot the statistic (-plot).")] bool plot = false,
        [Description("The graph device of -plot (-graph): data, bedgraph or none.")] string? graph = null,
        [Description("The report format (-rformat, default table): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq.")] string? reportFormat = null,
        [Description("Report qualifiers as typed on the command line: -raccshow -rdesshow -rusashow -norscoreshow -norstrandshow -rmaxall N -rmaxseq N (-rformat is the parameter reportFormat).")] string? reportQualifiers = null,
        [Description("Graph qualifiers as typed on the command line: -gtitle T -gsubtitle T -gxtitle T -gytitle T -gdesc T -goutfile NAME -gdirectory DIR.")] string? graphQualifiers = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.AddRange(["-sequence", SequenceSpec(null, sequence, "sequence.dat")]);
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (window.HasValue)
            args.AddRange(["-window", window.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (step.HasValue)
            args.AddRange(["-step", step.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(dataFile))
        {
            all["etcode.dat"] = dataFile;
            args.AddRange(["-datafile", "etcode.dat"]);
        }
        if (plot)
            args.Add("-plot");
        if (!string.IsNullOrEmpty(graph))
            args.AddRange(["-graph", graph]);
        if (!string.IsNullOrEmpty(reportFormat))
            args.AddRange(["-rformat", reportFormat]);
        AddQualifierTokens(args, reportQualifiers);
        AddQualifierTokens(args, graphQualifiers);
        if (plot && string.IsNullOrEmpty(graph))
            args.AddRange(["-graph", "data"]);
        args.AddRange(["-outfile", "out.txt"]);
        return RunEmbossMultiProgram("tcode", args, all, nameof(sequence), "out.txt");
    }
    [McpServerTool(Name = "checktrans", Title = "Genome — EMBOSS checktrans (STOP codons and ORFs of a protein)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 checktrans: for every protein sequence the stretches between its stops (*) that are at least orfml residues long: the report out.txt lists their positions and the number of stops, the sequences name_N are returned in the file seqout.txt (sequence format outputFormat) and the protein feature table in the file features.txt (feature format featureFormat, default gff3); with addlast (default) a stop is forced at the end (it is not counted). Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.")]
    public static EmbossProgramToolResult Checktrans(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The minimum length of an ORF to report (-orfml, default 100).")] int? orfml = null,
        [Description("Force the sequence to end with a stop (-addlast, default true).")] bool addlast = true,
        [Description("The sequence output format (-osformat, default fasta).")] string? outputFormat = null,
        [Description("The feature output format (-offormat, default gff3): gff, gff2, gff3, bed, pir, swiss, dasgff, draw, debug, refseqp.")] string? featureFormat = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.AddRange(["-sequence", SequenceSpec(null, sequence, "sequence.dat")]);
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (orfml.HasValue)
            args.AddRange(["-orfml", orfml.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!addlast)
            args.Add("-noaddlast");
        if (!string.IsNullOrEmpty(outputFormat))
            args.AddRange(["-osformat", outputFormat]);
        if (!string.IsNullOrEmpty(featureFormat))
            args.AddRange(["-offormat", featureFormat]);
        args.AddRange(["-outfile", "out.txt", "-outseq", "seqout.txt", "-outfeat", "features.txt"]);
        return RunEmbossMultiProgram("checktrans", args, all, nameof(sequence), "out.txt");
    }
    [McpServerTool(Name = "plotorf", Title = "Genome — EMBOSS plotorf (potential ORFs of a nucleotide sequence)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 plotorf: for each of the six frames the stretches from a start codon (start, default ATG) to the next stop codon of the frame (stop, default TAA,TAG,TGA; a stretch that reaches the end of the sequence ends there) as the filled rectangles of six plots; graph: data writes the six data files plotorf1.dat ... plotorf6.dat (returned in files), bedgraph one track file, none draws nothing (the devices that draw with plplot are refused). Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.")]
    public static EmbossProgramToolResult Plotorf(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The start codons, a comma list (-start, default ATG).")] string? start = null,
        [Description("The stop codons, a comma list (-stop, default TAA,TAG,TGA).")] string? stop = null,
        [Description("The graph device (-graph): data (default), bedgraph or none.")] string graph = "data",
        [Description("Graph qualifiers as typed on the command line: -gtitle T -gsubtitle T -gxtitle T -gytitle T -gdesc T -goutfile NAME -gdirectory DIR.")] string? graphQualifiers = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.AddRange(["-sequence", SequenceSpec(null, sequence, "sequence.dat")]);
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (!string.IsNullOrEmpty(start))
            args.AddRange(["-start", start]);
        if (!string.IsNullOrEmpty(stop))
            args.AddRange(["-stop", stop]);
        if (!string.IsNullOrEmpty(graph))
            args.AddRange(["-graph", graph]);
        AddQualifierTokens(args, graphQualifiers);
        return RunEmbossMultiProgram("plotorf", args, all, nameof(sequence), "stdout");
    }
}
