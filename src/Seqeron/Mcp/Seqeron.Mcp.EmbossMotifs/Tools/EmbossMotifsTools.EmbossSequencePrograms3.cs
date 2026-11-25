// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossMotifs, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossMotifs.Tools;

public partial class EmbossMotifsTools
{
    [McpServerTool(Name = "marscan", Title = "Genome — EMBOSS marscan (matrix/scaffold recognition signatures)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 marscan: the MAR/SAR recognition signatures: the 16 bp site awwrtaannwwgnnnc (one mismatch allowed) or its reverse gnnncwwnnttaywwt and the 8 bp site aataayaa or ttrttatt within 200 bp of each other; of a cluster of signatures the pair with the least distance is reported as the report out.txt (default format gff: the feature misc_signal with the positions of the two sites; every report format of EMBOSS, the report qualifiers of reportQualifiers). Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.")]
    public static EmbossProgramToolResult Marscan(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The report format (-rformat, default gff): simple, embl, genbank, gff, gff3, pir, swiss, dasgff, debug, listfile, dbmotif, diffseq, draw, restrict, excel, feattable, motif, nametable, regions, seqtable, srs, summary, table, tagseq.")] string? reportFormat = null,
        [Description("Report qualifiers as typed on the command line: -raccshow -rdesshow -rusashow -norscoreshow -norstrandshow -rmaxall N -rmaxseq N (-rformat is the parameter reportFormat).")] string? reportQualifiers = null,
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
        if (!string.IsNullOrEmpty(reportFormat))
            args.AddRange(["-rformat", reportFormat]);
        AddQualifierTokens(args, reportQualifiers);
        args.AddRange(["-outfile", "out.txt"]);
        return RunEmbossMultiProgram("marscan", args, all, nameof(sequence), "out.txt");
    }
}
