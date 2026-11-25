// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

public partial class EmbossPairwiseAlignmentTools
{
    [McpServerTool(Name = "wordfinder", Title = "Alignment — EMBOSS wordfinder (local alignments of large sequences by word matches)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 wordfinder: the target sequences (bsequence) are indexed by their words of wordlen residues, every query sequence (asequence) is looked up in the index, the diagonal with the most word residues fixes the start points and a banded Smith-Waterman alignment of width cells (gap penalties gapopen / gapextend, matrix datafile: an EMBOSS matrix name such as EDNAFULL or EBLOSUM62, or the text of a matrix file) is written in the alignment format alignFormat (every -aformat of EMBOSS, default simple, with the sub-header lines Wordscore and Alignlength) as out.txt; the matches and alignments excluded by limitmatch / lowmatch / limitalign / lowalign and the matches per target are written to err.txt. Equals the installed program on 300 runs (every input format of seqret, ranges, the report and graph options, the messages of the library; standard output, error stream, exit status and written files byte for byte). B09 WP63, F282.")]
    public static EmbossProgramToolResult Wordfinder(
        [Description("The text of the query sequence file (-asequence: any input format of seqret; one or more sequences), or a sequence as letters.")] string querySequence,
        [Description("The text of the target sequence file (-bsequence), or a sequence as letters.")] string targetSequence,
        [Description("The gap opening penalty (-gapopen, 0 to 1000, default 30).")] float? gapOpen = null,
        [Description("The gap extension penalty (-gapextend, 0 to 10, default 1.5).")] float? gapExtend = null,
        [Description("The scoring matrix (-datafile): the text of a matrix file.")] string? dataFile = null,
        [Description("The scoring matrix as an EMBOSS data file name (-datafile): EDNAFULL (nucleotides, default), EBLOSUM62 (proteins, default), EPAM250 ...")] string? matrix = null,
        [Description("The width of the alignment band (-width, at least 1, default 16).")] int? width = null,
        [Description("The word length of the initial matching (-wordlen, at least 3, default 6).")] int? wordLength = null,
        [Description("The maximum match score, 0 for no limit (-limitmatch).")] int? limitMatch = null,
        [Description("The maximum alignment length, 0 for no limit (-limitalign).")] int? limitAlign = null,
        [Description("The minimum match score (-lowmatch).")] int? lowMatch = null,
        [Description("The minimum alignment length (-lowalign).")] int? lowAlign = null,
        [Description("The alignment format (-aformat, default simple): simple, pair, fasta, msf, clustal, mega, nexus, phylip, selex, treecon, markx0 ... markx10, match, multiple, score, srs, tcoffee, unknown, sam, ...")] string? alignFormat = null,
        [Description("The residues per line of the alignment (-awidth).")] int? alignWidth = null,
        [Description("Alignment qualifiers as typed on the command line: -aaccshow -adesshow -ausashow -aglobal.")] string? alignQualifiers = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null)
    {
        if (string.IsNullOrWhiteSpace(querySequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(querySequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["asequence.dat"] = querySequence;
        args.AddRange(["-asequence", SequenceSpec(null, querySequence, "asequence.dat")]);
        all["bsequence.dat"] = targetSequence;
        args.AddRange(["-bsequence", SequenceSpec(null, targetSequence, "bsequence.dat")]);
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (gapOpen.HasValue)
            args.AddRange(["-gapopen", gapOpen.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (gapExtend.HasValue)
            args.AddRange(["-gapextend", gapExtend.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(dataFile))
        {
            all["matrix.txt"] = dataFile;
            args.AddRange(["-datafile", "matrix.txt"]);
        }
        if (!string.IsNullOrEmpty(matrix))
            args.AddRange(["-datafile", matrix]);
        if (width.HasValue)
            args.AddRange(["-width", width.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (wordLength.HasValue)
            args.AddRange(["-wordlen", wordLength.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (limitMatch.HasValue)
            args.AddRange(["-limitmatch", limitMatch.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (limitAlign.HasValue)
            args.AddRange(["-limitalign", limitAlign.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (lowMatch.HasValue)
            args.AddRange(["-lowmatch", lowMatch.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (lowAlign.HasValue)
            args.AddRange(["-lowalign", lowAlign.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(alignFormat))
            args.AddRange(["-aformat", alignFormat]);
        if (alignWidth.HasValue)
            args.AddRange(["-awidth", alignWidth.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        AddQualifierTokens(args, alignQualifiers);
        args.AddRange(["-outfile", "out.txt", "-errfile", "err.txt"]);
        return RunEmbossMultiProgram("wordfinder", args, all, nameof(querySequence), "out.txt");
    }
}
