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
    [McpServerTool(Name = "getorf", Title = "Genome — EMBOSS getorf (open reading frames)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 getorf: the regions between STOP codons or between START and STOP codons of the six frames (the circular sequence and the nucleotides around the codons too) are returned as protein or nucleotide sequences named name_1, name_2, ... with the header [start - end] (REVERSE SENSE) and the description of the sequence; the scan is Translator.FindOrfsGetorfProgram (the model of getorf_FindORFs over the EMBOSS genetic code files). Equals the installed program on 800 runs (every input format of seqret, ranges, the output formats, the messages of the library; byte for byte). B09 WP61, F280.")]
    public static EmbossProgramToolResult Getorf(
        [Description("The text of the sequence file (any input format of seqret), or a sequence as letters.")] string sequence,
        [Description("The genetic code (-table): 0 standard with AUG start only, 1 .. 6, 9 .. 16, 21 .. 23 (the EMBOSS EGC files).")] string? table = null,
        [Description("The minimum size of an ORF in nucleotides (-minsize, default 30).")] int? minSize = null,
        [Description("The maximum size of an ORF in nucleotides (-maxsize, default 1000000).")] int? maxSize = null,
        [Description("The type of output (-find): 0 protein between STOP codons, 1 protein between START and STOP, 2 nucleotides between STOP codons, 3 nucleotides between START and STOP, 4 nucleotides around START codons, 5 around initial STOP codons, 6 around ending STOP codons.")] string? find = null,
        [Description("A START codon at the beginning of a protein is an M (-methionine, default true).")] bool methionine = true,
        [Description("The sequence is circular (-circular).")] bool circular = false,
        [Description("ORFs of the reverse complement are found as well (-reverse, default true).")] bool reverse = true,
        [Description("The number of nucleotides on each side of the codon for find 4 .. 6 (-flanking, default 100).")] int? flanking = null,
        [Description("The input format (-sformat1); default: the program finds the format.")] string? sequenceFormat = null,
        [Description("Sequence input qualifiers as typed after the sequence on the command line: -sbegin1 -send1 -sreverse1, -slower1 -supper1 -snucleotide1 -sprotein1 -scircular1, -sid1 name, -sdbname1 db, -squick1, -iquery1 query, -ioffset1 bytes.")] string? qualifiers = null,
        [Description("The sequence output format (-osformat): fasta (default), embl, genbank, gcg, nbrf, ig, raw, msf, clustal, phylip, nexus, ... every format of seqret.")] string outputFormat = "fasta",
        CancellationToken cancellationToken = default,
        [Description("More files by name (file name -> text) that the program may read (a file named by a USA, @file).")] IReadOnlyDictionary<string, string>? files = null,
        [Description("The USA of the sequence input (-sequence) as the program reads it, instead of the text: [format::]file[begin:end:r], file:id, @list, asis::SEQUENCE.")] string? sequenceValue = null)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("The sequence cannot be empty", nameof(sequence));
        cancellationToken.ThrowIfCancellationRequested();
        var all = new Dictionary<string, string>(files ?? new Dictionary<string, string>());
        var args = new List<string>();
        all["sequence.dat"] = sequence;
        args.Add("-sequence");
        args.Add(SequenceSpec(sequenceValue, sequence, "sequence.dat"));
        if (!string.IsNullOrWhiteSpace(sequenceFormat))
            args.AddRange(["-sformat1", sequenceFormat]);
        AddQualifierTokens(args, qualifiers);
        if (!string.IsNullOrEmpty(table))
            args.AddRange(["-table", table]);
        if (minSize.HasValue)
            args.AddRange(["-minsize", minSize.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (maxSize.HasValue)
            args.AddRange(["-maxsize", maxSize.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        if (!string.IsNullOrEmpty(find))
            args.AddRange(["-find", find]);
        if (!methionine)
            args.Add("-nomethionine");
        if (circular)
            args.Add("-circular");
        if (!reverse)
            args.Add("-noreverse");
        if (flanking.HasValue)
            args.AddRange(["-flanking", flanking.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)]);
        return RunEmbossSequenceProgram("getorf", args, all, outputFormat, nameof(sequence));
    }
}
