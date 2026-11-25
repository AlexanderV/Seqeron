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
    [McpServerTool(Name = "est2genome", Title = "Alignment — EMBOSS est2genome (spliced EST / genomic alignment)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 est2genome (Richard Mott's est_genome, embest.c): the spliced alignment of EST sequences (a FASTA text, every record is aligned) to an unspliced genomic sequence (the first record of a FASTA text). A local alignment with introns: matches +match, mismatches -mismatch, a base deleted from either sequence -gapPenalty, an intron of any length -intronPenalty, or -splicePenalty when it starts after a donor and ends on an acceptor site (usesplice, default true: GT/AG for a forward gene, CT/AC for a reversed gene). mode 'both' (default) aligns the forward and the reverse-complemented EST with forward-gene splice sites, then repeats the better one with reversed-gene splice sites and reports the best (the Note line tells which: forward est / reversed est, forward / REVERSED GENE); 'forward' and 'reverse' (case-insensitive, unique abbreviations) align one strand. An alignment is written when its score is at least minScore: the Exon lines (score, percent identity, genomic and EST coordinates, names, EST description), the +Intron / -Intron / ?Intron lines, the Span line and, after a blank line, the Segment lines; with align the alignment itself (introns as their first and last five bases and their length, > forward / < reversed splice direction) in blocks of width columns and 'Alignment Score'. reverse prints the coordinates of a reversed EST on its original strand; best = false (-nobest) also writes the other alignment (with a single strand mode the program crashes: refused with InvalidOperationException); space is the megabytes of the matrix above which the linear-space divide and conquer is used (raised to 1.5e-6 x genome length with a warning in messages); shuffle > 0 sets minScore to the best score + 1 (the program's shuffle does not change the sequence); seed has no effect; width of 0 or less with align never ends in the program (refused). Negative penalties are used as their absolute values. sequenceQualifiers: EMBOSS sequence qualifiers after the files (-sbegin1 5 -send1 90 -sreverse1 for the ESTs, -sbegin2 ... for the genomic sequence). Byte for byte equal to est2genome on 330 runs (forward / reverse genes, ESTs on either strand, every option, the linear-space recursion, the crashes).")]
    public static Est2genomeToolResult Est2genome(
        [Description("The EST sequences as FASTA text (-estsequence); every record is aligned.")] string estFasta,
        [Description("The genomic sequence as FASTA text (-genomesequence); the first record is used.")] string genomeFasta,
        [Description("Score for matching two bases (-match, default 1).")] int match = 1,
        [Description("Cost for mismatching two bases (-mismatch, default 1).")] int mismatch = 1,
        [Description("Cost for deleting a single base in either sequence, excluding introns (-gappenalty, default 2).")] int gapPenalty = 2,
        [Description("Cost for an intron, independent of length (-intronpenalty, default 40).")] int intronPenalty = 40,
        [Description("Cost for an intron between a donor and an acceptor site (-splicepenalty, default 20).")] int splicePenalty = 20,
        [Description("Exclude alignments with scores below this threshold (-minscore, default 30).")] int minScore = 30,
        [Description("Report the coordinates of a reversed EST on the original strand (-reverse, default false).")] bool reverse = false,
        [Description("Use donor and acceptor splice sites (-usesplice, default true).")] bool useSplice = true,
        [Description("Comparison mode (-mode): 'both' (default), 'forward' or 'reverse'.")] string mode = "both",
        [Description("Print only the best alignment (-best, default true; false = -nobest).")] bool best = true,
        [Description("Space threshold in megabytes for the linear-space recursion (-space, default 10.0).")] float space = 10.0f,
        [Description("Shuffle count (-shuffle, default 0).")] int shuffle = 0,
        [Description("Random number seed (-seed, default 20825; no effect on the output).")] int seed = 20825,
        [Description("Show the alignment (-align, default false).")] bool align = false,
        [Description("Alignment width (-width, default 50).")] int width = 50,
        [Description("EMBOSS sequence qualifiers as typed after the files: -sbegin1 5 -send1 90 -sreverse1 -slower1 -supper1 (input 1 = the ESTs) and -sbegin2 ... (input 2 = the genomic sequence). " + EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null)
    {
        if (string.IsNullOrWhiteSpace(estFasta))
            throw new ArgumentException("The EST sequences cannot be empty", nameof(estFasta));
        if (string.IsNullOrWhiteSpace(genomeFasta))
            throw new ArgumentException("The genomic sequence cannot be empty", nameof(genomeFasta));
        var o = EmbossEst2genomePrograms.Run(new EmbossEst2genomeRequest(estFasta, genomeFasta, match, mismatch, gapPenalty, intronPenalty, splicePenalty, minScore, reverse, useSplice, mode, best, space, shuffle, seed, align, width, sequenceQualifiers));
        return new Est2genomeToolResult(o.Text, o.Messages);
    }
}
