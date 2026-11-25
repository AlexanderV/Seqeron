// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.MultipleAlignment, licensed under the terms stated in LICENSE (this directory).
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Alignment.Tools;
using static Seqeron.Mcp.Alignment.Tools.AlignmentTools;

namespace Seqeron.Mcp.MultipleAlignment.Tools;

public partial class MultipleAlignmentTools
{
    [McpServerTool(Name = "multiple_align_probcons", Title = "MSA — ProbCons 1.12 (port of the program)", ReadOnly = true)]
    [Description("ProbCons 1.12's multiple alignment of protein sequences (Do, Mahabhashyam, Brudno & Batzoglou 2005, Genome Res 15:330; probcons in.fa), ported step by step from the program's source: pair-HMM posterior match probabilities of every pair (forward / backward in single precision with the program's log-sum approximation; a match state and two pairs of insertion states; emission tables of BLOCKS 5.0 / BLOSUM62, transitions trained on BAliBASE 2.0; upper and lower case alike, other characters at the table default), `consistency` rounds of the consistency transformation P'(x_i-y_j) = 1/|S| sum_z sum_k P(x_i-z_k) P(z_k-y_j) on the sparse matrices (cells >= 0.01), an expected-accuracy guide tree (UPGMA-like: merge the largest value, average weighted by the merged value), progressive maximum-expected-accuracy alignment of profiles, and `iterativeRefinement` rounds of randomized-partitioning refinement (the bipartitions come from glibc rand() at its default seed 1, as in the program). Optional: `preTraining` rounds of EM re-estimation of the transition probabilities on the input itself, `alignmentOrder` (rows in alignment order), the cosmetic gap penalties gapOpenPenalty / gapExtensionPenalty (<= 0; in the program they overflow a heap buffer whenever the two profiles differ in size, so only equal-size merges are verified), `cutoff` (subtracted from the posterior before the profile alignments), `parameters` (a parameter file as text, probcons -p) and `annotate` (the -annot column scores, 0 - 200). Letter case is kept. VERIFICATION: stdout, standard error (parameter set, guide tree, -v lines) and written files equal the installed Debian probcons 1.12-14 on 407 / 407 committed command lines (scripts/alignment/generate_probcons_oracle.py), see probcons_command. totalScore is the natural BLOSUM62 affine sum of pairs (the program reports no score).")]
    public static ProbConsResultDto MultipleAlignProbCons(
        [Description("Protein sequences (letters, '*'; at least one). Letters outside the 20 amino acids have the program's default emissions.")] string[] sequences,
        [Description("Sequence headers (the first word is the name in the guide tree and Clustal output); omit for S1, S2, ...")] string[]? names = null,
        [Description("probcons -c: passes of consistency transformation, 0 - 5 (default 2).")] int consistency = 2,
        [Description("probcons -ir: passes of iterative refinement, 0 - 1000 (default 100).")] int iterativeRefinement = 100,
        [Description("probcons -pre: rounds of EM pre-training, 0 - 20 (default 0).")] int preTraining = 0,
        [Description("probcons -a: report the sequences in alignment order instead of input order.")] bool alignmentOrder = false,
        [Description("probcons -go: cosmetic gap-open penalty (<= 0; 0 together with a zero extension penalty = none).")] double gapOpenPenalty = 0,
        [Description("probcons -ge: cosmetic gap-extension penalty (<= 0).")] double gapExtensionPenalty = 0,
        [Description("probcons -co: value in [0, 1] subtracted from the posterior probabilities before the profile alignments.")] double cutoff = 0,
        [Description("probcons -p: the text of a parameter file (initial distribution, gap open, gap extension, alphabet, lower triangle of the match emissions, single emissions); omit for the built-in parameters.")] string? parameters = null,
        [Description("probcons -annot: also return the per-column annotation score (0 - 200: twice the mean posterior of the aligned residue pairs of the column, in percent).")] bool annotate = false,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        ArgumentNullException.ThrowIfNull(sequences);
        var options = new global::Seqeron.Genomics.TCoffee.ProbConsOptions
        {
            ConsistencyReps = consistency,
            IterativeRefinementReps = iterativeRefinement,
            PreTrainingReps = preTraining,
            AlignmentOrder = alignmentOrder,
            GapOpenPenalty = (float)gapOpenPenalty,
            GapContinuePenalty = (float)gapExtensionPenalty,
            Cutoff = (float)cutoff,
            Parameters = parameters == null ? null : global::Seqeron.Genomics.TCoffee.ProbConsParameters.Parse(parameters),
            Names = names,
            Annotate = annotate,
            CancellationToken = cancellationToken,
        };
        var run = global::Seqeron.Genomics.TCoffee.ProbConsAligner.ProbConsRun(sequences, options);
        var a = run.Alignment;
        return new ProbConsResultDto(a.AlignedSequences, a.Consensus, a.TotalScore, a.ExactTotalScore, run.Order, run.GuideTree, run.Annotation);
    }
    [McpServerTool(Name = "probcons_command", Title = "ProbCons 1.12 — command line (options, MFA / MSF / Clustal input files)", ReadOnly = true)]
    [Description("ProbCons 1.12 as a function of its command line and the text of the files it names (ProbConsProgram.Run): -clustalw, -c / --consistency, -ir / --iterative-refinement, -pre / --pre-training, -pairs, -viterbi, -v / --verbose, -annot FILE, -t / --train FILE, -e / --emissions, -p / --paramfile FILE, -a / --alignment-order, -go / --gap-open, -ge / --gap-extension, -co / --cutoff, with the program's number reading (strtol base 0, strtod prefix), ranges and error texts; input files are FASTA (any of '.', '-' and white space inside the sequence are removed; an empty sequence ends the file), MSF or Clustal alignments (ParseMSF with its quirks: letters upper-cased, gaps kept - the program then overflows a buffer, so gapped alignment input is not verified). Output: the alignment as FASTA (60 columns) or Clustal-like text with the conservation line (stdout), the pairwise alignment files of -pairs / -viterbi (header1-header2.fasta / .aln), the -annot column scores (setw(4) per line), the trained parameter file of -t (%.10f) - and the program's standard-error text byte for byte (heading, parameter set at %.10g, 'Loading sequence file', 'Alignment tree: ...', -v progress). Exit code 1 = the program exits with an error (unknown option, bad value, missing file, unknown character, -e without -t); with no input file the Debian build prints 'Zero sequences delivered to function DoAlign()' and exits with 0. VERIFICATION: exit code, stdout, stderr and files equal the installed Debian probcons 1.12-14 on 407 / 407 committed command lines.")]
    public static ProbConsCommandResultDto ProbConsCommand(
        [Description("The command-line arguments of probcons, one per element (for example [\"-clustalw\", \"-c\", \"1\", \"in.fa\"]).")] string[] arguments,
        [Description("The text of the files the arguments name: file name -> text (sequence files and the -p parameter file).")] Dictionary<string, string>? files = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (arguments == null)
            throw new ArgumentException("arguments must not be null.", nameof(arguments));
        var result = global::Seqeron.Genomics.TCoffee.ProbConsProgram.Run(arguments, files);
        return new ProbConsCommandResultDto(result.ExitCode, result.Stdout, result.Stderr, new Dictionary<string, string>(result.Files));
    }
}
