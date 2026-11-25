// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using Seqeron.Genomics.Kalign;
using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.MultipleAlignment.Tools;

using Seqeron.Mcp.Alignment.Tools;
using static Seqeron.Mcp.Alignment.Tools.AlignmentTools;

// B09 A2-20: Kalign 3.3.5 (port of the program, verified against the installed binary); delegates to KalignAligner.KalignRun / KalignProgram.Run.
public partial class MultipleAlignmentTools
{
    [McpServerTool(Name = "multiple_align_kalign", Title = "MSA — Kalign 3.3.5 (port of the program)", ReadOnly = true)]
    [Description("Kalign 3.3.5's multiple alignment (Lassmann 2019, Bioinformatics 36:1928; kalign -i in.fa -o out.fa), ported step by step from the program's source: the alphabet (DNA / protein) is detected from the letter counts (non-letters are gaps and are removed, sequences without residues are dropped, letter case is kept); distances = the bit-parallel edit distance of Myers (1999) of the shorter sequence (its first 255 residues) against the best window of the longer one, on the 13-letter reduced protein alphabet or ACGTN; guide tree = UPGMA (average linkage + 0.001) for fewer than 100 sequences, otherwise bisecting 2-means on the distances to 32 anchors; progressive alignment of profiles with the affine three-state dynamic program in Hirschberg's linear space (position-specific gap penalties, terminal gap extension), matrices RNA (the default for nucleotides: gap open 217, extension 39.4, terminal 292.6), dna (match 5 / mismatch -4, 8 / 6 / 0), CorBLOSUM66_13plus (default for proteins: 5.5 / 2 / 1) or GON250 (type divergent: 55 / 8 / 4). Quirks of 3.3.5 that the port reproduces and the oracle confirms: the terminal gap extension option is applied only together with a gap extension option (tgpe alone is ignored; gpe alone sets it to -1); the gap-open penalty is not charged for columns inserted into a profile; with threads = 1 the serial recursion is used (no shorter-first swap), with more threads the OpenMP recursion. VERIFICATION: rows equal the installed Debian kalign 1:3.3.5-1 (AVX2 build) on 526 / 526 committed runs (DNA / RNA / protein, every type and gap-penalty option, 1 / 2 / 4 threads, sequences of more than 500 columns, 100 - 1000 sequences, refusals). The bit-parallel distance of the portable (non-AVX2) build limits the pattern to 1 024 instead of 255 residues; sequences longer than 255 residues can therefore differ between builds. totalScore is the natural affine sum of pairs (BLOSUM62 10 / 0.2 for proteins, EDNAFULL 10 / 0.5 for nucleotides; 0 / NaN for thousands of sequences).")]
    public static KalignResultDto MultipleAlignKalign(
        [Description("Sequences to align (at least two with residues): letters (the alphabet of the detected type), any other character is a gap and is removed.")] string[] sequences,
        [Description("Kalign --type: 'auto' (default: RNA matrix for nucleotides, CorBLOSUM66_13plus for proteins), 'rna', 'dna', 'protein', 'divergent' (GON250, proteins only) or 'dna-internal' (library only: dna with terminal gap extension 8; the program's command line cannot select it). A type that does not fit the detected alphabet is refused like the program does.")] string type = "auto",
        [Description("Kalign --gpo: gap open penalty (positive = cost); omit or negative = the default of the type.")] double? gapOpen = null,
        [Description("Kalign --gpe: gap extension penalty; omit or negative = the default of the type.")] double? gapExtend = null,
        [Description("Kalign --tgpe: terminal gap extension penalty; applied only when gapExtend is also given (a quirk of 3.3.5: the program tests --gpe before applying --tgpe).")] double? terminalGapExtend = null,
        [Description("Kalign -n: number of threads (at least 1, default 4). Only 1 versus more than 1 matters: the one-thread program uses the serial recursion, which differs in the order of sequence pairs.")] int threads = 4,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        ArgumentNullException.ThrowIfNull(sequences);
        var options = new global::Seqeron.Genomics.Kalign.KalignOptions
        {
            Type = ParseKalignType(type),
            GapOpen = (float)(gapOpen ?? -1.0),
            GapExtend = (float)(gapExtend ?? -1.0),
            TerminalGapExtend = (float)(terminalGapExtend ?? -1.0),
            Threads = threads,
            CancellationToken = cancellationToken,
        };
        var run = global::Seqeron.Genomics.Kalign.KalignAligner.KalignRun(sequences, options);
        var a = run.Alignment;
        return new KalignResultDto(a.AlignedSequences, a.Consensus, a.TotalScore, a.ExactTotalScore, run.KeptIndices, run.Alphabet, run.GuideTree.Select(t => new[] { t.A, t.B, t.C }).ToArray());
    }

    private static global::Seqeron.Genomics.Kalign.KalignType ParseKalignType(string? type) => (type ?? "auto").Trim().ToLowerInvariant() switch
    {
        "auto" or "" => global::Seqeron.Genomics.Kalign.KalignType.Undefined,
        "rna" => global::Seqeron.Genomics.Kalign.KalignType.Rna,
        "dna" => global::Seqeron.Genomics.Kalign.KalignType.Dna,
        "dna-internal" => global::Seqeron.Genomics.Kalign.KalignType.DnaInternal,
        "protein" => global::Seqeron.Genomics.Kalign.KalignType.Protein,
        "divergent" => global::Seqeron.Genomics.Kalign.KalignType.ProteinDivergent,
        _ => throw new ArgumentException($"type must be 'auto', 'rna', 'dna', 'dna-internal', 'protein' or 'divergent', got '{type}'.", nameof(type)),
    };

    [McpServerTool(Name = "kalign_command", Title = "Kalign 3.3.5 — command line (type, gap penalties, threads, input files)", ReadOnly = true)]
    [Description("Kalign 3.3.5 as a function of its command line and the text of the files it names (KalignProgram.Run): the arguments are those of the program (-i / --input / --infile / --in, -o / --output / --outfile / --out, -f / --format fasta|msf|clu, --type rna|dna|internal|protein|divergent, --gpo, --gpe, --tgpe, -n / --nthreads, -q, --set, -h, -v / -V, --showw; getopt_long_only forms: one or two dashes, --name=value, unambiguous abbreviations; non-option arguments are input files, several are merged). The type word is matched like the program's strstr chain (rna, dna, internal, protein, divergent in this order: 'internal' contains 'rna' and therefore means RNA). The input files are read as FASTA, MSF or Clustal (detected from the first 100 lines like the program: MSF and Clustal rows are re-aligned after their gaps are removed); a standard input given as text is read before the named files. The alignment is written as FASTA with 60 residues per line (default), as MSF (GCG checksums, a date line taken from the clock; names longer than 256 characters are cut) or as Clustal (60-column blocks, 'Kalign (3.3.5) multiple sequence alignment' header) to the -o file or to standard output; without -q standard output starts with the banner and the program's log lines (detected alphabet, sequences read, CPU time of each stage, the warnings about aligned / unaligned input, '-n 1' alignment steps). Exit code 1 = the program stops with an error (unknown / ambiguous option, bad type word or format, fewer than one thread, a missing file, input of two alphabets, fewer than two sequences, a type that does not fit the sequences); the program segfaults on a single sequence or an MSF without Name: lines, the port reports exit code 1; letters outside the alphabet of the detected type leave uninitialised codes in the program and are refused. VERIFICATION: stdout / output file and exit status equal the installed Debian kalign 1:3.3.5-1 on 83 / 83 committed command lines (scripts/alignment/generate_kalign_command_oracle.py) and 251 / 251 more (formats, MSF / Clustal inputs, several inputs, non-quiet output under a time shim; scripts/alignment/generate_kalign_formats_oracle.py).")]
    public static KalignCommandResultDto KalignCommand(
        [Description("The command-line arguments of kalign, one per element (for example [\"-i\", \"in.fa\", \"-o\", \"out.fa\", \"--type\", \"dna\", \"-q\"]).")] string[] arguments,
        [Description("The text of the files the arguments name: file name -> FASTA / MSF / Clustal text.")] Dictionary<string, string>? files = null,
        [Description("Text of standard input (the program reads a standard input that is not a terminal before the named files; null = a terminal).")] string? stdin = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (arguments == null)
            throw new ArgumentException("arguments must not be null.", nameof(arguments));
        var result = global::Seqeron.Genomics.Kalign.KalignProgram.Run(arguments, files, stdin);
        return new KalignCommandResultDto(result.ExitCode, result.Stdout, new Dictionary<string, string>(result.Files), [.. result.Messages]);
    }
}
