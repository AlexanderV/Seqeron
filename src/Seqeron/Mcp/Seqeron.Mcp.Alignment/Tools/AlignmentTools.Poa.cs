using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-26: POA 2.0 (port of the program, verified against the installed binary); delegates to SequenceAligner.PoaRun / PoaProgram.Run.
public partial class AlignmentTools
{
    [McpServerTool(Name = "multiple_align_poa", Title = "MSA — POA 2.0 partial-order alignment (port of the program)", ReadOnly = true)]
    [Description("POA 2.0's partial-order multiple alignment (Lee, Grasso & Sharlow 2002, Bioinformatics 18:452; poa -read_fasta in.fa -clustal out.aln matrix), ported step by step from the program's source: every sequence is aligned to the growing partial-order graph by dynamic programming over the graph (each node scored against each residue, predecessors of both graphs, the program's affine gap model with a truncation / decay length and a free initial gap in local mode, tie-breaking by the smallest end cell) and fused into it (aligned identical letters share a node, other aligned letters form an alignment ring, the ring is kept contiguous in the graph's linear order). Default is iterative local alignment in input order; `global` aligns end to end, `progressive` merges the pairs in descending order of their alignment scores (or of the `pairScores` text of name1 name2 score lines), `fuseAll` fuses identical letters on rings, `heaviestBundle` adds the heaviest-bundle consensus sequences CONSENS0, CONSENS1, ... (a sequence joins a bundle when at least `heaviestBundleMinimum` of its residues lie on the path), `preserveSequenceOrder` restores the input order after progressive merging. The scoring matrix is a POA matrix file as text (symbol line, integer rows, GAP-PENALTIES=open ext1 ext2, GAP-PENALTIES-X, GAP-TRUNCATION-LENGTH, GAP-DECAY-LENGTH); omitted = BLOSUM80 with gap penalties 12 / 6 / 6 (the blosum80.mat of the distribution). Letters that are not symbols of the matrix become its first symbol (so lower-case letters need `case` = 'upper' for the protein matrix). Returns the rows, the program's PO file (the graph), the FASTA-PIR text and the CLUSTAL text. VERIFICATION: rows and all three files equal the installed Debian poa 2.0+20060928 on 359 / 359 committed command lines (protein BLOSUM80, nucleotide matrices with truncation / decay gap models, every option; scripts/alignment/generate_poa_oracle.py), see poa_command. Not ported: MSA / PO file input (-read_msa, -read_msa2, -read_msa_list, -subset, -remove).")]
    public static PoaResultDto MultipleAlignPoa(
        [Description("Sequences to align (at least one non-empty): residues; white space is removed.")] string[] sequences,
        [Description("Sequence names (the first word of a FASTA header, at most 31 characters are kept); omit for S1, S2, ...")] string[]? names = null,
        [Description("POA matrix file as text; omit for the BLOSUM80 matrix of the distribution (gap penalties 12 / 6 / 6).")] string? matrix = null,
        [Description("poa -do_global: global alignment (default: local alignment with free end gaps).")] bool global = false,
        [Description("poa -do_progressive: merge sequences / partial orders in descending order of their pairwise alignment scores (default: iterative in input order).")] bool progressive = false,
        [Description("poa -read_pairscores: the text of a pair-score file (lines 'name1 name2 score'); the merge order follows these scores.")] string? pairScores = null,
        [Description("poa -fuse_all: fuse identical letters on align rings.")] bool fuseAll = false,
        [Description("poa -hb: heaviest bundling: add the consensus sequences CONSENS0, CONSENS1, ... to the alignment.")] bool heaviestBundle = false,
        [Description("poa -hbmin: minimum fraction of a sequence's residues on a bundle's path for the sequence to join it (default 0.9).")] double heaviestBundleMinimum = 0.9,
        [Description("poa -best: restrict the PIR text to the sequences of the heaviest bundle (bundle 0; every alignment column is still written).")] bool best = false,
        [Description("poa -preserve_seqorder: write the sequences in their input order (progressive merging reorders them).")] bool preserveSequenceOrder = false,
        [Description("Letter case of the input: 'keep' (default), 'lower' (poa -tolower) or 'upper' (poa -toupper).")] string @case = "keep",
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        ArgumentNullException.ThrowIfNull(sequences);
        var seqs = new List<global::Seqeron.Genomics.Alignment.PoaSequence>();
        for (int i = 0; i < sequences.Length; i++)
        {
            if (sequences[i] == null)
                throw new ArgumentException("A sequence is null.", nameof(sequences));
            seqs.Add(new global::Seqeron.Genomics.Alignment.PoaSequence(names != null && i < names.Length ? names[i] : "S" + (i + 1), null, sequences[i]));
        }

        var options = new global::Seqeron.Genomics.Alignment.PoaOptions
        {
            Global = global,
            Progressive = progressive,
            PairScores = pairScores,
            FuseAll = fuseAll,
            HeaviestBundle = heaviestBundle,
            HeaviestBundleMinimum = (float)heaviestBundleMinimum,
            Best = best,
            PreserveSequenceOrder = preserveSequenceOrder,
            Case = ParsePoaCase(@case),
            CancellationToken = cancellationToken,
        };
        var m = matrix == null ? null : global::Seqeron.Genomics.Alignment.PoaScoreMatrix.Parse(matrix);
        var run = global::Seqeron.Genomics.Alignment.SequenceAligner.PoaRun(seqs, m, options);
        return new PoaResultDto(run.Names, run.Rows, run.Order, run.Po, run.Pir, run.Clustal, run.Stdout);
    }

    private static global::Seqeron.Genomics.Alignment.PoaCase ParsePoaCase(string? value) => (value ?? "keep").Trim().ToLowerInvariant() switch
    {
        "keep" or "" => global::Seqeron.Genomics.Alignment.PoaCase.Keep,
        "lower" => global::Seqeron.Genomics.Alignment.PoaCase.Lower,
        "upper" => global::Seqeron.Genomics.Alignment.PoaCase.Upper,
        _ => throw new ArgumentException($"case must be 'keep', 'lower' or 'upper', got '{value}'.", nameof(value)),
    };

    [McpServerTool(Name = "poa_command", Title = "POA 2.0 — command line (FASTA / matrix / pair-score files)", ReadOnly = true)]
    [Description("POA 2.0 as a function of its command line and the text of the files it names (PoaProgram.Run): poa [-read_fasta FILE] [-read_msa FILE] [-read_msa2 FILE] [-read_msa_list FILE] [-subset FILE] [-subset2 FILE] [-remove FILE] [-remove2 FILE] [-do_global] [-do_progressive] [-read_pairscores FILE] [-fuse_all] [-hb] [-hbmin X] [-best] [-preserve_seqorder] [-tolower] [-toupper] [-pir FILE] [-clustal FILE] [-po FILE] [-printmatrix LETTERS] [-v] MATRIXFILE, with the program's argument loop (every argument is tested against every flag; unknown flags are ignored; the first argument that does not start with a dash is the matrix file). Output: the PIR / CLUSTAL / PO files named by -pir / -clustal / -po, standard output (-v gap penalty tables, -printmatrix, the permutation line of -preserve_seqorder). Exit code 1 = the program stops with an error (no matrix or unreadable matrix, no sequence file, no sequences, conflicting flags), 255 = called without arguments (usage). The alignment progress lines of standard error are not reproduced (the program prints uninitialised %ld arguments there); the '...Read N sequences...' / '...Wrote N sequences...' lines are. The MSA / PO inputs (-read_msa, -read_msa2, -read_msa_list; the format is detected from the first line: a PO file, a FASTA-PIR alignment or a CLUSTAL alignment; rows are read in 510-character pieces like the program's fgets) and the filters -subset / -subset2 / -remove / -remove2 (SOURCENAME= lines) are ported; an alignment row name beyond 31 characters, a PO file whose links / sources / rings are inconsistent and a row without residues stop the run (exit 1 / 255) where the program aborts, crashes or loops. VERIFICATION: exit code, standard output, the Read / Wrote lines and the three output files equal the installed Debian poa 2.0+20060928 on 359 / 359 FASTA command lines and 687 / 687 MSA / PO input command lines (scripts/alignment/generate_poa_msa_oracle.py).")]
    public static PoaCommandResultDto PoaCommand(
        [Description("The command-line arguments of poa, one per element (for example [\"-read_fasta\", \"in.fa\", \"-clustal\", \"out.aln\", \"blosum80.mat\"]).")] string[] arguments,
        [Description("The text of the files the arguments name: file name -> text (FASTA / MSA / PO file, subset list, matrix file, pair-score file).")] Dictionary<string, string>? files = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (arguments == null)
            throw new ArgumentException("arguments must not be null.", nameof(arguments));
        var result = global::Seqeron.Genomics.Alignment.PoaProgram.Run(arguments, files);
        return new PoaCommandResultDto(result.ExitCode, result.Stdout, result.Stderr, new Dictionary<string, string>(result.Files));
    }
}
