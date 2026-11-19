using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-14b / A2-14c / A2-14d: ClustalW 2.1 as a command line (guide-tree files, -OUTORDER, output formats, secondary-structure gap-penalty masks); delegates to ClustalWProgram.Run.
public partial class AlignmentTools
{
    [McpServerTool(Name = "clustalw_command", Title = "ClustalW 2.1 — command line (alignment, profiles, guide trees, output formats, structure masks)", ReadOnly = true)]
    [Description("ClustalW 2.1 as a function of its command line and the text of the files it names (ClustalWProgram.Run): the arguments are those of the program (-INFILE=x.fa -TYPE=DNA|PROTEIN -ALIGN -OUTPUT=CLUSTAL|GCG|GDE|PHYLIP|PIR|NEXUS|FASTA -OUTORDER=INPUT|ALIGNED -CASE -SEQNOS -SEQNO_RANGE -RANGE=m,n -OUTFILE=name, -NEWTREE=name / -USETREE=name, -CONVERT (the input as it is in the chosen format), -PROFILE1=a -PROFILE2=b -PROFILE | -SEQUENCES with -NEWTREE1/2 / -USETREE1/2, -NOSECSTR1/2 -SECSTROUT=STRUCTURE|MASK|BOTH|NONE -HELIXGAP -STRANDGAP -LOOPGAP -TERMINALGAP -HELIXENDIN -HELIXENDOUT -STRANDENDIN -STRANDENDOUT, and every alignment parameter: -GAPOPEN -GAPEXT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE -KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -SCORE -PWMATRIX -PWGAPOPEN -PWGAPEXT -MATRIX -DNAMATRIX -PWDNAMATRIX -TRANSWEIGHT -NOPGAP -NOHGAP -NOVGAP -HGAPRESIDUES -GAPDIST -CLUSTERING=NJ|UPGMA -ITERATION=TREE|ALIGNMENT|NONE -NUMITER); option names are case-insensitive, values of enumerations are matched like the program does. `files` gives the text of every file the line names (input sequences / profiles in FASTA, CLUSTAL or GDE format; guide-tree files; matrix files); the result holds the files the program writes (alignment file in the chosen format with the default name <input base>.aln / .msf / .gde / .phy / .pir / .nxs / .fasta or -OUTFILE, guide-tree files <base>.dnd or the -NEWTREE name). Guide tree files: the neighbour-joining tree written by printPhylipTree (and the UPGMA tree), the user tree read with Tree::readTree (rooted or unrooted, with or without distances; the program's refusals - not starting with '(', wrong number of leaves, unknown name, unrooted without distances, trees that send the program into undefined behaviour - are tool errors). Row order: -OUTORDER=ALIGNED (the program's default) lists the rows in the order they entered the progressive alignment. Secondary-structure masks: !SS_ / !GM_ lines of CLUSTAL profile files and \"SS_ / \"GM_ entries of GDE files are read as the program does and change the gap penalties of -PROFILE / -SEQUENCES. VERIFICATION: guide trees and row order on 3 300 sets, the seven alignment file formats byte for byte on 2 611 sets (301 of them -CONVERT), structure masks on 1 200 profile pairs, all against the installed clustalw 2.1. Not implemented (registered in the audit): -TREE, -BOOTSTRAP, -STATS, -PIM and the input formats PIR, MSF, RSF, EMBL, PHYLIP, NEXUS.")]
    public static ClustalWCommandResultDto ClustalWCommand(
        [Description("The command-line arguments of clustalw, one per element (for example [\"-INFILE=x.fa\", \"-TYPE=DNA\", \"-ALIGN\", \"-OUTPUT=GCG\"]); a first argument without '-' or '/' is the input file.")] string[] arguments,
        [Description("The text of the files the arguments name: file name -> content (input sequences or profiles, -USETREE files, matrix files).")] Dictionary<string, string>? files = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (arguments == null || arguments.Length == 0)
            throw new ArgumentException("arguments must contain at least one command-line argument.", nameof(arguments));
        var result = global::Seqeron.Genomics.Alignment.ClustalWProgram.Run(arguments, files ?? new Dictionary<string, string>());
        return new ClustalWCommandResultDto(new Dictionary<string, string>(result.Files), [.. result.Messages]);
    }
}
