namespace Seqeron.Mcp.Analysis.Models;

/// <summary>Result of <c>nucmer</c>: the text of the <c>.delta</c> file (or of the <c>.cluster</c> file with <c>delta = false</c>) and what the programs write to standard error (and, with <c>coords</c>, the text of the coords file); <c>Log</c> is the whole standard error of the script, <c>Files</c> the files it leaves behind by name (<c>&lt;prefix&gt;.delta</c> ...).</summary>
public record NucmerToolResult(string? Delta, string? Cluster, string Messages, string? Coords = null, string Log = "", Dictionary<string, string>? Files = null);

/// <summary>The text a MUMmer show-* / delta-filter program writes to standard output, byte for byte.</summary>
public record MummerTextResult(string Text);

/// <summary>Result of <c>show_coords</c>, <c>show_snps</c>, <c>show_aligns</c> and <c>show_diff</c>: standard output and standard error.</summary>
public record MummerProgramResult(string Text, string Messages);

/// <summary>Result of <c>show_tiling</c>: standard output and the files of -t / -u / -p.</summary>
public record ShowTilingToolResult(string Text, string? ContigList, string? Unusable, string? PseudoMolecule);

/// <summary>Result of <c>dnadiff</c>: the files of the program (the lines of the report that perl's hash order shuffles are in a fixed order).</summary>
public record DnadiffToolResult(string Report, string Delta, string Delta1, string DeltaM, string Coords1, string CoordsM, string Snps, string ReferenceDiff, string QueryDiff, string? UnalignedReference, string? UnalignedQuery);

/// <summary>Result of <c>est2genome</c>: the <c>-outfile</c> text and what the program writes to standard error.</summary>
public record Est2genomeToolResult(string Text, string Messages);

/// <summary>Result of <c>promer</c>: the text of the <c>.delta</c> file (or of the <c>.cluster</c> file with <c>delta = false</c>), the optional coords file, what the programs write to standard error and the temporary files of the script (<c>.aaref</c>, <c>.aaqry</c>, <c>.mgaps</c>), which the script deletes.</summary>
public record PromerToolResult(string? Delta, string? Cluster, string Messages, string? Coords, string AminoReference, string AminoQuery, string Mgaps);

/// <summary>Result of <c>mummer_annotate</c>: standard output and the file <c>witherrors.gaps</c> the program writes.</summary>
public record MummerAnnotateToolResult(string Text, string WithErrors);

/// <summary>Result of <c>run_mummer1</c>: the files of the script (<c>.out</c>, <c>.gaps</c>, <c>.align</c>, <c>.errorsgaps</c>).</summary>
public record RunMummer1ToolResult(string Out, string Gaps, string Align, string ErrorsGaps);

/// <summary>Result of <c>combine_mums</c>: standard output, the errors file (<c>-W</c>), standard error and the files of the program by name (the errors file and the <c>.cover</c> files).</summary>
public record CombineMumsToolResult(string Text, string ErrorsFile, string Messages, Dictionary<string, string> Files);

/// <summary>Result of <c>run_mummer3</c>: the files of the script (<c>.out</c>, <c>.gaps</c>, <c>.align</c>, <c>.errorsgaps</c>), its standard output and its standard error.</summary>
public record RunMummer3ToolResult(string Out, string Gaps, string Align, string ErrorsGaps, string Stdout, string Messages);

/// <summary>Result of <c>repeat_match</c>: standard output, standard error and the exit status of the program.</summary>
public record RepeatMatchToolResult(string Text, string Messages, int ExitCode);

/// <summary>Result of <c>exact_tandems</c>: standard output (the table of tandem repeats, or the two progress lines of the installed script), standard error, the exit status and the files left behind.</summary>
public record ExactTandemsToolResult(string Text, string Messages, int ExitCode, string[] FilesLeft);

/// <summary>Result of <c>mapview</c>: the FIG files (or, with pdf / ps, the temporary FIG files the script feeds to fig2dev), the fig2dev command lines and output names, standard output, standard error, the exit status and the side file of an mgaps input.</summary>
public record MapviewToolResult(Dictionary<string, string> Figs, string[] Commands, string[] Outputs, string Stdout, string Messages, int ExitCode, Dictionary<string, string> Side);

/// <summary>Result of <c>mummerplot</c>: the files of the script (<c>.fplot</c>, <c>.rplot</c>, <c>.hplot</c>, <c>.gp</c>, <c>.filter</c>), standard output, standard error, the exit status and the gnuplot command line.</summary>
public record MummerplotToolResult(Dictionary<string, string> Files, string Stdout, string Messages, int ExitCode, string? GnuplotCommand);
