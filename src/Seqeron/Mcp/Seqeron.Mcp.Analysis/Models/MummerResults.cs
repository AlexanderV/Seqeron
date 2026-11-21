namespace Seqeron.Mcp.Analysis.Models;

/// <summary>Result of <c>nucmer</c>: the text of the <c>.delta</c> file (or of the <c>.cluster</c> file with <c>delta = false</c>) and what the programs write to standard error (and, with <c>coords</c>, the text of the coords file).</summary>
public record NucmerToolResult(string? Delta, string? Cluster, string Messages, string? Coords = null);

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
