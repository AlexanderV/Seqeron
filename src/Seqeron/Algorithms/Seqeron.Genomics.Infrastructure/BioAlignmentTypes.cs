using System.Text.Json.Nodes;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The alignment file formats of Biopython's <c>Bio.Align</c> handled by <c>AlignmentFormats</c> (the multiple-alignment text formats).</summary>
public enum AlignmentFileFormat
{
    /// <summary>Aligned FASTA (<c>fasta</c>).</summary>
    Fasta,
    /// <summary>A2M: upper case = match, lower case = insertion, dots and dashes = gaps (<c>a2m</c>).</summary>
    A2m,
    /// <summary>ClustalW / Clustal Omega alignment (<c>clustal</c>), 50 columns per block.</summary>
    Clustal,
    /// <summary>PHYLIP: 10-character names, one line per sequence (<c>phylip</c>).</summary>
    Phylip,
    /// <summary>NEXUS data block (<c>nexus</c>).</summary>
    Nexus,
    /// <summary>Stockholm 1.0 with <c>#=GF</c> / <c>#=GS</c> / <c>#=GR</c> / <c>#=GC</c> annotation (<c>stockholm</c>).</summary>
    Stockholm,
}

/// <summary>
/// One sequence of a <see cref="BioAlignment"/> (the Biopython <c>SeqRecord</c> of an alignment row): identifier, description and the PRINTED row (the aligned part of the sequence with <c>-</c> for a gap).
/// </summary>
public sealed class BioAlignmentRecord
{
    /// <summary>The identifier (<c>record.id</c>).</summary>
    public string Id { get; init; } = "";

    /// <summary>The description (<c>record.description</c>).</summary>
    public string Description { get; init; } = "";

    /// <summary>The aligned row, <c>-</c> for a gap.</summary>
    public string Row { get; init; } = "";

    /// <summary>The sequence annotations (<c>record.annotations</c>: <c>molecule_type</c> for NEXUS, <c>organism</c> / <c>accession</c> ... for Stockholm <c>#=GS</c>), in order.</summary>
    public JsonObject Annotations { get; init; } = [];

    /// <summary>The letter annotations (<c>record.letter_annotations</c>: Stockholm <c>#=GR</c>, one letter per residue of the row), in order.</summary>
    public JsonObject LetterAnnotations { get; init; } = [];

    /// <summary>Database cross-references (<c>record.dbxrefs</c>: Stockholm <c>#=GS ... DR</c>).</summary>
    public IReadOnlyList<string> DbXrefs { get; init; } = [];
}

/// <summary>
/// A multiple alignment as <c>Bio.Align.Alignment</c> holds it, reduced to what the text formats carry: the records with their printed rows, the alignment annotations (<c>alignment.annotations</c>: Stockholm <c>#=GF</c>) and
/// the column annotations (<c>alignment.column_annotations</c>: A2M <c>state</c>, Clustal <c>clustal_consensus</c>, Stockholm <c>#=GC</c>), the latter two as JSON objects so that their values keep Biopython's types.
/// </summary>
public sealed class BioAlignment
{
    /// <summary>The sequences, one per row.</summary>
    public IReadOnlyList<BioAlignmentRecord> Records { get; init; } = [];

    /// <summary><c>alignment.annotations</c> (Stockholm: <c>identifier</c>, <c>accession</c>, <c>definition</c>, <c>author</c> (list), <c>comment</c>, <c>references</c> ...).</summary>
    public JsonObject Annotations { get; init; } = [];

    /// <summary><c>alignment.column_annotations</c>: one character per column.</summary>
    public JsonObject ColumnAnnotations { get; init; } = [];

    /// <summary>The printed rows.</summary>
    public IReadOnlyList<string> Rows => [.. Records.Select(r => r.Row)];
}

/// <summary>Options of <c>AlignmentFormats.Write</c>.</summary>
public sealed record AlignmentWriteOptions
{
    /// <summary>The program of the Clustal header (<c>Biopython</c>).</summary>
    public string ClustalProgram { get; init; } = "Biopython";

    /// <summary>The version in the Clustal header (the file reports Biopython's version: <c>1.88</c> for the reference).</summary>
    public string ClustalVersion { get; init; } = "1.88";

    /// <summary>NEXUS: <c>true</c> interleaved (70 columns per block), <c>false</c> one line per taxon, <c>null</c> interleaved above 1000 columns (Biopython's default).</summary>
    public bool? NexusInterleave { get; init; }
}
