using System.Text.Json.Nodes;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The coordinate formats of Biopython's <c>Bio.Align</c> that describe pairwise alignments by their block coordinates (<c>AlignmentCoordinateFormats</c>).</summary>
public enum CoordinateAlignmentFormat
{
    /// <summary>UCSC BED (3 to 12 columns, <c>bed</c>).</summary>
    Bed,
    /// <summary>UCSC PSL / PSLX (<c>psl</c>).</summary>
    Psl,
    /// <summary>UCSC chain (<c>chain</c>).</summary>
    Chain,
    /// <summary>Multiple Alignment Format (<c>maf</c>), any number of sequences.</summary>
    Maf,
    /// <summary>Mauve XMFA (<c>mauve</c>), any number of sequences.</summary>
    Mauve,
    /// <summary>Sequence Alignment/Map (<c>sam</c>), one line per alignment.</summary>
    Sam,
    /// <summary>Exonerate vulgar / cigar lines (<c>exonerate</c>).</summary>
    Exonerate,
    /// <summary>UCSC bigBed (<c>bigbed</c>, binary).</summary>
    BigBed,
    /// <summary>UCSC bigMaf (<c>bigmaf</c>, binary).</summary>
    BigMaf,
    /// <summary>UCSC bigPsl (<c>bigpsl</c>, binary).</summary>
    BigPsl,
    /// <summary>GCG MSF (<c>msf</c>, reader of row alignments).</summary>
    Msf,
    /// <summary>EMBOSS pairwise / markx output (<c>emboss</c>, reader).</summary>
    Emboss,
    /// <summary>HHsearch / HHblits result (<c>hhr</c>, reader).</summary>
    Hhr,
    /// <summary>BLAST tabular output (<c>tabular</c>, reader).</summary>
    Tabular,
}

/// <summary>
/// One sequence of a <see cref="BioCoordAlignment"/> (a Biopython <c>SeqRecord</c>): its identifier, its length and the part of its residues that is known (a complete sequence has one part at 0, a sequence
/// whose contents are unknown, <c>Seq(None, length=n)</c>, has none; a PSLX file defines the residues of the aligned blocks only).
/// </summary>
public sealed class BioCoordSequence
{
    /// <summary>The identifier (<c>record.id</c>); null = <c>None</c>.</summary>
    public string? Id { get; init; }

    /// <summary>The length of the sequence; a negative value = the sequence object is <c>None</c> (a SAM target without a header).</summary>
    public long Length { get; init; }

    /// <summary>Known residues: start position -> text; null = the contents are undefined.</summary>
    public SortedDictionary<long, string>? Parts { get; init; }

    /// <summary>The description (<c>record.description</c>).</summary>
    public string Description { get; init; } = "";

    /// <summary>The annotations of the record (<c>record.annotations</c>: MAF <c>quality</c>, <c>leftStatus</c>, <c>leftCount</c>, <c>rightStatus</c>, <c>rightCount</c>).</summary>
    public JsonObject Annotations { get; init; } = [];

    /// <summary>The letter annotations of the record (<c>record.letter_annotations</c>: SAM <c>phred_quality</c>).</summary>
    public JsonObject LetterAnnotations { get; init; } = [];

    /// <summary>The residues <c>[start, end)</c>, or null when any of them is undefined (Biopython's <c>UndefinedSequenceError</c>).</summary>
    public string? Slice(long start, long end)
    {
        if (Parts is null)
            return null;
        if (end <= start)
            return "";
        foreach (var (s, text) in Parts)
        {
            if (s <= start && end <= s + text.Length)
                return text.Substring((int)(start - s), (int)(end - start));
        }

        return null;
    }
}

/// <summary>A pairwise alignment given by its block coordinates (Biopython's <c>Alignment</c>: <c>sequences</c>, <c>coordinates</c> and attributes such as <c>score</c>).</summary>
public sealed class BioCoordAlignment
{
    /// <summary>The sequences: target first, query second.</summary>
    public IReadOnlyList<BioCoordSequence> Sequences { get; init; } = [];

    /// <summary>The coordinates, one array per sequence (<c>alignment.coordinates</c>): a position decreasing along the alignment is a sequence on the reverse strand.</summary>
    public long[][] Coordinates { get; init; } = [];

    /// <summary>True for an alignment without coordinates (<c>alignment.coordinates is None</c>, an unmapped SAM line): <see cref="Coordinates"/> is empty and the first sequence is a placeholder.</summary>
    public bool NoCoordinates { get; init; }

    /// <summary>The attributes of the alignment object: <c>score</c>, <c>matches</c>, <c>misMatches</c>, <c>repMatches</c>, <c>nCount</c>, <c>thickStart</c>, <c>thickEnd</c>, <c>itemRgb</c>, <c>annotations</c> (<c>id</c> of a chain).</summary>
    public JsonObject Attributes { get; init; } = [];
}

/// <summary>Options of <c>AlignmentCoordinateFormats.Write</c>.</summary>
public sealed record CoordinateWriteOptions
{
    /// <summary>BED: the number of columns, 3 to 12 (default 12).</summary>
    public int BedN { get; init; } = 12;

    /// <summary>PSL: write the five header lines (default true).</summary>
    public bool Header { get; init; } = true;

    /// <summary>PSL: <c>null</c>, <c>lower</c> or <c>upper</c>: repeat-masked residues are counted in <c>repMatches</c>.</summary>
    public string? Mask { get; init; }

    /// <summary>PSL: the wildcard residue counted in <c>nCount</c> (default <c>N</c>; null = none).</summary>
    public string? Wildcard { get; init; } = "N";

    /// <summary>MAF: the metadata of the file header (<c>alignments.metadata</c>: <c>MAF Version</c>, <c>Scoring</c>, <c>Program</c>, <c>Comments</c>, <c>name</c>, <c>description</c>, <c>frames</c>, <c>mafDot</c>, <c>visibility</c>, <c>speciesOrder</c>); null = <c>{"MAF Version": "1"}</c>.</summary>
    public JsonObject? Metadata { get; init; }

    /// <summary>Mauve: the identifiers of the sequences (the file names, or the entry numbers of a combined file as strings); sequences are numbered by their index in this list.</summary>
    public IReadOnlyList<string>? Identifiers { get; init; }

    /// <summary>SAM: calculate the <c>MD</c> tag (<c>md=True</c>).</summary>
    public bool Md { get; init; }

    /// <summary>SAM: the target records of the header (<c>alignments.targets</c>).</summary>
    public IReadOnlyList<BioCoordSequence>? Targets { get; init; }

    /// <summary>Exonerate: <c>vulgar</c> (default) or <c>cigar</c>.</summary>
    public string ExonerateFormat { get; init; } = "vulgar";
}

/// <summary>The result of <c>AlignmentCoordinateFormats.ParseFile</c>: the alignments and the file-level attributes of the iterator (<c>alignments.metadata</c>, <c>alignments.targets</c>, ...).</summary>
public sealed class CoordinateAlignmentFile
{
    /// <summary>The alignments in file order.</summary>
    public IReadOnlyList<BioCoordAlignment> Alignments { get; init; } = [];

    /// <summary>The iterator's <c>metadata</c> dictionary.</summary>
    public JsonObject Metadata { get; init; } = [];

    /// <summary>The iterator's <c>targets</c> (SAM header @SQ records); empty when the format has none.</summary>
    public IReadOnlyList<BioCoordSequence> Targets { get; init; } = [];
}
