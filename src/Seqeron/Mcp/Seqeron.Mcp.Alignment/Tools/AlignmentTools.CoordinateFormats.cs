using System.ComponentModel;
using System.Text.Json.Nodes;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-27c / A2-27d: the Biopython Bio.Align coordinate formats (bed, psl, chain, maf, mauve, sam, exonerate; readers msf, emboss, hhr, tabular) as MCP tools; they delegate to AlignmentCoordinateFormats.
public partial class AlignmentTools
{
    [McpServerTool(Name = "write_coordinate_alignment", Title = "Alignment - Write alignments as BED / PSL / chain / MAF / Mauve / SAM / Exonerate (Bio.Align)", ReadOnly = true)]
    [Description("Writes alignments given by their block coordinates in a Biopython 1.88 Bio.Align coordinate format, byte for byte as Bio.Align.write does: bed (3 to 12 columns, bedN), psl (header, mask lower / upper for repeat matches, wildcard N), chain, maf (any number of sequences, score / pass / empty-region annotations, quality and i lines of the records, file header metadata), mauve (xmfa: header from the metadata and identifiers, one block per alignment), sam (header from metadata {\"HD\": {\"VN\": \"1.6\"}, \"PG\": [...]} and targets, MD tag with md = true, soft / hard clips, qualities, tags in attrs.annotations, attrs.operations), exonerate (vulgar or cigar lines with the header and footer). An alignment is a JSON object {\"sequences\": [{\"id\": \"chr1\", \"length\": 1000, \"parts\": {\"0\": \"ACGT...\"} or null}, {query}], \"coordinates\": [[target positions], [query positions]], \"attrs\": {\"score\": 42, \"matches\": 10, ...}} (target first, query second; a position that decreases along the alignment is a sequence on the reverse strand; parts = the known residues by start position, null = unknown contents; PSL counts matches / mismatches from the residues when they are known, otherwise every aligned base is a match). VERIFICATION: written text equal to Biopython 1.88 on 1 001 / 1 001 committed cases (scripts/alignment/generate_biopython_coord_oracle.py).")]
    public static CoordinateAlignmentResult WriteCoordinateAlignment(
        [Description("The alignments, each a JSON object as described (a JSON array of objects).")] string[] alignmentsJson,
        [Description("Format: bed, psl, chain, maf, mauve, sam or exonerate.")] string format = "psl",
        [Description("BED: number of columns, 3 to 12.")] int bedN = 12,
        [Description("PSL: write the five header lines.")] bool header = true,
        [Description("PSL: null, 'lower' or 'upper': repeat-masked residues of the target are counted as repMatches.")] string? mask = null,
        [Description("PSL: the wildcard residue counted in nCount (default N; empty = none).")] string? wildcard = "N",
        [Description("Mauve: the identifiers of the sequences (the file names, or the entry numbers of a combined file as strings: '0', '1', ...); sequences are numbered by their index.")] string[]? identifiers = null,
        [Description("MAF / Mauve: the metadata of the file header as a JSON object (MAF: {\"MAF Version\": \"1\", \"Scoring\": ..., \"Program\": ..., \"Comments\": [...], \"name\", \"description\", \"mafDot\", \"visibility\", \"speciesOrder\": [...]}, default {\"MAF Version\": \"1\"}; Mauve: {\"FormatVersion\": \"Mauve1\", \"File\": \"genome.fa\" (one combined file), \"BackboneFile\": ...}).")] string? metadataJson = null,
        [Description("SAM: calculate the MD tag (needs the contents of the target).")] bool md = false,
        [Description("SAM: the target records of the header (@SQ lines) as a JSON array of {\"id\": ..., \"length\": ..., \"description\": ..., \"annotations\": {...}} objects.")] string? targetsJson = null,
        [Description("Exonerate: vulgar or cigar.")] string exonerateFormat = "vulgar")
    {
        ArgumentNullException.ThrowIfNull(alignmentsJson);
        if (!global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.TryParseFormat(format, out var fmt) || fmt is not (CoordinateAlignmentFormat.Bed or CoordinateAlignmentFormat.Psl or CoordinateAlignmentFormat.Chain or CoordinateAlignmentFormat.Maf or CoordinateAlignmentFormat.Mauve or CoordinateAlignmentFormat.Sam or CoordinateAlignmentFormat.Exonerate))
            throw new ArgumentException($"Unknown coordinate format '{format}'. Supported writers: bed, psl, chain, maf, mauve, sam, exonerate.", nameof(format));
        var aligns = alignmentsJson.Select(j => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FromJson(JsonNode.Parse(j) ?? throw new ArgumentException("An alignment is not JSON."))).ToList();
        string text = global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.Write(aligns, fmt, new global::Seqeron.Genomics.Infrastructure.CoordinateWriteOptions { BedN = bedN, Header = header, Mask = mask, Wildcard = string.IsNullOrEmpty(wildcard) ? null : wildcard, Metadata = string.IsNullOrWhiteSpace(metadataJson) ? null : JsonNode.Parse(metadataJson) as JsonObject, Identifiers = identifiers, Md = md, ExonerateFormat = exonerateFormat, Targets = string.IsNullOrWhiteSpace(targetsJson) ? null : [.. (JsonNode.Parse(targetsJson) as JsonArray ?? throw new ArgumentException("targetsJson must be a JSON array.", nameof(targetsJson))).Select(t => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FromJson(new JsonObject { ["sequences"] = new JsonArray(t!.DeepClone()), ["coordinates"] = new JsonArray(), ["attrs"] = new JsonObject() }).Sequences[0])] });
        return new CoordinateAlignmentResult(text, global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FormatName(fmt), aligns.Count);
    }

    [McpServerTool(Name = "parse_coordinate_alignment", Title = "Alignment - Read BED / PSL / chain / MAF / Mauve / SAM / Exonerate / MSF / EMBOSS / HHR / tabular alignments (Bio.Align)", ReadOnly = true)]
    [Description("Reads the alignments of a BED, PSL / PSLX (nucleotides, and a protein aligned to translated DNA with its CDS feature), chain, MAF, Mauve xmfa, SAM (header metadata and targets, MD tag, hard / soft clips, qualities, tags), Exonerate vulgar / cigar, GCG MSF, EMBOSS pair / srspair / simple (needle, water, stretcher), HHsearch hhr or BLAST / FASTA tabular (-outfmt 7, -m 8CB / 8CC) file as Bio.Align.parse does and returns them as JSON objects {\"sequences\", \"coordinates\", \"attrs\"} (see write_coordinate_alignment): sequences with their lengths and the residues a PSLX file defines, the block coordinates (a decreasing coordinate = reverse strand), score / matches / misMatches / repMatches / nCount / thickStart / thickEnd / itemRgb / annotations.id as the format carries them. Text Biopython refuses (inconsistent columns, block counts, starts and ends) is refused with the same condition; the file-level values (the iterator's metadata, the SAM targets, the parser warnings) are returned in metadataJson / targetsJson / warnings. VERIFICATION: equal to Biopython 1.88 on 2 331 / 2 331 committed texts (incl. EMBOSS 6.6.0, HH-suite 3.3.0, BLAST+ and fasta36 outputs).")]
    public static CoordinateAlignmentResult ParseCoordinateAlignment(
        [Description("The text of the file.")] string text,
        [Description("Format: bed, psl, chain, maf, mauve, sam, exonerate, msf, emboss, hhr or tabular.")] string format = "psl")
    {
        ArgumentNullException.ThrowIfNull(text);
        if (!global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.TryParseFormat(format, out var fmt) || fmt is CoordinateAlignmentFormat.BigBed or CoordinateAlignmentFormat.BigMaf or CoordinateAlignmentFormat.BigPsl)
            throw new ArgumentException($"Unknown coordinate format '{format}'. Supported readers: bed, psl, chain, maf, mauve, sam, exonerate, msf, emboss, hhr, tabular.", nameof(format));
        var file = global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ParseFile(text, fmt);
        var arr = new JsonArray(file.Alignments.Select(global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ToJson).ToArray());
        var targets = new JsonArray(file.Targets.Select(t => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ToJson(new global::Seqeron.Genomics.Infrastructure.BioCoordAlignment { Sequences = [t], Coordinates = [] })["sequences"]![0]!.DeepClone()).ToArray());
        return new CoordinateAlignmentResult(arr.ToJsonString(), global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FormatName(fmt), file.Alignments.Count, file.Metadata.ToJsonString(), targets.ToJsonString(), [.. file.Warnings]);
    }
}
