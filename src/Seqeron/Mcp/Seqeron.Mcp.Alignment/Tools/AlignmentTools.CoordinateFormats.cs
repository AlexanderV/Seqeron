using System.ComponentModel;
using System.Text.Json.Nodes;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-27c: the Biopython Bio.Align coordinate formats (bed, psl, chain) as MCP tools; they delegate to AlignmentCoordinateFormats.
public partial class AlignmentTools
{
    [McpServerTool(Name = "write_coordinate_alignment", Title = "Alignment - Write pairwise alignments as BED / PSL / chain (Bio.Align)", ReadOnly = true)]
    [Description("Writes alignments given by their block coordinates in a Biopython 1.88 Bio.Align coordinate format, byte for byte as Bio.Align.write does: bed (3 to 12 columns, bedN), psl (header, mask lower / upper for repeat matches, wildcard N), chain, maf (any number of sequences, score / pass / empty-region annotations, quality and i lines of the records, file header metadata), mauve (xmfa: header from the metadata and identifiers, one block per alignment). An alignment is a JSON object {\"sequences\": [{\"id\": \"chr1\", \"length\": 1000, \"parts\": {\"0\": \"ACGT...\"} or null}, {query}], \"coordinates\": [[target positions], [query positions]], \"attrs\": {\"score\": 42, \"matches\": 10, ...}} (target first, query second; a position that decreases along the alignment is a sequence on the reverse strand; parts = the known residues by start position, null = unknown contents; PSL counts matches / mismatches from the residues when they are known, otherwise every aligned base is a match). VERIFICATION: written text equal to Biopython 1.88 on 620 / 620 committed cases (scripts/alignment/generate_biopython_coord_oracle.py).")]
    public static CoordinateAlignmentResult WriteCoordinateAlignment(
        [Description("The alignments, each a JSON object as described (a JSON array of objects).")] string[] alignmentsJson,
        [Description("Format: bed, psl, chain, maf or mauve.")] string format = "psl",
        [Description("BED: number of columns, 3 to 12.")] int bedN = 12,
        [Description("PSL: write the five header lines.")] bool header = true,
        [Description("PSL: null, 'lower' or 'upper': repeat-masked residues of the target are counted as repMatches.")] string? mask = null,
        [Description("PSL: the wildcard residue counted in nCount (default N; empty = none).")] string? wildcard = "N",
        [Description("Mauve: the identifiers of the sequences (the file names, or the entry numbers of a combined file as strings: '0', '1', ...); sequences are numbered by their index.")] string[]? identifiers = null,
        [Description("MAF / Mauve: the metadata of the file header as a JSON object (MAF: {\"MAF Version\": \"1\", \"Scoring\": ..., \"Program\": ..., \"Comments\": [...], \"name\", \"description\", \"mafDot\", \"visibility\", \"speciesOrder\": [...]}, default {\"MAF Version\": \"1\"}; Mauve: {\"FormatVersion\": \"Mauve1\", \"File\": \"genome.fa\" (one combined file), \"BackboneFile\": ...}).")] string? metadataJson = null)
    {
        ArgumentNullException.ThrowIfNull(alignmentsJson);
        if (!global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.TryParseFormat(format, out var fmt))
            throw new ArgumentException($"Unknown coordinate format '{format}'. Supported: bed, psl, chain, maf, mauve.", nameof(format));
        var aligns = alignmentsJson.Select(j => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FromJson(JsonNode.Parse(j) ?? throw new ArgumentException("An alignment is not JSON."))).ToList();
        string text = global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.Write(aligns, fmt, new global::Seqeron.Genomics.Infrastructure.CoordinateWriteOptions { BedN = bedN, Header = header, Mask = mask, Wildcard = string.IsNullOrEmpty(wildcard) ? null : wildcard, Metadata = string.IsNullOrWhiteSpace(metadataJson) ? null : JsonNode.Parse(metadataJson) as JsonObject, Identifiers = identifiers });
        return new CoordinateAlignmentResult(text, global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FormatName(fmt), aligns.Count);
    }

    [McpServerTool(Name = "parse_coordinate_alignment", Title = "Alignment - Read BED / PSL / chain alignments (Bio.Align)", ReadOnly = true)]
    [Description("Reads the alignments of a BED, PSL / PSLX, chain, MAF or Mauve xmfa file as Bio.Align.parse does and returns them as JSON objects {\"sequences\", \"coordinates\", \"attrs\"} (see write_coordinate_alignment): sequences with their lengths and the residues a PSLX file defines, the block coordinates (a decreasing coordinate = reverse strand), score / matches / misMatches / repMatches / nCount / thickStart / thickEnd / itemRgb / annotations.id as the format carries them. Text Biopython refuses (inconsistent columns, block counts, starts and ends) is refused with the same condition; a PSLX line of a protein aligned to translated DNA is not read. VERIFICATION: equal to Biopython 1.88 on 860 / 860 committed texts.")]
    public static CoordinateAlignmentResult ParseCoordinateAlignment(
        [Description("The text of the file.")] string text,
        [Description("Format: bed, psl, chain, maf or mauve.")] string format = "psl")
    {
        ArgumentNullException.ThrowIfNull(text);
        if (!global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.TryParseFormat(format, out var fmt))
            throw new ArgumentException($"Unknown coordinate format '{format}'. Supported: bed, psl, chain, maf, mauve.", nameof(format));
        var aligns = global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.Parse(text, fmt);
        var arr = new JsonArray(aligns.Select(global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ToJson).ToArray());
        return new CoordinateAlignmentResult(arr.ToJsonString(), global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FormatName(fmt), aligns.Count);
    }
}
