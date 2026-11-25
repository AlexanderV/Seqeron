using System.ComponentModel;
using System.Text.Json.Nodes;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-27d: the binary Biopython Bio.Align formats bigbed / bigmaf / bigpsl (bytes as base64) and the T-Coffee score_ascii evaluation; they delegate to AlignmentCoordinateFormats / SequenceAligner.TCoffeeEvaluate.
public partial class AlignmentTools
{
    private static CoordinateAlignmentFormat BigFormat(string format)
    {
        if (!global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.TryParseFormat(format, out var fmt) || fmt is not (CoordinateAlignmentFormat.BigBed or CoordinateAlignmentFormat.BigMaf or CoordinateAlignmentFormat.BigPsl))
            throw new ArgumentException($"Unknown binary format '{format}'. Supported: bigbed, bigmaf, bigpsl.", nameof(format));
        return fmt;
    }

    private static BioCoordSequence TargetFromJson(JsonNode t) => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FromJson(new JsonObject { ["sequences"] = new JsonArray(t.DeepClone()), ["coordinates"] = new JsonArray(), ["attrs"] = new JsonObject() }).Sequences[0];

    [McpServerTool(Name = "write_big_alignment", Title = "Alignment - Write alignments as bigBed / bigMaf / bigPsl (Bio.Align, base64)", ReadOnly = true)]
    [Description("Writes alignments given by their block coordinates (see write_coordinate_alignment for the JSON object) as a UCSC binary file, byte for byte as Biopython 1.88 Bio.Align.write does, and returns the bytes as base64: bigbed (BED3 - BED12 plus custom autoSQL fields, B+ tree of the chromosomes, R-tree, zoom levels, optional extra string indices), bigmaf (the MAF block of every alignment as an lstring field of the reference chromosome; the first row is named reference.chromosome) and bigpsl (BED12 + the 13 PSL fields: match counts from the residues, other-chromosome fields of the query, oCDS from the CDS feature of the query, oSequence, seqType from the molecule_type annotation; every alignment needs attrs.annotations). The data blocks are zlib-compressed (the compressed bytes depend on the zlib build: with compress = true the .NET deflate stream is used, which Biopython's zlib bytes need not equal; the structure, offsets and decompressed blocks are equal). VERIFICATION: with the compressed blocks of Biopython's own zlib fed back, 216 bigbed, 231 bigmaf and 328 bigpsl committed files equal byte for byte (scripts/alignment/generate_biopython_coord_oracle.py).")]
    public static BigAlignmentResult WriteBigAlignment(
        [Description("The alignments, each a JSON object (target first; for bigmaf all rows of the MAF block).")] string[] alignmentsJson,
        [Description("The chromosomes in the order of the alignments, as a JSON array of {\"id\": \"chr1\", \"length\": 1000} objects (bigmaf: the first id is reference.chromosome).")] string targetsJson,
        [Description("Format: bigbed, bigmaf or bigpsl.")] string format = "bigbed",
        [Description("bigbed: number of BED columns, 3 to 12.")] int bedN = 12,
        [Description("Compress the data blocks with zlib.")] bool compress = true,
        [Description("bigbed: number of items per data block.")] int itemsPerSlot = 512,
        [Description("bigbed: number of items per R-tree / B+ tree node.")] int blockSize = 256,
        [Description("bigbed: names of extra string columns to index.")] string[]? extraIndex = null,
        [Description("bigbed: a custom autoSQL declaration as JSON {\"name\", \"comment\", \"fields\": [[type, name, comment], ...]} (the BED columns first); the extra fields are read from attrs.annotations.")] string? declarationJson = null,
        [Description("bigpsl: write the CDS feature of the query in NCBI style.")] bool cds = false,
        [Description("bigpsl: include the query sequence.")] bool fa = false,
        [Description("bigpsl: null, 'lower' or 'upper': repeat-masked target residues are counted as repMatch.")] string? mask = null,
        [Description("bigpsl: the code of the wildcard letter counted in nCount (78 = N); null counts nothing, as Biopython's default string 'N' never equals a byte.")] int? wildcard = null)
    {
        ArgumentNullException.ThrowIfNull(alignmentsJson);
        ArgumentException.ThrowIfNullOrWhiteSpace(targetsJson);
        var fmt = BigFormat(format);
        var aligns = alignmentsJson.Select(j => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FromJson(JsonNode.Parse(j) ?? throw new ArgumentException("An alignment is not JSON."))).ToList();
        var targets = (JsonNode.Parse(targetsJson) as JsonArray ?? throw new ArgumentException("targetsJson must be a JSON array.", nameof(targetsJson))).Select(t => TargetFromJson(t!)).ToList();
        IReadOnlyList<BigBedField>? fields = null;
        string declName = "bed", declComment = "Browser Extensible Data";
        if (!string.IsNullOrWhiteSpace(declarationJson))
        {
            var d = JsonNode.Parse(declarationJson) as JsonObject ?? throw new ArgumentException("declarationJson must be a JSON object.", nameof(declarationJson));
            declName = d["name"]!.GetValue<string>();
            declComment = d["comment"]!.GetValue<string>();
            fields = [.. d["fields"]!.AsArray().Select(f => new BigBedField(f![0]!.GetValue<string>(), f[1]!.GetValue<string>(), f[2]!.GetValue<string>()))];
        }

        var options = new BigBedWriteOptions
        {
            BedN = bedN,
            Compress = compress,
            ItemsPerSlot = itemsPerSlot,
            BlockSize = blockSize,
            ExtraIndex = extraIndex ?? [],
            DeclarationFields = fields,
            DeclarationName = declName,
            DeclarationComment = declComment,
            Cds = cds,
            Fa = fa,
            Mask = mask,
            Wildcard = wildcard,
        };
        byte[] bytes = fmt switch
        {
            CoordinateAlignmentFormat.BigMaf => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.WriteBigMaf(aligns, targets, options),
            CoordinateAlignmentFormat.BigPsl => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.WriteBigPsl(aligns, targets, options),
            _ => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.WriteBigBed(aligns, targets, options),
        };
        return new BigAlignmentResult(Convert.ToBase64String(bytes), global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FormatName(fmt), bytes.Length, aligns.Count);
    }

    [McpServerTool(Name = "parse_big_alignment", Title = "Alignment - Read bigBed / bigMaf / bigPsl alignments (Bio.Align, base64)", ReadOnly = true)]
    [Description("Reads a UCSC binary alignment file (bigBed, bigMaf or bigPsl; base64) as Bio.Align.parse does and returns the alignments as JSON objects {\"sequences\", \"coordinates\", \"attrs\"} (see write_coordinate_alignment), the chromosomes (targets) and the autoSQL declaration {name, comment, fields}. bigbed: BED3 - BED12 columns plus the custom fields as attrs.annotations; bigmaf: the alignments of the MAF blocks (rows of the block, score / pass / empty annotations, quality / status annotations of the records; the first chromosome is named reference.chromosome); bigpsl: score, thickStart, thickEnd, itemRgb, matches, misMatches, repMatches, nCount, the query record with its sequence, molecule_type and CDS feature. A file Biopython refuses (not bigBed, truncated, inconsistent counts or coordinates) is refused with the same condition. VERIFICATION: equal to Biopython 1.88 on 229 bigbed, 963 bigmaf and 518 bigpsl committed files (damaged blocks and fields included).")]
    public static BigAlignmentParseResult ParseBigAlignment(
        [Description("The bytes of the file as base64.")] string base64,
        [Description("Format: bigbed, bigmaf or bigpsl.")] string format = "bigbed")
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(base64);
        var fmt = BigFormat(format);
        byte[] bytes;
        try
        {
            bytes = Convert.FromBase64String(base64);
        }
        catch (FormatException e)
        {
            throw new ArgumentException("base64 is not valid base64.", nameof(base64), e);
        }

        var file = fmt switch
        {
            CoordinateAlignmentFormat.BigMaf => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ParseBigMaf(bytes),
            CoordinateAlignmentFormat.BigPsl => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ParseBigPsl(bytes),
            _ => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ParseBigBed(bytes),
        };
        var arr = new JsonArray(file.Alignments.Select(global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ToJson).ToArray());
        var targets = new JsonArray(file.Targets.Select(t => global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ToJson(new BioCoordAlignment { Sequences = [t], Coordinates = [] })["sequences"]![0]!.DeepClone()).ToArray());
        var d = file.Declaration!;
        var decl = new JsonObject { ["name"] = d.Name, ["comment"] = d.Comment, ["fields"] = new JsonArray(d.Fields.Select(f => (JsonNode)new JsonArray(f.AsType, f.Name, f.Comment)).ToArray()) };
        return new BigAlignmentParseResult(arr.ToJsonString(), targets.ToJsonString(), decl.ToJsonString(), global::Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FormatName(fmt), file.Alignments.Count);
    }
}
