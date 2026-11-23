using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 WP67 (F307, R14-10b): Bio.Align.CodonAligner and the legacy Bio.Align.MultipleSeqAlignment of Biopython 1.88 as MCP tools.
public partial class AlignmentTools
{
    private const string MsaRecordsDescription = "The rows as a JSON array of records {\"id\": ..., \"name\": ..., \"description\": ..., \"seq\": \"ACGT-\", \"annotations\": {...}, \"letter_annotations\": {\"key\": \"string or array, one entry per letter\"}, \"dbxrefs\": [...]}; only seq is required (id / name / description default to <unknown id> / <unknown name> / <unknown description>).";

    private static SeqRecord RecordFromJson(JsonNode? node, string parameter)
    {
        var o = node as JsonObject ?? throw new ArgumentException("A record is a JSON object {\"id\", \"seq\", ...}.", parameter);
        string seq = o["seq"]?.GetValue<string>() ?? throw new ArgumentException("A record needs a seq.", parameter);
        return new SeqRecord(seq, o["id"]?.GetValue<string>(), o["name"]?.GetValue<string>(), o["description"]?.GetValue<string>(),
            o["annotations"] as JsonObject, o["letter_annotations"] as JsonObject, (o["dbxrefs"] as JsonArray)?.Select(x => x!.GetValue<string>()));
    }

    private static JsonObject RecordToJson(SeqRecord r) => new()
    {
        ["id"] = r.Id, ["name"] = r.Name, ["description"] = r.Description, ["seq"] = r.Seq, ["annotations"] = r.Annotations.DeepClone(),
        ["letter_annotations"] = r.LetterAnnotations.DeepClone(), ["dbxrefs"] = new JsonArray([.. r.DbXrefs.Select(x => (JsonNode)JsonValue.Create(x)!)]),
    };

    private static JsonObject MsaToJson(MultipleSeqAlignment m) => new()
    {
        ["records"] = new JsonArray([.. m.Select(r => (JsonNode)RecordToJson(r))]), ["annotations"] = m.Annotations.DeepClone(), ["column_annotations"] = m.ColumnAnnotations,
    };

    private static MultipleSeqAlignment MsaFromParts(string? recordsJson, string? annotationsJson, string? columnAnnotationsJson, string parameter)
    {
        var records = JsonNode.Parse(Need(recordsJson, parameter)) as JsonArray ?? throw new ArgumentException("The records are a JSON array of record objects.", parameter);
        return new MultipleSeqAlignment(records.Select(r => RecordFromJson(r, parameter)),
            string.IsNullOrWhiteSpace(annotationsJson) ? null : JsonNode.Parse(annotationsJson) as JsonObject ?? throw new ArgumentException("annotationsJson must be a JSON object.", nameof(annotationsJson)),
            string.IsNullOrWhiteSpace(columnAnnotationsJson) ? null : JsonNode.Parse(columnAnnotationsJson) as JsonObject ?? throw new ArgumentException("columnAnnotationsJson must be a JSON object.", nameof(columnAnnotationsJson)));
    }

    [McpServerTool(Name = "codon_aligner", Title = "Alignment — Biopython CodonAligner (protein to nucleotide with frame shifts)", ReadOnly = true)]
    [Description("Bio.Align.CodonAligner of Biopython 1.88: aligns a nucleotide sequence to its protein by dynamic programming over the three translated frames of the nucleotide sequence (genetic code by NCBI table number; DNA or RNA, a codon with a letter outside IUPAC is an error, ambiguous codons translate to X). Every residue is aligned to one codon; between two codons the nucleotide sequence may advance by 1, 2, 4 or 5 instead of 3 (a frame shift) at frameshiftMinusTwoScore / MinusOne / PlusOne / PlusTwo (default -3 each, or all four with frameshiftScore); leading and trailing nucleotides are free; residues and translated codons score matchScore / mismatchScore (defaults 1 / 0), 0 against the wildcard letter (default X). Returns the optimal score and every alignment of it as the vertices (target = residue positions, query = nucleotide positions; a frame shift is two vertices at the same residue position, the query going back for -2 / -1 and skipping for +1 / +2). Biopython's path iterator reads memory its traceback never wrote when the optimal alignment skips leading residues (it then lists one-vertex alignments for ever, with a heap overflow): those inputs are refused with an error, as the score() is defined for all inputs (scoreOnly = true). VERIFICATION: score, aligned score and every alignment in the order of the iterator equal to Biopython 1.88 on 1 400 of 1 500 recorded inputs (the other 100 are the refused ones; scripts/alignment/generate_codon_aligner_oracle.py).")]
    public static CodonAlignerResult CodonAlign(
        [Description("The protein sequence (amino-acid letters).")] string protein,
        [Description("The nucleotide sequence (DNA or RNA).")] string nucleotides,
        [Description("NCBI genetic code table number (default 1, the standard code).")] int tableNumber = 1,
        [Description("Score of a residue equal to the translated codon (default 1).")] double matchScore = 1,
        [Description("Score of a residue different from the translated codon (default 0).")] double mismatchScore = 0,
        [Description("A single wildcard letter that scores 0 against every letter (default X).")] string wildcard = "X",
        [Description("Score of all four frame shifts (default -3); the four below override it.")] double frameshiftScore = -3,
        [Description("Score of the shift -2 (the next codon overlaps the previous one by two nucleotides).")] double? frameshiftMinusTwoScore = null,
        [Description("Score of the shift -1.")] double? frameshiftMinusOneScore = null,
        [Description("Score of the shift +1 (one nucleotide is skipped between two codons).")] double? frameshiftPlusOneScore = null,
        [Description("Score of the shift +2.")] double? frameshiftPlusTwoScore = null,
        [Description("Return at most this many alignments (default 100).")] int maxAlignments = 100,
        [Description("true: only the score (defined for every input).")] bool scoreOnly = false)
    {
        ArgumentNullException.ThrowIfNull(protein);
        ArgumentNullException.ThrowIfNull(nucleotides);
        if (string.IsNullOrEmpty(wildcard) || wildcard.Length != 1)
            throw new ArgumentException("wildcard must be a single character.", nameof(wildcard));
        if (maxAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(maxAlignments), "maxAlignments must be at least 1.");
        var aligner = new CodonAligner(tableNumber)
        {
            MatchScore = matchScore,
            MismatchScore = mismatchScore,
            Wildcard = wildcard[0],
            FrameshiftScore = frameshiftScore,
        };
        if (frameshiftMinusTwoScore is { } a)
            aligner.FrameshiftMinusTwoScore = a;
        if (frameshiftMinusOneScore is { } b)
            aligner.FrameshiftMinusOneScore = b;
        if (frameshiftPlusOneScore is { } c)
            aligner.FrameshiftPlusOneScore = c;
        if (frameshiftPlusTwoScore is { } d)
            aligner.FrameshiftPlusTwoScore = d;
        if (scoreOnly)
            return new CodonAlignerResult(aligner.Score(protein, nucleotides), 0, []);
        var found = aligner.Align(protein, nucleotides);
        return new CodonAlignerResult(found.Score, found.Count, [.. found.Take(maxAlignments).Select(x => new CodonAlignmentDto(x.Coordinates))]);
    }

    [McpServerTool(Name = "multiple_seq_alignment", Title = "Alignment — legacy Biopython MultipleSeqAlignment (rows, columns, slices, add, sort, substitutions)", ReadOnly = true)]
    [Description("One operation of the legacy Bio.Align.MultipleSeqAlignment of Biopython 1.88 on an alignment of records (rows of one length with identifier, name, description, annotations, per-letter annotations, cross-references; alignment annotations and per-column annotations): summary (str(alignment): at most 20 rows, 50 letters), length (the number of rows), alignment_length, row (alignment[i], negative from the end), rows (alignment[a:b:c]: a sub-alignment; the column annotations are kept when all rows are), row_columns (alignment[i, a:b:c]: the part of a row), cell (alignment[i, j]), column (alignment[a:b:c, j]: the letters of a column), slice (alignment[rows, columns]: a sub-alignment with sliced column annotations), add (alignment + otherJson: columns appended pairwise, common annotations kept, per-column annotations concatenated), sort (by id, length, sequence, description or gaps; reverse), substitutions (the weighted letter-pair counts: the 'weight' annotations multiply, gaps are skipped, the matrix is symmetric with the counts shared between both sides), alignment (the modern Alignment object: gap-free sequences and block coordinates; Biopython 1.88 empties the letter annotations of the source rows in this call, so the returned rows keep theirs and the source is not returned), append / extend (a record / records of the alignment length), delete_row / delete_rows, column_annotations (the restricted dictionary: set_column_annotation adds one entry of the length the dictionary was created with). Slices are Python slices ('1:3', ':', '::-1', '-2:'). Errors are Biopython's (Sequences must all be the same length; list index out of range; index out of range; the per-column and per-letter annotation messages). Feature objects of records and file output through AlignIO are not part of the tool (use write_alignment for the text formats). VERIFICATION: results and error messages equal to Biopython 1.88 on 7 030 operations of 1 200 recorded alignments (scripts/alignment/generate_multiple_seq_alignment_oracle.py).")]
    public static MultipleSeqAlignmentResult MultipleSeqAlignmentOperation(
        [Description("The operation: summary, length, alignment_length, row, rows, row_columns, cell, column, slice, add, sort, substitutions, alignment, append, extend, delete_row, delete_rows, column_annotations, set_column_annotation.")] string operation,
        [Description(MsaRecordsDescription)] string recordsJson,
        [Description("Annotations of the whole alignment (JSON object).")] string? annotationsJson = null,
        [Description("Per-column annotations (JSON object: key -> string or array with one entry per column).")] string? columnAnnotationsJson = null,
        [Description("Row index (row, row_columns, cell, delete_row); negative counts from the end.")] int row = 0,
        [Description("Row slice 'a:b:c' (rows, column, slice, delete_rows; default ':').")] string? rows = null,
        [Description("Column slice 'a:b:c' (row_columns, slice; default ':').")] string? columns = null,
        [Description("Column index (cell, column); negative counts from the end.")] int column = 0,
        [Description("add: the second alignment as a JSON object {\"records\": [...], \"annotations\": {...}, \"column_annotations\": {...}}.")] string? otherJson = null,
        [Description("append: one record (JSON object); extend: a JSON array of records.")] string? newRecordsJson = null,
        [Description("sort key: id (default), length, sequence, description or gaps.")] string sortKey = "id",
        [Description("sort: descending.")] bool reverse = false,
        [Description("set_column_annotation: the key.")] string? annotationKey = null,
        [Description("set_column_annotation: the value (JSON: a string or an array).")] string? annotationValueJson = null)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(operation);
        string op = operation.Trim().ToLowerInvariant();
        var m = MsaFromParts(recordsJson, annotationsJson, columnAnnotationsJson, nameof(recordsJson));
        static string JStr(string v) => JsonValue.Create(v)!.ToJsonString();
        MultipleSeqAlignmentResult Alignment(MultipleSeqAlignment a) => new(op, MsaToJson(a).ToJsonString(), a.ToString());
        switch (op)
        {
            case "summary": return new(op, JStr(m.ToString()), m.ToString());
            case "length": return new(op, m.Count.ToString(CultureInfo.InvariantCulture));
            case "alignment_length": return new(op, m.GetAlignmentLength().ToString(CultureInfo.InvariantCulture));
            case "row": return new(op, RecordToJson(m.Row(row)).ToJsonString());
            case "rows": return Alignment(m.Rows(ParseSlice(rows)));
            case "row_columns": return new(op, RecordToJson(m.Row(row, ParseSlice(columns))).ToJsonString());
            case "cell": return new(op, JStr(m.Cell(row, column)));
            case "column": return new(op, JStr(m.Column(ParseSlice(rows), column)));
            case "slice": return Alignment(m.Slice(ParseSlice(rows), ParseSlice(columns)));
            case "add":
            {
                var other = JsonNode.Parse(Need(otherJson, nameof(otherJson))) as JsonObject ?? throw new ArgumentException("otherJson must be a JSON object.", nameof(otherJson));
                return Alignment(m.Add(MsaFromParts(other["records"]?.ToJsonString() ?? "[]", other["annotations"]?.ToJsonString(), other["column_annotations"]?.ToJsonString(), nameof(otherJson))));
            }

            case "sort":
            {
                Func<SeqRecord, IComparable>? key = sortKey.Trim().ToLowerInvariant() switch
                {
                    "id" => null,
                    "length" => r => r.Length,
                    "sequence" => r => r.Seq,
                    "description" => r => r.Description,
                    "gaps" => r => r.Seq.Count(c => c == '-'),
                    _ => throw new ArgumentException("sortKey: id, length, sequence, description or gaps.", nameof(sortKey)),
                };
                m.Sort(key, reverse);
                return Alignment(m);
            }

            case "substitutions":
            {
                var s = m.Substitutions();
                var values = new JsonArray();
                for (int i = 0; i < s.Letters.Length; i++)
                    values.Add((JsonNode)new JsonArray([.. Enumerable.Range(0, s.Letters.Length).Select(j => (JsonNode)JsonValue.Create(s.Values[i, j])!)]));
                return new(op, new JsonObject { ["letters"] = s.Letters, ["values"] = values }.ToJsonString());
            }

            case "alignment":
                return new(op, AlignmentToJson(m.ToCoordinateAlignment(keepSourceLetterAnnotations: true)));
            case "append":
                m.Append(RecordFromJson(JsonNode.Parse(Need(newRecordsJson, nameof(newRecordsJson))), nameof(newRecordsJson)));
                return Alignment(m);
            case "extend":
                m.Extend((JsonNode.Parse(Need(newRecordsJson, nameof(newRecordsJson))) as JsonArray ?? throw new ArgumentException("newRecordsJson must be a JSON array of records.", nameof(newRecordsJson))).Select(r => RecordFromJson(r, nameof(newRecordsJson))));
                return Alignment(m);
            case "delete_row":
                m.RemoveRow(row);
                return Alignment(m);
            case "delete_rows":
                m.RemoveRows(ParseSlice(rows));
                return Alignment(m);
            case "column_annotations":
                return new(op, m.ColumnAnnotations.ToJsonString());
            case "set_column_annotation":
                m.SetColumnAnnotation(Need(annotationKey, nameof(annotationKey)), JsonNode.Parse(Need(annotationValueJson, nameof(annotationValueJson))));
                return Alignment(m);
            default:
                throw new ArgumentException($"Unknown operation '{operation}'.", nameof(operation));
        }
    }
}
