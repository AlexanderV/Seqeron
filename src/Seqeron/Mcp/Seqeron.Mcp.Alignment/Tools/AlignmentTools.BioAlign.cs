using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 R14-10 (F303): Bio.Align of Biopython 1.88 as MCP tools: the PairwiseAligner (global / local / fogsaa, wildcard, strand), the counts of an alignment and the operations of the Alignment object.
public partial class AlignmentTools
{
    private const string AlignmentJsonDescription = "An alignment as JSON: {\"sequences\": [\"ACGT...\", ...] (or objects {\"id\": ..., \"length\": ..., \"parts\": {\"0\": \"ACGT...\"}}), \"coordinates\": [[target positions], [query positions], ...]} (one coordinate row per sequence; a position that decreases along the alignment is a sequence on the reverse strand), as parse_coordinate_alignment returns it.";

    private static BioCoordAlignment ParseAlignmentJson(string json, string parameter)
    {
        if (string.IsNullOrWhiteSpace(json))
            throw new ArgumentException("The alignment JSON is empty.", parameter);
        var root = JsonNode.Parse(json) ?? throw new ArgumentException("The alignment is not JSON.", parameter);
        if (root["sequences"] is JsonArray seqs)
        {
            for (int i = 0; i < seqs.Count; i++)
            {
                if (seqs[i] is JsonValue v && v.TryGetValue<string>(out var s))
                    seqs[i] = new JsonObject { ["id"] = null, ["length"] = (long)s.Length, ["parts"] = new JsonObject { ["0"] = s } };
            }
        }

        return AlignmentCoordinateFormats.FromJson(root);
    }

    private static List<BioCoordAlignment> ParseMany(string[]? json, string name) => [.. (json ?? throw new ArgumentException($"{name} is required for this operation.", name)).Select(j => ParseAlignmentJson(j, name))];

    private static string Need(string? value, string name) => string.IsNullOrWhiteSpace(value) ? throw new ArgumentException($"{name} is required for this operation.", name) : value;

    private static string AlignmentToJson(BioCoordAlignment a) => AlignmentCoordinateFormats.ToJson(a).ToJsonString();

    private static PySlice ParseSlice(string? text)
    {
        if (string.IsNullOrWhiteSpace(text))
            return PySlice.All;
        var parts = text.Split(':');
        if (parts.Length is < 2 or > 3)
            throw new ArgumentException($"'{text}' is not a slice 'start:stop[:step]'.", nameof(text));
        int? Number(string t) => string.IsNullOrWhiteSpace(t) ? null : int.Parse(t.Trim(), CultureInfo.InvariantCulture);
        return new PySlice(Number(parts[0]), Number(parts[1]), parts.Length == 3 ? Number(parts[2]) : null);
    }

    private static int[] ParseIntList(string text) => [.. (JsonNode.Parse(text) as JsonArray ?? throw new ArgumentException("A list of integers is needed, e.g. [1, 2, -1].", nameof(text))).Select(x => x!.GetValue<int>())];

    private static PairwiseAligner BuildAligner(string? mode, double? matchScore, double? mismatchScore, string? matrix, string? wildcard, string? scoring, string? gapScoresJson, string? gapFunctionsJson = null)
    {
        var aligner = new PairwiseAligner(string.IsNullOrWhiteSpace(scoring) ? null : scoring.Trim().ToLowerInvariant());
        if (!string.IsNullOrWhiteSpace(mode))
        {
            aligner.Mode = mode.Trim().ToLowerInvariant() switch
            {
                "global" => PairwiseAlignerMode.Global,
                "local" => PairwiseAlignerMode.Local,
                "fogsaa" => PairwiseAlignerMode.Fogsaa,
                _ => throw new ArgumentException($"Unknown mode '{mode}': global, local or fogsaa.", nameof(mode)),
            };
        }

        if (matchScore is { } ms)
            aligner.MatchScore = ms;
        if (mismatchScore is { } mm)
            aligner.MismatchScore = mm;
        if (!string.IsNullOrWhiteSpace(matrix))
            aligner.SubstitutionMatrix = ResolveMatrix(matrix);
        if (!string.IsNullOrEmpty(wildcard))
        {
            if (wildcard.Length != 1)
                throw new ArgumentException("wildcard must be a single character.", nameof(wildcard));
            aligner.Wildcard = wildcard[0];
        }

        if (!string.IsNullOrWhiteSpace(gapScoresJson))
        {
            var gaps = JsonNode.Parse(gapScoresJson) as JsonObject ?? throw new ArgumentException("gapScoresJson must be a JSON object {\"open_gap_score\": -2, \"extend_gap_score\": -1, ...}.", nameof(gapScoresJson));
            foreach (var (key, value) in gaps)
                aligner.SetGapScore(key, value!.GetValue<double>());
        }

        if (!string.IsNullOrWhiteSpace(gapFunctionsJson))
        {
            // B09 WP67 (F307): gap-score functions of (position, length) by family name: {"gap_score": {"kind": "lin", "p": [2, 0.5]}, "insertion_score": ..., "deletion_score": ...}
            var functions = JsonNode.Parse(gapFunctionsJson) as JsonObject ?? throw new ArgumentException("gapFunctionsJson must be a JSON object {\"gap_score\": {\"kind\": \"lin\", \"p\": [2, 0.5]}}.", nameof(gapFunctionsJson));
            foreach (var (key, value) in functions)
            {
                var spec = value as JsonObject ?? throw new ArgumentException($"gapFunctionsJson: '{key}' must be an object {{\"kind\": ..., \"p\": [...]}}.", nameof(gapFunctionsJson));
                string kind = spec["kind"]?.GetValue<string>() ?? throw new ArgumentException($"gapFunctionsJson: '{key}' has no kind.", nameof(gapFunctionsJson));
                double[] p = spec["p"] is JsonArray arr ? [.. arr.Select(x => x!.GetValue<double>())] : [];
                aligner.SetGapFunction(key, GapScoreFunctions.Create(kind, p));
            }
        }

        return aligner;
    }

    private const string AlignerParamNote = "gapScoresJson: a JSON object of Biopython gap-score properties (open_gap_score, extend_gap_score, gap_score, left_gap_score, right_open_gap_score, end_gap_score, target_gap_score, query_extend_gap_score, open_internal_insertion_score ... every name of the property tree; insertion = a gap in the target (sequence A), deletion = a gap in the query); the default of every gap score is -1.";

    [McpServerTool(Name = "pairwise_aligner", Title = "Pairwise — Biopython PairwiseAligner (global / local / FOGSAA, wildcard, strand)", ReadOnly = true)]
    [Description("Bio.Align.PairwiseAligner of Biopython 1.88 on two sequences: the optimal score, the algorithm Biopython selects (Needleman-Wunsch, Gotoh, Smith-Waterman, FOGSAA) and every alignment of the optimal score (target = sequenceA, query = sequenceB; block coordinates, printed rows, the Bio.Align text). mode global / local / fogsaa (the Fast Optimal Global Sequence Alignment Algorithm: one optimal alignment found by branch and bound, with Biopython's warnings when match < mismatch or a gap score exceeds the mismatch score), matchScore / mismatchScore or a substitution matrix (name or text), the wildcard letter (scores 0 against every letter, scalar scoring only), the presets scoring = blastn / megablast / blastp, the twelve gap scores by name, strand '-' (sequenceA against the reverse complement of sequenceB; the query coordinates decrease). " + AlignerParamNote + " Letters are compared case-insensitively. Gap-score functions of (position, length) (gapFunctionsJson: Biopython's Waterman-Smith-Beyer algorithm for gap_score / insertion_score / deletion_score, strands + and -; the functions are named families because a tool cannot receive code) and gap scores above zero (Needleman-Wunsch, Gotoh, Smith-Waterman and Waterman-Smith-Beyer with a positive score, local alignments included) are supported; FOGSAA ignores gap functions, as in Biopython. VERIFICATION: score, algorithm name, warnings and the alignments (as a set; FOGSAA the alignment itself) equal to Biopython 1.88 on 698 / 698 recorded configurations (scripts/alignment/generate_pairwise_aligner_oracle.py) and on 1 501 / 1 501 configurations with gap-score functions or gap scores above zero (scripts/alignment/generate_pairwise_wsb_oracle.py; including the RuntimeError the reference raises on one local alignment).")]
    public static PairwiseAlignerResult PairwiseAlignerAlign(
        [Description("The target sequence (any letters).")] string sequenceA,
        [Description("The query sequence (any letters).")] string sequenceB,
        [Description("global (default), local or fogsaa.")] string mode = "global",
        [Description("Score of an identical pair (default 1; ignored with a matrix).")] double matchScore = 1,
        [Description("Score of a different pair (default 0; ignored with a matrix).")] double mismatchScore = 0,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("A single wildcard letter that scores 0 against every letter (scalar scoring only).")] string? wildcard = null,
        [Description("Scoring preset: blastn, megablast or blastp (matrix and gap scores of Bio.Align); individual settings below override it.")] string? scoring = null,
        [Description("Gap scores: a JSON object of Biopython property names to numbers (see the description).")] string? gapScoresJson = null,
        [Description("Gap-score functions of (position, length) (Biopython's Waterman-Smith-Beyer algorithm): a JSON object gap_score / insertion_score / deletion_score -> {\"kind\": lin | log | sq | pos | step | tab | zero | bonus, \"p\": [a, b, c]}: lin -(a + b*length), log -(a + b*ln(length + 1)), sq -(a + b*length^2), pos -(a + b*position + c*length), step -a up to c positions and -b beyond, tab -(((position*7 + length*3) mod b)*0.5 + a), zero 0, bonus a*length - b.")] string? gapFunctionsJson = null,
        [Description("+ (default) or -: align sequenceA with the reverse complement of sequenceB.")] string strand = "+",
        [Description("Return at most this many alignments (default 100; the number of optimal alignments can be exponential).")] int maxAlignments = 100,
        [Description("true: only the score (no alignments; memory O(n) for large inputs).")] bool scoreOnly = false,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequenceA, nameof(sequenceA));
        RequireSequence(sequenceB, nameof(sequenceB));
        if (strand is not ("+" or "-"))
            throw new ArgumentException("strand must be '+' or '-'.", nameof(strand));
        if (maxAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(maxAlignments), "maxAlignments must be at least 1.");
        bool presetScalars = !string.IsNullOrWhiteSpace(scoring);
        var aligner = BuildAligner(mode, presetScalars && matchScore == 1 ? null : matchScore, presetScalars && mismatchScore == 0 ? null : mismatchScore, matrix, wildcard, scoring, gapScoresJson, gapFunctionsJson);
        if (scoreOnly)
            return new PairwiseAlignerResult(aligner.Score(sequenceA, sequenceB, strand[0]), aligner.Algorithm, [.. aligner.Warnings], 0, []);
        var found = aligner.Align(sequenceA, sequenceB, strand[0]);
        var list = found.Take(maxAlignments).Select(a =>
        {
            var rows = AlignmentOperations.Rows(a);
            double s = a.Attributes["score"] is { } n ? n.GetValue<double>() : found.Score;
            return new PairwiseAlignmentDto(rows[0], rows[1], a.Coordinates, s, AlignmentOperations.Format(a, aligner.SubstitutionMatrix));
        }).ToArray();
        return new PairwiseAlignerResult(found.Score, aligner.Algorithm, [.. aligner.Warnings], found.Count, list);
    }

    [McpServerTool(Name = "alignment_counts", Title = "Alignment — Counts (Bio.Align AlignmentCounts)", ReadOnly = true)]
    [Description("The counts of an alignment of any number of sequences as Bio.Align's alignment.counts(scoring): aligned letters, identities, mismatches, positives (with a substitution matrix: letters with a positive score), and the gaps split by side (left / internal / right), direction (insertion = a gap in the target, deletion = a gap in the query) and opened / extended; the counts add up over every pair of sequences. Scoring: none (letter equality); wildcard (letters equal to it are neither identities nor mismatches); matrix (name or text; adds the substitution score); or an aligner (matchScore / mismatchScore / matrix / wildcard / gapScoresJson as in pairwise_aligner) which adds the gap score and the total score. " + AlignerParamNote + " VERIFICATION: equal to Biopython 1.88 on 260 alignments x 8 scorings (every count and the printed text).")]
    public static AlignmentCountsResult AlignmentCountsOf(
        [Description(AlignmentJsonDescription)] string alignmentJson,
        [Description("A single wildcard letter (neither identity nor mismatch).")] string? wildcard = null,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("Aligner scoring: score of an identical pair (giving any aligner setting scores the alignment like PairwiseAligner).")] double? matchScore = null,
        [Description("Aligner scoring: score of a different pair.")] double? mismatchScore = null,
        [Description("Aligner scoring: a JSON object of Biopython gap-score properties (see the description).")] string? gapScoresJson = null)
    {
        var a = ParseAlignmentJson(alignmentJson, nameof(alignmentJson));
        AlignmentCounts c;
        if (matchScore is not null || mismatchScore is not null || !string.IsNullOrWhiteSpace(gapScoresJson))
            c = AlignmentOperations.Counts(a, BuildAligner(null, matchScore, mismatchScore, matrix, wildcard, null, gapScoresJson));
        else if (!string.IsNullOrWhiteSpace(matrix))
            c = AlignmentOperations.Counts(a, ResolveMatrix(matrix)!);
        else if (!string.IsNullOrEmpty(wildcard))
            c = AlignmentOperations.Counts(a, wildcard.Length == 1 ? wildcard[0] : throw new ArgumentException("wildcard must be a single character.", nameof(wildcard)));
        else
            c = AlignmentOperations.Counts(a);
        return new AlignmentCountsResult(
            c.Aligned, c.Identities, c.Mismatches, c.Positives, c.Gaps, c.OpenGaps, c.ExtendGaps, c.LeftGaps, c.InternalGaps, c.RightGaps,
            c.OpenLeftInsertions, c.ExtendLeftInsertions, c.OpenLeftDeletions, c.ExtendLeftDeletions,
            c.OpenInternalInsertions, c.ExtendInternalInsertions, c.OpenInternalDeletions, c.ExtendInternalDeletions,
            c.OpenRightInsertions, c.ExtendRightInsertions, c.OpenRightDeletions, c.ExtendRightDeletions,
            c.GapScore, c.SubstitutionScore, c.Score, c.ToString());
    }

    [McpServerTool(Name = "alignment_operation", Title = "Alignment — Bio.Align Alignment operations (map, slice, indices, frequencies, ...)", ReadOnly = true)]
    [Description("One operation of Bio.Align's Alignment object on an alignment given by its block coordinates (see parse_coordinate_alignment): shape (rows, columns), rows (the printed rows, alignment[i]), row (alignment[row, columns] with a slice or a list), column (alignment[rows, col]), letter (alignment[row, col]), slice (alignment[rows], alignment[rows, columns] with a slice 'a:b:c' or a list [c1, c2, ...]: a new alignment), aligned (pairwise: the aligned blocks), indices (the index of the letter in every column, -1 for a gap), inverse_indices (the column of every letter, -1 if not aligned), frequencies (the weight of each letter / gap per column; the 'weight' annotation of a record), substitutions (the letter-pair counts), reverse_complement, sort (by identifier or residues; reverse), add (alignment + alignment2), map (alignmentJson.map(alignment2Json): reads on a transcript on a genome), mapall (alignmentJson = the multiple alignment, alignmentsJson = one alignment per row), from_same_reference (alignmentsJson: pairwise alignments with the same first sequence combined), parse_printed (rows: gapped rows -> sequences and coordinates), text (str(alignment): the Biopython pretty print, with an optional matrix for the ':' marks), equal / compare (alignment vs alignment2). Slices are Python slices ('1:3', ':', '::-1', '-2:'); negative indices count from the end. VERIFICATION: equal to Biopython 1.88 on 572 recorded cases (every operation, 12 index expressions per alignment).")]
    public static AlignmentOperationResult AlignmentOperationOn(
        [Description("The operation: shape, rows, row, column, letter, slice, aligned, indices, inverse_indices, frequencies, substitutions, reverse_complement, sort, add, map, mapall, from_same_reference, parse_printed, text, equal, compare.")] string operation,
        [Description(AlignmentJsonDescription + " Not needed for from_same_reference and parse_printed.")] string? alignmentJson = null,
        [Description("The second alignment (add, map, equal, compare).")] string? alignment2Json = null,
        [Description("A JSON array of alignments (mapall: one per row of alignmentJson; from_same_reference: the pairwise alignments).")] string[]? alignmentsJson = null,
        [Description("Row index (row, letter) or row slice 'a:b:c' (column, slice; default ':').")] string? rows = null,
        [Description("Column index (column, letter), column slice 'a:b:c' or a JSON list [1, 2] (row, slice).")] string? columns = null,
        [Description("sort: sort descending. text: nothing.")] bool reverse = false,
        [Description("Optional substitution matrix (name or text) for text: ':' marks positive pairs.")] string? matrix = null,
        [Description("parse_printed: the gapped rows.")] string[]? printedRows = null)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(operation);
        string op = operation.Trim().ToLowerInvariant();
        BioCoordAlignment A() => ParseAlignmentJson(Need(alignmentJson, nameof(alignmentJson)), nameof(alignmentJson));
        BioCoordAlignment B() => ParseAlignmentJson(Need(alignment2Json, nameof(alignment2Json)), nameof(alignment2Json));
        List<BioCoordAlignment> Many() => ParseMany(alignmentsJson, nameof(alignmentsJson));
        static string JStr(string v) => JsonValue.Create(v)!.ToJsonString();
        static string Strings(IEnumerable<string> v) => new JsonArray([.. v.Select(x => (JsonNode)JsonValue.Create(x)!)]).ToJsonString();
        static JsonArray Longs(IEnumerable<long> v) => new([.. v.Select(x => (JsonNode)JsonValue.Create(x)!)]);
        static string Ints(IEnumerable<int[]> v) => new JsonArray([.. v.Select(r => (JsonNode)new JsonArray([.. r.Select(x => (JsonNode)JsonValue.Create(x)!)]))]).ToJsonString();
        switch (op)
        {
            case "shape":
                var (r, c) = AlignmentOperations.Shape(A());
                return new AlignmentOperationResult(op, $"[{r.ToString(CultureInfo.InvariantCulture)}, {c.ToString(CultureInfo.InvariantCulture)}]");
            case "rows":
                return new AlignmentOperationResult(op, Strings(AlignmentOperations.Rows(A())));
            case "row":
            {
                var a = A();
                int row = int.Parse(rows ?? "0", CultureInfo.InvariantCulture);
                if (string.IsNullOrWhiteSpace(columns))
                    return new AlignmentOperationResult(op, JStr(AlignmentOperations.Row(a, row)));
                return new AlignmentOperationResult(op, JStr(columns.TrimStart().StartsWith('[') ? AlignmentOperations.Row(a, row, ParseIntList(columns)) : AlignmentOperations.Row(a, row, ParseSlice(columns))));
            }

            case "column":
                return new AlignmentOperationResult(op, JStr(AlignmentOperations.Column(A(), ParseSlice(rows), int.Parse(columns ?? "0", CultureInfo.InvariantCulture))));
            case "letter":
                return new AlignmentOperationResult(op, JStr(AlignmentOperations.Letter(A(), int.Parse(rows ?? "0", CultureInfo.InvariantCulture), int.Parse(columns ?? "0", CultureInfo.InvariantCulture)).ToString()));
            case "slice":
            {
                var a = A();
                var rowSlice = ParseSlice(rows);
                if (string.IsNullOrWhiteSpace(columns))
                    return new AlignmentOperationResult(op, AlignmentToJson(AlignmentOperations.Slice(a, rowSlice)));
                return new AlignmentOperationResult(op, AlignmentToJson(columns.TrimStart().StartsWith('[') ? AlignmentOperations.Slice(a, rowSlice, ParseIntList(columns)) : AlignmentOperations.Slice(a, rowSlice, ParseSlice(columns))));
            }

            case "aligned":
                return new AlignmentOperationResult(op, new JsonArray([.. AlignmentOperations.Aligned(A()).Select(row => (JsonNode)new JsonArray([.. row.Select(seg => (JsonNode)Longs(seg))]))]).ToJsonString());
            case "indices":
                return new AlignmentOperationResult(op, Ints(AlignmentOperations.Indices(A())));
            case "inverse_indices":
                return new AlignmentOperationResult(op, Ints(AlignmentOperations.InverseIndices(A())));
            case "frequencies":
            {
                var freq = AlignmentOperations.Frequencies(A());
                var o = new JsonObject();
                foreach (var (letter, values) in freq)
                    o[letter.ToString()] = new JsonArray([.. values.Select(x => (JsonNode)JsonValue.Create(x)!)]);
                return new AlignmentOperationResult(op, o.ToJsonString());
            }

            case "substitutions":
            {
                var s = AlignmentOperations.Substitutions(A());
                var values = new JsonArray();
                for (int i = 0; i < s.Letters.Length; i++)
                    values.Add((JsonNode)new JsonArray([.. Enumerable.Range(0, s.Letters.Length).Select(j => (JsonNode)JsonValue.Create(s.Values[i, j])!)]));
                return new AlignmentOperationResult(op, new JsonObject { ["letters"] = s.Letters, ["values"] = values }.ToJsonString());
            }

            case "reverse_complement":
                return new AlignmentOperationResult(op, AlignmentToJson(AlignmentOperations.ReverseComplement(A())));
            case "sort":
                return new AlignmentOperationResult(op, AlignmentToJson(AlignmentOperations.Sort(A(), reverse: reverse)));
            case "add":
                return new AlignmentOperationResult(op, AlignmentToJson(AlignmentOperations.Add(A(), B())));
            case "map":
                return new AlignmentOperationResult(op, AlignmentToJson(AlignmentOperations.Map(A(), B())));
            case "mapall":
                return new AlignmentOperationResult(op, AlignmentToJson(AlignmentOperations.MapAll(A(), Many())));
            case "from_same_reference":
                return new AlignmentOperationResult(op, AlignmentToJson(AlignmentOperations.FromAlignmentsWithSameReference(Many())));
            case "parse_printed":
            {
                var (sequences, coordinates) = AlignmentOperations.ParsePrintedAlignment(printedRows ?? throw new ArgumentException("printedRows is required.", nameof(printedRows)));
                var o = new JsonObject
                {
                    ["sequences"] = new JsonArray([.. sequences.Select(x => (JsonNode)JsonValue.Create(x)!)]),
                    ["coordinates"] = new JsonArray([.. coordinates.Select(row => (JsonNode)Longs(row))]),
                };
                return new AlignmentOperationResult(op, o.ToJsonString());
            }

            case "text":
            {
                string text = AlignmentOperations.Format(A(), ResolveMatrix(matrix));
                return new AlignmentOperationResult(op, JStr(text), text);
            }

            case "equal":
                return new AlignmentOperationResult(op, AlignmentOperations.AreEqual(A(), B()) ? "true" : "false");
            case "compare":
                return new AlignmentOperationResult(op, AlignmentOperations.Compare(A(), B()).ToString(CultureInfo.InvariantCulture));
            default:
                throw new ArgumentException($"Unknown operation '{operation}'.", nameof(operation));
        }
    }
}
