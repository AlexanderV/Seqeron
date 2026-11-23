using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
using ModelContextProtocol.Server;
using global::Seqeron.Genomics.Alignment;
using global::Seqeron.Genomics.Analysis;

namespace Seqeron.Mcp.Analysis.Tools;

// B09 WP67 (F308, R14-13): Bio.Align.AlignInfo.SummaryInfo of Biopython 1.79 and the conservation scores of Capra and Singh (2007) as MCP tools.
public partial class AnalysisTools
{
    private const string MsaRecordsDescription = "The rows as a JSON array of records {\"id\": ..., \"seq\": \"ACGT-\", \"annotations\": {\"weight\": 0.5}}; only seq is required (id defaults to <unknown id>); the 'weight' annotation of a row weights it (default 1.0).";

    private static string NeedText(string? value, string name) => string.IsNullOrWhiteSpace(value) ? throw new ArgumentException($"{name} is required for this operation.", name) : value;

    private static MultipleSeqAlignment MsaFromParts(string? recordsJson, string parameter)
    {
        var records = JsonNode.Parse(NeedText(recordsJson, parameter)) as JsonArray ?? throw new ArgumentException("The records are a JSON array of record objects.", parameter);
        return new MultipleSeqAlignment(records.Select(r =>
        {
            var o = r as JsonObject ?? throw new ArgumentException("A record is a JSON object {\"id\", \"seq\", ...}.", parameter);
            return new SeqRecord(o["seq"]?.GetValue<string>() ?? throw new ArgumentException("A record needs a seq.", parameter), o["id"]?.GetValue<string>(), o["name"]?.GetValue<string>(), o["description"]?.GetValue<string>(),
                o["annotations"] as JsonObject, o["letter_annotations"] as JsonObject, (o["dbxrefs"] as JsonArray)?.Select(x => x!.GetValue<string>()));
        }));
    }

    private static string[] ParseStringArray(string? json, string parameter) =>
        string.IsNullOrWhiteSpace(json) ? [] : [.. (JsonNode.Parse(json) as JsonArray ?? throw new ArgumentException("A JSON array of strings is needed, e.g. [\"-\", \"*\"].", parameter)).Select(x => x!.GetValue<string>())];

    [McpServerTool(Name = "alignment_summary_info", Title = "Alignment — Biopython AlignInfo SummaryInfo (consensus, PSSM, information content, replacements)", ReadOnly = true)]
    [Description("One member of Bio.Align.AlignInfo.SummaryInfo of Biopython 1.79 (the last release that has the class; 1.88 keeps only get_column) on an alignment of records: dumb_consensus (per column the unique most common residue when its share of the residues, gaps not counted, reaches threshold, else the ambiguous letter; requireMultiple), gap_consensus (the same with the gap symbols counted), replacement_dictionary (letters: for every ordered pair of letters the summed product of the 'weight' annotations of the two rows that carry them in one column, the earlier row's letter first), pos_specific_score_matrix (the PSSM: per column the weighted count of every letter, with the consensus or axisSeq as the axis; charsToIgnore, the gap symbol is always ignored; the table text is str(pssm)), information_content (sum over the columns start..end-1 of the sum of f ln(f / e) / ln(logBase) with f the weighted frequency of the letter, gaps and charsToIgnore not counted, e from eFreqTableJson, pseudoCount k: f = (count + e k) / (total + k); the per-column values are ic_vector), ic_vector (the vector after information_content), print_info_content (position, letter of row repRecord, value with three decimals; needs eFreqTableJson), get_column. Errors and their messages are Biopython's (information_content without an expected-frequency table fails as it does there: TypeError; a negative pseudo count, a letter missing in the table, a range outside the first row, logBase 1 ...). VERIFICATION: the values of Biopython's own test suite (Tests/test_AlignInfo.py: ic 7.32029999423075, 133.061475107, 7.546, the printed PSSMs) and hand-derived values; the recorded run of 6 199 calls on the real class is the user-run package handoff/R14-13-aligninfo-biopython179 (its README says whether it has been integrated).")]
    public static AlignmentSummaryInfoResult AlignmentSummaryInfoOperation(
        [Description("The operation: dumb_consensus, gap_consensus, replacement_dictionary, pos_specific_score_matrix, information_content, ic_vector, print_info_content, get_column.")] string operation,
        [Description(MsaRecordsDescription)] string recordsJson,
        [Description("consensus: share the most common residue needs (default 0.7).")] double threshold = 0.7,
        [Description("consensus: the letter of a column without a clear consensus (default X).")] string ambiguous = "X",
        [Description("consensus: a column with one residue only is ambiguous.")] bool requireMultiple = false,
        [Description("replacement_dictionary: the letters to count (required).")] string? letters = null,
        [Description("pos_specific_score_matrix: the axis sequence (default the dumb consensus).")] string? axisSeq = null,
        [Description("pos_specific_score_matrix / information_content: a JSON array of letters that are not counted, e.g. [\"-\", \"*\"].")] string? charsToIgnoreJson = null,
        [Description("information_content: first column (0-based, default 0).")] int start = 0,
        [Description("information_content: end column, exclusive (default the length of the first row).")] int? end = null,
        [Description("information_content: expected frequency of each letter as a JSON object {\"A\": 0.25, \"C\": 0.25, ...} (gap symbols have none).")] string? eFreqTableJson = null,
        [Description("information_content: logarithm base (default 2: bits).")] double logBase = 2,
        [Description("information_content: pseudo count k (default 0).")] double pseudoCount = 0,
        [Description("get_column: the column (negative from the end).")] int col = 0,
        [Description("print_info_content: the row whose letters label the positions (default 0).")] int repRecord = 0)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(operation);
        string op = operation.Trim().ToLowerInvariant();
        var m = MsaFromParts(recordsJson, nameof(recordsJson));
        var summary = new SummaryInfo(m);
        var ignore = ParseStringArray(charsToIgnoreJson, nameof(charsToIgnoreJson));
        IReadOnlyDictionary<char, double>? table = null;
        if (!string.IsNullOrWhiteSpace(eFreqTableJson))
        {
            var o = JsonNode.Parse(eFreqTableJson) as JsonObject ?? throw new ArgumentException("eFreqTableJson must be a JSON object {\"A\": 0.25, ...}.", nameof(eFreqTableJson));
            if (o.Any(kv => kv.Key.Length != 1))
                throw new ArgumentException("Every key of eFreqTableJson is a single letter.", nameof(eFreqTableJson));
            table = o.ToDictionary(kv => kv.Key[0], kv => kv.Value!.GetValue<double>());
        }

        static string JStr(string v) => JsonValue.Create(v)!.ToJsonString();
        static JsonNode Num(double v) => JsonValue.Create(v)!;
        switch (op)
        {
            case "dumb_consensus": return new(op, JStr(summary.DumbConsensus(threshold, ambiguous, requireMultiple)));
            case "gap_consensus": return new(op, JStr(summary.GapConsensus(threshold, ambiguous, requireMultiple)));
            case "replacement_dictionary":
            {
                var rep = summary.ReplacementDictionary(letters);
                return new(op, new JsonArray([.. rep.Select(kv => (JsonNode)new JsonArray(JsonValue.Create(kv.Key.First.ToString()), JsonValue.Create(kv.Key.Second.ToString()), Num(kv.Value)))]).ToJsonString());
            }

            case "pos_specific_score_matrix":
            {
                var list = ignore.ToList();
                var pssm = summary.PosSpecificScoreMatrix(axisSeq, list);
                var rows = new JsonArray();
                for (int i = 0; i < pssm.Count; i++)
                {
                    var counts = new JsonObject();
                    foreach (var (letter, v) in pssm[i])
                        counts[letter.ToString()] = v;
                    rows.Add((JsonNode)new JsonObject { ["letter"] = pssm.GetResidue(i), ["counts"] = counts });
                }

                string text = pssm.ToString();
                return new(op, new JsonObject { ["rows"] = rows, ["text"] = text }.ToJsonString(), text);
            }

            case "information_content":
            {
                double v = summary.InformationContent(start, end, table, logBase, ignore, pseudoCount);
                return new(op, new JsonObject { ["value"] = v, ["ic_vector"] = new JsonArray([.. summary.IcVector.Select(Num)]) }.ToJsonString());
            }

            case "ic_vector":
                summary.InformationContent(start, end, table, logBase, ignore, pseudoCount);
                return new(op, new JsonArray([.. summary.IcVector.Select(Num)]).ToJsonString());
            case "print_info_content":
            {
                if (table is not null)
                    summary.InformationContent(start, end, table, logBase, ignore, pseudoCount);
                string text = summary.PrintInfoContent(repRecord);
                return new(op, JStr(text), text);
            }

            case "get_column":
                return new(op, JStr(summary.GetColumn(col)));
            default:
                throw new ArgumentException($"Unknown operation '{operation}'.", nameof(operation));
        }
    }

    [McpServerTool(Name = "msa_conservation", Title = "Alignment — per-column conservation scores (Capra and Singh 2007)", ReadOnly = true)]
    [Description("The conservation score of every column of a protein alignment as in the program of Capra and Singh (2007, Bioinformatics 23:1875): method js_divergence (Jensen-Shannon divergence in bits of the column distribution from the background, 0 for a column that looks like the background, close to 1 for a conserved column of a rare residue; the default method), shannon_entropy (the information score 1 - h with h = sum f ln f / ln(min(21, rows)): 1 conserved, 0 maximal entropy), property_entropy (the same over the groups AVLIMC, FWYH, STNQ, KR, DE, GP and the gap), relative_entropy (Kullback-Leibler divergence from the background / ln 20) or sum_of_pairs (the weighted mean BLOSUM62 score of all pairs of residues of the column). Frequencies carry the pseudo count 1e-7; gapPenalty (default true) multiplies each score by the weighted fraction of non-gap rows; sequenceWeights weights the rows by the position-based weights of Henikoff and Henikoff (1994); windowSize smooths the scores over that many columns on either side (the score of a column 0.5 of itself + 0.5 of the mean of its non-negative neighbours); background blosum62 (default) or swissprot, or a JSON array of 20 / 21 numbers in the order A R N D C Q E G H I L K M F P S T W Y V (-). Rows are upper-cased and read as the program's readers do (B = D, Z = Q, X and . = gap); a letter that is not one of the 20 amino acids is ignored by the entropies and is an error for sum_of_pairs and the weights. VERIFICATION: equal to an independent implementation of the formulas of the paper (every score of the unit tests to 1e-13) and to hand-derived values (the program itself is a Python 2 CGI script that is not run here).")]
    public static MsaConservationResult MsaConservationScores(
        [Description("At least one aligned protein sequence; equal lengths, gaps '-'.")] string[] sequences,
        [Description("js_divergence (default), shannon_entropy, property_entropy, relative_entropy or sum_of_pairs.")] string method = "js_divergence",
        [Description("Background distribution: blosum62 (default), swissprot, or a JSON array of 20 / 21 numbers.")] string background = "blosum62",
        [Description("Weight the rows (Henikoff and Henikoff 1994).")] bool sequenceWeights = false,
        [Description("Multiply every score by the weighted fraction of non-gap rows (default true).")] bool gapPenalty = true,
        [Description("Window size on either side (default 0: no smoothing).")] int windowSize = 0)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required", nameof(sequences));
        var m = method.Trim().ToLowerInvariant() switch
        {
            "js_divergence" or "jsd" => ConservationMethod.JsDivergence,
            "shannon_entropy" or "shannon" => ConservationMethod.ShannonEntropy,
            "property_entropy" => ConservationMethod.PropertyEntropy,
            "relative_entropy" => ConservationMethod.RelativeEntropy,
            "sum_of_pairs" => ConservationMethod.SumOfPairs,
            _ => throw new ArgumentException("method: js_divergence, shannon_entropy, property_entropy, relative_entropy or sum_of_pairs.", nameof(method)),
        };
        IReadOnlyList<double> bg = background.Trim().ToLowerInvariant() switch
        {
            "blosum62" or "" => ConservationScores.BlosumBackground,
            "swissprot" => ConservationScores.SwissProtBackground,
            var json when json.StartsWith('[') => [.. (JsonNode.Parse(background) as JsonArray ?? throw new ArgumentException("background must be a JSON array of numbers.", nameof(background))).Select(x => x!.GetValue<double>())],
            _ => throw new ArgumentException("background: blosum62, swissprot or a JSON array of 20 / 21 numbers.", nameof(background)),
        };
        var scores = global::Seqeron.Genomics.Analysis.MultipleAlignmentStatistics.ConservationProfile(sequences, m, bg, null, sequenceWeights, gapPenalty, windowSize);
        var rows = sequences.Select(r => r.ToUpperInvariant().Replace('B', 'D').Replace('Z', 'Q').Replace('X', '-').Replace('.', '-')).ToList();
        double[] weights = sequenceWeights ? ConservationScores.SequenceWeights(rows) : [];
        string[] columns = [.. Enumerable.Range(0, rows[0].Length).Select(i => new string([.. rows.Select(r => r[i])]))];
        return new MsaConservationResult(m.ToString(), windowSize, scores, weights, columns);
    }
}
