using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Analysis.Tools;

// EMBOSS program ports on the Analysis server (B09 A3-8 / A3-12 / A3-13); every tool delegates to GenomicAnalyzer.
public partial class AnalysisTools
{
    [McpServerTool(Name = "find_tandem_repeats_emboss", Title = "Genome — EMBOSS etandem / equicktandem", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 tandem-repeat programs ported statement by statement (diagonal-score models, not the exact-run finder find_tandem_repeats). program etandem (default): for every unit size minRepeat..maxRepeat a running per-column consensus score (+1 when the new base is the column's most frequent base, -1 otherwise; N counts as a fifth base, or as a plain -1 with mismatch) extends candidates in steps of the unit size; candidates need a best score >= threshold and, unless uniform, a consensus with two different bases; overlapping candidates lose to the highest score (ties: smaller unit). Reports score, size, count, identity = 100 * (score + n * (copies + 1)) / (2 * n * copies) and the lower-case consensus (phased as printed). program equicktandem: the quick pre-screen - the sequence compared with itself shifted by every gap 1..maxRepeat with a +1 / -1 diagonal score; reported stretches are masked for larger gaps; no consensus / identity. Hits are in the program's report order (etandem: highest score first). Non-ACGT symbols count as N (etandem) / one 'other' class (equicktandem). position is 0-based (EMBOSS Start - 1). Verified against etandem and equicktandem 6.6.0 on 1700 + 1300 random and planted sequences.")]
    public static FindEmbossTandemRepeatsResult FindTandemRepeatsEmboss(
        [Description("Nucleotide sequence: A/C/G/T plus IUPAC codes such as N (case-insensitive).")] string sequence,
        [Description("etandem (default) | equicktandem (case-insensitive).")] string program = "etandem",
        [Description("etandem -minrepeat: smallest unit size (>= 2, default 10 as in the program). Ignored by equicktandem.")] int minRepeat = 10,
        [Description("Largest unit size (-maxrepeat). Default: etandem = minRepeat (as in the program), equicktandem = 600.")] int? maxRepeat = null,
        [Description("Minimum score (-threshold, default 20; equicktandem needs >= 1).")] int threshold = 20,
        [Description("etandem -mismatch: count N as a mismatch instead of a fifth base (default false).")] bool mismatch = false,
        [Description("etandem -uniform: also report single-base consensus repeats (default false).")] bool uniform = false,
        CancellationToken cancellationToken = default)
    {
        var nucleotides = RequireIupacDna(sequence, nameof(sequence));
        string which = (program ?? "etandem").Trim().ToLowerInvariant();
        IReadOnlyList<EmbossTandemRepeat> found;
        switch (which)
        {
            case "etandem":
                found = GenomicAnalyzer.FindTandemRepeatsEtandem(nucleotides, minRepeat, maxRepeat ?? minRepeat, threshold, mismatch, uniform, cancellationToken);
                break;
            case "equicktandem":
                found = GenomicAnalyzer.FindTandemRepeatsEquicktandem(nucleotides, maxRepeat ?? 600, threshold, cancellationToken);
                break;
            default:
                throw new ArgumentException("program must be etandem or equicktandem", nameof(program));
        }

        bool detailed = which == "etandem";
        var items = found
            .Select(h => new EmbossTandemItem(h.Position, h.Length, h.Score, h.RepeatSize, h.RepeatCount,
                detailed ? Math.Round((double)h.Identity, 1) : null, h.Consensus))
            .ToArray();
        return new FindEmbossTandemRepeatsResult(which, items);
    }

    [McpServerTool(Name = "seqmatchall", Title = "Genome — EMBOSS seqmatchall (all-against-all word matches)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 seqmatchall: the all-against-all comparison of a set of sequences by exact word matches. For every pair i < j with length(i) > wordSize and length(j) >= wordSize every maximal exact match of at least wordSize symbols is reported (forward strand only; overlapping matches and repeated positions are all listed). Nucleotide or protein sequences (any letters, compared literally and case-insensitively); the program's word table leaves out words around its unknown symbol (N for a nucleotide sequence, X otherwise) with a scan that can keep a word that holds an N after a skipped one - reproduced exactly. Positions are 0-based (the program prints 1-based starts). Matches are listed per pair, by position (the program lists them in hash order: compare as sets). Verified against seqmatchall 6.6.0 on 330 random sets (2-5 sequences, nucleotide / N / lower case / protein, word sizes 2-12, 700+ pairs); on A/C/G/T input it equals find_common_regions mode MaximalExactMatches (MUMmer -maxmatch).")]
    public static SeqMatchAllResult SeqMatchAll(
        [Description("The sequences (at least one; nucleotide or protein letters).")] string[] sequences,
        [Description("Word size (seqmatchall -wordsize, default 4, >= 1).")] int wordSize = 4,
        CancellationToken cancellationToken = default)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required", nameof(sequences));
        if (sequences.Any(string.IsNullOrEmpty))
            throw new ArgumentException("Sequences cannot be null or empty", nameof(sequences));
        var items = GenomicAnalyzer.SeqMatchAll(sequences, wordSize, cancellationToken)
            .Select(m => new WordMatchItem(m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length))
            .ToArray();
        return new SeqMatchAllResult(items);
    }
}
