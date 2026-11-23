using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 R14-11 (F304): the score API of parasail 1.3.4 (score-only with end coordinates, statistics without a traceback, tables, last row and column; one query against many references).
public partial class AlignmentTools
{
    private const long ParasailMaxTableCells = 4_000_000;

    private static int[][] ToJagged(int[,] a)
    {
        var rows = new int[a.GetLength(0)][];
        for (int i = 0; i < rows.Length; i++)
        {
            rows[i] = new int[a.GetLength(1)];
            for (int j = 0; j < rows[i].Length; j++)
                rows[i][j] = a[i, j];
        }

        return rows;
    }

    [McpServerTool(Name = "parasail_score", Title = "Pairwise — parasail score API (score, end coordinates, statistics, tables)", ReadOnly = true)]
    [Description("The score API of parasail 1.3.4 for one query against one or many references: the optimal score with its end coordinates (endQuery, endReference: 0-based index of the last aligned residue), the statistics of the optimal alignment WITHOUT a traceback (matches, similar = positive pairs, length: parasail nw_stats / sw_stats / sg_stats), the full score tables (table; with statistics also the matches / similar / length tables) and the last row and column of the score matrix (rowcol), for global (function nw), local (sw) and every semi-global variant (sg = all end gaps free, sg_qb / sg_qe / sg_db / sg_de, sg_qx, sg_dx and the combinations sg_qb_db, sg_qb_de, sg_qe_db, sg_qe_de: q = query, d = database = reference; b = beginning, e = end, x = both ends). Gaps are parasail's affine penalties: a gap of length k costs open + (k - 1) * extend (both >= 0; BLAST-style 11 / 1 is open 11, extend 1). The program keeps O(min(length)) cells (the table needs m * n) and the query is prepared once for the whole list of references (profile reuse); parallel runs the references on all cores. Scoring: matrix (name or text; integer scores) or match / mismatch over the letters of the sequences. VERIFICATION: equal to parasail 1.3.4 (nw, sw and the eleven sg variants x plain, _stats, _table, _rowcol, _stats_table, _stats_rowcol) on 520 / 520 recorded pairs.")]
    public static ParasailScoreResult ParasailScore(
        [Description("The query (parasail s1).")] string query,
        [Description("The references (parasail s2, the database sequences), at least one.")] string[] references,
        [Description("parasail function: nw (global), sw (local) or a semi-global variant sg, sg_qb, sg_qe, sg_db, sg_de, sg_qx, sg_dx, sg_qb_db, sg_qb_de, sg_qe_db, sg_qe_de (default sw).")] string function = "sw",
        [Description("Penalty of the first gap position (>= 0; default 11).")] int open = 11,
        [Description("Penalty of each further gap position (>= 0; default 1).")] int extend = 1,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description("Score of an identical pair when no matrix is given (default 2).")] int match = 2,
        [Description("Score of a different pair when no matrix is given (default -3).")] int mismatch = -3,
        [Description("true: matches, similar and length of the optimal alignment (the *_stats functions).")] bool statistics = false,
        [Description("true: the score tables (and, with statistics, the matches / similar / length tables); at most 4 000 000 cells per reference.")] bool table = false,
        [Description("true: the last row and the last column of the score matrix (the *_rowcol functions).")] bool rowcol = false,
        [Description("true: score the references in parallel (the results keep their order).")] bool parallel = false,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(query, nameof(query));
        if (references is null || references.Length == 0)
            throw new ArgumentException("At least one reference is required.", nameof(references));
        foreach (string r in references)
            RequireSequence(r, nameof(references));
        string fn = (function ?? "sw").Trim().ToLowerInvariant();
        ParasailMode mode;
        var ends = default(ParasailEndGaps);
        if (fn == "nw")
            mode = ParasailMode.Global;
        else if (fn == "sw")
            mode = ParasailMode.Local;
        else if (ParasailEndGaps.TryParse(fn, out ends))
            mode = ParasailMode.SemiGlobal;
        else
            throw new ArgumentException($"Unknown function '{function}': nw, sw, sg, sg_qb, sg_qe, sg_db, sg_de, sg_qx, sg_dx, sg_qb_db, sg_qb_de, sg_qe_db or sg_qe_de.", nameof(function));
        var sub = ResolveMatrix(matrix);
        var pm = sub is not null ? ParasailMatrix.From(sub) : ParasailMatrix.Create(references.Prepend(query), match, mismatch);
        if (table && references.Any(r => (long)query.Length * r.Length > ParasailMaxTableCells))
            throw new ArgumentException($"A table needs query length x reference length <= {ParasailMaxTableCells} cells; use statistics / rowcol or shorter sequences.", nameof(table));
        var profile = new ParasailProfile(query, pm);
        var output = new ParasailOutput(statistics, table, rowcol);
        var results = profile.ScoreBatch(references, open, extend, mode, ends, output, parallel);
        var items = results.Select(r => new ParasailScoreItem(
            r.Score, r.EndQuery, r.EndReference, r.Matches, r.Similar, r.Length,
            r.Tables is { } t ? ToJagged(t.Score) : null, r.Tables?.Matches is { } tm ? ToJagged(tm) : null, r.Tables?.Similar is { } ts ? ToJagged(ts) : null, r.Tables?.Length is { } tl ? ToJagged(tl) : null,
            r.RowCol?.ScoreRow, r.RowCol?.ScoreColumn, r.RowCol?.MatchesRow, r.RowCol?.MatchesColumn, r.RowCol?.SimilarRow, r.RowCol?.SimilarColumn, r.RowCol?.LengthRow, r.RowCol?.LengthColumn)).ToArray();
        return new ParasailScoreResult(fn, items.Length, items);
    }
}
