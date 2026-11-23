using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>The kind of dynamic program of a <see cref="ParasailAligner"/> call (parasail <c>nw</c>, <c>sw</c>, <c>sg</c> and the <c>sg_*</c> variants).</summary>
public enum ParasailMode
{
    /// <summary>Global (<c>parasail_nw</c>): the score of the whole of both sequences; the end is (length - 1, length - 1).</summary>
    Global,

    /// <summary>Local (<c>parasail_sw</c>): the best-scoring pair of substrings; the end is the cell with the best score of the smallest reference index (then the smallest query index), as parasail's loop picks it.</summary>
    Local,

    /// <summary>Semi-global (<c>parasail_sg_flags</c>): end gaps of the sequences that <see cref="ParasailEndGaps"/> frees cost nothing; the end is the best cell of the last row / column.</summary>
    SemiGlobal,
}

/// <summary>
/// The end gaps of a semi-global call that cost nothing (parasail <c>sg_flags</c>: <c>s1_beg</c>, <c>s1_end</c>, <c>s2_beg</c>, <c>s2_end</c>; sequence 1 is the query, sequence 2 the reference). The names of the library's variants:
/// <c>sg</c> (all free), <c>sg_qb</c> / <c>sg_qe</c> / <c>sg_db</c> / <c>sg_de</c> (query / database begin / end free), <c>sg_qx</c> (both query ends), <c>sg_dx</c> (both database ends), and the combinations <c>sg_qb_db</c> ... <c>sg_qe_de</c>.
/// </summary>
/// <param name="S1Begin">The leading residues of the query are free (<c>sg_qb</c>).</param>
/// <param name="S1End">The trailing residues of the query are free (<c>sg_qe</c>).</param>
/// <param name="S2Begin">The leading residues of the reference are free (<c>sg_db</c>).</param>
/// <param name="S2End">The trailing residues of the reference are free (<c>sg_de</c>).</param>
public readonly record struct ParasailEndGaps(bool S1Begin = false, bool S1End = false, bool S2Begin = false, bool S2End = false)
{
    /// <summary>All four end gaps free (<c>sg</c>).</summary>
    public static ParasailEndGaps All => new(true, true, true, true);

    /// <summary>The flags of a parasail function name (<c>sg</c>, <c>sg_qb</c>, <c>sg_qx</c>, <c>sg_dx</c>, <c>sg_qb_de</c>, ...); false when the name is not a semi-global variant.</summary>
    public static bool TryParse(string? name, out ParasailEndGaps flags)
    {
        flags = default;
        if (name is null)
            return false;
        string n = name.Trim().ToLowerInvariant();
        if (n == "sg")
        {
            flags = All;
            return true;
        }

        if (!n.StartsWith("sg_", StringComparison.Ordinal))
            return false;
        bool qb = false, qe = false, db = false, de = false;
        foreach (string part in n[3..].Split('_'))
        {
            switch (part)
            {
                case "qb": qb = true; break;
                case "qe": qe = true; break;
                case "db": db = true; break;
                case "de": de = true; break;
                case "qx": qb = qe = true; break;
                case "dx": db = de = true; break;
                default: return false;
            }
        }

        flags = new ParasailEndGaps(qb, qe, db, de);
        return true;
    }
}

/// <summary>A substitution matrix of integer scores for a <see cref="ParasailAligner"/> (parasail <c>parasail_matrix_t</c>): the letters of the alphabet and their pair scores. Letters outside the alphabet are refused.</summary>
public sealed class ParasailMatrix
{
    private readonly int[,] _scores;
    private readonly int[] _index = new int[128];

    private ParasailMatrix(string name, string alphabet, int[,] scores)
    {
        Name = name;
        Alphabet = alphabet;
        _scores = scores;
        Array.Fill(_index, -1);
        for (int i = 0; i < alphabet.Length; i++)
            _index[alphabet[i]] = i;
    }

    /// <summary>The name of the matrix.</summary>
    public string Name { get; }

    /// <summary>The letters (upper case).</summary>
    public string Alphabet { get; }

    /// <summary>The matrix of <paramref name="matrix"/> (every score must be an integer).</summary>
    /// <exception cref="ArgumentException">The matrix has fractional scores.</exception>
    public static ParasailMatrix From(SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        if (matrix.IsFractional)
            throw new ArgumentException($"The matrix {matrix.Name} has fractional scores; the parasail-style score API takes integer scores.", nameof(matrix));
        string alphabet = matrix.Alphabet.ToUpperInvariant();
        var scores = new int[alphabet.Length, alphabet.Length];
        for (int i = 0; i < alphabet.Length; i++)
            for (int j = 0; j < alphabet.Length; j++)
                scores[i, j] = matrix.TryGetScore(alphabet[i], alphabet[j], out int s) ? s : 0;
        return new ParasailMatrix(matrix.Name, alphabet, scores);
    }

    /// <summary>The identity matrix of <paramref name="alphabet"/>: <paramref name="match"/> on the diagonal, <paramref name="mismatch"/> elsewhere (parasail <c>matrix_create</c>).</summary>
    public static ParasailMatrix Create(string alphabet, int match, int mismatch)
    {
        ArgumentException.ThrowIfNullOrEmpty(alphabet);
        string a = new([.. alphabet.ToUpperInvariant().Distinct()]);
        var scores = new int[a.Length, a.Length];
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < a.Length; j++)
                scores[i, j] = i == j ? match : mismatch;
        return new ParasailMatrix($"{match}/{mismatch}", a, scores);
    }

    /// <summary>The identity matrix over the letters of the sequences given (and nothing else).</summary>
    public static ParasailMatrix Create(IEnumerable<string> sequences, int match, int mismatch) => Create(string.Concat(sequences), match, mismatch);

    internal int Code(char c)
    {
        char u = char.ToUpperInvariant(c);
        int index = u < 128 ? _index[u] : -1;
        return index >= 0 ? index : throw new ArgumentException($"The letter '{c}' is not in the alphabet '{Alphabet}' of the matrix {Name}.");
    }

    internal int[] Encode(string sequence)
    {
        var codes = new int[sequence.Length];
        for (int i = 0; i < codes.Length; i++)
            codes[i] = Code(sequence[i]);
        return codes;
    }

    internal int Score(int a, int b) => _scores[a, b];
}

/// <summary>The tables of a <see cref="ParasailAligner"/> call that asked for them: <c>[i, j]</c> is the cell of query residue i and reference residue j (0-based), the values of <c>score_table</c> and, with statistics, <c>matches_table</c>, <c>similar_table</c>, <c>length_table</c>.</summary>
/// <param name="Score">The score of every cell.</param>
/// <param name="Matches">The identities of the alignment ending in each cell (statistics only).</param>
/// <param name="Similar">The positive pairs of that alignment.</param>
/// <param name="Length">The columns of that alignment.</param>
public sealed record ParasailTables(int[,] Score, int[,]? Matches, int[,]? Similar, int[,]? Length);

/// <summary>The last row and last column of the score matrix (<c>score_row</c> / <c>score_col</c> of the <c>_rowcol</c> functions) with their statistics.</summary>
/// <param name="ScoreRow">The scores of the last query row, per reference residue.</param>
/// <param name="ScoreColumn">The scores of the last reference column, per query residue.</param>
/// <param name="MatchesRow">Identities along the row (statistics only).</param>
/// <param name="MatchesColumn">Identities along the column.</param>
/// <param name="SimilarRow">Positive pairs along the row.</param>
/// <param name="SimilarColumn">Positive pairs along the column.</param>
/// <param name="LengthRow">Alignment columns along the row.</param>
/// <param name="LengthColumn">Alignment columns along the column.</param>
public sealed record ParasailRowCol(int[] ScoreRow, int[] ScoreColumn, int[]? MatchesRow, int[]? MatchesColumn, int[]? SimilarRow, int[]? SimilarColumn, int[]? LengthRow, int[]? LengthColumn);

/// <summary>The result of a <see cref="ParasailAligner"/> call (parasail <c>result</c>): score and end coordinates, the statistics of the optimal alignment when asked, the tables or the last row and column.</summary>
/// <param name="Score">The optimal score.</param>
/// <param name="EndQuery">0-based index of the last query residue of the best alignment (parasail <c>end_query</c>).</param>
/// <param name="EndReference">0-based index of the last reference residue (<c>end_ref</c>).</param>
/// <param name="Matches">Identical aligned pairs of the alignment (<c>matches</c>; statistics only).</param>
/// <param name="Similar">Aligned pairs with a positive score (<c>similar</c>).</param>
/// <param name="Length">Columns of the alignment, gaps included (<c>length</c>).</param>
/// <param name="Tables">The tables (<c>_table</c> functions) or null.</param>
/// <param name="RowCol">The last row and column (<c>_rowcol</c> functions) or null.</param>
public sealed record ParasailResult(int Score, int EndQuery, int EndReference, int? Matches, int? Similar, int? Length, ParasailTables? Tables, ParasailRowCol? RowCol);

/// <summary>What a <see cref="ParasailAligner"/> call returns besides the score and the end coordinates.</summary>
/// <param name="Statistics">Matches, similar and length of the optimal alignment without a traceback (<c>*_stats</c>).</param>
/// <param name="Table">The tables of every cell (<c>*_table</c>; O(m·n) memory).</param>
/// <param name="RowCol">The last row and column (<c>*_rowcol</c>).</param>
public readonly record struct ParasailOutput(bool Statistics = false, bool Table = false, bool RowCol = false);

/// <summary>
/// The score-only API of parasail 1.3.4 (reference <c>nw</c> / <c>sw</c> / <c>sg_flags</c> of the library, scalar code) on the library's own sequences: the optimal score with its end coordinates (<c>end_query</c>, <c>end_ref</c>), the statistics of the optimal alignment
/// without a traceback (<c>matches</c>, <c>similar</c>, <c>length</c>: <c>nw_stats</c>, <c>sw_stats</c>, <c>sg_stats</c>), the full score tables and the last row and column (<c>_table</c>, <c>_rowcol</c>), for global, local and every semi-global variant with affine gaps.
/// Gap costs are parasail's: a gap of length k costs <c>open + (k - 1) * gap</c> (both penalties &gt;= 0). The program keeps O(min(m, n)) cells (the longer sequence is the outer loop; ties are resolved exactly as the original loop order would)
/// unless a table is asked for. <see cref="ParasailProfile"/> keeps a query for many references.
/// </summary>
public static class ParasailAligner
{
    private const int NegInf = int.MinValue / 2;

    /// <summary>Runs one parasail-style program.</summary>
    /// <param name="query">Sequence 1 (parasail <c>s1</c>).</param>
    /// <param name="reference">Sequence 2 (<c>s2</c>, the database sequence).</param>
    /// <param name="open">Penalty of the first gap position (&gt;= 0).</param>
    /// <param name="gap">Penalty of each further gap position (&gt;= 0).</param>
    /// <param name="matrix">The scoring matrix.</param>
    /// <param name="mode">Global, local or semi-global.</param>
    /// <param name="ends">The free end gaps of a semi-global alignment.</param>
    /// <param name="output">What to compute besides the score.</param>
    /// <exception cref="ArgumentException">An empty sequence, a negative penalty, a letter outside the matrix alphabet.</exception>
    public static ParasailResult Run(string query, string reference, int open, int gap, ParasailMatrix matrix, ParasailMode mode, ParasailEndGaps ends = default, ParasailOutput output = default)
    {
        ArgumentNullException.ThrowIfNull(query);
        ArgumentNullException.ThrowIfNull(reference);
        ArgumentNullException.ThrowIfNull(matrix);
        return Run(matrix.Encode(query), matrix.Encode(reference), open, gap, matrix, mode, ends, output);
    }

    /// <summary><c>parasail_nw</c> (<c>nw_stats</c>, <c>nw_table</c> ... with <paramref name="output"/>).</summary>
    public static ParasailResult Nw(string query, string reference, int open, int gap, ParasailMatrix matrix, ParasailOutput output = default) =>
        Run(query, reference, open, gap, matrix, ParasailMode.Global, default, output);

    /// <summary><c>parasail_sw</c> (<c>sw_stats</c>, <c>sw_table</c> ... with <paramref name="output"/>).</summary>
    public static ParasailResult Sw(string query, string reference, int open, int gap, ParasailMatrix matrix, ParasailOutput output = default) =>
        Run(query, reference, open, gap, matrix, ParasailMode.Local, default, output);

    /// <summary><c>parasail_sg</c> and its variants by name (<c>sg</c>, <c>sg_qx</c>, <c>sg_dx</c>, <c>sg_qb_de</c> ...).</summary>
    /// <exception cref="ArgumentException">The name is not a semi-global variant.</exception>
    public static ParasailResult Sg(string variant, string query, string reference, int open, int gap, ParasailMatrix matrix, ParasailOutput output = default) =>
        Run(query, reference, open, gap, matrix, ParasailMode.SemiGlobal, ParasailEndGaps.TryParse(variant, out var flags) ? flags : throw new ArgumentException($"'{variant}' is not a semi-global variant (sg, sg_qb, sg_qe, sg_db, sg_de, sg_qx, sg_dx, sg_qb_db, ...).", nameof(variant)), output);

    internal static ParasailResult Run(int[] s1, int[] s2, int open, int gap, ParasailMatrix mx, ParasailMode mode, ParasailEndGaps ends, ParasailOutput output)
    {
        int n1 = s1.Length, n2 = s2.Length;
        if (n1 == 0 || n2 == 0)
            throw new ArgumentException("Both sequences must be non-empty.");
        if (open < 0 || gap < 0)
            throw new ArgumentException("The gap penalties must not be negative.");
        bool stats = output.Statistics, table = output.Table, rowcol = output.RowCol;
        bool local = mode == ParasailMode.Local, sg = mode == ParasailMode.SemiGlobal;
        if (!sg)
            ends = default;
        // the shorter sequence is the inner dimension (O(min(m, n)) cells): the outer loop runs over the longer one; a table needs the original orientation
        bool transposed = n1 < n2 && !table;
        int rows = transposed ? n2 : n1, cols = transposed ? n1 : n2;
        bool outerBeg = transposed ? ends.S2Begin : ends.S1Begin, innerBeg = transposed ? ends.S1Begin : ends.S2Begin;

        var h = new int[cols + 1];
        var f = new int[cols + 1];
        int[] hm = stats ? new int[cols + 1] : [], hs = stats ? new int[cols + 1] : [], hl = stats ? new int[cols + 1] : [];
        int[] fm = stats ? new int[cols + 1] : [], fs = stats ? new int[cols + 1] : [], fl = stats ? new int[cols + 1] : [];
        int[,]? tScore = table ? new int[n1, n2] : null;
        int[,]? tM = table && stats ? new int[n1, n2] : null, tS = table && stats ? new int[n1, n2] : null, tL = table && stats ? new int[n1, n2] : null;
        int[]? rowScore = rowcol ? new int[n2] : null, colScore = rowcol ? new int[n1] : null;
        int[]? rowM = rowcol && stats ? new int[n2] : null, colM = rowcol && stats ? new int[n1] : null;
        int[]? rowS = rowcol && stats ? new int[n2] : null, colS = rowcol && stats ? new int[n1] : null;
        int[]? rowL = rowcol && stats ? new int[n2] : null, colL = rowcol && stats ? new int[n1] : null;

        h[0] = 0;
        f[0] = NegInf;
        for (int c = 1; c <= cols; c++)
        {
            h[c] = local || innerBeg ? 0 : -open - (c - 1) * gap;
            f[c] = NegInf;
        }

        // the best end cell: score, original coordinates (i, j), statistics, and the rank that orders ties among the candidates
        int bestScore = NegInf, bestI = 0, bestJ = 0, bestM = 0, bestS = 0, bestL = 0, bestRank = int.MaxValue;
        long bestKey = long.MaxValue;
        bool bestSet = false;
        bool s1End = sg && ends.S1End, s2End = sg && ends.S2End;

        void Consider(int score, int i, int j, int m, int s, int l, int rank, long key)
        {
            bool better;
            if (!bestSet)
                better = true;
            else if (score != bestScore)
                better = score > bestScore;
            else
                better = rank < bestRank || (rank == bestRank && key < bestKey);
            if (!better)
                return;
            bestSet = true;
            (bestScore, bestI, bestJ, bestM, bestS, bestL, bestRank, bestKey) = (score, i, j, m, s, l, rank, key);
        }

        for (int r = 1; r <= rows; r++)
        {
            int nh = h[0], nhm = stats ? hm[0] : 0, nhs = stats ? hs[0] : 0, nhl = stats ? hl[0] : 0;
            int wh = local || outerBeg ? 0 : -open - (r - 1) * gap;
            int whm = 0, whs = 0, whl = 0;
            int e = NegInf, em = 0, es = 0, el = 0;
            h[0] = wh;
            if (stats)
            {
                hm[0] = 0;
                hs[0] = 0;
                hl[0] = 0;
            }

            for (int c = 1; c <= cols; c++)
            {
                int i = transposed ? c - 1 : r - 1, j = transposed ? r - 1 : c - 1;
                int nwh = nh, nwm = nhm, nws = nhs, nwl = nhl;
                nh = h[c];
                if (stats)
                {
                    nhm = hm[c];
                    nhs = hs[c];
                    nhl = hl[c];
                }

                int fOpen = nh - open, fExt = f[c] - gap;
                f[c] = Math.Max(fOpen, fExt);
                int eOpen = wh - open, eExt = e - gap;
                e = Math.Max(eOpen, eExt);
                int sub = mx.Score(s1[i], s2[j]);
                int hDiag = nwh + sub;
                wh = local ? Math.Max(hDiag, 0) : hDiag;
                wh = Math.Max(wh, e);
                wh = Math.Max(wh, f[c]);
                if (stats)
                {
                    if (fOpen > fExt)
                    {
                        fm[c] = nhm;
                        fs[c] = nhs;
                        fl[c] = nhl;
                    }

                    fl[c] += 1;
                    if (eOpen > eExt)
                    {
                        em = whm;
                        es = whs;
                        el = whl;
                    }

                    el += 1;
                    // tie order of the original loop: the diagonal, then the vertical gap (a gap in the reference), then the horizontal one; transposed, the roles of E and F change places
                    if (wh == hDiag)
                    {
                        whm = nwm + (s1[i] == s2[j] ? 1 : 0);
                        whs = nws + (sub > 0 ? 1 : 0);
                        whl = nwl + 1;
                    }
                    else if (wh == (transposed ? e : f[c]))
                    {
                        (whm, whs, whl) = transposed ? (em, es, el) : (fm[c], fs[c], fl[c]);
                    }
                    else
                    {
                        (whm, whs, whl) = transposed ? (fm[c], fs[c], fl[c]) : (em, es, el);
                    }

                    if (local && wh <= 0)
                    {
                        wh = 0;
                        whm = whs = whl = 0;
                    }

                    hm[c] = whm;
                    hs[c] = whs;
                    hl[c] = whl;
                }
                else if (local && wh <= 0)
                {
                    wh = 0;
                }

                h[c] = wh;
                if (table)
                {
                    tScore![i, j] = wh;
                    if (stats)
                    {
                        tM![i, j] = whm;
                        tS![i, j] = whs;
                        tL![i, j] = whl;
                    }
                }

                if (local)
                {
                    Consider(wh, i, j, whm, whs, whl, 0, (long)j * n1 + i);
                }
                else if (s1End || s2End)
                {
                    bool lastColumn = j == n2 - 1, lastRow = i == n1 - 1;
                    if ((s1End && lastColumn) || (s2End && lastRow))
                    {
                        // both ends free: a cell of the last row (not the corner) ranks before the cells of the last column; one end free: the cells of that line in order
                        int rank = 1;
                        long key = s1End ? i : j;
                        if (s1End && s2End && lastRow && !lastColumn)
                        {
                            rank = 0;
                            key = j;
                        }

                        Consider(wh, i, j, whm, whs, whl, rank, key);
                    }
                }
            }

            if (rowcol)
            {
                // WH is the last cell of the outer row: the last column when the loop runs over the query, else the last row
                if (!transposed)
                {
                    colScore![r - 1] = wh;
                    if (stats)
                    {
                        colM![r - 1] = whm;
                        colS![r - 1] = whs;
                        colL![r - 1] = whl;
                    }
                }
                else
                {
                    rowScore![r - 1] = wh;
                    if (stats)
                    {
                        rowM![r - 1] = whm;
                        rowS![r - 1] = whs;
                        rowL![r - 1] = whl;
                    }
                }
            }
        }

        if (rowcol)
        {
            for (int c = 1; c <= cols; c++)
            {
                if (!transposed)
                {
                    rowScore![c - 1] = h[c];
                    if (stats)
                    {
                        rowM![c - 1] = hm[c];
                        rowS![c - 1] = hs[c];
                        rowL![c - 1] = hl[c];
                    }
                }
                else
                {
                    colScore![c - 1] = h[c];
                    if (stats)
                    {
                        colM![c - 1] = hm[c];
                        colS![c - 1] = hs[c];
                        colL![c - 1] = hl[c];
                    }
                }
            }
        }

        int score, endQuery, endReference, matches, similar, length;
        if (local || (sg && (s1End || s2End)))
        {
            (score, endQuery, endReference, matches, similar, length) = bestSet ? (bestScore, bestI, bestJ, bestM, bestS, bestL) : (0, 0, 0, 0, 0, 0);
        }
        else
        {
            score = h[cols];
            endQuery = n1 - 1;
            endReference = n2 - 1;
            (matches, similar, length) = stats ? (hm[cols], hs[cols], hl[cols]) : (0, 0, 0);
        }

        return new ParasailResult(
            score, endQuery, endReference,
            stats ? matches : null, stats ? similar : null, stats ? length : null,
            table ? new ParasailTables(tScore!, tM, tS, tL) : null,
            rowcol ? new ParasailRowCol(rowScore!, colScore!, rowM, colM, rowS, colS, rowL, colL) : null);
    }
}

/// <summary>
/// A query prepared for many references (the profile of parasail's 1-query-versus-database use): the query and the scoring matrix are encoded once, <see cref="Score"/> and <see cref="ScoreBatch"/> then run the program against
/// each reference with the same results as <see cref="ParasailAligner.Run(string, string, int, int, ParasailMatrix, ParasailMode, ParasailEndGaps, ParasailOutput)"/>.
/// </summary>
public sealed class ParasailProfile
{
    private readonly int[] _query;

    /// <summary>Prepares <paramref name="query"/> for <paramref name="matrix"/>.</summary>
    /// <exception cref="ArgumentException">The query is empty or has a letter outside the matrix alphabet.</exception>
    public ParasailProfile(string query, ParasailMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(query);
        ArgumentNullException.ThrowIfNull(matrix);
        if (query.Length == 0)
            throw new ArgumentException("The query must not be empty.", nameof(query));
        Matrix = matrix;
        Query = query;
        _query = matrix.Encode(query);
    }

    /// <summary>The query.</summary>
    public string Query { get; }

    /// <summary>The scoring matrix.</summary>
    public ParasailMatrix Matrix { get; }

    /// <summary>Runs the program of <paramref name="mode"/> against one reference.</summary>
    public ParasailResult Score(string reference, int open, int gap, ParasailMode mode, ParasailEndGaps ends = default, ParasailOutput output = default)
    {
        ArgumentNullException.ThrowIfNull(reference);
        return ParasailAligner.Run(_query, Matrix.Encode(reference), open, gap, Matrix, mode, ends, output);
    }

    /// <summary>Runs the program against every reference (in order; <paramref name="parallel"/> uses all cores, the results keep the order).</summary>
    public IReadOnlyList<ParasailResult> ScoreBatch(IReadOnlyList<string> references, int open, int gap, ParasailMode mode, ParasailEndGaps ends = default, ParasailOutput output = default, bool parallel = false)
    {
        ArgumentNullException.ThrowIfNull(references);
        var results = new ParasailResult[references.Count];
        if (parallel)
            Parallel.For(0, references.Count, k => results[k] = Score(references[k], open, gap, mode, ends, output));
        else
            for (int k = 0; k < results.Length; k++)
                results[k] = Score(references[k], open, gap, mode, ends, output);
        return results;
    }
}
