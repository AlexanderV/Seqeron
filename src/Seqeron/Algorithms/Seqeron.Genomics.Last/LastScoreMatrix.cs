// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last: ScoreMatrix.cc, LambdaCalculator.cc, mcf_substitution_matrix_stats.cc, cbrc_linalg.cc), Copyright 2008 - 2019 Martin C. Frith, 2015 Yutaro Konta; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Last;

/// <summary>The substitution score matrix of LAST (<c>ScoreMatrix.cc</c>): the matrix file, the 66 x 66 tables by letter code, the scores of the ambiguity codes.</summary>
internal sealed class LastScoreMatrix
{
    /// <summary><c>INF</c> of <c>ScoreMatrixRow.hh</c>: the (negated) score of the delimiter symbol.</summary>
    public const int Inf = 1073692800;

    private static readonly string[] NtAmbiguities = ["MAC", "SCG", "KGT", "WTA", "RAG", "YCT", "BCGT", "DAGT", "HACT", "VACG", "NACGT"];

    public string RowSymbols { get; private set; } = "";

    public string ColSymbols { get; private set; } = "";

    public List<List<int>> Cells { get; private set; } = [];

    public int[][] CaseSensitive { get; } = NewTable();

    public int[][] CaseInsensitive { get; } = NewTable();

    public int MinScore { get; private set; }

    public int MaxScore { get; private set; }

    private static int[][] NewTable() => Enumerable.Range(0, LastAlphabet.Capacity).Select(_ => new int[LastAlphabet.Capacity]).ToArray();

    /// <summary><c>ScoreMatrix::setMatchMismatch</c>.</summary>
    public void SetMatchMismatch(int matchScore, int mismatchCost, string symbols)
    {
        RowSymbols = symbols;
        ColSymbols = symbols;
        Cells = [];
        for (int i = 0; i < symbols.Length; i++)
        {
            var row = Enumerable.Repeat(-mismatchCost, symbols.Length).ToList();
            row[i] = matchScore;
            Cells.Add(row);
        }
    }

    /// <summary><c>operator&gt;&gt;</c>: the header line of column symbols, then one line per row (symbol, scores).</summary>
    public void FromString(string text)
    {
        var colSymbols = new StringBuilder();
        var rowSymbols = new StringBuilder();
        var cells = new List<List<int>>();
        bool failed = false;
        foreach (string line in text.Split('\n'))
        {
            string[] words = line.Split([' ', '\t', '\r', '\f', '\v'], StringSplitOptions.RemoveEmptyEntries);
            if (words.Length == 0)
                continue;
            if (colSymbols.Length == 0)
            {
                if (words[0][0] == '#')
                    continue;
                if (words[0].Length == 3)
                    throw new LastProgramException("codon score matrices are not ported");
                foreach (string w in words)
                {
                    if (w.Length != 1)
                        failed = true;
                    colSymbols.Append(w);
                }
            }
            else
            {
                int numberOfColumns = colSymbols.Length;
                var row = new List<int>();
                bool ok = true;
                for (int i = 0; i < numberOfColumns; i++)
                {
                    if (i + 1 < words.Length && int.TryParse(words[i + 1], NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int score))
                        row.Add(score);
                    else
                    {
                        ok = false;
                        row.Add(0);
                    }
                }

                if (words[0].Length == 1 && ok)
                {
                    rowSymbols.Append(words[0]);
                    cells.Add(row);
                    if (words.Length > numberOfColumns + 1)
                        throw new LastProgramException("letter frequencies in a score matrix are not ported");
                }
                else
                {
                    // a line of column letter frequencies (not ported), or the end of the matrix
                    throw new LastProgramException("letter frequencies in a score matrix are not ported");
                }
            }
        }

        if (failed || cells.Count == 0)
            throw new LastProgramException("can't read the score matrix");
        RowSymbols = rowSymbols.ToString();
        ColSymbols = colSymbols.ToString();
        Cells = cells;
    }

    /// <summary><c>ScoreMatrix::init</c>.</summary>
    public void Init(byte[] symbolToIndex)
    {
        RowSymbols = RowSymbols.ToUpperInvariant();
        ColSymbols = ColSymbols.ToUpperInvariant();
        MinScore = MaxScore = Cells[0][0];
        foreach (var row in Cells)
        {
            foreach (int s in row)
            {
                MinScore = Math.Min(MinScore, s);
                MaxScore = Math.Max(MaxScore, s);
            }
        }

        for (int i = 0; i < LastAlphabet.Capacity; i++)
        {
            Array.Fill(CaseSensitive[i], MinScore);
            Array.Fill(CaseInsensitive[i], MinScore);
        }

        for (int i = 0; i < Cells.Count; i++)
        {
            for (int j = 0; j < Cells[i].Count; j++)
            {
                (int iu, int il) = UpperAndLower(RowSymbols[i], symbolToIndex);
                (int ju, int jl) = UpperAndLower(ColSymbols[j], symbolToIndex);
                int c = Cells[i][j];
                CaseSensitive[iu][jl] = Math.Min(c, 0);
                CaseSensitive[il][ju] = Math.Min(c, 0);
                CaseSensitive[il][jl] = Math.Min(c, 0);
                CaseSensitive[iu][ju] = c;
                CaseInsensitive[iu][ju] = c;
                CaseInsensitive[iu][jl] = c;
                CaseInsensitive[il][ju] = c;
                CaseInsensitive[il][jl] = c;
            }
        }

        int d = symbolToIndex[' '];
        for (int i = 0; i < LastAlphabet.Capacity; i++)
        {
            CaseSensitive[d][i] = -Inf;
            CaseSensitive[i][d] = -Inf;
            CaseInsensitive[d][i] = -Inf;
            CaseInsensitive[i][d] = -Inf;
        }
    }

    private static (int Upper, int Lower) UpperAndLower(char symbol, byte[] symbolToIndex)
    {
        int upper = symbolToIndex[symbol & 0xFF];
        int lower = symbolToIndex[char.ToLowerInvariant(symbol) & 0xFF];
        if (upper >= LastAlphabet.Capacity || lower >= LastAlphabet.Capacity)
            throw new LastProgramException("bad letter in score matrix: " + symbol);
        return (upper, lower);
    }

    /// <summary><c>ScoreMatrix::writeCommented</c>.</summary>
    public string WriteCommented()
    {
        var sb = new StringBuilder();
        int colWidth = Cells[0].Count < 20 ? 3 : 2;
        sb.Append("# ").Append(' ', 1);
        for (int i = 0; i < Cells[0].Count; i++)
            sb.Append(' ').Append(' ', colWidth - 1).Append(ColSymbols[i]);
        sb.Append('\n');
        for (int i = 0; i < Cells.Count; i++)
        {
            sb.Append("# ").Append(RowSymbols[i]);
            foreach (int c in Cells[i])
                sb.Append(' ').Append(c.ToString(CultureInfo.InvariantCulture).PadLeft(colWidth));
            sb.Append('\n');
        }

        return sb.ToString();
    }

    /// <summary><c>ScoreMatrix::addAmbiguousScores</c> for DNA (the codes M S K W R Y B D H V, and N with the options of <c>-X</c>).</summary>
    public void AddAmbiguousScores(bool isFullyAmbiguousRow, bool isFullyAmbiguousCol, byte[] symbolToIndex, double scale, double[] rowSymbolProbs, double[] colSymbolProbs)
    {
        int[][] fast = CaseInsensitive;
        int n = NtAmbiguities.Length;
        int numberOfAmbiguousRows = n - 1 + (isFullyAmbiguousRow ? 1 : 0);
        int numberOfAmbiguousCols = n - 1 + (isFullyAmbiguousCol ? 1 : 0);
        for (int k = 0; k < numberOfAmbiguousCols; k++)
        {
            char ambiguousSymbol = NtAmbiguities[k][0];
            if (ColSymbols.Contains(ambiguousSymbol))
                continue;
            ColSymbols = string.Concat(ColSymbols, ambiguousSymbol);
            for (int i = 0; i < Cells.Count; i++)
            {
                string rSymbols = AmbiguityList(NtAmbiguities, numberOfAmbiguousRows, RowSymbols[i]);
                string cSymbols = NtAmbiguities[k][1..];
                Cells[i].Add(JointScore(symbolToIndex, fast, scale, rowSymbolProbs, colSymbolProbs, rSymbols, cSymbols));
            }
        }

        for (int k = 0; k < numberOfAmbiguousRows; k++)
        {
            char ambiguousSymbol = NtAmbiguities[k][0];
            if (RowSymbols.Contains(ambiguousSymbol))
                continue;
            RowSymbols = string.Concat(RowSymbols, ambiguousSymbol);
            var row = new List<int>();
            for (int j = 0; j < Cells[0].Count; j++)
            {
                string rSymbols = NtAmbiguities[k][1..];
                string cSymbols = AmbiguityList(NtAmbiguities, numberOfAmbiguousCols, ColSymbols[j]);
                row.Add(JointScore(symbolToIndex, fast, scale, rowSymbolProbs, colSymbolProbs, rSymbols, cSymbols));
            }

            Cells.Add(row);
        }
    }

    private static string AmbiguityList(string[] ambiguities, int count, char symbol)
    {
        for (int i = 0; i < count; i++)
        {
            if (ambiguities[i][0] == symbol)
                return ambiguities[i][1..];
        }

        return symbol.ToString();
    }

    private static double SymbolProbSum(bool isOneSymbol, byte[] symbolToIndex, string symbols, double[] probs)
    {
        if (isOneSymbol)
            return 1;
        double p = 0;
        foreach (char c in symbols)
            p += probs[symbolToIndex[c]];
        return p;
    }

    private static int JointScore(byte[] symbolToIndex, int[][] fastMatrix, double scale, double[] rowSymbolProbs, double[] colSymbolProbs, string rSymbols, string cSymbols)
    {
        bool isOneRowSymbol = rSymbols.Length == 1;
        bool isOneColSymbol = cSymbols.Length == 1;
        double p = 0;
        foreach (char i in rSymbols)
        {
            foreach (char j in cSymbols)
            {
                int x = symbolToIndex[i];
                int y = symbolToIndex[j];
                double r = isOneRowSymbol ? 1 : rowSymbolProbs[x];
                double c = isOneColSymbol ? 1 : colSymbolProbs[y];
                p += r * c * Math.Exp(scale * fastMatrix[x][y]);
            }
        }

        double rowProbSum = SymbolProbSum(isOneRowSymbol, symbolToIndex, rSymbols, rowSymbolProbs);
        double colProbSum = SymbolProbSum(isOneColSymbol, symbolToIndex, cSymbols, colSymbolProbs);
        return NearestInt(Math.Log(p / (rowProbSum * colProbSum)) / scale);
    }

    /// <summary><c>nearestInt</c> of <c>qualityScoreUtil.hh</c>.</summary>
    public static int NearestInt(double x) => x > 0 ? (int)(0.5 + x) : -(int)(0.5 - x);
}

/// <summary><c>mcf::SubstitutionMatrixStats::calcUnbiased</c> for a matrix that is not one of the named protein matrices: lambda and the letter probabilities from <c>LambdaCalculator</c>.</summary>
internal sealed class LastMatrixStats
{
    public double Lambda { get; private set; } = -1;

    public double[] LetterProbs1 { get; private set; } = [];

    public double[] LetterProbs2 { get; private set; } = [];

    public bool IsBad => Lambda < 0;

    /// <summary><c>calcUnbiased</c>: the generator of <c>rand()</c> starts fresh in every process (glibc, seed 1).</summary>
    public void CalcUnbiased(int[][] scoreMatrix, int size)
    {
        var rng = new GlibcRandom(1);
        var mat = new double[size][];
        for (int i = 0; i < size; i++)
        {
            mat[i] = new double[size];
            for (int j = 0; j < size; j++)
                mat[i][j] = scoreMatrix[i][j];
        }

        var p1 = new double[size];
        var p2 = new double[size];
        Lambda = CalculateLambda(mat, size, p1, p2, 1000, 100, 1e-6, rng);
        if (IsBad)
            return;
        LetterProbs1 = p1;
        LetterProbs2 = p2;
    }

    /// <summary><c>flipDnaStrands</c> is not needed for the forward strand.</summary>
    private static double RoundToFewDigits(double x) => double.Parse(x.ToString("G6", CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);

    private static bool CalculateInvSum(double[][] matrix, int size, double tau, out double invSum, double[][] tmpMat, double[] tmpVec)
    {
        invSum = 0;
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                tmpMat[i][j] = Math.Exp(tau * matrix[i][j]);
        Array.Fill(tmpVec, 1.0, 0, size);
        if (!LinalgSolve(tmpMat, tmpVec, size))
            return false;
        double sum = 0;
        for (int i = 0; i < size; i++)
            sum += tmpVec[i];
        invSum = sum;
        return true;
    }

    /// <summary><c>linalgSolve</c> (Gauss-Jordan with partial pivoting); false = "singular matrix".</summary>
    private static bool LinalgSolve(double[][] m, double[] v, int s)
    {
        for (int k = 0; k < s; k++)
        {
            int iMax = k;
            for (int i = k; i < s; i++)
                if (Math.Abs(m[i][k]) > Math.Abs(m[iMax][k]))
                    iMax = i;
            if (iMax > k)
            {
                (m[k], m[iMax]) = (m[iMax], m[k]);
                (v[k], v[iMax]) = (v[iMax], v[k]);
            }

            if (m[k][k] == 0.0)
                return false;
            for (int i = 0; i < s; i++)
            {
                if (i == k)
                    continue;
                double q = m[i][k] / m[k][k];
                for (int j = k; j < s; j++)
                    m[i][j] -= q * m[k][j];
                v[i] -= q * v[k];
            }
        }

        for (int k = 0; k < s; k++)
            v[k] /= m[k][k];
        return true;
    }

    private static bool FindUpperBound(double[][] matrix, int size, out double ub)
    {
        ub = 0;
        double rMaxMin = double.MaxValue, cMaxMin = double.MaxValue;
        int lr = 0, lc = 0;
        for (int i = 0; i < size; i++)
        {
            double rMax = -double.MaxValue, rMin = double.MaxValue;
            for (int j = 0; j < size; j++)
            {
                if (matrix[i][j] > rMax)
                    rMax = matrix[i][j];
                if (matrix[i][j] < rMin)
                    rMin = matrix[i][j];
            }

            if (rMax == 0 && rMin == 0)
                lr++;
            else if (rMax <= 0 || rMin >= 0)
                return false;
            else if (rMax < rMaxMin)
                rMaxMin = rMax;
        }

        for (int j = 0; j < size; j++)
        {
            double cMax = -double.MaxValue, cMin = double.MaxValue;
            for (int i = 0; i < size; i++)
            {
                if (matrix[i][j] > cMax)
                    cMax = matrix[i][j];
                if (matrix[i][j] < cMin)
                    cMin = matrix[i][j];
            }

            if (cMax == 0 && cMin == 0)
                lc++;
            else if (cMax <= 0 || cMin >= 0)
                return false;
            else if (cMax < cMaxMin)
                cMaxMin = cMax;
        }

        if (lr == size)
            return false;
        ub = rMaxMin > cMaxMin ? 1.1 * Math.Log(1.0 * (size - lr)) / rMaxMin : 1.1 * Math.Log(1.0 * (size - lc)) / cMaxMin;
        return true;
    }

    private static bool CheckLambda(double[][] matrix, double lambda, int size, double[] letprob1, double[] letprob2, double[][] tmpMat)
    {
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                tmpMat[i][j] = Math.Exp(lambda * matrix[i][j]);
        Array.Fill(letprob2, 1.0, 0, size);
        if (!LinalgSolve(tmpMat, letprob2, size))
            throw new LastProgramException("singular matrix");
        for (int i = 0; i < size; i++)
        {
            double p = letprob2[i];
            if (p < 0 || p > 1)
                return false;
            letprob2[i] = RoundToFewDigits(p);
        }

        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                tmpMat[i][j] = Math.Exp(lambda * matrix[j][i]);
        Array.Fill(letprob1, 1.0, 0, size);
        if (!LinalgSolve(tmpMat, letprob1, size))
            throw new LastProgramException("singular matrix");
        for (int j = 0; j < size; j++)
        {
            double q = letprob1[j];
            if (q < 0 || q > 1)
                return false;
            letprob1[j] = RoundToFewDigits(q);
        }

        return true;
    }

    private static bool BinarySearch(double[][] matrix, int size, double lb, double ub, double[] letprob1, double[] letprob2, out double lambda, int maxIter, GlibcRandom rng)
    {
        lambda = 0;
        double l = 0, r = 0, lSum = 0, rSum = 0;
        int iter = 0;
        var tmpMat = new double[size][];
        for (int i = 0; i < size; i++)
            tmpMat[i] = new double[size];
        double[] tmpVec = letprob2;
        while (iter < maxIter && (l >= r || (lSum < 1.0 && rSum < 1.0) || (lSum > 1.0 && rSum > 1.0)))
        {
            l = lb + (ub - lb) * rng.Next() / 2147483647;
            r = lb + (ub - lb) * rng.Next() / 2147483647;
            if (!CalculateInvSum(matrix, size, l, out lSum, tmpMat, tmpVec) || !CalculateInvSum(matrix, size, r, out rSum, tmpMat, tmpVec))
            {
                l = 0;
                r = 0;
            }

            iter++;
        }

        if (iter >= maxIter)
            return false;
        while (lSum != 1.0 && rSum != 1.0 && (l + r) / 2.0 != l && (l + r) / 2.0 != r)
        {
            double mid = (l + r) / 2.0;
            if (!CalculateInvSum(matrix, size, mid, out double midSum, tmpMat, tmpVec))
                return false;
            if (Math.Abs(midSum) >= double.MaxValue)
                return false;
            if ((lSum < 1.0 && midSum >= 1.0) || (lSum > 1.0 && midSum <= 1.0))
            {
                r = mid;
                rSum = midSum;
            }
            else if ((rSum < 1.0 && midSum >= 1.0) || (rSum > 1.0 && midSum <= 1.0))
            {
                l = mid;
                lSum = midSum;
            }
            else
            {
                return false;
            }
        }

        if (Math.Abs(lSum - 1.0) < Math.Abs(rSum - 1.0))
        {
            if (CheckLambda(matrix, l, size, letprob1, letprob2, tmpMat))
            {
                lambda = l;
                return true;
            }

            return false;
        }

        if (CheckLambda(matrix, r, size, letprob1, letprob2, tmpMat))
        {
            lambda = r;
            return true;
        }

        return false;
    }

    private static double CalculateLambda(double[][] matrix, int size, double[] letprob1, double[] letprob2, int maxIter, int maxBoundarySearchIter, double lbRatio, GlibcRandom rng)
    {
        if (!FindUpperBound(matrix, size, out double ub))
            return -1;
        double lb = ub * lbRatio;
        double lambda = -1;
        int iter = 0;
        bool flag = false;
        while (!flag && iter < maxIter)
        {
            flag = BinarySearch(matrix, size, lb, ub, letprob1, letprob2, out double found, maxBoundarySearchIter, rng);
            if (flag)
                lambda = found;
            iter++;
        }

        return lambda;
    }
}
