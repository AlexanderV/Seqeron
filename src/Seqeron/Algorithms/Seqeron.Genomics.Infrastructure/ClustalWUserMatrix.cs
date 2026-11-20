using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Infrastructure;

// User-defined substitution matrices of ClustalW 2.1 (audit A1-10d): the files given with -MATRIX=file / -PWMATRIX=file (and -DNAMATRIX / -PWDNAMATRIX), read by
// SubMatrix::readUserMatrix (one matrix) and SubMatrix::readMatrixSeries (a CLUSTAL_SERIES file whose "MATRIX llimit ulimit file" lines name one matrix file per
// percent-identity range). The reader is a step-by-step port (clustalw 2.1+lgpl-7, src/substitutionMatrix/SubMatrix.cpp, sources.debian.org), quirks included:
//
//   * comment lines start with '#', blank lines (only white space) are skipped; the FIRST other line is the header: its alphabetic characters are the residues, in
//     column order (any other character, e.g. the '*' of the NCBI / EMBOSS files, is not a residue);
//   * the rows follow: an optional leading label token (a token whose first character is a letter is skipped), then the values of the row. Only the lower triangle
//     (columns <= row) is read, at least as many tokens as residues must be present, exactly (residues + 1) rows are required (the NCBI / EMBOSS files have one
//     more row than letters: the '*' row, which is read and never used; a file without such a row is rejected by ClustalW as "wrong format"), reading stops at a
//     line that starts with '#' or '!';
//   * a row whose first value contains a '.' is a floating-point row: every value of that row is multiplied by 10 (so 0.5 is 5) and truncated to a short;
//     otherwise the values are truncated integers (the scale is decided per row);
//   * a letter that is not one of "ABCDEFGHIKLMNOPQRSTUVWXYZ" (ClustalW's residue codes: no J, upper case only) gets no code ("residue ... not recognised", only a warning
//     in ClustalW): its row is skipped when the matrix is built, so the stored values are then read in a shifted order - reproduced here (Cross = -1, Cells()).
/// <summary>A user-defined ClustalW 2.1 substitution matrix (<c>-MATRIX=file</c> / <c>-PWMATRIX=file</c> / <c>-DNAMATRIX</c> / <c>-PWDNAMATRIX</c>, or one entry of a CLUSTAL_SERIES file), read exactly as ClustalW's <c>SubMatrix::readUserMatrix</c> does (quirks included, see the notes above).</summary>
public sealed class ClustalWUserMatrix
{
    /// <summary><c>UserParameters::aminoAcidCodes</c>: ClustalW's residue codes (the gap '-' last). Letters outside it (J, lower case) are not residues.</summary>
    public const string ResidueCodes = "ABCDEFGHIKLMNOPQRSTUVWXYZ-";

    /// <summary>Largest number of residues of a matrix (<c>UserParameters::maxAA</c> + 1 = 25 residue codes).</summary>
    public const int MaxResidues = 25;

    private readonly short[] _values;
    private readonly int[] _cross;

    private ClustalWUserMatrix(string name, string residues, int[] cross, short[] values)
    {
        Name = name;
        Residues = residues;
        _cross = cross;
        _values = values;
    }

    /// <summary>Name of the matrix (the file name, or a label chosen by the caller).</summary>
    public string Name { get; }

    /// <summary>The residues of the header line, in column order (alphabetic characters only).</summary>
    public string Residues { get; }

    /// <summary>For each header position the index of its residue in <see cref="ResidueCodes"/> (ClustalW's <c>Xref</c>), -1 for a letter ClustalW does not know.</summary>
    public IReadOnlyList<int> Cross => _cross;

    /// <summary>The values as ClustalW stores them (<c>Matrix userMat</c>): per file row the lower-triangle cells whose column residue is known, truncated to shorts (floating-point rows already x10).</summary>
    public IReadOnlyList<short> Values => _values;

    /// <summary>Number of residues of the header that ClustalW knows.</summary>
    public int RecognisedResidues => _cross.Count(c => c >= 0);

    /// <summary>
    /// The cells <c>SubMatrix::getMatrix</c> consumes from the stored values, in its order (row <c>i</c>, columns <c>j &lt;= i</c>, both residues known): the position
    /// pairs of the header and the stored value. For a header without unknown letters this is exactly the lower triangle of the file.
    /// </summary>
    public IEnumerable<(int Row, int Column, short Value)> Cells()
    {
        int ix = 0;
        for (int i = 0; i < _cross.Length; i++)
        {
            if (_cross[i] < 0)
                continue;
            for (int j = 0; j <= i; j++)
            {
                if (_cross[j] < 0)
                    continue;
                if (ix >= _values.Length)
                    throw new FormatException($"Matrix '{Name}' holds fewer values than its residues need.");
                yield return (i, j, _values[ix++]);
            }
        }
    }

    /// <summary>
    /// The canonical symmetric <see cref="SubstitutionMatrix"/> of the file (alphabet = <see cref="Residues"/>; cell (i, j) = cell (j, i) = the stored lower-triangle
    /// value, i.e. the integer ClustalW reads: floating-point rows x10, truncated; the upper triangle of the file is ignored, as in ClustalW).
    /// </summary>
    /// <exception cref="ArgumentException">A header letter is unknown to ClustalW or occurs twice.</exception>
    public SubstitutionMatrix ToSubstitutionMatrix()
    {
        int k = Residues.Length;
        if (_cross.Any(c => c < 0))
            throw new ArgumentException($"Matrix '{Name}' has residues ClustalW does not know; it has no canonical form.");
        var m = new int[k, k];
        foreach (var (row, col, value) in Cells())
        {
            m[row, col] = value;
            m[col, row] = value;
        }

        return new SubstitutionMatrix(Name, Residues, m);
    }

    /// <summary>
    /// The matrix as ClustalW would read it from a file holding <paramref name="matrix"/> (the symbols of the alphabet that are ClustalW residue codes, in alphabet
    /// order, e.g. without the '*'; lower triangle, integers). Lets a canonical <see cref="SubstitutionMatrix"/> be used as <c>-MATRIX</c> / <c>-PWMATRIX</c>.
    /// </summary>
    /// <exception cref="ArgumentException">No more than two of the symbols are ClustalW residues, or more than 25.</exception>
    public static ClustalWUserMatrix FromSubstitutionMatrix(SubstitutionMatrix matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        var letters = new StringBuilder();
        foreach (char c in matrix.Alphabet)
        {
            if (c != '-' && ResidueCodes.Contains(c))
                letters.Append(c);
        }

        string residues = letters.ToString();
        if (residues.Length < 3)
            throw new ArgumentException("A ClustalW matrix needs at least three residues.", nameof(matrix));
        if (residues.Length > MaxResidues)
            throw new ArgumentException($"A ClustalW matrix has at most {MaxResidues} residues.", nameof(matrix));

        var cross = residues.Select(c => ResidueCodes.IndexOf(c)).ToArray();
        var values = new List<short>();
        for (int i = 0; i < residues.Length; i++)
            for (int j = 0; j <= i; j++)
                values.Add(checked((short)matrix.GetScore(residues[i], residues[j])));
        return new ClustalWUserMatrix(matrix.Name, residues, cross, values.ToArray());
    }

    /// <summary>Reads a matrix file (<c>SubMatrix::readUserMatrix</c>; the text of the file, '\n' separated).</summary>
    /// <param name="text">The content of the matrix file.</param>
    /// <param name="name">Name of the matrix (default "USER").</param>
    /// <exception cref="FormatException">ClustalW reports "wrong format in matrix" / "in ... - single matrix expected" / "too many entries" (or would crash).</exception>
    public static ClustalWUserMatrix Parse(string text, string name = "USER")
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(name);
        var lines = SplitLines(text, name);
        int pos = 0;
        string? header = null;
        while (pos < lines.Count)
        {
            string line = lines[pos++];
            if (IsCommentLine(line))
                continue;
            if (line.StartsWith("CLUSTAL_SERIES", StringComparison.Ordinal))
                throw new FormatException($"In {name} - single matrix expected.");
            header = line;
            break;
        }

        int k = header is null ? 0 : header.Count(IsAsciiAlpha);
        if (k == 0)
            throw new FormatException($"Wrong format in matrix {name}.");
        if (k > MaxResidues)
            throw new FormatException($"Matrix {name} has {k} residues; ClustalW's amino-acid codes allow {MaxResidues}.");

        var residues = new string(header!.Where(IsAsciiAlpha).ToArray());
        var cross = new int[k];
        int recognised = 0;
        for (int i = 0; i < k; i++)
        {
            int j = ResidueCodes.IndexOf(residues[i]);
            cross[i] = j;
            if (j >= 0)
                recognised++;
        }
        var values = new List<short>();
        int ix = 0;
        while (pos < lines.Count)
        {
            string line = lines[pos++];
            if (line.Length > 0 && line[0] == '\n')
                continue;
            if (line.Length > 0 && (line[0] == '#' || line[0] == '!'))
                break;

            var args = Tokenise(line, k + 2);
            if (args.Count < recognised)
                throw new FormatException($"Wrong format in matrix {name}.");

            int farg = args[0].Length > 0 && IsAsciiAlpha(args[0][0]) ? 1 : 0;
            double scale = args[farg].Contains('.') ? 10.0 : 1.0;
            for (int i = 0; i <= ix; i++)
            {
                if (i >= k || cross[i] < 0)
                    continue;   // Xref is -1 beyond the header (the '*' column) and for a letter ClustalW does not know
                if (i + farg >= args.Count)
                    throw new FormatException($"Wrong format in matrix {name}: row {ix + 1} has too few values (ClustalW reads past the row).");
                double f = Atof(args[i + farg]);
                values.Add(unchecked((short)(int)(f * scale)));
            }

            ix++;
        }

        if (ix != k + 1)
            throw new FormatException($"Wrong format in matrix {name}: {ix} rows for {k} residues (ClustalW needs one row more than letters, the '*' row of the NCBI / EMBOSS files).");

        return new ClustalWUserMatrix(name, residues, cross, values.ToArray());
    }

    /// <summary>Reads a matrix file from disk (<see cref="Parse"/>; the name is the file name).</summary>
    public static ClustalWUserMatrix ReadFile(string path)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        return Parse(File.ReadAllText(path), path);
    }

    // ---- shared helpers of the file readers -------------------------------------------------------------------------------------------------------

    internal static bool IsAsciiAlpha(char c) => c is (>= 'A' and <= 'Z') or (>= 'a' and <= 'z');

    /// <summary><c>SubMatrix::commentline</c>: a line that starts with '#' or holds only white space.</summary>
    internal static bool IsCommentLine(string line)
    {
        if (line.Length > 0 && line[0] == '#')
            return true;
        foreach (char c in line)
        {
            if (c == '\n')
                break;
            if (!(c is ' ' or '\t' or '\v' or '\f' or '\r'))
                return false;
        }

        return true;
    }

    /// <summary>The lines of a file as <c>fgets</c> returns them (each keeps its '\n'); ClustalW reads at most 1023 characters per line, longer lines are rejected here.</summary>
    internal static List<string> SplitLines(string text, string name)
    {
        var lines = new List<string>();
        int start = 0;
        while (start < text.Length)
        {
            int nl = text.IndexOf('\n', start);
            int end = nl < 0 ? text.Length : nl + 1;
            if (end - start > 1023)
                throw new FormatException($"A line of {name} is longer than the 1023 characters ClustalW reads per line.");
            lines.Add(text.Substring(start, end - start));
            start = end;
        }

        return lines;
    }

    /// <summary><c>SubMatrix::getArgs</c>: <c>strtok</c> on " \t\n", at most <paramref name="max"/> tokens.</summary>
    internal static List<string> Tokenise(string line, int max)
    {
        var tokens = new List<string>();
        int i = 0;
        while (tokens.Count < max)
        {
            while (i < line.Length && line[i] is ' ' or '\t' or '\n')
                i++;
            if (i >= line.Length)
                break;
            int s = i;
            while (i < line.Length && !(line[i] is ' ' or '\t' or '\n'))
                i++;
            tokens.Add(line.Substring(s, i - s));
        }

        return tokens;
    }

    /// <summary>C <c>atof</c>: the longest numeric prefix of the token (0 when there is none).</summary>
    internal static double Atof(string token)
    {
        int i = 0;
        int n = token.Length;
        if (i < n && (token[i] is '+' or '-'))
            i++;
        int digitsStart = i;
        while (i < n && token[i] is >= '0' and <= '9')
            i++;
        int intDigits = i - digitsStart;
        int fracDigits = 0;
        if (i < n && token[i] == '.')
        {
            int f = i + 1;
            while (f < n && token[f] is >= '0' and <= '9')
                f++;
            fracDigits = f - (i + 1);
            if (intDigits + fracDigits > 0)
                i = f;
        }

        if (intDigits + fracDigits == 0)
            return 0.0;
        if (i < n && (token[i] is 'e' or 'E'))
        {
            int e = i + 1;
            if (e < n && (token[e] is '+' or '-'))
                e++;
            int expStart = e;
            while (e < n && token[e] is >= '0' and <= '9')
                e++;
            if (e > expStart)
                i = e;
        }

        return double.Parse(token.AsSpan(0, i), NumberStyles.Float, CultureInfo.InvariantCulture);
    }
}

/// <summary>One line <c>MATRIX llimit ulimit file</c> of a ClustalW matrix series: the matrix used for profile pairs whose percent identity lies in [llimit, ulimit].</summary>
/// <param name="LowerLimit">Lower percent-identity limit (0-100).</param>
/// <param name="UpperLimit">Upper limit (greater than <paramref name="LowerLimit"/>, at most 100).</param>
/// <param name="Matrix">The matrix.</param>
public sealed record ClustalWMatrixSeriesEntry(int LowerLimit, int UpperLimit, ClustalWUserMatrix Matrix);

/// <summary>What the file of <c>clustalw -MATRIX=file</c> holds: one matrix, or a series of matrices selected by the percent identity of the profile pair.</summary>
/// <param name="Matrix">The single matrix (<c>null</c> for a series).</param>
/// <param name="Series">The series entries (<c>null</c> for a single matrix); at most 10 (ClustalW's <c>MAXMAT</c>).</param>
public sealed record ClustalWMatrixFile(ClustalWUserMatrix? Matrix, IReadOnlyList<ClustalWMatrixSeriesEntry>? Series)
{
    /// <summary>ClustalW's <c>MAXMAT</c>: the most matrices of a series that are read.</summary>
    public const int MaxSeriesMatrices = 10;

    /// <summary>
    /// <c>SubMatrix::readMatrixSeries</c>: the first non-comment line decides: a line starting with <c>CLUSTAL_SERIES</c> opens a series (every following line starting
    /// with <c>MATRIX</c> reads "llimit ulimit file", the other lines are ignored), anything else is a single matrix (<see cref="ClustalWUserMatrix.Parse"/>).
    /// </summary>
    /// <param name="text">Content of the file given to <c>-MATRIX</c>.</param>
    /// <param name="readMatrixFile">Returns the content of a matrix file named in the series (ClustalW opens the name relative to its working directory); required for a series.</param>
    /// <param name="name">Name of the file (for messages and the single matrix).</param>
    /// <exception cref="FormatException">The file is malformed in a way ClustalW rejects.</exception>
    public static ClustalWMatrixFile Parse(string text, Func<string, string>? readMatrixFile = null, string name = "USER")
    {
        ArgumentNullException.ThrowIfNull(text);
        var lines = ClustalWUserMatrix.SplitLines(text, name);
        int pos = 0;
        bool series = false;
        while (pos < lines.Count)
        {
            string line = lines[pos++];
            if (ClustalWUserMatrix.IsCommentLine(line))
                continue;
            series = line.StartsWith("CLUSTAL_SERIES", StringComparison.Ordinal);
            break;
        }

        if (!series)
            return new ClustalWMatrixFile(ClustalWUserMatrix.Parse(text, name), null);

        if (readMatrixFile is null)
            throw new ArgumentException("A matrix series names matrix files; a reader for them is required.", nameof(readMatrixFile));

        var entries = new List<ClustalWMatrixSeriesEntry>();
        while (pos < lines.Count)
        {
            string line = lines[pos++];
            if (ClustalWUserMatrix.IsCommentLine(line) || !line.StartsWith("MATRIX", StringComparison.Ordinal))
                continue;

            if (!TryScan(line.Substring(6), out int lower, out int upper, out string file))
                throw new FormatException($"Bad format in file {name}.");
            if (lower < 0 || lower > 100 || upper < 0 || upper > 100)
                throw new FormatException($"Bad format in file {name}.");
            if (upper <= lower)
                throw new FormatException($"In file {name}: lower limit is greater than upper ({lower}-{upper}).");

            entries.Add(new ClustalWMatrixSeriesEntry(lower, upper, ClustalWUserMatrix.Parse(readMatrixFile(file), file)));
            if (entries.Count >= MaxSeriesMatrices)
                break;   // "The first 10 have been read in and will be used."
        }

        if (entries.Count == 0)
            throw new FormatException($"The series file {name} names no matrix (ClustalW would then use an uninitialised matrix).");
        return new ClustalWMatrixFile(null, entries);
    }

    /// <summary>Reads the file given to <c>-MATRIX</c>; the matrix files of a series are opened relative to the working directory, as ClustalW does.</summary>
    public static ClustalWMatrixFile ReadFile(string path)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        return Parse(File.ReadAllText(path), File.ReadAllText, path);
    }

    // sscanf(s, "%d %d %s"): two decimal integers and a white-space-free word
    private static bool TryScan(string s, out int a, out int b, out string word)
    {
        a = b = 0;
        word = "";
        int i = 0;
        if (!ScanInt(s, ref i, out a) || !ScanInt(s, ref i, out b))
            return false;
        while (i < s.Length && char.IsWhiteSpace(s[i]))
            i++;
        int st = i;
        while (i < s.Length && !char.IsWhiteSpace(s[i]))
            i++;
        if (i == st)
            return false;
        word = s.Substring(st, i - st);
        return true;
    }

    private static bool ScanInt(string s, ref int i, out int value)
    {
        value = 0;
        while (i < s.Length && char.IsWhiteSpace(s[i]))
            i++;
        int st = i;
        if (i < s.Length && (s[i] is '+' or '-'))
            i++;
        int digits = i;
        while (i < s.Length && s[i] is >= '0' and <= '9')
            i++;
        if (i == digits)
            return false;
        return int.TryParse(s.AsSpan(st, i - st), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out value);
    }
}
