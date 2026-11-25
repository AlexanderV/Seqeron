using System.Globalization;

namespace Seqeron.Genomics.Infrastructure;

public sealed partial class SubstitutionMatrix
{
    /// <summary>
    /// Parses a matrix from text (R14-8): the NCBI / EMBOSS / Biopython data-file format shared by <c>BLOSUM62</c>, <c>EBLOSUM62</c>, <c>NUC.4.4</c>, <c>EDNAFULL</c>, <c>PAM250</c> and every file of Biopython's
    /// <c>substitution_matrices/data</c>: lines starting with <c>#</c> are comments, the first other line lists the column symbols separated by white space, each following line is one row: its symbol, then one
    /// real score per column. Rows may come in the order of the columns or in another one (a row is placed by its symbol); every column symbol needs a row; symbols are upper-cased; scores may be
    /// integers or decimals (at most 6 decimals: a fractional matrix, <see cref="IsFractional"/>).
    /// </summary>
    /// <param name="text">The matrix text.</param>
    /// <param name="name">The name of the matrix (default <c>CUSTOM</c>); a <c>#</c> comment does not name it.</param>
    /// <exception cref="ArgumentException">An empty text, a symbol that is not one ASCII character, a duplicated symbol, a row with the wrong number of scores, a missing row or a score that is not a number.</exception>
    public static SubstitutionMatrix Parse(string text, string name = "CUSTOM")
    {
        ArgumentNullException.ThrowIfNull(text);
        var lines = text.Replace("\r\n", "\n").Replace('\r', '\n').Split('\n').Where(l => l.Trim().Length > 0 && !l.TrimStart().StartsWith('#')).ToList();
        if (lines.Count < 2)
            throw new ArgumentException("A matrix needs a header line of symbols and at least one row.", nameof(text));
        var header = lines[0].Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        foreach (string symbol in header)
        {
            if (symbol.Length != 1 || symbol[0] >= 128)
                throw new ArgumentException($"Matrix symbol '{symbol}' is not a single ASCII character (a codon matrix is CodonSubstitutionMatrix).", nameof(text));
        }

        string alphabet = string.Concat(header).ToUpperInvariant();
        int n = alphabet.Length;
        if (lines.Count - 1 != n)
            throw new ArgumentException($"The header lists {n} symbols but the text has {lines.Count - 1} rows.", nameof(text));
        var values = new decimal[n, n];
        var seen = new bool[n];
        for (int r = 1; r < lines.Count; r++)
        {
            var tokens = lines[r].Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            if (tokens[0].Length != 1)
                throw new ArgumentException($"Row {r}: the first token '{tokens[0]}' is not a single symbol.", nameof(text));
            int row = alphabet.IndexOf(char.ToUpperInvariant(tokens[0][0]), StringComparison.Ordinal);
            if (row < 0 || seen[row])
                throw new ArgumentException($"Row {r}: symbol '{tokens[0]}' is not a (new) column symbol.", nameof(text));
            seen[row] = true;
            if (tokens.Length - 1 != n)
                throw new ArgumentException($"Row '{tokens[0]}': {tokens.Length - 1} scores for {n} columns.", nameof(text));
            for (int c = 0; c < n; c++)
            {
                if (!decimal.TryParse(tokens[c + 1], NumberStyles.Float, CultureInfo.InvariantCulture, out decimal v))
                    throw new ArgumentException($"Row '{tokens[0]}', column '{header[c]}': '{tokens[c + 1]}' is not a number.", nameof(text));
                values[row, c] = v;
            }
        }

        bool allIntegral = values.Cast<decimal>().All(v => v == decimal.Truncate(v));
        if (allIntegral)
        {
            var ints = new int[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    ints[i, j] = checked((int)values[i, j]);
            }

            return new SubstitutionMatrix(name, alphabet, ints);
        }

        return new SubstitutionMatrix(name, alphabet, values);
    }

    /// <summary>The text of the matrix in the shared data-file format (<see cref="Parse"/> reads it back: the same cells).</summary>
    public string ToText()
    {
        var sb = new System.Text.StringBuilder();
        sb.Append("# ").Append(Name).Append('\n');
        sb.Append(' ');
        foreach (char c in Alphabet)
            sb.Append("  ").Append(c);
        sb.Append('\n');
        for (int i = 0; i < Alphabet.Length; i++)
        {
            sb.Append(Alphabet[i]);
            for (int j = 0; j < Alphabet.Length; j++)
                sb.Append(' ').Append(ValueAt(i, j).ToString(CultureInfo.InvariantCulture));
            sb.Append('\n');
        }

        return sb.ToString();
    }
}

/// <summary>
/// A substitution matrix over codons (R14-8): the SCHNEIDER empirical codon substitution matrix of Biopython 1.88 (Schneider, Cannarozzi &amp; Gonnet 2005, BMC Bioinformatics 6:134), 64 x 64 real scores indexed by the
/// three-letter codon (A, C, G, T; case-insensitive; U is read as T). The residue-pair matrices of <see cref="SubstitutionMatrix"/> are single-character; a codon is three characters, hence a type of its own.
/// </summary>
public sealed partial class CodonSubstitutionMatrix
{
    private readonly string[] _codons;
    private readonly decimal[,] _values;
    private readonly Dictionary<string, int> _index;

    /// <summary>The matrix name.</summary>
    public string Name { get; }

    /// <summary>The citation / description line of the data file.</summary>
    public string Description { get; }

    /// <summary>The codons in row / column order (64 for SCHNEIDER).</summary>
    public IReadOnlyList<string> Codons => _codons;

    private CodonSubstitutionMatrix(string name, string description, string[] codons, decimal[,] values)
    {
        Name = name;
        Description = description;
        _codons = codons;
        _values = values;
        _index = codons.Select((c, i) => (c, i)).ToDictionary(x => x.c, x => x.i, StringComparer.Ordinal);
    }

    private static readonly Lazy<CodonSubstitutionMatrix> LazySchneider = new(() =>
    {
        var codons = SchneiderCodons.Split(' ');
        int n = codons.Length;
        var v = new decimal[n, n];
        int k = 0;
        foreach (var token in SchneiderScores.Split(' '))
        {
            v[k / n, k % n] = decimal.Parse(token, NumberStyles.Float, CultureInfo.InvariantCulture);
            k++;
        }

        return new CodonSubstitutionMatrix("SCHNEIDER", SchneiderDescription, codons, v);
    });

    /// <summary>The SCHNEIDER codon matrix (Biopython 1.88 <c>substitution_matrices.load('SCHNEIDER')</c>, every cell equal).</summary>
    public static CodonSubstitutionMatrix Schneider => LazySchneider.Value;

    /// <summary>The score of a pair of codons; false when a codon is not three of A / C / G / T (or U).</summary>
    public bool TryGetValue(string codon1, string codon2, out decimal value)
    {
        value = 0;
        if (codon1 is null || codon2 is null || !_index.TryGetValue(Normalize(codon1), out int i) || !_index.TryGetValue(Normalize(codon2), out int j))
            return false;
        value = _values[i, j];
        return true;
    }

    /// <summary>The score of a pair of codons.</summary>
    /// <exception cref="ArgumentException">A codon is not three of A / C / G / T.</exception>
    public decimal GetValue(string codon1, string codon2) =>
        TryGetValue(codon1, codon2, out decimal v) ? v : throw new ArgumentException($"Codon pair '{codon1}','{codon2}' is not in matrix {Name}.");

    public static string Normalize(string codon) => codon.ToUpperInvariant().Replace('U', 'T');
}
