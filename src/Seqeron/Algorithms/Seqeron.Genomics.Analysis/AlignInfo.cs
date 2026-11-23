#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// <c>Bio.Align.AlignInfo.PSSM</c> of Biopython 1.79 (the last release that has it; 1.88 keeps only the deprecated <c>get_column</c>): the position-specific score matrix of an alignment, one row per alignment column with the letter of the
/// axis sequence and the (weighted) count of every residue. <see cref="ToString"/> is <c>str(pssm)</c>.
/// </summary>
public sealed class Pssm
{
    private readonly List<(string Letter, List<KeyValuePair<char, double>> Counts)> _rows;

    internal Pssm(List<(string Letter, List<KeyValuePair<char, double>> Counts)> rows)
    {
        _rows = rows;
    }

    /// <summary><c>len(pssm.pssm)</c>: the number of positions.</summary>
    public int Count => _rows.Count;

    /// <summary><c>pssm[pos]</c>: the counts of the residues at position <paramref name="pos"/> (insertion order = the letters of the alignment in sorted order).</summary>
    /// <exception cref="ArgumentOutOfRangeException">Python's <c>IndexError: list index out of range</c>.</exception>
    public IReadOnlyDictionary<char, double> this[int pos]
    {
        get
        {
            var row = _rows[Position(pos)];
            return row.Counts.ToDictionary(kv => kv.Key, kv => kv.Value);
        }
    }

    /// <summary><c>pssm.get_residue(pos)</c>: the axis letter at <paramref name="pos"/>.</summary>
    public string GetResidue(int pos) => _rows[Position(pos)].Letter;

    private int Position(int pos)
    {
        int i = pos < 0 ? pos + _rows.Count : pos;
        if (i < 0 || i >= _rows.Count)
            throw new ArgumentOutOfRangeException(nameof(pos), "list index out of range");
        return i;
    }

    /// <summary><c>str(pssm)</c>: see <see cref="Render"/>.</summary>
    public override string ToString() => Render();

    /// <summary><c>str(pssm)</c>: a header of the residues (sorted), then one line per position, the axis letter and the counts with one decimal.</summary>
    /// <exception cref="InvalidOperationException">An empty matrix (Python raises <c>IndexError: list index out of range</c>).</exception>
    public string Render()
    {
        if (_rows.Count == 0)
            throw new InvalidOperationException("list index out of range");
        var residues = _rows[0].Counts.Select(kv => kv.Key).OrderBy(c => c).ToList();
        var sb = new StringBuilder(" ");
        foreach (char r in residues)
            sb.Append("   ").Append(r);
        sb.Append('\n');
        foreach (var (letter, counts) in _rows)
        {
            sb.Append(letter).Append(' ');
            foreach (char r in residues)
                sb.Append(' ').Append(CFormat.Fixed(counts.First(kv => kv.Key == r).Value, 1));
            sb.Append('\n');
        }

        return sb.ToString();
    }
}

/// <summary>
/// <c>Bio.Align.AlignInfo.SummaryInfo</c> of Biopython 1.79 (B09 WP67, F308, R14-13; the last release of Biopython that has the class: 1.88 keeps the deprecated <c>get_column</c> only): the summary information of a
/// <see cref="MultipleSeqAlignment"/>. Members: <see cref="DumbConsensus"/>, <see cref="GapConsensus"/>, <see cref="ReplacementDictionary"/>, <see cref="PosSpecificScoreMatrix"/> (a <see cref="Pssm"/>),
/// <see cref="InformationContent"/> with <see cref="IcVector"/> (the <c>ic_vector</c> attribute), <see cref="GetColumn"/> and <see cref="PrintInfoContent"/> (<c>print_info_content</c>). Every method is a port of the Python code with its errors
/// (messages of the original; <c>ValueError</c> and <c>TypeError</c> are <see cref="ArgumentException"/> / <see cref="InvalidOperationException"/>, <c>KeyError</c> is <see cref="KeyNotFoundException"/>,
/// <c>ZeroDivisionError</c> is <see cref="DivideByZeroException"/>). The residues are the letters of the rows as they are (no case folding); the sequence weights are the <c>weight</c> annotations of the rows (default 1.0).
/// </summary>
public sealed class SummaryInfo
{
    private const char Gap = '-';

    /// <summary>Creates the summary of an alignment.</summary>
    public SummaryInfo(MultipleSeqAlignment alignment)
    {
        Alignment = alignment ?? throw new ArgumentNullException(nameof(alignment));
    }

    /// <summary><c>summary.alignment</c>.</summary>
    public MultipleSeqAlignment Alignment { get; }

    /// <summary><c>summary.ic_vector</c>: the information content of every column of the last <see cref="InformationContent"/> call (empty before the first).</summary>
    public IReadOnlyList<double> IcVector { get; private set; } = [];

    private static double Weight(SeqRecord r) => r.Annotations.TryGetPropertyValue("weight", out var w) && w is not null ? w.GetValue<double>() : 1.0;

    // dumb_consensus / gap_consensus differ in the letters they count
    private string Consensus(double threshold, string ambiguous, bool requireMultiple, bool countGaps)
    {
        var consensus = new StringBuilder();
        int conLen = Alignment.GetAlignmentLength();
        for (int n = 0; n < conLen; n++)
        {
            var atoms = new List<KeyValuePair<char, int>>();
            int numAtoms = 0;
            foreach (var record in Alignment)
            {
                if (n >= record.Length)
                    continue;
                char c = record.Seq[n];
                if (!countGaps && (c == '-' || c == '.'))
                    continue;
                int at = atoms.FindIndex(a => a.Key == c);
                if (at < 0)
                    atoms.Add(new KeyValuePair<char, int>(c, 1));
                else
                    atoms[at] = new KeyValuePair<char, int>(c, atoms[at].Value + 1);
                numAtoms++;
            }

            var maxAtoms = new List<char>();
            int maxSize = 0;
            foreach (var (atom, count) in atoms)
            {
                if (count > maxSize)
                {
                    maxAtoms = [atom];
                    maxSize = count;
                }
                else if (count == maxSize)
                {
                    maxAtoms.Add(atom);
                }
            }

            if (requireMultiple && numAtoms == 1)
                consensus.Append(ambiguous);
            else if (maxAtoms.Count == 1 && (double)maxSize / numAtoms >= threshold)
                consensus.Append(maxAtoms[0]);
            else
                consensus.Append(ambiguous);
        }

        return consensus.ToString();
    }

    /// <summary>
    /// <c>summary.dumb_consensus(threshold, ambiguous, require_multiple)</c>: per column the most common residue (gaps <c>-</c> and <c>.</c> not counted) when it is unique and its share of the counted residues is at least
    /// <paramref name="threshold"/>, else <paramref name="ambiguous"/>; with <paramref name="requireMultiple"/> a column with one residue only is ambiguous.
    /// </summary>
    public string DumbConsensus(double threshold = 0.7, string ambiguous = "X", bool requireMultiple = false) => Consensus(threshold, ambiguous, requireMultiple, countGaps: false);

    /// <summary><c>summary.gap_consensus(...)</c>: as <see cref="DumbConsensus"/> but the gap symbols are counted like residues (a column of mostly gaps gives a gap).</summary>
    public string GapConsensus(double threshold = 0.7, string ambiguous = "X", bool requireMultiple = false) => Consensus(threshold, ambiguous, requireMultiple, countGaps: true);

    /// <summary>
    /// <c>summary.replacement_dictionary(skip_chars, letters)</c>: for every ordered pair of <paramref name="letters"/> the summed product of the weights of the two rows that carry them in one column (every pair of rows once, the letter of the
    /// EARLIER row first); the keys are all pairs in the order of the letters. The <c>weight</c> annotations of the rows weight the pairs.
    /// </summary>
    /// <param name="letters">The letters to count (required: Python's default <c>None</c> raises <c>TypeError</c>).</param>
    /// <param name="skipChars">Not used; anything other than null is a <see cref="ArgumentException"/> (Python's <c>ValueError</c>).</param>
    public System.Collections.Generic.OrderedDictionary<(char First, char Second), double> ReplacementDictionary(string? letters = null, string? skipChars = null)
    {
        if (skipChars is not null)
            throw new ArgumentException("argument skip_chars has been deprecated; instead, please use 'letters' to specify the characters you want to include");
        if (letters is null)
            throw new InvalidOperationException("'NoneType' object is not iterable");
        var rep = new System.Collections.Generic.OrderedDictionary<(char, char), double>();
        foreach (char a in letters)
            foreach (char b in letters)
                rep[(a, b)] = 0;
        for (int r1 = 0; r1 < Alignment.Count; r1++)
        {
            for (int r2 = r1 + 1; r2 < Alignment.Count; r2++)
            {
                string s1 = Alignment[r1].Seq, s2 = Alignment[r2].Seq;
                double w1 = Weight(Alignment[r1]), w2 = Weight(Alignment[r2]);
                int len = Math.Min(s1.Length, s2.Length);
                for (int c = 0; c < len; c++)
                {
                    if (letters.Contains(s1[c]) && letters.Contains(s2[c]))
                        rep[(s1[c], s2[c])] += w1 * w2;
                }
            }
        }

        return rep;
    }

    // _get_all_letters: the letters of all rows, sorted by code point
    private string AllLetters() => new([.. new SortedSet<char>(Alignment.SelectMany(r => r.Seq))]);

    private static string RemoveAll(string letters, string piece) => piece.Length == 0 ? letters : letters.Replace(piece, "", StringComparison.Ordinal);

    /// <summary>
    /// <c>summary.pos_specific_score_matrix(axis_seq, chars_to_ignore)</c>: per alignment column the (weighted) count of every letter of the alignment (the sorted letters minus the gap symbol and <paramref name="charsToIgnore"/>); the axis
    /// letters are <paramref name="axisSeq"/> or, when null or empty, <see cref="DumbConsensus"/>. As in Biopython the gap symbol <c>-</c> is appended to <paramref name="charsToIgnore"/> in place.
    /// </summary>
    /// <exception cref="InvalidOperationException">An alignment without letters (Python's <c>AssertionError</c>), or an axis of another length than the alignment (<c>AssertionError</c>).</exception>
    /// <exception cref="ArgumentException">A residue that is not among the counted letters ("Residue X not found", Python's <c>ValueError</c>).</exception>
    public Pssm PosSpecificScoreMatrix(string? axisSeq = null, List<string>? charsToIgnore = null)
    {
        string allLetters = AllLetters();
        if (allLetters.Length == 0)
            throw new InvalidOperationException("AssertionError");
        charsToIgnore ??= [];
        charsToIgnore.Add("-");
        foreach (string c in charsToIgnore)
            allLetters = RemoveAll(allLetters, c);
        string left;
        if (!string.IsNullOrEmpty(axisSeq))
        {
            left = axisSeq;
            if (axisSeq.Length != Alignment.GetAlignmentLength())
                throw new InvalidOperationException("AssertionError");
        }
        else
        {
            left = DumbConsensus();
        }

        var info = new List<(string, List<KeyValuePair<char, double>>)>();
        for (int n = 0; n < left.Length; n++)
        {
            var score = allLetters.Select(l => new KeyValuePair<char, double>(l, 0)).ToList();
            foreach (var record in Alignment)
            {
                if (n >= record.Length)
                    continue;
                string residue = record.Seq[n].ToString();
                if (charsToIgnore.Contains(residue))
                    continue;
                double weight = Weight(record);
                int at = score.FindIndex(kv => kv.Key == residue[0]);
                if (at < 0)
                    throw new ArgumentException($"Residue {residue} not found");
                score[at] = new KeyValuePair<char, double>(residue[0], score[at].Value + weight);
            }

            info.Add((left[n].ToString(), score));
        }

        return new Pssm(info);
    }

    /// <summary>
    /// <c>summary.information_content(start, end, e_freq_table, log_base, chars_to_ignore, pseudo_count)</c>: the sum over the columns <paramref name="start"/> .. <paramref name="end"/> - 1 of
    /// <c>sum over the letters of f * log(f / e) / log(log_base)</c> with f the weighted frequency of the letter in the column (the gap symbol and <paramref name="charsToIgnore"/> not counted; the gap symbol, when counted, only
    /// enters the denominator) and e the expected frequency of <paramref name="eFreqTable"/>; with <paramref name="pseudoCount"/> k f is (count + e * k) / (total + k). <see cref="IcVector"/> is set to the per-column values.
    /// Biopython needs an expected-frequency table: without one the call fails with a <c>TypeError</c> as soon as a column has a letter (reproduced as <see cref="InvalidOperationException"/>).
    /// </summary>
    /// <param name="start">First column (0-based).</param>
    /// <param name="end">End column (exclusive; default the length of the first row).</param>
    /// <param name="eFreqTable">Expected frequency of every letter (not of the gap symbol).</param>
    /// <param name="logBase">Base of the logarithm (default 2: bits).</param>
    /// <param name="charsToIgnore">Letters that are not counted.</param>
    /// <param name="pseudoCount">Pseudo count k (default 0).</param>
    /// <exception cref="ArgumentException">Start or end outside the first row, a negative pseudo count, a letter missing in the table, a residue that is not among the counted letters.</exception>
    public double InformationContent(int start = 0, int? end = null, IReadOnlyDictionary<char, double>? eFreqTable = null, double logBase = 2, IEnumerable<string>? charsToIgnore = null, double pseudoCount = 0)
    {
        if (Alignment.Count == 0)
            throw new ArgumentOutOfRangeException(nameof(end), "list index out of range");
        int firstLength = Alignment[0].Length;
        int endColumn = end ?? firstLength;
        var ignore = charsToIgnore is null ? [] : charsToIgnore.ToList();
        if (start < 0 || endColumn > firstLength)
            throw new ArgumentException($"Start ({start}) and end ({endColumn}) are not in the range 0 to {firstLength}");
        string allLetters = AllLetters();
        foreach (string c in ignore)
            allLetters = RemoveAll(allLetters, c);
        var table = eFreqTable is { Count: > 0 } ? eFreqTable : null;

        var infoContent = new List<double>();
        for (int n = start; n < endColumn; n++)
        {
            var freq = LetterFreqs(n, allLetters, ignore, pseudoCount, table, out bool integers);
            infoContent.Add(ColumnInfoContent(freq, table, logBase, integers));
        }

        double total = 0;
        foreach (double v in infoContent)
            total += v;
        IcVector = infoContent;
        return total;
    }

    // _get_letter_freqs; the dictionary keeps the letters in the sorted order and holds integers (0) until it is divided: the flag tells which
    private List<KeyValuePair<char, double>> LetterFreqs(int residueNum, string letters, List<string> toIgnore, double pseudoCount, IReadOnlyDictionary<char, double>? eFreqTable, out bool integers)
    {
        integers = false;
        var freq = letters.Select(l => new KeyValuePair<char, double>(l, 0)).ToList();
        double totalCount = 0;
        if (pseudoCount < 0)
            throw new ArgumentException($"Positive value required for pseudo_count, {PyNumber(pseudoCount)} provided");
        foreach (var record in Alignment)
        {
            if (residueNum >= record.Length)
                throw new ArgumentOutOfRangeException(nameof(residueNum), "index out of range");
            string residue = record.Seq[residueNum].ToString();
            if (toIgnore.Contains(residue))
                continue;
            double weight = Weight(record);
            int at = freq.FindIndex(kv => kv.Key == residue[0]);
            if (at < 0)
                throw new ArgumentException($"Residue {residue} not found in letters {letters}");
            freq[at] = new KeyValuePair<char, double>(residue[0], freq[at].Value + weight);
            totalCount += weight;
        }

        if (eFreqTable is not null)
        {
            foreach (var (key, _) in freq)
                if (key != Gap && !eFreqTable.ContainsKey(key))
                    throw new ArgumentException($"{key} not found in expected frequency table");
        }

        if (totalCount == 0)
        {
            integers = true;                    // the dictionary still holds the integer 0 of dict.fromkeys
            return freq;
        }

        for (int i = 0; i < freq.Count; i++)
        {
            char letter = freq[i].Key;
            if (pseudoCount != 0 && eFreqTable is not null)
            {
                if (!eFreqTable.TryGetValue(letter, out double adjust))
                    throw new KeyNotFoundException($"'{letter}'");
                freq[i] = new KeyValuePair<char, double>(letter, (freq[i].Value + (adjust * pseudoCount)) / (totalCount + pseudoCount));
            }
            else
            {
                freq[i] = new KeyValuePair<char, double>(letter, freq[i].Value / totalCount);
            }
        }

        return freq;
    }

    // _get_column_info_content: random_expected is always None in Biopython, so a column without an expected-frequency table fails
    private static double ColumnInfoContent(List<KeyValuePair<char, double>> observed, IReadOnlyDictionary<char, double>? eFreqTable, double logBase, bool integers)
    {
        if (eFreqTable is not null)
        {
            foreach (var (key, _) in observed)
                if (key != Gap && !eFreqTable.ContainsKey(key))
                    throw new ArgumentException($"Frequency table provided does not contain observed letter {key}");
        }

        double total = 0.0;
        foreach (var (letter, obs) in observed)
        {
            double inner = 0.0;
            if (letter != Gap)
            {
                if (eFreqTable is null)
                    throw new InvalidOperationException($"unsupported operand type(s) for /: '{(integers ? "int" : "float")}' and 'NoneType'");
                double expected = eFreqTable[letter];
                if (expected == 0)
                    throw new DivideByZeroException("float division by zero");
                inner = obs / expected;
            }

            if (inner > 0)
            {
                double denominator = Math.Log(logBase);
                if (logBase <= 0)
                    throw new ArgumentException("math domain error");
                if (denominator == 0)
                    throw new DivideByZeroException("float division by zero");
                total += obs * Math.Log(inner) / denominator;
            }
        }

        return total;
    }

    private static string PyNumber(double v) => v == Math.Floor(v) && Math.Abs(v) < 1e15 ? ((long)v).ToString(CultureInfo.InvariantCulture) : v.ToString("R", CultureInfo.InvariantCulture);

    /// <summary><c>summary.get_column(col)</c>: the letters of an alignment column.</summary>
    public string GetColumn(int col) => Alignment.Column(col);

    /// <summary>
    /// <c>print_info_content(summary_info, fout, rep_record)</c>: three columns, the position, the letter of the representative row and the value of <see cref="IcVector"/> with three decimals; when the vector is empty
    /// <see cref="InformationContent"/> is called without a table first (as Biopython does).
    /// </summary>
    public string PrintInfoContent(int repRecord = 0)
    {
        if (IcVector.Count == 0)
            InformationContent();
        string rep = Alignment[repRecord].Seq;
        var sb = new StringBuilder();
        for (int pos = 0; pos < IcVector.Count; pos++)
            sb.Append(pos.ToString(CultureInfo.InvariantCulture)).Append(' ').Append(rep[pos]).Append(' ').Append(CFormat.Fixed(IcVector[pos], 3)).Append('\n');
        return sb.ToString();
    }
}
