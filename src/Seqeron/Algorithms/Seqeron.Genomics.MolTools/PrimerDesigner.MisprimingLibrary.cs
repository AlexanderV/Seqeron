using System.Globalization;

namespace Seqeron.Genomics.MolTools;

// Primer3 mispriming library (PRIMER_MISPRIMING_LIBRARY): libprimer3.c oligo_repeat_library_mispriming /
// pair_repeat_sim, p3_seq_lib.c (add_seq_to_seq_lib, parse_seq_name, upcase_and_check_char,
// reverse_complement_seq_lib), dpal.c (set_dpal_args, dpal_set_ambiguity_code_matrix,
// _dpal_long_nopath_maxgap1_local_end). Bit-exact to primer3-py 2.3.1 design_primers(misprime_lib=…).
public static partial class PrimerDesigner
{
    /// <summary>Primer3 default PRIMER_MAX_LIBRARY_MISPRIMING (<c>pr_set_default_global_args_1</c>: p_args.max_repeat_compl = 12.00).</summary>
    public const double Primer3MaxLibraryMispriming = 12.0;

    /// <summary>Primer3 default PRIMER_PAIR_MAX_LIBRARY_MISPRIMING (<c>pair_repeat_compl</c> = 24.00).</summary>
    public const double Primer3PairMaxLibraryMispriming = 24.0;

    /// <summary>
    /// Primer3 library mispriming score of one primer (PRIMER_LEFT/RIGHT_n_LIBRARY_MISPRIMING "score, name";
    /// <c>libprimer3.c</c> <c>oligo_repeat_library_mispriming</c>): for every library entry i (the caller's entries
    /// followed by their "reverse …" reverse complements) w_i = weight_i · <c>align</c>(primer, y_i), where a left
    /// (forward) primer is aligned with entry i and a right (reverse) primer — given 5′→3′ — with the reverse
    /// complement of entry i. Alignment: dpal with Primer3's primer-picking scores (+1 match, −1 mismatch, −0.25 against
    /// N, −2 per single-base gap, max gap 1, floored at 0) anchored at the primer's 3′ end (<c>DPAL_LOCAL_END</c>). By
    /// default (PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS = 0, Primer3's default) an IUPAC code in an entry can never align
    /// and N scores −0.25 — and, as in Primer3, a right primer is then aligned with <c>DPAL_LOCAL</c> (not end-anchored);
    /// with <paramref name="ambiguityCodesConsensus"/> = true (1) an IUPAC code scores +1 against every base it contains
    /// (<c>dpal_set_ambiguity_code_matrix</c>; N matches anything) and both primers use <c>DPAL_LOCAL_END</c>. An entry shorter than 3 nt scores its length. The reported entry is
    /// Primer3's <c>repeat_sim.max</c>: the first entry whose w exceeds the integer part of the running maximum
    /// (score 0 → the first entry).
    /// </summary>
    /// <param name="primer">Primer, 5′→3′ (case-insensitive).</param>
    /// <param name="isForward">True for a left (forward) primer, false for a right (reverse) primer.</param>
    /// <param name="library">The mispriming library.</param>
    /// <param name="ambiguityCodesConsensus">PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS (default false = 0, Primer3's default).</param>
    /// <returns>The library mispriming score and entry name.</returns>
    /// <exception cref="ArgumentNullException">Null primer or library.</exception>
    /// <exception cref="ArgumentException">Empty primer.</exception>
    public static LibraryMisprimingScore CalculateLibraryMispriming(
        string primer, bool isForward, PrimerMisprimingLibrary library, bool ambiguityCodesConsensus = false)
    {
        ArgumentNullException.ThrowIfNull(primer);
        ArgumentNullException.ThrowIfNull(library);
        if (primer.Length == 0)
            throw new ArgumentException("Primer cannot be empty.", nameof(primer));
        var r = ComputeLibraryMispriming(primer.ToUpperInvariant(), isForward, library, ambiguityCodesConsensus,
            double.PositiveInfinity);
        return new LibraryMisprimingScore(r.MaxScore, r.Name);
    }

    // Per-primer repeat_sim (scores of every library entry, Primer3's max index) and the max_repeat_compl decision.
    internal sealed class LibraryMispriming
    {
        public required double[] Scores { get; init; }
        public required int MaxIndex { get; init; }
        public required string Name { get; init; }
        public required bool Exceeds { get; init; }
        public double MaxScore => Scores[MaxIndex];
    }

    // oligo_repeat_library_mispriming. Every entry is scored (Primer3 stops at the first entry over the limit only
    // for a primer it then rejects; the accept/reject decision is the same).
    internal static LibraryMispriming ComputeLibraryMispriming(
        string primer, bool isLeft, PrimerMisprimingLibrary lib, bool consensus, double maxLibraryMispriming)
    {
        // max_lib_compl = (short) p_args.max_repeat_compl (C truncation toward zero, then 16-bit wrap).
        int maxLibCompl = double.IsPositiveInfinity(maxLibraryMispriming)
            ? int.MaxValue
            : unchecked((short)(int)maxLibraryMispriming);
        int n = lib.Names.Count;
        var scores = new double[n];
        int maxIndex = 0, max = 0;
        bool exceeds = false;
        for (int i = 0; i < n; i++)
        {
            // Left: align(s, seqs[i], local_end[_ambig]); right: align(s_r, rev_compl_seqs[i], consensus ?
            // local_end_ambig : local) — Primer3 uses the non-anchored LOCAL alignment there.
            var matrix = consensus ? DpalAmbiguityMatrix : DpalPrimerMatrix;
            double a = isLeft
                ? AlignLibrary(primer, lib.SequenceAt(i), matrix, localEnd: true)
                : AlignLibrary(primer, lib.ReverseComplementAt(i), matrix, localEnd: consensus);
            double w = lib.WeightAt(i) * a;
            if (w > short.MaxValue || w < short.MinValue)
                throw new InvalidOperationException("Out of range error occured calculating match to repeat library (Primer3).");
            scores[i] = w;
            if (w > max)
            {
                max = (int)w;
                maxIndex = i;
            }
            if (w > maxLibCompl)
                exceeds = true;
        }
        return new LibraryMispriming { Scores = scores, MaxIndex = maxIndex, Name = lib.Names[maxIndex], Exceeds = exceeds };
    }

    // pair_repeat_sim: max over entries of (int)(left.score[i] + right.score[i]) (starting at 0), with the first entry
    // reaching the maximum.
    internal static (int Score, string Name) PairLibraryMispriming(LibraryMispriming left, LibraryMispriming right,
        PrimerMisprimingLibrary lib)
    {
        int max = 0;
        string name = lib.Names[0];
        for (int i = 0; i < lib.Names.Count; i++)
        {
            int w = (int)(left.Scores[i] + right.Scores[i]);
            if (w > max)
            {
                max = w;
                name = lib.Names[i];
            }
        }
        return (max, name);
    }

    // libprimer3.c align() for DPAL_LOCAL / DPAL_LOCAL_END: a second sequence shorter than 3 scores its length;
    // otherwise the dpal score / 100, floored at 0 (dpal returns 0 for an empty first sequence).
    private static double AlignLibrary(string x, string y, int[] matrix, bool localEnd)
    {
        if (y.Length < 3)
            return y.Length;
        if (x.Length == 0)
            return 0;
        int smax = localEnd ? DpalLocalEnd(x, y, matrix) : DpalLocalFast(x, y, matrix: matrix);
        return smax < 0 ? 0.0 : smax / 100.0;
    }

    // dpal.c _dpal_long_nopath_maxgap1_local_end (score only, max gap 1, cells floored at 0; the optimum is taken in
    // the last row only, i.e. the alignment ends at the 3′-terminal base of X). Line-by-line port, including the
    // C behaviour for |X| ≤ 2 (row 1 reads X[1], the NUL terminator for |X| = 1, whose ssm row is INT_MIN; the
    // "last row" pass then recomputes row |X| − 1 from rows 0 and 1).
    private static int DpalLocalEnd(string x, string y, int[] m)
    {
        int xlen = x.Length, ylen = y.Length;
        const int gap = DpalGap;
        var s0 = new int[ylen];
        var s1 = new int[ylen];
        var s2 = new int[ylen];
        int smax = 0, score, a;
        char X(int i) => i < xlen ? x[i] : '\0';

        for (int j = 0; j < ylen; j++)
        {
            score = MatrixSsm(m, x[0], y[j]);
            if (score < 0) score = 0;
            s0[j] = score;
        }

        score = MatrixSsm(m, X(1), y[0]);
        if (score < 0) score = 0;
        s1[0] = score;
        for (int j = 1; j < ylen; j++)
        {
            score = s0[j - 1];
            if (j > 1 && (a = s0[j - 2] + gap) > score) score = a;
            score += MatrixSsm(m, X(1), y[j]);
            if (score < 0) score = 0;
            s1[j] = score;
        }

        for (int i = 2; i < xlen - 1; i++)
        {
            score = MatrixSsm(m, x[i], y[0]);
            if (score < 0) score = 0;
            s2[0] = score;
            score = s1[0];
            if ((a = s0[0] + gap) > score) score = a;
            score += MatrixSsm(m, x[i], y[1]);
            if (score < 0) score = 0;
            s2[1] = score;
            for (int j = 2; j < ylen; j++)
            {
                score = s0[j - 1];
                if ((a = s1[j - 2]) > score) score = a;
                score += gap;
                if ((a = s1[j - 1]) > score) score = a;
                score += MatrixSsm(m, x[i], y[j]);
                if (score < 0) score = 0;
                s2[j] = score;
            }
            (s0, s1, s2) = (s1, s2, s0);
        }

        // Last row (i = xlen − 1): the only row that updates the optimum (its cells are not needed afterwards; a
        // negative cell, floored at 0 in C, cannot raise the optimum either).
        int last = xlen - 1;
        score = MatrixSsm(m, x[last], y[0]);
        if (score > smax) smax = score;
        score = s1[0];
        if ((a = s0[0] + gap) > score) score = a;
        score += MatrixSsm(m, x[last], y[1]);
        if (score > smax) smax = score;
        for (int j = 2; j < ylen; j++)
        {
            score = s0[j - 1];
            if ((a = s1[j - 2]) > score) score = a;
            score += gap;
            if ((a = s1[j - 1]) > score) score = a;
            score += MatrixSsm(m, x[last], y[j]);
            if (score > smax) smax = score;
        }
        return smax;
    }

    // 7-bit dpal substitution matrix lookup; any other character has ssm INT_MIN (dpal_ssm is indexed by unsigned char
    // and only ACGTN / IUPAC rows are set).
    private static int MatrixSsm(int[] m, char x, char y) => x < 128 && y < 128 ? m[(x << 7) | y] : int.MinValue;

    // dpal.c set_dpal_args: A/C/G/T/N only (N −25 against everything, including N), +100 match, −100 mismatch; every
    // other character INT_MIN.
    private static readonly int[] DpalPrimerMatrix = BuildDpalPrimerMatrix();

    // dpal_set_ambiguity_code_matrix applied to set_dpal_args (Primer3 create_dpal_arg_holder local_ambig /
    // local_end_ambig): for IUPAC codes B D H V R Y K M S W N the score against a code or base is the maximum score
    // between the bases they represent.
    private static readonly int[] DpalAmbiguityMatrix = BuildDpalAmbiguityMatrix();

    private static int[] BuildDpalPrimerMatrix()
    {
        var m = new int[128 * 128];
        const string legal = "ACGTN";
        for (int i = 0; i < 128; i++)
            for (int j = 0; j < 128; j++)
            {
                int v = int.MinValue;
                if (i != 0 && j != 0 && legal.Contains((char)i) && legal.Contains((char)j))
                {
                    if (i == 'N' || j == 'N') v = DpalN;
                    else v = i == j ? DpalMatch : DpalMismatch;
                }
                m[(i << 7) | j] = v;
            }
        return m;
    }

    private static int[] BuildDpalAmbiguityMatrix()
    {
        var m = BuildDpalPrimerMatrix();
        const string ambCodes = "BDHVRYKMSWN", allBases = "ACGT";
        foreach (char c1 in ambCodes)
        {
            string bases1 = XlateAmbiguityCode(c1);
            foreach (char c2 in ambCodes)
            {
                string bases2 = XlateAmbiguityCode(c2);
                int extreme = int.MinValue;
                foreach (char b1 in bases1)
                    foreach (char b2 in bases2)
                        extreme = Math.Max(extreme, m[(b1 << 7) | b2]);
                m[(c1 << 7) | c2] = extreme;
            }
            foreach (char b2 in allBases)
            {
                int extreme = int.MinValue;
                foreach (char b1 in bases1)
                    extreme = Math.Max(extreme, m[(b1 << 7) | b2]);
                m[(c1 << 7) | b2] = extreme;
                m[(b2 << 7) | c1] = extreme;
            }
        }
        return m;
    }

    // dpal.c xlate_ambiguity_code.
    private static string XlateAmbiguityCode(char c) => c switch
    {
        'N' => "ACGT",
        'B' => "CGT",
        'D' => "AGT",
        'H' => "ACT",
        'V' => "ACG",
        'R' => "AG",
        'Y' => "CT",
        'K' => "GT",
        'M' => "AC",
        'S' => "CG",
        'W' => "AT",
        _ => throw new ArgumentOutOfRangeException(nameof(c)),
    };
}

/// <summary>
/// Primer3 library mispriming score of a primer or pair: the score (PRIMER_*_LIBRARY_MISPRIMING value) and the name of
/// the library entry it refers to ("reverse &lt;name&gt;" for an entry's reverse complement).
/// </summary>
/// <param name="Score">Weighted dpal score.</param>
/// <param name="Name">Library entry name.</param>
public readonly record struct LibraryMisprimingScore(double Score, string Name);

/// <summary>
/// A Primer3 mispriming library (PRIMER_MISPRIMING_LIBRARY; primer3-py <c>misprime_lib</c> name → sequence dictionary),
/// built as Primer3's <c>p3_seq_lib.c</c> does (<c>add_seq_to_seq_lib</c>, <c>reverse_complement_seq_lib</c>):
/// <list type="bullet">
/// <item>Each name may carry a weight after a <c>*</c> (e.g. <c>"Alu*2.5"</c>, <c>parse_seq_name</c>: the number
/// following the first <c>*</c>, read as C <c>strtod</c>; no <c>*</c> → weight 1; a missing number, a negative weight or
/// one above 100 is illegal). The full name, including the weight suffix, is reported.</item>
/// <item>Sequences are upper-cased; spaces, tabs, CR and LF are removed; A/C/G/T/N and the IUPAC codes
/// B D H V R Y K M S W are kept; any other character becomes N and is reported in <see cref="Warnings"/>
/// (<c>upcase_and_check_char</c>). An empty sequence is illegal.</item>
/// <item>After the caller's n entries Primer3 appends n entries "reverse &lt;name&gt;" holding their reverse
/// complements (IUPAC-aware, <c>p3_reverse_complement</c>) with the same weights; <see cref="Names"/> lists all 2n in
/// that order, which is the order Primer3 scans (and reports ties in).</item>
/// </list>
/// </summary>
public sealed class PrimerMisprimingLibrary
{
    private const double MaxLibraryWeight = 100.0; // PR_MAX_LIBRARY_WT
    private readonly string[] _seqs;
    private readonly double[] _weights;

    /// <summary>Builds a library from (name, sequence) entries, in order (primer3-py <c>misprime_lib</c>).</summary>
    /// <param name="entries">Library entries: name (optionally with a <c>*weight</c> suffix) → sequence.</param>
    /// <exception cref="ArgumentNullException">Null entries, name or sequence.</exception>
    /// <exception cref="ArgumentException">Empty sequence or illegal weight (Primer3 <c>add_seq_to_seq_lib</c>).</exception>
    public PrimerMisprimingLibrary(IEnumerable<KeyValuePair<string, string>> entries)
    {
        ArgumentNullException.ThrowIfNull(entries);
        var names = new List<string>();
        var seqs = new List<string>();
        var weights = new List<double>();
        var warnings = new List<string>();
        foreach (var (name, seq) in entries)
        {
            if (name is null || seq is null)
                throw new ArgumentNullException(nameof(entries), "Library entry names and sequences cannot be null.");
            double w = ParseWeight(name);
            if (w < 0)
                throw new ArgumentException($"Illegal weight in mispriming library entry '{name}' (Primer3: must be 0–100).", nameof(entries));
            if (seq.Length == 0)
                throw new ArgumentException($"Empty sequence in mispriming library entry '{name}'.", nameof(entries));
            var (clean, offender) = UpcaseAndCheck(seq);
            if (offender != '\0')
                warnings.Add($"Unrecognized character ({offender}) in mispriming library, entry {name}");
            names.Add(name);
            seqs.Add(clean);
            weights.Add(w);
        }

        int n = names.Count;
        var allNames = new string[2 * n];
        _seqs = new string[2 * n];
        _weights = new double[2 * n];
        for (int i = 0; i < n; i++)
        {
            allNames[i] = names[i];
            _seqs[i] = seqs[i];
            _weights[i] = weights[i];
            allNames[n + i] = "reverse " + names[i];
            _seqs[n + i] = DnaSequence.GetReverseComplementString(seqs[i]);
            _weights[n + i] = weights[i];
        }
        Count = n;
        Names = Array.AsReadOnly(allNames);
        Warnings = warnings.AsReadOnly();
    }

    /// <summary>Number of caller-supplied entries (Primer3 scans twice as many, see <see cref="Names"/>).</summary>
    public int Count { get; }

    /// <summary>All 2·<see cref="Count"/> entry names in Primer3's scan order: the caller's, then "reverse &lt;name&gt;".</summary>
    public IReadOnlyList<string> Names { get; }

    /// <summary>Primer3 warnings (unrecognized sequence characters replaced by N).</summary>
    public IReadOnlyList<string> Warnings { get; }

    /// <summary>The (cleaned) sequence of entry <paramref name="index"/> in <see cref="Names"/> order.</summary>
    public string SequenceAt(int index) => _seqs[index];

    /// <summary>The weight of entry <paramref name="index"/> in <see cref="Names"/> order.</summary>
    public double WeightAt(int index) => _weights[index];

    // Primer3 rev_compl_seqs[i]: the reverse complement of entry i (= entry i ± n).
    internal string ReverseComplementAt(int index) => _seqs[index < Count ? index + Count : index - Count];

    // p3_seq_lib.c parse_seq_name: the strtod value after the first '*' (no '*' → 1; no number or > 100 → −1).
    private static double ParseWeight(string name)
    {
        int star = name.IndexOf('*');
        if (star < 0)
            return 1;
        if (!TryStrtod(name, star + 1, out double v))
            return -1;
        return v > MaxLibraryWeight ? -1 : v;
    }

    // C strtod prefix parse (leading white space, sign, decimal digits with optional point and exponent, INF/INFINITY,
    // NAN); false when no number is found. (Hexadecimal floats are read as their leading "0".)
    private static bool TryStrtod(string s, int start, out double value)
    {
        value = 0;
        int i = start;
        while (i < s.Length && s[i] is ' ' or '\t' or '\n' or '\v' or '\f' or '\r') i++;
        int numStart = i;
        if (i < s.Length && s[i] is '+' or '-') i++;
        int afterSign = i;
        string rest = s.Substring(afterSign);
        bool negative = afterSign > numStart && s[numStart] == '-';
        if (rest.StartsWith("inf", StringComparison.OrdinalIgnoreCase))
        {
            value = negative ? double.NegativeInfinity : double.PositiveInfinity;
            return true;
        }
        if (rest.StartsWith("nan", StringComparison.OrdinalIgnoreCase))
        {
            value = double.NaN;
            return true;
        }
        int digits = 0;
        while (i < s.Length && char.IsAsciiDigit(s[i])) { i++; digits++; }
        if (i < s.Length && s[i] == '.')
        {
            i++;
            while (i < s.Length && char.IsAsciiDigit(s[i])) { i++; digits++; }
        }
        if (digits == 0)
            return false;
        if (i < s.Length && s[i] is 'e' or 'E')
        {
            int j = i + 1;
            if (j < s.Length && s[j] is '+' or '-') j++;
            if (j < s.Length && char.IsAsciiDigit(s[j]))
            {
                while (j < s.Length && char.IsAsciiDigit(s[j])) j++;
                i = j;
            }
        }
        value = double.Parse(s.AsSpan(numStart, i - numStart), NumberStyles.Float, CultureInfo.InvariantCulture);
        return true;
    }

    // p3_seq_lib.c upcase_and_check_char.
    private static (string Clean, char Offender) UpcaseAndCheck(string s)
    {
        var sb = new System.Text.StringBuilder(s.Length);
        char offender = '\0';
        foreach (char c in s)
        {
            switch (c)
            {
                case '\n' or ' ' or '\t' or '\r':
                    break;
                case 'a' or 'g' or 'c' or 't' or 'n' or 'A' or 'G' or 'C' or 'T' or 'N'
                    or 'b' or 'B' or 'd' or 'D' or 'h' or 'H' or 'v' or 'V' or 'r' or 'R'
                    or 'y' or 'Y' or 'k' or 'K' or 'm' or 'M' or 's' or 'S' or 'w' or 'W':
                    sb.Append(char.ToUpperInvariant(c));
                    break;
                default:
                    if (offender == '\0') offender = c;
                    sb.Append('N');
                    break;
            }
        }
        return (sb.ToString(), offender);
    }
}
