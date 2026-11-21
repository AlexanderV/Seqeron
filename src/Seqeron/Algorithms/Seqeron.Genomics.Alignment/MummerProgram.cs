using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>Which maximal matches <c>mummer</c> reports (<c>-mum</c>, <c>-mumreference</c> = <c>-mumcand</c>, <c>-maxmatch</c>).</summary>
public enum MummerMatchMode
{
    /// <summary><c>-mumreference</c> / <c>-mumcand</c> (the default): maximal matches unique in the reference.</summary>
    MumReference = 0,

    /// <summary><c>-mum</c>: maximal matches unique in the reference and in the query.</summary>
    Mum = 1,

    /// <summary><c>-maxmatch</c>: every maximal match.</summary>
    MaxMatch = 2,
}

/// <summary>The switches of the <c>mummer</c> command (MUMmer 3.23 / Kurtz <c>maxmat3</c>).</summary>
public sealed record MummerOptions
{
    /// <summary>The matching algorithm (<c>-mum</c>, <c>-mumreference</c>, <c>-maxmatch</c>); default <see cref="MummerMatchMode.MumReference"/>.</summary>
    public MummerMatchMode Mode { get; init; } = MummerMatchMode.MumReference;

    /// <summary><c>-n</c>: match only the characters a, c, g, t (any case); every other character of the reference and of the query is replaced by a character that matches nothing.</summary>
    public bool MatchNucleotidesOnly { get; init; }

    /// <summary><c>-l</c>: the minimum length of a match (default 20, must be positive).</summary>
    public long MinLength { get; init; } = 20;

    /// <summary><c>-b</c>: forward and reverse-complement matches.</summary>
    public bool Both { get; init; }

    /// <summary><c>-r</c>: only reverse-complement matches.</summary>
    public bool OnlyReverse { get; init; }

    /// <summary><c>-s</c>: show the matching substring under every match.</summary>
    public bool ShowString { get; init; }

    /// <summary><c>-c</c>: report the query position of a reverse-complement match relative to the original query (needs <c>-b</c> or <c>-r</c>).</summary>
    public bool ShowReversePositions { get; init; }

    /// <summary><c>-F</c>: force the 4 column format (reference id, reference position, query position, length) also for a single reference sequence.</summary>
    public bool FourColumn { get; init; }

    /// <summary><c>-L</c>: show the length of the query sequence on its header line.</summary>
    public bool ShowSequenceLengths { get; init; }

    /// <summary>The reference file name shown in the standard error lines of the program.</summary>
    public string ReferenceName { get; init; } = "reference.fa";

    /// <summary>The names of the query files (stderr lines only); missing names are <c>query.fa</c>.</summary>
    public IReadOnlyList<string>? QueryNames { get; init; }
}

/// <summary>One maximal match: 0-based start in the (concatenated) reference, 0-based start in the (possibly reverse-complemented) query and the length.</summary>
public readonly record struct MummerMatch(int ReferenceStart, int QueryStart, int Length);

/// <summary>A multiple FASTA file as <c>scanmultiplefastafile</c> (<c>maxmatinp.c</c>) reads it: the lower-cased residues of all records joined by one separator byte (255), the first white-space delimited word of every header.</summary>
internal sealed class MummerMultiSeq
{
    public const byte Separator = 255;
    public const byte SubjectReplacement = 252;   // WILDCARD - 2: replaces a non-ACGT residue of the reference with -n
    public const byte QueryReplacement = 251;     // WILDCARD - 3: replaces a non-ACGT residue of the query with -n

    public byte[] Text = [];
    public List<string> Names = [];
    public List<int> DescriptionLengths = [];      // DESCRIPTIONLENGTH: the word plus its terminating newline when the header continues
    public List<int> SeparatorPositions = [];

    public int SequenceCount => Names.Count;

    /// <summary>Scans the file text like the program (every <c>&gt;</c> outside a description starts a record; the residues before the first header belong to the first record).</summary>
    /// <exception cref="FormatException">The file is empty, has no record or only empty records (the program's error messages).</exception>
    public static MummerMultiSeq Scan(string fileName, string input, byte replacement)
    {
        if (input.Length == 0)
            throw new FormatException("mummer: cannot open file \"" + fileName + "\" or file \"" + fileName + "\" is empty");
        var ms = new MummerMultiSeq();
        var text = new List<byte>(input.Length);
        bool indesc = false, copy = false;
        var name = new StringBuilder();
        var seen = false;
        bool terminated = false;
        foreach (char raw in input)
        {
            byte c = (byte)raw;
            if (indesc)
            {
                if (copy)
                {
                    if (IsSpace(c))
                    {
                        copy = false;
                        terminated = true;
                    }
                    else
                        name.Append((char)c);
                }

                if (c == '\n')
                    indesc = false;
            }
            else if (c == '>')
            {
                if (seen)
                    ms.EndDescription(name, terminated);
                if (ms.Count > 0)
                {
                    ms.SeparatorPositions.Add(text.Count);
                    text.Add(Separator);
                }

                seen = true;
                ms.Count++;
                name.Clear();
                terminated = false;
                indesc = true;
                copy = true;
            }
            else if (!IsSpace(c))
            {
                byte l = ToLower(c);
                if (replacement != 0 && l is not ((byte)'a' or (byte)'c' or (byte)'g' or (byte)'t'))
                    l = replacement;
                text.Add(l);
            }
        }

        if (seen)
            ms.EndDescription(name, terminated);
        if (ms.Count == 0)
            throw new FormatException("mummer: no sequences in multiple fasta file");
        if (text.Count == 0)
            throw new FormatException("mummer: empty sequence in multiple fasta file");
        ms.Text = text.ToArray();
        return ms;
    }

    private int Count;

    private void EndDescription(StringBuilder name, bool terminated)
    {
        Names.Add(name.ToString());
        DescriptionLengths.Add(name.Length + (terminated ? 1 : 0));
    }

    /// <summary>The (sequence number, offset) of a position of the joined text (<c>pos2pospair</c>).</summary>
    public (int Sequence, int Offset) PositionPair(int position)
    {
        int k = 0;
        while (k < SeparatorPositions.Count && position > SeparatorPositions[k])
            k++;
        return (k, k == 0 ? position : position - SeparatorPositions[k - 1] - 1);
    }

    /// <summary>The sequences of the joined text (start, length) in order.</summary>
    public IEnumerable<(int Start, int Length)> Sequences()
    {
        int start = 0;
        foreach (int sep in SeparatorPositions)
        {
            yield return (start, sep - start);
            start = sep + 1;
        }

        yield return (start, Text.Length - start);
    }

    internal static bool IsSpace(byte c) => c is (byte)' ' or (byte)'\t' or (byte)'\n' or (byte)'\v' or (byte)'\f' or (byte)'\r';

    internal static byte ToLower(byte c) => c is >= (byte)'A' and <= (byte)'Z' ? (byte)(c + 32) : c;
}

/// <summary>
/// The suffix array of a joined multiple sequence text with the terminator ordered after every byte (the order in which the suffix tree of <c>mummer</c> enumerates its leaves), the inverse array and the
/// longest-common-prefix array, and the maximal-match enumeration of <c>findmaxmatches</c> / <c>findmumcandidates</c> (<c>findmaxmat.c</c>, <c>findmumcand.c</c>) over it.
/// </summary>
internal sealed class MummerIndex
{
    private readonly byte[] _text;
    private readonly int[] _sa;
    private readonly int[] _isa;
    private readonly int[] _lcp;     // _lcp[k] = lcp (suffix _sa[k-1], suffix _sa[k])

    public MummerIndex(byte[] text)
    {
        _text = text;
        int n = text.Length;
        _sa = new int[n];
        _isa = new int[n];
        _lcp = new int[n];
        BuildSuffixArray();
        for (int k = 0; k < n; k++)
            _isa[_sa[k]] = k;
        int h = 0;
        for (int i = 0; i < n; i++)
        {
            int r = _isa[i];
            if (r > 0)
            {
                int j = _sa[r - 1];
                while (i + h < n && j + h < n && text[i + h] == text[j + h])
                    h++;
                _lcp[r] = h;
                if (h > 0)
                    h--;
            }
            else
                h = 0;
        }
    }

    private void BuildSuffixArray()
    {
        int n = _text.Length;
        var rank = new int[n];
        var next = new int[n];
        var keys = new long[n];
        for (int i = 0; i < n; i++)
        {
            _sa[i] = i;
            rank[i] = _text[i];
        }

        int k = 1;
        while (true)
        {
            for (int i = 0; i < n; i++)
                keys[i] = ((long)rank[i] << 32) | (uint)(i + k < n ? rank[i + k] : int.MaxValue);
            // the keys are copied because Array.Sort permutes them together with the items
            var sortKeys = (long[])keys.Clone();
            for (int i = 0; i < n; i++)
                _sa[i] = i;
            Array.Sort(sortKeys, _sa);
            next[_sa[0]] = 0;
            for (int i = 1; i < n; i++)
                next[_sa[i]] = next[_sa[i - 1]] + (sortKeys[i] != sortKeys[i - 1] ? 1 : 0);
            Array.Copy(next, rank, n);
            if (rank[_sa[n - 1]] == n - 1 || k >= n)
                break;
            k <<= 1;
        }
    }

    private int CharAt(int k, int depth)
    {
        int p = _sa[k] + depth;
        return p < _text.Length ? _text[p] : 256;
    }

    /// <summary>Narrows the suffix interval [lo, hi] (whose suffixes agree on the first <paramref name="depth"/> bytes) to those whose next byte is <paramref name="c"/>; false when there is none.</summary>
    private bool Narrow(ref int lo, ref int hi, int depth, byte c)
    {
        int a = lo, b = hi + 1;
        while (a < b)
        {
            int m = (a + b) >> 1;
            if (CharAt(m, depth) < c)
                a = m + 1;
            else
                b = m;
        }

        int newLo = a;
        b = hi + 1;
        while (a < b)
        {
            int m = (a + b) >> 1;
            if (CharAt(m, depth) <= c)
                a = m + 1;
            else
                b = m;
        }

        if (a == newLo)
            return false;
        lo = newLo;
        hi = a - 1;
        return true;
    }

    /// <summary>Every maximal match of length at least <paramref name="minLength"/> between the text and <paramref name="query"/>, in the order <c>mummer</c> prints them.</summary>
    public List<MummerMatch> FindMatches(byte[] query, int minLength, MummerMatchMode mode, CancellationToken cancellationToken = default)
    {
        var result = new List<MummerMatch>();
        int n = _text.Length, m = query.Length;
        if (m < minLength)
            return result;
        bool maxMatch = mode == MummerMatchMode.MaxMatch;
        int lo = 0, hi = n - 1, depth = 0;
        for (int i = 0; i + minLength <= m; i++)
        {
            if ((i & 0x3FF) == 0)
                cancellationToken.ThrowIfCancellationRequested();
            // extend the match of the suffix i: the interval of the prefix of length depth is known
            int limit = maxMatch ? minLength : int.MaxValue;
            while (depth < limit && i + depth < m)
            {
                if (!Narrow(ref lo, ref hi, depth, query[i + depth]))
                    break;
                depth++;
                if (lo == hi && !maxMatch)
                {
                    int j = _sa[lo];
                    while (i + depth < m && j + depth < n && _text[j + depth] == query[i + depth])
                        depth++;
                    break;
                }
            }

            if (depth >= minLength)
            {
                if (maxMatch)
                {
                    for (int k = lo; k <= hi; k++)
                    {
                        int j = _sa[k];
                        if (j == 0 || i == 0 || _text[j - 1] != query[i - 1])
                        {
                            int l = depth;
                            while (i + l < m && j + l < n && _text[j + l] == query[i + l])
                                l++;
                            result.Add(new MummerMatch(j, i, l));
                        }
                    }
                }
                else if (lo == hi)
                {
                    int j = _sa[lo];
                    if (j == 0 || i == 0 || _text[j - 1] != query[i - 1])
                        result.Add(new MummerMatch(j, i, depth));
                }
            }

            // the interval of the suffix i + 1: follow the suffix link
            if (depth <= 1)
            {
                lo = 0;
                hi = n - 1;
                depth = 0;
            }
            else
            {
                int r = _isa[_sa[lo] + 1];
                depth--;
                lo = r;
                hi = r;
                while (lo > 0 && _lcp[lo] >= depth)
                    lo--;
                while (hi < n - 1 && _lcp[hi + 1] >= depth)
                    hi++;
            }
        }

        if (mode != MummerMatchMode.Mum)
            return result;
        return UniqueInQuery(result);
    }

    /// <summary><c>mumuniqueinquery</c> (<c>cleanMUMcand.c</c>): of the matches unique in the reference, the ones that are also unique in the query; sorted by increasing reference start and decreasing length.</summary>
    private static List<MummerMatch> UniqueInQuery(List<MummerMatch> candidates)
    {
        var output = new List<MummerMatch>();
        if (candidates.Count == 0)
            return output;
        var c = candidates.OrderBy(x => x.ReferenceStart).ThenByDescending(x => x.Length).ToList();
        long dbright = 0;
        bool ignorePrevious = false;
        for (int k = 0; k < c.Count; k++)
        {
            bool ignoreCurrent = false;
            long currentRight = (long)c[k].ReferenceStart + c[k].Length - 1;
            if (dbright > currentRight)
                ignoreCurrent = true;
            else if (dbright == currentRight)
            {
                ignoreCurrent = true;
                if (!ignorePrevious && c[k - 1].ReferenceStart == c[k].ReferenceStart)
                    ignorePrevious = true;
            }
            else
                dbright = currentRight;
            if (k > 0 && !ignorePrevious)
                output.Add(c[k - 1]);
            ignorePrevious = ignoreCurrent;
        }

        if (!ignorePrevious)
            output.Add(c[^1]);
        return output;
    }
}

/// <summary>
/// The <c>mummer</c> command of MUMmer 3.23 (Kurtz' <c>maxmat3</c>: <c>src/kurtz/mm3src</c> <c>maxmatopt.c</c>, <c>maxmatinp.c</c>, <c>procmaxmat.c</c>, <c>findmaxmat.c</c>, <c>findmumcand.c</c>): the maximal matches between a
/// reference file and one or more query files, printed per query record as a <c>&gt; id</c> header line and one line per match (reference position [and id], query position, length), for the forward and / or the
/// reverse-complement strand of the query.
/// </summary>
public static class MummerProgram
{
    /// <summary>Runs <c>mummer [options] reference query...</c>.</summary>
    /// <param name="referenceFasta">The text of the reference file (one or more records; they are joined by a separator no match can span).</param>
    /// <param name="queryFastas">The texts of the query files (at most 32), matched in order.</param>
    /// <param name="options">The switches (null = defaults).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>Standard output (matches) and the deterministic lines of standard error (the timing lines <c># CONSTRUCTIONTIME</c>, <c># COMPLETETIME</c>, <c># SPACE</c> are left out).</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The program refuses the switches (<c>-mum</c> with <c>-maxmatch</c>, <c>-b</c> with <c>-r</c>, <c>-c</c> without <c>-b</c> / <c>-r</c>, a non-positive <c>-l</c>, more than 32 query files, no query file).</exception>
    /// <exception cref="FormatException">An input file is empty or has no / only empty records.</exception>
    public static MummerTextOutput Run(string referenceFasta, IReadOnlyList<string> queryFastas, MummerOptions? options = null, CancellationToken cancellationToken = default)
    {
        ArgumentNullException.ThrowIfNull(referenceFasta);
        ArgumentNullException.ThrowIfNull(queryFastas);
        var o = options ?? new MummerOptions();
        if (o.Both && o.OnlyReverse)
            throw new ArgumentException("mummer: option -b and option -r exclude each other");
        if (o.ShowReversePositions && !o.Both && !o.OnlyReverse)
            throw new ArgumentException("mummer: option -c requires either option -b or option -r");
        if (o.MinLength <= 0)
            throw new ArgumentException("mummer: argument " + o.MinLength.ToString(CultureInfo.InvariantCulture) + " for option -l is not a positive integer");
        if (queryFastas.Count == 0)
            throw new ArgumentException("mummer: missing file arguments");
        if (queryFastas.Count > 32)
            throw new ArgumentException("mummer: too many query files, maximal number is 32");
        int minLength = (int)Math.Min(o.MinLength, int.MaxValue);
        bool forward = !o.OnlyReverse, reverse = o.Both || o.OnlyReverse;

        var messages = new StringBuilder();
        var reference = MummerMultiSeq.Scan(o.ReferenceName, referenceFasta, o.MatchNucleotidesOnly ? MummerMultiSeq.SubjectReplacement : (byte)0);
        messages.Append("# reading input file \"").Append(o.ReferenceName).Append("\" of length ").Append(reference.Text.Length.ToString(CultureInfo.InvariantCulture)).Append('\n');
        messages.Append("# construct suffix tree for sequence of length ").Append(reference.Text.Length.ToString(CultureInfo.InvariantCulture)).Append('\n');
        messages.Append("# (maximum reference length is 536870908)\n# (maximum query length is 4294967295)\n");
        var index = new MummerIndex(reference.Text);
        int maxDescription = reference.DescriptionLengths.Max();
        bool singleColumn = reference.SequenceCount == 1 && !o.FourColumn;

        var sb = new StringBuilder();
        for (int f = 0; f < queryFastas.Count; f++)
        {
            string qname = o.QueryNames != null && f < o.QueryNames.Count ? o.QueryNames[f] : "query.fa";
            var queries = MummerMultiSeq.Scan(qname, queryFastas[f], o.MatchNucleotidesOnly ? MummerMultiSeq.QueryReplacement : (byte)0);
            messages.Append("# reading input file \"").Append(qname).Append("\" of length ").Append(queries.Text.Length.ToString(CultureInfo.InvariantCulture)).Append('\n');
            messages.Append("# matching query-file \"").Append(qname).Append("\"\n# against subject-file \"").Append(o.ReferenceName).Append("\"\n");
            int number = 0;
            foreach (var (start, length) in queries.Sequences())
            {
                cancellationToken.ThrowIfCancellationRequested();
                var query = new byte[length];
                Array.Copy(queries.Text, start, query, 0, length);
                for (int strand = 0; strand < 2; strand++)
                {
                    bool isReverse = strand == 1;
                    if (isReverse ? !reverse : !forward)
                        continue;
                    sb.Append("> ").Append(queries.Names[number]);
                    if (isReverse)
                        sb.Append(" Reverse");
                    if (o.ShowSequenceLengths)
                        sb.Append("  Len = ").Append(length.ToString(CultureInfo.InvariantCulture));
                    sb.Append('\n');
                    byte[] q = isReverse ? ReverseComplement(query) : query;
                    foreach (var match in index.FindMatches(q, minLength, o.Mode, cancellationToken))
                    {
                        if (singleColumn)
                            sb.Append((match.ReferenceStart + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)).Append("  ");
                        else
                        {
                            var (seq, offset) = reference.PositionPair(match.ReferenceStart);
                            sb.Append("  ").Append(reference.Names[seq]).Append(' ', maxDescription - reference.DescriptionLengths[seq])
                                .Append("  ").Append((offset + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)).Append("  ");
                        }

                        long qpos = isReverse && o.ShowReversePositions ? length - match.QueryStart : match.QueryStart + 1L;
                        sb.Append(qpos.ToString(CultureInfo.InvariantCulture).PadLeft(8)).Append("  ")
                            .Append(match.Length.ToString(CultureInfo.InvariantCulture).PadLeft(8)).Append('\n');
                        if (o.ShowString)
                        {
                            for (int k = 0; k < match.Length; k++)
                                sb.Append((char)reference.Text[match.ReferenceStart + k]);
                            sb.Append('\n');
                        }
                    }
                }

                number++;
            }
        }

        return new MummerTextOutput(sb.ToString(), messages.ToString());
    }

    /// <summary><c>wccSequence</c> of <c>procmaxmat.c</c>: the Watson-Crick complement read backwards; the replacement characters of <c>-n</c> are kept, an unknown character becomes <c>n</c>.</summary>
    internal static byte[] ReverseComplement(byte[] s)
    {
        var r = new byte[s.Length];
        for (int i = 0; i < s.Length; i++)
            r[s.Length - 1 - i] = Complement(s[i]);
        return r;
    }

    private static byte Complement(byte c)
    {
        if (c >= MummerMultiSeq.QueryReplacement)
            return c;
        return (char)c switch
        {
            'a' => (byte)'t', 'c' => (byte)'g', 'g' => (byte)'c', 't' => (byte)'a', 'r' => (byte)'y', 'y' => (byte)'r', 's' => (byte)'s', 'w' => (byte)'w',
            'm' => (byte)'k', 'k' => (byte)'m', 'b' => (byte)'v', 'd' => (byte)'h', 'h' => (byte)'d', 'v' => (byte)'b', _ => (byte)'n',
        };
    }
}
