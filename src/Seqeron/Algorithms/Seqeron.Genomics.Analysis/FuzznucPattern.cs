using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// A parsed EMBOSS <c>fuzznuc</c> nucleotide pattern (the <c>-pattern</c> language): IUPAC letters, character classes <c>[ACG]</c>, exclusions
/// <c>{AG}</c>, repeats <c>N(3)</c> and length ranges <c>A(2,4)</c>, optional <c>-</c> separators, the anchors <c>&lt;</c> (start of the sequence) and
/// <c>&gt;</c> (end of the sequence) and an optional terminating period — e.g. <c>[CG](5)TG{A}N(1,5)C</c>.
/// </summary>
/// <remarks>
/// <para>
/// <b>Grammar</b> (EMBOSS 6.6.0 <c>fuzznuc -help</c>, <c>-pattern</c>): <c>pattern = ['&lt;'] element { ['-'] element } ['&gt;'] ['.']</c>;
/// <c>element = atom ['(' n [',' m] ')']</c>; <c>atom = IUPAC letter | '[' bases ']' | '{' bases '}'</c> with bases drawn from A, C, G, T (the binary rejects IUPAC
/// ambiguity letters inside brackets: "Nesting '[' in [] not allowed"). A single repeat <c>(n)</c> needs n ≥ 1; a range <c>(n,m)</c> needs 0 ≤ n ≤ m, m ≥ 1, and n = 0 makes the
/// element optional (<c>CA(0,2)G</c> matches <c>CG</c>, <c>CAG</c> and <c>CAAG</c>). Case is ignored. <c>{ACGT}</c> (no base allowed) is legal and matches nothing. Anything else
/// — an empty <c>[]</c> / <c>{}</c>, an unbalanced bracket, <c>(2,1)</c>, <c>(0)</c>, <c>(3,)</c>, an anchor in the middle, a letter that is not an IUPAC DNA code (U included) —
/// is rejected with an <see cref="ArgumentException"/> naming the position. (The binary accepts several malformed spellings in unspecified ways — <c>A&lt;C</c> reads as <c>AC</c>,
/// <c>{}</c> matches everything, <c>A(0,0)</c> reports empty hits; those are not part of the documented language and are not reproduced.)
/// </para>
/// <para>
/// <b>Compilation.</b> A pattern is a sequence of elements, each a set of allowed bases (a 4-bit mask over A, C, G, T) repeated between <see cref="PatternElement.MinRepeat"/> and
/// <see cref="PatternElement.MaxRepeat"/> times. Every non-empty base set is exactly one of the 15 IUPAC codes, so a pattern expands (<see cref="Expand"/>) into a finite set of
/// fixed-length IUPAC strings — one per choice of the repeat counts — which the canonical matchers (<see cref="ApproximateMatcher"/> /
/// <see cref="MotifFinder.FindDegenerateMotif(string, string, CancellationToken)"/>) scan; <see cref="GenomicAnalyzer.FindPatternHits(DnaSequence, string, int, bool, CancellationToken)"/>
/// unites their hits. Nothing here scans a sequence.
/// </para>
/// </remarks>
public sealed class FuzznucPattern
{
    /// <summary>Largest number of repeat-count combinations an <see cref="Expand"/> may enumerate (a pattern with more is rejected).</summary>
    public const int MaxCombinations = 100_000;

    private const string Bases = "ACGT";
    private const string IupacByPreference = "ACGTRYSWKMBDHVN"; // the 15 IUPAC DNA codes

    private readonly PatternElement[] _elements;

    private FuzznucPattern(string source, bool anchoredAtStart, bool anchoredAtEnd, PatternElement[] elements)
    {
        Source = source;
        AnchoredAtStart = anchoredAtStart;
        AnchoredAtEnd = anchoredAtEnd;
        _elements = elements;
    }

    /// <summary>The pattern text as given (trimmed).</summary>
    public string Source { get; }

    /// <summary>The pattern starts with <c>&lt;</c>: a hit must start at the first base of the strand it is searched on.</summary>
    public bool AnchoredAtStart { get; }

    /// <summary>The pattern ends with <c>&gt;</c>: a hit must end at the last base of the strand it is searched on.</summary>
    public bool AnchoredAtEnd { get; }

    /// <summary>The elements in pattern order.</summary>
    public IReadOnlyList<PatternElement> Elements => _elements;

    /// <summary>Length of the shortest hit (sum of the minimum repeats).</summary>
    public int MinLength => _elements.Sum(e => e.MinRepeat);

    /// <summary>Length of the longest hit (sum of the maximum repeats).</summary>
    public int MaxLength => _elements.Sum(e => e.MaxRepeat);

    /// <summary>True when every element has a fixed repeat count (all hits have the same length).</summary>
    public bool IsFixedLength => _elements.All(e => e.MinRepeat == e.MaxRepeat);

    /// <summary>
    /// Parses a pattern.
    /// </summary>
    /// <param name="pattern">The pattern text (see the class remarks for the grammar).</param>
    /// <exception cref="ArgumentNullException"><paramref name="pattern"/> is null.</exception>
    /// <exception cref="ArgumentException">The pattern is empty or violates the grammar; the message names the position.</exception>
    public static FuzznucPattern Parse(string pattern)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        return ParseCore(pattern, out string? error) ?? throw new ArgumentException(error, nameof(pattern));
    }

    /// <summary>Parses a pattern without throwing.</summary>
    /// <param name="pattern">The pattern text.</param>
    /// <param name="result">The parsed pattern, or null.</param>
    /// <param name="error">The reason, or null on success.</param>
    /// <returns>True when the pattern is valid.</returns>
    public static bool TryParse(string? pattern, out FuzznucPattern? result, out string? error)
    {
        if (pattern is null)
        {
            result = null;
            error = "The pattern is null.";
            return false;
        }

        result = ParseCore(pattern, out error);
        return result is not null;
    }

    private static FuzznucPattern? ParseCore(string pattern, out string? error)
    {
        error = null;
        string text = pattern.Trim();
        int i = 0;
        bool start = false, end = false;
        var elements = new List<PatternElement>();

        if (text.Length > 0 && text[0] == '<')
        {
            start = true;
            i++;
        }

        while (i < text.Length)
        {
            char c = text[i];
            if (c == '-')
            {
                i++;
                continue;
            }

            if (c == '.')
            {
                if (text[(i + 1)..].Trim().Length > 0)
                {
                    error = $"Illegal pattern '{pattern}': the period that ends a pattern must be last (position {i}).";
                    return null;
                }

                break;
            }

            if (c == '>')
            {
                string rest = text[(i + 1)..].TrimStart();
                if (rest.Length > 0 && rest != ".")
                {
                    error = $"Illegal pattern '{pattern}': '>' (end of sequence) must end the pattern (position {i}).";
                    return null;
                }

                end = true;
                break;
            }

            if (c == '<')
            {
                error = $"Illegal pattern '{pattern}': '<' (start of sequence) must begin the pattern (position {i}).";
                return null;
            }

            int mask;
            if (c is '[' or '{')
            {
                char close = c == '[' ? ']' : '}';
                int closeAt = text.IndexOf(close, i + 1);
                if (closeAt < 0)
                {
                    error = $"Illegal pattern '{pattern}': missing '{close}' for the '{c}' at position {i}.";
                    return null;
                }

                int listed = 0;
                for (int j = i + 1; j < closeAt; j++)
                {
                    int bit = Bases.IndexOf(char.ToUpperInvariant(text[j]));
                    if (bit < 0)
                    {
                        error = $"Illegal pattern '{pattern}': '{text[j]}' at position {j} is not one of A, C, G, T (a class or exclusion lists plain bases; ambiguity letters cannot be nested).";
                        return null;
                    }

                    listed |= 1 << bit;
                }

                if (closeAt == i + 1)
                {
                    error = $"Illegal pattern '{pattern}': empty '{c}{close}' at position {i}.";
                    return null;
                }

                mask = c == '[' ? listed : 0xF & ~listed;
                i = closeAt + 1;
            }
            else if (IupacByPreference.Contains(char.ToUpperInvariant(c)))
            {
                mask = MaskOf(char.ToUpperInvariant(c));
                i++;
            }
            else
            {
                error = $"Illegal pattern '{pattern}': '{c}' at position {i} is not an IUPAC nucleotide code, '[', '{{', '<', '>' or '-'.";
                return null;
            }

            int min = 1, max = 1;
            bool isRange = false;
            if (i < text.Length && text[i] == '(')
            {
                int closeAt = text.IndexOf(')', i + 1);
                if (closeAt < 0)
                {
                    error = $"Illegal pattern '{pattern}': missing ')' for the '(' at position {i}.";
                    return null;
                }

                string[] parts = text.Substring(i + 1, closeAt - i - 1).Split(',');
                if (parts.Length > 2 || !TryCount(parts[0], out min) || (parts.Length == 2 && !TryCount(parts[1], out max)))
                {
                    error = $"Illegal pattern '{pattern}': bad repeat '({text.Substring(i + 1, closeAt - i - 1)})' at position {i}; use (n) or (n,m) with whole numbers.";
                    return null;
                }

                if (parts.Length == 1)
                {
                    max = min;
                    if (min < 1)
                    {
                        error = $"Illegal pattern '{pattern}': a single repeat count must be at least 1 (position {i}).";
                        return null;
                    }
                }
                else if (min > max || max < 1)
                {
                    error = $"Illegal pattern '{pattern}': the repeat range ({min},{max}) at position {i} must satisfy 0 <= n <= m and m >= 1.";
                    return null;
                }

                isRange = parts.Length == 2;
                i = closeAt + 1;
            }

            elements.Add(new PatternElement(mask, min, max, isRange));
        }

        if (elements.Count == 0)
        {
            error = $"Illegal pattern '{pattern}': it has no element.";
            return null;
        }

        return new FuzznucPattern(text, start, end, elements.ToArray());

        static bool TryCount(string s, out int value)
            => int.TryParse(s.Trim(), System.Globalization.NumberStyles.None, System.Globalization.CultureInfo.InvariantCulture, out value) && value <= 10_000;
    }

    // Base set of an IUPAC DNA letter from the canonical NC-IUB table (IupacHelper.MatchesIupac), bit i = Bases[i].
    private static int MaskOf(char code)
    {
        int mask = 0;
        for (int i = 0; i < Bases.Length; i++)
            if (IupacHelper.MatchesIupac(Bases[i], code))
                mask |= 1 << i;
        return mask;
    }

    // The IUPAC letter whose base set is `mask` (1..15); '\0' for the empty set.
    internal static char CodeOf(int mask)
    {
        foreach (char code in IupacByPreference)
            if (MaskOf(code) == mask)
                return code;
        return '\0';
    }

    /// <summary>
    /// The pattern as read on the other strand: the elements in reverse order with every base set complemented (A↔T, C↔G) and the anchors swapped
    /// (<c>&lt;</c> becomes <c>&gt;</c>), which is what <c>fuzznuc -complement</c> searches — verified on the binary.
    /// </summary>
    public FuzznucPattern ReverseComplement()
    {
        var reversed = _elements.Reverse().Select(e => e with { BaseSet = ComplementMask(e.BaseSet) }).ToArray();
        return new FuzznucPattern(ToString(reversed, AnchoredAtEnd, AnchoredAtStart), AnchoredAtEnd, AnchoredAtStart, reversed);
    }

    private static int ComplementMask(int mask)
        => ((mask & 1) << 3) | ((mask & 8) >> 3) | ((mask & 2) << 1) | ((mask & 4) >> 1); // A(1)<->T(8), C(2)<->G(4)

    /// <summary>
    /// Every fixed-length IUPAC string the pattern stands for, one per distinct choice of the repeat counts (duplicates removed), ordered by length then ordinal. A hit of the
    /// pattern is a window equal to one of them (a base matches a letter when it is in the letter's base set). Empty strings (every optional element dropped) are not hits and are
    /// omitted; a pattern with an empty base set (<c>{ACGT}</c>) in a mandatory position expands to nothing.
    /// </summary>
    /// <exception cref="ArgumentException">The repeat ranges allow more than <see cref="MaxCombinations"/> combinations.</exception>
    public IReadOnlyList<string> Expand()
        => ExpandWithExtras().Select(e => e.Iupac).Distinct(StringComparer.Ordinal).ToList();

    /// <summary>
    /// <see cref="Expand"/> with, for each expansion, the positions that are <b>extra repeats</b> — the copies of an element <c>X(n,m)</c> beyond its first n. fuzznuc's
    /// mismatch search (EMBOSS <c>embpat.c</c> <c>patOUBrute</c>, used for every pattern that has a range when <c>-pmismatch</c> &gt; 0) lets a mismatch fall only on the mandatory
    /// part of the pattern; the extra repeats of a range must match exactly. The same string may appear with different extra-position sets (which element supplied the extras);
    /// each distinct (string, extras) pair is listed once, ordered by length, string, extras.
    /// </summary>
    /// <exception cref="ArgumentException">The repeat ranges allow more than <see cref="MaxCombinations"/> combinations.</exception>
    public IReadOnlyList<PatternExpansion> ExpandWithExtras()
    {
        double combinations = 1;
        foreach (var e in _elements)
            combinations *= e.MaxRepeat - e.MinRepeat + 1;
        if (combinations > MaxCombinations)
            throw new ArgumentException($"The pattern '{Source}' has {combinations:N0} repeat combinations (limit {MaxCombinations:N0}); narrow its ranges.");

        var results = new Dictionary<string, PatternExpansion>(StringComparer.Ordinal);
        Recurse(0, string.Empty, [], false);
        return results.Values.OrderBy(e => e.Iupac.Length).ThenBy(e => e.Iupac, StringComparer.Ordinal)
            .ThenBy(e => string.Join(',', e.ExtraPositions), StringComparer.Ordinal).ThenBy(e => e.EndsWithEmptyElement).ToList();

        void Recurse(int index, string prefix, List<int> extras, bool lastEmpty)
        {
            if (index == _elements.Length)
            {
                if (prefix.Length > 0)
                    results.TryAdd(prefix + "|" + string.Join(',', extras) + "|" + lastEmpty, new PatternExpansion(prefix, extras.ToArray(), lastEmpty));
                return;
            }

            var element = _elements[index];
            char code = CodeOf(element.BaseSet);
            for (int count = element.MinRepeat; count <= element.MaxRepeat; count++)
            {
                if (count > 0 && code == '\0')
                    break; // an element with no allowed base cannot be repeated into a hit
                var next = new List<int>(extras);
                for (int extra = element.MinRepeat; extra < count; extra++)
                    next.Add(prefix.Length + extra);
                Recurse(index + 1, code == '\0' ? prefix : prefix + new string(code, count), next, count == 0);
            }
        }
    }

    /// <summary>The normalised pattern text: upper case, classes and exclusions in A, C, G, T order, no separators.</summary>
    public override string ToString() => ToString(_elements, AnchoredAtStart, AnchoredAtEnd);

    private static string ToString(PatternElement[] elements, bool start, bool end)
    {
        var sb = new StringBuilder();
        if (start)
            sb.Append('<');
        foreach (var e in elements)
        {
            char code = CodeOf(e.BaseSet);
            string atom;
            if (e.BaseSet == 0)
                atom = "{ACGT}";
            else if (code is 'A' or 'C' or 'G' or 'T' or 'N')
                atom = code.ToString();
            else
                atom = "[" + string.Concat(Bases.Where((_, i) => (e.BaseSet >> i & 1) != 0)) + "]";
            sb.Append(atom);
            if (e.MinRepeat != 1 || e.MaxRepeat != 1)
                sb.Append(e.MinRepeat == e.MaxRepeat ? $"({e.MinRepeat})" : $"({e.MinRepeat},{e.MaxRepeat})");
        }

        if (end)
            sb.Append('>');
        return sb.ToString();
    }

    /// <summary>
    /// Parses a fuzznuc pattern file (<c>-pattern @file</c>): one pattern per line, or FASTA style (<c>&gt;name</c> then the pattern). A file whose first
    /// non-blank line starts with <c>&gt;</c> is read as FASTA: the name is the header up to the first blank, the pattern the following lines joined. Otherwise every non-blank
    /// line is a pattern. Unnamed patterns are called <c>pattern1</c>, <c>pattern2</c>, … by their 1-based position (fuzznuc's report names; with <c>-pname</c> the prefix is replaced).
    /// </summary>
    /// <param name="text">The file content.</param>
    /// <param name="defaultName">Prefix of the generated names (fuzznuc <c>-pname</c>, default <c>pattern</c>).</param>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    /// <exception cref="ArgumentException">A pattern is invalid (the message names its line).</exception>
    public static IReadOnlyList<NamedPattern> ParsePatternFile(string text, string defaultName = "pattern")
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(defaultName);
        var lines = text.Split('\n').Select(l => l.Trim()).Where(l => l.Length > 0).ToList();
        var result = new List<NamedPattern>();
        if (lines.Count == 0)
            return result;

        if (lines[0].StartsWith('>'))
        {
            string? name = null;
            var body = new StringBuilder();
            foreach (var line in lines.Append(">"))
            {
                if (line.StartsWith('>'))
                {
                    if (body.Length > 0 || name is not null)
                        result.Add(new NamedPattern(name ?? $"{defaultName}{result.Count + 1}", Parse(body.ToString())));
                    string header = line[1..].Trim();
                    name = header.Length == 0 ? null : string.Concat(header.TakeWhile(ch => !char.IsWhiteSpace(ch)));
                    body.Clear();
                }
                else
                {
                    body.Append(line);
                }
            }

            return result;
        }

        foreach (var line in lines)
            result.Add(new NamedPattern($"{defaultName}{result.Count + 1}", Parse(line)));
        return result;
    }
}

/// <summary>
/// One element of a <see cref="FuzznucPattern"/>: a set of allowed bases (bit 0 = A, bit 1 = C, bit 2 = G, bit 3 = T) repeated
/// <see cref="MinRepeat"/> … <see cref="MaxRepeat"/> times.
/// </summary>
/// <param name="BaseSet">Allowed bases as a 4-bit mask (A = 1, C = 2, G = 4, T = 8); 0 = no base allowed.</param>
/// <param name="MinRepeat">Fewest repeats (0 makes the element optional).</param>
/// <param name="MaxRepeat">Most repeats.</param>
/// <param name="IsRange">Written as a range <c>(n,m)</c> (even when n = m) — fuzznuc treats every pattern with such an element as a "range pattern" (regex / brute-force search).</param>
public readonly record struct PatternElement(int BaseSet, int MinRepeat, int MaxRepeat, bool IsRange = false)
{
    /// <summary>The allowed bases in A, C, G, T order.</summary>
    public string Bases
    {
        get
        {
            int set = BaseSet;
            return string.Concat("ACGT".Where((_, i) => (set >> i & 1) != 0));
        }
    }

    /// <summary>The IUPAC code with exactly this base set, or '\0' for the empty set.</summary>
    public char IupacCode => FuzznucPattern.CodeOf(BaseSet);
}

/// <summary>A fixed-length expansion of a <see cref="FuzznucPattern"/>: its IUPAC string and the positions of the extra repeats (see <see cref="FuzznucPattern.ExpandWithExtras"/>).</summary>
/// <param name="Iupac">The expansion as IUPAC letters.</param>
/// <param name="ExtraPositions">0-based positions within the expansion that are extra repeats of a range.</param>
/// <param name="EndsWithEmptyElement">The last pattern element contributes no base (an optional element, <c>X(0,m)</c>, chosen zero times).</param>
public readonly record struct PatternExpansion(string Iupac, IReadOnlyList<int> ExtraPositions, bool EndsWithEmptyElement = false);

/// <summary>A pattern with the name fuzznuc reports for it (<c>pattern1</c>, or the FASTA header).</summary>
/// <param name="Name">Pattern name.</param>
/// <param name="Pattern">The parsed pattern.</param>
public readonly record struct NamedPattern(string Name, FuzznucPattern Pattern);
