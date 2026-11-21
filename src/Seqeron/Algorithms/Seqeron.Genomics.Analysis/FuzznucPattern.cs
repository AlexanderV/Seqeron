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
/// fixed-length IUPAC strings — one per choice of the repeat counts — which the canonical matchers (<see cref="Seqeron.Genomics.Alignment.ApproximateMatcher"/> /
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
            var kind = PatternElementKind.Letter;
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
                kind = c == '[' ? PatternElementKind.Class : PatternElementKind.Exclusion;
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

            elements.Add(new PatternElement(mask, min, max, isRange, kind));
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
        var results = new Dictionary<string, PatternExpansion>(StringComparer.Ordinal);
        foreach (var path in ExpandPaths())
        {
            // an element with no allowed base cannot be repeated into a hit
            if (path.Elements.Any(index => _elements[index].BaseSet == 0))
                continue;
            string iupac = string.Concat(path.Elements.Select(index => CodeOf(_elements[index].BaseSet)));
            results.TryAdd(iupac + "|" + string.Join(',', path.ExtraPositions) + "|" + path.EndsWithEmptyElement,
                new PatternExpansion(iupac, path.ExtraPositions, path.EndsWithEmptyElement));
        }

        return results.Values.OrderBy(e => e.Iupac.Length).ThenBy(e => e.Iupac, StringComparer.Ordinal)
            .ThenBy(e => string.Join(',', e.ExtraPositions), StringComparer.Ordinal).ThenBy(e => e.EndsWithEmptyElement).ToList();
    }

    /// <summary>
    /// Every way the pattern can be unrolled into a fixed-length hit, one entry per choice of the repeat counts, as the <b>element</b> (index into <see cref="Elements"/>) that stands at each
    /// position — unlike <see cref="ExpandWithExtras"/> this keeps the kind of the element (letter, class, exclusion), which the literal sequence reading
    /// (<see cref="Matches(int, char)"/>) needs. Paths without a base are omitted; an element with an empty base set (<c>{ACGT}</c>) can still stand in a path.
    /// </summary>
    /// <param name="includeEmpty">Also the path of a pattern whose elements are ALL optional with every count 0 (the empty hit that <c>fuzznuc</c> lists as <c>end = start - 1</c>, A5-8c3b F217); default false.</param>
    /// <exception cref="ArgumentException">The repeat ranges allow more than <see cref="MaxCombinations"/> combinations.</exception>
    public IReadOnlyList<PatternPath> ExpandPaths(bool includeEmpty = false)
    {
        double combinations = 1;
        foreach (var e in _elements)
            combinations *= e.MaxRepeat - e.MinRepeat + 1;
        if (combinations > MaxCombinations)
            throw new ArgumentException(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"The pattern '{Source}' has {combinations:N0} repeat combinations (limit {MaxCombinations:N0}); narrow its ranges."));

        var paths = new List<PatternPath>();
        Recurse(0, [], [], false);
        return paths;

        void Recurse(int index, List<int> positions, List<int> extras, bool lastEmpty)
        {
            if (index == _elements.Length)
            {
                if (positions.Count > 0 || includeEmpty)
                    paths.Add(new PatternPath(positions.ToArray(), extras.ToArray(), lastEmpty));
                return;
            }

            var element = _elements[index];
            for (int count = element.MinRepeat; count <= element.MaxRepeat; count++)
            {
                var nextPositions = new List<int>(positions);
                var nextExtras = new List<int>(extras);
                for (int copy = 0; copy < count; copy++)
                {
                    if (copy >= element.MinRepeat)
                        nextExtras.Add(nextPositions.Count);
                    nextPositions.Add(index);
                }

                Recurse(index + 1, nextPositions, nextExtras, count == 0);
            }
        }
    }

    /// <summary>
    /// Literal reading of a sequence symbol by element <paramref name="elementIndex"/> (<c>fuzznuc</c> with an ambiguity code in the SEQUENCE, measured on 6.6.0, 1 350 oracle rows):
    /// <list type="bullet">
    /// <item>a pattern IUPAC <b>letter</b> matches a symbol whose base set lies within the letter's own set (R matches A, G and R; the sequence's R does not match the letters A, M or Y; <c>N</c> matches every symbol, X included);</item>
    /// <item>a <b>class</b> <c>[ACG]</c> matches the plain bases listed (a sequence ambiguity code never matches a class);</item>
    /// <item>an <b>exclusion</b> <c>{AG}</c> matches every symbol not listed (an ambiguity code and X included; <c>{ACGT}</c> matches exactly the non-A/C/G/T symbols).</item>
    /// </list>
    /// U reads as T; a symbol that is not an IUPAC letter (X) matches only <c>N</c> and an exclusion.
    /// </summary>
    /// <param name="elementIndex">Index into <see cref="Elements"/>.</param>
    /// <param name="symbol">The sequence symbol (any case).</param>
    public bool Matches(int elementIndex, char symbol)
    {
        var e = _elements[elementIndex];
        char upper = char.ToUpperInvariant(symbol);
        if (upper == 'U')
            upper = 'T';
        bool iupac = IupacByPreference.Contains(upper);
        int bits = iupac ? MaskOf(upper) : 0;
        switch (e.Kind)
        {
            case PatternElementKind.Class:
                return Bases.Contains(upper) && (bits & e.BaseSet) != 0;
            case PatternElementKind.Exclusion:
                {
                    bool plain = Bases.Contains(upper);
                    int listed = 0xF & ~e.BaseSet;
                    return !(plain && (bits & listed) != 0);
                }
            default:
                if (CodeOf(e.BaseSet) == 'N')
                    return true;
                return bits != 0 && (bits & ~e.BaseSet) == 0;
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
    /// Parses a fuzznuc pattern file (<c>-pattern @file</c>) with <see cref="EmbossPatternList.Read"/> (the program's <c>ajPatlistSeqRead</c>): one pattern per line, or FASTA style (a <c>&gt;name</c> header, optionally
    /// with <c>&lt;mismatch=n&gt;</c>, then the pattern over one or several lines); the first line decides unless the format is given. The name of a FASTA header is the whole header text; an unnamed pattern is called
    /// <c>pattern&lt;n&gt;</c> (the lines of a simple file: <c>-pname</c> or <c>pattern</c> followed by the line number, blank lines counted); each pattern carries its mismatches.
    /// </summary>
    /// <param name="text">The file content.</param>
    /// <param name="defaultName"><c>-pname</c> (null = not given).</param>
    /// <param name="mismatches"><c>-pmismatch</c>: the mismatches of a pattern whose header has no <c>&lt;mismatch=n&gt;</c>.</param>
    /// <param name="format"><c>-pformat</c>: <c>simple</c> or <c>fasta</c> (null = decided by the first line).</param>
    /// <exception cref="ArgumentNullException"><paramref name="text"/> is null.</exception>
    /// <exception cref="ArgumentException">A pattern is invalid (the message names its line).</exception>
    public static IReadOnlyList<NamedPattern> ParsePatternFile(string text, string? defaultName = null, int mismatches = 0, string? format = null)
    {
        ArgumentNullException.ThrowIfNull(text);
        var result = new List<NamedPattern>();
        foreach (var p in EmbossPatternList.Read("@file", defaultName, format, mismatches, _ => text))
        {
            string trimmed = p.Pattern.Trim();
            if (trimmed.Length > 0)
                result.Add(new NamedPattern(p.Name, Parse(trimmed), p.Mismatch));
        }

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
/// <param name="Kind">How the element was written: an IUPAC letter, a class <c>[..]</c> or an exclusion <c>{..}</c> (the three read an ambiguity code in the SEQUENCE differently, see <see cref="FuzznucPattern.Matches(int, char)"/>; for a plain A/C/G/T sequence they are the same set).
/// For an exclusion <paramref name="BaseSet"/> holds the bases that are allowed (the complement of the listed ones).</param>
public readonly record struct PatternElement(int BaseSet, int MinRepeat, int MaxRepeat, bool IsRange = false, PatternElementKind Kind = PatternElementKind.Letter)
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

/// <summary>How a <see cref="PatternElement"/> was written in the pattern.</summary>
public enum PatternElementKind
{
    /// <summary>An IUPAC nucleotide letter (<c>A</c>, <c>R</c>, <c>N</c>, ...).</summary>
    Letter = 0,

    /// <summary>A class <c>[ACG]</c> of plain bases.</summary>
    Class = 1,

    /// <summary>An exclusion <c>{AG}</c> of plain bases.</summary>
    Exclusion = 2,
}

/// <summary>One unrolling of a <see cref="FuzznucPattern"/> (<see cref="FuzznucPattern.ExpandPaths"/>): the element that stands at each position, the positions that are extra repeats, and whether the last element chose zero repeats.</summary>
/// <param name="Elements">Index into <see cref="FuzznucPattern.Elements"/> for every position of the hit.</param>
/// <param name="ExtraPositions">0-based positions that are extra repeats of a range (see <see cref="FuzznucPattern.ExpandWithExtras"/>).</param>
/// <param name="EndsWithEmptyElement">The last pattern element contributes no base.</param>
public readonly record struct PatternPath(IReadOnlyList<int> Elements, IReadOnlyList<int> ExtraPositions, bool EndsWithEmptyElement);

/// <summary>A fixed-length expansion of a <see cref="FuzznucPattern"/>: its IUPAC string and the positions of the extra repeats (see <see cref="FuzznucPattern.ExpandWithExtras"/>).</summary>
/// <param name="Iupac">The expansion as IUPAC letters.</param>
/// <param name="ExtraPositions">0-based positions within the expansion that are extra repeats of a range.</param>
/// <param name="EndsWithEmptyElement">The last pattern element contributes no base (an optional element, <c>X(0,m)</c>, chosen zero times).</param>
public readonly record struct PatternExpansion(string Iupac, IReadOnlyList<int> ExtraPositions, bool EndsWithEmptyElement = false);

/// <summary>A pattern with the name fuzznuc reports for it (<c>pattern1</c>, or the FASTA header).</summary>
/// <param name="Name">Pattern name.</param>
/// <param name="Pattern">The parsed pattern.</param>
/// <param name="Mismatch">The mismatches this pattern allows (a pattern file header's <c>&lt;mismatch=n&gt;</c> or <c>-pmismatch</c>); null = the caller's default.</param>
public readonly record struct NamedPattern(string Name, FuzznucPattern Pattern, int? Mismatch = null);
