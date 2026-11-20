using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// A parsed EMBOSS <c>fuzzpro</c> / <c>fuzztran</c> protein pattern (the <c>-pattern</c> language, PROSITE-like): amino-acid letters, <c>x</c> (any residue), classes <c>[ILV]</c>,
/// exclusions <c>{P}</c>, repeats <c>x(3)</c> and length ranges <c>x(2,4)</c>, optional <c>-</c> separators, the anchors <c>&lt;</c> (N terminus) / <c>&gt;</c> (C terminus) and an optional
/// terminating period — e.g. <c>[DE](2)HS{P}x(2)Px(2,4)C</c>.
/// </summary>
/// <remarks>
/// <para>
/// Grammar and letter semantics as run on EMBOSS 6.6.0 (<c>fuzzpro -pattern</c>, 900 oracle rows): the letters are case-insensitive; <c>x</c> matches every residue; the ambiguity letters
/// <c>B</c> = {D, N}, <c>Z</c> = {E, Q} and <c>J</c> = {I, L} stand for their pairs (and never for themselves), every other letter (also <c>O</c>, <c>U</c>) only for itself;
/// a class / exclusion lists plain letters — <c>B</c>, <c>Z</c>, <c>J</c>, <c>X</c> inside <c>[]</c> / <c>{}</c> are rejected ("Nesting '[' in [] not allowed" in the binary), as is
/// any non-letter. A single repeat <c>(n)</c> needs n ≥ 1; a range <c>(n,m)</c> needs 0 ≤ n ≤ m and m ≥ 1 (<c>(0,2)</c> makes the element optional). The anchors must begin / end the
/// pattern. Malformed spellings that the binary accepts in unspecified ways (<c>[]</c>, <c>{}</c> matching everything, <c>A&lt;</c> reading as <c>A</c>) are rejected.
/// </para>
/// <para>
/// A pattern expands into a finite set of fixed-length element-set sequences, one per choice of the repeat counts; <see cref="GenomicAnalyzer.FindProteinPatternHits(string, string, CancellationToken)"/>
/// scans a sequence with each and unites the windows, which is what the binary reports (every matching window length of a range, overlapping hits included; a start-anchored range
/// pattern reports only its longest match, the regular-expression class of <c>embpat.c</c>, exactly as <see cref="FuzznucPattern"/>).
/// </para>
/// </remarks>
public sealed class ProteinPattern
{
    /// <summary>Largest number of repeat-count combinations an expansion may enumerate (a pattern with more is rejected).</summary>
    public const int MaxCombinations = 100_000;

    private const uint AnyResidue = 0xFFFFFFFFu;

    private readonly ProteinPatternElement[] _elements;

    private ProteinPattern(string source, bool anchoredAtStart, bool anchoredAtEnd, ProteinPatternElement[] elements)
    {
        Source = source;
        AnchoredAtStart = anchoredAtStart;
        AnchoredAtEnd = anchoredAtEnd;
        _elements = elements;
    }

    /// <summary>The pattern text as given (trimmed).</summary>
    public string Source { get; }

    /// <summary>The pattern starts with <c>&lt;</c>: a hit must start at the first residue.</summary>
    public bool AnchoredAtStart { get; }

    /// <summary>The pattern ends with <c>&gt;</c>: a hit must end at the last residue.</summary>
    public bool AnchoredAtEnd { get; }

    /// <summary>The elements in pattern order.</summary>
    public IReadOnlyList<ProteinPatternElement> Elements => _elements;

    /// <summary>True when an element is written as a range <c>(n,m)</c>.</summary>
    internal bool HasRange => _elements.Any(e => e.IsRange);

    /// <summary>Parses a pattern.</summary>
    /// <param name="pattern">The pattern text.</param>
    /// <exception cref="ArgumentNullException"><paramref name="pattern"/> is null.</exception>
    /// <exception cref="ArgumentException">The pattern is empty or violates the grammar; the message names the position.</exception>
    public static ProteinPattern Parse(string pattern)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        return ParseCore(pattern, out string? error) ?? throw new ArgumentException(error, nameof(pattern));
    }

    /// <summary>Parses a pattern without throwing.</summary>
    /// <param name="pattern">The pattern text.</param>
    /// <param name="result">The parsed pattern, or null.</param>
    /// <param name="error">The reason, or null on success.</param>
    /// <returns>True when the pattern is valid.</returns>
    public static bool TryParse(string? pattern, out ProteinPattern? result, out string? error)
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

    private static ProteinPattern? ParseCore(string pattern, out string? error)
    {
        error = null;
        string text = pattern.Trim();
        int i = 0;
        bool start = false, end = false;
        var elements = new List<ProteinPatternElement>();
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
                    error = $"Illegal pattern '{pattern}': '>' (C terminus) must end the pattern (position {i}).";
                    return null;
                }

                end = true;
                break;
            }

            if (c == '<')
            {
                error = $"Illegal pattern '{pattern}': '<' (N terminus) must begin the pattern (position {i}).";
                return null;
            }

            uint set;
            if (c is '[' or '{')
            {
                char close = c == '[' ? ']' : '}';
                int closeAt = text.IndexOf(close, i + 1);
                if (closeAt < 0)
                {
                    error = $"Illegal pattern '{pattern}': missing '{close}' for the '{c}' at position {i}.";
                    return null;
                }

                if (closeAt == i + 1)
                {
                    error = $"Illegal pattern '{pattern}': empty '{c}{close}' at position {i}.";
                    return null;
                }

                uint listed = 0;
                for (int j = i + 1; j < closeAt; j++)
                {
                    char u = char.ToUpperInvariant(text[j]);
                    if (u is < 'A' or > 'Z')
                    {
                        error = $"Illegal pattern '{pattern}': non-letter '{text[j]}' at position {j} inside '{c}{close}'.";
                        return null;
                    }

                    if (u is 'B' or 'Z' or 'J' or 'X')
                    {
                        error = $"Illegal pattern '{pattern}': the ambiguity letter '{text[j]}' at position {j} cannot be nested inside '{c}{close}'.";
                        return null;
                    }

                    listed |= Bit(u);
                }

                set = c == '[' ? listed : AnyResidue & ~listed;
                i = closeAt + 1;
            }
            else
            {
                char u = char.ToUpperInvariant(c);
                if (u is < 'A' or > 'Z')
                {
                    error = $"Illegal pattern '{pattern}': '{c}' at position {i} is not an amino-acid letter, '[', '{{', '<', '>' or '-'.";
                    return null;
                }

                set = u switch
                {
                    'X' => AnyResidue,
                    'B' => Bit('D') | Bit('N'),
                    'Z' => Bit('E') | Bit('Q'),
                    'J' => Bit('I') | Bit('L'),
                    _ => Bit(u),
                };
                i++;
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

            elements.Add(new ProteinPatternElement(set, min, max, isRange));
        }

        if (elements.Count == 0)
        {
            error = $"Illegal pattern '{pattern}': it has no element.";
            return null;
        }

        return new ProteinPattern(text, start, end, elements.ToArray());

        static bool TryCount(string s, out int value)
            => int.TryParse(s.Trim(), System.Globalization.NumberStyles.None, System.Globalization.CultureInfo.InvariantCulture, out value) && value <= 10_000;
    }

    private static uint Bit(char upperLetter) => 1u << (upperLetter - 'A');

    /// <summary>Bit of a residue in an element set: bit 0..25 for the letters A..Z, bit 31 for every other character (a stop <c>*</c>, a gap).</summary>
    internal static uint ResidueBit(char residue)
    {
        char u = char.ToUpperInvariant(residue);
        return u is >= 'A' and <= 'Z' ? Bit(u) : 1u << 31;
    }

    /// <summary>Every fixed-length sequence of element sets the pattern stands for (one per choice of the repeat counts, duplicates removed); empty sequences are omitted.</summary>
    /// <exception cref="ArgumentException">The repeat ranges allow more than <see cref="MaxCombinations"/> combinations.</exception>
    internal IReadOnlyList<uint[]> Expand()
    {
        double combinations = 1;
        foreach (var e in _elements)
            combinations *= e.MaxRepeat - e.MinRepeat + 1;
        if (combinations > MaxCombinations)
            throw new ArgumentException($"The pattern '{Source}' has {combinations:N0} repeat combinations (limit {MaxCombinations:N0}); narrow its ranges.");

        var results = new Dictionary<string, uint[]>(StringComparer.Ordinal);
        Recurse(0, new List<uint>());
        return results.Values.ToList();

        void Recurse(int index, List<uint> prefix)
        {
            if (index == _elements.Length)
            {
                if (prefix.Count > 0)
                    results.TryAdd(string.Join(",", prefix), prefix.ToArray());
                return;
            }

            var element = _elements[index];
            for (int count = element.MinRepeat; count <= element.MaxRepeat; count++)
            {
                var next = new List<uint>(prefix);
                for (int r = 0; r < count; r++)
                    next.Add(element.Set);
                Recurse(index + 1, next);
            }
        }
    }

    /// <summary>The normalised pattern text: upper case, classes and exclusions in alphabetical order, no separators.</summary>
    public override string ToString()
    {
        var sb = new StringBuilder();
        if (AnchoredAtStart) sb.Append('<');
        foreach (var e in _elements)
            sb.Append(e);
        if (AnchoredAtEnd) sb.Append('>');
        return sb.ToString();
    }
}

/// <summary>
/// One element of a <see cref="ProteinPattern"/>: a set of residues repeated between <see cref="MinRepeat"/> and <see cref="MaxRepeat"/> times.
/// </summary>
/// <param name="Set">Residue bit set (bit 0..25 = A..Z, bit 31 = any other character); <c>x</c> and an exclusion include bit 31.</param>
/// <param name="MinRepeat">Fewest repeats (0 for an optional element).</param>
/// <param name="MaxRepeat">Most repeats.</param>
/// <param name="IsRange">Written as a range <c>(n,m)</c> (even when n = m) — the binary then uses its regular-expression search class.</param>
public readonly record struct ProteinPatternElement(uint Set, int MinRepeat, int MaxRepeat, bool IsRange = false)
{
    /// <summary>The letters of the set, A..Z order.</summary>
    public string Letters
    {
        get
        {
            var sb = new StringBuilder();
            for (int i = 0; i < 26; i++)
                if ((Set >> i & 1) != 0)
                    sb.Append((char)('A' + i));
            return sb.ToString();
        }
    }

    /// <inheritdoc />
    public override string ToString()
    {
        string atom = Atom();
        if (MinRepeat == 1 && MaxRepeat == 1 && !IsRange)
            return atom;
        return IsRange ? $"{atom}({MinRepeat},{MaxRepeat})" : $"{atom}({MinRepeat})";
    }

    private string Atom()
    {
        if (Set == 0xFFFFFFFFu)
            return "x";
        string letters = Letters;
        if (letters.Length == 1 && (Set >> 31 & 1) == 0)
            return letters;
        return (Set >> 31 & 1) != 0 ? "{" + Complement() + "}" : "[" + letters + "]";
    }

    private string Complement()
    {
        var sb = new StringBuilder();
        for (int i = 0; i < 26; i++)
            if ((Set >> i & 1) == 0)
                sb.Append((char)('A' + i));
        return sb.ToString();
    }
}
