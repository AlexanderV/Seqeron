using System.Runtime.CompilerServices;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// Canonical case folding and T/U spelling of sequence text, shared by every layer
/// (Infrastructure, Core, MolTools, MCP servers).
/// </summary>
/// <remarks>
/// Folding is ASCII-only: unlike <see cref="char.ToUpperInvariant(char)"/> it never maps a non-ASCII
/// character onto an ASCII sequence symbol (U+017F 'ſ' → 'S', U+0131 'ı' → 'I'), so such a character
/// stays invalid instead of silently reading as an IUPAC code / amino acid (review 2026-09,
/// TRANS-PROT-001; EMBOSS pepstats skips it).
/// </remarks>
public static class SequenceCase
{
    /// <summary>Upper-cases an ASCII letter; every other character is returned unchanged.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static char ToUpperAscii(char c) => char.IsAsciiLetterLower(c) ? (char)(c - 32) : c;

    /// <summary>
    /// String form of <see cref="ToUpperAscii(char)"/>: upper-cases ASCII letters only and returns
    /// <paramref name="s"/> itself when it has no ASCII lower-case letter. Length-preserving.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="s"/> is null.</exception>
    public static string ToUpperAscii(string s)
    {
        ArgumentNullException.ThrowIfNull(s);
        if (s.AsSpan().IndexOfAnyInRange('a', 'z') < 0)
            return s;
        return string.Create(s.Length, s, static (dest, src) =>
        {
            for (int i = 0; i < src.Length; i++)
                dest[i] = ToUpperAscii(src[i]);
        });
    }

    /// <summary>
    /// Span form of <see cref="ToUpperAscii(char)"/>: writes the ASCII-upper-cased
    /// <paramref name="source"/> into <paramref name="destination"/> without allocating.
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="destination"/> is shorter than <paramref name="source"/>.</exception>
    public static void ToUpperAscii(ReadOnlySpan<char> source, Span<char> destination)
    {
        if (destination.Length < source.Length)
            throw new ArgumentException("Destination is shorter than source.", nameof(destination));
        for (int i = 0; i < source.Length; i++)
            destination[i] = ToUpperAscii(source[i]);
    }

    /// <summary>Upper-case DNA spelling: <see cref="ToUpperAscii(string)"/>, then U read as T.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="s"/> is null.</exception>
    public static string ToUpperDna(string s) => ToUpperAscii(s).Replace('U', 'T');

    /// <summary>Upper-case RNA spelling: <see cref="ToUpperAscii(string)"/>, then T read as U.</summary>
    /// <exception cref="ArgumentNullException"><paramref name="s"/> is null.</exception>
    public static string ToUpperRna(string s) => ToUpperAscii(s).Replace('T', 'U');
}
