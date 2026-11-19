namespace SuffixTree.Esa;

/// <summary>
/// How much of a SESA container <see cref="SuffixIndexFactory.Open"/> checks before it hands the index out (ADR-003 section 6). Every level includes the cheaper ones.
/// </summary>
public enum EsaValidationLevel
{
    /// <summary>
    /// The default: the header page and the section table (magic, version, endian probe, header CRC-32C, file size, every section inside the file, aligned and not overlapping, the element counts and
    /// widths, the required sections and their agreement with the header, the alphabet and the symbol counts) and the layout of the packed text (65 bytes). O(sections), no pass over the payloads: a cold open
    /// takes milliseconds. A body that was damaged after the header was written is not noticed here; the <see cref="EsaAccessMode.Checked"/> mode keeps such a file memory safe.
    /// </summary>
    Header = 0,

    /// <summary>Also the CRC-32C of every section payload that carries one: one sequential pass over the file (hardware CRC where the CPU has it, a table path with the same value elsewhere).</summary>
    Sections = 1,

    /// <summary>
    /// Also the structure of the arrays, in one streaming O(n) pass with bounded working memory: the suffix array is a permutation of 0..n with <c>SA[0] = n</c> and its adjacent suffixes are in strictly ascending
    /// order, the LCP array equals the common prefixes of the adjacent suffixes, the child table equals the one the LCP array defines, the header statistics equal the recomputed ones, the Psi array (when present)
    /// equals <c>ISA[SA[i] + 1]</c>, every packed symbol is below sigma and the symbol counts equal a recount.
    /// </summary>
    Structure = 2,
}

/// <summary>Whether the reads of the opened index are checked against the bounds of their section (ADR-003 section 6).</summary>
public enum EsaAccessMode
{
    /// <summary>The default: every index that is read from the suffix array, the child table, the LCP array and its overflow list, and every symbol position, is validated against its section before it is used, and every walk is bounded; a body that does not make sense is an <see cref="System.IO.InvalidDataException"/>, never a read outside the mapping, a hang or a wrong answer that came from a read outside a section.</summary>
    Checked = 0,

    /// <summary>Opt-in for files that the process wrote itself or verified at the <see cref="EsaValidationLevel.Structure"/> level: the checks of the hot loops (the pattern descent, locate, the child walk) are left out. A corrupt file may crash the process or return garbage.</summary>
    Trusted = 1,
}

/// <summary>The options of <see cref="SuffixIndexFactory.Open"/> and <see cref="EsaSuffixTreeFactory.Open"/> for a SESA container; a v6 persistent file is opened as before and ignores them.</summary>
public sealed record EsaOpenOptions
{
    /// <summary>How much is validated; the default is <see cref="EsaValidationLevel.Header"/>.</summary>
    public EsaValidationLevel Validation { get; init; } = EsaValidationLevel.Header;

    /// <summary>How the reads are checked; the default is <see cref="EsaAccessMode.Checked"/>.</summary>
    public EsaAccessMode Access { get; init; } = EsaAccessMode.Checked;

    /// <summary>At the <see cref="EsaValidationLevel.Structure"/> level also compare <c>TEXT_SHA256</c> with the SHA-256 of the text (a header that records zeros is refused only when this is set and the header says a hash exists). Ignored at the other levels.</summary>
    public bool VerifyTextSha256 { get; init; }
}
