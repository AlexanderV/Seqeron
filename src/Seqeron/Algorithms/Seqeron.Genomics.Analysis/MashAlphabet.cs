using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The alphabet of a Mash sketch (<c>mash sketch -a | -z TEXT | -Z</c>; <c>Sketch.cpp</c> <c>setAlphabetFromString</c>, <c>sketchParameterSetup.cpp</c>, marbl/Mash v2.3): the set of characters a k-mer may hold, the case rule and the
/// consequences for the sketch: a window with a character outside the alphabet is skipped; the k-mer hash is 64-bit when <c>alphabetSize^k &gt; 2^32</c> (<see cref="Uses64BitHash"/>); the k-mer space of the p-value is
/// <c>alphabetSize^k</c> (<see cref="KmerSpace"/>); <c>-a</c> and <c>-z</c> imply <c>-n</c> (<see cref="ImpliesNoncanonical"/>) and <c>-a</c> implies <c>-k 9</c> (<see cref="ProteinKmerSize"/>) unless a k is given.
/// A5-10b2 (B09 F214).
/// </summary>
public sealed class MashAlphabet : IEquatable<MashAlphabet>
{
    /// <summary><c>alphabetNucleotide</c> of <c>Sketch.h</c>.</summary>
    public const string NucleotideCharacters = "ACGT";

    /// <summary><c>alphabetProtein</c> of <c>Sketch.h</c> (A-Z except BJOUXZ).</summary>
    public const string ProteinCharacters = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>The k-mer size <c>-a</c> sets when <c>-k</c> is not given.</summary>
    public const int ProteinKmerSize = 9;

    private readonly bool[] _table = new bool[256];

    private MashAlphabet(string characters, bool preserveCase, bool impliesNoncanonical, bool foldTable = false)
    {
        Characters = characters;
        PreserveCase = preserveCase;
        ImpliesNoncanonical = impliesNoncanonical;
        foreach (char c in characters)
            _table[foldTable && c is >= 'a' and <= 'z' ? (char)(c - 32) : c] = true;
    }

    /// <summary>
    /// The alphabet as <c>mash screen</c> / <c>taxscreen</c> use it (<c>CommandScreen.cpp</c>): the table is built by <c>setAlphabetFromString</c> BEFORE <c>parameters.preserveCase</c> is set from the sketch, so its letters are upper-cased
    /// whatever the sketch's case flag; the mixture itself is upper-cased only when the sketch does not preserve the case (a <c>-z acgt -Z</c> sketch screens nothing: its lower-case letters are not in the table of the mixture).
    /// </summary>
    internal MashAlphabet ScreenView() => new(Characters, PreserveCase, ImpliesNoncanonical, true);

    /// <summary>The characters of the alphabet in ascending order (what <c>Sketch::getAlphabetAsString</c> returns and the sketch file stores).</summary>
    public string Characters { get; }

    /// <summary><c>-Z</c>: the case of the letters is kept (the alphabet and the k-mers are not upper-cased; a letter whose case is not in the alphabet is skipped).</summary>
    public bool PreserveCase { get; }

    /// <summary>True for the alphabets chosen by <c>-a</c> or <c>-z</c>, which imply <c>-n</c> (strand preserved).</summary>
    public bool ImpliesNoncanonical { get; }

    /// <summary>The number of characters (<c>alphabetSize</c>).</summary>
    public int Size => Characters.Length;

    /// <summary>The nucleotide alphabet <c>ACGT</c> (the default, canonical k-mers), with the case rule of <c>-Z</c> when <paramref name="preserveCase"/>.</summary>
    public static MashAlphabet Nucleotide(bool preserveCase = false) => new(NucleotideCharacters, preserveCase, false);

    /// <summary>The amino-acid alphabet of <c>-a</c> (implies <c>-n</c>).</summary>
    public static MashAlphabet Protein(bool preserveCase = false) => new(ProteinCharacters, preserveCase, true);

    /// <summary>
    /// <c>setAlphabetFromString</c> for <c>-z TEXT</c> (implies <c>-n</c>): every character of the text, upper-cased unless <paramref name="preserveCase"/>, once.
    /// </summary>
    /// <exception cref="ArgumentException">A character is not ASCII (the program indexes its alphabet table with the signed <c>char</c>: a byte of 0x80 or more writes outside the table into the neighbouring fields of its parameter structure and has no defined result).</exception>
    public static MashAlphabet FromString(string characters, bool preserveCase = false)
        => Build(characters, preserveCase, true);

    /// <summary>
    /// The alphabet the options of <c>mash sketch</c> / <c>dist</c> / <c>triangle</c> choose (<c>sketchParameterSetup.cpp</c>): <c>-a</c> (protein) takes precedence over <c>-z</c>, else the nucleotides; <c>-Z</c> is <paramref name="preserveCase"/>.
    /// </summary>
    public static MashAlphabet FromOptions(bool protein = false, string? characters = null, bool preserveCase = false)
    {
        if (protein)
            return Protein(preserveCase);
        if (characters is not null)
            return FromString(characters, preserveCase);
        return Nucleotide(preserveCase);
    }

    /// <summary>The alphabet a sketch file stores (<c>alphabet</c> text and <c>preserveCase</c>; a file without an alphabet is a nucleotide sketch): <c>noncanonical</c> is the file's own flag.</summary>
    public static MashAlphabet FromSketchFile(string? characters, bool preserveCase)
        => Build(characters ?? NucleotideCharacters, preserveCase, false);

    private static MashAlphabet Build(string text, bool preserveCase, bool implies)
    {
        ArgumentNullException.ThrowIfNull(text);
        var set = new SortedSet<char>();
        foreach (char c in text)
        {
            if (c > 127)
                throw new ArgumentException($"The character U+{(int)c:X4} is not ASCII: Mash indexes its alphabet table with a signed char, so a byte of 0x80 or more has no defined meaning.", nameof(text));
            char u = !preserveCase && c is >= 'a' and <= 'z' ? (char)(c - 32) : c;
            set.Add(u);
        }

        var sb = new StringBuilder();
        foreach (char c in set)
            sb.Append(c);
        return new MashAlphabet(sb.ToString(), preserveCase, implies);
    }

    /// <summary>True for the default: the nucleotides <c>ACGT</c> with the case ignored (the fast path of the sketch).</summary>
    public bool IsDefaultNucleotide => !PreserveCase && Characters == NucleotideCharacters;

    /// <summary>True when the mixture of <c>mash screen</c> / <c>taxscreen</c> is 6-frame translated (<c>alphabet == alphabetProtein</c>).</summary>
    public bool IsProtein => Characters == ProteinCharacters;

    /// <summary>Whether the character is in the alphabet (after the case rule was applied to it by the caller).</summary>
    public bool Contains(char c) => c < 256 && _table[c];

    /// <summary><c>pow(alphabetSize, k)</c>: the k-mer space of the p-value (<c>Sketch::kmerSpace</c>).</summary>
    public double KmerSpace(int k) => Math.Pow(Size, k);

    /// <summary><c>use64 = pow(alphabetSize, k) &gt; pow(2, 32)</c> (<c>setAlphabetFromString</c>).</summary>
    public bool Uses64BitHash(int k) => KmerSpace(k) > 4294967296.0;

    /// <summary>The k-mer size of the options: <paramref name="requested"/>, or 9 for <c>-a</c> and 21 otherwise when no <c>-k</c> was given.</summary>
    public static int ResolveKmerSize(int? requested, bool protein) => requested ?? (protein ? ProteinKmerSize : 21);

    /// <inheritdoc />
    public bool Equals(MashAlphabet? other) => other is not null && Characters == other.Characters && PreserveCase == other.PreserveCase;

    /// <inheritdoc />
    public override bool Equals(object? obj) => Equals(obj as MashAlphabet);

    /// <inheritdoc />
    public override int GetHashCode() => HashCode.Combine(Characters, PreserveCase);

    /// <inheritdoc />
    public override string ToString() => Characters;

    // ------------------------------------------------------------------------------------------------ hashing of the alphabet k-mers

    // The 256 bytes of the installed Debian mash 2.3 (x86-64) that start 193 bytes before the table `complement` of Sketch.cpp ("Array from 0..25 of DNA complement of A..Z"): reverseComplement() indexes it with
    // (int) src - 'A', which for a character other than A-Z (a '*' separator of the screen, a lower-case letter with -Z, a digit ...) is out of the table: the program reads the bytes around it (read-only data: the same in every run).
    // Index 193 + 0 is 'A'.
    private static ReadOnlySpan<byte> ComplementWithNeighbours =>
    [
        0x67, 0x9E, 0x49, 0x6C, 0x9A, 0x27, 0x02, 0xFB, 0xDC, 0x33, 0x34, 0xA4, 0x46, 0x7A, 0xB0, 0x32, 0x18, 0x3C, 0xFF, 0x6A, 0xC4, 0xF0, 0xFF, 0xC8, 0xB9, 0x67, 0x04, 0x50, 0xC9, 0xDF, 0x1B, 0x43,
        0x34, 0x2B, 0x43, 0x52, 0xB6, 0x2B, 0xF1, 0x67, 0xE3, 0x1B, 0x4C, 0x7F, 0x42, 0x6E, 0x00, 0x4C, 0x22, 0x89, 0x5A, 0x7E, 0x2E, 0xA5, 0x9A, 0xF9, 0x96, 0x2A, 0x45, 0xEB, 0x0B, 0x39, 0x7C, 0xD8,
        0x2F, 0xFB, 0xE1, 0xB2, 0x74, 0x24, 0xFD, 0x2E, 0x22, 0x0C, 0xF6, 0x57, 0xF3, 0x1E, 0xCB, 0x0F, 0x44, 0x0F, 0x03, 0xBE, 0x8B, 0x29, 0xDC, 0x04, 0x67, 0x2F, 0xD1, 0x44, 0x11, 0x55, 0x13, 0x8B,
        0x94, 0xE9, 0xD7, 0x8F, 0x6D, 0x14, 0xA0, 0x11, 0x1C, 0x2F, 0x59, 0xD1, 0xAD, 0x2E, 0x71, 0x3C, 0xFB, 0x2F, 0x64, 0x77, 0xFC, 0x8C, 0xCE, 0xC4, 0xF9, 0xB9, 0x2F, 0x31, 0x31, 0x79, 0xDD, 0xB0,
        0x08, 0xE7, 0xA6, 0x8F, 0x31, 0x3D, 0xD2, 0x40, 0xC0, 0xA7, 0x9A, 0x58, 0xC0, 0x75, 0xC0, 0xA5, 0x0C, 0x72, 0xB1, 0x74, 0xF8, 0xD5, 0x14, 0xF9, 0x0C, 0x80, 0x32, 0x4D, 0x78, 0xBC, 0xFE, 0x8C,
        0x4E, 0x75, 0xF5, 0x69, 0xC5, 0x91, 0xA0, 0xB2, 0xCD, 0xB4, 0x16, 0xC0, 0x2C, 0x21, 0x44, 0x5F, 0x5C, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x54, 0x56, 0x47, 0x48, 0x4E, 0x4E, 0x43, 0x44, 0x4E, 0x4E, 0x4D, 0x4E, 0x4B, 0x4E, 0x4E, 0x4E, 0x4E, 0x59, 0x53, 0x41, 0x41, 0x42, 0x57, 0x4E, 0x52, 0x4E, 0x00, 0x00, 0x56, 0xD1, 0xFE,
        0xFF, 0x4C, 0xD1, 0xFE, 0xFF, 0x41, 0xD1, 0xFE, 0xFF, 0x36, 0xD1, 0xFE, 0xFF, 0x2B, 0xD1, 0xFE, 0xFF, 0x20, 0xD1, 0xFE, 0xFF, 0x15, 0xD1, 0xFE, 0xFF, 0x0A, 0xD1, 0xFE, 0xFF, 0x1F, 0xE9, 0xFE,
    ];

    /// <summary><c>complement[(int) c - 'A']</c> of <c>Sketch.cpp</c> for the byte <paramref name="c"/> (a signed <c>char</c> in the program): the DNA complement of A-Z, any other value as the program reads it out of the table.</summary>
    internal static byte Complement(int c)
    {
        int signedChar = c > 127 ? c - 256 : c;
        return ComplementWithNeighbours[193 + signedChar - 'A'];
    }

    /// <summary>
    /// Streams the hash (h1 of MurmurHash3_x64_128 over the k-mer's bytes) of every window of <paramref name="seq"/> that holds only characters of the alphabet, as <c>addMinHashes</c> / <c>hashSequence</c> take them: the sequence is upper-cased
    /// unless <c>-Z</c> (<c>seq[i] &gt; 96 &amp;&amp; seq[i] &lt; 123</c>), a k-mer is the forward one, or with <paramref name="canonical"/> the smaller (<c>memcmp(fwd, rev) &lt;= 0</c>) of itself and its reverse complement.
    /// A character that is not ASCII is outside every alphabet here (the program's table lookup has no defined result for it).
    /// </summary>
    internal static void StreamHashes(string seq, int k, bool canonical, uint seed, MashAlphabet alphabet, Action<ulong> sink, CancellationToken cancellationToken)
    {
        int n = seq.Length;
        if (k > n)
            return;
        var s = new byte[n];
        var ok = new bool[n];
        for (int i = 0; i < n; i++)
        {
            char c = seq[i];
            if (c > 127)
            {
                s[i] = (byte)'?';
                continue;
            }

            if (!alphabet.PreserveCase && c is >= 'a' and <= 'z')
                c = (char)(c - 32);
            s[i] = (byte)c;
            ok[i] = alphabet._table[c];
        }

        var rev = canonical ? new byte[k] : null;
        int run = 0;
        for (int i = 0; i < n; i++)
        {
            if ((i & 4095) == 0)
                cancellationToken.ThrowIfCancellationRequested();
            if (!ok[i])
            {
                run = 0;
                continue;
            }

            if (++run < k)
                continue;
            run = k;
            var fwd = new ReadOnlySpan<byte>(s, i - k + 1, k);
            ReadOnlySpan<byte> kmer = fwd;
            if (canonical)
            {
                for (int j = 0; j < k; j++)
                    rev![j] = Complement(fwd[k - 1 - j]);
                if (rev.AsSpan().SequenceCompareTo(fwd) < 0)
                    kmer = rev;
            }

            sink(KmerAnalyzer.MurmurHash3X64_128(kmer, seed).H1);
        }
    }
}
