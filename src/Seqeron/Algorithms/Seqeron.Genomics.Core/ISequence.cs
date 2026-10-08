namespace Seqeron.Genomics.Core;

/// <summary>
/// Base interface for all biological sequences.
/// </summary>
public interface ISequence
{
    /// <summary>
    /// Gets the sequence data as a string.
    /// </summary>
    string Sequence { get; }

    /// <summary>
    /// Gets the length of the sequence.
    /// </summary>
    int Length { get; }

    /// <summary>
    /// Gets the sequence type.
    /// </summary>
    SequenceType Type { get; }

    /// <summary>
    /// Gets the valid alphabet for this sequence type.
    /// </summary>
    IReadOnlySet<char> Alphabet { get; }

    /// <summary>
    /// Gets the character at the specified position.
    /// </summary>
    char this[int index] { get; }

    /// <summary>
    /// Gets a subsequence.
    /// </summary>
    ISequence Subsequence(int start, int length);

    /// <summary>
    /// Validates the sequence against the alphabet.
    /// </summary>
    bool IsValid();

    /// <summary>
    /// Gets the complement of the sequence (for nucleotides).
    /// </summary>
    ISequence? GetComplement();

    /// <summary>
    /// Gets the reverse of the sequence.
    /// </summary>
    ISequence GetReverse();

    /// <summary>
    /// Gets the reverse complement (for nucleotides).
    /// </summary>
    ISequence? GetReverseComplement();
}

/// <summary>
/// Sequence type enumeration.
/// </summary>
public enum SequenceType
{
    Dna,
    Rna,
    Protein,
    IupacDna,
    IupacRna,
    Quality
}

/// <summary>
/// Base class for biological sequences with common functionality.
/// </summary>
public abstract class SequenceBase : ISequence
{
    protected readonly string _sequence;

    /// <remarks>
    /// Only ASCII letters are upper-cased: culture-invariant Unicode case mapping would turn the
    /// non-IUPAC character U+017F 'ſ' into 'S' and let it pass <see cref="IsValid"/>
    /// (scikit-bio <c>DNA("ſ", lowercase=True)</c> rejects it).
    /// </remarks>
    protected SequenceBase(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        _sequence = string.Create(sequence.Length, sequence, static (dest, src) =>
        {
            for (int i = 0; i < src.Length; i++)
                dest[i] = SequenceExtensions.ToUpperAscii(src[i]);
        });
    }

    public string Sequence => _sequence;
    public int Length => _sequence.Length;
    public abstract SequenceType Type { get; }
    public abstract IReadOnlySet<char> Alphabet { get; }

    public char this[int index] => _sequence[index];

    public abstract ISequence Subsequence(int start, int length);

    public virtual bool IsValid()
    {
        foreach (char c in _sequence)
        {
            if (!Alphabet.Contains(c))
                return false;
        }
        return true;
    }

    public abstract ISequence? GetComplement();
    public abstract ISequence GetReverse();
    public abstract ISequence? GetReverseComplement();

    public override string ToString() => _sequence;
    public override int GetHashCode() => _sequence.GetHashCode();
    public override bool Equals(object? obj) => obj is SequenceBase other && _sequence == other._sequence;
}

/// <summary>
/// IUPAC DNA sequence with ambiguity codes support.
/// Supports: A, C, G, T, N (any), R (A/G), Y (C/T), W (A/T), S (G/C), K (G/T), M (A/C),
/// B (C/G/T), D (A/G/T), H (A/C/T), V (A/C/G)
/// </summary>
public class IupacDnaSequence : SequenceBase
{
    private static readonly HashSet<char> _alphabet = new()
    {
        'A', 'C', 'G', 'T', 'U',  // Standard bases
        'N',                      // Any base
        'R', 'Y',                 // Purine (A/G), Pyrimidine (C/T)
        'W', 'S',                 // Weak (A/T), Strong (G/C)
        'K', 'M',                 // Keto (G/T), Amino (A/C)
        'B', 'D', 'H', 'V',       // 3-base ambiguity
        '-', '.'                  // Gaps
    };

    /// <summary>
    /// Complement of one IUPAC symbol: symbols in <see cref="_alphabet"/> are complemented by the
    /// canonical <see cref="SequenceExtensions.GetComplementBase(char)"/> (IUPAC NC-IUB 1984 table;
    /// gaps '-'/'.' map to themselves); any other character becomes 'N'.
    /// </summary>
    private static char ComplementSymbol(char c)
        => _alphabet.Contains(c) ? SequenceExtensions.GetComplementBase(c) : 'N';

    private static readonly Dictionary<char, char[]> _expansions = new()
    {
        ['A'] = new[] { 'A' },
        ['C'] = new[] { 'C' },
        ['G'] = new[] { 'G' },
        ['T'] = new[] { 'T' },
        ['U'] = new[] { 'T' },
        ['N'] = new[] { 'A', 'C', 'G', 'T' },
        ['R'] = new[] { 'A', 'G' },
        ['Y'] = new[] { 'C', 'T' },
        ['W'] = new[] { 'A', 'T' },
        ['S'] = new[] { 'G', 'C' },
        ['K'] = new[] { 'G', 'T' },
        ['M'] = new[] { 'A', 'C' },
        ['B'] = new[] { 'C', 'G', 'T' },
        ['D'] = new[] { 'A', 'G', 'T' },
        ['H'] = new[] { 'A', 'C', 'T' },
        ['V'] = new[] { 'A', 'C', 'G' }
    };

    public IupacDnaSequence(string sequence) : base(sequence) { }

    public override SequenceType Type => SequenceType.IupacDna;
    public override IReadOnlySet<char> Alphabet => _alphabet;

    public override ISequence Subsequence(int start, int length)
        => new IupacDnaSequence(_sequence.Substring(start, length));

    public override ISequence? GetComplement()
    {
        var complement = new char[Length];
        for (int i = 0; i < Length; i++)
        {
            complement[i] = ComplementSymbol(_sequence[i]);
        }
        return new IupacDnaSequence(new string(complement));
    }

    public override ISequence GetReverse()
    {
        var reversed = new char[Length];
        for (int i = 0; i < Length; i++)
            reversed[i] = _sequence[Length - 1 - i];
        return new IupacDnaSequence(new string(reversed));
    }

    public override ISequence? GetReverseComplement()
    {
        var result = new char[Length];
        for (int i = 0; i < Length; i++)
        {
            result[i] = ComplementSymbol(_sequence[Length - 1 - i]);
        }
        return new IupacDnaSequence(new string(result));
    }

    /// <summary>
    /// Expands IUPAC code to possible bases.
    /// </summary>
    public static char[] ExpandCode(char iupacCode)
    {
        return _expansions.TryGetValue(char.ToUpperInvariant(iupacCode), out var bases)
            ? bases
            : new[] { 'N' };
    }

    /// <summary>
    /// Gets IUPAC code from a set of bases.
    /// </summary>
    /// <remarks>
    /// Inverse of the NC-IUB (1984) degenerate map (Biopython <c>IUPACData.ambiguous_dna_values</c> /
    /// <c>ambiguous_rna_values</c>): ASCII letters are case-folded and RNA U is treated as T, consistent
    /// with <see cref="ExpandCode"/> ('U' → T), so {A, U} → 'W' (Biopython <c>ambiguous_rna_values['W'] = "AU"</c>).
    /// Any set containing a non-base symbol, or an empty set, yields 'N'.
    /// </remarks>
    public static char GetIupacCode(IEnumerable<char> bases)
    {
        var baseSet = new HashSet<char>(bases.Select(b =>
        {
            char u = SequenceExtensions.ToUpperAscii(b);
            return u == 'U' ? 'T' : u;
        }));

        if (baseSet.SetEquals(new[] { 'A' })) return 'A';
        if (baseSet.SetEquals(new[] { 'C' })) return 'C';
        if (baseSet.SetEquals(new[] { 'G' })) return 'G';
        if (baseSet.SetEquals(new[] { 'T' })) return 'T';
        if (baseSet.SetEquals(new[] { 'A', 'G' })) return 'R';
        if (baseSet.SetEquals(new[] { 'C', 'T' })) return 'Y';
        if (baseSet.SetEquals(new[] { 'A', 'T' })) return 'W';
        if (baseSet.SetEquals(new[] { 'G', 'C' })) return 'S';
        if (baseSet.SetEquals(new[] { 'G', 'T' })) return 'K';
        if (baseSet.SetEquals(new[] { 'A', 'C' })) return 'M';
        if (baseSet.SetEquals(new[] { 'C', 'G', 'T' })) return 'B';
        if (baseSet.SetEquals(new[] { 'A', 'G', 'T' })) return 'D';
        if (baseSet.SetEquals(new[] { 'A', 'C', 'T' })) return 'H';
        if (baseSet.SetEquals(new[] { 'A', 'C', 'G' })) return 'V';

        return 'N';
    }

    /// <summary>
    /// Checks if pattern matches sequence at position (with IUPAC wildcards).
    /// </summary>
    public bool MatchesAt(string pattern, int position)
    {
        if (position < 0 || position + pattern.Length > Length)
            return false;

        for (int i = 0; i < pattern.Length; i++)
        {
            char p = char.ToUpperInvariant(pattern[i]);
            char s = _sequence[position + i];

            if (!CodesMatch(p, s))
                return false;
        }

        return true;
    }

    /// <summary>
    /// Checks if two IUPAC codes can represent the same base.
    /// </summary>
    /// <remarks>
    /// True when the NC-IUB (1984) base sets of the two codes intersect (Biopython
    /// <c>IUPACData.ambiguous_dna_values</c>; scikit-bio <c>DNA.degenerate_map</c>). ASCII letters are
    /// case-folded first (as in <see cref="ExpandCode"/> and <see cref="MatchesAt"/>), so ('r', 'a') matches;
    /// a symbol outside the table only matches itself.
    /// </remarks>
    public static bool CodesMatch(char code1, char code2)
    {
        code1 = SequenceExtensions.ToUpperAscii(code1);
        code2 = SequenceExtensions.ToUpperAscii(code2);
        var bases1 = _expansions.TryGetValue(code1, out var b1) ? b1 : new[] { code1 };
        var bases2 = _expansions.TryGetValue(code2, out var b2) ? b2 : new[] { code2 };

        return bases1.Intersect(bases2).Any();
    }

    /// <summary>
    /// Finds all positions where pattern matches (with IUPAC support).
    /// </summary>
    public IEnumerable<int> FindPattern(string pattern)
    {
        for (int i = 0; i <= Length - pattern.Length; i++)
        {
            if (MatchesAt(pattern, i))
                yield return i;
        }
    }

    /// <summary>
    /// Calculates ambiguity level (1.0 = no ambiguity, 0.0 = all N).
    /// </summary>
    public double GetAmbiguityLevel()
    {
        if (Length == 0) return 1.0;

        int unambiguous = _sequence.Count(c => c == 'A' || c == 'C' || c == 'G' || c == 'T');
        return unambiguous / (double)Length;
    }

    /// <summary>
    /// Generates all possible concrete sequences from IUPAC sequence.
    /// Warning: exponential complexity for many ambiguous positions.
    /// </summary>
    public IEnumerable<string> ExpandAll(int maxResults = 1000)
    {
        var results = new List<string> { "" };

        foreach (char c in _sequence)
        {
            var expansions = ExpandCode(c);
            var newResults = new List<string>();

            foreach (var prefix in results)
            {
                foreach (var b in expansions)
                {
                    newResults.Add(prefix + b);
                    if (newResults.Count >= maxResults)
                    {
                        foreach (var r in newResults)
                            yield return r;
                        yield break;
                    }
                }
            }

            results = newResults;
        }

        foreach (var r in results)
            yield return r;
    }
}

/// <summary>
/// Sequence with per-base quality scores (FASTQ).
/// </summary>
public class QualitySequence : SequenceBase
{
    private static readonly HashSet<char> _alphabet = new()
    {
        'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n'
    };

    private readonly byte[] _qualities;

    public QualitySequence(string sequence, byte[] qualities) : base(sequence)
    {
        if (qualities.Length != sequence.Length)
            throw new ArgumentException("Quality array must match sequence length");

        _qualities = qualities;
    }

    /// <summary>
    /// Creates a quality sequence from an ASCII-encoded FASTQ quality string, Q = ord(c) − <paramref name="phredOffset"/>.
    /// </summary>
    /// <remarks>
    /// FASTQ (Cock et al. 2010, NAR 38:1767): the quality string has exactly one printable-ASCII character
    /// (33–126) per base; Sanger/Phred+33 covers Q0–Q93, Illumina 1.3+/Phred+64 Q0–Q62. As in Biopython
    /// <c>Bio.SeqIO.QualityIO</c>, a length mismatch or a character outside [offset, 126] is an error
    /// (previously the string was silently truncated/zero-padded and low characters clamped to Q0).
    /// </remarks>
    /// <exception cref="ArgumentException">Length mismatch or a character outside [<paramref name="phredOffset"/>, 126].</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="phredOffset"/> outside [33, 126].</exception>
    public QualitySequence(string sequence, string qualityString, int phredOffset = 33)
        : base(sequence)
    {
        ArgumentNullException.ThrowIfNull(qualityString);
        ValidatePhredOffset(phredOffset);
        if (qualityString.Length != sequence.Length)
            throw new ArgumentException(
                $"Quality string length ({qualityString.Length}) must match sequence length ({sequence.Length}).",
                nameof(qualityString));

        _qualities = new byte[sequence.Length];
        for (int i = 0; i < qualityString.Length; i++)
        {
            char c = qualityString[i];
            if (c < phredOffset || c > MaxQualityChar)
                throw new ArgumentException(
                    $"Invalid quality character '{c}' (0x{(int)c:X2}) at position {i}: must be in [{phredOffset}, {MaxQualityChar}].",
                    nameof(qualityString));
            _qualities[i] = (byte)(c - phredOffset);
        }
    }

    /// <summary>Highest printable-ASCII quality character ('~', 126) in any FASTQ variant (Cock et al. 2010).</summary>
    private const int MaxQualityChar = 126;

    /// <summary>Highest Phred score representable in Sanger/Phred+33 FASTQ (126 − 33).</summary>
    public const byte MaxSangerPhred = 93;

    private static void ValidatePhredOffset(int phredOffset)
    {
        if (phredOffset < 33 || phredOffset > MaxQualityChar)
            throw new ArgumentOutOfRangeException(nameof(phredOffset), phredOffset,
                "Phred offset must be a printable-ASCII code in [33, 126] (33 = Sanger, 64 = Illumina 1.3+).");
    }

    public override SequenceType Type => SequenceType.Quality;
    public override IReadOnlySet<char> Alphabet => _alphabet;

    /// <summary>
    /// Gets the quality scores.
    /// </summary>
    public IReadOnlyList<byte> Qualities => _qualities;

    /// <summary>
    /// Gets quality at position.
    /// </summary>
    public byte GetQuality(int index) => _qualities[index];

    /// <summary>
    /// Gets mean quality score.
    /// </summary>
    public double MeanQuality => _qualities.Average(q => (double)q);

    /// <summary>
    /// Gets the quality string (Phred+<paramref name="phredOffset"/> encoding, default Sanger Phred+33).
    /// </summary>
    /// <remarks>
    /// Scores above the encoding's ceiling (126 − offset: Q93 for Phred+33, Q62 for Phred+64) are capped at
    /// '~' so the output stays printable ASCII, as Biopython <c>_get_sanger_quality_str</c> /
    /// <c>_get_illumina_quality_str</c> do ("Data loss - max PHRED quality 93 in Sanger FASTQ").
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="phredOffset"/> outside [33, 126].</exception>
    public string GetQualityString(int phredOffset = 33)
    {
        ValidatePhredOffset(phredOffset);
        int maxQ = MaxQualityChar - phredOffset;
        var chars = new char[_qualities.Length];
        for (int i = 0; i < _qualities.Length; i++)
        {
            chars[i] = (char)(Math.Min(_qualities[i], maxQ) + phredOffset);
        }
        return new string(chars);
    }

    public override ISequence Subsequence(int start, int length)
    {
        var subQual = new byte[length];
        Array.Copy(_qualities, start, subQual, 0, length);
        return new QualitySequence(_sequence.Substring(start, length), subQual);
    }

    /// <summary>
    /// Complement of one read base via the canonical <see cref="SequenceExtensions.GetComplementBase(char)"/>
    /// (IUPAC NC-IUB 1984 table). Characters that are not IUPAC nucleotide symbols become 'N'.
    /// The self-complementary IUPAC codes S, W and N are genuine complements (Biopython
    /// <c>complement("ASWRN")</c> = "TSWYN"), so "unchanged" is not by itself treated as unknown.
    /// </summary>
    private static char ComplementReadBase(char c)
    {
        char comp = SequenceExtensions.GetComplementBase(c);
        return comp == c && c is not ('S' or 'W' or 'N') ? 'N' : comp;
    }

    public override ISequence? GetComplement()
    {
        var complement = new char[Length];
        for (int i = 0; i < Length; i++)
        {
            complement[i] = ComplementReadBase(_sequence[i]);
        }
        return new QualitySequence(new string(complement), _qualities);
    }

    public override ISequence GetReverse()
    {
        var reversed = new char[Length];
        var revQual = new byte[Length];
        for (int i = 0; i < Length; i++)
        {
            reversed[i] = _sequence[Length - 1 - i];
            revQual[i] = _qualities[Length - 1 - i];
        }
        return new QualitySequence(new string(reversed), revQual);
    }

    public override ISequence? GetReverseComplement()
    {
        var result = new char[Length];
        var revQual = new byte[Length];
        for (int i = 0; i < Length; i++)
        {
            char c = _sequence[Length - 1 - i];
            result[i] = ComplementReadBase(c);
            revQual[i] = _qualities[Length - 1 - i];
        }
        return new QualitySequence(new string(result), revQual);
    }

    /// <summary>
    /// Trims low quality bases from ends.
    /// </summary>
    public QualitySequence TrimByQuality(byte minQuality = 20)
    {
        int start = 0;
        int end = Length - 1;

        while (start < Length && _qualities[start] < minQuality)
            start++;

        while (end > start && _qualities[end] < minQuality)
            end--;

        if (start > end)
            return new QualitySequence("", Array.Empty<byte>());

        int len = end - start + 1;
        var trimmedQual = new byte[len];
        Array.Copy(_qualities, start, trimmedQual, 0, len);

        return new QualitySequence(_sequence.Substring(start, len), trimmedQual);
    }

    /// <summary>
    /// Masks low quality bases with N.
    /// </summary>
    public QualitySequence MaskLowQuality(byte minQuality = 20)
    {
        var masked = new char[Length];
        for (int i = 0; i < Length; i++)
        {
            masked[i] = _qualities[i] >= minQuality ? _sequence[i] : 'N';
        }
        return new QualitySequence(new string(masked), _qualities);
    }

    /// <summary>
    /// Gets error probability from Phred score.
    /// </summary>
    public static double PhredToErrorProbability(byte phred)
        => Math.Pow(10, -phred / 10.0);

    /// <summary>
    /// Gets Phred score from error probability: Q = round(−10·log10 p) (Ewing &amp; Green 1998; Cock et al. 2010),
    /// capped at Q93 (<see cref="MaxSangerPhred"/>, the Sanger FASTQ maximum); p = 0 → Q93.
    /// </summary>
    /// <remarks>
    /// Rounds to the nearest integer like Biopython <c>_get_sanger_quality_str</c> (29.99 → 30, 9.55 → 10):
    /// p = 0.2 → Q7, p = 0.0011 → Q30 (ties to even, as Python <c>round</c>). Previously the value was truncated (Q6, Q29).
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="errorProb"/> is NaN or outside [0, 1].</exception>
    public static byte ErrorProbabilityToPhred(double errorProb)
    {
        if (double.IsNaN(errorProb) || errorProb < 0 || errorProb > 1)
            throw new ArgumentOutOfRangeException(nameof(errorProb), errorProb, "Error probability must be in [0, 1].");
        if (errorProb == 0)
            return MaxSangerPhred;
        return (byte)Math.Min(MaxSangerPhred, Math.Round(-10 * Math.Log10(errorProb)));
    }

    /// <summary>
    /// Calculates expected number of errors.
    /// </summary>
    public double ExpectedErrors()
        => _qualities.Sum(q => PhredToErrorProbability(q));
}
