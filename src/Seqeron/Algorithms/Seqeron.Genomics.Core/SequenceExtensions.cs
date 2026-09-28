using System.Runtime.CompilerServices;

namespace Seqeron.Genomics.Core;

/// <summary>
/// High-performance extension methods for sequence operations.
/// Provides Span-based overloads and CancellationToken support.
/// </summary>
public static class SequenceExtensions
{
    #region Span-based GC Content

    /// <summary>
    /// Calculates GC content as a percentage (0-100).
    /// Formula: (G + C) / (A + T + G + C + U) × 100
    /// Non-nucleotide characters are excluded from both numerator and denominator.
    /// Returns 0 for empty sequences or sequences with no valid nucleotides.
    /// </summary>
    /// <remarks>
    /// Matches the Wikipedia GC-content formula. Equals Biopython <c>gc_fraction(seq, "remove")</c> for any
    /// sequence without the IUPAC codes S/W; unlike Biopython, S (G|C) and W (A|T) are excluded here
    /// like every other ambiguity code — use <see cref="CalculateGcFraction(ReadOnlySpan{char},GcAmbiguityMode)"/>
    /// with <see cref="GcAmbiguityMode.Remove"/> for exact Biopython parity.
    /// Valid nucleotides: A, T, G, C, U (case-insensitive).
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double CalculateGcContent(this ReadOnlySpan<char> sequence)
    {
        // Single source of truth: percentage = fraction × 100 (bit-identical to (gc / valid) × 100).
        return sequence.CalculateGcFraction() * 100;
    }

    /// <summary>
    /// Counts the G/C nucleotides (numerator) and the valid A/T/G/C/U nucleotides (denominator)
    /// used by <see cref="CalculateGcFraction(ReadOnlySpan{char})"/> and
    /// <see cref="CalculateGcContent(ReadOnlySpan{char})"/> (case-insensitive; every other
    /// character, including IUPAC ambiguity codes, gaps and N, is excluded from both counts).
    /// Canonical counting primitive for callers that must report the raw counts.
    /// </summary>
    /// <param name="sequence">Nucleotide sequence.</param>
    /// <returns>(<c>GcCount</c> = #G + #C, <c>ValidCount</c> = #A + #T + #G + #C + #U).</returns>
    public static (int GcCount, int ValidCount) CountGcAndValidNucleotides(this ReadOnlySpan<char> sequence)
    {
        int gcCount = 0;
        int validCount = 0;
        for (int i = 0; i < sequence.Length; i++)
        {
            switch (sequence[i])
            {
                case 'G':
                case 'g':
                case 'C':
                case 'c':
                    gcCount++;
                    validCount++;
                    break;
                case 'A':
                case 'a':
                case 'T':
                case 't':
                case 'U':
                case 'u':
                    validCount++;
                    break;
            }
        }

        return (gcCount, validCount);
    }

    /// <summary>
    /// Calculates GC content for a string using Span optimization.
    /// Returns percentage (0-100).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double CalculateGcContentFast(this string sequence)
    {
        return sequence.AsSpan().CalculateGcContent();
    }

    /// <summary>
    /// Calculates GC content as a fraction (0-1).
    /// Formula: (G + C) / (A + T + G + C + U)
    /// Non-nucleotide characters are excluded from both numerator and denominator.
    /// Returns 0 for empty sequences or sequences with no valid nucleotides.
    /// </summary>
    /// <remarks>
    /// Matches the Wikipedia GC-content formula. Equals Biopython <c>gc_fraction(seq, "remove")</c> for any
    /// sequence without the IUPAC codes S/W; unlike Biopython, S (G|C) and W (A|T) are excluded here
    /// like every other ambiguity code — use <see cref="CalculateGcFraction(ReadOnlySpan{char},GcAmbiguityMode)"/>
    /// with <see cref="GcAmbiguityMode.Remove"/> for exact Biopython parity.
    /// Valid nucleotides: A, T, G, C, U (case-insensitive).
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double CalculateGcFraction(this ReadOnlySpan<char> sequence)
    {
        var (gcCount, validCount) = sequence.CountGcAndValidNucleotides();
        return validCount == 0 ? 0 : (double)gcCount / validCount;
    }

    /// <summary>
    /// Calculates GC content as a fraction (0-1) for a string using Span optimization.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double CalculateGcFractionFast(this string sequence)
    {
        return sequence.AsSpan().CalculateGcFraction();
    }

    #endregion

    #region Biopython-compatible GC fraction (IUPAC ambiguity modes)

    /// <summary>
    /// Strategy for handling IUPAC ambiguity codes when computing GC content, mirroring the
    /// <c>ambiguous</c> parameter of Biopython <c>Bio.SeqUtils.gc_fraction</c>.
    /// </summary>
    /// <remarks>
    /// Source: Biopython <c>Bio/SeqUtils/__init__.py</c> (master) <c>gc_fraction(seq, ambiguous="remove")</c>,
    /// retrieved 2026-06-23 from
    /// https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py .
    /// "Ambiguous nucleotides in this context are those different from ATCGSWU (S is G or C, and W is A or T)."
    /// </remarks>
    public enum GcAmbiguityMode
    {
        /// <summary>
        /// Biopython default. Numerator counts C, G and S; denominator counts only A, T, G, C, S, W, U.
        /// All other ambiguity codes (B,D,H,K,M,N,R,V,X,Y) are excluded from both numerator and denominator.
        /// </summary>
        Remove,

        /// <summary>
        /// Numerator counts only C, G and S; denominator is the FULL sequence length (every character),
        /// matching Biopython <c>ambiguous="ignore"</c>.
        /// </summary>
        Ignore,

        /// <summary>
        /// Ambiguity codes contribute their mean GC value (e.g. N/X = 0.5, V/B = 2/3, H/D = 1/3) to the
        /// numerator; denominator is the full sequence length. Matches Biopython <c>ambiguous="weighted"</c>.
        /// </summary>
        Weighted
    }

    // Weighted GC contribution of each IUPAC code, verbatim from Biopython _gc_values
    // (Bio/SeqUtils/__init__.py, master; retrieved 2026-06-23). Used only by GcAmbiguityMode.Weighted.
    // S=1 (G|C), W=0 (A|T); M/R/Y/K/X/N=0.5; V/B=2/3; H/D=1/3.
    private const double WeightedV = 2.0 / 3.0;
    private const double WeightedB = 2.0 / 3.0;
    private const double WeightedH = 1.0 / 3.0;
    private const double WeightedD = 1.0 / 3.0;
    private const double WeightedHalf = 0.5; // M, R, Y, K, X, N

    /// <summary>
    /// Calculates GC content as a fraction (0-1) using the IUPAC-ambiguity handling convention of
    /// Biopython <c>Bio.SeqUtils.gc_fraction(seq, ambiguous=...)</c>. This is an opt-in compatibility
    /// overload; the parameterless <see cref="CalculateGcFraction(ReadOnlySpan{char})"/> keeps the
    /// existing default behaviour (A/T/G/C/U only).
    /// </summary>
    /// <param name="sequence">Nucleotide sequence (case-insensitive).</param>
    /// <param name="mode">Ambiguity-handling mode; see <see cref="GcAmbiguityMode"/>.</param>
    /// <returns>GC fraction in [0,1]; 0 for an empty sequence or zero-length denominator.</returns>
    /// <remarks>
    /// Source: Biopython <c>gc_fraction</c> (master), retrieved 2026-06-23. For <c>Remove</c>:
    /// <c>gc = count(C,G,S)</c>, <c>length = gc + count(A,T,W,U)</c>. For <c>Ignore</c>:
    /// <c>gc = count(C,G,S)</c>, <c>length = len(seq)</c>. For <c>Weighted</c>: <c>gc</c> additionally
    /// adds <c>count(x)·_gc_values[x]</c> for x in BDHKMNRVXY, <c>length = len(seq)</c>.
    /// </remarks>
    public static double CalculateGcFraction(this ReadOnlySpan<char> sequence, GcAmbiguityMode mode)
    {
        if (sequence.IsEmpty) return 0;

        double gc = 0;        // numerator (may be fractional in Weighted mode)
        int strongWeakCount = 0; // A,T,G,C,S,W,U — the denominator for Remove mode
        int totalLength = sequence.Length; // denominator for Ignore/Weighted modes

        for (int i = 0; i < sequence.Length; i++)
        {
            switch (char.ToUpperInvariant(sequence[i]))
            {
                case 'G':
                case 'C':
                case 'S': // S = strong (G|C) -> full GC per Biopython "CGScgs"
                    gc += 1.0;
                    strongWeakCount++;
                    break;
                case 'A':
                case 'T':
                case 'U':
                case 'W': // W = weak (A|T) -> denominator only
                    strongWeakCount++;
                    break;
                default:
                    if (mode == GcAmbiguityMode.Weighted)
                        gc += WeightedGcValue(char.ToUpperInvariant(sequence[i]));
                    break;
            }
        }

        double length = mode == GcAmbiguityMode.Remove ? strongWeakCount : totalLength;
        return length > 0 ? gc / length : 0;
    }

    /// <summary>
    /// String overload of <see cref="CalculateGcFraction(ReadOnlySpan{char},GcAmbiguityMode)"/>.
    /// Opt-in Biopython <c>gc_fraction</c> compatibility; null/empty returns 0.
    /// </summary>
    public static double CalculateGcFraction(this string sequence, GcAmbiguityMode mode)
    {
        return string.IsNullOrEmpty(sequence) ? 0 : sequence.AsSpan().CalculateGcFraction(mode);
    }

    private static double WeightedGcValue(char upper) => upper switch
    {
        'V' => WeightedV,
        'B' => WeightedB,
        'H' => WeightedH,
        'D' => WeightedD,
        'M' or 'R' or 'Y' or 'K' or 'X' or 'N' => WeightedHalf,
        _ => 0.0 // any character not in the IUPAC weighted set contributes nothing
    };

    #endregion

    #region DNA Complement Core

    /// <summary>
    /// Gets the complement of a single DNA nucleotide (IUPAC-complete).
    /// Standard: A ↔ T, C ↔ G, U → A.
    /// Ambiguity codes per IUPAC NC-IUB 1984: R ↔ Y, S ↔ S, W ↔ W, K ↔ M, B ↔ V, D ↔ H, N ↔ N.
    /// Case-insensitive input, always returns uppercase.
    /// Non-IUPAC characters (gaps, etc.) pass through unchanged.
    /// </summary>
    /// <remarks>
    /// Source: Wikipedia Nucleic acid notation — IUPAC complement table.
    /// Cross-verified with Biopython complement() examples.
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static char GetComplementBase(char nucleotide) => nucleotide switch
    {
        'A' or 'a' => 'T',
        'T' or 't' => 'A',
        'G' or 'g' => 'C',
        'C' or 'c' => 'G',
        'U' or 'u' => 'A', // RNA: U pairs with A
        'R' or 'r' => 'Y', // Purine (A|G) → Pyrimidine (C|T)
        'Y' or 'y' => 'R', // Pyrimidine (C|T) → Purine (A|G)
        'S' or 's' => 'S', // Strong (C|G) → Strong (C|G)
        'W' or 'w' => 'W', // Weak (A|T) → Weak (A|T)
        'K' or 'k' => 'M', // Keto (G|T) → Amino (A|C)
        'M' or 'm' => 'K', // Amino (A|C) → Keto (G|T)
        'B' or 'b' => 'V', // Not A (C|G|T) → Not T (A|C|G)
        'D' or 'd' => 'H', // Not C (A|G|T) → Not G (A|C|T)
        'H' or 'h' => 'D', // Not G (A|C|T) → Not C (A|G|T)
        'V' or 'v' => 'B', // Not T (A|C|G) → Not A (C|G|T)
        'N' or 'n' => 'N', // Any nucleotide → Any nucleotide
        _ => nucleotide     // Non-IUPAC (gaps, etc.) pass through
    };

    /// <summary>
    /// Gets the complement of a single RNA nucleotide (IUPAC-complete).
    /// Standard: A ↔ U, C ↔ G, T → A (T accepted, pairs with A).
    /// Ambiguity codes per IUPAC NC-IUB 1984 (emitting RNA alphabet, U not T):
    /// R ↔ Y, S ↔ S, W ↔ W, K ↔ M, B ↔ V, D ↔ H, N ↔ N.
    /// Case-insensitive input, always returns uppercase for recognized bases.
    /// Non-IUPAC characters (gaps, etc.) pass through unchanged.
    /// </summary>
    /// <remarks>
    /// Matches Biopython complement_rna() behavior. The ambiguity complements are
    /// identical to the DNA path; only the base alphabet differs (T → U).
    /// Cross-verified against Biopython: complement_rna("ACGTUacgtuXYZxyz") → "UGCAAugcaaXRZxrz"
    /// (recognized bases uppercased here; unknowns pass through verbatim).
    /// See: https://biopython.org/docs/latest/api/Bio.Seq.html
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static char GetRnaComplementBase(char nucleotide) => nucleotide switch
    {
        'A' or 'a' => 'U',
        'U' or 'u' => 'A',
        'G' or 'g' => 'C',
        'C' or 'c' => 'G',
        'T' or 't' => 'A', // T pairs with A (Biopython: complement_rna("T") → "A")
        'R' or 'r' => 'Y', // Purine (A|G) → Pyrimidine (C|U)
        'Y' or 'y' => 'R', // Pyrimidine (C|U) → Purine (A|G)
        'S' or 's' => 'S', // Strong (G|C) → Strong (G|C)
        'W' or 'w' => 'W', // Weak (A|U) → Weak (A|U)
        'K' or 'k' => 'M', // Keto (G|U) → Amino (A|C)
        'M' or 'm' => 'K', // Amino (A|C) → Keto (G|U)
        'B' or 'b' => 'V', // Not A (C|G|U) → Not U (A|C|G)
        'D' or 'd' => 'H', // Not C (A|G|U) → Not G (A|C|U)
        'H' or 'h' => 'D', // Not G (A|C|U) → Not C (A|G|U)
        'V' or 'v' => 'B', // Not U (A|C|G) → Not A (C|G|U)
        'N' or 'n' => 'N', // Any nucleotide → Any nucleotide
        _ => nucleotide    // Non-IUPAC (gaps, etc.) pass through unchanged
    };

    #endregion

    #region Span-based Complement

    /// <summary>
    /// Gets the DNA complement into a destination span.
    /// </summary>
    /// <returns>True if successful, false if destination is too small.</returns>
    public static bool TryGetComplement(this ReadOnlySpan<char> sequence, Span<char> destination)
    {
        if (destination.Length < sequence.Length)
            return false;

        for (int i = 0; i < sequence.Length; i++)
        {
            destination[i] = GetComplementBase(sequence[i]);
        }

        return true;
    }

    /// <summary>
    /// Gets the DNA reverse complement into a destination span.
    /// </summary>
    public static bool TryGetReverseComplement(this ReadOnlySpan<char> sequence, Span<char> destination)
    {
        if (destination.Length < sequence.Length)
            return false;

        for (int i = 0; i < sequence.Length; i++)
        {
            destination[i] = GetComplementBase(sequence[sequence.Length - 1 - i]);
        }

        return true;
    }

    #endregion

    #region Span-based K-mer Operations

    /// <summary>
    /// Counts k-mers using span-based iteration (memory efficient).
    /// </summary>
    public static Dictionary<string, int> CountKmersSpan(this ReadOnlySpan<char> sequence, int k)
    {
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        var counts = new Dictionary<string, int>();

        if (sequence.Length < k)
            return counts;

        for (int i = 0; i <= sequence.Length - k; i++)
        {
            var kmer = sequence.Slice(i, k);
            var kmerStr = new string(kmer).ToUpperInvariant();

            if (!counts.TryAdd(kmerStr, 1))
                counts[kmerStr]++;
        }

        return counts;
    }

    /// <summary>
    /// Enumerates k-mers without allocating strings (yields spans).
    /// Use with caution - spans are only valid during enumeration.
    /// </summary>
    public static KmerEnumerator EnumerateKmers(this ReadOnlySpan<char> sequence, int k)
    {
        return new KmerEnumerator(sequence, k);
    }

    #endregion

    #region Span-based Hamming Distance

    /// <summary>
    /// Calculates Hamming distance between two spans of equal length.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int HammingDistance(this ReadOnlySpan<char> s1, ReadOnlySpan<char> s2)
    {
        if (s1.Length != s2.Length)
            throw new ArgumentException("Spans must have equal length for Hamming distance.");

        int distance = 0;
        for (int i = 0; i < s1.Length; i++)
        {
            if (char.ToUpperInvariant(s1[i]) != char.ToUpperInvariant(s2[i]))
                distance++;
        }

        return distance;
    }

    #endregion

    #region Validation

    /// <summary>
    /// Validates that a span contains only unambiguous DNA nucleotides A, C, G, T (case-insensitive).
    /// </summary>
    /// <remarks>
    /// Alphabet = IUPAC-IUB (1970) DNA bases = Biopython <c>IUPACData.unambiguous_dna_letters</c> ("GATC").
    /// IUPAC ambiguity codes (N, R, Y, …), U, gaps, whitespace and any other character are rejected;
    /// use <see cref="IsValidIupacDna"/> to accept incompletely specified bases. Only ASCII
    /// lower-case letters fold to upper case. An empty span is valid (vacuous truth; Biopython:
    /// "Zero-length sequences are always considered to be defined").
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool IsValidDna(this ReadOnlySpan<char> sequence)
        => IndexOfInvalidDna(sequence) < 0;

    /// <summary>
    /// Validates that a span contains only unambiguous RNA nucleotides A, C, G, U (case-insensitive).
    /// </summary>
    /// <remarks>
    /// Alphabet = IUPAC-IUB (1970) RNA bases = Biopython <c>IUPACData.unambiguous_rna_letters</c> ("GAUC").
    /// Same conventions as <see cref="IsValidDna"/> (T is rejected; empty span is valid).
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool IsValidRna(this ReadOnlySpan<char> sequence)
        => IndexOfInvalidRna(sequence) < 0;

    /// <summary>
    /// Returns the index of the first character that is not an unambiguous DNA nucleotide
    /// (A, C, G, T; ASCII case-insensitive), or -1 when every character is valid.
    /// Canonical predicate behind <see cref="IsValidDna"/> and <see cref="DnaSequence"/> validation.
    /// </summary>
    public static int IndexOfInvalidDna(this ReadOnlySpan<char> sequence)
    {
        for (int i = 0; i < sequence.Length; i++)
        {
            if (sequence[i] is not ('A' or 'C' or 'G' or 'T' or 'a' or 'c' or 'g' or 't'))
                return i;
        }
        return -1;
    }

    /// <summary>
    /// Returns the index of the first character that is not an unambiguous RNA nucleotide
    /// (A, C, G, U; ASCII case-insensitive), or -1 when every character is valid.
    /// Canonical predicate behind <see cref="IsValidRna"/> and <see cref="RnaSequence"/> validation.
    /// </summary>
    public static int IndexOfInvalidRna(this ReadOnlySpan<char> sequence)
    {
        for (int i = 0; i < sequence.Length; i++)
        {
            if (sequence[i] is not ('A' or 'C' or 'G' or 'U' or 'a' or 'c' or 'g' or 'u'))
                return i;
        }
        return -1;
    }

    /// <summary>
    /// Validates that a span contains only IUPAC DNA nucleotide codes — the 4 bases A, C, G, T plus the
    /// 11 NC-IUB (1984) incompletely-specified-base codes R, Y, S, W, K, M, B, D, H, V, N
    /// (ASCII case-insensitive).
    /// </summary>
    /// <remarks>
    /// Alphabet = Biopython <c>IUPACData.ambiguous_dna_letters</c> ("GATCRYWSMKHBVDN") = scikit-bio
    /// <c>DNA.definite_chars ∪ DNA.degenerate_chars</c>. U, gap symbols ('-', '.'), X and any other
    /// character are rejected (NC-IUB 1984 defines no gap symbol). Empty span is valid.
    /// </remarks>
    public static bool IsValidIupacDna(this ReadOnlySpan<char> sequence)
    {
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = ToUpperAscii(sequence[i]);
            if (c == 'U' || !IupacHelper.IsNucleotideCode(c))
                return false;
        }
        return true;
    }

    /// <summary>
    /// Validates that a span contains only IUPAC RNA nucleotide codes — A, C, G, U plus
    /// R, Y, S, W, K, M, B, D, H, V, N (ASCII case-insensitive).
    /// </summary>
    /// <remarks>
    /// Alphabet = Biopython <c>IUPACData.ambiguous_rna_letters</c> ("GAUCRYWSMKHBVDN") = scikit-bio
    /// <c>RNA.definite_chars ∪ RNA.degenerate_chars</c>. T, gaps and any other character are rejected.
    /// Empty span is valid.
    /// </remarks>
    public static bool IsValidIupacRna(this ReadOnlySpan<char> sequence)
    {
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = ToUpperAscii(sequence[i]);
            if (c == 'T' || !(c == 'U' || IupacHelper.IsNucleotideCode(c)))
                return false;
        }
        return true;
    }

    /// <summary>
    /// Upper-cases ASCII letters only. Unlike <see cref="char.ToUpperInvariant(char)"/> this never maps a
    /// non-ASCII character onto an ASCII sequence symbol (e.g. U+017F 'ſ' → 'S').
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static char ToUpperAscii(char c) => char.IsAsciiLetterLower(c) ? (char)(c - 32) : c;

    #endregion
}

/// <summary>
/// Enumerator for iterating k-mers as spans without string allocation.
/// </summary>
public ref struct KmerEnumerator
{
    private readonly ReadOnlySpan<char> _sequence;
    private readonly int _k;
    private int _position;

    internal KmerEnumerator(ReadOnlySpan<char> sequence, int k)
    {
        _sequence = sequence;
        _k = k;
        _position = -1;
    }

    public ReadOnlySpan<char> Current => _sequence.Slice(_position, _k);

    public bool MoveNext()
    {
        _position++;
        return _position <= _sequence.Length - _k;
    }

    public KmerEnumerator GetEnumerator() => this;
}
