using System;
using System.Buffers;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Threading;

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
    /// Matches Wikipedia GC-content formula and Biopython gc_fraction ("remove" mode).
    /// Valid nucleotides: A, T, G, C, U (case-insensitive).
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double CalculateGcContent(this ReadOnlySpan<char> sequence)
    {
        if (sequence.IsEmpty) return 0;

        int gcCount = 0;
        int validCount = 0;
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = sequence[i];
            switch (c)
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

        return validCount == 0 ? 0 : (double)gcCount / validCount * 100;
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
    /// Matches Wikipedia GC-content formula and Biopython gc_fraction ("remove" mode).
    /// Valid nucleotides: A, T, G, C, U (case-insensitive).
    /// </remarks>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double CalculateGcFraction(this ReadOnlySpan<char> sequence)
    {
        if (sequence.IsEmpty) return 0;

        int gcCount = 0;
        int validCount = 0;
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = sequence[i];
            switch (c)
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

    #region DNA Complement Core

    /// <summary>
    /// Gets the complement of a single DNA nucleotide.
    /// A ↔ T, C ↔ G (case-insensitive, returns uppercase).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static char GetComplementBase(char nucleotide) => nucleotide switch
    {
        'A' or 'a' => 'T',
        'T' or 't' => 'A',
        'G' or 'g' => 'C',
        'C' or 'c' => 'G',
        'U' or 'u' => 'A', // RNA support
        _ => nucleotide
    };

    /// <summary>
    /// Gets the complement of a single RNA nucleotide.
    /// A ↔ U, C ↔ G, T → A (case-insensitive, returns uppercase).
    /// Unknown bases are returned unchanged.
    /// </summary>
    /// <remarks>
    /// Matches Biopython complement_rna() behavior:
    /// - T is accepted and complements to A (T pairs with A regardless of context).
    /// - Unknown/non-nucleotide characters pass through unchanged.
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
        _ => nucleotide    // Unknown bases returned unchanged (Biopython behavior)
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
        var counts = new Dictionary<string, int>();

        if (sequence.Length < k || k <= 0)
            return counts;

        for (int i = 0; i <= sequence.Length - k; i++)
        {
            var kmer = sequence.Slice(i, k);
            var kmerStr = new string(kmer);

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
    /// Validates if a span contains only valid DNA characters.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool IsValidDna(this ReadOnlySpan<char> sequence)
    {
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = char.ToUpperInvariant(sequence[i]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                return false;
        }
        return true;
    }

    /// <summary>
    /// Validates if a span contains only valid RNA characters.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool IsValidRna(this ReadOnlySpan<char> sequence)
    {
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = char.ToUpperInvariant(sequence[i]);
            if (c != 'A' && c != 'C' && c != 'G' && c != 'U')
                return false;
        }
        return true;
    }

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
