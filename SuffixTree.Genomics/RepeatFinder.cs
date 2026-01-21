using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics;

/// <summary>
/// Finds various types of repeats in DNA sequences including microsatellites (STRs),
/// minisatellites (VNTRs), inverted repeats, and direct repeats.
/// </summary>
public static class RepeatFinder
{
    #region Microsatellite (STR) Detection

    /// <summary>
    /// Finds microsatellites (Short Tandem Repeats) in a DNA sequence.
    /// STRs are 1-6 bp motifs repeated consecutively.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minUnitLength">Minimum repeat unit length (default: 1).</param>
    /// <param name="maxUnitLength">Maximum repeat unit length (default: 6).</param>
    /// <param name="minRepeats">Minimum number of repeats to report (default: 3).</param>
    /// <returns>Collection of microsatellite repeats found.</returns>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        DnaSequence sequence,
        int minUnitLength = 1,
        int maxUnitLength = 6,
        int minRepeats = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minUnitLength < 1) throw new ArgumentOutOfRangeException(nameof(minUnitLength));
        if (maxUnitLength < minUnitLength) throw new ArgumentOutOfRangeException(nameof(maxUnitLength));
        if (minRepeats < 2) throw new ArgumentOutOfRangeException(nameof(minRepeats));

        return FindMicrosatellitesCore(sequence.Sequence, minUnitLength, maxUnitLength, minRepeats);
    }

    /// <summary>
    /// Finds microsatellites in a raw sequence string.
    /// </summary>
    public static IEnumerable<MicrosatelliteResult> FindMicrosatellites(
        string sequence,
        int minUnitLength = 1,
        int maxUnitLength = 6,
        int minRepeats = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var result in FindMicrosatellitesCore(sequence.ToUpperInvariant(), minUnitLength, maxUnitLength, minRepeats))
            yield return result;
    }

    private static IEnumerable<MicrosatelliteResult> FindMicrosatellitesCore(
        string seq,
        int minUnitLength,
        int maxUnitLength,
        int minRepeats)
    {
        var reported = new HashSet<(int Start, int End)>();

        for (int unitLen = minUnitLength; unitLen <= maxUnitLength; unitLen++)
        {
            for (int i = 0; i <= seq.Length - unitLen * minRepeats; i++)
            {
                string unit = seq.Substring(i, unitLen);

                // Skip if unit is just repetition of smaller unit
                if (IsRedundantUnit(unit))
                    continue;

                int repeats = 1;
                int j = i + unitLen;

                while (j + unitLen <= seq.Length && seq.Substring(j, unitLen) == unit)
                {
                    repeats++;
                    j += unitLen;
                }

                if (repeats >= minRepeats)
                {
                    int end = i + (repeats * unitLen) - 1;
                    
                    // Avoid reporting overlapping/contained repeats
                    if (!reported.Any(r => r.Start <= i && r.End >= end))
                    {
                        reported.Add((i, end));
                        yield return new MicrosatelliteResult(
                            Position: i,
                            RepeatUnit: unit,
                            RepeatCount: repeats,
                            TotalLength: repeats * unitLen,
                            RepeatType: ClassifyRepeatType(unit));
                    }
                }
            }
        }
    }

    private static bool IsRedundantUnit(string unit)
    {
        if (unit.Length <= 1) return false;

        // Check if unit is made of smaller repeating pattern
        for (int subLen = 1; subLen < unit.Length; subLen++)
        {
            if (unit.Length % subLen != 0) continue;

            string subUnit = unit.Substring(0, subLen);
            bool isRedundant = true;

            for (int i = subLen; i < unit.Length; i += subLen)
            {
                if (unit.Substring(i, subLen) != subUnit)
                {
                    isRedundant = false;
                    break;
                }
            }

            if (isRedundant) return true;
        }

        return false;
    }

    private static RepeatType ClassifyRepeatType(string unit)
    {
        return unit.Length switch
        {
            1 => RepeatType.Mononucleotide,
            2 => RepeatType.Dinucleotide,
            3 => RepeatType.Trinucleotide,
            4 => RepeatType.Tetranucleotide,
            5 => RepeatType.Pentanucleotide,
            6 => RepeatType.Hexanucleotide,
            _ => RepeatType.Complex
        };
    }

    #endregion

    #region Inverted Repeat Detection

    /// <summary>
    /// Finds inverted repeats (sequences that are reverse complements of each other).
    /// These can form hairpin/stem-loop structures.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minArmLength">Minimum length of each arm (default: 4).</param>
    /// <param name="maxLoopLength">Maximum loop length between arms (default: 50).</param>
    /// <param name="minLoopLength">Minimum loop length (default: 3).</param>
    /// <returns>Collection of inverted repeats found.</returns>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        DnaSequence sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minArmLength < 2) throw new ArgumentOutOfRangeException(nameof(minArmLength));
        if (minLoopLength < 0) throw new ArgumentOutOfRangeException(nameof(minLoopLength));

        return FindInvertedRepeatsCore(sequence.Sequence, minArmLength, maxLoopLength, minLoopLength);
    }

    /// <summary>
    /// Finds inverted repeats in a raw sequence string.
    /// </summary>
    public static IEnumerable<InvertedRepeatResult> FindInvertedRepeats(
        string sequence,
        int minArmLength = 4,
        int maxLoopLength = 50,
        int minLoopLength = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var result in FindInvertedRepeatsCore(sequence.ToUpperInvariant(), minArmLength, maxLoopLength, minLoopLength))
            yield return result;
    }

    private static IEnumerable<InvertedRepeatResult> FindInvertedRepeatsCore(
        string seq,
        int minArmLength,
        int maxLoopLength,
        int minLoopLength)
    {
        var reported = new HashSet<(int, int, int)>();

        for (int i = 0; i <= seq.Length - 2 * minArmLength - minLoopLength; i++)
        {
            for (int armLen = minArmLength; i + armLen <= seq.Length; armLen++)
            {
                string leftArm = seq.Substring(i, armLen);
                string leftArmRevComp = GetReverseComplement(leftArm);

                // Search for right arm
                int minJ = i + armLen + minLoopLength;
                int maxJ = Math.Min(i + armLen + maxLoopLength, seq.Length - armLen);

                for (int j = minJ; j <= maxJ; j++)
                {
                    if (j + armLen > seq.Length) break;

                    string rightArm = seq.Substring(j, armLen);

                    if (rightArm == leftArmRevComp)
                    {
                        int loopLength = j - (i + armLen);
                        string loop = loopLength > 0 ? seq.Substring(i + armLen, loopLength) : "";

                        var key = (i, j, armLen);
                        if (!reported.Contains(key))
                        {
                            reported.Add(key);
                            yield return new InvertedRepeatResult(
                                LeftArmStart: i,
                                RightArmStart: j,
                                ArmLength: armLen,
                                LoopLength: loopLength,
                                LeftArm: leftArm,
                                RightArm: rightArm,
                                Loop: loop,
                                CanFormHairpin: loopLength >= 3);
                        }
                    }
                }
            }
        }
    }

    private static string GetReverseComplement(string sequence)
    {
        var complement = new char[sequence.Length];
        for (int i = 0; i < sequence.Length; i++)
        {
            complement[sequence.Length - 1 - i] = sequence[i] switch
            {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => sequence[i]
            };
        }
        return new string(complement);
    }

    #endregion

    #region Direct Repeat Detection

    /// <summary>
    /// Finds direct repeats (identical sequences appearing multiple times).
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum repeat length (default: 5).</param>
    /// <param name="maxLength">Maximum repeat length (default: 50).</param>
    /// <param name="minSpacing">Minimum spacing between repeats (default: 1).</param>
    /// <returns>Collection of direct repeats found.</returns>
    public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
        DnaSequence sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minLength < 2) throw new ArgumentOutOfRangeException(nameof(minLength));
        if (maxLength < minLength) throw new ArgumentOutOfRangeException(nameof(maxLength));

        return FindDirectRepeatsCore(sequence.Sequence, minLength, maxLength, minSpacing);
    }

    /// <summary>
    /// Finds direct repeats in a raw sequence string.
    /// </summary>
    public static IEnumerable<DirectRepeatResult> FindDirectRepeats(
        string sequence,
        int minLength = 5,
        int maxLength = 50,
        int minSpacing = 1)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var result in FindDirectRepeatsCore(sequence.ToUpperInvariant(), minLength, maxLength, minSpacing))
            yield return result;
    }

    private static IEnumerable<DirectRepeatResult> FindDirectRepeatsCore(
        string seq,
        int minLength,
        int maxLength,
        int minSpacing)
    {
        var reported = new HashSet<(int, int, int)>();

        for (int len = minLength; len <= maxLength; len++)
        {
            for (int i = 0; i <= seq.Length - len * 2 - minSpacing; i++)
            {
                string repeat = seq.Substring(i, len);

                // Search for second occurrence
                for (int j = i + len + minSpacing; j <= seq.Length - len; j++)
                {
                    if (seq.Substring(j, len) == repeat)
                    {
                        var key = (i, j, len);
                        if (!reported.Contains(key))
                        {
                            reported.Add(key);
                            yield return new DirectRepeatResult(
                                FirstPosition: i,
                                SecondPosition: j,
                                RepeatSequence: repeat,
                                Length: len,
                                Spacing: j - i - len);
                        }
                    }
                }
            }
        }
    }

    #endregion

    #region Tandem Repeat Summary

    /// <summary>
    /// Gets a summary of all tandem repeats in a sequence.
    /// </summary>
    public static TandemRepeatSummary GetTandemRepeatSummary(
        DnaSequence sequence,
        int minRepeats = 3)
    {
        ArgumentNullException.ThrowIfNull(sequence);

        var microsatellites = FindMicrosatellites(sequence, 1, 6, minRepeats).ToList();

        var byType = microsatellites
            .GroupBy(m => m.RepeatType)
            .ToDictionary(g => g.Key, g => g.ToList());

        int totalBases = microsatellites.Sum(m => m.TotalLength);
        double percentageOfSequence = sequence.Length > 0
            ? (double)totalBases / sequence.Length * 100
            : 0;

        return new TandemRepeatSummary(
            TotalRepeats: microsatellites.Count,
            TotalRepeatBases: totalBases,
            PercentageOfSequence: percentageOfSequence,
            MononucleotideRepeats: byType.GetValueOrDefault(RepeatType.Mononucleotide)?.Count ?? 0,
            DinucleotideRepeats: byType.GetValueOrDefault(RepeatType.Dinucleotide)?.Count ?? 0,
            TrinucleotideRepeats: byType.GetValueOrDefault(RepeatType.Trinucleotide)?.Count ?? 0,
            TetranucleotideRepeats: byType.GetValueOrDefault(RepeatType.Tetranucleotide)?.Count ?? 0,
            LongestRepeat: microsatellites.OrderByDescending(m => m.TotalLength).FirstOrDefault(),
            MostFrequentUnit: microsatellites
                .GroupBy(m => m.RepeatUnit)
                .OrderByDescending(g => g.Count())
                .FirstOrDefault()?.Key);
    }

    #endregion

    #region Palindrome Detection

    /// <summary>
    /// Finds palindromic sequences (sequences that read the same 5' to 3' on both strands).
    /// These are recognition sites for many restriction enzymes.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="minLength">Minimum palindrome length (default: 4, must be even).</param>
    /// <param name="maxLength">Maximum palindrome length (default: 12).</param>
    /// <returns>Collection of palindromes found.</returns>
    public static IEnumerable<PalindromeResult> FindPalindromes(
        DnaSequence sequence,
        int minLength = 4,
        int maxLength = 12)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (minLength < 4 || minLength % 2 != 0) 
            throw new ArgumentOutOfRangeException(nameof(minLength), "Must be even and >= 4");
        if (maxLength < minLength) 
            throw new ArgumentOutOfRangeException(nameof(maxLength));

        return FindPalindromesCore(sequence.Sequence, minLength, maxLength);
    }

    /// <summary>
    /// Finds palindromes in a raw sequence string.
    /// </summary>
    public static IEnumerable<PalindromeResult> FindPalindromes(
        string sequence,
        int minLength = 4,
        int maxLength = 12)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var result in FindPalindromesCore(sequence.ToUpperInvariant(), minLength, maxLength))
            yield return result;
    }

    private static IEnumerable<PalindromeResult> FindPalindromesCore(
        string seq,
        int minLength,
        int maxLength)
    {
        for (int len = minLength; len <= maxLength; len += 2) // Palindromes must be even length
        {
            for (int i = 0; i <= seq.Length - len; i++)
            {
                string candidate = seq.Substring(i, len);
                string revComp = GetReverseComplement(candidate);

                if (candidate == revComp)
                {
                    yield return new PalindromeResult(
                        Position: i,
                        Sequence: candidate,
                        Length: len);
                }
            }
        }
    }

    #endregion
}

/// <summary>
/// Type of repeat unit.
/// </summary>
public enum RepeatType
{
    Mononucleotide,  // A, T, G, C
    Dinucleotide,    // AT, GC, CA, etc.
    Trinucleotide,   // CAG, CGG, etc.
    Tetranucleotide, // GATA, AAAT, etc.
    Pentanucleotide, // AAAAT, etc.
    Hexanucleotide,  // AAAAAG, etc.
    Complex          // Longer than 6bp
}

/// <summary>
/// Result of microsatellite (STR) detection.
/// </summary>
public readonly record struct MicrosatelliteResult(
    int Position,
    string RepeatUnit,
    int RepeatCount,
    int TotalLength,
    RepeatType RepeatType)
{
    /// <summary>
    /// Gets the full repeat sequence.
    /// </summary>
    public string FullSequence => string.Concat(Enumerable.Repeat(RepeatUnit, RepeatCount));
}

/// <summary>
/// Result of inverted repeat detection.
/// </summary>
public readonly record struct InvertedRepeatResult(
    int LeftArmStart,
    int RightArmStart,
    int ArmLength,
    int LoopLength,
    string LeftArm,
    string RightArm,
    string Loop,
    bool CanFormHairpin)
{
    /// <summary>
    /// Total length of the inverted repeat structure.
    /// </summary>
    public int TotalLength => 2 * ArmLength + LoopLength;
}

/// <summary>
/// Result of direct repeat detection.
/// </summary>
public readonly record struct DirectRepeatResult(
    int FirstPosition,
    int SecondPosition,
    string RepeatSequence,
    int Length,
    int Spacing);

/// <summary>
/// Result of palindrome detection.
/// </summary>
public readonly record struct PalindromeResult(
    int Position,
    string Sequence,
    int Length);

/// <summary>
/// Summary of tandem repeats in a sequence.
/// </summary>
public readonly record struct TandemRepeatSummary(
    int TotalRepeats,
    int TotalRepeatBases,
    double PercentageOfSequence,
    int MononucleotideRepeats,
    int DinucleotideRepeats,
    int TrinucleotideRepeats,
    int TetranucleotideRepeats,
    MicrosatelliteResult? LongestRepeat,
    string? MostFrequentUnit);
