using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuffixTree.Genomics;

/// <summary>
/// Provides algorithms for RNA secondary structure prediction and analysis.
/// Includes stem-loop detection, hairpin finding, and free energy calculations.
/// </summary>
public static class RnaSecondaryStructure
{
    #region Records and Types

    /// <summary>
    /// Represents a base pair in RNA secondary structure.
    /// </summary>
    public readonly record struct BasePair(int Position1, int Position2, char Base1, char Base2, BasePairType Type);

    /// <summary>
    /// Types of base pairs.
    /// </summary>
    public enum BasePairType
    {
        WatsonCrick,  // A-U, G-C
        Wobble,       // G-U
        NonCanonical  // Other pairs
    }

    /// <summary>
    /// Represents a stem structure (double-stranded region).
    /// </summary>
    public readonly record struct Stem(
        int Start5Prime,
        int End5Prime,
        int Start3Prime,
        int End3Prime,
        int Length,
        IReadOnlyList<BasePair> BasePairs,
        double FreeEnergy);

    /// <summary>
    /// Represents a loop structure.
    /// </summary>
    public readonly record struct Loop(
        LoopType Type,
        int Start,
        int End,
        int Size,
        string Sequence);

    /// <summary>
    /// Types of loops in RNA structure.
    /// </summary>
    public enum LoopType
    {
        Hairpin,      // Terminal loop
        Internal,     // Internal loop
        Bulge,        // Bulge loop (asymmetric)
        MultiLoop,    // Multi-branch loop
        External      // External (unpaired ends)
    }

    /// <summary>
    /// Represents a stem-loop (hairpin) structure.
    /// </summary>
    public readonly record struct StemLoop(
        int Start,
        int End,
        Stem Stem,
        Loop Loop,
        double TotalFreeEnergy,
        string DotBracketNotation);

    /// <summary>
    /// Represents a pseudoknot structure.
    /// </summary>
    public readonly record struct Pseudoknot(
        int Start1, int End1,
        int Start2, int End2,
        IReadOnlyList<BasePair> CrossingPairs);

    /// <summary>
    /// Complete secondary structure prediction result.
    /// </summary>
    public readonly record struct SecondaryStructure(
        string Sequence,
        string DotBracket,
        IReadOnlyList<BasePair> BasePairs,
        IReadOnlyList<StemLoop> StemLoops,
        IReadOnlyList<Pseudoknot> Pseudoknots,
        double MinimumFreeEnergy);

    #endregion

    #region Free Energy Parameters

    // Turner 2004 nearest-neighbor parameters (kcal/mol at 37°C)
    // Stacking energies for Watson-Crick pairs
    private static readonly Dictionary<string, double> StackingEnergies = new()
    {
        // Format: 5'-XY-3' / 3'-X'Y'-5' where X pairs with X', Y pairs with Y'
        { "AA/UU", -0.9 }, { "AU/UA", -1.1 }, { "UA/AU", -1.3 }, { "CU/GA", -2.1 },
        { "CA/GU", -2.1 }, { "GU/CA", -2.2 }, { "CG/GC", -2.4 }, { "UG/AC", -1.4 },
        { "AG/UC", -2.1 }, { "GA/CU", -2.4 }, { "GG/CC", -3.3 }, { "GC/CG", -3.4 },
        { "AC/UG", -2.2 }, { "UC/AG", -2.4 }, { "CC/GG", -3.3 }, { "UU/AA", -0.9 },
        // G-U wobble pairs
        { "GU/UG", -1.3 }, { "UG/GU", -1.3 },
    };

    // Hairpin loop initiation energies
    private static readonly Dictionary<int, double> HairpinLoopEnergies = new()
    {
        { 3, 5.4 }, { 4, 5.6 }, { 5, 5.7 }, { 6, 5.4 }, { 7, 6.0 },
        { 8, 5.5 }, { 9, 6.4 }, { 10, 6.5 }, { 12, 6.7 }, { 14, 7.0 },
        { 16, 7.2 }, { 18, 7.4 }, { 20, 7.5 }, { 25, 7.8 }, { 30, 8.0 }
    };

    // Special hairpin loop bonuses
    private static readonly Dictionary<string, double> SpecialHairpinLoops = new()
    {
        { "GAAA", -3.0 },  // GNRA tetraloop
        { "GCAA", -3.0 },
        { "GGAA", -3.0 },
        { "GUAA", -3.0 },
        { "UUCG", -3.0 },  // UNCG tetraloop
        { "UACG", -3.0 },
        { "UGCG", -3.0 },
        { "UCCG", -3.0 },
        { "CUUG", -2.0 },  // CUYG tetraloop
        { "CCUG", -2.0 },
    };

    #endregion

    #region Base Pairing

    /// <summary>
    /// Determines if two bases can form a pair.
    /// </summary>
    public static bool CanPair(char base1, char base2)
    {
        return GetBasePairType(base1, base2) != null;
    }

    /// <summary>
    /// Gets the type of base pair, or null if bases cannot pair.
    /// </summary>
    public static BasePairType? GetBasePairType(char base1, char base2)
    {
        char b1 = char.ToUpperInvariant(base1);
        char b2 = char.ToUpperInvariant(base2);

        // Watson-Crick pairs
        if ((b1 == 'A' && b2 == 'U') || (b1 == 'U' && b2 == 'A') ||
            (b1 == 'G' && b2 == 'C') || (b1 == 'C' && b2 == 'G'))
        {
            return BasePairType.WatsonCrick;
        }

        // Wobble pairs
        if ((b1 == 'G' && b2 == 'U') || (b1 == 'U' && b2 == 'G'))
        {
            return BasePairType.Wobble;
        }

        return null;
    }

    /// <summary>
    /// Gets the complement of a base (RNA).
    /// </summary>
    public static char GetComplement(char base_)
    {
        return char.ToUpperInvariant(base_) switch
        {
            'A' => 'U',
            'U' => 'A',
            'G' => 'C',
            'C' => 'G',
            _ => 'N'
        };
    }

    #endregion

    #region Stem-Loop Finding

    /// <summary>
    /// Finds all potential stem-loop structures in an RNA sequence.
    /// </summary>
    public static IEnumerable<StemLoop> FindStemLoops(
        string rnaSequence,
        int minStemLength = 3,
        int minLoopSize = 3,
        int maxLoopSize = 10,
        bool allowWobble = true)
    {
        if (string.IsNullOrEmpty(rnaSequence) || rnaSequence.Length < minStemLength * 2 + minLoopSize)
            yield break;

        string upper = rnaSequence.ToUpperInvariant();

        // Scan for potential hairpin loops
        for (int loopStart = minStemLength; loopStart <= upper.Length - minStemLength - minLoopSize; loopStart++)
        {
            for (int loopSize = minLoopSize; loopSize <= Math.Min(maxLoopSize, upper.Length - loopStart - minStemLength); loopSize++)
            {
                int loopEnd = loopStart + loopSize - 1;

                // Try to extend stem on both sides
                var stemLoop = TryBuildStemLoop(upper, loopStart, loopEnd, minStemLength, allowWobble);
                if (stemLoop != null)
                {
                    yield return stemLoop.Value;
                }
            }
        }
    }

    private static StemLoop? TryBuildStemLoop(string sequence, int loopStart, int loopEnd, int minStemLength, bool allowWobble)
    {
        var basePairs = new List<BasePair>();
        int stemLength = 0;

        int left = loopStart - 1;
        int right = loopEnd + 1;

        // Extend stem
        while (left >= 0 && right < sequence.Length)
        {
            var pairType = GetBasePairType(sequence[left], sequence[right]);

            if (pairType == null)
                break;

            if (pairType == BasePairType.Wobble && !allowWobble)
                break;

            basePairs.Add(new BasePair(left, right, sequence[left], sequence[right], pairType.Value));
            stemLength++;
            left--;
            right++;
        }

        if (stemLength < minStemLength)
            return null;

        // Build result
        int stemStart = left + 1;
        int stemEnd5 = loopStart - 1;
        int stemStart3 = loopEnd + 1;
        int stemEnd = right - 1;

        basePairs.Reverse(); // Order from 5' to 3'

        var stem = new Stem(
            Start5Prime: stemStart,
            End5Prime: stemEnd5,
            Start3Prime: stemStart3,
            End3Prime: stemEnd,
            Length: stemLength,
            BasePairs: basePairs,
            FreeEnergy: CalculateStemEnergy(sequence, basePairs));

        string loopSeq = sequence.Substring(loopStart, loopEnd - loopStart + 1);
        var loop = new Loop(
            Type: LoopType.Hairpin,
            Start: loopStart,
            End: loopEnd,
            Size: loopSeq.Length,
            Sequence: loopSeq);

        double loopEnergy = CalculateHairpinLoopEnergy(loopSeq, sequence[stemEnd5], sequence[stemStart3]);
        double totalEnergy = stem.FreeEnergy + loopEnergy;

        string dotBracket = GenerateDotBracket(sequence.Length, basePairs, stemStart, stemEnd);

        return new StemLoop(
            Start: stemStart,
            End: stemEnd,
            Stem: stem,
            Loop: loop,
            TotalFreeEnergy: totalEnergy,
            DotBracketNotation: dotBracket);
    }

    #endregion

    #region Energy Calculations

    /// <summary>
    /// Calculates the free energy of a stem region.
    /// </summary>
    public static double CalculateStemEnergy(string sequence, IReadOnlyList<BasePair> basePairs)
    {
        if (basePairs.Count < 2)
            return 0;

        double energy = 0;

        // Sum stacking energies
        for (int i = 0; i < basePairs.Count - 1; i++)
        {
            var pair1 = basePairs[i];
            var pair2 = basePairs[i + 1];

            string stack = $"{pair1.Base1}{pair2.Base1}/{pair1.Base2}{pair2.Base2}";
            energy += StackingEnergies.GetValueOrDefault(stack, -1.5); // Default stacking
        }

        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Calculates the free energy of a hairpin loop.
    /// </summary>
    public static double CalculateHairpinLoopEnergy(string loopSequence, char closingBase5, char closingBase3)
    {
        int size = loopSequence.Length;

        // Base loop initiation energy
        double energy = 0;
        if (HairpinLoopEnergies.TryGetValue(size, out double loopEnergy))
        {
            energy = loopEnergy;
        }
        else if (size > 30)
        {
            // Jacobson-Stockmayer extrapolation
            energy = 8.0 + 1.75 * 1.987 * 310.15 / 1000.0 * Math.Log((double)size / 30.0);
        }
        else
        {
            // Interpolate
            var keys = HairpinLoopEnergies.Keys.OrderBy(k => k).ToList();
            int lower = keys.LastOrDefault(k => k <= size);
            int upper = keys.FirstOrDefault(k => k >= size);
            if (lower > 0 && upper > 0)
            {
                double lowerE = HairpinLoopEnergies[lower];
                double upperE = HairpinLoopEnergies[upper];
                energy = lowerE + (upperE - lowerE) * (size - lower) / (upper - lower);
            }
            else
            {
                energy = 5.5; // Default
            }
        }

        // Tetraloop bonus
        if (size == 4 && SpecialHairpinLoops.TryGetValue(loopSequence.ToUpperInvariant(), out double bonus))
        {
            energy += bonus;
        }

        // Terminal mismatch penalty (simplified)
        char firstLoop = loopSequence.Length > 0 ? loopSequence[0] : 'N';
        char lastLoop = loopSequence.Length > 0 ? loopSequence[^1] : 'N';

        // All-C loop penalty
        if (loopSequence.All(c => char.ToUpperInvariant(c) == 'C'))
        {
            energy += 0.3 * size;
        }

        return Math.Round(energy, 2);
    }

    /// <summary>
    /// Calculates the minimum free energy (MFE) of an RNA sequence using dynamic programming.
    /// Simplified Nussinov-like algorithm.
    /// </summary>
    public static double CalculateMinimumFreeEnergy(string rnaSequence, int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(rnaSequence) || rnaSequence.Length < minLoopSize + 2)
            return 0;

        string seq = rnaSequence.ToUpperInvariant();
        int n = seq.Length;

        // DP table: dp[i,j] = MFE for subsequence from i to j
        var dp = new double[n, n];

        // Initialize
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                dp[i, j] = 0;

        // Fill DP table
        for (int length = minLoopSize + 2; length <= n; length++)
        {
            for (int i = 0; i <= n - length; i++)
            {
                int j = i + length - 1;

                // Option 1: j is unpaired
                dp[i, j] = dp[i, j - 1];

                // Option 2: j pairs with some k
                for (int k = i; k < j - minLoopSize; k++)
                {
                    if (CanPair(seq[k], seq[j]))
                    {
                        double pairEnergy = GetPairEnergy(seq[k], seq[j]);
                        double left = k > i ? dp[i, k - 1] : 0;
                        double enclosed = k + 1 < j ? dp[k + 1, j - 1] : 0;

                        double total = pairEnergy + left + enclosed;
                        dp[i, j] = Math.Min(dp[i, j], total);
                    }
                }
            }
        }

        return Math.Round(dp[0, n - 1], 2);
    }

    private static double GetPairEnergy(char base1, char base2)
    {
        var pairType = GetBasePairType(base1, base2);
        return pairType switch
        {
            BasePairType.WatsonCrick => -2.0,
            BasePairType.Wobble => -1.0,
            _ => 0
        };
    }

    #endregion

    #region Structure Prediction

    /// <summary>
    /// Predicts the secondary structure of an RNA sequence.
    /// </summary>
    public static SecondaryStructure PredictStructure(
        string rnaSequence,
        int minStemLength = 3,
        int minLoopSize = 3,
        int maxLoopSize = 10)
    {
        if (string.IsNullOrEmpty(rnaSequence))
        {
            return new SecondaryStructure(
                "", "", new List<BasePair>(), new List<StemLoop>(),
                new List<Pseudoknot>(), 0);
        }

        string seq = rnaSequence.ToUpperInvariant();

        // Find stem-loops
        var stemLoops = FindStemLoops(seq, minStemLength, minLoopSize, maxLoopSize)
            .OrderBy(sl => sl.TotalFreeEnergy)
            .ToList();

        // Select non-overlapping structures greedily by energy
        var selectedStemLoops = SelectNonOverlapping(stemLoops);

        // Collect all base pairs
        var allBasePairs = selectedStemLoops
            .SelectMany(sl => sl.Stem.BasePairs)
            .OrderBy(bp => bp.Position1)
            .ToList();

        // Generate dot-bracket notation
        string dotBracket = GenerateFullDotBracket(seq.Length, allBasePairs);

        // Calculate total MFE
        double mfe = selectedStemLoops.Sum(sl => sl.TotalFreeEnergy);

        // Detect pseudoknots
        var pseudoknots = DetectPseudoknots(allBasePairs).ToList();

        return new SecondaryStructure(
            Sequence: seq,
            DotBracket: dotBracket,
            BasePairs: allBasePairs,
            StemLoops: selectedStemLoops,
            Pseudoknots: pseudoknots,
            MinimumFreeEnergy: mfe);
    }

    private static List<StemLoop> SelectNonOverlapping(List<StemLoop> stemLoops)
    {
        var selected = new List<StemLoop>();
        var used = new HashSet<int>();

        foreach (var sl in stemLoops)
        {
            bool overlaps = false;
            for (int pos = sl.Start; pos <= sl.End; pos++)
            {
                if (used.Contains(pos))
                {
                    overlaps = true;
                    break;
                }
            }

            if (!overlaps)
            {
                selected.Add(sl);
                for (int pos = sl.Start; pos <= sl.End; pos++)
                {
                    used.Add(pos);
                }
            }
        }

        return selected;
    }

    #endregion

    #region Pseudoknot Detection

    /// <summary>
    /// Detects pseudoknots in a set of base pairs.
    /// A pseudoknot occurs when pairs cross: i < i' < j < j' for pairs (i,j) and (i',j').
    /// </summary>
    public static IEnumerable<Pseudoknot> DetectPseudoknots(IReadOnlyList<BasePair> basePairs)
    {
        var crossingGroups = new List<List<BasePair>>();

        for (int i = 0; i < basePairs.Count; i++)
        {
            for (int j = i + 1; j < basePairs.Count; j++)
            {
                var bp1 = basePairs[i];
                var bp2 = basePairs[j];

                // Check for crossing: i < i' < j < j'
                int i1 = Math.Min(bp1.Position1, bp1.Position2);
                int j1 = Math.Max(bp1.Position1, bp1.Position2);
                int i2 = Math.Min(bp2.Position1, bp2.Position2);
                int j2 = Math.Max(bp2.Position1, bp2.Position2);

                if (i1 < i2 && i2 < j1 && j1 < j2)
                {
                    // Found crossing pairs
                    yield return new Pseudoknot(
                        Start1: i1,
                        End1: j1,
                        Start2: i2,
                        End2: j2,
                        CrossingPairs: new List<BasePair> { bp1, bp2 });
                }
            }
        }
    }

    #endregion

    #region Dot-Bracket Notation

    /// <summary>
    /// Generates dot-bracket notation for a structure.
    /// </summary>
    private static string GenerateDotBracket(int length, IReadOnlyList<BasePair> basePairs, int start, int end)
    {
        var notation = new char[end - start + 1];
        for (int i = 0; i < notation.Length; i++)
            notation[i] = '.';

        foreach (var bp in basePairs)
        {
            if (bp.Position1 >= start && bp.Position1 <= end)
                notation[bp.Position1 - start] = '(';
            if (bp.Position2 >= start && bp.Position2 <= end)
                notation[bp.Position2 - start] = ')';
        }

        return new string(notation);
    }

    private static string GenerateFullDotBracket(int length, IReadOnlyList<BasePair> basePairs)
    {
        var notation = new char[length];
        for (int i = 0; i < length; i++)
            notation[i] = '.';

        foreach (var bp in basePairs)
        {
            int left = Math.Min(bp.Position1, bp.Position2);
            int right = Math.Max(bp.Position1, bp.Position2);
            notation[left] = '(';
            notation[right] = ')';
        }

        return new string(notation);
    }

    /// <summary>
    /// Parses dot-bracket notation to extract base pairs.
    /// </summary>
    public static IEnumerable<(int Position1, int Position2)> ParseDotBracket(string dotBracket)
    {
        var stack = new Stack<int>();

        for (int i = 0; i < dotBracket.Length; i++)
        {
            switch (dotBracket[i])
            {
                case '(':
                case '[':
                case '{':
                case '<':
                    stack.Push(i);
                    break;

                case ')':
                case ']':
                case '}':
                case '>':
                    if (stack.Count > 0)
                    {
                        int j = stack.Pop();
                        yield return (j, i);
                    }
                    break;
            }
        }
    }

    /// <summary>
    /// Validates dot-bracket notation.
    /// </summary>
    public static bool ValidateDotBracket(string dotBracket)
    {
        int count = 0;
        foreach (char c in dotBracket)
        {
            if (c == '(' || c == '[' || c == '{' || c == '<')
                count++;
            else if (c == ')' || c == ']' || c == '}' || c == '>')
                count--;

            if (count < 0)
                return false;
        }
        return count == 0;
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Finds inverted repeats (potential stem regions).
    /// </summary>
    public static IEnumerable<(int Start1, int End1, int Start2, int End2, int Length)> FindInvertedRepeats(
        string sequence,
        int minLength = 4,
        int minSpacing = 3,
        int maxSpacing = 100)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minLength * 2 + minSpacing)
            yield break;

        string upper = sequence.ToUpperInvariant();

        for (int i = 0; i <= upper.Length - minLength * 2 - minSpacing; i++)
        {
            for (int spacing = minSpacing; spacing <= Math.Min(maxSpacing, upper.Length - i - minLength * 2); spacing++)
            {
                int j = i + minLength + spacing;

                // Check for complementary regions
                int matchLen = 0;
                while (i + matchLen < j && j + matchLen < upper.Length)
                {
                    char expected = GetComplement(upper[j + matchLen]);
                    if (upper[i + matchLen] != expected)
                        break;
                    matchLen++;
                }

                if (matchLen >= minLength)
                {
                    yield return (i, i + matchLen - 1, j, j + matchLen - 1, matchLen);
                    // Skip overlapping matches
                    i += matchLen - 1;
                    break;
                }
            }
        }
    }

    /// <summary>
    /// Calculates the probability of a structure forming based on partition function (simplified).
    /// </summary>
    public static double CalculateStructureProbability(double structureEnergy, double ensembleEnergy, double temperature = 310.15)
    {
        const double R = 1.987; // cal/(mol·K)
        double RT = R * temperature / 1000.0; // kcal/mol

        double boltzmann = Math.Exp(-structureEnergy / RT);
        double partition = Math.Exp(-ensembleEnergy / RT);

        return partition > 0 ? boltzmann / partition : 0;
    }

    /// <summary>
    /// Generates a random RNA sequence.
    /// </summary>
    public static string GenerateRandomRna(int length, double gcContent = 0.5)
    {
        var random = new Random();
        var sb = new StringBuilder(length);

        for (int i = 0; i < length; i++)
        {
            double r = random.NextDouble();
            if (r < gcContent / 2)
                sb.Append('G');
            else if (r < gcContent)
                sb.Append('C');
            else if (r < gcContent + (1 - gcContent) / 2)
                sb.Append('A');
            else
                sb.Append('U');
        }

        return sb.ToString();
    }

    #endregion
}
