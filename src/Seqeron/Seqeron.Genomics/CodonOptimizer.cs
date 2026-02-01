using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics;

/// <summary>
/// Provides algorithms for codon optimization and sequence design for heterologous expression.
/// </summary>
public static class CodonOptimizer
{
    #region Records and Types

    /// <summary>
    /// Represents an organism's codon usage table.
    /// </summary>
    public readonly record struct CodonUsageTable(
        string OrganismName,
        IReadOnlyDictionary<string, double> CodonFrequencies,
        IReadOnlyDictionary<string, string> CodonToAminoAcid);

    /// <summary>
    /// Result of codon optimization.
    /// </summary>
    public readonly record struct OptimizationResult(
        string OriginalSequence,
        string OptimizedSequence,
        string ProteinSequence,
        double OriginalCAI,
        double OptimizedCAI,
        double GcContentOriginal,
        double GcContentOptimized,
        int ChangedCodons,
        IReadOnlyList<(int Position, string Original, string Optimized)> Changes);

    /// <summary>
    /// Optimization strategy options.
    /// </summary>
    public enum OptimizationStrategy
    {
        MaximizeCAI,           // Use most frequent codons
        BalancedOptimization,  // Balance CAI with other factors
        HarmonizeExpression,   // Match host codon usage distribution
        MinimizeSecondary,     // Avoid mRNA secondary structures
        AvoidRareCodeons       // Only replace rare codons
    }

    #endregion

    #region Standard Genetic Code

    private static readonly Dictionary<string, string> StandardGeneticCode = new()
    {
        // Phenylalanine
        { "UUU", "F" }, { "UUC", "F" },
        // Leucine
        { "UUA", "L" }, { "UUG", "L" }, { "CUU", "L" }, { "CUC", "L" }, { "CUA", "L" }, { "CUG", "L" },
        // Isoleucine
        { "AUU", "I" }, { "AUC", "I" }, { "AUA", "I" },
        // Methionine (Start)
        { "AUG", "M" },
        // Valine
        { "GUU", "V" }, { "GUC", "V" }, { "GUA", "V" }, { "GUG", "V" },
        // Serine
        { "UCU", "S" }, { "UCC", "S" }, { "UCA", "S" }, { "UCG", "S" }, { "AGU", "S" }, { "AGC", "S" },
        // Proline
        { "CCU", "P" }, { "CCC", "P" }, { "CCA", "P" }, { "CCG", "P" },
        // Threonine
        { "ACU", "T" }, { "ACC", "T" }, { "ACA", "T" }, { "ACG", "T" },
        // Alanine
        { "GCU", "A" }, { "GCC", "A" }, { "GCA", "A" }, { "GCG", "A" },
        // Tyrosine
        { "UAU", "Y" }, { "UAC", "Y" },
        // Stop codons
        { "UAA", "*" }, { "UAG", "*" }, { "UGA", "*" },
        // Histidine
        { "CAU", "H" }, { "CAC", "H" },
        // Glutamine
        { "CAA", "Q" }, { "CAG", "Q" },
        // Asparagine
        { "AAU", "N" }, { "AAC", "N" },
        // Lysine
        { "AAA", "K" }, { "AAG", "K" },
        // Aspartic acid
        { "GAU", "D" }, { "GAC", "D" },
        // Glutamic acid
        { "GAA", "E" }, { "GAG", "E" },
        // Cysteine
        { "UGU", "C" }, { "UGC", "C" },
        // Tryptophan
        { "UGG", "W" },
        // Arginine
        { "CGU", "R" }, { "CGC", "R" }, { "CGA", "R" }, { "CGG", "R" }, { "AGA", "R" }, { "AGG", "R" },
        // Glycine
        { "GGU", "G" }, { "GGC", "G" }, { "GGA", "G" }, { "GGG", "G" }
    };

    private static readonly Dictionary<string, List<string>> AminoAcidToCodons;

    static CodonOptimizer()
    {
        AminoAcidToCodons = StandardGeneticCode
            .GroupBy(kv => kv.Value)
            .ToDictionary(g => g.Key, g => g.Select(kv => kv.Key).ToList());
    }

    #endregion

    #region Predefined Codon Usage Tables

    /// <summary>
    /// E. coli K12 codon usage frequencies.
    /// </summary>
    public static readonly CodonUsageTable EColiK12 = new(
        "Escherichia coli K12",
        new Dictionary<string, double>
        {
            // High frequency codons for E. coli
            { "UUU", 0.58 }, { "UUC", 0.42 },
            { "UUA", 0.14 }, { "UUG", 0.13 }, { "CUU", 0.12 }, { "CUC", 0.10 }, { "CUA", 0.04 }, { "CUG", 0.47 },
            { "AUU", 0.49 }, { "AUC", 0.39 }, { "AUA", 0.11 },
            { "AUG", 1.00 },
            { "GUU", 0.28 }, { "GUC", 0.20 }, { "GUA", 0.17 }, { "GUG", 0.35 },
            { "UCU", 0.17 }, { "UCC", 0.15 }, { "UCA", 0.14 }, { "UCG", 0.14 }, { "AGU", 0.16 }, { "AGC", 0.25 },
            { "CCU", 0.18 }, { "CCC", 0.13 }, { "CCA", 0.20 }, { "CCG", 0.49 },
            { "ACU", 0.19 }, { "ACC", 0.40 }, { "ACA", 0.17 }, { "ACG", 0.25 },
            { "GCU", 0.18 }, { "GCC", 0.26 }, { "GCA", 0.23 }, { "GCG", 0.33 },
            { "UAU", 0.59 }, { "UAC", 0.41 },
            { "UAA", 0.61 }, { "UAG", 0.09 }, { "UGA", 0.30 },
            { "CAU", 0.57 }, { "CAC", 0.43 },
            { "CAA", 0.34 }, { "CAG", 0.66 },
            { "AAU", 0.49 }, { "AAC", 0.51 },
            { "AAA", 0.74 }, { "AAG", 0.26 },
            { "GAU", 0.63 }, { "GAC", 0.37 },
            { "GAA", 0.68 }, { "GAG", 0.32 },
            { "UGU", 0.46 }, { "UGC", 0.54 },
            { "UGG", 1.00 },
            { "CGU", 0.36 }, { "CGC", 0.36 }, { "CGA", 0.07 }, { "CGG", 0.11 }, { "AGA", 0.07 }, { "AGG", 0.04 },
            { "GGU", 0.35 }, { "GGC", 0.37 }, { "GGA", 0.13 }, { "GGG", 0.15 }
        },
        StandardGeneticCode);

    /// <summary>
    /// Saccharomyces cerevisiae (yeast) codon usage frequencies.
    /// </summary>
    public static readonly CodonUsageTable Yeast = new(
        "Saccharomyces cerevisiae",
        new Dictionary<string, double>
        {
            { "UUU", 0.59 }, { "UUC", 0.41 },
            { "UUA", 0.28 }, { "UUG", 0.29 }, { "CUU", 0.13 }, { "CUC", 0.06 }, { "CUA", 0.14 }, { "CUG", 0.11 },
            { "AUU", 0.46 }, { "AUC", 0.26 }, { "AUA", 0.27 },
            { "AUG", 1.00 },
            { "GUU", 0.39 }, { "GUC", 0.21 }, { "GUA", 0.21 }, { "GUG", 0.19 },
            { "UCU", 0.26 }, { "UCC", 0.16 }, { "UCA", 0.21 }, { "UCG", 0.10 }, { "AGU", 0.16 }, { "AGC", 0.11 },
            { "CCU", 0.31 }, { "CCC", 0.15 }, { "CCA", 0.42 }, { "CCG", 0.12 },
            { "ACU", 0.35 }, { "ACC", 0.22 }, { "ACA", 0.30 }, { "ACG", 0.14 },
            { "GCU", 0.38 }, { "GCC", 0.22 }, { "GCA", 0.29 }, { "GCG", 0.11 },
            { "UAU", 0.56 }, { "UAC", 0.44 },
            { "UAA", 0.48 }, { "UAG", 0.24 }, { "UGA", 0.29 },
            { "CAU", 0.64 }, { "CAC", 0.36 },
            { "CAA", 0.69 }, { "CAG", 0.31 },
            { "AAU", 0.59 }, { "AAC", 0.41 },
            { "AAA", 0.58 }, { "AAG", 0.42 },
            { "GAU", 0.65 }, { "GAC", 0.35 },
            { "GAA", 0.70 }, { "GAG", 0.30 },
            { "UGU", 0.63 }, { "UGC", 0.37 },
            { "UGG", 1.00 },
            { "CGU", 0.14 }, { "CGC", 0.06 }, { "CGA", 0.07 }, { "CGG", 0.04 }, { "AGA", 0.48 }, { "AGG", 0.21 },
            { "GGU", 0.47 }, { "GGC", 0.19 }, { "GGA", 0.22 }, { "GGG", 0.12 }
        },
        StandardGeneticCode);

    /// <summary>
    /// Human codon usage frequencies.
    /// </summary>
    public static readonly CodonUsageTable Human = new(
        "Homo sapiens",
        new Dictionary<string, double>
        {
            { "UUU", 0.45 }, { "UUC", 0.55 },
            { "UUA", 0.07 }, { "UUG", 0.13 }, { "CUU", 0.13 }, { "CUC", 0.20 }, { "CUA", 0.07 }, { "CUG", 0.41 },
            { "AUU", 0.36 }, { "AUC", 0.48 }, { "AUA", 0.16 },
            { "AUG", 1.00 },
            { "GUU", 0.18 }, { "GUC", 0.24 }, { "GUA", 0.11 }, { "GUG", 0.47 },
            { "UCU", 0.18 }, { "UCC", 0.22 }, { "UCA", 0.15 }, { "UCG", 0.06 }, { "AGU", 0.15 }, { "AGC", 0.24 },
            { "CCU", 0.28 }, { "CCC", 0.33 }, { "CCA", 0.27 }, { "CCG", 0.11 },
            { "ACU", 0.24 }, { "ACC", 0.36 }, { "ACA", 0.28 }, { "ACG", 0.12 },
            { "GCU", 0.26 }, { "GCC", 0.40 }, { "GCA", 0.23 }, { "GCG", 0.11 },
            { "UAU", 0.43 }, { "UAC", 0.57 },
            { "UAA", 0.28 }, { "UAG", 0.20 }, { "UGA", 0.52 },
            { "CAU", 0.41 }, { "CAC", 0.59 },
            { "CAA", 0.25 }, { "CAG", 0.75 },
            { "AAU", 0.46 }, { "AAC", 0.54 },
            { "AAA", 0.42 }, { "AAG", 0.58 },
            { "GAU", 0.46 }, { "GAC", 0.54 },
            { "GAA", 0.42 }, { "GAG", 0.58 },
            { "UGU", 0.45 }, { "UGC", 0.55 },
            { "UGG", 1.00 },
            { "CGU", 0.08 }, { "CGC", 0.19 }, { "CGA", 0.11 }, { "CGG", 0.21 }, { "AGA", 0.20 }, { "AGG", 0.20 },
            { "GGU", 0.16 }, { "GGC", 0.34 }, { "GGA", 0.25 }, { "GGG", 0.25 }
        },
        StandardGeneticCode);

    #endregion

    #region Codon Optimization

    /// <summary>
    /// Optimizes a coding sequence for expression in a target organism.
    /// </summary>
    public static OptimizationResult OptimizeSequence(
        string codingSequence,
        CodonUsageTable targetOrganism,
        OptimizationStrategy strategy = OptimizationStrategy.BalancedOptimization,
        double gcTargetMin = 0.40,
        double gcTargetMax = 0.60)
    {
        if (string.IsNullOrEmpty(codingSequence))
        {
            return new OptimizationResult("", "", "", 0, 0, 0, 0, 0, new List<(int, string, string)>());
        }

        string rna = codingSequence.ToUpperInvariant().Replace('T', 'U');

        if (rna.Length % 3 != 0)
        {
            // Trim to complete codons
            rna = rna.Substring(0, (rna.Length / 3) * 3);
        }

        var originalCodons = SplitIntoCodons(rna);
        var optimizedCodons = new List<string>();
        var changes = new List<(int Position, string Original, string Optimized)>();

        double originalCAI = CalculateCAI(rna, targetOrganism);
        var proteinBuilder = new StringBuilder();

        for (int i = 0; i < originalCodons.Count; i++)
        {
            string codon = originalCodons[i];
            string aminoAcid = TranslateCodon(codon);
            proteinBuilder.Append(aminoAcid);

            if (aminoAcid == "*")
            {
                // Keep stop codon
                optimizedCodons.Add(codon);
                continue;
            }

            string optimizedCodon = SelectOptimalCodon(aminoAcid, codon, targetOrganism, strategy);

            if (optimizedCodon != codon)
            {
                changes.Add((i * 3, codon, optimizedCodon));
            }

            optimizedCodons.Add(optimizedCodon);
        }

        string optimizedSequence = string.Join("", optimizedCodons);

        // Apply GC content balancing if needed
        if (strategy == OptimizationStrategy.BalancedOptimization)
        {
            optimizedSequence = BalanceGcContent(optimizedSequence, originalCodons, targetOrganism, gcTargetMin, gcTargetMax);
        }

        double optimizedCAI = CalculateCAI(optimizedSequence, targetOrganism);

        return new OptimizationResult(
            OriginalSequence: rna,
            OptimizedSequence: optimizedSequence,
            ProteinSequence: proteinBuilder.ToString(),
            OriginalCAI: originalCAI,
            OptimizedCAI: optimizedCAI,
            GcContentOriginal: CalculateGcContent(rna),
            GcContentOptimized: CalculateGcContent(optimizedSequence),
            ChangedCodons: changes.Count,
            Changes: changes);
    }

    private static string SelectOptimalCodon(string aminoAcid, string currentCodon, CodonUsageTable table, OptimizationStrategy strategy)
    {
        if (!AminoAcidToCodons.TryGetValue(aminoAcid, out var synonymousCodons))
            return currentCodon;

        if (synonymousCodons.Count == 1)
            return synonymousCodons[0];

        switch (strategy)
        {
            case OptimizationStrategy.MaximizeCAI:
                return synonymousCodons
                    .OrderByDescending(c => table.CodonFrequencies.GetValueOrDefault(c, 0))
                    .First();

            case OptimizationStrategy.AvoidRareCodeons:
                double currentFreq = table.CodonFrequencies.GetValueOrDefault(currentCodon, 0);
                if (currentFreq < 0.15)
                {
                    return synonymousCodons
                        .Where(c => table.CodonFrequencies.GetValueOrDefault(c, 0) >= 0.2)
                        .OrderByDescending(c => table.CodonFrequencies.GetValueOrDefault(c, 0))
                        .FirstOrDefault() ?? currentCodon;
                }
                return currentCodon;

            case OptimizationStrategy.HarmonizeExpression:
                // Use weighted random selection based on frequencies
                return SelectWeightedCodon(synonymousCodons, table);

            case OptimizationStrategy.BalancedOptimization:
            default:
                // Use codons with frequency >= 0.2, prefer highest
                var goodCodons = synonymousCodons
                    .Where(c => table.CodonFrequencies.GetValueOrDefault(c, 0) >= 0.15)
                    .OrderByDescending(c => table.CodonFrequencies.GetValueOrDefault(c, 0))
                    .ToList();
                return goodCodons.Count > 0 ? goodCodons[0] : currentCodon;
        }
    }

    private static string SelectWeightedCodon(List<string> codons, CodonUsageTable table)
    {
        var random = new Random();
        double totalWeight = codons.Sum(c => table.CodonFrequencies.GetValueOrDefault(c, 0.01));
        double r = random.NextDouble() * totalWeight;

        double cumulative = 0;
        foreach (var codon in codons)
        {
            cumulative += table.CodonFrequencies.GetValueOrDefault(codon, 0.01);
            if (r <= cumulative)
                return codon;
        }

        return codons[0];
    }

    private static string BalanceGcContent(string sequence, List<string> originalCodons, CodonUsageTable table, double minGc, double maxGc)
    {
        double currentGc = CalculateGcContent(sequence);

        if (currentGc >= minGc && currentGc <= maxGc)
            return sequence;

        var codons = SplitIntoCodons(sequence);
        bool needMoreGc = currentGc < minGc;

        for (int i = 0; i < codons.Count && (currentGc < minGc || currentGc > maxGc); i++)
        {
            string aminoAcid = TranslateCodon(codons[i]);
            if (aminoAcid == "*") continue;

            if (!AminoAcidToCodons.TryGetValue(aminoAcid, out var alternatives))
                continue;

            // Find alternative with appropriate GC content
            var sorted = needMoreGc
                ? alternatives.OrderByDescending(c => GetCodonGcContent(c))
                : alternatives.OrderBy(c => GetCodonGcContent(c));

            foreach (var alt in sorted)
            {
                if (table.CodonFrequencies.GetValueOrDefault(alt, 0) >= 0.1)
                {
                    codons[i] = alt;
                    break;
                }
            }

            currentGc = CalculateGcContent(string.Join("", codons));
        }

        return string.Join("", codons);
    }

    #endregion

    #region CAI Calculation

    /// <summary>
    /// Calculates the Codon Adaptation Index (CAI) for a sequence.
    /// </summary>
    public static double CalculateCAI(string codingSequence, CodonUsageTable table)
    {
        if (string.IsNullOrEmpty(codingSequence))
            return 0;

        string rna = codingSequence.ToUpperInvariant().Replace('T', 'U');
        var codons = SplitIntoCodons(rna);

        if (codons.Count == 0)
            return 0;

        double logSum = 0;
        int count = 0;

        foreach (var codon in codons)
        {
            string aminoAcid = TranslateCodon(codon);
            if (aminoAcid == "*") continue;

            double w = CalculateRelativeAdaptiveness(codon, aminoAcid, table);
            if (w > 0)
            {
                logSum += Math.Log(w);
                count++;
            }
        }

        return count > 0 ? Math.Exp(logSum / count) : 0;
    }

    private static double CalculateRelativeAdaptiveness(string codon, string aminoAcid, CodonUsageTable table)
    {
        if (!AminoAcidToCodons.TryGetValue(aminoAcid, out var synonymousCodons))
            return 1;

        double codonFreq = table.CodonFrequencies.GetValueOrDefault(codon, 0.01);
        double maxFreq = synonymousCodons.Max(c => table.CodonFrequencies.GetValueOrDefault(c, 0.01));

        return maxFreq > 0 ? codonFreq / maxFreq : 1;
    }

    #endregion

    #region Sequence Modification

    /// <summary>
    /// Removes restriction enzyme recognition sites from a sequence while preserving the protein.
    /// </summary>
    public static string RemoveRestrictionSites(string codingSequence, IEnumerable<string> restrictionSites, CodonUsageTable table)
    {
        if (string.IsNullOrEmpty(codingSequence))
            return "";

        string rna = codingSequence.ToUpperInvariant().Replace('T', 'U');
        var codons = SplitIntoCodons(rna);

        foreach (var site in restrictionSites)
        {
            string siteRna = site.ToUpperInvariant().Replace('T', 'U');
            string current = string.Join("", codons);

            while (current.Contains(siteRna))
            {
                int pos = current.IndexOf(siteRna, StringComparison.Ordinal);
                int codonIdx = pos / 3;

                // Try to change one of the codons overlapping the site
                for (int i = codonIdx; i <= Math.Min(codonIdx + 2, codons.Count - 1); i++)
                {
                    string aa = TranslateCodon(codons[i]);
                    if (aa == "*") continue;

                    if (AminoAcidToCodons.TryGetValue(aa, out var alts))
                    {
                        foreach (var alt in alts.Where(a => a != codons[i]))
                        {
                            var testCodons = new List<string>(codons) { [i] = alt };
                            string testSeq = string.Join("", testCodons);
                            if (!testSeq.Contains(siteRna))
                            {
                                codons[i] = alt;
                                break;
                            }
                        }
                    }
                }

                current = string.Join("", codons);

                // Prevent infinite loop
                if (current.Contains(siteRna))
                    break;
            }
        }

        return string.Join("", codons);
    }

    /// <summary>
    /// Reduces mRNA secondary structure by avoiding self-complementary regions.
    /// </summary>
    public static string ReduceSecondaryStructure(string codingSequence, CodonUsageTable table, int windowSize = 40)
    {
        if (string.IsNullOrEmpty(codingSequence) || codingSequence.Length < windowSize)
            return codingSequence;

        string rna = codingSequence.ToUpperInvariant().Replace('T', 'U');
        var codons = SplitIntoCodons(rna);

        for (int i = 0; i < codons.Count - windowSize / 3; i++)
        {
            string window = string.Join("", codons.Skip(i).Take(windowSize / 3 + 1));
            double structureScore = CalculateLocalStructure(window);

            if (structureScore > 0.5) // High structure propensity
            {
                // Try to reduce by changing codons
                for (int j = i; j < Math.Min(i + windowSize / 3, codons.Count); j++)
                {
                    string aa = TranslateCodon(codons[j]);
                    if (aa == "*") continue;

                    if (AminoAcidToCodons.TryGetValue(aa, out var alts))
                    {
                        string bestAlt = codons[j];
                        double bestScore = structureScore;

                        foreach (var alt in alts)
                        {
                            var testCodons = new List<string>(codons) { [j] = alt };
                            string testWindow = string.Join("", testCodons.Skip(i).Take(windowSize / 3 + 1));
                            double testScore = CalculateLocalStructure(testWindow);

                            if (testScore < bestScore && table.CodonFrequencies.GetValueOrDefault(alt, 0) >= 0.1)
                            {
                                bestScore = testScore;
                                bestAlt = alt;
                            }
                        }

                        codons[j] = bestAlt;
                    }
                }
            }
        }

        return string.Join("", codons);
    }

    private static double CalculateLocalStructure(string sequence)
    {
        int complementaryPairs = 0;
        int n = sequence.Length;

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 4; j < n; j++)
            {
                if (AreComplementary(sequence[i], sequence[j]))
                    complementaryPairs++;
            }
        }

        double maxPairs = (n * (n - 4)) / 2.0;
        return maxPairs > 0 ? complementaryPairs / maxPairs : 0;
    }

    private static bool AreComplementary(char b1, char b2)
    {
        return (b1 == 'A' && b2 == 'U') || (b1 == 'U' && b2 == 'A') ||
               (b1 == 'G' && b2 == 'C') || (b1 == 'C' && b2 == 'G');
    }

    #endregion

    #region Analysis Functions

    /// <summary>
    /// Analyzes rare codon usage in a sequence.
    /// </summary>
    public static IEnumerable<(int Position, string Codon, string AminoAcid, double Frequency)> FindRareCodons(
        string codingSequence,
        CodonUsageTable table,
        double threshold = 0.15)
    {
        if (string.IsNullOrEmpty(codingSequence))
            yield break;

        string rna = codingSequence.ToUpperInvariant().Replace('T', 'U');
        var codons = SplitIntoCodons(rna);

        for (int i = 0; i < codons.Count; i++)
        {
            double freq = table.CodonFrequencies.GetValueOrDefault(codons[i], 0);
            if (freq < threshold)
            {
                string aa = TranslateCodon(codons[i]);
                yield return (i * 3, codons[i], aa, freq);
            }
        }
    }

    /// <summary>
    /// Calculates codon frequency distribution for a sequence.
    /// </summary>
    public static Dictionary<string, int> CalculateCodonUsage(string codingSequence)
    {
        var usage = new Dictionary<string, int>();

        if (string.IsNullOrEmpty(codingSequence))
            return usage;

        string rna = codingSequence.ToUpperInvariant().Replace('T', 'U');
        var codons = SplitIntoCodons(rna);

        foreach (var codon in codons)
        {
            if (!usage.ContainsKey(codon))
                usage[codon] = 0;
            usage[codon]++;
        }

        return usage;
    }

    /// <summary>
    /// Compares codon usage between two sequences.
    /// </summary>
    public static double CompareCodonUsage(string sequence1, string sequence2)
    {
        var usage1 = CalculateCodonUsage(sequence1);
        var usage2 = CalculateCodonUsage(sequence2);

        var allCodons = usage1.Keys.Union(usage2.Keys).ToList();
        if (allCodons.Count == 0)
            return 0;

        int total1 = usage1.Values.Sum();
        int total2 = usage2.Values.Sum();

        if (total1 == 0 || total2 == 0)
            return 0;

        double correlation = 0;
        foreach (var codon in allCodons)
        {
            double freq1 = usage1.GetValueOrDefault(codon, 0) / (double)total1;
            double freq2 = usage2.GetValueOrDefault(codon, 0) / (double)total2;
            correlation += Math.Abs(freq1 - freq2);
        }

        return 1 - (correlation / 2);
    }

    #endregion

    #region Utility Methods

    private static List<string> SplitIntoCodons(string sequence)
    {
        var codons = new List<string>();
        for (int i = 0; i + 2 < sequence.Length; i += 3)
        {
            codons.Add(sequence.Substring(i, 3));
        }
        return codons;
    }

    private static string TranslateCodon(string codon)
    {
        return StandardGeneticCode.GetValueOrDefault(codon, "X");
    }

    private static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcFractionFast();

    private static double GetCodonGcContent(string codon)
    {
        int gc = codon.Count(c => c == 'G' || c == 'C');
        return gc / 3.0;
    }

    /// <summary>
    /// Creates a custom codon usage table from a reference sequence.
    /// </summary>
    public static CodonUsageTable CreateCodonTableFromSequence(string referenceSequence, string organismName)
    {
        var usage = CalculateCodonUsage(referenceSequence);
        var frequencies = new Dictionary<string, double>();

        // Group by amino acid and calculate relative frequencies
        var byAminoAcid = usage
            .Where(kv => StandardGeneticCode.ContainsKey(kv.Key))
            .GroupBy(kv => StandardGeneticCode[kv.Key]);

        foreach (var group in byAminoAcid)
        {
            int total = group.Sum(g => g.Value);
            foreach (var codon in group)
            {
                frequencies[codon.Key] = total > 0 ? (double)codon.Value / total : 0;
            }
        }

        return new CodonUsageTable(organismName, frequencies, StandardGeneticCode);
    }

    #endregion
}
