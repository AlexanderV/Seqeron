using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.MolTools;

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
    /// One window of the %MinMax codon-usage profile of Clarke &amp; Clark (2008).
    /// </summary>
    /// <param name="WindowStartCodon">0-based codon index of the first codon in the window.</param>
    /// <param name="PercentMinMax">
    /// Signed %MinMax value for the window: positive values are %Max (a run of predominantly
    /// common codons), negative values are %Min (a run of predominantly rare codons), and 0
    /// means the window's codon usage equals the per-amino-acid average. A value of −100
    /// is a window encoded entirely with the rarest synonymous codons; +100 is one encoded
    /// entirely with the most common synonymous codons. (Clarke &amp; Clark 2008.)
    /// </param>
    public readonly record struct MinMaxWindow(int WindowStartCodon, double PercentMinMax);

    /// <summary>
    /// One rare-codon cluster (RCC) detected by the Sherlocc rule of Chartier et&#160;al. (2012).
    /// </summary>
    /// <param name="StartCodon">0-based codon index of the first codon in the cluster window.</param>
    /// <param name="EndCodon">0-based codon index of the last codon in the cluster window (inclusive).</param>
    /// <param name="RareCount">Number of rare ("pause") codons inside the window.</param>
    public readonly record struct RareCodonCluster(int StartCodon, int EndCodon, int RareCount);

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

    /// <summary>
    /// Amino acids encoded by a single codon in the standard genetic code
    /// (Methionine/AUG and Tryptophan/UGG). Their relative adaptiveness w is always 1
    /// regardless of codon usage bias, so Sharp &amp; Li (1987) / Jansen et al. (2003)
    /// exclude them from CAI to avoid skewing the geometric mean.
    /// Derived from <see cref="AminoAcidToCodons"/> (groups of size 1), not hard-coded.
    /// </summary>
    private static readonly HashSet<string> SingleCodonAminoAcids;

    static CodonOptimizer()
    {
        AminoAcidToCodons = StandardGeneticCode
            .GroupBy(kv => kv.Value)
            .ToDictionary(g => g.Key, g => g.Select(kv => kv.Key).ToList());

        // Single-codon (non-degenerate) amino acids: those with exactly one synonymous
        // codon. In the standard genetic code these are Met ("M", AUG) and Trp ("W", UGG).
        // Stop ("*") has 3 codons and is handled separately. Source: Sharp & Li (1987).
        SingleCodonAminoAcids = AminoAcidToCodons
            .Where(kv => kv.Key != "*" && kv.Value.Count == 1)
            .Select(kv => kv.Key)
            .ToHashSet();
    }

    #endregion

    #region Predefined Codon Usage Tables

    /// <summary>
    /// E. coli K12 codon usage frequencies (relative fraction per amino acid).
    /// Source: Kazusa Codon Usage Database, species=316407 (E. coli K-12 substr. W3110, 4332 CDS).
    /// URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407
    /// </summary>
    public static readonly CodonUsageTable EColiK12 = new(
        "Escherichia coli K12",
        new Dictionary<string, double>
        {
            // Phenylalanine (F)
            { "UUU", 0.57 }, { "UUC", 0.43 },
            // Leucine (L)
            { "UUA", 0.13 }, { "UUG", 0.13 }, { "CUU", 0.10 }, { "CUC", 0.10 }, { "CUA", 0.04 }, { "CUG", 0.50 },
            // Isoleucine (I)
            { "AUU", 0.51 }, { "AUC", 0.42 }, { "AUA", 0.07 },
            // Methionine (M)
            { "AUG", 1.00 },
            // Valine (V)
            { "GUU", 0.26 }, { "GUC", 0.22 }, { "GUA", 0.15 }, { "GUG", 0.37 },
            // Serine (S)
            { "UCU", 0.15 }, { "UCC", 0.15 }, { "UCA", 0.12 }, { "UCG", 0.15 }, { "AGU", 0.15 }, { "AGC", 0.28 },
            // Proline (P)
            { "CCU", 0.16 }, { "CCC", 0.12 }, { "CCA", 0.19 }, { "CCG", 0.53 },
            // Threonine (T)
            { "ACU", 0.16 }, { "ACC", 0.44 }, { "ACA", 0.13 }, { "ACG", 0.27 },
            // Alanine (A)
            { "GCU", 0.16 }, { "GCC", 0.27 }, { "GCA", 0.21 }, { "GCG", 0.36 },
            // Tyrosine (Y)
            { "UAU", 0.57 }, { "UAC", 0.43 },
            // Stop (*)
            { "UAA", 0.64 }, { "UAG", 0.07 }, { "UGA", 0.29 },
            // Histidine (H)
            { "CAU", 0.57 }, { "CAC", 0.43 },
            // Glutamine (Q)
            { "CAA", 0.35 }, { "CAG", 0.65 },
            // Asparagine (N)
            { "AAU", 0.45 }, { "AAC", 0.55 },
            // Lysine (K)
            { "AAA", 0.76 }, { "AAG", 0.24 },
            // Aspartic acid (D)
            { "GAU", 0.63 }, { "GAC", 0.37 },
            // Glutamic acid (E)
            { "GAA", 0.69 }, { "GAG", 0.31 },
            // Cysteine (C)
            { "UGU", 0.44 }, { "UGC", 0.56 },
            // Tryptophan (W)
            { "UGG", 1.00 },
            // Arginine (R)
            { "CGU", 0.38 }, { "CGC", 0.40 }, { "CGA", 0.06 }, { "CGG", 0.10 }, { "AGA", 0.04 }, { "AGG", 0.02 },
            // Glycine (G)
            { "GGU", 0.34 }, { "GGC", 0.41 }, { "GGA", 0.11 }, { "GGG", 0.15 }
        },
        StandardGeneticCode);

    /// <summary>
    /// Saccharomyces cerevisiae (yeast) codon usage frequencies (relative fraction per amino acid).
    /// Source: Kazusa Codon Usage Database, species=4932.
    /// URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932
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
            { "UAA", 0.47 }, { "UAG", 0.23 }, { "UGA", 0.30 },
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
    /// Human codon usage frequencies (relative fraction per amino acid).
    /// Source: Kazusa Codon Usage Database, species=9606.
    /// URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
    /// </summary>
    public static readonly CodonUsageTable Human = new(
        "Homo sapiens",
        new Dictionary<string, double>
        {
            { "UUU", 0.46 }, { "UUC", 0.54 },
            { "UUA", 0.08 }, { "UUG", 0.13 }, { "CUU", 0.13 }, { "CUC", 0.20 }, { "CUA", 0.07 }, { "CUG", 0.40 },
            { "AUU", 0.36 }, { "AUC", 0.47 }, { "AUA", 0.17 },
            { "AUG", 1.00 },
            { "GUU", 0.18 }, { "GUC", 0.24 }, { "GUA", 0.12 }, { "GUG", 0.46 },
            { "UCU", 0.19 }, { "UCC", 0.22 }, { "UCA", 0.15 }, { "UCG", 0.05 }, { "AGU", 0.15 }, { "AGC", 0.24 },
            { "CCU", 0.29 }, { "CCC", 0.32 }, { "CCA", 0.28 }, { "CCG", 0.11 },
            { "ACU", 0.25 }, { "ACC", 0.36 }, { "ACA", 0.28 }, { "ACG", 0.11 },
            { "GCU", 0.27 }, { "GCC", 0.40 }, { "GCA", 0.23 }, { "GCG", 0.11 },
            { "UAU", 0.44 }, { "UAC", 0.56 },
            { "UAA", 0.30 }, { "UAG", 0.24 }, { "UGA", 0.47 },
            { "CAU", 0.42 }, { "CAC", 0.58 },
            { "CAA", 0.27 }, { "CAG", 0.73 },
            { "AAU", 0.47 }, { "AAC", 0.53 },
            { "AAA", 0.43 }, { "AAG", 0.57 },
            { "GAU", 0.46 }, { "GAC", 0.54 },
            { "GAA", 0.42 }, { "GAG", 0.58 },
            { "UGU", 0.46 }, { "UGC", 0.54 },
            { "UGG", 1.00 },
            { "CGU", 0.08 }, { "CGC", 0.18 }, { "CGA", 0.11 }, { "CGG", 0.20 }, { "AGA", 0.21 }, { "AGG", 0.21 },
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
        double gcTargetMax = 0.60,
        double rareCodonThreshold = 0.15)
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

            string optimizedCodon = SelectOptimalCodon(aminoAcid, codon, targetOrganism, strategy, rareCodonThreshold);

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

            // Rebuild changes to reflect GC balancing modifications
            changes.Clear();
            var finalCodons = SplitIntoCodons(optimizedSequence);
            for (int i = 0; i < originalCodons.Count && i < finalCodons.Count; i++)
            {
                if (originalCodons[i] != finalCodons[i])
                    changes.Add((i * 3, originalCodons[i], finalCodons[i]));
            }
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

    private static string SelectOptimalCodon(string aminoAcid, string currentCodon, CodonUsageTable table, OptimizationStrategy strategy, double rareCodonThreshold)
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
                if (currentFreq < rareCodonThreshold)
                {
                    return synonymousCodons
                        .Where(c => table.CodonFrequencies.GetValueOrDefault(c, 0) >= rareCodonThreshold)
                        .OrderByDescending(c => table.CodonFrequencies.GetValueOrDefault(c, 0))
                        .FirstOrDefault() ?? currentCodon;
                }
                return currentCodon;

            case OptimizationStrategy.HarmonizeExpression:
                // Use weighted random selection based on frequencies
                return SelectWeightedCodon(synonymousCodons, table);

            case OptimizationStrategy.BalancedOptimization:
            default:
                var goodCodons = synonymousCodons
                    .Where(c => table.CodonFrequencies.GetValueOrDefault(c, 0) >= rareCodonThreshold)
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
    /// Calculates the Codon Adaptation Index (CAI) for a sequence
    /// (Sharp &amp; Li 1987, <c>CAI = (∏ w_i)^(1/L)</c>, the geometric mean of the relative
    /// adaptiveness <c>w_i = f_i / max(f_j)</c> over the gene's codons; stop codons excluded).
    /// </summary>
    /// <param name="codingSequence">Coding sequence (DNA or RNA; case-insensitive).</param>
    /// <param name="table">Reference codon usage table.</param>
    /// <param name="excludeSingleCodonAminoAcids">
    /// When <see langword="true"/>, codons of amino acids that have a single codon in the
    /// standard genetic code (Met/AUG, Trp/UGG) are excluded from the geometric mean, as the
    /// original Sharp &amp; Li (1987) definition prescribes and Jansen et al. (2003) reiterate:
    /// "codon families containing a single codon (e.g. AUG and UGG …) should be excluded in
    /// computing CAI" because their w is always 1 regardless of bias. Default <see langword="false"/>
    /// preserves the historical inclusive behaviour (these codons counted with w = 1.0).
    /// </param>
    public static double CalculateCAI(string codingSequence, CodonUsageTable table, bool excludeSingleCodonAminoAcids = false)
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

            // Per Sharp & Li (1987): single-codon amino acids (Met/AUG, Trp/UGG) are excluded
            // from CAI when requested, since their w is always 1 and would skew the geometric mean.
            if (excludeSingleCodonAminoAcids && SingleCodonAminoAcids.Contains(aminoAcid)) continue;

            double w = CalculateRelativeAdaptiveness(codon, aminoAcid, table);
            if (double.IsNaN(w)) continue; // No frequency data for this AA in table

            logSum += Math.Log(w);
            count++;
        }

        return count > 0 ? Math.Exp(logSum / count) : 0;
    }

    private static double CalculateRelativeAdaptiveness(string codon, string aminoAcid, CodonUsageTable table)
    {
        if (!AminoAcidToCodons.TryGetValue(aminoAcid, out var synonymousCodons))
            return double.NaN; // Not a standard amino acid — no adaptiveness data

        double codonFreq = table.CodonFrequencies.GetValueOrDefault(codon, 0);
        double maxFreq = synonymousCodons.Max(c => table.CodonFrequencies.GetValueOrDefault(c, 0));

        if (maxFreq <= 0)
            return double.NaN; // No frequency data for this amino acid in the table

        // Clamp to 1e-6 to avoid ln(0) when codon is absent from an incomplete custom table
        // but other synonymous codons are present (maxFreq > 0, codonFreq = 0).
        // Sharp & Li (1987) did not encounter this case (complete reference sets),
        // but real-world partial tables may have gaps.
        return Math.Max(codonFreq / maxFreq, 1e-6);
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

    // Default sliding-window width (in codons) for the %MinMax profile.
    // "%MinMax results are typically averaged over an 18-codon sliding window."
    // Clarke TF, Clark PL (2008) "Rare Codons Cluster", PLoS ONE 3(10):e3412.
    private const int DefaultMinMaxWindowCodons = 18;

    // Sherlocc rare-codon-cluster (RCC) detection parameters.
    // "a seven position-wide window ... containing at least four pause positions out of seven."
    // Chartier M, Gaudreault F, Najmanovich R (2012) Bioinformatics 28(11):1438-1445,
    // doi:10.1093/bioinformatics/bts149.
    private const int DefaultClusterWindowCodons = 7;
    private const int DefaultClusterMinRareCodons = 4;

    /// <summary>
    /// Computes the %MinMax codon-usage profile of Clarke &amp; Clark (2008) over a sliding window.
    /// </summary>
    /// <remarks>
    /// For each amino acid <c>i</c> with <c>n</c> synonymous codons, let <c>Xij</c> be the usage
    /// frequency of the codon actually used, <c>Xmax,i</c> / <c>Xmin,i</c> the usage frequencies of
    /// the most / least common synonymous codon, and <c>Xavg,i</c> the arithmetic mean of the
    /// synonymous codon frequencies. Over a window of <paramref name="windowSize"/> codons:
    /// if Σ Xij &gt; Σ Xavg,i the window yields %Max = Σ(Xij − Xavg,i) / Σ(Xmax,i − Xavg,i) × 100
    /// (returned as a positive value); if Σ Xij &lt; Σ Xavg,i it yields %Min =
    /// Σ(Xavg,i − Xij) / Σ(Xavg,i − Xmin,i) × 100 (returned as a negative value). Codons of
    /// single-codon amino acids and unknown / stop codons (no synonymous spread) contribute 0 to
    /// both numerator and denominator. Source: Clarke &amp; Clark (2008), PLoS ONE 3(10):e3412.
    /// </remarks>
    /// <param name="codingSequence">DNA or RNA coding sequence (T is normalised to U).</param>
    /// <param name="table">Reference codon-usage table (per-amino-acid relative fractions).</param>
    /// <param name="windowSize">Sliding-window width in codons (default 18, per Clarke &amp; Clark 2008).</param>
    /// <returns>
    /// One <see cref="MinMaxWindow"/> per window position (codon indices
    /// 0 .. codonCount − windowSize). Empty if the sequence has fewer than
    /// <paramref name="windowSize"/> complete codons.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">windowSize &lt; 1.</exception>
    public static IReadOnlyList<MinMaxWindow> CalculateMinMaxProfile(
        string codingSequence,
        CodonUsageTable table,
        int windowSize = DefaultMinMaxWindowCodons)
    {
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1 codon.");

        var profile = new List<MinMaxWindow>();
        if (string.IsNullOrEmpty(codingSequence))
            return profile;

        string rna = codingSequence.ToUpperInvariant().Replace('T', 'U');
        var codons = SplitIntoCodons(rna);
        if (codons.Count < windowSize)
            return profile;

        // Per-codon (Xij), per-family average (Xavg), max (Xmax) and min (Xmin) frequencies.
        var xij = new double[codons.Count];
        var xavg = new double[codons.Count];
        var xmax = new double[codons.Count];
        var xmin = new double[codons.Count];
        for (int i = 0; i < codons.Count; i++)
        {
            string codon = codons[i];
            xij[i] = table.CodonFrequencies.GetValueOrDefault(codon, 0);
            string aa = TranslateCodon(codon);

            if (AminoAcidToCodons.TryGetValue(aa, out var synonyms) && synonyms.Count > 0)
            {
                double sum = 0, max = double.MinValue, min = double.MaxValue;
                foreach (var syn in synonyms)
                {
                    double f = table.CodonFrequencies.GetValueOrDefault(syn, 0);
                    sum += f;
                    if (f > max) max = f;
                    if (f < min) min = f;
                }
                xavg[i] = sum / synonyms.Count;
                xmax[i] = max;
                xmin[i] = min;
            }
            else
            {
                // Unknown codon: no synonymous family — contributes nothing to either side.
                xavg[i] = xij[i];
                xmax[i] = xij[i];
                xmin[i] = xij[i];
            }
        }

        for (int start = 0; start + windowSize <= codons.Count; start++)
        {
            double sumXij = 0, sumXavg = 0, sumMaxDelta = 0, sumMinDelta = 0;
            for (int k = start; k < start + windowSize; k++)
            {
                sumXij += xij[k];
                sumXavg += xavg[k];
                sumMaxDelta += xmax[k] - xavg[k];
                sumMinDelta += xavg[k] - xmin[k];
            }

            double percent;
            if (sumXij > sumXavg)
            {
                // %Max — positive value.
                percent = sumMaxDelta > 0 ? (sumXij - sumXavg) / sumMaxDelta * 100.0 : 0.0;
            }
            else if (sumXij < sumXavg)
            {
                // %Min — returned as a negative value (rare-codon side).
                percent = sumMinDelta > 0 ? -((sumXavg - sumXij) / sumMinDelta * 100.0) : 0.0;
            }
            else
            {
                percent = 0.0;
            }

            profile.Add(new MinMaxWindow(start, percent));
        }

        return profile;
    }

    /// <summary>
    /// Detects rare-codon clusters (RCCs) using the Sherlocc rule of Chartier et&#160;al. (2012):
    /// a window of <paramref name="windowSize"/> codons is a cluster when it contains at least
    /// <paramref name="minRareCodons"/> rare ("pause") codons.
    /// </summary>
    /// <remarks>
    /// A codon is "rare"/"pause" when its usage frequency in <paramref name="table"/> is strictly
    /// below <paramref name="rareThreshold"/> — the same per-codon criterion as
    /// <see cref="FindRareCodons"/>. Overlapping windows are merged into maximal clusters so a long
    /// rare run is reported once. This is opt-in; <see cref="FindRareCodons"/> (per-codon) is
    /// unchanged. Defaults reproduce the published Sherlocc rule "a seven position-wide window …
    /// containing at least four pause positions out of seven" (Chartier et&#160;al. 2012,
    /// doi:10.1093/bioinformatics/bts149).
    /// </remarks>
    /// <param name="codingSequence">DNA or RNA coding sequence (T is normalised to U).</param>
    /// <param name="table">Reference codon-usage table (per-amino-acid relative fractions).</param>
    /// <param name="rareThreshold">Per-codon rare-frequency cutoff (default 0.15, strict &lt;).</param>
    /// <param name="windowSize">Cluster window width in codons (default 7, per Sherlocc).</param>
    /// <param name="minRareCodons">Minimum rare codons in a window to call a cluster (default 4, per Sherlocc).</param>
    /// <returns>Maximal, non-overlapping <see cref="RareCodonCluster"/> regions in codon-index order.</returns>
    /// <exception cref="ArgumentOutOfRangeException">windowSize &lt; 1 or minRareCodons &lt; 1.</exception>
    public static IReadOnlyList<RareCodonCluster> FindRareCodonClusters(
        string codingSequence,
        CodonUsageTable table,
        double rareThreshold = 0.15,
        int windowSize = DefaultClusterWindowCodons,
        int minRareCodons = DefaultClusterMinRareCodons)
    {
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1 codon.");
        if (minRareCodons < 1)
            throw new ArgumentOutOfRangeException(nameof(minRareCodons), minRareCodons, "Minimum rare codons must be at least 1.");

        var clusters = new List<RareCodonCluster>();
        if (string.IsNullOrEmpty(codingSequence))
            return clusters;

        string rna = codingSequence.ToUpperInvariant().Replace('T', 'U');
        var codons = SplitIntoCodons(rna);
        if (codons.Count < windowSize)
            return clusters;

        // Mark each codon as rare (pause) when its table frequency is strictly below the threshold.
        var isRare = new bool[codons.Count];
        for (int i = 0; i < codons.Count; i++)
            isRare[i] = table.CodonFrequencies.GetValueOrDefault(codons[i], 0) < rareThreshold;

        int? mergedStart = null;
        int mergedEnd = -1;
        int windowRare = 0;
        for (int start = 0; start + windowSize <= codons.Count; start++)
        {
            if (start == 0)
            {
                for (int k = 0; k < windowSize; k++)
                    if (isRare[k]) windowRare++;
            }
            else
            {
                if (isRare[start - 1]) windowRare--;
                if (isRare[start + windowSize - 1]) windowRare++;
            }

            if (windowRare >= minRareCodons)
            {
                int end = start + windowSize - 1;
                if (mergedStart is null)
                {
                    mergedStart = start;
                    mergedEnd = end;
                }
                else if (start <= mergedEnd + 1)
                {
                    // Overlapping or adjacent qualifying window — extend the current cluster.
                    mergedEnd = end;
                }
                else
                {
                    clusters.Add(BuildCluster(mergedStart.Value, mergedEnd, isRare));
                    mergedStart = start;
                    mergedEnd = end;
                }
            }
        }

        if (mergedStart is not null)
            clusters.Add(BuildCluster(mergedStart.Value, mergedEnd, isRare));

        return clusters;
    }

    private static RareCodonCluster BuildCluster(int start, int end, bool[] isRare)
    {
        int count = 0;
        for (int k = start; k <= end; k++)
            if (isRare[k]) count++;
        return new RareCodonCluster(start, end, count);
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
