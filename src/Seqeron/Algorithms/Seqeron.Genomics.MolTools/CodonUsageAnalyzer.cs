using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Analyzes codon usage patterns in coding sequences.
/// Useful for studying codon bias, gene expression optimization, and evolutionary analysis.
/// </summary>
public static class CodonUsageAnalyzer
{
    #region Codon Usage Tables

    /// <summary>
    /// Counts codon occurrences in a coding sequence.
    /// </summary>
    /// <param name="sequence">Coding DNA sequence (must be multiple of 3).</param>
    /// <returns>Dictionary of codon counts.</returns>
    public static Dictionary<string, int> CountCodons(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CountCodonsCore(sequence.Sequence);
    }

    /// <summary>
    /// Counts codon occurrences in a raw sequence string.
    /// </summary>
    public static Dictionary<string, int> CountCodons(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, int>();

        return CountCodonsCore(sequence.ToUpperInvariant());
    }

    private static Dictionary<string, int> CountCodonsCore(string seq)
    {
        var counts = new Dictionary<string, int>();

        for (int i = 0; i + 3 <= seq.Length; i += 3)
        {
            string codon = seq.Substring(i, 3);
            if (IsValidCodon(codon))
            {
                counts.TryGetValue(codon, out int count);
                counts[codon] = count + 1;
            }
        }

        return counts;
    }

    private static bool IsValidCodon(string codon)
    {
        return codon.Length == 3 && codon.All(c => c is 'A' or 'C' or 'G' or 'T');
    }

    #endregion

    #region RSCU (Relative Synonymous Codon Usage)

    /// <summary>
    /// Calculates Relative Synonymous Codon Usage (RSCU).
    /// RSCU = (observed frequency) / (expected frequency if all synonymous codons are used equally)
    /// RSCU = 1 means no bias, > 1 means over-represented, < 1 means under-represented.
    /// </summary>
    public static Dictionary<string, double> CalculateRscu(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateRscuCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates RSCU from a raw sequence string.
    /// </summary>
    public static Dictionary<string, double> CalculateRscu(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return new Dictionary<string, double>();

        return CalculateRscuCore(sequence.ToUpperInvariant());
    }

    private static Dictionary<string, double> CalculateRscuCore(string seq)
    {
        var counts = CountCodonsCore(seq);
        var rscu = new Dictionary<string, double>();

        // Group codons by amino acid
        foreach (var aaGroup in CodonToAminoAcid.GroupBy(kv => kv.Value))
        {
            var synonymousCodons = aaGroup.Select(kv => kv.Key).ToList();
            int totalCount = synonymousCodons.Sum(c => counts.GetValueOrDefault(c, 0));
            int numSynonymous = synonymousCodons.Count;

            foreach (var codon in synonymousCodons)
            {
                int observed = counts.GetValueOrDefault(codon, 0);
                double expected = (double)totalCount / numSynonymous;

                rscu[codon] = expected > 0 ? observed / expected : 0;
            }
        }

        return rscu;
    }

    #endregion

    #region CAI (Codon Adaptation Index)

    /// <summary>
    /// Calculates Codon Adaptation Index (CAI) using a reference codon table.
    /// CAI measures how well codon usage matches highly expressed genes.
    /// Range: 0-1, where 1 means optimal codon usage.
    /// </summary>
    /// <param name="sequence">Coding sequence to analyze.</param>
    /// <param name="referenceRscu">RSCU values from reference set (e.g., highly expressed genes).</param>
    public static double CalculateCai(DnaSequence sequence, Dictionary<string, double> referenceRscu)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(referenceRscu);

        return CalculateCaiCore(sequence.Sequence, referenceRscu);
    }

    /// <summary>
    /// Calculates CAI from a raw sequence string.
    /// </summary>
    public static double CalculateCai(string sequence, Dictionary<string, double> referenceRscu)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateCaiCore(sequence.ToUpperInvariant(), referenceRscu);
    }

    private static double CalculateCaiCore(string seq, Dictionary<string, double> referenceRscu)
    {
        // Calculate relative adaptiveness (w) for each codon
        var relativeAdaptiveness = new Dictionary<string, double>();

        foreach (var aaGroup in CodonToAminoAcid.GroupBy(kv => kv.Value))
        {
            var synonymousCodons = aaGroup.Select(kv => kv.Key).ToList();
            double maxRscu = synonymousCodons.Max(c => referenceRscu.GetValueOrDefault(c, 0));

            foreach (var codon in synonymousCodons)
            {
                double rscu = referenceRscu.GetValueOrDefault(codon, 0);
                relativeAdaptiveness[codon] = maxRscu > 0 ? rscu / maxRscu : 0;
            }
        }

        // Calculate CAI as geometric mean of relative adaptiveness values
        double logSum = 0;
        int codonCount = 0;

        for (int i = 0; i + 3 <= seq.Length; i += 3)
        {
            string codon = seq.Substring(i, 3);
            if (IsValidCodon(codon) && relativeAdaptiveness.TryGetValue(codon, out double w) && w > 0)
            {
                logSum += Math.Log(w);
                codonCount++;
            }
        }

        return codonCount > 0 ? Math.Exp(logSum / codonCount) : 0;
    }

    /// <summary>
    /// Gets a reference RSCU table for E. coli highly expressed genes.
    /// </summary>
    public static Dictionary<string, double> EColiOptimalCodons => new()
    {
        // Preferred codons in E. coli highly expressed genes
        ["TTT"] = 0.30,
        ["TTC"] = 1.70,
        ["TTA"] = 0.07,
        ["TTG"] = 0.10,
        ["CTT"] = 0.10,
        ["CTC"] = 0.07,
        ["CTA"] = 0.03,
        ["CTG"] = 5.63,
        ["ATT"] = 0.45,
        ["ATC"] = 2.55,
        ["ATA"] = 0.01,
        ["ATG"] = 1.00,
        ["GTT"] = 2.03,
        ["GTC"] = 0.15,
        ["GTA"] = 0.83,
        ["GTG"] = 0.99,
        ["TCT"] = 1.65,
        ["TCC"] = 1.32,
        ["TCA"] = 0.10,
        ["TCG"] = 0.10,
        ["CCT"] = 0.12,
        ["CCC"] = 0.07,
        ["CCA"] = 0.32,
        ["CCG"] = 3.49,
        ["ACT"] = 1.43,
        ["ACC"] = 2.29,
        ["ACA"] = 0.15,
        ["ACG"] = 0.12,
        ["GCT"] = 1.51,
        ["GCC"] = 0.30,
        ["GCA"] = 1.07,
        ["GCG"] = 1.12,
        ["TAT"] = 0.44,
        ["TAC"] = 1.56,
        ["TAA"] = 1.64,
        ["TAG"] = 0.00,
        ["CAT"] = 0.40,
        ["CAC"] = 1.60,
        ["CAA"] = 0.14,
        ["CAG"] = 1.86,
        ["AAT"] = 0.10,
        ["AAC"] = 1.90,
        ["AAA"] = 1.52,
        ["AAG"] = 0.48,
        ["GAT"] = 0.56,
        ["GAC"] = 1.44,
        ["GAA"] = 1.48,
        ["GAG"] = 0.52,
        ["TGT"] = 0.50,
        ["TGC"] = 1.50,
        ["TGA"] = 1.36,
        ["TGG"] = 1.00,
        ["CGT"] = 4.88,
        ["CGC"] = 0.88,
        ["CGA"] = 0.04,
        ["CGG"] = 0.04,
        ["AGT"] = 0.13,
        ["AGC"] = 2.70,
        ["AGA"] = 0.04,
        ["AGG"] = 0.04,
        ["GGT"] = 2.76,
        ["GGC"] = 1.20,
        ["GGA"] = 0.04,
        ["GGG"] = 0.04
    };

    /// <summary>
    /// Gets a reference RSCU table for human highly expressed genes.
    /// </summary>
    public static Dictionary<string, double> HumanOptimalCodons => new()
    {
        ["TTT"] = 0.87,
        ["TTC"] = 1.13,
        ["TTA"] = 0.43,
        ["TTG"] = 0.77,
        ["CTT"] = 0.78,
        ["CTC"] = 1.17,
        ["CTA"] = 0.43,
        ["CTG"] = 2.41,
        ["ATT"] = 1.08,
        ["ATC"] = 1.41,
        ["ATA"] = 0.51,
        ["ATG"] = 1.00,
        ["GTT"] = 0.72,
        ["GTC"] = 0.95,
        ["GTA"] = 0.47,
        ["GTG"] = 1.86,
        ["TCT"] = 1.14,
        ["TCC"] = 1.32,
        ["TCA"] = 0.90,
        ["TCG"] = 0.33,
        ["CCT"] = 1.16,
        ["CCC"] = 1.29,
        ["CCA"] = 1.09,
        ["CCG"] = 0.45,
        ["ACT"] = 0.99,
        ["ACC"] = 1.41,
        ["ACA"] = 1.14,
        ["ACG"] = 0.45,
        ["GCT"] = 1.08,
        ["GCC"] = 1.60,
        ["GCA"] = 0.90,
        ["GCG"] = 0.42,
        ["TAT"] = 0.88,
        ["TAC"] = 1.12,
        ["TAA"] = 1.00,
        ["TAG"] = 0.80,
        ["CAT"] = 0.84,
        ["CAC"] = 1.16,
        ["CAA"] = 0.54,
        ["CAG"] = 1.46,
        ["AAT"] = 0.94,
        ["AAC"] = 1.06,
        ["AAA"] = 0.86,
        ["AAG"] = 1.14,
        ["GAT"] = 0.92,
        ["GAC"] = 1.08,
        ["GAA"] = 0.84,
        ["GAG"] = 1.16,
        ["TGT"] = 0.92,
        ["TGC"] = 1.08,
        ["TGA"] = 1.20,
        ["TGG"] = 1.00,
        ["CGT"] = 0.48,
        ["CGC"] = 1.08,
        ["CGA"] = 0.66,
        ["CGG"] = 1.20,
        ["AGT"] = 0.90,
        ["AGC"] = 1.41,
        ["AGA"] = 1.26,
        ["AGG"] = 1.32,
        ["GGT"] = 0.64,
        ["GGC"] = 1.36,
        ["GGA"] = 1.00,
        ["GGG"] = 1.00
    };

    #endregion

    #region Effective Number of Codons (ENC)

    /// <summary>
    /// Calculates Effective Number of Codons (ENC/Nc).
    /// ENC ranges from 20 (extreme bias - one codon per amino acid) to 61 (no bias).
    /// </summary>
    public static double CalculateEnc(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateEncCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates ENC from a raw sequence string.
    /// </summary>
    public static double CalculateEnc(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        return CalculateEncCore(sequence.ToUpperInvariant());
    }

    private static double CalculateEncCore(string seq)
    {
        var counts = CountCodonsCore(seq);

        // Group by amino acid degeneracy
        var fValues = new Dictionary<int, List<double>>();

        foreach (var aaGroup in CodonToAminoAcid.GroupBy(kv => kv.Value))
        {
            var synonymousCodons = aaGroup.Select(kv => kv.Key).ToList();
            int degeneracy = synonymousCodons.Count;

            if (degeneracy == 1) continue; // Skip Met and Trp

            var codonCounts = synonymousCodons.Select(c => counts.GetValueOrDefault(c, 0)).ToList();
            int n = codonCounts.Sum();

            if (n <= 1) continue;

            // Calculate F for this amino acid
            double sumPiSquared = codonCounts.Sum(ni => (double)ni * ni);
            double f = (n * sumPiSquared - 1) / (n - 1) / n;

            if (!fValues.ContainsKey(degeneracy))
                fValues[degeneracy] = new List<double>();

            fValues[degeneracy].Add(f);
        }

        // Calculate average F for each degeneracy class
        double enc = 0;

        // 2-fold degenerate (9 amino acids)
        if (fValues.ContainsKey(2) && fValues[2].Count > 0)
            enc += 9 / fValues[2].Average();
        else
            enc += 9;

        // 3-fold degenerate (1 amino acid - Ile)
        if (fValues.ContainsKey(3) && fValues[3].Count > 0)
            enc += 1 / fValues[3].Average();
        else
            enc += 1;

        // 4-fold degenerate (5 amino acids)
        if (fValues.ContainsKey(4) && fValues[4].Count > 0)
            enc += 5 / fValues[4].Average();
        else
            enc += 5;

        // 6-fold degenerate (3 amino acids)
        if (fValues.ContainsKey(6) && fValues[6].Count > 0)
            enc += 3 / fValues[6].Average();
        else
            enc += 3;

        // Add Met (1) + Trp (1) + Stop (3) = 5 single codons are fixed
        enc += 2;

        return Math.Min(61, Math.Max(20, enc));
    }

    #endregion

    #region Codon Usage Statistics

    /// <summary>
    /// Gets comprehensive codon usage statistics.
    /// </summary>
    public static CodonUsageStatistics GetStatistics(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return GetStatisticsCore(sequence.Sequence);
    }

    /// <summary>
    /// Gets codon usage statistics from a raw sequence string.
    /// </summary>
    public static CodonUsageStatistics GetStatistics(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return new CodonUsageStatistics(
                new Dictionary<string, int>(),
                new Dictionary<string, double>(),
                0, 0, 0, 0, 0, 0);

        return GetStatisticsCore(sequence.ToUpperInvariant());
    }

    private static CodonUsageStatistics GetStatisticsCore(string seq)
    {
        var counts = CountCodonsCore(seq);
        var rscu = CalculateRscuCore(seq);
        double enc = CalculateEncCore(seq);

        int totalCodons = counts.Values.Sum();

        // GC content at different codon positions
        int gc1 = 0, gc2 = 0, gc3 = 0;
        int positionCount = 0;

        for (int i = 0; i + 3 <= seq.Length; i += 3)
        {
            string codon = seq.Substring(i, 3);
            if (IsValidCodon(codon))
            {
                gc1 += IsGC(codon[0]) ? 1 : 0;
                gc2 += IsGC(codon[1]) ? 1 : 0;
                gc3 += IsGC(codon[2]) ? 1 : 0;
                positionCount++;
            }
        }

        double gc1Percent = positionCount > 0 ? (double)gc1 / positionCount * 100 : 0;
        double gc2Percent = positionCount > 0 ? (double)gc2 / positionCount * 100 : 0;
        double gc3Percent = positionCount > 0 ? (double)gc3 / positionCount * 100 : 0;
        double gc3s = gc3Percent; // GC3s (synonymous third position GC)

        return new CodonUsageStatistics(
            CodonCounts: counts,
            Rscu: rscu,
            Enc: enc,
            TotalCodons: totalCodons,
            Gc1: gc1Percent,
            Gc2: gc2Percent,
            Gc3: gc3Percent,
            Gc3s: gc3s);
    }

    private static bool IsGC(char c) => c is 'G' or 'C';

    #endregion

    #region Codon Table

    private static readonly Dictionary<string, char> CodonToAminoAcid = new()
    {
        ["TTT"] = 'F',
        ["TTC"] = 'F',
        ["TTA"] = 'L',
        ["TTG"] = 'L',
        ["CTT"] = 'L',
        ["CTC"] = 'L',
        ["CTA"] = 'L',
        ["CTG"] = 'L',
        ["ATT"] = 'I',
        ["ATC"] = 'I',
        ["ATA"] = 'I',
        ["ATG"] = 'M',
        ["GTT"] = 'V',
        ["GTC"] = 'V',
        ["GTA"] = 'V',
        ["GTG"] = 'V',
        ["TCT"] = 'S',
        ["TCC"] = 'S',
        ["TCA"] = 'S',
        ["TCG"] = 'S',
        ["CCT"] = 'P',
        ["CCC"] = 'P',
        ["CCA"] = 'P',
        ["CCG"] = 'P',
        ["ACT"] = 'T',
        ["ACC"] = 'T',
        ["ACA"] = 'T',
        ["ACG"] = 'T',
        ["GCT"] = 'A',
        ["GCC"] = 'A',
        ["GCA"] = 'A',
        ["GCG"] = 'A',
        ["TAT"] = 'Y',
        ["TAC"] = 'Y',
        ["TAA"] = '*',
        ["TAG"] = '*',
        ["CAT"] = 'H',
        ["CAC"] = 'H',
        ["CAA"] = 'Q',
        ["CAG"] = 'Q',
        ["AAT"] = 'N',
        ["AAC"] = 'N',
        ["AAA"] = 'K',
        ["AAG"] = 'K',
        ["GAT"] = 'D',
        ["GAC"] = 'D',
        ["GAA"] = 'E',
        ["GAG"] = 'E',
        ["TGT"] = 'C',
        ["TGC"] = 'C',
        ["TGA"] = '*',
        ["TGG"] = 'W',
        ["CGT"] = 'R',
        ["CGC"] = 'R',
        ["CGA"] = 'R',
        ["CGG"] = 'R',
        ["AGT"] = 'S',
        ["AGC"] = 'S',
        ["AGA"] = 'R',
        ["AGG"] = 'R',
        ["GGT"] = 'G',
        ["GGC"] = 'G',
        ["GGA"] = 'G',
        ["GGG"] = 'G'
    };

    #endregion
}

/// <summary>
/// Comprehensive codon usage statistics.
/// </summary>
public readonly record struct CodonUsageStatistics(
    IReadOnlyDictionary<string, int> CodonCounts,
    IReadOnlyDictionary<string, double> Rscu,
    double Enc,
    int TotalCodons,
    double Gc1,
    double Gc2,
    double Gc3,
    double Gc3s)
{
    /// <summary>
    /// Gets the overall GC content of the coding sequence.
    /// </summary>
    public double OverallGc => (Gc1 + Gc2 + Gc3) / 3;
}
