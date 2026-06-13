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
    /// For codon j of an amino acid with n synonymous codons and observed counts x,
    /// RSCU = x_j / ((1/n) * sum_k x_k) = (n * x_j) / sum_k x_k.
    /// RSCU = 1 means no bias, &gt; 1 means over-represented, &lt; 1 means under-represented.
    /// Definition per Sharp, Tuohy &amp; Mosurski (1986), Nucleic Acids Res. 14(13):5125-5143.
    /// Single-codon families (Met, Trp) always yield RSCU = 1 when present.
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

    // --- Constants per Wright F. (1990) Gene 87(1):23–29, as reproduced verbatim in
    //     Fuglsang A. (2004) BBRC 317:957–964 (Eqs. 1–5a). ---

    // Number of synonymous-codon amino acids in each degeneracy class of the
    // standard (NCBI table 1) genetic code, used as the numerators of Wright Eq. (3).
    private const int TwoFoldAminoAcidCount = 9;   // 9 doublets (His, Gln, …)
    private const int ThreeFoldAminoAcidCount = 1; // 1 triplet  (Ile)
    private const int FourFoldAminoAcidCount = 5;  // 5 quartets (Ala, Gly, Pro, Thr, Val)
    private const int SixFoldAminoAcidCount = 3;   // 3 sextets  (Leu, Ser, Arg)

    // The two single-codon amino acids Met (ATG) and Trp (TGG) each contribute exactly
    // one effective codon; this is the constant "2" in Wright Eq. (3).
    private const double SingleCodonAminoAcidContribution = 2.0;

    // Wright Eq. (3): if Nc exceeds 61 it is re-adjusted down to 61 (the maximum number
    // of sense codons in the standard genetic code).
    private const double MaxEffectiveCodons = 61.0;

    // Structural lower bound: every degeneracy class collapsed to one codon gives Nc = 20
    // (the extreme-bias limit stated by Wright/Fuglsang).
    private const double MinEffectiveCodons = 20.0;

    /// <summary>
    /// Calculates the Effective Number of Codons (ENC / Nc) per Wright (1990).
    /// Nc ranges from 20 (extreme bias — one codon per amino acid) to 61 (no bias).
    /// </summary>
    public static double CalculateEnc(DnaSequence sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return CalculateEncCore(sequence.Sequence);
    }

    /// <summary>
    /// Calculates ENC from a raw sequence string. Null/empty returns 0.
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

        // Wright Eq. (1): per-amino-acid codon homozygosity, grouped by degeneracy class.
        var fByDegeneracy = new Dictionary<int, List<double>>();

        foreach (var aaGroup in CodonToAminoAcid.GroupBy(kv => kv.Value))
        {
            // Exclude stop codons ('*'); they are not amino acids and not counted in Nc.
            if (aaGroup.Key == '*') continue;

            var synonymousCodons = aaGroup.Select(kv => kv.Key).ToList();
            int degeneracy = synonymousCodons.Count;

            if (degeneracy == 1) continue; // Met / Trp handled by SingleCodonAminoAcidContribution.

            int n = synonymousCodons.Sum(c => counts.GetValueOrDefault(c, 0));
            if (n <= 1) continue; // F̂ undefined for n ≤ 1 (denominator n − 1); Fuglsang 2004.

            // Wright Eq. (1): F̂ = (n·Σ p_i² − 1)/(n − 1), p_i = n_i/n.
            double sumPSquared = 0;
            foreach (var codon in synonymousCodons)
            {
                double p = (double)counts.GetValueOrDefault(codon, 0) / n;
                sumPSquared += p * p;
            }
            double f = (n * sumPSquared - 1) / (n - 1);

            if (!fByDegeneracy.TryGetValue(degeneracy, out var list))
            {
                list = new List<double>();
                fByDegeneracy[degeneracy] = list;
            }
            list.Add(f);
        }

        // Wright Eq. (4): the class average F̂ substitutes for any amino acid that cannot
        // be estimated within the same degeneracy class.
        double? f2 = AverageOrNull(fByDegeneracy, 2);
        double? f3 = AverageOrNull(fByDegeneracy, 3);
        double? f4 = AverageOrNull(fByDegeneracy, 4);
        double? f6 = AverageOrNull(fByDegeneracy, 6);

        // Wright Eq. (5a): when isoleucine (the only 3-fold amino acid) cannot be
        // estimated, F̂₃ = (F̂₂ + F̂₄)/2.
        if (f3 is null && f2 is not null && f4 is not null)
            f3 = (f2.Value + f4.Value) / 2.0;

        // Wright Eq. (3): Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆.
        // A class with no estimable F̂ (and, for Ile, no Eq. 5a fallback) contributes its
        // full codon count, i.e. all its codons are assumed effectively present.
        double enc = SingleCodonAminoAcidContribution
            + ClassContribution(TwoFoldAminoAcidCount, f2)
            + ClassContribution(ThreeFoldAminoAcidCount, f3)
            + ClassContribution(FourFoldAminoAcidCount, f4)
            + ClassContribution(SixFoldAminoAcidCount, f6);

        return Math.Min(MaxEffectiveCodons, Math.Max(MinEffectiveCodons, enc));
    }

    private static double? AverageOrNull(Dictionary<int, List<double>> fByDegeneracy, int degeneracy)
        => fByDegeneracy.TryGetValue(degeneracy, out var list) && list.Count > 0
            ? list.Average()
            : null;

    private static double ClassContribution(int aminoAcidCount, double? averageF)
        // No estimable homozygosity for the class ⇒ assume all codons of every amino acid
        // in the class are effectively in use (contribution equals the codon count).
        => averageF is double f && f > 0 ? aminoAcidCount / f : aminoAcidCount;

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
