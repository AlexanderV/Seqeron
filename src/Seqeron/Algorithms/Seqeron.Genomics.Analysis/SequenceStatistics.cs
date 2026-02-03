using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Provides comprehensive statistical analysis of biological sequences.
/// Includes composition, thermodynamics, and physical properties.
/// </summary>
public static class SequenceStatistics
{
    #region Basic Composition

    /// <summary>
    /// DNA/RNA nucleotide composition statistics.
    /// </summary>
    public readonly record struct NucleotideComposition(
        int Length,
        int CountA,
        int CountT,
        int CountG,
        int CountC,
        int CountU,
        int CountN,
        int CountOther,
        double GcContent,
        double AtContent,
        double GcSkew,
        double AtSkew);

    /// <summary>
    /// Amino acid composition statistics.
    /// </summary>
    public readonly record struct AminoAcidComposition(
        int Length,
        IReadOnlyDictionary<char, int> Counts,
        double MolecularWeight,
        double IsoelectricPoint,
        double Hydrophobicity,
        double ChargedResidueRatio,
        double AromaticResidueRatio);

    /// <summary>
    /// Calculates nucleotide composition of a DNA/RNA sequence.
    /// </summary>
    public static NucleotideComposition CalculateNucleotideComposition(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new NucleotideComposition(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }

        int a = 0, t = 0, g = 0, c = 0, u = 0, n = 0, other = 0;

        foreach (char ch in sequence.ToUpperInvariant())
        {
            switch (ch)
            {
                case 'A': a++; break;
                case 'T': t++; break;
                case 'G': g++; break;
                case 'C': c++; break;
                case 'U': u++; break;
                case 'N': n++; break;
                default: other++; break;
            }
        }

        int total = a + t + g + c + u;
        int gc = g + c;
        int at = a + t + u;

        double gcContent = total > 0 ? (double)gc / total : 0;
        double atContent = total > 0 ? (double)at / total : 0;
        double gcSkew = (g + c) > 0 ? (double)(g - c) / (g + c) : 0;
        double atSkew = (a + t) > 0 ? (double)(a - t) / (a + t) : 0;

        return new NucleotideComposition(
            Length: sequence.Length,
            CountA: a,
            CountT: t,
            CountG: g,
            CountC: c,
            CountU: u,
            CountN: n,
            CountOther: other,
            GcContent: gcContent,
            AtContent: atContent,
            GcSkew: gcSkew,
            AtSkew: atSkew);
    }

    /// <summary>
    /// Calculates amino acid composition of a protein sequence.
    /// </summary>
    public static AminoAcidComposition CalculateAminoAcidComposition(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new AminoAcidComposition(0, new Dictionary<char, int>(), 0, 7.0, 0, 0, 0);
        }

        var counts = new Dictionary<char, int>();
        foreach (char ch in sequence.ToUpperInvariant())
        {
            if (char.IsLetter(ch))
            {
                counts[ch] = counts.GetValueOrDefault(ch) + 1;
            }
        }

        int length = counts.Values.Sum();
        double mw = CalculateMolecularWeight(sequence);
        double pi = CalculateIsoelectricPoint(sequence);
        double hydro = CalculateHydrophobicity(sequence);

        // Charged residues: D, E (negative), K, R, H (positive)
        int charged = counts.GetValueOrDefault('D') + counts.GetValueOrDefault('E') +
                     counts.GetValueOrDefault('K') + counts.GetValueOrDefault('R') +
                     counts.GetValueOrDefault('H');
        double chargedRatio = length > 0 ? (double)charged / length : 0;

        // Aromatic residues: F, Y, W
        int aromatic = counts.GetValueOrDefault('F') + counts.GetValueOrDefault('Y') +
                      counts.GetValueOrDefault('W');
        double aromaticRatio = length > 0 ? (double)aromatic / length : 0;

        return new AminoAcidComposition(
            Length: length,
            Counts: counts,
            MolecularWeight: mw,
            IsoelectricPoint: pi,
            Hydrophobicity: hydro,
            ChargedResidueRatio: chargedRatio,
            AromaticResidueRatio: aromaticRatio);
    }

    #endregion

    #region Molecular Weight

    // Amino acid molecular weights (Da) - average isotopic mass
    private static readonly Dictionary<char, double> AminoAcidWeights = new()
    {
        { 'A', 89.09 },   { 'R', 174.20 }, { 'N', 132.12 }, { 'D', 133.10 },
        { 'C', 121.16 },  { 'E', 147.13 }, { 'Q', 146.15 }, { 'G', 75.07 },
        { 'H', 155.16 },  { 'I', 131.18 }, { 'L', 131.18 }, { 'K', 146.19 },
        { 'M', 149.21 },  { 'F', 165.19 }, { 'P', 115.13 }, { 'S', 105.09 },
        { 'T', 119.12 },  { 'W', 204.23 }, { 'Y', 181.19 }, { 'V', 117.15 },
        { 'U', 168.05 },  { 'O', 255.31 }, { 'B', 132.61 }, { 'Z', 146.64 },
        { 'X', 110.0 }    // Average for unknown
    };

    /// <summary>
    /// Calculates molecular weight of a protein sequence (Da).
    /// </summary>
    public static double CalculateMolecularWeight(string proteinSequence)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            return 0;

        double weight = 18.015; // Water molecule (for peptide bond formation)

        foreach (char aa in proteinSequence.ToUpperInvariant())
        {
            if (AminoAcidWeights.TryGetValue(aa, out double aaWeight))
            {
                weight += aaWeight - 18.015; // Subtract water for each peptide bond
            }
        }

        return weight + 18.015; // Add back water for N and C terminus
    }

    /// <summary>
    /// Calculates molecular weight of a DNA/RNA sequence.
    /// </summary>
    public static double CalculateNucleotideMolecularWeight(string sequence, bool isDna = true)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        // Average molecular weights for nucleotides (as monophosphates)
        double aWeight = isDna ? 331.2 : 347.2;
        double tWeight = 322.2;  // DNA only
        double uWeight = 324.2;  // RNA only
        double gWeight = isDna ? 347.2 : 363.2;
        double cWeight = isDna ? 307.2 : 323.2;

        double weight = 0;
        foreach (char ch in sequence.ToUpperInvariant())
        {
            weight += ch switch
            {
                'A' => aWeight,
                'T' => tWeight,
                'U' => uWeight,
                'G' => gWeight,
                'C' => cWeight,
                _ => 330.0 // Average
            };
        }

        return weight;
    }

    #endregion

    #region Isoelectric Point

    // pKa values for ionizable groups
    private static readonly Dictionary<char, (double pKa, int charge)> IonizableGroups = new()
    {
        { 'D', (3.9, -1) },  // Aspartic acid
        { 'E', (4.1, -1) },  // Glutamic acid
        { 'C', (8.3, -1) },  // Cysteine
        { 'Y', (10.1, -1) }, // Tyrosine
        { 'H', (6.0, 1) },   // Histidine
        { 'K', (10.5, 1) },  // Lysine
        { 'R', (12.5, 1) }   // Arginine
    };

    /// <summary>
    /// Calculates the isoelectric point (pI) of a protein.
    /// </summary>
    public static double CalculateIsoelectricPoint(string proteinSequence)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            return 7.0;

        // Count ionizable residues
        var counts = new Dictionary<char, int>();
        foreach (char aa in proteinSequence.ToUpperInvariant())
        {
            if (IonizableGroups.ContainsKey(aa))
                counts[aa] = counts.GetValueOrDefault(aa) + 1;
        }

        // Binary search for pH where net charge = 0
        double pHLow = 0.0;
        double pHHigh = 14.0;
        double pH = 7.0;
        const double precision = 0.01;

        // N-terminus pKa = 9.6, C-terminus pKa = 2.3
        const double nTermPka = 9.6;
        const double cTermPka = 2.3;

        while (pHHigh - pHLow > precision)
        {
            pH = (pHLow + pHHigh) / 2.0;

            // Calculate net charge at this pH
            double charge = 0;

            // N-terminus (positive)
            charge += 1.0 / (1.0 + Math.Pow(10, pH - nTermPka));

            // C-terminus (negative)
            charge -= 1.0 / (1.0 + Math.Pow(10, cTermPka - pH));

            // Side chains
            foreach (var (aa, count) in counts)
            {
                var (pKa, baseCharge) = IonizableGroups[aa];
                if (baseCharge > 0)
                {
                    // Positive residue
                    charge += count * 1.0 / (1.0 + Math.Pow(10, pH - pKa));
                }
                else
                {
                    // Negative residue
                    charge -= count * 1.0 / (1.0 + Math.Pow(10, pKa - pH));
                }
            }

            if (charge > 0)
                pHLow = pH;
            else
                pHHigh = pH;
        }

        return Math.Round(pH, 2);
    }

    #endregion

    #region Hydrophobicity

    // Kyte-Doolittle hydrophobicity scale
    private static readonly Dictionary<char, double> HydrophobicityScale = new()
    {
        { 'A', 1.8 },  { 'R', -4.5 }, { 'N', -3.5 }, { 'D', -3.5 },
        { 'C', 2.5 },  { 'E', -3.5 }, { 'Q', -3.5 }, { 'G', -0.4 },
        { 'H', -3.2 }, { 'I', 4.5 },  { 'L', 3.8 },  { 'K', -3.9 },
        { 'M', 1.9 },  { 'F', 2.8 },  { 'P', -1.6 }, { 'S', -0.8 },
        { 'T', -0.7 }, { 'W', -0.9 }, { 'Y', -1.3 }, { 'V', 4.2 }
    };

    /// <summary>
    /// Calculates the grand average of hydropathy (GRAVY) index.
    /// </summary>
    public static double CalculateHydrophobicity(string proteinSequence)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            return 0;

        double sum = 0;
        int count = 0;

        foreach (char aa in proteinSequence.ToUpperInvariant())
        {
            if (HydrophobicityScale.TryGetValue(aa, out double value))
            {
                sum += value;
                count++;
            }
        }

        return count > 0 ? sum / count : 0;
    }

    /// <summary>
    /// Calculates hydrophobicity profile using a sliding window.
    /// </summary>
    public static IEnumerable<double> CalculateHydrophobicityProfile(
        string proteinSequence,
        int windowSize = 9)
    {
        if (string.IsNullOrEmpty(proteinSequence) || windowSize > proteinSequence.Length)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        for (int i = 0; i <= upper.Length - windowSize; i++)
        {
            double sum = 0;
            for (int j = 0; j < windowSize; j++)
            {
                if (HydrophobicityScale.TryGetValue(upper[i + j], out double value))
                    sum += value;
            }
            yield return sum / windowSize;
        }
    }

    #endregion

    #region DNA Thermodynamics

    // Nearest-neighbor thermodynamic parameters (kcal/mol for ΔH and ΔS)
    private static readonly Dictionary<string, (double dH, double dS)> NearestNeighborParams = new()
    {
        { "AA", (-7.9, -22.2) }, { "TT", (-7.9, -22.2) },
        { "AT", (-7.2, -20.4) },
        { "TA", (-7.2, -21.3) },
        { "CA", (-8.5, -22.7) }, { "TG", (-8.5, -22.7) },
        { "GT", (-8.4, -22.4) }, { "AC", (-8.4, -22.4) },
        { "CT", (-7.8, -21.0) }, { "AG", (-7.8, -21.0) },
        { "GA", (-8.2, -22.2) }, { "TC", (-8.2, -22.2) },
        { "CG", (-10.6, -27.2) },
        { "GC", (-9.8, -24.4) },
        { "GG", (-8.0, -19.9) }, { "CC", (-8.0, -19.9) }
    };

    /// <summary>
    /// DNA thermodynamic properties.
    /// </summary>
    public readonly record struct ThermodynamicProperties(
        double DeltaH,
        double DeltaS,
        double DeltaG,
        double MeltingTemperature);

    /// <summary>
    /// Calculates thermodynamic properties of a DNA duplex.
    /// </summary>
    public static ThermodynamicProperties CalculateThermodynamics(
        string dnaSequence,
        double naConcentration = 0.05, // 50 mM
        double primerConcentration = 0.00000025) // 250 nM
    {
        if (string.IsNullOrEmpty(dnaSequence) || dnaSequence.Length < 2)
            return new ThermodynamicProperties(0, 0, 0, 0);

        string upper = dnaSequence.ToUpperInvariant();

        // Calculate ΔH and ΔS using nearest-neighbor method
        double dH = 0;
        double dS = 0;

        // Initiation parameters
        if (upper[0] == 'G' || upper[0] == 'C')
        {
            dH += 0.1;
            dS += -2.8;
        }
        else
        {
            dH += 2.3;
            dS += 4.1;
        }

        // Sum nearest-neighbor contributions
        for (int i = 0; i < upper.Length - 1; i++)
        {
            string dinuc = upper.Substring(i, 2);
            if (NearestNeighborParams.TryGetValue(dinuc, out var param))
            {
                dH += param.dH;
                dS += param.dS;
            }
        }

        // Salt correction for ΔS
        double saltCorrection = 0.368 * (upper.Length - 1) * Math.Log(naConcentration);
        dS += saltCorrection;

        // Calculate ΔG at 37°C (310.15 K)
        double dG = dH - (310.15 * dS / 1000.0);

        // Calculate Tm
        // Tm = ΔH / (ΔS + R * ln(Ct/4))
        // R = 1.987 cal/(mol·K)
        double R = 1.987;
        double tm = (dH * 1000) / (dS + R * Math.Log(primerConcentration / 4.0)) - 273.15;

        return new ThermodynamicProperties(
            DeltaH: Math.Round(dH, 2),
            DeltaS: Math.Round(dS, 2),
            DeltaG: Math.Round(dG, 2),
            MeltingTemperature: Math.Round(tm, 1));
    }

    /// <summary>
    /// Calculates simple melting temperature using Wallace rule or GC formula.
    /// </summary>
    public static double CalculateMeltingTemperature(string dnaSequence, bool useWallaceRule = true)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            return 0;

        var comp = CalculateNucleotideComposition(dnaSequence);

        if (useWallaceRule && dnaSequence.Length < ThermoConstants.WallaceMaxLength)
        {
            // Wallace rule for short oligos: Tm = 2(A+T) + 4(G+C)
            return ThermoConstants.CalculateWallaceTm(
                comp.CountA + comp.CountT,
                comp.CountG + comp.CountC);
        }
        else
        {
            // GC formula (Marmur-Doty)
            int total = comp.CountA + comp.CountT + comp.CountG + comp.CountC;
            if (total == 0) return 0;
            return ThermoConstants.CalculateMarmurDotyTm(comp.CountG + comp.CountC, total);
        }
    }

    #endregion

    #region Sequence Patterns

    /// <summary>
    /// Calculates dinucleotide frequencies.
    /// </summary>
    public static IReadOnlyDictionary<string, double> CalculateDinucleotideFrequencies(string sequence)
    {
        var counts = new Dictionary<string, int>();
        var freq = new Dictionary<string, double>();

        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return freq;

        string upper = sequence.ToUpperInvariant();

        // Count all dinucleotides
        int total = 0;
        for (int i = 0; i < upper.Length - 1; i++)
        {
            string dinuc = upper.Substring(i, 2);
            if (dinuc.All(c => "ATGCU".Contains(c)))
            {
                counts[dinuc] = counts.GetValueOrDefault(dinuc) + 1;
                total++;
            }
        }

        // Convert to frequencies
        foreach (var (dinuc, count) in counts)
        {
            freq[dinuc] = (double)count / total;
        }

        return freq;
    }

    /// <summary>
    /// Calculates observed vs expected dinucleotide ratios (CpG ratio, etc.).
    /// </summary>
    public static IReadOnlyDictionary<string, double> CalculateDinucleotideRatios(string sequence)
    {
        var ratios = new Dictionary<string, double>();

        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return ratios;

        var comp = CalculateNucleotideComposition(sequence);
        var dinucFreq = CalculateDinucleotideFrequencies(sequence);

        int total = comp.CountA + comp.CountT + comp.CountG + comp.CountC + comp.CountU;
        if (total == 0) return ratios;

        // Single nucleotide frequencies
        var singleFreq = new Dictionary<char, double>
        {
            { 'A', (double)comp.CountA / total },
            { 'T', (double)comp.CountT / total },
            { 'G', (double)comp.CountG / total },
            { 'C', (double)comp.CountC / total },
            { 'U', (double)comp.CountU / total }
        };

        // Calculate observed/expected ratios
        foreach (var (dinuc, observed) in dinucFreq)
        {
            double expected = singleFreq.GetValueOrDefault(dinuc[0]) *
                             singleFreq.GetValueOrDefault(dinuc[1]);
            ratios[dinuc] = expected > 0 ? observed / expected : 0;
        }

        return ratios;
    }

    /// <summary>
    /// Calculates codon usage frequencies.
    /// </summary>
    public static IReadOnlyDictionary<string, double> CalculateCodonFrequencies(
        string dnaSequence,
        int readingFrame = 0)
    {
        var counts = new Dictionary<string, int>();
        var freq = new Dictionary<string, double>();

        if (string.IsNullOrEmpty(dnaSequence) || dnaSequence.Length < 3)
            return freq;

        string upper = dnaSequence.ToUpperInvariant();
        int total = 0;

        for (int i = readingFrame; i <= upper.Length - 3; i += 3)
        {
            string codon = upper.Substring(i, 3);
            if (codon.All(c => "ATGC".Contains(c)))
            {
                counts[codon] = counts.GetValueOrDefault(codon) + 1;
                total++;
            }
        }

        foreach (var (codon, count) in counts)
        {
            freq[codon] = (double)count / total;
        }

        return freq;
    }

    #endregion

    #region Entropy and Complexity

    /// <summary>
    /// Calculates Shannon entropy of a sequence.
    /// </summary>
    public static double CalculateShannonEntropy(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        var counts = new Dictionary<char, int>();
        int total = 0;

        foreach (char ch in sequence.ToUpperInvariant())
        {
            if (char.IsLetter(ch))
            {
                counts[ch] = counts.GetValueOrDefault(ch) + 1;
                total++;
            }
        }

        if (total == 0) return 0;

        double entropy = 0;
        foreach (int count in counts.Values)
        {
            double freq = (double)count / total;
            if (freq > 0)
            {
                entropy -= freq * Math.Log2(freq);
            }
        }

        return entropy;
    }

    /// <summary>
    /// Calculates linguistic complexity of a sequence.
    /// </summary>
    public static double CalculateLinguisticComplexity(string sequence, int maxK = 6)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        string upper = sequence.ToUpperInvariant();
        int n = upper.Length;
        double totalRatio = 0;
        int kCount = 0;

        for (int k = 1; k <= Math.Min(maxK, n); k++)
        {
            var observedKmers = new HashSet<string>();
            for (int i = 0; i <= n - k; i++)
            {
                observedKmers.Add(upper.Substring(i, k));
            }

            // Maximum possible k-mers
            int maxPossible = Math.Min((int)Math.Pow(4, k), n - k + 1);
            if (maxPossible > 0)
            {
                totalRatio += (double)observedKmers.Count / maxPossible;
                kCount++;
            }
        }

        return kCount > 0 ? totalRatio / kCount : 0;
    }

    #endregion

    #region Protein Secondary Structure

    // Chou-Fasman propensity parameters
    private static readonly Dictionary<char, (double helix, double sheet, double turn)> SecondaryStructurePropensity = new()
    {
        { 'A', (1.42, 0.83, 0.66) }, { 'R', (0.98, 0.93, 0.95) },
        { 'N', (0.67, 0.89, 1.56) }, { 'D', (1.01, 0.54, 1.46) },
        { 'C', (0.70, 1.19, 1.19) }, { 'E', (1.51, 0.37, 0.74) },
        { 'Q', (1.11, 1.10, 0.98) }, { 'G', (0.57, 0.75, 1.56) },
        { 'H', (1.00, 0.87, 0.95) }, { 'I', (1.08, 1.60, 0.47) },
        { 'L', (1.21, 1.30, 0.59) }, { 'K', (1.16, 0.74, 1.01) },
        { 'M', (1.45, 1.05, 0.60) }, { 'F', (1.13, 1.38, 0.60) },
        { 'P', (0.57, 0.55, 1.52) }, { 'S', (0.77, 0.75, 1.43) },
        { 'T', (0.83, 1.19, 0.96) }, { 'W', (1.08, 1.37, 0.96) },
        { 'Y', (0.69, 1.47, 1.14) }, { 'V', (1.06, 1.70, 0.50) }
    };

    /// <summary>
    /// Predicts secondary structure propensities using Chou-Fasman parameters.
    /// </summary>
    public static IEnumerable<(double Helix, double Sheet, double Turn)> PredictSecondaryStructure(
        string proteinSequence,
        int windowSize = 7)
    {
        if (string.IsNullOrEmpty(proteinSequence) || windowSize > proteinSequence.Length)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        for (int i = 0; i <= upper.Length - windowSize; i++)
        {
            double helixSum = 0, sheetSum = 0, turnSum = 0;
            int count = 0;

            for (int j = 0; j < windowSize; j++)
            {
                if (SecondaryStructurePropensity.TryGetValue(upper[i + j], out var prop))
                {
                    helixSum += prop.helix;
                    sheetSum += prop.sheet;
                    turnSum += prop.turn;
                    count++;
                }
            }

            if (count > 0)
            {
                yield return (helixSum / count, sheetSum / count, turnSum / count);
            }
        }
    }

    #endregion

    #region Sequence Windows

    /// <summary>
    /// Calculates GC content in sliding windows.
    /// </summary>
    public static IEnumerable<double> CalculateGcContentProfile(
        string sequence,
        int windowSize = 100,
        int stepSize = 1)
    {
        if (string.IsNullOrEmpty(sequence) || windowSize > sequence.Length)
            yield break;

        string upper = sequence.ToUpperInvariant();

        for (int i = 0; i <= upper.Length - windowSize; i += stepSize)
        {
            int gc = 0;
            int total = 0;

            for (int j = 0; j < windowSize; j++)
            {
                char ch = upper[i + j];
                if (ch == 'G' || ch == 'C')
                {
                    gc++;
                    total++;
                }
                else if (ch == 'A' || ch == 'T' || ch == 'U')
                {
                    total++;
                }
            }

            yield return total > 0 ? (double)gc / total : 0;
        }
    }

    /// <summary>
    /// Calculates entropy in sliding windows.
    /// </summary>
    public static IEnumerable<double> CalculateEntropyProfile(
        string sequence,
        int windowSize = 50,
        int stepSize = 1)
    {
        if (string.IsNullOrEmpty(sequence) || windowSize > sequence.Length)
            yield break;

        for (int i = 0; i <= sequence.Length - windowSize; i += stepSize)
        {
            string window = sequence.Substring(i, windowSize);
            yield return CalculateShannonEntropy(window);
        }
    }

    #endregion

    #region Summary Statistics

    /// <summary>
    /// Comprehensive sequence statistics.
    /// </summary>
    public readonly record struct SequenceSummary(
        int Length,
        double GcContent,
        double Entropy,
        double Complexity,
        double MeltingTemperature,
        IReadOnlyDictionary<char, int> Composition);

    /// <summary>
    /// Generates comprehensive summary statistics for a DNA/RNA sequence.
    /// </summary>
    public static SequenceSummary SummarizeNucleotideSequence(string sequence)
    {
        var comp = CalculateNucleotideComposition(sequence);
        double entropy = CalculateShannonEntropy(sequence);
        double complexity = CalculateLinguisticComplexity(sequence);
        double tm = CalculateMeltingTemperature(sequence, useWallaceRule: sequence.Length < 14);

        var composition = new Dictionary<char, int>
        {
            { 'A', comp.CountA },
            { 'T', comp.CountT },
            { 'G', comp.CountG },
            { 'C', comp.CountC },
            { 'U', comp.CountU },
            { 'N', comp.CountN }
        };

        return new SequenceSummary(
            Length: comp.Length,
            GcContent: comp.GcContent,
            Entropy: entropy,
            Complexity: complexity,
            MeltingTemperature: tm,
            Composition: composition);
    }

    #endregion
}
