using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics;

/// <summary>
/// Designs hybridization probes for various applications (FISH, microarray, Northern blot, etc.).
/// </summary>
public static class ProbeDesigner
{
    #region Records

    /// <summary>
    /// Probe design parameters.
    /// </summary>
    public readonly record struct ProbeParameters(
        int MinLength,
        int MaxLength,
        double MinTm,
        double MaxTm,
        double MinGc,
        double MaxGc,
        int MaxHomopolymer,
        bool AvoidSecondaryStructure,
        double MaxSelfComplementarity);

    /// <summary>
    /// Default probe parameters for different applications.
    /// </summary>
    public static class Defaults
    {
        public static ProbeParameters Microarray => new(
            MinLength: 50, MaxLength: 70,
            MinTm: 75, MaxTm: 85,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 5,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.3);

        public static ProbeParameters FISH => new(
            MinLength: 200, MaxLength: 500,
            MinTm: 70, MaxTm: 90,
            MinGc: 0.35, MaxGc: 0.65,
            MaxHomopolymer: 8,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 0.4);

        public static ProbeParameters NorthernBlot => new(
            MinLength: 100, MaxLength: 300,
            MinTm: 65, MaxTm: 80,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 6,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.35);

        public static ProbeParameters qPCR => new(
            MinLength: 20, MaxLength: 30,
            MinTm: 68, MaxTm: 72,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 4,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.25);

        public static ProbeParameters SouthernBlot => new(
            MinLength: 150, MaxLength: 500,
            MinTm: 65, MaxTm: 75,
            MinGc: 0.35, MaxGc: 0.65,
            MaxHomopolymer: 7,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 0.4);
    }

    /// <summary>
    /// Designed probe.
    /// </summary>
    public readonly record struct Probe(
        string Sequence,
        int Start,
        int End,
        double Tm,
        double GcContent,
        double Score,
        ProbeType Type,
        IReadOnlyList<string> Warnings);

    /// <summary>
    /// Probe set for tiling.
    /// </summary>
    public readonly record struct TilingProbeSet(
        IReadOnlyList<Probe> Probes,
        int Coverage,
        double MeanTm,
        double TmRange);

    /// <summary>
    /// Probe validation result.
    /// </summary>
    public readonly record struct ProbeValidation(
        bool IsValid,
        double SpecificityScore,
        int OffTargetHits,
        double SelfComplementarity,
        bool HasSecondaryStructure,
        IReadOnlyList<string> Issues);

    /// <summary>
    /// Probe type.
    /// </summary>
    public enum ProbeType
    {
        Standard,
        Tiling,
        Antisense,
        LNA, // Locked Nucleic Acid
        MolecularBeacon
    }

    #endregion

    #region Probe Design

    /// <summary>
    /// Designs probes for a target sequence.
    /// </summary>
    public static IEnumerable<Probe> DesignProbes(
        string targetSequence,
        ProbeParameters? parameters = null,
        int maxProbes = 10)
    {
        var param = parameters ?? Defaults.Microarray;

        if (string.IsNullOrEmpty(targetSequence) || targetSequence.Length < param.MinLength)
            yield break;

        targetSequence = targetSequence.ToUpperInvariant();
        var candidates = new List<(Probe Probe, double Score)>();

        // Scan for candidate probes
        for (int length = param.MinLength; length <= param.MaxLength && length <= targetSequence.Length; length++)
        {
            for (int start = 0; start <= targetSequence.Length - length; start++)
            {
                string probeSeq = targetSequence.Substring(start, length);
                var probe = EvaluateProbe(probeSeq, start, param);

                if (probe.HasValue)
                {
                    candidates.Add((probe.Value, probe.Value.Score));
                }
            }
        }

        // Return top probes sorted by score
        foreach (var candidate in candidates.OrderByDescending(c => c.Score).Take(maxProbes))
        {
            yield return candidate.Probe;
        }
    }

    /// <summary>
    /// Evaluates a potential probe sequence.
    /// </summary>
    private static Probe? EvaluateProbe(string sequence, int start, ProbeParameters param)
    {
        var warnings = new List<string>();
        double score = 1.0;

        // Check GC content
        double gc = CalculateGcContent(sequence);
        if (gc < param.MinGc || gc > param.MaxGc)
        {
            score -= 0.3;
            warnings.Add($"GC content {gc:P0} outside range");
        }

        // Calculate Tm
        double tm = CalculateTm(sequence);
        if (tm < param.MinTm || tm > param.MaxTm)
        {
            score -= 0.3;
            warnings.Add($"Tm {tm:F1}°C outside range");
        }

        // Check homopolymers
        int maxHomopolymer = GetMaxHomopolymerLength(sequence);
        if (maxHomopolymer > param.MaxHomopolymer)
        {
            score -= 0.2;
            warnings.Add($"Homopolymer run of {maxHomopolymer}");
        }

        // Check self-complementarity
        double selfComp = CalculateSelfComplementarity(sequence);
        if (selfComp > param.MaxSelfComplementarity)
        {
            score -= 0.2;
            warnings.Add($"High self-complementarity {selfComp:P0}");
        }

        // Check for secondary structure potential
        if (param.AvoidSecondaryStructure)
        {
            bool hasStructure = HasSecondaryStructurePotential(sequence);
            if (hasStructure)
            {
                score -= 0.15;
                warnings.Add("Potential secondary structure");
            }
        }

        // Check for repeats
        if (HasSimpleRepeats(sequence))
        {
            score -= 0.1;
            warnings.Add("Contains simple repeats");
        }

        // Penalize extreme positions
        double positionPenalty = 0;
        if (sequence.StartsWith("G") || sequence.StartsWith("C"))
            positionPenalty += 0.02;
        if (sequence.EndsWith("G") || sequence.EndsWith("C"))
            positionPenalty += 0.02;
        score -= positionPenalty;

        if (score <= 0)
            return null;

        return new Probe(
            sequence,
            start,
            start + sequence.Length - 1,
            tm,
            gc,
            Math.Max(0, score),
            ProbeType.Standard,
            warnings);
    }

    /// <summary>
    /// Designs tiling probes to cover entire sequence.
    /// </summary>
    public static TilingProbeSet DesignTilingProbes(
        string targetSequence,
        int probeLength = 60,
        int overlap = 20,
        ProbeParameters? parameters = null)
    {
        var param = parameters ?? Defaults.Microarray with
        {
            MinLength = probeLength,
            MaxLength = probeLength
        };

        targetSequence = targetSequence.ToUpperInvariant();
        var probes = new List<Probe>();
        int step = probeLength - overlap;

        for (int start = 0; start <= targetSequence.Length - probeLength; start += step)
        {
            string probeSeq = targetSequence.Substring(start, probeLength);
            var probe = EvaluateProbe(probeSeq, start, param);

            if (probe.HasValue)
            {
                probes.Add(probe.Value with { Type = ProbeType.Tiling });
            }
            else
            {
                // Add with warnings for coverage
                double tm = CalculateTm(probeSeq);
                double gc = CalculateGcContent(probeSeq);
                probes.Add(new Probe(
                    probeSeq, start, start + probeLength - 1,
                    tm, gc, 0.3, ProbeType.Tiling,
                    new List<string> { "Suboptimal probe, included for coverage" }));
            }
        }

        // Calculate coverage
        int covered = 0;
        var coveredPositions = new bool[targetSequence.Length];
        foreach (var probe in probes)
        {
            for (int i = probe.Start; i <= probe.End && i < targetSequence.Length; i++)
            {
                if (!coveredPositions[i])
                {
                    coveredPositions[i] = true;
                    covered++;
                }
            }
        }

        double meanTm = probes.Average(p => p.Tm);
        double tmRange = probes.Max(p => p.Tm) - probes.Min(p => p.Tm);

        return new TilingProbeSet(probes, covered, meanTm, tmRange);
    }

    /// <summary>
    /// Designs antisense probe for RNA detection.
    /// </summary>
    public static IEnumerable<Probe> DesignAntisenseProbes(
        string mRnaSequence,
        ProbeParameters? parameters = null,
        int maxProbes = 5)
    {
        // Get reverse complement for antisense probes
        string antisense = DnaSequence.GetReverseComplementString(mRnaSequence);

        foreach (var probe in DesignProbes(antisense, parameters, maxProbes))
        {
            yield return probe with { Type = ProbeType.Antisense };
        }
    }

    /// <summary>
    /// Designs molecular beacon probe (hairpin structure for real-time detection).
    /// </summary>
    public static Probe? DesignMolecularBeacon(
        string targetSequence,
        int probeLength = 25,
        int stemLength = 5)
    {
        if (targetSequence.Length < probeLength)
            return null;

        targetSequence = targetSequence.ToUpperInvariant();

        // Find best region in target
        double bestScore = 0;
        string? bestLoop = null;
        int bestStart = 0;

        int loopLength = probeLength;
        for (int start = 0; start <= targetSequence.Length - loopLength; start++)
        {
            string loop = targetSequence.Substring(start, loopLength);
            double gc = CalculateGcContent(loop);
            double tm = CalculateTm(loop);

            double score = 1.0;
            if (gc < 0.40 || gc > 0.60) score -= 0.2;
            if (tm < 55 || tm > 65) score -= 0.2;
            if (GetMaxHomopolymerLength(loop) > 4) score -= 0.2;

            if (score > bestScore)
            {
                bestScore = score;
                bestLoop = loop;
                bestStart = start;
            }
        }

        if (bestLoop == null)
            return null;

        // Add stem sequences (GC-rich for stability)
        string stem5 = new string('G', stemLength / 2) + new string('C', stemLength - stemLength / 2);
        string stem3 = DnaSequence.GetReverseComplementString(stem5);

        string beaconSequence = stem5 + bestLoop + stem3;
        double beaconTm = CalculateTm(bestLoop); // Loop Tm is target Tm

        return new Probe(
            beaconSequence,
            bestStart,
            bestStart + loopLength - 1,
            beaconTm,
            CalculateGcContent(beaconSequence),
            bestScore,
            ProbeType.MolecularBeacon,
            new List<string> { $"Stem: {stemLength}bp, Loop: {loopLength}bp" });
    }

    #endregion

    #region Probe Validation

    /// <summary>
    /// Validates probe against a genome/transcriptome.
    /// </summary>
    public static ProbeValidation ValidateProbe(
        string probeSequence,
        IEnumerable<string> referenceSequences,
        int maxMismatches = 3)
    {
        probeSequence = probeSequence.ToUpperInvariant();
        var issues = new List<string>();
        int offTargetHits = 0;

        // Check off-target hits
        foreach (var reference in referenceSequences)
        {
            var hits = FindApproximateMatches(reference.ToUpperInvariant(), probeSequence, maxMismatches);
            offTargetHits += hits.Count();
        }

        if (offTargetHits > 1)
        {
            issues.Add($"{offTargetHits} potential off-target sites");
        }

        // Check self-complementarity
        double selfComp = CalculateSelfComplementarity(probeSequence);
        if (selfComp > 0.3)
        {
            issues.Add($"Self-complementarity: {selfComp:P0}");
        }

        // Check secondary structure
        bool hasStructure = HasSecondaryStructurePotential(probeSequence);
        if (hasStructure)
        {
            issues.Add("Potential secondary structure formation");
        }

        // Calculate specificity score
        double specificity = offTargetHits <= 1 ? 1.0 : 1.0 / offTargetHits;

        bool isValid = issues.Count == 0 || (offTargetHits <= 1 && selfComp <= 0.4);

        return new ProbeValidation(
            isValid,
            specificity,
            offTargetHits,
            selfComp,
            hasStructure,
            issues);
    }

    /// <summary>
    /// Checks probe specificity using suffix tree (fast).
    /// </summary>
    public static double CheckSpecificity(
        string probeSequence,
        global::SuffixTree.ISuffixTree genomeIndex)
    {
        probeSequence = probeSequence.ToUpperInvariant();

        // Check if probe sequence exists in genome
        var positions = genomeIndex.FindAllOccurrences(probeSequence);
        int hitCount = positions.Count;

        if (hitCount == 0)
            return 0; // Probe doesn't match target

        if (hitCount == 1)
            return 1.0; // Unique match

        // Multiple hits reduce specificity
        return 1.0 / hitCount;
    }

    #endregion

    #region Oligo Analysis

    /// <summary>
    /// Analyzes oligonucleotide properties.
    /// </summary>
    public static (double Tm, double GcContent, double MolecularWeight, double ExtinctionCoefficient)
        AnalyzeOligo(string sequence)
    {
        sequence = sequence.ToUpperInvariant();

        double tm = CalculateTm(sequence);
        double gc = CalculateGcContent(sequence);
        double mw = CalculateMolecularWeight(sequence);
        double extinction = CalculateExtinctionCoefficient(sequence);

        return (tm, gc, mw, extinction);
    }

    /// <summary>
    /// Calculates molecular weight.
    /// </summary>
    public static double CalculateMolecularWeight(string sequence)
    {
        // Average molecular weights of nucleotides
        double weight = 0;
        foreach (char c in sequence.ToUpperInvariant())
        {
            weight += c switch
            {
                'A' => 331.2,
                'C' => 307.2,
                'G' => 347.2,
                'T' => 322.2,
                'U' => 308.2,
                _ => 330.0
            };
        }

        // Subtract water for each phosphodiester bond
        weight -= (sequence.Length - 1) * 18.0;

        return weight;
    }

    /// <summary>
    /// Calculates extinction coefficient at 260nm.
    /// </summary>
    public static double CalculateExtinctionCoefficient(string sequence)
    {
        // Nearest-neighbor method (simplified)
        double coefficient = 0;
        sequence = sequence.ToUpperInvariant();

        // Individual nucleotide contributions
        foreach (char c in sequence)
        {
            coefficient += c switch
            {
                'A' => 15400,
                'C' => 7400,
                'G' => 11500,
                'T' => 8700,
                'U' => 9900,
                _ => 10000
            };
        }

        return coefficient;
    }

    /// <summary>
    /// Calculates concentration from absorbance.
    /// </summary>
    public static double CalculateConcentration(
        double absorbance260,
        double extinctionCoefficient,
        double pathLength = 1.0)
    {
        // Beer-Lambert law: A = εcl
        // c = A / (ε * l)
        return absorbance260 / (extinctionCoefficient * pathLength) * 1e6; // µM
    }

    #endregion

    #region Helper Methods

    private static double CalculateGcContent(string sequence)
    {
        int gc = sequence.Count(c => c == 'G' || c == 'C');
        return sequence.Length > 0 ? gc / (double)sequence.Length : 0;
    }

    private static double CalculateTm(string sequence)
    {
        int length = sequence.Length;

        if (length < 14)
        {
            // Wallace rule for short oligos
            int at = sequence.Count(c => c == 'A' || c == 'T');
            int gc = sequence.Count(c => c == 'G' || c == 'C');
            return 2 * at + 4 * gc;
        }
        else
        {
            // Salt-adjusted formula
            double gc = CalculateGcContent(sequence);
            return 81.5 + 16.6 * Math.Log10(0.05) + 41 * gc - 600.0 / length;
        }
    }

    private static int GetMaxHomopolymerLength(string sequence)
    {
        int maxRun = 1;
        int currentRun = 1;

        for (int i = 1; i < sequence.Length; i++)
        {
            if (sequence[i] == sequence[i - 1])
            {
                currentRun++;
                maxRun = Math.Max(maxRun, currentRun);
            }
            else
            {
                currentRun = 1;
            }
        }

        return maxRun;
    }

    private static double CalculateSelfComplementarity(string sequence)
    {
        string revComp = DnaSequence.GetReverseComplementString(sequence);
        int matches = 0;

        for (int i = 0; i < sequence.Length; i++)
        {
            if (sequence[i] == revComp[i])
                matches++;
        }

        return matches / (double)sequence.Length;
    }

    private static bool HasSecondaryStructurePotential(string sequence)
    {
        // Check for inverted repeats that could form hairpins
        int halfLen = sequence.Length / 2;

        for (int stemLen = 4; stemLen <= halfLen; stemLen++)
        {
            for (int i = 0; i <= sequence.Length - stemLen * 2 - 3; i++)
            {
                string left = sequence.Substring(i, stemLen);
                string right = sequence.Substring(i + stemLen + 3, stemLen);
                string rightRC = DnaSequence.GetReverseComplementString(right);

                int matches = 0;
                for (int j = 0; j < stemLen; j++)
                {
                    if (left[j] == rightRC[j])
                        matches++;
                }

                if (matches >= stemLen * 0.8)
                    return true;
            }
        }

        return false;
    }

    private static bool HasSimpleRepeats(string sequence)
    {
        // Check for di/tri-nucleotide repeats
        for (int unitLen = 2; unitLen <= 3; unitLen++)
        {
            for (int i = 0; i <= sequence.Length - unitLen * 4; i++)
            {
                string unit = sequence.Substring(i, unitLen);
                int repeats = 1;

                for (int j = i + unitLen; j <= sequence.Length - unitLen; j += unitLen)
                {
                    if (sequence.Substring(j, unitLen) == unit)
                        repeats++;
                    else
                        break;
                }

                if (repeats >= 4)
                    return true;
            }
        }

        return false;
    }

    private static IEnumerable<int> FindApproximateMatches(
        string text, string pattern, int maxMismatches)
    {
        for (int i = 0; i <= text.Length - pattern.Length; i++)
        {
            int mismatches = 0;
            for (int j = 0; j < pattern.Length && mismatches <= maxMismatches; j++)
            {
                if (text[i + j] != pattern[j])
                    mismatches++;
            }

            if (mismatches <= maxMismatches)
                yield return i;
        }
    }

    #endregion
}
