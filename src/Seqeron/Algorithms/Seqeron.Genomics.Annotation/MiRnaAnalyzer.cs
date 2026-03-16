using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Provides algorithms for microRNA analysis, target prediction, and seed matching.
/// </summary>
public static class MiRnaAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Represents a microRNA.
    /// </summary>
    public readonly record struct MiRna(
        string Name,
        string Sequence,
        string SeedSequence,
        int SeedStart,
        int SeedEnd);

    /// <summary>
    /// Represents a predicted miRNA target site.
    /// </summary>
    public readonly record struct TargetSite(
        int Start,
        int End,
        string TargetSequence,
        string MiRnaName,
        TargetSiteType Type,
        int SeedMatchLength,
        double Score,
        double FreeEnergy,
        string Alignment);

    /// <summary>
    /// Types of miRNA target sites.
    /// </summary>
    public enum TargetSiteType
    {
        Seed8mer,      // 8mer: perfect seed + A at position 1
        Seed7merM8,    // 7mer-m8: perfect seed match positions 2-8
        Seed7merA1,    // 7mer-A1: positions 2-7 + A at position 1
        Seed6mer,      // 6mer: positions 2-7
        Offset6mer,    // offset 6mer: positions 3-8
        Supplementary, // 3' supplementary pairing
        Centered       // centered site (positions 4-15)
    }

    /// <summary>
    /// Represents a potential precursor miRNA (pre-miRNA) hairpin.
    /// </summary>
    public readonly record struct PreMiRna(
        int Start,
        int End,
        string Sequence,
        string MatureSequence,
        string StarSequence,
        string Structure,
        double FreeEnergy);

    /// <summary>
    /// miRNA-mRNA duplex with alignment details.
    /// </summary>
    public readonly record struct MiRnaDuplex(
        string MiRnaSequence,
        string TargetSequence,
        string AlignmentString,
        int Matches,
        int Mismatches,
        int GUWobbles,
        int Gaps,
        double FreeEnergy);

    /// <summary>
    /// Result of comparing two miRNA seed regions.
    /// </summary>
    public readonly record struct SeedComparison(
        int Matches,
        int Mismatches,
        bool IsSameFamily);

    #endregion

    #region Seed Matching

    /// <summary>
    /// Extracts the seed region from a miRNA sequence (positions 2-8).
    /// </summary>
    public static string GetSeedSequence(string miRnaSequence)
    {
        if (string.IsNullOrEmpty(miRnaSequence) || miRnaSequence.Length < 8)
            return "";

        return miRnaSequence.Substring(1, 7).ToUpperInvariant();
    }

    /// <summary>
    /// Creates a MiRna record from a sequence.
    /// </summary>
    public static MiRna CreateMiRna(string name, string sequence)
    {
        string upper = sequence.ToUpperInvariant().Replace('T', 'U');
        string seed = GetSeedSequence(upper);

        return new MiRna(
            Name: name,
            Sequence: upper,
            SeedSequence: seed,
            SeedStart: 1,
            SeedEnd: 7);
    }

    /// <summary>
    /// Compares the seed regions of two miRNAs, returning the number of matches,
    /// mismatches (Hamming distance), and whether they belong to the same seed family.
    /// </summary>
    public static SeedComparison CompareSeedRegions(MiRna mirna1, MiRna mirna2)
    {
        string seed1 = mirna1.SeedSequence;
        string seed2 = mirna2.SeedSequence;

        if (string.IsNullOrEmpty(seed1) || string.IsNullOrEmpty(seed2))
            return new SeedComparison(Matches: 0, Mismatches: 0, IsSameFamily: false);

        int length = Math.Min(seed1.Length, seed2.Length);
        int matches = 0;
        int mismatches = 0;

        for (int i = 0; i < length; i++)
        {
            if (seed1[i] == seed2[i])
                matches++;
            else
                mismatches++;
        }

        // Account for length differences (if seeds have different lengths)
        mismatches += Math.Abs(seed1.Length - seed2.Length);

        bool isSameFamily = seed1 == seed2;

        return new SeedComparison(Matches: matches, Mismatches: mismatches, IsSameFamily: isSameFamily);
    }

    /// <summary>
    /// Finds all potential target sites for a miRNA in an mRNA sequence.
    /// Scans for the 6mer core (RC of miRNA positions 2-7), then extends to classify
    /// site types per Bartel (2009) and TargetScan conventions.
    /// </summary>
    public static IEnumerable<TargetSite> FindTargetSites(
        string mRnaSequence,
        MiRna miRna,
        double minScore = 0.5)
    {
        if (string.IsNullOrEmpty(mRnaSequence) || string.IsNullOrEmpty(miRna.Sequence))
            yield break;

        string mrna = mRnaSequence.ToUpperInvariant().Replace('T', 'U');
        string mirna = miRna.Sequence;

        string seedRC = GetReverseComplement(miRna.SeedSequence);
        if (seedRC.Length < 7)
            yield break;

        // seedRC layout (7 chars): [RC of pos8] [RC of pos7] ... [RC of pos2]
        // 6mer core = seedRC[1..7] = RC of miRNA positions 2-7 (6 chars)
        // offset 6mer pattern = seedRC[0..6] = RC of miRNA positions 3-8 (6 chars)
        string sixmerCore = seedRC.Substring(1, 6);
        string offset6Pat = seedRC.Substring(0, 6);

        // Track positions covered by 6mer-core-based sites to suppress overlapping offset 6mer
        var coveredPositions = new HashSet<int>();

        // Pass 1: Find 6mer core matches → classify as 8mer/7mer-m8/7mer-A1/6mer
        for (int i = 0; i <= mrna.Length - 6; i++)
        {
            if (mrna.Substring(i, 6) != sixmerCore)
                continue;

            // Check for position 8 match upstream: mrna[i-1] == seedRC[0]
            bool hasPos8 = i > 0 && mrna[i - 1] == seedRC[0];
            // Check for A opposite miRNA position 1 downstream: mrna[i+6] == 'A'
            bool hasA1 = i + 6 < mrna.Length && mrna[i + 6] == 'A';

            TargetSiteType type;
            int siteStart;
            int siteLength;
            int seedMatchLen;

            if (hasPos8 && hasA1)
            {
                // 8mer: match to positions 2-8 + A opposite position 1
                type = TargetSiteType.Seed8mer;
                siteStart = i - 1;
                siteLength = 8;
                seedMatchLen = 8;
            }
            else if (hasPos8)
            {
                // 7mer-m8: match to positions 2-8
                type = TargetSiteType.Seed7merM8;
                siteStart = i - 1;
                siteLength = 7;
                seedMatchLen = 7;
            }
            else if (hasA1)
            {
                // 7mer-A1: match to positions 2-7 + A opposite position 1
                type = TargetSiteType.Seed7merA1;
                siteStart = i;
                siteLength = 7;
                seedMatchLen = 7;
            }
            else
            {
                // 6mer: match to positions 2-7 only
                type = TargetSiteType.Seed6mer;
                siteStart = i;
                siteLength = 6;
                seedMatchLen = 6;
            }

            var site = CreateTargetSite(mrna, siteStart, siteLength, mirna, miRna.Name, type, seedMatchLen);
            if (site.Score >= minScore)
            {
                for (int j = siteStart; j < siteStart + siteLength; j++)
                    coveredPositions.Add(j);
                yield return site;
            }
        }

        // Pass 2: Find offset 6mer matches (positions 3-8) not overlapping higher-priority sites
        for (int i = 0; i <= mrna.Length - 6; i++)
        {
            if (mrna.Substring(i, 6) != offset6Pat)
                continue;

            // Skip if this position overlaps with a site already found
            bool overlaps = false;
            for (int j = i; j < i + 6; j++)
            {
                if (coveredPositions.Contains(j))
                {
                    overlaps = true;
                    break;
                }
            }
            if (overlaps) continue;

            // Also skip if full seedRC matches here (would have been found in pass 1 at i+1)
            if (i + 7 <= mrna.Length && mrna[i + 6] == seedRC[6])
                continue;

            var site = CreateTargetSite(mrna, i, 6, mirna, miRna.Name, TargetSiteType.Offset6mer, 6);
            if (site.Score >= minScore)
            {
                yield return site;
            }
        }
    }

    private static TargetSite CreateTargetSite(string mrna, int pos, int length, string mirna, string mirnaName, TargetSiteType type, int seedMatchLength)
    {
        // Extend alignment for full miRNA
        int extendedLength = Math.Min(mirna.Length, mrna.Length - pos);
        string targetSeq = mrna.Substring(pos, extendedLength);

        var duplex = AlignMiRnaToTarget(mirna, targetSeq);
        double score = CalculateTargetScore(type, duplex);

        return new TargetSite(
            Start: pos,
            End: pos + length - 1,
            TargetSequence: targetSeq,
            MiRnaName: mirnaName,
            Type: type,
            SeedMatchLength: seedMatchLength,
            Score: score,
            FreeEnergy: duplex.FreeEnergy,
            Alignment: duplex.AlignmentString);
    }

    #endregion

    #region Sequence Utilities

    /// <summary>
    /// Gets the reverse complement of an RNA sequence.
    /// </summary>
    public static string GetReverseComplement(string rnaSequence)
    {
        if (string.IsNullOrEmpty(rnaSequence))
            return "";

        var sb = new StringBuilder(rnaSequence.Length);

        for (int i = rnaSequence.Length - 1; i >= 0; i--)
        {
            char c = char.ToUpperInvariant(rnaSequence[i]);
            char complement = c switch
            {
                'A' => 'U',
                'U' => 'A',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => 'N'
            };
            sb.Append(complement);
        }

        return sb.ToString();
    }

    /// <summary>
    /// Checks if two bases can pair (Watson-Crick or G-U wobble).
    /// </summary>
    public static bool CanPair(char base1, char base2)
    {
        char b1 = char.ToUpperInvariant(base1);
        char b2 = char.ToUpperInvariant(base2);

        return (b1 == 'A' && b2 == 'U') || (b1 == 'U' && b2 == 'A') ||
               (b1 == 'G' && b2 == 'C') || (b1 == 'C' && b2 == 'G') ||
               (b1 == 'G' && b2 == 'U') || (b1 == 'U' && b2 == 'G');
    }

    /// <summary>
    /// Checks if a pair is a G-U wobble.
    /// </summary>
    public static bool IsWobblePair(char base1, char base2)
    {
        char b1 = char.ToUpperInvariant(base1);
        char b2 = char.ToUpperInvariant(base2);

        return (b1 == 'G' && b2 == 'U') || (b1 == 'U' && b2 == 'G');
    }

    #endregion

    #region Alignment and Scoring

    /// <summary>
    /// Aligns a miRNA to a target sequence.
    /// </summary>
    public static MiRnaDuplex AlignMiRnaToTarget(string miRnaSequence, string targetSequence)
    {
        if (string.IsNullOrEmpty(miRnaSequence) || string.IsNullOrEmpty(targetSequence))
        {
            return new MiRnaDuplex("", "", "", 0, 0, 0, 0, 0);
        }

        string mirna = miRnaSequence.ToUpperInvariant().Replace('T', 'U');
        string target = targetSequence.ToUpperInvariant().Replace('T', 'U');

        // Reverse target for 5'-3' alignment (miRNA binds 3' UTR in reverse orientation)
        var alignmentChars = new char[Math.Min(mirna.Length, target.Length)];
        int matches = 0, mismatches = 0, wobbles = 0;

        for (int i = 0; i < alignmentChars.Length; i++)
        {
            char m = mirna[i];
            // Target is read in reverse (5' to 3' of target = 3' to 5' of miRNA binding)
            int targetIdx = target.Length - 1 - i;
            if (targetIdx < 0) break;

            char t = target[targetIdx];

            if (CanPair(m, t))
            {
                if (IsWobblePair(m, t))
                {
                    alignmentChars[i] = ':';
                    wobbles++;
                }
                else
                {
                    alignmentChars[i] = '|';
                    matches++;
                }
            }
            else
            {
                alignmentChars[i] = ' ';
                mismatches++;
            }
        }

        double freeEnergy = CalculateDuplexEnergy(mirna, target, alignmentChars);

        return new MiRnaDuplex(
            MiRnaSequence: mirna,
            TargetSequence: target,
            AlignmentString: new string(alignmentChars),
            Matches: matches,
            Mismatches: mismatches,
            GUWobbles: wobbles,
            Gaps: 0,
            FreeEnergy: freeEnergy);
    }

    private static double CalculateDuplexEnergy(string mirna, string target, char[] alignment)
    {
        // Simplified energy calculation
        double energy = 0;

        for (int i = 0; i < alignment.Length - 1; i++)
        {
            if (alignment[i] == '|' && i + 1 < alignment.Length && alignment[i + 1] == '|')
            {
                energy -= 2.0; // Stacking bonus
            }
            else if (alignment[i] == '|')
            {
                energy -= 1.0; // Single match
            }
            else if (alignment[i] == ':')
            {
                energy -= 0.5; // Wobble
            }
            else
            {
                energy += 0.5; // Mismatch penalty
            }
        }

        return energy;
    }

    private static double CalculateTargetScore(TargetSiteType type, MiRnaDuplex duplex)
    {
        // Base scores proportional to Grimson et al. (2007) site-type efficacy weights:
        // 8mer=0.310, 7mer-m8=0.161, 7mer-A1=0.099; normalized to [0,1] with 8mer=1.0
        // 6mer and offset 6mer estimated from Agarwal et al. (2015).
        double baseScore = type switch
        {
            TargetSiteType.Seed8mer => 1.0,
            TargetSiteType.Seed7merM8 => 0.52,
            TargetSiteType.Seed7merA1 => 0.32,
            TargetSiteType.Seed6mer => 0.15,
            TargetSiteType.Offset6mer => 0.10,
            TargetSiteType.Supplementary => 0.05,
            TargetSiteType.Centered => 0.03,
            _ => 0.01
        };

        // Bonus for 3' supplementary pairing (extra matches beyond seed)
        if (duplex.Matches > 10)
        {
            baseScore += 0.05;
        }

        // Penalty for mismatches in alignment
        baseScore -= duplex.Mismatches * 0.01;

        return Math.Max(0, Math.Min(1, baseScore));
    }

    #endregion

    #region Pre-miRNA Energy Parameters

    // Turner 2004 nearest-neighbor stacking energies (kcal/mol at 37°C)
    // Source: NNDB — https://rna.urmc.rochester.edu/NNDB/turner04/
    // Format: 5'-XY-3' / 3'-X'Y'-5' where X pairs with X', Y pairs with Y'
    private static readonly Dictionary<string, double> StackingEnergies = new()
    {
        // Watson-Crick stacking (16 entries)
        { "AA/UU", -0.93 }, { "UU/AA", -0.93 },
        { "AU/UA", -1.10 },
        { "UA/AU", -1.33 },
        { "CU/GA", -2.08 }, { "AG/UC", -2.08 },
        { "CA/GU", -2.11 }, { "UG/AC", -2.11 },
        { "GU/CA", -2.24 }, { "AC/UG", -2.24 },
        { "GA/CU", -2.35 }, { "UC/AG", -2.35 },
        { "CG/GC", -2.36 },
        { "GG/CC", -3.26 }, { "CC/GG", -3.26 },
        { "GC/CG", -3.42 },
        // GU wobble stacking (20 entries)
        { "AG/UU", -0.55 }, { "UU/GA", -0.55 },
        { "UG/AU", -1.00 }, { "UA/GU", -1.00 },
        { "GA/UU", -1.27 }, { "UU/AG", -1.27 },
        { "AU/UG", -1.36 }, { "GU/UA", -1.36 },
        { "CG/GU", -1.41 }, { "UG/GC", -1.41 },
        { "GG/CU", -1.53 }, { "UC/GG", -1.53 },
        { "CU/GG", -2.11 }, { "GG/UC", -2.11 },
        { "GU/CG", -2.51 }, { "GC/UG", -2.51 },
        { "GG/UU", -0.50 }, { "UU/GG", -0.50 },
        { "UG/GU", +0.30 },
        { "GU/UG", +1.29 },
    };

    // Hairpin loop initiation energies (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/loop.txt
    // Sizes 3-9 experimentally determined; 10-30 extrapolated via
    // ΔG°(n) = ΔG°(9) + 1.75·R·T·ln(n/9)
    private static readonly Dictionary<int, double> HairpinLoopInitEnergies = new()
    {
        {  3, 5.4 }, {  4, 5.6 }, {  5, 5.7 }, {  6, 5.4 }, {  7, 6.0 },
        {  8, 5.5 }, {  9, 6.4 }, { 10, 6.5 }, { 11, 6.6 }, { 12, 6.7 },
        { 13, 6.8 }, { 14, 6.9 }, { 15, 6.9 }, { 16, 7.0 }, { 17, 7.1 },
        { 18, 7.1 }, { 19, 7.2 }, { 20, 7.2 }, { 21, 7.3 }, { 22, 7.3 },
        { 23, 7.4 }, { 24, 7.4 }, { 25, 7.5 }, { 26, 7.5 }, { 27, 7.5 },
        { 28, 7.6 }, { 29, 7.6 }, { 30, 7.7 }
    };

    // Terminal mismatch stacking energies (kcal/mol at 37°C)
    // Source: NNDB — rna.urmc.rochester.edu/NNDB/turner04/tm-parameters.html
    // Key = closing5' + firstMismatch(5'side) + lastMismatch(3'side) + closing3'
    private static readonly Dictionary<string, double> TerminalMismatchEnergies = new()
    {
        // Closing pair A-U
        { "AAAU", -0.8 }, { "AACU", -1.0 }, { "AAGU", -0.8 }, { "AAUU", -1.0 },
        { "ACAU", -0.6 }, { "ACCU", -0.7 }, { "ACGU", -0.6 }, { "ACUU", -0.7 },
        { "AGAU", -0.8 }, { "AGCU", -1.0 }, { "AGGU", -0.8 }, { "AGUU", -1.0 },
        { "AUAU", -0.6 }, { "AUCU", -0.8 }, { "AUGU", -0.6 }, { "AUUU", -0.8 },
        // Closing pair C-G
        { "CAAG", -1.5 }, { "CACG", -1.5 }, { "CAGG", -1.4 }, { "CAUG", -1.5 },
        { "CCAG", -1.0 }, { "CCCG", -1.1 }, { "CCGG", -1.0 }, { "CCUG", -0.8 },
        { "CGAG", -1.4 }, { "CGCG", -1.5 }, { "CGGG", -1.6 }, { "CGUG", -1.5 },
        { "CUAG", -1.0 }, { "CUCG", -1.4 }, { "CUGG", -1.0 }, { "CUUG", -1.2 },
        // Closing pair G-C
        { "GAAC", -1.1 }, { "GACC", -1.5 }, { "GAGC", -1.3 }, { "GAUC", -1.5 },
        { "GCAC", -1.1 }, { "GCCC", -0.7 }, { "GCGC", -1.1 }, { "GCUC", -0.5 },
        { "GGAC", -1.6 }, { "GGCC", -1.5 }, { "GGGC", -1.4 }, { "GGUC", -1.5 },
        { "GUAC", -1.1 }, { "GUCC", -1.0 }, { "GUGC", -1.1 }, { "GUUC", -0.7 },
        // Closing pair G-U
        { "GAAU", -0.3 }, { "GACU", -1.0 }, { "GAGU", -0.8 }, { "GAUU", -1.0 },
        { "GCAU", -0.6 }, { "GCCU", -0.7 }, { "GCGU", -0.6 }, { "GCUU", -0.7 },
        { "GGAU", -0.6 }, { "GGCU", -1.0 }, { "GGGU", -0.8 }, { "GGUU", -1.0 },
        { "GUAU", -0.6 }, { "GUCU", -0.8 }, { "GUGU", -0.6 }, { "GUUU", -0.6 },
        // Closing pair U-A
        { "UAAA", -1.0 }, { "UACA", -0.8 }, { "UAGA", -1.1 }, { "UAUA", -0.8 },
        { "UCAA", -0.7 }, { "UCCA", -0.6 }, { "UCGA", -0.7 }, { "UCUA", -0.5 },
        { "UGAA", -1.1 }, { "UGCA", -0.8 }, { "UGGA", -1.2 }, { "UGUA", -0.8 },
        { "UUAA", -0.7 }, { "UUCA", -0.6 }, { "UUGA", -0.7 }, { "UUUA", -0.5 },
        // Closing pair U-G
        { "UAAG", -1.0 }, { "UACG", -0.8 }, { "UAGG", -1.1 }, { "UAUG", -0.8 },
        { "UCAG", -0.7 }, { "UCCG", -0.6 }, { "UCGG", -0.7 }, { "UCUG", -0.5 },
        { "UGAG", -0.5 }, { "UGCG", -0.8 }, { "UGGG", -0.8 }, { "UGUG", -0.8 },
        { "UUAG", -0.7 }, { "UUCG", -0.6 }, { "UUGG", -0.7 }, { "UUUG", -0.5 },
    };

    // Terminal AU/GU penalty (kcal/mol at 37°C)
    // Source: NNDB — "Per AU end" on wc-parameters.html, "Per GU end" on gu-parameters.html
    private const double TerminalAU_GU_Penalty = 0.45;

    /// <summary>
    /// Calculates hairpin free energy using Turner 2004 nearest-neighbor parameters.
    /// Components: stacking energies + loop initiation + terminal mismatch + AU/GU penalties.
    /// </summary>
    private static double CalculateHairpinEnergy(string seq, int stemLength)
    {
        int n = seq.Length;
        int loopSize = n - 2 * stemLength;
        double energy = 0.0;

        // 1. Sum stacking energies for consecutive base pairs in the stem
        for (int i = 0; i < stemLength - 1; i++)
        {
            string key = $"{seq[i]}{seq[i + 1]}/{seq[n - 1 - i]}{seq[n - 2 - i]}";
            if (StackingEnergies.TryGetValue(key, out double stackE))
                energy += stackE;
        }

        // 2. Hairpin loop initiation
        if (HairpinLoopInitEnergies.TryGetValue(loopSize, out double loopE))
        {
            energy += loopE;
        }
        else if (loopSize > 30)
        {
            // Extrapolation: ΔG°(n) = ΔG°(9) + 1.75·R·T·ln(n/9)
            // R = 0.001987 kcal/(mol·K), T = 310.15 K → 1.75·R·T ≈ 1.079
            energy += 6.4 + 1.079 * Math.Log((double)loopSize / 9);
        }

        // 3. Terminal mismatch (for loops ≥ 4 nt)
        if (loopSize >= 4)
        {
            string tmKey = $"{seq[stemLength - 1]}{seq[stemLength]}{seq[n - stemLength - 1]}{seq[n - stemLength]}";
            if (TerminalMismatchEnergies.TryGetValue(tmKey, out double tmE))
                energy += tmE;
        }

        // 4. Terminal AU/GU penalties (outer and closing pairs)
        if (IsTerminalPenaltyPair(seq[0], seq[n - 1]))
            energy += TerminalAU_GU_Penalty;
        if (IsTerminalPenaltyPair(seq[stemLength - 1], seq[n - stemLength]))
            energy += TerminalAU_GU_Penalty;

        return energy;
    }

    private static bool IsTerminalPenaltyPair(char b1, char b2) =>
        (b1 == 'A' && b2 == 'U') || (b1 == 'U' && b2 == 'A') ||
        (b1 == 'G' && b2 == 'U') || (b1 == 'U' && b2 == 'G');

    #endregion

    #region Pre-miRNA Prediction

    /// <summary>
    /// Finds potential pre-miRNA hairpin structures in a sequence.
    /// </summary>
    public static IEnumerable<PreMiRna> FindPreMiRnaHairpins(
        string sequence,
        int minHairpinLength = 55,
        int maxHairpinLength = 120,
        int matureLength = 22)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minHairpinLength)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // Scan for potential hairpins
        for (int i = 0; i <= upper.Length - minHairpinLength; i++)
        {
            for (int len = minHairpinLength; len <= Math.Min(maxHairpinLength, upper.Length - i); len++)
            {
                string candidate = upper.Substring(i, len);

                // Check for hairpin structure
                var hairpin = AnalyzeHairpin(candidate, matureLength);
                if (hairpin != null)
                {
                    yield return new PreMiRna(
                        Start: i,
                        End: i + len - 1,
                        Sequence: candidate,
                        MatureSequence: hairpin.Value.Mature,
                        StarSequence: hairpin.Value.Star,
                        Structure: hairpin.Value.Structure,
                        FreeEnergy: hairpin.Value.Energy);
                }
            }
        }
    }

    private static (string Mature, string Star, string Structure, double Energy)? AnalyzeHairpin(string sequence, int matureLength)
    {
        int n = sequence.Length;
        if (n < 55) return null;

        // Find stem region by looking for complementary ends
        int stemLength = 0;
        int maxStem = Math.Min(n / 2 - 5, 35); // Leave room for loop

        for (int j = 0; j < maxStem; j++)
        {
            if (CanPair(sequence[j], sequence[n - 1 - j]))
            {
                stemLength++;
            }
            else
            {
                break;
            }
        }

        if (stemLength < 18) // Pre-miRNA needs ~18+ bp stem
            return null;

        // Check loop size (should be 3-15 nt)
        int loopSize = n - 2 * stemLength;
        if (loopSize < 3 || loopSize > 25)
            return null;

        // Extract mature miRNA (typically from 5' arm)
        int matureStart = 0;
        int matureEnd = Math.Min(matureLength, stemLength);
        string mature = sequence.Substring(matureStart, matureEnd);

        // Star sequence from 3' arm
        int starStart = n - matureEnd;
        string star = sequence.Substring(starStart, matureEnd);

        // Generate dot-bracket structure
        var structure = new char[n];
        for (int j = 0; j < n; j++)
        {
            if (j < stemLength)
                structure[j] = '(';
            else if (j >= n - stemLength)
                structure[j] = ')';
            else
                structure[j] = '.';
        }

        // Calculate free energy using Turner 2004 nearest-neighbor parameters
        double energy = CalculateHairpinEnergy(sequence, stemLength);

        return (mature, star, new string(structure), energy);
    }

    #endregion

    #region Target Context Analysis

    /// <summary>
    /// Analyzes the context around a target site (AU content, position in 3'UTR, etc.).
    /// </summary>
    public static (double AuContent, bool NearStart, bool NearEnd, double ContextScore) AnalyzeTargetContext(
        string mRnaSequence,
        int targetStart,
        int targetEnd,
        int contextWindow = 30)
    {
        if (string.IsNullOrEmpty(mRnaSequence))
            return (0, false, false, 0);

        string mrna = mRnaSequence.ToUpperInvariant();

        // Get context window
        int windowStart = Math.Max(0, targetStart - contextWindow);
        int windowEnd = Math.Min(mrna.Length, targetEnd + contextWindow);
        string context = mrna.Substring(windowStart, windowEnd - windowStart);

        // Calculate AU content
        int auCount = context.Count(c => c == 'A' || c == 'U');
        double auContent = (double)auCount / context.Length;

        // Check position
        bool nearStart = targetStart < mrna.Length * 0.15;
        bool nearEnd = targetEnd > mrna.Length * 0.85;

        // Calculate context score (higher AU content near target is favorable)
        double contextScore = auContent * 0.5;

        // Bonus for not being at very end
        if (!nearEnd && !nearStart)
        {
            contextScore += 0.3;
        }

        return (auContent, nearStart, nearEnd, Math.Min(1.0, contextScore));
    }

    /// <summary>
    /// Checks for site accessibility based on local structure.
    /// </summary>
    public static double CalculateSiteAccessibility(string mRnaSequence, int siteStart, int siteEnd)
    {
        if (string.IsNullOrEmpty(mRnaSequence) || siteStart < 0 || siteEnd >= mRnaSequence.Length)
            return 0;

        // Simplified accessibility: check for self-complementarity in local region
        int windowSize = 50;
        int start = Math.Max(0, siteStart - windowSize);
        int end = Math.Min(mRnaSequence.Length, siteEnd + windowSize);

        string window = mRnaSequence.Substring(start, end - start).ToUpperInvariant();

        // Count potential base pairs in the window (indicating structure)
        int structureScore = 0;
        for (int i = 0; i < window.Length; i++)
        {
            for (int j = i + 4; j < window.Length; j++)
            {
                if (CanPair(window[i], window[j]) && !IsWobblePair(window[i], window[j]))
                {
                    structureScore++;
                }
            }
        }

        // Higher structure score = less accessible
        double maxPairs = (window.Length * (window.Length - 4)) / 2;
        double structureDensity = structureScore / Math.Max(1, maxPairs);

        return Math.Max(0, 1.0 - structureDensity * 10);
    }

    #endregion

    #region miRNA Family Analysis

    /// <summary>
    /// Groups miRNAs by their seed sequence family.
    /// </summary>
    public static IEnumerable<(string SeedFamily, IReadOnlyList<MiRna> Members)> GroupBySeedFamily(IEnumerable<MiRna> miRnas)
    {
        return miRnas
            .GroupBy(m => m.SeedSequence)
            .Select(g => (g.Key, (IReadOnlyList<MiRna>)g.ToList()));
    }

    /// <summary>
    /// Finds miRNAs with similar seed sequences.
    /// </summary>
    public static IEnumerable<MiRna> FindSimilarMiRnas(MiRna query, IEnumerable<MiRna> database, int maxMismatches = 1)
    {
        string querySeed = query.SeedSequence;

        foreach (var mirna in database)
        {
            if (mirna.Name == query.Name)
                continue;

            int mismatches = 0;
            for (int i = 0; i < Math.Min(querySeed.Length, mirna.SeedSequence.Length); i++)
            {
                if (querySeed[i] != mirna.SeedSequence[i])
                    mismatches++;
            }

            if (mismatches <= maxMismatches)
            {
                yield return mirna;
            }
        }
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Calculates the GC content of a sequence.
    /// </summary>
    public static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcFractionFast();

    /// <summary>
    /// Generates all possible seed sequences for a given miRNA.
    /// </summary>
    public static IEnumerable<string> GenerateSeedVariants(string seedSequence, bool includeWobble = true)
    {
        yield return seedSequence;

        // Generate single-nucleotide variants at each position
        string bases = includeWobble ? "ACGU" : "ACGU";

        for (int i = 0; i < seedSequence.Length; i++)
        {
            foreach (char b in bases)
            {
                if (b != seedSequence[i])
                {
                    char[] variant = seedSequence.ToCharArray();
                    variant[i] = b;
                    yield return new string(variant);
                }
            }
        }
    }

    #endregion
}
