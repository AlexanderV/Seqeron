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

        // Get reverse complement of miRNA seed for target matching
        string seedRC = GetReverseComplement(miRna.SeedSequence);

        // Scan mRNA for seed matches
        for (int i = 0; i <= mrna.Length - 6; i++)
        {
            // Check for different seed match types
            var site = CheckSeedMatch(mrna, i, mirna, miRna.Name, seedRC);
            if (site != null && site.Value.Score >= minScore)
            {
                yield return site.Value;
            }
        }
    }

    private static TargetSite? CheckSeedMatch(string mrna, int pos, string mirna, string mirnaName, string seedRC)
    {
        if (pos + 8 > mrna.Length)
            return null;

        string target8 = mrna.Substring(pos, Math.Min(8, mrna.Length - pos));

        // Check for 8mer (A at position opposite to miRNA pos 1, then perfect seed)
        if (target8.Length >= 8)
        {
            string seedMatch = target8.Substring(0, 7);
            if (seedMatch == seedRC && target8[7] == 'A')
            {
                return CreateTargetSite(mrna, pos, 8, mirna, mirnaName, TargetSiteType.Seed8mer, 8);
            }
        }

        // Check for 7mer-m8 (positions 2-8 of miRNA, i.e., seed + position 8)
        if (target8.Length >= 7)
        {
            string m8Match = target8.Substring(0, 7);
            if (m8Match == seedRC)
            {
                return CreateTargetSite(mrna, pos, 7, mirna, mirnaName, TargetSiteType.Seed7merM8, 7);
            }
        }

        // Check for 7mer-A1 (positions 2-7 + A opposite position 1)
        if (target8.Length >= 7)
        {
            string a1Match = target8.Substring(0, 6);
            if (a1Match == seedRC.Substring(0, 6) && pos > 0 && mrna[pos - 1] == 'A')
            {
                return CreateTargetSite(mrna, pos - 1, 7, mirna, mirnaName, TargetSiteType.Seed7merA1, 7);
            }
        }

        // Check for 6mer (positions 2-7)
        if (target8.Length >= 6)
        {
            string match6 = target8.Substring(0, 6);
            if (match6 == seedRC.Substring(0, 6))
            {
                return CreateTargetSite(mrna, pos, 6, mirna, mirnaName, TargetSiteType.Seed6mer, 6);
            }
        }

        // Check for offset 6mer (positions 3-8)
        if (target8.Length >= 6 && seedRC.Length >= 7)
        {
            string offset6 = target8.Substring(0, 6);
            if (offset6 == seedRC.Substring(1, 6))
            {
                return CreateTargetSite(mrna, pos, 6, mirna, mirnaName, TargetSiteType.Offset6mer, 6);
            }
        }

        return null;
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
        double baseScore = type switch
        {
            TargetSiteType.Seed8mer => 1.0,
            TargetSiteType.Seed7merM8 => 0.9,
            TargetSiteType.Seed7merA1 => 0.85,
            TargetSiteType.Seed6mer => 0.7,
            TargetSiteType.Offset6mer => 0.6,
            TargetSiteType.Supplementary => 0.5,
            TargetSiteType.Centered => 0.4,
            _ => 0.3
        };

        // Bonus for 3' supplementary pairing
        if (duplex.Matches > 10)
        {
            baseScore += 0.1;
        }

        // Penalty for mismatches in seed
        baseScore -= duplex.Mismatches * 0.05;

        return Math.Max(0, Math.Min(1, baseScore));
    }

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

        // Calculate approximate free energy
        double energy = -stemLength * 1.5 + loopSize * 0.5;

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
