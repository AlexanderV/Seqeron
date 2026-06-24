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
    /// Inputs are case-insensitive and DNA T is treated as RNA U.
    /// </summary>
    public static bool CanPair(char base1, char base2)
    {
        char b1 = NormalizeBase(base1);
        char b2 = NormalizeBase(base2);

        return (b1 == 'A' && b2 == 'U') || (b1 == 'U' && b2 == 'A') ||
               (b1 == 'G' && b2 == 'C') || (b1 == 'C' && b2 == 'G') ||
               (b1 == 'G' && b2 == 'U') || (b1 == 'U' && b2 == 'G');
    }

    /// <summary>
    /// Checks if a pair is a G-U wobble.
    /// Inputs are case-insensitive and DNA T is treated as RNA U.
    /// </summary>
    public static bool IsWobblePair(char base1, char base2)
    {
        char b1 = NormalizeBase(base1);
        char b2 = NormalizeBase(base2);

        return (b1 == 'G' && b2 == 'U') || (b1 == 'U' && b2 == 'G');
    }

    /// <summary>
    /// Normalizes a nucleotide to the uppercase RNA alphabet (DNA T → RNA U).
    /// </summary>
    private static char NormalizeBase(char b)
    {
        char u = char.ToUpperInvariant(b);
        return u == 'T' ? 'U' : u;
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

        // miRNA pairs antiparallel to the target: miRNA 5'→3' index i pairs with the target
        // base read 3'→5', i.e. target[target.Length-1-i].
        // Source: Lewis BP et al. (2005) Cell 120:15–20 — targets are the reverse complement
        // of the miRNA seed; Watson-Crick (A-U, C-G) + G:U wobble (Crick 1966).
        var alignmentChars = new char[Math.Min(mirna.Length, target.Length)];
        int matches = 0, mismatches = 0, wobbles = 0;

        for (int i = 0; i < alignmentChars.Length; i++)
        {
            char m = mirna[i];
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

    /// <summary>
    /// Estimates the miRNA:target duplex folding free energy (kcal/mol, 37 °C) by summing
    /// Turner 2004 nearest-neighbor stacking energies over consecutive Watson-Crick / G:U
    /// paired positions of the antiparallel duplex. Simplified: only stacking terms are
    /// summed (no loop, bulge, or coaxial terms); see algorithm doc §5.3. All stacking
    /// values are the NNDB Turner 2004 set in <see cref="StackingEnergies"/>.
    /// </summary>
    private static double CalculateDuplexEnergy(string mirna, string target, char[] alignment)
    {
        double energy = 0;

        for (int i = 0; i < alignment.Length - 1; i++)
        {
            bool pairedI = alignment[i] == '|' || alignment[i] == ':';
            bool pairedNext = alignment[i + 1] == '|' || alignment[i + 1] == ':';

            if (!pairedI || !pairedNext)
                continue;

            // Antiparallel nearest-neighbor stack for miRNA positions i, i+1:
            // key = 5'-[m_i][m_{i+1}]-3' / 3'-[t_i][t_{i+1}]-5', where t_i = target[Len-1-i].
            char mi = mirna[i];
            char mn = mirna[i + 1];
            char ti = target[target.Length - 1 - i];
            char tn = target[target.Length - 1 - (i + 1)];

            string key = $"{mi}{mn}/{ti}{tn}";
            if (StackingEnergies.TryGetValue(key, out double stackE))
                energy += stackE;
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

    #region TargetScan context++ scoring (Agarwal et al. 2015) — opt-in

    /// <summary>
    /// Breakdown of the locally-computable TargetScan context++ score (Agarwal et al. 2015)
    /// for a single seed-matched target site. The context++ score is the SUM of per-feature
    /// contributions of a multiple-linear-regression model fit SEPARATELY per site type
    /// (8mer / 7mer-m8 / 7mer-A1 / 6mer). Only the features computable from the miRNA and the
    /// local 3'UTR sequence are realised here; full-transcript features are reported as
    /// <see cref="OmittedFeatures"/>. The <see cref="ContextScorePartial"/> is therefore a
    /// PARTIAL context++ score (an upper bound, since most omitted feature coefficients are
    /// negative), not the published headline CS. See algorithm doc §5.3.
    /// </summary>
    /// <param name="SiteType">The seed-match site type (must be one of the four canonical types).</param>
    /// <param name="Intercept">Site-type intercept term (raw coefficient; Agarwal_2015_parameters.txt).</param>
    /// <param name="LocalAuContribution">coeff(Local_AU) × scaled(local-AU fraction).</param>
    /// <param name="SRna1Contribution">Sum of sRNA1A/C/G indicator contributions (0 if miRNA nt1 = U).</param>
    /// <param name="SRna8Contribution">Sum of sRNA8A/C/G indicator contributions (0 if miRNA nt8 = U).</param>
    /// <param name="Site8Contribution">Sum of Site8A/C/G contributions (only 7mer-A1 / 6mer; else 0).</param>
    /// <param name="ContextScorePartial">Sum of all realised contributions above.</param>
    /// <param name="OmittedFeatures">Full-transcript context++ features NOT included (honest residual).</param>
    public readonly record struct ContextPlusPlusScore(
        TargetSiteType SiteType,
        double Intercept,
        double LocalAuContribution,
        double SRna1Contribution,
        double SRna8Contribution,
        double Site8Contribution,
        double ContextScorePartial,
        IReadOnlyList<string> OmittedFeatures);

    // ── Agarwal et al. (2015) context++ fitted coefficients ──────────────────────────────
    // Source (retrieved verbatim this session): TargetScan distribution parameter file
    //   Agarwal_2015_parameters.txt
    //   https://raw.githubusercontent.com/nsoranzo/targetscan/main/Agarwal_2015_parameters.txt
    // Computation/scaling logic (retrieved verbatim this session): targetscan_70_context_scores.pl
    //   getAgarwalContribution(), getLocalAU_contribution(), get_sRNA1_8_contributions(),
    //   getSite8_contribution()
    //   https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_context_scores.pl
    // Peer-reviewed model: Agarwal V, Bell GW, Nam JW, Bartel DP (2015) eLife 4:e05005,
    //   doi:10.7554/eLife.05005 — "one [regression] model for each of the four site types".
    // Parameter-file column order is 8mer / 7mer-m8 / 7mer-A1 / 6mer.

    // Intercept row (raw coeff, no scaling): siteType2siteTypeContribution in the perl.
    private const double CtxIntercept8mer = -0.589;
    private const double CtxIntercept7merM8 = -0.224;
    private const double CtxIntercept7merA1 = -0.195;
    private const double CtxIntercept6mer = -0.079;

    // Local_AU row: coeff, then min/max for min-max scaling of the local-AU fraction.
    private const double CtxLocalAuCoeff8mer = -0.254, CtxLocalAuMin8mer = 0.308, CtxLocalAuMax8mer = 0.814;
    private const double CtxLocalAuCoeff7merM8 = -0.177, CtxLocalAuMin7merM8 = 0.277, CtxLocalAuMax7merM8 = 0.782;
    private const double CtxLocalAuCoeff7merA1 = -0.075, CtxLocalAuMin7merA1 = 0.342, CtxLocalAuMax7merA1 = 0.801;
    private const double CtxLocalAuCoeff6mer = -0.040, CtxLocalAuMin6mer = 0.295, CtxLocalAuMax6mer = 0.772;

    // sRNA position-1 identity indicators (binary, used raw): coeff per site type.
    private const double CtxSRna1A8mer = -0.018, CtxSRna1A7merM8 = 0.010, CtxSRna1A7merA1 = -0.025, CtxSRna1A6mer = -0.002;
    private const double CtxSRna1C8mer = -0.021, CtxSRna1C7merM8 = 0.014, CtxSRna1C7merA1 = -0.021, CtxSRna1C6mer = 0.004;
    private const double CtxSRna1G8mer = 0.060, CtxSRna1G7merM8 = 0.062, CtxSRna1G7merA1 = 0.030, CtxSRna1G6mer = 0.018;

    // sRNA position-8 identity indicators (binary, used raw): coeff per site type.
    private const double CtxSRna8A8mer = 0.022, CtxSRna8A7merM8 = 0.004, CtxSRna8A7merA1 = -0.049, CtxSRna8A6mer = -0.015;
    private const double CtxSRna8C8mer = 0.012, CtxSRna8C7merM8 = -0.031, CtxSRna8C7merA1 = 0.033, CtxSRna8C6mer = 0.016;
    private const double CtxSRna8G8mer = 0.015, CtxSRna8G7merM8 = -0.008, CtxSRna8G7merA1 = -0.017, CtxSRna8G6mer = 0.006;

    // Site8 identity indicators (binary, used raw): defined ONLY for 7mer-A1 and 6mer
    // (the perl computes them only for siteType 1 = 7mer-A1 and 4 = 6mer).
    private const double CtxSite8A7merA1 = 0.000, CtxSite8A6mer = -0.002;
    private const double CtxSite8C7merA1 = 0.036, CtxSite8C6mer = 0.015;
    private const double CtxSite8G7merA1 = 0.015, CtxSite8G6mer = 0.012;

    // Number of flanking nucleotides used for the local-AU feature (getLocalAU_contribution
    // extracts 30 nt up- and downstream of the site).
    private const int LocalAuFlankLength = 30;

    /// <summary>
    /// Computes the locally-computable TargetScan context++ score (Agarwal et al. 2015) for one
    /// seed-matched target site, as an OPT-IN alternative to the default Grimson-proportional
    /// <c>Score</c> on <see cref="TargetSite"/>. Realises only the four features that depend solely
    /// on the miRNA and the local 3'UTR sequence — site-type intercept, local-AU, sRNA position-1/8
    /// nucleotide identity, and (for 7mer-A1 / 6mer) target site-position-8 identity — each computed
    /// and scaled exactly as in <c>targetscan_70_context_scores.pl</c>. Full-transcript features
    /// (3' supplementary pairing, SPS, TA, minimum distance, structural accessibility, PCT, 3'UTR
    /// and ORF length, ORF-8mer and offset-6mer counts) are NOT computed and are listed in
    /// <see cref="ContextPlusPlusScore.OmittedFeatures"/>; the returned score is therefore a partial
    /// context++ score, not the published headline CS.
    /// </summary>
    /// <param name="mRnaSequence">The mRNA / 3'UTR sequence the site was found in (RNA or DNA; T→U; case-insensitive).</param>
    /// <param name="miRna">The miRNA (its <c>Sequence</c> supplies nt1 and nt8).</param>
    /// <param name="site">A target site whose <c>Type</c> is one of the four canonical seed-match types.</param>
    /// <returns>The per-feature context++ contribution breakdown and their partial sum.</returns>
    /// <exception cref="ArgumentException">Thrown when <paramref name="site"/> is not one of the four canonical seed-match site types (8mer / 7mer-m8 / 7mer-A1 / 6mer).</exception>
    public static ContextPlusPlusScore ScoreTargetSiteContextPlusPlus(
        string mRnaSequence,
        MiRna miRna,
        TargetSite site)
    {
        if (site.Type is not (TargetSiteType.Seed8mer or TargetSiteType.Seed7merM8
            or TargetSiteType.Seed7merA1 or TargetSiteType.Seed6mer))
        {
            throw new ArgumentException(
                "context++ scoring is defined only for the four canonical seed-match site types " +
                "(8mer, 7mer-m8, 7mer-A1, 6mer) — Agarwal et al. (2015).", nameof(site));
        }

        string mrna = string.IsNullOrEmpty(mRnaSequence)
            ? ""
            : mRnaSequence.ToUpperInvariant().Replace('T', 'U');
        string mirna = string.IsNullOrEmpty(miRna.Sequence)
            ? ""
            : miRna.Sequence.ToUpperInvariant().Replace('T', 'U');

        TargetSiteType type = site.Type;

        double intercept = type switch
        {
            TargetSiteType.Seed8mer => CtxIntercept8mer,
            TargetSiteType.Seed7merM8 => CtxIntercept7merM8,
            TargetSiteType.Seed7merA1 => CtxIntercept7merA1,
            _ => CtxIntercept6mer
        };

        double localAu = LocalAuContribution(mrna, site.Start, site.End, type);
        double sRna1 = SRna1Contribution(mirna, type);
        double sRna8 = SRna8Contribution(mirna, type);
        double site8 = Site8Contribution(mrna, site.Start, type);

        double partial = intercept + localAu + sRna1 + sRna8 + site8;

        return new ContextPlusPlusScore(
            SiteType: type,
            Intercept: intercept,
            LocalAuContribution: localAu,
            SRna1Contribution: sRna1,
            SRna8Contribution: sRna8,
            Site8Contribution: site8,
            ContextScorePartial: partial,
            OmittedFeatures: OmittedContextPlusPlusFeatures);
    }

    /// <summary>Full-transcript context++ features that this local implementation does NOT compute.</summary>
    private static readonly IReadOnlyList<string> OmittedContextPlusPlusFeatures = new[]
    {
        "3P_score (3' supplementary pairing)",
        "SPS (predicted seed-pairing stability)",
        "TA_3UTR (target-site abundance)",
        "Min_dist (minimum distance to nearest 3'UTR end)",
        "SA (structural accessibility)",
        "PCT (probability of conserved targeting)",
        "Len_3UTR (3'UTR length)",
        "Len_ORF (ORF length)",
        "ORF8m (ORF 8mer count)",
        "Off6m (offset-6mer count in 3'UTR)"
    };

    // Local_AU: faithful port of getLocalAU_contribution. The local-AU fraction is the
    // position-weighted A/U content of the 30 nt up- and downstream of the site; it is then
    // min-max scaled and multiplied by the site-type Local_AU coefficient. Upstream positions
    // are weighted from the base immediately 5' of the site (i=0) outward; downstream from the
    // base immediately 3' of the site. Weighting offsets are +1 vs +2 per site type, exactly as
    // in the perl (8mer/7mer-m8 favour the upstream by one rank; 8mer/7mer-A1 the downstream).
    private static double LocalAuContribution(string mrna, int siteStart, int siteEnd, TargetSiteType type)
    {
        // Perl uses 1-based utrStart/utrEnd; here Start/End are 0-based inclusive site coordinates.
        // utrUp = up to 30 nt ending at the position immediately before the site (siteStart-1).
        // utrDown = up to 30 nt beginning at the position immediately after the site (siteEnd+1).
        double scoreSum = 0.0;
        double maxRaw = 0.0;

        // Upstream: walk from siteStart-1 backwards (i = 0 at the adjacent base).
        for (int i = 0; i < LocalAuFlankLength; i++)
        {
            int idx = siteStart - 1 - i;
            if (idx < 0) break;
            // 8mer (Seed8mer) and 7mer-m8 use 1/(i+1); 7mer-A1 and 6mer use 1/(i+2).
            double weight = (type is TargetSiteType.Seed8mer or TargetSiteType.Seed7merM8)
                ? 1.0 / (i + 1)
                : 1.0 / (i + 2);
            char b = mrna[idx];
            if (b == 'A' || b == 'U')
                scoreSum += weight;
            maxRaw += weight;
        }

        // Downstream: walk from siteEnd+1 forwards (i = 0 at the adjacent base).
        for (int i = 0; i < LocalAuFlankLength; i++)
        {
            int idx = siteEnd + 1 + i;
            if (idx >= mrna.Length) break;
            // 8mer and 7mer-A1 use 1/(i+2); 7mer-m8 and 6mer use 1/(i+1).
            double weight = (type is TargetSiteType.Seed8mer or TargetSiteType.Seed7merA1)
                ? 1.0 / (i + 2)
                : 1.0 / (i + 1);
            char b = mrna[idx];
            if (b == 'A' || b == 'U')
                scoreSum += weight;
            maxRaw += weight;
        }

        if (maxRaw == 0.0)
            return 0.0;

        double fraction = scoreSum / maxRaw;

        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxLocalAuCoeff8mer, CtxLocalAuMin8mer, CtxLocalAuMax8mer),
            TargetSiteType.Seed7merM8 => (CtxLocalAuCoeff7merM8, CtxLocalAuMin7merM8, CtxLocalAuMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxLocalAuCoeff7merA1, CtxLocalAuMin7merA1, CtxLocalAuMax7merA1),
            _ => (CtxLocalAuCoeff6mer, CtxLocalAuMin6mer, CtxLocalAuMax6mer)
        };

        double scaled = (fraction - min) / (max - min);
        return coeff * scaled;
    }

    // sRNA position-1 indicators: contributions are 0 when miRNA nt1 is U (perl: only computed
    // when sRNA1_nt ne "U"). Otherwise the indicator for the actual nt (A/C/G) is 1, others 0.
    private static double SRna1Contribution(string mirna, TargetSiteType type)
    {
        if (mirna.Length < 1) return 0.0;
        char nt1 = mirna[0];
        if (nt1 == 'U') return 0.0;

        (double a, double c, double g) = type switch
        {
            TargetSiteType.Seed8mer => (CtxSRna1A8mer, CtxSRna1C8mer, CtxSRna1G8mer),
            TargetSiteType.Seed7merM8 => (CtxSRna1A7merM8, CtxSRna1C7merM8, CtxSRna1G7merM8),
            TargetSiteType.Seed7merA1 => (CtxSRna1A7merA1, CtxSRna1C7merA1, CtxSRna1G7merA1),
            _ => (CtxSRna1A6mer, CtxSRna1C6mer, CtxSRna1G6mer)
        };

        return nt1 switch { 'A' => a, 'C' => c, 'G' => g, _ => 0.0 };
    }

    // sRNA position-8 indicators: 0 when miRNA nt8 is U; else indicator for the actual A/C/G.
    private static double SRna8Contribution(string mirna, TargetSiteType type)
    {
        if (mirna.Length < 8) return 0.0;
        char nt8 = mirna[7];
        if (nt8 == 'U') return 0.0;

        (double a, double c, double g) = type switch
        {
            TargetSiteType.Seed8mer => (CtxSRna8A8mer, CtxSRna8C8mer, CtxSRna8G8mer),
            TargetSiteType.Seed7merM8 => (CtxSRna8A7merM8, CtxSRna8C7merM8, CtxSRna8G7merM8),
            TargetSiteType.Seed7merA1 => (CtxSRna8A7merA1, CtxSRna8C7merA1, CtxSRna8G7merA1),
            _ => (CtxSRna8A6mer, CtxSRna8C6mer, CtxSRna8G6mer)
        };

        return nt8 switch { 'A' => a, 'C' => c, 'G' => g, _ => 0.0 };
    }

    // Site8 indicators: only defined for 7mer-A1 and 6mer (the perl computes them only for those
    // site types). The relevant target base is the nucleotide opposite miRNA position 8, i.e. the
    // base immediately 5' of the 6mer core. For 7mer-A1 / 6mer the 6mer core starts at site.Start,
    // so the position-8 base is mrna[site.Start - 1]. 0 when that base is U or out of range.
    private static double Site8Contribution(string mrna, int siteStart, TargetSiteType type)
    {
        if (type is not (TargetSiteType.Seed7merA1 or TargetSiteType.Seed6mer))
            return 0.0;

        int idx = siteStart - 1;
        if (idx < 0 || idx >= mrna.Length) return 0.0;
        char sitePos8 = mrna[idx];
        if (sitePos8 == 'U') return 0.0;

        (double a, double c, double g) = type == TargetSiteType.Seed7merA1
            ? (CtxSite8A7merA1, CtxSite8C7merA1, CtxSite8G7merA1)
            : (CtxSite8A6mer, CtxSite8C6mer, CtxSite8G6mer);

        return sitePos8 switch { 'A' => a, 'C' => c, 'G' => g, _ => 0.0 };
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

    #region Pre-miRNA Prediction — MFE-structure-based (opt-in)

    /// <summary>
    /// A pre-miRNA hairpin candidate assessed from the REAL minimum-free-energy (MFE) secondary
    /// structure produced by the validated Zuker–Stiegler folder
    /// (<see cref="Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure"/>,
    /// RNA-STRUCT-001), rather than the consecutive-pairing heuristic of
    /// <see cref="FindPreMiRnaHairpins"/>.
    /// </summary>
    /// <param name="Start">0-based start of the candidate window in the input sequence.</param>
    /// <param name="End">0-based inclusive end of the candidate window in the input sequence.</param>
    /// <param name="Sequence">The folded candidate (upper-cased; T read as U).</param>
    /// <param name="DotBracket">The folded MFE dot-bracket structure of the candidate.</param>
    /// <param name="StemBasePairs">Number of base pairs in the single dominant hairpin stem.</param>
    /// <param name="TerminalLoopStart">0-based start (within the candidate) of the terminal loop.</param>
    /// <param name="TerminalLoopSize">Length of the single terminal (apical) loop, in nucleotides.</param>
    /// <param name="FreeEnergy">
    /// ΔG° (kcal/mol) of the folded structure — by construction equal to
    /// <see cref="Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMinimumFreeEnergy"/>
    /// for the same candidate (negative for a stable hairpin).
    /// </param>
    /// <param name="Amfe">
    /// Adjusted minimal folding free energy: AMFE = 100·|MFE| / length (kcal/mol per 100 nt).
    /// </param>
    /// <param name="Mfei">
    /// Minimal folding free energy index: MFEI = AMFE / (G+C)% (Zhang et al. 2006).
    /// </param>
    public readonly record struct PreMiRnaMfe(
        int Start,
        int End,
        string Sequence,
        string DotBracket,
        int StemBasePairs,
        int TerminalLoopStart,
        int TerminalLoopSize,
        double FreeEnergy,
        double Amfe,
        double Mfei);

    // Pre-miRNA structural thresholds.
    // Stem: a putative miRNA is embedded in one arm of a fold-back hairpin with >=16 complementary
    // bases to the opposite arm — Ambros et al. (2003), RNA 9:277-279.
    private const int MfeMinStemBasePairs = 16;
    // Terminal (apical) loop: pre-miRNA loops are ~3-15 nt, up to ~25 — Bartel (2004), Cell 116:281.
    private const int MfeMinTerminalLoop = 3;
    private const int MfeMaxTerminalLoop = 25;
    // Per-100-nt normalisation factor for AMFE — Zhang et al. (2006), Cell Mol Life Sci 63:246-254:
    // "AMFE means the MFE of a RNA sequence with 100 nt in length" = MFE/length * 100.
    private const double AmfePer100Nt = 100.0;
    // MFEI discriminative threshold: "MFEI of miRNA precursors is >= 0.85, remarkably higher than
    // other RNAs" — Zhang et al. (2006). Used as the default MFE-fold acceptance cutoff.
    private const double DefaultMinMfei = 0.85;

    /// <summary>
    /// Computes the minimal folding free energy index (MFEI) of Zhang et al. (2006),
    /// Cell Mol Life Sci 63:246-254: MFEI = AMFE / (G+C)% where AMFE = 100·|MFE| / length.
    /// MFE is taken as a magnitude (|ΔG°|) so that the published "MFEI > 0.85" criterion applies
    /// directly to the negative ΔG° returned by the Turner-2004 folder.
    /// </summary>
    /// <param name="freeEnergy">ΔG° (kcal/mol), as returned by the folder (typically negative).</param>
    /// <param name="length">Sequence length in nucleotides (must be &gt; 0).</param>
    /// <param name="gcPercent">G+C content as a percentage in (0,100].</param>
    /// <returns>MFEI (dimensionless). Returns 0 when length or GC% is non-positive.</returns>
    public static double CalculateMfeIndex(double freeEnergy, int length, double gcPercent)
    {
        if (length <= 0 || gcPercent <= 0) return 0.0;
        double amfe = AmfePer100Nt * Math.Abs(freeEnergy) / length;
        return amfe / gcPercent;
    }

    /// <summary>
    /// Finds pre-miRNA hairpins by folding each candidate window with the validated Zuker–Stiegler
    /// MFE engine (<see cref="Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure"/>,
    /// RNA-STRUCT-001) and deriving the hairpin features (single terminal loop, paired-stem count,
    /// ΔG°, AMFE, MFEI) from the ACTUAL MFE dot-bracket structure — an opt-in alternative to the
    /// consecutive-pairing heuristic <see cref="FindPreMiRnaHairpins"/> (which remains the default).
    /// </summary>
    /// <remarks>
    /// Acceptance criteria (all from primary sources):
    /// <list type="bullet">
    /// <item>The MFE structure is a SINGLE dominant hairpin: exactly one contiguous terminal loop
    /// (one maximal run of unpaired bases flanked by paired bases), with no branching.</item>
    /// <item>Stem base pairs ≥ <see cref="MfeMinStemBasePairs"/> (16) — Ambros et al. (2003).</item>
    /// <item>Terminal loop ∈ [<see cref="MfeMinTerminalLoop"/>, <see cref="MfeMaxTerminalLoop"/>]
    /// (3–25 nt) — Bartel (2004).</item>
    /// <item>MFEI ≥ <paramref name="minMfei"/> (default 0.85) — Zhang et al. (2006).</item>
    /// </list>
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence (RNA or DNA; T read as U); case-insensitive.</param>
    /// <param name="minHairpinLength">Minimum candidate window length (default 55 nt).</param>
    /// <param name="maxHairpinLength">Maximum candidate window length (default 120 nt).</param>
    /// <param name="minMfei">Minimum MFEI to accept a candidate (default 0.85, Zhang 2006).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    /// <returns>Accepted pre-miRNA hairpins, each with features from the real MFE structure.</returns>
    public static IEnumerable<PreMiRnaMfe> FindPreMiRnaHairpinsByMfe(
        string sequence,
        int minHairpinLength = 55,
        int maxHairpinLength = 120,
        double minMfei = DefaultMinMfei,
        int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minHairpinLength)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        for (int i = 0; i <= upper.Length - minHairpinLength; i++)
        {
            for (int len = minHairpinLength; len <= Math.Min(maxHairpinLength, upper.Length - i); len++)
            {
                string candidate = upper.Substring(i, len);
                var assessed = AssessHairpinByMfe(candidate, minMfei, minLoopSize);
                if (assessed != null)
                {
                    var a = assessed.Value;
                    yield return new PreMiRnaMfe(
                        Start: i,
                        End: i + len - 1,
                        Sequence: candidate,
                        DotBracket: a.DotBracket,
                        StemBasePairs: a.StemBasePairs,
                        TerminalLoopStart: a.TerminalLoopStart,
                        TerminalLoopSize: a.TerminalLoopSize,
                        FreeEnergy: a.FreeEnergy,
                        Amfe: a.Amfe,
                        Mfei: a.Mfei);
                }
            }
        }
    }

    /// <summary>
    /// Folds a single candidate with the MFE engine and, if its MFE structure is a single dominant
    /// hairpin meeting the pre-miRNA structural thresholds, returns the derived features; otherwise
    /// returns <c>null</c>. Public so callers can assess one known sequence without the window scan.
    /// </summary>
    /// <param name="candidate">Candidate sequence (RNA or DNA; T read as U); case-insensitive.</param>
    /// <param name="minMfei">Minimum MFEI to accept (default 0.85, Zhang 2006).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    public static PreMiRnaMfe? AssessHairpinByMfe(
        string candidate,
        double minMfei = DefaultMinMfei,
        int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(candidate))
            return null;

        var mfe = Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure(candidate, minLoopSize);
        string structure = mfe.DotBracket;
        int n = structure.Length;
        if (n == 0)
            return null;

        // A single dominant hairpin: the MFE structure must have exactly one terminal loop
        // (one maximal run of '.' bounded on the left by '(' and on the right by ')'), all pairs
        // nested with no branching. Detect this directly from the dot-bracket.
        if (!TryDescribeSingleHairpin(structure, out int stemBp, out int loopStart, out int loopSize))
            return null;

        if (stemBp < MfeMinStemBasePairs)
            return null;
        if (loopSize < MfeMinTerminalLoop || loopSize > MfeMaxTerminalLoop)
            return null;

        double gcPercent = GcPercent(mfe.Sequence);
        double amfe = AmfePer100Nt * Math.Abs(mfe.FreeEnergy) / n;
        double mfei = CalculateMfeIndex(mfe.FreeEnergy, n, gcPercent);
        if (mfei < minMfei)
            return null;

        return new PreMiRnaMfe(
            Start: 0,
            End: n - 1,
            Sequence: mfe.Sequence,
            DotBracket: structure,
            StemBasePairs: stemBp,
            TerminalLoopStart: loopStart,
            TerminalLoopSize: loopSize,
            FreeEnergy: mfe.FreeEnergy,
            Amfe: amfe,
            Mfei: mfei);
    }

    /// <summary>
    /// Returns true iff <paramref name="dotBracket"/> describes a SINGLE hairpin: a sequence of
    /// '(' (the 5' stem, possibly interrupted by internal/bulge '.'), exactly one apical run of
    /// '.' (the terminal loop), then a sequence of ')' (the 3' stem). No multiloop / branching:
    /// once a ')' is seen no further '(' may appear. Counts stem base pairs and the terminal loop.
    /// </summary>
    private static bool TryDescribeSingleHairpin(
        string dotBracket, out int stemBasePairs, out int terminalLoopStart, out int terminalLoopSize)
    {
        stemBasePairs = 0;
        terminalLoopStart = -1;
        terminalLoopSize = 0;

        int opens = 0, closes = 0;
        bool seenClose = false;
        // The terminal loop is the run of dots immediately after the last '(' and before the
        // first ')'. Track the index just after the last '(' and the index of the first ')'.
        int lastOpenIndex = -1;
        int firstCloseIndex = -1;

        for (int k = 0; k < dotBracket.Length; k++)
        {
            char c = dotBracket[k];
            if (c == '(')
            {
                if (seenClose) return false; // a '(' after a ')' ⇒ branched / not a single hairpin
                opens++;
                lastOpenIndex = k;
            }
            else if (c == ')')
            {
                if (!seenClose) { seenClose = true; firstCloseIndex = k; }
                closes++;
            }
            else if (c != '.')
            {
                return false; // unexpected symbol
            }
        }

        if (opens == 0 || opens != closes) return false; // need a balanced, non-empty stem
        if (lastOpenIndex < 0 || firstCloseIndex < 0 || firstCloseIndex <= lastOpenIndex) return false;

        // Terminal loop = the dots strictly between the last '(' and the first ')'.
        terminalLoopStart = lastOpenIndex + 1;
        terminalLoopSize = firstCloseIndex - lastOpenIndex - 1;
        // Every position between them must be a dot (a single apical loop, no inner stem).
        for (int k = terminalLoopStart; k < firstCloseIndex; k++)
            if (dotBracket[k] != '.') return false;

        stemBasePairs = opens; // each '(' pairs with one ')'
        return true;
    }

    /// <summary>G+C content of a sequence as a percentage in [0,100].</summary>
    private static double GcPercent(string sequence)
    {
        if (string.IsNullOrEmpty(sequence)) return 0.0;
        int gc = 0;
        foreach (char c in sequence)
        {
            char u = char.ToUpperInvariant(c);
            if (u == 'G' || u == 'C') gc++;
        }
        return AmfePer100Nt * gc / sequence.Length;
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
