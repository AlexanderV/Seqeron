using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Phylogenetics;

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
    /// <param name="SaContribution">coeff(SA) × scaled(log10 of the unpaired probability of the 14-nt window centred on the seed match), from the Turner-2004 McCaskill partition function. 0 (and listed as omitted) when the window does not fit the local 3'UTR context.</param>
    /// <param name="ThreePrimePairingContribution">coeff(3P_score) × scaled(3' supplementary-pairing raw score). Computed from miRNA + 3'UTR.</param>
    /// <param name="MinDistContribution">coeff(Min_dist) × scaled(log10 distance to nearest 3'UTR end). Computed from the 3'UTR.</param>
    /// <param name="Len3UtrContribution">coeff(Len_3UTR) × scaled(log10 3'UTR length). Computed from the 3'UTR.</param>
    /// <param name="Off6mContribution">coeff(Off6m) × offset-6mer count in the 3'UTR (used raw). Computed from miRNA + 3'UTR.</param>
    /// <param name="SpsContribution">coeff(SPS) × scaled(caller-supplied SPS). 0 (and listed as omitted) unless an SPS value is supplied.</param>
    /// <param name="TaContribution">coeff(TA_3UTR) × scaled(caller-supplied TA). 0 (and listed as omitted) unless a TA value is supplied.</param>
    /// <param name="LenOrfContribution">coeff(Len_ORF) × scaled(log10 caller-supplied ORF length). 0 (and listed as omitted) unless an ORF length is supplied.</param>
    /// <param name="Orf8mContribution">coeff(ORF8m) × caller-supplied ORF-8mer count (used raw). 0 (and listed as omitted) unless an ORF-8mer count is supplied.</param>
    /// <param name="PctContribution">coeff(PCT) × scaled(PCT value). PCT is computed from the Friedman et al. (2009) branch-length score (Bls) via the published sigmoid; 0 (and listed as omitted) unless a <see cref="ContextPlusPlusInputs.Conservation"/> input is supplied.</param>
    /// <param name="BranchLengthScore">The Friedman et al. (2009) branch-length score (Bls) computed from the supplied conservation tree, or 0 when no conservation input was supplied.</param>
    /// <param name="Pct">The probability of conserved targeting (PCT) computed from <paramref name="BranchLengthScore"/> via the published sigmoid, or 0 when no conservation input was supplied.</param>
    /// <param name="ContextScorePartial">Sum of all realised contributions above.</param>
    /// <param name="OmittedFeatures">context++ features NOT included for this call (honest residual; depends on which optional inputs were supplied).</param>
    public readonly record struct ContextPlusPlusScore(
        TargetSiteType SiteType,
        double Intercept,
        double LocalAuContribution,
        double SRna1Contribution,
        double SRna8Contribution,
        double Site8Contribution,
        double SaContribution,
        double ThreePrimePairingContribution,
        double MinDistContribution,
        double Len3UtrContribution,
        double Off6mContribution,
        double SpsContribution,
        double TaContribution,
        double LenOrfContribution,
        double Orf8mContribution,
        double PctContribution,
        double BranchLengthScore,
        double Pct,
        double ContextScorePartial,
        IReadOnlyList<string> OmittedFeatures);

    /// <summary>
    /// Optional, caller-supplied context++ feature inputs that the library cannot derive from the
    /// miRNA and 3'UTR sequence alone, but whose contribution is computed FAITHFULLY (verbatim
    /// Agarwal et al. 2015 coefficient × scaling) when the caller provides the value.
    /// All are <c>null</c> by default; a null input means the feature stays an honest residual
    /// (reported in <see cref="ContextPlusPlusScore.OmittedFeatures"/>) and contributes 0.
    /// </summary>
    /// <param name="Sps">Predicted seed-pairing stability (Garcia et al. 2011 table value, kcal/mol; min-max scaled per site type). Data-blocked: requires the Garcia seed-region table.</param>
    /// <param name="Ta">Target-site abundance (transcriptome-wide log10 site count; min-max scaled). Data-blocked: requires a transcriptome.</param>
    /// <param name="OrfLength">ORF (CDS) length in nucleotides (log10-transformed, then min-max scaled). Requires the transcript's ORF.</param>
    /// <param name="Orf8mCount">Count of cognate 8mer sites in the ORF (used raw, no scaling). Requires the transcript's ORF.</param>
    /// <param name="Conservation">Multi-species conservation input (a phylogenetic tree + the set of species in which the site is conserved + the published per-site-type PCT sigmoid parameters). When supplied, the library computes the Friedman et al. (2009) branch-length score (Bls) and the resulting PCT value, then applies the bundled Agarwal et al. (2015) PCT coefficient. <c>null</c> ⇒ PCT stays an honest residual.</param>
    public readonly record struct ContextPlusPlusInputs(
        double? Sps = null,
        double? Ta = null,
        double? OrfLength = null,
        int? Orf8mCount = null,
        PctConservation? Conservation = null);

    // ── PCT (probability of conserved targeting) — Friedman et al. (2009) Genome Res 19:92 ──
    // Bls definition (Methods, retrieved verbatim this session): "The conservation of a given
    // sequence (e.g., an 8mer miR-1 site in a particular 3'UTR) was then assessed by summing the
    // total branch length in the phylogenetic tree connecting the subset of species having the
    // sequence perfectly aligned." (doi:10.1101/gr.082701.108; PMC2612969)
    // PCT definition (Methods): "the P_CT was defined as E[(S − B)/S]" — a signal-to-background
    // correction (targetscan.org/docs/pct.html: PCT ≈ (S/B − 1)/(S/B), or near zero for S/B < 1).
    // TargetScan parameterises the Bls→PCT relationship as a logistic (sigmoid) curve
    // (targetscan_70_BL_PCT.pl, calculatePCTthisBL, retrieved verbatim this session):
    //   $pct = $b0 + ( $b1 / (1 + $eConstant ** ( (0 - $b2) * $BL + $b3))); if ($pct < 0) { $pct = 0; }
    // i.e. PCT(Bls) = b0 + b1 / (1 + e^(−b2·Bls + b3)), truncated at 0.
    // The per-miRNA-family b0..b3 are in TargetScan's COMPILED, citation-required data tables and
    // are NOT published as numbers in Friedman 2009 → they are CALLER-SUPPLIED here (not bundled,
    // not invented). The library bundles only the published EQUATION + the Agarwal PCT coefficient.

    /// <summary>
    /// The four published logistic (sigmoid) parameters that map a branch-length score (Bls) to a
    /// PCT value for one miRNA family + site type, per <c>targetscan_70_BL_PCT.pl</c>
    /// (<c>calculatePCTthisBL</c>): <c>PCT(Bls) = B0 + B1 / (1 + e^(−B2·Bls + B3))</c>, truncated at 0.
    /// </summary>
    /// <remarks>
    /// These coefficients are fitted per miRNA family in TargetScan's compiled, citation-required
    /// data tables (8mer/7mer-m8/7mer-1a <c>PCT_parameters.txt</c>); Friedman et al. (2009) does not
    /// publish them as numbers. They are therefore caller-supplied (not bundled with the library),
    /// matching how SPS / TA / ORF inputs are caller-supplied.
    /// </remarks>
    /// <param name="B0">Offset term (b0).</param>
    /// <param name="B1">Scale term (b1).</param>
    /// <param name="B2">Bls slope inside the logistic (b2).</param>
    /// <param name="B3">Logistic offset (b3).</param>
    public readonly record struct PctSigmoidParameters(double B0, double B1, double B2, double B3);

    /// <summary>
    /// Multi-species conservation evidence for one target site: the phylogenetic tree relating the
    /// aligned species, the set of species in which the seed match is conserved (perfectly aligned),
    /// and the published per-site-type PCT sigmoid parameters. The library derives the Friedman
    /// et al. (2009) branch-length score (Bls) from <see cref="Tree"/> + <see cref="SpeciesWithSite"/>
    /// and the PCT value from <see cref="SigmoidParameters"/>.
    /// </summary>
    /// <param name="Tree">A phylogenetic tree (<see cref="PhylogeneticAnalyzer.PhyloNode"/>; e.g. parsed from Newick) whose LEAF names are species identifiers and whose <c>BranchLength</c> values are the edge lengths.</param>
    /// <param name="SpeciesWithSite">The species (leaf names) in which the seed match is conserved (perfectly aligned). Bls = total branch length of the minimal subtree connecting these species (Friedman 2009).</param>
    /// <param name="SigmoidParameters">The caller-supplied published per-site-type b0..b3 (see <see cref="PctSigmoidParameters"/>).</param>
    public readonly record struct PctConservation(
        PhylogeneticAnalyzer.PhyloNode Tree,
        IReadOnlyCollection<string> SpeciesWithSite,
        PctSigmoidParameters SigmoidParameters);

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

    // 3P_score row (3' supplementary pairing): coeff, then min/max for min-max scaling of the raw
    // pairing score. Agarwal_2015_parameters.txt: min=1, max=3.5 for all four site types.
    private const double Ctx3PCoeff8mer = -0.040, Ctx3PCoeff7merM8 = -0.055, Ctx3PCoeff7merA1 = -0.060, Ctx3PCoeff6mer = -0.024;
    private const double Ctx3PMin = 1.0, Ctx3PMax = 3.5;

    // Min_dist row: coeff, then min/max for min-max scaling of log10(distance to nearest 3'UTR end).
    private const double CtxMinDistCoeff8mer = 0.118, CtxMinDistCoeff7merM8 = 0.056, CtxMinDistCoeff7merA1 = 0.045, CtxMinDistCoeff6mer = 0.036;
    private const double CtxMinDistMin8mer = 1.415, CtxMinDistMax8mer = 3.113;
    private const double CtxMinDistMin7merM8 = 1.491, CtxMinDistMax7merM8 = 3.096;
    private const double CtxMinDistMin7merA1 = 1.431, CtxMinDistMax7merA1 = 3.117;
    private const double CtxMinDistMin6mer = 1.477, CtxMinDistMax6mer = 3.106;

    // Len_3UTR row: coeff, then min/max for min-max scaling of log10(3'UTR length).
    private const double CtxLen3UtrCoeff8mer = 0.310, CtxLen3UtrCoeff7merM8 = 0.154, CtxLen3UtrCoeff7merA1 = 0.129, CtxLen3UtrCoeff6mer = 0.045;
    private const double CtxLen3UtrMin8mer = 2.392, CtxLen3UtrMax8mer = 3.637;
    private const double CtxLen3UtrMin7merM8 = 2.409, CtxLen3UtrMax7merM8 = 3.615;
    private const double CtxLen3UtrMin7merA1 = 2.413, CtxLen3UtrMax7merA1 = 3.630;
    private const double CtxLen3UtrMin6mer = 2.405, CtxLen3UtrMax6mer = 3.620;

    // Off6m row (offset-6mer count): coeff only — used raw (NOT in the min-max regex of getAgarwalContribution).
    private const double CtxOff6mCoeff8mer = -0.020, CtxOff6mCoeff7merM8 = -0.011, CtxOff6mCoeff7merA1 = -0.020, CtxOff6mCoeff6mer = -0.010;

    // SPS row (seed-pairing stability; caller-supplied): coeff, then min/max for min-max scaling.
    private const double CtxSpsCoeff8mer = 0.210, CtxSpsCoeff7merM8 = 0.135, CtxSpsCoeff7merA1 = 0.095, CtxSpsCoeff6mer = 0.035;
    private const double CtxSpsMin8mer = -11.13, CtxSpsMax8mer = -5.52;
    private const double CtxSpsMin7merM8 = -11.13, CtxSpsMax7merM8 = -5.49;
    private const double CtxSpsMin7merA1 = -8.41, CtxSpsMax7merA1 = -3.33;
    private const double CtxSpsMin6mer = -8.57, CtxSpsMax6mer = -3.33;

    // TA_3UTR row (target-site abundance; caller-supplied): coeff, then min/max for min-max scaling.
    private const double CtxTaCoeff8mer = 0.222, CtxTaCoeff7merM8 = 0.139, CtxTaCoeff7merA1 = 0.117, CtxTaCoeff6mer = 0.058;
    private const double CtxTaMin8mer = 3.113, CtxTaMax8mer = 3.865;
    private const double CtxTaMin7merM8 = 3.067, CtxTaMax7merM8 = 3.887;
    private const double CtxTaMin7merA1 = 3.145, CtxTaMax7merA1 = 3.887;
    private const double CtxTaMin6mer = 3.113, CtxTaMax6mer = 3.887;

    // Len_ORF row (ORF length; caller-supplied): coeff, then min/max for min-max scaling of log10(ORF length).
    private const double CtxLenOrfCoeff8mer = 0.205, CtxLenOrfCoeff7merM8 = 0.100, CtxLenOrfCoeff7merA1 = 0.063, CtxLenOrfCoeff6mer = 0.029;
    private const double CtxLenOrfMin8mer = 2.788, CtxLenOrfMax8mer = 3.753;
    private const double CtxLenOrfMin7merM8 = 2.773, CtxLenOrfMax7merM8 = 3.729;
    private const double CtxLenOrfMin7merA1 = 2.773, CtxLenOrfMax7merA1 = 3.730;
    private const double CtxLenOrfMin6mer = 2.775, CtxLenOrfMax6mer = 3.731;

    // ORF8m row (ORF 8mer count; caller-supplied): coeff only — used raw (NOT min-max scaled).
    private const double CtxOrf8mCoeff8mer = -0.118, CtxOrf8mCoeff7merM8 = -0.044, CtxOrf8mCoeff7merA1 = -0.058, CtxOrf8mCoeff6mer = -0.060;

    // PCT row (probability of conserved targeting; Friedman et al. 2009): coeff, then min/max for
    // min-max scaling of the PCT value. PCT enters getAgarwalContribution exactly like the other
    // scaled features (targetscan_70_context_scores.pl getPCT_contribution → getAgarwalContribution,
    // PCT matches the min-max regex branch). Verbatim from Agarwal_2015_parameters.txt (PCT row),
    // column order 8mer / 7mer-m8 / 7mer-A1 / 6mer (coeff×4, min×4, max×4):
    //   PCT  -0.103  -0.048  -0.048  0.005  0  0  0  0  0.816  0.364  0.449  0.193
    private const double CtxPctCoeff8mer = -0.103, CtxPctMin8mer = 0.0, CtxPctMax8mer = 0.816;
    private const double CtxPctCoeff7merM8 = -0.048, CtxPctMin7merM8 = 0.0, CtxPctMax7merM8 = 0.364;
    private const double CtxPctCoeff7merA1 = -0.048, CtxPctMin7merA1 = 0.0, CtxPctMax7merA1 = 0.449;
    private const double CtxPctCoeff6mer = 0.005, CtxPctMin6mer = 0.0, CtxPctMax6mer = 0.193;

    // SA row (structural accessibility): coeff, then min/max for min-max scaling of log10(plfold).
    // Verbatim from Agarwal_2015_parameters.txt (SA row), column order 8mer / 7mer-m8 / 7mer-A1 / 6mer:
    //   SA  -0.115  -0.134  -0.077  -0.028  -4.356  -5.218  -4.23  -5.082  -0.661  -0.725  -0.588  -0.666
    private const double CtxSaCoeff8mer = -0.115, CtxSaMin8mer = -4.356, CtxSaMax8mer = -0.661;
    private const double CtxSaCoeff7merM8 = -0.134, CtxSaMin7merM8 = -5.218, CtxSaMax7merM8 = -0.725;
    private const double CtxSaCoeff7merA1 = -0.077, CtxSaMin7merA1 = -4.23, CtxSaMax7merA1 = -0.588;
    private const double CtxSaCoeff6mer = -0.028, CtxSaMin6mer = -5.082, CtxSaMax6mer = -0.666;

    // SA structural-accessibility window parameters (targetscan_70_context_scores.pl).
    //   getSA_contribution reads column index 13 of the RNAplfold _lunp row 7 nt downstream of the
    //   seed-match start → the 14th unpaired-probability column → the probability that the 14-nt
    //   stretch ENDING at that row is unpaired (RNAplfold man page; Bernhart et al. 2006). The
    //   eLife paper (Agarwal 2015, Fig 4A): "the log10 value of the unpaired probability for a
    //   14-nt region centered on the match to miRNA nucleotides 7 and 8".
    private const int SaUnpairedWindowLength = 14;   // RNAplfold -u column 14 used by getSA
    private const int SaRowOffsetDownstream = 7;      // grep -A7 → the row 7 nt 3' of utrStart
    // RNAplfold local-fold parameters used by runRNAplfold_all_UTRs: "RNAplfold -L 40 -W 80 -u 20".
    private const int SaPlfoldWindowSize = 80;        // -W (sliding-window averaging size)
    private const int SaPlfoldMaxSpan = 40;           // -L (maximum base-pair span)

    // Number of flanking nucleotides used for the local-AU feature (getLocalAU_contribution
    // extracts 30 nt up- and downstream of the site).
    private const int LocalAuFlankLength = 30;

    // 3' supplementary-pairing alignment constants (targetscan_70_context_scores.pl):
    //   $DESIRED_UTR_ALIGNMENT_LENGTH = 23; site-type-specific UTR start offset is 16 nt
    //   upstream of the seed match.
    private const int DesiredUtrAlignmentLength = 23;
    private const int ThreePrimeUtrStartOffset = 16;

    /// <summary>
    /// Computes the locally-computable TargetScan context++ score (Agarwal et al. 2015) for one
    /// seed-matched target site, as an OPT-IN alternative to the default Grimson-proportional
    /// <c>Score</c> on <see cref="TargetSite"/>. Realises only the four features that depend solely
    /// on the miRNA and the local 3'UTR sequence — site-type intercept, local-AU, sRNA position-1/8
    /// nucleotide identity, and (for 7mer-A1 / 6mer) target site-position-8 identity — each computed
    /// and scaled exactly as in <c>targetscan_70_context_scores.pl</c>. Full-transcript features
    /// 3' supplementary pairing (3P_score), minimum distance to the nearest 3'UTR end (Min_dist),
    /// 3'UTR length (Len_3UTR), and offset-6mer count (Off6m) — each computed and scaled exactly as
    /// in <c>targetscan_70_context_scores.pl</c>. Features that depend on data the library cannot
    /// derive from the miRNA and 3'UTR (SPS, TA, ORF length, ORF-8mer count) are computed faithfully
    /// only when the caller supplies them via <paramref name="inputs"/>; otherwise they contribute 0
    /// and are reported in <see cref="ContextPlusPlusScore.OmittedFeatures"/>. Structural
    /// accessibility (SA) is now computed from the Turner-2004 McCaskill partition function
    /// (<see cref="RnaSecondaryStructure.CalculateRegionUnpairedProbability"/>) — the 14-nt-window
    /// unpaired probability, log10-transformed and min-max scaled exactly as in
    /// <c>getSA_contribution</c> — and is omitted only when that window does not fit the local
    /// 3'UTR context. Conservation (PCT) is computed when the caller supplies a
    /// <see cref="ContextPlusPlusInputs.Conservation"/> input (a tree + the species in which the
    /// site is conserved + the published per-site-type sigmoid parameters): the library derives the
    /// Friedman et al. (2009) branch-length score (Bls), maps it to a PCT value via the published
    /// sigmoid, and applies the bundled Agarwal PCT coefficient. Otherwise PCT stays an honest
    /// residual. The returned score is therefore still a PARTIAL context++ score, not the published
    /// headline CS.
    /// </summary>
    /// <param name="mRnaSequence">The mRNA / 3'UTR sequence the site was found in (RNA or DNA; T→U; case-insensitive). Treated as the 3'UTR for the length / Min_dist / Off6m features.</param>
    /// <param name="miRna">The miRNA (its <c>Sequence</c> supplies nt1, nt8, the seed and the 3' supplementary region).</param>
    /// <param name="site">A target site whose <c>Type</c> is one of the four canonical seed-match types.</param>
    /// <param name="inputs">Optional caller-supplied feature values (SPS / TA / ORF length / ORF-8mer count) that cannot be derived from sequence alone.</param>
    /// <returns>The per-feature context++ contribution breakdown and their partial sum.</returns>
    /// <exception cref="ArgumentException">Thrown when <paramref name="site"/> is not one of the four canonical seed-match site types (8mer / 7mer-m8 / 7mer-A1 / 6mer).</exception>
    public static ContextPlusPlusScore ScoreTargetSiteContextPlusPlus(
        string mRnaSequence,
        MiRna miRna,
        TargetSite site,
        ContextPlusPlusInputs inputs = default)
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
        double sa = SaContribution(mrna, site.Start, type, out bool saIncluded);
        double threeP = ThreePrimePairingContribution(mrna, mirna, site.Start, site.End, type);
        double minDist = MinDistContribution(mrna, site.Start, site.End, type);
        double len3Utr = Len3UtrContribution(mrna, type);
        double off6m = Off6mContribution(mrna, mirna, type);

        double sps = inputs.Sps is double spsRaw ? SpsContribution(spsRaw, type) : 0.0;
        double ta = inputs.Ta is double taRaw ? TaContribution(taRaw, type) : 0.0;
        double lenOrf = inputs.OrfLength is double orfLen ? LenOrfContribution(orfLen, type) : 0.0;
        double orf8m = inputs.Orf8mCount is int orf8mCount ? Orf8mContribution(orf8mCount, type) : 0.0;

        double bls = 0.0, pct = 0.0, pctContribution = 0.0;
        if (inputs.Conservation is PctConservation conservation)
        {
            bls = ComputeBranchLengthScore(conservation.Tree, conservation.SpeciesWithSite);
            pct = PctFromBranchLength(bls, conservation.SigmoidParameters);
            pctContribution = PctContribution(pct, type);
        }

        double partial = intercept + localAu + sRna1 + sRna8 + site8 + sa
            + threeP + minDist + len3Utr + off6m
            + sps + ta + lenOrf + orf8m + pctContribution;

        return new ContextPlusPlusScore(
            SiteType: type,
            Intercept: intercept,
            LocalAuContribution: localAu,
            SRna1Contribution: sRna1,
            SRna8Contribution: sRna8,
            Site8Contribution: site8,
            SaContribution: sa,
            ThreePrimePairingContribution: threeP,
            MinDistContribution: minDist,
            Len3UtrContribution: len3Utr,
            Off6mContribution: off6m,
            SpsContribution: sps,
            TaContribution: ta,
            LenOrfContribution: lenOrf,
            Orf8mContribution: orf8m,
            PctContribution: pctContribution,
            BranchLengthScore: bls,
            Pct: pct,
            ContextScorePartial: partial,
            OmittedFeatures: BuildOmittedFeatures(inputs, saIncluded));
    }

    /// <summary>
    /// context++ features that remain residual for a given call. SA (structural accessibility) is
    /// now computed from the Turner-2004 McCaskill partition function and is omitted only when the
    /// 14-nt window does not fit the local 3'UTR context. PCT (conservation) is always omitted
    /// (needs a multi-species alignment). SPS, TA, ORF length and ORF-8mer count are omitted only
    /// when the caller did not supply them.
    /// </summary>
    private static IReadOnlyList<string> BuildOmittedFeatures(ContextPlusPlusInputs inputs, bool saIncluded)
    {
        var omitted = new List<string>();
        if (!saIncluded)
            omitted.Add("SA (structural accessibility — window does not fit the local 3'UTR context)");
        if (inputs.Conservation is null)
            omitted.Add("PCT (probability of conserved targeting — requires a multi-species alignment + tree)");
        if (inputs.Sps is null) omitted.Add("SPS (predicted seed-pairing stability — requires Garcia et al. 2011 seed-region table)");
        if (inputs.Ta is null) omitted.Add("TA_3UTR (target-site abundance — requires a transcriptome)");
        if (inputs.OrfLength is null) omitted.Add("Len_ORF (ORF length — requires the transcript ORF)");
        if (inputs.Orf8mCount is null) omitted.Add("ORF8m (ORF 8mer count — requires the transcript ORF)");
        return omitted;
    }

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

    // ── SA: structural accessibility (getSA_contribution) ─────────────────────────────────
    // Faithful port of getSA_contribution. TargetScan runs RNAplfold (-L 40 -W 80 -u 20) on each
    // UTR to produce a _lunp file, then for a site reads the probability that the 14-nt stretch
    // ENDING 7 nt downstream of the seed-match start is unpaired (column index 13 of `cut -f 2-`,
    // i.e. the 14th unpaired-probability column → L = 14), takes log10, and min-max scales it by
    // the SA coefficient. Here the unpaired probability is computed EXACTLY from the Turner-2004
    // McCaskill partition function (RnaSecondaryStructure.CalculateRegionUnpairedProbability),
    // over a local fold window of up to W = 80 nt centred on the SA window (max base-pair span
    // L = 40), mirroring RNAplfold's local-folding intent. When the 14-nt window cannot be placed
    // within the 3'UTR (too close to an end), SA is reported as omitted and contributes 0 —
    // matching the perl's "missing _lunp row → plfold = 0 → return 0".
    // Sources: targetscan_70_context_scores.pl getSA_contribution / runRNAplfold_all_UTRs;
    // RNAplfold man page + Bernhart et al. (2006) Bioinformatics 22:614; Agarwal et al. (2015)
    // eLife 4:e05005 Fig 4A ("log10 value of the unpaired probability for a 14-nt region centered
    // on the match to miRNA nucleotides 7 and 8").
    private static double SaContribution(string mrna, int siteStart, TargetSiteType type, out bool included)
    {
        included = false;

        // 1-based UTR start of the seed match; perl decrements it for 7mer-A1 (siteType 1) so the
        // positioning matches 8mer / 7mer-m8 sites (the legacy "|| siteType == 5" is dead code).
        int utrStart1 = siteStart + 1;
        if (type == TargetSiteType.Seed7merA1)
            utrStart1--;

        // The _lunp row read is 7 nt downstream of utrStart; that row's L=14 column is the
        // probability that the 14-nt stretch ENDING at that row is unpaired. Convert to a 0-based
        // window-end index.
        int windowEnd0 = (utrStart1 + SaRowOffsetDownstream) - 1;
        int windowStart0 = windowEnd0 - SaUnpairedWindowLength + 1;
        if (windowStart0 < 0 || windowEnd0 >= mrna.Length)
            return 0.0; // window does not fit → SA omitted (perl: plfold missing → 0)

        // Local fold context: up to W = 80 nt centred on the window, clamped to the UTR. RNAplfold
        // averages over length-W windows; folding this local context captures the local
        // accessibility of the 14-nt window (a base can only pair within ±L = 40 nt).
        int contextStart = Math.Max(0, windowEnd0 - (SaPlfoldWindowSize - SaUnpairedWindowLength) / 2 - SaUnpairedWindowLength + 1);
        contextStart = Math.Min(contextStart, windowStart0);
        int contextEnd = Math.Min(mrna.Length - 1, contextStart + SaPlfoldWindowSize - 1);
        contextStart = Math.Max(0, contextEnd - SaPlfoldWindowSize + 1);
        string context = mrna.Substring(contextStart, contextEnd - contextStart + 1);
        int localWindowEnd = windowEnd0 - contextStart;

        double plfold = RnaSecondaryStructure.CalculateRegionUnpairedProbability(
            context, localWindowEnd, SaUnpairedWindowLength);

        // log10(plfold); perl returns 0 when plfold is not a nonzero number (isNonzeroNumber).
        double log10Plfold = (plfold > 0) ? Math.Log10(plfold) : 0.0;

        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxSaCoeff8mer, CtxSaMin8mer, CtxSaMax8mer),
            TargetSiteType.Seed7merM8 => (CtxSaCoeff7merM8, CtxSaMin7merM8, CtxSaMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxSaCoeff7merA1, CtxSaMin7merA1, CtxSaMax7merA1),
            _ => (CtxSaCoeff6mer, CtxSaMin6mer, CtxSaMax6mer)
        };

        included = true;
        return ScaledContribution(log10Plfold, coeff, min, max);
    }

    // ── Generic min-max scaling (getAgarwalContribution) ──────────────────────────────────
    // Perl: scaledScore = (raw - min) / (max - min); contribution = coeff × scaledScore.
    // NOT clamped to [0,1] (scaled values < 0 or > 1 are kept, exactly as in the perl).
    private static double ScaledContribution(double raw, double coeff, double min, double max)
        => max == min ? coeff * raw : coeff * ((raw - min) / (max - min));

    // ── 3' supplementary pairing (3P_score) ───────────────────────────────────────────────
    // Faithful port of get3primePairingContribution, fed by extractSubseqForAlignment +
    // modifySubseqForAlignment (targetscan_70_context_scores.pl). The raw score is the maximum,
    // over all single-gap offsets in both the UTR and the miRNA, of: 1.0 per contiguous-stretch
    // base pair in the 3' "seed-supplementary" window (offset-adjusted positions 4..7) and 0.5
    // elsewhere, summed only over runs of ≥2 consecutive pairs, minus max(0,(offset-2)/2). The
    // raw score is then min-max scaled (min=1, max=3.5) and multiplied by the 3P_score coeff.
    private static readonly int[] ThreePrimeUtrStart  = { 0, 8, 8, 8, 8 }; // index by perl siteType 1..4
    private static readonly int[] ThreePrimeMiRnaStart = { 0, 7, 8, 8, 8 };
    private static readonly int[] ThreePrimeOverhang   = { 0, 1, 0, 0, 0 };

    private static double ThreePrimePairingContribution(string mrna, string mirna, int siteStart, int siteEnd, TargetSiteType type)
    {
        double raw = ThreePrimePairingRaw(mrna, mirna, siteStart, siteEnd, type);
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (Ctx3PCoeff8mer, Ctx3PMin, Ctx3PMax),
            TargetSiteType.Seed7merM8 => (Ctx3PCoeff7merM8, Ctx3PMin, Ctx3PMax),
            TargetSiteType.Seed7merA1 => (Ctx3PCoeff7merA1, Ctx3PMin, Ctx3PMax),
            _ => (Ctx3PCoeff6mer, Ctx3PMin, Ctx3PMax)
        };
        return ScaledContribution(raw, coeff, min, max);
    }

    // Perl site-type integers: 1 = 7mer-A1, 2 = 7mer-m8, 3 = 8mer, 4 = 6mer.
    private static int PerlSiteType(TargetSiteType type) => type switch
    {
        TargetSiteType.Seed7merA1 => 1,
        TargetSiteType.Seed7merM8 => 2,
        TargetSiteType.Seed8mer => 3,
        _ => 4
    };

    // Returns the raw 3P pairing score (the value getAgarwalContribution receives for "3P_score").
    private static double ThreePrimePairingRaw(string mrna, string mirna, int siteStart, int siteEnd, TargetSiteType type)
    {
        if (mrna.Length == 0 || mirna.Length == 0) return 0.0;

        int perlType = PerlSiteType(type);
        // C# 0-based inclusive Start/End ↔ perl 1-based utrStart/utrEnd.
        int utrStart = siteStart + 1;
        int utrEnd = siteEnd + 1;

        string subseq = ExtractSubseqForAlignment(mrna, utrStart, utrEnd, perlType);
        (string finalUtr, string mirnaForAlign) = ModifySubseqForAlignment(mirna, subseq, perlType);
        return ThreePrimePairingScore(perlType, finalUtr, mirnaForAlign);
    }

    // Port of extractSubseqForAlignment (siteType < 5 path only; the four canonical types).
    private static string ExtractSubseqForAlignment(string utr, int utrStart, int utrEnd, int perlType)
    {
        int realStart = utrStart - ThreePrimeUtrStartOffset;
        if (realStart < 0) realStart = 0;

        // real_end = utrEnd + 1 for siteType 1, 2, 4; else utrEnd.
        int realEnd = (perlType == 1 || perlType == 2 || perlType == 4) ? utrEnd + 1 : utrEnd;
        if (realStart >= realEnd) realStart = 0;

        int length = realEnd - realStart;
        // substr in 0-based: perl coordinates are 1-based, so the 0-based offset is realStart-1...
        // but the perl applies substr($seq, $real_start, $length) with $real_start already a
        // sequence offset (it is utrStart-16 where utrStart is 1-based, so it doubles as 0-based
        // start of the window). We mirror the perl arithmetic exactly using 0-based realStart.
        int begin = realStart;
        if (begin < 0) begin = 0;
        int avail = Math.Max(0, utr.Length - begin);
        int take = Math.Min(length, avail);
        string sub = take > 0 ? utr.Substring(begin, take) : "";

        // Pad leading N's up to DESIRED_UTR_ALIGNMENT_LENGTH (perl loop adds 23-len+1 N's).
        if (DesiredUtrAlignmentLength > sub.Length)
        {
            int nCount = DesiredUtrAlignmentLength - sub.Length + 1;
            sub = new string('N', nCount) + sub;
        }

        // Site at the very end of the UTR: perl appends trailing N's and converts a leading "NN".
        if (utr.Length < (begin + length))
        {
            int lengthDiff = begin + length - utr.Length;
            for (int i = 0; i < lengthDiff; i++)
            {
                sub += "N";
                if (sub.StartsWith("NN", StringComparison.Ordinal))
                    sub = "  " + sub.Substring(2);
            }
        }
        return sub;
    }

    // Port of modifySubseqForAlignment.
    private static (string finalUtr, string mirnaForAlign) ModifySubseqForAlignment(string matureMiRna, string subseq, int perlType)
    {
        int spacer1Length = matureMiRna.Length - DesiredUtrAlignmentLength;
        string spacer1 = spacer1Length > 0 ? new string(' ', spacer1Length) : "";

        int spacer2Length = spacer1Length < 0 ? -spacer1Length : 0;
        var spacer2 = new StringBuilder(new string(' ', spacer2Length));

        if (DesiredUtrAlignmentLength > subseq.Length)
        {
            int extra = DesiredUtrAlignmentLength - subseq.Length + 1;
            spacer2.Append(' ', extra);
        }

        // 7mer-A1 (1) and 8mer (3): chop one char off spacer2.
        if ((perlType == 1 || perlType == 3) && spacer2.Length > 0)
            spacer2.Length -= 1;

        string finalUtr = spacer1 + subseq;
        char[] rev = matureMiRna.ToCharArray();
        Array.Reverse(rev);
        string mirnaForAlign = spacer2.ToString() + new string(rev);
        return (finalUtr, mirnaForAlign);
    }

    // Port of get3primePairingContribution — returns the raw best score only.
    private static double ThreePrimePairingScore(int perlType, string utrIn, string mirnaIn)
    {
        // Normalise: strip spaces/newlines, uppercase, T→U.
        string utr = NormaliseForPairing(utrIn);
        string mirna = NormaliseForPairing(mirnaIn);

        // Reverse both (perl: $utr = reverse $utr; $mirna = reverse $mirna).
        utr = Reverse(utr);
        mirna = Reverse(mirna);

        int utrStart = ThreePrimeUtrStart[perlType];
        int mirnaStart = ThreePrimeMiRnaStart[perlType];
        int overhang = ThreePrimeOverhang[perlType];

        string utrNum = utrStart < utr.Length ? utr.Substring(utrStart) : "";
        string mirnaNum = mirnaStart < mirna.Length ? mirna.Substring(mirnaStart) : "";

        int maxScoreLen = Math.Max(utrNum.Length, mirnaNum.Length);

        int[] UTR = ToBaseCodes(utrNum);
        int[] MIRNA = ToBaseCodes(mirnaNum);

        double bestScore = double.NegativeInfinity;

        for (int offset = 0; offset < maxScoreLen; offset++)
        {
            // "top" gap orientation: MIRNA shifted by offset.
            double topScore = PairingRunScore(UTR, MIRNA, offset, overhang, gapOnTop: true);
            topScore -= Math.Max(0.0, (offset - 2) / 2.0);
            if (topScore > bestScore) bestScore = topScore;

            // "bottom" gap orientation: UTR shifted by offset.
            double bottomScore = PairingRunScore(UTR, MIRNA, offset, overhang, gapOnTop: false);
            bottomScore -= Math.Max(0.0, (offset - 2) / 2.0);
            if (bottomScore > bestScore) bestScore = bottomScore;
        }

        return double.IsNegativeInfinity(bestScore) ? 0.0 : bestScore;
    }

    // One alignment offset: scores only runs of ≥2 consecutive base pairs, taking the best such run.
    private static double PairingRunScore(int[] UTR, int[] MIRNA, int offset, int overhang, bool gapOnTop)
    {
        double score = 0.0;     // best committed run score for this offset
        double tempscore = 0.0; // current run score
        int prevmatch = 0;      // current run length

        int limit = gapOnTop
            ? Math.Min(MIRNA.Length - 1 - offset, UTR.Length - 1)
            : Math.Min(UTR.Length - 1 - offset, MIRNA.Length - 1);

        for (int i = 0; i <= limit; i++)
        {
            int u = gapOnTop ? UTR[i] : UTR[i + offset];
            int m = gapOnTop ? MIRNA[i + offset] : MIRNA[i];
            // Watson-Crick A:U/U:A → product 1×2 or 2×1 = 2; G:C/C:G → 3×4 or 4×3 = 12.
            bool isPair = (u * m == 2) || (u * m == 12);

            if (isPair)
            {
                int posCheck = gapOnTop ? (i + offset - overhang) : (i - overhang);
                if (prevmatch == 0) tempscore = 0.0;
                tempscore += (posCheck >= 4 && posCheck <= 7) ? 1.0 : 0.5;
                prevmatch++;
            }
            else
            {
                if (prevmatch >= 2 && tempscore > score) score = tempscore;
                tempscore = 0.0;
                prevmatch = 0;
            }
        }
        if (prevmatch >= 2 && tempscore > score) score = tempscore;
        return score;
    }

    private static string NormaliseForPairing(string s)
    {
        var sb = new StringBuilder(s.Length);
        foreach (char ch in s)
        {
            if (ch == ' ' || ch == '\n' || ch == '\r') continue;
            char u = char.ToUpperInvariant(ch);
            sb.Append(u == 'T' ? 'U' : u);
        }
        return sb.ToString();
    }

    private static string Reverse(string s)
    {
        char[] c = s.ToCharArray();
        Array.Reverse(c);
        return new string(c);
    }

    // tr/AUCGN/12345/ — any other character maps to 0 (cannot pair).
    private static int[] ToBaseCodes(string s)
    {
        var codes = new int[s.Length];
        for (int i = 0; i < s.Length; i++)
            codes[i] = s[i] switch { 'A' => 1, 'U' => 2, 'C' => 3, 'G' => 4, 'N' => 5, _ => 0 };
        return codes;
    }

    // ── Min_dist ───────────────────────────────────────────────────────────────────────────
    // Port of getMinDist_weighted_contribution for the single-isoform case (AIR_end = UTR length).
    // distTo5 = siteStart-1 (perl 1-based) ; distTo3 = UTRlen - siteEnd ; take the smaller, log10,
    // then min-max scale by the site-type Min_dist coefficients.
    private static double MinDistContribution(string mrna, int siteStart, int siteEnd, TargetSiteType type)
    {
        if (mrna.Length == 0) return 0.0;
        // perl 1-based site coordinates and UTR length (= AIR_end).
        int perlStart = siteStart + 1;
        int perlEnd = siteEnd + 1;
        int airEnd = mrna.Length;
        if (airEnd < perlEnd) return 0.0; // UTR shorter than the site (defensive)

        int distTo5 = perlStart - 1;
        int distTo3 = airEnd - perlEnd;
        int nearest = Math.Min(distTo5, distTo3);

        double log10 = nearest > 0 ? Math.Log10(nearest) : 0.0;

        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxMinDistCoeff8mer, CtxMinDistMin8mer, CtxMinDistMax8mer),
            TargetSiteType.Seed7merM8 => (CtxMinDistCoeff7merM8, CtxMinDistMin7merM8, CtxMinDistMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxMinDistCoeff7merA1, CtxMinDistMin7merA1, CtxMinDistMax7merA1),
            _ => (CtxMinDistCoeff6mer, CtxMinDistMin6mer, CtxMinDistMax6mer)
        };
        return ScaledContribution(log10, coeff, min, max);
    }

    // ── Len_3UTR ──────────────────────────────────────────────────────────────────────────
    // Port of get_len3UTR_weighted_contribution (single-isoform: UTR length = AIR_end = mrna.Length).
    private static double Len3UtrContribution(string mrna, TargetSiteType type)
    {
        int utrLength = mrna.Length;
        double log10 = utrLength > 0 ? Math.Log10(utrLength) : 0.0;

        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxLen3UtrCoeff8mer, CtxLen3UtrMin8mer, CtxLen3UtrMax8mer),
            TargetSiteType.Seed7merM8 => (CtxLen3UtrCoeff7merM8, CtxLen3UtrMin7merM8, CtxLen3UtrMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxLen3UtrCoeff7merA1, CtxLen3UtrMin7merA1, CtxLen3UtrMax7merA1),
            _ => (CtxLen3UtrCoeff6mer, CtxLen3UtrMin6mer, CtxLen3UtrMax6mer)
        };
        return ScaledContribution(log10, coeff, min, max);
    }

    // ── Off6m ─────────────────────────────────────────────────────────────────────────────
    // Port of getOffset6merSites + getOffset6mer_weighted_contribution (single-isoform). The
    // offset-6mer target pattern is the first 6 nt of reverse_complement(seed region 2-8); Off6m is
    // the count of (case-insensitive, overlapping) occurrences of that 6mer in the 3'UTR. Used RAW
    // (Off6m is not in the min-max regex of getAgarwalContribution).
    private static double Off6mContribution(string mrna, string mirna, TargetSiteType type)
    {
        if (mrna.Length == 0 || mirna.Length < 8) return 0.0;

        string seedRegion = mirna.Substring(1, 7);          // miRNA nt 2-8
        string revComp = ReverseComplementRna(seedRegion);  // 7 nt
        string pattern = revComp.Substring(0, 6);           // first 6 nt

        int count = 0;
        for (int i = 0; i + 6 <= mrna.Length; i++)
        {
            if (string.CompareOrdinal(mrna, i, pattern, 0, 6) == 0)
                count++;
        }

        double coeff = type switch
        {
            TargetSiteType.Seed8mer => CtxOff6mCoeff8mer,
            TargetSiteType.Seed7merM8 => CtxOff6mCoeff7merM8,
            TargetSiteType.Seed7merA1 => CtxOff6mCoeff7merA1,
            _ => CtxOff6mCoeff6mer
        };
        return coeff * count; // used raw (no scaling)
    }

    // reverse_complement over RNA (A↔U, G↔C); mirrors the perl reverse_complement after T→U.
    private static string ReverseComplementRna(string seq)
    {
        var c = new char[seq.Length];
        for (int i = 0; i < seq.Length; i++)
        {
            char b = seq[seq.Length - 1 - i];
            c[i] = b switch { 'A' => 'U', 'U' => 'A', 'T' => 'A', 'G' => 'C', 'C' => 'G', _ => b };
        }
        return new string(c);
    }

    // ── SPS / TA / Len_ORF / ORF8m (caller-supplied) ─────────────────────────────────────────
    private static double SpsContribution(double sps, TargetSiteType type)
    {
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxSpsCoeff8mer, CtxSpsMin8mer, CtxSpsMax8mer),
            TargetSiteType.Seed7merM8 => (CtxSpsCoeff7merM8, CtxSpsMin7merM8, CtxSpsMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxSpsCoeff7merA1, CtxSpsMin7merA1, CtxSpsMax7merA1),
            _ => (CtxSpsCoeff6mer, CtxSpsMin6mer, CtxSpsMax6mer)
        };
        return ScaledContribution(sps, coeff, min, max);
    }

    private static double TaContribution(double ta, TargetSiteType type)
    {
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxTaCoeff8mer, CtxTaMin8mer, CtxTaMax8mer),
            TargetSiteType.Seed7merM8 => (CtxTaCoeff7merM8, CtxTaMin7merM8, CtxTaMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxTaCoeff7merA1, CtxTaMin7merA1, CtxTaMax7merA1),
            _ => (CtxTaCoeff6mer, CtxTaMin6mer, CtxTaMax6mer)
        };
        return ScaledContribution(ta, coeff, min, max);
    }

    private static double LenOrfContribution(double orfLength, TargetSiteType type)
    {
        double log10 = orfLength > 0 ? Math.Log10(orfLength) : 0.0;
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxLenOrfCoeff8mer, CtxLenOrfMin8mer, CtxLenOrfMax8mer),
            TargetSiteType.Seed7merM8 => (CtxLenOrfCoeff7merM8, CtxLenOrfMin7merM8, CtxLenOrfMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxLenOrfCoeff7merA1, CtxLenOrfMin7merA1, CtxLenOrfMax7merA1),
            _ => (CtxLenOrfCoeff6mer, CtxLenOrfMin6mer, CtxLenOrfMax6mer)
        };
        return ScaledContribution(log10, coeff, min, max);
    }

    private static double Orf8mContribution(int count, TargetSiteType type)
    {
        double coeff = type switch
        {
            TargetSiteType.Seed8mer => CtxOrf8mCoeff8mer,
            TargetSiteType.Seed7merM8 => CtxOrf8mCoeff7merM8,
            TargetSiteType.Seed7merA1 => CtxOrf8mCoeff7merA1,
            _ => CtxOrf8mCoeff6mer
        };
        return coeff * count; // used raw (no scaling)
    }

    // ── PCT (probability of conserved targeting) — Friedman et al. (2009) ───────────────────

    /// <summary>
    /// Computes the Friedman et al. (2009) <b>branch-length score (Bls)</b> for a seed match: the
    /// total branch length of the minimal subtree of <paramref name="tree"/> that connects the
    /// <paramref name="speciesWithSite"/> (the species in which the site is perfectly aligned).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Per Friedman et al. (2009) Genome Res 19:92 (Methods): conservation "was then assessed by
    /// summing the total branch length in the phylogenetic tree connecting the subset of species
    /// having the sequence perfectly aligned." An edge belongs to that connecting subtree iff at
    /// least one species-with-site lies on each side of the edge (i.e. the edge is on a path
    /// between two such species). A single species (or none) yields Bls = 0.
    /// </para>
    /// <para>
    /// Leaf matching uses the leaf <see cref="PhylogeneticAnalyzer.PhyloNode.Name"/>; only species
    /// present in <paramref name="speciesWithSite"/> that also appear as a leaf of the tree count.
    /// </para>
    /// </remarks>
    /// <param name="tree">A phylogenetic tree whose leaf names are species identifiers and whose <c>BranchLength</c> values are edge lengths (e.g. parsed via <see cref="PhylogeneticAnalyzer.ParseNewick"/>).</param>
    /// <param name="speciesWithSite">Species (leaf names) in which the site is conserved.</param>
    /// <returns>The branch-length score (≥ 0).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="tree"/> or <paramref name="speciesWithSite"/> is null.</exception>
    public static double ComputeBranchLengthScore(
        PhylogeneticAnalyzer.PhyloNode tree,
        IReadOnlyCollection<string> speciesWithSite)
    {
        ArgumentNullException.ThrowIfNull(tree);
        ArgumentNullException.ThrowIfNull(speciesWithSite);

        var present = new HashSet<string>(speciesWithSite, StringComparer.Ordinal);
        if (present.Count == 0)
            return 0.0;

        double total = 0.0;
        AccumulateConnectingBranchLength(tree, present, out _, ref total);
        return total;
    }

    // Post-order walk. Returns, via <paramref name="subtreeHasSite"/>, whether the subtree rooted at
    // <paramref name="node"/> contains ≥1 species-with-site. A child edge is added to the connecting
    // subtree iff the child's subtree contains a species-with-site AND there is at least one
    // species-with-site outside that subtree (so the edge lies on a connecting path). Branch length
    // of the root edge itself is never counted (the root has no parent edge within the tree).
    private static void AccumulateConnectingBranchLength(
        PhylogeneticAnalyzer.PhyloNode node,
        HashSet<string> present,
        out bool subtreeHasSite,
        ref double total)
    {
        if (node.IsLeaf)
        {
            subtreeHasSite = present.Contains(node.Name);
            return;
        }

        bool anyChildHasSite = false;
        foreach (var child in node.Children)
        {
            AccumulateConnectingBranchLength(child, present, out bool childHasSite, ref total);
            if (childHasSite)
            {
                anyChildHasSite = true;
                // The edge from node→child is on a connecting path iff there is a species-with-site
                // both inside (childHasSite) and outside the child's subtree. "Outside" means the
                // total present count exceeds the count within the child's subtree. We detect this
                // structurally below after the loop; for the common case we add the edge whenever
                // the child subtree does not contain ALL present species.
                if (CountPresentInSubtree(child, present) < present.Count)
                    total += child.BranchLength;
            }
        }

        subtreeHasSite = anyChildHasSite;
    }

    private static int CountPresentInSubtree(PhylogeneticAnalyzer.PhyloNode node, HashSet<string> present)
    {
        if (node.IsLeaf)
            return present.Contains(node.Name) ? 1 : 0;
        int sum = 0;
        foreach (var child in node.Children)
            sum += CountPresentInSubtree(child, present);
        return sum;
    }

    /// <summary>
    /// Maps a branch-length score (Bls) to a PCT value using the published TargetScan logistic
    /// relationship (<c>targetscan_70_BL_PCT.pl</c>, <c>calculatePCTthisBL</c>):
    /// <c>PCT(Bls) = B0 + B1 / (1 + e^(−B2·Bls + B3))</c>, truncated at 0.
    /// </summary>
    /// <param name="branchLengthScore">The Friedman et al. (2009) branch-length score.</param>
    /// <param name="p">The caller-supplied published per-site-type sigmoid parameters.</param>
    /// <returns>PCT ∈ [0, ∞) (in practice [0, 1]); negative values are truncated to 0 as in the reference.</returns>
    public static double PctFromBranchLength(double branchLengthScore, PctSigmoidParameters p)
    {
        // targetscan_70_BL_PCT.pl: $pct = $b0 + ( $b1 / (1 + e ** ( (0 - $b2) * $BL + $b3)));
        double pct = p.B0 + p.B1 / (1.0 + Math.Exp(-p.B2 * branchLengthScore + p.B3));
        // targetscan_70_BL_PCT.pl: if ($pct < 0) { $pct = "0.0"; }
        return pct < 0.0 ? 0.0 : pct;
    }

    // PCT enters getAgarwalContribution as a min-max-scaled feature × the site-type PCT coefficient
    // (targetscan_70_context_scores.pl getPCT_contribution → getAgarwalContribution; PCT matches the
    // min-max regex branch alongside TA_3UTR/SPS/Local_AU/3P_score/SA/Len_ORF/Len_3UTR/Min_dist).
    private static double PctContribution(double pct, TargetSiteType type)
    {
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxPctCoeff8mer, CtxPctMin8mer, CtxPctMax8mer),
            TargetSiteType.Seed7merM8 => (CtxPctCoeff7merM8, CtxPctMin7merM8, CtxPctMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxPctCoeff7merA1, CtxPctMin7merA1, CtxPctMax7merA1),
            _ => (CtxPctCoeff6mer, CtxPctMin6mer, CtxPctMax6mer)
        };
        return ScaledContribution(pct, coeff, min, max);
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

    #region Drosha/Dicer cleavage-site prediction (Han 2006 / Park 2011 measuring rules) — opt-in

    /// <summary>
    /// Predicted Drosha (basal) and Dicer (apical) cleavage coordinates of a pri-/pre-miRNA hairpin and
    /// the resulting mature-miRNA (5p) and miRNA* (3p) spans, derived from the PUBLISHED measuring rules.
    /// All coordinates are 0-based indices into the input <see cref="Sequence"/>; spans are inclusive.
    /// </summary>
    /// <remarks>
    /// Measuring rules (each constant cited at its declaration):
    /// <list type="bullet">
    /// <item>Drosha cleaves ~11 bp from the basal ssRNA-dsRNA junction (Han et al. 2006).</item>
    /// <item>Dicer cleaves ~22 nt from the Drosha-generated 5' end — the 5' counting rule (Park et al. 2011).</item>
    /// <item>Each RNase III (Drosha, Dicer) cut leaves a 2-nt 3' overhang (Lee et al. 2003; Han et al. 2006).</item>
    /// </list>
    /// </remarks>
    /// <param name="BasalJunction">0-based index of the first paired base of the 5' arm (the basal ssRNA-dsRNA junction).</param>
    /// <param name="DroshaCut5Prime">0-based index of the 5' end of the 5p mature = Drosha cut on the 5' arm (basal junction + 11 bp).</param>
    /// <param name="DroshaCut3Prime">0-based index of the base on the 3' arm cut by Drosha (3' end of the 3p mature); 2 nt 3' of the 5p-paired position (the 2-nt 3' overhang).</param>
    /// <param name="MatureStart">0-based start of the 5p mature (= <see cref="DroshaCut5Prime"/>).</param>
    /// <param name="MatureEnd">0-based inclusive end of the 5p mature (= MatureStart + 22 − 1; Dicer 5' counting rule).</param>
    /// <param name="StarStart">0-based start of the 3p (miRNA*) span on the 3' arm.</param>
    /// <param name="StarEnd">0-based inclusive end of the 3p (miRNA*) span (= <see cref="DroshaCut3Prime"/>).</param>
    /// <param name="MatureSequence">The predicted 5p mature subsequence of <see cref="Sequence"/>.</param>
    /// <param name="StarSequence">The predicted 3p (miRNA*) subsequence of <see cref="Sequence"/>.</param>
    /// <param name="Sequence">The upper-cased (T→U) sequence the coordinates index into.</param>
    /// <param name="ThreePrimeOverhang">Length of the RNase III 2-nt 3' overhang produced at each cut.</param>
    /// <param name="HasCnncMotif">True iff a CNNC motif is present 16–18 nt 3' of the Drosha basal cut (Auyeung et al. 2013) — an optional confidence flag, not required.</param>
    public readonly record struct DroshaDicerCleavage(
        int BasalJunction,
        int DroshaCut5Prime,
        int DroshaCut3Prime,
        int MatureStart,
        int MatureEnd,
        int StarStart,
        int StarEnd,
        string MatureSequence,
        string StarSequence,
        string Sequence,
        int ThreePrimeOverhang,
        bool HasCnncMotif);

    // Drosha basal-junction "ruler": Drosha cleaves ~11 bp (≈ one helical turn) from the basal
    // stem-ssRNA junction — Han J et al. (2006), Cell 125:887-901, doi:10.1016/j.cell.2006.03.043:
    // "The cleavage site is determined mainly by the distance (approximately 11 bp) from the
    // stem-ssRNA junction."
    private const int DroshaCutBpFromBasalJunction = 11;

    // Dicer 5'-counting "ruler": Dicer cleaves ~22 nt from the Drosha-generated 5' end — Park JE et al.
    // (2011), Nature 475:201-205, doi:10.1038/nature10198: "the cleavage site determined mainly by the
    // distance (∼22 nucleotides) from the 5' end (5' counting rule)." This also fixes the mature length.
    private const int DicerCutNtFrom5PrimeEnd = 22;

    // RNase III staggered cut leaves a 2-nt 3' overhang (both Drosha and Dicer) — Lee Y et al. (2003);
    // Han et al. (2006). "Cleavage by RNase III domains results in 2-nt 3'-overhang end."
    private const int RNaseIII3PrimeOverhang = 2;

    // CNNC motif is positioned 16-18 nt downstream (3') of the Drosha basal cut — Auyeung VC et al.
    // (2013), Cell 152:844-858, doi:10.1016/j.cell.2013.01.031. Used only as an optional confidence flag.
    private const int CnncMinNtDownstreamOfDroshaCut = 16;
    private const int CnncMaxNtDownstreamOfDroshaCut = 18;
    private const int CnncMotifLength = 4; // C-N-N-C

    /// <summary>
    /// Predicts the Drosha (basal) and Dicer (apical) cleavage coordinates and the resulting 5p mature
    /// / 3p (miRNA*) spans of a pri-/pre-miRNA hairpin from the PUBLISHED distance ("ruler") rules, as
    /// an OPT-IN method (the existing detection methods and their defaults are unchanged):
    /// Drosha cuts ~11 bp from the basal ssRNA-dsRNA junction (Han et al. 2006); Dicer cuts ~22 nt from
    /// the Drosha-generated 5' end — the 5' counting rule (Park et al. 2011); each RNase III cut leaves
    /// a 2-nt 3' overhang (Lee et al. 2003; Han et al. 2006). It does NOT use a trained classifier.
    /// </summary>
    /// <param name="sequence">The pri-/pre-miRNA nucleotide sequence (RNA or DNA; T read as U; case-insensitive).</param>
    /// <param name="basalJunction">
    /// 0-based index of the basal ssRNA-dsRNA junction = the first paired base of the 5' arm of the stem
    /// (the position from which the ~11 bp Drosha ruler is measured). For a pri-miRNA this is the first
    /// base after the flanking basal ssRNA segment; for a miRBase-trimmed pre-miRNA whose annotated stem
    /// base coincides with the Drosha cut, pass the index of the stem base.
    /// </param>
    /// <returns>
    /// The predicted cleavage coordinates and mature/star spans, or <c>null</c> when the sequence is
    /// null/empty, the junction is out of range, or the hairpin is too short to place the predicted cuts.
    /// </returns>
    public static DroshaDicerCleavage? PredictDroshaDicerCleavage(string sequence, int basalJunction)
    {
        if (string.IsNullOrEmpty(sequence))
            return null;
        if (basalJunction < 0 || basalJunction >= sequence.Length)
            return null;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // Drosha cut on the 5' arm = basal junction + 11 bp (Han 2006). This is the 5' end of the 5p
        // mature (the base of the pre-miRNA on the 5' arm).
        int droshaCut5Prime = basalJunction + DroshaCutBpFromBasalJunction;
        if (droshaCut5Prime >= upper.Length)
            return null;

        // 5p mature: Dicer 5' counting rule fixes the mature length at ~22 nt (Park 2011), so the
        // mature spans [droshaCut5Prime, droshaCut5Prime + 22 - 1].
        int matureStart = droshaCut5Prime;
        int matureEnd = matureStart + DicerCutNtFrom5PrimeEnd - 1;
        if (matureEnd >= upper.Length)
            return null;

        // The 3' arm is cut staggered by 2 nt (RNase III 2-nt 3' overhang). The 3p mature is the
        // ~22-nt span on the 3' arm whose 3' end (the Drosha-generated end on the 3' arm) sits 2 nt
        // 3' of the position paired with the 5p mature's 5' base. Modelled by the conventional duplex
        // geometry: the 3p mature 5' end pairs ~2 nt inside the 5p mature 3' end (the 2-nt 3' overhang
        // at the apical Dicer cut), and the 3p span has the same ~22-nt length.
        int starEnd = matureEnd + RNaseIII3PrimeOverhang; // 3p 3' end (Drosha-generated end on 3' arm)
        int starStart = starEnd - DicerCutNtFrom5PrimeEnd + 1;
        if (starEnd >= upper.Length || starStart < 0)
            return null;

        int droshaCut3Prime = starEnd;

        string matureSeq = upper.Substring(matureStart, matureEnd - matureStart + 1);
        string starSeq = upper.Substring(starStart, starEnd - starStart + 1);

        bool hasCnnc = HasCnncMotifDownstream(upper, droshaCut5Prime);

        return new DroshaDicerCleavage(
            BasalJunction: basalJunction,
            DroshaCut5Prime: droshaCut5Prime,
            DroshaCut3Prime: droshaCut3Prime,
            MatureStart: matureStart,
            MatureEnd: matureEnd,
            StarStart: starStart,
            StarEnd: starEnd,
            MatureSequence: matureSeq,
            StarSequence: starSeq,
            Sequence: upper,
            ThreePrimeOverhang: RNaseIII3PrimeOverhang,
            HasCnncMotif: hasCnnc);
    }

    /// <summary>
    /// Returns true iff a CNNC motif (C, any, any, C) occurs within the 16–18 nt window 3' of the
    /// Drosha basal cut (Auyeung et al. 2013). Optional confidence signal only.
    /// </summary>
    private static bool HasCnncMotifDownstream(string seq, int droshaCut5Prime)
    {
        // Window: positions [cut + 16, cut + 18] start a CNNC (4-mer) downstream of the basal cut.
        for (int offset = CnncMinNtDownstreamOfDroshaCut; offset <= CnncMaxNtDownstreamOfDroshaCut; offset++)
        {
            int start = droshaCut5Prime + offset;
            if (start < 0 || start + CnncMotifLength > seq.Length)
                continue;
            if (seq[start] == 'C' && seq[start + CnncMotifLength - 1] == 'C')
                return true;
        }
        return false;
    }

    #endregion

    #region Pre-miRNA natural-vs-background classifier (trained logistic regression) — opt-in

    /// <summary>
    /// Structure/sequence features of a candidate pre-miRNA hairpin, computed entirely from the
    /// REAL minimum-free-energy (MFE) secondary structure of the validated Zuker–Stiegler folder
    /// (<see cref="RnaSecondaryStructure.CalculateMfeStructure"/>, RNA-STRUCT-001) plus the base
    /// composition. These are the inputs to the trained natural-vs-background logistic-regression
    /// classifier (<see cref="ClassifyPreMiRna"/>). Each feature is published:
    /// <list type="bullet">
    /// <item><b>MFE</b> — minimum free energy ΔG° (kcal/mol), the discriminative thermodynamic signal
    /// of Bonnet et al. (2004), Bioinformatics 20(17):2911-2917.</item>
    /// <item><b>AMFE / MFEI</b> — adjusted MFE = 100·|ΔG°|/length and MFEI = AMFE/(G+C)% of Zhang et
    /// al. (2006), Cell Mol Life Sci 63:246-254.</item>
    /// <item><b>PairedFraction</b> — fraction of bases that are paired in the MFE structure (the
    /// base-pairing-propensity / %paired structural feature family of microPred, Batuwita &amp;
    /// Palade (2009), Bioinformatics 25(8):989-995).</item>
    /// <item><b>GcContent, StemBasePairs, LoopSize</b> — composition and stem/loop geometry, the
    /// sequence/structure features of microPred (Batuwita &amp; Palade 2009) and Xue et al. (2005),
    /// BMC Bioinformatics 6:310.</item>
    /// </list>
    /// </summary>
    /// <param name="FreeEnergy">ΔG° of the MFE structure (kcal/mol; negative for a stable hairpin).</param>
    /// <param name="Amfe">Adjusted MFE = 100·|ΔG°|/length (Zhang 2006).</param>
    /// <param name="Mfei">MFEI = AMFE/(G+C)% (Zhang 2006).</param>
    /// <param name="GcContent">G+C content as a fraction in [0,1].</param>
    /// <param name="PairedFraction">Fraction of bases paired in the MFE structure, in [0,1].</param>
    /// <param name="StemBasePairs">Number of base pairs in the MFE structure.</param>
    /// <param name="LoopSize">Largest single run of unpaired bases enclosed by the stem (apical loop).</param>
    /// <param name="Length">Sequence length in nucleotides.</param>
    public readonly record struct PreMiRnaFeatures(
        double FreeEnergy,
        double Amfe,
        double Mfei,
        double GcContent,
        double PairedFraction,
        int StemBasePairs,
        int LoopSize,
        int Length);

    /// <summary>
    /// Result of the trained natural-vs-background pre-miRNA classifier
    /// (<see cref="ClassifyPreMiRna"/>): the extracted features, the logistic probability that the
    /// candidate is a genuine (natural) pre-miRNA hairpin rather than a background/shuffled sequence,
    /// and the resulting boolean call at the chosen decision threshold.
    /// </summary>
    /// <param name="Features">The structure/sequence features fed to the model.</param>
    /// <param name="NaturalProbability">P(natural) ∈ [0,1] from the logistic-regression model.</param>
    /// <param name="IsNatural">True iff <see cref="NaturalProbability"/> ≥ the decision threshold.</param>
    public readonly record struct PreMiRnaClassification(
        PreMiRnaFeatures Features,
        double NaturalProbability,
        bool IsNatural);

    // ── Trained logistic-regression coefficients ──────────────────────────────────────────
    // Model: logistic regression P(natural) = sigmoid(b0 + Σ b_j · z_j), where z_j is the j-th
    // feature standardised by the training-set mean/std (textbook logistic regression fit by
    // batch gradient ascent on the log-likelihood — Hastie, Tibshirani & Friedman, "The Elements
    // of Statistical Learning" 2nd ed. §4.4.1, doi:10.1007/978-0-387-84858-7).
    //
    // TRAINING DATA (public domain): positives = 13 real human pre-miRNA hairpin precursors from
    // miRBase (https://mirbase.org, Public Domain Dedication — https://mirbase.org/download/CURRENT/LICENSE/);
    // negatives = dinucleotide-preserving shuffles of the positives (Altschul & Erickson 1985
    // Eulerian-walk shuffle; the standard pre-miRNA-classifier background convention of Bonnet
    // et al. 2004 / Xue et al. 2005 / Batuwita & Palade 2009). 4 shuffles per positive → 52
    // negatives; 65 examples total. Fixed RNG seed 20060101 and a fixed deterministic train/test
    // split (first 70% by example index = train, last 30% = held-out test) make the fit and the
    // held-out metric fully reproducible. NO GPL miRDeep2 code was used; only the published method.
    // Coefficients regenerated by the offline trainer ScratchPreMiRnaTrainer (see Evidence doc §Training).
    //
    // Feature order: [FreeEnergy, Amfe, Mfei, GcContent, PairedFraction]. (StemBasePairs / LoopSize /
    // Length are extracted and reported but not model inputs — they are collinear with MFE/PairedFraction
    // and length-confounded; the five thermodynamic/composition/pairing features are the model.)
    private const double PreMiRnaLogRegBias = -4.340788257901692;
    private static readonly double[] PreMiRnaLogRegWeights =
        { -0.6299061497891402, 2.4689410267337104, 2.5770103435217253, -0.3315361036977469, 2.42751798356881 };
    private static readonly double[] PreMiRnaFeatureMean =
        { -22.988666666666667, 29.297338207740115, 0.6157639504286563, 0.477134229110014, 0.6153156920227864 };
    private static readonly double[] PreMiRnaFeatureStd =
        { 8.373410721510746, 10.37385177977526, 0.2084972456212268, 0.061249038538481466, 0.10661272255434924 };

    // Decision threshold: 0.5 is the natural Bayes cutoff for a balanced two-class logistic model
    // (Hastie, Tibshirani & Friedman §4.4). Exposed as the default of ClassifyPreMiRna.
    private const double DefaultNaturalThreshold = 0.5;

    /// <summary>
    /// Extracts the published structure/sequence features of a candidate pre-miRNA hairpin from its
    /// REAL minimum-free-energy secondary structure (folded once with the validated Zuker–Stiegler
    /// engine, RNA-STRUCT-001). The features are the inputs to <see cref="ClassifyPreMiRna"/>.
    /// </summary>
    /// <param name="sequence">Candidate nucleotide sequence (RNA or DNA; T read as U; case-insensitive).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    /// <returns>The extracted features, or <c>null</c> when the sequence is null/empty.</returns>
    public static PreMiRnaFeatures? ExtractPreMiRnaFeatures(string sequence, int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            return null;

        var mfe = RnaSecondaryStructure.CalculateMfeStructure(sequence, minLoopSize);
        string seq = mfe.Sequence;
        int n = seq.Length;
        if (n == 0)
            return null;

        int bp = mfe.BasePairs.Count;
        double pairedFraction = n > 0 ? (2.0 * bp) / n : 0.0; // each pair occupies two bases
        double gcFraction = GcPercent(seq) / AmfePer100Nt;    // GcPercent returns a percentage
        double amfe = AmfePer100Nt * Math.Abs(mfe.FreeEnergy) / n;
        double mfei = CalculateMfeIndex(mfe.FreeEnergy, n, GcPercent(seq));
        int loopSize = LargestEnclosedLoop(mfe.DotBracket);

        return new PreMiRnaFeatures(
            FreeEnergy: mfe.FreeEnergy,
            Amfe: amfe,
            Mfei: mfei,
            GcContent: gcFraction,
            PairedFraction: pairedFraction,
            StemBasePairs: bp,
            LoopSize: loopSize,
            Length: n);
    }

    /// <summary>
    /// Classifies a candidate pre-miRNA hairpin as a genuine (natural) pre-miRNA or a
    /// background/shuffled sequence using the bundled trained logistic-regression model over the
    /// published structure/sequence features (see <see cref="PreMiRnaFeatures"/>). This is an OPT-IN
    /// method — the existing detection methods (<see cref="FindPreMiRnaHairpins"/>,
    /// <see cref="FindPreMiRnaHairpinsByMfe"/>) and their defaults are unchanged. The model was
    /// trained from PUBLIC-DOMAIN miRBase precursors against dinucleotide-shuffled negatives
    /// (Bonnet et al. 2004 convention); NO trained GPL miRDeep2 code or weights are used.
    /// </summary>
    /// <param name="sequence">Candidate nucleotide sequence (RNA or DNA; T read as U; case-insensitive).</param>
    /// <param name="threshold">Decision threshold on P(natural) (default 0.5, the Bayes cutoff).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    /// <returns>The features, P(natural), and the boolean call; <c>null</c> for null/empty input.</returns>
    public static PreMiRnaClassification? ClassifyPreMiRna(
        string sequence,
        double threshold = DefaultNaturalThreshold,
        int minLoopSize = 3)
    {
        var features = ExtractPreMiRnaFeatures(sequence, minLoopSize);
        if (features is null)
            return null;

        double p = ScorePreMiRnaFeatures(features.Value);
        return new PreMiRnaClassification(features.Value, p, p >= threshold);
    }

    /// <summary>
    /// Evaluates the bundled trained logistic-regression model on a feature vector, returning
    /// P(natural) ∈ [0,1]. Exposed so callers can score features they extracted themselves.
    /// </summary>
    public static double ScorePreMiRnaFeatures(PreMiRnaFeatures features)
    {
        double[] raw =
        {
            features.FreeEnergy,
            features.Amfe,
            features.Mfei,
            features.GcContent,
            features.PairedFraction
        };

        double z = PreMiRnaLogRegBias;
        for (int j = 0; j < raw.Length; j++)
        {
            double std = PreMiRnaFeatureStd[j] == 0 ? 1.0 : PreMiRnaFeatureStd[j];
            double standardised = (raw[j] - PreMiRnaFeatureMean[j]) / std;
            z += PreMiRnaLogRegWeights[j] * standardised;
        }

        return Sigmoid(z);
    }

    /// <summary>Logistic (sigmoid) link 1/(1+e^-z) — Hastie, Tibshirani &amp; Friedman §4.4.1.</summary>
    private static double Sigmoid(double z) => 1.0 / (1.0 + Math.Exp(-z));

    /// <summary>
    /// Largest single run of unpaired bases ('.') that is enclosed within the structure (i.e. lies
    /// after the first '(' and before the last ')') — the apical loop of a hairpin. Returns 0 when
    /// there is no enclosed unpaired run.
    /// </summary>
    private static int LargestEnclosedLoop(string dotBracket)
    {
        int firstOpen = dotBracket.IndexOf('(');
        int lastClose = dotBracket.LastIndexOf(')');
        if (firstOpen < 0 || lastClose < 0 || lastClose <= firstOpen)
            return 0;

        int best = 0, run = 0;
        for (int k = firstOpen + 1; k < lastClose; k++)
        {
            if (dotBracket[k] == '.') { run++; if (run > best) best = run; }
            else run = 0;
        }
        return best;
    }

    /// <summary>
    /// Generates a dinucleotide-composition-preserving shuffle of <paramref name="sequence"/> using
    /// the Altschul &amp; Erickson (1985) Eulerian-walk algorithm: a random Eulerian trail through the
    /// directed multigraph whose edges are the input's adjacent-nucleotide pairs reproduces a sequence
    /// with the SAME first nucleotide, last nucleotide, and exact dinucleotide (doublet) counts. This
    /// is the standard pre-miRNA-classifier background convention (Bonnet et al. 2004; Xue et al. 2005;
    /// Batuwita &amp; Palade 2009): a length/composition-matched negative that preserves local
    /// composition so that discrimination must come from STRUCTURE, not base content.
    /// </summary>
    /// <param name="sequence">Input nucleotide sequence (upper/lowercase; characters preserved as-is except case is normalised to upper, T→U).</param>
    /// <param name="random">RNG (supply a fixed seed for reproducibility).</param>
    /// <returns>A shuffled sequence with identical dinucleotide counts; the input unchanged for length &lt; 2.</returns>
    public static string DinucleotideShuffle(string sequence, Random random)
    {
        if (random is null) throw new ArgumentNullException(nameof(random));
        if (string.IsNullOrEmpty(sequence)) return sequence ?? "";

        string s = sequence.ToUpperInvariant().Replace('T', 'U');
        int n = s.Length;
        if (n < 2) return s;

        char first = s[0];
        char last = s[n - 1];

        // Build the edge multiset: for each adjacent pair s[i]s[i+1], an edge s[i] -> s[i+1].
        var outEdges = new Dictionary<char, List<char>>();
        for (int i = 0; i < n - 1; i++)
        {
            if (!outEdges.TryGetValue(s[i], out var list))
            {
                list = new List<char>();
                outEdges[s[i]] = list;
            }
            list.Add(s[i + 1]);
        }

        var vertices = new List<char>(outEdges.Keys);

        // Altschul-Erickson: repeatedly (a) pick, for every vertex except `last`, one random outgoing
        // edge to form a "last-edge" arborescence rooted at `last`; the chosen last-edges must NOT
        // form a cycle (every vertex must reach `last`). (b) Shuffle the remaining edges. (c) Walk an
        // Eulerian trail from `first`, always taking the chosen last-edge of a vertex LAST. Retry the
        // arborescence draw until it is acyclic (a connected Eulerian multigraph guarantees success).
        for (int attempt = 0; attempt < 1000; attempt++)
        {
            var lastEdge = new Dictionary<char, char>();
            bool ok = true;
            foreach (char v in vertices)
            {
                if (v == last) continue;
                var list = outEdges[v];
                if (list.Count == 0) { ok = false; break; }
                lastEdge[v] = list[random.Next(list.Count)];
            }
            if (!ok) { if (vertices.Count == 1) break; else continue; }

            // Verify the chosen last-edges form a tree into `last` (no cycle that avoids `last`).
            if (!LastEdgesReachRoot(lastEdge, last, vertices))
                continue;

            // Build a per-vertex ordering of outgoing edges: random order, but the chosen last-edge
            // placed at the END so the Eulerian walk uses it last.
            var ordered = new Dictionary<char, List<char>>();
            foreach (char v in vertices)
            {
                var src = new List<char>(outEdges[v]);
                FisherYatesShuffle(src, random);
                if (lastEdge.TryGetValue(v, out char le))
                {
                    int idx = src.IndexOf(le);
                    if (idx >= 0) { src.RemoveAt(idx); src.Add(le); }
                }
                ordered[v] = src;
            }

            // Walk the Eulerian trail from `first`, consuming each vertex's edges in `ordered` order.
            var pos = new Dictionary<char, int>();
            foreach (char v in vertices) pos[v] = 0;
            var sb = new StringBuilder(n);
            sb.Append(first);
            char cur = first;
            for (int step = 0; step < n - 1; step++)
            {
                if (!ordered.TryGetValue(cur, out var list) || pos[cur] >= list.Count)
                {
                    ok = false;
                    break;
                }
                char nxt = list[pos[cur]++];
                sb.Append(nxt);
                cur = nxt;
            }
            if (ok && sb.Length == n)
                return sb.ToString();
        }

        // Fallback (extremely rare / degenerate composition): return the original (still composition-matched).
        return s;
    }

    // True iff, following the single chosen last-edge from every non-root vertex, every vertex
    // reaches `root` (i.e. the last-edges form an arborescence into root with no avoidant cycle).
    private static bool LastEdgesReachRoot(Dictionary<char, char> lastEdge, char root, List<char> vertices)
    {
        foreach (char v in vertices)
        {
            if (v == root) continue;
            char cur = v;
            int guard = 0;
            var seen = new HashSet<char>();
            while (cur != root)
            {
                if (!seen.Add(cur)) return false;           // cycle not through root
                if (!lastEdge.TryGetValue(cur, out cur)) return false; // dead end
                if (++guard > vertices.Count + 1) return false;
            }
        }
        return true;
    }

    private static void FisherYatesShuffle<T>(IList<T> list, Random random)
    {
        for (int i = list.Count - 1; i > 0; i--)
        {
            int j = random.Next(i + 1);
            (list[i], list[j]) = (list[j], list[i]);
        }
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
