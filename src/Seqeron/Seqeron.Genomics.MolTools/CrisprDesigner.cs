using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Designs CRISPR guide RNAs (gRNAs) and identifies PAM sequences.
/// Supports Cas9 (SpCas9, SaCas9) and Cas12a (Cpf1) systems.
/// </summary>
public static class CrisprDesigner
{
    #region PAM Definitions

    /// <summary>
    /// Gets the PAM sequence for a specific CRISPR system.
    /// </summary>
    public static CrisprSystem GetSystem(CrisprSystemType type) => type switch
    {
        CrisprSystemType.SpCas9 => new CrisprSystem("SpCas9", "NGG", 20, true, "Streptococcus pyogenes Cas9"),
        CrisprSystemType.SpCas9_NAG => new CrisprSystem("SpCas9-NAG", "NAG", 20, true, "SpCas9 with NAG PAM (lower activity)"),
        CrisprSystemType.SaCas9 => new CrisprSystem("SaCas9", "NNGRRT", 21, true, "Staphylococcus aureus Cas9"),
        CrisprSystemType.Cas12a => new CrisprSystem("Cas12a/Cpf1", "TTTV", 23, false, "Cas12a (Cpf1) - PAM before target"),
        CrisprSystemType.AsCas12a => new CrisprSystem("AsCas12a", "TTTV", 23, false, "Acidaminococcus sp. Cas12a"),
        CrisprSystemType.LbCas12a => new CrisprSystem("LbCas12a", "TTTV", 24, false, "Lachnospiraceae bacterium Cas12a"),
        CrisprSystemType.CasX => new CrisprSystem("CasX", "TTCN", 20, false, "CasX - compact Cas protein"),
        _ => throw new ArgumentException($"Unknown CRISPR system: {type}")
    };

    #endregion

    #region PAM Finding

    /// <summary>
    /// Finds all PAM sites in a sequence for a given CRISPR system.
    /// </summary>
    /// <param name="sequence">DNA sequence to search.</param>
    /// <param name="systemType">Type of CRISPR system.</param>
    /// <returns>Collection of PAM sites with their positions and orientations.</returns>
    public static IEnumerable<PamSite> FindPamSites(
        DnaSequence sequence,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindPamSitesCore(sequence.Sequence, GetSystem(systemType));
    }

    /// <summary>
    /// Finds all PAM sites in a raw sequence string.
    /// </summary>
    public static IEnumerable<PamSite> FindPamSites(
        string sequence,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        foreach (var site in FindPamSitesCore(sequence.ToUpperInvariant(), GetSystem(systemType)))
            yield return site;
    }

    private static IEnumerable<PamSite> FindPamSitesCore(string seq, CrisprSystem system)
    {
        string pamPattern = system.PamSequence;
        int guideLength = system.GuideLength;
        bool pamAfterTarget = system.PamAfterTarget;

        // Search forward strand
        for (int i = 0; i <= seq.Length - pamPattern.Length; i++)
        {
            if (MatchesPam(seq, i, pamPattern))
            {
                int targetStart, targetEnd;
                if (pamAfterTarget)
                {
                    // PAM is after target (Cas9): target-PAM
                    targetStart = i - guideLength;
                    targetEnd = i - 1;
                }
                else
                {
                    // PAM is before target (Cas12a): PAM-target
                    targetStart = i + pamPattern.Length;
                    targetEnd = targetStart + guideLength - 1;
                }

                if (targetStart >= 0 && targetEnd < seq.Length)
                {
                    string target = seq.Substring(targetStart, guideLength);
                    yield return new PamSite(
                        Position: i,
                        PamSequence: seq.Substring(i, pamPattern.Length),
                        TargetSequence: target,
                        TargetStart: targetStart,
                        IsForwardStrand: true,
                        System: system);
                }
            }
        }

        // Search reverse strand
        string revComp = DnaSequence.GetReverseComplementString(seq);
        for (int i = 0; i <= revComp.Length - pamPattern.Length; i++)
        {
            if (MatchesPam(revComp, i, pamPattern))
            {
                int revTargetStart, revTargetEnd;
                if (pamAfterTarget)
                {
                    revTargetStart = i - guideLength;
                    revTargetEnd = i - 1;
                }
                else
                {
                    revTargetStart = i + pamPattern.Length;
                    revTargetEnd = revTargetStart + guideLength - 1;
                }

                if (revTargetStart >= 0 && revTargetEnd < revComp.Length)
                {
                    // Convert position back to forward strand coordinates
                    int forwardPos = seq.Length - i - pamPattern.Length;
                    string target = revComp.Substring(revTargetStart, guideLength);

                    yield return new PamSite(
                        Position: forwardPos,
                        PamSequence: DnaSequence.GetReverseComplementString(revComp.Substring(i, pamPattern.Length)),
                        TargetSequence: target,
                        TargetStart: revTargetStart,
                        IsForwardStrand: false,
                        System: system);
                }
            }
        }
    }

    private static bool MatchesPam(string seq, int position, string pamPattern)
    {
        if (position + pamPattern.Length > seq.Length)
            return false;

        for (int i = 0; i < pamPattern.Length; i++)
        {
            char seqChar = seq[position + i];
            char pamChar = pamPattern[i];

            if (!MatchesIupac(seqChar, pamChar))
                return false;
        }

        return true;
    }

    private static bool MatchesIupac(char nucleotide, char iupacCode) => IupacHelper.MatchesIupac(nucleotide, iupacCode);

    #endregion

    #region Guide RNA Design

    /// <summary>
    /// Designs guide RNAs for a target region.
    /// </summary>
    /// <param name="sequence">DNA sequence containing the target region.</param>
    /// <param name="regionStart">Start of the region to target (0-based).</param>
    /// <param name="regionEnd">End of the region to target (0-based, inclusive).</param>
    /// <param name="systemType">Type of CRISPR system.</param>
    /// <param name="parameters">Optional design parameters.</param>
    /// <returns>Ranked list of guide RNA candidates.</returns>
    public static IEnumerable<GuideRnaCandidate> DesignGuideRnas(
        DnaSequence sequence,
        int regionStart,
        int regionEnd,
        CrisprSystemType systemType = CrisprSystemType.SpCas9,
        GuideRnaParameters? parameters = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (regionStart < 0 || regionStart >= sequence.Length)
            throw new ArgumentOutOfRangeException(nameof(regionStart));
        if (regionEnd < regionStart || regionEnd >= sequence.Length)
            throw new ArgumentOutOfRangeException(nameof(regionEnd));

        var effectiveParams = parameters ?? GuideRnaParameters.Default;
        var system = GetSystem(systemType);

        var pamSites = FindPamSitesCore(sequence.Sequence, system)
            .Where(p => IsInRegion(p, regionStart, regionEnd, system))
            .ToList();

        foreach (var pamSite in pamSites)
        {
            var candidate = EvaluateGuideRna(pamSite, sequence.Sequence, parameters, system);
            if (candidate.Score >= effectiveParams.MinScore)
                yield return candidate;
        }
    }

    /// <summary>
    /// Evaluates a single guide RNA sequence.
    /// </summary>
    public static GuideRnaCandidate EvaluateGuideRna(
        string guideSequence,
        CrisprSystemType systemType = CrisprSystemType.SpCas9,
        GuideRnaParameters? parameters = null)
    {
        if (string.IsNullOrEmpty(guideSequence))
            throw new ArgumentNullException(nameof(guideSequence));

        var effectiveParams = parameters ?? GuideRnaParameters.Default;
        var system = GetSystem(systemType);
        var seq = guideSequence.ToUpperInvariant();

        // Calculate GC content
        double gcContent = CalculateGcContent(seq);

        // Check for polyT (transcription terminator)
        bool hasPolyT = HasPolyT(seq, 4);

        // Calculate self-complementarity score
        double selfCompScore = CalculateSelfComplementarity(seq);

        // Check seed region (last 12 bp for Cas9)
        string seedRegion = system.PamAfterTarget
            ? seq.Substring(Math.Max(0, seq.Length - 12))
            : seq.Substring(0, Math.Min(12, seq.Length));
        double seedGc = CalculateGcContent(seedRegion);

        // Calculate overall score
        var issues = new List<string>();
        double score = 100;

        // GC content penalty
        if (gcContent < effectiveParams.MinGcContent)
        {
            score -= (effectiveParams.MinGcContent - gcContent) * 2;
            issues.Add($"Low GC content ({gcContent:F1}%)");
        }
        else if (gcContent > effectiveParams.MaxGcContent)
        {
            score -= (gcContent - effectiveParams.MaxGcContent) * 2;
            issues.Add($"High GC content ({gcContent:F1}%)");
        }

        // PolyT penalty
        if (hasPolyT)
        {
            score -= 20;
            issues.Add("Contains TTTT (potential Pol III terminator)");
        }

        // Self-complementarity penalty
        if (selfCompScore > 0.3)
        {
            score -= selfCompScore * 30;
            issues.Add("High self-complementarity");
        }

        // Seed region GC penalty
        if (seedGc < 30 || seedGc > 70)
        {
            score -= 10;
            issues.Add($"Suboptimal seed region GC ({seedGc:F1}%)");
        }

        // Check for restriction sites
        bool hasRestrictionSite = HasCommonRestrictionSite(seq);
        if (hasRestrictionSite)
        {
            score -= 5;
            issues.Add("Contains common restriction site");
        }

        return new GuideRnaCandidate(
            Sequence: seq,
            Position: -1,
            IsForwardStrand: true,
            GcContent: gcContent,
            SeedGcContent: seedGc,
            HasPolyT: hasPolyT,
            SelfComplementarityScore: selfCompScore,
            Score: Math.Max(0, score),
            Issues: issues,
            System: system);
    }

    private static GuideRnaCandidate EvaluateGuideRna(
        PamSite pamSite,
        string fullSequence,
        GuideRnaParameters? parameters,
        CrisprSystem system)
    {
        var systemType = system.Name switch
        {
            "SpCas9" => CrisprSystemType.SpCas9,
            "SaCas9" => CrisprSystemType.SaCas9,
            "Cas12a/Cpf1" => CrisprSystemType.Cas12a,
            _ => CrisprSystemType.SpCas9
        };

        var candidate = EvaluateGuideRna(pamSite.TargetSequence, systemType, parameters);

        return candidate with
        {
            Position = pamSite.TargetStart,
            IsForwardStrand = pamSite.IsForwardStrand
        };
    }

    private static bool IsInRegion(PamSite pamSite, int regionStart, int regionEnd, CrisprSystem system)
    {
        // Check if the cut site would be within the target region
        int cutSite = system.PamAfterTarget
            ? pamSite.Position - 3 // Cas9 cuts 3bp upstream of PAM
            : pamSite.Position + system.PamSequence.Length + 18; // Cas12a cuts downstream

        return cutSite >= regionStart && cutSite <= regionEnd;
    }

    #endregion

    #region Off-Target Analysis

    /// <summary>
    /// Predicts potential off-target sites for a guide RNA.
    /// </summary>
    /// <param name="guideSequence">The guide RNA sequence.</param>
    /// <param name="genome">The genome/sequence to search for off-targets.</param>
    /// <param name="maxMismatches">Maximum number of mismatches allowed.</param>
    /// <param name="systemType">Type of CRISPR system.</param>
    /// <returns>Collection of potential off-target sites.</returns>
    public static IEnumerable<OffTargetSite> FindOffTargets(
        string guideSequence,
        DnaSequence genome,
        int maxMismatches = 3,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
    {
        if (string.IsNullOrEmpty(guideSequence))
            throw new ArgumentNullException(nameof(guideSequence));
        ArgumentNullException.ThrowIfNull(genome);
        if (maxMismatches < 0 || maxMismatches > 5)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches));

        var system = GetSystem(systemType);
        var guide = guideSequence.ToUpperInvariant();

        // Find all PAM sites
        var pamSites = FindPamSitesCore(genome.Sequence, system);

        foreach (var pamSite in pamSites)
        {
            int mismatches = CountMismatches(guide, pamSite.TargetSequence);

            if (mismatches > 0 && mismatches <= maxMismatches)
            {
                // Calculate off-target score based on mismatch positions
                double score = CalculateOffTargetScore(guide, pamSite.TargetSequence, system);

                yield return new OffTargetSite(
                    Position: pamSite.Position,
                    Sequence: pamSite.TargetSequence,
                    Mismatches: mismatches,
                    MismatchPositions: GetMismatchPositions(guide, pamSite.TargetSequence),
                    IsForwardStrand: pamSite.IsForwardStrand,
                    OffTargetScore: score);
            }
        }
    }

    /// <summary>
    /// Calculates specificity score for a guide RNA (higher = more specific).
    /// </summary>
    public static double CalculateSpecificityScore(
        string guideSequence,
        DnaSequence genome,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
    {
        var offTargets = FindOffTargets(guideSequence, genome, 4, systemType).ToList();

        if (offTargets.Count == 0)
            return 100.0;

        // Score decreases with number and quality of off-targets
        double totalPenalty = offTargets.Sum(ot => ot.OffTargetScore);
        return Math.Max(0, 100 - totalPenalty);
    }

    private static int CountMismatches(string seq1, string seq2)
    {
        if (seq1.Length != seq2.Length)
            return int.MaxValue;

        int mismatches = 0;
        for (int i = 0; i < seq1.Length; i++)
        {
            if (seq1[i] != seq2[i])
                mismatches++;
        }
        return mismatches;
    }

    private static IReadOnlyList<int> GetMismatchPositions(string guide, string target)
    {
        var positions = new List<int>();
        int len = Math.Min(guide.Length, target.Length);

        for (int i = 0; i < len; i++)
        {
            if (guide[i] != target[i])
                positions.Add(i);
        }

        return positions;
    }

    private static double CalculateOffTargetScore(string guide, string target, CrisprSystem system)
    {
        double score = 0;
        int seedStart = system.PamAfterTarget ? guide.Length - 12 : 0;
        int seedEnd = system.PamAfterTarget ? guide.Length : 12;

        for (int i = 0; i < Math.Min(guide.Length, target.Length); i++)
        {
            if (guide[i] != target[i])
            {
                // Mismatches in seed region are more tolerated (lower off-target activity)
                // but still concerning
                bool inSeed = i >= seedStart && i < seedEnd;
                score += inSeed ? 5 : 2;
            }
        }

        return score;
    }

    #endregion

    #region Helper Methods

    private static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcContentFast();

    private static bool HasPolyT(string sequence, int length)
    {
        return sequence.Contains(new string('T', length));
    }

    private static double CalculateSelfComplementarity(string sequence)
    {
        string revComp = DnaSequence.GetReverseComplementString(sequence);
        int matches = 0;
        int total = sequence.Length;

        // Check for complementary regions
        for (int offset = 0; offset < sequence.Length; offset++)
        {
            for (int i = 0; i + offset < sequence.Length; i++)
            {
                if (sequence[i] == revComp[i + offset])
                    matches++;
            }
        }

        return (double)matches / (total * total);
    }

    private static bool HasCommonRestrictionSite(string sequence)
    {
        string[] commonSites = { "GAATTC", "GGATCC", "AAGCTT", "CTGCAG", "GCGGCCGC" };
        return commonSites.Any(site => sequence.Contains(site));
    }

    #endregion
}

/// <summary>
/// Type of CRISPR system.
/// </summary>
public enum CrisprSystemType
{
    /// <summary>Streptococcus pyogenes Cas9 (NGG PAM)</summary>
    SpCas9,
    /// <summary>SpCas9 with NAG PAM (lower efficiency)</summary>
    SpCas9_NAG,
    /// <summary>Staphylococcus aureus Cas9 (NNGRRT PAM)</summary>
    SaCas9,
    /// <summary>Cas12a/Cpf1 (TTTV PAM)</summary>
    Cas12a,
    /// <summary>Acidaminococcus sp. Cas12a</summary>
    AsCas12a,
    /// <summary>Lachnospiraceae bacterium Cas12a</summary>
    LbCas12a,
    /// <summary>CasX (TTCN PAM)</summary>
    CasX
}

/// <summary>
/// Represents a CRISPR system with its characteristics.
/// </summary>
public sealed record CrisprSystem(
    string Name,
    string PamSequence,
    int GuideLength,
    bool PamAfterTarget,
    string Description);

/// <summary>
/// Represents a PAM site in a sequence.
/// </summary>
public sealed record PamSite(
    int Position,
    string PamSequence,
    string TargetSequence,
    int TargetStart,
    bool IsForwardStrand,
    CrisprSystem System);

/// <summary>
/// Parameters for guide RNA design.
/// </summary>
public readonly record struct GuideRnaParameters(
    double MinGcContent,
    double MaxGcContent,
    double MinScore,
    bool AvoidPolyT,
    bool CheckSelfComplementarity)
{
    /// <summary>
    /// Default parameters for guide RNA design.
    /// </summary>
    public static GuideRnaParameters Default => new(
        MinGcContent: 40,
        MaxGcContent: 70,
        MinScore: 50,
        AvoidPolyT: true,
        CheckSelfComplementarity: true);
}

/// <summary>
/// A guide RNA candidate with quality metrics.
/// </summary>
public sealed record GuideRnaCandidate(
    string Sequence,
    int Position,
    bool IsForwardStrand,
    double GcContent,
    double SeedGcContent,
    bool HasPolyT,
    double SelfComplementarityScore,
    double Score,
    IReadOnlyList<string> Issues,
    CrisprSystem System)
{
    /// <summary>
    /// Gets the guide RNA sequence with the standard scaffold.
    /// </summary>
    public string FullGuideRna => Sequence + "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC";
}

/// <summary>
/// Represents a potential off-target site.
/// </summary>
public sealed record OffTargetSite(
    int Position,
    string Sequence,
    int Mismatches,
    IReadOnlyList<int> MismatchPositions,
    bool IsForwardStrand,
    double OffTargetScore);
