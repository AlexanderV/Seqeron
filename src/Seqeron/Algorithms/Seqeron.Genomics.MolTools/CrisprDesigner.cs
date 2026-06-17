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

        // Check seed region (last 10 bp for Cas9)
        // Evidence: Addgene - seed sequence is 8-10 bases at 3' end; using upper bound (10bp)
        string seedRegion = system.PamAfterTarget
            ? seq.Substring(Math.Max(0, seq.Length - 10))
            : seq.Substring(0, Math.Min(10, seq.Length));
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
        if (seedGc < 30 || seedGc > 80)
        {
            score -= 5;
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

        if (guideSequence.Length != system.GuideLength)
            throw new ArgumentException(
                $"Guide length ({guideSequence.Length}) does not match expected length ({system.GuideLength}) for {system.Name}.",
                nameof(guideSequence));

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
        // Seed region: PAM-proximal 12bp (Hsu et al. 2013: 8-12bp; using 12bp as conservative upper bound)
        int seedStart = system.PamAfterTarget ? guide.Length - 12 : 0;
        int seedEnd = system.PamAfterTarget ? guide.Length : 12;

        for (int i = 0; i < Math.Min(guide.Length, target.Length); i++)
        {
            if (guide[i] != target[i])
            {
                // Seed-region (PAM-proximal) mismatches are LESS tolerated by Cas9
                // (Hsu et al. 2013), so they receive a higher penalty here.
                bool inSeed = i >= seedStart && i < seedEnd;
                score += inSeed ? 5 : 2;
            }
        }

        return score;
    }

    #endregion

    #region Doench 2014 On-Target Score (Rule Set 1)

    // ---------------------------------------------------------------------------------------------
    // Doench et al. 2014 "Rule Set 1" on-target efficacy model (PUBLISHED, fully reproducible).
    //
    // Source: Doench, Hartenian, Graham, et al. "Rational design of highly active sgRNAs for
    //   CRISPR-Cas9-mediated gene inactivation." Nat Biotechnol 32, 1262-1267 (2014).
    //   PMID 25184501. doi:10.1038/nbt.3026.
    //
    // Coefficients transcribed verbatim from the reference implementation distributed with the
    //   CRISPOR tool (Haeussler et al. 2016, Genome Biol 17:148), file `doenchScore.py`:
    //   https://github.com/maximilianh/crisporWebsite/blob/master/doenchScore.py
    //   (this is the original on_target model published by Doench 2014; the same coefficient set
    //   used by Azimuth's `model_comparison` baseline and by CRISPOR/CRISPRscan).
    //
    // The model is a logistic-regression linear model over a fixed 30-nt context window:
    //     [4 nt upstream] + [20 nt protospacer] + [3 nt PAM (NGG)] + [3 nt downstream] = 30 nt.
    // It is the dot product of the feature vector with the published weights, passed through a
    // logistic sigmoid. Score in (0, 1); higher = more active. We expose it on a 0-100 scale.
    //
    // The position indices in the weight table are 0-based offsets directly into the 30-mer
    // (this matches the reference: subSeq = seq[pos : pos+len(modelSeq)]). The GC-count term
    // is computed over the 20-nt protospacer = seq[4:24].
    //
    // Cross-checks reproduced this session (independent Python run of the reference coefficients):
    //   "TATAGCTGCGATCTGAGGTAGGGAGGGACC" -> 0.7130893 (reference example value 0.713089368437)
    //   "TCCGCACCTGTCACGGTCGGGGCTTGGCGC" -> 0.0189838 (reference example value 0.0189838463593)
    // ---------------------------------------------------------------------------------------------

    private const double DoenchIntercept = 0.59763615;
    private const double DoenchGcHigh = -0.1665878;
    private const double DoenchGcLow = -0.2026259;

    /// <summary>
    /// (0-based 30-mer offset, sub-sequence, weight) feature table for the Doench 2014 Rule Set 1
    /// linear model. Single-nucleotide features (length-1) and dinucleotide features (length-2).
    /// Transcribed verbatim from the published reference implementation (see region comment).
    /// </summary>
    private static readonly (int Pos, string Seq, double Weight)[] DoenchParams =
    {
        (1, "G", -0.2753771), (2, "A", -0.3238875), (2, "C", 0.17212887), (3, "C", -0.1006662),
        (4, "C", -0.2018029), (4, "G", 0.24595663), (5, "A", 0.03644004), (5, "C", 0.09837684),
        (6, "C", -0.7411813), (6, "G", -0.3932644), (11, "A", -0.466099), (14, "A", 0.08537695),
        (14, "C", -0.013814), (15, "A", 0.27262051), (15, "C", -0.1190226), (15, "T", -0.2859442),
        (16, "A", 0.09745459), (16, "G", -0.1755462), (17, "C", -0.3457955), (17, "G", -0.6780964),
        (18, "A", 0.22508903), (18, "C", -0.5077941), (19, "G", -0.4173736), (19, "T", -0.054307),
        (20, "G", 0.37989937), (20, "T", -0.0907126), (21, "C", 0.05782332), (21, "T", -0.5305673),
        (22, "T", -0.8770074), (23, "C", -0.8762358), (23, "G", 0.27891626), (23, "T", -0.4031022),
        (24, "A", -0.0773007), (24, "C", 0.28793562), (24, "T", -0.2216372), (27, "G", -0.6890167),
        (27, "T", 0.11787758), (28, "C", -0.1604453), (29, "G", 0.38634258), (1, "GT", -0.6257787),
        (4, "GC", 0.30004332), (5, "AA", -0.8348362), (5, "TA", 0.76062777), (6, "GG", -0.4908167),
        (11, "GG", -1.5169074), (11, "TA", 0.7092612), (11, "TC", 0.49629861), (11, "TT", -0.5868739),
        (12, "GG", -0.3345637), (13, "GA", 0.76384993), (13, "GC", -0.5370252), (16, "TG", -0.7981461),
        (18, "GG", -0.6668087), (18, "TC", 0.35318325), (19, "CC", 0.74807209), (19, "TG", -0.3672668),
        (20, "AC", 0.56820913), (20, "CG", 0.32907207), (20, "GA", -0.8364568), (20, "GG", -0.7822076),
        (21, "TC", -1.029693), (22, "CG", 0.85619782), (22, "CT", -0.4632077), (23, "AA", -0.5794924),
        (23, "AG", 0.64907554), (24, "AG", -0.0773007), (24, "CG", 0.28793562), (24, "TG", -0.2216372),
        (26, "GT", 0.11787758), (28, "GG", -0.69774)
    };

    /// <summary>
    /// Calculates the Doench et al. 2014 "Rule Set 1" on-target efficacy score for an SpCas9 guide,
    /// expressed on a 0-100 scale (higher = predicted more active).
    /// </summary>
    /// <param name="context30Mer">
    /// The 30-nt sequence context required by the model:
    /// 4 nt upstream + 20 nt protospacer + 3 nt PAM (must be N<c>GG</c>) + 3 nt downstream.
    /// Case-insensitive; must contain only A/C/G/T.
    /// </param>
    /// <returns>The predicted on-target activity in the range [0, 100].</returns>
    /// <remarks>
    /// This is the published, exactly-reproducible linear model (logistic regression).
    /// It is NOT Doench "Rule Set 2" / Azimuth, which is a gradient-boosted-tree model that cannot
    /// be reproduced from published coefficients without the trained model file.
    /// Source: Doench et al. 2014, Nat Biotechnol 32:1262 (PMID 25184501); coefficients per the
    /// reference implementation in CRISPOR's <c>doenchScore.py</c>.
    /// </remarks>
    public static double CalculateOnTargetDoench2014(string context30Mer)
    {
        if (string.IsNullOrEmpty(context30Mer))
            throw new ArgumentNullException(nameof(context30Mer));

        var seq = context30Mer.ToUpperInvariant();

        if (seq.Length != 30)
            throw new ArgumentException(
                $"Doench 2014 requires a 30-nt context (4 upstream + 20 protospacer + 3 PAM + 3 downstream); got {seq.Length}.",
                nameof(context30Mer));

        foreach (var c in seq)
        {
            if (c is not ('A' or 'C' or 'G' or 'T'))
                throw new ArgumentException(
                    $"Doench 2014 context must contain only A/C/G/T; found '{c}'.",
                    nameof(context30Mer));
        }

        // PAM is at offsets 25-26 (the GG of NGG) within the 30-mer.
        if (seq[25] != 'G' || seq[26] != 'G')
            throw new ArgumentException(
                "Doench 2014 expects an SpCas9 NGG PAM at offsets 25-26 of the 30-nt context.",
                nameof(context30Mer));

        double score = DoenchIntercept;

        // GC-count term over the 20-nt protospacer = offsets [4, 24).
        int gcCount = 0;
        for (int i = 4; i < 24; i++)
        {
            if (seq[i] is 'G' or 'C')
                gcCount++;
        }
        double gcWeight = gcCount <= 10 ? DoenchGcLow : DoenchGcHigh;
        score += Math.Abs(10 - gcCount) * gcWeight;

        // Position-specific single- and di-nucleotide features.
        foreach (var (pos, modelSeq, weight) in DoenchParams)
        {
            if (string.CompareOrdinal(seq, pos, modelSeq, 0, modelSeq.Length) == 0)
                score += weight;
        }

        double probability = 1.0 / (1.0 + Math.Exp(-score));
        return probability * 100.0;
    }

    #endregion

    #region MIT / Hsu 2013 Off-Target Score

    // ---------------------------------------------------------------------------------------------
    // MIT / Hsu 2013 off-target specificity score (PUBLISHED, fully reproducible).
    //
    // Source: Hsu, Scott, Weinstein, et al. "DNA targeting specificity of RNA-guided Cas9
    //   nucleases." Nat Biotechnol 31, 827-832 (2013). PMID 23873081. doi:10.1038/nbt.2647.
    //   Scoring scheme as published on the (now-retired) crispr.mit.edu "about" page and
    //   transcribed in CRISPOR's `crispor.py` (calcHitScore / calcMitGuideScore).
    //   https://github.com/maximilianh/crisporWebsite/blob/master/crispor.py
    //
    // Single-hit score for one candidate off-target with mismatches relative to the 20-nt guide:
    //     score1 = product over mismatched positions i of (1 - W[i])
    //     score2 = 1 / ( ((19 - meanPairwiseMismatchDistance) / 19) * 4 + 1 )   [only if >=2 mm]
    //     score3 = 1 / (nmm^2)                                                   [only if >=1 mm]
    //     hitScore = score1 * score2 * score3 * 100
    //   where W is the published 20-position mismatch-penalty weight vector (PAM-distal -> proximal),
    //   meanPairwiseMismatchDistance is the mean of consecutive inter-mismatch gaps, and nmm is the
    //   number of mismatches. (score2/score3 special-cased to 1 below their thresholds, matching the
    //   reference; a 0-mismatch on-target therefore scores exactly 100.)
    //
    // Aggregate (guide-level) MIT specificity:
    //     guideScore = 100 / (100 + sum of all single-hit scores) * 100
    //   i.e. 100 with no off-targets; decreasing as off-target hits accumulate.
    //
    // Published W vector (20 positions, index 0 = PAM-distal .. index 19 = PAM-proximal):
    //   [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851,
    //    0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
    //
    // Cross-checks reproduced this session (independent Python run of the reference):
    //   perfect match -> 100; single mm @pos19 -> 100*(1-0.583)=41.7; @pos5 -> 60.5;
    //   mm @{5,15} -> score1=0.10406, score2=0.34545.., score3=0.25 -> 0.8987;
    //   aggregate of one 60.5 hit -> 100/(100+60.5)*100 = 62.305296.
    // ---------------------------------------------------------------------------------------------

    /// <summary>
    /// The published 20-position MIT/Hsu 2013 single-nucleotide mismatch-penalty weight vector.
    /// Index 0 is the PAM-distal (5') end of the 20-nt protospacer; index 19 is PAM-proximal (3').
    /// </summary>
    private static readonly double[] MitHitScoreWeights =
    {
        0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079,
        0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583
    };

    /// <summary>
    /// Calculates the MIT / Hsu 2013 single-hit off-target score for a candidate off-target site
    /// against a 20-nt guide. Returns a value in [0, 100]; 100 means an exact match (on-target),
    /// lower values mean the site is a weaker (less likely cut) off-target.
    /// </summary>
    /// <param name="guide20">The 20-nt guide/protospacer (PAM-distal first), A/C/G/T only.</param>
    /// <param name="offTarget20">The 20-nt candidate off-target protospacer, same length/orientation.</param>
    /// <returns>The single-hit MIT/Hsu score in [0, 100].</returns>
    /// <remarks>
    /// Source: Hsu et al. 2013, Nat Biotechnol 31:827 (PMID 23873081); formula per crispr.mit.edu
    /// as transcribed in CRISPOR's <c>calcHitScore</c>.
    /// </remarks>
    public static double CalculateMitHitScore(string guide20, string offTarget20)
    {
        if (string.IsNullOrEmpty(guide20))
            throw new ArgumentNullException(nameof(guide20));
        if (string.IsNullOrEmpty(offTarget20))
            throw new ArgumentNullException(nameof(offTarget20));
        if (guide20.Length != 20 || offTarget20.Length != 20)
            throw new ArgumentException("MIT/Hsu score requires two 20-nt sequences.");

        var g = guide20.ToUpperInvariant();
        var o = offTarget20.ToUpperInvariant();

        const int maxDist = 19;

        var dists = new List<int>();
        int mmCount = 0;
        int lastMmPos = -1;
        double score1 = 1.0;

        for (int pos = 0; pos < 20; pos++)
        {
            if (g[pos] != o[pos])
            {
                mmCount++;
                if (lastMmPos != -1)
                    dists.Add(pos - lastMmPos);
                score1 *= 1.0 - MitHitScoreWeights[pos];
                lastMmPos = pos;
            }
        }

        double score2;
        if (mmCount < 2)
        {
            score2 = 1.0;
        }
        else
        {
            double avgDist = (double)dists.Sum() / dists.Count;
            score2 = 1.0 / (((maxDist - avgDist) / maxDist) * 4.0 + 1.0);
        }

        double score3 = mmCount == 0 ? 1.0 : 1.0 / (mmCount * (double)mmCount);

        return score1 * score2 * score3 * 100.0;
    }

    /// <summary>
    /// Calculates the aggregate MIT / Hsu 2013 guide specificity score from a set of single-hit
    /// off-target scores: <c>100 / (100 + Σ hitScores) × 100</c>. Returns 100 when there are no
    /// off-target hits and decreases toward 0 as off-target burden grows.
    /// </summary>
    /// <param name="singleHitScores">The single-hit MIT/Hsu scores of all candidate off-targets.</param>
    /// <returns>The aggregate guide specificity in [0, 100].</returns>
    /// <remarks>
    /// Source: Hsu et al. 2013; aggregate "Guide score" per crispr.mit.edu (CRISPOR
    /// <c>calcMitGuideScore</c>). The reference rounds to an integer; this method returns the
    /// unrounded value so callers may round as they wish.
    /// </remarks>
    public static double CalculateMitSpecificityScore(IEnumerable<double> singleHitScores)
    {
        ArgumentNullException.ThrowIfNull(singleHitScores);
        double sum = singleHitScores.Sum();
        return 100.0 / (100.0 + sum) * 100.0;
    }

    /// <summary>
    /// Calculates the aggregate MIT / Hsu 2013 specificity score for a guide against a genome by
    /// enumerating off-target sites (via <see cref="FindOffTargets"/>) and combining their MIT/Hsu
    /// single-hit scores. Exact on-target matches are excluded (they are not off-targets).
    /// </summary>
    /// <param name="guideSequence">The 20-nt SpCas9 guide sequence.</param>
    /// <param name="genome">The genome/sequence to scan for off-targets.</param>
    /// <param name="maxMismatches">Maximum mismatches to consider (1-5).</param>
    /// <param name="systemType">CRISPR system (SpCas9 by default).</param>
    /// <returns>The aggregate MIT/Hsu specificity in [0, 100]; 100 when no off-targets are found.</returns>
    public static double CalculateMitSpecificityScore(
        string guideSequence,
        DnaSequence genome,
        int maxMismatches = 4,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
    {
        var offTargets = FindOffTargets(guideSequence, genome, maxMismatches, systemType).ToList();
        if (offTargets.Count == 0)
            return 100.0;

        var guide = guideSequence.ToUpperInvariant();
        var hitScores = offTargets
            .Where(ot => ot.Sequence.Length == 20 && guide.Length == 20)
            .Select(ot => CalculateMitHitScore(guide, ot.Sequence.ToUpperInvariant()));

        return CalculateMitSpecificityScore(hitScores);
    }

    #endregion

    #region CFD (Cutting Frequency Determination) Off-Target Score — Doench 2016

    // ---------------------------------------------------------------------------------------------
    // CFD (Cutting Frequency Determination) off-target score (PUBLISHED, fully reproducible from the
    // shipped matrices).
    //
    // Source: Doench, Fusi, Sullender, Hegde, et al. "Optimized sgRNA design to maximize activity and
    //   minimize off-target effects of CRISPR-Cas9." Nat Biotechnol 34, 184-191 (2016).
    //   PMID 26780180. doi:10.1038/nbt.3437. CFD is defined in the paper and Supplementary Table 19
    //   (per-position, per-mismatch-type percent-activity matrix) plus a PAM-activity table.
    //
    // Canonical reference implementation: `cfd-score-calculator.py` distributed by John Doench's lab
    //   and redistributed in CRISPOR (maximilianh/crisporWebsite, CFD_Scoring/) and bm2-lab/iGWOS
    //   (CFD/). The matrices are shipped as the Python pickles `mismatch_score.pkl` (240 entries =
    //   12 mismatch types x 20 positions) and `pam_scores.pkl` (16 NGG-region dinucleotides).
    //
    // The matrices below were obtained this session by decoding those authoritative pickles to text
    // and were cross-checked element-by-element across TWO independent repositories
    // (maximilianh/crisporWebsite and bm2-lab/iGWOS): 240/240 mismatch values and 16/16 PAM values
    // are identical between the two sources (zero diffs). They are reproduced here verbatim at full
    // double precision (the exact bit patterns decoded from the pickle).
    //
    // ALGORITHM (verbatim from `calc_cfd` in cfd-score-calculator.py):
    //   score = 1
    //   for i, off_base in enumerate(offTarget20):           # i = 0..19, 5' -> 3'
    //       if guide[i] == off_base: score *= 1
    //       else:  key = 'r' + guide[i]            (guide base, T written as U)
    //                  + ':d' + complement(off_base)  (base on the non-target DNA strand)
    //                  + ',' + str(i + 1)            (1-based position, 1 = position 0)
    //              score *= mm_scores[key]
    //   score *= pam_scores[ pam[-2:] ]              # last two PAM nt; PAM position 1 (N) -> 1
    //
    // ORIENTATION (pinned from the source, NOT from this code): the loop enumerates the 20-nt
    //   protospacer 5'->3', so position 1 (index 0) is the 5' / PAM-DISTAL end and position 20 is the
    //   3' / PAM-PROXIMAL (seed) end. (`sg = off[:-3]`, `pam = off[-2:]` => the protospacer precedes
    //   the PAM, so the high index is adjacent to the PAM.) Getting this backwards is the classic CFD
    //   bug; the orientation guard test asserts rC:dT,1 = 1.0 vs rC:dT,20 = 0.5 and fails on reversal.
    //
    // KEY CONVENTION (pinned from the source): 'rX' is the GUIDE (RNA) base with T written as U; 'dY'
    //   is the COMPLEMENT of the off-target protospacer base (i.e. the base on the off-target's
    //   non-target DNA strand that pairs against the guide). A position only contributes a penalty
    //   when guide[i] != offTarget[i]; matched positions contribute 1.0. Perfect match + GG PAM -> 1.0.
    //
    // CONTRACT: this implementation scores a 20-nt guide vs a 20-nt off-target protospacer plus the
    //   off-target's PAM, all A/C/G/T only (case-insensitive). The PAM score uses the last two PAM
    //   nucleotides; a 2-nt PAM string is taken as-is, a 3-nt PAM uses its last two. Insertions /
    //   deletions and non-ACGT bases are NOT supported (CFD is undefined for them) and throw.
    //
    // Cross-checks reproduced this session (independent Python re-derivation from the decoded pickle,
    // NOT from the C# arrays; the first is the published iGWOS doctest oracle):
    //   guide "GGGGGGGGGGGGGGGGGGGG" vs off "GGGGGGGGGGGGGGGGGAAA" + GG  -> 0.4635989007074176
    //     (= rG:dT,18 x rG:dT,19 x rG:dT,20 = 0.692307692 x 0.714285714 x 0.9375)
    //   perfect match + GG -> 1.0 ; perfect match + GA -> 0.069444 ; perfect + AG -> 0.259259
    //   guide "GACGCATAAAGATGAGACGC": off A@pos1 (rG:dT,1) -> 0.9 ; off A@pos20 (rC:dT,20) -> 0.5 ;
    //     off A@{1,20} -> 0.45 (product) ; off T@pos16 (rG:dA,16) -> 0.0
    // ---------------------------------------------------------------------------------------------

    /// <summary>
    /// CFD per-position mismatch percent-activity matrix (Doench 2016), keyed by mismatch type
    /// <c>"rX:dY"</c> (X = guide/RNA base with T as U; Y = complement of the off-target base), each
    /// value an array indexed 0..19 by 0-based protospacer position (index 0 = 5'/PAM-distal,
    /// index 19 = 3'/PAM-proximal). Verbatim from the authoritative <c>mismatch_score.pkl</c>.
    /// </summary>
    private static readonly Dictionary<string, double[]> CfdMismatchScores = new()
    {
        ["rA:dA"] = new[] { 1.0, 0.727272727, 0.705882353, 0.636363636, 0.363636364, 0.7142857140000001, 0.4375, 0.428571429, 0.6, 0.882352941, 0.307692308, 0.333333333, 0.3, 0.533333333, 0.2, 0.0, 0.133333333, 0.5, 0.538461538, 0.6 },
        ["rA:dC"] = new[] { 1.0, 0.8, 0.611111111, 0.625, 0.72, 0.7142857140000001, 0.705882353, 0.7333333329999999, 0.666666667, 0.5555555560000001, 0.65, 0.7222222220000001, 0.6521739129999999, 0.46666666700000003, 0.65, 0.192307692, 0.176470588, 0.4, 0.375, 0.764705882 },
        ["rA:dG"] = new[] { 0.857142857, 0.7857142859999999, 0.428571429, 0.352941176, 0.5, 0.454545455, 0.4375, 0.428571429, 0.571428571, 0.333333333, 0.4, 0.263157895, 0.21052631600000002, 0.214285714, 0.272727273, 0.0, 0.176470588, 0.19047619, 0.20689655199999998, 0.22727272699999998 },
        ["rC:dA"] = new[] { 1.0, 0.9090909090000001, 0.6875, 0.8, 0.636363636, 0.9285714290000001, 0.8125, 0.875, 0.875, 0.9411764709999999, 0.307692308, 0.538461538, 0.7, 0.7333333329999999, 0.066666667, 0.307692308, 0.46666666700000003, 0.642857143, 0.46153846200000004, 0.3 },
        ["rC:dC"] = new[] { 0.913043478, 0.695652174, 0.5, 0.5, 0.6, 0.5, 0.470588235, 0.642857143, 0.6190476189999999, 0.38888888899999996, 0.25, 0.444444444, 0.13636363599999998, 0.0, 0.05, 0.153846154, 0.058823529000000006, 0.133333333, 0.125, 0.058823529000000006 },
        ["rC:dT"] = new[] { 1.0, 0.727272727, 0.8666666670000001, 0.842105263, 0.571428571, 0.9285714290000001, 0.75, 0.65, 0.857142857, 0.8666666670000001, 0.75, 0.7142857140000001, 0.384615385, 0.35, 0.222222222, 1.0, 0.46666666700000003, 0.538461538, 0.428571429, 0.5 },
        ["rG:dA"] = new[] { 1.0, 0.636363636, 0.5, 0.363636364, 0.3, 0.666666667, 0.571428571, 0.625, 0.533333333, 0.8125, 0.384615385, 0.384615385, 0.3, 0.26666666699999997, 0.14285714300000002, 0.0, 0.25, 0.666666667, 0.666666667, 0.7 },
        ["rG:dG"] = new[] { 0.7142857140000001, 0.692307692, 0.384615385, 0.529411765, 0.7857142859999999, 0.681818182, 0.6875, 0.615384615, 0.538461538, 0.4, 0.428571429, 0.529411765, 0.42105263200000004, 0.428571429, 0.272727273, 0.0, 0.235294118, 0.47619047600000003, 0.448275862, 0.428571429 },
        ["rG:dT"] = new[] { 0.9, 0.846153846, 0.75, 0.9, 0.8666666670000001, 1.0, 1.0, 1.0, 0.642857143, 0.933333333, 1.0, 0.933333333, 0.923076923, 0.75, 0.9411764709999999, 1.0, 0.933333333, 0.692307692, 0.7142857140000001, 0.9375 },
        ["rU:dC"] = new[] { 0.956521739, 0.84, 0.5, 0.625, 0.64, 0.571428571, 0.588235294, 0.7333333329999999, 0.6190476189999999, 0.5, 0.4, 0.5, 0.260869565, 0.0, 0.05, 0.346153846, 0.117647059, 0.333333333, 0.25, 0.176470588 },
        ["rU:dG"] = new[] { 0.857142857, 0.857142857, 0.428571429, 0.647058824, 1.0, 0.9090909090000001, 0.6875, 1.0, 0.923076923, 0.533333333, 0.666666667, 0.947368421, 0.7894736840000001, 0.28571428600000004, 0.272727273, 0.666666667, 0.705882353, 0.428571429, 0.275862069, 0.090909091 },
        ["rU:dT"] = new[] { 1.0, 0.846153846, 0.7142857140000001, 0.47619047600000003, 0.5, 0.8666666670000001, 0.875, 0.8, 0.9285714290000001, 0.857142857, 0.75, 0.8, 0.692307692, 0.6190476189999999, 0.578947368, 0.9090909090000001, 0.533333333, 0.666666667, 0.28571428600000004, 0.5625 }
    };

    /// <summary>
    /// CFD PAM-activity table (Doench 2016) for the last two PAM nucleotides (the "GG" region of NGG).
    /// Verbatim from the authoritative <c>pam_scores.pkl</c>; canonical GG = 1.0. Position 1 of the
    /// PAM is N and contributes 1 (not part of this table).
    /// </summary>
    private static readonly Dictionary<string, double> CfdPamScores = new()
    {
        ["AA"] = 0.0,
        ["AC"] = 0.0,
        ["AG"] = 0.25925925899999996,
        ["AT"] = 0.0,
        ["CA"] = 0.0,
        ["CC"] = 0.0,
        ["CG"] = 0.107142857,
        ["CT"] = 0.0,
        ["GA"] = 0.06944444400000001,
        ["GC"] = 0.022222222000000003,
        ["GG"] = 1.0,
        ["GT"] = 0.016129031999999998,
        ["TA"] = 0.0,
        ["TC"] = 0.0,
        ["TG"] = 0.038961038999999996,
        ["TT"] = 0.0
    };

    private static char CfdComplement(char b) => b switch
    {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => throw new ArgumentException($"CFD requires A/C/G/T sequences; found '{b}'.")
    };

    /// <summary>
    /// Calculates the Doench 2016 CFD (Cutting Frequency Determination) off-target score for a 20-nt
    /// guide RNA against a 20-nt off-target protospacer with the off-target's PAM. The score is the
    /// product over the 20 protospacer positions of the per-position mismatch penalty (1.0 where the
    /// guide:off-target base pair matches) times the PAM-activity score for the off-target's PAM.
    /// A perfect match against a canonical NGG PAM scores exactly 1.0; weaker off-targets score lower.
    /// </summary>
    /// <param name="sgRna20">
    /// The 20-nt guide/protospacer, 5' (PAM-distal) first. A/C/G/T only, case-insensitive (T is
    /// treated as the RNA base U internally, per the model).
    /// </param>
    /// <param name="offTarget20">
    /// The 20-nt candidate off-target protospacer, same length and orientation as <paramref name="sgRna20"/>.
    /// </param>
    /// <param name="offTargetPam">
    /// The off-target's PAM. Only the last two nucleotides are scored (the N of NGG always contributes 1);
    /// pass either a 2-nt PAM (e.g. <c>"GG"</c>) or a 3-nt PAM (e.g. <c>"AGG"</c>). A/C/G/T only.
    /// </param>
    /// <returns>The CFD off-target score in [0, 1]; higher = more likely an active off-target.</returns>
    /// <remarks>
    /// Source: Doench et al. 2016, Nat Biotechnol 34:184 (PMID 26780180); matrices from the
    /// authoritative <c>mismatch_score.pkl</c>/<c>pam_scores.pkl</c> (cross-checked across CRISPOR and
    /// iGWOS). Insertions/deletions and non-ACGT bases are unsupported and throw; CFD is defined only
    /// for a 20-nt:20-nt protospacer comparison plus the off-target PAM.
    /// </remarks>
    public static double CalculateCfdScore(string sgRna20, string offTarget20, string offTargetPam)
    {
        if (string.IsNullOrEmpty(sgRna20))
            throw new ArgumentNullException(nameof(sgRna20));
        if (string.IsNullOrEmpty(offTarget20))
            throw new ArgumentNullException(nameof(offTarget20));
        if (string.IsNullOrEmpty(offTargetPam))
            throw new ArgumentNullException(nameof(offTargetPam));

        if (sgRna20.Length != 20)
            throw new ArgumentException($"CFD requires a 20-nt guide; got {sgRna20.Length}.", nameof(sgRna20));
        if (offTarget20.Length != 20)
            throw new ArgumentException($"CFD requires a 20-nt off-target protospacer; got {offTarget20.Length}.", nameof(offTarget20));

        var guide = sgRna20.ToUpperInvariant();
        var off = offTarget20.ToUpperInvariant();
        var pam = offTargetPam.ToUpperInvariant();

        // PAM: score the last two nucleotides (the N of NGG always contributes 1).
        if (pam.Length is not (2 or 3))
            throw new ArgumentException($"CFD PAM must be 2 or 3 nt; got {pam.Length}.", nameof(offTargetPam));
        string pamKey = pam.Length == 3 ? pam.Substring(1) : pam;
        if (!CfdPamScores.TryGetValue(pamKey, out double pamScore))
            throw new ArgumentException($"CFD PAM '{offTargetPam}' must contain only A/C/G/T.", nameof(offTargetPam));

        double score = 1.0;

        for (int i = 0; i < 20; i++)
        {
            char g = guide[i];
            char o = off[i];

            if (g is not ('A' or 'C' or 'G' or 'T'))
                throw new ArgumentException($"CFD requires A/C/G/T sequences; found '{g}'.", nameof(sgRna20));
            if (o is not ('A' or 'C' or 'G' or 'T'))
                throw new ArgumentException($"CFD requires A/C/G/T sequences; found '{o}'.", nameof(offTarget20));

            if (g == o)
                continue; // matched position contributes 1.0

            // rX = guide base (T written as U); dY = complement of the off-target base.
            char rBase = g == 'T' ? 'U' : g;
            char dBase = CfdComplement(o);
            string key = $"r{rBase}:d{dBase}";

            // Position is 1-based 1..20 (= index i+1); index 0 = 5'/PAM-distal end.
            score *= CfdMismatchScores[key][i];
        }

        return score * pamScore;
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
/// <param name="Position">PAM start coordinate, always expressed on the forward strand.</param>
/// <param name="PamSequence">The matched PAM sequence (forward-strand orientation).</param>
/// <param name="TargetSequence">The guide/protospacer sequence sliced from the strand on which the PAM was found.</param>
/// <param name="TargetStart">
/// Start index of <see cref="TargetSequence"/> on the strand the hit was found on.
/// For forward-strand hits (<see cref="IsForwardStrand"/> == true) this is a forward-strand
/// index. For reverse-strand hits it is an index into the reverse-complement string (used to
/// slice <see cref="TargetSequence"/>), NOT a forward-strand coordinate — unlike
/// <see cref="Position"/>, which is always forward-strand.
/// </param>
/// <param name="IsForwardStrand">True if the PAM was found on the forward strand; false for the reverse strand.</param>
/// <param name="System">The CRISPR system whose PAM/guide-length parameters produced this site.</param>
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
