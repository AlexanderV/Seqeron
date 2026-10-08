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
                    // Convert position back to forward strand coordinates.
                    // Convention (CRISPOR crispor.py findAllPams/flankSeqIter): coordinates are
                    // always forward-strand, sequences are always read on the protospacer strand.
                    int forwardPos = seq.Length - i - pamPattern.Length;

                    // Forward-strand start of the protospacer. On the reverse strand the guide
                    // lies on the opposite side of the PAM in forward coordinates:
                    //   Cas9 (PAM 3' of guide)   -> guide starts just after the PAM;
                    //   Cas12a (PAM 5' of guide) -> guide starts guideLength bases before it.
                    int forwardTargetStart = pamAfterTarget
                        ? forwardPos + pamPattern.Length
                        : forwardPos - guideLength;

                    string target = revComp.Substring(revTargetStart, guideLength);

                    yield return new PamSite(
                        Position: forwardPos,
                        PamSequence: revComp.Substring(i, pamPattern.Length),
                        TargetSequence: target,
                        TargetStart: forwardTargetStart,
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
    /// Designs guide RNAs for a target region, ranked best-first.
    /// </summary>
    /// <param name="sequence">DNA sequence containing the target region.</param>
    /// <param name="regionStart">Start of the region to target (0-based).</param>
    /// <param name="regionEnd">End of the region to target (0-based, inclusive).</param>
    /// <param name="systemType">Type of CRISPR system.</param>
    /// <param name="parameters">Optional design parameters.</param>
    /// <returns>
    /// Guide RNA candidates whose predicted cleavage site (<see cref="GetCutSite"/>) falls inside
    /// <paramref name="regionStart"/>..<paramref name="regionEnd"/> and whose
    /// <see cref="GuideRnaCandidate.Score"/> reaches <see cref="GuideRnaParameters.MinScore"/>,
    /// ordered by <see cref="GuideRnaParameters.Ranking"/> (best first) and then by position /
    /// forward-strand-first for a deterministic order.
    /// </returns>
    /// <remarks>
    /// Ranking follows the reference tool CRISPOR (<c>crispor.py</c>, <c>mergeGuideInfo</c>:
    /// <c>guideData.sort(reverse=True, key=…)</c> — descending by the selected score column).
    /// </remarks>
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
        return DesignGuideRnasCore(sequence, regionStart, regionEnd, systemType, parameters);
    }

    private static IEnumerable<GuideRnaCandidate> DesignGuideRnasCore(
        DnaSequence sequence, int regionStart, int regionEnd, CrisprSystemType systemType, GuideRnaParameters? parameters)
    {
        var effectiveParams = parameters ?? GuideRnaParameters.Default;
        var system = GetSystem(systemType);

        var candidates = FindPamSitesCore(sequence.Sequence, system)
            .Where(p => IsInRegion(p, regionStart, regionEnd))
            .Select(p => EvaluateGuideRna(p, sequence.Sequence, effectiveParams, system))
            .Where(c => c.Score >= effectiveParams.MinScore)
            .Where(c => effectiveParams.MaxSelfComplementaryStems is not int max
                        || c.SelfComplementaryStems <= max);

        // CRISPOR sorts the guide table descending by the selected score column; ties are broken
        // deterministically here (position, then forward strand first) so the output is stable.
        var ranked = effectiveParams.Ranking == GuideRnaRanking.OnTargetRuleSet2
            ? candidates.OrderByDescending(c => c.OnTargetScore ?? double.NegativeInfinity)
            : candidates.OrderByDescending(c => c.Score);

        return ranked.ThenBy(c => c.Position).ThenByDescending(c => c.IsForwardStrand);
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

        return EvaluateGuideRnaCore(guideSequence, GetSystem(systemType), parameters ?? GuideRnaParameters.Default);
    }

    private static GuideRnaCandidate EvaluateGuideRnaCore(
        string guideSequence,
        CrisprSystem system,
        GuideRnaParameters effectiveParams)
    {
        var seq = guideSequence.ToUpperInvariant();

        // Calculate GC content
        double gcContent = CalculateGcContent(seq);

        // Check for polyT (transcription terminator)
        bool hasPolyT = HasPolyT(seq, 4);

        // Calculate self-complementarity score
        double selfCompScore = CalculateSelfComplementarity(seq);

        // CHOPCHOP's published self-complementarity measure (4-bp stems); reported always, filtered
        // only when the caller sets a limit (the reference's own default is "report, do not filter").
        int selfCompStems = CountSelfComplementaryStems(seq, 4, effectiveParams.SelfComplementarityBackboneRegions);

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

        // PolyT penalty (opt-out via GuideRnaParameters.AvoidPolyT)
        if (hasPolyT && effectiveParams.AvoidPolyT)
        {
            score -= 20;
            issues.Add("Contains TTTT (potential Pol III terminator)");
        }

        // Self-complementarity penalty (opt-out via GuideRnaParameters.CheckSelfComplementarity)
        if (selfCompScore > 0.3 && effectiveParams.CheckSelfComplementarity)
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

        // Graf et al. 2019 inefficiency motifs. CRISPOR reports these as a warning for NGG systems
        // ("Inefficient") without altering any score, so no deduction is applied here either.
        var grafMotif = system.PamSequence == "NGG" ? GetGrafMotif(seq) : GrafMotifType.None;
        if (grafMotif == GrafMotifType.TtMotif)
            issues.Add("Graf 2019 TT-motif (inefficient in Pol III-based expression)");
        else if (grafMotif == GrafMotifType.GccMotif)
            issues.Add("Graf 2019 GCC-motif (generally inefficient)");

        // CHOPCHOP filters guides above filterSelfCompMax; the limit is opt-in here (null = the
        // reference's default of no filter), and the overrun is reported as an issue.
        if (effectiveParams.MaxSelfComplementaryStems is int maxStems && selfCompStems > maxStems)
            issues.Add($"{selfCompStems} self-complementary 4-bp stems (limit {maxStems})");

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
            System: system)
        {
            GrafMotif = grafMotif,
            SelfComplementaryStems = selfCompStems,
        };
    }

    private static GuideRnaCandidate EvaluateGuideRna(
        PamSite pamSite,
        string fullSequence,
        GuideRnaParameters effectiveParams,
        CrisprSystem system)
    {
        // The site's own system is carried through: a name-based remap would silently replace the
        // SpCas9-NAG / AsCas12a / LbCas12a / CasX definitions (and therefore the seed-region
        // orientation and the reported guide length) with SpCas9's.
        var candidate = EvaluateGuideRnaCore(pamSite.TargetSequence, pamSite.System, effectiveParams);

        string? context30Mer = TryGetRuleSet2Context(pamSite, fullSequence);

        return candidate with
        {
            Position = pamSite.TargetStart,
            IsForwardStrand = pamSite.IsForwardStrand,
            Context30Mer = context30Mer,
            OnTargetScore = context30Mer is null ? null : AzimuthRuleSet2.Score(context30Mer),
        };
    }

    /// <summary>
    /// Returns the forward-strand coordinate of the base immediately 3' of the predicted cleavage
    /// position of <paramref name="pamSite"/> on the PAM-bearing strand, for both strands.
    /// </summary>
    /// <param name="pamSite">A PAM site produced by <see cref="FindPamSites(DnaSequence, CrisprSystemType)"/>.</param>
    /// <returns>A 0-based forward-strand coordinate (may fall outside the sequence for sites at the ends).</returns>
    /// <remarks>
    /// Conventions follow the reference tool CRISPOR (<c>crispor.py</c>):
    /// <list type="bullet">
    /// <item><description>Cas9-type systems (PAM 3' of the protospacer) cut bluntly 3 bp 5' of the PAM
    /// ("the expected cleavage position located -3bp 5' of the PAM site"; the cut marker spans the three
    /// protospacer bases next to the PAM — <c>startFt = start - 3</c> on '+', <c>ftSeq + "---"</c> on '-').
    /// The returned coordinate is the third of those bases: <c>Position - 3</c> (forward) and
    /// <c>Position + pamLength + 2</c> (reverse).</description></item>
    /// <item><description>Cas12a/Cpf1-type systems (PAM 5' of the protospacer) make a staggered cut
    /// "after the 18th base on the non-targeted strand which has the TTTV PAM motif" (Zetsche et al.
    /// 2015, Cell 163:759, Fig. 3, as described by CRISPOR). The returned coordinate is the 19th
    /// protospacer base: <c>Position + pamLength + 18</c> (forward) and <c>Position - 19</c> (reverse).
    /// </description></item>
    /// </list>
    /// </remarks>
    public static int GetCutSite(PamSite pamSite)
    {
        ArgumentNullException.ThrowIfNull(pamSite);
        int pamLength = pamSite.PamSequence.Length;

        if (pamSite.System.PamAfterTarget)
            return pamSite.IsForwardStrand ? pamSite.Position - 3 : pamSite.Position + pamLength + 2;

        return pamSite.IsForwardStrand ? pamSite.Position + pamLength + 18 : pamSite.Position - 19;
    }

    private static bool IsInRegion(PamSite pamSite, int regionStart, int regionEnd)
    {
        int cutSite = GetCutSite(pamSite);
        return cutSite >= regionStart && cutSite <= regionEnd;
    }

    /// <summary>
    /// Extracts the 30-nt Rule Set 2 / Azimuth context of a PAM site from the surrounding sequence:
    /// 4 nt 5' of the protospacer + the 20-nt protospacer + the 3-nt NGG PAM + 3 nt 3' of the PAM,
    /// all read on the protospacer strand (CRISPOR <c>crisporEffScores.calcAllScores</c> builds exactly
    /// this window with <c>trimSeqs(seqs, -24, 6)</c> around the PAM start). Returns <c>null</c> when the
    /// system is not a 20-nt NGG system, the flanks do not fit, or the window is not all A/C/G/T
    /// (the reference scores such windows as "can't do Ns").
    /// </summary>
    private static string? TryGetRuleSet2Context(PamSite pamSite, string fullSequence)
    {
        var system = pamSite.System;
        if (system.PamSequence != "NGG" || system.GuideLength != 20 || !system.PamAfterTarget)
            return null;

        // Forward-strand span of [4 nt upstream | 20 nt guide | 3 nt PAM | 3 nt downstream]:
        // on the forward strand the upstream flank precedes the guide, on the reverse strand it follows it.
        int windowStart = pamSite.IsForwardStrand ? pamSite.TargetStart - 4 : pamSite.TargetStart - 6;
        if (windowStart < 0 || windowStart + 30 > fullSequence.Length)
            return null;

        string window = fullSequence.Substring(windowStart, 30);
        if (!pamSite.IsForwardStrand)
            window = DnaSequence.GetReverseComplementString(window);

        foreach (var c in window)
        {
            if (c is not ('A' or 'C' or 'G' or 'T'))
                return null;
        }

        return window;
    }

    /// <summary>
    /// Classifies a guide against the two inefficiency motifs described by Graf et al. (2019),
    /// Cell Reports 26:1098–1103 — the "TT-motif" and the "GCC-motif" at the 3' (PAM-proximal) end
    /// of an SpCas9 guide. Guides carrying either motif were markedly less efficient.
    /// </summary>
    /// <param name="guideSequence">The protospacer, 5'→3' (A/C/G/T, case-insensitive).</param>
    /// <returns>The motif found, or <see cref="GrafMotifType.None"/>.</returns>
    /// <remarks>
    /// Transcribed from the reference implementation <c>getGrafType</c> in CRISPOR's
    /// <c>crisporEffScores.py</c>: a guide has the TT-motif when it ends in <c>TTC</c> or <c>TTT</c>,
    /// or when its last four bases consist only of T and C with at least two Ts, or when its last four
    /// bases contain <c>TT</c> plus at least three Ts or at least one C; it has the GCC-motif when it
    /// ends in <c>[AGT]GCC</c> or in <c>GCCT</c>. CRISPOR applies the check only to NGG systems.
    /// </remarks>
    public static GrafMotifType GetGrafMotif(string guideSequence)
    {
        if (string.IsNullOrEmpty(guideSequence))
            throw new ArgumentNullException(nameof(guideSequence));

        var seq = guideSequence.ToUpperInvariant();
        if (seq.EndsWith("TTC", StringComparison.Ordinal) || seq.EndsWith("TTT", StringComparison.Ordinal))
            return GrafMotifType.TtMotif;

        string suffix = seq.Length >= 4 ? seq[^4..] : seq;
        int tCount = 0, cCount = 0;
        bool onlyTc = suffix.Length > 0;
        bool hasT = false, hasC = false;
        foreach (var c in suffix)
        {
            if (c == 'T') { tCount++; hasT = true; }
            else if (c == 'C') { cCount++; hasC = true; }
            else onlyTc = false;
        }

        // Reference: set(suffix) == set(["T", "C"]) and suffix.count("T") >= 2
        if (onlyTc && hasT && hasC && tCount >= 2)
            return GrafMotifType.TtMotif;

        // Reference: "TT" in suffix and (suffix.count("T") >= 3 or suffix.count("C") >= 1)
        if (suffix.Contains("TT", StringComparison.Ordinal) && (tCount >= 3 || cCount >= 1))
            return GrafMotifType.TtMotif;

        if (seq.EndsWith("GCC", StringComparison.Ordinal) && suffix.Length == 4 && suffix[0] is 'A' or 'G' or 'T')
            return GrafMotifType.GccMotif;

        if (seq.EndsWith("GCCT", StringComparison.Ordinal))
            return GrafMotifType.GccMotif;

        return GrafMotifType.None;
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
        return FindOffTargetsCore(guideSequence, genome, maxMismatches, systemType);
    }

    private static IEnumerable<OffTargetSite> FindOffTargetsCore(
        string guideSequence, DnaSequence genome, int maxMismatches, CrisprSystemType systemType)
    {
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

        // GC-count term over the 20-nt protospacer = offsets [4, 24); counted with the canonical
        // GC primitive (SequenceExtensions.CountGcAndValidNucleotides) over that fixed window.
        int gcCount = seq.AsSpan(4, 20).CountGcAndValidNucleotides().GcCount;
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

    #region Doench 2016 On-Target Score (Rule Set 2 / Azimuth)

    /// <summary>
    /// Calculates the Doench et al. 2016 "Rule Set 2" / Azimuth on-target efficacy score for an SpCas9
    /// guide, using sequence context only (the standard "Azimuth score" reported by CRISPOR).
    /// </summary>
    /// <param name="context30Mer">
    /// The 30-nt sequence context required by the model:
    /// 4 nt upstream + 20 nt protospacer + 3 nt PAM (must be N<c>GG</c>) + 3 nt downstream.
    /// Case-insensitive; must contain only A/C/G/T.
    /// </param>
    /// <returns>The predicted on-target activity, conventionally in the range [0, 1] (higher = more active).</returns>
    /// <remarks>
    /// Rule Set 2 is a trained gradient-boosted-tree model (not a published linear formula like
    /// <see cref="CalculateOnTargetDoench2014"/>). This uses a faithful, sklearn-free reconstruction of
    /// Microsoft Research's Azimuth <c>V3_model_nopos</c>; the score reproduces
    /// <c>azimuth.model_comparison.predict</c> for that model.
    /// Source: Doench et al. 2016, Nat Biotechnol 34:184 (PMID 26780180); trained model from
    /// https://github.com/MicrosoftResearch/Azimuth (BSD-3-Clause). See
    /// scripts/azimuth/extract_azimuth_model.py for the full provenance and reproduction.
    /// </remarks>
    public static double CalculateOnTargetRuleSet2(string context30Mer)
        => AzimuthRuleSet2.Score(context30Mer);

    /// <summary>
    /// Calculates the Doench et al. 2016 "Rule Set 2" / Azimuth on-target score using the gene-context
    /// (full) model, which additionally incorporates where the cut falls within the target protein.
    /// </summary>
    /// <param name="context30Mer">The 30-nt context (see <see cref="CalculateOnTargetRuleSet2(string)"/>).</param>
    /// <param name="aminoAcidCutPosition">The amino-acid position of the cut site within the target protein.</param>
    /// <param name="percentPeptide">
    /// The cut position as a percentage [0, 100] along the coding sequence (how far into the protein).
    /// </param>
    /// <returns>The predicted on-target activity, conventionally in the range [0, 1].</returns>
    /// <remarks>
    /// Uses Azimuth's <c>V3_model_full</c>. When gene context is unavailable, prefer the sequence-only
    /// <see cref="CalculateOnTargetRuleSet2(string)"/> overload. See the remarks there for provenance.
    /// </remarks>
    public static double CalculateOnTargetRuleSet2(string context30Mer, int aminoAcidCutPosition, double percentPeptide)
        => AzimuthRuleSet2.Score(context30Mer, aminoAcidCutPosition, percentPeptide);

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

    /// <summary>
    /// The standard sgRNA scaffold region CHOPCHOP offers for its self-complementarity check
    /// ("standard backbone" <c>AGGCTAGTCCGT</c>), to be passed to
    /// <see cref="CountSelfComplementaryStems"/> on the same strand as the guide.
    /// </summary>
    public const string StandardSgRnaBackboneRegion = "AGGCTAGTCCGT";

    /// <summary>
    /// Counts the self-complementary stems of a guide RNA: the number of <paramref name="stemLength"/>-mer
    /// windows of the protospacer that have at least 50% GC and can pair either with a downstream part of
    /// the guide itself or with one of the supplied scaffold (backbone) regions. Guide secondary structure
    /// of this kind impedes sgRNA activity.
    /// </summary>
    /// <param name="guideSequence">The protospacer 5'→3' (without the PAM, as in the reference).</param>
    /// <param name="stemLength">Stem length; the reference's <c>STEM_LEN</c> is 4.</param>
    /// <param name="backboneRegions">
    /// Optional scaffold regions written on the same strand as the guide (CHOPCHOP's <c>-BB</c> argument,
    /// e.g. <see cref="StandardSgRnaBackboneRegion"/>); they are reverse-complemented internally, exactly
    /// as CHOPCHOP does before scoring. <c>null</c> or empty = guide-internal stems only, which is the
    /// reference tool's own default.
    /// </param>
    /// <returns>The number of self-complementary stems (0 when the guide is no longer than a stem).</returns>
    /// <remarks>
    /// Transcribed from the reference implementation <c>selfComp</c> / <c>calcSelfComplementarity</c> in
    /// CHOPCHOP's <c>chopchop.py</c> (Labun et al., "CHOPCHOP v2/v3", Nucleic Acids Res 44:W272 (2016) /
    /// 47:W171 (2019)): with <c>rvs = revComp(guide)</c> and <c>L = len − STEM_LEN − 1</c>, every window
    /// <c>guide[i..i+STEM_LEN)</c> whose GC fraction is ≥ 0.5 counts once when it occurs in
    /// <c>rvs[0..L−i)</c> or inside a backbone region. The PAM is excluded ("Do not include PAM motif in
    /// folding calculations"). CHOPCHOP's own CLI default is <c>filterSelfCompMax = -1</c>, i.e. report
    /// without filtering; <see cref="GuideRnaParameters.MaxSelfComplementaryStems"/> is the filter.
    /// </remarks>
    public static int CountSelfComplementaryStems(
        string guideSequence,
        int stemLength = 4,
        IEnumerable<string>? backboneRegions = null)
    {
        if (string.IsNullOrEmpty(guideSequence))
            throw new ArgumentNullException(nameof(guideSequence));
        if (stemLength < 1)
            throw new ArgumentOutOfRangeException(nameof(stemLength), stemLength, "Stem length must be >= 1.");

        var fwd = guideSequence.ToUpperInvariant();
        if (fwd.Length <= stemLength)
            return 0;

        string rvs = DnaSequence.GetReverseComplementString(fwd);
        int l = fwd.Length - stemLength - 1;

        // CHOPCHOP reverse-complements the backbone regions it is given before matching against them.
        string[] backbones = backboneRegions is null
            ? Array.Empty<string>()
            : backboneRegions
                .Select(b => DnaSequence.GetReverseComplementString(b.ToUpperInvariant()))
                .ToArray();

        int folding = 0;
        for (int i = 0; i < fwd.Length - stemLength; i++)
        {
            // gccontent(stem) >= 0.5  <=>  2 * gcCount >= stemLength
            if (fwd.AsSpan(i, stemLength).CountGcAndValidNucleotides().GcCount * 2 < stemLength)
                continue;

            string stem = fwd.Substring(i, stemLength);
            bool pairs = rvs.AsSpan(0, l - i).IndexOf(stem.AsSpan()) >= 0;
            if (!pairs)
            {
                foreach (var backbone in backbones)
                {
                    if (backbone.Contains(stem, StringComparison.Ordinal)) { pairs = true; break; }
                }
            }
            if (pairs)
                folding++;
        }

        return folding;
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

    /// <summary>
    /// Recognition sequences of the common cloning enzymes screened by <see cref="EvaluateGuideRna"/>,
    /// resolved by name from the canonical enzyme table in <see cref="RestrictionAnalyzer.Enzymes"/>
    /// rather than duplicated as literals (EcoRI, BamHI, HindIII, PstI, NotI).
    /// </summary>
    private static readonly string[] CommonCloningSites = ResolveCommonCloningSites();

    private static string[] ResolveCommonCloningSites()
    {
        string[] names = { "EcoRI", "BamHI", "HindIII", "PstI", "NotI" };
        var sites = new string[names.Length];
        for (int i = 0; i < names.Length; i++)
        {
            sites[i] = RestrictionAnalyzer.GetEnzyme(names[i])?.RecognitionSequence
                ?? throw new InvalidOperationException(
                    $"Restriction enzyme '{names[i]}' is missing from the canonical enzyme table.");
        }
        return sites;
    }

    private static bool HasCommonRestrictionSite(string sequence)
    {
        foreach (var site in CommonCloningSites)
        {
            if (sequence.Contains(site, StringComparison.Ordinal))
                return true;
        }
        return false;
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
/// <remarks>
/// Coordinate/orientation convention (matches the CRISPOR reference implementation,
/// <c>crispor.py</c> <c>findAllPams</c> + <c>flankSeqIter</c>): <b>coordinates are always
/// forward-strand, 0-based</b>; <b>sequences are always read on the protospacer strand</b>
/// (the strand the PAM was matched on), so <see cref="PamSequence"/> always satisfies the
/// system's PAM motif under IUPAC matching and <see cref="TargetSequence"/> is the guide as
/// it would be ordered.
/// </remarks>
/// <param name="Position">PAM start coordinate, always expressed on the forward strand (0-based).</param>
/// <param name="PamSequence">
/// The matched PAM, read 5'→3' on the strand it was found on. For reverse-strand hits this is
/// the reverse complement of the forward-strand bases at <see cref="Position"/> (e.g. an SpCas9
/// hit reported as <c>TGG</c> reads <c>CCA</c> on the forward strand).
/// </param>
/// <param name="TargetSequence">The guide/protospacer sequence read 5'→3' on the strand the PAM was found on.</param>
/// <param name="TargetStart">
/// Forward-strand, 0-based start coordinate of the protospacer (the leftmost of the
/// <see cref="CrisprSystem.GuideLength"/> bases covered by <see cref="TargetSequence"/>), for
/// both strands. For reverse-strand hits <see cref="TargetSequence"/> is therefore the reverse
/// complement of the forward bases <c>[TargetStart, TargetStart + GuideLength)</c>.
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
/// Criterion by which <see cref="CrisprDesigner.DesignGuideRnas"/> ranks its output.
/// </summary>
public enum GuideRnaRanking
{
    /// <summary>Descending <see cref="GuideRnaCandidate.Score"/> (the composition-based quality score). Default.</summary>
    Score = 0,
    /// <summary>
    /// Descending <see cref="GuideRnaCandidate.OnTargetScore"/> (Doench 2016 Rule Set 2 / Azimuth);
    /// candidates with no 30-nt context rank last.
    /// </summary>
    OnTargetRuleSet2 = 1
}

/// <summary>
/// One of the two guide-inefficiency motifs described by Graf et al. (2019), Cell Reports 26:1098–1103.
/// </summary>
public enum GrafMotifType
{
    /// <summary>Neither motif present.</summary>
    None = 0,
    /// <summary>The "TT-motif" at the 3' end (poor Pol III-driven sgRNA expression).</summary>
    TtMotif = 1,
    /// <summary>The "GCC-motif" at the 3' end (<c>[AGT]GCC</c> / <c>GCCT</c>; generally inefficient).</summary>
    GccMotif = 2
}

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
    /// Ranking criterion applied by <see cref="CrisprDesigner.DesignGuideRnas"/>; defaults to
    /// <see cref="GuideRnaRanking.Score"/>.
    /// </summary>
    public GuideRnaRanking Ranking { get; init; }

    /// <summary>
    /// Maximum number of self-complementary 4-bp stems
    /// (<see cref="CrisprDesigner.CountSelfComplementaryStems"/>) a candidate may carry; guides above
    /// it are dropped by <see cref="CrisprDesigner.DesignGuideRnas"/> and flagged as an issue by
    /// <see cref="CrisprDesigner.EvaluateGuideRna"/>. <c>null</c> = no filter, which is the reference
    /// tool's own default (CHOPCHOP <c>filterSelfCompMax = -1</c>).
    /// </summary>
    public int? MaxSelfComplementaryStems { get; init; }

    /// <summary>
    /// Optional sgRNA scaffold regions, on the same strand as the guide, that the
    /// self-complementarity stem count may pair against (CHOPCHOP's <c>-BB/--backbone</c>; e.g.
    /// <see cref="CrisprDesigner.StandardSgRnaBackboneRegion"/>). <c>null</c> = guide-internal stems
    /// only, the reference tool's default.
    /// </summary>
    public IReadOnlyList<string>? SelfComplementarityBackboneRegions { get; init; }

    /// <summary>
    /// Default parameters for guide RNA design. The 40–70% GC acceptance window is the same one the
    /// reference tool CHOPCHOP uses (<c>chopchop.py</c>: <c>GC_LOW = 40</c>, <c>GC_HIGH = 70</c>;
    /// Labun et al., Nucleic Acids Res 44:W272 (2016) / 47:W171 (2019)).
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

    /// <summary>
    /// The 30-nt Rule Set 2 / Azimuth context of this guide read on the protospacer strand
    /// (4 nt 5' flank + 20-nt protospacer + 3-nt NGG PAM + 3 nt 3' flank), or <c>null</c> when it is
    /// unavailable: the guide was scored standalone, the system is not a 20-nt NGG system, the flanks
    /// do not fit inside the input sequence, or the window contains a non-A/C/G/T base.
    /// </summary>
    public string? Context30Mer { get; init; }

    /// <summary>
    /// The Doench et al. (2016) "Rule Set 2" / Azimuth on-target efficacy score of this guide
    /// (conventionally in [0, 1]; higher = predicted more active), or <c>null</c> when
    /// <see cref="Context30Mer"/> is unavailable. This is the published, trained on-target model
    /// (the "efficiency score" CRISPOR reports per guide), independent of the composition-based
    /// <see cref="Score"/>. See <see cref="CrisprDesigner.CalculateOnTargetRuleSet2(string)"/>.
    /// </summary>
    public double? OnTargetScore { get; init; }

    /// <summary>
    /// The Graf et al. (2019) inefficiency motif carried by this guide, if any
    /// (only evaluated for NGG systems, as in CRISPOR). See <see cref="CrisprDesigner.GetGrafMotif"/>.
    /// </summary>
    public GrafMotifType GrafMotif { get; init; }

    /// <summary>
    /// The number of self-complementary 4-bp stems of this guide, per CHOPCHOP's published
    /// self-complementarity measure. See <see cref="CrisprDesigner.CountSelfComplementaryStems"/>;
    /// <see cref="GuideRnaParameters.MaxSelfComplementaryStems"/> turns it into a filter. This is the
    /// sourced measure; <see cref="SelfComplementarityScore"/> is the library's older normalised
    /// complementary-pair fraction, kept for backward compatibility.
    /// </summary>
    public int SelfComplementaryStems { get; init; }
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
