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
        CrisprSystemType.SpCas9_NGA => new CrisprSystem("SpCas9-NGA", "NGA", 20, true, "SpCas9 with NGA PAM (alternative off-target PAM of NGG, Kleinstiver 2015)"),
        CrisprSystemType.SaCas9_NNGRRN => new CrisprSystem("SaCas9-NNGRRN", "NNGRRN", 21, true, "SaCas9 with relaxed NNGRRN PAM (alternative off-target PAM of NNGRRT)"),
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

        foreach (var site in FindPamSitesCore(SequenceExtensions.ToUpperAscii(sequence)!, GetSystem(systemType)))
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
    /// <param name="guideSequence">The protospacer / spacer (PAM excluded), A/C/G/T.</param>
    /// <param name="systemType">The CRISPR system, which fixes the seed side and the PAM-dependent checks.</param>
    /// <param name="parameters">Thresholds; <see cref="GuideRnaParameters.Default"/> when <c>null</c>.</param>
    /// <returns>
    /// The candidate with its metrics and a composite <see cref="GuideRnaCandidate.Score"/>.
    /// </returns>
    /// <remarks>
    /// <para>
    /// <b>The composite <see cref="GuideRnaCandidate.Score"/> is CHOPCHOP's own penalty sum, shown as a
    /// quality number.</b> CHOPCHOP (<c>chopchop.py</c>) ranks guides by an ascending <i>penalty</i>
    /// (<c>sort_CRISPR_guides</c> = <c>sorted(guides, key=attrgetter('score'))</c>) to which two composition terms
    /// contribute: <c>self.score += self.folding * SCORE['FOLDING']</c> (<c>SCORE['FOLDING'] = 1</c>, the 4-bp stem
    /// count of <see cref="CountSelfComplementaryStems"/>) and, in <c>calcGCContent</c>,
    /// <c>if GCcontent &gt; GC_HIGH or GCcontent &lt; GC_LOW: self.score += SCORE['CRISPR_BAD_GC']</c>
    /// (<c>GC_LOW = 40</c>, <c>GC_HIGH = 70</c>, <c>SCORE['CRISPR_BAD_GC'] = 300</c>, a flat amount, not proportional
    /// to the excess). This method returns <c>clamp(100 − penalty, 0, 100)</c> of that sum, i.e. <c>100 − stems</c>
    /// plus a drop to 0 for a guide outside the GC window (the 300 exceeds the whole 0–100 range); the
    /// clamp and the "100 minus" are this library's presentation of the penalty, the two weights are CHOPCHOP's.
    /// </para>
    /// <para>
    /// <b>Two deviations from CHOPCHOP's defaults, both switchable.</b> CHOPCHOP leaves the GC term off in every mode
    /// (<c>CRISPR_DEFAULT["SCORE_GC"] = False</c>, with the comment "this is already scored in many models!"; the only
    /// CLI flag, <c>-c</c>, is <c>store_false</c>, so the command line cannot switch it on); here
    /// <see cref="GuideRnaParameters.ScoreGcContent"/> is <c>true</c> by default because this composite has no
    /// efficiency model to cover GC, and <c>false</c> reproduces CHOPCHOP's default. The stem term is on, as in CHOPCHOP,
    /// and <see cref="GuideRnaParameters.CheckSelfComplementarity"/> = <c>false</c> is its <c>-SC/--noScoreSelfComp</c>.
    /// Weights, window edges (inclusive: exactly 40 % and 70 % are not penalised) and the GC formula
    /// (G+C over the protospacer length, PAM excluded) were checked against <c>chopchop.py</c> and are locked by
    /// oracle tests.
    /// </para>
    /// <para>
    /// <b>Report-only observations (no effect on the score).</b> <c>TTTT</c> (CRISPOR <c>crispor.py</c> prints "It cannot be
    /// transcribed with a U6 or U3 promoter" as a warning, never a score change; CHOPCHOP has no such term), a seed GC
    /// outside 30–80 %, a common cloning site, the Graf 2019 motifs, and a GC outside the window are listed in
    /// <see cref="GuideRnaCandidate.Issues"/> only. The seed-GC bounds and the cloning-enzyme list are this
    /// library's own convenience flags with no published source. Use <see cref="GuideRnaCandidate.OnTargetScore"/>
    /// (Rule Set 2) for a published efficacy estimate.
    /// </para>
    /// <para>
    /// <see cref="GuideRnaCandidate.SelfComplementarityScore"/> is the stem count divided by the number of 4-bp
    /// windows CHOPCHOP examines; it is derived from, and consistent with, <see cref="GuideRnaCandidate.SelfComplementaryStems"/>.
    /// </para>
    /// </remarks>
    public static GuideRnaCandidate EvaluateGuideRna(
        string guideSequence,
        CrisprSystemType systemType = CrisprSystemType.SpCas9,
        GuideRnaParameters? parameters = null)
    {
        if (string.IsNullOrEmpty(guideSequence))
            throw new ArgumentNullException(nameof(guideSequence));

        return EvaluateGuideRnaCore(
            NormalizeNucleotides(guideSequence, nameof(guideSequence)),
            GetSystem(systemType),
            parameters ?? GuideRnaParameters.Default);
    }

    /// <summary>
    /// Folds a raw-string nucleotide argument to upper case (ASCII only, so U+017F 'ſ' is not read as 'S'),
    /// reads <c>U</c> as <c>T</c> (an RNA spacer is its DNA-equivalent), and rejects every other character,
    /// IUPAC ambiguity codes included, with the same <see cref="ArgumentException"/> the file's other raw-string
    /// guards raise (the CFD / SaCas9 / MIT scorers, and the <see cref="DnaSequence"/> alphabet they mirror).
    /// </summary>
    private static string NormalizeNucleotides(string value, string paramName)
    {
        var chars = new char[value.Length];
        for (int i = 0; i < value.Length; i++)
        {
            char c = SequenceExtensions.ToUpperAscii(value[i]);
            if (c == 'U')
                c = 'T';
            else if (c is not ('A' or 'C' or 'G' or 'T'))
                throw new ArgumentException(
                    $"Expected A/C/G/T (U is read as T) but found '{value[i]}' at index {i}.", paramName);
            chars[i] = c;
        }
        return new string(chars);
    }

    private static GuideRnaCandidate EvaluateGuideRnaCore(
        string guideSequence,
        CrisprSystem system,
        GuideRnaParameters effectiveParams)
    {
        var seq = SequenceExtensions.ToUpperAscii(guideSequence)!;

        // Calculate GC content
        double gcContent = CalculateGcContent(seq);

        // Check for polyT (transcription terminator)
        bool hasPolyT = HasPolyT(seq, 4);

        // CHOPCHOP's published self-complementarity measure (4-bp stems); reported always, filtered
        // only when the caller sets a limit (the reference's own default is "report, do not filter").
        int selfCompStems = CountSelfComplementaryStems(seq, 4, effectiveParams.SelfComplementarityBackboneRegions);

        // Check seed region (last 10 bp for Cas9)
        // Evidence: Addgene CRISPR Guide, "Basic Mechanisms of CRISPR"
        // (https://www.addgene.org/guides/crispr/, fetched 2026-10-14): "the seed sequence (8–10 bases
        // at the 3′ end of the gRNA targeting sequence)"; using the upper bound (10 bp). The quote is
        // about Cas9; the PAM-before-target (Cas12a) branch mirrors it at the PAM-proximal 5' end.
        string seedRegion = system.PamAfterTarget
            ? seq.Substring(Math.Max(0, seq.Length - 10))
            : seq.Substring(0, Math.Min(10, seq.Length));
        double seedGc = CalculateGcContent(seedRegion);

        // Score = clamp(100 - penalty, 0, 100), where `penalty` is CHOPCHOP's ascending penalty sum
        // (chopchop.py: Guide.calcSelfComplementarity `self.score += self.folding * SCORE['FOLDING']`,
        // Guide.calcGCContent `self.score += SCORE['CRISPR_BAD_GC']` outside GC_LOW..GC_HIGH).
        // Everything else is reported in `issues` only (CRISPOR treats TTTT / GC extremes as warnings).
        var issues = new List<string>();
        double penalty = 0;

        // GC window (CHOPCHOP: GC_LOW = 40, GC_HIGH = 70, strict comparisons, flat penalty)
        bool gcOutside = gcContent < effectiveParams.MinGcContent || gcContent > effectiveParams.MaxGcContent;
        if (gcOutside)
        {
            issues.Add(gcContent < effectiveParams.MinGcContent
                ? $"Low GC content ({gcContent:F1}%)"
                : $"High GC content ({gcContent:F1}%)");
            if (effectiveParams.ScoreGcContent)
                penalty += BadGcContentPenalty;
        }

        // TTTT: report only (CRISPOR crispor.py: warning, no score change; CHOPCHOP has no such term).
        if (hasPolyT && effectiveParams.AvoidPolyT)
            issues.Add("Contains TTTT (potential Pol III terminator)");

        // Self-complementarity penalty (opt-out via GuideRnaParameters.CheckSelfComplementarity = -SC).
        bool stemPenalty = selfCompStems > 0 && effectiveParams.CheckSelfComplementarity;
        if (stemPenalty)
            penalty += selfCompStems * SelfComplementarityStemPenalty;

        // Seed region GC: report only (the 30-80 % bounds are not from any tool).
        if (seedGc < 30 || seedGc > 80)
            issues.Add($"Suboptimal seed region GC ({seedGc:F1}%)");

        // Common cloning site: report only (not from any tool).
        bool hasRestrictionSite = HasCommonRestrictionSite(seq);
        if (hasRestrictionSite)
            issues.Add("Contains common restriction site");

        // Graf et al. 2019 inefficiency motifs. CRISPOR reports these as a warning for NGG systems
        // ("Inefficient") without altering any score, so no deduction is applied here either.
        var grafMotif = system.PamSequence == "NGG" ? GetGrafMotif(seq) : GrafMotifType.None;
        if (grafMotif == GrafMotifType.TtMotif)
            issues.Add("Graf 2019 TT-motif (inefficient in Pol III-based expression)");
        else if (grafMotif == GrafMotifType.GccMotif)
            issues.Add("Graf 2019 GCC-motif (generally inefficient)");

        // Reported after the other issues so the order of the long-standing ones is unchanged.
        if (stemPenalty)
            issues.Add($"{selfCompStems} self-complementary 4-bp stems (-{selfCompStems * SelfComplementarityStemPenalty:0.##})");

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
            SelfComplementarityScore: SelfComplementarityWindowFraction(seq.Length, selfCompStems),
            Score: Math.Clamp(100 - penalty, 0, 100),
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

        var seq = NormalizeNucleotides(guideSequence, nameof(guideSequence));
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
    /// Predicts potential off-target sites for a guide RNA: enumerates every PAM site of
    /// <paramref name="systemType"/> on both strands of <paramref name="genome"/> and returns those
    /// whose guide-length protospacer differs from the guide by 1..<paramref name="maxMismatches"/>
    /// bases. Exact (0-mismatch) matches are treated as on-targets and are never returned.
    /// </summary>
    /// <param name="guideSequence">The guide RNA sequence; its length must match the system's guide length.</param>
    /// <param name="genome">The genome/sequence to search for off-targets.</param>
    /// <param name="maxMismatches">Maximum number of mismatches allowed, 0-5 (0 yields nothing).</param>
    /// <param name="systemType">Type of CRISPR system.</param>
    /// <returns>
    /// Potential off-target sites in a deterministic order: all forward-strand PAM sites by ascending
    /// position, then all reverse-strand sites. Substitutions only — bulges/indels are out of scope.
    /// </returns>
    /// <remarks>
    /// <see cref="OffTargetSite.OffTargetScore"/> is this library's own seed-weighted penalty (see
    /// <see cref="CalculateHeuristicSpecificityScore"/>), NOT a published per-hit model. For published per-hit
    /// scores use <see cref="CalculateMitHitScore"/> (MIT/Hsu 2013) or
    /// <see cref="CalculateCfdScore"/> (Doench 2016 CFD).
    /// <para>
    /// Recorded deviation on the on-target skip: the reference drops only the 0-mismatch hits that
    /// fall inside the caller-supplied on-target locus (<c>annotateOfftargets</c>, crispor.py
    /// L1920-1929: <c>chrom==inChrom and start &gt;= inStart and end &lt;= inEnd</c>; in its
    /// "auto-ontarget" mode, where no locus is given, just the first 0-mismatch hit it meets). Any
    /// further perfect match elsewhere in the genome IS scored there, as an off-target worth
    /// <c>calcHitScore = 100</c>. This method (whose contract is the 1..max-mismatch enumeration) takes no
    /// on-target locus and drops every exact match genome-wide; a guide with a second perfect match
    /// therefore looks more specific here than in the reference. The CRISPOR-faithful aggregates
    /// (<see cref="CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType, int?, int?)"/>,
    /// <see cref="CalculateMitSpecificityScoreWithAlternativePams"/>,
    /// <see cref="CalculateCfdSpecificityScore(string, DnaSequence, int, CrisprSystemType, CfdAggregation, int?, int?)"/>,
    /// <see cref="CalculateCfdSpecificityScoreWithAlternativePams"/>) accept an optional on-target interval
    /// that reproduces the reference's skip.
    /// </para>
    /// </remarks>
    public static IEnumerable<OffTargetSite> FindOffTargets(
        string guideSequence,
        DnaSequence genome,
        int maxMismatches = 3,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
    {
        return FindOffTargetsCore(
            ValidateScanArguments(guideSequence, genome, maxMismatches, systemType), genome, maxMismatches, systemType);
    }

    /// <summary>
    /// The guards every off-target entry point applies, in this order and with these exception types
    /// (null/empty guide, null genome, <c>maxMismatches</c> outside 0-5, guide length, then the guide
    /// alphabet); returns the normalised guide.
    /// </summary>
    private static string ValidateScanArguments(
        string guideSequence, DnaSequence genome, int maxMismatches, CrisprSystemType systemType)
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
        return NormalizeNucleotides(guideSequence, nameof(guideSequence));
    }

    /// <summary>
    /// Validates the optional on-target interval of the CRISPOR-faithful aggregates: both bounds or neither;
    /// 0 &lt;= start &lt; end. Returns <c>null</c> when neither is given (every exact match is dropped).
    /// </summary>
    private static (int Start, int End)? ResolveOnTargetInterval(int? onTargetStart, int? onTargetEnd)
    {
        if (onTargetStart is null && onTargetEnd is null)
            return null;
        if (onTargetStart is null || onTargetEnd is null)
            throw new ArgumentException(
                "The on-target interval needs both bounds (onTargetStart and onTargetEnd) or neither.",
                onTargetStart is null ? nameof(onTargetStart) : nameof(onTargetEnd));
        if (onTargetStart < 0)
            throw new ArgumentOutOfRangeException(nameof(onTargetStart), "The on-target interval start must be >= 0.");
        if (onTargetEnd <= onTargetStart)
            throw new ArgumentOutOfRangeException(nameof(onTargetEnd),
                "The on-target interval end must be greater than its start (0-based, half-open).");
        return (onTargetStart.Value, onTargetEnd.Value);
    }

    /// <summary>
    /// <see cref="FindOffTargetsWithPam"/> for the CRISPOR-faithful aggregates: the same guards and scan,
    /// plus the optional on-target interval (see <see cref="ScanOffTargetsCore"/>).
    /// </summary>
    private static List<PamAnnotatedOffTarget> CollectSystemHits(
        string guideSequence, DnaSequence genome, int maxMismatches, CrisprSystemType systemType,
        (int Start, int End)? onTarget) =>
        ScanOffTargetsCore(
                ValidateScanArguments(guideSequence, genome, maxMismatches, systemType),
                genome, maxMismatches, systemType, onTarget)
            .Select(static hit => new PamAnnotatedOffTarget(hit.Site, hit.Pam))
            .ToList();

    private static IEnumerable<OffTargetSite> FindOffTargetsCore(
        string guideSequence, DnaSequence genome, int maxMismatches, CrisprSystemType systemType) =>
        FindOffTargetsWithPamCore(guideSequence, genome, maxMismatches, systemType)
            .Select(static hit => hit.Site);

    /// <summary>
    /// The off-target scan of <see cref="FindOffTargetsCore"/>, additionally yielding each site's own
    /// PAM as read on the protospacer strand (<see cref="PamSite.PamSequence"/>). Scores that need the
    /// off-target PAM (CFD) use this so that neither the PAM window nor its strand handling is
    /// re-derived from coordinates.
    /// </summary>
    private static IEnumerable<(OffTargetSite Site, string Pam)> FindOffTargetsWithPamCore(
        string guideSequence, DnaSequence genome, int maxMismatches, CrisprSystemType systemType) =>
        ScanOffTargetsCore(guideSequence, genome, maxMismatches, systemType, onTarget: null)
            .Select(static hit => (hit.Site, hit.Pam));

    /// <summary>
    /// The scan behind every off-target entry point. With <paramref name="onTarget"/> <c>null</c> a site is
    /// reported when it has <c>1..maxMismatches</c> mismatches (every exact match is dropped genome-wide:
    /// the documented default). With an on-target interval, exact (0-mismatch) sites are reported too,
    /// except those whose footprint lies inside the interval - CRISPOR's skip, <c>annotateOfftargets</c>
    /// crispor.py L1920-1929: <c>editDist==0 and chrom==inChrom and start &gt;= inStart and end &lt;= inEnd</c>.
    /// The footprint is the protospacer together with its PAM, as the reference's <c>start</c>/<c>end</c>
    /// (a 23-mer for SpCas9; <c>parsePos</c> defaults <c>end</c> to <c>start+23</c>), on the forward strand,
    /// 0-based half-open, for either strand's hit (the reference never compares the strand).
    /// </summary>
    private static IEnumerable<(OffTargetSite Site, string Pam)> ScanOffTargetsCore(
        string guideSequence, DnaSequence genome, int maxMismatches, CrisprSystemType systemType,
        (int Start, int End)? onTarget)
    {
        var system = GetSystem(systemType);

        var guide = SequenceExtensions.ToUpperAscii(guideSequence)!;

        // Find all PAM sites
        var pamSites = FindPamSitesCore(genome.Sequence, system);

        foreach (var pamSite in pamSites)
        {
            int mismatches = CountMismatches(guide, pamSite.TargetSequence);

            if (mismatches == 0)
            {
                if (onTarget is not { } locus)
                    continue;
                int footprintStart = Math.Min(pamSite.TargetStart, pamSite.Position);
                int footprintEnd = Math.Max(pamSite.TargetStart + system.GuideLength,
                    pamSite.Position + pamSite.PamSequence.Length);
                if (footprintStart >= locus.Start && footprintEnd <= locus.End)
                    continue; // the on-target: not an off-target
            }
            else if (mismatches > maxMismatches)
            {
                continue;
            }

            // Calculate off-target score based on mismatch positions
            double score = CalculateOffTargetScore(guide, pamSite.TargetSequence, system);

            yield return (new OffTargetSite(
                Position: pamSite.Position,
                Sequence: pamSite.TargetSequence,
                Mismatches: mismatches,
                MismatchPositions: GetMismatchPositions(guide, pamSite.TargetSequence),
                IsForwardStrand: pamSite.IsForwardStrand,
                OffTargetScore: score), pamSite.PamSequence);
        }
    }

    /// <summary>
    /// Finds potential off-target sites exactly as <see cref="FindOffTargets"/> does, additionally
    /// reporting each site's own PAM as read on the protospacer strand. Use this when a score needs the
    /// off-target's PAM: the CFD model (<see cref="CalculateCfdScore"/>) and the alternative-PAM-aware
    /// MIT aggregate
    /// (<see cref="CalculateMitSpecificityScore(string, IEnumerable{PamAnnotatedOffTarget}, CrisprSystemType)"/>),
    /// whose downweight is defined on the PAM the hit actually carries.
    /// </summary>
    /// <param name="guideSequence">The guide sequence (must match the system's guide length).</param>
    /// <param name="genome">The genome/sequence to scan.</param>
    /// <param name="maxMismatches">Maximum mismatches to consider (0-5).</param>
    /// <param name="systemType">The CRISPR system whose PAM is scanned.</param>
    /// <returns>
    /// The same sites <see cref="FindOffTargets"/> yields, in the same deterministic order, each paired
    /// with its PAM. The PAM always matches <paramref name="systemType"/>'s motif — one call scans one
    /// system — so an alternative-PAM aggregate combines the results of several calls.
    /// </returns>
    /// <exception cref="ArgumentNullException">Thrown when the guide is null/empty or the genome is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="maxMismatches"/> is outside 0-5.</exception>
    /// <exception cref="ArgumentException">Thrown when the guide length does not match the system's.</exception>
    public static IEnumerable<PamAnnotatedOffTarget> FindOffTargetsWithPam(
        string guideSequence,
        DnaSequence genome,
        int maxMismatches = 3,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
    {
        // Same guards, same order, same exception types as FindOffTargets.
        return FindOffTargetsWithPamCore(
                ValidateScanArguments(guideSequence, genome, maxMismatches, systemType), genome, maxMismatches, systemType)
            .Select(static hit => new PamAnnotatedOffTarget(hit.Site, hit.Pam));
    }

    /// <summary>
    /// Calculates THIS LIBRARY'S OWN HEURISTIC specificity score for a guide RNA
    /// (higher = more specific):
    /// <c>max(0, 100 − Σ over off-targets of (5 per seed-region mismatch + 2 per distal mismatch))</c>,
    /// over the off-targets found by <see cref="FindOffTargets"/> with a fixed cap of 4 mismatches.
    /// The method name carries the word "heuristic" deliberately: no publication defines this score
    /// (see the remarks), while the two sibling aggregates in this class —
    /// <see cref="CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType, int?, int?)"/> and
    /// <see cref="CalculateCfdSpecificityScore(string, DnaSequence, int, CrisprSystemType, CfdAggregation, int?, int?)"/>
    /// — are reference-verified published models.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>THIS IS THIS LIBRARY'S OWN HEURISTIC, NOT A PUBLISHED MODEL.</b> Hsu et al. 2013
    /// (Nat Biotechnol 31:827) motivates only the *shape* of the weighting — that mismatches in the
    /// PAM-proximal seed region are less tolerated than PAM-distal ones. The numbers 5, 2, the 12-nt
    /// seed width, the linear sum and the subtraction from 100 are not from Hsu 2013 or from any other
    /// publication; no source was found for them, and none is claimed. Consequences: the score is
    /// linear in the mismatch count (the published models are multiplicative), it ignores which bases
    /// are involved and the off-target's own PAM, and it saturates at 0 after ~20 seed mismatches'
    /// worth of off-target burden, so two guides with very different real off-target profiles can
    /// share a score of 0. Treat it as a coarse triage number only.
    /// </para>
    /// <para>
    /// For a published, reference-verified model use one of the aggregates in this same class:
    /// <see cref="CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType, int?, int?)"/> (MIT/Hsu
    /// 2013 guide score, <c>100 / (100 + Σ hitScores) × 100</c> over per-hit
    /// <see cref="CalculateMitHitScore"/> values) or
    /// <see cref="CalculateCfdSpecificityScore(string, DnaSequence, int, CrisprSystemType, CfdAggregation, int?, int?)"/>
    /// (Doench 2016 CFD, aggregated as CRISPOR's <c>annotateOfftargets</c> does). Both reproduce the
    /// reference tool numerically; this method does not.
    /// </para>
    /// </remarks>
    /// <param name="guideSequence">The guide sequence; its length must match the system's guide length.</param>
    /// <param name="genome">The genome/sequence to scan for off-targets.</param>
    /// <param name="systemType">CRISPR system (SpCas9 by default).</param>
    /// <returns>The heuristic specificity in [0, 100]; 100 when no off-targets are found.</returns>
    public static double CalculateHeuristicSpecificityScore(
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

    /// <summary>
    /// Backward-compatible alias for
    /// <see cref="CalculateHeuristicSpecificityScore(string, DnaSequence, CrisprSystemType)"/>, which
    /// is the name to prefer: it says in the signature that this score is this library's own
    /// heuristic and not a published model. Kept so that existing callers (and the
    /// <c>crispr_specificity_score</c> MCP tool name, which is part of the published catalog) keep
    /// working unchanged; it forwards with no behaviour of its own and returns the identical value.
    /// Not marked <c>[Obsolete]</c> on purpose — the heuristic is a legitimate triage number when
    /// its contract is understood, and the repository builds with warnings as errors.
    /// </summary>
    /// <param name="guideSequence">See the primary overload.</param>
    /// <param name="genome">See the primary overload.</param>
    /// <param name="systemType">See the primary overload.</param>
    /// <returns>Exactly <see cref="CalculateHeuristicSpecificityScore(string, DnaSequence, CrisprSystemType)"/>.</returns>
    public static double CalculateSpecificityScore(
        string guideSequence,
        DnaSequence genome,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
        => CalculateHeuristicSpecificityScore(guideSequence, genome, systemType);

    /// <summary>
    /// Counts guide-to-target mismatches. Equal-length inputs are delegated to the canonical
    /// case-insensitive <see cref="SequenceExtensions.HammingDistance(ReadOnlySpan{char}, ReadOnlySpan{char})"/>;
    /// unequal lengths keep this type's historical <see cref="int.MaxValue"/> "not comparable"
    /// sentinel (the canonical throws), so callers that pre-validate the length are unaffected.
    /// </summary>
    private static int CountMismatches(string seq1, string seq2) =>
        seq1.Length != seq2.Length
            ? int.MaxValue
            : seq1.AsSpan().HammingDistance(seq2.AsSpan());

    /// <summary>
    /// Lists the 0-based positions at which <paramref name="guide"/> and <paramref name="target"/>
    /// differ, using the same case-insensitive ASCII comparison as the canonical
    /// <see cref="SequenceExtensions.HammingDistance(ReadOnlySpan{char}, ReadOnlySpan{char})"/>
    /// so that <c>MismatchPositions.Count == Mismatches</c> holds for every equal-length input
    /// (TestSpec INV-03). Unequal lengths are not comparable and yield an empty list, matching the
    /// <see cref="int.MaxValue"/> sentinel of <see cref="CountMismatches"/>; the off-target scan never
    /// produces them, as guide and protospacer are both of the system's guide length.
    /// </summary>
    private static IReadOnlyList<int> GetMismatchPositions(string guide, string target)
    {
        if (guide.Length != target.Length)
            return Array.Empty<int>();

        var positions = new List<int>();
        for (int i = 0; i < guide.Length; i++)
        {
            if (SequenceExtensions.ToUpperAscii(guide[i]) != SequenceExtensions.ToUpperAscii(target[i]))
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

        // Uses the same mismatch set as the reported Mismatches / MismatchPositions (one case-folding
        // rule for all three), so the penalty can never disagree with the reported mismatch count.
        foreach (int i in GetMismatchPositions(guide, target))
        {
            // Seed-region (PAM-proximal) mismatches are LESS tolerated by Cas9
            // (Hsu et al. 2013), so they receive a higher penalty here.
            bool inSeed = i >= seedStart && i < seedEnd;
            score += inSeed ? 5 : 2;
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
    //   "TATAGCTGCGATCTGAGGTAGGGAGGGACC" -> 0.7130893469547245 (the literal the reference *prints* as
    //       its own "expected result", 0.713089368437, is stale by ~2e-8 relative to the table it
    //       ships; the table is authoritative and is what both sides compute)
    //   "TCCGCACCTGTCACGGTCGGGGCTTGGCGC" -> 0.018983843193336125 (reference printed 0.0189838463593)
    //   Bulk: 477 contexts (400 pseudo-random NGG 30-mers, seed 20261014, + 77 `Context30Mer`
    //       values from DesignGuideRnas) scored by the reference loaded straight out of the
    //       downloaded `doenchScore.py` text -> 0 mismatches, max |delta| exactly 0 (bit-equal).
    //       Both GC branches and the gcCount == 10 boundary are covered (265 <= 10, 212 > 10).
    //
    // The reference's own CRISPOR pipeline reports this score as `int(100 * sigmoid)`
    //   (`crisporEffScores.py`, calcDoenchScores), i.e. truncated to a whole percent. This method
    //   returns the UNROUNDED sigmoid x 100, so it is the same model at full precision; a CRISPOR
    //   column value is `Math.Floor` of what this returns.
    //
    // NGG guard: the bare reference function does not check the PAM (CRISPOR only ever calls it for
    //   SpCas9/NGG guides), but the published feature table itself presumes NGG - every dinucleotide
    //   feature at offset 24 ends in G ((24,"AG"/"CG"/"TG")) and the offset-26 feature starts with G
    //   ((26,"GT")) - so scoring a non-NGG context would silently leave the model's domain. The
    //   guard is this library's explicit contract, matching AzimuthRuleSet2.Validate.
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
    /// <exception cref="ArgumentNullException"><paramref name="context30Mer"/> is <c>null</c>.</exception>
    /// <exception cref="ArgumentException">
    /// <paramref name="context30Mer"/> is not 30 characters long (the empty string included), contains
    /// a character other than A/C/G/T, or does not carry an NGG PAM at offsets 25-26.
    /// </exception>
    /// <remarks>
    /// This is the published, exactly-reproducible linear model (logistic regression):
    /// intercept + a GC-count term over the 20-nt protospacer (two published coefficients, one for
    /// GC &lt;= 10 and one for GC &gt; 10) + 70 position-specific single- and di-nucleotide weights,
    /// through a logistic sigmoid.
    /// It is NOT Doench "Rule Set 2" / Azimuth, which is a gradient-boosted-tree model that cannot
    /// be reproduced from published coefficients without the trained model file.
    /// Source: Doench et al. 2014, Nat Biotechnol 32:1262 (PMID 25184501); coefficients per the
    /// reference implementation in CRISPOR's <c>doenchScore.py</c>
    /// (https://github.com/maximilianh/crisporWebsite/blob/master/doenchScore.py, repository root),
    /// re-verified byte-for-byte and numerically against that file (see the region comment above).
    /// The reference's own CRISPOR pipeline reports the same score truncated to a whole percent
    /// (<c>int(100 * sigmoid)</c>, <c>crisporEffScores.py</c>); this method returns it unrounded.
    /// </remarks>
    public static double CalculateOnTargetDoench2014(string context30Mer)
    {
        // A null reference is a missing argument; an empty string is a wrong *value*, and the length
        // check below reports it the same way as every other malformed context. The reference draws
        // the same line (`assert isinstance(seq[0], basestring)` vs `assert len(seq[0]) > 0`,
        // azimuth/model_comparison.py L523-L524). Identical fix and rationale as WP1's
        // AzimuthRuleSet2.Validate("") (B08 F36); the test that locked ArgumentNullException for ""
        // was updated with the reason. The MCP tool is unaffected: it rejects null/empty itself,
        // as ArgumentException, before the library is reached.
        ArgumentNullException.ThrowIfNull(context30Mer);

        var seq = SequenceExtensions.ToUpperAscii(context30Mer)!;

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
    /// <param name="aminoAcidCutPosition">
    /// The amino-acid position of the cut site within the target protein. Passed to the model as-is
    /// (the reference does not treat any value of it specially).
    /// </param>
    /// <param name="percentPeptide">
    /// The cut position as a percentage [0, 100] along the coding sequence (how far into the protein),
    /// or <c>-1</c> for Azimuth's "no gene position" sentinel.
    /// </param>
    /// <returns>The predicted on-target activity, conventionally in the range [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="percentPeptide"/> is neither in [0, 100] nor the <c>-1</c> sentinel.
    /// </exception>
    /// <remarks>
    /// Uses Azimuth's <c>V3_model_full</c>, except for <paramref name="percentPeptide"/> == <c>-1</c>:
    /// the reference answers a request carrying that sentinel with the sequence-only
    /// <c>V3_model_nopos</c> model (<c>azimuth/model_comparison.py</c> L539 selects
    /// <c>V3_model_nopos.pickle</c> when <c>np.any(percent_peptide == -1)</c>, and L563 fills the
    /// gene-position frame with -1 rather than the caller's values), so this method returns exactly
    /// <see cref="CalculateOnTargetRuleSet2(string)"/> in that case.
    /// When gene context is unavailable, prefer the sequence-only
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

        var g = NormalizeNucleotides(guide20, nameof(guide20));
        var o = NormalizeNucleotides(offTarget20, nameof(offTarget20));

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
    /// <param name="guideSequence">The guide sequence (20 nt for SpCas9; see the length note below).</param>
    /// <param name="genome">The genome/sequence to scan for off-targets.</param>
    /// <param name="maxMismatches">Maximum mismatches to consider (0-5, as in <see cref="FindOffTargets"/>).</param>
    /// <param name="systemType">CRISPR system (SpCas9 by default). Must be a PAM-after-target (Cas9-class) system.</param>
    /// <param name="onTargetStart">
    /// Optional on-target interval, start (0-based, inclusive, forward-strand coordinates of
    /// <paramref name="genome"/>; both bounds or neither). See the remark on the on-target interval.
    /// </param>
    /// <param name="onTargetEnd">Optional on-target interval, end (exclusive); see <paramref name="onTargetStart"/>.</param>
    /// <returns>The aggregate MIT/Hsu specificity in [0, 100]; 100 when no off-targets are found.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="systemType"/> is a Cas12a/Cpf1-class (PAM-before-target) system, for
    /// which the reference provides no off-target model, or when the guide is shorter than 19 nt.
    /// </exception>
    /// <remarks>
    /// The per-hit model is chosen by system, exactly as the reference's <c>annotateOfftargets</c> does:
    /// a SaCas9 PAM (<c>NNGRRT</c>/<c>NNNRRT</c>) with a 21-nt guide is scored with the SaCas9-specific
    /// model of Tycko et al. 2018 (<see cref="CalculateSaCas9HitScore"/>, the reference's
    /// <c>calcSaHitScore</c>); every other PAM-after-target system uses the Hsu
    /// <see cref="CalculateMitHitScore"/>. A non-SaCas9 guide longer than 20 nt falls back to
    /// <c>calcHitScore</c>'s own last-20 fudge, and a 19-nt guide to its 'A'-pad.
    /// <para>
    /// Recorded deviation: this aggregate sums only off-targets of the one scanned system's PAM, whereas
    /// CRISPOR runs a single wider scan and additionally folds in alternative-PAM hits (an NGG scan that
    /// turns up a non-GG site) with a 0.2 weight. That weight cannot fire on a single-system scan — an
    /// SpCas9 scan yields only NGG sites and an SpCas9-NAG scan only NAG sites, so the reference's own
    /// condition (<c>pam == "NGG" and otSeq[-2:] != "GG"</c>) is never true here. For the reference's
    /// behaviour use
    /// <see cref="CalculateMitSpecificityScoreWithAlternativePams"/> (scans the main PAM and the
    /// alternative PAMs this library offers as systems) or
    /// <see cref="CalculateMitSpecificityScore(string, IEnumerable{PamAnnotatedOffTarget}, CrisprSystemType)"/>
    /// (aggregates any PAM-annotated hit set). This overload's behaviour is unchanged.
    /// </para>
    /// <para>
    /// <b>On-target interval (B08 round 3 G7).</b> By default every exact (0-mismatch) match in the genome is
    /// dropped, so a guide with a second perfect site looks more specific than in CRISPOR. The reference
    /// (<c>annotateOfftargets</c>, crispor.py L1920-1929) drops only 0-mismatch hits inside the caller-supplied
    /// on-target interval (<c>chrom==inChrom and start &gt;= inStart and end &lt;= inEnd</c>) and scores any other
    /// perfect match as an off-target worth <c>calcHitScore = 100</c> (CFD: the PAM score). Passing
    /// <paramref name="onTargetStart"/> / <paramref name="onTargetEnd"/> reproduces that: exact matches whose
    /// footprint - protospacer plus PAM, on the forward strand, 0-based half-open, the coordinates of
    /// <see cref="PamSite"/> - lies inside the interval are skipped (either strand, the reference never compares
    /// it), all others are aggregated. The reference's "auto-ontarget" mode (no interval: the first 0-mismatch hit it
    /// meets, in the aligner's output order) is deliberately not offered: that order is unspecified, and under
    /// CFD or with alternative PAMs the hits differ in score, so no deterministic equivalent exists.
    /// </para>
    /// </remarks>
    public static double CalculateMitSpecificityScore(
        string guideSequence,
        DnaSequence genome,
        int maxMismatches = 4,
        CrisprSystemType systemType = CrisprSystemType.SpCas9,
        int? onTargetStart = null,
        int? onTargetEnd = null)
    {
        EnsureMitModelApplies(systemType, nameof(systemType));
        var onTarget = ResolveOnTargetInterval(onTargetStart, onTargetEnd);

        var offTargets = CollectSystemHits(guideSequence, genome, maxMismatches, systemType, onTarget)
            .Select(static hit => hit.Site).ToList();
        if (offTargets.Count == 0)
            return 100.0;

        var guide = NormalizeNucleotides(guideSequence, nameof(guideSequence));

        // The reference picks the per-hit model by system, not one model for everything
        // (annotateOfftargets): a SaCas9 PAM uses the Tycko 2018 model, everything else calcHitScore.
        bool useSaModel = IsSaCas9Pam(GetSystem(systemType).PamSequence) && guide.Length == SaCas9GuideLength;
        var hitScores = offTargets
            .Select(ot => useSaModel
                ? CalculateSaCas9HitScore(guide, SequenceExtensions.ToUpperAscii(ot.Sequence)!)
                : CalculateMitHitScoreFudged(guide, SequenceExtensions.ToUpperAscii(ot.Sequence)!));

        return CalculateMitSpecificityScore(hitScores);
    }

    /// <summary>
    /// The factor the reference applies to an alternative-PAM off-target's MIT/Hsu single-hit score:
    /// <c>mitScore = mitScore * 0.2</c> (crispor.py <c>annotateOfftargets</c>, L1947-1958). Its own
    /// comment states the rationale: "this is a heuristic based on the guideSeq data where alternative
    /// PAMs represent only ~10% of all cleaveage events. We divide the MIT score by 5 to make sure that
    /// these off-targets are not ranked among the top but still appear in the list somewhat".
    /// </summary>
    private const double AlternativePamMitWeight = 0.2;

    /// <summary>
    /// The alternative off-target PAM motifs the reference will accept for a system's main PAM
    /// (crispor.py <c>offtargetPams</c>, L243-253), or an empty list when it lists none for that PAM.
    /// </summary>
    /// <remarks>
    /// Verbatim from the reference table, with its own citations: <c>NGG → NAG, NGA</c> ("MIT and eCrisp
    /// do that, they use the motif NGG + NAG, we add one more, based on the guideSeq results in Tsai et
    /// al, Nat Biot 2014"; "The NGA -&gt; NGG rule was described by Kleinstiver...Young 2015 'Improved
    /// Cas9 Specificity...'"), <c>NGA → NGG</c>, <c>NNGRRT → NNGRRN</c> ("NNGTRRT rule for S. aureus is
    /// in the new protocol 'SaCas9 User manual'"), <c>TTTV → TTTN</c>. The reference also requires
    /// "the length of the alternate PAM has to be the same as the original PAM". Entries of the table
    /// whose main PAM no <see cref="CrisprSystemType"/> uses (<c>NGK</c>, <c>ATTN</c>, <c>TTYN</c>) are
    /// not represented here.
    /// <para>
    /// Only the MIT downweight
    /// (<see cref="CalculateMitSpecificityScore(string, IEnumerable{PamAnnotatedOffTarget}, CrisprSystemType)"/>)
    /// is defined on these motifs; the reference's pipeline-level cutoff for alternative-PAM hits
    /// (<c>ALTPAMMINSCORE = 1.0</c>, L310-318, applied in the bowtie/BWA pre-filter at L4092-4133, not
    /// in the scoring function) is the optional <c>minAlternativePamHitScore</c> of
    /// <see cref="CalculateMitSpecificityScoreWithAlternativePams"/> (default 1.0, as there); the
    /// aggregate overload over hand-supplied hits does not filter.
    /// </para>
    /// </remarks>
    /// <param name="systemType">The CRISPR system whose main PAM is looked up.</param>
    /// <returns>The alternative PAM motifs (IUPAC), possibly empty; never null.</returns>
    public static IReadOnlyList<string> GetAlternativeOffTargetPams(CrisprSystemType systemType) =>
        GetSystem(systemType).PamSequence switch
        {
            "NGG" => ["NAG", "NGA"],
            "NGA" => ["NGG"],
            "NNGRRT" => ["NNGRRN"],
            "TTTV" => ["TTTN"],
            _ => []
        };

    /// <summary>
    /// Calculates the aggregate MIT / Hsu 2013 specificity score from off-targets that each carry their
    /// own PAM, applying the reference's alternative-PAM downweight: an off-target of an NGG (SpCas9)
    /// guide whose own PAM does not end in <c>GG</c> contributes only
    /// <c>0.2 ×</c> its single-hit score. This is the overload to use when the off-targets come from
    /// MORE THAN ONE PAM scan — e.g. an SpCas9 (<c>NGG</c>) scan combined with an SpCas9-NAG
    /// (<c>NAG</c>) scan, which is how this library reproduces CRISPOR's single wider search.
    /// </summary>
    /// <param name="guideSequence">The guide sequence (must match <paramref name="systemType"/>'s guide length).</param>
    /// <param name="offTargets">
    /// The off-targets to aggregate, each with the PAM it actually carries (as
    /// <see cref="FindOffTargetsWithPam"/> reports it). May mix PAMs freely; the main-vs-alternative
    /// decision is taken per hit from its own PAM, never from which scan produced it.
    /// </param>
    /// <param name="systemType">
    /// The guide's OWN system, i.e. the one whose PAM is the "main" PAM — not the system an individual
    /// hit was found with. Must be a PAM-after-target (Cas9-class) system.
    /// </param>
    /// <returns>The aggregate MIT/Hsu specificity in [0, 100]; 100 for an empty off-target set.</returns>
    /// <exception cref="ArgumentNullException">Thrown when the guide is null/empty or <paramref name="offTargets"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="systemType"/> is a Cas12a/Cpf1-class (PAM-before-target) system, when
    /// the guide length does not match the system's, or when an element is null / carries an empty PAM
    /// or a protospacer of a different length than the guide.
    /// </exception>
    /// <remarks>
    /// Source: CRISPOR <c>crispor.py</c> <c>annotateOfftargets</c>, L1947-1958 — after the per-hit MIT
    /// score is computed, <c>if pam=="NGG" and otSeq[-2:]!="GG": mitScore = mitScore * 0.2</c>. Three
    /// properties of that rule are reproduced exactly:
    /// <list type="bullet">
    /// <item><description>It is gated on the GUIDE's main PAM being literally <c>NGG</c>. An
    /// SpCas9-NAG guide (main PAM <c>NAG</c>) gets no downweight on any hit, because the reference's
    /// condition tests the searched PAM pattern, not the hit's.</description></item>
    /// <item><description>The test is on the hit's last two PAM bases, NOT on membership of
    /// <see cref="GetAlternativeOffTargetPams"/>: any non-<c>*GG</c> PAM (<c>NAG</c>, <c>NGA</c>, or
    /// anything else a caller supplies) is downweighted. This method therefore applies it to the
    /// concatenation <c>protospacer + PAM</c>, which is the reference's own <c>otSeq</c>.</description></item>
    /// <item><description>Only the MIT score is downweighted. The reference computes the hit's CFD
    /// score afterwards, from the full <c>otSeq</c> including its PAM, with no 0.2 factor — so
    /// <see cref="CalculateCfdSpecificityScore(IEnumerable{double}, CfdAggregation)"/> must NOT be fed
    /// downweighted values. (The CFD PAM table already penalises non-NGG PAMs on its own.)</description></item>
    /// </list>
    /// The per-hit model is chosen by system exactly as in
    /// <see cref="CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType, int?, int?)"/> (SaCas9 →
    /// Tycko 2018, everything else → Hsu with <c>calcHitScore</c>'s length handling), and the
    /// guide-level aggregate is the unrounded <c>calcMitGuideScore</c>.
    /// <para>
    /// Deliberate scope: the reference also pre-filters alternative-PAM hits with
    /// <c>ALTPAMMINSCORE = 1.0</c> in its search pipeline — that cutoff is available on
    /// <see cref="CalculateMitSpecificityScoreWithAlternativePams"/>, not applied to hand-supplied hits.
    /// </para>
    /// </remarks>
    public static double CalculateMitSpecificityScore(
        string guideSequence,
        IEnumerable<PamAnnotatedOffTarget> offTargets,
        CrisprSystemType systemType = CrisprSystemType.SpCas9)
    {
        if (string.IsNullOrEmpty(guideSequence))
            throw new ArgumentNullException(nameof(guideSequence));
        ArgumentNullException.ThrowIfNull(offTargets);
        EnsureMitModelApplies(systemType, nameof(systemType));

        var system = GetSystem(systemType);
        if (guideSequence.Length != system.GuideLength)
            throw new ArgumentException(
                $"Guide length ({guideSequence.Length}) does not match expected length ({system.GuideLength}) for {system.Name}.",
                nameof(guideSequence));

        var guide = NormalizeNucleotides(guideSequence, nameof(guideSequence));
        bool useSaModel = IsSaCas9Pam(system.PamSequence) && guide.Length == SaCas9GuideLength;
        // The reference gates the downweight on the SEARCHED pattern being exactly "NGG".
        bool mainPamIsNgg = string.Equals(system.PamSequence, "NGG", StringComparison.Ordinal);

        var hitScores = new List<double>();
        foreach (var hit in offTargets)
        {
            if (hit is null)
                throw new ArgumentException("Off-target list contains a null element.", nameof(offTargets));
            if (string.IsNullOrEmpty(hit.Pam))
                throw new ArgumentException(
                    "Every off-target must carry its own PAM; the alternative-PAM rule is defined on it.",
                    nameof(offTargets));

            string offTarget = SequenceExtensions.ToUpperAscii(hit.Site.Sequence)!;
            if (offTarget.Length != guide.Length)
                throw new ArgumentException(
                    $"Off-target protospacer length ({offTarget.Length}) does not match the guide ({guide.Length}).",
                    nameof(offTargets));

            double hitScore = useSaModel
                ? CalculateSaCas9HitScore(guide, offTarget)
                : CalculateMitHitScoreFudged(guide, offTarget);

            // crispor.py L1955: if pam=="NGG" and otSeq[-2:]!="GG": mitScore = mitScore * 0.2.
            // otSeq is the protospacer WITH its PAM, so the test is written on that concatenation here
            // rather than on the PAM alone — identical for the 3-nt PAMs of every NGG system, and the
            // reference's expression verbatim for any other PAM length a caller supplies.
            string otSeq = offTarget + SequenceExtensions.ToUpperAscii(hit.Pam)!;
            if (mainPamIsNgg && !otSeq.EndsWith("GG", StringComparison.Ordinal))
                hitScore *= AlternativePamMitWeight;

            hitScores.Add(hitScore);
        }

        return CalculateMitSpecificityScore(hitScores);
    }

    /// <summary>
    /// Calculates the aggregate MIT / Hsu 2013 specificity score for a guide against a genome over the
    /// guide's own PAM <b>and</b> the reference's alternative PAMs for it
    /// (<see cref="GetAlternativeOffTargetPams"/>), downweighting every alternative-PAM hit by 0.2 as
    /// CRISPOR does. This is the closest analogue of CRISPOR's guide score: its search accepts both the
    /// main and the alternative PAMs in one pass, whereas
    /// <see cref="CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType, int?, int?)"/> sees only
    /// the one system's PAM and therefore reports a guide as more specific than the reference does.
    /// </summary>
    /// <param name="guideSequence">The guide sequence (must match the system's guide length).</param>
    /// <param name="genome">The genome/sequence to scan for off-targets.</param>
    /// <param name="maxMismatches">Maximum mismatches to consider (0-5, as in <see cref="FindOffTargets"/>).</param>
    /// <param name="systemType">The guide's own system (SpCas9 by default); must be PAM-after-target.</param>
    /// <param name="minAlternativePamHitScore">
    /// Minimum unweighted hit score (0-100 scale) an alternative-PAM hit needs to be counted — the
    /// reference's <c>ALTPAMMINSCORE</c> (default <c>1.0</c>, as there; <c>0</c> = no filtering).
    /// </param>
    /// <param name="onTargetStart">
    /// Optional on-target interval, start (0-based, inclusive, forward-strand coordinates of
    /// <paramref name="genome"/>; both bounds or neither). See the remark on the on-target interval.
    /// </param>
    /// <param name="onTargetEnd">Optional on-target interval, end (exclusive); see <paramref name="onTargetStart"/>.</param>
    /// <returns>The aggregate MIT/Hsu specificity in [0, 100]; 100 when no off-target is found at all.</returns>
    /// <exception cref="ArgumentNullException">Thrown when the guide is null/empty or the genome is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="maxMismatches"/> is outside 0-5 or the cutoff is negative/NaN.</exception>
    /// <exception cref="ArgumentException">
    /// Thrown for a PAM-before-target system, or when the guide length does not match the system's.
    /// </exception>
    /// <remarks>
    /// Coverage, stated exactly: an alternative PAM can only be scanned here if some
    /// <see cref="CrisprSystemType"/> has it as its main PAM and the SAME guide length. SpCas9 (<c>NGG</c>)
    /// is covered for <c>NAG</c> (<see cref="CrisprSystemType.SpCas9_NAG"/>) and <c>NGA</c>
    /// (<see cref="CrisprSystemType.SpCas9_NGA"/>); SaCas9 (<c>NNGRRT</c>) for <c>NNGRRN</c>
    /// (<see cref="CrisprSystemType.SaCas9_NNGRRN"/>, 21-nt guide). A site matching both the main and a
    /// wider alternative pattern (every <c>NNGRRT</c> site is also <c>NNGRRN</c>) is one off-target, as
    /// in the reference's single search, and is counted once. The 0.2 weight is the reference's
    /// <c>pam=="NGG"</c>-gated rule, so it applies to NGG guides only; SaCas9 <c>NNGRRN</c> hits enter
    /// unweighted, exactly as <c>annotateOfftargets</c> does. Cas12a-class (<c>TTTV</c> → <c>TTTN</c>)
    /// has no PAM-after-target alternative system and throws here as in the genome overload. The PAMs
    /// actually scanned are <see cref="GetAlternativeOffTargetPams"/> intersected with that system set.
    /// <para>
    /// <paramref name="minAlternativePamHitScore"/> is the reference's <c>ALTPAMMINSCORE</c> /
    /// <c>--minAltPamScore</c> (crispor.py L318, L8439-8441), applied as L4108-4122 does: a hit whose PAM
    /// is in the alternative-PAM list is dropped when its unweighted hit score is below the cutoff
    /// (<c>if otScore &lt; minScore: continue</c>). The default is the reference's own
    /// <c>1.0</c> (B08 round 3 G6; it was <c>0</c> = no filtering before, which is still available by passing
    /// <c>0</c>). Because the reference tests
    /// the alternative list first, with SaCas9 the cutoff also applies to main-PAM (<c>NNGRRT</c>) hits.
    /// </para>
    /// <para>
    /// The optional on-target interval works as in
    /// <see cref="CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType, int?, int?)"/>
    /// and is applied to the main-PAM scan and to every alternative-PAM scan alike (the reference skips a
    /// 0-mismatch hit in the interval whatever its PAM).
    /// </para>
    /// </remarks>
    public static double CalculateMitSpecificityScoreWithAlternativePams(
        string guideSequence,
        DnaSequence genome,
        int maxMismatches = 4,
        CrisprSystemType systemType = CrisprSystemType.SpCas9,
        double minAlternativePamHitScore = 1.0,
        int? onTargetStart = null,
        int? onTargetEnd = null)
    {
        EnsureMitModelApplies(systemType, nameof(systemType));
        var onTarget = ResolveOnTargetInterval(onTargetStart, onTargetEnd);
        var hits = CollectAlternativePamHits(
            guideSequence, genome, maxMismatches, systemType, minAlternativePamHitScore, onTarget);

        return CalculateMitSpecificityScore(guideSequence, hits, systemType);
    }

    /// <summary>
    /// The hit set of CRISPOR's wider search, shared by the MIT and the CFD alternative-PAM aggregates:
    /// the guide's main-PAM off-targets plus those on every alternative PAM that is scannable
    /// (<see cref="GetAlternativeOffTargetPams"/>, <see cref="AlternativePamSystems"/>), a site matching
    /// both the main and a wider alternative pattern counted once, and the <c>ALTPAMMINSCORE</c>
    /// pre-filter (crispor.py L4104-4122) applied to the UNWEIGHTED MIT/Hsu (or SaCas9) hit score for
    /// both aggregates, because the reference filters the hit list once, before either score is computed.
    /// </summary>
    private static List<PamAnnotatedOffTarget> CollectAlternativePamHits(
        string guideSequence, DnaSequence genome, int maxMismatches, CrisprSystemType systemType,
        double minAlternativePamHitScore, (int Start, int End)? onTarget)
    {
        if (double.IsNaN(minAlternativePamHitScore) || minAlternativePamHitScore < 0)
            throw new ArgumentOutOfRangeException(nameof(minAlternativePamHitScore),
                "The alternative-PAM hit-score cutoff must be a non-negative number.");

        var hits = CollectSystemHits(guideSequence, genome, maxMismatches, systemType, onTarget);

        var main = GetSystem(systemType);
        var altPams = GetAlternativeOffTargetPams(systemType);
        foreach (string altPam in altPams)
        {
            foreach (var altSystem in AlternativePamSystems(altPam, main.GuideLength))
            {
                foreach (var hit in CollectSystemHits(guideSequence, genome, maxMismatches, altSystem, onTarget))
                {
                    // CRISPOR runs ONE search, so a site matching both the main and a wider alternative
                    // pattern (SaCas9 NNGRRT inside NNGRRN) is a single off-target; the main scan has it.
                    if (!MatchesPam(SequenceExtensions.ToUpperAscii(hit.Pam)!, 0, main.PamSequence))
                        hits.Add(hit);
                }
            }
        }

        if (minAlternativePamHitScore > 0 && altPams.Count > 0)
        {
            // crispor.py L4108-4122: a hit whose PAM is in the expanded alternative-PAM list is kept only
            // if its UNWEIGHTED hit score (calcHitScore / calcSaHitScore, before the 0.2 factor) is
            // >= ALTPAMMINSCORE ("if otScore < minScore: ... continue"). The membership test is on the
            // alternative list FIRST, exactly as the reference orders it, so for a PAM wider than the
            // main one (SaCas9 NNGRRN contains NNGRRT) main-PAM hits are subject to the cutoff as well.
            var guide = NormalizeNucleotides(guideSequence, nameof(guideSequence));
            bool useSaModel = IsSaCas9Pam(main.PamSequence) && guide.Length == SaCas9GuideLength;
            hits = hits.Where(hit =>
            {
                string pam = SequenceExtensions.ToUpperAscii(hit.Pam)!;
                if (!altPams.Any(alt => MatchesPam(pam, 0, alt)))
                    return true;
                string ot = SequenceExtensions.ToUpperAscii(hit.Site.Sequence)!;
                double score = useSaModel ? CalculateSaCas9HitScore(guide, ot) : CalculateMitHitScoreFudged(guide, ot);
                return score >= minAlternativePamHitScore;
            }).ToList();
        }

        return hits;
    }

    /// <summary>
    /// The <see cref="CrisprSystemType"/> values whose main PAM is <paramref name="altPam"/> and whose
    /// guide length is <paramref name="guideLength"/> — the only way an alternative PAM can be scanned
    /// by this library, whose off-target search takes one system at a time. The reference's own
    /// constraint that an alternative PAM has the same length as the main one is inherited from its
    /// table; the guide-length match is required additionally so the same guide can be used for the scan.
    /// </summary>
    private static IEnumerable<CrisprSystemType> AlternativePamSystems(string altPam, int guideLength)
    {
        foreach (var type in Enum.GetValues<CrisprSystemType>())
        {
            var system = GetSystem(type);
            if (system.PamAfterTarget
                && system.GuideLength == guideLength
                && string.Equals(system.PamSequence, altPam, StringComparison.Ordinal))
            {
                yield return type;
            }
        }
    }

    /// <summary>
    /// Throws when the MIT/Hsu off-target model does not apply to <paramref name="systemType"/>.
    /// The published model is SpCas9-specific, and the reference's CODE reports no MIT (or CFD)
    /// off-target score at all for Cpf1/Cas12a-class systems: in <c>annotateOfftargets</c> every hit
    /// gets <c>mitScore = 0.0</c> and <c>cfdScore = 0.0</c> when <c>pamIsCpf1(pam)</c>, and the
    /// guide-level aggregates are the sentinels <c>guideScore = -1</c> / <c>guideCfdScore = -1</c>.
    /// Refusing loudly is the honest translation of those sentinels: silently returning a perfect
    /// specificity of 100 for a guide whose off-targets were never scored would be unsafe, and
    /// returning the sentinel −1 would break this method's documented [0, 100] range.
    /// <para>
    /// Note on the reference's own prose: its Cpf1 banner (<c>showPamWarning</c>) says "there is no
    /// off-target ranking algorithm available in the literature, to our knowledge. We use Hsu and CFD
    /// scores below for off-target ranking, but they were developed for spCas9." That second sentence
    /// is not what the code does — the scoring path zeroes both scores for Cpf1 PAMs — so the code,
    /// not the banner, is taken as the reference behaviour here.
    /// </para>
    /// <para>
    /// Divergence, recorded deliberately: the reference's <c>pamIsCpf1</c> list (crispor.py L889-891)
    /// does not contain CasX's <c>TTCN</c> — that has its own <c>pamIsCasX</c> (L897-899) — so CRISPOR
    /// routes a CasX guide into the <c>else</c> branch and scores it with <c>calcHitScore</c>. That path
    /// is demonstrably mis-oriented for a PAM-before system. <c>pamIsCasX</c> sets
    /// <c>pamIsFirst = True</c> (L483-487), so <c>concatGuideAndPam</c> builds the scored string as
    /// <c>pamSeq + guideSeq</c> (L5543-5548) with the PAM at the 5' end, yet
    /// <c>annotateOfftargets</c> strips the PAM with <c>guideSeq[:len(guideSeq)-len(pam)]</c>
    /// (L1932-1933), which removes bases from the 3' end. The string that reaches <c>calcHitScore</c>
    /// therefore still carries the 4-nt PAM at its 5' end and has lost the 3'-most 4 protospacer bases
    /// — and <c>hitScoreM</c> is indexed PAM-proximal-last, so every position weight lands on the wrong
    /// base. This library refuses all PAM-before systems rather than reproduce that.
    /// </para>
    /// </summary>
    private static void EnsureMitModelApplies(CrisprSystemType systemType, string paramName)
    {
        var system = GetSystem(systemType);
        if (!system.PamAfterTarget)
            throw new ArgumentException(
                $"The MIT/Hsu 2013 off-target model is SpCas9-specific and is not defined for " +
                $"{system.Name} (PAM-before-target / Cas12a-class). CRISPOR reports no MIT or CFD " +
                $"off-target score for these systems. Score such sites with a system-specific model instead.",
                paramName);
    }

    /// <summary>
    /// Applies the reference <c>calcHitScore</c>'s own guide-length handling before scoring: longer
    /// than 20 nt uses the last 20 ("S. aureus requires 21bp long guides. We fudge by using only last
    /// 20bp", crispor.py L2056-2059), exactly 19 nt is left-padded with 'A' ("for 19bp guides, we fudge
    /// a little, but first pos has no weight anyways", L2060-2063). Guides shorter than 19 nt are
    /// refused: the reference does have a branch for them (<c>matrixStart = 20 - len</c>,
    /// <c>maxDist = len - 1</c>, L2064-2067) but its own comment reads "for shorter guides, I'm not sure
    /// if this score makes sense anymore, we force things".
    /// <para>
    /// Reachability, stated so the branches are not mistaken for live parity: with this library's
    /// current <see cref="CrisprSystem"/> set only the <c>== 20</c> branch runs. The five
    /// PAM-after-target systems are SpCas9, SpCas9-NAG and SpCas9-NGA (20 nt each) and SaCas9 and
    /// SaCas9-NNGRRN (21 nt each), and the SaCas9 ones are
    /// routed to the Tycko 2018 model (<see cref="CalculateSaCas9HitScore"/>) by the caller, exactly as
    /// <c>annotateOfftargets</c> routes it, so the <c>&gt; 20</c> fudge is never taken; no system has a
    /// 19-nt guide and <see cref="FindOffTargets"/> rejects a guide that does not match its system's
    /// length, so the <c>== 19</c> pad is never taken either. Both are kept because the reference
    /// applies them (the 19-nt pad in <c>annotateOfftargets</c> itself, L1934-1936), so a future
    /// non-20-nt Cas9-class system is scored rather than silently dropped.
    /// </para>
    /// </summary>
    private static double CalculateMitHitScoreFudged(string guide, string offTarget)
    {
        if (guide.Length != offTarget.Length)
            throw new ArgumentException("MIT/Hsu scoring requires equal-length sequences.", nameof(offTarget));

        if (guide.Length > 20)
            return CalculateMitHitScore(guide[^20..], offTarget[^20..]);
        if (guide.Length == 20)
            return CalculateMitHitScore(guide, offTarget);
        if (guide.Length == 19)
            return CalculateMitHitScore("A" + guide, "A" + offTarget);

        throw new ArgumentException(
            $"The MIT/Hsu off-target model is not defined for guides shorter than 19 nt (got {guide.Length}).",
            nameof(guide));
    }

    #endregion

    #region SaCas9 Off-Target Score — Tycko et al. 2018

    // ---------------------------------------------------------------------------------------------
    // SaCas9-specific off-target activity model (PUBLISHED, fully reproducible).
    //
    // Source: Tycko J, Wainberg M, Marinov GK, et al. "Pairwise library screen systematically
    //   interrogates Staphylococcus aureus Cas9 specificity in human cells."
    //   Nat Commun 9, 2962 (2018). doi:10.1038/s41467-018-05391-2.
    //   (the paper CRISPOR's calcSaHitScore cites: "saCas9 offtarget scoring from Tycko et al,
    //    https://www.nature.com/articles/s41467-018-05391-2")
    //
    // Reference implementation: CRISPOR's `calcSaHitScore` (crispor.py L2025-2043) delegates to
    //   `bin/src/pairwise-library-screen/predictSingle.py` (`SaCas9Scorer`), whose trained
    //   coefficients ship as `models/rstan_mult_<guideLen>_coeffs.txt`. Both were downloaded this
    //   session from raw.githubusercontent.com/maximilianh/crisporWebsite/master/ and the table below
    //   is a verbatim transcription of the `Median` column of `rstan_mult_21_coeffs.txt`
    //   (252 rows = 21 positions x 12 mismatch types).
    //
    // ALGORITHM (verbatim from `SaCas9Scorer.calcScore`):
    //   predActivity = 1
    //   for i, (dna, rna) in enumerate(zip(targetSeq, guideSeq)):
    //       if i == 0: continue          # "Mismatches at the 5' position are not truly determined
    //                                   #  by the model"
    //       if dna != rna: predActivity *= coeff[ 'd'+complement(dna) + ':r' + rnaBase(rna) ][i]
    //   calcSaHitScore = 100 * predActivity   # "to be compatible with the MIT score, has to be in
    //                                         #  the range 0-100 for the MIT aggregate guide score"
    //
    // KEY CONVENTION (pinned from the source): `mm_dict` maps the raw pair (targetBase, guideBase) to
    //   the coefficient id, e.g. ('A','C') -> "dT:rC" and ('A','T') -> "dT:rU". So 'd' is the
    //   COMPLEMENT of the off-target (DNA) base — the base on the non-target strand that pairs against
    //   the guide — and 'r' is the GUIDE (RNA) base with T written as U. Same convention as CFD, only
    //   written in the opposite order ("dX:rY" vs CFD's "rY:dX").
    //
    // POSITION CONVENTION (pinned from the source): the shipped tables index `Position` from the
    //   3' / PAM-proximal end, and `SaCas9Scorer.__init__` remaps it with
    //   `guideLenDict[guide_len - Position + 1]`, i.e. to a 1-based index from the 5' end. The arrays
    //   below are already in that remapped order: index 0 = 5'-most base (never read, see above),
    //   index 20 = PAM-proximal base. Mis-ordering is caught by the orientation guard test.
    //
    // Cross-checks reproduced this session (independent Python run of predictSingle.py's own
    // semantics against the downloaded coefficient table, i.e. the reference's OWN asserted values):
    //   calcScore("GCAACCACAAACCCACGAGGG", "GCAACCACAAACCCACGAGGG") = 1.0
    //   calcScore("GCAACCACAAACCCACGAGGG", "ACAAACACATACCCACAAGGA") = 0.041293256150042444
    //   calcScore("GGGTGAGTGAGTGTGTGCGTG", "GAGCGAGTGGGTGTGTGCGTG") = 0.13441983202750601
    //   => calcSaHitScore (x100) = 100.0, 4.1293256150042446, 13.4419832027506
    // ---------------------------------------------------------------------------------------------

    /// <summary>The guide length of the SaCas9 system, and the only length this model is shipped for here.</summary>
    private const int SaCas9GuideLength = 21;

    /// <summary>
    /// Trained SaCas9 mismatch coefficients (Tycko et al. 2018), keyed by mismatch type
    /// <c>"dX:rY"</c> (X = complement of the off-target base; Y = guide base with T as U), each value
    /// an array indexed 0..20 by 0-based protospacer position (index 0 = 5' end, index 20 =
    /// 3'/PAM-proximal). Verbatim <c>Median</c> column of the reference's <c>rstan_mult_21_coeffs.txt</c>.
    /// </summary>
    private static readonly Dictionary<string, double[]> SaCas9MismatchCoefficients = new()
    {
        ["dT:rC"] = new[] { 0.500182335832017, 0.546734106002513, 0.547691877828188, 0.552964194291369, 0.538528380008684, 0.325610215099261, 0.576468070360727, 0.127378670053041, 0.633439692869058, 0.615658340132625, 0.487865019225753, 0.42689462494337, 0.102397715418986, 0.200257694187704, 0.136499862475113, 0.0493188093803229, 0.136539686900136, 0.0927354706350068, 0.130261909041357, 0.0773508496492739, 0.105693733152788 },
        ["dT:rG"] = new[] { 0.889800314629899, 0.992970435672823, 0.732601715905163, 0.980427129883399, 0.685646931401168, 0.844922889454443, 0.667118034664953, 0.683173864400938, 0.555513095453494, 0.520509986734207, 0.457991759790709, 0.612354705406714, 0.37824405562458, 0.395337246155066, 0.27104610293669, 0.451540654839527, 0.443502321156338, 0.584714691584839, 0.436629040331668, 0.22354371179448, 0.445155052527266 },
        ["dT:rU"] = new[] { 0.507135821278318, 0.657829492961583, 0.636416873002814, 0.363267257971301, 0.576799741012485, 0.642009817365257, 0.567712358758302, 0.599538413152022, 0.506919342187349, 0.646750045500538, 0.455951859040721, 0.456786357841618, 0.267574282087781, 0.12373852758353, 0.0789972001603137, 0.10911606478632, 0.199571974143948, 0.16959064954525, 0.296554470273688, 0.0601305506668423, 0.0862142526234551 },
        ["dG:rA"] = new[] { 0.507585860878205, 0.522967414281177, 0.64726672914593, 0.178476371392183, 0.589826846023653, 0.445331765504143, 0.685794931847784, 0.325409595348562, 0.64243781038661, 0.224012343553474, 0.217821657175021, 0.0642029910217309, 0.0615709293577342, 0.0162972303713068, 0.0226050940568139, 0.0597797023253567, 0.050358810530371, 0.0775368431850832, 0.0152863394236489, 0.1432549909342, 0.0426130973647857 },
        ["dG:rG"] = new[] { 0.482584240469274, 0.643567704594203, 0.681666229047031, 0.416784632985585, 0.495035619290045, 0.531756470458885, 0.676004709924651, 0.458277575190559, 0.575785506970191, 0.444440796592595, 0.387918368049248, 0.216680572369354, 0.228139311830686, 0.0287977020564966, 0.128405660784272, 0.191384782904834, 0.0357634305414073, 0.0367789926238243, 0.070695309584881, 0.158363763886502, 0.104525045604611 },
        ["dG:rU"] = new[] { 0.49520183259973, 0.618178833577651, 0.581015956193127, 0.48725485300852, 0.942657171400989, 0.540642115587474, 0.873378872726365, 0.678015697554753, 0.58682811500121, 0.686634819793447, 0.431002446808513, 0.261650340816763, 0.317651503071752, 0.254249722211379, 0.187839704652641, 0.0979027817488908, 0.236699959495542, 0.113816114093481, 0.0573652386585105, 0.157638078782083, 0.179929057922201 },
        ["dC:rA"] = new[] { 0.497493947658051, 0.46259280684568, 0.643038591122457, 0.372212865767711, 0.680737138305581, 0.552943486205603, 0.584975615068689, 0.582110823409603, 0.495378197462189, 0.27782470178935, 0.307880431933795, 0.299394927832334, 0.141962480894423, 0.0292199017744031, 0.0195124051391844, 0.140212678240853, 0.107550670622112, 0.0533316085854886, 0.075377597399379, 0.0460670515595503, 0.109313834975141 },
        ["dC:rC"] = new[] { 0.502985937567947, 0.15914318618166, 0.274409757488698, 0.289442475421799, 0.300672177947638, 0.368351594652677, 0.349878390379607, 0.103481082335252, 0.529851864152853, 0.43201684115079, 0.155692052708749, 0.189216008110395, 0.0529315895573772, 0.039536033512892, 0.0232766388217564, 0.0191474738800195, 0.025852157724115, 0.0148089478006535, 0.0356756456096336, 0.0154733404301811, 0.0449961433951072 },
        ["dC:rU"] = new[] { 0.503187216328694, 0.490930802795924, 0.366150894868735, 0.518288493112618, 0.796643503923613, 0.530248054175443, 0.645199638303918, 0.51997749852277, 0.326200207604272, 0.406797059100196, 0.189198509598849, 0.23819640621691, 0.0970060141425968, 0.0302548103084032, 0.0307495242144697, 0.0522509311690574, 0.0575108872220017, 0.0400655724089505, 0.0554196638437215, 0.0535733539298569, 0.0805457311950328 },
        ["dA:rA"] = new[] { 0.497855206412811, 0.624360922159452, 0.698428615250848, 0.277678942354815, 0.38993398748401, 0.422521983677337, 0.557841348776102, 0.104010305290341, 0.441491391101617, 0.388385789116283, 0.507615988748957, 0.0630783302069826, 0.0388099660522349, 0.0460124390086707, 0.0961446136331311, 0.111029408600699, 0.056782864029314, 0.093036719073791, 0.059627368313025, 0.218904216245295, 0.0638608590837734 },
        ["dA:rC"] = new[] { 0.499950169216313, 0.640457481900362, 0.60640601604755, 0.629270731339457, 0.390709495232679, 0.280235604556182, 0.717284514586225, 0.41657664122309, 0.576259924027269, 0.646131400505293, 0.555381735098051, 0.251895624458276, 0.101076231082529, 0.346568359640956, 0.178707425151124, 0.117130325394167, 0.117154009003686, 0.0979767483566834, 0.136926154119441, 0.191997864741009, 0.043097315653776 },
        ["dA:rG"] = new[] { 0.763319435742999, 0.703936911913447, 0.475908917669472, 0.168587534931865, 0.515932076149009, 0.503183576570521, 0.565247609715247, 0.48471892642459, 0.56059950527931, 0.46065379532187, 0.374430944727598, 0.23572096660813, 0.0794278765493476, 0.0553640525834663, 0.0409451753782538, 0.222056916803022, 0.0586170001847033, 0.0832194891199737, 0.0372902912891327, 0.122380573189007, 0.177723354855625 },
    };

    /// <summary>
    /// Returns true for the PAM motifs the reference treats as SaCas9 (<c>pamIsSaCas9</c>:
    /// <c>NNGRRT</c> / <c>NNNRRT</c>, matched on the part before any <c>-</c> suffix). This library's
    /// <see cref="CrisprSystemType.SaCas9_NNGRRN"/> (the reference's alternative PAM of <c>NNGRRT</c>,
    /// never a main PAM there) is a SaCas9 system too and takes the same Tycko model.
    /// </summary>
    private static bool IsSaCas9Pam(string pam)
    {
        int dash = pam.IndexOf('-', StringComparison.Ordinal);
        string head = dash < 0 ? pam : pam[..dash];
        return head is "NNGRRT" or "NNNRRT" or "NNGRRN";
    }

    /// <summary>
    /// Calculates the SaCas9-specific single-hit off-target activity score of Tycko et al. 2018
    /// (Nat Commun 9:2962) for a 21-nt SaCas9 guide against a 21-nt off-target protospacer, on the
    /// 0-100 scale the reference uses so the value can enter the MIT aggregate
    /// (<see cref="CalculateMitSpecificityScore(IEnumerable{double})"/>). A perfect match scores
    /// exactly 100; each mismatch multiplies the activity by that position's trained coefficient.
    /// </summary>
    /// <param name="saGuide21">The 21-nt SaCas9 guide (PAM excluded), 5' first. A/C/G/T, case-insensitive.</param>
    /// <param name="offTarget21">The 21-nt off-target protospacer (PAM excluded), same orientation.</param>
    /// <returns>The single-hit score in [0, 100]; higher = more likely an active off-target.</returns>
    /// <remarks>
    /// Source: the trained model distributed with CRISPOR as
    /// <c>bin/src/pairwise-library-screen/predictSingle.py</c> + <c>models/rstan_mult_21_coeffs.txt</c>,
    /// applied by <c>crispor.py</c>'s <c>calcSaHitScore</c> as <c>100 × calcScore(...)</c>.
    /// Two faithful details of the reference are reproduced here: a mismatch at the 5'-most position
    /// contributes nothing ("Mismatches at the 5' position are not truly determined by the model"), and
    /// the model does not support gaps/bulges (they are refused).
    /// <para>
    /// Scope: the reference ships coefficient tables for guide lengths 20-24; only the 21-nt table is
    /// transcribed here, because 21 nt is the guide length of this library's only SaCas9 system
    /// (<see cref="CrisprSystemType.SaCas9"/>). Other lengths are refused rather than approximated.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentNullException">Either sequence is null or empty.</exception>
    /// <exception cref="ArgumentException">
    /// Either sequence is not 21 nt, or contains a character outside A/C/G/T (including the gap
    /// character the reference explicitly refuses).
    /// </exception>
    public static double CalculateSaCas9HitScore(string saGuide21, string offTarget21)
    {
        if (string.IsNullOrEmpty(saGuide21))
            throw new ArgumentNullException(nameof(saGuide21));
        if (string.IsNullOrEmpty(offTarget21))
            throw new ArgumentNullException(nameof(offTarget21));
        if (saGuide21.Length != SaCas9GuideLength)
            throw new ArgumentException(
                $"The Tycko 2018 SaCas9 model is transcribed here for a {SaCas9GuideLength}-nt guide; got {saGuide21.Length}.",
                nameof(saGuide21));
        if (offTarget21.Length != SaCas9GuideLength)
            throw new ArgumentException(
                $"The Tycko 2018 SaCas9 model requires a {SaCas9GuideLength}-nt off-target protospacer; got {offTarget21.Length}.",
                nameof(offTarget21));

        string guide = SequenceExtensions.ToUpperAscii(saGuide21)!;
        string off = SequenceExtensions.ToUpperAscii(offTarget21)!;

        double predictedActivity = 1.0;

        // The reference skips index 0 but still validates every base, so invalid input is refused
        // rather than silently scored. Unlike the MIT hit score, this scorer keeps its own strict
        // A/C/G/T check (no U -> T reading), as the CFD scorer does (its reference screens [^ATCG]).
        for (int i = 0; i < SaCas9GuideLength; i++)
        {
            char rna = guide[i];
            char dna = off[i];

            if (rna is not ('A' or 'C' or 'G' or 'T'))
                throw new ArgumentException($"The SaCas9 model requires A/C/G/T sequences; found '{rna}'.", nameof(saGuide21));
            if (dna is not ('A' or 'C' or 'G' or 'T'))
                throw new ArgumentException($"The SaCas9 model requires A/C/G/T sequences; found '{dna}'.", nameof(offTarget21));

            // "Mismatches at the 5' position are not truly determined by the model" (predictSingle.py).
            if (i == 0 || rna == dna)
                continue;

            // 'd' = complement of the off-target DNA base; 'r' = guide base with T written as U.
            string key = $"d{CfdComplement(dna)}:r{(rna == 'T' ? 'U' : rna)}";
            predictedActivity *= SaCas9MismatchCoefficients[key][i];
        }

        // "to be compatible with the MIT score, has to be in the range 0-100" (calcSaHitScore).
        return 100.0 * predictedActivity;
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

    /// <summary>
    /// Strict A/C/G/T complement for the CFD key construction. Deliberately NOT the canonical
    /// <see cref="SequenceExtensions.GetComplementBase(char)"/>: that one is IUPAC-complete, passes
    /// non-IUPAC characters through unchanged and maps U→A, whereas CFD must REFUSE anything outside
    /// A/C/G/T — the reference <c>calcCfdScore</c> screens its inputs with <c>re.search('[^ATCG]')</c>
    /// and returns the sentinel -1 rather than scoring them. Delegating would turn an ambiguity code
    /// into either an opaque missing-key failure or (for U) a silently wrong score, so this strict
    /// validator is kept by design.
    /// </summary>
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

        var guide = SequenceExtensions.ToUpperAscii(sgRna20)!;
        var off = SequenceExtensions.ToUpperAscii(offTarget20)!;
        var pam = SequenceExtensions.ToUpperAscii(offTargetPam)!;

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

    /// <summary>
    /// Aggregates per-off-target CFD scores into a single guide-level CFD specificity score, using the
    /// same <c>100 / (100 + Σ hits) × 100</c> form as the MIT aggregate — which is exactly what the
    /// reference does by default (<c>calcMitGuideScore(sum(cfdScores))</c>) — over the RAW CFD values
    /// in [0, 1], i.e. the CFD scores are NOT rescaled to 0-100 first.
    /// </summary>
    /// <param name="cfdScores">The per-off-target CFD scores (each in [0, 1], as returned by <see cref="CalculateCfdScore"/>).</param>
    /// <param name="aggregation">Which published aggregation to use; see <see cref="CfdAggregation"/>.</param>
    /// <returns>The aggregate CFD specificity; 100 for an empty off-target set under the default aggregation.</returns>
    /// <remarks>
    /// Source: CRISPOR <c>crispor.py</c> <c>annotateOfftargets</c>, which computes the guide's CFD
    /// specificity as <c>calcMitGuideScore(sum(cfdScores))</c> unless the <c>doCfdFix</c> flag (URL
    /// argument <c>fixCfd=1</c>) is set. Under that flag the reference applies TWO changes together:
    /// <c>calcCfdScore</c> rescales every per-hit CFD by ×100, and the aggregate becomes
    /// <c>calcCfdGuideScore(sum) = 100 × 100 / sum</c>. Over a sum of RAW per-hit CFD values the two
    /// factors therefore partly cancel: <c>100 × 100 / (100 × Σ raw) = 100 / Σ raw</c>, which is what
    /// <see cref="CfdAggregation.ParkinsonNormalized"/> computes here (this method always takes raw
    /// per-hit CFD values in [0, 1], for both aggregations). The reference rounds the default aggregate
    /// to an integer (<c>int(round(score × 100))</c>); this method returns the unrounded value so
    /// callers may round as they wish.
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown for <see cref="CfdAggregation.ParkinsonNormalized"/> when the CFD sum is 0: that formula
    /// divides by the sum and the reference is undefined there (it would raise <c>ZeroDivisionError</c>).
    /// </exception>
    public static double CalculateCfdSpecificityScore(
        IEnumerable<double> cfdScores,
        CfdAggregation aggregation = CfdAggregation.MitStyle)
    {
        ArgumentNullException.ThrowIfNull(cfdScores);
        double sum = cfdScores.Sum();

        return aggregation switch
        {
            CfdAggregation.MitStyle => 100.0 / (100.0 + sum) * 100.0,
            // The reference's doCfdFix path is 100 x 100 / Σ(100 x rawCfd) = 100 / Σ rawCfd: its per-hit
            // x100 rescaling cancels one of the two factors of 100. Inputs here are always raw CFD.
            CfdAggregation.ParkinsonNormalized when sum > 0 => 100.0 / sum,
            CfdAggregation.ParkinsonNormalized => throw new ArgumentOutOfRangeException(
                nameof(cfdScores),
                "The ParkinsonNormalized CFD aggregate (100 / sum of the raw CFD scores) is undefined for a CFD sum of 0."),
            _ => throw new ArgumentOutOfRangeException(nameof(aggregation))
        };
    }

    /// <summary>
    /// Calculates the aggregate CFD specificity score for a guide against a genome by enumerating
    /// off-target sites (via <see cref="FindOffTargets"/>) and combining their CFD scores. Exact
    /// on-target matches are excluded (they are not off-targets).
    /// </summary>
    /// <param name="guideSequence">The 20-nt guide sequence.</param>
    /// <param name="genome">The genome/sequence to scan for off-targets.</param>
    /// <param name="maxMismatches">Maximum mismatches to consider (0-5, as in <see cref="FindOffTargets"/>).</param>
    /// <param name="systemType">CRISPR system; must be a 20-nt PAM-after-target (Cas9-class) system.</param>
    /// <param name="aggregation">Which published aggregation to use.</param>
    /// <param name="onTargetStart">
    /// Optional on-target interval, start (0-based, inclusive, forward-strand coordinates of
    /// <paramref name="genome"/>; both bounds or neither). See the remark on the on-target interval.
    /// </param>
    /// <param name="onTargetEnd">Optional on-target interval, end (exclusive); see <paramref name="onTargetStart"/>.</param>
    /// <returns>The aggregate CFD specificity; 100 when no off-targets are found (default aggregation).</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="maxMismatches"/> is outside 0-5 (the same bound
    /// <see cref="FindOffTargets"/> enforces), or when <paramref name="aggregation"/> is
    /// <see cref="CfdAggregation.ParkinsonNormalized"/> and the CFD sum is 0.
    /// </exception>
    /// <remarks>
    /// Each off-target's CFD score needs that site's own PAM, which the scan reports on the protospacer
    /// strand. CFD is defined only for a 20-nt protospacer, so systems with a different guide length are
    /// refused rather than fudged (unlike the MIT model, the reference provides no length fallback for
    /// CFD); that refusal also covers SaCas9, for which the reference itself emits the per-hit CFD
    /// sentinel (<c>cfdScore = -1</c>, crispor.py L1943-1945). Note what the reference's guide-level
    /// value then is for SaCas9: every hit is dropped by its own <c>if cfdScore != -1</c> filter, so
    /// <c>guideCfdScore = calcMitGuideScore(sum([])) = 100</c> — a vacuous perfect score rather than a
    /// sentinel (unlike Cpf1, which does get <c>guideCfdScore = -1</c>). Refusing is the honest
    /// translation of that: a 100 for a guide whose off-targets were never CFD-scored would be unsafe.
    /// <para>
    /// Off-targets whose protospacer is not pure A/C/G/T are skipped rather than scored, mirroring the
    /// reference (<c>calcCfdScore</c> screens both sequences with <c>re.search('[^ATCG]')</c> and
    /// returns the sentinel -1; <c>annotateOfftargets</c> appends the hit only <c>if cfdScore != -1</c>).
    /// That guard is <b>defensive only and currently unreachable</b> from this overload: a
    /// <see cref="DnaSequence"/> is validated A/C/G/T at construction, so no ambiguous base can ever
    /// reach the scan. It is kept so the reference's exclusion rule is encoded rather than implied, and
    /// so the method stays correct if the genome type ever admits ambiguity codes.
    /// </para>
    /// <para>
    /// The optional on-target interval works as in
    /// <see cref="CalculateMitSpecificityScore(string, DnaSequence, int, CrisprSystemType, int?, int?)"/>
    /// (an exact off-target contributes its PAM score, <c>calcCfdScore</c> of a perfect match).
    /// </para>
    /// </remarks>
    public static double CalculateCfdSpecificityScore(
        string guideSequence,
        DnaSequence genome,
        int maxMismatches = 4,
        CrisprSystemType systemType = CrisprSystemType.SpCas9,
        CfdAggregation aggregation = CfdAggregation.MitStyle,
        int? onTargetStart = null,
        int? onTargetEnd = null)
    {
        var guide = ValidateCfdGenomeArguments(guideSequence, genome, maxMismatches, systemType);
        var onTarget = ResolveOnTargetInterval(onTargetStart, onTargetEnd);

        var cfdScores = new List<double>();
        foreach (var (site, pam) in ScanOffTargetsCore(guide, genome, maxMismatches, systemType, onTarget))
        {
            string offTarget = SequenceExtensions.ToUpperAscii(site.Sequence)!;

            // The reference screens both sequences for [^ATCG] and returns -1, which annotateOfftargets
            // then leaves out of the CFD sum. Defensive only: DnaSequence validates A/C/G/T at
            // construction (SequenceExtensions.IndexOfInvalidDna), so this never fires today.
            if (!offTarget.AsSpan().IsValidDna())
                continue;

            cfdScores.Add(CalculateCfdScore(guide, offTarget, pam));
        }

        return CalculateCfdSpecificityScore(cfdScores, aggregation);
    }

    /// <summary>
    /// The argument guards shared by the genome CFD aggregates (same order and exception types as
    /// <see cref="FindOffTargets"/>, which they cannot delegate to because CFD needs each hit's PAM);
    /// returns the normalised guide.
    /// </summary>
    private static string ValidateCfdGenomeArguments(
        string guideSequence, DnaSequence genome, int maxMismatches, CrisprSystemType systemType)
    {
        // Validate in the same order and with the same exception types as FindOffTargets, which this
        // overload cannot delegate to (it needs each hit's PAM, so it drives the PAM-annotated scan
        // directly and must therefore repeat FindOffTargets' own guards).
        if (string.IsNullOrEmpty(guideSequence))
            throw new ArgumentNullException(nameof(guideSequence));
        ArgumentNullException.ThrowIfNull(genome);
        if (maxMismatches < 0 || maxMismatches > 5)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches));

        EnsureMitModelApplies(systemType, nameof(systemType));
        var system = GetSystem(systemType);
        if (system.GuideLength != 20)
            throw new ArgumentException(
                $"The CFD model is defined only for a 20-nt protospacer; {system.Name} uses {system.GuideLength} nt.",
                nameof(systemType));

        if (guideSequence.Length != system.GuideLength)
            throw new ArgumentException(
                $"Guide length ({guideSequence.Length}) does not match expected length ({system.GuideLength}) for {system.Name}.",
                nameof(guideSequence));
        return NormalizeNucleotides(guideSequence, nameof(guideSequence));
    }

    /// <summary>
    /// Calculates the aggregate CFD specificity score for a guide against a genome over the guide's own
    /// PAM <b>and</b> the reference's alternative PAMs for it (<see cref="GetAlternativeOffTargetPams"/>:
    /// <c>NAG</c> and <c>NGA</c> for SpCas9). This is CRISPOR's guide CFD score: <c>annotateOfftargets</c>
    /// runs <c>calcCfdScore(guideSeq, otSeq)</c> on every hit of its one wider search, so a guide whose
    /// only off-target sits on <c>NAG</c> is not "perfectly specific" as
    /// <see cref="CalculateCfdSpecificityScore(string, DnaSequence, int, CrisprSystemType, CfdAggregation, int?, int?)"/>
    /// (main PAM only) reports it.
    /// </summary>
    /// <param name="guideSequence">The 20-nt guide sequence.</param>
    /// <param name="genome">The genome/sequence to scan for off-targets.</param>
    /// <param name="maxMismatches">Maximum mismatches to consider (0-5, as in <see cref="FindOffTargets"/>).</param>
    /// <param name="systemType">The guide's own system; a 20-nt PAM-after-target (Cas9-class) system.</param>
    /// <param name="aggregation">Which published aggregation to use.</param>
    /// <param name="minAlternativePamHitScore">
    /// The reference's <c>ALTPAMMINSCORE</c> (default <c>1.0</c>; <c>0</c> = no filtering), a minimum
    /// <b>MIT/Hsu</b> hit score (0-100 scale) for an alternative-PAM hit — see the remarks.
    /// </param>
    /// <param name="onTargetStart">
    /// Optional on-target interval, start (0-based, inclusive, forward-strand coordinates of
    /// <paramref name="genome"/>; both bounds or neither). See the remark on the on-target interval.
    /// </param>
    /// <param name="onTargetEnd">Optional on-target interval, end (exclusive); see <paramref name="onTargetStart"/>.</param>
    /// <returns>The aggregate CFD specificity; 100 when no off-targets are found (default aggregation).</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="maxMismatches"/> is outside 0-5, the cutoff is negative/NaN, or the
    /// <see cref="CfdAggregation.ParkinsonNormalized"/> CFD sum is 0.
    /// </exception>
    /// <remarks>
    /// Source: crispor.py. The hit list is pre-filtered once (L4085-4135): a hit whose PAM is in the
    /// expanded alternative list needs <c>calcHitScore(guideSeq, tSeqNoPam) &gt;= ALTPAMMINSCORE</c> —
    /// the MIT/Hsu hit score, not the CFD score, even though this aggregate later sums CFD values — and
    /// <c>annotateOfftargets</c> then appends <c>calcCfdScore(guideSeq, otSeq)</c> for each surviving
    /// hit. <b>No 0.2 factor applies on the CFD path</b>: the reference multiplies only the MIT score
    /// (<c>if pam=="NGG" and otSeq[-2:]!="GG": mitScore = mitScore * 0.2</c>, L1955), and takes
    /// <c>cfdScore</c> afterwards from the full <c>otSeq</c> including its PAM (the CFD PAM table, e.g.
    /// <c>AG</c> = 0.2593, is what penalises the non-NGG PAM). The hit collection (alternative-PAM list,
    /// scannable systems, main/alternative dedup, cutoff) is the one of
    /// <see cref="CalculateMitSpecificityScoreWithAlternativePams"/>. The CFD model needs a 20-nt
    /// protospacer, so SaCas9 is refused exactly as in the main-PAM CFD overload. The optional on-target
    /// interval works as in <see cref="CalculateMitSpecificityScoreWithAlternativePams"/>.
    /// </remarks>
    public static double CalculateCfdSpecificityScoreWithAlternativePams(
        string guideSequence,
        DnaSequence genome,
        int maxMismatches = 4,
        CrisprSystemType systemType = CrisprSystemType.SpCas9,
        CfdAggregation aggregation = CfdAggregation.MitStyle,
        double minAlternativePamHitScore = 1.0,
        int? onTargetStart = null,
        int? onTargetEnd = null)
    {
        var guide = ValidateCfdGenomeArguments(guideSequence, genome, maxMismatches, systemType);
        var onTarget = ResolveOnTargetInterval(onTargetStart, onTargetEnd);

        var cfdScores = new List<double>();
        foreach (var hit in CollectAlternativePamHits(
                     guide, genome, maxMismatches, systemType, minAlternativePamHitScore, onTarget))
        {
            string offTarget = SequenceExtensions.ToUpperAscii(hit.Site.Sequence)!;
            if (!offTarget.AsSpan().IsValidDna())
                continue; // the reference's [^ATCG] -> -1 -> left out (defensive, as in the main-PAM overload)

            cfdScores.Add(CalculateCfdScore(guide, offTarget, hit.Pam));
        }

        return CalculateCfdSpecificityScore(cfdScores, aggregation);
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

        var fwd = NormalizeNucleotides(guideSequence, nameof(guideSequence));
        if (fwd.Length <= stemLength)
            return 0;

        string rvs = DnaSequence.GetReverseComplementString(fwd);
        int l = fwd.Length - stemLength - 1;

        // CHOPCHOP reverse-complements the backbone regions it is given before matching against them.
        var backbones = new List<string>();
        if (backboneRegions is not null)
        {
            foreach (var region in backboneRegions)
            {
                if (region is null)
                    throw new ArgumentException("Backbone regions contain a null element.", nameof(backboneRegions));
                backbones.Add(DnaSequence.GetReverseComplementString(
                    NormalizeNucleotides(region, nameof(backboneRegions))));
            }
        }

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

    /// <summary>
    /// Points deducted from <see cref="GuideRnaCandidate.Score"/> per self-complementary 4-bp stem: CHOPCHOP's
    /// <c>SCORE['FOLDING'] = 1</c> (<c>chopchop.py</c>, <c>SCORE</c> dictionary).
    /// </summary>
    public const double SelfComplementarityStemPenalty = 1.0;

    /// <summary>
    /// CHOPCHOP's flat penalty for a guide whose GC content is outside the GC window:
    /// <c>SCORE['CRISPR_BAD_GC'] = 300</c> (<c>chopchop.py</c>, <c>SCORE</c> dictionary), added once in
    /// <c>calcGCContent</c> when <c>GCcontent &gt; GC_HIGH or GCcontent &lt; GC_LOW</c>. It exceeds the 0–100 range of
    /// <see cref="GuideRnaCandidate.Score"/>, so a guide outside the window scores 0 when
    /// <see cref="GuideRnaParameters.ScoreGcContent"/> is on.
    /// </summary>
    public const double BadGcContentPenalty = 300.0;

    /// <summary>
    /// The fraction of the 4-bp windows CHOPCHOP examines (<c>len − 4</c>, its <c>range(0, len(fwd) − STEM_LEN)</c>)
    /// that carry a self-complementary stem: <c>stems / (len − 4)</c>, in [0, 1]. The numerator is CHOPCHOP's
    /// <c>folding</c>; only the division by its own loop bound is this library's.
    /// </summary>
    private static double SelfComplementarityWindowFraction(int length, int stems)
    {
        int windows = length - 4;
        return windows > 0 ? (double)stems / windows : 0.0;
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
    CasX,
    /// <summary>
    /// SpCas9 with NGA PAM — the reference's second alternative off-target PAM of <c>NGG</c> and the
    /// main PAM of its <c>NGA → NGG</c> table row. Appended last so existing numeric values do not change.
    /// </summary>
    SpCas9_NGA,
    /// <summary>
    /// SaCas9 with the relaxed <c>NNGRRN</c> PAM — the reference's alternative off-target PAM of
    /// <c>NNGRRT</c> (21-nt guide). Appended last so existing numeric values do not change.
    /// </summary>
    SaCas9_NNGRRN
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
/// Selects how per-off-target CFD scores are aggregated into a guide-level CFD specificity score.
/// </summary>
/// <remarks>
/// Source: CRISPOR <c>crispor.py</c>. Both forms are present in the reference; the default there is
/// <see cref="MitStyle"/>, with <see cref="ParkinsonNormalized"/> selected by the <c>doCfdFix</c> flag
/// (URL argument <c>fixCfd=1</c>).
/// </remarks>
public enum CfdAggregation
{
    /// <summary>
    /// The reference default: <c>100 / (100 + Σ CFD) × 100</c> over the raw CFD values in [0, 1]
    /// (<c>calcMitGuideScore(sum(cfdScores))</c>). Yields 100 for an empty off-target set.
    /// </summary>
    MitStyle = 0,

    /// <summary>
    /// The <c>doCfdFix</c> variant, which makes CFD more comparable to the MIT score's magnitude:
    /// <c>100 / Σ CFD</c> over the raw CFD values in [0, 1] that this API always takes. The reference
    /// writes it as <c>calcCfdGuideScore(sum) = 100 × 100 / sum</c> ("suggested by Nicholas
    /// Parkinson") but under the same <c>doCfdFix</c> flag it also rescales every per-hit CFD by ×100
    /// inside <c>calcCfdScore</c>, so <c>100 × 100 / Σ(100 × raw) = 100 / Σ raw</c> — writing
    /// <c>100 × 100 / Σ raw</c> here would be a factor-100 error. Undefined for a CFD sum of 0 and not
    /// bounded by 100.
    /// </summary>
    ParkinsonNormalized = 1
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
/// <param name="MinGcContent">Lower edge of the GC window in percent (CHOPCHOP <c>GC_LOW</c> = 40).</param>
/// <param name="MaxGcContent">Upper edge of the GC window in percent (CHOPCHOP <c>GC_HIGH</c> = 70).</param>
/// <param name="MinScore">
/// Minimum <see cref="GuideRnaCandidate.Score"/> for <see cref="CrisprDesigner.DesignGuideRnas"/> to return a guide.
/// <b>Not from any tool</b> (CHOPCHOP and CRISPOR have no score threshold; their filters are the GC window
/// <c>filterGCmin/filterGCmax</c> and <c>filterSelfCompMax</c>, the latter being <see cref="MaxSelfComplementaryStems"/>).
/// With <see cref="ScoreGcContent"/> on, the score is <c>100 − stems</c> inside the window and 0 outside, so for a 20-nt
/// guide (at most 16 stems) every threshold from 1 to 84 selects exactly the guides inside the GC window; the default 50
/// is a value in that range, kept so existing callers see the same behaviour.
/// </param>
/// <param name="AvoidPolyT">Report <c>TTTT</c> in the issues list (no score effect).</param>
/// <param name="CheckSelfComplementarity">Apply CHOPCHOP's stem penalty (<c>-SC</c> off when <c>false</c>).</param>
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
    /// Whether a GC content outside <see cref="MinGcContent"/>..<see cref="MaxGcContent"/> adds CHOPCHOP's flat
    /// <see cref="CrisprDesigner.BadGcContentPenalty"/> (300) to the penalty, i.e. drops
    /// <see cref="GuideRnaCandidate.Score"/> to 0. Default <c>true</c>. CHOPCHOP itself ships this term <b>off</b> in
    /// every mode (<c>SCORE_GC = False</c>), so <c>false</c> reproduces CHOPCHOP's default; a GC outside the window is
    /// then only reported in <see cref="GuideRnaCandidate.Issues"/>.
    /// </summary>
    public bool ScoreGcContent { get; init; } = true;

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
/// <param name="Sequence">The protospacer / spacer, upper case, PAM excluded.</param>
/// <param name="Position">0-based position in the scanned sequence, or <c>-1</c> for a standalone evaluation.</param>
/// <param name="IsForwardStrand">Whether the guide was read from the forward strand.</param>
/// <param name="GcContent">GC content of the guide in percent.</param>
/// <param name="SeedGcContent">GC content of the 10-nt PAM-proximal seed region in percent.</param>
/// <param name="HasPolyT">Whether the guide contains <c>TTTT</c> (Pol III terminator).</param>
/// <param name="SelfComplementarityScore">
/// <see cref="SelfComplementaryStems"/> divided by the number of 4-bp windows CHOPCHOP examines (<c>length − 4</c>),
/// in [0, 1]. Derived from the published stem count (the numerator is CHOPCHOP's <c>folding</c>; the normalisation by
/// its own loop bound is this library's); it replaces the earlier unsourced complementary-pair fraction
/// (B08 round 3 G13). Not thermodynamic and not a probability.
/// </param>
/// <param name="Score">
/// Composite 0–100 quality number: <c>clamp(100 − penalty, 0, 100)</c> where <c>penalty</c> is CHOPCHOP's own sum of
/// <see cref="CrisprDesigner.SelfComplementarityStemPenalty"/> (1) per stem and, when
/// <see cref="GuideRnaParameters.ScoreGcContent"/>, <see cref="CrisprDesigner.BadGcContentPenalty"/> (300) for a GC
/// outside the window. TTTT, seed GC and cloning sites are reported in <paramref name="Issues"/> only. See
/// <see cref="CrisprDesigner.EvaluateGuideRna(string, CrisprSystemType, GuideRnaParameters?)"/>.
/// Use <see cref="OnTargetScore"/> (Rule Set 2) for a published efficacy estimate.
/// </param>
/// <param name="Issues">Human-readable reasons the score was reduced or the guide flagged.</param>
/// <param name="System">The CRISPR system the guide was evaluated for.</param>
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
    /// sourced measure and the one that drives <see cref="Score"/>
    /// (<see cref="SelfComplementarityScore"/> is this count divided by <c>length − 4</c>).
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

/// <summary>
/// An off-target site together with the PAM it carries, as read on the protospacer strand
/// (<see cref="PamSite.PamSequence"/>). Scores that are defined on the off-target's own PAM take this:
/// the CFD model and the alternative-PAM-aware MIT aggregate
/// (<see cref="CrisprDesigner.CalculateMitSpecificityScore(string, IEnumerable{PamAnnotatedOffTarget}, CrisprSystemType)"/>).
/// </summary>
/// <param name="Site">The off-target site.</param>
/// <param name="Pam">The PAM sequence at that site, on the protospacer strand.</param>
public sealed record PamAnnotatedOffTarget(OffTargetSite Site, string Pam);
