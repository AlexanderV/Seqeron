using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.MolTools;

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
            MinLength: 50, MaxLength: 60,
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
    /// A single gapped (Smith-Waterman) hit of a probe against a reference sequence.
    /// Unlike the ungapped Hamming scan used by <see cref="ValidateProbe"/>, a hit may span an
    /// insertion or deletion (see <see cref="HasGaps"/>), so off-target sites reachable only through
    /// an indel are detected (Altschul et al. 1990; Smith &amp; Waterman 1981).
    /// </summary>
    /// <param name="ReferenceIndex">Zero-based index of the reference sequence (in the supplied order) the hit was found in.</param>
    /// <param name="Start">Zero-based start position of the aligned region within the reference.</param>
    /// <param name="End">Zero-based end position (inclusive) of the aligned region within the reference.</param>
    /// <param name="Identity">
    /// Fraction of identical aligned columns over the probe length
    /// (identical columns / probe length), in [0, 1]. Computed exactly as in the library's
    /// reused ANI identity convention (Goris et al. 2007; pyani <c>ani_pid = ani_alnids / qlen</c>).
    /// </param>
    /// <param name="Coverage">Fraction of the probe covered by ungapped aligned columns (ungapped columns / probe length), in [0, 1].</param>
    /// <param name="HasGaps">True when the alignment contains at least one gap (insertion or deletion) — i.e. the hit is reachable only with an indel.</param>
    /// <param name="AlignedProbe">The probe side of the local alignment (may contain '-').</param>
    /// <param name="AlignedReference">The reference side of the local alignment (may contain '-').</param>
    public readonly record struct GappedProbeHit(
        int ReferenceIndex,
        int Start,
        int End,
        double Identity,
        double Coverage,
        bool HasGaps,
        string AlignedProbe,
        string AlignedReference);

    /// <summary>
    /// Result of the opt-in gapped (Smith-Waterman) off-target scan, separating the intended
    /// on-target match from genuine off-target hits.
    /// </summary>
    /// <param name="OnTargetHits">
    /// Perfect, ungapped, full-coverage exact matches (identity = 1.0, coverage = 1.0, no gaps).
    /// The first such hit is the probe's intended on-target binding site; any additional perfect
    /// exact matches are reported here too (a probe matching several identical sites is itself a
    /// specificity concern) but the first one is excluded from <see cref="OffTargetHits"/>.
    /// </param>
    /// <param name="OffTargetHits">
    /// Genuine off-target hits: every hit at or above the identity threshold that is NOT the single
    /// intended on-target site. This includes imperfect (mismatched) hits, indel-containing hits, and
    /// any extra perfect repeats. This is the corrected count that no longer pools the on-target match
    /// with off-targets (cf. <see cref="ProbeValidation.OffTargetHits"/>).
    /// </param>
    /// <param name="MinIdentity">The identity threshold used to call a hit (see <see cref="ScanOffTargetsGapped"/>).</param>
    public readonly record struct GappedSpecificityResult(
        IReadOnlyList<GappedProbeHit> OnTargetHits,
        IReadOnlyList<GappedProbeHit> OffTargetHits,
        double MinIdentity)
    {
        /// <summary>Number of genuine off-target hits (excludes the intended on-target match).</summary>
        public int OffTargetCount => OffTargetHits.Count;

        /// <summary>True when exactly one on-target site and no off-target hits were found.</summary>
        public bool IsSpecific => OnTargetHits.Count == 1 && OffTargetHits.Count == 0;
    }

    /// <summary>
    /// Karlin–Altschul statistics of an off-target alignment hit: the raw alignment score's
    /// statistical significance expressed as a bit score and an expectation (E) value.
    /// </summary>
    /// <remarks>
    /// Karlin &amp; Altschul (1990, PNAS 87:2264); Altschul et al. (1990, J Mol Biol 215:403).
    /// </remarks>
    /// <param name="RawScore">The raw alignment score S (sum of substitution/gap scores).</param>
    /// <param name="Lambda">
    /// The Karlin–Altschul scale parameter λ — the unique positive root of
    /// Σ_{i,j} p_i p_j e^{λ s_ij} = 1 for the scoring scheme and base frequencies.
    /// </param>
    /// <param name="K">The Karlin–Altschul search-space scale parameter K (supplied by the caller).</param>
    /// <param name="BitScore">The normalized bit score S' = (λS − ln K) / ln 2.</param>
    /// <param name="EValue">
    /// The expected number of distinct alignments scoring ≥ S by chance: E = K·m·n·e^{−λS} = m·n·2^{−S'}.
    /// </param>
    /// <param name="QueryLength">The query (probe) length m used in the search space.</param>
    /// <param name="DatabaseLength">The database (reference) length n used in the search space.</param>
    public readonly record struct KarlinAltschulStatistics(
        double RawScore,
        double Lambda,
        double K,
        double BitScore,
        double EValue,
        int QueryLength,
        long DatabaseLength);

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

    /// <summary>
    /// Result of evaluating a candidate hydrolysis (TaqMan) probe against the
    /// Applied Biosystems / Thermo Fisher TaqMan probe-design guidelines.
    /// Each boolean records whether one published rule is satisfied; <see cref="PassesAll"/>
    /// is the conjunction of every rule.
    /// </summary>
    /// <remarks>
    /// Sources (retrieved 2026-06-24):
    /// PREMIER Biosoft, "TaqMan probe design tips" (http://www.premierbiosoft.com/tech_notes/TaqMan.html);
    /// Thermo Fisher / Applied Biosystems "Designing a TaqMan Gene Expression Assay" and
    /// "TaqMan MGB Probe and Primer Sets" application notes.
    /// </remarks>
    public readonly record struct TaqManProbeEvaluation(
        string Sequence,
        bool NoGuanineAt5Prime,
        bool MoreCytosineThanGuanine,
        bool NoRunOfFourOrMoreG,
        bool GcContentInRange,
        bool LengthInRange,
        bool ProbeTmAbovePrimer,
        double Tm,
        double GcContent,
        int CytosineCount,
        int GuanineCount,
        bool PassesAll,
        IReadOnlyList<string> Violations);

    #endregion

    #region TaqMan (hydrolysis-probe) design — opt-in

    // --- Published TaqMan probe-design thresholds (Applied Biosystems / Thermo Fisher; PREMIER Biosoft) ---

    // Probe length range. PREMIER Biosoft: "TaqMan probes consist of a 18-22 bp oligonucleotide probe".
    // (IDT / Thermo state 18-30; we use the tighter ABI/PREMIER 18-22 default, configurable.)
    private const int TaqManMinLength = 18;
    private const int TaqManMaxLength = 22;

    // G+C content range. PREMIER Biosoft: "The G+C content should ideally be 30-80%".
    private const double TaqManMinGc = 0.30;
    private const double TaqManMaxGc = 0.80;

    // No runs of identical nucleotides, "especially four or more consecutive Gs" (PREMIER Biosoft / ABI).
    private const int TaqManMaxGuanineRun = 4;

    // Probe Tm should be ~10 °C higher than the primer Tm so the probe binds before Taq extends.
    // PREMIER Biosoft / Thermo Fisher: "TaqMan probe Tm should be 10 °C higher than the Primer Tm".
    private const double TaqManProbeTmDeltaAbovePrimer = 10.0;

    /// <summary>
    /// Evaluates a single candidate probe sequence against the published TaqMan
    /// (5'-nuclease hydrolysis probe) design guidelines. This is an opt-in chemistry-specific
    /// check; the generic <see cref="DesignProbes(string, ProbeParameters?, int)"/> designer is unchanged.
    /// </summary>
    /// <param name="probeSequence">The candidate probe sequence (5'→3').</param>
    /// <param name="primerTm">
    /// Melting temperature (°C) of the amplification primers. The probe Tm must be at least
    /// <c>primerTm + 10 °C</c> (Applied Biosystems / Thermo Fisher). Pass <see langword="null"/>
    /// to skip the probe-Tm-vs-primer gate (it is then reported as satisfied).
    /// </param>
    /// <param name="minLength">Minimum probe length (default 18, per ABI/PREMIER Biosoft).</param>
    /// <param name="maxLength">Maximum probe length (default 22, per ABI/PREMIER Biosoft).</param>
    /// <returns>A <see cref="TaqManProbeEvaluation"/> recording each rule outcome.</returns>
    public static TaqManProbeEvaluation EvaluateTaqManProbe(
        string probeSequence,
        double? primerTm = null,
        int minLength = TaqManMinLength,
        int maxLength = TaqManMaxLength)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);

        string seq = probeSequence.ToUpperInvariant();
        var violations = new List<string>();

        // Rule 1: no G at the 5' end (a 5' G adjacent to the reporter dye quenches
        // reporter fluorescence even after cleavage).
        bool noGuanineAt5Prime = seq.Length > 0 && seq[0] != 'G';
        if (!noGuanineAt5Prime)
            violations.Add("5' end is a guanine (quenches the reporter dye even after cleavage)");

        // Rule 2: more Cs than Gs in the probe sequence.
        int cCount = seq.Count(c => c == 'C');
        int gCount = seq.Count(c => c == 'G');
        bool moreCThanG = cCount > gCount;
        if (!moreCThanG)
            violations.Add($"not more C than G (C={cCount}, G={gCount})");

        // Rule 3: no run of four or more consecutive Gs.
        int maxGRun = GetMaxGuanineRunLength(seq);
        bool noRunOf4G = maxGRun < TaqManMaxGuanineRun;
        if (!noRunOf4G)
            violations.Add($"run of {maxGRun} consecutive Gs (>= {TaqManMaxGuanineRun})");

        // Rule 4: G+C content within 30-80%.
        double gc = CalculateGcContent(seq);
        bool gcInRange = gc >= TaqManMinGc && gc <= TaqManMaxGc;
        if (!gcInRange)
            violations.Add($"G+C content {gc:P0} outside {TaqManMinGc:P0}-{TaqManMaxGc:P0}");

        // Rule 5: probe length within range.
        bool lengthInRange = seq.Length >= minLength && seq.Length <= maxLength;
        if (!lengthInRange)
            violations.Add($"length {seq.Length} outside {minLength}-{maxLength} nt");

        // Rule 6: probe Tm at least ~10 °C above the primer Tm (when a primer Tm is supplied).
        double tm = CalculateTm(seq);
        bool probeTmAbovePrimer = !primerTm.HasValue || tm >= primerTm.Value + TaqManProbeTmDeltaAbovePrimer;
        if (!probeTmAbovePrimer)
            violations.Add(
                $"probe Tm {tm:F1} °C is not >= primer Tm {primerTm!.Value:F1} + {TaqManProbeTmDeltaAbovePrimer:F0} °C");

        bool passesAll = noGuanineAt5Prime && moreCThanG && noRunOf4G
            && gcInRange && lengthInRange && probeTmAbovePrimer;

        return new TaqManProbeEvaluation(
            seq,
            noGuanineAt5Prime,
            moreCThanG,
            noRunOf4G,
            gcInRange,
            lengthInRange,
            probeTmAbovePrimer,
            tm,
            gc,
            cCount,
            gCount,
            passesAll,
            violations);
    }

    /// <summary>
    /// Chooses the better TaqMan probe strand. Per Applied Biosystems / Thermo Fisher,
    /// the probe should be designed on the strand with <b>more Cs than Gs</b>; if a guanine
    /// occurs at the 5' end of one strand, the complement (antisense) strand should be used.
    /// Returns the sense (given) strand or its reverse complement, whichever better satisfies the rules.
    /// </summary>
    /// <param name="senseStrand">The candidate probe sequence on the sense strand (5'→3').</param>
    /// <param name="primerTm">Optional primer Tm for the probe-Tm gate (see <see cref="EvaluateTaqManProbe"/>).</param>
    /// <returns>
    /// A tuple of the chosen probe sequence (5'→3'), whether it is the reverse-complement
    /// (antisense) strand, and the evaluation of the chosen strand.
    /// </returns>
    public static (string Probe, bool IsReverseComplement, TaqManProbeEvaluation Evaluation) SelectTaqManStrand(
        string senseStrand,
        double? primerTm = null)
    {
        ArgumentNullException.ThrowIfNull(senseStrand);

        string sense = senseStrand.ToUpperInvariant();
        string antisense = DnaSequence.GetReverseComplementString(sense);

        var senseEval = EvaluateTaqManProbe(sense, primerTm);
        var antisenseEval = EvaluateTaqManProbe(antisense, primerTm);

        // Prefer a strand that passes all rules; otherwise prefer the one satisfying the two
        // hard reporter-dye rules (no 5'-G, then more C than G); finally fall back to the sense strand.
        if (senseEval.PassesAll && !antisenseEval.PassesAll)
            return (sense, false, senseEval);
        if (antisenseEval.PassesAll && !senseEval.PassesAll)
            return (antisense, true, antisenseEval);

        int senseRank = RankTaqManStrand(senseEval);
        int antisenseRank = RankTaqManStrand(antisenseEval);

        return antisenseRank > senseRank
            ? (antisense, true, antisenseEval)
            : (sense, false, senseEval);
    }

    // Ranks a strand: the no-5'-G rule then the more-C-than-G rule are the chemistry-critical
    // reporter-dye constraints, weighted above the remaining quality rules.
    private static int RankTaqManStrand(TaqManProbeEvaluation e)
    {
        int rank = 0;
        if (e.NoGuanineAt5Prime) rank += 4;
        if (e.MoreCytosineThanGuanine) rank += 2;
        if (e.NoRunOfFourOrMoreG) rank += 1;
        if (e.GcContentInRange) rank += 1;
        if (e.ProbeTmAbovePrimer) rank += 1;
        return rank;
    }

    private static int GetMaxGuanineRunLength(string sequence)
    {
        int maxRun = 0;
        int currentRun = 0;

        foreach (char c in sequence)
        {
            if (c == 'G')
            {
                currentRun++;
                maxRun = Math.Max(maxRun, currentRun);
            }
            else
            {
                currentRun = 0;
            }
        }

        return maxRun;
    }

    #endregion

    #region MGB (minor-groove binder) design rules — opt-in (qualitative; quantitative ΔTm is a residual)

    // --- Citable 3'-MGB probe-design rules (Kutyavin et al. 2000, Nucleic Acids Res 28(2):655-661) ---
    //
    // Kutyavin et al. (2000) established that conjugating a minor-groove binder (MGB) to the 3' end
    // of a DNA probe greatly stabilises the duplex, so MGB probes are designed SHORTER than
    // unmodified probes: "for MGB probes this length variation is narrowed to a range of 12-20mers"
    // (a 12mer MGB has ~the same Tm as a 27mer unmodified probe). The MGB is attached at the 3' end
    // ("3'-MGB-ODNs are easier to prepare … MGB-modified solid supports and automated DNA synthesis
    // can be used"). The QUANTITATIVE MGB ΔTm is empirical with no published closed-form model
    // (the stabilisation "varies by sequence"; A+T-rich MGB sites gain more than G+C-rich ones), so
    // only these qualitative DESIGN rules are implemented here; the quantitative MGB ΔTm is left as
    // an honest residual (see docs/Validation/LIMITATIONS.md).
    // Source (retrieved 2026-06-24): Kutyavin IV et al. (2000) Nucleic Acids Res 28(2):655-661,
    //   https://doi.org/10.1093/nar/28.2.655.

    // Kutyavin (2000): MGB-probe length range narrowed to 12-20mers.
    private const int MgbMinLength = 12;
    private const int MgbMaxLength = 20;

    /// <summary>
    /// Result of checking a candidate probe against the citable 3'-MGB (minor-groove binder)
    /// design rules of Kutyavin et al. (2000). These are <b>qualitative</b> design-rule checks
    /// only; the quantitative MGB ΔTm is empirical (no published formula) and is not computed.
    /// </summary>
    /// <param name="Sequence">The (upper-cased) probe sequence checked.</param>
    /// <param name="LengthInMgbRange">True when the probe length is within the MGB 12–20mer window.</param>
    /// <param name="Length">The probe length in nucleotides.</param>
    /// <param name="MgbAttachmentEnd">The recommended MGB attachment end (always 3', per Kutyavin 2000).</param>
    /// <param name="Guidance">Human-readable design guidance / any rule violations.</param>
    public readonly record struct MgbProbeDesign(
        string Sequence,
        bool LengthInMgbRange,
        int Length,
        string MgbAttachmentEnd,
        IReadOnlyList<string> Guidance);

    /// <summary>
    /// Evaluates a candidate probe against the citable 3'-MGB (minor-groove binder) probe-design
    /// rules of Kutyavin et al. (2000): the MGB is attached at the <b>3' end</b>, and MGB probes are
    /// designed <b>shorter</b> (12–20mer) than unmodified probes. This is an opt-in, qualitative
    /// design-rule check; the generic <see cref="DesignProbes(string, ProbeParameters?, int)"/>
    /// designer and all defaults are unchanged. The <b>quantitative</b> MGB ΔTm is empirical (no
    /// published closed-form model in Kutyavin 2000) and is deliberately NOT computed.
    /// </summary>
    /// <param name="probeSequence">The candidate probe sequence (5'→3').</param>
    /// <returns>An <see cref="MgbProbeDesign"/> recording the length-window outcome and 3'-MGB guidance.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="probeSequence"/> is null.</exception>
    public static MgbProbeDesign EvaluateMgbProbeDesign(string probeSequence)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);

        string seq = probeSequence.ToUpperInvariant();
        var guidance = new List<string>();

        bool lengthInRange = seq.Length >= MgbMinLength && seq.Length <= MgbMaxLength;
        if (!lengthInRange)
            guidance.Add(
                $"length {seq.Length} outside the MGB {MgbMinLength}-{MgbMaxLength}mer window (Kutyavin 2000)");

        // 3'-MGB placement guidance (always applies).
        guidance.Add("attach the minor-groove binder at the 3' end (Kutyavin 2000)");

        return new MgbProbeDesign(
            seq,
            lengthInRange,
            seq.Length,
            MgbAttachmentEnd: "3'",
            guidance);
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

        // Use optimized evaluation with prefix sums for O(1) GC lookup
        var candidates = DesignProbesOptimized(targetSequence, param, maxProbes);

        foreach (var probe in candidates)
        {
            yield return probe;
        }
    }

    /// <summary>
    /// Designs probes with genome-wide specificity check using suffix tree.
    /// O(n × m) for probe generation + O(m) per specificity check.
    /// </summary>
    /// <param name="targetSequence">Target sequence to design probes for.</param>
    /// <param name="genomeIndex">Pre-built suffix tree index for the genome (enables O(m) specificity lookup).</param>
    /// <param name="parameters">Probe design parameters.</param>
    /// <param name="maxProbes">Maximum number of probes to return.</param>
    /// <param name="requireUnique">If true, only return probes unique in the genome.</param>
    public static IEnumerable<Probe> DesignProbes(
        string targetSequence,
        global::SuffixTree.ISuffixTree genomeIndex,
        ProbeParameters? parameters = null,
        int maxProbes = 10,
        bool requireUnique = true)
    {
        var param = parameters ?? Defaults.Microarray;

        if (string.IsNullOrEmpty(targetSequence) || targetSequence.Length < param.MinLength)
            yield break;

        targetSequence = targetSequence.ToUpperInvariant();

        // Get candidates using optimized method
        var candidates = DesignProbesOptimized(targetSequence, param, maxProbes * 5); // Get more candidates for filtering

        int returned = 0;
        foreach (var probe in candidates)
        {
            if (returned >= maxProbes)
                yield break;

            // Fast O(m) specificity check using suffix tree
            double specificity = CheckSpecificity(probe.Sequence, genomeIndex);

            if (requireUnique && specificity < 1.0)
                continue; // Skip non-unique probes

            // Boost score based on specificity
            var adjustedProbe = probe with
            {
                Score = probe.Score * specificity
            };

            returned++;
            yield return adjustedProbe;
        }
    }

    /// <summary>
    /// Optimized probe design using prefix sums for O(1) GC content calculation.
    /// Total complexity: O(n × m) where n = sequence length, m = length range.
    /// </summary>
    private static List<Probe> DesignProbesOptimized(
        string targetSequence,
        ProbeParameters param,
        int maxProbes)
    {
        int n = targetSequence.Length;

        // Precompute GC prefix sums for O(1) GC content queries
        // gcPrefixSum[i] = count of G/C in sequence[0..i-1]
        int[] gcPrefixSum = new int[n + 1];
        for (int i = 0; i < n; i++)
        {
            char c = targetSequence[i];
            gcPrefixSum[i + 1] = gcPrefixSum[i] + (c == 'G' || c == 'C' ? 1 : 0);
        }

        var candidates = new List<(Probe Probe, double Score)>();

        // Scan for candidate probes
        for (int length = param.MinLength; length <= param.MaxLength && length <= n; length++)
        {
            for (int start = 0; start <= n - length; start++)
            {
                // O(1) GC content using prefix sums
                int gcCount = gcPrefixSum[start + length] - gcPrefixSum[start];
                double gc = (double)gcCount / length;

                // Early rejection based on GC - saves expensive substring operations
                if (gc < param.MinGc - 0.1 || gc > param.MaxGc + 0.1)
                    continue;

                string probeSeq = targetSequence.Substring(start, length);
                var probe = EvaluateProbeWithGc(probeSeq, start, param, gc);

                if (probe.HasValue)
                {
                    candidates.Add((probe.Value, probe.Value.Score));
                }
            }
        }

        // Return top probes sorted by score
        return candidates
            .OrderByDescending(c => c.Score)
            .Take(maxProbes)
            .Select(c => c.Probe)
            .ToList();
    }

    /// <summary>
    /// Evaluates probe with pre-calculated GC content (avoids redundant calculation).
    /// </summary>
    private static Probe? EvaluateProbeWithGc(string sequence, int start, ProbeParameters param, double gc)
    {
        var warnings = new List<string>();
        double score = 1.0;

        // GC already calculated - just check bounds
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
    /// Evaluates a potential probe sequence.
    /// </summary>
    private static Probe? EvaluateProbe(string sequence, int start, ProbeParameters param)
    {
        double gc = CalculateGcContent(sequence);
        return EvaluateProbeWithGc(sequence, start, param, gc);
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
    /// <param name="probeSequence">Probe sequence to validate.</param>
    /// <param name="referenceSequences">Reference sequences to search for off-target hits.</param>
    /// <param name="maxMismatches">Maximum allowed mismatches for approximate matching.
    /// Default: 3, based on CRISPR/Cas9 off-target tolerance of 3-5 bp mismatches
    /// per 20nt guide (Hsu et al. 2013, Fu et al. 2013).</param>
    /// <param name="selfComplementarityThreshold">Threshold above which self-complementarity
    /// generates a warning. Default: 0.3 (Microarray default). For random DNA the expected
    /// self-complementarity is ~0.25; values above this threshold indicate elevated
    /// palindromic character that may cause probe secondary structure.</param>
    public static ProbeValidation ValidateProbe(
        string probeSequence,
        IEnumerable<string> referenceSequences,
        int maxMismatches = 3,
        double selfComplementarityThreshold = 0.3)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(referenceSequences);

        probeSequence = probeSequence.ToUpperInvariant();
        var issues = new List<string>();

        // Empty probe is a degenerate input — cannot hybridize specifically
        if (probeSequence.Length == 0)
        {
            return new ProbeValidation(
                IsValid: false,
                SpecificityScore: 0.0,
                OffTargetHits: 0,
                SelfComplementarity: 0.0,
                HasSecondaryStructure: false,
                Issues: new List<string> { "Empty probe sequence" });
        }

        int offTargetHits = 0;

        // Check off-target hits via approximate matching
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
        if (selfComp > selfComplementarityThreshold)
        {
            issues.Add($"Self-complementarity: {selfComp:P0}");
        }

        // Check secondary structure
        bool hasStructure = HasSecondaryStructurePotential(probeSequence);
        if (hasStructure)
        {
            issues.Add("Potential secondary structure formation");
        }

        // Calculate specificity score:
        // 0 hits → 0.0 (probe doesn't hybridize to target — useless)
        // 1 hit  → 1.0 (unique match — ideal specificity)
        // N hits → 1.0/N (specificity decreases with cross-hybridization)
        double specificity;
        if (offTargetHits == 0)
            specificity = 0.0;
        else if (offTargetHits == 1)
            specificity = 1.0;
        else
            specificity = 1.0 / offTargetHits;

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

    // --- Gapped off-target scan thresholds (sourced) ---

    // Default minimum identity (identical aligned columns / probe length) to call an off-target.
    // Kane et al. (2000, Nucleic Acids Res. 28(22):4552-4557): "for a given oligonucleotide probe
    // any 'non-target' transcripts (cDNAs) >75% similar over the [...] target may show
    // cross-hybridization." 0.75 is therefore the empirically-grounded similarity threshold above
    // which cross-hybridization (an off-target) becomes a concern; callers may override it.
    private const double DefaultOffTargetMinIdentity = 0.75;

    // Extra reference length scanned per window beyond the probe length, to allow a few indels in
    // the local alignment (a gap shifts the reference frame). Two extra bases lets the Smith-Waterman
    // local alignment absorb short insertions/deletions while keeping each window O(probeLen) wide.
    private const int GappedScanGapAllowance = 2;

    // Gap character emitted by SequenceAligner in its aligned-output strings.
    private const char AlignmentGapChar = '-';

    /// <summary>
    /// Opt-in <b>gapped</b> off-target scan using the library's validated Smith-Waterman local
    /// aligner (<see cref="SequenceAligner.LocalAlign(string, string, ScoringMatrix?)"/>). Unlike the
    /// ungapped Hamming-distance scan in <see cref="ValidateProbe"/> (which only tolerates
    /// substitutions in a fixed-length window), this finds off-target sites reachable through
    /// insertions or deletions — the "BLAST-grade" improvement (gapped local alignment handles
    /// indels the ungapped scan misses; Altschul et al. 1990; Smith &amp; Waterman 1981).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The default ungapped <see cref="ValidateProbe"/> behaviour is unchanged; this is a separate,
    /// additive entry point. It also corrects the on/off-target pooling of
    /// <see cref="ProbeValidation.OffTargetHits"/>: the single intended on-target site (the perfect,
    /// ungapped, full-coverage exact match) is reported separately and is NOT counted as an off-target.
    /// </para>
    /// <para>
    /// This is an exhaustive sliding Smith-Waterman scan (O(g · n · m) over reference length g and
    /// probe length n·m), not a seeded BLAST index over a whole genome; for genome-scale search a
    /// k-mer/seed index would be required.
    /// </para>
    /// </remarks>
    /// <param name="probeSequence">Probe sequence to scan (5'→3'). Null throws; empty yields no hits.</param>
    /// <param name="referenceSequences">Reference sequences to scan for hits. Null throws.</param>
    /// <param name="minIdentity">
    /// Minimum alignment identity (identical aligned columns / probe length) to call a hit.
    /// Default 0.75 per Kane et al. (2000): non-targets &gt;75% similar over the probe may cross-hybridize.
    /// </param>
    /// <param name="scoring">
    /// Scoring matrix for the local alignment. Defaults to <see cref="SequenceAligner.BlastDna"/>
    /// (+2/-3, gap -2) — the same BLAST-style DNA scoring reused for the library's gapped ANI alignment.
    /// </param>
    /// <returns>A <see cref="GappedSpecificityResult"/> separating on-target from off-target hits.</returns>
    public static GappedSpecificityResult ScanOffTargetsGapped(
        string probeSequence,
        IEnumerable<string> referenceSequences,
        double minIdentity = DefaultOffTargetMinIdentity,
        ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(referenceSequences);

        var onTarget = new List<GappedProbeHit>();
        var offTarget = new List<GappedProbeHit>();

        string probe = probeSequence.ToUpperInvariant();
        if (probe.Length == 0)
            return new GappedSpecificityResult(onTarget, offTarget, minIdentity);

        var matrix = scoring ?? SequenceAligner.BlastDna;
        bool onTargetClaimed = false;

        int refIndex = 0;
        foreach (var reference in referenceSequences)
        {
            string text = (reference ?? string.Empty).ToUpperInvariant();
            foreach (var hit in ScanReferenceGapped(probe, text, refIndex, minIdentity, matrix))
            {
                // The intended on-target is the (first) perfect ungapped full-coverage exact match.
                bool isPerfectExact = !hit.HasGaps
                    && hit.Identity >= 1.0
                    && hit.Coverage >= 1.0;

                if (isPerfectExact && !onTargetClaimed)
                {
                    onTargetClaimed = true;
                    onTarget.Add(hit);
                }
                else if (isPerfectExact)
                {
                    // An additional perfect repeat: report it as an on-target-class match for
                    // visibility, but it is still a genuine extra binding site → counts as off-target.
                    onTarget.Add(hit);
                    offTarget.Add(hit);
                }
                else
                {
                    offTarget.Add(hit);
                }
            }

            refIndex++;
        }

        return new GappedSpecificityResult(onTarget, offTarget, minIdentity);
    }

    /// <summary>
    /// Slides a Smith-Waterman local alignment of <paramref name="probe"/> across
    /// <paramref name="text"/> and returns one best non-overlapping hit per distinct site whose
    /// identity is at least <paramref name="minIdentity"/>.
    /// </summary>
    private static List<GappedProbeHit> ScanReferenceGapped(
        string probe, string text, int refIndex, double minIdentity, ScoringMatrix matrix)
    {
        int probeLen = probe.Length;
        int windowLen = probeLen + GappedScanGapAllowance;
        var raw = new List<GappedProbeHit>();

        for (int start = 0; start < text.Length; start++)
        {
            int span = Math.Min(windowLen, text.Length - start);
            if (span <= 0)
                break;

            string window = text.Substring(start, span);
            AlignmentResult aln = SequenceAligner.LocalAlign(probe, window, matrix);

            // SequenceAligner.LocalAlign aligns (sequence1 = probe, sequence2 = window):
            // AlignedSequence1 is the probe side, AlignedSequence2 the reference side.
            string ap = aln.AlignedSequence1;
            string ar = aln.AlignedSequence2;
            if (ap.Length == 0)
                continue;

            int identical = 0;
            int ungapped = 0;
            bool hasGap = false;
            for (int k = 0; k < ap.Length; k++)
            {
                char c1 = ap[k];
                char c2 = ar[k];
                if (c1 == AlignmentGapChar || c2 == AlignmentGapChar)
                {
                    hasGap = true;
                    continue;
                }

                ungapped++;
                if (c1 == c2)
                    identical++;
            }

            double identity = (double)identical / probeLen;
            if (identity < minIdentity)
                continue;

            double coverage = (double)ungapped / probeLen;
            int absStart = start + aln.StartPosition2;
            int absEnd = start + aln.EndPosition2;

            raw.Add(new GappedProbeHit(
                refIndex, absStart, absEnd, identity, coverage, hasGap, ap, ar));
        }

        // Greedy best-per-site selection: take highest-identity (then highest-coverage, then
        // leftmost) hits first, accepting a hit only if its reference span does not overlap an
        // already-accepted one. Overlapping windows that re-detect the same site collapse to one hit.
        raw.Sort((a, b) =>
        {
            int byIdentity = b.Identity.CompareTo(a.Identity);
            if (byIdentity != 0) return byIdentity;
            int byCoverage = b.Coverage.CompareTo(a.Coverage);
            if (byCoverage != 0) return byCoverage;
            return a.Start.CompareTo(b.Start);
        });

        var accepted = new List<GappedProbeHit>();
        foreach (var hit in raw)
        {
            bool overlaps = accepted.Any(a => !(hit.End < a.Start || hit.Start > a.End));
            if (!overlaps)
                accepted.Add(hit);
        }

        accepted.Sort((a, b) => a.Start.CompareTo(b.Start));
        return accepted;
    }

    // --- Karlin–Altschul statistics (opt-in) ---
    //
    // Karlin & Altschul (1990, PNAS 87:2264) / Altschul et al. (1990, J Mol Biol 215:403):
    //   E = K·m·n·e^{−λS}                       (expected HSPs with score ≥ S by chance)
    //   S' = (λS − ln K) / ln 2                  (normalized "bit" score)
    //   E = m·n·2^{−S'}                          (E in terms of the bit score)
    //   λ is the unique positive root of  Σ_{i,j} p_i p_j e^{λ s_ij} = 1.
    // Verbatim formulas retrieved 2026-06-24 from the NCBI BLAST course
    //   "The Statistics of Sequence Similarity Scores" (Altschul),
    //   https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html, and from Durand,
    //   "BLAST (Karlin–Altschul) Statistics" (CMU 03-711, citing Karlin & Altschul 1990 and
    //   Altschul et al. 1990), http://www.cs.cmu.edu/~durand/03-711/2011/Lectures/Blast-informationContent-2011.pdf.
    // The theory requires a scoring scheme whose expected per-pair score is negative and that
    // has at least one positive score (Altschul, ibid.); otherwise λ is undefined.

    // Uniform nucleotide background frequency p_i = 0.25 for the four bases A,C,G,T
    // (the standard assumption when computing the +1/−3 nucleotide λ ≈ 1.374, NCBI blastn).
    private const double UniformBaseFrequency = 0.25;

    // Bisection bounds/iterations for solving Σ p_i p_j e^{λ s_ij} = 1. The function is strictly
    // increasing in λ once it crosses 1 (the expected score is negative, so it starts < 1 and a
    // positive score makes it diverge), so a simple bisection converges to the unique positive root.
    private const double LambdaSearchUpperBound = 100.0;
    private const int LambdaBisectionIterations = 200;

    /// <summary>
    /// Computes the Karlin–Altschul scale parameter λ for a match/mismatch nucleotide scoring
    /// scheme under uniform (0.25) base frequencies, by solving the defining equation
    /// Σ_{i,j} p_i p_j e^{λ s_ij} = 1 numerically (bisection on the unique positive root).
    /// </summary>
    /// <remarks>
    /// With four equiprobable bases and a simple match/mismatch matrix, the 16 ordered pairs split
    /// into 4 matches (probability 4·0.25² = 0.25) and 12 mismatches (probability 0.75), so the
    /// equation reduces to 0.25·e^{λ·match} + 0.75·e^{λ·mismatch} = 1. For the BLAST +1/−3 scheme
    /// this yields λ ≈ 1.374 (NCBI blastn). Karlin &amp; Altschul (1990); Altschul et al. (1990).
    /// </remarks>
    /// <param name="match">Match score (must be &gt; 0 — the required positive score).</param>
    /// <param name="mismatch">Mismatch score (must be &lt; 0).</param>
    /// <param name="baseFrequency">
    /// Per-base background frequency (default 0.25, uniform). The four bases are assumed equiprobable
    /// at this value; the four base frequencies must sum to 1, i.e. <paramref name="baseFrequency"/> = 0.25.
    /// </param>
    /// <returns>The positive λ solving the Karlin–Altschul equation.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when the scheme cannot define λ: the match score is not positive, the mismatch score is
    /// not negative, or the expected per-pair score is not negative (the theory's preconditions).
    /// </exception>
    public static double ComputeLambdaNucleotide(
        int match,
        int mismatch,
        double baseFrequency = UniformBaseFrequency)
    {
        // Karlin–Altschul preconditions: at least one positive score, and negative expected score.
        if (match <= 0)
            throw new ArgumentOutOfRangeException(nameof(match),
                "Karlin–Altschul λ is undefined: the scoring scheme must have at least one positive score.");
        if (mismatch >= 0)
            throw new ArgumentOutOfRangeException(nameof(mismatch),
                "Karlin–Altschul λ is undefined: the mismatch score must be negative.");

        // p(match) = 4 · p² (the four identical ordered pairs); p(mismatch) = 1 − p(match).
        double pMatch = 4.0 * baseFrequency * baseFrequency;
        double pMismatch = 1.0 - pMatch;

        // Expected per-pair score must be negative for the theory to hold.
        double expectedScore = pMatch * match + pMismatch * mismatch;
        if (expectedScore >= 0)
            throw new ArgumentOutOfRangeException(nameof(mismatch),
                "Karlin–Altschul λ is undefined: the expected per-pair score must be negative.");

        // f(λ) = p(match)·e^{λ·match} + p(mismatch)·e^{λ·mismatch} − 1.
        // f(0) = 0; f'(0) = expectedScore < 0 so f dips below 0 for small λ>0, then a positive
        // match score drives e^{λ·match} → ∞, so f crosses 0 exactly once at the positive root.
        static double F(double lambda, double pM, int m, double pMm, int mm)
            => pM * Math.Exp(lambda * m) + pMm * Math.Exp(lambda * mm) - 1.0;

        double lo = 0.0;
        double hi = LambdaSearchUpperBound;
        for (int i = 0; i < LambdaBisectionIterations; i++)
        {
            double mid = 0.5 * (lo + hi);
            if (F(mid, pMatch, match, pMismatch, mismatch) > 0.0)
                hi = mid;
            else
                lo = mid;
        }

        return 0.5 * (lo + hi);
    }

    /// <summary>
    /// Computes the Karlin–Altschul bit score and E-value for an off-target alignment hit's raw
    /// score, given the search-space dimensions and the scoring scheme.
    /// </summary>
    /// <remarks>
    /// <para>
    /// E = K·m·n·e^{−λS}, S' = (λS − ln K) / ln 2, E = m·n·2^{−S'} (Karlin &amp; Altschul 1990;
    /// Altschul et al. 1990). λ is computed from <paramref name="scoring"/> by
    /// <see cref="ComputeLambdaNucleotide"/>; K is supplied by the caller (the closed form requires
    /// the score-lattice machinery of Karlin–Altschul; for the BLAST +1/−3 nucleotide scheme the
    /// published value is K ≈ 0.711, NCBI blastn).
    /// </para>
    /// <para>This is additive and opt-in; <see cref="ScanOffTargetsGapped"/> and its defaults are unchanged.</para>
    /// </remarks>
    /// <param name="rawScore">The raw alignment score S of the hit.</param>
    /// <param name="queryLength">Query (probe) length m (&gt; 0).</param>
    /// <param name="databaseLength">Database (reference) length n (&gt; 0).</param>
    /// <param name="scoring">
    /// The scoring scheme. Its <see cref="ScoringMatrix.Match"/>/<see cref="ScoringMatrix.Mismatch"/>
    /// determine λ. Defaults to <see cref="SequenceAligner.BlastDna"/> (+2/−3). Pass a +1/−3 matrix to
    /// reproduce the published λ ≈ 1.374.
    /// </param>
    /// <param name="k">The Karlin–Altschul K parameter (default 0.711, the published nucleotide value).</param>
    /// <param name="baseFrequency">Per-base background frequency for λ (default 0.25, uniform).</param>
    /// <returns>The <see cref="KarlinAltschulStatistics"/> for the hit.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="scoring"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown for non-positive lengths or K, or a scheme for which λ is undefined.</exception>
    public static KarlinAltschulStatistics ComputeKarlinAltschul(
        double rawScore,
        int queryLength,
        long databaseLength,
        ScoringMatrix? scoring = null,
        double k = DefaultNucleotideK,
        double baseFrequency = UniformBaseFrequency)
    {
        var matrix = scoring ?? SequenceAligner.BlastDna;
        ArgumentNullException.ThrowIfNull(matrix);
        if (queryLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(queryLength), "Query length m must be positive.");
        if (databaseLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(databaseLength), "Database length n must be positive.");
        if (k <= 0)
            throw new ArgumentOutOfRangeException(nameof(k), "K must be positive.");

        double lambda = ComputeLambdaNucleotide(matrix.Match, matrix.Mismatch, baseFrequency);

        // S' = (λS − ln K) / ln 2  (Altschul et al. 1990).
        double bitScore = (lambda * rawScore - Math.Log(k)) / Math.Log(2.0);

        // E = K·m·n·e^{−λS}  (Karlin & Altschul 1990).
        double eValue = k * queryLength * databaseLength * Math.Exp(-lambda * rawScore);

        return new KarlinAltschulStatistics(
            rawScore, lambda, k, bitScore, eValue, queryLength, databaseLength);
    }

    // Published Karlin–Altschul K for the BLAST +1/−3 nucleotide scheme (NCBI blastn reports
    // Lambda ≈ 1.37, K ≈ 0.711 for match=1/mismatch=−3). K's full closed form needs the
    // score-probability lattice/geometric-spacing machinery of Karlin & Altschul (1990); it is
    // therefore exposed as a caller parameter, defaulted to this published value.
    private const double DefaultNucleotideK = 0.711;

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

    private static double CalculateGcContent(string sequence) =>
        sequence.Length > 0 ? sequence.CalculateGcFractionFast() : 0;

    private static double CalculateTm(string sequence)
    {
        int length = sequence.Length;

        if (length < ThermoConstants.WallaceMaxLength)
        {
            // Wallace rule for short oligos
            int at = sequence.Count(c => c == 'A' || c == 'T');
            int gc = sequence.Count(c => c == 'G' || c == 'C');
            return ThermoConstants.CalculateWallaceTm(at, gc);
        }
        else
        {
            // Salt-adjusted formula
            double gc = CalculateGcContent(sequence);
            return ThermoConstants.CalculateSaltAdjustedTm(gc, length);
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
