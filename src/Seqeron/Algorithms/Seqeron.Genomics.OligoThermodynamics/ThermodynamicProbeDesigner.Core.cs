// SPDX-License-Identifier: GPL-2.0-or-later
// Port of Primer3 (thal.c ntthal, dpal.c, masker.c, libprimer3; primer3-org/primer3 and primer3-py 2.3.1; GNU General Public License v2 or later).
// Part of Seqeron.Genomics.OligoThermodynamics, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.OligoThermodynamics;

/// <summary>Thermodynamic (ProbeDesigner) members of the Primer3 ntthal family; see the licence notes of this assembly.</summary>
public static partial class ThermodynamicProbeDesigner
{
    /// <summary>
    /// Probe design parameters.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Tm.</b> <see cref="MinTm"/>/<see cref="MaxTm"/> are compared with the probe Tm computed by
    /// Primer3's <c>seqtm</c> (<see cref="ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3"/>: SantaLucia 1998
    /// nearest-neighbour Tm with the SantaLucia salt correction for ≤ 36 nt, Primer3 <c>long_seq_tm</c>
    /// 81.5 + 16.6·log10([Mon]_eq) + 0.41·%GC − 600/N for longer probes) at the hybridization conditions
    /// <see cref="DnaConcentrationNanomolar"/>, <see cref="MonovalentMillimolar"/>, <see cref="DivalentMillimolar"/>,
    /// <see cref="DntpMillimolar"/>. The defaults are Primer3's hybridization-probe (internal-oligo) conditions
    /// PRIMER_INTERNAL_DNA_CONC = 50 nM, PRIMER_INTERNAL_SALT_MONOVALENT = 50 mM, PRIMER_INTERNAL_SALT_DIVALENT = 0,
    /// PRIMER_INTERNAL_DNTP_CONC = 0 (<c>libprimer3.cc</c> <c>pr_set_default_global_args</c>). The nearest-neighbour
    /// length limit is <see cref="MaxNearestNeighborLength"/> (Primer3 <c>nn_max_len</c>, default 36). The presets in
    /// <see cref="Defaults"/> state their own conditions; each preset's Tm window is reachable on this scale for its
    /// length and G+C window (audit round 3, A3-10).
    /// </para>
    /// <para>
    /// <b>Self-structure.</b> With <see cref="StructureScreen"/> = <see cref="ProbeStructureScreen.Thermodynamic"/>
    /// (default) a probe of ≤ <see cref="ThermodynamicScreenMaxLength"/> (default 60) nt made of A/C/G/T is screened as Primer3 screens a hybridization probe
    /// (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 1): ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin
    /// Tm (<see cref="ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure"/>) must not exceed <see cref="MaxStructureTm"/>
    /// (PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH = 47 °C). The self-dimer limit replaces
    /// <see cref="MaxSelfComplementarity"/>; the hairpin limit is applied when <see cref="AvoidSecondaryStructure"/>.
    /// Probes longer than <see cref="ThermodynamicScreenMaxLength"/> (thal.c <c>THAL_MAX_ALIGN</c>, default 60; a
    /// larger value is an opt-in that runs the unchanged ntthal recursions on longer probes), probes with non-ACGT
    /// bases and <see cref="ProbeStructureScreen.Heuristic"/> use the fallback screens: the self-dimer criterion is Primer3's
    /// alignment-mode internal-oligo screen (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0, <c>oligo_compl</c>): the dpal
    /// <c>self_any</c> (<see cref="ThermodynamicPrimerDesigner.CalculatePrimerSelfAnyComplementarity"/>) and <c>self_end</c>
    /// (<see cref="PrimerDesigner.CalculatePrimerSelfEndComplementarity"/>) scores, which have no length limit, must not
    /// exceed <see cref="MaxSelfAny"/> / <see cref="MaxSelfEnd"/> (PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END = 12.00);
    /// the hairpin criterion is the canonical sequence-only stem-loop screen
    /// <see cref="PrimerDesigner.HasHairpinPotential"/> (an exactly complementary stem of ≥ 4 bp — library convention,
    /// no published source — closing a loop of ≥ 3 nt, thal.c <c>min_hrpn_loop</c> = 3; audit round 5, A5-2). It is a
    /// conservative presence test, not a stability estimate: almost every random probe of &gt; 60 nt contains such a
    /// stem. It therefore only lowers the additive library score in <see cref="DesignProbes(string, ProbeParameters?, int)"/>
    /// (with <see cref="AvoidSecondaryStructure"/>) and is reported, but does not decide validity, in
    /// <see cref="ValidateProbe"/> (audit round 5, A5-5). With the default THAL_MAX_ALIGN = 60 a &gt; 60-nt probe gets no thermodynamic hairpin; raise
    /// <see cref="ThermodynamicScreenMaxLength"/> for the ntthal hairpin of longer probes (audit round 3, A3-9).
    /// <see cref="MaxSelfComplementarity"/> (position-wise fold-back fraction limit) is
    /// no longer used by any screen; it is kept for source compatibility.
    /// </para>
    /// </remarks>
    public readonly record struct ProbeParameters(
        int MinLength,
        int MaxLength,
        double MinTm,
        double MaxTm,
        double MinGc,
        double MaxGc,
        int MaxHomopolymer,
        bool AvoidSecondaryStructure,
        double MaxSelfComplementarity)
    {
        /// <summary>Probe (oligo) concentration in nM for the Tm (Primer3 PRIMER_INTERNAL_DNA_CONC, default 50).</summary>
        public double DnaConcentrationNanomolar { get; init; } = PrimerDesigner.Primer3InternalDnaConcentrationNanomolar;

        /// <summary>Monovalent cation concentration in mM (Primer3 PRIMER_INTERNAL_SALT_MONOVALENT, default 50).</summary>
        public double MonovalentMillimolar { get; init; } = PrimerDesigner.Primer3InternalMonovalentMillimolar;

        /// <summary>Mg²⁺ concentration in mM (Primer3 PRIMER_INTERNAL_SALT_DIVALENT, default 0).</summary>
        public double DivalentMillimolar { get; init; } = PrimerDesigner.Primer3InternalDivalentMillimolar;

        /// <summary>dNTP concentration in mM (Primer3 PRIMER_INTERNAL_DNTP_CONC, default 0).</summary>
        public double DntpMillimolar { get; init; } = PrimerDesigner.Primer3InternalDntpMillimolar;

        /// <summary>
        /// Longest probe (nt) whose Tm is the SantaLucia 1998 nearest-neighbour Tm; longer probes get Primer3's
        /// <c>long_seq_tm</c> (Primer3 <c>seqtm</c> <c>nn_max_len</c>; default <see cref="PrimerDesigner.Primer3MaxNnTmLength"/>
        /// = 36, Primer3's MAX_NN_TM_LENGTH; ≥ 0).
        /// </summary>
        public int MaxNearestNeighborLength { get; init; } = PrimerDesigner.Primer3MaxNnTmLength;

        /// <summary>Self-structure screen (default <see cref="ProbeStructureScreen.Thermodynamic"/>).</summary>
        public ProbeStructureScreen StructureScreen { get; init; } = ProbeStructureScreen.Thermodynamic;

        /// <summary>
        /// Maximum ntthal self-dimer / 3′ self-dimer / hairpin Tm (°C) of the thermodynamic screen
        /// (Primer3 PRIMER_INTERNAL_MAX_SELF_ANY_TH = _SELF_END_TH = _HAIRPIN_TH = 47 °C).
        /// </summary>
        public double MaxStructureTm { get; init; } = PrimerDesigner.Primer3MaxStructureTm;

        /// <summary>
        /// Longest probe (nt) screened by the <see cref="ProbeStructureScreen.Thermodynamic"/> ntthal screen — thal.h
        /// <c>THAL_MAX_ALIGN</c>. Default <see cref="ThermodynamicPrimerDesigner.NtthalMaxAlignLength"/> = 60 (Primer3 / primer3-py;
        /// unchanged behaviour). THAL_MAX_ALIGN is only a compile-time guard of thal.c, so a larger value (opt-in,
        /// ≤ <see cref="ThermodynamicPrimerDesigner.NtthalMaxSequenceLength"/>) screens longer ACGT probes with the unchanged ntthal
        /// self-dimer / 3′ self-dimer / hairpin recursions (= ntthal compiled with <c>-DTHAL_MAX_ALIGN=…</c>) instead of
        /// the fallback screens. Cost grows as O(n²·30²) per probe and screen — keep it to the probe lengths in use.
        /// Also the THAL_MAX_ALIGN of the <see cref="AssessCrossHybridization"/> / <see cref="ValidateProbe"/> non-target site
        /// duplex Tm (audit round 3, A3-27). Outside 60–10 000 → <see cref="ArgumentOutOfRangeException"/>.
        /// </summary>
        public int ThermodynamicScreenMaxLength { get; init; } = ThermodynamicPrimerDesigner.NtthalMaxAlignLength;

        /// <summary>
        /// Maximum Primer3 alignment-mode <c>self_any</c> of the fallback self-dimer screen
        /// (PRIMER_INTERNAL_MAX_SELF_ANY, Primer3 default 12.00; flagged when strictly greater).
        /// </summary>
        public double MaxSelfAny { get; init; } = PrimerDesigner.Primer3InternalMaxSelfComplementarity;

        /// <summary>
        /// Maximum Primer3 alignment-mode <c>self_end</c> of the fallback self-dimer screen
        /// (PRIMER_INTERNAL_MAX_SELF_END, Primer3 default 12.00; flagged when strictly greater).
        /// </summary>
        public double MaxSelfEnd { get; init; } = PrimerDesigner.Primer3InternalMaxSelfComplementarity;

        /// <summary>
        /// Ranking of <see cref="DesignProbes(string, ProbeParameters?, int)"/> (default
        /// <see cref="ProbeRanking.AdditiveScore"/>, the library heuristic). <see cref="ProbeRanking.Primer3Penalty"/>
        /// ranks by Primer3's internal-oligo objective <c>p_obj_fn</c> (PRIMER_INTERNAL_n_PENALTY) around
        /// <see cref="OptTm"/> / <see cref="OptLength"/> — the sourced ranking.
        /// </summary>
        public ProbeRanking Ranking { get; init; } = ProbeRanking.AdditiveScore;

        /// <summary>
        /// Optimum probe Tm (°C) of the <see cref="ProbeRanking.Primer3Penalty"/> ranking (Primer3 PRIMER_INTERNAL_OPT_TM,
        /// default 60 °C, <c>pr_set_default_global_args</c> <c>o_args.opt_tm</c>). With that ranking it must lie in
        /// [<see cref="MinTm"/>, <see cref="MaxTm"/>] (Primer3 <c>_pr_data_control</c>: "Optimum internal oligo Tm lower
        /// than minimum or higher than maximum"); unused by the additive score.
        /// </summary>
        public double OptTm { get; init; } = ThermodynamicPrimerDesigner.DefaultPrimer3Optima.OptTm;

        /// <summary>
        /// Optimum probe length (nt) of the <see cref="ProbeRanking.Primer3Penalty"/> ranking (Primer3 PRIMER_INTERNAL_OPT_SIZE,
        /// default 20, <c>o_args.opt_size</c>). With that ranking it must lie in [<see cref="MinLength"/>, <see cref="MaxLength"/>]
        /// (Primer3 <c>_pr_data_control</c>: "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE &gt; MAX_SIZE" / "&lt; MIN_SIZE");
        /// unused by the additive score.
        /// </summary>
        public int OptLength { get; init; } = ThermodynamicPrimerDesigner.DefaultPrimer3Optima.OptSize;
    }

    // ProbeParameters.ThermodynamicScreenMaxLength: THAL_MAX_ALIGN override range 60 (Primer3) … 10 000 (THAL_MAX_SEQ).
    private static void ValidateThermodynamicScreenMaxLength(ProbeParameters param, string paramName)
    {
        if (param.ThermodynamicScreenMaxLength < ThermodynamicPrimerDesigner.NtthalMaxAlignLength
            || param.ThermodynamicScreenMaxLength > ThermodynamicPrimerDesigner.NtthalMaxSequenceLength)
            throw new ArgumentOutOfRangeException(paramName,
                "ThermodynamicScreenMaxLength (THAL_MAX_ALIGN) must be in 60..10000.");
    }

    /// <summary>
    /// Designs a molecular beacon (Tyagi &amp; Kramer 1996): a target-complementary loop of
    /// <paramref name="probeLength"/> bases flanked by complementary GC-rich arms
    /// (5′ arm = ⌊stem/2⌋ G + remaining C, 3′ arm = its reverse complement).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The loop is the window with the best score: −0.2 each for GC outside 40–60 %, loop Tm outside the
    /// Tm window and a homopolymer run &gt; 4. The loop (probe–target) Tm is Primer3 <c>seqtm</c> at the Primer3
    /// hybridization-probe conditions (50 nM, 50 mM monovalent, no Mg²⁺/dNTP; <see cref="ProbeParameters"/>).
    /// With <paramref name="detectionTemperatureCelsius"/> = T the Tm window is the published
    /// [T + 7, T + 10] °C (molecular-beacon design rules of Tyagi &amp; Kramer / Marras et al.: the probe Tm and
    /// the stem Tm 7–10 °C above the detection temperature) and the beacon's stem-loop (hairpin) Tm is
    /// checked against T + 7 °C; without it the library's default window 55–65 °C is used.
    /// </para>
    /// <para>
    /// Warnings carry the stem/loop sizes and, for a beacon of ≤ <c>maxAlignLength</c> (default 60) nt, the ntthal hairpin Tm of the whole
    /// beacon (the stem melting temperature; <see cref="ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(string, double, double, double)"/>
    /// at the same conditions). <see cref="Probe.Tm"/> is the loop (probe–target) Tm.
    /// </para>
    /// </remarks>
    /// <param name="targetSequence">Target sequence (the loop is a window of it).</param>
    /// <param name="probeLength">Loop length (Tyagi &amp; Kramer: 15–30 nt; default 25).</param>
    /// <param name="stemLength">Arm length in bp (Tyagi &amp; Kramer: 5–7 bp; default 5).</param>
    /// <param name="detectionTemperatureCelsius">Optional detection (annealing) temperature T in °C.</param>
    /// <param name="maxAlignLength">ntthal <c>THAL_MAX_ALIGN</c> of the stem-loop hairpin Tm
    /// (<see cref="ThermodynamicPrimerDesigner.NtthalMaxAlignLength"/> = 60, Primer3's compile-time default, …
    /// <see cref="ThermodynamicPrimerDesigner.NtthalMaxSequenceLength"/> = 10 000): beacons up to this length get the ntthal hairpin
    /// Tm (opt-in for beacons longer than 60 nt; = thal.c compiled with <c>-DTHAL_MAX_ALIGN=…</c>; audit round 3, A3-27).</param>
    /// <exception cref="ArgumentNullException"><paramref name="targetSequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="probeLength"/> or <paramref name="stemLength"/> ≤ 0
    /// (a beacon needs a loop and two arms; the Tyagi &amp; Kramer ranges are recommendations and are not enforced),
    /// <paramref name="stemLength"/> so large that the beacon (2·stem + loop) cannot be built,
    /// or <paramref name="maxAlignLength"/> outside 60–10 000.</exception>
    public static Probe? DesignMolecularBeacon(
        string targetSequence,
        int probeLength = 25,
        int stemLength = 5,
        double? detectionTemperatureCelsius = null,
        int maxAlignLength = ThermodynamicPrimerDesigner.NtthalMaxAlignLength)
    {
        ArgumentNullException.ThrowIfNull(targetSequence);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(probeLength);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(stemLength);
        // The beacon (two arms + loop) must fit in a string (heavy tier, HEAVY-2: stemLength = int.MaxValue overflowed
        // into an OutOfMemoryException).
        if (2L * stemLength + probeLength > Array.MaxLength)
            throw new ArgumentOutOfRangeException(nameof(stemLength), stemLength, "The beacon (2·stemLength + probeLength) is too long.");
        if (maxAlignLength < ThermodynamicPrimerDesigner.NtthalMaxAlignLength || maxAlignLength > ThermodynamicPrimerDesigner.NtthalMaxSequenceLength)
            throw new ArgumentOutOfRangeException(nameof(maxAlignLength), maxAlignLength,
                "THAL_MAX_ALIGN must be in 60..10000.");
        if (targetSequence.Length < probeLength)
            return null;

        targetSequence = targetSequence.ToUpperInvariant();
        var conditions = Primer3ProbeConditions; // only the (default Primer3 internal-oligo) conditions are used

        double minLoopTm = detectionTemperatureCelsius is { } t0 ? t0 + BeaconMinTmAboveDetection : 55;
        double maxLoopTm = detectionTemperatureCelsius is { } t1 ? t1 + BeaconMaxTmAboveDetection : 65;

        // Find best region in target
        double bestScore = 0;
        string? bestLoop = null;
        int bestStart = 0;
        double bestTm = 0;

        int loopLength = probeLength;
        for (int start = 0; start <= targetSequence.Length - loopLength; start++)
        {
            string loop = targetSequence.Substring(start, loopLength);
            double gc = loop.CalculateGcFractionFast();
            double tmRaw = CalculateProbeTm(loop, conditions);
            double tm = double.IsNaN(tmRaw) ? 0.0 : tmRaw;

            double score = 1.0;
            if (gc < 0.40 || gc > 0.60) score -= 0.2;
            if (double.IsNaN(tmRaw) || tm < minLoopTm || tm > maxLoopTm) score -= 0.2;
            if (PrimerDesigner.FindLongestHomopolymer(loop) > 4) score -= 0.2;

            if (score > bestScore)
            {
                bestScore = score;
                bestLoop = loop;
                bestStart = start;
                bestTm = tm;
            }
        }

        if (bestLoop == null)
            return null;

        // Add stem sequences (GC-rich for stability)
        string stem5 = new string('G', stemLength / 2) + new string('C', stemLength - stemLength / 2);
        string stem3 = DnaSequence.GetReverseComplementString(stem5);

        string beaconSequence = stem5 + bestLoop + stem3;
        var warnings = new List<string> { $"Stem: {stemLength}bp, Loop: {loopLength}bp" };

        if (beaconSequence.Length <= maxAlignLength)
        {
            var hp = ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(
                beaconSequence,
                conditions.MonovalentMillimolar / 1000.0,
                conditions.DivalentMillimolar / 1000.0,
                conditions.DntpMillimolar / 1000.0,
                ThermodynamicPrimerDesigner.NtthalDefaultTemperatureCelsius,
                ThermodynamicPrimerDesigner.NtthalDefaultMaxLoop,
                maxAlignLength);
            if (hp is { } h)
            {
                warnings.Add($"Stem-loop (hairpin) Tm {h.TmCelsius:F1}°C (ntthal)");
                if (detectionTemperatureCelsius is { } t && h.TmCelsius < t + BeaconMinTmAboveDetection)
                    warnings.Add(
                        $"Stem-loop Tm {h.TmCelsius:F1}°C is less than {BeaconMinTmAboveDetection:0} °C above the detection temperature {t:0.#}°C");
            }
        }

        if (detectionTemperatureCelsius is { } td && bestTm < td + BeaconMinTmAboveDetection)
            warnings.Add(
                $"Probe Tm {bestTm:F1}°C is less than {BeaconMinTmAboveDetection:0} °C above the detection temperature {td:0.#}°C");

        return new Probe(
            beaconSequence,
            bestStart,
            bestStart + loopLength - 1,
            bestTm,
            beaconSequence.CalculateGcFractionFast(),
            bestScore,
            ProbeType.MolecularBeacon,
            warnings);
    }

    /// <summary>
    /// Kane et al. (2000) cross-hybridization criteria of a probe against known non-target sequences: for each
    /// non-target strand the probe may cross-hybridize when (a) its overall identity — identical columns of the
    /// best local alignment ÷ probe length — is strictly above <paramref name="maxIdentity"/> (Kane: &gt; 75 %
    /// similar over the probe), or (b) the longest stretch identical between probe and non-target is strictly
    /// longer than <paramref name="maxContiguousMatch"/> (Kane: &gt; 15 contiguous bases).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The local alignment is the canonical Smith–Waterman–Gotoh aligner
    /// <see cref="SequenceAligner.LocalAlignAffine(string, string, ScoringMatrix?)"/> with BLAST+ blastn scoring
    /// (<see cref="SequenceAligner.BlastDna"/>: +2/−3, gap existence 5, extension 2) by default — the reported
    /// optimal alignment is one of those Biopython's local <c>PairwiseAligner</c> enumerates. Non-target strands
    /// longer than 4096 nt are aligned in chunks overlapping by the longest reference span a positive-scoring
    /// local alignment can have (m + ⌈m·match/|extend|⌉), so the best score equals the whole-strand score; among
    /// equal-scoring chunks the first is reported. The longest contiguous match is the longest common substring
    /// from the canonical suffix tree (<c>SuffixTree.LongestCommonSubstringInfo</c>). Comparison is
    /// case-insensitive and literal (an N matches only N).
    /// </para>
    /// <para>
    /// With <paramref name="bothStrands"/> (default) each non-target is assessed as given and as its reverse
    /// complement, since a double-stranded non-target (genomic DNA, ds cDNA) offers the probe both strands;
    /// pass false for single-stranded non-targets given in the probe's sense.
    /// </para>
    /// </remarks>
    /// <param name="probeSequence">Probe sequence (non-empty).</param>
    /// <param name="nonTargetSequences">Known non-target sequences (null entries are treated as empty).</param>
    /// <param name="maxIdentity">Identity threshold in [0, 1] (default 0.75).</param>
    /// <param name="maxContiguousMatch">Contiguous-identity threshold in nt (≥ 0; default 15).</param>
    /// <param name="bothStrands">Assess the reverse complement of each non-target too (default true).</param>
    /// <param name="scoring">Local-alignment scoring (default <see cref="SequenceAligner.BlastDna"/>, affine gaps).</param>
    /// <param name="conditions">Hybridization conditions for the site duplex Tm
    /// (<see cref="CrossHybridizationAssessment.DuplexTm"/>; default Primer3 probe conditions 50 nM / 50 mM / 0 / 0).
    /// Its <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> is the ntthal <c>THAL_MAX_ALIGN</c> of the
    /// duplex (default 60; opt-in up to 10 000 for longer probes and sites — audit round 3, A3-27).</param>
    /// <param name="maxDuplexTm">Optional OligoArray-style specificity threshold (°C): a strand whose site duplex Tm
    /// is strictly above it is flagged (<see cref="CrossHybridizationAssessment.ExceedsDuplexTmThreshold"/>); null
    /// (default) applies only the Kane criteria. The threshold is assay-specific (Rouillard et al. 2003: user-set).</param>
    /// <returns>One assessment per non-target strand: index order, forward strand before reverse complement.</returns>
    /// <exception cref="ArgumentNullException">A null probe or non-target collection.</exception>
    /// <exception cref="ArgumentException">An empty probe.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxIdentity"/> outside [0, 1], a negative
    /// <paramref name="maxContiguousMatch"/>, or a <paramref name="conditions"/>
    /// <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> outside 60–10 000.</exception>
    public static IReadOnlyList<CrossHybridizationAssessment> AssessCrossHybridization(
        string probeSequence,
        IEnumerable<string> nonTargetSequences,
        double maxIdentity = KaneMaxIdentity,
        int maxContiguousMatch = KaneMaxContiguousMatch,
        bool bothStrands = true,
        ScoringMatrix? scoring = null,
        ProbeParameters? conditions = null,
        double? maxDuplexTm = null)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(nonTargetSequences);
        if (probeSequence.Length == 0)
            throw new ArgumentException("Probe sequence cannot be empty.", nameof(probeSequence));
        if (double.IsNaN(maxIdentity) || maxIdentity < 0 || maxIdentity > 1)
            throw new ArgumentOutOfRangeException(nameof(maxIdentity), "Identity threshold must be in [0, 1].");
        if (maxContiguousMatch < 0)
            throw new ArgumentOutOfRangeException(nameof(maxContiguousMatch), "Contiguous-match threshold cannot be negative.");

        if (conditions is { } stated)
            ValidateThermodynamicScreenMaxLength(stated, nameof(conditions));

        var matrix = scoring ?? SequenceAligner.BlastDna;
        var cond = conditions ?? Primer3ProbeConditions;
        string probe = probeSequence.ToUpperInvariant();
        // thal.c thal_check_errors: at least one strand ≤ THAL_MAX_ALIGN (ThermodynamicScreenMaxLength, default 60)
        // and neither longer than THAL_MAX_SEQ (10 000).
        int maxAlign = cond.ThermodynamicScreenMaxLength;
        bool duplexComputable = probe.Length <= ThermodynamicPrimerDesigner.NtthalMaxSequenceLength
            && probe.All(c => c is 'A' or 'C' or 'G' or 'T');
        var probeTree = global::SuffixTree.SuffixTree.Build(probe);
        var result = new List<CrossHybridizationAssessment>();

        int index = 0;
        foreach (var nonTarget in nonTargetSequences)
        {
            string forward = (nonTarget ?? string.Empty).ToUpperInvariant();
            result.Add(AssessStrand(index, false, forward));
            if (bothStrands)
                result.Add(AssessStrand(index, true, DnaSequence.GetReverseComplementString(forward)));
            index++;
        }

        return result;

        CrossHybridizationAssessment AssessStrand(int idx, bool reverse, string strand)
        {
            var (score, identical, siteStart, siteEnd) = BestLocalAlignment(probe, strand, matrix);
            int contiguous = strand.Length == 0 ? 0 : probeTree.LongestCommonSubstringInfo(strand).Substring.Length;
            double identity = (double)identical / probe.Length;
            double? duplexTm = null;
            if (duplexComputable && siteStart >= 0)
            {
                string site = strand.Substring(siteStart, siteEnd - siteStart + 1);
                if (Math.Min(probe.Length, site.Length) <= maxAlign
                    && site.Length <= ThermodynamicPrimerDesigner.NtthalMaxSequenceLength
                    && site.All(c => c is 'A' or 'C' or 'G' or 'T'))
                {
                    // The probe hybridizes to the strand complementary to the site (primer3-py calc_heterodimer).
                    var d = ThermodynamicPrimerDesigner.CalculateDimerThermodynamicsNtthal(
                        probe, DnaSequence.GetReverseComplementString(site), ThermodynamicPrimerDesigner.NtthalAlignmentMode.Any,
                        cond.MonovalentMillimolar / 1000.0, cond.DivalentMillimolar / 1000.0, cond.DntpMillimolar / 1000.0,
                        cond.DnaConcentrationNanomolar * 1e-9, ThermodynamicPrimerDesigner.NtthalDefaultTemperatureCelsius,
                        ThermodynamicPrimerDesigner.NtthalDefaultMaxLoop, maxAlign);
                    duplexTm = d?.TmCelsius ?? 0.0;
                }
            }

            return new CrossHybridizationAssessment(
                idx, reverse, identity, identical, score, contiguous,
                identity > maxIdentity, contiguous > maxContiguousMatch)
            {
                SiteStart = siteStart,
                SiteEnd = siteEnd,
                DuplexTm = duplexTm,
                ExceedsDuplexTmThreshold = maxDuplexTm is double limit && duplexTm is double tm && tm > limit,
            };
        }
    }

    /// <summary>
    /// Default probe parameters for different applications.
    /// </summary>
    /// <remarks>
    /// Every preset's Tm window is reachable on its own Tm scale and conditions for its length and G+C window
    /// (audit round 3, A3-10; primer3-py 2.3.1 <c>calc_tm</c> witnesses in the PROBE-DESIGN-001 tests). Primer3
    /// <c>long_seq_tm</c> (the Tm of every probe &gt; <see cref="ProbeParameters.MaxNearestNeighborLength"/>) at the
    /// 50 mM Primer3 probe conditions spans 71.25–85.35 °C over the FISH window (200–500 nt, G+C 35–65 %),
    /// 70.30–82.50 °C over the Northern window (100–300 nt, G+C 40–60 %) and 70.32–85.35 °C over the Southern window
    /// (150–500 nt, G+C 35–65 %); their Tm windows are library conventions (no published Tm window for these long
    /// probes was found — set them for the protocol's hybridization buffer). The qPCR window (68–70 °C, Primer Express)
    /// is reachable for 20–30-nt probes with G+C 30–80 % (nearest-neighbour Tm at the Primer3 probe conditions).
    /// </remarks>
    public static class Defaults
    {
        /// <summary>OligoArray 2.0 microarray Tm conditions: [Na⁺] = 1 M, oligo 1 µM (Rouillard, Zuker &amp; Gulari 2003).</summary>
        public const double OligoArrayMonovalentMillimolar = 1000.0;

        /// <summary>OligoArray 2.0 oligo concentration for the Tm, 1 µM = 1000 nM (Rouillard, Zuker &amp; Gulari 2003).</summary>
        public const double OligoArrayDnaConcentrationNanomolar = 1000.0;

        /// <summary>
        /// Long-oligo microarray probe preset (50–60 nt: Kane et al. 2000 50-mers, Agilent 60-mers; G+C 40–60 %) with
        /// OligoArray 2.0's Tm (Rouillard, Zuker &amp; Gulari 2003, NAR 31:3057): nearest-neighbour Tm over the whole
        /// oligo at [Na⁺] = 1 M and 1 µM oligo, Tm window 82–90 °C (the paper's microarray design setting). Here the
        /// nearest-neighbour Tm is Primer3 <c>oligotm</c> (SantaLucia 1998 unified parameters, C/4 for a
        /// non-self-complementary oligo) with <see cref="ProbeParameters.MaxNearestNeighborLength"/> = 60, no Mg²⁺/dNTP.
        /// The previous window 75–85 °C at the 50 mM Primer3 conditions was unreachable for 50–60-mers with G+C ≤ 60 %
        /// (<c>long_seq_tm</c> maximum 74.50 °C). The structure screen runs at the same conditions.
        /// </summary>
        public static ProbeParameters Microarray => new(
            MinLength: 50, MaxLength: 60,
            MinTm: 82, MaxTm: 90,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 5,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.3)
        {
            MonovalentMillimolar = OligoArrayMonovalentMillimolar,
            DnaConcentrationNanomolar = OligoArrayDnaConcentrationNanomolar,
            DivalentMillimolar = 0,
            DntpMillimolar = 0,
            MaxNearestNeighborLength = 60
        };

        /// <summary>FISH probe preset (200–500 nt, G+C 35–65 %, Tm 70–90 °C at the Primer3 probe conditions; library
        /// convention, reachable: <c>long_seq_tm</c> 71.25–85.35 °C over the window).</summary>
        public static ProbeParameters FISH => new(
            MinLength: 200, MaxLength: 500,
            MinTm: 70, MaxTm: 90,
            MinGc: 0.35, MaxGc: 0.65,
            MaxHomopolymer: 8,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 0.4);

        /// <summary>Northern-blot probe preset (100–300 nt, G+C 40–60 %, Tm 65–80 °C at the Primer3 probe conditions;
        /// library convention, reachable: <c>long_seq_tm</c> 70.30–82.50 °C over the window).</summary>
        public static ProbeParameters NorthernBlot => new(
            MinLength: 100, MaxLength: 300,
            MinTm: 65, MaxTm: 80,
            MinGc: 0.40, MaxGc: 0.60,
            MaxHomopolymer: 6,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.35);

        /// <summary>
        /// Hydrolysis (TaqMan) qPCR probe preset: Tm 68–70 °C and G+C 30–80 % (Applied Biosystems
        /// Primer Express probe guidelines: "Tm should be 68 to 70 °C", "Keep G-C content in the 30-80%
        /// range"; the Primer Express Tm scale itself is proprietary — here the Tm is Primer3 <c>seqtm</c>
        /// at <see cref="ProbeParameters.MonovalentMillimolar"/> etc.); length 20–30 nt within the
        /// 18–30 nt hydrolysis-probe range (IDT). Use <see cref="EvaluateTaqManProbe"/> for the
        /// chemistry-specific rules (no 5′ G, more C than G, no GGGG, probe Tm ≥ primer Tm + 10 °C).
        /// </summary>
        public static ProbeParameters qPCR => new(
            MinLength: 20, MaxLength: 30,
            MinTm: 68, MaxTm: 70,
            MinGc: 0.30, MaxGc: 0.80,
            MaxHomopolymer: 4,
            AvoidSecondaryStructure: true,
            MaxSelfComplementarity: 0.25);

        /// <summary>Southern-blot probe preset (150–500 nt, G+C 35–65 %, Tm 65–75 °C at the Primer3 probe conditions;
        /// library convention, reachable: <c>long_seq_tm</c> 70.32–85.35 °C over the window).</summary>
        public static ProbeParameters SouthernBlot => new(
            MinLength: 150, MaxLength: 500,
            MinTm: 65, MaxTm: 75,
            MinGc: 0.35, MaxGc: 0.65,
            MaxHomopolymer: 7,
            AvoidSecondaryStructure: false,
            MaxSelfComplementarity: 0.4);
    }

    /// <summary>
    /// Designs probes for a target sequence: every window of <see cref="ProbeParameters.MinLength"/>..<see cref="ProbeParameters.MaxLength"/>
    /// (GC within ±0.1 of the window) is scored 1.0 minus penalties — GC outside range 0.3, Tm (Primer3 <c>seqtm</c> at the
    /// parameters' conditions) outside range 0.3, homopolymer run above the maximum 0.2, self-complementarity 0.2,
    /// secondary structure 0.15 (only with <see cref="ProbeParameters.AvoidSecondaryStructure"/>), di-/trinucleotide
    /// microsatellite ≥ 4 copies 0.1, G/C at each terminus 0.02 (see <see cref="ProbeParameters"/> for the self-structure
    /// screens). Probes with score &gt; 0 are returned by score descending, ties by (length, start) ascending. This
    /// additive score (the default, <see cref="ProbeRanking.AdditiveScore"/>) is a library ranking heuristic: its penalty
    /// values have no published source. With <see cref="ProbeParameters.Ranking"/> = <see cref="ProbeRanking.Primer3Penalty"/>
    /// the same candidates (score &gt; 0) are instead ranked by Primer3's internal-oligo objective <c>p_obj_fn</c>
    /// (PRIMER_INTERNAL_n_PENALTY = |Tm − <see cref="ProbeParameters.OptTm"/>| + |length − <see cref="ProbeParameters.OptLength"/>|
    /// with Primer3's default internal-oligo weights, reported in <see cref="Probe.Primer3Penalty"/>) in Primer3's
    /// <c>primer_rec_comp</c> order (penalty ascending, start descending, length ascending; a probe whose Tm is not
    /// computable ranks last) — the sourced ranking. <see cref="DesignProbesPrimer3"/> is Primer3's complete picker
    /// (Primer3's acceptance limits as well as its ranking).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">With the <see cref="ProbeRanking.Primer3Penalty"/> ranking:
    /// <see cref="ProbeParameters.OptLength"/> outside [MinLength, MaxLength] or <see cref="ProbeParameters.OptTm"/> outside
    /// [MinTm, MaxTm] (Primer3 <c>_pr_data_control</c>); an undefined <see cref="ProbeParameters.Ranking"/> value;
    /// <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> outside 60–10 000; <see cref="ProbeParameters.MinLength"/> &lt; 1 (Primer3 <c>_pr_data_control</c>: PRIMER_INTERNAL_MIN_SIZE must be ≥ 1), a negative <see cref="ProbeParameters.MaxNearestNeighborLength"/>, or illegal reaction conditions (monovalent salt or DNA concentration not &gt; 0, divalent salt or dNTP negative, or non-finite — Primer3 <c>_pr_data_control</c>; audit heavy tier, HEAVY-2).</exception>
    public static IEnumerable<Probe> DesignProbes(
        string targetSequence,
        ProbeParameters? parameters = null,
        int maxProbes = 10)
    {
        var param = parameters ?? Defaults.Microarray;
        ValidateRanking(param, nameof(parameters));
        ValidateDesignParameters(param, nameof(parameters), checkLengths: true);
        return DesignProbesIterator(targetSequence, param, maxProbes);
    }

    private static IEnumerable<Probe> DesignProbesIterator(string targetSequence, ProbeParameters param, int maxProbes)
    {
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
    /// Designs probes with a genome-wide specificity check through a suffix tree: the candidates of
    /// <see cref="DesignProbes(string, ProbeParameters?, int)"/> are walked in their ranking order and each one's
    /// <see cref="CheckSpecificity(string, global::SuffixTree.ISuffixTree, bool)"/> value (1 / occurrences in the index)
    /// is applied — with <paramref name="requireUnique"/> a probe occurring more than once is dropped (the remaining
    /// probes keep their score and order), otherwise <see cref="Probe.Score"/> is multiplied by the specificity and the
    /// probes are re-ranked on the scaled score (stable: equal scores keep the documented tie order; the
    /// <see cref="ProbeRanking.Primer3Penalty"/> order is unaffected because the specificity does not enter the penalty).
    /// <b>Every</b> candidate is considered, lazily: the walk stops as soon as <paramref name="maxProbes"/> probes are
    /// produced, so the common case costs no more than the base design plus O(m) per inspected candidate.
    /// O(n × m) for probe generation + O(m) per specificity check.
    /// </summary>
    /// <param name="targetSequence">Target sequence to design probes for.</param>
    /// <param name="genomeIndex">Pre-built suffix tree index for the genome (enables O(m) specificity lookup).</param>
    /// <param name="parameters">Probe design parameters.</param>
    /// <param name="maxProbes">Maximum number of probes to return (none when ≤ 0).</param>
    /// <param name="requireUnique">If true, only return probes unique in the genome.</param>
    /// <param name="bothStrands">Count the probe's reverse-complement occurrences in the index too (passed to
    /// <see cref="CheckSpecificity(string, global::SuffixTree.ISuffixTree, bool)"/>: a probe whose reverse complement
    /// occurs in a double-stranded genome binds the other strand there; a reverse-palindromic probe is counted once).
    /// Default false: the indexed strand only (the behaviour before audit round 4, B07 F60) — pass true when the index
    /// holds one strand of a double-stranded genome.</param>
    /// <remarks>Candidates come in the order of <see cref="ProbeParameters.Ranking"/> (see
    /// <see cref="DesignProbes(string, ProbeParameters?, int)"/>), after the specificity scaling when
    /// <paramref name="requireUnique"/> is false. Before audit round 4 (B07 F59) only the top
    /// <paramref name="maxProbes"/> × 5 candidates were inspected and the scaled scores were not re-ranked.</remarks>
    /// <exception cref="ArgumentOutOfRangeException">Invalid Primer3-ranking optima or other invalid parameters (see
    /// <see cref="DesignProbes(string, ProbeParameters?, int)"/>).</exception>
    public static IEnumerable<Probe> DesignProbes(
        string targetSequence,
        global::SuffixTree.ISuffixTree genomeIndex,
        ProbeParameters? parameters = null,
        int maxProbes = 10,
        bool requireUnique = true,
        bool bothStrands = false)
    {
        var param = parameters ?? Defaults.Microarray;
        ValidateRanking(param, nameof(parameters));
        ValidateDesignParameters(param, nameof(parameters), checkLengths: true);
        return DesignProbesIterator(targetSequence, genomeIndex, param, maxProbes, requireUnique, bothStrands);
    }

    private static IEnumerable<Probe> DesignProbesIterator(
        string targetSequence,
        global::SuffixTree.ISuffixTree genomeIndex,
        ProbeParameters param,
        int maxProbes,
        bool requireUnique,
        bool bothStrands)
    {
        if (string.IsNullOrEmpty(targetSequence) || targetSequence.Length < param.MinLength || maxProbes <= 0)
            yield break;

        targetSequence = targetSequence.ToUpperInvariant();

        // Every candidate is considered, lazily in ranking order, until maxProbes survive (audit round 4, A4-1, F59):
        // requireUnique drops a probe with specificity < 1 (its score is unchanged, so the order is the base ranking);
        // otherwise the score is scaled by the specificity (never raised) and EnumerateRankedProbes re-ranks on it.
        Func<Probe, Probe?> applySpecificity = requireUnique
            ? probe => CheckSpecificity(probe.Sequence, genomeIndex, bothStrands) < 1.0 ? null : probe
            : probe => probe with { Score = probe.Score * CheckSpecificity(probe.Sequence, genomeIndex, bothStrands) };

        int returned = 0;
        foreach (var probe in EnumerateRankedProbes(targetSequence, param, applySpecificity))
        {
            yield return probe;
            if (++returned >= maxProbes)
                yield break;
        }
    }

    // Primer3 _pr_data_control checks of the optima used by the Primer3Penalty ranking.
    private static void ValidateRanking(ProbeParameters param, string paramName)
    {
        if (!Enum.IsDefined(param.Ranking))
            throw new ArgumentOutOfRangeException(paramName, $"Undefined probe ranking {param.Ranking}.");
        if (param.Ranking != ProbeRanking.Primer3Penalty)
            return;
        if (param.OptLength > param.MaxLength)
            throw new ArgumentOutOfRangeException(paramName,
                "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE > MAX_SIZE (Primer3 _pr_data_control): OptLength must not exceed MaxLength.");
        if (param.OptLength < param.MinLength)
            throw new ArgumentOutOfRangeException(paramName,
                "PRIMER_INTERNAL_{OPT,DEFAULT}_SIZE < MIN_SIZE (Primer3 _pr_data_control): OptLength must not be below MinLength.");
        if (!(param.OptTm >= param.MinTm && param.OptTm <= param.MaxTm))
            throw new ArgumentOutOfRangeException(paramName,
                "Optimum internal oligo Tm lower than minimum or higher than maximum (Primer3 _pr_data_control): OptTm must lie in [MinTm, MaxTm].");
    }

    // Primer3 p_obj_fn, internal-oligo branch (OT_INTL) with Primer3's default internal-oligo weights (Tm and size
    // weights 1, every other weight 0) around the parameters' optima; null when the Tm is not computable.
    private static double? ComputePrimer3ProbePenalty(double tm, int length, double gcFraction, ProbeParameters param) =>
        double.IsNaN(tm)
            ? null
            : ThermodynamicPrimerDesigner.CalculatePrimer3Penalty(
                new Primer3PenaltyInputs(tm, length, 100.0 * gcFraction),
                ThermodynamicPrimerDesigner.DefaultPrimer3Weights,
                new Primer3Optima(param.OptTm, param.OptLength, null));

    private static ProbeBase EvaluateProbeBase(string sequence, int start, int order, ProbeParameters param, double gc)
    {
        double tmRaw = CalculateProbeTm(sequence, param);
        bool tmComputable = !double.IsNaN(tmRaw);
        double tm = tmComputable ? tmRaw : 0.0;

        bool gcOutside = gc < param.MinGc || gc > param.MaxGc;
        bool tmOutside = !tmComputable || tm < param.MinTm || tm > param.MaxTm;
        int homopolymer = PrimerDesigner.FindLongestHomopolymer(sequence);
        bool repeats = HasSimpleRepeats(sequence);

        double positionPenalty = 0;
        if (sequence.StartsWith('G') || sequence.StartsWith('C'))
            positionPenalty += 0.02;
        if (sequence.EndsWith('G') || sequence.EndsWith('C'))
            positionPenalty += 0.02;

        double score = 1.0;
        if (gcOutside) score -= 0.3;
        if (tmOutside) score -= 0.3;
        if (homopolymer > param.MaxHomopolymer) score -= 0.2;
        if (repeats) score -= 0.1;
        score -= positionPenalty;

        double? primer3Penalty = param.Ranking == ProbeRanking.Primer3Penalty
            ? ComputePrimer3ProbePenalty(tmRaw, sequence.Length, gc, param)
            : null;

        return new ProbeBase(sequence, start, order, gc, tm, tmComputable, gcOutside, tmOutside,
            homopolymer, repeats, positionPenalty, score, primer3Penalty);
    }

    // Adds the self-structure screens to a base evaluation (penalties in the original order:
    // GC, Tm, homopolymer, self-complementarity, secondary structure, repeats, terminal G/C).
    private static Probe? FinishProbe(
        ProbeBase b, ProbeParameters param,
        Dictionary<string, (bool, string?, bool, string?)>? structureCache = null)
    {
        var warnings = new List<string>();
        double score = 1.0;

        if (b.GcOutside)
        {
            score -= 0.3;
            warnings.Add($"GC content {b.Gc:P0} outside range");
        }

        if (b.TmOutside)
        {
            score -= 0.3;
            warnings.Add(b.TmComputable
                ? $"Tm {b.Tm:F1}°C outside range"
                : "Tm not computable: probe contains a non-ACGT base");
        }

        if (b.Homopolymer > param.MaxHomopolymer)
        {
            score -= 0.2;
            warnings.Add($"Homopolymer run of {b.Homopolymer}");
        }

        (bool, string?, bool, string?) screen;
        if (structureCache is null || !structureCache.TryGetValue(b.Sequence, out screen))
        {
            screen = EvaluateSelfStructure(b.Sequence, param);
            structureCache?.Add(b.Sequence, screen);
        }
        var (selfComp, selfCompWarning, structure, structureWarning) = screen;
        if (selfComp)
        {
            score -= 0.2;
            warnings.Add(selfCompWarning!);
        }

        if (param.AvoidSecondaryStructure && structure)
        {
            score -= 0.15;
            warnings.Add(structureWarning!);
        }

        if (b.HasRepeats)
        {
            score -= 0.1;
            warnings.Add("Contains simple repeats");
        }

        score -= b.PositionPenalty;

        if (score <= 0)
            return null;

        return new Probe(
            b.Sequence,
            b.Start,
            b.Start + b.Sequence.Length - 1,
            b.Tm,
            b.Gc,
            Math.Max(0, score),
            ProbeType.Standard,
            warnings)
        {
            Primer3Penalty = b.Primer3Penalty,
        };
    }

    // Self-structure screens (see ProbeParameters remarks). Returns (self-complementarity flag + warning,
    // secondary-structure flag + warning); the caller applies the secondary-structure flag only when
    // AvoidSecondaryStructure is set.
    private static (bool SelfComp, string? SelfCompWarning, bool Structure, string? StructureWarning)
        EvaluateSelfStructure(string sequence, ProbeParameters param)
    {
        if (ComputeThermodynamicSelfStructure(sequence, param) is { } v)
        {
            {
                double max = param.MaxStructureTm;
                bool selfComp = v.SelfAnyTh > max || v.SelfEndTh > max;
                bool hairpin = v.HairpinTh > max;
                return (
                    selfComp,
                    selfComp ? $"Self-dimer Tm {Math.Max(v.SelfAnyTh, v.SelfEndTh):F1}°C exceeds {max:0.##}°C (ntthal)" : null,
                    hairpin,
                    hairpin ? $"Hairpin Tm {v.HairpinTh:F1}°C exceeds {max:0.##}°C (ntthal)" : null);
            }
        }

        var (selfDimer, selfDimerWarning) = AlignmentSelfDimerScreen(sequence, param);
        // The stem flag is used only with AvoidSecondaryStructure (FinishProbe); skip the O(n³) scan otherwise.
        bool stem = param.AvoidSecondaryStructure && HasSecondaryStructurePotential(sequence);
        return (
            selfDimer,
            selfDimerWarning,
            stem,
            stem ? "Potential secondary structure" : null);
    }

    // Fallback self-dimer criterion: Primer3 alignment-mode internal-oligo oligo_compl (dpal self_any, self_end;
    // PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END, default 12.00), no length limit. The decision only needs the
    // threshold comparison (early-exit dpal), so the warning names the exceeded limit; ValidateProbe reports the
    // exact values.
    private static (bool Flag, string? Warning) AlignmentSelfDimerScreen(string sequence, ProbeParameters param)
    {
        if (!ThermodynamicPrimerDesigner.ExceedsPrimer3SelfComplementarity(sequence, param.MaxSelfAny, param.MaxSelfEnd, out bool any))
            return (false, null);
        return (true, any
            ? string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Self-complementarity: Primer3 self_any exceeds {param.MaxSelfAny:0.00}")
            : string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Self-complementarity: Primer3 self_end exceeds {param.MaxSelfEnd:0.00}"));
    }

    private static (bool Flag, string? Warning) AlignmentSelfDimerScreen(double selfAny, double selfEnd, ProbeParameters param)
    {
        if (selfAny > param.MaxSelfAny)
            return (true, string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"Self-complementarity: Primer3 self_any {selfAny:0.00} exceeds {param.MaxSelfAny:0.00}"));
        if (selfEnd > param.MaxSelfEnd)
            return (true, string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"Self-complementarity: Primer3 self_end {selfEnd:0.00} exceeds {param.MaxSelfEnd:0.00}"));
        return (false, null);
    }

    // Primer3 thermodynamic self-structure of a probe (ntthal ANY / END1 self-dimer and hairpin Tm at the
    // parameters' conditions), or null when the fallback screen applies (Heuristic screen, longer than
    // ThermodynamicScreenMaxLength = THAL_MAX_ALIGN (default 60), or a non-ACGT base).
    private static PrimerDesigner.Primer3OligoStructure? ComputeThermodynamicSelfStructure(
        string sequence, ProbeParameters param)
    {
        if (param.StructureScreen != ProbeStructureScreen.Thermodynamic || sequence.Length > param.ThermodynamicScreenMaxLength)
            return null;

        return ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure(
            sequence, param.MonovalentMillimolar, param.DivalentMillimolar, param.DntpMillimolar,
            param.DnaConcentrationNanomolar, param.ThermodynamicScreenMaxLength);
    }

    /// <summary>
    /// The top <paramref name="maxProbes"/> of <see cref="EnumerateRankedProbes"/>: every window of the configured
    /// length range (prefix-sum GC) is ranked; the self-structure screens only lower a score, so they are evaluated
    /// lazily and the scan stops once no remaining candidate can enter the top <paramref name="maxProbes"/>; the
    /// result equals an exhaustive evaluation: probes with score &gt; 0, ordered by score descending, ties by
    /// (length, start) ascending (or by the Primer3 penalty with <see cref="ProbeRanking.Primer3Penalty"/>).
    /// </summary>
    private static List<Probe> DesignProbesOptimized(
        string targetSequence,
        ProbeParameters param,
        int maxProbes) =>
        maxProbes <= 0
            ? new List<Probe>()
            : EnumerateRankedProbes(targetSequence, param).Take(maxProbes).ToList();

    /// <summary>
    /// Lazily yields every candidate (score &gt; 0 after the self-structure screens, then <paramref name="adjust"/>;
    /// a null from <paramref name="adjust"/> drops the probe) in the ranking order of <see cref="ProbeParameters.Ranking"/>:
    /// <see cref="ProbeRanking.AdditiveScore"/> — (adjusted) score descending, ties by enumeration order (length, start)
    /// ascending; <see cref="ProbeRanking.Primer3Penalty"/> — Primer3 <c>primer_rec_comp</c> on the penalty (which
    /// <paramref name="adjust"/> does not change). <paramref name="adjust"/> must never raise a score: the additive order is
    /// then produced by a lazy merge — candidates are finished in descending base-score order (an upper bound of their final
    /// score) and a finished probe is released once it ranks before the bound of every unfinished candidate — so a prefix
    /// of the output costs only the candidates it needs, and the full output equals an exhaustive evaluation + stable sort.
    /// </summary>
    private static IEnumerable<Probe> EnumerateRankedProbes(
        string targetSequence,
        ProbeParameters param,
        Func<Probe, Probe?>? adjust = null)
    {
        var bases = EvaluateCandidateBases(targetSequence, param);
        var structureCache = new Dictionary<string, (bool, string?, bool, string?)>(StringComparer.Ordinal);

        if (param.Ranking == ProbeRanking.Primer3Penalty)
        {
            // The penalty does not depend on the self-structure screens: sort by primer_rec_comp (a non-computable
            // penalty ranks last) and finish lazily.
            bases.Sort((a, b) => ComparePrimer3Rank(
                a.Primer3Penalty ?? double.PositiveInfinity, a.Start, a.Sequence.Length,
                b.Primer3Penalty ?? double.PositiveInfinity, b.Start, b.Sequence.Length));
            foreach (var b in bases)
            {
                if (FinishProbe(b, param, structureCache) is { } p && (adjust is null ? p : adjust(p)) is { } q)
                    yield return q;
            }
            yield break;
        }

        // Rank key: score descending, then enumeration order ascending (= stable sort of the eager scan).
        static int CompareKeys((double Score, int Order) a, (double Score, int Order) b)
        {
            int c = b.Score.CompareTo(a.Score);
            return c != 0 ? c : a.Order.CompareTo(b.Order);
        }

        var pending = new PriorityQueue<Probe, (double Score, int Order)>(
            Comparer<(double Score, int Order)>.Create(CompareKeys));
        // Stable: equal base scores keep the enumeration order, so the bound (BaseScore, Order) is non-decreasing in rank.
        foreach (var b in bases.OrderByDescending(x => x.BaseScore))
        {
            // Structure screens and adjust never raise a score: a pending probe that ranks before this candidate's best
            // possible key ranks before every later candidate too.
            while (pending.TryPeek(out var head, out var key) && CompareKeys(key, (b.BaseScore, b.Order)) < 0)
            {
                pending.Dequeue();
                yield return head;
            }

            if (FinishProbe(b, param, structureCache) is { } p && (adjust is null ? p : adjust(p)) is { } q)
                pending.Enqueue(q, (q.Score, b.Order));
        }

        while (pending.TryDequeue(out var rest, out _))
            yield return rest;
    }

    // Every window of the configured length range (prefix-sum GC) with base score > 0, in enumeration order
    // (length ascending, then start ascending).
    private static List<ProbeBase> EvaluateCandidateBases(string targetSequence, ProbeParameters param)
    {
        int n = targetSequence.Length;

        // gcPrefixSum[i] / validPrefixSum[i] = count of G/C / of valid nucleotides (A/C/G/T/U) in sequence[0..i-1],
        // classified by the canonical SequenceExtensions.CountGcAndValidNucleotides, so a window's GC fraction equals
        // CalculateGcFractionFast(window) (= Primer3 gc_and_n_content: G+C over the non-N bases) in O(1).
        int[] gcPrefixSum = new int[n + 1];
        int[] validPrefixSum = new int[n + 1];
        for (int i = 0; i < n; i++)
        {
            var (gcBase, validBase) = targetSequence.AsSpan(i, 1).CountGcAndValidNucleotides();
            gcPrefixSum[i + 1] = gcPrefixSum[i] + gcBase;
            validPrefixSum[i + 1] = validPrefixSum[i] + validBase;
        }

        var bases = new List<ProbeBase>();
        int order = 0;
        for (int length = param.MinLength; length <= param.MaxLength && length <= n; length++)
        {
            for (int start = 0; start <= n - length; start++)
            {
                int gcCount = gcPrefixSum[start + length] - gcPrefixSum[start];
                int validCount = validPrefixSum[start + length] - validPrefixSum[start];
                double gc = validCount == 0 ? 0 : (double)gcCount / validCount;

                // Early rejection far outside the GC window.
                if (gc < param.MinGc - 0.1 || gc > param.MaxGc + 0.1)
                    continue;

                var b = EvaluateProbeBase(targetSequence.Substring(start, length), start, order++, param, gc);
                if (b.BaseScore > 0)
                    bases.Add(b);
            }
        }
        return bases;
    }

    /// <summary>
    /// Designs tiling probes that cover the whole target: windows every <c>probeLength − overlap</c> bases plus,
    /// when the last of them stops short of the 3′ end, one window anchored at the target end.
    /// </summary>
    /// <remarks>
    /// Windows of <paramref name="probeLength"/> start at <c>0, step, 2·step, …</c> (<c>step = probeLength − overlap</c>)
    /// while they fit; if the last one ends before the target end, a final window starts at
    /// <c>targetSequence.Length − probeLength</c> so the 3′ tail is covered — the end-anchored probe of CATCH
    /// (Metsky et al. 2019, <c>candidate_probes.make_candidate_probes_from_sequence</c>: "There are bases on the right
    /// that were never covered, so add another probe for this"). It is added exactly when
    /// <c>(Length − probeLength) mod step ≠ 0</c> (CATCH's own test, <c>len(seq) mod stride ≠ 0</c>, can duplicate the
    /// last window or miss a tail, e.g. 120 nt / 60 / stride 40); it overlaps its predecessor by more than
    /// <paramref name="overlap"/>. With <paramref name="overlap"/> ≥ 0 every base is covered
    /// (<see cref="TilingProbeSet.Coverage"/> = target length); a negative overlap leaves gaps between windows and
    /// <see cref="TilingProbeSet.Coverage"/> counts the positions actually covered. Each window is scored
    /// like <see cref="DesignProbes(string, ProbeParameters?, int)"/> (Tm = Primer3 <c>seqtm</c> at the
    /// parameters' conditions, thermodynamic self-structure screen for ≤ ThermodynamicScreenMaxLength nt, default 60). A window whose score is
    /// ≤ 0 is still emitted (score 0.3, warning "Suboptimal probe…") so coverage is preserved. Windows keep their
    /// tiling order (start ascending); with <see cref="ProbeParameters.Ranking"/> = <see cref="ProbeRanking.Primer3Penalty"/> each scored
    /// window also carries <see cref="Probe.Primer3Penalty"/>.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="targetSequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="probeLength"/> ≤ 0 or longer than the target
    /// (CATCH rejects a sequence shorter than the probe length), or <paramref name="overlap"/> ≥
    /// <paramref name="probeLength"/> (the step would not advance); invalid <paramref name="parameters"/> (Primer3-ranking optima, ThermodynamicScreenMaxLength outside 60–10 000, a negative MaxNearestNeighborLength or illegal reaction conditions, as for <see cref="DesignProbes(string, ProbeParameters?, int)"/>; its MinLength / MaxLength are not used).</exception>
    public static TilingProbeSet DesignTilingProbes(
        string targetSequence,
        int probeLength = 60,
        int overlap = 20,
        ProbeParameters? parameters = null)
    {
        ArgumentNullException.ThrowIfNull(targetSequence);
        if (probeLength <= 0)
            throw new ArgumentOutOfRangeException(nameof(probeLength), probeLength, "Probe length must be positive.");
        if (overlap >= probeLength)
            throw new ArgumentOutOfRangeException(nameof(overlap), overlap,
                $"Overlap must be less than the probe length ({probeLength}) so the tiling step is positive.");
        if (probeLength > targetSequence.Length)
            throw new ArgumentOutOfRangeException(nameof(probeLength), probeLength,
                $"Probe length must not exceed the target length ({targetSequence.Length}).");

        var param = parameters ?? Defaults.Microarray with
        {
            MinLength = probeLength,
            MaxLength = probeLength
        };
        ValidateRanking(param, nameof(parameters));
        ValidateDesignParameters(param, nameof(parameters), checkLengths: false);

        targetSequence = targetSequence.ToUpperInvariant();
        var probes = new List<Probe>();
        // long: a very negative overlap (gapped tiling) must not overflow the step (audit round 7, A7-3).
        long step = (long)probeLength - overlap;
        int lastStart = targetSequence.Length - probeLength;

        for (long start = 0; start <= lastStart; start += step)
            probes.Add(TilingWindow(targetSequence, (int)start, probeLength, param));

        // End-anchored window (CATCH) when the regular grid stops short of the 3' end.
        if (lastStart % step != 0)
            probes.Add(TilingWindow(targetSequence, lastStart, probeLength, param));

        // Covered positions: union of the windows (sorted by start).
        int covered = 0, coveredEnd = 0;
        foreach (var probe in probes)
        {
            covered += Math.Max(0, probe.End + 1 - Math.Max(probe.Start, coveredEnd));
            coveredEnd = Math.Max(coveredEnd, probe.End + 1);
        }

        double meanTm = probes.Average(p => p.Tm);
        double tmRange = probes.Max(p => p.Tm) - probes.Min(p => p.Tm);

        return new TilingProbeSet(probes, covered, meanTm, tmRange);
    }

    private static Probe TilingWindow(string targetSequence, int start, int probeLength, ProbeParameters param)
    {
        string probeSeq = targetSequence.Substring(start, probeLength);
        // One evaluation (all screens, eager); a window scoring ≤ 0 is still emitted for coverage with the same Tm
        // (0 when not computable) and GC, instead of recomputing them.
        var b = EvaluateProbeBase(probeSeq, start, 0, param, probeSeq.CalculateGcFractionFast());
        if (FinishProbe(b, param) is { } probe)
            return probe with { Type = ProbeType.Tiling };

        return new Probe(
            probeSeq, start, start + probeLength - 1,
            b.Tm, b.Gc, 0.3, ProbeType.Tiling,
            new List<string> { "Suboptimal probe, included for coverage" });
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
    /// Validates a probe: ungapped k-mismatch hit count over the references, self-structure screen, and
    /// (optionally) the Kane et al. (2000) cross-hybridization criteria against known non-target sequences.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Hits.</b> Every reference is scanned with the canonical ungapped k-mismatch (Hamming) matcher
    /// <see cref="ApproximateMatcher.FindWithMismatches(string, string, int)"/> (case-insensitive, overlapping; the
    /// strand given only by default, and with <paramref name="bothStrands"/> also for the probe's reverse complement —
    /// the probe's sites on the other strand of a double-stranded reference, as blastn <c>-strand both</c> and
    /// <see cref="CheckSpecificity(string, global::SuffixTree.ISuffixTree, bool)"/> with <c>bothStrands</c>; a reference
    /// position hit in both orientations — a reverse-palindromic probe, or a reverse-palindromic site — is one site,
    /// counted once, as <see cref="CheckSpecificity(string, global::SuffixTree.ISuffixTree, bool)"/> counts a
    /// palindromic probe once); <see cref="ProbeValidation.OffTargetHits"/> is the total and includes the intended site.
    /// Each hit is judged by the Kane et al. (2000) criteria on its ungapped diagonal — identity (L − mismatches) / L
    /// over the probe length &gt; <paramref name="maxNonTargetIdentity"/> or a run of identical positions longer than
    /// <paramref name="maxContiguousMatch"/> (<see cref="ProbeValidation.CrossHybridizingHits"/>); more than one such
    /// site records an issue (the probe can bind somewhere besides its intended site). With the default 3 mismatches
    /// every hit of a probe ≥ 13 nt meets the identity criterion ((L − 3) / L &gt; 0.75), so the decision equals the
    /// former "more than one hit" rule there. <see cref="ProbeValidation.SpecificityScore"/> = 1/N for N ≥ 1 hits
    /// (0 for none) is a library-defined uniqueness score (the share of the probe's N candidate binding sites
    /// taken by one site), not a published specificity metric (no published 1/N score was found — Kane 2000,
    /// OligoArray 2.0, Li &amp; Stormo 2001 and Primer3's PRIMER_INTERNAL_MAX_LIBRARY_MISHYB judge a site by
    /// identity / free energy / alignment score, not by a hit count); it is reported only and does not enter
    /// <see cref="ProbeValidation.IsValid"/>.
    /// </para>
    /// <para>
    /// <b>Self-structure.</b> The same screen as <see cref="DesignProbes(string, ProbeParameters?, int)"/>: with
    /// <see cref="ProbeStructureScreen.Thermodynamic"/> (default) an A/C/G/T probe of ≤ <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> (default 60) nt is screened as Primer3
    /// screens a hybridization probe — ntthal self-dimer (ANY), 3′ self-dimer (END1) and hairpin Tm
    /// (<see cref="ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure"/>, primer3-py <c>calc_homodimer</c> /
    /// <c>calc_end_stability</c> / <c>calc_hairpin</c> parity) at the <paramref name="conditions"/> salt/oligo
    /// concentrations must not exceed <see cref="ProbeParameters.MaxStructureTm"/> (PRIMER_INTERNAL_MAX_SELF_ANY_TH
    /// = _SELF_END_TH = _HAIRPIN_TH = 47 °C). Probes longer than <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/>
    /// (thal.c THAL_MAX_ALIGN, default 60; opt-in larger values run ntthal on longer probes), non-ACGT probes and
    /// <see cref="ProbeStructureScreen.Heuristic"/> use the fallback screens: Primer3 alignment-mode internal-oligo
    /// self-complementarity (dpal <c>self_any</c> / <c>self_end</c> &gt; <see cref="ProbeParameters.MaxSelfAny"/> /
    /// <see cref="ProbeParameters.MaxSelfEnd"/>, PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END = 12.00, no length limit)
    /// and the sequence-only stem-loop screen <see cref="PrimerDesigner.HasHairpinPotential"/> (an exactly complementary
    /// stem of ≥ 4 bp — library convention — closing a loop of ≥ 3 nt, thal.c <c>min_hrpn_loop</c>; almost every random
    /// probe of &gt; 60 nt is flagged, so raise <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> for a hairpin
    /// Tm decision on long probes). The fallback stem-loop flag is reported in
    /// <see cref="ProbeValidation.HasSecondaryStructure"/> and <see cref="ProbeValidation.Warnings"/> but is not an issue:
    /// validity rests on sourced criteria only (audit round 3, A3-12), and a stem ≥ 4 bp is a library convention
    /// (audit round 5, A5-5 — a correction of the validity rule; before it, this flag made almost every &gt; 60-nt probe
    /// invalid). The ntthal hairpin Tm limit (PRIMER_INTERNAL_MAX_HAIRPIN_TH) and the fallback alignment-mode self_any /
    /// self_end limits (PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END) are issues. The alignment-mode
    /// values are reported for every probe in <see cref="ProbeValidation.SelfAny"/> / <see cref="ProbeValidation.SelfEnd"/>;
    /// <see cref="ProbeValidation.SelfComplementarity"/> (fold-back fraction) is a library metric only.
    /// </para>
    /// <para>
    /// <b>Cross-hybridization.</b> When <paramref name="nonTargetSequences"/> is given, every non-target (both
    /// strands) is assessed with <see cref="AssessCrossHybridization"/> (Kane et al. 2000: overall identity &gt; 75 %
    /// or a contiguous identical stretch &gt; 15 nt → the probe may cross-hybridize); each cross-hybridizing
    /// non-target strand records an issue.
    /// </para>
    /// <para><see cref="ProbeValidation.IsValid"/> is true when no issue was recorded
    /// (<see cref="ProbeValidation.Warnings"/> do not count).</para>
    /// </remarks>
    /// <param name="probeSequence">Probe sequence to validate (case-insensitive). Null throws; empty → invalid result.</param>
    /// <param name="referenceSequences">Reference sequences scanned for ungapped hits (target included).</param>
    /// <param name="maxMismatches">Search radius of the ungapped site scan (≥ 0; default 3 — a library convention kept for
    /// compatibility, no published hybridization source: the figure is the CRISPR guide mismatch tolerance). Whether a
    /// found site is a problem is decided by the Kane criteria; to find every ungapped site that meets the Kane identity
    /// criterion pass ⌈L/4⌉ − 1 (the largest d with (L − d) / L &gt; 0.75).</param>
    /// <param name="selfComplementarityThreshold">Former fold-back-fraction limit of the fallback screen; kept for
    /// source compatibility and copied into <see cref="ProbeParameters.MaxSelfComplementarity"/>, but no longer used by
    /// any screen (the fallback self-dimer limits are <see cref="ProbeParameters.MaxSelfAny"/> / <see cref="ProbeParameters.MaxSelfEnd"/>)
    /// (default 0.3, the Microarray preset's <see cref="ProbeParameters.MaxSelfComplementarity"/>).</param>
    /// <param name="conditions">Hybridization conditions and structure-screen settings (default: Primer3 probe
    /// conditions 50 nM / 50 mM / 0 Mg²⁺ / 0 dNTP, thermodynamic screen, 47 °C — the settings of
    /// <see cref="Defaults.Microarray"/> at those conditions); its <see cref="ProbeParameters.MaxSelfComplementarity"/> is replaced by
    /// <paramref name="selfComplementarityThreshold"/>. Used: the salt / dNTP / oligo concentrations (ntthal screen and
    /// non-target duplex Tm; Primer3 <c>_pr_data_control</c> legality — monovalent and oligo &gt; 0, Mg²⁺ and dNTP ≥ 0,
    /// else <see cref="ArgumentOutOfRangeException"/>), <see cref="ProbeParameters.StructureScreen"/>,
    /// <see cref="ProbeParameters.MaxStructureTm"/>, <see cref="ProbeParameters.MaxSelfAny"/> / <see cref="ProbeParameters.MaxSelfEnd"/>,
    /// <see cref="ProbeParameters.ThermodynamicScreenMaxLength"/> (THAL_MAX_ALIGN of the screen and of the non-target duplex
    /// Tm; 60–10 000, else <see cref="ArgumentOutOfRangeException"/>).</param>
    /// <param name="nonTargetSequences">Optional known non-target sequences for the Kane assessment.</param>
    /// <param name="maxNonTargetIdentity">Kane identity threshold in [0, 1] (default 0.75; flagged when strictly above),
    /// for the reference sites and the non-targets; outside [0, 1] or NaN → <see cref="ArgumentOutOfRangeException"/>.</param>
    /// <param name="maxContiguousMatch">Kane contiguous-identity threshold in nt (≥ 0, default 15; flagged when strictly
    /// longer), for the reference sites and the non-targets; negative → <see cref="ArgumentOutOfRangeException"/>.</param>
    /// <param name="maxDuplexTm">Optional OligoArray-style off-target duplex-Tm threshold (°C; see
    /// <see cref="AssessCrossHybridization"/>); null = Kane criteria only.</param>
    /// <param name="bothStrands">Also scan the references for the probe's reverse complement (default false: the
    /// references' given strand only, the behaviour before audit round 4, B07 F60). A site found in both orientations
    /// at the same reference position is counted once and is cross-hybridizing when either orientation meets a Kane
    /// criterion. The non-target assessment is two-stranded regardless (<see cref="AssessCrossHybridization"/>).</param>
    public static ProbeValidation ValidateProbe(
        string probeSequence,
        IEnumerable<string> referenceSequences,
        int maxMismatches = 3,
        double selfComplementarityThreshold = 0.3,
        ProbeParameters? conditions = null,
        IEnumerable<string>? nonTargetSequences = null,
        double maxNonTargetIdentity = KaneMaxIdentity,
        int maxContiguousMatch = KaneMaxContiguousMatch,
        double? maxDuplexTm = null,
        bool bothStrands = false)
    {
        ArgumentNullException.ThrowIfNull(probeSequence);
        ArgumentNullException.ThrowIfNull(referenceSequences);
        if (double.IsNaN(maxNonTargetIdentity) || maxNonTargetIdentity < 0 || maxNonTargetIdentity > 1)
            throw new ArgumentOutOfRangeException(nameof(maxNonTargetIdentity), "Identity threshold must be in [0, 1].");
        if (maxContiguousMatch < 0)
            throw new ArgumentOutOfRangeException(nameof(maxContiguousMatch), "Contiguous-match threshold cannot be negative.");
        if (conditions is { } stated)
        {
            PrimerDesigner.ValidatePrimer3Conditions(stated.MonovalentMillimolar, stated.DivalentMillimolar,
                stated.DntpMillimolar, stated.DnaConcentrationNanomolar, nameof(conditions));
            ValidateThermodynamicScreenMaxLength(stated, nameof(conditions));
        }

        probeSequence = probeSequence.ToUpperInvariant();
        var issues = new List<string>();
        var warnings = new List<string>();

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
        int crossHybridizingHits = 0;

        // Canonical ungapped k-mismatch (Hamming) scan; case-insensitive, overlapping hits included. Each hit is a
        // candidate binding site; it counts as a cross-hybridizing site when it meets a Kane et al. (2000) criterion.
        // With bothStrands the probe's reverse complement is scanned too (the probe's sites on the other strand of a
        // double-stranded reference, as blastn -strand both / CheckSpecificity(bothStrands)); a reference position hit
        // in both orientations (a reverse-palindromic probe, or a site that is itself a reverse palindrome) is one
        // site, cross-hybridizing when either orientation meets a Kane criterion (audit round 4, A4-2, F60).
        string[] patterns = [probeSequence];
        if (bothStrands)
        {
            string reverseComplement = DnaSequence.GetReverseComplementString(probeSequence);
            if (!string.Equals(reverseComplement, probeSequence, StringComparison.Ordinal))
                patterns = [probeSequence, reverseComplement];
        }
        foreach (var reference in referenceSequences)
        {
            var sites = new Dictionary<int, bool>();
            foreach (string pattern in patterns)
            {
                foreach (var hit in ApproximateMatcher.FindWithMismatches(reference, pattern, maxMismatches))
                {
                    bool kane = UngappedSiteMeetsKaneCriteria(hit.MismatchPositions, probeSequence.Length, maxNonTargetIdentity, maxContiguousMatch);
                    sites[hit.Position] = sites.TryGetValue(hit.Position, out bool seen) ? seen || kane : kane;
                }
            }
            offTargetHits += sites.Count;
            crossHybridizingHits += sites.Values.Count(kane => kane);
        }

        // More than one site meeting the Kane criteria = the probe can bind somewhere besides its intended site.
        if (crossHybridizingHits > 1)
        {
            issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                $"{crossHybridizingHits} potential off-target sites (Kane 2000: identity > {maxNonTargetIdentity * 100:0.##}% or > {maxContiguousMatch} contiguous identical nt)"));
        }

        // Self-structure: the DesignProbes screen (Primer3 ntthal for ACGT probes ≤ ThermodynamicScreenMaxLength; otherwise Primer3
        // alignment-mode self_any / self_end + the PrimerDesigner.HasHairpinPotential stem-loop screen).
        var param = (conditions ?? Primer3ProbeConditions) with { MaxSelfComplementarity = selfComplementarityThreshold };
        double selfComp = CalculateSelfComplementarity(probeSequence);
        double alnSelfAny = ThermodynamicPrimerDesigner.CalculatePrimerSelfAnyComplementarity(probeSequence);
        double alnSelfEnd = PrimerDesigner.CalculatePrimerSelfEndComplementarity(probeSequence);
        var thermo = ComputeThermodynamicSelfStructure(probeSequence, param);
        bool selfCompIssue;
        bool hasStructure;
        if (thermo is { } t)
        {
            double max = param.MaxStructureTm;
            selfCompIssue = t.SelfAnyTh > max || t.SelfEndTh > max;
            hasStructure = t.HairpinTh > max;
            if (selfCompIssue)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Self-complementarity: ntthal self-dimer Tm {Math.Max(t.SelfAnyTh, t.SelfEndTh):F1}°C exceeds {max:0.##}°C"));
            if (hasStructure)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"Potential secondary structure formation: ntthal hairpin Tm {t.HairpinTh:F1}°C exceeds {max:0.##}°C"));
        }
        else
        {
            (selfCompIssue, string? selfCompWarning) = AlignmentSelfDimerScreen(alnSelfAny, alnSelfEnd, param);
            hasStructure = HasSecondaryStructurePotential(probeSequence);
            if (selfCompIssue)
                issues.Add(selfCompWarning!);
            // The sequence-only stem-loop screen (stem >= 4 bp: library convention) is a presence test that flags almost
            // every random probe of > 60 nt; per the validity rule it is reported (HasSecondaryStructure, Warnings) but
            // does not decide IsValid (audit round 5, A5-5). The ntthal hairpin Tm above (Primer3
            // PRIMER_INTERNAL_MAX_HAIRPIN_TH) decides validity wherever it runs.
            if (hasStructure)
                warnings.Add("Potential secondary structure formation");
        }

        // Kane et al. (2000) cross-hybridization criteria against known non-targets (optional).
        IReadOnlyList<CrossHybridizationAssessment> cross = Array.Empty<CrossHybridizationAssessment>();
        if (nonTargetSequences is not null)
        {
            cross = AssessCrossHybridization(probeSequence, nonTargetSequences, maxNonTargetIdentity, maxContiguousMatch,
                conditions: param, maxDuplexTm: maxDuplexTm);
            foreach (var c in cross.Where(c => c.CrossHybridizes))
            {
                issues.Add($"Cross-hybridization risk with non-target {c.NonTargetIndex}"
                    + (c.ReverseComplementStrand ? " (reverse complement)" : "")
                    + string.Create(System.Globalization.CultureInfo.InvariantCulture, $": identity {c.Identity * 100:F0}%, longest contiguous match {c.LongestContiguousMatch} nt (Kane 2000)")
                    + (c.ExceedsDuplexTmThreshold
                        ? string.Create(System.Globalization.CultureInfo.InvariantCulture, $", site duplex Tm {c.DuplexTm:F1}°C > {maxDuplexTm:0.##}°C")
                        : ""));
            }
        }

        // Library uniqueness score: 0 hits → 0, N ≥ 1 hits → 1/N.
        double specificity = offTargetHits == 0 ? 0.0 : 1.0 / offTargetHits;

        return new ProbeValidation(
            issues.Count == 0,
            specificity,
            offTargetHits,
            selfComp,
            hasStructure,
            issues)
        {
            ThermodynamicScreen = thermo is not null,
            SelfDimerTm = thermo?.SelfAnyTh,
            SelfEndDimerTm = thermo?.SelfEndTh,
            HairpinTm = thermo?.HairpinTh,
            SelfAny = alnSelfAny,
            SelfEnd = alnSelfEnd,
            CrossHybridization = cross,
            CrossHybridizingHits = crossHybridizingHits,
            Warnings = warnings,
        };
    }

    // Probe Tm = Primer3 seqtm at the parameters' hybridization conditions (NaN when not computable).
    private static double CalculateProbeTm(string sequence, ProbeParameters param) =>
        ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(
            sequence,
            param.DnaConcentrationNanomolar,
            param.MonovalentMillimolar,
            param.DivalentMillimolar,
            param.DntpMillimolar,
            param.MaxNearestNeighborLength);

    // Primer3 hybridization-probe (internal-oligo) conditions 50 nM / 50 mM / 0 Mg²⁺ / 0 dNTP, nearest-neighbour Tm up to
    // 36 nt, with the remaining settings of the Microarray preset (thermodynamic screen, 47 °C): the default conditions of
    // ValidateProbe, AssessCrossHybridization, DesignMolecularBeacon and AnalyzeOligo (the Microarray preset's own
    // conditions are OligoArray's 1 M / 1 µM since A3-10).
    private static ProbeParameters Primer3ProbeConditions => Defaults.Microarray with
    {
        DnaConcentrationNanomolar = PrimerDesigner.Primer3InternalDnaConcentrationNanomolar,
        MonovalentMillimolar = PrimerDesigner.Primer3InternalMonovalentMillimolar,
        DivalentMillimolar = PrimerDesigner.Primer3InternalDivalentMillimolar,
        DntpMillimolar = PrimerDesigner.Primer3InternalDntpMillimolar,
        MaxNearestNeighborLength = PrimerDesigner.Primer3MaxNnTmLength
    };

    // Eager ProbeParameters checks of the designers (heavy tier, HEAVY-2): MinLength < 1 crashed the window scan
    // (IndexOutOfRange) and illegal conditions surfaced lazily from the Tm with a foreign parameter name. Primer3
    // _pr_data_control: PRIMER_INTERNAL_MIN_SIZE must be >= 1 (checked only where MinLength is used), PRIMER_INTERNAL_SALT_
    // MONOVALENT / DNA_CONC > 0, SALT_DIVALENT / DNTP_CONC >= 0 (canonical ValidatePrimer3Conditions); the Tm's
    // max_nn_length must be >= 0; plus the ranking optima and THAL_MAX_ALIGN range.
    private static void ValidateDesignParameters(ProbeParameters param, string paramName, bool checkLengths)
    {
        if (checkLengths && param.MinLength < 1)
            throw new ArgumentOutOfRangeException(paramName,
                "PRIMER_INTERNAL_MIN_SIZE must be >= 1 (Primer3 _pr_data_control): MinLength must be at least 1.");
        if (param.MaxNearestNeighborLength < 0)
            throw new ArgumentOutOfRangeException(paramName, "MaxNearestNeighborLength (max_nn_length) cannot be negative.");
        PrimerDesigner.ValidatePrimer3Conditions(param.MonovalentMillimolar, param.DivalentMillimolar,
            param.DntpMillimolar, param.DnaConcentrationNanomolar, paramName);
        ValidateThermodynamicScreenMaxLength(param, paramName);
    }

    // Cheap per-candidate properties (everything except the self-structure screens). The ranking score
    // is rebuilt from these in the original penalty order, so lazily adding the structure penalties gives
    // bit-identical scores to an eager evaluation.
    private readonly record struct ProbeBase(
        string Sequence,
        int Start,
        int Order,
        double Gc,
        double Tm,
        bool TmComputable,
        bool GcOutside,
        bool TmOutside,
        int Homopolymer,
        bool HasRepeats,
        double PositionPenalty,
        double BaseScore,
        double? Primer3Penalty);

    private static PrimerDesigner.Primer3OligoStructure ComputePrimer3ProbeStructure(string oligo, Primer3ProbeSettings s) =>
        ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure(
            oligo, s.MonovalentMillimolar, s.DivalentMillimolar, s.DntpMillimolar, s.DnaConcentrationNanomolar)!.Value;

    // Tyagi & Kramer (1996) / Marras et al. design rule: the probe–target hybrid Tm and the stem
    // (hairpin) Tm should both be 7–10 °C above the detection (PCR annealing) temperature.
    private const double BeaconMinTmAboveDetection = 7.0;

    private const double BeaconMaxTmAboveDetection = 10.0;

    // Kane et al. (2000) criteria applied to an ungapped k-mismatch site of a probe of length L with mismatches at the
    // given (ascending, probe-relative) positions: identity (L − d) / L over the probe length > maxIdentity, or the
    // longest run of identical positions > maxContiguousMatch (the same measures AssessCrossHybridization takes from the
    // best local alignment / longest common substring, here on the ungapped diagonal of the hit).
    private static bool UngappedSiteMeetsKaneCriteria(
        IReadOnlyList<int> mismatchPositions, int probeLength, double maxIdentity, int maxContiguousMatch)
    {
        double identity = (double)(probeLength - mismatchPositions.Count) / probeLength;
        if (identity > maxIdentity)
            return true;
        int longest = 0, previous = -1;
        foreach (int position in mismatchPositions)
        {
            longest = Math.Max(longest, position - previous - 1);
            previous = position;
        }
        longest = Math.Max(longest, probeLength - previous - 1);
        return longest > maxContiguousMatch;
    }

    // --- Kane et al. (2000) cross-hybridization criteria ---
    // Kane MD et al. (2000) Nucleic Acids Res 28(22):4552-4557 (abstract): non-target transcripts ">75% similar
    // over the 50 base target may show cross-hybridization", and a non-target region must not include a stretch
    // of identical sequence ">15 contiguous bases"; summarised by later probe-design pipelines as "a probe is
    // likely to cross-hybridize with a nontarget if overall sequence identity is > 75% or if there is a
    // contiguous match > 15 bp".
    private const double KaneMaxIdentity = 0.75;

    private const int KaneMaxContiguousMatch = 15;

    /// <summary>
    /// Analyzes oligonucleotide properties: Tm (Primer3 <c>seqtm</c> at the Primer3 hybridization-probe
    /// conditions 50 nM / 50 mM monovalent / no Mg²⁺ / no dNTP — SantaLucia 1998 nearest-neighbour for
    /// ≤ 36 nt, <c>long_seq_tm</c> above; <c>NaN</c> for fewer than 2 bases or any non-ACGT base, e.g. RNA),
    /// G+C fraction, single-stranded molecular weight (<see cref="CalculateMolecularWeight(string)"/>) and the
    /// mononucleotide-sum ε260 (<see cref="CalculateExtinctionCoefficient(string)"/>).
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static (double Tm, double GcContent, double MolecularWeight, double ExtinctionCoefficient)
        AnalyzeOligo(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        sequence = sequence.ToUpperInvariant();

        double tm = CalculateProbeTm(sequence, Primer3ProbeConditions);
        double gc = sequence.CalculateGcFractionFast();
        double mw = CalculateMolecularWeight(sequence);
        double extinction = CalculateExtinctionCoefficient(sequence);

        return (tm, gc, mw, extinction);
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

    // Fallback hairpin criterion (> ThermodynamicScreenMaxLength nt, non-ACGT, ProbeStructureScreen.Heuristic): the
    // canonical sequence-only stem-loop screen PrimerDesigner.HasHairpinPotential with its defaults — an exactly
    // complementary stem of ≥ 4 bp (library convention, F49) closing a loop of ≥ 3 nt (thal.c min_hrpn_loop = 3) —
    // in place of the former private scan that tested a loop of exactly 3 nt with an unsourced ≥ 80 % stem match
    // (audit round 5, A5-2).
    private static bool HasSecondaryStructurePotential(string sequence) => PrimerDesigner.HasHairpinPotential(sequence);

    // Simple-sequence repeat: a di- or trinucleotide microsatellite of ≥ 4 complete copies, found by the
    // canonical MISA-convention scanner (primitive units only, A/C/G/T only).
    private static bool HasSimpleRepeats(string sequence) =>
        Seqeron.Genomics.Analysis.RepeatFinder.FindMicrosatellites(sequence, 2, 3, 4).Any();

    /// <summary>
    /// Primer3 internal-oligo candidate list (<c>make_internal_oligo_list</c> → <c>pick_primer_range</c>
    /// → <c>calc_and_check_oligo_features</c>, <c>OT_INTL</c>, thermodynamic mode) over
    /// <c>seq[regionStart, regionEnd)</c>, in Primer3's enumeration order (3′ end descending, then length
    /// ascending), with <see cref="Primer3Probe.Penalty"/> = <c>p_obj_fn</c> (internal-oligo weights).
    /// With <paramref name="screenStructure"/> = true (PRIMER_TASK=pick_hyb_probe_only, a <c>primer_list</c>
    /// output) the ntthal self-any / self-end / hairpin limits are applied while enumerating, exactly as
    /// Primer3 does for list output (a too-stable self-dimer is a "five-prime problem" ending the extension);
    /// with false (internal oligo for a primer pair, <c>primer_pairs</c> output) Primer3 postpones them to
    /// <c>choose_internal_oligo</c>: the structure fields are then <see cref="double.NaN"/> and the caller
    /// screens the oligo it picks (<see cref="PassesPrimer3ProbeStructure"/>).
    /// </summary>
    internal static List<Primer3Probe> EnumeratePrimer3InternalOligos(
        string seq, int regionStart, int regionEnd, Primer3ProbeSettings s, bool screenStructure,
        IReadOnlyList<int>? quality = null, string? lowercase = null)
    {
        // Primer3 o_args.weights defaults, with PRIMER_INTERNAL_WT_GC_PERCENT_GT/_LT and PRIMER_INTERNAL_OPT_GC_PERCENT.
        // PRIMER_INTERNAL_WT_LIBRARY_MISHYB is the o_args repeat_sim weight.
        var weights = ThermodynamicPrimerDesigner.DefaultPrimer3Weights with
        {
            GcGt = s.WeightGcPercentGt,
            GcLt = s.WeightGcPercentLt,
            LibraryMispriming = s.WeightLibraryMishyb,
            BoundGt = s.WeightBoundGt,
            BoundLt = s.WeightBoundLt,
            SequenceQuality = s.WeightSequenceQuality,
        };
        // PRIMER_INTERNAL_MISHYB_LIBRARY: scored while enumerating for list output (three_conditions) or when weighted
        // (calc_and_check_oligo_features), otherwise postponed to choose_internal_oligo.
        var library = s.ActiveMishybLibrary;
        bool scoreLibrary = library is not null && (screenStructure || s.WeightLibraryMishyb != 0);
        var optima = new Primer3Optima(s.OptTm, s.OptSize, s.OptGcPercent)
        {
            OptBound = s.OptBound,
        };
        double monovalentEq = PrimerDesigner.Primer3MonovalentEquivalent(s.MonovalentMillimolar, s.DivalentMillimolar, s.DntpMillimolar);
        var accepted = new List<Primer3Probe>();

        // pick_primer_range: for every 3' end, oligos of increasing length (5' extensions). A failure that no
        // 5' extension can cure — too many Ns, poly-X, self-any (Primer3 five_prime_problem bits) — ends the
        // extension loop for that 3' end; checks run in calc_and_check_oligo_features order and stop at the
        // first failure (Ns, GC, poly-X, Tm, self-any, self-end, hairpin).
        for (int end = regionEnd - 1; end >= regionStart + s.MinSize - 1; end--)
        {
            // PRIMER_LOWERCASE_MASKING (is_lowercase_masked, before every other check): a lower-case 3' base rejects every
            // oligo ending here (lowercase = the case-preserving template, null when masking is off).
            if (lowercase is not null && PrimerDesigner.IsLowercaseMaskedBase(lowercase[end]))
                continue;
            for (int len = s.MinSize; len <= s.MaxSize; len++)
            {
                int start = end - len + 1;
                if (start < regionStart)
                    break;

                string oligo = seq.Substring(start, len);
                int gc = 0, ns = 0;
                bool other = false;
                foreach (char c in oligo)
                {
                    if (c is 'G' or 'C') gc++;
                    else if (c == 'N') ns++;
                    else if (c is not ('A' or 'T')) other = true;
                }
                if (ns > 0)
                    break; // OP_TOO_MANY_NS (PRIMER_INTERNAL_MAX_NS_ACCEPTED = 0): five-prime problem

                double gcPercent = 100.0 * gc / len;
                if (gcPercent < s.MinGcPercent || gcPercent > s.MaxGcPercent)
                    continue;
                // sequence_quality_is_ok (OT_INTL: the minimum quality only), after the GC check: OP_LOW_SEQUENCE_QUALITY
                // is a five-prime problem (every 5' extension keeps the low base).
                int? minQuality = null;
                if (quality is not null)
                {
                    minQuality = ThermodynamicPrimerDesigner.CalculateSequenceQualityPrimer3(quality, start, len, true, s.QualityRangeMax).Min;
                    if (minQuality < s.MinQuality)
                        break;
                }
                if (PrimerDesigner.FindLongestHomopolymer(oligo) > s.MaxPolyX)
                    break; // OP_HIGH_POLY_X: five-prime problem
                if (other)
                    continue; // non-ACGT, non-N symbol: no Tm (seqtm error → low Tm)

                var (tm, bound) = ThermodynamicPrimerDesigner.Primer3SeqTm(oligo, s.DnaConcentrationNanomolar, monovalentEq, s.AnnealingTemperature);
                if (tm < s.MinTm || tm > s.MaxTm)
                    continue;
                // Fraction bound (PRIMER_ANNEALING_TEMP > 0; not a five-prime problem).
                if (s.AnnealingTemperature > 0.0 && (bound < s.MinBound || bound > s.MaxBound))
                    continue;

                double selfAny = double.NaN, selfEnd = double.NaN, hairpin = double.NaN;
                double? alnAny = null, alnEnd = null;
                if (screenStructure && !s.ThermodynamicOligoAlignment)
                {
                    // oligo_compl (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0): self_any, then self_end; no hairpin.
                    alnAny = ThermodynamicPrimerDesigner.CalculatePrimerSelfAnyComplementarity(oligo);
                    if (alnAny > s.MaxSelfAny)
                        break; // OP_HIGH_SELF_ANY: five-prime problem
                    alnEnd = PrimerDesigner.CalculatePrimerSelfEndComplementarity(oligo);
                    if (alnEnd > s.MaxSelfEnd)
                        continue;
                }
                else if (screenStructure)
                {
                    var st = ComputePrimer3ProbeStructure(oligo, s);
                    if (st.SelfAnyTh > s.MaxSelfAnyTh)
                        break; // OP_HIGH_SELF_ANY: five-prime problem
                    if (st.SelfEndTh > s.MaxSelfEndTh || st.HairpinTh > s.MaxHairpinTh)
                        continue;
                    (selfAny, selfEnd, hairpin) = (st.SelfAnyTh, st.SelfEndTh, st.HairpinTh);
                }

                // oligo_repeat_library_mispriming (OT_INTL), after the structure checks.
                ThermodynamicPrimerDesigner.LibraryMispriming? lib = null;
                if (scoreLibrary)
                {
                    lib = ComputePrimer3ProbeLibrary(oligo, s, library!);
                    if (lib.Exceeds)
                        break; // OP_HIGH_SIM_TO_NON_TEMPLATE_SEQ: five-prime problem
                }

                // p_obj_fn OT_INTL: no end_stability term (end_oligodg is computed for primers only).
                double penalty = ThermodynamicPrimerDesigner.CalculatePrimer3Penalty(
                    new Primer3PenaltyInputs(tm, len, gcPercent)
                    {
                        LibraryMispriming = lib?.MaxScore ?? 0.0,
                        Bound = bound,
                        SequenceQuality = minQuality,
                        QualityRangeMax = s.QualityRangeMax,
                    }, weights, optima);
                accepted.Add(new Primer3Probe(oligo, start, len, tm, gcPercent, selfAny, selfEnd, hairpin, penalty)
                {
                    SelfAny = alnAny,
                    SelfEnd = alnEnd,
                    LibraryMishyb = lib?.MaxScore,
                    LibraryMishybName = lib?.Name,
                    Bound = s.AnnealingTemperature > 0.0 && bound != ThermodynamicPrimerDesigner.Primer3OligoTmError ? bound : null,
                    MinSequenceQuality = minQuality,
                });
            }
        }
        return accepted;
    }

    /// <summary>
    /// Primer3 <c>choose_internal_oligo</c> postponed checks of one internal oligo (<c>oligo_compl_thermod</c>:
    /// self-any and self-end Tm, then <c>oligo_hairpin</c>; alignment mode: <c>oligo_compl</c> dpal self_any /
    /// self_end; then, unless already scored while enumerating, <c>oligo_repeat_library_mispriming</c> against the
    /// mishyb library) at the settings' conditions; returns the oligo with its structure (and library) values filled
    /// in, or <c>null</c> when a limit is exceeded.
    /// </summary>
    internal static Primer3Probe? PassesPrimer3ProbeStructure(Primer3Probe probe, Primer3ProbeSettings s)
    {
        if (!s.ThermodynamicOligoAlignment)
        {
            // oligo_compl with o_args (alignment mode): self_any, then self_end.
            double any = ThermodynamicPrimerDesigner.CalculatePrimerSelfAnyComplementarity(probe.Sequence);
            if (any > s.MaxSelfAny)
                return null;
            double end = PrimerDesigner.CalculatePrimerSelfEndComplementarity(probe.Sequence);
            if (end > s.MaxSelfEnd)
                return null;
            probe = probe with { SelfAny = any, SelfEnd = end };
        }
        else
        {
            var st = ComputePrimer3ProbeStructure(probe.Sequence, s);
            if (st.SelfAnyTh > s.MaxSelfAnyTh || st.SelfEndTh > s.MaxSelfEndTh || st.HairpinTh > s.MaxHairpinTh)
                return null;
            probe = probe with { SelfAnyTh = st.SelfAnyTh, SelfEndTh = st.SelfEndTh, HairpinTh = st.HairpinTh };
        }
        if (s.ActiveMishybLibrary is { } library && probe.LibraryMishyb is null)
        {
            var lib = ComputePrimer3ProbeLibrary(probe.Sequence, s, library);
            if (lib.Exceeds)
                return null;
            probe = probe with { LibraryMishyb = lib.MaxScore, LibraryMishybName = lib.Name };
        }
        return probe;
    }

    /// <summary>
    /// Picks hybridization probes exactly as Primer3 does for <c>PRIMER_TASK=pick_hyb_probe_only</c>
    /// (Rozen &amp; Skaletsky 2000; Untergasser et al. 2012; <c>libprimer3.cc</c> <c>make_internal_oligo_list</c>
    /// → <c>pick_primer_range</c> → <c>calc_and_check_oligo_features</c> → <c>p_obj_fn</c>, default
    /// thermodynamic mode): every window of <see cref="Primer3ProbeSettings.MinSize"/>..<see cref="Primer3ProbeSettings.MaxSize"/>
    /// bases made only of A/C/G/T (PRIMER_INTERNAL_MAX_NS_ACCEPTED = 0) is accepted when its G+C % is within
    /// [MinGcPercent, MaxGcPercent], its longest mononucleotide run ≤ MaxPolyX, its Tm (Primer3 <c>seqtm</c> with
    /// MAX_NN_TM_LENGTH = 36,
    /// <see cref="ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3"/>, at the settings' conditions) within
    /// [MinTm, MaxTm], and its ntthal self-dimer, 3′ self-dimer and hairpin Tm
    /// (<see cref="ThermodynamicPrimerDesigner.CalculatePrimer3OligoStructure"/>) ≤ MaxSelfAnyTh / MaxSelfEndTh / MaxHairpinTh.
    /// As in Primer3, windows are enumerated per 3′ end with growing length, and an N, a poly-X run or a too-stable
    /// self-dimer (Primer3's "five-prime problems") stops the extension for that 3′ end, so a longer window that
    /// would pass is not considered. Accepted probes are ranked by the Primer3 penalty |Tm − OptTm| + |length − OptSize|
    /// (<see cref="ThermodynamicPrimerDesigner.CalculatePrimer3Penalty"/> with Primer3's default internal-oligo weights) and
    /// ordered as Primer3's <c>primer_rec_comp</c>: penalty ascending, then start descending, then length
    /// ascending. With <see cref="Primer3ProbeSettings.ThermodynamicOligoAlignment"/> = false (Primer3
    /// PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 0) the structure limits are Primer3's alignment-mode ones instead: dpal
    /// <c>self_any</c> &gt; <see cref="Primer3ProbeSettings.MaxSelfAny"/> (a five-prime problem) or <c>self_end</c> &gt;
    /// <see cref="Primer3ProbeSettings.MaxSelfEnd"/> (both 12.00 by default) reject the window, no hairpin check, and the
    /// values are reported in <see cref="Primer3Probe.SelfAny"/>/<see cref="Primer3Probe.SelfEnd"/>.
    /// Verified against primer3-py 2.3.1 <c>design_primers</c> (PRIMER_INTERNAL_n_* values, both modes).
    /// </summary>
    /// <param name="template">Template sequence (case-insensitive); probes are picked on this strand.</param>
    /// <param name="settings">Picker settings (default: Primer3 defaults).</param>
    /// <param name="numReturn">PRIMER_NUM_RETURN (default 5; at least 1, as Primer3's <c>_pr_data_control</c> requires).</param>
    /// <param name="sequenceQuality">SEQUENCE_QUALITY: one integer quality per template base (null or empty = none). With it
    /// every oligo gets <see cref="Primer3Probe.MinSequenceQuality"/>, is checked against
    /// <see cref="Primer3ProbeSettings.MinQuality"/> and weighted by <see cref="Primer3ProbeSettings.WeightSequenceQuality"/>;
    /// its length must equal the template length and its values lie in [<see cref="Primer3ProbeSettings.QualityRangeMin"/>,
    /// <see cref="Primer3ProbeSettings.QualityRangeMax"/>] (Primer3 <c>_pr_data_control</c>).</param>
    /// <returns>Up to <paramref name="numReturn"/> probes, best first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Invalid sizes (MinSize &lt; 1, MaxSize &lt; MinSize,
    /// MaxSize &gt; 36), an invalid mishyb limit, or <paramref name="numReturn"/> &lt; 1 (Primer3 "PRIMER_NUM_RETURN &lt; 1").</exception>
    /// <exception cref="ArgumentException">A GC weight without <see cref="Primer3ProbeSettings.OptGcPercent"/> or a mishyb weight
    /// without a mishyb library (Primer3 <c>_pr_data_control</c>).</exception>
    public static IReadOnlyList<Primer3Probe> DesignProbesPrimer3(
        string template,
        Primer3ProbeSettings? settings = null,
        int numReturn = 5,
        IReadOnlyList<int>? sequenceQuality = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        var s = settings ?? new Primer3ProbeSettings();
        ValidatePrimer3ProbeSettings(s, nameof(settings));
        if (numReturn < 1)
            throw new ArgumentOutOfRangeException(nameof(numReturn), "PRIMER_NUM_RETURN < 1 (Primer3 _pr_data_control).");
        ThermodynamicPrimerDesigner.ValidatePrimer3Quality(sequenceQuality, template.Length, 0, s.MinQuality, s.QualityRangeMin,
            s.QualityRangeMax, 0.0, s.WeightSequenceQuality, nameof(sequenceQuality));

        var accepted = EnumeratePrimer3InternalOligos(template.ToUpperInvariant(), 0, template.Length, s, screenStructure: true,
            sequenceQuality is { Count: > 0 } ? sequenceQuality : null, s.LowercaseMasking ? template : null);

        // primer_rec_comp: quality ascending, then start descending, then length ascending.
        accepted.Sort((a, b) => ComparePrimer3Rank(a.Penalty, a.Start, a.Length, b.Penalty, b.Start, b.Length));

        return accepted.Count > numReturn ? accepted.GetRange(0, numReturn) : accepted;
    }

    private static ThermodynamicPrimerDesigner.LibraryMispriming ComputePrimer3ProbeLibrary(
        string oligo, Primer3ProbeSettings s, PrimerMisprimingLibrary library) =>
        ThermodynamicPrimerDesigner.ComputeLibraryMispriming(oligo, true, library, s.LibraryAmbiguityCodesConsensus,
            s.MaxLibraryMishyb, isInternal: true);

    // Primer3 primer_rec_comp: penalty (quality) ascending, then start descending, then length ascending.
    private static int ComparePrimer3Rank(double penaltyA, int startA, int lengthA, double penaltyB, int startB, int lengthB)
    {
        int c = penaltyA.CompareTo(penaltyB);
        if (c != 0) return c;
        c = startB.CompareTo(startA);
        return c != 0 ? c : lengthA.CompareTo(lengthB);
    }
}
