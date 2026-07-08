namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region ctDNA analysis (ONCO-CTDNA-001)

    /// <summary>
    /// Computes the probability of detecting at least one circulating-tumor-DNA (ctDNA) molecule in a
    /// plasma cfDNA sample under the Poisson sampling model. With <paramref name="genomeEquivalents"/> n
    /// sequenced haploid genome equivalents, mutant allele fraction <paramref name="mutantAlleleFraction"/> d,
    /// and <paramref name="reporterCount"/> k independent tumour reporters, the expected number of mutant
    /// molecules is λ = n·d·k and the detection probability is
    /// <code>
    /// p = 1 − e^(−n·d·k)
    /// </code>
    /// Source: US Patent 11,085,084 B2 ("the probability of observing a single tumor reporter in cfDNA
    /// follows a Poisson distribution with mean λ = n × d … x = 1 − e^(−nd) … for k independent reporters
    /// p = 1 − e^(−ndk)"); Avanzini et al. (2020), <i>Science Advances</i> 6(50):eabc4308.
    /// </summary>
    /// <param name="genomeEquivalents">Number of sequenced haploid genome equivalents n (≥ 0).</param>
    /// <param name="mutantAlleleFraction">Mutant allele fraction d ∈ [0, 1] (the ctDNA detection-limit fraction).</param>
    /// <param name="reporterCount">Number of independent tumour reporters k (≥ 1). Default 1.</param>
    /// <returns>Detection probability p = 1 − e^(−n·d·k) ∈ [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="genomeEquivalents"/> &lt; 0, <paramref name="mutantAlleleFraction"/> ∉ [0, 1] or not finite,
    /// or <paramref name="reporterCount"/> &lt; 1.
    /// </exception>
    public static double CtDnaDetectionProbability(
        int genomeEquivalents, double mutantAlleleFraction, int reporterCount = 1)
    {
        if (genomeEquivalents < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(genomeEquivalents), genomeEquivalents, "Genome equivalents (n) cannot be negative.");
        }

        if (double.IsNaN(mutantAlleleFraction) || mutantAlleleFraction < 0.0 || mutantAlleleFraction > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mutantAlleleFraction), mutantAlleleFraction, "Mutant allele fraction (d) must be in [0, 1].");
        }

        if (reporterCount < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(reporterCount), reporterCount, "Reporter count (k) must be at least 1.");
        }

        // λ = n·d·k; p = 1 − e^(−λ). λ = 0 ⇒ p = 0; large λ ⇒ p → 1 (never exceeds 1).
        double lambda = (double)genomeEquivalents * mutantAlleleFraction * reporterCount;
        return 1.0 - Math.Exp(-lambda); // p = 1 − e^(−λ).
    }

    /// <summary>
    /// Expected number of mutant molecules λ = n·d·k for a ctDNA sample. Source: US Patent 11,085,084 B2 /
    /// Avanzini et al. (2020) — Poisson mean λ = n × d (× k reporters); corroborated by the worked example of
    /// Pessoa et al. (2023): n = 15,000 genome equivalents at d = 0.001 (0.1% VAF) ⇒ λ = 15 mutant molecules.
    /// </summary>
    /// <param name="genomeEquivalents">Sequenced haploid genome equivalents n (≥ 0).</param>
    /// <param name="mutantAlleleFraction">Mutant allele fraction d ∈ [0, 1].</param>
    /// <param name="reporterCount">Independent tumour reporters k (≥ 1). Default 1.</param>
    /// <returns>Expected mutant-molecule count λ = n·d·k.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Same domain limits as <see cref="CtDnaDetectionProbability"/>.</exception>
    public static double ExpectedMutantMolecules(
        int genomeEquivalents, double mutantAlleleFraction, int reporterCount = 1)
    {
        if (genomeEquivalents < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(genomeEquivalents), genomeEquivalents, "Genome equivalents (n) cannot be negative.");
        }

        if (double.IsNaN(mutantAlleleFraction) || mutantAlleleFraction < 0.0 || mutantAlleleFraction > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mutantAlleleFraction), mutantAlleleFraction, "Mutant allele fraction (d) must be in [0, 1].");
        }

        if (reporterCount < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(reporterCount), reporterCount, "Reporter count (k) must be at least 1.");
        }

        return (double)genomeEquivalents * mutantAlleleFraction * reporterCount;
    }

    /// <summary>
    /// Decides whether ctDNA is detectable: at least one mutant molecule must be expected (λ = n·d·k ≥ 1)
    /// AND the Poisson detection probability must reach <paramref name="minDetectionProbability"/>. Source:
    /// the Poisson detection model (US Patent 11,085,084 B2; Avanzini et al. 2020) gives the probability;
    /// Newman et al. (2014), <i>Nat. Med.</i> 20(5):548–554 establish a validated detection range
    /// (0.025%–10% allele fraction) at the conventional 95% operating point. The λ ≥ 1 requirement enforces
    /// the physical limit that no mutant molecule can be observed when fewer than one is expected.
    /// </summary>
    /// <param name="genomeEquivalents">Sequenced haploid genome equivalents n (≥ 0).</param>
    /// <param name="mutantAlleleFraction">Mutant allele fraction d ∈ [0, 1].</param>
    /// <param name="reporterCount">Independent tumour reporters k (≥ 1). Default 1.</param>
    /// <param name="minDetectionProbability">Minimum p to call detected; default <see cref="DefaultCtDnaDetectionProbability"/> (0.95).</param>
    /// <returns><c>true</c> if λ ≥ 1 and p ≥ <paramref name="minDetectionProbability"/>; otherwise <c>false</c>.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Domain limits of <see cref="CtDnaDetectionProbability"/>, or <paramref name="minDetectionProbability"/> ∉ (0, 1].
    /// </exception>
    public static bool IsCtDnaDetected(
        int genomeEquivalents,
        double mutantAlleleFraction,
        int reporterCount = 1,
        double minDetectionProbability = DefaultCtDnaDetectionProbability)
    {
        if (double.IsNaN(minDetectionProbability) || minDetectionProbability <= 0.0 || minDetectionProbability > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minDetectionProbability), minDetectionProbability,
                "Minimum detection probability must be in (0, 1].");
        }

        double lambda = ExpectedMutantMolecules(genomeEquivalents, mutantAlleleFraction, reporterCount);
        if (lambda < MinExpectedMutantMolecules)
        {
            return false;
        }

        double probability = 1.0 - Math.Exp(-lambda);
        return probability >= minDetectionProbability;
    }

    /// <summary>
    /// Estimates the ctDNA tumour fraction from clonal heterozygous somatic SNVs observed in plasma at
    /// copy-neutral diploid loci. For such a variant the expected VAF is half the tumour-derived fraction
    /// (v = TF/2), so tumour fraction = 2 · (mean VAF). The mean is taken over the supplied variants'
    /// plasma VAFs (alt / total reads). Source: Antonello et al. (2024), CNAqc, <i>Genome Biology</i> 25:38
    /// (m = 1, n_tot = 2 special case of v = m·π / [2(1−π) + π·n_tot] gives v = π/2). The result is clamped
    /// to [0, 1] since a fraction cannot exceed 1.
    /// </summary>
    /// <param name="variants">Clonal heterozygous somatic SNVs observed in plasma (each VAF ≤ 0.5).</param>
    /// <returns>Estimated ctDNA tumour fraction ∈ [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="variants"/> is empty (tumour fraction undefined).</exception>
    /// <exception cref="ArgumentOutOfRangeException">A variant has invalid read counts or a VAF &gt; 0.5.</exception>
    public static double CalculateTumorFraction(IEnumerable<VariantObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        double vafSum = 0.0;
        int count = 0;
        foreach (VariantObservation variant in variants)
        {
            double vaf = CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
            if (vaf > 0.5)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(variants), vaf,
                    "A clonal heterozygous SNV cannot have VAF > 0.5 under the diploid model; locus is not copy-neutral diploid heterozygous.");
            }

            vafSum += vaf;
            count++;
        }

        if (count == 0)
        {
            throw new ArgumentException("Cannot estimate tumour fraction from an empty variant set.", nameof(variants));
        }

        double meanVaf = vafSum / count;
        double tumorFraction = TumorFractionFromVafFactor * meanVaf;
        return Math.Min(tumorFraction, 1.0); // a fraction cannot exceed 1.
    }

    /// <summary>
    /// Computes the mean plasma variant allele fraction across ctDNA reporters: the arithmetic mean of each
    /// variant's VAF = alt reads / total reads. Source: Newman et al. (2014), <i>Nat. Med.</i> 20(5):548–554 —
    /// ctDNA level is summarized as the fraction of mutant molecules across SNV/indel reporters.
    /// </summary>
    /// <param name="variants">Observed ctDNA reporters with plasma read counts.</param>
    /// <returns>Mean variant allele fraction ∈ [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="variants"/> is empty (mean undefined).</exception>
    public static double CalculateMeanVaf(IEnumerable<VariantObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        double vafSum = 0.0;
        int count = 0;
        foreach (VariantObservation variant in variants)
        {
            vafSum += CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
            count++;
        }

        if (count == 0)
        {
            throw new ArgumentException("Cannot compute mean VAF from an empty variant set.", nameof(variants));
        }

        return vafSum / count;
    }

    /// <summary>
    /// Converts a cfDNA input mass in nanograms to the number of haploid genome equivalents, using the
    /// standard conversion 3.3 pg per haploid genome (Devonshire et al. 2014, PMC4182654), i.e.
    /// GE = ng · 1000 / 3.3 (≈ 303 GE per ng; Alcaide et al. 2020, <i>Sci. Rep.</i> 10:12564). This gives the
    /// sampling depth n used in the Poisson detection model <see cref="CtDnaDetectionProbability"/>.
    /// </summary>
    /// <param name="cfDnaNanograms">cfDNA input mass in nanograms (≥ 0).</param>
    /// <returns>Haploid genome equivalents (a continuous count; not rounded).</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="cfDnaNanograms"/> is negative or not finite.</exception>
    public static double HaploidGenomeEquivalents(double cfDnaNanograms)
    {
        if (double.IsNaN(cfDnaNanograms) || double.IsInfinity(cfDnaNanograms) || cfDnaNanograms < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(cfDnaNanograms), cfDnaNanograms, "cfDNA mass (ng) must be a finite non-negative value.");
        }

        return cfDnaNanograms * PicogramsPerNanogram / PicogramsPerHaploidGenome;
    }

    #endregion


    #region Minimal/Molecular Residual Disease (ONCO-MRD-001)

    /// <summary>
    /// Minimum number of tracked patient-specific variants that must be detected in plasma for a
    /// tumour-informed MRD assay to call the sample ctDNA-positive (MRD-positive). Source: Reinert et al.
    /// (2019), <i>JAMA Oncology</i> 5(8):1124–1131, as quoted in the tumour-informed ctDNA MRD review
    /// (PMC9265001, Table 1): "Plasma samples with at least two tumor-specific SNVs are defined as
    /// ctDNA-positive." (Signatera personalised 16-plex assay.)
    /// </summary>
    public const int DefaultMrdPositivityThreshold = 2;

    /// <summary>
    /// Minimum number of alternate (mutant) supporting reads at a tracked locus for that variant to be
    /// counted as detected in plasma. A variant with no supporting reads contributes no MRD signal.
    /// Source: Wan et al. (2020), <i>Sci. Transl. Med.</i> 12(548):eaaz8084 — per-locus signal must exceed
    /// background; the publishable universal read-count cutoff is error-model specific, so the minimal
    /// presence rule (≥ 1 alt read) is used by default and is configurable. The panel-level ≥2 calling rule
    /// (<see cref="DefaultMrdPositivityThreshold"/>) is independent of this value.
    /// </summary>
    public const int DefaultMrdMinSupportingReads = 1;

    /// <summary>
    /// A single patient-specific tumour marker tracked in plasma for MRD, with its plasma read evidence.
    /// </summary>
    /// <param name="Chromosome">Contig / chromosome identifier of the tracked locus.</param>
    /// <param name="Position">1-based reference position of the tracked variant.</param>
    /// <param name="ReferenceAllele">Reference allele.</param>
    /// <param name="AlternateAllele">Tumour-specific alternate allele tracked at this locus.</param>
    /// <param name="PlasmaAltReads">Alternate (mutant) supporting reads observed in plasma at this locus (≥ 0).</param>
    /// <param name="PlasmaTotalReads">Total covering reads in plasma at this locus (≥ 0).</param>
    public readonly record struct TumorMarker(
        string Chromosome,
        int Position,
        string ReferenceAllele,
        string AlternateAllele,
        int PlasmaAltReads,
        int PlasmaTotalReads);

    /// <summary>Binary MRD call for a plasma timepoint.</summary>
    public enum MrdStatus
    {
        /// <summary>Fewer than the positivity threshold of tracked variants detected: no residual disease signal.</summary>
        Negative,

        /// <summary>At least the positivity threshold of tracked variants detected: molecular residual disease present.</summary>
        Positive
    }

    /// <summary>
    /// Result of a tumour-informed MRD assessment of one plasma sample against a patient-specific marker panel.
    /// </summary>
    /// <param name="Status">MRD-positive when <see cref="DetectedVariantCount"/> ≥ positivity threshold; else negative.</param>
    /// <param name="DetectedVariantCount">Number of tracked variants detected (alt reads ≥ minSupportingReads).</param>
    /// <param name="TrackedVariantCount">Total number of tracked variants in the panel.</param>
    /// <param name="IntegratedMutantAlleleFraction">
    /// Depth-weighted (read-pooled) mean plasma VAF across tracked loci: Σ alt / Σ total (INVAR IMAF). 0 when no reads.
    /// </param>
    /// <param name="DetectionProbability">
    /// Panel-level Poisson probability of detecting ≥ 1 ctDNA molecule, p = 1 − e^(−n·f·m), at the observed IMAF.
    /// </param>
    public readonly record struct MrdResult(
        MrdStatus Status,
        int DetectedVariantCount,
        int TrackedVariantCount,
        double IntegratedMutantAlleleFraction,
        double DetectionProbability);

    /// <summary>Per-timepoint MRD call in a longitudinal series.</summary>
    /// <param name="TimepointIndex">0-based index of the timepoint in input order.</param>
    /// <param name="Result">The MRD assessment for that timepoint's plasma sample.</param>
    public readonly record struct MrdTimepoint(int TimepointIndex, MrdResult Result);

    /// <summary>Longitudinal MRD tracking across ordered plasma timepoints.</summary>
    /// <param name="Timepoints">Per-timepoint MRD results, preserving input order.</param>
    /// <param name="FirstPositiveIndex">0-based index of the earliest MRD-positive timepoint, or −1 if none.</param>
    public readonly record struct MrdLongitudinalResult(
        IReadOnlyList<MrdTimepoint> Timepoints,
        int FirstPositiveIndex);

    /// <summary>
    /// Reports whether a tracked tumour marker is detected in plasma: its alternate (mutant) supporting-read
    /// count is at or above <paramref name="minSupportingReads"/>. Source: Wan et al. (2020) — a locus
    /// contributes ctDNA signal only when supported by mutant reads above background.
    /// </summary>
    /// <param name="marker">The tracked marker with its plasma read evidence.</param>
    /// <param name="minSupportingReads">Minimum alt reads to count as detected (default <see cref="DefaultMrdMinSupportingReads"/>).</param>
    /// <returns><c>true</c> if <c>PlasmaAltReads ≥ minSupportingReads</c>; otherwise <c>false</c>.</returns>
    public static bool IsVariantDetected(TumorMarker marker, int minSupportingReads = DefaultMrdMinSupportingReads)
    {
        if (minSupportingReads < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minSupportingReads), minSupportingReads, "Minimum supporting reads must be at least 1.");
        }

        return marker.PlasmaAltReads >= minSupportingReads;
    }

    /// <summary>
    /// Tumour-informed minimal/molecular residual disease (MRD) detection. Given a panel of patient-specific
    /// somatic markers (selected from the tumour and tracked in plasma), counts how many are detected and
    /// calls the sample MRD-positive when the number of detected variants reaches
    /// <paramref name="positivityThreshold"/> (default 2 of up to 16). Source: Reinert et al. (2019),
    /// <i>JAMA Oncology</i> 5(8):1124–1131 / Signatera analytical-validation white paper — "Plasma samples
    /// with at least two tumor-specific SNVs are defined as ctDNA-positive." The integrated mutant allele
    /// fraction (IMAF) is the depth-weighted mean plasma VAF across loci (Wan et al. 2020), and the
    /// panel-level Poisson detection probability p = 1 − e^(−n·f·m) reuses
    /// <see cref="CtDnaDetectionProbability"/> with m = number of tracked markers (Signatera white paper,
    /// Figure 2; Avanzini et al. 2020).
    /// </summary>
    /// <param name="tumorMarkers">Patient-specific tracked markers with plasma read evidence (non-empty).</param>
    /// <param name="positivityThreshold">Minimum detected variants to call positive (default <see cref="DefaultMrdPositivityThreshold"/>).</param>
    /// <param name="minSupportingReads">Minimum alt reads for a marker to count as detected (default <see cref="DefaultMrdMinSupportingReads"/>).</param>
    /// <param name="genomeEquivalents">Sequenced haploid genome equivalents n for the panel Poisson p (≥ 0; default 0 ⇒ p = 0).</param>
    /// <returns>The MRD call with detected count, tracked count, IMAF, and panel detection probability.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="tumorMarkers"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="tumorMarkers"/> is empty (no panel to interrogate).</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold or <paramref name="genomeEquivalents"/> is out of range.</exception>
    public static MrdResult DetectMRD(
        IEnumerable<TumorMarker> tumorMarkers,
        int positivityThreshold = DefaultMrdPositivityThreshold,
        int minSupportingReads = DefaultMrdMinSupportingReads,
        int genomeEquivalents = 0)
    {
        ArgumentNullException.ThrowIfNull(tumorMarkers);

        if (positivityThreshold < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(positivityThreshold), positivityThreshold, "Positivity threshold must be at least 1.");
        }

        if (minSupportingReads < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(minSupportingReads), minSupportingReads, "Minimum supporting reads must be at least 1.");
        }

        if (genomeEquivalents < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(genomeEquivalents), genomeEquivalents, "Genome equivalents (n) cannot be negative.");
        }

        int trackedCount = 0;
        int detectedCount = 0;
        long altReadSum = 0;
        long totalReadSum = 0;
        foreach (TumorMarker marker in tumorMarkers)
        {
            trackedCount++;
            altReadSum += Math.Max(0, marker.PlasmaAltReads);
            totalReadSum += Math.Max(0, marker.PlasmaTotalReads);
            if (IsVariantDetected(marker, minSupportingReads))
            {
                detectedCount++;
            }
        }

        if (trackedCount == 0)
        {
            throw new ArgumentException(
                "Cannot assess MRD against an empty marker panel.", nameof(tumorMarkers));
        }

        // INVAR integrated mutant allele fraction: depth-weighted (read-pooled) plasma VAF across loci.
        double imaf = totalReadSum == 0 ? 0.0 : (double)altReadSum / totalReadSum;

        // Panel-level Poisson p = 1 - e^(-n*f*m) with f = IMAF, m = tracked markers (Signatera Fig 2).
        double detectionProbability = CtDnaDetectionProbability(genomeEquivalents, imaf, trackedCount);

        MrdStatus status = detectedCount >= positivityThreshold ? MrdStatus.Positive : MrdStatus.Negative;

        return new MrdResult(status, detectedCount, trackedCount, imaf, detectionProbability);
    }

    /// <summary>
    /// Longitudinal MRD tracking: applies <see cref="DetectMRD"/> to each timepoint's plasma marker panel,
    /// preserving input order, and reports the earliest MRD-positive timepoint (serial surveillance of
    /// residual disease, Reinert et al. 2019). Each timepoint is an independent plasma marker panel.
    /// </summary>
    /// <param name="timepoints">Ordered plasma marker panels, one per timepoint.</param>
    /// <param name="positivityThreshold">Minimum detected variants to call positive (default <see cref="DefaultMrdPositivityThreshold"/>).</param>
    /// <param name="minSupportingReads">Minimum alt reads for a marker to count as detected (default <see cref="DefaultMrdMinSupportingReads"/>).</param>
    /// <param name="genomeEquivalents">Sequenced haploid genome equivalents n for the panel Poisson p (≥ 0; default 0).</param>
    /// <returns>Per-timepoint MRD results and the 0-based index of the first positive timepoint (−1 if none).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="timepoints"/> is null, or a timepoint panel is null.</exception>
    /// <exception cref="ArgumentException">A timepoint panel is empty.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold or <paramref name="genomeEquivalents"/> is out of range.</exception>
    public static MrdLongitudinalResult TrackVariantsOverTime(
        IEnumerable<IEnumerable<TumorMarker>> timepoints,
        int positivityThreshold = DefaultMrdPositivityThreshold,
        int minSupportingReads = DefaultMrdMinSupportingReads,
        int genomeEquivalents = 0)
    {
        ArgumentNullException.ThrowIfNull(timepoints);

        var results = new List<MrdTimepoint>();
        int firstPositiveIndex = -1;
        int index = 0;
        foreach (IEnumerable<TumorMarker> panel in timepoints)
        {
            ArgumentNullException.ThrowIfNull(panel);
            MrdResult result = DetectMRD(panel, positivityThreshold, minSupportingReads, genomeEquivalents);
            if (firstPositiveIndex < 0 && result.Status == MrdStatus.Positive)
            {
                firstPositiveIndex = index;
            }

            results.Add(new MrdTimepoint(index, result));
            index++;
        }

        return new MrdLongitudinalResult(results, firstPositiveIndex);
    }

    // ---------------------------------------------------------------------------------------------
    // INVAR-style background-subtracted, tumour-AF-weighted ctDNA signal estimation (ONCO-MRD-001).
    //
    // Faithfully reproduces the core, caller-reproducible part of the INVAR pipeline (Wan et al. 2020,
    // Sci. Transl. Med. 12(548):eaaz8084; reference implementation INVAR2, nrlab-CRUK/INVAR2,
    // R/4_detection/generalisedLikelihoodRatioTest.R and R/shared/detectionFunctions.R):
    //  (a) per-locus / per-context BACKGROUND error subtraction (caller-supplied background AF per locus);
    //  (b) tumour-allele-fraction-weighted aggregation (IMAFv2) and a generalised-likelihood-ratio (GLRT)
    //      detection statistic whose mixture model weights each locus by its tumour AF vs background.
    // This no-size variant takes an already-cleaned background model per locus. The full INVAR pipeline's
    // fragment-length (size) weighting, patient-specific outlier suppression, locus-noise filtering and
    // control-derived background-error estimation are implemented separately below (see
    // EstimateInvarSignalWithSize, SuppressOutlierLoci, EstimateLocusBackground / PassesBothStrandsFilter).
    // ---------------------------------------------------------------------------------------------

    /// <summary>
    /// Initial value of the per-sample ctDNA fraction <c>p</c> for the INVAR EM maximum-likelihood
    /// estimator. Source: INVAR2 <c>estimate_p_EM</c> (R/shared/detectionFunctions.R) — <c>initial_p = 0.01</c>.
    /// </summary>
    private const double InvarEmInitialP = 0.01;

    /// <summary>
    /// Number of expectation-maximisation iterations used to estimate the ctDNA fraction <c>p</c>.
    /// Source: INVAR2 <c>estimate_p_EM</c> (R/shared/detectionFunctions.R) — <c>iterations = 200</c>.
    /// </summary>
    private const int InvarEmIterations = 200;

    /// <summary>
    /// A tracked patient-specific locus for INVAR-style ctDNA signal estimation: its plasma read evidence,
    /// the tumour allele fraction of the variant, and a caller-supplied background (non-reference) error rate.
    /// </summary>
    /// <param name="PlasmaAltReads">Mutant (alternate) supporting reads observed in plasma at this locus (≥ 0).</param>
    /// <param name="PlasmaTotalReads">Total covering reads in plasma at this locus (≥ 0).</param>
    /// <param name="TumorAlleleFraction">
    /// Tumour allele fraction of the tracked variant, <c>AF</c> in INVAR (0 &lt; AF ≤ 1). Loci with higher
    /// tumour AF carry more ctDNA signal and are weighted more strongly in the likelihood (Wan et al. 2020).
    /// </param>
    /// <param name="BackgroundErrorRate">
    /// Caller-supplied per-locus / per-trinucleotide-context background (non-reference) error rate <c>e</c>
    /// in INVAR (0 ≤ e &lt; 1), estimated from control plasma samples at the same loci. Subtracted from the
    /// observed plasma signal and used as the null read-error rate in the likelihood model.
    /// </param>
    public readonly record struct InvarLocus(
        int PlasmaAltReads,
        int PlasmaTotalReads,
        double TumorAlleleFraction,
        double BackgroundErrorRate);

    /// <summary>
    /// Result of an INVAR-style background-subtracted, tumour-AF-weighted ctDNA signal estimate for one sample.
    /// </summary>
    /// <param name="IntegratedMutantAlleleFractionV2">
    /// Background-subtracted, depth-weighted aggregate tumour fraction (INVAR2 IMAFv2): the depth-weighted mean
    /// over loci of <c>max(0, locusVAF − backgroundRate)</c>. ≥ 0; equals 0 when no locus exceeds background.
    /// </param>
    /// <param name="EstimatedTumorFraction">
    /// Maximum-likelihood ctDNA fraction <c>p̂</c> from the INVAR EM estimator under the AF-weighted mixture model.
    /// </param>
    /// <param name="LikelihoodRatio">
    /// Generalised-likelihood-ratio detection statistic: <c>logL(p̂) − logL(p = 0)</c> (per-locus mean scaled,
    /// as in INVAR2). Larger ⇒ stronger ctDNA evidence; ≈ 0 for a pure-background sample.
    /// </param>
    /// <param name="Detected">
    /// <c>true</c> when <see cref="LikelihoodRatio"/> reaches <paramref name="detectionThreshold"/> AND at least
    /// one mutant read is present; otherwise <c>false</c>.
    /// </param>
    /// <param name="LocusCount">Number of informative loci used (tumour AF &gt; 0).</param>
    public readonly record struct InvarSignalResult(
        double IntegratedMutantAlleleFractionV2,
        double EstimatedTumorFraction,
        double LikelihoodRatio,
        bool Detected,
        int LocusCount);

    /// <summary>
    /// Computes the INVAR2 background-subtracted, depth-weighted integrated mutant allele fraction (IMAFv2):
    /// the depth-weighted mean over loci of <c>max(0, plasmaVAF − backgroundRate)</c>. Source: INVAR2
    /// <c>calculateIMAFv2</c> (R/4_detection/generalisedLikelihoodRatioTest.R) — per-context
    /// <c>MEAN_AF.BS = pmax(0, MEAN_AF − BACKGROUND_AF)</c> then
    /// <c>weighted.mean(MEAN_AF.BS, TOTAL_DP)</c>.
    /// </summary>
    /// <param name="loci">The tracked informative loci (non-null). Loci with 0 total reads contribute 0 weight.</param>
    /// <returns>The background-subtracted depth-weighted tumour fraction (≥ 0); 0 when no covering reads.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    public static double IntegratedMutantAlleleFractionV2(IEnumerable<InvarLocus> loci)
    {
        ArgumentNullException.ThrowIfNull(loci);

        double weightedSum = 0.0;
        long totalDepth = 0;
        foreach (InvarLocus locus in loci)
        {
            int total = Math.Max(0, locus.PlasmaTotalReads);
            if (total == 0)
            {
                continue;
            }

            double vaf = (double)Math.Max(0, locus.PlasmaAltReads) / total;

            // Per-locus/per-context background subtraction: signal above background only (INVAR2 pmax(0, .)).
            double backgroundSubtracted = Math.Max(0.0, vaf - locus.BackgroundErrorRate);

            // Depth-weighted aggregation (INVAR2 weighted.mean(., TOTAL_DP)).
            weightedSum += backgroundSubtracted * total;
            totalDepth += total;
        }

        return totalDepth == 0 ? 0.0 : weightedSum / totalDepth;
    }

    /// <summary>
    /// INVAR-style estimate of residual ctDNA signal from a panel of tracked loci, each carrying plasma read
    /// evidence, a tumour allele fraction, and a caller-supplied per-locus background error rate. Performs
    /// (a) per-locus background subtraction, (b) tumour-AF-weighted aggregation (IMAFv2), (c) maximum-likelihood
    /// estimation of the ctDNA fraction <c>p̂</c> and (d) a generalised-likelihood-ratio detection statistic,
    /// faithfully following INVAR2 (Wan et al. 2020; nrlab-CRUK/INVAR2 detectionFunctions.R, no-size variant).
    ///
    /// <para>Mixture model (per read at a locus with tumour AF <c>AF</c>, background <c>e</c>, ctDNA fraction
    /// <c>p</c>): the probability a read is mutant is <c>q = AF·(1−e)·p + (1−AF)·e·p + e·(1−p)</c>; loci with
    /// higher <c>AF</c> and lower <c>e</c> contribute more signal, so the estimate is signal-to-noise weighted.</para>
    /// </summary>
    /// <param name="loci">Tracked informative loci with plasma evidence, tumour AF and background rate (non-empty).</param>
    /// <param name="detectionThreshold">
    /// Minimum generalised-likelihood-ratio statistic to call the sample ctDNA-positive (≥ 0; default 0 ⇒
    /// any positive evidence with a mutant read is detected). Larger values trade sensitivity for specificity.
    /// </param>
    /// <returns>The INVAR signal estimate: IMAFv2, ML ctDNA fraction, likelihood-ratio statistic and detection call.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="loci"/> has no informative locus (tumour AF &gt; 0).</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="detectionThreshold"/> is negative, or a locus has an out-of-range tumour AF
    /// (must be 0 &lt; AF ≤ 1 to be informative) or background rate (must be 0 ≤ e &lt; 1).
    /// </exception>
    public static InvarSignalResult EstimateInvarSignal(
        IEnumerable<InvarLocus> loci,
        double detectionThreshold = 0.0)
    {
        ArgumentNullException.ThrowIfNull(loci);

        if (detectionThreshold < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(detectionThreshold), detectionThreshold, "Detection threshold cannot be negative.");
        }

        var materialised = loci as IReadOnlyCollection<InvarLocus> ?? loci.ToList();

        // IMAFv2 is computed over all covered loci (background-subtracted, depth-weighted).
        double imafV2 = IntegratedMutantAlleleFractionV2(materialised);

        // Build per-informative-locus vectors for the likelihood model (INVAR uses one row per molecule,
        // i.e. depth = total covering reads, M = mutant reads). Only loci with tumour AF > 0 are informative.
        var altReads = new List<double>();
        var totalReads = new List<double>();
        var tumorAf = new List<double>();
        var background = new List<double>();
        bool anyMutantRead = false;
        foreach (InvarLocus locus in materialised)
        {
            if (locus.TumorAlleleFraction < 0.0 || locus.TumorAlleleFraction > 1.0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(loci), locus.TumorAlleleFraction, "Tumour allele fraction must be in [0, 1].");
            }

            if (locus.BackgroundErrorRate < 0.0 || locus.BackgroundErrorRate >= 1.0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(loci), locus.BackgroundErrorRate, "Background error rate must be in [0, 1).");
            }

            // INVAR keeps only loci with tumour AF > 0 (filter(TUMOUR_AF > 0)).
            if (locus.TumorAlleleFraction <= 0.0)
            {
                continue;
            }

            int total = Math.Max(0, locus.PlasmaTotalReads);
            if (total == 0)
            {
                continue;
            }

            int alt = Math.Clamp(locus.PlasmaAltReads, 0, total);

            // INVAR guards a zero background by flooring it to one expected error in the locus depth
            // (BACKGROUND_AF = ifelse(BACKGROUND_AF > 0, BACKGROUND_AF, 1 / BACKGROUND_DP)), so log(e) is finite.
            double e = locus.BackgroundErrorRate > 0.0 ? locus.BackgroundErrorRate : 1.0 / total;

            altReads.Add(alt);
            totalReads.Add(total);
            tumorAf.Add(locus.TumorAlleleFraction);
            background.Add(e);
            if (alt > 0)
            {
                anyMutantRead = true;
            }
        }

        if (altReads.Count == 0)
        {
            throw new ArgumentException(
                "No informative locus (tumour AF > 0 with covering reads) to estimate INVAR signal.", nameof(loci));
        }

        double pMle = EstimateCtDnaFractionEm(altReads, totalReads, tumorAf, background);
        double nullLogLik = InvarLogLikelihood(altReads, totalReads, tumorAf, background, 0.0);
        double altLogLik = InvarLogLikelihood(altReads, totalReads, tumorAf, background, pMle);
        double likelihoodRatio = altLogLik - nullLogLik;

        bool detected = anyMutantRead && likelihoodRatio >= detectionThreshold;

        return new InvarSignalResult(imafV2, pMle, likelihoodRatio, detected, altReads.Count);
    }

    /// <summary>
    /// EM maximum-likelihood estimate of the ctDNA fraction <c>p</c> under the INVAR mixture model.
    /// Source: INVAR2 <c>estimate_p_EM</c> (R/shared/detectionFunctions.R):
    /// <c>g = AF·(1−e) + (1−AF)·e</c>;
    /// E-step <c>Z0 = (1−g)·p / ((1−g)·p + (1−e)·(1−p))</c>, <c>Z1 = g·p / (g·p + e·(1−p))</c>;
    /// M-step <c>p = Σ(M·Z1 + (R−M)·Z0) / ΣR</c>.
    /// </summary>
    private static double EstimateCtDnaFractionEm(
        IReadOnlyList<double> m,
        IReadOnlyList<double> r,
        IReadOnlyList<double> af,
        IReadOnlyList<double> e)
    {
        double p = InvarEmInitialP;
        for (int iter = 0; iter < InvarEmIterations; iter++)
        {
            double numerator = 0.0;
            double denominator = 0.0;
            for (int i = 0; i < m.Count; i++)
            {
                double g = (af[i] * (1.0 - e[i])) + ((1.0 - af[i]) * e[i]);
                double z0 = ((1.0 - g) * p) / (((1.0 - g) * p) + ((1.0 - e[i]) * (1.0 - p)));
                double z1 = (g * p) / ((g * p) + (e[i] * (1.0 - p)));
                numerator += (m[i] * z1) + ((r[i] - m[i]) * z0);
                denominator += r[i];
            }

            p = denominator == 0.0 ? 0.0 : numerator / denominator;
        }

        return p;
    }

    /// <summary>
    /// Per-locus-mean log-likelihood of the sample given a ctDNA fraction <c>p</c> under the INVAR mixture model.
    /// Source: INVAR2 <c>calc_log_likelihood</c> (R/shared/detectionFunctions.R):
    /// <c>q = AF·(1−e)·p + (1−AF)·e·p + e·(1−p)</c>;
    /// <c>logL = Σ[ lchoose(R, M) + M·log(q) + (R−M)·log(1−q) ] / length(R)</c>.
    /// </summary>
    private static double InvarLogLikelihood(
        IReadOnlyList<double> m,
        IReadOnlyList<double> r,
        IReadOnlyList<double> af,
        IReadOnlyList<double> e,
        double p)
    {
        double sum = 0.0;
        for (int i = 0; i < m.Count; i++)
        {
            double q = (af[i] * (1.0 - e[i]) * p) + ((1.0 - af[i]) * e[i] * p) + (e[i] * (1.0 - p));

            // Clamp q strictly inside (0, 1) so the logs stay finite at the boundaries (p = 0, e = 0 ⇒ q = 0).
            q = Math.Clamp(q, double.Epsilon, 1.0 - double.Epsilon);

            double lchoose = LogChoose(r[i], m[i]);
            sum += lchoose + (m[i] * Math.Log(q)) + ((r[i] - m[i]) * Math.Log(1.0 - q));
        }

        return sum / m.Count;
    }

    /// <summary>
    /// Natural log of the binomial coefficient C(n, k) via log-gamma, matching R's <c>lchoose(n, k)</c>:
    /// <c>lgamma(n+1) − lgamma(k+1) − lgamma(n−k+1)</c>.
    /// </summary>
    private static double LogChoose(double n, double k)
    {
        if (k < 0.0 || k > n)
        {
            return double.NegativeInfinity;
        }

        return LogGamma(n + 1.0) - LogGamma(k + 1.0) - LogGamma(n - k + 1.0);
    }

    // ---------------------------------------------------------------------------------------------
    // INVAR fragment-size weighting, patient-specific outlier suppression, locus-noise filtering
    // and control-derived background-error estimation (ONCO-MRD-001 — closing the residual).
    //
    // Ports of the open-source INVAR2 pipeline (nrlab-CRUK/INVAR2; Wan et al. 2020, Sci. Transl. Med.
    // 12(548):eaaz8084):
    //  (1) size-weighted GLRT — calc_likelihood_ratio_with_RL / calc_log_likelihood_with_RL /
    //      estimate_p_EM_with_RL (R/shared/detectionFunctions.R);
    //  (2) outlier suppression — repolish (R/3_outlier_suppression/outlierSuppression.R);
    //  (3) locus-noise filtering and (4) background-error estimation from control samples —
    //      createLociErrorRateTable (R/1_parse/onTargetErrorRatesAndFilter.R).
    // ---------------------------------------------------------------------------------------------

    /// <summary>
    /// Default outlier-suppression threshold <c>α</c>: a tracked locus is flagged as a patient-specific
    /// outlier when its one-sided binomial tail probability under the sample ctDNA estimate falls at or below
    /// <c>α / (number of loci)</c> (Bonferroni). Source: INVAR2 <c>outlierSuppression.R</c> option
    /// <c>--outlier-suppression</c> default <c>0.05</c>.
    /// </summary>
    public const double InvarDefaultOutlierSuppression = 0.05;

    /// <summary>
    /// Default maximum per-locus allele fraction for a control/plasma locus to be used when estimating the
    /// sample ctDNA fraction during outlier suppression (loci above this are assumed to carry real signal and
    /// are excluded from the null estimate). Source: INVAR2 <c>--allele-frequency-threshold</c> default <c>0.01</c>.
    /// </summary>
    public const double InvarDefaultAlleleFrequencyThreshold = 0.01;

    /// <summary>
    /// Default maximum mutant-read count for a locus to contribute to the null ctDNA estimate during outlier
    /// suppression. Source: INVAR2 <c>--maximum-mutant-reads</c> default <c>10</c>.
    /// </summary>
    public const int InvarDefaultMaximumMutantReads = 10;

    /// <summary>
    /// Default maximum fraction of control samples that may carry any alt read at a locus before the locus is
    /// blacklisted as recurrently noisy. Source: INVAR2 <c>--control-proportion</c> default <c>0.1</c>
    /// (<c>createLociErrorRateTable</c>: <c>N_SAMPLES_WITH_SIGNAL / N_SAMPLES &lt; proportion_of_controls</c>).
    /// </summary>
    public const double InvarDefaultControlProportion = 0.1;

    /// <summary>
    /// Default maximum mean control background allele fraction for a locus to pass the locus-noise filter.
    /// Source: INVAR2 <c>--max-background-allele-frequency</c> default <c>0.01</c>
    /// (<c>createLociErrorRateTable</c>: <c>BACKGROUND_AF &lt; max_background_mean_AF</c>).
    /// </summary>
    public const double InvarDefaultMaxBackgroundAlleleFrequency = 0.01;

    /// <summary>
    /// One molecule (read/fragment) covering a tracked locus, used by the INVAR fragment-size-weighted GLRT.
    /// Each molecule is mutant or wild-type and carries its cfDNA fragment length, the tumour allele fraction of
    /// the tracked variant, and the per-locus background error rate. Source: INVAR2
    /// <c>calculateLikelihoodRatioForSampleWithSize</c> (one row per molecule, <c>DP = 1</c>).
    /// </summary>
    /// <param name="IsMutant"><c>true</c> when this molecule supports the tumour (alternate) allele.</param>
    /// <param name="FragmentLength">cfDNA fragment length in bp (&gt; 0).</param>
    /// <param name="TumorAlleleFraction">Tumour allele fraction <c>AF</c> of the tracked variant (0 &lt; AF ≤ 1).</param>
    /// <param name="BackgroundErrorRate">Per-locus background (non-reference) error rate <c>e</c> (0 ≤ e &lt; 1).</param>
    public readonly record struct InvarMolecule(
        bool IsMutant,
        int FragmentLength,
        double TumorAlleleFraction,
        double BackgroundErrorRate);

    /// <summary>
    /// An empirical cfDNA fragment-size profile: for each fragment length, the probability that a molecule of
    /// that length is drawn from the distribution (mutant/tumour or wild-type/normal). Tumour-derived cfDNA is
    /// shorter, so the mutant profile is enriched at short lengths. Source: INVAR2
    /// <c>estimate_real_length_probability</c> / <c>sizeCharacterisation.R</c> (per-size <c>PROPORTION = COUNT/TOTAL</c>).
    /// </summary>
    public sealed class FragmentSizeProfile
    {
        private readonly IReadOnlyDictionary<int, double> _mutantProbability;
        private readonly IReadOnlyDictionary<int, double> _normalProbability;
        private readonly double _uniformProbability;

        /// <summary>
        /// Builds a size profile from per-length molecule counts for mutant (tumour) and wild-type (normal)
        /// fragments. Each profile is normalised to a probability mass (INVAR2 <c>PROPORTION = COUNT/TOTAL</c>).
        /// Lengths absent from a profile fall back to the uniform probability <c>1/(maxLength−minLength+1)</c>
        /// over the supplied length range (INVAR2 <c>onlyWeighMutants</c> fallback).
        /// </summary>
        /// <param name="mutantCounts">Per-fragment-length molecule counts for tumour-derived (mutant) reads.</param>
        /// <param name="normalCounts">Per-fragment-length molecule counts for wild-type (normal) reads.</param>
        /// <param name="minLength">Minimum fragment length in the size window (INVAR2 default 60).</param>
        /// <param name="maxLength">Maximum fragment length in the size window (INVAR2 default 300).</param>
        public FragmentSizeProfile(
            IReadOnlyDictionary<int, int> mutantCounts,
            IReadOnlyDictionary<int, int> normalCounts,
            int minLength = InvarDefaultMinFragmentLength,
            int maxLength = InvarDefaultMaxFragmentLength)
        {
            ArgumentNullException.ThrowIfNull(mutantCounts);
            ArgumentNullException.ThrowIfNull(normalCounts);
            if (minLength <= 0 || maxLength < minLength)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(maxLength), maxLength, "Require 0 < minLength <= maxLength for the size window.");
            }

            _mutantProbability = Normalise(mutantCounts);
            _normalProbability = Normalise(normalCounts);

            // Uniform fall-back over the size window (INVAR2: 1/((max-min)+1)).
            _uniformProbability = 1.0 / ((maxLength - minLength) + 1);
        }

        /// <summary>
        /// Private constructor for precomputed per-length probability masses (used by the KDE factory).
        /// </summary>
        private FragmentSizeProfile(
            IReadOnlyDictionary<int, double> mutantProbability,
            IReadOnlyDictionary<int, double> normalProbability,
            double uniformProbability)
        {
            _mutantProbability = mutantProbability;
            _normalProbability = normalProbability;
            _uniformProbability = uniformProbability;
        }

        /// <summary>
        /// Builds an <b>opt-in KDE-smoothed</b> size profile: the per-length signal-vs-noise weights are obtained
        /// from a <b>Gaussian kernel density estimate</b> of the fragment-length distribution rather than the raw
        /// discrete <c>COUNT/TOTAL</c> histogram. This matches INVAR2's <c>estimate_real_length_probability</c>,
        /// which calls R's <c>density()</c> (Gaussian kernel) on the per-length counts and integrates the smoothed
        /// density over each integer-length bin <c>[L−0.5, L+0.5]</c>; tumour-derived cfDNA is shorter, so the
        /// smoothed mutant profile up-weights short fragments without the sparse-bin noise of the raw histogram.
        /// <para>
        /// Kernel estimator (Silverman 1986, eq. 2.2a): for observed lengths <c>xᵢ</c> with weights
        /// <c>wᵢ = countᵢ / Σcount</c> and bandwidth <c>h</c>, the smoothed density is
        /// <c>f̂(t) = Σᵢ wᵢ · φ((t − xᵢ)/h) / h</c> with <c>φ</c> the standard-normal pdf (the Gaussian kernel,
        /// which integrates to 1 — Silverman eq. 2.2). The per-length probability mass is the analytic integral
        /// over the bin, <c>P(L) = Σᵢ wᵢ · [Φ((L+0.5 − xᵢ)/h) − Φ((L−0.5 − xᵢ)/h)]</c>, renormalised over the
        /// integer support <c>[minLength, maxLength]</c> so the masses sum to 1.
        /// </para>
        /// </summary>
        /// <param name="mutantCounts">Per-fragment-length molecule counts for tumour-derived (mutant) reads.</param>
        /// <param name="normalCounts">Per-fragment-length molecule counts for wild-type (normal) reads.</param>
        /// <param name="bandwidth">
        /// Explicit Gaussian-kernel bandwidth <c>h</c> (&gt; 0) in bp, or <c>null</c> (default) to choose it by
        /// Silverman's rule of thumb (see <paramref name="bandwidthAdjust"/>).
        /// </param>
        /// <param name="bandwidthAdjust">
        /// Multiplier applied to the auto-selected bandwidth (R <c>density(adjust=…)</c>); the bandwidth used is
        /// <c>bandwidthAdjust · h_Silverman</c>. Ignored when <paramref name="bandwidth"/> is supplied. Must be &gt; 0;
        /// default <c>1.0</c> (INVAR2 uses <c>0.03</c>).
        /// </param>
        /// <param name="minLength">Minimum fragment length in the size window (INVAR2 default 60).</param>
        /// <param name="maxLength">Maximum fragment length in the size window (INVAR2 default 300).</param>
        /// <returns>A KDE-smoothed <see cref="FragmentSizeProfile"/>.</returns>
        /// <exception cref="ArgumentNullException">Either count map is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Invalid window, bandwidth, or adjust.</exception>
        public static FragmentSizeProfile FromKernelDensity(
            IReadOnlyDictionary<int, int> mutantCounts,
            IReadOnlyDictionary<int, int> normalCounts,
            double? bandwidth = null,
            double bandwidthAdjust = 1.0,
            int minLength = InvarDefaultMinFragmentLength,
            int maxLength = InvarDefaultMaxFragmentLength)
        {
            ArgumentNullException.ThrowIfNull(mutantCounts);
            ArgumentNullException.ThrowIfNull(normalCounts);
            if (minLength <= 0 || maxLength < minLength)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(maxLength), maxLength, "Require 0 < minLength <= maxLength for the size window.");
            }

            if (bandwidth is <= 0.0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(bandwidth), bandwidth, "Bandwidth must be positive.");
            }

            if (bandwidthAdjust <= 0.0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(bandwidthAdjust), bandwidthAdjust, "Bandwidth adjust must be positive.");
            }

            IReadOnlyDictionary<int, double> mutant =
                SmoothToProbabilities(mutantCounts, bandwidth, bandwidthAdjust, minLength, maxLength);
            IReadOnlyDictionary<int, double> normal =
                SmoothToProbabilities(normalCounts, bandwidth, bandwidthAdjust, minLength, maxLength);

            double uniform = 1.0 / ((maxLength - minLength) + 1);
            return new FragmentSizeProfile(mutant, normal, uniform);
        }

        /// <summary>
        /// Smooths per-length counts into a per-integer-length probability mass via a weighted Gaussian KDE,
        /// integrated analytically over each integer bin and renormalised over the support. Empty / single-point
        /// inputs return an empty map (so the caller's uniform fall-back applies, mirroring INVAR2's guard
        /// <c>if (length(counts) &gt; 1)</c> in <c>estimate_real_length_probability</c>).
        /// </summary>
        private static IReadOnlyDictionary<int, double> SmoothToProbabilities(
            IReadOnlyDictionary<int, int> counts,
            double? bandwidth,
            double bandwidthAdjust,
            int minLength,
            int maxLength)
        {
            // Collect observed lengths with positive counts.
            var lengths = new List<int>(counts.Count);
            var weights = new List<double>(counts.Count);
            long total = 0;
            foreach (KeyValuePair<int, int> kv in counts)
            {
                int c = Math.Max(0, kv.Value);
                if (c > 0)
                {
                    lengths.Add(kv.Key);
                    weights.Add(c);
                    total += c;
                }
            }

            var result = new Dictionary<int, double>();

            // INVAR2 guard: need more than one distinct point to estimate a density; otherwise fall back to uniform.
            if (lengths.Count <= 1 || total == 0)
            {
                return result;
            }

            // Weights normalise to a probability mass (INVAR2: weights <- counts / sum(counts)).
            for (int i = 0; i < weights.Count; i++)
            {
                weights[i] /= total;
            }

            double h = bandwidth ?? (bandwidthAdjust * SilvermanBandwidth(lengths));

            // Analytic integral of the weighted Gaussian KDE over each integer bin [L-0.5, L+0.5]:
            //   P(L) = Σ wᵢ · [Φ((L+0.5 − xᵢ)/h) − Φ((L−0.5 − xᵢ)/h)]   (Φ = standard-normal CDF).
            double sum = 0.0;
            for (int len = minLength; len <= maxLength; len++)
            {
                double mass = 0.0;
                for (int i = 0; i < lengths.Count; i++)
                {
                    double upper = StandardNormalCdf(((len + 0.5) - lengths[i]) / h);
                    double lower = StandardNormalCdf(((len - 0.5) - lengths[i]) / h);
                    mass += weights[i] * (upper - lower);
                }

                // Store every in-support length (even a ~0 tail mass) so the smoothed profile is self-contained and
                // the uniform fall-back never applies inside the support — that fall-back is the DISCRETE profile's
                // behaviour for unobserved lengths, not the KDE's (whose smoothed mass is defined everywhere).
                result[len] = mass;
                sum += mass;
            }

            // Renormalise over the integer support so the per-length masses sum to 1.
            if (sum > 0.0)
            {
                var keys = new List<int>(result.Keys);
                foreach (int len in keys)
                {
                    result[len] /= sum;
                }
            }

            return result;
        }

        /// <summary>
        /// Silverman's rule-of-thumb bandwidth for a Gaussian kernel (Silverman 1986, eq. 3.31; R
        /// <c>bw.nrd0</c>): <c>h = 0.9 · min(σ̂, IQR/1.34) · n^(−1/5)</c> over the per-read sample obtained by
        /// expanding each observed length by its weight. The robust scale <c>min(σ̂, IQR/1.34)</c> falls back to
        /// <c>σ̂</c>, then <c>|x₁|</c>, then <c>1</c> when the quartiles coincide (R <c>bw.nrd0</c> guard).
        /// </summary>
        private static double SilvermanBandwidth(IReadOnlyList<int> distinctLengths)
        {
            // Use the distinct observed lengths as the sample (matching R's density(x): bw.nrd0 is computed on the
            // supplied length vector and ignores the kernel weights). At least two points are guaranteed here.
            int n = distinctLengths.Count;
            double mean = 0.0;
            foreach (int v in distinctLengths)
            {
                mean += v;
            }

            mean /= n;

            double ss = 0.0;
            foreach (int v in distinctLengths)
            {
                double d = v - mean;
                ss += d * d;
            }

            // Sample standard deviation (R sd() uses the n−1 divisor).
            double sd = Math.Sqrt(ss / (n - 1));

            var sorted = new List<int>(distinctLengths);
            sorted.Sort();
            double iqr = Quantile(sorted, 0.75) - Quantile(sorted, 0.25);

            // Robust scale with R bw.nrd0's fall-back chain (1.34 normalises the IQR of a Gaussian to its σ).
            double lo = Math.Min(sd, iqr / 1.34);
            if (lo <= 0.0)
            {
                double loFallback = Math.Abs(sorted[0]) > 0.0 ? Math.Abs(sorted[0]) : 1.0;
                lo = sd > 0.0 ? sd : loFallback;
            }

            return 0.9 * lo * Math.Pow(n, -0.2);
        }

        /// <summary>Type-7 (R default) linear-interpolation quantile of a sorted sample.</summary>
        private static double Quantile(IReadOnlyList<int> sorted, double probability)
        {
            int n = sorted.Count;
            if (n == 1)
            {
                return sorted[0];
            }

            // R type 7: h = (n − 1)·p; interpolate between order statistics h_floor and h_floor+1.
            double h = (n - 1) * probability;
            int lo = (int)Math.Floor(h);
            int hi = Math.Min(lo + 1, n - 1);
            return sorted[lo] + ((h - lo) * (sorted[hi] - sorted[lo]));
        }

        /// <summary>Probability that a tumour-derived (mutant) molecule has the given fragment length.</summary>
        public double MutantProbability(int fragmentLength) =>
            _mutantProbability.TryGetValue(fragmentLength, out double p) ? p : _uniformProbability;

        /// <summary>Probability that a wild-type (normal) molecule has the given fragment length.</summary>
        public double NormalProbability(int fragmentLength) =>
            _normalProbability.TryGetValue(fragmentLength, out double p) ? p : _uniformProbability;

        private static IReadOnlyDictionary<int, double> Normalise(IReadOnlyDictionary<int, int> counts)
        {
            long total = 0;
            foreach (KeyValuePair<int, int> kv in counts)
            {
                total += Math.Max(0, kv.Value);
            }

            var result = new Dictionary<int, double>(counts.Count);
            if (total == 0)
            {
                return result;
            }

            foreach (KeyValuePair<int, int> kv in counts)
            {
                result[kv.Key] = Math.Max(0, kv.Value) / (double)total;
            }

            return result;
        }

        /// <summary>
        /// Standard-normal cumulative distribution function <c>Φ(z) = ½·[1 + erf(z/√2)]</c>. Used to integrate the
        /// Gaussian kernel analytically over each fragment-length bin (the Gaussian kernel's antiderivative is
        /// <c>Φ</c>). Reference: <c>Φ(z) = ½[1 + erf(z/√2)]</c>.
        /// </summary>
        private static double StandardNormalCdf(double z) => 0.5 * (1.0 + Erf(z / Sqrt2));

        /// <summary>1/√2, used in the Φ(z) = ½[1 + erf(z/√2)] identity.</summary>
        private static readonly double Sqrt2 = Math.Sqrt(2.0);

        /// <summary>2/√π, the leading constant of the error-function series erf(x) = (2/√π)·Σ …</summary>
        private static readonly double TwoOverSqrtPi = 2.0 / Math.Sqrt(Math.PI);

        /// <summary>1/√π, the leading constant of the complementary-error continued fraction.</summary>
        private static readonly double InvSqrtPi = 1.0 / Math.Sqrt(Math.PI);

        /// <summary>
        /// |x| beyond which the Maclaurin series loses precision; above it the complementary continued fraction is
        /// used instead. Both are exact (double-precision-limited); 2.0 keeps each well-conditioned.
        /// </summary>
        private const double ErfSeriesBound = 2.0;

        /// <summary>
        /// Gauss error function, exact to double precision over the whole real line via two first-principles
        /// expansions (no tabulated approximation coefficients):
        /// <list type="bullet">
        /// <item>|x| ≤ 2 — the everywhere-convergent Maclaurin series
        /// <c>erf(x) = (2/√π)·Σ_{k≥0} (−1)ᵏ·x^(2k+1) / (k!·(2k+1))</c>, with each term derived from the previous by
        /// the ratio <c>t_k = t_{k−1}·(−x²)·(2k−1)/(k·(2k+1))</c>.</item>
        /// <item>|x| &gt; 2 — <c>erf(x) = sign(x)·(1 − erfc(|x|))</c> with <c>erfc</c> from its classical
        /// continued fraction <c>erfc(x) = (e^{−x²}/√π)·1/(x + ½/(x + 1/(x + 3⁄2/(x + 2/(x + …)))))</c>, which
        /// converges rapidly for <c>x ≥ 2</c>.</item>
        /// </list>
        /// </summary>
        private static double Erf(double x)
        {
            if (x == 0.0)
            {
                return 0.0;
            }

            double ax = Math.Abs(x);
            if (ax <= ErfSeriesBound)
            {
                double x2 = x * x;
                double term = x;          // k = 0 term: x / (0! · 1) = x.
                double sum = x;
                for (int k = 1; k < 200; k++)
                {
                    term *= -x2 * ((2 * k) - 1) / (k * ((2 * k) + 1));
                    double previous = sum;
                    sum += term;
                    if (sum == previous)
                    {
                        break; // Converged to machine precision.
                    }
                }

                return TwoOverSqrtPi * sum;
            }

            // Lentz evaluation of the continued fraction erfc(x) = (e^{-x²}/√π) / (x + a1/(x + a2/(x + …)))
            // with partial numerators a_k = k/2. Exact (double-precision-limited) for x ≥ 2.
            double cf = ax;
            for (int k = 60; k >= 1; k--)
            {
                cf = ax + ((k / 2.0) / cf);
            }

            double erfc = Math.Exp(-ax * ax) * InvSqrtPi / cf;
            double erf = 1.0 - erfc;
            return x < 0.0 ? -erf : erf;
        }
    }

    /// <summary>Default minimum cfDNA fragment length for the INVAR size window. Source: INVAR2
    /// <c>--minimum-fragment-length</c> default <c>60</c>.</summary>
    public const int InvarDefaultMinFragmentLength = 60;

    /// <summary>Default maximum cfDNA fragment length for the INVAR size window. Source: INVAR2
    /// <c>--maximum-fragment-length</c> default <c>300</c>.</summary>
    public const int InvarDefaultMaxFragmentLength = 300;

    /// <summary>
    /// INVAR fragment-size-weighted generalised-likelihood-ratio statistic. Each molecule contributes a size
    /// likelihood: a read of length L is weighted by the ratio of its probability under the tumour size profile
    /// (P1) to the normal size profile (P0), so that short, tumour-like fragments carry more ctDNA evidence.
    /// Source: INVAR2 <c>calc_likelihood_ratio_with_RL</c> / <c>calc_log_likelihood_with_RL</c> /
    /// <c>estimate_p_EM_with_RL</c> (R/shared/detectionFunctions.R):
    /// per molecule <c>L0 = (1−e)·P0·(1−p) + (1−g)·P1·p</c>, <c>L1 = e·P0·(1−p) + g·P1·p</c> with
    /// <c>g = AF·(1−e) + (1−AF)·e</c>; <c>logL = Σ[ M·log(L1) + (1−M)·log(L0) ] / n</c>; the EM and the
    /// detection statistic <c>LR = logL(p̂) − logL(0)</c> mirror the no-size variant.
    /// </summary>
    /// <param name="molecules">Per-molecule observations covering the tracked loci (non-empty, AF &gt; 0).</param>
    /// <param name="sizeProfile">Mutant/normal cfDNA fragment-size profiles (non-null).</param>
    /// <param name="detectionThreshold">Minimum LR to call ctDNA-positive (≥ 0; default 0).</param>
    /// <returns>The size-weighted INVAR signal estimate (ML <c>p̂</c>, LR, detection call).</returns>
    /// <exception cref="ArgumentNullException">Any argument is null.</exception>
    /// <exception cref="ArgumentException">No informative molecule (tumour AF &gt; 0).</exception>
    /// <exception cref="ArgumentOutOfRangeException">Negative threshold or out-of-range AF / background.</exception>
    public static InvarSignalResult EstimateInvarSignalWithSize(
        IEnumerable<InvarMolecule> molecules,
        FragmentSizeProfile sizeProfile,
        double detectionThreshold = 0.0)
    {
        ArgumentNullException.ThrowIfNull(molecules);
        ArgumentNullException.ThrowIfNull(sizeProfile);
        if (detectionThreshold < 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(detectionThreshold), detectionThreshold, "Detection threshold cannot be negative.");
        }

        var af = new List<double>();
        var background = new List<double>();
        var mutant = new List<double>();
        var p0 = new List<double>();
        var p1 = new List<double>();
        bool anyMutant = false;
        foreach (InvarMolecule mol in molecules)
        {
            if (mol.TumorAlleleFraction < 0.0 || mol.TumorAlleleFraction > 1.0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(molecules), mol.TumorAlleleFraction, "Tumour allele fraction must be in [0, 1].");
            }

            if (mol.BackgroundErrorRate < 0.0 || mol.BackgroundErrorRate >= 1.0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(molecules), mol.BackgroundErrorRate, "Background error rate must be in [0, 1).");
            }

            // Only informative molecules (tumour AF > 0) contribute (INVAR filter(TUMOUR_AF > 0)).
            if (mol.TumorAlleleFraction <= 0.0)
            {
                continue;
            }

            af.Add(mol.TumorAlleleFraction);
            // Background floored to a tiny positive value so log(L) stays finite (INVAR 1/BACKGROUND_DP guard).
            background.Add(mol.BackgroundErrorRate > 0.0 ? mol.BackgroundErrorRate : 1.0 / Math.Max(1, af.Count));
            mutant.Add(mol.IsMutant ? 1.0 : 0.0);
            p0.Add(sizeProfile.NormalProbability(mol.FragmentLength));
            p1.Add(sizeProfile.MutantProbability(mol.FragmentLength));
            if (mol.IsMutant)
            {
                anyMutant = true;
            }
        }

        if (af.Count == 0)
        {
            throw new ArgumentException(
                "No informative molecule (tumour AF > 0) to estimate INVAR signal.", nameof(molecules));
        }

        double pMle = EstimateCtDnaFractionEmWithRl(mutant, af, background, p0, p1);
        double nullLogLik = InvarLogLikelihoodWithRl(mutant, af, background, p0, p1, 0.0);
        double altLogLik = InvarLogLikelihoodWithRl(mutant, af, background, p0, p1, pMle);
        double likelihoodRatio = altLogLik - nullLogLik;

        bool detected = anyMutant && likelihoodRatio >= detectionThreshold;
        return new InvarSignalResult(double.NaN, pMle, likelihoodRatio, detected, af.Count);
    }

    /// <summary>
    /// EM estimate of the ctDNA fraction <c>p</c> under the INVAR fragment-size-weighted mixture (per molecule).
    /// Source: INVAR2 <c>estimate_p_EM_with_RL</c> (R/shared/detectionFunctions.R):
    /// E-step <c>Z0 = (1−g)·P1·p / ((1−g)·P1·p + (1−e)·P0·(1−p))</c>,
    /// <c>Z1 = g·P1·p / (g·P1·p + e·P0·(1−p))</c>; M-step <c>p = Σ(M·Z1 + (1−M)·Z0) / n</c>.
    /// </summary>
    private static double EstimateCtDnaFractionEmWithRl(
        IReadOnlyList<double> m,
        IReadOnlyList<double> af,
        IReadOnlyList<double> e,
        IReadOnlyList<double> p0,
        IReadOnlyList<double> p1)
    {
        double p = InvarEmInitialP;
        for (int iter = 0; iter < InvarEmIterations; iter++)
        {
            double numerator = 0.0;
            for (int i = 0; i < m.Count; i++)
            {
                double g = (af[i] * (1.0 - e[i])) + ((1.0 - af[i]) * e[i]);
                double intNorm = (1.0 - g) * p1[i] * p;
                double z0 = intNorm / (intNorm + ((1.0 - e[i]) * p0[i] * (1.0 - p)));
                double intMut = g * p1[i] * p;
                double z1 = intMut / (intMut + (e[i] * p0[i] * (1.0 - p)));
                numerator += (m[i] * z1) + ((1.0 - m[i]) * z0);
            }

            p = numerator / m.Count;
        }

        return p;
    }

    /// <summary>
    /// Per-molecule-mean log-likelihood of the sample under the INVAR fragment-size-weighted mixture.
    /// Source: INVAR2 <c>calc_log_likelihood_with_RL</c> (R/shared/detectionFunctions.R):
    /// <c>L0 = (1−e)·P0·(1−p) + (1−g)·P1·p</c>, <c>L1 = e·P0·(1−p) + g·P1·p</c>;
    /// <c>logL = Σ[ M·log(L1) + (1−M)·log(L0) ] / n</c>.
    /// </summary>
    private static double InvarLogLikelihoodWithRl(
        IReadOnlyList<double> m,
        IReadOnlyList<double> af,
        IReadOnlyList<double> e,
        IReadOnlyList<double> p0,
        IReadOnlyList<double> p1,
        double p)
    {
        double sum = 0.0;
        for (int i = 0; i < m.Count; i++)
        {
            double g = (af[i] * (1.0 - e[i])) + ((1.0 - af[i]) * e[i]);
            double l0 = ((1.0 - e[i]) * p0[i] * (1.0 - p)) + ((1.0 - g) * p1[i] * p);
            double l1 = (e[i] * p0[i] * (1.0 - p)) + (g * p1[i] * p);

            // Keep the logs finite at the boundaries (p = 0, profile zero).
            l0 = Math.Max(l0, double.Epsilon);
            l1 = Math.Max(l1, double.Epsilon);
            sum += (m[i] * Math.Log(l1)) + ((1.0 - m[i]) * Math.Log(l0));
        }

        return sum / m.Count;
    }

    /// <summary>
    /// Result of INVAR patient-specific outlier suppression for one tracked locus.
    /// </summary>
    /// <param name="Locus">The tracked locus.</param>
    /// <param name="IsOutlier">
    /// <c>true</c> when the locus carries more mutant signal than the sample-wide ctDNA estimate explains —
    /// its one-sided binomial tail probability is ≤ the Bonferroni outlier threshold (INVAR <c>OUTLIER.PASS = FALSE</c>).
    /// </param>
    /// <param name="BinomialTailProbability">
    /// One-sided binomial probability <c>P(X ≥ observed mutant reads)</c> under the sample ctDNA estimate
    /// (INVAR2 <c>binom.test(..., alternative = "greater")</c>).
    /// </param>
    public readonly record struct InvarOutlierResult(
        InvarLocus Locus,
        bool IsOutlier,
        double BinomialTailProbability);

    /// <summary>
    /// INVAR patient-specific outlier suppression. Estimates the sample ctDNA fraction from the loci that are
    /// consistent with the null (low AF, few mutant reads), then flags any locus whose mutant-read count is a
    /// one-sided binomial outlier under that estimate, using a Bonferroni-corrected threshold over the loci.
    /// Source: INVAR2 <c>repolish</c> (R/3_outlier_suppression/outlierSuppression.R):
    /// <c>P_THRESHOLD = outlierSuppression / n_distinct(loci)</c>;
    /// <c>P_ESTIMATE = max(estimate_p_EM(...), weighted.mean(AF, TUMOUR_AF))</c>;
    /// <c>OUTLIER.PASS = binom.test(mutReads, DP, P_ESTIMATE, "greater")$p.value &gt; P_THRESHOLD</c>.
    /// </summary>
    /// <param name="loci">Tracked informative loci (non-empty).</param>
    /// <param name="outlierSuppression">Family-wise outlier threshold α (default 0.05).</param>
    /// <param name="alleleFrequencyThreshold">Max per-locus AF used in the null estimate (default 0.01).</param>
    /// <param name="maximumMutantReads">Max mutant reads for a locus in the null estimate (default 10).</param>
    /// <returns>Per-locus outlier verdicts in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="loci"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="loci"/> is empty.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold is out of range.</exception>
    public static IReadOnlyList<InvarOutlierResult> SuppressOutlierLoci(
        IEnumerable<InvarLocus> loci,
        double outlierSuppression = InvarDefaultOutlierSuppression,
        double alleleFrequencyThreshold = InvarDefaultAlleleFrequencyThreshold,
        int maximumMutantReads = InvarDefaultMaximumMutantReads)
    {
        ArgumentNullException.ThrowIfNull(loci);
        if (outlierSuppression <= 0.0 || outlierSuppression > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(outlierSuppression), outlierSuppression, "Outlier suppression α must be in (0, 1].");
        }

        if (alleleFrequencyThreshold <= 0.0 || alleleFrequencyThreshold > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(alleleFrequencyThreshold), alleleFrequencyThreshold, "AF threshold must be in (0, 1].");
        }

        if (maximumMutantReads < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(maximumMutantReads), maximumMutantReads, "Maximum mutant reads cannot be negative.");
        }

        var materialised = loci as IReadOnlyList<InvarLocus> ?? loci.ToList();
        if (materialised.Count == 0)
        {
            throw new ArgumentException("Cannot suppress outliers in an empty locus set.", nameof(loci));
        }

        // Null ctDNA estimate from loci consistent with background (INVAR repolish forPEstimate filter).
        var nullM = new List<double>();
        var nullR = new List<double>();
        var nullAf = new List<double>();
        var nullBg = new List<double>();
        double afWeightSum = 0.0;
        double afWeight = 0.0;
        foreach (InvarLocus locus in materialised)
        {
            int total = Math.Max(0, locus.PlasmaTotalReads);
            if (total == 0 || locus.TumorAlleleFraction <= 0.0)
            {
                continue;
            }

            int alt = Math.Clamp(locus.PlasmaAltReads, 0, total);
            double vaf = (double)alt / total;
            if (vaf > alleleFrequencyThreshold || alt > maximumMutantReads)
            {
                continue;
            }

            double e = locus.BackgroundErrorRate > 0.0 ? locus.BackgroundErrorRate : 1.0 / total;
            nullM.Add(alt);
            nullR.Add(total);
            nullAf.Add(locus.TumorAlleleFraction);
            nullBg.Add(e);

            // weighted.mean(AF, TUMOUR_AF): weight the observed VAF by tumour AF.
            afWeightSum += vaf * locus.TumorAlleleFraction;
            afWeight += locus.TumorAlleleFraction;
        }

        double pEstimate = 0.0;
        if (nullM.Count > 0)
        {
            double emP = EstimateCtDnaFractionEm(nullM, nullR, nullAf, nullBg);
            double weightedMean = afWeight > 0.0 ? afWeightSum / afWeight : 0.0;
            pEstimate = Math.Max(emP, weightedMean);
        }

        // Bonferroni threshold over the distinct loci (INVAR P_THRESHOLD).
        double pThreshold = outlierSuppression / materialised.Count;

        var results = new List<InvarOutlierResult>(materialised.Count);
        foreach (InvarLocus locus in materialised)
        {
            int total = Math.Max(0, locus.PlasmaTotalReads);
            int alt = total == 0 ? 0 : Math.Clamp(locus.PlasmaAltReads, 0, total);
            double tail = BinomialUpperTail(alt, total, pEstimate);
            bool isOutlier = tail <= pThreshold;
            results.Add(new InvarOutlierResult(locus, isOutlier, tail));
        }

        return results;
    }

    /// <summary>
    /// One-sided upper-tail binomial probability <c>P(X ≥ x)</c> for <c>x</c> successes in <c>n</c> trials with
    /// success probability <c>p</c>, matching R's <c>binom.test(x, n, p, alternative = "greater")$p.value</c>.
    /// Returns 1 for <c>x ≤ 0</c> (INVAR <c>ifelse(x &lt;= 0, 1, ...)</c>).
    /// </summary>
    private static double BinomialUpperTail(int x, int n, double p)
    {
        if (x <= 0 || n <= 0)
        {
            return 1.0;
        }

        if (x > n)
        {
            return 0.0;
        }

        // Numerically stable: sum exp(lchoose + i·log p + (n−i)·log(1−p)) over i = x..n.
        if (p <= 0.0)
        {
            return 0.0; // No mutant reads expected under p = 0, yet x ≥ 1 ⇒ tail mass 0.
        }

        if (p >= 1.0)
        {
            return 1.0;
        }

        double logP = Math.Log(p);
        double log1MinusP = Math.Log(1.0 - p);
        double tail = 0.0;
        for (int i = x; i <= n; i++)
        {
            tail += Math.Exp(LogChoose(n, i) + (i * logP) + ((n - i) * log1MinusP));
        }

        return Math.Min(1.0, tail);
    }

    /// <summary>
    /// Per-control-sample observation of a tracked locus, used to estimate the per-locus background error rate
    /// and the locus-noise (blacklist) flag from a panel of control (non-cancer) plasma samples.
    /// Source: INVAR2 <c>createLociErrorRateTable</c> (R/1_parse/onTargetErrorRatesAndFilter.R).
    /// </summary>
    /// <param name="ControlSampleId">Identifier of the control sample this observation comes from.</param>
    /// <param name="AltForwardReads">Alt reads on the forward strand (<c>ALT_F</c>, ≥ 0).</param>
    /// <param name="AltReverseReads">Alt reads on the reverse strand (<c>ALT_R</c>, ≥ 0).</param>
    /// <param name="TotalReads">Total covering reads at the locus in this sample (<c>DP</c>, ≥ 0).</param>
    public readonly record struct ControlLocusObservation(
        string ControlSampleId,
        int AltForwardReads,
        int AltReverseReads,
        int TotalReads);

    /// <summary>
    /// Estimated background error model and locus-noise verdict for one tracked locus, derived from control samples.
    /// </summary>
    /// <param name="BackgroundErrorRate">
    /// Pooled control background allele fraction <c>BACKGROUND_AF = Σ(ALT_F+ALT_R) / Σ DP</c> across control samples.
    /// </param>
    /// <param name="ControlSampleCount">Number of control samples observed at this locus (<c>N_SAMPLES</c>).</param>
    /// <param name="ControlSamplesWithSignal">Number of control samples with ≥ 1 alt read (<c>N_SAMPLES_WITH_SIGNAL</c>).</param>
    /// <param name="LocusNoisePass">
    /// <c>true</c> when the locus is NOT recurrently noisy: signal appears in fewer than
    /// <c>controlProportion</c> of control samples AND the pooled background AF is below
    /// <c>maxBackgroundAlleleFrequency</c> (INVAR <c>LOCUS_NOISE.PASS</c>).
    /// </param>
    public readonly record struct LocusErrorRate(
        double BackgroundErrorRate,
        int ControlSampleCount,
        int ControlSamplesWithSignal,
        bool LocusNoisePass);

    /// <summary>
    /// Estimates a tracked locus's per-locus background error rate from a panel of control (non-cancer) plasma
    /// samples and computes the INVAR locus-noise (blacklist) verdict. The background error is the pooled control
    /// allele fraction; a locus fails the noise filter when alt reads recur in too many control samples or the
    /// pooled background AF is too high. Source: INVAR2 <c>createLociErrorRateTable</c>
    /// (R/1_parse/onTargetErrorRatesAndFilter.R): <c>BACKGROUND_AF = Σ(ALT_F+ALT_R)/Σ DP</c>;
    /// <c>LOCUS_NOISE.PASS = (N_SAMPLES_WITH_SIGNAL / N_SAMPLES) &lt; proportion_of_controls AND
    /// BACKGROUND_AF &lt; max_background_mean_AF</c>.
    /// </summary>
    /// <param name="controlObservations">Per-control-sample observations of the locus (non-empty).</param>
    /// <param name="controlProportion">Max fraction of control samples with signal before blacklisting (default 0.1).</param>
    /// <param name="maxBackgroundAlleleFrequency">Max pooled control background AF (default 0.01).</param>
    /// <returns>The estimated background rate and locus-noise verdict.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="controlObservations"/> is null.</exception>
    /// <exception cref="ArgumentException">No control observation supplied.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold is out of range.</exception>
    public static LocusErrorRate EstimateLocusBackground(
        IEnumerable<ControlLocusObservation> controlObservations,
        double controlProportion = InvarDefaultControlProportion,
        double maxBackgroundAlleleFrequency = InvarDefaultMaxBackgroundAlleleFrequency)
    {
        ArgumentNullException.ThrowIfNull(controlObservations);
        if (controlProportion <= 0.0 || controlProportion > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(controlProportion), controlProportion, "Control proportion must be in (0, 1].");
        }

        if (maxBackgroundAlleleFrequency <= 0.0 || maxBackgroundAlleleFrequency > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(maxBackgroundAlleleFrequency),
                maxBackgroundAlleleFrequency,
                "Maximum background allele frequency must be in (0, 1].");
        }

        long altSum = 0;
        long depthSum = 0;
        int sampleCount = 0;
        int samplesWithSignal = 0;
        foreach (ControlLocusObservation obs in controlObservations)
        {
            int alt = Math.Max(0, obs.AltForwardReads) + Math.Max(0, obs.AltReverseReads);
            int depth = Math.Max(0, obs.TotalReads);
            altSum += alt;
            depthSum += depth;
            sampleCount++;
            if (alt > 0)
            {
                samplesWithSignal++;
            }
        }

        if (sampleCount == 0)
        {
            throw new ArgumentException(
                "At least one control observation is required to estimate background.", nameof(controlObservations));
        }

        double backgroundAf = depthSum == 0 ? 0.0 : (double)altSum / depthSum;
        double signalFraction = (double)samplesWithSignal / sampleCount;
        bool locusNoisePass = signalFraction < controlProportion && backgroundAf < maxBackgroundAlleleFrequency;

        return new LocusErrorRate(backgroundAf, sampleCount, samplesWithSignal, locusNoisePass);
    }

    /// <summary>
    /// INVAR both-strands filter for a tracked locus: a locus passes when the alt allele is observed on BOTH the
    /// forward and reverse strands, or when there is no alt signal at all. Source: INVAR2
    /// <c>onTargetErrorRatesAndFilter.R</c> — <c>BOTH_STRANDS.PASS = ALT_F &gt; 0 &amp; ALT_R &gt; 0 | AF == 0</c>.
    /// </summary>
    /// <param name="altForwardReads">Alt reads on the forward strand (≥ 0).</param>
    /// <param name="altReverseReads">Alt reads on the reverse strand (≥ 0).</param>
    /// <returns><c>true</c> when the locus passes the both-strands filter.</returns>
    public static bool PassesBothStrandsFilter(int altForwardReads, int altReverseReads)
    {
        int f = Math.Max(0, altForwardReads);
        int r = Math.Max(0, altReverseReads);

        // AF == 0 ⟺ no alt reads on either strand.
        if (f == 0 && r == 0)
        {
            return true;
        }

        return f > 0 && r > 0;
    }

    #endregion

}
