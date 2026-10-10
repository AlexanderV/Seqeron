namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region CallSomaticMutations

    /// <summary>
    /// Classifies each tumor variant as somatic, germline, or not-detected by comparing the tumor and
    /// matched-normal allele frequencies. A variant is <see cref="SomaticStatus.Somatic"/> when its tumor
    /// VAF is at or above <paramref name="tumorVafThreshold"/> (present in tumor) and its matched-normal
    /// VAF is at or below <paramref name="normalVafThreshold"/> (absent from the normal, i.e. the normal is
    /// homozygous reference). This realizes the somatic state S = {(f_t, f_n): f_t ≠ f_n} restricted to a
    /// ref/ref normal genotype (Saunders et al. 2012; Kim et al. 2018).
    /// <para>
    /// This is a deterministic VAF-threshold rule (tumor limit of detection τ_t, normal absence ceiling τ_n),
    /// not a probabilistic caller: it does not weigh read depth or base quality. The published Mutect2
    /// likelihood model (TLOD/NLOD) is available as
    /// <see cref="CallSomaticMutationsMutect2(IEnumerable{VariantObservation}, int, bool, double, double)"/>.
    /// </para>
    /// </summary>
    /// <param name="variants">Tumor variant observations with matched-normal read evidence.</param>
    /// <param name="tumorVafThreshold">Minimum tumor VAF for presence (default <see cref="DefaultTumorVafThreshold"/>).</param>
    /// <param name="normalVafThreshold">Maximum normal VAF for absence (default <see cref="DefaultNormalVafThreshold"/>).</param>
    /// <returns>One <see cref="SomaticCall"/> per input variant, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A threshold is outside [0, 1].</exception>
    public static IReadOnlyList<SomaticCall> CallSomaticMutations(
        IEnumerable<VariantObservation> variants,
        double tumorVafThreshold = DefaultTumorVafThreshold,
        double normalVafThreshold = DefaultNormalVafThreshold)
    {
        ArgumentNullException.ThrowIfNull(variants);
        ValidateThreshold(tumorVafThreshold, nameof(tumorVafThreshold));
        ValidateThreshold(normalVafThreshold, nameof(normalVafThreshold));

        var calls = new List<SomaticCall>();
        foreach (var variant in variants)
        {
            calls.Add(Classify(variant, tumorVafThreshold, normalVafThreshold));
        }

        return calls;
    }

    /// <summary>
    /// Classifies a single variant. See <see cref="CallSomaticMutations(IEnumerable{VariantObservation}, double, double)"/>.
    /// </summary>
    public static SomaticCall Classify(
        VariantObservation variant,
        double tumorVafThreshold = DefaultTumorVafThreshold,
        double normalVafThreshold = DefaultNormalVafThreshold)
    {
        ValidateThreshold(tumorVafThreshold, nameof(tumorVafThreshold));
        ValidateThreshold(normalVafThreshold, nameof(normalVafThreshold));

        double tumorVaf = CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
        double normalVaf = CalculateVaf(variant.NormalAltReads, variant.NormalTotalReads);

        SomaticStatus status;
        if (tumorVaf < tumorVafThreshold)
        {
            // Not present in the tumor above the detection limit (Yan et al. 2021).
            status = SomaticStatus.NotDetected;
        }
        else if (normalVaf <= normalVafThreshold)
        {
            // Present in tumor, absent in normal (ref/ref) → somatic (Saunders et al. 2012).
            status = SomaticStatus.Somatic;
        }
        else
        {
            // Present in both tumor and normal → germline (Benjamin et al. 2019, Mutect2).
            status = SomaticStatus.Germline;
        }

        double score = status == SomaticStatus.Somatic
            ? CalculateSomaticScore(tumorVaf, normalVaf)
            : 0.0;

        return new SomaticCall(variant, tumorVaf, normalVaf, status, score);
    }

    #endregion


    #region FilterGermlineVariants

    /// <summary>
    /// Removes germline variants, returning only somatic calls. A variant is filtered out as germline
    /// when its matched-normal VAF exceeds <paramref name="normalVafThreshold"/>, mirroring Mutect2's rule
    /// of skipping variants clearly present in the matched normal (Benjamin et al. 2019). Variants not
    /// detected in the tumor (below <paramref name="tumorVafThreshold"/>) are also excluded.
    /// </summary>
    /// <param name="variants">Tumor variant observations with matched-normal read evidence.</param>
    /// <param name="tumorVafThreshold">Minimum tumor VAF for presence (default <see cref="DefaultTumorVafThreshold"/>).</param>
    /// <param name="normalVafThreshold">Maximum normal VAF for absence (default <see cref="DefaultNormalVafThreshold"/>).</param>
    /// <returns>The subset of variants classified as <see cref="SomaticStatus.Somatic"/>, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    public static IReadOnlyList<SomaticCall> FilterGermlineVariants(
        IEnumerable<VariantObservation> variants,
        double tumorVafThreshold = DefaultTumorVafThreshold,
        double normalVafThreshold = DefaultNormalVafThreshold)
    {
        ArgumentNullException.ThrowIfNull(variants);

        return CallSomaticMutations(variants, tumorVafThreshold, normalVafThreshold)
            .Where(c => c.Status == SomaticStatus.Somatic)
            .ToList();
    }

    #endregion


    #region CalculateSomaticScore

    /// <summary>
    /// Computes a deterministic somatic-confidence score in [0, 1] from a variant's read evidence. The
    /// score increases with tumor presence (f_t) and with the separation f_t − f_n between tumor and normal
    /// allele frequencies, reflecting the somatic criterion that f_t ≠ f_n with the normal at ref/ref
    /// (Saunders et al. 2012). When the normal VAF equals or exceeds the tumor VAF the score is 0
    /// (no somatic evidence).
    /// </summary>
    /// <param name="variant">The variant observation to score.</param>
    /// <returns>Somatic confidence in [0, 1].</returns>
    public static double CalculateSomaticScore(VariantObservation variant)
    {
        double tumorVaf = CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
        double normalVaf = CalculateVaf(variant.NormalAltReads, variant.NormalTotalReads);
        return CalculateSomaticScore(tumorVaf, normalVaf);
    }

    /// <summary>
    /// Computes the somatic-confidence score from tumor and normal allele frequencies:
    /// score = max(0, f_t − f_n), bounded in [0, 1] because both VAFs are in [0, 1].
    /// </summary>
    private static double CalculateSomaticScore(double tumorVaf, double normalVaf)
    {
        // Separation between tumor and normal allele frequencies; 0 when the normal carries the allele
        // at or above the tumor level (no somatic signal). Bounded in [0, 1] since both VAFs are in [0, 1].
        double separation = tumorVaf - normalVaf;
        return separation > 0.0 ? separation : 0.0;
    }

    #endregion


    #region Mutect2 somatic likelihoods model (count-based pileup form)

    /// <summary>
    /// Mutect2 default <c>--tumor-lod-to-emit</c>: an alt allele is emitted only when its tumor log10 odds
    /// (TLOD) exceeds 3.0. Source: GATK <c>M2ArgumentCollection.DEFAULT_EMISSION_LOG_10_ODDS = 3.0</c>
    /// (broadinstitute/gatk, tools/walkers/mutect/M2ArgumentCollection.java).
    /// </summary>
    public const double DefaultMutect2TumorLog10OddsThreshold = 3.0;

    /// <summary>
    /// Mutect2 default <c>--normal-lod</c>: with a matched normal, an allele is genotyped as somatic only
    /// when the normal log10 odds of hom-ref versus het (NLOD) exceeds 2.2; otherwise it is skipped as
    /// germline. Source: GATK <c>M2ArgumentCollection.DEFAULT_NORMAL_LOG_10_ODDS = 2.2</c>.
    /// </summary>
    public const double DefaultMutect2NormalLog10OddsThreshold = 2.2;

    /// <summary>
    /// Mutect2 default <c>--phred-scaled-global-read-mismapping-rate</c> (45): each read's likelihood under
    /// any allele is floored at (best allele likelihood) × 10^(−45/10). Source: GATK
    /// <c>LikelihoodEngineArgumentCollection.phredScaledGlobalReadMismappingRate = 45</c>, applied by
    /// <c>SomaticGenotypingEngine</c> through <c>AlleleLikelihoods.normalizeLikelihoods</c>.
    /// </summary>
    private const int Mutect2GlobalReadMismappingRatePhred = 45;

    /// <summary>
    /// Pair-HMM tri-state correction: a mismatching base of quality Q has likelihood ε/3 (ε = 10^(−Q/10)),
    /// the error spread over the three other bases. Source: GATK <c>PairHMM.TRISTATE_CORRECTION = 3.0</c>.
    /// </summary>
    private const double PairHmmTristateCorrection = 3.0;

    /// <summary>
    /// Convergence threshold of the Dirichlet mean-field iteration (L1 change / total &lt; 0.001).
    /// Source: GATK <c>SomaticLikelihoodsEngine.CONVERGENCE_THRESHOLD = 0.001</c>.
    /// </summary>
    private const double Mutect2ConvergenceThreshold = 0.001;

    /// <summary>Responsibility below which a likelihood term is dropped (GATK <c>NEGLIGIBLE_RESPONSIBILITY = 1e-10</c>).</summary>
    private const double Mutect2NegligibleResponsibility = 1.0e-10;

    /// <summary>x·ln x is treated as 0 below this x (GATK <c>SomaticLikelihoodsEngine.xLogx</c> cut-off 1e-8).</summary>
    private const double Mutect2EntropyCutoff = 1.0e-8;

    /// <summary>
    /// Result of the Mutect2 somatic likelihoods model for one variant.
    /// </summary>
    /// <param name="Variant">The variant that was classified.</param>
    /// <param name="TumorVaf">Empirical tumor VAF f_t = altReads / totalReads.</param>
    /// <param name="NormalVaf">Empirical matched-normal VAF f_n (0 when uncovered).</param>
    /// <param name="TumorLog10Odds">Mutect2 TLOD: log10 evidence ratio of {ref, alt} vs {ref} in the tumor.</param>
    /// <param name="NormalLog10Odds">
    /// Mutect2 NLOD: log10 odds that the normal is hom-ref rather than het (NaN when no matched normal).
    /// </param>
    /// <param name="Status">Somatic (TLOD &gt; τ_T and, with a normal, NLOD &gt; τ_N); Germline (TLOD &gt; τ_T, NLOD ≤ τ_N); NotDetected (TLOD ≤ τ_T).</param>
    public readonly record struct Mutect2SomaticCall(
        VariantObservation Variant,
        double TumorVaf,
        double NormalVaf,
        double TumorLog10Odds,
        double NormalLog10Odds,
        SomaticStatus Status);

    /// <summary>
    /// Classifies tumor variants with the GATK Mutect2 somatic likelihoods model evaluated on a biallelic
    /// pileup (Benjamin et al. 2019, <c>docs/mutect/mutect.tex</c>; <c>SomaticGenotypingEngine</c>):
    /// <list type="number">
    /// <item><description>Read likelihoods from base quality Q (Pair-HMM form): a read matching an allele has
    /// likelihood 1−ε, a mismatching read ε/3, ε = 10^(−Q/10); each read's worse likelihood is floored at the
    /// best × 10^(−4.5) (global mismapping rate Q45).</description></item>
    /// <item><description>TLOD = [ln P(R | {ref, alt}) − ln P(R | {ref})] / ln 10, where the evidence is the
    /// variational (mean-field Dirichlet, flat prior α = (1, 1)) model evidence of
    /// <c>SomaticLikelihoodsEngine.logEvidence</c>.</description></item>
    /// <item><description>NLOD = [Σ_r ln ℓ_r,ref − Σ_r (ln(ℓ_r,ref + ℓ_r,alt) + ln ½)] / ln 10 — hom-ref vs
    /// het in the normal (<c>diploidAltLogOdds</c>).</description></item>
    /// <item><description>Emit when TLOD &gt; <paramref name="tumorLog10OddsThreshold"/>; with a matched normal
    /// the allele is somatic only if NLOD &gt; <paramref name="normalLog10OddsThreshold"/>, otherwise it is
    /// skipped as germline. Tumor-only (<paramref name="hasMatchedNormal"/> = false) skips the normal test.</description></item>
    /// </list>
    /// Read counts are the only evidence available in <see cref="VariantObservation"/>, so every read is
    /// assigned the same base quality <paramref name="baseQuality"/> and ref reads = total − alt (biallelic
    /// site). Because the variant base is the only position at which the two haplotypes differ, the shared
    /// likelihood of the remaining read bases cancels in both log-odds, so the count form reproduces
    /// Mutect2's per-read computation exactly for uniform base quality. Not modelled here (they act on the
    /// read level or in <c>FilterMutectCalls</c>): local assembly, fragment merging / PCR-quality capping of
    /// overlapping mates, germline-resource population priors and the post-call filters.
    /// </summary>
    /// <param name="variants">Tumor variant observations with matched-normal read evidence.</param>
    /// <param name="baseQuality">Phred base quality Q (≥ 1) assigned to every read.</param>
    /// <param name="hasMatchedNormal">
    /// True for tumor–normal mode (the NLOD test is applied even when the normal has zero coverage, where
    /// NLOD = 0 and the allele is skipped as in Mutect2); false for tumor-only mode (ℓ_n = 1).
    /// </param>
    /// <param name="tumorLog10OddsThreshold">TLOD emission threshold (default 3.0).</param>
    /// <param name="normalLog10OddsThreshold">NLOD threshold (default 2.2).</param>
    /// <returns>One <see cref="Mutect2SomaticCall"/> per input variant, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="baseQuality"/> &lt; 1, a threshold is NaN/infinite, or a variant has invalid read counts.
    /// </exception>
    public static IReadOnlyList<Mutect2SomaticCall> CallSomaticMutationsMutect2(
        IEnumerable<VariantObservation> variants,
        int baseQuality,
        bool hasMatchedNormal = true,
        double tumorLog10OddsThreshold = DefaultMutect2TumorLog10OddsThreshold,
        double normalLog10OddsThreshold = DefaultMutect2NormalLog10OddsThreshold)
    {
        ArgumentNullException.ThrowIfNull(variants);
        ValidateBaseQuality(baseQuality);
        ValidateFiniteThreshold(tumorLog10OddsThreshold, nameof(tumorLog10OddsThreshold));
        ValidateFiniteThreshold(normalLog10OddsThreshold, nameof(normalLog10OddsThreshold));

        var calls = new List<Mutect2SomaticCall>();
        foreach (VariantObservation variant in variants)
        {
            double tumorVaf = CalculateVaf(variant.TumorAltReads, variant.TumorTotalReads);
            double normalVaf = CalculateVaf(variant.NormalAltReads, variant.NormalTotalReads);

            double tlod = CalculateMutect2TumorLog10Odds(
                variant.TumorTotalReads - variant.TumorAltReads, variant.TumorAltReads, baseQuality);
            double nlod = hasMatchedNormal
                ? CalculateMutect2NormalLog10Odds(
                    variant.NormalTotalReads - variant.NormalAltReads, variant.NormalAltReads, baseQuality)
                : double.NaN;

            SomaticStatus status;
            if (!(tlod > tumorLog10OddsThreshold))
            {
                status = SomaticStatus.NotDetected;
            }
            else if (!hasMatchedNormal || nlod > normalLog10OddsThreshold)
            {
                status = SomaticStatus.Somatic;
            }
            else
            {
                status = SomaticStatus.Germline;
            }

            calls.Add(new Mutect2SomaticCall(variant, tumorVaf, normalVaf, tlod, nlod, status));
        }

        return calls;
    }

    /// <summary>
    /// Mutect2 tumor log10 odds (TLOD) that an alt allele exists, for a biallelic pileup of
    /// <paramref name="refReads"/> ref and <paramref name="altReads"/> alt reads at uniform base quality:
    /// the log10 ratio of the variational model evidence with alleles {ref, alt} to that with {ref} only
    /// (GATK <c>SomaticGenotypingEngine.somaticLogOdds</c> / <c>SomaticLikelihoodsEngine.logEvidence</c>,
    /// flat Dirichlet prior (1, 1)). Returns 0 for an empty pileup (GATK: evidence 0 with no reads).
    /// </summary>
    /// <param name="refReads">Reads supporting the reference allele (≥ 0).</param>
    /// <param name="altReads">Reads supporting the alternate allele (≥ 0).</param>
    /// <param name="baseQuality">Phred base quality Q (≥ 1).</param>
    /// <returns>TLOD in log10 units.</returns>
    public static double CalculateMutect2TumorLog10Odds(int refReads, int altReads, int baseQuality)
    {
        ValidateReadCounts(refReads, altReads);
        ValidateBaseQuality(baseQuality);
        if (refReads + altReads == 0)
        {
            return 0.0;
        }

        Mutect2ReadLogLikelihoods(baseQuality, out double logMatch, out double logMismatch);

        // Two read classes: ref reads (ln ℓ_ref, ln ℓ_alt) = (match, mismatch); alt reads = (mismatch, match).
        double logEvidenceBoth = Mutect2LogEvidenceBiallelic(refReads, altReads, logMatch, logMismatch);

        // Ref-only allele set: a single allele with pseudocount 1 ⇒ g(α) = g(β) = 0 and every responsibility 1,
        // so the evidence is Σ_r ln ℓ_r,ref.
        double logEvidenceRefOnly = refReads * logMatch + altReads * logMismatch;

        return (logEvidenceBoth - logEvidenceRefOnly) / Math.Log(10.0);
    }

    /// <summary>
    /// Mutect2 normal log10 odds (NLOD) that the normal is homozygous reference rather than heterozygous:
    /// NLOD = [Σ_r ln ℓ_r,ref − Σ_r (logsumexp(ln ℓ_r,ref, ln ℓ_r,alt) + ln ½)] / ln 10
    /// (GATK <c>SomaticGenotypingEngine.diploidAltLogOdds</c>). Returns 0 for an empty pileup.
    /// </summary>
    /// <param name="refReads">Normal reads supporting the reference allele (≥ 0).</param>
    /// <param name="altReads">Normal reads supporting the alternate allele (≥ 0).</param>
    /// <param name="baseQuality">Phred base quality Q (≥ 1).</param>
    /// <returns>NLOD in log10 units.</returns>
    public static double CalculateMutect2NormalLog10Odds(int refReads, int altReads, int baseQuality)
    {
        ValidateReadCounts(refReads, altReads);
        ValidateBaseQuality(baseQuality);

        Mutect2ReadLogLikelihoods(baseQuality, out double logMatch, out double logMismatch);

        double logHalf = Math.Log(0.5);
        double hetPerRead = LogSumExp(logMatch, logMismatch) + logHalf; // symmetric for ref and alt reads
        double homRef = refReads * logMatch + altReads * logMismatch;
        double het = (refReads + altReads) * hetPerRead;
        return (homRef - het) / Math.Log(10.0);
    }

    /// <summary>
    /// Per-read log-likelihoods of the matching and mismatching allele at base quality Q:
    /// ln(1−ε) and ln(ε/3), with the mismatching value floored at ln(1−ε) − 4.5·ln 10 (global mismapping Q45).
    /// </summary>
    private static void Mutect2ReadLogLikelihoods(int baseQuality, out double logMatch, out double logMismatch)
    {
        double epsilon = Math.Pow(10.0, -baseQuality / 10.0);
        double match = Math.Log(1.0 - epsilon);
        double mismatch = Math.Log(epsilon / PairHmmTristateCorrection);

        // AlleleLikelihoods.normalizeLikelihoods: every allele's likelihood ≥ best + ln(10^(−45/10)).
        double best = Math.Max(match, mismatch);
        double cap = best - Mutect2GlobalReadMismappingRatePhred / 10.0 * Math.Log(10.0);
        logMatch = Math.Max(match, cap);
        logMismatch = Math.Max(mismatch, cap);
    }

    /// <summary>
    /// Variational log model evidence ln P(R | {ref, alt}) of GATK <c>SomaticLikelihoodsEngine.logEvidence</c>
    /// for two read classes (n_ref ref reads, n_alt alt reads) with a flat Dirichlet prior α = (1, 1):
    /// iterate β = α + Σ_r z̄_r with z̄_ra ∝ exp(ψ(β_a) − ψ(Σβ))·ℓ_ra from β = (1, 1) until
    /// ‖β_new − β‖₁ / Σβ_new &lt; 0.001, then
    /// ln P = g(α) − g(β) + Σ_r Σ_a z̄_ra (ln ℓ_ra − ln z̄_ra), g(ω) = lnΓ(Σω) − Σ lnΓ(ω_a).
    /// Reads of one class share identical responsibilities, so the per-read sums are class count × value.
    /// </summary>
    private static double Mutect2LogEvidenceBiallelic(int refReads, int altReads, double logMatch, double logMismatch)
    {
        // Allele order (ref, alt). Ref-read column (logMatch, logMismatch); alt-read column (logMismatch, logMatch).
        double[] prior = { 1.0, 1.0 };
        double[] posterior = { 1.0, 1.0 };
        bool converged = false;
        while (!converged)
        {
            Mutect2Responsibilities(posterior, logMatch, logMismatch, out double refReadRef, out double refReadAlt,
                out double altReadRef, out double altReadAlt);
            double[] updated =
            {
                refReads * refReadRef + altReads * altReadRef + prior[0],
                refReads * refReadAlt + altReads * altReadAlt + prior[1],
            };
            double distance = Math.Abs(posterior[0] - updated[0]) + Math.Abs(posterior[1] - updated[1]);
            converged = distance / (updated[0] + updated[1]) < Mutect2ConvergenceThreshold;
            posterior = updated;
        }

        Mutect2Responsibilities(posterior, logMatch, logMismatch, out double rr, out double ra,
            out double ar, out double aa);

        double priorContribution = LogDirichletNormalization(prior);
        double posteriorContribution = -LogDirichletNormalization(posterior);
        double refReadTerm = LikelihoodTerm(logMatch, rr) + LikelihoodTerm(logMismatch, ra) - XLogX(rr) - XLogX(ra);
        double altReadTerm = LikelihoodTerm(logMismatch, ar) + LikelihoodTerm(logMatch, aa) - XLogX(ar) - XLogX(aa);
        return priorContribution + posteriorContribution + refReads * refReadTerm + altReads * altReadTerm;
    }

    /// <summary>
    /// Mean-field responsibilities z̄_ra = softmax_a(ψ(β_a) − ψ(Σβ) + ln ℓ_ra) for a ref read and an alt read.
    /// </summary>
    private static void Mutect2Responsibilities(
        double[] posterior, double logMatch, double logMismatch,
        out double refReadRef, out double refReadAlt, out double altReadRef, out double altReadAlt)
    {
        double digammaSum = Digamma(posterior[0] + posterior[1]);
        double logWeightRef = Digamma(posterior[0]) - digammaSum;
        double logWeightAlt = Digamma(posterior[1]) - digammaSum;

        (refReadRef, refReadAlt) = Softmax2(logWeightRef + logMatch, logWeightAlt + logMismatch);
        (altReadRef, altReadAlt) = Softmax2(logWeightRef + logMismatch, logWeightAlt + logMatch);
    }

    private static (double First, double Second) Softmax2(double a, double b)
    {
        double norm = LogSumExp(a, b);
        return (Math.Exp(a - norm), Math.Exp(b - norm));
    }

    /// <summary>ln ℓ · z̄, dropped when z̄ is negligible (GATK <c>likelihoodsContribution</c>).</summary>
    private static double LikelihoodTerm(double logLikelihood, double responsibility)
        => responsibility < Mutect2NegligibleResponsibility ? 0.0 : logLikelihood * responsibility;

    /// <summary>x ln x with the GATK cut-off (0 for x &lt; 1e-8).</summary>
    private static double XLogX(double x) => x < Mutect2EntropyCutoff ? 0.0 : x * Math.Log(x);

    /// <summary>Dirichlet log normalization g(ω) = lnΓ(Σω) − Σ lnΓ(ω_a).</summary>
    private static double LogDirichletNormalization(double[] parameters)
    {
        double sum = 0.0;
        double sumLogGamma = 0.0;
        foreach (double p in parameters)
        {
            sum += p;
            sumLogGamma += LogGamma(p);
        }

        return LogGamma(sum) - sumLogGamma;
    }

    private static double LogSumExp(double a, double b)
    {
        double max = Math.Max(a, b);
        return max + Math.Log(Math.Exp(a - max) + Math.Exp(b - max));
    }

    /// <summary>
    /// Digamma ψ(x) for x &gt; 0: upward recurrence ψ(x) = ψ(x+1) − 1/x to x ≥ 6, then the asymptotic series
    /// ψ(x) ≈ ln x − 1/(2x) − 1/(12x²) + 1/(120x⁴) − 1/(252x⁶) + 1/(240x⁸) − 1/(132x¹⁰)
    /// (Abramowitz &amp; Stegun 6.3.18); absolute error &lt; 1e-12 in this range.
    /// </summary>
    private static double Digamma(double x)
    {
        double result = 0.0;
        while (x < 6.0)
        {
            result -= 1.0 / x;
            x += 1.0;
        }

        double inv = 1.0 / x;
        double inv2 = inv * inv;
        result += Math.Log(x) - 0.5 * inv
            - inv2 * (1.0 / 12.0 - inv2 * (1.0 / 120.0 - inv2 * (1.0 / 252.0 - inv2 * (1.0 / 240.0 - inv2 / 132.0))));
        return result;
    }

    private static void ValidateBaseQuality(int baseQuality)
    {
        if (baseQuality < 1)
        {
            throw new ArgumentOutOfRangeException(nameof(baseQuality), baseQuality, "Phred base quality must be at least 1.");
        }
    }

    private static void ValidateReadCounts(int refReads, int altReads)
    {
        if (refReads < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(refReads), refReads, "Ref read count cannot be negative.");
        }

        if (altReads < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(altReads), altReads, "Alt read count cannot be negative.");
        }
    }

    private static void ValidateFiniteThreshold(double value, string paramName)
    {
        if (!double.IsFinite(value))
        {
            throw new ArgumentOutOfRangeException(paramName, value, "Log-odds threshold must be finite.");
        }
    }

    #endregion


    #region Variant Allele Frequency

    /// <summary>
    /// Computes the empirical variant allele frequency (VAF) at a locus as the fraction of covering reads
    /// that support the alternate allele: VAF = altReads / totalReads. This is the model-free allele
    /// fraction derived from per-allele read depths (GATK <c>AD</c> field: alt AD / Σ AD; samtools read
    /// counts), distinct from Mutect2's Bayesian <c>AF</c> estimate. A site with no coverage
    /// (totalReads == 0) returns 0, since an uncovered site provides no evidence of the allele.
    /// </summary>
    /// <param name="altReads">Reads supporting the alternate allele (≥ 0, ≤ <paramref name="totalReads"/>).</param>
    /// <param name="totalReads">Total covering reads at the locus (≥ 0).</param>
    /// <returns>VAF in [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">A count is negative, or altReads &gt; totalReads.</exception>
    public static double CalculateVAF(int altReads, int totalReads) => CalculateVaf(altReads, totalReads);

    /// <summary>
    /// Computes the empirical VAF and its Wilson score confidence interval for the underlying allele
    /// proportion. The Wilson score interval (Wilson 1927) for n trials with p̂ = altReads/totalReads is:
    /// <code>
    /// center = (p̂ + z²/(2n)) / (1 + z²/n)
    /// margin = (z / (1 + z²/n)) · √( p̂(1−p̂)/n + z²/(4n²) )
    /// interval = center ± margin
    /// </code>
    /// where z is the standard-normal quantile for the requested confidence (z = 1.96 for 95%). The
    /// interval is bounded within [0, 1] (no overshoot) and has non-zero width even at p̂ = 0 or 1,
    /// unlike the Wald interval. Source: Wilson E.B. (1927), JASA 22(158):209–212, via the Binomial
    /// proportion confidence interval specification.
    /// </summary>
    /// <param name="altReads">Reads supporting the alternate allele (≥ 0, ≤ <paramref name="totalReads"/>).</param>
    /// <param name="totalReads">Total covering reads at the locus (&gt; 0 for a defined interval).</param>
    /// <param name="confidence">Two-sided confidence level in (0, 1); default 0.95 (z = 1.96).</param>
    /// <returns>The VAF point estimate with its Wilson score interval.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// A read count is invalid (see <see cref="CalculateVAF"/>), totalReads == 0, or confidence is outside (0, 1)
    /// at a level other than the supported 0.95.
    /// </exception>
    public static VafConfidenceInterval CalculateVAFConfidenceInterval(
        int altReads,
        int totalReads,
        double confidence = DefaultVafConfidence)
    {
        double vaf = CalculateVaf(altReads, totalReads);

        if (totalReads == 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(totalReads), "A confidence interval is undefined with zero coverage.");
        }

        double z = ZScoreFor(confidence);

        double n = totalReads;
        double pHat = vaf;
        double z2 = z * z;
        double denominator = 1.0 + z2 / n;
        double center = (pHat + z2 / (2.0 * n)) / denominator;
        double margin = (z / denominator) * Math.Sqrt(pHat * (1.0 - pHat) / n + z2 / (4.0 * n * n));

        // The Wilson interval is mathematically within [0, 1]; clamp guards against floating-point drift
        // at the exact boundaries p̂ = 0 (lower = 0) and p̂ = 1 (upper = 1).
        double lower = Math.Max(0.0, center - margin);
        double upper = Math.Min(1.0, center + margin);

        return new VafConfidenceInterval(vaf, lower, upper, confidence);
    }

    /// <summary>
    /// Adjusts an observed VAF for tumor purity and tumor-segment ploidy, recovering the per-tumour-copy
    /// mutant fraction (multiplicity × cancer cell fraction). Inverting the CNAqc expected-VAF relation
    /// v = (m·π) / (2(1−π) + π·n_tot) gives:
    /// <code>
    /// adjusted = vaf · (2(1−π) + π·ploidy) / π
    /// </code>
    /// where π = purity, ploidy = tumor total copy number n_tot, and the normal contribution is fixed at
    /// 2 (autosomal diploid). For a heterozygous somatic SNV in a diploid tumor (ploidy = 2) this reduces
    /// to adjusted = vaf / (π/2): e.g. observed VAF 0.4 at purity 0.8 ⇒ 1.0. Source: CNAqc
    /// (Genome Biology 2024); Tarabichi et al. (2017), PMC5538405.
    /// </summary>
    /// <param name="vaf">Observed VAF in [0, 1].</param>
    /// <param name="purity">Tumor purity π in (0, 1].</param>
    /// <param name="ploidy">Tumor total copy number n_tot at the locus (&gt; 0); 2 for a diploid segment.</param>
    /// <returns>The purity/ploidy-corrected mutant fraction (multiplicity × CCF).</returns>
    /// <exception cref="ArgumentOutOfRangeException">vaf ∉ [0, 1], purity ∉ (0, 1], or ploidy ≤ 0.</exception>
    public static double AdjustVAFForPurity(double vaf, double purity, double ploidy)
    {
        if (double.IsNaN(vaf) || vaf < 0.0 || vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(vaf), vaf, "VAF must be in the range [0, 1].");
        }

        if (double.IsNaN(purity) || purity <= 0.0 || purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(purity), purity, "Purity must be in the range (0, 1]; correction divides by purity.");
        }

        if (double.IsNaN(ploidy) || ploidy <= 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(ploidy), ploidy, "Ploidy (tumor total copy number) must be positive.");
        }

        // Weighted average total copies per cell: 2(1−π) from normal cells + π·ploidy from tumor cells.
        double averageCopiesPerCell = NormalDiploidCopyNumber * (1.0 - purity) + purity * ploidy;
        return vaf * averageCopiesPerCell / purity;
    }

    /// <summary>
    /// Estimates tumor purity ρ from the variant allele frequencies of clonal, heterozygous somatic SNVs
    /// at copy-neutral diploid loci. For such a variant (multiplicity m = 1, total copy number n_tot = 2),
    /// the CNAqc expected-VAF relation v = m·π / [2(1−π) + π·n_tot] reduces to v = π/2, so the per-variant
    /// purity is ρ = 2·v (Antonello et al. 2024, <i>Genome Biology</i> 25:38; CNAqc reports purity 60% ⇔ VAF 30%).
    /// The per-variant estimates are aggregated by their median, which is robust to subclonal / outlier VAFs.
    /// </summary>
    /// <param name="variants">Observed clonal heterozygous somatic SNVs (each contributes ρ = 2·VAF).</param>
    /// <returns>Estimated tumor purity ρ ∈ [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="variants"/> is empty (purity is undefined).</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// A variant has invalid read counts, or a VAF &gt; 0.5 (which would imply purity &gt; 1 under the diploid model).
    /// </exception>
    public static double EstimatePurityFromVAF(IEnumerable<VariantObservation> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var purities = new List<double>();
        foreach (VariantObservation variant in variants)
        {
            double vaf = CalculateVAF(variant.TumorAltReads, variant.TumorTotalReads);
            purities.Add(EstimatePurityFromVaf(vaf));
        }

        if (purities.Count == 0)
        {
            throw new ArgumentException("Cannot estimate purity from an empty variant set.", nameof(variants));
        }

        return Median(purities);
    }

    /// <summary>
    /// Estimates tumor purity ρ from a single clonal heterozygous somatic SNV VAF at a copy-neutral diploid
    /// locus using the closed form ρ = 2·v (the m = 1, n_tot = 2 special case of the CNAqc expected-VAF
    /// relation; Antonello et al. 2024, <i>Genome Biology</i> 25:38).
    /// </summary>
    /// <param name="vaf">Observed VAF v ∈ [0, 0.5] of the clonal heterozygous SNV.</param>
    /// <returns>Estimated tumor purity ρ = 2·v ∈ [0, 1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// vaf ∉ [0, 1], or vaf &gt; 0.5 (which would imply purity &gt; 1 under the diploid heterozygous model).
    /// </exception>
    public static double EstimatePurityFromVaf(double vaf)
    {
        if (double.IsNaN(vaf) || vaf < 0.0 || vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(vaf), vaf, "VAF must be in the range [0, 1].");
        }

        // ρ = 2·v; a heterozygous SNV in a diploid tumour cannot exceed VAF 0.5, so ρ > 1 is impossible input.
        double purity = HeterozygousDiploidPurityFactor * vaf;
        if (purity > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(vaf), vaf,
                "VAF > 0.5 implies purity > 1 under the heterozygous diploid model; the locus is not copy-neutral diploid heterozygous.");
        }

        return purity;
    }

    /// <summary>
    /// Estimates tumor purity ρ from clonal somatic mutations with known allele-specific copy-number state by
    /// inverting the CNAqc expected-VAF relation v = m·π / [2(1−π) + π·n_tot] for π:
    /// <code>
    /// π = 2·v / [m + v·(2 − n_tot)]
    /// </code>
    /// where m is the mutation multiplicity and n_tot the tumour total copy number (Antonello et al. 2024,
    /// <i>Genome Biology</i> 25:38). The per-variant estimates are aggregated by their median.
    /// </summary>
    /// <param name="variants">Clonal somatic mutations with VAF, multiplicity m, and tumour total copy number n_tot.</param>
    /// <returns>Estimated tumor purity ρ ∈ [0, 1].</returns>
    /// <exception cref="ArgumentNullException"><paramref name="variants"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="variants"/> is empty (purity is undefined).</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// A variant has vaf ∉ [0, 1], multiplicity &lt; 1, n_tot &lt; 1, multiplicity &gt; n_tot, or yields a purity
    /// outside [0, 1] (beyond a machine-ε rounding tolerance of (n_tot + 4)·ε, which is clamped to 1).
    /// </exception>
    public static double EstimatePurity(IEnumerable<PurityVariant> variants)
    {
        ArgumentNullException.ThrowIfNull(variants);

        var purities = new List<double>();
        foreach (PurityVariant variant in variants)
        {
            purities.Add(EstimatePurityFromAlleleSpecificVaf(variant));
        }

        if (purities.Count == 0)
        {
            throw new ArgumentException("Cannot estimate purity from an empty variant set.", nameof(variants));
        }

        return Median(purities);
    }

    /// <summary>
    /// Inverts the CNAqc expected-VAF relation for a single allele-specific variant:
    /// π = 2·v / [m + v·(2 − n_tot)].
    /// </summary>
    private const double PurityInversionMachineEpsilon = 2.220446049250313e-16; // IEEE-754 double machine ε (2^-52)

    private static double EstimatePurityFromAlleleSpecificVaf(in PurityVariant variant)
    {
        if (double.IsNaN(variant.Vaf) || variant.Vaf < 0.0 || variant.Vaf > 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Vaf, "VAF must be in the range [0, 1].");
        }

        if (variant.Multiplicity < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Multiplicity, "Mutation multiplicity m must be at least 1.");
        }

        if (variant.TumorTotalCopyNumber < 1)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.TumorTotalCopyNumber, "Tumour total copy number n_tot must be at least 1.");
        }

        // CNAqc expectations_generalised enumerates multiplicities only over 1..Major ≤ n_tot.
        if (variant.Multiplicity > variant.TumorTotalCopyNumber)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Multiplicity,
                "Multiplicity m cannot exceed the tumour total copy number n_tot (CNAqc: m ∈ 1..Major).");
        }

        // π = 2v / [m + v(2 − n_tot)], the algebraic inverse of v = mπ / [2(1−π) + π·n_tot].
        double denominator = variant.Multiplicity + variant.Vaf * (NormalDiploidCopyNumber - variant.TumorTotalCopyNumber);
        if (denominator <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Vaf,
                "The (VAF, multiplicity, copy-number) combination does not correspond to a purity in [0, 1].");
        }

        double purity = NormalDiploidCopyNumber * variant.Vaf / denominator;

        // An exact clonal peak at π = 1 (v = m/n_tot) can round to 1 + a few ulp; the inversion error is
        // ≤ 0.47·(n_tot + 4)·ε for n_tot ≤ 2000, so accept up to (n_tot + 4)·ε and clamp to 1.
        double tolerance = (variant.TumorTotalCopyNumber + 4) * PurityInversionMachineEpsilon;
        if (purity < 0.0 || purity > 1.0 + tolerance)
        {
            throw new ArgumentOutOfRangeException(
                nameof(variant), variant.Vaf,
                "The (VAF, multiplicity, copy-number) combination yields a purity outside [0, 1].");
        }

        return Math.Min(1.0, purity);
    }

    /// <summary>Median of a non-empty list of values (lower-mid average for even counts). Does not mutate the input.</summary>
    private static double Median(List<double> values)
    {
        double[] sorted = values.ToArray();
        Array.Sort(sorted);
        int n = sorted.Length;
        int mid = n / 2;
        return (n % 2 == 1) ? sorted[mid] : 0.5 * (sorted[mid - 1] + sorted[mid]);
    }

    #endregion


    #region Helpers

    /// <summary>
    /// Maps a two-sided confidence level to its standard-normal quantile z. Only the source-cited 95%
    /// level (z = 1.96, Wilson 1927) is supported; other levels would require additional cited z values.
    /// </summary>
    private static double ZScoreFor(double confidence)
    {
        if (double.IsNaN(confidence) || confidence <= 0.0 || confidence >= 1.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(confidence), confidence, "Confidence must be in the open interval (0, 1).");
        }

        if (Math.Abs(confidence - DefaultVafConfidence) < 1e-12)
        {
            return ZScore95;
        }

        throw new ArgumentOutOfRangeException(
            nameof(confidence), confidence,
            "Only the 0.95 confidence level (z = 1.96) is supported by the cited source.");
    }

    /// <summary>
    /// Computes a variant allele frequency f = altReads / totalReads. Returns 0 when there is no coverage
    /// (totalReads == 0) so an uncovered site is treated as the allele being absent.
    /// </summary>
    private static double CalculateVaf(int altReads, int totalReads)
    {
        if (altReads < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(altReads), "Alt read count cannot be negative.");
        }

        if (totalReads < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(totalReads), "Total read count cannot be negative.");
        }

        if (altReads > totalReads)
        {
            throw new ArgumentOutOfRangeException(nameof(altReads), "Alt read count cannot exceed total read count.");
        }

        return totalReads == 0 ? 0.0 : (double)altReads / totalReads;
    }

    private static void ValidateThreshold(double value, string paramName)
    {
        if (double.IsNaN(value) || value < 0.0 || value > 1.0)
        {
            throw new ArgumentOutOfRangeException(paramName, value, "Threshold must be in the range [0, 1].");
        }
    }

    #endregion

}
