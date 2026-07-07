using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Oncology;

public static partial class OncologyAnalyzer
{
    #region Tumor Mutational Burden (TMB)

    /// <summary>
    /// FDA-approved tumor-mutational-burden-high (TMB-H) cutoff for pembrolizumab, in mutations per
    /// megabase. A tumor is TMB-High when its TMB is at or above this value (the threshold is inclusive).
    /// Source: Marcus L et al. (2021), FDA Approval Summary: Pembrolizumab for the Treatment of Tumor
    /// Mutational Burden–High Solid Tumors, Clin Cancer Res 27(17):4685–4689 — pembrolizumab approved
    /// (June 16, 2020) for solid tumors with "TMB ≥10 mutations/megabase (mut/Mb)", companion diagnostic
    /// FoundationOne CDx. Value = 10.0 mut/Mb.
    /// </summary>
    public const double TmbHighThreshold = 10.0;

    /// <summary>
    /// TMB-High classification result. Only the FDA/F1CDx TMB-High cutoff (≥ 10 mut/Mb) is source-backed,
    /// so the classification is two-tier (High vs Low); intermediate research cut-points (e.g. tumor-type
    /// specific 6/20 boundaries) are not implemented because no authoritative source defines them.
    /// Source: Marcus et al. (2021).
    /// </summary>
    public enum TmbStatus
    {
        /// <summary>TMB below the FDA TMB-High cutoff (TMB &lt; 10 mut/Mb).</summary>
        Low,

        /// <summary>TMB at or above the FDA TMB-High cutoff (TMB ≥ 10 mut/Mb).</summary>
        High
    }

    /// <summary>
    /// Computes the tumor mutational burden (TMB) as the number of somatic coding mutations per megabase of
    /// sequenced region: TMB = <paramref name="mutationCount"/> / <paramref name="targetRegionMb"/>, in
    /// mutations/megabase (mut/Mb). Source: Chalmers ZR et al. (2017), Genome Medicine 9:34 — "TMB was
    /// defined as the number of somatic, coding, base substitution, and indel mutations per megabase of
    /// genome examined" (the FoundationOne 315-gene panel denominator is 1.1 Mb of coding genome). The
    /// caller supplies the already-filtered somatic mutation count (germline / known-driver filtering is the
    /// upstream somatic caller's responsibility, per Chalmers 2017's pre-count filtering and ONCO-SOMATIC-001).
    /// </summary>
    /// <param name="mutationCount">Number of counted somatic mutations (≥ 0).</param>
    /// <param name="targetRegionMb">Size of the sequenced coding region in megabases (&gt; 0, finite).</param>
    /// <returns>TMB in mutations per megabase (≥ 0).</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="mutationCount"/> is negative, or <paramref name="targetRegionMb"/> is not a finite
    /// value greater than 0 (TMB is undefined when the megabase denominator is 0).
    /// </exception>
    public static double CalculateTMB(int mutationCount, double targetRegionMb)
    {
        if (mutationCount < 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(mutationCount), mutationCount, "Mutation count cannot be negative.");
        }

        if (double.IsNaN(targetRegionMb) || double.IsInfinity(targetRegionMb) || targetRegionMb <= 0.0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(targetRegionMb), targetRegionMb,
                "Target region size (Mb) must be a finite value greater than 0; TMB is undefined at 0 Mb.");
        }

        // TMB = mutations / megabase (Chalmers et al. 2017).
        return mutationCount / targetRegionMb;
    }

    /// <summary>
    /// Computes TMB from a set of classified somatic calls (e.g. from <see cref="CallSomaticMutations"/>):
    /// counts the calls with <see cref="SomaticStatus.Somatic"/> status and divides by
    /// <paramref name="targetRegionMb"/>. Germline and not-detected calls are excluded, matching the TMB
    /// definition of counting only somatic mutations (Chalmers et al. 2017; Friends of Cancer Research TMB
    /// Harmonization, Merino et al. 2020). Thin wrapper over
    /// <see cref="CalculateTMB(int, double)"/>.
    /// </summary>
    /// <param name="calls">Classified somatic calls; only <see cref="SomaticStatus.Somatic"/> are counted.</param>
    /// <param name="targetRegionMb">Size of the sequenced coding region in megabases (&gt; 0, finite).</param>
    /// <returns>TMB in mutations per megabase (≥ 0).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="calls"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="targetRegionMb"/> is not finite and &gt; 0.</exception>
    public static double CalculateTMB(IEnumerable<SomaticCall> calls, double targetRegionMb)
    {
        ArgumentNullException.ThrowIfNull(calls);

        int somaticCount = calls.Count(c => c.Status == SomaticStatus.Somatic);
        return CalculateTMB(somaticCount, targetRegionMb);
    }

    /// <summary>
    /// Classifies a TMB value (mut/Mb) as <see cref="TmbStatus.High"/> when it is at or above the FDA
    /// TMB-High cutoff (≥ <see cref="TmbHighThreshold"/> = 10 mut/Mb; the boundary is inclusive), otherwise
    /// <see cref="TmbStatus.Low"/>. Source: Marcus et al. (2021), FDA Approval Summary — pembrolizumab for
    /// "TMB ≥10 mut/Mb" solid tumors (companion diagnostic FoundationOne CDx).
    /// </summary>
    /// <param name="tmb">Tumor mutational burden in mutations per megabase (≥ 0, finite).</param>
    /// <returns><see cref="TmbStatus.High"/> if tmb ≥ 10, else <see cref="TmbStatus.Low"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="tmb"/> is negative or not finite.</exception>
    public static TmbStatus ClassifyTMB(double tmb)
    {
        if (double.IsNaN(tmb) || double.IsInfinity(tmb) || tmb < 0.0)
        {
            throw new ArgumentOutOfRangeException(nameof(tmb), tmb, "TMB must be a finite value ≥ 0.");
        }

        return tmb >= TmbHighThreshold ? TmbStatus.High : TmbStatus.Low;
    }

    #endregion


    #region Microsatellite Instability (MSI)

    /// <summary>
    /// MSI score cutoff (as a fraction in [0,1]) at or above which a sample is microsatellite-instability-high
    /// (MSI-H). Source: niu-lab/msisensor2 README — "The recommended msi score cutoff value is 20%
    /// (msi high: msi score &gt;= 20%)", where the msi score is "number of msi sites / all valid sites". The
    /// boundary is inclusive. Value = 0.20 (20%).
    /// </summary>
    public const double MsiHighScoreThreshold = 0.20;

    /// <summary>
    /// Minimum number of unstable markers (out of the validated 5-marker reference panel) for a tumor to be
    /// classified MSI-H under the NCI/Bethesda categorical criteria. Source: Boland CR et al. (1998), Cancer
    /// Res 58(22):5248–5257 — a tumor is MSI-H "if two or more of the five markers show instability". Value = 2.
    /// </summary>
    public const int BethesdaMsiHighMarkerCount = 2;

    /// <summary>
    /// Number of unstable markers for MSI-L under the NCI/Bethesda criteria — "only one of the five markers
    /// shows instability". Source: Boland et al. (1998). Value = 1.
    /// </summary>
    public const int BethesdaMsiLowMarkerCount = 1;

    /// <summary>
    /// Microsatellite-instability status of a sample.
    /// <list type="bullet">
    /// <item><description><b>MSS</b> — microsatellite stable.</description></item>
    /// <item><description><b>MSI_Low</b> — low-frequency MSI (Bethesda: exactly 1 of 5 markers unstable).</description></item>
    /// <item><description><b>MSI_High</b> — high-frequency MSI (Bethesda ≥2/5; MSIsensor2 score ≥20%).</description></item>
    /// </list>
    /// Sources: Boland et al. (1998) (categorical MSS/MSI-L/MSI-H); niu-lab/msisensor2 (continuous score ≥20% → MSI-H).
    /// </summary>
    public enum MsiStatus
    {
        /// <summary>Microsatellite stable (no instability above the calling threshold).</summary>
        MSS,

        /// <summary>Low-frequency microsatellite instability (Bethesda: exactly one marker unstable).</summary>
        MSI_Low,

        /// <summary>High-frequency microsatellite instability (Bethesda ≥2/5 markers; MSIsensor2 score ≥20%).</summary>
        MSI_High
    }

    /// <summary>Result of an end-to-end MSI determination over a set of evaluated microsatellite loci.</summary>
    /// <param name="UnstableLoci">Number of loci called unstable (somatic indel in tumor vs normal).</param>
    /// <param name="TotalLoci">Total number of valid evaluated microsatellite loci.</param>
    /// <param name="Score">MSI score = <paramref name="UnstableLoci"/> / <paramref name="TotalLoci"/>, in [0,1].</param>
    /// <param name="Status">MSI-H (score ≥ 20%) or MSS, per the MSIsensor2 continuous-score cutoff.</param>
    public readonly record struct MsiResult(int UnstableLoci, int TotalLoci, double Score, MsiStatus Status);

    /// <summary>
    /// Computes the MSI score as the fraction of unstable microsatellite loci among all valid evaluated loci:
    /// score = <paramref name="unstableLoci"/> / <paramref name="totalLoci"/>. Source: niu-lab/msisensor2 —
    /// "the msi score (number of msi sites / all valid sites)"; Niu et al. (2014), Bioinformatics 30(7):1015 —
    /// the MSI score is the percentage of microsatellite sites with a somatic indel. Returned as a fraction in
    /// [0,1] (multiply by 100 for a percentage).
    /// </summary>
    /// <param name="unstableLoci">Number of loci called unstable (0 ≤ unstableLoci ≤ totalLoci).</param>
    /// <param name="totalLoci">Total number of valid evaluated loci (&gt; 0).</param>
    /// <returns>MSI score in [0,1].</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="totalLoci"/> ≤ 0 (score undefined with no valid loci), <paramref name="unstableLoci"/>
    /// is negative, or <paramref name="unstableLoci"/> &gt; <paramref name="totalLoci"/>.
    /// </exception>
    public static double CalculateMSIScore(int unstableLoci, int totalLoci)
    {
        if (totalLoci <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(totalLoci), totalLoci,
                "Total valid loci must be > 0; the MSI score is undefined with no evaluable loci.");
        }

        if (unstableLoci < 0 || unstableLoci > totalLoci)
        {
            throw new ArgumentOutOfRangeException(
                nameof(unstableLoci), unstableLoci,
                "Unstable loci must satisfy 0 ≤ unstableLoci ≤ totalLoci.");
        }

        // MSI score = unstable loci / valid loci (MSIsensor2; Niu et al. 2014).
        return (double)unstableLoci / totalLoci;
    }

    /// <summary>
    /// Classifies a continuous MSI score (fraction in [0,1]) as <see cref="MsiStatus.MSI_High"/> when it is at
    /// or above the MSIsensor2 cutoff (≥ <see cref="MsiHighScoreThreshold"/> = 20%; boundary inclusive),
    /// otherwise <see cref="MsiStatus.MSS"/>. Source: niu-lab/msisensor2 README — "msi high: msi score &gt;= 20%".
    /// MSIsensor2 defines only a binary MSI-H cutoff on the continuous score, so no MSI-L band is applied here;
    /// MSI-L is a marker-count concept handled by <see cref="ClassifyBethesdaPanel"/>.
    /// </summary>
    /// <param name="score">MSI score as a fraction in [0,1].</param>
    /// <returns><see cref="MsiStatus.MSI_High"/> if score ≥ 0.20, else <see cref="MsiStatus.MSS"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="score"/> is not a finite value in [0,1].</exception>
    public static MsiStatus ClassifyMSIStatus(double score)
    {
        if (double.IsNaN(score) || double.IsInfinity(score) || score < 0.0 || score > 1.0)
        {
            throw new ArgumentOutOfRangeException(nameof(score), score, "MSI score must be a finite value in [0,1].");
        }

        return score >= MsiHighScoreThreshold ? MsiStatus.MSI_High : MsiStatus.MSS;
    }

    /// <summary>
    /// Classifies a sample under the NCI/Bethesda categorical criteria from the number of unstable markers in a
    /// fixed reference panel: <see cref="MsiStatus.MSI_High"/> when ≥ 2 markers are unstable,
    /// <see cref="MsiStatus.MSI_Low"/> when exactly 1 is unstable, and <see cref="MsiStatus.MSS"/> when none is.
    /// Source: Boland CR et al. (1998), Cancer Res 58(22):5248–5257 — MSI-H "if two or more of the five markers
    /// show instability", MSI-L "if only one of the five markers shows instability", MSS otherwise.
    /// </summary>
    /// <param name="unstableMarkers">Number of unstable markers (0 ≤ unstableMarkers ≤ totalMarkers).</param>
    /// <param name="totalMarkers">Total markers evaluated in the panel (&gt; 0; classically 5).</param>
    /// <returns>The Bethesda MSI status.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// <paramref name="totalMarkers"/> ≤ 0, <paramref name="unstableMarkers"/> negative, or
    /// <paramref name="unstableMarkers"/> &gt; <paramref name="totalMarkers"/>.
    /// </exception>
    public static MsiStatus ClassifyBethesdaPanel(int unstableMarkers, int totalMarkers)
    {
        if (totalMarkers <= 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(totalMarkers), totalMarkers, "Total markers must be > 0.");
        }

        if (unstableMarkers < 0 || unstableMarkers > totalMarkers)
        {
            throw new ArgumentOutOfRangeException(
                nameof(unstableMarkers), unstableMarkers,
                "Unstable markers must satisfy 0 ≤ unstableMarkers ≤ totalMarkers.");
        }

        // Boland et al. (1998): ≥2 unstable → MSI-H; exactly 1 → MSI-L; 0 → MSS.
        if (unstableMarkers >= BethesdaMsiHighMarkerCount)
        {
            return MsiStatus.MSI_High;
        }

        return unstableMarkers == BethesdaMsiLowMarkerCount ? MsiStatus.MSI_Low : MsiStatus.MSS;
    }

    /// <summary>
    /// End-to-end MSI determination from per-locus stability calls: counts the unstable loci, computes the MSI
    /// score (<see cref="CalculateMSIScore"/>), and classifies it with the MSIsensor2 continuous-score cutoff
    /// (<see cref="ClassifyMSIStatus"/>). Each element of <paramref name="locusUnstableFlags"/> is one valid
    /// evaluated microsatellite locus: <c>true</c> = unstable (somatic indel), <c>false</c> = stable. This
    /// matches the MSIsensor pipeline where per-locus instability (chi-square tumor-vs-normal length
    /// distributions) is determined upstream and the sample score is the fraction of unstable loci
    /// (Niu et al. 2014; niu-lab/msisensor2).
    /// </summary>
    /// <param name="locusUnstableFlags">Per-locus stability flags (one entry per valid evaluated locus).</param>
    /// <returns>The unstable/total counts, the MSI score, and the MSI status.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="locusUnstableFlags"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The sequence is empty (no valid loci to evaluate).</exception>
    public static MsiResult DetectMSI(IEnumerable<bool> locusUnstableFlags)
    {
        ArgumentNullException.ThrowIfNull(locusUnstableFlags);

        int total = 0;
        int unstable = 0;
        foreach (bool isUnstable in locusUnstableFlags)
        {
            total++;
            if (isUnstable)
            {
                unstable++;
            }
        }

        if (total == 0)
        {
            throw new ArgumentOutOfRangeException(
                nameof(locusUnstableFlags), "At least one valid locus is required; the MSI score is undefined for an empty set.");
        }

        double score = CalculateMSIScore(unstable, total);
        return new MsiResult(unstable, total, score, ClassifyMSIStatus(score));
    }

    #endregion

}
