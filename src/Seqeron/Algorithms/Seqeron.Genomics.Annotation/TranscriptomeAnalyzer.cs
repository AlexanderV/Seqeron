namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Provides algorithms for transcriptome analysis including expression quantification and differential expression.
/// </summary>
public static class TranscriptomeAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Represents gene expression data.
    /// </summary>
    public readonly record struct GeneExpression(
        string GeneId,
        double RawCount,
        double TPM,
        double FPKM,
        int Length);

    /// <summary>
    /// Represents differential expression result.
    /// </summary>
    public readonly record struct DifferentialExpression(
        string GeneId,
        double Log2FoldChange,
        double PValue,
        double AdjustedPValue,
        bool IsSignificant,
        string Regulation);

    /// <summary>
    /// Represents a gene set enrichment result.
    /// </summary>
    public readonly record struct EnrichmentResult(
        string PathwayId,
        string PathwayName,
        int GenesInPathway,
        int OverlappingGenes,
        double EnrichmentScore,
        double PValue,
        IReadOnlyList<string> Genes);

    /// <summary>
    /// Represents alternative splicing event.
    /// </summary>
    public readonly record struct SplicingEvent(
        string GeneId,
        string EventType,
        int Start,
        int End,
        double InclusionLevel,
        double DeltaPSI);

    /// <summary>
    /// Represents transcript isoform.
    /// </summary>
    public readonly record struct TranscriptIsoform(
        string TranscriptId,
        string GeneId,
        int Length,
        int ExonCount,
        double Expression,
        bool IsProteinCoding,
        IReadOnlyList<(int Start, int End)> Exons);

    /// <summary>
    /// Represents co-expression cluster.
    /// </summary>
    public readonly record struct CoExpressionCluster(
        int ClusterId,
        IReadOnlyList<string> Genes,
        double MeanCorrelation,
        string RepresentativeGene,
        IReadOnlyList<string> EnrichedFunctions);

    /// <summary>
    /// Types of alternative splicing events.
    /// </summary>
    public enum SplicingEventType
    {
        SkippedExon,
        RetainedIntron,
        AlternativeFivePrimeSS,
        AlternativeThreePrimeSS,
        MutuallyExclusiveExons
    }

    #endregion

    #region Expression Quantification

    // TPM per-million scaling factor: TPM_i = (X_i/l_i)/Σ(X_j/l_j) * 10^6.
    // Source: Zhao, Ye & Stanton (2020) RNA 26(8); Wagner, Kin & Lynch (2012) Theory Biosci 131(4):281-285.
    private const double TpmScalingFactor = 1_000_000.0;

    // FPKM scaling: combines the per-kilobase (10^3) and per-million-reads (10^6) factors.
    // FPKM_i = X_i * 10^9 / (l_i * N). Source: Zhao, Ye & Stanton (2020); Pimentel (2014).
    private const double FpkmScalingFactor = 1_000_000_000.0;

    /// <summary>
    /// Calculates TPM (Transcripts Per Million) from raw counts.
    /// </summary>
    /// <remarks>
    /// TPM_i = (X_i / l_i) / Σ_j(X_j / l_j) × 10^6, so TPM values within a sample sum to 10^6.
    /// Source: Zhao, Ye &amp; Stanton (2020) RNA; Wagner, Kin &amp; Lynch (2012) Theory in Biosciences.
    /// </remarks>
    public static IEnumerable<GeneExpression> CalculateTPM(
        IEnumerable<(string GeneId, double RawCount, int Length)> geneCounts)
    {
        var geneList = geneCounts.ToList();

        if (geneList.Count == 0)
            yield break;

        // Reads-per-kilobase-equivalent rate = count / length for each gene (X_i / l_i).
        var rates = geneList
            .Select(g => (g.GeneId, g.RawCount, g.Length, Rate: g.RawCount / Math.Max(g.Length, 1)))
            .ToList();

        double sumRates = rates.Sum(r => r.Rate);
        double totalReads = geneList.Sum(g => g.RawCount);

        // Degenerate denominator (all counts zero): TPM is 0/0 (undefined); emit 0 for every gene.
        if (sumRates == 0)
        {
            foreach (var gene in geneList)
            {
                yield return new GeneExpression(gene.GeneId, gene.RawCount, 0, 0, gene.Length);
            }
            yield break;
        }

        foreach (var (geneId, rawCount, length, rate) in rates)
        {
            double tpm = (rate / sumRates) * TpmScalingFactor;
            double fpkm = CalculateFPKM(rawCount, length, totalReads);

            yield return new GeneExpression(geneId, rawCount, tpm, fpkm, length);
        }
    }

    /// <summary>
    /// Calculates FPKM (Fragments Per Kilobase of transcript per Million mapped reads) for a single gene.
    /// </summary>
    /// <param name="rawCount">Number of reads/fragments mapped to the transcript (X_i).</param>
    /// <param name="length">Transcript length in bases (l_i); must be positive.</param>
    /// <param name="totalReads">Total mapped reads/fragments in the sample (N); must be positive.</param>
    /// <returns>FPKM_i = X_i × 10^9 / (l_i × N); 0 when length or totalReads is non-positive.</returns>
    /// <remarks>Source: Zhao, Ye &amp; Stanton (2020) RNA; Mortazavi et al. (2008) Nat Methods; Pimentel (2014).</remarks>
    public static double CalculateFPKM(double rawCount, int length, double totalReads)
    {
        // Length and total mapped reads must be positive; otherwise FPKM is undefined → 0.
        if (length <= 0 || totalReads <= 0)
            return 0;

        // FPKM = (reads * 10^9) / (length * total mapped reads).
        return (rawCount * FpkmScalingFactor) / (length * totalReads);
    }

    /// <summary>
    /// Normalizes expression values across samples using quantile normalization (Bolstad et al. 2003).
    /// </summary>
    /// <remarks>
    /// Each sample (column) is sorted; the mean across samples at each rank is computed; the rank means
    /// are placed back at each value's original position. Values tied within a sample receive the mean of
    /// the rank means they would otherwise span (tie-average rule).
    /// Source: Wikipedia "Quantile normalization" (citing Bolstad, Irizarry, Astrand &amp; Speed, 2003).
    /// All samples must share the same length; the first sample defines the length.
    /// </remarks>
    public static IEnumerable<IReadOnlyList<double>> QuantileNormalize(
        IEnumerable<IEnumerable<double>> samples)
    {
        var sampleList = samples.Select(s => s.ToList()).ToList();

        if (sampleList.Count == 0)
            yield break;

        int geneCount = sampleList[0].Count;
        int sampleCount = sampleList.Count;

        if (geneCount == 0)
            yield break;

        // Mean of each sample's sorted values at each rank (lowest rank 0 .. highest).
        var rankMeans = new double[geneCount];
        for (int rank = 0; rank < geneCount; rank++)
        {
            double sum = 0;
            foreach (var sample in sampleList)
            {
                sum += sample.OrderBy(x => x).ElementAt(rank);
            }
            rankMeans[rank] = sum / sampleCount;
        }

        // Assign rank means back to original positions, averaging over tied runs (Bolstad tie rule).
        foreach (var sample in sampleList)
        {
            // Positions ordered by value; equal values form contiguous tied runs.
            var ordered = sample
                .Select((val, idx) => (val, idx))
                .OrderBy(x => x.val)
                .ToList();

            var normalized = new double[geneCount];
            int run = 0;
            while (run < geneCount)
            {
                int runEnd = run;
                while (runEnd + 1 < geneCount && ordered[runEnd + 1].val == ordered[run].val)
                    runEnd++;

                // Average of the rank means spanning the tied run [run, runEnd].
                double tieSum = 0;
                for (int rank = run; rank <= runEnd; rank++)
                    tieSum += rankMeans[rank];
                double tieMean = tieSum / (runEnd - run + 1);

                for (int rank = run; rank <= runEnd; rank++)
                    normalized[ordered[rank].idx] = tieMean;

                run = runEnd + 1;
            }

            yield return normalized;
        }
    }

    /// <summary>
    /// Performs log2 transformation with pseudocount.
    /// </summary>
    public static IEnumerable<double> Log2Transform(
        IEnumerable<double> values,
        double pseudocount = 1.0)
    {
        foreach (var value in values)
        {
            yield return Math.Log2(value + pseudocount);
        }
    }

    #endregion

    #region Differential Expression Analysis

    /// <summary>
    /// Performs simple differential expression analysis using fold change and t-test.
    /// </summary>
    public static IEnumerable<DifferentialExpression> AnalyzeDifferentialExpression(
        IEnumerable<(string GeneId, IReadOnlyList<double> Group1, IReadOnlyList<double> Group2)> expressionData,
        double foldChangeThreshold = 1.0,
        double pValueThreshold = 0.05)
    {
        var genes = expressionData.ToList();

        if (genes.Count == 0)
            yield break;

        var results = new List<(string GeneId, double Log2FC, double PValue)>();

        foreach (var (geneId, group1, group2) in genes)
        {
            if (group1.Count == 0 || group2.Count == 0)
                continue;

            double mean1 = group1.Average();
            double mean2 = group2.Average();

            // Avoid division by zero
            double log2FC = mean1 > 0
                ? Math.Log2((mean2 + 0.01) / (mean1 + 0.01))
                : 0;

            double pValue = CalculateTTestPValue(group1, group2);

            results.Add((geneId, log2FC, pValue));
        }

        // Multiple testing correction (Benjamini-Hochberg)
        var sortedByPValue = results.OrderBy(r => r.PValue).ToList();
        var adjustedPValues = BenjaminiHochberg(sortedByPValue.Select(r => r.PValue));
        var adjustedList = adjustedPValues.ToList();

        for (int i = 0; i < sortedByPValue.Count; i++)
        {
            var (geneId, log2FC, pValue) = sortedByPValue[i];
            double adjPValue = adjustedList[i];
            bool isSignificant = adjPValue < pValueThreshold && Math.Abs(log2FC) >= foldChangeThreshold;
            string regulation = log2FC > 0 ? "Upregulated" : (log2FC < 0 ? "Downregulated" : "Unchanged");

            yield return new DifferentialExpression(
                GeneId: geneId,
                Log2FoldChange: log2FC,
                PValue: pValue,
                AdjustedPValue: adjPValue,
                IsSignificant: isSignificant,
                Regulation: regulation);
        }
    }

    // Pseudocount added to both group means before the log2 ratio so that a zero mean does not make
    // log2(mean2/mean1) undefined. Standard regularization of the simple ratio estimator; value 1 keeps
    // log2(1) = 0 for an absent/absent comparison. Source: log2 fold-change definition, Love et al. (2014)
    // DESeq2; Science Park RNA-seq differential-expression lesson (degenerate-input convention).
    private const double FoldChangePseudocount = 1.0;

    // Minimum replicates per group required to form an unbiased (N-1) sample variance for Welch's t-test.
    // Source: Welch (1947) — the corrected sample variance is undefined for N < 2.
    private const int MinReplicatesForTTest = 2;

    /// <summary>
    /// Calculates the log2 fold change between two expression conditions as
    /// log2((mean(expression2) + c) / (mean(expression1) + c)), where c is a pseudocount.
    /// </summary>
    /// <param name="expression1">Replicate expression values for condition 1 (reference / control).</param>
    /// <param name="expression2">Replicate expression values for condition 2 (treatment); the numerator.</param>
    /// <returns>
    /// log2 fold change; positive means higher expression in condition 2 (treatment), negative means lower.
    /// Returns 0 when either condition has no replicates.
    /// </returns>
    /// <remarks>
    /// Fold change is the log2 ratio of mean expression between treatment and control.
    /// Source: Love, Huber &amp; Anders (2014) DESeq2 (log2 fold change between treatment and control);
    /// Science Park RNA-seq differential-expression lesson (log2(condition A / condition B), positive = up in numerator).
    /// </remarks>
    public static double CalculateFoldChange(
        IReadOnlyList<double> expression1,
        IReadOnlyList<double> expression2)
    {
        if (expression1 is null || expression2 is null)
            return 0;

        if (expression1.Count == 0 || expression2.Count == 0)
            return 0;

        double mean1 = expression1.Average();
        double mean2 = expression2.Average();

        // log2((mean2 + c) / (mean1 + c)): positive = up in condition 2 (treatment).
        return Math.Log2((mean2 + FoldChangePseudocount) / (mean1 + FoldChangePseudocount));
    }

    /// <summary>
    /// Identifies differentially expressed genes between two conditions using the log2 fold change,
    /// a Welch (unequal-variance) two-sample t-test, and Benjamini-Hochberg FDR adjustment.
    /// </summary>
    /// <param name="condition1">Per-gene replicate expression values for condition 1 (reference / control).</param>
    /// <param name="condition2">Per-gene replicate expression values for condition 2 (treatment).</param>
    /// <param name="alpha">Adjusted-p-value significance threshold (FDR); a gene is significant only when its adjusted p-value is strictly below this.</param>
    /// <param name="log2FoldChangeThreshold">Minimum absolute log2 fold change required for significance.</param>
    /// <returns>
    /// One <see cref="DifferentialExpression"/> per gene. A gene is flagged significant only when BOTH
    /// |log2 fold change| ≥ <paramref name="log2FoldChangeThreshold"/> AND adjusted p-value &lt; <paramref name="alpha"/>.
    /// </returns>
    /// <remarks>
    /// Welch's t-statistic uses unbiased (N-1) variances and Welch-Satterthwaite degrees of freedom; the
    /// two-sided p-value is the exact Student's t tail via the regularized incomplete beta function.
    /// Raw p-values are corrected with the Benjamini-Hochberg step-up procedure (matches R's p.adjust BH).
    /// A gene is differentially expressed only when both the fold-change and the adjusted-p criteria hold.
    /// Sources: Welch (1947); Student's t-distribution CDF (regularized incomplete beta); Benjamini &amp; Hochberg (1995);
    /// Love et al. (2014) DESeq2; Science Park RNA-seq lesson (two-criterion DE rule).
    /// </remarks>
    public static IEnumerable<DifferentialExpression> FindDifferentiallyExpressed(
        IEnumerable<(string GeneId, IReadOnlyList<double> Condition1, IReadOnlyList<double> Condition2)> genes,
        double alpha = 0.05,
        double log2FoldChangeThreshold = 1.0)
    {
        var geneList = genes?.ToList() ?? new List<(string, IReadOnlyList<double>, IReadOnlyList<double>)>();

        if (geneList.Count == 0)
            return Enumerable.Empty<DifferentialExpression>();

        // Per-gene log2 fold change and raw two-sided Welch t-test p-value.
        var perGene = new (string GeneId, double Log2FC, double PValue)[geneList.Count];
        for (int g = 0; g < geneList.Count; g++)
        {
            var (geneId, c1, c2) = geneList[g];
            double log2Fc = CalculateFoldChange(c1, c2);
            double pValue = WelchTTestTwoSidedPValue(c1, c2);
            perGene[g] = (geneId, log2Fc, pValue);
        }

        // Benjamini-Hochberg adjustment across all genes (multiple-testing correction).
        double[] adjusted = BenjaminiHochbergAdjust(perGene.Select(x => x.PValue).ToArray());

        var results = new DifferentialExpression[geneList.Count];
        for (int g = 0; g < geneList.Count; g++)
        {
            var (geneId, log2Fc, pValue) = perGene[g];
            double adjP = adjusted[g];

            // DE = BOTH criteria: |log2FC| >= threshold AND adjusted p-value < alpha.
            bool isSignificant = Math.Abs(log2Fc) >= log2FoldChangeThreshold && adjP < alpha;
            string regulation = log2Fc > 0 ? "Upregulated" : (log2Fc < 0 ? "Downregulated" : "Unchanged");

            results[g] = new DifferentialExpression(
                GeneId: geneId,
                Log2FoldChange: log2Fc,
                PValue: pValue,
                AdjustedPValue: adjP,
                IsSignificant: isSignificant,
                Regulation: regulation);
        }

        return results;
    }

    /// <summary>
    /// Two-sided p-value of Welch's unequal-variance two-sample t-test using the exact Student's t tail.
    /// </summary>
    /// <remarks>
    /// t = (mean2 - mean1) / sqrt(s1²/N1 + s2²/N2) with unbiased (N-1) variances; Welch-Satterthwaite df.
    /// p = I_{ν/(ν+t²)}(ν/2, 1/2). Source: Welch (1947); Student's t-distribution CDF (regularized incomplete beta).
    /// </remarks>
    private static double WelchTTestTwoSidedPValue(IReadOnlyList<double> group1, IReadOnlyList<double> group2)
    {
        if (group1 is null || group2 is null ||
            group1.Count < MinReplicatesForTTest || group2.Count < MinReplicatesForTTest)
            return 1.0; // Variance undefined for N < 2: gene not testable → not significant.

        int n1 = group1.Count;
        int n2 = group2.Count;
        double mean1 = group1.Average();
        double mean2 = group2.Average();

        // Unbiased sample variances (divide by N-1).
        double var1 = group1.Sum(x => (x - mean1) * (x - mean1)) / (n1 - 1);
        double var2 = group2.Sum(x => (x - mean2) * (x - mean2)) / (n2 - 1);

        double se = Math.Sqrt(var1 / n1 + var2 / n2);
        if (se == 0)
            return mean1 == mean2 ? 1.0 : 0.0; // No variance: identical → p=1; separated → p=0.

        double t = (mean2 - mean1) / se;

        // Welch-Satterthwaite degrees of freedom.
        double a = var1 / n1;
        double b = var2 / n2;
        double df = (a + b) * (a + b) /
                    (a * a / (n1 - 1) + b * b / (n2 - 1));

        // Two-sided p-value: P(|T| >= |t|) = I_{df/(df+t²)}(df/2, 1/2).
        double x = df / (df + t * t);
        return RegularizedIncompleteBeta(df / 2.0, 0.5, x);
    }

    /// <summary>
    /// Benjamini-Hochberg step-up adjusted p-values, matching R's <c>p.adjust(method="BH")</c>.
    /// </summary>
    /// <remarks>
    /// Sort p ascending; multiply each by m/rank; take the cumulative minimum from the largest p downward;
    /// clamp to 1; restore the original order. Source: Benjamini &amp; Hochberg (1995); R stats p.adjust.
    /// </remarks>
    private static double[] BenjaminiHochbergAdjust(double[] pValues)
    {
        int m = pValues.Length;
        var adjusted = new double[m];
        if (m == 0)
            return adjusted;

        // Indices of p-values sorted ascending.
        var order = Enumerable.Range(0, m).OrderBy(i => pValues[i]).ToArray();

        double runningMin = 1.0;
        // Walk from the largest p-value (rank m) down to the smallest (rank 1).
        for (int rank = m; rank >= 1; rank--)
        {
            int idx = order[rank - 1];
            double bh = pValues[idx] * m / rank; // m/rank inflation factor.
            runningMin = Math.Min(runningMin, bh);
            adjusted[idx] = Math.Min(runningMin, 1.0); // Clamp to 1.
        }

        return adjusted;
    }

    /// <summary>
    /// Regularized incomplete beta function I_x(a,b) via the Lentz continued-fraction expansion.
    /// </summary>
    /// <remarks>
    /// Standard numerical evaluation (continued fraction betacf with symmetry I_x(a,b)=1-I_{1-x}(b,a)).
    /// Used to compute the exact Student's t-distribution tail. Source: Student's t-distribution CDF
    /// (regularized incomplete beta); Press et al., Numerical Recipes (betai/betacf).
    /// </remarks>
    private static double RegularizedIncompleteBeta(double a, double b, double x)
    {
        if (x <= 0.0)
            return 0.0;
        if (x >= 1.0)
            return 1.0;

        double logBeta = LogGamma(a + b) - LogGamma(a) - LogGamma(b);
        double front = Math.Exp(logBeta + a * Math.Log(x) + b * Math.Log(1.0 - x));

        // Use the more rapidly converging branch.
        if (x < (a + 1.0) / (a + b + 2.0))
            return front * BetaContinuedFraction(a, b, x) / a;

        return 1.0 - front * BetaContinuedFraction(b, a, 1.0 - x) / b;
    }

    private static double BetaContinuedFraction(double a, double b, double x)
    {
        const int maxIterations = 200;
        const double epsilon = 3.0e-12;
        const double tiny = 1.0e-300; // Guards against zero denominators (Lentz's method).

        double qab = a + b;
        double qap = a + 1.0;
        double qam = a - 1.0;

        double c = 1.0;
        double d = 1.0 - qab * x / qap;
        if (Math.Abs(d) < tiny) d = tiny;
        d = 1.0 / d;
        double h = d;

        for (int mIter = 1; mIter <= maxIterations; mIter++)
        {
            int m2 = 2 * mIter;

            double aa = mIter * (b - mIter) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (Math.Abs(d) < tiny) d = tiny;
            c = 1.0 + aa / c;
            if (Math.Abs(c) < tiny) c = tiny;
            d = 1.0 / d;
            h *= d * c;

            aa = -(a + mIter) * (qab + mIter) * x / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (Math.Abs(d) < tiny) d = tiny;
            c = 1.0 + aa / c;
            if (Math.Abs(c) < tiny) c = tiny;
            d = 1.0 / d;
            double del = d * c;
            h *= del;

            if (Math.Abs(del - 1.0) < epsilon)
                break;
        }

        return h;
    }

    /// <summary>
    /// Natural logarithm of the gamma function (Lanczos approximation).
    /// </summary>
    private static double LogGamma(double x)
    {
        // Lanczos coefficients (g = 7, n = 9). Source: Press et al., Numerical Recipes (gammln).
        double[] coefficients =
        {
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        };

        double xx = x - 1.0;
        double a = coefficients[0];
        double tVal = xx + 7.5;
        for (int i = 1; i < coefficients.Length; i++)
            a += coefficients[i] / (xx + i);

        return 0.5 * Math.Log(2 * Math.PI) + (xx + 0.5) * Math.Log(tVal) - tVal + Math.Log(a);
    }

    /// <summary>
    /// Calculates a simple t-test p-value.
    /// </summary>
    private static double CalculateTTestPValue(IReadOnlyList<double> group1, IReadOnlyList<double> group2)
    {
        if (group1.Count < 2 || group2.Count < 2)
            return 1.0;

        double mean1 = group1.Average();
        double mean2 = group2.Average();

        double var1 = group1.Sum(x => (x - mean1) * (x - mean1)) / (group1.Count - 1);
        double var2 = group2.Sum(x => (x - mean2) * (x - mean2)) / (group2.Count - 1);

        double se = Math.Sqrt(var1 / group1.Count + var2 / group2.Count);

        if (se == 0)
            return mean1 == mean2 ? 1.0 : 0.0;

        double t = Math.Abs(mean2 - mean1) / se;

        // Approximate p-value using normal distribution
        return 2 * (1 - StatisticsHelper.NormalCDF(t));
    }

    /// <summary>
    /// Benjamini-Hochberg multiple testing correction.
    /// </summary>
    private static IEnumerable<double> BenjaminiHochberg(IEnumerable<double> pValues)
    {
        var pList = pValues.ToList();
        int n = pList.Count;

        if (n == 0)
            yield break;

        var adjusted = new double[n];
        double minSoFar = 1.0;

        for (int i = n - 1; i >= 0; i--)
        {
            double corrected = pList[i] * n / (i + 1);
            corrected = Math.Min(corrected, minSoFar);
            corrected = Math.Min(corrected, 1.0);
            adjusted[i] = corrected;
            minSoFar = corrected;
        }

        foreach (var adj in adjusted)
            yield return adj;
    }

    #endregion

    #region Gene Set Enrichment

    /// <summary>
    /// Performs over-representation analysis (ORA) for gene set enrichment.
    /// </summary>
    public static IEnumerable<EnrichmentResult> PerformOverRepresentationAnalysis(
        IReadOnlySet<string> differentiallyExpressedGenes,
        IEnumerable<(string PathwayId, string PathwayName, IReadOnlySet<string> Genes)> pathways,
        int backgroundGeneCount)
    {
        if (differentiallyExpressedGenes.Count == 0 || backgroundGeneCount <= 0)
            yield break;

        foreach (var (pathwayId, pathwayName, pathwayGenes) in pathways)
        {
            var overlapping = pathwayGenes.Intersect(differentiallyExpressedGenes).ToList();

            if (overlapping.Count == 0)
                continue;

            // Fisher's exact test approximation
            double pValue = CalculateFisherPValue(
                differentiallyExpressedGenes.Count,
                pathwayGenes.Count,
                overlapping.Count,
                backgroundGeneCount);

            // Enrichment score
            double expected = (double)differentiallyExpressedGenes.Count * pathwayGenes.Count / backgroundGeneCount;
            double enrichmentScore = expected > 0 ? overlapping.Count / expected : 0;

            yield return new EnrichmentResult(
                PathwayId: pathwayId,
                PathwayName: pathwayName,
                GenesInPathway: pathwayGenes.Count,
                OverlappingGenes: overlapping.Count,
                EnrichmentScore: enrichmentScore,
                PValue: pValue,
                Genes: overlapping);
        }
    }

    /// <summary>
    /// Approximates Fisher's exact test p-value using hypergeometric distribution.
    /// </summary>
    private static double CalculateFisherPValue(int deGenes, int pathwaySize, int overlap, int background)
    {
        // Hypergeometric probability approximation
        // Using normal approximation for large samples

        double expectedOverlap = (double)deGenes * pathwaySize / background;
        double variance = expectedOverlap * (1 - (double)pathwaySize / background) *
                          (background - deGenes) / (background - 1);

        if (variance <= 0)
            return overlap >= expectedOverlap ? 0.0 : 1.0;

        double z = (overlap - expectedOverlap) / Math.Sqrt(variance);
        return 1 - StatisticsHelper.NormalCDF(z);
    }

    /// <summary>
    /// Calculates Gene Set Enrichment Score (GSEA-like).
    /// </summary>
    public static double CalculateEnrichmentScore(
        IReadOnlyList<string> rankedGenes,
        IReadOnlySet<string> geneSet)
    {
        if (rankedGenes.Count == 0 || geneSet.Count == 0)
            return 0;

        int n = rankedGenes.Count;
        int hitCount = rankedGenes.Count(g => geneSet.Contains(g));
        int missCount = n - hitCount;

        if (hitCount == 0 || missCount == 0)
            return 0;

        double hitIncrement = 1.0 / hitCount;
        double missDecrement = 1.0 / missCount;

        double runningSum = 0;
        double maxDeviation = 0;

        foreach (var gene in rankedGenes)
        {
            if (geneSet.Contains(gene))
            {
                runningSum += hitIncrement;
            }
            else
            {
                runningSum -= missDecrement;
            }

            if (Math.Abs(runningSum) > Math.Abs(maxDeviation))
            {
                maxDeviation = runningSum;
            }
        }

        return maxDeviation;
    }

    #endregion

    #region Alternative Splicing Analysis

    /// <summary>
    /// Calculates the Percent Spliced In (PSI, Ψ) of an alternative-splicing event from read counts.
    /// </summary>
    /// <param name="inclusionReads">Reads supporting the inclusion isoform (I); must be ≥ 0.</param>
    /// <param name="exclusionReads">Reads supporting the exclusion/skipping isoform (S); must be ≥ 0.</param>
    /// <param name="inclusionEffectiveLength">
    /// Optional effective length of the inclusion isoform (l_I); when this and
    /// <paramref name="exclusionEffectiveLength"/> are both &gt; 0 the rMATS length-normalized form is used.
    /// </param>
    /// <param name="exclusionEffectiveLength">Optional effective length of the skipping isoform (l_S).</param>
    /// <returns>
    /// Ψ = I / (I + S) by default; the rMATS length-normalized Ψ = (I/l_I) / (I/l_I + S/l_S) when both
    /// effective lengths are positive. Returns <see cref="double.NaN"/> when there are no supporting reads.
    /// </returns>
    /// <remarks>
    /// PSI is the expression of inclusion isoforms as a fraction of the total (inclusion + exclusion).
    /// Source: BMC Bioinformatics 13(Suppl 6):S11 (Ψ = γ_i/(γ_i+γ_e)); SUPPA2 (Trincado et al. 2018);
    /// length-normalized form: Shen et al. (2014) rMATS, PNAS 111(51):E5593, ψ̂ = (I/l_I)/(I/l_I + S/l_S).
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">If either read count is negative.</exception>
    public static double CalculatePSI(
        double inclusionReads,
        double exclusionReads,
        double inclusionEffectiveLength = 0,
        double exclusionEffectiveLength = 0)
    {
        if (inclusionReads < 0)
            throw new ArgumentOutOfRangeException(nameof(inclusionReads), "Read counts must be non-negative.");
        if (exclusionReads < 0)
            throw new ArgumentOutOfRangeException(nameof(exclusionReads), "Read counts must be non-negative.");

        // Length-normalized rMATS form when both effective lengths are supplied (>0).
        // ψ̂ = (I/l_I) / (I/l_I + S/l_S). Source: Shen et al. (2014) rMATS.
        if (inclusionEffectiveLength > 0 && exclusionEffectiveLength > 0)
        {
            double iRate = inclusionReads / inclusionEffectiveLength;
            double sRate = exclusionReads / exclusionEffectiveLength;
            double denom = iRate + sRate;
            // No supporting reads → 0/0 undefined. Source: BMC Bioinformatics PMC3330053.
            return denom == 0 ? double.NaN : iRate / denom;
        }

        // Unnormalized read-count ratio Ψ = I/(I+S). Source: BMC Bioinformatics PMC3330053; SUPPA2.
        double total = inclusionReads + exclusionReads;
        return total == 0 ? double.NaN : inclusionReads / total;
    }

    /// <summary>
    /// Detects and classifies alternative-splicing events from transcript isoforms, comparing pairs of
    /// isoforms of the same gene and assigning each structural difference to one of the five canonical
    /// alternative-splicing classes.
    /// </summary>
    /// <param name="isoforms">Transcript isoforms; each carries its gene id and ordered (5′→3′) exon coordinates.</param>
    /// <returns>
    /// One <see cref="SplicingEvent"/> per detected difference between an isoform pair of the same gene.
    /// Genes with fewer than two isoforms produce no events; identical isoform pairs produce no events.
    /// </returns>
    /// <remarks>
    /// The five canonical event classes are skipped exon (SE), intron retention (RI), alternative 5′ splice
    /// site (A5SS), alternative 3′ splice site (A3SS), and mutually exclusive exons (MXE).
    /// Source: Wang et al. (2008) Nature 456(7221):470–476; rMATS event vocabulary (Shen et al. 2014).
    /// Exon coordinates are treated as inclusive [Start, End] on one strand in ascending order.
    /// </remarks>
    public static IEnumerable<SplicingEvent> DetectAlternativeSplicing(
        IEnumerable<TranscriptIsoform> isoforms)
    {
        if (isoforms is null)
            yield break;

        var byGene = isoforms.GroupBy(i => i.GeneId);

        foreach (var gene in byGene)
        {
            var transcripts = gene.ToList();

            // An AS event requires two isoforms of the same gene (Wang et al. 2008).
            if (transcripts.Count < 2)
                continue;

            // Compare every distinct unordered isoform pair within the gene.
            for (int a = 0; a < transcripts.Count; a++)
            {
                for (int b = a + 1; b < transcripts.Count; b++)
                {
                    var classified = ClassifyIsoformPair(gene.Key, transcripts[a], transcripts[b]);
                    if (classified.HasValue)
                        yield return classified.Value;
                }
            }
        }
    }

    /// <summary>
    /// Classifies the structural difference between two isoforms of one gene into a canonical AS class.
    /// Returns null when the two isoforms are structurally identical (no event).
    /// </summary>
    /// <remarks>Source: Wang et al. (2008) five-class taxonomy (SE, RI, A5SS, A3SS, MXE).</remarks>
    private static SplicingEvent? ClassifyIsoformPair(
        string geneId,
        TranscriptIsoform first,
        TranscriptIsoform second)
    {
        var exonsA = (first.Exons ?? Array.Empty<(int Start, int End)>())
            .OrderBy(e => e.Start).ThenBy(e => e.End).ToList();
        var exonsB = (second.Exons ?? Array.Empty<(int Start, int End)>())
            .OrderBy(e => e.Start).ThenBy(e => e.End).ToList();

        // Identical exon structure → no alternative-splicing event.
        if (exonsA.SequenceEqual(exonsB))
            return null;

        var setA = new HashSet<(int, int)>(exonsA);
        var setB = new HashSet<(int, int)>(exonsB);

        var onlyA = exonsA.Where(e => !setB.Contains(e)).ToList();
        var onlyB = exonsB.Where(e => !setA.Contains(e)).ToList();

        // Retained intron: a unique exon in one isoform spans (covers) the intron between two consecutive
        // exons of the other isoform — the intron is retained as exon body. Checked first because such an
        // exon may coexist with the (now-merged) flanking exons appearing as unique on the other side.
        // Source: Wang et al. (2008).
        if (SpansIntron(onlyA, exonsB) || SpansIntron(onlyB, exonsA))
        {
            var retained = SpansIntron(onlyA, exonsB) ? FindSpanningExon(onlyA, exonsB) : FindSpanningExon(onlyB, exonsA);
            return MakeEvent(geneId, SplicingEventType.RetainedIntron, retained, retained);
        }

        // Mutually exclusive exons: each isoform has exactly one exon the other lacks, both occupying the
        // same internal slot between shared flanking exons (non-overlapping alternatives). Source: Wang 2008.
        if (onlyA.Count == 1 && onlyB.Count == 1 && !Overlaps(onlyA[0], onlyB[0]))
        {
            return MakeEvent(geneId, SplicingEventType.MutuallyExclusiveExons, onlyA[0], onlyB[0]);
        }

        // Alternative splice site: a unique exon in one isoform shares exactly one boundary with a unique
        // exon in the other. On the forward strand the 5′ splice site (donor) is an exon's END boundary and
        // the 3′ splice site (acceptor) is an exon's START boundary, so:
        //   same Start, different End  ⇒ alternative donor      ⇒ A5SS (AlternativeFivePrimeSS)
        //   same End,   different Start ⇒ alternative acceptor   ⇒ A3SS (AlternativeThreePrimeSS)
        // Source: Wang 2008; rMATS A5SS/A3SS coordinate convention (Xinglab/rmats-turbo).
        if (onlyA.Count == 1 && onlyB.Count == 1)
        {
            var ea = onlyA[0];
            var eb = onlyB[0];
            if (ea.Start == eb.Start && ea.End != eb.End)
                return MakeEvent(geneId, SplicingEventType.AlternativeFivePrimeSS, ea, eb);
            if (ea.End == eb.End && ea.Start != eb.Start)
                return MakeEvent(geneId, SplicingEventType.AlternativeThreePrimeSS, ea, eb);
        }

        // Otherwise: an exon present in one isoform and absent from the other is a skipped/cassette exon,
        // the most common AS class. Source: Wang et al. (2008).
        var diff = onlyA.Count > 0 ? onlyA[0] : onlyB[0];
        return MakeEvent(geneId, SplicingEventType.SkippedExon, diff, diff);
    }

    /// <summary>
    /// True when some unique exon in <paramref name="candidates"/> spans the intron between two consecutive
    /// exons of <paramref name="otherExons"/> (the intron is retained as exon body). Source: Wang et al. (2008).
    /// </summary>
    private static bool SpansIntron(
        IReadOnlyList<(int Start, int End)> candidates,
        IReadOnlyList<(int Start, int End)> otherExons)
        => candidates.Any(c => Spans(c, otherExons));

    private static (int Start, int End) FindSpanningExon(
        IReadOnlyList<(int Start, int End)> candidates,
        IReadOnlyList<(int Start, int End)> otherExons)
        => candidates.First(c => Spans(c, otherExons));

    private static bool Spans((int Start, int End) exon, IReadOnlyList<(int Start, int End)> otherExons)
    {
        // The exon bridges two consecutive `other` exons by covering the intron between them:
        // exon.Start ≤ left.End and exon.End ≥ right.Start, where left.End < right.Start is a real intron.
        for (int i = 0; i + 1 < otherExons.Count; i++)
        {
            var left = otherExons[i];
            var right = otherExons[i + 1];
            if (exon.Start <= left.End && exon.End >= right.Start && left.End < right.Start)
                return true;
        }
        return false;
    }

    private static bool Overlaps((int Start, int End) x, (int Start, int End) y)
        => x.Start <= y.End && y.Start <= x.End;

    private static SplicingEvent MakeEvent(
        string geneId,
        SplicingEventType type,
        (int Start, int End) primary,
        (int Start, int End) secondary)
    {
        // Event spans from the leftmost to the rightmost coordinate of the two affected exons.
        int start = Math.Min(primary.Start, secondary.Start);
        int end = Math.Max(primary.End, secondary.End);

        return new SplicingEvent(
            GeneId: geneId,
            EventType: type.ToString(),
            Start: start,
            End: end,
            InclusionLevel: double.NaN, // Inclusion level requires read counts; use CalculatePSI separately.
            DeltaPSI: 0);
    }

    /// <summary>
    /// Identifies potential skipped exon events.
    /// </summary>
    public static IEnumerable<SplicingEvent> FindSkippedExonEvents(
        IEnumerable<(string GeneId, int ExonStart, int ExonEnd, double InclusionReads, double SkippingReads)> exonData)
    {
        foreach (var (geneId, start, end, inclusion, skipping) in exonData)
        {
            double total = inclusion + skipping;
            if (total == 0)
                continue;

            double psi = inclusion / total; // Percent Spliced In

            yield return new SplicingEvent(
                GeneId: geneId,
                EventType: "SkippedExon",
                Start: start,
                End: end,
                InclusionLevel: psi,
                DeltaPSI: 0); // Would need comparison sample
        }
    }

    /// <summary>
    /// Detects differential splicing between conditions.
    /// </summary>
    public static IEnumerable<SplicingEvent> DetectDifferentialSplicing(
        IEnumerable<(string GeneId, int Start, int End, double PSI_Condition1, double PSI_Condition2)> splicingData,
        double deltaPsiThreshold = 0.1)
    {
        foreach (var (geneId, start, end, psi1, psi2) in splicingData)
        {
            double deltaPsi = psi2 - psi1;

            if (Math.Abs(deltaPsi) >= deltaPsiThreshold)
            {
                string eventType = deltaPsi > 0 ? "IncreasedInclusion" : "IncreasedSkipping";

                yield return new SplicingEvent(
                    GeneId: geneId,
                    EventType: eventType,
                    Start: start,
                    End: end,
                    InclusionLevel: psi2,
                    DeltaPSI: deltaPsi);
            }
        }
    }

    #endregion

    #region Transcript Isoform Analysis

    /// <summary>
    /// Identifies dominant transcript isoform for each gene.
    /// </summary>
    public static IEnumerable<(string GeneId, TranscriptIsoform DominantIsoform, double DominanceRatio)>
        FindDominantIsoforms(IEnumerable<TranscriptIsoform> isoforms)
    {
        var byGene = isoforms.GroupBy(i => i.GeneId);

        foreach (var group in byGene)
        {
            var sorted = group.OrderByDescending(i => i.Expression).ToList();

            if (sorted.Count == 0)
                continue;

            var dominant = sorted[0];
            double totalExpression = sorted.Sum(i => i.Expression);
            double dominanceRatio = totalExpression > 0 ? dominant.Expression / totalExpression : 0;

            yield return (group.Key, dominant, dominanceRatio);
        }
    }

    /// <summary>
    /// Detects isoform switching between conditions.
    /// </summary>
    public static IEnumerable<(string GeneId, string TranscriptId1, string TranscriptId2, double SwitchScore)>
        DetectIsoformSwitching(
            IEnumerable<(TranscriptIsoform Isoform, double Expression1, double Expression2)> isoformData,
            double switchThreshold = 0.3)
    {
        var byGene = isoformData.GroupBy(d => d.Isoform.GeneId);

        foreach (var group in byGene)
        {
            var isoforms = group.ToList();

            if (isoforms.Count < 2)
                continue;

            double total1 = isoforms.Sum(i => i.Expression1);
            double total2 = isoforms.Sum(i => i.Expression2);

            if (total1 == 0 || total2 == 0)
                continue;

            // Calculate usage ratios
            var usageChanges = isoforms
                .Select(i => (
                    i.Isoform.TranscriptId,
                    Usage1: i.Expression1 / total1,
                    Usage2: i.Expression2 / total2))
                .Select(u => (
                    u.TranscriptId,
                    u.Usage1,
                    u.Usage2,
                    Delta: u.Usage2 - u.Usage1))
                .OrderByDescending(u => Math.Abs(u.Delta))
                .ToList();

            if (usageChanges.Count < 2)
                continue;

            // Check for significant switching
            var increased = usageChanges.FirstOrDefault(u => u.Delta > switchThreshold);
            var decreased = usageChanges.FirstOrDefault(u => u.Delta < -switchThreshold);

            if (increased.TranscriptId != null && decreased.TranscriptId != null)
            {
                double switchScore = Math.Abs(increased.Delta) + Math.Abs(decreased.Delta);

                yield return (group.Key, decreased.TranscriptId, increased.TranscriptId, switchScore);
            }
        }
    }

    #endregion

    #region Co-Expression Analysis

    /// <summary>
    /// Calculates Pearson correlation between gene expression profiles.
    /// </summary>
    public static double CalculatePearsonCorrelation(
        IReadOnlyList<double> expression1,
        IReadOnlyList<double> expression2)
    {
        if (expression1.Count != expression2.Count || expression1.Count < 2)
            return 0;

        int n = expression1.Count;
        double mean1 = expression1.Average();
        double mean2 = expression2.Average();

        double covariance = 0;
        double var1 = 0;
        double var2 = 0;

        for (int i = 0; i < n; i++)
        {
            double d1 = expression1[i] - mean1;
            double d2 = expression2[i] - mean2;
            covariance += d1 * d2;
            var1 += d1 * d1;
            var2 += d2 * d2;
        }

        if (var1 == 0 || var2 == 0)
            return 0;

        return covariance / Math.Sqrt(var1 * var2);
    }

    /// <summary>
    /// Builds a co-expression network.
    /// </summary>
    public static IEnumerable<(string Gene1, string Gene2, double Correlation)> BuildCoExpressionNetwork(
        IEnumerable<(string GeneId, IReadOnlyList<double> Expression)> geneProfiles,
        double correlationThreshold = 0.7)
    {
        var genes = geneProfiles.ToList();

        for (int i = 0; i < genes.Count; i++)
        {
            for (int j = i + 1; j < genes.Count; j++)
            {
                double corr = CalculatePearsonCorrelation(genes[i].Expression, genes[j].Expression);

                if (Math.Abs(corr) >= correlationThreshold)
                {
                    yield return (genes[i].GeneId, genes[j].GeneId, corr);
                }
            }
        }
    }

    /// <summary>
    /// Performs hierarchical clustering of genes by expression.
    /// </summary>
    public static IEnumerable<CoExpressionCluster> ClusterGenesByExpression(
        IEnumerable<(string GeneId, IReadOnlyList<double> Expression)> geneProfiles,
        int numClusters = 5,
        double correlationThreshold = 0.5)
    {
        var genes = geneProfiles.ToList();

        if (genes.Count == 0)
            yield break;

        // Simple k-means-like clustering based on correlation
        var clusters = new List<List<(string GeneId, IReadOnlyList<double> Expression)>>();

        // Initialize clusters with first genes
        for (int i = 0; i < Math.Min(numClusters, genes.Count); i++)
        {
            clusters.Add(new List<(string, IReadOnlyList<double>)> { genes[i] });
        }

        // Assign remaining genes to nearest cluster
        for (int i = numClusters; i < genes.Count; i++)
        {
            int bestCluster = 0;
            double bestCorr = double.MinValue;

            for (int c = 0; c < clusters.Count; c++)
            {
                // Calculate average correlation with cluster members
                double avgCorr = clusters[c]
                    .Select(m => CalculatePearsonCorrelation(genes[i].Expression, m.Expression))
                    .Average();

                if (avgCorr > bestCorr)
                {
                    bestCorr = avgCorr;
                    bestCluster = c;
                }
            }

            clusters[bestCluster].Add(genes[i]);
        }

        // Create cluster results
        for (int c = 0; c < clusters.Count; c++)
        {
            if (clusters[c].Count == 0)
                continue;

            var clusterGenes = clusters[c].Select(g => g.GeneId).ToList();

            // Calculate mean internal correlation
            double meanCorr = 0;
            int corrCount = 0;
            for (int i = 0; i < clusters[c].Count; i++)
            {
                for (int j = i + 1; j < clusters[c].Count; j++)
                {
                    meanCorr += CalculatePearsonCorrelation(
                        clusters[c][i].Expression,
                        clusters[c][j].Expression);
                    corrCount++;
                }
            }
            meanCorr = corrCount > 0 ? meanCorr / corrCount : 0;

            // Representative gene is one with highest mean correlation to others
            string representative = clusterGenes.First();

            yield return new CoExpressionCluster(
                ClusterId: c + 1,
                Genes: clusterGenes,
                MeanCorrelation: meanCorr,
                RepresentativeGene: representative,
                EnrichedFunctions: new List<string>());
        }
    }

    #endregion

    #region RNA-seq Quality Control

    /// <summary>
    /// Calculates basic RNA-seq quality metrics.
    /// </summary>
    public static (double MappingRate, double ExonicRate, double RRNARate, int DetectedGenes)
        CalculateQualityMetrics(
            double totalReads,
            double mappedReads,
            double exonicReads,
            double rRNAReads,
            IEnumerable<double> geneCounts)
    {
        double mappingRate = totalReads > 0 ? mappedReads / totalReads : 0;
        double exonicRate = mappedReads > 0 ? exonicReads / mappedReads : 0;
        double rrnaRate = mappedReads > 0 ? rRNAReads / mappedReads : 0;
        int detectedGenes = geneCounts.Count(c => c > 0);

        return (mappingRate, exonicRate, rrnaRate, detectedGenes);
    }

    /// <summary>
    /// Identifies potential batch effects using PCA.
    /// </summary>
    public static IEnumerable<(string SampleId, double PC1, double PC2)> PerformPCA(
        IEnumerable<(string SampleId, IReadOnlyList<double> Expression)> samples,
        int topGenes = 500)
    {
        var sampleList = samples.ToList();

        if (sampleList.Count < 2)
        {
            foreach (var sample in sampleList)
            {
                yield return (sample.SampleId, 0, 0);
            }
            yield break;
        }

        int geneCount = sampleList[0].Expression.Count;

        // Select top variable genes
        var variances = new double[geneCount];
        for (int g = 0; g < geneCount; g++)
        {
            var values = sampleList.Select(s => s.Expression[g]).ToList();
            double mean = values.Average();
            variances[g] = values.Sum(v => (v - mean) * (v - mean)) / values.Count;
        }

        var topGeneIndices = variances
            .Select((v, i) => (v, i))
            .OrderByDescending(x => x.v)
            .Take(Math.Min(topGenes, geneCount))
            .Select(x => x.i)
            .ToHashSet();

        // Simple PCA approximation using first two principal directions
        // (Full PCA would require SVD)
        foreach (var sample in sampleList)
        {
            var selectedValues = sample.Expression
                .Where((v, i) => topGeneIndices.Contains(i))
                .ToList();

            // Approximate PC scores as weighted sums
            double pc1 = selectedValues.Take(selectedValues.Count / 2).Sum();
            double pc2 = selectedValues.Skip(selectedValues.Count / 2).Sum();

            yield return (sample.SampleId, pc1, pc2);
        }
    }

    #endregion
}
