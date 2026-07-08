namespace Seqeron.Genomics.Population;

/// <summary>
/// Provides algorithms for population genetics analysis.
/// </summary>
public static class PopulationGeneticsAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Represents a genetic variant with allele frequencies.
    /// </summary>
    public readonly record struct Variant(
        string Id,
        string Chromosome,
        int Position,
        string ReferenceAllele,
        string AlternateAllele,
        double AlleleFrequency,
        int SampleCount);

    /// <summary>
    /// Represents population diversity statistics.
    /// </summary>
    public readonly record struct DiversityStatistics(
        double NucleotideDiversity,
        double WattersonTheta,
        double TajimasD,
        int SegregratingSites,
        int SampleSize,
        double HeterozygosityObserved,
        double HeterozygosityExpected);

    /// <summary>
    /// Represents F-statistics for population structure.
    /// </summary>
    public readonly record struct FStatistics(
        double Fst,
        double Fis,
        double Fit,
        string Population1,
        string Population2);

    /// <summary>
    /// Represents Hardy-Weinberg equilibrium test result.
    /// </summary>
    public readonly record struct HardyWeinbergResult(
        string VariantId,
        int ObservedAA,
        int ObservedAa,
        int Observedaa,
        double ExpectedAA,
        double ExpectedAa,
        double Expectedaa,
        double ChiSquare,
        double PValue,
        bool InEquilibrium);

    /// <summary>
    /// Represents linkage disequilibrium between two variants.
    /// </summary>
    public readonly record struct LinkageDisequilibrium(
        string Variant1,
        string Variant2,
        double DPrime,
        double RSquared,
        double Distance);

    /// <summary>
    /// Represents a haplotype block.
    /// </summary>
    public readonly record struct HaplotypeBlock(
        int Start,
        int End,
        IReadOnlyList<string> Variants,
        IReadOnlyList<(string Haplotype, double Frequency)> Haplotypes);

    /// <summary>
    /// Represents selection scan result.
    /// </summary>
    public readonly record struct SelectionSignal(
        string Region,
        int Start,
        int End,
        double Score,
        string TestType,
        double PValue,
        string Interpretation);

    /// <summary>
    /// Holds the unstandardized integrated haplotype score for a focal SNP together with
    /// its component integrated EHH areas and the derived allele frequency used for
    /// frequency-binned standardization (Voight et al. 2006).
    /// </summary>
    /// <param name="UnstandardizedIHS">ln(iHH_A / iHH_D) per Voight et al. (2006).</param>
    /// <param name="IhhAncestral">Integrated EHH (area under EHH curve) for the ancestral allele.</param>
    /// <param name="IhhDerived">Integrated EHH (area under EHH curve) for the derived allele.</param>
    /// <param name="DerivedAlleleFrequency">Frequency of the derived (core) allele in the sample.</param>
    public readonly record struct IhsResult(
        double UnstandardizedIHS,
        double IhhAncestral,
        double IhhDerived,
        double DerivedAlleleFrequency);

    /// <summary>
    /// Summarizes a genome-wide selection scan window: the fraction of SNPs whose
    /// absolute standardized iHS exceeds the extreme threshold (|iHS| &gt; 2) per
    /// Voight et al. (2006).
    /// </summary>
    public readonly record struct SelectionScanWindow(
        int WindowIndex,
        int SnpCount,
        int ExtremeCount,
        double ProportionExtreme);

    /// <summary>
    /// Represents ancestry proportion.
    /// </summary>
    public readonly record struct AncestryProportion(
        string IndividualId,
        IReadOnlyDictionary<string, double> Proportions);

    #endregion

    #region Allele Frequency Calculations

    /// <summary>
    /// Calculates allele frequencies from genotype counts.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when any genotype count is negative.
    /// </exception>
    public static (double MajorFreq, double MinorFreq) CalculateAlleleFrequencies(
        int homozygousMajor,
        int heterozygous,
        int homozygousMinor)
    {
        if (homozygousMajor < 0)
            throw new ArgumentOutOfRangeException(nameof(homozygousMajor), homozygousMajor, "Genotype count cannot be negative.");
        if (heterozygous < 0)
            throw new ArgumentOutOfRangeException(nameof(heterozygous), heterozygous, "Genotype count cannot be negative.");
        if (homozygousMinor < 0)
            throw new ArgumentOutOfRangeException(nameof(homozygousMinor), homozygousMinor, "Genotype count cannot be negative.");

        int totalAlleles = 2 * (homozygousMajor + heterozygous + homozygousMinor);

        if (totalAlleles == 0)
            return (0, 0);

        int majorAlleles = 2 * homozygousMajor + heterozygous;
        int minorAlleles = 2 * homozygousMinor + heterozygous;

        return ((double)majorAlleles / totalAlleles, (double)minorAlleles / totalAlleles);
    }

    /// <summary>
    /// Calculates minor allele frequency (MAF) from genotypes.
    /// </summary>
    public static double CalculateMAF(IEnumerable<int> genotypes)
    {
        var genotypeList = genotypes.ToList();

        if (genotypeList.Count == 0)
            return 0;

        // Genotypes: 0 = homozygous ref, 1 = heterozygous, 2 = homozygous alt
        int totalAlleles = genotypeList.Count * 2;
        int altAlleles = genotypeList.Sum();

        double altFreq = (double)altAlleles / totalAlleles;
        return Math.Min(altFreq, 1 - altFreq);
    }

    /// <summary>
    /// Filters variants by minor allele frequency.
    /// </summary>
    public static IEnumerable<Variant> FilterByMAF(
        IEnumerable<Variant> variants,
        double minMAF = 0.01,
        double maxMAF = 0.5)
    {
        foreach (var variant in variants)
        {
            double maf = Math.Min(variant.AlleleFrequency, 1 - variant.AlleleFrequency);

            if (maf >= minMAF && maf <= maxMAF)
            {
                yield return variant;
            }
        }
    }

    #endregion

    #region Diversity Statistics

    /// <summary>
    /// Calculates nucleotide diversity (π).
    /// </summary>
    public static double CalculateNucleotideDiversity(
        IEnumerable<IReadOnlyList<char>> sequences)
    {
        var seqList = sequences.ToList();

        if (seqList.Count < 2)
            return 0;

        int n = seqList.Count;
        int length = seqList[0].Count;
        double totalDiff = 0;
        int comparisons = 0;

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                int diffs = 0;
                for (int k = 0; k < length; k++)
                {
                    if (seqList[i][k] != seqList[j][k])
                        diffs++;
                }
                totalDiff += diffs;
                comparisons++;
            }
        }

        return totalDiff / (comparisons * length);
    }

    /// <summary>
    /// Calculates Watterson's theta estimator.
    /// </summary>
    public static double CalculateWattersonTheta(int segregatingSites, int sampleSize, int sequenceLength)
    {
        if (sampleSize < 2 || sequenceLength <= 0)
            return 0;

        // Harmonic number a1
        double a1 = 0;
        for (int i = 1; i < sampleSize; i++)
        {
            a1 += 1.0 / i;
        }

        return (double)segregatingSites / (a1 * sequenceLength);
    }

    /// <summary>
    /// Calculates Tajima's D statistic.
    /// Formula (Tajima 1989, Wikipedia):
    ///   D = (k̂ − S/a₁) / √(e₁·S + e₂·S·(S−1))
    /// where k̂ = average number of pairwise differences (NOT per-site).
    /// </summary>
    /// <param name="averagePairwiseDifferences">
    /// k̂ — average number of nucleotide differences per pair of sequences.
    /// Equals π (per-site nucleotide diversity) × L (sequence length).
    /// </param>
    /// <param name="segregatingSites">S — number of polymorphic positions.</param>
    /// <param name="sampleSize">n — number of sequences in the sample (requires n ≥ 3).</param>
    public static double CalculateTajimasD(
        double averagePairwiseDifferences,
        int segregatingSites,
        int sampleSize)
    {
        if (segregatingSites == 0 || sampleSize < 3)
            return 0;

        int n = sampleSize;

        // Calculate harmonic numbers
        double a1 = 0, a2 = 0;
        for (int i = 1; i < n; i++)
        {
            a1 += 1.0 / i;
            a2 += 1.0 / (i * i);
        }

        // Watterson estimate of expected pairwise differences: S / a₁
        double wattersonEstimate = (double)segregatingSites / a1;

        // Calculate constants (Tajima 1989)
        double b1 = (n + 1.0) / (3 * (n - 1));
        double b2 = 2.0 * (n * n + n + 3) / (9 * n * (n - 1));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);

        // Calculate variance of d = k̂ − S/a₁
        double variance = e1 * segregatingSites + e2 * segregatingSites * (segregatingSites - 1);

        if (variance <= 0)
            return 0;

        // Tajima's D = d / √V̂(d)
        double d = averagePairwiseDifferences - wattersonEstimate;
        return d / Math.Sqrt(variance);
    }

    /// <summary>
    /// Calculates comprehensive diversity statistics.
    /// </summary>
    public static DiversityStatistics CalculateDiversityStatistics(
        IEnumerable<IReadOnlyList<char>> sequences)
    {
        var seqList = sequences.ToList();

        if (seqList.Count < 2)
        {
            return new DiversityStatistics(0, 0, 0, 0, seqList.Count, 0, 0);
        }

        int n = seqList.Count;
        int length = seqList[0].Count;

        // Count segregating sites
        int segregatingSites = 0;
        for (int pos = 0; pos < length; pos++)
        {
            char first = seqList[0][pos];
            if (seqList.Any(s => s[pos] != first))
                segregatingSites++;
        }

        double pi = CalculateNucleotideDiversity(seqList);
        double theta = CalculateWattersonTheta(segregatingSites, n, length);

        // k̂ = average pairwise differences (unnormalized) = π × L
        double kHat = pi * length;
        double tajD = CalculateTajimasD(kHat, segregatingSites, n);

        // Calculate heterozygosity
        double hetObs = CalculateObservedHeterozygosity(seqList);
        double hetExp = CalculateExpectedHeterozygosity(seqList);

        return new DiversityStatistics(
            NucleotideDiversity: pi,
            WattersonTheta: theta,
            TajimasD: tajD,
            SegregratingSites: segregatingSites,
            SampleSize: n,
            HeterozygosityObserved: hetObs,
            HeterozygosityExpected: hetExp);
    }

    /// <summary>
    /// Nei's (1978) unbiased gene diversity per site, serving as the haploid analogue
    /// of observed heterozygosity. Formula: (n/(n-1)) × (1 − Σp_i²) averaged over sites.
    /// Source: Nei M. (1978) "Estimation of average heterozygosity and genetic distance
    /// from a small number of individuals", Genetics 89(3):583–590.
    /// </summary>
    private static double CalculateObservedHeterozygosity(List<IReadOnlyList<char>> sequences)
    {
        if (sequences.Count < 2)
            return 0;

        int n = sequences.Count;
        int length = sequences[0].Count;
        double totalHet = 0;

        for (int pos = 0; pos < length; pos++)
        {
            var alleleCounts = sequences
                .GroupBy(s => s[pos])
                .ToDictionary(g => g.Key, g => g.Count());

            double sumPiSquared = alleleCounts.Values
                .Select(c => (double)c / n)
                .Select(p => p * p)
                .Sum();

            // Nei's unbiased estimator: n/(n-1) × (1 - Σp²)
            totalHet += (double)n / (n - 1) * (1 - sumPiSquared);
        }

        return totalHet / length;
    }

    /// <summary>
    /// Basic gene diversity per site (expected heterozygosity).
    /// Formula: (1 − Σp_i²) averaged over sites.
    /// Source: Wikipedia — Zygosity, "Heterozygosity in population genetics".
    /// H_e = 1 − Σ f_i² where f_i = allele frequency.
    /// </summary>
    private static double CalculateExpectedHeterozygosity(List<IReadOnlyList<char>> sequences)
    {
        if (sequences.Count < 2)
            return 0;

        int length = sequences[0].Count;
        double totalHet = 0;

        for (int pos = 0; pos < length; pos++)
        {
            var alleleCounts = sequences
                .GroupBy(s => s[pos])
                .ToDictionary(g => g.Key, g => g.Count());

            int n = sequences.Count;
            double sumPiSquared = alleleCounts.Values
                .Select(c => (double)c / n)
                .Select(p => p * p)
                .Sum();

            totalHet += 1 - sumPiSquared;
        }

        return totalHet / length;
    }

    #endregion

    #region Hardy-Weinberg Equilibrium

    /// <summary>
    /// Tests Hardy-Weinberg equilibrium for a variant.
    /// </summary>
    public static HardyWeinbergResult TestHardyWeinberg(
        string variantId,
        int observedAA,
        int observedAa,
        int observedaa,
        double significanceLevel = 0.05)
    {
        int n = observedAA + observedAa + observedaa;

        if (n == 0)
        {
            return new HardyWeinbergResult(variantId, 0, 0, 0, 0, 0, 0, 0, 1, true);
        }

        // Calculate allele frequencies
        double p = (2.0 * observedAA + observedAa) / (2.0 * n);
        double q = 1 - p;

        // Expected counts under HWE
        double expectedAA = p * p * n;
        double expectedAa = 2 * p * q * n;
        double expectedaa = q * q * n;

        // Chi-square test
        double chiSquare = 0;

        if (expectedAA > 0)
            chiSquare += Math.Pow(observedAA - expectedAA, 2) / expectedAA;
        if (expectedAa > 0)
            chiSquare += Math.Pow(observedAa - expectedAa, 2) / expectedAa;
        if (expectedaa > 0)
            chiSquare += Math.Pow(observedaa - expectedaa, 2) / expectedaa;

        // P-value (1 degree of freedom)
        double pValue = 1 - ChiSquareCDF(chiSquare, 1);

        return new HardyWeinbergResult(
            VariantId: variantId,
            ObservedAA: observedAA,
            ObservedAa: observedAa,
            Observedaa: observedaa,
            ExpectedAA: expectedAA,
            ExpectedAa: expectedAa,
            Expectedaa: expectedaa,
            ChiSquare: chiSquare,
            PValue: pValue,
            InEquilibrium: pValue >= significanceLevel);
    }

    private static double ChiSquareCDF(double x, int df)
    {
        if (x < 0)
            return 0;

        if (x == 0)
            return 0;

        // For df=1, use a more stable formula
        // CDF = 2 * Phi(sqrt(x)) - 1 where Phi is standard normal CDF
        // Or use regularized incomplete gamma function
        return RegularizedGammaP(df / 2.0, x / 2.0);
    }

    /// <summary>
    /// Computes the regularized lower incomplete gamma function P(a, x).
    /// P(a, x) = γ(a, x) / Γ(a) where γ is the lower incomplete gamma function.
    /// </summary>
    private static double RegularizedGammaP(double a, double x)
    {
        if (x < 0 || a <= 0)
            return 0;

        if (x == 0)
            return 0;

        // Use series expansion for small x, continued fraction for large x
        if (x < a + 1)
        {
            return GammaSeriesP(a, x);
        }
        else
        {
            return 1.0 - GammaContinuedFractionQ(a, x);
        }
    }

    /// <summary>
    /// Series expansion for regularized incomplete gamma P(a, x).
    /// </summary>
    private static double GammaSeriesP(double a, double x)
    {
        const int maxIterations = 200;
        const double epsilon = 1e-14;

        double logGammaA = LogGamma(a);
        double sum = 1.0 / a;
        double term = sum;

        for (int n = 1; n < maxIterations; n++)
        {
            term *= x / (a + n);
            sum += term;
            if (Math.Abs(term) < Math.Abs(sum) * epsilon)
                break;
        }

        return sum * Math.Exp(-x + a * Math.Log(x) - logGammaA);
    }

    /// <summary>
    /// Continued fraction for regularized incomplete gamma Q(a, x) = 1 - P(a, x).
    /// </summary>
    private static double GammaContinuedFractionQ(double a, double x)
    {
        const int maxIterations = 200;
        const double epsilon = 1e-14;
        const double tiny = 1e-30;

        double logGammaA = LogGamma(a);
        double b = x + 1.0 - a;
        double c = 1.0 / tiny;
        double d = 1.0 / b;
        double h = d;

        for (int n = 1; n < maxIterations; n++)
        {
            double an = -n * (n - a);
            b += 2.0;
            d = an * d + b;
            if (Math.Abs(d) < tiny) d = tiny;
            c = b + an / c;
            if (Math.Abs(c) < tiny) c = tiny;
            d = 1.0 / d;
            double delta = d * c;
            h *= delta;
            if (Math.Abs(delta - 1.0) < epsilon)
                break;
        }

        return Math.Exp(-x + a * Math.Log(x) - logGammaA) * h;
    }

    /// <summary>
    /// Computes ln(Γ(x)) using Lanczos approximation.
    /// </summary>
    private static double LogGamma(double x)
    {
        if (x <= 0)
            return double.PositiveInfinity;

        // Lanczos approximation coefficients
        double[] c =
        {
            76.18009172947146,
            -86.50532032941677,
            24.01409824083091,
            -1.231739572450155,
            0.1208650973866179e-2,
            -0.5395239384953e-5
        };

        double y = x;
        double tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.Log(tmp);

        double ser = 1.000000000190015;
        for (int j = 0; j < 6; j++)
        {
            y += 1;
            ser += c[j] / y;
        }

        return -tmp + Math.Log(2.5066282746310005 * ser / x);
    }

    #endregion

    #region Population Structure (F-statistics)

    /// <summary>
    /// Calculates Wright's variance-based Fst between populations.
    /// Formula: Fst = σ²_S / p̄(1-p̄) where σ²_S is the weighted variance of allele frequencies
    /// among subpopulations and p̄(1-p̄) is the expected heterozygosity.
    /// Source: Wright (1965) Evolution 19:395-420; Wikipedia: Fixation index §Definition.
    /// </summary>
    public static double CalculateFst(
        IEnumerable<(double AlleleFreq, int SampleSize)> population1,
        IEnumerable<(double AlleleFreq, int SampleSize)> population2)
    {
        var pop1 = population1.ToList();
        var pop2 = population2.ToList();

        if (pop1.Count == 0 || pop2.Count == 0)
            return 0;

        if (pop1.Count != pop2.Count)
            throw new ArgumentException(
                $"The two populations' per-locus allele frequency counts must match; " +
                $"got {pop1.Count} and {pop2.Count}.",
                nameof(population2));

        double numerator = 0;
        double denominator = 0;

        for (int i = 0; i < pop1.Count; i++)
        {
            double p1 = pop1[i].AlleleFreq;
            double p2 = pop2[i].AlleleFreq;
            int n1 = pop1[i].SampleSize;
            int n2 = pop2[i].SampleSize;

            double pBar = (n1 * p1 + n2 * p2) / (n1 + n2);
            double variance = ((p1 - pBar) * (p1 - pBar) * n1 +
                               (p2 - pBar) * (p2 - pBar) * n2) / (n1 + n2);

            double het = pBar * (1 - pBar);

            numerator += variance;
            denominator += het;
        }

        return denominator > 0 ? numerator / denominator : 0;
    }

    /// <summary>
    /// Calculates pairwise Fst matrix for multiple populations.
    /// </summary>
    public static double[,] CalculatePairwiseFst(
        IEnumerable<(string PopulationId, IReadOnlyList<(double AlleleFreq, int SampleSize)> Variants)> populations)
    {
        var popList = populations.ToList();
        int n = popList.Count;
        var fstMatrix = new double[n, n];

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double fst = CalculateFst(popList[i].Variants, popList[j].Variants);
                fstMatrix[i, j] = fst;
                fstMatrix[j, i] = fst;
            }
        }

        return fstMatrix;
    }

    /// <summary>
    /// Calculates F-statistics (Fis, Fit, Fst).
    /// </summary>
    public static FStatistics CalculateFStatistics(
        string pop1Name,
        string pop2Name,
        IEnumerable<(int HetObs1, int N1, int HetObs2, int N2, double AlleleFreq1, double AlleleFreq2)> variantData)
    {
        var data = variantData.ToList();

        if (data.Count == 0)
            return new FStatistics(0, 0, 0, pop1Name, pop2Name);

        double totalHetObs = 0;
        double totalHetExp = 0;
        double totalHetTotal = 0;
        int totalN = 0;

        foreach (var (hetObs1, n1, hetObs2, n2, p1, p2) in data)
        {
            double pBar = (n1 * p1 + n2 * p2) / (n1 + n2);

            totalHetObs += hetObs1 + hetObs2;
            totalHetExp += 2 * p1 * (1 - p1) * n1 + 2 * p2 * (1 - p2) * n2;
            totalHetTotal += 2 * pBar * (1 - pBar) * (n1 + n2);
            totalN += n1 + n2;
        }

        double hi = totalN > 0 ? totalHetObs / totalN : 0;
        double hs = totalN > 0 ? totalHetExp / totalN : 0;
        double ht = totalN > 0 ? totalHetTotal / totalN : 0;

        double fis = hs > 0 ? 1 - hi / hs : 0;
        double fit = ht > 0 ? 1 - hi / ht : 0;
        double fst = ht > 0 ? 1 - hs / ht : 0;

        return new FStatistics(
            Fst: fst,
            Fis: fis,
            Fit: fit,
            Population1: pop1Name,
            Population2: pop2Name);
    }

    #endregion

    #region Linkage Disequilibrium

    /// <summary>
    /// Calculates linkage disequilibrium between two variants.
    /// 
    /// r² is computed as the squared Pearson correlation of genotype values (0, 1, 2).
    /// From Wikipedia (LD for diploid frequencies): the diploid correlation R_AB
    /// equals the haplotype-level r_AB (Hill &amp; Robertson 1968, Wright 1933).
    /// 
    /// D' uses D estimated from the diploid genotype covariance: D = Cov(X₁,X₂)/2,
    /// then normalized per Lewontin (1964): D' = D / D_max, clamped to [0, 1].
    /// </summary>
    public static LinkageDisequilibrium CalculateLD(
        string variant1Id,
        string variant2Id,
        IEnumerable<(int Geno1, int Geno2)> genotypes,
        int distance)
    {
        var genoList = genotypes.ToList();

        if (genoList.Count == 0)
        {
            return new LinkageDisequilibrium(variant1Id, variant2Id, 0, 0, distance);
        }

        int n = genoList.Count;

        // Calculate allele frequencies (genotype 0=AA, 1=AB, 2=BB)
        // Allele frequency p = frequency of B allele
        double p1 = genoList.Sum(g => g.Geno1) / (2.0 * n);
        double p2 = genoList.Sum(g => g.Geno2) / (2.0 * n);
        double q1 = 1 - p1;
        double q2 = 1 - p2;

        // r² from genotypes: squared Pearson correlation of genotype values (0, 1, 2).
        // From Wikipedia (LD for diploid frequencies): the diploid correlation R_AB
        // equals the haplotype correlation r_AB, so r² = Cor(X₁,X₂)².
        double mean1 = genoList.Average(g => (double)g.Geno1);
        double mean2 = genoList.Average(g => (double)g.Geno2);
        double cov = genoList.Sum(g => (g.Geno1 - mean1) * (g.Geno2 - mean2)) / n;
        double var1 = genoList.Sum(g => Math.Pow(g.Geno1 - mean1, 2)) / n;
        double var2 = genoList.Sum(g => Math.Pow(g.Geno2 - mean2, 2)) / n;

        double rSquared = (var1 > 0 && var2 > 0) ? (cov * cov) / (var1 * var2) : 0;
        rSquared = Math.Clamp(rSquared, 0.0, 1.0);  // squared correlation ∈ [0,1]; clamp FP rounding noise (mirrors the D' clamp below)

        // D estimated from diploid genotype covariance.
        // From Wikipedia (LD for diploid frequencies): Cov_diploid(X₁,X₂) = 2D
        // with 0/1/2 encoding, therefore D = Cov(X₁,X₂) / 2.
        double d = cov / 2;

        double dMax = d >= 0
            ? Math.Min(p1 * q2, q1 * p2)
            : Math.Min(p1 * p2, q1 * q2);

        double dPrime = dMax > 1e-10 ? Math.Abs(d) / dMax : 0;
        dPrime = Math.Min(dPrime, 1.0);  // Ensure D' ≤ 1

        return new LinkageDisequilibrium(
            Variant1: variant1Id,
            Variant2: variant2Id,
            DPrime: dPrime,
            RSquared: rSquared,
            Distance: distance);
    }

    /// <summary>
    /// Identifies haplotype blocks using adjacent-pair r² threshold.
    /// Simplified Gabriel et al. (2002) method: consecutive variants with r² ≥ threshold form a block.
    /// </summary>
    public static IEnumerable<HaplotypeBlock> FindHaplotypeBlocks(
        IEnumerable<(string VariantId, int Position, IReadOnlyList<int> Genotypes)> variants,
        double ldThreshold = 0.7)
    {
        var variantList = variants.OrderBy(v => v.Position).ToList();

        if (variantList.Count < 2)
            yield break;

        int blockStart = variantList[0].Position;
        var blockVariants = new List<string> { variantList[0].VariantId };

        for (int i = 1; i < variantList.Count; i++)
        {
            var prev = variantList[i - 1];
            var curr = variantList[i];

            // Calculate LD
            var genoPairs = prev.Genotypes
                .Zip(curr.Genotypes, (g1, g2) => (g1, g2))
                .ToList();

            var ld = CalculateLD(
                prev.VariantId,
                curr.VariantId,
                genoPairs,
                curr.Position - prev.Position);

            if (ld.RSquared >= ldThreshold)
            {
                blockVariants.Add(curr.VariantId);
            }
            else
            {
                // End current block
                if (blockVariants.Count >= 2)
                {
                    yield return new HaplotypeBlock(
                        Start: blockStart,
                        End: prev.Position,
                        Variants: blockVariants.ToList(),
                        Haplotypes: new List<(string, double)>());
                }

                // Start new block
                blockStart = curr.Position;
                blockVariants = new List<string> { curr.VariantId };
            }
        }

        // Final block
        if (blockVariants.Count >= 2)
        {
            yield return new HaplotypeBlock(
                Start: blockStart,
                End: variantList[^1].Position,
                Variants: blockVariants.ToList(),
                Haplotypes: new List<(string, double)>());
        }
    }

    #endregion

    #region Selection Tests

    /// <summary>
    /// Calculates integrated haplotype score (iHS).
    /// </summary>
    public static double CalculateIHS(
        IReadOnlyList<double> ehh0,
        IReadOnlyList<double> ehh1,
        IReadOnlyList<int> positions)
    {
        if (ehh0.Count != ehh1.Count || ehh0.Count != positions.Count || ehh0.Count < 2)
            return 0;

        // Integrate EHH for ancestral (0) and derived (1) alleles
        double ihh0 = 0, ihh1 = 0;

        for (int i = 1; i < positions.Count; i++)
        {
            double dist = positions[i] - positions[i - 1];
            ihh0 += (ehh0[i - 1] + ehh0[i]) / 2 * dist;
            ihh1 += (ehh1[i - 1] + ehh1[i]) / 2 * dist;
        }

        if (ihh0 <= 0 || ihh1 <= 0)
            return 0;

        return Math.Log(ihh1 / ihh0);
    }

    /// <summary>
    /// Scans for selection signals using multiple tests.
    /// </summary>
    public static IEnumerable<SelectionSignal> ScanForSelection(
        IEnumerable<(string Region, int Start, int End, double TajimaD, double Fst, double IHS)> regions,
        double tajimaDThreshold = -2.0,
        double fstThreshold = 0.25,
        double ihsThreshold = 2.0)
    {
        foreach (var (region, start, end, tajD, fst, ihs) in regions)
        {
            // Positive selection signals
            if (tajD < tajimaDThreshold)
            {
                yield return new SelectionSignal(
                    Region: region,
                    Start: start,
                    End: end,
                    Score: tajD,
                    TestType: "TajimasD",
                    PValue: EstimateSelectionPValue(tajD, "TajimasD"),
                    Interpretation: "Possible positive/purifying selection (excess rare variants)");
            }

            if (fst > fstThreshold)
            {
                yield return new SelectionSignal(
                    Region: region,
                    Start: start,
                    End: end,
                    Score: fst,
                    TestType: "Fst",
                    PValue: EstimateSelectionPValue(fst, "Fst"),
                    Interpretation: "Possible local adaptation (high differentiation)");
            }

            if (Math.Abs(ihs) > ihsThreshold)
            {
                yield return new SelectionSignal(
                    Region: region,
                    Start: start,
                    End: end,
                    Score: ihs,
                    TestType: "iHS",
                    PValue: EstimateSelectionPValue(ihs, "iHS"),
                    Interpretation: ihs > 0
                        ? "Positive selection on derived allele"
                        : "Positive selection on ancestral allele");
            }
        }
    }

    private static double EstimateSelectionPValue(double score, string testType)
    {
        // Simplified p-value estimation using normal approximation
        double z = testType switch
        {
            "TajimasD" => Math.Abs(score),
            "Fst" => score * 10, // Scale Fst
            "iHS" => Math.Abs(score),
            _ => Math.Abs(score)
        };

        return 2 * (1 - StatisticsHelper.NormalCDF(z));
    }

    // ----- Canonical iHS pipeline (Voight et al. 2006; Sabeti et al. 2002; Szpiech & Hernandez 2014) -----

    // EHH integration is truncated where EHH first drops below this cutoff, per Voight et al.
    // (2006) Materials and Methods ("nearest points ... where the EHH drops below 0.05") and the
    // rehh default parameter limehh = 0.05 (Gautier et al. 2017).
    private const double EhhIntegrationCutoff = 0.05;

    // A standardized iHS is "extreme" when its magnitude exceeds 2; genome-wide selection is
    // quantified by the proportion of SNPs with |iHS| > 2 (Voight et al. 2006, Materials and
    // Methods: "quantified by the proportion of SNPs with |iHS| > 2").
    private const double ExtremeIhsThreshold = 2.0;

    /// <summary>
    /// Computes Extended Haplotype Homozygosity (EHH) for the chromosomes carrying a given core
    /// allele, extended over a marker window. EHH is the probability that two randomly chosen
    /// core-carrying chromosomes are identical over the window:
    /// EHH = Σ_h C(n_h, 2) / C(n_c, 2), where n_h is the count of each distinct extended haplotype
    /// and n_c the number of core-carrying chromosomes (Sabeti et al. 2002; Szpiech &amp; Hernandez
    /// 2014, Eq. 3).
    /// </summary>
    /// <param name="extendedHaplotypes">
    /// The extended haplotype string (alleles over the window) for each core-carrying chromosome.
    /// </param>
    /// <returns>
    /// EHH in [0, 1]. Returns 1 for a single chromosome (a sample of one is trivially homozygous);
    /// returns 0 for an empty sample.
    /// </returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="extendedHaplotypes"/> is null.</exception>
    public static double CalculateEhh(IReadOnlyList<string> extendedHaplotypes)
    {
        ArgumentNullException.ThrowIfNull(extendedHaplotypes);

        int nc = extendedHaplotypes.Count;
        if (nc == 0)
            return 0;
        if (nc == 1)
            return 1;

        var counts = new Dictionary<string, int>(StringComparer.Ordinal);
        foreach (string hap in extendedHaplotypes)
            counts[hap] = counts.TryGetValue(hap, out int c) ? c + 1 : 1;

        double numerator = 0;
        foreach (int nh in counts.Values)
            numerator += Choose2(nh);

        return numerator / Choose2(nc);
    }

    /// <summary>
    /// Computes the unstandardized integrated haplotype score (iHS) at a focal SNP from phased
    /// haplotypes. EHH is tracked outward in both directions for the ancestral (0) and derived (1)
    /// core alleles, the EHH curves are integrated by the trapezoidal rule against marker positions
    /// (truncating where EHH first drops below 0.05), and the score is
    /// unstandardized iHS = ln(iHH_A / iHH_D) (Voight et al. 2006).
    /// </summary>
    /// <param name="haplotypes">
    /// Phased haplotypes; each string holds one allele per marker ('0' = ancestral, '1' = derived).
    /// All strings must have the same length, equal to <paramref name="positions"/>.Count.
    /// </param>
    /// <param name="positions">Chromosomal (or genetic) positions of the markers, strictly increasing.</param>
    /// <param name="coreIndex">Index of the focal SNP within each haplotype.</param>
    /// <returns>
    /// The iHS result. <see cref="IhsResult.UnstandardizedIHS"/> follows the Voight et al. (2006)
    /// sign convention ln(iHH_A / iHH_D): negative values indicate unusually long haplotypes on the
    /// derived allele (candidate positive selection on the derived allele).
    /// </returns>
    /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
    /// <exception cref="ArgumentException">Thrown when haplotype lengths are inconsistent or the core is monomorphic.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="coreIndex"/> is out of range.</exception>
    public static IhsResult CalculateIHS(
        IReadOnlyList<string> haplotypes,
        IReadOnlyList<int> positions,
        int coreIndex)
    {
        ArgumentNullException.ThrowIfNull(haplotypes);
        ArgumentNullException.ThrowIfNull(positions);
        if (haplotypes.Count == 0)
            throw new ArgumentException("At least one haplotype is required.", nameof(haplotypes));
        if (coreIndex < 0 || coreIndex >= positions.Count)
            throw new ArgumentOutOfRangeException(nameof(coreIndex));

        int markerCount = positions.Count;
        foreach (string hap in haplotypes)
        {
            if (hap is null)
                throw new ArgumentException("Haplotype strings cannot be null.", nameof(haplotypes));
            if (hap.Length != markerCount)
                throw new ArgumentException("Each haplotype must have one allele per position.", nameof(haplotypes));
        }

        var ancestral = new List<int>();
        var derived = new List<int>();
        for (int h = 0; h < haplotypes.Count; h++)
        {
            char allele = haplotypes[h][coreIndex];
            if (allele == '0')
                ancestral.Add(h);
            else if (allele == '1')
                derived.Add(h);
            else
                throw new ArgumentException("Core alleles must be '0' (ancestral) or '1' (derived).", nameof(haplotypes));
        }

        if (ancestral.Count == 0 || derived.Count == 0)
            throw new ArgumentException("The focal SNP must be polymorphic (both ancestral and derived alleles present).", nameof(haplotypes));

        double ihhAncestral = IntegrateEhh(haplotypes, positions, coreIndex, ancestral);
        double ihhDerived = IntegrateEhh(haplotypes, positions, coreIndex, derived);

        double derivedFreq = (double)derived.Count / (ancestral.Count + derived.Count);
        double unstandardized = (ihhAncestral <= 0 || ihhDerived <= 0)
            ? 0
            : Math.Log(ihhAncestral / ihhDerived);

        return new IhsResult(unstandardized, ihhAncestral, ihhDerived, derivedFreq);
    }

    /// <summary>
    /// Standardizes unstandardized iHS scores within derived-allele-frequency bins so the result is
    /// approximately standard normal and comparable across allele frequencies:
    /// iHS = (x − E_p[x]) / SD_p[x], where the expectation and standard deviation are taken over the
    /// empirical distribution of SNPs whose derived allele frequency p falls in the same bin
    /// (Voight et al. 2006).
    /// </summary>
    /// <param name="scores">Unstandardized iHS values with their derived allele frequencies.</param>
    /// <param name="binCount">Number of equal-width frequency bins over (0, 1). Default 20 (bin width 0.05).</param>
    /// <returns>Standardized iHS values aligned with the input order.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="scores"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="binCount"/> &lt; 1.</exception>
    public static IReadOnlyList<double> StandardizeIHS(
        IReadOnlyList<(double Unstandardized, double DerivedAlleleFrequency)> scores,
        int binCount = 20)
    {
        ArgumentNullException.ThrowIfNull(scores);
        if (binCount < 1)
            throw new ArgumentOutOfRangeException(nameof(binCount), binCount, "Bin count must be at least 1.");

        // Group SNP indices by frequency bin.
        var bins = new Dictionary<int, List<int>>();
        for (int i = 0; i < scores.Count; i++)
        {
            double p = scores[i].DerivedAlleleFrequency;
            int bin = (int)(p * binCount);
            if (bin >= binCount)
                bin = binCount - 1; // include p == 1.0 in the top bin
            if (bin < 0)
                bin = 0;
            if (!bins.TryGetValue(bin, out var list))
                bins[bin] = list = new List<int>();
            list.Add(i);
        }

        var result = new double[scores.Count];
        foreach (var indices in bins.Values)
        {
            double mean = 0;
            foreach (int idx in indices)
                mean += scores[idx].Unstandardized;
            mean /= indices.Count;

            double variance = 0;
            foreach (int idx in indices)
            {
                double d = scores[idx].Unstandardized - mean;
                variance += d * d;
            }
            // Sample standard deviation of the empirical bin distribution.
            double sd = indices.Count > 1 ? Math.Sqrt(variance / (indices.Count - 1)) : 0;

            foreach (int idx in indices)
                result[idx] = sd > 0 ? (scores[idx].Unstandardized - mean) / sd : 0;
        }

        return result;
    }

    /// <summary>
    /// Performs a genome-wide selection scan over fixed windows of consecutive SNPs, quantifying
    /// each window by the proportion of SNPs with |iHS| &gt; 2 (Voight et al. 2006, Materials and
    /// Methods). SNPs are assumed to be supplied in genomic order.
    /// </summary>
    /// <param name="standardizedScores">Standardized iHS value for each SNP, in genomic order.</param>
    /// <param name="windowSize">Number of consecutive SNPs per window (Voight et al. used 50).</param>
    /// <returns>One <see cref="SelectionScanWindow"/> per non-overlapping window.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="standardizedScores"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="windowSize"/> &lt; 1.</exception>
    public static IEnumerable<SelectionScanWindow> ScanForSelection(
        IReadOnlyList<double> standardizedScores,
        int windowSize = 50)
    {
        ArgumentNullException.ThrowIfNull(standardizedScores);
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1.");

        int windowIndex = 0;
        for (int start = 0; start < standardizedScores.Count; start += windowSize)
        {
            int end = Math.Min(start + windowSize, standardizedScores.Count);
            int snpCount = end - start;
            int extreme = 0;
            for (int i = start; i < end; i++)
            {
                if (Math.Abs(standardizedScores[i]) > ExtremeIhsThreshold)
                    extreme++;
            }

            yield return new SelectionScanWindow(
                WindowIndex: windowIndex++,
                SnpCount: snpCount,
                ExtremeCount: extreme,
                ProportionExtreme: (double)extreme / snpCount);
        }
    }

    /// <summary>
    /// Integrates the EHH decay curve for one core-allele subset outward in both directions from the
    /// focal SNP using the trapezoidal rule, summing the areas to the left and right and truncating
    /// each direction where EHH first drops below the 0.05 cutoff (Voight et al. 2006; Szpiech &amp;
    /// Hernandez 2014, Eq. 4/7).
    /// </summary>
    private static double IntegrateEhh(
        IReadOnlyList<string> haplotypes,
        IReadOnlyList<int> positions,
        int coreIndex,
        IReadOnlyList<int> coreChromosomes)
    {
        double area = 0;
        area += IntegrateDirection(haplotypes, positions, coreIndex, coreChromosomes, step: +1);
        area += IntegrateDirection(haplotypes, positions, coreIndex, coreChromosomes, step: -1);
        return area;
    }

    private static double IntegrateDirection(
        IReadOnlyList<string> haplotypes,
        IReadOnlyList<int> positions,
        int coreIndex,
        IReadOnlyList<int> coreChromosomes,
        int step)
    {
        // At the core itself every core-carrying chromosome is identical, so EHH = 1.
        double previousEhh = 1;
        int previousPos = positions[coreIndex];
        double area = 0;

        for (int marker = coreIndex + step; marker >= 0 && marker < positions.Count; marker += step)
        {
            int lo = Math.Min(coreIndex, marker);
            int hi = Math.Max(coreIndex, marker);
            var window = new List<string>(coreChromosomes.Count);
            foreach (int chrom in coreChromosomes)
                window.Add(haplotypes[chrom].Substring(lo, hi - lo + 1));

            double ehh = CalculateEhh(window);
            int pos = positions[marker];

            // Add the trapezoid up to this marker, then stop once EHH has decayed below the cutoff.
            area += (previousEhh + ehh) / 2 * Math.Abs(pos - previousPos);

            if (ehh < EhhIntegrationCutoff)
                break;

            previousEhh = ehh;
            previousPos = pos;
        }

        return area;
    }

    /// <summary>
    /// Number of unordered pairs C(n, 2) = n(n − 1) / 2.
    /// </summary>
    private static double Choose2(int n) => n * (n - 1) / 2.0;

    #endregion

    #region Ancestry Analysis

    // ADMIXTURE default stopping criterion ε for the convergence rule (Equation 5):
    // declare convergence once the log-likelihood gain falls below ε.
    // Alexander, Novembre & Lange (2009), Genome Research 19(9):1655–1664, Eq. 5
    // ("we choose ε = 10⁻⁴ as the default stopping criterion in ADMIXTURE").
    private const double AncestryLogLikelihoodTolerance = 1e-4;

    /// <summary>
    /// Estimates individual ancestry (admixture) proportions by maximum likelihood given
    /// <b>fixed</b> reference-population allele frequencies, i.e. the supervised / projection
    /// mode of ADMIXTURE. For each individual the ancestry vector q (one fraction per reference
    /// population) is found with the FRAPPE expectation-maximization update on the binomial
    /// admixture log-likelihood of Alexander, Novembre &amp; Lange (2009).
    /// </summary>
    /// <param name="individuals">
    /// Individuals to estimate. <c>Genotypes[j]</c> is the observed number of copies of allele 1
    /// at SNP <c>j</c> (0, 1, or 2); any other value is treated as missing and that SNP is skipped
    /// for the individual. The genotype length must equal the reference-panel SNP count.
    /// </param>
    /// <param name="referencePops">
    /// Reference populations with fixed allele-1 frequencies (one per SNP, in [0,1]); all panels
    /// must cover the same SNP set in the same order.
    /// </param>
    /// <param name="maxIterations">Maximum EM iterations per individual (default 100).</param>
    /// <returns>
    /// One <see cref="AncestryProportion"/> per valid individual, whose <c>Proportions</c> are
    /// keyed by reference-population id and sum to 1.
    /// </returns>
    /// <remarks>
    /// Model (Alexander et al. 2009): q_ik is the fraction of individual i's genome from
    /// population k; f_kj is the allele-1 frequency in population k at SNP j. With F fixed, the
    /// EM update (their Equation 4) is
    /// q_ik^{n+1} = (1/2J) Σ_j [ g_ij·a_ijk + (2−g_ij)·b_ijk ], where
    /// a_ijk = q_ik f_kj / Σ_m q_im f_mj and b_ijk = q_ik (1−f_kj) / Σ_m q_im (1−f_mj).
    /// Iteration stops early once the log-likelihood gain (Equation 2) falls below
    /// <see cref="AncestryLogLikelihoodTolerance"/> (their Equation 5).
    /// </remarks>
    public static IEnumerable<AncestryProportion> EstimateAncestry(
        IEnumerable<(string IndividualId, IReadOnlyList<int> Genotypes)> individuals,
        IEnumerable<(string PopulationId, IReadOnlyList<double> AlleleFrequencies)> referencePops,
        int maxIterations = 100)
    {
        var indList = individuals.ToList();
        var refList = referencePops.ToList();

        if (indList.Count == 0 || refList.Count == 0)
            yield break;

        int k = refList.Count;
        int m = refList[0].AlleleFrequencies.Count;

        foreach (var (indId, genotypes) in indList)
        {
            if (genotypes.Count != m)
                continue;

            double[] proportions = EstimateIndividualAncestry(genotypes, refList, k, m, maxIterations);

            var propDict = new Dictionary<string, double>(k);
            for (int pop = 0; pop < k; pop++)
                propDict[refList[pop].PopulationId] = proportions[pop];

            yield return new AncestryProportion(indId, propDict);
        }
    }

    /// <summary>
    /// Runs the FRAPPE EM (Alexander et al. 2009, Eq. 4) for one individual with the reference
    /// allele frequencies held fixed, returning ancestry fractions that sum to 1.
    /// </summary>
    private static double[] EstimateIndividualAncestry(
        IReadOnlyList<int> genotypes,
        List<(string PopulationId, IReadOnlyList<double> AlleleFrequencies)> refList,
        int k,
        int m,
        int maxIterations)
    {
        // Initialize q uniformly over the K reference populations.
        var q = new double[k];
        for (int pop = 0; pop < k; pop++)
            q[pop] = 1.0 / k;

        double previousLogLik = AncestryLogLikelihood(genotypes, refList, q, k, m);

        for (int iter = 0; iter < maxIterations; iter++)
        {
            var next = new double[k];

            // EM update q_ik^{n+1} = (1/2J) Σ_j [ g_ij a_ijk + (2-g_ij) b_ijk ]  (Eq. 4),
            // where J counts only the informative (non-missing) SNPs for this individual.
            int informativeSnps = 0;
            for (int snp = 0; snp < m; snp++)
            {
                int g = genotypes[snp];
                if (g < 0 || g > 2)
                    continue; // missing genotype contributes no likelihood term (Eq. 2)
                informativeSnps++;

                // Admixed allele-1 / allele-2 frequencies for this individual: Σ_m q_im f_mj.
                double mixAllele1 = 0;
                double mixAllele2 = 0;
                for (int pop = 0; pop < k; pop++)
                {
                    double f = refList[pop].AlleleFrequencies[snp];
                    mixAllele1 += q[pop] * f;
                    mixAllele2 += q[pop] * (1 - f);
                }

                for (int pop = 0; pop < k; pop++)
                {
                    double f = refList[pop].AlleleFrequencies[snp];
                    double a = mixAllele1 > 0 ? q[pop] * f / mixAllele1 : 0;
                    double b = mixAllele2 > 0 ? q[pop] * (1 - f) / mixAllele2 : 0;
                    next[pop] += g * a + (2 - g) * b;
                }
            }

            if (informativeSnps == 0)
                break; // nothing informative: keep the uniform prior

            // Divide by 2J (J = informative SNP count). The EM preserves Σ_k q_ik = 1 exactly.
            double normalizer = 2.0 * informativeSnps;
            for (int pop = 0; pop < k; pop++)
                q[pop] = next[pop] / normalizer;

            double logLik = AncestryLogLikelihood(genotypes, refList, q, k, m);
            if (logLik - previousLogLik < AncestryLogLikelihoodTolerance)
                break; // converged per Eq. 5
            previousLogLik = logLik;
        }

        return q;
    }

    /// <summary>
    /// Binomial admixture log-likelihood for one individual given fixed allele frequencies
    /// (Alexander et al. 2009, Equation 2), summed over informative SNPs.
    /// </summary>
    private static double AncestryLogLikelihood(
        IReadOnlyList<int> genotypes,
        List<(string PopulationId, IReadOnlyList<double> AlleleFrequencies)> refList,
        double[] q,
        int k,
        int m)
    {
        double logLik = 0;
        for (int snp = 0; snp < m; snp++)
        {
            int g = genotypes[snp];
            if (g < 0 || g > 2)
                continue;

            double mixAllele1 = 0;
            double mixAllele2 = 0;
            for (int pop = 0; pop < k; pop++)
            {
                double f = refList[pop].AlleleFrequencies[snp];
                mixAllele1 += q[pop] * f;
                mixAllele2 += q[pop] * (1 - f);
            }

            if (g > 0 && mixAllele1 > 0)
                logLik += g * Math.Log(mixAllele1);
            if (g < 2 && mixAllele2 > 0)
                logLik += (2 - g) * Math.Log(mixAllele2);
        }
        return logLik;
    }

    #endregion

    #region Inbreeding

    // Default minimum number of homozygous SNPs a run must contain to be retained.
    // PLINK 1.9 --homozyg-snp default is 100 consecutive SNPs (Chang et al. 2015,
    // PLINK 1.9 "Runs of homozygosity" documentation, --homozyg-snp).
    private const int DefaultRohMinSnps = 100;

    // Default minimum physical length (in base pairs) for a retained run.
    // PLINK 1.9 --homozyg-kb default is 1000 kb = 1,000,000 bp (Chang et al. 2015).
    private const int DefaultRohMinLengthBp = 1_000_000;

    // Default maximum number of heterozygous (opposite) genotypes tolerated inside a
    // single run before it is broken. The consecutive-runs method of Marras et al.
    // (2015) allows a small number of opposite genotypes (maxOppRun) to account for
    // genotyping error; PLINK's scanning window likewise permits one heterozygous call
    // per window (--homozyg-window-het default 1).
    private const int DefaultRohMaxHeterozygotes = 1;

    // Default maximum physical gap (in base pairs) between two consecutive SNPs that may
    // still belong to the same run. PLINK 1.9 --homozyg-gap default is 1000 kb =
    // 1,000,000 bp (Chang et al. 2015): SNPs farther apart than this break the run.
    private const int DefaultRohMaxGapBp = 1_000_000;

    /// <summary>
    /// Calculates the genomic inbreeding coefficient F_ROH from runs of homozygosity:
    /// F_ROH = (Σ L_ROH) / L_AUTO, the total length of an individual's runs of homozygosity
    /// (above a chosen minimum length) divided by the length of the autosomal genome covered
    /// by SNPs (McQuillan et al. 2008, Eq. for F_roh = ΣL_roh / L_auto).
    /// </summary>
    /// <param name="rohSegments">
    /// Half-open ROH intervals <c>[Start, End)</c>; each contributes length <c>End − Start</c>.
    /// </param>
    /// <param name="genomeLength">
    /// L_AUTO — the SNP-covered autosomal genome length in the same units as the segments
    /// (base pairs). Must be positive.
    /// </param>
    /// <returns>
    /// F_ROH in [0, 1] when segments lie within the genome. Returns 0 when
    /// <paramref name="genomeLength"/> is non-positive (no defined denominator).
    /// </returns>
    public static double CalculateInbreedingFromROH(
        IEnumerable<(int Start, int End)> rohSegments,
        int genomeLength)
    {
        ArgumentNullException.ThrowIfNull(rohSegments);

        if (genomeLength <= 0)
            return 0;

        long totalROH = rohSegments.Sum(r => (long)(r.End - r.Start));
        return (double)totalROH / genomeLength;
    }

    /// <summary>
    /// Identifies runs of homozygosity (ROH) with the window-free consecutive-runs method
    /// of Marras et al. (2015): the genome is scanned SNP by SNP in ascending position
    /// order; a candidate run is extended while it contains at most
    /// <paramref name="maxHeterozygotes"/> opposite (heterozygous) genotypes and no gap
    /// between consecutive SNPs exceeds <paramref name="maxGap"/>. A run that violates either
    /// limit is closed (ending at the last SNP that still satisfied them), and a new run
    /// starts at the breaking SNP. A closed run is reported only when it contains at least
    /// <paramref name="minSnps"/> homozygous SNPs and spans at least <paramref name="minLength"/>
    /// base pairs (PLINK 1.9 --homozyg-snp / --homozyg-kb thresholds; Chang et al. 2015).
    /// </summary>
    /// <param name="genotypes">
    /// Per-SNP <c>(Position, Genotype)</c> pairs, where genotype 0 = homozygous reference and
    /// 2 = homozygous alternate (both homozygous) and 1 = heterozygous (opposite). Positions
    /// need not be pre-sorted; they are ordered ascending internally.
    /// </param>
    /// <param name="minSnps">Minimum number of SNPs in a retained run (default 100; PLINK --homozyg-snp).</param>
    /// <param name="minLength">Minimum physical length in bp of a retained run (default 1,000,000; PLINK --homozyg-kb 1000).</param>
    /// <param name="maxHeterozygotes">Maximum opposite (heterozygous) genotypes tolerated within a run (default 1).</param>
    /// <param name="maxGap">Maximum bp gap between consecutive SNPs in one run (default 1,000,000; PLINK --homozyg-gap 1000).</param>
    /// <returns>
    /// The retained runs as <c>(Start, End, SnpCount)</c>, where Start/End are the positions of
    /// the first and last SNP of the run (inclusive) and SnpCount is the number of SNPs in it,
    /// emitted in ascending Start order.
    /// </returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="genotypes"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="minSnps"/> &lt; 1, or any of <paramref name="minLength"/>,
    /// <paramref name="maxHeterozygotes"/>, <paramref name="maxGap"/> is negative.
    /// </exception>
    public static IEnumerable<(int Start, int End, int SnpCount)> FindROH(
        IEnumerable<(int Position, int Genotype)> genotypes,
        int minSnps = DefaultRohMinSnps,
        int minLength = DefaultRohMinLengthBp,
        int maxHeterozygotes = DefaultRohMaxHeterozygotes,
        int maxGap = DefaultRohMaxGapBp)
    {
        ArgumentNullException.ThrowIfNull(genotypes);
        if (minSnps < 1)
            throw new ArgumentOutOfRangeException(nameof(minSnps), minSnps, "Minimum SNP count must be at least 1.");
        if (minLength < 0)
            throw new ArgumentOutOfRangeException(nameof(minLength), minLength, "Minimum length cannot be negative.");
        if (maxHeterozygotes < 0)
            throw new ArgumentOutOfRangeException(nameof(maxHeterozygotes), maxHeterozygotes, "Maximum heterozygotes cannot be negative.");
        if (maxGap < 0)
            throw new ArgumentOutOfRangeException(nameof(maxGap), maxGap, "Maximum gap cannot be negative.");

        return FindROHIterator(genotypes, minSnps, minLength, maxHeterozygotes, maxGap);
    }

    private static IEnumerable<(int Start, int End, int SnpCount)> FindROHIterator(
        IEnumerable<(int Position, int Genotype)> genotypes,
        int minSnps,
        int minLength,
        int maxHeterozygotes,
        int maxGap)
    {
        var genoList = genotypes.OrderBy(g => g.Position).ToList();
        if (genoList.Count == 0)
            yield break;

        // State of the currently growing run. The run is reported on the closed interval
        // [runStartIndex .. lastHomIndex]: a run must begin and end on a homozygous SNP, so
        // trailing tolerated heterozygotes are not part of the emitted run. snpCountAtLastHom
        // is the SNP count of that emitted interval (interior tolerated heterozygotes included).
        int runStartIndex = 0;     // first SNP of the current run (homozygous)
        int snpCount = 0;          // SNPs scanned into the current run so far
        int hetCount = 0;          // opposite genotypes inside the current run so far
        int lastHomIndex = 0;      // index of the most recent homozygous SNP in the run
        int snpCountAtLastHom = 0; // snpCount as of lastHomIndex (the emitted interval size)

        bool QualifiesForEmit(out (int, int, int) run)
        {
            int start = genoList[runStartIndex].Position;
            int end = genoList[lastHomIndex].Position;
            run = (start, end, snpCountAtLastHom);
            return snpCountAtLastHom >= minSnps && end - start >= minLength;
        }

        for (int i = 0; i < genoList.Count; i++)
        {
            var (pos, geno) = genoList[i];
            bool isHet = geno == 1;

            // A gap larger than maxGap breaks the run before this SNP is considered.
            bool gapBreak = snpCount > 0 && pos - genoList[i - 1].Position > maxGap;
            // Adding this SNP would push opposite-genotype count past the tolerance.
            bool hetBreak = isHet && hetCount + 1 > maxHeterozygotes;

            if (snpCount > 0 && (gapBreak || hetBreak))
            {
                if (QualifiesForEmit(out var run))
                    yield return run;

                // Restart a fresh run at the breaking SNP. A heterozygous breaker cannot seed a
                // homozygous run, so begin counting from the next homozygous SNP instead.
                runStartIndex = i;
                snpCount = isHet ? 0 : 1;
                hetCount = 0;
                lastHomIndex = i;
                snpCountAtLastHom = isHet ? 0 : 1;
                continue;
            }

            if (snpCount == 0)
            {
                if (isHet)
                    continue; // skip leading heterozygous SNPs until a run can start
                runStartIndex = i;
            }

            snpCount++;
            if (isHet)
            {
                hetCount++;
            }
            else
            {
                lastHomIndex = i;
                snpCountAtLastHom = snpCount;
            }
        }

        if (snpCount > 0 && QualifiesForEmit(out var finalRun))
            yield return finalRun;
    }

    #endregion
}
