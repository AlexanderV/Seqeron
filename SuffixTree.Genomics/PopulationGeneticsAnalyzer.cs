using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics;

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
    public static (double MajorFreq, double MinorFreq) CalculateAlleleFrequencies(
        int homozygousMajor,
        int heterozygous,
        int homozygousMinor)
    {
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
    /// </summary>
    public static double CalculateTajimasD(
        double nucleotideDiversity,
        double wattersonTheta,
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

        // Calculate constants
        double b1 = (n + 1.0) / (3 * (n - 1));
        double b2 = 2.0 * (n * n + n + 3) / (9 * n * (n - 1));
        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);

        // Calculate variance
        double variance = e1 * segregatingSites + e2 * segregatingSites * (segregatingSites - 1);

        if (variance <= 0)
            return 0;

        // Tajima's D
        double d = nucleotideDiversity - wattersonTheta;
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
        double tajD = CalculateTajimasD(pi, theta, segregatingSites, n);

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

    private static double CalculateObservedHeterozygosity(List<IReadOnlyList<char>> sequences)
    {
        if (sequences.Count < 2)
            return 0;

        int length = sequences[0].Count;
        int hetSites = 0;

        for (int pos = 0; pos < length; pos++)
        {
            var alleles = sequences.Select(s => s[pos]).Distinct().ToList();
            if (alleles.Count > 1)
                hetSites++;
        }

        return (double)hetSites / length;
    }

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

        // Approximation using incomplete gamma function
        return LowerIncompleteGamma(df / 2.0, x / 2.0) / Gamma(df / 2.0);
    }

    private static double LowerIncompleteGamma(double a, double x)
    {
        if (x < 0)
            return 0;

        double sum = 0;
        double term = 1.0 / a;
        sum = term;

        for (int n = 1; n < 100; n++)
        {
            term *= x / (a + n);
            sum += term;
            if (Math.Abs(term) < 1e-10)
                break;
        }

        return Math.Pow(x, a) * Math.Exp(-x) * sum;
    }

    private static double Gamma(double x)
    {
        // Stirling approximation for gamma function
        if (x < 0.5)
            return Math.PI / (Math.Sin(Math.PI * x) * Gamma(1 - x));

        x -= 1;
        double[] g = { 1.0, 0.5772156649015329, -0.6558780715202538, -0.0420026350340952 };
        double result = g[0];

        for (int i = 1; i < g.Length; i++)
            result += g[i] * Math.Pow(x, i);

        return Math.Sqrt(2 * Math.PI / x) * Math.Pow(x / Math.E, x) * result;
    }

    #endregion

    #region Population Structure (F-statistics)

    /// <summary>
    /// Calculates Weir and Cockerham's Fst between populations.
    /// </summary>
    public static double CalculateFst(
        IEnumerable<(double AlleleFreq, int SampleSize)> population1,
        IEnumerable<(double AlleleFreq, int SampleSize)> population2)
    {
        var pop1 = population1.ToList();
        var pop2 = population2.ToList();

        if (pop1.Count == 0 || pop2.Count == 0)
            return 0;

        double numerator = 0;
        double denominator = 0;

        for (int i = 0; i < Math.Min(pop1.Count, pop2.Count); i++)
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

        // Calculate allele frequencies
        double p1 = genoList.Sum(g => g.Geno1) / (2.0 * genoList.Count);
        double p2 = genoList.Sum(g => g.Geno2) / (2.0 * genoList.Count);
        double q1 = 1 - p1;
        double q2 = 1 - p2;

        // Calculate haplotype frequencies (assuming random phase)
        double p11 = 0;
        foreach (var (geno1, geno2) in genoList)
        {
            // Estimate haplotype from genotypes
            if (geno1 == 2 && geno2 == 2) p11 += 1;
            else if (geno1 == 2 && geno2 == 1) p11 += 0.5;
            else if (geno1 == 1 && geno2 == 2) p11 += 0.5;
            else if (geno1 == 1 && geno2 == 1) p11 += 0.25;
        }
        p11 /= genoList.Count;

        // Calculate D
        double d = p11 - p1 * p2;

        // Calculate D'
        double dMax = d >= 0
            ? Math.Min(p1 * q2, q1 * p2)
            : Math.Min(p1 * p2, q1 * q2);

        double dPrime = dMax != 0 ? Math.Abs(d) / dMax : 0;

        // Calculate r²
        double denom = p1 * q1 * p2 * q2;
        double rSquared = denom > 0 ? (d * d) / denom : 0;

        return new LinkageDisequilibrium(
            Variant1: variant1Id,
            Variant2: variant2Id,
            DPrime: dPrime,
            RSquared: rSquared,
            Distance: distance);
    }

    /// <summary>
    /// Identifies haplotype blocks using the four-gamete test.
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
                End: variantList.Last().Position,
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

        return 2 * (1 - NormalCDF(z));
    }

    private static double NormalCDF(double x)
    {
        return 0.5 * (1 + Erf(x / Math.Sqrt(2)));
    }

    private static double Erf(double x)
    {
        double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
        double a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;

        int sign = x < 0 ? -1 : 1;
        x = Math.Abs(x);

        double t = 1.0 / (1.0 + p * x);
        double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);

        return sign * y;
    }

    #endregion

    #region Ancestry Analysis

    /// <summary>
    /// Estimates ancestry proportions using a simplified ADMIXTURE-like approach.
    /// </summary>
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

            // Initialize proportions uniformly
            var proportions = new double[k];
            for (int i = 0; i < k; i++)
                proportions[i] = 1.0 / k;

            // EM-like optimization (simplified)
            for (int iter = 0; iter < maxIterations; iter++)
            {
                var newProportions = new double[k];

                for (int snp = 0; snp < m; snp++)
                {
                    int geno = genotypes[snp];

                    // Calculate likelihood for each population
                    var likelihoods = new double[k];
                    double totalLik = 0;

                    for (int pop = 0; pop < k; pop++)
                    {
                        double p = refList[pop].AlleleFrequencies[snp];

                        // Genotype probability under HWE
                        double prob = geno switch
                        {
                            0 => (1 - p) * (1 - p),
                            1 => 2 * p * (1 - p),
                            2 => p * p,
                            _ => 0.25
                        };

                        likelihoods[pop] = proportions[pop] * prob;
                        totalLik += likelihoods[pop];
                    }

                    if (totalLik > 0)
                    {
                        for (int pop = 0; pop < k; pop++)
                        {
                            newProportions[pop] += likelihoods[pop] / totalLik;
                        }
                    }
                }

                // Normalize
                double sum = newProportions.Sum();
                for (int pop = 0; pop < k; pop++)
                {
                    proportions[pop] = sum > 0 ? newProportions[pop] / sum : 1.0 / k;
                }
            }

            // Create result
            var propDict = new Dictionary<string, double>();
            for (int pop = 0; pop < k; pop++)
            {
                propDict[refList[pop].PopulationId] = proportions[pop];
            }

            yield return new AncestryProportion(indId, propDict);
        }
    }

    #endregion

    #region Inbreeding

    /// <summary>
    /// Calculates inbreeding coefficient from runs of homozygosity.
    /// </summary>
    public static double CalculateInbreedingFromROH(
        IEnumerable<(int Start, int End)> rohSegments,
        int genomeLength)
    {
        if (genomeLength <= 0)
            return 0;

        long totalROH = rohSegments.Sum(r => (long)(r.End - r.Start));
        return (double)totalROH / genomeLength;
    }

    /// <summary>
    /// Identifies runs of homozygosity (ROH).
    /// </summary>
    public static IEnumerable<(int Start, int End, int SnpCount)> FindROH(
        IEnumerable<(int Position, int Genotype)> genotypes,
        int minSnps = 50,
        int minLength = 1_000_000,
        int maxHeterozygotes = 1)
    {
        var genoList = genotypes.OrderBy(g => g.Position).ToList();

        if (genoList.Count < minSnps)
            yield break;

        int segmentStart = genoList[0].Position;
        int snpCount = 0;
        int hetCount = 0;

        for (int i = 0; i < genoList.Count; i++)
        {
            var (pos, geno) = genoList[i];

            if (geno == 1) // Heterozygous
            {
                hetCount++;
            }

            snpCount++;

            if (hetCount > maxHeterozygotes)
            {
                // End segment
                if (snpCount >= minSnps)
                {
                    int segmentEnd = genoList[i - 1].Position;
                    if (segmentEnd - segmentStart >= minLength)
                    {
                        yield return (segmentStart, segmentEnd, snpCount - 1);
                    }
                }

                // Start new segment
                segmentStart = pos;
                snpCount = 1;
                hetCount = geno == 1 ? 1 : 0;
            }
        }

        // Check final segment
        if (snpCount >= minSnps)
        {
            int segmentEnd = genoList.Last().Position;
            if (segmentEnd - segmentStart >= minLength)
            {
                yield return (segmentStart, segmentEnd, snpCount);
            }
        }
    }

    #endregion
}
