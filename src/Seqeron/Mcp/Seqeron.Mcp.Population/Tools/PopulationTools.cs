using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Population;

namespace Seqeron.Mcp.Population.Tools;

/// <summary>
/// MCP tools for population genetics operations.
/// </summary>
[McpServerToolType]
public class PopulationTools
{
    #region PopulationGeneticsAnalyzer

    // -------- Allele Frequency --------

    [McpServerTool(Name = "allele_frequencies", Title = "Population — Allele Frequencies", ReadOnly = true)]
    [Description("Calculate major and minor allele frequencies from diploid genotype counts (homozygous-major, heterozygous, homozygous-minor).")]
    public static AlleleFrequenciesResult AlleleFrequencies(
        [Description("Count of homozygous-major genotypes (AA).")] int homozygousMajor,
        [Description("Count of heterozygous genotypes (Aa).")] int heterozygous,
        [Description("Count of homozygous-minor genotypes (aa).")] int homozygousMinor)
    {
        var (major, minor) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor, heterozygous, homozygousMinor);
        return new AlleleFrequenciesResult(major, minor);
    }

    [McpServerTool(Name = "minor_allele_frequency", Title = "Population — Minor Allele Frequency", ReadOnly = true)]
    [Description("Compute minor allele frequency (MAF) from a vector of diploid genotypes encoded as 0 (hom-ref), 1 (het), or 2 (hom-alt).")]
    public static MinorAlleleFrequencyResult MinorAlleleFrequency(
        [Description("Genotype vector with values in {0,1,2}.")] int[] genotypes)
    {
        var maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);
        return new MinorAlleleFrequencyResult(maf);
    }

    [McpServerTool(Name = "filter_variants_by_maf", Title = "Population — Filter Variants by MAF", ReadOnly = true)]
    [Description("Filter variants whose minor allele frequency lies within [minMAF, maxMAF].")]
    public static FilterVariantsByMafResult FilterVariantsByMaf(
        [Description("Variants to filter.")] VariantItem[] variants,
        [Description("Inclusive lower bound on MAF (default 0.01).")] double minMAF = 0.01,
        [Description("Inclusive upper bound on MAF (default 0.5).")] double maxMAF = 0.5)
    {
        var source = variants.Select(v => new PopulationGeneticsAnalyzer.Variant(
            v.Id, v.Chromosome, v.Position, v.ReferenceAllele, v.AlternateAllele,
            v.AlleleFrequency, v.SampleCount));

        var filtered = PopulationGeneticsAnalyzer
            .FilterByMAF(source, minMAF, maxMAF)
            .Select(v => new VariantItem(
                v.Id, v.Chromosome, v.Position, v.ReferenceAllele, v.AlternateAllele,
                v.AlleleFrequency, v.SampleCount))
            .ToList();

        return new FilterVariantsByMafResult(filtered);
    }

    // -------- Diversity Statistics --------

    [McpServerTool(Name = "nucleotide_diversity", Title = "Population — Nucleotide Diversity (π)", ReadOnly = true)]
    [Description("Compute nucleotide diversity π — average pairwise per-site differences across an aligned set of equal-length sequences.")]
    public static NucleotideDiversityResult NucleotideDiversity(
        [Description("Aligned sequences of equal length.")] string[] sequences)
    {
        var seqs = sequences.Select(s => (IReadOnlyList<char>)s.ToCharArray()).ToList();
        var pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(seqs);
        return new NucleotideDiversityResult(pi);
    }

    [McpServerTool(Name = "wattersons_theta", Title = "Population — Watterson's θ", ReadOnly = true)]
    [Description("Compute Watterson's θ estimator from segregating site count, sample size, and sequence length.")]
    public static WattersonsThetaResult WattersonsTheta(
        [Description("Number of segregating (polymorphic) sites S.")] int segregatingSites,
        [Description("Sample size n (number of sequences).")] int sampleSize,
        [Description("Length of analysed sequence in nucleotides.")] int sequenceLength)
    {
        var theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(
            segregatingSites, sampleSize, sequenceLength);
        return new WattersonsThetaResult(theta);
    }

    [McpServerTool(Name = "tajimas_d", Title = "Population — Tajima's D", ReadOnly = true)]
    [Description("Compute Tajima's D (Tajima 1989) from average pairwise differences k̂, segregating sites S, and sample size n (n ≥ 3 required).")]
    public static TajimasDResult TajimasD(
        [Description("k̂ — average number of pairwise differences between sequences (NOT per-site).")] double averagePairwiseDifferences,
        [Description("Number of segregating sites S.")] int segregatingSites,
        [Description("Sample size n.")] int sampleSize)
    {
        var d = PopulationGeneticsAnalyzer.CalculateTajimasD(
            averagePairwiseDifferences, segregatingSites, sampleSize);
        return new TajimasDResult(d);
    }

    [McpServerTool(Name = "diversity_statistics", Title = "Population — Diversity Statistics", ReadOnly = true)]
    [Description("Compute combined diversity statistics (π, Watterson's θ, Tajima's D, segregating sites, observed/expected heterozygosity) from aligned sequences.")]
    public static DiversityStatisticsResult DiversityStatistics(
        [Description("Aligned sequences of equal length.")] string[] sequences)
    {
        var seqs = sequences.Select(s => (IReadOnlyList<char>)s.ToCharArray()).ToList();
        var s = PopulationGeneticsAnalyzer.CalculateDiversityStatistics(seqs);
        return new DiversityStatisticsResult(
            NucleotideDiversity: s.NucleotideDiversity,
            WattersonTheta: s.WattersonTheta,
            TajimasD: s.TajimasD,
            SegregatingSites: s.SegregratingSites,
            SampleSize: s.SampleSize,
            HeterozygosityObserved: s.HeterozygosityObserved,
            HeterozygosityExpected: s.HeterozygosityExpected);
    }

    // -------- Hardy-Weinberg --------

    [McpServerTool(Name = "hardy_weinberg_test", Title = "Population — Hardy-Weinberg Test", ReadOnly = true)]
    [Description("Test Hardy-Weinberg equilibrium for a biallelic variant via chi-square goodness-of-fit (1 df) on observed AA/Aa/aa counts.")]
    public static HardyWeinbergTestResult HardyWeinbergTest(
        [Description("Variant identifier.")] string variantId,
        [Description("Observed homozygous-AA count.")] int observedAA,
        [Description("Observed heterozygous-Aa count.")] int observedAa,
        [Description("Observed homozygous-aa count.")] int observedaa,
        [Description("Significance level α (default 0.05).")] double significanceLevel = 0.05)
    {
        var r = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            variantId, observedAA, observedAa, observedaa, significanceLevel);
        return new HardyWeinbergTestResult(
            r.VariantId, r.ObservedAA, r.ObservedAa, r.Observedaa,
            r.ExpectedAA, r.ExpectedAa, r.Expectedaa,
            r.ChiSquare, r.PValue, r.InEquilibrium);
    }

    // -------- F-statistics --------

    [McpServerTool(Name = "fst", Title = "Population — Wright's Fst", ReadOnly = true)]
    [Description("Compute Wright's variance-based Fst between two populations from per-variant allele frequencies and sample sizes (Wright 1965).")]
    public static FstResult Fst(
        [Description("Per-variant {alleleFreq, sampleSize} for population 1.")] AlleleItem[] population1,
        [Description("Per-variant {alleleFreq, sampleSize} for population 2.")] AlleleItem[] population2)
    {
        var p1 = population1.Select(a => (a.AlleleFreq, a.SampleSize));
        var p2 = population2.Select(a => (a.AlleleFreq, a.SampleSize));
        var fst = PopulationGeneticsAnalyzer.CalculateFst(p1, p2);
        return new FstResult(fst);
    }

    [McpServerTool(Name = "pairwise_fst", Title = "Population — Pairwise Fst Matrix", ReadOnly = true)]
    [Description("Compute pairwise Fst for a set of populations; returns the population id list and the symmetric Fst matrix.")]
    public static PairwiseFstResult PairwiseFst(
        [Description("Populations with their per-variant allele frequencies.")] PopulationItem[] populations)
    {
        var input = populations.Select(p =>
            (p.PopulationId,
             (IReadOnlyList<(double AlleleFreq, int SampleSize)>)p.Variants
                 .Select(v => (v.AlleleFreq, v.SampleSize))
                 .ToList()));

        var matrix = PopulationGeneticsAnalyzer.CalculatePairwiseFst(input);

        int n = populations.Length;
        var rows = new List<IReadOnlyList<double>>(n);
        for (int i = 0; i < n; i++)
        {
            var row = new double[n];
            for (int j = 0; j < n; j++)
                row[j] = matrix[i, j];
            rows.Add(row);
        }

        var ids = populations.Select(p => p.PopulationId).ToList();
        return new PairwiseFstResult(ids, rows);
    }

    [McpServerTool(Name = "f_statistics", Title = "Population — F-statistics (Fis, Fit, Fst)", ReadOnly = true)]
    [Description("Compute Wright's F-statistics (Fis, Fit, Fst) between two populations from per-variant heterozygosities and allele frequencies.")]
    public static FStatisticsResult FStatistics(
        [Description("Population 1 name.")] string pop1Name,
        [Description("Population 2 name.")] string pop2Name,
        [Description("Per-variant data: observed heterozygote counts and sample sizes for each population, plus allele frequencies.")] VariantDataItem[] variantData)
    {
        var data = variantData.Select(v =>
            (v.HetObs1, v.N1, v.HetObs2, v.N2, v.AlleleFreq1, v.AlleleFreq2));
        var r = PopulationGeneticsAnalyzer.CalculateFStatistics(pop1Name, pop2Name, data);
        return new FStatisticsResult(r.Fst, r.Fis, r.Fit, r.Population1, r.Population2);
    }

    // -------- Linkage Disequilibrium --------

    [McpServerTool(Name = "linkage_disequilibrium", Title = "Population — Linkage Disequilibrium", ReadOnly = true)]
    [Description("Compute pairwise linkage disequilibrium (D' and r²) between two variants from diploid genotype pairs (0/1/2 encoding). Uses Hill & Robertson (1968) and Lewontin (1964) normalizations.")]
    public static LinkageDisequilibriumResult LinkageDisequilibrium(
        [Description("Variant 1 identifier.")] string variant1Id,
        [Description("Variant 2 identifier.")] string variant2Id,
        [Description("Per-individual genotype pairs at the two variants (0/1/2).")] GenotypePairItem[] genotypes,
        [Description("Genomic distance between variants (e.g., bp).")] int distance)
    {
        var pairs = genotypes.Select(g => (g.Geno1, g.Geno2));
        var ld = PopulationGeneticsAnalyzer.CalculateLD(variant1Id, variant2Id, pairs, distance);
        return new LinkageDisequilibriumResult(
            ld.Variant1, ld.Variant2, ld.DPrime, ld.RSquared, ld.Distance);
    }

    [McpServerTool(Name = "haplotype_blocks", Title = "Population — Haplotype Blocks", ReadOnly = true)]
    [Description("Detect haplotype blocks via adjacent-pair r² ≥ ldThreshold (simplified Gabriel et al. 2002).")]
    public static HaplotypeBlocksResult HaplotypeBlocks(
        [Description("Variants with positions and per-individual genotypes (0/1/2).")] VariantGenotypesItem[] variants,
        [Description("Minimum r² to extend a block (default 0.7).")] double ldThreshold = 0.7)
    {
        var input = variants.Select(v => (v.VariantId, v.Position, (IReadOnlyList<int>)v.Genotypes.ToArray()));
        var blocks = PopulationGeneticsAnalyzer
            .FindHaplotypeBlocks(input, ldThreshold)
            .Select(b => new HaplotypeBlockItem(
                b.Start,
                b.End,
                b.Variants.ToList(),
                b.Haplotypes.Select(h => new HaplotypeFrequencyItem(h.Haplotype, h.Frequency)).ToList()))
            .ToList();
        return new HaplotypeBlocksResult(blocks);
    }

    // -------- Selection --------

    [McpServerTool(Name = "integrated_haplotype_score", Title = "Population — Integrated Haplotype Score (iHS)", ReadOnly = true)]
    [Description("Compute integrated haplotype score iHS = ln(iHH₁ / iHH₀) by trapezoidal integration of EHH curves over genomic positions.")]
    public static IhsResult IntegratedHaplotypeScore(
        [Description("EHH values for the ancestral (0) allele, aligned with positions.")] double[] ehh0,
        [Description("EHH values for the derived (1) allele, aligned with positions.")] double[] ehh1,
        [Description("Genomic positions corresponding to ehh0/ehh1 entries.")] int[] positions)
    {
        var ihs = PopulationGeneticsAnalyzer.CalculateIHS(ehh0, ehh1, positions);
        return new IhsResult(ihs);
    }

    [McpServerTool(Name = "scan_selection_signals", Title = "Population — Scan Selection Signals", ReadOnly = true)]
    [Description("Scan genomic regions for selection signals using Tajima's D, Fst, and iHS thresholds. Emits one signal per crossing test per region.")]
    public static ScanSelectionSignalsResult ScanSelectionSignals(
        [Description("Regions with precomputed Tajima's D, Fst, and iHS scores.")] RegionItem[] regions,
        [Description("Tajima's D threshold for negative-D selection signal (default -2.0).")] double tajimaDThreshold = -2.0,
        [Description("Fst threshold for differentiation signal (default 0.25).")] double fstThreshold = 0.25,
        [Description("|iHS| threshold for haplotype-based selection signal (default 2.0).")] double ihsThreshold = 2.0)
    {
        var input = regions.Select(r => (r.Region, r.Start, r.End, r.TajimaD, r.Fst, r.IHS));
        var signals = PopulationGeneticsAnalyzer
            .ScanForSelection(input, tajimaDThreshold, fstThreshold, ihsThreshold)
            .Select(s => new SelectionSignalItem(
                s.Region, s.Start, s.End, s.Score, s.TestType, s.PValue, s.Interpretation))
            .ToList();
        return new ScanSelectionSignalsResult(signals);
    }

    [McpServerTool(Name = "estimate_ancestry", Title = "Population — Estimate Ancestry", ReadOnly = true)]
    [Description("Estimate per-individual ancestry proportions across reference populations via a simplified ADMIXTURE-like EM procedure.")]
    public static EstimateAncestryResult EstimateAncestry(
        [Description("Individuals with per-SNP genotypes (0/1/2).")] IndividualItem[] individuals,
        [Description("Reference populations with per-SNP allele frequencies (alt allele).")] RefPopItem[] referencePops,
        [Description("Maximum EM iterations (default 100).")] int maxIterations = 100)
    {
        var ind = individuals.Select(i => (i.IndividualId, (IReadOnlyList<int>)i.Genotypes.ToArray()));
        var refs = referencePops.Select(r => (r.PopulationId, (IReadOnlyList<double>)r.AlleleFrequencies.ToArray()));
        var items = PopulationGeneticsAnalyzer
            .EstimateAncestry(ind, refs, maxIterations)
            .Select(a => new AncestryProportionItem(
                a.IndividualId,
                a.Proportions.ToDictionary(kv => kv.Key, kv => kv.Value)))
            .ToList();
        return new EstimateAncestryResult(items);
    }

    // -------- Inbreeding / ROH --------

    [McpServerTool(Name = "inbreeding_from_roh", Title = "Population — Inbreeding from ROH", ReadOnly = true)]
    [Description("Estimate the genomic inbreeding coefficient F_ROH as Σ(ROH lengths) / genomeLength.")]
    public static InbreedingFromRohResult InbreedingFromRoh(
        [Description("Runs of homozygosity as {start, end} segments.")] RohSegmentItem[] rohSegments,
        [Description("Total assayed genome length.")] int genomeLength)
    {
        var segs = rohSegments.Select(s => (s.Start, s.End));
        var f = PopulationGeneticsAnalyzer.CalculateInbreedingFromROH(segs, genomeLength);
        return new InbreedingFromRohResult(f);
    }

    [McpServerTool(Name = "runs_of_homozygosity", Title = "Population — Runs of Homozygosity", ReadOnly = true)]
    [Description("Identify runs of homozygosity (ROH) from per-SNP genotype calls (0/1/2) using minimum SNP count, minimum length, and maximum heterozygote tolerance.")]
    public static RunsOfHomozygosityResult RunsOfHomozygosity(
        [Description("Per-SNP {position, genotype} entries (genotype: 0=hom-ref, 1=het, 2=hom-alt).")] GenotypePositionItem[] genotypes,
        [Description("Minimum number of SNPs in a run (default 50).")] int minSnps = 50,
        [Description("Minimum run length in bp (default 1,000,000).")] int minLength = 1_000_000,
        [Description("Maximum tolerated heterozygous calls per run (default 1).")] int maxHeterozygotes = 1)
    {
        var input = genotypes.Select(g => (g.Position, g.Genotype));
        var items = PopulationGeneticsAnalyzer
            .FindROH(input, minSnps, minLength, maxHeterozygotes)
            .Select(r => new RohRegionItem(r.Start, r.End, r.SnpCount))
            .ToList();
        return new RunsOfHomozygosityResult(items);
    }

    #endregion
}
