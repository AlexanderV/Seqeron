using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Mcp.Metagenomics.Tools;

/// <summary>
/// MCP tools for metagenomics analysis operations.
/// </summary>
[McpServerToolType]
public class MetagenomicsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private MetagenomicsTools() { }

    #region MetagenomicsAnalyzer

    /// <summary>
    /// Classify metagenomic reads with the Kraken k-mer / LCA algorithm: collect each read's
    /// canonical k-mer hits, build the classification tree over the hit taxa weighted by k-mer
    /// count, and assign the leaf of the maximum-scoring root-to-leaf path (LCA-of-leaves on ties).
    /// Reads with no hits are assigned the root / unclassified taxon (Wood &amp; Salzberg 2014).
    /// </summary>
    [McpServerTool(Name = "classify_reads", Title = "Metagenomics — Classify Reads", ReadOnly = true)]
    [Description("Classify metagenomic reads with the Kraken k-mer/LCA algorithm against a canonical-k-mer→taxon-id database and a taxonomy tree. Returns the assigned taxon, RTL score, lineage and C/Q confidence per read.")]
    public static ClassifyReadsResult ClassifyReads(
        [Description("Reads to classify (id + nucleotide sequence).")] IReadOnlyList<ReadInput> reads,
        [Description("Flattened canonical-k-mer→taxon-id database. Keys must be canonical k-mers (lex-min of forward / reverse complement); values are taxon ids in the taxonomy tree.")] IReadOnlyList<KmerDatabaseEntry> kmerDatabase,
        [Description("Taxonomy tree nodes (id, name, rank, parent id); the root is self-parented.")] IReadOnlyList<TaxonNodeInput> taxonomy,
        [Description("k-mer length (default 31, per Kraken).")] int k = 31)
    {
        var dict = new Dictionary<string, int>(kmerDatabase.Count);
        foreach (var entry in kmerDatabase)
            dict[entry.Kmer] = entry.TaxonId;

        var tree = ToTaxonomyTree(taxonomy);
        var inputs = reads.Select(r => (r.Id, r.Sequence));
        var classifications = MetagenomicsAnalyzer.ClassifyReads(inputs, dict, tree, k).ToList();
        return new ClassifyReadsResult(classifications);
    }

    /// <summary>
    /// Build a Kraken canonical-k-mer→taxon-id database from labeled reference sequences:
    /// each shared k-mer is mapped to the lowest common ancestor of its owning taxa.
    /// </summary>
    [McpServerTool(Name = "build_kmer_database", Title = "Metagenomics — Build K-mer Database", ReadOnly = true)]
    [Description("Build a Kraken canonical-k-mer→taxon-id database from labeled reference sequences. Skips ambiguous nucleotides; a k-mer shared by several taxa maps to their lowest common ancestor in the taxonomy tree.")]
    public static BuildKmerDatabaseResult BuildKmerDatabase(
        [Description("Reference sequences (taxon id + nucleotide sequence).")] IReadOnlyList<ReferenceGenomeInput> referenceGenomes,
        [Description("Taxonomy tree nodes (id, name, rank, parent id); the root is self-parented.")] IReadOnlyList<TaxonNodeInput> taxonomy,
        [Description("k-mer length (default 31).")] int k = 31)
    {
        var tree = ToTaxonomyTree(taxonomy);
        var inputs = referenceGenomes.Select(g => (g.TaxonId, g.Sequence));
        var dict = MetagenomicsAnalyzer.BuildKmerDatabase(inputs, tree, k);
        var entries = dict.Select(kv => new KmerDatabaseEntry(kv.Key, kv.Value)).ToList();
        return new BuildKmerDatabaseResult(entries, entries.Count);
    }

    /// <summary>
    /// Aggregate per-read classifications into a sample-level taxonomic profile
    /// (per-rank abundances + Shannon / Simpson diversity at species level).
    /// </summary>
    [McpServerTool(Name = "taxonomic_profile", Title = "Metagenomics — Taxonomic Profile", ReadOnly = true)]
    [Description("Aggregate per-read classifications into rank-wise abundances plus Shannon/Simpson diversity (computed at species level).")]
    public static TaxonomicProfileResult TaxonomicProfile(
        [Description("Per-read classifications produced by classify_reads.")]
        IReadOnlyList<MetagenomicsAnalyzer.TaxonomicClassification> classifications)
    {
        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);
        return new TaxonomicProfileResult(
            KingdomAbundance: ToAbundanceList(profile.KingdomAbundance),
            PhylumAbundance: ToAbundanceList(profile.PhylumAbundance),
            GenusAbundance: ToAbundanceList(profile.GenusAbundance),
            SpeciesAbundance: ToAbundanceList(profile.SpeciesAbundance),
            ShannonDiversity: profile.ShannonDiversity,
            SimpsonDiversity: profile.SimpsonDiversity,
            TotalReads: profile.TotalReads,
            ClassifiedReads: profile.ClassifiedReads);
    }

    /// <summary>
    /// Compute alpha-diversity metrics (Shannon, Simpson, inverse Simpson, Chao1,
    /// observed species, Pielou's evenness) for a single sample.
    /// </summary>
    [McpServerTool(Name = "alpha_diversity", Title = "Metagenomics — Alpha Diversity", ReadOnly = true)]
    [Description("Compute Shannon, Simpson, inverse Simpson, Chao1 (Chao 1984), observed species, and Pielou's evenness for a sample's abundance vector.")]
    public static MetagenomicsAnalyzer.AlphaDiversity AlphaDiversity(
        [Description("Per-species abundances. For Chao1, integer count data is required; otherwise Chao1 falls back to S_obs.")]
        IReadOnlyList<AbundanceItem> abundances)
    {
        var dict = ToAbundanceDict(abundances);
        return MetagenomicsAnalyzer.CalculateAlphaDiversity(dict);
    }

    /// <summary>
    /// Compute beta-diversity metrics (Bray–Curtis, Jaccard) between two samples.
    /// </summary>
    [McpServerTool(Name = "beta_diversity", Title = "Metagenomics — Beta Diversity", ReadOnly = true)]
    [Description("Compute Bray–Curtis and Jaccard distances plus shared/unique species counts between two samples. UniFrac is reported as 0 (no phylogenetic tree input).")]
    public static MetagenomicsAnalyzer.BetaDiversity BetaDiversity(
        [Description("Name/identifier of sample 1.")] string sample1Name,
        [Description("Sample 1 abundance vector.")] IReadOnlyList<AbundanceItem> sample1,
        [Description("Name/identifier of sample 2.")] string sample2Name,
        [Description("Sample 2 abundance vector.")] IReadOnlyList<AbundanceItem> sample2)
    {
        return MetagenomicsAnalyzer.CalculateBetaDiversity(
            sample1Name,
            ToAbundanceDict(sample1),
            sample2Name,
            ToAbundanceDict(sample2));
    }

    /// <summary>
    /// Bin contigs into MAGs using k-means on GC, normalized coverage, and
    /// tetranucleotide frequency (Pearson distance, Teeling 2004 TETRA);
    /// reports completeness and GC-variance contamination (Parks 2014).
    /// </summary>
    [McpServerTool(Name = "bin_contigs", Title = "Metagenomics — Bin Contigs (MAG)", ReadOnly = true)]
    [Description("Cluster contigs into metagenome-assembled genomes using k-means over GC content, normalized coverage and tetranucleotide frequency (TETRA, Teeling 2004). Reports completeness and GC-variance contamination (Parks 2014).")]
    public static BinContigsResult BinContigs(
        [Description("Contigs with sequences and per-contig coverage.")] IReadOnlyList<ContigInput> contigs,
        [Description("Maximum number of bins / k for k-means (default 10).")] int numBins = 10,
        [Description("Minimum total bin length to be reported, in bp (default 500000).")] double minBinSize = 500000,
        [Description("Expected genome size in bp for completeness estimation (default 4_000_000).")] double expectedGenomeSize = 4_000_000)
    {
        var inputs = contigs.Select(c => (c.ContigId, c.Sequence, c.Coverage));
        var bins = MetagenomicsAnalyzer.BinContigs(inputs, numBins, minBinSize, expectedGenomeSize).ToList();
        return new BinContigsResult(bins);
    }

    /// <summary>
    /// Predict functional annotations for proteins via motif containment against
    /// a function database (KO / pathway / COG category).
    /// </summary>
    [McpServerTool(Name = "predict_functions", Title = "Metagenomics — Predict Functions", ReadOnly = true)]
    [Description("Predict functional annotations for proteins by motif containment against a function database, returning function/pathway/KO and inferred COG category.")]
    public static PredictFunctionsResult PredictFunctions(
        [Description("Predicted proteins (gene id + amino-acid sequence).")] IReadOnlyList<ProteinInput> proteins,
        [Description("Function database entries (motif key + function/pathway/KO).")] IReadOnlyList<FunctionDatabaseEntry> functionDatabase)
    {
        var dict = new Dictionary<string, (string Function, string Pathway, string Ko)>(functionDatabase.Count);
        foreach (var e in functionDatabase)
            dict[e.Motif] = (e.Function, e.Pathway, e.Ko);

        var inputs = proteins.Select(p => (p.GeneId, p.ProteinSequence));
        var annotations = MetagenomicsAnalyzer.PredictFunctions(inputs, dict).ToList();
        return new PredictFunctionsResult(annotations);
    }

    /// <summary>
    /// Compute functional richness / Shannon functional diversity over a set of
    /// annotations and the per-pathway hit counts.
    /// </summary>
    [McpServerTool(Name = "functional_diversity", Title = "Metagenomics — Functional Diversity", ReadOnly = true)]
    [Description("Compute functional richness, Shannon functional diversity, and per-pathway hit counts from a set of functional annotations.")]
    public static FunctionalDiversityResult FunctionalDiversity(
        [Description("Functional annotations (e.g. from predict_functions).")]
        IReadOnlyList<MetagenomicsAnalyzer.FunctionalAnnotation> annotations)
    {
        var (richness, diversity, pathwayCounts) =
            MetagenomicsAnalyzer.CalculateFunctionalDiversity(annotations);

        var counts = pathwayCounts
            .Select(kv => new PathwayCountItem(kv.Key, kv.Value))
            .ToList();

        return new FunctionalDiversityResult(richness, diversity, counts);
    }

    /// <summary>
    /// Search genes for antibiotic-resistance markers via motif containment
    /// against a resistance database.
    /// </summary>
    [McpServerTool(Name = "find_resistance_genes", Title = "Metagenomics — Resistance Gene Search", ReadOnly = true)]
    [Description("Search nucleotide genes for antibiotic-resistance markers via motif containment against a resistance database; reports the gene id, hit name, antibiotic class, and identity.")]
    public static FindResistanceGenesResult FindResistanceGenes(
        [Description("Genes to scan (gene id + nucleotide sequence).")] IReadOnlyList<GeneInput> genes,
        [Description("Resistance database entries (motif + name + antibiotic class).")] IReadOnlyList<ResistanceDatabaseEntry> resistanceDatabase)
    {
        var dict = new Dictionary<string, (string Name, string AntibioticClass)>(resistanceDatabase.Count);
        foreach (var e in resistanceDatabase)
            dict[e.Motif] = (e.Name, e.AntibioticClass);

        var inputs = genes.Select(g => (g.GeneId, g.Sequence));
        var hits = MetagenomicsAnalyzer.FindResistanceGenes(inputs, dict)
            .Select(t => new ResistanceGeneItem(t.GeneId, t.ResistanceGene, t.AntibioticClass, t.Identity))
            .ToList();
        return new FindResistanceGenesResult(hits);
    }

    /// <summary>
    /// Test for differential taxon abundance between two condition groups using
    /// Welch's t-test on per-sample abundances; flags hits with |log2FC| &gt; 1
    /// and p &lt; threshold.
    /// </summary>
    [McpServerTool(Name = "differential_abundance", Title = "Metagenomics — Differential Abundance", ReadOnly = true)]
    [Description("Welch's t-test for differential taxon abundance between two condition groups. Reports log2 fold-change, p-value, and a significance flag (p < threshold AND |log2FC| > 1).")]
    public static DifferentialAbundanceResult DifferentialAbundance(
        [Description("Per-sample abundance vectors for condition 1.")] IReadOnlyList<AbundanceSample> condition1Samples,
        [Description("Per-sample abundance vectors for condition 2.")] IReadOnlyList<AbundanceSample> condition2Samples,
        [Description("p-value significance threshold (default 0.05).")] double pValueThreshold = 0.05)
    {
        var c1 = condition1Samples.Select(s => (IReadOnlyDictionary<string, double>)ToAbundanceDict(s.Items));
        var c2 = condition2Samples.Select(s => (IReadOnlyDictionary<string, double>)ToAbundanceDict(s.Items));

        var items = MetagenomicsAnalyzer.DifferentialAbundance(c1, c2, pValueThreshold)
            .Select(t => new DifferentialAbundanceItem(t.Taxon, t.FoldChange, t.PValue, t.Significant))
            .ToList();
        return new DifferentialAbundanceResult(items);
    }

    #endregion

    #region PanGenomeAnalyzer

    /// <summary>
    /// Construct a pan-genome (core / accessory / unique partition + statistics
    /// + open vs closed classification, Tettelin et al. 2005) from a set of genomes.
    /// </summary>
    [McpServerTool(Name = "construct_pangenome", Title = "Pan-genome — Construct", ReadOnly = true)]
    [Description("Construct a pan-genome (core/accessory/unique partition, genome fluidity, open-vs-closed classification per Tettelin 2005) from a set of genomes.")]
    public static PanGenomeResultDto ConstructPanGenome(
        [Description("Genomes (id + ordered list of genes).")] IReadOnlyList<GenomeInput> genomes,
        [Description("Sequence-identity threshold for ortholog clustering (default 0.9).")] double identityThreshold = 0.9,
        [Description("Fraction of genomes a cluster must appear in to be 'core' (default 0.99).")] double coreFraction = 0.99)
    {
        var dict = ToGenomesDict(genomes);
        var result = PanGenomeAnalyzer.ConstructPanGenome(dict, identityThreshold, coreFraction);

        var genomeToGenes = result.GenomeToGenes
            .Select(kv => new GenomeToGenesItem(kv.Key, kv.Value.ToList()))
            .ToList();

        var statsDto = new PanGenomeStatisticsDto(
            result.Statistics.TotalGenomes,
            result.Statistics.TotalGenes,
            result.Statistics.CoreGeneCount,
            result.Statistics.AccessoryGeneCount,
            result.Statistics.UniqueGeneCount,
            result.Statistics.CoreFraction,
            result.Statistics.GenomeFluidity,
            result.Statistics.Type.ToString());

        return new PanGenomeResultDto(
            result.CoreGenes.ToList(),
            result.AccessoryGenes.ToList(),
            result.UniqueGenes.ToList(),
            genomeToGenes,
            statsDto);
    }

    /// <summary>
    /// Cluster genes from multiple genomes into ortholog groups using a 7-mer
    /// Jaccard similarity threshold.
    /// </summary>
    [McpServerTool(Name = "cluster_genes", Title = "Pan-genome — Cluster Genes", ReadOnly = true)]
    [Description("Cluster genes from multiple genomes into ortholog groups using 7-mer Jaccard similarity. Returns clusters with member gene/genome ids, average identity, and a representative consensus sequence.")]
    public static ClusterGenesResult ClusterGenes(
        [Description("Genomes (id + ordered list of genes).")] IReadOnlyList<GenomeInput> genomes,
        [Description("Sequence-identity threshold for cluster membership (default 0.9).")] double identityThreshold = 0.9)
    {
        var dict = ToGenomesDict(genomes);
        var clusters = PanGenomeAnalyzer.ClusterGenes(dict, identityThreshold).ToList();
        return new ClusterGenesResult(clusters);
    }

    /// <summary>
    /// Build a per-genome gene presence/absence matrix from a set of clusters.
    /// </summary>
    [McpServerTool(Name = "gene_presence_absence_matrix", Title = "Pan-genome — Presence/Absence Matrix", ReadOnly = true)]
    [Description("Build a per-genome gene presence/absence matrix against a list of gene clusters. Each row carries flattened (cluster id → present) entries.")]
    public static PresenceAbsenceMatrixResult GenePresenceAbsenceMatrix(
        [Description("Genomes (id + ordered list of genes).")] IReadOnlyList<GenomeInput> genomes,
        [Description("Gene clusters (e.g. from cluster_genes).")]
        IReadOnlyList<PanGenomeAnalyzer.GeneCluster> clusters)
    {
        var dict = ToGenomesDict(genomes);
        var rows = PanGenomeAnalyzer.CreatePresenceAbsenceMatrix(dict, clusters)
            .Select(r => new GenePresenceRowDto(
                r.GenomeId,
                r.GenePresence.Select(kv => new GenePresenceItem(kv.Key, kv.Value)).ToList(),
                r.TotalGenes,
                r.PresentGenes))
            .ToList();
        return new PresenceAbsenceMatrixResult(rows);
    }

    /// <summary>
    /// Fit Heaps' law to the pan-genome new-gene-discovery curve via permuted
    /// presence/absence orderings (Tettelin et al. 2008; micropan heaps()).
    /// </summary>
    [McpServerTool(Name = "fit_heaps_law", Title = "Pan-genome — Fit Heaps' Law", ReadOnly = true)]
    [Description("Fit Heaps' law n(N) = Intercept · N^(-Alpha) to the permuted new-gene-discovery curve (Tettelin 2008; micropan heaps). Pan-genome is open when Alpha < 1.")]
    public static HeapsLawFitDto FitHeapsLaw(
        [Description("Genomes (id + ordered list of genes).")] IReadOnlyList<GenomeInput> genomes,
        [Description("Sequence-identity threshold for ortholog clustering (default 0.9).")] double identityThreshold = 0.9,
        [Description("Number of random genome-order permutations to pool over (default 100).")] int permutations = 100)
    {
        var dict = ToGenomesDict(genomes);
        var fit = PanGenomeAnalyzer.FitHeapsLaw(dict, identityThreshold, permutations);
        return new HeapsLawFitDto(fit.Intercept, fit.Alpha, fit.IsOpen);
    }

    /// <summary>
    /// Filter clusters down to those present in at least <paramref name="threshold"/>·
    /// <paramref name="totalGenomes"/> genomes (the "core" set).
    /// </summary>
    [McpServerTool(Name = "core_gene_clusters", Title = "Pan-genome — Core Gene Clusters", ReadOnly = true)]
    [Description("Filter gene clusters down to the core set: those present in at least floor(threshold * totalGenomes) genomes.")]
    public static CoreGeneClustersResult CoreGeneClusters(
        [Description("All gene clusters.")]
        IReadOnlyList<PanGenomeAnalyzer.GeneCluster> clusters,
        [Description("Total number of genomes in the analysis.")] int totalGenomes,
        [Description("Core fraction threshold (default 0.99).")] double threshold = 0.99)
    {
        var core = PanGenomeAnalyzer.GetCoreGeneClusters(clusters, totalGenomes, threshold).ToList();
        return new CoreGeneClustersResult(core);
    }

    /// <summary>
    /// Concatenate a single genome's representative sequences for the supplied
    /// core clusters into a per-genome alignment block.
    /// </summary>
    [McpServerTool(Name = "core_genome_alignment", Title = "Pan-genome — Core Genome Alignment", ReadOnly = true)]
    [Description("Concatenate a single genome's representative sequences for the supplied core clusters into a per-genome alignment block.")]
    public static CoreGenomeAlignmentResult CoreGenomeAlignment(
        [Description("Genomes (id + ordered list of genes).")] IReadOnlyList<GenomeInput> genomes,
        [Description("Core gene clusters (e.g. from core_gene_clusters).")]
        IReadOnlyList<PanGenomeAnalyzer.GeneCluster> coreClusters,
        [Description("Genome id to extract the alignment block for.")] string genomeId)
    {
        var dict = ToGenomesDict(genomes);
        var alignment = PanGenomeAnalyzer.CreateCoreGenomeAlignment(dict, coreClusters, genomeId);
        return new CoreGenomeAlignmentResult(alignment);
    }

    /// <summary>
    /// Summarise accessory clusters (present in &gt;1 but &lt; all genomes) with
    /// genome membership and frequency.
    /// </summary>
    [McpServerTool(Name = "accessory_genes", Title = "Pan-genome — Accessory Genes", ReadOnly = true)]
    [Description("Summarise accessory clusters (present in >1 but not all genomes) with their genome membership and frequency = genomeCount / totalGenomes.")]
    public static AccessoryGenesResult AccessoryGenes(
        [Description("All gene clusters.")]
        IReadOnlyList<PanGenomeAnalyzer.GeneCluster> clusters,
        [Description("Total number of genomes in the analysis.")] int totalGenomes)
    {
        var items = PanGenomeAnalyzer.AnalyzeAccessoryGenes(clusters, totalGenomes)
            .Select(t => new AccessoryGeneItem(t.ClusterId, t.GenomesWithGene.ToList(), t.Frequency))
            .ToList();
        return new AccessoryGenesResult(items);
    }

    /// <summary>
    /// List, per genome, the cluster ids that occur only in that genome
    /// (singleton accessory genes).
    /// </summary>
    [McpServerTool(Name = "find_genome_specific_genes", Title = "Pan-genome — Genome-Specific Genes", ReadOnly = true)]
    [Description("For each genome, list the cluster ids that occur only in that genome (singleton accessory clusters).")]
    public static GenomeSpecificGenesResult FindGenomeSpecificGenes(
        [Description("Genomes (id + ordered list of genes).")] IReadOnlyList<GenomeInput> genomes,
        [Description("Gene clusters covering the supplied genomes.")]
        IReadOnlyList<PanGenomeAnalyzer.GeneCluster> clusters)
    {
        var dict = ToGenomesDict(genomes);
        var items = PanGenomeAnalyzer.FindGenomeSpecificGenes(dict, clusters)
            .Select(t => new GenomeSpecificGeneItem(t.GenomeId, t.UniqueGeneIds.ToList()))
            .ToList();
        return new GenomeSpecificGenesResult(items);
    }

    /// <summary>
    /// Pick phylogenetic-marker clusters from the core set: single-copy core clusters
    /// (present in all <paramref name="totalGenomes"/> genomes with exactly one gene per
    /// genome) that contain at least one parsimony-informative site, ranked by descending
    /// parsimony-informative-site count and capped at <paramref name="maxMarkers"/>
    /// (Ding et al. 2018, panX; Page et al. 2015, Roary).
    /// </summary>
    [McpServerTool(Name = "select_phylogenetic_markers", Title = "Pan-genome — Select Phylogenetic Markers", ReadOnly = true)]
    [Description("Pick phylogenetic markers: single-copy core clusters (present in all genomes with exactly one gene each) with >= 1 parsimony-informative site, ranked by descending parsimony-informative-site count, capped at maxMarkers (panX/Roary).")]
    public static SelectPhylogeneticMarkersResult SelectPhylogeneticMarkers(
        [Description("Genomes (id + ordered list of genes), used to recover each cluster's member sequences.")]
        IReadOnlyList<GenomeInput> genomes,
        [Description("Core gene clusters to filter (e.g. from core_gene_clusters).")]
        IReadOnlyList<PanGenomeAnalyzer.GeneCluster> coreClusters,
        [Description("Total number of genomes in the analysis.")] int totalGenomes,
        [Description("Maximum number of markers to return (default 100).")] int maxMarkers = 100)
    {
        var dict = ToGenomesDict(genomes);
        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(dict, coreClusters, totalGenomes, maxMarkers).ToList();
        return new SelectPhylogeneticMarkersResult(markers);
    }

    #endregion

    #region Helpers

    private static Seqeron.Genomics.Metagenomics.TaxonomyTree ToTaxonomyTree(
        IReadOnlyList<TaxonNodeInput> nodes)
    {
        var taxa = nodes.Select(n =>
            new Seqeron.Genomics.Metagenomics.TaxonNode(n.Id, n.Name, n.Rank, n.ParentId));
        return new Seqeron.Genomics.Metagenomics.TaxonomyTree(taxa);
    }

    private static List<AbundanceItem> ToAbundanceList(IReadOnlyDictionary<string, double> dict)
    {
        var list = new List<AbundanceItem>(dict.Count);
        foreach (var kv in dict)
            list.Add(new AbundanceItem(kv.Key, kv.Value));
        return list;
    }

    private static Dictionary<string, double> ToAbundanceDict(IReadOnlyList<AbundanceItem> items)
    {
        var dict = new Dictionary<string, double>(items.Count);
        foreach (var item in items)
            dict[item.Name] = item.Fraction;
        return dict;
    }

    private static Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> ToGenomesDict(
        IReadOnlyList<GenomeInput> genomes)
    {
        var dict = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>(genomes.Count);
        foreach (var g in genomes)
        {
            var genes = g.Genes.Select(x => (x.GeneId, x.Sequence)).ToList();
            dict[g.GenomeId] = genes;
        }
        return dict;
    }

    #endregion
}
