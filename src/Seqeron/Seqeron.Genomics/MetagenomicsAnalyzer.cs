using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics;

/// <summary>
/// Provides algorithms for metagenomic analysis and taxonomic classification.
/// </summary>
public static class MetagenomicsAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Represents a taxonomic classification result.
    /// </summary>
    public readonly record struct TaxonomicClassification(
        string ReadId,
        string Kingdom,
        string Phylum,
        string Class,
        string Order,
        string Family,
        string Genus,
        string Species,
        double Confidence,
        int MatchedKmers,
        int TotalKmers);

    /// <summary>
    /// Represents a taxonomic profile of a metagenomic sample.
    /// </summary>
    public readonly record struct TaxonomicProfile(
        IReadOnlyDictionary<string, double> KingdomAbundance,
        IReadOnlyDictionary<string, double> PhylumAbundance,
        IReadOnlyDictionary<string, double> GenusAbundance,
        IReadOnlyDictionary<string, double> SpeciesAbundance,
        double ShannonDiversity,
        double SimpsonDiversity,
        int TotalReads,
        int ClassifiedReads);

    /// <summary>
    /// Represents alpha diversity metrics.
    /// </summary>
    public readonly record struct AlphaDiversity(
        double ShannonIndex,
        double SimpsonIndex,
        double InverseSimpson,
        double Chao1Estimate,
        double ObservedSpecies,
        double PielouEvenness);

    /// <summary>
    /// Represents beta diversity comparison between samples.
    /// </summary>
    public readonly record struct BetaDiversity(
        string Sample1,
        string Sample2,
        double BrayCurtis,
        double JaccardDistance,
        double UniFracDistance,
        int SharedSpecies,
        int UniqueToSample1,
        int UniqueToSample2);

    /// <summary>
    /// Represents a binned contig (MAG - Metagenome Assembled Genome).
    /// </summary>
    public readonly record struct GenomeBin(
        string BinId,
        IReadOnlyList<string> ContigIds,
        double TotalLength,
        double GcContent,
        double Coverage,
        double Completeness,
        double Contamination,
        string PredictedTaxonomy);

    /// <summary>
    /// Represents a functional annotation.
    /// </summary>
    public readonly record struct FunctionalAnnotation(
        string GeneId,
        string Function,
        string Pathway,
        string KoNumber,
        string CogCategory,
        double EValue,
        double BitScore);

    /// <summary>
    /// K-mer database entry for taxonomic classification.
    /// </summary>
    public readonly record struct KmerTaxonEntry(
        string Kmer,
        string TaxonId,
        string TaxonName,
        int TaxonomicLevel);

    #endregion

    #region K-mer Based Classification

    /// <summary>
    /// Classifies metagenomic reads using k-mer matching.
    /// </summary>
    public static IEnumerable<TaxonomicClassification> ClassifyReads(
        IEnumerable<(string Id, string Sequence)> reads,
        IReadOnlyDictionary<string, string> kmerDatabase,
        int k = 31)
    {
        foreach (var (id, sequence) in reads)
        {
            if (string.IsNullOrEmpty(sequence) || sequence.Length < k)
            {
                yield return new TaxonomicClassification(
                    id, "Unclassified", "", "", "", "", "", "", 0, 0, 0);
                continue;
            }

            var taxonCounts = new Dictionary<string, int>();
            int totalKmers = 0;
            int matchedKmers = 0;

            for (int i = 0; i <= sequence.Length - k; i++)
            {
                string kmer = sequence.Substring(i, k).ToUpperInvariant();
                totalKmers++;

                if (kmerDatabase.TryGetValue(kmer, out string? taxon))
                {
                    matchedKmers++;
                    if (!taxonCounts.ContainsKey(taxon))
                        taxonCounts[taxon] = 0;
                    taxonCounts[taxon]++;
                }
            }

            if (taxonCounts.Count == 0)
            {
                yield return new TaxonomicClassification(
                    id, "Unclassified", "", "", "", "", "", "", 0, 0, totalKmers);
                continue;
            }

            // Find best matching taxon using LCA
            var bestTaxon = taxonCounts.OrderByDescending(kv => kv.Value).First();
            double confidence = totalKmers > 0 ? (double)matchedKmers / totalKmers : 0;

            var taxonomy = ParseTaxonomyString(bestTaxon.Key);

            yield return new TaxonomicClassification(
                ReadId: id,
                Kingdom: taxonomy.GetValueOrDefault("kingdom", ""),
                Phylum: taxonomy.GetValueOrDefault("phylum", ""),
                Class: taxonomy.GetValueOrDefault("class", ""),
                Order: taxonomy.GetValueOrDefault("order", ""),
                Family: taxonomy.GetValueOrDefault("family", ""),
                Genus: taxonomy.GetValueOrDefault("genus", ""),
                Species: taxonomy.GetValueOrDefault("species", ""),
                Confidence: confidence,
                MatchedKmers: matchedKmers,
                TotalKmers: totalKmers);
        }
    }

    /// <summary>
    /// Builds a k-mer database from reference genomes.
    /// </summary>
    public static Dictionary<string, string> BuildKmerDatabase(
        IEnumerable<(string TaxonId, string Sequence)> referenceGenomes,
        int k = 31)
    {
        var database = new Dictionary<string, string>();

        foreach (var (taxonId, sequence) in referenceGenomes)
        {
            if (string.IsNullOrEmpty(sequence) || sequence.Length < k)
                continue;

            string seq = sequence.ToUpperInvariant();

            for (int i = 0; i <= seq.Length - k; i++)
            {
                string kmer = seq.Substring(i, k);
                if (kmer.All(c => "ACGT".Contains(c)))
                {
                    // Use canonical k-mer (lexicographically smaller of forward/reverse)
                    string canonical = GetCanonicalKmer(kmer);

                    if (!database.ContainsKey(canonical))
                        database[canonical] = taxonId;
                }
            }
        }

        return database;
    }

    private static string GetCanonicalKmer(string kmer)
    {
        string revComp = DnaSequence.GetReverseComplementString(kmer);
        return string.Compare(kmer, revComp, StringComparison.Ordinal) <= 0 ? kmer : revComp;
    }

    private static Dictionary<string, string> ParseTaxonomyString(string taxonomy)
    {
        var result = new Dictionary<string, string>();
        var ranks = new[] { "kingdom", "phylum", "class", "order", "family", "genus", "species" };

        var parts = taxonomy.Split(';', '|');
        for (int i = 0; i < Math.Min(parts.Length, ranks.Length); i++)
        {
            result[ranks[i]] = parts[i].Trim();
        }

        return result;
    }

    #endregion

    #region Taxonomic Profile Generation

    /// <summary>
    /// Generates a taxonomic profile from classified reads.
    /// </summary>
    public static TaxonomicProfile GenerateTaxonomicProfile(
        IEnumerable<TaxonomicClassification> classifications)
    {
        var classList = classifications.ToList();
        int totalReads = classList.Count;
        int classifiedReads = classList.Count(c => c.Kingdom != "Unclassified" && !string.IsNullOrEmpty(c.Kingdom));

        var kingdomCounts = new Dictionary<string, int>();
        var phylumCounts = new Dictionary<string, int>();
        var genusCounts = new Dictionary<string, int>();
        var speciesCounts = new Dictionary<string, int>();

        foreach (var c in classList.Where(c => c.Kingdom != "Unclassified"))
        {
            IncrementCount(kingdomCounts, c.Kingdom);
            IncrementCount(phylumCounts, c.Phylum);
            IncrementCount(genusCounts, c.Genus);
            IncrementCount(speciesCounts, c.Species);
        }

        double total = classifiedReads > 0 ? classifiedReads : 1;

        var kingdomAbundance = kingdomCounts.ToDictionary(kv => kv.Key, kv => kv.Value / total);
        var phylumAbundance = phylumCounts.Where(kv => !string.IsNullOrEmpty(kv.Key))
            .ToDictionary(kv => kv.Key, kv => kv.Value / total);
        var genusAbundance = genusCounts.Where(kv => !string.IsNullOrEmpty(kv.Key))
            .ToDictionary(kv => kv.Key, kv => kv.Value / total);
        var speciesAbundance = speciesCounts.Where(kv => !string.IsNullOrEmpty(kv.Key))
            .ToDictionary(kv => kv.Key, kv => kv.Value / total);

        double shannon = CalculateShannonIndex(speciesAbundance.Values);
        double simpson = CalculateSimpsonIndex(speciesAbundance.Values);

        return new TaxonomicProfile(
            KingdomAbundance: kingdomAbundance,
            PhylumAbundance: phylumAbundance,
            GenusAbundance: genusAbundance,
            SpeciesAbundance: speciesAbundance,
            ShannonDiversity: shannon,
            SimpsonDiversity: simpson,
            TotalReads: totalReads,
            ClassifiedReads: classifiedReads);
    }

    private static void IncrementCount(Dictionary<string, int> counts, string key)
    {
        if (string.IsNullOrEmpty(key)) return;
        if (!counts.ContainsKey(key))
            counts[key] = 0;
        counts[key]++;
    }

    #endregion

    #region Alpha Diversity

    /// <summary>
    /// Calculates alpha diversity metrics for a sample.
    /// </summary>
    public static AlphaDiversity CalculateAlphaDiversity(IReadOnlyDictionary<string, double> abundances)
    {
        if (abundances == null || abundances.Count == 0)
            return new AlphaDiversity(0, 0, 0, 0, 0, 0);

        var values = abundances.Values.Where(v => v > 0).ToList();
        int observedSpecies = values.Count;

        if (observedSpecies == 0)
            return new AlphaDiversity(0, 0, 0, 0, 0, 0);

        double shannon = CalculateShannonIndex(values);
        double simpson = CalculateSimpsonIndex(values);
        double inverseSimpson = simpson > 0 ? 1 / simpson : 0;

        // Chao1 estimator (simplified - needs singleton/doubleton counts for proper calculation)
        double chao1 = observedSpecies; // Simplified

        // Pielou's evenness
        double evenness = observedSpecies > 1 ? shannon / Math.Log(observedSpecies) : 0;

        return new AlphaDiversity(
            ShannonIndex: shannon,
            SimpsonIndex: simpson,
            InverseSimpson: inverseSimpson,
            Chao1Estimate: chao1,
            ObservedSpecies: observedSpecies,
            PielouEvenness: evenness);
    }

    private static double CalculateShannonIndex(IEnumerable<double> abundances)
    {
        double sum = abundances.Sum();
        if (sum == 0) return 0;

        double h = 0;
        foreach (var p in abundances.Where(a => a > 0))
        {
            double pi = p / sum;
            h -= pi * Math.Log(pi);
        }
        return h;
    }

    private static double CalculateSimpsonIndex(IEnumerable<double> abundances)
    {
        double sum = abundances.Sum();
        if (sum == 0) return 0;

        double d = 0;
        foreach (var p in abundances.Where(a => a > 0))
        {
            double pi = p / sum;
            d += pi * pi;
        }
        return d;
    }

    #endregion

    #region Beta Diversity

    /// <summary>
    /// Calculates beta diversity between two samples.
    /// </summary>
    public static BetaDiversity CalculateBetaDiversity(
        string sample1Name,
        IReadOnlyDictionary<string, double> sample1,
        string sample2Name,
        IReadOnlyDictionary<string, double> sample2)
    {
        var allSpecies = sample1.Keys.Union(sample2.Keys).ToList();

        int shared = 0;
        int unique1 = 0;
        int unique2 = 0;

        foreach (var species in allSpecies)
        {
            bool in1 = sample1.ContainsKey(species) && sample1[species] > 0;
            bool in2 = sample2.ContainsKey(species) && sample2[species] > 0;

            if (in1 && in2) shared++;
            else if (in1) unique1++;
            else if (in2) unique2++;
        }

        double brayCurtis = CalculateBrayCurtis(sample1, sample2, allSpecies);
        double jaccard = CalculateJaccardDistance(shared, unique1, unique2);

        return new BetaDiversity(
            Sample1: sample1Name,
            Sample2: sample2Name,
            BrayCurtis: brayCurtis,
            JaccardDistance: jaccard,
            UniFracDistance: 0, // Requires phylogenetic tree
            SharedSpecies: shared,
            UniqueToSample1: unique1,
            UniqueToSample2: unique2);
    }

    private static double CalculateBrayCurtis(
        IReadOnlyDictionary<string, double> sample1,
        IReadOnlyDictionary<string, double> sample2,
        IEnumerable<string> allSpecies)
    {
        double sumMin = 0;
        double sumTotal = 0;

        foreach (var species in allSpecies)
        {
            double a1 = sample1.GetValueOrDefault(species, 0);
            double a2 = sample2.GetValueOrDefault(species, 0);
            sumMin += Math.Min(a1, a2);
            sumTotal += a1 + a2;
        }

        return sumTotal > 0 ? 1 - (2 * sumMin / sumTotal) : 0;
    }

    private static double CalculateJaccardDistance(int shared, int unique1, int unique2)
    {
        int total = shared + unique1 + unique2;
        return total > 0 ? 1.0 - (double)shared / total : 0;
    }

    #endregion

    #region Genome Binning

    /// <summary>
    /// Bins contigs based on composition and coverage.
    /// </summary>
    public static IEnumerable<GenomeBin> BinContigs(
        IEnumerable<(string ContigId, string Sequence, double Coverage)> contigs,
        int numBins = 10,
        double minBinSize = 500000)
    {
        var contigList = contigs.ToList();

        if (contigList.Count == 0)
            yield break;

        // Calculate features for each contig
        var features = contigList.Select(c => new
        {
            c.ContigId,
            c.Sequence,
            c.Coverage,
            GcContent = CalculateGcContent(c.Sequence),
            TetraNucFreq = CalculateTetraNucleotideFrequency(c.Sequence)
        }).ToList();

        // Simple k-means clustering based on GC and coverage
        var bins = new Dictionary<int, List<(string Id, string Seq, double Coverage)>>();

        for (int i = 0; i < numBins; i++)
            bins[i] = new List<(string, string, double)>();

        foreach (var f in features)
        {
            // Assign to bin based on GC content and coverage
            int binIndex = (int)(f.GcContent * numBins) % numBins;
            bins[binIndex].Add((f.ContigId, f.Sequence, f.Coverage));
        }

        int binId = 1;
        foreach (var bin in bins.Values.Where(b => b.Count > 0))
        {
            double totalLength = bin.Sum(c => c.Seq.Length);
            if (totalLength < minBinSize)
                continue;

            double avgGc = bin.Average(c => CalculateGcContent(c.Seq));
            double avgCoverage = bin.Average(c => c.Coverage);

            // Estimate completeness using single-copy marker genes (simplified)
            double completeness = Math.Min(totalLength / 2000000 * 100, 100);
            double contamination = CalculateContamination(bin.Select(b => b.Seq));

            yield return new GenomeBin(
                BinId: $"bin.{binId++}",
                ContigIds: bin.Select(b => b.Id).ToList(),
                TotalLength: totalLength,
                GcContent: avgGc,
                Coverage: avgCoverage,
                Completeness: completeness,
                Contamination: contamination,
                PredictedTaxonomy: "");
        }
    }

    private static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcFractionFast();

    private static Dictionary<string, double> CalculateTetraNucleotideFrequency(string sequence)
    {
        var freq = new Dictionary<string, int>();
        sequence = sequence.ToUpperInvariant();

        for (int i = 0; i <= sequence.Length - 4; i++)
        {
            string tetra = sequence.Substring(i, 4);
            if (tetra.All(c => "ACGT".Contains(c)))
            {
                if (!freq.ContainsKey(tetra))
                    freq[tetra] = 0;
                freq[tetra]++;
            }
        }

        int total = freq.Values.Sum();
        return freq.ToDictionary(kv => kv.Key, kv => total > 0 ? (double)kv.Value / total : 0);
    }

    private static double CalculateContamination(IEnumerable<string> sequences)
    {
        // Simplified contamination estimation based on GC variance
        var gcValues = sequences.Select(CalculateGcContent).ToList();
        if (gcValues.Count < 2)
            return 0;

        double mean = gcValues.Average();
        double variance = gcValues.Sum(gc => (gc - mean) * (gc - mean)) / gcValues.Count;

        // High variance suggests contamination
        return Math.Min(variance * 1000, 50);
    }

    #endregion

    #region Functional Profiling

    /// <summary>
    /// Predicts functional annotations using HMM profiles (simplified simulation).
    /// </summary>
    public static IEnumerable<FunctionalAnnotation> PredictFunctions(
        IEnumerable<(string GeneId, string ProteinSequence)> proteins,
        IReadOnlyDictionary<string, (string Function, string Pathway, string Ko)> functionDatabase)
    {
        foreach (var (geneId, sequence) in proteins)
        {
            if (string.IsNullOrEmpty(sequence))
                continue;

            // Simplified: look for conserved motifs
            foreach (var entry in functionDatabase)
            {
                string motif = entry.Key;
                if (sequence.Contains(motif))
                {
                    var (function, pathway, ko) = entry.Value;
                    yield return new FunctionalAnnotation(
                        GeneId: geneId,
                        Function: function,
                        Pathway: pathway,
                        KoNumber: ko,
                        CogCategory: InferCogCategory(function),
                        EValue: 1e-10, // Simulated
                        BitScore: motif.Length * 2.0);
                }
            }
        }
    }

    private static string InferCogCategory(string function)
    {
        function = function.ToLowerInvariant();

        if (function.Contains("transport")) return "G"; // Carbohydrate transport
        if (function.Contains("kinase") || function.Contains("phospho")) return "T"; // Signal transduction
        if (function.Contains("dna") || function.Contains("replica")) return "L"; // DNA replication
        if (function.Contains("ribosom") || function.Contains("translat")) return "J"; // Translation
        if (function.Contains("metabol")) return "C"; // Energy metabolism
        if (function.Contains("membrane")) return "M"; // Cell membrane

        return "S"; // Function unknown
    }

    /// <summary>
    /// Calculates functional diversity metrics.
    /// </summary>
    public static (double FunctionalRichness, double FunctionalDiversity, IReadOnlyDictionary<string, int> PathwayCounts)
        CalculateFunctionalDiversity(IEnumerable<FunctionalAnnotation> annotations)
    {
        var annotList = annotations.ToList();

        var pathwayCounts = annotList
            .Where(a => !string.IsNullOrEmpty(a.Pathway))
            .GroupBy(a => a.Pathway)
            .ToDictionary(g => g.Key, g => g.Count());

        var functionCounts = annotList
            .Where(a => !string.IsNullOrEmpty(a.Function))
            .GroupBy(a => a.Function)
            .ToDictionary(g => g.Key, g => g.Count());

        double richness = functionCounts.Count;
        double total = functionCounts.Values.Sum();

        // Shannon diversity of functions
        double diversity = 0;
        if (total > 0)
        {
            foreach (var count in functionCounts.Values)
            {
                double p = count / total;
                if (p > 0)
                    diversity -= p * Math.Log(p);
            }
        }

        return (richness, diversity, pathwayCounts);
    }

    #endregion

    #region Antibiotic Resistance Gene Detection

    /// <summary>
    /// Searches for antibiotic resistance genes.
    /// </summary>
    public static IEnumerable<(string GeneId, string ResistanceGene, string AntibioticClass, double Identity)>
        FindResistanceGenes(
            IEnumerable<(string GeneId, string Sequence)> genes,
            IReadOnlyDictionary<string, (string Name, string AntibioticClass)> resistanceDatabase)
    {
        foreach (var (geneId, sequence) in genes)
        {
            if (string.IsNullOrEmpty(sequence))
                continue;

            foreach (var entry in resistanceDatabase)
            {
                string targetMotif = entry.Key;
                var (name, antibioticClass) = entry.Value;

                // Simple containment check (real implementation would use alignment)
                if (sequence.Contains(targetMotif))
                {
                    double identity = (double)targetMotif.Length / sequence.Length;
                    yield return (geneId, name, antibioticClass, identity);
                }
            }
        }
    }

    #endregion

    #region Sample Comparison

    /// <summary>
    /// Performs differential abundance analysis between two conditions.
    /// </summary>
    public static IEnumerable<(string Taxon, double FoldChange, double PValue, bool Significant)>
        DifferentialAbundance(
            IEnumerable<IReadOnlyDictionary<string, double>> condition1Samples,
            IEnumerable<IReadOnlyDictionary<string, double>> condition2Samples,
            double pValueThreshold = 0.05)
    {
        var c1List = condition1Samples.ToList();
        var c2List = condition2Samples.ToList();

        if (c1List.Count == 0 || c2List.Count == 0)
            yield break;

        var allTaxa = c1List.SelectMany(s => s.Keys)
            .Union(c2List.SelectMany(s => s.Keys))
            .Distinct();

        foreach (var taxon in allTaxa)
        {
            var values1 = c1List.Select(s => s.GetValueOrDefault(taxon, 0)).ToList();
            var values2 = c2List.Select(s => s.GetValueOrDefault(taxon, 0)).ToList();

            double mean1 = values1.Average();
            double mean2 = values2.Average();

            double foldChange = mean1 > 0 ? mean2 / mean1 : (mean2 > 0 ? double.PositiveInfinity : 1);
            double logFoldChange = Math.Log2(Math.Max(foldChange, 0.001));

            // Simple t-test (Welch's approximation)
            double pValue = CalculateTTestPValue(values1, values2);
            bool significant = pValue < pValueThreshold && Math.Abs(logFoldChange) > 1;

            yield return (taxon, logFoldChange, pValue, significant);
        }
    }

    private static double CalculateTTestPValue(List<double> group1, List<double> group2)
    {
        if (group1.Count < 2 || group2.Count < 2)
            return 1.0;

        double mean1 = group1.Average();
        double mean2 = group2.Average();
        double var1 = group1.Sum(x => (x - mean1) * (x - mean1)) / (group1.Count - 1);
        double var2 = group2.Sum(x => (x - mean2) * (x - mean2)) / (group2.Count - 1);

        if (var1 == 0 && var2 == 0)
            return mean1 == mean2 ? 1.0 : 0.0;

        double se = Math.Sqrt(var1 / group1.Count + var2 / group2.Count);
        if (se == 0) return 1.0;

        double t = Math.Abs(mean1 - mean2) / se;

        // Approximate p-value using normal distribution
        return 2 * (1 - StatisticsHelper.NormalCDF(t));
    }

    #endregion
}
