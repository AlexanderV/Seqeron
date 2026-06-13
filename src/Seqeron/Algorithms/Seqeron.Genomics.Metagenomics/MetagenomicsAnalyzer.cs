using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Metagenomics;

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
    /// Result of a pathway over-representation (enrichment) test.
    /// <paramref name="PValue"/> is the right-tail hypergeometric probability
    /// P(X ≥ Overlap): the chance of seeing at least this many query genes in the
    /// pathway by random sampling from the background.
    /// </summary>
    public readonly record struct PathwayEnrichment(
        string Pathway,
        int Overlap,
        int PathwaySize,
        int QuerySize,
        int BackgroundSize,
        double PValue);

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
    /// Classifies metagenomic reads using k-mer matching with a flat best-hit rule.
    /// </summary>
    /// <remarks>
    /// Each read is assigned to the single taxon with the highest count of matching (canonical)
    /// k-mers. This is a BEST-HIT (highest k-mer count) classifier — it is NOT an LCA / Kraken
    /// weighted root-to-leaf classifier: there is no taxonomy tree, no lowest-common-ancestor
    /// resolution, and no per-rank weighting. Each database k-mer maps to exactly one taxon
    /// (see <see cref="BuildKmerDatabase"/>). Ties (two or more taxa with the same best k-mer count)
    /// resolve to an arbitrary best-count taxon determined by dictionary/enumeration ordering.
    /// Confidence is simply matched-best-taxon k-mers / total non-ambiguous k-mers.
    /// </remarks>
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

            for (int i = 0; i <= sequence.Length - k; i++)
            {
                string kmer = sequence.Substring(i, k).ToUpperInvariant();

                // Skip k-mers containing ambiguous nucleotides (Kraken-style filtering)
                if (!kmer.All(c => "ACGT".Contains(c)))
                    continue;

                // Use canonical k-mer (strand-independent) for database lookup
                string canonical = GetCanonicalKmer(kmer);
                totalKmers++;

                if (kmerDatabase.TryGetValue(canonical, out string? taxon))
                {
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

            // Best-hit rule: classify to the taxon with the most k-mer hits.
            // NOT an LCA — ties resolve to an arbitrary best-count taxon (enumeration order).
            var bestTaxon = taxonCounts.OrderByDescending(kv => kv.Value).First();
            // C = k-mers supporting the chosen classification, Q = non-ambiguous k-mers
            int matchedKmers = bestTaxon.Value;
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
    /// Builds a flat k-mer → taxon database from reference genomes.
    /// </summary>
    /// <remarks>
    /// Each canonical k-mer is mapped to exactly ONE taxon — the FIRST reference genome in which
    /// it is encountered (subsequent occurrences in other taxa are ignored). This is NOT a Kraken
    /// LCA database: there is no taxonomy tree and no lowest-common-ancestor assignment for k-mers
    /// shared across taxa. The resulting database supports only the flat best-hit classification
    /// performed by <see cref="ClassifyReads"/>.
    /// </remarks>
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
        var classified = classList.Where(c => c.Kingdom != "Unclassified" && !string.IsNullOrEmpty(c.Kingdom)).ToList();
        int classifiedReads = classified.Count;

        var kingdomCounts = new Dictionary<string, int>();
        var phylumCounts = new Dictionary<string, int>();
        var genusCounts = new Dictionary<string, int>();
        var speciesCounts = new Dictionary<string, int>();

        foreach (var c in classified)
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

        // Chao1 estimator — Chao (1984)
        double chao1 = CalculateChao1(values, observedSpecies);

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

    /// <summary>
    /// Chao1 richness estimator — Chao (1984).
    /// S_Chao1 = S_obs + f1²/(2·f2) when f2 > 0;
    /// S_Chao1 = S_obs + f1·(f1−1)/2 when f2 = 0 (bias-corrected).
    /// f1 = singletons (species with count = 1), f2 = doubletons (count = 2).
    /// Requires integer count data; for proportional data, returns S_obs.
    /// </summary>
    private static double CalculateChao1(List<double> values, int observedSpecies)
    {
        const double intTolerance = 1e-9;
        bool isCountData = values.All(v => Math.Abs(v - Math.Round(v)) < intTolerance);
        if (!isCountData)
            return observedSpecies;

        int f1 = values.Count(v => Math.Abs(v - 1.0) < intTolerance);
        int f2 = values.Count(v => Math.Abs(v - 2.0) < intTolerance);

        if (f1 == 0)
            return observedSpecies;

        if (f2 > 0)
            return observedSpecies + (double)(f1 * f1) / (2 * f2);

        // Bias-corrected form when f2 = 0
        return observedSpecies + (double)(f1 * (f1 - 1)) / 2;
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
    /// Bins contigs using k-means clustering on compositional and coverage features.
    /// Features: GC content, tetranucleotide frequency (per TETRA/Teeling 2004), coverage.
    /// Quality: completeness from bin length / expected genome size,
    ///          contamination from within-bin GC standard deviation normalized by theoretical max.
    /// </summary>
    /// <param name="contigs">Input contigs with sequences and coverage.</param>
    /// <param name="numBins">Maximum number of bins (k for k-means).</param>
    /// <param name="minBinSize">Minimum total length for a bin to be reported.</param>
    /// <param name="expectedGenomeSize">Expected genome size in bp for completeness estimation.</param>
    public static IEnumerable<GenomeBin> BinContigs(
        IEnumerable<(string ContigId, string Sequence, double Coverage)> contigs,
        int numBins = 10,
        double minBinSize = 500000,
        double expectedGenomeSize = 4_000_000)
    {
        var contigList = contigs.ToList();

        if (contigList.Count == 0)
            yield break;

        // Calculate features for each contig: GC, coverage, TNF
        var features = contigList.Select(c => new ContigFeatures(
            GcContent: CalculateGcContent(c.Sequence),
            Coverage: c.Coverage,
            Tnf: CalculateTetraNucleotideFrequency(c.Sequence)
        )).ToList();

        // Normalize coverage to [0, 1] for distance computation
        double maxCov = features.Max(f => f.Coverage);
        if (maxCov <= 0) maxCov = 1;
        var normalized = features.Select(f => f with { Coverage = f.Coverage / maxCov }).ToList();

        // K-means clustering on [GC, normalized coverage, TNF Pearson distance]
        int k = Math.Min(numBins, contigList.Count);
        var assignments = KMeansCluster(normalized, k);

        // Group contigs by cluster
        var clusters = new Dictionary<int, List<int>>();
        for (int i = 0; i < assignments.Length; i++)
        {
            if (!clusters.ContainsKey(assignments[i]))
                clusters[assignments[i]] = new List<int>();
            clusters[assignments[i]].Add(i);
        }

        int binId = 1;
        foreach (var cluster in clusters.Values)
        {
            double totalLength = cluster.Sum(i => contigList[i].Sequence.Length);
            if (totalLength < minBinSize)
                continue;

            double avgGc = cluster.Average(i => features[i].GcContent);
            double avgCoverage = cluster.Average(i => contigList[i].Coverage);
            double completeness = Math.Min(totalLength / expectedGenomeSize * 100, 100);
            double contamination = CalculateContamination(cluster.Select(i => features[i].GcContent));

            yield return new GenomeBin(
                BinId: $"bin.{binId++}",
                ContigIds: cluster.Select(i => contigList[i].ContigId).ToList(),
                TotalLength: totalLength,
                GcContent: avgGc,
                Coverage: avgCoverage,
                Completeness: completeness,
                Contamination: contamination,
                PredictedTaxonomy: "");
        }
    }

    private readonly record struct ContigFeatures(
        double GcContent,
        double Coverage,
        Dictionary<string, double> Tnf);

    private static int[] KMeansCluster(List<ContigFeatures> features, int k, int maxIterations = 50)
    {
        int n = features.Count;
        if (n == 0 || k <= 0)
            return Array.Empty<int>();

        k = Math.Min(k, n);
        var assignments = new int[n];

        // Initialize centroids by spreading across GC-sorted data (deterministic, diverse)
        var sortedIndices = Enumerable.Range(0, n).OrderBy(i => features[i].GcContent).ToList();
        var centroids = new ContigFeatures[k];
        for (int c = 0; c < k; c++)
        {
            int idx = sortedIndices[c * n / k];
            centroids[c] = features[idx];
        }

        for (int iter = 0; iter < maxIterations; iter++)
        {
            bool changed = false;

            // Assignment step: assign each point to nearest centroid
            for (int i = 0; i < n; i++)
            {
                double minDist = double.MaxValue;
                int best = 0;
                for (int c = 0; c < k; c++)
                {
                    double dist = CompositeDistance(features[i], centroids[c]);
                    if (dist < minDist)
                    {
                        minDist = dist;
                        best = c;
                    }
                }
                if (assignments[i] != best)
                {
                    assignments[i] = best;
                    changed = true;
                }
            }

            if (!changed)
                break;

            // Update step: recompute centroids
            for (int c = 0; c < k; c++)
            {
                var members = Enumerable.Range(0, n).Where(i => assignments[i] == c).ToList();
                if (members.Count == 0)
                    continue;

                double avgGc = members.Average(i => features[i].GcContent);
                double avgCov = members.Average(i => features[i].Coverage);

                // Average TNF vector across cluster members
                var avgTnf = new Dictionary<string, double>();
                var allKeys = members.SelectMany(i => features[i].Tnf.Keys).Distinct().ToList();
                foreach (var key in allKeys)
                    avgTnf[key] = members.Average(i => features[i].Tnf.GetValueOrDefault(key, 0));

                centroids[c] = new ContigFeatures(avgGc, avgCov, avgTnf);
            }
        }

        return assignments;
    }

    /// <summary>
    /// Composite distance: GC difference + coverage difference + TNF Pearson distance.
    /// Per literature, binning uses GC content, coverage, and tetranucleotide frequencies
    /// (Wikipedia: Binning metagenomics; Teeling 2004 TETRA).
    /// </summary>
    private static double CompositeDistance(ContigFeatures a, ContigFeatures b)
    {
        double gcDist = Math.Abs(a.GcContent - b.GcContent);
        double covDist = Math.Abs(a.Coverage - b.Coverage);
        double tnfDist = TnfPearsonDistance(a.Tnf, b.Tnf);
        return gcDist + covDist + tnfDist;
    }

    /// <summary>
    /// TNF distance based on Pearson correlation, per TETRA methodology (Teeling 2004).
    /// Returns (1 - r) / 2, normalized to [0, 1].
    /// </summary>
    private static double TnfPearsonDistance(Dictionary<string, double> a, Dictionary<string, double> b)
    {
        var allKeys = a.Keys.Union(b.Keys).ToList();
        if (allKeys.Count == 0) return 1.0;

        double[] va = allKeys.Select(k => a.GetValueOrDefault(k, 0)).ToArray();
        double[] vb = allKeys.Select(k => b.GetValueOrDefault(k, 0)).ToArray();

        double meanA = va.Average(), meanB = vb.Average();
        double sumAB = 0, sumA2 = 0, sumB2 = 0;
        for (int i = 0; i < va.Length; i++)
        {
            double da = va[i] - meanA, db = vb[i] - meanB;
            sumAB += da * db;
            sumA2 += da * da;
            sumB2 += db * db;
        }

        double denom = Math.Sqrt(sumA2) * Math.Sqrt(sumB2);
        double r = denom > 0 ? sumAB / denom : 0;
        return (1.0 - r) / 2.0;
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

    /// <summary>
    /// Estimates contamination from within-bin GC standard deviation, normalized by
    /// the theoretical maximum (0.5) for values in [0, 1].
    /// Low GC variance indicates a single taxonomic source (Parks et al. 2014).
    /// </summary>
    private static double CalculateContamination(IEnumerable<double> gcValues)
    {
        var gcList = gcValues.ToList();
        if (gcList.Count < 2)
            return 0;

        double mean = gcList.Average();
        double variance = gcList.Sum(gc => (gc - mean) * (gc - mean)) / gcList.Count;
        double stdDev = Math.Sqrt(variance);

        // Normalize by theoretical max std dev (0.5) for values in [0, 1]
        return Math.Min(stdDev / 0.5 * 100, 100);
    }

    #endregion

    #region Functional Profiling

    // Ungapped Karlin-Altschul parameters for the BLOSUM62 scoring system, taken from
    // the first (ungapped) row of `blosum62_values` in NCBI's BLAST core blast_stat.c.
    // Source: NCBI C BLAST Toolkit, algo/blast/core/blast_stat.c.
    private const double Blosum62UngappedLambda = 0.3176;
    private const double Blosum62UngappedK = 0.134;

    // ln(2): the bit-score normalization divides by ln 2 so the score is expressed in bits.
    // Source: Altschul et al., NCBI BLAST tutorial, "The Statistics of Sequence Similarity Scores".
    private static readonly double Ln2 = Math.Log(2.0);

    // BLOSUM62 diagonal (self-match) substitution scores per amino acid.
    // An exact (ungapped) self-alignment of a protein segment scores the sum of these
    // entries over its residues. Source: NCBI BLOSUM62 matrix file (blast/matrices/BLOSUM62).
    private static readonly IReadOnlyDictionary<char, int> Blosum62Diagonal = new Dictionary<char, int>
    {
        ['A'] = 4, ['R'] = 5, ['N'] = 6, ['D'] = 6, ['C'] = 9, ['Q'] = 5, ['E'] = 5,
        ['G'] = 6, ['H'] = 8, ['I'] = 4, ['L'] = 4, ['K'] = 5, ['M'] = 5, ['F'] = 6,
        ['P'] = 7, ['S'] = 4, ['T'] = 5, ['W'] = 11, ['Y'] = 7, ['V'] = 4,
    };

    /// <summary>
    /// Predicts functional annotations by homology-based annotation transfer: each query
    /// protein is matched against a signature database; for every database signature that
    /// occurs exactly in the protein, the ungapped BLOSUM62 raw score of the matched
    /// segment is converted into a bit score and an E-value using Karlin-Altschul
    /// statistics, and the annotation (function, pathway, KO) of the single best hit
    /// (lowest E-value) is transferred to the gene.
    /// </summary>
    /// <remarks>
    /// Bit score: S' = (λ·S − ln K) / ln 2; E-value: E = K·m·n·e^(−λ·S), equivalently
    /// E = m·n·2^(−S'), with the ungapped BLOSUM62 λ = 0.3176, K = 0.134 (NCBI blast_stat.c)
    /// and the matched-segment raw score S summed from BLOSUM62 diagonal scores. m and n are
    /// the lengths of the query protein and the matched signature. Source: Altschul et al.,
    /// NCBI BLAST tutorial.
    /// </remarks>
    /// <param name="proteins">Query genes as (GeneId, ProteinSequence) pairs (single-letter amino acids).</param>
    /// <param name="functionDatabase">Maps a signature (subsequence) to its (Function, Pathway, KO).</param>
    /// <returns>One best-hit <see cref="FunctionalAnnotation"/> per gene that matched at least one signature.</returns>
    /// <exception cref="ArgumentNullException">If <paramref name="proteins"/> or <paramref name="functionDatabase"/> is null.</exception>
    public static IEnumerable<FunctionalAnnotation> PredictFunctions(
        IEnumerable<(string GeneId, string ProteinSequence)> proteins,
        IReadOnlyDictionary<string, (string Function, string Pathway, string Ko)> functionDatabase)
    {
        if (proteins is null) throw new ArgumentNullException(nameof(proteins));
        if (functionDatabase is null) throw new ArgumentNullException(nameof(functionDatabase));

        return Iterate();

        IEnumerable<FunctionalAnnotation> Iterate()
        {
            foreach (var (geneId, sequence) in proteins)
            {
                if (string.IsNullOrWhiteSpace(sequence))
                    continue;

                FunctionalAnnotation? best = null;

                foreach (var entry in functionDatabase)
                {
                    string signature = entry.Key;
                    if (string.IsNullOrEmpty(signature) || !sequence.Contains(signature, StringComparison.Ordinal))
                        continue;

                    int rawScore = Blosum62SelfScore(signature);
                    double bitScore = FunctionalBitScore(rawScore);
                    // m·n is the search space: query length × matched-signature length.
                    double eValue = ExpectedValue(rawScore, sequence.Length, signature.Length);

                    var (function, pathway, ko) = entry.Value;
                    var candidate = new FunctionalAnnotation(
                        GeneId: geneId,
                        Function: function,
                        Pathway: pathway,
                        KoNumber: ko,
                        CogCategory: InferCogCategory(function),
                        EValue: eValue,
                        BitScore: bitScore);

                    // Best hit = lowest E-value (most significant). Source: BLAST E-value ranks significance.
                    if (best is null || candidate.EValue < best.Value.EValue)
                        best = candidate;
                }

                if (best is not null)
                    yield return best.Value;
            }
        }
    }

    /// <summary>
    /// Ungapped BLOSUM62 raw self-alignment score of a protein segment: the sum of the
    /// BLOSUM62 diagonal scores over its residues. Unknown residues contribute 0.
    /// Source: NCBI BLOSUM62 matrix file.
    /// </summary>
    public static int Blosum62SelfScore(string segment)
    {
        if (string.IsNullOrEmpty(segment))
            return 0;

        int score = 0;
        foreach (char residue in segment)
        {
            if (Blosum62Diagonal.TryGetValue(char.ToUpperInvariant(residue), out int s))
                score += s;
        }

        return score;
    }

    /// <summary>
    /// Bit (normalized) score from a raw ungapped BLOSUM62 alignment score:
    /// S' = (λ·S − ln K) / ln 2. Source: Altschul et al., NCBI BLAST tutorial.
    /// </summary>
    public static double FunctionalBitScore(double rawScore)
        => (Blosum62UngappedLambda * rawScore - Math.Log(Blosum62UngappedK)) / Ln2;

    /// <summary>
    /// Karlin-Altschul E-value of an ungapped BLOSUM62 alignment with raw score
    /// <paramref name="rawScore"/> over a search space of size m·n:
    /// E = K·m·n·e^(−λ·S). Source: Altschul et al., NCBI BLAST tutorial.
    /// </summary>
    /// <param name="rawScore">Raw alignment score S.</param>
    /// <param name="queryLength">Query length m.</param>
    /// <param name="subjectLength">Matched-subject length n.</param>
    public static double ExpectedValue(double rawScore, int queryLength, int subjectLength)
        => Blosum62UngappedK * queryLength * subjectLength
           * Math.Exp(-Blosum62UngappedLambda * rawScore);

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
    /// Pathway over-representation (enrichment) analysis using the hypergeometric test.
    /// For each pathway, computes the right-tail probability of observing at least as many
    /// query genes in the pathway as were observed, by random sampling from the background:
    /// P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M, n−i) / C(N, n), where N = background size,
    /// M = pathway size, n = query size, x = overlap. Source: PNNL Proteomics Data Analysis
    /// §8.2 Over-Representation Analysis.
    /// </summary>
    /// <param name="queryGenes">The gene set of interest (e.g. predicted/differential genes).</param>
    /// <param name="pathwayDatabase">Maps pathway id → its member genes.</param>
    /// <param name="backgroundGenes">
    /// The background gene universe. If null or empty, the union of all pathway members is
    /// used as the background.
    /// </param>
    /// <returns>One <see cref="PathwayEnrichment"/> per pathway, ascending by p-value.</returns>
    /// <exception cref="ArgumentNullException">If <paramref name="queryGenes"/> or <paramref name="pathwayDatabase"/> is null.</exception>
    public static IReadOnlyList<PathwayEnrichment> FindPathwayEnrichment(
        IEnumerable<string> queryGenes,
        IReadOnlyDictionary<string, IReadOnlyCollection<string>> pathwayDatabase,
        IEnumerable<string>? backgroundGenes = null)
    {
        if (queryGenes is null) throw new ArgumentNullException(nameof(queryGenes));
        if (pathwayDatabase is null) throw new ArgumentNullException(nameof(pathwayDatabase));

        var query = new HashSet<string>(queryGenes, StringComparer.Ordinal);

        var background = backgroundGenes is not null
            ? new HashSet<string>(backgroundGenes, StringComparer.Ordinal)
            : new HashSet<string>(StringComparer.Ordinal);
        if (background.Count == 0)
        {
            foreach (var members in pathwayDatabase.Values)
                foreach (var gene in members)
                    background.Add(gene);
        }
        // The query is part of the universe being sampled from.
        background.UnionWith(query);

        int n = query.Count;            // sample size (interesting genes)
        int bigN = background.Count;    // population size (background)

        var results = new List<PathwayEnrichment>(pathwayDatabase.Count);
        foreach (var (pathwayId, rawMembers) in pathwayDatabase)
        {
            // Count pathway members and overlap relative to the background universe.
            var members = new HashSet<string>(rawMembers, StringComparer.Ordinal);
            members.IntersectWith(background);

            int bigM = members.Count;   // successes in population (pathway size)
            int x = 0;                  // observed successes in sample (overlap)
            foreach (var gene in query)
                if (members.Contains(gene)) x++;

            double pValue = HypergeometricUpperTail(x, bigN, bigM, n);
            results.Add(new PathwayEnrichment(pathwayId, x, bigM, n, bigN, pValue));
        }

        results.Sort((a, b) => a.PValue.CompareTo(b.PValue));
        return results;
    }

    /// <summary>
    /// Right-tail hypergeometric probability P(X ≥ x) for drawing <paramref name="n"/> items
    /// from a population of <paramref name="bigN"/> containing <paramref name="bigM"/> successes:
    /// P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M, n−i) / C(N, n).
    /// Computed in log-space (log-Gamma) for numerical stability. Source: PNNL ORA §8.2.
    /// </summary>
    public static double HypergeometricUpperTail(int x, int bigN, int bigM, int n)
    {
        // Degenerate / empty-sum cases: no over-representation possible ⇒ p = 1.
        if (x <= 0 || bigM <= 0 || n <= 0 || bigN <= 0)
            return 1.0;

        double logDenom = LogChoose(bigN, n);
        double cumulative = 0.0; // Σ_{i=0}^{x-1} P(X = i)
        for (int i = 0; i < x; i++)
        {
            // P(X = i) is 0 when the partial table is infeasible (LogChoose → −∞).
            double logTerm = LogChoose(bigM, i) + LogChoose(bigN - bigM, n - i) - logDenom;
            if (!double.IsNegativeInfinity(logTerm))
                cumulative += Math.Exp(logTerm);
        }

        return Math.Clamp(1.0 - cumulative, 0.0, 1.0);
    }

    /// <summary>Log of the binomial coefficient C(n, k) via log-Gamma; −∞ when k ∉ [0, n].</summary>
    private static double LogChoose(int n, int k)
    {
        if (k < 0 || k > n)
            return double.NegativeInfinity;

        return LogGamma(n + 1) - LogGamma(k + 1) - LogGamma(n - k + 1);
    }

    /// <summary>Lanczos approximation of ln Γ(x) for x &gt; 0.</summary>
    private static double LogGamma(double x)
    {
        if (x <= 0)
            return double.PositiveInfinity;

        // Lanczos approximation coefficients.
        double[] c =
        {
            76.18009172947146,
            -86.50532032941677,
            24.01409824083091,
            -1.231739572450155,
            0.1208650973866179e-2,
            -0.5395239384953e-5,
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

    // ResFinder-style acquired-resistance-gene detection thresholds.
    // ResFinder reports the best-matching database gene found by a sequence search and
    // applies a percent-identity (%ID) cutoff and a coverage (length) cutoff.
    // Zankari et al. (2012) JAC 67(11):2640-2644: the web service default %ID is 100% and
    // is user-selectable; a gene must "cover at least 2/5 of the length of the resistance
    // gene in the database". The 90% identity / 60% coverage pair is the documented
    // ResFinder web-service operating point (Zankari et al. 2017, JAC 72(10):2764-2768).
    // Source URLs recorded in docs/Evidence/META-RESIST-001-Evidence.md.

    /// <summary>Default minimum percent identity (0–1) for calling a resistance gene present.</summary>
    public const double DefaultResistanceIdentityThreshold = 0.90;   // ResFinder web service default %ID.

    /// <summary>Default minimum coverage (0–1) of the reference gene length.</summary>
    public const double DefaultResistanceCoverageThreshold = 0.60;   // ResFinder min coverage (Zankari 2012/2017).

    /// <summary>
    /// A detected antibiotic-resistance gene: the best-matching reference gene for a contig
    /// that meets the identity and coverage thresholds.
    /// </summary>
    /// <param name="ContigId">Identifier of the query contig/gene that was searched.</param>
    /// <param name="ResistanceGene">Name of the matched reference resistance gene.</param>
    /// <param name="AntibioticClass">Antibiotic class conferred by the matched gene.</param>
    /// <param name="PercentIdentity">
    /// BLAST-style percent identity (0–1): identical positions divided by the gapless
    /// alignment length of the best ungapped match.
    /// </param>
    /// <param name="Coverage">
    /// Fraction (0–1) of the reference gene length spanned by the best ungapped match.
    /// </param>
    public readonly record struct ResistanceHit(
        string ContigId,
        string ResistanceGene,
        string AntibioticClass,
        double PercentIdentity,
        double Coverage);

    /// <summary>
    /// Detects acquired antibiotic-resistance genes in assembled contigs against a
    /// caller-supplied reference database, following the ResFinder methodology: for each
    /// reference gene the best ungapped alignment to the contig is located, percent identity
    /// and reference coverage are computed, and the reference gene is reported only when it
    /// passes both thresholds. For each contig the single best-matching reference gene is
    /// reported (highest identity, ties broken by higher coverage), mirroring ResFinder's
    /// "best-matching gene" output.
    /// </summary>
    /// <param name="contigs">Assembled contigs (id + nucleotide sequence) to screen.</param>
    /// <param name="referenceGenes">
    /// Reference resistance genes (id, full nucleotide sequence, gene name, antibiotic class).
    /// The caller supplies the curated database (e.g. ResFinder/CARD); no gene list is hard-coded.
    /// </param>
    /// <param name="identityThreshold">Minimum percent identity (0–1). Defaults to the ResFinder web-service value.</param>
    /// <param name="coverageThreshold">Minimum reference coverage (0–1). Defaults to the ResFinder value.</param>
    /// <returns>One <see cref="ResistanceHit"/> per contig that has a passing best match.</returns>
    /// <exception cref="ArgumentNullException">If <paramref name="contigs"/> or <paramref name="referenceGenes"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">If a threshold is outside [0, 1].</exception>
    public static IEnumerable<ResistanceHit> FindAntibioticResistanceGenes(
        IEnumerable<(string ContigId, string Sequence)> contigs,
        IEnumerable<(string GeneId, string Sequence, string Name, string AntibioticClass)> referenceGenes,
        double identityThreshold = DefaultResistanceIdentityThreshold,
        double coverageThreshold = DefaultResistanceCoverageThreshold)
    {
        if (contigs is null) throw new ArgumentNullException(nameof(contigs));
        if (referenceGenes is null) throw new ArgumentNullException(nameof(referenceGenes));
        if (identityThreshold < 0.0 || identityThreshold > 1.0)
            throw new ArgumentOutOfRangeException(nameof(identityThreshold), identityThreshold, "Identity threshold must be in [0, 1].");
        if (coverageThreshold < 0.0 || coverageThreshold > 1.0)
            throw new ArgumentOutOfRangeException(nameof(coverageThreshold), coverageThreshold, "Coverage threshold must be in [0, 1].");

        var refList = referenceGenes
            .Where(r => !string.IsNullOrEmpty(r.Sequence))
            .ToList();

        return Iterate();

        IEnumerable<ResistanceHit> Iterate()
        {
            foreach (var (contigId, sequence) in contigs)
            {
                if (string.IsNullOrEmpty(sequence))
                    continue;

                ResistanceHit? best = null;

                foreach (var reference in refList)
                {
                    var (identity, coverage) = BestUngappedMatch(sequence, reference.Sequence);

                    // ResFinder: report only matches passing both the %ID and coverage cutoffs.
                    if (identity < identityThreshold || coverage < coverageThreshold)
                        continue;

                    var candidate = new ResistanceHit(
                        ContigId: contigId,
                        ResistanceGene: reference.Name,
                        AntibioticClass: reference.AntibioticClass,
                        PercentIdentity: identity,
                        Coverage: coverage);

                    // Best-matching gene: highest identity, ties broken by greater coverage.
                    if (best is null
                        || candidate.PercentIdentity > best.Value.PercentIdentity
                        || (candidate.PercentIdentity == best.Value.PercentIdentity
                            && candidate.Coverage > best.Value.Coverage))
                    {
                        best = candidate;
                    }
                }

                if (best is not null)
                    yield return best.Value;
            }
        }
    }

    /// <summary>
    /// Finds the best ungapped (gapless) alignment of <paramref name="reference"/> within
    /// <paramref name="contig"/> and returns its BLAST-style percent identity and reference
    /// coverage. The reference is slid across the contig at every offset (including offsets
    /// where it overhangs an end, so partially assembled genes on a contig edge are scored);
    /// for each offset identical positions are counted over the overlapping window. The
    /// returned identity = identical positions / window length (gaps are not introduced, so
    /// the gapless alignment length equals the overlap length); coverage = window length /
    /// reference length. The offset maximizing identical positions is chosen.
    /// </summary>
    private static (double Identity, double Coverage) BestUngappedMatch(string contig, string reference)
    {
        int n = contig.Length;
        int m = reference.Length;
        int bestMatches = 0;
        int bestWindow = 0;

        // Offsets from -(m-1) (reference overhangs the left end) to (n-1) (overhangs the right end).
        for (int offset = -(m - 1); offset <= n - 1; offset++)
        {
            int start = Math.Max(0, offset);
            int end = Math.Min(n, offset + m);     // exclusive on the contig
            int window = end - start;
            if (window <= 0)
                continue;

            int matches = 0;
            for (int i = 0; i < window; i++)
            {
                // contig[start + i] aligns to reference[(start - offset) + i]
                if (contig[start + i] == reference[start - offset + i])
                    matches++;
            }

            // Prefer the window with the most identical positions; ties favour the longer window.
            if (matches > bestMatches || (matches == bestMatches && window > bestWindow))
            {
                bestMatches = matches;
                bestWindow = window;
            }
        }

        if (bestWindow == 0)
            return (0.0, 0.0);

        double identity = (double)bestMatches / bestWindow;   // identical positions / gapless alignment length
        double coverage = (double)bestWindow / m;             // fraction of reference gene length covered
        return (identity, coverage);
    }

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
