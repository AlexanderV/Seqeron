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
    /// Represents the taxonomic classification of a single read by the Kraken-style
    /// k-mer / lowest-common-ancestor classifier.
    /// </summary>
    /// <param name="ReadId">Input read identifier.</param>
    /// <param name="TaxonId">
    /// Assigned taxon id — the leaf of the maximum-scoring root-to-leaf (RTL) path in the read's
    /// classification tree (Wood &amp; Salzberg 2014). Reads with no k-mer hits are assigned the
    /// root / unclassified taxon (<see cref="TaxonomyTree.RootId"/>).
    /// </param>
    /// <param name="TaxonName">Name of the assigned taxon, from the taxonomy tree.</param>
    /// <param name="Rank">Rank label of the assigned taxon, from the taxonomy tree.</param>
    /// <param name="RtlScore">
    /// Score of the winning RTL path: the sum of node weights (k-mer hit counts) along the path
    /// from the root to the assigned leaf.
    /// </param>
    /// <param name="Confidence">
    /// Kraken C/Q confidence: C = number of k-mers mapped to a taxon in the clade rooted at the
    /// assigned label, Q = number of non-ambiguous k-mers queried. 0 when there are no hits.
    /// </param>
    /// <param name="MatchedKmers">C — k-mers supporting the assigned clade.</param>
    /// <param name="TotalKmers">Q — non-ambiguous k-mers queried against the database.</param>
    /// <param name="Kingdom">Kingdom on the assigned taxon's lineage, or empty.</param>
    /// <param name="Phylum">Phylum on the assigned taxon's lineage, or empty.</param>
    /// <param name="Class">Class on the assigned taxon's lineage, or empty.</param>
    /// <param name="Order">Order on the assigned taxon's lineage, or empty.</param>
    /// <param name="Family">Family on the assigned taxon's lineage, or empty.</param>
    /// <param name="Genus">Genus on the assigned taxon's lineage, or empty.</param>
    /// <param name="Species">Species on the assigned taxon's lineage, or empty.</param>
    public readonly record struct TaxonomicClassification(
        string ReadId,
        int TaxonId,
        string TaxonName,
        string Rank,
        int RtlScore,
        double Confidence,
        int MatchedKmers,
        int TotalKmers,
        string Kingdom,
        string Phylum,
        string Class,
        string Order,
        string Family,
        string Genus,
        string Species);

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

    #region K-mer / LCA Taxonomic Classification (Kraken)

    // Kraken (Wood & Salzberg 2014, Genome Biology 15:R46) taxonomic classifier.
    //
    // Database build: "a database that contains records consisting of a k-mer and the LCA of all
    // organisms whose genomes contain that k-mer". As reference sequences are processed, "if a
    // k-mer from a sequence has had its LCA value previously set, then the LCA of the stored value
    // and the current sequence's taxon is calculated" and stored.
    //
    // Per-read classification: the read's hit taxa "and their ancestors in the taxonomy tree form
    // what we term the classification tree ... Each node in the classification tree is weighted with
    // the number of k-mers in K(S) that mapped to the taxon associated with that node. ... each
    // root-to-leaf (RTL) path in the classification tree is scored by calculating the sum of all
    // node weights along the path. The maximum scoring RTL path in the classification tree is the
    // classification path [and its leaf is the assigned label]." Tie-break: "if there are multiple
    // maximally scoring paths, the LCA of all those paths' leaves is selected."
    //
    // Unclassified: "Sequences for which none of the k-mers in K(S) are found in any genome are left
    // unclassified by this algorithm" — reported here as the root taxon (TaxonomyTree.RootId).
    //
    // Confidence reuses Kraken 2's C/Q score: C = k-mers mapped to a taxon in the clade rooted at the
    // assigned label, Q = non-ambiguous k-mers queried.

    private const string AcgtAlphabet = "ACGT";

    /// <summary>
    /// Classifies metagenomic reads with the Kraken k-mer / LCA algorithm: for each read it
    /// collects the (canonical) k-mer hits, builds the classification tree over the hit taxa and
    /// their ancestors weighted by k-mer count, finds the maximum-scoring root-to-leaf (RTL) path,
    /// and assigns the leaf of that path (LCA-of-leaves on ties). Reads with no hits are assigned
    /// the root / unclassified taxon.
    /// </summary>
    /// <param name="reads">Reads to classify (id + nucleotide sequence).</param>
    /// <param name="kmerDatabase">
    /// Canonical-k-mer → taxon-id database (see <see cref="BuildKmerDatabase"/>). Keys are canonical
    /// k-mers; values are taxon ids present in <paramref name="taxonomy"/>.
    /// </param>
    /// <param name="taxonomy">Taxonomy tree providing parent chains and the LCA operation.</param>
    /// <param name="k">K-mer length (default 31, per Kraken). Must be positive.</param>
    /// <returns>One <see cref="TaxonomicClassification"/> per input read, in input order.</returns>
    /// <exception cref="ArgumentNullException">A required argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> is not positive.</exception>
    public static IEnumerable<TaxonomicClassification> ClassifyReads(
        IEnumerable<(string Id, string Sequence)> reads,
        IReadOnlyDictionary<string, int> kmerDatabase,
        TaxonomyTree taxonomy,
        int k = 31)
    {
        ArgumentNullException.ThrowIfNull(reads);
        ArgumentNullException.ThrowIfNull(kmerDatabase);
        ArgumentNullException.ThrowIfNull(taxonomy);
        if (k <= 0) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be positive.");

        return Iterate();

        IEnumerable<TaxonomicClassification> Iterate()
        {
            foreach (var (id, sequence) in reads)
                yield return ClassifyRead(id, sequence, kmerDatabase, taxonomy, k);
        }
    }

    private static TaxonomicClassification ClassifyRead(
        string id,
        string sequence,
        IReadOnlyDictionary<string, int> kmerDatabase,
        TaxonomyTree taxonomy,
        int k)
    {
        // Collect per-taxon k-mer hit counts (the leaf weights of the classification tree) and Q.
        var hitCounts = new Dictionary<int, int>();
        int totalKmers = 0; // Q = non-ambiguous k-mers queried

        if (!string.IsNullOrEmpty(sequence) && sequence.Length >= k)
        {
            string seq = sequence.ToUpperInvariant();
            for (int i = 0; i <= seq.Length - k; i++)
            {
                string kmer = seq.Substring(i, k);
                if (!IsAcgtOnly(kmer))
                    continue; // skip ambiguous k-mers (not counted in Q)

                totalKmers++;
                string canonical = GetCanonicalKmer(kmer);
                if (kmerDatabase.TryGetValue(canonical, out int taxon) && taxonomy.Contains(taxon))
                {
                    hitCounts.TryGetValue(taxon, out int c);
                    hitCounts[taxon] = c + 1;
                }
            }
        }

        if (hitCounts.Count == 0)
            return BuildResult(id, taxonomy.Root, 0, 0, totalKmers, taxonomy);

        // Classification tree node weights: weight[node] = k-mers mapped to that exact taxon
        // (= hitCounts). The RTL score of a leaf L is Σ weight over taxa lying on the path Root..L.
        //
        // Leaves of the classification tree = hit taxa that are not a proper ancestor of another hit
        // taxon. (A hit taxon that is an ancestor of another hit taxon is an internal node, never a
        // leaf, since the deeper hit extends the path.) Score each candidate leaf as the sum of
        // hit-weights over taxa lying on its root path.
        var hitTaxa = hitCounts.Keys.ToList();
        var leaves = new List<int>();
        foreach (int candidate in hitTaxa)
        {
            bool isAncestorOfAnother = hitTaxa.Any(other =>
                other != candidate && taxonomy.IsAncestorOf(candidate, other));
            if (!isAncestorOfAnother)
                leaves.Add(candidate);
        }

        int bestScore = int.MinValue;
        var bestLeaves = new List<int>();
        foreach (int leaf in leaves)
        {
            int score = 0;
            foreach (int node in taxonomy.GetPathToRoot(leaf))
                if (hitCounts.TryGetValue(node, out int w))
                    score += w;

            if (score > bestScore)
            {
                bestScore = score;
                bestLeaves.Clear();
                bestLeaves.Add(leaf);
            }
            else if (score == bestScore)
            {
                bestLeaves.Add(leaf);
            }
        }

        // Tie-break: LCA of all maximally-scoring leaves.
        int assigned = bestLeaves.Count == 1 ? bestLeaves[0] : taxonomy.Lca(bestLeaves);

        // C = k-mers mapped to any taxon in the clade rooted at the assigned label.
        int cladeKmers = 0;
        foreach (var (taxon, weight) in hitCounts)
            if (taxonomy.IsAncestorOf(assigned, taxon))
                cladeKmers += weight;

        return BuildResult(id, assigned, bestScore, cladeKmers, totalKmers, taxonomy);
    }

    private static TaxonomicClassification BuildResult(
        string id, int assignedTaxon, int rtlScore, int cladeKmers, int totalKmers, TaxonomyTree taxonomy)
    {
        var node = taxonomy.GetNode(assignedTaxon);
        double confidence = totalKmers > 0 ? (double)cladeKmers / totalKmers : 0.0;
        var ranks = ExtractRankLineage(assignedTaxon, taxonomy);

        return new TaxonomicClassification(
            ReadId: id,
            TaxonId: assignedTaxon,
            TaxonName: node.Name,
            Rank: node.Rank,
            RtlScore: rtlScore,
            Confidence: confidence,
            MatchedKmers: cladeKmers,
            TotalKmers: totalKmers,
            Kingdom: ranks.GetValueOrDefault("kingdom", ""),
            Phylum: ranks.GetValueOrDefault("phylum", ""),
            Class: ranks.GetValueOrDefault("class", ""),
            Order: ranks.GetValueOrDefault("order", ""),
            Family: ranks.GetValueOrDefault("family", ""),
            Genus: ranks.GetValueOrDefault("genus", ""),
            Species: ranks.GetValueOrDefault("species", ""));
    }

    /// <summary>
    /// Reads the seven standard ranks off the assigned taxon's lineage (root path), keyed by the
    /// lowercased rank label on each node, for compatibility with the downstream profile.
    /// </summary>
    private static Dictionary<string, string> ExtractRankLineage(int taxonId, TaxonomyTree taxonomy)
    {
        var result = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        foreach (int node in taxonomy.GetPathToRoot(taxonId))
        {
            var n = taxonomy.GetNode(node);
            string rank = (n.Rank ?? "").ToLowerInvariant();
            switch (rank)
            {
                case "kingdom":
                case "domain":
                    result["kingdom"] = n.Name; break;
                case "phylum":
                case "class":
                case "order":
                case "family":
                case "genus":
                case "species":
                    result[rank] = n.Name; break;
            }
        }
        return result;
    }

    /// <summary>
    /// Builds a Kraken-style canonical-k-mer → taxon database from labeled reference sequences.
    /// Each canonical k-mer is mapped to the <em>lowest common ancestor</em> of the taxa of all
    /// references that contain it: the first reference sets the k-mer's taxon, and every subsequent
    /// reference replaces the stored value with its LCA against that reference's taxon.
    /// </summary>
    /// <param name="referenceSequences">
    /// Labeled references (taxon id + nucleotide sequence). The taxon id must exist in
    /// <paramref name="taxonomy"/>.
    /// </param>
    /// <param name="taxonomy">Taxonomy tree used to compute the LCA of shared k-mers.</param>
    /// <param name="k">K-mer length (default 31). Must be positive.</param>
    /// <returns>A canonical-k-mer → taxon-id database for <see cref="ClassifyReads"/>.</returns>
    /// <exception cref="ArgumentNullException">A required argument is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="k"/> is not positive.</exception>
    /// <exception cref="KeyNotFoundException">A reference's taxon id is not in <paramref name="taxonomy"/>.</exception>
    public static Dictionary<string, int> BuildKmerDatabase(
        IEnumerable<(int TaxonId, string Sequence)> referenceSequences,
        TaxonomyTree taxonomy,
        int k = 31)
    {
        ArgumentNullException.ThrowIfNull(referenceSequences);
        ArgumentNullException.ThrowIfNull(taxonomy);
        if (k <= 0) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be positive.");

        var database = new Dictionary<string, int>();

        foreach (var (taxonId, sequence) in referenceSequences)
        {
            if (string.IsNullOrEmpty(sequence) || sequence.Length < k)
                continue;
            if (!taxonomy.Contains(taxonId))
                throw new KeyNotFoundException($"Reference taxon id {taxonId} is not in the taxonomy tree.");

            string seq = sequence.ToUpperInvariant();
            for (int i = 0; i <= seq.Length - k; i++)
            {
                string kmer = seq.Substring(i, k);
                if (!IsAcgtOnly(kmer))
                    continue;

                string canonical = GetCanonicalKmer(kmer);
                if (database.TryGetValue(canonical, out int existing))
                    database[canonical] = taxonomy.Lca(existing, taxonId); // collapse to LCA
                else
                    database[canonical] = taxonId;
            }
        }

        return database;
    }

    private static bool IsAcgtOnly(string kmer)
    {
        foreach (char c in kmer)
            if (AcgtAlphabet.IndexOf(c) < 0)
                return false;
        return true;
    }

    private static string GetCanonicalKmer(string kmer)
    {
        string revComp = DnaSequence.GetReverseComplementString(kmer);
        return string.Compare(kmer, revComp, StringComparison.Ordinal) <= 0 ? kmer : revComp;
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

    // Number of distinct tetranucleotides over the {A,C,G,T} alphabet (4^4),
    // the full TETRA signature dimensionality (Teeling et al. 2004).
    private const int TetranucleotideAlphabetSize = 256;

    /// <summary>
    /// Computes the TETRA z-score tetranucleotide signature of a DNA sequence
    /// (Teeling et al. 2004, BMC Bioinformatics 5:163; Environ Microbiol 6(9):938–947).
    /// <para>
    /// This is the <b>opt-in, Markov-corrected</b> signature: for each of the 256 tetranucleotides
    /// the observed count is compared to the value <b>predicted by a maximal-order (2nd-order)
    /// Markov model</b> from the constituent di-/trinucleotide composition, and the divergence is
    /// expressed as a z-score using the Schbath (1997) variance approximation:
    /// </para>
    /// <list type="bullet">
    /// <item><description>Expected count E(n1n2n3n4) = N(n1n2n3)·N(n2n3n4) / N(n2n3).</description></item>
    /// <item><description>Variance var(n1n2n3n4) = E·[N(n2n3)−N(n1n2n3)]·[N(n2n3)−N(n2n3n4)] / N(n2n3)².</description></item>
    /// <item><description>z(n1n2n3n4) = (N(n1n2n3n4) − E) / √var.</description></item>
    /// </list>
    /// <para>
    /// As in TETRA, the sequence is first extended by its reverse complement so the signature is
    /// strand-symmetric; counts are taken over overlapping words on the extended strand. A
    /// tetranucleotide whose denominator N(n2n3)=0 or whose variance ≤0 receives z=0 (no evidence
    /// of over-/under-representation). The 256-component z-score vectors of two sequences are then
    /// compared with <see cref="TetranucleotideZScoreCorrelation"/> (Pearson correlation).
    /// </para>
    /// <para>
    /// This does <b>not</b> replace the default raw-frequency path used internally by
    /// <see cref="BinContigs"/>; it is provided for callers who want the full TETRA z-score signature.
    /// </para>
    /// </summary>
    /// <param name="sequence">DNA sequence (case-insensitive; non-ACGT characters are skipped).</param>
    /// <returns>
    /// A 256-entry map keyed by the ACGT tetranucleotide, with the Teeling z-score for each.
    /// Because of the reverse-complement extension, an input of ≥2 usable ACGT bases already yields
    /// a ≥4-nt extended strand and a non-trivial signature; null, empty, or single-base input
    /// produces an all-zero 256-entry map.
    /// </returns>
    public static IReadOnlyDictionary<string, double> CalculateTetranucleotideZScores(string sequence)
    {
        var z = new Dictionary<string, double>(TetranucleotideAlphabetSize);

        // TETRA extends the sequence by its reverse complement to make the signature
        // strand-symmetric (compensates leading/lagging-strand tetranucleotide skew).
        string extended = ExtendWithReverseComplement(sequence);

        var n4 = new Dictionary<string, int>();
        var n3 = new Dictionary<string, int>();
        var n2 = new Dictionary<string, int>();
        CountOligonucleotides(extended, n4, n3, n2);

        foreach (string tetra in EnumerateTetranucleotides())
        {
            z[tetra] = TetranucleotideZScore(tetra, n4, n3, n2);
        }

        return z;
    }

    /// <summary>
    /// Pearson correlation of the TETRA z-score signatures of two DNA sequences
    /// (Teeling et al. 2004) — the tetranucleotide correlation coefficient.
    /// Identical signatures correlate to 1.0; unrelated compositions correlate near 0.
    /// </summary>
    /// <param name="sequenceA">First DNA sequence.</param>
    /// <param name="sequenceB">Second DNA sequence.</param>
    /// <returns>
    /// Pearson r ∈ [−1, 1] over the aligned 256-component z-score vectors; 0 when either vector
    /// has zero variance (no usable signal).
    /// </returns>
    public static double TetranucleotideZScoreCorrelation(string sequenceA, string sequenceB)
    {
        var za = CalculateTetranucleotideZScores(sequenceA);
        var zb = CalculateTetranucleotideZScores(sequenceB);

        // Both maps share the same fixed 256-key ordering (EnumerateTetranucleotides),
        // so iterating one key set aligns the vectors component-wise.
        double[] va = za.Values.ToArray();
        double[] vb = EnumerateTetranucleotides().Select(t => zb[t]).ToArray();

        return PearsonCorrelation(va, vb);
    }

    /// <summary>
    /// Teeling/Schbath z-score for one tetranucleotide from the maximal-order (2nd-order) Markov
    /// expected count and its variance approximation. Returns 0 when the denominator N(n2n3)=0
    /// or the variance is non-positive (no evidence of over-/under-representation).
    /// </summary>
    private static double TetranucleotideZScore(
        string tetra,
        Dictionary<string, int> n4,
        Dictionary<string, int> n3,
        Dictionary<string, int> n2)
    {
        int observed = n4.GetValueOrDefault(tetra, 0);
        int prefix3 = n3.GetValueOrDefault(tetra.Substring(0, 3), 0);   // N(n1n2n3)
        int suffix3 = n3.GetValueOrDefault(tetra.Substring(1, 3), 0);   // N(n2n3n4)
        int middle2 = n2.GetValueOrDefault(tetra.Substring(1, 2), 0);   // N(n2n3)

        if (middle2 == 0)
            return 0.0;

        // Expected count: E(n1n2n3n4) = N(n1n2n3)·N(n2n3n4) / N(n2n3)   (Teeling 2004)
        double expected = (double)prefix3 * suffix3 / middle2;

        // Variance (Schbath 1997 approximation):
        // var = E·[N(n2n3)−N(n1n2n3)]·[N(n2n3)−N(n2n3n4)] / N(n2n3)²
        double variance = expected
            * (middle2 - prefix3) * (double)(middle2 - suffix3)
            / ((double)middle2 * middle2);

        if (variance <= 0)
            return 0.0;

        // z(n1n2n3n4) = (N(n1n2n3n4) − E) / √var
        return (observed - expected) / Math.Sqrt(variance);
    }

    /// <summary>
    /// Appends the reverse complement of the ACGT-filtered, upper-cased sequence (TETRA
    /// strand-symmetric extension). Non-ACGT characters are dropped before extension.
    /// </summary>
    private static string ExtendWithReverseComplement(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return string.Empty;

        var forward = new StringBuilder(sequence.Length);
        foreach (char c in sequence)
        {
            char u = char.ToUpperInvariant(c);
            if (u is 'A' or 'C' or 'G' or 'T')
                forward.Append(u);
        }

        var rc = new StringBuilder(forward.Length);
        for (int i = forward.Length - 1; i >= 0; i--)
        {
            rc.Append(forward[i] switch
            {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => 'N'
            });
        }

        return forward.Append(rc).ToString();
    }

    /// <summary>
    /// Counts overlapping di-, tri-, and tetranucleotide words over an ACGT-only string in a
    /// single pass. Words containing any non-ACGT character are skipped.
    /// </summary>
    private static void CountOligonucleotides(
        string s,
        Dictionary<string, int> n4,
        Dictionary<string, int> n3,
        Dictionary<string, int> n2)
    {
        for (int i = 0; i + 2 <= s.Length; i++)
        {
            if (i + 2 <= s.Length) Increment(n2, s.Substring(i, 2));
            if (i + 3 <= s.Length) Increment(n3, s.Substring(i, 3));
            if (i + 4 <= s.Length) Increment(n4, s.Substring(i, 4));
        }
    }

    private static void Increment(Dictionary<string, int> map, string key)
    {
        if (!key.All(c => c is 'A' or 'C' or 'G' or 'T'))
            return;
        map[key] = map.GetValueOrDefault(key, 0) + 1;
    }

    private static readonly char[] DnaBases = { 'A', 'C', 'G', 'T' };

    /// <summary>Enumerates all 256 ACGT tetranucleotides in a fixed lexicographic order.</summary>
    private static IEnumerable<string> EnumerateTetranucleotides()
    {
        foreach (char a in DnaBases)
            foreach (char b in DnaBases)
                foreach (char c in DnaBases)
                    foreach (char d in DnaBases)
                        yield return new string(new[] { a, b, c, d });
    }

    /// <summary>Pearson product-moment correlation of two equal-length vectors; 0 if either is constant.</summary>
    private static double PearsonCorrelation(double[] a, double[] b)
    {
        if (a.Length == 0 || a.Length != b.Length)
            return 0.0;

        double meanA = a.Average(), meanB = b.Average();
        double sumAB = 0, sumA2 = 0, sumB2 = 0;
        for (int i = 0; i < a.Length; i++)
        {
            double da = a[i] - meanA, db = b[i] - meanB;
            sumAB += da * db;
            sumA2 += da * da;
            sumB2 += db * db;
        }

        double denom = Math.Sqrt(sumA2) * Math.Sqrt(sumB2);
        return denom > 0 ? sumAB / denom : 0.0;
    }

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
    /// P(X ≥ x) = Σ_{i=x}^{min(n,M)} C(M,i)·C(N−M, n−i) / C(N, n)
    ///          = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M, n−i) / C(N, n).
    /// The upper tail is summed directly (i = x … min(n, M)) rather than as 1 − lower-tail, to
    /// avoid catastrophic cancellation when the tail probability is tiny (e.g. 7.88e-8 for the
    /// PNNL §8.2 example): subtracting a sum ≈ 1 from 1.0 would destroy the small result.
    /// Each PMF term is computed in log-space (log-Gamma) for numerical stability.
    /// Source: PNNL ORA §8.2; Boyle et al. (2004), Bioinformatics 20(18):3710-3715.
    /// </summary>
    public static double HypergeometricUpperTail(int x, int bigN, int bigM, int n)
    {
        // Degenerate / empty-sum cases: no over-representation possible ⇒ p = 1.
        if (x <= 0 || bigM <= 0 || n <= 0 || bigN <= 0)
            return 1.0;

        double logDenom = LogChoose(bigN, n);
        double tail = 0.0; // Σ_{i=x}^{min(n,M)} P(X = i)
        int upper = Math.Min(n, bigM); // beyond this, C(M,i) = 0 (cannot draw more successes than exist)
        for (int i = x; i <= upper; i++)
        {
            // P(X = i) is 0 when the partial table is infeasible (LogChoose → −∞).
            double logTerm = LogChoose(bigM, i) + LogChoose(bigN - bigM, n - i) - logDenom;
            if (!double.IsNegativeInfinity(logTerm))
                tail += Math.Exp(logTerm);
        }

        return Math.Clamp(tail, 0.0, 1.0);
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
    /// reference length. The offset maximizing identical positions is chosen; on an equal
    /// match count the higher-identity (shorter) window is preferred so the alignment is never
    /// padded with flanking mismatches (a padded window would lower identity and could spuriously
    /// fail the identity threshold even when a perfect HSP exists).
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

            // Prefer the window with the most identical positions; on an equal match count
            // favour the shorter (higher-identity) window so the chosen alignment is not padded
            // with mismatching flanks. bestWindow == 0 marks the unset state.
            if (matches > bestMatches
                || (matches == bestMatches && bestWindow != 0 && window < bestWindow))
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

    #region Significant Taxa Detection (Mann–Whitney U)

    // Normal-approximation constants for the Mann–Whitney U statistic.
    // Mean m_U = n1*n2/2 and variance n1*n2*(n1+n2+1)/12 per Mann & Whitney (1947),
    // Ann. Math. Statist. 18(1):50–60 (see docs/Evidence/META-TAXA-001-Evidence.md).
    private const double MannWhitneyVarianceDivisor = 12.0;
    // Continuity correction subtracts 0.5 from |U − m_U| (SciPy mannwhitneyu, use_continuity=True default).
    private const double ContinuityCorrection = 0.5;

    /// <summary>
    /// Result of a single taxon's two-group differential-abundance test.
    /// </summary>
    /// <param name="Taxon">Taxon identifier.</param>
    /// <param name="U">Mann–Whitney U statistic (the larger of U1, U2 used for the z-score).</param>
    /// <param name="Z">Normal-approximation z-score, <c>(|U − m_U| − cc) / σ_U</c>.</param>
    /// <param name="PValue">Two-tailed asymptotic p-value, <c>2·SF(|z|)</c>, clamped to [0,1].</param>
    /// <param name="Significant"><c>PValue &lt; pThreshold</c>.</param>
    public readonly record struct SignificantTaxon(
        string Taxon,
        double U,
        double Z,
        double PValue,
        bool Significant);

    /// <summary>
    /// Result of the Mann–Whitney U (Wilcoxon rank-sum) test under the normal approximation.
    /// </summary>
    /// <param name="U1">U statistic for <c>group1</c>: <c>R1 − n1(n1+1)/2</c>.</param>
    /// <param name="U2">U statistic for <c>group2</c>: <c>n1·n2 − U1</c>.</param>
    /// <param name="Z">z-score from the larger U: <c>(|U − m_U| − cc) / σ_U</c>.</param>
    /// <param name="PValue">Two-tailed asymptotic p-value <c>2·SF(|z|)</c>, clamped to [0,1].</param>
    public readonly record struct MannWhitneyResult(double U1, double U2, double Z, double PValue);

    /// <summary>
    /// Computes the Mann–Whitney U (Wilcoxon rank-sum) test between two independent samples using
    /// the asymptotic normal approximation with midrank tie handling.
    /// </summary>
    /// <remarks>
    /// U1 = R1 − n1(n1+1)/2 where R1 is the rank sum of <paramref name="group1"/> in the pooled
    /// midrank ranking; U2 = n1·n2 − U1. The z-score uses m_U = n1·n2/2 and the tie-corrected
    /// σ_U = sqrt(n1·n2·(n1+n2+1)/12 − n1·n2·Σ(t³−t)/(12·n·(n−1))). Mann &amp; Whitney (1947).
    /// </remarks>
    /// <param name="group1">First sample's observations (e.g., abundances). Must be non-empty.</param>
    /// <param name="group2">Second sample's observations. Must be non-empty.</param>
    /// <param name="useContinuityCorrection">Subtract 0.5 from |U − m_U| (default true, matching SciPy).</param>
    /// <returns>U1, U2, z, and the two-tailed asymptotic p-value.</returns>
    /// <exception cref="ArgumentNullException">A group is null.</exception>
    /// <exception cref="ArgumentException">A group is empty.</exception>
    public static MannWhitneyResult MannWhitneyU(
        IReadOnlyList<double> group1,
        IReadOnlyList<double> group2,
        bool useContinuityCorrection = true)
    {
        ArgumentNullException.ThrowIfNull(group1);
        ArgumentNullException.ThrowIfNull(group2);
        if (group1.Count == 0 || group2.Count == 0)
            throw new ArgumentException("Both groups must contain at least one observation.");

        int n1 = group1.Count;
        int n2 = group2.Count;
        int n = n1 + n2;

        // Pool both samples and assign midranks (ties share the average of their positions).
        var pooled = new (double Value, int Group)[n];
        for (int i = 0; i < n1; i++) pooled[i] = (group1[i], 1);
        for (int j = 0; j < n2; j++) pooled[n1 + j] = (group2[j], 2);
        Array.Sort(pooled, (a, b) => a.Value.CompareTo(b.Value));

        var ranks = new double[n];
        double tieTermSum = 0; // Σ (t_k³ − t_k) over tie groups
        int idx = 0;
        while (idx < n)
        {
            int start = idx;
            while (idx < n && pooled[idx].Value == pooled[start].Value) idx++;
            int tieCount = idx - start;
            // Ranks are 1-based; midrank = average of the (start+1 .. idx) positions.
            double midRank = (start + 1 + idx) / 2.0;
            for (int k = start; k < idx; k++) ranks[k] = midRank;
            if (tieCount > 1)
            {
                double t = tieCount;
                tieTermSum += t * t * t - t;
            }
        }

        double r1 = 0;
        for (int k = 0; k < n; k++)
            if (pooled[k].Group == 1) r1 += ranks[k];

        // U1 = R1 − n1(n1+1)/2 ; U2 = n1·n2 − U1  (Mann & Whitney 1947).
        double u1 = r1 - (double)n1 * (n1 + 1) / 2.0;
        double nProduct = (double)n1 * n2;
        double u2 = nProduct - u1;

        double meanU = nProduct / 2.0;
        // Tie-corrected variance: n1·n2/12 · [(n+1) − Σ(t³−t)/(n(n−1))].
        double variance = n > 1
            ? nProduct / MannWhitneyVarianceDivisor * ((n + 1) - tieTermSum / ((double)n * (n - 1)))
            : 0.0;

        double z;
        double pValue;
        if (variance <= 0)
        {
            // Degenerate: all observations tied → σ = 0, no evidence against H0.
            z = 0.0;
            pValue = 1.0;
        }
        else
        {
            double sigma = Math.Sqrt(variance);
            double uForZ = Math.Max(u1, u2);
            double distance = Math.Abs(uForZ - meanU);
            if (useContinuityCorrection)
                distance = Math.Max(0.0, distance - ContinuityCorrection);
            z = distance / sigma;
            // Two-tailed: 2·SF(|z|) = 2·(1 − Φ(|z|)). Clamp to [0,1] for tiny numerical overshoot.
            pValue = Math.Clamp(2.0 * (1.0 - StatisticsHelper.NormalCDF(z)), 0.0, 1.0);
        }

        return new MannWhitneyResult(u1, u2, z, pValue);
    }

    /// <summary>
    /// Identifies taxa whose abundances differ significantly between two sample groups using a
    /// per-taxon Mann–Whitney U (Wilcoxon rank-sum) test, the standard non-parametric approach for
    /// differential abundance in metagenomics (Xia &amp; Sun 2017).
    /// </summary>
    /// <param name="profiles">Per-sample taxon→abundance maps. A taxon absent in a profile counts as abundance 0.</param>
    /// <param name="groups">Group label (1 or 2) for each profile, aligned by index with <paramref name="profiles"/>.</param>
    /// <param name="pThreshold">Significance threshold; a taxon is significant when its p-value is below it.</param>
    /// <param name="useContinuityCorrection">Passed to <see cref="MannWhitneyU"/> (default true).</param>
    /// <returns>One <see cref="SignificantTaxon"/> per taxon observed in any profile, ordered by ascending p-value.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="profiles"/> or <paramref name="groups"/> is null.</exception>
    /// <exception cref="ArgumentException">Counts mismatch, fewer than two groups present, or a profile lacks a 1/2 label.</exception>
    public static IReadOnlyList<SignificantTaxon> FindSignificantTaxa(
        IReadOnlyList<IReadOnlyDictionary<string, double>> profiles,
        IReadOnlyList<int> groups,
        double pThreshold = 0.05,
        bool useContinuityCorrection = true)
    {
        ArgumentNullException.ThrowIfNull(profiles);
        ArgumentNullException.ThrowIfNull(groups);
        if (profiles.Count != groups.Count)
            throw new ArgumentException("profiles and groups must have the same length.");
        if (profiles.Count == 0)
            return Array.Empty<SignificantTaxon>();

        var group1Indices = new List<int>();
        var group2Indices = new List<int>();
        for (int i = 0; i < groups.Count; i++)
        {
            if (groups[i] == 1) group1Indices.Add(i);
            else if (groups[i] == 2) group2Indices.Add(i);
            else throw new ArgumentException($"Group label at index {i} must be 1 or 2.");
        }

        if (group1Indices.Count == 0 || group2Indices.Count == 0)
            throw new ArgumentException("Both group 1 and group 2 must contain at least one profile.");

        var allTaxa = new SortedSet<string>(StringComparer.Ordinal);
        foreach (var profile in profiles)
            foreach (var taxon in profile.Keys)
                allTaxa.Add(taxon);

        var results = new List<SignificantTaxon>(allTaxa.Count);
        foreach (var taxon in allTaxa)
        {
            var g1 = group1Indices.Select(i => profiles[i].GetValueOrDefault(taxon, 0.0)).ToList();
            var g2 = group2Indices.Select(i => profiles[i].GetValueOrDefault(taxon, 0.0)).ToList();

            var mw = MannWhitneyU(g1, g2, useContinuityCorrection);
            double u = Math.Max(mw.U1, mw.U2);
            results.Add(new SignificantTaxon(taxon, u, mw.Z, mw.PValue, mw.PValue < pThreshold));
        }

        // Deterministic ordering: ascending p-value, then taxon name for stable ties.
        results.Sort((a, b) =>
        {
            int cmp = a.PValue.CompareTo(b.PValue);
            return cmp != 0 ? cmp : string.CompareOrdinal(a.Taxon, b.Taxon);
        });
        return results;
    }

    #endregion

    #region CheckM-style marker-gene completeness / contamination (opt-in; META-BIN-001)

    /// <summary>
    /// A collocated single-copy marker set <c>s</c> — one of the marker sets <c>M</c> over which
    /// CheckM averages completeness and contamination (Parks et al. 2015). The marker IDs are
    /// opaque identifiers (e.g. Pfam accessions) keyed against a per-bin marker-count map.
    /// </summary>
    /// <param name="MarkerIds">The marker identifiers in this collocated set (set <c>s</c>).</param>
    public readonly record struct MarkerSet(IReadOnlyList<string> MarkerIds);

    /// <summary>
    /// CheckM-style marker-gene quality of a bin: completeness and contamination (percentages),
    /// computed from collocated single-copy marker sets per Parks et al. (2015).
    /// </summary>
    /// <param name="Completeness">Completeness in % (0–100).</param>
    /// <param name="Contamination">Contamination in % (≥ 0; can exceed 100 with heavy duplication).</param>
    /// <param name="MarkerSetCount">Number of collocated marker sets |M| used.</param>
    /// <param name="MarkerCount">Total number of markers across all sets.</param>
    /// <param name="MarkersPresent">Distinct markers found at least once.</param>
    public readonly record struct BinMarkerQuality(
        double Completeness,
        double Contamination,
        int MarkerSetCount,
        int MarkerCount,
        int MarkersPresent);

    /// <summary>
    /// A loaded single-copy marker HMM: a Plan7 profile plus its marker identifier and the
    /// per-sequence bit-score threshold used to call the marker present.
    /// </summary>
    /// <param name="MarkerId">Marker identifier (e.g. Pfam accession "PF00318").</param>
    /// <param name="Hmm">The parsed Plan7 profile HMM.</param>
    /// <param name="BitScoreThreshold">Bit-score gate for a hit (the Pfam GA1 gathering threshold).</param>
    public readonly record struct MarkerHmm(string MarkerId, Plan7ProfileHmm Hmm, double BitScoreThreshold);

    // ln 2, for converting Plan7 natural-log (nat) log-odds scores to bits.
    // Source: Eddy (2011) PLoS Comput Biol 7:e1002195 — HMMER reports log-odds in bits = nats / ln 2.
    private const double LogOddsBitsPerNat = 0.69314718055994530941723212145818;

    // Fallback per-sequence bit-score gate when a profile has no GA1 gathering threshold.
    // Conservative default mirroring ProteinMotifFinder.FindDomainsByHmm.
    private const double DefaultMarkerBitScoreThreshold = 10.0;

    // The bundled CC0 universal single-copy ribosomal-protein marker HMMs (see Resources/README.md).
    // Each is treated as its own singleton collocated marker set in the default set (ribosomal
    // proteins of distinct families are not collocated as one operon-set here — see algorithm doc).
    // Provenance: EMBL-EBI InterPro Pfam HMM API, 2026-06-25; licence CC0 (public domain).
    private static readonly (string Resource, string MarkerId)[] BundledRibosomalMarkers =
    {
        ("PF00318_Ribosomal_S2.hmm",  "PF00318"),
        ("PF00177_Ribosomal_S7.hmm",  "PF00177"),
        ("PF00410_Ribosomal_S8.hmm",  "PF00410"),
        ("PF00380_Ribosomal_S9.hmm",  "PF00380"),
        ("PF00338_Ribosomal_S10.hmm", "PF00338"),
        ("PF00411_Ribosomal_S11.hmm", "PF00411"),
        ("PF00203_Ribosomal_S19.hmm", "PF00203"),
        ("PF00687_Ribosomal_L1.hmm",  "PF00687"),
        ("PF00297_Ribosomal_L3.hmm",  "PF00297"),
    };

    // The Pfam-defined members of the GTDB bac120 (Bacteria) domain-level universal single-copy
    // marker set (Parks et al. 2018 Nat Biotechnol 36:996; GTDB-Tk bac120.marker_info). bac120 is
    // 120 markers (6 Pfam + 114 TIGRFAM); only the 6 Pfam markers are CC0 and bundled. The 114
    // TIGRFAM-defined markers (CC BY-SA 4.0) are NOT bundled — supply them via LoadMarkerHmms.
    // Provenance + licence: see Resources/README.md.
    private static readonly (string Resource, string MarkerId)[] BundledBacterialPfamMarkers =
    {
        ("PF00380_Ribosomal_S9.hmm", "PF00380"),
        ("PF00410_Ribosomal_S8.hmm", "PF00410"),
        ("PF00466_Ribosomal_L10.hmm", "PF00466"),
        ("PF01025_GrpE.hmm",          "PF01025"),
        ("PF02576_DUF150.hmm",        "PF02576"),
        ("PF03726_PNPase.hmm",        "PF03726"),
    };

    // The Pfam-defined members of the GTDB ar122 (Archaea) domain-level universal single-copy
    // marker set (Parks et al. 2018; GTDB-Tk ar122.marker_info). ar122 is 122 markers (35 Pfam +
    // 87 TIGRFAM); only the 35 Pfam markers are CC0 and bundled. The 87 TIGRFAM-defined markers
    // (CC BY-SA 4.0) are NOT bundled — supply them via LoadMarkerHmms. See Resources/README.md.
    private static readonly (string Resource, string MarkerId)[] BundledArchaealPfamMarkers =
    {
        ("PF00368_HMGCoA_red.hmm",      "PF00368"),
        ("PF00410_Ribosomal_S8.hmm",    "PF00410"),
        ("PF00466_Ribosomal_L10.hmm",   "PF00466"),
        ("PF00687_Ribosomal_L1.hmm",    "PF00687"),
        ("PF00827_Ribosomal_L15e.hmm",  "PF00827"),
        ("PF00900_Ribosomal_S4e.hmm",   "PF00900"),
        ("PF01000_RNA_pol_A_bac.hmm",   "PF01000"),
        ("PF01015_Ribosomal_S3Ae.hmm",  "PF01015"),
        ("PF01090_Ribosomal_S19e.hmm",  "PF01090"),
        ("PF01092_Ribosomal_S6e.hmm",   "PF01092"),
        ("PF01157_Ribosomal_L21e.hmm",  "PF01157"),
        ("PF01191_RNA_pol_Rpb5_C.hmm",  "PF01191"),
        ("PF01194_RNA_pol_N.hmm",       "PF01194"),
        ("PF01198_Ribosomal_L31e.hmm",  "PF01198"),
        ("PF01200_Ribosomal_S28e.hmm",  "PF01200"),
        ("PF01269_Fibrillarin.hmm",     "PF01269"),
        ("PF01280_Ribosomal_L19e.hmm",  "PF01280"),
        ("PF01282_Ribosomal_S24e.hmm",  "PF01282"),
        ("PF01496_V_ATPase_I.hmm",      "PF01496"),
        ("PF01655_Ribosomal_L32e.hmm",  "PF01655"),
        ("PF01798_Nop.hmm",             "PF01798"),
        ("PF01864_DUF46.hmm",           "PF01864"),
        ("PF01866_Diphthamide_syn.hmm", "PF01866"),
        ("PF01868_UPF0086.hmm",         "PF01868"),
        ("PF01984_dsDNA_bind.hmm",      "PF01984"),
        ("PF01990_ATPsynt_F.hmm",       "PF01990"),
        ("PF02006_DUF137.hmm",          "PF02006"),
        ("PF02978_SRP_SPB.hmm",         "PF02978"),
        ("PF03874_RNA_pol_Rpb4.hmm",    "PF03874"),
        ("PF04019_DUF359.hmm",          "PF04019"),
        ("PF04104_DNA_primase_lrg.hmm", "PF04104"),
        ("PF04919_DUF655.hmm",          "PF04919"),
        ("PF07541_EIF_2_alpha.hmm",     "PF07541"),
        ("PF13656_RNA_pol_L_2.hmm",     "PF13656"),
        ("PF13685_FeADH_2.hmm",         "PF13685"),
    };

    /// <summary>
    /// Computes CheckM-style genome completeness and contamination from collocated single-copy
    /// marker SETS and a per-marker copy-count map, using the exact CheckM formula
    /// (Parks et al. 2015, <i>Genome Res</i> 25:1043; reference implementation
    /// <c>MarkerSet.genomeCheck</c>).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Let <c>M</c> be the set of collocated marker sets and <c>G_M</c> the markers identified in
    /// the genome. For each set <c>s ∈ M</c>:
    /// </para>
    /// <list type="bullet">
    /// <item><description><c>present_s = |s ∩ G_M|</c> — markers in <c>s</c> found ≥ 1 time.</description></item>
    /// <item><description><c>multiCopy_s = Σ_{g ∈ s} C_g</c>, where <c>C_g = N − 1</c> for a marker
    /// found <c>N ≥ 1</c> times and <c>0</c> for a missing marker.</description></item>
    /// </list>
    /// <para>
    /// Completeness = <c>100 · (1/|M|) · Σ_{s∈M} present_s/|s|</c>;
    /// Contamination = <c>100 · (1/|M|) · Σ_{s∈M} multiCopy_s/|s|</c>.
    /// </para>
    /// </remarks>
    /// <param name="markerSets">The collocated marker sets <c>M</c>. Empty sets are ignored.</param>
    /// <param name="markerCounts">
    /// Map from marker id to the number of copies found in the bin. Missing keys (or values ≤ 0)
    /// are treated as absent.
    /// </param>
    /// <returns>The completeness/contamination percentages and marker tallies.</returns>
    /// <exception cref="ArgumentNullException">If an argument is null.</exception>
    public static BinMarkerQuality EstimateBinQualityFromMarkerCounts(
        IEnumerable<MarkerSet> markerSets,
        IReadOnlyDictionary<string, int> markerCounts)
    {
        ArgumentNullException.ThrowIfNull(markerSets);
        ArgumentNullException.ThrowIfNull(markerCounts);

        double comp = 0.0;
        double cont = 0.0;
        int usedSets = 0;
        int totalMarkers = 0;
        var distinctPresent = new HashSet<string>(StringComparer.Ordinal);

        foreach (var set in markerSets)
        {
            var ids = set.MarkerIds;
            if (ids is null || ids.Count == 0)
                continue; // empty marker set contributes nothing and would divide by zero

            int present = 0;
            int multiCopy = 0;
            foreach (var marker in ids)
            {
                int count = markerCounts.TryGetValue(marker, out var c) && c > 0 ? c : 0;
                if (count >= 1)
                {
                    present++;                       // C_g counted once toward |s ∩ G_M|
                    distinctPresent.Add(marker);
                }
                if (count > 1)
                    multiCopy += count - 1;           // C_g = N - 1 for a duplicated marker
            }

            comp += (double)present / ids.Count;
            cont += (double)multiCopy / ids.Count;
            usedSets++;
            totalMarkers += ids.Count;
        }

        // |M| = 0 ⇒ no information; completeness and contamination are both 0 (avoid div-by-zero).
        double completeness = usedSets > 0 ? 100.0 * comp / usedSets : 0.0;
        double contamination = usedSets > 0 ? 100.0 * cont / usedSets : 0.0;

        return new BinMarkerQuality(
            Completeness: completeness,
            Contamination: contamination,
            MarkerSetCount: usedSets,
            MarkerCount: totalMarkers,
            MarkersPresent: distinctPresent.Count);
    }

    /// <summary>
    /// Counts, for each marker HMM, how many of a bin's predicted proteins it detects — i.e. the
    /// copy number of that marker in the bin. A protein is a hit for a marker when its Plan7
    /// Viterbi log-odds bit score against the marker's profile meets the marker's bit-score
    /// threshold (the Pfam GA1 gathering threshold). Reuses the <see cref="Plan7ProfileHmm"/> engine.
    /// </summary>
    /// <param name="proteins">The bin's predicted protein sequences (single-letter amino acids).</param>
    /// <param name="markerHmms">The single-copy marker HMMs to score.</param>
    /// <returns>Map from marker id to the number of proteins detected as that marker.</returns>
    /// <exception cref="ArgumentNullException">If an argument is null.</exception>
    public static IReadOnlyDictionary<string, int> DetectMarkers(
        IEnumerable<string> proteins,
        IEnumerable<MarkerHmm> markerHmms)
    {
        ArgumentNullException.ThrowIfNull(proteins);
        ArgumentNullException.ThrowIfNull(markerHmms);

        var hmms = markerHmms.ToList();
        var proteinList = proteins.Where(p => !string.IsNullOrEmpty(p)).ToList();

        var counts = new Dictionary<string, int>(StringComparer.Ordinal);
        foreach (var marker in hmms)
        {
            int copies = 0;
            foreach (var protein in proteinList)
            {
                double nats = marker.Hmm.ViterbiScore(protein);
                double bits = double.IsNegativeInfinity(nats) ? double.NegativeInfinity : nats / LogOddsBitsPerNat;
                if (bits >= marker.BitScoreThreshold)
                    copies++;
            }
            // A marker id may appear once; if duplicated in the input, accumulate.
            counts[marker.MarkerId] = counts.TryGetValue(marker.MarkerId, out var prev) ? prev + copies : copies;
        }
        return counts;
    }

    /// <summary>
    /// Detects single-copy markers in a bin's proteins with the Plan7 HMM engine and computes
    /// CheckM-style completeness and contamination (Parks et al. 2015) over the supplied collocated
    /// marker sets. Opt-in; the TNF/proxy <see cref="BinContigs"/> path and its defaults are unchanged.
    /// </summary>
    /// <param name="proteins">The bin's predicted protein sequences.</param>
    /// <param name="markerSets">The collocated marker sets <c>M</c> (ids must match the HMM marker ids).</param>
    /// <param name="markerHmms">The marker HMMs used to detect the markers.</param>
    /// <returns>The bin's completeness/contamination quality.</returns>
    /// <exception cref="ArgumentNullException">If an argument is null.</exception>
    public static BinMarkerQuality EstimateBinQualityFromMarkers(
        IEnumerable<string> proteins,
        IEnumerable<MarkerSet> markerSets,
        IEnumerable<MarkerHmm> markerHmms)
    {
        ArgumentNullException.ThrowIfNull(proteins);
        ArgumentNullException.ThrowIfNull(markerSets);
        ArgumentNullException.ThrowIfNull(markerHmms);

        var counts = DetectMarkers(proteins, markerHmms);
        return EstimateBinQualityFromMarkerCounts(markerSets, counts);
    }

    /// <summary>
    /// Loads the bundled CC0 universal single-copy ribosomal-protein marker HMMs (9 Pfam families:
    /// S2, S7, S8, S9, S10, S11, S19, L1, L3). Each marker's bit-score gate is its Pfam GA1
    /// gathering threshold (falling back to a conservative default if the profile has none).
    /// </summary>
    /// <remarks>
    /// This is a SMALL universal set, not CheckM's full lineage-specific marker DB. For lineage-
    /// specific marker sets + tree-based placement, supply your own marker HMMs via
    /// <see cref="LoadMarkerHmms(IEnumerable{System.IO.TextReader}, System.Collections.Generic.IReadOnlyDictionary{string, double})"/>.
    /// Provenance + CC0 licence: see <c>Resources/README.md</c>.
    /// </remarks>
    /// <returns>The bundled marker HMMs.</returns>
    public static IReadOnlyList<MarkerHmm> LoadBundledRibosomalMarkerHmms()
        => LoadEmbeddedMarkerHmms(BundledRibosomalMarkers);

    /// <summary>
    /// Loads the bundled Pfam-defined markers of the GTDB <c>bac120</c> Bacteria domain-level
    /// universal single-copy marker set (Parks et al. 2018, <i>Nat Biotechnol</i> 36:996; GTDB-Tk).
    /// bac120 has 120 markers (6 Pfam + 114 TIGRFAM); these are the 6 CC0 Pfam markers
    /// (PF00380, PF00410, PF00466, PF01025, PF02576, PF03726).
    /// </summary>
    /// <remarks>
    /// The 114 TIGRFAM-defined bac120 markers are licensed CC BY-SA 4.0 (not public domain) and are
    /// NOT bundled; supply them via
    /// <see cref="LoadMarkerHmms(IEnumerable{System.IO.TextReader}, System.Collections.Generic.IReadOnlyDictionary{string, double})"/>.
    /// Each marker's bit-score gate is its Pfam GA1 gathering threshold. Provenance + CC0 licence:
    /// see <c>Resources/README.md</c>.
    /// </remarks>
    /// <returns>The 6 bundled bac120 Pfam marker HMMs.</returns>
    public static IReadOnlyList<MarkerHmm> LoadBundledBacterialMarkerHmms()
        => LoadEmbeddedMarkerHmms(BundledBacterialPfamMarkers);

    /// <summary>
    /// Loads the bundled Pfam-defined markers of the GTDB <c>ar122</c> Archaea domain-level
    /// universal single-copy marker set (Parks et al. 2018, <i>Nat Biotechnol</i> 36:996; GTDB-Tk).
    /// ar122 has 122 markers (35 Pfam + 87 TIGRFAM); these are the 35 CC0 Pfam markers.
    /// </summary>
    /// <remarks>
    /// The 87 TIGRFAM-defined ar122 markers are licensed CC BY-SA 4.0 (not public domain) and are
    /// NOT bundled; supply them via
    /// <see cref="LoadMarkerHmms(IEnumerable{System.IO.TextReader}, System.Collections.Generic.IReadOnlyDictionary{string, double})"/>.
    /// Each marker's bit-score gate is its Pfam GA1 gathering threshold. Provenance + CC0 licence:
    /// see <c>Resources/README.md</c>.
    /// </remarks>
    /// <returns>The 35 bundled ar122 Pfam marker HMMs.</returns>
    public static IReadOnlyList<MarkerHmm> LoadBundledArchaealMarkerHmms()
        => LoadEmbeddedMarkerHmms(BundledArchaealPfamMarkers);

    // Shared loader for an embedded-resource marker table: parses each profile with the Plan7
    // engine and gates on its Pfam GA1 gathering threshold (conservative default if absent).
    private static IReadOnlyList<MarkerHmm> LoadEmbeddedMarkerHmms(
        (string Resource, string MarkerId)[] markers)
    {
        var result = new List<MarkerHmm>(markers.Length);
        var asm = typeof(MetagenomicsAnalyzer).Assembly;
        foreach (var (resource, markerId) in markers)
        {
            string resourceName = $"Seqeron.Genomics.Metagenomics.Resources.{resource}";
            using var stream = asm.GetManifestResourceStream(resourceName)
                ?? throw new InvalidOperationException($"Embedded marker HMM not found: {resourceName}");
            using var reader = new System.IO.StreamReader(stream);
            var hmm = Plan7ProfileHmm.Parse(reader);
            double threshold = hmm.GatheringThreshold ?? DefaultMarkerBitScoreThreshold;
            result.Add(new MarkerHmm(markerId, hmm, threshold));
        }
        return result;
    }

    /// <summary>
    /// Builds the default marker SETS for the bundled universal ribosomal markers: one singleton
    /// collocated set per Pfam family (these distinct ribosomal families are scored independently;
    /// CheckM's operon-based collocation grouping requires the full lineage DB — see algorithm doc).
    /// </summary>
    /// <returns>One singleton <see cref="MarkerSet"/> per bundled marker id.</returns>
    public static IReadOnlyList<MarkerSet> BundledRibosomalMarkerSets()
        => BundledRibosomalMarkers
            .Select(m => new MarkerSet(new[] { m.MarkerId }))
            .ToList();

    /// <summary>
    /// Builds the default marker SETS for the bundled bac120 Bacteria Pfam markers: one singleton
    /// collocated set per Pfam family. CheckM's operon-based collocation grouping requires the full
    /// lineage DB — see the algorithm doc — so each universal family is scored independently here.
    /// </summary>
    /// <returns>One singleton <see cref="MarkerSet"/> per bundled bac120 Pfam marker id (6 sets).</returns>
    public static IReadOnlyList<MarkerSet> BundledBacterialMarkerSets()
        => BundledBacterialPfamMarkers
            .Select(m => new MarkerSet(new[] { m.MarkerId }))
            .ToList();

    /// <summary>
    /// Builds the default marker SETS for the bundled ar122 Archaea Pfam markers: one singleton
    /// collocated set per Pfam family. CheckM's operon-based collocation grouping requires the full
    /// lineage DB — see the algorithm doc — so each universal family is scored independently here.
    /// </summary>
    /// <returns>One singleton <see cref="MarkerSet"/> per bundled ar122 Pfam marker id (35 sets).</returns>
    public static IReadOnlyList<MarkerSet> BundledArchaealMarkerSets()
        => BundledArchaealPfamMarkers
            .Select(m => new MarkerSet(new[] { m.MarkerId }))
            .ToList();

    /// <summary>
    /// Caller-supplied marker-HMM loader: parses HMMER3/f ASCII profiles from the given readers
    /// (e.g. lineage-specific Pfam/TIGRFAM marker HMMs from the full CheckM data) into
    /// <see cref="MarkerHmm"/> entries keyed by each profile's accession (falling back to its name).
    /// Lets users with the full CheckM marker database supply their own markers.
    /// </summary>
    /// <param name="hmmReaders">Readers over HMMER3/f ASCII <c>.hmm</c> profiles.</param>
    /// <param name="bitScoreThresholds">
    /// Optional per-marker-id bit-score gates. When absent for a marker, the profile's GA1
    /// gathering threshold is used, falling back to a conservative default.
    /// </param>
    /// <returns>The loaded marker HMMs, in reader order.</returns>
    /// <exception cref="ArgumentNullException">If <paramref name="hmmReaders"/> is null.</exception>
    public static IReadOnlyList<MarkerHmm> LoadMarkerHmms(
        IEnumerable<System.IO.TextReader> hmmReaders,
        IReadOnlyDictionary<string, double>? bitScoreThresholds = null)
    {
        ArgumentNullException.ThrowIfNull(hmmReaders);

        var result = new List<MarkerHmm>();
        foreach (var reader in hmmReaders)
        {
            var hmm = Plan7ProfileHmm.Parse(reader);
            string markerId = !string.IsNullOrEmpty(hmm.Accession) ? hmm.Accession : hmm.Name;
            double threshold =
                bitScoreThresholds is not null && bitScoreThresholds.TryGetValue(markerId, out var t)
                    ? t
                    : hmm.GatheringThreshold ?? DefaultMarkerBitScoreThreshold;
            result.Add(new MarkerHmm(markerId, hmm, threshold));
        }
        return result;
    }

    #endregion
}
