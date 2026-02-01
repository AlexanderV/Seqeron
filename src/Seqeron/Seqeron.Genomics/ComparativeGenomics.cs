using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics;

/// <summary>
/// Provides comparative genomics algorithms for analyzing relationships between genomes.
/// Includes synteny detection, genome rearrangements, and ortholog identification.
/// </summary>
public static class ComparativeGenomics
{
    /// <summary>
    /// Represents a syntenic block between two genomes.
    /// </summary>
    public readonly record struct SyntenicBlock(
        string Genome1Id,
        int Start1,
        int End1,
        string Genome2Id,
        int Start2,
        int End2,
        bool IsInverted,
        int GeneCount,
        double Identity);

    /// <summary>
    /// Represents an orthologous gene pair.
    /// </summary>
    public readonly record struct OrthologPair(
        string Gene1Id,
        string Gene2Id,
        double Identity,
        double Coverage,
        int AlignmentLength);

    /// <summary>
    /// Represents a genome rearrangement event.
    /// </summary>
    public readonly record struct RearrangementEvent(
        RearrangementType Type,
        string GenomeId,
        int Position,
        int Length,
        string? TargetPosition = null);

    /// <summary>
    /// Types of genome rearrangements.
    /// </summary>
    public enum RearrangementType
    {
        Inversion,
        Translocation,
        Deletion,
        Insertion,
        Duplication,
        Transposition
    }

    /// <summary>
    /// Represents a gene for comparative analysis.
    /// </summary>
    public readonly record struct Gene(
        string Id,
        string GenomeId,
        int Start,
        int End,
        char Strand,
        string? Sequence = null);

    /// <summary>
    /// Result of comparative genome analysis.
    /// </summary>
    public readonly record struct ComparisonResult(
        IReadOnlyList<SyntenicBlock> SyntenicBlocks,
        IReadOnlyList<OrthologPair> Orthologs,
        IReadOnlyList<RearrangementEvent> Rearrangements,
        double OverallSynteny,
        int ConservedGenes,
        int GenomeSpecificGenes1,
        int GenomeSpecificGenes2);

    /// <summary>
    /// Finds syntenic blocks between two genomes using gene order.
    /// </summary>
    public static IEnumerable<SyntenicBlock> FindSyntenicBlocks(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        IReadOnlyDictionary<string, string> orthologMap,
        int minBlockSize = 3,
        int maxGap = 5)
    {
        if (genome1Genes.Count == 0 || genome2Genes.Count == 0)
            yield break;

        // Build position maps
        var genome2Positions = new Dictionary<string, int>();
        for (int i = 0; i < genome2Genes.Count; i++)
        {
            genome2Positions[genome2Genes[i].Id] = i;
        }

        // Find runs of collinear orthologs
        var anchors = new List<(int pos1, int pos2, bool strand)>();

        for (int i = 0; i < genome1Genes.Count; i++)
        {
            var gene1 = genome1Genes[i];
            if (orthologMap.TryGetValue(gene1.Id, out string? ortholog) &&
                genome2Positions.TryGetValue(ortholog, out int pos2))
            {
                var gene2 = genome2Genes[pos2];
                bool sameStrand = gene1.Strand == gene2.Strand;
                anchors.Add((i, pos2, sameStrand));
            }
        }

        if (anchors.Count < minBlockSize)
            yield break;

        // Find collinear chains
        var blocks = FindCollinearChains(anchors, minBlockSize, maxGap);

        foreach (var block in blocks)
        {
            if (block.Count < minBlockSize) continue;

            int start1 = genome1Genes[block.Min(a => a.pos1)].Start;
            int end1 = genome1Genes[block.Max(a => a.pos1)].End;
            int start2 = genome2Genes[block.Min(a => a.pos2)].Start;
            int end2 = genome2Genes[block.Max(a => a.pos2)].End;

            bool isInverted = block.First().pos2 > block.Last().pos2;

            yield return new SyntenicBlock(
                Genome1Id: genome1Genes[0].GenomeId,
                Start1: start1,
                End1: end1,
                Genome2Id: genome2Genes[0].GenomeId,
                Start2: Math.Min(start2, end2),
                End2: Math.Max(start2, end2),
                IsInverted: isInverted,
                GeneCount: block.Count,
                Identity: 1.0); // Could compute actual identity
        }
    }

    private static List<List<(int pos1, int pos2, bool strand)>> FindCollinearChains(
        List<(int pos1, int pos2, bool strand)> anchors,
        int minSize,
        int maxGap)
    {
        var chains = new List<List<(int pos1, int pos2, bool strand)>>();

        // Sort by position in genome 1
        anchors = anchors.OrderBy(a => a.pos1).ToList();

        var currentChain = new List<(int pos1, int pos2, bool strand)>();
        int? lastPos2 = null;
        bool? direction = null;

        foreach (var anchor in anchors)
        {
            if (lastPos2 == null)
            {
                currentChain.Add(anchor);
                lastPos2 = anchor.pos2;
                continue;
            }

            int diff = anchor.pos2 - lastPos2.Value;
            bool currentDir = diff > 0;

            // Check if this anchor continues the chain
            bool continuesChain = Math.Abs(diff) <= maxGap &&
                                  Math.Abs(diff) >= 1 &&
                                  (direction == null || direction == currentDir);

            if (continuesChain)
            {
                currentChain.Add(anchor);
                direction = currentDir;
            }
            else
            {
                // Save current chain and start new one
                if (currentChain.Count >= minSize)
                {
                    chains.Add(new List<(int pos1, int pos2, bool strand)>(currentChain));
                }
                currentChain.Clear();
                currentChain.Add(anchor);
                direction = null;
            }

            lastPos2 = anchor.pos2;
        }

        if (currentChain.Count >= minSize)
        {
            chains.Add(currentChain);
        }

        return chains;
    }

    /// <summary>
    /// Identifies orthologous gene pairs using sequence similarity.
    /// </summary>
    public static IEnumerable<OrthologPair> FindOrthologs(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        double minIdentity = 0.3,
        double minCoverage = 0.5)
    {
        foreach (var gene1 in genome1Genes.Where(g => !string.IsNullOrEmpty(g.Sequence)))
        {
            OrthologPair? bestMatch = null;

            foreach (var gene2 in genome2Genes.Where(g => !string.IsNullOrEmpty(g.Sequence)))
            {
                var (identity, coverage, alignLen) = CalculateSequenceSimilarity(
                    gene1.Sequence!, gene2.Sequence!);

                if (identity >= minIdentity && coverage >= minCoverage)
                {
                    if (bestMatch == null || identity > bestMatch.Value.Identity)
                    {
                        bestMatch = new OrthologPair(
                            Gene1Id: gene1.Id,
                            Gene2Id: gene2.Id,
                            Identity: identity,
                            Coverage: coverage,
                            AlignmentLength: alignLen);
                    }
                }
            }

            if (bestMatch != null)
                yield return bestMatch.Value;
        }
    }

    /// <summary>
    /// Finds reciprocal best hits (RBH) for ortholog identification.
    /// </summary>
    public static IEnumerable<OrthologPair> FindReciprocalBestHits(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        double minIdentity = 0.3)
    {
        // Best hits from genome1 to genome2
        var bestHits1To2 = new Dictionary<string, (string gene2, double score)>();
        // Best hits from genome2 to genome1
        var bestHits2To1 = new Dictionary<string, (string gene1, double score)>();

        // Find best hits in both directions
        foreach (var gene1 in genome1Genes.Where(g => !string.IsNullOrEmpty(g.Sequence)))
        {
            double bestScore = 0;
            string? bestMatch = null;

            foreach (var gene2 in genome2Genes.Where(g => !string.IsNullOrEmpty(g.Sequence)))
            {
                var (identity, coverage, _) = CalculateSequenceSimilarity(
                    gene1.Sequence!, gene2.Sequence!);
                double score = identity * coverage;

                if (score > bestScore && identity >= minIdentity)
                {
                    bestScore = score;
                    bestMatch = gene2.Id;
                }
            }

            if (bestMatch != null)
                bestHits1To2[gene1.Id] = (bestMatch, bestScore);
        }

        foreach (var gene2 in genome2Genes.Where(g => !string.IsNullOrEmpty(g.Sequence)))
        {
            double bestScore = 0;
            string? bestMatch = null;

            foreach (var gene1 in genome1Genes.Where(g => !string.IsNullOrEmpty(g.Sequence)))
            {
                var (identity, coverage, _) = CalculateSequenceSimilarity(
                    gene2.Sequence!, gene1.Sequence!);
                double score = identity * coverage;

                if (score > bestScore && identity >= minIdentity)
                {
                    bestScore = score;
                    bestMatch = gene1.Id;
                }
            }

            if (bestMatch != null)
                bestHits2To1[gene2.Id] = (bestMatch, bestScore);
        }

        // Find reciprocal best hits
        foreach (var (gene1, (gene2, score)) in bestHits1To2)
        {
            if (bestHits2To1.TryGetValue(gene2, out var reverse) && reverse.gene1 == gene1)
            {
                yield return new OrthologPair(
                    Gene1Id: gene1,
                    Gene2Id: gene2,
                    Identity: score,
                    Coverage: 1.0,
                    AlignmentLength: 0);
            }
        }
    }

    private static (double identity, double coverage, int alignLength) CalculateSequenceSimilarity(
        string seq1, string seq2)
    {
        // Simple k-mer based similarity (faster than full alignment)
        const int k = 5;

        if (seq1.Length < k || seq2.Length < k)
            return (0, 0, 0);

        var kmers1 = new HashSet<string>();
        for (int i = 0; i <= seq1.Length - k; i++)
            kmers1.Add(seq1.Substring(i, k).ToUpperInvariant());

        var kmers2 = new HashSet<string>();
        for (int i = 0; i <= seq2.Length - k; i++)
            kmers2.Add(seq2.Substring(i, k).ToUpperInvariant());

        int shared = kmers1.Intersect(kmers2).Count();
        int total = kmers1.Union(kmers2).Count();

        double identity = total > 0 ? (double)shared / total : 0;
        double coverage = (double)shared / Math.Max(kmers1.Count, kmers2.Count);
        int alignLen = Math.Min(seq1.Length, seq2.Length);

        return (identity, coverage, alignLen);
    }

    /// <summary>
    /// Detects genome rearrangements by comparing gene orders.
    /// </summary>
    public static IEnumerable<RearrangementEvent> DetectRearrangements(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        IReadOnlyDictionary<string, string> orthologMap)
    {
        // Build position and strand maps
        var gene1Positions = genome1Genes
            .Select((g, i) => (g, i))
            .ToDictionary(x => x.g.Id, x => (pos: x.i, strand: x.g.Strand));

        var gene2Positions = genome2Genes
            .Select((g, i) => (g, i))
            .ToDictionary(x => x.g.Id, x => (pos: x.i, strand: x.g.Strand));

        // Map genome1 positions to genome2 positions via orthologs
        var mappedPositions = new List<(int pos1, int pos2, bool inverted)>();

        foreach (var gene1 in genome1Genes)
        {
            if (orthologMap.TryGetValue(gene1.Id, out string? ortholog) &&
                gene2Positions.TryGetValue(ortholog, out var gene2Info))
            {
                bool inverted = gene1.Strand != genome2Genes[gene2Info.pos].Strand;
                mappedPositions.Add((gene1Positions[gene1.Id].pos, gene2Info.pos, inverted));
            }
        }

        if (mappedPositions.Count < 2)
            yield break;

        mappedPositions = mappedPositions.OrderBy(p => p.pos1).ToList();

        // Detect inversions
        for (int i = 0; i < mappedPositions.Count - 1; i++)
        {
            var current = mappedPositions[i];
            var next = mappedPositions[i + 1];

            // Check for inversion (position decreases or strand changes)
            if (current.pos2 > next.pos2 || current.inverted != next.inverted)
            {
                yield return new RearrangementEvent(
                    Type: RearrangementType.Inversion,
                    GenomeId: genome1Genes[0].GenomeId,
                    Position: genome1Genes[current.pos1].Start,
                    Length: genome1Genes[next.pos1].End - genome1Genes[current.pos1].Start);
            }

            // Check for large gaps (potential deletions/insertions)
            int gap1 = next.pos1 - current.pos1;
            int gap2 = Math.Abs(next.pos2 - current.pos2);

            if (gap1 > 1 && gap2 <= 1)
            {
                // Gap in genome1 but not genome2 - potential deletion in genome2
                yield return new RearrangementEvent(
                    Type: RearrangementType.Deletion,
                    GenomeId: genome2Genes[0].GenomeId,
                    Position: genome2Genes[current.pos2].End,
                    Length: gap1 - 1);
            }
            else if (gap2 > 1 && gap1 <= 1)
            {
                // Gap in genome2 but not genome1 - potential insertion in genome2
                yield return new RearrangementEvent(
                    Type: RearrangementType.Insertion,
                    GenomeId: genome2Genes[0].GenomeId,
                    Position: genome2Genes[current.pos2].End,
                    Length: gap2 - 1);
            }
        }
    }

    /// <summary>
    /// Performs comprehensive comparative analysis between two genomes.
    /// </summary>
    public static ComparisonResult CompareGenomes(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        double minOrthologIdentity = 0.3,
        int minSyntenicBlockSize = 3)
    {
        // Find orthologs
        var orthologs = FindReciprocalBestHits(genome1Genes, genome2Genes, minOrthologIdentity).ToList();

        // Build ortholog map
        var orthologMap = orthologs.ToDictionary(o => o.Gene1Id, o => o.Gene2Id);

        // Find syntenic blocks
        var syntenicBlocks = FindSyntenicBlocks(
            genome1Genes, genome2Genes, orthologMap, minSyntenicBlockSize).ToList();

        // Detect rearrangements
        var rearrangements = DetectRearrangements(genome1Genes, genome2Genes, orthologMap).ToList();

        // Calculate statistics
        var orthologGenes1 = new HashSet<string>(orthologs.Select(o => o.Gene1Id));
        var orthologGenes2 = new HashSet<string>(orthologs.Select(o => o.Gene2Id));

        int specific1 = genome1Genes.Count(g => !orthologGenes1.Contains(g.Id));
        int specific2 = genome2Genes.Count(g => !orthologGenes2.Contains(g.Id));

        double synteny = syntenicBlocks.Count > 0
            ? (double)syntenicBlocks.Sum(b => b.GeneCount) / Math.Min(genome1Genes.Count, genome2Genes.Count)
            : 0;

        return new ComparisonResult(
            SyntenicBlocks: syntenicBlocks,
            Orthologs: orthologs,
            Rearrangements: rearrangements,
            OverallSynteny: Math.Min(1.0, synteny),
            ConservedGenes: orthologs.Count,
            GenomeSpecificGenes1: specific1,
            GenomeSpecificGenes2: specific2);
    }

    /// <summary>
    /// Calculates the reversal distance between two gene orders.
    /// Uses breakpoint-based approximation for signed permutations.
    /// </summary>
    public static int CalculateReversalDistance(
        IReadOnlyList<int> permutation1,
        IReadOnlyList<int> permutation2)
    {
        if (permutation1.Count != permutation2.Count)
            throw new ArgumentException("Permutations must have the same length");

        int n = permutation1.Count;
        if (n <= 1) return 0;

        // Convert to relative permutation
        var positionMap = new Dictionary<int, int>();
        for (int i = 0; i < n; i++)
            positionMap[permutation2[i]] = i;

        var relative = permutation1.Select(x => positionMap[x]).ToList();

        // Count breakpoints: positions where consecutive elements are not adjacent
        // Include boundaries as implicit breakpoints
        int breakpoints = 0;

        // Check if first element is in correct position (extended permutation with 0 at start)
        if (relative[0] != 0)
            breakpoints++;

        // Check adjacency between consecutive elements
        for (int i = 0; i < n - 1; i++)
        {
            if (Math.Abs(relative[i + 1] - relative[i]) != 1)
                breakpoints++;
        }

        // Check if last element is in correct position (extended permutation with n at end)
        if (relative[n - 1] != n - 1)
            breakpoints++;

        // Lower bound on reversal distance (breakpoint distance / 2)
        return (breakpoints + 1) / 2;
    }

    /// <summary>
    /// Finds conserved gene clusters (genes that appear together in multiple genomes).
    /// </summary>
    public static IEnumerable<IReadOnlyList<string>> FindConservedClusters(
        IReadOnlyList<IReadOnlyList<Gene>> genomes,
        IReadOnlyDictionary<string, string> orthologGroups,
        int minClusterSize = 3,
        int maxGap = 2)
    {
        if (genomes.Count < 2)
            yield break;

        // Find gene clusters in first genome
        var genome1 = genomes[0];
        var clusters = new List<List<string>>();

        for (int start = 0; start < genome1.Count; start++)
        {
            var cluster = new List<string>();

            for (int i = start; i < Math.Min(start + 20, genome1.Count); i++)
            {
                if (orthologGroups.ContainsKey(genome1[i].Id))
                {
                    cluster.Add(orthologGroups[genome1[i].Id]);
                }
            }

            if (cluster.Count >= minClusterSize)
            {
                clusters.Add(cluster);
            }
        }

        // Check which clusters are conserved in other genomes
        foreach (var cluster in clusters)
        {
            bool conserved = true;

            for (int g = 1; g < genomes.Count && conserved; g++)
            {
                var genome = genomes[g];
                var geneGroups = genome
                    .Where(gene => orthologGroups.ContainsKey(gene.Id))
                    .Select(gene => orthologGroups[gene.Id])
                    .ToList();

                // Check if cluster genes appear within maxGap of each other
                var clusterSet = new HashSet<string>(cluster);
                int found = 0;
                int lastFoundPos = -maxGap - 1;

                for (int i = 0; i < geneGroups.Count; i++)
                {
                    if (clusterSet.Contains(geneGroups[i]))
                    {
                        if (i - lastFoundPos <= maxGap + 1)
                        {
                            found++;
                        }
                        lastFoundPos = i;
                    }
                }

                conserved = found >= minClusterSize;
            }

            if (conserved)
            {
                yield return cluster;
            }
        }
    }

    /// <summary>
    /// Calculates Average Nucleotide Identity (ANI) between two genomes.
    /// </summary>
    public static double CalculateANI(
        string genome1Sequence,
        string genome2Sequence,
        int fragmentSize = 1000,
        double minFragmentIdentity = 0.7)
    {
        if (string.IsNullOrEmpty(genome1Sequence) || string.IsNullOrEmpty(genome2Sequence))
            return 0;

        var identities = new List<double>();

        // Fragment genome1 and compare to genome2
        for (int i = 0; i <= genome1Sequence.Length - fragmentSize; i += fragmentSize / 2)
        {
            string fragment = genome1Sequence.Substring(i, fragmentSize);
            double bestIdentity = FindBestFragmentMatch(fragment, genome2Sequence);

            if (bestIdentity >= minFragmentIdentity)
            {
                identities.Add(bestIdentity);
            }
        }

        return identities.Count > 0 ? identities.Average() : 0;
    }

    private static double FindBestFragmentMatch(string fragment, string genome)
    {
        // Use SuffixTree for efficient longest common substring search
        var suffixTree = global::SuffixTree.SuffixTree.Build(genome.ToUpperInvariant());
        string lcs = suffixTree.LongestCommonSubstring(fragment.ToUpperInvariant());

        // Calculate identity based on LCS length relative to fragment length
        double identity = fragment.Length > 0 ? (double)lcs.Length / fragment.Length : 0;

        return Math.Min(identity, 1.0);
    }

    /// <summary>
    /// Generates a dot plot comparison between two sequences.
    /// Uses SuffixTree for efficient O(m+k) word matching.
    /// </summary>
    public static IEnumerable<(int x, int y)> GenerateDotPlot(
        string sequence1,
        string sequence2,
        int wordSize = 10,
        int stepSize = 1)
    {
        if (string.IsNullOrEmpty(sequence1) || string.IsNullOrEmpty(sequence2))
            yield break;

        // Build SuffixTree on sequence2 for efficient pattern matching
        var suffixTree = global::SuffixTree.SuffixTree.Build(sequence2.ToUpperInvariant());

        // Find matching words from sequence1 in sequence2
        for (int i = 0; i <= sequence1.Length - wordSize; i += stepSize)
        {
            string word = sequence1.Substring(i, wordSize).ToUpperInvariant();
            var positions = suffixTree.FindAllOccurrences(word);

            foreach (int j in positions)
            {
                yield return (i, j);
            }
        }
    }
}
