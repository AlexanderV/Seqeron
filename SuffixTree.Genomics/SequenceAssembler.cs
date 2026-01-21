using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace SuffixTree.Genomics;

/// <summary>
/// Provides genome assembly algorithms for constructing contigs from sequence reads.
/// Supports overlap-layout-consensus and de Bruijn graph approaches.
/// </summary>
public static class SequenceAssembler
{
    /// <summary>
    /// Result of sequence assembly.
    /// </summary>
    public readonly record struct AssemblyResult(
        IReadOnlyList<string> Contigs,
        int TotalReads,
        int AssembledReads,
        double N50,
        int LongestContig,
        int TotalLength);

    /// <summary>
    /// Represents an overlap between two reads.
    /// </summary>
    public readonly record struct Overlap(
        int ReadIndex1,
        int ReadIndex2,
        int OverlapLength,
        int Position1,
        int Position2);

    /// <summary>
    /// Assembly parameters.
    /// </summary>
    public readonly record struct AssemblyParameters(
        int MinOverlap = 20,
        double MinIdentity = 0.9,
        int KmerSize = 31,
        int MinContigLength = 100);

    /// <summary>
    /// Assembles reads using overlap-layout-consensus approach.
    /// </summary>
    public static AssemblyResult AssembleOLC(
        IReadOnlyList<string> reads,
        AssemblyParameters? parameters = null)
    {
        var param = parameters ?? new AssemblyParameters();

        if (reads == null || reads.Count == 0)
            return new AssemblyResult(Array.Empty<string>(), 0, 0, 0, 0, 0);

        // Step 1: Find overlaps
        var overlaps = FindAllOverlaps(reads, param.MinOverlap, param.MinIdentity);

        // Step 2: Build overlap graph and find layout
        var contigs = BuildContigsFromOverlaps(reads, overlaps, param);

        // Filter by minimum length
        contigs = contigs.Where(c => c.Length >= param.MinContigLength).ToList();

        // Calculate statistics
        var stats = CalculateStats(contigs, reads.Count);

        return stats;
    }

    /// <summary>
    /// Assembles reads using de Bruijn graph approach.
    /// </summary>
    public static AssemblyResult AssembleDeBruijn(
        IReadOnlyList<string> reads,
        AssemblyParameters? parameters = null)
    {
        var param = parameters ?? new AssemblyParameters();

        if (reads == null || reads.Count == 0)
            return new AssemblyResult(Array.Empty<string>(), 0, 0, 0, 0, 0);

        int k = param.KmerSize;

        // Build de Bruijn graph
        var graph = BuildDeBruijnGraph(reads, k);

        // Find Eulerian paths / contigs
        var contigs = FindContigs(graph, k);

        // Filter by minimum length
        contigs = contigs.Where(c => c.Length >= param.MinContigLength).ToList();

        return CalculateStats(contigs, reads.Count);
    }

    /// <summary>
    /// Finds all overlaps between reads above minimum threshold.
    /// </summary>
    public static IReadOnlyList<Overlap> FindAllOverlaps(
        IReadOnlyList<string> reads,
        int minOverlap = 20,
        double minIdentity = 0.9)
    {
        var overlaps = new List<Overlap>();

        for (int i = 0; i < reads.Count; i++)
        {
            for (int j = 0; j < reads.Count; j++)
            {
                if (i == j) continue;

                var overlap = FindOverlap(reads[i], reads[j], minOverlap, minIdentity);
                if (overlap.HasValue)
                {
                    overlaps.Add(new Overlap(i, j, overlap.Value.length,
                        overlap.Value.pos1, overlap.Value.pos2));
                }
            }
        }

        return overlaps;
    }

    /// <summary>
    /// Finds suffix-prefix overlap between two sequences.
    /// </summary>
    public static (int length, int pos1, int pos2)? FindOverlap(
        string seq1, string seq2,
        int minOverlap = 20,
        double minIdentity = 0.9)
    {
        // Check if suffix of seq1 overlaps with prefix of seq2
        int maxPossible = Math.Min(seq1.Length, seq2.Length);

        for (int overlapLen = maxPossible; overlapLen >= minOverlap; overlapLen--)
        {
            string suffix = seq1.Substring(seq1.Length - overlapLen);
            string prefix = seq2.Substring(0, overlapLen);

            double identity = CalculateIdentity(suffix, prefix);
            if (identity >= minIdentity)
            {
                return (overlapLen, seq1.Length - overlapLen, 0);
            }
        }

        return null;
    }

    /// <summary>
    /// Calculates sequence identity between two strings of equal length.
    /// </summary>
    public static double CalculateIdentity(string seq1, string seq2)
    {
        if (seq1.Length != seq2.Length) return 0;
        if (seq1.Length == 0) return 1;

        int matches = 0;
        for (int i = 0; i < seq1.Length; i++)
        {
            if (char.ToUpperInvariant(seq1[i]) == char.ToUpperInvariant(seq2[i]))
                matches++;
        }

        return (double)matches / seq1.Length;
    }

    private static List<string> BuildContigsFromOverlaps(
        IReadOnlyList<string> reads,
        IReadOnlyList<Overlap> overlaps,
        AssemblyParameters param)
    {
        var contigs = new List<string>();
        var used = new HashSet<int>();

        // Build adjacency: for each read, best successor
        var bestSuccessor = new Dictionary<int, (int next, int overlap)>();
        var hasPredecessor = new HashSet<int>();

        foreach (var ov in overlaps.OrderByDescending(o => o.OverlapLength))
        {
            int r1 = ov.ReadIndex1;
            int r2 = ov.ReadIndex2;

            if (!bestSuccessor.ContainsKey(r1))
            {
                bestSuccessor[r1] = (r2, ov.OverlapLength);
                hasPredecessor.Add(r2);
            }
        }

        // Find starting reads (no predecessor in best-overlap chain)
        var starters = new List<int>();
        for (int i = 0; i < reads.Count; i++)
        {
            if (!hasPredecessor.Contains(i))
                starters.Add(i);
        }

        // Build contigs from each starter
        foreach (int start in starters)
        {
            if (used.Contains(start)) continue;

            var sb = new StringBuilder();
            int current = start;

            while (current != -1 && !used.Contains(current))
            {
                used.Add(current);
                string read = reads[current];

                if (sb.Length == 0)
                {
                    sb.Append(read);
                }
                else if (bestSuccessor.TryGetValue(current, out var info))
                {
                    // Already added previous, now extend
                }

                if (bestSuccessor.TryGetValue(current, out var next))
                {
                    int overlap = next.overlap;
                    string nextRead = reads[next.next];
                    if (!used.Contains(next.next))
                    {
                        sb.Append(nextRead.Substring(overlap));
                    }
                    current = next.next;
                }
                else
                {
                    current = -1;
                }
            }

            if (sb.Length > 0)
                contigs.Add(sb.ToString());
        }

        // Add unused reads as singleton contigs
        for (int i = 0; i < reads.Count; i++)
        {
            if (!used.Contains(i) && reads[i].Length >= param.MinContigLength)
            {
                contigs.Add(reads[i]);
            }
        }

        return contigs;
    }

    private static Dictionary<string, List<string>> BuildDeBruijnGraph(
        IReadOnlyList<string> reads, int k)
    {
        var graph = new Dictionary<string, List<string>>();

        foreach (string read in reads)
        {
            for (int i = 0; i <= read.Length - k; i++)
            {
                string kmer = read.Substring(i, k);
                string prefix = kmer.Substring(0, k - 1);
                string suffix = kmer.Substring(1);

                if (!graph.ContainsKey(prefix))
                    graph[prefix] = new List<string>();

                graph[prefix].Add(suffix);
            }
        }

        return graph;
    }

    private static List<string> FindContigs(
        Dictionary<string, List<string>> graph, int k)
    {
        var contigs = new List<string>();

        // Count in-degrees and out-degrees
        var outDegree = new Dictionary<string, int>();
        var inDegree = new Dictionary<string, int>();

        foreach (var kvp in graph)
        {
            outDegree[kvp.Key] = kvp.Value.Count;
            foreach (string target in kvp.Value)
            {
                inDegree[target] = inDegree.GetValueOrDefault(target, 0) + 1;
            }
        }

        // Find all nodes
        var allNodes = new HashSet<string>(graph.Keys);
        foreach (var targets in graph.Values)
            foreach (var t in targets)
                allNodes.Add(t);

        // Find starting points (out > in or sources)
        var startNodes = new List<string>();
        foreach (string node in allNodes)
        {
            int outD = outDegree.GetValueOrDefault(node, 0);
            int inD = inDegree.GetValueOrDefault(node, 0);

            if (outD > inD || (outD > 0 && inD == 0))
                startNodes.Add(node);
        }

        if (startNodes.Count == 0 && graph.Count > 0)
            startNodes.Add(graph.Keys.First());

        // Trace paths from starting nodes
        var usedEdges = new Dictionary<string, HashSet<int>>();

        foreach (string start in startNodes)
        {
            var contig = TraceContig(graph, start, usedEdges);
            if (contig.Length >= k)
                contigs.Add(contig);
        }

        return contigs;
    }

    private static string TraceContig(
        Dictionary<string, List<string>> graph,
        string start,
        Dictionary<string, HashSet<int>> usedEdges)
    {
        var sb = new StringBuilder(start);
        string current = start;

        while (graph.TryGetValue(current, out var neighbors) && neighbors.Count > 0)
        {
            if (!usedEdges.ContainsKey(current))
                usedEdges[current] = new HashSet<int>();

            // Find unused edge
            int edgeIndex = -1;
            for (int i = 0; i < neighbors.Count; i++)
            {
                if (!usedEdges[current].Contains(i))
                {
                    edgeIndex = i;
                    break;
                }
            }

            if (edgeIndex == -1) break;

            usedEdges[current].Add(edgeIndex);
            string next = neighbors[edgeIndex];
            sb.Append(next[^1]); // Append last character
            current = next;
        }

        return sb.ToString();
    }

    /// <summary>
    /// Calculates N50 and other assembly statistics.
    /// </summary>
    public static AssemblyResult CalculateStats(
        IReadOnlyList<string> contigs, int totalReads)
    {
        if (contigs.Count == 0)
            return new AssemblyResult(contigs, totalReads, 0, 0, 0, 0);

        var sortedLengths = contigs.Select(c => c.Length).OrderDescending().ToList();
        int totalLength = sortedLengths.Sum();
        int halfLength = totalLength / 2;

        // Calculate N50
        int cumulative = 0;
        double n50 = 0;
        foreach (int len in sortedLengths)
        {
            cumulative += len;
            if (cumulative >= halfLength)
            {
                n50 = len;
                break;
            }
        }

        return new AssemblyResult(
            contigs,
            totalReads,
            totalReads, // Simplified: assume all reads used
            n50,
            sortedLengths.First(),
            totalLength);
    }

    /// <summary>
    /// Merges two overlapping contigs.
    /// </summary>
    public static string MergeContigs(string contig1, string contig2, int overlapLength)
    {
        if (overlapLength <= 0 || overlapLength > Math.Min(contig1.Length, contig2.Length))
            return contig1 + contig2;

        return contig1 + contig2.Substring(overlapLength);
    }

    /// <summary>
    /// Scaffolds contigs using paired-end information.
    /// </summary>
    public static IReadOnlyList<string> Scaffold(
        IReadOnlyList<string> contigs,
        IReadOnlyList<(int contig1, int contig2, int gapSize)> links,
        char gapCharacter = 'N')
    {
        if (contigs.Count == 0) return contigs;

        var scaffolds = new List<string>();
        var used = new HashSet<int>();

        // Group links by first contig
        var linkMap = links.GroupBy(l => l.contig1)
            .ToDictionary(g => g.Key, g => g.ToList());

        for (int i = 0; i < contigs.Count; i++)
        {
            if (used.Contains(i)) continue;

            var sb = new StringBuilder(contigs[i]);
            used.Add(i);

            // Follow links
            int current = i;
            while (linkMap.TryGetValue(current, out var nextLinks))
            {
                var link = nextLinks.FirstOrDefault(l => !used.Contains(l.contig2));
                if (link.contig2 == 0 && link.gapSize == 0 && used.Contains(0))
                    break;

                if (!used.Contains(link.contig2))
                {
                    // Add gap
                    sb.Append(new string(gapCharacter, Math.Max(1, link.gapSize)));
                    sb.Append(contigs[link.contig2]);
                    used.Add(link.contig2);
                    current = link.contig2;
                }
                else
                {
                    break;
                }
            }

            scaffolds.Add(sb.ToString());
        }

        return scaffolds;
    }

    /// <summary>
    /// Calculates coverage depth at each position of the reference.
    /// </summary>
    public static int[] CalculateCoverage(string reference, IReadOnlyList<string> reads, int minOverlap = 20)
    {
        var coverage = new int[reference.Length];

        foreach (string read in reads)
        {
            // Find where read maps to reference
            int pos = FindBestAlignment(reference, read, minOverlap);
            if (pos >= 0)
            {
                for (int i = pos; i < pos + read.Length && i < reference.Length; i++)
                {
                    coverage[i]++;
                }
            }
        }

        return coverage;
    }

    private static int FindBestAlignment(string reference, string read, int minOverlap)
    {
        int bestPos = -1;
        int bestScore = minOverlap - 1;

        for (int pos = 0; pos <= reference.Length - read.Length; pos++)
        {
            int matches = 0;
            for (int i = 0; i < read.Length; i++)
            {
                if (char.ToUpperInvariant(reference[pos + i]) ==
                    char.ToUpperInvariant(read[i]))
                    matches++;
            }

            if (matches > bestScore)
            {
                bestScore = matches;
                bestPos = pos;
            }
        }

        return bestPos;
    }

    /// <summary>
    /// Computes the consensus sequence from multiple aligned reads.
    /// </summary>
    public static string ComputeConsensus(IReadOnlyList<string> alignedReads)
    {
        if (alignedReads.Count == 0) return "";

        int length = alignedReads[0].Length;
        var sb = new StringBuilder();

        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int>();
            foreach (string read in alignedReads)
            {
                if (pos < read.Length)
                {
                    char c = char.ToUpperInvariant(read[pos]);
                    if (c != '-' && c != 'N')
                    {
                        counts[c] = counts.GetValueOrDefault(c, 0) + 1;
                    }
                }
            }

            if (counts.Count > 0)
            {
                char consensus = counts.MaxBy(kvp => kvp.Value).Key;
                sb.Append(consensus);
            }
            else
            {
                sb.Append('N');
            }
        }

        return sb.ToString();
    }

    /// <summary>
    /// Quality trims reads based on quality scores.
    /// </summary>
    public static IReadOnlyList<string> QualityTrimReads(
        IReadOnlyList<(string sequence, string quality)> reads,
        int minQuality = 20,
        int minLength = 50)
    {
        var trimmed = new List<string>();

        foreach (var (sequence, quality) in reads)
        {
            int start = 0;
            int end = sequence.Length;

            // Trim from start
            while (start < end && (quality[start] - 33) < minQuality)
                start++;

            // Trim from end
            while (end > start && (quality[end - 1] - 33) < minQuality)
                end--;

            if (end - start >= minLength)
            {
                trimmed.Add(sequence.Substring(start, end - start));
            }
        }

        return trimmed;
    }

    /// <summary>
    /// Error corrects reads using k-mer frequency analysis.
    /// </summary>
    public static IReadOnlyList<string> ErrorCorrectReads(
        IReadOnlyList<string> reads,
        int kmerSize = 21,
        int minKmerFrequency = 3)
    {
        // Build k-mer frequency table
        var kmerCounts = new Dictionary<string, int>();
        foreach (string read in reads)
        {
            for (int i = 0; i <= read.Length - kmerSize; i++)
            {
                string kmer = read.Substring(i, kmerSize).ToUpperInvariant();
                kmerCounts[kmer] = kmerCounts.GetValueOrDefault(kmer, 0) + 1;
            }
        }

        var corrected = new List<string>();

        foreach (string read in reads)
        {
            var sb = new StringBuilder(read);

            for (int i = 0; i <= read.Length - kmerSize; i++)
            {
                string kmer = read.Substring(i, kmerSize).ToUpperInvariant();

                if (kmerCounts.GetValueOrDefault(kmer, 0) < minKmerFrequency)
                {
                    // Try to correct by substituting middle base
                    int midPos = i + kmerSize / 2;
                    char original = char.ToUpperInvariant(sb[midPos]);

                    foreach (char replacement in new[] { 'A', 'C', 'G', 'T' })
                    {
                        if (replacement == original) continue;

                        sb[midPos] = replacement;
                        string newKmer = sb.ToString().Substring(i, kmerSize).ToUpperInvariant();

                        if (kmerCounts.GetValueOrDefault(newKmer, 0) >= minKmerFrequency)
                        {
                            break; // Keep correction
                        }
                        else
                        {
                            sb[midPos] = original; // Revert
                        }
                    }
                }
            }

            corrected.Add(sb.ToString());
        }

        return corrected;
    }
}
