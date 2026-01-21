using System;
using System.Collections.Generic;
using System.Linq;

namespace SuffixTree.Genomics;

/// <summary>
/// Provides genome assembly quality assessment algorithms.
/// Includes N50/L50, contiguity metrics, completeness analysis, and repeat detection.
/// </summary>
public static class GenomeAssemblyAnalyzer
{
    #region Records

    /// <summary>
    /// Represents an assembled sequence (contig or scaffold).
    /// </summary>
    public readonly record struct AssembledSequence(
        string Id,
        string Sequence,
        int Length,
        double GcContent,
        int GapCount,
        int TotalGapLength);

    /// <summary>
    /// Assembly statistics summary.
    /// </summary>
    public readonly record struct AssemblyStatistics(
        int TotalSequences,
        long TotalLength,
        long TotalLengthNoGaps,
        int N50,
        int L50,
        int N90,
        int L90,
        int LargestContig,
        int SmallestContig,
        double MeanLength,
        double MedianLength,
        double GcContent,
        int TotalGaps,
        long TotalGapLength,
        double GapPercentage);

    /// <summary>
    /// Nx/Lx statistics for various thresholds.
    /// </summary>
    public readonly record struct NxStatistics(
        int Threshold,
        int Nx,
        int Lx,
        long CumulativeLength);

    /// <summary>
    /// Gap information in assembly.
    /// </summary>
    public readonly record struct GapInfo(
        string SequenceId,
        int Start,
        int End,
        int Length,
        string GapType);

    /// <summary>
    /// Repeat annotation.
    /// </summary>
    public readonly record struct RepeatAnnotation(
        string SequenceId,
        int Start,
        int End,
        string RepeatClass,
        string RepeatFamily,
        double DivergencePercent,
        char Strand);

    /// <summary>
    /// BUSCO-like completeness result.
    /// </summary>
    public readonly record struct CompletenessResult(
        int TotalGenes,
        int Complete,
        int CompleteSingleCopy,
        int CompleteDuplicated,
        int Fragmented,
        int Missing,
        double CompletenessPercent,
        double DuplicationPercent);

    /// <summary>
    /// Scaffold gap structure.
    /// </summary>
    public readonly record struct ScaffoldStructure(
        string ScaffoldId,
        IReadOnlyList<(string ContigId, int Start, int End)> Contigs,
        IReadOnlyList<GapInfo> Gaps,
        int TotalLength,
        int ContigLength,
        int GapLength);

    /// <summary>
    /// Assembly comparison result.
    /// </summary>
    public readonly record struct AssemblyComparison(
        string Assembly1Name,
        string Assembly2Name,
        double AlignedFraction1,
        double AlignedFraction2,
        int Breakpoints,
        int Inversions,
        int Translocations,
        double SequenceIdentity);

    #endregion

    #region Basic Statistics

    /// <summary>
    /// Calculates comprehensive assembly statistics.
    /// </summary>
    public static AssemblyStatistics CalculateStatistics(
        IEnumerable<(string Id, string Sequence)> sequences)
    {
        var seqList = sequences.Select(s => ParseSequence(s.Id, s.Sequence)).ToList();

        if (seqList.Count == 0)
        {
            return new AssemblyStatistics(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }

        var lengths = seqList.Select(s => s.Length).OrderByDescending(l => l).ToList();
        long totalLength = lengths.Sum(l => (long)l);
        long totalNoGaps = seqList.Sum(s => (long)(s.Length - s.TotalGapLength));

        // Calculate N50/L50
        var n50Stats = CalculateNx(lengths, totalLength, 50);
        var n90Stats = CalculateNx(lengths, totalLength, 90);

        // GC content
        long gcCount = seqList.Sum(s => CountGC(s.Sequence));
        long totalBases = seqList.Sum(s => CountBases(s.Sequence));
        double gcContent = totalBases > 0 ? gcCount / (double)totalBases : 0;

        // Gap statistics
        int totalGaps = seqList.Sum(s => s.GapCount);
        long totalGapLength = seqList.Sum(s => (long)s.TotalGapLength);
        double gapPercentage = totalLength > 0 ? totalGapLength * 100.0 / totalLength : 0;

        return new AssemblyStatistics(
            seqList.Count,
            totalLength,
            totalNoGaps,
            n50Stats.Nx,
            n50Stats.Lx,
            n90Stats.Nx,
            n90Stats.Lx,
            lengths[0],
            lengths[^1],
            totalLength / (double)seqList.Count,
            lengths[lengths.Count / 2],
            gcContent,
            totalGaps,
            totalGapLength,
            gapPercentage);
    }

    /// <summary>
    /// Parses sequence and extracts metrics.
    /// </summary>
    private static AssembledSequence ParseSequence(string id, string sequence)
    {
        int gapCount = 0;
        int totalGapLength = 0;
        int currentGapLength = 0;

        foreach (char c in sequence)
        {
            if (c == 'N' || c == 'n')
            {
                currentGapLength++;
            }
            else if (currentGapLength > 0)
            {
                gapCount++;
                totalGapLength += currentGapLength;
                currentGapLength = 0;
            }
        }

        if (currentGapLength > 0)
        {
            gapCount++;
            totalGapLength += currentGapLength;
        }

        long gcCount = CountGC(sequence);
        long bases = CountBases(sequence);
        double gc = bases > 0 ? gcCount / (double)bases : 0;

        return new AssembledSequence(id, sequence, sequence.Length, gc, gapCount, totalGapLength);
    }

    /// <summary>
    /// Counts GC bases in sequence.
    /// </summary>
    private static long CountGC(string sequence)
    {
        return sequence.Count(c => c == 'G' || c == 'g' || c == 'C' || c == 'c');
    }

    /// <summary>
    /// Counts non-N bases.
    /// </summary>
    private static long CountBases(string sequence)
    {
        return sequence.Count(c => c != 'N' && c != 'n');
    }

    /// <summary>
    /// Calculates Nx and Lx statistics.
    /// </summary>
    public static NxStatistics CalculateNx(
        IReadOnlyList<int> sortedLengths,
        long totalLength,
        int threshold)
    {
        if (sortedLengths.Count == 0 || totalLength == 0)
            return new NxStatistics(threshold, 0, 0, 0);

        long targetLength = (long)(totalLength * threshold / 100.0);
        long cumulative = 0;
        int count = 0;

        foreach (int length in sortedLengths)
        {
            cumulative += length;
            count++;

            if (cumulative >= targetLength)
            {
                return new NxStatistics(threshold, length, count, cumulative);
            }
        }

        return new NxStatistics(threshold, sortedLengths[^1], count, cumulative);
    }

    /// <summary>
    /// Calculates Nx statistics for multiple thresholds.
    /// </summary>
    public static IEnumerable<NxStatistics> CalculateNxCurve(
        IEnumerable<int> lengths,
        params int[] thresholds)
    {
        var sorted = lengths.OrderByDescending(l => l).ToList();
        long total = sorted.Sum(l => (long)l);

        if (thresholds.Length == 0)
            thresholds = new[] { 10, 20, 30, 40, 50, 60, 70, 80, 90 };

        foreach (int t in thresholds.OrderBy(t => t))
        {
            yield return CalculateNx(sorted, total, t);
        }
    }

    /// <summary>
    /// Calculates auN (area under Nx curve) - assembly contiguity metric.
    /// </summary>
    public static double CalculateAuN(IEnumerable<int> lengths)
    {
        var sorted = lengths.OrderByDescending(l => l).ToList();
        if (sorted.Count == 0) return 0;

        long total = sorted.Sum(l => (long)l);
        if (total == 0) return 0;

        // auN is the sum of (length * length) / total
        double auN = sorted.Sum(l => (double)l * l) / total;
        return auN;
    }

    #endregion

    #region Gap Analysis

    /// <summary>
    /// Identifies all gaps in assembly sequences.
    /// </summary>
    public static IEnumerable<GapInfo> FindGaps(
        IEnumerable<(string Id, string Sequence)> sequences,
        int minGapLength = 1)
    {
        foreach (var (id, sequence) in sequences)
        {
            int gapStart = -1;

            for (int i = 0; i <= sequence.Length; i++)
            {
                bool isN = i < sequence.Length && (sequence[i] == 'N' || sequence[i] == 'n');

                if (isN && gapStart < 0)
                {
                    gapStart = i;
                }
                else if (!isN && gapStart >= 0)
                {
                    int length = i - gapStart;
                    if (length >= minGapLength)
                    {
                        string gapType = ClassifyGap(length);
                        yield return new GapInfo(id, gapStart, i - 1, length, gapType);
                    }
                    gapStart = -1;
                }
            }
        }
    }

    /// <summary>
    /// Classifies gap by length.
    /// </summary>
    private static string ClassifyGap(int length)
    {
        return length switch
        {
            < 10 => "Short",
            < 100 => "Medium",
            < 1000 => "Long",
            _ => "Scaffold"
        };
    }

    /// <summary>
    /// Analyzes gap distribution in assembly.
    /// </summary>
    public static (int Count, double MeanLength, double MedianLength, int MaxLength, IReadOnlyDictionary<string, int> TypeCounts)
        AnalyzeGapDistribution(IEnumerable<GapInfo> gaps)
    {
        var gapList = gaps.ToList();

        if (gapList.Count == 0)
            return (0, 0, 0, 0, new Dictionary<string, int>());

        var lengths = gapList.Select(g => g.Length).OrderBy(l => l).ToList();
        var typeCounts = gapList.GroupBy(g => g.GapType)
                                .ToDictionary(g => g.Key, g => g.Count());

        return (
            gapList.Count,
            lengths.Average(),
            lengths[lengths.Count / 2],
            lengths.Max(),
            typeCounts);
    }

    #endregion

    #region Scaffold Analysis

    /// <summary>
    /// Analyzes scaffold structure (contigs and gaps).
    /// </summary>
    public static IEnumerable<ScaffoldStructure> AnalyzeScaffolds(
        IEnumerable<(string Id, string Sequence)> scaffolds,
        int minGapLength = 10)
    {
        foreach (var (id, sequence) in scaffolds)
        {
            var contigs = new List<(string, int, int)>();
            var gaps = new List<GapInfo>();

            int contigStart = 0;
            int contigId = 1;
            int gapStart = -1;

            for (int i = 0; i <= sequence.Length; i++)
            {
                bool isN = i < sequence.Length && (sequence[i] == 'N' || sequence[i] == 'n');

                if (isN && gapStart < 0)
                {
                    // End of contig, start of gap
                    if (i > contigStart)
                    {
                        contigs.Add(($"{id}_contig{contigId}", contigStart, i - 1));
                        contigId++;
                    }
                    gapStart = i;
                }
                else if (!isN && gapStart >= 0)
                {
                    // End of gap
                    int gapLength = i - gapStart;
                    if (gapLength >= minGapLength)
                    {
                        gaps.Add(new GapInfo(id, gapStart, i - 1, gapLength, ClassifyGap(gapLength)));
                    }
                    contigStart = i;
                    gapStart = -1;
                }
            }

            // Handle final contig
            if (gapStart < 0 && contigStart < sequence.Length)
            {
                contigs.Add(($"{id}_contig{contigId}", contigStart, sequence.Length - 1));
            }

            int totalContigLength = contigs.Sum(c => c.Item3 - c.Item2 + 1);
            int totalGapLen = gaps.Sum(g => g.Length);

            yield return new ScaffoldStructure(
                id, contigs, gaps,
                sequence.Length, totalContigLength, totalGapLen);
        }
    }

    /// <summary>
    /// Extracts contigs from scaffolds.
    /// </summary>
    public static IEnumerable<(string Id, string Sequence)> ExtractContigs(
        IEnumerable<(string Id, string Sequence)> scaffolds,
        int minContigLength = 200)
    {
        foreach (var (id, sequence) in scaffolds)
        {
            var contigs = new List<(int Start, int End)>();
            int contigStart = -1;

            for (int i = 0; i <= sequence.Length; i++)
            {
                bool isN = i < sequence.Length && (sequence[i] == 'N' || sequence[i] == 'n');

                if (!isN && contigStart < 0)
                {
                    contigStart = i;
                }
                else if (isN && contigStart >= 0)
                {
                    if (i - contigStart >= minContigLength)
                    {
                        contigs.Add((contigStart, i));
                    }
                    contigStart = -1;
                }
            }

            if (contigStart >= 0 && sequence.Length - contigStart >= minContigLength)
            {
                contigs.Add((contigStart, sequence.Length));
            }

            int contigNum = 1;
            foreach (var (start, end) in contigs)
            {
                yield return ($"{id}_contig{contigNum++}", sequence[start..end]);
            }
        }
    }

    #endregion

    #region Completeness Assessment

    /// <summary>
    /// Assesses assembly completeness using marker genes (BUSCO-like).
    /// </summary>
    public static CompletenessResult AssessCompleteness(
        IEnumerable<(string Id, string Sequence)> assembly,
        IEnumerable<(string GeneId, string Sequence)> markerGenes,
        double identityThreshold = 0.9,
        double coverageThreshold = 0.9)
    {
        var assemblySeqs = assembly.ToList();
        var markers = markerGenes.ToList();

        if (markers.Count == 0)
            return new CompletenessResult(0, 0, 0, 0, 0, 0, 0, 0);

        int complete = 0;
        int singleCopy = 0;
        int duplicated = 0;
        int fragmented = 0;
        int missing = 0;

        foreach (var (geneId, geneSeq) in markers)
        {
            var hits = FindMarkerHits(geneSeq, assemblySeqs, identityThreshold, coverageThreshold);

            if (hits.Count == 0)
            {
                missing++;
            }
            else if (hits.Any(h => h.Coverage >= coverageThreshold))
            {
                complete++;
                var completeHits = hits.Count(h => h.Coverage >= coverageThreshold);
                if (completeHits == 1)
                    singleCopy++;
                else
                    duplicated++;
            }
            else
            {
                fragmented++;
            }
        }

        double completeness = markers.Count > 0 ? complete * 100.0 / markers.Count : 0;
        double duplication = complete > 0 ? duplicated * 100.0 / complete : 0;

        return new CompletenessResult(
            markers.Count,
            complete,
            singleCopy,
            duplicated,
            fragmented,
            missing,
            completeness,
            duplication);
    }

    /// <summary>
    /// Finds marker gene hits in assembly.
    /// </summary>
    private static List<(string SeqId, double Identity, double Coverage)> FindMarkerHits(
        string markerSeq,
        List<(string Id, string Sequence)> assembly,
        double identityThreshold,
        double coverageThreshold)
    {
        var hits = new List<(string, double, double)>();

        // Simplified k-mer based search
        int kmerSize = Math.Min(31, markerSeq.Length / 3);
        if (kmerSize < 11) kmerSize = 11;

        var markerKmers = new HashSet<string>();
        for (int i = 0; i <= markerSeq.Length - kmerSize; i++)
        {
            markerKmers.Add(markerSeq.Substring(i, kmerSize).ToUpperInvariant());
        }

        foreach (var (seqId, sequence) in assembly)
        {
            int matchingKmers = 0;
            var seqUpper = sequence.ToUpperInvariant();

            for (int i = 0; i <= seqUpper.Length - kmerSize; i++)
            {
                if (markerKmers.Contains(seqUpper.Substring(i, kmerSize)))
                    matchingKmers++;
            }

            double coverage = markerKmers.Count > 0 ?
                matchingKmers / (double)markerKmers.Count : 0;

            if (coverage >= coverageThreshold * 0.5) // Relaxed for fragmented
            {
                // Estimate identity from k-mer matches
                double identity = Math.Min(1.0, coverage * 1.1);
                if (identity >= identityThreshold * 0.8)
                {
                    hits.Add((seqId, identity, coverage));
                }
            }
        }

        return hits;
    }

    /// <summary>
    /// Estimates genome completeness by k-mer spectrum analysis.
    /// </summary>
    public static (double Completeness, double ErrorRate, long EstimatedGenomeSize)
        EstimateCompletenessFromKmers(
            IEnumerable<(string Kmer, int Count)> kmerSpectrum,
            int expectedCoverage = 0)
    {
        var spectrum = kmerSpectrum.ToList();

        if (spectrum.Count == 0)
            return (0, 0, 0);

        // Find coverage peak (mode)
        var countDistribution = spectrum
            .Where(k => k.Count > 1) // Exclude singletons (likely errors)
            .GroupBy(k => k.Count)
            .Select(g => (Count: g.Key, Frequency: g.Count()))
            .OrderByDescending(g => g.Frequency)
            .ToList();

        if (countDistribution.Count == 0)
            return (0, 1.0, 0);

        int peakCoverage = expectedCoverage > 0 ?
            expectedCoverage :
            countDistribution.First().Count;

        // Calculate error rate from singleton ratio
        int singletons = spectrum.Count(k => k.Count == 1);
        double errorRate = singletons / (double)spectrum.Count;

        // Estimate genome size
        long totalKmers = spectrum.Sum(k => (long)k.Count);
        long estimatedSize = peakCoverage > 0 ? totalKmers / peakCoverage : 0;

        // Estimate completeness
        int solidKmers = spectrum.Count(k => k.Count >= peakCoverage / 2);
        double completeness = solidKmers / (double)Math.Max(1, spectrum.Count - singletons);

        return (completeness, errorRate, estimatedSize);
    }

    #endregion

    #region Repeat Analysis

    /// <summary>
    /// Identifies repetitive sequences using k-mer frequency.
    /// </summary>
    public static IEnumerable<(string SequenceId, int Start, int End, int Copies)> FindRepetitiveRegions(
        IEnumerable<(string Id, string Sequence)> sequences,
        int kmerSize = 15,
        int minCopies = 3,
        int windowSize = 100)
    {
        // Build k-mer frequency table
        var kmerCounts = new Dictionary<string, int>();

        foreach (var (id, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();
            for (int i = 0; i <= upper.Length - kmerSize; i++)
            {
                string kmer = upper.Substring(i, kmerSize);
                if (!kmer.Contains('N'))
                {
                    kmerCounts[kmer] = kmerCounts.GetValueOrDefault(kmer) + 1;
                }
            }
        }

        // Find regions with high k-mer frequency
        foreach (var (id, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();
            int regionStart = -1;
            int maxCopies = 0;

            for (int i = 0; i <= upper.Length - kmerSize; i += windowSize / 10)
            {
                // Check window
                int windowEnd = Math.Min(i + windowSize, upper.Length);
                int highCopyKmers = 0;
                int localMaxCopies = 0;

                for (int j = i; j <= windowEnd - kmerSize; j++)
                {
                    string kmer = upper.Substring(j, kmerSize);
                    if (!kmer.Contains('N') && kmerCounts.TryGetValue(kmer, out int count))
                    {
                        if (count >= minCopies)
                        {
                            highCopyKmers++;
                            localMaxCopies = Math.Max(localMaxCopies, count);
                        }
                    }
                }

                bool isRepetitive = highCopyKmers > windowSize / (kmerSize * 2);

                if (isRepetitive && regionStart < 0)
                {
                    regionStart = i;
                    maxCopies = localMaxCopies;
                }
                else if (!isRepetitive && regionStart >= 0)
                {
                    yield return (id, regionStart, i - 1, maxCopies);
                    regionStart = -1;
                    maxCopies = 0;
                }
                else if (isRepetitive)
                {
                    maxCopies = Math.Max(maxCopies, localMaxCopies);
                }
            }

            if (regionStart >= 0)
            {
                yield return (id, regionStart, upper.Length - 1, maxCopies);
            }
        }
    }

    /// <summary>
    /// Identifies tandem repeats.
    /// </summary>
    public static IEnumerable<(string SequenceId, int Start, int End, string Unit, int Copies, double Purity)>
        FindTandemRepeats(
            IEnumerable<(string Id, string Sequence)> sequences,
            int minUnitLength = 2,
            int maxUnitLength = 50,
            int minCopies = 3)
    {
        foreach (var (id, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();

            for (int i = 0; i < upper.Length - minUnitLength * minCopies; i++)
            {
                for (int unitLen = minUnitLength; unitLen <= maxUnitLength && i + unitLen * minCopies <= upper.Length; unitLen++)
                {
                    string unit = upper.Substring(i, unitLen);
                    if (unit.Contains('N')) continue;

                    // Count copies
                    int copies = 1;
                    int matches = 0;
                    int pos = i + unitLen;

                    while (pos + unitLen <= upper.Length)
                    {
                        int unitMatches = 0;
                        for (int j = 0; j < unitLen; j++)
                        {
                            if (upper[pos + j] == unit[j])
                                unitMatches++;
                        }

                        if (unitMatches >= unitLen * 0.8) // 80% match
                        {
                            copies++;
                            matches += unitMatches;
                            pos += unitLen;
                        }
                        else
                        {
                            break;
                        }
                    }

                    if (copies >= minCopies)
                    {
                        int totalLength = copies * unitLen;
                        double purity = (unitLen + matches) / (double)(totalLength);

                        yield return (id, i, i + totalLength - 1, unit, copies, purity);

                        // Skip past this repeat
                        i += totalLength - 1;
                        break;
                    }
                }
            }
        }
    }

    /// <summary>
    /// Calculates repeat content statistics.
    /// </summary>
    public static (long TotalRepeatLength, double RepeatPercentage, IReadOnlyDictionary<string, long> RepeatClassLengths)
        CalculateRepeatContent(
            IEnumerable<RepeatAnnotation> repeats,
            long genomeLength)
    {
        var repeatList = repeats.ToList();

        if (repeatList.Count == 0 || genomeLength == 0)
            return (0, 0, new Dictionary<string, long>());

        long totalLength = repeatList.Sum(r => (long)(r.End - r.Start + 1));
        double percentage = totalLength * 100.0 / genomeLength;

        var classLengths = repeatList
            .GroupBy(r => r.RepeatClass)
            .ToDictionary(g => g.Key, g => g.Sum(r => (long)(r.End - r.Start + 1)));

        return (totalLength, percentage, classLengths);
    }

    #endregion

    #region Assembly Comparison

    /// <summary>
    /// Compares two assemblies using k-mer content.
    /// </summary>
    public static AssemblyComparison CompareAssemblies(
        IEnumerable<(string Id, string Sequence)> assembly1,
        IEnumerable<(string Id, string Sequence)> assembly2,
        string name1 = "Assembly1",
        string name2 = "Assembly2",
        int kmerSize = 21)
    {
        var kmers1 = ExtractKmers(assembly1, kmerSize);
        var kmers2 = ExtractKmers(assembly2, kmerSize);

        int shared = kmers1.Keys.Count(k => kmers2.ContainsKey(k));

        double aligned1 = kmers1.Count > 0 ? shared / (double)kmers1.Count : 0;
        double aligned2 = kmers2.Count > 0 ? shared / (double)kmers2.Count : 0;

        // Estimate identity from shared k-mers
        double identity = (aligned1 + aligned2) / 2;

        return new AssemblyComparison(
            name1, name2,
            aligned1, aligned2,
            0, 0, 0, identity);
    }

    /// <summary>
    /// Extracts k-mers from assembly.
    /// </summary>
    private static Dictionary<string, int> ExtractKmers(
        IEnumerable<(string Id, string Sequence)> sequences,
        int kmerSize)
    {
        var kmers = new Dictionary<string, int>();

        foreach (var (_, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();
            for (int i = 0; i <= upper.Length - kmerSize; i++)
            {
                string kmer = upper.Substring(i, kmerSize);
                if (!kmer.Contains('N'))
                {
                    kmers[kmer] = kmers.GetValueOrDefault(kmer) + 1;
                }
            }
        }

        return kmers;
    }

    /// <summary>
    /// Finds syntenic blocks between assemblies.
    /// </summary>
    public static IEnumerable<(string Seq1, int Start1, int End1, string Seq2, int Start2, int End2, bool IsInverted)>
        FindSyntenicBlocks(
            IEnumerable<(string Id, string Sequence)> assembly1,
            IEnumerable<(string Id, string Sequence)> assembly2,
            int minBlockSize = 1000,
            int kmerSize = 21)
    {
        var asm1 = assembly1.ToList();
        var asm2 = assembly2.ToList();

        // Build k-mer index for assembly2
        var kmerIndex = new Dictionary<string, List<(string SeqId, int Position)>>();

        foreach (var (seqId, sequence) in asm2)
        {
            var upper = sequence.ToUpperInvariant();
            for (int i = 0; i <= upper.Length - kmerSize; i++)
            {
                string kmer = upper.Substring(i, kmerSize);
                if (!kmer.Contains('N'))
                {
                    if (!kmerIndex.ContainsKey(kmer))
                        kmerIndex[kmer] = new List<(string, int)>();
                    kmerIndex[kmer].Add((seqId, i));
                }
            }
        }

        // Find matching regions
        foreach (var (seq1Id, seq1) in asm1)
        {
            var upper1 = seq1.ToUpperInvariant();
            var anchors = new List<(int Pos1, string Seq2, int Pos2)>();

            for (int i = 0; i <= upper1.Length - kmerSize; i += kmerSize)
            {
                string kmer = upper1.Substring(i, kmerSize);
                if (!kmer.Contains('N') && kmerIndex.TryGetValue(kmer, out var hits))
                {
                    foreach (var (seq2Id, pos2) in hits)
                    {
                        anchors.Add((i, seq2Id, pos2));
                    }
                }
            }

            // Cluster anchors into blocks
            var blocks = ClusterAnchors(anchors, minBlockSize);

            foreach (var block in blocks)
            {
                yield return (seq1Id, block.Start1, block.End1,
                             block.Seq2, block.Start2, block.End2, block.IsInverted);
            }
        }
    }

    /// <summary>
    /// Clusters k-mer anchors into syntenic blocks.
    /// </summary>
    private static List<(int Start1, int End1, string Seq2, int Start2, int End2, bool IsInverted)>
        ClusterAnchors(
            List<(int Pos1, string Seq2, int Pos2)> anchors,
            int minBlockSize)
    {
        var result = new List<(int, int, string, int, int, bool)>();

        if (anchors.Count == 0)
            return result;

        // Group by target sequence
        var bySeq2 = anchors.GroupBy(a => a.Seq2);

        foreach (var group in bySeq2)
        {
            var sorted = group.OrderBy(a => a.Pos1).ToList();

            int blockStart1 = sorted[0].Pos1;
            int blockStart2 = sorted[0].Pos2;
            int lastPos1 = sorted[0].Pos1;
            int lastPos2 = sorted[0].Pos2;
            bool? isInverted = null;

            for (int i = 1; i < sorted.Count; i++)
            {
                int gap1 = sorted[i].Pos1 - lastPos1;
                int gap2 = sorted[i].Pos2 - lastPos2;

                bool currentInverted = gap2 < 0;

                if (isInverted == null)
                    isInverted = currentInverted;

                // Check if continuous
                if (gap1 > minBlockSize || Math.Abs(gap2) > minBlockSize || currentInverted != isInverted)
                {
                    // End current block
                    if (lastPos1 - blockStart1 >= minBlockSize)
                    {
                        result.Add((blockStart1, lastPos1, group.Key,
                                   Math.Min(blockStart2, lastPos2),
                                   Math.Max(blockStart2, lastPos2),
                                   isInverted ?? false));
                    }

                    blockStart1 = sorted[i].Pos1;
                    blockStart2 = sorted[i].Pos2;
                    isInverted = null;
                }

                lastPos1 = sorted[i].Pos1;
                lastPos2 = sorted[i].Pos2;
            }

            // Handle last block
            if (lastPos1 - blockStart1 >= minBlockSize)
            {
                result.Add((blockStart1, lastPos1, group.Key,
                           Math.Min(blockStart2, lastPos2),
                           Math.Max(blockStart2, lastPos2),
                           isInverted ?? false));
            }
        }

        return result;
    }

    #endregion

    #region Quality Assessment

    /// <summary>
    /// Calculates per-base quality metrics.
    /// </summary>
    public static IEnumerable<(string SequenceId, int Position, int WindowSize, double GcContent, int NCount, double Complexity)>
        CalculateLocalQuality(
            IEnumerable<(string Id, string Sequence)> sequences,
            int windowSize = 1000)
    {
        foreach (var (id, sequence) in sequences)
        {
            var upper = sequence.ToUpperInvariant();

            for (int i = 0; i < upper.Length; i += windowSize / 2)
            {
                int end = Math.Min(i + windowSize, upper.Length);
                string window = upper[i..end];

                int gc = window.Count(c => c == 'G' || c == 'C');
                int nCount = window.Count(c => c == 'N');
                int bases = window.Length - nCount;

                double gcContent = bases > 0 ? gc / (double)bases : 0;
                double complexity = CalculateLinguisticComplexity(window);

                yield return (id, i, end - i, gcContent, nCount, complexity);
            }
        }
    }

    /// <summary>
    /// Calculates linguistic complexity (ratio of distinct k-mers to possible k-mers).
    /// </summary>
    private static double CalculateLinguisticComplexity(string sequence, int kmerSize = 4)
    {
        if (sequence.Length < kmerSize)
            return 0;

        var kmers = new HashSet<string>();
        int validKmers = 0;

        for (int i = 0; i <= sequence.Length - kmerSize; i++)
        {
            string kmer = sequence.Substring(i, kmerSize);
            if (!kmer.Contains('N'))
            {
                kmers.Add(kmer);
                validKmers++;
            }
        }

        if (validKmers == 0)
            return 0;

        // Maximum possible k-mers for DNA
        int maxKmers = (int)Math.Pow(4, kmerSize);
        int possibleKmers = Math.Min(validKmers, maxKmers);

        return kmers.Count / (double)possibleKmers;
    }

    /// <summary>
    /// Identifies potentially misassembled regions.
    /// </summary>
    public static IEnumerable<(string SequenceId, int Start, int End, string Reason, double Score)>
        FindSuspiciousRegions(
            IEnumerable<(string Id, string Sequence)> sequences,
            double gcDeviation = 0.15,
            double minComplexity = 0.3)
    {
        // First pass: calculate global GC
        var seqList = sequences.ToList();
        long totalGc = 0;
        long totalBases = 0;

        foreach (var (_, sequence) in seqList)
        {
            totalGc += sequence.Count(c => c == 'G' || c == 'g' || c == 'C' || c == 'c');
            totalBases += sequence.Count(c => c != 'N' && c != 'n');
        }

        double globalGc = totalBases > 0 ? totalGc / (double)totalBases : 0.5;

        // Second pass: find anomalies
        foreach (var (id, sequence) in seqList)
        {
            var quality = CalculateLocalQuality(
                new[] { (id, sequence) }, windowSize: 500).ToList();

            int regionStart = -1;
            string? reason = null;
            double maxScore = 0;

            foreach (var (_, pos, size, gc, nCount, complexity) in quality)
            {
                bool suspicious = false;
                double score = 0;
                string? currentReason = null;

                // High N content
                if (nCount > size * 0.1)
                {
                    suspicious = true;
                    score = Math.Max(score, nCount / (double)size);
                    currentReason = "High N content";
                }

                // GC deviation
                if (Math.Abs(gc - globalGc) > gcDeviation)
                {
                    suspicious = true;
                    score = Math.Max(score, Math.Abs(gc - globalGc) / gcDeviation);
                    currentReason = currentReason == null ? "GC deviation" : currentReason + "; GC deviation";
                }

                // Low complexity
                if (complexity < minComplexity && nCount < size * 0.5)
                {
                    suspicious = true;
                    score = Math.Max(score, 1 - complexity / minComplexity);
                    currentReason = currentReason == null ? "Low complexity" : currentReason + "; Low complexity";
                }

                if (suspicious && regionStart < 0)
                {
                    regionStart = pos;
                    reason = currentReason;
                    maxScore = score;
                }
                else if (!suspicious && regionStart >= 0)
                {
                    yield return (id, regionStart, pos - 1, reason!, maxScore);
                    regionStart = -1;
                    reason = null;
                    maxScore = 0;
                }
                else if (suspicious)
                {
                    maxScore = Math.Max(maxScore, score);
                    if (currentReason != null && !reason!.Contains(currentReason))
                        reason += "; " + currentReason;
                }
            }

            if (regionStart >= 0)
            {
                yield return (id, regionStart, sequence.Length - 1, reason!, maxScore);
            }
        }
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Filters sequences by length.
    /// </summary>
    public static IEnumerable<(string Id, string Sequence)> FilterByLength(
        IEnumerable<(string Id, string Sequence)> sequences,
        int minLength = 0,
        int maxLength = int.MaxValue)
    {
        return sequences.Where(s => s.Sequence.Length >= minLength && s.Sequence.Length <= maxLength);
    }

    /// <summary>
    /// Sorts sequences by length descending.
    /// </summary>
    public static IEnumerable<(string Id, string Sequence)> SortByLength(
        IEnumerable<(string Id, string Sequence)> sequences,
        bool descending = true)
    {
        return descending ?
            sequences.OrderByDescending(s => s.Sequence.Length) :
            sequences.OrderBy(s => s.Sequence.Length);
    }

    /// <summary>
    /// Calculates sequence length distribution.
    /// </summary>
    public static IReadOnlyDictionary<string, int> CalculateLengthDistribution(
        IEnumerable<int> lengths,
        params int[] bins)
    {
        if (bins.Length == 0)
            bins = new[] { 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000 };

        var distribution = bins.ToDictionary(b => $"<{b}", _ => 0);
        distribution[$">={bins[^1]}"] = 0;

        foreach (int length in lengths)
        {
            bool found = false;
            foreach (int bin in bins)
            {
                if (length < bin)
                {
                    distribution[$"<{bin}"]++;
                    found = true;
                    break;
                }
            }
            if (!found)
                distribution[$">={bins[^1]}"]++;
        }

        return distribution;
    }

    #endregion
}
