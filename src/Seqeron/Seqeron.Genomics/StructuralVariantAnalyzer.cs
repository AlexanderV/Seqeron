using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics;

/// <summary>
/// Provides algorithms for structural variant detection and analysis.
/// </summary>
public static class StructuralVariantAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Types of structural variants.
    /// </summary>
    public enum SVType
    {
        Deletion,
        Duplication,
        Inversion,
        Insertion,
        Translocation,
        ComplexRearrangement,
        CopyNumberVariation
    }

    /// <summary>
    /// Represents a structural variant.
    /// </summary>
    public readonly record struct StructuralVariant(
        string Id,
        string Chromosome,
        int Start,
        int End,
        SVType Type,
        int Length,
        double Quality,
        int SupportingReads,
        string? InsertedSequence);

    /// <summary>
    /// Represents a breakpoint.
    /// </summary>
    public readonly record struct Breakpoint(
        string Chromosome1,
        int Position1,
        char Strand1,
        string Chromosome2,
        int Position2,
        char Strand2,
        int SupportingReads,
        double Quality);

    /// <summary>
    /// Represents copy number segment.
    /// </summary>
    public readonly record struct CopyNumberSegment(
        string Chromosome,
        int Start,
        int End,
        double LogRatio,
        int CopyNumber,
        double BAlleleFrequency,
        int ProbeCount);

    /// <summary>
    /// Represents a read pair signature.
    /// </summary>
    public readonly record struct ReadPairSignature(
        string ReadId,
        string Chromosome1,
        int Position1,
        char Strand1,
        string Chromosome2,
        int Position2,
        char Strand2,
        int InsertSize,
        bool IsDiscordant);

    /// <summary>
    /// Represents a split read.
    /// </summary>
    public readonly record struct SplitRead(
        string ReadId,
        string Chromosome,
        int PrimaryPosition,
        int SupplementaryPosition,
        int ClipLength,
        string ClippedSequence);

    /// <summary>
    /// Represents SV annotation.
    /// </summary>
    public readonly record struct SVAnnotation(
        string SVId,
        IReadOnlyList<string> AffectedGenes,
        IReadOnlyList<string> AffectedExons,
        string FunctionalImpact,
        double PopulationFrequency,
        bool IsPathogenic);

    #endregion

    #region Read Pair Analysis

    /// <summary>
    /// Identifies discordant read pairs.
    /// </summary>
    public static IEnumerable<ReadPairSignature> FindDiscordantPairs(
        IEnumerable<(string ReadId, string Chr1, int Pos1, char Strand1, string Chr2, int Pos2, char Strand2, int InsertSize)> readPairs,
        int expectedInsertSize = 400,
        int insertSizeStdDev = 50,
        int maxInsertSize = 10000)
    {
        double lowerBound = expectedInsertSize - 3 * insertSizeStdDev;
        double upperBound = expectedInsertSize + 3 * insertSizeStdDev;

        foreach (var (readId, chr1, pos1, strand1, chr2, pos2, strand2, insertSize) in readPairs)
        {
            bool isDiscordant = false;
            string reason = "";

            // Check for interchromosomal
            if (chr1 != chr2)
            {
                isDiscordant = true;
                reason = "Interchromosomal";
            }
            // Check insert size
            else if (insertSize < lowerBound || insertSize > upperBound)
            {
                isDiscordant = true;
                reason = insertSize > upperBound ? "LargeInsert" : "SmallInsert";
            }
            // Check orientation (expecting FR for standard library)
            else if (!((strand1 == '+' && strand2 == '-') || (strand1 == '-' && strand2 == '+')))
            {
                isDiscordant = true;
                reason = "AbnormalOrientation";
            }

            if (isDiscordant || insertSize > maxInsertSize)
            {
                yield return new ReadPairSignature(
                    ReadId: readId,
                    Chromosome1: chr1,
                    Position1: pos1,
                    Strand1: strand1,
                    Chromosome2: chr2,
                    Position2: pos2,
                    Strand2: strand2,
                    InsertSize: insertSize,
                    IsDiscordant: true);
            }
        }
    }

    /// <summary>
    /// Clusters discordant read pairs into SV candidates.
    /// </summary>
    public static IEnumerable<StructuralVariant> ClusterDiscordantPairs(
        IEnumerable<ReadPairSignature> discordantPairs,
        int clusterDistance = 500,
        int minSupport = 3)
    {
        var pairs = discordantPairs.OrderBy(p => p.Chromosome1).ThenBy(p => p.Position1).ToList();

        if (pairs.Count == 0)
            yield break;

        var clusters = new List<List<ReadPairSignature>>();
        var currentCluster = new List<ReadPairSignature> { pairs[0] };

        for (int i = 1; i < pairs.Count; i++)
        {
            var prev = pairs[i - 1];
            var curr = pairs[i];

            bool sameCluster = prev.Chromosome1 == curr.Chromosome1 &&
                               prev.Chromosome2 == curr.Chromosome2 &&
                               Math.Abs(curr.Position1 - prev.Position1) <= clusterDistance &&
                               Math.Abs(curr.Position2 - prev.Position2) <= clusterDistance;

            if (sameCluster)
            {
                currentCluster.Add(curr);
            }
            else
            {
                if (currentCluster.Count >= minSupport)
                {
                    clusters.Add(currentCluster);
                }
                currentCluster = new List<ReadPairSignature> { curr };
            }
        }

        if (currentCluster.Count >= minSupport)
        {
            clusters.Add(currentCluster);
        }

        int svId = 1;
        foreach (var cluster in clusters)
        {
            var sv = CreateSVFromCluster(cluster, svId++);
            if (sv != null)
            {
                yield return sv.Value;
            }
        }
    }

    private static StructuralVariant? CreateSVFromCluster(List<ReadPairSignature> cluster, int id)
    {
        if (cluster.Count == 0)
            return null;

        var first = cluster[0];
        int start = cluster.Min(p => p.Position1);
        int end = cluster.Max(p => p.Position2);

        SVType type;
        if (first.Chromosome1 != first.Chromosome2)
        {
            type = SVType.Translocation;
        }
        else if (first.Strand1 == first.Strand2)
        {
            type = SVType.Inversion;
        }
        else if (end - start > 10000)
        {
            type = SVType.Deletion;
        }
        else
        {
            type = SVType.Duplication;
        }

        return new StructuralVariant(
            Id: $"SV{id}",
            Chromosome: first.Chromosome1,
            Start: start,
            End: end,
            Type: type,
            Length: Math.Abs(end - start),
            Quality: Math.Min(cluster.Count * 10.0, 100.0),
            SupportingReads: cluster.Count,
            InsertedSequence: null);
    }

    #endregion

    #region Split Read Analysis

    /// <summary>
    /// Identifies split reads from soft-clipped alignments.
    /// </summary>
    public static IEnumerable<SplitRead> FindSplitReads(
        IEnumerable<(string ReadId, string Chromosome, int Position, string Cigar, string Sequence)> alignments,
        int minClipLength = 20)
    {
        foreach (var (readId, chromosome, position, cigar, sequence) in alignments)
        {
            var clips = ParseSoftClips(cigar);

            foreach (var (clipPos, clipLen, isLeft) in clips)
            {
                if (clipLen >= minClipLength)
                {
                    int clipStart = isLeft ? 0 : sequence.Length - clipLen;
                    string clippedSeq = sequence.Substring(clipStart, clipLen);

                    int suppPos = isLeft ? position : position + GetAlignedLength(cigar);

                    yield return new SplitRead(
                        ReadId: readId,
                        Chromosome: chromosome,
                        PrimaryPosition: position,
                        SupplementaryPosition: suppPos,
                        ClipLength: clipLen,
                        ClippedSequence: clippedSeq);
                }
            }
        }
    }

    private static List<(int Position, int Length, bool IsLeft)> ParseSoftClips(string cigar)
    {
        var clips = new List<(int, int, bool)>();
        int pos = 0;
        int numStart = 0;

        for (int i = 0; i < cigar.Length; i++)
        {
            if (char.IsDigit(cigar[i]))
            {
                if (numStart < 0)
                    numStart = i;
            }
            else
            {
                if (numStart >= 0)
                {
                    int len = int.Parse(cigar.Substring(numStart, i - numStart));

                    if (cigar[i] == 'S')
                    {
                        clips.Add((pos, len, pos == 0));
                    }
                    else if (cigar[i] == 'M' || cigar[i] == 'D' || cigar[i] == 'N')
                    {
                        pos += len;
                    }
                }
                numStart = -1;
            }
        }

        return clips;
    }

    private static int GetAlignedLength(string cigar)
    {
        int length = 0;
        int numStart = 0;

        for (int i = 0; i < cigar.Length; i++)
        {
            if (char.IsDigit(cigar[i]))
            {
                if (numStart == 0 || !char.IsDigit(cigar[numStart]))
                    numStart = i;
            }
            else
            {
                if (i > numStart && char.IsDigit(cigar[numStart]))
                {
                    int len = int.Parse(cigar.Substring(numStart, i - numStart));
                    if (cigar[i] == 'M' || cigar[i] == 'D' || cigar[i] == 'N')
                    {
                        length += len;
                    }
                }
                numStart = i + 1;
            }
        }

        return length;
    }

    /// <summary>
    /// Clusters split reads to identify breakpoints.
    /// </summary>
    public static IEnumerable<Breakpoint> ClusterSplitReads(
        IEnumerable<SplitRead> splitReads,
        int clusterDistance = 10,
        int minSupport = 2)
    {
        var reads = splitReads.OrderBy(r => r.Chromosome).ThenBy(r => r.PrimaryPosition).ToList();

        if (reads.Count == 0)
            yield break;

        var clusters = new List<List<SplitRead>>();
        var currentCluster = new List<SplitRead> { reads[0] };

        for (int i = 1; i < reads.Count; i++)
        {
            var prev = reads[i - 1];
            var curr = reads[i];

            bool sameCluster = prev.Chromosome == curr.Chromosome &&
                               Math.Abs(curr.PrimaryPosition - prev.PrimaryPosition) <= clusterDistance;

            if (sameCluster)
            {
                currentCluster.Add(curr);
            }
            else
            {
                if (currentCluster.Count >= minSupport)
                {
                    clusters.Add(currentCluster);
                }
                currentCluster = new List<SplitRead> { curr };
            }
        }

        if (currentCluster.Count >= minSupport)
        {
            clusters.Add(currentCluster);
        }

        foreach (var cluster in clusters)
        {
            int pos = (int)cluster.Average(r => r.PrimaryPosition);
            int suppPos = (int)cluster.Average(r => r.SupplementaryPosition);

            yield return new Breakpoint(
                Chromosome1: cluster[0].Chromosome,
                Position1: pos,
                Strand1: '+',
                Chromosome2: cluster[0].Chromosome,
                Position2: suppPos,
                Strand2: '-',
                SupportingReads: cluster.Count,
                Quality: Math.Min(cluster.Count * 15.0, 100.0));
        }
    }

    #endregion

    #region Copy Number Analysis

    /// <summary>
    /// Segments copy number data using circular binary segmentation-like approach.
    /// </summary>
    public static IEnumerable<CopyNumberSegment> SegmentCopyNumber(
        IEnumerable<(string Chromosome, int Position, double LogRatio, double BAF)> probes,
        double changeThreshold = 0.3,
        int minProbes = 5)
    {
        var probeList = probes.OrderBy(p => p.Chromosome).ThenBy(p => p.Position).ToList();

        if (probeList.Count == 0)
            yield break;

        var currentChrom = probeList[0].Chromosome;
        int segStart = probeList[0].Position;
        var segmentProbes = new List<(int Position, double LogRatio, double BAF)>
        {
            (probeList[0].Position, probeList[0].LogRatio, probeList[0].BAF)
        };

        for (int i = 1; i < probeList.Count; i++)
        {
            var probe = probeList[i];

            bool newSegment = probe.Chromosome != currentChrom;

            if (!newSegment && segmentProbes.Count >= minProbes)
            {
                double currentMean = segmentProbes.Average(p => p.LogRatio);
                if (Math.Abs(probe.LogRatio - currentMean) > changeThreshold)
                {
                    newSegment = true;
                }
            }

            if (newSegment)
            {
                if (segmentProbes.Count >= minProbes)
                {
                    yield return CreateSegment(currentChrom, segStart, segmentProbes);
                }

                currentChrom = probe.Chromosome;
                segStart = probe.Position;
                segmentProbes.Clear();
            }

            segmentProbes.Add((probe.Position, probe.LogRatio, probe.BAF));
        }

        if (segmentProbes.Count >= minProbes)
        {
            yield return CreateSegment(currentChrom, segStart, segmentProbes);
        }
    }

    private static CopyNumberSegment CreateSegment(
        string chromosome,
        int start,
        List<(int Position, double LogRatio, double BAF)> probes)
    {
        double meanLogR = probes.Average(p => p.LogRatio);
        double meanBAF = probes.Average(p => p.BAF);

        // Estimate copy number from log ratio (assuming diploid baseline)
        int copyNumber = (int)Math.Round(2 * Math.Pow(2, meanLogR));
        copyNumber = Math.Max(0, Math.Min(copyNumber, 10));

        return new CopyNumberSegment(
            Chromosome: chromosome,
            Start: start,
            End: probes.Last().Position,
            LogRatio: meanLogR,
            CopyNumber: copyNumber,
            BAlleleFrequency: meanBAF,
            ProbeCount: probes.Count);
    }

    /// <summary>
    /// Identifies copy number variants from segments.
    /// </summary>
    public static IEnumerable<StructuralVariant> IdentifyCNVs(
        IEnumerable<CopyNumberSegment> segments,
        int normalCopyNumber = 2,
        int minLength = 10000)
    {
        int id = 1;

        foreach (var segment in segments)
        {
            if (segment.CopyNumber == normalCopyNumber)
                continue;

            if (segment.End - segment.Start < minLength)
                continue;

            SVType type = segment.CopyNumber < normalCopyNumber
                ? SVType.Deletion
                : SVType.Duplication;

            yield return new StructuralVariant(
                Id: $"CNV{id++}",
                Chromosome: segment.Chromosome,
                Start: segment.Start,
                End: segment.End,
                Type: type,
                Length: segment.End - segment.Start,
                Quality: Math.Abs(segment.LogRatio) * 50,
                SupportingReads: segment.ProbeCount,
                InsertedSequence: null);
        }
    }

    #endregion

    #region SV Merging and Filtering

    /// <summary>
    /// Merges overlapping structural variants.
    /// </summary>
    public static IEnumerable<StructuralVariant> MergeOverlappingSVs(
        IEnumerable<StructuralVariant> variants,
        double overlapFraction = 0.5)
    {
        var svList = variants.OrderBy(v => v.Chromosome).ThenBy(v => v.Start).ToList();

        if (svList.Count == 0)
            yield break;

        var merged = new List<StructuralVariant>();
        var current = svList[0];

        for (int i = 1; i < svList.Count; i++)
        {
            var next = svList[i];

            if (current.Chromosome == next.Chromosome &&
                current.Type == next.Type &&
                CalculateOverlap(current, next) >= overlapFraction)
            {
                // Merge
                current = new StructuralVariant(
                    Id: current.Id,
                    Chromosome: current.Chromosome,
                    Start: Math.Min(current.Start, next.Start),
                    End: Math.Max(current.End, next.End),
                    Type: current.Type,
                    Length: Math.Max(current.End, next.End) - Math.Min(current.Start, next.Start),
                    Quality: Math.Max(current.Quality, next.Quality),
                    SupportingReads: current.SupportingReads + next.SupportingReads,
                    InsertedSequence: current.InsertedSequence ?? next.InsertedSequence);
            }
            else
            {
                yield return current;
                current = next;
            }
        }

        yield return current;
    }

    private static double CalculateOverlap(StructuralVariant sv1, StructuralVariant sv2)
    {
        int overlapStart = Math.Max(sv1.Start, sv2.Start);
        int overlapEnd = Math.Min(sv1.End, sv2.End);

        if (overlapEnd <= overlapStart)
            return 0;

        int overlapLen = overlapEnd - overlapStart;
        int minLen = Math.Min(sv1.Length, sv2.Length);

        return minLen > 0 ? (double)overlapLen / minLen : 0;
    }

    /// <summary>
    /// Filters structural variants by quality and support.
    /// </summary>
    public static IEnumerable<StructuralVariant> FilterSVs(
        IEnumerable<StructuralVariant> variants,
        double minQuality = 20,
        int minSupport = 2,
        int minLength = 50,
        int maxLength = 100_000_000)
    {
        return variants.Where(v =>
            v.Quality >= minQuality &&
            v.SupportingReads >= minSupport &&
            v.Length >= minLength &&
            v.Length <= maxLength);
    }

    #endregion

    #region SV Annotation

    /// <summary>
    /// Annotates structural variants with gene information.
    /// </summary>
    public static IEnumerable<SVAnnotation> AnnotateSVs(
        IEnumerable<StructuralVariant> variants,
        IEnumerable<(string GeneId, string Chromosome, int Start, int End, IReadOnlyList<(int Start, int End)> Exons)> genes)
    {
        var geneList = genes.ToList();

        foreach (var sv in variants)
        {
            var affectedGenes = new List<string>();
            var affectedExons = new List<string>();

            foreach (var gene in geneList.Where(g => g.Chromosome == sv.Chromosome))
            {
                // Check if SV overlaps gene
                if (sv.Start <= gene.End && sv.End >= gene.Start)
                {
                    affectedGenes.Add(gene.GeneId);

                    // Check affected exons
                    for (int i = 0; i < gene.Exons.Count; i++)
                    {
                        var exon = gene.Exons[i];
                        if (sv.Start <= exon.End && sv.End >= exon.Start)
                        {
                            affectedExons.Add($"{gene.GeneId}:exon{i + 1}");
                        }
                    }
                }
            }

            string impact = DetermineImpact(sv, affectedGenes.Count, affectedExons.Count);

            yield return new SVAnnotation(
                SVId: sv.Id,
                AffectedGenes: affectedGenes,
                AffectedExons: affectedExons,
                FunctionalImpact: impact,
                PopulationFrequency: 0, // Would need population database
                IsPathogenic: impact == "HIGH" || impact == "MODERATE");
        }
    }

    private static string DetermineImpact(StructuralVariant sv, int geneCount, int exonCount)
    {
        if (exonCount > 0)
        {
            return sv.Type switch
            {
                SVType.Deletion => "HIGH",
                SVType.Duplication => "MODERATE",
                SVType.Inversion => "HIGH",
                SVType.Translocation => "HIGH",
                _ => "MODERATE"
            };
        }

        if (geneCount > 0)
        {
            return "MODIFIER";
        }

        return "LOW";
    }

    #endregion

    #region SV Genotyping

    /// <summary>
    /// Genotypes a structural variant in a sample.
    /// </summary>
    public static (string Genotype, double Quality) GenotypeSV(
        StructuralVariant sv,
        int refReads,
        int altReads,
        int totalReads)
    {
        if (totalReads == 0)
            return ("./.", 0);

        double altFraction = (double)altReads / totalReads;
        double refFraction = (double)refReads / totalReads;

        string genotype;
        double quality;

        if (altFraction < 0.1)
        {
            genotype = "0/0"; // Homozygous reference
            quality = refReads * 3.0;
        }
        else if (altFraction > 0.9)
        {
            genotype = "1/1"; // Homozygous alternate
            quality = altReads * 3.0;
        }
        else if (altFraction >= 0.3 && altFraction <= 0.7)
        {
            genotype = "0/1"; // Heterozygous
            quality = (refReads + altReads) * 2.0;
        }
        else
        {
            genotype = "0/1"; // Likely heterozygous
            quality = (refReads + altReads) * 1.5;
        }

        return (genotype, Math.Min(quality, 99));
    }

    #endregion

    #region Breakpoint Assembly

    /// <summary>
    /// Assembles breakpoint junction sequence from split reads.
    /// </summary>
    public static string? AssembleBreakpointSequence(
        IEnumerable<SplitRead> splitReads,
        int minOverlap = 10)
    {
        var reads = splitReads.OrderBy(r => r.ClipLength).ToList();

        if (reads.Count == 0)
            return null;

        // Simple assembly: use longest clipped sequence
        var longest = reads.MaxBy(r => r.ClipLength);
        return longest.ClippedSequence;
    }

    /// <summary>
    /// Identifies microhomology at breakpoint junctions.
    /// </summary>
    public static (int MicrohomologyLength, string Sequence) FindMicrohomology(
        string leftFlank,
        string rightFlank,
        int maxLength = 20)
    {
        if (string.IsNullOrEmpty(leftFlank) || string.IsNullOrEmpty(rightFlank))
            return (0, "");

        leftFlank = leftFlank.ToUpperInvariant();
        rightFlank = rightFlank.ToUpperInvariant();

        int maxMH = Math.Min(maxLength, Math.Min(leftFlank.Length, rightFlank.Length));

        for (int len = maxMH; len >= 1; len--)
        {
            string leftEnd = leftFlank.Substring(leftFlank.Length - len);
            string rightStart = rightFlank.Substring(0, len);

            if (leftEnd == rightStart)
            {
                return (len, leftEnd);
            }
        }

        return (0, "");
    }

    #endregion
}
