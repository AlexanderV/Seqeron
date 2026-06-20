using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Provides algorithms for structural variant detection and analysis.
/// </summary>
public static class StructuralVariantAnalyzer
{
    #region Constants

    // Default anomaly cutoff in units of insert-size standard deviation. A read pair is flagged
    // discordant-by-span when its insert size falls outside mean ± cutoff·sd.
    // Source: BreakDancer README (-c, default 3); Fan et al. 2014 (Curr Protoc Bioinformatics) "3 s.d.".
    private const double DefaultInsertSizeCutoffSd = 3.0;

    // Default minimum number of supporting read pairs required to report an SV.
    // Source: BreakDancer README (-r minimum number of read pairs to establish a connection, default 2).
    private const int DefaultMinSupport = 2;

    // Default position tolerance (bases) within which split-read junctions are clustered into one
    // breakpoint. Source: ClipCrop (Suzuki et al. 2011) — breakpoints are "clustered within 5-base
    // differences"; mapping imprecision spreads a single breakpoint across nearby positions.
    private const int DefaultBreakpointClusterTolerance = 5;

    // Default minimum number of supporting split reads required to report a breakpoint.
    // Source: SoftSearch (Hart et al. 2013) — "A putative breakpoint is defined when there is at
    // least x soft-clipped reads beginning at position y"; the configurable minimum, set here to 2
    // to match the sibling read-pair support default (BreakDancer -r = 2).
    private const int DefaultBreakpointMinSupport = 2;

    // Diploid copy-number baseline (reference ploidy). CNVkit converts a log2 ratio to absolute copy
    // number as CN = ploidy · 2^log2; for a diploid genome "the absolute copy number is calculated as
    // 2 * 2^(log2 value)". Source: CNVkit calling docs / cnvlib/call.py `_log2_ratio_to_absolute_pure`.
    private const int DiploidPloidy = 2;

    // Default read-depth window size (positions per window). Read depth is summarised over
    // non-overlapping fixed-size windows. Yoon et al. (2009) used 100-bp windows; the size is a
    // parameter and 100 is the literature default. Source: Yoon et al. 2009, Genome Research 19(9).
    private const int DefaultCnvWindowSize = 100;

    #endregion

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
    /// Identifies discordant (anomalous) read pairs from paired-end mappings.
    /// </summary>
    /// <remarks>
    /// A pair is anomalous when it shows "unexpected separation distances or orientation"
    /// (BreakDancer, Chen et al. 2009): interchromosomal mapping, an insert size outside
    /// mean ± <paramref name="cutoffSd"/>·sd, or an orientation that is neither FR nor RF.
    /// </remarks>
    /// <param name="readPairs">Mapped read pairs (mate coordinates, strands, insert size).</param>
    /// <param name="expectedInsertSize">Mean library insert size.</param>
    /// <param name="insertSizeStdDev">Library insert-size standard deviation.</param>
    /// <param name="cutoffSd">Anomaly cutoff in units of standard deviation (BreakDancer -c, default 3).</param>
    /// <param name="maxInsertSize">Hard upper bound above which a pair is always anomalous.</param>
    public static IEnumerable<ReadPairSignature> FindDiscordantPairs(
        IEnumerable<(string ReadId, string Chr1, int Pos1, char Strand1, string Chr2, int Pos2, char Strand2, int InsertSize)> readPairs,
        int expectedInsertSize = 400,
        int insertSizeStdDev = 50,
        double cutoffSd = DefaultInsertSizeCutoffSd,
        int maxInsertSize = 10000)
    {
        ArgumentNullException.ThrowIfNull(readPairs);

        // Anomaly bounds: a pair is discordant-by-span when its insert size is strictly
        // outside [mean − c·sd, mean + c·sd] (BreakDancer README: bounds = mean ± c·std).
        double lowerBound = expectedInsertSize - cutoffSd * insertSizeStdDev;
        double upperBound = expectedInsertSize + cutoffSd * insertSizeStdDev;

        foreach (var (readId, chr1, pos1, strand1, chr2, pos2, strand2, insertSize) in readPairs)
        {
            bool isDiscordant =
                // Interchromosomal mapping is a linking/translocation signature (Medvedev et al. 2009).
                chr1 != chr2
                // Span outside the cutoff is a deletion (larger) or insertion (smaller) signature.
                || insertSize < lowerBound
                || insertSize > upperBound
                // Concordant orientation is FR or RF (mates pointing inward); anything else is abnormal.
                || !IsConcordantOrientation(strand1, strand2);

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
    /// Returns true when two mates have the concordant forward-reverse (FR, inward-facing)
    /// orientation for a standard short-insert library: the upstream mate on the '+' strand and the
    /// downstream mate on the '−' strand, so the reads point towards one another.
    /// </summary>
    /// <remarks>
    /// Standard Illumina paired-end libraries yield FR pairs (SAM proper-pair FLAG 0x02). Every other
    /// orientation is a discordant signature: FF/RR (same strand) supports an inversion, and RF
    /// (reverse-forward, outward-facing / "everted") supports a tandem duplication — DELLY, LUMPY,
    /// Manta and SVXplorer all read an FR cluster as a deletion candidate and an RF cluster as a
    /// duplication candidate. RF is "proper" only for opposite-orientation mate-pair libraries, not for
    /// the short-insert library modelled here. Source: cureffi.org / BWA proper-pair convention
    /// ("RF, FF or RR … that's a problem"); Rausch et al. 2012 (DELLY); SVXplorer (Kumar et al. 2020).
    /// </remarks>
    private static bool IsConcordantOrientation(char strand1, char strand2) =>
        strand1 == '+' && strand2 == '-';

    /// <summary>
    /// Classifies a discordant read-pair signature into a structural-variant type using the
    /// paired-end mapping (PEM) signatures of Medvedev, Stanciu &amp; Brudno (2009).
    /// </summary>
    /// <remarks>
    /// Decision order (source-traceable):
    /// <list type="number">
    /// <item>mates on different chromosomes → <see cref="SVType.Translocation"/> (linking/CTX signature; ASSUMPTION A1: chromosome difference takes precedence over orientation);</item>
    /// <item>same chromosome, mates on the same strand → <see cref="SVType.Inversion"/> (flipped orientation);</item>
    /// <item>same chromosome, reverse-forward (RF, outward-facing / everted) mates → <see cref="SVType.Duplication"/> (tandem-duplication signature; DELLY, LUMPY, Manta, SVXplorer read an RF cluster as a duplication);</item>
    /// <item>same chromosome, forward-reverse (FR) mates, span &gt; mean + c·sd → <see cref="SVType.Deletion"/> (span larger than insert size);</item>
    /// <item>same chromosome, forward-reverse (FR) mates, span &lt; mean − c·sd → <see cref="SVType.Insertion"/> (span smaller than insert size);</item>
    /// <item>otherwise → <see cref="SVType.ComplexRearrangement"/> (anomalous but not matching a basic signature).</item>
    /// </list>
    /// </remarks>
    /// <param name="pair">A read-pair signature (typically one flagged by <see cref="FindDiscordantPairs"/>).</param>
    /// <param name="expectedInsertSize">Mean library insert size.</param>
    /// <param name="insertSizeStdDev">Library insert-size standard deviation.</param>
    /// <param name="cutoffSd">Anomaly cutoff in units of standard deviation (default 3).</param>
    public static SVType ClassifySV(
        ReadPairSignature pair,
        int expectedInsertSize = 400,
        int insertSizeStdDev = 50,
        double cutoffSd = DefaultInsertSizeCutoffSd)
    {
        // Interchromosomal: linking signature → translocation (Medvedev et al. 2009; BreakDancer CTX).
        if (pair.Chromosome1 != pair.Chromosome2)
            return SVType.Translocation;

        // Same chromosome, same-strand mates: one read's orientation is flipped → inversion.
        if (pair.Strand1 == pair.Strand2)
            return SVType.Inversion;

        // Same chromosome, reverse-forward (RF, outward-facing / everted): the mates have swapped
        // relative order but kept opposite strands → tandem-duplication signature (DELLY, LUMPY,
        // Manta and SVXplorer all read an RF cluster as a duplication candidate). The data model
        // stores the upstream mate first, so RF is strand1 '−', strand2 '+'.
        if (pair.Strand1 == '-' && pair.Strand2 == '+')
            return SVType.Duplication;

        double upperBound = expectedInsertSize + cutoffSd * insertSizeStdDev;
        double lowerBound = expectedInsertSize - cutoffSd * insertSizeStdDev;

        // Forward-reverse (FR) mates from here on: span larger than the insert size → deletion
        // (Medvedev et al. 2009).
        if (pair.InsertSize > upperBound)
            return SVType.Deletion;

        // Span smaller than the insert size → insertion (Medvedev et al. 2009).
        if (pair.InsertSize < lowerBound)
            return SVType.Insertion;

        // Anomalous (e.g. flagged by the maxInsertSize guard) but no basic signature matched.
        return SVType.ComplexRearrangement;
    }

    /// <summary>
    /// Detects structural variants from paired-end mappings: flags discordant pairs, clusters
    /// nearby pairs supporting the same event, and reports an SV for each cluster meeting the
    /// minimum read-pair support. This is the canonical SV-detection entry point.
    /// </summary>
    /// <param name="readPairs">Mapped read pairs (mate coordinates, strands, insert size).</param>
    /// <param name="expectedInsertSize">Mean library insert size.</param>
    /// <param name="insertSizeStdDev">Library insert-size standard deviation.</param>
    /// <param name="cutoffSd">Anomaly cutoff in units of standard deviation (default 3).</param>
    /// <param name="clusterDistance">Maximum coordinate gap to keep adjacent discordant pairs in one cluster.</param>
    /// <param name="minSupport">Minimum supporting read pairs to report an SV (BreakDancer -r, default 2).</param>
    public static IEnumerable<StructuralVariant> DetectSVs(
        IEnumerable<(string ReadId, string Chr1, int Pos1, char Strand1, string Chr2, int Pos2, char Strand2, int InsertSize)> readPairs,
        int expectedInsertSize = 400,
        int insertSizeStdDev = 50,
        double cutoffSd = DefaultInsertSizeCutoffSd,
        int clusterDistance = 500,
        int minSupport = DefaultMinSupport)
    {
        ArgumentNullException.ThrowIfNull(readPairs);

        var discordant = FindDiscordantPairs(
            readPairs, expectedInsertSize, insertSizeStdDev, cutoffSd);

        return ClusterDiscordantPairs(
            discordant, clusterDistance, minSupport,
            expectedInsertSize, insertSizeStdDev, cutoffSd);
    }

    /// <summary>
    /// Clusters discordant read pairs into SV candidates.
    /// </summary>
    public static IEnumerable<StructuralVariant> ClusterDiscordantPairs(
        IEnumerable<ReadPairSignature> discordantPairs,
        int clusterDistance = 500,
        int minSupport = DefaultMinSupport,
        int expectedInsertSize = 400,
        int insertSizeStdDev = 50,
        double cutoffSd = DefaultInsertSizeCutoffSd)
    {
        ArgumentNullException.ThrowIfNull(discordantPairs);

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
            var sv = CreateSVFromCluster(cluster, svId++, expectedInsertSize, insertSizeStdDev, cutoffSd);
            if (sv != null)
            {
                yield return sv.Value;
            }
        }
    }

    private static StructuralVariant? CreateSVFromCluster(
        List<ReadPairSignature> cluster,
        int id,
        int expectedInsertSize,
        int insertSizeStdDev,
        double cutoffSd)
    {
        if (cluster.Count == 0)
            return null;

        var first = cluster[0];
        int start = cluster.Min(p => p.Position1);
        int end = cluster.Max(p => p.Position2);

        // SV type comes from the evidence-based PEM signature of the representative pair
        // (Medvedev, Stanciu & Brudno 2009), not from an ad-hoc size threshold.
        SVType type = ClassifySV(first, expectedInsertSize, insertSizeStdDev, cutoffSd);

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

    /// <summary>
    /// Detects structural-variant breakpoints from split (soft-clipped) reads. Each split read
    /// contributes its aligned/clipped junction coordinate (the single-base breakpoint, i.e. the
    /// "marginal point between a clipped sequence and matched sequence" of Suzuki et al. 2011);
    /// reads on the same chromosome whose junctions agree within <paramref name="clusterTolerance"/>
    /// are grouped, and a breakpoint is reported for every group meeting <paramref name="minSupport"/>.
    /// This is the canonical breakpoint-detection entry point.
    /// </summary>
    /// <remarks>
    /// Method (source-traceable):
    /// <list type="number">
    /// <item>The per-read breakpoint coordinate is <see cref="SplitRead.SupplementaryPosition"/>, the
    /// junction between aligned and clipped portions (ClipCrop, Suzuki et al. 2011; single-base
    /// resolution per Tattini et al. 2015).</item>
    /// <item>Junctions are sorted by chromosome then position and split into clusters whenever the gap
    /// to the previous junction exceeds <paramref name="clusterTolerance"/> or the chromosome changes
    /// (ClipCrop "clustered within 5-base differences"; SAM POS is per-contig).</item>
    /// <item>A cluster is reported only when it contains at least <paramref name="minSupport"/> reads,
    /// and its support count is the cluster size (SoftSearch, Hart et al. 2013: "at least x
    /// soft-clipped reads beginning at position y").</item>
    /// </list>
    /// The reported coordinate is the rounded mean of the member junctions (ASSUMPTION A1: the sources
    /// fix the junction and the tolerance window but not the summary statistic; the mean stays inside
    /// the same ≤ tolerance neighbourhood and does not change cluster membership or support).
    /// </remarks>
    /// <param name="splitReads">Split reads, each carrying its anchored position and junction coordinate.</param>
    /// <param name="clusterTolerance">Maximum junction-position gap (bases) to keep adjacent reads in one breakpoint (ClipCrop, default 5).</param>
    /// <param name="minSupport">Minimum supporting split reads to report a breakpoint (SoftSearch, default 2).</param>
    /// <exception cref="ArgumentNullException"><paramref name="splitReads"/> is null.</exception>
    public static IEnumerable<Breakpoint> FindBreakpoints(
        IEnumerable<SplitRead> splitReads,
        int clusterTolerance = DefaultBreakpointClusterTolerance,
        int minSupport = DefaultBreakpointMinSupport)
    {
        ArgumentNullException.ThrowIfNull(splitReads);

        return FindBreakpointsIterator(splitReads, clusterTolerance, minSupport);
    }

    private static IEnumerable<Breakpoint> FindBreakpointsIterator(
        IEnumerable<SplitRead> splitReads,
        int clusterTolerance,
        int minSupport)
    {
        // Sort by chromosome then by the junction coordinate so a single linear scan can group
        // adjacent junctions (ClipCrop: breakpoints are "sorted and clustered").
        var reads = splitReads
            .OrderBy(r => r.Chromosome, StringComparer.Ordinal)
            .ThenBy(r => r.SupplementaryPosition)
            .ToList();

        if (reads.Count == 0)
            yield break;

        var currentCluster = new List<SplitRead> { reads[0] };

        for (int i = 1; i < reads.Count; i++)
        {
            var prev = reads[i - 1];
            var curr = reads[i];

            // Same breakpoint iff same chromosome and junction within the tolerance window.
            // The gap is computed in 64-bit width: junctions are reference coordinates that
            // may span the full Int32 range, and an Int32 subtraction (or Math.Abs of
            // int.MinValue) would overflow and throw on extreme/opposite-sign coordinates.
            bool sameCluster =
                prev.Chromosome == curr.Chromosome &&
                Math.Abs((long)curr.SupplementaryPosition - prev.SupplementaryPosition) <= clusterTolerance;

            if (sameCluster)
            {
                currentCluster.Add(curr);
            }
            else
            {
                if (currentCluster.Count >= minSupport)
                    yield return CreateBreakpoint(currentCluster);

                currentCluster = new List<SplitRead> { curr };
            }
        }

        if (currentCluster.Count >= minSupport)
            yield return CreateBreakpoint(currentCluster);
    }

    private static Breakpoint CreateBreakpoint(List<SplitRead> cluster)
    {
        // Reported coordinate = rounded mean of member junctions (ASSUMPTION A1; sub-tolerance only).
        int position = (int)Math.Round(cluster.Average(r => (double)r.SupplementaryPosition));

        return new Breakpoint(
            Chromosome1: cluster[0].Chromosome,
            Position1: position,
            Strand1: '+',
            Chromosome2: cluster[0].Chromosome,
            Position2: position,
            Strand2: '-',
            SupportingReads: cluster.Count,
            Quality: Math.Min(cluster.Count * 15.0, 100.0));
    }

    /// <summary>
    /// Refines a candidate breakpoint region to its consensus single-base junction using the split
    /// reads whose junction falls inside the region. Returns the most frequently observed junction
    /// coordinate (the mode; ties broken by the rounded mean of the modal coordinates), or
    /// <see langword="null"/> when no split read supports the region.
    /// </summary>
    /// <remarks>
    /// Source-traceable: the breakpoint is the aligned/clipped junction (Suzuki et al. 2011); the
    /// junction shared by the most reads is the best single-base estimate (SoftSearch, Hart et al.
    /// 2013: a breakpoint is the position where soft-clipped reads accumulate). The region bounds are
    /// inclusive.
    /// </remarks>
    /// <param name="chromosome">Chromosome of the candidate region.</param>
    /// <param name="regionStart">Inclusive start of the candidate region (reference coordinate).</param>
    /// <param name="regionEnd">Inclusive end of the candidate region (reference coordinate).</param>
    /// <param name="splitReads">Split reads to draw junction support from.</param>
    /// <returns>The consensus junction coordinate, or null if no read's junction lies in the region.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="splitReads"/> is null.</exception>
    public static int? RefineBreakpoint(
        string chromosome,
        int regionStart,
        int regionEnd,
        IEnumerable<SplitRead> splitReads)
    {
        ArgumentNullException.ThrowIfNull(splitReads);

        var junctions = splitReads
            .Where(r => r.Chromosome == chromosome &&
                        r.SupplementaryPosition >= regionStart &&
                        r.SupplementaryPosition <= regionEnd)
            .Select(r => r.SupplementaryPosition)
            .ToList();

        if (junctions.Count == 0)
            return null;

        // Consensus = the junction supported by the most reads (mode). The position where clipped
        // reads accumulate is the breakpoint (SoftSearch, Hart et al. 2013).
        int maxCount = junctions
            .GroupBy(p => p)
            .Max(g => g.Count());

        var modalPositions = junctions
            .GroupBy(p => p)
            .Where(g => g.Count() == maxCount)
            .Select(g => g.Key)
            .ToList();

        // Single mode → that coordinate; tie → rounded mean of the tied modal coordinates (they are
        // already within the same support neighbourhood).
        return modalPositions.Count == 1
            ? modalPositions[0]
            : (int)Math.Round(modalPositions.Average(p => (double)p));
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

    #region Read-Depth CNV Detection (SV-CNV-001)

    /// <summary>
    /// Detects copy-number variation from per-position read depth using the read-depth-of-coverage
    /// model. The depth track is summarised into non-overlapping windows of
    /// <paramref name="windowSize"/> positions; each window's mean read depth is converted to a
    /// log2 ratio against a reference depth and then to an integer copy number. This is the canonical
    /// CNV entry point.
    /// </summary>
    /// <remarks>
    /// Source-traceable method:
    /// <list type="number">
    /// <item>Read depth is "counting the number of mapped reads in [fixed-size] windows" and is a
    /// quantitative measure of copy number — there is "a linear relationship between coverage and
    /// copy number" (Yoon et al. 2009, Genome Research 19(9):1586–1592). A trailing partial window
    /// is dropped so every reported window has exactly <paramref name="windowSize"/> positions.</item>
    /// <item>The reference (diploid baseline) depth is, when not supplied, the overall median of the
    /// non-zero window means, mirroring the overall median <c>m</c> of the Yoon GC-correction
    /// equation <c>r' = r·m/m_GC</c> (ASSUMPTION A1).</item>
    /// <item>The log2 ratio is <c>log2(windowMeanRD / referenceRD)</c> and the integer copy number is
    /// <c>round(ploidy · 2^log2)</c> — CNVkit <c>_log2_ratio_to_absolute_pure</c>: <c>n = r·2^v</c>;
    /// for a diploid genome "the absolute copy number is calculated as 2 * 2^(log2 value)".</item>
    /// <item>A window with zero read depth has an undefined ratio (<c>log2 0 = −∞</c>) and is reported
    /// as a no-call (excluded), not a finite call (Yoon RD = read count; CNVkit treats unusable
    /// signal as a no-call). Copy number is clamped to be non-negative (CNVkit <c>max(0, n)</c>).</item>
    /// </list>
    /// </remarks>
    /// <param name="depthData">Per-position read depth (mapped read count) along one contig, 0-based.</param>
    /// <param name="windowSize">Number of positions per non-overlapping window (Yoon et al. used 100).</param>
    /// <param name="referenceDepth">
    /// Reference (copy-number-neutral) mean read depth used as the log2 anchor. When
    /// <see langword="null"/>, the overall median of non-zero window means is used (ASSUMPTION A1).
    /// </param>
    /// <param name="chromosome">Contig label recorded on each emitted segment.</param>
    /// <returns>One <see cref="CopyNumberSegment"/> per full, non-zero-depth window, in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="depthData"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> is not positive.</exception>
    public static IEnumerable<CopyNumberSegment> DetectCNV(
        IReadOnlyList<int> depthData,
        int windowSize = DefaultCnvWindowSize,
        double? referenceDepth = null,
        string chromosome = "chr1")
    {
        ArgumentNullException.ThrowIfNull(depthData);
        ArgumentOutOfRangeException.ThrowIfNegativeOrZero(windowSize);

        return DetectCnvIterator(depthData, windowSize, referenceDepth, chromosome);
    }

    private static IEnumerable<CopyNumberSegment> DetectCnvIterator(
        IReadOnlyList<int> depthData,
        int windowSize,
        double? referenceDepth,
        string chromosome)
    {
        int windowCount = depthData.Count / windowSize;
        if (windowCount == 0)
            yield break;

        // Mean read depth per non-overlapping full window (Yoon et al. 2009: read counts per window).
        var windowMeans = new double[windowCount];
        for (int w = 0; w < windowCount; w++)
        {
            long sum = 0;
            int baseIndex = w * windowSize;
            for (int j = 0; j < windowSize; j++)
                sum += depthData[baseIndex + j];
            windowMeans[w] = (double)sum / windowSize;
        }

        // Reference depth: explicit value, or the overall median of non-zero window means
        // (Yoon overall-median baseline m; ASSUMPTION A1). A zero/absent reference means no call.
        double reference = referenceDepth ?? OverallMedianNonZero(windowMeans);
        if (reference <= 0)
            yield break;

        for (int w = 0; w < windowCount; w++)
        {
            double meanRd = windowMeans[w];

            // Zero-depth window: log2(0) is undefined (−∞) → no-call, not a finite call
            // (Yoon RD = read count; CNVkit no-call for unusable signal).
            if (meanRd <= 0)
                continue;

            double logRatio = Math.Log2(meanRd / reference);
            int copyNumber = LogRatioToCopyNumber(logRatio);

            yield return new CopyNumberSegment(
                Chromosome: chromosome,
                Start: w * windowSize,
                End: w * windowSize + windowSize - 1,
                LogRatio: logRatio,
                CopyNumber: copyNumber,
                BAlleleFrequency: double.NaN, // BAF requires allele-specific data, not provided here.
                ProbeCount: windowSize);
        }
    }

    /// <summary>
    /// Converts a sequence of per-window log2 ratios into copy-number segments, merging maximal runs
    /// of consecutive windows that share the same integer copy number. This is the segmentation
    /// variant of read-depth CNV calling: it consumes log2 ratios directly (e.g. from an external
    /// reference profile) rather than raw depth.
    /// </summary>
    /// <remarks>
    /// Each log2 ratio is converted with the same rule as <see cref="DetectCNV"/>
    /// (<c>CN = round(2 · 2^log2)</c>, CNVkit <c>_log2_ratio_to_absolute_pure</c>); adjacent windows
    /// with equal copy number are merged into one segment whose <see cref="CopyNumberSegment.LogRatio"/>
    /// is the mean of its members and whose <see cref="CopyNumberSegment.ProbeCount"/> is the number of
    /// merged windows. NaN log2 ratios are dropped (no-call). Coordinates are window indices.
    /// </remarks>
    /// <param name="logRatios">Per-window log2 ratios (window i covers coordinate i).</param>
    /// <param name="chromosome">Contig label recorded on each emitted segment.</param>
    /// <returns>Merged copy-number segments in input order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="logRatios"/> is null.</exception>
    public static IEnumerable<CopyNumberSegment> SegmentCopyNumber(
        IEnumerable<double> logRatios,
        string chromosome = "chr1")
    {
        ArgumentNullException.ThrowIfNull(logRatios);

        return SegmentCopyNumberIterator(logRatios, chromosome);
    }

    private static IEnumerable<CopyNumberSegment> SegmentCopyNumberIterator(
        IEnumerable<double> logRatios,
        string chromosome)
    {
        int index = 0;
        int? runCopyNumber = null;
        int runStart = 0;
        var runLogRatios = new List<double>();

        foreach (double logRatio in logRatios)
        {
            int position = index++;

            // NaN ratio → no usable signal → no-call window; it breaks the current run.
            if (double.IsNaN(logRatio))
            {
                if (runCopyNumber.HasValue)
                {
                    yield return BuildSegment(chromosome, runStart, position - 1, runLogRatios, runCopyNumber.Value);
                    runCopyNumber = null;
                    runLogRatios = new List<double>();
                }
                continue;
            }

            int copyNumber = LogRatioToCopyNumber(logRatio);

            if (runCopyNumber == copyNumber)
            {
                runLogRatios.Add(logRatio);
            }
            else
            {
                if (runCopyNumber.HasValue)
                    yield return BuildSegment(chromosome, runStart, position - 1, runLogRatios, runCopyNumber.Value);

                runCopyNumber = copyNumber;
                runStart = position;
                runLogRatios = new List<double> { logRatio };
            }
        }

        if (runCopyNumber.HasValue)
            yield return BuildSegment(chromosome, runStart, index - 1, runLogRatios, runCopyNumber.Value);
    }

    private static CopyNumberSegment BuildSegment(
        string chromosome, int start, int end, List<double> logRatios, int copyNumber) =>
        new(
            Chromosome: chromosome,
            Start: start,
            End: end,
            LogRatio: logRatios.Average(),
            CopyNumber: copyNumber,
            BAlleleFrequency: double.NaN,
            ProbeCount: logRatios.Count);

    /// <summary>
    /// Converts a log2 ratio to an integer copy number for a diploid genome:
    /// <c>round(ploidy · 2^log2)</c>, clamped to be non-negative.
    /// </summary>
    /// <remarks>
    /// CNVkit <c>_log2_ratio_to_absolute_pure</c>: <c>n = r · 2^v</c> (here r = ploidy = 2), rounded to
    /// the nearest integer; copy number is physically ≥ 0 (CNVkit <c>max(0, n)</c>). CNVkit's
    /// <c>do_call</c> rounds with NumPy's <c>ndarray.round()</c>, which uses round-half-to-even
    /// (banker's rounding: numpy.round(0.5)=0, numpy.round(2.5)=2), so the conversion here uses
    /// <see cref="MidpointRounding.ToEven"/> to match the reference implementation exactly at
    /// half-integer copy numbers.
    /// </remarks>
    private static int LogRatioToCopyNumber(double logRatio)
    {
        double copies = DiploidPloidy * Math.Pow(2, logRatio);
        int rounded = (int)Math.Round(copies, MidpointRounding.ToEven);
        return Math.Max(0, rounded);
    }

    /// <summary>
    /// Returns the median of the strictly positive (non-zero) values, or 0 when none are positive.
    /// </summary>
    /// <remarks>
    /// The reference baseline is the overall median <c>m</c> of windows (Yoon et al. 2009 GC-correction
    /// equation); zero-depth windows are no-calls and excluded from the baseline.
    /// </remarks>
    private static double OverallMedianNonZero(IReadOnlyList<double> values)
    {
        var positive = values.Where(v => v > 0).OrderBy(v => v).ToList();
        if (positive.Count == 0)
            return 0;

        int mid = positive.Count / 2;
        return positive.Count % 2 == 1
            ? positive[mid]
            : (positive[mid - 1] + positive[mid]) / 2.0;
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
