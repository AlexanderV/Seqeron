using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics;

/// <summary>
/// Provides chromosome-level analysis algorithms.
/// Includes karyotyping, centromere/telomere detection, synteny analysis, and aneuploidy detection.
/// </summary>
public static class ChromosomeAnalyzer
{
    #region Constants

    /// <summary>
    /// Human telomere repeat sequence.
    /// </summary>
    public const string HumanTelomereRepeat = "TTAGGG";

    /// <summary>
    /// Human centromere alpha-satellite consensus.
    /// </summary>
    public const string AlphaSatelliteConsensus = "AATGAATATTTCTTTTATGTTCCTTAAAGTAGAAATGTCAAGAATATGTTAAGCCTTAAATG";

    #endregion

    #region Records

    /// <summary>
    /// Represents a chromosome.
    /// </summary>
    public readonly record struct Chromosome(
        string Name,
        long Length,
        int? CentromereStart,
        int? CentromereEnd,
        int? TelomereStartLength,
        int? TelomereEndLength,
        double GcContent,
        string? CytogeneticBand);

    /// <summary>
    /// Karyotype information.
    /// </summary>
    public readonly record struct Karyotype(
        int TotalChromosomes,
        int AutosomeCount,
        IReadOnlyList<string> SexChromosomes,
        long TotalGenomeSize,
        double MeanChromosomeLength,
        int PloidyLevel,
        bool HasAneuploidy,
        IReadOnlyList<string> Abnormalities);

    /// <summary>
    /// Cytogenetic band (G-band).
    /// </summary>
    public readonly record struct CytogeneticBand(
        string Chromosome,
        int Start,
        int End,
        string Name,
        string Stain,
        double GcContent,
        double GeneDensity);

    /// <summary>
    /// Telomere analysis result.
    /// </summary>
    public readonly record struct TelomereResult(
        string Chromosome,
        bool Has5PrimeTelomere,
        int TelomereLength5Prime,
        bool Has3PrimeTelomere,
        int TelomereLength3Prime,
        double RepeatPurity5Prime,
        double RepeatPurity3Prime,
        bool IsCriticallyShort);

    /// <summary>
    /// Centromere analysis result.
    /// </summary>
    public readonly record struct CentromereResult(
        string Chromosome,
        int? Start,
        int? End,
        int Length,
        string CentromereType,
        double AlphaSatelliteContent,
        bool IsAcrocentric);

    /// <summary>
    /// Synteny block between species.
    /// </summary>
    public readonly record struct SyntenyBlock(
        string Species1Chromosome,
        int Species1Start,
        int Species1End,
        string Species2Chromosome,
        int Species2Start,
        int Species2End,
        char Strand,
        int GeneCount,
        double SequenceIdentity);

    /// <summary>
    /// Chromosomal rearrangement.
    /// </summary>
    public readonly record struct ChromosomalRearrangement(
        string Type,
        string Chromosome1,
        int Position1,
        string? Chromosome2,
        int? Position2,
        int? Size,
        string? Description);

    /// <summary>
    /// Copy number state.
    /// </summary>
    public readonly record struct CopyNumberState(
        string Chromosome,
        int Start,
        int End,
        int CopyNumber,
        double LogRatio,
        double Confidence);

    #endregion

    #region Karyotype Analysis

    /// <summary>
    /// Analyzes karyotype from chromosome data.
    /// </summary>
    public static Karyotype AnalyzeKaryotype(
        IEnumerable<(string Name, long Length, bool IsSexChromosome)> chromosomes,
        int expectedPloidyLevel = 2)
    {
        var chromList = chromosomes.ToList();

        if (chromList.Count == 0)
        {
            return new Karyotype(0, 0, new List<string>(), 0, 0, 0, false, new List<string>());
        }

        var sexChroms = chromList.Where(c => c.IsSexChromosome).Select(c => c.Name).ToList();
        var autosomes = chromList.Where(c => !c.IsSexChromosome).ToList();

        long totalSize = chromList.Sum(c => c.Length);
        double meanLength = totalSize / (double)chromList.Count;

        // Detect aneuploidy (simplified - looks for missing or extra chromosomes)
        var abnormalities = new List<string>();
        bool hasAneuploidy = false;

        // Group autosomes by base name (e.g., "chr1" from "chr1_1", "chr1_2")
        var autosomeGroups = autosomes
            .GroupBy(c => GetChromosomeBaseName(c.Name))
            .ToList();

        foreach (var group in autosomeGroups)
        {
            int count = group.Count();
            if (count != expectedPloidyLevel)
            {
                hasAneuploidy = true;
                if (count < expectedPloidyLevel)
                    abnormalities.Add($"Monosomy {group.Key}");
                else if (count > expectedPloidyLevel)
                    abnormalities.Add($"Trisomy {group.Key}");
            }
        }

        return new Karyotype(
            chromList.Count,
            autosomes.Count,
            sexChroms,
            totalSize,
            meanLength,
            expectedPloidyLevel,
            hasAneuploidy,
            abnormalities);
    }

    /// <summary>
    /// Gets base chromosome name (strips copy suffixes).
    /// </summary>
    private static string GetChromosomeBaseName(string name)
    {
        // Handle formats like "chr1_1", "chr1_2" or "chr1a", "chr1b"
        int underscoreIdx = name.LastIndexOf('_');
        if (underscoreIdx > 0 && underscoreIdx < name.Length - 1)
        {
            if (int.TryParse(name[(underscoreIdx + 1)..], out _))
                return name[..underscoreIdx];
        }

        return name;
    }

    /// <summary>
    /// Detects ploidy level from read depth.
    /// </summary>
    public static (int PloidyLevel, double Confidence) DetectPloidy(
        IEnumerable<double> normalizedDepths,
        double expectedDiploidDepth = 1.0)
    {
        var depths = normalizedDepths.ToList();

        if (depths.Count == 0)
            return (2, 0);

        double medianDepth = depths.OrderBy(d => d).ElementAt(depths.Count / 2);
        double ratio = medianDepth / expectedDiploidDepth;

        // Determine ploidy
        int ploidy = (int)Math.Round(ratio * 2);
        ploidy = Math.Max(1, Math.Min(8, ploidy)); // Limit to reasonable range

        // Calculate confidence based on how close to integer ploidy
        double fractionalPart = Math.Abs(ratio * 2 - ploidy);
        double confidence = 1.0 - fractionalPart * 2;

        return (ploidy, Math.Max(0, confidence));
    }

    #endregion

    #region Telomere Analysis

    /// <summary>
    /// Analyzes telomeres at chromosome ends.
    /// </summary>
    public static TelomereResult AnalyzeTelomeres(
        string chromosomeName,
        string sequence,
        string telomereRepeat = "TTAGGG",
        int searchLength = 10000,
        int minTelomereLength = 500,
        int criticalLength = 3000)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new TelomereResult(chromosomeName, false, 0, false, 0, 0, 0, true);
        }

        sequence = sequence.ToUpperInvariant();
        telomereRepeat = telomereRepeat.ToUpperInvariant();
        string telomereRepeatRC = DnaSequence.GetReverseComplementString(telomereRepeat);

        // Analyze 5' end (should have CCCTAA repeats = reverse complement)
        int search5End = Math.Min(searchLength, sequence.Length);
        var (length5, purity5) = MeasureTelomereLength(
            sequence[..search5End], telomereRepeatRC, fromEnd: false);

        // Analyze 3' end (should have TTAGGG repeats)
        int search3Start = Math.Max(0, sequence.Length - searchLength);
        var (length3, purity3) = MeasureTelomereLength(
            sequence[search3Start..], telomereRepeat, fromEnd: true);

        bool has5Prime = length5 >= minTelomereLength;
        bool has3Prime = length3 >= minTelomereLength;
        bool isCritical = (has5Prime && length5 < criticalLength) ||
                          (has3Prime && length3 < criticalLength);

        return new TelomereResult(
            chromosomeName,
            has5Prime, length5,
            has3Prime, length3,
            purity5, purity3,
            isCritical);
    }

    /// <summary>
    /// Measures telomere length and repeat purity.
    /// </summary>
    private static (int Length, double Purity) MeasureTelomereLength(
        string region,
        string repeatUnit,
        bool fromEnd)
    {
        int repeatLen = repeatUnit.Length;
        if (region.Length < repeatLen)
            return (0, 0);

        int telomereLength = 0;
        int matchingBases = 0;
        int totalBases = 0;

        int start = fromEnd ? region.Length - repeatLen : 0;
        int step = fromEnd ? -repeatLen : repeatLen;

        while (true)
        {
            if (start < 0 || start + repeatLen > region.Length)
                break;

            string window = region.Substring(start, repeatLen);
            int matches = 0;

            for (int i = 0; i < repeatLen; i++)
            {
                if (window[i] == repeatUnit[i])
                    matches++;
            }

            double similarity = matches / (double)repeatLen;

            if (similarity >= 0.7) // Allow some divergence
            {
                telomereLength += repeatLen;
                matchingBases += matches;
                totalBases += repeatLen;
                start += step;
            }
            else
            {
                break;
            }
        }

        double purity = totalBases > 0 ? matchingBases / (double)totalBases : 0;
        return (telomereLength, purity);
    }

    /// <summary>
    /// Estimates telomere length from qPCR T/S ratio.
    /// </summary>
    public static double EstimateTelomereLengthFromTSRatio(
        double tsRatio,
        double referenceRatio = 1.0,
        double referenceLength = 7000)
    {
        // T/S ratio is proportional to telomere length
        return referenceLength * tsRatio / referenceRatio;
    }

    #endregion

    #region Centromere Analysis

    /// <summary>
    /// Analyzes centromere region.
    /// </summary>
    public static CentromereResult AnalyzeCentromere(
        string chromosomeName,
        string sequence,
        int windowSize = 100000,
        double minAlphaSatelliteContent = 0.3)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new CentromereResult(chromosomeName, null, null, 0, "Unknown", 0, false);
        }

        sequence = sequence.ToUpperInvariant();

        // Scan for regions with high repetitive content and low GC variability
        int? centStart = null;
        int? centEnd = null;
        double maxScore = 0;

        for (int i = 0; i < sequence.Length - windowSize; i += windowSize / 4)
        {
            int end = Math.Min(i + windowSize, sequence.Length);
            string window = sequence[i..end];

            // Check for alpha-satellite-like content
            double repeatContent = EstimateRepeatContent(window);
            double gcVariability = CalculateGcVariability(window, 1000);

            // Centromeres have high repeat content and low GC variability
            double score = repeatContent * (1 - gcVariability);

            if (score > maxScore && repeatContent > minAlphaSatelliteContent)
            {
                maxScore = score;
                centStart = i;
                centEnd = end;
            }
        }

        // Extend centromere boundaries
        if (centStart.HasValue && centEnd.HasValue)
        {
            // Extend left
            while (centStart > windowSize / 2)
            {
                string window = sequence[(centStart.Value - windowSize / 2)..centStart.Value];
                if (EstimateRepeatContent(window) >= minAlphaSatelliteContent * 0.7)
                    centStart -= windowSize / 2;
                else
                    break;
            }

            // Extend right
            while (centEnd < sequence.Length - windowSize / 2)
            {
                string window = sequence[centEnd.Value..(centEnd.Value + windowSize / 2)];
                if (EstimateRepeatContent(window) >= minAlphaSatelliteContent * 0.7)
                    centEnd += windowSize / 2;
                else
                    break;
            }
        }

        int length = centStart.HasValue && centEnd.HasValue ? centEnd.Value - centStart.Value : 0;

        // Determine centromere type
        string centType = DetermineCentromereType(sequence.Length, centStart, centEnd);
        bool isAcrocentric = centType == "Acrocentric";

        return new CentromereResult(
            chromosomeName,
            centStart,
            centEnd,
            length,
            centType,
            maxScore,
            isAcrocentric);
    }

    /// <summary>
    /// Estimates repeat content using k-mer frequency.
    /// </summary>
    private static double EstimateRepeatContent(string sequence, int kmerSize = 15)
    {
        if (sequence.Length < kmerSize * 2)
            return 0;

        var kmerCounts = new Dictionary<string, int>();

        for (int i = 0; i <= sequence.Length - kmerSize; i++)
        {
            string kmer = sequence.Substring(i, kmerSize);
            if (!kmer.Contains('N'))
            {
                kmerCounts[kmer] = kmerCounts.GetValueOrDefault(kmer) + 1;
            }
        }

        if (kmerCounts.Count == 0)
            return 0;

        // Count k-mers appearing more than once
        int repeatedKmers = kmerCounts.Values.Count(c => c > 1);
        int totalRepeatInstances = kmerCounts.Values.Where(c => c > 1).Sum();

        return totalRepeatInstances / (double)(sequence.Length - kmerSize + 1);
    }

    /// <summary>
    /// Calculates GC content variability.
    /// </summary>
    private static double CalculateGcVariability(string sequence, int windowSize)
    {
        var gcValues = new List<double>();

        for (int i = 0; i < sequence.Length - windowSize; i += windowSize)
        {
            string window = sequence.Substring(i, windowSize);
            gcValues.Add(window.CalculateGcFractionFast());
        }

        if (gcValues.Count < 2)
            return 0;

        double mean = gcValues.Average();
        double variance = gcValues.Sum(v => (v - mean) * (v - mean)) / gcValues.Count;

        return Math.Sqrt(variance);
    }

    /// <summary>
    /// Determines centromere type based on position.
    /// </summary>
    private static string DetermineCentromereType(int chromosomeLength, int? centStart, int? centEnd)
    {
        if (!centStart.HasValue || !centEnd.HasValue)
            return "Unknown";

        int centMid = (centStart.Value + centEnd.Value) / 2;
        double position = centMid / (double)chromosomeLength;

        return position switch
        {
            < 0.15 => "Acrocentric",
            < 0.35 => "Submetacentric",
            < 0.65 => "Metacentric",
            < 0.85 => "Submetacentric",
            _ => "Acrocentric"
        };
    }

    #endregion

    #region Cytogenetic Bands

    /// <summary>
    /// Predicts G-band pattern from sequence.
    /// </summary>
    public static IEnumerable<CytogeneticBand> PredictGBands(
        string chromosomeName,
        string sequence,
        int bandSize = 5000000,
        double darkBandGcThreshold = 0.37,
        double lightBandGcThreshold = 0.45)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        sequence = sequence.ToUpperInvariant();
        int bandNumber = 1;
        int arm = 1; // p arm = 1, q arm = 2

        for (int i = 0; i < sequence.Length; i += bandSize)
        {
            int end = Math.Min(i + bandSize, sequence.Length);
            string region = sequence[i..end];

            // Calculate GC content
            int total = region.Count(c => c != 'N');
            double gcContent = total > 0 ? region.CalculateGcFractionFast() : 0.5;

            // Determine stain type
            string stain;
            if (gcContent < darkBandGcThreshold)
                stain = "gpos100"; // Dark band
            else if (gcContent < lightBandGcThreshold)
                stain = "gpos50";  // Medium band
            else
                stain = "gneg";    // Light band

            // Estimate gene density (simplified - AT-rich regions have lower gene density)
            double geneDensity = gcContent * 2; // Simplified correlation

            string bandName = $"{chromosomeName}{(arm == 1 ? "p" : "q")}{bandNumber}";

            yield return new CytogeneticBand(
                chromosomeName,
                i,
                end - 1,
                bandName,
                stain,
                gcContent,
                geneDensity);

            bandNumber++;

            // Switch arms at midpoint (simplified)
            if (i + bandSize >= sequence.Length / 2 && arm == 1)
            {
                arm = 2;
                bandNumber = 1;
            }
        }
    }

    /// <summary>
    /// Identifies heterochromatin regions.
    /// </summary>
    public static IEnumerable<(int Start, int End, string Type)> FindHeterochromatinRegions(
        string sequence,
        int windowSize = 100000,
        double minRepeatContent = 0.5)
    {
        if (string.IsNullOrEmpty(sequence))
            yield break;

        sequence = sequence.ToUpperInvariant();
        int? regionStart = null;

        for (int i = 0; i < sequence.Length - windowSize; i += windowSize / 2)
        {
            string window = sequence.Substring(i, windowSize);
            double repeatContent = EstimateRepeatContent(window);

            if (repeatContent >= minRepeatContent)
            {
                if (!regionStart.HasValue)
                    regionStart = i;
            }
            else if (regionStart.HasValue)
            {
                string type = DetermineHeterochromatinType(sequence, regionStart.Value, i);
                yield return (regionStart.Value, i - 1, type);
                regionStart = null;
            }
        }

        if (regionStart.HasValue)
        {
            yield return (regionStart.Value, sequence.Length - 1, "Constitutive");
        }
    }

    /// <summary>
    /// Determines heterochromatin type.
    /// </summary>
    private static string DetermineHeterochromatinType(string sequence, int start, int end)
    {
        // Check position
        double position = (start + end) / 2.0 / sequence.Length;

        if (position < 0.05 || position > 0.95)
            return "Telomeric";
        if (position > 0.45 && position < 0.55)
            return "Centromeric";

        return "Constitutive";
    }

    #endregion

    #region Synteny Analysis

    /// <summary>
    /// Identifies synteny blocks between two genomes.
    /// </summary>
    public static IEnumerable<SyntenyBlock> FindSyntenyBlocks(
        IEnumerable<(string Chr1, int Start1, int End1, string Gene1,
                    string Chr2, int Start2, int End2, string Gene2)> orthologPairs,
        int minGenes = 3,
        int maxGap = 10)
    {
        var pairs = orthologPairs.ToList();

        if (pairs.Count < minGenes)
            yield break;

        // Group by chromosome pairs
        var chromPairs = pairs.GroupBy(p => (p.Chr1, p.Chr2));

        foreach (var group in chromPairs)
        {
            // Sort by position in first genome
            var sorted = group.OrderBy(p => p.Start1).ToList();

            // Find collinear runs
            int blockStart = 0;
            bool isForward = true;

            for (int i = 1; i < sorted.Count; i++)
            {
                var prev = sorted[i - 1];
                var curr = sorted[i];

                // Check if positions are collinear
                bool collinear = false;
                bool currentForward = curr.Start2 > prev.End2;

                if (i == 1)
                    isForward = currentForward;

                if (currentForward == isForward)
                {
                    int gap1 = curr.Start1 - prev.End1;
                    int gap2 = Math.Abs(curr.Start2 - prev.End2);

                    if (gap1 <= maxGap * 1000000 && gap2 <= maxGap * 1000000)
                        collinear = true;
                }

                if (!collinear || i == sorted.Count - 1)
                {
                    int blockEnd = collinear ? i : i - 1;
                    int geneCount = blockEnd - blockStart + 1;

                    if (geneCount >= minGenes)
                    {
                        var blockGenes = sorted.Skip(blockStart).Take(geneCount).ToList();
                        var first = blockGenes.First();
                        var last = blockGenes.Last();

                        yield return new SyntenyBlock(
                            group.Key.Chr1,
                            first.Start1,
                            last.End1,
                            group.Key.Chr2,
                            Math.Min(first.Start2, last.Start2),
                            Math.Max(first.End2, last.End2),
                            isForward ? '+' : '-',
                            geneCount,
                            0.9); // Placeholder identity
                    }

                    blockStart = i;
                    if (i < sorted.Count - 1)
                        isForward = sorted[i + 1].Start2 > curr.End2;
                }
            }
        }
    }

    /// <summary>
    /// Detects chromosomal rearrangements from synteny blocks.
    /// </summary>
    public static IEnumerable<ChromosomalRearrangement> DetectRearrangements(
        IEnumerable<SyntenyBlock> syntenyBlocks)
    {
        var blocks = syntenyBlocks.OrderBy(b => b.Species1Chromosome)
                                  .ThenBy(b => b.Species1Start)
                                  .ToList();

        for (int i = 0; i < blocks.Count - 1; i++)
        {
            var current = blocks[i];
            var next = blocks[i + 1];

            // Same chromosome in species1
            if (current.Species1Chromosome == next.Species1Chromosome)
            {
                // Check for inversion
                if (current.Species2Chromosome == next.Species2Chromosome &&
                    current.Strand != next.Strand)
                {
                    yield return new ChromosomalRearrangement(
                        "Inversion",
                        current.Species1Chromosome,
                        current.Species1End,
                        null,
                        next.Species1Start,
                        next.Species1Start - current.Species1End,
                        $"Inversion between {current.Species1Chromosome}:{current.Species1End}-{next.Species1Start}");
                }

                // Check for translocation
                if (current.Species2Chromosome != next.Species2Chromosome)
                {
                    yield return new ChromosomalRearrangement(
                        "Translocation",
                        current.Species1Chromosome,
                        current.Species1End,
                        next.Species2Chromosome,
                        next.Species2Start,
                        null,
                        $"Translocation from {current.Species2Chromosome} to {next.Species2Chromosome}");
                }
            }
        }
    }

    #endregion

    #region Aneuploidy Detection

    /// <summary>
    /// Detects aneuploidy from read depth data.
    /// </summary>
    public static IEnumerable<CopyNumberState> DetectAneuploidy(
        IEnumerable<(string Chromosome, int Position, double Depth)> depthData,
        double medianDepth,
        int binSize = 1000000)
    {
        var data = depthData.ToList();

        if (data.Count == 0 || medianDepth <= 0)
            yield break;

        // Group by chromosome
        var byChrom = data.GroupBy(d => d.Chromosome);

        foreach (var chromGroup in byChrom)
        {
            // Bin the data
            var bins = chromGroup
                .GroupBy(d => d.Position / binSize)
                .OrderBy(g => g.Key)
                .ToList();

            foreach (var bin in bins)
            {
                double meanDepth = bin.Average(d => d.Depth);
                double logRatio = Math.Log2(meanDepth / medianDepth);

                // Determine copy number
                int copyNumber = (int)Math.Round(Math.Pow(2, logRatio) * 2);
                copyNumber = Math.Max(0, Math.Min(10, copyNumber));

                // Calculate confidence
                double expected = copyNumber / 2.0;
                double observed = Math.Pow(2, logRatio);
                double confidence = 1.0 - Math.Min(1.0, Math.Abs(expected - observed));

                yield return new CopyNumberState(
                    chromGroup.Key,
                    bin.Key * binSize,
                    (bin.Key + 1) * binSize - 1,
                    copyNumber,
                    logRatio,
                    confidence);
            }
        }
    }

    /// <summary>
    /// Identifies whole chromosome aneuploidy.
    /// </summary>
    public static IEnumerable<(string Chromosome, int CopyNumber, string Type)> IdentifyWholeChromosomeAneuploidy(
        IEnumerable<CopyNumberState> copyNumberStates,
        double minFraction = 0.8)
    {
        var states = copyNumberStates.ToList();

        var byChrom = states.GroupBy(s => s.Chromosome);

        foreach (var chromGroup in byChrom)
        {
            var cnCounts = chromGroup
                .GroupBy(s => s.CopyNumber)
                .Select(g => (CopyNumber: g.Key, Fraction: g.Count() / (double)chromGroup.Count()))
                .OrderByDescending(g => g.Fraction)
                .ToList();

            if (cnCounts.Count > 0)
            {
                var dominant = cnCounts.First();

                if (dominant.Fraction >= minFraction && dominant.CopyNumber != 2)
                {
                    string type = dominant.CopyNumber switch
                    {
                        0 => "Nullisomy",
                        1 => "Monosomy",
                        3 => "Trisomy",
                        4 => "Tetrasomy",
                        _ => $"Copy number = {dominant.CopyNumber}"
                    };

                    yield return (chromGroup.Key, dominant.CopyNumber, type);
                }
            }
        }
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Calculates chromosome arm ratio (p/q).
    /// </summary>
    public static double CalculateArmRatio(int centromerePosition, int chromosomeLength)
    {
        if (centromerePosition <= 0 || chromosomeLength <= 0)
            return 0;

        int pArmLength = centromerePosition;
        int qArmLength = chromosomeLength - centromerePosition;

        return qArmLength > 0 ? pArmLength / (double)qArmLength : 0;
    }

    /// <summary>
    /// Classifies chromosome by arm ratio.
    /// </summary>
    public static string ClassifyChromosomeByArmRatio(double armRatio)
    {
        return armRatio switch
        {
            >= 0.9 and <= 1.1 => "Metacentric",
            >= 0.5 and < 0.9 => "Submetacentric",
            >= 0.2 and < 0.5 => "Acrocentric",
            < 0.2 => "Telocentric",
            > 1.1 and <= 2.0 => "Submetacentric",
            > 2.0 and <= 5.0 => "Acrocentric",
            > 5.0 => "Telocentric",
            _ => "Unknown"
        };
    }

    /// <summary>
    /// Estimates chromosome age from telomere length.
    /// </summary>
    public static double EstimateCellDivisionsFromTelomereLength(
        int currentLength,
        int birthLength = 15000,
        int lossPerDivision = 50)
    {
        if (lossPerDivision <= 0)
            return 0;

        int lost = birthLength - currentLength;
        return Math.Max(0, lost / (double)lossPerDivision);
    }

    #endregion
}
