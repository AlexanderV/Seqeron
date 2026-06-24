using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Chromosome;

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

    /// <summary>
    /// Length in base pairs of the human alpha-satellite (alphoid) monomer — the fundamental
    /// tandemly-repeated unit of human centromeric DNA.
    /// Source: Willard HF (1985); Waye JS, Willard HF (1987); review Hartley G, O'Neill RJ (2019),
    /// Alpha satellite DNA biology (PMC6121732): "Alpha satellite DNA is composed of fundamental
    /// 171bp monomeric repeat units."
    /// </summary>
    public const int AlphaSatelliteMonomerLength = 171;

    /// <summary>
    /// Canonical 17-bp CENP-B box consensus motif, written with IUPAC ambiguity codes:
    /// 5'-YTTCGTTGGAARCGGGA-3' (Y = C/T, R = A/G).
    /// Source: Masumoto H, Masukata H, Muro Y, Nozaki N, Okazaki T (1989), J Cell Biol 109(4):1963-1973;
    /// consensus as reported in PMC6121732 ("5'-T/CTCGTTGGAAA/GCGGGA-3'") and PMC4843215
    /// ("YTTCGTTGGAARCGGGA"). CENP-B binds this 17-bp sequence within a subset of alpha-satellite monomers.
    /// </summary>
    public const string CenpBBoxConsensus = "YTTCGTTGGAARCGGGA";

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
    /// Alpha-satellite (alphoid) specific detection result.
    /// Unlike <see cref="CentromereResult.AlphaSatelliteContent"/> (a generic tandem-repeat-density
    /// score), this captures alpha-satellite-specific molecular signatures: a ~171 bp tandem
    /// periodicity and CENP-B box occurrences.
    /// </summary>
    /// <param name="IsAlphaSatellite">True when both molecular signatures are met:
    /// strong ~171 bp tandem periodicity AND the sequence is AT-rich (alpha satellite is AT-rich).</param>
    /// <param name="PeriodicityScore">Self-similarity (fraction of bases identical to the base
    /// <see cref="AlphaSatelliteMonomerLength"/> positions upstream), in [0,1]. ~1 for a clean
    /// tandem array, ~0.25 for random DNA.</param>
    /// <param name="BestPeriod">The repeat period (bp) with the highest self-similarity within the
    /// searched tolerance window around 171 bp, or 0 if none searched.</param>
    /// <param name="AtContent">Fraction of A/T bases in [0,1].</param>
    /// <param name="CenpBBoxCount">Number of CENP-B box (17-bp) matches found.</param>
    public readonly record struct AlphaSatelliteResult(
        bool IsAlphaSatellite,
        double PeriodicityScore,
        int BestPeriod,
        double AtContent,
        int CenpBBoxCount);

    /// <summary>
    /// Higher-order repeat (HOR) structure of an alpha-satellite array.
    /// An alpha-satellite HOR is a block of <see cref="MonomersPerUnit"/> distinct ~171 bp monomers
    /// that is itself tandemly repeated with high identity between copies (inter-HOR), while the
    /// monomers WITHIN one unit are much more divergent (intra-HOR). Source: McNulty SM, Sullivan BA
    /// (2018), "Alpha satellite DNA biology" (PMC6121732): "A defined number of individual monomers …
    /// that are 50–70% identical in sequence are arranged tandemly to form a HOR unit"; "HOR within a
    /// given array are 97–100% identical".
    /// </summary>
    /// <param name="HasHigherOrderStructure">True when a HOR period ≥ 2 monomers was detected (the
    /// array is organised into multi-monomer HOR units), false for a purely monomeric / single-monomer
    /// homogeneous array.</param>
    /// <param name="MonomersPerUnit">The HOR period: number of monomers per HOR unit (the smallest
    /// block size k such that monomers k apart are inter-HOR identical across the array). 1 when there
    /// is no multi-monomer HOR organisation.</param>
    /// <param name="HorUnitLengthBp">Length of one HOR unit in base pairs
    /// (<see cref="MonomersPerUnit"/> × 171 bp).</param>
    /// <param name="HorCopyNumber">Number of complete HOR units tiling the analysed monomers
    /// (⌊monomer count / <see cref="MonomersPerUnit"/>⌋).</param>
    /// <param name="MonomerCount">Number of ~171 bp monomers the array was split into.</param>
    /// <param name="MeanInterHorIdentity">Mean percent identity between monomers at the same position
    /// in different HOR copies (i and i+period). NaN when there are no inter-HOR pairs.</param>
    /// <param name="MeanIntraHorIdentity">Mean percent identity between distinct monomers within one
    /// HOR unit. NaN when the unit has a single monomer.</param>
    public readonly record struct HorResult(
        bool HasHigherOrderStructure,
        int MonomersPerUnit,
        int HorUnitLengthBp,
        int HorCopyNumber,
        int MonomerCount,
        double MeanInterHorIdentity,
        double MeanIntraHorIdentity);

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
                abnormalities.Add($"{GetAneuploidyTerm(count)} {group.Key}");
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
    /// Returns standard cytogenetic aneuploidy term for a given copy count.
    /// Per ISCN / standard nomenclature (Wikipedia: Aneuploidy).
    /// </summary>
    private static string GetAneuploidyTerm(int copyCount) => copyCount switch
    {
        0 => "Nullisomy",
        1 => "Monosomy",
        2 => "Disomy",
        3 => "Trisomy",
        4 => "Tetrasomy",
        5 => "Pentasomy",
        _ => $"Polysomy ({copyCount} copies)"
    };

    /// <summary>
    /// Gets base chromosome name (strips copy suffixes).
    /// </summary>
    private static string GetChromosomeBaseName(string name)
    {
        // Strips a trailing "_N" integer copy suffix, e.g. "chr1_2" -> "chr1".
        // Only numeric "_N" suffixes are removed; letter suffixes (e.g. "chr1a") are left intact.
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

        var sorted = depths.OrderBy(d => d).ToList();
        double medianDepth = sorted.Count % 2 == 1
            ? sorted[sorted.Count / 2]
            : (sorted[sorted.Count / 2 - 1] + sorted[sorted.Count / 2]) / 2.0;
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
    /// Determines centromere type based on arm ratio per Levan et al. (1964).
    /// Source: Levan A, Fredga K, Sandberg AA. "Nomenclature for centromeric position
    /// on chromosomes". Hereditas. 1964;52(2):201-220.
    /// </summary>
    private static string DetermineCentromereType(int chromosomeLength, int? centStart, int? centEnd)
    {
        if (!centStart.HasValue || !centEnd.HasValue)
            return "Unknown";

        int centMid = (centStart.Value + centEnd.Value) / 2;
        int pArm = Math.Min(centMid, chromosomeLength - centMid);
        int qArm = Math.Max(centMid, chromosomeLength - centMid);

        if (pArm == 0)
            return "Telocentric";

        double armRatio = (double)qArm / pArm;

        return armRatio switch
        {
            <= 1.7 => "Metacentric",
            <= 3.0 => "Submetacentric",
            < 7.0 => "Subtelocentric",
            _ => "Acrocentric"
        };
    }

    #endregion

    #region Alpha-Satellite-Specific Detection

    // Sourced parameters for alpha-satellite-specific detection.
    //
    // Monomer period: alpha satellite is a tandem repeat of a ~171 bp monomer
    //   (Willard 1985; Waye & Willard 1987; review PMC6121732). See AlphaSatelliteMonomerLength.
    // Period tolerance: monomers diverge and indels occur, so the observed tandem period varies
    //   slightly around 171 bp; we scan a small window around the canonical length and take the best.
    // AT-richness: alpha satellite is described as an "AT-rich 171-bp alphoid monomer" (PMC6121732).
    //   A balanced AT content is 0.5; alpha satellite sits above it, so we require AT > 0.5.
    // Periodicity threshold: monomers within an array share 50-70% identity (PMC6121732), so
    //   base-level self-similarity at the monomer period for a genuine tandem array is well above
    //   the ~0.25 expected for random DNA. We require periodicity >= 0.50 (the lower bound of the
    //   reported 50-70% monomer identity).

    /// <summary>Half-width (bp) of the period search window scanned around the 171 bp monomer length.</summary>
    private const int MonomerPeriodTolerance = 5;

    /// <summary>Minimum base-level self-similarity at the monomer period to call a tandem array
    /// (lower bound of the 50-70% intra-array monomer identity, PMC6121732).</summary>
    private const double MinPeriodicityScore = 0.50;

    /// <summary>Minimum AT fraction for the AT-rich alpha-satellite signature (above the 0.5 balance point).</summary>
    private const double MinAlphaSatelliteAtContent = 0.50;

    /// <summary>Length (bp) of the CENP-B box motif (Masumoto et al. 1989).</summary>
    private const int CenpBBoxLength = 17;

    /// <summary>
    /// Detects alpha-satellite (alphoid)-specific signal in a sequence: a ~171 bp tandem periodicity
    /// combined with AT-richness, plus a count of CENP-B box occurrences. This is alpha-satellite
    /// SPECIFIC, in contrast to <see cref="AnalyzeCentromere"/> whose
    /// <see cref="CentromereResult.AlphaSatelliteContent"/> is a generic tandem-repeat-density score.
    /// </summary>
    /// <param name="sequence">DNA sequence to test (case-insensitive).</param>
    /// <returns>An <see cref="AlphaSatelliteResult"/>. For an empty/too-short sequence, returns a
    /// no-signal result (not alpha-satellite, zero scores).</returns>
    public static AlphaSatelliteResult DetectAlphaSatellite(string sequence)
    {
        // Need at least two monomers plus the tolerance to measure periodicity at the monomer period.
        int minLength = AlphaSatelliteMonomerLength + MonomerPeriodTolerance + 1;
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minLength)
            return new AlphaSatelliteResult(false, 0, 0, 0, 0);

        sequence = sequence.ToUpperInvariant();

        // 1) AT content (alpha satellite is AT-rich).
        int atCount = 0;
        int acgtCount = 0;
        foreach (char c in sequence)
        {
            if (c is 'A' or 'T') { atCount++; acgtCount++; }
            else if (c is 'C' or 'G') { acgtCount++; }
        }
        double atContent = acgtCount > 0 ? atCount / (double)acgtCount : 0;

        // 2) Tandem periodicity: scan periods around the 171 bp monomer length and keep the best
        //    base-level self-similarity (fraction of positions identical to the base `period` upstream).
        double bestScore = 0;
        int bestPeriod = 0;
        int lowPeriod = AlphaSatelliteMonomerLength - MonomerPeriodTolerance;
        int highPeriod = AlphaSatelliteMonomerLength + MonomerPeriodTolerance;

        for (int period = lowPeriod; period <= highPeriod; period++)
        {
            if (period >= sequence.Length)
                break;

            int matches = 0;
            int comparisons = sequence.Length - period;
            for (int i = period; i < sequence.Length; i++)
            {
                if (sequence[i] == sequence[i - period])
                    matches++;
            }

            double score = comparisons > 0 ? matches / (double)comparisons : 0;
            if (score > bestScore)
            {
                bestScore = score;
                bestPeriod = period;
            }
        }

        // 3) CENP-B box occurrences.
        int cenpBCount = CountCenpBBoxes(sequence);

        bool isAlphaSatellite = bestScore >= MinPeriodicityScore && atContent > MinAlphaSatelliteAtContent;

        return new AlphaSatelliteResult(isAlphaSatellite, bestScore, bestPeriod, atContent, cenpBCount);
    }

    /// <summary>
    /// Finds the start positions (0-based) of all forward-strand CENP-B box matches in a sequence.
    /// The CENP-B box is the 17-bp consensus 5'-YTTCGTTGGAARCGGGA-3' (Masumoto et al. 1989), where
    /// Y = C/T and R = A/G; all other consensus positions must match exactly.
    /// </summary>
    /// <param name="sequence">DNA sequence to scan (case-insensitive).</param>
    /// <returns>0-based start indices of each match, in ascending order.</returns>
    public static IReadOnlyList<int> FindCenpBBoxes(string sequence)
    {
        var positions = new List<int>();
        if (string.IsNullOrEmpty(sequence) || sequence.Length < CenpBBoxLength)
            return positions;

        string upper = sequence.ToUpperInvariant();
        for (int i = 0; i <= upper.Length - CenpBBoxLength; i++)
        {
            if (MatchesIupac(upper, i, CenpBBoxConsensus))
                positions.Add(i);
        }

        return positions;
    }

    // --- Higher-order repeat (HOR) structure detection (opt-in; additive to monomer detection) ---
    //
    // Sourced thresholds. An alpha-satellite HOR is a block of N ~171 bp monomers tandemly repeated;
    // copies of the same HOR position are near-identical (inter-HOR) while distinct monomers within a
    // unit are divergent (intra-HOR):
    //   - Inter-HOR identity: "HOR within a given array are 97–100% identical"; HOR copies "differ in
    //     sequence by only a few percent"; "mutual sequence divergence of <5%". (McNulty & Sullivan 2018,
    //     PMC6121732; Rosandić et al. 2024, PMC11050224; Alkan et al. 2007.) We require >= 95% identity
    //     between monomers k apart to accept k as a HOR period (the conservative <5%-divergence bound).
    //   - Intra-HOR monomer identity: monomers within a unit "are 50–70% identical" / "differ in
    //     sequence by 10–40%" (PMC6121732). These are NOT inter-HOR identical, so a true HOR period is
    //     the SMALLEST k for which the k-periodic identity clears the inter-HOR bar — at k=1 (adjacent
    //     monomers) a multi-monomer HOR fails the bar because adjacent monomers are only 50–70% identical.
    //   - HOR period definition: "HOR unit length is determined by where the next monomer shows nearly
    //     total sequence identity to the first monomer in the HOR" (PMC6121732) — i.e. monomer-level
    //     identity periodicity, exactly the k we search for.

    /// <summary>Minimum percent identity between two monomers k apart for k to be accepted as a HOR
    /// period — the conservative inter-HOR bound (HOR copies differ by &lt;5%, i.e. are 97–100%
    /// identical; PMC6121732 / Rosandić 2024 / Alkan 2007).</summary>
    private const double InterHorMinIdentityPercent = 95.0;

    /// <summary>Fraction of k-periodic monomer pairs that must clear
    /// <see cref="InterHorMinIdentityPercent"/> for k to count as a consistent HOR period. A HOR array
    /// is highly homogeneous, so the periodicity must hold across (essentially all of) the array, not
    /// at a single pair.</summary>
    private const double HorPeriodConsistencyFraction = 0.90;

    /// <summary>
    /// Detects higher-order repeat (HOR) structure in an alpha-satellite array. The array is split
    /// into consecutive ~171 bp monomers; the monomer-vs-monomer identity is computed with the
    /// library aligner (<see cref="SequenceAligner.GlobalAlign(string,string,ScoringMatrix?)"/> +
    /// <see cref="SequenceAligner.CalculateStatistics"/>); the HOR period is the SMALLEST block size
    /// k (≥ 1) such that monomers k apart are inter-HOR identical (≥ 95%) consistently across the
    /// array. This is opt-in and additive: it does not change <see cref="DetectAlphaSatellite"/>,
    /// <see cref="AnalyzeCentromere"/>, or the Levan classification.
    /// </summary>
    /// <param name="sequence">Alpha-satellite array (case-insensitive). Should already be alpha-satellite;
    /// this method does not re-test the alphoid signature.</param>
    /// <param name="monomerLength">Monomer length used to split the array (default 171 bp, the alphoid
    /// monomer; Willard 1985 / PMC6121732).</param>
    /// <returns>An <see cref="HorResult"/>. For fewer than two monomers, returns a no-structure result
    /// (period 1, copy number = monomer count, NaN identities).</returns>
    public static HorResult DetectHigherOrderRepeat(string sequence, int monomerLength = AlphaSatelliteMonomerLength)
    {
        if (monomerLength < 1)
            throw new ArgumentOutOfRangeException(nameof(monomerLength), "Monomer length must be positive.");

        if (string.IsNullOrEmpty(sequence))
            return new HorResult(false, 1, monomerLength, 0, 0, double.NaN, double.NaN);

        string upper = sequence.ToUpperInvariant();

        // 1) Split the array into consecutive full monomers (trailing partial monomer is ignored).
        int monomerCount = upper.Length / monomerLength;
        if (monomerCount < 2)
            return new HorResult(false, 1, monomerLength, monomerCount, monomerCount, double.NaN, double.NaN);

        var monomers = new string[monomerCount];
        for (int i = 0; i < monomerCount; i++)
            monomers[i] = upper.Substring(i * monomerLength, monomerLength);

        // 2) Pairwise monomer identity (percent), via the library global aligner. Cached on demand so
        //    each (i,j) pair is aligned at most once.
        var identityCache = new Dictionary<(int, int), double>();
        double Identity(int a, int b)
        {
            if (a == b) return 100.0;
            var key = a < b ? (a, b) : (b, a);
            if (identityCache.TryGetValue(key, out double cached))
                return cached;
            var alignment = SequenceAligner.GlobalAlign(monomers[key.Item1], monomers[key.Item2]);
            double id = SequenceAligner.CalculateStatistics(alignment).Identity;
            identityCache[key] = id;
            return id;
        }

        // 3) Find the smallest period k (1 <= k <= monomerCount/2) whose k-periodic monomer pairs are
        //    inter-HOR identical consistently across the array. k = 1 means a homogeneous single-monomer
        //    repeat (no multi-monomer HOR); k >= 2 is a genuine HOR period.
        int detectedPeriod = 0;
        int maxPeriod = monomerCount / 2;
        for (int k = 1; k <= maxPeriod; k++)
        {
            int pairs = 0, consistent = 0;
            for (int i = 0; i + k < monomerCount; i++)
            {
                pairs++;
                if (Identity(i, i + k) >= InterHorMinIdentityPercent)
                    consistent++;
            }

            if (pairs > 0 && consistent >= HorPeriodConsistencyFraction * pairs)
            {
                detectedPeriod = k;
                break;
            }
        }

        // 4) No periodic high-identity structure at any k: the monomers are mutually divergent and not
        //    HOR-organised. Report no structure (period 1 = treat each monomer as its own unit).
        if (detectedPeriod == 0)
        {
            double intraNoHor = MeanPairwiseIdentity(monomers.Length, 0, monomers.Length, Identity);
            return new HorResult(false, 1, monomerLength, monomerCount, monomerCount, double.NaN, intraNoHor);
        }

        int copyNumber = monomerCount / detectedPeriod;
        int unitLengthBp = detectedPeriod * monomerLength;

        // 5) Inter-HOR identity: mean identity between monomers at the same unit position in different
        //    copies (i and i+period across the array).
        double interSum = 0; int interN = 0;
        for (int i = 0; i + detectedPeriod < monomerCount; i++)
        {
            interSum += Identity(i, i + detectedPeriod);
            interN++;
        }
        double meanInter = interN > 0 ? interSum / interN : double.NaN;

        // 6) Intra-HOR identity: mean identity between the distinct monomers WITHIN the first unit.
        double meanIntra = detectedPeriod >= 2
            ? MeanPairwiseIdentity(detectedPeriod, 0, detectedPeriod, Identity)
            : double.NaN;

        bool hasHor = detectedPeriod >= 2;
        return new HorResult(hasHor, detectedPeriod, unitLengthBp, copyNumber, monomerCount, meanInter, meanIntra);
    }

    /// <summary>
    /// Mean of all distinct unordered pairwise identities among monomer indices [start, end).
    /// Returns NaN when fewer than two indices are in range.
    /// </summary>
    private static double MeanPairwiseIdentity(int count, int start, int end, Func<int, int, double> identity)
    {
        double sum = 0; int n = 0;
        for (int i = start; i < end && i < count; i++)
            for (int j = i + 1; j < end && j < count; j++)
            {
                sum += identity(i, j);
                n++;
            }

        return n > 0 ? sum / n : double.NaN;
    }

    /// <summary>Counts CENP-B box matches (see <see cref="FindCenpBBoxes"/>).</summary>
    private static int CountCenpBBoxes(string upperSequence)
    {
        int count = 0;
        for (int i = 0; i <= upperSequence.Length - CenpBBoxLength; i++)
        {
            if (MatchesIupac(upperSequence, i, CenpBBoxConsensus))
                count++;
        }

        return count;
    }

    /// <summary>
    /// Returns true if the window of <paramref name="upperSequence"/> starting at
    /// <paramref name="start"/> matches the IUPAC <paramref name="pattern"/>. Only the ambiguity
    /// codes occurring in the CENP-B box consensus (Y = C/T, R = A/G) are supported.
    /// </summary>
    private static bool MatchesIupac(string upperSequence, int start, string pattern)
    {
        for (int j = 0; j < pattern.Length; j++)
        {
            char c = upperSequence[start + j];
            char p = pattern[j];
            bool ok = p switch
            {
                'Y' => c is 'C' or 'T',
                'R' => c is 'A' or 'G',
                _ => c == p
            };
            if (!ok)
                return false;
        }

        return true;
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
                            double.NaN); // Not computable from coordinate-only input
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

                // Check for deletion (same species2 chromosome, same strand, asymmetric gap)
                // Per Wikipedia (Chromosomal rearrangement): deletion = segment is removed
                if (current.Species2Chromosome == next.Species2Chromosome &&
                    current.Strand == next.Strand)
                {
                    int gap1 = next.Species1Start - current.Species1End;
                    int gap2 = current.Strand == '+'
                        ? next.Species2Start - current.Species2End
                        : current.Species2Start - next.Species2End;

                    if (gap1 > 0 && gap2 >= 0 && gap1 > gap2 * 2)
                    {
                        yield return new ChromosomalRearrangement(
                            "Deletion",
                            current.Species1Chromosome,
                            current.Species1End,
                            null,
                            next.Species1Start,
                            gap1 - Math.Max(0, gap2),
                            $"Deletion in species 2: {current.Species1Chromosome}:{current.Species1End}-{next.Species1Start}");
                    }
                }
            }
        }

        // Detect duplications: overlapping species 1 regions mapping to different species 2 locations
        // Per Wikipedia (Chromosomal rearrangement): duplication = segment is copied
        for (int i = 0; i < blocks.Count; i++)
        {
            for (int j = i + 1; j < blocks.Count; j++)
            {
                if (blocks[i].Species1Chromosome == blocks[j].Species1Chromosome)
                {
                    bool overlaps = blocks[i].Species1Start < blocks[j].Species1End &&
                                   blocks[j].Species1Start < blocks[i].Species1End;

                    if (overlaps)
                    {
                        bool sameTarget = blocks[i].Species2Chromosome == blocks[j].Species2Chromosome &&
                            blocks[i].Species2Start == blocks[j].Species2Start &&
                            blocks[i].Species2End == blocks[j].Species2End;

                        if (!sameTarget)
                        {
                            int overlapStart = Math.Max(blocks[i].Species1Start, blocks[j].Species1Start);
                            int overlapEnd = Math.Min(blocks[i].Species1End, blocks[j].Species1End);

                            yield return new ChromosomalRearrangement(
                                "Duplication",
                                blocks[i].Species1Chromosome,
                                overlapStart,
                                blocks[j].Species2Chromosome,
                                blocks[j].Species2Start,
                                overlapEnd - overlapStart,
                                $"Duplication: {blocks[i].Species1Chromosome}:{overlapStart}-{overlapEnd} maps to multiple locations");
                        }
                    }
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
                        5 => "Pentasomy",
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
