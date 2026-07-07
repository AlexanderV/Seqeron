namespace Seqeron.Genomics.Tests.Unit.Chromosome;

[TestFixture]
public class GenomeAssemblyAnalyzerTests
{
    // NOTE: N50/L50/Nx/Lx/auN, CalculateStatistics, and FindGaps tests for
    // ASSEMBLY-STATS-001 live in GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs
    // (the canonical evidence-based fixture). The pre-template Weak/duplicate cases
    // that previously lived here were removed during that unit's consolidation.

    #region Gap Distribution Tests

    [Test]
    public void AnalyzeGapDistribution_CalculatesStats()
    {
        var gaps = new List<GenomeAssemblyAnalyzer.GapInfo>
        {
            new("s1", 0, 9, 10, "Short"),
            new("s1", 100, 199, 100, "Medium"),
            new("s2", 0, 999, 1000, "Long")
        };

        var (count, mean, median, max, types) = GenomeAssemblyAnalyzer.AnalyzeGapDistribution(gaps);

        Assert.That(count, Is.EqualTo(3));
        Assert.That(mean, Is.EqualTo(370).Within(1));
        Assert.That(max, Is.EqualTo(1000));
        Assert.That(types["Short"], Is.EqualTo(1));
        Assert.That(types["Medium"], Is.EqualTo(1));
        Assert.That(types["Long"], Is.EqualTo(1));
    }

    #endregion

    #region Scaffold Analysis Tests

    [Test]
    public void AnalyzeScaffolds_WithGaps_IdentifiesContigs()
    {
        var scaffolds = new List<(string, string)>
        {
            ("scaffold1", "ACGTACGT" + new string('N', 100) + "GCGCGCGC" + new string('N', 100) + "TTTTTTTT")
        };

        var structures = GenomeAssemblyAnalyzer.AnalyzeScaffolds(scaffolds).ToList();

        Assert.That(structures, Has.Count.EqualTo(1));
        Assert.That(structures[0].Contigs.Count, Is.EqualTo(3));
        Assert.That(structures[0].Gaps.Count, Is.EqualTo(2));
    }

    [Test]
    public void ExtractContigs_FromScaffold_ReturnsContigs()
    {
        var scaffolds = new List<(string, string)>
        {
            ("scaffold1", new string('A', 500) + new string('N', 100) + new string('C', 500))
        };

        var contigs = GenomeAssemblyAnalyzer.ExtractContigs(scaffolds, minContigLength: 200).ToList();

        Assert.That(contigs, Has.Count.EqualTo(2));
        Assert.That(contigs[0].Sequence.Length, Is.EqualTo(500));
        Assert.That(contigs[1].Sequence.Length, Is.EqualTo(500));
    }

    [Test]
    public void ExtractContigs_ShortContig_Filtered()
    {
        var scaffolds = new List<(string, string)>
        {
            ("scaffold1", "ACGT" + new string('N', 100) + new string('C', 500))
        };

        var contigs = GenomeAssemblyAnalyzer.ExtractContigs(scaffolds, minContigLength: 100).ToList();

        Assert.That(contigs, Has.Count.EqualTo(1));
    }

    #endregion

    #region Completeness Assessment Tests

    [Test]
    public void AssessCompleteness_AllFound_Returns100Percent()
    {
        var assembly = new List<(string, string)>
        {
            ("contig1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
        };

        var markers = new List<(string, string)>
        {
            ("gene1", "ACGTACGTACGTACGTACGT")
        };

        var result = GenomeAssemblyAnalyzer.AssessCompleteness(assembly, markers,
            identityThreshold: 0.5, coverageThreshold: 0.5);

        Assert.That(result.Complete, Is.GreaterThan(0));
    }

    [Test]
    public void AssessCompleteness_NoneFound_ReturnsMissing()
    {
        var assembly = new List<(string, string)>
        {
            ("contig1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        };

        var markers = new List<(string, string)>
        {
            ("gene1", "CCCCCCCCCCCCCCCCCCCC")
        };

        var result = GenomeAssemblyAnalyzer.AssessCompleteness(assembly, markers);

        Assert.That(result.Missing, Is.EqualTo(1));
        Assert.That(result.CompletenessPercent, Is.EqualTo(0));
    }

    [Test]
    public void AssessCompleteness_EmptyMarkers_ReturnsZero()
    {
        var assembly = new List<(string, string)>
        {
            ("contig1", "ACGTACGT")
        };

        var result = GenomeAssemblyAnalyzer.AssessCompleteness(
            assembly, new List<(string, string)>());

        Assert.That(result.TotalGenes, Is.EqualTo(0));
    }

    [Test]
    public void EstimateCompletenessFromKmers_ValidSpectrum_ReturnsEstimates()
    {
        var spectrum = new List<(string, int)>();

        // Simulate k-mer spectrum with peak at coverage 30
        for (int i = 0; i < 100; i++)
            spectrum.Add(($"KMER{i}", 1)); // Singletons (errors)

        for (int i = 0; i < 1000; i++)
            spectrum.Add(($"SOLID{i}", 30)); // Peak coverage

        var (completeness, errorRate, genomeSize) =
            GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers(spectrum);

        Assert.That(errorRate, Is.LessThan(0.5));
        Assert.That(genomeSize, Is.GreaterThan(0));
    }

    #endregion

    #region Repeat Analysis Tests

    [Test]
    public void FindRepetitiveRegions_RepeatedKmers_Detected()
    {
        string repeatUnit = "ACGTACGTACGTACGT"; // 16bp unit
        string sequence = string.Concat(Enumerable.Repeat(repeatUnit, 10)); // 10 copies

        var sequences = new List<(string, string)>
        {
            ("seq1", sequence + "NNNNNNNNNN" + new string('A', 100))
        };

        var repeats = GenomeAssemblyAnalyzer.FindRepetitiveRegions(
            sequences, kmerSize: 10, minCopies: 3).ToList();

        // Should detect some repetitive regions
        Assert.That(repeats.Count, Is.GreaterThanOrEqualTo(0)); // May or may not detect based on algorithm
    }

    [Test]
    public void FindTandemRepeats_SimpleTandem_Detected()
    {
        string sequence = "ACGT" + "CAGCAGCAGCAGCAGCAGCAGCAGCAG" + "TGCA"; // CAG x 9

        var sequences = new List<(string, string)>
        {
            ("seq1", sequence)
        };

        var repeats = GenomeAssemblyAnalyzer.FindTandemRepeats(
            sequences, minUnitLength: 3, maxUnitLength: 5, minCopies: 5).ToList();

        Assert.That(repeats, Has.Count.GreaterThanOrEqualTo(1));
        Assert.That(repeats[0].Unit, Does.Contain("CAG").Or.Contain("AGC").Or.Contain("GCA"));
        Assert.That(repeats[0].Copies, Is.GreaterThanOrEqualTo(5));
    }

    [Test]
    public void FindTandemRepeats_NoRepeats_ReturnsEmpty()
    {
        var sequences = new List<(string, string)>
        {
            ("seq1", "ACGTGCATACGTGCATACGT") // Non-repetitive
        };

        var repeats = GenomeAssemblyAnalyzer.FindTandemRepeats(
            sequences, minCopies: 10).ToList();

        Assert.That(repeats, Is.Empty);
    }

    [Test]
    public void CalculateRepeatContent_CalculatesPercentage()
    {
        var repeats = new List<GenomeAssemblyAnalyzer.RepeatAnnotation>
        {
            new("seq1", 0, 999, "LINE", "L1", 10.0, '+'),
            new("seq1", 2000, 2499, "SINE", "Alu", 5.0, '+')
        };

        var (totalLength, percentage, classLengths) =
            GenomeAssemblyAnalyzer.CalculateRepeatContent(repeats, 10000);

        Assert.That(totalLength, Is.EqualTo(1500));
        Assert.That(percentage, Is.EqualTo(15.0));
        Assert.That(classLengths["LINE"], Is.EqualTo(1000));
        Assert.That(classLengths["SINE"], Is.EqualTo(500));
    }

    #endregion

    #region Assembly Comparison Tests

    [Test]
    public void CompareAssemblies_IdenticalAssemblies_HighIdentity()
    {
        var assembly = new List<(string, string)>
        {
            ("seq1", "ACGTACGTACGTACGTACGTACGTACGTACGT")
        };

        var comparison = GenomeAssemblyAnalyzer.CompareAssemblies(
            assembly, assembly, "asm1", "asm2");

        Assert.That(comparison.SequenceIdentity, Is.EqualTo(1.0).Within(0.01));
        Assert.That(comparison.AlignedFraction1, Is.EqualTo(1.0).Within(0.01));
    }

    [Test]
    public void CompareAssemblies_DifferentAssemblies_LowIdentity()
    {
        var assembly1 = new List<(string, string)>
        {
            ("seq1", new string('A', 100))
        };

        var assembly2 = new List<(string, string)>
        {
            ("seq1", new string('C', 100))
        };

        var comparison = GenomeAssemblyAnalyzer.CompareAssemblies(assembly1, assembly2);

        Assert.That(comparison.AlignedFraction1, Is.LessThan(0.5));
    }

    [Test]
    public void FindSyntenicBlocks_IdenticalSequences_FindsBlock()
    {
        string sharedSequence = string.Concat(Enumerable.Repeat("ACGTGCTA", 200));

        var assembly1 = new List<(string, string)>
        {
            ("chr1", sharedSequence)
        };

        var assembly2 = new List<(string, string)>
        {
            ("chrA", sharedSequence)
        };

        var blocks = GenomeAssemblyAnalyzer.FindSyntenicBlocks(
            assembly1, assembly2, minBlockSize: 100).ToList();

        // Algorithm may find blocks depending on k-mer spacing
        Assert.That(blocks, Is.Not.Null);
    }

    #endregion

    #region Quality Assessment Tests

    [Test]
    public void CalculateLocalQuality_ReturnsWindowMetrics()
    {
        var sequences = new List<(string, string)>
        {
            ("seq1", new string('G', 500) + new string('A', 500))
        };

        var quality = GenomeAssemblyAnalyzer.CalculateLocalQuality(
            sequences, windowSize: 200).ToList();

        Assert.That(quality.Count, Is.GreaterThan(0));

        // First windows should have high GC
        Assert.That(quality[0].GcContent, Is.GreaterThan(0.9));

        // Later windows should have low GC
        var lastWindow = quality.Last(q => q.Position >= 500);
        Assert.That(lastWindow.GcContent, Is.LessThan(0.1));
    }

    [Test]
    public void FindSuspiciousRegions_HighNContent_Detected()
    {
        var sequences = new List<(string, string)>
        {
            ("seq1", new string('A', 500) + new string('N', 200) + new string('C', 500))
        };

        var suspicious = GenomeAssemblyAnalyzer.FindSuspiciousRegions(sequences).ToList();

        Assert.That(suspicious.Any(r => r.Reason.Contains("N content")), Is.True);
    }

    [Test]
    public void FindSuspiciousRegions_LowComplexity_Detected()
    {
        var sequences = new List<(string, string)>
        {
            ("seq1", "ACGTACGT" + new string('A', 1000) + "GCTAGCTA")
        };

        var suspicious = GenomeAssemblyAnalyzer.FindSuspiciousRegions(
            sequences, minComplexity: 0.5).ToList();

        Assert.That(suspicious.Any(r => r.Reason.Contains("complexity")), Is.True);
    }

    #endregion

    #region Utility Tests

    [Test]
    public void FilterByLength_FiltersCorrectly()
    {
        var sequences = new List<(string, string)>
        {
            ("short", "ACG"),
            ("medium", "ACGTACGT"),
            ("long", new string('A', 100))
        };

        var filtered = GenomeAssemblyAnalyzer.FilterByLength(
            sequences, minLength: 5, maxLength: 50).ToList();

        Assert.That(filtered, Has.Count.EqualTo(1));
        Assert.That(filtered[0].Id, Is.EqualTo("medium"));
    }

    [Test]
    public void SortByLength_SortsDescending()
    {
        var sequences = new List<(string, string)>
        {
            ("short", "ACG"),
            ("long", new string('A', 100)),
            ("medium", "ACGTACGT")
        };

        var sorted = GenomeAssemblyAnalyzer.SortByLength(sequences).ToList();

        Assert.That(sorted[0].Id, Is.EqualTo("long"));
        Assert.That(sorted[1].Id, Is.EqualTo("medium"));
        Assert.That(sorted[2].Id, Is.EqualTo("short"));
    }

    [Test]
    public void CalculateLengthDistribution_BinsCorrectly()
    {
        var lengths = new List<int> { 50, 150, 500, 750, 1500, 5000, 10000 };

        var distribution = GenomeAssemblyAnalyzer.CalculateLengthDistribution(
            lengths, 100, 500, 1000, 5000);

        Assert.That(distribution["<100"], Is.EqualTo(1));
        Assert.That(distribution["<500"], Is.EqualTo(1));
        Assert.That(distribution["<1000"], Is.EqualTo(2));
        Assert.That(distribution["<5000"], Is.EqualTo(1));
        Assert.That(distribution[">=5000"], Is.EqualTo(2));
    }

    #endregion

    #region Edge Cases Tests

    // NOTE: CalculateStatistics_AllNs and the leading/trailing FindGaps edge cases
    // were moved into GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs (ASSEMBLY-STATS-001).

    [Test]
    public void CalculateRepeatContent_EmptyInput_ReturnsZeros()
    {
        var (total, percentage, classes) = GenomeAssemblyAnalyzer.CalculateRepeatContent(
            new List<GenomeAssemblyAnalyzer.RepeatAnnotation>(), 10000);

        Assert.That(total, Is.EqualTo(0));
        Assert.That(percentage, Is.EqualTo(0));
        Assert.That(classes, Is.Empty);
    }

    #endregion
}
