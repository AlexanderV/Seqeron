namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for chromosome analysis.
/// Verifies telomere detection, centromere analysis, and karyotype.
///
/// Test Units: CHROM-TELO-001, CHROM-CENT-001, CHROM-KARYO-001, CHROM-SYNT-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Chromosome")]
public class ChromosomeSnapshotTests
{
    [Test]
    public Task AnalyzeTelomeres_KnownSequence_MatchesSnapshot()
    {
        // Build a sequence with telomere repeats at both ends
        string telomere = string.Concat(Enumerable.Repeat("TTAGGG", 200));
        string body = new('N', 5000);
        string chromosome = telomere + body + string.Concat(Enumerable.Repeat("CCCTAA", 200));

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chr1", chromosome,
            searchLength: 2000, minTelomereLength: 100, criticalLength: 500);

        return Verify(new
        {
            result.Chromosome,
            result.Has5PrimeTelomere,
            result.TelomereLength5Prime,
            result.Has3PrimeTelomere,
            result.TelomereLength3Prime,
            RepeatPurity5 = Math.Round(result.RepeatPurity5Prime, 4),
            RepeatPurity3 = Math.Round(result.RepeatPurity3Prime, 4),
            result.IsCriticallyShort
        });
    }

    [Test]
    public Task AnalyzeKaryotype_KnownChromosomes_MatchesSnapshot()
    {
        var chromosomes = new (string Name, long Length, bool IsSexChromosome)[]
        {
            ("chr1", 248956422, false),
            ("chr2", 242193529, false),
            ("chr3", 198295559, false),
            ("chrX", 156040895, true),
            ("chrY", 57227415, true)
        };

        var karyotype = ChromosomeAnalyzer.AnalyzeKaryotype(chromosomes, expectedPloidyLevel: 2);
        return Verify(new
        {
            karyotype.TotalChromosomes,
            karyotype.AutosomeCount,
            SexChromosomes = karyotype.SexChromosomes,
            karyotype.PloidyLevel,
            karyotype.TotalGenomeSize,
            karyotype.MeanChromosomeLength,
            karyotype.HasAneuploidy
        });
    }

    [Test]
    public Task FindSyntenyBlocks_KnownOrthologPairs_MatchesSnapshot()
    {
        var orthologPairs = new[]
        {
            ("chr1", 100, 200, "GeneA", "chr2", 300, 400, "GeneA'"),
            ("chr1", 250, 350, "GeneB", "chr2", 450, 550, "GeneB'"),
            ("chr1", 400, 500, "GeneC", "chr2", 600, 700, "GeneC'"),
            ("chr1", 550, 650, "GeneD", "chr2", 750, 850, "GeneD'"),
            ("chr3", 100, 200, "GeneE", "chr4", 100, 200, "GeneE'"),
            ("chr3", 250, 350, "GeneF", "chr4", 250, 350, "GeneF'"),
            ("chr3", 400, 500, "GeneG", "chr4", 400, 500, "GeneG'")
        }.Select(t => (Chr1: t.Item1, Start1: t.Item2, End1: t.Item3, Gene1: t.Item4,
                       Chr2: t.Item5, Start2: t.Item6, End2: t.Item7, Gene2: t.Item8));

        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(orthologPairs, minGenes: 3)
            .Select(b => new { b.Species1Chromosome, b.Species2Chromosome, b.GeneCount })
            .ToList();

        return Verify(new { BlockCount = blocks.Count, Blocks = blocks });
    }

    [Test]
    public Task DetectAneuploidy_KnownDepthData_MatchesSnapshot()
    {
        // Simulate depth data: chr1 normal (2N), chr2 trisomy (3N)
        var depthData = new List<(string Chromosome, int Position, double Depth)>();
        for (int i = 0; i < 10; i++)
        {
            depthData.Add(("chr1", i * 1000000, 30.0 + (i % 3))); // ~normal
            depthData.Add(("chr2", i * 1000000, 45.0 + (i % 3))); // ~1.5x = trisomy
        }

        var states = ChromosomeAnalyzer.DetectAneuploidy(depthData, medianDepth: 30.0, binSize: 2000000)
            .Select(s => new { s.Chromosome, s.CopyNumber, s.LogRatio, s.Confidence })
            .ToList();

        return Verify(new { States = states });
    }

    [Test]
    public Task AnalyzeCentromere_KnownSequence_MatchesSnapshot()
    {
        // Build a sequence with alpha-satellite repeats in center
        string flank = new string('A', 5000);
        string alphaRepeat = string.Concat(Enumerable.Repeat("AATGGAAATGAAATCCAACT", 300)); // ~171bp alpha satellite unit approximation
        string chromosome = flank + alphaRepeat + flank;

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chr1", chromosome,
            windowSize: 1000, minAlphaSatelliteContent: 0.1);

        return Verify(new
        {
            result.Chromosome,
            result.Start,
            result.End,
            result.Length,
            result.CentromereType,
            AlphaSatellite = Math.Round(result.AlphaSatelliteContent, 4),
            result.IsAcrocentric
        });
    }
}
