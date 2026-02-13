namespace Seqeron.Genomics.Tests.Snapshots;

/// <summary>
/// Snapshot tests for metagenomics analysis.
/// Verifies alpha/beta diversity and contig binning output stability.
///
/// Test Units: META-ALPHA-001, META-BETA-001, META-BIN-001
/// </summary>
[TestFixture]
[Category("Snapshot")]
[Category("Metagenomics")]
public class MetagenomicsSnapshotTests
{
    [Test]
    public Task AlphaDiversity_KnownAbundances_MatchesSnapshot()
    {
        var abundances = new Dictionary<string, double>
        {
            ["E_coli"] = 0.35,
            ["B_subtilis"] = 0.25,
            ["S_aureus"] = 0.20,
            ["P_aeruginosa"] = 0.12,
            ["L_monocytogenes"] = 0.08
        };
        var alpha = MetagenomicsAnalyzer.CalculateAlphaDiversity(abundances);
        return Verify(new
        {
            Shannon = Math.Round(alpha.ShannonIndex, 4),
            Simpson = Math.Round(alpha.SimpsonIndex, 4),
            InvSimpson = Math.Round(alpha.InverseSimpson, 4),
            Chao1 = Math.Round(alpha.Chao1Estimate, 4),
            alpha.ObservedSpecies,
            Evenness = Math.Round(alpha.PielouEvenness, 4)
        });
    }

    [Test]
    public Task BetaDiversity_TwoSamples_MatchesSnapshot()
    {
        var sample1 = new Dictionary<string, double>
        {
            ["E_coli"] = 0.40,
            ["B_subtilis"] = 0.30,
            ["S_aureus"] = 0.20,
            ["P_aeruginosa"] = 0.10
        };
        var sample2 = new Dictionary<string, double>
        {
            ["E_coli"] = 0.10,
            ["B_subtilis"] = 0.15,
            ["S_aureus"] = 0.50,
            ["L_monocytogenes"] = 0.25
        };
        var beta = MetagenomicsAnalyzer.CalculateBetaDiversity("Soil", sample1, "Water", sample2);
        return Verify(new
        {
            beta.Sample1,
            beta.Sample2,
            BrayCurtis = Math.Round(beta.BrayCurtis, 4),
            Jaccard = Math.Round(beta.JaccardDistance, 4),
            beta.SharedSpecies,
            beta.UniqueToSample1,
            beta.UniqueToSample2
        });
    }

    [Test]
    public Task ClassifyReads_KnownDatabase_MatchesSnapshot()
    {
        // Build a minimal k-mer database
        var kmerDb = new Dictionary<string, string>
        {
            ["ACGTACGTACGTACGTACGTACGTACGTACGT"] = "E_coli",
            ["TTTTAAAACCCCGGGGTTTTAAAACCCCGGGG"] = "B_subtilis"
        };
        var reads = new[]
        {
            ("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            ("read2", "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"),
            ("read3", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
        };
        var classifications = MetagenomicsAnalyzer.ClassifyReads(reads, kmerDb, k: 31)
            .Select(c => new { c.ReadId, c.Confidence, c.MatchedKmers })
            .ToList();

        return Verify(new { Classifications = classifications });
    }

    [Test]
    public Task GenerateTaxonomicProfile_MatchesSnapshot()
    {
        // Create mock classifications using positional construction
        var classifications = new MetagenomicsAnalyzer.TaxonomicClassification[]
        {
            new("read1", "Bacteria", "Proteobacteria",
                "Gammaproteobacteria", "Enterobacterales",
                "Enterobacteriaceae", "Escherichia", "E_coli",
                0.95, 10, 11),
            new("read2", "Bacteria", "Firmicutes",
                "Bacilli", "Bacillales",
                "Bacillaceae", "Bacillus", "B_subtilis",
                0.90, 8, 10)
        };

        var profile = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifications);

        return Verify(new
        {
            profile.TotalReads,
            profile.ClassifiedReads,
            Shannon = Math.Round(profile.ShannonDiversity, 4),
            Simpson = Math.Round(profile.SimpsonDiversity, 4)
        });
    }

    [Test]
    public Task BinContigs_KnownContigs_MatchesSnapshot()
    {
        var contigs = Enumerable.Range(1, 20).Select(i =>
        {
            var seq = string.Concat(Enumerable.Repeat(i % 2 == 0 ? "ACGT" : "GCGC", 500));
            return ($"contig_{i}", seq, (double)(10 + i * 2));
        });

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 3, minBinSize: 100)
            .Select(b => new { b.BinId, ContigCount = b.ContigIds.Count, b.TotalLength, GcContent = Math.Round(b.GcContent, 4) })
            .ToList();

        return Verify(new { Bins = bins });
    }
}
