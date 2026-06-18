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
        // Tiny taxonomy: root → E_coli, B_subtilis.
        var tree = new TaxonomyTree(new[]
        {
            new TaxonNode(1, "root", "root", 1),
            new TaxonNode(100, "E_coli", "species", 1),
            new TaxonNode(200, "B_subtilis", "species", 1),
        });

        // Build the canonical-k-mer → taxon-id database (k=31) from reference sequences so the
        // k-mer lengths are guaranteed consistent with the reads.
        var kmerDb = MetagenomicsAnalyzer.BuildKmerDatabase(new[]
        {
            (100, "ACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            (200, "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"),
        }, tree, k: 31);

        var reads = new[]
        {
            ("read1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
            ("read2", "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT"),
            ("read3", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
        };
        var classifications = MetagenomicsAnalyzer.ClassifyReads(reads, kmerDb, tree, k: 31)
            .Select(c => new { c.ReadId, c.TaxonId, c.TaxonName, c.RtlScore, c.Confidence, c.MatchedKmers, c.TotalKmers })
            .ToList();

        return Verify(new { Classifications = classifications });
    }

    [Test]
    public Task GenerateTaxonomicProfile_MatchesSnapshot()
    {
        // Create mock classifications using positional construction
        var classifications = new MetagenomicsAnalyzer.TaxonomicClassification[]
        {
            new("read1", TaxonId: 100, TaxonName: "E_coli", Rank: "species", RtlScore: 10,
                Confidence: 0.95, MatchedKmers: 10, TotalKmers: 11,
                Kingdom: "Bacteria", Phylum: "Proteobacteria",
                Class: "Gammaproteobacteria", Order: "Enterobacterales",
                Family: "Enterobacteriaceae", Genus: "Escherichia", Species: "E_coli"),
            new("read2", TaxonId: 200, TaxonName: "B_subtilis", Rank: "species", RtlScore: 8,
                Confidence: 0.90, MatchedKmers: 8, TotalKmers: 10,
                Kingdom: "Bacteria", Phylum: "Firmicutes",
                Class: "Bacilli", Order: "Bacillales",
                Family: "Bacillaceae", Genus: "Bacillus", Species: "B_subtilis")
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
