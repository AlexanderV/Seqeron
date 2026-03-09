using NUnit.Framework;
using Seqeron.Genomics;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Comprehensive tests for MetagenomicsAnalyzer.BinContigs (Genome Binning).
/// 
/// Test Unit: META-BIN-001
/// Algorithm: Metagenomic Genome Binning (MAG assembly)
/// 
/// Evidence Sources:
/// - Wikipedia: Binning (metagenomics)
/// - Teeling et al. (2004): TETRA tetranucleotide analysis
/// - Parks et al. (2014): CheckM quality metrics
/// - Maguire et al. (2020): MAG binning limitations
/// 
/// Key Concepts:
/// - Compositional binning based on GC content, tetranucleotide frequencies, and coverage
/// - K-means clustering on composite distance (GC + coverage + TNF Pearson)
/// - Quality metrics: Completeness (length / expected genome size) and Contamination (GC std dev)
/// - Minimum bin size filtering
/// </summary>
[TestFixture]
[Category("Metagenomics")]
[Category("META-BIN-001")]
public class MetagenomicsAnalyzer_GenomeBinning_Tests
{
    #region Test Data Helpers

    /// <summary>
    /// Creates a high-GC contig (approximately 100% G+C).
    /// Used to test GC-based separation.
    /// </summary>
    private static string CreateHighGcSequence(int length)
    {
        // Alternating G and C for 100% GC content
        var sb = new System.Text.StringBuilder(length);
        for (int i = 0; i < length; i++)
        {
            sb.Append(i % 2 == 0 ? 'G' : 'C');
        }
        return sb.ToString();
    }

    /// <summary>
    /// Creates a low-GC contig (approximately 0% G+C).
    /// Used to test GC-based separation.
    /// </summary>
    private static string CreateLowGcSequence(int length)
    {
        // Alternating A and T for 0% GC content
        var sb = new System.Text.StringBuilder(length);
        for (int i = 0; i < length; i++)
        {
            sb.Append(i % 2 == 0 ? 'A' : 'T');
        }
        return sb.ToString();
    }

    /// <summary>
    /// Creates a sequence with approximately 50% GC content.
    /// </summary>
    private static string CreateMidGcSequence(int length)
    {
        // Repeating GCTA for ~50% GC
        var sb = new System.Text.StringBuilder(length);
        string pattern = "GCTA";
        for (int i = 0; i < length; i++)
        {
            sb.Append(pattern[i % 4]);
        }
        return sb.ToString();
    }

    #endregion

    #region Must Tests - Empty and Null Handling

    /// <summary>
    /// M1: Empty input must return empty result.
    /// Source: Standard edge case, implementation contract.
    /// </summary>
    [Test]
    public void BinContigs_EmptyInput_ReturnsEmpty()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();

        var bins = MetagenomicsAnalyzer.BinContigs(contigs).ToList();

        Assert.That(bins, Is.Empty);
    }

    /// <summary>
    /// M2: Single contig below minimum size returns empty.
    /// Source: minBinSize parameter contract.
    /// </summary>
    [Test]
    public void BinContigs_SingleContigBelowMinSize_ReturnsEmpty()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("contig1", CreateMidGcSequence(1000), 10.0) // 1000 bp < 500000 bp default
        };

        var bins = MetagenomicsAnalyzer.BinContigs(contigs).ToList();

        Assert.That(bins, Is.Empty);
    }

    #endregion

    #region Must Tests - Basic Functionality

    /// <summary>
    /// M3: Valid contigs return non-empty bins.
    /// Source: Wikipedia - Binning produces MAGs from contigs.
    /// 50 contigs × 20kb = 1MB total, well above minBinSize=100kb.
    /// </summary>
    [Test]
    public void BinContigs_ValidContigs_ReturnsNonEmptyBins()
    {
        // Create contigs large enough to pass minimum size filter
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();
        for (int i = 0; i < 50; i++)
        {
            contigs.Add(($"contig_{i}", CreateHighGcSequence(20000), 10.0));
        }

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 5, minBinSize: 100000).ToList();

        Assert.That(bins, Is.Not.Empty,
            "50 contigs × 20kb = 1MB total with 100kb minBinSize must produce at least one bin");
    }

    /// <summary>
    /// M4: All bin IDs must be unique.
    /// Source: Data integrity invariant.
    /// </summary>
    [Test]
    public void BinContigs_BinIds_AreUnique()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();
        // Create diverse contigs to generate multiple bins - 50 high GC + 50 low GC
        for (int i = 0; i < 100; i++)
        {
            string seq = i < 50 ? CreateHighGcSequence(15000) : CreateLowGcSequence(15000);
            contigs.Add(($"contig_{i}", seq, 10.0 + i));
        }

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 10, minBinSize: 50000).ToList();

        // 100 contigs × 15kb = 1.5Mb total, minBinSize=50kb → MUST produce bins
        Assert.That(bins, Is.Not.Empty,
            "100 diverse 15kb contigs with 50kb min bin size MUST produce at least one bin");

        var binIds = bins.Select(b => b.BinId).ToList();
        Assert.That(binIds.Distinct().Count(), Is.EqualTo(binIds.Count),
            "All bin IDs must be unique");
    }

    #endregion

    #region Must Tests - Quality Metric Ranges

    /// <summary>
    /// M5: Completeness must be in valid range [0, 100] and equal min(totalLength / expectedGenomeSize × 100, 100).
    /// Source: CheckM standard - Completeness is percentage.
    /// Formula: completeness = min(bin.TotalLength / expectedGenomeSize × 100, 100)
    /// Default expectedGenomeSize = 4,000,000 bp.
    /// </summary>
    [Test]
    public void BinContigs_Completeness_InValidRange()
    {
        // 10 contigs × 200kb = 2MB total; expectedGenomeSize = 4MB (default)
        // numBins=1 forces all contigs into a single bin for deterministic verification
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();
        for (int i = 0; i < 10; i++)
        {
            contigs.Add(($"contig_{i}", CreateMidGcSequence(200000), 15.0));
        }

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 1, minBinSize: 100000).ToList();

        Assert.That(bins, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(bins[0].Completeness, Is.GreaterThanOrEqualTo(0));
            Assert.That(bins[0].Completeness, Is.LessThanOrEqualTo(100));
            // 2,000,000 / 4,000,000 × 100 = 50.0
            Assert.That(bins[0].Completeness, Is.EqualTo(50.0),
                "Completeness = min(totalLength / expectedGenomeSize × 100, 100) = 2M/4M × 100 = 50.0");
        });
    }

    /// <summary>
    /// M6: Contamination must be in valid range [0, 100].
    /// Source: CheckM standard - Contamination is percentage.
    /// Formula: contamination = min(stddev(GC values) / 0.5 × 100, 100)
    /// Theoretical max stddev for values in [0,1] is 0.5 (half at 0, half at 1).
    /// </summary>
    [Test]
    public void BinContigs_Contamination_InValidRange()
    {
        // Scenario 1: Uniform GC → stddev = 0 → contamination = 0
        var uniformContigs = new List<(string ContigId, string Sequence, double Coverage)>();
        for (int i = 0; i < 10; i++)
            uniformContigs.Add(($"contig_{i}", CreateMidGcSequence(200000), 15.0));

        var uniformBins = MetagenomicsAnalyzer.BinContigs(uniformContigs, numBins: 1, minBinSize: 100000).ToList();

        Assert.That(uniformBins, Has.Count.EqualTo(1));
        Assert.That(uniformBins[0].Contamination, Is.EqualTo(0.0),
            "All contigs same GC → zero within-bin GC stddev → contamination = 0");

        // Scenario 2: Maximum GC variance → stddev = 0.5 → contamination = 100
        // Equal mix of GC=1.0 (GCGCGC...) and GC=0.0 (ATATAT...)
        var mixedContigs = new List<(string ContigId, string Sequence, double Coverage)>();
        for (int i = 0; i < 10; i++)
            mixedContigs.Add(($"high_{i}", CreateHighGcSequence(200000), 15.0));
        for (int i = 0; i < 10; i++)
            mixedContigs.Add(($"low_{i}", CreateLowGcSequence(200000), 15.0));

        var mixedBins = MetagenomicsAnalyzer.BinContigs(mixedContigs, numBins: 1, minBinSize: 100000).ToList();

        Assert.That(mixedBins, Has.Count.EqualTo(1));
        Assert.That(mixedBins[0].Contamination, Is.EqualTo(100.0),
            "Equal mix of GC=1.0 and GC=0.0 → stddev=0.5 (theoretical max) → contamination=100");
    }

    /// <summary>
    /// M7: GC content must be in valid range [0, 1] and match known values for pure bins.
    /// Source: Mathematical definition - GC = (G+C)/(A+T+G+C).
    /// Pure high-GC bin (all GCGCGC...) → GC = 1.0; pure low-GC bin (all ATATAT...) → GC = 0.0.
    /// </summary>
    [Test]
    public void BinContigs_GcContent_InValidRange()
    {
        // Use extreme GC populations that k-means will separate into pure bins
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();
        for (int i = 0; i < 20; i++)
            contigs.Add(($"highgc_{i}", CreateHighGcSequence(25000), 10.0));
        for (int i = 0; i < 20; i++)
            contigs.Add(($"lowgc_{i}", CreateLowGcSequence(25000), 10.0));

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 10, minBinSize: 50000).ToList();

        Assert.Multiple(() =>
        {
            foreach (var bin in bins)
            {
                // Range invariant: GC ∈ [0, 1]
                Assert.That(bin.GcContent, Is.GreaterThanOrEqualTo(0),
                    $"Bin {bin.BinId} GC must be >= 0");
                Assert.That(bin.GcContent, Is.LessThanOrEqualTo(1),
                    $"Bin {bin.BinId} GC must be <= 1");

                // Exact value check for pure bins
                bool allHigh = bin.ContigIds.All(id => id.StartsWith("highgc_"));
                bool allLow = bin.ContigIds.All(id => id.StartsWith("lowgc_"));

                if (allHigh)
                    Assert.That(bin.GcContent, Is.EqualTo(1.0),
                        $"Bin {bin.BinId}: pure high-GC bin (all GCGCGC...) must have GC = 1.0");
                else if (allLow)
                    Assert.That(bin.GcContent, Is.EqualTo(0.0),
                        $"Bin {bin.BinId}: pure low-GC bin (all ATATAT...) must have GC = 0.0");
            }
        });
    }

    #endregion

    #region Must Tests - Data Integrity

    /// <summary>
    /// M8: Each bin's TotalLength must equal the sum of its contig lengths.
    /// Source: Mathematical invariant — bin.TotalLength = Σ len(contig) for contig ∈ bin.
    /// </summary>
    [Test]
    public void BinContigs_TotalLength_EqualsContigLengthSum()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("c1", CreateHighGcSequence(200000), 10.0),
            ("c2", CreateHighGcSequence(150000), 10.0),
            ("c3", CreateHighGcSequence(250000), 10.0)
        };

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 3, minBinSize: 100000).ToList();

        var contigLengths = contigs.ToDictionary(c => c.ContigId, c => (double)c.Sequence.Length);

        Assert.Multiple(() =>
        {
            foreach (var bin in bins)
            {
                double expectedLength = bin.ContigIds.Sum(id => contigLengths[id]);
                Assert.That(bin.TotalLength, Is.EqualTo(expectedLength),
                    $"Bin {bin.BinId}: TotalLength must equal sum of its contig lengths");
            }

            // Total binned length cannot exceed input
            double totalBinLength = bins.Sum(b => b.TotalLength);
            double totalContigLength = contigs.Sum(c => c.Sequence.Length);
            Assert.That(totalBinLength, Is.LessThanOrEqualTo(totalContigLength),
                "Total binned length cannot exceed input length");
        });
    }

    /// <summary>
    /// M9: Minimum bin size filtering works correctly.
    /// Source: Implementation specification.
    /// </summary>
    [Test]
    public void BinContigs_MinBinSize_FiltersSmallBins()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("large1", CreateHighGcSequence(300000), 10.0),
            ("large2", CreateHighGcSequence(300000), 10.0),
            ("small1", CreateLowGcSequence(10000), 10.0)  // This alone would be below threshold
        };

        double minBinSize = 200000;
        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 5, minBinSize: minBinSize).ToList();

        Assert.Multiple(() =>
        {
            foreach (var bin in bins)
            {
                Assert.That(bin.TotalLength, Is.GreaterThanOrEqualTo(minBinSize),
                    $"Bin {bin.BinId} should be >= minBinSize ({minBinSize})");
            }
        });
    }

    /// <summary>
    /// M11: Contig IDs are preserved in bins; bins are disjoint.
    /// Source: Data integrity requirement.
    /// </summary>
    [Test]
    public void BinContigs_ContigIds_PreservedInBins()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("contig_alpha", CreateHighGcSequence(200000), 10.0),
            ("contig_beta", CreateHighGcSequence(200000), 12.0),
            ("contig_gamma", CreateHighGcSequence(200000), 11.0)
        };

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 3, minBinSize: 100000).ToList();

        var allBinnedContigIds = bins.SelectMany(b => b.ContigIds).ToList();
        var inputContigIds = contigs.Select(c => c.ContigId).ToHashSet();

        Assert.Multiple(() =>
        {
            // All binned contigs must come from input
            Assert.That(allBinnedContigIds.All(id => inputContigIds.Contains(id)), Is.True,
                "All binned contig IDs must come from input");

            // Bins must be disjoint: no contig appears in more than one bin
            Assert.That(allBinnedContigIds.Count, Is.EqualTo(allBinnedContigIds.Distinct().Count()),
                "Each contig must appear in at most one bin (bins are disjoint)");
        });
    }

    /// <summary>
    /// M12: Bin coverage equals the arithmetic mean of its constituent contig coverages.
    /// Source: Mathematical definition — bin coverage = Σcoverage / n for contigs in bin.
    /// </summary>
    [Test]
    public void BinContigs_Coverage_PreservedInBins()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("c1", CreateHighGcSequence(200000), 20.0),
            ("c2", CreateHighGcSequence(200000), 30.0),
            ("c3", CreateHighGcSequence(200000), 25.0)
        };

        // numBins=1 forces all contigs into a single bin
        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 1, minBinSize: 100000).ToList();

        Assert.That(bins, Has.Count.EqualTo(1),
            "numBins=1 with 600kb total (> 100kb minBinSize) must produce exactly one bin");

        // Coverage = average(20.0, 30.0, 25.0) = 75.0 / 3 = 25.0
        Assert.That(bins[0].Coverage, Is.EqualTo(25.0).Within(1e-10),
            "Bin coverage must equal the arithmetic mean of its contig coverages");
    }

    #endregion

    #region Must Tests - GC-Based Separation

    /// <summary>
    /// M10: High GC and low GC contigs must separate into different bins.
    /// Source: Wikipedia - GC content is primary compositional feature.
    /// Source: TETRA paper - Compositional signals distinguish genomes.
    /// K-means on GC + coverage + TNF Pearson distance must separate extreme GC populations.
    /// </summary>
    [Test]
    public void BinContigs_HighGcVsLowGc_SeparatesIntoDifferentBins()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();

        // High GC contigs (100% GC)
        for (int i = 0; i < 30; i++)
            contigs.Add(($"highgc_{i}", CreateHighGcSequence(25000), 10.0));

        // Low GC contigs (0% GC)
        for (int i = 0; i < 30; i++)
            contigs.Add(($"lowgc_{i}", CreateLowGcSequence(25000), 10.0));

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 10, minBinSize: 100000).ToList();

        Assert.That(bins.Count, Is.GreaterThanOrEqualTo(2),
            "Extreme GC difference (0% vs 100%) must produce at least 2 bins");

        // No bin may mix high-GC and low-GC contigs
        Assert.Multiple(() =>
        {
            foreach (var bin in bins)
            {
                bool hasHigh = bin.ContigIds.Any(id => id.StartsWith("highgc_"));
                bool hasLow = bin.ContigIds.Any(id => id.StartsWith("lowgc_"));
                Assert.That(hasHigh && hasLow, Is.False,
                    $"Bin {bin.BinId} mixes high-GC and low-GC contigs — must separate extreme GC");
            }
        });
    }

    #endregion

    #region Should Tests

    /// <summary>
    /// S1: Diverse input creates multiple bins.
    /// Source: Wikipedia - Binning separates different organisms.
    /// </summary>
    [Test]
    public void BinContigs_MultipleBins_FromDiverseInput()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();

        // Create contigs with varying GC content to encourage multiple bins
        for (int i = 0; i < 100; i++)
        {
            // Vary GC content by mixing sequences
            string seq = i < 33 ? CreateHighGcSequence(20000) :
                         i < 66 ? CreateMidGcSequence(20000) :
                                  CreateLowGcSequence(20000);
            contigs.Add(($"contig_{i}", seq, 10.0 + i % 20));
        }

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 10, minBinSize: 50000).ToList();

        // With diverse input (3 distinct GC groups), should get multiple bins
        Assert.That(bins.Count, Is.GreaterThanOrEqualTo(2),
            "Diverse input with 3 GC groups should produce at least 2 bins");
    }

    /// <summary>
    /// S2: NumBins parameter limits maximum bin count.
    /// Source: K-means with k clusters produces at most k non-empty groups.
    /// </summary>
    [Test]
    public void BinContigs_NumBins_AffectsMaximumBinCount()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();
        for (int i = 0; i < 100; i++)
        {
            contigs.Add(($"contig_{i}", CreateMidGcSequence(20000), 10.0 + i));
        }

        var binsSmall = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 3, minBinSize: 50000).ToList();
        var binsLarge = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 20, minBinSize: 50000).ToList();

        Assert.Multiple(() =>
        {
            // 100 contigs × 20kb = 2MB total > 50kb → must produce bins
            Assert.That(binsSmall, Is.Not.Empty, "Must produce bins from 2MB of contigs");
            Assert.That(binsLarge, Is.Not.Empty, "Must produce bins from 2MB of contigs");

            // K-means with k clusters produces at most k non-empty groups
            Assert.That(binsSmall.Count, Is.LessThanOrEqualTo(3),
                "Bin count must not exceed numBins=3");
            Assert.That(binsLarge.Count, Is.LessThanOrEqualTo(20),
                "Bin count must not exceed numBins=20");
        });
    }

    /// <summary>
    /// S3: Large dataset completes in reasonable time.
    /// Source: Practical requirement for metagenomics data.
    /// </summary>
    [Test]
    [CancelAfter(30000)] // 30 second timeout
    public void BinContigs_LargeDataset_CompletesWithinTimeout()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();

        // Create 500 contigs of 10kb each = 5MB total
        for (int i = 0; i < 500; i++)
        {
            contigs.Add(($"contig_{i}", CreateMidGcSequence(10000), 10.0 + i % 50));
        }

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 20, minBinSize: 100000).ToList();

        Assert.That(bins, Is.Not.Null, "Should complete without timeout");
    }

    #endregion

    #region Edge Cases

    /// <summary>
    /// C1: All contigs with same GC, coverage, and TNF should cluster together.
    /// Source: K-means with zero inter-point distance converges to minimal clusters.
    /// </summary>
    [Test]
    public void BinContigs_AllSameGc_ClustersTogether()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>();
        for (int i = 0; i < 50; i++)
            contigs.Add(($"contig_{i}", CreateMidGcSequence(15000), 10.0));

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 5, minBinSize: 50000).ToList();

        // 50 × 15kb = 750kb, minBinSize=50kb → must produce at least 1 bin
        Assert.That(bins.Count, Is.GreaterThanOrEqualTo(1),
            "Uniform contigs totaling 750kb with 50kb threshold must produce at least one bin");

        // Uniform features should concentrate contigs in few bins, not scatter
        Assert.That(bins.Count, Is.LessThanOrEqualTo(3),
            "Uniform GC/coverage/TNF input should not be split into many bins");
    }

    /// <summary>
    /// Very short sequences handled gracefully.
    /// </summary>
    [Test]
    public void BinContigs_VeryShortSequences_HandledGracefully()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("short1", "ACGT", 10.0),
            ("short2", "GC", 10.0),
            ("short3", "AT", 10.0)
        };

        // Should not throw, just return empty due to size filter
        Assert.DoesNotThrow(() =>
        {
            var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 3, minBinSize: 100).ToList();
        });
    }

    /// <summary>
    /// Zero coverage contigs handled correctly.
    /// </summary>
    [Test]
    public void BinContigs_ZeroCoverage_HandledGracefully()
    {
        var contigs = new List<(string ContigId, string Sequence, double Coverage)>
        {
            ("c1", CreateMidGcSequence(200000), 0.0),
            ("c2", CreateMidGcSequence(200000), 0.0)
        };

        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 3, minBinSize: 100000).ToList();

        Assert.That(bins, Is.Not.Null, "Should handle zero coverage");
    }

    #endregion
}
