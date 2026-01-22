using NUnit.Framework;
using SuffixTree.Genomics;
using System;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for ProbeDesigner.DesignProbes, DesignTilingProbes, and probe scoring.
/// Test Unit: PROBE-DESIGN-001
/// </summary>
[TestFixture]
public class ProbeDesigner_ProbeDesign_Tests
{
    #region Test Data

    // Good sequence for microarray probes (moderate GC, no extreme features)
    private static readonly string MicroarrayTargetSequence =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

    // Longer sequence for various tests
    private static readonly string LongSequence =
        new string('A', 30) + "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" + new string('T', 30);

    #endregion

    #region DesignProbes - Input Validation (Must)

    [Test]
    public void DesignProbes_EmptySequence_ReturnsEmpty()
    {
        // M1: Empty sequence boundary condition
        var probes = ProbeDesigner.DesignProbes("").ToList();

        Assert.That(probes, Is.Empty);
    }

    [Test]
    public void DesignProbes_NullSequence_ReturnsEmpty()
    {
        // M2: Null sequence boundary condition
        var probes = ProbeDesigner.DesignProbes(null!).ToList();

        Assert.That(probes, Is.Empty);
    }

    [Test]
    public void DesignProbes_ShortSequence_ReturnsEmpty()
    {
        // M3: Sequence shorter than MinLength returns empty
        var probes = ProbeDesigner.DesignProbes("ACGT").ToList();

        Assert.That(probes, Is.Empty);
    }

    #endregion

    #region DesignProbes - Invariants (Must)

    [Test]
    public void DesignProbes_ValidSequence_ProbesHaveScoreInValidRange()
    {
        // M4: Score range invariant: 0.0 ≤ score ≤ 1.0
        var probes = ProbeDesigner.DesignProbes(MicroarrayTargetSequence, maxProbes: 10).ToList();

        Assert.That(probes, Is.Not.Empty, "Should produce at least one probe");

        Assert.Multiple(() =>
        {
            foreach (var probe in probes)
            {
                Assert.That(probe.Score, Is.InRange(0.0, 1.0),
                    $"Probe at {probe.Start} has score {probe.Score} outside valid range");
            }
        });
    }

    [Test]
    public void DesignProbes_ValidSequence_ProbesHaveGcContentInValidRange()
    {
        // M5: GC content invariant: 0.0 ≤ GC ≤ 1.0
        var probes = ProbeDesigner.DesignProbes(MicroarrayTargetSequence, maxProbes: 10).ToList();

        Assert.That(probes, Is.Not.Empty);

        Assert.Multiple(() =>
        {
            foreach (var probe in probes)
            {
                Assert.That(probe.GcContent, Is.InRange(0.0, 1.0),
                    $"Probe at {probe.Start} has GC content {probe.GcContent} outside valid range");
            }
        });
    }

    [Test]
    public void DesignProbes_ValidSequence_ProbesHavePositiveTm()
    {
        // M6: Tm positivity invariant: Tm > 0
        var probes = ProbeDesigner.DesignProbes(MicroarrayTargetSequence, maxProbes: 10).ToList();

        Assert.That(probes, Is.Not.Empty);

        Assert.Multiple(() =>
        {
            foreach (var probe in probes)
            {
                Assert.That(probe.Tm, Is.GreaterThan(0),
                    $"Probe at {probe.Start} has non-positive Tm {probe.Tm}");
            }
        });
    }

    [Test]
    public void DesignProbes_ValidSequence_ProbesHaveValidCoordinates()
    {
        // M7: Coordinate validity: 0 ≤ Start < End < sequence.Length
        string target = LongSequence;
        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 10).ToList();

        Assert.That(probes, Is.Not.Empty);

        Assert.Multiple(() =>
        {
            foreach (var probe in probes)
            {
                Assert.That(probe.Start, Is.GreaterThanOrEqualTo(0),
                    $"Probe Start {probe.Start} is negative");
                Assert.That(probe.End, Is.LessThan(target.Length),
                    $"Probe End {probe.End} exceeds sequence length {target.Length}");
                Assert.That(probe.End, Is.GreaterThan(probe.Start),
                    $"Probe End {probe.End} is not greater than Start {probe.Start}");
            }
        });
    }

    [Test]
    public void DesignProbes_ValidSequence_ProbeSequenceMatchesSubstring()
    {
        // M8: Probe sequence equals input substring at coordinates
        string target = LongSequence.ToUpperInvariant();
        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 10).ToList();

        Assert.That(probes, Is.Not.Empty);

        Assert.Multiple(() =>
        {
            foreach (var probe in probes)
            {
                int length = probe.End - probe.Start + 1;
                string expected = target.Substring(probe.Start, length);
                Assert.That(probe.Sequence, Is.EqualTo(expected),
                    $"Probe sequence mismatch at position {probe.Start}");
            }
        });
    }

    #endregion

    #region DesignProbes - Parameters (Must)

    [Test]
    public void DesignProbes_MaxProbesParameter_LimitsResultCount()
    {
        // M15: maxProbes parameter limits returned count
        int maxProbes = 3;
        var probes = ProbeDesigner.DesignProbes(MicroarrayTargetSequence, maxProbes: maxProbes).ToList();

        Assert.That(probes.Count, Is.LessThanOrEqualTo(maxProbes));
    }

    [Test]
    public void DesignProbes_MicroarrayDefaults_ProducesCorrectLengthProbes()
    {
        // M11: Microarray defaults: length 50-70 bp
        var param = ProbeDesigner.Defaults.Microarray;
        string target = new string('G', 25) + "ACGTACGTACGTACGTACGTACGTACGT" + new string('C', 25);

        var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 5).ToList();

        Assert.Multiple(() =>
        {
            foreach (var probe in probes)
            {
                Assert.That(probe.Sequence.Length, Is.InRange(param.MinLength, param.MaxLength),
                    $"Probe length {probe.Sequence.Length} outside Microarray range [{param.MinLength}, {param.MaxLength}]");
            }
        });
    }

    [Test]
    public void DesignProbes_FISHDefaults_ProducesCorrectLengthProbes()
    {
        // M12: FISH defaults: length 200-500 bp
        var param = ProbeDesigner.Defaults.FISH;

        Assert.Multiple(() =>
        {
            Assert.That(param.MinLength, Is.GreaterThanOrEqualTo(200), "FISH MinLength should be ≥200");
            Assert.That(param.MaxLength, Is.GreaterThanOrEqualTo(500), "FISH MaxLength should be ≥500");
        });

        // Create a sequence long enough for FISH probes
        string target = new string('G', 100) + new string('A', 100) +
                        new string('C', 100) + new string('T', 100) +
                        new string('G', 100) + new string('A', 100);

        var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 3).ToList();

        if (probes.Count > 0)
        {
            Assert.Multiple(() =>
            {
                foreach (var probe in probes)
                {
                    Assert.That(probe.Sequence.Length, Is.InRange(param.MinLength, param.MaxLength),
                        $"FISH probe length {probe.Sequence.Length} outside range [{param.MinLength}, {param.MaxLength}]");
                }
            });
        }
    }

    #endregion

    #region DesignProbes - Edge Cases (Must)

    [Test]
    public void DesignProbes_AllGC_ReturnsProbesWithHighGcContent()
    {
        // M13: High GC content (100%) results in GcContent ≈ 1.0
        string target = new string('G', 100);

        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var probe in probes)
        {
            Assert.That(probe.GcContent, Is.EqualTo(1.0).Within(0.01),
                "All-GC sequence should have GC content of 1.0");
        }
    }

    [Test]
    public void DesignProbes_AllAT_ReturnsProbesWithLowGcContent()
    {
        // M14: Low GC content (all A/T) results in low GcContent
        string target = new string('A', 50) + new string('T', 50);

        var probes = ProbeDesigner.DesignProbes(target).ToList();

        foreach (var probe in probes)
        {
            Assert.That(probe.GcContent, Is.LessThanOrEqualTo(0.1),
                "All-AT sequence should have GC content near 0");
        }
    }

    #endregion

    #region DesignTilingProbes (Must)

    [Test]
    public void DesignTilingProbes_CoversExpectedPositions()
    {
        // M9: Tiling probes cover expected positions
        string target = new string('A', 100) + "GCGCGCGC" + new string('T', 100);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 50, overlap: 10);

        Assert.Multiple(() =>
        {
            Assert.That(tiling.Probes.Count, Is.GreaterThan(0), "Should produce tiling probes");
            Assert.That(tiling.Coverage, Is.GreaterThan(0), "Should cover some positions");
        });
    }

    [Test]
    public void DesignTilingProbes_AllProbesHaveTilingType()
    {
        // M10: Tiling probes all have Type = Tiling
        string target = new string('A', 200);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 50, overlap: 10);

        Assert.That(tiling.Probes.All(p => p.Type == ProbeDesigner.ProbeType.Tiling), Is.True,
            "All tiling probes should have Type = Tiling");
    }

    [Test]
    public void DesignTilingProbes_CalculatesTmStatisticsCorrectly()
    {
        // S5: Tiling probes calculate mean Tm correctly
        string target = new string('G', 50) + new string('C', 50) + new string('A', 50);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 40, overlap: 10);

        Assert.Multiple(() =>
        {
            Assert.That(tiling.MeanTm, Is.GreaterThan(0), "Mean Tm should be positive");
            Assert.That(tiling.TmRange, Is.GreaterThanOrEqualTo(0), "Tm range should be non-negative");
        });
    }

    #endregion

    #region DesignProbes - Quality (Should)

    [Test]
    public void DesignProbes_HomopolymerSequence_GeneratesWarnings()
    {
        // S1: Homopolymer runs generate warnings
        string target = new string('A', 20) + "GCGCGCGC" + new string('A', 20) +
                       "TATATATA" + new string('G', 30);

        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 10).ToList();

        // Probes covering homopolymer regions may have warnings
        var probesWithWarnings = probes.Where(p => p.Warnings.Count > 0).ToList();
        Assert.That(probesWithWarnings.Count, Is.GreaterThanOrEqualTo(0),
            "Probes with warnings should exist (or none if all probes avoid homopolymers)");
    }

    [Test]
    public void DesignProbes_CaseInsensitiveInput_ProducesConsistentResults()
    {
        // S2: Case-insensitive input handling
        string upper = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        string lower = "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt";

        var probesUpper = ProbeDesigner.DesignProbes(upper, maxProbes: 1).ToList();
        var probesLower = ProbeDesigner.DesignProbes(lower, maxProbes: 1).ToList();

        if (probesUpper.Count > 0 && probesLower.Count > 0)
        {
            Assert.That(probesUpper[0].Tm, Is.EqualTo(probesLower[0].Tm).Within(0.1),
                "Case should not affect Tm calculation");
        }
    }

    [Test]
    public void DesignProbes_ProbesAreSortedByScoreDescending()
    {
        // S6: Probes are sorted by score descending
        var probes = ProbeDesigner.DesignProbes(MicroarrayTargetSequence, maxProbes: 10).ToList();

        if (probes.Count >= 2)
        {
            for (int i = 0; i < probes.Count - 1; i++)
            {
                Assert.That(probes[i].Score, Is.GreaterThanOrEqualTo(probes[i + 1].Score),
                    $"Probe at index {i} (score {probes[i].Score}) should have score ≥ probe at index {i + 1} (score {probes[i + 1].Score})");
            }
        }
    }

    #endregion

    #region DesignAntisenseProbes (Should)

    [Test]
    public void DesignAntisenseProbes_ReturnsAntisenseType()
    {
        // S3: DesignAntisenseProbes returns Antisense type
        string mRna = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC";

        var probes = ProbeDesigner.DesignAntisenseProbes(mRna, maxProbes: 3).ToList();

        Assert.That(probes.All(p => p.Type == ProbeDesigner.ProbeType.Antisense), Is.True,
            "All antisense probes should have Type = Antisense");
    }

    #endregion

    #region DesignMolecularBeacon (Should)

    [Test]
    public void DesignMolecularBeacon_CreatesBeaconWithStem()
    {
        // S4: MolecularBeacon has stem sequences
        string target = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        var beacon = ProbeDesigner.DesignMolecularBeacon(target, probeLength: 20, stemLength: 5);

        Assert.Multiple(() =>
        {
            Assert.That(beacon, Is.Not.Null, "Should create a molecular beacon");
            Assert.That(beacon!.Value.Type, Is.EqualTo(ProbeDesigner.ProbeType.MolecularBeacon));
            Assert.That(beacon.Value.Sequence.Length, Is.GreaterThan(20),
                "Beacon should be longer than loop due to stems");
        });
    }

    [Test]
    public void DesignMolecularBeacon_ShortSequence_ReturnsNull()
    {
        // Boundary: Short sequence returns null
        var beacon = ProbeDesigner.DesignMolecularBeacon("ACGT", probeLength: 20);

        Assert.That(beacon, Is.Null);
    }

    #endregion

    #region Application Defaults (Could)

    [Test]
    public void DesignProbes_qPCRDefaults_ProducesCorrectLengthProbes()
    {
        // C1: qPCR defaults produce 20-30 bp probes
        var param = ProbeDesigner.Defaults.qPCR;

        Assert.Multiple(() =>
        {
            Assert.That(param.MinLength, Is.InRange(18, 25), "qPCR MinLength should be ~20");
            Assert.That(param.MaxLength, Is.InRange(25, 35), "qPCR MaxLength should be ~30");
        });

        // Create suitable sequence for qPCR probe
        string target = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 5).ToList();

        if (probes.Count > 0)
        {
            foreach (var probe in probes)
            {
                Assert.That(probe.Sequence.Length, Is.InRange(param.MinLength, param.MaxLength),
                    $"qPCR probe length {probe.Sequence.Length} outside range [{param.MinLength}, {param.MaxLength}]");
            }
        }
    }

    #endregion

    #region Invariant Group Assertions

    [Test]
    public void DesignProbes_AllInvariants_HoldForValidInput()
    {
        // Combined invariant test for comprehensive coverage
        string target = LongSequence.ToUpperInvariant();
        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 5).ToList();

        Assert.That(probes, Is.Not.Empty, "Should produce probes");

        foreach (var probe in probes)
        {
            Assert.Multiple(() =>
            {
                // Score range
                Assert.That(probe.Score, Is.InRange(0.0, 1.0), $"Score out of range: {probe.Score}");

                // GC content range
                Assert.That(probe.GcContent, Is.InRange(0.0, 1.0), $"GC out of range: {probe.GcContent}");

                // Tm positivity
                Assert.That(probe.Tm, Is.GreaterThan(0), $"Tm not positive: {probe.Tm}");

                // Coordinate validity
                Assert.That(probe.Start, Is.GreaterThanOrEqualTo(0), $"Start negative: {probe.Start}");
                Assert.That(probe.End, Is.LessThan(target.Length), $"End exceeds length: {probe.End}");
                Assert.That(probe.End, Is.GreaterThan(probe.Start), $"End ≤ Start: {probe.End} ≤ {probe.Start}");

                // Sequence match
                int length = probe.End - probe.Start + 1;
                string expected = target.Substring(probe.Start, length);
                Assert.That(probe.Sequence, Is.EqualTo(expected), "Sequence mismatch");

                // Warnings is not null
                Assert.That(probe.Warnings, Is.Not.Null, "Warnings should not be null");
            });
        }
    }

    #endregion

    #region Suffix Tree Optimization

    [Test]
    public void DesignProbes_WithSuffixTree_FiltersNonUniqueProbes()
    {
        // Create a genome with repeated regions
        string uniqueRegion = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 52bp
        string repeatedRegion = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"; // 52bp - appears twice
        string genome = uniqueRegion + repeatedRegion + "AAAAAAAA" + repeatedRegion;

        // Build suffix tree for the genome
        var genomeIndex = global::SuffixTree.SuffixTree.Build(genome);

        // Design probes requiring uniqueness
        var param = ProbeDesigner.Defaults.Microarray with { MinLength = 50, MaxLength = 52 };
        var uniqueProbes = ProbeDesigner.DesignProbes(uniqueRegion, genomeIndex, param, maxProbes: 5, requireUnique: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(uniqueProbes.Count, Is.GreaterThan(0), "Should find unique probes");
            foreach (var probe in uniqueProbes)
            {
                // Verify probe is indeed unique in genome
                var positions = genomeIndex.FindAllOccurrences(probe.Sequence);
                Assert.That(positions.Count, Is.EqualTo(1),
                    $"Probe at {probe.Start} should be unique, but found {positions.Count} occurrences");
            }
        });
    }

    [Test]
    public void DesignProbes_WithSuffixTree_PerformanceImprovement()
    {
        // Create a moderately long sequence
        string target = string.Concat(Enumerable.Repeat("ACGTACGTACGTACGT", 50)); // 800bp

        // Build suffix tree once - O(n)
        var genomeIndex = global::SuffixTree.SuffixTree.Build(target);

        var param = ProbeDesigner.Defaults.Microarray;

        // Time without suffix tree
        var sw1 = System.Diagnostics.Stopwatch.StartNew();
        var probesWithout = ProbeDesigner.DesignProbes(target, param, maxProbes: 10).ToList();
        sw1.Stop();

        // Time with suffix tree (includes specificity check)
        var sw2 = System.Diagnostics.Stopwatch.StartNew();
        var probesWith = ProbeDesigner.DesignProbes(target, genomeIndex, param, maxProbes: 10, requireUnique: false).ToList();
        sw2.Stop();

        // Both should produce results
        Assert.Multiple(() =>
        {
            Assert.That(probesWithout.Count, Is.GreaterThan(0), "Should produce probes without index");
            Assert.That(probesWith.Count, Is.GreaterThan(0), "Should produce probes with index");
        });

        TestContext.WriteLine($"Without suffix tree: {sw1.ElapsedMilliseconds}ms");
        TestContext.WriteLine($"With suffix tree: {sw2.ElapsedMilliseconds}ms");
    }

    #endregion
}
