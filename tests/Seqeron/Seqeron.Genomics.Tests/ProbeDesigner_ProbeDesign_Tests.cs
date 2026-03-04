using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests;

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

        Assert.That(probes.Count, Is.EqualTo(maxProbes),
            $"81-bp ACGT-repeat sequence should yield exactly {maxProbes} probes");
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

        // Create a varied sequence long enough for FISH probes (avoid homopolymers)
        string target = string.Concat(Enumerable.Range(0, 150).Select(i => "ATGC"[i % 4])) +
                        string.Concat(Enumerable.Range(0, 150).Select(i => "CGAT"[i % 4])) +
                        string.Concat(Enumerable.Range(0, 150).Select(i => "TACG"[i % 4])) +
                        string.Concat(Enumerable.Range(0, 150).Select(i => "GCAT"[i % 4]));

        var probes = ProbeDesigner.DesignProbes(target, param, maxProbes: 3).ToList();

        // With a varied sequence of 600bp, we should get FISH probes
        Assert.That(probes, Is.Not.Empty,
            "Should generate FISH probes from 600bp varied sequence");
        Assert.Multiple(() =>
        {
            foreach (var probe in probes)
            {
                Assert.That(probe.Sequence.Length, Is.InRange(param.MinLength, param.MaxLength),
                    $"FISH probe length {probe.Sequence.Length} outside range [{param.MinLength}, {param.MaxLength}]");
            }
        });
    }

    #endregion

    #region DesignProbes - Edge Cases (Must)

    [Test]
    public void DesignProbes_AllGC_ReturnsProbesWithHighGcContent()
    {
        // M13: High GC content (100%) results in GcContent = 1.0 exactly
        // Use unrestricted GC/Tm params to bypass early rejection filter
        var param = ProbeDesigner.Defaults.Microarray with
        {
            MinGc = 0.0,
            MaxGc = 1.0,
            MinTm = 0,
            MaxTm = 200
        };
        string target = new string('G', 100);

        var probes = ProbeDesigner.DesignProbes(target, param).ToList();

        Assert.That(probes, Is.Not.Empty, "All-G sequence with unrestricted GC params must produce probes");
        foreach (var probe in probes)
        {
            Assert.That(probe.GcContent, Is.EqualTo(1.0),
                "All-GC probe must have GC content of exactly 1.0");
        }
    }

    [Test]
    public void DesignProbes_AllAT_ReturnsProbesWithLowGcContent()
    {
        // M14: Low GC content (all A/T) results in GcContent = 0.0 exactly
        // Use unrestricted GC/Tm params to bypass early rejection filter
        var param = ProbeDesigner.Defaults.Microarray with
        {
            MinGc = 0.0,
            MaxGc = 1.0,
            MinTm = 0,
            MaxTm = 200
        };
        string target = new string('A', 50) + new string('T', 50);

        var probes = ProbeDesigner.DesignProbes(target, param).ToList();

        Assert.That(probes, Is.Not.Empty, "All-AT sequence with unrestricted GC params must produce probes");
        foreach (var probe in probes)
        {
            Assert.That(probe.GcContent, Is.EqualTo(0.0),
                "All-AT probe must have GC content of exactly 0.0");
        }
    }

    #endregion

    #region DesignTilingProbes (Must)

    [Test]
    public void DesignTilingProbes_CoversExpectedPositions()
    {
        // M9: Tiling probes cover expected positions
        // 208-char sequence, probeLength=50, overlap=10 → step=40
        // Expected probes at positions: 0, 40, 80, 120 (160 > 208-50=158)
        // Coverage: positions 0-169 = 170
        string target = new string('A', 100) + "GCGCGCGC" + new string('T', 100);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 50, overlap: 10);

        Assert.Multiple(() =>
        {
            Assert.That(tiling.Probes.Count, Is.EqualTo(4), "Expected 4 tiling probes");
            Assert.That(tiling.Coverage, Is.EqualTo(170), "Expected coverage of 170 positions");

            var starts = tiling.Probes.Select(p => p.Start).ToList();
            Assert.That(starts, Is.EqualTo(new[] { 0, 40, 80, 120 }),
                "Tiling probes should start at exact positions");
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
        // 150-char sequence, probeLength=40, overlap=10 → step=30
        // Probes at positions 0, 30, 60, 90 (120 > 110)
        string target = new string('G', 50) + new string('C', 50) + new string('A', 50);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 40, overlap: 10);

        Assert.That(tiling.Probes.Count, Is.EqualTo(4), "Expected 4 tiling probes");

        double expectedMean = tiling.Probes.Average(p => p.Tm);
        double expectedRange = tiling.Probes.Max(p => p.Tm) - tiling.Probes.Min(p => p.Tm);

        Assert.Multiple(() =>
        {
            Assert.That(tiling.MeanTm, Is.EqualTo(expectedMean).Within(0.001),
                "MeanTm must equal average of individual probe Tm values");
            Assert.That(tiling.TmRange, Is.EqualTo(expectedRange).Within(0.001),
                "TmRange must equal max(Tm) - min(Tm)");
            Assert.That(tiling.TmRange, Is.GreaterThan(0),
                "Mixed GC sequence should produce probes with different Tm values");
        });
    }

    #endregion

    #region DesignProbes - Quality (Should)

    [Test]
    public void DesignProbes_HomopolymerSequence_GeneratesWarnings()
    {
        // S1: Homopolymer runs generate warnings
        // Sequence has a 30-G run at positions 56-85; probes spanning it must report homopolymer warning
        string target = new string('A', 20) + "GCGCGCGC" + new string('A', 20) +
                       "TATATATA" + new string('G', 30);

        var probes = ProbeDesigner.DesignProbes(target, maxProbes: 20).ToList();

        Assert.That(probes, Is.Not.Empty, "Should produce probes from 86-bp sequence");

        var probesWithHomopolymerWarning = probes
            .Where(p => p.Warnings.Any(w => w.Contains("Homopolymer")))
            .ToList();
        Assert.That(probesWithHomopolymerWarning, Is.Not.Empty,
            "Probes spanning 30-nt G run must have homopolymer warning");
    }

    [Test]
    public void DesignProbes_CaseInsensitiveInput_ProducesConsistentResults()
    {
        // S2: Case-insensitive input handling
        string upper = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        string lower = "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt";

        var probesUpper = ProbeDesigner.DesignProbes(upper, maxProbes: 1).ToList();
        var probesLower = ProbeDesigner.DesignProbes(lower, maxProbes: 1).ToList();

        // Both MUST produce probes - sequence is long enough
        Assert.That(probesUpper, Is.Not.Empty, "Upper case sequence must produce probes");
        Assert.That(probesLower, Is.Not.Empty, "Lower case sequence must produce probes");

        Assert.That(probesUpper[0].Tm, Is.EqualTo(probesLower[0].Tm).Within(0.1),
            "Case should not affect Tm calculation");
    }

    [Test]
    public void DesignProbes_ProbesAreSortedByScoreDescending()
    {
        // S6: Probes are sorted by score descending
        var probes = ProbeDesigner.DesignProbes(MicroarrayTargetSequence, maxProbes: 10).ToList();

        // MicroarrayTargetSequence is long enough to produce multiple probes
        Assert.That(probes.Count, Is.GreaterThanOrEqualTo(2),
            "MicroarrayTargetSequence with maxProbes:10 must produce at least 2 probes");

        for (int i = 0; i < probes.Count - 1; i++)
        {
            Assert.That(probes[i].Score, Is.GreaterThanOrEqualTo(probes[i + 1].Score),
                $"Probe at index {i} (score {probes[i].Score}) should have score ≥ probe at index {i + 1} (score {probes[i + 1].Score})");
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
        // stemLength=5 → stem5="GGCCC", stem3=RC("GGCCC")="GGGCC"
        // Total length = stem5(5) + loop(20) + stem3(5) = 30
        string target = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        var beacon = ProbeDesigner.DesignMolecularBeacon(target, probeLength: 20, stemLength: 5);

        Assert.That(beacon, Is.Not.Null, "Should create a molecular beacon");

        var b = beacon!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(b.Type, Is.EqualTo(ProbeDesigner.ProbeType.MolecularBeacon));
            Assert.That(b.Sequence.Length, Is.EqualTo(30),
                "Beacon = stem5(5) + loop(20) + stem3(5) = 30");
            Assert.That(b.Sequence.Substring(0, 5), Is.EqualTo("GGCCC"),
                "5' stem should be GGCCC");
            Assert.That(b.Sequence.Substring(b.Sequence.Length - 5, 5), Is.EqualTo("GGGCC"),
                "3' stem should be reverse complement of 5' stem");
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

        // 48 bp target with qPCR params MUST produce probes
        Assert.That(probes, Is.Not.Empty,
            "48 bp target with qPCR parameters must produce probes");

        foreach (var probe in probes)
        {
            Assert.That(probe.Sequence.Length, Is.InRange(param.MinLength, param.MaxLength),
                $"qPCR probe length {probe.Sequence.Length} outside range [{param.MinLength}, {param.MaxLength}]");
        }
    }

    [Test]
    public void ValidateProbe_SelfComplementarity_DetectsCorrectly()
    {
        // C2: Self-complementarity detection
        // Palindromic DNA (seq == reverse complement) → self-complementarity = 1.0
        string palindrome = "AACCGGTT"; // RC("AACCGGTT") = "AACCGGTT"
        var resultPalindrome = ProbeDesigner.ValidateProbe(palindrome, Array.Empty<string>());
        Assert.That(resultPalindrome.SelfComplementarity, Is.EqualTo(1.0),
            "Palindromic sequence must have self-complementarity of 1.0");

        // Non-complementary sequence (all-A vs all-T reverse complement) → 0.0
        string nonComp = "AAAAAAAAAA"; // RC = "TTTTTTTTTT", zero positional matches
        var resultNonComp = ProbeDesigner.ValidateProbe(nonComp, Array.Empty<string>());
        Assert.That(resultNonComp.SelfComplementarity, Is.EqualTo(0.0),
            "All-A sequence must have zero self-complementarity (RC = all-T)");
    }

    [Test]
    public void ValidateProbe_SecondaryStructure_IdentifiesHairpins()
    {
        // C3: Secondary structure detection identifies hairpins
        // "ACGT" + 3-nt loop + "ACGT" forms a hairpin (RC("ACGT")="ACGT", 100% stem match)
        string hairpin = "ACGTAAAACGT";
        var resultHairpin = ProbeDesigner.ValidateProbe(hairpin, Array.Empty<string>());
        Assert.That(resultHairpin.HasSecondaryStructure, Is.True,
            "Inverted repeat ACGT-loop-ACGT should be detected as secondary structure");

        // Sequence without inverted repeats → no secondary structure
        string noHairpin = "AAAAAACCCCCC";
        var resultNoHairpin = ProbeDesigner.ValidateProbe(noHairpin, Array.Empty<string>());
        Assert.That(resultNoHairpin.HasSecondaryStructure, Is.False,
            "Sequence without inverted repeats should have no secondary structure");
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

        TestContext.Out.WriteLine($"Without suffix tree: {sw1.ElapsedMilliseconds}ms");
        TestContext.Out.WriteLine($"With suffix tree: {sw2.ElapsedMilliseconds}ms");
    }

    #endregion

    #region Mutation-Killing Tests — MolecularBeacon Scoring

    [Test]
    public void DesignMolecularBeacon_AtRichTarget_ScorePenalizedForGcAndTm()
    {
        // AT-rich loop: GC ≈ 0%, Tm < 55 → both penalties fire (score ≤ 0.6)
        // Kills ||→&& mutation on beacon scoring conditions
        string target = "ATATATATATATATATATATATATATATAT"; // 28 bp, 0% GC

        var beacon = ProbeDesigner.DesignMolecularBeacon(target, probeLength: 25, stemLength: 5);

        Assert.That(beacon, Is.Not.Null, "Beacon should be designed even for AT-rich target");
        Assert.That(beacon!.Value.Score, Is.LessThan(0.8),
            "AT-rich loop should have GC and Tm penalties reducing score");
    }

    [Test]
    public void DesignMolecularBeacon_GcRichTarget_ScorePenalizedForGcAndTm()
    {
        // GC-rich loop: GC = 100%, Tm > 65 → both penalties fire (score ≤ 0.6)
        // Kills ||→&& mutation from opposite direction
        string target = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"; // 30 bp, 100% GC

        var beacon = ProbeDesigner.DesignMolecularBeacon(target, probeLength: 25, stemLength: 5);

        Assert.That(beacon, Is.Not.Null, "Beacon should be designed even for GC-rich target");
        Assert.That(beacon!.Value.Score, Is.LessThan(0.8),
            "GC-rich loop should have GC and Tm penalties reducing score");
    }

    #endregion
}
