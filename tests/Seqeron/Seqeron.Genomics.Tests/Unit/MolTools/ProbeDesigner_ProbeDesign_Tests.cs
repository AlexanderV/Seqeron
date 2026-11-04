namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for ProbeDesigner.DesignProbes, DesignTilingProbes, and probe scoring.
/// Test Unit: PROBE-DESIGN-001
/// </summary>
[TestFixture]
public class ProbeDesigner_ProbeDesign_Tests
{
    #region Test Data

    // Good sequence for microarray probes (moderate GC, no extreme features)
    private const string MicroarrayTargetSequence =
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
        // M11: Microarray defaults: length 50-60 bp
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

    [Test]
    public void DesignProbes_WindowWithN_GcIsFractionOfNonNBases_MatchesPrimer3()
    {
        // M16 (B07 audit round 3, A3-19): GC of a window with N bases = G+C over the non-N bases, as Primer3
        // libprimer3.cc gc_and_n_content (100·num_gc/num_gcat, N excluded) and the canonical CalculateGcFractionFast.
        // primer3-py 2.3.1 design_primers(pick_hyb_probe_only, PRIMER_INTERNAL_MAX_NS_ACCEPTED = 2) on this 20-mer:
        // PRIMER_INTERNAL_0_GC_PERCENT = 50.0 (9 G/C of 18 non-N bases; the former G+C/length gave 0.45).
        const string target = "GACNTGAAGCNCTTAGCAAC";
        var param = new ProbeDesigner.ProbeParameters(
            MinLength: 20, MaxLength: 20, MinTm: -1000, MaxTm: 1000, MinGc: 0.0, MaxGc: 1.0,
            MaxHomopolymer: 10, AvoidSecondaryStructure: false, MaxSelfComplementarity: 1.0);

        var probes = ProbeDesigner.DesignProbes(target, param).ToList();
        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 20, overlap: 0, parameters: param);

        Assert.That(probes, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(probes[0].GcContent, Is.EqualTo(0.5));
            Assert.That(tiling.Probes.Single().GcContent, Is.EqualTo(0.5), "tiling path (eager) agrees");
        });
    }

    #endregion

    #region DesignTilingProbes (Must)

    [Test]
    public void DesignTilingProbes_CoversExpectedPositions()
    {
        // M9: Tiling probes cover expected positions
        // 208-char sequence, probeLength=50, overlap=10 → step=40
        // Grid starts 0, 40, 80, 120 (160 > 208-50=158) end at 169; (158 mod 40 = 38 ≠ 0) → end-anchored window at
        // 158 (CATCH make_candidate_probes_from_sequence) covers 158..207. Coverage: positions 0-207 = 208.
        // (B07 audit round 7, A7-2: previously 4 probes / coverage 170 — the 3' tail 170..207 was never tiled.)
        string target = new string('A', 100) + "GCGCGCGC" + new string('T', 100);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 50, overlap: 10);

        Assert.Multiple(() =>
        {
            Assert.That(tiling.Probes.Count, Is.EqualTo(5), "Expected 4 grid probes + 1 end-anchored probe");
            Assert.That(tiling.Coverage, Is.EqualTo(208), "Expected coverage of all 208 positions");

            var starts = tiling.Probes.Select(p => p.Start).ToList();
            Assert.That(starts, Is.EqualTo(new[] { 0, 40, 80, 120, 158 }),
                "Tiling probes should start at exact positions");
            Assert.That(tiling.Probes[^1].Sequence, Is.EqualTo(target.Substring(158, 50)));
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
        // Grid probes at positions 0, 30, 60, 90 (120 > 110) + end-anchored probe at 110 (110 mod 30 ≠ 0; A7-2)
        string target = new string('G', 50) + new string('C', 50) + new string('A', 50);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 40, overlap: 10);

        Assert.That(tiling.Probes.Count, Is.EqualTo(5), "Expected 5 tiling probes");

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

    // ── B07 audit round 7, A7-2: argument guards, end-anchored window, truthful coverage ─────────────────────────

    [Test]
    public void DesignTilingProbes_TargetShorterThanProbe_ThrowsArgumentOutOfRange()
    {
        // Repro 1: 30-nt target with the default 60-nt probe used to throw InvalidOperationException
        // ("Sequence contains no elements", Average over no windows). CATCH rejects such a sequence (ValueError).
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignTilingProbes(new string('A', 30)));
        Assert.That(ex!.ParamName, Is.EqualTo("probeLength"));
        Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignTilingProbes("", probeLength: 20, overlap: 0));
        Assert.Throws<ArgumentNullException>(() => ProbeDesigner.DesignTilingProbes(null!, probeLength: 20, overlap: 0));
    }

    [TestCase(0, 0, "probeLength")]
    [TestCase(-5, -10, "probeLength")]
    [TestCase(20, 20, "overlap")]
    [TestCase(20, 25, "overlap")]
    public void DesignTilingProbes_NonAdvancingStep_ThrowsInsteadOfLooping(int probeLength, int overlap, string param)
    {
        // Repro 2: probeLength ≤ 0 or overlap ≥ probeLength (step ≤ 0) used to loop forever.
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() =>
            ProbeDesigner.DesignTilingProbes(new string('A', 100), probeLength, overlap));
        Assert.That(ex!.ParamName, Is.EqualTo(param));
    }

    [Test]
    public void DesignTilingProbes_110nt_60_20_TilesThe3PrimeTail()
    {
        // Repro 3: 110 nt, probe 60, overlap 20 → grid starts 0, 40 (80 > 50) cover 0..99 only; the last 10 nt were
        // never covered yet the set claimed full coverage. CATCH anchors a final probe at 110 − 60 = 50.
        string target = string.Concat(Enumerable.Repeat("ACGTTGCAAG", 10)) + "GATCCGATCA"; // 110 nt
        Assert.That(target, Has.Length.EqualTo(110));
        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 60, overlap: 20);

        Assert.Multiple(() =>
        {
            Assert.That(tiling.Probes.Select(p => p.Start), Is.EqualTo(new[] { 0, 40, 50 }));
            Assert.That(tiling.Probes[^1].End, Is.EqualTo(109));
            Assert.That(tiling.Probes[^1].Sequence, Is.EqualTo(target[50..]));
            Assert.That(tiling.Coverage, Is.EqualTo(IndependentCoverage(110, new[] { 0, 40, 50 }, 60)).And.EqualTo(110));
        });
    }

    // Starts = the grid 0, step, … ≤ L − P, plus L − P iff (L − P) mod step ≠ 0; Coverage = independently counted
    // union of [start, start + P). Includes CATCH's own corner cases: 100/60/40 (CATCH duplicates start 40) and
    // 120/60/40 (CATCH's len mod stride = 0 test misses the 100..119 tail), P = L, step 1, and negative overlap (gaps).
    [TestCase(100, 60, 20, new[] { 0, 40 })]
    [TestCase(120, 60, 20, new[] { 0, 40, 60 })]
    [TestCase(60, 60, 20, new[] { 0 })]
    [TestCase(61, 60, 20, new[] { 0, 1 })]
    [TestCase(25, 10, 9, new[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 })]
    [TestCase(22, 10, 4, new[] { 0, 6, 12 })]
    [TestCase(30, 10, 4, new[] { 0, 6, 12, 18, 20 })]
    [TestCase(45, 10, -5, new[] { 0, 15, 30, 35 })]
    [TestCase(55, 10, -5, new[] { 0, 15, 30, 45 })]
    public void DesignTilingProbes_StartsAndCoverage_MatchIndependentTiling(int length, int probeLength, int overlap, int[] expectedStarts)
    {
        string target = string.Concat(Enumerable.Repeat("GATTACAGCG", (length + 9) / 10))[..length];
        var param = new ProbeDesigner.ProbeParameters(
            MinLength: probeLength, MaxLength: probeLength, MinTm: -1000, MaxTm: 1000, MinGc: 0.0, MaxGc: 1.0,
            MaxHomopolymer: 10, AvoidSecondaryStructure: false, MaxSelfComplementarity: 1.0);

        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength, overlap, param);

        Assert.Multiple(() =>
        {
            Assert.That(tiling.Probes.Select(p => p.Start), Is.EqualTo(expectedStarts));
            Assert.That(tiling.Probes.Select(p => p.Sequence), Is.EqualTo(expectedStarts.Select(s => target.Substring(s, probeLength))));
            Assert.That(tiling.Probes.Select(p => p.End), Is.EqualTo(expectedStarts.Select(s => s + probeLength - 1)));
            Assert.That(tiling.Coverage, Is.EqualTo(IndependentCoverage(length, expectedStarts, probeLength)));
            if (overlap >= 0)
                Assert.That(tiling.Coverage, Is.EqualTo(length), "non-negative overlap covers every base");
        });
    }

    private static int IndependentCoverage(int length, int[] starts, int probeLength)
    {
        var covered = new bool[length];
        foreach (int s in starts)
            for (int i = s; i < s + probeLength; i++)
                covered[i] = true;
        return covered.Count(c => c);
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

    #region B07 audit round 7, A7-3 + A7-4: degenerate / null arguments

    private const string BeaconTarget = "ATGCGTACGTTAGCCGATCGATCGGCTAGCTAGGATCCGATCGTAGCTAGCATCGACTGAC";

    // A7-3: probeLength 0 returned a stem-only beacon (empty loop), stemLength 0 a loop with no arms, negatives threw a
    // raw string AOORE ('length' / 'count'). A beacon needs a loop and two arms (Tyagi & Kramer 1996); the 15–30 nt /
    // 5–7 bp ranges are documented recommendations, so only impossible values (≤ 0) are rejected.
    [TestCase(0, 5, "probeLength")]
    [TestCase(-1, 5, "probeLength")]
    [TestCase(int.MinValue, 5, "probeLength")]
    [TestCase(25, 0, "stemLength")]
    [TestCase(25, -1, "stemLength")]
    [TestCase(25, int.MinValue, "stemLength")]
    public void DesignMolecularBeacon_NonPositiveLoopOrStem_ThrowsArgumentOutOfRange(int probeLength, int stemLength, string param)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => ProbeDesigner.DesignMolecularBeacon(BeaconTarget, probeLength, stemLength));
        Assert.That(ex!.ParamName, Is.EqualTo(param));
    }

    [TestCase(1, 1)]
    [TestCase(10, 3)]
    [TestCase(35, 9)]
    public void DesignMolecularBeacon_PositiveOutsideRecommendedRanges_StillDesigned(int probeLength, int stemLength)
    {
        // Recommendations are not requirements: the smallest legal beacon still has a loop of probeLength and two arms.
        var beacon = ProbeDesigner.DesignMolecularBeacon(BeaconTarget, probeLength, stemLength);
        Assert.That(beacon, Is.Not.Null);
        Assert.That(beacon!.Value.Sequence, Has.Length.EqualTo(probeLength + 2 * stemLength));
    }

    [Test]
    public void DesignMolecularBeacon_Null_ThrowsArgumentNull()
    {
        // A7-4: was NullReferenceException.
        var ex = Assert.Throws<ArgumentNullException>(() => ProbeDesigner.DesignMolecularBeacon(null!));
        Assert.That(ex!.ParamName, Is.EqualTo("targetSequence"));
    }

    [Test]
    public void AnalyzeOligo_Null_ThrowsArgumentNull()
    {
        // A7-4: was NullReferenceException (ProbeDesigner's design/evaluate siblings throw ArgumentNullException).
        var ex = Assert.Throws<ArgumentNullException>(() => ProbeDesigner.AnalyzeOligo(null!));
        Assert.That(ex!.ParamName, Is.EqualTo("sequence"));
    }

    [TestCase(null)]
    [TestCase("")]
    public void CalculateExtinctionCoefficient_NullOrEmpty_ReturnsZero(string? sequence)
    {
        // A7-4: null was NullReferenceException; the siblings CalculateExtinctionCoefficientNearestNeighbor and
        // CalculateMolecularWeight return 0 for null / empty, and so does the mononucleotide sum of no bases.
        Assert.That(ProbeDesigner.CalculateExtinctionCoefficient(sequence!), Is.EqualTo(0));
        Assert.That(ProbeDesigner.CalculateExtinctionCoefficientNearestNeighbor(sequence!), Is.EqualTo(0));
    }

    [Test]
    public void CalculateExtinctionCoefficient_NonDegenerate_Unchanged()
    {
        // Mononucleotide sum: A 15400 + C 7400 + G 11500 + T 8700 + U 9900 + N 10000.
        Assert.That(ProbeDesigner.CalculateExtinctionCoefficient("acgtUN"), Is.EqualTo(62900));
    }

    [Test]
    public void DesignTilingProbes_HugeNegativeOverlap_NoStepOverflow()
    {
        // Sweep: overlap = int.MinValue overflowed step = probeLength − overlap to a negative value and crashed in
        // Substring. A negative overlap is gapped tiling; with a step beyond the target only window 0 and the
        // end-anchored window remain.
        string target = new string('A', 40) + new string('G', 40);
        var tiling = ProbeDesigner.DesignTilingProbes(target, probeLength: 20, overlap: int.MinValue);
        Assert.That(tiling.Probes.Select(p => p.Start), Is.EqualTo(new[] { 0, 60 }));
        Assert.That(tiling.Coverage, Is.EqualTo(40));
    }

    #endregion

    #region B07 heavy tier (HEAVY-2, F75): ProbeParameters / condition guards found by the guard fuzz

    private const string GuardTarget = "ATGCGTACGTTAGCCGATCGATCGGCTAGCTAGGATCCGATCGTAGCTAGCATCGACTGACTTAGGCA";

    private static IEnumerable<TestCaseData> IllegalDesignParameters()
    {
        var q = ProbeDesigner.Defaults.qPCR;
        // MinLength < 1 crashed the window scan with IndexOutOfRangeException (Primer3: PRIMER_INTERNAL_MIN_SIZE ≥ 1).
        yield return new TestCaseData(q with { MinLength = 0, MaxLength = 3 }).SetName("MinLength 0");
        yield return new TestCaseData(q with { MinLength = -5, MaxLength = 3 }).SetName("MinLength -5");
        // Illegal conditions surfaced lazily from the Tm with a foreign ParamName (maxNearestNeighborLength / naConcentration).
        yield return new TestCaseData(q with { MaxNearestNeighborLength = -1 }).SetName("MaxNearestNeighborLength -1");
        yield return new TestCaseData(q with { MonovalentMillimolar = -1 }).SetName("monovalent -1");
        yield return new TestCaseData(q with { DnaConcentrationNanomolar = 0 }).SetName("DNA 0");
        yield return new TestCaseData(q with { DivalentMillimolar = double.NaN }).SetName("divalent NaN");
        yield return new TestCaseData(q with { DntpMillimolar = double.PositiveInfinity }).SetName("dNTP +inf");
    }

    [TestCaseSource(nameof(IllegalDesignParameters))]
    public void DesignProbes_IllegalParameters_ThrowEagerlyWithParameterName(ProbeDesigner.ProbeParameters p)
    {
        // Eager: the exception is raised by the call itself, before enumeration.
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbes(GuardTarget, p))!.ParamName,
                Is.EqualTo("parameters"));
            var index = global::SuffixTree.SuffixTree.Build(GuardTarget);
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignProbes(GuardTarget, index, p))!.ParamName,
                Is.EqualTo("parameters"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(() => ProbeDesigner.DesignAntisenseProbes(GuardTarget, p).ToList())!.ParamName,
                Is.EqualTo("parameters"));
        });
    }

    [Test]
    public void DesignTilingProbes_IllegalConditions_Throw_UnusedLengthsIgnored()
    {
        var q = ProbeDesigner.Defaults.qPCR;
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(
                () => ProbeDesigner.DesignTilingProbes(GuardTarget, 20, 5, q with { MonovalentMillimolar = 0 }))!.ParamName,
                Is.EqualTo("parameters"));
            Assert.That(Assert.Throws<ArgumentOutOfRangeException>(
                () => ProbeDesigner.DesignTilingProbes(GuardTarget, 20, 5, q with { MaxNearestNeighborLength = -1 }))!.ParamName,
                Is.EqualTo("parameters"));
            // Tiling windows have probeLength nt; the parameters' MinLength / MaxLength are not used, so 0 is accepted.
            var set = ProbeDesigner.DesignTilingProbes(GuardTarget, 20, 5, q with { MinLength = 0, MaxLength = 0 });
            Assert.That(set.Coverage, Is.EqualTo(GuardTarget.Length));
        });
    }

    [Test]
    public void DesignProbes_LegalPresets_Unchanged()
    {
        // The new guards accept every preset (OligoArray 1 M / 1 µM included) and the defaults.
        foreach (var p in new[] { ProbeDesigner.Defaults.qPCR, ProbeDesigner.Defaults.Microarray with { StructureScreen = ProbeDesigner.ProbeStructureScreen.Heuristic } })
            Assert.DoesNotThrow(() => ProbeDesigner.DesignProbes(GuardTarget, p, 3).ToList());
        Assert.DoesNotThrow(() => ProbeDesigner.DesignProbes(GuardTarget, (ProbeDesigner.ProbeParameters?)null, 1).ToList());
    }

    [TestCase(0.0, 50.0, 0.0, 0.0, "dnaConcentrationNanomolar")]
    [TestCase(double.PositiveInfinity, 50.0, 0.0, 0.0, "dnaConcentrationNanomolar")]
    [TestCase(50.0, 0.0, 0.0, 0.0, "monovalentMillimolar")]
    [TestCase(50.0, double.PositiveInfinity, 0.0, 0.0, "monovalentMillimolar")]
    [TestCase(50.0, 50.0, -1.0, 0.0, "divalentMillimolar")]
    [TestCase(50.0, 50.0, double.NaN, 0.0, "divalentMillimolar")]
    [TestCase(50.0, 50.0, 0.0, double.PositiveInfinity, "dntpMillimolar")]
    public void EvaluateTaqManProbe_IllegalConditions_ThrowWithOwnParameterName(double dna, double mv, double dv, double dntp, string param)
    {
        // Before: an infinite salt reached ThermoConstants and threw with ParamName "naConcentration".
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => ProbeDesigner.EvaluateTaqManProbe("CCATCACCCTACATCACC", null, 18, 22, dna, mv, dv, dntp));
        Assert.That(ex!.ParamName, Is.EqualTo(param));
    }

    [Test]
    public void DesignMolecularBeacon_StemTooLongForAString_ThrowsArgumentOutOfRange()
    {
        // Before: stemLength = int.MaxValue → OutOfMemoryException building the arms.
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => ProbeDesigner.DesignMolecularBeacon(BeaconTarget, 25, int.MaxValue));
        Assert.That(ex!.ParamName, Is.EqualTo("stemLength"));
    }

    #endregion
}
