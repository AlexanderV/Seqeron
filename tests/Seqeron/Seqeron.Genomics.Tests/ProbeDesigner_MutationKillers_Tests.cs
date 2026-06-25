using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Targeted mutation-killing tests for ProbeDesigner.cs (checklist 04 rows 24/25:
/// PROBE-DESIGN-001, PROBE-VALID-001).
///
/// The canonical suite asserts invariant RANGES (score∈[0,1], GC∈[0,1]); it leaves the
/// exact formulas un-pinned, so Stryker arithmetic/relational/logical mutants survived
/// (37.6% baseline). These tests transcribe each documented formula/rule as ground truth
/// and assert the production value EQUALS it, so a mutated operator diverges and is killed.
/// </summary>
[TestFixture]
public class ProbeDesigner_MutationKillers_Tests
{
    // ── Oligo analysis: molecular weight (Σ residue mass − (n−1)·H2O) ────────────────

    [Test]
    [TestCase("ACGT", 1253.8)]   // 331.2+307.2+347.2+322.2 − 3·18
    [TestCase("A", 331.2)]       // single base, no phosphodiester subtraction
    [TestCase("AA", 644.4)]      // 662.4 − 1·18
    [TestCase("U", 308.2)]       // uracil branch
    public void CalculateMolecularWeight_MatchesResidueSumMinusWater(string seq, double expected)
    {
        ProbeDesigner.CalculateMolecularWeight(seq).Should().BeApproximately(expected, 1e-6);
    }

    [Test]
    public void CalculateMolecularWeight_WaterSubtractionScalesWithLength()
    {
        // (n−1)·18 term: AAA = 3·331.2 − 2·18 = 993.6 − 36 = 957.6 (kills the (n−1)·18 mutants)
        ProbeDesigner.CalculateMolecularWeight("AAA").Should().BeApproximately(957.6, 1e-6);
    }

    // ── Oligo analysis: extinction coefficient (Σ per-base ε260) ─────────────────────

    [Test]
    [TestCase("ACGT", 43000.0)] // 15400+7400+11500+8700
    [TestCase("AA", 30800.0)]   // 2·15400
    [TestCase("U", 9900.0)]     // uracil branch
    public void CalculateExtinctionCoefficient_MatchesPerBaseSum(string seq, double expected)
    {
        ProbeDesigner.CalculateExtinctionCoefficient(seq).Should().BeApproximately(expected, 1e-6);
    }

    // ── Oligo analysis: concentration (Beer–Lambert c = A/(ε·l)·1e6) ─────────────────

    [Test]
    [TestCase(1.0, 10000.0, 1.0, 100.0)]  // 1/10000·1e6
    [TestCase(3.0, 10000.0, 1.0, 300.0)]  // absorbance scales numerator
    [TestCase(1.0, 10000.0, 2.0, 50.0)]   // path length doubles denominator
    [TestCase(1.0, 20000.0, 1.0, 50.0)]   // extinction doubles denominator
    public void CalculateConcentration_MatchesBeerLambert(double a260, double eps, double path, double expected)
    {
        ProbeDesigner.CalculateConcentration(a260, eps, path).Should().BeApproximately(expected, 1e-6);
    }

    // ── AnalyzeOligo: Tm (Wallace for <14 nt, salt-adjusted for ≥14 nt) + GC fraction ─

    [Test]
    [TestCase("ACGT", 2, 2)] // at=2, gc=2
    [TestCase("AAAA", 4, 0)]
    [TestCase("GCGC", 0, 4)]
    public void AnalyzeOligo_ShortOligoTm_UsesWallaceRuleOnExactAtGcCounts(string seq, int at, int gc)
    {
        var result = ProbeDesigner.AnalyzeOligo(seq);

        // Pins the A/T and G/C counting predicates: a mis-count changes the Wallace Tm.
        result.Tm.Should().BeApproximately(ThermoConstants.CalculateWallaceTm(at, gc), 1e-6);
    }

    [Test]
    public void AnalyzeOligo_LongOligoTm_UsesSaltAdjustedFormula()
    {
        // 20 nt (≥ WallaceMaxLength) → salt-adjusted branch on the GC fraction.
        const string seq = "GCGCGCATATGCGCGCATAT"; // 20 nt, gc = 12/20 = 0.6
        var result = ProbeDesigner.AnalyzeOligo(seq);

        result.Tm.Should().BeApproximately(ThermoConstants.CalculateSaltAdjustedTm(0.6, seq.Length), 1e-6);
    }

    [Test]
    [TestCase("ACGT", 0.5)]
    [TestCase("GGGG", 1.0)]
    [TestCase("AAAA", 0.0)]
    [TestCase("GGGGCCCCAA", 0.8)]
    public void AnalyzeOligo_GcContent_MatchesFraction(string seq, double expected)
    {
        ProbeDesigner.AnalyzeOligo(seq).GcContent.Should().BeApproximately(expected, 1e-9);
    }

    // ── ValidateProbe: specificity score (0 hits→0, 1→1, N→1/N) + off-target rule ────

    private const string Probe8 = "ACGTTGCA"; // not self-complementary, no 2° structure at 8 nt

    [Test]
    public void ValidateProbe_NoMatch_SpecificityZero()
    {
        var v = ProbeDesigner.ValidateProbe(Probe8, new[] { "TTTTTTTTTTTT" }, maxMismatches: 0);

        v.OffTargetHits.Should().Be(0);
        v.SpecificityScore.Should().Be(0.0, "0 hits ⇒ probe does not hybridize ⇒ specificity 0");
    }

    [Test]
    public void ValidateProbe_UniqueMatch_SpecificityOne()
    {
        var v = ProbeDesigner.ValidateProbe(Probe8, new[] { "GGGG" + Probe8 + "GGGG" }, maxMismatches: 0);

        v.OffTargetHits.Should().Be(1);
        v.SpecificityScore.Should().Be(1.0);
        v.Issues.Should().NotContain(i => i.Contains("off-target"));
    }

    [Test]
    public void ValidateProbe_TwoMatches_SpecificityHalf_AndOffTargetIssue()
    {
        var refs = new[] { "GGGG" + Probe8 + "GGGG", "AAAA" + Probe8 + "AAAA" };
        var v = ProbeDesigner.ValidateProbe(Probe8, refs, maxMismatches: 0);

        v.OffTargetHits.Should().Be(2);
        v.SpecificityScore.Should().BeApproximately(0.5, 1e-9, "N hits ⇒ 1/N");
        v.Issues.Should().Contain(i => i.Contains("off-target"), "more than one hit is flagged");
    }

    [Test]
    public void ValidateProbe_EmptyProbe_IsInvalidDegenerate()
    {
        var v = ProbeDesigner.ValidateProbe("", new[] { "ACGT" });

        v.IsValid.Should().BeFalse();
        v.SpecificityScore.Should().Be(0.0);
        v.Issues.Should().Contain(i => i.Contains("Empty"));
    }

    // ── Self-complementarity = (# positions where base == its reverse-complement)/len ─

    [Test]
    [TestCase("ACGT", 1.0)]  // revcomp == ACGT → all 4 positions match
    [TestCase("AAAA", 0.0)]  // revcomp == TTTT → 0 matches
    [TestCase("ACGA", 0.5)]  // revcomp == TCGT → C,G match → 2/4
    public void ValidateProbe_SelfComplementarity_MatchesPositionwiseFraction(string seq, double expected)
    {
        var v = ProbeDesigner.ValidateProbe(seq, System.Array.Empty<string>(), maxMismatches: 0);

        v.SelfComplementarity.Should().BeApproximately(expected, 1e-9);
    }

    // ── DesignTilingProbes: step = probeLength − overlap, full coverage, Tm stats ─────

    [Test]
    public void DesignTilingProbes_ProbeStartsAndCoverage_AreExact()
    {
        const string target = "GCGCGCGCGCATATATATATGC"; // 22 nt
        var set = ProbeDesigner.DesignTilingProbes(target, probeLength: 10, overlap: 4);

        // step = 10 − 4 = 6 → starts 0,6,12 (start ≤ 22 − 10 = 12)
        set.Probes.Select(p => p.Start).Should().Equal(0, 6, 12);
        set.Probes.Should().OnlyContain(p => p.Type == ProbeDesigner.ProbeType.Tiling);
        // windows [0..9],[6..15],[12..21] union = positions 0..21 = 22 covered
        set.Coverage.Should().Be(22);
        // End coordinate = start + probeLength − 1
        set.Probes[0].End.Should().Be(9);
        set.Probes[2].End.Should().Be(21);
    }

    [Test]
    public void DesignTilingProbes_TmStatistics_MatchProbeAggregates()
    {
        const string target = "GCGCGCGCGCATATATATATGCATGCATGC"; // 30 nt
        var set = ProbeDesigner.DesignTilingProbes(target, probeLength: 10, overlap: 4);

        double expectedMean = set.Probes.Average(p => p.Tm);
        double expectedRange = set.Probes.Max(p => p.Tm) - set.Probes.Min(p => p.Tm);

        set.MeanTm.Should().BeApproximately(expectedMean, 1e-9);
        set.TmRange.Should().BeApproximately(expectedRange, 1e-9, "range = max − min");
    }

    // ── DesignMolecularBeacon: stem5 = G^(stem/2)·C^(stem−stem/2), stem3 = revcomp(stem5) ─

    [Test]
    public void DesignMolecularBeacon_StemConstruction_IsExact()
    {
        // Uniform target → every loop window scores equally → first (start 0) wins (strict >).
        string target = string.Concat(Enumerable.Repeat("ACGT", 10)); // 40 nt
        var beacon = ProbeDesigner.DesignMolecularBeacon(target, probeLength: 25, stemLength: 5);

        beacon.Should().NotBeNull();
        // stemLength 5 → stem5 = "GG"+"CCC" = "GGCCC"; stem3 = revcomp("GGCCC") = "GGGCC"
        beacon!.Value.Sequence.Should().StartWith("GGCCC").And.EndWith("GGGCC");
        beacon.Value.Type.Should().Be(ProbeDesigner.ProbeType.MolecularBeacon);
        beacon.Value.Start.Should().Be(0, "uniform target makes the first window the best by strict >");
        // Beacon body length = stem5 + loop(25) + stem3 = 5 + 25 + 5 = 35
        beacon.Value.Sequence.Length.Should().Be(35);
    }

    [Test]
    public void DesignMolecularBeacon_TargetShorterThanProbe_ReturnsNull()
    {
        ProbeDesigner.DesignMolecularBeacon("ACGT", probeLength: 25).Should().BeNull();
    }

    // ── EvaluateProbeWithGc scoring: exact penalty bookkeeping via DesignProbes ───────

    private static ProbeDesigner.ProbeParameters Tiny(
        double minGc = 0.0, double maxGc = 1.0, double minTm = 0, double maxTm = 1000,
        int maxHomo = 100, bool avoidStructure = false, double maxSelfComp = 1.0) =>
        new(MinLength: 4, MaxLength: 4, MinTm: minTm, MaxTm: maxTm,
            MinGc: minGc, MaxGc: maxGc, MaxHomopolymer: maxHomo,
            AvoidSecondaryStructure: avoidStructure, MaxSelfComplementarity: maxSelfComp);

    [Test]
    public void EvaluateProbe_NoPenalties_ScoreIsOne()
    {
        // "ACGT": GC 0.5 in range, Tm in range, no homopolymer issue, starts A / ends T → no 3'/5' GC penalty.
        var probes = ProbeDesigner.DesignProbes("ACGT", Tiny()).ToList();

        probes.Should().ContainSingle();
        probes[0].Score.Should().BeApproximately(1.0, 1e-9);
    }

    [Test]
    public void EvaluateProbe_GcEndsPenalty_SubtractsPointZeroTwoEach()
    {
        // "GACG": starts G (+0.02) and ends G (+0.02) → score 1.0 − 0.04 = 0.96.
        var probes = ProbeDesigner.DesignProbes("GACG", Tiny()).ToList();

        probes.Should().ContainSingle();
        probes[0].Score.Should().BeApproximately(0.96, 1e-9, "G/C at both termini each cost 0.02");
    }

    [Test]
    public void EvaluateProbe_GcOutOfRange_SubtractsPointThree()
    {
        // GC=0.5 with MinGc 0.55 (within the 0.1 early-reject slack, so still evaluated) → −0.3 penalty.
        // "ACGT" starts A / ends T → no terminus penalty, so score = 1.0 − 0.3 = 0.7.
        var probes = ProbeDesigner.DesignProbes("ACGT", Tiny(minGc: 0.55)).ToList();

        probes.Should().ContainSingle();
        probes[0].Score.Should().BeApproximately(0.7, 1e-9);
        probes[0].Warnings.Should().Contain(w => w.Contains("GC content"));
    }

    [Test]
    public void EvaluateProbe_TmOutOfRange_SubtractsPointThree()
    {
        // Force Tm below MinTm (MinTm = 999) → −0.3. "ACGT": no GC/terminus penalty → 0.7.
        var probes = ProbeDesigner.DesignProbes("ACGT", Tiny(minTm: 999)).ToList();

        probes.Should().ContainSingle();
        probes[0].Score.Should().BeApproximately(0.7, 1e-9);
        probes[0].Warnings.Should().Contain(w => w.Contains("Tm"));
    }

    [Test]
    public void DesignProbes_TargetShorterThanMinLength_ReturnsEmpty()
    {
        ProbeDesigner.DesignProbes("ACG", Tiny()).Should().BeEmpty();
    }

    [Test]
    public void DesignProbes_GcBeyondEarlyRejectSlack_RejectsCandidate()
    {
        // GC=0.5, MinGc 0.65 ⇒ 0.5 < 0.65 − 0.1 = 0.55 ⇒ early-rejected ⇒ no probe emitted.
        ProbeDesigner.DesignProbes("ACGT", Tiny(minGc: 0.65)).Should().BeEmpty();
    }
}
