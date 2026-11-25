// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the MolTools area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("MolTools")]
public class MolToolsCombinatorialTests
{

    /// <summary>Permissive structural filters so only the length/GC/Tm windows can gate acceptance.</summary>
    private static PrimerParameters WithWindows(int minLen, int maxLen, double gcLo, double gcHi, double tmLo, double tmHi) =>
        ThermodynamicPrimerDesigner.DefaultParameters with
        {
            MinLength = minLen,
            MaxLength = maxLen,
            MinGcContent = gcLo,
            MaxGcContent = gcHi,
            MinTm = tmLo,
            MaxTm = tmHi,
            MaxHomopolymer = 1000,
            MaxDinucleotideRepeats = 1000,
            Check3PrimeStability = false,
            Avoid3PrimeGC = false,
        };

    [Test, Combinatorial]
    public void PrimerDesign_AcceptanceIsConjunctionOfLengthGcTmWindows(
        [Values(15, 20, 22)] int minLen,
        [Values(18, 20, 25)] int maxLen,
        [Values(0, 1, 2)] int gcWin,
        [Values(0, 1, 2)] int tmWin)
    {
        var (gcLo, gcHi) = GcWindows[gcWin];
        var (tmLo, tmHi) = TmWindows[tmWin];

        var cand = ThermodynamicPrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true,
            WithWindows(minLen, maxLen, gcLo, gcHi, tmLo, tmHi));

        double gc = PrimerDesigner.CalculateGcContent(CleanPrimer20);
        double tm = ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(CleanPrimer20);
        int len = CleanPrimer20.Length;

        bool lenOk = len >= minLen && len <= maxLen;
        bool gcOk = gc >= gcLo && gc <= gcHi;
        bool tmOk = tm >= tmLo && tm <= tmHi;

        cand.IsValid.Should().Be(lenOk && gcOk && tmOk,
            $"accept ⟺ len∈[{minLen},{maxLen}] ∧ GC∈[{gcLo},{gcHi}] ∧ Tm∈[{tmLo},{tmHi}]");

        // Every violated window contributes exactly its own diagnostic; satisfied windows stay silent.
        cand.Issues.Any(i => i.StartsWith("Length")).Should().Be(!lenOk);
        cand.Issues.Any(i => i.StartsWith("GC content")).Should().Be(!gcOk);
        cand.Issues.Any(i => i.StartsWith("Tm")).Should().Be(!tmOk);
    }

    /// <summary>
    /// Interaction witness: from a config that accepts the primer, narrowing ANY
    /// single window past the primer's measured value flips acceptance to reject —
    /// confirming each axis is an independent necessary condition.
    /// </summary>
    [Test]
    public void PrimerDesign_EachWindow_IndependentlyGatesAcceptance()
    {
        var wide = WithWindows(1, 100, 0, 100, 0, 200);
        ThermodynamicPrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true, wide).IsValid
            .Should().BeTrue("a clean primer is accepted under an all-permissive config");

        double gc = PrimerDesigner.CalculateGcContent(CleanPrimer20);
        double tm = ThermodynamicPrimerDesigner.CalculateMeltingTemperaturePrimer3(CleanPrimer20);

        ThermodynamicPrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true, wide with { MinGcContent = gc + 5 })
            .IsValid.Should().BeFalse("a GC floor above the primer's GC excludes it");
        ThermodynamicPrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true, wide with { MaxTm = tm - 1 })
            .IsValid.Should().BeFalse("a Tm ceiling below the primer's Tm excludes it");
        ThermodynamicPrimerDesigner.EvaluatePrimer(CleanPrimer20, 0, true, wide with { MaxLength = CleanPrimer20.Length - 1 })
            .IsValid.Should().BeFalse("a length ceiling below the primer's length excludes it");
    }

    /// <summary>
    /// Worked example through the end-to-end design pipeline: any pair
    /// <see cref="ThermodynamicPrimerDesigner.DesignPrimers"/> returns must honour every
    /// configured window on BOTH primers (acceptance is enforced during selection).
    /// </summary>
    [Test]
    public void PrimerDesign_DesignedPair_HonoursEveryWindow()
    {
        var template = new DnaSequence(DiverseDna(600));
        // Wide Tm window on the Primer3-default (SantaLucia NN) scale for 18–25-mers.
        var param = ThermodynamicPrimerDesigner.DefaultParameters with
        {
            MinLength = 18, MaxLength = 25, MinGcContent = 35, MaxGcContent = 65, MinTm = 45, MaxTm = 65,
        };

        var result = ThermodynamicPrimerDesigner.DesignPrimers(template, 280, 320, param);

        foreach (var p in new[] { result.Forward, result.Reverse }.Where(p => p is not null))
        {
            p!.Length.Should().BeInRange(param.MinLength, param.MaxLength);
            p.GcContent.Should().BeInRange(param.MinGcContent, param.MaxGcContent);
            p.MeltingTemperature.Should().BeInRange(param.MinTm, param.MaxTm);
        }
    }

    private static ThermodynamicProbeDesigner.ProbeParameters MkProbeParam(
        int minLen, int maxLen, double gcLo, double gcHi, double tmLo, double tmHi) =>
        new(MinLength: minLen, MaxLength: maxLen, MinTm: tmLo, MaxTm: tmHi,
            MinGc: gcLo, MaxGc: gcHi, MaxHomopolymer: 100,
            AvoidSecondaryStructure: false, MaxSelfComplementarity: 1.0);

    [Test, Combinatorial]
    public void ProbeDesign_ReturnedProbes_HonourLengthAndGcWindows(
        [Values(18, 20, 22)] int minLen,
        [Values(24, 28, 32)] int maxLen,
        [Values(0, 1, 2)] int gcWin,
        [Values(0, 1, 2)] int tmWin)
    {
        var (gcLo, gcHi) = ProbeGcWindows[gcWin];
        var (tmLo, tmHi) = ProbeTmWindows[tmWin];

        var probes = ThermodynamicProbeDesigner.DesignProbes(
            ProbeTarget, MkProbeParam(minLen, maxLen, gcLo, gcHi, tmLo, tmHi), maxProbes: 20).ToList();

        foreach (var p in probes)
        {
            p.Sequence.Length.Should().BeInRange(minLen, maxLen, "probe length is hard-bounded by the window");
            p.GcContent.Should().BeInRange(gcLo - 0.1 - 1e-9, gcHi + 0.1 + 1e-9,
                "GC honours the window within the ±0.1 early-rejection slack");
        }
    }

    /// <summary>
    /// Interaction witness: the Tm window is the soft scoring axis. With every other
    /// term held neutral, moving the Tm band from "includes all probe Tms" to "excludes
    /// them" lowers each common probe's score by exactly the 0.3 Tm penalty.
    /// </summary>
    [Test]
    public void ProbeDesign_TmWindow_ContributesExactScoringPenalty()
    {
        var tmIn = MkProbeParam(20, 20, 0.0, 1.0, 0, 200);    // GC + Tm both always in range
        var tmOut = tmIn with { MinTm = 200, MaxTm = 300 };    // only the Tm term now fires

        var inByStart = ThermodynamicProbeDesigner.DesignProbes(ProbeTarget, tmIn, maxProbes: 50).ToDictionary(p => p.Start);
        var outByStart = ThermodynamicProbeDesigner.DesignProbes(ProbeTarget, tmOut, maxProbes: 50).ToDictionary(p => p.Start);

        var common = inByStart.Keys.Intersect(outByStart.Keys).ToList();
        common.Should().NotBeEmpty("a uniform Tm penalty preserves ranking, so the same probes are returned");
        foreach (int start in common)
            (inByStart[start].Score - outByStart[start].Score).Should().BeApproximately(0.3, 1e-6,
                "an out-of-Tm probe loses exactly the 0.3 Tm penalty");
    }

    /// <summary>
    /// Interaction witness: the GC band carries a ±0.1 early-rejection slack. A 65%-GC
    /// 20-mer is accepted when 0.65 ≤ MaxGc+0.1 (band 0.40–0.60) but rejected once it
    /// exceeds the slack (band 0.40–0.50 ⇒ ceiling 0.60).
    /// </summary>
    [Test]
    public void ProbeDesign_GcWindow_AppliesTenPercentSlack()
    {
        const string gc65 = "GCGCGCGCGCGCGAATATAT"; // 20 nt, 13 GC = 0.65

        ThermodynamicProbeDesigner.DesignProbes(gc65, MkProbeParam(20, 20, 0.40, 0.60, 0, 200), maxProbes: 5)
            .Should().ContainSingle("0.65 is within MaxGc(0.60)+0.1 slack");
        ThermodynamicProbeDesigner.DesignProbes(gc65, MkProbeParam(20, 20, 0.40, 0.50, 0, 200), maxProbes: 5)
            .Should().BeEmpty("0.65 exceeds MaxGc(0.50)+0.1 slack");
    }

    /// <summary>
    /// Worked example through the qPCR-defaults pipeline: every returned probe respects
    /// the strict length window and the result set is capped at maxProbes.
    /// </summary>
    [Test]
    public void ProbeDesign_QpcrDefaults_RespectLengthWindowAndCap()
    {
        var param = ThermodynamicProbeDesigner.Defaults.qPCR; // 20–30 nt
        var probes = ThermodynamicProbeDesigner.DesignProbes(ProbeTarget, param, maxProbes: 10).ToList();

        probes.Count.Should().BeLessThanOrEqualTo(10);
        foreach (var p in probes)
            p.Sequence.Length.Should().BeInRange(param.MinLength, param.MaxLength);
    }

    [Test, Combinatorial]
    public void ProbeValid_SpecificityAndIssues_FollowValidationRules(
        [Values(0, 1, 3)] int offTargetCount,
        [Values(0, 1, 2)] int probeIdx,
        [Values(0.25, 0.40, 0.60)] double selfCompThreshold)
    {
        var (probe, expSelfComp, expSelfAny, expHairpin) = ValidationProbes[probeIdx];
        string reference = BuildOffTargetReference(probe, offTargetCount);

        // Fallback self-structure screen: the self-complementarity criterion is Primer3 alignment-mode self_any /
        // self_end > PRIMER_INTERNAL_MAX_SELF_ANY/_END 12.00 (audit round 2, A6); the fold-back fraction and the legacy
        // selfComplementarityThreshold no longer gate the issue.
        var v = ThermodynamicProbeDesigner.ValidateProbe(probe, new[] { reference }, maxMismatches: 0,
            selfComplementarityThreshold: selfCompThreshold,
            conditions: ThermodynamicProbeDesigner.Defaults.Microarray with { StructureScreen = ProbeDesigner.ProbeStructureScreen.Heuristic });

        // Specificity invariants (#4/#5/#6) — depend solely on off-target multiplicity.
        double expSpec = offTargetCount == 0 ? 0.0 : offTargetCount == 1 ? 1.0 : 1.0 / offTargetCount;
        v.OffTargetHits.Should().Be(offTargetCount);
        v.SpecificityScore.Should().BeApproximately(expSpec, 1e-9);

        // Self-complementarity measured exactly and bounded to [0,1] (#2).
        v.SelfComplementarity.Should().BeApproximately(expSelfComp, 1e-9);
        v.SelfAny.Should().Be(expSelfAny);
        v.SelfComplementarity.Should().BeInRange(0.0, 1.0);
        v.HasSecondaryStructure.Should().Be(expHairpin, "fallback hairpin = PrimerDesigner.HasHairpinPotential");

        // IsValid = no recorded issue (off-target multiplicity, self-complementarity). The fallback stem-loop flag
        // (library convention) is a warning, not an issue (audit round 5, A5-5).
        bool offIssue = offTargetCount > 1;
        bool selfIssue = expSelfAny > PrimerDesigner.Primer3InternalMaxSelfComplementarity;
        int issueCount = (offIssue ? 1 : 0) + (selfIssue ? 1 : 0);
        bool expectedValid = issueCount == 0;

        v.IsValid.Should().Be(expectedValid);
        v.Issues.Should().HaveCount(issueCount);
        v.Issues.Any(i => i.Contains("off-target")).Should().Be(offIssue);
        v.Issues.Any(i => i.StartsWith("Self-complementarity")).Should().Be(selfIssue);
        v.Issues.Should().NotContain("Potential secondary structure formation");
        v.Warnings.Contains("Potential secondary structure formation").Should().Be(expHairpin);
    }

    /// <summary>
    /// Interaction witness: at exact matching the suffix-tree CheckSpecificity score
    /// agrees with ValidateProbe's specificity for every probe × hit-count — the two
    /// specificity routes are consistent.
    /// </summary>
    [Test]
    public void ProbeValid_CheckSpecificity_AgreesWithValidateProbe()
    {
        foreach (var (probe, _, _, _) in ValidationProbes)
            foreach (int k in new[] { 0, 1, 3 })
            {
                string reference = BuildOffTargetReference(probe, k);
                double viaTree = ProbeDesigner.CheckSpecificity(probe, global::SuffixTree.SuffixTree.Build(reference));
                double viaValidate = ThermodynamicProbeDesigner.ValidateProbe(probe, new[] { reference }, maxMismatches: 0).SpecificityScore;
                viaTree.Should().BeApproximately(viaValidate, 1e-9, $"both score {k} exact hits identically");
            }
    }

    /// <summary>
    /// Interaction witness: the mismatch budget governs off-target sensitivity — a site
    /// carrying one substitution is missed under strict matching but found once one
    /// mismatch is permitted. — spec S4.
    /// </summary>
    [Test]
    public void ProbeValid_MaxMismatches_GovernsOffTargetSensitivity()
    {
        string probe = ValidationProbes[1].Seq;
        char[] mutated = probe.ToCharArray();
        mutated[10] = mutated[10] == 'A' ? 'C' : 'A';                 // exactly one substitution
        string reference = new string('C', 30) + new string(mutated) + new string('C', 30);

        ThermodynamicProbeDesigner.ValidateProbe(probe, new[] { reference }, maxMismatches: 0).OffTargetHits
            .Should().Be(0, "strict matching misses a 1-mismatch site");
        ThermodynamicProbeDesigner.ValidateProbe(probe, new[] { reference }, maxMismatches: 1).OffTargetHits
            .Should().Be(1, "allowing one mismatch finds it");
    }

    [Test, Combinatorial]
    public void Hairpin_FoldsPlantedStem_AcrossStemLoopAndEngine(
        [Values(3, 4, 5)] int stem,
        [Values(4, 5, 6)] int loop,
        [Values(HairpinEngine.NnFolder, HairpinEngine.Ntthal)] HairpinEngine engine)
    {
        string seq = Hairpin(stem, loop);

        if (engine == HairpinEngine.NnFolder)
        {
            var hp = PrimerDesigner.FindMostStableHairpin(seq);
            hp.Should().NotBeNull("the planted G·C stem must fold a hairpin");
            hp!.Value.StemLength.Should().Be(stem, "the recovered stem is the planted G·C stem length");
            hp.Value.LoopSize.Should().Be(loop, "the recovered loop is the planted poly-T loop length");
            // Stem monotonicity (the robust NN law: more stem stacks ⇒ lower ΔG37). The ΔG37 SIGN is
            // not asserted — a short stem under a large loop can be net-unstable (ΔG37 > 0).
            if (stem > 3)
                PrimerDesigner.FindMostStableHairpin(Hairpin(3, loop))!.Value.DeltaG37
                    .Should().BeGreaterThan(hp.Value.DeltaG37, "a longer stem lowers ΔG37 (more stacks)");
        }
        else
        {
            var th = ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(seq);
            th.Should().NotBeNull("the ntthal model folds the planted hairpin");
            th!.Value.BasePairs.Should().Be(stem - 1, "ntthal counts the stem's NN stacks (= base pairs − 1)");
            double.IsNaN(th.Value.DeltaG37).Should().BeFalse("ntthal yields a finite ΔG37");
            if (stem > 3)
                ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal(Hairpin(3, loop))!.Value.DeltaG37
                    .Should().BeGreaterThan(th.Value.DeltaG37, "a longer stem lowers ΔG37 under ntthal too");
        }
    }

    /// <summary>
    /// Interaction witness (specialLoop axis): the ntthal model's bundled special-loop bonus makes a
    /// stabilised GNRA tetraloop hairpin at least as stable as the same stem with a generic loop.
    /// </summary>
    [Test]
    public void Hairpin_NtthalSpecialLoop_NoLessStableThanGenericLoop()
    {
        var special = ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGCGAAAGCCC"); // GNRA (GAAA) tetraloop
        var generic = ThermodynamicPrimerDesigner.CalculateHairpinThermodynamicsNtthal("GGGCTTTTGCCC"); // generic poly-T loop
        special.Should().NotBeNull();
        generic.Should().NotBeNull();
        special!.Value.DeltaG37.Should().BeLessThanOrEqualTo(generic!.Value.DeltaG37 + 1e-9,
            "a recognised special tetraloop is no less stable than a generic loop of the same length");
    }

    /// <summary>Deterministic well-mixed ACGT sequence (LCG) so PAMs/sites occur at many positions.</summary>
    private static string DiverseDna(int n, uint seed = 0x12345678u)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-DESIGN-001 — Primer-candidate acceptance (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 22.
    // Dimensions: minLen(3) × maxLen(3) × gcRange(3) × tmRange(3). Grid 3⁴ = 81.
    //
    // Model (Primer3 / Untergasser 2012 acceptance semantics): a candidate primer
    // is accepted iff it lies inside EVERY configured window simultaneously —
    // length ∈ [MinLength, MaxLength], GC% ∈ [MinGcContent, MaxGcContent] and
    // Tm ∈ [MinTm, MaxTm] (plus structural filters held constant here). Acceptance
    // is therefore the logical CONJUNCTION of the per-axis membership tests, and a
    // violation of any single window must surface its own diagnostic.
    //
    // The combinatorial point: the four windows interact multiplicatively. The probe
    // primer (20 nt, 50% GC, Primer3-default Tm 56.43 °C) is placed so each axis straddles its three
    // windows — minLen {15,20,22} and maxLen {18,20,25} bracket length 20 on both
    // sides, the GC windows straddle 50%, the Tm windows straddle 56.43 °C — so every
    // axis genuinely flips acceptance and the AND is exercised across the grid.
    // ═══════════════════════════════════════════════════════════════════════

    // 20-mer, 50% GC, Primer3-default Tm 56.4298 °C (primer3-py 2.3.1 calc_tm), hairpin-free,
    // max homopolymer 2 (verified independently).
    private const string CleanPrimer20 = "GACGCTGTCTGAGACTAGAA";

    private static readonly (double Lo, double Hi)[] GcWindows = { (30, 45), (45, 55), (55, 70) };
    private static readonly (double Lo, double Hi)[] TmWindows = { (45, 55), (55, 65), (57, 67) };

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROBE-DESIGN-001 — Hybridization-probe design (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 24.
    // Dimensions: minLen(3) × maxLen(3) × tmRange(3) × gcRange(3). Grid 3⁴ = 81.
    //
    // Model (microarray/qPCR probe selection): ThermodynamicProbeDesigner.DesignProbes scans every
    // window in [MinLength, MaxLength], keeps those whose GC fraction lies within the
    // configured band (with a ±0.1 early-rejection slack), scores each (penalising GC,
    // Tm, homopolymer, self-complementarity, structure deviations) and returns the top
    // probes. Two of the four windows are HARD constraints on any returned probe —
    // length ∈ [MinLength, MaxLength] (loop bound) and GC ∈ [MinGc−0.1, MaxGc+0.1]
    // (early rejection) — whereas the Tm window is a SOFT scoring term (−0.3 penalty),
    // so out-of-Tm probes can still be returned. Tm itself is the salt-adjusted formula
    // 81.5 + 16.6·log10[Na⁺] + 41·GC − 600/N (≥14 nt) — Howley 1979 / SantaLucia.
    //
    // The combinatorial point: the four windows interact during the scan. The grid
    // asserts the HARD invariants hold for EVERY returned probe under EVERY (len×len×
    // gc×tm) configuration — i.e. the search never emits a probe outside the length or
    // GC bands regardless of how the other windows are set. The soft Tm axis and the GC
    // ±0.1 slack are pinned by dedicated witnesses below.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string ProbeTarget = DiverseDna(400);
    private static readonly (double Lo, double Hi)[] ProbeGcWindows = { (0.30, 0.50), (0.40, 0.60), (0.50, 0.70) };
    private static readonly (double Lo, double Hi)[] ProbeTmWindows = { (50, 65), (55, 70), (60, 80) };

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROBE-VALID-001 — Probe validation (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 25.
    // Spec: tests/TestSpecs/PROBE-VALID-001.md — canonical methods ValidateProbe &
    //       CheckSpecificity.
    //
    // Axis mapping (documented deviation): the canonical validators have NO gc/tm
    // windows — those belong to design (PROBE-DESIGN-001). Probe VALIDATION judges a
    // probe against a reference for (a) cross-hybridization / off-target specificity and
    // (b) self-complementarity / secondary structure. So the genuine validation axes are
    // off-target multiplicity, the probe's self-complementarity, and the self-comp
    // threshold (selfCompMax) — the grid uses offTargetCount × selfCompProbe ×
    // selfCompThreshold.
    //
    // Model (microarray cross-hybridization; off-target editing): specificity is a strict
    // function of exact hit count — 0 hits ⇒ 0.0 (probe doesn't bind), 1 ⇒ 1.0 (unique),
    // N ⇒ 1/N (invariants #4-6). A probe is reported invalid when it accumulates issues
    // (>1 off-target site, self-comp above threshold, secondary structure) UNLESS the
    // lenient clause holds (≤1 hit AND self-comp ≤ 0.4). The combinatorial point: the
    // IsValid decision composes the off-target and self-complementarity checks. Since audit
    // round 2 (B07 A6) the fallback self-complementarity issue is Primer3's alignment-mode
    // self_any/self_end > 12.00; the legacy selfCompThreshold axis must have no effect.
    // ═══════════════════════════════════════════════════════════════════════

    // selfComp = fraction of positions Watson-Crick-paired with the mirror position.
    // SelfAny = Primer3 alignment-mode self_any (dpal.c + align(), compiled): 0.00, 10.00, 20.00 (self_end equal).
    // Hairpin = the fallback stem-loop screen PrimerDesigner.HasHairpinPotential (exact ≥ 4-bp stem, loop ≥ 3 nt;
    // audit round 5, A5-2): poly-A has no stem; GCGC (pos 2) · GCGC (pos 16) and ACGT (pos 0) · ACGT (pos 8) are exact
    // 4-bp stems with ≥ 3-nt loops. (The former private scan tested a loop of exactly 3 nt only and missed both.)
    private static readonly (string Seq, double SelfComp, double SelfAny, bool Hairpin)[] ValidationProbes =
    {
        ("AAAAAAAAAAAAAAAAAAAA", 0.0, 0.0, false),
        ("TGGCGCGGGGTAACGCGCGC", 0.5, 10.0, true),
        ("ACGTACGTACGTACGTACGT", 1.0, 20.0, true),
    };

    /// <summary>Reference holding exactly <paramref name="k"/> exact copies of the probe, C-padded, G-spaced.</summary>
    private static string BuildOffTargetReference(string probe, int k)
    {
        string core = k == 0 ? "" : string.Join(new string('G', 30), Enumerable.Repeat(probe, k));
        return new string('C', 40) + core + new string('C', 40);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-HAIRPIN-001 — Unimolecular hairpin folder / secondary-structure Tm (MolTools)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 241.
    // Spec: tests/TestSpecs/PRIMER-HAIRPIN-001.md
    //       (PrimerDesigner.FindMostStableHairpin / CalculateHairpinThermodynamicsNtthal). ADVANCED §10.
    //
    // Sources: SantaLucia & Hicks (2004) NN hairpin model (stem stacks + loop ΔG); primer3 ntthal.
    //
    // Model: the most stable hairpin's ΔG37 is the NN sum over the stem stacks plus a loop ΔG penalty.
    //   A planted G^stem · T^loop · C^stem hairpin folds a stem of exactly `stem` base pairs around a
    //   `loop`-nt loop. More stem stacks lower ΔG37; a longer loop raises (destabilises) ΔG37.
    //
    // Dimensions: stemLen(3) × loopLen(3) × specialLoop(2). Grid 3×3×2 = 18 (exhaustive).
    //
    // Axis mapping (documented): specialLoop is realised as the FOLDER axis — the plain NN folder
    // (FindMostStableHairpin, no special-loop bonus) vs the primer3 ntthal model
    // (CalculateHairpinThermodynamicsNtthal, which bundles the special tri/tetraloop ΔG bonuses).
    // Both fold the same planted stem; the bonus itself is pinned by a dedicated witness.
    // ═══════════════════════════════════════════════════════════════════════

    public enum HairpinEngine { NnFolder, Ntthal }

    private static string Hairpin(int stem, int loop) => new string('G', stem) + new string('T', loop) + new string('C', stem);
}
