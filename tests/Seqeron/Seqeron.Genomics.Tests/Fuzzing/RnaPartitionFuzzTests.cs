using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using PartitionFunctionResult = Seqeron.Genomics.Analysis.RnaSecondaryStructure.PartitionFunctionResult;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for RNA-PARTITION-001 — the McCASKILL EQUILIBRIUM PARTITION FUNCTION.
/// Unit under test:
///   RnaSecondaryStructure.CalculatePartitionFunction(string, double basePairEnergy, double temperature)
///   → PartitionFunctionResult(double PartitionFunction Z, IReadOnlyDictionary&lt;(int,int),double&gt; P[i,j])
/// in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs — the McCaskill
/// (1990) inside recursion Q/Qᵇ giving Z = Q[1,n] in O(n³)/O(n²), followed by the outside
/// recursion P[i,j] = Qᵇ[i,j]·O[i,j]/Z over the simplified fixed-per-pair energy model.
/// Checklist: docs/checklists/03_FUZZING.md, row 154 (strategy code BE).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Degenerate, boundary and out-of-domain inputs must NEVER fail in an undisciplined
/// way: no IndexOutOfRange off the n×n DP tables (Q, Qᵇ, outside) on length-0/1 spans;
/// no DivideByZero in P = Qᵇ·O/Z when Z is degenerate; no NaN/Inf in Z or any probability
/// (notably at the temperature edge, where RT→0 would blow up exp(−βE)); no FALSE
/// probability mass on a sequence that cannot pair; and no hang / super-cubic blow-up on a
/// long homopolymer. Every input resolves to a well-defined, theory-correct ensemble OR a
/// documented validation outcome (null → ArgumentNullException; T ≤ 0 → ArgumentOutOfRange).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Scope vs the sibling RNA fuzz rows
/// ───────────────────────────────────────────────────────────────────────────
/// Row 71 (RNA-STRUCT-001) covers greedy PredictStructure; rows 72/73 the base-pair /
/// energy surfaces; rows 149–153 the dot-bracket / hairpin / inverted-repeat / base-pair /
/// MFE-fold surfaces. THIS row is the McCaskill PARTITION FUNCTION — the DP SUM over the
/// whole Boltzmann ensemble (Z + every base-pair probability), distinct from the single-
/// optimum MFE fold (RNA-MFE-001) and from the scalar single-structure Boltzmann weight
/// CalculateStructureProbability (a different method — see the dedicated temperature-edge
/// test below that pins which method is RNA-PARTITION-001). BE corners stress the DP tables:
/// empty, single base, all-unpaired homopolymer, sub-hairpin spans.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy code BE = Boundary Exploitation (граничні значення: 0, -1, empty)
/// — docs/checklists/03_FUZZING.md §Description. Row-154 targets:
///   "empty, single base, all-unpaired".
/// ───────────────────────────────────────────────────────────────────────────
///
/// The DOCUMENTED contract under test (docs/algorithms/RnaStructure/RNA_Partition_Function.md):
///   • §3.3 / §6.1: null → ArgumentNullException; T ≤ 0 → ArgumentOutOfRangeException;
///     empty sequence → Z = 1 with an EMPTY probability map (only the empty structure);
///     a sequence with no admissible pair ("AAAA") → Z = 1; pair span ≤ 3 only ("GC") → Z = 1.
///   • INV-01 (§2.4): Z ≥ 1 always — the empty/all-unpaired structure contributes exp(0)=1.
///   • INV-02 (§2.4): every P[i,j] ∈ [0,1] — a sum of Boltzmann probabilities of structures
///     containing (i,j). Per-column the pair-probabilities into a base sum to ≤ 1 (a base can
///     be paired in at most one pair per structure), so each row/column sum of P ≤ 1.
///   • INV-03 (§2.4): with E_bp = 0, Z equals the number of admissible structures (every
///     Boltzmann weight is exp(0)=1). §7.1 worked example: Z("GGGAAACCC", E_bp=0) = 20 and
///     P[0,8] = 6/20 = 0.30.
///   • INV-04 (§2.4): Z strictly increases as E_bp decreases (each pair weight grows).
///   • §2.2 / §5.2: P uses the FULL outside recursion (external + enclosing terms), so a
///     nestable pair receives strictly more than its external-only share.
///   • Keys are 0-based with i &lt; j; only A-U, G-C, G-U pair; T is treated as U.
///
/// Complexity / hang-safety (§4.3): O(n³) time, O(n²) space. The "all-unpaired" /
/// "extreme length" homopolymers are kept MODEST (n ≤ 200) and guarded with [CancelAfter]
/// so a regression that introduced a hang or non-terminating DP would FAIL, not wedge the
/// suite. The known [Explicit] performance baseline (random n=300, §7.1) is NOT invoked here.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaPartitionFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed LOCALLY so generated fuzz inputs are reproducible.</summary>
    private static string RandomRna(int length, int seed)
    {
        const string bases = "ACGU";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>A homopolymer of a single base — no two positions can ever pair (all-unpaired).</summary>
    private static string Homopolymer(char b, int length) => new string(b, length);

    /// <summary>
    /// The WELL-FORMED-RESULT helper for the partition-function ensemble. Pins every
    /// documented invariant that must hold for ANY input (INV-01, INV-02 and the ≤ 1
    /// per-base mass bound): Z is finite and ≥ 1; every P[i,j] is a finite value in [0,1]
    /// with valid 0-based keys (i &lt; j, in range); and the total probability mass entering
    /// any single base is ≤ 1 (a base pairs in at most one pair per structure, so its
    /// summed pairing probability cannot exceed 1).
    /// </summary>
    private static void AssertWellFormedEnsemble(PartitionFunctionResult r, int n, string because)
    {
        double z = r.PartitionFunction;
        double.IsNaN(z).Should().BeFalse($"Z must never be NaN ({because})");
        double.IsInfinity(z).Should().BeFalse($"Z must never be ±∞ ({because})");
        z.Should().BeGreaterThanOrEqualTo(1.0, $"INV-01: Z ≥ 1 — the empty structure always weighs exp(0)=1 ({because})");

        r.BasePairProbabilities.Should().NotBeNull($"the probability map is never null ({because})");

        // Per-base accumulated pairing probability — must not exceed 1 for any base.
        var perBase = new double[Math.Max(n, 0)];

        foreach (var kv in r.BasePairProbabilities)
        {
            (int i, int j) = kv.Key;
            double p = kv.Value;

            i.Should().BeLessThan(j, $"keys are ordered i < j ({because})");
            i.Should().BeGreaterThanOrEqualTo(0, $"keys are 0-based and in range ({because})");
            j.Should().BeLessThan(n, $"keys are 0-based and in range ({because})");

            double.IsNaN(p).Should().BeFalse($"P[{i},{j}] must never be NaN ({because})");
            double.IsInfinity(p).Should().BeFalse($"P[{i},{j}] must never be ±∞ ({because})");
            p.Should().BeInRange(-1e-12, 1.0 + 1e-9,
                $"INV-02: P[{i},{j}] is a probability in [0,1] ({because})");

            perBase[i] += p;
            perBase[j] += p;
        }

        for (int b = 0; b < perBase.Length; b++)
            perBase[b].Should().BeLessThanOrEqualTo(1.0 + 1e-9,
                $"total pairing mass into base {b} cannot exceed 1 (a base pairs at most once per structure) ({because})");
    }

    /// <summary>Sum of all base-pair probabilities — the expected number of base pairs in the ensemble.</summary>
    private static double TotalPairMass(PartitionFunctionResult r) =>
        r.BasePairProbabilities.Values.Sum();

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-PARTITION-001 — McCaskill partition function : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region RNA-PARTITION-001 — McCaskill partition function

    #region BE — Boundary: empty sequence → Z = 1, empty probability map, no DP-table crash

    /// <summary>
    /// BE — "empty": the empty sequence. The contract short-circuits "" to Z = 1 with an
    /// EMPTY probability map BEFORE allocating any DP buffer, so there is no n×n table to
    /// index and no IndexOutOfRange / DivideByZero can occur (RNA_Partition_Function.md
    /// §3.3, §6.1, INV-01). We pin BOTH the exact value (Z = 1, no pairs) and the no-throw
    /// guarantee, and that pathological energy / temperature parameters do not perturb the
    /// empty boundary (the empty structure is the only structure regardless of weights).
    /// </summary>
    [Test]
    public void Empty_ReturnsZOneEmptyProbabilities_NoCrash()
    {
        var act = () => RnaSecondaryStructure.CalculatePartitionFunction("");
        PartitionFunctionResult r = act.Should().NotThrow("empty input is a documented no-op, not an error").Subject;

        r.PartitionFunction.Should().Be(1.0, "empty sequence: only the empty structure exists, Z = 1 (§6.1, INV-01)");
        r.BasePairProbabilities.Should().BeEmpty("no positions exist, so no pair can form (§3.3)");
        AssertWellFormedEnsemble(r, 0, "empty sequence");

        // Pathological energy / temperature on the empty boundary still yields Z = 1, no pairs.
        foreach (double e in new[] { -100.0, 0.0, 100.0 })
        {
            var rr = RnaSecondaryStructure.CalculatePartitionFunction("", basePairEnergy: e);
            rr.PartitionFunction.Should().Be(1.0, $"empty sequence is Z = 1 for any E_bp (here {e})");
            rr.BasePairProbabilities.Should().BeEmpty();
        }
        RnaSecondaryStructure.CalculatePartitionFunction("", temperature: 1e-9)
            .PartitionFunction.Should().Be(1.0, "empty sequence is Z = 1 at any positive temperature");
    }

    #endregion

    #region BE — Boundary: single base → no pair, Z = 1, empty probabilities, no 1×1-table crash

    /// <summary>
    /// BE — "single base": a length-1 sequence. A base cannot pair with itself and there is
    /// no second position, so the ONLY structure is the open chain → Z = 1 with an EMPTY
    /// probability map. The danger this targets is an IndexOutOfRange on a 1×1 DP table or an
    /// off-by-one in the GetQ(i+1, j-1) empty-interval base case. All four bases are pinned so
    /// the result does not depend on which single residue is supplied (RNA_Partition_Function.md
    /// §6.1, INV-01).
    /// </summary>
    [Test]
    public void SingleBase_AllResidues_ReturnsZOneNoPairs()
    {
        foreach (char b in "ACGU")
        {
            var act = () => RnaSecondaryStructure.CalculatePartitionFunction(b.ToString());
            PartitionFunctionResult r = act.Should().NotThrow($"a single '{b}' is a documented no-op").Subject;

            r.PartitionFunction.Should().Be(1.0, $"a single '{b}' has no partner — only the open chain (Z = 1)");
            r.BasePairProbabilities.Should().BeEmpty($"a single '{b}' admits no pair");
            AssertWellFormedEnsemble(r, 1, $"single base '{b}'");
        }
    }

    /// <summary>
    /// BE — sub-hairpin spans: lengths short enough that no pair (i,j) can enclose the
    /// minimum 3-nt hairpin loop (need j − i &gt; 3). "GC"/"GCG"/"GCGC" are perfectly
    /// complementary yet too short to fold, so Z = 1 and the map is empty — a guard against
    /// counting steric-ally forbidden pairs (RNA_Partition_Function.md §6.1, ASM-03).
    /// </summary>
    [Test]
    public void SubHairpinSpans_CannotPair_ReturnsZOne()
    {
        foreach (string tooShort in new[] { "GC", "GCG", "GCGC", "CG", "AU", "GCGCG" })
        {
            PartitionFunctionResult r = RnaSecondaryStructure.CalculatePartitionFunction(tooShort);
            r.PartitionFunction.Should().Be(1.0,
                $"'{tooShort}' cannot enclose a 3-nt hairpin loop (span ≤ 3), so only the open chain exists");
            r.BasePairProbabilities.Should().BeEmpty($"'{tooShort}' has no admissible pair (min hairpin loop = 3)");
            AssertWellFormedEnsemble(r, tooShort.Length, $"sub-hairpin '{tooShort}'");
        }
    }

    #endregion

    #region BE — Boundary: all-unpaired (homopolymer / non-complementary) → Z = 1, no false mass

    /// <summary>
    /// BE — "all-unpaired": a homopolymer (e.g. "AAAA…"). No two positions are complementary,
    /// so NO pair is admissible and the only structure is the open chain → Z = 1 exactly, with
    /// an EMPTY probability map and ZERO total pair mass. This is the core BE target: a guard
    /// against FALSE probability mass and against a DivideByZero / NaN when the probability loop
    /// runs over a degenerate (pair-free) ensemble (RNA_Partition_Function.md §6.1, INV-01).
    /// All-A, all-C, all-G, all-U are pinned across several lengths; [CancelAfter] proves the
    /// DP terminates on the homopolymer.
    /// </summary>
    [Test]
    [CancelAfter(20_000)]
    public void AllUnpaired_Homopolymer_ZIsOneNoMass(CancellationToken _)
    {
        foreach (char b in "ACGU")
        {
            foreach (int len in new[] { 2, 4, 8, 16, 50 })
            {
                string seq = Homopolymer(b, len);
                PartitionFunctionResult r = RnaSecondaryStructure.CalculatePartitionFunction(seq);

                r.PartitionFunction.Should().Be(1.0,
                    $"homopolymer '{b}'×{len} has no complementary pair — only the open chain (Z = 1)");
                r.BasePairProbabilities.Should().BeEmpty(
                    $"homopolymer '{b}'×{len} admits no pair, so NO false probability mass (INV-01)");
                TotalPairMass(r).Should().Be(0.0, $"zero expected pairs in '{b}'×{len}");
                AssertWellFormedEnsemble(r, len, $"homopolymer '{b}'×{len}");
            }
        }
    }

    /// <summary>
    /// BE — all-unpaired via a NON-complementary mixed sequence (not a homopolymer): a run of
    /// only A and G has no Watson-Crick / wobble partner among {A,G} (A pairs only U; G pairs
    /// C and U), so it still cannot fold → Z = 1, no pairs. Confirms the "no false mass"
    /// guarantee is about pairing ADMISSIBILITY, not literal character identity.
    /// </summary>
    [Test]
    public void AllUnpaired_NonComplementaryMix_ZIsOne()
    {
        foreach (string seq in new[] { "AGAGAGAG", "AAGGAAGG", "GAGAGAGAGA", "AGAAGAAGAA" })
        {
            PartitionFunctionResult r = RnaSecondaryStructure.CalculatePartitionFunction(seq);
            r.PartitionFunction.Should().Be(1.0,
                $"'{seq}' has only A/G — no admissible pair among them, so Z = 1 (open chain only)");
            r.BasePairProbabilities.Should().BeEmpty($"'{seq}' admits no pair");
            AssertWellFormedEnsemble(r, seq.Length, $"non-complementary mix '{seq}'");
        }
    }

    #endregion

    #region BE — Boundary: temperature / energy edges → no NaN, no DivideByZero

    /// <summary>
    /// BE — temperature edge (physically meaningful regime). Across the biological / physical
    /// range and well beyond it the ensemble is fully well-formed: Z finite ≥ 1, every
    /// P ∈ [0,1], no NaN/Inf leaking into P = Qᵇ·O/Z. (Distinct from CalculateStructureProbability,
    /// which is a SEPARATE scalar method — this is the McCaskill DP, RNA-PARTITION-001.)
    /// </summary>
    [Test]
    [CancelAfter(20_000)]
    public void TemperatureEdges_FoldableSequence_WellFormed(CancellationToken _)
    {
        const string seq = "GGGGAAAACCCC";
        foreach (double t in new[] { 5.0, 37.0, 100.0, 273.15, 310.15, 373.15, 1000.0, 1e6 })
        {
            var act = () => RnaSecondaryStructure.CalculatePartitionFunction(seq, temperature: t);
            PartitionFunctionResult r = act.Should().NotThrow($"positive temperature {t} K is valid").Subject;
            AssertWellFormedEnsemble(r, seq.Length, $"temperature {t} K");
        }
    }

    /// <summary>
    /// BE — sub-Kelvin temperature edge (the genuine overflow hazard). As T → 0⁺, RT → 0 and
    /// the per-pair Boltzmann weight exp(−βE_bp) for a stabilising pair overflows IEEE double to
    /// +∞ — inherent to a non-log-space partition function (the doc's only precondition is
    /// T &gt; 0, and §2.4 notes Z scales with RT). The fuzzing contract here is the DISCIPLINED
    /// outcome: no exception, and critically NO NaN leaking into any probability via P = Qᵇ·O/Z
    /// (the source guards the probability loop with <c>if (z &gt; 0)</c>, so an overflow to +∞
    /// yields an EMPTY probability map rather than ∞/∞ = NaN). Z may be a well-defined +∞ here,
    /// but it must never be NaN and must remain ≥ 1.
    /// </summary>
    [Test]
    [CancelAfter(20_000)]
    public void TemperatureEdges_SubKelvin_NoNaN_NoThrow(CancellationToken _)
    {
        const string seq = "GGGGAAAACCCC";
        foreach (double t in new[] { 1e-3, 0.01, 0.1, 0.5 })
        {
            var act = () => RnaSecondaryStructure.CalculatePartitionFunction(seq, temperature: t);
            PartitionFunctionResult r = act.Should().NotThrow($"positive temperature {t} K is accepted").Subject;

            double z = r.PartitionFunction;
            double.IsNaN(z).Should().BeFalse($"Z must never be NaN even at {t} K (overflow → well-defined +∞, not NaN)");
            z.Should().BeGreaterThanOrEqualTo(1.0, $"INV-01: Z ≥ 1 at {t} K");

            foreach (var kv in r.BasePairProbabilities.Values)
                double.IsNaN(kv).Should().BeFalse($"no probability may be NaN at {t} K (∞/∞ guard)");
        }
    }

    /// <summary>
    /// BE — energy edge: E_bp = 0 (every Boltzmann weight is exp(0)=1) and extreme magnitudes.
    /// At E_bp = 0, Z must equal the COUNT of admissible structures (INV-03). A large positive
    /// E_bp (destabilising) drives pair weights → 0, so the ensemble collapses toward the open
    /// chain: Z → 1 and pair mass → 0. A large negative E_bp must still yield finite, in-range
    /// probabilities. None of these may produce NaN/Inf (RNA_Partition_Function.md §2.4).
    /// </summary>
    [Test]
    public void EnergyEdges_NoNaN_AndZeroEnergyCountsStructures()
    {
        const string seq = "GGGAAACCC";

        // E_bp = 0 → Z is the structure count; documented worked example = 20.
        PartitionFunctionResult zero = RnaSecondaryStructure.CalculatePartitionFunction(seq, basePairEnergy: 0.0);
        zero.PartitionFunction.Should().BeApproximately(20.0, 1e-9,
            "INV-03 / §7.1: with E_bp = 0, Z('GGGAAACCC') = number of admissible structures = 20");
        AssertWellFormedEnsemble(zero, seq.Length, "E_bp = 0");

        // Strongly destabilising → ensemble collapses to the open chain (Z → 1, mass → 0).
        PartitionFunctionResult hot = RnaSecondaryStructure.CalculatePartitionFunction(seq, basePairEnergy: 100.0);
        hot.PartitionFunction.Should().BeApproximately(1.0, 1e-6,
            "a strongly destabilising E_bp drives every pair weight → 0, so Z → 1");
        TotalPairMass(hot).Should().BeApproximately(0.0, 1e-6, "no pairing mass when pairs are forbidden by energy");
        AssertWellFormedEnsemble(hot, seq.Length, "E_bp = +100");

        // Strongly stabilising → still finite and in-range (no NaN/Inf).
        PartitionFunctionResult cold = RnaSecondaryStructure.CalculatePartitionFunction(seq, basePairEnergy: -10.0);
        AssertWellFormedEnsemble(cold, seq.Length, "E_bp = -10");
    }

    #endregion

    #region BE — Validation: null sequence, non-positive temperature

    /// <summary>
    /// BE — documented validation outcomes (RNA_Partition_Function.md §3.3, §6.1): a null
    /// sequence is an ArgumentNullException; a non-positive temperature (0 or negative Kelvin
    /// is physically meaningless, would make RT = 0 → exp(±∞) → NaN) is an
    /// ArgumentOutOfRangeException. These are explicit contracts, NOT silent NaN leakage.
    /// </summary>
    [Test]
    public void Validation_NullAndNonPositiveTemperature_Throw()
    {
        var nullAct = () => RnaSecondaryStructure.CalculatePartitionFunction(null!);
        nullAct.Should().Throw<ArgumentNullException>("null sequence is a contract violation (§3.3)");

        foreach (double t in new[] { 0.0, -1.0, -310.15, double.NegativeInfinity })
        {
            var act = () => RnaSecondaryStructure.CalculatePartitionFunction("GGGAAACCC", temperature: t);
            act.Should().Throw<ArgumentOutOfRangeException>(
                $"temperature {t} K is non-positive — rejected, never a NaN leak (§3.3, §6.1)");
        }
    }

    #endregion

    #region POSITIVE sanity — a foldable sequence yields Z > 1 and real probability mass

    /// <summary>
    /// POSITIVE sanity (anti-rubber-stamp): a clearly FOLDABLE sequence "GGGGAAAACCCC" (a
    /// 4-bp G·C stem closing a 4-nt loop) must have Z &gt; 1 and assign STRICTLY POSITIVE
    /// base-pair probability (all in [0,1]) to the plausible stem pairs (the nested G·C pairs
    /// flanking the AAAA loop), proving the DP actually counts structures rather than returning
    /// the degenerate Z = 1. We assert at least one stem pair has meaningful mass and that the
    /// most-probable stem pair is one of the canonical (i, n−1−i) G·C closures. We do NOT
    /// hard-code Vienna numbers (simplified model, ASM-01); we assert the documented INVARIANTS.
    /// </summary>
    [Test]
    [CancelAfter(20_000)]
    public void Foldable_HasZGreaterThanOneAndStemMass(CancellationToken _)
    {
        const string seq = "GGGGAAAACCCC"; // n = 12; G(0..3) complementary to C(8..11)
        int n = seq.Length;

        PartitionFunctionResult r = RnaSecondaryStructure.CalculatePartitionFunction(seq);
        AssertWellFormedEnsemble(r, n, "foldable GGGGAAAACCCC");

        r.PartitionFunction.Should().BeGreaterThan(1.0,
            "a foldable stem-loop has at least one pairing structure beyond the open chain → Z > 1");
        r.BasePairProbabilities.Should().NotBeEmpty("a foldable sequence admits base pairs");
        TotalPairMass(r).Should().BeGreaterThan(0.0, "the foldable ensemble carries real pairing mass");

        // The canonical stem pairs are the symmetric G·C closures (i, n-1-i) for i = 0..3.
        var canonicalStem = new[] { (0, 11), (1, 10), (2, 9), (3, 8) };
        foreach (var pair in canonicalStem)
        {
            r.BasePairProbabilities.Should().ContainKey(pair,
                $"the complementary G·C closure {pair} is an admissible stem pair");
            r.BasePairProbabilities[pair].Should().BeGreaterThan(0.0,
                $"the stem pair {pair} occurs in part of the ensemble → strictly positive probability");
        }

        // The single most probable pair must be a genuine stem closure, not noise.
        var top = r.BasePairProbabilities.OrderByDescending(kv => kv.Value).First();
        canonicalStem.Should().Contain(top.Key,
            "the highest-probability base pair is one of the canonical G·C stem closures");
    }

    /// <summary>
    /// POSITIVE sanity — documented WORKED EXAMPLE (§7.1): with E_bp = 0 every weight is 1,
    /// so Z("GGGAAACCC") is exactly the number of admissible non-crossing structures = 20, and
    /// the outermost pair (0,8) appears in 6 of them → P[0,8] = 6/20 = 0.30. Pinning this exact
    /// rational value verifies BOTH the inside count AND the FULL outside recursion (the
    /// external-term-only formula would under-report (0,8); §2.2, §5.2).
    /// </summary>
    [Test]
    public void WorkedExample_GGGAAACCC_ZeroEnergy_MatchesDoc()
    {
        PartitionFunctionResult r = RnaSecondaryStructure.CalculatePartitionFunction("GGGAAACCC", basePairEnergy: 0.0);

        r.PartitionFunction.Should().BeApproximately(20.0, 1e-9,
            "§7.1: Z('GGGAAACCC', E_bp=0) = 20 admissible structures");
        r.BasePairProbabilities.Should().ContainKey((0, 8));
        r.BasePairProbabilities[(0, 8)].Should().BeApproximately(0.30, 1e-9,
            "§7.1: 6 of 20 structures contain (0,8) → P[0,8] = 6/20 = 0.30 (full outside recursion)");
        AssertWellFormedEnsemble(r, 9, "worked example GGGAAACCC");
    }

    /// <summary>
    /// INV-04 (§2.4): Z strictly increases as E_bp decreases — each pair weight exp(−βE_bp)
    /// grows and Z is a positive combination of them. A monotone sweep of decreasing E_bp on a
    /// foldable sequence must give strictly increasing Z (and always ≥ 1).
    /// </summary>
    [Test]
    public void Monotonicity_ZIncreasesAsEnergyDecreases()
    {
        const string seq = "GGGGAAAACCCC";
        double[] energies = { 2.0, 1.0, 0.0, -1.0, -2.0, -3.0 }; // strictly decreasing
        double prevZ = double.NegativeInfinity;
        foreach (double e in energies)
        {
            double z = RnaSecondaryStructure.CalculatePartitionFunction(seq, basePairEnergy: e).PartitionFunction;
            z.Should().BeGreaterThanOrEqualTo(1.0, "INV-01 holds at every energy");
            z.Should().BeGreaterThan(prevZ, $"INV-04: Z strictly increases as E_bp decreases (at E_bp = {e})");
            prevZ = z;
        }
    }

    #endregion

    #region BE — Robustness sweep: random fuzz inputs across lengths never violate invariants

    /// <summary>
    /// BE robustness sweep: a deterministic stream of random RNAs across a range of lengths
    /// (including the boundary lengths 0/1/2 and a modest homopolymer-comparable upper bound)
    /// must NEVER violate the documented invariants — no throw, Z finite ≥ 1, every P ∈ [0,1],
    /// per-base mass ≤ 1, valid keys. [CancelAfter] proves O(n³) termination on every draw.
    /// This is the broad anti-fragility net behind the targeted BE cases above.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void RandomSweep_AllLengths_NeverViolatesInvariants(CancellationToken _)
    {
        foreach (int len in new[] { 0, 1, 2, 3, 5, 8, 13, 21, 40, 80, 120 })
        {
            for (int seed = 0; seed < 8; seed++)
            {
                string seq = RandomRna(len, seed * 131 + len);
                var act = () => RnaSecondaryStructure.CalculatePartitionFunction(seq);
                PartitionFunctionResult r =
                    act.Should().NotThrow($"len={len} seed={seed} '{seq}' must not throw").Subject;
                AssertWellFormedEnsemble(r, len, $"random len={len} seed={seed} '{seq}'");
            }
        }
    }

    #endregion

    #endregion
}
