// RNA-HAIRPIN-001 — Hairpin Loop and Stem Free-Energy Calculation (Turner 2004)
// Fuzz tests (strategy BE = Boundary Exploitation).
// Algorithm doc: docs/algorithms/RnaStructure/Hairpin_Energy_Calculation.md
// Canonical tests: tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_HairpinEnergy_Tests.cs (RNA-HAIRPIN-001)
// Evidence: docs/Evidence/RNA-HAIRPIN-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-HAIRPIN-001.md
// Source: Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH (2004).
//         Proc. Natl. Acad. Sci. USA 101:7287-7292. doi:10.1073/pnas.0401799101
//         Parameters from NNDB Turner 2004 (rna.urmc.rochester.edu/NNDB/turner04).

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for RNA-HAIRPIN-001 — the Turner 2004 nearest-neighbor free-energy
/// calculation of the two elementary stem-loop motifs:
///   • <see cref="RnaSecondaryStructure.CalculateHairpinLoopEnergy(string,char,char,bool)"/>
///     — ΔG°37 of a hairpin LOOP (terminal single-stranded loop closed by a base pair).
///   • <see cref="RnaSecondaryStructure.CalculateStemEnergy(string,IReadOnlyList{BasePair})"/>
///     — ΔG°37 of a STEM (double-stranded helix): Σ of the P−1 nearest-neighbor stacking
///     terms over P stacked base pairs, plus a +0.45 per-AU/GU-end penalty.
/// Both live in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs.
///
/// This is the DISTINCT energy unit (RNA-HAIRPIN-001, checklist row 150). The hairpin/stem
/// MOTIF DETECTOR <c>FindStemLoops</c> is the separate unit RNA-STEMLOOP-001 (row 72,
/// covered by RnaStructureFuzzTests), so this file is scoped strictly to the ENERGY
/// contract — in particular the loop&lt;minLoop, no-stem and empty boundary facets.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs and asserts the code
/// NEVER fails in an undisciplined way: no hang, no unhandled runtime exception
/// (IndexOutOfRange off an empty/short loop string or empty base-pair list,
/// ArgumentOutOfRange from internal indexing, DivideByZero, NaN/Inf), and no nonsense
/// (non-finite) output. Every input must resolve EITHER to a well-defined, theory-correct
/// ΔG°37 value OR to a documented, intentional outcome.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy BE = Boundary Exploitation — targets "loop&lt;minLoop, no stem, empty"
/// ───────────────────────────────────────────────────────────────────────────
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE), row 150.
///
///   • loop &lt; minLoop — the steric minimum hairpin loop is 3 nt. The nearest-neighbor
///       rules PROHIBIT hairpin loops with fewer than 3 nucleotides; loops &lt; 3 nt are
///       sterically impossible and do not form (Hairpin_Energy_Calculation.md §2.4 INV-02,
///       §6.1; NNDB hairpin.html). The implementation returns a deliberately PROHIBITIVE
///       exact sentinel of 100.0 so a downstream optimizer can never select such a loop —
///       NOT a normal low energy, never a crash on the 0/1/2-nt boundary.
///
///   • no stem — a helix of FEWER THAN 2 base pairs has P−1 = 0 stacking terms. An EMPTY
///       base-pair list therefore yields exactly 0 (Hairpin_Energy_Calculation.md §2.4
///       INV-04, §3.3, §6.1); a SINGLE base pair has 0 stacking terms too (only a possible
///       terminal-end penalty remains). No index walk off an empty/one-element list.
///
///   • empty — an EMPTY loop string is size 0 &lt; 3, so it is prohibited (100.0), never an
///       IndexOutOfRange on loop[0]/loop[^1]; an EMPTY stem (no pairs) is 0. Both degenerate
///       inputs are well-defined, never DivideByZero / IndexOutOfRange.
///
/// Watched failure modes: IndexOutOfRange on a too-short loop string, off-by-one at the
/// minLoop=3 boundary (a 2-nt loop wrongly scored, or a 3-nt loop wrongly prohibited),
/// false sub-minimal "hairpin" energies, and an index walk off an empty base-pair list.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Theory-correct contract asserted (Hairpin_Energy_Calculation.md §2.4, §3.3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   • INV-02 — loops &lt; 3 nt return the EXACT prohibitive sentinel 100.0 (not merely "large").
///   • The 3-nt loop is the FIRST sterically-allowed loop: it returns a finite, NON-prohibitive
///     value (initiation(3) 5.4 [+ all-C penalty]) — the off-by-one boundary is exact.
///   • INV-04 — CalculateStemEnergy of 0 base pairs is exactly 0; a single pair sums 0 stacks.
///   • INV-01 — deterministic in (loop, closing pair, special-GU flag) / (sequence, pairs).
///   • Every returned ΔG°37 is FINITE (no NaN/Inf) for every boundary / fuzzed input.
///   • POSITIVE sanity: a real, stable hairpin (a complementary stem + a ≥3-nt loop) scores a
///     NEGATIVE total ΔG, and a sub-minimal (≤2-nt) loop in the same context is rejected.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaHairpinFuzzTests
{
    private const double Tol = 1e-9;

    /// <summary>The documented steric minimum hairpin loop size (NNDB Turner 2004).</summary>
    private const int MinLoop = 3;

    /// <summary>The prohibitive sentinel returned for loops &lt; 3 nt (INV-02).</summary>
    private const double ProhibitiveSentinel = 100.0;

    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static Random Rng(int seed) => new(seed);

    private static string RandomLoop(Random rng, int length)
    {
        const string bases = "ACGU";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    private static char RandomBase(Random rng)
    {
        const string bases = "ACGU";
        return bases[rng.Next(bases.Length)];
    }

    /// <summary>
    /// Asserts a hairpin-loop ΔG°37 is well-formed for the given loop size: it is FINITE
    /// (never NaN/Inf) and respects the steric floor — a loop &lt; 3 nt MUST be the exact
    /// prohibitive sentinel (INV-02), and a loop ≥ 3 nt MUST be a normal (non-prohibitive)
    /// finite value. This is the loop-side "well-formed result" check.
    /// </summary>
    private static void AssertWellFormedHairpinEnergy(double dg, int loopSize)
    {
        double.IsNaN(dg).Should().BeFalse("hairpin ΔG°37 must never be NaN");
        double.IsInfinity(dg).Should().BeFalse("hairpin ΔG°37 must never be infinite");

        if (loopSize < MinLoop)
        {
            dg.Should().Be(ProhibitiveSentinel,
                "loops < 3 nt are sterically prohibited and return the exact 100.0 sentinel (INV-02)");
        }
        else
        {
            dg.Should().BeLessThan(ProhibitiveSentinel,
                "a sterically-allowed loop (≥ 3 nt) must score a normal, non-prohibitive energy");
        }
    }

    /// <summary>
    /// Asserts a stem ΔG°37 is well-formed: FINITE, and for &lt; 2 base pairs there are
    /// zero stacking terms so the magnitude can only come from the (bounded) terminal-end
    /// penalty — in particular an empty list is exactly 0 (INV-04).
    /// </summary>
    private static void AssertWellFormedStemEnergy(double dg, int pairCount)
    {
        double.IsNaN(dg).Should().BeFalse("stem ΔG°37 must never be NaN");
        double.IsInfinity(dg).Should().BeFalse("stem ΔG°37 must never be infinite");

        if (pairCount == 0)
            dg.Should().Be(0.0, "a stem of 0 base pairs has no stacking and no terminal pairs (INV-04)");
    }

    #endregion

    #region RNA-HAIRPIN-001 — Hairpin Loop and Stem Free Energy (Turner 2004)

    // ─────────────────────────────────────────────────────────────────────────
    // BE: loop < minLoop — the steric floor (loops < 3 nt prohibited → exact 100.0)
    // ─────────────────────────────────────────────────────────────────────────

    // A 2-nt loop (the boundary just below the minimum) is prohibited: exact sentinel,
    // never an off-by-one normal score, never an IndexOutOfRange on loop[0]/loop[^1].
    [Test]
    public void HairpinLoop_TwoNtLoop_IsProhibited_ExactSentinel()
    {
        double dg = CalculateHairpinLoopEnergy("AA", 'G', 'C');
        dg.Should().Be(ProhibitiveSentinel,
            "a 2-nt loop is below the 3-nt steric floor and returns the exact prohibitive 100.0 (INV-02)");
    }

    // A 1-nt loop is likewise prohibited.
    [Test]
    public void HairpinLoop_OneNtLoop_IsProhibited_ExactSentinel()
    {
        double dg = CalculateHairpinLoopEnergy("A", 'G', 'C');
        dg.Should().Be(ProhibitiveSentinel,
            "a 1-nt loop is below the 3-nt steric floor and returns the exact prohibitive 100.0 (INV-02)");
    }

    // Off-by-one boundary, both sides: 0,1,2-nt loops are ALL prohibited (exact 100.0),
    // and the FIRST allowed size (3 nt) is a normal, non-prohibitive finite value.
    // Closing pair varied to keep the assertion sequence-robust at the boundary.
    [Test]
    public void HairpinLoop_MinLoopBoundary_SubMinimalProhibited_ThreeNtAllowed()
    {
        var rng = Rng(150_001);
        for (int trial = 0; trial < 64; trial++)
        {
            char c5 = RandomBase(rng);
            char c3 = RandomBase(rng);

            for (int size = 0; size <= 2; size++)
            {
                string loop = RandomLoop(rng, size);
                double sub = CalculateHairpinLoopEnergy(loop, c5, c3);
                sub.Should().Be(ProhibitiveSentinel,
                    $"a {size}-nt loop is sub-minimal and must be the exact prohibitive 100.0 sentinel (INV-02)");
            }

            string threeNt = RandomLoop(rng, MinLoop);
            double allowed = CalculateHairpinLoopEnergy(threeNt, c5, c3);
            AssertWellFormedHairpinEnergy(allowed, MinLoop);
        }
    }

    // Fuzz across the whole low-size range: for any random loop string of size 0..2 and any
    // closing pair, the result is the exact sentinel and never crashes; size ≥ 3 is finite
    // and non-prohibitive — the steric floor holds uniformly, no false sub-minimal hairpin.
    [Test]
    public void HairpinLoop_FuzzedShortLoops_NeverFalseSubMinimalHairpin()
    {
        var rng = Rng(150_002);
        for (int trial = 0; trial < 400; trial++)
        {
            int size = rng.Next(0, 6); // 0..5 — straddles the 3-nt floor
            string loop = RandomLoop(rng, size);
            char c5 = RandomBase(rng);
            char c3 = RandomBase(rng);

            double dg = 0;
            Action act = () => dg = CalculateHairpinLoopEnergy(loop, c5, c3);
            act.Should().NotThrow("a fuzzed short loop must never crash the energy calculation");

            AssertWellFormedHairpinEnergy(dg, size);
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE: empty — empty loop string / empty stem
    // ─────────────────────────────────────────────────────────────────────────

    // An EMPTY loop string is size 0 < 3 → prohibited; must NOT IndexOutOfRange on loop[0].
    [Test]
    public void HairpinLoop_EmptyLoopString_IsProhibited_NoCrash()
    {
        double dg = 0;
        Action act = () => dg = CalculateHairpinLoopEnergy("", 'G', 'C');
        act.Should().NotThrow("an empty loop string must be handled, not throw IndexOutOfRange on loop[0]");
        dg.Should().Be(ProhibitiveSentinel,
            "an empty (size-0) loop is below the 3-nt floor → exact prohibitive 100.0 (INV-02)");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE: no stem — empty / single-pair base-pair list (P−1 = 0 stacks)
    // ─────────────────────────────────────────────────────────────────────────

    // An EMPTY base-pair list has no stacking and no terminal pair → exactly 0 (INV-04).
    // Must not index off an empty list (basePairs[0] / basePairs[^1]).
    [Test]
    public void StemEnergy_EmptyBasePairs_IsZero_NoCrash()
    {
        double dg = 0;
        Action act = () => dg = CalculateStemEnergy("ACGU", new List<BasePair>());
        act.Should().NotThrow("an empty base-pair list must not index off the empty list");
        dg.Should().Be(0.0, "a stem of 0 base pairs has P−1 = 0 stacks and no terminal pairs → 0 (INV-04)");
    }

    // A SINGLE base pair has P−1 = 0 stacking terms; the only contribution is the bounded
    // terminal-end penalty. Result must be finite and equal to its end-penalty contribution:
    // a G-C single pair (no AU/GU end) → 0; an A-U single pair (two AU ends, same pair) →
    // +0.45 twice = +0.90. This isolates the no-stacking ("no stem") path.
    [Test]
    public void StemEnergy_SingleBasePair_NoStacks_OnlyEndPenalty()
    {
        var gc = new List<BasePair> { new(0, 3, 'G', 'C', BasePairType.WatsonCrick) };
        double dgGc = CalculateStemEnergy("GCGC", gc);
        dgGc.Should().Be(0.0,
            "a single G-C pair has 0 stacking terms and no AU/GU end → ΔG = 0");

        var au = new List<BasePair> { new(0, 3, 'A', 'U', BasePairType.WatsonCrick) };
        double dgAu = CalculateStemEnergy("AUAU", au);
        dgAu.Should().Be(0.90,
            "a single A-U pair has 0 stacking terms; both helix ends are this AU pair → +0.45 twice = +0.90");
    }

    // Fuzz the no-stem / minimal-stem corner: a single random canonical pair always yields a
    // FINITE energy that comes only from the bounded end penalty (never NaN/Inf, never a crash).
    [Test]
    public void StemEnergy_FuzzedSinglePair_FiniteEndPenaltyOnly()
    {
        (char b1, char b2, BasePairType t)[] canonical =
        {
            ('A', 'U', BasePairType.WatsonCrick), ('U', 'A', BasePairType.WatsonCrick),
            ('G', 'C', BasePairType.WatsonCrick), ('C', 'G', BasePairType.WatsonCrick),
            ('G', 'U', BasePairType.Wobble),      ('U', 'G', BasePairType.Wobble),
        };

        var rng = Rng(150_003);
        for (int trial = 0; trial < 200; trial++)
        {
            var (b1, b2, t) = canonical[rng.Next(canonical.Length)];
            var pairs = new List<BasePair> { new(0, 3, b1, b2, t) };

            double dg = 0;
            Action act = () => dg = CalculateStemEnergy($"{b1}xx{b2}", pairs);
            act.Should().NotThrow("a single-pair stem must never crash");
            AssertWellFormedStemEnergy(dg, pairCount: 1);
            // P−1 = 0 stacks: the magnitude is bounded by the at-most-two end penalties.
            Math.Abs(dg).Should().BeLessThanOrEqualTo(0.90 + Tol,
                "a single pair contributes 0 stacks; only ≤ two +0.45 end penalties remain");
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Robustness sweep — fuzzed loops + closing pairs never crash, always finite
    // ─────────────────────────────────────────────────────────────────────────

    // Broad fuzz across loop sizes (incl. the 3-nt floor, special-loop sizes 4/6 and long
    // loops hitting the log-extrapolation branch) and all closing-pair combinations: the
    // calculation is total — every input returns a finite, steric-floor-respecting value.
    [Test]
    public void HairpinLoop_BroadFuzz_AlwaysFiniteAndRespectsFloor()
    {
        var rng = Rng(150_004);
        for (int trial = 0; trial < 600; trial++)
        {
            int size = rng.Next(0, 40); // includes < 3 (prohibited) and > 30 (extrapolation)
            string loop = RandomLoop(rng, size);
            char c5 = RandomBase(rng);
            char c3 = RandomBase(rng);
            bool guFlag = rng.Next(2) == 0;

            double dg = 0;
            Action act = () => dg = CalculateHairpinLoopEnergy(loop, c5, c3, guFlag);
            act.Should().NotThrow("a fuzzed hairpin loop must never crash");
            AssertWellFormedHairpinEnergy(dg, size);
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Determinism (INV-01)
    // ─────────────────────────────────────────────────────────────────────────

    [Test]
    public void HairpinLoop_FuzzedInputs_AreDeterministic()
    {
        var rng = Rng(150_005);
        for (int trial = 0; trial < 200; trial++)
        {
            int size = rng.Next(0, 12);
            string loop = RandomLoop(rng, size);
            char c5 = RandomBase(rng);
            char c3 = RandomBase(rng);
            bool guFlag = rng.Next(2) == 0;

            double a = CalculateHairpinLoopEnergy(loop, c5, c3, guFlag);
            double b = CalculateHairpinLoopEnergy(loop, c5, c3, guFlag);
            a.Should().Be(b, "hairpin energy must be deterministic in its inputs (INV-01)");
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // POSITIVE sanity — a real hairpin is detected/stable; a sub-minimal loop is rejected
    // ─────────────────────────────────────────────────────────────────────────

    // The canonical worked hairpin "GGGGAAAACCCC": a 4-bp G-C stem (GGGG / CCCC) enclosing a
    // 4-nt AAAA loop closed by a G-C pair. FindStemLoops (RNA-STEMLOOP-001) detects it with the
    // documented stem span and loop size, and — the RNA-HAIRPIN-001 contract under test — the
    // hairpin LOOP and STEM energies it carries are both well-formed, the loop is non-prohibitive,
    // and the TOTAL ΔG is NEGATIVE (the strong G-C stem dominates → a stable, real hairpin).
    [Test]
    public void PositiveSanity_ClearHairpin_DetectedWithNegativeStableEnergy()
    {
        const string seq = "GGGGAAAACCCC";
        var loops = FindStemLoops(seq, minStemLength: 3, minLoopSize: 3, maxLoopSize: 10).ToList();

        loops.Should().NotBeEmpty("GGGGAAAACCCC has a clear 4-bp G-C stem enclosing a 4-nt loop");

        var best = loops.OrderBy(sl => sl.TotalFreeEnergy).First();
        best.Loop.Type.Should().Be(LoopType.Hairpin, "the detected motif is a terminal hairpin loop");
        best.Loop.Size.Should().BeGreaterThanOrEqualTo(MinLoop, "the loop must respect the 3-nt steric floor");

        // The hairpin LOOP energy recomputed directly from the detected loop + closing pair is the
        // well-formed, non-prohibitive RNA-HAIRPIN-001 value.
        char closing5 = seq[best.Stem.End5Prime];
        char closing3 = seq[best.Stem.Start3Prime];
        double loopEnergy = CalculateHairpinLoopEnergy(best.Loop.Sequence, closing5, closing3);
        AssertWellFormedHairpinEnergy(loopEnergy, best.Loop.Size);

        // The STEM energy is well-formed; the strong G-C helix is stabilising (negative).
        double stemEnergy = CalculateStemEnergy(seq, best.Stem.BasePairs);
        AssertWellFormedStemEnergy(stemEnergy, best.Stem.BasePairs.Count);
        stemEnergy.Should().BeLessThan(0.0, "a multi-bp G-C stem is stabilising (negative ΔG)");

        // The full stem-loop ΔG is negative: a stable, real hairpin.
        best.TotalFreeEnergy.Should().BeLessThan(0.0,
            "a clear G-C-stemmed hairpin with a ≥3-nt loop is net stabilising (negative total ΔG)");
    }

    // The complement of the positive case: shrink the loop below the steric floor (2-nt AA loop
    // with the SAME G-C closing context) and the hairpin loop energy is REJECTED (prohibitive),
    // proving the boundary discriminates a real hairpin from a sterically-impossible one.
    [Test]
    public void PositiveSanity_SubMinimalLoop_RejectedAgainstAllowedThreeNt()
    {
        double twoNt   = CalculateHairpinLoopEnergy("AA",  'G', 'C'); // below floor
        double threeNt = CalculateHairpinLoopEnergy("AAA", 'G', 'C'); // first allowed

        twoNt.Should().Be(ProhibitiveSentinel,
            "a 2-nt loop is sterically impossible → exact prohibitive 100.0 (INV-02)");
        threeNt.Should().BeLessThan(ProhibitiveSentinel,
            "a 3-nt loop is the first allowed size → a normal, non-prohibitive energy");
        threeNt.Should().BeLessThan(twoNt,
            "the allowed 3-nt loop is far more favourable than the prohibited 2-nt loop");
    }

    #endregion
}
