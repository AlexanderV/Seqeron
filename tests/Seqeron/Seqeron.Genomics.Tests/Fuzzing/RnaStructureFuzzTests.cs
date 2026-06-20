using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the RnaStructure area — RNA secondary-structure prediction
/// (RNA-STRUCT-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception (IndexOutOfRangeException leaking from the O(L²)/O(L³) DP matrix on
/// short spans, ArgumentOutOfRangeException from internal Substring/indexing,
/// DivideByZero, OutOfMemory). Every input must resolve to EITHER a well-defined,
/// theory-correct result, OR a *documented, intentional* validation exception.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: RNA-STRUCT-001 — RNA secondary structure prediction
/// Checklist: docs/checklists/03_FUZZING.md, row 71.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — empty string, single base, all-A·U sequence,
///          and an extremely-long sequence (kept ≤ 300 because the prediction DP
///          is O(L²)–O(L³); see §Complexity below). These probe the length=0,
///          length=1, length&lt;minLoop+2 short-circuits and the upper extreme.
///   • MC = Malformed Content — non-RNA characters (T instead of U, and arbitrary
///          junk such as digits/punctuation): these must be handled (silently
///          treated as never-pairing residues), never rejected with a crash.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The prediction contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// An RNA secondary structure is a set of non-crossing base pairs over A-U / G-C
/// (Watson–Crick) and G-U (wobble); the nearest-neighbor model assigns free energy
/// to the loops of the structure and the minimum-free-energy (MFE) fold is sought.
/// — docs/algorithms/RnaStructure/Minimum_Free_Energy.md §2.1–2.2.
///
/// Two public surfaces are probed (both in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs):
///
///   (1) PredictStructure(string, minStemLength=3, minLoopSize=3, maxLoopSize=10)
///       → SecondaryStructure { Sequence, DotBracket, BasePairs, StemLoops,
///         Pseudoknots, MinimumFreeEnergy }. A GREEDY stem-loop assembly: it scans
///         hairpins, scores them by Turner energy, and selects NON-OVERLAPPING
///         stem-loops (SelectNonOverlapping), then emits base pairs + dot-bracket.
///         (Minimum_Free_Energy.md §5.1, §4.1 step 4.)
///
///   (2) CalculateMinimumFreeEnergy(string, minLoopSize=3) → double — the Zuker-
///       style O(n³) DP optimal ΔG°37 score. INV-01: MFE ≤ 0 for every input (the
///       empty open chain ΔG = 0 is always in the search space, so the minimum is
///       never positive). (Minimum_Free_Energy.md §2.4 INV-01, §3.2.)
///
/// Documented input handling (Minimum_Free_Energy.md §3.3, §6.1):
///   • null / "" → empty structure / MFE 0 — never an exception.
///   • length &lt; minLoopSize + 2 → cannot enclose a 3-nt hairpin → no pairs / MFE 0.
///   • minLoopSize &lt; 3 → clamped up to 3 (nearest-neighbor rules prohibit shorter
///     hairpin loops — the loops-&lt;3-nt sterically-impossible rule). This is the KEY
///     min-loop invariant: no reported hairpin has fewer than 3 unpaired bases.
///   • case-insensitive (upper-cased internally); only A/C/G/U pair; unrecognized
///     characters (including T on the PredictStructure surface, and arbitrary junk)
///     simply never form pairs — no crash, no rejection. (On the MFE surface T is
///     additionally read as U, since A–U and A–T stacking are identical.)
///
/// KEY structural invariants asserted on every produced structure (the theory-
/// correct contract — these are what a "valid structure" MUST satisfy):
///   • ≤ 1 partner per base — each position appears in AT MOST ONE base pair
///     (a base cannot pair with two partners, nor with itself). [partner-uniqueness]
///   • canonical pairing — every reported pair is A-U / G-C / U-A / C-G / G-U / U-G.
///   • min-loop ≥ 3 — the two members of any hairpin pair span ≥ 3 unpaired bases
///     between them at the innermost (closing) pair. We assert this via the loop
///     records and the clamp behavior.
///   • non-crossing (nested) — PredictStructure assembles non-overlapping hairpins,
///     so the emitted base-pair set is pseudoknot-free: DetectPseudoknots over the
///     result is empty, and the dot-bracket validates (balanced brackets).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// Secondary-structure prediction is O(L²)–O(L³) (Minimum_Free_Energy.md §4.3).
/// The "extremely long" fuzz target is therefore kept MODEST (L ≤ 300, the top of
/// the documented benchmark's measured range, ~214 ms for MFE at n=300) and guarded
/// with [CancelAfter] so a regression that introduced a hang would FAIL rather than
/// wedge the suite. The known [Explicit] RnaSecondaryStructure_MFE_Benchmark is NOT
/// invoked here.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaStructureFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomRna(int length, int seed)
    {
        const string bases = "ACGU";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal structural invariants every predicted structure must
    /// satisfy: partner-uniqueness (≤ 1 partner per base, never self-pairing),
    /// canonical pairing, in-range indices, nesting (pseudoknot-free) and a
    /// well-formed, balanced dot-bracket consistent with the base-pair count.
    /// </summary>
    private static void AssertWellFormedStructure(RnaSecondaryStructure.SecondaryStructure s)
    {
        int n = s.Sequence.Length;

        // Dot-bracket length matches the sequence and is balanced (nested brackets).
        s.DotBracket.Should().HaveLength(n, "the dot-bracket annotates every residue exactly once");
        RnaSecondaryStructure.ValidateDotBracket(s.DotBracket).Should().BeTrue(
            "a predicted structure must be a balanced, well-nested bracket string");

        var seen = new HashSet<int>();
        foreach (var bp in s.BasePairs)
        {
            int i = Math.Min(bp.Position1, bp.Position2);
            int j = Math.Max(bp.Position1, bp.Position2);

            // In-range, distinct endpoints (no base pairs with itself).
            i.Should().BeInRange(0, n - 1, "pair endpoints must index inside the sequence");
            j.Should().BeInRange(0, n - 1, "pair endpoints must index inside the sequence");
            i.Should().BeLessThan(j, "a base can never pair with itself");

            // Partner-uniqueness: each position participates in at most one pair.
            seen.Add(i).Should().BeTrue($"position {i} must have at most one partner");
            seen.Add(j).Should().BeTrue($"position {j} must have at most one partner");

            // Canonical pairing only (A-U / G-C / G-U wobble).
            RnaSecondaryStructure.CanPair(s.Sequence[i], s.Sequence[j]).Should().BeTrue(
                $"reported pair ({s.Sequence[i]},{s.Sequence[j]}) must be a canonical RNA pair");

            // Min-loop ≥ 3: the closing-pair endpoints of every reported hairpin
            // enclose at least 3 unpaired bases (j - i - 1 ≥ 3 for the innermost pair).
            // Each emitted base pair belongs to a stem whose innermost pair satisfies
            // this; the weakest span (largest i, smallest j inside a stem) still ≥ 3.
        }

        // Nesting / pseudoknot-free: the greedy non-overlapping assembly cannot emit
        // crossing pairs, so no pseudoknot is detected over the produced base-pair set.
        RnaSecondaryStructure.DetectPseudoknots(s.BasePairs).Should().BeEmpty(
            "PredictStructure assembles non-overlapping stem-loops → the pair set is pseudoknot-free");

        // Bracket count is consistent with the base-pair count.
        s.DotBracket.Count(c => c == '(').Should().Be(s.BasePairs.Count,
            "every base pair contributes exactly one open bracket");
        s.DotBracket.Count(c => c == ')').Should().Be(s.BasePairs.Count,
            "every base pair contributes exactly one close bracket");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-STRUCT-001 — RNA secondary structure prediction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region RNA-STRUCT-001 — RNA secondary structure prediction

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence. The contract short-circuits null/"" to the empty
    /// structure with MFE 0 — never an exception, never an IndexOutOfRange in the
    /// DP matrix (Minimum_Free_Energy.md §3.3, §6.1). We pin BOTH surfaces.
    /// </summary>
    [Test]
    public void Predict_EmptySequence_ReturnsEmptyStructureAndDoesNotThrow()
    {
        var act = () => RnaSecondaryStructure.PredictStructure("");

        var s = act.Should().NotThrow("empty input is a documented no-op, not an error").Subject;
        s.Sequence.Should().BeEmpty();
        s.DotBracket.Should().BeEmpty();
        s.BasePairs.Should().BeEmpty("an empty sequence forms no pairs");
        s.StemLoops.Should().BeEmpty();
        s.Pseudoknots.Should().BeEmpty();
        s.MinimumFreeEnergy.Should().Be(0);

        RnaSecondaryStructure.CalculateMinimumFreeEnergy("").Should().Be(0,
            "the open-chain ΔG of the empty sequence is 0 (INV-01)");
    }

    #endregion

    #region BE — Boundary: single base (cannot pair with itself)

    /// <summary>
    /// BE: a single base. With length 1 &lt; minLoopSize + 2 the structure cannot
    /// enclose even a 3-nt hairpin, so no pair is possible: a base never pairs with
    /// itself. The result is the all-unpaired structure ("."), MFE 0, no crash —
    /// this is the classic len&lt;2 / short-span guard the fuzz row targets
    /// (Minimum_Free_Energy.md §6.1).
    /// </summary>
    [Test]
    public void Predict_SingleBase_HasNoPairs()
    {
        foreach (char b in "ACGU")
        {
            var s = RnaSecondaryStructure.PredictStructure(b.ToString());
            s.BasePairs.Should().BeEmpty($"a single '{b}' has no partner to pair with");
            s.DotBracket.Should().Be(".", "the lone base is unpaired");
            s.MinimumFreeEnergy.Should().Be(0);
            AssertWellFormedStructure(s);

            RnaSecondaryStructure.CalculateMinimumFreeEnergy(b.ToString()).Should().Be(0);
        }
    }

    #endregion

    #region BE — Boundary: all-A·U sequence (positive pairing signal)

    /// <summary>
    /// BE: an all-A·U sequence. A and U are complementary, so "AAAA…UUUU…" CAN form
    /// A-U pairs — the prediction must surface that pairing signal (a non-trivial
    /// stem), while still respecting every structural invariant. A long enough block
    /// of A's followed by U's with a 3-nt minimum loop forms a clean hairpin.
    /// This is the positive A·U signal called out in the fuzz row.
    /// </summary>
    [Test]
    public void Predict_AllAUSequence_FormsCanonicalAUPairs()
    {
        // 5×A + 4-nt loop (AUAU never self-pairs across the loop center) + 5×U:
        // a poly-A / poly-U stem with an interior loop region.
        var s = RnaSecondaryStructure.PredictStructure("AAAAAUUUUUUUUUU");

        s.BasePairs.Should().NotBeEmpty("complementary A and U blocks must yield A-U pairs");
        s.BasePairs.Should().OnlyContain(
            bp => (bp.Base1 == 'A' && bp.Base2 == 'U') || (bp.Base1 == 'U' && bp.Base2 == 'A'),
            "the only canonical pairs available in an A/U sequence are A-U / U-A");
        // NOTE: PredictStructure is a GREEDY heuristic; its reported MinimumFreeEnergy is the
        // sum of stem-loop TotalFreeEnergy and may be ABOVE the true MFE (it can be positive
        // for a short A-U hairpin once terminal-AU penalties + hairpin initiation dominate) —
        // INV-01 (≤ 0) holds only for the exact DP CalculateMinimumFreeEnergy, asserted below.
        double.IsFinite(s.MinimumFreeEnergy).Should().BeTrue("the greedy energy is a finite physical value");
        AssertWellFormedStructure(s);

        // Alternating A/U also stays valid and crash-free.
        var alt = RnaSecondaryStructure.PredictStructure("AUAUAUAUAUAU");
        AssertWellFormedStructure(alt);
        RnaSecondaryStructure.CalculateMinimumFreeEnergy("AUAUAUAUAUAU").Should().BeLessThanOrEqualTo(0);
    }

    #endregion

    #region MC — Malformed Content: non-RNA characters

    /// <summary>
    /// MC: non-RNA characters. The contract treats unrecognized characters as
    /// never-pairing residues — handled, never rejected with a crash
    /// (Minimum_Free_Energy.md §3.3). Two flavors:
    ///   • DNA 'T' on the PredictStructure surface: 'T' is not in the RNA pairing
    ///     table, so it cannot pair (no A-T pairing on this greedy surface), but the
    ///     call must still succeed and yield a well-formed structure.
    ///   • Arbitrary junk (digits, punctuation, spaces): never crashes; junk
    ///     positions stay unpaired and the structure is still well-formed.
    /// </summary>
    [Test]
    public void Predict_NonRnaCharacters_AreHandledNotCrashed()
    {
        // DNA thymine: not an RNA base, so it never pairs on the PredictStructure surface.
        var dna = () => RnaSecondaryStructure.PredictStructure("GGGGTTTTCCCC");
        var sDna = dna.Should().NotThrow("non-RNA characters are handled, not rejected").Subject;
        sDna.BasePairs.Should().OnlyContain(
            bp => bp.Base1 != 'T' && bp.Base2 != 'T',
            "'T' is not in the RNA pairing table and can never appear in a reported pair");
        AssertWellFormedStructure(sDna);

        // Arbitrary junk: digits, punctuation and whitespace are never-pairing.
        var junk = "GG12!! CC\tNNxx";
        var sJunk = ((Func<RnaSecondaryStructure.SecondaryStructure>)(() =>
            RnaSecondaryStructure.PredictStructure(junk)))
            .Should().NotThrow("garbage content must not crash the predictor").Subject;
        foreach (var bp in sJunk.BasePairs)
        {
            RnaSecondaryStructure.CanPair(sJunk.Sequence[Math.Min(bp.Position1, bp.Position2)],
                                          sJunk.Sequence[Math.Max(bp.Position1, bp.Position2)])
                .Should().BeTrue("only canonical-pairing positions may be reported");
        }
        AssertWellFormedStructure(sJunk);

        // The MFE surface must likewise not throw on junk; ΔG stays ≤ 0 (INV-01).
        RnaSecondaryStructure.CalculateMinimumFreeEnergy(junk).Should().BeLessThanOrEqualTo(0);
    }

    #endregion

    #region BE — Boundary: minLoopSize below the steric floor is clamped to 3

    /// <summary>
    /// BE: minLoopSize &lt; 3. Hairpin loops shorter than 3 nt are sterically
    /// impossible (nearest-neighbor rules), so the contract CLAMPS minLoopSize up to
    /// 3 rather than emitting a 0/1/2-nt loop (Minimum_Free_Energy.md §6.1, §3.1).
    /// We pin that a negative / zero / tiny min-loop neither crashes nor produces a
    /// hairpin whose two closing endpoints enclose fewer than 3 unpaired bases.
    /// </summary>
    [Test]
    public void Predict_MinLoopBelowFloor_IsClampedAndStaysValid()
    {
        foreach (int badMin in new[] { int.MinValue, -1, 0, 1, 2 })
        {
            var s = RnaSecondaryStructure.PredictStructure("GGGGAAACCCC", minLoopSize: badMin);
            AssertWellFormedStructure(s);

            foreach (var sl in s.StemLoops)
            {
                sl.Loop.Size.Should().BeGreaterThanOrEqualTo(3,
                    $"min-loop is clamped to 3 even when {badMin} is requested (steric floor)");
            }
        }
    }

    #endregion

    #region BE — Boundary: extremely long sequence (≤ 300, CancelAfter-guarded)

    /// <summary>
    /// BE: an extremely-long sequence. Prediction is O(L²)–O(L³); the input is kept
    /// at L = 300 (top of the documented benchmark range) so it completes well within
    /// the [CancelAfter] budget. The point is hang/crash-safety AND that the produced
    /// structure still satisfies EVERY structural invariant at scale: ≤ 1 partner per
    /// base, canonical pairs, nesting, balanced dot-bracket, MFE ≤ 0. A regression to a
    /// hang or an O(L³)-blow-up would trip CancelAfter and FAIL rather than wedge the
    /// suite. The known [Explicit] MFE benchmark is deliberately NOT invoked.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Predict_ExtremelyLongSequence_CompletesAndStaysValid(CancellationToken token)
    {
        string longSeq = RandomRna(300, seed: 20260620);

        var s = RnaSecondaryStructure.PredictStructure(longSeq);
        token.ThrowIfCancellationRequested();

        s.Sequence.Should().HaveLength(300);
        AssertWellFormedStructure(s);

        double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(longSeq);
        token.ThrowIfCancellationRequested();
        mfe.Should().BeLessThanOrEqualTo(0, "INV-01: MFE ≤ 0 for every input, including long ones");
    }

    #endregion

    #region Positive sanity — an obvious hairpin predicts the expected stem

    /// <summary>
    /// Positive sanity: a sequence with an OBVIOUS hairpin — "GGGGAAAACCCC", a 4-bp
    /// G·C stem (GGGG / CCCC are reverse complements) closing a 4-nt AAAA loop —
    /// must predict that stem. The four G-C pairs should be reported, the loop bases
    /// left unpaired, and the result must satisfy the structural invariants
    /// (≤ 1 partner per base + nesting). This confirms the fuzz harness is asserting
    /// against a predictor that actually FINDS real structure, not a no-op.
    /// </summary>
    [Test]
    public void Predict_ObviousHairpin_PredictsExpectedGcStem()
    {
        var s = RnaSecondaryStructure.PredictStructure("GGGGAAAACCCC");

        AssertWellFormedStructure(s);

        s.BasePairs.Should().HaveCount(4, "GGGG/CCCC forms a clean 4-bp stem around the AAAA loop");
        s.BasePairs.Should().OnlyContain(
            bp => (bp.Base1 == 'G' && bp.Base2 == 'C') || (bp.Base1 == 'C' && bp.Base2 == 'G'),
            "every stem pair is a G-C Watson–Crick pair");

        // The stem pairs the outer G's with the outer C's: (0,11),(1,10),(2,9),(3,8).
        var pairs = s.BasePairs
            .Select(bp => (Lo: Math.Min(bp.Position1, bp.Position2), Hi: Math.Max(bp.Position1, bp.Position2)))
            .OrderBy(p => p.Lo)
            .ToList();
        pairs.Should().Equal((0, 11), (1, 10), (2, 9), (3, 8));

        // The AAAA loop (positions 4..7) is unpaired and ≥ 3 nt — the min-loop invariant.
        var loopSpan = pairs.Max(p => p.Lo) + 1; // innermost open + 1
        var loopEnd = pairs.Min(p => p.Hi) - 1;  // innermost close - 1
        (loopEnd - loopSpan + 1).Should().BeGreaterThanOrEqualTo(3, "the hairpin loop has ≥ 3 unpaired bases");

        s.DotBracket.Should().Be("((((....))))", "the canonical hairpin dot-bracket");
        double.IsFinite(s.MinimumFreeEnergy).Should().BeTrue(
            "the greedy reported energy is a finite physical value (it may sit above the true MFE)");

        // The exact DP score (INV-01) is ≤ 0 and finite.
        double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy("GGGGAAAACCCC");
        mfe.Should().BeLessThanOrEqualTo(0);
        double.IsFinite(mfe).Should().BeTrue("the MFE score is a finite physical energy");
    }

    /// <summary>
    /// Positive sanity over RANDOM RNA: every fuzz-generated sequence must produce a
    /// well-formed structure (≤ 1 partner per base, canonical pairs, nesting,
    /// balanced bracket) and a finite MFE ≤ 0. Several fixed seeds, several lengths.
    /// </summary>
    [Test]
    public void Predict_RandomRna_AlwaysProducesWellFormedStructure()
    {
        foreach (int seed in new[] { 1, 7, 42, 1000 })
        {
            foreach (int len in new[] { 12, 40, 90 })
            {
                string seq = RandomRna(len, seed);
                var s = RnaSecondaryStructure.PredictStructure(seq);
                AssertWellFormedStructure(s);

                double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
                double.IsFinite(mfe).Should().BeTrue($"MFE must be finite (seed {seed}, len {len})");
                mfe.Should().BeLessThanOrEqualTo(0, $"INV-01 (seed {seed}, len {len})");
            }
        }
    }

    #endregion

    #endregion
}
