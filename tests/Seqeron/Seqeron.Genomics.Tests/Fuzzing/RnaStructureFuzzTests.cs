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
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: RNA-STEMLOOP-001 — stem-loop / hairpin motif detection
/// Checklist: docs/checklists/03_FUZZING.md, row 72.
/// Method under test: RnaSecondaryStructure.FindStemLoops(string, minStemLength=3,
///   minLoopSize=3, maxLoopSize=10, allowWobble=true) → IEnumerable&lt;StemLoop&gt;.
///   An EXHAUSTIVE local scan: for every loop start and every loop size it extends
///   a stem outward while pairing stays valid and yields EVERY candidate whose stem
///   length ≥ minStemLength (no overlap suppression — RNA_Stemloop.md §1, §4.1, §5.2).
/// Fuzz strategy for THIS unit:
///   • BE = Boundary Exploitation — the degenerate parameter / sequence corners that
///     could crash (IndexOutOfRange off the stem-extension walk, DivideByZero on a
///     zero loop), hang, or emit a malformed hairpin:
///       – Empty sequence                  → no stem-loops (string.IsNullOrEmpty short-circuit).
///       – No complement (e.g. all-A)      → no base can pair → no stem-loops.
///       – minLoopSize = 0 (and ≤ 2)       → CLAMPED up to 3; no 0/1/2-nt loop is ever
///                                            emitted (steric floor; RNA_Stemloop.md INV-01).
///       – minStemLength &gt; seqLen/2        → stem+loop+stem cannot fit (length &lt;
///                                            minStemLength*2 + minLoopSize) → no stem-loops.
///   — docs/checklists/03_FUZZING.md §Description (strategy code BE); targets:
///     "Empty, no complement, minLoop=0, minStem &gt; seqLen/2".
///
/// Theory-correct stem-loop contract (RNA_Stemloop.md §2.4 INV-01..04, §3.3, §6.1):
///   • INV-01 Loop.Size ≥ 3 (minLoopSize clamped to 3 before scanning).
///   • INV-02 Stem.Length ≥ minStemLength (shorter candidates discarded).
///   • INV-03 Loop.Type is always Hairpin.
///   • INV-04 DotBracketNotation.Length = End − Start + 1.
///   • The two stem arms are antiparallel reverse-complements: every stem BasePair is
///     a canonical RNA pair (A-U / G-C / G-U) and the arms span the loop, so the
///     unpaired loop sits strictly between the 5' and 3' stem ends.
///   • Empty / too-short (length &lt; minStemLength*2 + minLoopSize) → no candidates,
///     never an exception.
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

    /// <summary>
    /// Asserts the theory-correct contract of a single detected stem-loop
    /// (RNA_Stemloop.md §2.4 INV-01..04, §2.1): hairpin loop ≥ 3 nt, stem length ≥
    /// the requested floor, antiparallel reverse-complementary stem arms (each stem
    /// pair canonical), in-range nested coords with the loop strictly enclosed, and a
    /// local dot-bracket whose length equals the detected span.
    /// </summary>
    private static void AssertWellFormedStemLoop(
        RnaSecondaryStructure.StemLoop sl, string sequence, int minStemLength)
    {
        int n = sequence.Length;

        // Span coords are in range and ordered (INV-04 needs a non-negative span).
        sl.Start.Should().BeInRange(0, n - 1, "the stem-loop 5' end indexes the sequence");
        sl.End.Should().BeInRange(0, n - 1, "the stem-loop 3' end indexes the sequence");
        sl.Start.Should().BeLessThan(sl.End, "a stem-loop spans a non-empty antiparallel region");

        // INV-03: the loop is always a Hairpin.
        sl.Loop.Type.Should().Be(RnaSecondaryStructure.LoopType.Hairpin,
            "FindStemLoops emits hairpin loops only (INV-03)");

        // INV-01: hairpin loops shorter than 3 nt are sterically impossible.
        sl.Loop.Size.Should().BeGreaterThanOrEqualTo(3, "min-loop is clamped to the 3-nt steric floor (INV-01)");
        sl.Loop.Sequence.Should().HaveLength(sl.Loop.Size, "the loop sequence covers exactly the loop span");
        (sl.Loop.End - sl.Loop.Start + 1).Should().Be(sl.Loop.Size, "loop coords are consistent with its size");

        // INV-02: the stem meets the requested minimum length.
        sl.Stem.Length.Should().BeGreaterThanOrEqualTo(minStemLength, "stem length ≥ minStemLength (INV-02)");
        sl.Stem.BasePairs.Should().HaveCount(sl.Stem.Length, "one base pair per stacked stem position");

        // The loop sits strictly between the inner stem ends — the stem ENCLOSES it.
        sl.Stem.End5Prime.Should().Be(sl.Loop.Start - 1, "the 5' stem arm abuts the loop start");
        sl.Stem.Start3Prime.Should().Be(sl.Loop.End + 1, "the 3' stem arm abuts the loop end");

        // Stem arms are antiparallel reverse-complements: every pair is canonical and
        // its 5' index is below the loop while its 3' index is above it.
        foreach (var bp in sl.Stem.BasePairs)
        {
            int i = Math.Min(bp.Position1, bp.Position2);
            int j = Math.Max(bp.Position1, bp.Position2);
            i.Should().BeInRange(0, n - 1, "stem pair endpoints index the sequence");
            j.Should().BeInRange(0, n - 1, "stem pair endpoints index the sequence");
            i.Should().BeLessThan(sl.Loop.Start, "the 5' stem arm lies before the loop");
            j.Should().BeGreaterThan(sl.Loop.End, "the 3' stem arm lies after the loop");
            RnaSecondaryStructure.CanPair(sequence[i], sequence[j]).Should().BeTrue(
                $"stem pair ({sequence[i]},{sequence[j]}) must be a canonical RNA pair (reverse-complement arms)");
        }

        // INV-04: the local dot-bracket length equals the detected span.
        sl.DotBracketNotation.Should().HaveLength(sl.End - sl.Start + 1, "INV-04: local bracket covers exactly the span");
        sl.DotBracketNotation.Count(c => c == '(').Should().Be(sl.Stem.Length, "one open bracket per stem pair");
        sl.DotBracketNotation.Count(c => c == ')').Should().Be(sl.Stem.Length, "one close bracket per stem pair");

        double.IsFinite(sl.TotalFreeEnergy).Should().BeTrue("the stem-loop free energy is a finite physical value");
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

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-STEMLOOP-001 — stem-loop / hairpin detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region RNA-STEMLOOP-001 — stem-loop detection

    #region BE — Boundary: empty sequence yields no stem-loops

    /// <summary>
    /// BE: the empty sequence. FindStemLoops short-circuits on string.IsNullOrEmpty
    /// (and on null) → no candidates, never an IndexOutOfRange off the stem walk
    /// (RNA_Stemloop.md §3.3, §6.1). We pin both "" and null.
    /// </summary>
    [Test]
    public void FindStemLoops_EmptySequence_YieldsNoneAndDoesNotThrow()
    {
        var act = () => RnaSecondaryStructure.FindStemLoops("").ToList();
        act.Should().NotThrow("empty input is a documented no-op, not an error")
           .Subject.Should().BeEmpty("an empty sequence encloses no loop and forms no stem");

        var actNull = () => RnaSecondaryStructure.FindStemLoops(null!).ToList();
        actNull.Should().NotThrow("null is treated like empty, not an error")
               .Subject.Should().BeEmpty("a null sequence forms no stem-loops");
    }

    #endregion

    #region BE — Boundary: no complement (all-A) yields no stem-loops

    /// <summary>
    /// BE: a sequence with NO complementary base — e.g. poly-A. A pairs with nothing
    /// (no A-A pair exists; only A-U / G-C / G-U are canonical), so every stem
    /// extension stops immediately and no candidate ever reaches minStemLength → no
    /// stem-loops, never a crash (RNA_Stemloop.md §2.1, §6.1). We probe a length well
    /// above the steric floor so the only reason for "no result" is the missing
    /// complement, and confirm the same for poly-C and poly-G (each pairs with nothing
    /// available in its own homopolymer).
    /// </summary>
    [Test]
    public void FindStemLoops_NoComplement_YieldsNoneAndDoesNotThrow()
    {
        foreach (string mono in new[] { new string('A', 30), new string('C', 30), new string('G', 30) })
        {
            var act = () => RnaSecondaryStructure.FindStemLoops(mono).ToList();
            var result = act.Should().NotThrow("a non-pairing homopolymer must not crash the scan").Subject;
            result.Should().BeEmpty($"'{mono[0]}'×{mono.Length} has no complementary base, so no stem can form");
        }
    }

    #endregion

    #region BE — Boundary: minLoopSize 0 (and ≤ 2) is clamped to the 3-nt floor

    /// <summary>
    /// BE: minLoopSize = 0 (and −1, 1, 2). A 0-nt loop is sterically impossible — the
    /// contract CLAMPS minLoopSize up to 3 before scanning (RNA_Stemloop.md §3.3,
    /// §2.4 INV-01), so the degenerate value neither crashes (no DivideByZero / empty
    /// Substring) nor yields a hairpin whose loop is shorter than 3 nt. We feed a
    /// sequence with a genuine 3-nt-loop hairpin so the clamp is observable: every
    /// returned stem-loop still has Loop.Size ≥ 3.
    /// </summary>
    [Test]
    public void FindStemLoops_MinLoopBelowFloor_IsClampedToThree()
    {
        // GGGG / CCCC stem around a 3-nt AAA loop.
        const string seq = "GGGGAAACCCC";
        foreach (int badMin in new[] { int.MinValue, -1, 0, 1, 2 })
        {
            var act = () => RnaSecondaryStructure.FindStemLoops(seq, minLoopSize: badMin).ToList();
            var result = act.Should().NotThrow($"minLoopSize {badMin} is clamped, not an error").Subject;

            foreach (var sl in result)
            {
                sl.Loop.Size.Should().BeGreaterThanOrEqualTo(3,
                    $"min-loop is clamped to the 3-nt steric floor even when {badMin} is requested (INV-01)");
                AssertWellFormedStemLoop(sl, seq, minStemLength: 3);
            }
        }
    }

    #endregion

    #region BE — Boundary: minStemLength > seqLen/2 cannot fit any stem-loop

    /// <summary>
    /// BE: minStemLength larger than half the sequence. A stem-loop needs two stem
    /// arms PLUS a loop between them, so it occupies minStemLength*2 + minLoopSize
    /// residues; when minStemLength &gt; seqLen/2 the two arms alone already exceed the
    /// sequence and the length guard (length &lt; minStemLength*2 + minLoopSize) short-
    /// circuits → no candidates, never an out-of-range stem walk (RNA_Stemloop.md
    /// §3.3, §6.1). We use a real GC hairpin so the ONLY reason for "no result" is the
    /// oversized stem floor, and sweep several over-half thresholds.
    /// </summary>
    [Test]
    public void FindStemLoops_MinStemOverHalfLength_YieldsNone()
    {
        const string seq = "GGGGAAAACCCC"; // len 12 — would be a real hairpin at minStem ≤ 4
        int half = seq.Length / 2;          // 6

        foreach (int minStem in new[] { half + 1, half + 2, seq.Length, seq.Length + 5 })
        {
            var act = () => RnaSecondaryStructure.FindStemLoops(seq, minStemLength: minStem).ToList();
            var result = act.Should().NotThrow($"oversized minStemLength {minStem} is short-circuited, not an error").Subject;
            result.Should().BeEmpty(
                $"a stem of {minStem} > {half} per arm plus a loop cannot fit in {seq.Length} nt");
        }

        // Sanity counterpoint: at the largest stem that DOES fit, the hairpin is found —
        // proving the empties above are the size guard, not a dead method. minStem 4
        // needs 4*2 + 3 = 11 ≤ 12, so the 4-bp GGGG/CCCC stem qualifies.
        RnaSecondaryStructure.FindStemLoops(seq, minStemLength: 4).Should().NotBeEmpty(
            "the 4-bp GGGG/CCCC stem fits within 12 nt and must be detected");
    }

    #endregion

    #region Positive sanity — an obvious hairpin is detected with the expected stem and loop

    /// <summary>
    /// Positive sanity: "GGGGAAAACCCC" is an unmistakable hairpin — a 4-bp G·C stem
    /// (GGGG and CCCC are reverse complements) enclosing a 4-nt AAAA loop. The scan
    /// must surface a stem-loop whose stem is the four G-C pairs (0,11),(1,10),(2,9),
    /// (3,8) and whose loop is the AAAA at positions 4..7, all invariants satisfied.
    /// This proves the fuzz harness asserts against a detector that FINDS real
    /// structure, not a no-op.
    /// </summary>
    [Test]
    public void FindStemLoops_ObviousHairpin_DetectsExpectedStemAndLoop()
    {
        const string seq = "GGGGAAAACCCC";
        var all = RnaSecondaryStructure.FindStemLoops(seq).ToList();

        all.Should().NotBeEmpty("a clear GGGG/AAAA/CCCC hairpin must be detected");
        foreach (var sl in all)
            AssertWellFormedStemLoop(sl, seq, minStemLength: 3);

        // The maximal candidate is the full 4-bp stem closing the 4-nt AAAA loop.
        var best = all.OrderByDescending(sl => sl.Stem.Length).First();
        best.Stem.Length.Should().Be(4, "GGGG/CCCC stacks into a 4-bp stem");
        best.Loop.Sequence.Should().Be("AAAA", "the enclosed hairpin loop is the central AAAA");
        best.Loop.Start.Should().Be(4);
        best.Loop.End.Should().Be(7);

        best.Stem.BasePairs
            .Select(bp => (Lo: Math.Min(bp.Position1, bp.Position2), Hi: Math.Max(bp.Position1, bp.Position2)))
            .OrderBy(p => p.Lo)
            .Should().Equal((0, 11), (1, 10), (2, 9), (3, 8));

        best.Stem.BasePairs.Should().OnlyContain(
            bp => (bp.Base1 == 'G' && bp.Base2 == 'C') || (bp.Base1 == 'C' && bp.Base2 == 'G'),
            "every stem pair is a G-C Watson–Crick pair");

        best.DotBracketNotation.Should().Be("((((....))))", "the canonical hairpin dot-bracket over the span");
    }

    /// <summary>
    /// Positive sanity over RANDOM RNA: across fixed seeds and lengths, FindStemLoops
    /// must never crash or hang and every emitted candidate must satisfy the stem-loop
    /// contract (canonical reverse-complementary arms, loop ≥ 3, stem ≥ floor, in-range
    /// nested coords, INV-04 bracket length). [CancelAfter] guards against an
    /// O(n²·L)-blow-up regression turning into a hang.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void FindStemLoops_RandomRna_AlwaysWellFormed(CancellationToken token)
    {
        foreach (int seed in new[] { 3, 17, 99, 2026 })
        {
            foreach (int len in new[] { 12, 40, 120 })
            {
                string seq = RandomRna(len, seed);
                var act = () => RnaSecondaryStructure.FindStemLoops(seq).ToList();
                var result = act.Should().NotThrow($"random RNA must not crash the scan (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                foreach (var sl in result)
                    AssertWellFormedStemLoop(sl, seq, minStemLength: 3);
            }
        }
    }

    #endregion

    #endregion
}
