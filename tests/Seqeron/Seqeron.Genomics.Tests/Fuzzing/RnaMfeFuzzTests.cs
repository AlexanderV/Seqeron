using System;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for RNA-MFE-001 — the MINIMUM FREE ENERGY (MFE) folding DP.
/// Unit under test: RnaSecondaryStructure.CalculateMinimumFreeEnergy(string, int minLoopSize)
/// in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs — a Zuker–
/// Stiegler O(n³) dynamic program over the loop decomposition, scored with the Turner 2004
/// nearest-neighbor parameters, returning the optimal ΔG°37 in kcal/mol (rounded to 2 dp).
/// Checklist: docs/checklists/03_FUZZING.md, row 152 (strategy code BE).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Degenerate, boundary and out-of-domain inputs must NEVER fail in an undisciplined
/// way: no IndexOutOfRange off the n×n DP table on length-0/1 spans, no DivideByZero,
/// no NaN/Inf energy, no false base pairs in a sequence that cannot pair, and no
/// hang / super-cubic blow-up on a long homopolymer. Every input must resolve to a
/// well-defined, theory-correct score OR a documented validation outcome.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Scope vs the sibling RNA fuzz rows
/// ───────────────────────────────────────────────────────────────────────────
/// Row 71 (RNA-STRUCT-001) covers the greedy PredictStructure surface; row 73
/// (RNA-ENERGY-001) covers the per-loop / Boltzmann energy surfaces. THIS row is the
/// MFE DP FOLD itself — CalculateMinimumFreeEnergy — focused on the BE corners that
/// stress the DP table: empty, single base, all-A (no pairs), homopolymer. See
/// RnaStructureFuzzTests.cs for the structure / stem-loop / loop-energy rows.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy code BE = Boundary Exploitation (граничні значення: 0, -1, empty)
/// — docs/checklists/03_FUZZING.md §Description. Row-152 targets:
///   "empty, single base, all-A (no pairs), homopolymer".
/// ───────────────────────────────────────────────────────────────────────────
///
/// The DP contract under test (docs/algorithms/RnaStructure/Minimum_Free_Energy.md):
///   • §3.3 / §6.1: null / "" → 0; length &lt; minLoopSize + 2 → 0 (cannot enclose a
///     3-nt hairpin); minLoopSize &lt; 3 → clamped up to 3; only A/C/G/U pair (WC + G-U
///     wobble); unrecognized characters never pair; no exception for unfoldable input.
///   • §6.1 edge case: a homopolymer (e.g. "AAAAAAAA") has NO possible base pair, so the
///     only structure is the open chain with ΔG = 0 → MFE = 0 (no false pairs).
///   • INV-01 (§2.4): MFE ≤ 0 for every input — the open-chain (ΔG = 0) structure is
///     always in the search space, so the optimum is never positive.
///   • INV-02 (§2.4): MFE is non-increasing under suffix extension —
///     MFE(prefix) ≥ MFE(prefix + suffix) — extension only adds folding options.
///   • INV-03 (§2.4): the score is deterministic — repeated evaluation is identical.
///   • Worked example (§7.1, NNDB Turner 2004 Hairpin Example 1): "CACAAAAAAAUGUG"
///     folds to −1.41 (NNDB tabulated −1.4); Example 2 "CACAGAAAGUGUG" → −1.91.
///
/// Complexity / hang-safety (§4.3): the fold is O(n³) with MAXLOOP = 30. The long
/// homopolymer / random "extreme length" targets are kept MODEST (n ≤ 300, the top of
/// the documented benchmark range, ~214 ms at n=300) and guarded with [CancelAfter] so
/// a regression that introduced a hang or a non-terminating DP would FAIL, not wedge the
/// suite. The known [Explicit] RnaSecondaryStructure_MFE_Benchmark is NOT invoked here.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaMfeFuzzTests
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

    /// <summary>
    /// The well-formed-result helper for the MFE SCORE surface: an MFE is always a finite
    /// physical energy and, by INV-01, never positive. (The DP returns a single double; the
    /// "structure well-formedness" contract belongs to PredictStructure / RNA-STRUCT-001.)
    /// </summary>
    private static void AssertWellFormedMfe(double mfe, string because)
    {
        double.IsNaN(mfe).Should().BeFalse($"MFE must never be NaN ({because})");
        double.IsInfinity(mfe).Should().BeFalse($"MFE must never be ±∞ ({because})");
        mfe.Should().BeLessThanOrEqualTo(0, $"INV-01: MFE ≤ 0 for every input ({because})");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-MFE-001 — minimum free energy folding (Zuker DP) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region RNA-MFE-001 — minimum free energy folding

    #region BE — Boundary: empty / null sequence → 0, no DP-table crash

    /// <summary>
    /// BE — "empty": the empty (and null) sequence. The contract short-circuits null/""
    /// to ΔG = 0 BEFORE allocating any DP buffer, so there is no n×n table to index and
    /// no IndexOutOfRange / DivideByZero can occur (Minimum_Free_Energy.md §3.3, §6.1).
    /// We pin BOTH the value (exactly 0) and the no-throw guarantee, and that the default
    /// vs clamped minLoopSize make no difference on this empty boundary.
    /// </summary>
    [Test]
    public void Mfe_EmptyOrNull_ReturnsZeroAndDoesNotThrow()
    {
        var actEmpty = () => RnaSecondaryStructure.CalculateMinimumFreeEnergy("");
        double e = actEmpty.Should().NotThrow("empty input is a documented no-op, not an error").Subject;
        e.Should().Be(0, "the open-chain ΔG of the empty sequence is exactly 0 (INV-01 / §6.1)");

        var actNull = () => RnaSecondaryStructure.CalculateMinimumFreeEnergy(null!);
        actNull.Should().NotThrow("null is treated like empty, never an exception")
               .Subject.Should().Be(0);

        // A pathological minLoopSize on the empty boundary still short-circuits to 0.
        RnaSecondaryStructure.CalculateMinimumFreeEnergy("", minLoopSize: 0).Should().Be(0);
        RnaSecondaryStructure.CalculateMinimumFreeEnergy("", minLoopSize: int.MinValue).Should().Be(0);
    }

    #endregion

    #region BE — Boundary: single base → no pair possible, 0, no 1×1-table crash

    /// <summary>
    /// BE — "single base": a length-1 sequence. With 1 &lt; minLoopSize + 2 the DP cannot
    /// enclose even a 3-nt hairpin and a base can never pair with itself, so the MFE is
    /// exactly 0 with NO false pair. The danger this targets is an IndexOutOfRange on a
    /// 1×1 DP table or an off-by-one in the short-span guard; the contract short-circuits
    /// length &lt; minLoopSize + 2 to 0 (Minimum_Free_Energy.md §3.3, §6.1). All four bases
    /// are pinned so the result does not depend on which single residue is supplied.
    /// </summary>
    [Test]
    public void Mfe_SingleBase_ReturnsZeroAndDoesNotThrow()
    {
        foreach (char b in "ACGU")
        {
            var act = () => RnaSecondaryStructure.CalculateMinimumFreeEnergy(b.ToString());
            double mfe = act.Should().NotThrow($"a single '{b}' is a documented no-op, not an error").Subject;
            mfe.Should().Be(0, $"a single '{b}' has no partner — the only structure is the open chain (ΔG = 0)");
        }

        // The lengths strictly below the steric floor (minLoopSize + 2 = 5) likewise fold to 0:
        // no pair (i,j) can enclose a 3-nt hairpin, so the DP never records a favorable pairing.
        foreach (string tooShort in new[] { "GC", "GCG", "GCGC" })
        {
            RnaSecondaryStructure.CalculateMinimumFreeEnergy(tooShort).Should().Be(0,
                $"'{tooShort}' is shorter than minLoopSize + 2 and cannot enclose a 3-nt hairpin → MFE 0");
        }
    }

    #endregion

    #region BE — Boundary: all-A (no pairs) → exactly 0, no FALSE pairs

    /// <summary>
    /// BE — "all-A (no pairs)": poly-A has no complementary base anywhere (A pairs with
    /// nothing — only A-U / G-C / G-U are canonical), so the ONLY admissible structure is
    /// the open chain and the MFE is EXACTLY 0 (Minimum_Free_Energy.md §6.1 homopolymer
    /// edge case). The failure this catches is a DP that records a FALSE pair (a negative
    /// energy where none is physically possible) — the score must be 0, not below 0. We
    /// sweep several lengths well above the steric floor so the only reason for "0" is the
    /// missing complement, not the short-span short-circuit.
    /// </summary>
    [Test]
    public void Mfe_AllA_NoPairsPossible_ReturnsExactlyZero()
    {
        foreach (int len in new[] { 5, 8, 20, 50 })
        {
            string polyA = new string('A', len);
            double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(polyA);
            mfe.Should().Be(0,
                $"A×{len} has no complementary base — the open chain (ΔG = 0) is the only fold; " +
                "a negative score would mean a FALSE pair was recorded");
        }
    }

    #endregion

    #region BE — Boundary: homopolymers (every base) → no self-pair, exactly 0

    /// <summary>
    /// BE — "homopolymer": no single nucleotide pairs with itself — there is no A-A, C-C,
    /// G-G or U-U canonical pair — so EVERY homopolymer folds to the open chain with MFE
    /// exactly 0. (Wobble G-U needs two DIFFERENT residues, so a pure G or pure U run also
    /// has no pair.) This generalizes the all-A target to all four bases and pins that the
    /// no-pair result holds regardless of which residue dominates and at lengths far above
    /// the steric floor — never a false pair, never a crash. A modestly long homopolymer
    /// (n = 200) is included under [CancelAfter] to prove the DP TERMINATES (no hang /
    /// non-termination on a degenerate all-same input) and still returns exactly 0.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Mfe_Homopolymers_HaveNoSelfPair_ReturnZero(CancellationToken token)
    {
        foreach (char b in "ACGU")
        {
            foreach (int len in new[] { 6, 15, 40 })
            {
                string mono = new string(b, len);
                var act = () => RnaSecondaryStructure.CalculateMinimumFreeEnergy(mono);
                double mfe = act.Should().NotThrow($"a non-pairing '{b}'×{len} must not crash the DP").Subject;
                mfe.Should().Be(0,
                    $"'{b}'×{len} has no self-pair (no A-A/C-C/G-G/U-U pair) → MFE 0, no false pairing");
            }
        }

        // A modestly long homopolymer must TERMINATE (O(n³) DP, no hang) and still be 0.
        string longMono = new string('G', 200);
        double longMfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(longMono);
        token.ThrowIfCancellationRequested();
        longMfe.Should().Be(0, "a long pure-G run has no pair — the DP terminates with the open-chain 0");
    }

    #endregion

    #region BE — Boundary: minLoopSize below the 3-nt steric floor is clamped

    /// <summary>
    /// BE: a minLoopSize below the 3-nt steric floor (including int.MinValue and 0).
    /// Hairpin loops shorter than 3 nt are sterically impossible, so the DP CLAMPS
    /// minLoopSize up to 3 before folding (Minimum_Free_Energy.md §3.1, §6.1) rather than
    /// admitting a 0/1/2-nt loop. The degenerate value must neither crash nor over-stabilize:
    /// the score with any sub-floor minLoopSize must EQUAL the score at the clamped floor of 3.
    /// </summary>
    [Test]
    public void Mfe_MinLoopBelowFloor_IsClampedToThree()
    {
        const string seq = "GGGGAAACCCC"; // a genuine hairpin with a 3-nt AAA loop
        double atFloor = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq, minLoopSize: 3);

        foreach (int badMin in new[] { int.MinValue, -1, 0, 1, 2 })
        {
            var act = () => RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq, minLoopSize: badMin);
            double mfe = act.Should().NotThrow($"minLoopSize {badMin} is clamped, not an error").Subject;
            mfe.Should().Be(atFloor,
                $"minLoopSize {badMin} clamps to the 3-nt floor — identical fold to minLoopSize = 3");
            AssertWellFormedMfe(mfe, $"clamped minLoopSize {badMin}");
        }
    }

    #endregion

    #region INV-01 / INV-02 / INV-03 — never positive, monotone, deterministic

    /// <summary>
    /// INV-01: MFE ≤ 0 for EVERY input. The open-chain (ΔG = 0) structure is always in the
    /// search space, so the optimum can never be positive (Minimum_Free_Energy.md §2.4).
    /// We assert it over a spread of deterministic and fuzz-generated sequences — including
    /// foldable, non-pairing, junk and mixed-case inputs — and confirm the score is finite.
    /// </summary>
    [Test]
    public void Mfe_NeverPositive_OverManyInputs()
    {
        foreach (string seq in new[]
        {
            "GGGGAAAACCCC", "CACAAAAAAAUGUG", "CACAGAAAGUGUG",
            "AUAUAUAUAUAU", "AAAAAAAA", "CCCCCCCC", "gggaaaccc",
            "GG12!! CC\tNNxx", "GGGGUUUUCCCC"
        })
        {
            double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
            AssertWellFormedMfe(mfe, $"input \"{seq}\"");
        }

        foreach (int seed in new[] { 1, 7, 42, 1000 })
            foreach (int len in new[] { 8, 25, 60 })
            {
                string seq = RandomRna(len, seed);
                AssertWellFormedMfe(RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq),
                    $"random seed {seed} len {len}");
            }
    }

    /// <summary>
    /// INV-02: MFE is non-increasing under suffix extension — MFE(prefix) ≥ MFE(prefix +
    /// suffix), because extending the sequence only ADDS folding options and the optimum
    /// over a superset cannot be larger (Minimum_Free_Energy.md §2.4). We grow a hairpin
    /// sequence one suffix at a time and assert the score never increases.
    /// </summary>
    [Test]
    public void Mfe_NonIncreasingUnderSuffixExtension()
    {
        string baseSeq = RandomRna(20, seed: 20260620);
        double prev = RnaSecondaryStructure.CalculateMinimumFreeEnergy(baseSeq);

        var rng = new Random(20260621);
        const string bases = "ACGU";
        string grown = baseSeq;
        for (int step = 0; step < 12; step++)
        {
            grown += bases[rng.Next(bases.Length)];
            double now = RnaSecondaryStructure.CalculateMinimumFreeEnergy(grown);
            now.Should().BeLessThanOrEqualTo(prev + 1e-9,
                $"INV-02: extending \"{baseSeq}\" (step {step}) only adds folding options → MFE cannot rise");
            prev = now;
        }
    }

    /// <summary>
    /// INV-03: the MFE score is deterministic — the DP has no randomness, so repeated
    /// evaluation of the same sequence yields the identical value (Minimum_Free_Energy.md
    /// §2.4). We re-evaluate several sequences (including the long homopolymer and a long
    /// random fold) and assert bit-for-bit equality across repeated calls.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Mfe_IsDeterministic_AcrossRepeatedCalls(CancellationToken token)
    {
        foreach (string seq in new[]
        {
            "CACAAAAAAAUGUG", "GGGGAAAACCCC", new string('A', 100), RandomRna(120, seed: 555)
        })
        {
            double first = RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq);
            token.ThrowIfCancellationRequested();
            for (int rep = 0; rep < 3; rep++)
                RnaSecondaryStructure.CalculateMinimumFreeEnergy(seq).Should().Be(first,
                    "INV-03: the DP is deterministic — identical input yields the identical score");
        }
    }

    #endregion

    #region BE — extreme length, CancelAfter-guarded: terminates, finite, ≤ 0

    /// <summary>
    /// BE — "extremely long" companion to the homopolymer target: a random length-300
    /// fold (top of the documented benchmark range, ~214 ms) must TERMINATE within the
    /// [CancelAfter] budget — proving the O(n³) DP does not hang or blow up — and return a
    /// finite MFE ≤ 0 (INV-01). A regression to a hang trips CancelAfter and FAILS rather
    /// than wedging the suite.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void Mfe_ExtremelyLongSequence_TerminatesFiniteAndNonPositive(CancellationToken token)
    {
        string longSeq = RandomRna(300, seed: 20260620);
        double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(longSeq);
        token.ThrowIfCancellationRequested();
        AssertWellFormedMfe(mfe, "random length-300 fold");
    }

    #endregion

    #region Positive sanity — a foldable hairpin yields the documented NEGATIVE MFE

    /// <summary>
    /// Positive sanity (NNDB Turner 2004 Hairpin worked examples, Minimum_Free_Energy.md
    /// §7.1 / RNA-MFE-001 M1·M2): the fold must surface real, NEGATIVE stability with the
    /// EXACT documented values — proving the fuzz harness asserts against a DP that actually
    /// folds, not a no-op that always returns 0.
    ///   • "CACAAAAAAAUGUG" → −1.41 (NNDB Example 1; tabulated −1.4).
    ///   • "CACAGAAAGUGUG"  → −1.91 (NNDB Example 2; GG first-mismatch bonus).
    /// We also pin the no-pair counterpoint: a non-pairing sequence of the SAME length is
    /// exactly 0 — confirming the negative scores are genuine folding, not a constant.
    /// </summary>
    [Test]
    public void Mfe_FoldableHairpins_ReturnDocumentedNegativeEnergy()
    {
        double ex1 = RnaSecondaryStructure.CalculateMinimumFreeEnergy("CACAAAAAAAUGUG");
        ex1.Should().BeApproximately(-1.41, 1e-2, "NNDB Turner 2004 Hairpin Example 1 (tabulated −1.4)");

        double ex2 = RnaSecondaryStructure.CalculateMinimumFreeEnergy("CACAGAAAGUGUG");
        ex2.Should().BeApproximately(-1.91, 1e-2, "NNDB Turner 2004 Hairpin Example 2 (GG first mismatch, −1.9)");

        // No-pair counterpoint at the same length as Example 1 → exactly 0 (no false fold).
        RnaSecondaryStructure.CalculateMinimumFreeEnergy(new string('A', 14)).Should().Be(0,
            "an unfoldable run of the same length scores exactly 0 — the negatives above are real folding");
    }

    /// <summary>
    /// Positive sanity — a clean G·C hairpin "GGGGAAAACCCC" (a 4-bp GGGG/CCCC stem over a
    /// 4-nt AAAA loop) is a strong fold, so its DP MFE must be finite and clearly negative
    /// (G-C stacking dominates). This is the foldable-vs-no-pair contrast called out by the
    /// row: a foldable sequence yields a negative energy; a no-pair one yields 0.
    /// </summary>
    [Test]
    public void Mfe_ObviousGcHairpin_HasFiniteNegativeEnergy()
    {
        double mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy("GGGGAAAACCCC");
        double.IsFinite(mfe).Should().BeTrue("a foldable hairpin yields a finite physical energy");
        mfe.Should().BeLessThan(0, "a 4-bp G·C stem is a stable fold → strictly negative MFE");
    }

    #endregion

    #endregion
}
