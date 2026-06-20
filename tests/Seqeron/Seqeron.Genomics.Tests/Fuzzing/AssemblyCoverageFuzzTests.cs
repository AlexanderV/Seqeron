using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Assembly area — Coverage (Depth) Calculation (ASSEMBLY-COVER-001),
/// the per-base sequencing-depth array
/// <see cref="SequenceAssembler.CalculateCoverage(string, IReadOnlyList{string}, int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop
/// (the O(r·n·m) placement+increment must always terminate), no state corruption,
/// no nonsense output (a depth array whose length differs from the reference, a
/// negative depth, a non-finite derived mean), and no *unhandled* runtime
/// exception — in particular NO DivideByZero when a CALLER derives the mean from a
/// zero-length reference or an empty/no-mapped-read input, and NO IndexOutOfRange
/// when a read would extend past the reference end. Every input must resolve to
/// EITHER a well-defined, theory-correct depth array OR a *documented, intentional*
/// validation exception (ArgumentNullException for a null reference / null reads —
/// contract §3.3, §6.1). A raw runtime exception, a hang, a wrong-length array, a
/// negative count, or an IndexOutOfRange on an over-long read is a bug, not a
/// passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-COVER-001 — Coverage (Depth) Calculation
/// Checklist: docs/checklists/03_FUZZING.md, row 142.
/// Algorithm doc: docs/algorithms/Assembly/Coverage_Calculation.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row:
///          – NO READS: empty reads list → an all-zero depth array of length
///            reference.Length (§3.3, §6.1, INV-05). The CALLER-derived mean
///            (Σ depth / G) is 0/G = 0 — NO DivideByZero, and the array itself is
///            all zeros (no read contributes any depth).
///          – ZERO-LENGTH REF: empty reference → an EMPTY (length-0) depth array
///            (§3.3, §6.1). A caller deriving the mean as Σ depth / G hits the
///            length-0 denominator; the unit must still return cleanly (empty
///            array, no crash) and a caller guarding G == 0 sees 0/0 — the test
///            asserts the documented empty-array result and that no exception
///            escapes the unit.
///          – SINGLE READ: one placed read → depth 1 over its [p, p+L) span and 0
///            elsewhere (INV-04), Σ depth = read length, caller mean = L/G; a read
///            that fails to place (below minOverlap, or longer than the reference)
///            contributes 0 everywhere (INV-05), with NO IndexOutOfRange.
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Coverage_Calculation.md §2.2, §2.4, §3)
/// ───────────────────────────────────────────────────────────────────────────
/// For a reference of length G and placed reads, depth[i] = #{ reads r : pᵣ ≤ i &lt;
/// pᵣ + Lᵣ } (§2.2). Invariants: INV-01 output length == G; INV-02 depth[i] ≥ 0
/// (a count of reads); INV-03 Σ depth[i] == total mapped bases; INV-04 a read at p
/// of length L increments exactly [p, min(p+L, G)); INV-05 a read that fails to
/// place (best match &lt; minOverlap, or longer than G) adds 0 everywhere (§2.4).
/// Summary statistics are CALLER-derived from the array: average depth C = Σ depth
/// / G, breadth = #{ i : depth[i] ≥ 1 } / G (§2.2, §5.2). Edge cases (§3.3, §6.1):
/// empty reference → empty array; empty reads → all-zero array of length G; null
/// reference or reads → ArgumentNullException. Placement is case-insensitive,
/// ungapped best-match with a minOverlap floor (default 20); ties keep the leftmost
/// (§4.2). Worked example (§7.1): reference "ACGTTGCAAT" (G=10) with reads "ACGTT",
/// "TTGCA", "GCAAT" placed at 0/3/5 → depth [1,1,1,2,2,2,2,2,1,1], Σ=15, mean 1.5,
/// breadth 1.0.
///   SequenceAssembler.CalculateCoverage(
///       string reference, IReadOnlyList&lt;string&gt; reads, int minOverlap = 20)
///   → int[]
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyCoverageFuzzTests
{
    // Documented default placement floor (Coverage_Calculation.md §3.1).
    private const int DefaultMinOverlap = 20;

    #region Helpers

    private static readonly char[] DnaAlphabet = { 'A', 'C', 'G', 'T' };

    /// <summary>A random DNA reference of the given length over {A,C,G,T}.</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(DnaAlphabet[rng.Next(DnaAlphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Asserts a depth array is WELL-FORMED per the documented contract regardless
    /// of the (possibly degenerate) input:
    ///   INV-01 length == reference length;
    ///   INV-02 every element ≥ 0 (a count of reads, never negative);
    ///   the caller-derived mean (Σ depth / G) is finite for G &gt; 0 (no overflow /
    ///   NaN from valid counts), and is left to the caller to GUARD when G == 0.
    /// </summary>
    private static void AssertWellFormed(int[] depth, string reference)
    {
        depth.Should().NotBeNull("the unit always returns an array (never null)");
        depth.Should().HaveCount(reference.Length, "INV-01: depth length == reference length");
        depth.Should().OnlyContain(d => d >= 0, "INV-02: depth is a count of reads, never negative");

        if (reference.Length > 0)
        {
            double mean = depth.Sum() / (double)reference.Length;
            double.IsFinite(mean).Should().BeTrue("the caller-derived mean Σ depth / G is finite for G > 0");
            mean.Should().BeGreaterThanOrEqualTo(0.0, "mean coverage is non-negative");
        }
    }

    /// <summary>
    /// Independent oracle for the documented depth rule (Coverage_Calculation.md
    /// §2.2 / §4.1) restricted to EXACT placements: each read is located at its
    /// single leftmost exact occurrence (a special case of "best ungapped match" in
    /// which the read matches perfectly, isolating the source-defined counting rule
    /// per §5.4 deviation #1), then increments [p, min(p+L, G)). Reads with no exact
    /// occurrence are treated as unplaced (contribute 0, INV-05). Used to cross-check
    /// the unit on inputs built from exact substrings of the reference.
    /// </summary>
    private static int[] OracleExact(string reference, IReadOnlyList<string> reads)
    {
        var depth = new int[reference.Length];
        foreach (string read in reads)
        {
            if (read.Length == 0 || read.Length > reference.Length) continue;
            int p = reference.IndexOf(read, StringComparison.OrdinalIgnoreCase);
            if (p < 0) continue;
            for (int i = p; i < p + read.Length && i < reference.Length; i++)
                depth[i]++;
        }
        return depth;
    }

    #endregion

    #region ASSEMBLY-COVER-001 — Coverage (Depth) Calculation (BE: no reads, zero-length ref, single read)

    #region Positive sanity — hand-computed documented depth, mean and breadth

    // Documented worked example (§7.1): reference "ACGTTGCAAT" (G=10) with three unique
    // 5-mers placed at 0/3/5 → depth [1,1,1,2,2,2,2,2,1,1], Σ=15, mean 1.5, breadth 1.0.
    [Test]
    public void CalculateCoverage_WorkedExample_DocumentedPerBaseDepthMeanBreadth()
    {
        const string reference = "ACGTTGCAAT"; // length 10
        var reads = new[] { "ACGTT", "TTGCA", "GCAAT" }; // placed at 0, 3, 5

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        depth.Should().Equal(new[] { 1, 1, 1, 2, 2, 2, 2, 2, 1, 1 }, "documented per-base depth (§7.1)");
        AssertWellFormed(depth, reference);

        depth.Sum().Should().Be(15, "INV-03: Σ depth == total mapped bases (5+5+5)");
        (depth.Sum() / (double)reference.Length).Should().Be(1.5, "mean = Σ depth / G = 15/10 (§7.1)");
        (depth.Count(d => d >= 1) / (double)reference.Length).Should().Be(1.0, "breadth = covered/G = 10/10 (§7.1)");
    }

    // A read tiling a reference exactly once → depth 1 everywhere, mean 1.0, breadth 1.0.
    [Test]
    public void CalculateCoverage_SingleFullLengthRead_DepthOneEverywhere()
    {
        const string reference = "ACGTACGTAC"; // G = 10
        var reads = new[] { reference };        // full-length read at position 0

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 10);

        depth.Should().OnlyContain(d => d == 1, "a full-length read covers every position exactly once (INV-04)");
        AssertWellFormed(depth, reference);
        depth.Sum().Should().Be(10, "INV-03: Σ depth == read length == G");
    }

    // Non-overlapping tiling reads → depth 1 in each tile, exact partition of the reference.
    [Test]
    public void CalculateCoverage_NonOverlappingTiles_DepthOneAcrossPartition()
    {
        const string reference = "ACGTTGCAAT"; // G = 10
        var reads = new[] { "ACGTT", "GCAAT" }; // [0,5) and [5,10)

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        depth.Should().Equal(new[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, "two disjoint tiles cover every base once");
        depth.Sum().Should().Be(10, "INV-03: Σ depth == 5 + 5");
        AssertWellFormed(depth, reference);
    }

    #endregion

    #region BE — Boundary: no reads (empty reads list ⇒ all-zero array, caller mean 0/G, no DivideByZero)

    // Empty reads list → all-zero depth array of length G (§3.3, §6.1, INV-05).
    [Test]
    public void CalculateCoverage_NoReads_AllZeroArrayOfReferenceLength()
    {
        const string reference = "ACGTACGTAC"; // G = 10

        int[] depth = SequenceAssembler.CalculateCoverage(reference, Array.Empty<string>());

        depth.Should().HaveCount(10, "INV-01: length == reference length even with no reads");
        depth.Should().OnlyContain(d => d == 0, "no mapped read contributes any depth (INV-05)");
        AssertWellFormed(depth, reference);
    }

    // No mapped reads → caller deriving the mean Σ depth / G gets 0/G = 0, NOT a DivideByZero.
    [Test]
    public void CalculateCoverage_NoReads_CallerMeanIsZero_NoDivideByZero()
    {
        const string reference = "ACGTACGTAC"; // G = 10 > 0

        int[] depth = SequenceAssembler.CalculateCoverage(reference, Array.Empty<string>());

        Action deriveMean = () =>
        {
            double mean = depth.Sum() / (double)reference.Length;
            mean.Should().Be(0.0, "0 mapped bases / G = 0");
        };
        deriveMean.Should().NotThrow("a populated reference with no reads divides by G, not zero");
    }

    // Reads present but ALL fail to place (below minOverlap) → all-zero array (INV-05).
    [Test]
    public void CalculateCoverage_ReadsAllBelowMinOverlap_AllZero()
    {
        const string reference = "ACGTACGTACGTACGTACGTACGT"; // G = 24
        // Reads share fewer than minOverlap (20) matching chars anywhere → unplaced.
        var reads = new[] { "TTTTT", "GGGGG", "CCCCC" };

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 20);

        depth.Should().OnlyContain(d => d == 0, "reads below minOverlap add 0 everywhere (INV-05)");
        AssertWellFormed(depth, reference);
    }

    #endregion

    #region BE — Boundary: zero-length reference (empty ref ⇒ empty array, no crash, no DivideByZero in unit)

    // Empty reference → empty (length-0) depth array (§3.3, §6.1). The unit itself does
    // not divide; a caller deriving the mean must guard G == 0 (0/0). No exception escapes.
    [Test]
    public void CalculateCoverage_EmptyReference_ReturnsEmptyArray()
    {
        int[] depth = SequenceAssembler.CalculateCoverage(string.Empty, new[] { "ACGT", "TTGCA" });

        depth.Should().BeEmpty("an empty reference yields an empty depth array (§3.3, §6.1)");
    }

    // Empty reference + empty reads → still an empty array; no DivideByZero on the
    // length-0 denominator inside the unit (the unit returns before any division).
    [Test]
    public void CalculateCoverage_EmptyReferenceNoReads_EmptyArray_NoCrash()
    {
        Action act = () =>
        {
            int[] depth = SequenceAssembler.CalculateCoverage(string.Empty, Array.Empty<string>());
            depth.Should().BeEmpty("zero-length reference → empty array regardless of reads");
        };

        act.Should().NotThrow("the unit must not divide by the length-0 reference denominator");
    }

    // A caller deriving the mean from a zero-length reference must GUARD G == 0; the unit
    // hands back an empty array (Σ = 0), and the documented summary is undefined/guarded.
    [Test]
    public void CalculateCoverage_EmptyReference_CallerMustGuardZeroLengthDenominator()
    {
        int[] depth = SequenceAssembler.CalculateCoverage(string.Empty, new[] { "ACGT" });

        depth.Sum().Should().Be(0, "no positions exist to accumulate depth");
        // Documented summary mean = Σ depth / G with G == 0 is 0/0; the CALLER guards it.
        double meanGuarded = depth.Length == 0 ? 0.0 : depth.Sum() / (double)depth.Length;
        meanGuarded.Should().Be(0.0, "a G==0 guard yields a defined mean of 0 (no DivideByZero)");
    }

    #endregion

    #region BE — Boundary: single read (depth 1 over its span, 0 elsewhere, caller mean L/G; unplaced ⇒ 0)

    // One read placed in the interior → depth 1 over exactly its span, 0 elsewhere (INV-04).
    [Test]
    public void CalculateCoverage_SingleInteriorRead_DepthOneOverSpanZeroElsewhere()
    {
        const string reference = "AAACGTGTTT"; // G = 10
        var reads = new[] { "CGTGT" };          // exact occurrence at position 3 → [3,8)

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        depth.Should().Equal(new[] { 0, 0, 0, 1, 1, 1, 1, 1, 0, 0 }, "single read covers exactly [3,8) (INV-04)");
        depth.Sum().Should().Be(5, "INV-03: Σ depth == read length 5");
        (depth.Sum() / (double)reference.Length).Should().Be(0.5, "mean = L/G = 5/10");
        AssertWellFormed(depth, reference);
    }

    // One read at the reference END → depth 1 over its trailing span, no IndexOutOfRange.
    [Test]
    public void CalculateCoverage_SingleReadAtReferenceEnd_NoIndexOutOfRange()
    {
        const string reference = "ACGTTGCAAT"; // G = 10
        var reads = new[] { "GCAAT" };          // occurrence at position 5 → [5,10)

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        depth.Should().Equal(new[] { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 }, "trailing read covers [5,10), no overrun (INV-04)");
        AssertWellFormed(depth, reference);
    }

    // A SINGLE read LONGER than the reference cannot be placed → all-zero, NO IndexOutOfRange (§6.1, INV-05).
    [Test]
    public void CalculateCoverage_SingleReadLongerThanReference_AllZeroNoOverrun()
    {
        const string reference = "ACGTT"; // G = 5
        var reads = new[] { "ACGTTGCAAT" };  // length 10 > 5 → cannot be placed

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        depth.Should().HaveCount(5, "INV-01: array length == reference length");
        depth.Should().OnlyContain(d => d == 0, "an over-long read cannot place → 0 everywhere (INV-05)");
        AssertWellFormed(depth, reference);
    }

    // A single read whose best match is below minOverlap → unplaced, all-zero (INV-05).
    [Test]
    public void CalculateCoverage_SingleReadBelowMinOverlap_AllZero()
    {
        const string reference = "ACGTACGTACGTACGTACGTACGT"; // G = 24
        var reads = new[] { "TTTTTTTTTTTTTTTTTTTT" };          // 20 T's: best match < 20 → unplaced

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 20);

        depth.Should().OnlyContain(d => d == 0, "a read below minOverlap adds 0 everywhere (INV-05)");
        AssertWellFormed(depth, reference);
    }

    // A single read matching case-insensitively places and counts (lowercase read, §3.3).
    [Test]
    public void CalculateCoverage_SingleLowercaseRead_PlacesCaseInsensitively()
    {
        const string reference = "ACGTTGCAAT"; // G = 10
        var reads = new[] { "acgtt" };           // lowercase exact match at 0 → [0,5)

        int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: 5);

        depth.Should().Equal(new[] { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 }, "case-insensitive placement (§3.3)");
        AssertWellFormed(depth, reference);
    }

    #endregion

    #region BE — Validation: null reference / null reads ⇒ documented ArgumentNullException (§3.3, §6.1)

    [Test]
    public void CalculateCoverage_NullReference_Throws()
    {
        Action act = () => SequenceAssembler.CalculateCoverage(null!, new[] { "ACGT" });

        act.Should().Throw<ArgumentNullException>("null reference is the documented validation contract (§3.3)");
    }

    [Test]
    public void CalculateCoverage_NullReads_Throws()
    {
        Action act = () => SequenceAssembler.CalculateCoverage("ACGT", null!);

        act.Should().Throw<ArgumentNullException>("null reads list is the documented validation contract (§3.3)");
    }

    #endregion

    #region BE — Broad fuzz: degenerate refs/reads never crash; exact-substring reads match the oracle

    // Fuzz the boundary axes (G ∈ {0,1,…}, #reads ∈ {0,…}) directly: the unit must NEVER
    // throw on a non-null reference + non-null reads, always returns a well-formed array,
    // and a caller deriving the mean with a G==0 guard never divides by zero.
    [Test]
    [CancelAfter(30_000)]
    public void CalculateCoverage_FuzzedBoundaryShapes_NeverThrows_WellFormed()
    {
        var rng = new Random(142_001);
        for (int trial = 0; trial < 1500; trial++)
        {
            int g = rng.Next(0, 30);              // includes the zero-length reference boundary
            string reference = RandomDna(rng, g);

            int nReads = rng.Next(0, 6);          // includes the no-reads boundary
            var reads = new List<string>(nReads);
            for (int r = 0; r < nReads; r++)
                reads.Add(RandomDna(rng, rng.Next(0, 12))); // includes empty / over-long reads

            int minOverlap = rng.Next(0, 8);

            int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap);

            AssertWellFormed(depth, reference);

            // Caller-side summary with the documented G==0 guard never divides by zero.
            Action deriveSummaries = () =>
            {
                double mean = depth.Length == 0 ? 0.0 : depth.Sum() / (double)depth.Length;
                double breadth = depth.Length == 0 ? 0.0 : depth.Count(d => d >= 1) / (double)depth.Length;
                double.IsFinite(mean).Should().BeTrue();
                breadth.Should().BeInRange(0.0, 1.0);
            };
            deriveSummaries.Should().NotThrow("guarded summaries never DivideByZero on G==0 / no reads");
        }
    }

    // Build reads as EXACT substrings of a random reference so placement is deterministic
    // (each read at its leftmost exact occurrence): the unit's depth array must equal the
    // independent exact-placement oracle, and Σ depth == total placed bases (INV-03).
    [Test]
    [CancelAfter(30_000)]
    public void CalculateCoverage_ExactSubstringReads_MatchesExactOracle()
    {
        var rng = new Random(142_002);
        for (int trial = 0; trial < 800; trial++)
        {
            int g = rng.Next(1, 30);
            string reference = RandomDna(rng, g);

            int nReads = rng.Next(0, 6);
            var reads = new List<string>(nReads);
            int minLenReadSeen = int.MaxValue;
            for (int r = 0; r < nReads; r++)
            {
                int len = rng.Next(1, g + 1);          // 1..G, fits in the reference
                int start = rng.Next(0, g - len + 1);
                string read = reference.Substring(start, len);
                reads.Add(read);
                minLenReadSeen = Math.Min(minLenReadSeen, len);
            }

            // minOverlap ≤ every read length so each exact read can place.
            int minOverlap = nReads == 0 ? 1 : Math.Min(minLenReadSeen, 1);

            int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap);

            AssertWellFormed(depth, reference);
            int[] expected = OracleExact(reference, reads);
            depth.Should().Equal(expected, "exact-substring placement matches the documented counting rule (§2.2)");

            // INV-03 sanity on the oracle's deterministic placements.
            depth.Sum().Should().Be(expected.Sum(), "INV-03: Σ depth == total mapped bases");
        }
    }

    #endregion

    #endregion
}
