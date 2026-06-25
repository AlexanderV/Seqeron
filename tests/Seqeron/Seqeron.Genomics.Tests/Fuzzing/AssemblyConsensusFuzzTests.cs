using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Assembly area — Consensus Computation (ASSEMBLY-CONSENSUS-001),
/// the column-wise majority/threshold consensus caller
/// <see cref="SequenceAssembler.ComputeConsensus(IReadOnlyList{string}, double, char)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop
/// (the O(L×R) two-pass-per-column reduction must always terminate), no state
/// corruption, no nonsense output (a consensus whose length differs from the
/// longest read, an emitted gap '-'/'.', a committed residue that is not the unique
/// strict-max above threshold), and no *unhandled* runtime exception — in
/// particular NO DivideByZero on an empty or all-gap (zero non-gap-depth) column,
/// and NO nondeterministic tie-break. Every input must resolve to EITHER a
/// well-defined, theory-correct result OR a *documented, intentional* validation
/// exception (ArgumentNullException for a null read list — contract §3.3, §6.1).
/// A raw runtime exception, a hang, a wrong-length consensus, an emitted gap, or an
/// order-dependent tie winner is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-CONSENSUS-001 — Consensus Computation
/// Checklist: docs/checklists/03_FUZZING.md, row 140.
/// Algorithm doc: docs/algorithms/Assembly/Consensus_Computation.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row:
///          – SINGLE READ: one read (and a 1-deep alignment) → the consensus equals
///            that read (uppercased, gaps→ambiguous), with NO DivideByZero / variance
///            crash on a depth-1 column (§6.1, INV-02: 1/1 = 1.0 ≥ threshold).
///          – CONFLICTING READS: reads disagreeing at a column → the DOCUMENTED,
///            DETERMINISTIC outcome: the unique strict-max residue iff its frequency
///            ≥ threshold, else the ambiguous symbol; a TIE for the max emits the
///            ambiguous symbol regardless of read order (INV-02/INV-03 — not
///            order-dependent garbage, never an arbitrary MaxBy winner).
///          – EMPTY: no reads → the documented empty string "" (no columns), with NO
///            DivideByZero on a zero-depth column (§6.1, §3.3).
///          – ALL-N: every read all-'N' → every column commits 'N' (N is a plain
///            residue here, not a gap), no crash; the all-GAP analogue ('-'/'.') →
///            every column emits the ambiguous symbol with NO DivideByZero (INV-04,
///            INV-05).
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Consensus_Computation.md §2.4, §3)
/// ───────────────────────────────────────────────────────────────────────────
/// For each alignment column n (0 ≤ n &lt; L, L = longest read), tally only NON-GAP
/// residues ('-' and '.' skipped, INV-04); num_atoms = number of such residues.
/// Let max_size be the maximum count and max_atoms the residues achieving it. A
/// residue is committed iff it is the UNIQUE max (len(max_atoms)==1, INV-03) AND
/// max_size / num_atoms ≥ threshold (INV-02); otherwise the ambiguous symbol is
/// emitted. Residues are uppercased. The consensus length equals L (INV-01). An
/// all-gap / empty column emits ambiguous with NO division by zero (INV-05). Empty
/// read list → "" (§6.1); null list → ArgumentNullException (§3.3, §6.1).
/// Defaults: threshold = 0.5 (simple-majority plurality), ambiguous = 'N'.
///   SequenceAssembler.ComputeConsensus(
///       IReadOnlyList&lt;string&gt; alignedReads, double threshold = 0.5, char ambiguous = 'N')
///   → string
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyConsensusFuzzTests
{
    // Documented defaults (Consensus_Computation.md §3.1, §4.2).
    private const double DefaultThreshold = 0.5;
    private const char DefaultAmbiguous = 'N';
    private const char GapDash = '-';
    private const char GapDot = '.';

    #region Helpers

    /// <summary>
    /// Independent oracle implementing the documented decision rule verbatim
    /// (Consensus_Computation.md §2.4 / §4.1): per column, tally non-gap residues,
    /// commit the UNIQUE strict-max residue iff its frequency ≥ threshold, else
    /// emit <paramref name="ambiguous"/>. Used to cross-check the unit under test on
    /// fuzzed inputs without re-using its implementation.
    /// </summary>
    private static string Oracle(IReadOnlyList<string> reads, double threshold, char ambiguous)
    {
        int length = reads.Count == 0 ? 0 : reads.Max(r => r.Length);
        var sb = new StringBuilder(length);

        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int>();
            int numAtoms = 0;
            foreach (string r in reads)
            {
                if (pos >= r.Length) continue;
                char c = char.ToUpperInvariant(r[pos]);
                if (c == GapDash || c == GapDot) continue;
                counts[c] = counts.GetValueOrDefault(c, 0) + 1;
                numAtoms++;
            }

            int maxSize = counts.Count == 0 ? 0 : counts.Values.Max();
            int tiedAtMax = counts.Values.Count(v => v == maxSize);
            char maxResidue = counts.Where(kvp => kvp.Value == maxSize)
                                    .Select(kvp => kvp.Key)
                                    .FirstOrDefault(ambiguous);

            bool committed = maxSize > 0
                             && tiedAtMax == 1
                             && (double)maxSize / numAtoms >= threshold;

            sb.Append(committed ? maxResidue : ambiguous);
        }

        return sb.ToString();
    }

    /// <summary>
    /// Asserts a consensus is WELL-FORMED per the documented contract regardless of
    /// the (possibly degenerate) input:
    ///   INV-01 length == longest read (== L; 0 for an empty list);
    ///   INV-04 no gap character '-'/'.' ever appears in the output;
    ///   every emitted character is either the ambiguous symbol or an UPPERCASE
    ///   residue that actually occurs (uppercased) somewhere in that column
    ///   (committed residues are drawn from the column's non-gap residues, §3.2).
    /// </summary>
    private static void AssertWellFormed(
        string consensus, IReadOnlyList<string> reads, char ambiguous)
    {
        int length = reads.Count == 0 ? 0 : reads.Max(r => r.Length);
        consensus.Should().HaveLength(length, "INV-01: consensus length == longest read (alignment length)");

        for (int pos = 0; pos < consensus.Length; pos++)
        {
            char c = consensus[pos];
            c.Should().NotBe(GapDash, "INV-04: gap '-' never appears in the consensus");
            c.Should().NotBe(GapDot, "INV-04: gap '.' never appears in the consensus");

            if (c == ambiguous) continue;

            // A committed (non-ambiguous) residue must be an uppercased residue
            // present in this column among the supplied reads.
            var columnResidues = reads
                .Where(r => pos < r.Length)
                .Select(r => char.ToUpperInvariant(r[pos]))
                .Where(ch => ch != GapDash && ch != GapDot)
                .ToHashSet();
            columnResidues.Should().Contain(c,
                "a committed residue must occur (uppercased) in its column (§3.2)");
        }
    }

    private static readonly char[] DnaAlphabet = { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', '-', '.', 'N' };

    /// <summary>A fuzzed pre-aligned read of the given length over the DNA + gap + N alphabet.</summary>
    private static string RandomRead(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(DnaAlphabet[rng.Next(DnaAlphabet.Length)]);
        return sb.ToString();
    }

    #endregion

    #region ASSEMBLY-CONSENSUS-001 — Consensus Computation (BE: single read, conflicting reads, empty, all-N)

    #region Positive sanity — hand-computed documented consensus

    // Mostly-agreeing reads: every column has a clear ≥ 0.5 majority (hand-computed).
    //   col0 {A,A,A,T} A=3/4=0.75 ≥ 0.5 → A
    //   col1 {C,C,C,C} C=4/4=1.0      → C
    //   col2 {G,G,A,G} G=3/4=0.75     → G
    //   col3 {T,T,T,T} T=4/4=1.0      → T
    [Test]
    public void ComputeConsensus_MostlyAgreeingReads_DocumentedMajority()
    {
        var reads = new[] { "ACGT", "ACGT", "ACAT", "TCGT" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        consensus.Should().Be("ACGT", "every column has a unique strict-max residue with frequency ≥ 0.5");
        AssertWellFormed(consensus, reads, DefaultAmbiguous);
    }

    // Documented worked example (§7.1): threshold flips a 2/3 column from committed to ambiguous.
    [Test]
    public void ComputeConsensus_WorkedExample_ThresholdControlsColumn0()
    {
        var reads = new[] { "ACGT", "ACGT", "TCGT" };

        SequenceAssembler.ComputeConsensus(reads, threshold: 0.5)
            .Should().Be("ACGT", "col0 A=2/3≈0.667 ≥ 0.5 → A (§7.1)");
        SequenceAssembler.ComputeConsensus(reads, threshold: 0.7)
            .Should().Be("NCGT", "col0 A=2/3≈0.667 < 0.7 → N (§7.1)");
    }

    // A clear conflict resolves to the documented strict-max winner (deterministic, not MaxBy garbage).
    [Test]
    public void ComputeConsensus_ClearConflict_ResolvesToStrictMaxWinner()
    {
        // col0 {G,G,G,A,C} G=3/5=0.6 ≥ 0.5, unique max → G.
        var reads = new[] { "G", "G", "G", "A", "C" };

        SequenceAssembler.ComputeConsensus(reads)
            .Should().Be("G", "unique strict-max residue above threshold is committed (INV-02)");
    }

    #endregion

    #region BE — Boundary: single read (depth-1 column ⇒ that read, no DivideByZero)

    // One read → consensus equals that read (1/1 = 1.0 ≥ threshold for every non-gap column);
    // no DivideByZero / variance crash on a depth-1 column.
    [Test]
    public void ComputeConsensus_SingleRead_EqualsThatRead()
    {
        var reads = new[] { "ACGTACGTAC" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        consensus.Should().Be("ACGTACGTAC", "a single read commits every non-gap column (1/1 ≥ threshold)");
        AssertWellFormed(consensus, reads, DefaultAmbiguous);
    }

    // Single read with lowercase + gaps: residues are uppercased; gap columns become ambiguous;
    // length still equals the read length (INV-01, INV-04).
    [Test]
    public void ComputeConsensus_SingleReadWithGapsAndCase_UppercasedGapsAmbiguous()
    {
        var reads = new[] { "ac-gt.A" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        consensus.Should().Be("ACNGTNA", "uppercase residues, gap columns ('-','.') → ambiguous 'N'");
        AssertWellFormed(consensus, reads, DefaultAmbiguous);
    }

    // Fuzz: a single random read always yields a same-length consensus, never throws.
    [Test]
    public void ComputeConsensus_SingleRandomRead_NeverThrows_WellFormed()
    {
        var rng = new Random(140_001);
        for (int trial = 0; trial < 400; trial++)
        {
            int len = rng.Next(0, 40);
            var reads = new[] { RandomRead(rng, len) };

            string consensus = SequenceAssembler.ComputeConsensus(reads);

            AssertWellFormed(consensus, reads, DefaultAmbiguous);
            consensus.Should().Be(Oracle(reads, DefaultThreshold, DefaultAmbiguous),
                "single-read consensus matches the documented decision rule");
        }
    }

    #endregion

    #region BE — Boundary: conflicting reads (deterministic majority / tie → ambiguous, order-independent)

    // A 50/50 tie for the maximum count → ambiguous symbol (INV-03), NOT an arbitrary winner.
    [Test]
    public void ComputeConsensus_TwoWayTie_EmitsAmbiguous()
    {
        // col0 {A,T} both count 1 → tie → N. col1 {C,C} → C.
        var reads = new[] { "AC", "TC" };

        SequenceAssembler.ComputeConsensus(reads)
            .Should().Be("NC", "a tie for the max count emits the ambiguous symbol (INV-03)");
    }

    // Three-way tie (A=C=G each 1 of 3) → ambiguous; deterministic.
    [Test]
    public void ComputeConsensus_ThreeWayTie_EmitsAmbiguous()
    {
        var reads = new[] { "A", "C", "G" };

        SequenceAssembler.ComputeConsensus(reads)
            .Should().Be("N", "three residues tied for the max count → ambiguous (INV-03)");
    }

    // Tie-break MUST be order-independent: every permutation of the SAME conflicting reads
    // yields the IDENTICAL consensus (not order-dependent garbage, not arbitrary MaxBy).
    [Test]
    public void ComputeConsensus_ConflictingReads_DeterministicAcrossPermutations()
    {
        var rng = new Random(140_002);
        for (int trial = 0; trial < 200; trial++)
        {
            // Build a column-conflicting read set of equal-length reads.
            int depth = rng.Next(2, 7);
            int width = rng.Next(1, 8);
            var reads = Enumerable.Range(0, depth).Select(_ => RandomRead(rng, width)).ToList();

            string baseline = SequenceAssembler.ComputeConsensus(reads);

            // Shuffle the read order several times; consensus must be invariant.
            for (int shuffle = 0; shuffle < 4; shuffle++)
            {
                var permuted = reads.OrderBy(_ => rng.Next()).ToList();
                SequenceAssembler.ComputeConsensus(permuted)
                    .Should().Be(baseline, "consensus is order-independent (deterministic tie-break, INV-03)");
            }

            AssertWellFormed(baseline, reads, DefaultAmbiguous);
        }
    }

    // Sub-threshold plurality (no residue reaches the cut-off) → ambiguous (INV-02).
    [Test]
    public void ComputeConsensus_SubThresholdPlurality_EmitsAmbiguous()
    {
        // col0 {A,A,C,G,T} A=2/5=0.4 < 0.5, unique max but below threshold → N.
        var reads = new[] { "A", "A", "C", "G", "T" };

        SequenceAssembler.ComputeConsensus(reads)
            .Should().Be("N", "a unique max below the threshold emits ambiguous (INV-02)");
    }

    // Custom ambiguous symbol (e.g. protein 'X') is honoured on no-consensus columns.
    [Test]
    public void ComputeConsensus_CustomAmbiguousSymbol_UsedOnNoConsensus()
    {
        var reads = new[] { "AC", "TC" }; // col0 tie → ambiguous, col1 C.

        SequenceAssembler.ComputeConsensus(reads, threshold: 0.5, ambiguous: 'X')
            .Should().Be("XC", "the supplied ambiguous symbol replaces 'N' (§3.1)");
    }

    #endregion

    #region BE — Boundary: empty (no reads / empty reads ⇒ no DivideByZero on zero-depth columns)

    // No reads → empty string, no columns, NO DivideByZero (§6.1).
    [Test]
    public void ComputeConsensus_EmptyReadList_ReturnsEmptyString()
    {
        SequenceAssembler.ComputeConsensus(Array.Empty<string>())
            .Should().BeEmpty("an empty read list has no columns → \"\" (§6.1)");
    }

    // Reads that are all empty strings → longest read is 0 → empty consensus, no crash.
    [Test]
    public void ComputeConsensus_AllEmptyReads_ReturnsEmptyString()
    {
        var reads = new[] { "", "", "" };

        SequenceAssembler.ComputeConsensus(reads)
            .Should().BeEmpty("longest read length 0 ⇒ empty consensus, no zero-depth division");
    }

    // Null read list → documented ArgumentNullException (§3.3, §6.1).
    [Test]
    public void ComputeConsensus_NullReadList_Throws()
    {
        Action act = () => SequenceAssembler.ComputeConsensus(null!);

        act.Should().Throw<ArgumentNullException>("null read list is the documented validation contract (§3.3)");
    }

    // Ragged reads with a zero-depth interior column: a column covered by no read (all reads
    // shorter, or all gaps) yields the ambiguous symbol with NO DivideByZero (INV-05).
    [Test]
    public void ComputeConsensus_AllGapColumn_AmbiguousNoDivideByZero()
    {
        // col0 {A,A} → A; col1 {-,.} all gaps, num_atoms=0 → ambiguous, NO /0 (INV-05).
        var reads = new[] { "A-", "A." };

        SequenceAssembler.ComputeConsensus(reads)
            .Should().Be("AN", "an all-gap (zero non-gap-depth) column emits ambiguous, no /0 (INV-05)");
    }

    // Ragged reads: the consensus spans the LONGEST read; columns past shorter reads have
    // reduced depth but never divide by zero (INV-01, INV-05).
    [Test]
    public void ComputeConsensus_RaggedReads_LengthEqualsLongest()
    {
        var reads = new[] { "ACGTACG", "AC", "ACG" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        consensus.Should().HaveLength(7, "INV-01: consensus length == longest read");
        AssertWellFormed(consensus, reads, DefaultAmbiguous);
        // col0 {A,A,A}→A col1 {C,C,C}→C col2 {G,G}→G col3..6 {single} → from the long read.
        consensus.Should().Be("ACGTACG", "ragged tail columns commit the sole covering read");
    }

    #endregion

    #region BE — Boundary: all-N reads (N is a plain residue ⇒ N consensus; all-gap ⇒ ambiguous)

    // Every read all-'N': N is a normal residue (NOT a gap), so every column commits 'N'.
    [Test]
    public void ComputeConsensus_AllNReads_CommitsN()
    {
        var reads = new[] { "NNNN", "NNNN", "NNNN" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        consensus.Should().Be("NNNN", "N is a plain residue: unique max at frequency 1.0 → committed N");
        AssertWellFormed(consensus, reads, DefaultAmbiguous);
    }

    // A single all-N read: 1/1 ≥ threshold → all N, no crash on depth-1 N columns.
    [Test]
    public void ComputeConsensus_SingleAllNRead_CommitsN()
    {
        var reads = new[] { "NNNNN" };

        SequenceAssembler.ComputeConsensus(reads)
            .Should().Be("NNNNN", "a single all-N read commits N at every column");
    }

    // Every read all-GAP ('-'/'.'): every column has zero non-gap depth → all ambiguous, NO /0.
    [Test]
    public void ComputeConsensus_AllGapReads_AllAmbiguous_NoDivideByZero()
    {
        var reads = new[] { "----", "....", "-.-." };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        consensus.Should().Be("NNNN", "all-gap columns commit nothing → ambiguous, no /0 (INV-05)");
        AssertWellFormed(consensus, reads, DefaultAmbiguous);
    }

    #endregion

    #region BE — Broad fuzz: random aligned reads never crash, always match the documented rule

    [Test]
    [CancelAfter(30_000)]
    public void ComputeConsensus_RandomAlignedReads_NeverThrows_MatchesOracle()
    {
        var rng = new Random(140_003);
        for (int trial = 0; trial < 1500; trial++)
        {
            int depth = rng.Next(0, 9);
            int width = rng.Next(0, 30);
            var reads = Enumerable.Range(0, depth).Select(_ => RandomRead(rng, width)).ToList();
            double threshold = rng.NextDouble(); // 0..1 inclusive-of-0
            char ambiguous = rng.Next(2) == 0 ? 'N' : 'X';

            string consensus = SequenceAssembler.ComputeConsensus(reads, threshold, ambiguous);

            AssertWellFormed(consensus, reads, ambiguous);
            consensus.Should().Be(Oracle(reads, threshold, ambiguous),
                "the unit matches the documented decision rule under fuzzed depth/width/threshold");
        }
    }

    // Boundary thresholds 0.0 and 1.0: at threshold 0, ANY unique-max residue commits (≥ 0 always
    // true); at threshold 1.0, only unanimous non-gap columns commit. No crash either way.
    [Test]
    public void ComputeConsensus_BoundaryThresholds_BehaveDeterministically()
    {
        var reads = new[] { "AAAC", "AAGC", "AATC" };

        // threshold 0.0: every column's unique max commits (col1 A=3/3, col2 tie? A=1,G=1,T=1 tie → N).
        // col0 A=3/3→A; col1 A=3/3→A; col2 {A,G,T} tie → N even at 0.0 (tie ⇒ ambiguous, INV-03); col3 C→C.
        SequenceAssembler.ComputeConsensus(reads, threshold: 0.0)
            .Should().Be("AANC", "threshold 0 still emits ambiguous on a tie (INV-03 dominates)");

        // threshold 1.0: only unanimous non-gap columns commit; col0/col1/col3 unanimous, col2 not.
        SequenceAssembler.ComputeConsensus(reads, threshold: 1.0)
            .Should().Be("AANC", "threshold 1.0 commits only unanimous columns");
    }

    #endregion

    #endregion
}
