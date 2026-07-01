using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Matching area — IUPAC-Degenerate Consensus Generation
/// (MOTIF-GENERATE-001), the threshold-then-encode consensus caller
/// <see cref="MotifFinder.GenerateConsensus(IEnumerable{string})"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop
/// (the O(n×m) per-column tally must always terminate), no state corruption, and
/// no nonsense output — a consensus whose length differs from the column count,
/// an emitted symbol outside the 15 IUPAC codes, or a non-deterministic tie/threshold
/// winner. The unit forms a per-column COUNT MATRIX (the A/C/G/T tally for each
/// column of the supplied equal-length rows) and emits one symbol per column, so the
/// degenerate boundaries of a counts→consensus generator are exactly: EMPTY COUNTS
/// (no columns, or a column whose A/C/G/T total is zero — no DivideByZero, no
/// NullReference), SINGLE COLUMN (a one-column matrix → a length-1 consensus), and
/// TIES (a column where two+ bases share the maximum count — the documented
/// threshold/fallback rule must resolve them DETERMINISTICALLY). Every input must
/// resolve to EITHER a well-defined, theory-correct result OR the single documented
/// validation exception (ArgumentNullException for a null collection — §3.3, §6.1).
/// A raw runtime exception (DivideByZero on a zero-total column, NullReference /
/// IndexOutOfRange on empty/ragged input), a hang, a wrong-length consensus, an
/// out-of-alphabet symbol, or an order-dependent winner is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MOTIF-GENERATE-001 — IUPAC-Degenerate Consensus Generation
/// Checklist: docs/checklists/03_FUZZING.md, row 171.
/// Algorithm doc: docs/algorithms/Pattern_Matching/IUPAC_Degenerate_Consensus.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the row
///     ("empty counts, single column, ties"):
///       – EMPTY COUNTS: an empty collection → the documented "" (§3.3, §6.1); rows
///         that are all empty → a 0-column matrix → "" with NO DivideByZero on a
///         zero-total column; a column tallied entirely from non-ACGT residues has a
///         zero A/C/G/T total and MUST NOT divide-by-zero or null-deref — it falls
///         back to the most-frequent (zero) base 'A' deterministically (§5.2).
///       – SINGLE COLUMN: a one-column matrix (width-1 rows) → a length-1 consensus
///         (INV-01); a unanimous single column → that base (INV-02), no crash.
///       – TIES: a column with two+ equal-max bases → the documented rule: each base
///         strictly over θ·n (θ = 0.25) enters the IUPAC set and the set is encoded
///         via the NC-IUB table (INV-04); a column where NO base passes (e.g. four
///         equal bases each at 25 %) falls back to the most-frequent base with an
///         ALPHABETICAL tie-break — deterministic, never N (§5.2, §6.1).
/// — docs/checklists/03_FUZZING.md §Description (BE = граничні значення 0/-1/MaxInt/empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (IUPAC_Degenerate_Consensus.md §2.2, §3, §5.2, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// Column count = the FIRST row's length (INV-01). For each column, tally A,C,G,T
/// over the uppercased rows (non-ACGT residues are IGNORED, not rejected — §3.3);
/// shorter rows simply contribute fewer counts at trailing columns. Retain the set
/// B = { b : count(b) > θ·n }, θ = 0.25, n = number of rows, STRICT '>' (INV-05).
/// If B is non-empty, emit IUPAC(B) via the NC-IUB 1984 set→symbol table (§2.2,
/// INV-04); otherwise emit the single most-frequent base with an alphabetical
/// (A&lt;C&lt;G&lt;T) tie-break (§5.2). Output ⊆ the 15 IUPAC symbols (INV-03). Null
/// collection → ArgumentNullException; empty collection → "" (§3.3).
///   MotifFinder.GenerateConsensus(IEnumerable&lt;string&gt;) → string
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MotifGenerateFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    /// <summary>The 15 IUPAC symbols that are the only legal output characters (INV-03).</summary>
    private const string IupacSymbols = "ACGTRYSWKMBDHVN";

    private const double Theta = 0.25;

    #region Helpers

    /// <summary>
    /// Independent oracle implementing the documented decision rule verbatim
    /// (IUPAC_Degenerate_Consensus.md §2.2, §4.1, §5.2) without re-using the unit:
    /// per column, tally A,C,G,T over the uppercased rows (ignoring non-ACGT and
    /// rows too short to reach the column); keep bases whose count is STRICTLY greater
    /// than θ·n; if any pass, map the set to its NC-IUB symbol; otherwise emit the
    /// most-frequent base with an alphabetical tie-break.
    /// </summary>
    private static string Oracle(IReadOnlyList<string> rows)
    {
        int n = rows.Count;
        if (n == 0) return "";
        int width = rows[0].Length;
        var sb = new StringBuilder(width);
        double threshold = n * Theta;

        for (int col = 0; col < width; col++)
        {
            var counts = new int[Alphabet.Length]; // A,C,G,T
            foreach (string r in rows)
            {
                if (col >= r.Length) continue;            // shorter row contributes nothing here
                int idx = Array.IndexOf(Alphabet, char.ToUpperInvariant(r[col]));
                if (idx >= 0) counts[idx]++;              // non-ACGT ignored
            }

            string passing = string.Concat(
                Enumerable.Range(0, Alphabet.Length)
                          .Where(b => counts[b] > threshold) // strict '>'
                          .Select(b => Alphabet[b]));        // already in A,C,G,T order

            if (passing.Length > 0)
            {
                sb.Append(NcIub(passing));
            }
            else
            {
                // No base passes → most-frequent base, alphabetical tie-break (first max in A→C→G→T).
                int best = 0;
                for (int b = 1; b < counts.Length; b++)
                    if (counts[b] > counts[best]) best = b; // strict '>' keeps earliest on tie
                sb.Append(Alphabet[best]);
            }
        }

        return sb.ToString();
    }

    /// <summary>NC-IUB 1984 set→symbol map for a base set written in A&lt;C&lt;G&lt;T order (§2.2).</summary>
    private static char NcIub(string bases) => bases switch
    {
        "A" => 'A', "C" => 'C', "G" => 'G', "T" => 'T',
        "AG" => 'R', "CT" => 'Y', "CG" => 'S', "AT" => 'W', "GT" => 'K', "AC" => 'M',
        "CGT" => 'B', "AGT" => 'D', "ACT" => 'H', "ACG" => 'V',
        _ => 'N',
    };

    /// <summary>
    /// Asserts a consensus is WELL-FORMED per the documented contract regardless of the
    /// (possibly degenerate) input: INV-01 length == column count (== first row's length,
    /// or 0 for an empty matrix); INV-03 every output symbol is one of the 15 IUPAC codes.
    /// </summary>
    private static void AssertWellFormed(string consensus, IReadOnlyList<string> rows)
    {
        int width = rows.Count == 0 ? 0 : rows[0].Length;
        consensus.Should().HaveLength(width, "INV-01: consensus length == column count (first row's length)");
        foreach (char c in consensus)
            IupacSymbols.Should().Contain(c.ToString(),
                "INV-03: every output symbol is one of the 15 IUPAC codes");
    }

    private static string RandomRow(Random rng, int width, char[] alphabet)
    {
        var sb = new StringBuilder(width);
        for (int i = 0; i < width; i++)
            sb.Append(alphabet[rng.Next(alphabet.Length)]);
        return sb.ToString();
    }

    private static readonly char[] AlphabetMixedCase = { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' };

    #endregion

    #region MOTIF-GENERATE-001 — IUPAC-Degenerate Consensus Generation (BE: empty counts, single column, ties)

    #region Positive sanity — hand-computed documented consensus

    // Documented worked example (§7.1): col0 {A,G} (each 1 > 2×0.25=0.5) → R; cols 1-3 unanimous → TGC.
    [Test]
    public void Generate_DocumentedSample_ProducesRTGC()
    {
        var rows = new[] { "ATGC", "GTGC" };

        string consensus = MotifFinder.GenerateConsensus(rows);

        consensus.Should().Be("RTGC", "documented degenerate-consensus worked example (§7.1)");
        AssertWellFormed(consensus, rows);
        consensus.Should().Be(Oracle(rows), "matches the documented decision rule");
    }

    // Numerical walk-through (§7.1): ["C","G","T"], n=3, θ·n=0.75; each count 1 > 0.75 → {C,G,T} → B.
    [Test]
    public void Generate_DocumentedThreeBaseColumn_ProducesB()
    {
        var rows = new[] { "C", "G", "T" };

        MotifFinder.GenerateConsensus(rows)
            .Should().Be("B", "{C,G,T} each strictly over θ·n=0.75 → NC-IUB symbol B (§7.1)");
    }

    // A clear unanimous matrix reproduces the input base at every column (INV-02 singleton sets).
    [Test]
    public void Generate_UnanimousColumns_ReproduceTheBase()
    {
        var rows = new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT" };

        MotifFinder.GenerateConsensus(rows)
            .Should().Be("ACGTACGT", "every column is a unanimous singleton set → that base (INV-02)");
    }

    #endregion

    #region BE — empty counts (empty collection, all-empty rows, zero-total column)

    // Empty collection → documented "" with no columns, NO DivideByZero / NullReference (§3.3, §6.1).
    [Test]
    public void Generate_EmptyCollection_ReturnsEmptyString()
    {
        MotifFinder.GenerateConsensus(Array.Empty<string>())
            .Should().BeEmpty("an empty collection has no columns → \"\" (§3.3, §6.1)");
    }

    // Rows that are all empty strings → first-row width 0 → empty consensus, no zero-row division.
    [Test]
    public void Generate_AllEmptyRows_ReturnsEmptyString()
    {
        MotifFinder.GenerateConsensus(new[] { "", "", "" })
            .Should().BeEmpty("first row width 0 ⇒ zero columns ⇒ empty consensus, no DivideByZero");
    }

    // Null collection → documented ArgumentNullException (§3.3, §6.1).
    [Test]
    public void Generate_NullCollection_Throws()
    {
        Action act = () => MotifFinder.GenerateConsensus(null!);

        act.Should().Throw<ArgumentNullException>("null collection is the documented validation contract (§3.3)");
    }

    // A column whose A/C/G/T total is ZERO (all residues non-ACGT, which are ignored) must NOT
    // DivideByZero or null-deref: with no passing base it falls back to most-frequent (all-zero) → 'A'.
    [Test]
    public void Generate_ZeroTotalColumn_NoDivideByZero_FallsBackToA()
    {
        var rows = new[] { "N", "-", "X" }; // none counted → zero A/C/G/T total in the single column

        Action act = () => MotifFinder.GenerateConsensus(rows);

        act.Should().NotThrow("a zero-total column must not DivideByZero / null-deref (§5.2, §6.1)");
        string consensus = MotifFinder.GenerateConsensus(rows);
        consensus.Should().Be("A", "no base passes a positive threshold ⇒ most-frequent (zero) base, alphabetical → A");
        AssertWellFormed(consensus, rows);
    }

    // Fuzz: a column built ENTIRELY from non-ACGT residues never crashes and always yields a
    // legal IUPAC symbol (the deterministic all-zero fallback).
    [Test]
    [CancelAfter(30_000)]
    public void Generate_RandomNonAcgtOnlyColumns_NeverThrows_LegalSymbol()
    {
        const string nonAcgt = "NRYSWKMBDHV-.* 0123456789xz";
        var rng = new Random(171_001);
        for (int trial = 0; trial < 500; trial++)
        {
            int depth = rng.Next(1, 7);
            int width = rng.Next(1, 12);
            var rows = Enumerable.Range(0, depth)
                .Select(_ => RandomRow(rng, width, nonAcgt.ToCharArray()))
                .ToList();

            string consensus = MotifFinder.GenerateConsensus(rows);

            AssertWellFormed(consensus, rows);
            consensus.Should().Be(Oracle(rows), "matches the documented zero-total fallback rule");
        }
    }

    #endregion

    #region BE — single column (one-column matrix → length-1 consensus)

    // A one-column matrix → a length-1 consensus equal to that column's argmax (no crash).
    [Test]
    public void Generate_SingleColumn_LengthOneArgmax()
    {
        // {A,A,G}: n=3, θ·n=0.75; A=2>0.75 passes, G=1>0.75 passes → {A,G} → R.
        MotifFinder.GenerateConsensus(new[] { "A", "A", "G" })
            .Should().Be("R", "single column {A,G} both over θ·n → R (length 1)");

        // {C}: singleton → C (length 1).
        MotifFinder.GenerateConsensus(new[] { "C" })
            .Should().Be("C", "a single-row single-column matrix → that base");
    }

    // Fuzz: any single-column matrix yields a length-1 consensus, never crashes, matches the oracle.
    [Test]
    [CancelAfter(30_000)]
    public void Generate_RandomSingleColumn_AlwaysLengthOne_MatchesOracle()
    {
        var rng = new Random(171_002);
        for (int trial = 0; trial < 800; trial++)
        {
            int depth = rng.Next(1, 12);
            var rows = Enumerable.Range(0, depth)
                .Select(_ => RandomRow(rng, 1, AlphabetMixedCase))
                .ToList();

            string consensus = MotifFinder.GenerateConsensus(rows);

            consensus.Should().HaveLength(1, "a one-column matrix yields a length-1 consensus (INV-01)");
            AssertWellFormed(consensus, rows);
            consensus.Should().Be(Oracle(rows), "matches the documented decision rule");
        }
    }

    #endregion

    #region BE — ties (equal-max columns: threshold IUPAC set or alphabetical fallback)

    // Two-base ties over the threshold encode the NC-IUB symbol for the pair (INV-04), and the
    // result is discriminating (a single-base consensus would be wrong here).
    [TestCase("A", "G", 'R')] // {A,G} → R
    [TestCase("C", "T", 'Y')] // {C,T} → Y
    [TestCase("C", "G", 'S')] // {C,G} → S
    [TestCase("A", "T", 'W')] // {A,T} → W
    [TestCase("G", "T", 'K')] // {G,T} → K
    [TestCase("A", "C", 'M')] // {A,C} → M
    public void Generate_TwoBaseTie_EncodesIupacPair(string a, string b, char expected)
    {
        // n=2, θ·n=0.5; each base count 1 > 0.5 ⇒ both enter the set ⇒ the pair's IUPAC symbol.
        MotifFinder.GenerateConsensus(new[] { a, b })
            .Should().Be(expected.ToString(), "an equal-max pair over θ·n encodes its NC-IUB symbol (INV-04)");
    }

    // Three-base ties over the threshold encode the triple's NC-IUB symbol.
    [TestCase("C", "G", "T", 'B')]
    [TestCase("A", "G", "T", 'D')]
    [TestCase("A", "C", "T", 'H')]
    [TestCase("A", "C", "G", 'V')]
    public void Generate_ThreeBaseTie_EncodesIupacTriple(string a, string b, string c, char expected)
    {
        // n=3, θ·n=0.75; each count 1 > 0.75 ⇒ all three enter the set.
        MotifFinder.GenerateConsensus(new[] { a, b, c })
            .Should().Be(expected.ToString(), "an equal-max triple over θ·n encodes its NC-IUB symbol (INV-04)");
    }

    // Four equal bases (each at exactly 25 %): NONE passes the strict θ·n boundary, so the
    // documented fallback emits the most-frequent base with an ALPHABETICAL tie-break → 'A',
    // NEVER 'N' (§5.2, §6.1). This is the canonical non-trivial tie corner of this unit.
    [Test]
    public void Generate_FourEqualBases_FallsBackToA_NotN()
    {
        MotifFinder.GenerateConsensus(new[] { "A", "C", "G", "T" })
            .Should().Be("A", "no base exceeds θ·n=1.0 (strict '>') ⇒ most-frequent alphabetical fallback → A (§5.2)");
    }

    // A base sitting EXACTLY at the threshold is EXCLUDED (strict '>', INV-05): n=4, θ·n=1.0;
    // a base with count 1 (= 25 %) does not pass — only the base over the threshold survives.
    [Test]
    public void Generate_BaseExactlyAtThreshold_Excluded()
    {
        // col0: A=3, C=1 (n=4, θ·n=1.0). A=3>1 passes; C=1 not > 1 ⇒ excluded ⇒ {A} → A.
        MotifFinder.GenerateConsensus(new[] { "A", "A", "A", "C" })
            .Should().Be("A", "C at exactly 25 % is excluded by the strict '>' boundary (INV-05)");
    }

    // Determinism: every permutation of the SAME rows yields the IDENTICAL consensus — the
    // threshold/fallback tie-break is order-independent (INV-04/§5.2; row-order invariant).
    [Test]
    public void Generate_RowOrderIndependent_Deterministic()
    {
        var rng = new Random(171_003);
        for (int trial = 0; trial < 300; trial++)
        {
            int depth = rng.Next(2, 8);
            int width = rng.Next(1, 10);
            var rows = Enumerable.Range(0, depth)
                .Select(_ => RandomRow(rng, width, Alphabet))
                .ToList();

            string baseline = MotifFinder.GenerateConsensus(rows);

            for (int shuffle = 0; shuffle < 4; shuffle++)
            {
                var permuted = rows.OrderBy(_ => rng.Next()).ToList();
                MotifFinder.GenerateConsensus(permuted)
                    .Should().Be(baseline, "consensus is row-order-independent / deterministic");
            }
        }
    }

    #endregion

    #region BE — broad fuzz: random matrices never crash, match the documented rule

    [Test]
    [CancelAfter(30_000)]
    public void Generate_RandomEqualWidthMatrices_NeverThrows_MatchesOracle()
    {
        var rng = new Random(171_004);
        for (int trial = 0; trial < 1500; trial++)
        {
            int depth = rng.Next(1, 10);
            int width = rng.Next(0, 30);
            var rows = Enumerable.Range(0, depth)
                .Select(_ => RandomRow(rng, width, AlphabetMixedCase))
                .ToList();

            string consensus = MotifFinder.GenerateConsensus(rows);

            AssertWellFormed(consensus, rows);
            consensus.Should().Be(Oracle(rows),
                "the unit matches the documented decision rule under fuzzed depth/width/case");
        }
    }

    // Ragged matrices (column count taken from the FIRST row, shorter rows contribute nothing
    // at trailing columns — §5.2): the unit must NOT IndexOutOfRange, and the consensus length
    // stays the first row's length.
    [Test]
    [CancelAfter(30_000)]
    public void Generate_RandomRaggedMatrices_NoIndexOutOfRange_MatchesOracle()
    {
        var rng = new Random(171_005);
        for (int trial = 0; trial < 1000; trial++)
        {
            int depth = rng.Next(1, 8);
            var rows = Enumerable.Range(0, depth)
                .Select(_ => RandomRow(rng, rng.Next(0, 20), Alphabet))
                .ToList();

            Action act = () => MotifFinder.GenerateConsensus(rows);
            act.Should().NotThrow<IndexOutOfRangeException>(
                "shorter rows contribute nothing at trailing columns, never indexed out of range (§5.2)");

            string consensus = MotifFinder.GenerateConsensus(rows);
            AssertWellFormed(consensus, rows);
            consensus.Should().Be(Oracle(rows), "matches the documented ragged-row tally rule");
        }
    }

    #endregion

    #endregion
}
