using System.Text;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Matching area — Consensus from a Multiple Alignment (MOTIF-CONS-001),
/// the column-wise most-frequent ("majority") consensus caller
/// <see cref="MotifFinder.CreateConsensusFromAlignment(IEnumerable{string})"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop
/// (the O(n×m) per-column tally must always terminate), no state corruption, no
/// nonsense output (a consensus whose length differs from the alignment width, an
/// emitted symbol outside the documented alphabet {A,C,G,T}, a non-deterministic
/// tie winner), and no *unhandled* runtime exception — in particular NO
/// IndexOutOfRange when rows have unequal lengths, and NO DivideByZero /
/// NullReference on an empty alignment. Every input must resolve to EITHER a
/// well-defined, theory-correct result OR a *documented, intentional* validation
/// exception (ArgumentNullException for a null collection; ArgumentException for
/// unequal-length rows or a non-ACGT character — contract §3.3, §6.1). A raw
/// runtime exception, a hang, a wrong-length consensus, an out-of-alphabet symbol,
/// or an order-dependent tie winner is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MOTIF-CONS-001 — Consensus from a Multiple Alignment
/// Checklist: docs/checklists/03_FUZZING.md, row 169.
/// Algorithm doc: docs/algorithms/Pattern_Matching/Consensus_From_Alignment.md
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the row:
///          – UNEQUAL ROW LENGTHS: rows of differing lengths are NOT silently indexed
///            (no IndexOutOfRange at a high column over a short row); the documented
///            contract is a thrown ArgumentException (§3.3, §6.1).
///          – EMPTY ALIGNMENT: zero rows → the documented empty string "" with NO
///            DivideByZero / NullReference on a zero-row tally (§3.3, §6.1).
///          – SINGLE ROW: one row → the consensus equals that row exactly (uppercased),
///            each column's only base being its own maximum, no crash (§6.1, INV-04).
///   • MC = Malformed Content — out-of-domain residues: any character outside
///          {A,C,G,T} after uppercasing → the documented ArgumentException, never a
///          silent emit of garbage and never an out-of-alphabet output symbol (§3.3).
/// — docs/checklists/03_FUZZING.md §Description (BE = граничні значення 0/-1/MaxInt/empty;
///   MC = Malformed Content / невалідний контент).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Consensus_From_Alignment.md §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// For equal-length aligned DNA strings, build the 4×n profile (A,C,G,T counts per
/// column) and emit, per column, the base with the maximum count (INV-02), scanning
/// the alphabet A→C→G→T with a strict '>' so ties resolve to the alphabetically
/// EARLIEST base (INV-03). Inputs are uppercased (case-insensitive). The consensus
/// length equals the common row width (INV-01); identical/single rows reproduce that
/// row (INV-04); the result is fully deterministic (INV-05). Validation (§3.3):
/// null collection → ArgumentNullException; empty collection → ""; unequal-length
/// rows → ArgumentException; any non-ACGT character (after uppercasing) →
/// ArgumentException.
///   MotifFinder.CreateConsensusFromAlignment(IEnumerable&lt;string&gt;) → string
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MotifConsensusFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>
    /// Independent oracle implementing the documented decision rule verbatim
    /// (Consensus_From_Alignment.md §4.1): per column, tally A,C,G,T over the
    /// uppercased rows and emit the base with the maximum count, scanning A→C→G→T so
    /// ties resolve to the alphabetically-earliest base. Assumes pre-validated,
    /// equal-length, all-ACGT rows. Used to cross-check the unit on fuzzed inputs
    /// without re-using its implementation.
    /// </summary>
    private static string Oracle(IReadOnlyList<string> rows)
    {
        if (rows.Count == 0) return "";
        int width = rows[0].Length;
        var sb = new StringBuilder(width);

        for (int col = 0; col < width; col++)
        {
            var counts = new int[Alphabet.Length];
            foreach (string r in rows)
                counts[Array.IndexOf(Alphabet, char.ToUpperInvariant(r[col]))]++;

            int best = 0;
            for (int b = 1; b < counts.Length; b++)
                if (counts[b] > counts[best]) best = b; // strict '>' ⇒ alphabetically-earliest tie winner
            sb.Append(Alphabet[best]);
        }

        return sb.ToString();
    }

    /// <summary>
    /// Asserts a consensus is WELL-FORMED per the documented contract regardless of
    /// the (possibly degenerate) equal-length input:
    ///   INV-01 length == common row width (== 0 for an empty alignment);
    ///   output alphabet ⊆ {A,C,G,T} (§3.2, §6.2) — never an out-of-alphabet symbol;
    ///   INV-02 each emitted base actually ATTAINS the maximum count in its column.
    /// </summary>
    private static void AssertWellFormed(string consensus, IReadOnlyList<string> rows)
    {
        int width = rows.Count == 0 ? 0 : rows[0].Length;
        consensus.Should().HaveLength(width, "INV-01: consensus length == alignment width");

        foreach (char c in consensus)
            Alphabet.Should().Contain(c, "every output symbol is from the documented alphabet {A,C,G,T} (§3.2)");

        for (int col = 0; col < consensus.Length; col++)
        {
            var counts = new int[Alphabet.Length];
            foreach (string r in rows)
                counts[Array.IndexOf(Alphabet, char.ToUpperInvariant(r[col]))]++;

            int max = counts.Max();
            int idx = Array.IndexOf(Alphabet, consensus[col]);
            counts[idx].Should().Be(max, "INV-02: the emitted base attains the column maximum count");
        }
    }

    /// <summary>An aligned row of the given width over a chosen alphabet.</summary>
    private static string RandomRow(Random rng, int width, char[] alphabet)
    {
        var sb = new StringBuilder(width);
        for (int i = 0; i < width; i++)
            sb.Append(alphabet[rng.Next(alphabet.Length)]);
        return sb.ToString();
    }

    private static readonly char[] AlphabetMixedCase = { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't' };

    #endregion

    #region MOTIF-CONS-001 — Consensus from a Multiple Alignment (BE: unequal rows, empty, single; MC: invalid chars)

    #region Positive sanity — hand-computed documented consensus

    // Documented Rosalind CONS worked example (§7.1, §136): profile column maxima
    //   A=5 1 0 0 5 5 0 0 / C=0 0 1 4 2 0 6 1 / G=1 1 6 3 0 1 0 0 / T=1 5 0 0 0 1 1 6
    //   → 5→A 5→T 6→G 4→C 5→A 5→A 6→C 6→T  ⇒  "ATGCAACT".
    [Test]
    public void CreateConsensus_RosalindSample_DocumentedConsensus()
    {
        var rows = new[]
        {
            "ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC",
            "TTGGAACT", "ATGCCATT", "ATGGCACT",
        };

        string consensus = MotifFinder.CreateConsensusFromAlignment(rows);

        consensus.Should().Be("ATGCAACT", "documented Rosalind CONS column maxima (§7.1)");
        AssertWellFormed(consensus, rows);
    }

    // A clear per-column majority resolves to the strict-max base (deterministic).
    [Test]
    public void CreateConsensus_ClearMajority_ResolvesToStrictMaxPerColumn()
    {
        // col0 A=3,T=1→A; col1 C=4→C; col2 G=3,A=1→G; col3 T=4→T.
        var rows = new[] { "ACGT", "ACGT", "ACAT", "TCGT" };

        MotifFinder.CreateConsensusFromAlignment(rows)
            .Should().Be("ACGT", "each column's unique strict-max base is committed (INV-02)");
    }

    // Tie column resolves to the alphabetically-EARLIEST base (INV-03), and this is
    // discriminating: a reversed tie-break would return the later base.
    [Test]
    public void CreateConsensus_TieColumn_ResolvesAlphabeticallyEarliest()
    {
        // col0 {A,G} tie 1/1 → A (not G); col1 {T,T} → T.   col0 {G,T} tie → G (not T).
        MotifFinder.CreateConsensusFromAlignment(new[] { "AT", "GT" })
            .Should().Be("AT", "tie (A,G) resolves to alphabetically-earliest A (INV-03)");
        MotifFinder.CreateConsensusFromAlignment(new[] { "GT", "TT" })
            .Should().Be("GT", "tie (G,T) resolves to alphabetically-earliest G (INV-03)");
    }

    #endregion

    #region BE — Boundary: single row (consensus == that row, uppercased, no crash)

    // One row → consensus equals that row exactly (each column's only base is its max).
    [Test]
    public void CreateConsensus_SingleRow_EqualsThatRow()
    {
        var rows = new[] { "ACGTACGTAC" };

        string consensus = MotifFinder.CreateConsensusFromAlignment(rows);

        consensus.Should().Be("ACGTACGTAC", "a single row commits its own base at every column (INV-04)");
        AssertWellFormed(consensus, rows);
    }

    // Single lowercase row → uppercased consensus equal to the row (case-insensitive, §3.3).
    [Test]
    public void CreateConsensus_SingleLowercaseRow_UppercasedRow()
    {
        MotifFinder.CreateConsensusFromAlignment(new[] { "acgtac" })
            .Should().Be("ACGTAC", "input is uppercased before processing (case-insensitive contract)");
    }

    // Fuzz: a single random row always yields itself (uppercased), never throws.
    [Test]
    public void CreateConsensus_SingleRandomRow_NeverThrows_EqualsUppercasedRow()
    {
        var rng = new Random(169_001);
        for (int trial = 0; trial < 400; trial++)
        {
            int width = rng.Next(0, 40);
            string row = RandomRow(rng, width, AlphabetMixedCase);
            var rows = new[] { row };

            string consensus = MotifFinder.CreateConsensusFromAlignment(rows);

            consensus.Should().Be(row.ToUpperInvariant(), "single-row consensus == uppercased row (INV-04)");
            AssertWellFormed(consensus, rows);
            consensus.Should().Be(Oracle(rows), "matches the documented decision rule");
        }
    }

    #endregion

    #region BE — Boundary: empty alignment (no rows ⇒ "", no DivideByZero / NullReference)

    // No rows → empty string, no columns, NO DivideByZero / NullReference (§3.3, §6.1).
    [Test]
    public void CreateConsensus_EmptyAlignment_ReturnsEmptyString()
    {
        MotifFinder.CreateConsensusFromAlignment(Array.Empty<string>())
            .Should().BeEmpty("an empty collection has no columns → \"\" (§3.3, §6.1)");
    }

    // Rows that are all empty strings → common width 0 → empty consensus, no crash.
    [Test]
    public void CreateConsensus_AllEmptyRows_ReturnsEmptyString()
    {
        MotifFinder.CreateConsensusFromAlignment(new[] { "", "", "" })
            .Should().BeEmpty("common width 0 ⇒ empty consensus, no zero-row division");
    }

    // Null collection → documented ArgumentNullException (§3.3, §6.1).
    [Test]
    public void CreateConsensus_NullCollection_Throws()
    {
        Action act = () => MotifFinder.CreateConsensusFromAlignment(null!);

        act.Should().Throw<ArgumentNullException>("null collection is the documented validation contract (§3.3)");
    }

    #endregion

    #region BE — Boundary: unequal row lengths (documented ArgumentException, no IndexOutOfRange)

    // Rows of differing lengths → documented ArgumentException, NOT IndexOutOfRange when a
    // high column is indexed over a short row (§3.3, §6.1).
    [Test]
    public void CreateConsensus_UnequalRowLengths_ThrowsArgumentException()
    {
        var rows = new[] { "ACGT", "AC", "ACGTACGT" };

        Action act = () => MotifFinder.CreateConsensusFromAlignment(rows);

        act.Should().Throw<ArgumentException>(
            "unequal-length rows are rejected by validation, never indexed out of range (§3.3)")
            .Which.Message.Should().Contain("length");
    }

    // A short row AFTER a long first row: the failure must surface as ArgumentException
    // (the validation guard), never an IndexOutOfRange while tallying high columns.
    [Test]
    public void CreateConsensus_LongFirstThenShortRow_ThrowsArgumentExceptionNotIndexOutOfRange()
    {
        var rows = new[] { "ACGTACGT", "AC" };

        Action act = () => MotifFinder.CreateConsensusFromAlignment(rows);

        act.Should().Throw<ArgumentException>(
            "the validation guard fires before any high-column indexing of the short row");
        act.Should().NotThrow<IndexOutOfRangeException>();
    }

    // Fuzz: any ragged row set throws ArgumentException — never IndexOutOfRange, never garbage.
    [Test]
    [CancelAfter(30_000)]
    public void CreateConsensus_RandomRaggedRows_AlwaysThrowsArgumentException_NoIndexOutOfRange()
    {
        var rng = new Random(169_002);
        for (int trial = 0; trial < 600; trial++)
        {
            int depth = rng.Next(2, 8);
            var rows = new List<string>(depth);
            for (int i = 0; i < depth; i++)
                rows.Add(RandomRow(rng, rng.Next(0, 25), Alphabet));

            // Force at least one length mismatch.
            if (rows.Select(r => r.Length).Distinct().Count() == 1)
                rows[rng.Next(depth)] += "A";

            Action act = () => MotifFinder.CreateConsensusFromAlignment(rows);

            act.Should().Throw<ArgumentException>(
                "ragged rows are rejected, never indexed out of range");
            act.Should().NotThrow<IndexOutOfRangeException>();
        }
    }

    #endregion

    #region MC — Malformed Content: non-ACGT residues (documented ArgumentException)

    // A non-ACGT character (gap, IUPAC code, digit, space) → documented ArgumentException,
    // never a silent emit or an out-of-alphabet output symbol (§3.3, §6.2).
    [TestCase("ACGT", "ANGT")]   // 'N' IUPAC code not accepted by this most-frequent variant
    [TestCase("ACGT", "AC-T")]   // gap
    [TestCase("ACGT", "AC1T")]   // digit
    [TestCase("ACGT", "AC T")]   // space
    [TestCase("ACGT", "ACGU")]   // RNA uracil
    public void CreateConsensus_NonAcgtCharacter_ThrowsArgumentException(string ok, string bad)
    {
        var rows = new[] { ok, bad };

        Action act = () => MotifFinder.CreateConsensusFromAlignment(rows);

        act.Should().Throw<ArgumentException>(
            "any character outside {A,C,G,T} after uppercasing is rejected (§3.3)");
    }

    // Lowercase ACGT is VALID (case-insensitive), so it must NOT be treated as malformed.
    [Test]
    public void CreateConsensus_LowercaseAcgt_NotMalformed()
    {
        MotifFinder.CreateConsensusFromAlignment(new[] { "acgt", "ACGT" })
            .Should().Be("ACGT", "lowercase ACGT is uppercased, not rejected (§3.3)");
    }

    // Fuzz: injecting a single out-of-domain character into otherwise-valid equal-length
    // rows always throws ArgumentException — never returns an out-of-alphabet symbol.
    [Test]
    [CancelAfter(30_000)]
    public void CreateConsensus_RandomInvalidChar_AlwaysThrowsArgumentException()
    {
        const string invalidPool = "NRYSWKMBDHVU-.* 0123456789xz";
        var rng = new Random(169_003);
        for (int trial = 0; trial < 600; trial++)
        {
            int depth = rng.Next(1, 6);
            int width = rng.Next(1, 20);
            var rows = Enumerable.Range(0, depth)
                .Select(_ => RandomRow(rng, width, Alphabet))
                .ToArray();

            // Corrupt one position in one row with a guaranteed-invalid character.
            int badRow = rng.Next(depth);
            int badCol = rng.Next(width);
            char bad = invalidPool[rng.Next(invalidPool.Length)];
            var chars = rows[badRow].ToCharArray();
            chars[badCol] = bad;
            rows[badRow] = new string(chars);

            Action act = () => MotifFinder.CreateConsensusFromAlignment(rows);

            act.Should().Throw<ArgumentException>(
                $"row {badRow} col {badCol} contains invalid '{bad}' (§3.3)");
        }
    }

    #endregion

    #region BE — Broad fuzz: random equal-length valid alignments never crash, match the documented rule

    [Test]
    [CancelAfter(30_000)]
    public void CreateConsensus_RandomValidAlignments_NeverThrows_MatchesOracle()
    {
        var rng = new Random(169_004);
        for (int trial = 0; trial < 1500; trial++)
        {
            int depth = rng.Next(1, 9);
            int width = rng.Next(0, 30);
            var rows = Enumerable.Range(0, depth)
                .Select(_ => RandomRow(rng, width, AlphabetMixedCase))
                .ToList();

            string consensus = MotifFinder.CreateConsensusFromAlignment(rows);

            AssertWellFormed(consensus, rows);
            consensus.Should().Be(Oracle(rows),
                "the unit matches the documented decision rule under fuzzed depth/width/case");
        }
    }

    // Determinism / row-order independence: every permutation of the SAME rows yields the
    // IDENTICAL consensus (INV-05; alphabetical tie-break is not order-dependent).
    [Test]
    public void CreateConsensus_RowOrderIndependent_Deterministic()
    {
        var rng = new Random(169_005);
        for (int trial = 0; trial < 300; trial++)
        {
            int depth = rng.Next(2, 7);
            int width = rng.Next(1, 10);
            var rows = Enumerable.Range(0, depth)
                .Select(_ => RandomRow(rng, width, Alphabet))
                .ToList();

            string baseline = MotifFinder.CreateConsensusFromAlignment(rows);

            for (int shuffle = 0; shuffle < 4; shuffle++)
            {
                var permuted = rows.OrderBy(_ => rng.Next()).ToList();
                MotifFinder.CreateConsensusFromAlignment(permuted)
                    .Should().Be(baseline, "consensus is order-independent / deterministic (INV-05)");
            }
        }
    }

    #endregion

    #endregion
}
