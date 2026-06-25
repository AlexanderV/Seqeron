using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Analysis area — genomic repeat enumeration (GENOMIC-REPEAT-001), the suffix-tree
/// repeat finder <see cref="GenomicAnalyzer.FindRepeats(DnaSequence, int)"/>: for a DNA sequence it
/// enumerates every DISTINCT substring that occurs at least twice (overlapping occurrences allowed)
/// and has length ≥ a given <c>minLength</c>, each returned as a <see cref="RepeatInfo"/> carrying the
/// repeated substring and its ascending 0-based occurrence positions. The classical suffix-tree
/// characterisation: a substring occurring k times is the path label of an internal node with k leaves
/// below it; a repeat is therefore a NON-EMPTY substring occurring ≥ 2 times.
/// — docs/algorithms/Repeat_Analysis/Repeat_Detection.md §2.1–2.2, §3, §4.1.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate and boundary sequences to the unit and asserts the code NEVER fails in an
/// undisciplined way: no crash / IndexOutOfRange on the empty sequence, a single base, a sequence with
/// NO repeat (all-distinct content), or a sequence that is ENTIRELY one repeat (homopolymer / "ACAC…"),
/// and NO hang / quadratic blow-up on dense-repeat input (guarded with [CancelAfter]). It also asserts
/// the code NEVER produces nonsense: never a length-0 "repeat" (even at minLength ≤ 0 — the empty string
/// is not a repeat, §5.4 #3), never a repeat occurring only once (INV-01: Count ≥ 2), never a repeat
/// shorter than max(1, minLength) (INV-05), never an out-of-bounds or non-start Position (INV-03), never
/// occurrences that are not identical copies of the reported substring, never a duplicated distinct
/// substring, never a non-deterministic result. The subject is a <see cref="DnaSequence"/> (uppercased
/// and validated to {A,C,G,T} at construction), so out-of-domain residues / null bytes / unicode cannot
/// reach the scan; for valid input the method itself throws no validation exception (§3.3). A raw
/// runtime exception, a hang, a length-0 / single-occurrence / out-of-bounds repeat, a missed or
/// fabricated repeat, or a non-deterministic output is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: GENOMIC-REPEAT-001 — repeat enumeration (FindRepeats: all distinct substrings occurring ≥ 2×,
/// length ≥ minLength), the Analysis-area GenomicAnalyzer repeat facade.
/// Checklist: docs/checklists/03_FUZZING.md, row 178 (BE — "no repeat, full repeat, minLen edge").
/// Algorithm doc: docs/algorithms/Repeat_Analysis/Repeat_Detection.md (Test Unit ID GENOMIC-REPEAT-001).
/// Distinct from rows 13–17 (REP-* units on the standalone RepeatFinder class): THIS unit is the
/// suffix-tree repeat finder on <see cref="GenomicAnalyzer"/> in the Analysis area (rows 175–178).
///
/// Fuzz strategy exercised for THIS unit (BE = Boundary Exploitation — граничні значення: 0, -1,
/// MaxInt, empty — docs/checklists/03_FUZZING.md §Description), mapped to the row's three targets:
///   • NO REPEAT — a subject whose content is all-distinct so NO substring of length ≥ minLength recurs
///     (e.g. "ACGT", or any string of distinct characters): empty enumeration, NO crash, NO fabricated
///     repeat (§6.1 "no repeat → empty"). The empty string and a single base are the EMPTY/atom boundary.
///   • FULL REPEAT — a subject that is ENTIRELY a repeat: a homopolymer ("AAAA…") whose maximal repeat
///     is the length-(n-1) prefix with overlapping occurrences {0,1}, or a period-2 tile ("ACAC…").
///     Every shorter qualifying substring is also enumerated. No quadratic hang on the dense case
///     (kept modest, guarded by [CancelAfter]) (§6.1 "overlapping run AAAAAAAAAA → AAAAAAAAA@{0,1}").
///   • minLen EDGE — minLength ∈ {int.MinValue, -1, 0, 1, =length, >length, int.MaxValue}: at the
///     ≤ 0 boundary the effective minimum clamps to 1 so NO length-0 repeat is ever emitted and the
///     scan terminates (§5.4 #3); at minLength = length only full-length repeats survive (overlapping
///     homopolymer corner); at minLength > length the result is empty; no off-by-one at the boundary
///     and no infinite loop / IndexOutOfRange (§6.1 "minLength > all repeats → empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Repeat_Detection.md §2.4 INV-01..06, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   INV-01 every returned repeat occurs ≥ 2 times (Count ≥ 2; corresponds to an internal node);
///   INV-02 RepeatInfo.Length == Sequence.Length (string depth = substring length);
///   INV-03 each Positions[i] is a true 0-based start of Sequence in the text;
///   INV-05 FindRepeats returns only substrings with Length ≥ max(1, minLength) (§5.4 #3);
///   INV-06 Positions is sorted ascending.
///   Empty / no-repeat subject → empty enumeration; no validation exception for valid DnaSequence (§3.3).
///   GenomicAnalyzer.FindRepeats(DnaSequence, int minLength) → IEnumerable&lt;RepeatInfo&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class GenomicRepeatFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>A random ACGT string of the given length (length 0 ⇒ empty string).</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(Alphabet[rng.Next(Alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent oracle for FindRepeats, built directly from the documented model (§2.1–2.2, §3.2),
    /// NOT from the unit: the set of every DISTINCT substring of <paramref name="text"/> that occurs at
    /// least twice (overlapping occurrences counted) and has length ≥ max(1, minLength), keyed by the
    /// substring and mapped to its ascending list of all 0-based start positions. A repeat is a
    /// non-empty substring (the empty string is excluded even when minLength ≤ 0, §5.4 #3).
    /// </summary>
    private static Dictionary<string, List<int>> Oracle(string text, int minLength)
    {
        int effectiveMin = Math.Max(1, minLength);
        var occ = new Dictionary<string, List<int>>(StringComparer.Ordinal);

        for (int len = effectiveMin; len <= text.Length; len++)
        {
            for (int start = 0; start + len <= text.Length; start++)
            {
                string sub = text.Substring(start, len);
                if (!occ.TryGetValue(sub, out var list))
                {
                    list = new List<int>();
                    occ[sub] = list;
                }
                list.Add(start);
            }
        }

        // Keep only substrings occurring >= 2 times (INV-01). Positions are already ascending by
        // construction (outer loop scans start left-to-right).
        return occ.Where(kv => kv.Value.Count >= 2)
                  .ToDictionary(kv => kv.Key, kv => kv.Value, StringComparer.Ordinal);
    }

    /// <summary>
    /// Asserts a single <see cref="RepeatInfo"/> is WELL-FORMED against the documented invariants and the
    /// text it was found in: non-empty substring of length ≥ max(1, minLength) (INV-05, §5.4 #3), Count
    /// ≥ 2 (INV-01), Length == Sequence.Length (INV-02), Positions ascending & distinct (INV-06), and
    /// every position is an in-bounds 0-based start where the text spells exactly the reported substring
    /// (INV-03 — occurrences are identical copies, none fabricated/shifted).
    /// </summary>
    private static void AssertRepeatWellFormed(RepeatInfo r, string text, int minLength)
    {
        int effectiveMin = Math.Max(1, minLength);

        r.IsEmpty.Should().BeFalse("a yielded repeat is never the None/empty value");
        r.Sequence.Should().NotBeNullOrEmpty("a repeat is a NON-EMPTY substring (§5.4 #3)");
        r.Length.Should().Be(r.Sequence.Length, "Length is the substring length (INV-02)");
        r.Length.Should().BeGreaterThanOrEqualTo(effectiveMin,
            "repeat length ≥ max(1, minLength); no length-0 repeat even at minLength ≤ 0 (INV-05, §5.4 #3)");

        r.Count.Should().Be(r.Positions.Count, "Count == Positions.Count");
        r.Count.Should().BeGreaterThanOrEqualTo(2, "a repeat occurs ≥ 2 times (INV-01)");

        r.Positions.Should().BeInAscendingOrder("Positions are sorted ascending (INV-06)");
        r.Positions.Should().OnlyHaveUniqueItems("each occurrence start is distinct");

        foreach (int pos in r.Positions)
        {
            pos.Should().BeGreaterThanOrEqualTo(0, "Position is a 0-based offset (INV-03)");
            (pos + r.Length).Should().BeLessThanOrEqualTo(text.Length,
                "the occurrence is in-bounds (no read past the end) (INV-03)");
            text.Substring(pos, r.Length).Should().Be(r.Sequence,
                "the text spells exactly the reported substring at Position (identical copies, INV-03)");
        }
    }

    /// <summary>
    /// Runs FindRepeats, asserts every yielded repeat is well-formed, asserts the reported (substring →
    /// positions) map EXACTLY equals the independent oracle (no missed repeat, no fabricated repeat, no
    /// duplicate distinct substring, positions match the full overlapping occurrence set), and asserts
    /// determinism (two enumerations agree). Returns the materialised result.
    /// </summary>
    private static List<RepeatInfo> AssertWellFormed(string text, int minLength)
    {
        var dna = new DnaSequence(text);
        var result = GenomicAnalyzer.FindRepeats(dna, minLength).ToList();

        foreach (var r in result)
            AssertRepeatWellFormed(r, dna.Sequence, minLength);

        // No distinct substring is yielded twice.
        result.Select(r => r.Sequence).Should()
            .OnlyHaveUniqueItems("each distinct repeated substring is reported exactly once");

        // Exact set + positions equivalence with the spec oracle.
        var actual = result.ToDictionary(r => r.Sequence, r => r.Positions.ToList(), StringComparer.Ordinal);
        var expected = Oracle(dna.Sequence, minLength);

        actual.Keys.Should().BeEquivalentTo(expected.Keys,
            "the repeat SET must exactly match the spec oracle (none missed, none fabricated)");
        foreach (var kv in expected)
            actual[kv.Key].Should().Equal(kv.Value,
                $"all overlapping 0-based occurrences of '{kv.Key}' must be reported, ascending (INV-03/06)");

        // Determinism: re-enumerate and compare.
        var again = GenomicAnalyzer.FindRepeats(dna, minLength)
            .ToDictionary(r => r.Sequence, r => r.Positions.ToList(), StringComparer.Ordinal);
        again.Keys.Should().BeEquivalentTo(actual.Keys, "enumeration is deterministic");

        return result;
    }

    #endregion

    #region GENOMIC-REPEAT-001 — repeat enumeration (BE: no repeat, full repeat, minLen edge)

    // ─── Positive sanity: a known repeat is found at the documented coordinates ─────────────────

    [Test]
    public void PositiveSanity_KnownRepeat_FoundAtDocumentedCoordinates()
    {
        // Doc §7.1 worked example: ATCGATCGA has the LRS "ATCGA" at {0,4}. FindRepeats(.,5) must therefore
        // report exactly that one length-≥5 repeat (the only substring occurring ≥ 2× with length ≥ 5).
        var dna = new DnaSequence("ATCGATCGA");
        var len5 = GenomicAnalyzer.FindRepeats(dna, minLength: 5).ToList();

        len5.Should().ContainSingle("only 'ATCGA' occurs ≥ 2× with length ≥ 5 (doc §7.1)");
        len5[0].Sequence.Should().Be("ATCGA");
        len5[0].Positions.Should().Equal(0, 4);
        len5[0].Length.Should().Be(5);
        len5[0].Count.Should().Be(2);

        AssertWellFormed("ATCGATCGA", minLength: 5);
        // And full oracle equivalence at a smaller minLength (many shorter repeats too).
        AssertWellFormed("ATCGATCGA", minLength: 1);
    }

    [Test]
    public void PositiveSanity_DocExample_AllRepeatsAtMinLength3()
    {
        // Doc §7.1: FindRepeats("ACGTACGTTTTTACGT", 3) reports exactly these 8 distinct substrings.
        var dna = new DnaSequence("ACGTACGTTTTTACGT");
        var result = AssertWellFormed("ACGTACGTTTTTACGT", minLength: 3);

        var map = result.ToDictionary(r => r.Sequence, r => r.Positions.ToList(), StringComparer.Ordinal);
        map.Should().HaveCount(8, "the documented full set has 8 distinct repeats (doc §7.1)");
        map["ACG"].Should().Equal(0, 4, 12);
        map["ACGT"].Should().Equal(0, 4, 12);
        map["CGT"].Should().Equal(1, 5, 13);
        map["TAC"].Should().Equal(3, 11);
        map["TACG"].Should().Equal(3, 11);
        map["TACGT"].Should().Equal(3, 11);
        map["TTT"].Should().Equal(7, 8, 9);
        map["TTTT"].Should().Equal(7, 8);
    }

    // ─── NO REPEAT: all-distinct content ⇒ empty enumeration, no false repeat ───────────────────

    [Test]
    public void NoRepeat_Empty_NoCrash()
    {
        AssertWellFormed(string.Empty, minLength: 1).Should().BeEmpty("ε has no substring occurring twice (§6.1)");
        AssertWellFormed(string.Empty, minLength: 0).Should().BeEmpty("empty subject, minLength 0 ⇒ still empty");
    }

    [Test]
    public void NoRepeat_SingleBase_NoCrash()
    {
        foreach (var s in new[] { "A", "C", "G", "T" })
            AssertWellFormed(s, minLength: 1).Should().BeEmpty($"'{s}' is a single character; nothing recurs");
    }

    [Test]
    public void NoRepeat_AllDistinctContent_Empty()
    {
        // "ACGT" — every length-≥1 substring is unique (no character repeats), so no repeat at all (§6.1).
        AssertWellFormed("ACGT", minLength: 1).Should().BeEmpty("ACGT has no repeated substring (§6.1)");
        // "ACGTACG" repeats single chars but not for minLength large enough; at minLength 4 nothing recurs.
        AssertWellFormed("ACGTACG", minLength: 4).Should().BeEmpty("no substring of length ≥ 4 recurs in ACGTACG");
    }

    [Test]
    public void NoRepeat_RandomLongDistinctRun_NoForwardRepeatAtThatLength()
    {
        // A random sequence's LRS length grows ~logarithmically; for length n, minLength = n (the whole
        // string) can never recur, so FindRepeats(.,n) is always empty — a clean "no repeat" boundary.
        var rng = new Random(178_001);
        for (int t = 0; t < 300; t++)
        {
            string s = RandomDna(rng, rng.Next(1, 40));
            AssertWellFormed(s, minLength: s.Length).Should()
                .BeEmpty("a string is never a repeated substring of itself (length-n substring occurs once)");
        }
    }

    // ─── FULL REPEAT: the whole sequence is one repeat (homopolymer / period-2 tile) ────────────

    [Test]
    public void FullRepeat_Homopolymer_MaximalRepeatWithOverlappingOccurrences()
    {
        // "AAAAAAAAAA" (n=10): per doc §6.1 the LRS is "AAAAAAAAA" (length 9) at positions {0,1}
        // (overlap allowed). FindRepeats(.,9) must report exactly that single maximal repeat.
        var dna = new DnaSequence("AAAAAAAAAA");
        var top = GenomicAnalyzer.FindRepeats(dna, minLength: 9).ToList();

        top.Should().ContainSingle("only the length-9 prefix recurs (overlapping) (§6.1)");
        top[0].Sequence.Should().Be("AAAAAAAAA");
        top[0].Positions.Should().Equal(0, 1);
        top[0].Length.Should().Be(9);

        // And the full enumeration matches the oracle: every length 1..9 'A'-run recurs (overlapping).
        var all = AssertWellFormed("AAAAAAAAAA", minLength: 1);
        all.Should().HaveCount(9, "lengths 1..9 each give one distinct 'A'-run repeat");
    }

    [Test]
    public void FullRepeat_Period2Tile_AllPrefixesRecur()
    {
        // "ACACACAC" (n=8) is entirely a repeat of "AC". The maximal repeated substrings are the
        // length-(n-2) windows; the full enumeration is verified against the oracle (no hang, no
        // off-by-one). The unit "AC" itself recurs at {0,2,4,6}.
        var result = AssertWellFormed("ACACACAC", minLength: 1);
        var map = result.ToDictionary(r => r.Sequence, r => r.Positions.ToList(), StringComparer.Ordinal);
        map.Should().ContainKey("AC");
        map["AC"].Should().Equal(new[] { 0, 2, 4, 6 }, "the period-2 unit recurs at every even offset (overlap allowed)");
    }

    [Test]
    [CancelAfter(60_000)]
    public void FullRepeat_DenseHomopolymer_NoQuadraticHang()
    {
        // Dense full-repeat input is the worst case for the O(n²) enumeration (§4.3). Keep modest and
        // guard with [CancelAfter]: a long homopolymer must complete and stay well-formed, not hang.
        foreach (int n in new[] { 50, 120, 200 })
        {
            string s = new string('A', n);
            var result = AssertWellFormed(s, minLength: Math.Max(1, n - 5));
            // The top minLength window (n - 5 .. n - 1) gives a handful of overlapping 'A'-runs.
            result.Should().NotBeEmpty($"a length-{n} homopolymer has long repeated runs");
            result.Should().OnlyContain(r => r.Sequence.All(c => c == 'A'),
                "every repeat in a homopolymer is an 'A'-run");
        }
    }

    [Test]
    [CancelAfter(60_000)]
    public void FullRepeat_DenseTandemUnit_NoQuadraticHang()
    {
        // A long "ACGT"-tiled sequence is also dense in repeats; verify it completes and matches oracle.
        foreach (int copies in new[] { 20, 40 })
        {
            string s = string.Concat(Enumerable.Repeat("ACGT", copies));
            // Use a high minLength so the oracle stays small but the unit must still scan the dense text.
            AssertWellFormed(s, minLength: s.Length - 6);
        }
    }

    // ─── minLen EDGE: 0 / 1 / =length / >length / extreme ints (§5.4 #3, §6.1) ──────────────────

    [Test]
    public void MinLen_ZeroAndNegative_ClampToOne_NoLengthZeroRepeat_NoInfiniteLoop()
    {
        // minLength ≤ 0 must clamp to an effective minimum of 1: NEVER a length-0 "repeat", NEVER an
        // infinite loop (§5.4 #3). A subject with a known repeat at minLength {int.MinValue, -1, 0, 1}
        // must yield the SAME well-formed result as minLength 1, all with Length ≥ 1.
        string text = "ACGTACGT"; // "ACGT" recurs at {0,4}; "ACG"@{0,4}; "CGT"@{1,5}; single chars too
        var baseline = AssertWellFormed(text, minLength: 1)
            .ToDictionary(r => r.Sequence, r => r.Positions.ToList(), StringComparer.Ordinal);

        foreach (int min in new[] { int.MinValue, -1, 0 })
        {
            var result = AssertWellFormed(text, minLength: min);
            result.Should().OnlyContain(r => r.Length >= 1, "no length-0 repeat at minLength ≤ 0 (§5.4 #3)");
            result.ToDictionary(r => r.Sequence, r => r.Positions.ToList(), StringComparer.Ordinal)
                .Should().BeEquivalentTo(baseline, "minLength ≤ 0 behaves as minLength 1 (clamp, §5.4 #3)");
        }
    }

    [Test]
    public void MinLen_EqualsSequenceLength_OnlyFullLengthRepeatsSurvive()
    {
        // minLength = n: only a length-n substring could qualify, but a length-n substring occurs once,
        // so a non-repeating subject yields empty (off-by-one at the boundary check).
        AssertWellFormed("ACGTACGT", minLength: 8).Should()
            .BeEmpty("the length-8 substring is the whole string; it occurs once (no length-n repeat)");

        // Overlapping homopolymer corner: for "AAAA" (n=4), minLength = 4 ⇒ "AAAA" occurs once ⇒ empty,
        // but minLength = 3 ⇒ "AAA"@{0,1} recurs (overlap). Exercises the =length / =length-1 boundary.
        AssertWellFormed("AAAA", minLength: 4).Should().BeEmpty("'AAAA' length 4 occurs once");
        var three = AssertWellFormed("AAAA", minLength: 3);
        three.Should().ContainSingle();
        three[0].Sequence.Should().Be("AAA");
        three[0].Positions.Should().Equal(0, 1);
    }

    [Test]
    public void MinLen_GreaterThanSequenceLength_Empty_NoCrash()
    {
        foreach (var text in new[] { "", "A", "ACGTACGT", "AAAAAAAA" })
        {
            int n = text.Length;
            AssertWellFormed(text, minLength: n + 1).Should()
                .BeEmpty($"minLength {n + 1} exceeds |'{text}'|; no substring that long exists (§6.1)");
            GenomicAnalyzer.FindRepeats(new DnaSequence(text), int.MaxValue).ToList()
                .Should().BeEmpty("int.MaxValue minLength filters everything, no crash/overflow");
        }
    }

    [Test]
    public void MinLen_BoundaryValues_NoCrash_RespectFilter()
    {
        // A single subject swept across the boundary minLength values; each result well-formed and
        // monotonically shrinking as minLength rises (a higher floor only removes shorter repeats).
        string text = "ACGTACGTACGT";
        int prev = int.MaxValue;
        foreach (int min in new[] { int.MinValue, -1, 0, 1, 2, 3, 4, 8, 11, 12, 13 })
        {
            var result = AssertWellFormed(text, min);
            int count = result.Count;
            if (min > 1) // for min ≤ 1 the effective floor is 1 (clamp); compare only the rising tail.
                count.Should().BeLessThanOrEqualTo(prev, "raising minLength never adds repeats (monotone filter)");
            prev = count;
        }
    }

    // ─── Broad random fuzz: well-formedness + oracle equivalence over many shapes ────────────────

    [Test]
    [CancelAfter(90_000)]
    public void RandomSequences_WellFormed_AndMatchOracle()
    {
        var rng = new Random(178_002);
        for (int t = 0; t < 1200; t++)
        {
            string text = RandomDna(rng, rng.Next(0, 60));
            int bucket = t % 5;
            int min = bucket switch
            {
                0 => 1,
                1 => 2,
                2 => 3,
                3 => rng.Next(-3, 10),
                _ => rng.Next(0, Math.Max(1, text.Length) + 2),
            };
            AssertWellFormed(text, min);
        }
    }

    [Test]
    [CancelAfter(90_000)]
    public void RandomSeededRepeats_AlwaysFound_AndWellFormed()
    {
        // Embed a guaranteed repeated unit twice (with random filler between) so the result is non-empty
        // often, then verify full well-formedness and oracle equivalence. The embedded unit (or a
        // qualifying prefix of it) must appear among the reported repeats.
        var rng = new Random(178_003);
        for (int t = 0; t < 600; t++)
        {
            string unit = RandomDna(rng, rng.Next(2, 6));
            string filler = RandomDna(rng, rng.Next(0, 8));
            string text = unit + filler + unit;

            var result = AssertWellFormed(text, minLength: unit.Length);
            // The embedded unit occurs ≥ 2 times and has length == minLength, so it MUST be reported.
            result.Select(r => r.Sequence).Should()
                .Contain(unit, "an embedded length-≥minLength substring occurring twice is reported (INV-01)");
        }
    }

    #endregion
}
