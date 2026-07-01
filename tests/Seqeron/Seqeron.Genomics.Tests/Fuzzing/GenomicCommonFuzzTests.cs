using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Analysis area — Common Region Detection (GENOMIC-COMMON-001),
/// the exact, deterministic longest-common-substring (LCS) facade
/// <see cref="GenomicAnalyzer.FindLongestCommonRegion(DnaSequence, DnaSequence)"/> and its
/// per-start-position helper <see cref="GenomicAnalyzer.FindCommonRegions(DnaSequence, DnaSequence, int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain sequence pairs / parameters to the
/// unit and asserts the code NEVER fails in an undisciplined way: no hang/infinite loop (the
/// O(n+m) suffix-tree streaming and the O(n + m·log m) per-start binary search must always
/// terminate), no state corruption, no nonsense output (a reported region that is NOT a
/// contiguous substring of BOTH inputs, a region whose reported 0-based positions don't spell
/// it, a "common" region fabricated when nothing is shared, a region shorter than max(1,
/// minLength), a NON-DETERMINISTIC result), no NullReference/crash on a single-char or empty
/// pair, and no DivideByZero anywhere (there are no ratios in this unit — but the binary-search
/// bounds must not underflow/overflow on the empty / single-char boundary). Every input must
/// resolve to a well-defined, theory-correct result; this unit declares NO validation exception
/// (empty input is accepted and yields CommonRegion.None / an empty enumeration — §3.3, §6.1).
/// A raw runtime exception, a hang, a false common region, a wrong position, a region longer
/// than the true LCS, or a non-deterministic output is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: GENOMIC-COMMON-001 — Common Region Detection (longest common SUBSTRING, contiguous)
/// Checklist: docs/checklists/03_FUZZING.md, row 175.
/// Algorithm doc: docs/algorithms/Sequence_Comparison/Common_Region_Detection.md
/// Distinct from row 164 PROTMOTIF-COMMON-001 (ProteinMotifFinder, protein motif collection)
/// and row 173 MOTIF-SHARED-001 (MotifFinder multi-sequence quorum word enumeration): THIS unit
/// is the two-sequence DNA longest-common-substring facade on GenomicAnalyzer.
///
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the row:
///       – single input: the minimal pair shapes — one sequence empty, or a single-character
///         sequence on one/both sides. Empty pair → CommonRegion.None / empty enumeration, NO
///         NullReference/crash (§6.1, INV-05). A single shared base → length-1 region at valid
///         0-based positions (S4, INV-01).
///       – disjoint: sequences sharing NO contiguous character (e.g. all-A vs all-C) → an
///         EMPTY result, NO false common region fabricated (INV-05, §6.1).
///       – identical: seq1 == seq2 → the WHOLE sequence at positions 0/0 (a string is a
///         substring of itself, §6.1); the result is DETERMINISTIC across repeated calls
///         (INV-01, INV-02).
///       – minLength < 1 (0, −1, int.MinValue): treated as 1 — no empty / length-0 regions
///         (§6.1, INV-04). MaxInt minLength → no region long enough → empty, no crash.
/// — docs/checklists/03_FUZZING.md §Description (BE = граничні значення: 0, -1, MaxInt, empty).
///
/// Note on Malformed Content / Injection: each operand is a <see cref="DnaSequence"/>, which is
/// uppercased and validated to the {A,C,G,T} alphabet at construction, so out-of-domain residues,
/// null bytes and unicode cannot reach this method; this is therefore a pure boundary (BE) row
/// over the pair shape (single / disjoint / identical / empty) and the integer minLength, exactly
/// as the checklist row specifies.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Common_Region_Detection.md §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// Given two DNA sequences, the longest common SUBSTRING is the longest CONTIGUOUS string that
/// is a substring of both (NOT a gapped subsequence — §2.1). FindLongestCommonRegion returns a
/// CommonRegion(Sequence, PositionInFirst, PositionInSecond) where:
///   INV-01 the returned substring is contiguous and occurs in both at the reported 0-based
///          positions: sequence1[PositionInFirst .. +Length] == sequence2[PositionInSecond ..
///          +Length] == Sequence;
///   INV-02 no common contiguous substring strictly longer than the returned one exists;
///   INV-03 on a length tie the substring first found in sequence2 is returned (deterministic);
///   INV-05 empty input or no shared character → CommonRegion.None (empty, length 0, positions −1).
/// FindCommonRegions(minLength): for each start position in sequence2, the SINGLE longest common
/// substring of length ≥ max(1,minLength) that also occurs in sequence1, distinct substrings
/// reported once with the first occurrence in sequence1 and the start in sequence2 (INV-04 —
/// right-maximal-per-start, deduplicated, NOT every common substring; minLength < 1 ≡ 1).
///   GenomicAnalyzer.FindLongestCommonRegion(DnaSequence, DnaSequence) → CommonRegion
///   GenomicAnalyzer.FindCommonRegions(DnaSequence, DnaSequence, int minLength) → IEnumerable&lt;CommonRegion&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class GenomicCommonFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>A random ACGT string of the given length.</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(Alphabet[rng.Next(Alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent O(n·m) dynamic-programming oracle for the longest common SUBSTRING length of
    /// two strings (the textbook contiguous LCS recurrence, §2.2). Returns the maximal length and
    /// the full set of distinct maximal-length common substrings — the unit must return ONE of
    /// these (its documented first-in-sequence2 tie-break, INV-03), and that representative must
    /// be exactly this length (INV-02). Built from the spec, not the unit.
    /// </summary>
    private static (int Length, HashSet<string> Candidates) OracleLcs(string a, string b)
    {
        var candidates = new HashSet<string>();
        int best = 0;
        if (a.Length == 0 || b.Length == 0)
            return (0, candidates);

        // dp[i,j] = length of the common substring ending at a[i-1], b[j-1].
        var prev = new int[b.Length + 1];
        var cur = new int[b.Length + 1];
        for (int i = 1; i <= a.Length; i++)
        {
            for (int j = 1; j <= b.Length; j++)
            {
                if (a[i - 1] == b[j - 1])
                {
                    cur[j] = prev[j - 1] + 1;
                    if (cur[j] > best)
                    {
                        best = cur[j];
                        candidates.Clear();
                        candidates.Add(a.Substring(i - cur[j], cur[j]));
                    }
                    else if (cur[j] == best)
                    {
                        candidates.Add(a.Substring(i - cur[j], cur[j]));
                    }
                }
                else
                {
                    cur[j] = 0;
                }
            }
            (prev, cur) = (cur, prev);
            Array.Clear(cur, 0, cur.Length);
        }
        return (best, candidates);
    }

    /// <summary>
    /// Asserts a FindLongestCommonRegion result is WELL-FORMED per the documented contract,
    /// cross-checked against the independent DP oracle:
    ///   • when the oracle LCS length is 0 (empty input / no shared char) the region is
    ///     CommonRegion.None — empty, length 0, positions −1 (INV-05);
    ///   • otherwise the reported Length equals the oracle's maximal length (INV-02 — never
    ///     longer, never shorter than the true LCS);
    ///   • the reported Sequence is one of the oracle's maximal-length candidates;
    ///   • the reported substring GENUINELY occurs at PositionInFirst in seq1 and
    ///     PositionInSecond in seq2 — contiguous, exact (INV-01).
    /// </summary>
    private static void AssertWellFormedRegion(CommonRegion region, string a, string b)
    {
        var (oracleLen, candidates) = OracleLcs(a, b);

        if (oracleLen == 0)
        {
            region.IsEmpty.Should().BeTrue("no shared character ⇒ CommonRegion.None (INV-05)");
            region.Sequence.Should().BeEmpty();
            region.Length.Should().Be(0);
            region.PositionInFirst.Should().Be(-1, "None reports position −1 (§3.2, INV-05)");
            region.PositionInSecond.Should().Be(-1);
            return;
        }

        region.IsEmpty.Should().BeFalse("a shared substring of length {0} exists (INV-05)", oracleLen);
        region.Length.Should().Be(oracleLen, "reported length is the maximal common-substring length (INV-02)");
        region.Sequence.Should().HaveLength(oracleLen, "Length == |Sequence|");
        candidates.Should().Contain(region.Sequence,
            "the reported substring is a maximal-length common substring (INV-02/INV-03)");

        // INV-01: positions are valid and spell the substring in BOTH sequences.
        region.PositionInFirst.Should().BeInRange(0, a.Length - region.Length,
            "PositionInFirst is a valid 0-based start in sequence1 (INV-01)");
        region.PositionInSecond.Should().BeInRange(0, b.Length - region.Length,
            "PositionInSecond is a valid 0-based start in sequence2 (INV-01)");
        a.Substring(region.PositionInFirst, region.Length).Should().Be(region.Sequence,
            "the substring genuinely occurs at PositionInFirst in sequence1 (INV-01)");
        b.Substring(region.PositionInSecond, region.Length).Should().Be(region.Sequence,
            "the substring genuinely occurs at PositionInSecond in sequence2 (INV-01)");
    }

    /// <summary>
    /// Asserts a FindCommonRegions result is WELL-FORMED per the documented contract (INV-04):
    ///   • distinct Sequence values (one record per substring);
    ///   • every region has length ≥ max(1, minLength) and is non-empty (no length-0 region);
    ///   • every region GENUINELY occurs (contiguous, exact) at its reported PositionInFirst in
    ///     sequence1 and PositionInSecond (its start) in sequence2;
    ///   • the start position in sequence2 indeed begins that region (right-maximal-per-start);
    ///   • each reported substring is a true common contiguous substring of both inputs.
    /// </summary>
    private static void AssertWellFormedRegions(
        IReadOnlyList<CommonRegion> regions, string a, string b, int minLength)
    {
        int effectiveMin = Math.Max(1, minLength);
        regions.Select(r => r.Sequence).Should().OnlyHaveUniqueItems("distinct substrings reported once (INV-04)");

        foreach (var r in regions)
        {
            r.IsEmpty.Should().BeFalse("a region is a non-empty contiguous substring (INV-04)");
            r.Length.Should().BeGreaterThanOrEqualTo(effectiveMin,
                "every region has length ≥ max(1, minLength) (INV-04)");

            r.PositionInFirst.Should().BeInRange(0, a.Length - r.Length, "valid start in sequence1");
            r.PositionInSecond.Should().BeInRange(0, b.Length - r.Length, "valid start in sequence2");
            a.Substring(r.PositionInFirst, r.Length).Should().Be(r.Sequence,
                "region genuinely occurs at PositionInFirst in sequence1 (INV-04)");
            b.Substring(r.PositionInSecond, r.Length).Should().Be(r.Sequence,
                "region begins at its reported start in sequence2 (INV-04, right-maximal-per-start)");
        }
    }

    #endregion

    #region GENOMIC-COMMON-001 — Common Region Detection (BE: single input, disjoint, identical)

    #region Positive sanity — hand-computed documented LCS

    // Documented worked example (§7.1): ACGTACGT vs TTACGTGG → TACGT, length 5, pos 3 / 1.
    [Test]
    public void FindLongestCommonRegion_DocumentedWorkedExample_ReturnsTacgt()
    {
        var a = new DnaSequence("ACGTACGT");
        var b = new DnaSequence("TTACGTGG");

        CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(a, b);

        r.Sequence.Should().Be("TACGT", "documented LCS (§7.1)");
        r.Length.Should().Be(5);
        r.PositionInFirst.Should().Be(3, "TACGT starts at index 3 of ACGTACGT (INV-01)");
        r.PositionInSecond.Should().Be(1, "TACGT starts at index 1 of TTACGTGG (INV-01)");
        AssertWellFormedRegion(r, "ACGTACGT", "TTACGTGG");
    }

    // A planted common block ("GATTACA") is reported when present in both flanks.
    [Test]
    public void FindLongestCommonRegion_PlantedCommonBlock_IsReported()
    {
        var a = new DnaSequence("CCCCGATTACATTTT");
        var b = new DnaSequence("AAGATTACAGGGGGG");

        CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(a, b);

        r.Sequence.Should().Be("GATTACA", "the planted 7-mer is the longest shared contiguous block");
        AssertWellFormedRegion(r, "CCCCGATTACATTTT", "AAGATTACAGGGGGG");
    }

    #endregion

    #region BE — Boundary: single input (empty / single-char pair → no crash)

    // Empty first sequence: only the empty string qualifies → CommonRegion.None, no crash (§6.1).
    [Test]
    public void FindLongestCommonRegion_EmptyFirst_ReturnsNone_NoCrash()
    {
        var a = new DnaSequence("");
        var b = new DnaSequence("ACGTACGT");

        Action act = () => GenomicAnalyzer.FindLongestCommonRegion(a, b);
        act.Should().NotThrow("empty input is a documented boundary, not an error (§3.3, §6.1)");

        CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(a, b);
        r.IsEmpty.Should().BeTrue("empty sequence shares nothing → CommonRegion.None (INV-05)");
        r.PositionInFirst.Should().Be(-1);
        r.PositionInSecond.Should().Be(-1);
    }

    // Both sequences empty → CommonRegion.None, no crash, no DivideByZero in the LCS streaming.
    [Test]
    public void FindLongestCommonRegion_BothEmpty_ReturnsNone_NoCrash()
    {
        var empty = new DnaSequence("");

        Action act = () => GenomicAnalyzer.FindLongestCommonRegion(empty, empty);
        act.Should().NotThrow("the empty/empty pair is a documented boundary (§6.1)");

        GenomicAnalyzer.FindLongestCommonRegion(empty, empty).IsEmpty.Should().BeTrue("nothing shared (INV-05)");
    }

    // Single shared base — minimal non-empty LCS. "A" vs "TTTAT" → A at seq1 pos 0, seq2 pos 3 (S4).
    [Test]
    public void FindLongestCommonRegion_SingleCharOverlap_ReturnsLengthOneRegion()
    {
        var a = new DnaSequence("A");
        var b = new DnaSequence("TTTAT");

        CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(a, b);

        r.Sequence.Should().Be("A", "only the single base A is shared");
        r.Length.Should().Be(1);
        r.PositionInFirst.Should().Be(0);
        r.PositionInSecond.Should().Be(3, "first A in TTTAT is at index 3 (INV-01)");
        AssertWellFormedRegion(r, "A", "TTTAT");
    }

    // Single-char pair sharing NOTHING ("A" vs "C") → CommonRegion.None, no crash.
    [Test]
    public void FindLongestCommonRegion_SingleCharDisjoint_ReturnsNone()
    {
        var r = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence("A"), new DnaSequence("C"));
        r.IsEmpty.Should().BeTrue("A and C share no character (INV-05)");
    }

    // Fuzz: one operand empty, the other random → never throws, always CommonRegion.None.
    [Test]
    [CancelAfter(30_000)]
    public void FindLongestCommonRegion_OneEmptyOperand_RandomOther_AlwaysNone_NeverThrows()
    {
        var rng = new Random(175_001);
        var empty = new DnaSequence("");
        for (int trial = 0; trial < 600; trial++)
        {
            string s = RandomDna(rng, rng.Next(0, 40));
            var other = new DnaSequence(s);

            CommonRegion left = GenomicAnalyzer.FindLongestCommonRegion(empty, other);
            CommonRegion right = GenomicAnalyzer.FindLongestCommonRegion(other, empty);

            left.IsEmpty.Should().BeTrue("empty-on-left shares nothing (INV-05)");
            right.IsEmpty.Should().BeTrue("empty-on-right shares nothing (INV-05)");
            AssertWellFormedRegion(left, "", s);
            AssertWellFormedRegion(right, s, "");
        }
    }

    // Fuzz: single-character first operand over random second → length ≤ 1, never throws.
    [Test]
    [CancelAfter(30_000)]
    public void FindLongestCommonRegion_SingleCharFirst_RandomSecond_NeverThrows_MatchesOracle()
    {
        var rng = new Random(175_002);
        for (int trial = 0; trial < 600; trial++)
        {
            string a = Alphabet[rng.Next(Alphabet.Length)].ToString();
            string b = RandomDna(rng, rng.Next(0, 30));

            CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence(a), new DnaSequence(b));

            r.Length.Should().BeLessThanOrEqualTo(1, "a single-base sequence can share at most one base (INV-02)");
            AssertWellFormedRegion(r, a, b);
        }
    }

    #endregion

    #region BE — Boundary: disjoint inputs (no common character → None, no false common region)

    // All-A vs all-C: disjoint character sets → CommonRegion.None, no false region.
    [Test]
    public void FindLongestCommonRegion_DisjointHomopolymers_ReturnsNone()
    {
        var r = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence("AAAAAAAA"), new DnaSequence("CCCCCCCC"));
        r.IsEmpty.Should().BeTrue("AAAA… and CCCC… share no contiguous substring (INV-05)");
        r.Length.Should().Be(0);
    }

    // FindCommonRegions on disjoint inputs → empty enumeration, no false common region.
    [Test]
    public void FindCommonRegions_DisjointHomopolymers_ReturnsEmpty()
    {
        var regions = GenomicAnalyzer
            .FindCommonRegions(new DnaSequence("GGGGGGGG"), new DnaSequence("TTTTTTTT"), minLength: 1)
            .ToList();
        regions.Should().BeEmpty("GGGG… and TTTT… share no substring (INV-04)");
    }

    // Fuzz: operands drawn from disjoint sub-alphabets ({A,C} vs {G,T}) → no shared character at
    // all → always CommonRegion.None / empty regions, never a fabricated common region.
    [Test]
    [CancelAfter(30_000)]
    public void FindLongestCommonRegion_DisjointAlphabets_NeverEmitsFalseCommonRegion()
    {
        var rng = new Random(175_003);
        char[] left = { 'A', 'C' };
        char[] right = { 'G', 'T' };
        for (int trial = 0; trial < 600; trial++)
        {
            var sa = new StringBuilder();
            int la = rng.Next(0, 30);
            for (int i = 0; i < la; i++) sa.Append(left[rng.Next(left.Length)]);
            var sb = new StringBuilder();
            int lb = rng.Next(0, 30);
            for (int i = 0; i < lb; i++) sb.Append(right[rng.Next(right.Length)]);
            string a = sa.ToString(), b = sb.ToString();

            CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence(a), new DnaSequence(b));
            r.IsEmpty.Should().BeTrue("disjoint alphabets share no character ⇒ no common region (INV-05)");

            var regions = GenomicAnalyzer.FindCommonRegions(new DnaSequence(a), new DnaSequence(b), 1).ToList();
            regions.Should().BeEmpty("no shared substring ⇒ empty enumeration (INV-04)");
            AssertWellFormedRegion(r, a, b);
        }
    }

    #endregion

    #region BE — Boundary: identical inputs (whole sequence at 0/0, deterministic)

    // Identical sequences: the whole sequence is the LCS at positions 0/0 (a string is a
    // substring of itself, §6.1).
    [Test]
    public void FindLongestCommonRegion_Identical_ReturnsWholeSequenceAtZeroZero()
    {
        const string seq = "ACGTACGTTT";
        var a = new DnaSequence(seq);
        var b = new DnaSequence(seq);

        CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(a, b);

        r.Sequence.Should().Be(seq, "a string is a substring of itself → whole sequence (§6.1)");
        r.Length.Should().Be(seq.Length);
        r.PositionInFirst.Should().Be(0);
        r.PositionInSecond.Should().Be(0);
        AssertWellFormedRegion(r, seq, seq);
    }

    // Determinism: identical inputs produce the SAME region across repeated calls.
    [Test]
    public void FindLongestCommonRegion_Identical_DeterministicAcrossCalls()
    {
        var a = new DnaSequence("ACGTACGTACGT");
        var b = new DnaSequence("ACGTACGTACGT");

        CommonRegion first = GenomicAnalyzer.FindLongestCommonRegion(a, b);
        CommonRegion second = GenomicAnalyzer.FindLongestCommonRegion(a, b);

        second.Sequence.Should().Be(first.Sequence, "result is deterministic (INV-03)");
        second.PositionInFirst.Should().Be(first.PositionInFirst);
        second.PositionInSecond.Should().Be(first.PositionInSecond);
    }

    // FindCommonRegions on identical inputs: for each start i in seq2 the longest match is the
    // whole suffix seq[i..]; the deduplicated set is exactly {seq, seq[1..], …} of length ≥ min.
    [Test]
    public void FindCommonRegions_Identical_ReturnsRightMaximalSuffixesAtEachStart()
    {
        const string seq = "ACGT";
        var a = new DnaSequence(seq);
        var b = new DnaSequence(seq);

        var regions = GenomicAnalyzer.FindCommonRegions(a, b, minLength: 1).ToList();

        regions.Select(r => r.Sequence).Should().BeEquivalentTo(
            new[] { "ACGT", "CGT", "GT", "T" },
            "for identical inputs the longest match at each start i is the suffix seq[i..] (INV-04)");
        AssertWellFormedRegions(regions, seq, seq, 1);
    }

    // Fuzz: identical random pairs → LCS is the whole sequence at 0/0 (or None when empty),
    // deterministic, never throws.
    [Test]
    [CancelAfter(30_000)]
    public void FindLongestCommonRegion_IdenticalRandom_WholeSequence_Deterministic()
    {
        var rng = new Random(175_004);
        for (int trial = 0; trial < 600; trial++)
        {
            string s = RandomDna(rng, rng.Next(0, 40));
            var a = new DnaSequence(s);
            var b = new DnaSequence(s);

            CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(a, b);

            if (s.Length == 0)
            {
                r.IsEmpty.Should().BeTrue("empty identical pair shares nothing (INV-05)");
            }
            else
            {
                r.Length.Should().Be(s.Length, "identical sequences share the whole sequence (INV-02)");
                r.Sequence.Should().Be(s);
                r.PositionInFirst.Should().Be(0);
                r.PositionInSecond.Should().Be(0);
            }

            // Determinism cross-check.
            GenomicAnalyzer.FindLongestCommonRegion(a, b).Sequence.Should().Be(r.Sequence, "deterministic (INV-03)");
            AssertWellFormedRegion(r, s, s);
        }
    }

    #endregion

    #region BE — Boundary: minLength guards (0, −1, int.MinValue ≡ 1; int.MaxValue → empty)

    // minLength below 1 is treated as 1 — no empty / length-0 region, same set as minLength=1.
    [Test]
    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void FindCommonRegions_MinLengthBelowOne_TreatedAsOne_NoEmptyRegion(int minLength)
    {
        var a = new DnaSequence("ACGT");
        var b = new DnaSequence("ACGT");

        var regions = GenomicAnalyzer.FindCommonRegions(a, b, minLength).ToList();
        var baseline = GenomicAnalyzer.FindCommonRegions(a, b, 1).ToList();

        regions.Should().OnlyContain(r => r.Length >= 1, "minLength < 1 ≡ 1, no length-0 region (INV-04)");
        regions.Select(r => r.Sequence).Should().BeEquivalentTo(
            baseline.Select(r => r.Sequence), "minLength < 1 behaves exactly like minLength == 1 (§6.1)");
    }

    // int.MaxValue minLength: no region can be that long → empty enumeration, no overflow/crash.
    [Test]
    [CancelAfter(15_000)]
    public void FindCommonRegions_MinLengthMaxInt_ReturnsEmpty_NoOverflow()
    {
        var a = new DnaSequence("ACGTACGT");
        var b = new DnaSequence("ACGTACGT");

        Action act = () => GenomicAnalyzer.FindCommonRegions(a, b, int.MaxValue).ToList();
        act.Should().NotThrow("an unreachable minLength is a boundary, not an error (BE: MaxInt)");

        GenomicAnalyzer.FindCommonRegions(a, b, int.MaxValue)
            .Should().BeEmpty("no common substring is int.MaxValue long ⇒ empty");
    }

    #endregion

    #region BE — Broad fuzz: random pairs / minLength never crash, match the documented contract

    // FindLongestCommonRegion over random pairs: never throws, the reported length equals the
    // independent DP-oracle LCS length, and the region is well-formed in both operands.
    [Test]
    [CancelAfter(60_000)]
    public void FindLongestCommonRegion_RandomPairs_NeverThrows_MatchesOracle()
    {
        var rng = new Random(175_005);
        for (int trial = 0; trial < 1500; trial++)
        {
            string a = RandomDna(rng, rng.Next(0, 35));
            string b = RandomDna(rng, rng.Next(0, 35));

            CommonRegion r = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence(a), new DnaSequence(b));

            AssertWellFormedRegion(r, a, b);
        }
    }

    // FindCommonRegions over random pairs and random minLength: never throws/hangs, every reported
    // region is a genuine common contiguous substring of length ≥ max(1, minLength).
    [Test]
    [CancelAfter(60_000)]
    public void FindCommonRegions_RandomPairsAndMinLength_NeverThrows_WellFormed()
    {
        var rng = new Random(175_006);
        for (int trial = 0; trial < 1500; trial++)
        {
            string a = RandomDna(rng, rng.Next(0, 35));
            string b = RandomDna(rng, rng.Next(0, 35));
            int minLength = rng.Next(-3, 12);

            var regions = GenomicAnalyzer.FindCommonRegions(new DnaSequence(a), new DnaSequence(b), minLength).ToList();

            AssertWellFormedRegions(regions, a, b, minLength);

            // Cross-check against the LCS oracle: the single longest region (if any) cannot exceed
            // the true LCS length, and when a common substring ≥ max(1,minLength) exists, at least
            // one region must be reported.
            var (oracleLen, _) = OracleLcs(a, b);
            int effMin = Math.Max(1, minLength);
            if (regions.Count > 0)
                regions.Max(r => r.Length).Should().BeLessThanOrEqualTo(oracleLen,
                    "no reported region exceeds the true longest common substring (INV-02/INV-04)");
            if (oracleLen >= effMin)
                regions.Should().NotBeEmpty(
                    "a common substring of length ≥ max(1,minLength) exists ⇒ at least one region (INV-04)");
        }
    }

    // Determinism over random pairs: the unit's longest-region result is stable across calls.
    [Test]
    [CancelAfter(30_000)]
    public void FindLongestCommonRegion_RandomPairs_Deterministic()
    {
        var rng = new Random(175_007);
        for (int trial = 0; trial < 400; trial++)
        {
            var a = new DnaSequence(RandomDna(rng, rng.Next(1, 30)));
            var b = new DnaSequence(RandomDna(rng, rng.Next(1, 30)));

            CommonRegion first = GenomicAnalyzer.FindLongestCommonRegion(a, b);
            CommonRegion second = GenomicAnalyzer.FindLongestCommonRegion(a, b);

            second.Sequence.Should().Be(first.Sequence, "result set is deterministic (INV-03)");
            second.PositionInFirst.Should().Be(first.PositionInFirst);
            second.PositionInSecond.Should().Be(first.PositionInSecond);
        }
    }

    #endregion

    #endregion
}
