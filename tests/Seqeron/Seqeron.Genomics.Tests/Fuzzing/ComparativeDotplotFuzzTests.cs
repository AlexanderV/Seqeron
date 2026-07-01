using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Comparative-Genomics dot plot — the word-match (k-tuple)
/// dot matrix of two sequences (COMPGEN-DOTPLOT-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang / quadratic
/// blow-up on highly repetitive input, no IndexOutOfRange / negative-length
/// Substring when a sequence is shorter than the word size, no DivideByZero on a
/// zero-size matrix, and no nonsense output (a dot whose coordinate is out of
/// bounds, or whose word does not actually match). Every input must resolve to
/// EITHER a well-defined, theory-correct dot set, OR a *documented, intentional*
/// validation exception (ArgumentOutOfRangeException for a non-positive
/// word/step size). A raw runtime exception, a hang, or an out-of-bounds dot is a
/// bug, not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-DOTPLOT-001 — Word-match (k-tuple) dot plot
/// Checklist: docs/checklists/03_FUZZING.md, row 134.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation (граничні значення: empty, 1-base, k-boundary) —
///          the degenerate dot-plot boundaries called out in the checklist row:
///          EMPTY (one/both sequences empty), SINGLE BASE (length-1 input, or a
///          sequence shorter than the word size, so no length-k word is
///          extractable), PALINDROME (a string-palindrome vs itself → the
///          documented main + anti-diagonal symmetric match pattern), and
///          REPEAT-RICH (a tandem-repeat sequence vs itself → many off-diagonal
///          dots at the repeat period, which must terminate, not hang).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The dot-plot contract under test (Dot_Plot_Generation.md)
/// ───────────────────────────────────────────────────────────────────────────
/// GenerateDotPlot reports the EMBOSS dottup word-match relation
///   D = { (i, j) : A[i..i+w−1] = B[j..j+w−1] }, 0 ≤ i ≤ n−w, 0 ≤ j ≤ m−w,
/// where A = sequence1 (→ x), B = sequence2 (→ y), w = wordSize, equality is
/// character-by-character and CASE-INSENSITIVE (both upper-cased), and i advances
/// in steps of stepSize along A (Gibbs &amp; McIntyre 1970; Rice et al. 2000 dottup).
/// The API entry under test is
///   ComparativeGenomics.GenerateDotPlot(
///       string sequence1, string sequence2,
///       int wordSize = 10, int stepSize = 1)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs
///    lines 1169–1207; the empty/short guards at lines 1189–1192).
///
/// THE DOCUMENTED INVARIANTS (Dot_Plot_Generation.md §2.4):
///   • INV-01: (x, y) ∈ output ⇔ A[x..x+w−1] = B[y..y+w−1] (case-insensitive).
///   • INV-02: self-comparison (A = B) contains every (i, i), 0 ≤ i ≤ n−w (the
///             FULL main diagonal — a word always matches itself).
///   • INV-03: empty output when either sequence is null/empty OR shorter than w.
///   • INV-04: all overlapping occurrences are reported; x is a multiple of stepSize.
/// The documented edge cases (§3.3 / §6.1):
///   • null / empty / shorter-than-w sequence → empty result (no exception);
///   • disjoint alphabets → empty result;
///   • wordSize ≤ 0 or stepSize ≤ 0 → ArgumentOutOfRangeException;
///   • highly repetitive input → output size is O(n·m) worst case but the
///     enumeration must terminate (§4.3, §6.2).
///
/// The four BE checklist targets map to these documented behaviours:
///   • empty       → empty dot plot, NO DivideByZero / IndexOutOfRange on a
///                   zero-size matrix (INV-03; lines 1189–1190).
///   • single base → a length-1 (or any length &lt; w) sequence yields NO dot for
///                   w &gt; 1 (no length-w word can start, INV-03; line 1191), and
///                   the documented single-residue rule for w = 1; NEVER a
///                   negative-length Substring or IndexOutOfRange.
///   • palindrome  → a STRING-palindrome P (P reads the same forwards and
///                   backwards) compared to ITSELF at w = 1 yields BOTH the main
///                   diagonal (i,i) (INV-02) AND the anti-diagonal (i, n−1−i),
///                   because P[i] == P[n−1−i]; the dot set is symmetric under
///                   (x,y) ↔ (y,x) for any self-comparison.
///   • repeat-rich → a tandem repeat (period p) vs itself emits dots on every
///                   diagonal offset that is a multiple of p; the enumeration
///                   must TERMINATE (no quadratic hang) under [CancelAfter].
/// A positive-sanity test pins the documented golden vector (§7.1) and the full
/// main diagonal of a self-comparison, plus a known shared subword at the correct
/// (x, y) offsets.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeDotplotFuzzTests
{
    #region Test Data and Helpers

    /// <summary>Collects the dot plot into a set of (x, y) coordinate pairs.</summary>
    private static HashSet<(int x, int y)> Dots(string s1, string s2, int wordSize, int stepSize = 1)
        => ComparativeGenomics.GenerateDotPlot(s1, s2, wordSize, stepSize).ToHashSet();

    /// <summary>
    /// A well-formed dot plot satisfies the documented structural contract for ANY input
    /// (Dot_Plot_Generation.md §2.4, §3.2):
    ///   • every x is a valid 0-based word START in sequence1: 0 ≤ x ≤ |s1| − w,
    ///   • every y is a valid 0-based word START in sequence2: 0 ≤ y ≤ |s2| − w,
    ///   • every x is a multiple of stepSize (INV-04),
    ///   • the word actually matches, case-insensitively, at (x, y) (INV-01 — the dot is real),
    ///   • coordinates are unique within a single enumeration (no duplicate (x,y)).
    /// This is the always-on safety net: no dot may point outside the matrix or to a non-match.
    /// </summary>
    private static void AssertWellFormed(
        IReadOnlyList<(int x, int y)> dots, string s1, string s2, int wordSize, int stepSize = 1)
    {
        var seen = new HashSet<(int, int)>();
        foreach (var (x, y) in dots)
        {
            x.Should().BeInRange(0, s1.Length - wordSize,
                "x is a valid word start in sequence1: 0 ≤ x ≤ |s1|−w (§3.2)");
            y.Should().BeInRange(0, s2.Length - wordSize,
                "y is a valid word start in sequence2: 0 ≤ y ≤ |s2|−w (§3.2)");
            (x % stepSize).Should().Be(0, "x must be a multiple of stepSize (INV-04)");

            // INV-01: the dot is real — the length-w words actually match case-insensitively.
            string w1 = s1.Substring(x, wordSize).ToUpperInvariant();
            string w2 = s2.Substring(y, wordSize).ToUpperInvariant();
            w1.Should().Be(w2, "a dot at (x,y) means A[x..x+w−1] == B[y..y+w−1] (INV-01)");

            seen.Add((x, y)).Should().BeTrue("each (x,y) coordinate appears at most once per enumeration");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  COMPGEN-DOTPLOT-001 — Word-match dot plot : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region COMPGEN-DOTPLOT-001 — Word-match (k-tuple) dot plot

    #region BE — Boundary: empty (zero-size matrix)

    /// <summary>
    /// BE: both sequences empty is the floor of the size axis — a 0×0 matrix. The documented result
    /// is the empty dot plot (INV-03), with no DivideByZero / IndexOutOfRange on the zero-size matrix
    /// (the IsNullOrEmpty guard, ComparativeGenomics.cs line 1189, short-circuits before any indexing).
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_BothEmpty_EmptyNoCrash()
    {
        List<(int x, int y)> dots = null!;
        FluentActions.Invoking(() => dots = Dots("", "", wordSize: 1).ToList())
            .Should().NotThrow("a 0×0 matrix must yield a clean empty dot plot, never a crash");

        dots.Should().BeEmpty("no length-w word can be formed from empty input (INV-03)");
        AssertWellFormed(dots, "", "", wordSize: 1);
    }

    /// <summary>
    /// BE: ONE empty sequence (the other non-empty) is the asymmetric floor — a 0×m / n×0 matrix.
    /// No word can be located in the empty side, so the result is empty (INV-03). Both argument
    /// orders are probed: an empty sequence1 hits the slide-loop's empty guard, an empty sequence2
    /// hits the suffix-tree-over-empty path; neither may throw.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_OneEmpty_EmptyNoCrash()
    {
        List<(int x, int y)> fwd = null!, rev = null!;
        FluentActions.Invoking(() => fwd = Dots("ACGTACGT", "", wordSize: 2).ToList())
            .Should().NotThrow("non-empty vs empty must not crash on the empty side");
        FluentActions.Invoking(() => rev = Dots("", "ACGTACGT", wordSize: 2).ToList())
            .Should().NotThrow("empty vs non-empty must be equally safe");

        fwd.Should().BeEmpty("no word can match into an empty sequence2 (INV-03)");
        rev.Should().BeEmpty("no word can start in an empty sequence1 (INV-03)");
    }

    #endregion

    #region BE — Boundary: single base / shorter than word size

    /// <summary>
    /// BE: a length-1 sequence with the default wordSize 10 is the canonical "shorter than k" corner
    /// — the case most likely to throw a negative-length Substring (Substring(0, 10) on a 1-char
    /// string) or IndexOutOfRange. The documented behaviour is an empty result (INV-03; the
    /// `Length &lt; wordSize` guard, line 1191), NEVER an exception.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_SingleBaseDefaultWord_EmptyNoIndexOutOfRange()
    {
        List<(int x, int y)> dots = null!;
        FluentActions.Invoking(() => dots = ComparativeGenomics.GenerateDotPlot("A", "A").ToList())
            .Should().NotThrow("a 1-base sequence shorter than the default word size must not index past the end");

        dots.Should().BeEmpty("no length-10 word can start in a 1-base sequence (INV-03)");
    }

    /// <summary>
    /// BE: both sequences exactly one base shorter than the word size — the off-by-one boundary of the
    /// short-circuit guard. |s| = w − 1 means n − w = −1, so the slide loop must not execute even once
    /// (a negative upper bound). Documented empty result, no negative-length Substring.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_OneShorterThanWord_EmptyNoNegativeSubstring()
    {
        // |"ACGT"| = 4, wordSize 5 → n − w = −1 on both axes.
        List<(int x, int y)> dots = null!;
        FluentActions.Invoking(() => dots = Dots("ACGT", "ACGT", wordSize: 5).ToList())
            .Should().NotThrow("|s| = w−1 must not produce a negative-length Substring");

        dots.Should().BeEmpty("no length-5 word can start when both sequences are length 4 (INV-03)");
    }

    /// <summary>
    /// BE: a single base vs itself at wordSize 1 is the smallest NON-degenerate dot plot — a 1×1
    /// matrix with exactly one cell. The documented single-residue rule places a dot iff the
    /// characters match (§2.2, w=1), so "A" vs "A" → exactly {(0,0)} and "A" vs "C" → empty. Pins
    /// that the boundary of the minimal valid input is the correct single dot, not empty/crash.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_SingleBaseVsItself_OneDotWhenEqual()
    {
        var match = Dots("A", "A", wordSize: 1);
        var mismatch = Dots("A", "C", wordSize: 1);

        AssertWellFormed(match.ToList(), "A", "A", wordSize: 1);
        match.Should().BeEquivalentTo(new[] { (0, 0) },
            "a 1×1 matrix of equal residues has exactly one dot at (0,0) (single-residue rule §2.2)");
        mismatch.Should().BeEmpty("a dot is placed only on an exact match (§2.2, disjoint → empty §6.1)");
    }

    #endregion

    #region BE — Boundary: palindrome (main + anti-diagonal symmetry)

    /// <summary>
    /// BE: a STRING-palindrome P (P reads the same forwards and backwards) compared to ITSELF at
    /// wordSize 1 is the canonical symmetric dot plot. Because P[i] == P[n−1−i], the dot set contains
    /// BOTH the main diagonal (i, i) (INV-02 — a residue always matches itself) AND the anti-diagonal
    /// (i, n−1−i) (Wikipedia: a palindrome shows as an anti-diagonal in a self dot plot). We pin the
    /// EXACT match set for a clean palindrome with distinct mirror-paired residues, so a missing
    /// anti-diagonal or a wrong coordinate convention fails.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_StringPalindromeSelf_MainAndAntiDiagonal()
    {
        // ABCBA: P[0]=P[4]='A', P[1]=P[3]='B', P[2]='C'. Distinct residues except the mirror pairs,
        // so the ONLY equal-character pairs are the main diagonal and the anti-diagonal.
        const string p = "ABCBA";
        int n = p.Length;
        var dots = Dots(p, p, wordSize: 1);

        AssertWellFormed(dots.ToList(), p, p, wordSize: 1);

        var expected = new HashSet<(int, int)>();
        for (int i = 0; i < n; i++)
        {
            expected.Add((i, i));            // main diagonal: P[i] == P[i]
            expected.Add((i, n - 1 - i));    // anti-diagonal: P[i] == P[n−1−i] (palindrome)
        }
        dots.Should().BeEquivalentTo(expected,
            "a string-palindrome self dot plot is exactly the main diagonal ∪ anti-diagonal (INV-02 + Wikipedia palindrome pattern)");
    }

    /// <summary>
    /// BE: symmetry property of ANY self-comparison — the dot plot of A vs A is invariant under
    /// transposition (x,y) ↔ (y,x), because A[x..]==A[y..] ⇔ A[y..]==A[x..] (INV-01 is a symmetric
    /// relation when both axes are the same sequence). This is the structural symmetry the palindrome
    /// case is a special case of; we assert it directly on a palindrome at a word size &gt; 1.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_PalindromeSelf_DotSetIsTransposeSymmetric()
    {
        // A reverse-complement-free string palindrome over a 3-letter alphabet, length 7.
        const string p = "GATCTAG"; // reversed == "GATCTAG"? check: reverse = GATCTAG → palindrome.
        var dots = Dots(p, p, wordSize: 2);

        AssertWellFormed(dots.ToList(), p, p, wordSize: 2);

        var transposed = dots.Select(d => (d.y, d.x)).ToHashSet();
        dots.Should().BeEquivalentTo(transposed,
            "a self-comparison dot set is symmetric under (x,y) ↔ (y,x) (INV-01 over a single sequence)");
        // The main diagonal must be present at w=2: every (i,i) for 0 ≤ i ≤ n−2.
        for (int i = 0; i <= p.Length - 2; i++)
            dots.Should().Contain((i, i), "the full main diagonal is present in any self-comparison (INV-02)");
    }

    #endregion

    #region BE — Boundary: repeat-rich (off-diagonal flood, must terminate)

    /// <summary>
    /// BE: a tandem repeat (period 2) vs itself is the worst case for output size — the documented
    /// O(n·m) blow-up (§4.3, §6.2). "ATATAT…" has the word "AT" starting at every even index and "TA"
    /// at every odd index, so the dots form off-diagonal lines at every even diagonal offset. The
    /// enumeration must TERMINATE (no quadratic hang) under a tight cancel budget, and every dot must
    /// be well-formed. We also pin the EXACT documented count to prove no member is dropped/duplicated.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_TandemRepeatSelf_TerminatesWithExactCount()
    {
        // (AT)×10 = 20 chars; wordSize 2 → word starts 0..18 (19 of them). A start at index i is
        // "AT" when i is even, "TA" when i is odd. In s the word "AT" occurs at even j (0,2,…,18 →
        // 10 positions) and "TA" at odd j (1,3,…,17 → 9 positions). Even starts: i ∈ {0,2,…,18} → 10
        // starts, each matching the 10 "AT" occurrences. Odd starts: i ∈ {1,3,…,17} → 9 starts, each
        // matching the 9 "TA" occurrences. Total = 10×10 + 9×9 = 181.
        string s = string.Concat(Enumerable.Repeat("AT", 10)); // length 20
        const int w = 2;

        List<(int x, int y)> dots = null!;
        FluentActions.Invoking(() => dots = Dots(s, s, w).ToList())
            .Should().NotThrow("a repeat-rich self dot plot must terminate, not hang");

        AssertWellFormed(dots, s, s, w);

        dots.Should().HaveCount(181,
            "(AT)×10 self dot plot at w=2: 10 even starts ×10 \"AT\" + 9 odd starts ×9 \"TA\" = 181 dots (§4.3 O(n·m))");

        // The repeat structure: dots lie on diagonals whose offset (y − x) is EVEN (same period parity).
        dots.Should().OnlyContain(d => (d.y - d.x) % 2 == 0,
            "a period-2 tandem repeat matches only at even diagonal offsets (the repeat period)");
    }

    /// <summary>
    /// BE: a homopolymer (period 1) vs itself is the absolute worst case — every word equals every
    /// other word, so the dot plot is the COMPLETE matrix of (n−w+1)² dots (§4.3 "self-comparison of a
    /// homopolymer"). The hazard is a quadratic hang; sizes are kept modest with a cancel budget. We
    /// pin the exact full-matrix count and that no dot escapes the matrix.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_HomopolymerSelf_FullMatrixTerminates()
    {
        string s = new string('A', 30); // 30 A's
        const int w = 5;
        int starts = s.Length - w + 1;   // 26 word starts per axis

        List<(int x, int y)> dots = null!;
        FluentActions.Invoking(() => dots = Dots(s, s, w).ToList())
            .Should().NotThrow("a homopolymer self dot plot is the dense worst case but must terminate");

        AssertWellFormed(dots, s, s, w);
        dots.Should().HaveCount(starts * starts,
            "every length-w word of a homopolymer equals every other → the complete (n−w+1)² matrix (§4.3)");
    }

    #endregion

    #region Positive sanity — the documented contract on hand-built examples

    /// <summary>
    /// Positive sanity (Dot_Plot_Generation.md §7.1 golden vector): AGCGT (x) vs AT (y) at wordSize 1.
    /// The only equal characters are the A's at (0,0) and the T's at (4,1), so the dot set is exactly
    /// {(0,0),(4,1)}. This is the load-bearing correctness check distinguishing a theory-correct
    /// word-match enumeration from a wrong coordinate convention or a first-occurrence-only impl.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_GoldenVector_HuttleyExample()
    {
        var dots = Dots("AGCGT", "AT", wordSize: 1);

        AssertWellFormed(dots.ToList(), "AGCGT", "AT", wordSize: 1);
        dots.Should().BeEquivalentTo(new[] { (0, 0), (4, 1) },
            "AGCGT vs AT at w=1: A==A at (0,0) and T==T at (4,1) — the documented golden vector (§7.1)");
    }

    /// <summary>
    /// Positive sanity — full main diagonal (INV-02). A self-comparison of a sequence with all-distinct
    /// residues at wordSize 1 is EXACTLY the main diagonal {(i,i)} and nothing else: every word matches
    /// itself, and distinct residues forbid any off-diagonal dot. Pins that the diagonal is complete and
    /// no spurious off-diagonal dot is invented.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_SelfDistinctResidues_ExactMainDiagonal()
    {
        const string s = "ACGT";
        var dots = Dots(s, s, wordSize: 1);

        AssertWellFormed(dots.ToList(), s, s, wordSize: 1);
        var expected = Enumerable.Range(0, s.Length).Select(i => (i, i)).ToHashSet();
        dots.Should().BeEquivalentTo(expected,
            "ACGT self-comparison is exactly the full main diagonal (i,i) (INV-02), no off-diagonal noise");
    }

    /// <summary>
    /// Positive sanity — a KNOWN shared subword appears at the correct (x, y) offsets. The word
    /// "CATG" occurs once in each input, at x=3 in sequence1 and y=2 in sequence2, so the only dot is
    /// (3, 2). Pins the exact coordinate convention (x = start in seq1, y = start in seq2) — a
    /// transposed implementation would report (2, 3).
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_KnownSharedSubword_CorrectOffset()
    {
        //          0123456
        const string s1 = "GGGCATGTTT"; // "CATG" at index 3
        const string s2 = "TTCATGAA";   // "CATG" at index 2
        var dots = Dots(s1, s2, wordSize: 4);

        AssertWellFormed(dots.ToList(), s1, s2, wordSize: 4);
        dots.Should().BeEquivalentTo(new[] { (3, 2) },
            "the shared word CATG starts at x=3 in seq1 and y=2 in seq2 → single dot (3,2); x↔seq1, y↔seq2");
    }

    /// <summary>
    /// Positive sanity — case-insensitivity (§3.3). Lower-case input must match its upper-case
    /// counterpart because both axes are upper-cased before matching. "acgt" vs "ACGT" at wordSize 1
    /// yields the same full main diagonal as the all-upper case, with no extra or missing dots.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void GenerateDotPlot_MixedCase_IsCaseInsensitive()
    {
        var dots = Dots("acgt", "ACGT", wordSize: 1);

        AssertWellFormed(dots.ToList(), "acgt", "ACGT", wordSize: 1);
        var expected = Enumerable.Range(0, 4).Select(i => (i, i)).ToHashSet();
        dots.Should().BeEquivalentTo(expected,
            "matching upper-cases both sequences → case-insensitive equality (§3.3); diagonal preserved");
    }

    /// <summary>
    /// BE/Validation: a non-positive word or step size is undefined for a word-match dot plot and is
    /// rejected EAGERLY with ArgumentOutOfRangeException (§3.3, §6.1; ComparativeGenomics.cs lines
    /// 1175–1178) — before enumeration begins, not deferred to first MoveNext. Pins the documented
    /// validation boundary on both parameters and both non-positive directions (0 and negative).
    /// </summary>
    [Test]
    public void GenerateDotPlot_NonPositiveWordOrStep_ThrowsEagerly()
    {
        // Eager throw: calling the method (without enumerating) must already throw.
        FluentActions.Invoking(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", wordSize: 0))
            .Should().Throw<ArgumentOutOfRangeException>("wordSize 0 is undefined (§3.3)")
            .Which.ParamName.Should().Be("wordSize");
        FluentActions.Invoking(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", wordSize: -1))
            .Should().Throw<ArgumentOutOfRangeException>("negative wordSize is undefined (§3.3)")
            .Which.ParamName.Should().Be("wordSize");
        FluentActions.Invoking(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", wordSize: 2, stepSize: 0))
            .Should().Throw<ArgumentOutOfRangeException>("stepSize 0 is undefined (§3.3)")
            .Which.ParamName.Should().Be("stepSize");
        FluentActions.Invoking(() => ComparativeGenomics.GenerateDotPlot("ACGT", "ACGT", wordSize: 2, stepSize: -1))
            .Should().Throw<ArgumentOutOfRangeException>("negative stepSize is undefined (§3.3)")
            .Which.ParamName.Should().Be("stepSize");
    }

    /// <summary>
    /// Fuzz sweep — across many random sequence pairs of varied lengths (including empty and
    /// shorter-than-w), alphabets, repeat content, word sizes and step sizes, GenerateDotPlot must
    /// ALWAYS terminate, never throw, and return a well-formed dot set: every coordinate in-bounds,
    /// every x a multiple of stepSize (INV-04), and every dot a real case-insensitive word match
    /// (INV-01). As a partition floor we also pin that a self-comparison ALWAYS contains the full main
    /// diagonal (INV-02). The exact dot SET is asserted only in the dedicated golden tests over
    /// hand-built sequences; here the well-formed net + the self-diagonal floor are the theory-correct,
    /// input-independent checks.
    /// </summary>
    [Test]
    [CancelAfter(120000)]
    public void GenerateDotPlot_RandomPairs_AlwaysWellFormedWithSelfDiagonal()
    {
        const string alphabet = "ACGT";
        for (int seed = 0; seed < 80; seed++)
        {
            var rng = new Random(seed);
            int len1 = rng.Next(0, 25);
            int len2 = rng.Next(0, 25);
            int wordSize = rng.Next(1, 6);
            int stepSize = rng.Next(1, 4);

            // Occasionally bias toward repeat-rich content to stress the off-diagonal flood.
            string s1 = RandomSeq(rng, len1, alphabet, repeatBias: seed % 3 == 0);
            string s2 = RandomSeq(rng, len2, alphabet, repeatBias: seed % 3 == 1);

            List<(int x, int y)> dots = null!;
            FluentActions.Invoking(() => dots = Dots(s1, s2, wordSize, stepSize).ToList())
                .Should().NotThrow($"random pair must never crash (seed {seed})");
            AssertWellFormed(dots, s1, s2, wordSize, stepSize);

            // INV-02 floor: a self-comparison at stepSize 1 contains every (i,i) for valid i.
            if (s1.Length >= wordSize)
            {
                var self = Dots(s1, s1, wordSize, stepSize: 1).ToList();
                AssertWellFormed(self, s1, s1, wordSize, stepSize: 1);
                for (int i = 0; i <= s1.Length - wordSize; i++)
                    self.Should().Contain((i, i),
                        $"self-comparison always contains the full main diagonal (INV-02, seed {seed})");
            }
        }
    }

    /// <summary>Builds a random sequence; when <paramref name="repeatBias"/> is set, tiles a short motif to force repeats.</summary>
    private static string RandomSeq(Random rng, int len, string alphabet, bool repeatBias)
    {
        if (len == 0) return string.Empty;
        if (repeatBias)
        {
            int period = rng.Next(1, 4);
            var motif = new char[period];
            for (int i = 0; i < period; i++) motif[i] = alphabet[rng.Next(alphabet.Length)];
            var sb = new System.Text.StringBuilder(len);
            for (int i = 0; i < len; i++) sb.Append(motif[i % period]);
            return sb.ToString();
        }
        var chars = new char[len];
        for (int i = 0; i < len; i++) chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    #endregion

    #endregion
}
