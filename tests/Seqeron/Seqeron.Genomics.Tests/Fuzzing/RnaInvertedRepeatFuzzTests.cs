// RNA-INVERT-001 — RNA Inverted-Repeat (potential stem) detection.
// Fuzz tests (strategy BE = Boundary Exploitation).
// Algorithm doc: docs/algorithms/RNA_Secondary_Structure/Inverted_Repeats.md
// Canonical tests: tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_FindInvertedRepeats_Tests.cs (RNA-INVERT-001)
// Evidence: docs/Evidence/RNA-INVERT-001-Evidence.md
// Source: RnaSecondaryStructure.FindInvertedRepeats — RnaSecondaryStructure.cs.
//         Alamro et al. (2021) IUPACpal, BMC Bioinformatics 22:51. doi:10.1186/s12859-021-03983-2
//         (perfect, k=0 antiparallel reverse-complement W G W̄ᴿ pattern; no G-U wobble — ASM-01).

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for RNA-INVERT-001 — <see cref="RnaSecondaryStructure.FindInvertedRepeats(string,int,int,int)"/>,
/// detection of perfect inverted repeats in an RNA/DNA sequence: a left arm <c>W</c>, an
/// intervening loop <c>G</c>, and a right arm equal to the <b>reverse complement</b> of <c>W</c>
/// (the <c>W G W̄ᴿ</c> pattern, |G| ≥ 0) — the sequence signature of a potential hairpin stem.
/// Lives in src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs.
///
/// This is the RNA-context inverted-repeat unit (checklist row 151). It is DISTINCT from the
/// DNA repeat-analysis unit REP-INV-001 (row 15, <c>RepeatFinder.FindInvertedRepeats</c>): this
/// detector uses the RNA complement table (<c>GetRnaComplementBase</c>: A⟷U, G⟷C, T→A, IUPAC
/// degenerate codes mapped, non-IUPAC passed through) and reports MAXIMAL, NON-OVERLAPPING perfect
/// stems via outward antiparallel extension over a bounded loop window
/// (Inverted_Repeats.md §1, §4).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs and asserts the code NEVER fails
/// in an undisciplined way: no hang, no unhandled runtime exception (IndexOutOfRange / negative
/// substring length off a single/short sequence, ArgumentOutOfRange from internal indexing), and
/// no nonsense output (a "repeat" whose arms are not actually reverse complements, out of bounds,
/// or shorter than minLength). Every input must resolve EITHER to well-formed, theory-correct
/// inverted repeats OR to a documented, intentional outcome (an empty result).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Strategy BE = Boundary Exploitation — targets "palindrome, no complementarity, single base"
/// ───────────────────────────────────────────────────────────────────────────
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE), row 151.
///
///   • palindrome — a perfectly self-complementary stem (the canonical POSITIVE). The doc's
///       worked example <c>UUACGAAAAAACGUAA</c> (Inverted_Repeats.md §7.1) MUST be detected with
///       the EXACT documented coordinates (Start1 0, End1 4, Start2 11, End2 15, Length 5). The
///       limiting case |G| = 0 (a true palindrome with no spacer, e.g. GGGGCCCC) is also checked
///       at the minSpacing = 0 boundary — arm spans exact, complementarity exact, no off-by-one.
///
///   • no complementarity — a sequence with NO self-complementary region (e.g. an all-A
///       homopolymer: A's reverse complement is U, which never occurs) yields an EMPTY result,
///       no crash, no false repeat (Inverted_Repeats.md §6.1 "no complement possibilities").
///
///   • single base — a length-1 (and length-0/empty/null) sequence cannot hold two arms plus a
///       loop, so it yields an EMPTY result. The guard <c>n &lt; 2·minLength + minSpacing</c> must
///       fire BEFORE any indexing — never an IndexOutOfRange or a negative substring length on the
///       degenerate single/short input (Inverted_Repeats.md §3.3, §6.1).
///
/// Watched failure modes: IndexOutOfRange / negative-length substring on a single/short sequence;
/// off-by-one in the reported arm spans (Length, Start/End); a false repeat reported where no real
/// reverse-complement pairing exists; G-U wobble wrongly accepted (ASM-01 forbids it — G's
/// complement is C, not U).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Theory-correct contract asserted (Inverted_Repeats.md §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   • INV-01 — right arm = reverse complement of left arm:
///       complement(seq[Start2+Length-1-k]) == seq[Start1+k] for all k ∈ [0,Length).
///   • INV-02 — equal arm lengths: End1-Start1+1 == Length == End2-Start2+1.
///   • INV-03 — loop length Start2-End1-1 ∈ [minSpacing, maxSpacing].
///   • INV-04 — Length ≥ minLength.
///   • INV-05 — Start1 ≤ End1 < Start2 ≤ End2 (arms disjoint, left precedes right), all in [0,n).
///   • ASM-01 — strict Watson-Crick/IUPAC complement, NO G-U wobble.
///   • Degenerate / boundary inputs (null/empty/single/short, bad parameters) → empty, never throw.
///   • POSITIVE sanity: a known palindrome/stem is detected with documented arm coordinates;
///     a non-complementary sequence yields none.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaInvertedRepeatFuzzTests
{
    /// <summary>Documented defaults (Inverted_Repeats.md §3.1).</summary>
    private const int DefaultMinLength = 4;
    private const int DefaultMinSpacing = 3;
    private const int DefaultMaxSpacing = 100;

    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static Random Rng(int seed) => new(seed);

    private static string RandomRna(Random rng, int length)
    {
        const string bases = "ACGU";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Reverse complement under the RNA complement table — the same map the detector uses for
    /// arm pairing (<c>GetComplement</c> ⇒ <c>GetRnaComplementBase</c>). Used to build canonical
    /// positives and to independently re-verify reported arms.
    /// </summary>
    private static string ReverseComplement(string s)
    {
        var chars = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            chars[i] = GetComplement(s[s.Length - 1 - i]);
        return new string(chars);
    }

    /// <summary>
    /// Asserts ONE reported inverted repeat is well-formed against the whole documented contract
    /// (Inverted_Repeats.md §2.4 INV-01..05): in-bounds disjoint arms (INV-05), equal arm lengths
    /// matching Length (INV-02), Length ≥ minLength (INV-04), loop length within [minSpacing,
    /// maxSpacing] (INV-03), and EVERY position pairing antiparallel under the strict complement
    /// with NO G-U wobble (INV-01 / ASM-01).
    /// </summary>
    private static void AssertWellFormedRepeat(
        string seq,
        (int Start1, int End1, int Start2, int End2, int Length) r,
        int minLength,
        int minSpacing,
        int maxSpacing)
    {
        int n = seq.Length;

        // INV-05 — ordering and in-bounds.
        r.Start1.Should().BeGreaterThanOrEqualTo(0, "left arm start must be in bounds (INV-05)");
        r.Start1.Should().BeLessThanOrEqualTo(r.End1, "left arm start ≤ end (INV-05)");
        r.End1.Should().BeLessThan(r.Start2, "left arm must end before the right arm begins (INV-05)");
        r.Start2.Should().BeLessThanOrEqualTo(r.End2, "right arm start ≤ end (INV-05)");
        r.End2.Should().BeLessThan(n, "right arm end must be in bounds (INV-05)");

        // INV-02 — equal arm lengths, both equal to Length.
        (r.End1 - r.Start1 + 1).Should().Be(r.Length, "left arm length must equal Length (INV-02)");
        (r.End2 - r.Start2 + 1).Should().Be(r.Length, "right arm length must equal Length (INV-02)");

        // INV-04 — minimum arm length.
        r.Length.Should().BeGreaterThanOrEqualTo(minLength, "arm length must be ≥ minLength (INV-04)");

        // INV-03 — loop length within the window.
        int loop = r.Start2 - r.End1 - 1;
        loop.Should().BeInRange(minSpacing, maxSpacing, "loop length must be in [minSpacing, maxSpacing] (INV-03)");

        // INV-01 / ASM-01 — strict antiparallel reverse-complement pairing, no wobble.
        for (int k = 0; k < r.Length; k++)
        {
            char left = seq[r.Start1 + k];
            char right = seq[r.Start2 + r.Length - 1 - k];
            GetComplement(right).Should().Be(left,
                $"position k={k} must pair antiparallel under the strict complement (INV-01/ASM-01); " +
                $"left='{left}' right='{right}'");
        }
    }

    private static void AssertAllWellFormed(
        string seq,
        IReadOnlyList<(int Start1, int End1, int Start2, int End2, int Length)> repeats,
        int minLength = DefaultMinLength,
        int minSpacing = DefaultMinSpacing,
        int maxSpacing = DefaultMaxSpacing)
    {
        foreach (var r in repeats)
            AssertWellFormedRepeat(seq, r, minLength, minSpacing, maxSpacing);
    }

    #endregion

    #region RNA-INVERT-001 — RNA Inverted-Repeat Detection (perfect W G W̄ᴿ)

    // ─────────────────────────────────────────────────────────────────────────
    // BE: palindrome — the canonical POSITIVE (self-complementary stem)
    // ─────────────────────────────────────────────────────────────────────────

    // The documented worked example (Inverted_Repeats.md §7.1) must be detected with the EXACT
    // coordinates. This pins arm spans, length and complementarity against off-by-one drift.
    [Test]
    public void Palindrome_DocWorkedExample_DetectedWithExactCoordinates()
    {
        var irs = FindInvertedRepeats("UUACGAAAAAACGUAA").ToList();

        irs.Should().ContainSingle("the worked example holds exactly one maximal inverted repeat");
        irs[0].Should().Be((0, 4, 11, 15, 5),
            "left arm [0,4]=UUACG, right arm [11,15]=CGUAA = reverse complement of UUACG (Inverted_Repeats.md §7.1)");
        AssertAllWellFormed("UUACGAAAAAACGUAA", irs);
    }

    // |G| = 0 boundary: a true palindrome with NO spacer (GGGGCCCC) at minSpacing = 0. The two
    // 4-nt arms abut; the right arm is the reverse complement of the left. Exact arm spans —
    // this exercises the smallest-loop boundary without an off-by-one into the loop window.
    [Test]
    public void Palindrome_ZeroSpacer_DetectedAtMinSpacingBoundary()
    {
        // GGGG | (no loop) | CCCC. rev-comp(GGGG) = CCCC.
        var irs = FindInvertedRepeats("GGGGCCCC", minLength: 4, minSpacing: 0, maxSpacing: 100).ToList();

        irs.Should().ContainSingle("a perfect palindrome with no spacer is a single maximal inverted repeat");
        irs[0].Should().Be((0, 3, 4, 7, 4), "arms [0,3]=GGGG and [4,7]=CCCC abut with a zero-length loop");
        AssertAllWellFormed("GGGGCCCC", irs, minLength: 4, minSpacing: 0, maxSpacing: 100);
    }

    // A classic stem (8 nt) + a real loop, built so the reverse complement is unambiguous.
    // The detected stem must span the full complementary arms with the loop in between.
    [Test]
    public void Palindrome_StemPlusLoop_DetectedWithCorrectArmsAndLoop()
    {
        // stem GCGUACGC (rev-comp = GCGUACGC, a self-complementary palindromic stem) | loop AAAAA |
        const string stem = "GCGUACGC";
        const string loop = "AAAAA";
        string seq = stem + loop + ReverseComplement(stem);

        var irs = FindInvertedRepeats(seq, minLength: 4, minSpacing: 3, maxSpacing: 100).ToList();

        irs.Should().ContainSingle();
        var r = irs[0];
        r.Start1.Should().Be(0);
        r.End1.Should().Be(stem.Length - 1);
        r.Length.Should().Be(stem.Length);
        r.Start2.Should().Be(stem.Length + loop.Length);
        r.End2.Should().Be(seq.Length - 1);
        (r.Start2 - r.End1 - 1).Should().Be(loop.Length, "the reported loop must equal the inserted spacer length");
        AssertAllWellFormed(seq, irs);
    }

    // Fuzz: a synthesized stem + loop with a RANDOM stem and a random loop ALWAYS yields a
    // well-formed repeat whose left arm covers the synthesized stem; the detector must find at
    // least the planted stem and never report a malformed/false arm.
    [Test]
    public void Palindrome_RandomSyntheticStems_AlwaysDetectedWellFormed()
    {
        var rng = Rng(151_001);
        for (int trial = 0; trial < 200; trial++)
        {
            int armLen = rng.Next(4, 9);     // ≥ minLength
            int loopLen = rng.Next(3, 12);   // within default window
            string left = RandomRna(rng, armLen);
            string loop = RandomRna(rng, loopLen);
            string seq = left + loop + ReverseComplement(left);

            var irs = FindInvertedRepeats(seq).ToList();

            // A planted perfect stem GUARANTEES at least one inverted repeat exists; the detector
            // must report a non-empty set. (We do NOT pin the planted arm's exact coordinates: the
            // detector reports MAXIMAL, NON-OVERLAPPING stems longest-first, so a longer coincidental
            // stem elsewhere may legitimately crowd out the planted one — INV greedy selection, §4.)
            irs.Should().NotBeEmpty($"a planted stem of {armLen} nt + {loopLen}-nt loop must be detected: {seq}");
            AssertAllWellFormed(seq, irs);

            // Every reported right arm must be the exact reverse complement of its left arm —
            // re-verified independently of the detector's internal pairing loop.
            foreach (var r in irs)
            {
                string detectedLeft = seq.Substring(r.Start1, r.Length);
                string detectedRight = seq.Substring(r.Start2, r.Length);
                detectedRight.Should().Be(ReverseComplement(detectedLeft),
                    $"the reported right arm must be the exact reverse complement of the left arm in '{seq}'");
            }
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE: no complementarity — empty result, no false repeat, no crash
    // ─────────────────────────────────────────────────────────────────────────

    // An all-A homopolymer has NO self-complementary region (A's complement is U, absent), so no
    // inverted repeat exists — the documented "no complement possibilities" case (§6.1).
    [Test]
    public void NoComplementarity_AllAdenine_YieldsEmpty()
    {
        var irs = FindInvertedRepeats(new string('A', 64)).ToList();
        irs.Should().BeEmpty("an all-A sequence has no reverse-complement region, so no inverted repeat (§6.1)");
    }

    // A parallel DIRECT repeat (left == right, 5'→3') is NOT an inverted repeat: the right arm
    // must be the *reverse complement*, not identical (Inverted_Repeats.md §6.1). The arm here is
    // deliberately NOT self-reverse-complementary (rev-comp("AAGCCAA") = "UUGGCUU" ≠ "AAGCCAA"),
    // so the planted direct repeat can never masquerade as an inverted repeat (unlike a palindromic
    // arm, which legitimately is both).
    [Test]
    public void NoComplementarity_DirectRepeat_NotReportedAsInvertedRepeat()
    {
        const string arm = "AAGCCAA";
        arm.Should().NotBe(ReverseComplement(arm), "test premise: the arm must not be self-reverse-complementary");
        const string seq = "AAGCCAA" + "UUUUUUU" + "AAGCCAA"; // direct repeat across a loop

        var irs = FindInvertedRepeats(seq).ToList();

        irs.Should().NotContain(r => seq.Substring(r.Start1, r.Length) == seq.Substring(r.Start2, r.Length),
            "a direct (identical, non-palindromic) repeat is not an inverted repeat — only reverse complements qualify (§6.1)");
        // Whatever (if anything) is reported must still be a genuine reverse-complement repeat.
        AssertAllWellFormed(seq, irs);
    }

    // G-U wobble is NOT accepted (ASM-01): rev-comp(GGGG)=CCCC, NOT UUUU. A GGGG…UUUU layout has
    // no strict Watson-Crick stem, so no inverted repeat is reported there.
    [Test]
    public void NoComplementarity_GuWobble_NotAccepted()
    {
        // If wobble were (wrongly) allowed, GGGG/UUUU would pair. Strict complement forbids it.
        const string seq = "GGGGAAAUUUU";
        var irs = FindInvertedRepeats(seq).ToList();

        irs.Should().NotContain(r =>
                seq.Substring(r.Start1, r.Length).All(c => c == 'G') &&
                seq.Substring(r.Start2, r.Length).All(c => c == 'U'),
            "G-U wobble is not a Watson-Crick complement (ASM-01) — GGGG/UUUU must not pair");
        AssertAllWellFormed(seq, irs);
    }

    // Fuzz: random sequences never produce a malformed or false repeat. Whatever is reported must
    // pass the full INV-01..05 / ASM-01 check; nothing reported is fine too.
    [Test]
    [CancelAfter(30_000)]
    public void NoComplementarity_RandomSequences_NeverReportMalformedOrFalseRepeat()
    {
        var rng = Rng(151_002);
        for (int trial = 0; trial < 400; trial++)
        {
            int len = rng.Next(0, 60);
            string seq = RandomRna(rng, len);

            List<(int, int, int, int, int)> irs = null!;
            Action act = () => irs = FindInvertedRepeats(seq).ToList();
            act.Should().NotThrow($"random RNA must never crash the detector: '{seq}'");

            AssertAllWellFormed(seq, irs);
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE: single base — and the rest of the too-short / degenerate boundary
    // ─────────────────────────────────────────────────────────────────────────

    // A SINGLE base cannot hold two arms + a loop: empty result, NO IndexOutOfRange / negative
    // substring length — the n < 2·minLength + minSpacing guard fires before any indexing.
    [Test]
    public void SingleBase_YieldsEmpty_NoCrash()
    {
        foreach (var seq in new[] { "A", "C", "G", "U", "T", "N", "X" })
        {
            List<(int, int, int, int, int)> irs = null!;
            Action act = () => irs = FindInvertedRepeats(seq).ToList();
            act.Should().NotThrow($"a single-base sequence must not crash: '{seq}'");
            irs.Should().BeEmpty($"'{seq}' is too short to hold two arms and a loop (§3.3, §6.1)");
        }
    }

    // Null / empty are documented to yield an empty result with no throw (§3.3).
    [Test]
    public void NullAndEmpty_YieldEmpty_NoCrash()
    {
        List<(int, int, int, int, int)> fromNull = null!;
        Action actNull = () => fromNull = FindInvertedRepeats(null!).ToList();
        actNull.Should().NotThrow("null input must yield an empty result, not throw (§3.3)");
        fromNull.Should().BeEmpty();

        FindInvertedRepeats(string.Empty).ToList().Should().BeEmpty("empty input yields empty (§3.3)");
    }

    // Sweep every length from 0 up to (and just past) the minimum holdable size for the defaults.
    // For the defaults the smallest sequence that CAN hold a repeat is 2·4+3 = 11; everything
    // below 11 must be empty and never crash, regardless of content.
    [Test]
    public void TooShort_BelowMinimumHoldableLength_AlwaysEmpty_NoCrash()
    {
        var rng = Rng(151_003);
        const int minHoldable = 2 * DefaultMinLength + DefaultMinSpacing; // 11
        for (int len = 0; len < minHoldable; len++)
        {
            for (int trial = 0; trial < 24; trial++)
            {
                string seq = RandomRna(rng, len);
                List<(int, int, int, int, int)> irs = null!;
                Action act = () => irs = FindInvertedRepeats(seq).ToList();
                act.Should().NotThrow($"len={len} '{seq}' must not crash the detector");
                irs.Should().BeEmpty($"len={len} < {minHoldable} cannot hold two {DefaultMinLength}-nt arms + a {DefaultMinSpacing}-nt loop");
            }
        }
    }

    // Even a perfectly self-complementary sequence shorter than the holdable minimum is empty:
    // a 1-arm "palindrome" (e.g. "GC") cannot form a two-arm stem.
    [Test]
    public void TooShort_PerfectButTooSmall_YieldsEmpty()
    {
        // GC is self-complementary (rev-comp(GC)=GC) but only 2 nt — far below 2·4+3.
        FindInvertedRepeats("GC").Should().BeEmpty("a 2-nt self-complementary string cannot hold two min-length arms");
        // Even with minLength=1, minSpacing=0 the holdable minimum is 2 — "G" (len 1) is still empty.
        FindInvertedRepeats("G", minLength: 1, minSpacing: 0).Should().BeEmpty("a single base cannot form two arms even at minLength=1");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // BE: degenerate parameters — bad arguments yield empty, never throw (§3.1, §3.3)
    // ─────────────────────────────────────────────────────────────────────────

    [Test]
    public void DegenerateParameters_OutOfRange_YieldEmpty_NoCrash()
    {
        const string seq = "UUACGAAAAAACGUAA"; // a sequence WITH a real repeat at defaults

        FindInvertedRepeats(seq, minLength: 0).Should().BeEmpty("minLength < 1 → empty (§3.1)");
        FindInvertedRepeats(seq, minLength: -5).Should().BeEmpty("negative minLength → empty (§3.1)");
        FindInvertedRepeats(seq, minSpacing: -1).Should().BeEmpty("minSpacing < 0 → empty (§3.1)");
        FindInvertedRepeats(seq, minSpacing: 5, maxSpacing: 4).Should().BeEmpty("maxSpacing < minSpacing → empty (§3.1)");
    }

    // minLength larger than half the sequence cannot fit two arms: empty, no crash, no negative
    // index (the REP-INV-001 "minLen > seqLen/2" boundary, mirrored for the RNA unit).
    [Test]
    public void DegenerateParameters_MinLengthExceedsHalf_YieldsEmpty()
    {
        const string seq = "GGGGCCCC"; // length 8; a min arm of 5 cannot fit two non-overlapping arms
        FindInvertedRepeats(seq, minLength: 5, minSpacing: 0).Should().BeEmpty(
            "two arms of 5 nt cannot fit in an 8-nt sequence (n < 2·minLength)");
        FindInvertedRepeats(seq, minLength: 100, minSpacing: 0).Should().BeEmpty(
            "an absurdly large minLength yields empty, not a crash");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Determinism / lazy-evaluation discipline
    // ─────────────────────────────────────────────────────────────────────────

    // Repeated enumeration of the same input yields identical results (pure function).
    [Test]
    public void Determinism_RepeatedCalls_AreIdentical()
    {
        var rng = Rng(151_004);
        for (int trial = 0; trial < 50; trial++)
        {
            string seq = RandomRna(rng, rng.Next(12, 40));
            var a = FindInvertedRepeats(seq).ToList();
            var b = FindInvertedRepeats(seq).ToList();
            b.Should().Equal(a, $"FindInvertedRepeats must be deterministic for '{seq}'");
            AssertAllWellFormed(seq, a);
        }
    }

    #endregion
}
