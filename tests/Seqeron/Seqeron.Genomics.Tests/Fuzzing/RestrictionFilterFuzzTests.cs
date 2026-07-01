using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using FluentAssertions.Execution;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for RESTR-FILTER-001 — Restriction Enzyme Filtering
/// (RestrictionAnalyzer.GetEnzymesByCutLength / GetBluntCutters / GetStickyCutters).
/// Checklist: docs/checklists/03_FUZZING.md, row 224 (MolTools area, strategy BE).
/// Algorithm doc: docs/algorithms/MolTools/Restriction_Enzyme_Filtering.md.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed, out-of-domain and boundary inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no crash, no hang/infinite
/// loop, no NaN/Infinity, no state corruption, and no *unhandled* runtime exception
/// leaking from internal logic. Every input must resolve to a well-defined,
/// theory-correct result. For this unit the "fuzz surface" is unusual: the filters
/// take NO sequence input — they are pure set operations over a fixed, curated
/// built-in enzyme library (doc §1, §2.4, §3.3). The only attacker-controlled inputs
/// are the integer length bounds of GetEnzymesByCutLength. So Boundary Exploitation
/// (BE = 0, -1, int.MaxValue, int.MinValue, empty interval) is applied to those
/// integers, and the blunt/sticky/length predicate contract is pinned exactly.
///
/// The doc (§3.3) is explicit: all four methods are TOTAL — they never throw and
/// never return null. A non-positive or inverted range yields an empty sequence,
/// not an exception. Anything else is a bug.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The filter contract under test — Restriction_Enzyme_Filtering.md §2.2/§2.4
/// ───────────────────────────────────────────────────────────────────────────
/// For an enzyme e with recognition sequence r(e) and per-strand cut offsets
/// cf(e)/cr(e) measured from the 5' end of the site:
///   • len(e) = |r(e)|.
///   • e is BLUNT iff cf(e) = cr(e) (center/symmetric cut → both strands terminate
///     in a base pair); otherwise e is STICKY (staggered cut → 5'/3' overhang). [§2.2]
///   • ByLength(min,max) = { e : min ≤ len(e) ≤ max }, BOTH bounds inclusive. [INV-03]
///   • ByLength(L) = ByLength(L,L). [INV-04]
///   • min > max ⇒ ByLength(min,max) = ∅. [INV-05]
///   • Blunt ∩ Sticky = ∅ and Blunt ∪ Sticky = Library (total partition). [INV-01]
///
/// ───────────────────────────────────────────────────────────────────────────
/// Boundaries from the checklist row hint ("no sites, all-pass, all-fail")
/// ───────────────────────────────────────────────────────────────────────────
///   • "no sites"  → the length predicate matches NOTHING (e.g. [9,10] above the
///                   8-nt max for undivided sites, or non-positive bounds) ⇒ empty.
///   • "all-pass"  → the predicate matches EVERY candidate (a range spanning the
///                   whole library, including via int.MinValue..int.MaxValue) ⇒
///                   the entire library is returned.
///   • "all-fail"  → no candidate satisfies the predicate ⇒ empty result.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Independently-derived expected values (NOT echoed off the code's arrays)
/// ───────────────────────────────────────────────────────────────────────────
/// The library is enumerated against the doc (§2.1, §4.2, §6.1, §7.1) and primary
/// sources (Wikipedia "Restriction enzyme"/"Sticky and blunt ends"/"List of …
/// cutting sites"; NEB/REBASE). Recognition-length distribution: 4-nt and 6-nt and
/// 8-nt undivided Type-II sites, plus exactly one 13-nt interrupted palindrome SfiI
/// (GGCCNNNN^NGGCC, PMC548270) which lies OUTSIDE the 4–8 nt range (doc §6.1).
/// Therefore ByLength(4,8) = full library minus SfiI (doc §6.1 / §7.1).
/// The blunt (center-cut) enzymes, each cut site cross-checked against REBASE/NEB:
///   AluI AG^CT, DpnI Gm6A^TC, EcoRV GAT^ATC, HaeIII GG^CC, HincII GTY^RAC,
///   RsaI GT^AC, ScaI AGT^ACT, SmaI CCC^GGG, StuI AGG^CCT, SwaI ATTT^AAAT.
/// Representative sticky producers: EcoRI G^AATTC (5' overhang), KpnI GGTAC^C and
/// PstI CTGCA^G (3' overhangs). These ground the positive-sanity assertions below;
/// a wrong cut position for ANY one enzyme breaks the exact-set checks.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RestrictionFilterFuzzTests
{
    #region Helpers

    private const int Seed = 224_0001; // local-only deterministic seed for the fuzz sweep

    /// <summary>
    /// The blunt (center-cut, cf == cr) enzymes, derived independently from REBASE/NEB/
    /// Wikipedia — NOT read off the code's cut-position arrays. Used to pin the exact
    /// blunt/sticky partition.
    /// </summary>
    private static readonly string[] ExpectedBlunt =
    {
        "AluI", "DpnI", "EcoRV", "HaeIII", "HincII", "RsaI", "ScaI", "SmaI", "StuI", "SwaI"
    };

    private static IReadOnlyList<RestrictionEnzyme> Library =>
        RestrictionAnalyzer.Enzymes.Values.ToList();

    private static HashSet<string> Names(IEnumerable<RestrictionEnzyme> enzymes) =>
        enzymes.Select(e => e.Name).ToHashSet(StringComparer.Ordinal);

    /// <summary>
    /// Every returned enzyme must be a real library member with a coherent record:
    /// a non-empty recognition sequence, a recognition length equal to that string's
    /// length, and a blunt flag consistent with cf == cr. No corruption / phantom rows.
    /// </summary>
    private static void AssertWellFormed(IEnumerable<RestrictionEnzyme> result)
    {
        var libraryNames = Names(Library);
        foreach (var e in result)
        {
            e.Should().NotBeNull();
            e.Name.Should().NotBeNullOrEmpty();
            libraryNames.Should().Contain(e.Name, "filters may only return enzymes from the built-in library");
            e.RecognitionSequence.Should().NotBeNullOrEmpty();
            e.RecognitionLength.Should().Be(e.RecognitionSequence.Length,
                "RecognitionLength is the recognition-string length (doc §2.2)");
            e.IsBluntEnd.Should().Be(e.CutPositionForward == e.CutPositionReverse,
                "blunt iff the two strands are cut at the same offset (doc §2.2, §5.2)");
        }
    }

    #endregion

    #region RESTR-FILTER-001 — Restriction enzyme filtering

    #region Positive sanity — hand-checkable worked example (doc §7.1, §6.1)

    // Reproduces the doc §7.1 worked example: GetBluntCutters() returns exactly the
    // center-cut producers; GetEnzymesByCutLength(6,6) returns only 6-cutters; and
    // GetEnzymesByCutLength(4,8) returns every undivided site, i.e. the whole library
    // minus the single 13-nt interrupted palindrome SfiI (doc §6.1).
    [Test]
    public void WorkedExample_BluntSet_SixCutters_AndUndividedRange_MatchDoc()
    {
        var blunt = Names(RestrictionAnalyzer.GetBluntCutters());
        var sixCutters = RestrictionAnalyzer.GetEnzymesByCutLength(6, 6).ToList();
        var undivided = RestrictionAnalyzer.GetEnzymesByCutLength(4, 8).ToList();

        using (new AssertionScope())
        {
            // Blunt = exactly the externally-sourced center-cut set (no overhang producer leaks in).
            blunt.Should().BeEquivalentTo(ExpectedBlunt,
                "GetBluntCutters must return exactly the documented center-cut (blunt) enzymes (doc §2.2, INV-02)");

            // 6-cutters: every result is 6 nt and includes the canonical examples.
            sixCutters.Should().OnlyContain(e => e.RecognitionLength == 6);
            Names(sixCutters).Should().Contain(new[] { "EcoRI", "BamHI", "PstI" });
            Names(sixCutters).Should().NotContain("AluI");  // 4-cutter
            Names(sixCutters).Should().NotContain("NotI");  // 8-cutter

            // [4,8] = the full library except the lone 13-nt interrupted palindrome SfiI (doc §6.1).
            undivided.Should().HaveCount(Library.Count - 1);
            Names(undivided).Should().NotContain("SfiI",
                "SfiI (GGCCNNNN^NGGCC, length 13) lies outside the 4–8 nt undivided range (doc §6.1)");
            undivided.Should().OnlyContain(e => e.RecognitionLength >= 4 && e.RecognitionLength <= 8);
        }
    }

    // Sticky set is the exact complement of blunt: contains the documented overhang
    // producers and excludes every blunt cutter (doc §2.2, §4.2).
    [Test]
    public void WorkedExample_StickySet_IsComplementOfBlunt()
    {
        var sticky = Names(RestrictionAnalyzer.GetStickyCutters());

        using (new AssertionScope())
        {
            sticky.Should().Contain(new[] { "EcoRI", "KpnI", "PstI" },
                "EcoRI (5' overhang), KpnI/PstI (3' overhangs) are staggered cuts → sticky");
            sticky.Should().NotContain(ExpectedBlunt,
                "no blunt (center-cut) enzyme may appear in the sticky set");
        }
    }

    #endregion

    #region BE — "all-pass": a range spanning the whole library returns everything

    // A range whose bounds bracket every recognition length (1..int.MaxValue, and the
    // extreme int.MinValue..int.MaxValue) is the "all-pass" boundary: the predicate
    // is satisfied by EVERY enzyme, so the entire library — SfiI included — is returned.
    [Test]
    [CancelAfter(20000)]
    public void ByLength_AllPassRange_ReturnsEntireLibrary()
    {
        var widest = RestrictionAnalyzer.GetEnzymesByCutLength(1, int.MaxValue).ToList();
        var extreme = RestrictionAnalyzer.GetEnzymesByCutLength(int.MinValue, int.MaxValue).ToList();

        using (new AssertionScope())
        {
            AssertWellFormed(widest);
            AssertWellFormed(extreme);

            // all-pass ⇒ count == library size, and the 13-nt SfiI is now INCLUDED
            // (unlike the [4,8] undivided range) because no upper bound excludes it.
            widest.Should().HaveCount(Library.Count);
            Names(widest).Should().BeEquivalentTo(Names(Library));
            Names(widest).Should().Contain("SfiI", "an unbounded upper limit includes the 13-nt site");

            extreme.Should().HaveCount(Library.Count,
                "int.MinValue..int.MaxValue must not overflow or drop any enzyme");
            Names(extreme).Should().BeEquivalentTo(Names(Library));
        }
    }

    #endregion

    #region BE — "all-fail" / "no sites": empty-result boundaries

    // "no sites": a range entirely above the 8-nt maximum for undivided sites returns
    // nothing — but only because SfiI (13 nt) does not fall in [9,12]. Two empty cases:
    //   [9,10]  — above the undivided max, below SfiI's 13 (doc §6.1)
    //   [0,0]   — no site has length 0
    //   [-1,0]  — non-positive bounds (doc §6.1)
    // and the int.MaxValue/int.MinValue degenerate single-points.
    [Test]
    [CancelAfter(20000)]
    public void ByLength_NoMatchingLength_ReturnsEmpty()
    {
        using (new AssertionScope())
        {
            RestrictionAnalyzer.GetEnzymesByCutLength(9, 10).Should().BeEmpty(
                "[9,10] is above the 8-nt undivided max and below SfiI's 13 nt (doc §6.1)");
            RestrictionAnalyzer.GetEnzymesByCutLength(0, 0).Should().BeEmpty(
                "no recognition site has length 0");
            RestrictionAnalyzer.GetEnzymesByCutLength(-1, 0).Should().BeEmpty(
                "non-positive bounds match nothing (doc §6.1)");
            RestrictionAnalyzer.GetEnzymesByCutLength(int.MaxValue, int.MaxValue).Should().BeEmpty(
                "no site is int.MaxValue nt long");
            RestrictionAnalyzer.GetEnzymesByCutLength(int.MinValue, int.MinValue).Should().BeEmpty(
                "no site is int.MinValue nt long");
            RestrictionAnalyzer.GetEnzymesByCutLength(3).Should().BeEmpty(
                "the smallest recognition site is 4 nt; no 3-cutter exists in the library");
        }
    }

    // INV-05: an inverted range (min > max), including extreme inversions, describes an
    // empty interval ⇒ empty result, never an exception (doc §3.3, §6.1).
    [Test]
    [CancelAfter(20000)]
    public void ByLength_InvertedRange_ReturnsEmpty_NeverThrows()
    {
        using (new AssertionScope())
        {
            FluentActions.Invoking(() => RestrictionAnalyzer.GetEnzymesByCutLength(8, 4).ToList())
                .Should().NotThrow();
            RestrictionAnalyzer.GetEnzymesByCutLength(8, 4).Should().BeEmpty();
            RestrictionAnalyzer.GetEnzymesByCutLength(int.MaxValue, int.MinValue).Should().BeEmpty(
                "the maximally inverted range is still an empty interval");
            RestrictionAnalyzer.GetEnzymesByCutLength(6, 5).Should().BeEmpty();
        }
    }

    #endregion

    #region BE — single-length boundaries pin exact per-length sets (INV-03/INV-04)

    // Exact per-length partition for the populated lengths 4/6/8, plus INV-04 equality of
    // the single-length overload with the equal-bounds range. Counts derived independently:
    // the library has 9 four-cutters, 24 six-cutters, 5 eight-cutters and 1 thirteen-cutter
    // (SfiI), summing to the full library; 4+6+8 alone exclude SfiI.
    [Test]
    [CancelAfter(20000)]
    public void ByLength_ExactLengths_PartitionTheUndividedSites()
    {
        var four = RestrictionAnalyzer.GetEnzymesByCutLength(4).ToList();
        var six = RestrictionAnalyzer.GetEnzymesByCutLength(6).ToList();
        var eight = RestrictionAnalyzer.GetEnzymesByCutLength(8).ToList();

        using (new AssertionScope())
        {
            four.Should().OnlyContain(e => e.RecognitionLength == 4);
            six.Should().OnlyContain(e => e.RecognitionLength == 6);
            eight.Should().OnlyContain(e => e.RecognitionLength == 8);

            Names(four).Should().Contain(new[] { "AluI", "HaeIII", "TaqI" });
            Names(six).Should().Contain(new[] { "EcoRI", "BamHI", "PstI" });
            Names(eight).Should().Contain(new[] { "NotI", "PacI", "AscI" });

            // 4+6+8 together = every undivided site = [4,8] = library minus SfiI.
            (four.Count + six.Count + eight.Count).Should().Be(Library.Count - 1,
                "the only library enzyme outside lengths {4,6,8} is the 13-nt SfiI");

            // INV-04: single-length overload == equal-bounds range overload.
            foreach (int length in new[] { 4, 6, 8 })
            {
                Names(RestrictionAnalyzer.GetEnzymesByCutLength(length))
                    .Should().BeEquivalentTo(Names(RestrictionAnalyzer.GetEnzymesByCutLength(length, length)),
                        $"ByLength({length}) must equal ByLength({length},{length}) (INV-04)");
            }
        }
    }

    #endregion

    #region BE — blunt/sticky total partition (INV-01) and totality (§3.3)

    // INV-01: blunt and sticky are disjoint and together cover the whole library; the
    // blunt set is exactly the externally-sourced center-cut set. Both filters are total
    // (never null, never throw) per §3.3 — the "all-pass on one side, all-fail on the
    // other" structure of a Boolean split.
    [Test]
    [CancelAfter(20000)]
    public void BluntSticky_TotalPartition_AndTotality()
    {
        IEnumerable<RestrictionEnzyme>? blunt = null;
        IEnumerable<RestrictionEnzyme>? sticky = null;

        using (new AssertionScope())
        {
            FluentActions.Invoking(() =>
            {
                blunt = RestrictionAnalyzer.GetBluntCutters();
                sticky = RestrictionAnalyzer.GetStickyCutters();
            }).Should().NotThrow();

            blunt.Should().NotBeNull();
            sticky.Should().NotBeNull();

            var bluntNames = Names(blunt!);
            var stickyNames = Names(sticky!);
            var libraryNames = Names(Library);

            AssertWellFormed(blunt!);
            AssertWellFormed(sticky!);

            bluntNames.Should().BeEquivalentTo(ExpectedBlunt);
            bluntNames.Overlaps(stickyNames).Should().BeFalse("blunt ∩ sticky = ∅ (INV-01)");
            (bluntNames.Count + stickyNames.Count).Should().Be(libraryNames.Count,
                "blunt + sticky counts sum to the full library (total partition, INV-01)");
            bluntNames.Union(stickyNames).Should().BeEquivalentTo(libraryNames,
                "blunt ∪ sticky = library (INV-01)");
        }
    }

    #endregion

    #region Randomized BE sweep — never crash / hang / corrupt over fuzzed int bounds

    // Boundary sweep over random integer bounds drawn from the BE pool
    // {int.MinValue, -1, 0, 1, 3..14, int.MaxValue}. For EVERY (min,max) pair the
    // filter must: not throw, not hang, and return a result that exactly equals the
    // independently-recomputed predicate { e : min ≤ len(e) ≤ max } over the library —
    // including the empty set when the interval is empty/inverted. This pins the real
    // algorithmic contract on fuzzed input, not merely "doesn't crash".
    [Test]
    [CancelAfter(60000)]
    public void ByLength_RandomizedBoundarySweep_MatchesIndependentPredicate()
    {
        var rng = new Random(Seed);
        int[] pool = { int.MinValue, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, int.MaxValue };
        var library = Library;

        for (int iter = 0; iter < 4000; iter++)
        {
            int min = pool[rng.Next(pool.Length)];
            int max = pool[rng.Next(pool.Length)];

            List<RestrictionEnzyme> actual = null!;
            Action act = () => actual = RestrictionAnalyzer.GetEnzymesByCutLength(min, max).ToList();
            act.Should().NotThrow($"ByLength({min},{max}) is total (doc §3.3)");

            AssertWellFormed(actual);

            // Independent recomputation of the inclusive-range predicate over the library.
            // long arithmetic avoids any int overflow when comparing against int.MinValue/MaxValue.
            var expected = library
                .Where(e => (long)e.RecognitionLength >= min && (long)e.RecognitionLength <= max)
                .Select(e => e.Name)
                .ToHashSet(StringComparer.Ordinal);

            Names(actual).Should().BeEquivalentTo(expected,
                $"ByLength({min},{max}) must equal {{ e : {min} ≤ len(e) ≤ {max} }} (INV-03/INV-05)");

            if (min > max)
                actual.Should().BeEmpty($"inverted range [{min},{max}] is empty (INV-05)");
        }
    }

    #endregion

    #endregion
}
