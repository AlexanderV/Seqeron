using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the RnaStructure area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: RNA-STRUCT-001 — secondary-structure prediction (RnaStructure).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 71.
///
/// API under test (RnaSecondaryStructure.PredictStructure):
///   Finds complementary stem-loops, greedily selects a non-overlapping set by energy, and
///   reports the resulting base pairs (i,j), the dot-bracket string and the total MFE.
///   Base pairing requires Watson–Crick / wobble complementarity, a minimum stem length and
///   a loop within [minLoopSize, maxLoopSize].
///
/// Relations (derived from complementary pairing, NOT from output):
///   • COMP (empty ⇒ 0 pairs): an empty sequence — and any sequence with no complementary
///          run able to close a stem (e.g. a homopolymer) — yields no base pairs.
///   • MON  (more complementary bases ⇒ ≥ pairs): lengthening the complementary arms of a
///          hairpin can only add base pairs, so the pair count is non-decreasing in arm length.
///   • INV  (non-pairing 3' tail ⇒ existing pairs preserved): appending a tail that cannot
///          pair (a homopolymer) at the 3' end leaves every base pair of the original hairpin
///          intact at the same positions.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class RnaStructureMetamorphicTests
{
    #region Helpers

    /// <summary>The set of base-pair index pairs (i,j) predicted for a sequence.</summary>
    private static HashSet<(int I, int J)> PairSet(string rna) =>
        RnaSecondaryStructure.PredictStructure(rna)
            .BasePairs.Select(bp => (bp.Position1, bp.Position2))
            .ToHashSet();

    /// <summary>A perfect GC hairpin: an n-base G arm, a 4-base poly-A loop, an n-base C arm.</summary>
    private static string GcHairpin(int armLength) =>
        new string('G', armLength) + "AAAA" + new string('C', armLength);

    #endregion

    #region COMP — no complementarity ⇒ no base pairs

    [Test]
    [Description("COMP: an empty sequence yields no base pairs, an empty dot-bracket and zero MFE.")]
    public void PredictStructure_EmptySequence_HasNoPairs()
    {
        var result = RnaSecondaryStructure.PredictStructure("");

        result.BasePairs.Should().BeEmpty(because: "there are no bases to pair");
        result.DotBracket.Should().BeEmpty(because: "an empty sequence has an empty structure string");
        result.MinimumFreeEnergy.Should().Be(0, because: "no pairs contribute any stacking energy");
    }

    [Test]
    [Description("COMP: a homopolymer cannot form any complementary stem, so no base pairs are predicted.")]
    public void PredictStructure_Homopolymer_HasNoPairs()
    {
        foreach (var homopolymer in new[] { "AAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCC" })
            RnaSecondaryStructure.PredictStructure(homopolymer)
                .BasePairs.Should().BeEmpty(
                    because: $"'{homopolymer[0]}' cannot pair with itself, so no stem can close");
    }

    #endregion

    #region MON — longer complementary arms cannot reduce the pair count

    [Test]
    [Description("MON: extending the complementary arms of a hairpin can only add base pairs — the pair count is non-decreasing in arm length and strictly grows over the range.")]
    public void PredictStructure_LongerComplementaryArms_DoNotReducePairs()
    {
        int previous = -1;
        var counts = new List<int>();

        foreach (int arm in new[] { 3, 4, 5, 6, 7 })
        {
            int pairs = RnaSecondaryStructure.PredictStructure(GcHairpin(arm)).BasePairs.Count;
            counts.Add(pairs);

            pairs.Should().BeGreaterThanOrEqualTo(previous,
                because: $"a longer ({arm}-base) complementary arm keeps every pair the shorter arm could form and may add more");
            previous = pairs;
        }

        counts.Last().Should().BeGreaterThan(counts.First(),
            because: "the 7-base stem forms strictly more base pairs than the 3-base stem");
    }

    #endregion

    #region INV — a non-pairing 3' tail preserves existing pairs

    [Test]
    [Description("INV: appending a homopolymer tail (which cannot pair) at the 3' end leaves every base pair of the original hairpin intact at the same positions.")]
    public void PredictStructure_NonPairing3PrimeTail_PreservesExistingPairs()
    {
        const string core = "GGGGGAAAACCCCC"; // a 5-bp GC hairpin closing a poly-A loop
        var corePairs = PairSet(core);
        corePairs.Should().NotBeEmpty(because: "the core hairpin must actually fold for the test to be meaningful");

        // A poly-A tail at the 3' end cannot pair (no U to pair with) and does not shift the
        // hairpin's positions, so the original pairs must survive unchanged.
        foreach (var tail in new[] { "A", "AAAA", "AAAAAAAA" })
        {
            var extendedPairs = PairSet(core + tail);

            extendedPairs.IsSupersetOf(corePairs).Should().BeTrue(
                because: $"a non-pairing 3' tail of {tail.Length} A's cannot break or move the existing base pairs");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: RNA-STEMLOOP-001 — stem-loop (hairpin) detection (RnaStructure).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 72.
    //
    // API under test (RnaSecondaryStructure.FindStemLoops):
    //   For each loop window, extends a stem outward as long as the flanking bases are
    //   complementary (Watson–Crick or wobble), returning a StemLoop whose Stem.Length is the
    //   number of consecutive complementary pairs. The stem is built only from the arms; the
    //   loop bases lie strictly between the arms and are never part of the stem.
    //
    // Relations (derived from outward complementary extension, NOT from output):
    //   • COMP (no complement ⇒ no stem): a sequence with no complementary arms (a homopolymer,
    //          or like-base "arms") yields no stem-loops.
    //   • MON  (longer arms ⇒ longer stem): lengthening the complementary arms of a hairpin
    //          lengthens the maximal stem — its length is non-decreasing and strictly grows.
    //   • INV  (loop content irrelevant to stem pairing): with the arms and loop length fixed,
    //          changing the loop bases leaves the hairpin's stem base-pairs identical, since the
    //          stem is built from the arms alone.
    // ───────────────────────────────────────────────────────────────────────────

    #region RNA-STEMLOOP-001 — Helpers

    private static HashSet<(int I, int J)> StemPairSet(RnaSecondaryStructure.StemLoop sl) =>
        sl.Stem.BasePairs.Select(bp => (bp.Position1, bp.Position2)).ToHashSet();

    #endregion

    #region RNA-STEMLOOP-001 COMP — no complementary arms ⇒ no stem-loop

    [Test]
    [Description("COMP: a sequence with no complementary arms (a homopolymer or like-base arms) produces no stem-loops.")]
    public void FindStemLoops_NoComplementarity_ReturnsNothing()
    {
        foreach (var seq in new[] { "AAAAAAAAAAAA", "GGGGAAAGGGG" })
            RnaSecondaryStructure.FindStemLoops(seq).Should().BeEmpty(
                because: $"'{seq}' has no two complementary arms able to close a stem");
    }

    #endregion

    #region RNA-STEMLOOP-001 MON — longer complementary arms lengthen the maximal stem

    [Test]
    [Description("MON: lengthening a hairpin's complementary arms lengthens the maximal stem — non-decreasing in arm length and strictly larger over the range.")]
    public void FindStemLoops_LongerArms_LengthenMaximalStem()
    {
        int previous = -1;
        var maxima = new List<int>();

        foreach (int arm in new[] { 3, 4, 5, 6, 7 })
        {
            var stemLoops = RnaSecondaryStructure.FindStemLoops(GcHairpin(arm)).ToList();
            int maxStem = stemLoops.Count == 0 ? 0 : stemLoops.Max(sl => sl.Stem.Length);
            maxima.Add(maxStem);

            maxStem.Should().BeGreaterThanOrEqualTo(previous,
                because: $"a longer ({arm}-base) complementary arm cannot shorten the maximal stem");
            previous = maxStem;
        }

        maxima.Last().Should().BeGreaterThan(maxima.First(),
            because: "the 7-base arms yield a strictly longer maximal stem than the 3-base arms");
    }

    #endregion

    #region RNA-STEMLOOP-001 INV — loop content does not change the stem base-pairs

    [Test]
    [Description("INV: with the arms (5×G / 5×C) and loop length (4) fixed, changing the loop bases leaves the hairpin's stem base-pairs identical, because the stem is built from the arms alone.")]
    public void FindStemLoops_LoopContent_DoesNotChangeStemPairing()
    {
        // Arms G(0..4) / C(9..13) close a length-4 loop at 5..8 → stem pairs (i, 13-i), i=0..4.
        var expected = Enumerable.Range(0, 5).Select(i => (i, 13 - i)).ToHashSet();

        foreach (var loop in new[] { "AAAA", "GAAA", "UUCG", "CUUG" })
        {
            string seq = "GGGGG" + loop + "CCCCC";
            var stemLoops = RnaSecondaryStructure.FindStemLoops(seq).ToList();

            stemLoops.Any(sl => StemPairSet(sl).SetEquals(expected)).Should().BeTrue(
                because: $"the 5-bp G:C stem is determined by the arms, so loop '{loop}' yields the same stem base-pairs");
            stemLoops.Max(sl => sl.Stem.Length).Should().Be(5,
                because: $"the arms are 5 bases long, so the maximal stem is length 5 regardless of loop '{loop}'");
        }
    }

    #endregion
}
