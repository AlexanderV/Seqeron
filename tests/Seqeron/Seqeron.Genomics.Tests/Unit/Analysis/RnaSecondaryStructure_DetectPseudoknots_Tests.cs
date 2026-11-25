// RNA-PSEUDOKNOT-001 — Pseudoknot Detection
// Evidence: docs/Evidence/RNA-PSEUDOKNOT-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-PSEUDOKNOT-001.md
// Source: Antczak M et al. (2018). Bioinformatics 34(8):1304-1312 (crossing condition i<k<j<l,
//         "([)]" example); Smit S et al. (2008). RNA 14(3):410-416 (pseudoknot = crossing pairs);
//         biotite.structure.pseudoknots; Wikipedia "Pseudoknot" (Rivas & Eddy 1999).

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_DetectPseudoknots_Tests
{
    private static BasePair Pair(int p1, int p2) =>
        new(p1, p2, 'A', 'U', BasePairType.WatsonCrick);

    #region DetectPseudoknots

    // M1 — H-type crossing "([)]": pairs (0,2)+(1,3) satisfy i<k<j<l (0<1<2<3) → exactly one pseudoknot.
    // Evidence: Antczak (2018) crossing condition i<k<j<l and the "([)]" DBL example.
    [Test]
    public void DetectPseudoknots_HTypeCrossingPairs_ReturnsOnePseudoknot()
    {
        var pairs = new List<BasePair> { Pair(0, 2), Pair(1, 3) };

        var result = DetectPseudoknots(pairs).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1),
                "(0,2) and (1,3) cross (0<1<2<3) — exactly one pseudoknot per Antczak (2018)");
            Assert.That(result[0].CrossingPairs, Has.Count.EqualTo(2),
                "a pseudoknot is formed by the two crossing base pairs");
            Assert.That(result[0].CrossingPairs, Does.Contain(Pair(0, 2)),
                "the first crossing pair (0,2) must be carried in the result");
            Assert.That(result[0].CrossingPairs, Does.Contain(Pair(1, 3)),
                "the second crossing pair (1,3) must be carried in the result");
        });
    }

    // M2 — Nested control: (0,5) fully contains (1,4) (i<k<l<j) → not well nested is FALSE → zero pseudoknots.
    // Evidence: Antczak (2018) nested condition i<k<l<j; Wikipedia "not well nested".
    [Test]
    public void DetectPseudoknots_NestedPairs_ReturnsNone()
    {
        var pairs = new List<BasePair> { Pair(0, 5), Pair(1, 4) };

        var result = DetectPseudoknots(pairs).ToList();

        Assert.That(result, Is.Empty,
            "(0,5) fully contains (1,4): 0<1<4<5 is nested, not crossing — no pseudoknot");
    }

    // M3 — Disjoint control: (0,2) and (3,5) do not overlap (j<k) → zero pseudoknots.
    // Evidence: crossing requires range overlap; disjoint pairs are side-by-side.
    [Test]
    public void DetectPseudoknots_DisjointPairs_ReturnsNone()
    {
        var pairs = new List<BasePair> { Pair(0, 2), Pair(3, 5) };

        var result = DetectPseudoknots(pairs).ToList();

        Assert.That(result, Is.Empty,
            "(0,2) ends before (3,5) begins (2<3) — disjoint, no crossing");
    }

    // M4 — Endpoint normalization: same crossing as M1 but endpoints stored as (close,open).
    // Evidence: biotite — crossing defined on open<close, so pairs must be min/max normalized first.
    [Test]
    public void DetectPseudoknots_ReversedEndpoints_DetectsSameCrossing()
    {
        var pairs = new List<BasePair> { Pair(2, 0), Pair(3, 1) };

        var result = DetectPseudoknots(pairs).ToList();

        Assert.That(result, Has.Count.EqualTo(1),
            "(2,0) and (3,1) normalize to (0,2) and (1,3) which cross — must detect the same pseudoknot as M1");
    }

    // M5 — Reported coordinates: result exposes normalized, open-first crossing endpoints.
    // Evidence: crossing definition Start1<Start2<End1<End2 (i<k<j<l).
    [Test]
    public void DetectPseudoknots_CrossingPair_ReportsNormalizedCoordinates()
    {
        var pairs = new List<BasePair> { Pair(0, 2), Pair(1, 3) };

        var knot = DetectPseudoknots(pairs).Single();

        Assert.Multiple(() =>
        {
            Assert.That(knot.Start1, Is.EqualTo(0), "i = opening of the first-opening pair (0)");
            Assert.That(knot.Start2, Is.EqualTo(1), "k = opening of the second pair (1)");
            Assert.That(knot.End1, Is.EqualTo(2), "j = closing of the first-opening pair (2)");
            Assert.That(knot.End2, Is.EqualTo(3), "l = closing of the second pair (3)");
            Assert.That(knot.Start1 < knot.Start2 && knot.Start2 < knot.End1 && knot.End1 < knot.End2,
                Is.True, "reported coordinates must satisfy i<k<j<l");
        });
    }

    // S1 — Empty input → no pseudoknots (a pseudoknot needs two crossing pairs).
    // Evidence: derived from crossing definition; Smit (2008) (removal requires crossing pairs).
    [Test]
    public void DetectPseudoknots_EmptyInput_ReturnsNone()
    {
        var result = DetectPseudoknots(new List<BasePair>()).ToList();

        Assert.That(result, Is.Empty, "an empty pair set cannot contain a crossing");
    }

    // S2 — Single pair → no pseudoknots (one pair cannot cross anything).
    // Evidence: derived corner case (>=2 pairs needed to cross).
    [Test]
    public void DetectPseudoknots_SinglePair_ReturnsNone()
    {
        var result = DetectPseudoknots(new List<BasePair> { Pair(0, 4) }).ToList();

        Assert.That(result, Is.Empty, "a single base pair has nothing to cross with");
    }

    // S2b — Null input → no pseudoknots (documented contract, no exception).
    [Test]
    public void DetectPseudoknots_NullInput_ReturnsNone()
    {
        var result = DetectPseudoknots(null!).ToList();

        Assert.That(result, Is.Empty, "null pair set is treated as empty per the documented contract");
    }

    // S3 — Mixed set: nested + crossing + disjoint together → only the crossing relation is reported.
    // Evidence: INV-01..INV-03 applied independently to each pair-of-pairs.
    [Test]
    public void DetectPseudoknots_MixedPairs_ReportsOnlyCrossing()
    {
        // (0,9) nests (1,8); (2,3) is disjoint from all; (11,14) crosses (12,16).
        // All of {0,9},{1,8},{2,3} lie within [0,9] and so are disjoint from the crossing block (>=11).
        var pairs = new List<BasePair>
        {
            Pair(0, 9), Pair(1, 8), Pair(2, 3), Pair(11, 14), Pair(12, 16),
        };

        var result = DetectPseudoknots(pairs).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1),
                "only (11,14) and (12,16) cross (11<12<14<16); the nested and disjoint pairs do not");
            Assert.That(result[0].CrossingPairs, Does.Contain(Pair(11, 14)),
                "the crossing result must carry pair (11,14)");
            Assert.That(result[0].CrossingPairs, Does.Contain(Pair(12, 16)),
                "the crossing result must carry pair (12,16)");
        });
    }

    // S4 — Order independence: shuffling the input list yields the identical set of crossings.
    // Evidence: INV-05 — pure combinatorial scan is order-independent.
    [Test]
    public void DetectPseudoknots_ReorderedInput_ProducesSameCrossingSet()
    {
        var ordered = new List<BasePair> { Pair(0, 2), Pair(1, 3), Pair(5, 9), Pair(7, 11) };
        var shuffled = new List<BasePair> { Pair(7, 11), Pair(0, 2), Pair(5, 9), Pair(1, 3) };

        var fromOrdered = DetectPseudoknots(ordered)
            .Select(k => (k.Start1, k.End1, k.Start2, k.End2)).OrderBy(t => t).ToList();
        var fromShuffled = DetectPseudoknots(shuffled)
            .Select(k => (k.Start1, k.End1, k.Start2, k.End2)).OrderBy(t => t).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(fromOrdered, Has.Count.EqualTo(2),
                "two crossings present: (0,2)x(1,3) and (5,9)x(7,11)");
            Assert.That(fromShuffled, Is.EqualTo(fromOrdered),
                "detected crossing set must not depend on input ordering");
        });
    }

    // S5 — Three mutually-crossing pairs "([{)]}": (0,3),(1,4),(2,5). Each pair crosses each other
    // (0<1<3<4, 0<2<3<5, 1<2<4<5), so all C(3,2)=3 pairwise crossings are reported separately.
    // Evidence: Antczak (2018) crossing condition i<k<j<l applied to every pair-of-pairs; the
    // documented contract reports each crossing pair-of-pairs as its own Pseudoknot (no order grouping).
    [Test]
    public void DetectPseudoknots_ThreeMutuallyCrossingPairs_ReportsEachPairwiseCrossing()
    {
        var pairs = new List<BasePair> { Pair(0, 3), Pair(1, 4), Pair(2, 5) };

        var result = DetectPseudoknots(pairs).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(3),
                "([{)]}: all three pairs mutually cross — C(3,2)=3 pairwise crossings reported");
            Assert.That(result.All(k =>
                    k.Start1 < k.Start2 && k.Start2 < k.End1 && k.End1 < k.End2),
                Is.True, "each reported crossing satisfies i<k<j<l");
            var reported = result
                .Select(k => (k.Start1, k.End1, k.Start2, k.End2))
                .OrderBy(t => t).ToList();
            Assert.That(reported, Is.EqualTo(new List<(int, int, int, int)>
            {
                (0, 3, 1, 4), (0, 3, 2, 5), (1, 4, 2, 5),
            }), "the three reported crossings are (0,3)x(1,4), (0,3)x(2,5), (1,4)x(2,5)");
        });
    }

    // C1 — Property test (O(n^2) invariant): every reported pseudoknot's pairs satisfy i<k<j<l.
    // Evidence: INV-01 crossing condition. Deterministic seed for reproducibility.
    [Test]
    public void DetectPseudoknots_RandomPairSets_EveryResultSatisfiesCrossingCondition()
    {
        var rng = new System.Random(20260614); // fixed seed for determinism
        for (int trial = 0; trial < 200; trial++)
        {
            int n = rng.Next(2, 8);
            var pairs = new List<BasePair>();
            for (int p = 0; p < n; p++)
            {
                int x = rng.Next(0, 20);
                int y = rng.Next(0, 20);
                if (x == y) y = x + 1; // ensure two distinct positions
                pairs.Add(Pair(x, y));
            }

            foreach (var knot in DetectPseudoknots(pairs))
            {
                Assert.That(
                    knot.Start1 < knot.Start2 && knot.Start2 < knot.End1 && knot.End1 < knot.End2,
                    Is.True,
                    $"every reported pseudoknot must satisfy i<k<j<l; got " +
                    $"({knot.Start1},{knot.End1}) x ({knot.Start2},{knot.End2})");
            }
        }
    }

    #endregion

    #region Pseudoknot orders / dot-bracket-letter rendering (B12 review 2026-11-25)

    private static List<BasePair> Pairs(params (int, int)[] p) => p.Select(t => Pair(t.Item1, t.Item2)).ToList();

    // ToDotBracket = ViennaRNA 2.7.2 vrna_db_from_ptable (RNA.db_from_ptable on the same pair table).
    [TestCase(4, new[] { 0, 2, 1, 3 }, "([)]", TestName = "RnaPseudoknot001_ToDotBracket_HType_MatchesViennaRna")]
    [TestCase(6, new[] { 0, 3, 1, 4, 2, 5 }, "([{)]}", TestName = "RnaPseudoknot001_ToDotBracket_ThreeMutuallyCrossing_MatchesViennaRna")]
    [TestCase(6, new[] { 0, 4, 1, 3, 2, 5 }, "(([))]", TestName = "RnaPseudoknot001_ToDotBracket_BiotiteExample_MatchesViennaRna")]
    [TestCase(22, new[] { 0, 8, 1, 7, 4, 13, 5, 12, 10, 20, 11, 19 }, "((..[[.)).((]].....)).", TestName = "RnaPseudoknot001_ToDotBracket_KissingHairpins_MatchesViennaRna")]
    [TestCase(15, new[] { 0, 4, 1, 12, 2, 6, 3, 10, 8, 13, 9, 14 }, "([[{).].(<}.])>", TestName = "RnaPseudoknot001_ToDotBracket_FourFamilies_MatchesViennaRna")]
    [TestCase(12, new[] { 6, 0, 10, 2, 9, 3, 8, 4 }, "(.[[[.).]]].", TestName = "RnaPseudoknot001_ToDotBracket_ClosingFirstStorage_MatchesViennaRna")]
    [TestCase(5, new int[0], ".....", TestName = "RnaPseudoknot001_ToDotBracket_NoPairs_AllDots")]
    public void RnaPseudoknot001_ToDotBracket_MatchesViennaRnaDbFromPtable(int length, int[] flat, string expected)
    {
        var pairs = Enumerable.Range(0, flat.Length / 2).Select(t => (flat[2 * t], flat[2 * t + 1])).ToList();

        string db = ToDotBracket(length, pairs);

        Assert.Multiple(() =>
        {
            Assert.That(db, Is.EqualTo(expected), "RNA.db_from_ptable (ViennaRNA 2.7.2)");
            Assert.That(ParseDotBracket(db).Select(p => (p.Position1, p.Position2)).OrderBy(p => p.Position1),
                Is.EqualTo(pairs.Select(p => (Math.Min(p.Item1, p.Item2), Math.Max(p.Item1, p.Item2))).OrderBy(p => p.Item1)),
                "ParseDotBracket(ToDotBracket(P)) = P");
            Assert.That(ValidateDotBracket(db), Is.True);
        });
    }

    // 30 mutually crossing pairs use all 30 ViennaRNA families "([{<A..Z"; a 31st would need a 31st family:
    // ViennaRNA warns and drops it (RNA.db_from_ptable leaves positions 30 and 61 as '.'); we reject instead.
    [Test]
    public void RnaPseudoknot001_ToDotBracket_ThirtyFamiliesLimit()
    {
        var p30 = Enumerable.Range(0, 30).Select(i => (i, i + 30)).ToList();
        var p31 = Enumerable.Range(0, 31).Select(i => (i, i + 31)).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(ToDotBracket(60, p30),
                Is.EqualTo("([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ)]}>abcdefghijklmnopqrstuvwxyz"));
            Assert.Throws<ArgumentException>(() => ToDotBracket(62, p31));
        });
    }

    [Test]
    public void RnaPseudoknot001_ToDotBracket_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => ToDotBracket(5, new[] { (2, 2) }), "self pair");
            Assert.Throws<ArgumentException>(() => ToDotBracket(5, new[] { (-1, 3) }), "negative position");
            Assert.Throws<ArgumentException>(() => ToDotBracket(5, new[] { (0, 5) }), "position >= length");
            Assert.Throws<ArgumentException>(() => ToDotBracket(8, new[] { (0, 5), (5, 7) }), "base in two pairs");
            Assert.Throws<ArgumentException>(() => ToDotBracket(8, new[] { (0, 5), (0, 5) }), "duplicate pair");
            Assert.Throws<ArgumentNullException>(() => ToDotBracket(8, null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => ToDotBracket(-1, Array.Empty<(int, int)>()));
        });
    }

    // Orders: biotite structure.pseudoknots docstring example "ABCbac" → [0 0 1] and "(([))]";
    // other rows = exhaustive brute force (maximum nested subset per layer, ties → 5'-most pair kept).
    [TestCase(new[] { 0, 4, 1, 3, 2, 5 }, new[] { 0, 0, 1 }, "(([))]", TestName = "RnaPseudoknot001_Orders_BiotiteDocExample")]
    [TestCase(new[] { 0, 2, 1, 3 }, new[] { 0, 1 }, "([)]", TestName = "RnaPseudoknot001_Orders_HTypeEqualStems_5PrimeStemOrder0")]
    [TestCase(new[] { 0, 3, 1, 4, 2, 5 }, new[] { 0, 1, 2 }, "([{)]}", TestName = "RnaPseudoknot001_Orders_ThreeMutuallyCrossing")]
    [TestCase(new[] { 0, 8, 1, 7, 4, 13, 5, 12, 10, 20, 11, 19 }, new[] { 0, 0, 1, 1, 0, 0 }, "((..[[.)).((]].....))", TestName = "RnaPseudoknot001_Orders_KissingHairpins")]
    [TestCase(new[] { 0, 6, 2, 10, 3, 9, 4, 8 }, new[] { 1, 0, 0, 0 }, "[.(((.].)))", TestName = "RnaPseudoknot001_Orders_LargerStemKeptAtOrder0")]
    [TestCase(new[] { 0, 4, 1, 12, 2, 6, 3, 10, 8, 13, 9, 14 }, new[] { 0, 1, 1, 2, 0, 3 }, "([[{).].(<}.])>", TestName = "RnaPseudoknot001_Orders_IteratedMaximumNeedsFourLayers")]
    public void RnaPseudoknot001_AssignPseudoknotOrders_MatchesSourcedProcedure(int[] flat, int[] expectedOrders, string expectedDbl)
    {
        var pairs = Enumerable.Range(0, flat.Length / 2).Select(t => Pair(flat[2 * t], flat[2 * t + 1])).ToList();
        int length = flat.Max() + 1;

        var orders = AssignPseudoknotOrders(pairs);

        Assert.Multiple(() =>
        {
            Assert.That(orders, Is.EqualTo(expectedOrders));
            Assert.That(ToDotBracketByPseudoknotOrder(length, pairs), Is.EqualTo(expectedDbl));
        });
    }

    // Property (seeded): every order layer is nested, order 0 is a maximum nested subset (size cross-checked
    // against an independent brute force over all subsets), and the DBL string parses back to the pairs.
    [Test]
    public void RnaPseudoknot001_Orders_RandomPairSets_LayersNestedAndOrder0Maximum()
    {
        var rng = new Random(20261125);
        for (int iter = 0; iter < 300; iter++)
        {
            int p = rng.Next(0, 9);
            int length = 2 * p + rng.Next(0, 6);
            var pos = Enumerable.Range(0, length).OrderBy(_ => rng.Next()).ToList();
            var pairs = Enumerable.Range(0, p).Select(t => Pair(pos[2 * t], pos[2 * t + 1])).ToList();

            var orders = AssignPseudoknotOrders(pairs);
            for (int a = 0; a < p; a++)
                for (int b = a + 1; b < p; b++)
                    if (orders[a] == orders[b])
                        Assert.That(DetectPseudoknots(new[] { pairs[a], pairs[b] }), Is.Empty, "pairs of equal order never cross");

            int best = 0;
            for (int mask = 0; mask < (1 << p); mask++)
            {
                var subset = Enumerable.Range(0, p).Where(t => (mask & (1 << t)) != 0).Select(t => pairs[t]).ToList();
                if (subset.Count > best && !DetectPseudoknots(subset).Any()) best = subset.Count;
            }
            Assert.That(orders.Count(o => o == 0), Is.EqualTo(best), "order-0 layer is a maximum nested subset");

            string dbl = ToDotBracketByPseudoknotOrder(length, pairs);
            Assert.That(ParseDotBracket(dbl).OrderBy(x => x.Position1),
                Is.EqualTo(pairs.Select(x => (Math.Min(x.Position1, x.Position2), Math.Max(x.Position1, x.Position2))).OrderBy(x => x.Item1)));
        }
    }

    [Test]
    public void RnaPseudoknot001_Orders_InvalidInput_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AssignPseudoknotOrders(new List<BasePair>()), Is.Empty);
            Assert.Throws<ArgumentNullException>(() => AssignPseudoknotOrders(null!));
            Assert.Throws<ArgumentException>(() => AssignPseudoknotOrders(Pairs((3, 3))));
            Assert.Throws<ArgumentException>(() => AssignPseudoknotOrders(Pairs((0, 5), (5, 9))));
            Assert.Throws<ArgumentException>(() => AssignPseudoknotOrders(Pairs((-2, 5))));
            Assert.Throws<ArgumentException>(() => ToDotBracketByPseudoknotOrder(5, Pairs((0, 5))));
        });
    }

    #endregion
}
