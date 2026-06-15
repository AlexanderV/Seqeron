// META-PATHWAY-001 — Metabolic Pathway Enrichment (Over-Representation Analysis)
// Evidence: docs/Evidence/META-PATHWAY-001-Evidence.md
// TestSpec: tests/TestSpecs/META-PATHWAY-001.md
// Source: Boyle EI et al. (2004). Bioinformatics 20(18):3710-3715 (GO::TermFinder);
//         PNNL Proteomics Data Analysis in R/Bioconductor §8.2 Over-Representation Analysis.

using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class MetagenomicsAnalyzer_FindPathwayEnrichment_Tests
{
    private static IReadOnlyDictionary<string, IReadOnlyCollection<string>> Db(
        params (string Id, string[] Genes)[] entries)
        => entries.ToDictionary(e => e.Id, e => (IReadOnlyCollection<string>)e.Genes);

    #region HypergeometricUpperTail — exact probabilities (M1–M7)

    // M1 — PNNL §8.2 worked example: N=8000, M=400, n=100, x=20 → P(X≥20)=7.88e-8.
    [Test]
    public void HypergeometricUpperTail_PnnlWorkedExample_MatchesPublishedPValue()
    {
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(20, 8000, 400, 100);

        Assert.That(p, Is.EqualTo(7.884747232900224e-08).Within(1e-15),
            "PNNL §8.2: P(X≥20) with N=8000, M=400, n=100 equals 7.88e-8 (published value).");
    }

    // M2 — Hypergeometric symmetry: swapping M↔n leaves P(X≥x) unchanged (same as PNNL example).
    [Test]
    public void HypergeometricUpperTail_SwappedMandN_GivesSameProbability()
    {
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(20, 8000, 100, 400);

        Assert.That(p, Is.EqualTo(7.884747232900224e-08).Within(1e-15),
            "INV-04: hypergeometric is symmetric under M↔n; (8000,100,400) equals (8000,400,100).");
    }

    // M3 — All query genes in pathway: N=10, M=5, n=5, x=5 → P(X=5)=C(5,5)C(5,0)/C(10,5)=1/252.
    [Test]
    public void HypergeometricUpperTail_AllQueryInPathway_EqualsOneOver252()
    {
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(5, 10, 5, 5);

        Assert.That(p, Is.EqualTo(1.0 / 252.0).Within(1e-10),
            "P(X≥5)=P(X=5)=C(5,5)C(5,0)/C(10,5)=1/252 by exact rational evaluation.");
    }

    // M4 — Partial overlap: N=4, M=2, n=2, x=1 → P(X≥1)=1−P(X=0)=1−C(2,0)C(2,2)/C(4,2)=5/6.
    [Test]
    public void HypergeometricUpperTail_PartialOverlap_EqualsFiveSixths()
    {
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(1, 4, 2, 2);

        Assert.That(p, Is.EqualTo(5.0 / 6.0).Within(1e-10),
            "P(X≥1)=1−C(2,0)C(2,2)/C(4,2)=1−1/6=5/6 by exact rational evaluation.");
    }

    // M5 — x=0: empty upper sum ⇒ P=1 (no over-representation possible). PNNL §8.2 corner case.
    [Test]
    public void HypergeometricUpperTail_ZeroOverlap_EqualsOne()
    {
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(0, 10, 5, 5);

        Assert.That(p, Is.EqualTo(1.0).Within(1e-12),
            "INV-02: x=0 makes the upper sum empty, so P(X≥0)=1.");
    }

    // M6 — At least one: N=10, M=5, n=5, x=1 → P(X≥1)=1−C(5,0)C(5,5)/C(10,5)=1−1/252=251/252.
    [Test]
    public void HypergeometricUpperTail_AtLeastOne_Equals251Over252()
    {
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(1, 10, 5, 5);

        Assert.That(p, Is.EqualTo(251.0 / 252.0).Within(1e-10),
            "P(X≥1)=1−P(X=0)=1−C(5,0)C(5,5)/C(10,5)=251/252.");
    }

    // M7 — Degenerate population: any of N, M, n ≤ 0 ⇒ P=1 (no success drawable). INV-03.
    [Test]
    public void HypergeometricUpperTail_DegeneratePopulation_EqualsOne()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MetagenomicsAnalyzer.HypergeometricUpperTail(2, 0, 5, 5), Is.EqualTo(1.0),
                "N=0: empty population ⇒ p=1.");
            Assert.That(MetagenomicsAnalyzer.HypergeometricUpperTail(2, 10, 0, 5), Is.EqualTo(1.0),
                "M=0: no successes in population ⇒ p=1.");
            Assert.That(MetagenomicsAnalyzer.HypergeometricUpperTail(2, 10, 5, 0), Is.EqualTo(1.0),
                "n=0: no draws ⇒ p=1.");
        });
    }

    // M7b — Infeasible overlap (x > min(M,n)): cannot draw more successes than exist ⇒ P=0.
    // N=10, M=3, n=5, x=4 → upper tail is empty (i ranges x..min(n,M)=3 < 4) ⇒ tail = 0.
    // Cross-checked with SciPy: hypergeom.sf(3, 10, 3, 5) = 0.0.
    [Test]
    public void HypergeometricUpperTail_OverlapExceedsAvailableSuccesses_EqualsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MetagenomicsAnalyzer.HypergeometricUpperTail(4, 10, 3, 5), Is.EqualTo(0.0),
                "x=4 > M=3: no feasible draw of 4 successes from 3 ⇒ P(X≥4)=0 (SciPy sf=0).");
            Assert.That(MetagenomicsAnalyzer.HypergeometricUpperTail(4, 10, 5, 3), Is.EqualTo(0.0),
                "x=4 > n=3: cannot draw 4 successes in only 3 draws ⇒ P(X≥4)=0 (symmetry, SciPy sf=0).");
        });
    }

    #endregion

    #region FindPathwayEnrichment — end-to-end, sorting, fields (M8–M14, S1, S2)

    // M8 / M10 — End-to-end with explicit background: query{g1,g2,g3}, P1={g1..g5}, N=10.
    // Overlap=3, M=5, n=3, N=10 → P(X≥3)=1/12; all PathwayEnrichment fields equal the counts used.
    [Test]
    public void FindPathwayEnrichment_ExplicitBackground_ComputesExactPValueAndFields()
    {
        var query = new[] { "g1", "g2", "g3" };
        var db = Db(("P1", new[] { "g1", "g2", "g3", "g4", "g5" }));
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10" };

        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db, background).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Pathway, Is.EqualTo("P1"), "Pathway id passed through.");
            Assert.That(result.Overlap, Is.EqualTo(3), "All 3 query genes are pathway members (x=3).");
            Assert.That(result.PathwaySize, Is.EqualTo(5), "Pathway has 5 members within the background (M=5).");
            Assert.That(result.QuerySize, Is.EqualTo(3), "Distinct query genes (n=3).");
            Assert.That(result.BackgroundSize, Is.EqualTo(10), "Background universe size (N=10).");
            Assert.That(result.PValue, Is.EqualTo(1.0 / 12.0).Within(1e-10),
                "P(X≥3) with N=10,M=5,n=3 = 1/12 by exact evaluation.");
        });
    }

    // M9 — Two pathways, one enriched and one not: results returned ascending by p-value. INV-05.
    [Test]
    public void FindPathwayEnrichment_MultiplePathways_SortedAscendingByPValue()
    {
        var query = new[] { "g1", "g2", "g3" };
        var db = Db(
            ("Enriched", new[] { "g1", "g2", "g3", "g4", "g5" }),   // overlap 3
            ("NotEnriched", new[] { "g6", "g7", "g8" }));            // overlap 0 → p=1
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10" };

        var results = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db, background);

        Assert.Multiple(() =>
        {
            Assert.That(results[0].Pathway, Is.EqualTo("Enriched"), "Most significant pathway first.");
            Assert.That(results[1].Pathway, Is.EqualTo("NotEnriched"), "Least significant pathway last.");
            Assert.That(results[0].PValue, Is.LessThanOrEqualTo(results[1].PValue),
                "INV-05: results ascending by p-value.");
            Assert.That(results[1].PValue, Is.EqualTo(1.0).Within(1e-12),
                "Disjoint pathway has overlap 0 ⇒ p=1.");
        });
    }

    // M11 — Null query → ArgumentNullException.
    [Test]
    public void FindPathwayEnrichment_NullQuery_Throws()
        => Assert.Throws<ArgumentNullException>(
            () => MetagenomicsAnalyzer.FindPathwayEnrichment(null!, Db(("P1", new[] { "g1" }))),
            "Null query genes must throw ArgumentNullException.");

    // M12 — Null pathway database → ArgumentNullException.
    [Test]
    public void FindPathwayEnrichment_NullDatabase_Throws()
        => Assert.Throws<ArgumentNullException>(
            () => MetagenomicsAnalyzer.FindPathwayEnrichment(new[] { "g1" }, null!),
            "Null pathway database must throw ArgumentNullException.");

    // M13 — Empty pathway database → empty result list.
    [Test]
    public void FindPathwayEnrichment_EmptyDatabase_ReturnsEmpty()
    {
        var results = MetagenomicsAnalyzer.FindPathwayEnrichment(
            new[] { "g1" }, Db());

        Assert.That(results, Is.Empty, "No pathways to test ⇒ empty result list.");
    }

    // M14 — Query disjoint from pathway: overlap 0, p=1. INV-02.
    [Test]
    public void FindPathwayEnrichment_NoOverlap_PValueOne()
    {
        var query = new[] { "x1", "x2" };
        var db = Db(("P1", new[] { "g1", "g2", "g3" }));
        var background = new[] { "g1", "g2", "g3", "x1", "x2" };

        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db, background).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.Overlap, Is.EqualTo(0), "No query gene is a pathway member.");
            Assert.That(result.PValue, Is.EqualTo(1.0).Within(1e-12), "INV-02: overlap 0 ⇒ p=1.");
        });
    }

    // S1 — Default background (none supplied) = union(pathway members ∪ query); N derived from that.
    // query{g1,g2,g3}, P1={g1..g5}; union = {g1..g5} ⇒ N=5, M=5, n=3, x=3 → P(X≥3)=1.0.
    [Test]
    public void FindPathwayEnrichment_DefaultBackground_UsesUnionOfMembersAndQuery()
    {
        var query = new[] { "g1", "g2", "g3" };
        var db = Db(("P1", new[] { "g1", "g2", "g3", "g4", "g5" }));

        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.BackgroundSize, Is.EqualTo(5),
                "ASSUMPTION: default background = union(pathway ∪ query) = {g1..g5}, N=5.");
            Assert.That(result.PValue, Is.EqualTo(1.0).Within(1e-12),
                "With N=M=5, n=3, x=3 the whole population is the pathway ⇒ P(X≥3)=1.");
        });
    }

    // S2 — Duplicate query genes are counted once (set semantics): n=3, not 5.
    [Test]
    public void FindPathwayEnrichment_DuplicateQueryGenes_CountedOnce()
    {
        var query = new[] { "g1", "g1", "g2", "g3", "g3" };
        var db = Db(("P1", new[] { "g1", "g2", "g3", "g4", "g5" }));
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10" };

        var result = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db, background).Single();

        Assert.Multiple(() =>
        {
            Assert.That(result.QuerySize, Is.EqualTo(3), "Distinct query genes only (n=3).");
            Assert.That(result.PValue, Is.EqualTo(1.0 / 12.0).Within(1e-10),
                "Same as M8 after de-duplication: P(X≥3) with N=10,M=5,n=3 = 1/12.");
        });
    }

    #endregion

    #region Property test (C1)

    // C1 — INV-01: p-value is always within [0,1] over a range of feasible inputs.
    [Test]
    public void HypergeometricUpperTail_VariousInputs_PValueWithinUnitInterval()
    {
        Assert.Multiple(() =>
        {
            for (int bigN = 2; bigN <= 12; bigN++)
                for (int bigM = 1; bigM < bigN; bigM++)
                    for (int n = 1; n < bigN; n++)
                        for (int x = 0; x <= Math.Min(bigM, n); x++)
                        {
                            double p = MetagenomicsAnalyzer.HypergeometricUpperTail(x, bigN, bigM, n);
                            Assert.That(p, Is.InRange(0.0, 1.0),
                                $"INV-01: P must be in [0,1] for x={x},N={bigN},M={bigM},n={n}.");
                        }
        });
    }

    #endregion
}
