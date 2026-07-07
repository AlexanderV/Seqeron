namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Metagenomics pathway over-representation unit — hypergeometric
/// pathway-enrichment (ORA) via <see cref="MetagenomicsAnalyzer.FindPathwayEnrichment"/>
/// and its core right-tail probability <see cref="MetagenomicsAnalyzer.HypergeometricUpperTail"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// no NaN/Infinity leaking into a probability, and no *unhandled* runtime
/// exception (IndexOutOfRangeException, NullReferenceException,
/// DivideByZeroException, OverflowException, …). Every input must produce EITHER
/// a well-defined, theory-correct result, OR a *documented, intentional*
/// validation exception (ArgumentNullException). A raw runtime exception, a hang,
/// a p-value outside [0, 1], or a *fabricated* enrichment on degenerate input is a
/// bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-PATHWAY-001 — pathway over-representation analysis (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 195.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation (0, -1, MaxInt, empty)
///   — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Fuzz targets (checklist row 195): "no pathway genes, all-pathway, empty".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// FindPathwayEnrichment implements hypergeometric over-representation analysis:
/// for each caller-supplied pathway it counts the pathway members within the
/// background universe (M), the distinct query genes (n), the background size (N,
/// with the query unioned into it), and the query∩pathway overlap (x), then
/// returns the right-tail probability
///   P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M, n−i) / C(N, n)
/// computed in log-space via the Lanczos ln Γ, clamped to [0, 1]. Pathways are
/// returned ascending by p-value (most significant first).
///   — docs/algorithms/Metagenomics/Pathway_Enrichment_ORA.md §2.2, §4;
///     src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs
///     (FindPathwayEnrichment lines 1025–1068, HypergeometricUpperTail 1081–1099).
///
/// Documented invariants this fixture pins (Pathway_Enrichment_ORA.md §2.4):
///   • INV-01: 0 ≤ p-value ≤ 1 (it is a probability; clamped).
///   • INV-02: x = 0 ⇒ p-value = 1 (empty upper sum: 1 − 0 = 1). — §6.1.
///   • INV-03: N, M, or n ≤ 0 ⇒ p-value = 1 (no success drawable from a
///     degenerate population). — §6.1 (Degenerate N/M/n ≤ 0).
///   • INV-04: P(X ≥ x) invariant under swapping M ↔ n (hypergeometric symmetry).
///   • INV-05: pathways returned ascending by p-value.
///
/// Boundary / malformed-input handling fixed by the doc (§3.3, §6.1) and source,
/// which these fuzz tests pin so the contract can never silently drift:
///   • EMPTY (BE):
///       – empty pathway database → empty result list, no crash. — §6.1.
///       – empty query set → every pathway has overlap x = 0 ⇒ p = 1.0 (INV-02).
///       – empty/null background → background defaults to the union of all pathway
///         members (then the query is unioned in). — §3.3, §5.2, source 1035–1045.
///   • NO PATHWAY GENES (BE): a query whose genes are members of NO pathway ⇒
///     overlap x = 0 for every pathway ⇒ p = 1.0. Over-representation is NEVER
///     fabricated when nothing is shared. — INV-02, §6.1.
///   • ALL-PATHWAY (BE): a query equal to the FULL pathway membership (complete
///     enrichment) ⇒ x = M = n, the most significant attainable p-value (< 1, and
///     never below 0). The exact small-case value is hand-derivable. — §2.2.
///   • Argument validation: null queryGenes / pathwayDatabase → ArgumentNullException
///     (source 1030–1031), checked eagerly before any work. — §3.3, §6.1.
///   • Duplicate query genes counted once (HashSet ordinal). — §6.1, §3.1.
///
/// Positive sanity (worked examples, derived INDEPENDENTLY from the hypergeometric
/// formula, NOT echoed off the implementation):
///   • §7.2 API example: query {g1,g2,g3}; pathway P1={g1..g5}; background {g1..g10};
///     ⇒ N = 10, M = 5, n = 3, x = 3 ⇒ P(X≥3) = C(5,3)C(5,0)/C(10,3) = 10/120 = 1/12.
///   • §7.1 PNNL example: HypergeometricUpperTail(20, 8000, 400, 100) = 7.884747×10⁻⁸.
///   • all-pathway: query {g1,g2,g3} = full pathway P={g1,g2,g3}; background {g1..g6}
///     ⇒ N = 6, M = 3, n = 3, x = 3 ⇒ P(X≥3) = C(3,3)C(3,0)/C(6,3) = 1/20 = 0.05.
/// A genuine full-enrichment must therefore yield a small positive p-value, so a
/// passing "no crash" result cannot be a degenerate analyzer that returns 1.0 (or
/// nothing) for everything.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism
/// ───────────────────────────────────────────────────────────────────────────
/// All inputs are hand-built or generated from a LOCALLY fixed-seed
/// `new Random(seed)` (never a shared static Rng), so every run is reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MetagenomicsPathwayFuzzTests
{
    // Independently-derived hypergeometric constants (1 − Σ tail form), NOT echoed
    // off the implementation:
    //   §7.2 API:  C(5,3)·C(5,0)/C(10,3) = 10/120 = 1/12.
    private const double ApiExamplePValue = 1.0 / 12.0;
    //   §7.1 PNNL: HypergeometricUpperTail(20, 8000, 400, 100) = 7.884747×10⁻⁸.
    private const double PnnlPValue = 7.884747217109681e-08;
    //   all-pathway: C(3,3)·C(3,0)/C(6,3) = 1/20 = 0.05.
    private const double AllPathwayPValue = 1.0 / 20.0;

    private static Dictionary<string, IReadOnlyCollection<string>> BuildPathwayDatabase() =>
        new()
        {
            ["P1"] = new[] { "g1", "g2", "g3", "g4", "g5" },
        };

    #region META-PATHWAY-001 — hypergeometric pathway over-representation

    // ════════════════════════════════════════════════════════════════════════
    //  Positive sanity — worked examples must be reproduced EXACTLY.
    //  Guards against a degenerate analyzer (constant 1.0, or empty output) that
    //  would pass every boundary test below.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void HypergeometricUpperTail_ApiWorkedExample_EqualsOneTwelfth()
    {
        // §7.2: P(X≥3 | N=10, M=5, n=3) = C(5,3)·C(5,0)/C(10,3) = 10/120 = 1/12.
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(x: 3, bigN: 10, bigM: 5, n: 3);

        p.Should().BeApproximately(ApiExamplePValue, 1e-12,
            "single feasible term C(5,3)·C(5,0)/C(10,3) = 1/12 — §7.2 worked example");
    }

    [Test]
    public void HypergeometricUpperTail_PnnlWorkedExample_MatchesPublishedValue()
    {
        // §7.1: published right-tail value ≈ 7.88×10⁻⁸ for the large finite universe.
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(x: 20, bigN: 8000, bigM: 400, n: 100);

        p.Should().BeApproximately(PnnlPValue, 1e-13,
            "matches the published PNNL §8.2 right-tail value 7.884747×10⁻⁸ — §7.1");
        p.Should().BeGreaterThan(0.0).And.BeLessThan(1.0)
            .And.NotBe(double.NaN);
    }

    [Test]
    public void FindPathwayEnrichment_ApiWorkedExample_ReportsOneTwelfthAndCounts()
    {
        var query = new[] { "g1", "g2", "g3" };
        var pathways = BuildPathwayDatabase();
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10" };

        var results = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, background);

        var r = results.Should().ContainSingle().Subject;
        r.Pathway.Should().Be("P1");
        r.BackgroundSize.Should().Be(10, "background {g1..g10} unioned with query ⊆ it ⇒ N = 10");
        r.PathwaySize.Should().Be(5, "all 5 P1 members lie in the background ⇒ M = 5");
        r.QuerySize.Should().Be(3, "distinct query genes ⇒ n = 3");
        r.Overlap.Should().Be(3, "g1,g2,g3 are all in P1 ⇒ x = 3");
        r.PValue.Should().BeApproximately(ApiExamplePValue, 1e-12,
            "P(X≥3 | N=10,M=5,n=3) = 1/12 — §7.2 API example, derived from the formula not the code");
    }

    // ───────────────────────────────────────────────────────────────────────
    // INV-04: hypergeometric symmetry — P(X ≥ x) is invariant under M ↔ n.
    // Pinned with two independently-equal calls (not echoed off the code).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void HypergeometricUpperTail_IsSymmetricUnderSwapOfMAndN()
    {
        double a = MetagenomicsAnalyzer.HypergeometricUpperTail(x: 2, bigN: 10, bigM: 4, n: 6);
        double b = MetagenomicsAnalyzer.HypergeometricUpperTail(x: 2, bigN: 10, bigM: 6, n: 4);

        a.Should().BeApproximately(b, 1e-12, "P(X≥x) is symmetric under M ↔ n — INV-04");
        a.Should().BeInRange(0.0, 1.0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // INV-05: pathways ranked ascending by p-value (most significant first).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindPathwayEnrichment_ResultsAreSortedAscendingByPValue()
    {
        var query = new[] { "g1", "g2", "g3" };
        var pathways = new Dictionary<string, IReadOnlyCollection<string>>
        {
            ["enriched"] = new[] { "g1", "g2", "g3" },        // x = 3 (full overlap) ⇒ small p
            ["partial"]  = new[] { "g1", "g4", "g5" },        // x = 1 ⇒ larger p
            ["disjoint"] = new[] { "g6", "g7", "g8" },        // x = 0 ⇒ p = 1.0
        };
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10" };

        var results = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, background);

        results.Select(r => r.PValue).Should().BeInAscendingOrder("INV-05: ranked by significance");
        results[0].Pathway.Should().Be("enriched", "the fully-overlapping pathway is most significant");
        results.Should().OnlyContain(r => r.PValue >= 0.0 && r.PValue <= 1.0, "INV-01");
    }

    #endregion

    #region META-PATHWAY-001 — BE boundary: ALL-PATHWAY (full enrichment)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: ALL-PATHWAY (BE). The query equals the FULL pathway membership
    // (complete enrichment / 100% pathway completeness): x = M = n. The result
    // must be the most significant attainable p-value — small, strictly positive,
    // and ≤ 1 (never below 0). Exact hand-derived value for the small case.
    //   query {g1,g2,g3} = pathway P {g1,g2,g3}; background {g1..g6}
    //   ⇒ N=6, M=3, n=3, x=3 ⇒ P(X≥3) = C(3,3)·C(3,0)/C(6,3) = 1/20 = 0.05.
    // — Pathway_Enrichment_ORA.md §2.2.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindPathwayEnrichment_AllPathwayGenesInQuery_GivesExactSmallPValue()
    {
        var query = new[] { "g1", "g2", "g3" };
        var pathways = new Dictionary<string, IReadOnlyCollection<string>>
        {
            ["full"] = new[] { "g1", "g2", "g3" },
        };
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "g6" };

        var r = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, background)
            .Should().ContainSingle().Subject;

        r.Overlap.Should().Be(3, "the whole pathway is in the query ⇒ x = M = 3");
        r.PathwaySize.Should().Be(3);
        r.QuerySize.Should().Be(3);
        r.BackgroundSize.Should().Be(6);
        r.PValue.Should().BeApproximately(AllPathwayPValue, 1e-12,
            "complete enrichment P(X≥3 | N=6,M=3,n=3) = C(3,3)/C(6,3) = 1/20 — §2.2, derived from the formula");
        r.PValue.Should().BeGreaterThan(0.0).And.BeLessThan(1.0).And.NotBe(double.NaN);
    }

    // Direct full-enrichment via the core function: x = M = n with a small surplus
    // background. p = C(M,M)·C(N−M,0)/C(N,n) = 1/C(N,n), independently derivable.
    [Test]
    public void HypergeometricUpperTail_FullEnrichment_EqualsReciprocalOfChooseNn()
    {
        // N=6, M=3, n=3, x=3 ⇒ p = 1 / C(6,3) = 1/20.
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(x: 3, bigN: 6, bigM: 3, n: 3);

        p.Should().BeApproximately(1.0 / 20.0, 1e-12, "full enrichment ⇒ p = 1/C(N,n) = 1/20");
    }

    #endregion

    #region META-PATHWAY-001 — BE boundary: NO PATHWAY GENES (zero overlap)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: NO PATHWAY GENES (BE). A query whose genes are members of NO
    // pathway ⇒ overlap x = 0 for every pathway ⇒ p-value = 1.0. Enrichment is
    // never fabricated when nothing is shared. — INV-02, §6.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindPathwayEnrichment_QuerySharesNoPathwayGenes_PValueIsOne()
    {
        var query = new[] { "q1", "q2", "q3" };           // none of these are pathway members
        var pathways = BuildPathwayDatabase();             // P1 = g1..g5
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "q1", "q2", "q3" };

        var r = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, background)
            .Should().ContainSingle().Subject;

        r.Overlap.Should().Be(0, "no query gene is a pathway member ⇒ x = 0");
        r.PValue.Should().Be(1.0, "x = 0 ⇒ empty upper sum ⇒ p = 1.0 (INV-02, §6.1)");
    }

    // Core function: x = 0 must return exactly 1.0 regardless of N, M, n (INV-02).
    [Test]
    public void HypergeometricUpperTail_ZeroOverlap_ReturnsExactlyOne()
    {
        MetagenomicsAnalyzer.HypergeometricUpperTail(x: 0, bigN: 10, bigM: 5, n: 3)
            .Should().Be(1.0, "x = 0 ⇒ empty upper sum: 1 − 0 = 1 (INV-02)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // BE: degenerate population N, M, or n ≤ 0 (including the -1 boundary) ⇒
    // p-value = 1.0 — no success drawable from a degenerate population (INV-03).
    // Also exercises the int.MaxValue boundary on x (overlap larger than any
    // feasible population) which yields an empty feasible sum ⇒ clamps to ≤ 1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void HypergeometricUpperTail_DegenerateOrBoundaryArguments_StayInUnitInterval()
    {
        // INV-03 degenerate cases — all p = 1.0.
        MetagenomicsAnalyzer.HypergeometricUpperTail(x: 1, bigN: 0, bigM: 5, n: 3).Should().Be(1.0, "N = 0 (INV-03)");
        MetagenomicsAnalyzer.HypergeometricUpperTail(x: 1, bigN: 10, bigM: 0, n: 3).Should().Be(1.0, "M = 0 (INV-03)");
        MetagenomicsAnalyzer.HypergeometricUpperTail(x: 1, bigN: 10, bigM: 5, n: 0).Should().Be(1.0, "n = 0 (INV-03)");

        // -1 boundary (BE) — negative population/sample is degenerate ⇒ 1.0.
        MetagenomicsAnalyzer.HypergeometricUpperTail(x: 1, bigN: -1, bigM: 5, n: 3).Should().Be(1.0, "N = -1 (BE)");
        MetagenomicsAnalyzer.HypergeometricUpperTail(x: -1, bigN: 10, bigM: 5, n: 3).Should().Be(1.0, "x = -1 ⇒ x ≤ 0 (BE)");

        // int.MaxValue overlap (BE): infeasible (x > min(n,M)) ⇒ empty sum ⇒ p = 0,
        // still inside [0, 1], no overflow, no NaN.
        double p = MetagenomicsAnalyzer.HypergeometricUpperTail(x: int.MaxValue, bigN: 10, bigM: 5, n: 3);
        p.Should().BeInRange(0.0, 1.0).And.NotBe(double.NaN, "x ≫ min(n,M) ⇒ empty feasible sum, p clamped (INV-01)");
    }

    #endregion

    #region META-PATHWAY-001 — BE boundary: EMPTY

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY PATHWAY DATABASE (BE). Nothing to test ⇒ empty result
    // list, no crash, no fabricated row. — §6.1 (Empty pathway database).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindPathwayEnrichment_EmptyPathwayDatabase_YieldsEmptyList()
    {
        var query = new[] { "g1", "g2", "g3" };
        var emptyDb = new Dictionary<string, IReadOnlyCollection<string>>();

        var results = MetagenomicsAnalyzer.FindPathwayEnrichment(query, emptyDb);

        results.Should().BeEmpty("an empty pathway database has nothing to test (§6.1)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY QUERY SET (BE). No interesting genes ⇒ every pathway
    // has overlap x = 0 ⇒ p = 1.0 (INV-02). No crash on the empty enumeration.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindPathwayEnrichment_EmptyQuery_AllPValuesAreOne()
    {
        var emptyQuery = Array.Empty<string>();
        var pathways = BuildPathwayDatabase();
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "g6" };

        var results = MetagenomicsAnalyzer.FindPathwayEnrichment(emptyQuery, pathways, background);

        results.Should().ContainSingle();
        results[0].QuerySize.Should().Be(0, "no query genes ⇒ n = 0");
        results[0].Overlap.Should().Be(0, "no query genes ⇒ x = 0");
        results[0].PValue.Should().Be(1.0, "x = 0 (also n = 0) ⇒ p = 1.0 (INV-02/INV-03)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY / NULL BACKGROUND (BE). When no background is supplied
    // (null or empty), it defaults to the union of all pathway members, then the
    // query is unioned in. — §3.3, §5.2, source 1035–1045.
    //   query {g1,g2,g3}; P1={g1..g5}; no background ⇒ background = {g1..g5}
    //   ⇒ N = 5, M = 5, n = 3, x = 3 ⇒ P(X≥3) = C(5,3)·C(0,0)/C(5,3) = 1.0.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindPathwayEnrichment_NullBackground_DefaultsToUnionOfPathwayMembers()
    {
        var query = new[] { "g1", "g2", "g3" };
        var pathways = BuildPathwayDatabase(); // P1 = g1..g5

        var r = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, backgroundGenes: null)
            .Should().ContainSingle().Subject;

        r.BackgroundSize.Should().Be(5, "default background = union of pathway members {g1..g5} ⇒ N = 5 (§5.2)");
        r.PathwaySize.Should().Be(5);
        r.Overlap.Should().Be(3);
        // N = M = 5, n = 3, x = 3 ⇒ only i = 3 term: C(5,3)·C(0,0)/C(5,3) = 1.0.
        r.PValue.Should().BeApproximately(1.0, 1e-12,
            "when background = pathway, every query member is a success ⇒ p = 1.0");
    }

    // Empty (non-null) background behaves identically to null: it also defaults
    // to the pathway-member union. — source line 1038 (background.Count == 0).
    [Test]
    public void FindPathwayEnrichment_EmptyBackground_DefaultsToUnionOfPathwayMembers()
    {
        var query = new[] { "g1", "g2", "g3" };
        var pathways = BuildPathwayDatabase();

        var r = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, backgroundGenes: Array.Empty<string>())
            .Should().ContainSingle().Subject;

        r.BackgroundSize.Should().Be(5, "an empty background also defaults to the pathway-member union (line 1038)");
    }

    #endregion

    #region META-PATHWAY-001 — duplicate query genes (set semantics) & arg validation

    // ───────────────────────────────────────────────────────────────────────
    // Duplicate query genes are counted once (HashSet, ordinal). — §6.1, §3.1.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindPathwayEnrichment_DuplicateQueryGenes_CountedOnce()
    {
        var query = new[] { "g1", "g1", "g1", "g2" }; // distinct = {g1, g2}
        var pathways = BuildPathwayDatabase();
        var background = new[] { "g1", "g2", "g3", "g4", "g5", "g6" };

        var r = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, background)
            .Should().ContainSingle().Subject;

        r.QuerySize.Should().Be(2, "duplicates collapse: distinct {g1,g2} ⇒ n = 2 (§6.1, set semantics)");
        r.Overlap.Should().Be(2, "both distinct query genes are P1 members ⇒ x = 2");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Argument validation (BE): null queryGenes / pathwayDatabase →
    // ArgumentNullException, checked eagerly. — §3.3, §6.1, source 1030–1031.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindPathwayEnrichment_NullArguments_ThrowDocumentedExceptions()
    {
        var query = new[] { "g1" };
        var pathways = BuildPathwayDatabase();

        Action nullQuery = () => MetagenomicsAnalyzer.FindPathwayEnrichment(null!, pathways);
        nullQuery.Should().Throw<ArgumentNullException>().WithParameterName("queryGenes");

        Action nullDb = () => MetagenomicsAnalyzer.FindPathwayEnrichment(query, null!);
        nullDb.Should().Throw<ArgumentNullException>().WithParameterName("pathwayDatabase");
    }

    #endregion

    #region META-PATHWAY-001 — randomized boundary sweep

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random boundary batch (BE) under a time budget.
    // A deterministic, locally-seeded generator builds many pathway databases and
    // queries spanning the boundaries (empty queries, disjoint genes, full overlap,
    // varying universe sizes). FindPathwayEnrichment must process every case
    // without crashing or hanging, and EVERY result must be well-formed:
    //   • p-value ∈ [0, 1], finite (no NaN/Infinity) — INV-01;
    //   • x = 0 ⇒ p = 1.0 — INV-02;
    //   • results sorted ascending by p-value — INV-05;
    //   • counts consistent: 0 ≤ Overlap ≤ min(PathwaySize, QuerySize),
    //     PathwaySize ≤ BackgroundSize, QuerySize ≤ BackgroundSize.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void FindPathwayEnrichment_RandomBoundaryBatch_NeverCrashesAndStaysWellFormed()
    {
        var rng = new Random(20260620); // locally fixed seed — deterministic

        for (int trial = 0; trial < 300; trial++)
        {
            int universe = rng.Next(0, 25); // includes 0 (empty universe)
            var genePool = Enumerable.Range(0, Math.Max(universe, 1)).Select(i => $"g{i}").ToArray();

            // Random query: sometimes empty, sometimes disjoint, sometimes overlapping.
            int qSize = rng.Next(0, universe + 1);
            var query = Enumerable.Range(0, qSize)
                .Select(_ => genePool[rng.Next(genePool.Length)])
                .ToArray();

            // 0..3 pathways, each a random subset (possibly empty / possibly the whole pool).
            int nPathways = rng.Next(0, 4);
            var pathways = new Dictionary<string, IReadOnlyCollection<string>>();
            for (int p = 0; p < nPathways; p++)
            {
                int mSize = rng.Next(0, universe + 1);
                var members = Enumerable.Range(0, mSize)
                    .Select(_ => genePool[rng.Next(genePool.Length)])
                    .ToArray();
                pathways[$"P{p}"] = members;
            }

            // Background: null (default), empty, or an explicit random universe.
            string[]? background = rng.Next(3) switch
            {
                0 => null,
                1 => Array.Empty<string>(),
                _ => genePool.Take(rng.Next(0, genePool.Length + 1)).ToArray(),
            };

            var results = MetagenomicsAnalyzer.FindPathwayEnrichment(query, pathways, background);

            results.Count.Should().Be(pathways.Count, "exactly one result per pathway");
            results.Select(r => r.PValue).Should().BeInAscendingOrder("INV-05");

            foreach (var r in results)
            {
                r.PValue.Should().BeInRange(0.0, 1.0, "p-value is a probability (INV-01)");
                double.IsNaN(r.PValue).Should().BeFalse("no NaN p-value on boundary input");
                double.IsInfinity(r.PValue).Should().BeFalse("no Infinity p-value on boundary input");

                r.Overlap.Should().BeGreaterThanOrEqualTo(0);
                r.Overlap.Should().BeLessThanOrEqualTo(Math.Min(r.PathwaySize, r.QuerySize),
                    "overlap cannot exceed pathway size or query size");
                r.PathwaySize.Should().BeLessThanOrEqualTo(r.BackgroundSize,
                    "pathway members are intersected with the background");
                r.QuerySize.Should().BeLessThanOrEqualTo(r.BackgroundSize,
                    "the query is unioned into the background");

                if (r.Overlap == 0)
                    r.PValue.Should().Be(1.0, "x = 0 ⇒ p = 1.0 (INV-02)");
            }
        }
    }

    #endregion
}
