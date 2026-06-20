using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Comparative-Genomics area — Genome Rearrangement Detection
/// (COMPGEN-REARR-001), the breakpoint-model rearrangement detector
/// <see cref="ComparativeGenomics.DetectRearrangements"/> (and its classifier
/// <see cref="ComparativeGenomics.ClassifyRearrangement"/>).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop,
/// no state corruption, no nonsense output (a negative breakpoint count, a count
/// outside the documented range, a spurious breakpoint on an identity order, a
/// missing breakpoint), and no *unhandled* runtime exception (IndexOutOfRange on a
/// single-gene or empty permutation, off-by-one in the adjacency walk, DivideByZero).
/// Every input must resolve to EITHER a well-defined, theory-correct result OR a
/// *documented, intentional* validation exception (ArgumentNullException on a null
/// argument — contract §3.3). A raw runtime exception, a hang, a wrong breakpoint
/// count, or a spurious event on an identity order is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-REARR-001 — Genome Rearrangement Detection (Breakpoints)
/// Checklist: docs/checklists/03_FUZZING.md, row 137.
/// Algorithm doc: docs/algorithms/Comparative_Genomics/Genome_Rearrangement_Detection.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: IDENTICAL ORDER (genome B == genome A's order → the
///          canonical no-change case, b(β) = 0 breakpoints — §2.4 INV-01), FULL
///          REVERSAL (genome B is the complete strand-negating reversal of A → the
///          two reversal boundaries, exactly 2 breakpoints — a single reversal
///          creates ≤ 2 breakpoints, the internal descending pairs (−(k+1),−k)
///          satisfy y = x+1 and are conserved — §2.2, §2.4 INV-02), and SINGLE
///          GENE (a length-1 permutation has NO internal adjacency → 0 breakpoints,
///          no IndexOutOfRange — §6.1 "Fewer than 2 mappable orthologs ⇒ 0 events").
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Scope relative to COMPGEN-REVERSAL-001 (row 138, CalculateReversalDistance)
/// ───────────────────────────────────────────────────────────────────────────
/// THIS row is the general breakpoint / rearrangement-boundary model: it COUNTS
/// disrupted adjacencies (DetectRearrangements) and classifies each boundary
/// (ClassifyRearrangement). The MINIMUM number of reversals (reversal *distance*,
/// Hannenhalli–Pevzner) is the separate unit COMPGEN-REVERSAL-001
/// (CalculateReversalDistance — §5.3 "Not implemented", §2.5) and is NOT exercised
/// here. The breakpoint distance is only a lower bound on the reversal distance
/// (d ≥ b/2 — §2.5), so the two units measure different quantities.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The breakpoint contract under test (Genome_Rearrangement_Detection.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Markers (genome-1 orthologs, via orthologMap) are read in genome-1 order and
/// relabelled to their genome-2 rank, signed by relative strand; the permutation
/// α is extended with sentinels 0 and n+1. A consecutive pair (x, y) of the
/// extended permutation is a BREAKPOINT iff it is not an identity adjacency —
/// equivalently iff y ≠ x + 1 (§2.2, §4.2 INV-02; the single signed test subsumes
/// both the (x,y) and (−y,−x) clauses because a reversal negates signs). The
/// breakpoint count b(α) equals the breakpoint distance d_BP = n − (common
/// adjacencies) and lies in [0, n+1] (§2.4 INV-04). b(β) = 0 for the identity
/// (§2.4 INV-01). Fewer than 2 mappable orthologs ⇒ 0 events (§6.1). Null
/// arguments ⇒ ArgumentNullException, validated eagerly before iteration (§3.3).
/// The classifier maps a sign-reversing boundary → Inversion and an
/// orientation-preserving discontinuity → Transposition (§2.4 INV-05).
///   ComparativeGenomics.DetectRearrangements(
///       IReadOnlyList&lt;Gene&gt; genome1Genes, IReadOnlyList&lt;Gene&gt; genome2Genes,
///       IReadOnlyDictionary&lt;string,string&gt; orthologMap)
///       → IEnumerable&lt;RearrangementEvent&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeRearrangementFuzzTests
{
    #region Helpers

    /// <summary>Genome-1 gene g{idx} with the given strand; coordinates spread so Start is distinct.</summary>
    private static ComparativeGenomics.Gene G1(int idx, char strand)
        => new($"g{idx}", "genome1", idx * 100, idx * 100 + 50, strand);

    /// <summary>Genome-2 reference gene h{idx} (forward strand unless overridden).</summary>
    private static ComparativeGenomics.Gene H(int idx, char strand = '+')
        => new($"h{idx}", "genome2", idx * 100, idx * 100 + 50, strand);

    /// <summary>Genome-2 reference as the identity h0..h(n-1), all '+'.</summary>
    private static List<ComparativeGenomics.Gene> Genome2Ref(int n)
        => Enumerable.Range(0, n).Select(i => H(i)).ToList();

    /// <summary>
    /// Materialise the events and assert every one is well-formed per the documented contract:
    /// the count is non-negative and within the documented maximum n+1 (§2.4 INV-04), and each
    /// event carries a non-negative Length and a parsable / classifiable signature.
    /// </summary>
    private static List<ComparativeGenomics.RearrangementEvent> AssertWellFormed(
        IEnumerable<ComparativeGenomics.RearrangementEvent> events, int markerCount)
    {
        var list = events.ToList();

        // The extended permutation of n markers has exactly n+1 internal pairs, so the breakpoint
        // count cannot exceed n+1 (INV-04). A permutation of <2 markers has no internal adjacency.
        int max = markerCount < 2 ? 0 : markerCount + 1;
        list.Count.Should().BeInRange(0, max,
            $"INV-04: breakpoint count must lie in [0, n+1] for n={markerCount} markers");

        foreach (var e in list)
        {
            e.Length.Should().BeGreaterThanOrEqualTo(0, "a breakpoint span |Δ| is non-negative");
            // Classification must never throw and must yield a defined enum value.
            var type = ComparativeGenomics.ClassifyRearrangement(e);
            type.Should().BeOneOf(
                new[]
                {
                    ComparativeGenomics.RearrangementType.Inversion,
                    ComparativeGenomics.RearrangementType.Transposition,
                },
                "every detected boundary classifies as Inversion or Transposition (INV-05)");
        }

        return list;
    }

    /// <summary>Builds an ortholog map gi → h{targetRanks[i]} for the given target permutation.</summary>
    private static Dictionary<string, string> Map(int[] targetRanks)
    {
        var map = new Dictionary<string, string>();
        for (int i = 0; i < targetRanks.Length; i++)
            map[$"g{i}"] = $"h{targetRanks[i]}";
        return map;
    }

    #endregion

    #region COMPGEN-REARR-001 — Genome Rearrangement Detection (BE: identical order, full reversal, single gene)

    #region BE — Boundary: identical order (b(β) = 0, the canonical no-change case — INV-01)

    // Identity order at every size: genome 1 reproduces genome 2's order verbatim, same strands.
    // Documented value: b(β) = 0 breakpoints (§2.4 INV-01). The canonical no-change baseline.
    [Test]
    public void DetectRearrangements_IdenticalOrder_AlwaysZeroBreakpoints()
    {
        foreach (int n in new[] { 2, 3, 5, 8, 13, 21 })
        {
            var genome2 = Genome2Ref(n);
            var genome1 = Enumerable.Range(0, n).Select(i => G1(i, '+')).ToList();
            var map = Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{i}");

            var events = AssertWellFormed(
                ComparativeGenomics.DetectRearrangements(genome1, genome2, map), n);

            events.Should().BeEmpty(
                $"INV-01: identical signed gene order has b(β)=0, no breakpoints (n={n})");
        }
    }

    // Identical order is robust to a uniform strand flip on BOTH genomes: relative strand is
    // unchanged, so the relabelled signs are all '+' and the order is still the identity ⇒ 0.
    [Test]
    public void DetectRearrangements_IdenticalOrderBothReverseStrand_StillZeroBreakpoints()
    {
        const int n = 6;
        var genome2 = Enumerable.Range(0, n).Select(i => H(i, '-')).ToList();
        var genome1 = Enumerable.Range(0, n).Select(i => G1(i, '-')).ToList();
        var map = Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{i}");

        var events = AssertWellFormed(
            ComparativeGenomics.DetectRearrangements(genome1, genome2, map), n);

        events.Should().BeEmpty(
            "relative strand is identical, so the permutation is still the identity ⇒ b(β)=0");
    }

    #endregion

    #region BE — Boundary: full reversal (the two reversal boundaries — INV-02)

    // A genuine biological full reversal of the identity (+1..+n) negates the signs of the whole
    // block, yielding (−n, −(n-1), …, −1). Extended (0, −n, …, −1, n+1): the only breakpoints are
    // the two flanks (0,−n) and (−1,n+1); every internal descending pair (−(k+1),−k) satisfies
    // y = x+1 and is conserved (§2.2, §2.4 INV-02). Documented value: EXACTLY 2 breakpoints.
    [Test]
    public void DetectRearrangements_FullReversal_ReturnsExactlyTwoBreakpoints()
    {
        foreach (int n in new[] { 2, 3, 4, 5, 8, 12 })
        {
            var genome2 = Genome2Ref(n);
            // genome-1 marker i maps to the (n-1-i)-th genome-2 gene, with flipped strand (a true
            // reversal re-orients the block) ⇒ relabelled permutation (−n, −(n-1), …, −1).
            var genome1 = Enumerable.Range(0, n).Select(i => G1(i, '-')).ToList();
            var map = new Dictionary<string, string>();
            for (int i = 0; i < n; i++) map[$"g{i}"] = $"h{n - 1 - i}";

            var events = AssertWellFormed(
                ComparativeGenomics.DetectRearrangements(genome1, genome2, map), n);

            events.Should().HaveCount(2,
                $"a full strand-negating reversal creates exactly the 2 flanking breakpoints (n={n})");
            events.Should().OnlyContain(
                e => ComparativeGenomics.ClassifyRearrangement(e)
                     == ComparativeGenomics.RearrangementType.Inversion,
                "both boundaries of a sign-negating reversal classify as Inversion (INV-05)");
        }
    }

    // An ORIENTATION-PRESERVING order reversal (descending ranks but '+' strand throughout) is the
    // worst case for the breakpoint metric: every internal pair (k+1, k) has y = k ≠ (k+1)+1, so
    // EVERY consecutive pair is a breakpoint. For n markers the extended permutation has n+1 pairs
    // ⇒ n+1 breakpoints (the INV-04 ceiling). Distinguishes the metric from a signed full reversal.
    [Test]
    public void DetectRearrangements_OrientationPreservingDescending_HitsBreakpointCeiling()
    {
        foreach (int n in new[] { 2, 3, 4, 6 })
        {
            var genome2 = Genome2Ref(n);
            var genome1 = Enumerable.Range(0, n).Select(i => G1(i, '+')).ToList();
            var map = new Dictionary<string, string>();
            for (int i = 0; i < n; i++) map[$"g{i}"] = $"h{n - 1 - i}";

            var events = AssertWellFormed(
                ComparativeGenomics.DetectRearrangements(genome1, genome2, map), n);

            events.Should().HaveCount(n + 1,
                $"a sign-stable descending order disrupts every adjacency ⇒ n+1 breakpoints (n={n})");
        }
    }

    #endregion

    #region BE — Boundary: single gene (length-1 permutation, no internal adjacency — §6.1)

    // A single mappable ortholog: no internal adjacency to evaluate ⇒ 0 events. Critically, the
    // length-1 permutation must not trigger IndexOutOfRange/DivideByZero in the adjacency walk.
    [Test]
    public void DetectRearrangements_SingleGene_ReturnsNoEventsAndDoesNotThrow()
    {
        var genome2 = Genome2Ref(3);
        var genome1 = new List<ComparativeGenomics.Gene> { G1(0, '+') };
        var map = new Dictionary<string, string> { ["g0"] = "h0" };

        Action act = () =>
        {
            var events = AssertWellFormed(
                ComparativeGenomics.DetectRearrangements(genome1, genome2, map), 1);
            events.Should().BeEmpty("a single marker has no internal adjacency ⇒ 0 breakpoints (§6.1)");
        };

        act.Should().NotThrow("a length-1 permutation must not crash the adjacency walk");
    }

    // Single gene against a single-gene reference (1×1), and against a many-gene reference: either
    // way only one anchor resolves, so there is no adjacency ⇒ 0 events, no crash.
    [Test]
    public void DetectRearrangements_SingleGeneVariousReferences_NoEventsNoCrash()
    {
        var oneToOne = ComparativeGenomics.DetectRearrangements(
            new List<ComparativeGenomics.Gene> { G1(0, '+') },
            new List<ComparativeGenomics.Gene> { H(0) },
            new Dictionary<string, string> { ["g0"] = "h0" });

        var oneToMany = ComparativeGenomics.DetectRearrangements(
            new List<ComparativeGenomics.Gene> { G1(0, '-') },
            Genome2Ref(5),
            new Dictionary<string, string> { ["g0"] = "h3" });

        AssertWellFormed(oneToOne, 1).Should().BeEmpty("1×1 has no internal adjacency ⇒ 0 events");
        AssertWellFormed(oneToMany, 1).Should().BeEmpty("a lone resolved anchor ⇒ 0 events");
    }

    // Only ONE anchor in the map resolves even though both genomes are multi-gene: the documented
    // "fewer than 2 mappable orthologs ⇒ 0 events" floor (§6.1), no off-by-one into an empty walk.
    [Test]
    public void DetectRearrangements_OnlyOneResolvableAnchor_ReturnsNoEvents()
    {
        var genome2 = Genome2Ref(4);
        var genome1 = Enumerable.Range(0, 4).Select(i => G1(i, '+')).ToList();
        // g0 resolves; g1..g3 point to non-existent genome-2 genes.
        var map = new Dictionary<string, string>
        {
            ["g0"] = "h0", ["g1"] = "hX", ["g2"] = "hY", ["g3"] = "hZ",
        };

        var events = AssertWellFormed(
            ComparativeGenomics.DetectRearrangements(genome1, genome2, map), 1);

        events.Should().BeEmpty("<2 mappable orthologs ⇒ no internal adjacency ⇒ 0 events (§6.1)");
    }

    #endregion

    #region BE — Boundary: empty genomes & null arguments (degenerate / documented exceptions — §3.3)

    // Empty genomes: 0 events, no exception (§6.1 "Empty genomes ⇒ 0 events").
    [Test]
    public void DetectRearrangements_EmptyGenomes_ReturnsNoEventsNoThrow()
    {
        Action act = () =>
        {
            var events = AssertWellFormed(
                ComparativeGenomics.DetectRearrangements(
                    new List<ComparativeGenomics.Gene>(),
                    new List<ComparativeGenomics.Gene>(),
                    new Dictionary<string, string>()),
                0);
            events.Should().BeEmpty("no markers ⇒ no breakpoints");
        };

        act.Should().NotThrow("empty genomes are a documented degenerate input, not an error");
    }

    // Each null argument throws ArgumentNullException eagerly (before iteration) — contract §3.3.
    [Test]
    public void DetectRearrangements_NullArguments_ThrowArgumentNullEagerly()
    {
        var genome = new List<ComparativeGenomics.Gene> { G1(0, '+'), G1(1, '+') };
        var map = new Dictionary<string, string>();

        // Eager: the exception fires on the call itself, not on enumeration.
        var nullG1 = () => ComparativeGenomics.DetectRearrangements(null!, Genome2Ref(2), map);
        var nullG2 = () => ComparativeGenomics.DetectRearrangements(genome, null!, map);
        var nullMap = () => ComparativeGenomics.DetectRearrangements(genome, Genome2Ref(2), null!);

        nullG1.Should().Throw<ArgumentNullException>("null genome1 is invalid (§3.3)");
        nullG2.Should().Throw<ArgumentNullException>("null genome2 is invalid (§3.3)");
        nullMap.Should().Throw<ArgumentNullException>("null orthologMap is invalid (§3.3)");
    }

    #endregion

    #region Positive sanity & metric discrimination (the documented breakpoint contract)

    // Hand-computed: the Hunter Lecture 16 worked example (−2,−3,+1,+6,−5,−4) has b(α)=6 (doc §7.1).
    // A positive sanity check that the metric counts breakpoints, not merely returns 0.
    [Test]
    public void DetectRearrangements_HunterWorkedExample_ReturnsSixBreakpoints()
    {
        var genome2 = Genome2Ref(6);
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            G1(0, '-'), // → h1 (rank 2), opposite strand ⇒ −2
            G1(1, '-'), // → h2 (rank 3) ⇒ −3
            G1(2, '+'), // → h0 (rank 1) ⇒ +1
            G1(3, '+'), // → h5 (rank 6) ⇒ +6
            G1(4, '-'), // → h4 (rank 5) ⇒ −5
            G1(5, '-'), // → h3 (rank 4) ⇒ −4
        };
        var map = new Dictionary<string, string>
        {
            ["g0"] = "h1", ["g1"] = "h2", ["g2"] = "h0",
            ["g3"] = "h5", ["g4"] = "h4", ["g5"] = "h3",
        };

        var events = AssertWellFormed(
            ComparativeGenomics.DetectRearrangements(genome1, genome2, map), 6);

        events.Should().HaveCount(6,
            "permutation (−2,−3,+1,+6,−5,−4) has exactly 6 breakpoints (Hunter Lecture 16, doc §7.1)");
    }

    // Hand-computed: one signed reversed block (+1,−4,−3,−2,+5) has exactly 2 breakpoints — the
    // two flanks of the reversed segment; the internal (−4,−3),(−3,−2) satisfy y=x+1 (doc §2.2).
    [Test]
    public void DetectRearrangements_SingleReversedBlock_ReturnsTwoBreakpoints()
    {
        var genome2 = Genome2Ref(5);
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            G1(0, '+'), // +1
            G1(1, '-'), // −4
            G1(2, '-'), // −3
            G1(3, '-'), // −2
            G1(4, '+'), // +5
        };
        var map = Map(new[] { 0, 3, 2, 1, 4 });

        var events = AssertWellFormed(
            ComparativeGenomics.DetectRearrangements(genome1, genome2, map), 5);

        events.Should().HaveCount(2,
            "one reversed block (+1,−4,−3,−2,+5) produces exactly 2 breakpoints");
    }

    // Metric discrimination: identity ⇒ 0, full signed reversal ⇒ 2, sign-stable reversal ⇒ n+1.
    // All on the SAME n so the three documented outcomes are visibly distinct (not all green-zero).
    [Test]
    public void DetectRearrangements_MetricDiscriminatesIdentityReversalAndCeiling()
    {
        const int n = 5;
        var genome2 = Genome2Ref(n);

        int identity = ComparativeGenomics.DetectRearrangements(
            Enumerable.Range(0, n).Select(i => G1(i, '+')).ToList(), genome2,
            Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{i}")).Count();

        int signedReversal = ComparativeGenomics.DetectRearrangements(
            Enumerable.Range(0, n).Select(i => G1(i, '-')).ToList(), genome2,
            Map(Enumerable.Range(0, n).Select(i => n - 1 - i).ToArray())).Count();

        int signStableReversal = ComparativeGenomics.DetectRearrangements(
            Enumerable.Range(0, n).Select(i => G1(i, '+')).ToList(), genome2,
            Map(Enumerable.Range(0, n).Select(i => n - 1 - i).ToArray())).Count();

        identity.Should().Be(0, "identity ⇒ 0 breakpoints (INV-01)");
        signedReversal.Should().Be(2, "a signed full reversal ⇒ 2 breakpoints (INV-02)");
        signStableReversal.Should().Be(n + 1, "a sign-stable descending order ⇒ n+1 breakpoints (INV-04 ceiling)");
    }

    #endregion

    #region Robustness — randomized permutations always well-formed (no hang, no crash)

    // Random signed permutations of varying size: the result must ALWAYS be well-formed
    // (count in [0, n+1], every event classifiable), regardless of the random order/strands.
    [Test]
    [CancelAfter(30000)]
    public void DetectRearrangements_RandomSignedPermutations_AlwaysWellFormed()
    {
        var rng = new Random(137); // locally-seeded, deterministic.

        for (int trial = 0; trial < 400; trial++)
        {
            int n = rng.Next(0, 16);
            var genome2 = Genome2Ref(n);

            // Random permutation of genome-2 ranks via a Fisher–Yates shuffle.
            var ranks = Enumerable.Range(0, n).ToArray();
            for (int i = n - 1; i > 0; i--)
            {
                int j = rng.Next(i + 1);
                (ranks[i], ranks[j]) = (ranks[j], ranks[i]);
            }

            var genome1 = Enumerable.Range(0, n)
                .Select(i => G1(i, rng.Next(2) == 0 ? '+' : '-')).ToList();
            var map = Map(ranks);

            var events = AssertWellFormed(
                ComparativeGenomics.DetectRearrangements(genome1, genome2, map), n);

            // INV-01 cross-check: if the random permutation happens to be the identity with all '+',
            // it must be 0 — but that is already covered by AssertWellFormed's range + the dedicated
            // identity tests; here we only assert the universal well-formedness above.
            events.Should().NotBeNull();
        }
    }

    #endregion

    #endregion
}
