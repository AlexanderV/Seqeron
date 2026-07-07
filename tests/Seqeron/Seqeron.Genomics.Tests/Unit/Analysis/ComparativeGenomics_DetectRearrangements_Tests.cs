// COMPGEN-REARR-001 — Genome Rearrangement Detection by Breakpoints
// Evidence: docs/Evidence/COMPGEN-REARR-001-Evidence.md
// TestSpec: tests/TestSpecs/COMPGEN-REARR-001.md
// Source: Bafna V, Pevzner PA (1998). Sorting by Transpositions. SIAM J. Discrete Math. 11(2):224-240.
//         Tannier E, Zheng C, Sankoff D (2009) breakpoint distance (PMC3887456).
//         Hunter College CompBio Lecture 16 (worked example b(alpha)=6; criterion y != x+1).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class ComparativeGenomics_DetectRearrangements_Tests
{
    #region Helpers

    // Genome-1 gene gi with the given strand.
    private static ComparativeGenomics.Gene G1(int idx, char strand)
        => new($"g{idx}", "genome1", idx * 100, idx * 100 + 50, strand);

    // Genome-2 reference gene hi (all forward strand unless overridden).
    private static ComparativeGenomics.Gene H(int idx, char strand = '+')
        => new($"h{idx}", "genome2", idx * 100, idx * 100 + 50, strand);

    // Build the genome-2 reference as h0..h(n-1), all '+'.
    private static List<ComparativeGenomics.Gene> Genome2Ref(int n)
        => Enumerable.Range(0, n).Select(i => H(i)).ToList();

    #endregion

    #region DetectRearrangements — MUST Tests

    // M1 — Hunter worked example: signed permutation (-2,-3,+1,+6,-5,-4) has exactly 6 breakpoints.
    // Source: Hunter Lecture 16 lines 287-289 (b(alpha)=6).
    [Test]
    public void DetectRearrangements_HunterWorkedExample_ReturnsSixBreakpoints()
    {
        // Arrange: genome2 = h0..h5 all '+'. Map gi to ranks {2,3,1,6,5,4} with signs {-,-,+,+,-,-}.
        // sign = (g1.Strand == h.Strand) ? + : -, and all h are '+', so a '-' g1 yields a negative value.
        var genome2 = Genome2Ref(6);
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            G1(0, '-'), // -> h1 (rank 2), opposite strand => -2
            G1(1, '-'), // -> h2 (rank 3) => -3
            G1(2, '+'), // -> h0 (rank 1) => +1
            G1(3, '+'), // -> h5 (rank 6) => +6
            G1(4, '-'), // -> h4 (rank 5) => -5
            G1(5, '-'), // -> h3 (rank 4) => -4
        };
        var map = new Dictionary<string, string>
        {
            ["g0"] = "h1", ["g1"] = "h2", ["g2"] = "h0",
            ["g3"] = "h5", ["g4"] = "h4", ["g5"] = "h3",
        };

        // Act
        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        // Assert: exactly 6 breakpoints, matching the published count for this permutation.
        Assert.That(events, Has.Count.EqualTo(6),
            "permutation (-2,-3,+1,+6,-5,-4) has exactly 6 breakpoints (Hunter Lecture 16)");
    }

    // M2 — Identity / collinear order yields 0 breakpoints. Source: b_beta(beta)=0 (Hunter).
    [Test]
    public void DetectRearrangements_IdenticalOrder_ReturnsNoEvents()
    {
        var genome2 = Genome2Ref(5);
        var genome1 = Enumerable.Range(0, 5).Select(i => G1(i, '+')).ToList();
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");

        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        Assert.That(events, Is.Empty, "identical signed gene order has b(beta)=0, no breakpoints");
    }

    // M3 — Single reversed block (+1,-4,-3,-2,+5) has exactly 2 breakpoints.
    // Source: reversal negates signs (Hunter line 256) + breakpoint criterion y != x+1.
    [Test]
    public void DetectRearrangements_SingleReversedBlock_ReturnsTwoBreakpoints()
    {
        // Arrange: ranks {1,4,3,2,5}; signs {+,-,-,-,+} => (+1,-4,-3,-2,+5).
        var genome2 = Genome2Ref(5);
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            G1(0, '+'), // -> h0 (rank 1) => +1
            G1(1, '-'), // -> h3 (rank 4) => -4
            G1(2, '-'), // -> h2 (rank 3) => -3
            G1(3, '-'), // -> h1 (rank 2) => -2
            G1(4, '+'), // -> h4 (rank 5) => +5
        };
        var map = new Dictionary<string, string>
        {
            ["g0"] = "h0", ["g1"] = "h3", ["g2"] = "h2", ["g3"] = "h1", ["g4"] = "h4",
        };

        // Act
        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        // Assert: only the two boundaries of the reversed block (+1,-4) and (-2,+5) are breakpoints.
        Assert.That(events, Has.Count.EqualTo(2),
            "one reversed block (+1,-4,-3,-2,+5) produces exactly 2 breakpoints");
    }

    // M4 — A sign-consecutive descending pair (-(k+1),-k) is NOT a breakpoint.
    // Source: Hunter "(-5,-4) is not a breakpoint since (4,5) appears in beta".
    [Test]
    public void DetectRearrangements_SignConsecutiveDescendingPair_NotCountedAsBreakpoint()
    {
        // Arrange: permutation (-2,-1,+3,+4,+5). Pair (-2,-1) satisfies y = x+1 (-1 = -2+1) => not a BP.
        // The only breakpoints are (0,-2) and (-1,+3). Internal (-2,-1) must be excluded.
        var genome2 = Genome2Ref(5);
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            G1(0, '-'), // -> h1 (rank 2) => -2
            G1(1, '-'), // -> h0 (rank 1) => -1
            G1(2, '+'), // -> h2 (rank 3) => +3
            G1(3, '+'), // -> h3 (rank 4) => +4
            G1(4, '+'), // -> h4 (rank 5) => +5
        };
        var map = new Dictionary<string, string>
        {
            ["g0"] = "h1", ["g1"] = "h0", ["g2"] = "h2", ["g3"] = "h3", ["g4"] = "h4",
        };

        // Act
        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        // Assert: exactly 2 breakpoints; the sign-consecutive (-2,-1) is excluded by the y=x+1 rule.
        Assert.That(events, Has.Count.EqualTo(2),
            "(-2,-1) satisfies y=x+1 and is NOT a breakpoint; only (0,-2) and (-1,+3) are");
    }

    // M5 — ClassifyRearrangement returns Inversion for a sign-flipped boundary.
    // Source: a reversal negates signs (Hunter line 256).
    [Test]
    public void ClassifyRearrangement_SignFlippedBoundary_ReturnsInversion()
    {
        // Arrange: detect on the single reversed block; every breakpoint there involves a negative
        // value (the reversed block is negative), so each classifies as an Inversion.
        var genome2 = Genome2Ref(5);
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            G1(0, '+'), G1(1, '-'), G1(2, '-'), G1(3, '-'), G1(4, '+'),
        };
        var map = new Dictionary<string, string>
        {
            ["g0"] = "h0", ["g1"] = "h3", ["g2"] = "h2", ["g3"] = "h1", ["g4"] = "h4",
        };
        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        // Act
        var types = events.Select(ComparativeGenomics.ClassifyRearrangement).ToList();

        // Assert: a reversal boundary (sign change / negative value) is classified Inversion.
        Assert.That(types, Is.All.EqualTo(ComparativeGenomics.RearrangementType.Inversion),
            "boundaries of a sign-negating reversed block classify as Inversion");
    }

    // M6 — ClassifyRearrangement returns Transposition for an orientation-preserving relocation.
    // Source: a transposition moves a block preserving orientation (Bafna & Pevzner 1998).
    [Test]
    public void ClassifyRearrangement_OrientationPreservingRelocation_ReturnsTransposition()
    {
        // Arrange: all-positive permutation (+1,+2,+4,+5,+3) — block {4,5} moved before 3,
        // no sign change. The breakpoints are positive-to-positive discontinuities.
        var genome2 = Genome2Ref(5);
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            G1(0, '+'), // +1
            G1(1, '+'), // +2
            G1(2, '+'), // +4
            G1(3, '+'), // +5
            G1(4, '+'), // +3
        };
        var map = new Dictionary<string, string>
        {
            ["g0"] = "h0", ["g1"] = "h1", ["g2"] = "h3", ["g3"] = "h4", ["g4"] = "h2",
        };
        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        // Act
        var types = events.Select(ComparativeGenomics.ClassifyRearrangement).ToList();

        // Assert: all boundaries are orientation-preserving (positive) => Transposition.
        Assert.Multiple(() =>
        {
            Assert.That(events, Is.Not.Empty, "an all-positive relocation produces breakpoints");
            Assert.That(types, Is.All.EqualTo(ComparativeGenomics.RearrangementType.Transposition),
                "orientation-preserving positive discontinuities classify as Transposition");
        });
    }

    // M7 — Fewer than 2 mappable orthologs => no events (no internal adjacency).
    [Test]
    public void DetectRearrangements_FewerThanTwoOrthologs_ReturnsNoEvents()
    {
        var genome2 = Genome2Ref(3);
        var genome1 = new List<ComparativeGenomics.Gene> { G1(0, '+'), G1(1, '+') };
        var map = new Dictionary<string, string> { ["g0"] = "h0" }; // only one anchor resolves

        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        Assert.That(events, Is.Empty, "a permutation of <2 markers has no internal adjacency, no breakpoints");
    }

    // M8 — Null genome1 => ArgumentNullException (eager validation).
    [Test]
    public void DetectRearrangements_NullGenome1_Throws()
    {
        var genome2 = Genome2Ref(3);
        var map = new Dictionary<string, string>();

        Assert.Throws<ArgumentNullException>(
            () => ComparativeGenomics.DetectRearrangements(null!, genome2, map),
            "null genome1 must throw ArgumentNullException before iteration");
    }

    // M9 — Null ortholog map => ArgumentNullException (eager validation).
    [Test]
    public void DetectRearrangements_NullOrthologMap_Throws()
    {
        var genome1 = new List<ComparativeGenomics.Gene> { G1(0, '+'), G1(1, '+') };
        var genome2 = Genome2Ref(3);

        Assert.Throws<ArgumentNullException>(
            () => ComparativeGenomics.DetectRearrangements(genome1, genome2, null!),
            "null orthologMap must throw ArgumentNullException before iteration");
    }

    // M9b — Null genome2 => ArgumentNullException (eager validation; contract §3.3: any null arg throws).
    [Test]
    public void DetectRearrangements_NullGenome2_Throws()
    {
        var genome1 = new List<ComparativeGenomics.Gene> { G1(0, '+'), G1(1, '+') };
        var map = new Dictionary<string, string>();

        Assert.Throws<ArgumentNullException>(
            () => ComparativeGenomics.DetectRearrangements(genome1, null!, map),
            "null genome2 must throw ArgumentNullException before iteration");
    }

    // M10 — ClassifyRearrangement falls back to the stored Type when TargetPosition is absent/unparsable.
    // Source: documented fallback (ComparativeGenomics.cs ClassifyRearrangement: "otherwise trust the stored Type").
    [Test]
    public void ClassifyRearrangement_NoParsableTargetPosition_ReturnsStoredType()
    {
        Assert.Multiple(() =>
        {
            // No TargetPosition => the stored Type is returned verbatim.
            var noTarget = new ComparativeGenomics.RearrangementEvent(
                Type: ComparativeGenomics.RearrangementType.Transposition,
                GenomeId: "g", Position: 10, Length: 1, TargetPosition: null);
            Assert.That(ComparativeGenomics.ClassifyRearrangement(noTarget),
                Is.EqualTo(ComparativeGenomics.RearrangementType.Transposition),
                "null TargetPosition => return the stored Type");

            // Malformed TargetPosition (no signed pair) => still the stored Type.
            var malformed = new ComparativeGenomics.RearrangementEvent(
                Type: ComparativeGenomics.RearrangementType.Inversion,
                GenomeId: "g", Position: 10, Length: 1, TargetPosition: "not-a-pair");
            Assert.That(ComparativeGenomics.ClassifyRearrangement(malformed),
                Is.EqualTo(ComparativeGenomics.RearrangementType.Inversion),
                "unparsable TargetPosition => return the stored Type");
        });
    }

    #endregion

    #region DetectRearrangements — SHOULD Tests

    // S1 — An ortholog whose target gene is absent in genome2 is skipped; remaining anchors evaluated.
    [Test]
    public void DetectRearrangements_DanglingOrtholog_SkippedAndRemainingEvaluated()
    {
        // Arrange: identity order for g0..g4 plus a dangling g5 -> hX. The dangling anchor is dropped,
        // leaving the identity permutation => 0 breakpoints.
        var genome2 = Genome2Ref(5);
        var genome1 = Enumerable.Range(0, 6).Select(i => G1(i, '+')).ToList();
        var map = Enumerable.Range(0, 5).ToDictionary(i => $"g{i}", i => $"h{i}");
        map["g5"] = "hDoesNotExist";

        // Act
        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        // Assert: dangling anchor ignored; remaining 5 anchors are identity => no breakpoints.
        Assert.That(events, Is.Empty, "dangling ortholog skipped; remaining identity order has no breakpoints");
    }

    // S2 — Event count equals the breakpoint distance d_BP = n - (common adjacencies). Source: Tannier (PMC3887456).
    [Test]
    public void DetectRearrangements_EventCount_EqualsBreakpointDistance()
    {
        // Arrange: permutation (+1,-4,-3,-2,+5) (single reversed block). n=5.
        // Common adjacencies of the extended permutation = pairs with y=x+1:
        //   (0,+1),(-4,-3),(-3,-2),(-2... no),(+5,+6) -> conserved are (0,1),(-4,-3),(-3,-2),(5,6) = 4.
        // Total internal pairs = n+1 = 6; breakpoints = 6 - 4 = 2 = d_BP.
        var genome2 = Genome2Ref(5);
        var genome1 = new List<ComparativeGenomics.Gene>
        {
            G1(0, '+'), G1(1, '-'), G1(2, '-'), G1(3, '-'), G1(4, '+'),
        };
        var map = new Dictionary<string, string>
        {
            ["g0"] = "h0", ["g1"] = "h3", ["g2"] = "h2", ["g3"] = "h1", ["g4"] = "h4",
        };

        // Act
        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        // Assert: 2 breakpoints == d_BP for this permutation (n+1 internal pairs minus 4 conserved).
        Assert.That(events, Has.Count.EqualTo(2),
            "breakpoint event count equals d_BP = (n+1 internal pairs) - (conserved adjacencies)");
    }

    // S3 — Empty genomes => no events, no exception.
    [Test]
    public void DetectRearrangements_EmptyGenomes_ReturnsNoEvents()
    {
        var genome1 = new List<ComparativeGenomics.Gene>();
        var genome2 = new List<ComparativeGenomics.Gene>();
        var map = new Dictionary<string, string>();

        var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).ToList();

        Assert.That(events, Is.Empty, "no markers => no breakpoints, no exception");
    }

    #endregion

    #region DetectRearrangements — COULD Tests (properties)

    // C1 — Property: breakpoint event count is always in [0, n+1] for n markers.
    // Source: the extended permutation of n markers has n+1 internal pairs (INV-04).
    [Test]
    public void DetectRearrangements_AnyPermutation_EventCountWithinZeroToNPlusOne()
    {
        // Deterministic fixed permutations of varying size and shape.
        var cases = new[]
        {
            new[] { 0, 1, 2, 3 },       // identity (4)
            new[] { 3, 2, 1, 0 },       // full reverse (4)
            new[] { 1, 0, 3, 2 },       // adjacent swaps (4)
            new[] { 2, 0, 1, 4, 3 },    // mixed (5)
        };

        Assert.Multiple(() =>
        {
            foreach (var targetRanks in cases)
            {
                int n = targetRanks.Length;
                var genome2 = Genome2Ref(n);
                var genome1 = Enumerable.Range(0, n).Select(i => G1(i, '+')).ToList();
                var map = new Dictionary<string, string>();
                for (int i = 0; i < n; i++) map[$"g{i}"] = $"h{targetRanks[i]}";

                int count = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).Count();

                Assert.That(count, Is.InRange(0, n + 1),
                    $"INV-04: breakpoint count must lie in [0, n+1] for n={n}");
            }
        });
    }

    // C2 — Property: identical inputs always give 0 breakpoints regardless of size (INV-01).
    [Test]
    public void DetectRearrangements_IdenticalInputs_AlwaysZeroAcrossSizes()
    {
        Assert.Multiple(() =>
        {
            foreach (int n in new[] { 2, 3, 5, 8 })
            {
                var genome2 = Genome2Ref(n);
                var genome1 = Enumerable.Range(0, n).Select(i => G1(i, '+')).ToList();
                var map = Enumerable.Range(0, n).ToDictionary(i => $"g{i}", i => $"h{i}");

                int count = ComparativeGenomics.DetectRearrangements(genome1, genome2, map).Count();

                Assert.That(count, Is.Zero, $"INV-01: identity order has 0 breakpoints for n={n}");
            }
        });
    }

    #endregion
}
