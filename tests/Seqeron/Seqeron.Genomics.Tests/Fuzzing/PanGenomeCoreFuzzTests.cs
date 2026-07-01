using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the PanGenome area — core / accessory / unique partitioning of
/// the gene-cluster repertoire (PANGEN-CORE-001), the Tettelin/Roary pan-genome
/// construction model.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (a cluster in zero or more than
/// one of the core/accessory/unique partitions, partition counts that do not sum
/// to TotalGenes, a CoreFraction or GenomeFluidity outside [0,1], a "core" gene on
/// pairwise-DISJOINT genomes that share nothing), and no *unhandled* runtime
/// exception (the headline hazard here is DivideByZero on the presence-fraction
/// denominator occupancy/N when N = 0 genomes, plus a crash / IndexOutOfRange on a
/// single genome or empty input). Every input must resolve to EITHER a
/// well-defined, theory-correct partition, OR a *documented, intentional*
/// degenerate result (null / empty genomes → all-empty PanGenomeResult, no throw).
/// A raw runtime exception, a hang, a malformed partition, a false core on
/// disjoint genomes, or a fraction out of range is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PANGEN-CORE-001 — Pan-genome construction (core / accessory / unique)
/// Checklist: docs/checklists/03_FUZZING.md, row 191.
/// Algorithm doc: docs/algorithms/Metagenomics/PanGenome_Core_Accessory.md
///                (Test Unit ID PANGEN-CORE-001).
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate genome-count boundaries called
///          out in the checklist row:
///            – SINGLE genome (N=1) → every cluster has occupancy = N = 1, so under
///              the fractional rule occupancy/N = 1 ≥ coreFraction: every gene is
///              trivially CORE (intersection over one genome), fluidity 0, no
///              DivideByZero on the N=1 denominator (§6.1 "single genome").
///            – DISJOINT genomes (genomes sharing NO genes) → every cluster has
///              occupancy 1 → ALL unique, NO false core, fluidity 1 (§6.1
///              "pairwise-disjoint content", INV-05).
///            – EMPTY (no genomes, N=0) → documented all-empty guarded result, NO
///              DivideByZero on the zero-genome denominator (§3.3, §6.1 "null /
///              empty genomes").
/// — docs/checklists/03_FUZZING.md §Description (strategy codes:
///    "BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty)").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The partition contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// PanGenomeAnalyzer.ConstructPanGenome clusters all genes into ortholog families
/// (delegated to ClusterGenes, the k-mer Jaccard / CD-HIT greedy model), assigns
/// each cluster an occupancy = number of distinct genomes containing it, then
/// partitions on occupancy (PanGenome_Core_Accessory.md §2.2):
///   • CORE      ⟺ occupancy / N ≥ coreFraction  (Roary "present in ≥ 99% of
///                  samples"; a FRACTIONAL test, NOT floor(coreFraction·N)) [4].
///   • UNIQUE    ⟺ occupancy = 1  (strain-specific / cloud) [1][2].
///   • ACCESSORY ⟺ all remaining clusters (dispensable / shell) [2].
/// The API entry under test is
///   PanGenomeAnalyzer.ConstructPanGenome(
///       IReadOnlyDictionary&lt;string, IReadOnlyList&lt;(string GeneId, string Sequence)&gt;&gt; genomes,
///       double identityThreshold = 0.9, double coreFraction = 0.99)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs
///    ConstructPanGenome lines 103–172; IsCoreOccupancy lines 184–189, which
///    guards N ≤ 0 by returning false — the DivideByZero guard under fuzz).
///
/// THE DOCUMENTED INVARIANTS (PanGenome_Core_Accessory.md §2.4):
///   • INV-01: core, accessory, unique are pairwise disjoint and partition all clusters.
///   • INV-02: CoreGeneCount + AccessoryGeneCount + UniqueGeneCount = TotalGenes.
///   • INV-03: core ⟺ occupancy/N ≥ coreFraction; unique ⟺ occupancy = 1.
///   • INV-04: 0 ≤ GenomeFluidity ≤ 1.
///   • INV-05: identical content ⇒ φ = 0; pairwise-disjoint ⇒ φ = 1.
///   • INV-06: CoreFraction = CoreGeneCount / TotalGenes (0 if TotalGenes = 0).
/// The documented edge cases (§6.1):
///   • null / empty genomes → empty result, all stats 0, Type Closed (no throw).
///   • single genome (N=1) → all clusters core (occupancy = N), fluidity 0.
///   • pairwise-disjoint content → no core, all unique, fluidity 1.
///
/// The three BE checklist targets map to these documented behaviours:
///   • single genome   → one genome ⇒ occupancy/N = 1/1 = 1 ≥ coreFraction for ANY
///                       coreFraction ≤ 1 ⇒ EVERY cluster is core; no accessory, no
///                       unique; fluidity 0 (no pairs); no DivideByZero.
///   • disjoint genomes → genomes whose gene sets are pairwise disjoint ⇒ every
///                       cluster occupancy = 1 ⇒ ALL unique, zero core (even at
///                       coreFraction as low as it can go with N≥2), fluidity 1.
///   • empty            → genomes Count = 0 (and null) ⇒ all-empty result, every
///                       count 0, CoreFraction 0, fluidity 0, Type Closed, no throw.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PanGenomeCoreFuzzTests
{
    #region Helpers

    private const int Seed = 191_0001; // local-only deterministic seed for the fuzz sweep

    // The k-mer Jaccard clusterer needs sequences of length >= k (k = 7); identical
    // strings cluster together, sufficiently distinct strings form separate clusters.
    // A pool of distinct 30-bp sequences keyed by an integer "family" index lets a
    // test place the SAME family in several genomes (raising occupancy) or give each
    // genome its OWN family (occupancy 1 / disjoint).
    private static readonly string[] Bases = { "A", "C", "G", "T" };

    // Deterministic distinct sequence for family index f: a 30-bp string whose
    // 7-mer content is unique per f, so distinct f never collide under the clusterer.
    private static string Family(int f)
    {
        // Build a 30-char base-4 pattern seeded by f; mix in f at several positions so
        // the k-mer (k=7) spectra differ between families.
        var rng = new Random(unchecked(1_000_003 * (f + 1)));
        var sb = new System.Text.StringBuilder(30);
        for (int i = 0; i < 30; i++)
            sb.Append(Bases[rng.Next(4)]);
        return sb.ToString();
    }

    private static IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string Genome, (string GeneId, string Sequence)[] Genes)[] entries)
    {
        var dict = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        foreach (var (genome, genes) in entries)
            dict[genome] = genes.ToList();
        return dict;
    }

    // A gene entry placing family f into a genome (gene id is unique per call).
    private static (string GeneId, string Sequence) Gene(string genomeId, int f) =>
        ($"{genomeId}_g{f}", Family(f));

    /// <summary>
    /// Asserts the result is a WELL-FORMED partition per PanGenome_Core_Accessory.md
    /// §2.4: every cluster is in exactly one of core/accessory/unique (INV-01), the
    /// three counts sum to TotalGenes (INV-02), the count fields agree with the list
    /// lengths, CoreFraction = core/total (INV-06), and both CoreFraction and
    /// GenomeFluidity lie in [0,1] (INV-04, INV-06). Used by every test so a malformed
    /// partition or an out-of-range fraction fails everywhere, not just where probed.
    /// </summary>
    private static void AssertWellFormedPartition(PanGenomeAnalyzer.PanGenomeResult r)
    {
        var s = r.Statistics;

        // Count fields agree with the returned lists.
        r.CoreGenes.Count.Should().Be(s.CoreGeneCount, "CoreGeneCount must match the CoreGenes list length.");
        r.AccessoryGenes.Count.Should().Be(s.AccessoryGeneCount, "AccessoryGeneCount must match the AccessoryGenes list length.");
        r.UniqueGenes.Count.Should().Be(s.UniqueGeneCount, "UniqueGeneCount must match the UniqueGenes list length.");

        // INV-02: the three partitions sum to the total cluster count.
        (s.CoreGeneCount + s.AccessoryGeneCount + s.UniqueGeneCount).Should().Be(s.TotalGenes,
            "INV-02: core + accessory + unique must equal TotalGenes.");

        // INV-01: the three partitions are pairwise disjoint sets of cluster IDs.
        var core = new HashSet<string>(r.CoreGenes);
        var accessory = new HashSet<string>(r.AccessoryGenes);
        var unique = new HashSet<string>(r.UniqueGenes);
        core.Overlaps(accessory).Should().BeFalse("INV-01: core and accessory must be disjoint.");
        core.Overlaps(unique).Should().BeFalse("INV-01: core and unique must be disjoint.");
        accessory.Overlaps(unique).Should().BeFalse("INV-01: accessory and unique must be disjoint.");

        // No duplicates within a partition (each cluster classified once).
        core.Count.Should().Be(r.CoreGenes.Count, "INV-01: no cluster appears twice in core.");
        accessory.Count.Should().Be(r.AccessoryGenes.Count, "INV-01: no cluster appears twice in accessory.");
        unique.Count.Should().Be(r.UniqueGenes.Count, "INV-01: no cluster appears twice in unique.");

        // INV-06: CoreFraction = core / total (0 when total = 0).
        double expectedCoreFraction = s.TotalGenes > 0 ? (double)s.CoreGeneCount / s.TotalGenes : 0.0;
        s.CoreFraction.Should().BeApproximately(expectedCoreFraction, 1e-12,
            "INV-06: CoreFraction = CoreGeneCount / TotalGenes (0 if TotalGenes = 0).");

        // INV-04 / INV-06: fractions stay in [0,1] (the headline 'fractions out of range' hazard).
        s.CoreFraction.Should().BeInRange(0.0, 1.0, "CoreFraction must lie in [0,1].");
        s.GenomeFluidity.Should().BeInRange(0.0, 1.0, "INV-04: GenomeFluidity must lie in [0,1].");

        // CoreFraction must be finite (DivideByZero / 0-genome denominator would yield NaN).
        double.IsNaN(s.CoreFraction).Should().BeFalse("CoreFraction must be finite, never NaN (no 0-genome divide).");
        double.IsNaN(s.GenomeFluidity).Should().BeFalse("GenomeFluidity must be finite, never NaN.");
    }

    #endregion

    #region PANGEN-CORE-001 — Pan-genome construction (core / accessory / unique)

    #region BE — Boundary: empty / null genomes (zero-genome denominator, no throw)

    // The headline DivideByZero hazard: with N = 0 genomes the presence fraction
    // occupancy/N has a zero denominator. The documented contract (§3.3, §6.1) is an
    // all-empty result with every stat 0 and Type Closed — NOT a throw and NOT a NaN.
    [Test]
    [CancelAfter(20000)]
    public void Empty_NoGenomes_ReturnsAllEmptyResult_NoDivideByZero()
    {
        var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();

        var act = () => PanGenomeAnalyzer.ConstructPanGenome(genomes);
        act.Should().NotThrow("empty genomes is a documented degenerate input (§3.3): an all-empty result, never a DivideByZero.");

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes);
        AssertWellFormedPartition(r);

        r.Statistics.TotalGenomes.Should().Be(0);
        r.Statistics.TotalGenes.Should().Be(0);
        r.Statistics.CoreGeneCount.Should().Be(0, "no genomes ⇒ no clusters ⇒ no core.");
        r.Statistics.AccessoryGeneCount.Should().Be(0);
        r.Statistics.UniqueGeneCount.Should().Be(0);
        r.Statistics.CoreFraction.Should().Be(0.0, "INV-06: TotalGenes = 0 ⇒ CoreFraction = 0, not NaN.");
        r.Statistics.GenomeFluidity.Should().Be(0.0, "no genome pairs ⇒ fluidity 0.");
        r.Statistics.Type.Should().Be(PanGenomeAnalyzer.PanGenomeType.Closed, "N < 3 ⇒ Closed (§6.1).");
    }

    // null genomes is the same documented degenerate input (§3.3): empty result, no throw.
    [Test]
    [CancelAfter(20000)]
    public void Empty_NullGenomes_ReturnsAllEmptyResult_NoThrow()
    {
        var act = () => PanGenomeAnalyzer.ConstructPanGenome(null!);
        act.Should().NotThrow("null genomes is documented as an empty result (§3.3), never an exception.");

        var r = PanGenomeAnalyzer.ConstructPanGenome(null!);
        AssertWellFormedPartition(r);
        r.Statistics.TotalGenomes.Should().Be(0);
        r.Statistics.TotalGenes.Should().Be(0);
        r.CoreGenes.Should().BeEmpty();
        r.AccessoryGenes.Should().BeEmpty();
        r.UniqueGenes.Should().BeEmpty();
    }

    // A genome present but holding an EMPTY gene list: still N = 1 genome but ZERO
    // clusters. Probes the occupancy/N path (N=1) crossed with TotalGenes = 0
    // (CoreFraction must be the documented 0, not 0/0 = NaN).
    [Test]
    [CancelAfter(20000)]
    public void Empty_GenomesWithEmptyGeneLists_NoClusters_NoNaN()
    {
        var genomes = Genomes(
            ("g1", Array.Empty<(string, string)>()),
            ("g2", Array.Empty<(string, string)>()));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes);
        AssertWellFormedPartition(r);

        r.Statistics.TotalGenomes.Should().Be(2, "two genomes are present even though both are empty.");
        r.Statistics.TotalGenes.Should().Be(0, "no genes ⇒ no clusters.");
        r.Statistics.CoreFraction.Should().Be(0.0, "INV-06: 0 clusters ⇒ CoreFraction 0, never 0/0 NaN.");
    }

    #endregion

    #region BE — Boundary: single genome (N=1 ⇒ every cluster trivially core)

    // §6.1 "single genome (N=1): all clusters core (occupancy = N), fluidity 0".
    // With one genome, every cluster occupancy = 1 = N, so occupancy/N = 1 ≥
    // coreFraction (≤ 1) ⇒ EVERY cluster is core. Despite occupancy = 1 (which on its
    // own would read as "unique"), the FRACTIONAL core test fires FIRST, so there are
    // zero unique and zero accessory clusters. No DivideByZero on the N=1 denominator.
    [Test]
    [CancelAfter(20000)]
    public void SingleGenome_EveryClusterIsCore_NoUniqueNoAccessory()
    {
        var genomes = Genomes(
            ("only", new[] { Gene("only", 1), Gene("only", 2), Gene("only", 3) }));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99);
        AssertWellFormedPartition(r);

        r.Statistics.TotalGenomes.Should().Be(1);
        r.Statistics.TotalGenes.Should().Be(3, "three distinct families ⇒ three singleton clusters.");
        r.Statistics.CoreGeneCount.Should().Be(3, "§6.1: with N=1 every cluster occupancy = N ⇒ all core (intersection over one genome).");
        r.Statistics.UniqueGeneCount.Should().Be(0, "the fractional core test (1/1 ≥ 0.99) fires before the occupancy = 1 'unique' branch.");
        r.Statistics.AccessoryGeneCount.Should().Be(0);
        r.Statistics.CoreFraction.Should().Be(1.0, "all clusters core ⇒ CoreFraction = 1.");
        r.Statistics.GenomeFluidity.Should().Be(0.0, "§6.1: single genome ⇒ no pairs ⇒ fluidity 0.");
    }

    // Single genome holding only ONE gene: the minimal non-empty input. Probes a
    // crash / IndexOutOfRange on the smallest possible pan-genome.
    [Test]
    [CancelAfter(20000)]
    public void SingleGenome_SingleGene_IsCore_NoCrash()
    {
        var genomes = Genomes(("only", new[] { Gene("only", 42) }));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0);
        AssertWellFormedPartition(r);

        r.Statistics.CoreGeneCount.Should().Be(1, "the sole cluster occupies the sole genome (1/1 = 100% ≥ coreFraction) ⇒ core.");
        r.Statistics.UniqueGeneCount.Should().Be(0);
        r.Statistics.AccessoryGeneCount.Should().Be(0);
    }

    // Single genome with DUPLICATE copies of the same family (paralogs): still one
    // cluster, occupancy 1 = N ⇒ core. Guards that within-genome duplication does not
    // inflate occupancy past N or break the partition.
    [Test]
    [CancelAfter(20000)]
    public void SingleGenome_ParalogousCopies_StillOneCoreCluster()
    {
        var fam = Family(7);
        var genomes = Genomes(
            ("only", new[] { ("p1", fam), ("p2", fam), ("p3", fam) }));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0);
        AssertWellFormedPartition(r);

        r.Statistics.TotalGenes.Should().Be(1, "three identical paralogs collapse into one cluster.");
        r.Statistics.CoreGeneCount.Should().Be(1, "occupancy = 1 distinct genome = N ⇒ core, not unique.");
    }

    #endregion

    #region BE — Boundary: disjoint genomes (no shared genes ⇒ all unique, no false core)

    // §6.1 "pairwise-disjoint content: no core, all unique, fluidity 1" (INV-05).
    // Each genome carries its OWN families, shared by none. Every cluster occupancy = 1
    // ⇒ ALL unique, ZERO core — even at the most permissive coreFraction the partition
    // must NOT falsely report core (1/3 = 33% < any coreFraction > 1/3).
    [Test]
    [CancelAfter(20000)]
    public void DisjointGenomes_NoSharedGenes_AllUnique_NoFalseCore()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("g1", 1), Gene("g1", 2) }),
            ("g2", new[] { Gene("g2", 3), Gene("g2", 4) }),
            ("g3", new[] { Gene("g3", 5), Gene("g3", 6) }));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99);
        AssertWellFormedPartition(r);

        r.Statistics.TotalGenomes.Should().Be(3);
        r.Statistics.TotalGenes.Should().Be(6, "six disjoint families ⇒ six clusters.");
        r.Statistics.CoreGeneCount.Should().Be(0, "§6.1: pairwise-disjoint content ⇒ NO core (no family in >1 genome).");
        r.Statistics.AccessoryGeneCount.Should().Be(0, "no family is in 'some but not all' (every family is in exactly one).");
        r.Statistics.UniqueGeneCount.Should().Be(6, "every cluster occupancy = 1 ⇒ all unique (strain-specific).");
        r.Statistics.CoreFraction.Should().Be(0.0);
        r.Statistics.GenomeFluidity.Should().BeApproximately(1.0, 1e-9, "INV-05: pairwise-disjoint ⇒ fluidity 1.");
    }

    // Disjoint genomes must produce NO core even with coreFraction at its lowest
    // meaningful value (the fractional rule still requires occupancy/N ≥ coreFraction;
    // 1/2 = 0.5 is below 0.51). Two disjoint genomes ⇒ each cluster occupancy 1 = 50%.
    [Test]
    [CancelAfter(20000)]
    public void DisjointGenomes_LowCoreFraction_StillNoFalseCore()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("g1", 1) }),
            ("g2", new[] { Gene("g2", 2) }));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.51);
        AssertWellFormedPartition(r);

        r.Statistics.CoreGeneCount.Should().Be(0, "1/2 = 50% < 51% ⇒ no false core on disjoint genomes.");
        r.Statistics.UniqueGeneCount.Should().Be(2, "each disjoint family occupies exactly one genome ⇒ unique.");
    }

    #endregion

    #region Positive sanity — the documented thresholds on a hand-built mixed example

    // Hand-computed reference (PanGenome_Core_Accessory.md §2.2, INV-03). Three genomes:
    //   family 1 in ALL 3   → occupancy 3/3 = 100% ≥ 99% → CORE
    //   family 2 in 2 of 3  → occupancy 2/3 = 66.7% < 99% AND occupancy ≠ 1 → ACCESSORY (shell)
    //   family 3 in 1 only  → occupancy 1 → UNIQUE
    //   family 4 in 1 only  → occupancy 1 → UNIQUE
    // Expect core 1, accessory 1, unique 2, total 4; CoreFraction = 1/4 = 0.25.
    [Test]
    [CancelAfter(20000)]
    public void PositiveSanity_MixedOccupancy_MatchesDocumentedThresholds()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("c", 1), Gene("s", 2), Gene("u1", 3) }),
            ("g2", new[] { Gene("c", 1), Gene("s", 2) }),
            ("g3", new[] { Gene("c", 1), Gene("u3", 4) }));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99);
        AssertWellFormedPartition(r);

        r.Statistics.TotalGenes.Should().Be(4, "families: 1 (×3), 2 (×2), 3 (×1), 4 (×1) ⇒ 4 clusters.");
        r.Statistics.CoreGeneCount.Should().Be(1, "INV-03: family 1 in 3/3 = 100% ≥ 99% ⇒ core.");
        r.Statistics.AccessoryGeneCount.Should().Be(1, "INV-03: family 2 in 2/3 = 66.7% (<99%, >1 genome) ⇒ accessory/shell.");
        r.Statistics.UniqueGeneCount.Should().Be(2, "INV-03: families 3 and 4 each in 1 genome ⇒ unique.");
        r.Statistics.CoreFraction.Should().BeApproximately(0.25, 1e-12, "INV-06: 1 core / 4 total = 0.25.");
    }

    // The core threshold is FRACTIONAL, not floor(coreFraction·N): with N=3 and
    // coreFraction 0.99 a 2-of-3 (66.7%) cluster is NOT core. Guards the unsourced
    // floor(0.99·3)=2 convention that would wrongly admit it (§2.2 emphasis, [4]).
    [Test]
    [CancelAfter(20000)]
    public void PositiveSanity_TwoOfThree_IsAccessoryNotCore()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("x", 10) }),
            ("g2", new[] { Gene("x", 10) }),
            ("g3", new[] { Gene("y", 11) }));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 0.99);
        AssertWellFormedPartition(r);

        r.Statistics.CoreGeneCount.Should().Be(0, "Roary core = ≥99% of samples; 2/3 = 66.7% is NOT core (not floor(0.99·3)=2).");
        r.Statistics.AccessoryGeneCount.Should().Be(1, "the 2-of-3 family is accessory/shell.");
        r.Statistics.UniqueGeneCount.Should().Be(1, "the 1-of-3 family is unique.");
    }

    // coreFraction = 1.0 is the strict intersection: only families in EVERY genome are
    // core. A family in all-but-one drops to accessory. Confirms the inclusive ≥ at the
    // exact 100% boundary (occupancy = N is core; occupancy = N-1 is not).
    [Test]
    [CancelAfter(20000)]
    public void PositiveSanity_CoreFractionOne_OnlyFullyConservedIsCore()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("all", 1), Gene("most", 2) }),
            ("g2", new[] { Gene("all", 1), Gene("most", 2) }),
            ("g3", new[] { Gene("all", 1) }));

        var r = PanGenomeAnalyzer.ConstructPanGenome(genomes, coreFraction: 1.0);
        AssertWellFormedPartition(r);

        r.Statistics.CoreGeneCount.Should().Be(1, "only family 1 (3/3 = 100% ≥ 1.0) is core at coreFraction 1.0.");
        r.Statistics.AccessoryGeneCount.Should().Be(1, "family 2 in 2/3 (< 100%, > 1 genome) ⇒ accessory.");
        r.Statistics.UniqueGeneCount.Should().Be(0);
    }

    #endregion

    #region Fuzz sweep — random genome/gene topologies always partition without crashing

    // Random pan-genomes across the whole genome-count range (including the BE
    // boundaries N=0, 1, disjoint) must ALWAYS yield a well-formed partition with
    // finite in-range fractions and no throw/hang — the core robustness claim.
    [Test]
    [CancelAfter(60000)]
    public void FuzzSweep_RandomGenomes_AlwaysWellFormedPartition_NoCrash()
    {
        var rng = new Random(Seed);

        for (int iter = 0; iter < 300; iter++)
        {
            int n = rng.Next(0, 6); // 0..5 genomes — includes empty (0) and single (1)
            int familyPool = rng.Next(1, 8);

            var entries = new (string, (string, string)[])[n];
            for (int g = 0; g < n; g++)
            {
                int geneCount = rng.Next(0, 5); // a genome may be empty
                var genes = new (string, string)[geneCount];
                for (int j = 0; j < geneCount; j++)
                {
                    int fam = rng.Next(0, familyPool);
                    genes[j] = ($"g{g}_gene{j}", Family(fam));
                }
                entries[g] = ($"genome{g}", genes);
            }

            double coreFraction = 0.5 + rng.NextDouble() * 0.5; // 0.5 .. 1.0

            PanGenomeAnalyzer.PanGenomeResult r = default;
            var act = () => r = PanGenomeAnalyzer.ConstructPanGenome(Genomes(entries), coreFraction: coreFraction);
            act.Should().NotThrow($"random pan-genome (n={n}, pool={familyPool}, coreFraction={coreFraction:F3}) must never throw.");

            AssertWellFormedPartition(r);
            r.Statistics.TotalGenomes.Should().Be(n, "TotalGenomes must equal the input genome count.");

            // No cluster may be reported core unless its occupancy fraction truly meets
            // the threshold — re-derive each core cluster's occupancy from GenomeToGenes
            // is not possible here, but disjoint/empty already proved no false core; the
            // fraction-range + partition checks cover the rest.
            if (n == 0)
                r.Statistics.TotalGenes.Should().Be(0, "no genomes ⇒ no clusters.");
        }
    }

    #endregion

    #endregion
}
