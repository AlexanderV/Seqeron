namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Comparative-Genomics area — conserved gene clustering
/// (COMPGEN-CLUSTER-001), the common-interval gene-cluster model.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (a cluster smaller than the
/// documented minimum, a label that is not contiguous in every genome, a
/// non-deterministic ordering), and no *unhandled* runtime exception
/// (IndexOutOfRange on a 1-gene genome, DivideByZero on a single genome,
/// NullReferenceException on a gene with no ortholog group). Every input must
/// resolve to EITHER a well-defined, theory-correct result, OR a *documented,
/// intentional* validation exception (ArgumentNullException). A raw runtime
/// exception, a hang, a malformed cluster, or a non-deterministic result is a bug,
/// not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-CLUSTER-001 — Conserved gene clusters (common intervals)
/// Checklist: docs/checklists/03_FUZZING.md, row 132.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate clustering boundaries called out
///          in the checklist row: a single genome, genomes with no conserved genes
///          (no common interval), and identical genomes (every window conserved).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The clustering contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A conserved gene cluster is a COMMON INTERVAL of the genomes' ortholog-group
/// orderings (Bui-Xuan, Habib &amp; Paul 2013 Def. 1, citing Uno &amp; Yagiura 2000):
/// a set of ortholog-group labels that occupies a contiguous window (interval) in
/// EVERY genome — no foreign group may sit between members. Each genome is read as
/// the sequence of its genes' ortholog-group labels in chromosomal order; a gene
/// with no group is a window-breaking non-member. The API entry under test is
///   ComparativeGenomics.FindConservedClusters(
///       IReadOnlyList&lt;IReadOnlyList&lt;Gene&gt;&gt; genomes,
///       IReadOnlyDictionary&lt;string,string&gt; orthologGroups,
///       int minClusterSize = 3, int maxGap = 2)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs
///    lines 914–989; the private interval test IsIntervalOf lines 1000–1021).
/// Each returned cluster is the sorted list of its ortholog-group labels; the
/// collection is ordered by cluster size then lexicographically by joined labels.
///
/// THE DOCUMENTED INVARIANTS (Conserved_Gene_Clusters.md §2.4):
///   • INV-01: every reported cluster is a label set contiguous in EVERY genome.
///   • INV-02: every cluster has size ≥ minClusterSize AND ≥ 2 (interval needs i &lt; j).
///   • INV-03: a set contiguous in some but not all genomes is excluded.
///   • INV-04: fewer than 2 genomes ⇒ empty result (common interval is a K ≥ 2 notion).
///   • INV-05: output is deterministic and order-independent.
/// The documented edge cases (§6.1):
///   • &lt; 2 genomes → empty; gene without ortholog group → breaks windows, never a
///     member; repeated group label → any matching window counts; no conserved set
///     ≥ minClusterSize → empty; null genomes / orthologGroups → ArgumentNullException.
///
/// The three BE checklist targets map to these documented behaviours:
///   • single genome     → empty (INV-04): one genome is not a family, so the
///                          conserved-cluster question is vacuous. The hazard probed
///                          is a crash / DivideByZero / IndexOutOfRange on a 1-element
///                          clustering — the `genomes.Count &lt; 2` guard short-circuits
///                          to empty before any window scan.
///   • no conserved genes → empty (INV-03/§6.1): genomes that share no contiguous
///                          label set of size ≥ minClusterSize produce no common
///                          interval; nothing is falsely clustered, and no crash.
///   • identical genomes  → every window of the shared sequence is conserved (INV-01),
///                          deterministically, with no infinite loop (the worst-case
///                          quadratic window scan still terminates).
/// A positive-sanity test pins the documented golden vector (§7.1 / Example 1) and
/// that two tight groups with a clear separating foreign group reconstruct as
/// distinct clusters; identical genomes collapse to the full set of conserved windows.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeClusterFuzzTests
{
    #region Helpers

    /// <summary>
    /// Builds a genome (ordered gene list) from a string of single-character group
    /// labels. Gene ids are made unique per genome so the ortholog-group map keys do
    /// not collide across genomes; the map maps every gene id to its label.
    /// </summary>
    private static (IReadOnlyList<ComparativeGenomics.Gene> genome, Dictionary<string, string> map) BuildGenome(
        string labels, int genomeIndex)
    {
        var genes = new List<ComparativeGenomics.Gene>(labels.Length);
        var map = new Dictionary<string, string>(StringComparer.Ordinal);
        for (int i = 0; i < labels.Length; i++)
        {
            string geneId = $"g{genomeIndex}_{i}";
            genes.Add(new ComparativeGenomics.Gene(geneId, $"G{genomeIndex}", i * 100, i * 100 + 50, '+'));
            map[geneId] = labels[i].ToString();
        }
        return (genes, map);
    }

    /// <summary>
    /// Assembles several label-string genomes into the (genomes, combined-map) pair the
    /// API expects. Genomes are 0-indexed so per-genome gene ids stay unique.
    /// </summary>
    private static (List<IReadOnlyList<ComparativeGenomics.Gene>> genomes, Dictionary<string, string> map) Build(
        params string[] labelStrings)
    {
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>>(labelStrings.Length);
        var map = new Dictionary<string, string>(StringComparer.Ordinal);
        for (int g = 0; g < labelStrings.Length; g++)
        {
            var (genome, gmap) = BuildGenome(labelStrings[g], g);
            genomes.Add(genome);
            foreach (var kv in gmap) map[kv.Key] = kv.Value;
        }
        return (genomes, map);
    }

    /// <summary>Canonical key (ordinal-joined sorted labels) for one cluster.</summary>
    private static string Key(IReadOnlyList<string> cluster) => string.Join("", cluster);

    /// <summary>The set of canonical cluster keys returned by the unit.</summary>
    private static HashSet<string> KeySet(IEnumerable<IReadOnlyList<string>> clusters) =>
        clusters.Select(Key).ToHashSet(StringComparer.Ordinal);

    /// <summary>
    /// A well-formed clustering result satisfies the documented structural invariants:
    ///   • every cluster is non-empty and has size ≥ 2 (INV-02; the interval lower bound),
    ///   • every cluster's labels are sorted ordinally and distinct (output shape §3.2),
    ///   • the returned collection is sorted by size then joined-label order (INV-05 shape),
    ///   • no two clusters share the same canonical label set (sets are deduplicated).
    /// This is the partition-style sanity net: it does not by itself verify INV-01
    /// (contiguity in every genome) — the targeted tests pin that — but it catches a
    /// degenerate or malformed cluster (size &lt; 2, unsorted, duplicated).
    /// </summary>
    private static void AssertWellFormedClusters(IReadOnlyList<IReadOnlyList<string>> clusters, int effectiveMin)
    {
        var keys = new HashSet<string>(StringComparer.Ordinal);
        IReadOnlyList<string>? prev = null;
        foreach (var c in clusters)
        {
            c.Should().NotBeNull("a cluster must never be null");
            c.Count.Should().BeGreaterThanOrEqualTo(Math.Max(effectiveMin, 2),
                "every cluster has size ≥ max(minClusterSize, 2) (INV-02)");
            c.Should().BeInAscendingOrder(StringComparer.Ordinal, "cluster labels are ordinal-sorted (§3.2)");
            c.Should().OnlyHaveUniqueItems("a cluster is a SET of labels, no duplicates (§3.2)");

            keys.Add(Key(c)).Should().BeTrue("distinct cluster sets must not be reported twice (deduplicated)");

            if (prev is not null)
            {
                int sizeCmp = prev.Count.CompareTo(c.Count);
                sizeCmp.Should().BeLessThanOrEqualTo(0, "clusters are ordered by ascending size (INV-05)");
                if (sizeCmp == 0)
                    string.CompareOrdinal(Key(prev), Key(c)).Should().BeLessThan(0,
                        "equal-size clusters are ordered lexicographically by joined labels (INV-05)");
            }
            prev = c;
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  COMPGEN-CLUSTER-001 — Conserved gene clusters : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region COMPGEN-CLUSTER-001 — Conserved gene clusters (common intervals)

    #region BE — Boundary: single genome (family-of-one is vacuous)

    /// <summary>
    /// BE: a single genome is the degenerate floor of the genome-count axis. A common
    /// interval is a K ≥ 2 family notion (Conserved_Gene_Clusters.md §2.4 INV-04,
    /// §6.1): with one genome every window is trivially "common", so the question is
    /// vacuous and the documented answer is EMPTY. The hazard this probes is a crash,
    /// DivideByZero, or IndexOutOfRange on a 1-element clustering — the
    /// `genomes.Count &lt; 2` guard (ComparativeGenomics.cs lines 926–927) short-circuits
    /// to empty BEFORE any window scan, so a genome of any length is safe.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_SingleGenome_IsEmptyNoCrash()
    {
        var (genomes, map) = Build("abcde");

        List<IReadOnlyList<string>> clusters = null!;
        FluentActions.Invoking(() => clusters = ComparativeGenomics.FindConservedClusters(genomes, map).ToList())
            .Should().NotThrow("a single genome must short-circuit to empty, never crash on a family-of-one");

        clusters.Should().BeEmpty("with fewer than 2 genomes the common-interval question is vacuous → empty (INV-04)");
    }

    /// <summary>
    /// BE: a single 1-gene genome is the smallest possible non-empty single-genome
    /// input — the corner where a naive window scan could take a degenerate span. It
    /// must still hit the &lt; 2-genome guard and return empty with NO IndexOutOfRange on
    /// the length-1 label sequence.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_SingleSingleGeneGenome_IsEmptyNoIndexOutOfRange()
    {
        var (genomes, map) = Build("a");

        FluentActions.Invoking(() => ComparativeGenomics.FindConservedClusters(genomes, map).ToList())
            .Should().NotThrow("a 1-gene single genome must not take an over-long span");

        ComparativeGenomics.FindConservedClusters(genomes, map).Should().BeEmpty(
            "still fewer than 2 genomes → empty (INV-04)");
    }

    /// <summary>
    /// BE: the empty genome LIST is the floor below a single genome. Zero genomes is
    /// also fewer than two, so the result is empty — no exception, no division by a
    /// zero genome count.
    /// </summary>
    [Test]
    public void FindConservedClusters_EmptyGenomeList_IsEmptyNoCrash()
    {
        var genomes = new List<IReadOnlyList<ComparativeGenomics.Gene>>();
        var map = new Dictionary<string, string>();

        FluentActions.Invoking(() => ComparativeGenomics.FindConservedClusters(genomes, map).ToList())
            .Should().NotThrow("zero genomes is fewer than two → empty, never a DivideByZero on the genome count");

        ComparativeGenomics.FindConservedClusters(genomes, map).Should().BeEmpty(
            "zero genomes → empty (INV-04)");
    }

    #endregion

    #region BE — Boundary: no conserved genes (no common interval)

    /// <summary>
    /// BE: genomes that share NO contiguous label set of size ≥ minClusterSize have no
    /// common interval → empty, with nothing falsely clustered (INV-03,
    /// Conserved_Gene_Clusters.md §6.1 "no conserved set ≥ minClusterSize → empty").
    /// Two genomes built from completely disjoint group alphabets cannot share any
    /// label at all, let alone a contiguous set, so the result is empty.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_DisjointAlphabets_IsEmpty()
    {
        // Genome 0 uses {a,b,c}, genome 1 uses {x,y,z}: no shared label → no common interval.
        var (genomes, map) = Build("abc", "xyz");

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2).ToList();

        AssertWellFormedClusters(clusters, effectiveMin: 2);
        clusters.Should().BeEmpty("genomes sharing no labels share no contiguous set → no common interval (INV-03)");
    }

    /// <summary>
    /// BE: shared labels that are SCRAMBLED so that no set of size ≥ minClusterSize is
    /// contiguous in both genomes also yields empty — the discriminating "no conserved
    /// genes" case where the alphabets overlap but the arrangement defeats every
    /// candidate window. Here genome 1 reverses-with-interleave the order so that the
    /// only contiguous size-3 set of genome 0 ({a,b,c}) is broken by a foreign element
    /// in genome 1. This pins that overlap alone does NOT manufacture a cluster (INV-03).
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_OverlappingButNeverContiguous_IsEmpty()
    {
        // Genome 0: a b c   (the only size-3 window is {a,b,c})
        // Genome 1: a x b x c — a,b,c never form a foreign-free contiguous window.
        var (genomes, map) = Build("abc", "axbxc");

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3).ToList();

        AssertWellFormedClusters(clusters, effectiveMin: 3);
        clusters.Should().BeEmpty("the size-3 set {a,b,c} is contiguous in genome 0 but split by 'x' in genome 1 → excluded (INV-03)");
    }

    /// <summary>
    /// BE: genes with NO ortholog group are window-breaking non-members (§6.1). When
    /// every gene of one genome lacks a group, no label set can be an interval of it,
    /// so the result is empty — and crucially there is NO NullReferenceException on the
    /// missing-key lookup (the code substitutes a non-member sentinel, not null).
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_GenomeWithNoOrthologGroups_IsEmptyNoNullRef()
    {
        // Build two genomes, then DROP every map entry for genome 1 so its genes are ungrouped.
        var (genomes, map) = Build("abc", "abc");
        foreach (var gene in genomes[1])
            map.Remove(gene.Id);

        List<IReadOnlyList<string>> clusters = null!;
        FluentActions.Invoking(() => clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2).ToList())
            .Should().NotThrow("an ungrouped gene must be a window-breaking non-member, never a null deref");

        clusters.Should().BeEmpty("every gene of genome 1 is ungrouped → no label set is an interval of it → empty (§6.1)");
    }

    #endregion

    #region BE — Boundary: identical genomes (every window conserved)

    /// <summary>
    /// BE: identical genomes are the upper boundary — every contiguous window of the
    /// shared sequence is, by construction, a common interval. For the distinct-label
    /// permutation "abcd" replicated across two genomes the conserved clusters are
    /// exactly the windows of size ≥ minClusterSize: {a,b,c}, {b,c,d}, {a,b,c,d}. We pin
    /// the EXACT set (not "contains"), proving the upper boundary enumerates every
    /// window and invents none. The scan terminates (CancelAfter) — no infinite loop on
    /// the worst-case all-conserved input.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_IdenticalGenomes_AllWindowsConserved()
    {
        var (genomes, map) = Build("abcd", "abcd");

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3).ToList();

        AssertWellFormedClusters(clusters, effectiveMin: 3);
        KeySet(clusters).Should().BeEquivalentTo(new[] { "abc", "bcd", "abcd" },
            "for identical 'abcd' genomes every size-≥3 contiguous window is a common interval: {a,b,c}, {b,c,d}, {a,b,c,d}");
    }

    /// <summary>
    /// BE: identical genomes replicated across MANY copies must give the SAME conserved
    /// set as the two-copy case — adding identical genomes cannot add or remove a common
    /// interval (each is already an interval of the shared sequence). This pins INV-01
    /// at the upper boundary and that the K-genome verification loop does not spuriously
    /// drop a window. The all-identical input is also the canonical non-termination
    /// hazard (maximal candidate set); CancelAfter guards it.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_FiveIdenticalGenomes_SameAsTwoCopies()
    {
        var (two, twoMap) = Build("abcd", "abcd");
        var (five, fiveMap) = Build("abcd", "abcd", "abcd", "abcd", "abcd");

        var twoClusters = ComparativeGenomics.FindConservedClusters(two, twoMap, minClusterSize: 3).ToList();
        var fiveClusters = ComparativeGenomics.FindConservedClusters(five, fiveMap, minClusterSize: 3).ToList();

        AssertWellFormedClusters(fiveClusters, effectiveMin: 3);
        KeySet(fiveClusters).Should().BeEquivalentTo(KeySet(twoClusters),
            "identical genomes are intervals of one another, so adding copies leaves the conserved set unchanged (INV-01)");
    }

    /// <summary>
    /// BE: identical genomes collapse a REPEATED-LABEL sequence into one cluster per
    /// distinct label set, not one per position. For "aabbcc" replicated, the whole-set
    /// cluster {a,b,c} is reported once (cluster identity is the SET of labels;
    /// duplicate windows collapse, §5.2). Pins that the upper boundary deduplicates
    /// rather than emitting a cluster per occurrence.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_IdenticalGenomesRepeatedLabels_CollapseToDistinctSets()
    {
        var (genomes, map) = Build("aabbcc", "aabbcc");

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3).ToList();

        AssertWellFormedClusters(clusters, effectiveMin: 3);
        KeySet(clusters).Should().Contain("abc",
            "the full distinct set {a,b,c} is a common interval of the identical repeated sequence");
        KeySet(clusters).Should().OnlyHaveUniqueItems("each distinct label SET is reported once, not once per occurrence (§5.2)");
    }

    #endregion

    #region BE — Boundary: degenerate minClusterSize (interval lower bound)

    /// <summary>
    /// BE: a degenerate minClusterSize ≤ 0 is the floor of the cluster-size axis. The
    /// interval definition requires i &lt; j, so size is raised to ≥ 2 internally
    /// (`effectiveMin = Math.Max(minClusterSize, 2)`, ComparativeGenomics.cs line 929,
    /// INV-02). A request for size 0 or a negative size must therefore NOT produce a
    /// size-0/1 (single-element or empty) cluster and must not crash; for identical
    /// "abc" genomes the smallest valid clusters are the size-2 windows plus the full
    /// set: {a,b}, {b,c}, {a,b,c}.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_NonPositiveMinClusterSize_FloorsToTwo()
    {
        var (genomes, map) = Build("abc", "abc");

        foreach (int degenerate in new[] { int.MinValue, -5, 0, 1 })
        {
            var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: degenerate).ToList();

            AssertWellFormedClusters(clusters, effectiveMin: 2);
            KeySet(clusters).Should().BeEquivalentTo(new[] { "ab", "bc", "abc" },
                $"minClusterSize {degenerate} is floored to 2 (interval needs i<j, INV-02), so no size-0/1 cluster appears");
        }
    }

    /// <summary>
    /// BE: a minClusterSize LARGER than any genome cannot be satisfied — no window has
    /// that many distinct labels — so the result is empty (no conserved set ≥
    /// minClusterSize, §6.1). Pins the upper end of the size axis: an over-large
    /// threshold floors the result to empty rather than crashing or wrapping.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_MinClusterSizeAboveGenomeLength_IsEmpty()
    {
        var (genomes, map) = Build("abc", "abc");

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 99).ToList();

        clusters.Should().BeEmpty("no window has 99 distinct labels, so no cluster meets the threshold → empty (§6.1)");
    }

    #endregion

    #region BE — Validation: null contract

    /// <summary>
    /// BE: null genomes or orthologGroups is the one input the contract rejects with an
    /// exception (Conserved_Gene_Clusters.md §6.1; ComparativeGenomics.cs lines 920–921).
    /// We pin ArgumentNullException carrying the documented parameter name for each, so a
    /// null argument is a disciplined defensive throw, not an unhandled NRE deep in the
    /// scan.
    /// </summary>
    [Test]
    public void FindConservedClusters_NullArguments_ThrowArgumentNull()
    {
        var (genomes, map) = Build("abc", "abc");

        FluentActions.Invoking(() => ComparativeGenomics.FindConservedClusters(null!, map).ToList())
            .Should().Throw<ArgumentNullException>("null genomes is rejected defensively")
            .Which.ParamName.Should().Be("genomes");

        FluentActions.Invoking(() => ComparativeGenomics.FindConservedClusters(genomes, null!).ToList())
            .Should().Throw<ArgumentNullException>("null orthologGroups is rejected defensively")
            .Which.ParamName.Should().Be("orthologGroups");
    }

    #endregion

    #region Positive sanity — the documented metric on hand-built examples

    /// <summary>
    /// Positive sanity (Conserved_Gene_Clusters.md §7.1, Bui-Xuan et al. 2013 Example 1):
    /// pins the EXACT golden vector, not just "green". For P1 = (1 2 3 4 5 6 7) and
    /// P2 = (7 2 1 3 6 4 5) the common intervals of size ≥ 2 are EXACTLY
    ///   {1,2}, {1,2,3}, {3,4,5,6}, {4,5}, {4,5,6}, {1..6}, {1..7}.
    /// This is the load-bearing correctness check that the common-interval enumeration
    /// is theory-correct, distinguishing it from an over- or under-reporting clusterer.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_GoldenVectorExample1_ExactCommonIntervals()
    {
        // Labels "1".."7"; genome 0 = identity, genome 1 = the paper's permutation.
        var (genomes, map) = Build("1234567", "7213645");

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2).ToList();

        AssertWellFormedClusters(clusters, effectiveMin: 2);
        KeySet(clusters).Should().BeEquivalentTo(
            new[] { "12", "123", "3456", "45", "456", "123456", "1234567" },
            "the size-≥2 common intervals of P1=Id7, P2=(7 2 1 3 6 4 5) are exactly the paper's Example 1 set (§7.1)");
    }

    /// <summary>
    /// Positive sanity — two TIGHT conserved groups with a clear separation reconstruct
    /// as DISTINCT clusters, the core business meaning of conserved gene clustering. Two
    /// operons {a,b,c} and {d,e,f} are kept internally contiguous in both genomes, but a
    /// DIFFERENT foreign group is wedged between them in each genome ('y' in genome 0,
    /// 'z' in genome 1) and the operons are reordered. Because the separators differ,
    /// NO cross-operon union (e.g. {a,b,c} ∪ separator ∪ {d,e,f}) can be a foreign-free
    /// window in BOTH genomes, so the two operons are each reported but no spurious
    /// merged cluster — the "clear gap" reconstructs the two groups exactly (INV-03).
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_TwoTightGroupsWithGap_ReconstructsBothNoSpuriousMerge()
    {
        // Genome 0: a b c  y  d e f      operon1, separator y, operon2
        // Genome 1: d e f  z  a b c      operons swapped, DIFFERENT separator z
        var (genomes, map) = Build("abcydef", "defzabc");

        var clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 3).ToList();

        AssertWellFormedClusters(clusters, effectiveMin: 3);
        var keys = KeySet(clusters);
        keys.Should().Contain("abc", "operon {a,b,c} is contiguous (foreign-free) in both genomes → a conserved cluster");
        keys.Should().Contain("def", "operon {d,e,f} is contiguous in both genomes → a conserved cluster");
        keys.Should().BeEquivalentTo(new[] { "abc", "def" },
            "the separators differ ('y' vs 'z'), so no cross-operon set is a foreign-free window in both genomes → exactly the two operons are reported (INV-03)");
    }

    /// <summary>
    /// Positive sanity — genome ORDER independence (INV-05). Permuting the order of the
    /// input genomes must not change the conserved-cluster set: the common-interval
    /// notion is symmetric in the family. Pins that the result is determined by the
    /// genomes' contents, not their argument order.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_GenomeOrderIndependent()
    {
        var (forward, fmap) = Build("abcydef", "defzabc");
        var (reversed, rmap) = Build("defzabc", "abcydef");

        var forwardKeys = KeySet(ComparativeGenomics.FindConservedClusters(forward, fmap, minClusterSize: 3));
        var reversedKeys = KeySet(ComparativeGenomics.FindConservedClusters(reversed, rmap, minClusterSize: 3));

        reversedKeys.Should().BeEquivalentTo(forwardKeys,
            "the common-interval set is symmetric in the genome family → order-independent (INV-05)");
    }

    /// <summary>
    /// Positive sanity — DETERMINISM (INV-05). Repeated invocations on the same input
    /// must return byte-for-byte the same ordered cluster list, with no reliance on hash
    /// iteration order. Pins that the documented deterministic sort (by size then
    /// joined-label order) makes the output reproducible.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindConservedClusters_RepeatedRuns_AreDeterministic()
    {
        var (genomes, map) = Build("1234567", "7213645");

        List<string> Run() => ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: 2)
            .Select(Key).ToList();

        var first = Run();
        for (int i = 0; i < 5; i++)
            Run().Should().Equal(first, "the deterministic sort makes the ordered output identical across runs (INV-05)");
    }

    /// <summary>
    /// Fuzz sweep — across many random multi-genome families, varied lengths, alphabet
    /// sizes, minClusterSize values and ungrouped-gene rates, FindConservedClusters must
    /// ALWAYS terminate, never throw, and return only well-formed clusters (size ≥ 2,
    /// sorted, deduplicated, contiguous in EVERY genome). The contiguity check here is
    /// the strong net: we independently re-verify INV-01 for each reported cluster.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void FindConservedClusters_RandomFamilies_AlwaysWellFormedAndContiguous()
    {
        const string alphabet = "abcde";
        for (int seed = 0; seed < 60; seed++)
        {
            var rng = new Random(seed);
            int genomeCount = rng.Next(1, 5);          // includes the single-genome degenerate case
            int len = rng.Next(0, 9);                  // includes the empty genome
            int minSize = rng.Next(-1, 5);             // includes degenerate ≤ 0

            // All genomes are permutations/sequences over the same small alphabet so
            // conserved clusters can genuinely occur; some genes are randomly ungrouped.
            var labelStrings = new string[genomeCount];
            for (int g = 0; g < genomeCount; g++)
            {
                var sb = new char[len];
                for (int i = 0; i < len; i++) sb[i] = alphabet[rng.Next(alphabet.Length)];
                labelStrings[g] = new string(sb);
            }
            var (genomes, map) = Build(labelStrings);
            // Randomly drop a fraction of map entries → ungrouped (window-breaking) genes.
            foreach (var genome in genomes)
                foreach (var gene in genome)
                    if (rng.NextDouble() < 0.15) map.Remove(gene.Id);

            List<IReadOnlyList<string>> clusters = null!;
            FluentActions.Invoking(() => clusters = ComparativeGenomics.FindConservedClusters(genomes, map, minClusterSize: minSize).ToList())
                .Should().NotThrow($"random family must never crash (seed {seed})");

            int effectiveMin = Math.Max(minSize, 2);
            AssertWellFormedClusters(clusters, effectiveMin);

            // INV-01: re-verify each reported cluster is contiguous (foreign-free window) in EVERY genome.
            foreach (var cluster in clusters)
            {
                var set = cluster.ToHashSet(StringComparer.Ordinal);
                foreach (var genome in genomes)
                {
                    string[] labels = genome.Select(gn =>
                        map.TryGetValue(gn.Id, out var grp) ? grp : " __none__ ").ToArray();
                    IsContiguousIntervalOf(labels, set).Should().BeTrue(
                        $"reported cluster {{{Key(cluster)}}} must be a foreign-free contiguous window in every genome (INV-01, seed {seed})");
                }
            }
        }
    }

    /// <summary>
    /// Independent re-implementation of the common-interval test for the fuzz sweep's
    /// INV-01 cross-check: is <paramref name="set"/> the exact label set of some
    /// contiguous, foreign-free window of <paramref name="labels"/>? Mirrors the
    /// documented interval definition (set of all elements of a window; no foreign
    /// element) without reusing the unit under test.
    /// </summary>
    private static bool IsContiguousIntervalOf(string[] labels, HashSet<string> set)
    {
        for (int start = 0; start < labels.Length; start++)
        {
            if (!set.Contains(labels[start])) continue;
            var seen = new HashSet<string>(StringComparer.Ordinal);
            for (int end = start; end < labels.Length; end++)
            {
                if (!set.Contains(labels[end])) break; // foreign element → not this location
                seen.Add(labels[end]);
                if (seen.Count == set.Count) return true;
            }
        }
        return false;
    }

    #endregion

    #endregion
}
