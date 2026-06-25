using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology tumour clonal-phylogeny reconstruction area — ONCO-PHYLO-001.
/// The unit under test is the deterministic clonal-tree builder implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
/// <see cref="OncologyAnalyzer.ReconstructPhylogeny(IReadOnlyList{OncologyAnalyzer.CcfCluster}, double)"/>
/// — reconstructs a rooted clonal (tumour) phylogeny from per-sample cancer cell fraction (CCF)
/// clusters — together with its trunk/branch partition helpers
/// <see cref="OncologyAnalyzer.IdentifyTrunkMutations(OncologyAnalyzer.ClonalPhylogeny)"/> and
/// <see cref="OncologyAnalyzer.IdentifyBranchMutations(OncologyAnalyzer.ClonalPhylogeny)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang (tree building MUST
/// terminate — no cycles, no clone made its own ancestor), no NullReference on a
/// node with no parent, and no nonsense output. Every input must resolve to EITHER
/// a well-formed, theory-correct rooted tree OR a *documented, intentional* fault —
/// here an ArgumentNullException (null clusters / null CCF list), an ArgumentException
/// (ragged / empty / NaN / out-of-[0,1] CCF, or a duplicate id), or an
/// ArgumentOutOfRangeException (negative / NaN tolerance). The result is always a
/// SINGLE rooted spanning tree: each non-root cluster has exactly one parent, no
/// cycles, every node reachable from the synthetic normal root, and every edge
/// obeys the perfect-phylogeny constraints. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-PHYLO-001 — tumour clonal phylogeny reconstruction (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 114.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 114): "single clone, identical clones, no shared mutations".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// The CCF-cluster model is the cell-fraction analogue of mutation sets: under the
/// infinite-sites / perfect-phylogeny model an ancestor's CCF dominates its
/// descendant's in every sample (Tumor_Phylogeny_Reconstruction.md §2.2), so
/// "ancestor mutation-set ⊆ descendant" is realised as lineage precedence
/// (ancestor CCF ≥ descendant CCF − ε) plus the presence pattern (a parent absent
/// in a sample cannot parent a child present there). The BE targets map onto the
/// documented contract as follows:
///   • single clone     ⇔ ONE cluster ⇒ trivial tree: root → that cluster, no other
///                          edges, no crash; trunk = {cluster}, branches = {}
///                          (Tumor_Phylogeny_Reconstruction.md §6.1, §7.1);
///   • identical clones  ⇔ several clusters with the SAME CCF vector ⇒ deterministic
///                          handling: the sum rule (Eq. 5) forbids two equal-CCF
///                          siblings under one parent, so they chain; no infinite
///                          loop, no duplicate-parent corruption, no self-ancestry
///                          (§6.1 "Two equal-CCF clusters under one parent", INV-02);
///   • no shared muts    ⇔ clusters with DISJOINT presence (each private to a
///                          different sample) ⇒ independent lineages off the trunk;
///                          no false ancestry inferred between disjoint clones
///                          (presence pattern, §2.2 / §7.1).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (Tumor_Phylogeny_Reconstruction.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • Lineage precedence: every edge u→v has u.CCF[i] ≥ v.CCF[i] − ε ∀i, and
///     u.CCF[i] = 0 ⇒ v.CCF[i] = 0 (presence pattern)                  (INV-01, §2.2 Eq.2)
///   • Sum rule: per node, per sample, Σ_children v.CCF[i] ≤ u.CCF[i] + ε (INV-02, §2.2 Eq.5)
///   • Single rooted tree: every cluster has exactly one parent, no cycles (INV-03, §2.4)
///   • Trunk and Branch partition the clusters (disjoint, union = all)  (INV-04)
///   • Deterministic: identical input ⇒ identical edge set            (INV-05, §2.4)
///   • Synthetic root id is distinct from every cluster id; root.CCF = 1 ∀ sample (§3.2, §4.1)
///   • Empty cluster list ⇒ root-only tree (no edges, trunk = {}, branches = {}) (§6.1)
///   • null clusters / null CCF list ⇒ ArgumentNullException          (§3.3)
///   • ragged / empty / NaN / out-of-[0,1] CCF, duplicate id ⇒ ArgumentException (§3.3)
///   • negative / NaN tolerance ⇒ ArgumentOutOfRangeException         (§3.3)
///
/// No source bug was found; no test was weakened.
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// All tree-building tests carry [CancelAfter] so a non-terminating build fails loudly.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyClonalPhylogenyFuzzTests
{
    private const int HangGuardMs = 30_000;

    // ── Well-formed-tree assertion helper ────────────────────────────────────
    // Pin the documented tree-validity contract on EVERY reconstructed phylogeny so
    // a fuzz test cannot rubber-stamp a malformed tree (a missing/duplicate parent,
    // a cycle, a self-ancestry edge, an unreachable node, or an edge that violates
    // lineage precedence / the sum rule) green. This is the structural backbone the
    // BE targets all funnel through.
    private static void AssertWellFormedTree(
        ClonalPhylogeny tree,
        IReadOnlyList<CcfCluster> input,
        double tolerance)
    {
        tree.Clusters.Should().NotBeNull("a phylogeny always carries its cluster list");
        tree.Clusters.Count.Should().Be(input.Count, "every input cluster is a node (in input order)");
        tree.Edges.Should().NotBeNull("a phylogeny always carries its edge list");

        // The synthetic root id must be distinct from every real cluster id (§3.2, §4.1).
        var clusterIds = new HashSet<int>();
        foreach (CcfCluster c in input)
        {
            clusterIds.Add(c.Id);
        }

        clusterIds.Should().NotContain(tree.RootId,
            "the synthetic normal root id must be distinct from every cluster id (§3.2)");

        // One edge per non-root cluster: a spanning tree has exactly n edges (INV-03).
        tree.Edges.Count.Should().Be(input.Count,
            "a rooted spanning tree over n clusters has exactly n edges (INV-03)");

        // Each cluster has EXACTLY one parent — no missing parent (NullReference risk),
        // no duplicate parent (corruption). ParentOf must be non-null for every cluster.
        var parent = new Dictionary<int, int>();
        foreach (CcfCluster c in input)
        {
            int? p = tree.ParentOf(c.Id);
            p.Should().NotBeNull($"cluster {c.Id} must have exactly one parent (no orphan) (INV-03)");
            parent[c.Id] = p!.Value;

            // A clone may never be its own parent / ancestor.
            p!.Value.Should().NotBe(c.Id, $"cluster {c.Id} must not be its own parent (no self-ancestry)");
        }

        // No edge may name a child more than once (single-parent ⇒ no duplicate child).
        var childIds = tree.Edges.Select(e => e.ChildId).ToList();
        childIds.Should().OnlyHaveUniqueItems("each cluster is the child of exactly one edge (INV-03)");
        var childSet = new HashSet<int>(childIds);
        foreach (CcfCluster c in input)
        {
            childSet.Should().Contain(c.Id, $"cluster {c.Id} must appear as a child exactly once");
        }

        // Acyclicity + reachability: walk every node up to the root in a bounded number
        // of hops. A cycle (or self-ancestry) would loop forever; the bound makes that a
        // loud failure instead of a hang inside the assertion.
        int n = input.Count;
        foreach (CcfCluster c in input)
        {
            int hops = 0;
            int cur = c.Id;
            while (cur != tree.RootId)
            {
                parent.ContainsKey(cur).Should().BeTrue(
                    $"node {cur} on the path from {c.Id} must have a recorded parent (reachable, no broken chain)");
                cur = parent[cur];
                hops++;
                hops.Should().BeLessThanOrEqualTo(n + 1,
                    $"the path from cluster {c.Id} to the root must terminate (acyclic, no self-ancestry) (INV-03)");
            }
        }

        // Per-edge perfect-phylogeny constraints and the sum rule (INV-01, INV-02).
        var ccfById = new Dictionary<int, double[]>();
        int sampleCount = tree.SampleCount;
        double[] rootCcf = new double[sampleCount];
        Array.Fill(rootCcf, 1.0);
        ccfById[tree.RootId] = rootCcf;
        foreach (CcfCluster c in input)
        {
            ccfById[c.Id] = c.CcfPerSample.ToArray();
        }

        // INV-01: lineage precedence on every edge.
        foreach (ClonalEdge e in tree.Edges)
        {
            double[] pc = ccfById[e.ParentId];
            double[] cc = ccfById[e.ChildId];
            for (int i = 0; i < sampleCount; i++)
            {
                pc[i].Should().BeGreaterThanOrEqualTo(cc[i] - tolerance - 1e-12,
                    $"edge {e.ParentId}→{e.ChildId}: ancestor CCF ≥ descendant CCF − ε in sample {i} (INV-01)");

                if (pc[i] <= 0.0)
                {
                    cc[i].Should().Be(0.0,
                        $"edge {e.ParentId}→{e.ChildId}: a parent absent in sample {i} cannot parent a child present there (presence pattern, INV-01)");
                }
            }
        }

        // INV-02: per node, per sample, children CCF sum ≤ parent CCF + ε.
        // Scope: REAL cluster parents only. The synthetic root is the documented
        // spanning-tree backstop — when a noisy clone fits no admissible cluster parent
        // it is attached to the root even if that breaches the root's artificial budget
        // of 1.0, because every clone must have a parent (source §5.2, ReconstructPhylogeny
        // fallback comment; Popic 2015 spanning tree). The root is therefore exempt here.
        var childrenByParent = tree.Edges
            .Where(e => e.ParentId != tree.RootId)
            .GroupBy(e => e.ParentId)
            .ToDictionary(g => g.Key, g => g.Select(e => e.ChildId).ToList());
        foreach (KeyValuePair<int, List<int>> kv in childrenByParent)
        {
            double[] pc = ccfById[kv.Key];
            for (int i = 0; i < sampleCount; i++)
            {
                double sum = 0.0;
                foreach (int child in kv.Value)
                {
                    sum += ccfById[child][i];
                }

                sum.Should().BeLessThanOrEqualTo(pc[i] + tolerance + 1e-12,
                    $"node {kv.Key}: children CCF sum ≤ parent CCF + ε in sample {i} (sum rule, INV-02)");
            }
        }

        // INV-04: trunk and branch partition the clusters (disjoint, union = all).
        var trunk = new HashSet<int>(IdentifyTrunkMutations(tree));
        var branch = new HashSet<int>(IdentifyBranchMutations(tree));
        trunk.Overlaps(branch).Should().BeFalse("trunk and branch sets are disjoint (INV-04)");
        var union = new HashSet<int>(trunk);
        union.UnionWith(branch);
        union.Should().BeEquivalentTo(clusterIds, "trunk ∪ branch covers exactly the clusters (INV-04)");
    }

    private static CcfCluster Cluster(int id, params double[] ccf) => new(id, ccf);

    // ═════════════════════════════════════════════════════════════════════════
    #region ONCO-PHYLO-001 — clonal phylogeny reconstruction
    // ═════════════════════════════════════════════════════════════════════════

    // ── BE: single clone ─────────────────────────────────────────────────────
    // One cluster ⇒ trivial tree: root → cluster, exactly one edge, that cluster
    // is the trunk, no branches, no crash (§6.1 "Single cluster", §7.1).

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_SingleClone_FormsTrivialRootedTree()
    {
        var clusters = new[] { Cluster(7, 0.6) };

        ClonalPhylogeny tree = ReconstructPhylogeny(clusters);

        AssertWellFormedTree(tree, clusters, DefaultPhylogenyTolerance);
        tree.ParentOf(7).Should().Be(tree.RootId, "the only clone attaches directly to the synthetic root (§6.1)");
        tree.ChildrenOf(7).Should().BeEmpty("a single clone has no descendants");
        IdentifyTrunkMutations(tree).Should().Equal(new[] { 7 }, "the single clone is the whole trunk (§6.1)");
        IdentifyBranchMutations(tree).Should().BeEmpty("a single clone yields no subclonal branches (§6.1)");
    }

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_SingleClone_AnyInRangeCcfAndId_NoCrash([Range(0, 19)] int seed)
    {
        var rng = new Random(seed);
        int id = rng.Next(2) == 0 ? rng.Next(-1000, 1000) : int.MaxValue;
        int samples = rng.Next(1, 5);
        double[] ccf = new double[samples];
        for (int i = 0; i < samples; i++)
        {
            // Include the boundary values 0 and 1 explicitly.
            int pick = rng.Next(3);
            ccf[i] = pick == 0 ? 0.0 : pick == 1 ? 1.0 : rng.NextDouble();
        }

        var clusters = new[] { new CcfCluster(id, ccf) };

        ClonalPhylogeny tree = ReconstructPhylogeny(clusters);

        AssertWellFormedTree(tree, clusters, DefaultPhylogenyTolerance);
        tree.ParentOf(id).Should().Be(tree.RootId, "the single clone always hangs off the root regardless of its CCF (§6.1)");
        tree.RootId.Should().NotBe(id, "the synthetic root id must stay distinct from the clone id (§4.1)");
        tree.RootId.Should().BeLessThanOrEqualTo(-1, "the synthetic root id is conventionally negative (≤ -1) (§4.1)");
    }

    // ── BE: identical clones ─────────────────────────────────────────────────
    // Two+ clusters with the SAME CCF vector. The sum rule (Eq. 5) forbids two
    // equal-CCF siblings under one parent (their CCFs would sum to 2×, exceeding the
    // parent), so they CHAIN deterministically. No infinite loop, no duplicate
    // parent, no self-ancestry (§6.1, INV-02, INV-05).

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_TwoIdenticalClones_SumRuleForcesChainNotSiblings()
    {
        // Founder at 1.0 then two equal 0.6 subclones (single sample).
        var clusters = new[] { Cluster(1, 1.0), Cluster(2, 0.6), Cluster(3, 0.6) };

        ClonalPhylogeny tree = ReconstructPhylogeny(clusters);

        AssertWellFormedTree(tree, clusters, DefaultPhylogenyTolerance);
        // 0.6 + 0.6 = 1.2 > 1.0 ⇒ the two equal subclones cannot both sit under the founder.
        tree.ParentOf(2).Should().Be(1, "the first 0.6 subclone nests under the founder");
        tree.ParentOf(3).Should().Be(2, "the second equal subclone chains below the first (sum rule, §6.1)");
        tree.ChildrenOf(1).Should().HaveCount(1, "the founder admits only one of two equal-CCF children (INV-02)");
    }

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_ManyIdenticalClones_FormChainNeverCycle([Range(0, 19)] int seed)
    {
        var rng = new Random(seed);
        int copies = rng.Next(2, 8);
        double shared = 0.2 + (rng.NextDouble() * 0.6); // a single non-trivial CCF, well below 1
        var list = new List<CcfCluster> { Cluster(0, 1.0) }; // a founder so the chain has somewhere to start
        for (int k = 1; k <= copies; k++)
        {
            list.Add(Cluster(k, shared));
        }

        ClonalPhylogeny tree = ReconstructPhylogeny(list);

        // The whole point: identical clones must NOT crash, hang, or corrupt the tree.
        AssertWellFormedTree(tree, list, DefaultPhylogenyTolerance);

        // With identical CCFs the sum rule allows at most one such child per parent
        // (2×shared > shared), so the identical clones form a single chain: each node
        // (after the founder) has at most one identical-clone child.
        for (int k = 1; k <= copies; k++)
        {
            int identicalChildren = tree.ChildrenOf(k).Count(c => c >= 1);
            identicalChildren.Should().BeLessThanOrEqualTo(1,
                "the sum rule admits at most one equal-CCF child per node (no sibling explosion) (INV-02)");
        }
    }

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_AllClonesIdenticalIncludingFounder_NoSelfAncestry([Range(0, 19)] int seed)
    {
        var rng = new Random(seed);
        int copies = rng.Next(2, 6);
        int samples = rng.Next(1, 4);
        double[] ccf = new double[samples];
        for (int i = 0; i < samples; i++)
        {
            ccf[i] = 0.3 + (rng.NextDouble() * 0.5);
        }

        var list = new List<CcfCluster>();
        for (int k = 0; k < copies; k++)
        {
            list.Add(new CcfCluster(k, (double[])ccf.Clone()));
        }

        ClonalPhylogeny tree = ReconstructPhylogeny(list);

        // No founder dominates (all equal) yet the result must still be a valid spanning
        // tree: a single deepest chain off the root, no clone its own ancestor, no cycle.
        AssertWellFormedTree(tree, list, DefaultPhylogenyTolerance);
        tree.Edges.Should().HaveCount(copies, "exactly one edge per identical clone (spanning tree)");
    }

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_IdenticalClones_Deterministic([Range(0, 9)] int seed)
    {
        var rng = new Random(seed);
        int copies = rng.Next(2, 7);
        double shared = 0.25 + (rng.NextDouble() * 0.5);
        var list = new List<CcfCluster> { Cluster(0, 1.0) };
        for (int k = 1; k <= copies; k++)
        {
            list.Add(Cluster(k, shared));
        }

        ClonalPhylogeny a = ReconstructPhylogeny(list);
        ClonalPhylogeny b = ReconstructPhylogeny(list);

        a.Edges.Should().Equal(b.Edges, "identical input ⇒ identical edge set (INV-05)");
    }

    // ── BE: no shared mutations (disjoint presence) ──────────────────────────
    // Clusters that are each private to a DIFFERENT sample (disjoint presence
    // patterns) share no mutations. The presence pattern forbids one descending from
    // another, so they form independent lineages off the trunk — no FALSE ancestry
    // inferred between disjoint clones (§2.2, §7.1).

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_DisjointPrivateClones_FormIndependentLineagesOffTrunk()
    {
        // Trunk present in both samples; two clones each private to one sample (disjoint).
        var clusters = new[]
        {
            Cluster(1, 1.0, 1.0), // trunk
            Cluster(2, 0.6, 0.0), // private to sample 0
            Cluster(3, 0.0, 0.7), // private to sample 1
        };

        ClonalPhylogeny tree = ReconstructPhylogeny(clusters);

        AssertWellFormedTree(tree, clusters, DefaultPhylogenyTolerance);
        tree.ParentOf(2).Should().Be(1, "clone 2 (private to s0) cannot descend from clone 3 (s0=0) ⇒ attaches to trunk");
        tree.ParentOf(3).Should().Be(1, "clone 3 (private to s1) cannot descend from clone 2 (s1=0) ⇒ attaches to trunk");
        tree.ChildrenOf(2).Should().NotContain(3, "no false ancestry between disjoint clones (presence pattern)");
        tree.ChildrenOf(3).Should().NotContain(2, "no false ancestry between disjoint clones (presence pattern)");
    }

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_DisjointClonesNoTrunk_AllHangOffRoot()
    {
        // No clone present in every sample: each is private to one sample (no shared founder).
        var clusters = new[]
        {
            Cluster(10, 0.8, 0.0, 0.0),
            Cluster(20, 0.0, 0.7, 0.0),
            Cluster(30, 0.0, 0.0, 0.9),
        };

        ClonalPhylogeny tree = ReconstructPhylogeny(clusters);

        AssertWellFormedTree(tree, clusters, DefaultPhylogenyTolerance);
        // None can parent another (every pair is disjoint), so the root is the only valid
        // ancestor for all three ⇒ a star/forest off the root, no false nesting.
        tree.ParentOf(10).Should().Be(tree.RootId, "disjoint clone attaches to the root (no valid cluster ancestor)");
        tree.ParentOf(20).Should().Be(tree.RootId, "disjoint clone attaches to the root (no valid cluster ancestor)");
        tree.ParentOf(30).Should().Be(tree.RootId, "disjoint clone attaches to the root (no valid cluster ancestor)");
    }

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_RandomDisjointPrivateClones_NoFalseAncestry([Range(0, 19)] int seed)
    {
        var rng = new Random(seed);
        int samples = rng.Next(2, 6);
        var list = new List<CcfCluster>();
        // One private clone per sample, each present ONLY in its own sample (disjoint).
        for (int s = 0; s < samples; s++)
        {
            double[] ccf = new double[samples];
            ccf[s] = 0.2 + (rng.NextDouble() * 0.7);
            list.Add(new CcfCluster(s + 1, ccf));
        }

        ClonalPhylogeny tree = ReconstructPhylogeny(list);

        AssertWellFormedTree(tree, list, DefaultPhylogenyTolerance);
        // Disjoint clones cannot be ancestor/descendant of each other ⇒ none is a child of
        // another cluster; every clone hangs directly off the root.
        foreach (CcfCluster c in list)
        {
            tree.ParentOf(c.Id).Should().Be(tree.RootId,
                $"clone {c.Id} is private to its own sample and shares no presence with any other ⇒ root child (no false ancestry)");
        }
    }

    // ── BE: empty cohort boundary ────────────────────────────────────────────

    [Test]
    [CancelAfter(HangGuardMs)]
    public void Be_EmptyCohort_ReturnsRootOnlyTree()
    {
        ClonalPhylogeny tree = ReconstructPhylogeny(Array.Empty<CcfCluster>());

        tree.Clusters.Should().BeEmpty("an empty cohort has no nodes (§6.1)");
        tree.Edges.Should().BeEmpty("a root-only tree has no edges (§6.1)");
        IdentifyTrunkMutations(tree).Should().BeEmpty("no clusters ⇒ no trunk (§6.1)");
        IdentifyBranchMutations(tree).Should().BeEmpty("no clusters ⇒ no branches (§6.1)");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region ONCO-PHYLO-001 — documented validation faults (BE: malformed input)
    // ═════════════════════════════════════════════════════════════════════════

    [Test]
    public void Validation_NullClusters_Throws()
    {
        Action act = () => ReconstructPhylogeny(null!);
        act.Should().Throw<ArgumentNullException>("null clusters is a documented contract violation (§3.3)");
    }

    [Test]
    public void Validation_NullCcfList_Throws()
    {
        var clusters = new[] { new CcfCluster(1, null!) };
        Action act = () => ReconstructPhylogeny(clusters);
        act.Should().Throw<ArgumentNullException>("a null per-sample CCF list is rejected (§3.3)");
    }

    [Test]
    public void Validation_EmptyCcfList_Throws()
    {
        var clusters = new[] { new CcfCluster(1, Array.Empty<double>()) };
        Action act = () => ReconstructPhylogeny(clusters);
        act.Should().Throw<ArgumentException>("an empty CCF list (zero samples) is rejected (§3.3)");
    }

    [Test]
    public void Validation_RaggedCcfLists_Throws()
    {
        var clusters = new[] { Cluster(1, 1.0, 1.0), Cluster(2, 0.5) };
        Action act = () => ReconstructPhylogeny(clusters);
        act.Should().Throw<ArgumentException>("clusters must share one sample count (§3.3)");
    }

    [Test]
    public void Validation_DuplicateIds_Throws()
    {
        var clusters = new[] { Cluster(1, 1.0), Cluster(1, 0.5) };
        Action act = () => ReconstructPhylogeny(clusters);
        act.Should().Throw<ArgumentException>("cluster ids must be unique (§3.3)");
    }

    [Test]
    public void Validation_NaNCcf_Throws()
    {
        var clusters = new[] { Cluster(1, double.NaN) };
        Action act = () => ReconstructPhylogeny(clusters);
        act.Should().Throw<ArgumentException>("a NaN CCF is rejected, not silently handled (§3.3)");
    }

    [Test]
    public void Validation_CcfAboveOne_Throws()
    {
        var clusters = new[] { Cluster(1, 1.0001) };
        Action act = () => ReconstructPhylogeny(clusters);
        act.Should().Throw<ArgumentException>("a CCF above 1 is rejected, not clamped (§3.3)");
    }

    [Test]
    public void Validation_NegativeCcf_Throws()
    {
        var clusters = new[] { Cluster(1, -0.0001) };
        Action act = () => ReconstructPhylogeny(clusters);
        act.Should().Throw<ArgumentException>("a negative CCF is rejected (§3.3)");
    }

    [Test]
    public void Validation_NegativeTolerance_Throws()
    {
        var clusters = new[] { Cluster(1, 1.0) };
        Action act = () => ReconstructPhylogeny(clusters, tolerance: -0.01);
        act.Should().Throw<ArgumentOutOfRangeException>("ε must be non-negative (§3.3)");
    }

    [Test]
    public void Validation_NaNTolerance_Throws()
    {
        var clusters = new[] { Cluster(1, 1.0) };
        Action act = () => ReconstructPhylogeny(clusters, tolerance: double.NaN);
        act.Should().Throw<ArgumentOutOfRangeException>("ε must be a number (§3.3)");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region ONCO-PHYLO-001 — positive sanity (nesting & branching reconstruct correctly)
    // ═════════════════════════════════════════════════════════════════════════

    // A clear nested clonal structure (founder ⊃ subclone1 ⊃ subclone2 by CCF) must
    // reconstruct the correct linear parent/child chain.
    [Test]
    [CancelAfter(HangGuardMs)]
    public void Positive_NestedClonalStructure_ReconstructsLinearChain()
    {
        var clusters = new[] { Cluster(1, 1.0), Cluster(2, 0.6), Cluster(3, 0.3) };

        ClonalPhylogeny tree = ReconstructPhylogeny(clusters);

        AssertWellFormedTree(tree, clusters, DefaultPhylogenyTolerance);
        tree.ParentOf(1).Should().Be(tree.RootId, "founder (1.0) attaches to the normal root");
        tree.ParentOf(2).Should().Be(1, "subclone1 (0.6) nests under the founder (deepest valid ancestor)");
        tree.ParentOf(3).Should().Be(2, "subclone2 (0.3) nests under subclone1 (deepest valid ancestor)");
        IdentifyTrunkMutations(tree).Should().Equal(new[] { 1, 2, 3 }, "a pure chain is all trunk, no branch point");
        IdentifyBranchMutations(tree).Should().BeEmpty("no divergence ⇒ no branch clusters");
    }

    // A branching structure (founder with two divergent, sample-private subclones)
    // must reconstruct the correct branching with both subclones as siblings.
    [Test]
    [CancelAfter(HangGuardMs)]
    public void Positive_BranchingStructure_ReconstructsTwoSiblings()
    {
        var clusters = new[]
        {
            Cluster(1, 1.0, 1.0), // founder/trunk
            Cluster(2, 0.6, 0.0), // divergent subclone, sample 0
            Cluster(3, 0.0, 0.7), // divergent subclone, sample 1
        };

        ClonalPhylogeny tree = ReconstructPhylogeny(clusters);

        AssertWellFormedTree(tree, clusters, DefaultPhylogenyTolerance);
        tree.ParentOf(1).Should().Be(tree.RootId, "founder attaches to the root");
        tree.ChildrenOf(1).OrderBy(x => x).Should().Equal(new[] { 2, 3 }, "the founder has two divergent sibling subclones");
        IdentifyTrunkMutations(tree).Should().Equal(new[] { 1 }, "the trunk ends at the branch point (the founder)");
        IdentifyBranchMutations(tree).OrderBy(x => x).Should().Equal(new[] { 2, 3 }, "both divergent subclones are subclonal branches");
    }

    // Tolerance ε must admit an edge whose violation is within the margin, and the
    // stricter ε = 0 must reject the same noisy edge (documented behaviour, §4.2).
    [Test]
    [CancelAfter(HangGuardMs)]
    public void Positive_ToleranceAdmitsNearViolation_StrictRejects()
    {
        // A (total 1.0) is processed first; in sample 0 the descendant B (0.55) slightly
        // exceeds ancestor A (0.50) — a noise overshoot. With ε = 0.1, A is admissible
        // (0.50 ≥ 0.55 − 0.10) and deeper than the root.
        var clusters = new[] { Cluster(1, 0.50, 0.50), Cluster(2, 0.55, 0.0) };

        ClonalPhylogeny lenient = ReconstructPhylogeny(clusters, tolerance: 0.1);
        AssertWellFormedTree(lenient, clusters, 0.1);
        lenient.ParentOf(2).Should().Be(1,
            "with ε = 0.1 the 0.05 overshoot is within margin ⇒ clone 2 nests under clone 1");

        ClonalPhylogeny strict = ReconstructPhylogeny(clusters, tolerance: 0.0);
        AssertWellFormedTree(strict, clusters, 0.0);
        strict.ParentOf(2).Should().Be(strict.RootId,
            "with ε = 0 the 0.05 overshoot violates lineage precedence ⇒ clone 2 attaches to the root, not clone 1");
    }

    #endregion
}
