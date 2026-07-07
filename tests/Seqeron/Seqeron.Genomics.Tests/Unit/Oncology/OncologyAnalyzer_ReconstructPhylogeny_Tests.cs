// ONCO-PHYLO-001 — Tumor Phylogeny Reconstruction (clonal tree from CCF clusters)
// Evidence: docs/Evidence/ONCO-PHYLO-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-PHYLO-001.md
// Source: Popic V et al. (2015). Fast and scalable inference of multi-sample cancer lineages.
//         Genome Biology 16:91. https://doi.org/10.1186/s13059-015-0647-8
//         Zheng L et al. (2022). PICTograph. Bioinformatics 38(15):3677-3683.
//         https://doi.org/10.1093/bioinformatics/btac367
//
// Expected edges/relationships below are derived independently from the cited rules:
//   Lineage precedence (Eq.2): parent.CCF[i] >= child.CCF[i]-e and parent=0 => child=0 per sample.
//   Sum rule (Eq.5): sum over children of v.CCF[i] <= u.CCF[i]+e per node, per sample.
// They are NOT copied from the implementation output.

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_ReconstructPhylogeny_Tests
{
    private static OncologyAnalyzer.CcfCluster C(int id, params double[] ccf) =>
        new(id, ccf);

    private static (int parent, int child)[] EdgeTuples(OncologyAnalyzer.ClonalPhylogeny p) =>
        p.Edges.Select(e => (e.ParentId, e.ChildId)).OrderBy(t => t.Item1).ThenBy(t => t.Item2).ToArray();

    #region ReconstructPhylogeny

    // M1 — Linear chain: single sample A=1.0,B=0.6,C=0.3. Eq.2 (anc>=desc): each nests in the previous.
    // Deepest-valid-ancestor => root->A->B->C. Popic (2015) Eq.2.
    [Test]
    public void ReconstructPhylogeny_DescendingSingleSampleCcf_FormsLinearChain()
    {
        var clusters = new[] { C(1, 1.0), C(2, 0.6), C(3, 0.3) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        int root = p.RootId;
        Assert.Multiple(() =>
        {
            Assert.That(p.ParentOf(1), Is.EqualTo(root), "A (CCF 1.0) attaches to the normal root (1.0 >= 1.0)");
            Assert.That(p.ParentOf(2), Is.EqualTo(1), "B (0.6) nests under A: deepest valid ancestor, 1.0 >= 0.6");
            Assert.That(p.ParentOf(3), Is.EqualTo(2), "C (0.3) nests under B: deepest valid ancestor, 0.6 >= 0.3");
        });
    }

    // M2 — Branching: 2 samples A=[1,1] trunk, B=[0.6,0] private s1, C=[0,0.7] private s2.
    // B,C cannot be ancestor/descendant of each other (presence pattern, constraint 1) => both children of A.
    // Popic (2015) constraints (1)(2)(3).
    [Test]
    public void ReconstructPhylogeny_TwoPrivateSubclones_FormsBranchingTree()
    {
        var clusters = new[] { C(1, 1.0, 1.0), C(2, 0.6, 0.0), C(3, 0.0, 0.7) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        int root = p.RootId;
        Assert.Multiple(() =>
        {
            Assert.That(p.ParentOf(1), Is.EqualTo(root), "A (CCF 1 in both samples) is the trunk attached to root");
            Assert.That(p.ParentOf(2), Is.EqualTo(1), "B private to s1 attaches to A; cannot descend from C (C.s1=0 < B.s1)");
            Assert.That(p.ParentOf(3), Is.EqualTo(1), "C private to s2 attaches to A; cannot descend from B (B.s2=0 < C.s2)");
            Assert.That(p.ChildrenOf(1).OrderBy(x => x), Is.EqualTo(new[] { 2, 3 }), "A has the two sibling branches B and C");
        });
    }

    // M3 — Sum rule forces a chain: single sample A=1.0,B=0.6,C=0.6.
    // B and C cannot both be children of A (0.6+0.6=1.2 > 1.0, Eq.5) => C nests under B (0.6<=0.6).
    // Popic (2015) Eq.5.
    [Test]
    public void ReconstructPhylogeny_TwoEqualSubclones_SumRuleForcesChain()
    {
        var clusters = new[] { C(1, 1.0), C(2, 0.6), C(3, 0.6) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        int root = p.RootId;
        Assert.Multiple(() =>
        {
            Assert.That(p.ParentOf(1), Is.EqualTo(root), "A attaches to root");
            Assert.That(p.ParentOf(2), Is.EqualTo(1), "first 0.6 subclone nests under A");
            Assert.That(p.ParentOf(3), Is.EqualTo(2),
                "second 0.6 subclone cannot also be A's child (sum 1.2 > 1.0); it chains under B (Eq.5)");
        });
    }

    // M8 — Single cluster: A=1.0 attaches to root, is the only (trunk) node.
    [Test]
    public void ReconstructPhylogeny_SingleCluster_AttachesToRoot()
    {
        var clusters = new[] { C(7, 1.0) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        Assert.Multiple(() =>
        {
            Assert.That(p.Edges.Count, Is.EqualTo(1), "exactly one edge: root->cluster");
            Assert.That(p.ParentOf(7), Is.EqualTo(p.RootId), "the only cluster attaches to the synthetic root");
        });
    }

    // S1 — Empty input: root-only tree, no edges.
    [Test]
    public void ReconstructPhylogeny_EmptyInput_ReturnsRootOnlyTree()
    {
        OncologyAnalyzer.ClonalPhylogeny p =
            OncologyAnalyzer.ReconstructPhylogeny(Array.Empty<OncologyAnalyzer.CcfCluster>());

        Assert.Multiple(() =>
        {
            Assert.That(p.Edges, Is.Empty, "no clusters => no edges");
            Assert.That(p.Clusters, Is.Empty, "no clusters present");
        });
    }

    // S2 — Tolerance admits a noisy near-violation. A=[0.50,0.50] (total 1.0, processed first), B=[0.55,0.0].
    // In sample 1 the descendant B (0.55) slightly exceeds ancestor A (0.50) due to noise. With e=0.1, A is valid
    // (0.50 >= 0.55-0.10 = 0.45) and deeper than root => B nests under A. Popic 2015 Eq.2 (relaxed by e).
    [Test]
    public void ReconstructPhylogeny_ToleranceAdmitsNearViolation_NestsUnderParent()
    {
        var clusters = new[] { C(1, 0.50, 0.50), C(2, 0.55, 0.0) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters, tolerance: 0.1);

        Assert.That(p.ParentOf(2), Is.EqualTo(1),
            "with e=0.1, B.s1 (0.55) is admissible under A.s1 (0.50 >= 0.55-0.10); deepest valid ancestor is A, not root");
    }

    // S3 — Strict (e=0) rejects the noisy edge. Same input; A.s1 (0.50) >= B.s1 (0.55) is false, so A is not a
    // valid parent; the deepest valid candidate is the root (1.0 >= 0.55) => B attaches to root. Popic 2015 Eq.2 (strict).
    [Test]
    public void ReconstructPhylogeny_StrictRejectsNoisyEdge_AttachesToRoot()
    {
        var clusters = new[] { C(1, 0.50, 0.50), C(2, 0.55, 0.0) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        Assert.That(p.ParentOf(2), Is.EqualTo(p.RootId),
            "with e=0, B.s1 (0.55) > A.s1 (0.50) violates Eq.2 under A; the only valid parent is the root (1.0 >= 0.55)");
    }

    // C1 — Determinism: same input run twice yields identical edge set (INV-05).
    [Test]
    public void ReconstructPhylogeny_RunTwice_ProducesIdenticalEdges()
    {
        var clusters = new[] { C(1, 1.0, 1.0), C(2, 0.6, 0.0), C(3, 0.0, 0.7) };

        var first = EdgeTuples(OncologyAnalyzer.ReconstructPhylogeny(clusters));
        var second = EdgeTuples(OncologyAnalyzer.ReconstructPhylogeny(clusters));

        Assert.That(second, Is.EqualTo(first), "reconstruction is deterministic (INV-05)");
    }

    #endregion

    #region Invariants

    // M6 — INV-1: every edge has ancestor CCF >= descendant CCF per sample (Eq.2). Checked on M1 and M2 trees.
    [Test]
    public void ReconstructPhylogeny_EveryEdge_SatisfiesAncestorGeqDescendant()
    {
        var inputs = new[]
        {
            new[] { C(1, 1.0), C(2, 0.6), C(3, 0.3) },
            new[] { C(1, 1.0, 1.0), C(2, 0.6, 0.0), C(3, 0.0, 0.7) },
        };

        foreach (var clusters in inputs)
        {
            OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);
            var ccf = BuildCcfLookup(p);
            foreach (OncologyAnalyzer.ClonalEdge e in p.Edges)
            {
                double[] parent = ccf[e.ParentId];
                double[] child = ccf[e.ChildId];
                for (int i = 0; i < p.SampleCount; i++)
                {
                    Assert.That(parent[i], Is.GreaterThanOrEqualTo(child[i] - 1e-12),
                        $"INV-1: edge {e.ParentId}->{e.ChildId} sample {i}: ancestor CCF must be >= descendant CCF (Eq.2)");
                }
            }
        }
    }

    // M7 — INV-2 property: for random valid CCF inputs (fixed seed 42), per-node children CCF sum <= parent CCF (Eq.5).
    [Test]
    public void ReconstructPhylogeny_RandomValidInputs_SatisfySumRule()
    {
        var rng = new Random(42); // fixed documented seed for determinism
        for (int trial = 0; trial < 50; trial++)
        {
            int n = 2 + rng.Next(6);   // 2..7 clusters
            int k = 1 + rng.Next(3);   // 1..3 samples

            // Generate MODEL-CONSISTENT CCFs: per sample the cluster CCFs sum to <= 1, so a feasible tree
            // (at minimum a star under the root) always exists and the sum rule (Eq.5) is satisfiable. This is
            // the valid domain of INV-2; inputs that over-saturate a sample violate the perfect-phylogeny model.
            var perSample = new double[k][];
            for (int s = 0; s < k; s++)
            {
                double remaining = 1.0;
                var col = new double[n];
                for (int id = 0; id < n; id++)
                {
                    double share = Math.Round(rng.NextDouble() * remaining, 3);
                    col[id] = share;
                    remaining -= share;
                }

                perSample[s] = col;
            }

            var clusters = new List<OncologyAnalyzer.CcfCluster>(n);
            for (int id = 1; id <= n; id++)
            {
                var ccf = new double[k];
                for (int s = 0; s < k; s++)
                {
                    ccf[s] = perSample[s][id - 1];
                }

                clusters.Add(new OncologyAnalyzer.CcfCluster(id, ccf));
            }

            OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);
            var ccfLookup = BuildCcfLookup(p);

            // Verify sum rule at every node (root + each cluster).
            var nodes = new List<int> { p.RootId };
            nodes.AddRange(clusters.Select(c => c.Id));
            foreach (int node in nodes)
            {
                var children = p.ChildrenOf(node);
                for (int s = 0; s < p.SampleCount; s++)
                {
                    double sum = children.Sum(ch => ccfLookup[ch][s]);
                    Assert.That(sum, Is.LessThanOrEqualTo(ccfLookup[node][s] + 1e-9),
                        $"INV-2: trial {trial} node {node} sample {s}: children CCF sum must not exceed parent CCF (Eq.5)");
                }
            }
        }
    }

    private static Dictionary<int, double[]> BuildCcfLookup(OncologyAnalyzer.ClonalPhylogeny p)
    {
        double[] rootCcf = new double[p.SampleCount];
        Array.Fill(rootCcf, 1.0);
        var map = new Dictionary<int, double[]> { [p.RootId] = rootCcf };
        foreach (OncologyAnalyzer.CcfCluster c in p.Clusters)
        {
            map[c.Id] = c.CcfPerSample.ToArray();
        }

        return map;
    }

    #endregion

    #region Trunk / Branch identification

    // M4 + M5 — Trunk = {A}; Branches = {B, C} on the branching M2 tree. Popic (2015).
    [Test]
    public void IdentifyTrunkAndBranch_BranchingTree_TrunkIsCommonAncestor()
    {
        var clusters = new[] { C(1, 1.0, 1.0), C(2, 0.6, 0.0), C(3, 0.0, 0.7) };
        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        IReadOnlyList<int> trunk = OncologyAnalyzer.IdentifyTrunkMutations(p);
        IReadOnlyList<int> branches = OncologyAnalyzer.IdentifyBranchMutations(p);

        Assert.Multiple(() =>
        {
            Assert.That(trunk, Is.EqualTo(new[] { 1 }), "trunk = the single clonal ancestor A before the branch point");
            Assert.That(branches.OrderBy(x => x), Is.EqualTo(new[] { 2, 3 }), "branches = the two subclones B and C");
        });
    }

    // Trunk on a linear chain is the whole chain; no branches.
    [Test]
    public void IdentifyTrunkAndBranch_LinearChain_AllTrunkNoBranches()
    {
        var clusters = new[] { C(1, 1.0), C(2, 0.6), C(3, 0.3) };
        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IdentifyTrunkMutations(p), Is.EqualTo(new[] { 1, 2, 3 }),
                "a pure chain has every node on the trunk (single-child path from root)");
            Assert.That(OncologyAnalyzer.IdentifyBranchMutations(p), Is.Empty,
                "a pure chain has no subclonal branches");
        });
    }

    [Test]
    public void IdentifyTrunkAndBranch_EmptyTree_BothEmpty()
    {
        OncologyAnalyzer.ClonalPhylogeny p =
            OncologyAnalyzer.ReconstructPhylogeny(Array.Empty<OncologyAnalyzer.CcfCluster>());

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IdentifyTrunkMutations(p), Is.Empty, "no clusters => no trunk");
            Assert.That(OncologyAnalyzer.IdentifyBranchMutations(p), Is.Empty, "no clusters => no branches");
        });
    }

    #endregion

    #region Validation

    [Test]
    public void ReconstructPhylogeny_NullClusters_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => OncologyAnalyzer.ReconstructPhylogeny(null!),
            "null cluster list must throw ArgumentNullException");
    }

    [Test]
    public void ReconstructPhylogeny_CcfOutOfRange_Throws()
    {
        var clusters = new[] { C(1, 1.0), C(2, 1.5) };
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ReconstructPhylogeny(clusters),
            "CCF > 1 must throw ArgumentException");
    }

    [Test]
    public void ReconstructPhylogeny_NaNCcf_Throws()
    {
        var clusters = new[] { C(1, 1.0), C(2, double.NaN) };
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ReconstructPhylogeny(clusters),
            "NaN CCF must throw ArgumentException");
    }

    [Test]
    public void ReconstructPhylogeny_RaggedSampleCounts_Throws()
    {
        var clusters = new[] { C(1, 1.0, 1.0), C(2, 0.5) };
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ReconstructPhylogeny(clusters),
            "clusters with differing sample counts must throw ArgumentException");
    }

    [Test]
    public void ReconstructPhylogeny_DuplicateIds_Throws()
    {
        var clusters = new[] { C(1, 1.0), C(1, 0.5) };
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.ReconstructPhylogeny(clusters),
            "duplicate cluster ids must throw ArgumentException");
    }

    [Test]
    public void ReconstructPhylogeny_NegativeTolerance_Throws()
    {
        var clusters = new[] { C(1, 1.0) };
        Assert.Throws<ArgumentOutOfRangeException>(
            () => OncologyAnalyzer.ReconstructPhylogeny(clusters, tolerance: -0.1),
            "negative tolerance must throw ArgumentOutOfRangeException");
    }

    #endregion
}
