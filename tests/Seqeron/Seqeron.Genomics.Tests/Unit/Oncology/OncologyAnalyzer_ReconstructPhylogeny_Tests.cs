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
// They are NOT copied from the implementation output. Exact topologies / error scores / tree counts are the output of
// the original LICHeE code (github.com/viq854/lichee, LICHeE/release/lichee.jar: PHYNetwork + getLineageTrees +
// evaluateLineageTrees, cell-prevalence mode VAF_MAX = 1, VAF_ERROR_MARGIN = e) run on the same clusters (B24 F18).

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_ReconstructPhylogeny_Tests
{
    private static OncologyAnalyzer.CcfCluster C(int id, params double[] ccf) =>
        new(id, ccf);

    private static (int parent, int child)[] EdgeTuples(OncologyAnalyzer.ClonalPhylogeny p) =>
        p.Edges.Select(e => (e.ParentId, e.ChildId)).OrderBy(t => t.Item1).ThenBy(t => t.Item2).ToArray();

    #region ReconstructPhylogeny

    // M1 — single sample A=1.0,B=0.6,C=0.3. B can only hang under A (root: 1.0+0.6 > 1, Eq.5); C fits under A
    // (0.6+0.3 <= 1) or under B, so exactly two trees are valid. LICHeE (lichee.jar) enumerates both and ranks
    // A->{B,C} first (error score 0 for both; enumeration order breaks the tie).
    [Test]
    public void ReconstructPhylogeny_DescendingSingleSampleCcf_MatchesLicheeTopTree()
    {
        var clusters = new[] { C(1, 1.0), C(2, 0.6), C(3, 0.3) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        int root = p.RootId;
        Assert.Multiple(() =>
        {
            Assert.That(p.ParentOf(1), Is.EqualTo(root), "A (CCF 1.0) attaches to the normal root");
            Assert.That(p.ParentOf(2), Is.EqualTo(1), "B (0.6) cannot share the root with A (1.6 > 1, Eq.5)");
            Assert.That(p.ParentOf(3), Is.EqualTo(1), "LICHeE top tree places C under A (0.6+0.3 <= 1.0)");
            Assert.That(p.ValidTreeCount, Is.EqualTo(2), "LICHeE: 'Found 2 valid tree(s)'");
            Assert.That(p.ErrorScore, Is.EqualTo(0.0), "every sum rule holds exactly");
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
    // B and C cannot both be children of A (0.6+0.6=1.2 > 1.0, Eq.5) => one nests under the other (0.6<=0.6).
    // LICHeE checkAndAddEdge orients the equal-CCF pair later->earlier (err12 < err21 is false), so the single valid
    // tree is root->A->C->B (lichee.jar). Popic (2015) Eq.5.
    [Test]
    public void ReconstructPhylogeny_TwoEqualSubclones_SumRuleForcesChain()
    {
        var clusters = new[] { C(1, 1.0), C(2, 0.6), C(3, 0.6) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        int root = p.RootId;
        Assert.Multiple(() =>
        {
            Assert.That(p.ParentOf(1), Is.EqualTo(root), "A attaches to root");
            Assert.That(p.ParentOf(3), Is.EqualTo(1), "one 0.6 subclone nests under A");
            Assert.That(p.ParentOf(2), Is.EqualTo(3),
                "the other 0.6 subclone cannot also be A's child (sum 1.2 > 1.0); it chains below (Eq.5)");
            Assert.That(p.ValidTreeCount, Is.EqualTo(1), "LICHeE finds exactly one valid tree");
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

        Assert.Multiple(() =>
        {
            Assert.That(p.ParentOf(2), Is.EqualTo(1),
                "with e=0.1, B.s1 (0.55) is admissible under A.s1 (0.50 >= 0.55-0.10)");
            Assert.That(p.ErrorScore, Is.EqualTo(0.050000000000000044).Within(1e-15),
                "LICHeE error score sqrt((0.55-0.50)^2) (lichee.jar)");
        });
    }

    // S3 — Strict (e=0) rejects the noisy edge. Same input; A.s1 (0.50) < B.s1 (0.55) violates Eq.2 under A, so the
    // root is B's only admissible parent — but then the root's children sum to 0.50+0.55 = 1.05 > 1.0 in sample 1
    // (Eq.5). No valid tree exists; LICHeE (lichee.jar, also with the complete network) reports none. The former greedy
    // returned root->B anyway, violating Eq.5 (B24 F18). ε = 0 is explicit since the default became LICHeE's 0.1
    // (B24 F44); at the default lichee.jar (-e 0.1) admits A->B (0.55 <= 0.5 + 0.1): 1 tree, error 0.050000000000000044.
    [Test]
    public void ReconstructPhylogeny_StrictRejectsNoisyEdge_NoValidTree()
    {
        var clusters = new[] { C(1, 0.50, 0.50), C(2, 0.55, 0.0) };

        OncologyAnalyzer.ClonalPhylogeny atDefault = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.TryReconstructPhylogeny(clusters, out _, tolerance: 0.0), Is.False,
                "no spanning tree satisfies Eq.5 at e=0");
            Assert.Throws<InvalidOperationException>(() => OncologyAnalyzer.ReconstructPhylogeny(clusters, tolerance: 0.0),
                "ReconstructPhylogeny reports the absence of a valid tree instead of returning an invalid one");
            Assert.That(atDefault.ParentOf(1), Is.EqualTo(atDefault.RootId), "default e = 0.1 (lichee.jar -e 0.1)");
            Assert.That(atDefault.ParentOf(2), Is.EqualTo(1));
            Assert.That(atDefault.ValidTreeCount, Is.EqualTo(1));
            Assert.That(atDefault.ErrorScore, Is.EqualTo(0.050000000000000044));
        });
    }

    // F18 repro — e=0, A=[0.3,1.0], B=[0.2,0.2], C=[0,0.8]. Valid tree: A->{B,C} (s2: 0.8+0.2 = 1.0 <= 1.0).
    // The former greedy debited A's budget to 1.0-0.8 = 0.19999999999999996 < 0.2 and fell back to root->B
    // (root s2: 1.0+0.2 = 1.2 > 1, Eq.5 violated). LICHeE (lichee.jar): "default 1 0.0 A:root B:A C:A".
    [Test]
    public void ReconstructPhylogeny_ExactSumRuleEquality_MatchesLichee()
    {
        var clusters = new[] { C(1, 0.3, 1.0), C(2, 0.2, 0.2), C(3, 0.0, 0.8) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        Assert.Multiple(() =>
        {
            Assert.That(p.ParentOf(1), Is.EqualTo(p.RootId));
            Assert.That(p.ParentOf(2), Is.EqualTo(1), "B fits under A: 0.8+0.2 = 1.0 <= A.s2 = 1.0 (Eq.5)");
            Assert.That(p.ParentOf(3), Is.EqualTo(1));
            Assert.That(p.ValidTreeCount, Is.EqualTo(1));
            Assert.That(p.ErrorScore, Is.EqualTo(0.0));
        });
    }

    // Complete-network fallback — e=0.02, 3 samples: A=[0.482,0,0.443137], B=[0.519,0.796779,0.56],
    // C=[0.19,0.418589,0.24]. The default (level-adjacent) network admits no valid tree; with ALL_EDGES LICHeE finds one:
    // A,B under the root, C under B, error score sqrt((0.482+0.519-1)^2 + (0.443137+0.56-1)^2) = 0.003292532308117998
    // (lichee.jar output "complete 1 0.003292532308117998").
    [Test]
    public void ReconstructPhylogeny_DefaultNetworkInfeasible_UsesCompleteNetworkLikeLichee()
    {
        var clusters = new[]
        {
            C(1, 0.482, 0.0, 0.443137),
            C(2, 0.519, 0.796779, 0.56),
            C(3, 0.19, 0.418589, 0.24),
        };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters, tolerance: 0.02);

        Assert.Multiple(() =>
        {
            Assert.That(p.UsedCompleteNetwork, Is.True);
            Assert.That(p.ParentOf(1), Is.EqualTo(p.RootId));
            Assert.That(p.ParentOf(2), Is.EqualTo(p.RootId));
            Assert.That(p.ParentOf(3), Is.EqualTo(2));
            Assert.That(p.ValidTreeCount, Is.EqualTo(1));
            Assert.That(p.ErrorScore, Is.EqualTo(0.003292532308117998).Within(1e-15));
        });
    }

    // Enumeration — 4 private single-sample clusters 0.05..0.053: at the default e = 0.1 lichee.jar (-e 0.1) finds 24
    // valid trees (15 at e = 0) and ranks the chain root->D->C->B->A first in both cases (error 0; B24 F44 re-derivation).
    [Test]
    public void ReconstructPhylogeny_ManyValidTrees_CountAndTopTreeMatchLichee()
    {
        var clusters = new[] { C(1, 0.05), C(2, 0.051), C(3, 0.052), C(4, 0.053) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);
        OncologyAnalyzer.ClonalPhylogeny strict = OncologyAnalyzer.ReconstructPhylogeny(clusters, tolerance: 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(p.ValidTreeCount, Is.EqualTo(24));
            Assert.That(strict.ValidTreeCount, Is.EqualTo(15));
            Assert.That(EdgeTuples(strict), Is.EqualTo(EdgeTuples(p)));
            Assert.That(p.ErrorScore, Is.EqualTo(0.0));
            Assert.That(p.ParentOf(4), Is.EqualTo(p.RootId));
            Assert.That(p.ParentOf(3), Is.EqualTo(4));
            Assert.That(p.ParentOf(2), Is.EqualTo(3));
            Assert.That(p.ParentOf(1), Is.EqualTo(2));
        });
    }

    // F44 (A21) — the default ε is LICHeE's -e default 0.1, also in cell-prevalence mode (-cp only sets VAF_MAX =
    // MAX_ALLOWED_VAF = 1; LineageEngine option parsing, Parameters.VAF_ERROR_MARGIN = 0.1).
    [Test]
    public void DefaultPhylogenyTolerance_IsLicheeDefault()
    {
        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(new[] { C(1, 0.5, 0.5), C(2, 0.55, 0.0) });

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.DefaultPhylogenyTolerance, Is.EqualTo(0.1));
            Assert.That(p.Tolerance, Is.EqualTo(0.1));
        });
    }

    // F44 (A22) — grow-call cap (LICHeE Parameters.MAX_NUM_GROW_CALLS, here reduced through the internal test hook).
    // 7 clusters with presence profile 11, e = 0.05. lichee.jar driven by the harness with MAX_NUM_GROW_CALLS patched:
    //   cap 40  -> 0 trees (default and ALL_EDGES networks both capped, fixNetwork removes nothing) -> no tree;
    //   cap 80  -> 7 trees, error 0.07071067811865477, root->{1,7}, 1->{2,4}, 2->3, 4->5, 7->6;
    //   cap 320 -> 17 trees, error 0.050000000000000044, root->{1,7}, 1->{2,4}, 2->3, 4->5, 5->6;
    //   10^8    -> 176 trees, error 0, root->{1,6}, 1->{2,4}, 2->3, 4->5, 5->7.
    // LICHeE's frame returns at the cap while its callers keep looping, but no further recursion (hence no tree) can
    // happen, so the result is the best of the trees completed within the budget — as the port's immediate stop.
    [TestCase(40, 0, double.NaN, new int[0])]
    [TestCase(80, 7, 0.07071067811865477, new[] { -1, 1, 2, 1, 4, 7, -1 })]
    [TestCase(320, 17, 0.050000000000000044, new[] { -1, 1, 2, 1, 4, 5, -1 })]
    [TestCase(OncologyAnalyzer.MaxPhylogenyGrowCalls, 176, 0.0, new[] { -1, 1, 2, 1, 4, -1, 5 })]
    public void GrowCallCap_ReturnsBestTreeWithinBudgetLikeLichee(int cap, int trees, double error, int[] parents)
    {
        var clusters = new[]
        {
            C(1, 0.9, 0.8), C(2, 0.5, 0.45), C(3, 0.4, 0.3), C(4, 0.3, 0.35),
            C(5, 0.2, 0.25), C(6, 0.1, 0.15), C(7, 0.15, 0.1),
        };

        bool found = OncologyAnalyzer.TryReconstructPhylogeny(clusters, out OncologyAnalyzer.ClonalPhylogeny p, 0.05, cap);

        if (trees == 0)
        {
            Assert.That(found, Is.False, "lichee.jar: no tree within the budget (default and complete network)");
            return;
        }

        // parents[k] = parent of cluster k + 1 (-1 = root), from the jar's edge list.
        Assert.Multiple(() =>
        {
            Assert.That(found, Is.True);
            Assert.That(p.ValidTreeCount, Is.EqualTo(trees));
            Assert.That(p.ErrorScore, Is.EqualTo(error));
            Assert.That(Enumerable.Range(1, 7).Select(id => p.ParentOf(id) == p.RootId ? -1 : p.ParentOf(id)!.Value),
                Is.EqualTo(parents));
        });
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

    // Trunk = alterations present in all tumour cells (Werner et al. 2017 Sci Rep 7:44991: "alterations that are in the
    // trunk of the tree must be present in all cells of the tumour"). A single-child chain root->A(1.0)->B(0.5)->C(0.25)
    // keeps B and C subclonal: A-only cells (0.5) lack B. The former structural rule reported {A,B,C} (B24 F19).
    [Test]
    public void IdentifyTrunkAndBranch_SubclonalChain_OnlyClonalNodeIsTrunk()
    {
        var clusters = new[] { C(1, 1.0), C(2, 0.5), C(3, 0.25) };
        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogeny(clusters);

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IdentifyTrunkMutations(p), Is.EqualTo(new[] { 1 }),
                "only the CCF-1 cluster is present in every tumour cell");
            Assert.That(OncologyAnalyzer.IdentifyBranchMutations(p), Is.EqualTo(new[] { 2, 3 }),
                "CCF 0.5 / 0.25 clusters are subclonal");
        });
    }

    // Two clonal clusters (CCF 1 in every sample) form a two-node trunk; the trunk criterion uses the tree's e.
    [Test]
    public void IdentifyTrunkMutations_ClonalChainAndTolerance_FollowCcfCriterion()
    {
        var clonal = OncologyAnalyzer.ReconstructPhylogeny(new[] { C(1, 1.0, 1.0), C(2, 1.0, 1.0), C(3, 0.4, 0.0) });
        var nearClonalStrict = OncologyAnalyzer.ReconstructPhylogeny(new[] { C(1, 0.97), C(2, 0.5) }, tolerance: 0.0);
        var nearClonalNoisy = OncologyAnalyzer.ReconstructPhylogeny(new[] { C(1, 0.97), C(2, 0.5) }, tolerance: 0.05);

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.IdentifyTrunkMutations(clonal), Has.Count.EqualTo(2)
                .And.EquivalentTo(new[] { 1, 2 }), "both CCF-1 clusters are truncal");
            Assert.That(OncologyAnalyzer.IdentifyTrunkMutations(nearClonalStrict), Is.Empty,
                "with e=0 a CCF of 0.97 is subclonal");
            Assert.That(OncologyAnalyzer.IdentifyTrunkMutations(nearClonalNoisy), Is.EqualTo(new[] { 1 }),
                "with e=0.05, 0.97 >= 1-0.05 counts as clonal");
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
