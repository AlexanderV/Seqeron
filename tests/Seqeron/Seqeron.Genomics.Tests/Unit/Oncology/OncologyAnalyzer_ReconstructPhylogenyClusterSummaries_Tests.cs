// ONCO-PHYLO-001 — LICHeE per-cluster edge error margins (B24 F42) and fixNetwork (B24 F43)
// Evidence: docs/Evidence/ONCO-PHYLO-001-Evidence.md (§ F42/F43)
// TestSpec: tests/TestSpecs/ONCO-PHYLO-001.md
// Source: github.com/viq854/lichee (master 26c2a70) LICHeE/src/lineage/PHYNetwork.java getAAFErrorMargin /
//         checkAndAddEdge / fixNetwork, LineageEngine.buildLineage steps 4–6, AAFClusterer.Cluster
//         recomputeCentroidAndStdDev, SNVGroup robustness (MIN_ROBUST_CLUSTER_SUPPORT = 2).
//
// Every expected topology / tree count / error score / centroid / SD below is the output of the original LICHeE
// classes (LICHeE/release/lichee.jar, OpenJDK 21) driven by a harness that builds SNVGroups from the member CCF rows,
// computes each cluster with Cluster.recomputeCentroidAndStdDev and runs LineageEngine.buildLineage steps 4–6
// verbatim (VAF_MAX = 1, VAF_ERROR_MARGIN = e). Not copied from the implementation. Fixtures derived at e = 0 pass
// tolerance 0 explicitly since the library default became LICHeE's 0.1 (B24 F44); DefaultTolerance_* locks the
// lichee.jar -e 0.1 results.

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_ReconstructPhylogenyClusterSummaries_Tests
{
    /// <summary>Builds a cluster summary from member CCF rows (LICHeE recomputeCentroidAndStdDev).</summary>
    private static OncologyAnalyzer.CcfClusterSummary M(int id, params double[][] members) =>
        OncologyAnalyzer.CcfClusterSummary.FromMembers(id, members);

    private static double[] R(params double[] v) => v;

    private static OncologyAnalyzer.CcfCluster Centroid(OncologyAnalyzer.CcfClusterSummary s) => new(s.Id, s.CcfPerSample);

    // p1 (4 samples, e = 0.05): A = members {0.44}x4 / {0.56}x4 (centroid 0.5, SD 0.06), B = 2 x [0.46,0.46,0.46,0.57].
    private static OncologyAnalyzer.CcfClusterSummary[] P1() => new[]
    {
        M(1, R(0.44, 0.44, 0.44, 0.44), R(0.56, 0.56, 0.56, 0.56)),
        M(2, R(0.46, 0.46, 0.46, 0.57), R(0.46, 0.46, 0.46, 0.57)),
    };

    // t02150 (2 samples, e = 0).
    private static OncologyAnalyzer.CcfClusterSummary[] T02150() => new[]
    {
        M(1, R(0.32, 0), R(0.292, 0), R(0.394, 0), R(0.556, 0), R(0.312, 0)),
        M(2, R(0.539, 0.447), R(0.505, 0.469), R(0.57, 0.422)),
        M(3, R(0, 0.409), R(0, 0.515)),
    };

    // t01061 (3 samples, e = 0.1).
    private static OncologyAnalyzer.CcfClusterSummary[] T01061() => new[]
    {
        M(1, R(0.786, 0.439, 0.478), R(0.63, 0.568, 0.532)),
        M(2, R(0.999, 0.95, 0.973), R(0.922, 0.87, 0.892)),
        M(3, R(0.323, 0.4, 0.444), R(1.0, 0.189, 0.691), R(0.769, 0.881, 0.825), R(1.0, 0.521, 0.789)),
        M(4, R(0.863, 0.503, 0.292), R(0.787, 0.387, 0.593)),
    };

    // t00298 (2 samples, e = 0).
    private static OncologyAnalyzer.CcfClusterSummary[] T00298() => new[]
    {
        M(1, R(0, 0.729), R(0, 0.41), R(0, 0.677), R(0, 0.041), R(0, 0.55)),
        M(2, R(0.336, 0), R(0.438, 0), R(0.445, 0)),
        M(3, R(0.303, 0), R(0.298, 0), R(0.298, 0), R(0.319, 0)),
        M(4, R(0.444, 0.289), R(0.468, 0.311), R(0.359, 0.355)),
        M(5, R(0, 0.204), R(0, 0.001)),
        M(6, R(0, 0.15), R(0, 0.342), R(0, 0.266), R(0, 0.325)),
    };

    #region F42 — FromMembers = LICHeE Cluster.recomputeCentroidAndStdDev

    // Centroid = sum / n, SD = sqrt(sum (x - mean)^2 / n) (population SD), bit-identical to the Java doubles printed by
    // the harness ("# cluster 2 s0 mean 0.5379999999999999 sd 0.02654555832275271", ...).
    [Test]
    public void FromMembers_CentroidAndPopulationSd_BitIdenticalToLichee()
    {
        OncologyAnalyzer.CcfClusterSummary[] t = T02150();
        OncologyAnalyzer.CcfClusterSummary[] p = P1();

        Assert.Multiple(() =>
        {
            Assert.That(t[0].CcfPerSample[0], Is.EqualTo(0.3748));
            Assert.That(t[0].StdDevPerSample[0], Is.EqualTo(0.09695854784391114));
            Assert.That(t[0].MemberCount, Is.EqualTo(5));
            Assert.That(t[1].CcfPerSample[0], Is.EqualTo(0.5379999999999999));
            Assert.That(t[1].StdDevPerSample[0], Is.EqualTo(0.02654555832275271));
            Assert.That(t[1].CcfPerSample[1], Is.EqualTo(0.44599999999999995));
            Assert.That(t[1].StdDevPerSample[1], Is.EqualTo(0.01920069443188622));
            Assert.That(t[2].CcfPerSample[1], Is.EqualTo(0.46199999999999997));
            Assert.That(t[2].StdDevPerSample[1], Is.EqualTo(0.05300000000000002));
            Assert.That(t[2].CcfPerSample[0], Is.EqualTo(0.0), "absent sample: centroid 0");
            Assert.That(p[0].CcfPerSample[3], Is.EqualTo(0.5));
            Assert.That(p[0].StdDevPerSample[3], Is.EqualTo(0.060000000000000026));
            Assert.That(p[1].StdDevPerSample[3], Is.EqualTo(0.0));
        });
    }

    #endregion

    #region F42 — per-cluster margins change the LICHeE network (lichee.jar)

    // p1: static e = 0.05 admits only B->A (A at s4 0.5 < 0.57 - 0.05) -> root->B->A, 1 tree, error
    // 0.06928203230275505 (jar, SD 0). With margins se_A = 1.96*0.06/sqrt(2) = 0.0832 both directions pass and
    // checkAndAddEdge keeps A->B (one-sided excess 0.07 < 0.12); A->B breaks the sum rule at s4 (0.57 > 0.5 + 0.05),
    // the complete network adds nothing (same level) -> lichee.jar: "Found 0 valid tree(s)".
    [Test]
    public void ClusterMargins_FlipEdgeOrientation_NoValidTreeLikeLichee()
    {
        OncologyAnalyzer.CcfClusterSummary[] clusters = P1();

        bool found = OncologyAnalyzer.TryReconstructPhylogenyFromClusterSummaries(clusters, out _, tolerance: 0.05);
        OncologyAnalyzer.ClonalPhylogeny staticTree =
            OncologyAnalyzer.ReconstructPhylogeny(clusters.Select(Centroid).ToArray(), tolerance: 0.05);

        Assert.Multiple(() =>
        {
            Assert.That(found, Is.False, "lichee.jar: 0 valid trees with per-cluster margins");
            Assert.That(staticTree.ParentOf(2), Is.EqualTo(staticTree.RootId));
            Assert.That(staticTree.ParentOf(1), Is.EqualTo(2));
            Assert.That(staticTree.ErrorScore, Is.EqualTo(0.06928203230275505), "lichee.jar static (SD 0) error score");
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(clusters, 0.05),
                NUnit.Framework.Throws.InvalidOperationException);
        });
    }

    // t02150 (e = 0): lichee.jar with margins -> complete network, 2 valid trees, top tree root->{1,2,3}, error 0;
    // static (SD 0) -> default network, 1 tree, 2->1, root->{2,3}, error 0.
    [Test]
    public void ClusterMargins_T02150_MatchesLicheeCompleteNetworkTree()
    {
        OncologyAnalyzer.CcfClusterSummary[] clusters = T02150();

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(clusters, 0.0);
        OncologyAnalyzer.ClonalPhylogeny s = OncologyAnalyzer.ReconstructPhylogeny(clusters.Select(Centroid).ToArray(), 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(p.UsedCompleteNetwork, Is.True);
            Assert.That(p.ValidTreeCount, Is.EqualTo(2));
            Assert.That(p.ErrorScore, Is.EqualTo(0.0));
            Assert.That(p.ParentOf(1), Is.EqualTo(p.RootId));
            Assert.That(p.ParentOf(2), Is.EqualTo(p.RootId));
            Assert.That(p.ParentOf(3), Is.EqualTo(p.RootId));
            Assert.That(s.UsedCompleteNetwork, Is.False);
            Assert.That(s.ValidTreeCount, Is.EqualTo(1));
            Assert.That(s.ParentOf(1), Is.EqualTo(2));
        });
    }

    // t01061 (e = 0.1): lichee.jar with margins -> 0 valid trees; static (SD 0) -> 1 tree root->2->3->4->1,
    // error 0.10016236818286589.
    [Test]
    public void ClusterMargins_T01061_NoValidTreeLikeLichee()
    {
        OncologyAnalyzer.CcfClusterSummary[] clusters = T01061();

        bool found = OncologyAnalyzer.TryReconstructPhylogenyFromClusterSummaries(clusters, out _, tolerance: 0.1);
        OncologyAnalyzer.ClonalPhylogeny s =
            OncologyAnalyzer.ReconstructPhylogeny(clusters.Select(Centroid).ToArray(), tolerance: 0.1);

        Assert.Multiple(() =>
        {
            Assert.That(found, Is.False);
            Assert.That(s.ParentOf(2), Is.EqualTo(s.RootId));
            Assert.That(s.ParentOf(3), Is.EqualTo(2));
            Assert.That(s.ParentOf(4), Is.EqualTo(3));
            Assert.That(s.ParentOf(1), Is.EqualTo(4));
            Assert.That(s.ErrorScore, Is.EqualTo(0.10016236818286589));
        });
    }

    // t00298 (e = 0): lichee.jar with margins -> complete network, 28 valid trees, top tree root->{1,3,4,5}, 4->2,
    // 1->6, error 0; static (SD 0) -> 5 trees, root->{1,4}, 4->2, 2->3, 1->6, 6->5.
    [Test]
    public void ClusterMargins_T00298_TopTreeMatchesLichee()
    {
        OncologyAnalyzer.CcfClusterSummary[] clusters = T00298();

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(clusters, 0.0);
        OncologyAnalyzer.ClonalPhylogeny s = OncologyAnalyzer.ReconstructPhylogeny(clusters.Select(Centroid).ToArray(), 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(p.UsedCompleteNetwork, Is.True);
            Assert.That(p.ValidTreeCount, Is.EqualTo(28));
            Assert.That(p.ErrorScore, Is.EqualTo(0.0));
            Assert.That(new[] { 1, 3, 4, 5 }.Select(p.ParentOf), Is.All.EqualTo(p.RootId));
            Assert.That(p.ParentOf(2), Is.EqualTo(4));
            Assert.That(p.ParentOf(6), Is.EqualTo(1));
            Assert.That(s.ValidTreeCount, Is.EqualTo(5));
            Assert.That(s.ParentOf(3), Is.EqualTo(2));
            Assert.That(s.ParentOf(5), Is.EqualTo(6));
        });
    }

    #endregion

    // F44 — the same fixtures at the default e = 0.1 (lichee.jar -e 0.1): t02150 margins -> default network, 1 tree,
    // root->2, 2->{1,3}, error 0.016000000000000014 (static, the Java centroids as 2 identical members: same tree and
    // error); t00298
    // margins -> complete network, 44 trees, same top tree as at e = 0; f1 -> nothing removed, root->1->2, error
    // 0.050000000000000044.
    [Test]
    public void DefaultTolerance_ReDerivedFixtures_MatchLichee()
    {
        OncologyAnalyzer.ClonalPhylogeny t02150 = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(T02150());
        OncologyAnalyzer.ClonalPhylogeny t02150s = OncologyAnalyzer.ReconstructPhylogeny(T02150().Select(Centroid).ToArray());
        OncologyAnalyzer.ClonalPhylogeny t00298 = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(T00298());
        OncologyAnalyzer.ClonalPhylogeny f1 = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(
            new[] { M(1, R(0.5, 0.5)), M(2, R(0.55, 0)) });

        Assert.Multiple(() =>
        {
            Assert.That(t02150.UsedCompleteNetwork, Is.False);
            Assert.That(t02150.ValidTreeCount, Is.EqualTo(1));
            Assert.That(t02150.ErrorScore, Is.EqualTo(0.016000000000000014));
            Assert.That(new[] { 1, 2, 3 }.Select(t02150.ParentOf), Is.EqualTo(new int?[] { 2, t02150.RootId, 2 }));
            Assert.That(t02150s.ErrorScore, Is.EqualTo(0.016000000000000014));
            Assert.That(new[] { 1, 2, 3 }.Select(t02150s.ParentOf), Is.EqualTo(new int?[] { 2, t02150s.RootId, 2 }));
            Assert.That(t00298.UsedCompleteNetwork, Is.True);
            Assert.That(t00298.ValidTreeCount, Is.EqualTo(44));
            Assert.That(t00298.ErrorScore, Is.EqualTo(0.0));
            Assert.That(new[] { 1, 3, 4, 5 }.Select(t00298.ParentOf), Is.All.EqualTo(t00298.RootId));
            Assert.That(t00298.ParentOf(2), Is.EqualTo(4));
            Assert.That(t00298.ParentOf(6), Is.EqualTo(1));
            Assert.That(f1.RemovedClusterIds, Is.Empty);
            Assert.That(f1.ParentOf(1), Is.EqualTo(f1.RootId));
            Assert.That(f1.ParentOf(2), Is.EqualTo(1));
            Assert.That(f1.ErrorScore, Is.EqualTo(0.050000000000000044));
        });
    }

    #region F42 — static behaviour preserved

    // With every SD = 0 the margin is max(e, 0 + 0) = e (root: e + 0 = e), i.e. exactly the static network: the
    // F18 complete-network fixture gives the same tree and error score bit for bit.
    [Test]
    public void ZeroSd_IdenticalToStaticReconstruction()
    {
        double[][] ccf =
        {
            new[] { 0.482, 0.0, 0.443137 },
            new[] { 0.519, 0.796779, 0.56 },
            new[] { 0.19, 0.418589, 0.24 },
        };
        var summaries = ccf.Select((v, i) => new OncologyAnalyzer.CcfClusterSummary(i + 1, v, new double[3], 3)).ToArray();

        OncologyAnalyzer.ClonalPhylogeny a = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(summaries, 0.02);
        OncologyAnalyzer.ClonalPhylogeny b = OncologyAnalyzer.ReconstructPhylogeny(summaries.Select(Centroid).ToArray(), 0.02);

        Assert.Multiple(() =>
        {
            Assert.That(a.Edges, Is.EqualTo(b.Edges));
            Assert.That(a.ErrorScore, Is.EqualTo(b.ErrorScore));
            Assert.That(a.ErrorScore, Is.EqualTo(0.003292532308117998).Within(1e-15));
            Assert.That(a.ValidTreeCount, Is.EqualTo(b.ValidTreeCount));
            Assert.That(a.UsedCompleteNetwork, Is.EqualTo(b.UsedCompleteNetwork));
        });
    }

    // PHYNode.getStdDev returns 0 for a sample outside the cluster's presence profile, so an SD supplied where the
    // centroid is 0 has no effect.
    [Test]
    public void SdAtAbsentSample_Ignored()
    {
        OncologyAnalyzer.CcfClusterSummary[] baseClusters = T02150();
        OncologyAnalyzer.CcfClusterSummary[] noisy = baseClusters
            .Select(c => c with { StdDevPerSample = c.CcfPerSample.Select((v, i) => v > 0 ? c.StdDevPerSample[i] : 0.9).ToArray() })
            .ToArray();

        OncologyAnalyzer.ClonalPhylogeny a = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(baseClusters);
        OncologyAnalyzer.ClonalPhylogeny b = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(noisy);

        Assert.That(b.Edges, Is.EqualTo(a.Edges));
    }

    [Test]
    public void EmptyInput_ReturnsRootOnlyTree()
    {
        OncologyAnalyzer.ClonalPhylogeny p =
            OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(Array.Empty<OncologyAnalyzer.CcfClusterSummary>());

        Assert.That(p.Edges, Is.Empty);
    }

    #endregion

    #region F42 — validation

    [Test]
    public void Validation_BadSummaries_Throw()
    {
        double[] ccf = { 0.5, 0.4 };
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(null!), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(
                new[] { new OncologyAnalyzer.CcfClusterSummary(1, ccf, null!, 2) }), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(
                new[] { new OncologyAnalyzer.CcfClusterSummary(1, ccf, new[] { 0.1 }, 2) }), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(
                new[] { new OncologyAnalyzer.CcfClusterSummary(1, ccf, new[] { 0.1, -0.1 }, 2) }), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(
                new[] { new OncologyAnalyzer.CcfClusterSummary(1, ccf, new[] { 0.1, double.NaN }, 2) }), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(
                new[] { new OncologyAnalyzer.CcfClusterSummary(1, ccf, new[] { 0.1, 0.1 }, 0) }), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(
                new[] { new OncologyAnalyzer.CcfClusterSummary(1, ccf, new[] { 0.1, 0.1 }, 2) }, tolerance: -1),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => OncologyAnalyzer.CcfClusterSummary.FromMembers(1, Array.Empty<double[]>()), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OncologyAnalyzer.CcfClusterSummary.FromMembers(1, new[] { R(0.1, 0.2), R(0.3) }), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OncologyAnalyzer.CcfClusterSummary.FromMembers(1, new[] { R() }), NUnit.Framework.Throws.ArgumentException);
        });
    }

    #endregion

    #region F43 — fixNetwork (LineageEngine.buildLineage step 5) vs lichee.jar

    private static OncologyAnalyzer.CcfClusterSummary MR(int id, int robustMembers, params double[][] members) =>
        OncologyAnalyzer.CcfClusterSummary.FromMembers(id, members) with { RobustMemberCount = robustMembers };

    // f1 (= F18 S3 repro, ε 0): A=[0.5,0.5], B=[0.55,0], one member each (non-robust, size 1). No tree; fixNetwork
    // scans nodes in id order with a strict '<' on size, so the tie keeps A (node 1) -> A dropped; root->B.
    // lichee.jar: "removed [1]", 1 tree, error 0, edge -1 -> 2. Without fixNetwork: no tree (F18).
    [Test]
    public void FixNetwork_SizeTie_DropsFirstNodeLikeLichee()
    {
        var clusters = new[] { M(1, R(0.5, 0.5)), M(2, R(0.55, 0)) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(clusters, 0.0);
        bool keepAll = OncologyAnalyzer.TryReconstructPhylogenyFromClusterSummaries(
            clusters, out _, 0.0, removeNonRobustClusters: false);

        Assert.Multiple(() =>
        {
            Assert.That(p.RemovedClusterIds, Is.EqualTo(new[] { 1 }));
            Assert.That(p.Clusters.Select(c => c.Id), Is.EqualTo(new[] { 2 }));
            Assert.That(p.Edges, Is.EqualTo(new[] { new OncologyAnalyzer.ClonalEdge(p.RootId, 2) }));
            Assert.That(p.ValidTreeCount, Is.EqualTo(1));
            Assert.That(p.ErrorScore, Is.EqualTo(0.0));
            Assert.That(p.UsedCompleteNetwork, Is.False);
            Assert.That(keepAll, Is.False, "LICHeE without fixNetwork: no valid tree");
        });
    }

    // f2: A = 3 members, 1 robust (non-robust, size 3); B = 1 member (non-robust, size 1) -> the smallest (B) is
    // dropped. lichee.jar: "removed [2]", edge -1 -> 1.
    [Test]
    public void FixNetwork_DropsSmallestNonRobustCluster()
    {
        var clusters = new[]
        {
            MR(1, 1, R(0.5, 0.5), R(0.5, 0.5), R(0.5, 0.5)),
            M(2, R(0.55, 0)),
        };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(clusters, 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(p.RemovedClusterIds, Is.EqualTo(new[] { 2 }));
            Assert.That(p.ParentOf(1), Is.EqualTo(p.RootId));
            Assert.That(p.ValidTreeCount, Is.EqualTo(1));
        });
    }

    // f4 / f4r: LICHeE's --clustersFile input never marks a cluster robust (RobustMemberCount = 0): A (2 members) is
    // dropped before B (3) -> root->B (jar "removed [1]"). The same clusters with robust members (default) are both
    // robust: nothing to drop, complete network also infeasible -> jar "trees 0".
    [Test]
    public void FixNetwork_ClustersFileInputWithoutRobustMembers_MatchesLichee()
    {
        var nonRobust = new[]
        {
            MR(1, 0, R(0.5, 0.5), R(0.5, 0.5)),
            MR(2, 0, R(0.55, 0), R(0.55, 0), R(0.55, 0)),
        };
        var robust = nonRobust.Select(c => c with { RobustMemberCount = null }).ToArray();

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(nonRobust, 0.0);
        bool robustFound = OncologyAnalyzer.TryReconstructPhylogenyFromClusterSummaries(robust, out _, 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(p.RemovedClusterIds, Is.EqualTo(new[] { 1 }));
            Assert.That(p.ParentOf(2), Is.EqualTo(p.RootId));
            Assert.That(robustFound, Is.False);
        });
    }

    // f6: A robust (2 members) [0.5,0.5]; B=[0.55,0], C=[0,0.55] single members. Removing B (first size-1 node) still
    // leaves root children 0.5 + 0.55 > 1 in sample 2 -> C removed next -> root->A (jar "removed [2, 3]").
    [Test]
    public void FixNetwork_RemovesRepeatedlyUntilTreeFound()
    {
        var clusters = new[] { M(1, R(0.5, 0.5), R(0.5, 0.5)), M(2, R(0.55, 0)), M(3, R(0, 0.55)) };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(clusters, 0.0);

        Assert.Multiple(() =>
        {
            Assert.That(p.RemovedClusterIds, Is.EqualTo(new[] { 2, 3 }));
            Assert.That(p.Edges, Is.EqualTo(new[] { new OncologyAnalyzer.ClonalEdge(p.RootId, 1) }));
            Assert.That(OncologyAnalyzer.IdentifyBranchMutations(p), Is.EqualTo(new[] { 1 }),
                "removed clusters are neither trunk nor branch");
        });
    }

    // f3: the F18 complete-network fixture (each cluster 2 robust members) plus a non-robust single-member D =
    // [0.01,0.01,0.01]. fixNetwork drops D, the default network is still infeasible, so ALL_EDGES runs on the
    // reduced set (LICHeE rebuilds from the mutated groups): jar "removed [4]", complete network, root->{1,2}, 2->3,
    // error 0.003292532308117998.
    [Test]
    public void FixNetwork_ThenCompleteNetworkOnRemainingClusters()
    {
        var clusters = new[]
        {
            M(1, R(0.482, 0, 0.443137), R(0.482, 0, 0.443137)),
            M(2, R(0.519, 0.796779, 0.56), R(0.519, 0.796779, 0.56)),
            M(3, R(0.19, 0.418589, 0.24), R(0.19, 0.418589, 0.24)),
            M(4, R(0.01, 0.01, 0.01)),
        };

        OncologyAnalyzer.ClonalPhylogeny p = OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(clusters, 0.02);

        Assert.Multiple(() =>
        {
            Assert.That(p.RemovedClusterIds, Is.EqualTo(new[] { 4 }));
            Assert.That(p.UsedCompleteNetwork, Is.True);
            Assert.That(p.ParentOf(1), Is.EqualTo(p.RootId));
            Assert.That(p.ParentOf(2), Is.EqualTo(p.RootId));
            Assert.That(p.ParentOf(3), Is.EqualTo(2));
            Assert.That(p.ValidTreeCount, Is.EqualTo(1));
            Assert.That(p.ErrorScore, Is.EqualTo(0.003292532308117998));
        });
    }

    // f7: removals [2, 4] leave robust A and C conflicting; complete network infeasible too -> jar "trees 0".
    [Test]
    public void FixNetwork_NoTreeAfterRemovals_TryReturnsFalse()
    {
        var clusters = new[]
        {
            M(1, R(0.5, 0.5), R(0.5, 0.5)),
            M(2, R(0.55, 0)),
            M(3, R(0, 0.55), R(0, 0.55)),
            M(4, R(0, 0.3)),
        };

        bool found = OncologyAnalyzer.TryReconstructPhylogenyFromClusterSummaries(clusters, out _, 0.0);

        Assert.That(found, Is.False);
    }

    [Test]
    public void FixNetwork_RobustMemberCountOutOfRange_Throws()
    {
        var bad = new[] { MR(1, 3, R(0.5), R(0.4)) };
        var negative = new[] { MR(1, -1, R(0.5), R(0.4)) };

        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(bad), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OncologyAnalyzer.ReconstructPhylogenyFromClusterSummaries(negative), NUnit.Framework.Throws.ArgumentException);
        });
    }

    #endregion
}
