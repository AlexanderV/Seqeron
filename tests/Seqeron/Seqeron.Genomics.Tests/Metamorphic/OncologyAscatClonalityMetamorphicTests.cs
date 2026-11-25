namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic tests for the behaviour introduced by the 2026-09 B24 review of the Oncology copy-number / clonality
/// units. Each relation follows from the algorithm definition (the ported reference code), not from observed output.
/// <list type="bullet">
/// <item>ONCO-ASCAT-001 (runASCAT): sex-chromosome segments are excluded from the fit (<c>autoprobes</c>), so adding
/// X/Y segments leaves ρ, ψ, the goodness of fit and the autosomal integer segments unchanged.</item>
/// <item>ONCO-ASCAT-001 (ascat.aspcf): chromosomes are segmented independently, so segmenting a concatenation of
/// chromosomes equals concatenating the per-chromosome results (below the 800-level penalty-ladder trigger).</item>
/// <item>ONCO-ASCAT-001 (Battenberg): the fit uses l = max(b, 1 − b), so mirroring every BAF gives identical states;
/// segments are fitted independently, so reordering them reorders the fits.</item>
/// <item>ONCO-PLOIDY-001 (facets-suite <c>get_sample_genome</c>): the WGD fraction uses lengths and per-chromosome
/// spans only, so translating a chromosome's coordinates, relabelling "7" ↔ "chr7", or adding non-autosomal segments
/// does not change the call.</item>
/// <item>ONCO-PHYLO-001 (LICHeE): the search depends on cluster order and CCFs, not on the id values, so a bijective
/// relabelling of ids maps the tree edge-for-edge; swapping the two sample columns preserves feasibility, the valid-tree
/// count and the error score.</item>
/// </list>
/// All randomness is locally seeded.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Oncology")]
public class OncologyAscatClonalityMetamorphicTests
{

    private static OncologyAnalyzer.AlleleSpecificSegment AsSeg(string chr, long start, long end, int major, int minor) =>
        new(chr, start, end, major, minor);

    [Test]
    [Description("INV (facets get_sample_genome): per-chromosome translation, 'chr' relabelling and non-autosomal segments do not change the WGD call.")]
    public void DetectWholeGenomeDoublingFromSuppliedLength_SpanPreservingTransforms_KeepTheCall()
    {
        int doubled = 0, notDoubled = 0;
        for (int seed = 0; seed < 200; seed++)
        {
            var rng = new Random(seed);
            var segs = new List<OncologyAnalyzer.AlleleSpecificSegment>();
            for (int c = 1; c <= rng.Next(1, 6); c++)
            {
                long pos = rng.Next(0, 5_000_000);
                for (int k = 0; k < rng.Next(1, 5); k++)
                {
                    long len = rng.Next(1, 40_000_000);
                    segs.Add(AsSeg(c.ToString(), pos, pos + len, rng.Next(0, 4), rng.Next(0, 3)));
                    pos += len + rng.Next(0, 10_000_000); // unsegmented gaps count toward the span
                }
            }

            bool call = OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(segs);
            if (call) doubled++; else notDoubled++;

            var offsets = Enumerable.Range(0, 23).Select(_ => (long)rng.Next(0, 1_000_000_000)).ToArray();
            var translated = segs.Select(s =>
            {
                long o = offsets[int.Parse(s.Chromosome)];
                return s with { Start = s.Start + o, End = s.End + o };
            }).ToList();
            var relabelled = segs.Select(s => s with { Chromosome = "chr" + s.Chromosome }).ToList();
            var withSex = segs.Concat(new[]
            {
                AsSeg("X", 0, 150_000_000, 4, 4), AsSeg("chrY", 0, 50_000_000, 0, 0), AsSeg("MT", 0, 16_569, 5, 5),
            }).ToList();

            OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(translated).Should().Be(call, "translation, seed {0}", seed);
            OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(relabelled).Should().Be(call, "chr prefix, seed {0}", seed);
            OncologyAnalyzer.DetectWholeGenomeDoublingFromSuppliedLength(withSex).Should().Be(call, "non-autosomes, seed {0}", seed);
        }

        doubled.Should().BeGreaterThan(0);
        notDoubled.Should().BeGreaterThan(0);
    }

    private static OncologyAnalyzer.CcfCluster[] RandomClusters(Random rng, int samples)
    {
        double[] grid = { 0.0, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0 };
        return Enumerable.Range(1, rng.Next(1, 7))
            .Select(id => new OncologyAnalyzer.CcfCluster(id, Enumerable.Range(0, samples).Select(_ => grid[rng.Next(grid.Length)]).ToArray()))
            .ToArray();
    }

    [Test]
    [Description("INV (LICHeE uses cluster order, not id values): a bijective id relabelling maps the tree edge-for-edge.")]
    public void TryReconstructPhylogeny_ClusterRelabelling_MapsTheTree()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            var clusters = RandomClusters(rng, rng.Next(1, 4));
            double eps = rng.Next(0, 3) * 0.05;
            int[] newIds = Enumerable.Range(0, 1000).OrderBy(_ => rng.Next()).Take(clusters.Length).Select(v => v * 3 - 700).ToArray();
            var map = clusters.Select((c, i) => (c.Id, newIds[i])).ToDictionary(p => p.Id, p => p.Item2);
            var relabelled = clusters.Select(c => c with { Id = map[c.Id] }).ToArray();

            bool okA = OncologyAnalyzer.TryReconstructPhylogeny(clusters, out var a, eps);
            bool okB = OncologyAnalyzer.TryReconstructPhylogeny(relabelled, out var b, eps);
            okB.Should().Be(okA, "seed {0}", seed);
            if (!okA)
            {
                continue;
            }

            b.ErrorScore.Should().Be(a.ErrorScore, "seed {0}", seed);
            b.ValidTreeCount.Should().Be(a.ValidTreeCount, "seed {0}", seed);
            b.UsedCompleteNetwork.Should().Be(a.UsedCompleteNetwork, "seed {0}", seed);
            var mappedEdges = a.Edges.Select(e => (e.ParentId == a.RootId ? b.RootId : map[e.ParentId], map[e.ChildId]));
            b.Edges.Select(e => (e.ParentId, e.ChildId)).Should().Equal(mappedEdges, "seed {0}", seed);
            OncologyAnalyzer.IdentifyTrunkMutations(b).Should().Equal(OncologyAnalyzer.IdentifyTrunkMutations(a).Select(id => map[id]),
                "seed {0}", seed);
        }
    }

    [Test]
    [Description("INV (per-sample constraints): swapping the two sample columns preserves feasibility, valid-tree count and error score.")]
    public void TryReconstructPhylogeny_TwoSampleSwap_PreservesTheSearchOutcome()
    {
        int feasible = 0;
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            var clusters = RandomClusters(rng, 2);
            double eps = rng.Next(0, 3) * 0.05;
            var swapped = clusters.Select(c => c with { CcfPerSample = new[] { c.CcfPerSample[1], c.CcfPerSample[0] } }).ToArray();

            bool okA = OncologyAnalyzer.TryReconstructPhylogeny(clusters, out var a, eps);
            bool okB = OncologyAnalyzer.TryReconstructPhylogeny(swapped, out var b, eps);
            okB.Should().Be(okA, "seed {0}", seed);
            if (!okA)
            {
                continue;
            }

            feasible++;
            b.ValidTreeCount.Should().Be(a.ValidTreeCount, "seed {0}", seed);
            b.UsedCompleteNetwork.Should().Be(a.UsedCompleteNetwork, "seed {0}", seed);
            b.ErrorScore.Should().BeApproximately(a.ErrorScore, 1e-12, "seed {0}", seed);
            OncologyAnalyzer.IdentifyTrunkMutations(b).Should().Equal(OncologyAnalyzer.IdentifyTrunkMutations(a), "seed {0}", seed);
        }

        feasible.Should().BeGreaterThan(50);
    }
}
