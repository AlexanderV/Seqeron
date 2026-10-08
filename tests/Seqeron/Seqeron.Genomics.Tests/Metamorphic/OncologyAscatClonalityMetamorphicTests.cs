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
/// <item>ONCO-CCF-001 (Ckmeans.1d.dp): the optimum is computed on the sorted values, so permuting the input gives the
/// same centroids and the same (value, cluster) pairs.</item>
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
    private static OncologyAnalyzer.AlleleSpecificSegmentSummary Seg(string chr, long start, double logR, double baf, int probes) =>
        new(chr, start, start + 1_000_000, logR, baf, probes);

    private static List<OncologyAnalyzer.AlleleSpecificSegmentSummary> ModelGenome(Random rng, double rho, double psi, int n)
    {
        var segs = new List<OncologyAnalyzer.AlleleSpecificSegmentSummary>();
        for (int i = 0; i < n; i++)
        {
            int major = rng.Next(0, 5), minor = rng.Next(0, major + 1);
            if (major + minor == 0)
            {
                major = 1;
            }

            int total = major + minor;
            double logR = Math.Log2((2 * (1 - rho) + rho * total) / (2 * (1 - rho) + rho * psi)) + (rng.NextDouble() - 0.5) * 0.02;
            double baf = (1 - rho + rho * minor) / (2 * (1 - rho) + rho * total);
            segs.Add(Seg((1 + i % 22).ToString(), i * 2_000_000L, logR, baf, rng.Next(5, 300)));
        }

        return segs;
    }

    [Test]
    [Description("INV (ASCAT autoprobes): appending X/Y segments leaves ρ, ψ, GoF and the autosomal integer segments unchanged.")]
    public void FitPurityPloidy_AddingSexChromosomeSegments_LeavesTheFitUnchanged()
    {
        int found = 0;
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(seed);
            double rho = 0.3 + 0.05 * rng.Next(0, 14), psi = 1.6 + 0.05 * rng.Next(0, 60);
            var autosomal = ModelGenome(rng, rho, psi, rng.Next(3, 8));
            var withSex = autosomal.ToList();
            withSex.Add(Seg("X", 0, -0.8 + rng.NextDouble(), rng.NextDouble(), rng.Next(5, 500)));
            withSex.Add(Seg("chrY", 0, -1.0 + rng.NextDouble(), 0.1, rng.Next(5, 50)));

            bool okA = OncologyAnalyzer.TryFitPurityPloidy(autosomal, out var a);
            bool okB = OncologyAnalyzer.TryFitPurityPloidy(withSex, out var b);
            okB.Should().Be(okA, "sex-chromosome probes do not enter the ASCAT fit (seed {0})", seed);
            if (!okA)
            {
                continue;
            }

            found++;
            b.Purity.Should().Be(a.Purity, "seed {0}", seed);
            b.Psi.Should().Be(a.Psi, "seed {0}", seed);
            b.GoodnessOfFit.Should().Be(a.GoodnessOfFit, "seed {0}", seed);
            b.IsNonAberrant.Should().Be(a.IsNonAberrant, "seed {0}", seed);
            b.Segments.Take(autosomal.Count).Should().Equal(a.Segments, "autosomal integer segments unchanged (seed {0})", seed);
        }

        found.Should().BeGreaterThan(10, "most model genomes have an ASCAT optimum");
    }

    [Test]
    [Description("INV (ascat.aspcf per-chromosome loop): ASPCF(chr1 ++ chr2 ++ chr3) = ASPCF(chr1) ++ ASPCF(chr2) ++ ASPCF(chr3).")]
    public void SegmentAlleleSpecificAspcf_ChromosomesAreSegmentedIndependently()
    {
        for (int seed = 0; seed < 25; seed++)
        {
            var rng = new Random(seed);
            var blocks = new List<OncologyAnalyzer.AlleleSpecificLocus[]>();
            for (int c = 1; c <= 3; c++)
            {
                int n = rng.Next(1, 150);
                double level = rng.Next(-2, 3) * 0.4, shift = rng.Next(0, 4) * 0.1;
                var block = new OncologyAnalyzer.AlleleSpecificLocus[n];
                for (int i = 0; i < n; i++)
                {
                    if (i == n / 2)
                    {
                        level += 0.6;
                        shift = 0.3 - shift;
                    }

                    double baf = Math.Clamp(0.5 + (rng.Next(2) == 0 ? shift : -shift) + (rng.NextDouble() - 0.5) * 0.08, 0.0, 1.0);
                    block[i] = new OncologyAnalyzer.AlleleSpecificLocus(c.ToString(), i * 500L, level + (rng.NextDouble() - 0.5) * 0.3, baf);
                }

                blocks.Add(block);
            }

            var whole = OncologyAnalyzer.SegmentAlleleSpecificAspcf(blocks.SelectMany(b => b));
            var parts = blocks.SelectMany(b => OncologyAnalyzer.SegmentAlleleSpecificAspcf(b)).ToList();
            whole.Should().Equal(parts, "chromosomes are segmented independently (seed {0})", seed);
        }
    }

    [Test]
    [Description("INV (Battenberg l = max(b, 1 − b)): mirroring every BAF gives identical states; reordering segments reorders the fits.")]
    public void FitSubclonalCopyNumber_BafMirrorAndSegmentOrder_AreInvariant()
    {
        for (int seed = 0; seed < 50; seed++)
        {
            var rng = new Random(seed);
            double rho = 0.1 + 0.9 * rng.NextDouble(), psi = 1.5 + 3.5 * rng.NextDouble();
            var segs = Enumerable.Range(0, rng.Next(1, 12))
                .Select(i => Seg("1", i * 2_000_000L, (rng.NextDouble() - 0.5) * 3.0, rng.NextDouble(), 20))
                .ToArray();
            var mirrored = segs.Select(s => s with { MeanBAF = 1.0 - s.MeanBAF }).ToArray();

            var fits = OncologyAnalyzer.FitSubclonalCopyNumber(segs, rho, psi);
            var mirroredFits = OncologyAnalyzer.FitSubclonalCopyNumber(mirrored, rho, psi);
            for (int i = 0; i < segs.Length; i++)
            {
                mirroredFits[i].PrimaryState.Should().Be(fits[i].PrimaryState, "seed {0} segment {1}", seed, i);
                mirroredFits[i].SecondaryState.Should().Be(fits[i].SecondaryState, "seed {0} segment {1}", seed, i);
                mirroredFits[i].IsSubclonal.Should().Be(fits[i].IsSubclonal, "seed {0} segment {1}", seed, i);
            }

            int[] order = Enumerable.Range(0, segs.Length).OrderBy(_ => rng.Next()).ToArray();
            var reordered = OncologyAnalyzer.FitSubclonalCopyNumber(order.Select(i => segs[i]).ToArray(), rho, psi);
            reordered.Should().Equal(order.Select(i => fits[i]), "segments are fitted independently (seed {0})", seed);
        }
    }

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

    [Test]
    [Description("INV (Ckmeans.1d.dp works on sorted data): permuting the CCFs gives identical centroids and identical (value, cluster) pairs.")]
    public void ClusterCcfValues_InputPermutation_IsInvariant()
    {
        for (int seed = 0; seed < 300; seed++)
        {
            var rng = new Random(seed);
            int n = rng.Next(1, 40);
            double[] values = Enumerable.Range(0, n).Select(_ => rng.Next(0, 60) / 50.0 * 0.8).ToArray();
            int k = rng.Next(1, n + 1);
            double[] permuted = values.OrderBy(_ => rng.Next()).ToArray();

            var a = OncologyAnalyzer.ClusterCcfValues(values, k);
            var b = OncologyAnalyzer.ClusterCcfValues(permuted, k);

            b.Centroids.Should().Equal(a.Centroids, "seed {0}", seed);
            b.ClonalClusterIndex.Should().Be(a.ClonalClusterIndex, "seed {0}", seed);
            var pairsA = values.Select((v, i) => (v, a.Assignments[i])).OrderBy(p => p.v).ThenBy(p => p.Item2);
            var pairsB = permuted.Select((v, i) => (v, b.Assignments[i])).OrderBy(p => p.v).ThenBy(p => p.Item2);
            pairsB.Should().Equal(pairsA, "each value keeps its cluster (seed {0})", seed);
        }
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
