using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the behaviour introduced by the 2026-09 B24 review (ONCO-ASCAT-001, ONCO-CCF-001,
/// ONCO-PHYLO-001, ONCO-HETERO-001): the runASCAT / ASCAT-manual fit, the ascat.aspcf segmentation, the Battenberg
/// sub-clonal fit, the LICHeE trunk rule and the maftools MATH score.
///
/// Oracles are recomputed independently from the published equations (never routed through production helpers):
/// <list type="bullet">
/// <item>ASCAT model (Van Loo et al. 2010, PNAS 107:16910; ascat.runAscat.R): a clonal segment with integer allele copy
/// numbers (nA, nB) at purity ρ and ploidy ψ has logR r = log2((2(1−ρ)+ρ(nA+nB)) / (2(1−ρ)+ρψ)) (γ = 1) and BAF
/// b = (1−ρ+ρ·nB) / (2(1−ρ)+ρ(nA+nB)).</item>
/// <item>Battenberg (Nik-Zainal et al. 2012; Wedge-lab/battenberg R/fitcopynumber.R): the sub-clonal fraction τ solves
/// l = [1−ρ+ρ(τM₁+(1−τ)M₂)] / [2−2ρ+ρ(τT₁+(1−τ)T₂)] for the two states of the nearest edge (which differ by one copy
/// of one allele); a clonal call has |level(state) − l| &lt; maxdist = 0.01.</item>
/// <item>ascat.aspcf: each segment's logR level is the mean of the raw logR of its loci; its BAF is 0.5 + μ ∈ [0.5, 1].</item>
/// <item>Werner et al. 2017 (Sci Rep 7:44991): trunk clusters are present in all cells of every sample (CCF ≥ 1 − ε).</item>
/// <item>MATH (Mroz &amp; Rocco 2013; maftools mathScore.R) is a MAD/median ratio: invariant to input order and to
/// scaling every VAF by c &gt; 0.</item>
/// </list>
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Oncology")]
public class OncologyAscatClonalityProperties
{
    // -------------------------------------------------------------------------
    // ASCAT exact-model generators (Van Loo 2010 equations, γ = 1)
    // -------------------------------------------------------------------------

    private static double ModelLogR(double rho, double psi, int total) =>
        Math.Log2((2.0 * (1.0 - rho) + rho * total) / (2.0 * (1.0 - rho) + rho * psi));

    private static double ModelMinorBaf(double rho, int minor, int total) =>
        (1.0 - rho + rho * minor) / (2.0 * (1.0 - rho) + rho * total);

    private static Gen<(double Rho, double Psi, (int Major, int Minor, int Probes)[] Karyotype)> ExactAscatGenomeGen() =>
        from rhoPercent in Gen.Choose(20, 100)
        from psiTwentieths in Gen.Choose(30, 100) // ψ ∈ [1.5, 5.0] step 0.05
        from n in Gen.Choose(1, 6)
        from karyotype in (from major in Gen.Choose(0, 5)
                           from minor in Gen.Choose(0, 5)
                           from probes in Gen.Choose(1, 200)
                           select (Major: Math.Max(major, minor), Minor: Math.Min(major, minor), Probes: probes))
            .Where(k => k.Major + k.Minor >= 1)
            .ArrayOf(n)
        select (rhoPercent / 100.0, psiTwentieths / 20.0, karyotype);

    private static OncologyAnalyzer.AlleleSpecificSegmentSummary[] ToSummaries(
        double rho, double psi, (int Major, int Minor, int Probes)[] karyotype, bool upperBafOrientation)
    {
        var segments = new OncologyAnalyzer.AlleleSpecificSegmentSummary[karyotype.Length];
        long start = 0;
        for (int i = 0; i < karyotype.Length; i++)
        {
            (int major, int minor, int probes) = karyotype[i];
            double b = ModelMinorBaf(rho, minor, major + minor);
            segments[i] = new OncologyAnalyzer.AlleleSpecificSegmentSummary(
                (1 + i % 22).ToString(), start, start + 1_000_000, ModelLogR(rho, psi, major + minor),
                upperBafOrientation ? 1.0 - b : b, probes);
            start += 2_000_000;
        }

        return segments;
    }

    /// <summary>
    /// ASCAT manual fit (rho_manual/psi_manual) on an exact clonal genome generated from the ASCAT equations at (ρ, ψ):
    /// the distance to the integer lattice is 0, so GoF = 100 %, every segment is recovered as its generating
    /// (major, minor) integer pair, and the output ploidy is the probe-weighted mean total copy number
    /// (ascat.runAscat.R). Checked in both BAF orientations (b and 1 − b are the same heterozygous signal).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EvaluatePurityPloidy_ExactModelGenome_RecoversKaryotypeWithFullGoodnessOfFit()
    {
        return Prop.ForAll(ExactAscatGenomeGen().ToArbitrary(), g =>
        {
            foreach (bool upper in new[] { false, true })
            {
                var segments = ToSummaries(g.Rho, g.Psi, g.Karyotype, upper);
                var fit = OncologyAnalyzer.EvaluatePurityPloidy(segments, g.Rho, g.Psi);

                if (fit.GoodnessOfFit < 100.0 - 1e-6 || fit.GoodnessOfFit > 100.0 + 1e-9)
                {
                    return false.Label($"GoF {fit.GoodnessOfFit} ≠ 100 (ρ={g.Rho}, ψ={g.Psi}, upper={upper})");
                }

                for (int i = 0; i < segments.Length; i++)
                {
                    var s = fit.Segments[i];
                    if (s.MajorCopyNumber != g.Karyotype[i].Major || s.MinorCopyNumber != g.Karyotype[i].Minor)
                    {
                        return false.Label(
                            $"segment {i}: {s.MajorCopyNumber}:{s.MinorCopyNumber} ≠ {g.Karyotype[i].Major}:{g.Karyotype[i].Minor}");
                    }
                }

                double probes = g.Karyotype.Sum(k => (double)k.Probes);
                double ploidy = g.Karyotype.Sum(k => (double)(k.Major + k.Minor) * k.Probes) / probes;
                if (Math.Abs(fit.Ploidy - ploidy) > 1e-12 * Math.Max(1.0, ploidy) || fit.Purity != g.Rho || fit.Psi != g.Psi)
                {
                    return false.Label($"ploidy {fit.Ploidy} ≠ {ploidy} or (ρ, ψ) not echoed");
                }
            }

            return true.Label("ok");
        });
    }

    /// <summary>
    /// Consistency of the grid fit with the manual fit (both are runASCAT paths): whenever <c>TryFitPurityPloidy</c>
    /// finds an optimum, <c>EvaluatePurityPloidy</c> at the reported (ρ, ψ) emits the identical integer segments, output
    /// ploidy and non-aberrant flag, and — when the optimum lies at ρ ≤ 1 (a grid ρ &gt; 1 is reported as 1) — the same
    /// goodness of fit.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 40)]
    public Property TryFitPurityPloidy_AgreesWithEvaluateAtTheReportedOptimum()
    {
        return Prop.ForAll(ExactAscatGenomeGen().ToArbitrary(), g =>
        {
            var segments = ToSummaries(g.Rho, g.Psi, g.Karyotype, upperBafOrientation: false);
            if (!OncologyAnalyzer.TryFitPurityPloidy(segments, out var fit))
            {
                return true.Label("ASCAT: rho = NA");
            }

            var manual = OncologyAnalyzer.EvaluatePurityPloidy(segments, fit.Purity, fit.Psi);
            bool same = manual.Segments.SequenceEqual(fit.Segments) && manual.Ploidy == fit.Ploidy
                        && manual.IsNonAberrant == fit.IsNonAberrant;
            bool gofOk = fit.Purity >= 1.0 || Math.Abs(manual.GoodnessOfFit - fit.GoodnessOfFit) <= 1e-9;
            return (same && gofOk).Label(
                $"ρ={fit.Purity}, ψ={fit.Psi}: GoF grid {fit.GoodnessOfFit} vs manual {manual.GoodnessOfFit}, same={same}");
        });
    }

    // -------------------------------------------------------------------------
    // Battenberg sub-clonal fit
    // -------------------------------------------------------------------------

    private static double CornerLevel(double rho, int major, int minor) =>
        (1.0 - rho + rho * major) / (2.0 - 2.0 * rho + rho * (major + minor));

    private static Gen<(double Rho, double Psi, OncologyAnalyzer.AlleleSpecificSegmentSummary[] Segments)> BattenbergGen() =>
        from rhoPercent in Gen.Choose(5, 100)
        from psiTwentieths in Gen.Choose(30, 100)
        from n in Gen.Choose(1, 8)
        from segs in (from logRMilli in Gen.Choose(-1500, 1500)
                      from bafMilli in Gen.Choose(10, 990) // l = max(b, 1 − b) ≤ 0.99 keeps every edge equation regular
                      select new OncologyAnalyzer.AlleleSpecificSegmentSummary(
                          "1", 0, 1000, logRMilli / 1000.0, bafMilli / 1000.0, 10)).ArrayOf(n)
        select (rhoPercent / 100.0, psiTwentieths / 20.0, segs);

    /// <summary>
    /// Battenberg two-state model: a clonal segment has one state at fraction 1 whose BAF level lies within maxdist 0.01
    /// of l = max(b, 1 − b); a sub-clonal segment has two non-negative integer states on one edge of the copy-number square
    /// (differing by exactly one copy of one allele), fractions τ and 1 − τ summing to 1, and the τ-mixture reproduces
    /// the segment BAF l (the equation Battenberg solves for τ).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitSubclonalCopyNumber_StatesSatisfyBattenbergTwoStateModel()
    {
        return Prop.ForAll(BattenbergGen().ToArbitrary(), g =>
        {
            var fits = OncologyAnalyzer.FitSubclonalCopyNumber(g.Segments, g.Rho, g.Psi);
            if (fits.Count != g.Segments.Length)
            {
                return false.Label("one fit per segment");
            }

            for (int i = 0; i < fits.Count; i++)
            {
                var f = fits[i];
                double l = Math.Max(g.Segments[i].MeanBAF, 1.0 - g.Segments[i].MeanBAF);
                var s1 = f.PrimaryState;
                if (s1.MajorCopyNumber < 0 || s1.MinorCopyNumber < 0 || !f.Segment.Equals(g.Segments[i]))
                {
                    return false.Label($"segment {i}: negative CN or segment not echoed");
                }

                if (!f.IsSubclonal)
                {
                    bool clonalOk = f.SecondaryState is null && s1.CellFraction == 1.0
                                    && Math.Abs(CornerLevel(g.Rho, s1.MajorCopyNumber, s1.MinorCopyNumber) - l) < 0.01;
                    if (!clonalOk)
                    {
                        return false.Label($"segment {i}: clonal call off the maxdist rule");
                    }

                    continue;
                }

                if (f.SecondaryState is not { } s2)
                {
                    return false.Label($"segment {i}: sub-clonal without a second state");
                }

                int dMajor = Math.Abs(s1.MajorCopyNumber - s2.MajorCopyNumber);
                int dMinor = Math.Abs(s1.MinorCopyNumber - s2.MinorCopyNumber);
                bool oneEdge = s2.MinorCopyNumber >= 0 && dMajor + dMinor == 1;
                double tau = s1.CellFraction;
                bool fractionsSum = Math.Abs(tau + s2.CellFraction - 1.0) <= 1e-12 * Math.Max(1.0, Math.Abs(tau));

                double numerator = 1.0 - g.Rho + g.Rho * (tau * s1.MajorCopyNumber + (1.0 - tau) * s2.MajorCopyNumber);
                double denominator = 2.0 - 2.0 * g.Rho
                                     + g.Rho * (tau * s1.TotalCopyNumber + (1.0 - tau) * s2.TotalCopyNumber);
                bool reproducesBaf = Math.Abs(numerator - l * denominator)
                                     <= 1e-9 * (Math.Abs(numerator) + Math.Abs(l * denominator) + 1.0);

                if (!(oneEdge && fractionsSum && reproducesBaf))
                {
                    return false.Label(
                        $"segment {i}: states {s1.MajorCopyNumber}:{s1.MinorCopyNumber}/{s2.MajorCopyNumber}:{s2.MinorCopyNumber}, " +
                        $"τ={tau}, edge={oneEdge}, sum={fractionsSum}, baf={reproducesBaf}");
                }
            }

            return true.Label("ok");
        });
    }

    // -------------------------------------------------------------------------
    // ascat.aspcf segmentation
    // -------------------------------------------------------------------------

    private static Gen<OncologyAnalyzer.AlleleSpecificLocus[]> LociGen() =>
        from chromosomes in Gen.Choose(1, 3)
        from perChromosome in Gen.Choose(1, 120).ArrayOf(chromosomes)
        from seed in Gen.Choose(0, int.MaxValue)
        select BuildLoci(perChromosome, seed);

    private static OncologyAnalyzer.AlleleSpecificLocus[] BuildLoci(int[] perChromosome, int seed)
    {
        var rng = new Random(seed);
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        for (int c = 0; c < perChromosome.Length; c++)
        {
            double level = rng.Next(-3, 4) * 0.3;
            double imbalance = rng.Next(0, 4) * 0.1;
            for (int i = 0; i < perChromosome[c]; i++)
            {
                if (rng.NextDouble() < 0.05)
                {
                    level = rng.Next(-3, 4) * 0.3;
                    imbalance = rng.Next(0, 4) * 0.1;
                }

                double baf = Math.Clamp(0.5 + (rng.Next(2) == 0 ? imbalance : -imbalance) + (rng.NextDouble() - 0.5) * 0.1, 0.0, 1.0);
                loci.Add(new OncologyAnalyzer.AlleleSpecificLocus(
                    (c + 1).ToString(), i * 1000L, level + (rng.NextDouble() - 0.5) * 0.4, baf));
            }
        }

        return loci.ToArray();
    }

    /// <summary>
    /// ascat.aspcf output contract: the segments tile the loci in input order without gaps or overlaps (Σ LocusCount = n),
    /// never span two chromosomes, carry the first/last locus positions as Start/End, report as logR level the mean of the
    /// raw logR of their loci, and a mirrored BAF 0.5 + μ ∈ [0.5, 1].
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property SegmentAlleleSpecificAspcf_TilesLociWithRawLogRMeans()
    {
        return Prop.ForAll(LociGen().ToArbitrary(), loci =>
        {
            var segments = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci);
            int cursor = 0;
            foreach (var s in segments)
            {
                if (s.LocusCount < 1 || cursor + s.LocusCount > loci.Length)
                {
                    return false.Label($"segment overruns the loci at {cursor}");
                }

                var run = loci.Skip(cursor).Take(s.LocusCount).ToArray();
                bool oneChromosome = run.All(l => l.Chromosome == s.Chromosome);
                bool bounds = s.Start == run[0].Position && s.End == run[^1].Position;
                double mean = 0.0;
                foreach (var l in run)
                {
                    mean += l.LogR;
                }

                mean /= run.Length;
                bool logROk = Math.Abs(s.MeanLogR - mean) <= 1e-12;
                bool bafOk = s.MeanBAF >= 0.5 && s.MeanBAF <= 1.0;
                if (!(oneChromosome && bounds && logROk && bafOk))
                {
                    return false.Label($"segment at {cursor}: chrom={oneChromosome}, bounds={bounds}, logR={logROk}, baf={bafOk}");
                }

                cursor += s.LocusCount;
            }

            return (cursor == loci.Length).Label($"Σ LocusCount {cursor} ≠ {loci.Length}");
        });
    }

    // -------------------------------------------------------------------------
    // LICHeE trunk (Werner 2017) and MATH (maftools)
    // -------------------------------------------------------------------------

    private static Arbitrary<(OncologyAnalyzer.CcfCluster[] Clusters, double Tolerance)> PhylogenyArbitrary() =>
        (from samples in Gen.Choose(1, 3)
         from n in Gen.Choose(1, 6)
         from vectors in Gen.Elements(0, 250, 500, 750, 900, 970, 1000).Select(v => v / 1000.0).ArrayOf(samples).ArrayOf(n)
         from tolMilli in Gen.Choose(0, 100)
         select (vectors.Select((v, i) => new OncologyAnalyzer.CcfCluster(i + 1, v)).ToArray(), tolMilli / 1000.0))
        .ToArbitrary();

    /// <summary>
    /// Werner et al. 2017: every trunk cluster has CCF ≥ 1 − ε in every sample, and the trunk is maximal — the last trunk
    /// node (or the root, for an empty trunk) has no child that is clonal in every sample. Oracle recomputed from the
    /// input CCFs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IdentifyTrunkMutations_IsTheMaximalAllClonalRootPath()
    {
        return Prop.ForAll(PhylogenyArbitrary(), p =>
        {
            if (!OncologyAnalyzer.TryReconstructPhylogeny(p.Clusters, out var phylo, p.Tolerance))
            {
                return true.Label("no valid tree");
            }

            var ccf = p.Clusters.ToDictionary(c => c.Id, c => c.CcfPerSample);
            bool Clonal(int id) => ccf[id].All(v => v >= 1.0 - p.Tolerance);

            var trunk = OncologyAnalyzer.IdentifyTrunkMutations(phylo);
            if (!trunk.All(Clonal))
            {
                return false.Label($"trunk [{string.Join(",", trunk)}] contains a cluster with CCF < 1 − ε");
            }

            int last = trunk.Count == 0 ? phylo.RootId : trunk[^1];
            bool maximal = !phylo.ChildrenOf(last).Any(Clonal);
            return maximal.Label($"trunk [{string.Join(",", trunk)}] stops before an all-clonal child of {last}");
        });
    }

    private static Arbitrary<(double[] Values, int Seed, int ScalePermille)> VafArbitrary() =>
        (from n in Gen.Choose(1, 25)
         from values in Gen.Choose(1, 1000).Select(v => v / 1000.0).ArrayOf(n)
         from seed in Gen.Choose(0, int.MaxValue)
         from scale in Gen.Choose(50, 1000)
         select (values, seed, scale)).ToArbitrary();

    /// <summary>
    /// MATH depends only on the multiset of VAFs (both medians are order statistics): any permutation gives a
    /// bit-identical score.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateITH_IsPermutationInvariant()
    {
        return Prop.ForAll(VafArbitrary(), t =>
        {
            var rng = new Random(t.Seed);
            double[] shuffled = t.Values.OrderBy(_ => rng.Next()).ToArray();
            return (OncologyAnalyzer.CalculateITH(shuffled) == OncologyAnalyzer.CalculateITH(t.Values))
                .Label("MATH changed under a permutation of the VAFs");
        });
    }

    /// <summary>
    /// MATH = 100·1.4826·MAD/median is scale-free: multiplying every VAF by c ∈ (0, 1] (so VAFs stay in [0, 1])
    /// scales MAD and median alike and leaves the score unchanged (up to rounding).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateITH_IsInvariantToScalingAllVafs()
    {
        return Prop.ForAll(VafArbitrary(), t =>
        {
            double c = t.ScalePermille / 1000.0;
            double baseMath = OncologyAnalyzer.CalculateITH(t.Values);
            double scaledMath = OncologyAnalyzer.CalculateITH(t.Values.Select(v => v * c).ToArray());
            return (Math.Abs(scaledMath - baseMath) <= 1e-9 * Math.Max(1.0, baseMath))
                .Label($"MATH {baseMath} → {scaledMath} under scaling by {c}");
        });
    }
}
