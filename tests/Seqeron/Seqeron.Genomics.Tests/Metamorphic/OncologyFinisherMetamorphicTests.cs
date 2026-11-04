namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// FIN-B24 heavy tier — metamorphic relations for the behaviour added by the B24 finisher (F24–F66). Each relation follows
/// from the ported reference code (ASCAT, Battenberg, CNAqc, GISTIC2, Ckmeans.1d.dp, LICHeE, maftools), not from
/// observed output. All randomness is locally seeded.
/// <list type="bullet">
/// <item>ascat.aspcf with germline genotypes (F36/F59): the BAF of germline-homozygous loci is never read (ASCAT sets it
/// to NA); with every locus heterozygous the segmentation equals the heterozygous-only overload on one chromosome.</item>
/// <item>ascat.asmultipcf (F53): the joint cost Σ_tracks is symmetric in the samples, so permuting the samples permutes
/// the per-sample results.</item>
/// <item>runASCAT sex model (F38): XX is the default; under XY the X/Y segments are haploid (nB = 0) and the fit, which
/// uses autosomes only, is unchanged.</item>
/// <item>Battenberg <c>determine_copynumber</c> (F39/F41): maxdist ≥ 1 forces pval = 1 (clonal); constant SNP BAFs reproduce
/// the plain fit; the bootstrap is a deterministic function of the seed and does not change the fitted states.</item>
/// <item>Battenberg <c>segment.baf.phased</c> (F40): BAFphased is BAF or 1 − BAF; NaN BAFs are dropped.</item>
/// <item>ascat.metrics (F30): nMajor/nMinor are max/min (label-free); sex chromosomes are excluded.</item>
/// <item>ASCAT probe-weighted ploidy (F31): equal probe counts ≡ the bp-weighted Patchwork mean on equal-length segments.</item>
/// <item>UCSC cytoband arm table (F55): p-arm + q-arm = chromosome length; the centromere lies inside; lengths agree with
/// the autosome-length table.</item>
/// <item>GISTIC2 ziggurat (F50–F54): the four event kinds reconstruct each sample's capped copy number − 2 per marker
/// (amp/aod add, del/doa subtract), and sample order only relabels the events.</item>
/// </list>
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Oncology")]
public class OncologyFinisherMetamorphicTests
{
    // -------------------------------------------------------------------------
    // ascat.aspcf with germline genotypes (F36, F59)
    // -------------------------------------------------------------------------

    private static List<OncologyAnalyzer.AlleleSpecificLocus> RandomLoci(Random rng, int chromosomes, int perChromosome)
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        for (int c = 1; c <= chromosomes; c++)
        {
            int blocks = rng.Next(1, 4);
            int block = 0;
            for (int i = 0; i < perChromosome; i++)
            {
                if (i == perChromosome * (block + 1) / blocks)
                {
                    block++;
                }

                double level = 0.15 * ((c + block) % 4) - 0.2 + 0.013;
                double bafLevel = 0.1 * ((c * 3 + block) % 4);
                double logR = level + (rng.NextDouble() - 0.5) * 0.2;
                double baf = Math.Clamp(0.5 + (rng.Next(2) == 0 ? -1 : 1) * bafLevel + (rng.NextDouble() - 0.5) * 0.1, 0, 1);
                loci.Add(new OncologyAnalyzer.AlleleSpecificLocus(c.ToString(), 1000L * (i + 1), logR, baf));
            }
        }

        return loci;
    }

    [Test]
    [Description("INV (ascat.aspcf: Tumor_BAF[germline homozygous] <- NA): changing the BAF of homozygous loci (including to NaN) leaves the segmentation bit-identical.")]
    public void SegmentAlleleSpecificAspcf_HomozygousLocusBaf_IsNeverRead()
    {
        for (int seed = 0; seed < 12; seed++)
        {
            var rng = new Random(seed);
            var loci = RandomLoci(rng, 2, 60);
            var het = loci.Select(_ => rng.NextDouble() < 0.7).ToArray();
            var a = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, het);

            var altered = loci.Select((l, i) => het[i] ? l : l with { BAF = rng.Next(3) == 0 ? double.NaN : rng.NextDouble() }).ToList();
            var b = OncologyAnalyzer.SegmentAlleleSpecificAspcf(altered, het);

            Assert.Multiple(() =>
            {
                Assert.That(b.SegmentedLogR, Is.EqualTo(a.SegmentedLogR), $"logR, seed {seed}");
                Assert.That(b.SegmentedBaf, Is.EqualTo(a.SegmentedBaf), $"BAF, seed {seed}");
                Assert.That(b.Segments, Is.EqualTo(a.Segments), $"segments, seed {seed}");
            });
        }
    }

    [Test]
    [Description("EQ (ascat.aspcf docs: with every locus heterozygous the levels/BAF are those of the het-only overload, on one chromosome with no zero level): all-het germline-aware ≡ het-only segmentation.")]
    public void SegmentAlleleSpecificAspcf_AllHeterozygous_MatchesHetOnlyOverload()
    {
        for (int seed = 0; seed < 12; seed++)
        {
            var rng = new Random(100 + seed);
            var loci = RandomLoci(rng, 1, 80);
            var hetOnly = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci);
            var full = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, Enumerable.Repeat(true, loci.Count).ToArray());

            Assert.That(full.Segments.Count, Is.EqualTo(hetOnly.Count), $"segment count, seed {seed}");
            for (int s = 0; s < hetOnly.Count; s++)
            {
                var x = hetOnly[s];
                var y = full.Segments[s];
                Assert.Multiple(() =>
                {
                    Assert.That((y.Chromosome, y.Start, y.End, y.LocusCount), Is.EqualTo((x.Chromosome, x.Start, x.End, x.LocusCount)), $"seed {seed} seg {s}");
                    Assert.That(y.MeanLogR, Is.EqualTo(x.MeanLogR).Within(1e-12), $"logR seed {seed} seg {s}");
                    Assert.That(y.MeanBAF, Is.EqualTo(x.MeanBAF).Within(1e-12), $"BAF seed {seed} seg {s}");
                });
            }
        }
    }

    // -------------------------------------------------------------------------
    // ascat.asmultipcf (F53)
    // -------------------------------------------------------------------------

    [Test]
    [Description("INV (ascat.asmultipcf joint cost Σ over tracks is symmetric in the samples): reversing the sample order reverses the per-sample segmentations.")]
    public void SegmentAlleleSpecificAsMultiPcf_SampleOrderPermutation_PermutesResults()
    {
        for (int seed = 0; seed < 8; seed++)
        {
            var rng = new Random(200 + seed);
            int samplesCount = rng.Next(2, 4);
            var baseLoci = RandomLoci(rng, 2, 40);
            var samples = new List<IReadOnlyList<OncologyAnalyzer.AlleleSpecificLocus>>();
            for (int s = 0; s < samplesCount; s++)
            {
                samples.Add(baseLoci.Select(l => l with
                {
                    LogR = l.LogR * (0.6 + 0.2 * s) + (rng.NextDouble() - 0.5) * 0.05,
                    BAF = Math.Clamp(l.BAF + (rng.NextDouble() - 0.5) * 0.05, 0, 1),
                }).ToList());
            }

            var het = baseLoci.Select(_ => rng.NextDouble() < 0.8).ToArray();
            var forward = OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(samples, het);
            var reversed = OncologyAnalyzer.SegmentAlleleSpecificAsMultiPcf(samples.AsEnumerable().Reverse().ToList(), het);

            for (int s = 0; s < samplesCount; s++)
            {
                var a = forward[s];
                var b = reversed[samplesCount - 1 - s];
                Assert.That(b.Segments.Select(g => (g.Chromosome, g.Start, g.End, g.LocusCount)),
                    Is.EqualTo(a.Segments.Select(g => (g.Chromosome, g.Start, g.End, g.LocusCount))), $"breakpoints, seed {seed} sample {s}");
                for (int i = 0; i < a.SegmentedLogR.Count; i++)
                {
                    Assert.That(b.SegmentedLogR[i], Is.EqualTo(a.SegmentedLogR[i]).Within(1e-12), $"logR seed {seed} sample {s} probe {i}");
                    Assert.That(double.IsNaN(b.SegmentedBaf[i]) ? -1 : b.SegmentedBaf[i],
                        Is.EqualTo(double.IsNaN(a.SegmentedBaf[i]) ? -1 : a.SegmentedBaf[i]).Within(1e-12), $"BAF seed {seed} sample {s} probe {i}");
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // runASCAT sex model (F38)
    // -------------------------------------------------------------------------

    private static List<OncologyAnalyzer.AlleleSpecificSegmentSummary> ModelGenomeWithSexChromosomes(Random rng, double rho, double psi)
    {
        var segs = new List<OncologyAnalyzer.AlleleSpecificSegmentSummary>();
        string[] chrs = Enumerable.Range(1, 8).Select(i => i.ToString()).Concat(new[] { "X", "X", "Y" }).ToArray();
        for (int i = 0; i < chrs.Length; i++)
        {
            int major = rng.Next(1, 4), minor = rng.Next(0, major + 1);
            int total = major + minor;
            double logR = Math.Log2((2 * (1 - rho) + rho * total) / (2 * (1 - rho) + rho * psi)) + (rng.NextDouble() - 0.5) * 0.02;
            double baf = (1 - rho + rho * minor) / (2 * (1 - rho) + rho * total);
            segs.Add(new OncologyAnalyzer.AlleleSpecificSegmentSummary(chrs[i], i * 2_000_000L, i * 2_000_000L + 1_000_000, logR, baf, rng.Next(10, 200)));
        }

        return segs;
    }

    private static bool IsSex(string chr) => chr is "X" or "Y";

    [Test]
    [Description("EQ/INV (runASCAT haploidchrs): the XX model equals the gender-less default; XY keeps ρ, GoF and autosomal segments and emits nB = 0 on X/Y.")]
    public void EvaluatePurityPloidy_SexModel_FemaleIsDefault_MaleHaploidSexChromosomes()
    {
        for (int seed = 0; seed < 20; seed++)
        {
            var rng = new Random(300 + seed);
            double rho = 0.4 + 0.05 * rng.Next(0, 12), psi = 1.8 + 0.1 * rng.Next(0, 20);
            var segs = ModelGenomeWithSexChromosomes(rng, rho, psi);

            var dflt = OncologyAnalyzer.EvaluatePurityPloidy(segs, rho, psi);
            var female = OncologyAnalyzer.EvaluatePurityPloidy(segs, OncologyAnalyzer.AscatSexModel.Female, rho, psi);
            var male = OncologyAnalyzer.EvaluatePurityPloidy(segs, OncologyAnalyzer.AscatSexModel.Male, rho, psi);

            Assert.Multiple(() =>
            {
                Assert.That((female.Purity, female.Ploidy, female.GoodnessOfFit), Is.EqualTo((dflt.Purity, dflt.Ploidy, dflt.GoodnessOfFit)), $"seed {seed}");
                Assert.That(female.Segments, Is.EqualTo(dflt.Segments), $"female segments, seed {seed}");
                Assert.That((male.Purity, male.GoodnessOfFit), Is.EqualTo((female.Purity, female.GoodnessOfFit)), $"male fit, seed {seed}");
                Assert.That(male.Segments.Where(s => !IsSex(s.Chromosome)), Is.EqualTo(female.Segments.Where(s => !IsSex(s.Chromosome))), $"autosomes, seed {seed}");
                Assert.That(male.Segments.Where(s => IsSex(s.Chromosome)).Select(s => Math.Min(s.MajorCopyNumber, s.MinorCopyNumber)),
                    Is.All.EqualTo(0), $"haploid X/Y nB = 0, seed {seed}");
            });
        }
    }

    // -------------------------------------------------------------------------
    // Battenberg determine_copynumber (F39, F41) and segment.baf.phased (F40)
    // -------------------------------------------------------------------------

    private static List<OncologyAnalyzer.SubclonalSegmentSnpBafs> RandomSnpSegments(Random rng, double rho, double psi, int n, bool constantSnps)
    {
        var list = new List<OncologyAnalyzer.SubclonalSegmentSnpBafs>();
        for (int i = 0; i < n; i++)
        {
            int major = rng.Next(1, 4), minor = rng.Next(0, major + 1);
            double tau = rng.Next(2) == 0 ? 1.0 : 0.3 + 0.4 * rng.NextDouble();
            int major2 = major + 1;
            double total = tau * (major + minor) + (1 - tau) * (major2 + minor);
            double mtot = tau * major + (1 - tau) * major2;
            double logR = Math.Log2((2 * (1 - rho) + rho * total) / (2 * (1 - rho) + rho * psi));
            double baf = (1 - rho + rho * mtot) / (2 * (1 - rho) + rho * total);
            int snps = rng.Next(3, 30);
            var snpBafs = Enumerable.Range(0, snps)
                .Select(_ => constantSnps ? baf : Math.Clamp(baf + (rng.NextDouble() - 0.5) * 0.08, 0, 1)).ToArray();
            var summary = new OncologyAnalyzer.AlleleSpecificSegmentSummary(
                (1 + i % 22).ToString(), i * 2_000_000L, i * 2_000_000L + 1_000_000, logR, baf, snps);
            list.Add(new OncologyAnalyzer.SubclonalSegmentSnpBafs(summary, snpBafs));
        }

        return list;
    }

    [Test]
    [Description("EQ (determine_copynumber: sd(BAFke) == 0 ⇒ pval 0 ⇔ the plain distance rule): constant SNP BAFs reproduce FitSubclonalCopyNumber; maxdist = 1 ⇒ every pval = 1 and every segment clonal.")]
    public void FitSubclonalCopyNumberWithSnpTest_ConstantSnpsAndMaxDistShortcut()
    {
        for (int seed = 0; seed < 25; seed++)
        {
            var rng = new Random(400 + seed);
            double rho = 0.5 + 0.05 * rng.Next(0, 10), psi = 1.9 + 0.1 * rng.Next(0, 15);
            var constant = RandomSnpSegments(rng, rho, psi, 8, constantSnps: true);
            var plain = OncologyAnalyzer.FitSubclonalCopyNumber(constant.Select(s => s.Segment).ToList(), rho, psi);
            var tested = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(constant, rho, psi);
            Assert.That(tested.Select(t => t.Fit), Is.EqualTo(plain), $"constant SNP BAFs, seed {seed}");

            var noisy = RandomSnpSegments(rng, rho, psi, 8, constantSnps: false);
            var shortcut = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(noisy, rho, psi, maxBafDistance: 1.0);
            Assert.Multiple(() =>
            {
                Assert.That(shortcut.Select(t => t.PValue), Is.All.EqualTo(1.0), $"maxdist pval, seed {seed}");
                Assert.That(shortcut.Select(t => t.Fit.IsSubclonal), Is.All.False, $"maxdist clonal, seed {seed}");
            });
        }
    }

    [Test]
    [Description("DET (R set.seed stream): the Battenberg bootstrap is reproducible for a seed, and its fitted states equal the t-test fit (the bootstrap only adds solutions / CIs).")]
    public void FitSubclonalCopyNumberWithBootstrap_SeedReproducible_FitUnchanged()
    {
        for (int seed = 0; seed < 6; seed++)
        {
            var rng = new Random(500 + seed);
            double rho = 0.5 + 0.05 * rng.Next(0, 10), psi = 1.9 + 0.1 * rng.Next(0, 15);
            var segs = RandomSnpSegments(rng, rho, psi, 6, constantSnps: false);
            var a = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(segs, rho, psi, seed: 17 + seed, permutations: 200);
            var b = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(segs, rho, psi, seed: 17 + seed, permutations: 200);
            var tested = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(segs, rho, psi);

            for (int i = 0; i < segs.Count; i++)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(b[i].Solutions, Is.EqualTo(a[i].Solutions), $"seed {seed} segment {i}");
                    Assert.That(a[i].Fit, Is.EqualTo(tested[i].Fit), $"fit, seed {seed} segment {i}");
                    Assert.That(a[i].PValue, Is.EqualTo(tested[i].PValue), $"pval, seed {seed} segment {i}");
                });
            }
        }
    }

    [Test]
    [Description("DEF (segment.baf.phased: BAFphased = ifelse(BAFsegm > 0.5, BAF, 1 − BAF); NA BAFs dropped): one row per non-NaN SNP, BAFphased ∈ {BAF, 1 − BAF}, segment BAF ∈ [0, 1]; deterministic.")]
    public void SegmentPhasedBaf_PhasedIsBafOrMirror_NaNDropped_Deterministic()
    {
        for (int seed = 0; seed < 10; seed++)
        {
            var rng = new Random(600 + seed);
            var snps = new List<OncologyAnalyzer.PhasedBafSnp>();
            for (int c = 1; c <= 2; c++)
            {
                int n = rng.Next(60, 160);
                double level = 0.5 + 0.1 * rng.Next(0, 4);
                for (int i = 0; i < n; i++)
                {
                    double b = rng.Next(2) == 0 ? level : 1 - level;
                    b = Math.Clamp(b + (rng.NextDouble() - 0.5) * 0.1, 0, 1);
                    snps.Add(new OncologyAnalyzer.PhasedBafSnp(c.ToString(), 1000L * (i + 1), rng.Next(25) == 0 ? double.NaN : b));
                }
            }

            var a = OncologyAnalyzer.SegmentPhasedBaf(snps);
            var again = OncologyAnalyzer.SegmentPhasedBaf(snps);
            Assert.Multiple(() =>
            {
                Assert.That(again, Is.EqualTo(a), $"deterministic, seed {seed}");
                Assert.That(a.Count, Is.EqualTo(snps.Count(s => !double.IsNaN(s.Baf))), $"rows, seed {seed}");
                Assert.That(a.All(r => r.BafPhased == r.Baf || r.BafPhased == 1 - r.Baf), Is.True, $"phased, seed {seed}");
                Assert.That(a.All(r => r.BafSegment >= 0 && r.BafSegment <= 1), Is.True, $"segment BAF range, seed {seed}");
            });
        }
    }

    // -------------------------------------------------------------------------
    // ascat.metrics (F30), probe-weighted ploidy (F31)
    // -------------------------------------------------------------------------

    private static List<OncologyAnalyzer.AlleleSpecificSegment> RandomIntegerGenome(Random rng, int n)
    {
        var segs = new List<OncologyAnalyzer.AlleleSpecificSegment>();
        for (int i = 0; i < n; i++)
        {
            long start = i * 10_000_000L;
            long len = 1_000_000L * rng.Next(1, 9); // whole-Mb lengths ⇒ exact weight sums in any order
            segs.Add(new OncologyAnalyzer.AlleleSpecificSegment((1 + rng.Next(22)).ToString(), start, start + len, rng.Next(0, 6), rng.Next(0, 4)));
        }

        return segs;
    }

    [Test]
    [Description("INV (ascat.metrics: nMajor ≥ nMinor; autosomes only; tapply sums): allele-label swap, appended X/Y segments and segment permutation leave every metric unchanged; GI/LOH ∈ [0, 1]; WGD follows mode_majA.")]
    public void ComputeAscatGenomeMetrics_LabelSwapSexChromosomesPermutation_AreInvariant()
    {
        for (int seed = 0; seed < 60; seed++)
        {
            var rng = new Random(700 + seed);
            var segs = RandomIntegerGenome(rng, rng.Next(1, 25));
            var m = OncologyAnalyzer.ComputeAscatGenomeMetrics(segs);

            var swapped = segs.Select(s => s with { MajorCopyNumber = s.MinorCopyNumber, MinorCopyNumber = s.MajorCopyNumber }).ToList();
            var withSex = segs.Concat(new[]
            {
                new OncologyAnalyzer.AlleleSpecificSegment("X", 0, 150_000_000, 4, 0),
                new OncologyAnalyzer.AlleleSpecificSegment("chrY", 0, 50_000_000, 3, 3),
            }).ToList();
            var permuted = segs.OrderBy(_ => rng.Next()).ToList();

            var expectedWgd = m.ModeMajorAllele switch
            {
                1 => OncologyAnalyzer.AscatWgdStatus.NoWgd,
                2 => OncologyAnalyzer.AscatWgdStatus.Wgd,
                >= 3 and <= 5 => OncologyAnalyzer.AscatWgdStatus.WgdPlus,
                _ => OncologyAnalyzer.AscatWgdStatus.NotAvailable,
            };

            Assert.Multiple(() =>
            {
                Assert.That(OncologyAnalyzer.ComputeAscatGenomeMetrics(swapped), Is.EqualTo(m), $"label swap, seed {seed}");
                Assert.That(OncologyAnalyzer.ComputeAscatGenomeMetrics(withSex), Is.EqualTo(m), $"X/Y ignored, seed {seed}");
                Assert.That(OncologyAnalyzer.ComputeAscatGenomeMetrics(permuted), Is.EqualTo(m), $"permutation, seed {seed}");
                Assert.That(m.WgdStatus, Is.EqualTo(expectedWgd), $"WGD from mode_majA, seed {seed}");
                Assert.That(m.LossOfHeterozygosity, Is.InRange(0.0, 1.0));
                Assert.That(m.GenomicInstability is null == (expectedWgd == OncologyAnalyzer.AscatWgdStatus.NotAvailable), Is.True, $"GI NA iff WGD NA, seed {seed}");
                Assert.That(m.GenomicInstability ?? 0.5, Is.InRange(0.0, 1.0));
                Assert.That(m.ModeMinorAllele, Is.LessThanOrEqualTo(Math.Max(m.ModeMajorAllele, 5)));
            });
        }
    }

    [Test]
    [Description("EQ (ASCAT ψ = Σ CN·n / Σ n vs Patchwork Σ CN·L / Σ L): with equal segment lengths and equal probe counts the two ploidies coincide; multiplying every probe count by k leaves ψ unchanged.")]
    public void EstimatePloidy_ProbeWeighted_UniformProbesEqualLengths_MatchesBpWeighted()
    {
        for (int seed = 0; seed < 60; seed++)
        {
            var rng = new Random(800 + seed);
            int n = rng.Next(1, 30);
            var segs = Enumerable.Range(0, n)
                .Select(i => new OncologyAnalyzer.AlleleSpecificSegment((1 + i % 22).ToString(), i * 5_000_000L, i * 5_000_000L + 2_000_000, rng.Next(0, 6), rng.Next(0, 4)))
                .ToList();
            int probes = rng.Next(1, 500);
            var counts = Enumerable.Range(0, n).Select(_ => rng.Next(1, 300)).ToList();
            int k = rng.Next(2, 9);

            Assert.Multiple(() =>
            {
                Assert.That(OncologyAnalyzer.EstimatePloidy(segs, Enumerable.Repeat(probes, n)),
                    Is.EqualTo(OncologyAnalyzer.EstimatePloidy(segs)).Within(1e-12), $"uniform probes, seed {seed}");
                Assert.That(OncologyAnalyzer.EstimatePloidy(segs, counts.Select(c => c * k)),
                    Is.EqualTo(OncologyAnalyzer.EstimatePloidy(segs, counts)).Within(1e-12), $"probe scaling, seed {seed}");
            });
        }
    }

    // -------------------------------------------------------------------------
    // UCSC cytoband arm table (F55)
    // -------------------------------------------------------------------------

    [TestCase(OncologyAnalyzer.ReferenceGenome.GRCh38)]
    [TestCase(OncologyAnalyzer.ReferenceGenome.GRCh37)]
    [Description("DEF (UCSC cytoBand, GISTIC2 arm split): 0 < CentromereStart ≤ PArmEnd ≤ CentromereEnd ≤ Length, p + q = Length, and the autosome lengths equal GetAutosomeLengths (chrom.sizes).")]
    public void GetChromosomeArmLengths_ArmBookkeeping_IsConsistent(OncologyAnalyzer.ReferenceGenome genome)
    {
        var arms = OncologyAnalyzer.GetChromosomeArmLengths(genome);
        var autosomes = OncologyAnalyzer.GetAutosomeLengths(genome);

        Assert.That(arms.Select(a => a.Chromosome),
            Is.EqualTo(Enumerable.Range(1, 22).Select(i => i.ToString()).Concat(new[] { "X", "Y" })));
        foreach (var a in arms)
        {
            Assert.Multiple(() =>
            {
                Assert.That(a.CentromereStart, Is.GreaterThanOrEqualTo(0), a.Chromosome);
                Assert.That(a.CentromereStart, Is.LessThanOrEqualTo(a.PArmEnd), a.Chromosome);
                Assert.That(a.PArmEnd, Is.LessThanOrEqualTo(a.CentromereEnd), a.Chromosome);
                Assert.That(a.CentromereEnd, Is.LessThan(a.Length), a.Chromosome);
                Assert.That(a.PArmLength + a.QArmLength, Is.EqualTo(a.Length), a.Chromosome);
            });
        }

        Assert.That(arms.Take(22).Select(a => a.Length), Is.EqualTo(autosomes));
        Assert.That(arms.Take(22).Sum(a => a.Length), Is.EqualTo(OncologyAnalyzer.GetAutosomalGenomeLength(genome)));
    }

    // -------------------------------------------------------------------------
    // GISTIC2 ziggurat deconstruction (F50–F54)
    // -------------------------------------------------------------------------

    private static (List<OncologyAnalyzer.ZigguratChromosome> Layout, List<IReadOnlyList<OncologyAnalyzer.ZigguratSegment>> Samples) RandomCohort(Random rng)
    {
        var layout = new List<OncologyAnalyzer.ZigguratChromosome>();
        int chromosomes = rng.Next(1, 4);
        for (int c = 1; c <= chromosomes; c++)
        {
            layout.Add(new OncologyAnalyzer.ZigguratChromosome(c.ToString(), rng.Next(0, 15), rng.Next(1, 15)));
        }

        var samples = new List<IReadOnlyList<OncologyAnalyzer.ZigguratSegment>>();
        int n = rng.Next(1, 5);
        for (int s = 0; s < n; s++)
        {
            var segs = new List<OncologyAnalyzer.ZigguratSegment>();
            foreach (var chr in layout)
            {
                int markers = chr.PArmMarkerCount + chr.QArmMarkerCount;
                int pos = 1;
                while (pos <= markers)
                {
                    int end = Math.Min(markers, pos + rng.Next(0, 8));
                    double value = rng.Next(4) == 0 ? 0.0 : Math.Round((rng.NextDouble() - 0.5) * 3.4, 2);
                    segs.Add(new OncologyAnalyzer.ZigguratSegment(chr.Chromosome, pos, end, value));
                    pos = end + 1;
                }
            }

            samples.Add(segs);
        }

        return (layout, samples);
    }

    private static double Cn2(double log2) => Math.Pow(2, Math.Clamp(log2, -1.5, 1.5) + 1) - 2; // GISTIC2 cap 1.5, then 2^(x+1) − 2

    [Test]
    [Description("DEF (GISTIC2 make_final_Qs / reconstruct_genomes): per sample and marker, Σ amp + Σ aod − Σ del − Σ doa over the events covering it equals the capped copy number − 2 of the input segment.")]
    public void DeconstructZiggurat_EventsReconstructEverySample()
    {
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(900 + seed);
            var (layout, samples) = RandomCohort(rng);
            var result = OncologyAnalyzer.DeconstructZiggurat(layout, samples);

            for (int s = 0; s < samples.Count; s++)
            {
                foreach (var seg in samples[s])
                {
                    for (int marker = seg.StartMarker; marker <= seg.EndMarker; marker++)
                    {
                        double sum = result.Events
                            .Where(e => e.Sample == s && e.Chromosome == seg.Chromosome && e.StartMarker <= marker && marker <= e.EndMarker)
                            .Sum(e => e.Type is OncologyAnalyzer.ZigguratEventType.Amplification or OncologyAnalyzer.ZigguratEventType.AmplificationOverDeletion
                                ? e.Amplitude
                                : -e.Amplitude);
                        Assert.That(sum, Is.EqualTo(Cn2(seg.Value)).Within(1e-9),
                            $"seed {seed} sample {s} chr {seg.Chromosome} marker {marker}");
                    }
                }
            }
        }
    }

    [Test]
    [Description("INV (the cohort length × amplitude table is a sum over samples): reversing the sample order maps every event to the mirrored sample index.")]
    public void DeconstructZiggurat_SampleOrderPermutation_RelabelsEvents()
    {
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(1000 + seed);
            var (layout, samples) = RandomCohort(rng);
            int n = samples.Count;
            var forward = OncologyAnalyzer.DeconstructZiggurat(layout, samples);
            var reversed = OncologyAnalyzer.DeconstructZiggurat(layout, samples.AsEnumerable().Reverse().ToList());

            static string Key(OncologyAnalyzer.ZigguratEvent e) =>
                string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{e.Type}|{e.Chromosome}|{e.StartMarker}|{e.EndMarker}|{e.Amplitude:R}|{e.StartLevel:R}|{e.EndLevel:R}|{e.ArmFraction:R}|{e.ArmLevel:R}");

            for (int s = 0; s < n; s++)
            {
                var a = forward.ForSample(s).Select(Key).OrderBy(k => k, StringComparer.Ordinal);
                var b = reversed.ForSample(n - 1 - s).Select(Key).OrderBy(k => k, StringComparer.Ordinal);
                Assert.That(b, Is.EqualTo(a), $"seed {seed} sample {s}");
            }
        }
    }
}
