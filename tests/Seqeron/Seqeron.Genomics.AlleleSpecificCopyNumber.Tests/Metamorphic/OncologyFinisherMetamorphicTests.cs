// SPDX-License-Identifier: GPL-3.0-only
// Tests of the Port of ASCAT 3.2.0 and Battenberg 3.0.0; part of Seqeron.Genomics.AlleleSpecificCopyNumber.Tests, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the production project).
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

    private static List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> RandomLoci(Random rng, int chromosomes, int perChromosome)
    {
        var loci = new List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus>();
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
                loci.Add(new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus(c.ToString(), 1000L * (i + 1), logR, baf));
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
            var a = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAspcf(loci, het);

            var altered = loci.Select((l, i) => het[i] ? l : l with { BAF = rng.Next(3) == 0 ? double.NaN : rng.NextDouble() }).ToList();
            var b = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAspcf(altered, het);

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
            var hetOnly = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAspcf(loci);
            var full = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAspcf(loci, Enumerable.Repeat(true, loci.Count).ToArray());

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
            var samples = new List<IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus>>();
            for (int s = 0; s < samplesCount; s++)
            {
                samples.Add(baseLoci.Select(l => l with
                {
                    LogR = l.LogR * (0.6 + 0.2 * s) + (rng.NextDouble() - 0.5) * 0.05,
                    BAF = Math.Clamp(l.BAF + (rng.NextDouble() - 0.5) * 0.05, 0, 1),
                }).ToList());
            }

            var het = baseLoci.Select(_ => rng.NextDouble() < 0.8).ToArray();
            var forward = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAsMultiPcf(samples, het);
            var reversed = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAsMultiPcf(samples.AsEnumerable().Reverse().ToList(), het);

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

    private static List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary> ModelGenomeWithSexChromosomes(Random rng, double rho, double psi)
    {
        var segs = new List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary>();
        string[] chrs = Enumerable.Range(1, 8).Select(i => i.ToString()).Concat(new[] { "X", "X", "Y" }).ToArray();
        for (int i = 0; i < chrs.Length; i++)
        {
            int major = rng.Next(1, 4), minor = rng.Next(0, major + 1);
            int total = major + minor;
            double logR = Math.Log2((2 * (1 - rho) + rho * total) / (2 * (1 - rho) + rho * psi)) + (rng.NextDouble() - 0.5) * 0.02;
            double baf = (1 - rho + rho * minor) / (2 * (1 - rho) + rho * total);
            segs.Add(new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary(chrs[i], i * 2_000_000L, i * 2_000_000L + 1_000_000, logR, baf, rng.Next(10, 200)));
        }

        return segs;
    }

    [Test]
    [Description("EQ/INV (runASCAT haploidchrs): the XX model equals the gender-less default; XY keeps ρ, GoF and autosomal segments and emits nB = 0 on X/Y.")]
    public void EvaluatePurityPloidy_SexModel_FemaleIsDefault_MaleHaploidSexChromosomes()
    {
        for (int seed = 0; seed < 20; seed++)
        {
            var rng = new Random(300 + seed);
            double rho = 0.4 + 0.05 * rng.Next(0, 12), psi = 1.8 + 0.1 * rng.Next(0, 20);
            var segs = ModelGenomeWithSexChromosomes(rng, rho, psi);

            var dflt = AlleleSpecificCopyNumberAnalyzer.EvaluatePurityPloidy(segs, rho, psi);
            var female = AlleleSpecificCopyNumberAnalyzer.EvaluatePurityPloidy(segs, AlleleSpecificCopyNumberAnalyzer.AscatSexModel.Female, rho, psi);
            var male = AlleleSpecificCopyNumberAnalyzer.EvaluatePurityPloidy(segs, AlleleSpecificCopyNumberAnalyzer.AscatSexModel.Male, rho, psi);

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

    private static List<AlleleSpecificCopyNumberAnalyzer.SubclonalSegmentSnpBafs> RandomSnpSegments(Random rng, double rho, double psi, int n, bool constantSnps)
    {
        var list = new List<AlleleSpecificCopyNumberAnalyzer.SubclonalSegmentSnpBafs>();
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
            var summary = new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary(
                (1 + i % 22).ToString(), i * 2_000_000L, i * 2_000_000L + 1_000_000, logR, baf, snps);
            list.Add(new AlleleSpecificCopyNumberAnalyzer.SubclonalSegmentSnpBafs(summary, snpBafs));
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
            var plain = AlleleSpecificCopyNumberAnalyzer.FitSubclonalCopyNumber(constant.Select(s => s.Segment).ToList(), rho, psi);
            var tested = AlleleSpecificCopyNumberAnalyzer.FitSubclonalCopyNumberWithSnpTest(constant, rho, psi);
            Assert.That(tested.Select(t => t.Fit), Is.EqualTo(plain), $"constant SNP BAFs, seed {seed}");

            var noisy = RandomSnpSegments(rng, rho, psi, 8, constantSnps: false);
            var shortcut = AlleleSpecificCopyNumberAnalyzer.FitSubclonalCopyNumberWithSnpTest(noisy, rho, psi, maxBafDistance: 1.0);
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
            var a = AlleleSpecificCopyNumberAnalyzer.FitSubclonalCopyNumberWithBootstrap(segs, rho, psi, seed: 17 + seed, permutations: 200);
            var b = AlleleSpecificCopyNumberAnalyzer.FitSubclonalCopyNumberWithBootstrap(segs, rho, psi, seed: 17 + seed, permutations: 200);
            var tested = AlleleSpecificCopyNumberAnalyzer.FitSubclonalCopyNumberWithSnpTest(segs, rho, psi);

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
            var snps = new List<AlleleSpecificCopyNumberAnalyzer.PhasedBafSnp>();
            for (int c = 1; c <= 2; c++)
            {
                int n = rng.Next(60, 160);
                double level = 0.5 + 0.1 * rng.Next(0, 4);
                for (int i = 0; i < n; i++)
                {
                    double b = rng.Next(2) == 0 ? level : 1 - level;
                    b = Math.Clamp(b + (rng.NextDouble() - 0.5) * 0.1, 0, 1);
                    snps.Add(new AlleleSpecificCopyNumberAnalyzer.PhasedBafSnp(c.ToString(), 1000L * (i + 1), rng.Next(25) == 0 ? double.NaN : b));
                }
            }

            var a = AlleleSpecificCopyNumberAnalyzer.SegmentPhasedBaf(snps);
            var again = AlleleSpecificCopyNumberAnalyzer.SegmentPhasedBaf(snps);
            Assert.Multiple(() =>
            {
                Assert.That(again, Is.EqualTo(a), $"deterministic, seed {seed}");
                Assert.That(a.Count, Is.EqualTo(snps.Count(s => !double.IsNaN(s.Baf))), $"rows, seed {seed}");
                Assert.That(a.All(r => r.BafPhased == r.Baf || r.BafPhased == 1 - r.Baf), Is.True, $"phased, seed {seed}");
                Assert.That(a.All(r => r.BafSegment >= 0 && r.BafSegment <= 1), Is.True, $"segment BAF range, seed {seed}");
            });
        }
    }


    private static bool IsSex(string chr) => chr is "X" or "Y";
}
