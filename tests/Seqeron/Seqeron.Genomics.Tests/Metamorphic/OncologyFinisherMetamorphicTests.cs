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

    private static bool IsSex(string chr) => chr is "X" or "Y";

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
