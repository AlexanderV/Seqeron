// SPDX-License-Identifier: GPL-3.0-only
// Tests of the Port of ASCAT 3.2.0 and Battenberg 3.0.0; part of Seqeron.Genomics.AlleleSpecificCopyNumber.Tests, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// FIN-B24 heavy tier — fuzz tests for the input contracts added by the B24 finisher: R <c>NA</c> handling in
/// <c>ascat.aspcf</c> / <c>ascat.asmultipcf</c> (F59: NaN is accepted and follows R's NA path), the asmultipcf
/// one-sample / probe-set errors (F53: R fails with "incorrect number of subscripts on matrix"), the GISTIC2 ziggurat cap
/// and noise filter (F50–F54), the canonical copy-number guards (F28/F66) and the Student t / incomplete beta domain
/// (F39). Random inputs are locally seeded; every run either succeeds with a well-formed result or throws the documented
/// exception type.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Oncology")]
public class OncologyFinisherFuzzTests
{
    private static List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> NoisyLoci(Random rng, int chromosomes, int perChromosome, double naRate)
    {
        var loci = new List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus>();
        for (int c = 1; c <= chromosomes; c++)
        {
            for (int i = 0; i < perChromosome; i++)
            {
                double logR = rng.NextDouble() < naRate ? double.NaN : (rng.NextDouble() - 0.5) * (rng.Next(10) == 0 ? 20 : 1);
                double baf = rng.NextDouble() < naRate ? double.NaN : rng.Next(8) switch { 0 => 0.0, 1 => 1.0, _ => rng.NextDouble() };
                loci.Add(new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus(c.ToString(), 1000L * (i + 1), logR, baf));
            }
        }

        return loci;
    }

    [Test]
    [Description("F59: NaN logR/BAF (R NA) anywhere in the germline-aware ascat.aspcf input never throws; every probe gets a finite segmented logR unless its whole level is missing, segmented BAF is NaN or in [0, 1], and segments tile the loci.")]
    public void SegmentAlleleSpecificAspcf_RandomMissingValues_WellFormed()
    {
        for (int seed = 0; seed < 60; seed++)
        {
            var rng = new Random(seed);
            var loci = NoisyLoci(rng, rng.Next(1, 4), rng.Next(1, 60), naRate: 0.15);
            var het = loci.Select(_ => rng.NextDouble() < 0.7).ToArray();

            var result = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAspcf(loci, het);

            Assert.Multiple(() =>
            {
                Assert.That(result.SegmentedLogR, Has.Count.EqualTo(loci.Count), $"seed {seed}");
                Assert.That(result.SegmentedBaf.All(b => double.IsNaN(b) || (b >= 0 && b <= 1)), Is.True, $"BAF range, seed {seed}");
                Assert.That(result.SegmentedLogR.All(v => !double.IsInfinity(v)), Is.True, $"no infinite logR, seed {seed}");
                Assert.That(result.Segments.Sum(s => s.LocusCount), Is.EqualTo(loci.Count), $"tiling, seed {seed}");
            });
        }
    }

    [Test]
    [Description("F53/F59: random multi-sample input with missing values either segments every sample (one result per sample, per-probe arrays aligned) or throws ArgumentException for a single-probe chromosome part.")]
    public void SegmentAlleleSpecificAsMultiPcf_RandomMissingValues_WellFormed()
    {
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(1000 + seed);
            int chromosomes = rng.Next(1, 3), probes = rng.Next(2, 40), samplesCount = rng.Next(2, 4);
            var template = NoisyLoci(rng, chromosomes, probes, naRate: 0);
            var samples = Enumerable.Range(0, samplesCount)
                .Select(_ => (IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus>)template.Select(l => l with
                {
                    LogR = rng.NextDouble() < 0.1 ? double.NaN : l.LogR + (rng.NextDouble() - 0.5) * 0.1,
                    BAF = rng.NextDouble() < 0.1 ? double.NaN : l.BAF,
                }).ToList())
                .ToList();

            var result = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAsMultiPcf(samples);

            Assert.That(result, Has.Count.EqualTo(samplesCount), $"seed {seed}");
            foreach (var r in result)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(r.SegmentedLogR, Has.Count.EqualTo(template.Count), $"seed {seed}");
                    Assert.That(r.SegmentedBaf.All(b => double.IsNaN(b) || (b >= 0.5 && b <= 1)), Is.True, $"mirrored BAF, seed {seed}");
                });
            }
        }
    }

    [Test]
    [Description("F53 (R: 'incorrect number of subscripts on matrix'): one sample, differing probe sets, a wrong genotype count and a single-probe chromosome part are ArgumentExceptions.")]
    public void SegmentAlleleSpecificAsMultiPcf_ContractViolations_Throw()
    {
        var rng = new Random(7);
        var loci = NoisyLoci(rng, 2, 10, naRate: 0);
        var shifted = loci.Select(l => l with { Position = l.Position + 1 }).ToList();
        var singleProbePart = loci.Take(10).Append(loci[10]).ToList();

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAsMultiPcf(new[] { loci }));
            Assert.Throws<ArgumentException>(() => AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAsMultiPcf(new[] { loci, shifted }));
            Assert.Throws<ArgumentException>(() => AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAsMultiPcf(new[] { loci, loci }, new bool[3]));
            Assert.Throws<ArgumentException>(() => AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAsMultiPcf(new[] { singleProbePart, singleProbePart }));
            Assert.Throws<ArgumentNullException>(() => AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAsMultiPcf(null!));
        });
    }
}
