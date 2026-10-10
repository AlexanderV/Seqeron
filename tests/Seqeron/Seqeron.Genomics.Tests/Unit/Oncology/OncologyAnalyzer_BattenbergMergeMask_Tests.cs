// ONCO-ASCAT-001 — Battenberg merge_segments and mask_high_cn_segments (B24 F60)
// Evidence: docs/Evidence/ONCO-ASCAT-001-Evidence.md (§ F60)
// TestSpec: tests/TestSpecs/ONCO-ASCAT-001.md (§19)
// Source: Wedge-lab/battenberg master (57a8f7e) R/fitcopynumber.R merge_segments (calc_nmin/calc_nmaj, updateNeighbour,
//         updateAround, checkStatus, merge_seg; distance > 3e6, same clonal solution, same square + Welch t.test p < 0.05
//         on logR and BAFphased with > 10 values each) and mask_high_cn_segments; GenomicRanges 1.54.1 distance/width,
//         GenomeInfoDb rankSeqlevels (seqlevel order of makeGRangesFromDataFrame).
// Inputs and R expectations: BattenbergCallSubclonesData (G1–G4). Branches hit in R (verbose merge_segments):
//   G1 (opt 2 mean): too few values ×2, same clonal merge, no-significant-difference merge, significant ×2, different
//      squares, distance > 3 Mb; max_allowed_state 5 masks the (6,1) segment (39 of its 40 SNPs: half-open startpos < Position).
//   G2 (opt 3): chromosomes "10","2","X" reordered to 2, 10, X; cascade of merges (3 → 1 segments), LOH merge via t-test.
//   G3 (opt 1 median): chr3/chr1/chrX → chr1, chr3, chrX; t-test merge of two sub-clonal segments; distance > 3 Mb;
//      BAFseg 1 segment (nMaj 1001 from cn_upper_limit) masked at the default 250.
//   G4 (opt 3): median of the merged BAFphased is exactly 1 ⇒ mean used; seqlevels MT, GL000192.1, 2, Y,
//      1_gl000191_random, X, 10 → 2, 10, X, Y, MT, 1_gl000191_random, GL000192.1.

using static Seqeron.Genomics.Tests.Unit.Oncology.BattenbergCallSubclonesData;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_BattenbergMergeMask_Tests
{
    private static (List<OncologyAnalyzer.PhasedBafSegmentedSnp> Snps, List<OncologyAnalyzer.LogRProbe> LogR,
        IReadOnlyList<OncologyAnalyzer.BattenbergSegmentCall> Initial, OncologyAnalyzer.BattenbergSegmentMerge Merge) RunMerge(string name)
    {
        Genome g = Genomes[name];
        var (snps, logR) = Generate(g);
        var initial = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(
            OncologyAnalyzer.BuildBattenbergSegments(snps, logR), g.Rho, g.Psit, g.RSeed);
        var merge = OncologyAnalyzer.MergeBattenbergSegments(initial, snps, logR, g.Rho, g.Psit, bafOption: g.Option);
        return (snps, logR, initial, merge);
    }

    [TestCase("G1")]
    [TestCase("G2")]
    [TestCase("G3")]
    [TestCase("G4")]
    public void MergeBattenbergSegments_MatchesR(string name)
    {
        var (_, _, initial, merge) = RunMerge(name);
        Expectation e = Expected[name];
        Assert.Multiple(() =>
        {
            // The merge input is R's first determine_copynumber after set.seed(rseed).
            AssertRows(initial, e.Initial, $"{name} determine_copynumber #1");
            Assert.That(merge.Segments, Has.Count.EqualTo(e.Merged.Length), $"{name}: merged subclones rows.");
            for (int i = 0; i < Math.Min(merge.Segments.Count, e.Merged.Length); i++)
            {
                var a = merge.Segments[i];
                var x = e.Merged[i];
                string m = $"{name} merged row {i + 1}";
                Assert.That((a.Chromosome, a.Start, a.End), Is.EqualTo((x.Chr, x.Start, x.End)), m + " extent (GenomeInfoDb seqlevel order).");
                Assert.That(a.Baf, Is.EqualTo(x.Baf).Within(1e-12), m + " BAF.");
                Assert.That(a.LogR, Is.EqualTo(x.LogR).Within(1e-12), m + " LogR.");
                Assert.That((a.Call.Fit.PrimaryState.MajorCopyNumber, a.Call.Fit.PrimaryState.MinorCopyNumber), Is.EqualTo((x.Maj1, x.Min1)), m + " nMaj1_A/nMin1_A kept.");
                Assert.That(a.Call.Fit.PrimaryState.CellFraction, Is.EqualTo(x.Frac1).Within(1e-12), m + " frac1_A kept.");
            }

            var runs = Runs(merge.SegmentedSnps);
            Assert.That(runs, Has.Count.EqualTo(e.MergedRuns.Length), $"{name}: BAFsegmented runs.");
            for (int k = 0; k < Math.Min(runs.Count, e.MergedRuns.Length); k++)
            {
                var x = e.MergedRuns[k];
                Assert.That((runs[k].Chr, runs[k].First, runs[k].Last, runs[k].Count), Is.EqualTo((x.Chr, x.First, x.Last, x.Count)), $"{name} run {k + 1}.");
                Assert.That(runs[k].BafSeg, Is.EqualTo(x.BafSeg).Within(1e-12), $"{name} run {k + 1} BAFseg.");
            }
        });
    }

    [TestCase("G1")]
    [TestCase("G2")]
    [TestCase("G3")]
    [TestCase("G4")]
    public void MaskHighCopyNumberSegments_MatchesR(string name)
    {
        var (_, logR, _, merge) = RunMerge(name);
        Genome g = Genomes[name];
        Expectation e = Expected[name];
        // Copy-number states do not depend on the RNG stream (only the bootstrap columns do).
        var second = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(
            OncologyAnalyzer.BuildBattenbergSegments(merge.SegmentedSnps, logR), g.Rho, g.Psit, seed: 1);
        var mask = OncologyAnalyzer.MaskHighCopyNumberSegments(second, merge.SegmentedSnps, g.MaxState);
        Assert.Multiple(() =>
        {
            Assert.That((mask.MaskedCount, mask.MaskedSize), Is.EqualTo((e.MaskedCount, e.MaskedSize)), $"{name}: masked_count / masked_size.");
            Assert.That(mask.SegmentedSnps.Count(r => double.IsNaN(r.BafSegment)), Is.EqualTo(e.MaskedBafSegRows), $"{name}: masked BAFseg rows.");
            Assert.That(mask.Calls.Select(c => c.IsMasked), Is.EqualTo(e.Final.Select(r => double.IsNaN(r.Row[4]))), $"{name}: masked rows (nMaj1_A NA).");
            for (int i = 0; i < mask.Calls.Count; i++)
            {
                var c = mask.Calls[i];
                Assert.That(c.Fit.PrimaryState.CellFraction, Is.EqualTo(e.Final[i].Row[6]).Within(1e-12), $"{name} row {i + 1} frac1_A kept.");
                if (c.IsMasked && c.Solutions.Count > 0)
                {
                    Assert.That(c.Solutions[0].MajorCopyNumber1, Is.Null);
                }
            }
        });
    }

    [Test]
    public void MaskHighCopyNumberSegments_HalfOpenRange_KeepsFirstSnp()
    {
        var (_, logR, _, merge) = RunMerge("G1");
        var second = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(
            OncologyAnalyzer.BuildBattenbergSegments(merge.SegmentedSnps, logR), 0.7, 2.6, seed: 1);
        var mask = OncologyAnalyzer.MaskHighCopyNumberSegments(second, merge.SegmentedSnps, maxAllowedState: 5);
        var masked = mask.Calls.Single(c => c.IsMasked).Fit.Segment;
        var rows = mask.SegmentedSnps.Where(r => r.Chromosome == masked.Chromosome && r.Position >= masked.Start && r.Position <= masked.End).ToList();
        Assert.Multiple(() =>
        {
            Assert.That((masked.Chromosome, masked.Start, masked.End, masked.End - masked.Start), Is.EqualTo(("2", 6_120_000L, 6_900_000L, 780_000L)));
            Assert.That(double.IsNaN(rows[0].BafSegment), Is.False, "startpos < Position: the first SNP keeps BAFseg.");
            Assert.That(rows.Skip(1).All(r => double.IsNaN(r.BafSegment)), Is.True);
            Assert.That(OncologyAnalyzer.MaskHighCopyNumberSegments(second, merge.SegmentedSnps).MaskedCount, Is.Zero, "default 250: (6,1) kept.");
        });
    }

    // GenomeInfoDb::rankSeqlevels as used by makeGRangesFromDataFrame (R: s <- unique(x); s[rankSeqlevels(s)] <- s).
    [TestCase(new[] { "MT", "Y", "X", "10", "2", "GL000192.1", "hs37d5", "1_gl000191_random", "M" },
        new[] { "2", "10", "X", "Y", "M", "MT", "1_gl000191_random", "GL000192.1", "hs37d5" })]
    [TestCase(new[] { "chrX", "chr2", "chrUn_gl000220", "chr10", "chrM", "chr1_gl000191_random", "chr22" },
        new[] { "chr2", "chr10", "chr22", "chrX", "chrM", "chr1_gl000191_random", "chrUn_gl000220" })]
    [TestCase(new[] { "X", "IV", "II", "I", "M" }, new[] { "I", "II", "IV", "X", "M" })]
    [TestCase(new[] { "X", "Y", "1" }, new[] { "1", "X", "Y" })]
    [TestCase(new[] { "X", "I", "V" }, new[] { "I", "V", "X" })]
    [TestCase(new[] { "chr", "2b", "2a", "2A", "2", "3L", "3R", "Z", "W", "U", "Mito", "Xa", "Ya", "Uq", "MTx", "Wq", "Zq", "chrW", "b", "a" },
        new[] { "chrW", "chr", "2", "2A", "2a", "2b", "3L", "3R", "W", "Z", "U", "Wq", "Zq", "Xa", "Ya", "Uq", "Mito", "MTx", "a", "b" })]
    [TestCase(new[] { "ch1", "CH2", "CHR3", "chr4", "5" }, new[] { "CHR3", "chr4", "CH2", "ch1", "5" })]
    public void BattenbergSeqlevelOrder_MatchesGenomeInfoDb(string[] names, string[] expected)
    {
        Assert.That(OncologyAnalyzer.BattenbergSeqlevelOrder(names), Is.EqualTo(expected));
    }

    [Test]
    public void MergeBattenbergSegments_ConstantData_ThrowsLikeRTTest()
    {
        // Two adjacent segments in the same square with different solutions (clonal (2,1) vs sub-clonal) whose 12 logR and
        // 12 BAFphased values are all identical: R t.test stops with "data are essentially constant".
        var snps = new List<OncologyAnalyzer.PhasedBafSegmentedSnp>();
        var logR = new List<OncologyAnalyzer.LogRProbe>();
        for (int i = 0; i < 24; i++)
        {
            long pos = 1000 + i * 1000;
            snps.Add(new OncologyAnalyzer.PhasedBafSegmentedSnp("1", pos, 0.6, 0.6, i < 12 ? 0.6296 : 0.668));
            logR.Add(new OncologyAnalyzer.LogRProbe("1", pos, 0.1));
        }

        var calls = OncologyAnalyzer.FitSubclonalCopyNumberWithBootstrap(OncologyAnalyzer.BuildBattenbergSegments(snps, logR), 0.7, 2.6, 1);
        Assert.That(calls.Select(c => c.Fit.IsSubclonal), Is.EqualTo(new[] { false, true }));
        Assert.Throws<InvalidOperationException>(() => OncologyAnalyzer.MergeBattenbergSegments(calls, snps, logR, 0.7, 2.6));
    }

    [Test]
    public void MergeBattenbergSegments_InvalidInput_Throws()
    {
        var (snps, logR, initial, _) = RunMerge("G1");
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.MergeBattenbergSegments(null!, snps, logR, 0.7, 2.6));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.MergeBattenbergSegments(initial, null!, logR, 0.7, 2.6));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.MergeBattenbergSegments(initial, snps, null!, 0.7, 2.6));
            var inf = logR.ToList();
            inf[3] = inf[3] with { LogR = double.PositiveInfinity };
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.MergeBattenbergSegments(initial, snps, inf, 0.7, 2.6), "±∞ logR");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.MergeBattenbergSegments(initial, snps, logR.Where(p => p.Chromosome != "2").ToList(), 0.7, 2.6),
                "chromosome without logR probes (R stopifnot)");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.MergeBattenbergSegments(initial.Where(c => c.Fit.Segment.Chromosome != "2").ToList(), snps, logR, 0.7, 2.6),
                "chromosome without calls (R stopifnot)");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.MergeBattenbergSegments(initial, snps, logR, 0.7, 2.6, bafOption: (OncologyAnalyzer.BattenbergSegmentBafOption)4));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.MergeBattenbergSegments(initial, snps, logR, 0.0, 2.6));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.MaskHighCopyNumberSegments(null!, snps));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.MaskHighCopyNumberSegments(initial, null!));
        });
    }

    private static List<(string Chr, long First, long Last, int Count, double BafSeg)> Runs(IReadOnlyList<OncologyAnalyzer.PhasedBafSegmentedSnp> rows)
    {
        var runs = new List<(string, long, long, int, double)>();
        int from = 0;
        for (int i = 0; i < rows.Count; i++)
        {
            if (i < rows.Count - 1 && rows[i + 1].Chromosome == rows[i].Chromosome && rows[i + 1].BafSegment == rows[i].BafSegment)
            {
                continue;
            }

            runs.Add((rows[i].Chromosome, rows[from].Position, rows[i].Position, i - from + 1, rows[i].BafSegment));
            from = i + 1;
        }

        return runs;
    }
}
