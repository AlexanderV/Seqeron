// ONCO-ASCAT-001 — Battenberg phased-BAF segmentation (segment.baf.phased) feeding determine_copynumber (B24 F40)
// Evidence: docs/Evidence/ONCO-ASCAT-001-Evidence.md (§ F40)
// TestSpec: tests/TestSpecs/ONCO-ASCAT-001.md (§14)
// Source: Wedge-lab/battenberg master (57a8f7e), R/segmentation.R (segment.baf.phased, adjustSegmValues),
//         R/fastPCF.R (selectFastPcf, runFastPcf, runPcfSubset, filterMarkS4, PottsCompact, findEst, getMad),
//         R/fitcopynumber.R (determine_copynumber), R/orderEdges.R. Nik-Zainal S et al. (2012) Cell 149:994.
//
// Inputs are generated in-test by a deterministic integer LCG that is replicated verbatim in the R harness
// (x ← (1664525·x + 1013904223) mod 2^24, u = x / 2^24; all operations exact in IEEE double in both languages):
// a piecewise-constant BAF level per true segment, haplotype blocks of `block` SNPs randomly switched (BAF ↦ 1 − BAF),
// and noise sd·(u₁ + u₂ + u₃ − 1.5)·2, clamped to [0, 1]. Expected rows: Battenberg R functions sourced verbatim in
// R 4.3.3 (plots stubbed), run-length summary of the output (chromosome, first/last position, SNPs, BAFseg).

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_BattenbergPhasedSegmentation_Tests
{
    private static List<OncologyAnalyzer.PhasedBafSnp> Gen(
        string chr, int n, double[] levels, int[] bounds, int block, double seed, double sdv, long pos0, long step,
        int gapAfter = -1, long gap = 0)
    {
        double x = seed;
        double Next()
        {
            x = (1664525 * x + 1013904223) % 16777216;
            return x / 16777216;
        }

        var result = new List<OncologyAnalyzer.PhasedBafSnp>(n);
        bool flip = false;
        int seg = 1;
        for (int i = 1; i <= n; i++)
        {
            while (seg < bounds.Length && i > bounds[seg - 1])
            {
                seg++;
            }

            if ((i - 1) % block == 0)
            {
                flip = Next() < 0.5;
            }

            double noise = sdv * ((Next() + Next() + Next()) - 1.5) * 2;
            double b = levels[seg - 1] + noise;
            if (flip)
            {
                b = 1 - b;
            }

            long pos = pos0 + i * step + (gapAfter > 0 && i > gapAfter ? gap : 0);
            result.Add(new OncologyAnalyzer.PhasedBafSnp(chr, pos, Math.Min(1, Math.Max(0, b))));
        }

        return result;
    }

    internal static List<OncologyAnalyzer.PhasedBafSnp> Track(string name) => name switch
    {
        "t1" => Gen("1", 300, new[] { 0.5, 0.72, 0.6 }, new[] { 100, 200, 300 }, 15, 12345, 0.04, 1000, 1000),
        "t2" => Gen("2", 2000, new[] { 0.55, 0.8, 0.5, 0.67 }, new[] { 500, 900, 1500, 2000 }, 40, 777, 0.05, 50000, 1500, 1200, 4_000_000)
            .Concat(Gen("3", 40, new[] { 0.7 }, new[] { 40 }, 10, 99, 0.03, 100, 100)).ToList(),
        "t3" => Gen("4", 16000, new[] { 0.5, 0.75, 0.62, 0.9, 0.5 }, new[] { 3000, 7000, 9000, 12500, 16000 }, 60, 4242, 0.05, 10000, 200),
        "t4" => Gen("5", 400, new[] { 0.6, 0.8 }, new[] { 250, 400 }, 20, 31337, 0.04, 1000, 1000),
        "t5" => Gen("6", 1500, new[] { 0.55, 0.62, 0.58, 0.7 }, new[] { 300, 700, 1100, 1500 }, 7, 2024, 0.1, 1, 300),
        _ => throw new ArgumentOutOfRangeException(nameof(name)),
    };

    private static List<(string Chr, long First, long Last, int Count, double BafSeg)> Runs(
        IReadOnlyList<OncologyAnalyzer.PhasedBafSegmentedSnp> rows)
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

    private static void AssertRuns(
        IReadOnlyList<OncologyAnalyzer.PhasedBafSegmentedSnp> rows,
        (string Chr, long First, long Last, int Count, double BafSeg)[] expected,
        double sumPhased, double sumSeg, int n)
    {
        var runs = Runs(rows);
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.EqualTo(n), "Battenberg output rows.");
            Assert.That(runs, Has.Count.EqualTo(expected.Length), "Battenberg segments (BAFseg runs).");
            for (int k = 0; k < Math.Min(runs.Count, expected.Length); k++)
            {
                Assert.That((runs[k].Chr, runs[k].First, runs[k].Last, runs[k].Count),
                    Is.EqualTo((expected[k].Chr, expected[k].First, expected[k].Last, expected[k].Count)), $"segment {k + 1} extent.");
                Assert.That(runs[k].BafSeg, Is.EqualTo(expected[k].BafSeg).Within(1e-12), $"segment {k + 1} BAFseg.");
            }

            Assert.That(rows.Sum(r => r.BafPhased), Is.EqualTo(sumPhased).Within(1e-12 * sumPhased), "Σ BAFphased (haplotype-block correction).");
            Assert.That(rows.Sum(r => r.BafSegment), Is.EqualTo(sumSeg).Within(1e-12 * sumSeg), "Σ BAFseg (R sums in long double; relative 1e−12).");
        });
    }

    // t1: 300 SNPs (< 1000 ⇒ runFastPcf fractions 0.15/0.15), 15-SNP switched haplotype blocks, default options.
    [Test]
    public void SegmentPhasedBaf_Track1_MatchesBattenberg()
    {
        var rows = OncologyAnalyzer.SegmentPhasedBaf(Track("t1"));
        AssertRuns(rows, new[]
        {
            ("1", 2000L, 101000L, 100, 0.49788217067718499),
            ("1", 102000L, 201000L, 100, 0.71850197553634598),
            ("1", 202000L, 301000L, 100, 0.58956708192825302),
        }, 176.48445771217348, 180.59512281417841, 300);
    }

    // t2: chr2 2000 SNPs (fractions 0.12/0.05) with a 4 Mb gap (≥ 3 Mb ⇒ presegment break) + chr3 40 SNPs (< 50 ⇒ mean);
    // calc_seg_baf_option 3 (default) and 1 agree here (no median of 0/1), option 2 keeps the PCF means.
    [TestCase(OncologyAnalyzer.BattenbergSegmentBafOption.MedianUnlessExtreme, TestName = "Track2_Option3")]
    [TestCase(OncologyAnalyzer.BattenbergSegmentBafOption.Median, TestName = "Track2_Option1")]
    public void SegmentPhasedBaf_Track2_MatchesBattenberg(OncologyAnalyzer.BattenbergSegmentBafOption option)
    {
        var rows = OncologyAnalyzer.SegmentPhasedBaf(
            Track("t2"), new OncologyAnalyzer.BattenbergPhasedSegmentationOptions { SegmentBafOption = option });
        AssertRuns(rows, new[]
        {
            ("2", 51500L, 800000L, 500, 0.54674800336360896),
            ("2", 801500L, 1400000L, 400, 0.80187466144561803),
            ("2", 1401500L, 1850000L, 300, 0.50032104849815395),
            ("2", 5851500L, 6300000L, 300, 0.49780743718147302),
            ("2", 6301500L, 7050000L, 500, 0.66956448376178701),
            ("3", 200L, 4100L, 40, 0.50574477136135099),
        }, 1246.572572066784, 1248.5744446992874, 2040);
    }

    [Test]
    public void SegmentPhasedBaf_Track2_MeanOption_MatchesBattenberg()
    {
        var rows = OncologyAnalyzer.SegmentPhasedBaf(
            Track("t2"),
            new OncologyAnalyzer.BattenbergPhasedSegmentationOptions { SegmentBafOption = OncologyAnalyzer.BattenbergSegmentBafOption.Mean });
        AssertRuns(rows, new[]
        {
            ("2", 51500L, 800000L, 500, 0.54146806168556205),
            ("2", 801500L, 1400000L, 400, 0.80023094666004202),
            ("2", 1401500L, 1850000L, 300, 0.50063822877407005),
            ("2", 5851500L, 6300000L, 300, 0.50009172968069704),
            ("2", 6301500L, 7050000L, 500, 0.67059457948207901),
            ("3", 200L, 4100L, 40, 0.505747132062912),
        }, 1246.572572066784, 1246.572572066784, 2040);
    }

    // t3: 16 000 SNPs (≥ 15 000 ⇒ runPcfSubset, 5000-SNP windows advancing by 4000).
    [Test]
    public void SegmentPhasedBaf_Track3_PcfSubset_MatchesBattenberg()
    {
        var rows = OncologyAnalyzer.SegmentPhasedBaf(Track("t3"));
        AssertRuns(rows, new[]
        {
            ("4", 10200L, 610000L, 3000, 0.50033228099346105),
            ("4", 610200L, 1410000L, 4000, 0.74888144135475199),
            ("4", 1410200L, 1810000L, 2000, 0.61671185612678503),
            ("4", 1810200L, 2510000L, 3500, 0.89983910620212604),
            ("4", 2510200L, 3210000L, 3500, 0.50094892680644998),
        }, 10639.269587749242, 10632.704436182976, 16000);
    }

    // t4: prior breakpoints (SVs at 150 500 and 320 000) with gamma 5, kmin 5; and no_segmentation on the same track.
    [Test]
    public void SegmentPhasedBaf_Track4_PriorBreakpoints_MatchesBattenberg()
    {
        var rows = OncologyAnalyzer.SegmentPhasedBaf(Track("t4"), new OncologyAnalyzer.BattenbergPhasedSegmentationOptions
        {
            Gamma = 5,
            Kmin = 5,
            PriorBreakpoints = new[]
            {
                new OncologyAnalyzer.BattenbergPriorBreakpoint("5", 150500),
                new OncologyAnalyzer.BattenbergPriorBreakpoint("5", 320000),
            },
        });
        AssertRuns(rows, new[]
        {
            ("5", 2000L, 150000L, 149, 0.59599266052246103),
            ("5", 151000L, 251000L, 101, 0.60799378395080605),
            ("5", 252000L, 320000L, 69, 0.80668817043304497),
            ("5", 321000L, 401000L, 81, 0.80411475181579595),
        }, 270.88505104064944, 271.00505725383766, 400);
    }

    [Test]
    public void SegmentPhasedBaf_Track4_NoSegmentation_MatchesBattenberg()
    {
        var rows = OncologyAnalyzer.SegmentPhasedBaf(Track("t4"), new OncologyAnalyzer.BattenbergPhasedSegmentationOptions { NoSegmentation = true });
        AssertRuns(rows, new[] { ("5", 2000L, 401000L, 400, 0.63947307586669999) }, 270.88505104064944, 255.78923034668, 400);
    }

    // t5: noisy (sd 0.1), 7-SNP haplotype blocks, close levels — breakpoints placed by the penalised cost, not the truth.
    [Test]
    public void SegmentPhasedBaf_Track5_Noisy_MatchesBattenberg()
    {
        var rows = OncologyAnalyzer.SegmentPhasedBaf(Track("t5"));
        AssertRuns(rows, new[]
        {
            ("6", 301L, 100801L, 336, 0.50309495329856901),
            ("6", 101101L, 210901L, 367, 0.61098066091537495),
            ("6", 211201L, 329701L, 396, 0.50803442597389203),
            ("6", 330001L, 450001L, 401, 0.70338219404220603),
        }, 869.7522670245171, 876.50769936084771, 1500);
    }

    // Battenberg-like logR track: per chromosome-local SNP index i, lv[seg] + ((i mod 7) − 3)·0.01 at every SNP position.
    internal static List<OncologyAnalyzer.LogRProbe> LogRTrack(
        IReadOnlyList<OncologyAnalyzer.PhasedBafSegmentedSnp> rows, double[] lv, int[] bounds)
    {
        var probes = new List<OncologyAnalyzer.LogRProbe>();
        foreach (var chr in rows.Select(r => r.Chromosome).Distinct())
        {
            int seg = 1, i = 0;
            foreach (var r in rows.Where(r => r.Chromosome == chr))
            {
                i++;
                while (seg < bounds.Length && i > bounds[seg - 1])
                {
                    seg++;
                }

                probes.Add(new OncologyAnalyzer.LogRProbe(chr, r.Position, lv[seg - 1] + ((i % 7) - 3) * 0.01));
            }
        }

        return probes;
    }

    // End to end: segment.baf.phased → determine_copynumber (R: set.seed, determine_copynumber(BAFvals, LogRvals, ρ,
    // ψ_all, γ = 1, ctrans, ctrans, maxdist 0.01, siglevel 0.05, noperms 1000, cn_upper_limit 1000)). Rows:
    // (chr, startpos, endpos, LogR, pval, nMaj1_A, nMin1_A, frac1_A, nMaj2_A, nMin2_A, frac2_A); −1 = NA (clonal).
    private static readonly (string Chr, long Start, long End, double LogR, double PVal, int Maj1, int Min1, double Frac1, int Maj2, int Min2, double Frac2)[] E1 =
    {
        ("1", 2000, 101000, -0.00029999999999999997, 1.0, 1, 1, 1.0, -1, -1, double.NaN),
        ("1", 102000, 201000, 0.15010000000000001119, 6.2072657065460898e-36, 2, 0, 0.47709449087559203, 2, 1, 0.52290550912440803),
        ("1", 202000, 301000, -0.19950000000000001066, 5.0025818343502595e-13, 1, 1, 0.37649750519839498, 2, 1, 0.62350249480160502),
    };

    private static readonly (string Chr, long Start, long End, double LogR, double PVal, int Maj1, int Min1, double Frac1, int Maj2, int Min2, double Frac2)[] E2 =
    {
        ("2", 51500, 800000, -1.0006000000000000e-01, 9.5844256487078194e-46, 1, 1, 0.75731976285090297, 2, 1, 0.24268023714909701),
        ("2", 801500, 1400000, 3.0002499999999999e-01, 2.2067308446755200e-101, 3, 0, 0.39163557478371203, 3, 1, 0.60836442521628797),
        ("2", 1401500, 1850000, -3.3333333333333301e-05, 1.0, 1, 1, 1.0, -1, -1, double.NaN),
        ("2", 5851500, 6300000, 5.0000000000000003e-02, 1.0, 2, 2, 1.0, -1, -1, double.NaN),
        ("2", 6301500, 7050000, 2.0005999999999999e-01, 3.2600528882290903e-20, 2, 1, 0.79257756414276204, 3, 1, 0.20742243585723799),
        ("3", 200, 4100, -1.0000000000000001e-01, 1.0, 1, 1, 1.0, -1, -1, double.NaN),
    };

    [TestCase("e1")]
    [TestCase("e2")]
    public void BuildBattenbergSegments_FeedsSnpTest_MatchesDetermineCopynumber(string name)
    {
        IReadOnlyList<OncologyAnalyzer.PhasedBafSegmentedSnp> rows;
        List<OncologyAnalyzer.LogRProbe> logR;
        double rho, psit;
        var expected = name == "e1" ? E1 : E2;
        if (name == "e1")
        {
            rows = OncologyAnalyzer.SegmentPhasedBaf(Track("t1"));
            logR = LogRTrack(rows, new[] { 0.0, 0.15, -0.2 }, new[] { 100, 200, 300 });
            logR.Add(new OncologyAnalyzer.LogRProbe("1", 150500, double.PositiveInfinity)); // excluded (!is.infinite)
            (rho, psit) = (0.7, 2.6);
        }
        else
        {
            rows = OncologyAnalyzer.SegmentPhasedBaf(Track("t2"));
            logR = LogRTrack(rows, new[] { -0.1, 0.3, 0.0, 0.05, 0.2 }, new[] { 500, 900, 1200, 1500, 2000 });
            (rho, psit) = (0.85, 3.0);
        }

        var segments = OncologyAnalyzer.BuildBattenbergSegments(rows, logR);
        var fits = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(segments, rho, psit);

        Assert.That(fits, Has.Count.EqualTo(expected.Length), "determine_copynumber segments (switchpoints).");
        Assert.Multiple(() =>
        {
            for (int k = 0; k < expected.Length; k++)
            {
                var e = expected[k];
                var s = segments[k].Segment;
                var f = fits[k];
                Assert.That((s.Chromosome, s.Start, s.End), Is.EqualTo((e.Chr, e.Start, e.End)), $"seg {k + 1} chr/startpos/endpos.");
                Assert.That(s.MeanLogR, Is.EqualTo(e.LogR).Within(1e-12), $"seg {k + 1} LogR.");
                Assert.That(f.PValue, Is.EqualTo(e.PVal).Within(1e-9 * e.PVal), $"seg {k + 1} pval.");
                Assert.That((f.Fit.PrimaryState.MajorCopyNumber, f.Fit.PrimaryState.MinorCopyNumber), Is.EqualTo((e.Maj1, e.Min1)), $"seg {k + 1} nMaj1_A:nMin1_A.");
                Assert.That(f.Fit.PrimaryState.CellFraction, Is.EqualTo(e.Frac1).Within(1e-12), $"seg {k + 1} frac1_A.");
                if (e.Maj2 >= 0)
                {
                    Assert.That((f.Fit.SecondaryState!.Value.MajorCopyNumber, f.Fit.SecondaryState!.Value.MinorCopyNumber), Is.EqualTo((e.Maj2, e.Min2)), $"seg {k + 1} nMaj2_A:nMin2_A.");
                    Assert.That(f.Fit.SecondaryState!.Value.CellFraction, Is.EqualTo(e.Frac2).Within(1e-12), $"seg {k + 1} frac2_A.");
                }
                else
                {
                    Assert.That(f.Fit.SecondaryState, Is.Null, $"seg {k + 1} clonal.");
                }
            }
        });
    }

    // ---- Invariants / edge cases (Battenberg semantics) ----

    [Test]
    public void SegmentPhasedBaf_MissingBafDropped_AndSmallPresegmentUsesMean()
    {
        var snps = new List<OncologyAnalyzer.PhasedBafSnp>
        {
            new("X", 10, 0.2), new("X", 20, double.NaN), new("X", 30, 0.3), new("X", 40, 0.9),
        };
        var rows = OncologyAnalyzer.SegmentPhasedBaf(snps);
        // < 50 SNPs: BAFsegm = mean(BAF) = 0.4666… ≤ 0.5 ⇒ BAFphased = 1 − BAF = (0.8, 0.7, 0.1);
        // BAFphseg = mean = 0.5333…, option 3 ⇒ median 0.7 (not 0/1).
        Assert.Multiple(() =>
        {
            Assert.That(rows.Select(r => r.Position), Is.EqualTo(new long[] { 10, 30, 40 }));
            Assert.That(rows.Select(r => r.BafPhased), Is.EqualTo(new[] { 0.8, 0.7, 1 - 0.9 }).Within(1e-15));
            Assert.That(rows.Select(r => r.BafSegment), Is.All.EqualTo(0.7).Within(1e-15));
        });
    }

    [Test]
    public void SegmentPhasedBaf_ExtremeMedian_KeepsMeanUnderOption3()
    {
        var snps = new[] { new OncologyAnalyzer.PhasedBafSnp("1", 1, 1.0), new("1", 2, 1.0), new("1", 3, 0.7) };
        var opt3 = OncologyAnalyzer.SegmentPhasedBaf(snps);
        var opt1 = OncologyAnalyzer.SegmentPhasedBaf(snps, new OncologyAnalyzer.BattenbergPhasedSegmentationOptions { SegmentBafOption = OncologyAnalyzer.BattenbergSegmentBafOption.Median });
        Assert.Multiple(() =>
        {
            Assert.That(opt3.Select(r => r.BafSegment), Is.All.EqualTo(2.7 / 3).Within(1e-15), "median 1 ⇒ mean kept.");
            Assert.That(opt1.Select(r => r.BafSegment), Is.All.EqualTo(1.0), "option 1 ⇒ median.");
        });
    }

    [Test]
    public void SegmentPhasedBaf_ArgumentGuards()
    {
        var ok = new[] { new OncologyAnalyzer.PhasedBafSnp("1", 1, 0.5) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.SegmentPhasedBaf(null!));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.SegmentPhasedBaf(new[] { new OncologyAnalyzer.PhasedBafSnp("1", 1, 1.2) }));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.SegmentPhasedBaf(ok, new OncologyAnalyzer.BattenbergPhasedSegmentationOptions { Kmin = 0 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.SegmentPhasedBaf(ok, new OncologyAnalyzer.BattenbergPhasedSegmentationOptions { PhaseKmin = 15 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.SegmentPhasedBaf(ok, new OncologyAnalyzer.BattenbergPhasedSegmentationOptions { Gamma = double.NaN }));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.SegmentPhasedBaf(ok, new OncologyAnalyzer.BattenbergPhasedSegmentationOptions { SegmentBafOption = (OncologyAnalyzer.BattenbergSegmentBafOption)4 }));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.BuildBattenbergSegments(null!, Array.Empty<OncologyAnalyzer.LogRProbe>()));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.BuildBattenbergSegments(
                new[] { new OncologyAnalyzer.PhasedBafSegmentedSnp("1", 1, 0.5, 0.5, 1.5) }, Array.Empty<OncologyAnalyzer.LogRProbe>()));
        });
    }

    [Test]
    public void BuildBattenbergSegments_NoLogRProbes_LogRIsZero()
    {
        var rows = new[]
        {
            new OncologyAnalyzer.PhasedBafSegmentedSnp("1", 100, 0.6, 0.6, 0.65),
            new OncologyAnalyzer.PhasedBafSegmentedSnp("1", 200, 0.3, 0.7, 0.65),
            new OncologyAnalyzer.PhasedBafSegmentedSnp("2", 300, 0.7, 0.7, 0.65), // chromosome change ⇒ new segment
        };
        var segs = OncologyAnalyzer.BuildBattenbergSegments(rows, new[] { new OncologyAnalyzer.LogRProbe("1", 150, double.NaN) });
        Assert.Multiple(() =>
        {
            Assert.That(segs, Has.Count.EqualTo(2));
            Assert.That(segs[0].Segment.MeanLogR, Is.EqualTo(0.0), "is.na(LogR) ⇒ 0.");
            Assert.That(segs[0].PhasedSnpBafs, Is.EqualTo(new[] { 0.6, 0.7 }), "BAFke = BAFphased.");
            Assert.That((segs[0].Segment.Start, segs[0].Segment.End, segs[0].Segment.LocusCount), Is.EqualTo((100L, 200L, 2)));
        });
    }
}
