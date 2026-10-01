using Seqeron.Genomics.Tests.Properties;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier fuzz tests for the API added by the B04 completeness audit (docs/Validation/review-2026-09/B04.md,
/// F34–F65): random parameter sets and adversarial inputs (empty, all-N, homopolymers, IUPAC, lower case, U, long
/// periods) for the DUST engines, longdust, BBDuk / windowed / k-mer entropy, the TRF parameter API with random
/// apparent-size tables and its formatters, degenerate repeats in every reporting mode and the MISA definition parser /
/// scan. Every call must either succeed with well-formed output or throw the documented argument exception.
///
/// Test Units: SEQ-COMPLEX-DUST-001, SEQ-COMPLEX-001, SEQ-COMPLEX-WINDOW-001, SEQ-COMPLEX-KMER-001, REP-APPROX-001,
/// REP-DIRECT-001, REP-STR-001, REP-TANDEM-001.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Complexity")]
public class B04AuditFuzzTests
{
    private const int Seed = 1004;

    private static string Random(Random rng, string alphabet, int length) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    /// <summary>Adversarial nucleotide inputs (IUPAC DNA unless <paramref name="anySymbols"/>).</summary>
    private static string Adversarial(Random rng, bool anySymbols = false)
    {
        int n = rng.Next(0, 400);
        switch (rng.Next(9))
        {
            case 0: return string.Empty;
            case 1: return new string('N', n);
            case 2: return new string("ACGT"[rng.Next(4)], n);
            case 3: return string.Concat(Enumerable.Repeat(Random(rng, "ACGT", rng.Next(1, 7)), n / 3 + 1));
            case 4: return Random(rng, "ACGTRYKMSWBDHVN", n);
            case 5: return B04AuditProperties.Segmented(rng, n, "N").ToLowerInvariant();
            case 6: return anySymbols ? Random(rng, "ACGTUXacgtu-*N", n) : B04AuditProperties.Segmented(rng, n, "NNNN");
            case 7: return Random(rng, "ACGT", n);
            default: return B04AuditProperties.Segmented(rng, n, "RYN");
        }
    }

    private static bool WellFormed(IReadOnlyList<(int Start, int End)> iv, int n)
    {
        for (int i = 0; i < iv.Count; i++)
        {
            if (iv[i].Start < 0 || iv[i].End > n || iv[i].Start >= iv[i].End) return false;
            if (i > 0 && iv[i - 1].End > iv[i].Start) return false;
        }
        return true;
    }

    /// <summary>
    /// F35/F36/F53: both DUST engines on adversarial IUPAC input with random valid parameters return well-formed intervals
    /// and a mask of the input length; dustmasker rejects non-IUPAC symbols with ArgumentException and parameters outside
    /// its ranges (window 8–64, 10·T an integer 2–64) with ArgumentOutOfRangeException; linker 0 / 33 is rejected.
    /// </summary>
    [Test]
    public void DustEngines_AdversarialInput_WellFormedOrDocumentedException()
    {
        var rng = new Random(Seed);
        for (int trial = 0; trial < 400; trial++)
        {
            string s = Adversarial(rng);
            int w = rng.Next(8, 65), level = rng.Next(2, 65), linker = rng.Next(1, 33);
            foreach (var engine in new[] { DustEngine.Sdust, DustEngine.Dustmasker })
            {
                var iv = SequenceComplexity.FindLowComplexityIntervals(s, w, level / 10.0, linker, engine);
                Assert.That(WellFormed(iv, s.Length), Is.True, $"{engine} {s}");
                Assert.That(SequenceComplexity.MaskLowComplexity(s, w, level / 10.0, 'N', linker, trial % 2 == 0, engine), Has.Length.EqualTo(s.Length));
            }
        }

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGUACGT", 20, 2.0, 1, DustEngine.Dustmasker));
            Assert.Throws<ArgumentException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT*ACGT", 20, 2.0, 1, DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 7, 2.0, 1, DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 65, 2.0, 1, DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 64, 2.05, 1, DustEngine.Dustmasker));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 64, 2.0, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityIntervals("ACGT", 64, 2.0, 33));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateDustScore("ACGTACGT", 4));
        });
    }

    /// <summary>
    /// F34: longdust with random k (1–11), window, threshold, X-drop, start count, strand, approximate and GC options on
    /// adversarial input returns well-formed intervals and finite scores; invalid options throw. (k ≤ 11 keeps the fixture
    /// fast: like longdust.c, which memsets its 4^k-entry count table before every backward / forward scan
    /// (longdust.c:210, :228) and enumerates all 4^k k-mers for the GC classes (:39), the port costs O(4^k) per scan, so
    /// k = 12–14 takes seconds per call — parity with the reference, not a defect.)
    /// </summary>
    [Test]
    public void Longdust_RandomOptions_AdversarialInput_WellFormed()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 250; trial++)
        {
            string s = Adversarial(rng, anySymbols: true);
            int k = rng.Next(1, 12), w = rng.Next(1, 600), xdrop = rng.Next(0, 80), minStart = rng.Next(2, 6);
            double t = rng.Next(1, 30) / 10.0;
            double? gc = rng.Next(3) == 0 ? rng.Next(5, 96) / 100.0 : null;
            var iv = SequenceComplexity.FindLongdustRegions(s, k, w, t, xdrop, minStart, rng.Next(2) == 0, rng.Next(2) == 0, gc);
            Assert.That(WellFormed(iv, s.Length), Is.True, $"k={k} w={w} {s}");
            Assert.That(double.IsFinite(SequenceComplexity.CalculateLongdustScore(s, k, gc)), Is.True, s);
        }

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", k: 15));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", windowSize: 65535));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", threshold: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLongdustRegions("ACGT", minStartCount: 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateLongdustScore("ACGT", 7, 1.0));
        });
    }

    /// <summary>
    /// F37–F39, F54, F55: fixed-alphabet LC, windowed complexity, LCR, BBDuk and all k-mer entropy corrections on
    /// adversarial input (U, lower case, IUPAC, gaps) give finite values in their ranges; regions lie in the sequence.
    /// </summary>
    [Test]
    public void EntropyAndLcFamilies_AdversarialInput_FiniteAndInRange()
    {
        var rng = new Random(Seed + 2);
        var corrections = new[] { KmerEntropyCorrection.None, KmerEntropyCorrection.MillerMadow, KmerEntropyCorrection.Grassberger };
        for (int trial = 0; trial < 300; trial++)
        {
            string s = Adversarial(rng, anySymbols: true);
            int m = rng.Next(1, 20), w = rng.Next(1, 80), step = rng.Next(1, 20), k = rng.Next(1, 8);
            int distinct = s.ToUpperInvariant().Distinct().Count();
            double lc = SequenceComplexity.CalculateLinguisticComplexity(s, m, Math.Max(1, distinct) + rng.Next(0, 3));
            Assert.That(lc, Is.InRange(0.0, 1.0), s);
            foreach (var p in SequenceComplexity.CalculateWindowedComplexity(s, w, step, m))
            {
                Assert.That(p.LinguisticComplexity, Is.InRange(0.0, 1.0), s);
                Assert.That(p.ShannonEntropy, Is.InRange(0.0, 2.0 + 1e-12), s);
            }
            foreach (var r in SequenceComplexity.FindLowComplexityRegions(s, w, rng.Next(0, 21) / 10.0))
                Assert.That(r.Start >= 0 && r.End < s.Length && r.Length == r.End - r.Start + 1, Is.True, s);
            foreach (var r in SequenceComplexity.FindLowEntropyRegionsBbduk(s, rng.Next(0, 101) / 100.0, Math.Max(Math.Min(k, 5) + 1, w), Math.Min(k, 5)))
                Assert.That(r.Start >= 0 && r.End < s.Length && r.Length == r.End - r.Start + 1, Is.True, s);
            foreach (var c in corrections)
            {
                double v = SequenceComplexity.CalculateKmerEntropy(s, k, c);
                double norm = SequenceComplexity.CalculateKmerEntropy(s, k, c, normalize: true);
                Assert.That(double.IsFinite(v) && double.IsFinite(norm), Is.True, $"{c} {s}");
            }
        }

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => SequenceComplexity.CalculateLinguisticComplexity("ACGT", 3, 3));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateLinguisticComplexity("ACGT", 3, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowEntropyRegionsBbduk("ACGT", 1.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.FindLowComplexityRegions("ACGT", 4, double.NaN).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceComplexity.CalculateKmerEntropy("ACGT", 2, (KmerEntropyCorrection)9));
            Assert.Throws<ArgumentException>(() => SequenceComplexity.ParseKmerEntropyCorrection("chao-shen"));
        });
    }

    /// <summary>
    /// F41–F46, F56: TRF with random parameter sets, random valid apparent-size tables and adversarial input (all-N,
    /// homopolymers, IUPAC, lower case, a ≥ 1365-bp pattern at PI 20) never throws and keeps every result invariant;
    /// invalid parameter sets / tables are rejected by Validate.
    /// </summary>
    [Test]
    public void Trf_RandomParametersAndTables_AdversarialInput_KeepInvariants()
    {
        var rng = new Random(Seed + 3);
        int len = TandemRepeatsFinderParameters.ApparentSizeTableLength;
        for (int trial = 0; trial < 160; trial++)
        {
            string s = trial % 4 == 0 ? B04AuditProperties.TandemRich(rng, rng.Next(1, 600), withN: true) : Adversarial(rng);
            var table = new int[len];
            for (int d = 1; d < len; d++) table[d] = rng.Next(0, Math.Max(d, 20));
            var p = new TandemRepeatsFinderParameters
            {
                MatchWeight = rng.Next(1, 4), MismatchPenalty = rng.Next(1, 8), IndelPenalty = rng.Next(1, 8),
                MatchProbability = rng.Next(2) == 0 ? 75 : 80, IndelProbability = rng.Next(1, 101),
                MinScore = rng.Next(1, 100), MaxPeriod = rng.Next(1, 2001), MaxRepeatLength = rng.Next(1, 5000),
                EliminateRedundancy = rng.Next(3) != 0, FlankLength = rng.Next(0, 3) * 25,
                ApparentSizeTable = rng.Next(2) == 0 ? table : null,
            };
            var rows = RepeatFinder.FindApproximateTandemRepeats(s, p, out int count);
            foreach (var r in rows)
                Assert.That(B04AuditProperties.CheckTrfResult(s, p, r, count), Is.Empty, $"{p} {s}");
            Assert.That(RepeatFinder.MaskApproximateTandemRepeats(s, rows, softMask: true), Has.Length.EqualTo(s.Length));
        }

        // A long-period pattern (> 1365 bp) at PI 20 exercises the narrow band (TRF aborts above 150 band cells).
        string unit = Random(rng, "ACGT", 1400);
        string longPeriod = Random(rng, "ACGT", 50) + unit + unit + Random(rng, "ACGT", 50);
        var lp = TandemRepeatsFinderParameters.Recommended with { IndelProbability = 20, MaxPeriod = 2000 };
        var longRows = RepeatFinder.FindApproximateTandemRepeats(longPeriod, lp, out int longCount);
        foreach (var r in longRows)
            Assert.That(B04AuditProperties.CheckTrfResult(longPeriod, lp, r, longCount), Is.Empty);

        Assert.Multiple(() =>
        {
            var bad = new int[len];
            bad[100] = 100; // > max(d, 20) − 1
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RepeatFinder.FindApproximateTandemRepeats("ACGT", TandemRepeatsFinderParameters.Recommended with { ApparentSizeTable = bad }).ToList());
            Assert.Throws<ArgumentException>(() =>
                RepeatFinder.FindApproximateTandemRepeats("ACGT", TandemRepeatsFinderParameters.Recommended with { ApparentSizeTable = new int[5] }).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RepeatFinder.FindApproximateTandemRepeats("ACGT", TandemRepeatsFinderParameters.Recommended with { MatchProbability = 90 }).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RepeatFinder.FindApproximateTandemRepeats("ACGT", TandemRepeatsFinderParameters.Recommended with { MaxPeriod = 2001 }).ToList());
            Assert.Throws<ArgumentOutOfRangeException>(() => TandemRepeatsFinderParameters.ExactApparentSizeTable(70));
        });
    }

    /// <summary>
    /// F57/F64/F65: the .dat / -ngs / HTML table / summary / alignment-page formatters accept every result list the
    /// finder produces (including none, with and without an explicit OutputCount) and reject repeats outside the sequence.
    /// </summary>
    [Test]
    public void TrfFormatters_AnyFinderOutput_FormatsWithoutError()
    {
        var rng = new Random(Seed + 4);
        for (int trial = 0; trial < 80; trial++)
        {
            string s = B04AuditProperties.TandemRich(rng, rng.Next(1, 500), withN: trial % 2 == 0);
            var p = TandemRepeatsFinderParameters.Recommended with { FlankLength = trial % 3 == 0 ? 50 : 0, MatchProbability = trial % 2 == 0 ? 80 : 75 };
            var rows = RepeatFinder.FindApproximateTandemRepeats(s, p, out int count);
            Assert.That(RepeatFinder.FormatTrfDatLines(s, rows, "seq name", p), Does.Contain("Sequence: seq name"));
            Assert.That(RepeatFinder.FormatTrfDatLines(s, rows, "seq", p, TrfDatLayout.Ngs).Split('\n').Count(l => l.Length > 0),
                Is.EqualTo(rows.Count == 0 ? 0 : rows.Count + 1));
            Assert.That(RepeatFinder.FormatTrfHtmlTables(s, rows, "seq", p, "f"), Has.Count.EqualTo(Math.Max(1, (rows.Count + 119) / 120)));
            Assert.That(RepeatFinder.FormatTrfAlignmentPages(s, rows, "seq", p, "f", count), Is.Not.Empty);
            Assert.That(RepeatFinder.FormatTrfAlignmentPages(s, rows, "seq", p, "f"), Is.Not.Empty);
            Assert.That(RepeatFinder.FormatTrfHtmlSummary(new[] { ("seq", rows.Count), ("other", 0) }, p, "f").Html, Does.Contain("seq"));
        }

        var outside = new ApproximateTandemRepeatResult(5, 10, 2, 2, "AC", 5, 100, 0, 20);
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FormatTrfDatLines("ACACACAC", new[] { outside }, "x", TandemRepeatsFinderParameters.Recommended));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.MaskApproximateTandemRepeats("ACACACAC", new[] { outside }));
    }

    /// <summary>
    /// F47/F58/F59: degenerate repeats with random minLength / k / maxLength / minSpacing in every distance, strand,
    /// reporting and compatibility mode on adversarial input return rows inside the sequence that respect the filters.
    /// </summary>
    [Test]
    public void DegenerateRepeats_RandomOptions_RowsRespectFilters()
    {
        var rng = new Random(Seed + 5);
        for (int trial = 0; trial < 200; trial++)
        {
            string s = Adversarial(rng);
            if (s.Length > 150) s = s[..150];
            int k = rng.Next(1, 4), minLength = k + rng.Next(1, 12);
            int maxLength = rng.Next(2) == 0 ? int.MaxValue : minLength + rng.Next(0, 30);
            int minSpacing = rng.Next(3) == 0 ? int.MinValue : rng.Next(-20, 20);
            var distance = rng.Next(2) == 0 ? ApproximateRepeatDistance.Edit : ApproximateRepeatDistance.Hamming;
            bool pal = rng.Next(2) == 0;
            var reporting = rng.Next(2) == 0 ? DegenerateRepeatReporting.AllMaximal : DegenerateRepeatReporting.BestPerSeed;
            bool compat = rng.Next(2) == 0;
            var rows = RepeatFinder.FindDegenerateRepeats(s, minLength, k, distance, pal, maxLength, minSpacing, reporting, compat).ToList();
            foreach (var x in rows)
            {
                Assert.That(x.FirstPosition >= 0 && x.SecondPosition >= 0
                    && x.FirstPosition + x.FirstLength <= s.Length && x.SecondPosition + x.SecondLength <= s.Length
                    && x.FirstLength >= minLength && x.SecondLength >= minLength
                    && x.FirstLength <= maxLength && x.SecondLength <= maxLength
                    && x.Spacing >= minSpacing && x.Distance <= k && x.Distance >= 0, Is.True, $"{x} {s}");
            }
            var approx = RepeatFinder.FindApproximateDirectRepeats(s, minLength, k, maxLength, minSpacing, rng.Next(2) == 0, reporting, compat).ToList();
            Assert.That(approx.All(x => x.Mismatches <= k && x.Length >= minLength && x.SecondPosition + x.Length <= s.Length), Is.True, s);
        }
    }

    /// <summary>
    /// F48/F61: ParseMisaDefinition on random text either returns a valid map (sizes ≥ 1, minima ≥ 2, ordered) or throws
    /// ArgumentException; the scan accepts every valid map (unit sizes up to 30) on adversarial input in both modes.
    /// </summary>
    [Test]
    public void MisaDefinitionParserAndScan_RandomText_ValidMapOrArgumentException()
    {
        var rng = new Random(Seed + 6);
        const string chars = "0123456789--  ,\tdefx:()";
        int parsedCount = 0;
        for (int trial = 0; trial < 3000; trial++)
        {
            string text = Random(rng, chars, rng.Next(0, 30));
            IReadOnlyDictionary<int, int> map;
            try
            {
                map = RepeatFinder.ParseMisaDefinition(text);
            }
            catch (ArgumentException)
            {
                continue;
            }
            parsedCount++;
            Assert.That(map.Count > 0 && map.All(kv => kv.Key >= 1 && kv.Value >= 2) && map.Keys.SequenceEqual(map.Keys.OrderBy(x => x)), Is.True, text);
        }
        Assert.That(parsedCount, Is.GreaterThan(0));

        for (int trial = 0; trial < 150; trial++)
        {
            string s = Adversarial(rng, anySymbols: true);
            var map = Enumerable.Range(1, rng.Next(1, 6)).Select(_ => rng.Next(1, 31)).Distinct().ToDictionary(u => u, _ => rng.Next(2, 6));
            foreach (var mode in new[] { MicrosatelliteScanMode.MaximalRuns, MicrosatelliteScanMode.MisaRegex })
            {
                var ssrs = RepeatFinder.FindMicrosatellites(s, map, mode).ToList();
                Assert.That(ssrs.All(m => m.Position >= 0 && m.Position + m.TotalLength <= s.Length && map.ContainsKey(m.RepeatUnit.Length)
                    && m.RepeatCount >= map[m.RepeatUnit.Length]), Is.True, s);
                var summary = RepeatFinder.GetTandemRepeatSummary(s, map, mode);
                Assert.That(summary.TotalRepeats, Is.EqualTo(ssrs.Count), s);
                Assert.That(summary.PercentageOfSequence, Is.InRange(0.0, 100.0), s);
            }
        }
    }
}
