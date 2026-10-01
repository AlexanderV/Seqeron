using Seqeron.Genomics.Tests.Properties;

namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the behaviour added by the B04 completeness audit
/// (docs/Validation/review-2026-09/B04.md, F34–F65): longdust strand symmetry (F34), degenerate palindromic repeats
/// under reverse complement (F47), symbol-relabelling / case / U↔T invariance of the k-mer entropy corrections (F54),
/// fixed-alphabet LC, windowed LC and BBDuk entropy (F37–F39, F55), case invariance of the TRF parameter API with
/// its statistics and formatters (F41, F57, F64) and of the MISA scan mode compounds (F48). Fixed seeds; every relation
/// is exact.
///
/// Test Units: SEQ-COMPLEX-DUST-001, SEQ-COMPLEX-KMER-001, SEQ-COMPLEX-001, SEQ-COMPLEX-WINDOW-001, REP-DIRECT-001,
/// REP-APPROX-001, REP-STR-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Complexity")]
public class B04AuditMetamorphicTests
{
    private const int Seed = 20261001;

    private static char Complement(char c) => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c };

    private static string ReverseComplement(string s) => new(s.Reverse().Select(Complement).ToArray());

    // A bijection of {A, C, G, T} that is neither the identity nor the complement.
    private static string Relabel(string s) =>
        new(s.Select(c => c switch { 'A' => 'C', 'C' => 'G', 'G' => 'T', 'T' => 'A', _ => c }).ToArray());

    /// <summary>
    /// MR (F34): the default longdust result is the union of the forward scan and the reverse-complement scan mapped
    /// back, so the intervals of revcomp(S) are the mirror images [n − e, n − s) of those of S.
    /// </summary>
    [Test]
    public void Longdust_BothStrands_MirrorsUnderReverseComplement()
    {
        var rng = new Random(Seed);
        for (int trial = 0; trial < 120; trial++)
        {
            string s = B04AuditProperties.Segmented(rng, rng.Next(20, 900), trial % 3 == 0 ? "N" : "");
            int k = rng.Next(2, 8), w = rng.Next(20, 500);
            int n = s.Length;
            var forward = SequenceComplexity.FindLongdustRegions(s, k, w);
            var mirrored = SequenceComplexity.FindLongdustRegions(ReverseComplement(s), k, w)
                .Select(iv => (Start: n - iv.End, End: n - iv.Start)).OrderBy(iv => iv.Start).ToList();
            Assert.That(mirrored, Is.EqualTo(forward), $"k={k} w={w} s={s}");
        }
    }

    /// <summary>
    /// MR (F47): a palindromic (reverse-complement) degenerate repeat pairs instance A with revcomp(B); in revcomp(S)
    /// the same pair appears at the mirrored coordinates, so the set of unordered instance pairs (with distances) maps
    /// onto itself under p ↦ n − p − len. Checked for edit and Hamming distance, with wildcards. (A set, not a multiset:
    /// a pair with equal starts and different lengths is reported in both orientations, its mirror image only once.)
    /// </summary>
    [Test]
    public void DegenerateRepeats_Palindromic_MirrorUnderReverseComplement()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 120; trial++)
        {
            string s = B04AuditProperties.Segmented(rng, rng.Next(15, 90), trial % 2 == 0 ? "N" : "");
            int k = rng.Next(1, 3), minLength = k + rng.Next(2, 8);
            var distance = trial % 2 == 0 ? ApproximateRepeatDistance.Edit : ApproximateRepeatDistance.Hamming;
            int n = s.Length;

            static (int, int, int, int, int) Pair(int p1, int l1, int p2, int l2, int d) =>
                (p1, l1).CompareTo((p2, l2)) <= 0 ? (p1, l1, p2, l2, d) : (p2, l2, p1, l1, d);

            var direct = RepeatFinder.FindDegenerateRepeats(s, minLength, k, distance, true, int.MaxValue, int.MinValue)
                .Select(x => Pair(x.FirstPosition, x.FirstLength, x.SecondPosition, x.SecondLength, x.Distance))
                .Distinct().OrderBy(x => x).ToList();
            var mirrored = RepeatFinder.FindDegenerateRepeats(ReverseComplement(s), minLength, k, distance, true, int.MaxValue, int.MinValue)
                .Select(x => Pair(n - x.FirstPosition - x.FirstLength, x.FirstLength, n - x.SecondPosition - x.SecondLength, x.SecondLength, x.Distance))
                .Distinct().OrderBy(x => x).ToList();
            Assert.That(mirrored, Is.EqualTo(direct), $"{distance} k={k} l={minLength} s={s}");
        }
    }

    /// <summary>
    /// MR (F54): every k-mer entropy estimate (plug-in, Miller–Madow, Grassberger; raw and normalised) depends only on
    /// the multiset of k-mer counts, so it is unchanged by a bijective relabelling of the bases and by lower case.
    /// </summary>
    [Test]
    public void KmerEntropy_AllCorrections_InvariantUnderRelabellingAndCase()
    {
        var rng = new Random(Seed + 2);
        var corrections = new[] { KmerEntropyCorrection.None, KmerEntropyCorrection.MillerMadow, KmerEntropyCorrection.Grassberger };
        for (int trial = 0; trial < 300; trial++)
        {
            string s = B04AuditProperties.Segmented(rng, rng.Next(1, 200));
            int k = rng.Next(1, 7);
            foreach (var c in corrections)
            {
                foreach (bool normalize in new[] { false, true })
                {
                    double v = SequenceComplexity.CalculateKmerEntropy(s, k, c, normalize);
                    Assert.That(SequenceComplexity.CalculateKmerEntropy(Relabel(s), k, c, normalize), Is.EqualTo(v).Within(1e-12), $"{c} {normalize} {s}");
                    Assert.That(SequenceComplexity.CalculateKmerEntropy(s.ToLowerInvariant(), k, c, normalize), Is.EqualTo(v), $"{c} {normalize} {s}");
                }
            }
        }
    }

    /// <summary>
    /// MR (F37–F39, F55): fixed-alphabet LC, windowed LC / Shannon points and BBDuk entropy masking depend only on the
    /// equality pattern / counts of the bases, so a bijective relabelling leaves them unchanged; lower case and (for the
    /// BBDuk / windowed string paths) U for T give the same result.
    /// </summary>
    [Test]
    public void LinguisticWindowedAndBbduk_InvariantUnderRelabellingCaseAndU()
    {
        var rng = new Random(Seed + 3);
        for (int trial = 0; trial < 200; trial++)
        {
            string s = B04AuditProperties.Segmented(rng, rng.Next(1, 220));
            int m = rng.Next(1, 12), w = rng.Next(5, 50), step = rng.Next(1, 10), kb = rng.Next(1, 6);
            string rel = Relabel(s), variant = s.Replace('T', 'U').ToLowerInvariant();

            Assert.That(SequenceComplexity.CalculateLinguisticComplexity(rel, m, 4),
                Is.EqualTo(SequenceComplexity.CalculateLinguisticComplexity(s, m, 4)), s);

            var points = SequenceComplexity.CalculateWindowedComplexity(s, w, step, m).ToList();
            var relPoints = SequenceComplexity.CalculateWindowedComplexity(rel, w, step, m).ToList();
            // LC is an exact ratio of integer counts; the entropy sums the same terms in another order (≤ 1 ulp apart).
            Assert.That(relPoints.Select(p => (p.Position, p.WindowStart, p.WindowEnd, p.LinguisticComplexity)),
                Is.EqualTo(points.Select(p => (p.Position, p.WindowStart, p.WindowEnd, p.LinguisticComplexity))), s);
            Assert.That(relPoints.Zip(points).All(z => Math.Abs(z.First.ShannonEntropy - z.Second.ShannonEntropy) <= 1e-12), s);
            Assert.That(SequenceComplexity.CalculateWindowedComplexity(s.ToLowerInvariant(), w, step, m).ToList(), Is.EqualTo(points), s);

            double cutoff = rng.Next(0, 101) / 100.0;
            var bb = SequenceComplexity.FindLowEntropyRegionsBbduk(s, cutoff, w, kb).Select(r => (r.Start, r.End)).ToList();
            Assert.That(SequenceComplexity.FindLowEntropyRegionsBbduk(rel, cutoff, w, kb).Select(r => (r.Start, r.End)), Is.EqualTo(bb), s);
            Assert.That(SequenceComplexity.FindLowEntropyRegionsBbduk(variant, cutoff, w, kb).Select(r => (r.Start, r.End)), Is.EqualTo(bb), s);
        }
    }

    /// <summary>
    /// MR (F41, F57, F64): TRF upper-cases its input, so the parameter API returns identical results (statistics,
    /// alignment rows, flanks, output indices) for a lower-case copy, and the .dat lines and alignment pages are
    /// identical too.
    /// </summary>
    [Test]
    public void Trf_ParameterApiAndFormatters_InvariantUnderCase()
    {
        var rng = new Random(Seed + 4);
        var sets = new[]
        {
            TandemRepeatsFinderParameters.Recommended,
            TandemRepeatsFinderParameters.Recommended with { MatchProbability = 75, IndelProbability = 20, FlankLength = 50 },
            TandemRepeatsFinderParameters.Recommended with { MismatchPenalty = 3, IndelPenalty = 5, MinScore = 40, MaxPeriod = 200, EliminateRedundancy = false },
        };
        for (int trial = 0; trial < 60; trial++)
        {
            string s = B04AuditProperties.TandemRich(rng, rng.Next(30, 500), withN: trial % 2 == 0);
            var p = sets[trial % sets.Length];
            var upper = RepeatFinder.FindApproximateTandemRepeats(s, p, out int count);
            var lower = RepeatFinder.FindApproximateTandemRepeats(s.ToLowerInvariant(), p, out int lowerCount);
            Assert.That(lower, Is.EqualTo(upper), s);
            Assert.That(lowerCount, Is.EqualTo(count), s);
            Assert.That(RepeatFinder.FormatTrfDatLines(s.ToLowerInvariant(), lower, "x", p), Is.EqualTo(RepeatFinder.FormatTrfDatLines(s, upper, "x", p)), s);
            Assert.That(RepeatFinder.FormatTrfAlignmentPages(s.ToLowerInvariant(), lower, "x", p, "x", lowerCount),
                Is.EqualTo(RepeatFinder.FormatTrfAlignmentPages(s, upper, "x", p, "x", count)), s);
        }
    }

    /// <summary>
    /// MR (F48): misa.pl matches its SSR regex case-insensitively, so the MisaRegex SSR list and the compound
    /// microsatellites (start, end, type, components) do not change when the input is lower-cased; the same holds for the
    /// maximal-run mode.
    /// </summary>
    [Test]
    public void MisaScanAndCompounds_InvariantUnderCase()
    {
        var rng = new Random(Seed + 5);
        var maps = new IReadOnlyDictionary<int, int>[]
        {
            RepeatFinder.MisaDefaultMinRepeats,
            new Dictionary<int, int> { [1] = 3, [2] = 2, [3] = 2, [4] = 2, [5] = 2, [6] = 2 },
            new Dictionary<int, int> { [2] = 3, [7] = 2, [9] = 2, [10] = 2 },
        };
        for (int trial = 0; trial < 150; trial++)
        {
            string s = B04AuditProperties.SsrRich(rng, rng.Next(10, 300), dirty: false);
            var map = maps[trial % maps.Length];
            int maxInterruption = rng.Next(0, 101);
            foreach (var mode in new[] { MicrosatelliteScanMode.MisaRegex, MicrosatelliteScanMode.MaximalRuns })
            {
                Assert.That(RepeatFinder.FindMicrosatellites(s.ToLowerInvariant(), map, mode).ToList(),
                    Is.EqualTo(RepeatFinder.FindMicrosatellites(s, map, mode).ToList()), s);
                var up = RepeatFinder.FindCompoundMicrosatellites(s, map, maxInterruption, mode);
                var low = RepeatFinder.FindCompoundMicrosatellites(s.ToLowerInvariant(), map, maxInterruption, mode);
                Assert.That(low.Select(c => (c.Start, c.End, c.MisaType, c.Components.Count)),
                    Is.EqualTo(up.Select(c => (c.Start, c.End, c.MisaType, c.Components.Count))), s);
                Assert.That(low.SelectMany(c => c.Components), Is.EqualTo(up.SelectMany(c => c.Components)), s);
            }
        }
    }
}
